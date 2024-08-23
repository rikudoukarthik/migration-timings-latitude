##### data processing ----------

# reading, processing, modifying and setting up data for migration-timings-latitude specifically

data_process <- function(rawpath, senspath, data_obj) {
  
require(lubridate)
require(tidyverse)

# preimp <- c("CATEGORY","COMMON.NAME","SUBSPECIES.COMMON.NAME","OBSERVATION.COUNT","LOCALITY.ID",
#             "LOCALITY.TYPE","STATE","COUNTY","LATITUDE","LONGITUDE","OBSERVATION.DATE",
#             "TIME.OBSERVATIONS.STARTED","OBSERVER.ID","PROTOCOL.TYPE","DURATION.MINUTES",
#             "EFFORT.DISTANCE.KM","LOCALITY","NUMBER.OBSERVERS","ALL.SPECIES.REPORTED",
#             "GROUP.IDENTIFIER","SAMPLING.EVENT.IDENTIFIER","TRIP.COMMENTS")

# nms <- names(read.delim(rawpath, nrows = 1, sep = "\t", header = T, quote = "", 
#                         stringsAsFactors = F, na.strings = c(""," ", NA)))
# nms[!(nms %in% preimp)] <- "NULL"
# nms[nms %in% preimp] <- NA
# data <- read.delim(rawpath, colClasses = nms, sep = "\t", header = T, quote = "",
#                    stringsAsFactors = F, na.strings = c(""," ",NA)) 

# nms1 <- names(read.delim(senspath, nrows = 1, sep = "\t", header = T, quote = "", 
#                          stringsAsFactors = F, na.strings = c(""," ", NA)))
# nms1[!(nms1 %in% preimp)] <- "NULL"
# nms1[nms1 %in% preimp] <- NA
# senssp <- read.delim(senspath, colClasses = nms1, sep = "\t", header = T, quote = "",
#                      stringsAsFactors = F, na.strings = c(""," ",NA))


# data <- bind_rows(data, senssp) 


### adding useful columns 
data <- data_obj %>% 
  mutate(GROUP.ID = ifelse(is.na(GROUP.IDENTIFIER), SAMPLING.EVENT.IDENTIFIER, 
                           GROUP.IDENTIFIER), 
         OBSERVATION.DATE = as_date(OBSERVATION.DATE), 
         YEAR = year(OBSERVATION.DATE), 
         MONTH = month(OBSERVATION.DATE),
         DAY.M = day(OBSERVATION.DATE),
         DAY.Y = yday(OBSERVATION.DATE),
         WEEK.Y = week(OBSERVATION.DATE),
         M.YEAR = if_else(MONTH <= 5, YEAR-1, YEAR), # from 1st June to 31st May
         WEEK.MY = if_else(WEEK.Y > 21, WEEK.Y-21, 52-(21-WEEK.Y)),
         # 1 deg resolution (since 1deg ~ 100km)
         LATITUDE1 = round(LATITUDE, 0)) 
  

### removing duplicate lists
data <- data %>%
  group_by(GROUP.ID, COMMON.NAME) %>%
  arrange(SAMPLING.EVENT.IDENTIFIER) %>% 
  distinct(GROUP.ID, COMMON.NAME, .keep_all = TRUE) %>% 
  ungroup()


save(data, file = "data/data.RData")
  
return(data)

}


##### data filters ----------

# data file path
# path to list of group accounts to be filtered out
# path to classification of year-month as COVID categories


dataqual_filt <- function(data, groupaccspath, maxvel = 20, minsut = 2){
  
  # data object already in environment

  require(tidyverse)
  require(lubridate)
  
  ### list of group accounts to be filtered
  groupaccs <- read.csv(groupaccspath, 
                        na.strings = c(""," ",NA), quote = "", header = T, 
                        nrows = 401)  # excluding empty cells
  groupaccs <- groupaccs %>% 
    mutate(CATEGORY = case_when(GA.1 == 1 ~ "GA.1", 
                                GA.2 == 1 ~ "GA.2", 
                                TRUE ~ "NG"))
  filtGA <- groupaccs %>% filter(CATEGORY == "GA.1") %>% select(OBSERVER.ID)
  

  
  ### main data filtering ######
  
  data <- data %>% 
    filter(YEAR >= 2014) %>% # retaining only data from 2014 onward
    anti_join(filtGA) %>% # removing data from group accounts
    group_by(GROUP.ID) %>% 
    mutate(NO.SP = n_distinct(COMMON.NAME)) %>%
    ungroup() %>% 
    mutate(TIME.D = hour(as_datetime(paste(OBSERVATION.DATE,
                                           TIME.OBSERVATIONS.STARTED))),
           MIN = minute(as_datetime(paste(OBSERVATION.DATE,
                                          TIME.OBSERVATIONS.STARTED))),
           SPEED = EFFORT.DISTANCE.KM*60/DURATION.MINUTES, # kmph
           SUT = NO.SP*60/DURATION.MINUTES, # species per hour
           # calculate hour checklist ended
           END = floor((TIME.D*60 + MIN + DURATION.MINUTES)/60))
  

save(data, file = "data/data.RData")
  
return(data)
  
}



##### weekly reporting frequencies ----------

# calculating weekly reporting frequencies

calc_repfreq <- function(data){
  
  lat_week_totlists <- data %>% 
    filter(ALL.SPECIES.REPORTED == 1) %>% # complete lists only
    group_by(LAT.BOX, WEEK.MY) %>% 
    reframe(TOT.LISTS = n_distinct(GROUP.ID)) %>% 
    complete(WEEK.MY, nesting(LAT.BOX), fill = list(TOT.LISTS = 0)) |> 
    arrange(desc(LAT.BOX), WEEK.MY) |> 
    relocate(LAT.BOX)
  
  spec_lat_week <- data |> 
    distinct(COMMON.NAME, LAT.BOX, WEEK.MY) |> 
    # for each species, we want all latitude-week combos (only for relevant lat)
    group_by(COMMON.NAME) |> 
    left_join(lat_week_totlists, by = c("LAT.BOX", "WEEK.MY")) |> 
    complete(LAT.BOX, WEEK.MY, fill = list(TOT.LISTS = 0)) |> 
    ungroup() |> 
    arrange(COMMON.NAME, desc(LAT.BOX), WEEK.MY)
  
  data_repfreq <- spec_lat_week |> 
    left_join(data |> filter(ALL.SPECIES.REPORTED == 1), 
              by = c("COMMON.NAME", "LAT.BOX", "WEEK.MY")) |> 
    group_by(COMMON.NAME, LAT.BOX, WEEK.MY, TOT.LISTS) %>% 
    reframe(NO.LISTS = n_distinct(GROUP.ID, na.rm = TRUE)) |> 
    filter(TOT.LISTS > 0) |> 
    mutate(REP.FREQ = 100*NO.LISTS/TOT.LISTS)
  
  # saving data
  save(lat_week_totlists, spec_lat_week, data_repfreq,
        file = "data/data_repfreq.RData")

  return(data_repfreq)

}

# filter for species of interest & prepare winter and full datasets ----------

prep_data_spec <- function(data){

  # select species with certain migratory status acc. to SoIB
  soib_migrants <- read_csv("data/SoIB_mapping_2023.csv") |> 
    filter(Migratory.Status.Within.India %in% c("Winter Migrant",
                                                "Winter Migrant & Localized Summer Migrant"))
  
  # for GrWa, excluding data from Himalayan states and from WB eastwards, 
  # to focus on ssp. viridanus with extralimital breeding range
  
  data_full <- data %>% 
    filter(COMMON.NAME %in% soib_migrants$eBird.English.Name.2023) %>%
    # remove migratorily complicated subspecies
    filter(!(COMMON.NAME == "Greenish Warbler" & 
               STATE %in% c("West Bengal", "Sikkim", "Meghalaya", "Assam", "Arunachal Pradesh",
                            "Nagaland", "Manipur", "Mizoram", "Tripura",
                            "Jammu and Kashmir", "Himachal Pradesh", "Uttarakhand", 
                            "Uttar Pradesh", "Bihar"))) # these two states?
    
  # if present outside winter (as well as in winter), then ignore
  non_winter <- data_full |> 
    filter(WEEK.MY %in% 1:10) %>% 
    group_by(COMMON.NAME, LATITUDE1) |> 
    reframe(N.NONWIN.WEEKS = n_distinct(WEEK.MY)) |> 
    mutate(TRUE.NONWINTER = case_when(N.NONWIN.WEEKS >= 3 ~ TRUE,
                                      TRUE ~ FALSE)) %>% 
    distinct(COMMON.NAME, LATITUDE1, TRUE.NONWINTER)
  
  win_spec_info <- data_full |> 
    # winter dataset
    filter(MONTH %in% c(12, 1)) |> 
    # only consider latitudes where species "truly wintering"
    # i.e., present at least 6 (of 11) winter weeks
    group_by(COMMON.NAME, LATITUDE1) |> 
    reframe(N.WIN.WEEKS = n_distinct(WEEK.MY)) |> 
    mutate(TRUE.WINTER = case_when(N.WIN.WEEKS >= 6 ~ TRUE,
                                   TRUE ~ FALSE)) |> 
    arrange(COMMON.NAME, desc(LATITUDE1)) |> 
    group_by(COMMON.NAME) |> 
    # index of lat box for each species, for lat gap calculations below
    mutate(INDEX = row_number()) |> 
    # filter for true winterers
    # (also remove species present in winter & non-winter)
    left_join(non_winter, by = c("COMMON.NAME", "LATITUDE1")) %>% 
    filter(TRUE.WINTER == TRUE & 
             (TRUE.NONWINTER == FALSE | is.na(TRUE.NONWINTER))) |> 
    # filter out non-consecutive lat boxes for each species
    mutate(DIFF.LAG = INDEX - lag(INDEX, 1),
           DIFF.LEAD = -(INDEX - lead(INDEX, 1))) |> 
    # this allows for gaps but not gaps with only 1 valid row
    filter(DIFF.LAG == 1 | DIFF.LEAD == 1)  |> 
    # exclude species with large (> 1 deg) latitudanl gaps in winter range
    mutate(DIFF.LAG = INDEX - lag(INDEX, 1),
           DIFF.LEAD = -(INDEX - lead(INDEX, 1))) |> 
    mutate(DIFF.LAG = replace_na(DIFF.LAG, 0),
           DIFF.LEAD = replace_na(DIFF.LEAD, 0)) |> 
    mutate(MAX.GAP = max(max(DIFF.LAG), max(DIFF.LEAD)),
           SPEC.LARGE.GAP = case_when((MAX.GAP == 1) | 
                                          (n_distinct(LATITUDE1)/MAX.GAP >= 4) ~ FALSE,
                                      TRUE ~ TRUE)) |> 
    ungroup() |> 
    filter(SPEC.LARGE.GAP == FALSE) |> 
    distinct(COMMON.NAME, LATITUDE1, N.WIN.WEEKS) |>
    # exclude species having lat breadth below threshold
    group_by(COMMON.NAME) |> 
    mutate(N.LAT = n_distinct(LATITUDE1),
           LAT.MAX = max(LATITUDE1),
           LAT.MIN = min(LATITUDE1),
           LAT.MED = median(LATITUDE1),
           # +1 to account for both northern and southern limit boxes
           # (to match with N.LAT)
           LAT.BREADTH = (LAT.MAX - LAT.MIN) + 1) |> 
    ungroup() |> 
    # latitudinal breadth threshold 
    # (using no. and not breadth to account for gaps)
    filter(N.LAT >= 5)  |> 
    # rename lat column
    rename(LAT.BOX = LATITUDE1) |> 
    relocate(LAT.BOX, N.WIN.WEEKS, .after = last_col()) |> 
    # join useful SoIB columns
    left_join(soib_migrants |> 
                dplyr::select(eBird.English.Name.2023,
                              Order, Family, Diet.Guild, Habitat.Specialization, Migratory.Status.Within.India,
                              all_of(starts_with("SoIBv2."))),
              by = c("COMMON.NAME" = "eBird.English.Name.2023"))

  
  # prepare full dataset

  data_full <- data_full |> 
    rename(LAT.BOX = LATITUDE1) |> 
    right_join(win_spec_info |> 
                dplyr::select(COMMON.NAME, LAT.BOX, N.LAT, LAT.BREADTH, 
                              LAT.MIN, LAT.MAX, LAT.MED),
               by = c("COMMON.NAME", "LAT.BOX"))
  
  # calculate reporting frequencies for each species
  # per lat box per week
  data_repfreq <- calc_repfreq(data_full)
  
  # save data objects
  save(data_full, data_repfreq, win_spec_info,
        file = "data/dataforanalysis.RData")
  
}

# calculate synchronicity measures --------------------------------------------------

calc_sync_meas <- function(data_grouped) {
  
  # initial trials 
  
  # reframe(TAU = cor.test(LAT.BOX, DAY.MY, method = "kendall")$estimate) %>% 
  # reframe(PW.DIST.MEAN = stats::dist(DAY.MY) %>% as.numeric() %>% mean(),
  #         PW.DIST.SD = stats::dist(DAY.MY) %>% as.numeric() %>% sd()) %>% 
  # reframe(SLOPE = (first(LAT.SCALE) - last(LAT.SCALE))/(first(DAY.SCALE) - last(DAY.SCALE))) %>% 
  # R2 = (lm(DAY.MY ~ LAT.BOX, data = cur_group()) %>% summary())$r.squared,
  
  data_calc <- data_grouped %>% 
    reframe(SLOPE = summary(lm(DAY.MY ~ LAT.BOX, data = cur_group()))$coef[2, 1],
            SE = summary(lm(DAY.MY ~ LAT.BOX, data = cur_group()))$coef[2, 2]) %>% 
    mutate(
      SYNC.SLOPE = ((1/SLOPE) %>% abs() %>% log() %>% scale())[,1], # quickness
      SYNC.TIGHT = -log(SE) # uniformity
    ) 
  
  return(data_calc)
  
}


# plot lat on Y and day on X --------------------------------------------------------

gg_lat_day <- function(data, date_type, ill = TRUE, sync = FALSE) {
  
  col_theme <- if (date_type == "ARR") "#006994" else if (date_type == "DEP") "#b22222"
  
  text_y <- if (date_type == "ARR") 60 else if (date_type == "DEP") 260
  text_x <- if (date_type == "ARR") 15 else if (date_type == "DEP") 30

  
  data_filt <- data %>% 
    filter(DATE.TYPE == date_type) %>% 
    {if (sync == TRUE) {
      arrange(., SYNC.SLOPE, SYNC.TIGHT)
    } else {
      .
    }}
  
  
  plot <- data_filt %>% 
    # arrange in order of synchronicity
    {if (sync == TRUE) {
      mutate(.,
             COMMON.NAME = factor(COMMON.NAME, levels = unique(data_filt$COMMON.NAME)),
             TEXT.LABEL = glue(
               "Quick = {round(SYNC.SLOPE, 3)}\nTight = {round(SYNC.TIGHT, 3)}"
               ))
    } else {
      .
    }} %>% 
    {if (ill == TRUE) {
      filter(., COMMON.NAME %in% illust_species)
    } else {
      .
    }} %>% 
    ggplot(aes(x = LAT.BOX, y = DAY.MY)) + 
    theme(strip.placement = "outside") +
    geom_point(size = 2, alpha = 0.75, col = col_theme, fill = col_theme) +
    geom_smooth(method = "lm", se = TRUE, col = "black", fill = "black") +
    {if (sync == TRUE) {
      geom_text(mapping = aes(x = text_x, y = text_y, label = TEXT.LABEL), size = 3)
    }} +
    facet_wrap(~ COMMON.NAME, 
               nrow = if (ill == TRUE) 5 else 9) +
    scale_x_continuous(breaks = seq(0, 40, 5)) +
    scale_y_continuous(breaks = seq(0, 400, 30)) +
    coord_flip() 
    
  return(plot)
  
}

gg_lat_day_sync <- function(data, date_type, ill = TRUE) {
  
  col_theme <- if (date_type == "ARR") "#006994" else if (date_type == "DEP") "#b22222"
  
  text_y <- if (date_type == "ARR") 60 else if (date_type == "DEP") 260
  text_x <- if (date_type == "ARR") 15 else if (date_type == "DEP") 30

  
  data_filt <- data %>% 
    filter(DATE.TYPE == date_type) %>% 
    {if (date_type == "ARR") {
      arrange(., SYNC.IND)
    } else if (date_type == "DEP") {
      arrange(., desc(SYNC.IND))
    }}

  
  plot <- data_filt %>% 
    # arrange in order of synchronicity
    mutate(COMMON.NAME = factor(COMMON.NAME, levels = unique(data_filt$COMMON.NAME)),
           TEXT.LABEL = glue("Slope = {round(SYNC.IND, 3)}")) %>% 
    {if (ill == TRUE) {
      filter(., COMMON.NAME %in% illust_species)
    } else {
      .
    }} %>% 
    ggplot(aes(x = LAT.BOX, y = DAY.MY)) + 
    theme(strip.placement = "outside") +
    geom_point(size = 2, alpha = 0.75, col = col_theme, fill = col_theme) +
    geom_smooth(method = "lm", se = TRUE, col = "black", fill = "black") +
    geom_text(mapping = aes(x = text_x, y = text_y, label = TEXT.LABEL), size = 3) +
    facet_wrap(~ COMMON.NAME, 
               nrow = if (ill == TRUE) 5 else 9) +
    scale_x_continuous(breaks = seq(0, 40, 5)) +
    scale_y_continuous(breaks = seq(0, 400, 30)) +
    coord_flip() 
  
  return(plot)
  
}

# get relative weights from PCA -----------------------------------------------------

get_pca_prop <- function(data) {
  
  # variables should be rows and samples should be columns
  sync_mat <- data %>% 
    mutate(SPEC.MIG = glue("{DATE.TYPE}_{COMMON.NAME}")) %>% 
    dplyr::select(SPEC.MIG, SYNC.SLOPE, SYNC.TIGHT) %>% 
    pivot_longer(cols = c("SYNC.SLOPE", "SYNC.TIGHT"), 
                 names_to = "SYNC.VAR", values_to = "VALUES") %>% 
    pivot_wider(names_from = "SPEC.MIG", values_from = "VALUES") %>% 
    mutate(across(c(everything(), -SYNC.VAR), ~ as.numeric(.))) %>% 
    as.matrix()
  
  rownames(sync_mat) <- sync_mat[, 1]
  sync_mat <- sync_mat[, -1]
  class(sync_mat) <- "numeric"
  
  sync_pca <- prcomp(t(sync_mat), scale = TRUE)
  
  sync_var <- sync_pca$sdev^2
  sync_var_prop <- round(sync_var/sum(sync_var), 3)
  
  return(sync_var_prop)
  
}

