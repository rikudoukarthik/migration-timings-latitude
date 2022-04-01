##### data processing ----------

# reading, processing, modifying and setting up data for migration-timings-latitude specifically

data_process <- function(rawpath, senspath) {
  
require(lubridate)
require(tidyverse)

preimp <- c("CATEGORY","COMMON.NAME","SUBSPECIES.COMMON.NAME","OBSERVATION.COUNT", "LOCALITY.ID",
            "LOCALITY.TYPE","STATE","COUNTY","LATITUDE","LONGITUDE","OBSERVATION.DATE",
            "TIME.OBSERVATIONS.STARTED","OBSERVER.ID","PROTOCOL.TYPE","DURATION.MINUTES",
            "EFFORT.DISTANCE.KM","LOCALITY","NUMBER.OBSERVERS","ALL.SPECIES.REPORTED",
            "GROUP.IDENTIFIER","SAMPLING.EVENT.IDENTIFIER","TRIP.COMMENTS")

nms <- names(read.delim(rawpath, nrows = 1, sep = "\t", header = T, quote = "", 
                        stringsAsFactors = F, na.strings = c(""," ", NA)))
nms[!(nms %in% preimp)] <- "NULL"
nms[nms %in% preimp] <- NA
data <- read.delim(rawpath, colClasses = nms, sep = "\t", header = T, quote = "",
                   stringsAsFactors = F, na.strings = c(""," ",NA)) 

nms1 <- names(read.delim(senspath, nrows = 1, sep = "\t", header = T, quote = "", 
                         stringsAsFactors = F, na.strings = c(""," ", NA)))
nms1[!(nms1 %in% preimp)] <- "NULL"
nms1[nms1 %in% preimp] <- NA
senssp <- read.delim(senspath, colClasses = nms1, sep = "\t", header = T, quote = "",
                     stringsAsFactors = F, na.strings = c(""," ",NA))


data <- bind_rows(data, senssp) 


met_week <- function(dates) {
  require(lubridate)
  normal_year <- c((0:363 %/% 7 + 1), 52)
  leap_year   <- c(normal_year[1:59], 9, normal_year[60:365])
  year_day    <- yday(dates)
  return(ifelse(leap_year(dates), leap_year[year_day], normal_year[year_day])) 
}


### adding useful columns 
data <- data %>% 
  mutate(GROUP.ID = ifelse(is.na(GROUP.IDENTIFIER), SAMPLING.EVENT.IDENTIFIER, 
                           GROUP.IDENTIFIER), 
         OBSERVATION.DATE = as.Date(OBSERVATION.DATE), 
         YEAR = year(OBSERVATION.DATE), 
         MONTH = month(OBSERVATION.DATE),
         DAY.M = day(OBSERVATION.DATE),
         DAY.Y = yday(OBSERVATION.DATE),
         WEEK.Y = met_week(OBSERVATION.DATE),
         M.YEAR = if_else(DAY.Y <= 151, YEAR-1, YEAR), # from 1st June to 31st May
         WEEK.MY = if_else(WEEK.Y > 21, WEEK.Y-21, 52-(21-WEEK.Y)),
         # 1 deg resolution (since 1deg ~ 100km)
         LATITUDE1 = round(LATITUDE, 0)) 
  

### removing duplicate lists
data <- data %>%
  group_by(GROUP.ID, COMMON.NAME) %>%
  arrange(SAMPLING.EVENT.IDENTIFIER) %>% 
  slice(1) %>% ungroup()



assign("data", data, .GlobalEnv)

save(data, file = "data/data.RData")

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
  


  
  ### new observer data (to calculate no. of new observers metric) #######
  
  new_obsr_data <- data %>% 
    select(c("YEAR", "MONTH", "STATE", "SAMPLING.EVENT.IDENTIFIER",
             "LAST.EDITED.DATE", "OBSERVATION.DATE", "OBSERVER.ID")) %>% 
    mutate(LAST.EDITED.DATE = ymd_hms(LAST.EDITED.DATE)) %>% 
    group_by(OBSERVER.ID) %>% 
    arrange(LAST.EDITED.DATE) %>% 
    ungroup() %>% 
    distinct(OBSERVER.ID, .keep_all = TRUE) %>%
    mutate(YEAR = year(LAST.EDITED.DATE),
           MONTH = month(LAST.EDITED.DATE)) %>% 
    filter(YEAR >= 2014) %>% 
    rename(LE.YEAR = YEAR,
           LE.MONTH = MONTH) %>% 
    mutate(YEAR = year(OBSERVATION.DATE), 
           MONTH = month(OBSERVATION.DATE))
  
  # filtering
  new_obsr_data <- new_obsr_data %>% anti_join(filtGA) 
  save(new_obsr_data, file = "data/new_obsr_data.RData")
  
  
  
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
  
  # choose checklists without info on duration with 3 or fewer species
  temp <- data %>%
    filter(ALL.SPECIES.REPORTED == 1, PROTOCOL.TYPE != "Incidental") %>%
    group_by(GROUP.ID) %>% slice(1) %>%
    filter(NO.SP <= 3, is.na(DURATION.MINUTES)) %>%
    distinct(GROUP.ID)
  
  # exclude records based on various criteria 
  data <- data %>%
    mutate(ALL.SPECIES.REPORTED = 
             case_when(ALL.SPECIES.REPORTED == 1 & 
                         (GROUP.ID %in% temp | 
                            SPEED > maxvel |
                            (SUT < minsut & NO.SP <= 3) | 
                            PROTOCOL.TYPE == "Incidental" | 
                            (!is.na(TIME.D) & ((TIME.D <= 4 & END <= 4) | 
                                                 (TIME.D >= 20 & END <= 28)
                            )
                            )
                         ) ~ 0, 
                       ALL.SPECIES.REPORTED == 0 ~ 0,
                       TRUE ~ 1)) %>% 
    select(-SPEED, -SUT, -MIN, -END) %>% 
    filter(ALL.SPECIES.REPORTED == 1)
  
  assign("data", data, .GlobalEnv)
  
  save(data, file = "data/data.RData")
}



##### weekly reporting frequencies ----------

# calculating weekly reporting frequencies

sp_repfreq <- function(data){
  
dataRF <- data %>% 
  filter(ALL.SPECIES.REPORTED == 1) %>% # complete lists only
  group_by(LATITUDE1, WEEK.MY) %>% 
  mutate(TOT.LISTS = n_distinct(GROUP.ID)) %>% 
  group_by(LATITUDE1, WEEK.MY, TOT.LISTS, COMMON.NAME) %>% 
  summarise(NO.LISTS = n_distinct(GROUP.ID),
            REP.FREQ = NO.LISTS/TOT.LISTS)

assign("dataRF", dataRF, .GlobalEnv)

save(dataRF, file = "data/dataRF.RData")

}

##### select certain species and subspecies ----------

sp_select <- function(data){
  
  # for GrWa, excluding data from Himalayan states and from WB eastwards, 
  # to focus on ssp. viridanus with extralimital breeding range
  
  data0 <- data %>% 
    filter((CATEGORY == "species" & COMMON.NAME %in% c("Blyth's Reed Warbler", 
                                                       "Brown-breasted Flycatcher", 
                                                       "Greenish Warbler")) |
             # (CATEGORY == "issf" & str_detect(SUBSPECIES.COMMON.NAME, "viridanus")) |
             (CATEGORY == "issf" & COMMON.NAME %in% c("Blyth's Reed Warbler", 
                                                         "Brown-breasted Flycatcher"))) %>%
    filter(!(COMMON.NAME == "Greenish Warbler" & 
               STATE %in% c("West Bengal", "Sikkim", "Meghalaya", "Assam", "Arunachal Pradesh",
                            "Nagaland", "Manipur", "Mizoram", "Tripura",
                            "Jammu and Kashmir", "Himachal Pradesh", "Uttarakhand", 
                            "Uttar Pradesh", "Bihar"))) # these two states?
    
  
  assign("data0", data0, .GlobalEnv)
  
  return(data0)
  
  }