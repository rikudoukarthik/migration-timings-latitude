# packages ####
library(auk)
library(tidyverse)
library(lubridate)
library(patchwork)
library(MASS)


# task ####

# setting input and output files
f_in <- "ebd_IN_relDec-2020.txt"
f_out <- "ebd_filtered_BCItask.txt"

# importing data into R. takes ~2min on my system
dataimport <- f_in %>% auk_ebd() %>% 
  # setting filter for the three species
  auk_species(species = c("Greenish Warbler", 
                          "Blyth's Reed Warbler", 
                          "Brown-breasted Flycatcher")) %>% 
  # running filter
  auk_filter(file = f_out) %>% 
  # reading text file into R
  read_ebd() %>% glimpse()



# creating a function to convert dates to week-of-year, 
# instead of using lubridate::week() which results in 53 weeks 
# (starts week from Monday rather than 1st)
met_week <- function(dates) {
  normal_year <- c((0:363 %/% 7 + 1), 52)
  leap_year   <- c(normal_year[1:59], 9, normal_year[60:365])
  year_day    <- lubridate::yday(dates)
  
  return(ifelse(lubridate::leap_year(dates), leap_year[year_day], normal_year[year_day]))
}


# 129652 rows and 45 columns. Need to remove a lot of the irrelevant columns.

# not ideal to name object "data" but for ease
data <- dataimport %>% select(-c(1:5,9:24,28:29,31:37,39:44)) %>% 
  mutate(days = yday(observation_date), weeks = met_week(observation_date)) %>% 
  group_by(common_name) %>% 
  mutate(week26plus = ifelse(weeks>26, weeks-26, 52-(26-weeks))) 
# "week26plus" for convenience and intuitiveness and to find max-min


# Plot1: visualising change in reports with time and latitude
plot1 <- ggplot(data, aes(weeks, latitude)) + 
  scale_x_continuous(breaks=(seq(0,52,4))) +
  scale_y_continuous(breaks=(seq(0,40,5))) +
  labs(title="Reports of migrant species varying with timing and latitude") +
  geom_point(colour="dark blue", alpha=1/3) +
  geom_vline(xintercept=26) +
  annotate("text", x=24, y=10, label = "Jun", angle=90) +
  annotate("text", x=27, y=10, label = "Jul", angle=90) +
  facet_wrap(~common_name) + theme_bw()


# Plot2: visualising change in reports with winter-time and latitude
plot2 <- ggplot(data, aes(week26plus, latitude)) + 
  scale_x_continuous(breaks=(seq(0,52,4))) +
  scale_y_continuous(breaks=(seq(0,40,5))) +
  labs(title="Reports of migrant species varying with winter-timing and latitude",
       subtitle="X-axis represents number of weeks after the 26th week of a year") +
  geom_point(colour="dark red", alpha=1/3) +
  geom_vline(xintercept=26) +
  annotate("text", x=24, y=10, label = "Dec", angle=90) +
  annotate("text", x=27, y=10, label = "Jan", angle=90) +
  facet_wrap(~common_name) + theme_bw()



# investigating first arrivals and last departures
# ("dataX" objects from 2 to 6 are in the "trial-and-error code" section)

data7 <- data %>% 
  mutate(latitude=round(latitude), # since 1deg ~ 100km
         rounding = factor("1deg (~100km)")) %>% # as factor to help with facet
  group_by(common_name, rounding, latitude) %>% 
  summarise(arrival = min(week26plus), departure = max(week26plus)) %>% 
  pivot_longer(cols = arrival:departure,
               names_to = "date_type", values_to = "weeks26plus") %>% 
  mutate(date_type = factor(date_type, levels = c("departure","arrival")))

data8 <- data %>% 
  mutate(latitude=round(latitude,1), # since 0.1deg ~ 10km
         rounding = factor("0.1deg (~10km)")) %>% # as factor to help with facet
  group_by(common_name, rounding, latitude) %>% 
  summarise(arrival = min(week26plus), departure = max(week26plus)) %>% 
  pivot_longer(cols = arrival:departure,
               names_to = "date_type", values_to = "weeks26plus") %>% 
  mutate(date_type = factor(date_type, levels = c("departure","arrival")))

data9 <- bind_rows(data8, data7)


# Plot3: visualising dates at two latitude scales and across the three species
plot3 <- ggplot(data9) + 
  theme_bw() + theme(strip.placement = "outside") +
  geom_point(aes(x=weeks26plus, 
                 y=latitude,
                 colour=date_type),
             alpha=2/3) +
  scale_colour_manual(values=c("maroon","sky blue")) +
  geom_vline(aes(xintercept=26)) +
  geom_text(aes(x=26, y=10, label="Dec"), nudge_x=-2, angle=90) +
  geom_text(aes(x=26, y=10, label="Jan"), nudge_x=1, angle=90) +
  facet_grid(rows = vars(data9$rounding), 
             cols = vars(data9$common_name),
             switch = "y") +
  scale_y_continuous(breaks=(seq(0,40,5))) +
  ggtitle("Variation of first arrival and last departure dates with latitude") 


# Plot4: visualising latitude and arrival/departure dates with smooth lines
plot4 <- 
  ggplot(data7, aes(latitude, weeks26plus, colour=common_name)) + 
  scale_x_continuous(breaks=(seq(0,52,4))) +
  scale_y_continuous(limits=(c(0,55)),
                     breaks=(seq(0,55,5))) +
  labs(title="1deg (~100km)") +
  scale_colour_manual(values=c("#FFCC66","#996633","#CCCC66"),
                      aesthetics = c("colour","fill"),
                      labels = c("BRWa","BBFl","GrWa"),
                      name = "Species Code") +
  geom_point() +
  geom_smooth(orientation="x", span=5, aes(fill=common_name)) +
  geom_hline(aes(yintercept=26)) +
  geom_text(aes(y=26, x=10, label="Jan"), colour="black", nudge_y=2) +
  geom_text(aes(y=26, x=10, label="Dec"), colour="black", nudge_y=-2) +
  facet_wrap(~data7$date_type, dir="v", strip.position = "left") +
  theme_bw() + theme(plot.title = element_text(hjust=0.5),
                     axis.title = element_blank()) +
  ggplot(data8, aes(latitude, weeks26plus, colour=common_name)) + 
  scale_x_continuous(breaks=(seq(0,52,4))) +
  scale_y_continuous(limits=(c(0,55)),
                     breaks=(seq(0,55,5))) +
  labs(title="0.1deg (~10km)") +
  scale_colour_manual(values=c("#FFCC66","#996633","#CCCC66"),
                      aesthetics = c("colour","fill"),
                      labels = c("BRWa","BBFl","GrWa"),
                      name = "Species Code") +
  geom_point() +
  geom_smooth(orientation="x", span=0.5, aes(fill=common_name)) +
  geom_hline(aes(yintercept=26)) +
  geom_text(aes(y=26, x=10, label="Jan"), colour="black", nudge_y=2) +
  geom_text(aes(y=26, x=10, label="Dec"), colour="black", nudge_y=-2) +
  facet_wrap(~data8$date_type, dir="v", strip.position = "left") +
  theme_bw() + theme(plot.title = element_text(hjust=0.5),
                     strip.text = element_blank(),
                     axis.title = element_blank()) +
  plot_annotation(title = "Variation of first arrival and last departure dates with latitude",
                  subtitle = ~atop("Y-axis represents number of weeks after the 26th week of a year", 
                                   "Facet columns show two different scales of latitude"),
                  caption = "latitude",
                  theme = theme(plot.title = element_text(hjust=0.5, size=18),
                                plot.subtitle = element_text(hjust=0.5),
                                plot.caption = element_text(hjust=0.47, size=14))) +
  plot_layout(guides="collect") &
  theme(strip.placement = "outside") 






data7.b2 <- data %>% 
  mutate(latitude=round(latitude), 
         rounding = factor("1deg (~100km)")) %>% 
  group_by(common_name, rounding, latitude) %>% 
  summarise(n = n(),
            arrival = unname(quantile(week26plus, probs = 0.025, type=1)), 
            departure = unname(quantile(week26plus, probs = 0.975, type=1))) %>% 
  pivot_longer(cols = arrival:departure,
               names_to = "date_type", values_to = "weeks26plus") %>% 
  mutate(date_type = factor(date_type, levels = c("departure","arrival")))

# or it might be better to find middle ground and use 75% observations
data7.b3 <- data %>% 
  mutate(latitude=round(latitude), 
         rounding = factor("1deg (~100km)")) %>% 
  group_by(common_name, rounding, latitude) %>% 
  summarise(n = n(),
            arrival = unname(quantile(week26plus, probs = 0.125, type=1)), 
            departure = unname(quantile(week26plus, probs = 0.875, type=1))) %>% 
  pivot_longer(cols = arrival:departure,
               names_to = "date_type", values_to = "weeks26plus") %>% 
  mutate(date_type = factor(date_type, levels = c("departure","arrival")))

# however, judging from the distribution of points in Plot2 (original) and from the 
# values in column "n", I still feel it is better to go with 95%.


data8.b2 <- data %>% 
  mutate(latitude=round(latitude,1), 
         rounding = factor("0.1deg (~10km)")) %>% 
  group_by(common_name, rounding, latitude) %>% 
  summarise(n = n(),
            arrival = unname(quantile(week26plus, probs = 0.025, type=1)), 
            departure = unname(quantile(week26plus, probs = 0.975, type=1))) %>% 
  pivot_longer(cols = arrival:departure,
               names_to = "date_type", values_to = "weeks26plus") %>% 
  mutate(date_type = factor(date_type, levels = c("departure","arrival")))

data8.b3 <- data %>% 
  mutate(latitude=round(latitude,1), 
         rounding = factor("0.1deg (~10km)")) %>% 
  group_by(common_name, rounding, latitude) %>% 
  summarise(n = n(),
            arrival = unname(quantile(week26plus, probs = 0.125, type=1)), 
            departure = unname(quantile(week26plus, probs = 0.875, type=1))) %>% 
  pivot_longer(cols = arrival:departure,
               names_to = "date_type", values_to = "weeks26plus") %>% 
  mutate(date_type = factor(date_type, levels = c("departure","arrival")))

data9.b2 <- bind_rows(data8.b2, data7.b2)
data9.b3 <- bind_rows(data8.b3, data7.b3)

## ##


## visualisations ##

# Plot3.b: visualising dates at two latitude scales and across the three species
plot3.b2 <- ggplot(data9.b2) + 
  theme_bw() + theme(strip.placement = "outside") +
  geom_point(aes(x=weeks26plus, 
                 y=latitude,
                 colour=date_type),
             alpha=2/3) +
  scale_colour_manual(values=c("maroon","sky blue")) +
  geom_vline(aes(xintercept=26)) +
  geom_text(aes(x=26, y=10, label="Dec"), nudge_x=-2, angle=90) +
  geom_text(aes(x=26, y=10, label="Jan"), nudge_x=1, angle=90) +
  facet_grid(rows = vars(data9.b2$rounding), 
             cols = vars(data9.b2$common_name),
             switch = "y") +
  scale_y_continuous(breaks=(seq(0,40,5))) +
  ggtitle("Variation of first arrival and last departure dates (95%) with latitude")

plot3.b3 <- ggplot(data9.b3) + 
  theme_bw() + theme(strip.placement = "outside") +
  geom_point(aes(x=weeks26plus, 
                 y=latitude,
                 colour=date_type),
             alpha=2/3) +
  scale_colour_manual(values=c("maroon","sky blue")) +
  geom_vline(aes(xintercept=26)) +
  geom_text(aes(x=26, y=10, label="Dec"), nudge_x=-2, angle=90) +
  geom_text(aes(x=26, y=10, label="Jan"), nudge_x=1, angle=90) +
  facet_grid(rows = vars(data9.b3$rounding), 
             cols = vars(data9.b3$common_name),
             switch = "y") +
  scale_y_continuous(breaks=(seq(0,40,5))) +
  ggtitle("Variation of first arrival and last departure dates (75%) with latitude")

plot3.b <- plot3.b2 + plot3.b3 + plot_layout(guides="collect", nrow=2) 


# Plot4.b: visualising latitude and arrival/departure dates with smooth lines
plot4.b2 <- 
  ggplot(data7.b2, aes(latitude, weeks26plus, colour=common_name)) + 
  scale_x_continuous(breaks=(seq(0,52,2))) +
  scale_y_continuous(limits=(c(0,55)),
                     breaks=(seq(0,55,5))) +
  labs(title="1deg (~100km)") +
  scale_colour_manual(values=c("#FFCC66","#996633","#CCCC66"),
                      aesthetics = c("colour","fill"),
                      labels = c("BRWa","BBFl","GrWa"),
                      name = "Species Code") +
  geom_point() +
  geom_smooth(orientation="x", span=1, aes(fill=common_name)) +
  geom_hline(aes(yintercept=26)) +
  geom_text(aes(y=26, x=10, label="Jan"), colour="black", nudge_y=2) +
  geom_text(aes(y=26, x=10, label="Dec"), colour="black", nudge_y=-2) +
  facet_wrap(~data7.b2$date_type, dir="v", strip.position = "left") +
  theme_bw() + theme(plot.title = element_text(hjust=0.5),
                     axis.title = element_blank()) +
  ggplot(data8.b2, aes(latitude, weeks26plus, colour=common_name)) + 
  scale_x_continuous(breaks=(seq(0,52,2))) +
  scale_y_continuous(limits=(c(0,55)),
                     breaks=(seq(0,55,5))) +
  labs(title="0.1deg (~10km)") +
  scale_colour_manual(values=c("#FFCC66","#996633","#CCCC66"),
                      aesthetics = c("colour","fill"),
                      labels = c("BRWa","BBFl","GrWa"),
                      name = "Species Code") +
  geom_point() +
  geom_smooth(orientation="x", span=0.5, aes(fill=common_name)) +
  geom_hline(aes(yintercept=26)) +
  geom_text(aes(y=26, x=10, label="Jan"), colour="black", nudge_y=2) +
  geom_text(aes(y=26, x=10, label="Dec"), colour="black", nudge_y=-2) +
  facet_wrap(~data8.b2$date_type, dir="v", strip.position = "left") +
  theme_bw() + theme(plot.title = element_text(hjust=0.5),
                     strip.text = element_blank(),
                     axis.title = element_blank()) +
  plot_annotation(title = "95% observations",
                  caption = "latitude",
                  theme = theme(plot.title = element_text(size=16),
                                plot.caption = element_text(hjust=0.47, size=12))) +
  plot_layout(guides="collect") &
  theme(strip.placement = "outside") 


plot4.b3 <- 
  ggplot(data7.b3, aes(latitude, weeks26plus, colour=common_name)) + 
  scale_x_continuous(breaks=(seq(0,52,2))) +
  scale_y_continuous(limits=(c(0,55)),
                     breaks=(seq(0,55,5))) +
  labs(title="1deg (~100km)") +
  scale_colour_manual(values=c("#FFCC66","#996633","#CCCC66"),
                      aesthetics = c("colour","fill"),
                      labels = c("BRWa","BBFl","GrWa"),
                      name = "Species Code") +
  geom_point() +
  geom_smooth(orientation="x", span=1, aes(fill=common_name)) +
  geom_hline(aes(yintercept=26)) +
  geom_text(aes(y=26, x=10, label="Jan"), colour="black", nudge_y=2) +
  geom_text(aes(y=26, x=10, label="Dec"), colour="black", nudge_y=-2) +
  facet_wrap(~data7.b3$date_type, dir="v", strip.position = "left") +
  theme_bw() + theme(plot.title = element_text(hjust=0.5),
                     axis.title = element_blank()) +
  ggplot(data8.b3, aes(latitude, weeks26plus, colour=common_name)) + 
  scale_x_continuous(breaks=(seq(0,52,2))) +
  scale_y_continuous(limits=(c(0,55)),
                     breaks=(seq(0,55,5))) +
  labs(title="0.1deg (~10km)") +
  scale_colour_manual(values=c("#FFCC66","#996633","#CCCC66"),
                      aesthetics = c("colour","fill"),
                      labels = c("BRWa","BBFl","GrWa"),
                      name = "Species Code") +
  geom_point() +
  geom_smooth(orientation="x", span=0.5, aes(fill=common_name)) +
  geom_hline(aes(yintercept=26)) +
  geom_text(aes(y=26, x=10, label="Jan"), colour="black", nudge_y=2) +
  geom_text(aes(y=26, x=10, label="Dec"), colour="black", nudge_y=-2) +
  facet_wrap(~data8.b3$date_type, dir="v", strip.position = "left") +
  theme_bw() + theme(plot.title = element_text(hjust=0.5),
                     strip.text = element_blank(),
                     axis.title = element_blank()) +
  plot_annotation(title = "75% observations",
                  caption = "latitude",
                  theme = theme(plot.title = element_text(size=16),
                                plot.caption = element_text(hjust=0.47, size=12))) +
  plot_layout(guides="collect") &
  theme(strip.placement = "outside") 


plot4.b <- wrap_elements(plot4.b2) + wrap_elements(plot4.b3) +
  plot_annotation(title = "Variation of first arrival and last departure dates with latitude",
                  theme = theme(plot.title = element_text(hjust=0.5, size=20))) +
  plot_layout(guides="collect", nrow=2) 

## ##



### 2: most eBirded areas should not dominate graphs 


## modifying data ##

# removing duplicate observations of same latitude and week26plus value
datax <- data %>% group_by(common_name, latitude, week26plus) %>% 
  summarise(n=n(), weeks=weeks)

# getting kernel density estimates
dens.1 <- kde2d(datax$weeks, datax$latitude)
dens.2 <- kde2d(datax$week26plus, datax$latitude)

# creating a new data frame of that 2d density grid
grid.1 <- data.frame(with(dens.1, expand.grid(x,y)), as.vector(dens.1$z))
names(grid.1) <- c("xgr", "ygr", "zgr")

grid.2 <- data.frame(with(dens.2, expand.grid(x,y)), as.vector(dens.2$z))
names(grid.2) <- c("xgr", "ygr", "zgr")

# fitting a model
mod.1 <- loess(zgr~xgr*ygr, data=grid.1, span=0.5)
mod.2 <- loess(zgr~xgr*ygr, data=grid.2, span=0.5)

# applying the model to the original data to estimate density at that point
datax$pointdens.1 <- predict(mod.1, 
                             newdata=data.frame(xgr=datax$weeks, 
                                                ygr=datax$latitude))
datax$pointdens.2 <- predict(mod.2, 
                             newdata=data.frame(xgr=datax$week26plus, 
                                                ygr=datax$latitude))

## ##


## visualisations ##

# Plot1.b: visualising change in reports with time and latitude
plot1.b <- ggplot(datax, aes(weeks, latitude, colour=pointdens.1)) + 
  scale_x_continuous(breaks=(seq(0,52,4))) +
  scale_y_continuous(breaks=(seq(0,40,5))) +
  coord_trans(y="log10") +
  scale_colour_gradient(low="dark blue", high="#6495ED") +
  labs(title="Reports of migrant species varying with timing and latitude") +
  geom_jitter(size=0.8) +
  geom_vline(xintercept=26) +
  annotate("text", x=24, y=10, label = "Jun", angle=90) +
  annotate("text", x=27, y=10, label = "Jul", angle=90) +
  facet_wrap(~common_name) + theme_bw()

# Plot2.b: visualising change in reports with winter-time and latitude
plot2.b  <-  ggplot(datax, aes(week26plus, latitude, colour=pointdens.2)) + 
  scale_x_continuous(breaks=(seq(0,52,4))) +
  scale_y_continuous(breaks=(seq(0,40,5))) +
  coord_trans(y="log10") +
  scale_colour_gradient(low="dark red", high="#FA8072") +
  labs(title="Reports of migrant species varying with winter-timing and latitude",
       subtitle="X-axis represents number of weeks after the 26th week of a year") +
  geom_jitter(size=0.8) +
  geom_vline(xintercept=26) +
  annotate("text", x=24, y=10, label = "Dec", angle=90) +
  annotate("text", x=27, y=10, label = "Jan", angle=90) +
  facet_wrap(~common_name) + theme_bw()

## ##

