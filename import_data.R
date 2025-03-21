Cabauw_measure <- read.csv("../../Data/Data_Cabauw/Cabauw_measure.csv")
Cabauw_measure$Year = year(Cabauw_measure$Year)
Cabauw_measure$DateTime = ymd_hms(Cabauw_measure$DateTime)
Cabauw_measure = Cabauw_measure[Cabauw_measure$Year > 2000 & Cabauw_measure$Year < 2020,]


Cabauw_KNW <- read.csv("../../Data/Cabauw_KNW/modeldata.csv")
Cabauw_KNW$DateTime = ymd_hms(Cabauw_KNW$DateTime)
Cabauw_KNW$Year = year(Cabauw_KNW$Year)
Cabauw_KNW = Cabauw_KNW[Cabauw_KNW$Year > min(Cabauw_KNW$Year) & Cabauw_KNW$Year < max(Cabauw_KNW$Year),]


Cabauw_RACMO <- read.csv("../../Data/Data_RACMO/fulldata.csv")
Cabauw_RACMO$DateTime = ymd_hms(Cabauw_RACMO$DateTime)
Cabauw_RACMO = Cabauw_RACMO[Cabauw_RACMO$Year > min(Cabauw_RACMO$Year) & Cabauw_RACMO$Year < max(Cabauw_RACMO$Year),]


Schiphol_measure <- read.csv("../../Data/Schiphol_measure/dataSchiphol.csv")
Schiphol_measure$DateTime = ymd_hms(Schiphol_measure$DateTime)
Schiphol_measure = Schiphol_measure[Schiphol_measure$Year > min(Schiphol_measure$Year) & Schiphol_measure$Year < max(Schiphol_measure$Year),]


Schiphol_RACMO <- read.csv("../../Data/Data_RACMO/Schiphol/RACMO_Schiphol.csv")
Schiphol_RACMO$DateTime = ymd_hms(Schiphol_RACMO$DateTime)
Schiphol_RACMO = Schiphol_RACMO[Schiphol_RACMO$Year > min(Schiphol_RACMO$Year) & Schiphol_RACMO$Year < max(Schiphol_RACMO$Year),]




all_data = list(measure = Cabauw_measure,KNW = Cabauw_KNW, RACMO = Cabauw_RACMO)




all_data$RACMOh6 <- Cabauw_RACMO %>%
  mutate(
    DT = DateTime
  )%>%
  group_by(DateTime = date(DT) + hours(floor(hour(DT)/6)*6)) %>%
  summarise(
    Year = as.integer(mean(Year)),
    w10m = mean(w10m),
    F050 = mean(fh050),
    F100 = mean(fh100),
    F150 = mean(fh150),
    F200 = mean(fh200),
    F250 = mean(fh250),
    F300 = mean(fh300),
    w10max = max(w10max)
  )

all_data$RACMOh12 <- Cabauw_RACMO %>%
  mutate(
    DT = DateTime
  )%>%
  group_by(DateTime = date(DT) + hours(floor(hour(DT)/12)*12)) %>%
  summarise(
    Year = as.integer(mean(Year)),
    w10m = mean(w10m),
    F050 = mean(fh050),
    F100 = mean(fh100),
    F150 = mean(fh150),
    F200 = mean(fh200),
    F250 = mean(fh250),
    F300 = mean(fh300),
    w10max = max(w10max)
  )


all_data$RACMOh24 <- Cabauw_RACMO %>%
  mutate(
    DT = DateTime
  )%>%
  group_by(DateTime = date(DT)) %>%
  summarise(
    Year = as.integer(mean(Year)),
    w10m = mean(w10m),
    F050 = mean(fh050),
    F100 = mean(fh100),
    F150 = mean(fh150),
    F200 = mean(fh200),
    F250 = mean(fh250),
    F300 = mean(fh300),
    w10max = max(w10max)
  )
