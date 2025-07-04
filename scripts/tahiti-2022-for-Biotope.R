#------------------------------------------------------------------------------#
#### Symphonie des Iles ####
        
# cahiers de la fondation Biotope
# analyses des indices acoustiques
# Jean-Yves Barnagaud (jean-yves.barnagaud@ephe.psl.eu)
#------------------------------------------------------------------------------#

library(lubridate)
library(hms)
library(suncalc)
library(sf)
library(ggplot2)
library(ade4)
library(factoextra)
library(RColorBrewer)
library(patchwork)

## données ---------------------------------------------------------------------

# position des enregistreurs
xy <- read.csv2("data/spatial/positions_SM.csv")

wd <- "data/acoustic-indices/"
df.all <- data.frame()

# données des SM

for(k in c("Moorea","Tahiti_Papehue","Tahiti_controle")){

filename <- paste(wd,"indices_acou_",k,"__",sep="")

if(k =="Moorea"){
  vec <- c(1,2,3,5)
} else {
  vec <- 1:5
}

df.tp <- data.frame()

for (i in vec) {
  
# charger les données 
  
  tp0 <- read.csv2(
    paste(
      filename,
      i,
      ".csv",
      sep = ""
    )
  )

# START de chaque fichier (1 fichier = 1 période de 30')
  
  yr <- substring(tp0$FILE,
                  first = nchar(tp0$FILE) - 18,
                  last = nchar(tp0$FILE) - 15)
  mth <- substring(tp0$FILE,
                   first = nchar(tp0$FILE) - 14,
                   last = nchar(tp0$FILE) - 13)
  dy <- substring(tp0$FILE,
                  first = nchar(tp0$FILE) - 12,
                  last = nchar(tp0$FILE) - 11)
  
  hr <- substring(tp0$FILE,
                  first = nchar(tp0$FILE) - 9,
                  last = nchar(tp0$FILE) - 8)
  mn <- substring(tp0$FILE,
                  first = nchar(tp0$FILE) - 7,
                  last = nchar(tp0$FILE) - 6)
  sc <- substring(tp0$FILE,
                  first = nchar(tp0$FILE) - 5,
                  last = nchar(tp0$FILE) - 4)
  start <- ymd_hms(paste(yr, mth, dy, hr, mn, sc, sep = "-"), tz = "Pacific/Tahiti")
  
  tp0$START_30mn <- ymd_hms(start, tz = "Pacific/Tahiti")

# agréger les indices par périodes de 30' (devient la résolution temporelle)

tp30 <- aggregate(tp0[,1:11],by=list(tp0$FILE,tp0$START_30mn,tp0$SM),FUN="median")
colnames(tp30)[1:3] = c("FILE","START_30MN","SM")
tp30$SITE <- paste(k, i, sep = "_")

# supprimer 1er / dernier fichiers (annonces de début / fin)

tp30b <- subset(tp30,START_30MN!=max(tp30$START_30MN) & START_30MN!=min(tp30$START_30MN))
tp30b$AREA <- k 

df.tp <- rbind(df.tp,tp30b)

}

df.all <- rbind(df.all,df.tp)

}

df.all$AREA[which(df.all$AREA == "Moorea")] = "Moorea"
df.all$AREA[which(df.all$AREA == "Tahiti_Papehue")] = "Tahiti - Papehue"
df.all$AREA[which(df.all$AREA == "Tahiti_controle")] = "Tahiti - Arahurahu"
df.all$AREA <- factor(df.all$AREA)

# séparer jour / nuit

sun <- getSunlightTimes(
  date = as.Date("2022-04-05"),
  lat = -17.34,
  lon = -149.54,
  keep = c("sunrise", "sunset"),
  tz = "Pacific/Tahiti"
)

df.all$START_HOUR <- as_hms(df.all$START_30MN)


df.day <- subset(df.all,START_HOUR > as_hms(sun[1,"sunrise"]) & START_HOUR < as_hms(sun[1,"sunset"]))
df.night <- subset(df.all,START_HOUR < as_hms(sun[1,"sunrise"]) | START_HOUR > as_hms(sun[1,"sunset"]))

## 1 - caractérisation des paysages acoustiques --------------------------------

ind.sound <- c("BIOAC", "H", "ACI", "NDSI", "ADI", "NP")

# ACP sur soundscape jour
s.df.day <- df.day[, ind.sound]
s.day0 <- dudi.pca(s.df.day, scannf = F, nf = 5)
df.day$SITE <- factor(df.day$SITE)
s.day.bca <-bca(s.day0,fac = df.day$SITE,scannf = F, nf = 2)
site.day <- factor(c(rep("Moorea",4),rep("Tahiti - Papehue",5),rep("Tahiti - Arahurahu",5)))

s.class(s.day0$li, fac = df.day$AREA)
s.class(s.day.bca$li, fac = site.day)

# ACP sur soundscape nuit
s.df.night <- df.night[, ind.sound]
s.night0 <- dudi.pca(s.df.night, scannf = F, nf = 2)
df.night$SITE <- factor(df.night$SITE)
s.night.bca <-bca(s.night0,fac = df.night$SITE,scannf = F, nf = 2)
site.night <- factor(c(rep("Tahiti - Arahurahu",5),rep("Tahiti - Papehue",5)))

s.class(s.night0$li, fac = df.night$AREA)
s.class(s.night.bca$li, fac = site.night)

# figure
palette.poly.day <- c("#7e03a8", "#fdc527", "#e66c5c")
palette.poly.night <- c( "#fdc527", "#e66c5c")

p.day <- fviz_pca_biplot(
  s.day0,
  label = "var",
  habillage = df.day$AREA,
  axes = c(1, 2),
  addEllipses=TRUE,
  ellipse.level=0.95,
  palette = palette.poly.day, geom="point",
  alpha = 0.5,
  mean.point.size = 0,
  col.var = "black",
  ggtheme = theme_minimal(),
  title = "(A) Jour",
  legend.title = "Sites"
)
p.day

p.night <- fviz_pca_biplot(
  s.night0,
  label = "var",
  habillage = df.night$AREA,
  axes = c(1, 2),
  addEllipses=TRUE,
  ellipse.level=0.95,
  palette = palette.poly.night, geom="point",
  alpha = 0.5,
  mean.point.size = 0,
  col.var = "black",
  ggtheme = theme_minimal(),
  title = "(B) Nuit",
  legend.title = "Sites"
)
p.night

p.day / p.night
ggsave("outputs/ordination-sites.png",width = 7.5, height = 15)

## 2 - Analyse spatiale : cartographie -----------------------------------------

# la cartographie est faite sous qGIS par simplicité

avg.day <- aggregate(df.day[, c("BIOAC", "H", "ACI", "ADI", "NDSI", "NP")], by = list(df.day$SITE), FUN =
                       "mean")
cv.day <- aggregate(
  df.day[, c("BIOAC", "H", "ACI", "ADI", "NDSI", "NP")],
  by = list(df.day$SITE),
  FUN = function(x) {
    sd(x) / mean(x)
  }
)
colnames(cv.day)[-1] <- c("cv.BIOAC","cv.H","cv.ACI","cv.ADI","cv.NDSI","cv.NP")
synth.day <- merge(avg.day, cv.day, by = "Group.1")
colnames(synth.day)[1] <- "name"
xy.day <- merge(synth.day,xy,by = "name" , all = F)
  
avg.night <- aggregate(df.night[, c("BIOAC", "H", "ACI", "ADI", "NDSI", "NP")], by = list(df.night$SITE), FUN =
                         "mean")
cv.night <- aggregate(
  df.night[, c("BIOAC", "H", "ACI", "ADI", "NDSI", "NP")],
  by = list(df.night$SITE),
  FUN = function(x) {
    sd(x) / mean(x)
  }
)
colnames(cv.night)[-1] <- c("cv.BIOAC","cv.H","cv.ACI","cv.ADI","cv.NDSI","cv.NP")
synth.night <- merge(avg.night, cv.night, by = "Group.1")
colnames(synth.night)[1] <- "name"
xy.night <- merge(synth.night,xy,by = "name" , all = F)



## 3 - Séries temporelles ------------------------------------------------------

# Papehue 
st.tahiti.pa <- subset(df.all,AREA == "Tahiti - Papehue")

ndsi.t <- ggplot(st.tahiti.pa)+
  aes(x = START_30MN, y = NDSI)+
  geom_point()+
  geom_smooth(color =  "#e66c5c")+
  facet_wrap(~SITE)+
  theme_minimal()+
  labs(x = "Heure",y = "NDSI",title="")

bioac.t <- ggplot(st.tahiti.pa)+
  aes(x = START_30MN, y = BIOAC)+
  geom_point()+
  geom_smooth(color =  "steelblue")+
  facet_wrap(~SITE)+
  theme_minimal()+
  labs(x = "Heure",y = "BIOAC",title="")

# Arahurahu 
st.tahiti.ar <- subset(df.all,AREA == "Tahiti - Arahurahu")

ndsi.ar <- ggplot(st.tahiti.ar)+
  aes(x = START_30MN, y = NDSI)+
  geom_point()+
  geom_smooth(color =  "#e66c5c")+
  facet_wrap(~SITE)+
  theme_minimal()+
  labs(x = "Heure",y = "NDSI",title="")

bioac.ar <- ggplot(st.tahiti.ar)+
  aes(x = START_30MN, y = BIOAC)+
  geom_point()+
  geom_smooth(color =  "steelblue")+
  facet_wrap(~SITE)+
  theme_minimal()+
  labs(x = "Heure",y = "BIOAC",title="")

ndsi.t / bioac.t
ggsave("outputs/time_series_pa.png",width = 15, height = 7)

ndsi.ar / bioac.ar
ggsave("outputs/time_series_ar.png",width = 15, height = 7)



