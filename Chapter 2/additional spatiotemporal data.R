##### Daymet data (Gridded)
# https://github.com/khufkens/daymetr

library(daymetr)
library(raster)

files<-list.files(tempdir(),full.names = T)
files
file.remove(files)

lat<-35.97358
lon<- -79.10037
var<-"prcp"
year<-2017

download_daymet_ncss(location = c(lat+0.1,lon-0.1,lat-0.1,lon+0.1),
                     start = year,
                     end = year,
                     param = var)

file<-list.files(tempdir(),full.names = T,pattern = paste(var,"_daily_",year,sep=""))
file
daymet_grid <- brick(file)
plot(daymet_grid)

##### PhenoCam data
library(phenocamapi)
library(tidyverse)

rois <- get_rois() %>% 
  filter(lat>=25 & lat<=50  & lon < -60 & lon> -130) %>%  # focus on CONUS
  filter(roitype=="DB") %>% # deciduous
  mutate(site_years=as.numeric(site_years)) %>% 
  filter(site_years>5)
rois

ggplot() +
  geom_polygon( data=map_data("state"), aes(x=long, y=lat, group = group),color="white", fill="grey95" ) +
  geom_point(data=rois,aes(x=lon,y=lat),cex=2,color="blue")+
  theme_void()


phenocam_ts<- vector(mode = "list", length = nrow(rois))
for (i in 1:nrow(rois)) {
  ts <- get_pheno_ts(site = rois$site[i], vegType = rois$roitype[i], roiID = rois$sequence_number[i], type = '3day') %>% 
    dplyr::select(date,midday_gcc) %>% 
    mutate(lon=rois$lon[i],lat=rois$lat[i],site = rois$site[i])
  phenocam_ts[[i]]<-ts
}
phenocam_ts_all<-rbindlist(phenocam_ts) %>% 
  mutate(date=as.Date(date))
  # tidyr::complete(date = seq(min(date), max(date), by = "day")) # fill in NAs for missing data
head(phenocam_ts_all)
