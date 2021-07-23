## Make static ggplot maps
library(data.table)
library(tigris)
library(tidycensus)
library(sf)
library(ggplot2)
library(grid)
library(gridExtra)
## Hardcoded map theme for all the ggplot stuff (I just fiddle with legend.position depending on map... I've found it really hard to automate this, if you want the legend inside the plot).
map_theme <- theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size=16),
        legend.text = element_text(size=5),
        # legend.justification=c(0,0),
        legend.position=c(0.9,0.1),
        # legend.position = 'bottom',
        legend.box.just	= 'right',
        legend.box = 'vertical',
        panel.grid.major = element_line(color=NA),
        panel.grid.minor = element_line(color=NA),
        # legend.key.width=unit(1.5,"line"),
        # legend.key.height=unit(1,"line"),
        legend.background = element_rect(fill='white',color='black'),
        legend.title = element_text(size=0),
        panel.background = element_rect(fill = "lightblue",colour = NA))
# plot.margin=unit(c(-.1,-.1,-.1,-.1), "cm"))
## For this map, everything is cropped relative to the bounding box of "phila_tracts" here - 
## but this could be any sf object.
# phila_fips <- '42101'
# phila_tracts <- tracts(state=42,county=101,class='sf')
# phila_tracts <- st_transform(phila_tracts, crs = 4326)
make_map <- function(phila_tracts,vars,clean_var_names=NULL,custom_breaks=NULL,scale_name='',show_median=T) {
## Get all background shapefiles from tigris
message('Getting counties...')
all_counties <- counties(class='sf')
all_counties <- st_transform(all_counties, crs = st_crs(phila_tracts))
all_counties <- st_make_valid(all_counties)
get_overlaps <- function(x,y) {
  overlaps <- st_intersects(x, y)
  overlaps <- unlist(lapply(1:dim(overlaps)[1], function(i) ifelse(length((overlaps[[i]]))==1, TRUE, FALSE)))
  return(x[overlaps,])
}
these_counties_sf <- st_crop(all_counties, phila_tracts)
these_counties <- as.data.table(these_counties_sf)
these_counties[, fips := paste0(STATEFP,COUNTYFP)]
these_states <- unique(these_counties[, STATEFP])
these_counties <- these_counties[, fips]
bad_fips <- c('26099','06079')
these_counties <- these_counties[!(these_counties %in% bad_fips)]
## Roads
# message('Getting roads...')
# get_county_roads <- function(c) {
#   roads <- roads(state=substr(c,1,2), county=substr(c,3,5), class='sf')
#   return(roads)
# }
# roads <- do.call(rbind, lapply(these_counties, get_county_roads))
# roads <- roads[!is.na(roads$RTTYP),]
# roads <- roads[roads$RTTYP!='C', ]
# roads <- roads[roads$RTTYP!='M', ]
# roads <- st_transform(roads, crs=st_crs(phila_tracts))
# roads <- st_crop(roads, st_bbox(phila_tracts))
## Parks
message('Getting parks...')
get_county_parks <- function(s) {
  parks <- landmarks(class='sf', state=s, type='area')
  return(parks)
}
parks <- do.call(rbind, lapply(these_states, get_county_parks))
parks <- parks[grepl('Park|Cmtry', parks$FULLNAME),]
parks <- st_transform(parks, crs=st_crs(phila_tracts))
parks <- st_crop(parks, phila_tracts)
## Water
message('Getting water...')
get_county_water <- function(c) {
  water <- area_water(state=substr(c,1,2), county=substr(c,3,5), class='sf')
  return(water)
}
water <- do.call(rbind, lapply(these_counties, get_county_water))
water <- st_transform(water, crs=st_crs(phila_tracts))
water <- st_crop(water, st_bbox(phila_tracts))
st_erase <- function(x, y) suppressWarnings(st_difference(st_buffer(st_make_valid(x),0), st_buffer(st_make_valid(st_union(st_combine(y))),0)))
these_counties_sf <- st_erase(these_counties_sf, water)

## Make map
message('Making maps...')
map_var_map <- function(var) {
  clean_name <- var
  if(!is.null(clean_var_names)) {
    clean_name <- clean_var_names[raw==var, clean]
  }
  if(is.null(custom_breaks)) {
    if(show_median) limits_to_use <- c(min(phila_tracts[[var]],na.rm=T),median(phila_tracts[[var]],na.rm=T),max(phila_tracts[[var]],na.rm=T))
    if(!show_median) limits_to_use <- c(min(phila_tracts[[var]],na.rm=T),max(phila_tracts[[var]],na.rm=T))
    if(max(limits_to_use)<max(phila_tracts[!is.na(phila_tracts[[var]]),][[var]])) phila_tracts[phila_tracts[[var]]>max(limits_to_use) & !is.na(phila_tracts[[var]]),][[var]] <- max(limits_to_use)
    if(min(limits_to_use)>min(phila_tracts[!is.na(phila_tracts[[var]]),][[var]])) phila_tracts[phila_tracts[[var]]<min(limits_to_use) & !is.na(phila_tracts[[var]]),][[var]] <- min(limits_to_use)
    if(!show_median) labels_to_use <- c(paste0('<=',round(min(limits_to_use),2)),paste0('>=',round(max(limits_to_use),2)))
    if(show_median) labels_to_use <- c(paste0('<=',round(min(limits_to_use),2)),round(median(phila_tracts[[var]],na.rm=T),2),paste0('>=',round(max(limits_to_use),2)))
  }
  if(!is.null(custom_breaks)) {
    limits_to_use <- custom_breaks[[var]]
    if(max(limits_to_use)<max(phila_tracts[!is.na(phila_tracts[[var]]),][[var]])) phila_tracts[phila_tracts[[var]]>max(limits_to_use) & !is.na(phila_tracts[[var]]),][[var]] <- max(limits_to_use)
    if(min(limits_to_use)>min(phila_tracts[!is.na(phila_tracts[[var]]),][[var]])) phila_tracts[phila_tracts[[var]]<min(limits_to_use) & !is.na(phila_tracts[[var]]),][[var]] <- min(limits_to_use)
    labels_to_use <- c(paste0('<=',round(min(limits_to_use),2)),paste0('>=',round(max(limits_to_use),2)))
  }
  if(length(unique(limits_to_use))==1) {
    limits_to_use <- unique(limits_to_use)
    labels_to_use <- limits_to_use
  }
  gg1 <- ggplot() + 
    geom_sf(data=these_counties_sf,
            alpha=1,
            fill='white',
            color=NA,
            lwd=0.5,
            inherit.aes = FALSE) + 
    geom_sf(data=phila_tracts,
            aes(fill=get(var)),
            alpha=0.7,
            color=NA,
            lwd=0,
            inherit.aes = FALSE) + 
    geom_sf(data=parks,
            alpha=1,
            fill='#e5f5e0',
            color=NA,
            lwd=0,
            inherit.aes = FALSE) + 
    #geom_sf(data=water,
    #        alpha=1,
    #        fill='#ccd4e3',
    #        color=NA,
    #        lwd=0,
    #        inherit.aes = FALSE) +
    # geom_sf(data=roads,
    #         alpha=1,
    #         color='black',
    #         lwd=0.1) + 
    ggtitle(clean_name) + 
    scale_fill_viridis_c(name=scale_name,
                         option='plasma',
                         limits=c(min(limits_to_use),max(limits_to_use)),
                         breaks=limits_to_use,
                         labels=labels_to_use) + 
    map_theme + 
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    # guides(fill = guide_legend(order = 1), 
    #        shape = guide_legend(order = 2)) + 
    # guides(colourbar = guide_legend(override.aes = list(size = 1))) 
    theme(legend.key.size = unit(0.015, "npc"))
    # theme(legend.box="vertical",legend.box.just = "right")
  return(gg1)
}
all_maps <- lapply(vars, map_var_map)
# png(paste0('C:/Users/ngraetz/Documents/repos/covid/eco_apartheid/results/',city,'.png'),height=11,width=11, units='in',res=600)
# grid.arrange(grobs=all_maps,ncol=3)
# dev.off()
return(all_maps)
}

