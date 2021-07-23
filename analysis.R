library(microACS)
library(cowplot)
library(tigris)
library(tidycensus)
library(data.table)
library(ggplot2)
library(raster)
library(sf)
library(flextable)
library(officer)
library(relaimpo)

repo <- 'C:/Users/ncgra/Dropbox/Penn/repos/eco_apartheid/'
inputs_dir <- paste0(repo,'inputs/')
outputs_dir <- paste0(repo,'outputs/')

####################################################################################################################
## 1. Based on top 20 most populated MSAs, grab national county/tract shapefiles for analysis.
####################################################################################################################
metro_tracts <- function(metro_name, st, cb) {
  message(metro_name)
  ## First, identify which counties intersect the metro area 
  metro <- filter(cb, grepl(metro_name, NAME))
  metro_zips <- zips[metro,]
  ctcodes <- ct[metro,][,c('COUNTYFP','STATEFP')]
  ## Grab all tracts in these counties
  tr <- Reduce(rbind, lapply(1:dim(ctcodes)[1], function(x) tracts(state=ctcodes[x,]$STATEFP, county=ctcodes[x,]$COUNTYFP, cb=TRUE, class='sf')))
  ## Now, find out which tracts are within the metro area
  within <- st_within(tr, metro)
  within_lgl <- map_lgl(within, function(x) {
    if (length(x) == 1) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  })
  output <- tr[within_lgl,]
  output$metro <- metro_name
  ## Return metro tracts with all info
  return(output)
}
## Get national shapefiles
st <- states(cb=T, class='sf')
ct <- counties(cb=T, class='sf')
zips <- zctas(class='sf')
cb <- core_based_statistical_areas(cb=T,class='sf')
cb_pops <- as.data.table(get_estimates(geography='cbsa',product='population'))
## Process top 20 metro areas
top20_cb <- cb_pops[order(-value)]
top20_cb <- gsub(' Metro Area','',top20_cb[1:20, NAME])
cb <- cb[cb$NAME %in% top20_cb,]
if(file.exists(paste0(inputs_dir,'top20_tracts.RDS'))) all <- readRDS(paste0(inputs_dir,'top20_tracts.RDS'))
if(!file.exists(paste0(inputs_dir,'top20_tracts.RDS'))) {
  all <- Reduce(rbind, lapply(top20_cb, metro_tracts, st, cb))
  saveRDS(all, paste0(inputs_dir,'top20_tracts.RDS'))
}
message(paste0('Total tracts: ', dim(all)[1]))
table(all$metro)

####################################################################################################################
## 2. Based on the county/tract shapefiles, process all national layers at the tract-level.
##    - Eco-apartheid
##        - Medicaid gap in enrollment (microACS)
##        - Heat vulnerability, electricity as % income / local MODIS heat (microACS)
##        - Cold vulnerability, gas as % of income (microACS)
##        - Proportion essential workers (ACS 5-year)
##        - Rent as percent of income (ACS 5-year)
##        - Police killings, spatially weighted for spillover effects? (Mapping Police Violence)
##    - Life expectancy
####################################################################################################################

####################################################################################################################
## Merge microACS variables with MODIS
####################################################################################################################

all$county_fips <- paste0(all$STATEFP, all$COUNTYFP)
all_counties <- unique(all$county_fips)
all_counties <- all_counties[all_counties!='06037']
all_states <- unique(all$STATEFP)
all_states <- st[st$STATEFP %in% unique(all$STATEFP), ]$STUSPS

## Load 2017 MODIS and create indoor heat island variable
for(month in c('06','07','08')) {
  message(paste0('Extracting ', month))
  modis <- raster(paste0(inputs_dir,'lst_day_v6_mean_max_1m_2017_',month,'_00.tif'))
  modis <- crop(modis, all)
  vals <- extract(modis, all)
  tempc <- lapply(vals, mean, na.rm=T)
  tempc <- unlist(tempc)
  all[[paste0('temp',month)]] <- tempc * (9/5) + 32
}
all$heat <- (all$temp06 + all$temp07 + all$temp08) / 3

## Change back to data.table for everything else to be fast
all.dt <- as.data.table(all)
for(fun_name in c('Medicaid enrollment gap','Electricity burden','Gas burden')) {
  sae <- readRDS(paste0(inputs_dir,fun_name,'.RDS'))
  all.dt <- merge(all.dt,sae,by='GEOID',all.x=T)
}
all.dt[, heat_decile := cut(heat, quantile(heat, probs=seq(0,1,0.1), na.rm=T), labels=1:10, include.lowest=T)]
all.dt[, elec_decile := cut(`Electricity burden`, quantile(`Electricity burden`, probs=seq(0,1,0.1), na.rm=T), labels=rev(1:10), include.lowest=T)]
all.dt[, heat_vulnerability := as.numeric(as.character(heat_decile)) / as.numeric(as.character(elec_decile))]
all.dt[, heat_decile := as.numeric(as.character(heat_decile))]
all.dt[, elec_decile := as.numeric(as.character(elec_decile))]

####################################################################################################################
## ACS variables from tidycensus (saved to repo, but code below to regenerate if necessary)
####################################################################################################################

# acs_cb <- data.table(name=c('Percent non-white','Percent black','Essential worker rate','Housing cost'),
#                      table=c('B02001','B02001','S2401','DP04'),
#                      var=c('B02001_002','B02001_003',NA,NA),
#                      denom=c('B02001_001','B02001_001',NA,'DP04_0136E'))
# 
# ## Rent burden = proportion of renters spending more than 30% monthly income on rent (DP04_0141E + DP04_0142E / DP04_0136E)
# ## Essential workers = proportion of total labor force (S2401_C01_001) in essential occupations according to ACLU occupational categories
# ## Percent Black and percent non-white = proportion over total tract population  
# process_county_acs <- function(v,fips) {
#   ## Get ACS data from tidycensus
#   s <- substr(fips,1,2)
#   c <- substr(fips,3,5)
#   hhincome <- as.data.table(get_acs(geography='tract',year=2018,table=acs_cb[name==v, table],state=s,county=c,cache_table = T))
#   hhincome <- dcast(hhincome, GEOID ~ variable, value.var='estimate')
#   ## Process variables
#   if(v=='Percent non-white') hhincome[, estimate := 1 - (get(acs_cb[name==v, var]) / get(acs_cb[name==v, denom]))]
#   if(v=='Percent black') hhincome[, estimate := get(acs_cb[name==v, var]) / get(acs_cb[name==v, denom])]
#   if(v=='Essential worker rate') {
#     ## Occupation categories from ACLU
#     essential_vars <- c('S2401_C01_015','S2401_C01_031','S2401_C01_030','S2401_C01_032',
#                         'S2401_C01_036','S2401_C01_034','S2401_C01_035','S2401_C01_028',
#                         'S2401_C01_027','S2401_C01_024','S2401_C01_023','S2401_C01_019',
#                         'S2401_C01_025','S2401_C01_020')
#     hhincome[, estimate := rowSums(.SD), .SDcols = essential_vars]
#     hhincome[, estimate := estimate / S2401_C01_001]
#     hhincome[, total_pop := S2401_C01_001]
#     hhincome[estimate==0, estimate := NA]
#   }
#   if(v=='Housing cost') {
#     ## Numerator = (housing units with selected monthly owner costs w/mortage > 35% monthly income) + 
#     ##             (housing units with rent > 35% monthly income)
#     hhincome[, estimate := rowSums(.SD), .SDcols = c('DP04_0115','DP04_0142')]
#     ## Denominator = total housing units
#     hhincome[, estimate := estimate / DP04_0001] 
#   }
#   ## Save outcome variable and GEOID
#   hhincome <- hhincome[, c('GEOID','estimate')]
#   hhincome[is.nan(estimate), estimate := NA]
#   hhincome <- hhincome[!is.na(estimate),]
#   setnames(hhincome, 'estimate', v)
#   saveRDS(hhincome, paste0(inputs_dir,'acs/',v,'/',fips,'.RDS'))
#   return(NULL)
# }
# for(v in c('Percent non-white','Percent black','Essential worker rate','Housing cost')) {
#   dir.create(paste0(outputs_dir,'acs/',v,'/'))
#   for(f in unique(all.dt[, county_fips])) {
#     message(paste0(v,', ', f))
#     process_county_acs(v,fips=f)
#   }
# }
# ## Combine
# for(fun_name in c('Percent non-white','Percent black','Essential worker rate','Housing cost')) {
#   geo_list <- rbindlist(lapply(list.files(paste0(outputs_dir,'acs/',fun_name),full.names=T), readRDS))
#   saveRDS(geo_list, paste0(outputs_dir,fun_name,'.RDS'))
# }

####################################################################################################################
## Police killings 2013-2020 (saved to repo, but code below to regenerate if necessary)
####################################################################################################################
# police <- fread(paste0(inputs_dir,"MPVDatasetDownload.csv"))
# police[, outcome := 1]
# police <- police[, list(police_killings=sum(outcome)), by='Zipcode']
# ## Get zip code crosswalk
# cross <- fread(paste0(inputs_dir,'zip_tract_crosswalk.csv'))
# cross[, Zipcode := as.character(ZIP)]
# cross[nchar(Zipcode)==3, Zipcode := paste0('00',Zipcode)]
# cross[nchar(Zipcode)==4, Zipcode := paste0('0',Zipcode)]
# setnames(cross, 'TRACT', 'GEOID')
# cross[, GEOID := as.character(GEOID)]
# cross[nchar(GEOID)==10, GEOID := paste0('0',GEOID)]
# ## Merge to get police killings in any zip code overlapping tract
# all_tract_police <- merge(all.dt[,c('GEOID')], cross[,c('GEOID','Zipcode')], all.x=T)
# all_tract_police[, Zipcode := as.numeric(Zipcode)]
# all_tract_police <- merge(all_tract_police, police, by='Zipcode', all.x=T)
# all_tract_police[is.na(Zipcode), police_killings := 0]
# all_tract_police[is.na(police_killings), police_killings := 0]
# all_tract_police <- all_tract_police[, list(police_killings=sum(police_killings)), by='GEOID']
# saveRDS(all_tract_police, paste0(inputs_dir,'police_killings.RDS'))
## Make smoother version (sum of all nearby tracts)
all_tract_police <- readRDS(paste0(inputs_dir,'police_killings.RDS'))

###############################
## MERGE EVERYTHING
###############################
le <- fread(paste0(inputs_dir,"US_A.CSV"))
le[, GEOID := as.character(`Tract ID`)]
le[nchar(GEOID)==10, GEOID := paste0('0',GEOID)]
## Life expectancy
final <- merge(all.dt, le[,c('GEOID','e(0)')], by='GEOID', all.x=T)
## Eco-apartheid variables
for(ind in c('Percent non-white','Percent black','Essential worker rate','Housing cost','police_killings')) {
  d <- readRDS(paste0(inputs_dir,ind,'.RDS'))
  final <- merge(final, d, by='GEOID', all.x=T)
}
saveRDS(final, paste0(outputs_dir,'final_eco_data.RDS'))

####################################################################################################################
## 3. Combine all layers into "eco-apartheid index" (PCA, first component?)
####################################################################################################################

## Merge EJSCREEN (too big for Git, download here: https://gaftp.epa.gov/EJSCREEN/2018/)
final <- readRDS(paste0(outputs_dir,'final_eco_data.RDS'))
ej <- fread(paste0(inputs_dir,"EJSCREEN_2019_StatePctiles.csv"))
ej[, GEOID := as.character(ID)]
ej[nchar(GEOID)==11, GEOID := paste0('0',GEOID)]
ej[, GEOID := substr(GEOID,1,11)]
ej[, ej := as.numeric(CANCER)] ## P_CANCR_D2 = EJ index, CANCER = raw value from EPA NATA
ej[, ej_cancer := as.numeric(P_CANCR_D2)]
ej <- ej[, list(ej=mean(ej), ej_cancer=mean(ej_cancer)), by='GEOID']
final <- merge(final, ej, by='GEOID', all.x=T)

## PCA
pca_vars <- c('Essential worker rate','Housing cost','police_killings','Medicaid enrollment gap', 'heat_vulnerability','ej')
sapply(final[,pca_vars,with=F], function(x) sum(is.na(x)))
for(v in pca_vars) final <- final[!is.na(get(v)),]
pca_dt <- copy(final[,pca_vars,with=F])
pca_dt <- scale(pca_dt)
pca <- prcomp(pca_dt)
final[, eco_apartheid := pca$x[,1] * -1]
final <- final[!is.na(`e(0)`), ]

## Tables
setwd(outputs_dir)
pca_table <- as.data.table(round(pca$rotation,2))
pca_table[, variable := rownames(pca$rotation)]
setcolorder(pca_table,'variable')
ft <- flextable(pca_table)
doc <- read_docx() %>%
  body_add_flextable(value = ft, split = TRUE) %>%
  body_end_section_landscape() %>% # a landscape section is ending here
  print( target = paste0("PCA_",Sys.Date(),".docx" ))
pca_table <- as.data.table(round(summary(pca)$importance,2))
pca_table[, variable := rownames(summary(pca)$importance)]
setcolorder(pca_table,'variable')
ft <- flextable(pca_table)
doc <- read_docx() %>%
  body_add_flextable(value = ft, split = TRUE) %>%
  body_end_section_landscape() %>% # a landscape section is ending here
  print( target = paste0("PCA2_",Sys.Date(),".docx" ))

## Merge SVI to compare
svi <- fread(paste0(inputs_dir,"SVI2018_US.csv"))
svi <- svi[, c('FIPS','RPL_THEMES')]
setnames(svi, c('GEOID','SVI'))
svi[, GEOID := as.character(GEOID)]
svi[nchar(GEOID)==10, GEOID := paste0('0',GEOID)]
final <- merge(final, svi, by='GEOID', all.x=T)
cor(final[, c('eco_apartheid','e(0)','SVI')])['eco_apartheid','e(0)']

## Merge redlining to compare
if(!file.exists(paste0(inputs_dir,'final_redlining_scores.RDS'))) {
  all <- readRDS(paste0(inputs_dir,'top20_tracts.RDS'))
  rl <- st_read(paste0(inputs_dir,'holc_ad_data.shp'))
  rl$outcome[rl$holc_grade=='A'] <- 1
  rl$outcome[rl$holc_grade=='B'] <- 2
  rl$outcome[rl$holc_grade=='C'] <- 3
  rl$outcome[rl$holc_grade=='D'] <- 4
  rl <- st_transform(rl, st_crs(all))
  rl <- st_make_valid(rl)
  get_rl_score <- function(zip) {
    message(zip)
    ## Crop RL boundaries to tract.
    t <- all[all$GEOID==zip,]
    rl_t <- st_intersection(rl, t)
    ## Calculate geographic-volume weights for each RL category.
    rl_weights <- data.table(rl=0:4, w=rep(0,5))
    for(r in unique(rl_t$outcome)) {
      rl_weights[rl==r, w := (sum(st_area(rl_t[rl_t$outcome==r,])) / st_area(t))]
    }
    rl_weights[rl==0, w := 1 - sum(rl_weights[, w])]
    ## Calculate contemporary tract-level RL score based on geographic-volume weighting of historic RL boundaries.
    rl_score <- rl_weights[, weighted.mean(rl,w)]
    return(rl_score)
  }
  redlining_scores <- unlist(lapply(all$GEOID, get_rl_score))
  all$redlining <- redlining_scores
  all_redlining <- as.data.table(all[c('GEOID','redlining')])
  all_redlining[, geometry := NULL]
}
if(file.exists(paste0(inputs_dir,'final_redlining_scores.RDS'))) all_redlining <- readRDS(paste0(inputs_dir,'final_redlining_scores.RDS'))
final <- merge(final, all_redlining, by='GEOID')

## Merge Area Deprivation Index (block group level, average to tract level)
adi <- fread(paste0(inputs_dir,"us_bg.txt"))
adi[, GEOID := as.character(as.numeric(FIPS))]
adi[nchar(GEOID)==11, GEOID := paste0('0',GEOID)]
adi[, GEOID := substr(GEOID,1,11)]
adi[, ADI_NATRANK := as.numeric(ADI_NATRANK)]
adi <- adi[, list(ADI=mean(ADI_NATRANK)), by='GEOID']
final <- merge(final, adi, by='GEOID', all.x=T)

## Save final merged dataset
saveRDS(final, paste0(outputs_dir,'final_eco_data_PCA_',Sys.Date(),'.RDS'))

####################################################################################################################
## 4. Create comparison figures to other index (bivariate correlations with life expectancy?):
##    - Social Vulnerability Index (CDC)
##    - Redlining (HOLC)
####################################################################################################################

## Scatter plots
final <- readRDS(paste0(outputs_dir,'final_eco_data_PCA_2021-02-01.RDS'))
pdf(paste0(outputs_dir,'/Figure1.pdf'),height=12,width=14)
gg1 <- ggplot(data=final,
       aes(x=eco_apartheid,
           y=`e(0)`,
           fill=`Percent non-white`)) + 
  geom_point(shape=21,color='black',size=5) + 
  geom_smooth(method='lm',color='black',se=F) + 
  scale_fill_viridis_c(name='% non-white') + 
  guides(alpha=FALSE) +
  labs(y='Life expectancy', x='Eco-apartheid') + 
  lims(x=c(-3,7.5)) + 
  theme_bw()
gg2 <- ggplot(data=final,
              aes(x=eco_apartheid,
                  y=`e(0)`,
                  fill=`Percent non-white`,
                  alpha=`Percent non-white`)) +
  geom_point(shape=21,color='black',size=5) +
  geom_smooth(method='lm',color='black',se=F) +
  scale_fill_viridis_c(name='% non-white') +
  guides(alpha=FALSE) +
  labs(y='Life expectancy', x='Eco-apartheid') +
  lims(x=c(-3,7.5)) +
  facet_wrap(~metro,ncol=4) +
  theme_bw()
grid.arrange(grobs=list(gg1,gg2 + theme(legend.position = 'none')),ncol=1)
dev.off()

## Tables
get_diff <- function(m) {
  m1 <- lm(`e(0)`~eco_apartheid,data=final)
  preds <- data.table(eco_apartheid=quantile(final[,eco_apartheid], probs=c(0.05,0.95)),
                      q=c(5,95),
                      metro=m)
  race_cor <- round(cor(final[,c('eco_apartheid','Percent non-white')])[1,2],2)
  if(m!='all') {
    m1 <- lm(`e(0)`~eco_apartheid,data=final[metro==m,])
    preds <- data.table(eco_apartheid=quantile(final[metro==m,eco_apartheid], probs=c(0.05,0.95)),
                        q=c(5,95),
                        metro=m)
    race_cor <- round(cor(final[metro==m,c('eco_apartheid','Percent non-white')])[1,2],2)
  }
  diffs <- rnorm(10000,predict(m1, newdata=preds)[1],predict(m1, newdata=preds, se=T)$se.fit[1]) - rnorm(10000,predict(m1, newdata=preds)[2],predict(m1, newdata=preds, se=T)$se.fit[2])
  preds[, pred := predict(m1, newdata=preds)]
  preds[, upper := pred + predict(m1, newdata=preds, se=T)$se.fit*1.96]
  preds[, lower := pred - predict(m1, newdata=preds, se=T)$se.fit*1.96]
  preds <- dcast(preds, metro ~ q, value.var=c('pred','upper','lower'))
  preds[, pred_diff := mean(diffs)]
  preds[, lower_diff := quantile(diffs,0.05)]
  preds[, upper_diff := quantile(diffs,0.95)]
  preds[, N := dim(final)[1]]
  if(m!='all') preds[, N := dim(final[metro==m,])[1]]
  preds[, race_cor := race_cor]
  return(preds)
}
metro_preds <- rbindlist(lapply(c('all',unique(final[,metro])), get_diff))
metro_preds[, pred_5 := paste0(round(pred_5,1), '(',round(lower_5,1),'-',round(upper_5,1),')')]
metro_preds[, pred_95 := paste0(round(pred_95,1), '(',round(lower_95,1),'-',round(upper_95,1),')')]
metro_preds[, diff := paste0(round(pred_diff,1), '(',round(lower_diff,1),'-',round(upper_diff,1),')')]
metro_preds <- metro_preds[order(-pred_diff)]
metro_preds <- metro_preds[,c('metro','pred_5','pred_95','diff','race_cor','N')]
setwd(outputs_dir)
ft <- flextable(metro_preds)
doc <- read_docx() %>%
  body_add_flextable(value = ft, split = TRUE) %>%
  body_end_section_landscape() %>% # a landscape section is ending here
  print(target = paste0("LE_DIFF_TABLE_",Sys.Date(),".docx" ))

## Correlation plots
get_cor <- function(m,v) {
  c <- cor(final[metro==m, c(v,'e(0)'), with=F], use='complete.obs')[v,'e(0)']
  return(data.table(metro=m,cor=c))
}
all_cors <- cor(final[, c('eco_apartheid','SVI','redlining','ej_cancer','ADI','e(0)'), with=F], use='complete.obs')
all_cors <- data.table(metro='All tracts',variable=rownames(all_cors),value=all_cors[,'e(0)'])
all_cors <- all_cors[variable!='e(0)',]
all_cors[variable=='eco_apartheid', variable := 'Eco-apartheid']
all_cors[variable=='redlining', variable := 'HOLC']
all_cors[variable=='ej_cancer', variable := 'EJ']
all_cors[, value := round(value,2)]
eco_cors <- rbindlist(lapply(unique(final[,metro]), get_cor, v='eco_apartheid'))
setnames(eco_cors, 'cor', 'eco_corr')
svi_cors <- rbindlist(lapply(unique(final[,metro]), get_cor, v='SVI'))
setnames(svi_cors, 'cor', 'svi_corr')
red_cors <- rbindlist(lapply(unique(final[,metro]), get_cor, v='redlining'))
setnames(red_cors, 'cor', 'red_corr')
ej_cors <- rbindlist(lapply(unique(final[,metro]), get_cor, v='ej_cancer'))
setnames(ej_cors, 'cor', 'ej_corr')
adi_cors <- rbindlist(lapply(unique(final[,metro]), get_cor, v='ADI'))
setnames(adi_cors, 'cor', 'adi_corr')
cors <- merge(eco_cors, svi_cors, by='metro')
cors <- merge(cors, red_cors, by='metro')
cors <- merge(cors, ej_cors, by='metro')
cors <- merge(cors, adi_cors, by='metro')
cors <- melt(cors, id.vars='metro')
metro_order <- eco_cors[order(eco_corr)]
cors[, metro := factor(metro, levels=metro_order$metro)]
cors[, value := round(value,2)]
cors[variable=='eco_corr', variable := 'Eco-apartheid']
cors[variable=='svi_corr', variable := 'SVI']
cors[variable=='adi_corr', variable := 'ADI']
cors[variable=='red_corr', variable := 'HOLC']
cors[variable=='ej_corr', variable := 'EJ']
les <- copy(final[, c('metro','e(0)')])
les[, le_range := round(max(`e(0)`)-min(`e(0)`)), by='metro']
cors <- merge(cors, unique(les[,c('metro','le_range')]), by='metro')
cors[, variable := factor(variable, levels=c('Eco-apartheid','SVI','ADI','EJ','HOLC'))]
cors[, metro := paste0(metro,' (',le_range,')')]
metro_order <- cors[variable=='Eco-apartheid',]
metro_order <- metro_order[order(value)]
cors <- rbind(cors,all_cors, fill=T)
cors[, value_string := as.character(value)]
cors[nchar(value_string)==4, value_string := paste0(value,'0')]
cors[, metro := factor(metro, levels=c('All tracts', metro_order$metro))]
cors[, labsize := ifelse(metro=='All tracts', 'big', 'small')]
pdf(paste0(outputs_dir,'Figure_3.pdf'),height=8,width=12)
ggplot() + 
  geom_tile(data=cors,
            aes(x=variable,y=metro,fill=value)) + 
  geom_label(data=cors,
             aes(x=variable,y=metro,label=value_string,size=labsize)) + 
  scale_fill_viridis_c(direction=-1,name='Correlation w/\ntract-level life\nexpectancy',option='plasma') + 
  scale_size_manual(values=c(6,4),guide=F) + 
  labs(x='',y='') + 
  theme_bw() + 
  theme(legend.position = 'bottom',
        legend.key.width = unit(2,'cm'),
        axis.text = element_text(size=14),
        axis.title = element_text(size=18))
dev.off()

####################################################################################################################
## 5. Compare to COVID-19 peak reported death rate by zip code
##    - Chicago: weekly deaths beginning March 1
##      - https://data.cityofchicago.org/Health-Human-Services/COVID-19-Cases-Tests-and-Deaths-by-ZIP-Code/yhhz-zm2v/data
##    - Philadelphia: weekly deaths beginning June 4
##      - https://github.com/ambientpointcorp/covid19-philadelphia/tree/master/deaths_by_zipcode
##    - NYC
##      - Cumulative deaths by zip code: https://github.com/nychealth/coronavirus-data/blob/master/totals/data-by-modzcta.csv
##      - Weekly deaths by borough: https://github.com/nychealth/coronavirus-data/blob/master/trends/data-by-day.csv
##    - LA
##      - Cumulative deaths by "city" (1-3 zip codes): http://dashboard.publichealth.lacounty.gov/covid19_surveillance_dashboard/ 
####################################################################################################################

final <- readRDS(paste0(outputs_dir,'final_eco_data_PCA_2021-02-01.RDS'))
final[, geometry := NULL]
all <- readRDS(paste0(inputs_dir,'top20_tracts.RDS'))
final_sf <- merge(all,final,by=c('GEOID','metro'))
get_rl_score <- function(zip, covid, geo_var, eco_sf) {
  ## Crop RL boundaries to tract.
  message(zip)
  t <- covid[covid[[geo_var]]==zip,]
  t <- st_make_valid(t)
  rl_t <- st_intersection(eco_sf, t)
  ## Calculate geographic-volume weights for each tract.
  rl_weights <- data.table(rl=rl_t$GEOID,estimate=rl_t$eco_apartheid,w=0)
  for(r in unique(rl_t$GEOID)) {
    rl_weights[rl==r, w := as.numeric(sum(st_area(rl_t[rl_t$GEOID==r,])) / st_area(t))]
  }
  elec_mean <- rl_weights[, weighted.mean(estimate,w)]
  if(is.nan(elec_mean)) elec_mean <- 0
  message(paste0('   ',sum(rl_weights$w)))
  return(elec_mean)
}
## LA
la_cities <- st_read(paste0(inputs_dir,"Countywide_Statistical_Areas__CSA_.shp"))
la_covid <- fread(paste0(inputs_dir,"LA_County_Covid19_CSA_case_death_table_Jan26.csv"))
la_covid[!(geo_merge %in% la_cities$LABEL),]
la_covid[, LABEL := geo_merge]
la_covid[, unadj_death_rate100k := (deaths_final/population)*100000]
la_covid <- la_covid[population>1000,] ## Median is 16k, there are some very sparse tracts with 0 COVID deaths.
la_covid <- merge(la_cities, la_covid, by='LABEL')
la_covid <- st_transform(la_covid, st_crs(final_sf))
la_sf <- final_sf[final_sf$metro=='Los Angeles-Long Beach-Anaheim, CA',]
la_sf <- st_make_valid(la_sf)
eco_cities <- unlist(lapply(la_covid$LABEL, get_rl_score, covid=la_covid, geo_var='LABEL', eco_sf=la_sf))
la_covid$eco_apartheid <- eco_cities
## NYC
nyc_covid <- fread(paste0(inputs_dir,"data-by-modzcta_jan26.csv"))
nyc_covid[, unadj_death_rate100k := COVID_DEATH_RATE]
all_zips <- zctas(state='NY',class='sf',year=2010)
all_zips$MODIFIED_ZCTA <- as.numeric(all_zips$ZCTA5CE10)
nyc_covid <- merge(all_zips,nyc_covid,by='MODIFIED_ZCTA')
nyc_covid <- st_transform(nyc_covid, st_crs(final_sf))
nyc_sf <- final_sf[final_sf$metro=='New York-Newark-Jersey City, NY-NJ-PA',]
nyc_sf <- st_make_valid(nyc_sf)
eco_cities <- unlist(lapply(nyc_covid$MODIFIED_ZCTA, get_rl_score, covid=nyc_covid, geo_var='MODIFIED_ZCTA', eco_sf=nyc_sf))
nyc_covid$eco_apartheid <- eco_cities
## Chicago
chicago_covid <- fread(paste0(inputs_dir,"COVID-19_Cases__Tests__and_Deaths_by_ZIP_Code_Jan26.csv"))
chicago_covid <- chicago_covid[`Week End`=='11/07/2020',]
chicago_covid[, unadj_death_rate100k := `Death Rate - Cumulative`]
chicago_covid[, zip := `ZIP Code`]
all_zips <- zctas(state='IL',class='sf',year=2010)
all_zips$zip <- as.numeric(all_zips$ZCTA5CE10)
chicago_covid <- merge(all_zips,chicago_covid,by='zip')
chicago_covid <- st_transform(chicago_covid, st_crs(final_sf))
chicago_sf <- final_sf[final_sf$metro=='Chicago-Naperville-Elgin, IL-IN-WI',]
chicago_sf <- st_make_valid(chicago_sf)
eco_cities <- unlist(lapply(chicago_covid$zip, get_rl_score, covid=chicago_covid, geo_var='zip', eco_sf=chicago_sf))
chicago_covid$eco_apartheid <- eco_cities
## Plot all cities COVID19 * eco-apartheid index
la_covid <- as.data.table(la_covid)
la_covid[, metro := 'Los Angeles-Long Beach-Anaheim, CA']
nyc_covid <- as.data.table(nyc_covid)
nyc_covid[, metro := 'New York-Newark-Jersey City, NY-NJ-PA']
chicago_covid <- as.data.table(chicago_covid)
chicago_covid[, metro := 'Chicago-Naperville-Elgin, IL-IN-WI']
all_covid <- rbindlist(list(la_covid[,c('metro','unadj_death_rate100k','eco_apartheid')],
                            chicago_covid[,c('metro','unadj_death_rate100k','eco_apartheid')],
                            nyc_covid[,c('metro','unadj_death_rate100k','eco_apartheid')]))
for(m in unique(all_covid[,metro])) {
  m1 <- lm(unadj_death_rate100k~eco_apartheid, data=all_covid[metro==m,])
  all_covid[metro==m, slope := round(m1$coefficients[2],2)]
  all_covid[metro==m, pred_covid := predict(m1)]
}
slopes <- unique(all_covid[,c('metro','slope','eco_apartheid','pred_covid')])
slopes <- slopes[, max := max(eco_apartheid), by='metro']
slopes <- slopes[eco_apartheid==max,]
pdf(paste0(outputs_dir,'Figure_4.pdf'),height=8,width=12)
ggplot() + 
  geom_smooth(data=all_covid,
              aes(x=eco_apartheid,
                  y=unadj_death_rate100k,
                  fill=metro,color=metro),
              method='lm',size=1) +
  geom_point(data=all_covid,
             aes(x=eco_apartheid,
                 y=unadj_death_rate100k,
                 fill=metro),
             shape=21,color='black',size=5) +
  scale_fill_viridis_d() + 
  scale_color_viridis_d() + 
  guides(color=F) + 
  ylim(c(0,500)) +
  labs(x='Eco-apartheid',y='Cumulative crude death rate per 100k',fill='Metro') + 
  theme_bw() +
  theme(legend.position = 'bottom',
        axis.text = element_text(size=14),
        axis.title = element_text(size=18))
dev.off()

## Tables
get_diff <- function(m) {
  m1 <- lm(unadj_death_rate100k~eco_apartheid,data=all_covid[metro==m,])
  preds <- data.table(eco_apartheid=quantile(all_covid[metro==m,eco_apartheid], probs=c(0.1,0.9)),
                      q=c(5,95),
                      metro=m)
  preds[, pred := predict(m1, newdata=preds, se=T)$fit]
  preds[, upper := pred + predict(m1, newdata=preds, se=T)$se.fit*1.96]
  preds[, lower := pred - predict(m1, newdata=preds, se=T)$se.fit*1.96]
  return(preds)
}
metro_preds <- rbindlist(lapply(unique(all_covid[,metro]), get_diff))
metro_preds <- dcast(metro_preds, metro ~ q, value.var=c('pred','upper','lower'))
metro_preds[, diff := pred_95-pred_5]
metro_preds[order(diff)]
metro_preds[,paste0(round(pred_5,1), '(',round(lower_5,1),'-',round(upper_5,1),')')]
metro_preds[,paste0(round(pred_95,1), '(',round(lower_95,1),'-',round(upper_95,1),')')]

get_diff <- function(m) {
  m1 <- lm(unadj_death_rate100k~eco_apartheid,data=all_covid)
  preds <- data.table(eco_apartheid=quantile(all_covid[,eco_apartheid], probs=c(0.05,0.95)),
                      q=c(5,95),
                      metro=m)
  if(m!='all') {
    m1 <- lm(unadj_death_rate100k~eco_apartheid,data=all_covid[metro==m,])
    preds <- data.table(eco_apartheid=quantile(all_covid[metro==m,eco_apartheid], probs=c(0.05,0.95)),
                        q=c(5,95),
                        metro=m)
  }
  diffs <- rnorm(10000,predict(m1, newdata=preds)[1],predict(m1, newdata=preds, se=T)$se.fit[1]) - rnorm(10000,predict(m1, newdata=preds)[2],predict(m1, newdata=preds, se=T)$se.fit[2])
  preds[, pred := predict(m1, newdata=preds)]
  preds[, upper := pred + predict(m1, newdata=preds, se=T)$se.fit*1.96]
  preds[, lower := pred - predict(m1, newdata=preds, se=T)$se.fit*1.96]
  preds <- dcast(preds, metro ~ q, value.var=c('pred','upper','lower'))
  preds[, pred_diff := mean(diffs)]
  preds[, lower_diff := quantile(diffs,0.05)]
  preds[, upper_diff := quantile(diffs,0.95)]
  preds[, N := dim(all_covid)[1]]
  if(m!='all') preds[, N := dim(all_covid[metro==m,])[1]]
  preds[, slope := round(coef(m1)[['eco_apartheid']],1)]
  return(preds)
}
metro_preds <- rbindlist(lapply(c(unique(all_covid[,metro])), get_diff))
metro_preds[, pred_5 := paste0(round(pred_5,1), '(',round(lower_5,1),'-',round(upper_5,1),')')]
metro_preds[, pred_95 := paste0(round(pred_95,1), '(',round(lower_95,1),'-',round(upper_95,1),')')]
metro_preds[, diff := paste0(round(pred_diff,1), '(',round(lower_diff,1),'-',round(upper_diff,1),')')]
metro_preds <- metro_preds[order(-pred_diff)]
metro_preds <- metro_preds[,c('metro','slope','pred_5','pred_95','diff','N')]
setwd(outputs_dir)
ft <- flextable(metro_preds)
doc <- read_docx() %>%
  body_add_flextable(value = ft, split = TRUE) %>%
  body_end_section_landscape() %>% # a landscape section is ending here
  print( target = paste0("COVID_DIFF_TABLE_",Sys.Date(),".docx" ))

####################################################################################################################
## 6. Extra figures
####################################################################################################################

## CITY-SPECIFIC MAPS
all <- readRDS(paste0(inputs_dir,'top20_tracts.RDS'))
final <- readRDS(paste0(outputs_dir,'final_eco_data_PCA_2021-02-01.RDS'))
final_sf <- merge(all, final[,c('GEOID','eco_apartheid')], by='GEOID')
source(paste0(repo,'map_template.R'))
## Philly: 42101
## St Louis: 29189, 29510
## Detroit: 26163
get_map <- function(m) {
  message(m)
  if(m=='Philadelphia') target_counties <- '^42101'
  if(m=='St Louis') target_counties <- '^29189|^29510'
  if(m=='Detroit') target_counties <- '^26163'
  map <- make_map(final_sf[grepl(target_counties,final_sf$GEOID),], vars='eco_apartheid', scale_name = 'Eco-apartheid', custom_breaks = list('eco_apartheid'=c(-2.5,3.5)))
  scat <- ggplot(data=final[grepl(target_counties,final_sf$GEOID),],
                 aes(x=eco_apartheid,
                     y=`e(0)`,
                     fill=`Percent non-white`)) + 
    geom_point(shape=21,color='black',size=5) + 
    geom_smooth(method='lm',color='black',se=F) + 
    scale_fill_viridis_c(name='% non-white',guide=F) + 
    guides(alpha=FALSE) +
    labs(y='Life expectancy', x='Eco-apartheid') + 
    lims(y=c(60,90), x=c(-2.5,3.2)) + 
    theme_bw() + 
    theme(strip.text = element_text(size=10),
          legend.position = 'bottom',
          legend.key.height = unit(0.5,'in'),
          legend.key.width = unit(1,'in'))
  return(list(map[[1]]+ggtitle(NULL),scat))
}
all_maps <- lapply(c('Philadelphia','St Louis','Detroit'),get_map)
gLegend<-function(a.plot){
  if ("ggplot" %in% class(a.plot)) {
    tmp <- ggplot_gtable(ggplot_build(a.plot))
  } else if ("grob" %in% class(a.plot)) {
    tmp <- .gplot
  }
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
scat <- ggplot(data=final[grepl('^42101',final_sf$GEOID),],
               aes(x=eco_apartheid,
                   y=`e(0)`,
                   fill=`Percent non-white`*100)) + 
  geom_point(shape=21,color='black',size=5) + 
  geom_smooth(method='lm',color='black',se=F) + 
  scale_fill_viridis_c(name='% non-white') + 
  guides(alpha=FALSE) +
  labs(y='Life expectancy', x='Eco-apartheid') + 
  lims(y=c(60,90), x=c(-2.5,3)) + 
  theme_bw() + theme(legend.position = 'bottom', legend.key.width = unit(0.9,'in'))
leg <- gLegend(scat)
scat <- ggplot(data=final[grepl('^42101',final_sf$GEOID),],
               aes(x=eco_apartheid,
                   y=`e(0)`,
                   fill=eco_apartheid)) + 
  geom_point(shape=21,color='black',size=5) + 
  geom_smooth(method='lm',color='black',se=F) + 
  scale_fill_viridis_c(name='Eco-apartheid',option='plasma') + 
  guides(alpha=FALSE) +
  labs(y='Life expectancy', x='Eco-apartheid') + 
  lims(y=c(60,90), x=c(-2.5,3)) + 
  theme_bw() + theme(legend.position = 'bottom', legend.key.width = unit(0.9,'in'))
map_leg <- gLegend(scat)
pdf(paste0(outputs_dir,'Figure_2.pdf'),height=14,width=9)
grid.arrange(grobs=list(all_maps[[1]][[1]]+ggtitle('Philadelphia')+theme(legend.position = 'none',plot.title = element_text(hjust = 0.5)),all_maps[[1]][[2]],
                        all_maps[[2]][[1]]+ggtitle('St. Louis')+theme(legend.position = 'none',plot.title = element_text(hjust = 0.5)),all_maps[[2]][[2]],
                        all_maps[[3]][[1]]+ggtitle('Detroit')+theme(legend.position = 'none',plot.title = element_text(hjust = 0.5)),all_maps[[3]][[2]],
                        leg,map_leg),
             layout_matrix=rbind(c(1,1,2,2),
                                 c(1,1,2,2),
                                 c(1,1,2,2),
                                 c(1,1,2,2),
                                 c(3,3,4,4),
                                 c(3,3,4,4),
                                 c(3,3,4,4),
                                 c(3,3,4,4),
                                 c(5,5,6,6),
                                 c(5,5,6,6),
                                 c(5,5,6,6),
                                 c(5,5,6,6),
                                 c(7,7,7,7),
                                 c(8,8,8,8)))
dev.off()
