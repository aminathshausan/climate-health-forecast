# Author: Dr Aminath Shausan 

####################################
#this program  pre-processes data
###############################

.libPaths(<PUT PATH TO DIRECTORY>)

options(digits=10)
rm(list = ls())
set.seed(963258)

### load required libraries
library(sf)
library(sp)
library(dplyr)
library(ggplot2)
library(INLA)
########

####################################
## load required cleaned data
#####################################
## read map of AU
map <- read_sf(dsn = './data/shapeFile', layer ='STE_2021_AUST_GDA2020')
#remove Christmas and Cocos islands
map <- map %>%
  filter(!(STE_NAME21 %in%  c('Australian Capital Territory','Other Territories',
                              'Outside Australia'))) %>%
  arrange(STE_NAME21)  
map <- map %>%
  select(STE_CODE21, STE_NAME21, geometry)

## load the combined influenza and climate data  
data<- read.csv('./data/pathogen/combinedInflClimate.csv') ## 
unique(data$year)
data$date <- as.Date(data$date, "%Y-%m-%d")
## drop covid period (2020 to 2021)  
data.req <- data%>%
  filter(!(year %in% c(2020, 2021))) #%>%
print(unique(data.req$year))
### order data according to date and state
data.req <- data.req[order(data.req$date,data.req$state),] 
rownames(data.req)<-1:nrow(data.req)

## map region acronym to their actual name
region_lookup <- c("NSW" = "New South Wales",  
                   "NT"  = "Northern Territory",  
                   "QLD" = "Queensland",  
                   "SA"  = "South Australia",  
                   "TAS" = "Tasmania",  
                   "VIC" = "Victoria",  
                   "WA"  = "Western Australia")
data.req$STE_NAME21<- region_lookup[data.req$state]


# add levels of humidity
summary(data.req$RH)
data.req$RH_cat <- cut(data.req$RH,
                       breaks = quantile(data.req$RH, probs = c(0, 1/3, 2/3, 1)),
                       labels = c("low","med","high"),
                       include.lowest = TRUE)

data.req$RH_cat <- factor(data.req$RH_cat, levels = c("low","med","high"))
str(data.req)

###############################################################################
##.  define priors : Penalized PC prior
############################################################################
precision.prior <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))

###############################################################################
##.  define model fit function 
############################################################################
fit_model <- function(formula, data = df, family = "poisson", config = FALSE)
  
{
  model <- inla(formula = formula, data = data, family = family, offset = log(E),
                control.inla = list(strategy = 'adaptive'), 
                control.compute = list(dic = TRUE, waic=TRUE, config = config, 
                                       cpo = TRUE, return.marginals = FALSE),
                control.fixed = list(correlation.matrix = FALSE, # 
                                     prec.intercept = 1, prec = 1),
                control.predictor = list(link = 1, compute = TRUE), 
                #num.threads = 8,
                verbose = FALSE)
  return(model)
}

