##  
# Author: Dr Aminath Shausan 

####################################
#this program  performs sensitivity analysis after imputing missing weekly data
## time frame is months starting from  January 2008 to Dec 2022  
 
###############################


.libPaths("/Users/aminath/r-libraries")

options(digits=10)
rm(list = ls())
set.seed(963258)
#control.compute=list(save.memory=TRUE)

### load required libraries
# library(sf)
# library(sp)
library(dplyr)
library(ggplot2)
library(INLA)
########


## load required data
# source("./src/00_data_preprocess.R")

###########################
## Read and preprocess data
#####################################
## load the imputed combined influenza and climate data  
data<- read.csv('./data/pathogen/imputed_combinedInflClimate.csv') ## data for 2008 - 2022; 7 states
unique(data$year) ## covid period excluded
data$date <- as.Date(data$date, "%Y-%m-%d")

data.req <- data
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
                control.fixed = list(correlation.matrix = FALSE, #TRUE, 
                                     prec.intercept = 1, prec = 1),
                control.predictor = list(link = 1, compute = TRUE), 
                #num.threads = 8,
                verbose = FALSE)
  #model <- inla.rerun(model)
  return(model)
}


#############################################
## structure data for modelling
###########################################
### create lag 1 of temperature variable  
data.model <- data.req %>%
  group_by(state) %>%
  mutate(tasLag1 = dplyr::lag(tas, n =1, default = NULL)#,
  ) %>%
  ungroup()%>%
  as.data.frame()

#drop first 7 rows where lag1 values = NA
data.model <- data.model  %>% 
  filter(!is.na(tasLag1))

#add spatial, year, time ,  indexes
data.model <- data.model %>%
  mutate(
    ID.region = as.numeric(factor(state)),
    ID.year = as.numeric((factor(year)))
  )

data.model <- data.model %>%
  group_by(state) %>%
  arrange(date) %>%
  mutate(ID.time = as.numeric(row_number())) %>%
  ungroup() %>%
  as.data.frame()

## define model offset E = popn/100,000  
data.model$E <- data.model$popn/10^5
str(data.model)


################################################################
## Use the optimal model M3 (=f.m2)
########################################################

## rename required observed column
data.model <- rename(data.model, obs = observed_imp)

## M3: model with seasonality
f.m2 <- obs ~ 1 +  
  f(inla.group(tas, n=20), model = "rw2", constr = TRUE,
    scale.model = TRUE,  hyper = precision.prior) + 
  f(inla.group(tasLag1, n=20), model = "rw2", constr = TRUE,
    scale.model = TRUE,  hyper = precision.prior)   +  
  f(ID.time, replicate = ID.region, model = "ar1",  constr = TRUE,
    hyper = precision.prior) +
  f(ID.region, model = "iid", constr = TRUE,
    hyper = precision.prior)


#################################################################################
## ----Perform forecasts for various windows at lead times 1, 2, ..., 12
## select df with sliding window of 1 month  
## start first origin from dec 2017  to last origin Dec 2019 
## (this gives 13 forecasting origins) 
#################################################################################

## step 1: Define a function to get data for each lead time
getFldData = function(df, fld, nstep){

  t.indx = end.train.indx + nstep+(fld-1)  
  id.fit <- fld:t.indx
  id.frcs <- tail(id.fit, n=nstep)
  
  data.train <- df%>%
    filter(ID.time %in% id.fit)#%>%  
  
  data.train$obs[data.train$ID.time %in% id.frcs]<- NA 
  
  return(data.train)
}

## step 2: Define a function to get data for the required historic period  
getHistData = function(start_hist_year){
  
  data.window <- data.model %>%
    filter(year >= start_hist_year) ## 
  
  print(unique(data.window$year))
  print(length(unique(filter(data.window, year < 2018)$year)))  ## gives training window length

  ## step 2: re-assign temporal indices
  data.window <- data.window %>%
    group_by(state) %>%
    arrange(date) %>%
    mutate(ID.time = as.numeric(row_number()), 
          ID.year = as.numeric(row_number())) %>%
    ungroup() %>%
    as.data.frame()
  
  return(data.window)
}

## step 4: get required data for a given historic period
##change start_hist_year from 2008(10yrs of historic data), 2009, ..., 2015 (3 yrs of historic data)
rm(data.req2)
data.req2 <- getHistData(start_hist_year = 2008)

## step 3: get temporal index corresponding to Dec 2017 (the first forecast origin)
end.train.indx <- unique(filter(data.req2, year == 2017, month == 12)$ID.time)
print(end.train.indx)

## step 4: define a data frame to hold predictions at each origin
{
  rm(data.pred)
  data.pred <- data.req2[c("state","date","year","month", "obs", "popn", "ID.time")]
  data.pred$q2 <- NaN
} 

#######################
## Step 5:  fit model and save results for each forecast origin and a sliding window of 1 month
## do this step for each model separately
################
# use 1:25 origins 
for(i in 1:25) {
  
  rm(origin) 
  origin <- i ##change  
  print(sprintf("--- origin--: %d", origin))
  
  for(h in 1:12){ 
    set.seed(123+i+h)

    nstep = h ## change from 1 to  12 steps ahead
    print(sprintf("nstep: %d", nstep))

    rm(data.valid)
    data.valid <- getFldData(data.req2,fld = origin, nstep = nstep) ### change this according to number of steps to predict
    
  ### for sliding window, re-define time index
  data.valid <- data.valid %>%
    group_by(state) %>%
    arrange(date) %>%
    mutate(ID.time = as.numeric(row_number())
    ) %>%
    ungroup() %>%
    as.data.frame()
  
    rm(mdl)
    ## fit model
    mdl <- fit_model(formula = f.m2, data = data.valid, family = "poisson") ##"
    summary(mdl)

    t.indx = end.train.indx + nstep+(origin-1)
    print(sprintf("t_index: %d", t.indx))
    
    if (nstep ==1){
        data.pred$q2[data.pred$ID.time %in% origin: (t.indx)] <-  as.integer(mdl$summary.fitted.values[,"0.5quant"])
    }
    else {
      data.pred$q2[data.pred$ID.time == t.indx] <-  as.integer(tail(mdl$summary.fitted.values[,"0.5quant"],7))
    }
  }
  
  
  {
    ### save prediction for each forecast origin
    colname = paste(c("o", i), collapse = "")
    names(data.pred)[names(data.pred) == 'q2'] <- colname
  }
  
}

## save result in the form: (Eg; w10_m0_noClimate) for windowsize_m2_tas_tasLag1_nonlinr
### use this descriptions:
##    w10_m2_tas_tasLag1_nonlinr

write.csv(data.pred, file= "./results/sensitivity_imputed/w10_m2_tas_tasLag1_nonlinr.csv", row.names = F)  ## change name as required


###################################################################################
## Fit optimal model (M3 = f.m2) to data containing optimal historic window 
## that is data from 2013 to 2022
###################################################################################

rm(data.req2)
rm(data.req2)
data.req2 <- getHistData(start_hist_year = 2013)
# print(unique(data.req2$year))

data.req2 <- data.req2 %>%
  group_by(state) %>%
  arrange(date) %>%
  mutate(ID.time = as.numeric(row_number()))%>%
  ungroup() %>%
  as.data.frame()

rm(mdl)
mdl <- fit_model(f.m2, data = data.req2, family = "poisson")
summary(mdl)
print(c('waic,dic', mdl$waic$waic, mdl$dic$dic))

## save optimal model fit
save(mdl, file = "./results/model_fits/optimal_model_fit.RData")

################################################################

