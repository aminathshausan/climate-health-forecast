##  
# Author: Dr Aminath Shausan 

####################################
#this program  performs sensitivity analysis for seasonality
## time frame is months starting from  January 2008 to Dec 2022  
 
###############################

## load required data
source("./src/00_data_preprocess.R")

####################################
## Structure data to include additional variables
#####################################
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

## M3: model with seasonality
f.m2 <- obs ~ 1 +  
  f(inla.group(tas, n=20), model = "rw2", constr = TRUE,
    scale.model = TRUE,  hyper = precision.prior) + 
  f(inla.group(tasLag1, n=20), model = "rw2", constr = TRUE,
    scale.model = TRUE,  hyper = precision.prior)   +  
  f(month, replicate = ID.region, model = "rw1", cyclic = TRUE, constr = TRUE,
    scale.model = TRUE,  hyper = precision.prior) +
  f(ID.region, model = "iid", constr = TRUE,
    hyper = precision.prior)


#################################################################################
## ----Perform forecasts for various windows at lead times 1, 2, ..., 12
## select df with sliding window of 1 month  
## start first origin from dec 2017  to last origin Dec 2018 
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
# 
# origin <- 4
# fld <- origin
# nstep <-12
# rm(data.valid)
# data.valid <- getFldData(fld = origin, nstep = nstep) 

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
# use 1:13 origins 
for(i in 1:13) {
  
  rm(origin) 
  origin <- i ##change  
  print(sprintf("--- origin--: %d", origin))
  
  for(h in 1:1){
    set.seed(123+i+h)

    nstep = h ## change from 1 to  12 steps ahead
    print(sprintf("nstep: %d", nstep))

    rm(data.valid)
    data.valid <- getFldData(fld = origin, nstep = nstep) ### change this according to number of steps to predict
    
  # ### for sliding window, re-define time index
  # data.valid <- data.valid %>%
  #   group_by(state) %>%
  #   arrange(date) %>%
  #   mutate(ID.time = as.numeric(row_number())
  #   ) %>%
  #   ungroup() %>%
  #   as.data.frame()
  
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

write.csv(data.pred, file= "./results/sensitivity_covid/w10_m2_tas_tasLag1_nonlinr.csv", row.names = F)  ## change name as required


###################################################################################
## Fit optimal model (M3 = f.m2) to data containing optimal historic window 
## that is data from 2013 to 2022
###################################################################################


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
save(mdl, file = "./results/sensitivity_seasonality/optimal_model_fit.RData")

################################################################

