##
# Author: Dr Aminath Shausan 

####################################
#this program visualises data and model fits
##################################

## load required data
source("./src/00_data_preprocess.R")
########################################
# plot influenza incidence rate (per 100,000 population) heat maps (month and year) per region 
########################################
iir_sptemp <- 
  data.req %>% 
  group_by(year, month, STE_NAME21) %>%
  # calculate region level incidence rate
  summarise(cases = sum(obs),
            popn = sum(popn)) %>% 
  mutate(iir = cases / popn * 10^5) %>% 
  ggplot(aes(x = month, y = factor(year), fill = iir)) + 
  geom_raster() +
  ylab('Year') + 
  xlab('Month') +
  scale_fill_distiller(  name    = 'Incidence rate',  
                         transform = "log1p",
                         breaks = c(0, 30, 100, 400, 900),
                         labels = c(0, 30, 100, 400, 900),
                         palette =  "YlGnBu",  direction = 1) + 
  scale_x_continuous(breaks = c(1,4,7,10), labels = c('Jan', 'Apr', 'Jul', 'Oct')) +
   scale_y_discrete(breaks = levels(factor(data.req$year))[seq(1, length(levels(factor(data.req$year))), by = 2)]  )+
   theme_bw(base_size = 12, base_family = "sans", header_family = "sans") +
   facet_wrap(~ STE_NAME21, nrow = 3, ncol = 3, scales = "free")

iir_sptemp
ggsave('./figs/iir_sptemp.pdf', height = 20, width = 30, units = "cm")
  

############################################################
# plot  temperature (tas) heat map (month and year) per region 
#########################################################
tas_sptemp <- 
  data.req %>% 
  group_by(year, month, STE_NAME21) %>%
  summarise(tas = mean(tas), .groups = "drop") %>%
  ggplot(aes(x = month, y = factor(year), fill = tas)) + 
  geom_raster() +
  xlab("Month") +
  ylab("Year") +
  scale_fill_distiller(  name    = "Temperature (°C)",  palette =  "RdBu",  direction = -1)+ ## YlGn
  scale_x_continuous(
    breaks = c(1, 4, 7, 10),
    labels = c("Jan", "Apr", "Jul", "Oct")
  ) +
  scale_y_discrete(breaks = levels(factor(data.req$year))[seq(1, length(levels(factor(data.req$year))), by = 2)]  )+
  theme_bw(base_size = 12) +
  facet_wrap(~ STE_NAME21, nrow = 3, ncol = 3, scales = "free")

tas_sptemp
ggsave('./figs/tas_sptemp.pdf', height = 20, width = 30, units = "cm")

###################################################
# plot  humidity (RH) heat map (month and year) per region 
#################################################
rh_sptemp <- 
  data.req %>% 
  group_by(year, month, STE_NAME21) %>%
  summarise(tas = mean(RH), .groups = "drop") %>%
  ggplot(aes(x = month, y = factor(year), fill = tas)) + 
  geom_raster() +
  xlab("Month") +
  ylab("Year") +
  scale_fill_distiller(  name    = "Humidity",  palette =  "BrBG",  direction = -1)+ ##RdBu BrBG
  scale_x_continuous(
    breaks = c(1, 4, 7, 10),
    labels = c("Jan", "Apr", "Jul", "Oct")
  ) +
  scale_y_discrete(breaks = levels(factor(data.req$year))[seq(1, length(levels(factor(data.req$year))), by = 2)]  )+
  theme_bw(base_size = 12) +
  facet_wrap(~ STE_NAME21, nrow = 3, ncol = 3, scales = "free")
rh_sptemp
ggsave('./figs/RH_sptemp.pdf', height = 20, width = 30, units = "cm")



################################################################################
## plot results from optimal model (M3) ##
##############################################################################
## filter data from 2013 onward (this contains optimal model fit data)
## load optimal model fit (model 3 = f.m2)
load('./results/optimal_model_fit.RData')

rm(data.req2)
data.req2 <- data.req %>%
  filter(year >= 2013) ## change year from 2008(10yrs of historic data) to 20017(1yrs of historic)
print(unique(data.req2$year))

data.req2 <- data.req2 %>%
  group_by(state) %>%
  arrange(date) %>%
  mutate(ID.time = as.numeric(row_number()), 
         ID.year = as.numeric(year)) %>%
  ungroup() %>%
  as.data.frame()

### ------ plot a maps of spatial random effects  
# extract posterior mean estimates for the spatial random effect
print(names(mdl$summary.random))
post_mean_region <- mdl$summary.random$ID.region$mean
map$post_mean <- post_mean_region
mn <-min(map$post_mean)
mx <-max(map$post_mean)

space_effects <- ggplot() + 
  geom_sf(data = map, aes(fill = post_mean), lwd = 0, color = NA) +
  scale_fill_distiller(palette = "PRGn", direction = -1, 
                       limits = c(min(mn,-mx),max(mx,-mn))) +
  labs(fill = "Contribution to \n log(Incidence rate)") +
  geom_sf_text(data = map, aes(label = STE_NAME21), size = 4, color = "black", family = "sans")+
  theme_void(base_size = 12, base_family = "sans", header_family = "sans") 

space_effects
ggsave('./figs/spatial_effect.pdf', height = 20, width = 30, units = "cm")

### ----- plot long term temporal random effects  
length(unique(mdl$summary.random$ID.time$ID)) 
nrow(mdl$summary.random$ID.time)
data.plt <- data.req2 %>%
  mutate(q1 = mdl$summary.random$ID.time[,"0.025quant"],
         q2 = mdl$summary.random$ID.time[,"0.5quant"],
         q3 =mdl$summary.random$ID.time[,"0.975quant"],
  )


x_breaks  <- c(1, 13, 25, 37, 49, 61, 73, 85) ## , 96
x_labels  <- c("2013", "2014", "2015", "2016", "2017", "2018", "2019", "2022")
temporal_effects <- 
  # data.plt %>%
  ggplot(data.plt) + 
  geom_ribbon(aes(x = ID.time, ymin = q1, ymax = q3), 
              fill = "steelblue", alpha = 0.5) +   ## "cadetblue"
  geom_line(aes(x = ID.time, y = q2), col = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  xlab("Time") +
  ylab("Contribution to log(Incidence rate)") +
  scale_y_continuous() +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  theme_bw() + 
  facet_wrap(. ~ STE_NAME21, nrow = 3, ncol = 3, scales = "free")
temporal_effects

ggsave(paste0(getwd(),'/figs/temporal_effect.pdf'), height = 30, width = 25, units = "cm")

## ------- plot estimated temperature-influenza  curves
dim(mdl$summary.random$`inla.group(tas, n = 20`) 
# get smooth terms for temperature
tas_eff     <- mdl$summary.random$`inla.group(tas, n = 20)`
tasLag1_eff <- mdl$summary.random$`inla.group(tasLag1, n = 20)`

tas_eff_df <- tas_eff %>%
  mutate(RR_mean = exp(mean - mean(mean)),
         RR_lower = exp(`0.025quant` - mean(mean)),
         RR_upper = exp(`0.975quant` - mean(mean))
  ) 
tasLag1_eff_df <- tasLag1_eff %>%
  mutate(RR_mean = exp(mean - mean(mean)),
         RR_lower = exp(`0.025quant` - mean(mean)),
         RR_upper = exp(`0.975quant` - mean(mean))
  ) 

RR_tas_effect <- 
  ggplot(tas_eff_df) + 
  geom_ribbon(aes(x = ID, ymin = RR_lower, ymax = RR_upper), 
              fill = "grey", alpha = 0.5) +   ## "cadetblue"
  geom_line(aes(x = ID, y = RR_mean), col = "black") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey9") +
  labs(x="Temperature (°C)", y = "Relative risk") +
  # theme(panel.grid.major = element_line(colour = gray(0.5), linetype = "dashed",
  #                                       linewidth  = 0.1),panel.background = element_blank() )+ 
  theme(panel.background = element_blank() )+
  scale_y_continuous() +
  scale_x_continuous() +
  theme_bw()  

RR_tas_effect
ggsave(paste0(getwd(),'/figs/RR_tas.pdf'), height = 20, width = 20, units = "cm")

RR_tasLag1_effect <- 
  ggplot(tasLag1_eff_df) + 
  geom_ribbon(aes(x = ID, ymin = RR_lower, ymax = RR_upper), 
              fill = "grey", alpha = 0.5) +   ## "cadetblue"
  geom_line(aes(x = ID, y = RR_mean), col = "black") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey9") +
  labs(x="Temperature (°C)", y = "Relative risk") +
  # theme(panel.grid.major = element_line(colour = gray(0.5), linetype = "dashed",
  #                                       linewidth  = 0.1),panel.background = element_blank() )+ 
  theme(panel.background = element_blank() )+
  scale_y_continuous() +
  scale_x_continuous() +
  theme_bw()  

RR_tasLag1_effect

ggsave(paste0(getwd(),'/figs/RR_tasLag1.pdf'), height = 20, width = 20, units = "cm")


#############################################################
## Compare spatial effects exposure-response relationship effect from the original M3 model
## and the sensitivity models for M3 using window 5
##########################################################
## Load the original and sensitivity optimal models separately and rename mdl accordingly
mdl_senstvty <- mdl
post_mean_region_senstvty <- mdl_senstvty$summary.random$ID.region$mean ## spatial effect
rm(mdl)
mdl_original <- mdl
post_mean_region_original <- mdl_original$summary.random$ID.region$mean ## spatial effect


## ---- plot exposure-response relationships ----
dim(mdl_original$summary.random$`inla.group(tas, n = 20`) 
dim(mdl_senstvty$summary.random$`inla.group(tas, n = 20`) 

# get smooth terms for temperature
tas_eff_original     <- mdl_original$summary.random$`inla.group(tas, n = 20)`
tasLag1_eff_original <- mdl_original$summary.random$`inla.group(tasLag1, n = 20)`
tas_eff_sensitivity     <- mdl_senstvty$summary.random$`inla.group(tas, n = 20)`
tasLag1_eff_sensitivity <- mdl_senstvty$summary.random$`inla.group(tasLag1, n = 20)`

tas_eff_df_original <- tas_eff_original %>%
  mutate(RR_mean = exp(mean - mean(mean)),
         RR_lower = exp(`0.025quant` - mean(mean)),
         RR_upper = exp(`0.975quant` - mean(mean))
  ) 
tasLag1_eff_df_original <- tasLag1_eff_original %>%
  mutate(RR_mean = exp(mean - mean(mean)),
         RR_lower = exp(`0.025quant` - mean(mean)),
         RR_upper = exp(`0.975quant` - mean(mean))
  ) 

tas_eff_df_sensitivity <- tas_eff_sensitivity %>%
  mutate(RR_mean = exp(mean - mean(mean)),
         RR_lower = exp(`0.025quant` - mean(mean)),
         RR_upper = exp(`0.975quant` - mean(mean))
  ) 
tasLag1_eff_df_sensitivity <- tasLag1_eff_sensitivity %>%
  mutate(RR_mean = exp(mean - mean(mean)),
         RR_lower = exp(`0.025quant` - mean(mean)),
         RR_upper = exp(`0.975quant` - mean(mean))
  ) 

## plot tas effect 
RR_tas_effect_combined <-
  ggplot() +
  # Original model
  geom_ribbon(data = tas_eff_df_original,
              aes(x = ID, ymin = RR_lower, ymax = RR_upper, fill = 'M3'),
              alpha = 0.5) + ## fill = "grey",
  geom_line(data = tas_eff_df_original,
            aes(x = ID, y = RR_mean, colour = "M3",),  linewidth = 1 ) + 
  # Sensitivity model
  geom_ribbon(data = tas_eff_df_sensitivity,
              aes(x = ID, ymin = RR_lower, ymax = RR_upper, fill = "M3_season"), ## change name accordingly
              alpha = 0.5) +
  geom_line(data = tas_eff_df_sensitivity,
            aes(x = ID, y = RR_mean, colour = "M3_season"),      linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "black") +
  scale_colour_manual(name = '', values = c('M3' = 'grey',
                                            'M3_season'= 'steelblue'))+
  scale_fill_manual(name = '', values = c('M3' = 'grey',
                                          'M3_season' = 'steelblue')) +
  labs( x = "Temperature (°C)",  y = "Relative risk") +
  scale_y_continuous() +
  scale_x_continuous() +
  theme_bw() +
  theme(panel.background = element_blank(), 
        legend.position = 'inside',
        legend.position.inside = c(0.08, 0.98), 
        legend.justification = c('left', 'top'),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        text = element_text(family = "sans", size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
  )

RR_tas_effect_combined

ggsave(paste0(getwd(),'/figs/RR_tas_season_combined.pdf'), height = 20, width = 20, units = "cm")

## plot tasLag1 effect 
RR_tasLag1_effect_combined <-
  ggplot() +
  # Original model
  geom_ribbon(data = tasLag1_eff_df_original,
              aes(x = ID, ymin = RR_lower, ymax = RR_upper, fill = 'M3'),
              alpha = 0.5) + ##  
  geom_line(data = tasLag1_eff_df_original,
            aes(x = ID, y = RR_mean, colour = "M3",),  linewidth = 1 ) +  
  # Sensitivity model
  geom_ribbon(data = tasLag1_eff_df_sensitivity,
              aes(x = ID, ymin = RR_lower, ymax = RR_upper, fill = "M3_season"), ## change name accordingly
              alpha = 0.5) +
  geom_line(data = tasLag1_eff_df_sensitivity,
            aes(x = ID, y = RR_mean, colour = "M3_season"),      linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "black") +
  scale_colour_manual(name = '', values = c('M3' = 'grey',
                                            'M3_season'= 'steelblue'))+
  scale_fill_manual(name = '', values = c('M3' = 'grey',
                                          'M3_season' = 'steelblue')) +
  labs( x = "Temperature (°C)",  y = "Relative risk") +
  scale_y_continuous() +
  scale_x_continuous() +
  theme_bw() +
  theme(panel.background = element_blank(), 
        legend.position = 'inside',
        legend.position.inside = c(0.98, 0.98), 
        legend.justification = c('right', 'top'),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        text = element_text(family = "sans", size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
  )

RR_tasLag1_effect_combined

ggsave(paste0(getwd(),'/figs/RR_tasLag1_season_combined.pdf'), height = 20, width = 20, units = "cm")

