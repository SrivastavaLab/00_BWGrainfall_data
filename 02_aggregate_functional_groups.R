# This scripts reads in the "raw" data that was downloaded in 01_accessing_data.
# The primary goal of this script is to generate bromeliad-level summaries that 
# can be combined with data on final insect communities.  We create summaries of
# ibuttons, functional groups and taxonomic groups.


# reading in data ---------------------------------------------------------

library(bwgtools)
library(assertr)
library(tidyverse)

#before you can merge the functional groups and summarize with them you have to 
#get this info from the BWG, then also the traits. Read it in off a file, select
#what is necessary and THEN do this work
# "what is necessary" == everything not Tachet

final_inverts <- read_csv("Data/BWG_final_invertebrates.csv")
bwg_names <- read_csv("Data/BWG_final_invertebrate_traits.csv")

## merge with taxonomy and
#"canonical" traits
invert_traits <- final_inverts %>% 
  rename(bwg_name = species) %>% 
  left_join(bwg_names, by = "bwg_name") %>% 
  filter(realm == "aquatic")


#summarize functional groups
func_groups <- sum_func_groups(invert_traits,
                               grps = list(~site, 
                                           ~site_brom.id,
                                           ~predation,
                                           ~functional_group))



## At this point, we should have combined all species within trophic (predation) and functional groups. 
## This implies that there should be no duplication of either category within func_groups
## if that were NOT true, then there are duplicate rows. 
## select columns that must be unique, drop duplicates, check for decrease in nrow
distinct_rows <- func_groups %>% 
  ungroup %>% 
  select(site, site_brom.id, predation, functional_group) %>% 
  distinct %>% 
  nrow

func_groups <- func_groups %>% 
  assertr::verify(nrow(.) == distinct_rows)



## check a known problem:
func_groups %>% filter(site_brom.id == "cardoso_DEC25")

#' This table shows that there is a piercer in both predator and prey categories

# Abundance and biomass at functional and trophic levels ------------------------------------------------------------

# note that func_groups is calculated from invert_traits

#summarize the biomass of each functional group. Produce one row per bromeliad, 
#which 210 bromeliads. There should be NO missing biomass. If an animal has no
#functional group, it has "NA_bio"
func_bio_abd <- func_groups %>%
  ungroup %>%
  ## total_biomass should be present for every bromeliad
  assert(not_na, total_biomass, total_abundance, total_taxa) %>%
  select(site_brom.id, predation, functional_group, total_biomass, total_abundance, total_taxa) %>%
  gather(response, value, starts_with("total")) %>% 
  unite("pred_func_resp", predation, functional_group, response) %>% 
  spread(pred_func_resp, value, fill = 0) %>% 
  verify(nrow(.) == 210)

## total biomass for each trophic level. if an animal has no known trophic
## level, it is "NA_bio"
trophic_bio <-  func_bio_abd %>% glimpse %>% 
  mutate(predator_biomass = 
           predator_engulfer_total_biomass + 
           predator_piercer_total_biomass,
         prey_biomass     = 
           prey_filter.feeder_total_biomass + 
           prey_gatherer_total_biomass + 
           prey_piercer_total_biomass + 
           prey_scraper_total_biomass + 
           prey_shredder_total_biomass) %>% 
  mutate(predator_abundance = 
           predator_engulfer_total_abundance + 
           predator_piercer_total_abundance,
         prey_abundance     = 
           prey_filter.feeder_total_abundance + 
           prey_gatherer_total_abundance + 
           prey_piercer_total_abundance + 
           prey_scraper_total_abundance + 
           prey_shredder_total_abundance)

# taxonomic level biomass -------------------------------------------------

##begin taxonomic level biomass

sum_taxa_groups<-function (merged_data, grps = list(~site, ~site_brom.id, ~ord,~subclass,
                                                    ~family, ~subfamily, ~genus))
{
  merged_data %>% dplyr::group_by_(.dots = grps) %>% dplyr::summarize(total_abundance = sum(abundance),
                                                                      total_biomass = sum(biomass), total_taxa = n())
}

ord_groups <- sum_taxa_groups(invert_traits,
                              grps = list(~site, ~site_brom.id, ~ord, ~subclass,
                                          ~family, ~subfamily, ~genus))

check_rename_summarize_spread <- function(df, colname){
  df_ready <- df %>%
    ungroup %>%
    ## total_biomass should be present for every bromeliad
    assert(not_na, total_biomass) %>%
    select_(.dots = c("site_brom.id", colname, "total_biomass"))
  
  df_ready[[colname]] <- if_else(is.na(df_ready[[colname]]),
                                 paste0(colname, "_", df_ready[[colname]]),
                                 df_ready[[colname]])
  
  df_ready[[colname]] <-  paste0(df_ready[[colname]], "_bio")
  
  df_ready %>%
    group_by_(.dots = c("site_brom.id", colname)) %>%
    summarize(total_biomass = sum(total_biomass)) %>%
    spread_(key_col = colname, value_col = "total_biomass", fill = 0)
}

ord_bio <- check_rename_summarize_spread(ord_groups, "ord")

subclass_bio <- check_rename_summarize_spread(ord_groups, "subclass")

family_bio <- check_rename_summarize_spread(ord_groups, "family")

subfamily_bio <- check_rename_summarize_spread(ord_groups, "subfamily")

genus_bio <- check_rename_summarize_spread(ord_groups, "genus")

## merge all the bromeliad-level variables
## start with bromeliad.physical because it must be complete!



# summarize ibuttons ------------------------------------------------------

# read those ibuttons in
ibuttons <- read_csv("Data/BWGrainfall_long_ibuttons.csv")

# summarize ibuttons for combining with bromeliad-level data later
ibutton_data <- ibuttons %>%
  group_by(site, site_brom.id) %>%
  summarise(mean_max = mean(max.temp, na.rm = TRUE), mean_min = mean(min.temp, na.rm = TRUE),
            mean_mean = mean(mean.temp, na.rm = TRUE), sd_max = sd(max.temp, na.rm = TRUE),
            sd_min = sd(min.temp, na.rm = TRUE), sd_mean = sd(mean.temp, na.rm = TRUE),
            cv_max = 100*(sd_max/mean_max), cv_min = 100*(sd_min/mean_min),
            cv_mean = 100*(sd_mean/mean_mean)) %>%
  ungroup %>% 
  gather(variable, observed, 3:11) %>%
  replace_na(list(observed = NA_real_)) %>%
  select(-site) %>%
  spread(variable, observed, fill = NA_real_) %>%
  rename(max_temp = mean_max, min_temp = mean_min, mean_temp = mean_mean,
         sd_max_temp = sd_max, sd_min_temp = sd_min, sd_mean_temp = sd_mean,
         cv_max_temp = cv_max, cv_min_temp = cv_min, cv_mean_temp = cv_mean)

ibutton_data_cardoso_corrected <- ibutton_data %>% 
  mutate(site_brom.id = str_replace(site_brom.id, " [A-Z]$", ""),
         site_brom.id = str_replace(site_brom.id, "[A-Z]$", ""),
         site_brom.id = str_replace_all(site_brom.id, " ", ""))

# join everything together ------------------------------------------------

bromeliad_variables <- read_csv("Data/BWG_bromeliad_variables.csv")

fulldata  <-  bromeliad_variables %>%
  left_join(trophic_bio, by = "site_brom.id") %>%
  left_join(ord_bio, by = "site_brom.id")%>%
  left_join(subclass_bio, by = "site_brom.id")%>%
  left_join(family_bio, by = "site_brom.id")%>%
  left_join(subfamily_bio, by = "site_brom.id")%>%
  left_join(genus_bio, by = "site_brom.id")%>%
  left_join(ibutton_data_cardoso_corrected, by = "site_brom.id") #note ibutton data is 205 rows not 210

#predict a missing colombia maxvol measurement - also changed on rawdata
maxvolmod<-lm(log(maxvol)~leaf.number, data=subset(fulldata,site=="colombia"))
par(mfrow=c(2,2)); plot(maxvolmod) #+mean.diam+catchment.area excluded as not useful
summary(maxvolmod)#r2 is 30
fulldata$maxvol[fulldata$site_brom.id=="colombia_29"]<-exp(predict(maxvolmod, data.frame(leaf.number= 33))) #predicted missing colombia_29 maxvol is 564

glimpse(fulldata)

write_csv(fulldata, "Data/BWG_wide_functional_groups_ibuttons.csv")