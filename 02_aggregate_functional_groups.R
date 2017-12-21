# This scripts reads in the "raw" data that was downloaded in 01_accessing_data.
# The primary goal of this script is to generate bromeliad-level summaries that 
# can be combined with data on final insect communities.  We create summaries of
# ibuttons, functional groups and taxonomic groups.

## Some data will be held out for now. These data will appear in later versions
## of this same database. Lines were data is held out are marked with **HELD OUT**


# reading in data ---------------------------------------------------------

library(bwgtools)
library(assertr)
library(tidyverse)
library(stringr)

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


# Gustavo’s “Meso and Top” categorization ---------------------------------

meso_top <- readxl::read_excel("Meso & Top.xlsx", range = "A1:C64")
meso_top %>% tail

final_inverts


# need to combine meso_top with the final_inverts to generate Gustavo's biomass numbers.

# match names
meso_top_fixed <- meso_top %>% 
  rename(site = Site,
         species = BWGCode) %>% 
  # fix site mispellings
  mutate(site = tolower(site),
         site = if_else(site == "french g", "frenchguiana", site),
         site = if_else(site == "puerto r", "puertorico", site),
         site = if_else(site == "costa rica", "costarica", site))

unmatched_sites <- meso_top_fixed$site %>% unique %>% setdiff(final_inverts$site)

if(length(unmatched_sites) != 0) stop("there are unmatch sites!")

# should also be matches for species

unmatched_spp <- meso_top_fixed$species %>% unique %>% setdiff(final_inverts$species %>% unique)
if(length(unmatched_spp) != 0) stop("there are unmatched spp!")

meso_top_by_bromeliads <- final_inverts %>% 
  left_join(meso_top_fixed, by = c("site", "species")) %>% 
  drop_na(`Predator category`) %>% 
  select(site_brom.id, site, pred_cat = `Predator category`, biomass) %>% 
  group_by(site_brom.id, site, pred_cat) %>% 
  summarize(total_biomass = sum(biomass)) %>% 
  ungroup %>% 
  spread(pred_cat, total_biomass, fill = 0) %>% 
  rename(meso_predator_biomass = Meso, 
         top_predator_biomass = Top)


# join everything together ------------------------------------------------

bromeliad_variables <- read_csv("Data/BWG_bromeliad_variables.csv")

fulldata  <-  bromeliad_variables %>%
  left_join(trophic_bio, by = "site_brom.id") %>%
  left_join(ord_bio, by = "site_brom.id")%>%
  left_join(subclass_bio, by = "site_brom.id")%>%
  left_join(family_bio, by = "site_brom.id")%>%
  left_join(subfamily_bio, by = "site_brom.id")%>%
  left_join(genus_bio, by = "site_brom.id")%>%
  left_join(ibutton_data_cardoso_corrected, by = "site_brom.id") %>% #note ibutton data is 205 rows not 210
  # join on the "mesopredator" and "toppredator" biomasses
  left_join(meso_top_by_bromeliads) %>% 
  replace_na(list(meso_predator_biomass = 0, 
                  top_predator_biomass = 0))

#predict a missing colombia maxvol measurement - also changed on rawdata
maxvolmod<-lm(log(maxvol)~leaf.number, data=subset(fulldata,site=="colombia"))
par(mfrow=c(2,2)); plot(maxvolmod) #+mean.diam+catchment.area excluded as not useful
summary(maxvolmod)#r2 is 30
fulldata$maxvol[fulldata$site_brom.id=="colombia_29"]<-exp(predict(maxvolmod, data.frame(leaf.number= 33))) #predicted missing colombia_29 maxvol is 564

glimpse(fulldata)


# writing out data --------------------------------------------------------


# remove taxonomic groups under subfamily, remove variables ending _taxa or
# _abundance, and remove meso_predator_biomass and top_predator_biomass. We
# could save this csv as a different name, and I can just adjust my code in
# 01_datasets_for_paper1 to read in new name.

fulldata %>% 
  select(-(Aelosoma_bio:Wyeomyia_bio), 
         # no 'abundance' in paper 1 
         -contains("abundance"),
         - meso_predator_biomass, 
         - top_predator_biomass) %>% 
  write_csv("Data/BWG_wide_functional_groups_ibuttons.csv")

# take the long-format functional group data (abudance and biomass) and add ibuttons
# Start with the data that contains everything: all traits and all functional groups
invert_traits %>% 
  # summarize by grouping and summing over functional and taxonomic groups. throw in order for good measure
  sum_func_groups( grps = list(~site, 
                               ~site_brom.id,
                               ~predation,
                               ~functional_group, 
                               ~ord,
                               ~family)) %>% 
  # finally add in the temperature data
  left_join(ibutton_data_cardoso_corrected %>% select(site_brom.id, mean_temp, sd_mean_temp)) %>% 
  # **HELD OUT**
  # filter out data that will appear in later paper
  # "total abundance" "total taxa" (perhaps
  # "predation"),
  select(-total_abundance, -total_taxa, -predation) %>% 
  write_csv("Data/BWG_long_functional_groups_temp.csv")

