
# move this part back to paper 1
# filter out only aquatics
# filter(realm == "aquatic")


#summarize functional groups
func_groups <- sum_func_groups(invert_traits,
                               grps = list(~site, ~site_brom.id,
                                           ~predation, ~functional_group))



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


# Abundance and biomass at functional and trophic levels ------------------------------------------------------------


# summarize trophic groups
trophic_groups <- sum_func_groups(invert_traits,
                                  grps = list(~site, ~site_brom.id,
                                              ~predation))


# note that func_groups is calculated from invert_traits

#summarize the biomass of each functional group. Produce one row per bromeliad, 
#which 210 bromeliads. There should be NO missing biomass. If an animal has no
#functional group, it has "NA_bio"
func_bio <- func_groups %>%
  ungroup %>%
  ## total_biomass should be present for every bromeliad
  assert(not_na, total_biomass) %>%
  select(site_brom.id, functional_group, total_biomass) %>%
  # assert(not_na, functional_group) %>%
  group_by(site_brom.id, functional_group) %>%
  summarise_each(funs(sum)) %>%
  mutate(functional_group = if_else(is.na(functional_group), paste0("predation_", functional_group), functional_group),
         functional_group = paste0(functional_group, "_bio")) %>%
  spread(functional_group, total_biomass, fill = 0) %>% 
  verify(nrow(.) == 210)

func_abd <- func_groups %>%
  ungroup %>%
  ## total_abundance should be present for every bromeliad
  assert(not_na, total_abundance) %>%
  ## if it is a piercer, append predation to it
  mutate(functional_group = if_else(functional_group == "piercer", 
                                    paste0(predation, "_", functional_group),
                                    functional_group)) %>% 
  select(site_brom.id, functional_group, total_abundance) %>%
  mutate(functional_group = if_else(is.na(functional_group), paste0("predation_", functional_group), functional_group),
         functional_group = paste0(functional_group, "_abd")) %>%
  spread(functional_group, total_abundance, fill = 0) %>% 
  verify(nrow(.) == 210)


## total biomass for each trophic level. if an animal has no known trophic
## level, it is "NA_bio"
trophic_bio <- trophic_groups %>%
  ungroup %>%
  ## total_biomass should be present for every bromeliad
  assert(not_na, total_biomass) %>%
  select(site_brom.id, predation, total_biomass) %>%
  mutate(predation = if_else(is.na(predation), paste0("predation_", predation), predation),
         predation = paste0(predation, "_bio")) %>%
  spread(predation, total_biomass, fill = 0) %>% 
  verify(nrow(.) == 210)

#abundance of each trophic level
trophic_abd <- trophic_groups %>%
  ungroup %>%
  ## total_biomass should be present for every bromeliad
  assert(not_na, total_biomass) %>%
  select(site_brom.id, predation, total_abundance) %>%
  mutate(predation = if_else(is.na(predation), paste0("predation_", predation), predation),
         predation = paste0(predation, "_abd")) %>%
  spread(predation, total_abundance, fill = 0) %>% 
  verify(nrow(.) == 210)


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



fulldata  <-  bromeliad_variables %>%
  left_join(func_bio, by = "site_brom.id") %>%
  left_join(func_abd, by = "site_brom.id") %>%
  left_join(trophic_bio, by = "site_brom.id") %>%
  left_join(trophic_abd, by = "site_brom.id")%>%
  left_join(ord_bio, by = "site_brom.id")%>%
  left_join(subclass_bio, by = "site_brom.id")%>%
  left_join(family_bio, by = "site_brom.id")%>%
  left_join(subfamily_bio, by = "site_brom.id")%>%
  left_join(genus_bio, by = "site_brom.id")%>%
  left_join(brom_hydro, by = "site_brom.id")%>%
  left_join(ibutton_data, by = "site_brom.id") #note ibutton data is 205 rows not 210

#predict a missing colombia maxvol measurement - also changed on rawdata
maxvolmod<-lm(log(maxvol)~leaf.number, data=subset(fulldata,site=="colombia"))
par(mfrow=c(2,2)); plot(maxvolmod) #+mean.diam+catchment.area excluded as not useful
summary(maxvolmod)#r2 is 30
fulldata$maxvol[fulldata$site_brom.id=="colombia_29"]<-exp(predict(maxvolmod, data.frame(leaf.number= 33))) #predicted missing colombia_29 maxvol is 564

# setwd("C:/Users/Diane/Dropbox/BWG Drought experiment/Paper 1_thresholds/Data")

write_csv(fulldata, "Data/paper_1_data.csv")