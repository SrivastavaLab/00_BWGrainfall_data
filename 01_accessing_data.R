# setwd("C:/Users/Diane/Dropbox/BWG Drought experiment/Paper 1_thresholds")
source("../Rscripts/checkbwgversion.R")
rdrop2::drop_auth(cache = FALSE)

#to manually update bwgtools, run this, may need to install dependent packages too
# devtools::install_github("SrivastavaLab/bwgtools", dependencies = TRUE)
# devtools::install_github("Srivastavalab/bwgtools")

# Loading packages ------------------------------------

library(bwgtools)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(curl)
library(purrr)
library(lubridate)
library(assertr)

# obtain data -----------------------------------------

# site.info
site_info <- combine_tab(sheetname = "site.info")

# site.weather
site_weather <- combine_tab(sheetname = "site.weather")

# bromeliad.physical
bromeliad_physical <- combine_tab(sheetname = "bromeliad.physical")

# leaf.waterdepths
leaf_waterdepths <- combine_tab(sheetname = "leaf.waterdepths")

# terrestrial.taxa
#terrestrial_taxa <- combine_tab(sheetname = "terrestrial.taxa")

# bromeliad.terrestrial
bromeliad_terrestrial <- combine_tab(sheetname = "bromeliad.terrestrial")

# bromeliad.ibuttons
ibuttons <- combine_tab(sheetname = "bromeliad.ibuttons")

# bromeliad.initial.inverts
initial_inverts <- combine_tab(sheetname = "bromeliad.initial.inverts")

# bromeliad.final.inverts
final_inverts <- bwgtools::combine_tab(sheetname = "bromeliad.final.inverts")


### trait data -- from database
bwg_names <- bwgdata::bwg_get("species")

# checking final inverts, merging with bwg_names --------------------------------------------------

# checking for duplicate columns -- there should be 30 sites with 7 unique names
complete_sites <- final_inverts %>%
  select(site, site_brom.id) %>%
  distinct %>%
  group_by(site) %>%
  tally %>%
  assertr::verify(.$n == rep(30, 7))

stopifnot(identical(structure(list(site = c("argentina", "cardoso", "colombia", "costarica", 
                                            "frenchguiana", "macae", "puertorico"),
                                   n = c(30L, 30L, 30L, 
                                         30L, 30L, 30L, 30L)),
                              class = c("tbl_df", "tbl", "data.frame"
                              ),
                              row.names = c(NA, -7L),
                              .Names = c("site", "n")),
                    complete_sites))

final_inverts %>% 
  write_csv("Data/BWG_final_invertebrates.csv")

# bwg_names processing and preparation ------------------------------------



# that is a list and it causes
## pain in subsequent joins. fortunately we don't need it.
## Dropping it for now; future versions won't create this
## problem in the 1st place
if (inherits(bwg_names$names, "list")) {
  bwg_names <- bwg_names %>% 
    select(-names)
}

# fix the variable classes
bwg_names <- bwg_names %>% 
  mutate_all(funs(parse_guess)) 


## who are these pesky animals who are neither predator, nor prey, nor NA?
## (there should be none at all)
bwg_names %>% 
  assert(in_set(c("predator","prey")), predation) %>% 
  filter(!(predation %in% c("predator","prey", NA_character_))) %>% 
  {stopifnot(nrow(.) == 0)}


### merge with taxonomy and "canonical" traits
invert_traits <- final_inverts %>% 
  rename(bwg_name = species) %>% 
  left_join(bwg_names, by = "bwg_name")


# Hydrology -- removing leaky leaves ---------------------------------------------------

## using this very long data.frame, identify leaves which fail the above criteria. drop them
hydro_rogue_hunting <- hydro_variables(leaf_waterdepths, site_info, bromeliad_physical, dohydro = FALSE)


# 3. correlation (pearsons r) in water depths between lateral leaves <0.30

# calculated by simply spreading the leaf label into columns (so that
# measurements on the same date are on the same row) and calcualting the
# pairwise correlation coefficient
uncorrelated_leaves <- hydro_rogue_hunting %>%
  spread(leaf, depth) %>% 
  select(-date) %>% 
  nest(leafa, leafb) %>% 
  mutate(leaf_corr = map_dbl(data, ~ cor(.x$leafa, .x$leafb, use = "pairwise.complete.obs"))) %>% 
  select(-data) %>% 
  filter(leaf_corr < 0.3)

# Because CV and the mean are very related measurements, it makes sense to
# calculate them at the same time, for every leaf, then reshape the data to get a Rogue list
hydro_leaf_summaries <- hydro_rogue_hunting %>% 
  select(-date) %>% 
  ungroup %>% 
  nest(depth) %>% 
  mutate(leaf_summaries = map(data,
                              ~ data_frame(mean_depth = mean(.x$depth, na.rm = TRUE),
                                           sd_depth   = sd(.x$depth, na.rm = TRUE), 
                                           cv_depth   = sd_depth / mean_depth))) %>% 
  select(-data) %>% 
  unnest(leaf_summaries)

# 2. Same rogue leaf in (1) has 5x larger cv than other lateral leaf

# again, spreading leaves to make a pairwise comparison. Doing perhaps a weird 
# formula, where the difference between CVs is dividied by the smallest CV (i
# just made that up, does it look right?)
different_CV_leaves <- hydro_leaf_summaries %>% 
  select(-mean_depth, -sd_depth) %>% 
  spread(leaf, cv_depth) %>% 
  ## where is the difference in CV x5
  rowwise %>% 
  filter((abs(leafa - leafb) / min(leafa, leafb)) >= 3.5)


# 1. Mean water depth of rogue (i.e. leaky) leaf <7mm, other lateral leaf
# >7mm (a better rule might be #days with <5mm water and here 70% or more
#       days with <5mm water could be a good guide to rogue leaves)

# the simplest way to calculate this is just by using any to pull out bromeliads
# with suspiciously disparate mean depths
very_different_means <- hydro_leaf_summaries %>% 
  select(-sd_depth, -cv_depth) %>% 
  # where is at least one leaf <7mm?
  group_by(site_brom.id) %>% 
  filter(any(mean_depth < 7) & any(mean_depth > 7))
  
# out of curiosity let us calculate the second criteria
very_low_days <- hydro_rogue_hunting %>% 
  drop_na(depth) %>% 
  nest(date, depth) %>% 
  mutate(ndays_low = map_dbl(data, ~ sum(.x$depth <= 5) / nrow(.x))) %>% 
  filter(ndays_low >= 0.7)

# Are these different at all? 

very_low_days %>% anti_join(very_different_means) # present in very_low not in very_different

very_different_means %>% anti_join(very_low_days) # presente in very_different not very_low


# 4. Leaf starts expt dry (<5mm at least once in first 2 days)

# This is done in a two-step process, simply because the process is slow. First,
# filter days to within the first few days of that block (thanks to the grouping
# variable that is already applied, this minimum is calculated within each 
# bromeliad, which are nested within temporal block.). Then we drop NA because
# not everybody measured bromeliads, and some will start with NA.
non_na_first_days <- hydro_rogue_hunting %>% 
  filter(date <= min(date) + days(2)) %>% 
  drop_na(depth)

dry_right_away <- non_na_first_days %>% 
  ungroup %>% 
  filter(depth < 2) %>% 
  select(-date, -depth) %>% 
  distinct


# semi_join all these together, to create a Rogue's gallery
rogue_leaves <- very_different_means %>% 
  semi_join(dry_right_away) %>%
  semi_join(different_CV_leaves) %>% 
  semi_join(uncorrelated_leaves) %>% 
  ungroup


# anti_join hydro with these
leaf_waterdepths_no_rogue <- hydro_rogue_hunting %>% 
  anti_join(rogue_leaves, by = c("site_brom.id", "site", "trt.name", "leaf", "watered_first", "temporal.block"))
# run water_summary_calc on this list (possibly after rearranging? and nesting again?)

# water
leaf_water_nested <- leaf_waterdepths_no_rogue %>% 
  dplyr::arrange(site, date, trt.name, leaf) %>%
  nest(date, depth)

hydro_calculated <- leaf_water_nested %>% 
  mutate(water_variables = map(data, ~water_summary_calc(.x$depth)))


# Hydrology -- the final version of hydro variables ------------------------------------

hydro <- hydro_calculated %>% 
  select(-data) %>% 
  unnest(water_variables)

## summarize hydro:
mean_hydro <- hydro %>%
  ungroup %>% 
  select_(.dots = c("site_brom.id", "mean.depth", "max.depth", "min.depth", "amplitude", 
                    "net_fluc", "tot_fluc", "sd.depth", "cv.depth", "wetness", "prop.overflow.days", 
                    "prop.driedout.days", "long_dry", "long_wet", "n_driedout", "n_overflow")) %>% 
  group_by(site_brom.id) %>%
  summarise_each(funs(median(., na.rm = TRUE))) %>%
  replace_na(list(long_dry = 0, long_wet = 0, n_driedout = 0, n_overflow = 0))

when_last <- hydro %>%
  ungroup %>%
  select(site_brom.id, last_dry, last_wet) %>%
  group_by(site_brom.id) %>%
  summarise_each(funs(min(., na.rm = TRUE))) %>%
  replace_na(list(last_dry = 65, last_wet = 65))

brom_hydro <- mean_hydro %>% left_join(when_last)


# decomposition -----------------------------------------------------------

decompositon <- get_decomp(bromeliad_physical = bromeliad_physical)


# as yet unused variables -------------------------------------------------

# rainfall
# ibutton metrics
# initial_inverts


# Bromeliads --------------------------------------------------------------

bromeliad_variables <- bromeliad_physical %>%
  left_join(decompositon, by = c("site_brom.id", "trt.name")) %>% 
  verify(nrow(.) == 210)


# write data ------------------------------------------

write_csv(invert_traits, "Data/invert_data.csv")

