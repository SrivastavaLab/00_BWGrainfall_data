
library(ggplot2)
library(viridis)

# bromeliad.ibuttons
ibuttons <- combine_tab(sheetname = "bromeliad.ibuttons")


ibuttons %>% 
  filter(site != "puertorico") %>%
  ggplot(aes(x = date, y = mean.temp, ymin = min.temp, ymax = max.temp, group = site_brom.id)) +
  geom_pointrange(alpha = 0.2) + 
  facet_wrap(~site, scales = "free_x")


ibuttons %>% 
  filter(site != "puertorico") %>%
  ggplot(aes(x = date, y = mean.temp, colour = trt.name)) +
  geom_line(alpha = 0.7) + 
  facet_wrap(~site, scales = "free_x") + 
  scale_color_viridis(option = "A", discrete = TRUE)

# whaaa what happened to CR's treatments???

ibuttons %>% 
  filter(site != "puertorico") %>% 
  filter(site == "costarica")

# ugh they would have to be read in and reworked. 
