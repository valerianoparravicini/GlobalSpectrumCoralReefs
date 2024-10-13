library(tidyverse)
targets::tar_load(summary_fish_functions)

data <- summary_fish_functions 
data[data$location == "Kon\xe9" , "location"] <- "Kone"

data <- data %>%
  group_by(location, dataset_id, site) %>%
  drop_na(location) %>%
  mutate(time = (year - min(year))/(max(year) - min(year)))  %>%
  mutate(loc_site = paste(location, site, sep = "_"))


library(tidybayes)

####### cc#######
sub <- data %>%
  filter(!is.na(coral_cover)) %>%
  filter(nyear>4) %>%
  filter(!location == "KONIAMBO")

ggplot(sub) +
  geom_point(aes(x = year, y = coral_cover)) +
  facet_wrap(~loc_site, scales = "free")

ggplot(sub) +
  geom_point(aes(x = sqrt(coral_cover), y = log(outNP))) +
  geom_smooth(aes(x = sqrt(coral_cover), y = log(outNP)), method = "lm") +
  theme_bw()  
ggplot(sub) +
  geom_point(aes(x = sqrt(coral_cover), y = log(exNP))) +
  geom_smooth(aes(x = sqrt(coral_cover), y = log(exNP)), method = "lm") +
  theme_bw()  
ggplot(sub) +
  geom_point(aes(x = sqrt(coral_cover), y = log(exNP))) +
  geom_smooth(aes(x = sqrt(coral_cover), y = log(exNP), group = site, color = zone), method = "lm") +
  theme_bw() +
  facet_wrap(~location, scale = "free")

######## plots ############
data <- summary_fish_functions
subset <- data %>%
  filter(location %in% c("Moorea", "Tetiaroa", "MaheNWCarbonate", 
  "MaheWCarbonate", "CousinCarbonate", "CousinPatch", "Hoga", "Tubuai") | 
    site %in% c("Moorea", "Tetiaroa", "MaheNWCarbonate", 
                "MaheWCarbonate", "CousinCarbonate", "CousinPatch", "Hoga", "Tubuai"))
subset <- data %>%
  filter(location %in% c("Hoga"))

ggplot(subset) +
  geom_point(aes(x = year, y = log(outNP), color = zone)) +
  geom_smooth(aes(x = year, y = log(outNP), color = zone)) +
  facet_wrap(~site) 

corcov <- coral %>%
  filter(site %in% c("Tiahura34", "MaheNWCarbonate", "MaheWCarbonate", "CousinCarbonate", "CousinPatch")|
           location %in% c( "Tetiaroa", "Hoga", "Tubuai", "Montebellos")) 


ggplot(subset) +
  geom_point(aes(x = year, y = coral_cover, color = zone)) +
  geom_smooth(aes(x = year, y = coral_cover)) +
  facet_wrap(~site, scale = "free") 

########### cousin ################

logz <- function(x){
  x <- log(x)
  (x - mean(x))/sd(x)
}

cousin_c <- corcov %>%
  filter(site %in% c("CousinPatch")) %>%
  mutate(fase = case_when(year == 1994 ~ "before",
                          year %in% c(2005) ~ "during",
                          year %in% c(2011) ~ "post")) %>%
  drop_na(fase)

cousin_f <- subset %>%
  filter(site %in% c("CousinPatch"))%>%
  mutate(fase = case_when(year == 1994 ~ "before",
                          year %in% c(2005) ~ "during",
                          year %in% c(2011) ~ "post")) %>%
  drop_na(fase)

ggplot(cousin_f) +
  geom_boxplot(aes(x = fase, y = plank, color = location))

cousin_long <- cousin_f %>%
  ungroup() %>%
  select(fase, exP, exNP, herb, plank, prod) %>%
  pivot_longer(-c(fase)) %>%
  filter(value>0) %>%
  group_by(name) %>%
  mutate(ref = median(value[fase == "before"])) %>%
  mutate(value_scaled = value/ref)

ggplot(cousin_long) +
  geom_boxplot(aes(x = name, y = value_scaled, color = fase), outlier.alpha = 0) +
  coord_flip()

fit_cousin <- brm(value_scaled ~ fase*name, data = cousin_long, backend = "cmdstanr", family = "student")
summary(fit_cousin)


nd <- select(cousin_long, name, fase) %>% unique()
pred <- fitted(fit_cousin, newdata = nd, summary = F)[1:1000,]
pred <- pred %>% as.data.frame() %>% pivot_longer(1:15) %>%
  mutate(fun = rep(nd$name, 1000),
         fase = rep(nd$fase, 1000))

pred <- cbind(nd, fitted(fit_cousin, newdata = nd, summary = T)) %>%
  group_by(name) %>%
  mutate(ref = Estimate[fase == "before"])


library(fishualize)
library(tidybayes)
library(ggplot2)
library(tidyverse)

ggplot(cousin_c) +
  geom_point(aes(x = fase, y = coral_cover))

cousin_c2 <- cousin_c %>%
  group_by(year) %>%
  summarise(coral_cover = median(coral_cover))

pcous_cor <- ggplot(cousin_c2) +
  geom_smooth(aes(x = year, y = (coral_cover)), color = "black") +
  geom_point(aes(x = year, y = (coral_cover)), data = cousin_c, 
             position = position_jitter(width = 1)) +
  theme_classic() +
  scale_y_continuous(trans = "sqrt") +
  labs(y = "Coral cover", x = "Year")
pcous_cor

plot_cousin <- 
ggplot(pred) +
  geom_hline(yintercept = 1, linetype = 3)  +
  
  geom_point(aes(x = fase, y = Estimate/ref, color = name),
             position = position_dodge(width = 0.5), size = 3) +
  geom_line(aes(x = as.numeric(as.factor(fase)), y = Estimate/ref, color = name,
               linetype = name),
            position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = `Q2.5`/ref, ymax = `Q97.5`/ref, x = fase, color = name),
                alpha = 0.6,  position = position_dodge(width = 0.5), width = 0, size = 2) +
  theme_bw() +
  scale_color_fish_d() +
  scale_fill_fish_d() + 
  labs(x = "Fase", y = "Scaled function", color = "Function", linetype = "Function",
       title = "Cousin (shifted)") +
  theme(text = element_text(size = 14)) 
plot_cousin


cousin_c2 <- cousin_c %>%
  group_by(fase) %>%
  summarise(coral_cover = median(coral_cover))

cousin_f2 <- cousin_f %>%
  group_by(fase)  %>%
  summarise_if(is.numeric, median)



########### mahe ################
mahe_c <- corcov %>%
  filter(site %in% c("MaheNWCarbonate")) %>%
  mutate(fase = case_when(year == 1994 ~ "before",
                          year %in% c(2005)~ "during",
                          year == 2011 ~ "post")) %>%
  drop_na(fase)
mahe_f <- subset %>%
  filter(site %in% c("MaheNWCarbonate"))%>%
  mutate(fase = case_when(year == 1994 ~ "before",
                          year %in% c(2005) ~ "during",
                          year == 2011 ~ "post")) %>%
  drop_na(fase)

ggplot(mahe_f) +
  geom_boxplot(aes(x = fase, y = plank, color = location))

mahe_c2 <- mahe_c %>%
  group_by(year) %>%
  summarise(coral_cover = median(coral_cover))
mahe_f2 <- mahe_f %>%
  group_by(fase) %>%
  summarise_if(is.numeric, median)



mahe_long <- mahe_f %>%
  ungroup() %>%
  select(fase, exP, exNP, herb, plank, prod) %>%
  pivot_longer(-c(fase)) %>%
  filter(value>0) %>%
  #mutate_if(is.numeric, log) %>%
  group_by(name) %>%
  mutate(ref = median(value[fase == "before"])) %>%
  mutate(value_scaled = value/ref)

ggplot(mahe_long) +
  geom_boxplot(aes(x = name, y = value_scaled, color = fase), outlier.alpha = 0)

ggplot(mahe_long) +
  geom_boxplot(aes(x = fase, y = value_scaled, color = name), outlier.alpha = 0) 

library(brms)
fit_mahe <- brm(value_scaled ~ fase*name, data = mahe_long, 
                backend = "cmdstanr", family = "student")
summary(fit_mahe)
#marginal_effects(fit_cousin)


nd <- select(mahe_long, name, fase) %>% unique()
pred <- fitted(fit_mahe, newdata = nd, summary = F)[1:1000,]
pred <- pred %>% as.data.frame() %>% pivot_longer(1:15) %>%
  mutate(fun = rep(nd$name, 1000),
         fase = rep(nd$fase, 1000))

pred <- cbind(nd, fitted(fit_mahe, newdata = nd, summary = T)) %>%
  group_by(name) %>%
  mutate(ref = Estimate[fase == "before"]) 

pmahe_cor <- ggplot(mahe_c2) +
  geom_smooth(aes(x = year, y = (coral_cover)), color = "black") +
  geom_point(aes(x = year, y = (coral_cover)), data = mahe_c, 
             position = position_jitter(width = 1)) +
  theme_classic() +
  scale_y_continuous(trans = "sqrt") +
  labs(y = "Coral cover", x = "Year")


plot_mahe <- 
ggplot(pred) +
  geom_hline(yintercept = 1, linetype = 3)  +
  
  geom_point(aes(x = fase, y = Estimate/ref, color = name),
             position = position_dodge(width = 0.5), size = 3) +
  geom_line(aes(x = as.numeric(as.factor(fase)), y = Estimate/ref, color = name,
                linetype = name),
            position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = `Q2.5`/ref, ymax = `Q97.5`/ref, x = fase, color = name),
                alpha = 0.6,  position = position_dodge(width = 0.5), width = 0, size = 2) +
  
  theme_bw() +
  scale_color_fish_d() +
  scale_fill_fish_d() + 
  labs(x = "Fase", y = "Scaled function", color = "Function", linetype = "Function",
       title = "Mahe (Recovered)") +
  theme(text = element_text(size = 14)) +
  scale_y_continuous(limits = c(0,3))
plot_mahe

library(patchwork)
plot_sey <- pmahe_cor +  pcous_cor + plot_mahe + theme(legend.position = "none") + 
  plot_cousin  + plot_layout(ncol = 2, heights = c(1,4))
plot_sey
ggsave("output/ts_seychelles.png", width = 10, height = 5)


#############" Tetiaroa #########""
tet_c <- corcov %>%
  filter(location %in% c("Tetiaroa")) %>%
  mutate(fase = case_when(year == 2005 ~ "before",
                          year %in% c(2013)~ "during",
                          year == 2019 ~ "post")) %>%
  drop_na(fase)
tet_f <- subset %>%
  filter(location %in% c("Tetiaroa"))%>%
  mutate(fase = case_when(year == 2005 ~ "before",
                          year %in% c(2013)~ "during",
                          year == 2019 ~ "post"))%>%
  drop_na(fase)

ggplot(mahe_f) +
  geom_boxplot(aes(x = fase, y = plank, color = location))

tet_c2 <- tet_c %>%
  group_by(year) %>%
  summarise(coral_cover = median(coral_cover))
tet_f2 <- tet_f %>%
  group_by(fase) %>%
  summarise_if(is.numeric, median)



tet_long <-tet_f %>%
  ungroup() %>%
  select(fase, exP, exNP, herb, plank, prod) %>%
  pivot_longer(-c(fase)) %>%
  filter(value>0) %>%
  #mutate_if(is.numeric, log) %>%
  group_by(name) %>%
  mutate(ref = median(value[fase == "before"])) %>%
  mutate(value_scaled = value/ref)

ggplot(tet_long) +
  geom_boxplot(aes(x = name, y = value_scaled, color = fase), outlier.alpha = 0)

ggplot(tet_long) +
  geom_boxplot(aes(x = fase, y = value_scaled, color = name), outlier.alpha = 0) 

library(brms)
fit_tet <- brm(value_scaled ~ fase*name, data = tet_long, 
                backend = "cmdstanr", family = "student")
summary(fit_tet)
#marginal_effects(fit_cousin)


nd <- select(tet_long, name, fase) %>% unique()
pred <- fitted(fit_tet, newdata = nd, summary = F)[1:1000,]
pred <- pred %>% as.data.frame() %>% pivot_longer(1:15) %>%
  mutate(fun = rep(nd$name, 1000),
         fase = rep(nd$fase, 1000))

pred <- cbind(nd, fitted(fit_tet, newdata = nd, summary = T)) %>%
  group_by(name) %>%
  mutate(ref = Estimate[fase == "before"]) 

ptet_cor <- ggplot(tet_c2) +
  geom_smooth(aes(x = year, y = (coral_cover)), color = "black") +
  geom_point(aes(x = year, y = (coral_cover)), data = tet_c, 
             position = position_jitter(width = 1)) +
  theme_classic() +
  scale_y_continuous(trans = "sqrt") +
  labs(y = "Coral cover", x = "Year")


plot_tet <- 
  ggplot(pred) +
  geom_hline(yintercept = 1, linetype = 3)  +
  geom_point(aes(x = fase, y = Estimate/ref, color = name),
             position = position_dodge(width = 0.5), size = 3) +
  geom_line(aes(x = as.numeric(as.factor(fase)), y = Estimate/ref, color = name,
                linetype = name),
            position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = `Q2.5`/ref, ymax = `Q97.5`/ref, x = fase, color = name),
                alpha = 0.6,  position = position_dodge(width = 0.5), width = 0, size = 2) +
  
  theme_bw() +
  scale_color_fish_d() +
  scale_fill_fish_d() + 
  labs(x = "Fase", y = "Scaled function", color = "Function", linetype = "Function",
       title = "Tetiaroa (Recovered)") +
  theme(text = element_text(size = 14)) +
  scale_y_continuous(limits = c(-0.1,2))
plot_tet

library(patchwork)
 ptet_cor + plot_tet + theme(legend.position = "right") +
   plot_layout(ncol = 1, heights = c(1,4))

ggsave("output/ts_tet.png", width = 5, height = 5)


#############" Moorea #########""
moo_c <- corcov %>%
  filter(location %in% c("Moorea")) %>%
  mutate(fase = case_when(year == 2005 ~ "before",
                          year %in% c(2011)~ "during",
                          year == 2017 ~ "post")) %>%
  drop_na(fase)
moo_f <- subset %>%
  filter(location %in% c("Moorea"))%>%
  mutate(fase = case_when(year == 2005 ~ "before",
                          year %in% c(2011)~ "during",
                          year == 2017 ~ "post"))%>%
  drop_na(fase)

ggplot(moo_f) +
  geom_boxplot(aes(x = fase, y = coral_cover, color = location))

moo_c2 <- moo_c %>%
  group_by(year) %>%
  summarise(coral_cover = median(coral_cover))
moo_f2 <- moo_f %>%
  group_by(fase) %>%
  summarise_if(is.numeric, median)



moo_long <-moo_f %>%
  ungroup() %>%
  select(fase, exP, exNP, herb, plank, prod) %>%
  pivot_longer(-c(fase)) %>%
  filter(value>0) %>%
  #mutate_if(is.numeric, log) %>%
  group_by(name) %>%
  mutate(ref = median(value[fase == "before"])) %>%
  mutate(value_scaled = value/ref)


library(brms)
fit_moo <- brm(value_scaled ~ fase*name, data = moo_long, 
               backend = "cmdstanr", family = "student")
summary(fit_moo)
#marginal_effects(fit_cousin)


nd <- select(moo_long, name, fase) %>% unique()
pred <- fitted(fit_moo, newdata = nd, summary = F)[1:1000,]
pred <- pred %>% as.data.frame() %>% pivot_longer(1:15) %>%
  mutate(fun = rep(nd$name, 1000),
         fase = rep(nd$fase, 1000))

pred <- cbind(nd, fitted(fit_moo, newdata = nd, summary = T)) %>%
  group_by(name) %>%
  mutate(ref = Estimate[fase == "before"]) 

pmoo_cor <- ggplot(moo_c2) +
  geom_smooth(aes(x = year, y = (coral_cover)), color = "black") +
  geom_point(aes(x = year, y = (coral_cover)), data = moo_c, 
             position = position_jitter(width = 1)) +
  theme_classic() +
  scale_y_continuous(trans = "sqrt") +
  labs(y = "Coral cover", x = "Year")


plot_moo <- 
  ggplot(pred) +
   geom_hline(yintercept = 1, linetype = 3)  +
  geom_point(aes(x = fase, y = Estimate/ref, color = name),
             position = position_dodge(width = 0.5), size = 3) +
  geom_line(aes(x = as.numeric(as.factor(fase)), y = Estimate/ref, color = name,
                linetype = name),
            position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = `Q2.5`/ref, ymax = `Q97.5`/ref, x = fase, color = name),
                alpha = 0.6,  position = position_dodge(width = 0.5), width = 0, size = 2) +
  
  theme_bw() +
  scale_color_fish_d() +
  scale_fill_fish_d() + 
  labs(x = "Fase", y = "Scaled function", color = "Function", linetype = "Function",
       title = "Moorea (Recovered)") +
  theme(text = element_text(size = 14)) +
  scale_y_continuous(limits = c(-0.1,3))
plot_moo

library(patchwork)
plot_poly <- ptet_cor +pmoo_cor + plot_tet +  theme(legend.position = "none") + 
  plot_moo + theme(legend.position = "right") +
  plot_layout(ncol = 2, heights = c(1,4))
plot_poly
ggsave("output/ts_poly.png", width = 10, height = 5)


 pmahe_cor +  pcous_cor + plot_mahe  +theme(legend.position = "none") + 
  plot_cousin  + ptet_cor +pmoo_cor + plot_tet +  theme(legend.position = "none") + 
   plot_moo + theme(legend.position = "none") +
   plot_layout(ncol = 2, heights = c(1,4, 1, 4))

 
 ggsave("output/ts_cases.png", width = 10, height = 10)
 