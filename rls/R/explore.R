targets::tar_load(transect_fish_benthic)
targets::tar_load(rls_predictors)
library(tidyverse)
library(fishualize)

comb <- transect_fish_benthic %>%
  left_join(select(rls_predictors, SiteCode, Realm, min_sst_5year, max_sst_5year, median_sst_1year, gravtot3)) %>%
  mutate(sst = median_sst_1year) %>%
  drop_na(sst) %>%
  mutate(check_sst = 28 > min_sst_5year | 28 < max_sst_5year) %>%
  drop_na(sst) %>%
  group_by(Realm) %>%
  mutate(pristine = case_when(gravtot3 < quantile(gravtot3, 0.10) ~ 1, TRUE ~ 0),
         impact = case_when(gravtot3 > quantile(gravtot3, 0.9) ~ 1, TRUE ~ 0))
  # group_by(SiteCode) %>%
  # summarise_if(is.numeric, median, na.rm = T)


ggplot(comb) +
  geom_point(aes(x = pc1, y = pc2)) +
  geom_point(aes(x = pc1, y = pc2), color = "red", 
             data = comb[comb$herb/comb$biom > quantile(comb$herb/comb$biom, 0.9),]) +
  scale_color_fish()

sum(comb$check_sst)/nrow(comb)

ggplot(comb) +
  geom_point(aes(x = sqrt(turf + algae + cca), y = log(plank)))  +
  geom_smooth(aes(x = sqrt(turf + algae + cca), y = log(plank)), method = "lm")
ggplot(comb) +
  geom_point(aes(x = sqrt(turf + algae + cca), y = log(plank)))  +
  geom_smooth(aes(x = sqrt(turf + algae + cca), y = log(plank)), method = "lm")
  

ggplot(comb) +
  geom_point(aes(x = sqrt(rock), y = log(herb/biom)), method = "lm") +
  geom_smooth(aes(x = sqrt(rock), y = log(herb/biom)), method = "lm")

comb_sel <- comb %>%
  filter(herb>0, plank>0, Accretion_kgyr>0) %>%
  mutate(plank = residuals(lm(log(plank) ~  log(biom) + sst)), 
         herb = residuals(lm(log(herb) ~   log(biom) + sst)), 
         prod = residuals(lm(log(prod) ~  log(biom) + sst)),  
         pisc = residuals(lm(log(pisc) ~ log(biom) +  sst)),  
         exP = residuals(lm(log(exP) ~   log(biom) + sst)), 
         exN = residuals(lm(log(exN) ~  log(biom) + sst))) %>%
  mutate(calc_log = log(Calc_ghr*100/(100-(left + sand))), 
         orgg_log = log(OrgGrowth_gyr*100/(100-(left + sand))), 
         gpp_log = log(GPP_ghr*100/(100-(left + sand))),  
         accr_log = log(Accretion_kgyr*100/(100-(left + sand))), 
         space = (BranchSpace_cm3*100/(100-(left + sand))), 
         sbiom_log = log( BiomassStand_g*100/(100-(left + sand))),
         rug_log = log(Rugosity*100/(100-(left + sand))),
         turn = log(turn),  outNP = log(outNP)) 

comb_fish <- comb_sel %>%
  ungroup() %>%
  select(herb, prod, plank, exP, exNP) %>%
  drop_na(herb, prod, plank, exP, exNP)

comb_ben <- comb %>%
  ungroup() %>%
  select(BiomassStand_g, Calc_ghr, BranchSpace_cm3, Rugosity,
         GPP_ghr) %>%
  drop_na()
library(ade4)
dudi1 <- dudi.pca(comb_fish, scale = TRUE, center = TRUE, scan = FALSE, nf = 2)
dudi2 <- dudi.pca(comb_ben, scale = TRUE, center = TRUE, scan = FALSE, nf = 2)
coin1 <- coinertia(dudi1,dudi2, scan = FALSE, nf = 2)
coin1$aY

summary(coin1)
plot(coin1)

plot(coin1$lX)
plot(dudi1$li)

library(GGally)
#ggpairs(comb_sel)

library(corrplot)
corrplot(cor(comb_sel))

comb1 <-  comb %>%
  filter(herb>0, plank>0, Accretion_kgyr>0) %>%
  mutate(plank = residuals(lm(log(plank) ~  log(biom) + sst)), 
         herb = residuals(lm(log(herb) ~   log(biom) + sst)), 
         prod = residuals(lm(log(prod) ~  log(biom) + sst)),  
         pisc = residuals(lm(log(pisc) ~ log(biom) +  sst)),  
         exP = residuals(lm(log(exP) ~   log(biom) + sst)), 
         exN = residuals(lm(log(exN) ~  log(biom) + sst))) %>%
  mutate(calc_log = log(Calc_ghr*100/(100-(left + sand))), 
         orgg_log = log(OrgGrowth_gyr*100/(100-(left + sand))), 
         gpp_log = log(GPP_ghr*100/(100-(left + sand))),  
         accr_log = log(Accretion_kgyr*100/(100-(left + sand))), 
         space = (BranchSpace_cm3*100/(100-(left + sand))), 
         sbiom_log = log( BiomassStand_g*100/(100-(left + sand))),
         rug_log = log(Rugosity*100/(100-(left + sand))),
         turn = log(turn),  outNP = log(outNP)) %>%
  drop_na(herb, prod, plank, exP, exNP)

ggplot(comb) +
  geom_point(aes(x = log(Accretion_kgyr), y = log(plank)))

ggplot(comb) +
  geom_smooth(aes(x = log(algae), y = (herb)))
summary(comb_sel)
pca <- prcomp(log1p(comb_ben), scale = T, center = T) 
summary(pca)
biplot(pca)

pc <- as.data.frame(pca$x) %>%
  cbind(comb)
rot <- as.data.frame(pca$rotation)

rot1 <- ggplot(rot) +
  geom_vline(aes(xintercept = 0)) +
  geom_hline(aes(yintercept = 0)) +
  geom_point(aes(x = PC1, y = PC2), data = rot) +
  geom_label(aes(x = PC1, y = PC2, label = rownames(rot)), size = 5) +
  coord_equal() 
rot1

rot2 <- ggplot(rot) +
  geom_vline(aes(xintercept = 0)) +
  geom_hline(aes(yintercept = 0)) +
  geom_point(aes(x = PC3, y = PC2), data = rot) +
  geom_label(aes(x = PC3, y = PC2, label = rownames(rot)), size = 2) +
  coord_equal() 
rot2

library(patchwork)
library(fishualize)

ggplot(rot) +
  coord_equal() +
  geom_vline(aes(xintercept = 0)) +
  geom_hline(aes(yintercept = 0)) +
  geom_point(aes(x = PC1, y = PC3), data = rot) +
  geom_label(aes(x = PC1, y = PC3, label = rownames(rot))) +
  coord_equal() 



rot1 +
  ggplot(pc) +
  geom_point(aes(x = PC1, y = PC2, color = coral), data = pc, alpha = 0.4, size = 3) +
  geom_vline(aes(xintercept = 0)) +
  geom_hline(aes(yintercept = 0))+
  coord_equal() +
  xlim(c(-6,6)) +
  ylim(c(-6,6)) +
  scale_color_fish() +

  ggplot(pc, aes(x = PC1, y = PC2)) +
  geom_point(aes(x = PC1, y = PC2, color = Realm),
             data = pc, 
             alpha = 0.4, size = 2) +
  ggforce::geom_mark_hull(aes(x = PC1, y = PC2),
                          data = pc[pc$pristine == 1,],
                          concavity=100,radius=unit(.15,"cm"), color = "blue",
                          expand=unit(.15,"cm"),alpha=.15,size=1)+
  ggforce::geom_mark_hull(aes(x = PC1, y = PC2),
                          data = pc[pc$impact == 1,],
                          concavity=100,radius=unit(.15,"cm"), color = "red",
                          expand=unit(.15,"cm"),alpha=.15,size=1)+
  geom_vline(aes(xintercept = 0)) +
  geom_hline(aes(yintercept = 0))+
  coord_equal() +
  xlim(c(-6,6)) +
  ylim(c(-6,6)) +
  scale_color_fish_d() +
  
  ggplot(pc, aes(x = pc1, y = pc2)) +
  geom_point(aes(x = pc1, y = pc2),
             data = pc, 
             alpha = 0.4, size = 2) +
  
  ggforce::geom_mark_hull(aes(x = pc1, y = pc2),
                          data = pc[pc$gravtot3> quantile(pc$gravtot3, 0.95),],
                          concavity=100,radius=unit(.15,"cm"), color = "red",
                          expand=unit(.15,"cm"),alpha=.15,size=1)+
  ggforce::geom_mark_hull(aes(x = pc1, y = pc2),
                          data = pc[pc$gravtot3 < quantile(pc$gravtot3, 0.05),],
                          concavity=100,radius=unit(.15,"cm"), color = "blue",
                          expand=unit(.15,"cm"),alpha=.15,size=1)+
  geom_vline(aes(xintercept = 0)) +
  geom_hline(aes(yintercept = 0))+
  coord_equal() +
  scale_color_fish() & theme_bw()



comb <-  transect_fish_benthic %>%
  left_join(select(rls_predictors, SiteCode, min_sst_5year, max_sst_5year, median_sst_1year, gravtot3)) %>%
  mutate(sst = median_sst_1year) %>%
  drop_na(sst) %>%
  mutate(check_sst = 28 > min_sst_5year | 28 < max_sst_5year) %>%
  drop_na(sst) #%>%

ggplot(comb1) +
  geom_point(aes(x = pc1, y = pc2, color = plank)) +
  scale_color_fish()

ggplot(comb) +
  geom_point(aes(x = log(gravtot3), y = log(prod), color = (sst))) +
  geom_smooth(aes(x = log(gravtot3), y = log(prod)), method = "lm")  +
  scale_color_fish() +
  labs(x = "Human gravity", y = "log(Production)") +
  
  
  ggplot(comb) +
  geom_point(aes(x = log(gravtot3), y = log(turn), color = (sst))) +
  geom_smooth(aes(x = log(gravtot3), y = log(turn)), method = "lm") +
  scale_color_fish() +
  labs(x = "Human gravity", y = "log(turnover)")+
  
  
  ggplot(comb) +
  geom_point(aes(x = log(gravtot3), y = log(exNP), color = (sst))) +
  geom_smooth(aes(x = log(gravtot3), y = log(exNP)), method = "lm")  +
  scale_color_fish() +
  labs(x = "Human gravity", y = "log(flux N:P)")+
  
  ggplot(comb) +
  geom_point(aes(x = log(gravtot3), y = log(exP), color = (sst))) +
  geom_smooth(aes(x = log(gravtot3), y = log(exP)), method = "lm")  +
  scale_color_fish() +
  labs(x = "Human gravity", y = "log(P flux)")+
  
  ggplot(comb) +
  geom_point(aes(x = log(gravtot3), y = log(plank), color = sst)) +
  geom_smooth(aes(x = log(gravtot3), y = log(plank)), method = "lm")  +
  scale_color_fish() +
  labs(x = "Human gravity", y = "log(Planktivory) (gC/day)") +
  
  
  ggplot(comb) +
  geom_point(aes(x = log(gravtot3), y = log(herb), color = sst)) +
  geom_smooth(aes(x = log(gravtot3), y = log(herb)), method = "lm")  +
  scale_color_fish() +
  labs(x = "Human gravity", y = "log(Herbivory) (gC/day)") +
  
  
  ggplot(comb) +
  geom_point(aes(x = log(gravtot3), y = log(Calc_ghr), color = sst)) +
  geom_smooth(aes(x = log(gravtot3), y = log(Calc_ghr)), method = "lm")  +
  scale_color_fish() +
  labs(x = "Human gravity", y = "log(Calcification) (gC/day)") +
  
  
  
  ggplot(comb) +
  geom_point(aes(x = log(gravtot3), y = log(Accretion_kgyr), color = coral)) +
  geom_smooth(aes(x = log(gravtot3), y = log(Accretion_kgyr)), method = "lm")  +
  scale_color_fish() +
  labs(x = "Human gravity", y = "log(Accretion_kgyr)") +
  
  ggplot(comb) +
  geom_point(aes(x = log(gravtot3), y = log(GPP_ghr), color = coral)) +
  geom_smooth(aes(x = log(gravtot3), y = log(GPP_ghr)), method = "lm")  +
  scale_color_fish() +
  labs(x = "Human gravity", y = "log(Gpp)") +
  
  
  ggplot(comb) +
  geom_point(aes(x = log(gravtot3), y = log(BiomassStand_g), color = coral)) +
  geom_smooth(aes(x = log(gravtot3), y = log(BiomassStand_g)), method = "lm")  +
  scale_color_fish() +
  labs(x = "Human gravity", y = "log(BiomassStand_g)") 
  
  

ggplot(comb) +
  geom_point(aes(x = log(coral), y = log(prod), color = (sst))) +
  geom_smooth(aes(x = log(coral), y = log(prod)), method = "lm")  +
  scale_color_fish() +
  labs(x = "Human gravity", y = "log(Production)") +
  
  
  ggplot(comb) +
  geom_point(aes(x = log(coral), y = log(turn), color = (sst))) +
  geom_smooth(aes(x = log(coral), y = log(turn)), method = "lm") +
  scale_color_fish() +
  labs(x = "Human gravity", y = "log(turnover)")+
  
  
  ggplot(comb) +
  geom_point(aes(x = log(coral), y = log(exNP), color = (sst))) +
  geom_smooth(aes(x = log(coral), y = log(exNP)), method = "lm")  +
  scale_color_fish() +
  labs(x = "Human gravity", y = "log(flux N:P)")+
  
  ggplot(comb) +
  geom_point(aes(x = log(coral), y = log(exP), color = (sst))) +
  geom_smooth(aes(x = log(coral), y = log(exP)), method = "lm")  +
  scale_color_fish() +
  labs(x = "Human gravity", y = "log(P flux)")+
  
  ggplot(comb) +
  geom_point(aes(x = log(coral), y = log(plank), color = sst)) +
  geom_smooth(aes(x = log(coral), y = log(plank)), method = "lm")  +
  scale_color_fish() +
  labs(x = "Human gravity", y = "log(Planktivory) (gC/day)") +
  
  
  ggplot(comb) +
  geom_point(aes(x = log(coral), y = log(herb), color = sst)) +
  geom_smooth(aes(x = log(coral), y = log(herb)), method = "lm")  +
  scale_color_fish() +
  labs(x = "Human gravity", y = "log(Herbivory) (gC/day)") +
  
  
  ggplot(comb) +
  geom_point(aes(x = log(coral), y = log(Calc_ghr), color = sst)) +
  geom_smooth(aes(x = log(coral), y = log(Calc_ghr)), method = "lm")  +
  scale_color_fish() +
  labs(x = "Human gravity", y = "log(Calcification) (gC/day)") +
  
  
  
  ggplot(comb) +
  geom_point(aes(x = log(coral), y = log(Accretion_kgyr), color = coral)) +
  geom_smooth(aes(x = log(coral), y = log(Accretion_kgyr)), method = "lm")  +
  scale_color_fish() +
  labs(x = "Human gravity", y = "log(Accretion_kgyr)") +
  
  ggplot(comb) +
  geom_point(aes(x = log(coral), y = log(GPP_ghr), color = coral)) +
  geom_smooth(aes(x = log(coral), y = log(GPP_ghr)), method = "lm")  +
  scale_color_fish() +
  labs(x = "Human gravity", y = "log(Gpp)") +
  
  
  ggplot(comb) +
  geom_point(aes(x = log(coral), y = log(BiomassStand_g), color = coral)) +
  geom_smooth(aes(x = log(coral), y = log(BiomassStand_g)), method = "lm")  +
  scale_color_fish() +
  labs(x = "Human gravity", y = "log(BiomassStand_g)") 

  # ggplot(comb) +
  # geom_point(aes(x = sqrt(coral), y = log(prod), color = (sst))) +
  # geom_smooth(aes(x = sqrt(coral), y = log(prod)), method = "lm")  +
  # scale_color_fish() +
  # labs(x = "Coral", y = "log(Production)") +
  # 
  # 
  # ggplot(comb) +
  # geom_point(aes(x = sqrt(coral), y = log(turn), color = (sst))) +
  # geom_smooth(aes(x = sqrt(coral), y = log(turn)), method = "lm") +
  # scale_color_fish() +
  # labs(x = "Coral", y = "log(turnover)")+
  # 
  # 
  # ggplot(comb) +
  # geom_point(aes(x = sqrt(coral), y = log(exNP), color = (sst))) +
  # geom_smooth(aes(x = sqrt(coral), y = log(exNP)), method = "lm")  +
  # scale_color_fish() +
  # labs(x = "Coral", y = "log(flux N:P)")+
  # 
  # ggplot(comb) +
  # geom_point(aes(x = sqrt(coral), y = log(exP), color = (sst))) +
  # geom_smooth(aes(x = sqrt(coral), y = log(exP)), method = "lm")  +
  # scale_color_fish() +
  # labs(x = "Coral", y = "log(P flux)")+
  # 
  # ggplot(comb) +
  # geom_point(aes(x = sqrt(coral), y = log(plank), color = sst)) +
  # geom_smooth(aes(x = sqrt(coral), y = log(plank)), method = "lm")  +
  # scale_color_fish() +
  # labs(x = "Coral", y = "log(Planktivory) (gC/day)") +
  # 
  # 
  # ggplot(comb) +
  # geom_point(aes(x = sqrt(coral), y = log(herb), color = sst)) +
  # geom_smooth(aes(x = sqrt(coral), y = log(herb)), method = "lm")  +
  # scale_color_fish() +
  # labs(x = "Coral", y = "log(Herbivory) (gC/day)")

###### pca ipfc #########
cols <- function(){
  c("#00183aff","#c2c2c2ff","#e6e6e6ff","#6ac2e5ff",
    "#85ffbcff","#fffc5cff","#f6ab13ff","#df3416ff")
}
theme_ppt <- function(){
  theme_classic() +
    theme(rect = element_rect(fill = cols()[1], color = cols()[1]), 
          panel.background = element_rect(fill = cols()[1], color = cols()[1]), 
          plot.background = element_rect(fill = cols()[1], color = cols()[1]), 
          text = element_text(color = cols()[3], size = 18),
          axis.text = element_text(color = cols()[3], size = 14),
          axis.ticks = element_line(color = cols()[3]), 
          axis.line = element_line(colour = cols()[3], size = 1) 
    ) 
}


sub <- comb %>%
  ungroup() %>%
  dplyr::select(exN, exP, plank, herb, prod, turn, pisc, 
         GPP_ghr, OrgGrowth_gyr, Calc_ghr, Storage_kg, Rugosity, BranchSpace_cm3, Accretion_kgyr) %>%
  mutate_all(function(x){(log10(x + 1))})

# pca
pca <- prcomp(sub, scale = T, center = T) 

biplot(pca)

pc <- as.data.frame(pca$x) %>%
  cbind(comb)

rot <- as.data.frame(pca$rotation)


# plotting
rot1 <- ggplot(rot) +
  geom_point(aes(x = PC1, y = PC2),
             data = pc, 
             alpha = 0.4, size = 1, color = cols()[2]) +
  geom_vline(aes(xintercept = 0), color = "grey90") +
  geom_hline(aes(yintercept = 0), color = "grey90") +
  geom_segment(aes(x = 0, y = 0, xend = PC1*11, yend = PC2*11, color = rownames(rot)), 
               size = 1, data = rot) +
  geom_label(aes(x = PC1*11, y = PC2*11, color = rownames(rot),
                 label = c("ExN", "ExP", "Plank", "Herb", "Prod" ,
                           "Turn" , "Pisc" ,"GPP", "OrgGr", "Calc", "Stor", "Rug", 
                           "Space", "Accr" )), fill = "#00183aff",
size = 3) +
  coord_equal() +
  scale_color_manual(values = c("#BF961BFF", "#997100FF",   "#BE7F36FF","#003C74FF", "#003D64FF", "#F19564FF",
                                "#086986FF",  "#FA885BFF","#06A2BCFF", "#00C5E6FF", "#3FDBFBFF",  "#E8673BFF", "#D94B2BFF", "#A0F5F7FF")) +
  xlim(c(-6,6)) + ylim(c(-6,5)) +
  theme_ppt() +
  theme(legend.position = "none")
rot1

plot <- 
  rot1 +
  
  ggplot(comb) +
  geom_point(aes(x = coral, y = biom/1000), color = "grey90", alpha = 0.3) +
  geom_hline(yintercept = exp(mean(log(comb$biom/1000))), color = "grey95") +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(trans = "sqrt") +
  theme_ppt() +
  theme(axis.text = element_text(size = 10)) +
  labs(x = "Coral cover (%)", y = "Fish biomass (kg)")
plot  


  ggsave(filename = "plot_ppt.png", plot = plot, 
         width = 27, height = 12,
         units = "cm")

  
  
  plot <- ggplot(rot) +
    geom_point(aes(x = PC1, y = PC2, color = (gravtot3)),
               data = pc, 
               alpha = 0.4, size = 1) +
    geom_vline(aes(xintercept = 0), color = "grey90") +
    geom_hline(aes(yintercept = 0), color = "grey90") +
       xlim(c(-6,6)) + ylim(c(-6,5)) +
    scale_color_fish(option = "Trimma_lantana", trans = "log10") +
    theme_ppt() +
    labs(color = "Human gravity") +
    theme(legend.position = "bottom", legend.text = element_text(size = 12)) +
    guides(color = guide_colorbar(title.position = 'top', title.hjust = .5,
                                 barwidth = unit(15, 'lines'), 
                                 barheight = unit(.7, 'lines')))
  plot

  
  ggsave(filename = "plot_ppt_2.png", plot = plot, 
         width = 12, height = 12,
         units = "cm")
  
  fish(option = "Trimma_lantana", n = 5)
  
  plot <- ggplot(rot) +
    geom_point(aes(x = PC1, y = PC2, color = (gravtot3)),
               data = pc, 
               alpha = 0.4, size = 1) +
    geom_vline(aes(xintercept = 0), color = "grey90") +
    geom_hline(aes(yintercept = 0), color = "grey90") +
    xlim(c(-6,6)) + ylim(c(-6,5)) +
    scale_color_fish(option = "Trimma_lantana", trans = "log10") +
    theme_ppt() +
    labs(color = "Human gravity") +
    ggforce::geom_mark_hull(aes(x = PC1, y = PC2),
                            data = pc[pc$pristine == 1,],
                            concavity=100,radius=unit(.15,"cm"), color = "#4AB2B8FF" ,
                            expand=unit(.15,"cm"),alpha=.15,size=1)+
    ggforce::geom_mark_hull(aes(x = PC1, y = PC2),
                            data = pc[pc$impact == 1,],
                            concavity=100,radius=unit(.15,"cm"), color =  "#ED1C24FF",
                            expand=unit(.15,"cm"),alpha=.15,size=1)+
    
    theme(legend.position = "bottom", legend.text = element_text(size = 12)) +
    guides(color = guide_colorbar(title.position = 'top', title.hjust = .5,
                                  barwidth = unit(15, 'lines'), 
                                  barheight = unit(.7, 'lines')))
  plot
  
  
  ggsave(filename = "plot_ppt_3.png", plot = plot, 
         width = 12, height = 12,
         units = "cm")
  
  
  
  


