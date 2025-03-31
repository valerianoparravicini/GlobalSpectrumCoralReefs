#file load

rls_benthic_functions <- read_csv("../2_rls_benthos/output/RLS_benthic_functions.csv")

fish_functions <- read_csv("../1_rls_fish/output/transect_fish_functions.csv")

list_dataset_functions <- combine_benthic_fish_functions(rls_benthic_functions, fish_functions) 

list_pca_functions <- pca_functions(list_dataset_functions, benth=rls_benthic_functions, fish=fish_functions  )

figure_2(list_pca_functions , benth=rls_benthic_functions, fish=fish_functions  )

pred <- read_csv("data/RLS_predictors.csv")

#functions

combine_benthic_fish_functions <- function(rls_benthic_functions, fish_functions) {
  
  
  dataset_functions <- merge(rls_benthic_functions, fish_functions, by = "SurveyID")
  
  rls_functions <- dataset_functions 
  
  # Building dataframes
  dataset_functions <- data.frame(Survey_ID = rep(dataset_functions$SurveyID, 14),
                                  functions = c(dataset_functions$herb, dataset_functions$plank, dataset_functions$prod,
                                                dataset_functions$pisc, dataset_functions$exN, dataset_functions$exP,
                                                dataset_functions$turn, dataset_functions$GPP_ghr, dataset_functions$OrgGrowth_gyr,
                                                dataset_functions$Calc_ghr, dataset_functions$Storage_kg, dataset_functions$Accretion_kgyr,
                                                dataset_functions$Rugosity, dataset_functions$BranchSpace_cm3),
                                  Label = c(rep("HERBIVORY \n", length(dataset_functions$SurveyID)),
                                            rep("PLANKTIVORY \n", length(dataset_functions$SurveyID)),
                                            rep("FISH \nPRODUCTION", length(dataset_functions$SurveyID)),
                                            rep("PISCIVORY \n", length(dataset_functions$SurveyID)),
                                            rep("NITROGEN \nEXCRETION", length(dataset_functions$SurveyID)),
                                            rep("PHOSPHORUS \nEXCRETION", length(dataset_functions$SurveyID)),
                                            rep("FISH \nTURNOVER", length(dataset_functions$SurveyID)),
                                            rep("GROSS PRIMARY \nPRODUCTION", length(dataset_functions$SurveyID)),
                                            rep("ORGANIC \nGROWTH", length(dataset_functions$SurveyID)),
                                            rep("CALCIFICATION \n", length(dataset_functions$SurveyID)),
                                            rep("CARBONATE \nSTORAGE", length(dataset_functions$SurveyID)),
                                            rep("INORGANIC \nGROWTH", length(dataset_functions$SurveyID)),
                                            rep("RUGOSITY \n", length(dataset_functions$SurveyID)),
                                            rep("BRANCH \nSPACING", length(dataset_functions$SurveyID))),
                                  colors = c("#003C74FF", "#003D64FF", "#086986FF", "#06A2BCFF", "#00C5E6FF", "#3FDBFBFF", "#A0F5F7FF",
                                             "#BF961BFF", "#997100FF", "#BE7F36FF", "#F19564FF", "#FA885BFF", "#E8673BFF", "#D94B2BFF"))
  
  dataset_functions$Label <- factor(dataset_functions$Label, levels = c("HERBIVORY \n", "PLANKTIVORY \n", "FISH \nPRODUCTION", 
                                                                        "PISCIVORY \n", "NITROGEN \nEXCRETION", "PHOSPHORUS \nEXCRETION", 
                                                                        "FISH \nTURNOVER", "GROSS PRIMARY \nPRODUCTION", "ORGANIC \nGROWTH",
                                                                        "CALCIFICATION \n", "CARBONATE \nSTORAGE", "INORGANIC \nGROWTH", 
                                                                        "RUGOSITY \n", "BRANCH \nSPACING"))
  
  dataset_functions <- dataset_functions %>% mutate(functions_log10 = (functions+0.000001),
                                                    functions_log10_rounded = (functions+0.000001))
  
  dataset_functions %>% group_by(Label) %>% summarise (min = min((functions_log10)), max = max((functions_log10)))
  
  
 save(dataset_functions, file="output/dataset_functions.RData")
 save(rls_functions, file="output/rls_functions.RData")
  
}




figure_1 <- function(dataset_functions, rls_functions) {
  
  
  # Setting map parameters ####
  # warning messages due to an update from PROJ4 to PROJ6. In the meantime, the warning is just a nuisance and has no implications
  load(url("https://github.com/valentinitnelav/RandomScripts/blob/master/NaturalEarth.RData?raw=true"))
  PROJ <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" 
  NE_countries_rob  <- spTransform(NE_countries, CRSobj = PROJ) 
  NE_graticules_rob <- spTransform(NE_graticules, CRS(PROJ))
  NE_box_rob <- spTransform(NE_box, CRSobj = PROJ) 
  spatial_points_y <- SpatialPoints(cbind(lbl.Y$lon, lbl.Y$lat), proj4string = CRS("+proj=longlat +datum=WGS84"))
  prj.coord <- spTransform(spatial_points_y, CRS(PROJ))
  lbl.Y.prj <- cbind(prj.coord, lbl.Y) %>% data.frame() ; names(lbl.Y.prj)[1:2] <- c("X.prj","Y.prj") 
  spatial_points_x <- SpatialPoints(cbind(lbl.X$lon, lbl.X$lat), proj4string = CRS("+proj=longlat +datum=WGS84"))
  prj.coord <- spTransform(spatial_points_x, CRS(PROJ))
  lbl.X.prj <- cbind(prj.coord, lbl.X) %>% data.frame() ; names(lbl.X.prj)[1:2] <- c("X.prj","Y.prj") 
  lbl.X.prj <- lbl.X.prj[c(seq(1, nrow(lbl.X.prj), by = 2)), ]
  
  rls <- rls_functions
  
  spatial_points_rls <- SpatialPoints(cbind(rls$SiteLongitude , rls$SiteLatitude), proj4string = CRS("+proj=longlat +datum=WGS84"))
  data_Spatial_proj = spTransform(spatial_points_rls, CRS(PROJ)) %>% 
    as.data.frame() %>% rename(Long_Robin = coords.x1, Lat_Robin = coords.x2) 
  rls = rls %>% cbind(.,data_Spatial_proj)
  
  
  Figure_1A = ggplot(data = rls, aes(x = Long_Robin, y = Lat_Robin)) +
    geom_path(data=NE_graticules_rob, aes(long, lat, group=group), linetype="dotted", color="grey40", size = 0.25) +  
    geom_polygon(data=NE_countries_rob, aes(long,lat, group=group), colour="gray80", fill="gray80", linewidth = 0.25) +
      geom_polygon(data=NE_box_rob, aes(x=long, y=lat), colour="grey80", fill="transparent", size = 0.25) +
    
      
      
      geom_text(data = lbl.Y.prj, aes(x = coords.x1, y = coords.x2, label = lbl), color="grey80", size=2) +
      geom_text(data = lbl.X.prj, aes(x = coords.x1, y = coords.x2, label = lbl), color="grey80", size=2) +
      geom_point(show.legend = F, size = 3, fill = "grey40", shape = 21) +
      coord_sf(ylim = c(-4000000,4000000), expand = F) +
    
    annotate(geom = 'segment', 
             y = -4950000, 
             yend = -4950000, 
             x = -15150000, 
             xend = 15150000, 
             size = 0.25,
             color="grey80") +
    
    annotate(geom = 'segment', 
             y = 4950000, 
             yend = 4950000, 
             x = -15150000, 
             xend = 15150000, 
             size = 0.25,
             color="grey80") +
    
    theme(axis.line.x = element_line(size = 1)) +
      
      theme_void() + theme(legend.position="bottom")
  
  
                            
  
  
  # Figure 1B
  (Figure_1B <- (ggplot(data = dataset_functions, aes(x=functions_log10, fill = Label)) +
                   geom_histogram(color="black", linewidth = .1, bins = 15) + 
                   facet_wrap(~Label, ncol = 7, scales = "free_x") +
                   scale_x_log10(name = "", labels = label_comma()) +
                   scale_fill_manual(values = dataset_functions$colors) + theme_classic() +
                   theme(legend.position = "none", strip.text = element_text(hjust = 0),
                         strip.background = element_blank(),
                         panel.background = element_blank(),
                         panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(),
                         axis.line = element_line(colour = "black")) + scale_y_continuous(name = "number of reefs")))
  
  dataset_functions_split <- dataset_functions %>% group_by(Label) %>% group_split()
  
  Fig_1B = vector("list", length = 14)
  K = 1
  (Fig_1B[[K]] <- ggplot(data = dataset_functions_split[[K]], aes(x=functions_log10, fill = Label)) +
      ggtitle(unique(dataset_functions_split[[K]]$Label)) + 
      geom_histogram(color="black", linewidth = .1, bins = 15) + 
      scale_x_log10(name = expression(paste("g."* h^-1)), breaks = c(0.1, 1000)) +
      scale_y_continuous(name = "", limits = c(0,400), breaks = c(0, 300), expand = c(0,0)) +
      scale_fill_manual(values = dataset_functions_split[[K]]$colors) + theme_classic() +
      theme(legend.position = "none", strip.text = element_text(hjust = 0),
            strip.background = element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 8),
            axis.title.x = element_text(size = 10),
            plot.margin = unit(c(0.5,0.5,1,0.5), "cm"),
            plot.title = element_text(size = 10, face = "bold"),
            axis.line = element_line(colour = "black")))
  
  K = 2
  (Fig_1B[[K]] <- ggplot(data = dataset_functions_split[[K]], aes(x=functions_log10, fill = Label)) +
      ggtitle(unique(dataset_functions_split[[K]]$Label)) + 
      geom_histogram(color="black", linewidth = .1, bins = 15) + 
      scale_x_log10(name = expression(paste("g."* h^-1)), breaks = c(0.1, 1000)) +
      scale_y_continuous(name = "", limits = c(0,400), breaks = c(0, 300), expand = c(0,0)) +
      scale_fill_manual(values = dataset_functions_split[[K]]$colors) + theme_classic() +
      theme(legend.position = "none", strip.text = element_text(hjust = 0),
            strip.background = element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 8),
            axis.title.x = element_text(size = 10),
            plot.margin = unit(c(0.5,0.5,1,0.5), "cm"),
            plot.title = element_text(size = 10, face = "bold"),
            axis.line = element_line(colour = "black")))
  
  K = 3
  (Fig_1B[[K]] <- ggplot(data = dataset_functions_split[[K]], aes(x=functions_log10, fill = Label)) +
      ggtitle(unique(dataset_functions_split[[K]]$Label)) + 
      geom_histogram(color="black", linewidth = .1, bins = 15) + 
      scale_x_log10(name = expression(paste("kg."* yr^-1)), breaks = c(0.1, 100)) +
      scale_y_continuous(name = "", limits = c(0,400), breaks = c(0, 300), expand = c(0,0)) +
      scale_fill_manual(values = dataset_functions_split[[K]]$colors) + theme_classic() +
      theme(legend.position = "none", strip.text = element_text(hjust = 0),
            strip.background = element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 8),
            axis.title.x = element_text(size = 10),
            plot.margin = unit(c(0.5,0.5,1,0.5), "cm"),
            plot.title = element_text(size = 10, face = "bold"),
            axis.line = element_line(colour = "black")))
  
  K = 4
  (Fig_1B[[K]] <- ggplot(data = dataset_functions_split[[K]], aes(x=functions_log10, fill = Label)) +
      ggtitle(unique(dataset_functions_split[[K]]$Label)) + 
      geom_histogram(color="black", linewidth = .1, bins = 15) + 
      scale_x_log10(name = expression(paste("g."* h^-1)), breaks = c(0.1, 1000)) +
      scale_y_continuous(name = "", limits = c(0,400), breaks = c(0, 300), expand = c(0,0)) +
      scale_fill_manual(values = dataset_functions_split[[K]]$colors) + theme_classic() +
      theme(legend.position = "none", strip.text = element_text(hjust = 0),
            strip.background = element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 8),
            axis.title.x = element_text(size = 10),
            plot.margin = unit(c(0.5,0.5,1,0.5), "cm"),
            plot.title = element_text(size = 10, face = "bold"),
            axis.line = element_line(colour = "black")))
  
  K = 5
  (Fig_1B[[K]] <- ggplot(data = dataset_functions_split[[K]], aes(x=functions_log10, fill = Label)) +
      ggtitle(unique(dataset_functions_split[[K]]$Label)) + 
      geom_histogram(color="black", linewidth = .1, bins = 15) + 
      scale_x_log10(name = expression(paste("g."* h^-1)), breaks = c(0.1, 1000)) +
      scale_y_continuous(name = "", limits = c(0,400), breaks = c(0, 300), expand = c(0,0)) +
      scale_fill_manual(values = dataset_functions_split[[K]]$colors) + theme_classic() +
      theme(legend.position = "none", strip.text = element_text(hjust = 0),
            strip.background = element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 8),
            axis.title.x = element_text(size = 10),
            plot.margin = unit(c(0.5,0.5,1,0.5), "cm"),
            plot.title = element_text(size = 10, face = "bold"),
            axis.line = element_line(colour = "black")))
  
  K = 6
  (Fig_1B[[K]] <- ggplot(data = dataset_functions_split[[K]], aes(x=functions_log10, fill = Label)) +
      ggtitle(unique(dataset_functions_split[[K]]$Label)) + 
      geom_histogram(color="black", linewidth = .1, bins = 15) + 
      scale_x_log10(name = expression(paste("g."* h^-1)), breaks = c(0.1, 100)) +
      scale_y_continuous(name = "", limits = c(0,400), breaks = c(0, 300), expand = c(0,0)) +
      scale_fill_manual(values = dataset_functions_split[[K]]$colors) + theme_classic() +
      theme(legend.position = "none", strip.text = element_text(hjust = 0),
            strip.background = element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 8),
            axis.title.x = element_text(size = 10),
            plot.margin = unit(c(0.5,0.5,1,0.5), "cm"),
            plot.title = element_text(size = 10, face = "bold"),
            axis.line = element_line(colour = "black")))
  
  K = 7
  (Fig_1B[[K]] <- ggplot(data = dataset_functions_split[[K]], aes(x=functions_log10, fill = Label)) +
      ggtitle(unique(dataset_functions_split[[K]]$Label)) + 
      geom_histogram(color="black", linewidth = .1, bins = 15) + 
      scale_x_log10(name = expression(paste("g."* h^-1)), breaks = c(0.0001, 0.001)) +
      scale_y_continuous(name = "", limits = c(0,400), breaks = c(0, 300), expand = c(0,0)) +
      scale_fill_manual(values = dataset_functions_split[[K]]$colors) + theme_classic() +
      theme(legend.position = "none", strip.text = element_text(hjust = 0),
            strip.background = element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 8),
            axis.title.x = element_text(size = 10),
            plot.margin = unit(c(0.5,0.5,1,0.5), "cm"),
            plot.title = element_text(size = 10, face = "bold"),
            axis.line = element_line(colour = "black")))
  
  K = 8
  (Fig_1B[[K]] <- ggplot(data = dataset_functions_split[[K]], aes(x=functions_log10, fill = Label)) +
      ggtitle(unique(dataset_functions_split[[K]]$Label)) + 
      geom_histogram(color="black", linewidth = .1, bins = 15) + 
      scale_x_log10(name = expression(paste("g."* h^-1)), breaks = c(0.01, 10)) +
      scale_y_continuous(name = "", limits = c(0,400), breaks = c(0, 300), expand = c(0,0)) +
      scale_fill_manual(values = dataset_functions_split[[K]]$colors) + theme_classic() +
      theme(legend.position = "none", strip.text = element_text(hjust = 0),
            strip.background = element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 8),
            axis.title.x = element_text(size = 10),
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
            plot.title = element_text(size = 10, face = "bold"),
            axis.line = element_line(colour = "black")))
  
  K = 9
  (Fig_1B[[K]] <- ggplot(data = dataset_functions_split[[K]], aes(x=functions_log10, fill = Label)) +
      ggtitle(unique(dataset_functions_split[[K]]$Label)) + 
      geom_histogram(color="black", linewidth = .1, bins = 15) + 
      scale_x_log10(name = expression(paste("kg." * yr^-1)), breaks = c(10, 10000)) +
      scale_y_continuous(name = "", limits = c(0,400), breaks = c(0, 300), expand = c(0,0)) +
      scale_fill_manual(values = dataset_functions_split[[K]]$colors) + theme_classic() +
      theme(legend.position = "none", strip.text = element_text(hjust = 0),
            strip.background = element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 8),
            axis.title.x = element_text(size = 10),
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
            plot.title = element_text(size = 10, face = "bold"),
            axis.line = element_line(colour = "black")))
  
  K = 10
  (Fig_1B[[K]] <- ggplot(data = dataset_functions_split[[K]], aes(x=functions_log10, fill = Label)) +
      ggtitle(unique(dataset_functions_split[[K]]$Label)) + 
      geom_histogram(color="black", linewidth = .1, bins = 15) + 
      scale_x_log10(name = expression(paste("g." * m^-2 *"."* h^-1)), breaks = c(0.01, 10)) +
      scale_y_continuous(name = "", limits = c(0,400), breaks = c(0, 300), expand = c(0,0)) +
      scale_fill_manual(values = dataset_functions_split[[K]]$colors) + theme_classic() +
      theme(legend.position = "none", strip.text = element_text(hjust = 0),
            strip.background = element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 8),
            axis.title.x = element_text(size = 10),
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
            plot.title = element_text(size = 10, face = "bold"),
            axis.line = element_line(colour = "black")))
  
  K = 11
  (Fig_1B[[K]] <- ggplot(data = dataset_functions_split[[K]], aes(x=functions_log10, fill = Label)) +
      ggtitle(unique(dataset_functions_split[[K]]$Label)) + 
      geom_histogram(color="black", linewidth = .1, bins = 15) + 
      scale_x_log10(name = expression(paste("kg." * m^-2)), breaks = c(0.1, 100)) +
      scale_y_continuous(name = "", limits = c(0,400), breaks = c(0, 300), expand = c(0,0)) +
      scale_fill_manual(values = dataset_functions_split[[K]]$colors) + theme_classic() +
      theme(legend.position = "none", strip.text = element_text(hjust = 0),
            strip.background = element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 8),
            axis.title.x = element_text(size = 10),
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
            plot.title = element_text(size = 10, face = "bold"),
            axis.line = element_line(colour = "black")))
  
  K = 12
  (Fig_1B[[K]] <- ggplot(data = dataset_functions_split[[K]], aes(x=functions_log10, fill = Label)) +
      ggtitle(unique(dataset_functions_split[[K]]$Label)) + 
      geom_histogram(color="black", linewidth = .1, bins = 15) + 
      scale_x_log10(name = expression(paste("kg." * yr^-1)), breaks = c(0.1, 10)) +
      scale_y_continuous(name = "", limits = c(0,400), breaks = c(0, 300), expand = c(0,0)) +
      scale_fill_manual(values = dataset_functions_split[[K]]$colors) + theme_classic() +
      theme(legend.position = "none", strip.text = element_text(hjust = 0),
            strip.background = element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 8),
            axis.title.x = element_text(size = 10),
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
            plot.title = element_text(size = 10, face = "bold"),
            axis.line = element_line(colour = "black")))
  
  K = 13
  (Fig_1B[[K]] <- ggplot(data = dataset_functions_split[[K]], aes(x=functions_log10, fill = Label)) +
      ggtitle(unique(dataset_functions_split[[K]]$Label)) + 
      geom_histogram(color="black", linewidth = .1, bins = 15) + 
      scale_x_log10(name = "dimensionless", breaks = c(1, 5)) +
      scale_y_continuous(name = "", limits = c(0,400), breaks = c(0, 300), expand = c(0,0)) +
      scale_fill_manual(values = dataset_functions_split[[K]]$colors) + theme_classic() +
      theme(legend.position = "none", strip.text = element_text(hjust = 0),
            strip.background = element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 8),
            axis.title.x = element_text(size = 10),
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
            plot.title = element_text(size = 10, face = "bold"),
            axis.line = element_line(colour = "black")))
  
  K = 14
  (Fig_1B[[K]] <- ggplot(data = dataset_functions_split[[K]], aes(x=functions_log10, fill = Label)) +
      ggtitle(unique(dataset_functions_split[[K]]$Label)) + 
      geom_histogram(color="black", linewidth = .1, bins = 15) + 
      scale_x_log10(name = expression(cm^3), breaks = c(10, 10000)) +
      scale_y_continuous(name = "", limits = c(0,400), breaks = c(0, 300), expand = c(0,0)) +
      scale_fill_manual(values = dataset_functions_split[[K]]$colors) + theme_classic() +
      theme(legend.position = "none", strip.text = element_text(hjust = 0),
            strip.background = element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 8),
            axis.title.x = element_text(size = 10),
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
            plot.title = element_text(size = 10, face = "bold"),
            axis.line = element_line(colour = "black")))
  
  
  layout <- "
AAAAAAAAA
#########
#BCDEFGH#
#########
#ILMNOPQ#
"
  
  

  Fig_1 <- Figure_1A + 
    
  Fig_1B[[8]] + Fig_1B[[9]] + Fig_1B[[10]] + Fig_1B[[11]] + Fig_1B[[12]] + Fig_1B[[13]] + Fig_1B[[14]] +
  
  Fig_1B[[1]] + Fig_1B[[2]] + Fig_1B[[3]] + Fig_1B[[4]] + Fig_1B[[5]] + Fig_1B[[6]] + Fig_1B[[7]] +
  
  plot_layout(ncol = 1, nrow=5, heights = c(2.5,0.5, 1,0.2, 1), design = layout) 
  
  ggsave("output/Fig_1.png", width = 40, height = 20, units = "cm", dpi = 320)
  
}
  