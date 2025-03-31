

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
  
  return(list(dataset_functions, rls_functions))
  
}#eo combine_benthic_fish_functions




figure_1 <- function(list_dataset_functions) {
  
  dataset_functions <- list_dataset_functions[[1]]
  rls_functions <- list_dataset_functions[[2]]
  
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
             y = -4000000, 
             yend = -4000000, 
             x = -15850000, 
             xend = 15850000, 
             size = 0.5,
             color="grey80") +
    annotate(geom = 'segment', 
             y = 4000000, 
             yend = 4000000, 
             x = -15850000, 
             xend = 15850000, 
             size = 0.5,
             color="grey80") +
    
    theme(axis.line.x = element_line(linewidth = 1)) +
    
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
      geom_histogram(color="black", linewidth = .1, bins = 20) + 
      scale_x_log10(name = expression(paste("g."* h^-1)), breaks = c(0.1, 100)) +
      scale_y_continuous(name = "", limits = c(0,500), breaks = c(0, 2250,500), expand = c(0,0)) +
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
      geom_histogram(color="black", linewidth = .1, bins = 20) + 
      scale_x_log10(name = expression(paste("g."* h^-1)), breaks = c(0.1, 100)) +
      scale_y_continuous(name = "", limits = c(0,500), breaks = c(0, 250,500), expand = c(0,0)) +
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
      geom_histogram(color="black", linewidth = .1, bins = 20) + 
      scale_x_log10(name = expression(paste("kg."* yr^-1)), breaks = c(0.1, 100)) +
      scale_y_continuous(name = "", limits = c(0,500), breaks = c(0, 250,500), expand = c(0,0)) +
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
      geom_histogram(color="black", linewidth = .1, bins = 20) + 
      scale_x_log10(name = expression(paste("g."* h^-1)), breaks = c(0.1, 100)) +
      scale_y_continuous(name = "", limits = c(0,500), breaks = c(0, 250,500), expand = c(0,0)) +
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
      geom_histogram(color="black", linewidth = .1, bins = 20) + 
      scale_x_log10(name = expression(paste("g."* h^-1)), breaks = c(0.1, 100)) +
      scale_y_continuous(name = "", limits = c(0,500), breaks = c(0, 250,500), expand = c(0,0)) +
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
      geom_histogram(color="black", linewidth = .1, bins = 20) + 
      scale_x_log10(name = expression(paste("g."* h^-1)), breaks = c(0.1, 100)) +
      scale_y_continuous(name = "", limits = c(0,500), breaks = c(0, 250,500), expand = c(0,0)) +
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
      geom_histogram(color="black", linewidth = .1, bins = 20) + 
      scale_x_log10(name = expression(paste("g."* h^-1)), breaks = c(0, 0.001)) +
      scale_y_continuous(name = "", limits = c(0,500), breaks = c(0, 250,500), expand = c(0,0)) +
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
      geom_histogram(color="black", linewidth = .1, bins = 20) + 
      scale_x_log10(name = expression(paste("g."* h^-1)), breaks = c(0.01, 10)) +
      scale_y_continuous(name = "", limits = c(0,500), breaks = c(0, 250,500), expand = c(0,0)) +
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
      geom_histogram(color="black", linewidth = .1, bins = 20) + 
      scale_x_log10(name = expression(paste("kg." * yr^-1)), breaks = c(10, 10000)) +
      scale_y_continuous(name = "", limits = c(0,500), breaks = c(0, 250,500), expand = c(0,0)) +
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
      geom_histogram(color="black", linewidth = .1, bins = 20) + 
      scale_x_log10(name = expression(paste("g." * m^-2 *"."* h^-1)), breaks = c(0.01, 10)) +
      scale_y_continuous(name = "", limits = c(0,500), breaks = c(0, 250,500), expand = c(0,0)) +
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
      geom_histogram(color="black", linewidth = .1, bins = 20) + 
      scale_x_log10(name = expression(paste("kg." * m^-2)), breaks = c(0.1, 100)) +
      scale_y_continuous(name = "", limits = c(0,500), breaks = c(0, 250, 500), expand = c(0,0)) +
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
      geom_histogram(color="black", linewidth = .1, bins = 20) + 
      scale_x_log10(name = expression(paste("kg." * yr^-1)), breaks = c(0.1, 10)) +
      scale_y_continuous(name = "", limits = c(0,500), breaks = c(0, 250,500), expand = c(0,0)) +
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
      geom_histogram(color="black", linewidth = .1, bins = 20) + 
      scale_x_log10(name = "dimensionless", breaks = c(1, 5)) +
      scale_y_continuous(name = "", limits = c(0,500), breaks = c(0, 250,500), expand = c(0,0)) +
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
      geom_histogram(color="black", linewidth = .1, bins = 20) + 
      scale_x_log10(name = expression(cm^3), breaks = c(1, 10000)) +
      scale_y_continuous(name = "", limits = c(0,500), breaks = c(0, 250,500), expand = c(0,0)) +
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
  
}#eo fig 1



pca_functions <- function(list_dataset_functions, benth, fish) {

  dataset_functions <- list_dataset_functions[[1]]
  
  f.fun <- c("herb", "plank", "prod","pisc",  "exN", "exP", "turn")
  b.fun <- c("GPP_ghr", "OrgGrowth_gyr",  "Calc_ghr", "Storage_kg", "Accretion_kgyr",  "Rugosity", "BranchSpace_cm3") 
  
  f.cols <- fish(7, option = fish_palettes()[127])
  b.cols <- fish(7, option = "Holocentrus_adscensionis")
  
  
  labels <- rbind(data.frame(fun = f.fun, type="fish", lab=NA, cols=f.cols), data.frame(fun = b.fun, type="benthos", lab=NA, cols=b.cols))
  labels$lab[labels$fun=="herb"] <- "HERB"
  labels$lab[labels$fun=="plank"] <- "PLANK"
  labels$lab[labels$fun=="prod"] <- "PROD"
  labels$lab[labels$fun=="pisc"] <- "PISC"
  labels$lab[labels$fun=="exN"] <- "EX.N"
  labels$lab[labels$fun=="exP"] <- "EX.P"
  labels$lab[labels$fun=="turn"] <- "TURN"
  labels$lab[labels$fun=="GPP_ghr"] <- "GPP"
  labels$lab[labels$fun=="Rugosity"] <- "RUG"
  labels$lab[labels$fun=="Calc_ghr"] <- "CALC"
  labels$lab[labels$fun=="Storage_kg"] <- "STOR"
  labels$lab[labels$fun=="OrgGrowth_gyr"] <- "ORG.G"
  labels$lab[labels$fun=="BranchSpace_cm3"] <- "SPACE"
  labels$lab[labels$fun=="Accretion_kgyr"] <- "INORG.G"
  labels
  
  labels$n <- c(1:nrow(labels))
  
  length(unique(fish$SurveyID))
  length(unique(benth$SurveyID))
  
  sites <- data.frame(SurveyID = unique(c(benth$SurveyID, fish$SurveyID)))
  head(sites)
  length(unique(sites$SurveyID))
  
  sites$benthos <- benth$SurveyID[match(sites$SurveyID, benth$SurveyID)]
  sites$fish <- fish$SurveyID[match(sites$SurveyID, fish$SurveyID)]
  head(sites)
  
  all <- sites[complete.cases(sites),]
  nrow(all)
  
  all[, f.fun] <- fish[match(all$SurveyID, fish$SurveyID), f.fun]
  all[, b.fun] <- benth[match(all$SurveyID, fish$SurveyID), b.fun]
  
  all <- all[complete.cases(all),]
  nrow(all)
  
  rownames(all)<-all$SurveyID
  
  log.all <- all[,c(f.fun, b.fun)]+1
  #case_study <- case_study[,c(f.fun, b.fun)]+1
  #log.all <- rbind(log.all, case_study)
  #log.all[log.all==0] <- NA
  #log.all <- log.all[complete.cases(log.all),]
  log.all <- log(log.all)
  #nrow(log.all)
  
  pca_fun <- prcomp(log.all, center=T, scale=T)
  # biplot(pca)
  
  save(pca_fun, file="output/pca_functions.RData")
  
  return(list(pca_fun, all, labels))
  

  
}#eo pca_functions
  
  

figure_2 <- function(list_pca_functions, benth, fish) {
  
pca <- list_pca_functions[[1]]

all <- list_pca_functions[[2]]

labels <- list_pca_functions[[3]]


  
pc.df <- data.frame(pca$x[,1:6])
vecs <- data.frame(pca$rotation[,1:6])
vecs$labs <- labels$lab[match(rownames(vecs), labels$fun)]
vecs$type <- labels$type[match(rownames(vecs), labels$fun)]
vecs$n <- labels$n[match(rownames(vecs), labels$fun)]
vars <- round((pca$sdev^2 / sum(pca$sdev^2)), 3)*100

head(pc.df)

mult <- min((max(pc.df[,"PC2"]) - min(pc.df[,"PC2"])/(max(vecs[,"PC2"])-min(vecs[,"PC2"]))),(max(pc.df[,"PC1"]) - min(pc.df[,"PC1"])/(max(vecs[,"PC1"])-min(vecs[,"PC1"]))))
vecs$PC1b <- vecs$PC1 * mult * 0.7
vecs$PC2b <- vecs$PC2 * mult * 0.7
vecs$PC3b <- vecs$PC3# * mult * 0.7

cols <- labels$cols
names(cols) <- labels$lab

pc.df$SurveyID <- rownames(pc.df)
pc.df$biom <- fish$biom[match(pc.df$SurveyID, fish$SurveyID)]
pc.df$cover <- benth$coral[match(pc.df$SurveyID, benth$SurveyID)]

options(ggrepel.max.overlaps = Inf)

p1 <- ggplot()+
  geom_point(data=pc.df, aes(PC1, PC2, fill=log(biom)), shape=21, col="white", stroke=0.2, alpha=0.5)+
  geom_segment(data=vecs, aes(x=0, y=0, xend=PC1b, yend=PC2b, col=labs), arrow=arrow(length=unit(1,"mm")))+
  geom_label_repel(data=vecs, aes(PC1b, PC2b, label=n, colour=labs),fontface="bold",size=2.5, force=0.0005, label.padding=unit(0.2, "mm"), label.size=1)+
  geom_label_repel(data=vecs, aes(PC1b, PC2b, label=n),fontface="bold",size=2.5, force=0.0005, label.padding=unit(0.2, "mm"), label.size=0)+
  #geom_text_repel(data=vecs, aes(PC1b, PC2b, label=n),fontface="bold",size=3, force=0.0005)+
  guides(colour="none", fill="none")+
  scale_fill_distiller(palette="Greys", direction=1)+
  scale_colour_manual(values=cols)+
  theme_bw()+theme(panel.grid=element_blank(), legend.title=element_blank(), legend.key.width=unit(1, "mm"), legend.key.height=unit(1, "mm"), legend.position=c(0.16, 0.93), legend.background=element_blank())+
  labs(x=paste("PC1 (", vars[1], "%)", sep=""), y=paste("PC2 (", vars[2], "%)", sep=""))


p2 <- ggplot()+
  geom_point(data=pc.df, aes(PC2, PC3, fill=log(biom)), shape=21, col="white", stroke=0.2, alpha=0.5)+
  geom_segment(data=vecs, aes(x=0, y=0, xend=PC2*7, yend=PC3*4, col=labs), arrow=arrow(length=unit(1,"mm")))+
  geom_label_repel(data=vecs, aes(PC2*7, PC3*4, label=n, colour=labs),fontface="bold",size=2.5, force=0.0005, label.padding=unit(0.2, "mm"), label.size=1)+
  geom_label_repel(data=vecs, aes(PC2*7, PC3*4, label=n),fontface="bold",size=2.5, force=0.0005, label.padding=unit(0.2, "mm"), label.size=0)+
  guides(colour="none", fill="none")+
  scale_colour_manual(values=cols)+
  scale_fill_distiller(palette="Greys", direction=1)+
  theme_bw()+theme(panel.grid=element_blank())+
  labs(x=paste("PC2 (", vars[2], "%)", sep=""), y=paste("PC3 (", vars[3], "%)", sep=""))


p3 <- ggplot()+
  geom_point(data=pc.df, aes(PC1, PC4, fill=log(biom)), shape=21, col="white", stroke=0.2, alpha=0.5)+
  geom_segment(data=vecs, aes(x=0, y=0, xend=PC1*7, yend=PC4*4, col=labs), arrow=arrow(length=unit(1,"mm")))+
  #geom_text_repel(data=vecs, aes(PC1*8, PC4*5, label=n), col='black', fontface="bold",size=3, 
  geom_label_repel(data=vecs, aes(PC1*7, PC4*4, label=n, colour=labs), fontface="bold",size=2.5, force=0.0008, label.padding=unit(0.2, "mm"), label.size=1)+
  geom_label_repel(data=vecs, aes(PC1*7, PC4*4, label=n),fontface="bold",size=2.5, force=0.0008, label.padding=unit(0.2, "mm"), label.size=0)+
  guides(colour="none", fill="none")+
  scale_colour_manual(values=cols)+
  scale_fill_distiller(palette="Greys", direction=1)+
  theme_bw()+theme(panel.grid=element_blank())+
  labs(x=paste("PC1 (", vars[1], "%)", sep=""), y=paste("PC4 (", vars[4], "%)", sep=""))


vecs$lab2 <- paste("(", vecs$n, ") ", vecs$labs, sep="")
vecs$lab2 <- factor(vecs$lab2, levels=rev(vecs$lab2))

vecs2 <- melt(vecs[,c("PC1", "PC2", "PC3", "PC4", "labs", "n", "lab2")], id.var=c("labs", "n", "lab2"))
head(vecs2)

p4 <- ggplot(vecs2, aes(y=lab2, x=value))+
  geom_bar(stat="identity", aes(x=abs(value)), fill="white")+
  geom_bar(stat="identity", aes(x=-value), fill="white")+
  geom_bar(stat="identity", aes(fill=labs), col="black", size=0.1)+
  geom_vline(xintercept=0)+
  guides(fill="none")+
  scale_fill_manual(values=cols)+
  facet_wrap(~variable, nrow=1, scales="free_x")+
  xlab("Loading")+
  theme_classic()+theme(axis.title.y=element_blank(), strip.background=element_blank(), axis.text.x=element_text(angle=45, hjust=1))



all$biom <- fish$biom[match(all$SurveyID, fish$SurveyID)]
all$cover <- benth$coral[match(all$SurveyID, benth$SurveyID)]


p5 <- ggplot(all, aes(cover, biom/1000))+
  geom_point(col="grey", shape=21, aes(fill=log(biom)), alpha=0.5)+
  guides(fill="none")+
  geom_smooth(col="black", se=F, linetype="dashed", linewidth=0.5)+
  scale_y_log10()+
  #scale_x_sqrt()+
  scale_fill_distiller(palette="Greys", direction=1)+
  theme_classic()+
  labs(x="Coral cover (%)", y="Fish biomass (kg)")


Fig_2 <- plot_grid(
  plot_grid(p1, p2, p3, nrow=1, labels=c("a", "b", "c")), 
  plot_grid(p4, p5, rel_widths=c(1, 0.8), labels=c("d","e")), 
  ncol=1)


ggsave("output/Fig_2.png", width = 28, height = 17, units = "cm", dpi = 320)

  
}
  

##### --------> Prepare data for PC models

prepare_PC_data <- function(list_pca_functions, pred) {
  
  dat <- as.data.frame(list_pca_functions[[1]]$x)
  dat$SurveyID <- rownames(dat)
  
  ################################# 
  
  IDs <- dat$SurveyID
  
  env <- c("SurveyID", "mean_chl_5year", "mean_sss_5year", "mean_npp_5year", "mean_sst_5year", "mean_pH_1year_5year")
  
  environment <- lapply(IDs, function(x) {
    
    na.omit(pred[pred$SurveyID == x, which(colnames(pred) %in% env)])[1,]  
  })
  
  env_df <- do.call(rbind, environment) 
  env_df <- na.omit(env_df)
  
  dat <- dat[dat$SurveyID %in% env_df$SurveyID,]
  pred <- pred[pred$SurveyID %in% env_df$SurveyID,]
  
  pred$No.take.multizoned[pred$No.take.multizoned == "No take"] <- "No take"
  pred$No.take.multizoned[pred$No.take.multizoned == "No take "] <- "No take"
  pred$No.take.multizoned[pred$No.take.multizoned == "No Take"] <- "No take"
  pred$No.take.multizoned[pred$No.take.multizoned == "No take multizoned"] <- "No take"
  pred$No.take.multizoned[pred$No.take.multizoned == "Restricted take multizoned"] <- "Restricted take"
  pred$No.take.multizoned[pred$No.take.multizoned == "Restricted take"] <- "Restricted take"
  pred$No.take.multizoned[pred$No.take.multizoned == "Fishing"] <- "Fishing"
  pred$No.take.multizoned[is.na(pred$No.take.multizoned)] <- "Not protected"
  
  IDs <- dat$SurveyID
  
  DHW <- lapply(IDs, function(x) {
    
    c(na.omit(pred[pred$SurveyID == x,]$mean_DHW_5year)[1], na.omit(pred[pred$SurveyID == x,]$mean_DHW_1year[1]))  
  })
  
  names(DHW) <- IDs
  
  
  DHW_df <- do.call(rbind, DHW) 
  
  
  sites <- sapply(IDs, function(x) {
    
    pred[pred$SurveyID == x,]$SiteCode[1]
    
  })
  
  gravtot <- sapply(IDs, function(x) {
    
    pred[pred$SurveyID == x,]$gravtot2[1]
    
  })
  
  mpa <- sapply(IDs, function(x) {
    
    pred[pred$SurveyID == x,]$No.take.multizoned[1]
    
  })
  
  pc <- prcomp(env_df[,2:6],
               center = TRUE,
               scale = TRUE)
  
  dat$DHW5 <- DHW_df[,1]
  dat$DHW1 <- DHW_df[,2]
  dat$gravtot <- gravtot
  dat$SiteCode <- sites
  dat$mpa <- mpa
  dat$env1 <- pc$x[,1]
  dat$env2 <- pc$x[,2]
  
  dat
  
}


######## prepare data for benthic models

prepare_benthic_models <- function(pred, benthos) {
  
  ################################# BENTHOS
  
  IDsB <- benthos$SurveyID
  
  env <- c("SurveyID", "mean_chl_5year", "mean_sss_5year", "mean_npp_5year", "mean_sst_5year", "mean_pH_1year_5year")
  environment <- lapply(IDsB, function(x) {
    
    na.omit(pred[pred$SurveyID == x, which(colnames(pred) %in% env)])[1,]  
  })
  
  env_df <- do.call(rbind, environment) 
  env_df <- na.omit(env_df)
  
  benthos <- benthos[benthos$SurveyID %in% env_df$SurveyID,]
  #fish <- fish[fish$SurveyID %in% env_df$SurveyID,]
  pred <- pred[pred$SurveyID %in% env_df$SurveyID,]
  
  pred$No.take.multizoned[pred$No.take.multizoned == "No take"] <- "No take"
  pred$No.take.multizoned[pred$No.take.multizoned == "No take "] <- "No take"
  pred$No.take.multizoned[pred$No.take.multizoned == "No Take"] <- "No take"
  pred$No.take.multizoned[pred$No.take.multizoned == "No take multizoned"] <- "No take"
  pred$No.take.multizoned[pred$No.take.multizoned == "Restricted take multizoned"] <- "Restricted take"
  pred$No.take.multizoned[pred$No.take.multizoned == "Restricted take"] <- "Restricted take"
  pred$No.take.multizoned[pred$No.take.multizoned == "Fishing"] <- "Fishing"
  pred$No.take.multizoned[is.na(pred$No.take.multizoned)] <- "Not protected"
  
  IDsB <- benthos$SurveyID
  
  
  DHWB <- lapply(IDsB, function(x) {
    
    c(na.omit(pred[pred$SurveyID == x,]$mean_DHW_5year)[1], na.omit(pred[pred$SurveyID == x,]$mean_DHW_1year[1]))  
    
    
  })
  
  names(DHWB) <- IDsB
  
  
  DHWB_df <- do.call(rbind, DHWB) 
  
  
  sites <- sapply(IDsB, function(x) {
    
    pred[pred$SurveyID == x,]$SiteCode[1]
    
  })
  
  gravtot <- sapply(IDsB, function(x) {
    
    pred[pred$SurveyID == x,]$gravtot2[1]
    
  })
  
  mpa <- sapply(IDsB, function(x) {
    
    pred[pred$SurveyID == x,]$No.take.multizoned[1]
    
  })
  
  pc <- prcomp(env_df[,2:6],
               center = TRUE,
               scale. = TRUE)
  
  benthos$DHW5 <- DHWB_df[,1]
  benthos$DHW1 <- DHWB_df[,2]
  benthos$gravtot <- gravtot
  benthos$SiteCode <- sites
  benthos$mpa <- mpa
  benthos$env1 <- pc$x[,1]
  benthos$env2 <- pc$x[,2]
  
  benthos
  
}





##################################### Prepare data for fish models


prepare_fish_models <- function(pred, fish) {
  IDsF <- fish$SurveyID
  
  env <- c("SurveyID", "mean_chl_5year", "mean_sss_5year", "mean_npp_5year", "mean_sst_5year", "mean_pH_1year_5year")
  environment <- lapply(IDsF, function(x) {
    
    na.omit(pred[pred$SurveyID == x, which(colnames(pred) %in% env)])[1,]  
  })
  
  env_df <- do.call(rbind, environment) 
  env_df <- na.omit(env_df)
  
  fish <- fish[fish$SurveyID %in% env_df$SurveyID,]
  pred <- pred[pred$SurveyID %in% env_df$SurveyID,]
  
  
  
  pred$No.take.multizoned[pred$No.take.multizoned == "No take"] <- "No take"
  pred$No.take.multizoned[pred$No.take.multizoned == "No take "] <- "No take"
  pred$No.take.multizoned[pred$No.take.multizoned == "No Take"] <- "No take"
  pred$No.take.multizoned[pred$No.take.multizoned == "No take multizoned"] <- "No take"
  pred$No.take.multizoned[pred$No.take.multizoned == "Restricted take multizoned"] <- "Restricted take"
  pred$No.take.multizoned[pred$No.take.multizoned == "Restricted take"] <- "Restricted take"
  pred$No.take.multizoned[pred$No.take.multizoned == "Fishing"] <- "Fishing"
  pred$No.take.multizoned[is.na(pred$No.take.multizoned)] <- "Not protected"
  
  IDsF <- fish$SurveyID
  
  DHWF <- lapply(IDsF, function(x) {
    
    c(na.omit(pred[pred$SurveyID == x,]$mean_DHW_5year)[1], na.omit(pred[pred$SurveyID == x,]$mean_DHW_1year[1]))  
    
    
  })
  
  names(DHWF) <- IDsF
  
  DHWF_df <- do.call(rbind, DHWF) 
  
  
  sites <- sapply(IDsF, function(x) {
    
    pred[pred$SurveyID == x,]$SiteCode[1]
    
  })
  
  gravtot <- sapply(IDsF, function(x) {
    
    pred[pred$SurveyID == x,]$gravtot2[1]
    
  })
  
  mpa <- sapply(IDsF, function(x) {
    
    pred[pred$SurveyID == x,]$No.take.multizoned[1]
    
  })
  
  environment <- lapply(IDsF, function(x) {
    
    na.omit(pred[pred$SurveyID == x, which(colnames(pred) %in% env)])[1,]  
  })
  
  env_df <- do.call(rbind, environment) 
  env_df <- na.omit(env_df)
  
  pc <- prcomp(env_df[,2:6],
               center = TRUE,
               scale = TRUE)
  
  fish$DHW5 <- DHWF_df[,1]
  fish$DHW1 <- DHWF_df[,2]
  fish$SiteCode <- sites
  fish$gravtot <- gravtot
  fish$mpa <- mpa
  fish$env1 <- pc$x[,1]
  fish$env2 <- pc$x[,2]
  
  fish
}



###########  PC models env 

run_models_env <- function(dat_mod) {
  
  library(brms)
  
  pr <- get_prior(PC1 ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = dat_mod)
  pr$prior[which(pr$class == "b")] <- "normal(0,1)"
  
  fit_pc1_env <- brm(PC1 ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = dat_mod, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000, seed = 10, prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  
  pr <- get_prior(PC2 ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = dat_mod)
  pr$prior[which(pr$class == "b")] <- "normal(0,1)"
  
  fit_pc2_env <- brm(PC2 ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = dat_mod, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000, seed = 10, prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  
  pr <- get_prior(PC3 ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = dat_mod)
  pr$prior[which(pr$class == "b")] <- "normal(0,1)"
  
  fit_pc3_env <- brm(PC3 ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = dat_mod, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000, seed = 10, prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  
  pr <- get_prior(PC4 ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = dat_mod)
  pr$prior[which(pr$class == "b")] <- "normal(0,1)"
  
  fit_pc4_env <- brm(PC4 ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = dat_mod, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000, seed = 10, prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  
  fit_pc_env <- list(fit_pc1_env,fit_pc2_env, fit_pc3_env, fit_pc4_env)
  
  fit_pc_env
}


###########  PC models no_env 


run_models_noenv <- function(dat_mod) {
  
  library(brms)
  
  pr <- get_prior(PC1 ~ DHW1  + gravtot + mpa + (1|SiteCode), data = dat_mod)
  pr$prior[which(pr$class == "b")] <- "normal(0,1)"
  
  fit_pc1_noenv <- brm(PC1 ~ DHW1  + gravtot + mpa  + (1|SiteCode), data = dat_mod, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000, seed = 10, prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  
  pr <- get_prior(PC2 ~ DHW1  + gravtot + mpa + (1|SiteCode), data = dat_mod)
  pr$prior[which(pr$class == "b")] <- "normal(0,1)"
  
  fit_pc2_noenv <- brm(PC2 ~ DHW1  + gravtot + mpa  + (1|SiteCode), data = dat_mod, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000, seed = 10, prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  
  pr <- get_prior(PC3 ~ DHW1  + gravtot + mpa + (1|SiteCode), data = dat_mod)
  pr$prior[which(pr$class == "b")] <- "normal(0,1)"
  
  fit_pc3_noenv <- brm(PC3 ~ DHW1  + gravtot + mpa  + (1|SiteCode), data = dat_mod, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000, seed = 10, prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  
  pr <- get_prior(PC4 ~ DHW1  + gravtot + mpa + (1|SiteCode), data = dat_mod)
  pr$prior[which(pr$class == "b")] <- "normal(0,1)"
  
  fit_pc4_noenv <- brm(PC4 ~ DHW1  + gravtot + mpa  + (1|SiteCode), data = dat_mod, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000, seed = 10, prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  
  fit_pc_noenv <- list(fit_pc1_noenv,fit_pc2_noenv, fit_pc3_noenv, fit_pc4_noenv)
  
  fit_pc_noenv
  
}





ind_benthic_models_env <- function(benthos) {
  
  ###########  IND benthic models env
  
  pr <- get_prior(log1p(Calc_ghr) ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = benthos)
  pr$prior[which(pr$class == "b")] <- "normal(0,1)"
  
  #fit env
  
  #control = list(adapt_delta = 0.9, max_treedepth = 12),
  
  fit_cc <- brm(log1p(coral) ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = benthos, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000, seed = 10,  prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  fit_cal <- brm(log1p(Calc_ghr) ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = benthos, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000, seed = 10,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  fit_sto <- brm(log1p(Storage_kg) ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = benthos, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000, seed = 10,  prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  fit_rug <- brm(log1p(Rugosity) ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = benthos, backend = "cmdstanr", cores = 4, threads = 30, iter = 10000, seed = 10,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 15))
  fit_gpp <- brm(log1p(GPP_ghr) ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = benthos, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000, seed = 10,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  fit_brs <- brm(log1p(BranchSpace_cm3) ~ DHW1  + mpa + gravtot + env1 + env2 + (1|SiteCode), data = benthos, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000, prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  fit_ing <- brm(log1p(Accretion_kgyr) ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = benthos, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000,  prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  fit_org <- brm(log1p(OrgGrowth_gyr) ~ DHW1  + gravtot + mpa +env1 + env2 +  (1|SiteCode), data = benthos, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  
  
  fit_b_env <- list(fit_cc,fit_cal, fit_sto, fit_rug, fit_gpp, fit_brs, fit_org, fit_ing)
  
  fit_b_env
}



#fit no env

ind_benthic_models_noenv <- function(benthos) {
  
  pr <- get_prior(log1p(Calc_ghr) ~ DHW1  + gravtot + mpa + (1|SiteCode), data = benthos)
  pr$prior[which(pr$class == "b")] <- "normal(0,1)"
  
  fit_cc_noenv <- brm(log1p(coral) ~ DHW1  + gravtot + mpa + (1|SiteCode), data = benthos, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000, seed = 10,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  fit_cal_noenv <- brm(log1p(Calc_ghr) ~ DHW1  + gravtot + mpa + (1|SiteCode), data = benthos, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000, seed = 10,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  fit_sto_noenv <- brm(log1p(Storage_kg) ~ DHW1  + gravtot + mpa  + (1|SiteCode), data = benthos, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000, seed = 10,  prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  fit_rug_noenv <- brm(log1p(Rugosity) ~ DHW1  + gravtot + mpa  + (1|SiteCode), data = benthos, backend = "cmdstanr", cores = 4, threads = 30, iter = 10000, seed = 10,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 15))
  fit_gpp_noenv <- brm(log1p(GPP_ghr) ~ DHW1  + gravtot + mpa + (1|SiteCode), data = benthos, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000, seed = 10,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  fit_brs_noenv <- brm(log1p(BranchSpace_cm3) ~ DHW1  + mpa + gravtot  + (1|SiteCode), data = benthos, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  fit_ing_noenv <- brm(log1p(Accretion_kgyr) ~ DHW1  + gravtot + mpa + (1|SiteCode), data = benthos, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  fit_org_noenv <- brm(log1p(OrgGrowth_gyr) ~ DHW1  + gravtot + mpa +  (1|SiteCode), data = benthos, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  
  
  fit_b_noenv <- list(fit_cc_noenv,fit_cal_noenv, fit_sto_noenv, fit_rug_noenv, fit_gpp_noenv, fit_brs_noenv, fit_org_noenv, fit_ing_noenv)
  
  fit_b_noenv 
}



###########  fish models 

ind_fish_models_env <- function(fish) {
  
  library(brms)
  
  pr <- get_prior(log1p(herb) ~ DHW1  + gravtot + mpa + env1 + env2 +  (1|SiteCode), data = fish)
  pr$prior[1] <- "normal(0,1)"
  
  
  #fit env
  
  fit_her <- brm(log1p(herb) ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = fish, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000, seed = 10,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  fit_pla <- brm(log1p(plank) ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = fish, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000, seed = 10,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  fit_pis <- brm(log1p(pisc) ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = fish, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000, seed = 10,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  fit_pro <- brm(log1p(prod) ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = fish, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000, seed = 10,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  fit_exP <- brm(log1p(exP) ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = fish, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000, seed = 10,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  fit_exN <- brm(log1p(exN) ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = fish, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000, seed = 10,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  fit_turn <- brm(log1p(turn) ~ DHW1  + gravtot + mpa + env1 + env2 + (1|SiteCode), data = fish, backend = "cmdstanr", cores = 4, threads = 30, iter = 10000, seed = 10,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 15))
  
  
  fit_f_env <- list(fit_her, fit_pla, fit_pis, fit_pro, fit_exP, fit_exN, fit_turn)
  
  fit_f_env 
  
}


#fit no env

ind_fish_models_noenv <- function(fish) {
  
  pr <- get_prior(log1p(herb) ~ DHW1  + gravtot + mpa + (1|SiteCode), data = fish)
  pr$prior[1] <- "normal(0,1)"
  
  fit_her_noenv <- brm(log1p(herb) ~ DHW1  + gravtot + mpa  + (1|SiteCode), data = fish, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000, seed = 10,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  fit_pla_noenv <- brm(log1p(plank) ~ DHW1  + gravtot + mpa + (1|SiteCode), data = fish, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000, seed = 10,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  fit_pis_noenv <- brm(log1p(pisc) ~ DHW1  + gravtot + mpa  + (1|SiteCode), data = fish, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000, seed = 10,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  fit_pro_noenv <- brm(log1p(prod) ~ DHW1  + gravtot + mpa  + (1|SiteCode), data = fish, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000, seed = 10,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  fit_exP_noenv <- brm(log1p(exP) ~ DHW1  + gravtot + mpa  + (1|SiteCode), data = fish, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000, seed = 10,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  fit_exN_noenv <- brm(log1p(exN) ~ DHW1  + gravtot + mpa  + (1|SiteCode), data = fish, backend = "cmdstanr", cores = 4, threads = 30, iter = 5000, seed = 10,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 12))
  fit_turn_noenv <- brm(log1p(turn) ~ DHW1  + gravtot + mpa  + (1|SiteCode), data = fish, backend = "cmdstanr", cores = 4, threads = 30, iter = 10000, seed = 10,   prior = pr, control = list(adapt_delta = 0.9, max_treedepth = 15))
  
  
  fit_f_noenv <- list(fit_her_noenv, fit_pla_noenv, fit_pis_noenv, fit_pro_noenv, fit_exP_noenv, fit_exN_noenv, fit_turn_noenv)
  
  fit_f_noenv 
}








plot_pc <- function(labels2, models_env, models_noenv) {
  
  #load("output/fit_pcs.RData")
  #load("output/fit_pcs_noenv.RData")
  
  #fit_pcs <- fish_models_env
  #fit_pcs_noenv <- fish_models_noenv
  
  stat_interval <- function(..., .width = .95) {
    stat_summary(..., fun.data = ggdist::median_hdci,
                 geom = "linerange", fun.args = list(.width = .width))
  }
  stat_median <- function(...) {
    stat_summary(..., fun = median, geom = "point")
  }
  extract_named_vec <- function(df, col, names_from = "fun") {
    out <- df[[col]]
    names(out) <- df[[names_from]]
    out
  }
  ref <- labels2
  model_list <- list(
    `(1) PC1` = list(models_env[[1]], models_noenv[[1]]),
    `(2) PC2` = list(models_env[[2]], models_noenv[[2]]),
    `(3) PC3` = list(models_env[[3]], models_noenv[[3]]),
    `(4) PC4` = list(models_env[[4]], models_noenv[[4]])
  )
  sample_post <- purrr::map_dfr(model_list, function(x) {
    noenv <- x[[2]]
    x <- x[[1]]
    r2 <- brms::bayes_R2(x, summary = FALSE, re_formula = NA)
    r2_noenv <- brms::bayes_R2(noenv, summary = FALSE, re_formula = NA)
    brms::as_draws_df(x, variable = c("b_DHW1", "b_gravtot", "b_mpaNotake")) |>
      dplyr::mutate(
        `italic(R^2) * \" w/o env.\"` = .env$r2_noenv[, 1],
        `italic(R^2) * \" w/ env.\"` = .env$r2[, 1],
        b_gravtot = b_gravtot * 1e4
      ) |>
      dplyr::rename(
        DHW = b_DHW1, `\"Gravity (x 10\"^-4 * \")\"` = b_gravtot,
        `No-take` = b_mpaNotake
      ) |>
      dplyr::select(
        DHW, `\"Gravity (x 10\"^-4 * \")\"`, `No-take`,
        `italic(R^2) * \" w/ env.\"`,`italic(R^2) * \" w/o env.\"`
      ) |>
      tidyr::pivot_longer(
        cols = everything(), names_to = "predictor", values_to = "value"
      ) |>
      suppressWarnings()
  }, .id = "lab") |>
    dplyr::left_join(ref, by = join_by(lab)) |>
    dplyr::mutate(
      lab = factor(
        lab, levels = c("(1) PC1", "(2) PC2", "(3) PC3", "(4) PC4")
      ),
      predictor = factor(
        predictor, levels = c(
          "DHW", "\"Gravity (x 10\"^-4 * \")\"", "No-take",
          "italic(R^2) * \" w/o env.\"", "italic(R^2) * \" w/ env.\""
        )
      )
    )
  cols <- extract_named_vec(ref, "cols")
  labs <- extract_named_vec(ref, "lab")
  fig3c <- ggplot(data = sample_post) +
    aes(x = lab, y = value, colour = lab) +
    geom_hline(
      data = dplyr::filter(sample_post, !grepl("env\\.", predictor)),
      mapping = aes(yintercept = 0), linetype = 2
    ) +
    stat_interval(.width = .80, linewidth = 1.6, colour = "black") +
    stat_interval(.width = .95, linewidth = 0.6, colour = "black") +
    stat_interval(.width = .80, linewidth = 1.5) +
    stat_interval(.width = .95, linewidth = 0.5) +
    stat_median(size = 1.8, colour = "black") +
    stat_median(size = 1.7) +
    scale_colour_manual(values = cols) +
    coord_flip() +
    labs(x = NULL, y = NULL) + 
    labs(caption = "Medians with 95% (thin), 80% (thick) C.I.") +
    facet_grid(~ predictor, scales = "free_x", labeller = label_parsed) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      panel.spacing = unit(1.5, "cm", data = NULL),
      panel.background = element_rect(color = 'grey50')
    )
  #ggsave("output/fig3_c_pcs.png", fig3c, dev = "png", width = 10, height = 3.2)
  
  fig3c
}



################ FIG 3


figure_3 <- function(list_pca_functions, pred) {

  
pca <- list_pca_functions[[1]]  

f.fun <- c("herb", "plank", "prod","pisc",  "exN", "exP", "turn")
b.fun <- c("GPP_ghr", "OrgGrowth_gyr",  "Calc_ghr", "Storage_kg", "Accretion_kgyr",  "Rugosity", "BranchSpace_cm3") 

f.cols <- fish(7, option = fish_palettes()[127])
b.cols <- fish(7, option = "Holocentrus_adscensionis")

labels <- rbind(data.frame(fun = f.fun, type="fish", lab=NA, cols=f.cols), data.frame(fun = b.fun, type="benthos", lab=NA, cols=b.cols))
labels$lab[labels$fun=="herb"] <- "HERB"
labels$lab[labels$fun=="plank"] <- "PLANK"
labels$lab[labels$fun=="prod"] <- "PROD"
labels$lab[labels$fun=="pisc"] <- "PISC"
labels$lab[labels$fun=="exN"] <- "EX.N"
labels$lab[labels$fun=="exP"] <- "EX.P"
labels$lab[labels$fun=="turn"] <- "TURN"
labels$lab[labels$fun=="GPP_ghr"] <- "GPP"
labels$lab[labels$fun=="Rugosity"] <- "RUG"
labels$lab[labels$fun=="Calc_ghr"] <- "CALC"
labels$lab[labels$fun=="Storage_kg"] <- "STOR"
labels$lab[labels$fun=="OrgGrowth_gyr"] <- "ORG.G"
labels$lab[labels$fun=="BranchSpace_cm3"] <- "SPACE"
labels$lab[labels$fun=="Accretion_kgyr"] <- "INORG.G"

labels$n <- c(1:nrow(labels))



# biplot(pca)

pc.df <- data.frame(pca$x[,1:6])
vecs <- data.frame(pca$rotation[,1:6])
vecs$labs <- labels$lab[match(rownames(vecs), labels$fun)]
vecs$type <- labels$type[match(rownames(vecs), labels$fun)]
vecs$n <- labels$n[match(rownames(vecs), labels$fun)]
vars <- round((pca$sdev^2 / sum(pca$sdev^2)), 3)*100


mult <- min((max(pc.df[,"PC2"]) - min(pc.df[,"PC2"])/(max(vecs[,"PC2"])-min(vecs[,"PC2"]))),(max(pc.df[,"PC1"]) - min(pc.df[,"PC1"])/(max(vecs[,"PC1"])-min(vecs[,"PC1"]))))
vecs$PC1b <- vecs$PC1 * mult * 0.7
vecs$PC2b <- vecs$PC2 * mult * 0.7
vecs$PC3b <- vecs$PC3# * mult * 0.7

cols <- labels$cols
names(cols) <- labels$lab

# 2) add human predictors

pc.df$dhw <- pred$mean_DHW_1year[match(rownames(pc.df), pred$SurveyID)]
pc.df$gravity <- pred$gravtot2[match(rownames(pc.df), pred$SurveyID)]

pc.df <- na.omit(pc.df) 

#plot_grid(
#  ggplot(pc.df, aes(PC2, PC3, fill=sqrt(dhw)))+geom_point(shape=21)+scale_fill_viridis(option="magma"),
#  ggplot(pc.df, aes(PC1, PC4, fill=log(gravity)))+geom_point(shape=21)+scale_fill_viridis())


# 3)find upper/lower hulls

n_rows <- round((nrow(pc.df)/100)*10, 0)
row_n <- (nrow(pc.df)-n_rows):nrow(pc.df)

dhw_low <- cbind(pc.df[pc.df$dhw==0,], type="lowest")  #use all the zeros (around 280)
dhw_high <-cbind(pc.df[pc.df$dhw %in% sort(pc.df$dhw)[row_n], ], type="highest")
dhw_sub <- rbind(dhw_low, dhw_high)
dhw_hulls <- rbind(dhw_low[chull(dhw_low[,c("PC2", "PC3")]), ], dhw_high[chull(dhw_high[,c("PC2", "PC3")]), ])

grav_high <-cbind(pc.df[pc.df$gravity %in% sort(pc.df$gravity)[row_n], ], type="highest")
grav_low <-cbind(pc.df[pc.df$gravity %in% sort(pc.df$gravity, decreasing=T)[row_n], ], type="lowest")
grav_sub <- rbind(grav_low, grav_high)
grav_hulls <- rbind(grav_low[chull(grav_low[,c("PC1", "PC4")]), ], grav_high[chull(grav_high[,c("PC1", "PC4")]), ])

# 4) plot

display.brewer.pal(n = 8, name = 'Dark2')
brewer.pal(n = 8, name = "Dark2")
display.brewer.pal(n = 8, name = "Paired")


dhw_sub$type2 <- ifelse(dhw_sub$type=="highest", "High DHW", "Low DHW")

dhw_plot <- plot_grid(
  ggplot(data=dhw_sub)+geom_density(aes(x=PC2, fill=type), alpha=0.6, linewidth=0.2)+
    guides(fill="none")+
    xlim(c(-6,7))+
    theme_void()+theme(plot.margin=margin(0,0,0,0))+
    scale_fill_viridis(discrete = TRUE, option="cividis", direction=1, begin=0.8, end=0.2),
  #scale_fill_manual(values=c("#D95F02", "#1B9E77")),
  ggplot()+theme_void(),
  ggplot()+theme_void(),
  ggplot()+theme_void(),
  ggplot()+theme_void(),
  ggplot()+theme_void(),
  ggplot()+
    geom_point(data=pc.df, aes(PC2, PC3), col="grey", size=0.5)+
    stat_density_2d(data=dhw_sub, aes(PC2, PC3, col=type2, fill=type2), breaks=c(0.01), geom="polygon", alpha=0.5)+
    #geom_polygon(data=grav_hulls, aes(PC1, PC4, colour=type), fill=NA, linewidth=0.5, linetype="dotted")+
    geom_point(data=dhw_sub, aes(PC2, PC3,fill=type2), size=0.8, shape=21, stroke=0.2)+
    geom_segment(data=vecs, aes(x=0, y=0, xend=PC2*9, yend=PC3*5))+
    geom_label_repel(data=vecs, aes(PC2*9, PC3*5, label=n), fontface="bold",size=2.5, force=0.0008, label.padding=unit(0.3, "mm"))+
    xlim(c(-6,7))+
    ylim(c(-5,6))+
    # scale_fill_viridis(discrete = TRUE) +
    #scale_color_manual(values=c( "#D95F02","#1B9E77"))+
    #scale_fill_manual(values=c("#D95F02", "#1B9E77"))+
    scale_colour_viridis(discrete = TRUE, option="cividis", direction=1, begin=0.8, end=0.2)+
    scale_fill_viridis(discrete = TRUE, option="cividis", direction=1, begin=0.8, end=0.2)+
    #labs(x=paste("PC2 (", vars[2], "%)", sep=""), y=paste("PC3 (", vars[3], "%)", sep=""))+
    #scale_alpha(range=c(0.4, 0.6))+
    theme_bw()+theme(panel.grid=element_blank(), legend.title=element_blank(), legend.key.width=unit(2, "mm"), legend.key.height=unit(1, "mm"), legend.position=c(0.18, 0.92), legend.background=element_blank(), plot.background=element_blank()), 
  ggplot()+theme_void(),
  ggplot(data=dhw_sub)+geom_density(aes(x=PC3, fill=type), alpha=0.6, linewidth=0.2)+
    guides(fill="none")+
    xlim(c(-5,6))+
    coord_flip()+
    theme_void()+
    scale_fill_viridis(discrete = TRUE, option="cividis", direction=1, begin=0.8, end=0.2),
  #scale_fill_manual(values=c("#D95F02", "#1B9E77")),
  ncol=3, rel_heights=c(0.3,-0.12, 1), rel_widths=c(1, -0.12, 0.3), align="hv", axis="tblr")




grav_sub$type2 <- ifelse(grav_sub$type=="highest", "High Gravity", "Low Gravity")

grav_plot <- plot_grid(
  ggplot(data=grav_sub)+geom_density(aes(x=PC1, fill=type), alpha=0.6, linewidth=0.2)+
    guides(fill="none")+
    xlim(c(-11,6.5))+
    theme_void()+theme(plot.margin=margin(0,0,0,0))+
    scale_fill_viridis(discrete = TRUE,  direction=-1, begin=0, end=0.5),
  #scale_fill_manual(values=c("#D95F02", "#1B9E77")),
  ggplot()+theme_void(),
  ggplot()+theme_void(),
  ggplot()+theme_void(),
  ggplot()+theme_void(),
  ggplot()+theme_void(),
  ggplot()+
    geom_point(data=pc.df, aes(PC1, PC4), col="grey", size=0.5)+
    stat_density_2d(data=grav_sub, aes(PC1, PC4, col=type2, fill=type2), breaks=c(0.02), geom="polygon", alpha=0.5)+
    #geom_polygon(data=grav_hulls, aes(PC1, PC4, colour=type), fill=NA, linewidth=0.5, linetype="dotted")+
    geom_point(data=grav_sub, aes(PC1, PC4,fill=type2), size=0.8, shape=21, stroke=0.2)+
    geom_segment(data=vecs, aes(x=0, y=0, xend=PC1*9, yend=PC4*5))+
    geom_label_repel(data=vecs, aes(PC1*9, PC4*5, label=n), fontface="bold",size=2.5, force=0.0008, label.padding=unit(0.3, "mm"))+
    xlim(c(-11,6.5))+
    ylim(c(-3,5.5))+
    scale_colour_viridis(discrete = TRUE,  direction=-1, begin=0, end=0.5)+
    scale_fill_viridis(discrete = TRUE,  direction=-1, begin=0, end=0.5) +
    #scale_color_manual(values=c( "#D95F02","#1B9E77"))+
    #scale_fill_manual(values=c("#D95F02", "#1B9E77"))+
    #scale_alpha(range=c(0.4, 0.6))+
    #labs(x=paste("PC1 (", vars[1], "%)", sep=""), y=paste("PC4 (", vars[4], "%)", sep=""))+
    theme_bw()+theme(panel.grid=element_blank(), legend.title=element_blank(), legend.key.width=unit(2, "mm"), legend.key.height=unit(1, "mm"), legend.position=c(0.2, 0.92), legend.background=element_blank(), plot.background=element_blank()), 
  ggplot()+theme_void(),
  ggplot(data=grav_sub)+geom_density(aes(x=PC4, fill=type), alpha=0.6, linewidth=0.2)+
    guides(fill="none")+
    coord_flip()+
    xlim(c(-3,5.5))+
    theme_void()+
    scale_fill_viridis(discrete = TRUE,  direction=-1, begin=0, end=0.5),
  #scale_fill_manual(values=c("#D95F02", "#1B9E77")),
  ncol=3, rel_heights=c(0.3,-0.12, 1), rel_widths=c(1, -0.12, 0.3), align="hv", axis="tblr")

fig3ab <- plot_grid(NULL, dhw_plot, NULL, grav_plot, NULL, labels=c("", "a","", "b", ""), nrow=1, rel_widths = c(0.02,1,0.1,1,0.01))

#ggsave("output/fig3ab.png", fig_3, dev = "png", width = 10, height = 4.8)
fig3ab

}




combine_fig_3 <- function(fig_3ab, fig_3c) {
  
  fig_3cb <- plot_grid(fig_3c, NULL, rel_widths = c(1,0.08))
  fig3 <- plot_grid(
    fig_3ab, fig_3cb, 
    labels = c("", "c"), nrow = 2,
    rel_heights = c(1.3,0.7)
  )
  
  ggsave("output/Fig_3.png", fig3, width = 30, height = 20, units = "cm", dpi = 320)
  
  
  
}







  
