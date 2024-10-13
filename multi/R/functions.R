benthic_functions_estimation <- function(rls, traits, geo, pred, rls_info, rls_sites, indivs) {
  
  library(readxl) ; library(sp) ; library(ggforce) ; library(tidyverse); library(patchwork) ; library(scales)
  library(reshape2) ; library(cowplot) ; library(ggrepel) ; library(fishualize) ; library(viridis) ; library(ggcorrplot)
  
  ### AIM: CALCULATE BENTHIC FUNCTIONS USING SURVEY DATA
  rls$SurveyDate1 <- as.Date(rls$SurveyDate, format="%d/%m/%Y")
  rls <- rls[!rls$SurveyID==912340599,] # 0 % cover for all
  rls_cols <- rls_info$col
  names(rls_cols) <- rls_info$taxon
  
  rls$taxon_original <- rls$taxon
  rls$taxon <- ifelse(rls$taxon=="Rock" | rls$taxon=="coral rubble", "turf algae", as.character(rls$taxon))
  
  # macroalgae acts similar to fleshy? 
  rls$taxon <- ifelse(rls$taxon=="canopy forming macroalgae" | rls$taxon=="understory macroalgae", "fleshy algae", as.character(rls$taxon))
  
  ########################################## simulate sizes
  ##########################################################
  
  #  make df with rows = "individuals". give a size in cm
  #  simple model: assume surveys are 1 m2 quadrats. everything 10x10cm. 
  
  rls$N <- round(rls$cover) # % cover = number of "individuals"
  sites <- unique(rls$SurveyID)
  ntot<-length(unique(sites))
  coords <- expand.grid(c(1:10), c(1:10))  # for plotting
  
  surveys <- length(unique(indivs$SurveyID))
  surveys_dataset <- indivs %>% group_by(SurveyID) %>% group_split()
  
  test <- surveys_dataset[[2]] %>% group_by(taxon) %>% arrange(taxon, x, y)
  check_XY_difference <- function(x, y) { c(FALSE, abs(diff(x)) == 1 | abs(diff(y)) == 1) }
  
  # Sum the diam_cm values based on the condition within the same taxon
  surveys_size_arranged <- vector(mode = "list", length = surveys)
  for (i in 1:surveys) {
    surveys_size_arranged[[i]] <- surveys_dataset[[i]] %>%
      group_by(taxon) %>%
      mutate(group_id = cumsum(!check_XY_difference(x, y))) %>%
      group_by(taxon, group_id) %>%
      summarise(SurveyID = first(SurveyID),
                X = mean(x),
                Y = mean(y),
                diam_cm = sum(diam_cm)) }
  
  indivs = surveys_size_arranged %>% bind_rows()
  
  hemisphere<-function(d, tb, ...){
    r <- (d/2)
    SA<-(2*pi*r^2)
    TB<-SA*tb
    VOL<-((2/3)*pi*r^3)
    FULL <- VOL # No branch space
    data.frame(SA, TB, VOL, FULL)}
  
  flat<-function(d, th, tb, ...){
    r <- (d/2)
    SA<-pi*(r^2)
    TB<-SA*tb
    VOL<-SA*th
    FULL <- VOL # No branch space
    data.frame(SA, TB, VOL, FULL)}
  
  branches<-function(d, bh, br, bpa, tb, ...){ 
    r <- (d/2)
    SA<-(pi*(r^2)*(bpa*(2*pi*br*bh)+(pi*(br^2))))
    TB<-SA*tb
    VOL<-(pi*(r^2)*(bpa*(pi*(br^2)*bh)))
    FULL <- (pi*((d/2)^2)) * bh
    data.frame(SA,TB, VOL, FULL)}
  
  cone<-function(d, bh, th, tb, ...){
    r <- d/2
    #SA <- pi*r*(r+sqrt((bh^2)+(r^2)))
    SA <- pi * r * sqrt((r^2)+(bh^2)) #pi*r*l
    TB<-SA*tb
    VOL<-th*(SA) 
    FULL <- (pi*(r^2)) * bh
    data.frame(SA,TB, VOL, FULL)}
  
  ############################################# add geometry
  ##########################################################
  
  # use formula to get VOL, SA, FULL, TB
  # use growth rates to get geometry after 1 year
  
  growth <- subset(traits, trait_name=="Growth rate")
  growth <- aggregate(value~rls_group, growth, mean) 
  growth # means for now. / perhaps randomly sample from dist?  
  indivs$growth <- growth$value[match(indivs$taxon, growth$rls_group)]
  indivs$diam_cm_next <- indivs$diam_cm + (indivs$growth /10 )
  
  geo$bh[geo$rls == "Branching"] = 5
  
  # set formula
  indivs$fun <- rls_info$geometry[match(indivs$taxon, rls_info$taxon)]
  use <- indivs[complete.cases(indivs), ] # removes all without trait data 
  indivs = indivs %>% data.frame()
  use = use %>% data.frame()
  
  t0 <- NULL
  t1 <- NULL
  ntot <- nrow(use)
  
  for (i in 1:ntot){
    funct <- use[i, "fun"]
    gp <- use[i,"taxon"]
    pars <- geo[geo$rls==gp, c("br","bh", "bpa", "th", "tb")]
    pars$d <- use[i,"diam_cm"]
    t0 <- rbind(t0, unlist(do.call(Map, c(f= as.name(funct), pars))))
    pars$d <- use[i,"diam_cm_next"]
    t1 <- rbind(t1, unlist(do.call(Map, c(f= as.name(funct), pars))))
    print(round(which(c(1:ntot)==i)/ntot*100, 2))
  }
  colnames(t1) <- paste(colnames(t1), "_next", sep="")
  morphs <- cbind(use, t0, t1)
  
  # Rugosity corrected from Husband et al. 2021 (https://doi.org/10.1007/s00338-022-02247-6)
  morphs$Rugosity[morphs$taxon %in% c("Halimeda", "fleshy algae", "turf algae")] <- 1
  morphs$Rugosity[morphs$taxon %in% c("Encrusting", "coralline algae")] <- 1.6
  morphs$Rugosity[morphs$taxon == "Hemispherical"] <- 1.8
  morphs$Rugosity[morphs$taxon == "Tabular"] <- 2.0
  morphs$Rugosity[morphs$taxon == "Corymbose"] <- 2.3
  morphs$Rugosity[morphs$taxon == "Branching"] <- 3.0
  morphs$Rugosity[morphs$taxon == "Laminar"] <- 3.5
  
  # Rugosity per SurveyID
  morphs_Rugo <- morphs %>% mutate(Rugosity_reef = (Rugosity * diam_cm)) %>% group_by(SurveyID) %>% 
    summarise(Rugosity_reef = mean(Rugosity_reef)/10)
  
  # Merge Rugosity to previous data
  morphs <- merge(morphs, morphs_Rugo) %>% dplyr::select(-Rugosity) %>% rename(Rugosity = Rugosity_reef)
  
  # merge additional traits with morphs data 
  trait_avs <- aggregate(value~rls_group+trait_name, traits, mean) # or sample?
  
  ######################## (1) carbonate production & storage
  # change in volume and skeletal density
  sd_dat <- subset(trait_avs, trait_name=="Skeletal density")
  morphs$skel_den <- sd_dat$value[match(morphs$taxon, sd_dat$rls_group)]
  morphs$Storage_kg <- (morphs$VOL*morphs$skel_den )/1000 # (cm3 * gcm3 = g )
  morphs$Accretion_kgyr <- ((morphs$VOL_next-morphs$VOL)*morphs$skel_den)/1000 
  morphs$Storage_kg[is.na(morphs$Storage_kg)] <- 0  # others contribute 0
  morphs$Accretion_kgyr[is.na(morphs$Accretion_kgyr)] <- 0 # others contribute 0
  
  ######################## (2) branch space & rugosity
  morphs$BranchSpace_cm3<- (morphs$FULL - morphs$VOL) 
  remove <- c("fleshy algae", "turf algae", "Halimeda") # don't contribute?
  morphs$BranchSpace_cm3[morphs$taxon %in% remove ]<-NA
  morphs$BranchSpace_cm3[is.na(morphs$BranchSpace_cm3)] <- 0 # others contribute 0
  
  ######################## (3) photosynthesis & calcification
  gpp_dat <- subset(trait_avs, trait_name=="Gross photosynthesis")
  morphs$gpp <- gpp_dat$value[match(morphs$taxon, gpp_dat$rls_group)] #ugcm2hr
  morphs$GPP_ghr <- ( morphs$SA * morphs$gpp) / (1000^2)
  
  calc_dat <- subset(trait_avs, trait_name=="Calcification rate")
  morphs$calc <- calc_dat$value[match(morphs$taxon, calc_dat$rls_group)] #ugcm2hr
  morphs$Calc_ghr <- ( morphs$SA * morphs$calc) / (1000^2)
  morphs$Calc_ghr[is.na(morphs$Calc_ghr)] <- 0 # others contribute 0
  
  ######################## (4) organic production & storage
  # corals by SA algae by vol
  morphs$calc_status <- rls_info$calcification[match(morphs$taxon, rls_info$taxon)]
  morphs$OrgGrowth_gyr <- (morphs$TB_next - morphs$TB) / 1000 # corals
  morphs$OrgGrowth_gyr <- ifelse(morphs$calc_status=="non_calcifying", 
                                 (((morphs$VOL_next - morphs$VOL)*0.07) - 15) / 1000, morphs$OrgGrowth_gyr)
  # based on https://doi.org/10.1515/BOT.2002.063
  morphs$BiomassStand_g <- (morphs$TB) / 1000
  morphs$BiomassStand_g  <- ifelse(morphs$calc_status=="non_calcifying", 
                                   ((morphs$VOL*0.07) - 15) / 1000, morphs$BiomassStand_g) 
  # based on https://doi.org/10.1515/BOT.2002.063
  morphs$mean_sst <- pred$mean_sst_5year[match(morphs$SurveyID, pred$SurveyID)]
  sum(is.na(morphs$mean_sst)) 
  morphs[is.na(morphs$mean_sst),]
  morphs[is.na(morphs$mean_sst),"mean_sst"] <- 25 # Change later for solitary islands
  
  algae_temp_factor <- function(toptprime, erprime, ei, temp_c) {
    topt <- 288.15 + plogis(toptprime) * 30
    er <- ei * plogis(erprime)
    boltz <- 8.62e-5
    ei * exp(
      er * ((temp_c + 273.15) - topt) / ((temp_c + 273.15) * boltz * topt)
    ) /
      (ei - er * (1 - exp(
        ei * ((temp_c + 273.15) - topt) / ((temp_c + 273.15) * boltz * topt))
      ))
  }
  coral_temp_factor <- function(g_max = 1, rho, temp_shift_c, ln_sigma, temp_c) {
    g_max * exp( -1 * exp(rho * (temp_c - temp_shift_c)) -
                   exp(ln_sigma) * (temp_c - temp_shift_c)^2
    )
  }
  
  morphs$coral_sst_factor <- coral_temp_factor( g_max = 1, rho = 1.022, temp_shift_c = 30.197, ln_sigma = -4.518, temp_c = morphs$mean_sst)
  morphs$algae_sst_factor <- algae_temp_factor(toptprime = 0.21, erprime = -1.31, ei = 3, temp_c = morphs$mean_sst  )
  
  algae <- c("coralline algae", "fleshy algae", "Halimeda", "turf algae")
  morphs$coral_sst_factor[morphs$taxon %in% algae] <- 1
  morphs$algae_sst_factor[!morphs$taxon %in% algae] <- 1
  rates <- c("OrgGrowth_gyr", "GPP_ghr", "Accretion_kgyr", "Calc_ghr")
  
  morphs[,rates] <- morphs[,rates] * morphs$coral_sst_factor * morphs$algae_sst_factor
  
  ###################################### taxa contributions
  ##########################################################
  # first sum functions by site + taxon
  functions <- c("Storage_kg", "Accretion_kgyr", "GPP_ghr", "Calc_ghr", "BranchSpace_cm3", "OrgGrowth_gyr", "BiomassStand_g")
  fun_taxa <- aggregate(.~SurveyID + taxon, morphs[,c("SurveyID", "taxon", functions)], sum)
  
  plotdat <- melt(fun_taxa, id=c("SurveyID", "taxon"))
  plotdat <- plotdat[plotdat$SurveyID %in% unique(plotdat$SurveyID)[1:70],]
  taxa_plot <- ggplot(data=plotdat)+
    geom_bar(aes(y=as.factor(SurveyID), x=value, fill=taxon), stat="identity")+
    scale_fill_manual(values=rls_cols)+
    facet_wrap(~variable, scales="free_x", nrow=1)+
    theme_bw()+ylab("SurveyID")+
    theme(axis.text.y=element_text(size=2), strip.background=element_blank())
  
  ##################################### site-level functions
  ##########################################################
  sites     <- data.frame(SurveyID=unique(rls$SurveyID))
  #rls_sites <- get(load("benthic_functions_MM/data/data_RLS/RLS_sitesInfos.RData"))
  add_cols  <- c("SiteCode", "SiteLatitude", "SiteLongitude", "Country","SiteName", "Location", "Depth", "SurveyDate")
  sites[,add_cols] <- rls_sites[match(sites$SurveyID, rls_sites$SurveyID), add_cols]
  
  rls$broad_cat <- rls_info$broad_cat[match(rls$taxon, rls_info$taxon)]
  cover         <- aggregate(cover~SurveyID+broad_cat, rls, sum)
  cover         <- dcast(cover, SurveyID~broad_cat, value.var="cover")
  cover$total   <- rowSums(subset(cover, select=-c(SurveyID)))
  cover$flat    <- cover$total-cover$coral
  sites         <- merge(sites, cover, all.x=F, all.y=F)
  
  funs              <- aggregate(.~ SurveyID, morphs[,c("SurveyID",functions)], sum)
  sites[,functions] <- funs[match(sites$SurveyID, funs$SurveyID), functions]
  
  funs  <- merge(funs, morphs_Rugo, by = "SurveyID")  %>% rename(Rugosity = Rugosity_reef)
  sites <- merge(sites, morphs_Rugo, by = "SurveyID") %>% rename(Rugosity = Rugosity_reef)
  
  sites.all <- na.omit(sites)
  pca       <- prcomp(sqrt(sites.all[,c(functions, "Rugosity")]), scale=T, center=T)
  sites.all[,c("pc1", "pc2")] <- pca$x[,1:2]
  vecs      <- data.frame(pc1=pca$rotation[,1], pc2=pca$rotation[,2], 
                          labels = rownames(pca$rotation), pc1b=pca$rotation[,1], pc2b=pca$rotation[,2])
  vecs$pc2b[vecs$labels=="Calc_ghr"] <- vecs$pc2b[vecs$labels=="Calc_ghr"]+0.05
  vecs$pc2b[vecs$labels=="GPP_ghr"]  <- vecs$pc2b[vecs$labels=="GPP_ghr"]+0.05
  vars <- round((pca$sdev^2 / sum(pca$sdev^2)), 3)*100
  
  data_corrected_size = melt(sites.all[,c("SurveyID",functions, "Rugosity")], id.var="SurveyID")
  xlsx::write.xlsx(data_corrected_size, "output/data_corrected_size.xlsx")
  

}


plot_figure2 <- function(Fish_Functions, Benthos_Functions, case_study) {

  library(readxl) ; library(sp) ; library(ggforce) ; library(tidyverse); library(patchwork) ; library(scales)
  library(reshape2) ; library(cowplot) ; library(ggrepel) ; library(fishualize) ; library(viridis) ; library(ggcorrplot)
  

fish  <- Fish_Functions
benth <- Benthos_Functions
dataset_functions <- merge(Benthos_Functions, Fish_Functions, by = "SurveyID")

# Adding case study as the same format of RLS dataset from Mike
case_study <- case_study %>% 
  select(c(2:5, 7, 9, 11, 13, 15, 19, 22:26, 28, 31)) %>% drop_na() %>% 
  mutate(SurveyID = paste(Location, Site, Year, sep = "_")) %>% 
  select(-c(Year, Location, Site)) %>% column_to_rownames("SurveyID") %>% 
  select(8:9, 11, 10, 13, 12, 14, 3, 6, 4, 1, 2, 7, 5)
colnames(case_study) <- c("herb", "plank", "prod","pisc",  "exN", "exP", "turn", 
                          "GPP_ghr", "OrgGrowth_gyr",  "Calc_ghr", "Storage_kg", 
                          "Accretion_kgyr",  "Rugosity", "BranchSpace_cm3")

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
case_study <- case_study[,c(f.fun, b.fun)]+1
log.all <- rbind(log.all, case_study)
#log.all[log.all==0] <- NA
#log.all <- log.all[complete.cases(log.all),]
log.all <- log(log.all)
nrow(log.all)

pca <- prcomp(log.all, center=T, scale=T)
# biplot(pca)

# ggcor

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
p1

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
p2

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
p3

plot_grid(p1, p2, p3, nrow=1)

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
p4


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
p5

(Figure_2 <- plot_grid(
  plot_grid(p1, p2, p3, nrow=1, labels=c("a", "b", "c")), 
  plot_grid(p4, p5, rel_widths=c(1, 0.8), labels=c("d","e")), 
  ncol=1))

ggsave("output/fig2.png", Figure_2,  width = 10, height = 7)

}




plot_figure1 <- function(Benthos_Functions, Fish_Functions, rls, data_ts_coord) {
  
  
  ########################
  ####### FIGURE 1 #######
  ########################
  
  library(readxl) ; library(sp) ; library(ggforce) ; library(tidyverse); library(patchwork) ; library(scales)
  library(reshape2) ; library(cowplot) ; library(ggrepel) ; library(fishualize) ; library(viridis) ; library(ggcorrplot); library(sf)
  
  ### Loading projections
  load(url("https://github.com/valentinitnelav/RandomScripts/blob/master/NaturalEarth.RData?raw=true"))
  PROJ              <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" 
  
  ### Viz Function distribution across the world (Benthos and Fish)
  # Building dataframes
  
  dataset_functions <- merge(Benthos_Functions, Fish_Functions, by = "SurveyID")
  
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
  
  dataset_functions <- dataset_functions %>% mutate(functions_log10 = log10(functions+0.000001),
                                                    functions_log10_rounded = round(log10(functions+0.000001), 1))
  
  dataset_functions %>% group_by(Label) %>% summarise (min = min((functions_log10)), max = max((functions_log10)))
  
  # Figure 1B
  (Figure_1B <- (ggplot(data = dataset_functions, aes(x=functions, fill = Label)) +
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
  (Fig_1B[[K]] <- ggplot(data = dataset_functions_split[[K]], aes(x=functions, fill = Label)) +
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
            axis.text = element_text(size = 12),
            axis.title.x = element_text(size = 14),
            plot.margin = unit(c(0.5,0.5,1,0.5), "cm"),
            plot.title = element_text(size = 10, face = "bold"),
            axis.line = element_line(colour = "black")))
  
  K = 2
  (Fig_1B[[K]] <- ggplot(data = dataset_functions_split[[K]], aes(x=functions, fill = Label)) +
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
            axis.text = element_text(size = 12),
            axis.title.x = element_text(size = 14),
            plot.margin = unit(c(0.5,0.5,1,0.5), "cm"),
            plot.title = element_text(size = 10, face = "bold"),
            axis.line = element_line(colour = "black")))
  
  K = 3
  (Fig_1B[[K]] <- ggplot(data = dataset_functions_split[[K]], aes(x=functions, fill = Label)) +
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
            axis.text = element_text(size = 12),
            axis.title.x = element_text(size = 14),
            plot.margin = unit(c(0.5,0.5,1,0.5), "cm"),
            plot.title = element_text(size = 10, face = "bold"),
            axis.line = element_line(colour = "black")))
  
  K = 4
  (Fig_1B[[K]] <- ggplot(data = dataset_functions_split[[K]], aes(x=functions, fill = Label)) +
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
            axis.text = element_text(size = 12),
            axis.title.x = element_text(size = 14),
            plot.margin = unit(c(0.5,0.5,1,0.5), "cm"),
            plot.title = element_text(size = 10, face = "bold"),
            axis.line = element_line(colour = "black")))
  
  K = 5
  (Fig_1B[[K]] <- ggplot(data = dataset_functions_split[[K]], aes(x=functions, fill = Label)) +
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
            axis.text = element_text(size = 12),
            axis.title.x = element_text(size = 14),
            plot.margin = unit(c(0.5,0.5,1,0.5), "cm"),
            plot.title = element_text(size = 10, face = "bold"),
            axis.line = element_line(colour = "black")))
  
  K = 6
  (Fig_1B[[K]] <- ggplot(data = dataset_functions_split[[K]], aes(x=functions, fill = Label)) +
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
            axis.text = element_text(size = 12),
            axis.title.x = element_text(size = 14),
            plot.margin = unit(c(0.5,0.5,1,0.5), "cm"),
            plot.title = element_text(size = 10, face = "bold"),
            axis.line = element_line(colour = "black")))
  
  K = 7
  (Fig_1B[[K]] <- ggplot(data = dataset_functions_split[[K]], aes(x=functions, fill = Label)) +
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
            axis.text = element_text(size = 12),
            axis.title.x = element_text(size = 14),
            plot.margin = unit(c(0.5,0.5,1,0.5), "cm"),
            plot.title = element_text(size = 10, face = "bold"),
            axis.line = element_line(colour = "black")))
  
  K = 8
  (Fig_1B[[K]] <- ggplot(data = dataset_functions_split[[K]], aes(x=functions, fill = Label)) +
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
            axis.text = element_text(size = 12),
            axis.title.x = element_text(size = 14),
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
            plot.title = element_text(size = 10, face = "bold"),
            axis.line = element_line(colour = "black")))
  
  K = 9
  (Fig_1B[[K]] <- ggplot(data = dataset_functions_split[[K]], aes(x=functions, fill = Label)) +
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
            axis.text = element_text(size = 12),
            axis.title.x = element_text(size = 14),
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
            plot.title = element_text(size = 10, face = "bold"),
            axis.line = element_line(colour = "black")))
  
  K = 10
  (Fig_1B[[K]] <- ggplot(data = dataset_functions_split[[K]], aes(x=functions, fill = Label)) +
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
            axis.text = element_text(size = 12),
            axis.title.x = element_text(size = 14),
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
            plot.title = element_text(size = 10, face = "bold"),
            axis.line = element_line(colour = "black")))
  
  K = 11
  (Fig_1B[[K]] <- ggplot(data = dataset_functions_split[[K]], aes(x=functions, fill = Label)) +
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
            axis.text = element_text(size = 12),
            axis.title.x = element_text(size = 14),
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
            plot.title = element_text(size = 10, face = "bold"),
            axis.line = element_line(colour = "black")))
  
  K = 12
  (Fig_1B[[K]] <- ggplot(data = dataset_functions_split[[K]], aes(x=functions, fill = Label)) +
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
            axis.text = element_text(size = 12),
            axis.title.x = element_text(size = 14),
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
            plot.title = element_text(size = 10, face = "bold"),
            axis.line = element_line(colour = "black")))
  
  K = 13
  (Fig_1B[[K]] <- ggplot(data = dataset_functions_split[[K]], aes(x=functions, fill = Label)) +
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
            axis.text = element_text(size = 12),
            axis.title.x = element_text(size = 14),
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
            plot.title = element_text(size = 10, face = "bold"),
            axis.line = element_line(colour = "black")))
  
  K = 14
  (Fig_1B[[K]] <- ggplot(data = dataset_functions_split[[K]], aes(x=functions, fill = Label)) +
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
            axis.text = element_text(size = 12),
            axis.title.x = element_text(size = 14),
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
            plot.title = element_text(size = 10, face = "bold"),
            axis.line = element_line(colour = "black")))
  
  Figure_1B_top <- (Fig_1B[[1]] + Fig_1B[[2]] + Fig_1B[[3]] + Fig_1B[[4]] + 
                      Fig_1B[[5]] + Fig_1B[[6]] + Fig_1B[[7]]) + plot_layout(nrow = 1)
  Figure_1B_bot <- (Fig_1B[[8]] + Fig_1B[[9]] + Fig_1B[[10]] + Fig_1B[[11]] + Fig_1B[[12]] + 
                      Fig_1B[[13]] + Fig_1B[[14]]) + plot_layout(nrow = 1)
  (Figure_1B <- Figure_1B_top / Figure_1B_bot)
  
  # Setting map parameters ####
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
  
  spatial_points_rls <- SpatialPoints(cbind(rls$SiteLongitude , rls$SiteLatitude), proj4string = CRS("+proj=longlat +datum=WGS84"))
  data_Spatial_proj = spTransform(spatial_points_rls, CRS(PROJ)) %>% 
    as.data.frame() %>% rename(Long_Robin = coords.x1, Lat_Robin = coords.x2) 
  rls = rls %>% cbind(.,data_Spatial_proj)
  
  spatial_points_ts <- SpatialPoints(cbind(data_ts_coord$Long , data_ts_coord$Lat), proj4string = CRS("+proj=longlat +datum=WGS84"))
  data_ts_coord_proj <- spTransform(spatial_points_ts, CRS(PROJ)) %>% as.data.frame() %>% 
    rename(Long_Robin = coords.x1, Lat_Robin = coords.x2) 
  data_ts_coord = data_ts_coord %>% cbind(.,data_ts_coord_proj) %>% dplyr::filter(Lat <= 80)
  
  (Figure_1A = ggplot(data = rls, aes(x = Long_Robin, y = Lat_Robin)) +
      geom_polygon(data=NE_countries_rob, aes(long,lat, group=group), colour="gray80", fill="gray80", linewidth = 0.25) +
      geom_polygon(data=NE_box_rob, aes(x=long, y=lat), colour="grey50", fill="transparent", size = 0.25) +
      geom_path(data=NE_graticules_rob, aes(long, lat, group=group), linetype="dotted", color="grey50", size = 0.25) +
      geom_text(data = lbl.Y.prj, aes(x = coords.x1, y = coords.x2, label = lbl), color="grey50", size=2) +
      geom_text(data = lbl.X.prj, aes(x = coords.x1, y = coords.x2, label = lbl), color="grey50", size=2) +
      geom_point(show.legend = F, size = 2, fill = "#9966CC", shape = 21) +
      coord_fixed(ratio = 1) + theme_void() + theme(legend.position="bottom"))
  

  
  ggsave("output/fig1a.png", Figure_1A , dev = "png", width = 10, height = 4)
  ggsave("output/fig1b.png", Figure_1B , dev = "png", width = 10, height = 4)
  

}




time_series <- function(ts_benthos, ts_fish, Seychelles_data, Moorea_data, rls_info, traits, geo, predictors_ts, Fish_functions, case_study, Benthos_Functions, Fish_Functions, pred) {
  
  library(readxl) ; library(sp) ; library(ggforce) ; library(tidyverse); library(patchwork) ; library(scales)
  library(reshape2) ; library(cowplot) ; library(ggrepel) ; library(fishualize) ; library(viridis) ; library(ggcorrplot); library(sf)
  
  
  ###########################
  ####### TIME SERIES #######
  ###########################
  
  # Key sites – Checking data status
  ts_benthos_keysites <- ts_benthos %>% dplyr::filter(Location %in% c("Cousin", "Mahe", "Praslin", "St Anne", "Moorea", "Tetiaroa", "Hoga"))
  table(ts_benthos$Location)
  
  ts_fish_keysites <- ts_fish %>% dplyr::filter(Location %in% c("Cousin", "Mahe", "Praslin", "St Anne", "Moorea", "Tetiaroa", "Hoga"))
  table(ts_fish_keysites$Location) # Cousin and Mahe coming from Nick – Fish team already got the data.
  # Everything tracks
  
  # Look at benthos resolution
  table(ts_benthos_keysites$Genus, ts_benthos_keysites$Location) # Seychelle is poor
  
  #### 1) Seychelles
  ts_benthos_seychelles <- ts_benthos %>% 
    dplyr::filter(Location %in% c("Mahe", "Cousin", "Praslin", "St Anne"))
  
  ts_benthos_seychelles$rls_broad_gr = NA
  ts_benthos_seychelles$rls_broad_gr[ts_benthos_seychelles$Group != "NA"] = ts_benthos_seychelles$Group[ts_benthos_seychelles$Group != "NA"]
  ts_benthos_seychelles$rls_broad_gr[ts_benthos_seychelles$Group == "NA"] = paste("Hard_Coral", 
                                                                                  ts_benthos_seychelles$Shape[ts_benthos_seychelles$Group == "NA"],
                                                                                  sep = "_")
  unique(ts_benthos_seychelles$Category)
  
  # Define all combination
  all_combinations <- expand.grid(
    Location = unique(ts_benthos_seychelles$Location),
    Site = unique(ts_benthos_seychelles$Site),
    Year = unique(ts_benthos_seychelles$Year),
    rls_broad_gr = unique(ts_benthos_seychelles$rls_broad_gr))
  
  # Perform the summarization
  ts_benthos_seychelles %>%
    group_by(Country, Location, Site, Depth, Year, Replicate, rls_broad_gr) %>%
    summarise(Cover = mean(Cover, na.rm = TRUE)) %>%
    right_join(all_combinations, by = c("Location", "Site", "Year", "rls_broad_gr")) %>%
    ungroup() %>% group_by(Country, Location, Site,        Year, rls_broad_gr) %>% 
    ggplot(aes(y = Cover, x = as.character(Year), fill = rls_broad_gr, group = )) + 
    geom_bar(stat="identity", width = 0.5) +
    facet_wrap(Replicate~Site) # Some replicates are not going up to 100%
  
  # Define the cover for each site vs. replicate
  total_cover_seychelles <- ts_benthos_seychelles %>%
    group_by(Country, Location, Site, Depth, Year, Replicate, rls_broad_gr) %>%
    summarise(Cover = mean(Cover, na.rm = TRUE)) %>%
    right_join(all_combinations, by = c("Location", "Site", "Year", "rls_broad_gr")) %>%
    ungroup() %>% group_by(Country, Location, Site, Replicate, Year) %>% 
    summarise(Cover_tot = sum(Cover)) %>% 
    mutate(std_cover = 100/Cover_tot)
  
  # With cover up to 100%
  # Fish team focusing on 1994, 2005 and 2008
  ts_benthos_seychelles %>%
    group_by(Country, Location, Site, Depth, Year, Replicate, rls_broad_gr) %>%
    summarise(Cover = mean(Cover, na.rm = TRUE)) %>%
    right_join(all_combinations, by = c("Location", "Site", "Year", "rls_broad_gr")) %>%
    right_join(total_cover_seychelles, by = c("Country", "Location", "Site", "Year", "Replicate")) %>% 
    mutate(Cover_corrected = Cover * std_cover) %>% 
    ungroup() %>% 
    right_join(all_combinations, by = c("Location", "Site", "Year", "rls_broad_gr")) %>%
    dplyr::filter(Year <= 2008) %>% 
    ggplot(aes(y = Cover_corrected, x = as.character(Year), fill = rls_broad_gr)) + 
    geom_bar(stat="identity", width = 0.5) +
    facet_wrap(Site~Replicate) # Some replicates are not going up to 100%
  
  ######
  Seychelles_data <- Seychelles_data %>% dplyr::select(-c(Country, habitat, `struc complexity`))
  Seychelles_data_melt <- melt(Seychelles_data, id = c("Site","Replicate","year","Transect length")) %>% 
    mutate(value = as.numeric(value)) 
  Seychelles_data_melt$value[is.na(Seychelles_data_melt$value)] = 0
  
  data_seychelles_barplot <- Seychelles_data_melt %>% group_by(Site, Replicate, year, variable) %>% summarise(value = mean(value)) %>% 
    group_by(Site, year, variable) %>% summarise(value = mean(value)) %>% drop_na() %>% 
    dplyr::filter(variable != "...74", variable != "...75") %>% 
    ggplot(aes(y = value * 10, x = as.character(year), fill = variable)) + 
    geom_bar(stat="identity", width = 0.5) +
    facet_wrap(~Site) # 5 sites in 2008 not having a 100% cover resolution (be carefful)
  
  # Attribute rls group
  Seychelles_dataset_site <- ts_benthos_seychelles %>%
    group_by(Country, Location, Site, Depth, Year, Replicate, rls_broad_gr) %>%
    summarise(Cover = mean(Cover, na.rm = TRUE)) %>%
    right_join(all_combinations, by = c("Location", "Site", "Year", "rls_broad_gr")) %>%
    right_join(total_cover_seychelles, by = c("Country", "Location", "Site", "Year", "Replicate")) %>% 
    mutate(Cover_corrected = Cover * std_cover) %>% 
    ungroup() %>% 
    right_join(all_combinations, by = c("Location", "Site", "Year", "rls_broad_gr")) %>%
    dplyr::filter(Year <= 2008) %>% drop_na()
  
  #### 2) Indonesia
  
  Hoga_dataset <- ts_benthos %>% dplyr::filter(Country == "Indonesia") 
  
  Hoga_dataset$rls_broad_gr = NA
  Hoga_dataset$rls_broad_gr[Hoga_dataset$Group != "NA"] = Hoga_dataset$Group[Hoga_dataset$Group != "NA"]
  Hoga_dataset$rls_broad_gr[Hoga_dataset$Group == "NA" | Hoga_dataset$Category == "Hard living coral"] = 
    paste("Hard_Coral", Hoga_dataset$Shape[Hoga_dataset$Group == "NA" | Hoga_dataset$Category == "Hard living coral"], sep = "_")
  Hoga_dataset$rls_broad_gr[Hoga_dataset$rls_broad_gr == "Hard_Coral_NA"] = 
    paste("Hard_Coral", Hoga_dataset$Genus[Hoga_dataset$rls_broad_gr == "Hard_Coral_NA"], sep = "_")
  Hoga_dataset$rls_broad_gr[Hoga_dataset$rls_broad_gr == "Hard_Coral_NA"] = 
    Hoga_dataset$Category[Hoga_dataset$rls_broad_gr == "Hard_Coral_NA"]
  
  total_cover_indonesia <- Hoga_dataset %>% group_by(Location, Site, Zone, Replicate, Year, rls_broad_gr) %>% 
    summarise(Cover = sum(Cover, na.rm = TRUE)) %>%
    ungroup() %>% group_by(Location, Site, Zone, Replicate, Year) %>% 
    summarise(Cover_tot = sum(Cover)) %>% 
    mutate(std_cover = 100/Cover_tot) 
  
  # With cover up to 100%
  Hoga_dataset_no_replicate <- Hoga_dataset %>%
    group_by(Location, Site, Zone, Replicate, Year, rls_broad_gr) %>%
    summarise(Cover = sum(Cover, na.rm = TRUE)) %>%
    right_join(total_cover_indonesia, by = c("Location", "Site", "Zone", "Year", "Replicate")) %>% 
    mutate(Cover_corrected = Cover * std_cover) 
  #ungroup() %>% 
  #group_by(Location, Site, Zone, Year, rls_broad_gr) %>% 
  #summarise(Cover_corrected = sum(Cover_corrected, na.rm = TRUE)/3)
  
  # Some gaps to fill again
  Hoga_fill_gap <- Hoga_dataset_no_replicate %>% 
    dplyr::filter(Zone == "Slope") %>% 
    group_by(Location, Site, Replicate, Zone, Year) %>% 
    summarise(Cover_tot = sum(Cover_corrected)) %>% 
    mutate(std_cover = 100/Cover_tot)
  
  Hoga_dataset_site <- Hoga_dataset_no_replicate %>%
    group_by(Location, Site, Zone, Replicate, Year, rls_broad_gr) %>%
    summarise(Cover_corrected = sum(Cover_corrected, na.rm = TRUE)) %>%
    right_join(Hoga_fill_gap, by = c("Location", "Site", "Replicate", "Zone", "Year")) %>% 
    mutate(Cover_corrected = Cover_corrected * std_cover) %>% 
    ungroup() 
  
  # B3 & KDS
  
  Hoga_dataset_site$rls_broad_gr[Hoga_dataset_site$rls_broad_gr %in% c("Hard living coral", "Hard_Coral_Acropora", "Hard_Coral_Branching")] = "Branching"
  Hoga_dataset_site$rls_broad_gr[Hoga_dataset_site$rls_broad_gr %in% c("Hard_Coral_Digitate")] = "Corymbose"
  Hoga_dataset_site$rls_broad_gr[Hoga_dataset_site$rls_broad_gr %in% c("Hard_Coral_Encrusting", "Hard_Coral_Leptastrea", "Hard_Coral_Montipora", "Hard_Coral_Leptoseris")] = "Encrusting"
  Hoga_dataset_site$rls_broad_gr[Hoga_dataset_site$rls_broad_gr %in% c("Hard_Coral_Favia", "Hard_Coral_Free-living", "Hard_Coral_Leptastrea", "Hard_Coral_Massive", "Hard_Coral_Submassive")] = "Massive"
  Hoga_dataset_site$rls_broad_gr[Hoga_dataset_site$rls_broad_gr %in% c("Hard_Coral_Foliose", "Hard_Coral_Mycedium", "Hard_Coral_Pachyseris", "Hard_Coral_Pavona")] = "Laminar"
  Hoga_dataset_site$rls_broad_gr[Hoga_dataset_site$rls_broad_gr %in% c("Hard_Coral_Millepora", "Hard_Coral_Tabular")] = "Tabular"
  Hoga_dataset_site$rls_broad_gr[Hoga_dataset_site$rls_broad_gr %in% c("Actiniaria", "Alcyonacea", "Hard_Coral_Galaxea", "Hard_Coral_Galaxea", "Hard_Coral_Tubastraea", "Hard_Coral_Goniopora")] = "soft"
  Hoga_dataset_site$rls_broad_gr[Hoga_dataset_site$rls_broad_gr %in% c("Algae", "Hard_Coral_Caulerpa", "Hard_Coral_Dictyota")] = "fleshy algae"
  Hoga_dataset_site$rls_broad_gr[Hoga_dataset_site$rls_broad_gr %in% c("Hard_Coral_Turbinaria")] = "understory macroalgae"
  Hoga_dataset_site$rls_broad_gr[Hoga_dataset_site$rls_broad_gr %in% c("Hard_Coral_Halimeda")] = "Halimeda"
  Hoga_dataset_site$rls_broad_gr[Hoga_dataset_site$rls_broad_gr %in% c("Coralline algae")] = "coralline algae"
  Hoga_dataset_site$rls_broad_gr[Hoga_dataset_site$rls_broad_gr %in% c("Crinoidea", "Bivalvia", "Tunicata", "Zoantharia", "Porifera", "Other fauna", "others")] = "others"
  Hoga_dataset_site$rls_broad_gr[Hoga_dataset_site$rls_broad_gr %in% c("Rubble", "Rock", "Turf algae", "Hard dead coral")] = "Rock and Rubble"
  Hoga_dataset_site$rls_broad_gr[Hoga_dataset_site$rls_broad_gr %in% c("Sand", "Silt")] = "Sand"
  
  Hoga_dataset$rls_broad_gr = NA
  Hoga_dataset$rls_broad_gr[Hoga_dataset$Group != "NA"] = Hoga_dataset$Group[Hoga_dataset$Group != "NA"]
  Hoga_dataset$rls_broad_gr[Hoga_dataset$Group == "NA" | Hoga_dataset$Category == "Hard living coral"] = 
    paste("Hard_Coral", Hoga_dataset$Shape[Hoga_dataset$Group == "NA" | Hoga_dataset$Category == "Hard living coral"], sep = "_")
  Hoga_dataset$rls_broad_gr[Hoga_dataset$rls_broad_gr == "Hard_Coral_NA"] = 
    paste("Hard_Coral", Hoga_dataset$Genus[Hoga_dataset$rls_broad_gr == "Hard_Coral_NA"], sep = "_")
  Hoga_dataset$rls_broad_gr[Hoga_dataset$rls_broad_gr == "Hard_Coral_NA"] = 
    Hoga_dataset$Category[Hoga_dataset$rls_broad_gr == "Hard_Coral_NA"]
  
  total_cover_indonesia <- Hoga_dataset %>% group_by(Location, Site, Zone, Replicate, Year, rls_broad_gr) %>% 
    summarise(Cover = sum(Cover, na.rm = TRUE)) %>%
    ungroup() %>% group_by(Location, Site, Zone, Replicate, Year) %>% 
    summarise(Cover_tot = sum(Cover)) %>% 
    mutate(std_cover = 100/Cover_tot) 
  
  # With cover up to 100%
  Hoga_dataset_no_replicate <- Hoga_dataset %>%
    group_by(Location, Site, Zone, Replicate, Year, rls_broad_gr) %>%
    summarise(Cover = sum(Cover, na.rm = TRUE)) %>%
    right_join(total_cover_indonesia, by = c("Location", "Site", "Zone", "Year", "Replicate")) %>% 
    mutate(Cover_corrected = Cover * std_cover) 
  #ungroup() %>% 
  #group_by(Location, Site, Zone, Year, rls_broad_gr) %>% 
  #summarise(Cover_corrected = sum(Cover_corrected, na.rm = TRUE)/3)
  
  # Some gaps to fill again
  Hoga_fill_gap <- Hoga_dataset_no_replicate %>% 
    dplyr::filter(Zone == "Slope") %>% 
    group_by(Location, Site, Replicate, Zone, Year) %>% 
    summarise(Cover_tot = sum(Cover_corrected)) %>% 
    mutate(std_cover = 100/Cover_tot)
  
  Hoga_dataset_site <- Hoga_dataset_no_replicate %>%
    group_by(Location, Site, Zone, Replicate, Year, rls_broad_gr) %>%
    summarise(Cover_corrected = sum(Cover_corrected, na.rm = TRUE)) %>%
    right_join(Hoga_fill_gap, by = c("Location", "Site", "Replicate", "Zone", "Year")) %>% 
    mutate(Cover_corrected = Cover_corrected * std_cover) %>% 
    ungroup() 
  
  Moorea_benthos_ts <- Moorea_data %>% 
    group_by(Campaign, Year, Season, `Marine Area`, Habitat, Transect, Substrate) %>% 
    summarise(proportion = mean(proportion)) %>% ungroup() %>% 
    dplyr::filter(Habitat == "Outer slope", Season == "Mar") %>% 
    dplyr::select(-c("Campaign", "Season", "Habitat"))
  
  Seychelles_dataset_site %>% 
    ggplot(aes(y = Cover_corrected, x = as.character(Year), fill = rls_broad_gr)) + 
    geom_bar(stat="identity", width = 0.5) +
    facet_wrap(Site~Replicate) 
  
  # Export
  
  Moorea_benthos_ts$rls_broad_gr = NA
  Moorea_benthos_ts$rls_broad_gr[Moorea_benthos_ts$Substrate %in% c("Millepora", "Napopora", "Psammocora")] = "Branching"
  Moorea_benthos_ts$rls_broad_gr[Moorea_benthos_ts$Substrate %in% c("Pocillopora")] = "Corymbose"
  Moorea_benthos_ts$rls_broad_gr[Moorea_benthos_ts$Substrate %in% c("Montipora", "Synarea", "Gardineroseris", "Fungia", "Lobophylia", "Corail non identifié")] = "Encrusting"
  Moorea_benthos_ts$rls_broad_gr[Moorea_benthos_ts$Substrate %in% c("Porites", "Favia", "Leptastrea", "Astreopora", "Cyphastrea", "Acanthastrea", "Herpolitha", "Dipsastrea")] = "Massive"
  Moorea_benthos_ts$rls_broad_gr[Moorea_benthos_ts$Substrate %in% c("Pavona", "Leptoseris")] = "Laminar"
  Moorea_benthos_ts$rls_broad_gr[Moorea_benthos_ts$Substrate %in% c("Acropora")] = "Tabular"
  Moorea_benthos_ts$rls_broad_gr[Moorea_benthos_ts$Substrate %in% c("Boodlea")] = "fleshy algae"
  Moorea_benthos_ts$rls_broad_gr[Moorea_benthos_ts$Substrate %in% c("Macroalgae", "Turbinaria")] = "understory macroalgae"
  Moorea_benthos_ts$rls_broad_gr[Moorea_benthos_ts$Substrate %in% c("Halimeda")] = "Halimeda"
  Moorea_benthos_ts$rls_broad_gr[Moorea_benthos_ts$Substrate %in% c("Other")] = "others"
  Moorea_benthos_ts$rls_broad_gr[Moorea_benthos_ts$Substrate %in% c("Rubble", "Dead coral", "Pavement", "Astreopora bleach", "Pocillopora bleach", "Montipora bleach", "Acropora bleach")] = "Rock and Rubble"
  Moorea_benthos_ts$rls_broad_gr[Moorea_benthos_ts$Substrate %in% c("Sand")] = "Sand"
  
  Moorea_benthos_ts %>% 
    dplyr::filter(rls_broad_gr %in% c("Branching", "Corymbose", "Encrusting", "Massive", "Laminar", "Tabular")) %>% 
    ggplot(aes(y = proportion/3, x = as.character(Year), fill = rls_broad_gr)) + 
    geom_bar(stat="identity", width = 0.5) +
    facet_wrap(~`Marine Area`) 
  
  ############################
  ### CLEANING TIME SERIES ###
  ############################
  
  ts_benthos_indonesia  <- Hoga_dataset_site
  ts_benthos_seychelles <- Seychelles_dataset_site
  ts_benthos_polynesia  <- Moorea_benthos_ts
  
  # Organize each dataset in a single file
  ts_benthos_polynesia <- ts_benthos_polynesia %>% 
    group_by(Year, `Marine Area`, Transect, rls_broad_gr) %>% 
    summarise(proportion = sum(proportion) * 100)
  
  data_1 <- ts_benthos_polynesia %>% dplyr::select(Year, `Marine Area`, Transect, rls_broad_gr, proportion) %>% cbind(dataset = rep("Moorea", length(ts_benthos_polynesia$Year)))
  data_2 <- ts_benthos_seychelles %>% dplyr::select(Year, Site, Replicate, rls_broad_gr, Cover_corrected) %>% cbind(dataset = rep("Seychelles", length(ts_benthos_seychelles$Year)))
  data_3 <- ts_benthos_indonesia %>% dplyr::select(Year, Site, Replicate, rls_broad_gr, Cover_corrected) %>% cbind(dataset = rep("Hoga", length(ts_benthos_indonesia$Year)))
  
  colnames(data_1) <- c("Year", "Site", "Replicate", "taxon", "cover", "dataset")
  colnames(data_2) <- c("Year", "Site", "Replicate", "taxon", "cover", "dataset")
  colnames(data_3) <- c("Year", "Site", "Replicate", "taxon", "cover", "dataset")
  
  # Unique dataset
  ts_benthos_case_study <- rbind(data_1, data_2, data_3)
  
  ###
  
  rls_cols <- rls_info$col
  names(rls_cols) <- rls_info$taxon
  
  ########################################## edit categories
  ##########################################################
  
  # rock and rubble act more as turf?
  ts_benthos_case_study$taxon_original <- ts_benthos_case_study$taxon
  ts_benthos_case_study$taxon <- ifelse(ts_benthos_case_study$taxon=="Rock" | ts_benthos_case_study$taxon=="coral rubble", "turf algae", as.character(ts_benthos_case_study$taxon))
  
  # macroalgae acts similar to fleshy? 
  ts_benthos_case_study$taxon <- ifelse(ts_benthos_case_study$taxon=="canopy forming macroalgae" | ts_benthos_case_study$taxon=="understory macroalgae", "fleshy algae", as.character(ts_benthos_case_study$taxon))
  
  ts_benthos_case_study$N <- ceiling(ts_benthos_case_study$cover) # % cover = number of "individuals"
  sites <- unique(paste(ts_benthos_case_study$Year, ts_benthos_case_study$Site, ts_benthos_case_study$Replicate, sep = "_"))
  ts_benthos_case_study$SurveyID = paste(ts_benthos_case_study$Year, ts_benthos_case_study$Site, ts_benthos_case_study$Replicate, sep = "_")
  ntot   <-length(unique(sites))
  coords <- expand.grid(c(1:10), c(1:10))  # for plotting
  
  # Matching RLS previous dataset
  ts_benthos_case_study$taxon[ts_benthos_case_study$taxon == "Hard_Coral_Branching"] = "Branching"
  ts_benthos_case_study$taxon[ts_benthos_case_study$taxon == "Hard_Coral_Tabular"] = "Tabular"
  ts_benthos_case_study$taxon[ts_benthos_case_study$taxon == "Hard_Coral_Encrusting"] = "Encrusting"
  ts_benthos_case_study$taxon[ts_benthos_case_study$taxon %in% c("Hard_Coral_Massive", "Hard_Coral_Submassive", "Massive")] = "Hemispherical"
  ts_benthos_case_study$taxon[ts_benthos_case_study$taxon %in% c("Alcyonacea", "soft")] = "ahermatypic coral"
  
  indivs <- NULL
  i = sites[1]
  for(i in sites){
    site <- ts_benthos_case_study[ts_benthos_case_study$SurveyID==i  ,]
    inds <- rep(site$taxon, site$N)
    diam_cm <- 10
    xy <- coords[sample(nrow(coords)),][1:length(inds),] # randomise for plotting
    temp <- data.frame(SurveyID=i, taxon = inds, diam_cm, x=xy$Var1, y=xy$Var2)
    indivs  <- rbind(indivs , temp)
    print(round(which(sites==i)/ntot*100, 2)) 
  }
  
  indivs = indivs %>% drop_na()
  table(indivs$SurveyID)
  indivs$x <- as.numeric(indivs$x)
  indivs$y <- as.numeric(indivs$y)
  
  sizeplot <- ggplot(indivs[indivs$SurveyID %in% sample(sites, 10),])+
    geom_point(aes(x = x, y = y, col = taxon), shape = 15, size = 2) +
    theme_void()+
    theme(aspect.ratio=1, plot.title=element_text(hjust=0.5), 
          strip.text=element_text(size=6),
          legend.key.height=unit(2,"mm"), plot.margin=margin(20,20,20,20))+
    facet_wrap(~SurveyID, nrow=5)
  sizeplot
  
  # Cumulated size
  
  surveys <- length(unique(indivs$SurveyID))
  surveys_dataset <- indivs %>% group_by(SurveyID) %>% group_split()
  
  test <- surveys_dataset[[1]] %>% group_by(taxon) %>% arrange(taxon, x, y)
  check_XY_difference <- function(x, y) { c(FALSE, abs(diff(x)) == 1 | abs(diff(y)) == 1) }
  
  # Sum the diam_cm values based on the condition within the same taxon
  
  surveys_size_arranged <- vector(mode = "list", length = surveys)
  for (i in 1:surveys) {
    surveys_size_arranged[[i]] <- surveys_dataset[[i]] %>%
      group_by(taxon) %>%
      mutate(group_id = cumsum(!check_XY_difference(x, y))) %>%
      group_by(taxon, group_id) %>%
      summarise(SurveyID = first(SurveyID),
                X = mean(x),
                Y = mean(y),
                diam_cm = sum(diam_cm)) }
  
  indivs = surveys_size_arranged %>% bind_rows()
  
  hemisphere<-function(d, tb, ...){
    r <- (d/2)
    SA<-(2*pi*r^2)
    TB<-SA*tb
    VOL<-((2/3)*pi*r^3)
    FULL <- VOL # No branch space
    data.frame(SA, TB, VOL, FULL)}
  
  flat<-function(d, th, tb, ...){
    r <- (d/2)
    SA<-pi*(r^2)
    TB<-SA*tb
    VOL<-SA*th
    FULL <- VOL # No branch space
    data.frame(SA, TB, VOL, FULL)}
  
  branches<-function(d, bh, br, bpa, tb, ...){ 
    r <- (d/2)
    SA<-(pi*(r^2)*(bpa*(2*pi*br*bh)+(pi*(br^2))))
    TB<-SA*tb
    VOL<-(pi*(r^2)*(bpa*(pi*(br^2)*bh)))
    FULL <- (pi*((d/2)^2)) * bh
    data.frame(SA,TB, VOL, FULL)}
  
  cone<-function(d, bh, th, tb, ...){
    r <- d/2
    #SA <- pi*r*(r+sqrt((bh^2)+(r^2)))
    SA <- pi * r * sqrt((r^2)+(bh^2)) #pi*r*l
    TB<-SA*tb
    VOL<-th*(SA) 
    FULL <- (pi*(r^2)) * bh
    data.frame(SA,TB, VOL, FULL)}
  
  ############################################# add geometry
  ##########################################################
  
  # use formula to get VOL, SA, FULL, TB
  # use growth rates to get geometry after 1 year
  
  growth <- subset(traits, trait_name=="Growth rate")
  growth <- aggregate(value~rls_group, growth, mean) 
  growth # means for now. / perhaps randomly sample from dist?  
  indivs$growth <- growth$value[match(indivs$taxon, growth$rls_group)]
  indivs$diam_cm_next <- indivs$diam_cm + (indivs$growth /10 )
  head(indivs, 30)
  
  indivs$taxon <- ifelse(indivs$taxon=="Rock and Rubble" | indivs$taxon=="Rubble", "Rock and Rubble", as.character(indivs$taxon))
  indivs$taxon <- ifelse(indivs$taxon=="Macroalgae", "fleshy algae", as.character(indivs$taxon))
  
  
  # set formula
  indivs$fun <- rls_info$geometry[match(indivs$taxon, rls_info$taxon)]
  use <- indivs[complete.cases(indivs), ] # removes all without trait data 
  head(use, 30)
  
  use <- use %>% data.frame()
  
  t0 <- NULL
  t1 <- NULL
  ntot <- nrow(use)
  
  for (i in 1:ntot){
    funct <- use[i, "fun"]
    gp <- use[i,"taxon"]
    pars <- geo[geo$rls==gp, c("br","bh", "bpa", "th", "tb")]
    pars$d <- use[i,"diam_cm"]
    t0 <- rbind(t0, unlist(do.call(Map, c(f= as.name(funct), pars))))
    pars$d <- use[i,"diam_cm_next"]
    t1 <- rbind(t1, unlist(do.call(Map, c(f= as.name(funct), pars))))
    print(round(which(c(1:ntot)==i)/ntot*100, 2))
  }
  colnames(t1) <- paste(colnames(t1), "_next", sep="")
  morphs <- cbind(use, t0, t1)
  head(morphs)
  
  # Rugosity corrected from Husband et al. 2021 (https://doi.org/10.1007/s00338-022-02247-6)
  morphs$Rugosity[morphs$taxon %in% c("Halimeda", "fleshy algae", "turf algae")] <- 1
  morphs$Rugosity[morphs$taxon %in% c("Encrusting", "coralline algae")] <- 1.6
  morphs$Rugosity[morphs$taxon == "Hemispherical"] <- 1.8
  morphs$Rugosity[morphs$taxon == "Tabular"] <- 2.0
  morphs$Rugosity[morphs$taxon == "Corymbose"] <- 2.3
  morphs$Rugosity[morphs$taxon == "Branching"] <- 3.0
  morphs$Rugosity[morphs$taxon == "Laminar"] <- 3.5
  
  # Rugosity per SurveyID
  morphs_Rugo <- morphs %>% mutate(Rugosity_reef = (Rugosity * diam_cm)) %>% group_by(SurveyID) %>% 
    summarise(Rugosity_reef = mean(Rugosity_reef)/10)
  
  hist(morphs_Rugo$Rugosity_reef)
  
  # Merge Rugosity to previous data
  morphs <- merge(morphs, morphs_Rugo) %>% dplyr::select(-Rugosity) %>% rename(Rugosity = Rugosity_reef)
  
  # merge additional traits with morphs data 
  trait_avs <- aggregate(value~rls_group+trait_name, traits, mean) # or sample?
  trait_avs
  
  
  ######################## (1) carbonate production & storage
  
  # change in volume and skeletal density
  sd_dat <- subset(trait_avs, trait_name=="Skeletal density")
  morphs$skel_den <- sd_dat$value[match(morphs$taxon, sd_dat$rls_group)]
  morphs$Storage_kg <- (morphs$VOL*morphs$skel_den )/1000 # (cm3 * gcm3 = g )
  morphs$Accretion_kgyr <- ((morphs$VOL_next-morphs$VOL)*morphs$skel_den)/1000 
  morphs$Storage_kg[is.na(morphs$Storage_kg)] <- 0  # others contribute 0
  morphs$Accretion_kgyr[is.na(morphs$Accretion_kgyr)] <- 0 # others contribute 0
  
  ######################## (2) branch space & rugosity
  
  morphs$BranchSpace_cm3<- (morphs$FULL - morphs$VOL) 
  remove <- c("fleshy algae", "turf algae", "Halimeda") # don't contribute?
  morphs$BranchSpace_cm3[morphs$taxon %in% remove ]<-NA
  morphs$BranchSpace_cm3[is.na(morphs$BranchSpace_cm3)] <- 0 # others contribute 0
  
  ######################## (3) photosynthesis & calcification
  
  gpp_dat <- subset(trait_avs, trait_name=="Gross photosynthesis")
  morphs$gpp <- gpp_dat$value[match(morphs$taxon, gpp_dat$rls_group)] #ugcm2hr
  morphs$GPP_ghr <- ( morphs$SA * morphs$gpp) / (1000^2)
  
  calc_dat <- subset(trait_avs, trait_name=="Calcification rate")
  morphs$calc <- calc_dat$value[match(morphs$taxon, calc_dat$rls_group)] #ugcm2hr
  morphs$Calc_ghr <- ( morphs$SA * morphs$calc) / (1000^2)
  morphs$Calc_ghr[is.na(morphs$Calc_ghr)] <- 0 # others contribute 0
  head(morphs)
  
  ######################## (4) organic production & storage
  
  # corals by SA algae by vol
  morphs$calc_status <- rls_info$calcification[match(morphs$taxon, rls_info$taxon)]
  morphs$OrgGrowth_gyr <- (morphs$TB_next - morphs$TB) / 1000 # corals
  morphs$OrgGrowth_gyr <- ifelse(morphs$calc_status=="non_calcifying", 
                                 (((morphs$VOL_next - morphs$VOL)*0.07) - 15) / 1000, morphs$OrgGrowth_gyr)
  # based on https://doi.org/10.1515/BOT.2002.063
  # check halimeda? 
  
  morphs$BiomassStand_g <- (morphs$TB) / 1000
  morphs$BiomassStand_g  <- ifelse(morphs$calc_status=="non_calcifying", 
                                   ((morphs$VOL*0.07) - 15) / 1000, morphs$BiomassStand_g) 
  
  # based on https://doi.org/10.1515/BOT.2002.063
  
  pred <- predictors_ts %>% 
    dplyr::select(Year, Location, mean_sst) %>% 
    mutate(merge = paste(Location, Year, sep = "_"))
  
  dataset_sites_merging <- data.frame(Site = ts_benthos_case_study$Site, Location = ts_benthos_case_study$dataset) %>% distinct()
  morphs$Site = gsub("^\\d+_(.*)_\\d+$", "\\1", morphs$SurveyID)
  morphs$Year <- as.numeric(substr(morphs$SurveyID, 1, 4))
  
  morphs <- morphs %>% left_join(dataset_sites_merging) #%>% left_join(pred)
  morphs$merge <- paste(morphs$Location, morphs$Year, sep = "_")
  morphs <- morphs %>% left_join(pred, by = "merge") %>% dplyr::select(-merge)
  
  algae_temp_factor <- function(toptprime, erprime, ei, temp_c) {
    topt <- 288.15 + plogis(toptprime) * 30
    er <- ei * plogis(erprime)
    boltz <- 8.62e-5
    ei * exp(
      er * ((temp_c + 273.15) - topt) / ((temp_c + 273.15) * boltz * topt)
    ) /
      (ei - er * (1 - exp(
        ei * ((temp_c + 273.15) - topt) / ((temp_c + 273.15) * boltz * topt))
      ))
  }
  coral_temp_factor <- function(g_max = 1, rho, temp_shift_c, ln_sigma, temp_c) {
    g_max * exp( -1 * exp(rho * (temp_c - temp_shift_c)) -
                   exp(ln_sigma) * (temp_c - temp_shift_c)^2
    )
  }
  
  str(morphs)
  
  morphs$coral_sst_factor <- coral_temp_factor( g_max = 1, rho = 1.022, temp_shift_c = 30.197, ln_sigma = -4.518,    temp_c = morphs$mean_sst)
  
  morphs$algae_sst_factor <- algae_temp_factor(
    toptprime = 0.21, erprime = -1.31, ei = 3, temp_c = morphs$mean_sst  )
  
  algae <- c("coralline algae", "fleshy algae", "Halimeda", "turf algae")
  morphs$coral_sst_factor[morphs$taxon %in% algae] <- 1
  morphs$algae_sst_factor[!morphs$taxon %in% algae] <- 1
  rates <- c("OrgGrowth_gyr", "GPP_ghr", "Accretion_kgyr", "Calc_ghr")
  
  morphs[,rates] <- morphs[,rates] * morphs$coral_sst_factor * morphs$algae_sst_factor
  
  ###################################### taxa contributions
  ##########################################################
  # first sum functions by site + taxon
  
  functions <- c("Storage_kg", "Accretion_kgyr", "GPP_ghr", "Calc_ghr", "BranchSpace_cm3", "OrgGrowth_gyr", "BiomassStand_g")
  
  fun_taxa <- aggregate(.~SurveyID + taxon, morphs[,c("SurveyID", "taxon", functions)], sum)
  
  plotdat <- reshape2::melt(fun_taxa, id=c("SurveyID", "taxon"))
  plotdat <- plotdat[plotdat$SurveyID %in% unique(plotdat$SurveyID)[100:170],]
  taxa_plot <- ggplot(data=plotdat)+
    geom_bar(aes(y=as.factor(SurveyID), x=value, fill=taxon), stat="identity")+
    scale_fill_manual(values=rls_cols)+
    facet_wrap(~variable, scales="free_x", nrow=1)+
    theme_bw()+ylab("SurveyID")+
    theme(axis.text.y=element_text(size=2), strip.background=element_blank())
  taxa_plot
  
  funs <- merge(fun_taxa, morphs_Rugo, by = "SurveyID") %>% rename(Rugosity = Rugosity_reef)
  
  label_decomposition <- str_split(funs$SurveyID, fixed("_"))
  label_decomposition_Location = NA ; label_decomposition_Year = NA
  for (i in 1:length(label_decomposition)) {
    label_decomposition_Location[[i]] <- label_decomposition[[i]][2]
    label_decomposition_Year[[i]] <- label_decomposition[[i]][1]}
  funs$Location <- abind::abind(label_decomposition_Location)
  funs$Year <- abind::abind(label_decomposition_Year)
  
  funs = funs %>% rename(Site = Location) %>% left_join(dataset_sites_merging)
  
  ### Arrange Dataset for Figure_4 Benthos
  Moorea_function <- funs %>% dplyr::filter(Location %in% c("Moorea"), Year %in% c(2005, 2010, 2015)) %>% 
    mutate(Phase = c(rep("before", 202), rep("during", 56), rep("post", 168)))
  Hoga_function <- funs %>% dplyr::filter(Location %in% c("Hoga"), Year %in% c(2005, 2010, 2015)) %>% 
    mutate(Phase = c(rep("before", 72), rep("during", 64), rep("post", 73)))
  Seychelles_function <- funs %>% dplyr::filter(Location %in% c("Seychelles")) %>% 
    mutate(Phase = c(rep("before", 848), rep("during", 862), rep("post", 1043)))
  # Standardisation
  functions_ref <- rbind(Moorea_function %>% dplyr::filter(Phase == "before"),
                         Hoga_function %>% dplyr::filter(Phase == "before"),
                         Seychelles_function %>% dplyr::filter(Phase == "before")) %>% 
    dplyr::select(c("Year", "Location", "Site", "Storage_kg", "Accretion_kgyr", "GPP_ghr", "Calc_ghr", "BranchSpace_cm3", "OrgGrowth_gyr", "Rugosity")) %>% 
    group_by(Year, Location, Site) %>% summarise_all(mean)
  functions_ref <- data.frame(Location         = rep(functions_ref$Location, 7),
                              Site             = rep(functions_ref$Site, 7),
                              Year             = rep(functions_ref$Year, 7),
                              Functions_ref    = c(functions_ref$Storage_kg, functions_ref$Accretion_kgyr,
                                                   functions_ref$GPP_ghr, functions_ref$Calc_ghr,
                                                   functions_ref$BranchSpace_cm3, functions_ref$OrgGrowth_gyr,
                                                   functions_ref$Rugosity),
                              Functions_labels = c(rep("(11) STOR", length(functions_ref$Location)),
                                                   rep("(12) INORG. G", length(functions_ref$Location)), 
                                                   rep("(8) GPP", length(functions_ref$Location)), 
                                                   rep("(10) CALC", length(functions_ref$Location)), 
                                                   rep("(14) SPACE", length(functions_ref$Location)),
                                                   rep("(9) ORG.G", length(functions_ref$Location)), 
                                                   rep("(13) RUG", length(functions_ref$Location))))
  
  functions_benthos <- rbind(Moorea_function, Hoga_function, Seychelles_function) %>% data.frame()
  functions_benthos_df <- data.frame(Location         = rep(functions_benthos$Location, 7),
                                     Site             = rep(functions_benthos$Site, 7),
                                     Phase            = rep(functions_benthos$Phase, 7),
                                     Functions_values = c(functions_benthos$Storage_kg, functions_benthos$Accretion_kgyr,
                                                          functions_benthos$GPP_ghr, functions_benthos$Calc_ghr,
                                                          functions_benthos$BranchSpace_cm3, functions_benthos$OrgGrowth_gyr,
                                                          functions_benthos$Rugosity),
                                     Functions_labels = c(rep("(11) STOR", length(functions_benthos$Location)),
                                                          rep("(12) INORG. G", length(functions_benthos$Location)), 
                                                          rep("(8) GPP", length(functions_benthos$Location)), 
                                                          rep("(10) CALC", length(functions_benthos$Location)), 
                                                          rep("(14) SPACE", length(functions_benthos$Location)),
                                                          rep("(9) ORG.G", length(functions_benthos$Location)), 
                                                          rep("(13) RUG", length(functions_benthos$Location))))
  
  
  functions_benthos_df = functions_benthos_df %>% right_join(functions_ref, by = c("Location", "Site", "Functions_labels")) %>% 
    mutate(., standarized_function = Functions_values / Functions_ref)
  
  functions_benthos_df_plot <- functions_benthos_df %>% group_by(Location, Site, Phase, Functions_labels) %>% 
    summarise(standarized_function_avg = mean(standarized_function), standarized_function_sd = sd(standarized_function)) 
  # %>% dplyr::filter(Location != "Cousin Patch", Location != "Mahe NW Carbonate")
  
  functions_benthos_df_plot$Functions_labels <- factor(functions_benthos_df_plot$Functions_labels, 
                                                       levels = c("(8) GPP", "(9) ORG.G", "(10) CALC", "(11) STOR",
                                                                  "(12) INORG. G", "(13) RUG", "(14) SPACE"))
  
  # Final Figure
  Figure_4_benthos <- ggplot(functions_benthos_df_plot) +
    geom_hline(yintercept = 1, linetype = 3)  +
    geom_point(aes(x = Phase, y = standarized_function_avg, color = Functions_labels),
               position = position_dodge(width = 0.5), size = 3) +
    geom_line(aes(x = as.numeric(as.factor(Phase)), y = standarized_function_avg, color = Functions_labels,
                  linetype = Functions_labels),
              position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = standarized_function_avg - standarized_function_sd, 
                      ymax = standarized_function_avg + standarized_function_sd, 
                      x = Phase, color = Functions_labels),
                  alpha = 0.6,  position = position_dodge(width = 0.5), width = 0, linewidth = 2) +
    labs(x = "Phase", y = "Scaled function", color = "Function", linetype = "Function") +
    scale_color_manual(values = c("#BF961BFF", "#997100FF", "#BE7F36FF", "#F19564FF", "#FA885BFF", "#E8673BFF", "#D94B2BFF")) +
    scale_y_continuous(limits = c(-0.1,3)) + 
    facet_wrap(Location~Site) +
    theme_bw() +
    theme(strip.text = element_text(hjust = 0, size = 12),
          strip.background = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          plot.margin = unit(c(0.5,0.5,1,0.5), "cm"),
          plot.title = element_text(size = 10, face = "bold"),
          axis.line = element_line(colour = "black"))
  
  xlsx::write.xlsx(functions_benthos_df_plot %>% data.frame(), "output/functions_benthos_df_ts.xlsx")
  
  ### Add fish data from Nina
  Fish_Moorea     <- Fish_functions %>% dplyr::filter(dataset_id == "Moorea_MPA", zone == "forereef", year %in% c(2005, 2010, 2015))
  Fish_Seychelles <- Fish_functions %>% dplyr::filter(dataset_id == "SEYC", year %in% c(1994, 2005, 2008))
  Fish_Hoga       <- Fish_functions %>% dplyr::filter(dataset_id == "HOGA", zone == "Slope", year %in% c(2005, 2010, 2015))
  
  functions_benthos_final = functions_benthos %>% 
    dplyr::select("Year", "Site", "Location", "Site", "Storage_kg", "Accretion_kgyr", "GPP_ghr", "Calc_ghr", "BranchSpace_cm3", "OrgGrowth_gyr", "Rugosity") %>% 
    group_by(Year, Location, Site) %>%
    summarise_all(list(mean = mean, sd = sd)) 
  
  table(functions_benthos_final$Year, functions_benthos_final$Location)
  
  Fish_functions  <- rbind(Fish_Moorea, Fish_Seychelles, Fish_Hoga) %>% 
    dplyr::select("year", "site", "location", "biom", "herb", "plank", "pisc", "prod", "exP", "egP", "exN", "egN", "Ic", 
                  "turn", "egNP", "exNP", "outN", "outP", "outNP", "coral_cover") %>% 
    group_by(year, location, site) %>%
    summarise_all(list(mean = mean, sd = sd)) 
  
  unique(Fish_functions$site)
  unique(functions_benthos$Site)
  
  functions_benthos$Site <- gsub(" ", "", functions_benthos$Site)
  Fish_functions$site <- gsub(" ", "", Fish_functions$site)
  
  head(functions_benthos)
  
  benthos_final = functions_benthos %>% dplyr::select(-c(SurveyID, taxon)) %>% group_by(Location, Site, Year, Phase) %>% 
    summarise(Storage_avg = mean(Storage_kg), Storage_sd = sd(Storage_kg), 
              Accretion_kgyr_avg = mean(Accretion_kgyr), Accretion_kgyr_sd = sd(Accretion_kgyr), 
              GPP_ghr_avg = mean(GPP_ghr), GPP_ghr_sd = sd(GPP_ghr),
              Calc_ghr_avg = mean(Calc_ghr), Calc_ghr_sd = sd(Calc_ghr),
              BranchSpace_cm3_avg = mean(BranchSpace_cm3), BranchSpace_cm3_sd = sd(BranchSpace_cm3),
              OrgGrowth_gyr_avg = mean(OrgGrowth_gyr), OrgGrowth_gyr_sd = sd(OrgGrowth_gyr),
              BiomassStand_g_avg = mean(BiomassStand_g), BiomassStand_g_sd = sd(BiomassStand_g),
              Rugosity_avg = mean(Rugosity), Rugosity_sd = sd(Rugosity)) 
  
  Fish_functions$location[Fish_functions$location %in% c("Cousin", "Mahe", "Praslin", "Ste")] = "Seychelles"
  
  Fish_functions = Fish_functions %>% rename(Year = year, Location = location, Site = site)
  final_dataset = merge(benthos_final, Fish_functions, by = c("Year", "Location", "Site"), all = T) 
  
  ################
  ### FIGURE 4 ###
  ################
  
  # Adding case study as the same format of RLS dataset from Mike
  case_study <- case_study %>% 
    select(c(2:5, 7, 9, 11, 13, 15, 19, 22:26, 28, 31)) %>% drop_na() %>% 
    mutate(SurveyID = paste(Location, Site, Year, sep = "_")) %>% 
    select(-c(Year, Location, Site)) %>% column_to_rownames("SurveyID") %>% 
    select(8:9, 11, 10, 13, 12, 14, 3, 6, 4, 1, 2, 7, 5)
  colnames(case_study) <- c("herb", "plank", "prod","pisc",  "exN", "exP", "turn", 
                            "GPP_ghr", "OrgGrowth_gyr",  "Calc_ghr", "Storage_kg", 
                            "Accretion_kgyr", "Rugosity", "BranchSpace_cm3")
  
  head(fish)
  
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
  
  benth <- Benthos_Functions
  fish <- Fish_Functions
  
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
  
  colnames(fish)
  colnames(benth)
  
  all[, f.fun] <- fish[match(all$SurveyID, fish$SurveyID), f.fun]
  all[, b.fun] <- benth[match(all$SurveyID, fish$SurveyID), b.fun]
  
  all <- all[complete.cases(all),]
  nrow(all)
  
  rownames(all)<-all$SurveyID
  
  log.all <- all[,c(f.fun, b.fun)]+1
  case_study <- case_study[,c(f.fun, b.fun)]+1
  log.all <- rbind(log.all, case_study)
  #log.all[log.all==0] <- NA
  #log.all <- log.all[complete.cases(log.all),]
  log.all <- log(log.all)
  nrow(log.all)
  
  pca <- prcomp(log.all, center=T, scale=T)
  
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
  p1
  
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
  p2
  
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
  p3
  
  plot_grid(p1, p2, p3, nrow=1)
  
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
  p4
  
  
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
  p5
  
  Figure_2 <- plot_grid(
    plot_grid(p1, p2, p3, nrow=1, labels=c("a", "b", "c")), 
    plot_grid(p4, p5, rel_widths=c(1, 0.8), labels=c("d","e")), 
    ncol=1)
  
  pca_df <- pc.df %>% rownames_to_column("Label")
  pca_df_sites <- pca_df[c(1172:1396),] %>% 
    dplyr::select(-c(SurveyID, biom, cover))
  
  label_decomposition <- str_split(pca_df_sites$Label, fixed("_"))
  for (i in 1:length(label_decomposition)) {
    pca_df_sites$Location[i] <- label_decomposition[[i]][1]
    pca_df_sites$Site[i] <- label_decomposition[[i]][2]
    pca_df_sites$Year[i] <- as.numeric(label_decomposition[[i]][3])}
  
  pca_df_sites %>% group_by(Location) %>% summarise(min = min(Year), max = max(Year)) %>% mutate(diff = max - min)
  pca_df_sites %>% distinct(Location, Site) %>% group_by(Location) %>% summarise(n = n())
  
  (CS_1A <- pca_df_sites %>% dplyr::filter(Location == "Hoga") %>% 
      ggplot(aes(x = Year, y = PC1)) + 
      geom_line(aes(col = Site, group = Site), show.legend = F) + 
      geom_point(shape = 21, aes(fill = Site), color = "black", size = 0.75, show.legend = F) +
      scale_x_continuous(name = "", breaks = seq(1995, 2020, 5), limits = c(1994, 2020)) +
      scale_y_continuous(name = "PC1 (33.3%)", breaks = seq(-1.5, 4.5, 2), limits = c(-1.75, 4.5)) +
      theme_classic() +
      scale_color_manual(values = c("#8B6A06", "#B78C08","#E3AE09", "#FFD700", "#FFE55C", "#FFED8A")) +
      scale_fill_manual( values = c("#8B6A06", "#B78C08","#E3AE09", "#FFD700", "#FFE55C", "#FFED8A")) +
      ggtitle("Indonesia") +
      theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
            panel.grid = element_line(colour = NA),
            axis.text.x = element_text(size = 0),
            axis.text.y = element_text(size = 12),
            axis.title.x = element_text(size = 14, vjust = -3),
            axis.title.y = element_text(size = 14, vjust = 3),
            legend.title = element_text(size = 15),
            strip.text = element_text(size = 14),
            plot.title = element_text(size = 15, face = "bold"),
            plot.margin = unit(c(.5, .5, .5, .5), "cm")))
  
  (CS_1B <- pca_df_sites %>% dplyr::filter(Location == "Hoga") %>% 
      ggplot(aes(x = Year, y = PC2)) + 
      geom_line(aes(col = Site, group = Site), show.legend = F) + 
      geom_point(shape = 21, aes(fill = Site), color = "black", size = 0.75, show.legend = F) +
      scale_x_continuous(name = "", breaks = seq(1995, 2020, 5), limits = c(1994, 2020)) +
      scale_y_continuous(name = "PC2 (24.5%)", breaks = seq(-4, 4, 2), limits = c(-4, 4)) +
      theme_classic() +
      scale_color_manual(values = c("#8B6A06", "#B78C08","#E3AE09", "#FFD700", "#FFE55C", "#FFED8A")) +
      scale_fill_manual( values = c("#8B6A06", "#B78C08","#E3AE09", "#FFD700", "#FFE55C", "#FFED8A")) +
      theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
            panel.grid = element_line(colour = NA),
            axis.text.x = element_text(size = 0),
            axis.text.y = element_text(size = 12),
            axis.title.x = element_text(size = 14, vjust = -3),
            axis.title.y = element_text(size = 14, vjust = 3),
            legend.title = element_text(size = 15),
            strip.text = element_text(size = 14),
            plot.margin = unit(c(.5, .5, .5, .5), "cm")))
  
  (CS_1C <- pca_df_sites %>% dplyr::filter(Location == "Hoga") %>% 
      ggplot(aes(x = Year, y = PC3)) + 
      geom_line(aes(col = Site, group = Site), show.legend = F) + 
      geom_point(shape = 21, aes(fill = Site), color = "black", size = 0.75, show.legend = F) +
      scale_x_continuous(name = "", breaks = seq(1995, 2020, 5), limits = c(1994, 2020)) +
      scale_y_continuous(name = "PC3 (14.4%)", breaks = seq(-3.5, 2.5, 2), limits = c(-3.5, 2.5)) +
      theme_classic() +
      scale_color_manual(values = c("#8B6A06", "#B78C08","#E3AE09", "#FFD700", "#FFE55C", "#FFED8A")) +
      scale_fill_manual( values = c("#8B6A06", "#B78C08","#E3AE09", "#FFD700", "#FFE55C", "#FFED8A")) +
      theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
            panel.grid = element_line(colour = NA),
            axis.text.x = element_text(size = 0),
            axis.text.y = element_text(size = 12),
            axis.title.x = element_text(size = 14, vjust = -3),
            axis.title.y = element_text(size = 14, vjust = 3),
            legend.title = element_text(size = 15),
            strip.text = element_text(size = 14),
            plot.margin = unit(c(.5, .5, .5, .5), "cm")))
  
  (CS_1D <- pca_df_sites %>% dplyr::filter(Location == "Hoga") %>% 
      ggplot(aes(x = Year, y = PC4)) + 
      geom_line(aes(col = Site, group = Site), show.legend = F) + 
      geom_point(shape = 21, aes(fill = Site), color = "black", size = 0.75, show.legend = F) +
      scale_x_continuous(name = "", breaks = seq(1995, 2020, 5), limits = c(1994, 2020), labels = c("", 2000, "", 2010, "", 2020)) +
      scale_y_continuous(name = "PC4 (6.7%)", breaks = seq(-2, 2, 2), limits = c(-2, 2.25)) +
      theme_classic() +
      scale_color_manual(values = c("#8B6A06", "#B78C08","#E3AE09", "#FFD700", "#FFE55C", "#FFED8A")) +
      scale_fill_manual( values = c("#8B6A06", "#B78C08","#E3AE09", "#FFD700", "#FFE55C", "#FFED8A")) +
      theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
            panel.grid = element_line(colour = NA),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            axis.title.x = element_text(size = 14, vjust = -3),
            axis.title.y = element_text(size = 14, vjust = 3),
            legend.title = element_text(size = 15),
            strip.text = element_text(size = 14),
            plot.margin = unit(c(.5, .5, .5, .5), "cm")))
  
  (CS_2A <- pca_df_sites %>% dplyr::filter(Location == "Moorea") %>% 
      ggplot(aes(x = Year, y = PC1)) + 
      geom_line(aes(col = Site, group = Site), show.legend = F) + 
      geom_point(shape = 21, aes(fill = Site), color = "black", size = 0.75, show.legend = F) +
      scale_x_continuous(name = "", breaks = seq(1995, 2020, 5), limits = c(1994, 2020)) +
      scale_y_continuous(name = "", breaks = seq(-1.5, 4.5, 2), limits = c(-1.75, 4.5)) +
      theme_classic() +
      scale_color_manual(values = c("#002600", "#003300", "#004c00", "#006600", "#008000", "#198c19", 
                                    "#329932", "#66b266", "#7fbf7f", "#99cc99", "#b2d8b2")) +
      scale_fill_manual( values = c("#002600", "#003300", "#004c00", "#006600", "#008000", "#198c19", 
                                    "#329932", "#66b266", "#7fbf7f", "#99cc99", "#b2d8b2")) +
      ggtitle("French Polynesia") +
      theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
            panel.grid = element_line(colour = NA),
            axis.text.x = element_text(size = 0),
            axis.text.y = element_text(size = 0),
            axis.title.x = element_text(size = 14, vjust = -3),
            axis.title.y = element_text(size = 14, vjust = 3),
            legend.title = element_text(size = 15),
            strip.text = element_text(size = 14),
            plot.title = element_text(size = 15, face = "bold"),
            plot.margin = unit(c(.5, .5, .5, .5), "cm")))
  
  (CS_2B <- pca_df_sites %>% dplyr::filter(Location == "Moorea") %>% 
      ggplot(aes(x = Year, y = PC2)) + 
      geom_line(aes(col = Site, group = Site), show.legend = F) + 
      geom_point(shape = 21, aes(fill = Site), color = "black", size = 0.75, show.legend = F) +
      scale_x_continuous(name = "", breaks = seq(1995, 2020, 5), limits = c(1994, 2020)) +
      scale_y_continuous(name = "", breaks = seq(-4, 4, 2), limits = c(-4, 4)) +
      theme_classic() +
      scale_color_manual(values = c("#002600", "#003300", "#004c00", "#006600", "#008000", "#198c19", 
                                    "#329932", "#66b266", "#7fbf7f", "#99cc99", "#b2d8b2")) +
      scale_fill_manual( values = c("#002600", "#003300", "#004c00", "#006600", "#008000", "#198c19", 
                                    "#329932", "#66b266", "#7fbf7f", "#99cc99", "#b2d8b2")) +
      theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
            panel.grid = element_line(colour = NA),
            axis.text.x = element_text(size = 0),
            axis.text.y = element_text(size = 0),
            axis.title.x = element_text(size = 14, vjust = -3),
            axis.title.y = element_text(size = 14, vjust = 3),
            legend.title = element_text(size = 15),
            strip.text = element_text(size = 14),
            plot.margin = unit(c(.5, .5, .5, .5), "cm")))
  
  (CS_2C <- pca_df_sites %>% dplyr::filter(Location == "Moorea") %>% 
      ggplot(aes(x = Year, y = PC3)) + 
      geom_line(aes(col = Site, group = Site), show.legend = F) + 
      geom_point(shape = 21, aes(fill = Site), color = "black", size = 0.75, show.legend = F) +
      scale_x_continuous(name = "", breaks = seq(1995, 2020, 5), limits = c(1994, 2020)) +
      scale_y_continuous(name = "", breaks = seq(-3.5, 2.5, 2), limits = c(-3.5, 2.5)) +
      theme_classic() +
      scale_color_manual(values = c("#002600", "#003300", "#004c00", "#006600", "#008000", "#198c19", 
                                    "#329932", "#66b266", "#7fbf7f", "#99cc99", "#b2d8b2")) +
      scale_fill_manual( values = c("#002600", "#003300", "#004c00", "#006600", "#008000", "#198c19", 
                                    "#329932", "#66b266", "#7fbf7f", "#99cc99", "#b2d8b2")) +
      theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
            panel.grid = element_line(colour = NA),
            axis.text.x = element_text(size = 0),
            axis.text.y = element_text(size = 0),
            axis.title.x = element_text(size = 14, vjust = -3),
            axis.title.y = element_text(size = 14, vjust = 3),
            legend.title = element_text(size = 15),
            strip.text = element_text(size = 14),
            plot.margin = unit(c(.5, .5, .5, .5), "cm")))
  
  (CS_2D <- pca_df_sites %>% dplyr::filter(Location == "Moorea") %>% 
      ggplot(aes(x = Year, y = PC4)) + 
      geom_line(aes(col = Site, group = Site), show.legend = F) + 
      geom_point(shape = 21, aes(fill = Site), color = "black", size = 0.75, show.legend = F) +
      scale_x_continuous(name = "", breaks = seq(1995, 2020, 5), limits = c(1994, 2020), labels = c("", 2000, "", 2010, "", 2020)) +
      scale_y_continuous(name = "", breaks = seq(-2, 2, 2), limits = c(-2, 2.25)) +
      theme_classic() +
      scale_color_manual(values = c("#002600", "#003300", "#004c00", "#006600", "#008000", "#198c19", 
                                    "#329932", "#66b266", "#7fbf7f", "#99cc99", "#b2d8b2")) +
      scale_fill_manual( values = c("#002600", "#003300", "#004c00", "#006600", "#008000", "#198c19", 
                                    "#329932", "#66b266", "#7fbf7f", "#99cc99", "#b2d8b2")) +
      theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
            panel.grid = element_line(colour = NA),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 0),
            axis.title.x = element_text(size = 14, vjust = -3),
            axis.title.y = element_text(size = 14, vjust = 3),
            legend.title = element_text(size = 15),
            strip.text = element_text(size = 14),
            plot.margin = unit(c(.5, .5, .5, .5), "cm")))
  
  (CS_3A <- pca_df_sites %>% dplyr::filter(Location == "Seychelles") %>% 
      ggplot(aes(x = Year, y = PC1)) + 
      geom_line(aes(col = Site, group = Site), show.legend = F) + 
      geom_point(shape = 21, aes(fill = Site), color = "black", size = 0.75, show.legend = F) +
      scale_x_continuous(name = "", breaks = seq(1995, 2020, 5), limits = c(1994, 2020)) +
      scale_y_continuous(name = "", breaks = seq(-1.5, 4.5, 2), limits = c(-1.75, 4.5)) +
      ggtitle("Seychelles") +
      theme_classic() +
      scale_color_manual(values = c("#260818", "#2f1969", "#351c75", "#3a2868", "#3f1e5f", "#4a236f", "#54287f", "#5f2d8f",
                                    "#554a75", "#71639c", "#7f6faf", "#8e7cc3", "#c3b8de", "#d9d2e9", "#ead2de", "#e1c0d0",
                                    "#d9aec3", "#d5a6bd", "#ce95b3", "#c888a9", "#c27ba0")) +
      scale_fill_manual( values = c("#260818", "#2f1969", "#351c75", "#3a2868", "#3f1e5f", "#4a236f", "#54287f", "#5f2d8f",
                                    "#554a75", "#71639c", "#7f6faf", "#8e7cc3", "#c3b8de", "#d9d2e9", "#ead2de", "#e1c0d0",
                                    "#d9aec3", "#d5a6bd", "#ce95b3", "#c888a9", "#c27ba0")) +
      theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
            panel.grid = element_line(colour = NA),
            axis.text.x = element_text(size = 0),
            axis.text.y = element_text(size = 0),
            axis.title.x = element_text(size = 14, vjust = -3),
            axis.title.y = element_text(size = 14, vjust = 3),
            legend.title = element_text(size = 15),
            strip.text = element_text(size = 14),
            plot.title = element_text(size = 15, face = "bold"),
            plot.margin = unit(c(.5, .5, .5, .5), "cm")))
  
  (CS_3B <- pca_df_sites %>% dplyr::filter(Location == "Seychelles") %>% 
      ggplot(aes(x = Year, y = PC2)) + 
      geom_line(aes(col = Site, group = Site), show.legend = F) + 
      geom_point(shape = 21, aes(fill = Site), color = "black", size = 0.75, show.legend = F) +
      scale_x_continuous(name = "", breaks = seq(1995, 2020, 5), limits = c(1994, 2020)) +
      scale_y_continuous(name = "", breaks = seq(-4, 4, 2), limits = c(-4, 4)) +
      theme_classic() +
      scale_color_manual(values = c("#260818", "#2f1969", "#351c75", "#3a2868", "#3f1e5f", "#4a236f", "#54287f", "#5f2d8f",
                                    "#554a75", "#71639c", "#7f6faf", "#8e7cc3", "#c3b8de", "#d9d2e9", "#ead2de", "#e1c0d0",
                                    "#d9aec3", "#d5a6bd", "#ce95b3", "#c888a9", "#c27ba0")) +
      scale_fill_manual( values = c("#260818", "#2f1969", "#351c75", "#3a2868", "#3f1e5f", "#4a236f", "#54287f", "#5f2d8f",
                                    "#554a75", "#71639c", "#7f6faf", "#8e7cc3", "#c3b8de", "#d9d2e9", "#ead2de", "#e1c0d0",
                                    "#d9aec3", "#d5a6bd", "#ce95b3", "#c888a9", "#c27ba0")) +
      theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
            panel.grid = element_line(colour = NA),
            axis.text.x = element_text(size = 0),
            axis.text.y = element_text(size = 0),
            axis.title.x = element_text(size = 14, vjust = -3),
            axis.title.y = element_text(size = 14, vjust = 3),
            legend.title = element_text(size = 15),
            strip.text = element_text(size = 14),
            plot.margin = unit(c(.5, .5, .5, .5), "cm")))
  
  (CS_3C <- pca_df_sites %>% dplyr::filter(Location == "Seychelles") %>% 
      ggplot(aes(x = Year, y = PC3)) + 
      geom_line(aes(col = Site, group = Site), show.legend = F) + 
      geom_point(shape = 21, aes(fill = Site), color = "black", size = 0.75, show.legend = F) +
      scale_x_continuous(name = "", breaks = seq(1995, 2020, 5), limits = c(1994, 2020)) +
      scale_y_continuous(name = "", breaks = seq(-3.5, 2.5, 2), limits = c(-3.5, 2.5)) +
      theme_classic() +
      scale_color_manual(values = c("#260818", "#2f1969", "#351c75", "#3a2868", "#3f1e5f", "#4a236f", "#54287f", "#5f2d8f",
                                    "#554a75", "#71639c", "#7f6faf", "#8e7cc3", "#c3b8de", "#d9d2e9", "#ead2de", "#e1c0d0",
                                    "#d9aec3", "#d5a6bd", "#ce95b3", "#c888a9", "#c27ba0")) +
      scale_fill_manual( values = c("#260818", "#2f1969", "#351c75", "#3a2868", "#3f1e5f", "#4a236f", "#54287f", "#5f2d8f",
                                    "#554a75", "#71639c", "#7f6faf", "#8e7cc3", "#c3b8de", "#d9d2e9", "#ead2de", "#e1c0d0",
                                    "#d9aec3", "#d5a6bd", "#ce95b3", "#c888a9", "#c27ba0")) +
      theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
            panel.grid = element_line(colour = NA),
            axis.text.x = element_text(size = 0),
            axis.text.y = element_text(size = 0),
            axis.title.x = element_text(size = 14, vjust = -3),
            axis.title.y = element_text(size = 14, vjust = 3),
            legend.title = element_text(size = 15),
            strip.text = element_text(size = 14),
            plot.margin = unit(c(.5, .5, .5, .5), "cm")))
  
  (CS_3D <- pca_df_sites %>% dplyr::filter(Location == "Seychelles") %>% 
      ggplot(aes(x = Year, y = PC4)) + 
      geom_line(aes(col = Site, group = Site), show.legend = F) + 
      geom_point(shape = 21, aes(fill = Site), color = "black", size = 0.75, show.legend = F) +
      scale_x_continuous(name = "", breaks = seq(1995, 2020, 5), limits = c(1994, 2020), labels = c("", 2000, "", 2010, "", 2020)) +
      scale_y_continuous(name = "", breaks = seq(-2, 2, 2), limits = c(-2, 2.25)) +
      theme_classic() +
      scale_color_manual(values = c("#260818", "#2f1969", "#351c75", "#3a2868", "#3f1e5f", "#4a236f", "#54287f", "#5f2d8f",
                                    "#554a75", "#71639c", "#7f6faf", "#8e7cc3", "#c3b8de", "#d9d2e9", "#ead2de", "#e1c0d0",
                                    "#d9aec3", "#d5a6bd", "#ce95b3", "#c888a9", "#c27ba0")) +
      scale_fill_manual( values = c("#260818", "#2f1969", "#351c75", "#3a2868", "#3f1e5f", "#4a236f", "#54287f", "#5f2d8f",
                                    "#554a75", "#71639c", "#7f6faf", "#8e7cc3", "#c3b8de", "#d9d2e9", "#ead2de", "#e1c0d0",
                                    "#d9aec3", "#d5a6bd", "#ce95b3", "#c888a9", "#c27ba0")) +
      theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
            panel.grid = element_line(colour = NA),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 0),
            axis.title.x = element_text(size = 14, vjust = -3),
            axis.title.y = element_text(size = 14, vjust = 3),
            legend.title = element_text(size = 15),
            strip.text = element_text(size = 14),
            plot.margin = unit(c(.5, .5, .5, .5), "cm")))
  
  library(patchwork)
  (Panel_Figure_CS <- CS_1A + plot_spacer() + CS_2A + plot_spacer() + CS_3A + 
      CS_1B + plot_spacer() + CS_2B + plot_spacer() + CS_3B + 
      CS_1C + plot_spacer() + CS_2C + plot_spacer() + CS_3C + 
      CS_1D + plot_spacer() + CS_2D + plot_spacer() + CS_3D + 
      plot_layout(nrow = 4, widths = c(30, -3.75, 30, -3.75, 30)) )
  
  ## Request Vale
  colors = c("#8B6A06", "#B78C08", "#E3AE09", "#FFD700", "#FFE55C", "#FFED8A", "#002600", "#003300", "#004c00", "#006600", 
             "#008000", "#198c19", "#329932", "#66b266", "#7fbf7f", "#99cc99", "#b2d8b2", "#260818", "#2f1969", "#351c75", 
             "#3a2868", "#3f1e5f", "#4a236f", "#54287f", "#5f2d8f", "#554a75", "#71639c", "#7f6faf", "#8e7cc3", "#c3b8de", 
             "#d9d2e9", "#ead2de", "#e1c0d0", "#d9aec3", "#d5a6bd", "#ce95b3", "#c888a9", "#c27ba0")
  
  pca_ordered = pca_df_sites %>% mutate(Location = factor(Location, levels = c("Hoga", "Moorea", "Seychelles"))) %>% 
    arrange(Location) 
  PC1 = pca_ordered %>% mutate(Site = factor(Site, levels = unique(pca_ordered$Site))) %>% 
    ggplot(aes(x = Year, y = PC1)) + 
    geom_line(aes(col = Site, group = Site), show.legend = F) + 
    geom_point(shape = 21, aes(fill = Site), color = "black", size = 0.75, show.legend = F) +
    scale_x_continuous(name = "", breaks = seq(1995, 2020, 5), limits = c(1994, 2020), labels = c("", 2000, "", 2010, "", 2020)) +
    scale_y_reverse(name = "PC1 (33.3%)", limits = c(4.5, -1.75), breaks = seq(4.5,-1.5, -2)) +
    theme_classic() +
    scale_color_manual(values = colors) + scale_fill_manual( values = colors) +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
          panel.grid = element_line(colour = NA),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14, vjust = -3),
          axis.title.y = element_text(size = 14, vjust = 3),
          legend.title = element_text(size = 15),
          strip.text = element_text(size = 14),
          plot.title = element_text(size = 15, face = "bold"),
          plot.margin = unit(c(.5, .5, .5, .5), "cm"))
  
  PC2 = pca_ordered %>% mutate(Site = factor(Site, levels = unique(pca_ordered$Site))) %>% 
    ggplot(aes(x = Year, y = PC2)) + 
    geom_line(aes(col = Site, group = Site), show.legend = F) + 
    geom_point(shape = 21, aes(fill = Site), color = "black", size = 0.75, show.legend = F) +
    scale_x_continuous(name = "", breaks = seq(1995, 2020, 5), limits = c(1994, 2020), labels = c("", 2000, "", 2010, "", 2020)) +
    scale_y_continuous(name = "PC2 (24.5%)", breaks = seq(-4, 4, 2), limits = c(-4, 4)) +
    theme_classic() +
    scale_color_manual(values = colors) + scale_fill_manual( values = colors) +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
          panel.grid = element_line(colour = NA),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14, vjust = -3),
          axis.title.y = element_text(size = 14, vjust = 3),
          legend.title = element_text(size = 15),
          strip.text = element_text(size = 14),
          plot.margin = unit(c(.5, .5, .5, .5), "cm"))
  
  PC3 <- pca_ordered %>% mutate(Site = factor(Site, levels = unique(pca_ordered$Site))) %>% 
    ggplot(aes(x = Year, y = PC3)) + 
    geom_line(aes(col = Site, group = Site), show.legend = F) + 
    geom_point(shape = 21, aes(fill = Site), color = "black", size = 0.75, show.legend = F) +
    scale_x_continuous(name = "", breaks = seq(1995, 2020, 5), limits = c(1994, 2020), labels = c("", 2000, "", 2010, "", 2020)) +
    scale_y_continuous(name = "PC3 (14.4%)", breaks = seq(-3.5, 2.5, 2), limits = c(-3.5, 2.5)) +
    theme_classic() +
    scale_color_manual(values = colors) + scale_fill_manual( values = colors) +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
          panel.grid = element_line(colour = NA),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14, vjust = -3),
          axis.title.y = element_text(size = 14, vjust = 3),
          legend.title = element_text(size = 15),
          strip.text = element_text(size = 14),
          plot.margin = unit(c(.5, .5, .5, .5), "cm"))
  
  PC4 <- pca_ordered %>% mutate(Site = factor(Site, levels = unique(pca_ordered$Site))) %>% 
    ggplot(aes(x = Year, y = PC4)) + 
    geom_line(aes(col = Site, group = Site), show.legend = F) + 
    geom_point(shape = 21, aes(fill = Site), color = "black", size = 0.75, show.legend = F) +
    scale_x_continuous(name = "", breaks = seq(1995, 2020, 5), limits = c(1994, 2020), labels = c("", 2000, "", 2010, "", 2020)) +
    scale_y_continuous(name = "PC4 (6.7%)", breaks = seq(-2, 2, 2), limits = c(-2, 2.25)) +
    theme_classic() +
    scale_color_manual(values = colors) + scale_fill_manual( values = colors) +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
          panel.grid = element_line(colour = NA),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14, vjust = -3),
          axis.title.y = element_text(size = 14, vjust = 3),
          legend.title = element_text(size = 15),
          strip.text = element_text(size = 14),
          plot.margin = unit(c(.5, .5, .5, .5), "cm"))
  
  (Merged_Figure_CS <- PC1 + plot_spacer() + PC2 + plot_spacer() + PC3 + plot_spacer() + PC4 + 
      plot_layout(nrow = 1, widths = c(30, -1, 30, -1, 30, -1, 30)) )
  
  ggsave("output/fig4_a.png", Merged_Figure_CS, dev = "png", width = 10, height = 3)
  
  ### Correlations
  Moorea <- pca_df_sites %>% dplyr::filter(Location == "Moorea") %>% select(PC1, PC2, PC3, PC4)
  p.mat <- cor_pmat(Moorea)
  corr_Moorea <- round(cor(Moorea), 3)
  cor_3 = ggcorrplot(corr_Moorea, hc.order = F, type = 'lower', show.diag = F, method = "square",
                     colors = c("#E5E4E2", "#008000", "#004c00"), lab = TRUE) + 
    theme(legend.position = "none")
  Indonesia <- pca_df_sites %>% dplyr::filter(Location == "Hoga") %>% select(PC1, PC2, PC3, PC4)
  p.mat <- cor_pmat(Indonesia)
  corr_Indo <- round(cor(Indonesia), 3)
  cor_2 = ggcorrplot(corr_Indo, hc.order = F, type = 'lower', show.diag = F,
                     colors = c("#8F4B07", "#FFED8A", "#FFD700"), lab = TRUE) + 
    theme(legend.position = "none")
  Seychelles <- pca_df_sites %>% dplyr::filter(Location == "Seychelles") %>% 
    select(PC1, PC2, PC3, PC4)
  p.mat <- cor_pmat(Seychelles)
  corr_Sey <- round(cor(Seychelles), 3)
  cor_1 = ggcorrplot(corr_Sey, hc.order = F, type = 'lower', show.diag = F,
                     colors = c("#c27ba0", "#7f6faf", "#351c75"), lab = TRUE) + 
    theme(legend.position = "none")
  (correlations <- cor_1 + cor_2 + cor_3)
  
  ggsave("output/fig4_b.png", correlations , dev = "png", width = 10, height = 3)
  
  
  
}




