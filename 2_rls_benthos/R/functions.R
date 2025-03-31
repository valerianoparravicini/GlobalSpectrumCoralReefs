

### --->> 1 - import and clean rls

clean_rls_dataset <- function(rls_raw, rls_info) {
  
  ### AIM: CALCULATE BENTHIC FUNCTIONS USING SURVEY DATA
  rls_raw$SurveyDate1 <- as.Date(rls_raw$SurveyDate, format="%d/%m/%Y")
  rls_raw <- rls_raw[!rls_raw$SurveyID==912340599,] # 0 % cover for all
  rls_cols <- rls_info$col
  names(rls_cols) <- rls_info$taxon
  
  rls_raw$taxon_original <- rls_raw$taxon
  rls_raw$taxon <- ifelse(rls_raw$taxon=="Rock" | rls_raw$taxon=="coral rubble", "turf algae", as.character(rls_raw$taxon))
  
  # macroalgae acts similar to fleshy? 
  rls_raw$taxon <- ifelse(rls_raw$taxon=="canopy forming macroalgae" | rls_raw$taxon=="understory macroalgae", "fleshy algae", as.character(rls_raw$taxon))
  
  rls <- rls_raw
  
  return(rls)
  
}#clean_rls_dataset



### --->> 2 - simulate sizes // make df with rows = "individuals". give a size in cm // simple model: assume surveys are 1 m2 quadrats. everything 10x10cm.

sim_colony_size_rls <- function(rls) {
  
  rls$N <- round(rls$cover) # % cover = number of "individuals"
  head(rls)
  
  sites <- unique(rls$SurveyID)
  ntot<-length(unique(sites))
  coords <- expand.grid(c(1:10), c(1:10))  # for plotting
  
  
  
  indivs <- parallel::mclapply(sites, function(i) {
    
    set.seed(10)
    site <- rls[rls$SurveyID==i  ,]
    inds <- rep(site$taxon, site$N)
    diam_cm <- 10
    xy <- coords[sample(nrow(coords)),][1:length(inds),] # randomise for plotting
    data.frame(SurveyID=i, taxon = inds, diam_cm, x=xy$Var1, y=xy$Var2)
    
  }, mc.cores=50)
  
  return(do.call(rbind, indivs))
  
  
}#sim_colony_size_rls


### --->> 3 - creates surveys_size_arranged 

surveys_size_arranged <- function(indivs) {
  
  surveys <- length(unique(indivs$SurveyID))
  surveys_dataset <- indivs %>% group_by(SurveyID) %>% group_split()
  
  test <- surveys_dataset[[2]] %>% group_by(taxon) %>% arrange(taxon, x, y)
  check_XY_difference <- function(x, y) { c(FALSE, abs(diff(x)) == 1 | abs(diff(y)) == 1) }
  
  # Sum the diam_cm values based on the condition within the same taxon
  surveys_size_arranged_lp <- pbmcapply::pbmclapply(1:surveys, function(i) {
    
    surveys_dataset[[i]] %>%
      group_by(taxon) %>%
      mutate(group_id = cumsum(!check_XY_difference(x, y))) %>%
      group_by(taxon, group_id) %>%
      summarise(SurveyID = first(SurveyID),
                X = mean(x),
                Y = mean(y),
                diam_cm = sum(diam_cm)) }, mc.cores=50)
  
  indivs = surveys_size_arranged_lp %>% bind_rows()
  
  return(indivs)
  
}#surveys_size_arranged 



### --->> 4 - geometric shape functions (hemisphere, flat, branches, cone)

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



### --->> 5 - add geometry // use formula to get VOL, SA, FULL, TB// use growth rates to get geometry after 1 year

geometry_change_1_year <- function(traits, indivs, geo, rls_info) {
  
  growth <- subset(traits, trait_name=="Growth rate")
  growth <- aggregate(value~rls_group, growth, mean) 
  
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
  
  pb = txtProgressBar(min = 0, max = ntot, initial = 0, style = 3) 
  
  for (i in 1:ntot){
    funct <- use[i, "fun"]
    gp <- use[i,"taxon"]
    pars <- geo[geo$rls==gp, c("br","bh", "bpa", "th", "tb")]
    pars$d <- use[i,"diam_cm"]
    t0 <- rbind(t0, unlist(do.call(Map, c(f= as.name(funct), pars))))
    pars$d <- use[i,"diam_cm_next"]
    t1 <- rbind(t1, unlist(do.call(Map, c(f= as.name(funct), pars))))
    #print(round(which(c(1:ntot)==i)/ntot*100, 2))
    setTxtProgressBar(pb,i)
  }
  
  close(pb)
  
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
  
  return(list(morphs,  trait_avs, morphs_Rugo))
  
}#geometry_change_1_year



### --->> 6 - Temperature factors for algae and coral

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


### --->> 7 - Compute functions

compute_functions <- function(morphs_traits, pred, rls_info) {
  
  morphs <- morphs_traits[[1]]
  trait_avs <- morphs_traits[[2]]
  
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
  
  morphs$coral_sst_factor <- coral_temp_factor( g_max = 1, rho = 1.022, temp_shift_c = 30.197, ln_sigma = -4.518, temp_c = morphs$mean_sst)
  morphs$algae_sst_factor <- algae_temp_factor(toptprime = 0.21, erprime = -1.31, ei = 3, temp_c = morphs$mean_sst  )
  
  algae <- c("coralline algae", "fleshy algae", "Halimeda", "turf algae")
  morphs$coral_sst_factor[morphs$taxon %in% algae] <- 1
  morphs$algae_sst_factor[!morphs$taxon %in% algae] <- 1
  rates <- c("OrgGrowth_gyr", "GPP_ghr", "Accretion_kgyr", "Calc_ghr")
  
  morphs[,rates] <- morphs[,rates] * morphs$coral_sst_factor * morphs$algae_sst_factor
  
  return(morphs)
  
}




taxa_contribution_plot <- function(rls_info, morphs) {
  
  rls_cols <- rls_info$col
  names(rls_cols) <- rls_info$taxon
  
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
  
  taxa_plot 
  
}#taxa_contribution_plot



site_level_functions <- function(rls, rls_sites, rls_info, morphs, morphs_traits) {
  
  morphs_Rugo <- morphs_traits[[3]] 
  
  sites     <- data.frame(SurveyID=unique(rls$SurveyID))
  add_cols  <- c("SiteCode", "SiteLatitude", "SiteLongitude", "Country","SiteName", "Location", "Depth", "SurveyDate")
  sites[,add_cols] <- rls_sites[match(sites$SurveyID, rls_sites$SurveyID), add_cols]
  
  rls$broad_cat <- rls_info$broad_cat[match(rls$taxon, rls_info$taxon)]
  cover         <- aggregate(cover~SurveyID+broad_cat, rls, sum)
  cover         <- dcast(cover, SurveyID~broad_cat, value.var="cover")
  cover$total   <- rowSums(subset(cover, select=-c(SurveyID)))
  cover$flat    <- cover$total-cover$coral
  sites         <- merge(sites, cover, all.x=F, all.y=F)
  
  functions <- c("Storage_kg", "Accretion_kgyr", "GPP_ghr", "Calc_ghr", "BranchSpace_cm3", "OrgGrowth_gyr", "BiomassStand_g")
  funs              <- aggregate(.~ SurveyID, morphs[,c("SurveyID",functions)], sum)
  sites[,functions] <- funs[match(sites$SurveyID, funs$SurveyID), functions]
  
  funs  <- merge(funs, morphs_Rugo, by = "SurveyID")  %>% rename(Rugosity = Rugosity_reef)
  sites <- merge(sites, morphs_Rugo, by = "SurveyID") %>% rename(Rugosity = Rugosity_reef)
  
  sites.all <- na.omit(sites)
  
  #pca       <- prcomp(sqrt(sites.all[,c(functions, "Rugosity")]), scale=T, center=T)
  #sites.all[,c("pc1", "pc2")] <- pca$x[,1:2]
  #vecs      <- data.frame(pc1=pca$rotation[,1], pc2=pca$rotation[,2], 
  #                        labels = rownames(pca$rotation), pc1b=pca$rotation[,1], pc2b=pca$rotation[,2])
  #vecs$pc2b[vecs$labels=="Calc_ghr"] <- vecs$pc2b[vecs$labels=="Calc_ghr"]+0.05
  #vecs$pc2b[vecs$labels=="GPP_ghr"]  <- vecs$pc2b[vecs$labels=="GPP_ghr"]+0.05
  #vars <- round((pca$sdev^2 / sum(pca$sdev^2)), 3)*100
  
  data_corrected_size = melt(sites.all[,c("SurveyID",functions, "Rugosity")], id.var="SurveyID")
  
  write_csv(data_corrected_size, "output/data_corrected_size.csv")
  write_csv(sites, "output/RLS_benthic_functions.csv")
  
  return(data_corrected_size)
}


