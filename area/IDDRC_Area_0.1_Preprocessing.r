#CND14: Drop CPRedo,    #UNCD14: keep all    #CHOPD14: keep all (check if any same image name)
#CHOPD7: keep all   #CN D7: keep all   #UNC D7: keep all
#CHOP D35: keep all   #CN D35: check if near 0 images are removed: otherwise keep all   #UNC R35 Keep all. 
#Keep all D56

#also remove replicates with very few hCOs from UNC 

# INPUT FILES ############################################################################################
rds.here <- "/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/"

D7Area <- readRDS(paste0(rds.here, "IDDRCD7Area.rds"))
D14Area <- readRDS(paste0(rds.here, "IDDRCD14Area.rds"))
D35Area <- readRDS(paste0(rds.here, "IDDRCD35Area.rds"))
D56Area <- readRDS(paste0(rds.here, "IDDRCD56Area.rds"))

#All output files with READY appeneded to original name

#########################################################################################################
#D7
#Remove tails of distribution
Bottom2p5 <- quantile(D7Area$AreaMM2, 0.025)
Top2p5 <- quantile(D7Area$AreaMM2, 0.975)
D7Area <- filter(D7Area, AreaMM2 < Top2p5)
D7Area <- filter(D7Area, AreaMM2 > Bottom2p5)

saveRDS(D7Area, paste0(rds.here, "IDDRCD7AreaReady.rds"))

# D14 ##################################################################################
#remove CP redo, UNCrep7 and ONLY AKs meausrements for CNRep1
D14Area <- filter(D14Area, Initial != "CPReDO")
D14.CNR1 <- filter(D14Area, Experiment == "CNRep1")
D14.CNR1 <- filter(D14.CNR1, Initial == "AK")
D14Area <- filter(D14Area, Experiment != "CNRep1")
D14Area <- filter(D14Area, Experiment != "UNCRep7")

D14Area <- rbind(D14Area, D14.CNR1)

#Remove tails of distribution
Bottom2p5 <- quantile(D14Area$D14AreaMM2, 0.025)
Top2p5 <- quantile(D14Area$D14AreaMM2, 0.975)
D14Area <- filter(D14Area, D14AreaMM2 < Top2p5)
D14Area <- filter(D14Area, D14AreaMM2 > Bottom2p5)

D14Area <- dplyr::rename(D14Area, AreaMM2 = D14AreaMM2)

saveRDS(D14Area, paste0(rds.here, "IDDRCD14AreaReady.rds"))

# D35 ############################################################################
#Remove tails of distribution
Bottom2p5 <- quantile(D35Area$D35AreaMM2, 0.025)
Top2p5 <- quantile(D35Area$D35AreaMM2, 0.975)
D35Area <- filter(D35Area, D35AreaMM2 < Top2p5)
D35Area <- filter(D35Area, D35AreaMM2 > Bottom2p5)

#needs to trim some CHOP replicates
D35Area.CNR1 <- filter(D35Area, Experiment == "CNRep1")
D35Area.CNR1 <- filter(D35Area.CNR1, Initial == "MY")
D35Area <- filter(D35Area, Experiment != "CNRep1")
D35Area <- filter(D35Area, Experiment != "UNCRep7")

D35Area.CNR2 <- filter(D35Area, Experiment == "CNRep2")
#for CNR2, looks like LTD was just done twice, randomly select one
D35Area.CNR2.LTD <- filter(D35Area.CNR2, Initial == "LTD")
D35Area <- filter(D35Area, Experiment != "CNRep2")
#filter by image name, get list of area measurents in vector
#filter by randomly selected area measurment, join back to new df
CNR2.LTD.image <- c(unique(D35Area.CNR2.LTD$FileName_Original)) 
smalldf <- filter(D35Area.CNR2.LTD, FileName_Original == CNR2.LTD.image[1])
CNR2.LTD.AM <- c(smalldf$D35AreaMM2)
smalldf <- filter(smalldf, D35AreaMM2 == sample(D35AreaMM2,1))
CNR2.LTD.ok <- smalldf

for(i in 2:length(CNR2.LTD.image)){
  smalldf <- filter(D35Area.CNR2.LTD, FileName_Original == CNR2.LTD.image[i])
  CNR2.LTD.AM <- c(smalldf$D35AreaMM2)
  smalldf <- filter(smalldf, D35AreaMM2 == sample(D35AreaMM2,1))
  CNR2.LTD.ok <- rbind(smalldf,CNR2.LTD.ok)
}

D35Area.CHOPR7 <- filter(D35Area, Experiment == "CHOPRep7")
D35Area.CHOPR7 <- filter(D35Area.CHOPR7, Initial == "CP")
D35Area <- filter(D35Area, Experiment != "CHOPRep7")

D35Area.CHOPR5 <- filter(D35Area, Experiment == "CHOPRep5")
D35Area.CHOPR5 <- filter(D35Area.CHOPR5, Initial == "LTDAnn")
D35Area <- filter(D35Area, Experiment != "CHOPRep5")

D35Area.CHOPR2 <- filter(D35Area, Experiment == "CHOPRep2")
#for CHOPR2, looks like LTD was just done twice, randomly select one
D35Area.CHOPR2.LTD <- filter(D35Area.CHOPR2, Initial == "LTD")
D35Area <- filter(D35Area, Experiment != "CHOPRep2")
#filter by image name, get list of area measurents in vector
#filter by randomly selected area measurment, join back to new df
CHOPR2.LTD.image <- c(unique(D35Area.CHOPR2.LTD$FileName_Original)) 
smalldf <- filter(D35Area.CHOPR2.LTD, FileName_Original == CHOPR2.LTD.image[1])
CHOPR2.LTD.AM <- c(smalldf$D35AreaMM2)
smalldf <- filter(smalldf, D35AreaMM2 == sample(D35AreaMM2,1))
CHOPR2.LTD.ok <- smalldf

for(i in 2:length(CHOPR2.LTD.image)){
  smalldf <- filter(D35Area.CHOPR2.LTD, FileName_Original == CHOPR2.LTD.image[i])
  CHOPR2.LTD.AM <- c(smalldf$D35AreaMM2)
  smalldf <- filter(smalldf, D35AreaMM2 == sample(D35AreaMM2,1))
  CHOPR2.LTD.ok <- rbind(smalldf,CHOPR2.LTD.ok)
}

D35Area <- rbind(D35Area,CNR2.LTD.ok,D35Area.CNR1,D35Area.CHOPR7,D35Area.CHOPR5,CHOPR2.LTD.ok)

D35Area <- dplyr::rename(D35Area, AreaMM2 = D35AreaMM2)

saveRDS(D35Area, paste0(rds.here, "IDDRCD35AreaReady.rds"))

# D56 ##################################################################
#Remove tails of distribution
Bottom2p5 <- quantile(D56Area$D56AreaMM2, 0.025)
Top2p5 <- quantile(D56Area$D56AreaMM2, 0.975)
D56Area <- filter(D56Area, D56AreaMM2 < Top2p5)
D56Area <- filter(D56Area, D56AreaMM2 > Bottom2p5)

D56Area <- filter(D56Area, Experiment != "CNRep7")

D56Area <- dplyr::rename(D56Area, AreaMM2 = D56AreaMM2)

saveRDS(D56Area, paste0(rds.here, "IDDRCD56AreaReady.rds"))



