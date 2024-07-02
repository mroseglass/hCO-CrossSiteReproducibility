library("ggsignif")
library(tidyverse)
library(dplyr)
library(readxl)
library(plotrix)
setwd("/work/users/r/o/roseg")

# INPUT FILES ############################################################################################
rds.here <- "/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/"

D7Area <- readRDS(paste0(rds.here, "IDDRCD7Area.rds"))
D14Area <- readRDS(paste0(rds.here, "IDDRCD14Area.rds"))
D35Area <- readRDS(paste0(rds.here, "IDDRCD35Area.rds"))
D56Area <- readRDS(paste0(rds.here, "IDDRCD56Area.rds"))

#Also missing UNC Rep5 Day 56, so add that back in
UNCD56 <- read_csv("/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/UNC_D56_AreaAdd.csv")
UNCD56 <- UNCD56[,-7] 
UNCD56 <- mutate(UNCD56, Experiment = paste0(Site,Replicate)) %>%
  mutate(Mask = Initial)
D56Area <- rbind(UNCD56,D56Area)

# Replace CN images with registered images for each day
D7Area$Site <- gsub("CN","CN_First",D7Area$Site)
D14Area$Site <- gsub("CN","CN_First",D14Area$Site)
D35Area$Site <- gsub("CN","CN_First",D35Area$Site)
D56Area$Site <- gsub("CN","CN_First",D56Area$Site)

#Update Day 7 with registered images 
D7 <- read.csv("/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/1221CNRegisteredImages/D7/7_123456789_8bit_CN_RegisteredIntensity_Final_hCOs.csv")
images_discard <- read_excel("/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/1221CNRegisteredImages/D7/CN_Registered_D7_DoNotKeep.xlsx")
D7 <- anti_join(D7,images_discard, by = "FileName_Intensity")
#remove pictures of teh 12 well plate
D7 <- D7 %>% filter(grepl(')In',FileName_Intensity))
#remove debris
D7 <- filter(D7, AreaShape_Area > 50000)

#Get filtered output to look like rest of D7
D7 <- dplyr::select(D7,c(ObjectNumber,FileName_hCOD36,AreaShape_Area)) %>%
  dplyr::rename(ObjectID = ObjectNumber) %>%
  dplyr::rename(ImageName = FileName_hCOD36) %>% 
  mutate(uniqueID = paste(ImageName,ObjectID, sep = "_")) %>%
  mutate(Replicate = str_extract(ImageName,"CN R[1-9]")) %>%
  mutate(Site = str_extract(ImageName, "CN")) %>%
  mutate(AreaMM2 = AreaShape_Area/0.16/1.31/1000000)
D7$Replicate <- gsub("CN R","Rep",D7$Replicate)
Initial <- rep("Reg", length(D7$ImageName))
D7 <- cbind(D7,Initial)
D7 <- mutate(D7, Experiment = paste0(Site,Replicate))
D7Area <- rbind(D7,D7Area)

#Update Day 14 with registered images
list.cp.output <- list.files("/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/CN_Registered_NoDownSample_CPoutput/D14/")
i=1
CP.out <- read.csv(paste0("/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/CN_Registered_NoDownSample_CPoutput/D14/",list.cp.output[i]))
for (i in 2:length(list.cp.output)){
  CP.out.tmp <- read.csv(paste0("/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/CN_Registered_NoDownSample_CPoutput/D14/",list.cp.output[i]))
  CP.out <- rbind(CP.out,CP.out.tmp)
}

#filters for getting hCO size objects
CP.out.filter <- filter(CP.out, AreaShape_Area > 100000)
CP.out.filter <- filter(CP.out.filter, AreaShape_Area < 1000000)

#Get filtered output to look like rest of D14
CP_out_slim <- mutate(CP.out.filter, Area_ID = paste(ObjectNumber,FileName_Edges))  
CN_D14_ready <- dplyr::select(CP_out_slim, c(Area_ID,AreaShape_Area))
Site <- rep("CN", length(CN_D14_ready$Area_ID))
CN_D14_ready <- cbind(Site,CN_D14_ready)
CN_D14_ready <- mutate(CN_D14_ready, Replicate = str_extract(CN_D14_ready$Area_ID, "N R."))
CN_D14_ready$Replicate <- gsub("N R","Rep",CN_D14_ready$Replicate)
CN_D14_ready <- mutate(CN_D14_ready, Experiment = paste0(Site,Replicate))
Assignment <- rep("Registered", length(CN_D14_ready$Area_ID))
Initial <- rep("RG", length(CN_D14_ready$Area_ID))
CN_D14_ready <- cbind(CN_D14_ready,Assignment,Initial)
CN_D14_ready <- dplyr::rename(CN_D14_ready, FileName_D14Original = Area_ID)
CN_D14_ready <- dplyr::rename(CN_D14_ready, AreaOccupied_Area = AreaShape_Area)
CN_D14_ready <- mutate(CN_D14_ready, D14AreaMM2 = AreaOccupied_Area/0.16/1.31/1000000) 
D14Area <- rbind(CN_D14_ready,D14Area)

#Day 35 CN no downsamle registered images
list.cp.output <- list.files("/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/CN_Registered_NoDownSample_CPoutput/D35/")
i=1
CP.out <- read.csv(paste0("/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/CN_Registered_NoDownSample_CPoutput/D35/",list.cp.output[i]))
for (i in 2:length(list.cp.output)){
  CP.out.tmp <- read.csv(paste0("/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/CN_Registered_NoDownSample_CPoutput/D35/",list.cp.output[i]))
  CP.out <- rbind(CP.out,CP.out.tmp)
}

#filters for getting hCO size objects
CP.out.filter <- filter(CP.out, AreaShape_Area > 100000)
CP.out.filter <- filter(CP.out.filter, AreaShape_Area < 1000000)

D35 <- dplyr::select(CP.out.filter, c(FileName_hCOD36,AreaShape_Area)) %>%
  dplyr::rename(FileName_Original = FileName_hCOD36) %>%
  mutate(D35AreaMM2 = AreaShape_Area/0.16/1.31/1000000) %>%
  mutate(Replicate = str_extract(FileName_Original, "CNR[1-9]")) %>%
  mutate(Site = str_extract(FileName_Original, "CN"))
D35$Replicate <- gsub("CNR","Rep",D35$Replicate)
Initial <- rep("Reg", length(D35$FileName_Original))
D35 <- cbind(D35,Initial)
D35 <- mutate(D35, Experiment = paste0(Site,Replicate))
D35Area <- rbind(D35,D35Area)

#Day 56 CN no downsamle registered images
list.cp.output <- list.files("/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/1221CNRegisteredImages/D56/CND56_JusttheGoods/")
i=1
CP.out <- read.csv(paste0("/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/1221CNRegisteredImages/D56/CND56_JusttheGoods/",list.cp.output[i]))
for (i in 2:length(list.cp.output)){
  CP.out.tmp <- read.csv(paste0("/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/1221CNRegisteredImages/D56/CND56_JusttheGoods/",list.cp.output[i]))
  CP.out <- rbind(CP.out,CP.out.tmp)
}
#filter CP.out by images to keep
CP.out <- mutate(CP.out, UniqueID = paste(ObjectNumber,FileName_hCOD36, sep = "v"))
CP.out.filter <- read.csv("/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/1221CNRegisteredImages/D56/Keep_Edges_CN_D56_NoDownsample.csv") %>%
  mutate(UniqueID = paste(ObjectNumber,FileName_hCOD36, sep = "v")) %>%
  select(UniqueID)
CP.out.filter$UniqueID <- gsub("Edges.jpeg",".jpg", CP.out.filter$UniqueID)
#notworkingwell
CP.out.filter <- inner_join(CP.out.filter,CP.out, by = "UniqueID")
CP.out.filter <- select(CP.out.filter, c(FileName_hCOD36, AreaShape_Area))
Mask <- rep("Manual", length(CP.out.filter$FileName_hCOD36))
CP.out.filter <- cbind(Mask,CP.out.filter)

#seperately add the retrace files
list.cp.output <- list.files("/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/1221CNRegisteredImages/D56/CN_D56_Retrace/")
i=1
CP.out <- read.csv(paste0("/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/1221CNRegisteredImages/D56/CN_D56_Retrace/",list.cp.output[i]))
for (i in 2:length(list.cp.output)){
  CP.out.tmp <- read.csv(paste0("/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/1221CNRegisteredImages/D56/CN_D56_Retrace/",list.cp.output[i]))
  CP.out <- rbind(CP.out,CP.out.tmp)
}
CP.out <- filter(CP.out, AreaShape_Area > 3000) %>%
  filter(AreaShape_Area != "NaN") %>% 
  select(c(FileName_hCOD36, AreaShape_Area))
Mask <- rep("Manual", length(CP.out$FileName_hCOD36))
CP.out <- cbind(Mask,CP.out)

D56 <- rbind(CP.out.filter,CP.out)
D56 <-   dplyr::rename(D56,FileName_D56Original = FileName_hCOD36) %>%
  mutate(D56AreaMM2 = AreaShape_Area/0.16/1.31/1000000) %>%
  mutate(Replicate = str_extract(FileName_D56Original, "CN R[1-9]")) %>%
  mutate(Site = str_extract(FileName_D56Original, "CN"))
D56$Replicate <- gsub("CN R","Rep",D56$Replicate)
Initial <- rep("Reg", length(D56$FileName_D56Original))
D56 <- cbind(D56,Initial)
D56 <- mutate(D56, Experiment = paste0(Site,Replicate)) %>% dplyr::rename(AreaOccupied_D56 = AreaShape_Area)
D56Area <- rbind(D56,D56Area)

# Save Updated Area
saveRDS(D7Area, paste0(rds.here, "CNRegistered_D7Area.rds"))
saveRDS(D14Area, paste0(rds.here, "CNRegistered_D14Area.rds"))
saveRDS(D35Area,paste0(rds.here, "CNRegistered_D35Area.rds"))
saveRDS(D56Area, paste0(rds.here, "CNRegistered_D56Area.rds"))

#pull together scRNAseq related areas
area.scrnaseq.84 <- read_csv("/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/IDDRC_D84_scRNAseq_hCO_area.csv")
area.scrnaseq.84 <- dplyr::rename(area.scrnaseq.84,AreaOccupied_Area = AreaShape_Area) %>%
  dplyr::rename(Replicate = Rep) %>%
  filter(Site != "CN") %>%
  select(c(Site,Replicate,AreaOccupied_Area,FileName_hCOD36)) %>%
  mutate(AreaMM2 = AreaOccupied_Area/0.16/1.31/1000000)
area.scrnaseq.84$Replicate <- gsub("R","Rep", area.scrnaseq.84$Replicate)

cn.update <- read_csv("/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/1221CNRegisteredImages/CN_scRNAseq_D84/scRNAseq_RegisteredhCOD36_Trace_KeepersOnly.csv")
cn.update <- dplyr::select(cn.update, c(FileName_hCOD36,AreaShape_Area))
cn.update <-  dplyr::rename(cn.update,AreaOccupied_Area = AreaShape_Area) %>%
  mutate(AreaMM2 = AreaOccupied_Area/0.16/1.31/1000000) %>%
  mutate(Replicate = str_extract(FileName_hCOD36, "CN R[1-9]")) %>%
  mutate(Site = str_extract(FileName_hCOD36, "CN"))
cn.update$Replicate <- gsub("CN R","Rep",cn.update$Replicate)

area.scrnaseq.84 <- rbind(area.scrnaseq.84,cn.update) %>%
  dplyr::rename(FileName = FileName_hCOD36)
area.scrnaseq.84$Replicate <- gsub("Rep 1","Rep1", area.scrnaseq.84$Replicate)
Day <- rep("D84", length(area.scrnaseq.84$FileName))
area.scrnaseq.84 <- cbind(Day,area.scrnaseq.84)

#add in CN D14
area.scrnaseq.14 <- D14Area[grep("scRNA",D14Area$FileName_D14Original),] %>%
  dplyr::rename(FileName = FileName_D14Original) %>%
  dplyr::rename(AreaMM2 = D14AreaMM2)
UNCR5 <- c("UNC","IDDRC UNC Rep5 Plate2 CB scRNA  D14 4.8.22 0178.jpg",103128,"Rep5","UNCRep5","Added","RG",103128/0.16/1.31/1000000)
UNCR3 <- c("UNC","UNC_R3_P2_CoA scrnaseq d14 10.15.21 0048.jpg",120894,"Rep3","UNCRep3","Added","RG",120894/0.16/1.31/1000000)
CHOPR4 <- c("CHOP","20211018-EW38R4-scRNA-1.tif",167017,"Rep4","CHOPRep4","Added","RG",167017/0.16/1.31/1000000)
CHOPR5 <- c("CHOP","20220323-ew42r5-c6-scrna-1.tif",200894,"Rep5","CHOPRep5","Added","RG",200894/0.16/1.31/1000000)
CHOPR6 <- c("CHOP","20220324-ew42r6pl1-day14-c3-scrna-1.tif",207004,"Rep6","CHOPRep6","Added","RG",207004/0.16/1.31/1000000)
area.scrnaseq.14 <- rbind(area.scrnaseq.14,UNCR5,UNCR3,CHOPR6,CHOPR5,CHOPR4)
Day <- rep("D14", length(area.scrnaseq.14$FileName))
area.scrnaseq.14 <- cbind(Day,area.scrnaseq.14)

area.scrnaseq.14 <- select(area.scrnaseq.14,colnames(area.scrnaseq.84))
area.scrnaseq <- rbind(area.scrnaseq.14,area.scrnaseq.84)

write_csv(area.scrnaseq,"/work/users/r/o/roseg/IDDRC/IDDRCDatabase/CrossSectionalAreaRNAseq.csv")
#########################################################################################
ggplot(CP.out, aes(x=AreaShape_Area))+
  geom_histogram()

#filters for getting hCO size objects
#less than 4million
#larger than 100,000?
CP.out.filter <- filter(CP.out, AreaShape_Area > 100000)
CP.out.filter <- filter(CP.out.filter, AreaShape_Area < 1000000)




ggplot(CP.out, aes(x=AreaShape_Area))+
  geom_histogram()

ggplot(D14Area, aes(x=Site, y=D14AreaMM2)) +
  geom_boxplot(outlier.shape = NA, lwd=0.2) + 
  geom_point(aes(color=Replicate),position=position_jitterdodge(), alpha=0.9)+
  theme(axis.text.x = element_text(angle = 45)) +
  geom_signif(test = "t.test",comparisons = list(c("CN", "CHOP"),
                                                 c("CN", "UNC"),
                                                 c("UNC","CHOP")),
              y_position = c(2,3,4))

