####START####

####Set working directory under your home directory ('R'is my home directory)
setwd("~/R")

###Lake project###
#(6 Cheshire meres & 8 Welsh lakes)

###packages require in this script###
##The package(reshape) should after package(tidyverse), otherwise the rename funtion in package(reshape) will not work

#plyr
library(plyr)
#tidyverse
library(tidyverse)
#cowplot
library(cowplot)
#scales
library(scales)
#reshape
library(reshape)
library(reshape2)
#gtable
library(gtable)
#data.table
library(data.table)
#gridExtra
library(gridExtra)
#ggrepel
library(ggrepel)
#Using ggplot for the NMDS plot https://chrischizinski.github.io/rstats/vegan-ggplot2/
library(vegan)
#stat_cor
library (ggpubr)



###metaBEAT command###
####12S

#metaBEAT_global.py \
#-Q Querymap.txt \
#-R REFmap.txt \
#--cluster --clust_match 1 --clust_cov 3 \
#--blast --min_ident 1 \
#-m 12S -n 5 \
#-E -v \
#-@  haikuilee@gmail.com \
#-o 12S_LP-trim_30-merged-nonchimeras-cluster_1c3-blast-min_ident_1 &> log_non-chimeras

####Cytb

#metaBEAT_global.py \
#-Q Querymap.txt \
#-R REFmap.txt \
#--cluster --clust_match 1 --clust_cov 3 \
#--blast --min_ident 0.95 \
#-m Cytb -n 5 \
#-E -v \
#-@  haikuilee@gmail.com \
#-o Cytb_LP-trim_30-merged-nonchimeras-cluster_1c3-blast-min_ident_0.95 &> log_non-chimeras

###ID description###

####Add more ID information into the metaBEAT results
#Total_n (Total sample numbers from a site, if the number is same for all of sites, then no need to add)
#OSL_n (Offshore line sample numbers from a lake, if the number is are same for all of sites, then no need to add)
#SL_n (Shore line sample numbers from a lake, if the number is are same for all of sites, then no need to add)
#Lake (name of lake or mere)
#Project (Lake or Mere)
#Transect OR Location (ShoreLine or Off Shore line)	
#Lat(latitude)
#Log(longitude)


###Cytb####

#Read the csv file with ID information
LP_Cytb <- read.csv(file="Lake project/Cytb_LP_Rdataset_Dec_2017.csv", header = TRUE)
         

##add the reads of 'Salmo' into 'Salmo_trutta', due to only 'Salmo_trutta' under the genus 'Salmo'
##for Llyn_Cwellyn; Llyn_Ogwen; Llyn_Traffwll
##keep 'Salmo' in Llyn_Padarn, due to two species 'Salmo_trutta' and 'Salmo_salar' under the genus 'Salmo', need phylogenetic tree to solve this problem

for (i in 1:nrow(LP_Cytb)){
   if (LP_Cytb$Lake[i] == "Llyn_Cwellyn" | LP_Cytb$Lake[i] == "Llyn_Ogwen" | LP_Cytb$Lake[i] == "Llyn_Traffwll")
       {LP_Cytb$Salmo_trutta[i] <- LP_Cytb$Salmo_trutta[i]+LP_Cytb$Salmo[i]; 
                           LP_Cytb$Salmo[i]<- 0}
}


##add 'Salvelinus' into 'Salvelinus_alpinus' and delete Salvelinus
for (i in 1:nrow(LP_Cytb)){LP_Cytb$Salvelinus_alpinus[i]=LP_Cytb$Salvelinus_alpinus[i]+LP_Cytb$Salvelinus[i]; 
  LP_Cytb$Salvelinus[i]<-0}
LP_Cytb$Salvelinus <- NULL

##add 'Oncorhynchus' into 'Oncorhynchus_mykiss and delete Oncorhynchus
for (i in 1:nrow(LP_Cytb)){LP_Cytb$Oncorhynchus_mykiss[i]=LP_Cytb$Oncorhynchus_mykiss[i]+LP_Cytb$Oncorhynchus[i]; 
LP_Cytb$Oncorhynchus[i]<-0}

LP_Cytb$Oncorhynchus <- NULL

#Summary the the total read counts for each sample
LP_Cytb$SUM <- rowSums(LP_Cytb[12:35])

LP_Cytb_sample <- LP_Cytb[which(LP_Cytb$Transect != 'Control'),]
  

LP_Cytb_sample_Total <- aggregate(SUM~Lake+Total_n+Project+Locus, data = LP_Cytb_sample, FUN=sum)

#Reshape the dataset


#Using the 1:11, and 36 variables as ID

LP_Cytb_rs <- melt(LP_Cytb,id=c(1:11,36))

#Rename the variables

LP_Cytb_rs <- rename (LP_Cytb_rs,c("variable"="Species","value"="Reads"))



#Calculate each species reads percentage for each sample

LP_Cytb_rs$Ratio <- LP_Cytb_rs$Reads/LP_Cytb_rs$SUM

LP_Cytb_rs <-LP_Cytb_rs[which(LP_Cytb_rs$Reads >0),]

###First threshold####

###Avoid false positive approaches###
#Low-frequency noise threshold (0.0007) based on the compare results under different thresholds (0.0004, 0.0007, 0.0001)
#If you want to decrease the value of threshold, you need to re-read the original CSV file

###Threshold_Cytb###

Threshold_Cytb <- 0.0007



#if the Ratio less than the Threshold, then the reads change to zero


for (a in 1:nrow(LP_Cytb_rs)) {
  if (!is.na(LP_Cytb_rs$Ratio[a])) {
    if(LP_Cytb_rs$Ratio[a] < Threshold_Cytb){LP_Cytb_rs$Reads[a] <- 0} 
  }
}

###POS###

#Check the Control samples, then remove all of control samples, positive species and unassigned species

LP_Cytb_rs_control <- LP_Cytb_rs[which(LP_Cytb_rs$Transect == 'Control'),]

#Check the species read counts in Control samples

LP_Cytb_rs_control_filter <- LP_Cytb_rs_control[which(LP_Cytb_rs_control$Species != 'Astatotilapia_calliptera' 
                                                      & LP_Cytb_rs_control$Species != 'unassigned'),]


ggplot(LP_Cytb_rs_control_filter, aes(x = Species, y= Reads, fill = Control))+
  geom_bar(stat="identity",position="dodge")+facet_wrap(~Lake, ncol = 5, scales = "free")+
  labs(x="Species", y="Read counts", title= "Control samples_Cytb_T0.0007")+theme_bw()+
  theme(text=element_text(size=20),axis.text.x = element_text(angle =45, vjust = 1, hjust = 0.95),
        plot.title = element_text(size = rel(0.8),face = "bold", hjust = 0.5))




ggsave("Control samples_Cytb_T0.0007.jpeg",path = "Lake project/Figures_new/", width = 10, height = 8, units = "in", dpi=500)

###Samples###


LP_Cytb_rs_sample <- LP_Cytb_rs[which(LP_Cytb_rs$Transect != 'Control' 
                    & LP_Cytb_rs$Species != 'Astatotilapia_calliptera' 
                    & LP_Cytb_rs$Species != 'unassigned' 
                    & LP_Cytb_rs$Species != 'Cyprinidae'
                    & LP_Cytb_rs$Species != 'Percidae'
                    & LP_Cytb_rs$Species != 'Salmonidae'),]



###12S####
##no need to add any species reads from unassigned OTU BLAST

#Read the csv file with ID information
LP_12S <- read.csv(file="Lake project/12S_LP_Rdataset_Dec_2017.csv", header = TRUE)


#Summary the the total read counts for each sample
LP_12S$SUM <- rowSums(LP_12S[12:40])


LP_12S_sample <- LP_12S[which(LP_12S$Transect != 'Control'),]

LP_12S_sample_Total <- aggregate(SUM~Lake+Total_n+Project+Locus, data = LP_12S_sample, FUN=sum)

#Reshape the dataset


#Using the 1:11, and 41 variables as ID

LP_12S_rs <- melt(LP_12S,id=c(1:11,41))

#Rename the variables

LP_12S_rs <- rename (LP_12S_rs,c("variable"="Species","value"="Reads"))

LP_12S_rs <-LP_12S_rs[which(LP_12S_rs$Reads >0),]

#Calculate each species reads percentage for each sample

LP_12S_rs$Ratio <- LP_12S_rs$Reads/LP_12S_rs$SUM

###First threshold####

#Low-frequency noise threshold (0.003) based on the compare results under different thresholds (0.001,0.002,0.003,0.004,0.005)
#if the Ratio less than the Threshold, then the reads change to zero
Threshold_12S <- 0.003


for (b in 1:nrow(LP_12S_rs)) {
  if (!is.na(LP_12S_rs$Ratio[b])) {
    if(LP_12S_rs$Ratio[b] < Threshold_12S){LP_12S_rs$Reads[b] <- 0} 
  }
}

###Control samples###
#Check the Control samples, then remove all of control samples, positive species and unassigned species

LP_12S_rs_control <- LP_12S_rs[which(LP_12S_rs$Transect == 'Control'
                                     &LP_12S_rs$Control != 'POS1'),]


LP_12S_rs_control_filter <- LP_12S_rs_control[which(LP_12S_rs_control$Species != 'Astatotilapia_calliptera' 
                                                      & LP_12S_rs_control$Species != 'unassigned'
                                                    & LP_12S_rs_control$Species != 'Cyprinidae'),]


ggplot(LP_12S_rs_control_filter, aes(x = Species, y= Reads, fill = Control))+
  geom_bar(stat="identity",position="dodge")+facet_wrap(~Lake,ncol =5, scales = "free")+
  labs(x="Species", y="Read counts", title= "Control samples_12S_T0.003")+theme_bw()+
  theme(text=element_text(size=14),axis.text.x = element_text(angle =45, vjust = 1, hjust = 0.95),
        plot.title = element_text(size = rel(0.8),face = "bold", hjust = 0.5))


ggsave("Control samples_12S_noPOS_T0.003.jpeg",path = "Lake project/Figures_new/", width = 10, height = 8, units = "in", dpi=500)



LP_12S_rs_control_pos <- LP_12S_rs[which(LP_12S_rs$Transect == 'Control'),]


LP_12S_rs_control_filter_2 <- LP_12S_rs_control_pos[which(LP_12S_rs_control_pos$Species != 'Astatotilapia_calliptera' 
                                                    & LP_12S_rs_control_pos$Species != 'unassigned'
                                                    & LP_12S_rs_control_pos$Species != 'Cyprinidae'),]


ggplot(LP_12S_rs_control_filter_2, aes(x = Species, y= Reads, fill = Control))+
  geom_bar(stat="identity",position="dodge")+facet_wrap(~Lake,ncol =5, scales = "free")+
  labs(x="Species", y="Read counts", title= "Control samples_12S_T0.003")+theme_bw()+
  theme(text=element_text(size=14),axis.text.x = element_text(angle =45, vjust = 1, hjust = 0.95),
        plot.title = element_text(size = rel(0.8),face = "bold", hjust = 0.5))


ggsave("Control samples_12S_T0.003.jpeg",path = "Lake project/Figures_new/", width = 10, height = 8, units = "in", dpi=500)




###Samples###
LP_12S_rs_sample <- LP_12S_rs[which(LP_12S_rs$Transect != 'Control' 
                                    & LP_12S_rs$Species != 'Astatotilapia_calliptera' 
                                    & LP_12S_rs$Species != 'unassigned'),]



###Combine the Cytb and 12S sample dataset
LP_rs <- rbind(LP_Cytb_rs_sample,LP_12S_rs_sample)

LP_Sample_Total <- rbind(LP_Cytb_sample_Total,LP_12S_sample_Total) #Total reads of all samples from this lake

LP_rs_mere <- LP_rs[which(LP_rs$Project == 'Mere'),]


LP_rs_wlake <- LP_rs[which(LP_rs$Project == 'Lake'),]




####NMDS_12S####

names(LP_12S_rs_sample)

NMDS_12S <- LP_12S_rs_sample [c("SampleID","Lake","Transect","Species","Reads" )]

NMDS_12S <- NMDS_12S[which(NMDS_12S$Reads >0),]


levels(droplevels(NMDS_12S$Species))

NMDS_12S <- NMDS_12S[which(NMDS_12S$Species != 'Thymallus_thymallus'
                           &NMDS_12S$Species != 'Lampetra_fluviatilis'
                           &NMDS_12S$Species != 'Leuciscus_idus'
                           &NMDS_12S$Species != 'Blicca_bjoerkna'),]




NMDS_12S <- NMDS_12S[which(NMDS_12S$Lake !="Llyn_Cwellyn" |
                             NMDS_12S$Species != "Rutilus_rutilus"),]


NMDS_12S <- NMDS_12S[which(NMDS_12S$Lake !="Llyn_Ogwen" |
                             NMDS_12S$Species != "Abramis_brama"),]

NMDS_12S <- NMDS_12S[which(NMDS_12S$Lake !="Llyn_Ogwen" |
                             NMDS_12S$Species != "Cottus_gobio"),]


NMDS_12S <- NMDS_12S[which(NMDS_12S$Lake !="Llyn_Penrhyn" |
                             NMDS_12S$Species != "Esox_lucius"),]



NMDS_12S <- NMDS_12S[which(NMDS_12S$Lake !="Llyn_Traffwll" |
                             NMDS_12S$Species != "Cottus_gobio"),]



NMDS_12S <- NMDS_12S[which(NMDS_12S$Lake !="Betley_Mere" |
                             NMDS_12S$Species != "Gobio_gobio"),]





NMDS_12S$Species <-mapvalues(NMDS_12S$Species, c("Anguilla_anguilla","Abramis_brama","Carassius_carassius","Cyprinus_carpio", 
                                                 "Leucaspius_delineatus","Phoxinus_phoxinus","Rutilus_rutilus","Scardinius_erythrophthalmus",
                                                 "Tinca_tinca","Esox_lucius","Cottus_gobio","Gasterosteus_aculeatus",     
                                                 "Perca_fluviatilis","Salmo_trutta","Gobio_gobio","Pseudorasbora_parva",        
                                                 "Rhodeus_amarus","Squalius_cephalus","Alburnus_alburnus","Pungitius_pungitius",
                                                 "Oncorhynchus_mykiss","Salmo_salar","Salvelinus_alpinus"), 
                             c("EEL","BRE","CRU","CAR",
                               "SUN","MIN","ROA","RUD",
                               "TEN","PIK","BUL","3SS",
                               "PER","BTR","GUD","TMG",
                               "BIT","CHU","BLE","9SS",
                               "RTR","SAL","CHA"))


NMDS_12S$Lake <-mapvalues(NMDS_12S$Lake, c("Betley_Mere","Chapel_Mere","Fenemere","Kenfig_Pool","Llan_Bwch-llyn","Llangorse_Lake",  
                                           "Llyn_Cwellyn","Llyn_Ogwen","Llyn_Padarn","Llyn_Penrhyn","Llyn_Traffwll","Maer_Pool",       
                                           "Oss_Mere","Watch_Lane_Flash"), 
                          c("BET","CAM","FEN","KEN","LLB","LLG","CWE","OGW","PAD","PEN","TRA","MAP","OSS","WLF"))



NMDS_12S_wide <- reshape(NMDS_12S, idvar= c("SampleID","Lake","Transect"), timevar= "Species", direction = "wide")




names(NMDS_12S_wide)

NMDS_12S_wide <- rename (NMDS_12S_wide,c("Reads.EEL"="EEL","Reads.BRE"="BRE","Reads.BLE"="BLE","Reads.CRU"="CRU",
                                         "Reads.CAR"="CAR","Reads.GUD"="GUD","Reads.SUN"="SUN","Reads.MIN"="MIN",
                                         "Reads.TMG"="TMG","Reads.BIT"="BIT","Reads.ROA"="ROA","Reads.RUD"="RUD",
                                         "Reads.CHU"="CHU","Reads.TEN"="TEN","Reads.PIK"="PIK","Reads.BUL"="BUL","Reads.3SS"="3SS",
                                         "Reads.9SS"="9SS","Reads.RTR"="RTR","Reads.SAL"="SAL","Reads.BTR"="BTR","Reads.CHA"="CHA","Reads.PER"="PER"))


NMDS_12S_wide.rownames <- data.frame(NMDS_12S_wide[,-1], row.names=NMDS_12S_wide[,1])


NMDS_12S_wide.rownames[is.na(NMDS_12S_wide.rownames)] <- 0



NMDS_12S_wide.active <- NMDS_12S_wide.rownames[, 3:ncol(NMDS_12S_wide.rownames)]




NMDS_12S_data <- metaMDS(NMDS_12S_wide.active, distance="bray", k=2, trymax=1000)


#check stressplot
stressplot(NMDS_12S_data)






#Using the scores function from vegan to extract the site scores and convert to a data.frame

NMDS_12S_data.scores <- as.data.frame(scores(NMDS_12S_data))

# create a column of sampleID names, from the rownames of data.scores
NMDS_12S_data.scores$Rep <- rownames(NMDS_12S_data.scores)

#  add the Treatment variable created earlier
NMDS_12S_data.scores$Lake <- NMDS_12S_wide.rownames$Lake


NMDS_12S_data.scores$Transect <- NMDS_12S_wide.rownames$Transect

NMDS_12S_data.scores$Locus <- "sRNA"

#look at the data
head(NMDS_12S_data.scores)



#Using the scores function from vegan to extract the species scores and convert to a data.frame
NMDS_12S_species.scores <- as.data.frame(scores(NMDS_12S_data, "species"))

# create a column of species, from the rownames of species.scores
NMDS_12S_species.scores$Species <- rownames(NMDS_12S_species.scores)

NMDS_12S_species.scores$Locus <- "sRNA"

#look at the data
head(NMDS_12S_species.scores)



NMDS_12S_species.scores$Species <-mapvalues(NMDS_12S_species.scores$Species, c("X3SS","X9SS"), 
                                                                            c("3SS","9SS"))



#Using the chull function from vegan to extract the hull values and convert to a data.frame

NMDS_12S.PAD <- NMDS_12S_data.scores[NMDS_12S_data.scores$Lake == "PAD", ][chull(NMDS_12S_data.scores[NMDS_12S_data.scores$Lake == 
                                                                                                        "PAD", c("NMDS1", "NMDS2")]), ]
NMDS_12S.CWE <- NMDS_12S_data.scores[NMDS_12S_data.scores$Lake == "CWE", ][chull(NMDS_12S_data.scores[NMDS_12S_data.scores$Lake == 
                                                                                                        "CWE", c("NMDS1", "NMDS2")]), ]

NMDS_12S.LLG <- NMDS_12S_data.scores[NMDS_12S_data.scores$Lake == "LLG", ][chull(NMDS_12S_data.scores[NMDS_12S_data.scores$Lake == 
                                                                                                        "LLG", c("NMDS1", "NMDS2")]), ]

NMDS_12S.TRA <- NMDS_12S_data.scores[NMDS_12S_data.scores$Lake == "TRA", ][chull(NMDS_12S_data.scores[NMDS_12S_data.scores$Lake == 
                                                                                                        "TRA", c("NMDS1", "NMDS2")]), ]

NMDS_12S.OGW <- NMDS_12S_data.scores[NMDS_12S_data.scores$Lake == "OGW", ][chull(NMDS_12S_data.scores[NMDS_12S_data.scores$Lake == 
                                                                                                        "OGW", c("NMDS1", "NMDS2")]), ]

NMDS_12S.KEN <- NMDS_12S_data.scores[NMDS_12S_data.scores$Lake == "KEN", ][chull(NMDS_12S_data.scores[NMDS_12S_data.scores$Lake == 
                                                                                                        "KEN", c("NMDS1", "NMDS2")]), ]

NMDS_12S.PEN <- NMDS_12S_data.scores[NMDS_12S_data.scores$Lake == "PEN", ][chull(NMDS_12S_data.scores[NMDS_12S_data.scores$Lake == 
                                                                                                        "PEN", c("NMDS1", "NMDS2")]), ]

NMDS_12S.LLB <- NMDS_12S_data.scores[NMDS_12S_data.scores$Lake == "LLB", ][chull(NMDS_12S_data.scores[NMDS_12S_data.scores$Lake == 
                                                                                                        "LLB", c("NMDS1", "NMDS2")]), ]

NMDS_12S.WLF <- NMDS_12S_data.scores[NMDS_12S_data.scores$Lake == "WLF", ][chull(NMDS_12S_data.scores[NMDS_12S_data.scores$Lake == 
                                                                                                        "WLF", c("NMDS1", "NMDS2")]), ]

NMDS_12S.OSS <- NMDS_12S_data.scores[NMDS_12S_data.scores$Lake == "OSS", ][chull(NMDS_12S_data.scores[NMDS_12S_data.scores$Lake == 
                                                                                                        "OSS", c("NMDS1", "NMDS2")]), ]

NMDS_12S.FEN <- NMDS_12S_data.scores[NMDS_12S_data.scores$Lake == "FEN", ][chull(NMDS_12S_data.scores[NMDS_12S_data.scores$Lake == 
                                                                                                        "FEN", c("NMDS1", "NMDS2")]), ]

NMDS_12S.MAP <- NMDS_12S_data.scores[NMDS_12S_data.scores$Lake == "MAP", ][chull(NMDS_12S_data.scores[NMDS_12S_data.scores$Lake == 
                                                                                                        "MAP", c("NMDS1", "NMDS2")]), ]

NMDS_12S.BET <- NMDS_12S_data.scores[NMDS_12S_data.scores$Lake == "BET", ][chull(NMDS_12S_data.scores[NMDS_12S_data.scores$Lake == 
                                                                                                        "BET", c("NMDS1", "NMDS2")]), ]

NMDS_12S.CAM <- NMDS_12S_data.scores[NMDS_12S_data.scores$Lake == "CAM", ][chull(NMDS_12S_data.scores[NMDS_12S_data.scores$Lake == 
                                                                                                        "CAM", c("NMDS1", "NMDS2")]), ]




NMDS_12S_group.scores <- NMDS_12S_data.scores


NMDS_12S_group.scores$Lake <-mapvalues(NMDS_12S_group.scores$Lake, c("CWE","OGW","PAD","TRA","PEN","KEN","MAP","CAM","LLB","LLG","OSS","FEN","WLF","BET"), 
                                c("G1","G1","G1","G2","G2","G3","G3","G3","G3","G3","G3","G3","G4","G4"))


NMDS_12S.G1 <- NMDS_12S_group.scores[NMDS_12S_group.scores$Lake == "G1", ][chull(NMDS_12S_group.scores[NMDS_12S_group.scores$Lake == 
                                                                                                         "G1", c("NMDS1", "NMDS2")]), ]
NMDS_12S.G2 <- NMDS_12S_group.scores[NMDS_12S_group.scores$Lake == "G2", ][chull(NMDS_12S_group.scores[NMDS_12S_group.scores$Lake == 
                                                                                                         "G2", c("NMDS1", "NMDS2")]), ]
NMDS_12S.G3 <- NMDS_12S_group.scores[NMDS_12S_group.scores$Lake == "G3", ][chull(NMDS_12S_group.scores[NMDS_12S_group.scores$Lake == 
                                                                                                         "G3", c("NMDS1", "NMDS2")]), ]
NMDS_12S.G4 <- NMDS_12S_group.scores[NMDS_12S_group.scores$Lake == "G4", ][chull(NMDS_12S_group.scores[NMDS_12S_group.scores$Lake == 
                                                                                                         "G4", c("NMDS1", "NMDS2")]), ]

NMDS_12S_hull.group <- rbind(NMDS_12S.G1,NMDS_12S.G2,NMDS_12S.G3,NMDS_12S.G4)
NMDS_12S_hull.group$Locus <- "sRNA"


#combine sampling site data

NMDS_12S_hull.data <- rbind(NMDS_12S.PAD,NMDS_12S.CWE,NMDS_12S.LLG,NMDS_12S.TRA,NMDS_12S.OGW,NMDS_12S.KEN,NMDS_12S.PEN,
                            NMDS_12S.LLB,NMDS_12S.WLF,NMDS_12S.OSS,NMDS_12S.FEN,NMDS_12S.MAP,NMDS_12S.BET,NMDS_12S.CAM)  
NMDS_12S_hull.data$Locus <- "sRNA"


####NMDS_Cytb####



names(LP_Cytb_rs_sample)

NMDS_Cytb <- LP_Cytb_rs_sample [c("SampleID","Lake","Transect","Species","Reads" )]

NMDS_Cytb <- NMDS_Cytb[which(NMDS_Cytb$Reads >0),]


levels(droplevels(NMDS_Cytb$Species))

NMDS_Cytb <- NMDS_Cytb[which(NMDS_Cytb$Species != 'Thymallus_thymallus'
                           &NMDS_Cytb$Species != 'Lampetra_fluviatilis'
                           &NMDS_Cytb$Species != 'Leuciscus_idus'
                           &NMDS_Cytb$Species != 'Blicca_bjoerkna'
                           &NMDS_Cytb$Species != 'Salmo'),]





NMDS_Cytb <- NMDS_Cytb[which(NMDS_Cytb$Lake !="Llyn_Cwellyn" |
                             NMDS_Cytb$Species != "Rutilus_rutilus"),]


NMDS_Cytb <- NMDS_Cytb[which(NMDS_Cytb$Lake !="Llyn_Ogwen" |
                             NMDS_Cytb$Species != "Abramis_brama"),]

NMDS_Cytb <- NMDS_Cytb[which(NMDS_Cytb$Lake !="Llyn_Ogwen" |
                             NMDS_Cytb$Species != "Cottus_gobio"),]


NMDS_Cytb <- NMDS_Cytb[which(NMDS_Cytb$Lake !="Llyn_Penrhyn" |
                             NMDS_Cytb$Species != "Esox_lucius"),]


NMDS_Cytb <- NMDS_Cytb[which(NMDS_Cytb$Lake !="Llyn_Traffwll" |
                             NMDS_Cytb$Species != "Cottus_gobio"),]




NMDS_Cytb$Species <-mapvalues(NMDS_Cytb$Species, c("Anguilla_anguilla","Abramis_brama","Carassius_carassius","Cyprinus_carpio", 
                                                 "Leucaspius_delineatus","Phoxinus_phoxinus","Rutilus_rutilus","Scardinius_erythrophthalmus",
                                                 "Tinca_tinca","Esox_lucius","Cottus_gobio","Gasterosteus_aculeatus",     
                                                 "Perca_fluviatilis","Salmo_trutta","Carassius_auratus",        
                                                 "Oncorhynchus_mykiss","Salmo_salar","Salvelinus_alpinus"), 
                                                 c("EEL","BRE","CRU","CAR",
                                                   "SUN","MIN","ROA","RUD",
                                                   "TEN","PIK","BUL","3SS",
                                                   "PER","BTR","GOF",
                                                   "RTR","SAL","CHA"))


NMDS_Cytb$Lake <-mapvalues(NMDS_Cytb$Lake, c("Betley_Mere","Chapel_Mere","Fenemere","Kenfig_Pool","Llan_Bwch-llyn","Llangorse_Lake",  
                                             "Llyn_Cwellyn","Llyn_Ogwen","Llyn_Padarn","Llyn_Penrhyn","Llyn_Traffwll","Maer_Pool",       
                                             "Oss_Mere","Watch_Lane_Flash"), 
                                            c("BET","CAM","FEN","KEN","LLB","LLG","CWE","OGW","PAD","PEN","TRA","MAP","OSS","WLF"))
    


NMDS_Cytb_wide <- reshape(NMDS_Cytb, idvar= c("SampleID","Lake","Transect"), timevar= "Species", direction = "wide")

names(NMDS_Cytb_wide)

NMDS_Cytb_wide <- rename (NMDS_Cytb_wide,c( "Reads.EEL"="EEL","Reads.BRE"="BRE","Reads.CRU"="CRU",
                                            "Reads.CAR"="CAR","Reads.SUN"="SUN","Reads.MIN"="MIN",
                                            "Reads.GOF"="GOF","Reads.BIT"="BIT","Reads.ROA"="ROA","Reads.RUD"="RUD",
                                            "Reads.TEN"="TEN","Reads.PIK"="PIK","Reads.BUL"="BUL","Reads.3SS"="3SS",
                                            "Reads.RTR"="RTR","Reads.SAL"="SAL","Reads.BTR"="BTR","Reads.CHA"="CHA","Reads.PER"="PER"))


NMDS_Cytb_wide.rownames <- data.frame(NMDS_Cytb_wide[,-1], row.names=NMDS_Cytb_wide[,1])


NMDS_Cytb_wide.rownames[is.na(NMDS_Cytb_wide.rownames)] <- 0



NMDS_Cytb_wide.active <- NMDS_Cytb_wide.rownames[, 3:ncol(NMDS_Cytb_wide.rownames)]




NMDS_Cytb_data <- metaMDS(NMDS_Cytb_wide.active, distance="bray", k=2, trymax=1500)


#check stressplot
stressplot(NMDS_Cytb_data)





#Using the scores function from vegan to extract the site scores and convert to a data.frame

NMDS_Cytb_data.scores <- as.data.frame(scores(NMDS_Cytb_data))

# create a column of replicate names, from the rownames of data.scores
NMDS_Cytb_data.scores$Rep <- rownames(NMDS_Cytb_data.scores)

#  add the Treatment variable created earlier
NMDS_Cytb_data.scores$Lake <- NMDS_Cytb_wide.rownames$Lake

NMDS_Cytb_data.scores$Transect <- NMDS_Cytb_wide.rownames$Transect


NMDS_Cytb_data.scores$Locus <- "Cytb"

#look at the data
head(NMDS_Cytb_data.scores)



#Using the scores function from vegan to extract the species scores and convert to a data.frame
NMDS_Cytb_species.scores <- as.data.frame(scores(NMDS_Cytb_data, "species"))

# create a column of species, from the rownames of species.scores
NMDS_Cytb_species.scores$Species <- rownames(NMDS_Cytb_species.scores)

NMDS_Cytb_species.scores$Locus <- "Cytb"

#look at the data
head(NMDS_Cytb_species.scores)



NMDS_Cytb_species.scores$Species <-mapvalues(NMDS_Cytb_species.scores$Species, c("X3SS"), 
                                            c("3SS"))





#Using the chull function from vegan to extract the hull values and convert to a data.frame

NMDS_Cytb.PAD <- NMDS_Cytb_data.scores[NMDS_Cytb_data.scores$Lake == "PAD", ][chull(NMDS_Cytb_data.scores[NMDS_Cytb_data.scores$Lake == 
                                                                                                        "PAD", c("NMDS1", "NMDS2")]), ]
NMDS_Cytb.CWE <- NMDS_Cytb_data.scores[NMDS_Cytb_data.scores$Lake == "CWE", ][chull(NMDS_Cytb_data.scores[NMDS_Cytb_data.scores$Lake == 
                                                                                                        "CWE", c("NMDS1", "NMDS2")]), ]

NMDS_Cytb.LLG <- NMDS_Cytb_data.scores[NMDS_Cytb_data.scores$Lake == "LLG", ][chull(NMDS_Cytb_data.scores[NMDS_Cytb_data.scores$Lake == 
                                                                                                        "LLG", c("NMDS1", "NMDS2")]), ]

NMDS_Cytb.TRA <- NMDS_Cytb_data.scores[NMDS_Cytb_data.scores$Lake == "TRA", ][chull(NMDS_Cytb_data.scores[NMDS_Cytb_data.scores$Lake == 
                                                                                                        "TRA", c("NMDS1", "NMDS2")]), ]

NMDS_Cytb.OGW <- NMDS_Cytb_data.scores[NMDS_Cytb_data.scores$Lake == "OGW", ][chull(NMDS_Cytb_data.scores[NMDS_Cytb_data.scores$Lake == 
                                                                                                        "OGW", c("NMDS1", "NMDS2")]), ]

NMDS_Cytb.KEN <- NMDS_Cytb_data.scores[NMDS_Cytb_data.scores$Lake == "KEN", ][chull(NMDS_Cytb_data.scores[NMDS_Cytb_data.scores$Lake == 
                                                                                                        "KEN", c("NMDS1", "NMDS2")]), ]

NMDS_Cytb.PEN <- NMDS_Cytb_data.scores[NMDS_Cytb_data.scores$Lake == "PEN", ][chull(NMDS_Cytb_data.scores[NMDS_Cytb_data.scores$Lake == 
                                                                                                        "PEN", c("NMDS1", "NMDS2")]), ]

NMDS_Cytb.LLB <- NMDS_Cytb_data.scores[NMDS_Cytb_data.scores$Lake == "LLB", ][chull(NMDS_Cytb_data.scores[NMDS_Cytb_data.scores$Lake == 
                                                                                                        "LLB", c("NMDS1", "NMDS2")]), ]

NMDS_Cytb.WLF <- NMDS_Cytb_data.scores[NMDS_Cytb_data.scores$Lake == "WLF", ][chull(NMDS_Cytb_data.scores[NMDS_Cytb_data.scores$Lake == 
                                                                                                        "WLF", c("NMDS1", "NMDS2")]), ]

NMDS_Cytb.OSS <- NMDS_Cytb_data.scores[NMDS_Cytb_data.scores$Lake == "OSS", ][chull(NMDS_Cytb_data.scores[NMDS_Cytb_data.scores$Lake == 
                                                                                                      "OSS", c("NMDS1", "NMDS2")]), ]

NMDS_Cytb.FEN <- NMDS_Cytb_data.scores[NMDS_Cytb_data.scores$Lake == "FEN", ][chull(NMDS_Cytb_data.scores[NMDS_Cytb_data.scores$Lake == 
                                                                                                      "FEN", c("NMDS1", "NMDS2")]), ]

NMDS_Cytb.MAP <- NMDS_Cytb_data.scores[NMDS_Cytb_data.scores$Lake == "MAP", ][chull(NMDS_Cytb_data.scores[NMDS_Cytb_data.scores$Lake == 
                                                                                                      "MAP", c("NMDS1", "NMDS2")]), ]

NMDS_Cytb.BET <- NMDS_Cytb_data.scores[NMDS_Cytb_data.scores$Lake == "BET", ][chull(NMDS_Cytb_data.scores[NMDS_Cytb_data.scores$Lake == 
                                                                                                      "BET", c("NMDS1", "NMDS2")]), ]

NMDS_Cytb.CAM <- NMDS_Cytb_data.scores[NMDS_Cytb_data.scores$Lake == "CAM", ][chull(NMDS_Cytb_data.scores[NMDS_Cytb_data.scores$Lake == 
                                                                                                      "CAM", c("NMDS1", "NMDS2")]), ]




NMDS_Cytb_group.scores <- NMDS_Cytb_data.scores


NMDS_Cytb_group.scores$Lake <-mapvalues(NMDS_Cytb_group.scores$Lake, c("CWE","OGW","PAD","TRA","PEN","KEN","MAP","CAM","LLB","LLG","OSS","FEN","WLF","BET"), 
                                       c("G1","G1","G1","G2","G2","G3","G3","G3","G3","G3","G3","G3","G4","G4"))


NMDS_Cytb.G1 <- NMDS_Cytb_group.scores[NMDS_Cytb_group.scores$Lake == "G1", ][chull(NMDS_Cytb_group.scores[NMDS_Cytb_group.scores$Lake == 
                                                                                                         "G1", c("NMDS1", "NMDS2")]), ]
NMDS_Cytb.G2 <- NMDS_Cytb_group.scores[NMDS_Cytb_group.scores$Lake == "G2", ][chull(NMDS_Cytb_group.scores[NMDS_Cytb_group.scores$Lake == 
                                                                                                         "G2", c("NMDS1", "NMDS2")]), ]
NMDS_Cytb.G3 <- NMDS_Cytb_group.scores[NMDS_Cytb_group.scores$Lake == "G3", ][chull(NMDS_Cytb_group.scores[NMDS_Cytb_group.scores$Lake == 
                                                                                                         "G3", c("NMDS1", "NMDS2")]), ]
NMDS_Cytb.G4 <- NMDS_Cytb_group.scores[NMDS_Cytb_group.scores$Lake == "G4", ][chull(NMDS_Cytb_group.scores[NMDS_Cytb_group.scores$Lake == 
                                                                                                         "G4", c("NMDS1", "NMDS2")]), ]

NMDS_Cytb_hull.group <- rbind(NMDS_Cytb.G1,NMDS_Cytb.G2,NMDS_Cytb.G3,NMDS_Cytb.G4)
NMDS_Cytb_hull.group$Locus <- "Cytb"

#combine sampling site data

NMDS_Cytb_hull.data <- rbind(NMDS_Cytb.PAD,NMDS_Cytb.CWE,NMDS_Cytb.LLG,NMDS_Cytb.TRA,NMDS_Cytb.OGW,NMDS_Cytb.KEN,NMDS_Cytb.PEN,
                            NMDS_Cytb.LLB,NMDS_Cytb.WLF,NMDS_Cytb.OSS,NMDS_Cytb.FEN,NMDS_Cytb.MAP,NMDS_Cytb.BET,NMDS_Cytb.CAM)  
NMDS_Cytb_hull.data$Locus <- "Cytb"




NMDS_hull.data <- rbind(NMDS_12S_hull.data,NMDS_Cytb_hull.data)

NMDS_species.scores <- rbind(NMDS_12S_species.scores,NMDS_Cytb_species.scores)

NMDS_data.scores <- rbind(NMDS_12S_data.scores,NMDS_Cytb_data.scores)

NMDS_hull.group <- rbind(NMDS_12S_hull.group,NMDS_Cytb_hull.group)




NMDS_12S_hull.data$Lake <- factor(NMDS_12S_hull.data$Lake, levels = (c("CWE","PAD","OGW","PEN","TRA","KEN","LLB","LLG","MAP","CAM","OSS","FEN","WLF","BET")))


NMDS_12S_data.scores$Lake <- factor(NMDS_12S_data.scores$Lake, levels = (c("CWE","PAD","OGW","PEN","TRA","KEN","LLB","LLG","MAP","CAM","OSS","FEN","WLF","BET")))

NMDS_12S_data.scores$Transect <- factor(NMDS_12S_data.scores$Transect, levels = (c("SL","OSL")))



NMDS_Cytb_hull.data$Lake <- factor(NMDS_Cytb_hull.data$Lake, levels = (c("CWE","PAD","OGW","PEN","TRA","KEN","LLB","LLG","MAP","CAM","OSS","FEN","WLF","BET")))


NMDS_Cytb_data.scores$Lake <- factor(NMDS_Cytb_data.scores$Lake, levels = (c("CWE","PAD","OGW","PEN","TRA","KEN","LLB","LLG","MAP","CAM","OSS","FEN","WLF","BET")))

NMDS_Cytb_data.scores$Transect <- factor(NMDS_Cytb_data.scores$Transect, levels = (c("SL","OSL")))

colpalette <-c("#8cd6b5",
               "#4f2d46",
               "#638d48",
               "#c48943",
               "#669bb7",
               "#b74ace",
               "#cb9f9c",
               "#cfd058",
               "#74d656",
               "#a688c9",
               "#c34980",
               "#5541a7",
               "#4b4d32",
               "#c34837")





###ellipse based on sampling group



Type_label_1 <- data.frame(x=c(-3.1,-0.7,0.5,1),y=c(-0.5,0.5,1.1,-1.5),Locus = c("Cytb"), labs = c("Com1","Com2","Com3","Com4"))

Type_label_2 <- data.frame(x=c(-3.1,-1.2,1.5,-1),y=c(-0.5,1,1,-1.2),Locus = c("sRNA"), labs = c("Com1","Com2","Com3","Com4"))





NMDS_Cytb_FIG <- ggplot() +
                stat_ellipse(data=NMDS_Cytb_hull.group,aes(x=NMDS1,y=NMDS2,group=Lake),type = "norm",linetype = "longdash", level = 0.3, size=0.5) + # add the convex hulls
                geom_point(data=NMDS_Cytb_data.scores,aes(x=NMDS1,y=NMDS2,colour=Lake,shape=Transect),size=3.5) + # add the point markers
                geom_text_repel(data=NMDS_Cytb_species.scores,aes(x=NMDS1,y=NMDS2,label=Species),colour="black",size=6,fontface = "italic") +  # add the species labels
                geom_text(data=Type_label_1, aes(x=x, y=y,label=labs), size=7,colour="red")+
                scale_colour_manual(name="Sampling lake",values=colpalette) +
                scale_shape_manual(name="Location",values = c(16,17)) +
                guides(col = guide_legend(ncol = 7, byrow = TRUE,order = 1), shape = guide_legend(ncol = 1,byrow = TRUE))+
                theme_bw()+theme(text=element_text(size=20),legend.position="bottom",legend.title = element_text(size=15),legend.text = element_text(size=15),
                                 panel.grid.minor=element_blank(),panel.grid.major =element_blank(),strip.text=element_blank())


# NMDS_Cytb_FIG <- ggplot() + 
#             stat_ellipse(data=NMDS_Cytb_hull.group,aes(x=NMDS1,y=NMDS2,group=Lake),type = "norm",linetype = "longdash", level = 0.3, size=0.5) + # add the convex hulls
#             geom_point(data=NMDS_Cytb_data.scores,aes(x=NMDS1,y=NMDS2,colour=Lake),size=3.5) + # add the point markers
#             geom_text_repel(data=NMDS_Cytb_species.scores,aes(x=NMDS1,y=NMDS2,label=Species),colour="black",size=6,fontface = "italic") +  # add the species labels
#             geom_text(data=Type_label_1, aes(x=x, y=y,label=labs), size=7,colour="red")+
#             scale_colour_manual(name="Sampling lake",values=colpalette) +
#             scale_shape_manual(name="Location",values = c(16,17)) +
#             guides(col = guide_legend(ncol = 7, byrow = TRUE,order = 1), shape = guide_legend(ncol = 1,byrow = TRUE))+
#             theme_bw()+theme(text=element_text(size=20),legend.position="bottom",legend.title = element_text(size=15),legend.text = element_text(size=15),
#                              panel.grid.minor=element_blank(),panel.grid.major =element_blank(),strip.text=element_blank())

legend_NMDS <- get_legend(NMDS_Cytb_FIG)

NMDS_Cytb_FIG <- NMDS_Cytb_FIG + theme(legend.position="none")


NMDS_12S_FIG <- ggplot() + 
              stat_ellipse(data=NMDS_12S_hull.group,aes(x=NMDS1,y=NMDS2,group=Lake),type = "norm",linetype = "longdash", level = 0.3, size=0.5) + # add the convex hulls
              geom_point(data=NMDS_12S_data.scores,aes(x=NMDS1,y=NMDS2,colour=Lake,shape=Transect),size=3.5) + # add the point markers
              geom_text_repel(data=NMDS_12S_species.scores,aes(x=NMDS1,y=NMDS2,label=Species),colour="black",size=6,fontface = "italic") +  # add the species labels
              geom_text(data=Type_label_2, aes(x=x, y=y,label=labs), size=7,colour="red")+
              scale_colour_manual(name="Sampling lake",values=colpalette) +
              scale_shape_manual(name="Location",values = c(16,17)) +
              guides(col = guide_legend(ncol = 7, byrow = TRUE,order = 1), shape = guide_legend(ncol = 1, byrow = TRUE))+
              theme_bw()+theme(text=element_text(size=20),legend.position="None",legend.title = element_text(size=15),legend.text = element_text(size=15),
                               panel.grid.minor=element_blank(),panel.grid.major =element_blank(),strip.text=element_blank())




 



###ANOSIM_12S####



NMDS_12S_wide.rownames <- data.frame(NMDS_12S_wide[,-1], row.names=NMDS_12S_wide[,1])


NMDS_12S_wide.rownames[is.na(NMDS_12S_wide.rownames)] <- 0

NMDS_12S_Group <- NMDS_12S_wide.rownames

NMDS_12S_Group$Lake <-mapvalues(NMDS_12S_Group$Lake, c("CWE","OGW","PAD","TRA","PEN","KEN","MAP","CAM","LLB","LLG","OSS","FEN","WLF","BET"), 
                                                     c("G1","G1","G1","G2","G2","G3","G3","G3","G3","G3","G3","G3","G4","G4"))






#Global
NMDS_12S.dist <- vegdist(NMDS_12S_wide.active,method="bray")


sRNA.env <- droplevels(NMDS_12S_Group$Lake)

sRNA.env <- as.data.frame(sRNA.env)

sRNA.env<- rename (sRNA.env,c("sRNA.env"="Lake"))



sRNA.anosim<- anosim(NMDS_12S.dist, sRNA.env$Lake)
summary(sRNA.anosim)
plot(sRNA.anosim)





NMDS_12S_Group_sub <- NMDS_12S_Group[which(NMDS_12S_Group$Lake == "G2"
                                        | NMDS_12S_Group$Lake == "G3"),]



NMDS_12S_wide.active_sub <- NMDS_12S_Group_sub [, 3:ncol(NMDS_12S_Group_sub)]




NMDS_12S.dist_sub <- vegdist(NMDS_12S_wide.active_sub,method="bray")


sRNA.env <- droplevels(NMDS_12S_Group_sub$Lake)

sRNA.env <- as.data.frame(sRNA.env)

sRNA.env<- rename (sRNA.env,c("sRNA.env"="Lake"))



sRNA.anosim<- anosim(NMDS_12S.dist_sub, sRNA.env$Lake)
summary(sRNA.anosim)
plot(sRNA.anosim)





NMDS_12S_wide.rownames_G1 <- NMDS_12S_wide.rownames[which(NMDS_12S_wide.rownames$Lake == "CWE"
                                                        | NMDS_12S_wide.rownames$Lake == "OGW"
                                                        | NMDS_12S_wide.rownames$Lake == "PAD"),]



NMDS_12S_wide.active_G1 <- NMDS_12S_wide.rownames_G1 [, 3:ncol(NMDS_12S_wide.rownames_G1)]




NMDS_12S.dist_G1 <- vegdist(NMDS_12S_wide.active_G1,method="bray")


sRNA.env <- droplevels(NMDS_12S_wide.rownames_G1$Lake)

sRNA.env <- as.data.frame(sRNA.env)

sRNA.env<- rename (sRNA.env,c("sRNA.env"="Lake"))



sRNA.anosim<- anosim(NMDS_12S.dist_G1, sRNA.env$Lake)
summary(sRNA.anosim)
plot(sRNA.anosim)






NMDS_12S_wide.rownames_G2 <- NMDS_12S_wide.rownames[which(NMDS_12S_wide.rownames$Lake == "TRA"
                                                          | NMDS_12S_wide.rownames$Lake == "PEN"),]



NMDS_12S_wide.active_G2 <- NMDS_12S_wide.rownames_G2 [, 3:ncol(NMDS_12S_wide.rownames_G2)]




NMDS_12S.dist_G2 <- vegdist(NMDS_12S_wide.active_G2,method="bray")


sRNA.env <- droplevels(NMDS_12S_wide.rownames_G2$Lake)

sRNA.env <- as.data.frame(sRNA.env)

sRNA.env<- rename (sRNA.env,c("sRNA.env"="Lake"))



sRNA.anosim<- anosim(NMDS_12S.dist_G2, sRNA.env$Lake)
summary(sRNA.anosim)
plot(sRNA.anosim)







NMDS_12S_wide.rownames_G3 <- NMDS_12S_wide.rownames[which(NMDS_12S_wide.rownames$Lake == "KEN"
                                                          | NMDS_12S_wide.rownames$Lake == "MAP"
                                                          | NMDS_12S_wide.rownames$Lake == "CAM"
                                                          | NMDS_12S_wide.rownames$Lake == "LLB"
                                                          | NMDS_12S_wide.rownames$Lake == "LLG"
                                                          | NMDS_12S_wide.rownames$Lake == "OSS"
                                                          | NMDS_12S_wide.rownames$Lake == "FEN"),]




NMDS_12S_wide.active_G3 <- NMDS_12S_wide.rownames_G3 [, 3:ncol(NMDS_12S_wide.rownames_G3)]




NMDS_12S.dist_G3 <- vegdist(NMDS_12S_wide.active_G3,method="bray")


sRNA.env <- droplevels(NMDS_12S_wide.rownames_G3$Lake)

sRNA.env <- as.data.frame(sRNA.env)

sRNA.env<- rename (sRNA.env,c("sRNA.env"="Lake"))



sRNA.anosim<- anosim(NMDS_12S.dist_G3, sRNA.env$Lake)
summary(sRNA.anosim)
plot(sRNA.anosim)





NMDS_12S_wide.rownames_G4 <- NMDS_12S_wide.rownames[which(NMDS_12S_wide.rownames$Lake == "WLF"
                                                          | NMDS_12S_wide.rownames$Lake == "BET"),]




NMDS_12S_wide.active_G4 <- NMDS_12S_wide.rownames_G4 [, 3:ncol(NMDS_12S_wide.rownames_G4)]




NMDS_12S.dist_G4 <- vegdist(NMDS_12S_wide.active_G4,method="bray")


sRNA.env <- droplevels(NMDS_12S_wide.rownames_G4$Lake)

sRNA.env <- as.data.frame(sRNA.env)

sRNA.env<- rename (sRNA.env,c("sRNA.env"="Lake"))



sRNA.anosim<- anosim(NMDS_12S.dist_G4, sRNA.env$Lake)
summary(sRNA.anosim)
plot(sRNA.anosim)




###ANOSIM_Cytb####



NMDS_Cytb_wide.rownames <- data.frame(NMDS_Cytb_wide[,-1], row.names=NMDS_Cytb_wide[,1])


NMDS_Cytb_wide.rownames[is.na(NMDS_Cytb_wide.rownames)] <- 0

NMDS_Cytb_Group <- NMDS_Cytb_wide.rownames

NMDS_Cytb_Group$Lake <-mapvalues(NMDS_Cytb_Group$Lake, c("CWE","OGW","PAD","TRA","PEN","KEN","MAP","CAM","LLB","LLG","OSS","FEN","WLF","BET"), 
                                c("G1","G1","G1","G2","G2","G3","G3","G3","G3","G3","G3","G3","G4","G4"))






#Global
NMDS_Cytb.dist <- vegdist(NMDS_Cytb_wide.active,method="bray")


Cytb.env <- droplevels(NMDS_Cytb_Group$Lake)

Cytb.env <- as.data.frame(Cytb.env)

Cytb.env<- rename (Cytb.env,c("Cytb.env"="Lake"))



Cytb.anosim<- anosim(NMDS_Cytb.dist, Cytb.env$Lake)
summary(Cytb.anosim)
plot(Cytb.anosim)





NMDS_Cytb_Group_sub <- NMDS_Cytb_Group[which(NMDS_Cytb_Group$Lake == "G2"
                                           | NMDS_Cytb_Group$Lake == "G3"),]



NMDS_Cytb_wide.active_sub <- NMDS_Cytb_Group_sub [, 3:ncol(NMDS_Cytb_Group_sub)]




NMDS_Cytb.dist_sub <- vegdist(NMDS_Cytb_wide.active_sub,method="bray")


Cytb.env <- droplevels(NMDS_Cytb_Group_sub$Lake)

Cytb.env <- as.data.frame(Cytb.env)

Cytb.env<- rename (Cytb.env,c("Cytb.env"="Lake"))



Cytb.anosim<- anosim(NMDS_Cytb.dist_sub, Cytb.env$Lake)
summary(Cytb.anosim)
plot(Cytb.anosim)





NMDS_Cytb_wide.rownames_G1 <- NMDS_Cytb_wide.rownames[which(NMDS_Cytb_wide.rownames$Lake == "CWE"
                                                          | NMDS_Cytb_wide.rownames$Lake == "OGW"
                                                          | NMDS_Cytb_wide.rownames$Lake == "PAD"),]



NMDS_Cytb_wide.active_G1 <- NMDS_Cytb_wide.rownames_G1 [, 3:ncol(NMDS_Cytb_wide.rownames_G1)]




NMDS_Cytb.dist_G1 <- vegdist(NMDS_Cytb_wide.active_G1,method="bray")


Cytb.env <- droplevels(NMDS_Cytb_wide.rownames_G1$Lake)

Cytb.env <- as.data.frame(Cytb.env)

Cytb.env<- rename (Cytb.env,c("Cytb.env"="Lake"))



Cytb.anosim<- anosim(NMDS_Cytb.dist_G1, Cytb.env$Lake)
summary(Cytb.anosim)
plot(Cytb.anosim)






NMDS_Cytb_wide.rownames_G2 <- NMDS_Cytb_wide.rownames[which(NMDS_Cytb_wide.rownames$Lake == "TRA"
                                                          | NMDS_Cytb_wide.rownames$Lake == "PEN"),]



NMDS_Cytb_wide.active_G2 <- NMDS_Cytb_wide.rownames_G2 [, 3:ncol(NMDS_Cytb_wide.rownames_G2)]




NMDS_Cytb.dist_G2 <- vegdist(NMDS_Cytb_wide.active_G2,method="bray")


Cytb.env <- droplevels(NMDS_Cytb_wide.rownames_G2$Lake)

Cytb.env <- as.data.frame(Cytb.env)

Cytb.env<- rename (Cytb.env,c("Cytb.env"="Lake"))



Cytb.anosim<- anosim(NMDS_Cytb.dist_G2, Cytb.env$Lake)
summary(Cytb.anosim)
plot(Cytb.anosim)







NMDS_Cytb_wide.rownames_G3 <- NMDS_Cytb_wide.rownames[which(NMDS_Cytb_wide.rownames$Lake == "KEN"
                                                          | NMDS_Cytb_wide.rownames$Lake == "MAP"
                                                          | NMDS_Cytb_wide.rownames$Lake == "CAM"
                                                          | NMDS_Cytb_wide.rownames$Lake == "LLB"
                                                          | NMDS_Cytb_wide.rownames$Lake == "LLG"
                                                          | NMDS_Cytb_wide.rownames$Lake == "OSS"
                                                          | NMDS_Cytb_wide.rownames$Lake == "FEN"),]




NMDS_Cytb_wide.active_G3 <- NMDS_Cytb_wide.rownames_G3 [, 3:ncol(NMDS_Cytb_wide.rownames_G3)]




NMDS_Cytb.dist_G3 <- vegdist(NMDS_Cytb_wide.active_G3,method="bray")


Cytb.env <- droplevels(NMDS_Cytb_wide.rownames_G3$Lake)

Cytb.env <- as.data.frame(Cytb.env)

Cytb.env<- rename (Cytb.env,c("Cytb.env"="Lake"))



Cytb.anosim<- anosim(NMDS_Cytb.dist_G3, Cytb.env$Lake)
summary(Cytb.anosim)
plot(Cytb.anosim)





NMDS_Cytb_wide.rownames_G4 <- NMDS_Cytb_wide.rownames[which(NMDS_Cytb_wide.rownames$Lake == "WLF"
                                                            | NMDS_Cytb_wide.rownames$Lake == "BET"),]




NMDS_Cytb_wide.active_G4 <- NMDS_Cytb_wide.rownames_G4 [, 3:ncol(NMDS_Cytb_wide.rownames_G4)]




NMDS_Cytb.dist_G4 <- vegdist(NMDS_Cytb_wide.active_G4,method="bray")


Cytb.env <- droplevels(NMDS_Cytb_wide.rownames_G4$Lake)

Cytb.env <- as.data.frame(Cytb.env)

Cytb.env<- rename (Cytb.env,c("Cytb.env"="Lake"))



Cytb.anosim<- anosim(NMDS_Cytb.dist_G4, Cytb.env$Lake)
summary(Cytb.anosim)
plot(Cytb.anosim)




###Welsh lakes####


LP_rs_wlake <- LP_rs[which(LP_rs$Project == 'Lake'),]
LP_rs_wlake_SO <- LP_rs_wlake[which(LP_rs_wlake$Reads >0),]



LP_rs_wlake_SO2 = data.frame(Lake=factor(LP_rs_wlake_SO$Lake),Project=factor(LP_rs_wlake_SO$Project),Total_n=LP_rs_wlake_SO$Total_n,
                             Locus=factor(LP_rs_wlake_SO$Locus),Species=factor(LP_rs_wlake_SO$Species))



#summary the read counts based on species

LP_rs_wlake_SO_sum <- aggregate(Reads~Lake+Total_n+Project+Locus+Species, data = LP_rs_wlake_SO, FUN=sum)



LP_rs_wlake_SO_sum2 <- merge(LP_rs_wlake_SO_sum,LP_Sample_Total,by=c("Lake","Total_n","Project","Locus"))


#RRC relative read count 

#LP_rs_wlake_SO_sum2$RRC <- LP_rs_wlake_SO_sum2$Reads/LP_rs_wlake_SO_sum2$SUM


#LP_rs_wlake_SO_sum2_RRC <- LP_rs_wlake_SO_sum2[c(-6,-7)]


#RRC2 relative read count_proportion of this species read counts of the samples with this species

LP_rs_wlake_SO_SUM <- aggregate(SUM~Lake+Total_n+Project+Locus+Species, data = LP_rs_wlake_SO, FUN=sum)


LP_rs_wlake_SO_sum3 <- merge(LP_rs_wlake_SO_sum,LP_rs_wlake_SO_SUM,by=c("Lake","Total_n","Project","Locus","Species")) #total reads of the samples with this species

LP_rs_wlake_SO_sum3$RRC2 <- LP_rs_wlake_SO_sum3$Reads/LP_rs_wlake_SO_sum3$SUM




##Reshape the data based on Loci


LP_rs_wlake_loci =dcast(LP_rs_wlake_SO2, Lake+Total_n+Project+Species~Locus)

LP_rs_wlake_loci$Cytb_SO <- LP_rs_wlake_loci$Cytb/LP_rs_wlake_loci$Total_n

LP_rs_wlake_loci$"sRNA_SO" <- LP_rs_wlake_loci$`12S`/LP_rs_wlake_loci$Total_n

LP_rs_wlake_loci <- rename (LP_rs_wlake_loci,c("Cytb"="Cytb_N","12S"="sRNA_N"))


#LP_rs_wlake_SO_sum_loci <-dcast(LP_rs_wlake_SO_sum2_RRC, Lake+Total_n+Project+Species~Locus)

#LP_rs_wlake_SO_sum_loci <- rename (LP_rs_wlake_SO_sum_loci,c("Cytb"="Cytb_RRC","12S"="sRNA_RRC"))



LP_rs_wlake_SO_sum_loci2 <-dcast(LP_rs_wlake_SO_sum3, Lake+Total_n+Project+Species~Locus)

LP_rs_wlake_SO_sum_loci2 <- rename (LP_rs_wlake_SO_sum_loci2,c("Cytb"="Cytb_RRC2","12S"="sRNA_RRC2"))



LP_rs_wlake_SO_sum_loci3 <-dcast(LP_rs_wlake_SO_sum, Lake+Total_n+Project+Species~Locus)

LP_rs_wlake_SO_sum_loci3 <- rename (LP_rs_wlake_SO_sum_loci3,c("Cytb"="Cytb_RC","12S"="sRNA_RC"))


#merge SO and RRC
#Merge_wlake_loci <- merge(LP_rs_wlake_loci, LP_rs_wlake_SO_sum_loci, by=c("Lake","Total_n","Project","Species"))


#merge SO and RRC wuth RRC2 
Merge_wlake_loci <- merge(LP_rs_wlake_loci, LP_rs_wlake_SO_sum_loci2, by=c("Lake","Total_n","Project","Species"))

Merge_wlake_loci <- merge(Merge_wlake_loci, LP_rs_wlake_SO_sum_loci3, by=c("Lake","Total_n","Project","Species"))

Merge_wlake_loci <- Merge_wlake_loci[which(Merge_wlake_loci$Species != "Salmo"),] 

Merge_wlake_loci [is.na(Merge_wlake_loci )] <- 0

#rm two variables "Cytb_N","sRNA_N"
Merge_wlake_loci$Cytb_N <- Merge_wlake_loci$sRNA_N <- NULL


levels(droplevels(Merge_wlake_loci$Species))

#MAX_SO

for (i in 1:nrow(Merge_wlake_loci)) {
  if (Merge_wlake_loci$Cytb_SO[i] > Merge_wlake_loci$sRNA_SO[i]) {Merge_wlake_loci$MAX_SO[i] <- Merge_wlake_loci$Cytb_SO[i]} 
  else {Merge_wlake_loci$MAX_SO[i] <- Merge_wlake_loci$sRNA_SO[i]}
}


#MAX_RRC2
for (i in 1:nrow(Merge_wlake_loci)) {
  if (Merge_wlake_loci$Cytb_RRC2[i] > Merge_wlake_loci$sRNA_RRC2[i]) {Merge_wlake_loci$MAX_RRC2[i] <- Merge_wlake_loci$Cytb_RRC2[i]} 
  else {Merge_wlake_loci$MAX_RRC2[i] <- Merge_wlake_loci$sRNA_RRC2[i]}
}

#Number of loci

for (i in 1:nrow(Merge_wlake_loci)) {
  if (Merge_wlake_loci$Cytb_RRC[i] > 0 & Merge_wlake_loci$sRNA_RRC[i] > 0) {Merge_wlake_loci$locus_n[i] <- 2} 
  else {Merge_wlake_loci$locus_n[i] <- 1}
}




Tench <- read.csv (file="Lake project/Tench.csv", header = TRUE)



Merge_wlake_loci <- rbind (Merge_wlake_loci,Tench)

#Expeted DARFOR scale


Wlake_Expected <- read.csv(file="Lake project/Welsh_Lake_expected_DAFOR_Presence.csv", header = TRUE)

Wlake_Expected_EDNA <- merge(Merge_wlake_loci, Wlake_Expected, by=c("Lake","Species"))


#Site occupancy score

for (i in 1:nrow(Wlake_Expected_EDNA)) {
  if (Wlake_Expected_EDNA$MAX_SO[i] > 0.8) {Wlake_Expected_EDNA$DAFOR[i] <- 5} 
  else if (Wlake_Expected_EDNA$MAX_SO[i] > 0.6) {Wlake_Expected_EDNA$DAFOR[i] <- 4} 
  else if (Wlake_Expected_EDNA$MAX_SO[i] > 0.3) {Wlake_Expected_EDNA$DAFOR[i] <- 3} 
  else if (Wlake_Expected_EDNA$MAX_SO[i] > 0.1) {Wlake_Expected_EDNA$DAFOR[i] <- 2} 
  else if (Wlake_Expected_EDNA$MAX_SO[i] > 0) {Wlake_Expected_EDNA$DAFOR[i] <- 1} 
  else {Wlake_Expected_EDNA$DAFOR[i] <- 0} 
}



#Confidence score

for (i in 1:nrow(Wlake_Expected_EDNA)) {
  if (Wlake_Expected_EDNA$locus_n[i] ==1 & Wlake_Expected_EDNA$DAFOR[i] == 1 & Wlake_Expected_EDNA$MAX_RRC[i] > 0.5) {Wlake_Expected_EDNA$Presence[i] <- 5}
  else if (Wlake_Expected_EDNA$locus_n[i] ==1 & Wlake_Expected_EDNA$DAFOR[i] > 1 & Wlake_Expected_EDNA$MAX_RRC[i] > 0.1) {Wlake_Expected_EDNA$Presence[i] <- 5}
  else if (Wlake_Expected_EDNA$locus_n[i] ==2 & Wlake_Expected_EDNA$DAFOR[i] == 1 & Wlake_Expected_EDNA$MAX_RRC[i] > 0.1) {Wlake_Expected_EDNA$Presence[i] <- 5}
  else if (Wlake_Expected_EDNA$locus_n[i] ==2 & Wlake_Expected_EDNA$DAFOR[i] > 1 & Wlake_Expected_EDNA$MAX_RRC[i] > 0.01) {Wlake_Expected_EDNA$Presence[i] <- 5}
  else if (Wlake_Expected_EDNA$locus_n[i] ==1 & Wlake_Expected_EDNA$DAFOR[i] == 1 & Wlake_Expected_EDNA$MAX_RRC[i] > 0.1) {Wlake_Expected_EDNA$Presence[i] <- 3}
  else if (Wlake_Expected_EDNA$locus_n[i] ==1 & Wlake_Expected_EDNA$DAFOR[i] >= 1 & Wlake_Expected_EDNA$MAX_RRC[i] > 0.05) {Wlake_Expected_EDNA$Presence[i] <- 3}
  else if (Wlake_Expected_EDNA$locus_n[i] ==2 & Wlake_Expected_EDNA$DAFOR[i] == 1 & Wlake_Expected_EDNA$MAX_RRC[i] > 0.01) {Wlake_Expected_EDNA$Presence[i] <- 3}
  else if (Wlake_Expected_EDNA$locus_n[i] ==2 & Wlake_Expected_EDNA$DAFOR[i] >= 1 & Wlake_Expected_EDNA$MAX_RRC[i] > 0.005) {Wlake_Expected_EDNA$Presence[i] <- 3}
  else if (Wlake_Expected_EDNA$locus_n[i] ==1 & Wlake_Expected_EDNA$DAFOR[i] >= 1 & Wlake_Expected_EDNA$MAX_RRC[i] > 0.005) {Wlake_Expected_EDNA$Presence[i] <- 1}
  else if (Wlake_Expected_EDNA$DAFOR[i] > 1 & Wlake_Expected_EDNA$MAX_RRC[i] < 0.005) {Wlake_Expected_EDNA$Presence[i] <- 1}
  else if (Wlake_Expected_EDNA$DAFOR[i] == 1 & Wlake_Expected_EDNA$MAX_RRC[i] < 0.005) {Wlake_Expected_EDNA$Presence[i] <- -1}
  else {Wlake_Expected_EDNA$Presence[i] <- 0} 
}



#Abundance score, corrected score (i.e. confidence score × site occupancy score)

for (i in 1:nrow(Wlake_Expected_EDNA)) {
  if (Wlake_Expected_EDNA$Presence[i]*Wlake_Expected_EDNA$DAFOR[i]>=25) {Wlake_Expected_EDNA$DAFOR_c[i] <- 5}
  else if (Wlake_Expected_EDNA$Presence[i]*Wlake_Expected_EDNA$DAFOR[i]>=15) {Wlake_Expected_EDNA$DAFOR_c[i] <- 4}
  else if (Wlake_Expected_EDNA$Presence[i]*Wlake_Expected_EDNA$DAFOR[i]>=9) {Wlake_Expected_EDNA$DAFOR_c[i] <- 3}
  else if (Wlake_Expected_EDNA$Presence[i]*Wlake_Expected_EDNA$DAFOR[i]>=4) {Wlake_Expected_EDNA$DAFOR_c[i] <- 2}
  else if (Wlake_Expected_EDNA$Presence[i]*Wlake_Expected_EDNA$DAFOR[i]>=1) {Wlake_Expected_EDNA$DAFOR_c[i] <- 1}
  else if (Wlake_Expected_EDNA$Presence[i]*Wlake_Expected_EDNA$DAFOR[i]==0) {Wlake_Expected_EDNA$DAFOR_c[i] <- 0}
  else {Wlake_Expected_EDNA$DAFOR_c[i] <- 1} 
}




##exclude unlikely present species




Wlake_excl_absent <- Wlake_Expected_EDNA[which(Wlake_Expected_EDNA$Presence != -1),]

#reorder based on the lake
Wlake_excl_absent$Lake <- factor(Wlake_excl_absent$Lake, levels = (c("Llyn_Padarn","Llyn_Cwellyn","Llangorse_Lake",
                                                                     "Llyn_Traffwll","Llyn_Ogwen","Kenfig_Pool","Llyn_Penrhyn","Llan_Bwch-llyn")))


levels(droplevels(Wlake_excl_absent$Species))

Wlake_excl_absent$Species <-mapvalues(Wlake_excl_absent$Species, c("Anguilla_anguilla","Abramis_brama","Carassius_auratus","Phoxinus_phoxinus", 
                                                                   "Rutilus_rutilus","Scardinius_erythrophthalmus","Tinca_tinca","Esox_lucius",            
                                                                   "Cottus_gobio","Gasterosteus_aculeatus","Perca_fluviatilis","Oncorhynchus_mykiss",
                                                                   "Salmo_salar","Salmo_trutta","Salvelinus_alpinus","Alburnus_alburnus", 
                                                                   "Gobio_gobio","Leuciscus_idus","Pungitius_pungitius"), 
                                                                c("EEL","BRE","GOF","MIN",
                                                                  "ROA","RUD","TEN","PIK",
                                                                  "BUL","3SS","PER","RTR",
                                                                  "SAL","BTR","CHA","BLE",
                                                                  "GUD","IDE","9SS"))

#presence confidence plot

ggplot(Wlake_excl_absent, aes(x = Ex_Presence, y = Presence))+
  geom_point(stat="identity",size=3)+
  facet_wrap(~Lake,ncol =4)+
  scale_y_continuous(limits = c(-1, 6), breaks=seq(-1, 6, 1))+
  scale_x_continuous(limits = c(-1, 6), breaks=seq(-1, 6, 1))+
  labs(x="Ex_Presence", y="Presence")+theme_bw()+
  theme(text=element_text(size=15))+
  geom_text_repel(aes(label=Species),size=3)



ggplot(Wlake_excl_absent, aes(x = Ex_Presence, y = Presence))+
  geom_point(stat="identity",size=3)+
  geom_smooth(formula= y ~ x,method="lm", se=FALSE,lwd=0.5,na.rm=TRUE, colour='black')+
  scale_y_continuous(limits = c(-1, 6), breaks=seq(-1, 6, 1))+
  scale_x_continuous(limits = c(-1, 6), breaks=seq(-1, 6, 1))+
  labs(x="Ex_Presence", y="Presence")+theme_bw()+
  theme(text=element_text(size=15))+
  geom_text_repel(aes(label=Species),size=3)

ggplot(Wlake_excl_absent, aes(x = Ex_Presence, y = Presence))+
  geom_point(stat="identity",size=3)+
  scale_y_continuous(limits = c(-1, 6), breaks=seq(-1, 6, 1))+
  scale_x_continuous(limits = c(-1, 6), breaks=seq(-1, 6, 1))+
  geom_smooth(formula= y ~ x,method="lm", se=FALSE,lwd=0.5,na.rm=TRUE, colour='black')+
  labs(x="Ex_Presence", y="Presence")+theme_bw()+
  theme(text=element_text(size=15))


##spearman_test
rhofun<-function(x, y) {
  corr=(cor.test(x, y,
                 alternative="two.sided", method="spearman",exact = FALSE))}


rhofun(Wlake_excl_absent$Ex_Presence,Wlake_excl_absent$Presence)$statistic

rhofun(Wlake_excl_absent$Ex_Presence,Wlake_excl_absent$Presence)$p.value

rhofun(Wlake_excl_absent$Ex_Presence,Wlake_excl_absent$Presence)$estimate


Wlake_excl_absent$Lake <-mapvalues(Wlake_excl_absent$Lake, c("Kenfig_Pool","Llan_Bwch-llyn","Llangorse_Lake",  
                                                 "Llyn_Cwellyn","Llyn_Ogwen","Llyn_Padarn","Llyn_Penrhyn","Llyn_Traffwll"), 
                                                 c("KEN","LLB","LLG","CWE","OGW","PAD","PEN","TRA"))












####Mere####

LP_rs_mere <- LP_rs[which(LP_rs$Project == 'Mere'),]
LP_rs_mere_SO <- LP_rs_mere[which(LP_rs_mere$Reads >0),]

LP_rs_mere_SO2 = data.frame(Lake=factor(LP_rs_mere_SO$Lake),Project=factor(LP_rs_mere_SO$Project),Total_n=LP_rs_mere_SO$Total_n,
                            Locus=factor(LP_rs_mere_SO$Locus),Species=factor(LP_rs_mere_SO$Species))

#summary the read counts based on species

LP_rs_mere_SO_sum <- aggregate(Reads~Lake+Total_n+Project+Locus+Species, data = LP_rs_mere_SO, FUN=sum)



#LP_rs_mere_SO_sum2 <- merge(LP_rs_mere_SO_sum,LP_Sample_Total,by=c("Lake","Total_n","Project","Locus"))


#RRC relative read count 

#LP_rs_mere_SO_sum2$RRC <- LP_rs_mere_SO_sum2$Reads/LP_rs_mere_SO_sum2$SUM


#LP_rs_mere_SO_sum2_RRC <- LP_rs_mere_SO_sum2[c(-6,-7)]


#RRC2 relative read count_proportion of this species read counts of the samples with this species

LP_rs_mere_SO_SUM <- aggregate(SUM~Lake+Total_n+Project+Locus+Species, data = LP_rs_mere_SO, FUN=sum)


LP_rs_mere_SO_sum3 <- merge(LP_rs_mere_SO_sum,LP_rs_mere_SO_SUM,by=c("Lake","Total_n","Project","Locus","Species")) #total reads of the samples with this species

LP_rs_mere_SO_sum3$RRC2 <- LP_rs_mere_SO_sum3$Reads/LP_rs_mere_SO_sum3$SUM




##Reshape the data based on Loci


LP_rs_mere_loci =dcast(LP_rs_mere_SO2, Lake+Total_n+Project+Species~Locus)

LP_rs_mere_loci$Cytb_SO <- LP_rs_mere_loci$Cytb/LP_rs_mere_loci$Total_n

LP_rs_mere_loci$"sRNA_SO" <- LP_rs_mere_loci$`12S`/LP_rs_mere_loci$Total_n

LP_rs_mere_loci <- rename (LP_rs_mere_loci,c("Cytb"="Cytb_N","12S"="sRNA_N"))


#LP_rs_mere_SO_sum_loci <-dcast(LP_rs_mere_SO_sum2_RRC, Lake+Total_n+Project+Species~Locus)

#LP_rs_mere_SO_sum_loci <- rename (LP_rs_mere_SO_sum_loci,c("Cytb"="Cytb_RRC","12S"="sRNA_RRC"))



LP_rs_mere_SO_sum_loci2 <-dcast(LP_rs_mere_SO_sum3, Lake+Total_n+Project+Species~Locus)

LP_rs_mere_SO_sum_loci2 <- rename (LP_rs_mere_SO_sum_loci2,c("Cytb"="Cytb_RRC2","12S"="sRNA_RRC2"))


LP_rs_mere_SO_sum_loci3 <-dcast(LP_rs_mere_SO_sum, Lake+Total_n+Project+Species~Locus)

LP_rs_mere_SO_sum_loci3 <- rename (LP_rs_mere_SO_sum_loci3,c("Cytb"="Cytb_RC","12S"="sRNA_RC"))


#merge SO and RRC
#Merge_mere_loci <- merge(LP_rs_mere_loci, LP_rs_mere_SO_sum_loci, by=c("Lake","Total_n","Project","Species"))


#merge SO and RRC wuth RRC2 
Merge_mere_loci <- merge(LP_rs_mere_loci, LP_rs_mere_SO_sum_loci2, by=c("Lake","Total_n","Project","Species"))


Merge_mere_loci <- merge(Merge_mere_loci, LP_rs_mere_SO_sum_loci3, by=c("Lake","Total_n","Project","Species"))


Merge_mere_loci [is.na(Merge_mere_loci )] <- 0

#rm two variables "Cytb_N","sRNA_N"
Merge_mere_loci$Cytb_N <- Merge_mere_loci$sRNA_N <- NULL


levels(droplevels(Merge_mere_loci$Species))

#MAX_SO

for (i in 1:nrow(Merge_mere_loci)) {
  if (Merge_mere_loci$Cytb_SO[i] > Merge_mere_loci$sRNA_SO[i]) {Merge_mere_loci$MAX_SO[i] <- Merge_mere_loci$Cytb_SO[i]} 
  else {Merge_mere_loci$MAX_SO[i] <- Merge_mere_loci$sRNA_SO[i]}
}


#MAX_RRC2
for (i in 1:nrow(Merge_mere_loci)) {
  if (Merge_mere_loci$Cytb_RRC2[i] > Merge_mere_loci$sRNA_RRC2[i]) {Merge_mere_loci$MAX_RRC2[i] <- Merge_mere_loci$Cytb_RRC2[i]} 
  else {Merge_mere_loci$MAX_RRC2[i] <- Merge_mere_loci$sRNA_RRC2[i]}
}

#Number of loci

for (i in 1:nrow(Merge_mere_loci)) {
  if (Merge_mere_loci$Cytb_RRC[i] > 0 & Merge_mere_loci$sRNA_RRC[i] > 0) {Merge_mere_loci$locus_n[i] <- 2} 
  else {Merge_mere_loci$locus_n[i] <- 1}
}






#Expeted DARFOR scale


mere_Expected <- read.csv(file="Lake project/Mere_Lake_expected_Bio_Num_Presence.csv", header = TRUE)

mere_Expected_EDNA <- merge(Merge_mere_loci, mere_Expected, by=c("Lake","Species"))


#Site occupancy score

for (i in 1:nrow(mere_Expected_EDNA)) {
  if (mere_Expected_EDNA$MAX_SO[i] > 0.8) {mere_Expected_EDNA$DAFOR[i] <- 5} 
  else if (mere_Expected_EDNA$MAX_SO[i] > 0.6) {mere_Expected_EDNA$DAFOR[i] <- 4} 
  else if (mere_Expected_EDNA$MAX_SO[i] > 0.3) {mere_Expected_EDNA$DAFOR[i] <- 3} 
  else if (mere_Expected_EDNA$MAX_SO[i] > 0.1) {mere_Expected_EDNA$DAFOR[i] <- 2} 
  else if (mere_Expected_EDNA$MAX_SO[i] > 0) {mere_Expected_EDNA$DAFOR[i] <- 1} 
  else {mere_Expected_EDNA$DAFOR[i] <- 0} 
}



#Confidence score

for (i in 1:nrow(mere_Expected_EDNA)) {
  if (mere_Expected_EDNA$locus_n[i] ==1 & mere_Expected_EDNA$DAFOR[i] == 1 & mere_Expected_EDNA$MAX_RRC[i] > 0.5) {mere_Expected_EDNA$Presence[i] <- 5}
  else if (mere_Expected_EDNA$locus_n[i] ==1 & mere_Expected_EDNA$DAFOR[i] > 1 & mere_Expected_EDNA$MAX_RRC[i] > 0.1) {mere_Expected_EDNA$Presence[i] <- 5}
  else if (mere_Expected_EDNA$locus_n[i] ==2 & mere_Expected_EDNA$DAFOR[i] == 1 & mere_Expected_EDNA$MAX_RRC[i] > 0.1) {mere_Expected_EDNA$Presence[i] <- 5}
  else if (mere_Expected_EDNA$locus_n[i] ==2 & mere_Expected_EDNA$DAFOR[i] > 1 & mere_Expected_EDNA$MAX_RRC[i] > 0.01) {mere_Expected_EDNA$Presence[i] <- 5}
  else if (mere_Expected_EDNA$locus_n[i] ==1 & mere_Expected_EDNA$DAFOR[i] == 1 & mere_Expected_EDNA$MAX_RRC[i] > 0.1) {mere_Expected_EDNA$Presence[i] <- 3}
  else if (mere_Expected_EDNA$locus_n[i] ==1 & mere_Expected_EDNA$DAFOR[i] >= 1 & mere_Expected_EDNA$MAX_RRC[i] > 0.05) {mere_Expected_EDNA$Presence[i] <- 3}
  else if (mere_Expected_EDNA$locus_n[i] ==2 & mere_Expected_EDNA$DAFOR[i] == 1 & mere_Expected_EDNA$MAX_RRC[i] > 0.01) {mere_Expected_EDNA$Presence[i] <- 3}
  else if (mere_Expected_EDNA$locus_n[i] ==2 & mere_Expected_EDNA$DAFOR[i] >= 1 & mere_Expected_EDNA$MAX_RRC[i] > 0.005) {mere_Expected_EDNA$Presence[i] <- 3}
  else if (mere_Expected_EDNA$locus_n[i] ==1 & mere_Expected_EDNA$DAFOR[i] >= 1 & mere_Expected_EDNA$MAX_RRC[i] > 0.005) {mere_Expected_EDNA$Presence[i] <- 1}
  else if (mere_Expected_EDNA$DAFOR[i] > 1 & mere_Expected_EDNA$MAX_RRC[i] < 0.005) {mere_Expected_EDNA$Presence[i] <- 1}
  else if (mere_Expected_EDNA$DAFOR[i] == 1 & mere_Expected_EDNA$MAX_RRC[i] < 0.005) {mere_Expected_EDNA$Presence[i] <- -1}
  else {mere_Expected_EDNA$Presence[i] <- 0} 
}


#Abundance score, corrected score (i.e. confidence score × site occupancy score)

for (i in 1:nrow(mere_Expected_EDNA)) {
  if (mere_Expected_EDNA$Presence[i]*mere_Expected_EDNA$DAFOR[i]>=25) {mere_Expected_EDNA$DAFOR_c[i] <- 5}
  else if (mere_Expected_EDNA$Presence[i]*mere_Expected_EDNA$DAFOR[i]>=15) {mere_Expected_EDNA$DAFOR_c[i] <- 4}
  else if (mere_Expected_EDNA$Presence[i]*mere_Expected_EDNA$DAFOR[i]>=9) {mere_Expected_EDNA$DAFOR_c[i] <- 3}
  else if (mere_Expected_EDNA$Presence[i]*mere_Expected_EDNA$DAFOR[i]>=4) {mere_Expected_EDNA$DAFOR_c[i] <- 2}
  else if (mere_Expected_EDNA$Presence[i]*mere_Expected_EDNA$DAFOR[i]>=1) {mere_Expected_EDNA$DAFOR_c[i] <- 1}
  else if (mere_Expected_EDNA$Presence[i]*mere_Expected_EDNA$DAFOR[i]==0) {mere_Expected_EDNA$DAFOR_c[i] <- 0}
  else {mere_Expected_EDNA$DAFOR_c[i] <- 1} 
}




##exclude unlikely present species




mere_excl_absent <- mere_Expected_EDNA[which(mere_Expected_EDNA$Presence != -1),]

#reorder based on the lake
mere_excl_absent$Lake <- factor(mere_excl_absent$Lake, levels = (c("Watch_Lane_Flash","Oss_Mere","Fenemere",
                                                                   "Maer_Pool","Betley_Mere","Chapel_Mere")))


levels(droplevels(mere_excl_absent$Species))



mere_excl_absent$Species <-mapvalues(mere_excl_absent$Species, c("Anguilla_anguilla","Abramis_brama","Carassius_carassius","Cyprinus_carpio", 
                                                                 "Leucaspius_delineatus","Phoxinus_phoxinus","Rutilus_rutilus","Scardinius_erythrophthalmus",
                                                                 "Tinca_tinca","Esox_lucius","Cottus_gobio","Gasterosteus_aculeatus",     
                                                                 "Perca_fluviatilis","Salmo_trutta","Gobio_gobio","Pseudorasbora_parva",        
                                                                                         "Rhodeus_amarus","Squalius_cephalus"), 
                                                             c("EEL","BRE","CRU","CAR",
                                                               "SUN","MIN","ROA","RUD",
                                                               "TEN","PIK","BUL","3SS",
                                                               "PER","BTR","GUD","TMG",
                                                               "BIT","CHU"))



ggplot(mere_excl_absent, aes(x = Ex_Presence, y = Presence))+
  geom_point(stat="identity",size=3)+
  facet_wrap(~Lake,ncol =2)+
  scale_y_continuous(limits = c(-1, 5), breaks=seq(-1, 5, 1))+
  scale_x_continuous(limits = c(-1, 5), breaks=seq(-1, 5, 1))+
  labs(x="Ex_Presence", y="Presence")+theme_bw()+
  theme(text=element_text(size=15))+
  geom_text_repel(aes(label=Species),size=3)

ggplot(mere_excl_absent, aes(x = Ex_Presence, y = Presence))+
  geom_point(stat="identity",size=3)+
  geom_smooth(formula= y ~ x,method="lm", se=FALSE,lwd=0.5,na.rm=TRUE, colour='black')+
  scale_y_continuous(limits = c(-1, 5), breaks=seq(-1, 5, 1))+
  scale_x_continuous(limits = c(-1, 5), breaks=seq(-1, 5, 1))+
  labs(x="Ex_Presence", y="Presence")+theme_bw()+
  theme(text=element_text(size=15))+
  geom_text_repel(aes(label=Species),size=3)

ggplot(mere_excl_absent, aes(x = Ex_Presence, y = Presence))+
  geom_point(stat="identity",size=3)+
  scale_y_continuous(limits = c(-1, 5), breaks=seq(-1, 5, 1))+
  scale_x_continuous(limits = c(-1, 5), breaks=seq(-1, 5, 1))+
  geom_smooth(formula= y ~ x,method="lm", se=FALSE,lwd=0.5,na.rm=TRUE, colour='black')+
  labs(x="Ex_Presence", y="Presence")+theme_bw()+
  theme(text=element_text(size=15))






mere_excl_absent[is.na(mere_excl_absent)] <- 1#Replace the NA with 1, otherwise the NA data will be removed


mere_excl_absent$Lake <-mapvalues(mere_excl_absent$Lake, c("Betley_Mere","Chapel_Mere","Fenemere","Maer_Pool","Oss_Mere","Watch_Lane_Flash"), 
                                                        c("BET","CAM","FEN","MAP","OSS","WLF"))



##confidence####

Wlake_confidence <- Wlake_excl_absent[c("Lake","Species","Project","Ex_Presence","Presence")]

mere_confidence <- mere_excl_absent[c("Lake","Species","Project","Ex_Presence","Presence")]




confidence <- rbind (Wlake_confidence, mere_confidence)


rhofun(confidence$Ex_Presence,confidence$Presence)$statistic

rhofun(confidence$Ex_Presence,confidence$Presence)$p.value

rhofun(confidence$Ex_Presence,confidence$Presence)$estimate

rhofun(confidence$Ex_Presence,confidence$Presence)$parameter


ggplot(confidence, aes(x = Ex_Presence, y = Presence))+
  geom_point(stat="identity",size=3)+
  scale_y_continuous(limits = c(-1, 5), breaks=seq(-1, 5, 1))+
  scale_x_continuous(limits = c(-1, 5), breaks=seq(-1, 5, 1))+
  geom_smooth(formula= y ~ x,method="lm", se=FALSE,lwd=0.5,na.rm=TRUE, colour='black')+
  labs(x="Ex_Presence", y="Presence")+theme_bw()+
  theme(text=element_text(size=15))



##correlation####

Wlake_cor <- Wlake_excl_absent[c("Lake","Species","Total_n","Project","Cytb_SO","sRNA_SO","Cytb_RC","sRNA_RC")]

mere_cor <- mere_excl_absent[c("Lake","Species","Total_n","Project","Cytb_SO","sRNA_SO","Cytb_RC","sRNA_RC")]


correlation <- rbind (Wlake_cor, mere_cor)

correlation$Cytb_reads <- correlation$Cytb_RC/correlation$Total_n

correlation$sRNA_reads <- correlation$sRNA_RC/correlation$Total_n

for (i in 1:nrow(correlation)) {
  if (correlation$Cytb_reads[i] ==0) {correlation$Cytb_reads[i] <- 1}
}

for (i in 1:nrow(correlation)) {
  if (correlation$sRNA_reads[i] ==0) {correlation$sRNA_reads[i] <- 1}
}





correlation$Lake <- factor(correlation$Lake, levels = (c("CWE","PAD","OGW","PEN","TRA","KEN","LLB","LLG","MAP","CAM","OSS","FEN","WLF","BET")))






corrfig1 <- ggplot(correlation, aes(x = Cytb_SO, y = Cytb_reads, colour=Lake))+
              geom_point(stat="identity",size=3)+
              scale_x_continuous(limits = c(0, 1), breaks=seq(0, 1.0, 0.2))+
              scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                            labels = trans_format("log10", math_format(10^.x)))+
              scale_colour_manual(name="Sampling lake",values=colpalette) +
              geom_smooth(formula= y ~ x,method="lm", se=FALSE,lwd=1.0,na.rm=TRUE)+
              guides(col = guide_legend(ncol = 7, byrow = TRUE))+
              labs(x="Site occupancy (Cytb)", y="Averaged read counts (Cytb)")+theme_bw()+
              stat_cor(method = "pearson",label.sep = "; ",label.y.npc="centre",label.x=0.6,size=3, hjust = 0, show.legend = FALSE)+
              theme(text=element_text(size=18), legend.position = "bottom",legend.title = element_text(size=15),legend.text = element_text(size=15),
                    panel.grid.minor=element_blank(),panel.grid.major =element_blank())





legend_corrfig1 <- get_legend(corrfig1)

corrfig1 <- corrfig1 + theme(legend.position="none")+
             geom_text(aes(x=0, y=1e4, label= '(a)'),hjust = 0, size=6, colour = "black",parse = TRUE)


corrfig2 <- ggplot(correlation, aes(x = sRNA_SO, y = sRNA_reads, colour=Lake))+
            geom_point(stat="identity",size=3)+
            scale_x_continuous(limits = c(0, 1), breaks=seq(0, 1.0, 0.2))+
            scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                          labels = trans_format("log10", math_format(10^.x)))+
            scale_colour_manual(name="Sampling lake",values=colpalette) +
            geom_smooth(formula= y ~ x,method="lm", se=FALSE,lwd=1.0,na.rm=TRUE)+
            labs(x="Site occupancy (12S)", y="Averaged read counts (12S)")+theme_bw()+
            stat_cor(method = "pearson",label.sep = "; ",label.y.npc="centre",label.x=0.6,size=3, hjust = 0)+
            theme(text=element_text(size=18),legend.position = "none",
                  panel.grid.minor=element_blank(),panel.grid.major =element_blank())



corrfig2 <- corrfig2 +
  geom_text(aes(x=0, y=2e4, label= '(b)'),hjust = 0, size=6, colour = "black",parse = TRUE)




corrfig3 <- ggplot(correlation, aes(x = sRNA_SO, y = Cytb_SO, colour=Lake))+
            geom_point(stat="identity",size=3)+
            scale_x_continuous(limits = c(0, 1), breaks=seq(0, 1.0, 0.2))+
            scale_y_continuous(limits = c(0, 1), breaks=seq(0, 1.0, 0.2))+
            scale_colour_manual(name="Sampling lake",values=colpalette) +
            geom_smooth(formula= y ~ x,method="lm", se=FALSE,lwd=1,na.rm=TRUE)+
            labs(x="Site occupancy (12S)", y="Site occupancy (Cytb)")+theme_bw()+
            stat_cor(method = "pearson",label.sep = "; ",label.y.npc="centre",label.x=0.7,size=3, hjust = 0)+
            theme(text=element_text(size=18),legend.position = "none",
                  panel.grid.minor=element_blank(),panel.grid.major =element_blank())



corrfig3 <- corrfig3 +
  geom_text(aes(x=0, y=1.0, label= '(c)'),hjust = 0, size=6, colour = "black",parse = TRUE)




corrfig4 <- ggplot(correlation, aes(x = sRNA_reads, y = Cytb_reads,colour=Lake))+
            geom_point(stat="identity",size=3)+
            scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                          labels = trans_format("log10", math_format(10^.x)))+
            scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                          labels = trans_format("log10", math_format(10^.x)))+
            scale_colour_manual(name="Sampling lake",values=colpalette) +
            geom_smooth(formula= y ~ x,method="lm", se=FALSE,lwd=1.0,na.rm=TRUE)+
            labs(x="Averaged read counts (12S)", y="Averaged read counts (Cytb)")+theme_bw()+
            stat_cor(method = "pearson",label.sep = "; ",label.y.npc="centre",label.x=2.8,size=3, hjust = 0)+
            theme(text=element_text(size=18),legend.position = "none",
                  panel.grid.minor=element_blank(),panel.grid.major =element_blank())





corrfig4 <- corrfig4 +
  geom_text(aes(x=0.1, y=1e4, label= '(d)'),hjust = 0, size=6, colour = "black",parse = TRUE)



corr_fig1 <- grid.arrange(corrfig1,corrfig2,corrfig3,corrfig4, ncol = 2, nrow = 2,
                                  widths = c(2.7,2.7), heights = c(2.5,2.5))


corr_fig2 <- grid.arrange(corr_fig1, legend_corrfig1, nrow = 2, heights = c(2.5,0.3))


#FigS1.2_corr####

ggsave("FigS1.2_corr.jpeg",corr_fig2,path = "Lake project/Figures_new/",  width = 10, height = 10, units = "in", dpi=500)





#eDNA abudance estimating####


Wlake_excl_absent$Lake <- factor(Wlake_excl_absent$Lake, levels = (c("CWE","PAD","OGW","PEN","TRA","KEN","LLB","LLG")))




DAFOR_c_welsh <-ggplot(Wlake_excl_absent, aes(x = Ex_DAFOR, y = DAFOR_c))+
                geom_point(stat="identity",size=2.5)+
                facet_wrap(~Lake,ncol =4)+
                scale_y_continuous(limits = c(0, 5), breaks=seq(0, 5, 1))+
                scale_x_continuous(limits = c(0, 5), breaks=seq(0, 5, 1))+
                geom_smooth(formula= y ~ x,method="lm", se=FALSE,lwd=0.5,na.rm=TRUE, colour='black')+
                stat_cor(method = "spearman",label.sep = "\n",label.y=1,label.x=3,size=3.5,hjust = 0)+
                labs(x="Expected DAFOR scale", y="DAFOR scale based on eDNA data")+theme_bw()+
                theme(text=element_text(size=12),strip.background = element_blank(),strip.text=element_text(face="bold"),
                      panel.grid.minor=element_blank(),panel.grid.major =element_blank())+
                geom_text_repel(aes(label=Species),size=3,fontface ="italic")




DAFOR_c_len <- length(levels(Wlake_excl_absent$Lake))

DAFOR_c_vars <- data.frame(expand.grid(levels(Wlake_excl_absent$Lake)))

colnames(DAFOR_c_vars) <- c("Lake")


DAFOR_c_dat3 <- data.frame(x = rep(0, DAFOR_c_len), y = rep(4.9, DAFOR_c_len), 
                           DAFOR_c_vars, labs=paste0("(",letters[1:8],")"))




DAFOR_1<-DAFOR_c_welsh +
          geom_text(data=DAFOR_c_dat3, aes(x, y, label=labs),hjust = 0, size=3,fontface ="bold",parse = TRUE)





colpalette2 <-c("#8cd6b5",
               "#4f2d46",
               "#638d48",
               "#c48943",
               "#669bb7",
               "#b74ace",
               "#cb9f9c",
               "#cfd058")

# 
# 
# DAFOR_welsh <-ggplot(Wlake_excl_absent, aes(x = Ex_DAFOR, y = DAFOR_c, colour=Lake))+
#               geom_point(stat="identity",size=3,position = position_jitterdodge(dodge.width = 0.02))+
#               scale_y_continuous(limits = c(0, 5), breaks=seq(0, 5, 1))+
#               scale_x_continuous(limits = c(0, 5), breaks=seq(0, 5, 1))+
#               scale_colour_manual(name="Sampling site",values=colpalette2) +
#               geom_smooth(formula= y ~ x,method="lm", se=FALSE,lwd=1.0,na.rm=TRUE)+
#               labs(x="Expected DAFOR scale", y="DAFOR scale based on eDNA data")+theme_bw()+
#               stat_cor(method = "spearman",label.sep = "; ",label.y.npc="centre",label.x=2.8,size=5, hjust = 0)+
#               theme(text=element_text(size=18),legend.position = "none",
#                     panel.grid.minor=element_blank(),panel.grid.major =element_blank())
# 
# ggplot(Wlake_excl_absent, aes(x = Ex_DAFOR, y = DAFOR_c))+
#   geom_point(stat="identity",size=3)+
#   scale_y_continuous(limits = c(0, 5), breaks=seq(0, 5, 1))+
#   scale_x_continuous(limits = c(0, 5), breaks=seq(0, 5, 1))+
#   facet_wrap(~Lake,ncol =4)+
#   scale_colour_manual(name="Sampling site",values=colpalette2) +
#   geom_smooth(formula= y ~ x,method="lm", se=FALSE,lwd=1.0,na.rm=TRUE)+
#   labs(x="Expected DAFOR scale", y="DAFOR scale based on eDNA data")+theme_bw()+
#   stat_cor(method = "spearman",label.sep = "; ",label.y.npc="centre",label.x=2.8,size=5, hjust = 0)+
#   theme(text=element_text(size=18),legend.position = "none",
#         panel.grid.minor=element_blank(),panel.grid.major =element_blank())+
#   geom_text_repel(aes(label=Species),size=3)
# 
# 
# 
# ddply(Wlake_excl_absent, .(Lake), summarise, s=rhofun(Ex_DAFOR,DAFOR_c)$statistic,
#       cor.est=rhofun(Ex_DAFOR,DAFOR_c)$estimate,
#       pval=rhofun(Ex_DAFOR,DAFOR_c)$p.value,
#       alt=rhofun(Ex_DAFOR,DAFOR_c)$alternative)
# 
# 
# ddply(Wlake_excl_absent, .(Lake), summarise, s=corfun(Ex_DAFOR,DAFOR_c)$statistic,
#       cor.est=corfun(Ex_DAFOR,DAFOR_c)$estimate,
#       pval=corfun(Ex_DAFOR,DAFOR_c)$p.value,
#       alt=corfun(Ex_DAFOR,DAFOR_c)$alternative,
#       df=corfun(Ex_DAFOR,DAFOR_c)$parameter)
# 
# 
# 
# 
# DAFOR_welsh <- DAFOR_welsh + 
#   geom_text(aes(x=0, y=5, label= 'A'),hjust = 0, size=7, colour = "black",parse = TRUE)
# 
# 





mere_excl_absent$Lake <- factor(mere_excl_absent$Lake, levels = (c("MAP","CAM","OSS","FEN","WLF","BET")))




colpalette3 <-c("#74d656",
                "#a688c9",
                "#c34980",
                "#5541a7",
                "#4b4d32",
                "#c34837")





DAFOR_Ind<- ggplot(mere_excl_absent, aes(x = Ind_abu, y = DAFOR_c))+
            geom_point(stat="identity",size=2.5)+
            facet_wrap(~Lake,ncol =3)+
            scale_y_continuous(limits = c(0, 5), breaks=seq(0, 5, 1))+
            scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                          labels = trans_format("log10", math_format(10^.x)))+
            geom_smooth(formula= y ~ x,method="lm", se=FALSE,lwd=0.5,na.rm=TRUE, colour='black')+
            stat_cor(method = "spearman",label.sep = "\n",label.y=1,label.x=1.4,size=3.5,hjust = 0)+
            labs(x= expression(paste("Individual density (ind. ha"^"-1",")")), y="DAFOR scale based on eDNA data")+theme_bw()+
            theme(text=element_text(size=12),strip.background = element_blank(),strip.text=element_text(face="bold"),
                  panel.grid.minor=element_blank(),panel.grid.major =element_blank())+
            geom_text_repel(aes(label=Species),size=3,fontface ="italic")






DAFOR_Ind_len <- length(levels(mere_excl_absent$Lake))

DAFOR_Ind_vars <- data.frame(expand.grid(levels(mere_excl_absent$Lake)))

colnames(DAFOR_Ind_vars) <- c("Lake")


DAFOR_Ind_dat3 <- data.frame(x = rep(1, DAFOR_Ind_len), y = rep(4.9, DAFOR_Ind_len), 
                           DAFOR_Ind_vars,labs=paste0("(",letters[9:14],")"))



DAFOR_2<-DAFOR_Ind + geom_text(data=DAFOR_Ind_dat3, aes(x, y, label=labs),hjust = 0, size=3, fontface ="bold", parse = TRUE)



DAFOR_Bio<- ggplot(mere_excl_absent, aes(x = Bio_den, y = DAFOR_c))+
            geom_point(stat="identity",size=2.5)+
            facet_wrap(~Lake,ncol =3)+
            scale_y_continuous(limits = c(0, 5), breaks=seq(0, 5, 1))+
            scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                          labels = trans_format("log10", math_format(10^.x)))+
            geom_smooth(formula= y ~ x,method="lm", se=FALSE,lwd=0.5,na.rm=TRUE, colour='black')+
            stat_cor(method = "spearman",label.sep = "\n",label.y=1,label.x=1.2,size=3.5,hjust = 0)+
            labs(x= expression(paste("Biomass density (kg ha"^"-1",")")), y="DAFOR scale based on eDNA data")+theme_bw()+
            theme(text=element_text(size=12),strip.background = element_blank(),strip.text=element_text(face="bold"),
                  panel.grid.minor=element_blank(),panel.grid.major =element_blank())+
            geom_text_repel(aes(label=Species),size=3,fontface ="italic")



DAFOR_Bio_len <- length(levels(mere_excl_absent$Lake))

DAFOR_Bio_vars <- data.frame(expand.grid(levels(mere_excl_absent$Lake)))

colnames(DAFOR_Bio_vars) <- c("Lake")


DAFOR_Bio_dat3 <- data.frame(x = rep(0.9, DAFOR_Bio_len), y = rep(4.9, DAFOR_Bio_len), 
                             DAFOR_Bio_vars, labs=paste0("(",letters[1:6],")"))



DAFOR_3<-DAFOR_Bio + geom_text(data=DAFOR_Bio_dat3, aes(x, y, label=labs),hjust = 0, size=3, fontface ="bold", parse = TRUE)



#FigS1.4_Bio####
ggsave("FigS1.4_Bio.jpeg",DAFOR_3, path = "Lake project/Figures_new/",  width = 6.5, height = 4, units = "in", dpi=500)



# DAFOR_FNL <- grid.arrange(DAFOR_1,DAFOR_2,DAFOR_3, ncol = 1, nrow = 3,
#                           widths = c(2.7), heights = c(2.5,2.5,2.5))
# 
# ggsave("Fig3_DAFOR.jpeg",DAFOR_FNL, path = "Lake project/Figures_new/",  width = 10, height = 10, units = "in", dpi=500)



DAFOR_FNL <- grid.arrange(DAFOR_1,DAFOR_2, ncol = 1, nrow = 2,
                          widths = c(2.7), heights = c(2.0,2.0))



#Fig3_DAFOR####
ggsave("Fig3_DAFOR.jpeg",DAFOR_FNL, path = "Lake project/Figures_new/",  width = 8, height = 9, units = "in", dpi=500)

ggsave("Fig3_DAFOR.pdf",DAFOR_FNL, path = "Lake project/Figures_new/",  width = 9, height = 10, units = "in", dpi=500)

# DAFOR_Ind <-  ggplot(mere_excl_absent, aes(x = Ind_abu, y = DAFOR_c, colour=Lake))+
#               geom_point(stat="identity",size=3,position = position_jitterdodge(dodge.width = 0.02))+
#               scale_y_continuous(limits = c(0, 5), breaks=seq(0, 5, 1))+
#               scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
#                             labels = trans_format("log10", math_format(10^.x)))+
#               geom_smooth(formula= y ~ x,method="lm", se=FALSE,lwd=1,na.rm=TRUE)+
#               scale_colour_manual(name="Sampling site",values=colpalette3) +
#               labs(x="Individual density (ind./ha)", y="")+theme_bw()+
#               stat_cor(method = "spearman",label.sep = "; ",label.y.npc="centre",label.x=2.1,size=5, hjust = 0, show.legend = FALSE)+
#               theme(text=element_text(size=18),legend.position = "none",
#                     panel.grid.minor=element_blank(),panel.grid.major =element_blank(),axis.title.y=element_blank(), axis.ticks.y=element_blank(),axis.text.y=element_blank())
# 
# 
# 
# 
# 
# DAFOR_Ind <- DAFOR_Ind  + 
#   geom_text(aes(x=0.7, y=5, label= 'B'),hjust = 0, size=7, colour = "black",parse = TRUE)
# 
# 
# 
# DAFOR_fig1 <- grid.arrange(DAFOR_welsh, DAFOR_Ind, ncol = 2, nrow = 1,
#                            widths = c(2.7,2.7), heights = c(2.5))
# 
# 
# 
# DAFOR_fig2 <- grid.arrange(DAFOR_fig1, legend_corrfig1, nrow = 2,
#                            heights = c(2.5,0.3))
# 
# 
# #ggsave("Fig5_DAFOR.jpeg",DAFOR_fig2, path = "Lake project/Figures_new/",  width = 11, height = 6, units = "in", dpi=500)
# 
# 
# DAFOR_Bio <- ggplot(mere_excl_absent, aes(x = Bio_den, y = DAFOR_c, colour=Lake))+
#               geom_point(stat="identity",size=3,position = position_jitterdodge(dodge.width = 0.02))+
#               scale_y_continuous(limits = c(0, 5), breaks=seq(0, 5, 1))+
#               scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
#                             labels = trans_format("log10", math_format(10^.x)))+
#               geom_smooth(formula= y ~ x,method="lm", se=FALSE,lwd=1,na.rm=TRUE)+
#               scale_colour_manual(name="Sampling site",values=colpalette3) +
#               guides(col = guide_legend(ncol = 6, byrow = TRUE))+
#               labs(x="Biomass density (kg/ha)", y="DAFOR scale based on eDNA data")+theme_bw()+
#               stat_cor(method = "spearman",label.sep = "; ",label.y.npc="centre",label.x=1.1,size=5, hjust = 0, show.legend = FALSE)+
#               theme(text=element_text(size=13),legend.position = "bottom",legend.title = element_text(size=11),legend.text = element_text(size=11),
#                     panel.grid.minor=element_blank(),panel.grid.major =element_blank())
# 
#   
#   
# #ggsave("FigS3_DAFOR_bio.jpeg",DAFOR_Bio, path = "Lake project/Figures_new/",  width = 6, height = 6, units = "in", dpi=500)
# 
# 
# 
# 








###fish composition####


Wlake_Fish_SN <- dcast(LP_rs_wlake_SO2,  Lake+Total_n+Project+Locus~Species, fun.aggregate = length)

Wlake_Fish_SN.melt = melt(Wlake_Fish_SN, id= c("Lake","Total_n","Project","Locus"))


Wlake_Fish_SN.melt <- rename (Wlake_Fish_SN.melt,c("variable"="Species","value"="SN"))

Wlake_Fish_SN.melt <- Wlake_Fish_SN.melt[which(Wlake_Fish_SN.melt$SN >0 ),]


Wlake_fish_com  <- merge(LP_rs_wlake_SO_sum,Wlake_Fish_SN.melt,by=c("Lake","Total_n","Project","Species","Locus"))



Wlake_Expected_EDNA[which(Wlake_Expected_EDNA$Presence == -1),]

#excduding the species records based on the above results


Wlake_fish_com <- Wlake_fish_com[which(Wlake_fish_com$Lake !="Kenfig_Pool" |
                                         Wlake_fish_com$Species != "Leuciscus_idus"),]


Wlake_fish_com <- Wlake_fish_com[which(Wlake_fish_com$Lake !="Llyn_Cwellyn" |
                                         Wlake_fish_com$Species != "Rutilus_rutilus"),]


Wlake_fish_com <- Wlake_fish_com[which(Wlake_fish_com$Lake !="Llyn_Ogwen" |
                                         Wlake_fish_com$Species != "Abramis_brama"),]

Wlake_fish_com <- Wlake_fish_com[which(Wlake_fish_com$Lake !="Llyn_Ogwen" |
                                         Wlake_fish_com$Species != "Cottus_gobio"),]


Wlake_fish_com <- Wlake_fish_com[which(Wlake_fish_com$Lake !="Llyn_Penrhyn" |
                                         Wlake_fish_com$Species != "Esox_lucius"),]



Wlake_fish_com <- Wlake_fish_com[which(Wlake_fish_com$Lake !="Llyn_Traffwll" |
                                         Wlake_fish_com$Species != "Cottus_gobio"),]


Wlake_fish_com <- Wlake_fish_com[which(Wlake_fish_com$Lake !="Llyn_Padarn" |
                                         Wlake_fish_com$Species != "Salmo"),]




mere_Fish_SN <- dcast(LP_rs_mere_SO2,  Lake+Total_n+Project+Locus~Species, fun.aggregate = length)

mere_Fish_SN.melt = melt(mere_Fish_SN, id= c("Lake","Total_n","Project","Locus"))


mere_Fish_SN.melt <- rename (mere_Fish_SN.melt,c("variable"="Species","value"="SN"))

mere_Fish_SN.melt <- mere_Fish_SN.melt[which(mere_Fish_SN.melt$SN >0 ),]


mere_fish_com  <- merge(LP_rs_mere_SO_sum,mere_Fish_SN.melt,by=c("Lake","Total_n","Project","Species","Locus"))


#excduding the species records based on the above results

mere_Expected_EDNA[which(mere_Expected_EDNA$Presence == -1),]

mere_fish_com <- mere_fish_com[which(mere_fish_com$Lake !="Betley_Mere" |
                                       mere_fish_com$Species != "Gobio_gobio"),]





Fish_com <- rbind(Wlake_fish_com, mere_fish_com)

Fish_com <- ddply(Fish_com, .(Lake,Locus), mutate, Percent_reads = Reads/sum(Reads))

Fish_com <- ddply(Fish_com, .(Lake,Locus), mutate, Percent_SN = SN/sum(SN))




Fish_com$Lake <-mapvalues(Fish_com$Lake, c("Betley_Mere","Chapel_Mere","Fenemere","Kenfig_Pool","Llan_Bwch-llyn","Llangorse_Lake",  
                                           "Llyn_Cwellyn","Llyn_Ogwen","Llyn_Padarn","Llyn_Penrhyn","Llyn_Traffwll","Maer_Pool",       
                                           "Oss_Mere","Watch_Lane_Flash"), 
                                          c("BET","CAM","FEN","KEN","LLB","LLG","CWE","OGW","PAD","PEN","TRA","MAP","OSS","WLF"))



Fish_com$Lake <- factor(Fish_com$Lake, levels = (c("CWE","PAD","OGW","PEN","TRA","KEN","LLB","LLG","MAP","CAM","OSS","FEN","WLF","BET")))

Fish_com$Locus <- factor(Fish_com$Locus, levels = (c("12S","Cytb")))



Fish_com$Species <- factor(Fish_com$Species, levels = (c("Abramis_brama","Alburnus_alburnus","Anguilla_anguilla","Carassius_auratus","Carassius_carassius","Cottus_gobio",
                                                         "Cyprinus_carpio","Esox_lucius", "Gasterosteus_aculeatus","Gobio_gobio","Leucaspius_delineatus","Oncorhynchus_mykiss","Perca_fluviatilis",
                                                         "Phoxinus_phoxinus","Pseudorasbora_parva","Pungitius_pungitius","Rhodeus_amarus","Rutilus_rutilus",
                                                         "Salmo_salar","Salmo_trutta","Salvelinus_alpinus", "Scardinius_erythrophthalmus","Squalius_cephalus","Tinca_tinca")))

levels(droplevels(Fish_com$Species))


Fish_com$Species <-mapvalues(Fish_com$Species, c("Abramis_brama","Alburnus_alburnus","Anguilla_anguilla","Carassius_auratus","Carassius_carassius","Cottus_gobio",
                                                 "Cyprinus_carpio","Esox_lucius", "Gasterosteus_aculeatus","Gobio_gobio","Leucaspius_delineatus","Oncorhynchus_mykiss","Perca_fluviatilis",
                                                 "Phoxinus_phoxinus","Pseudorasbora_parva","Pungitius_pungitius","Rhodeus_amarus","Rutilus_rutilus",
                                                 "Salmo_salar","Salmo_trutta","Salvelinus_alpinus", "Scardinius_erythrophthalmus","Squalius_cephalus","Tinca_tinca"), 
                                               c("BRE","BLE","EEL","GOF","CRU","BUL",
                                                 "CAR","PIK","3SS","GUD","SUN","RTR","PER",
                                                 "MIN","TMG","9SS","BIT","ROA",
                                                 "SAL","BTR","CHA","RUD","CHU","TEN"))
                  

###http://tools.medialab.sciences-po.fr/iwanthue/

palette2 <- c("#d28939",
              "#7440ce",
              "#66d64f",
              "#ce4dc5",
              "#ded23d",
              "#05a496",
              "#a1cb52",
              "#7578aa",
              "#6cd6a0",
              "#d54b75",
              "#4b7f46",
              "#d84d35",
              "#8acccd",
              "#00d3ce",
              "#839fcb",
              "#867d2f",
              "#7d4374",
              "#d2cf92",
              "#074cf8",
              "#d9a395",
              "#533432",
              "#cf91c7",
              "#4a737a",
              "#009100")

ggplot(Fish_com,aes(x=Locus,y=Percent_reads,fill=Species))+
  geom_bar(stat="identity",position="stack",width = 0.8,colour="black")+facet_grid(Lake ~ ., switch = "both")+coord_flip()+
  scale_y_continuous(labels=percent,expand = c(0,0),limits = c(0, 1.01), breaks=seq(0, 1, 0.25))+
  scale_fill_manual(name="Species",values=palette2)+
  labs(x="Sampling lake",y="Proportion of read counts")+
  guides(fill = guide_legend(ncol = 8,byrow = TRUE, reverse= TRUE))+
  theme_bw()+theme(text=element_text(size=15),legend.position="bottom",legend.text = element_text(size=12),
                   strip.background = element_blank(),panel.border = element_blank(),strip.text.y=element_text(angle=180),
                   strip.placement= "outside",strip.text=element_text(face="bold"),axis.ticks.length = unit(0.2, "cm"),
                   legend.title = element_text(size=12),panel.grid.minor=element_blank(),panel.grid.major =element_blank())

###FigS1.3_Fish_com_reads####


ggsave("FigS1.3_Fish_com_reads.jpeg",path = "Lake project/Figures_new/",  width = 15, height = 8, units = "in", dpi=500)




####SN####



ggplot(Fish_com,aes(x=Locus,y=Percent_SN,fill=Species))+
  geom_bar(stat="identity",position="stack",width = 0.8,colour="black")+facet_grid(Lake ~ ., switch = "both")+coord_flip()+
  scale_y_continuous(labels=percent,expand = c(0,0),limits = c(0, 1.01), breaks=seq(0, 1, 0.25))+
  scale_fill_manual(name="Species", values=palette2)+
  labs(x="Sampling lake",y="Proportion of site occupancy")+
  guides(fill = guide_legend(ncol = 8,byrow = TRUE, reverse= TRUE))+
  theme_bw()+theme(text=element_text(size=15),legend.position="bottom",legend.text = element_text(size=12),
                  strip.background = element_blank(),panel.border = element_blank(),strip.text.y=element_text(angle=180),
                  strip.placement= "outside",strip.text=element_text(face="bold"),axis.ticks.length = unit(0.2, "cm"),
                  legend.title = element_text(size=12),panel.grid.minor=element_blank(),panel.grid.major =element_blank())

###Fig2_Fish_com####

ggsave("Fig2_Fish_com_SN.jpeg",path = "Lake project/Figures_new/",  width = 15, height = 8, units = "in", dpi=500)

ggsave("Fig2_Fish_com_SN.pdf",path = "Lake project/Figures_new/",  width = 15, height = 8, units = "in", dpi=500)


###Cluster with SO####

Fish_com$SO <- Fish_com$SN/Fish_com$Total_n
  
Fish_cluster <-  Fish_com[c(1,4,5,10)]




Fish_clu_cytb <- Fish_cluster[which(Fish_cluster$Locus == 'Cytb'),]
  
Fish_clu_cytb <- Fish_clu_cytb[c(-3)]


Fish_clu_cytb_wide <- reshape(Fish_clu_cytb, idvar= c("Lake"), timevar= "Species", direction = "wide")

Fish_clu_cytb_wide <- data.frame(Fish_clu_cytb_wide[,-1], row.names=Fish_clu_cytb_wide[,1])


Fish_clu_cytb_wide[is.na(Fish_clu_cytb_wide)] <- 0



library("cluster")
library("factoextra")


# Visualize



res.hc_cytb <- Fish_clu_cytb_wide %>%
  scale() %>%                    # Scale the data
  dist(method = "canberra") %>% # Compute dissimilarity matrix
  hclust(method = "ward.D2")     # Compute hierachical clustering


res.dist_Cytb <- dist(Fish_clu_cytb_wide, method = "canberra")
res.coph_Cytb <- cophenetic(res.hc_cytb)
# Correlation between cophenetic distance and the original distance
cor(res.dist_Cytb, res.coph_Cytb)

# Visualize using factoextra

Clu_Cytb <- fviz_dend(res.hc_cytb, k = 4, # Cut in three groups
                      cex = 1,lwd = 1, xlab = "Sampling lake", main="",ggtheme = theme_bw(base_size = 20),# label size
                      k_colors = c("#669bb7", "#638d48", "#5541a7","#ab641c"),
                      color_labels_by_k = TRUE, # color labels by groups
                      rect = TRUE # Add rectangle around groups
                      )+theme(panel.grid.minor=element_blank(),panel.grid.major =element_blank())


Clu_Cytb1 <-  ggdraw(Clu_Cytb)+draw_plot_label(c("Type 2","Type 1","Type 3","Type 3"), c(0.12,0.28,0.55,0.8), c(0.6,0.6,0.6,0.6), fontface = "plain",size = 18,colour="red")





Fish_clu_sRNA <- Fish_cluster[which(Fish_cluster$Locus == '12S'),]

Fish_clu_sRNA <- Fish_clu_sRNA[c(-3)]


Fish_clu_sRNA_wide <- reshape(Fish_clu_sRNA, idvar= c("Lake"), timevar= "Species", direction = "wide")

Fish_clu_sRNA_wide <- data.frame(Fish_clu_sRNA_wide[,-1], row.names=Fish_clu_sRNA_wide[,1])


Fish_clu_sRNA_wide[is.na(Fish_clu_sRNA_wide)] <- 0



# Visualize



res.hc_sRNA <- Fish_clu_sRNA_wide %>%
  scale() %>%                    # Scale the data
  dist(method = "canberra") %>% # Compute dissimilarity matrix
  hclust(method = "ward.D2")     # Compute hierachical clustering

?dist

res.dist_sRNA <- dist(Fish_clu_sRNA_wide, method = "canberra")
res.coph_sRNA <- cophenetic(res.hc_sRNA)
# Correlation between cophenetic distance and the original distance
cor(res.dist_sRNA, res.coph_sRNA)

# Visualize using factoextra

Clu_sRNA <- fviz_dend(res.hc_sRNA, k = 4, # Cut in three groups
                      cex = 1,lwd = 1, xlab = "Sampling lake", main = "", ggtheme = theme_bw(base_size = 20),# label size
                      k_colors = c("#638d48", "#ab641c","#669bb7", "#5541a7"),
                      color_labels_by_k = TRUE, # color labels by groups
                      rect = TRUE # Add rectangle around groups
                      )+theme(panel.grid.minor=element_blank(),panel.grid.major =element_blank())


Clu_sRNA1 <-  ggdraw(Clu_sRNA)+draw_plot_label(c("Type 1","Type 3","Type 2","Type 3"), c(0.12,0.28,0.4,0.7), c(0.67,0.67,0.67,0.67), fontface = "plain",size = 18,colour="red")



dend_cytb <- Fish_clu_cytb_wide %>%
  scale() %>%                    # Scale the data
  dist(method = "canberra") %>% # Compute dissimilarity matrix
  hclust(method = "ward.D2")%>% as.dendrogram
  
dend_sRNA <- Fish_clu_sRNA_wide %>%
  scale() %>%                    # Scale the data
  dist(method = "canberra") %>% # Compute dissimilarity matrix
  hclust(method = "ward.D2")%>% as.dendrogram


library(dendextend)

dend_list <- dendlist(dend_cytb, dend_sRNA)

cor.dendlist(dend_list, method = "cophenetic")

dendlist(dend_cytb, dend_sRNA) %>%
  untangle(method = "step1side") %>% # Find the best alignment layout
  tanglegram()                       # Draw the two dendrograms



####Fig4####

NMDS_Clu_fig1 <- grid.arrange(Clu_Cytb1,Clu_sRNA1,NMDS_Cytb_FIG,NMDS_12S_FIG, ncol = 2, nrow = 2,
                              widths = c(2.7,2.7), heights = c(2.2,2.5))


NMDS_Clu_fig2 <- grid.arrange(NMDS_Clu_fig1, legend_NMDS, nrow = 2, heights = c(2.5,0.2))



NMDS_Clu_final <- ggdraw(NMDS_Clu_fig2)+draw_plot_label(c("(a1)","(a2)","(b1)","(b2)"), c(0, 0, 0.49,0.49), c(0.99, 0.58, 0.99, 0.58), fontface = "plain",size = 15)


ggsave("Fig4_NMDS_Clu.jpeg",NMDS_Clu_final,path = "Lake project/Figures_new/",  width = 15, height = 12, units = "in", dpi=500)

ggsave("Fig4_NMDS_Clu.pdf",NMDS_Clu_final,path = "Lake project/Figures_new/",  width = 15, height = 12, units = "in", dpi=500)



####Sample-based rarefaction analyses####

#install package rich Rossi Diversity 2011, 3, 112-120
#requires installation of vegan
#requires permute
#requires boot

#load rich library
library(permute)
library(vegan)
library(rich)
library(boot)


#Sequencing depth#####

names(LP_rs)

LP_rs_rare <- LP_rs[c(1,11,13,14)]

LP_rs_rare<-LP_rs_rare[which(LP_rs_rare$Reads>0),]

LP_rs_rare_Cytb <- LP_rs_rare[which(LP_rs_rare$Locus=="Cytb"),]

names(LP_rs_rare_Cytb)

LP_rs_rare_Cytb_sum <- aggregate(Reads~Locus+Species, data = LP_rs_rare_Cytb, FUN=mean)


LP_rs_rare_Cytb_sum$Reads <- round(LP_rs_rare_Cytb_sum$Reads,0)

LP_rs_rare_Cytb_sum_b <- reshape(LP_rs_rare_Cytb_sum, idvar= c("Locus"), timevar= "Species", direction = "wide")

LP_rs_rare_Cytb_sum_b$Locus<-LP_rs_rare_Cytb_sum_b$Reads.Salmo <- NULL




LP_rs_rare_sRNA <- LP_rs_rare[which(LP_rs_rare$Locus=="12S"),]

names(LP_rs_rare_sRNA)

LP_rs_rare_sRNA_sum <- aggregate(Reads~Locus+Species, data = LP_rs_rare_sRNA, FUN=mean)


LP_rs_rare_sRNA_sum$Reads <- round(LP_rs_rare_sRNA_sum$Reads,0)

LP_rs_rare_sRNA_sum<- LP_rs_rare_sRNA_sum[which(LP_rs_rare_sRNA_sum$Reads>0),]


LP_rs_rare_sRNA_sum_b <- reshape(LP_rs_rare_sRNA_sum, idvar= c("Locus"), timevar= "Species", direction = "wide")

LP_rs_rare_sRNA_sum_b$Locus <- LP_rs_rare_sRNA_sum_b$Reads.Leuciscus_idus<-NULL




#FigS1.1_sequencing_depth#####


# Make plot


png("Lake project/Figures_new/FigS1.1_sequencing_depth.jpeg", width = 12, height = 6, units = 'in', res = 500, pointsize = 12)

layout(matrix(c(1,2), 1, 2, byrow = TRUE),
       heights=c(2)) #arrange the plots 


rarecurve(LP_rs_rare_Cytb_sum_b, step = 10, sample=7436, 
          xlab = "Number of sequencing read counts", ylab = "Number of species",
          label = FALSE,col="blue",lwd=3,cex.lab=1.5,yaxt="n")
axis(2, at=seq(0, 18, by=2) , labels= seq(0, 18, by=2))   


mtext('(a)', side = 3, line = 1, adj = 0.0, cex = 1.5) 

rarecurve(LP_rs_rare_sRNA_sum_b, step = 10, sample=5743, 
          xlab = "Number of sequencing read counts", ylab = "Number of species",
          label = FALSE,col="blue",lwd=3,cex.lab=1.5,yaxt="n")
axis(2, at=seq(0, 24, by=2) , labels= seq(0, 24, by=2))  


mtext('(b)', side = 3, line = 1, adj = 0.0, cex = 1.5) 

dev.off()






###### transform back


NMDS_12S_SL <- NMDS_12S[which(NMDS_12S$Transect=="SL"),]

NMDS_Cytb_SL <- NMDS_Cytb[which(NMDS_Cytb$Transect=="SL"),]


RFA_12S <- NMDS_12S_SL[c(1,2,4,5)]

RFA_Cytb <- NMDS_Cytb_SL[c(1,2,4,5)]


# CWE ----------------------------------------------------------------------
RFA_12S_CWE <- RFA_12S[which(RFA_12S$Lake == 'CWE'),]


RFA_12S_CWE_back <- reshape(RFA_12S_CWE, idvar= c("SampleID","Lake"), timevar= "Species", direction = "wide")


RFA_12S_CWE_back[nrow(RFA_12S_CWE_back) + 1,] = c("CWE14","CWE",0,0,0,0,0)  
  

  

RFA_12S_CWE_noID <- RFA_12S_CWE_back [c(3:ncol(RFA_12S_CWE_back))]

RFA_12S_CWE_noID [is.na(RFA_12S_CWE_noID)] <- 0

RFA_12S_CWE_noID <- mutate_all(RFA_12S_CWE_noID, function(x) as.numeric(as.character(x)))



test_CWE_12S <- rich(matrix=RFA_12S_CWE_noID, nrandom=499, verbose=TRUE)

test_CWE_12S$cr # observed cumulative species richness

test_CWE_12S$mr # observed mean value of species richness over the n samples

rarefaction_CWE_12S <-rarc(RFA_12S_CWE_noID, nrandom=499)

rarefaction_CWE_12S$out

RF_12S_CWE <- rarefaction_CWE_12S$out [c(1,6)]

RF_12S_CWE$Lake <- c('CWE')

RF_12S_CWE$Locus  <- c('sRNA')




RFA_Cytb_CWE <- RFA_Cytb[which(RFA_Cytb$Lake == 'CWE'),]


RFA_Cytb_CWE_back <- reshape(RFA_Cytb_CWE, idvar= c("SampleID","Lake"), timevar= "Species", direction = "wide")



RFA_Cytb_CWE_back[nrow(RFA_Cytb_CWE_back) + 1,] = c("CWE14","CWE",0,0,0,0)

RFA_Cytb_CWE_back[nrow(RFA_Cytb_CWE_back) + 1,] = c("CWE15","CWE",0,0,0,0) 

RFA_Cytb_CWE_back[nrow(RFA_Cytb_CWE_back) + 1,] = c("CWE19","CWE",0,0,0,0)  


RFA_Cytb_CWE_noID <- RFA_Cytb_CWE_back [c(3:ncol(RFA_Cytb_CWE_back))]

RFA_Cytb_CWE_noID [is.na(RFA_Cytb_CWE_noID)] <- 0

RFA_Cytb_CWE_noID <- mutate_all(RFA_Cytb_CWE_noID, function(x) as.numeric(as.character(x)))




test_CWE_Cytb <- rich(matrix=RFA_Cytb_CWE_noID, nrandom=499, verbose=TRUE)

test_CWE_Cytb$cr # observed cumulative species richness

test_CWE_Cytb$mr # observed mean value of species richness over the n samples

rarefaction_CWE_Cytb <-rarc(RFA_Cytb_CWE_noID, nrandom=499)

rarefaction_CWE_Cytb$out

RF_Cytb_CWE <- rarefaction_CWE_Cytb$out [c(1,6)]

RF_Cytb_CWE$Lake <- c('CWE')

RF_Cytb_CWE$Locus  <- c('Cytb')



RF_CWE <- rbind(RF_Cytb_CWE,RF_12S_CWE)


# OGW ----------------------------------------------------------------------
RFA_12S_OGW <- RFA_12S[which(RFA_12S$Lake == 'OGW'),]


RFA_12S_OGW_back <- reshape(RFA_12S_OGW, idvar= c("SampleID","Lake"), timevar= "Species", direction = "wide")

RFA_12S_OGW_noID <- RFA_12S_OGW_back [c(4:ncol(RFA_12S_OGW_back))]

RFA_12S_OGW_noID [is.na(RFA_12S_OGW_noID)] <- 0



test_OGW_12S <- rich(matrix=RFA_12S_OGW_noID, nrandom=499, verbose=TRUE)

test_OGW_12S$cr # observed cumulative species richness

test_OGW_12S$mr # observed mean value of species richness over the n samples

rarefaction_OGW_12S <-rarc(RFA_12S_OGW_noID, nrandom=499)

rarefaction_OGW_12S$out

RF_12S_OGW <- rarefaction_OGW_12S$out [c(1,6)]

RF_12S_OGW$Lake <- c('OGW')

RF_12S_OGW$Locus  <- c('sRNA')




RFA_Cytb_OGW <- RFA_Cytb[which(RFA_Cytb$Lake == 'OGW'),]


RFA_Cytb_OGW_back <- reshape(RFA_Cytb_OGW, idvar= c("SampleID","Lake"), timevar= "Species", direction = "wide")

RFA_Cytb_OGW_noID <- RFA_Cytb_OGW_back [c(3:ncol(RFA_Cytb_OGW_back))]

RFA_Cytb_OGW_noID [is.na(RFA_Cytb_OGW_noID)] <- 0



test_OGW_Cytb <- rich(matrix=RFA_Cytb_OGW_noID, nrandom=499, verbose=TRUE)

test_OGW_Cytb$cr # observed cumulative species richness

test_OGW_Cytb$mr # observed mean value of species richness over the n samples

rarefaction_OGW_Cytb <-rarc(RFA_Cytb_OGW_noID, nrandom=499)

rarefaction_OGW_Cytb$out

RF_Cytb_OGW <- rarefaction_OGW_Cytb$out [c(1,6)]

RF_Cytb_OGW$Lake <- c('OGW')

RF_Cytb_OGW$Locus  <- c('Cytb')



RF_OGW <- rbind(RF_Cytb_OGW,RF_12S_OGW)




# PAD ----------------------------------------------------------------------
RFA_12S_PAD <- RFA_12S[which(RFA_12S$Lake == 'PAD'),]


RFA_12S_PAD_back <- reshape(RFA_12S_PAD, idvar= c("SampleID","Lake"), timevar= "Species", direction = "wide")

RFA_12S_PAD_noID <- RFA_12S_PAD_back [c(3:ncol(RFA_12S_PAD_back))]

RFA_12S_PAD_noID [is.na(RFA_12S_PAD_noID)] <- 0



test_PAD_12S <- rich(matrix=RFA_12S_PAD_noID, nrandom=499, verbose=TRUE)

test_PAD_12S$cr # observed cumulative species richness

test_PAD_12S$mr # observed mean value of species richness over the n samples

rarefaction_PAD_12S <-rarc(RFA_12S_PAD_noID, nrandom=499)

rarefaction_PAD_12S$out

RF_12S_PAD <- rarefaction_PAD_12S$out [c(1,6)]

RF_12S_PAD$Lake <- c('PAD')

RF_12S_PAD$Locus  <- c('sRNA')




RFA_Cytb_PAD <- RFA_Cytb[which(RFA_Cytb$Lake == 'PAD'),]


RFA_Cytb_PAD_back <- reshape(RFA_Cytb_PAD, idvar= c("SampleID","Lake"), timevar= "Species", direction = "wide")



RFA_Cytb_PAD_noID <- RFA_Cytb_PAD_back [c(3:ncol(RFA_Cytb_PAD_back))]

RFA_Cytb_PAD_noID [is.na(RFA_Cytb_PAD_noID)] <- 0

RFA_Cytb_PAD_noID <- mutate_all(RFA_Cytb_PAD_noID, function(x) as.numeric(as.character(x)))


test_PAD_Cytb <- rich(matrix=RFA_Cytb_PAD_noID, nrandom=499, verbose=TRUE)

test_PAD_Cytb$cr # observed cumulative species richness

test_PAD_Cytb$mr # observed mean value of species richness over the n samples

rarefaction_PAD_Cytb <-rarc(RFA_Cytb_PAD_noID, nrandom=499)

rarefaction_PAD_Cytb$out

RF_Cytb_PAD <- rarefaction_PAD_Cytb$out [c(1,6)]

RF_Cytb_PAD$Lake <- c('PAD')

RF_Cytb_PAD$Locus  <- c('Cytb')



RF_PAD <- rbind(RF_Cytb_PAD,RF_12S_PAD)





# TRA ----------------------------------------------------------------------
RFA_12S_TRA <- RFA_12S[which(RFA_12S$Lake == 'TRA'),]


RFA_12S_TRA_back <- reshape(RFA_12S_TRA, idvar= c("SampleID","Lake"), timevar= "Species", direction = "wide")



RFA_12S_TRA_back[nrow(RFA_12S_TRA_back) + 1,] = c("TRA18","TRA",0,0,0,0,0,0,0)  




RFA_12S_TRA_noID <- RFA_12S_TRA_back [c(3:ncol(RFA_12S_TRA_back))]

RFA_12S_TRA_noID [is.na(RFA_12S_TRA_noID)] <- 0

RFA_12S_TRA_noID <- mutate_all(RFA_12S_TRA_noID, function(x) as.numeric(as.character(x)))




test_TRA_12S <- rich(matrix=RFA_12S_TRA_noID, nrandom=499, verbose=TRUE)

test_TRA_12S$cr # observed cumulative species richness

test_TRA_12S$mr # observed mean value of species richness over the n samples

rarefaction_TRA_12S <-rarc(RFA_12S_TRA_noID, nrandom=499)

rarefaction_TRA_12S$out

RF_12S_TRA <- rarefaction_TRA_12S$out [c(1,6)]

RF_12S_TRA$Lake <- c('TRA')

RF_12S_TRA$Locus  <- c('sRNA')




RFA_Cytb_TRA <- RFA_Cytb[which(RFA_Cytb$Lake == 'TRA'),]


RFA_Cytb_TRA_back <- reshape(RFA_Cytb_TRA, idvar= c("SampleID","Lake"), timevar= "Species", direction = "wide")




RFA_Cytb_TRA_noID <- RFA_Cytb_TRA_back [c(3:ncol(RFA_Cytb_TRA_back))]

RFA_Cytb_TRA_noID [is.na(RFA_Cytb_TRA_noID)] <- 0

RFA_Cytb_TRA_noID <- mutate_all(RFA_Cytb_TRA_noID, function(x) as.numeric(as.character(x)))


test_TRA_Cytb <- rich(matrix=RFA_Cytb_TRA_noID, nrandom=499, verbose=TRUE)

test_TRA_Cytb$cr # observed cumulative species richness

test_TRA_Cytb$mr # observed mean value of species richness over the n samples

rarefaction_TRA_Cytb <-rarc(RFA_Cytb_TRA_noID, nrandom=499)

rarefaction_TRA_Cytb$out

RF_Cytb_TRA <- rarefaction_TRA_Cytb$out [c(1,6)]

RF_Cytb_TRA$Lake <- c('TRA')

RF_Cytb_TRA$Locus  <- c('Cytb')


RF_TRA <- rbind(RF_Cytb_TRA,RF_12S_TRA)





# PEN ----------------------------------------------------------------------

RFA_12S_PEN <- RFA_12S[which(RFA_12S$Lake == 'PEN'),]


RFA_12S_PEN_back <- reshape(RFA_12S_PEN, idvar= c("SampleID","Lake"), timevar= "Species", direction = "wide")



RFA_12S_PEN_noID <- RFA_12S_PEN_back [c(3:ncol(RFA_12S_PEN_back))]

RFA_12S_PEN_noID [is.na(RFA_12S_PEN_noID)] <- 0

RFA_12S_PEN_noID <- mutate_all(RFA_12S_PEN_noID, function(x) as.numeric(as.character(x)))




test_PEN_12S <- rich(matrix=RFA_12S_PEN_noID, nrandom=499, verbose=TRUE)

test_PEN_12S$cr # observed cumulative species richness

test_PEN_12S$mr # observed mean value of species richness over the n samples

rarefaction_PEN_12S <-rarc(RFA_12S_PEN_noID, nrandom=499)

rarefaction_PEN_12S$out

RF_12S_PEN <- rarefaction_PEN_12S$out [c(1,6)]

RF_12S_PEN$Lake <- c('PEN')

RF_12S_PEN$Locus  <- c('sRNA')



RFA_Cytb_PEN <- RFA_Cytb[which(RFA_Cytb$Lake == 'PEN'),]


RFA_Cytb_PEN_back <- reshape(RFA_Cytb_PEN, idvar= c("SampleID","Lake"), timevar= "Species", direction = "wide")




RFA_Cytb_PEN_noID <- RFA_Cytb_PEN_back [c(3:ncol(RFA_Cytb_PEN_back))]

RFA_Cytb_PEN_noID [is.na(RFA_Cytb_PEN_noID)] <- 0

RFA_Cytb_PEN_noID <- mutate_all(RFA_Cytb_PEN_noID, function(x) as.numeric(as.character(x)))


test_PEN_Cytb <- rich(matrix=RFA_Cytb_PEN_noID, nrandom=499, verbose=TRUE)

test_PEN_Cytb$cr # observed cumulative species richness

test_PEN_Cytb$mr # observed mean value of species richness over the n samples

rarefaction_PEN_Cytb <-rarc(RFA_Cytb_PEN_noID, nrandom=499)

rarefaction_PEN_Cytb$out

RF_Cytb_PEN <- rarefaction_PEN_Cytb$out [c(1,6)]

RF_Cytb_PEN$Lake <- c('PEN')

RF_Cytb_PEN$Locus  <- c('Cytb')


RF_PEN <- rbind(RF_Cytb_PEN,RF_12S_PEN)




# KEN ----------------------------------------------------------------------

RFA_12S_KEN <- RFA_12S[which(RFA_12S$Lake == 'KEN'),]


RFA_12S_KEN_back <- reshape(RFA_12S_KEN, idvar= c("SampleID","Lake"), timevar= "Species", direction = "wide")



RFA_12S_KEN_noID <- RFA_12S_KEN_back [c(3:ncol(RFA_12S_KEN_back))]

RFA_12S_KEN_noID [is.na(RFA_12S_KEN_noID)] <- 0

RFA_12S_KEN_noID <- mutate_all(RFA_12S_KEN_noID, function(x) as.numeric(as.character(x)))




test_KEN_12S <- rich(matrix=RFA_12S_KEN_noID, nrandom=499, verbose=TRUE)

test_KEN_12S$cr # observed cumulative species richness

test_KEN_12S$mr # observed mean value of species richness over the n samples

rarefaction_KEN_12S <-rarc(RFA_12S_KEN_noID, nrandom=499)

rarefaction_KEN_12S$out

RF_12S_KEN <- rarefaction_KEN_12S$out [c(1,6)]

RF_12S_KEN$Lake <- c('KEN')

RF_12S_KEN$Locus  <- c('sRNA')



RFA_Cytb_KEN <- RFA_Cytb[which(RFA_Cytb$Lake == 'KEN'),]


RFA_Cytb_KEN_back <- reshape(RFA_Cytb_KEN, idvar= c("SampleID","Lake"), timevar= "Species", direction = "wide")




RFA_Cytb_KEN_noID <- RFA_Cytb_KEN_back [c(3:ncol(RFA_Cytb_KEN_back))]

RFA_Cytb_KEN_noID [is.na(RFA_Cytb_KEN_noID)] <- 0

RFA_Cytb_KEN_noID <- mutate_all(RFA_Cytb_KEN_noID, function(x) as.numeric(as.character(x)))


test_KEN_Cytb <- rich(matrix=RFA_Cytb_KEN_noID, nrandom=499, verbose=TRUE)

test_KEN_Cytb$cr # observed cumulative species richness

test_KEN_Cytb$mr # observed mean value of species richness over the n samples

rarefaction_KEN_Cytb <-rarc(RFA_Cytb_KEN_noID, nrandom=499)

rarefaction_KEN_Cytb$out

RF_Cytb_KEN <- rarefaction_KEN_Cytb$out [c(1,6)]

RF_Cytb_KEN$Lake <- c('KEN')

RF_Cytb_KEN$Locus  <- c('Cytb')


RF_KEN <- rbind(RF_Cytb_KEN,RF_12S_KEN)



# LLB ----------------------------------------------------------------------

RFA_12S_LLB <- RFA_12S[which(RFA_12S$Lake == 'LLB'),]


RFA_12S_LLB_back <- reshape(RFA_12S_LLB, idvar= c("SampleID","Lake"), timevar= "Species", direction = "wide")



RFA_12S_LLB_noID <- RFA_12S_LLB_back [c(3:ncol(RFA_12S_LLB_back))]

RFA_12S_LLB_noID [is.na(RFA_12S_LLB_noID)] <- 0

RFA_12S_LLB_noID <- mutate_all(RFA_12S_LLB_noID, function(x) as.numeric(as.character(x)))




test_LLB_12S <- rich(matrix=RFA_12S_LLB_noID, nrandom=499, verbose=TRUE)

test_LLB_12S$cr # observed cumulative species richness

test_LLB_12S$mr # observed mean value of species richness over the n samples

rarefaction_LLB_12S <-rarc(RFA_12S_LLB_noID, nrandom=499)

rarefaction_LLB_12S$out

RF_12S_LLB <- rarefaction_LLB_12S$out [c(1,6)]

RF_12S_LLB$Lake <- c('LLB')

RF_12S_LLB$Locus  <- c('sRNA')



RFA_Cytb_LLB <- RFA_Cytb[which(RFA_Cytb$Lake == 'LLB'),]


RFA_Cytb_LLB_back <- reshape(RFA_Cytb_LLB, idvar= c("SampleID","Lake"), timevar= "Species", direction = "wide")




RFA_Cytb_LLB_noID <- RFA_Cytb_LLB_back [c(3:ncol(RFA_Cytb_LLB_back))]

RFA_Cytb_LLB_noID [is.na(RFA_Cytb_LLB_noID)] <- 0

RFA_Cytb_LLB_noID <- mutate_all(RFA_Cytb_LLB_noID, function(x) as.numeric(as.character(x)))


test_LLB_Cytb <- rich(matrix=RFA_Cytb_LLB_noID, nrandom=499, verbose=TRUE)

test_LLB_Cytb$cr # observed cumulative species richness

test_LLB_Cytb$mr # observed mean value of species richness over the n samples

rarefaction_LLB_Cytb <-rarc(RFA_Cytb_LLB_noID, nrandom=499)

rarefaction_LLB_Cytb$out

RF_Cytb_LLB <- rarefaction_LLB_Cytb$out [c(1,6)]

RF_Cytb_LLB$Lake <- c('LLB')

RF_Cytb_LLB$Locus  <- c('Cytb')


RF_LLB <- rbind(RF_Cytb_LLB,RF_12S_LLB)





# LLG ----------------------------------------------------------------------

RFA_12S_LLG <- RFA_12S[which(RFA_12S$Lake == 'LLG'),]


RFA_12S_LLG_back <- reshape(RFA_12S_LLG, idvar= c("SampleID","Lake"), timevar= "Species", direction = "wide")



RFA_12S_LLG_noID <- RFA_12S_LLG_back [c(3:ncol(RFA_12S_LLG_back))]

RFA_12S_LLG_noID [is.na(RFA_12S_LLG_noID)] <- 0

RFA_12S_LLG_noID <- mutate_all(RFA_12S_LLG_noID, function(x) as.numeric(as.character(x)))




test_LLG_12S <- rich(matrix=RFA_12S_LLG_noID, nrandom=499, verbose=TRUE)

test_LLG_12S$cr # observed cumulative species richness

test_LLG_12S$mr # observed mean value of species richness over the n samples

rarefaction_LLG_12S <-rarc(RFA_12S_LLG_noID, nrandom=499)

rarefaction_LLG_12S$out

RF_12S_LLG <- rarefaction_LLG_12S$out [c(1,6)]

RF_12S_LLG$Lake <- c('LLG')

RF_12S_LLG$Locus  <- c('sRNA')



RFA_Cytb_LLG <- RFA_Cytb[which(RFA_Cytb$Lake == 'LLG'),]


RFA_Cytb_LLG_back <- reshape(RFA_Cytb_LLG, idvar= c("SampleID","Lake"), timevar= "Species", direction = "wide")




RFA_Cytb_LLG_noID <- RFA_Cytb_LLG_back [c(3:ncol(RFA_Cytb_LLG_back))]

RFA_Cytb_LLG_noID [is.na(RFA_Cytb_LLG_noID)] <- 0

RFA_Cytb_LLG_noID <- mutate_all(RFA_Cytb_LLG_noID, function(x) as.numeric(as.character(x)))


test_LLG_Cytb <- rich(matrix=RFA_Cytb_LLG_noID, nrandom=499, verbose=TRUE)

test_LLG_Cytb$cr # observed cumulative species richness

test_LLG_Cytb$mr # observed mean value of species richness over the n samples

rarefaction_LLG_Cytb <-rarc(RFA_Cytb_LLG_noID, nrandom=499)

rarefaction_LLG_Cytb$out

RF_Cytb_LLG <- rarefaction_LLG_Cytb$out [c(1,6)]

RF_Cytb_LLG$Lake <- c('LLG')

RF_Cytb_LLG$Locus  <- c('Cytb')


RF_LLG <- rbind(RF_Cytb_LLG,RF_12S_LLG)





# MAP ----------------------------------------------------------------------

RFA_12S_MAP <- RFA_12S[which(RFA_12S$Lake == 'MAP'),]


RFA_12S_MAP_back <- reshape(RFA_12S_MAP, idvar= c("SampleID","Lake"), timevar= "Species", direction = "wide")

RFA_12S_MAP_back[nrow(RFA_12S_MAP_back) + 1,] = c("MAP16","MAP",0,0,0,0)  


RFA_12S_MAP_noID <- RFA_12S_MAP_back [c(3:ncol(RFA_12S_MAP_back))]

RFA_12S_MAP_noID [is.na(RFA_12S_MAP_noID)] <- 0

RFA_12S_MAP_noID <- mutate_all(RFA_12S_MAP_noID, function(x) as.numeric(as.character(x)))




test_MAP_12S <- rich(matrix=RFA_12S_MAP_noID, nrandom=499, verbose=TRUE)

test_MAP_12S$cr # observed cumulative species richness

test_MAP_12S$mr # observed mean value of species richness over the n samples

rarefaction_MAP_12S <-rarc(RFA_12S_MAP_noID, nrandom=499)

rarefaction_MAP_12S$out

RF_12S_MAP <- rarefaction_MAP_12S$out [c(1,6)]

RF_12S_MAP$Lake <- c('MAP')

RF_12S_MAP$Locus  <- c('sRNA')



RFA_Cytb_MAP <- RFA_Cytb[which(RFA_Cytb$Lake == 'MAP'),]


RFA_Cytb_MAP_back <- reshape(RFA_Cytb_MAP, idvar= c("SampleID","Lake"), timevar= "Species", direction = "wide")



RFA_Cytb_MAP_noID <- RFA_Cytb_MAP_back [c(3:ncol(RFA_Cytb_MAP_back))]

RFA_Cytb_MAP_noID [is.na(RFA_Cytb_MAP_noID)] <- 0

RFA_Cytb_MAP_noID <- mutate_all(RFA_Cytb_MAP_noID, function(x) as.numeric(as.character(x)))


test_MAP_Cytb <- rich(matrix=RFA_Cytb_MAP_noID, nrandom=499, verbose=TRUE)

test_MAP_Cytb$cr # observed cumulative species richness

test_MAP_Cytb$mr # observed mean value of species richness over the n samples

rarefaction_MAP_Cytb <-rarc(RFA_Cytb_MAP_noID, nrandom=499)

rarefaction_MAP_Cytb$out

RF_Cytb_MAP <- rarefaction_MAP_Cytb$out [c(1,6)]

RF_Cytb_MAP$Lake <- c('MAP')

RF_Cytb_MAP$Locus  <- c('Cytb')


RF_MAP <- rbind(RF_Cytb_MAP,RF_12S_MAP)



# CAM ----------------------------------------------------------------------

RFA_12S_CAM <- RFA_12S[which(RFA_12S$Lake == 'CAM'),]


RFA_12S_CAM_back <- reshape(RFA_12S_CAM, idvar= c("SampleID","Lake"), timevar= "Species", direction = "wide")

RFA_12S_CAM_back[nrow(RFA_12S_CAM_back) + 1,] = c("CAM09","CAM",0,0,0,0)  


RFA_12S_CAM_noID <- RFA_12S_CAM_back [c(3:ncol(RFA_12S_CAM_back))]

RFA_12S_CAM_noID [is.na(RFA_12S_CAM_noID)] <- 0

RFA_12S_CAM_noID <- mutate_all(RFA_12S_CAM_noID, function(x) as.numeric(as.character(x)))




test_CAM_12S <- rich(matrix=RFA_12S_CAM_noID, nrandom=499, verbose=TRUE)

test_CAM_12S$cr # observed cumulative species richness

test_CAM_12S$mr # observed mean value of species richness over the n samples

rarefaction_CAM_12S <-rarc(RFA_12S_CAM_noID, nrandom=499)

rarefaction_CAM_12S$out

RF_12S_CAM <- rarefaction_CAM_12S$out [c(1,6)]

RF_12S_CAM$Lake <- c('CAM')

RF_12S_CAM$Locus  <- c('sRNA')



RFA_Cytb_CAM <- RFA_Cytb[which(RFA_Cytb$Lake == 'CAM'),]


RFA_Cytb_CAM_back <- reshape(RFA_Cytb_CAM, idvar= c("SampleID","Lake"), timevar= "Species", direction = "wide")

RFA_Cytb_CAM_back[nrow(RFA_Cytb_CAM_back) + 1,] = c("CAM07","CAM",0,0,0,0,0)  
RFA_Cytb_CAM_back[nrow(RFA_Cytb_CAM_back) + 1,] = c("CAM09","CAM",0,0,0,0,0)  



RFA_Cytb_CAM_noID <- RFA_Cytb_CAM_back [c(3:ncol(RFA_Cytb_CAM_back))]

RFA_Cytb_CAM_noID [is.na(RFA_Cytb_CAM_noID)] <- 0

RFA_Cytb_CAM_noID <- mutate_all(RFA_Cytb_CAM_noID, function(x) as.numeric(as.character(x)))


test_CAM_Cytb <- rich(matrix=RFA_Cytb_CAM_noID, nrandom=499, verbose=TRUE)

test_CAM_Cytb$cr # observed cumulative species richness

test_CAM_Cytb$mr # observed mean value of species richness over the n samples

rarefaction_CAM_Cytb <-rarc(RFA_Cytb_CAM_noID, nrandom=499)

rarefaction_CAM_Cytb$out

RF_Cytb_CAM <- rarefaction_CAM_Cytb$out [c(1,6)]

RF_Cytb_CAM$Lake <- c('CAM')

RF_Cytb_CAM$Locus  <- c('Cytb')


RF_CAM <- rbind(RF_Cytb_CAM,RF_12S_CAM)





# OSS ----------------------------------------------------------------------

RFA_12S_OSS <- RFA_12S[which(RFA_12S$Lake == 'OSS'),]


RFA_12S_OSS_back <- reshape(RFA_12S_OSS, idvar= c("SampleID","Lake"), timevar= "Species", direction = "wide")



RFA_12S_OSS_noID <- RFA_12S_OSS_back [c(3:ncol(RFA_12S_OSS_back))]

RFA_12S_OSS_noID [is.na(RFA_12S_OSS_noID)] <- 0

RFA_12S_OSS_noID <- mutate_all(RFA_12S_OSS_noID, function(x) as.numeric(as.character(x)))




test_OSS_12S <- rich(matrix=RFA_12S_OSS_noID, nrandom=499, verbose=TRUE)

test_OSS_12S$cr # observed cumulative species richness

test_OSS_12S$mr # observed mean value of species richness over the n samples

rarefaction_OSS_12S <-rarc(RFA_12S_OSS_noID, nrandom=499)

rarefaction_OSS_12S$out

RF_12S_OSS <- rarefaction_OSS_12S$out [c(1,6)]

RF_12S_OSS$Lake <- c('OSS')

RF_12S_OSS$Locus  <- c('sRNA')



RFA_Cytb_OSS <- RFA_Cytb[which(RFA_Cytb$Lake == 'OSS'),]


RFA_Cytb_OSS_back <- reshape(RFA_Cytb_OSS, idvar= c("SampleID","Lake"), timevar= "Species", direction = "wide")



RFA_Cytb_OSS_noID <- RFA_Cytb_OSS_back [c(3:ncol(RFA_Cytb_OSS_back))]

RFA_Cytb_OSS_noID [is.na(RFA_Cytb_OSS_noID)] <- 0

RFA_Cytb_OSS_noID <- mutate_all(RFA_Cytb_OSS_noID, function(x) as.numeric(as.character(x)))


test_OSS_Cytb <- rich(matrix=RFA_Cytb_OSS_noID, nrandom=499, verbose=TRUE)

test_OSS_Cytb$cr # observed cumulative species richness

test_OSS_Cytb$mr # observed mean value of species richness over the n samples

rarefaction_OSS_Cytb <-rarc(RFA_Cytb_OSS_noID, nrandom=499)

rarefaction_OSS_Cytb$out

RF_Cytb_OSS <- rarefaction_OSS_Cytb$out [c(1,6)]

RF_Cytb_OSS$Lake <- c('OSS')

RF_Cytb_OSS$Locus  <- c('Cytb')


RF_OSS <- rbind(RF_Cytb_OSS,RF_12S_OSS)




# FEN ----------------------------------------------------------------------

RFA_12S_FEN <- RFA_12S[which(RFA_12S$Lake == 'FEN'),]


RFA_12S_FEN_back <- reshape(RFA_12S_FEN, idvar= c("SampleID","Lake"), timevar= "Species", direction = "wide")



RFA_12S_FEN_noID <- RFA_12S_FEN_back [c(3:ncol(RFA_12S_FEN_back))]

RFA_12S_FEN_noID [is.na(RFA_12S_FEN_noID)] <- 0

RFA_12S_FEN_noID <- mutate_all(RFA_12S_FEN_noID, function(x) as.numeric(as.character(x)))




test_FEN_12S <- rich(matrix=RFA_12S_FEN_noID, nrandom=499, verbose=TRUE)

test_FEN_12S$cr # observed cumulative species richness

test_FEN_12S$mr # observed mean value of species richness over the n samples

rarefaction_FEN_12S <-rarc(RFA_12S_FEN_noID, nrandom=499)

rarefaction_FEN_12S$out

RF_12S_FEN <- rarefaction_FEN_12S$out [c(1,6)]

RF_12S_FEN$Lake <- c('FEN')

RF_12S_FEN$Locus  <- c('sRNA')



RFA_Cytb_FEN <- RFA_Cytb[which(RFA_Cytb$Lake == 'FEN'),]


RFA_Cytb_FEN_back <- reshape(RFA_Cytb_FEN, idvar= c("SampleID","Lake"), timevar= "Species", direction = "wide")



RFA_Cytb_FEN_noID <- RFA_Cytb_FEN_back [c(3:ncol(RFA_Cytb_FEN_back))]

RFA_Cytb_FEN_noID [is.na(RFA_Cytb_FEN_noID)] <- 0

RFA_Cytb_FEN_noID <- mutate_all(RFA_Cytb_FEN_noID, function(x) as.numeric(as.character(x)))


test_FEN_Cytb <- rich(matrix=RFA_Cytb_FEN_noID, nrandom=499, verbose=TRUE)

test_FEN_Cytb$cr # observed cumulative species richness

test_FEN_Cytb$mr # observed mean value of species richness over the n samples

rarefaction_FEN_Cytb <-rarc(RFA_Cytb_FEN_noID, nrandom=499)

rarefaction_FEN_Cytb$out

RF_Cytb_FEN <- rarefaction_FEN_Cytb$out [c(1,6)]

RF_Cytb_FEN$Lake <- c('FEN')

RF_Cytb_FEN$Locus  <- c('Cytb')


RF_FEN <- rbind(RF_Cytb_FEN,RF_12S_FEN)





# WLF ----------------------------------------------------------------------

RFA_12S_WLF <- RFA_12S[which(RFA_12S$Lake == 'WLF'),]


RFA_12S_WLF_back <- reshape(RFA_12S_WLF, idvar= c("SampleID","Lake"), timevar= "Species", direction = "wide")



RFA_12S_WLF_noID <- RFA_12S_WLF_back [c(3:ncol(RFA_12S_WLF_back))]

RFA_12S_WLF_noID [is.na(RFA_12S_WLF_noID)] <- 0

RFA_12S_WLF_noID <- mutate_all(RFA_12S_WLF_noID, function(x) as.numeric(as.character(x)))




test_WLF_12S <- rich(matrix=RFA_12S_WLF_noID, nrandom=499, verbose=TRUE)

test_WLF_12S$cr # observed cumulative species richness

test_WLF_12S$mr # observed mean value of species richness over the n samples

rarefaction_WLF_12S <-rarc(RFA_12S_WLF_noID, nrandom=499)

rarefaction_WLF_12S$out

RF_12S_WLF <- rarefaction_WLF_12S$out [c(1,6)]

RF_12S_WLF$Lake <- c('WLF')

RF_12S_WLF$Locus  <- c('sRNA')



RFA_Cytb_WLF <- RFA_Cytb[which(RFA_Cytb$Lake == 'WLF'),]


RFA_Cytb_WLF_back <- reshape(RFA_Cytb_WLF, idvar= c("SampleID","Lake"), timevar= "Species", direction = "wide")



RFA_Cytb_WLF_back[nrow(RFA_Cytb_WLF_back) + 1,] = c("WLF11","WLF",0,0,0,0,0,0,0,0,0)  




RFA_Cytb_WLF_noID <- RFA_Cytb_WLF_back [c(3:ncol(RFA_Cytb_WLF_back))]

RFA_Cytb_WLF_noID [is.na(RFA_Cytb_WLF_noID)] <- 0

RFA_Cytb_WLF_noID <- mutate_all(RFA_Cytb_WLF_noID, function(x) as.numeric(as.character(x)))


test_WLF_Cytb <- rich(matrix=RFA_Cytb_WLF_noID, nrandom=499, verbose=TRUE)

test_WLF_Cytb$cr # observed cumulative species richness

test_WLF_Cytb$mr # observed mean value of species richness over the n samples

rarefaction_WLF_Cytb <-rarc(RFA_Cytb_WLF_noID, nrandom=499)

rarefaction_WLF_Cytb$out

RF_Cytb_WLF <- rarefaction_WLF_Cytb$out [c(1,6)]

RF_Cytb_WLF$Lake <- c('WLF')

RF_Cytb_WLF$Locus  <- c('Cytb')


RF_WLF <- rbind(RF_Cytb_WLF,RF_12S_WLF)







# BET ----------------------------------------------------------------------

RFA_12S_BET <- RFA_12S[which(RFA_12S$Lake == 'BET'),]


RFA_12S_BET_back <- reshape(RFA_12S_BET, idvar= c("SampleID","Lake"), timevar= "Species", direction = "wide")

RFA_12S_BET_back[nrow(RFA_12S_BET_back) + 1,] = c("BET13","BET",0,0,0,0,0,0,0,0)  



RFA_12S_BET_noID <- RFA_12S_BET_back [c(3:ncol(RFA_12S_BET_back))]

RFA_12S_BET_noID [is.na(RFA_12S_BET_noID)] <- 0

RFA_12S_BET_noID <- mutate_all(RFA_12S_BET_noID, function(x) as.numeric(as.character(x)))




test_BET_12S <- rich(matrix=RFA_12S_BET_noID, nrandom=499, verbose=TRUE)

test_BET_12S$cr # observed cumulative species richness

test_BET_12S$mr # observed mean value of species richness over the n samples

rarefaction_BET_12S <-rarc(RFA_12S_BET_noID, nrandom=499)

rarefaction_BET_12S$out

RF_12S_BET <- rarefaction_BET_12S$out [c(1,6)]

RF_12S_BET$Lake <- c('BET')

RF_12S_BET$Locus  <- c('sRNA')



RFA_Cytb_BET <- RFA_Cytb[which(RFA_Cytb$Lake == 'BET'),]


RFA_Cytb_BET_back <- reshape(RFA_Cytb_BET, idvar= c("SampleID","Lake"), timevar= "Species", direction = "wide")



RFA_Cytb_BET_back[nrow(RFA_Cytb_BET_back) + 1,] = c("BET20","BET",0,0,0,0,0,0,0,0)  




RFA_Cytb_BET_noID <- RFA_Cytb_BET_back [c(3:ncol(RFA_Cytb_BET_back))]

RFA_Cytb_BET_noID [is.na(RFA_Cytb_BET_noID)] <- 0

RFA_Cytb_BET_noID <- mutate_all(RFA_Cytb_BET_noID, function(x) as.numeric(as.character(x)))


test_BET_Cytb <- rich(matrix=RFA_Cytb_BET_noID, nrandom=499, verbose=TRUE)

test_BET_Cytb$cr # observed cumulative species richness

test_BET_Cytb$mr # observed mean value of species richness over the n samples

rarefaction_BET_Cytb <-rarc(RFA_Cytb_BET_noID, nrandom=499)

rarefaction_BET_Cytb$out

RF_Cytb_BET <- rarefaction_BET_Cytb$out [c(1,6)]

RF_Cytb_BET$Lake <- c('BET')

RF_Cytb_BET$Locus  <- c('Cytb')


RF_BET <- rbind(RF_Cytb_BET,RF_12S_BET)


####RFA_FIG#####

RF <- rbind(RF_CWE,RF_OGW,RF_PAD,RF_TRA,RF_PEN,RF_KEN,RF_LLB,RF_LLG,RF_MAP,RF_CAM,RF_OSS,RF_FEN,RF_WLF,RF_BET)

RF$Lake <- factor(RF$Lake, levels = (c("CWE","PAD","OGW","PEN","TRA","KEN","LLB","LLG","MAP","CAM","OSS","FEN","WLF","BET")))



RF_line <- data.frame(Locus = c("Cytb","sRNA","Cytb","sRNA"), sample = c(7, 7,11,11))


RF_label <- data.frame(x=c(0,0),y=c(9,9),Locus = c("Cytb","sRNA"), labs = c("(a)","(b)"))





ggplot(RF, aes(x= sample, y= mean.richness,colour= Lake))+
geom_point(stat="identity",size=3)+
geom_text(data=RF_label, aes(x=x, y=y,label=labs), size=5,colour="black")+
facet_wrap(~Locus)+
scale_x_continuous(limits = c(0, 20), breaks=seq(0, 20, 2))+
scale_y_continuous(limits = c(0, 9), breaks=seq(0, 9, 1))+
labs(x="Number of shore samples", y="Number of species")+
scale_colour_manual(name="Sampling lake",values=colpalette) +
guides(col = guide_legend(ncol = 7, byrow = TRUE))+
theme_bw()+ geom_vline(aes(xintercept = sample), data = RF_line,linetype="dashed",size=1)+
theme(text=element_text(size=15),legend.position="bottom",legend.title = element_text(size=12),legend.text = element_text(size=12),
      panel.grid.minor=element_blank(),panel.grid.major =element_blank(),strip.text=element_blank())


##Fig5.3_RFA_SL####

ggsave("FigS5.3_RFA_SL.jpeg", path = "Lake project/Figures_new/",  width = 8, height = 6, units = "in", dpi=500)



####Transect#####

LP_rs_wlake_SO_tran = data.frame(Lake=factor(LP_rs_wlake_SO$Lake),Project=factor(LP_rs_wlake_SO$Project),Transect=factor(LP_rs_wlake_SO$Transect),
                             Locus=factor(LP_rs_wlake_SO$Locus),Species=factor(LP_rs_wlake_SO$Species))



wlake_tran_site =dcast(LP_rs_wlake_SO_tran, Project+Lake+Species~Locus+Transect)




wlake_tran_site_sub <- wlake_tran_site[which(wlake_tran_site$Lake == "Llyn_Cwellyn"|
                                             wlake_tran_site$Lake == "Llyn_Padarn"|
                                             wlake_tran_site$Lake == "Llangorse_Lake" ),]

Wlake_Expected_EDNA[which(Wlake_Expected_EDNA$Presence == -1),]

#excduding the species records based on the above results




wlake_tran_site_sub  <- wlake_tran_site_sub[which(wlake_tran_site_sub$Lake !="Llyn_Cwellyn" |
                                                   wlake_tran_site_sub$Species != "Rutilus_rutilus"),]



wlake_tran_site_sub <- wlake_tran_site_sub[which(wlake_tran_site_sub$Lake !="Llyn_Padarn" |
                                                    wlake_tran_site_sub$Species != "Salmo"),]



LP_rs_mere_SO_tran = data.frame(Lake=factor(LP_rs_mere_SO$Lake),Project=factor(LP_rs_mere_SO$Project),Transect=factor(LP_rs_mere_SO$Transect),
                                 Locus=factor(LP_rs_mere_SO$Locus),Species=factor(LP_rs_mere_SO$Species))



mere_tran_site =dcast(LP_rs_mere_SO_tran, Project+Lake+Species~Locus+Transect)


#excduding the species records based on the above results

mere_Expected_EDNA[which(mere_Expected_EDNA$Presence == -1),]

mere_tran_site <- mere_tran_site[which(mere_tran_site$Lake !="Betley_Mere" |
                                         mere_tran_site$Species != "Gobio_gobio"),]



LP_tran_site <- rbind (wlake_tran_site_sub,mere_tran_site)


LP_tran_site$Lake <-mapvalues(LP_tran_site$Lake, c("Betley_Mere","Chapel_Mere","Fenemere","Llangorse_Lake",  
                                                  "Llyn_Cwellyn","Llyn_Padarn","Maer_Pool","Oss_Mere","Watch_Lane_Flash"), 
                                                  c("BET","CAM","FEN","LLG","CWE","PAD","MAP","OSS","WLF"))



LP_tran_site$Lake <- factor(LP_tran_site$Lake, levels = (c("CWE","PAD","LLG","MAP","CAM","OSS","FEN","WLF","BET")))


colpalette4 <-c("#8cd6b5",
               "#4f2d46",
               "#cfd058",
               "#74d656",
               "#a688c9",
               "#c34980",
               "#5541a7",
               "#4b4d32",
               "#c34837")




tranfig1 <- ggplot(LP_tran_site, aes(x = Cytb_SL/10, y = Cytb_OSL/10, colour=Lake))+
            geom_point(stat="identity",size=3,position = position_jitterdodge(dodge.width = 0.02))+
            scale_y_continuous(limits = c(0, 1.0), breaks=seq(0, 1.0, 0.2))+
            scale_x_continuous(limits = c(0, 1.0), breaks=seq(0, 1.0, 0.2))+
            scale_colour_manual(name="Sampling site",values=colpalette4) +
            geom_smooth(formula= y ~ x,method="lm", se=FALSE,lwd=1.0,na.rm=TRUE)+
            guides(col = guide_legend(ncol = 9, byrow = TRUE))+
            labs(x="Site occupancy of shore samples (Cytb)", y="Site occupancy of offshore samples (Cytb)")+theme_bw()+
            stat_cor(method = "pearson",label.sep = "; ",label.y.npc="centre",label.x=0.75,size=3, hjust = 0, show.legend = FALSE)+
            theme(text=element_text(size=18), legend.position = "bottom",legend.title = element_text(size=15),legend.text = element_text(size=15),
                  panel.grid.minor=element_blank(),panel.grid.major =element_blank())




##Pearson_test
corfun<-function(x, y) {
  corr=(cor.test(x, y,
                 alternative="two.sided", method="pearson"))}

ddply(LP_tran_site, .(Lake), summarise, s=corfun(Cytb_SL,Cytb_OSL)$statistic,
      cor.est=corfun(Cytb_SL,Cytb_OSL)$estimate,
      pval=corfun(Cytb_SL,Cytb_OSL)$p.value,
      alt=corfun(Cytb_SL,Cytb_OSL)$alternative,
      df=corfun(Cytb_SL,Cytb_OSL)$parameter)

ddply(LP_tran_site, .(Lake), summarise, s=corfun(`12S_SL`,`12S_OSL`)$statistic,
      cor.est=corfun(`12S_SL`,`12S_OSL`)$estimate,
      pval=corfun(`12S_SL`,`12S_OSL`)$p.value,
      alt=corfun(`12S_SL`,`12S_OSL`)$alternative,
      df=corfun(`12S_SL`,`12S_OSL`)$parameter)


legend_tranfig1 <- get_legend(tranfig1)

tranfig1  <- tranfig1  + theme(legend.position="none")+
  geom_text(aes(x=0, y=1.0, label= '(a)'),hjust = 0, size=5, colour = "black")+
  geom_text(aes(x=0.18, y=0.9, label= 'CAR'),hjust = 0, size=5, colour = "black",fontface = "italic")+
  geom_text(aes(x=0.72, y=0.1, label= 'BRE'),hjust = 0, size=5, colour = "black",fontface = "italic")




tranfig2 <- ggplot(LP_tran_site, aes(x =  `12S_SL`/10, y = `12S_OSL`/10, colour=Lake))+
            geom_point(stat="identity",size=3,position = position_jitterdodge(dodge.width = 0.02))+
            scale_y_continuous(limits = c(0, 1.0), breaks=seq(0, 1.0, 0.2))+
            scale_x_continuous(limits = c(0, 1.0), breaks=seq(0, 1.0, 0.2))+
            scale_colour_manual(name="Sampling site",values=colpalette4) +
            geom_smooth(formula= y ~ x,method="lm", se=FALSE,lwd=1.0,na.rm=TRUE)+
            guides(col = guide_legend(ncol = 9, byrow = TRUE))+
            labs(x="Site occupancy of shore samples (12S)", y="Site occupancy of offshore samples (12S)")+theme_bw()+
            stat_cor(method = "pearson",label.sep = "; ",label.y.npc="centre",label.x=0.75,size=3, hjust = 0, show.legend = FALSE)+
            theme(text=element_text(size=18), legend.position = "bottom",legend.title = element_text(size=15),legend.text = element_text(size=15),
                  panel.grid.minor=element_blank(),panel.grid.major =element_blank())




tranfig2  <- tranfig2  + theme(legend.position="none")+
  geom_text(aes(x=0, y=1.0, label= '(b)'),hjust = 0, size=5, colour = "black",parse = TRUE)



TRAN_fig1 <- grid.arrange(tranfig1, tranfig2, ncol = 2, nrow = 1,
                           widths = c(2.7,2.7), heights = c(2.5))



TRAN_fig2 <- grid.arrange(TRAN_fig1, legend_tranfig1, nrow = 2,
                           heights = c(2.5,0.3))



####FigS5.1_Tran####
ggsave("FigS5.1_Tran.jpeg",TRAN_fig2, path = "Lake project/Figures_new/",  width = 11, height = 6.5, units = "in", dpi=500)



LP_tran_site$Species <- mapvalues(LP_tran_site$Species, c("Anguilla_anguilla","Abramis_brama","Carassius_carassius","Carassius_auratus","Cyprinus_carpio", 
                                                          "Leucaspius_delineatus","Phoxinus_phoxinus","Rutilus_rutilus","Scardinius_erythrophthalmus",
                                                          "Tinca_tinca","Esox_lucius","Cottus_gobio","Gasterosteus_aculeatus",     
                                                          "Perca_fluviatilis","Salmo_trutta","Gobio_gobio","Pseudorasbora_parva",        
                                                          "Rhodeus_amarus","Squalius_cephalus","Alburnus_alburnus","Pungitius_pungitius",
                                                          "Oncorhynchus_mykiss","Salmo_salar","Salvelinus_alpinus"), 
                                                        c("EEL","BRE","CRU","GOF","CAR",
                                                          "SUN","MIN","ROA","RUD",
                                                          "TEN","PIK","BUL","3SS",
                                                          "PER","BTR","GUD","TMG",
                                                          "BIT","CHU","BLE","9SS",
                                                          "RTR","SAL","CHA"))





LP_tran_site_Cytb <- LP_tran_site[which(LP_tran_site$Species != "GOF"&
                                         LP_tran_site$Species != "CRU"&
                                         LP_tran_site$Species != "GUD"&
                                         LP_tran_site$Species != "TMG"&
                                         LP_tran_site$Species != "BIT"&
                                         LP_tran_site$Species != "SUN"&
                                         LP_tran_site$Species != "CHU"),]



LP_tran_site_Cytb$Species <- factor(LP_tran_site_Cytb$Species, levels = (c("BRE","EEL","CAR","PIK","3SS","PER","MIN","ROA","SAL",
                                                                         "BTR","CHA","RUD","TEN","BUL")))




Tran_Sp_Cytb <-ggplot(LP_tran_site_Cytb, aes(x = Cytb_SL/10, y = Cytb_OSL/10))+
                geom_point(stat="identity",size=2.5)+
                facet_wrap(~Species,ncol =5)+
                scale_y_continuous(limits = c(0, 1.0), breaks=seq(0, 1.0, 0.2))+
                scale_x_continuous(limits = c(0, 1.0), breaks=seq(0, 1.0, 0.2))+
                geom_smooth(formula= y ~ x,method="lm", se=FALSE,lwd=0.5,na.rm=TRUE, colour='black')+
                stat_cor(method = "pearson",label.sep = "\n",label.x=0.5,label.y=0.2,size=3,hjust = 0)+
                labs(x="Site occupancy of shore samples (Cytb)", y="Site occupancy of offshore samples (Cytb)")+theme_bw()+
                theme(text=element_text(size=12),strip.background = element_blank(),strip.text=element_text(face="bold.italic"),
                      panel.grid.minor=element_blank(),panel.grid.major =element_blank())+
                geom_text_repel(aes(label=Lake),size=3)




Tran_Sp_len <- length(levels(LP_tran_site_Cytb$Species))

Tran_Sp_vars <- data.frame(expand.grid(levels(LP_tran_site_Cytb$Species)))

colnames(Tran_Sp_vars) <- c("Species")


Tran_Sp_dat3 <- data.frame(x = rep(0, Tran_Sp_len), y = rep(0.98, Tran_Sp_len), 
                           Tran_Sp_vars, labs=paste0("(",letters[1:14],1,")"))



Tran_Sp_1<- Tran_Sp_Cytb + geom_text(data=Tran_Sp_dat3, aes(x, y, label=labs),hjust = 0, size=3,fontface ="bold",parse = TRUE)





LP_tran_site_12S <- LP_tran_site[which(LP_tran_site$Species != "GOF"&
                                          LP_tran_site$Species != "CRU"&
                                          LP_tran_site$Species != "SUN"&
                                          LP_tran_site$Species != "TMG"&
                                          LP_tran_site$Species != "BUL"&
                                          LP_tran_site$Species != "CHU"&
                                          LP_tran_site$Species != "BIT"),]



LP_tran_site_12S$Species <- factor(LP_tran_site_12S$Species, levels = (c("BRE","EEL","CAR","PIK","3SS","PER","MIN","ROA","SAL",
                                                                           "BTR","CHA","RUD","TEN","GUD")))





Tran_Sp_sRNA <-ggplot(LP_tran_site_12S, aes(x = `12S_SL`/10, y = `12S_OSL`/10))+
                geom_point(stat="identity",size=2.5)+
                facet_wrap(~Species,ncol =5)+
                scale_y_continuous(limits = c(0, 1.0), breaks=seq(0, 1.0, 0.2))+
                scale_x_continuous(limits = c(0, 1.0), breaks=seq(0, 1.0, 0.2))+
                geom_smooth(formula= y ~ x,method="lm", se=FALSE,lwd=0.5,na.rm=TRUE, colour='black')+
                stat_cor(method = "pearson",label.sep = "\n",label.x=0.5,label.y=0.15,size=3,hjust = 0)+
                labs(x="Site occupancy of shore samples (12S)", y="Site occupancy of offshore samples (12S)")+theme_bw()+
                theme(text=element_text(size=12),strip.background = element_blank(),strip.text=element_text(face="bold.italic"),
                      panel.grid.minor=element_blank(),panel.grid.major =element_blank())+
                geom_text_repel(aes(label=Lake),size=3)




Tran_Sp_len <- length(levels(LP_tran_site_12S$Species))

Tran_Sp_vars <- data.frame(expand.grid(levels(LP_tran_site_12S$Species)))

colnames(Tran_Sp_vars) <- c("Species")


Tran_Sp_dat3 <- data.frame(x = rep(0, Tran_Sp_len), y = rep(0.98, Tran_Sp_len), 
                           Tran_Sp_vars, labs=paste0("(",letters[1:14],2,")"))



Tran_Sp_2<- Tran_Sp_sRNA + geom_text(data=Tran_Sp_dat3, aes(x, y, label=labs),hjust = 0, size=3,fontface ="bold",parse = TRUE)


Tran_Sp <- grid.arrange(Tran_Sp_1,Tran_Sp_2, ncol = 1, nrow = 2,
                          widths = c(2.7), heights = c(2.5,2.5))

####FigS5.2_Sp####

ggsave("FigS5.2_Sp.jpeg",Tran_Sp, path = "Lake project/Figures_new/",  width = 9, height = 11, units = "in", dpi=500)






species_richness <- data.frame(OSL=c(6,6,7,5,5,8,9,10,8),SL=c(6,7,8,5,5,7,9,10,9))


species_richness_rs <- melt(species_richness)

#Rename the variables

species_richness_rs <- rename (species_richness_rs,c("variable"="Location","value"="Number"))



t.test(Number ~ Location, data = species_richness_rs, paired = TRUE)



