# This code prepares the data for qSIP calculations, including calculating tube-level WADs for all OTUs and analyzing the shifts among tubes


#Set working directory and load libraries & scripts:

library(reshape2)
library(VennDiagram)
source("qSIP_repo/sample.vec.R")                      #sample.vec
source("qSIP_repo/WAD.func.R")                        #WAD.func
source("qSIP_repo/fit.norm.func.R")                   #fit.norm.func
source("qSIP_repo/boot.WAD.func.R")                   #boot.WAD.func
source("qSIP_repo/diff.wad.calc.R")                   #diff.wad.calc
source("qSIP_repo/boot.diff.wad.R")                   #boot.diff.wad
source("qSIP_repo/MW.calc.R")                         #MW.calc
source("qSIP_repo/MW.calc.Schildkraut.R")             #MW.calc.Schildkraut
source("qSIP_repo/comparison.message.R")              #comparison.message
source("qSIP_repo/ape.calc.R")                        #ape.calc
source("qSIP_repo/boot.diff.ape.R")                   #boot.diff.ape
source("qSIP_repo/r.calc.R")                          #r.calc
source("qSIP_repo/boot.diff.r.R")                     #boot.diff.r
source("qSIP_repo/boot.TUBE.func.R")                  #boot.TUBE.func
source("qSIP_repo/f.calc.R")                          #f.calc
source("qSIP_repo/boot.diff.f.R")                     #boot.diff.f
source("qSIP_repo/all.taxa.calcs.R")                  #all.taxa.calcs
source("qSIP_repo/id.reps.R")                         #id.reps
source("qSIP_repo/select.rep.R")                      #select.rep
source("qSIP_repo/explore.filter.taxa.R")             #explore.filter.taxa
source("qSIP_repo/filter.taxa.R")                     #filter.taxa
source("qSIP_repo/explore.filter.fractions.taxa.R")   #explore.filter.fractions.taxa
source("qSIP_repo/filter.fractions.taxa.R")           #filter.fractions.taxa  
source("qSIP_repo/WAD.by.taxon.func.R")               #WAD.by.taxon.func
source("qSIP_repo/SE.WAD.by.taxon.plot.R")            #SE.WAD.by.taxon.plot
source("qSIP_repo/find.unlabeled.correction.R")       #find.unlabeled.correction
source("qSIP_repo/find.labeled.correction.R")         #find.labeled.correction
source("qSIP_repo/td.pos.resid.R")                    #td.pos.resid
source("qSIP_repo/td.abs.resid.R")                    #td.abs.resid
source("qSIP_repo/bu.abs.resid.R")                    #bu.abs.resid
source("qSIP_repo/select.best.iteration.R")           #select.best.iteration
source("qSIP_repo/find.labeled.correction.plot.R")    #find.labeled.correction.plot
source("qSIP_repo/get.seq.taxa.nums.R")               #get.seq.taxa.nums
source("qSIP_repo/add.lab.WAD.corr.summary.R")        #add.lab.WAD.corr.summary
source("qSIP_repo/apply.unlabeled.correction.R")      #apply.unlabeled.correction
source("qSIP_repo/apply.labeled.correction.R")        #apply.labeled.correction

data.melted <- read.table("data.melted.csv", header = TRUE, sep = ",")

data.melted$trt.code<- gsub('Angelo_18', 'ang_18O', data.melted$Soil_Isotope)
data.melted$trt.code<- gsub('Angelo_16', 'ang_16O', data.melted$trt.code)

data.melted$trt.code<- gsub('Hopland_18', 'hop_18O', data.melted$trt.code)
data.melted$trt.code<- gsub('Hopland_16', 'hop_16O', data.melted$trt.code)

data.melted$trt.code<- gsub('Sedgwick_18', 'sed_18O', data.melted$trt.code)
data.melted$trt.code<- gsub('Sedgwick_16', 'sed_16O', data.melted$trt.code)

data.melted$Isotope_Treatment<- data.melted$trt.code
data.melted<- data.melted %>%
  select(-trt.code, -trt.code.1)

#Calculate corrected WADs for each tube and taxon:
WAD.by.taxon <- WAD.by.taxon.func(X=data.melted, vars=c("taxon", "Density_g_ml", "copies.ul", "Tube", "Soil_Isotope"))

#Look at the results:
names(WAD.by.taxon)          #names of the two data frames in the output list
head(WAD.by.taxon$obs.wads)  #looking at the head of the first data frame
head(WAD.by.taxon[[1]])      #another way to look at the head of the first data frame
WAD.by.taxon$reps.by.trt     #looking at the head of the second data frame
WAD.by.taxon[[2]]            #another way to look at the head of the second data frame


##### Perform the tube-level WAD corrections for ALL replicates of ALL treatments, and then proceed with filtering and standard qSIP analysis:

#Calculate the shift in WAD for the specified unlabeled replicates from their global mean using the taxa common to ALL replicates of ALL treatments:
#MF NOTE: I used global mean of taxa common to all replicates within each treatment (site) because there are so few taxa common to all treatments, and because the WAD shifts for each site look different
#Identify the appropriate shift for unlabeled WADs:
#Note: this function takes a while to run (i.e., ~_?_min for ~150 taxa common to all unlabeled and all labeled treatments):
#Include all unlabeled treatments and all labeled treatments:
system.time(unlab.WAD.corr.list.ang <- find.unlabeled.correction(LIST=WAD.by.taxon, unlab.tmts=c("Angelo_16"), lab.tmts=c("Angelo_18"), CI=0.90))
system.time(unlab.WAD.corr.list.hop <- find.unlabeled.correction(LIST=WAD.by.taxon, unlab.tmts=c("Hopland_16"), lab.tmts=c("Hopland_18"), CI=0.90))
system.time(unlab.WAD.corr.list.sed <- find.unlabeled.correction(LIST=WAD.by.taxon, unlab.tmts=c("Sedgwick_16"), lab.tmts=c("Sedgwick_18"), CI=0.90))

#Look at the results:
names(unlab.WAD.corr.list)                           #names of the two data frames in the output list
unlab.WAD.corr.list.hop$WAD.norm.fit.parms               #looking at the first object -- a data frame
unlab.WAD.corr.list.ang$WAD.norm.fit.parms
unlab.WAD.corr.list.sed$WAD.norm.fit.parms

unlab.WAD.corr.list.hop[[1]]                             #another way to look at the first object
unlab.WAD.corr.list.hop$corr.names                       #looking at the second object in the list -- a vector of the names of the corrected unlabeled replicates
unlab.WAD.corr.list[[2]]                             #another way to look at the second object
head(unlab.WAD.corr.list.hop$WAD.table.corr)             #looking at the head of the third object -- a data frame
head(unlab.WAD.corr.list[[3]])                       #another way to look at the head of the third object

#####################
#####################
#Calculate the shift in WAD for each of the labeled replicates (3 labeled treatments x 3 replicates = 9 total labeled replicates):
#1_AN_Exu_18O:
system.time(ang.18O.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.list.ang, lab.replicate="R1.Angelo_18", lab.names=c("R1.Angelo_18", "R2.Angelo_18", "R3.Angelo_18"), method="td.pos.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))
system.time(R2.ang.18O.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.list.ang, lab.replicate="R2.Angelo_18", lab.names=c("R1.Angelo_18", "R2.Angelo_18", "R3.Angelo_18"), method="td.pos.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))
system.time(R3.ang.18O.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.list.ang, lab.replicate="R3.Angelo_18", lab.names=c("R1.Angelo_18", "R2.Angelo_18", "R3.Angelo_18"), method="td.pos.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))

system.time(hop.18O.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.list.hop, lab.replicate="R1.Hopland_18", lab.names=c("R1.Hopland_18", "R2.Hopland_18", "R3.Hopland_18"), method="td.pos.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))
system.time(R2.hop.18O.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.list.hop, lab.replicate="R2.Hopland_18", lab.names=c("R1.Hopland_18", "R2.Hopland_18", "R3.Hopland_18"), method="td.pos.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))
system.time(R3.hop.18O.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.list.hop, lab.replicate="R3.Hopland_18", lab.names=c("R1.Hopland_18", "R2.Hopland_18", "R3.Hopland_18"), method="td.pos.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))

system.time(sed.18O.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.list.sed, lab.replicate="R1.Sedgwick_18", lab.names=c("R1.Sedgwick_18", "R2.Sedgwick_18", "R3.Sedgwick_18"), method="td.pos.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))
system.time(R2.sed.18O.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.list.sed, lab.replicate="R2.Sedgwick_18", lab.names=c("R1.Sedgwick_18", "R2.Sedgwick_18", "R3.Sedgwick_18"), method="td.pos.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))
system.time(R3.sed.18O.td.pos.resid.list <- find.labeled.correction(LIST=unlab.WAD.corr.list.sed, lab.replicate="R3.Sedgwick_18", lab.names=c("R1.Sedgwick_18", "R2.Sedgwick_18", "R3.Sedgwick_18"), method="td.pos.resid", unlab.SD.percentile=0.5, lab.SD.percentile=0.5, min.num.nongrowers=10))

#########################
########################
#Find the best iteration:
R1.ang.18O.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=ang.18O.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
find.labeled.correction.plot(find.labeled.correction.list=ang.18O.td.pos.resid.list, filename="R1_Ang_18O_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R1.ang.18O.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")

R2.ang.18O.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R2.ang.18O.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
find.labeled.correction.plot(find.labeled.correction.list=R2.ang.18O.td.pos.resid.list, filename="R2_Ang_18O_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R2.ang.18O.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")

R3.ang.18O.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R3.ang.18O.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
find.labeled.correction.plot(find.labeled.correction.list=R3.ang.18O.td.pos.resid.list, filename="R3_Ang_18O_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R2.ang.18O.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")

R1.hop.18O.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=hop.18O.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
find.labeled.correction.plot(find.labeled.correction.list=hop.18O.td.pos.resid.list, filename="R1_Ang_18O_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R1.hop.18O.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")

R2.hop.18O.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R2.hop.18O.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
find.labeled.correction.plot(find.labeled.correction.list=R2.hop.18O.td.pos.resid.list, filename="R2_Ang_18O_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R2.hop.18O.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")

R3.hop.18O.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R3.hop.18O.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
find.labeled.correction.plot(find.labeled.correction.list=R3.hop.18O.td.pos.resid.list, filename="R3_Ang_18O_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R3.hop.18O.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")

R1.sed.18O.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=sed.18O.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
find.labeled.correction.plot(find.labeled.correction.list=sed.18O.td.pos.resid.list, filename="R1_Ang_18O_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R1.sed.18O.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")

R2.sed.18O.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R2.sed.18O.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
find.labeled.correction.plot(find.labeled.correction.list=R2.sed.18O.td.pos.resid.list, filename="R2_Ang_18O_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R2.sed.18O.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")

R3.sed.18O.td.pos.resid.best.iteration.list <- select.best.iteration(putative.nongrower.metrics=R3.sed.18O.td.pos.resid.list$putative.nongrower.metrics, quantile.threshold=0.10)
find.labeled.correction.plot(find.labeled.correction.list=R3.sed.18O.td.pos.resid.list, filename="R3_Ang_18O_td_pos_resid_plots", path="qSIP_output/Figures/", highlight=R3.sed.18O.td.pos.resid.best.iteration.list$best.iteration, highlight.col="magenta")

#Correct the summarized and raw data for tube-level shifts in WAD identified above for ALL unlabeled and labeled replicates:
#Apply the tube-level shift to the already-summarized taxon-level labeled replicate WADs to 'correct' them and create a new summary list (analagous to 'unlab.WAD.corr.list') that includes those corrected values:

#First, define character vectors listing the set of labeled replicates and the names of the shifts corresponding to those replicates:
labeled.treatments.sed <- c("Sedgwick_18")
all.possible.labeled.replicates.sed <- as.vector(outer(c("R1", "R2", "R3"), labeled.treatments.sed, paste, sep="."))

labeled.treatments.2.sed <- c("sed_18O")
all.possible.labeled.replicates.sed <- as.vector(outer(c("R1", "R2", "R3"), labeled.treatments.2.sed, paste, sep="."))
actual.labeled.replicates.dots.sed <- gsub(pattern="\\_", replacement="\\.", x=all.possible.labeled.replicates.sed, perl=TRUE)
actual.labeled.replicates.dots.corrections.sed <- paste(actual.labeled.replicates.dots.sed, "td.pos.resid.best.iteration.list$norm.correction", sep=".")

#Evaluate the names of the corrections to get a vector of the corrections themselves:
actual.labeled.replicates.dots.corrections.values.sed <- unlist(lapply(as.list(actual.labeled.replicates.dots.corrections.sed), function(x) eval(parse(text=x))))

WAD.by.taxon$reps.by.trt.sed <- WAD.by.taxon$reps.by.trt %>%
  filter(Soil_Isotope == "Sedgwick_18" | Soil_Isotope == "Sedgwick_16")

#Now, apply the tube-level shifts for all of those replicates: lab resp to add, reps.by.trt
WAD.corr.list.sed <- add.lab.WAD.corr.summary(summary.list=unlab.WAD.corr.list.sed, reps.by.trt=WAD.by.taxon$reps.by.trt.sed, lab.reps.to.add=all.possible.labeled.replicates.sed, lab.shifts=actual.labeled.replicates.dots.corrections.values.sed, CI=0.90)

#Do the same for Hopland
labeled.treatments.hop <- c("Hopland_18")
all.possible.labeled.replicates.hop <- as.vector(outer(c("R1", "R2", "R3"), labeled.treatments.hop, paste, sep="."))

labeled.treatments.2.hop <- c("hop_18O")
all.possible.labeled.replicates.hop.2 <- as.vector(outer(c("R1", "R2", "R3"), labeled.treatments.2.hop, paste, sep="."))
actual.labeled.replicates.dots.hop <- gsub(pattern="\\_", replacement="\\.", x=all.possible.labeled.replicates.hop.2, perl=TRUE)
actual.labeled.replicates.dots.corrections.hop <- paste(actual.labeled.replicates.dots.hop, "td.pos.resid.best.iteration.list$norm.correction", sep=".")

#Evaluate the names of the corrections to get a vector of the corrections themselves:
actual.labeled.replicates.dots.corrections.values.hop <- unlist(lapply(as.list(actual.labeled.replicates.dots.corrections.hop), function(x) eval(parse(text=x))))

WAD.by.taxon$reps.by.trt.hop <- WAD.by.taxon$reps.by.trt %>%
  filter(Soil_Isotope == "Hopland_18" | Soil_Isotope == "Hopland_16")

#Now, apply the tube-level shifts for all of those replicates: lab resp to add, reps.by.trt
WAD.corr.list.hop <- add.lab.WAD.corr.summary(summary.list=unlab.WAD.corr.list.hop, reps.by.trt=WAD.by.taxon$reps.by.trt.hop, lab.reps.to.add=all.possible.labeled.replicates.hop, lab.shifts=actual.labeled.replicates.dots.corrections.values.hop, CI=0.90)

#Do the same for Angelo
labeled.treatments.ang <- c("Angelo_18")
all.possible.labeled.replicates.ang <- as.vector(outer(c("R1", "R2", "R3"), labeled.treatments.ang, paste, sep="."))

labeled.treatments.2.ang <- c("ang_18O")
all.possible.labeled.replicates.ang.2 <- as.vector(outer(c("R1", "R2", "R3"), labeled.treatments.2.ang, paste, sep="."))
actual.labeled.replicates.dots.ang <- gsub(pattern="\\_", replacement="\\.", x=all.possible.labeled.replicates.ang.2, perl=TRUE)
actual.labeled.replicates.dots.corrections.ang <- paste(actual.labeled.replicates.dots.ang, "td.pos.resid.best.iteration.list$norm.correction", sep=".")

#Evaluate the names of the corrections to get a vector of the corrections themselves:
actual.labeled.replicates.dots.corrections.values.ang <- unlist(lapply(as.list(actual.labeled.replicates.dots.corrections.ang), function(x) eval(parse(text=x))))

WAD.by.taxon$reps.by.trt.ang <- WAD.by.taxon$reps.by.trt %>%
  filter(Soil_Isotope == "Angelo_18" | Soil_Isotope == "Angelo_16")

#Now, apply the tube-level shifts for all of those replicates: lab resp to add, reps.by.trt
WAD.corr.list.ang <- add.lab.WAD.corr.summary(summary.list=unlab.WAD.corr.list.ang, reps.by.trt=WAD.by.taxon$reps.by.trt.ang, lab.reps.to.add=all.possible.labeled.replicates.ang, lab.shifts=actual.labeled.replicates.dots.corrections.values.ang, CI=0.90)

#Look at the results:
names(WAD.corr.list.ang)                           #names of the two data frames in the output list
WAD.corr.list.ang$WAD.norm.fit.parms               #looking at the first object -- a data frame
WAD.corr.list.ang$corr.names                       #looking at the second object in the list -- a vector of the names of the corrected unlabeled replicates
head(WAD.corr.list.ang$WAD.table.corr)             #looking at the head of the third object -- a data frame

#get raw data
data<- read.table("table1.csv", header=TRUE, sep=",") #sample data
data<-data %>% select(-X) #remove random column generated by loading file

#reorder columns
data<- data %>%
  select("Sample_ID","Isotope_Treatment", "Tube", "Fraction", "Soil_Isotope", "Density_g_ml", "DNA_ng_ul", "Total_16S_copies", everything())

data.sed<- data %>%
  filter(Soil_Isotope == "Sedgwick_16" | Soil_Isotope == "Sedgwick_18")

#Apply the tube-level shifts to 'correct' the unlabeled WADs in 'data' (corrects all unlabeled replicates at once):
data.corr.sed <- apply.unlabeled.correction(raw.data=data.sed, correction.table=unlab.WAD.corr.list$WAD.norm.fit.parms, reps.by.trt=WAD.by.taxon$reps.by.trt.sed, vars=c("Density_g_ml", "Tube", "Soil_Isotope"))

#Apply the tube-level shift to 'correct' the labeled WADs in 'data.corr' (corrects one labeled replicate at a time):
for (i in 1:length(all.possible.labeled.replicates)){
  data.corr.sed <- apply.labeled.correction(raw.data=data.sed, raw.data.corr=data.corr.sed, lab.replicate=all.possible.labeled.replicates[i], correction.value=eval(parse(text=actual.labeled.replicates.dots.corrections[i])), reps.by.trt=WAD.by.taxon$reps.by.trt.sed, vars=c("Density_g_ml", "Tube", "Soil_Isotope"))
}

#Re-calculate number of copies per uL, based on relative abundance and total number of copies per uL:
ncopies.corr.sed <- data.corr.sed$Total_16S_copies*data.corr.sed[,9:(ncol(data.corr.sed)-1)]

ncopies.corr.sed <- cbind(data.corr.sed[,1:8], ncopies.corr.sed)  # add first 15 columns of data.corr to ncopies.corr

dim(ncopies.corr.sed)


#Melt data.corr into long format by tube, sample, tmt, rep, fraction, DNA conc, and density;
#Do this for copies.ul and for relative abundance, which is just our data.corr file. Merge these to into 1 masterfile: data.corr.melted
ncopies.corr.sed.melted <- melt(ncopies.corr.sed, id=c("Sample_ID","Isotope_Treatment", "Tube", "Fraction", "Soil_Isotope", "Density_g_ml", "DNA_ng_ul", "Total_16S_copies"), variable.name="taxon", value.name="copies.ul")


rel.abundance.corr.sed.melted <- melt(data.corr.sed, id=c("Sample_ID","Isotope_Treatment", "Tube", "Fraction", "Soil_Isotope", "Density_g_ml", "DNA_ng_ul", "Total_16S_copies"), variable.name="taxon", value.name="rel.abundance")
data.corr.melted.sed <- merge(ncopies.corr.sed.melted, rel.abundance.corr.sed.melted)
head(data.corr.melted.sed)

#Merge taxa data and reorder data frame by taxon and SampleID and fraction:
#data.corr.melted <- merge(data.corr.melted, taxa.id)
data.corr.melted.sed <- data.corr.melted.sed[order(data.corr.melted.sed$taxon, data.corr.melted.sed$Sample_ID, data.corr.melted.sed$Fraction),]

row.names(data.corr.melted.sed) <- 1:dim(data.corr.melted.sed)[1]   #rename observations to be sequential
head(data.corr.melted.sed)

#Calculate corrected WADs for each tube and taxon:
WAD.by.taxon.corr.sed <- WAD.by.taxon.func(X=data.corr.melted.sed, vars=c("taxon", "Density_g_ml", "copies.ul", "Tube", "Isotope_Treatment"))

#Look at the results:
names(WAD.by.taxon.corr.sed)          #names of the two data frames in the output list
head(WAD.by.taxon.corr.sed$obs.wads)  #looking at the head of the first data frame
WAD.by.taxon.corr.sed$reps.by.trt     #looking at the head of the second data frame


write.csv(data.corr.melted.sed, "data.corr.melted.sed.csv")

Tcompare<- read.table("trt.comparisons.csv", header=TRUE, sep=",", colClasses=c("numeric","factor","character","character","character","numeric","character"))

#AP First remove taxa in tubes where the taxon appears in less than 4 fractions per tube
#AP first count the number of rows in a tube for each taxon where t.copies.ul does not equal 0 (i.e. how many fractions in each tube #does each taxon appear?)
data.melted.hopland$Tube<- as.factor(data.melted.hopland$Tube)
is.numeric(data.melted.hopland$Density_g_ml)

data.melted.heavyfilt <- data.corr.melted.sed %>%
  group_by(taxon, Tube) %>%
  summarise_at(vars(copies.ul), ~sum(. != 0))

head(data.melted.heavyfilt)

#AP rename column in test to "nonzero.count.t.copies.ul" - this signifies how many fractions contain nonzero copies of the taxon
names(data.melted.heavyfilt)[3] <- "nonzero.count.t.copies.ul."
head(data.melted.heavyfilt)

#AP join to data.melted.2 the new "data.melted.2.heavyfilt" data frame which adds the new column of nonzero count to the original data.melted
join.test <- left_join(data.corr.melted.sed, data.melted.heavyfilt)
View(join.test)

#AP remove rows from join.test where nonzero.count.t.copies.ul is less than 3 (i.e. remove rows of the data.melted for a taxon for each tube where the taxon is present in less than 3 fractions)
data.melted.heavyfilt <- subset(join.test, join.test$nonzero.count.t.copies.ul > 1)

#this dataframe contains taxa that appear in at least 4 replicates in each tube
#now use Ben's code to filter so that taxa appear in at least 2 out of 3 biological replicates
#create a single data.melted file for each site that contains Hopland, Angelo, and Sedwick separately

#Sedgwick
data.melted.sedgwick <- filter.taxa(DATA=data.melted.heavyfilt, trt.code.1= Tcompare$trt.code.1[3], trt.code.2= NULL, trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "Soil_Isotope"), min.reps=2)
data.melted.sedgwick <- filter.taxa(DATA=data.melted.sedgwick, trt.code.1=NULL, trt.code.2=Tcompare$trt.code.2[3], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "Soil_Isotope"), min.reps=2)

#Files are ready for all.taxa.calcs
set.seed(100)

system.time(all.comparisons.sedgwick <- all.taxa.calcs(X.all=data.melted.sedgwick, comparisons=Tcompare[3,], vars=c("taxon", "Density_g_ml", "copies.ul", "Tube", "Soil_Isotope", "DNA_ng_uL"), growth.model="exponential", prop.O.from.water=0.6, v.frac=50, copies.cell=6, pgC.cell=0.1, CI=0.95, draws=1000, tailed.test=1))
write.csv(all.comparisons.sedgwick, file = "qSIP_output/all.comparisons.sedgwick.itag.csv")

##################ANGELO
data.ang<- data %>%
  filter(Soil_Isotope == "Angelo_16" | Soil_Isotope == "Angelo_18")

#Apply the tube-level shifts to 'correct' the unlabeled WADs in 'data' (corrects all unlabeled replicates at once):
data.corr.ang <- apply.unlabeled.correction(raw.data=data.ang, correction.table=unlab.WAD.corr.list.ang$WAD.norm.fit.parms, reps.by.trt=WAD.by.taxon$reps.by.trt.ang, vars=c("Density_g_ml", "Tube", "Soil_Isotope"))


#Apply the tube-level shift to 'correct' the labeled WADs in 'data.corr' (corrects one labeled replicate at a time):
for (i in 1:length(all.possible.labeled.replicates.ang)){
  data.corr.ang <- apply.labeled.correction(raw.data=data.ang, raw.data.corr=data.corr.ang, lab.replicate=all.possible.labeled.replicates.ang[i], correction.value=eval(parse(text=actual.labeled.replicates.dots.corrections.ang[i])), reps.by.trt=WAD.by.taxon$reps.by.trt.ang, vars=c("Density_g_ml", "Tube", "Soil_Isotope"))
}

#Re-calculate number of copies per uL, baang on relative abundance and total number of copies per uL:
ncopies.corr.ang <- data.corr.ang$Total_16S_copies*data.corr.ang[,9:(ncol(data.corr.ang)-1)]

ncopies.corr.ang <- cbind(data.corr.ang[,1:8], ncopies.corr.ang)  # add first 15 columns of data.corr to ncopies.corr

dim(ncopies.corr.ang)


#Melt data.corr into long format by tube, sample, tmt, rep, fraction, DNA conc, and density;
#Do this for copies.ul and for relative abundance, which is just our data.corr file. Merge these to into 1 masterfile: data.corr.melted
ncopies.corr.ang.melted <- melt(ncopies.corr.ang, id=c("Sample_ID","Isotope_Treatment", "Tube", "Fraction", "Soil_Isotope", "Density_g_ml", "DNA_ng_ul", "Total_16S_copies"), variable.name="taxon", value.name="copies.ul")


rel.abundance.corr.ang.melted <- melt(data.corr.ang, id=c("Sample_ID","Isotope_Treatment", "Tube", "Fraction", "Soil_Isotope", "Density_g_ml", "DNA_ng_ul", "Total_16S_copies"), variable.name="taxon", value.name="rel.abundance")
data.corr.melted.ang <- merge(ncopies.corr.ang.melted, rel.abundance.corr.ang.melted)
head(data.corr.melted.ang)

#Merge taxa data and reorder data frame by taxon and SampleID and fraction:
#data.corr.melted <- merge(data.corr.melted, taxa.id)
data.corr.melted.ang <- data.corr.melted.ang[order(data.corr.melted.ang$taxon, data.corr.melted.ang$Sample_ID, data.corr.melted.ang$Fraction),]

row.names(data.corr.melted.ang) <- 1:dim(data.corr.melted.ang)[1]   #rename observations to be sequential
head(data.corr.melted.ang)

write.csv(data.corr.melted.ang, "data.corr.melted.ang.csv")

##################HOPLAND
data.hop<- data %>%
  filter(Soil_Isotope == "Hopland_16" | Soil_Isotope == "Hopland_18")

#Apply the tube-level shifts to 'correct' the unlabeled WADs in 'data' (corrects all unlabeled replicates at once):
data.corr.hop <- apply.unlabeled.correction(raw.data=data.hop, correction.table=unlab.WAD.corr.list.hop$WAD.norm.fit.parms, reps.by.trt=WAD.by.taxon$reps.by.trt.hop, vars=c("Density_g_ml", "Tube", "Soil_Isotope"))


#Apply the tube-level shift to 'correct' the labeled WADs in 'data.corr' (corrects one labeled replicate at a time):
for (i in 1:length(all.possible.labeled.replicates.hop)){
  data.corr.hop <- apply.labeled.correction(raw.data=data.hop, raw.data.corr=data.corr.hop, lab.replicate=all.possible.labeled.replicates.hop[i], correction.value=eval(parse(text=actual.labeled.replicates.dots.corrections.hop[i])), reps.by.trt=WAD.by.taxon$reps.by.trt.hop, vars=c("Density_g_ml", "Tube", "Soil_Isotope"))
}

#Re-calculate number of copies per uL, bahop on relative abundance and total number of copies per uL:
ncopies.corr.hop <- data.corr.hop$Total_16S_copies*data.corr.hop[,9:(ncol(data.corr.hop)-1)]

ncopies.corr.hop <- cbind(data.corr.hop[,1:8], ncopies.corr.hop)  # add first 15 columns of data.corr to ncopies.corr

dim(ncopies.corr.hop)


#Melt data.corr into long format by tube, sample, tmt, rep, fraction, DNA conc, and density;
#Do this for copies.ul and for relative abundance, which is just our data.corr file. Merge these to into 1 masterfile: data.corr.melted
ncopies.corr.hop.melted <- melt(ncopies.corr.hop, id=c("Sample_ID","Isotope_Treatment", "Tube", "Fraction", "Soil_Isotope", "Density_g_ml", "DNA_ng_ul", "Total_16S_copies"), variable.name="taxon", value.name="copies.ul")


rel.abundance.corr.hop.melted <- melt(data.corr.hop, id=c("Sample_ID","Isotope_Treatment", "Tube", "Fraction", "Soil_Isotope", "Density_g_ml", "DNA_ng_ul", "Total_16S_copies"), variable.name="taxon", value.name="rel.abundance")
data.corr.melted.hop <- merge(ncopies.corr.hop.melted, rel.abundance.corr.hop.melted)
head(data.corr.melted.hop)

#Merge taxa data and reorder data frame by taxon and SampleID and fraction:
#data.corr.melted <- merge(data.corr.melted, taxa.id)
data.corr.melted.hop <- data.corr.melted.hop[order(data.corr.melted.hop$taxon, data.corr.melted.hop$Sample_ID, data.corr.melted.hop$Fraction),]

row.names(data.corr.melted.hop) <- 1:dim(data.corr.melted.hop)[1]   #rename observations to be sequential
head(data.corr.melted.hop)

write.csv(data.corr.melted.hop, "data.corr.melted.hop.csv")

#AP First remove taxa in tubes where the taxon appears in less than 4 fractions per tube
#AP first count the number of rows in a tube for each taxon where t.copies.ul does not equal 0 (i.e. how many fractions in each tube #does each taxon appear?)
data.melted.corr.hop <- read.table("data.corr.melted.hop.csv", header = TRUE, sep = ",")
library(dplyr)

data.melted.heavyfilt <- data.melted.corr.hop %>%
  group_by(taxon, Tube) %>%
  summarise_at(vars(copies.ul), ~sum(. != 0))

head(data.melted.heavyfilt)

#AP rename column in test to "nonzero.count.t.copies.ul" - this signifies how many fractions contain nonzero copies of the taxon
names(data.melted.heavyfilt)[3] <- "nonzero.count.t.copies.ul."
head(data.melted.heavyfilt)

#AP join to data.melted.2 the new "data.melted.2.heavyfilt" data frame which adds the new column of nonzero count to the original data.melted
join.test <- left_join(data.melted.corr.hop, data.melted.heavyfilt)
View(join.test)

#AP remove rows from join.test where nonzero.count.t.copies.ul is less than 3 (i.e. remove rows of the data.melted for a taxon for each tube where the taxon is present in less than 3 fractions)
data.melted.heavyfilt <- subset(join.test, join.test$nonzero.count.t.copies.ul > 1)

#this dataframe contains taxa that appear in at least 4 replicates in each tube
#now use Ben's code to filter so that taxa appear in at least 2 out of 3 biological replicates
#create a single data.melted file for each site that contains Hopland, Angelo, and Sedwick separately

data.melted.hopland <- filter.taxa(DATA=data.melted.heavyfilt, trt.code.1= Tcompare$trt.code.1[1], trt.code.2= NULL, trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "Soil_Isotope"), min.reps=2)
data.melted.hopland <- filter.taxa(DATA=data.melted.hopland, trt.code.1=NULL, trt.code.2=Tcompare$trt.code.2[1], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "Soil_Isotope"), min.reps=2)
data.melted.hopland<- data.melted.hopland %>%
  select(-X)

#Files are ready for all.taxa.calcs
set.seed(100)
head(data.melted.hopland)

system.time(all.comparisons.hopland <- all.taxa.calcs(X.all=data.melted.hopland, comparisons=Tcompare[1,], vars=c("taxon", "Density_g_ml", "copies.ul", "Tube", "Soil_Isotope", "DNA_ng_ul"), growth.model="exponential", prop.O.from.water=0.6, v.frac=50, copies.cell=6, pgC.cell=0.1, CI=0.95, draws=1000, tailed.test=1))
write.csv(all.comparisons.sedgwick, file = "qSIP_output/all.comparisons.sedgwick.itag.csv")
