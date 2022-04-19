setwd("/Users/meganfoley/Documents/winterSIP/manuscript/ember.correction")

data.corr.melted.hop <- read.table("data.corr.melted.hop.txt", header = TRUE, sep = "\t")
Tcompare <- read.table("trt.comparisons.txt", header=TRUE, sep="\t", colClasses=c("numeric","factor","character","character","character","numeric","character"))

library(dplyr)
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

data.corr.melted.hop<- data.corr.melted.hop%>%
  select(-X)

data.corr.melted.hop$Isotope_Treatment<- paste(data.corr.melted.hop$Isotope_Treatment,"O",sep = "")
data.corr.melted.hop$Fraction<-as.factor(data.corr.melted.hop$Fraction)
data.corr.melted.hop$Tube<-as.factor(data.corr.melted.hop$Tube)
data.corr.melted.hop$Isotope_Treatment<-as.factor(data.corr.melted.hop$Isotope_Treatment)
data.corr.melted.hop$Soil_Isotope<-as.character(data.corr.melted.hop$Soil)
data.corr.melted.hop$Sample_ID<-as.character(data.corr.melted.hop$Sample_ID)
data.corr.melted.hop$taxon<-as.factor(data.corr.melted.hop$taxon)

#AP First remove taxa in tubes where the taxon appears in less than 4 fractions per tube
#AP first count the number of rows in a tube for each taxon where t.copies.ul does not equal 0 (i.e. how many fractions in each tube #does each taxon appear?)
data.melted.heavyfilt <- data.corr.melted.hop %>%
  group_by(taxon, Tube) %>%
  summarise_at(vars(copies.ul), ~sum(. != 0))

head(data.melted.heavyfilt)

#AP rename column in test to "nonzero.count.t.copies.ul" - this signifies how many fractions contain nonzero copies of the taxon
names(data.melted.heavyfilt)[3] <- "nonzero.count.t.copies.ul."
head(data.melted.heavyfilt)

#AP join to data.melted.2 the new "data.melted.2.heavyfilt" data frame which adds the new column of nonzero count to the original data.melted
join.test <- left_join(data.corr.melted.hop, data.melted.heavyfilt)

#AP remove rows from join.test where nonzero.count.t.copies.ul is less than 3 (i.e. remove rows of the data.melted for a taxon for each tube where the taxon is present in less than 3 fractions)
data.melted.heavyfilt <- subset(join.test, join.test$nonzero.count.t.copies.ul > 1)

#this dataframe contains taxa that appear in at least 4 replicates in each tube
#now use Ben's code to filter so that taxa appear in at least 2 out of 3 biological replicates
#create a single data.melted file for each site that contains Hopland, Angelo, and Sedwick separately

#Sedgwick
data.melted.hopland <- filter.taxa(DATA=data.melted.heavyfilt, trt.code.1= Tcompare$trt.code.1[1], trt.code.2= NULL, trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "Soil_Isotope"), min.reps=2)
data.melted.hopland <- filter.taxa(DATA=data.melted.hopland, trt.code.1=NULL, trt.code.2=Tcompare$trt.code.2[1], trt.refs=NULL, vars=c("taxon", "copies.ul", "Tube", "Soil_Isotope"), min.reps=2)

#Files are ready for all.taxa.calcs
set.seed(100)

system.time(all.comparisons.hopland <- all.taxa.calcs(X.all=data.melted.hopland, comparisons=Tcompare[1,], vars=c("taxon", "Density_g_ml", "copies.ul", "Tube", "Soil_Isotope", "DNA_ng_ul"), growth.model="exponential", prop.O.from.water=0.6, v.frac=50, copies.cell=6, pgC.cell=0.1, CI=0.95, draws=1000, tailed.test=1))
write.csv(all.comparisons.hopland, file = "qSIP_output/all.comparisons.hopland.csv")
