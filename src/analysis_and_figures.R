---
title: "It's a wormy world: Meta-analysis reveals several decades of change in the global abundance of
the parasitic nematodes Anisakis spp. and Pseudoterranova spp. in marine fishes and invertebrates"
author: "Evan Fiorenza and Chelsea Wood"
date: "updated 21 January 2020"
---
  
rm(list=ls())#Clear Workspace

#Load required Packages---------------

library(metafor)
library(ggplot2)
library(rfishbase)
library(optimx)
library(maptools)
library(maps)
library(lmerTest)

#Load Required Data files----------------
Dat = data.frame(read.table('Data/anisakid_data.csv', header = TRUE, sep=','))#Load the full anisakid dataset, including the publications that we excluded and didn't take data from 
FAO = data.frame(read.table('Data/anisakid_fao.csv', header = TRUE, sep=',')) #Load a dataset that contains the FAO area data obtained by matching the data points to FAO polygons in ArcGIS
Hosts = data.frame(read.table('Data/host_data.csv', header = TRUE, sep=','))#Load the host list to use to get data from Fishbase

#Data Manipulation--------------------
Dat2=Dat[which(Dat$Include_abstract==1),]#Trim the data to only include the publications that made it past the title and abstract screening stages
AllDat=Dat2[which(Dat2$Use_Data==''),]#Trim the data to only include publications that we included in the final data set

# Calculate the standardized abundance by either using the reported mean abundance or prevalence * intensity
AllDat$Standard_Abundance=(AllDat$Prevalence/100)*AllDat$Intensity_Mean
AllDat$Standard_Abundance[which(AllDat$Abundance_mean_per_fish>=0)]=AllDat$Abundance_mean_per_fish[which(AllDat$Abundance_mean_per_fish>=0)]

#Standardize year by taking the reported year for those that are in a single year or by taking the median if a range of years is reported
AllDat$Standard_Year=AllDat$Year_collection
AllDat$Standard_Year[which(AllDat$Year_collection_min>0)]=(AllDat$Year_collection_min[which(AllDat$Year_collection_min>0)]+AllDat$Year_collection_max[which(AllDat$Year_collection_min>0)])/2

#Standardize to the standard deviation of the abundance when an error measure is reported (Std Dev, Std Error, 95% Confidence interval)
AllDat$Standard_AbundanceStdDev=NA
AllDat$Standard_AbundanceStdDev[which(AllDat$Abundance_SD>=0)]=AllDat$Abundance_SD[which(AllDat$Abundance_SD>=0)]
AllDat$Standard_AbundanceStdDev[which(AllDat$Abundance_SE>=0)]=AllDat$Abundance_SE[which(AllDat$Abundance_SE>=0)]*sqrt(AllDat$Total_Fish_Examined[which(AllDat$Abundance_SE>=0)]) #Back calculate from Std Error to Std Deviation
AllDat$Standard_AbundanceStdDev[which(AllDat$Abundance_ConfidenceLevel95_Max>=0)]=(AllDat$Abundance_ConfidenceLevel95_Max[which(AllDat$Abundance_ConfidenceLevel95_Max>=0)]-AllDat$Abundance_mean_per_fish[which(AllDat$Abundance_ConfidenceLevel95_Max>=0)])/1.96*sqrt(AllDat$Total_Fish_Examined[which(AllDat$Abundance_ConfidenceLevel95_Max>=0)])#Back Calculate from the 95% CI to Std deviation
AllDat$Standard_AbundanceStdDev[which(AllDat$Abundance_Variance>=0)]=sqrt(AllDat$Abundance_Variance[which(AllDat$Abundance_Variance>=0)])

#Standardize to the standard deviation of the intensity when an error measure is reported (Std Dev, Std Error, 95% Confidence interval)
AllDat$Standard_IntensityStdDev=NA
AllDat$Standard_IntensityStdDev[which(AllDat$Intensity_SD>=0)]=AllDat$Intensity_SD[which(AllDat$Intensity_SD>=0)]
AllDat$Standard_IntensityStdDev[which(AllDat$Intensity_SE>=0)]=AllDat$Intensity_SE[which(AllDat$Intensity_SE>=0)]*sqrt(AllDat$Total_Fish_Examined[which(AllDat$Intensity_SE>=0)])
AllDat$Standard_IntensityStdDev[which(AllDat$Intensity_ConfidenceLevel95_Max>=0)]=(AllDat$Intensity_ConfidenceLevel95_Max[which(AllDat$Intensity_ConfidenceLevel95_Max>=0)]-AllDat$Intensity_Mean[which(AllDat$Intensity_ConfidenceLevel95_Max>=0)])*sqrt(AllDat$Total_Fish_Examined[which(AllDat$Intensity_ConfidenceLevel95_Max>=0)])/1.96

#Calculate an estimated Standard deviation when a range of abundance or intensities is reported
AllDat$Calc_StdDev=NA
abmax=which(AllDat$Abundance_Range_Max>=0) #find the subset of the data where the abundance range is reported
imax=which(AllDat$Intensity_Range_High>=0) #find the subset of the data where the intensity range is reported

#Create a function to optimize to find the standard deviation
#This function assumes that all the data follows a negative binomial distribution and that the max reported is the 95th quantile
STDNegBinom=function(size){
  abs(qnbinom(.95,size=size,mu=AllDat$Standard_Abundance[x])-AllDat$Abundance_Range_Max[x])
}

#Run the optimization function for abundance ranges
for(x in abmax){
  z=optim(1,STDNegBinom,lower=0,upper=10,method='Brent')#Optimize  
  AllDat$Calc_StdDev[x]=sqrt(AllDat$Standard_Abundance[x]+(((AllDat$Standard_Abundance[x])^2)/z$par)) #Calculates the standard deviation based on the formula for mean and size for a negative binomial distribution 
  
}

#Run the optimization function for intensity ranges
for(x in imax){
  z=optim(10,STDNegBinom,lower=.00000000001,upper=100,method='Brent')  
  AllDat$Calc_StdDev[x]=sqrt(AllDat$Standard_Abundance[x]+(((AllDat$Standard_Abundance[x])^2)/z$par))
  
}

#Standardize reported lengths to a single measure, either the reported mean or the median of a range for each type of length reported. When only 'length; is reported, assumed to be standard length
#For species where we can't convert to standard length, we use the reported fork or total length for correction
AllDat$Standard_Length=AllDat$Length_cm
AllDat$Standard_Length[which(AllDat$Length_min_cm>0)]=(AllDat$Length_min_cm[which(AllDat$Length_min_cm>0)]+AllDat$Length_max_cm[which(AllDat$Length_min_cm>0)])/2
AllDat$Standard_Length[which(AllDat$Length_median_cm>0)]=AllDat$Length_median_cm[which(AllDat$Length_median_cm>0)]
AllDat$Standard_TLength=AllDat$Length_TL_cm
AllDat$Standard_TLength[which(AllDat$Length_TL_cm_min>0)]=(AllDat$Length_TL_cm_min[which(AllDat$Length_TL_cm_min>0)]+AllDat$Length_TL_cm_max[which(AllDat$Length_TL_cm_min>0)])/2
AllDat$Standard_SLength=AllDat$Length_SL_cm
AllDat$Standard_SLength[which(AllDat$Length_SL_min_cm>0)]=(AllDat$Length_SL_min_cm[which(AllDat$Length_SL_min_cm>0)]+AllDat$Length_SL_max_cm[which(AllDat$Length_SL_min_cm>0)])/2
AllDat$Standard_FLength=AllDat$Length_FL_cm_mean
AllDat$Standard_FLength[which(AllDat$Length_FL_cm_min>0)]=(AllDat$Length_FL_cm_min[which(AllDat$Length_FL_cm_min>0)]+AllDat$Length_FL_cm_max[which(AllDat$Length_FL_cm_min>0)])/2
AllDat$Standard_DMLength=AllDat$Length_Dorsal_Mantle_cm
AllDat$Standard_DMLength[which(AllDat$Length_Dorsal_Mantle_Range_Min>0)]=(AllDat$Length_Dorsal_Mantle_Range_Min_cm[which(AllDat$Length_Dorsal_Mantle_Range_Min_cm>0)]+AllDat$Length_Dorsal_Mantle_Range_Max_cm[which(AllDat$Length_Dorsal_Mantle_Range_Min_cm>0)])/2
AllDat$Standard_PALength=AllDat$Length_PAL_cm


#Extract data from fishbase

#Extract and take the mean of the length-length conversions on fishbase
for (x in 1:length(Hosts$Scientific)){
  zzz=length_length(as.character(Hosts$Scientific[x]))
  Hosts$FLSLa[x]=mean(zzz$a[which(zzz$Length1=='SL' & zzz$Length2=='FL')])
  Hosts$FLSLb[x]=mean(zzz$b[which(zzz$Length1=='SL' & zzz$Length2=='FL')])
  Hosts$TLSLa[x]=mean(zzz$a[which(zzz$Length1=='SL' & zzz$Length2=='TL')])
  Hosts$TLSLb[x]=mean(zzz$b[which(zzz$Length1=='SL' & zzz$Length2=='TL')])
  Hosts$FLTLa[x]=mean(zzz$a[which(zzz$Length1=='TL' & zzz$Length2=='FL')])
  Hosts$FLTLb[x]=mean(zzz$b[which(zzz$Length1=='TL' & zzz$Length2=='FL')])
  Hosts$TLFLa[x]=mean(zzz$a[which(zzz$Length1=='FL' & zzz$Length2=='TL')])
  Hosts$TLFLb[x]=mean(zzz$b[which(zzz$Length1=='FL' & zzz$Length2=='TL')])
}


#Export a .csv file with all the host metadata from fishbase
#write.table(Hosts, "host_metadata.csv", col.names = TRUE, row.names = FALSE, sep = ',')

#March extracted data from the host metadata and fishbase to the anisakid dataset
AllDat$TrueHost=Hosts$Scientific[match(AllDat$Host,Hosts$DataSpecies)]#Correct spelling and scientific names to the accepted ones

#Append Length conversion coefficients to Anisakid data
AllDat$FLSLa=Hosts$FLSLa[match(AllDat$Host,Hosts$DataSpecies)]
AllDat$FLSLb=Hosts$FLSLb[match(AllDat$Host,Hosts$DataSpecies)]
AllDat$TLSLa=Hosts$TLSLa[match(AllDat$Host,Hosts$DataSpecies)]
AllDat$TLSLb=Hosts$TLSLb[match(AllDat$Host,Hosts$DataSpecies)]
AllDat$FLTLa=Hosts$FLTLa[match(AllDat$Host,Hosts$DataSpecies)]
AllDat$FLTLb=Hosts$FLTLb[match(AllDat$Host,Hosts$DataSpecies)]
AllDat$TLFLa=Hosts$TLFLa[match(AllDat$Host,Hosts$DataSpecies)]
AllDat$TLFLb=Hosts$TLFLb[match(AllDat$Host,Hosts$DataSpecies)]

#Correct all length measurements to Standard/Preanal/Dorsal-mantle length (species dependent, though all are essentially the length of the body)
AllDat$Standard_LengthCalc=AllDat$Standard_Length
AllDat$Standard_LengthCalcFL[is.na(AllDat$TLSLb)]=AllDat$Standard_TLength[is.na(AllDat$TLSLb)]*AllDat$TLFLb[is.na(AllDat$TLSLb)]+AllDat$TLFLa[is.na(AllDat$TLSLb)]
AllDat$Standard_LengthCalcTL[is.na(AllDat$FLSLb)]=AllDat$Standard_FLength[is.na(AllDat$FLSLb)]*AllDat$FLTLb[is.na(AllDat$FLSLb)]+AllDat$FLTLa[is.na(AllDat$FLSLb)]
AllDat$Standard_LengthCalc[which(AllDat$Standard_SLength>0)]=AllDat$Standard_SLength[which(AllDat$Standard_SLength>0)]
AllDat$Standard_LengthCalc[which(AllDat$Standard_TLength>0)]=AllDat$TLSLa[which(AllDat$Standard_TLength>0)]+AllDat$TLSLb[which(AllDat$Standard_TLength>0)]*AllDat$Standard_TLength[which(AllDat$Standard_TLength>0)]
AllDat$Standard_LengthCalc[which(AllDat$Standard_FLength>0)]=AllDat$FLSLa[which(AllDat$Standard_FLength>0)]+AllDat$FLSLb[which(AllDat$Standard_FLength>0)]*AllDat$Standard_FLength[which(AllDat$Standard_FLength>0)]
AllDat$Standard_LengthCalc[which(AllDat$Standard_PALength>0)]=AllDat$Length_PAL_cm[which(AllDat$Standard_PALength>0)]
AllDat$Standard_LengthCalc[which(AllDat$Standard_LengthCalcTL>0)]=AllDat$TLSLa[which(AllDat$Standard_LengthCalcTL>0)]+AllDat$TLSLb[which(AllDat$Standard_LengthCalcTL>0)]*AllDat$Standard_LengthCalcTL[which(AllDat$Standard_LengthCalcTL>0)]
AllDat$Standard_LengthCalc[which(AllDat$Standard_LengthCalcFL>0)]=AllDat$FLSLa[which(AllDat$Standard_LengthCalcFL>0)]+AllDat$FLSLb[which(AllDat$Standard_LengthCalcFL>0)]*AllDat$Standard_LengthCalcFL[which(AllDat$Standard_LengthCalcFL>0)]
AllDat$Standard_LengthCalc[which(AllDat$Standard_DMLength>0)]=AllDat$Standard_DMLength[which(AllDat$Standard_DMLength>0)]
AllDat$Standard_LengthCalc[which(AllDat$Standard_FLength>0 & is.na(AllDat$FLSLb)&is.na(AllDat$TLSLb))]=AllDat$Standard_FLength[which(AllDat$Standard_FLength>0 &is.na(AllDat$FLSLb)&is.na(AllDat$TLSLb))]
AllDat$Standard_LengthCalc[which(AllDat$Standard_TLength>0 & is.na(AllDat$FLSLb)&is.na(AllDat$TLSLb))]=AllDat$Standard_TLength[which(AllDat$Standard_TLength>0 &is.na(AllDat$FLSLb)&is.na(AllDat$TLSLb))]
AllDat$Standard_LengthCalc[which(AllDat$Standard_FLength>0 & is.na(AllDat$FLSLb)&is.na(AllDat$FLTLb))]=AllDat$Standard_FLength[which(AllDat$Standard_FLength>0 &is.na(AllDat$FLSLb)&is.na(AllDat$FLTLb))]
AllDat$Standard_LengthCalc[which(AllDat$Standard_TLength>0 & is.na(AllDat$TLSLb)&is.na(AllDat$TLFLb))]=AllDat$Standard_TLength[which(AllDat$Standard_TLength>0 &is.na(AllDat$TLSLb)&is.na(AllDat$TLFLb))]

#Calculate the standard deviation for abundance when using prevalence * intensity by propagating the standard deviation of the prevalence (n*p*(1-p)) and intensity
AllDat$Calc_StdIDev=NA
AllDat$Calc_StdIDev=AllDat$Standard_Abundance*sqrt((sqrt(AllDat$Prevalence/100*AllDat$Total_Fish_Examined*(1-AllDat$Prevalence/100))/AllDat$Prevalence/100)^(2)+(AllDat$Standard_IntensityStdDev/AllDat$Intensity_Mean)^(2))
AllDat$Calc_StdIDev[which(AllDat$Prevalence==0 & AllDat$Standard_IntensityStdDev)]=0

#Bring all the errors together 
AllDat$Standard_StandardDev=AllDat$Standard_AbundanceStdDev
AllDat$Standard_StandardDev[which(AllDat$Calc_StdDev>=0)]=AllDat$Calc_StdDev[which(AllDat$Calc_StdDev>=0)]
AllDat$Standard_StandardDev[which(AllDat$Calc_StdIDev>=0)]=AllDat$Calc_StdIDev[which(AllDat$Calc_StdIDev>=0)]
AllDat$Standard_StandardDev[which(AllDat$Standard_Abundance==0)]=0

AllDat2=AllDat[is.na(AllDat$Implicit),]

#count the total number of hosts in the dataset across both anisakid genera

head(AllDat2)
length(levels(AllDat2$Host[which(AllDat2$Wild_Farmed=='Wild')]))
#225 minus the ten that have no data = 215


#Subset the data to just Anisakis sp
Anisakis=AllDat2[grep('Anisakis',AllDat2$Parasite),]

#Aggregate the Anisakis data to the level of Anisakis sp
AnisakisAll2=aggregate(cbind(Standard_Abundance,Standard_StandardDev^2)~TrueHost+Method_Counting_parasites+Total_Fish_Examined+Latitude_decimaldegrees+Longitude_decimaldegrees+Standard_Year+Standard_LengthCalc+Portion_of_Body+ID+Wild_Farmed,data = Anisakis,FUN='sum',na.rm=F)
AnisakisAll2$Standard_StandardDev=sqrt(AnisakisAll2$V2) #Propagate error of the individual parasite species through to the aggregated (assuming independence between the species)

#Append FAO major regions to the Anisakis dataset by matching the lat and longs from the two data sets
AnisakisAll2$FAO=FAO$F_CODE[match(paste(AnisakisAll2$Latitude_decimaldegrees,AnisakisAll2$Longitude_decimaldegrees), paste(FAO$Latitude_dd,FAO$Longitude_dd))]
AnisakisAll2$FAO[is.na(AnisakisAll2$FAO)]=0
AnisakisAll2$FAO=as.factor(AnisakisAll2$FAO)#convert the numeric FAO zones to a factor


# Need to fix the messed-up FAO regions that come out as 0. These are points that border too closely on land.
# I did this by hand by printing out a csv, looking up lats and longs, and correcting 0s to appropriate numbers.

# write.csv(AnisakisAll2,"AnisakisAll2.csv")

# Once I fixed the 0s by hand, I read the updated file back in.

AnisakisAll2<-read.csv("data/processed_data_anisakis.csv",header=T,sep=",")



#Subset the data to just Pseudoterranova sp
Pseudoterranova=AllDat2[grep('Pseudoterranova',AllDat2$Parasite),]


#Aggregate the Pseudoterranova data to the level of Pseudoterranova sp
PseudoterranovaAll2=aggregate(cbind(Standard_Abundance,Standard_StandardDev^2)~TrueHost+Method_Counting_parasites+Total_Fish_Examined+Latitude_decimaldegrees+Longitude_decimaldegrees+Standard_Year+Standard_LengthCalc+Portion_of_Body+ID+Wild_Farmed,data = Pseudoterranova,FUN='sum',na.rm=F)
PseudoterranovaAll2$Standard_StandardDev=sqrt(PseudoterranovaAll2$V2) #Propagate error of the individual parasite species through to the aggregated (assuming independence between the species)

#Append FAO major regions to the Pseudoterranova dataset by matching the lat and longs from the two data sets
PseudoterranovaAll2$FAO=FAO$F_CODE[match(paste(PseudoterranovaAll2$Latitude_decimaldegrees,PseudoterranovaAll2$Longitude_decimaldegrees), paste(FAO$Latitude_dd,FAO$Longitude_dd))]
PseudoterranovaAll2$FAO[is.na(PseudoterranovaAll2$FAO)]=0
PseudoterranovaAll2$FAO=as.factor(PseudoterranovaAll2$FAO)#convert the numeric FAO zones to a factor


# Need to fix the messed-up FAO regions that come out as 0. These are points that border too closely on land.
# I did this by hand by printing out a csv, looking up lats and longs, and correcting 0s to appropriate numbers.

# write.csv(PseudoterranovaAll2,"PseudoterranovaAll2.csv")

# Once I fixed the 0s by hand, I read the updated file back in.

PseudoterranovaAll2<-read.csv("data/processed_data_pseudoterranova.csv",header=T,sep=",")





###### ANISAKIS ANALYSIS

AAll=escalc(measure='GEN', yi=(Standard_Abundance)^(1/4),vi=.25*Standard_StandardDev/(Standard_Abundance+1)*(Standard_Abundance^(1/4)),ni=Total_Fish_Examined,data=AnisakisAll2)

Amodtime.abun.sd=rma.mv(yi=yi,V=vi+1,mods=~Standard_Year+Standard_LengthCalc,random=list(~1|TrueHost/Portion_of_Body,~1|FAO,~1|Method_Counting_parasites,~1|ID),method='ML',data=AAll[which(AAll$Wild_Farmed=='Wild'),],control=list(optimizer='optim'))

Amodtime.abun.sd.null=rma.mv(yi=yi,V=vi+1,mods=~Standard_LengthCalc,random=list(~1|TrueHost/Portion_of_Body,~1|FAO,~1|Method_Counting_parasites,~1|ID),method='ML',data=AAll[which(AAll$Wild_Farmed=='Wild'),],control=list(optimizer='optim'))

Amodtime.outlier.removed=rma.mv(yi=yi,V=vi+1,mods=~Standard_Year+Standard_LengthCalc,random=list(~1|TrueHost/Portion_of_Body,~1|FAO,~1|Method_Counting_parasites,~1|ID),method='ML',data=AAll[which(AAll$Wild_Farmed=='Wild' & AAll$Standard_Year >1962),],control=list(optimizer='optim'))

summary(Amodtime.abun.sd)
summary(Amodtime.abun.sd.null)
summary(Amodtime.outlier.removed)

#count up the number of host species and FAO regions and the number of data points obtained by each method
#for the beginning of the results section
length(AAll$FAO[which(AAll$Wild_Farmed=='Wild')])
levels(as.factor(AAll$FAO[which(AAll$Wild_Farmed=='Wild')]))
length(levels(AAll$TrueHost[which(AAll$Wild_Farmed=='Wild')]))
tapply(AAll$Method_Counting_parasites[which(AAll$Wild_Farmed=='Wild')],AAll$Method_Counting_parasites[which(AAll$Wild_Farmed=='Wild')],length)

#calculate a pseudo-R2 for this model according to https://stackoverflow.com/questions/22356450/getting-r-squared-from-a-mixed-effects-multilevel-model-in-metafor
(sum(Amodtime.abun.sd.null$sigma2)-sum(Amodtime.abun.sd$sigma2))/sum(Amodtime.abun.sd.null$sigma2)


###### Test for species-level drivers in Anisakis
AnisHost=levels(AAll$TrueHost)
Host.Influence=data.frame(Host=AnisHost,TimeEffect=NA)

for (x in c(1:61,64:89,91:127,129:174,176:length(AnisHost))){
  print(x)
  Amodtime2.abun.sd=rma.mv(yi=yi,V=vi+1,mods=~Standard_Year+Standard_LengthCalc,random=list(~1|TrueHost/Portion_of_Body,~1|FAO,~1|Method_Counting_parasites,~1|ID),method='ML',data=AAll[which(AAll$Wild_Farmed=='Wild'& AAll$TrueHost!=AnisHost[x]),],control=list(optimizer='optim'))
  Host.Influence$TimeEffect[x]=Amodtime2.abun.sd$b[2]
  Host.Influence$p[x]=Amodtime2.abun.sd$pval[2]
  
}

Host.Influence$Diff=Host.Influence$TimeEffect-Amodtime.abun.sd$b[2]
Host.Influence[which(Host.Influence$Diff==max(Host.Influence$Diff,na.rm=T)),]
Host.Influence[which(Host.Influence$Diff==min(Host.Influence$Diff,na.rm=T)),]

#quantify the influence of Aphanopus carbo

Amodtime.noAcarbo=rma.mv(yi=yi,V=vi+1,mods=~Standard_Year+Standard_LengthCalc,random=list(~1|TrueHost/Portion_of_Body,~1|FAO,~1|Method_Counting_parasites,~1|ID),method='ML',data=AAll[which(AAll$Wild_Farmed=='Wild'& AAll$TrueHost!='Aphanopus carbo'),],control=list(optimizer='optim'))
summary(Amodtime.noAcarbo)

Amodtime.regular=rma.mv(yi=yi,V=vi+1,mods=~Standard_Year+Standard_LengthCalc,random=list(~1|TrueHost/Portion_of_Body,~1|FAO,~1|Method_Counting_parasites,~1|ID),method='ML',data=AAll[which(AAll$Wild_Farmed=='Wild'),],control=list(optimizer='optim'))
summary(Amodtime.regular)


###### Test for FAO-region drivers in Anisakis

AnisFAO=levels(as.factor(AAll$FAO))
FAO.Influence=data.frame(FAO=AnisFAO,TimeEffect=NA)
Amodtime2.abun.sd=rma.mv(yi=yi,V=vi+1,mods=~Standard_Year+Standard_LengthCalc,random=list(~1|TrueHost/Portion_of_Body,~1|FAO,~1|Method_Counting_parasites,~1|ID),method='ML',data=AAll[which(AAll$Wild_Farmed=='Wild'),],control=list(optimizer='optim'))

for (x in c(1:length(AnisFAO))){
  print(x)
  Amodtime2.abun.sd.FAO=rma.mv(yi=yi,V=vi+1,mods=~Standard_Year+Standard_LengthCalc,random=list(~1|TrueHost/Portion_of_Body,~1|FAO,~1|Method_Counting_parasites,~1|ID),method='ML',data=AAll[which(AAll$Wild_Farmed=='Wild'& AAll$FAO!=AnisFAO[x]),],control=list(optimizer='optim'))
  FAO.Influence$TimeEffect[x]=Amodtime2.abun.sd.FAO$b[2]
  FAO.Influence$p[x]=Amodtime2.abun.sd.FAO$pval[2]
  
}

FAO.Influence$Diff=FAO.Influence$TimeEffect-Amodtime2.abun.sd$b[2]
FAO.Influence[which(FAO.Influence$Diff==max(FAO.Influence$Diff,na.rm=T)),]
FAO.Influence[which(FAO.Influence$Diff==min(FAO.Influence$Diff,na.rm=T)),]

#result without region 27
Amodtime2.abun.sd.FAO=rma.mv(yi=yi,V=vi+1,mods=~Standard_Year+Standard_LengthCalc,random=list(~1|TrueHost/Portion_of_Body,~1|FAO,~1|Method_Counting_parasites,~1|ID),method='ML',data=AAll[which(AAll$Wild_Farmed=='Wild'& AAll$FAO!=AnisFAO[2]),],control=list(optimizer='optim'))

#does region 27 contain a lot of our records?
thing<-AAll[which(AAll$Wild_Farmed=='Wild'),]
thingy<-tapply(thing$FAO,thing$FAO,FUN=length)
as.data.frame(thingy)
sum(thingy)


###### Test for detection-technique drivers in Anisakis

AnisDetect=levels(as.factor(AAll$Method_Counting_parasites))

# Exclude any study that used digestion
Amodtime2.abun.sd.digest=rma.mv(yi=yi,V=vi+1,mods=~Standard_Year+Standard_LengthCalc,random=list(~1|TrueHost/Portion_of_Body,~1|FAO,~1|Method_Counting_parasites,~1|ID),method='ML',data=AAll[which(AAll$Wild_Farmed=='Wild'& AAll$Method_Counting_parasites!="Digestion" & AAll$Method_Counting_parasites!="Digestion and Microscopy" & AAll$Method_Counting_parasites!="Digestion and UV"),],control=list(optimizer='optim'))
summary(Amodtime2.abun.sd.digest)  

# Exclude any study that used UV
Amodtime2.abun.sd.uv=rma.mv(yi=yi,V=vi+1,mods=~Standard_Year+Standard_LengthCalc,random=list(~1|TrueHost/Portion_of_Body,~1|FAO,~1|Method_Counting_parasites,~1|ID),method='ML',data=AAll[which(AAll$Wild_Farmed=='Wild'& AAll$Method_Counting_parasites!="Digestion and UV" & AAll$Method_Counting_parasites!="Microscopyand UV" & AAll$Method_Counting_parasites!="UV"),],control=list(optimizer='optim'))
summary(Amodtime2.abun.sd.uv)  




###### PSEUDOTERRANOVA ANALYSIS
PAll=escalc(measure='MN', mi=Standard_Abundance^(1/4),sdi=.25*Standard_StandardDev/(Standard_Abundance+1)*(Standard_Abundance^(1/4)),ni=Total_Fish_Examined,data=PseudoterranovaAll2)

Pmodtime.abun.corsd=rma.mv(yi=yi,V=vi+1,mods=~Standard_Year+Standard_LengthCalc,random=list(~1|TrueHost/Portion_of_Body,~1|FAO,~1|Method_Counting_parasites,~1|ID),data=PAll[which(PAll$Wild_Farmed=='Wild'),],method = 'ML',control=list(optimizer='optim'))

Pmodtime.abun.null=rma.mv(yi=yi,V=vi+1,mods=~1,random=list(~1|TrueHost/Portion_of_Body,~1|FAO,~1|Method_Counting_parasites,~1|ID),data=PAll[which(PAll$Wild_Farmed=='Wild'),],method = 'ML',control=list(optimizer='optim'))

summary(Pmodtime.abun.corsd)
summary(Pmodtime.abun.null)

#calculate a pseudo-R2 for this model according to https://stackoverflow.com/questions/22356450/getting-r-squared-from-a-mixed-effects-multilevel-model-in-metafor
(sum(Pmodtime.abun.null$sigma2)-sum(Pmodtime.abun.corsd$sigma2)/sum(Pmodtime.abun.null$sigma2))

#count up the number of host species and FAO regions and the number of data points obtained by each method
#for the beginning of the results section
length(PAll[which(PAll$Wild_Farmed=='Wild'),])
levels(as.factor(PAll$FAO[which(PAll$Wild_Farmed=='Wild')]))
length(levels(PAll$TrueHost[which(PAll$Wild_Farmed=='Wild')]))
tapply(PAll$Method_Counting_parasites[which(PAll$Wild_Farmed=='Wild')],PAll$Method_Counting_parasites[which(PAll$Wild_Farmed=='Wild')],length)


###### Test for species-level drivers in Pseudoterranova

PseHost=levels(PAll$TrueHost)
PseHost.Influence=data.frame(Host=PseHost,TimeEffect=NA,p=NA)

for (x in c(1:31,32,33:52,53,54,56:62,63,64:72,73,74:83,85:length(PseHost))){
  print(x)
  Pmodtime2.abun.sd=rma.mv(yi=yi,V=vi+1,mods=~Standard_Year+Standard_LengthCalc,random=list(~1|TrueHost/Portion_of_Body,~1|FAO,~1|Method_Counting_parasites,~1|ID),method='ML',data=PAll[which(PAll$Wild_Farmed=='Wild'& PAll$TrueHost!=PseHost[x]),],control=list(optimizer='optim'))
  PseHost.Influence$TimeEffect[x]=Pmodtime2.abun.sd$b[2]
  PseHost.Influence$p[x]=Pmodtime2.abun.sd$pval[2]
}

PseHost.Influence$Diff=PseHost.Influence$TimeEffect-Pmodtime.abun.corsd$b[2]
PseHost.Influence[which(PseHost.Influence$Diff==max(PseHost.Influence$Diff,na.rm = T)),]
PseHost.Influence[which(PseHost.Influence$Diff==min(PseHost.Influence$Diff,na.rm = T)),]


###### Test for FAO-region drivers in Anisakis

PseFAO=levels(as.factor(PAll$FAO))
PseFAO.Influence=data.frame(FAO=PseFAO,TimeEffect=NA,p=NA)
Pmodtime2.abun.sd=rma.mv(yi=yi,V=vi+1,mods=~Standard_Year+Standard_LengthCalc,random=list(~1|TrueHost/Portion_of_Body,~1|FAO,~1|Method_Counting_parasites,~1|ID),method='ML',data=PAll[which(PAll$Wild_Farmed=='Wild'),],control=list(optimizer='optim'))

for (x in c(1:length(PseFAO))){
  print(x)
  Pmodtime2.abun.sd.FAO=rma.mv(yi=yi,V=vi+1,mods=~Standard_Year+Standard_LengthCalc,random=list(~1|TrueHost/Portion_of_Body,~1|FAO,~1|Method_Counting_parasites,~1|ID),method='ML',data=PAll[which(PAll$Wild_Farmed=='Wild'& PAll$FAO!=PseFAO[x]),],control=list(optimizer='optim'))
  PseFAO.Influence$TimeEffect[x]=Pmodtime2.abun.sd.FAO$b[2]
  PseFAO.Influence$p[x]=Pmodtime2.abun.sd.FAO$pval[2]
}

PseFAO.Influence$Diff=PseFAO.Influence$TimeEffect-Pmodtime2.abun.sd$b[2]
PseFAO.Influence[which(PseFAO.Influence$Diff==max(PseFAO.Influence$Diff,na.rm = T)),]
PseFAO.Influence[which(PseFAO.Influence$Diff==min(PseFAO.Influence$Diff,na.rm = T)),]


###### Test for detection-technique drivers in Pseudoterranova

PseDetect=levels(as.factor(PAll$Method_Counting_parasites))

# Exclude any study that used digestion
Pmodtime2.abun.sd.digest=rma.mv(yi=yi,V=vi+1,mods=~Standard_Year+Standard_LengthCalc,random=list(~1|TrueHost/Portion_of_Body,~1|FAO,~1|Method_Counting_parasites,~1|ID),method='ML',data=PAll[which(PAll$Wild_Farmed=='Wild'& PAll$Method_Counting_parasites!="Digestion and Microscopy"),],control=list(optimizer='optim'))
summary(Pmodtime2.abun.sd.digest)  

# Exclude any study that used UV
Pmodtime2.abun.sd.uv=rma.mv(yi=yi,V=vi+1,mods=~Standard_Year+Standard_LengthCalc,random=list(~1|TrueHost/Portion_of_Body,~1|FAO,~1|Method_Counting_parasites,~1|ID),method='ML',data=PAll[which(PAll$Wild_Farmed=='Wild'& PAll$Method_Counting_parasites!="Microscopyand UV"),],control=list(optimizer='optim'))
summary(Pmodtime2.abun.sd.uv) 







#### MAKE PLOTS


# FIGURE 2
mapWorld <- borders("world", colour="grey90", fill="grey90")#get some data to plot a world map 
Anismap=ggplot()+mapWorld+geom_bin2d(aes(x=Longitude_decimaldegrees,y=Latitude_decimaldegrees),binwidth=c(4,4),data = AnisakisAll2)+theme_minimal()+labs(y='',x='',fill='Count')+scale_fill_gradientn(colors = c("#00DEFF","#00A0FC","#0366FB","#2124A8","#093164","#093164"),limits=c(-5,55),breaks=c(0,10,20,30,40,50),name="Number of records")+ theme(panel.background = element_rect(fill = NA, colour = "black"))+scale_x_continuous(breaks = c(-180,-120,-60,0,60,120,180),limits=c(-180,180))+scale_y_continuous(breaks = c(-90,-60,-30,0,30,60,90),limits=c(-90,90))
Pseumap=ggplot()+mapWorld+geom_bin2d(aes(x=Longitude_decimaldegrees,y=Latitude_decimaldegrees),binwidth=c(4,4),data = PseudoterranovaAll2)+theme_minimal()+labs(y='',x='',fill='Count')+scale_fill_gradientn(colors = c("#00DEFF","#00A0FC","#0366FB","#2124A8","#093164","#093164"),limits=c(-5,55),breaks=c(0,10,20,30,40,50))+ theme(panel.background = element_rect(fill = NA, colour = "black"))+scale_x_continuous(breaks = c(-180,-120,-60,0,60,120,180),limits=c(-180,180))+scale_y_continuous(breaks = c(-90,-60,-30,0,30,60,90),limits=c(-90,90))


# FIGURE 3
A.Mod.Graph=data.frame(Standard_Year=1962:2015)
Apred=data.frame(predict.rma(Amodtime.abun.sd,newmods = cbind(A.Mod.Graph$Standard_Year,0)))
A.Mod.Graph$Standard_Abundance=(Apred$pred)^4
A.Mod.Graph$LCL=(Apred$ci.lb)^4
A.Mod.Graph$UCL=(Apred$ci.ub)^4
anisplot<-ggplot()+
  geom_point(size=2,data=AAll[which(AAll$Wild_Farmed=='Wild'),],aes(x=Standard_Year,y=(Standard_Abundance/Standard_LengthCalc),color=TrueHost,shape=Portion_of_Body))+
  theme_minimal()+ 
  scale_color_discrete(guide=F)+
  theme(legend.position = 'top',legend.text=element_text(size=16),legend.title=element_text(size=18))+
  geom_line(data=A.Mod.Graph,aes(x=Standard_Year,y=Standard_Abundance),color="#2171b5",size=2)+
  geom_ribbon(data=A.Mod.Graph,aes(x=Standard_Year,ymin=LCL,ymax=UCL),alpha=.5,fill="#B2DFEE")+
  ylab('')+
  xlab('')+
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"))+
  scale_shape_manual(name="Portion of host examined",labels=c('Alimentary tract','Filet','Viscera','Whole'),values=c(15,16,17,18))

#estimate magnitude of change over time by going into A.Mod.Graph and pulling out estimates for the appropriate year
A.Mod.Graph
anisakis_standardized_change<-1.035568e+00/3.655275e-03
anisakis_standardized_change

P.Mod.Graph=data.frame(Standard_Year=1978:2015)
Ppred=data.frame(predict.rma(Pmodtime.abun.corsd,newmods = cbind(P.Mod.Graph$Standard_Year,0)))
P.Mod.Graph$Standard_Abundance=(Ppred$pred)^4
P.Mod.Graph$LCL=(Ppred$ci.lb)^4
P.Mod.Graph$UCL=(Ppred$ci.ub)^4
pseuplot<-ggplot()+
  geom_point(data=PseudoterranovaAll2,aes(x=Standard_Year,y=Standard_Abundance/Standard_LengthCalc,color=TrueHost,shape=Portion_of_Body))+ 
  theme_minimal()+
  theme(legend.position='none')+
  geom_line(data=P.Mod.Graph,aes(x=Standard_Year,y=Standard_Abundance),color='#ff6961')+
  geom_ribbon(data=P.Mod.Graph,aes(x=Standard_Year,ymin=LCL,ymax=UCL),alpha=.2,fill='#ff6961')+
  ylab('')+
  xlab('Year')+
  theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"))+
  scale_shape_manual(name="Portion of Host Examined",labels=c('Alimentary tract','Filet','Viscera','Whole'),values=c(15,16,17,18))

library(cowplot)
ggdraw(plot=NULL,xlim=c(0,10),ylim=c(0,10))+
  draw_plot(anisplot,x=0,y=5,width=10,height=5)+
  draw_plot(pseuplot,x=2.77,y=0,width=7.1,height=5)+
  draw_text("Abundance (average number of worms per fish, length standardized)",x=0.15,y=0.4,size=18,hjust=0,vjust=0,angle=90)+
  draw_text("A",size=20,x=0.5,y=9.5)+
  draw_text("B",size=20,x=0.5,y=4.5)


# SUPPLEMENTARY FIGURE S1
All_HostInfluence=Host.Influence
randomeffectsA=ranef(Amodtime.abun.sd)
colnames(All_HostInfluence)=c('Host','Anisakis_TimeEffect','Anisakis_time_p','Anisakis_Host_Influence')
randomeffectsA$TrueHost$Host=row.names(randomeffectsA$TrueHost)
randomeffectsAnisakisHosts=randomeffectsA$TrueHost
All_HostInfluence$Anisakis_HostBurden=randomeffectsAnisakisHosts$intrcpt[match(All_HostInfluence$Host,randomeffectsAnisakisHosts$Host)]
All_HostInfluence$Anisakis_Host_Influence[is.na(All_HostInfluence$Anisakis_HostBurden)]=NA

write.csv(All_HostInfluence,'data/HostInfoTable.csv')

length(All_HostInfluence$Host)
thing0.111<-c(rep(1,93))
thing111.223<-(rep(2,92))
All_HostInfluence$group<-c(thing0.111,thing111.223)#Had to switch from rbind() to c() to get it to work properly
All_HostInfluence$Host<-gsub("Doryteuthis gahi\xca","Doryteuthis gahi",All_HostInfluence$Host)
All_HostInfluence$Host<-as.factor(All_HostInfluence$Host)
All_HostInfluence$Host=factor(All_HostInfluence$Host,levels=rev(levels(All_HostInfluence$Host)))

dev.off

first<-ggplot(aes(y=Host),data=All_HostInfluence[which(All_HostInfluence$group==1),])+
  geom_raster(aes(x='Anisakis Burden',fill=((Anisakis_HostBurden-mean(Anisakis_HostBurden,na.rm=T))/sd(Anisakis_HostBurden,na.rm=T))))+
  geom_raster(aes(x='Anisakis Influence',fill=-1*(Anisakis_Host_Influence-mean(Anisakis_Host_Influence,na.rm=T))/sd(Anisakis_Host_Influence,na.rm=T)))+
  scale_fill_gradient2(low="#d7191c",mid="white",high="#2c7bb6",midpoint=0,name='Relative effect',na.value='black')+
  theme_minimal()+
  theme(axis.title.x = element_blank(),legend.position = "none",text=element_text(size=8))

second<-ggplot(aes(y=Host),data=All_HostInfluence[which(All_HostInfluence$group==2),])+
  geom_raster(aes(x='Anisakis Burden',fill=((Anisakis_HostBurden-mean(Anisakis_HostBurden,na.rm=T))/sd(Anisakis_HostBurden,na.rm=T))))+
  geom_raster(aes(x='Anisakis Influence',fill=-1*(Anisakis_Host_Influence-mean(Anisakis_Host_Influence,na.rm=T))/sd(Anisakis_Host_Influence,na.rm=T)))+
  scale_fill_gradient2(low="#d7191c",mid="white",high="#2c7bb6",midpoint=0,name='Relative effect',na.value='black')+
  theme_minimal()+
  ylab("")+
  theme(axis.title.x = element_blank(),text=element_text(size=8))

library(cowplot)

ggdraw(plot=NULL,xlim=c(0,10),ylim=c(0,10))+
  draw_plot(first,x=0,y=0,width=5,height=10)+
  draw_plot(second,x=5,y=0,width=5,height=10)


# SUPPLEMENTARY FIGURE S2
All_FAOInfluence=FAO.Influence
colnames(All_FAOInfluence)=c('FAO_Region','Anisakis_TimeEffect','Anisakis_p','Anisakis_FAO_Influence')
All_FAOInfluence$Pseudoterranova_Influence=PseFAO.Influence$Diff[match(All_FAOInfluence$FAO,PseFAO.Influence$FAO)]
randomeffectsA$FAO$FAO=row.names(randomeffectsA$FAO)
randomeffectsAnisakisFAO=randomeffectsA$FAO
All_FAOInfluence$Anisakis_FAOBurden=randomeffectsAnisakisFAO$intrcpt[match(All_FAOInfluence$FAO,randomeffectsAnisakisFAO$FAO)]
All_FAOInfluence$Pseudoterranova_FAOBurden=randomeffectsP$FAO$intrcpt[match(All_FAOInfluence$FAO,row.names(randomeffectsP$FAO))]
All_FAOInfluence$Anisakis_FAO_Influence[is.na(All_FAOInfluence$Anisakis_FAO)]=NA
All_FAOInfluence$Pseudoterranova_Influence[is.na(All_FAOInfluence$Pseudoterranova_FAOBurden)]=NA
write.csv(All_FAOInfluence,'data/FAOInfoTable.csv')
All_FAOInfluence$FAO_Region<-factor(All_FAOInfluence$FAO_Region,levels=rev(levels(All_FAOInfluence$FAO_Region)))

ggplot(aes(y=FAO_Region),data=All_FAOInfluence)+
  geom_raster(aes(x='Anisakis Burden',fill=((Anisakis_FAOBurden-mean(Anisakis_FAOBurden,na.rm=T))/sd(Anisakis_FAOBurden,na.rm=T))))+
  geom_raster(aes(x='Anisakis Influence',fill=-1*((Anisakis_FAO_Influence-mean(Anisakis_FAO_Influence,na.rm=T))/sd(Anisakis_FAO_Influence,na.rm=T))))+
  geom_raster(aes(x='Pseudoterranova Burden',fill=((Pseudoterranova_FAOBurden-mean(Pseudoterranova_FAOBurden,na.rm=T))/sd(Pseudoterranova_FAOBurden,na.rm=T))))+
  geom_raster(aes(x='Pseudoterranova Influence',fill=-1*((Pseudoterranova_Influence-mean(Pseudoterranova_Influence,na.rm=T))/sd(Pseudoterranova_Influence,na.rm=T))))+
  scale_fill_gradient2(low="#d7191c",mid="white",high="#2c7bb6",midpoint=0,name='Relative effect',na.value='black')+
  scale_y_discrete()+
  ylab("FAO region")+
  theme_minimal()+
  theme(axis.title.x=element_blank())


# SUPPLEMENTARY FIGURE S3
All_HostInfluence=PseHost.Influence
randomeffectsP=ranef(Pmodtime.abun.corsd)
colnames(All_HostInfluence)=c('Host','Pse_TimeEffect','Pse_time_p','Pseudoterranova_Influence')
All_HostInfluence$Pseudoterranova_HostBurden=randomeffectsP$TrueHost$intrcpt[match(All_HostInfluence$Host,row.names(randomeffectsP$TrueHost))]
All_HostInfluence$Pseudoterranova_Influence[is.na(All_HostInfluence$Pseudoterranova_HostBurden)]=NA

length(All_HostInfluence$Host)
thing0.111<-c(rep(1,47))
thing111.223<-(rep(2,46))
All_HostInfluence$group<-c(thing0.111,thing111.223)#Had to switch from rbind() to c() to get it to work properly
All_HostInfluence$Host=factor(All_HostInfluence$Host,levels=rev(levels(All_HostInfluence$Host)))

dev.off()

first<-ggplot(aes(y=Host),data=All_HostInfluence[which(All_HostInfluence$group==1),])+
  geom_raster(aes(x='Pseudoterranova Burden',fill=(Pseudoterranova_HostBurden-mean(Pseudoterranova_HostBurden,na.rm=T))/sd(Pseudoterranova_HostBurden,na.rm=T)))+
  geom_raster(aes(x='Pseudoterranova Influence',fill=-1*(Pseudoterranova_Influence-mean(Pseudoterranova_Influence,na.rm=T))/sd(Pseudoterranova_Influence,na.rm=T)))+
  scale_fill_gradient2(low="#d7191c",mid="white",high="#2c7bb6",midpoint=0,name='Relative effect',na.value='black')+
  theme_minimal()+
  theme(axis.title.x = element_blank(),legend.position = "none",text=element_text(size=8))

second<-ggplot(aes(y=Host),data=All_HostInfluence[which(All_HostInfluence$group==2),])+
  geom_raster(aes(x='Pseudoterranova Burden',fill=(Pseudoterranova_HostBurden-mean(Pseudoterranova_HostBurden,na.rm=T))/sd(Pseudoterranova_HostBurden,na.rm=T)))+
  geom_raster(aes(x='Pseudoterranova Influence',fill=-1*(Pseudoterranova_Influence-mean(Pseudoterranova_Influence,na.rm=T))/sd(Pseudoterranova_Influence,na.rm=T)))+
  scale_fill_gradient2(low="#d7191c",mid="white",high="#2c7bb6",midpoint=0,name='Relative effect',na.value='black')+
  theme_minimal()+
  ylab("")+
  theme(axis.title.x = element_blank(),text=element_text(size=8))

ggdraw(plot=NULL,xlim=c(0,10),ylim=c(0,10))+
  draw_plot(first,x=0,y=0,width=5,height=10)+
  draw_plot(second,x=5,y=0,width=5,height=10)



