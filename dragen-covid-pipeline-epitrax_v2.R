#!/usr/bin/env Rscript
# Author: OLP/01/2022
# Usage: Getting collection dates and other metadata from EpiTrax for samples with LabIds and Customer Ids
# Parameters: Only one: Run Name i.e. Rscript / . . . ./dragen-covid-pipeline-epitrax-v2.R UT-A01290-220113

args = commandArgs(trailingOnly=TRUE)

library("dplyr", warn.conflicts = FALSE)
library("grid")
#library("gridBase")
library("rbin")#
library("tidyverse")
library("lubridate", warn.conflicts = FALSE)
library("tibble")

args = commandArgs(trailingOnly=TRUE)
# arg<-c('UT-A01290-220113');
s<- paste('/Volumes/NGS/Analysis/covidseq/',args[1],sep="/");
ss<-gsub("\\.*.*/","", s);date<-ymd(substr(ss, 11, 16))

# reading covidseq_summary

loc<-paste(paste('/Volumes/NGS/Analysis/covidseq',args[1],sep="/"),"sumcovidseq",sep ="/")
covidseq_summary<-read.csv(paste(loc,"covidseq_summary.csv", sep="/"));covidseq_summary<-covidseq_summary[,c(1,9)]
colnames(covidseq_summary)[2] <- 'Test results'
cat("1.- Reading covidseq_summary from",loc,"\n")

sampleNumber<-covidseq_summary$sample;
# cleaning
sampleNumber<-gsub("\\NTC.*","", sampleNumber);sampleNumber<-gsub("\\Pos.*","", sampleNumber);
sampleNumber<-gsub("\\Emp.*","", sampleNumber);sampleNumber<-gsub("\\Ext.*","", sampleNumber);
sampleNumber<-sampleNumber[sampleNumber !=""]
covidseq_summary<-merge(sampleNumber,covidseq_summary, by.x = 'x', by.y = 'sample'); # cleaning covidseq_summary
colnames(covidseq_summary) <- c("sample","Test results")


# Accessing EpiTrax files (NAS): detecting last file created (NGS......csv and All_ncovid.......csv)
tmpshot <- fileSnapshot("/Volumes/IDGenomics_NAS/NextStrain/EpiTrax/", pattern=glob2rx("All*.csv"))#
s<-rownames(tmpshot$info[which.max(tmpshot$info$mtime),]);
tmpshot <- fileSnapshot("/Volumes/IDGenomics_NAS/NextStrain/EpiTrax/", pattern=glob2rx("NGS*.csv"))#
s2<-rownames(tmpshot$info[which.max(tmpshot$info$mtime),]);
cat("2.- These files have been detected from EpiTrax export as the last created: ",s, " and ",s2,"\n")


locc<-paste("/Volumes/IDGenomics_NAS/NextStrain/EpiTrax/",s,sep="")
nNGS<-read.csv(locc)
cls <- function() cat(rep("\n",100));
NGS<-read.csv(paste("/Volumes/IDGenomics_NAS/NextStrain/EpiTrax/",s2,sep=""))



# Selecting usefull columns from EpiTrax files
# 1.  from All . . .NoCovid
nNGS_Covid<-nNGS[,c(2,3,4,6,8,10)]; # (lastName, firstName, DOB, CollectionDate, sampleNumer, summitterId)
# 2. from NGS_covid
NGS_Covid<-NGS[,c(2,3,4,6,8,11)]; # (lastName, firstName, DOB, CollectionDate, sampleNumer, summitterId)

sampleNumber_nNGS<-nNGS_Covid$sampleNumber;
sampleNumber_NGS<-NGS_Covid$sampleNumber;sampleNumber<-as.character(NGS_Covid$sampleNumber)
NGS_Covid <-subset(NGS_Covid, select=-c(sampleNumber));
NGS_Covid<-add_column(NGS_Covid, sampleNumber, .after = 'submitterId');

# merging the two tables (All_No... and NGS) to have the EpiTrax data
EpiTrax_data<-bind_rows(NGS_Covid, nNGS_Covid);


# searching info for samples with Lab ID: ## inner join
Matrix1<-merge(EpiTrax_data,covidseq_summary, by.x='sampleNumber', by.y ='sample')
Matrix1<-Matrix1[,c(2,3,4,6,1,7)]; # all info for sample swith LbId are in Matrix1
#######################################################################################
# searching info for samples without LabID, i.e. customer ids IHC, ...etc. # Left join

LMM<-NGS_Covid %>% filter(grepl('/',submitterId));
sampleNumber2<-gsub("\\.*.*/","", LMM$submitterId );
LMM<- cbind(LMM,sampleNumber2); #collumn with Customer id at the end from EpiTrax data
covidseq_summary<-covidseq_summary[-c(grep("^[[:digit:]]", covidseq_summary$sample)), ];
LMM<-NGS_Covid %>% filter(grepl('/',submitterId));
sampleNumber2<-gsub("\\.*.*/","", LMM$submitterId );
LMM<- cbind(LMM,sampleNumber2); #
Matrix2<-merge(covidseq_summary, LMM, by.x ='sample', by.y = 'sampleNumber2', all.x = TRUE)
x1<- Matrix2$sample;y1 <- Matrix2$sampleNumber;y111<-as.character(y1);cc<-coalesce(y111,x1);
Matrix2<- cbind(Matrix2,cc);Matrix2<-Matrix2[,c(3,4,5,8,9,2)];
colnames(Matrix2)[5] <- 'sampleNumber';

ngs_file<-bind_rows(Matrix1, Matrix2);
ngs_file$'Test results'[is.na(ngs_file$'Test results')]= "No able to be sequenced";ngs_file[is.na(ngs_file)]= ""
ngs_file2<-cbind(ngs_file, 'Name of facility performing the test'='UPHL', 'ordering facility' ='UHPL',
                 'Name or description of the test ordered'="SARS-CoV-2 Whole Genome Sequencing",
                 'lab test date'=ymd(date),'lab test status' ="preliminary",'facility' ="UPHL")
ngs_file<-ngs_file2[,c(1,2,3,7,9,4,5,10,11,6,12)]

# Eliminate possible duplicated based on the sampleNumber
ngs_file<-ngs_file[!duplicated(ngs_file$sampleNumber),]

# save ngs file
# =============================================================
ngs_file_name<-paste(paste('ngs',args[1],sep="_"),'csv',sep=".");
write.csv(ngs_file,paste(loc,ngs_file_name,sep="/"),row.names=FALSE)
cat("The file",ngs_file_name, "should be completed and saved at" ,loc ,"\n")