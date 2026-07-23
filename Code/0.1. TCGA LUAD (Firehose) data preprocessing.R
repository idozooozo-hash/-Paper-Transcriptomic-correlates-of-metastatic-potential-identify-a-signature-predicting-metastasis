########################################
## Directory
########################################
dir_fh="G:/My Drive/DB/Firehose/LUAD/" # The location where the TCGA(obtained from Firehose) data is stored
dir_tcga="G:/My Drive/DB/TCGA/mRNA_protein_coding/" # The location where the TCGA(obtained from Firehose) data is stored
dir_data="G:/My Drive/project/metastasis_potential_biomarker/manuscript/Revision(scientific_reports)/github/data/"



########################################
## Library
########################################
library(dplyr)
library(stringr)



########################################
## Load Data
########################################

## Sample information
load(paste0(dir_fh,"sinfo_merged.Rdata")) # sinfo
msinfo=sinfo
load(paste0(dir_fh,"sinfo_only_clin.Rdata")) # sinfo 
osinfo=sinfo
fsinfo=left_join(msinfo,osinfo)
rm(msinfo,osinfo)

## Read count
load(paste0(dir_fh,"read_count.Rdata")) # rcm

## TPM
load(paste0(dir_fh,"tpm.Rdata")) # tpm



########################################
## Dataset Preview
########################################
head(fsinfo$`patient.samples.sample.bcr_sample_barcode`)
# [1] "tcga-05-4245-01a" "tcga-05-4382-01a" "tcga-05-4384-01a" "tcga-05-4396-01a"
# [5] "tcga-05-4402-01a" "tcga-05-4405-01a"
head(fsinfo$`patient.bcr_patient_barcode`)
# [1] "tcga-05-4245" "tcga-05-4382" "tcga-05-4384" "tcga-05-4396" "tcga-05-4402"
# [6] "tcga-05-4405"

tpm[1:5,1:5]
#       TCGA-05-4244-01A-01R-1107-07 TCGA-05-4249-01A-01R-1107-07
# A1BG                     4.7565003                    6.9204708
# A1CF                     0.0000000                    0.4027222
# A2BP1                    1.4570164                    1.3839393
# A2LD1                    7.0927804                    6.4928610
# A2ML1                    0.4319973                    1.3839393
#       TCGA-05-4250-01A-01R-1107-07 TCGA-05-4382-01A-01R-1206-07
# A1BG                      5.696542                     7.198727
# A1CF                      0.000000                     0.000000
# A2BP1                     0.000000                     0.000000
# A2LD1                     7.249191                     6.821053
# A2ML1                     0.000000                     2.532591
#       TCGA-05-4384-01A-01R-1755-07
# A1BG                      7.004132
# A1CF                      0.000000
# A2BP1                     0.000000
# A2LD1                     6.468824
# A2ML1                     0.000000

rcm[1:5,1:5]
#       TCGA-05-4244-01A-01R-1107-07 TCGA-05-4244-01A-01R-1107-07.1
# A1BG                         74.57                   9.832852e-07
# A1CF                          0.00                   0.000000e+00
# A2BP1                         5.00                   3.549782e-08
# A2LD1                       388.18                   7.351586e-06
# A2ML1                         1.00                   5.431279e-09
#       TCGA-05-4244-01A-01R-1107-07.2 TCGA-05-4249-01A-01R-1107-07
# A1BG                              NA                       373.14
# A1CF                              NA                         1.00
# A2BP1                             NA                         5.00
# A2LD1                             NA                       276.63
# A2ML1                             NA                         5.00
#       TCGA-05-4249-01A-01R-1107-07.1
# A1BG                    6.370817e-06
# A1CF                    1.476505e-08
# A2BP1                   3.737857e-08
# A2LD1                   5.170925e-06
# A2ML1                   3.447295e-08



########################################
## sample information data
########################################

## only Tumor
nid=fsinfo$patient.bcr_patient_barcode[grep("Normal|normal",fsinfo$patient.samples.sample.sample_type)]
fsinfo=fsinfo[!(fsinfo$patient.bcr_patient_barcode %in% nid) , ]
sinfo=setNames(fsinfo[,c("patient.bcr_patient_barcode","patient.diagnosis")] , c("Sample","Cancer Type"))


## Age
sinfo$Age=as.numeric(fsinfo$patient.age_at_initial_pathologic_diagnosis)


## metastasis yes/no
# 참고 열 : patient.stage_event.tnm_categories.pathologic_categories.pathologic_m
minfo=fsinfo[ , colnames(fsinfo)[grep("new.*event_type",colnames(fsinfo)) ]]
yinfo=fsinfo[ , colnames(fsinfo)[setdiff(grep("new_tumor_event_after_initial_treatment",colnames(fsinfo)) , grep("days",colnames(fsinfo)))]]
mid=fsinfo$`patient.bcr_patient_barcode`[apply(minfo, 1, function(c) sum(grepl("metastasis",c)))>0 & apply(yinfo, 1, function(c) sum(grepl("yes",c)))>0] # meta id
nid=setdiff(fsinfo$`patient.bcr_patient_barcode`[apply(yinfo, 1, function(c) sum(grepl("no",c)))>0] , mid) # non meta id
sinfo$`Metastasis Ocurrence`=ifelse(fsinfo$patient.bcr_patient_barcode %in% mid , "Distant Metastasis" , ifelse(fsinfo$`patient.bcr_patient_barcode` %in% nid, "No Metastasis",NA))


## metastasis day info
ecols=colnames(fsinfo)[grep("new.*event_type",colnames(fsinfo)) ]
ecoll=apply(fsinfo[,ecols], 1, function(c) ecols[grepl("metastasis",c)])
dcoll=lapply(ecoll, function(col) ifelse(length(col)>0 , paste0(gsub("\\.new_neoplasm_event_types\\.new_neoplasm_event_type(-2)?$", "", col),".days_to_new_tumor_event_after_initial_treatment") , NA))
mdays=as.numeric(sapply(1:nrow(fsinfo), function(i) ifelse(!is.na(dcoll[[i]]), min(fsinfo[i,dcoll[[i]]],na.rm=T) , NA)))
fdays=as.numeric(apply(fsinfo[,grep("days_to_last_followup",colnames(fsinfo))] , 1, function(d) max(d,na.rm=T)))
sinfo$`Metastasis Free Survival`=ifelse(fsinfo$patient.bcr_patient_barcode %in% mid , mdays , ifelse(fsinfo$`patient.bcr_patient_barcode` %in% nid, fdays,NA))


## recurrence yes/no
rinfo=fsinfo[ , colnames(fsinfo)[grep("new.*event_type",colnames(fsinfo)) ]]
yinfo=fsinfo[ , colnames(fsinfo)[setdiff(grep("new_tumor_event_after_initial_treatment",colnames(fsinfo)) , grep("days",colnames(fsinfo)))]]
rid=fsinfo$`patient.bcr_patient_barcode`[apply(rinfo, 1, function(c) sum(grepl("recurrence",c)))>0 & apply(yinfo, 1, function(c) sum(grepl("yes",c)))>0] # meta id
nid=setdiff(fsinfo$`patient.bcr_patient_barcode`[apply(yinfo, 1, function(c) sum(grepl("no",c)))>0] , rid) # non meta id
sinfo$`Recurrence Ocurrence`=ifelse(fsinfo$patient.bcr_patient_barcode %in% rid , "Locoregional Recurrence" , ifelse(fsinfo$`patient.bcr_patient_barcode` %in% nid, "No Recurrence",NA))


## recurrence day info
ecols=colnames(fsinfo)[grep("new.*event_type",colnames(fsinfo)) ]
ecoll=apply(fsinfo[,ecols], 1, function(c) ecols[grepl("recurrence",c)])
dcoll=lapply(ecoll, function(col) ifelse(length(col)>0 , paste0(gsub("\\.new_neoplasm_event_types\\.new_neoplasm_event_type(-2)?$", "", col),".days_to_new_tumor_event_after_initial_treatment") , NA))
rdays=as.numeric(sapply(1:nrow(fsinfo), function(i) ifelse(!is.na(dcoll[[i]]), min(fsinfo[i,dcoll[[i]]],na.rm=T) , NA)))
fdays=as.numeric(apply(fsinfo[,grep("days_to_last_followup",colnames(fsinfo))] , 1, function(d) max(d,na.rm=T)))
sinfo$`Relapse Free Survival`=ifelse(fsinfo$patient.bcr_patient_barcode %in% rid , rdays , ifelse(fsinfo$`patient.bcr_patient_barcode` %in% nid, fdays,NA))


## survival yes/no
sinfo$`Vital Status`=str_to_title(fsinfo$patient.vital_status)
sdays=apply(fsinfo[ , grep("days_to_last_followup|days_to_last_known_alive",colnames(fsinfo))] , 1, function(c) max(c,na.rm=T))
ddays=apply(fsinfo[ , grep("days_to_last_followup|patient.days_to_death",colnames(fsinfo))] , 1, function(c) min(c,na.rm=T))
sinfo$`Overall Survival`=as.numeric(ifelse(sinfo$`Vital Status`=="Dead", ddays , ifelse(sinfo$`Vital Status`=="Alive" , sdays , NA)))


## stage
sinfo$`Pathologic Stage`=toupper(gsub("stage |a$|b$|c$","",fsinfo$patient.stage_event.pathologic_stage))
sinfo$`Clinical Stage`=toupper(gsub("stage |a$|b$|c$","",fsinfo$patient.omfs.omf.omf_staging.clinical_stage))


## drug
sinfo$`Drug Name`=str_to_title(fsinfo$patient.drugs.drug.drug_name)
sinfo$`Drug Type`=str_to_title(fsinfo$patient.drugs.drug.therapy_types.therapy_type)
sinfo$`Drug Response`=str_to_title(fsinfo$patient.drugs.drug.measure_of_response)


## radiation
rid=setdiff(fsinfo$`patient.bcr_patient_barcode`[apply(fsinfo[,c("patient.radiation_therapy","patient.follow_ups.follow_up.additional_radiation_therapy")], 1, function(c) sum(c=="yes"))>0],NA)
nid=setdiff(fsinfo$`patient.bcr_patient_barcode`[apply(fsinfo[,c("patient.radiation_therapy","patient.follow_ups.follow_up.additional_radiation_therapy")], 1, function(c) sum(c=="no"))>0],c(NA,rid))
sinfo$`Radiation Therapy`=ifelse(sinfo$Sample %in% rid, "Yes" , ifelse(sinfo$Sample %in% nid , "No" , NA))
sinfo$`Radiation Therapy Type`=ifelse(sinfo$Sample %in% rid, fsinfo$patient.radiations.radiation.radiation_type , NA)
sinfo$`Radiation Therapy Dosage`=as.numeric(ifelse(sinfo$Sample %in% rid, fsinfo$patient.radiations.radiation.radiation_dosage , NA))
sinfo$`Radiation Target`=ifelse(sinfo$Sample %in% rid, fsinfo$patient.radiations.radiation.anatomic_treatment_site , NA)
sinfo$`Radiation Response`=ifelse(sinfo$Sample %in% rid, fsinfo$patient.radiations.radiation.measure_of_response , NA)


## Race
sinfo$Race=fsinfo$patient.race_list.race


## Gender
sinfo$Gender=str_to_title(fsinfo$patient.gender)


## first collected sample
sum(duplicated(fsinfo$patient.bcr_patient_barcode)) # 0 
# ==> first collected tissue 구하는 과정은 생략

sinfo$Tissue=fsinfo$patient.samples.sample.bcr_sample_barcode



########################################
## Firehose expression data
########################################

## read count
rcm=rcm[,!grepl("\\.2$|\\.1$",colnames(rcm))]
colnames(rcm)=tolower(sapply(str_split(colnames(rcm),"-"), function(c) paste0(c[1:4],collapse="-")))
rcm=rcm[,colnames(rcm) %in% sinfo$Tissue]
colnames(rcm)=sinfo$Sample[match(colnames(rcm),sinfo$Tissue)]

dugs=unique(rownames(rcm)[duplicated(rownames(rcm))])
for(dug in dugs) {
    ercm=rcm[rownames(rcm)==dug,]    
    rcm=rbind(rcm[rownames(rcm)!=dug,] , matrix(apply(ercm,2,mean), nrow=1, dimnames=list(dug , colnames(ercm))))    
}


## tpm
tpm=tpm[,!grepl("\\.2$|\\.1$",colnames(tpm))]
colnames(tpm)=tolower(sapply(str_split(colnames(tpm),"-"), function(c) paste0(c[1:4],collapse="-")))
tpm=tpm[,colnames(tpm) %in% sinfo$Tissue]
colnames(tpm)=sinfo$Sample[match(colnames(tpm),sinfo$Tissue)]

dugs=unique(rownames(tpm)[duplicated(rownames(tpm))])
for(dug in dugs) {
    etpm=tpm[rownames(tpm)==dug,]    
    tpm=rbind(tpm[rownames(tpm)!=dug,] , matrix(apply(etpm,2,mean), nrow=1, dimnames=list(dug , colnames(etpm))))    
}



########################################
## Save
########################################

## Remain only common samples
csam=Reduce(intersect , list(sinfo$Sample,colnames(rcm),colnames(tpm)))
sinfo=sinfo[match(csam,sinfo$Sample),]
rcm=rcm[,match(csam,colnames(rcm))]
tpm=tpm[,match(csam,colnames(tpm))]


## save
# save(sinfo, file=paste0(dir_data,"0.1. sample information(TCGA-Firehose).Rdata"))
# save(rcm, file=paste0(dir_data,"0.1. read count(TCGA-Firehose).Rdata"))
# save(tpm, file=paste0(dir_data,"0.1. tpm(TCGA-Firehose).Rdata"))