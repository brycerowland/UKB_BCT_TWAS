library(data.table)
phe <- fread("/proj/yunligrp/UKBB_phen_29983/ukb9883.tab",header=T,data.table=F)
### pregnancy
names(phe)[grep(3140,names(phe))]
preg <- phe[!is.na(phe$f.3140.0.0) & phe$f.3140.0.0==1,1]

### HIV??
names(phe)[grep(23064,names(phe))]

names(phe)[grep(20003,names(phe))]
### cancer treatments/drugs: 1140869524 - 1140870284; 
#1141157448 epoetin alfa product
#1140870568         erythropoietin product
#1140870600         epoetin alfa
#1140870604         eprex 4000iu/1ml injection
#1140870618         epoetin beta
#1140870626         recormon 1000iu injection+diluent
#1141178816         aranesp 10micrograms/0.4ml prefilled syringe
#1141178858         darbepoetin alfa

drugs <- phe[,c(1,2439:2486)]
drugs.excl <- apply(drugs[,-1],1,function(x)any(x%in%c(1140869523:1140870285,1141157448,1140870568,1140870600,1140870604,1140870618,1140870626,1141178816,1141178858),na.rm=T))
drugs.excl <- drugs[drugs.excl,1]

### cancer code
names(phe)[grep(20001,names(phe))]

quit()
# [1] "f.20001.0.0" "f.20001.0.1" "f.20001.0.2" "f.20001.0.3" "f.20001.0.4"
# [6] "f.20001.0.5" "f.20001.1.0" "f.20001.1.1" "f.20001.1.2" "f.20001.1.3"
#[11] "f.20001.1.4" "f.20001.1.5" "f.20001.2.0" "f.20001.2.1" "f.20001.2.2"
#[16] "f.20001.2.3" "f.20001.2.4" "f.20001.2.5"

ccode <- c(1085,1074,1070,1063,1058,1056,1055,1053,1052,1051,1050,1048,1047)
c20001 <- phe[,c(1,2319:2324)]
c20001.excl <- apply(c20001[,-1],1,function(x)any(x%in%ccode,na.rm=T))
c20001.excl <- c20001[c20001.excl,1]

### ICD10 main/secondary diag
code19 <- fread("coding19.tsv",header=T,data.table=F)

icd10.main <- c("B20","B21","B22","B23","B24","C40","C41","C81","C82","C83","C84","C85","C86","C87","C88","C89","C90","C91","C92","C93","C94","C95","C96","D45","D46","D47","D55","D56","D57","D58","D59","D60","D61","D63",
              "D640","D641","D642","D643","D644","D65","D66","D67","D68","D69","D70","D71","D72","D73","D74","D75","D76","D77","D80","D81","D82","D83","D84","D85","D86","D87","D88","D89","K70","K71","K74",
              "R70","R71","R72","R73","R74","R75","R76","R77","R78","R79")
icd10 <- NULL

for (i in icd10.main) icd10 <- c(icd10,code19[grep(i,code19[,1]),1])
icd10 <- icd10[-grep("Block",icd10)]

#>  grep(41202,names(phe))
# [1] 3905 3906 3907 3908 3909 3910 3911 3912 3913 3914 3915 3916 3917 3918 3919
#[16] 3920 3921 3922 3923 3924 3925 3926 3927 3928 3929 3930 3931 3932 3933 3934
#[31] 3935 3936 3937 3938 3939 3940 3941 3942 3943 3944 3945 3946 3947 3948 3949
#[46] 3950 3951 3952 3953 3954 3955 3956 3957 3958 3959 3960 3961 3962 3963 3964
#[61] 3965 3966 3967 3968 3969 3970

c41202 <- phe[,c(1,3905:3970)]
c41202.excl <- apply(c41202[,-1],1,function(x)any(x%in%icd10,na.rm=T))
c41202.excl <- c41202[c41202.excl,1]

c41204 <- phe[,c(1,3999:4182)]
c41204.excl <- apply(c41204[,-1],1,function(x)any(x%in%icd10,na.rm=T))
c41204.excl <- c41204[c41204.excl,1]

### ICD9 main/secondary diag
code87 <- fread("coding87.tsv",header=T,data.table=F)

icd9.main <- c(170,200:208,2384,2385,2386,2387,286,287,288,289,282,283,284,571,790)
icd9 <- code87[c(1205:1215,1403:1460,1675:1678,2194:2259,2271:2337,3958:3972,7540:7550),1]

c41203 <- phe[,c(1,3971:3998)]
c41203.excl <- apply(c41203[,-1],1,function(x)any(x%in%icd9,na.rm=T))
c41203.excl <- c41203[c41203.excl,1]

c41205 <- phe[,c(1,4183:4212)]
c41205.excl <- apply(c41205[,-1],1,function(x)any(x%in%icd9,na.rm=T))
c41205.excl <- c41205[c41205.excl,1]

### operative procedures opcs4
code240 <- fread("coding240.tsv",header=T,data.table=F)
op <- code240[c(3605:3623,8434:8442),1]

c41200 <- phe[,c(1,3834:3882)]
c41200.excl <- apply(c41200[,-1],1,function(x)any(x%in%op,na.rm=T))
c41200.excl <- c41200[c41200.excl,1]

c41210 <- phe[,c(1,4261:4346)]
c41210.excl <- apply(c41210[,-1],1,function(x)any(x%in%op,na.rm=T))
c41210.excl <- c41210[c41210.excl,1]

excl <- unique(c(preg,drugs.excl,c20001.excl,c41202.excl,c41204.excl,c41203.excl,c41205.excl,c41200.excl,c41210.excl))
write(excl,"BCratios_exclusion.txt",ncol=1)

phe <- phe[!phe[,1]%in%excl,] ###466057

######
#                   African Any other Asian background
#                      3008                       1704
#Any other Black background Any other mixed background
#                       109                        970
#Any other white background     Asian or Asian British
#                     15279                         39
#               Bangladeshi     Black or Black British
#                       214                         22
#                   British                  Caribbean
#                    410906                       3931
#                   Chinese                Do not know
#                      1485                        203
#                    Indian                      Irish
#                      5468                      12174
#                     Mixed         Other ethnic group
#                        45                       4216
#                 Pakistani       Prefer not to answer
#                      1693                       1513
#                     White            White and Asian
#                       528                        791
#   White and Black African  White and Black Caribbean
#                       384                        557
#######

phe <- phe[(!is.na(phe$f.21000.0.0)) & (phe$f.21000.0.0%in%c("Any other white background","British","Irish","White")),] ###438887
phe <- phe[is.na(phe$f.22010),] ### poor heterozygosity/missingness 438458
excl1 <- read.table("/gpfs/gsfs7/users/NHLBI_UKBB/test_MPV/bolt.in_plink_but_not_imputed.FID_IID.968.txt")
phe <- phe[!phe[,1]%in%excl1,] ### 438458

### ID, RBC, HGB, HCT, MCV, MCH, MCHC, RDW, Reticulocyte count, PLT, MPV, PDW, WBC, LYM, MON, NEU, EOS, BAS, htrc, mrv, msv, sex, age, PC1-40
#bcr <- na.omit(phe[,c("f.eid","f.30010.0.0","f.30020.0.0","f.30030.0.0","f.30040.0.0","f.30050.0.0","f.30060.0.0","f.30070.0.0","f.30250.0.0",
#"f.30080.0.0","f.30100.0.0","f.30110.0.0","f.30000.0.0","f.30120.0.0","f.30130.0.0",
#"f.30140.0.0","f.30150.0.0","f.30160.0.0","f.30300.0.0","f.30260.0.0","f.30270.0.0","f.31.0.0","f.21022.0.0",
#"f.22009.0.1","f.22009.0.2","f.22009.0.3","f.22009.0.4","f.22009.0.5",
#"f.22009.0.6","f.22009.0.7","f.22009.0.8","f.22009.0.9","f.22009.0.10",
#"f.22009.0.11","f.22009.0.12","f.22009.0.13","f.22009.0.14","f.22009.0.15",
#"f.22009.0.16","f.22009.0.17","f.22009.0.18","f.22009.0.19","f.22009.0.20",
#"f.22009.0.21","f.22009.0.22","f.22009.0.23","f.22009.0.24","f.22009.0.25",
#"f.22009.0.26","f.22009.0.27","f.22009.0.28","f.22009.0.29","f.22009.0.30",
#"f.22009.0.31","f.22009.0.32","f.22009.0.33","f.22009.0.34","f.22009.0.35",
#"f.22009.0.36","f.22009.0.37","f.22009.0.38","f.22009.0.39","f.22009.0.40")])
#> dim(bcr)
#[1] 407478     63

#colnames(bcr) <- c("FID","rbc","hgb","hct","mcv","mch","mchc","rdw","rtc","plt","mpv","pdw","wbc","lym","mon","neu","eos","bas","hrtc","mrv","msv","sex","age",
#paste("pc",1:40,sep=""))
#bcr.cor <- cor(bcr[,2:21])

#library(gplots)

#tiff(paste("UKBB_BCR_cor.tiff",sep=""),width=1600,height=1600,res=184,compression="lzw")
#heatmap.2(bcr.cor,col=bluered, dendrogram="none", trace="none", Rowv=F, Colv=F, denscol="yellow", cexRow=1.2, cexCol=1.2, symm=T, symkey=T) 
#dev.off()
#####################################################################################

bcr <- phe[,c("f.eid","f.30010.0.0","f.30020.0.0","f.30030.0.0","f.30040.0.0","f.30050.0.0","f.30060.0.0","f.30070.0.0","f.30250.0.0",
"f.30080.0.0","f.30100.0.0","f.30110.0.0","f.30000.0.0","f.30120.0.0","f.30130.0.0",
"f.30140.0.0","f.30150.0.0","f.30160.0.0","f.30300.0.0","f.30260.0.0","f.30270.0.0","f.31.0.0","f.21003.0.0",
"f.22009.0.1","f.22009.0.2","f.22009.0.3","f.22009.0.4","f.22009.0.5",
"f.22009.0.6","f.22009.0.7","f.22009.0.8","f.22009.0.9","f.22009.0.10",
"f.22009.0.11","f.22009.0.12","f.22009.0.13","f.22009.0.14","f.22009.0.15",
"f.22009.0.16","f.22009.0.17","f.22009.0.18","f.22009.0.19","f.22009.0.20",
"f.22009.0.21","f.22009.0.22","f.22009.0.23","f.22009.0.24","f.22009.0.25",
"f.22009.0.26","f.22009.0.27","f.22009.0.28","f.22009.0.29","f.22009.0.30",
"f.22009.0.31","f.22009.0.32","f.22009.0.33","f.22009.0.34","f.22009.0.35",
"f.22009.0.36","f.22009.0.37","f.22009.0.38","f.22009.0.39","f.22009.0.40")]

colnames(bcr) <- c("FID","rbc","hgb","hct","mcv","mch","mchc","rdw","rtc","plt","mpv","pdw","wbc","lym","mon","neu","eos","bas","hrtc","mrv","msv","sex","age",
paste("pc",1:40,sep=""))

#bcr <- bcr[(!is.na(bcr$wbc) & bcr$wbc<=200) | (!is.na(bcr$hgb) & bcr$hgb<=20) | (!is.na(bcr$hct) & bcr$hct<=60) | (!is.na(bcr$plt) & bcr$plt<=1000) | (!is.na(bcr$mpv) & bcr$mpv<=quantile(bcr$mpv,0.96,na.rm=T)),]
bcr <- bcr[(!is.na(bcr$wbc) & bcr$wbc<=200) | (!is.na(bcr$hgb) & bcr$hgb<=20) | (!is.na(bcr$hct) & bcr$hct<=60) | (!is.na(bcr$plt) & bcr$plt<=1000),]
bcr <- bcr[,-c(5:7,19,21)]

#> dim(bcr)
#[1] 418366     58

#### remove those with 0 except for eos and bas that will be replaced with min/10
bcr <- bcr[!is.na(bcr$rtc) & bcr$rtc!=0 & !is.na(bcr$mon) & bcr$mon!=0 & !is.na(bcr$neu) & bcr$neu!=0,]
excl2 <- read.table("/gpfs/gsfs7/users/NHLBI_UKBB/test_MPV/bolt.in_plink_but_not_imputed.FID_IID.652.txt")
bcr <- bcr[!bcr[,1]%in%excl2,]

x.ind <- fread("/gpfs/gsfs7/users/NHLBI_UKBB/imputed/ukb28525_imp_chrX_v3_s486743.sample",header=T,data.table=F)
bcr <- bcr[bcr[,1]%in%x.ind[,1],]
bcr <- na.omit(bcr)

bcr[!is.na(bcr$eos) & bcr$eos==0,"eos"] <- 0.005
bcr[!is.na(bcr$bas) & bcr$bas==0,"bas"] <- 0.005
bcr[!is.na(bcr$rbc),"rbc"] <- bcr[!is.na(bcr$rbc),"rbc"]*1000
bcr[!is.na(bcr$rtc),"rtc"] <- bcr[!is.na(bcr$rtc),"rtc"]*1000

