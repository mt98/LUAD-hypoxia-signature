####This work looks at validating Brian's LUAD signature in multiple datasets, after that look at Sun and Shi signatures and Buffa signature



setwd("/Users/mkhan/Documents/PhD work/Brian RNA-seq scripts")
########################################################
#GSE8894
#HGU133plus2

load("GSE8894_survival.RData")
load("GSE8894_exprs.RData")
###only work on the lung adenocarcinoma
GSE8894_survival<-GSE8894_survival[GSE8894_survival$cell_type=="Adenocarcinoma",]

exprs<-GSE8894_expr[,intersect(colnames(GSE8894_expr),rownames(GSE8894_survival))]###this line is finding those samples which are common between patient samples clinical and expression data
centre_exprs<-sweep(exprs,1,apply(exprs,1,median))####median centering expression data
####The median centering of the 
p<-centre_exprs[original_seed_genes,]
GSE8894_yhat<- pamr.predict(train_model,p, train_threshold)
####Censor the data
GSE8894_survival$status   <-ifelse(GSE8894_survival$RFS_status=="non_recurrence",0,1)
censorship<-60
GSE8894_survival$rfs_month   <- as.numeric( as.character(GSE8894_survival$rfs_month) )
GSE8894_survival_event_censorship  <- ifelse(GSE8894_survival$rfs_month<= censorship &GSE8894_survival$status == 1 , 1 ,0  )
GSE8894_survival_time_censorship   <- ifelse( GSE8894_survival_event_censorship == 0 & GSE8894_survival$rfs_month>=censorship ,censorship ,GSE8894_survival$rfs_month )   
GSE8894_survival    <- cbind( GSE8894_survival, GSE8894_survival_event_censorship, GSE8894_survival_time_censorship )
colnames(GSE8894_survival)[ (ncol(GSE8894_survival)-1):ncol(GSE8894_survival) ] <- c("c_event","c_event_time")
####KM plot
res.cox <- coxph(Surv(c_event_time,c_event)~GSE8894_yhat, data=GSE8894_survival)
summary(res.cox)
fit <- survfit(Surv(c_event_time,c_event)~GSE8894_yhat, data=GSE8894_survival)
par(pty="s")
gplotmeantr<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Relapse free survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantr$plot <- gplotmeantr$plot + labs(
  title    = "GSE8894"        
)
gplotmeantr <- ggpar(
  gplotmeantr,
  font.title    = c(25, "bold"))
gplotmeantr$plot<-gplotmeantr$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantr$table <- ggpar(gplotmeantr$table,
                            font.title = list(size = 16))

###Look at the MVA###use the age, sex, stage 
#scp mqbpkmk3@csf3.itservices.manchester.ac.uk:~/GSE8894_mva.RData "/Users/mkhan/Documents/Brian RNA-seq scripts"
load("GSE8894_mva.RData")
###look at the GSE8894
gse8894<-merge(GSE8894_survival,gse,by="row.names")
res.cox <- coxph(Surv(c_event_time,c_event)~age, data=gse8894)
summary(res.cox)

res.cox <- coxph(Surv(c_event_time,c_event)~gender, data=gse8894)
summary(res.cox)

#GSE3141
#HGU133plus2
load( "GSE3141_survival.RData")
load("GSE3141_exprs.RData")

GSE3141_survival$status   <-ifelse(GSE3141_survival$Status=="alive",0,1)
censorship<-60
GSE3141_survival$Overall_survival_time    <- as.numeric( as.character(GSE3141_survival$Overall_survival_time ) )
GSE3141_survival_event_censorship  <- ifelse(GSE3141_survival$Overall_survival_time <= censorship & GSE3141_survival$status == 1 , 1 ,0  )
GSE3141_survival_time_censorship   <- ifelse( GSE3141_survival_event_censorship == 0 & GSE3141_survival$Overall_survival_time >=censorship ,censorship ,GSE3141_survival$Overall_survival_time )   
GSE3141_survival     <- cbind( GSE3141_survival, GSE3141_survival_event_censorship, GSE3141_survival_time_censorship )
colnames(GSE3141_survival)[ (ncol(GSE3141_survival)-1):ncol(GSE3141_survival) ] <- c("c_event","c_event_time")
GSE3141_survival$Cell.type<-as.character(GSE3141_survival$Cell.type)
GSE3141_survival<-GSE3141_survival[GSE3141_survival$Cell.type==" A",]


exprs<-GSE3141_expr[,intersect(colnames(GSE3141_expr),rownames(GSE3141_survival))]
centre_exprs<-sweep(exprs,1,apply(exprs,1,median))
sig_exprs<-centre_exprs[original_seed_genes,]
GSE3141_yhat<- pamr.predict(train_model,sig_exprs, train_threshold)



res.cox <- coxph(Surv(c_event_time,c_event)~GSE3141_yhat, data=GSE3141_survival)
summary(res.cox)
fit <- survfit(Surv(c_event_time,c_event)~GSE3141_yhat, data=GSE3141_survival)
par(pty="s")
gplotmeantro1<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantro1$plot <- gplotmeantro1$plot + labs(
  title    = "GSE3141"        
)
gplotmeantro1 <- ggpar(
  gplotmeantro1,
  font.title    = c(25, "bold"))
gplotmeantro1$plot<-gplotmeantro1$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantro1$table <- ggpar(gplotmeantro1$table,
                            font.title = list(size = 16))


###Can't find it in the MVA 
rm(list=c('exprs', 'GSE3141_survival',"GSE3141_yhat"))


####################################################
#GSE31210; all Adenocarcinoma, saved as GSE31210_survival object
#HGU133plus2

load("GSE31210_survival.RData")#GSE31210_survival
load("GSE31210_exprs.RData")


exprs<-exprs[,intersect(colnames(exprs),rownames(GSE31210_survival))]
centre_exprs<-sweep(exprs,1,apply(exprs,1,median))

sig_exprs<-madSummary133plus2(hgu133plus2.db,ad_seed_genes,centre_exprs)
rownames(sig_exprs)[which(rownames(sig_exprs)=='CCN5')]<-'WISP2'

GSE31210_yhat<- pamr.predict(train_model,sig_exprs, train_threshold)
res.cox<-coxph(Surv(censored_os_time,censored_os_status)~GSE31210_yhat, data=GSE31210_survival)
summary(res.cox)
survdiff(Surv(censored_os_time,censored_os_status)~GSE31210_yhat, data=GSE31210_survival)

fit <- survfit(Surv(censored_os_time,censored_os_status)~ GSE31210_yhat, data=GSE31210_survival)
par(pty="s")
gplotmeantro7<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantro7$plot <- gplotmeantro7$plot + labs(
  title    = "GSE31210"        
)
gplotmeantro7 <- ggpar(
  gplotmeantro7,
  font.title    = c(25, "bold"))
gplotmeantro7$plot<-gplotmeantro7$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantro7$table <- ggpar(gplotmeantro7$table,
                            font.title = list(size = 16))


GSE31210_fit<-survfit(Surv(censored_os_time,censored_os_status)~GSE31210_yhat, data=GSE31210_survival)
plotKM(GSE31210_fit,GSE31210_survival,"OS")
survdiff(Surv(censored_rfs_time,censored_rfs_status)~GSE31210_yhat, data=GSE31210_survival)
GSE31210_fit<-survfit(Surv(censored_rfs_time,censored_rfs_status)~GSE31210_yhat, data=GSE31210_survival)

fit <- survfit(Surv(censored_rfs_time,censored_rfs_status)~ GSE31210_yhat, data=GSE31210_survival)
par(pty="s")
gplotmeantr4<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Relapse free survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantr4$plot <- gplotmeantr4$plot + labs(
  title    = "GSE31210"        
)
gplotmeantr4 <- ggpar(
  gplotmeantr4,
  font.title    = c(25, "bold"))
gplotmeantr4$plot<-gplotmeantr4$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantr4$table <- ggpar(gplotmeantr4$table,
                            font.title = list(size = 16))


plotKM(GSE31210_fit,GSE31210_survival,"RFS")

rm(list=c('exprs', 'GSE31210_survival',"GSE31210_yhat"))



#############################################################

#GSE42127
#Illumina HWG-6 v3.0
#KM p = 0.09 for whole data set; KM p = 0.05 for patients with adjuvant chemotherapy


load( "GSE42127_LUAD_survival.RData")
load("GSE42127_exprs.RData")


#GSE42127_ad_survival<-GSE42127_ad_survival[which(GSE42127_ad_survival$had_adjuvant_chemo==T),]


exprs<-exprs[,intersect(colnames(exprs),rownames(GSE42127_ad_survival))]
center_exprs<-sweep(exprs,1,apply(exprs,1,median))
sig_exprs<-madSummaryAgilent(illuminaHumanv3.db,ad_seed_genes,center_exprs)

rownames(sig_exprs)[which(rownames(sig_exprs)=='CCN5')]<-'WISP2'


GSE42127_yhat<- pamr.predict(train_model,sig_exprs, train_threshold)
res.cox<-coxph(Surv(censored_os_time,censored_os_status)~GSE42127_yhat, data=GSE42127_ad_survival)
summary(res.cox)
survdiff(Surv(censored_os_time,censored_os_status)~GSE42127_yhat, data=GSE42127_ad_survival)#NS
GSE42127_fit<-survfit(Surv(censored_os_time,censored_os_status)~GSE42127_yhat, data=GSE42127_ad_survival)
plotKM(GSE42127_fit,GSE42127_ad_survival,"OS")

fit <- survfit(Surv(censored_os_time,censored_os_status)~ GSE42127_yhat, data=GSE42127_ad_survival)
par(pty="s")
gplotmeantro8<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantro8$plot <- gplotmeantro8$plot + labs(
  title    = "GSE42127"        
)
gplotmeantro8 <- ggpar(
  gplotmeantro8,
  font.title    = c(25, "bold"))
gplotmeantro8$plot<-gplotmeantro8$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantro8$table <- ggpar(gplotmeantro8$table,
                            font.title = list(size = 16))

rm(list=c('exprs', 'GSE42127_ad_survival',"GSE42127_yhat"))
###########################################################################################

#GSE19188
#HGU133plus2
load("GSE19188_LUAD_survival.RData")
load("GSE19188_exprs.RData")


exprs<-exprs[,intersect(colnames(exprs),rownames(GSE19188_ad_survival))]
centre_exprs<-sweep(exprs,1,apply(exprs,1,median))


sig_exprs<-madSummary133plus2(hgu133plus2.db,ad_seed_genes,centre_exprs)
rownames(sig_exprs)[which(rownames(sig_exprs)=='CCN5')]<-'WISP2'

GSE19188_yhat<- pamr.predict(train_model,sig_exprs, train_threshold)
res.cox<-coxph(Surv(censored_os_time,censored_os_status)~GSE19188_yhat, data=GSE19188_ad_survival)
summary(res.cox)
survdiff(Surv(censored_os_time,censored_os_status)~GSE19188_yhat, data=GSE19188_ad_survival)
GSE19188_fit<-survfit(Surv(censored_os_time,censored_os_status)~GSE19188_yhat, data=GSE19188_ad_survival)

fit <- survfit(Surv(censored_os_time,censored_os_status)~ GSE19188_yhat, data=GSE19188_ad_survival)
par(pty="s")
gplotmeantro9<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantro9$plot <- gplotmeantro9$plot + labs(
  title    = "GSE19188"        
)
gplotmeantro9 <- ggpar(
  gplotmeantro9,
  font.title    = c(25, "bold"))
gplotmeantro9$plot<-gplotmeantro9$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantro9$table <- ggpar(gplotmeantro9$table,
                            font.title = list(size = 16))

plotKM(GSE19188_fit,GSE19188_ad_survival,"OS")


rm(list=c('exprs', 'GSE19188_ad_survival',"GSE19188_yhat"))

##########################################################################################
#GSE72094
#Bespoke array, annotation from pData of GEOquery object

load( "GSE72094_LUAD_survival.RData"    )
load( "GSE72094_exprs.RData"    )

exprs<-exprs[,intersect(colnames(exprs),rownames(GSE72094_ad_survival))]
center_exprs<-sweep(exprs,1,apply(exprs,1,median))
sig_exprs<-center_exprs[intersect(original_seed_genes,rownames(center_exprs)),]


GSE72094_yhat<- pamr.predict(train_model,sig_exprs, train_threshold)
survdiff(Surv(censored_time,censored_status)~GSE72094_yhat, data=GSE72094_ad_survival)
GSE72094_fit<-survfit(Surv(censored_time,censored_status)~GSE72094_yhat, data=GSE72094_ad_survival)
res.cox <- coxph(Surv(censored_time,censored_status)~GSE72094_yhat, data= GSE72094_ad_survival)
summary(res.cox)
fit <- survfit(Surv(censored_time,censored_status)~GSE72094_yhat, data= GSE72094_ad_survival)
par(pty="s")
gplotmeantro2<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantro2$plot <- gplotmeantro2$plot + labs(
  title    = "GSE72094"        
)
gplotmeantro2 <- ggpar(
  gplotmeantro2,
  font.title    = c(25, "bold"))
gplotmeantro2$plot<-gplotmeantro2$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantro2$table <- ggpar(gplotmeantro2$table,
                            font.title = list(size = 16))
####MVA
#scp -r mqbpkmk3@csf3.itservices.manchester.ac.uk:~/GSE72094.rds "/Users/mkhan/Documents/Brian RNA-seq scripts"
GSE72094MVA<-readRDS(file="GSE72094.rds")
GSE72094_MVA<-merge(GSE72094_ad_survival,GSE72094MVA,by="row.names")


GSE72094_MVA$Age<-as.numeric(GSE72094_MVA$Age)
res.cox <- coxph(Surv(censored_time,censored_status)~Age, data=GSE72094_MVA)
summary(res.cox)
res.cox <- coxph(Surv(censored_time,censored_status)~GSE72094_MVA[,14], data=GSE72094_MVA)
summary(res.cox)

stage<-sapply(GSE72094_MVA[,13],function(x){
	if(x=="1"|x=="1A"|x=="1B"|x=="2A"|x=="2B"){y=1}
	else{y=2}})
GSE72094_MVA$stage<-as.factor(stage)
res.cox <- coxph(Surv(censored_time,censored_status)~as.factor(stage), data=GSE72094_MVA)
summary(res.cox)

GSE72094_MVA_edited<-GSE72094_MVA[,12][(!GSE72094_MVA=="Missing"),]
res.cox <- coxph(Surv(censored_time,censored_status)~GSE72094_MVA_edited[,12], data=GSE72094_MVA_edited)
summary(res.cox)
res.cox <- coxph(Surv(censored_time,censored_status)~GSE72094_MVA[,14]+stage+GSE72094_yhat, data=GSE72094_MVA)
summary(res.cox)
rm(list=c('center_exprs','exprs', 'GSE72094_ad_survival',"GSE72094_yhat"))


############################################################################################
#####GSE41271 

#gds <- getGEO("GSE41271",GSEMatrix=TRUE)
#gse41271<-exprs(gds[[1]])
#saveRDS(gse41271, file = "GSE41271.rds")
#gse41271_survival<-(pData(phenoData(gds[[1]]))
#saveRDS(gse41271_survival, file = "GSE41271_survival.rds")

#scp mqbpkmk3@csf3.itservices.manchester.ac.uk:~/GSE41271.rds "/Users/mkhan/Documents/Brian RNA-seq scripts"
#scp mqbpkmk3@csf3.itservices.manchester.ac.uk:~/GSE41271_survival.rds "/Users/mkhan/Documents/Brian RNA-seq scripts"

GSE41271<-readRDS(file = "GSE41271.rds")



###processing the survival data for GSE42127
#write.table(GSE41271_survival,file="GSE41271_survival.csv",sep=",",row.names=FALSE)
survivalGSE41271<-read.csv(file="GSE41271_survival.csv",sep=",",header=T)
####censoring data
survivalGSE41271$status<-ifelse(survivalGSE41271$vital.statistics.ch1=="A",0,1)
censorship<-60
survivalGSE41271$OS   <- as.numeric( as.character(survivalGSE41271$OS ) )
GSE41271_survival_event_censorship  <- ifelse(survivalGSE41271$OS <= censorship &survivalGSE41271$status == 1 , 1 ,0  )
GSE41271_survival_time_censorship   <- ifelse( GSE41271_survival_event_censorship == 0 & survivalGSE41271$OS>=censorship ,censorship ,
survivalGSE41271$OS )   
survivalGSE41271   <- cbind(survivalGSE41271, GSE41271_survival_event_censorship, GSE41271_survival_time_censorship )
colnames(survivalGSE41271)[ (ncol(survivalGSE41271)-1):ncol(survivalGSE41271) ] <- c("c_event","c_event_time")

survivalGSE41271$rstatus<-ifelse(survivalGSE41271$recurrence.ch1=="N",0,1)
censorship<-60
survivalGSE41271$RFS   <- as.numeric( as.character(survivalGSE41271$RFS ) )
GSE41271_rsurvival_event_censorship  <- ifelse(survivalGSE41271$RFS <= censorship &survivalGSE41271$rstatus == 1 , 1 ,0  )
GSE41271_rsurvival_time_censorship   <- ifelse( GSE41271_rsurvival_event_censorship == 0 & survivalGSE41271$RFS>=censorship ,censorship ,
survivalGSE41271$RFS )   
survivalGSE41271   <- cbind(survivalGSE41271, GSE41271_rsurvival_event_censorship, GSE41271_rsurvival_time_censorship )
colnames(survivalGSE41271)[ (ncol(survivalGSE41271)-1):ncol(survivalGSE41271) ] <- c("c_revent","c_revent_time")
GSE41271_survival<-survivalGSE41271[survivalGSE41271$histology.ch1=="Adenocarcinoma",]
rownames(GSE41271_survival)<-GSE41271_survival$geo_accession
exprs<-GSE41271[,intersect(colnames(GSE41271),rownames(GSE41271_survival))]
exprs<-log2(exprs)
centre_exprs<-sweep(exprs,1,apply(exprs,1,median))
sig_exprs<-madSummaryAgilent(illuminaHumanv3.db,ad_seed_genes,centre_exprs)
rownames(sig_exprs)[which(rownames(sig_exprs)=='CCN5')]<-'WISP2'
GSE42127_yhat<- pamr.predict(train_model,sig_exprs, train_threshold)

res.cox <- coxph(Surv(c_event_time,c_event)~GSE42127_yhat, data=survivalGSE41271 )
summary(res.cox)
fit <- survfit(Surv(c_event_time,c_event)~GSE42127_yhat, data=survivalGSE41271 )
par(pty="s")
gplotmeantro3<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantro3$plot <- gplotmeantro3$plot + labs(
  title    = "GSE41271"        
)
gplotmeantro3 <- ggpar(
  gplotmeantro3,
  font.title    = c(25, "bold"))
gplotmeantro3$plot<-gplotmeantro3$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantro3$table <- ggpar(gplotmeantro3$table,
                            font.title = list(size = 16))

####RFS
res.cox <- coxph(Surv(c_revent_time,c_revent)~GSE42127_yhat, data=survivalGSE41271 )
summary(res.cox)
fit <- survfit(Surv(c_revent_time,c_revent)~GSE42127_yhat, data=survivalGSE41271 )
par(pty="s")
gplotmeantr1<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Relapse free survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantr1$plot <- gplotmeantr1$plot + labs(
  title    = "GSE41271"        
)
gplotmeantr1 <- ggpar(
  gplotmeantr1,
  font.title    = c(25, "bold"))
gplotmeantr1$plot<-gplotmeantr1$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantr1$table <- ggpar(gplotmeantr1$table,
                            font.title = list(size = 16))

#####
#gds <- getGEO("GSE29013",GSEMatrix=TRUE)
#gse29013<-exprs(gds[[1]])
#saveRDS(gse29013, file = "GSE29013.rds")
#gse29013_survival<-(pData(phenoData(gds[[1]])))
#saveRDS(gse29013_survival, file = "GSE29013_survival.rds")

#scp mqbpkmk3@csf3.itservices.manchester.ac.uk:~/GSE29013.rds "/Users/mkhan/Documents/Brian RNA-seq scripts"
#scp mqbpkmk3@csf3.itservices.manchester.ac.uk:~/GSE29013_survival.rds "/Users/mkhan/Documents/Brian RNA-seq scripts"

GSE29013<-readRDS(file = "GSE29013.rds")
GSE29013_survival<-readRDS(file = "GSE29013_survival.rds")
GSE29013_survival$time_to_progression<-as.numeric(GSE29013_survival$time_to_progression) * 12
GSE29013_survival$death_time<-as.numeric(GSE29013_survival$death_time) * 12

####censoring data

censorship<-60
GSE29013_survival_event_censorship  <- ifelse(GSE29013_survival$death_time<= censorship &GSE29013_survival$death_event == 1 , 1 ,0  )
GSE29013_survival_time_censorship   <- ifelse(GSE29013_survival_event_censorship == 0 & GSE29013_survival$death_time>=censorship ,censorship ,
GSE29013_survival$death_time )   
GSE29013_survival  <- cbind(GSE29013_survival, GSE29013_survival_event_censorship, GSE29013_survival_time_censorship )
colnames(GSE29013_survival)[ (ncol(GSE29013_survival)-1):ncol(GSE29013_survival) ] <- c("c_event","c_event_time")

GSE29013_rsurvival_event_censorship  <- ifelse(GSE29013_survival$time_to_progression<= censorship &GSE29013_survival$progression == 1 , 1 ,0  )
GSE29013_rsurvival_time_censorship   <- ifelse(GSE29013_rsurvival_event_censorship == 0 & GSE29013_survival$time_to_progression>=censorship ,censorship ,
GSE29013_survival$time_to_progression)   
GSE29013_survival  <- cbind(GSE29013_survival, GSE29013_rsurvival_event_censorship, GSE29013_rsurvival_time_censorship )
colnames(GSE29013_survival)[ (ncol(GSE29013_survival)-1):ncol(GSE29013_survival) ] <- c("c_revent","c_revent_time")


GSE29013_survival<-GSE29013_survival[GSE29013_survival$histology=="Adenocarcinoma",]

exprs<-GSE29013[,intersect(colnames(GSE29013),rownames(GSE29013_survival))]
centre_exprs<-sweep(exprs,1,apply(exprs,1,median))
sig_exprs<-madSummary133plus2(hgu133plus2.db,ad_seed_genes,centre_exprs)
rownames(sig_exprs)[which(rownames(sig_exprs)=='CCN5')]<-'WISP2'
GSE29013_yhat<- pamr.predict(train_model,sig_exprs, train_threshold)

res.cox <- coxph(Surv(c_event_time,c_event)~GSE29013_yhat, data=GSE29013_survival )
summary(res.cox)
fit <- survfit(Surv(c_event_time,c_event)~GSE29013_yhat, data=GSE29013_survival )
par(pty="s")
gplotmeantro4<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantro4$plot <- gplotmeantro4$plot + labs(
  title    = "GSE29013"        
)
gplotmeantro4 <- ggpar(
  gplotmeantro4,
  font.title    = c(25, "bold"))
gplotmeantro4$plot<-gplotmeantro4$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantro4$table <- ggpar(gplotmeantro4$table,
                            font.title = list(size = 16))

####RFS
res.cox <- coxph(Surv(c_revent_time,c_revent)~GSE29013_yhat, data=GSE29013_survival )
summary(res.cox)
fit <- survfit(Surv(c_revent_time,c_revent)~GSE29013_yhat, data=GSE29013_survival )
par(pty="s")
gplotmeantr2<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Progression free survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantr2$plot <- gplotmeantr2$plot + labs(
  title    = "GSE29013"        
)
gplotmeantr2 <- ggpar(
  gplotmeantr2,
  font.title    = c(25, "bold"))
gplotmeantr2$plot<-gplotmeantr2$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantr2$table <- ggpar(gplotmeantr2$table,
                            font.title = list(size = 16))

#####


#gds <- getGEO("GSE50081",GSEMatrix=TRUE)
#gse50081<-exprs(gds[[1]])
#saveRDS(gse50081, file = "GSE50081.rds")
#gse50081_survival<-(pData(phenoData(gds[[1]])))
#saveRDS(gse50081_survival, file = "GSE50081_survival.rds")

#scp mqbpkmk3@csf3.itservices.manchester.ac.uk:~/GSE50081.rds "/Users/mkhan/Documents/Brian RNA-seq scripts"
#scp mqbpkmk3@csf3.itservices.manchester.ac.uk:~/GSE50081_survival.rds "/Users/mkhan/Documents/Brian RNA-seq scripts"
GSE50081<-readRDS(file = "GSE50081.rds")
GSE50081_survival<-readRDS(file = "GSE50081_survival.rds")

GSE50081_survival[,54]<-as.numeric(as.character(GSE50081_survival[,54]))*12
GSE50081_survival$Status<-ifelse(GSE50081_survival[,53]=="alive",0,1)
censorship<-60
GSE50081_survival_event_censorship  <- ifelse(GSE50081_survival[,54]<= censorship &GSE50081_survival$Status == 1 , 1 ,0  )
GSE50081_survival_time_censorship   <- ifelse(GSE50081_survival_event_censorship == 0 & GSE50081_survival[,54]>=censorship ,censorship ,
GSE50081_survival[,54])   
GSE50081_survival  <- cbind(GSE50081_survival, GSE50081_survival_event_censorship, GSE50081_survival_time_censorship )
colnames(GSE50081_survival)[ (ncol(GSE50081_survival)-1):ncol(GSE50081_survival) ] <- c("c_event","c_event_time")

GSE50081_survival<-GSE50081_survival[GSE50081_survival$histology=="adenocarcinoma",]
exprs<-GSE50081[,intersect(colnames(GSE50081),rownames(GSE50081_survival))]
centre_exprs<-sweep(exprs,1,apply(exprs,1,median))
sig_exprs<-madSummary133plus2(hgu133plus2.db,ad_seed_genes,centre_exprs)
rownames(sig_exprs)[which(rownames(sig_exprs)=='CCN5')]<-'WISP2'
GSE50081_yhat<- pamr.predict(train_model,sig_exprs, train_threshold)

res.cox <- coxph(Surv(c_event_time,c_event)~GSE50081_yhat, data=GSE50081_survival )
summary(res.cox)
fit <- survfit(Surv(c_event_time,c_event)~GSE50081_yhat, data=GSE50081_survival )
par(pty="s")
gplotmeantro5<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantro5$plot <- gplotmeantro5$plot + labs(
  title    = "GSE50081"        
)
gplotmeantro5 <- ggpar(
  gplotmeantro5,
  font.title    = c(25, "bold"))
gplotmeantro5$plot<-gplotmeantro5$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantro5$table <- ggpar(gplotmeantro5$table,
                            font.title = list(size = 16))



####MVA
res.cox <- coxph(Surv(c_event_time,c_event)~as.numeric(GSE50081_survival[,44]), data=GSE50081_survival )
summary(res.cox)

res.cox <- coxph(Surv(c_event_time,c_event)~GSE50081_survival[,50], data=GSE50081_survival )
summary(res.cox)

GSE50081_survival_edited<-GSE50081_survival[!GSE50081_survival[,51]=="Unable to determine",]
GSE50081_survival_edited$smoking<-ifelse(GSE50081_survival_edited$smoking=="Never",0,1)
res.cox <- coxph(Surv(c_event_time,c_event)~smoking, data=GSE50081_survival_edited )
summary(res.cox)

GSE50081_survival$Stage<-sapply(GSE50081_survival[,52],function(x){
	if(x=="1A"|x=="1B"){y=1} else{y=2}})


res.cox <- coxph(Surv(c_event_time,c_event)~Stage, data=GSE50081_survival )
summary(res.cox)
res.cox <- coxph(Surv(c_event_time,c_event)~Stage+GSE50081_yhat, data=GSE50081_survival )
summary(res.cox)



#gds <- getGEO("GSE37745",GSEMatrix=TRUE)
#gse37745<-exprs(gds[[1]])
#saveRDS(gse37745, file = "GSE37745.rds")
#gse37745_survival<-(pData(phenoData(gds[[1]])))
#saveRDS(gse37745_survival, file = "GSE37745_survival.rds")
#scp mqbpkmk3@csf3.itservices.manchester.ac.uk:~/GSE37745.rds "/Users/mkhan/Documents/Brian RNA-seq scripts"
#scp mqbpkmk3@csf3.itservices.manchester.ac.uk:~/GSE37745_survival.rds "/Users/mkhan/Documents/Brian RNA-seq scripts"
GSE37745<-readRDS(file = "GSE37745.rds")
GSE37745_survival<-readRDS(file = "GSE37745_survival.rds")
GSE37745_survival[,43]<-as.numeric(GSE37745_survival[,43]) / 12
GSE37745_survival$Status<-ifelse(GSE37745_survival[,45]=="no",0,1)
####censoring data

censorship<-60
GSE37745_survival_event_censorship  <- ifelse(GSE37745_survival[,43]<= censorship &GSE37745_survival$Status == 1 , 1 ,0  )
GSE37745_survival_time_censorship   <- ifelse(GSE37745_survival_event_censorship == 0 & GSE37745_survival[,43]>=censorship ,censorship ,
GSE37745_survival[,43] )   
GSE37745_survival  <- cbind(GSE37745_survival, GSE37745_survival_event_censorship, GSE37745_survival_time_censorship )
colnames(GSE37745_survival)[ (ncol(GSE37745_survival)-1):ncol(GSE37745_survival) ] <- c("c_event","c_event_time")
GSE37745_survival<-GSE37745_survival[GSE37745_survival$histology=="adeno",]

exprs<-GSE37745[,intersect(colnames(GSE37745),rownames(GSE37745_survival))]
centre_exprs<-sweep(exprs,1,apply(exprs,1,median))
sig_exprs<-madSummary133plus2(hgu133plus2.db,ad_seed_genes,centre_exprs)
rownames(sig_exprs)[which(rownames(sig_exprs)=='CCN5')]<-'WISP2'
GSE37745_yhat<- pamr.predict(train_model,sig_exprs, train_threshold)

res.cox <- coxph(Surv(c_event_time,c_event)~GSE37745_yhat, data=GSE37745_survival )
summary(res.cox)###confusing result 0.637 [0.30-1.35]



#####

#gds <- getGEO("GSE30219",GSEMatrix=TRUE)
#gse30219<-exprs(gds[[1]])
#saveRDS(gse30219, file = "GSE30219.rds")
#gse30219_survival<-(pData(phenoData(gds[[1]])))
#saveRDS(gse30219_survival, file = "GSE30219_survival.rds")
#scp mqbpkmk3@csf3.itservices.manchester.ac.uk:~/GSE30219.rds "/Users/mkhan/Documents/Brian RNA-seq scripts"
#scp mqbpkmk3@csf3.itservices.manchester.ac.uk:~/GSE30219_survival.rds "/Users/mkhan/Documents/Brian RNA-seq scripts"

#scp mqbpkmk3@csf3.itservices.manchester.ac.uk:~/GSE37745_survival.rds "/Users/mkhan/Documents/Brian RNA-seq scripts"
GSE30219<-readRDS(file = "GSE30219.rds")
GSE30219_survival<-readRDS(file = "GSE30219_survival.rds")
GSE30219_survival[,43]<-as.numeric(GSE30219_survival[,43]) 
GSE30219_survival<-GSE30219_survival[!is.na(GSE30219_survival[,43]),]
GSE30219_survival$Status<-ifelse(GSE30219_survival[,50]=="ALIVE",0,1)
####censoring data

censorship<-60
GSE30219_survival_event_censorship  <- ifelse(GSE30219_survival[,43]<= censorship &GSE30219_survival$Status == 1 , 1 ,0  )
GSE30219_survival_time_censorship   <- ifelse(GSE30219_survival_event_censorship == 0 & GSE30219_survival[,43]>=censorship ,censorship ,
GSE30219_survival[,43] )   
GSE30219_survival  <- cbind(GSE30219_survival, GSE30219_survival_event_censorship, GSE30219_survival_time_censorship )
colnames(GSE30219_survival)[ (ncol(GSE30219_survival)-1):ncol(GSE30219_survival) ] <- c("c_event","c_event_time")
GSE30219_survival<-GSE30219_survival[GSE30219_survival[,45]=="ADC",]

exprs<-GSE30219[,intersect(colnames(GSE30219),rownames(GSE30219_survival))]
centre_exprs<-sweep(exprs,1,apply(exprs,1,median))
sig_exprs<-madSummary133plus2(hgu133plus2.db,ad_seed_genes,centre_exprs)
rownames(sig_exprs)[which(rownames(sig_exprs)=='CCN5')]<-'WISP2'
GSE30219_yhat<- pamr.predict(train_model,sig_exprs, train_threshold)

res.cox <- coxph(Surv(c_event_time,c_event)~GSE30219_yhat, data=GSE30219_survival )
summary(res.cox)
fit <- survfit(Surv(c_event_time,c_event)~ GSE30219_yhat, data= GSE30219_survival  )
par(pty="s")
gplotmeantro6<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantro6$plot <- gplotmeantro6$plot + labs(
  title    = "GSE30219"        
)
gplotmeantro6 <- ggpar(
  gplotmeantro6,
  font.title    = c(25, "bold"))
gplotmeantro6$plot<-gplotmeantro6$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantro6$table <- ggpar(gplotmeantro6$table,
                            font.title = list(size = 16))




splots<-list(gplotmeantro1,gplotmeantro2,gplotmeantro3,gplotmeantro4,gplotmeantro5,gplotmeantro6,gplotmeantro7,gplotmeantro8,gplotmeantro9)
s<-arrange_ggsurvplots(splots,ncol=3,nrow=3)
ggsave("/Users/mkhan/Documents/Brian RNA-seq scripts/Figure 3.pdf",plot=s,height=18,width=18)

splots<-list(gplotmeantr,gplotmeantr1, gplotmeantr4)
s<-arrange_ggsurvplots(splots,ncol=2,nrow=2)
ggsave("/Users/mkhan/Documents/Brian RNA-seq scripts/Figure S_RFS.pdf",plot=s,height=12,width=12)

####GSE68465###did not get all the probesets
#gds <- getGEO("GSE68465",GSEMatrix=TRUE)
#gse68465<-exprs(gds[[1]])
#saveRDS(gse68465, file = "GSE68465.rds")
#gse68465_survival<-(pData(phenoData(gds[[1]])))
#saveRDS(gse68465_survival, file = "GSE68465_survival.rds")
#scp mqbpkmk3@csf3.itservices.manchester.ac.uk:~/GSE68465.rds "/Users/mkhan/Documents/Brian RNA-seq scripts"
#scp mqbpkmk3@csf3.itservices.manchester.ac.uk:~/GSE68465_survival.rds "/Users/mkhan/Documents/Brian RNA-seq scripts"
GSE68465<-readRDS(file = "GSE68465.rds")
GSE68465_survival<-readRDS(file = "GSE68465_survival.rds")
GSE68465_survival<-GSE68465_survival[GSE68465_survival$disease_state=="Lung Adenocarcinoma",]
###Survival

GSE68465_survival$Status<-ifelse(GSE68465_survival$vital_status=="Alive",0,1)
GSE68465_survival$months_to_last_contact_or_death<-as.numeric(GSE68465_survival$months_to_last_contact_or_death)
censorship<-60
GSE68465_survival_event_censorship  <- ifelse(GSE68465_survival$months_to_last_contact_or_death<= censorship &GSE68465_survival$Status == 1 , 1 ,0  )
GSE68465_survival_time_censorship   <- ifelse(GSE68465_survival_event_censorship == 0 & GSE68465_survival$months_to_last_contact_or_death>=censorship ,censorship ,
GSE68465_survival$months_to_last_contact_or_death )   
GSE68465_survival  <- cbind(GSE68465_survival,GSE68465_survival_event_censorship,GSE68465_survival_time_censorship )
colnames(GSE68465_survival)[ (ncol(GSE68465_survival)-1):ncol(GSE68465_survival) ] <- c("c_event","c_event_time")

exprs<-GSE68465[,intersect(colnames(GSE68465),rownames(GSE68465_survival))]
exprs<-log2(exprs)
centre_exprs<-sweep(exprs,1,apply(exprs,1,median))
sig_exprs<-madSummary133plus2(hgu133a2.db,ad_seed_genes,centre_exprs)###problem with annotations


####GSE14814###did not get matching probesets
#gds <- getGEO("GSE14814",GSEMatrix=TRUE)
#gse14814<-exprs(gds[[1]])
#saveRDS(gse14814, file = "GSE14814.rds")
#gse14814_survival<-(pData(phenoData(gds[[1]])))
#saveRDS(gse14814_survival, file = "GSE14814_survival.rds")
#scp mqbpkmk3@csf3.itservices.manchester.ac.uk:~/GSE14814.rds "/Users/mkhan/Documents/Brian RNA-seq scripts"
#scp mqbpkmk3@csf3.itservices.manchester.ac.uk:~/GSE14814_survival.rds "/Users/mkhan/Documents/Brian RNA-seq scripts"
GSE14814<-readRDS(file = "GSE14814.rds")
GSE14814_survival<-readRDS(file = "GSE14814_survival.rds")
GSE14814_survival<-GSE14814_survival[GSE14814_survival[,53]=="ADC",]
###Survival

GSE14814_survival$Status<-ifelse(GSE14814_survival$vital_status=="Alive",0,1)
GSE14814_survival$months_to_last_contact_or_death<-as.numeric(GSE14814_survival$months_to_last_contact_or_death)
censorship<-60
GSE14814_survival_event_censorship  <- ifelse(GSE14814_survival$months_to_last_contact_or_death<= censorship &GSE14814_survival$Status == 1 , 1 ,0  )
GSE14814_survival_time_censorship   <- ifelse(GSE14814_survival_event_censorship == 0 & GSE14814_survival$months_to_last_contact_or_death>=censorship ,censorship ,
GSE14814_survival$months_to_last_contact_or_death )   
GSE14814_survival  <- cbind(GSE14814_survival,GSE14814_survival_event_censorship,GSE14814_survival_time_censorship )
colnames(GSE14814_survival)[ (ncol(GSE14814_survival)-1):ncol(GSE14814_survival) ] <- c("c_event","c_event_time")

exprs<-GSE14814[,intersect(colnames(GSE14814),rownames(GSE14814_survival))]
centre_exprs<-sweep(exprs,1,apply(exprs,1,median))
sig_exprs<-madSummary133plus2(hgu133a2.db,ad_seed_genes[1:10],centre_exprs)###problem with annotations



###Look at the different LUAD signatures performance in the different cohorts individually

###Sun signature
LUAD_genes<-read.table("/Users/mkhan/Documents/PhD work/Brian RNA-seq scripts/luad.csv",sep=",",header=T)
rownames(LUAD_genes)<-LUAD_genes$Gene.symbols
load("GSE8894_survival.RData")
load("GSE8894_exprs.RData")
###only work on the lung adenocarcinoma
GSE8894_survival<-GSE8894_survival[GSE8894_survival$cell_type=="Adenocarcinoma",]

exprs<-GSE8894_expr[,intersect(colnames(GSE8894_expr),rownames(GSE8894_survival))]###this line is finding those samples which are common between patient samples clinical and 
####Censor the data
GSE8894_survival$status   <-ifelse(GSE8894_survival$RFS_status=="non_recurrence",0,1)
censorship<-60
GSE8894_survival$rfs_month   <- as.numeric( as.character(GSE8894_survival$rfs_month) )
GSE8894_survival_event_censorship  <- ifelse(GSE8894_survival$rfs_month<= censorship &GSE8894_survival$status == 1 , 1 ,0  )
GSE8894_survival_time_censorship   <- ifelse( GSE8894_survival_event_censorship == 0 & GSE8894_survival$rfs_month>=censorship ,censorship ,GSE8894_survival$rfs_month )   
GSE8894_survival    <- cbind( GSE8894_survival, GSE8894_survival_event_censorship, GSE8894_survival_time_censorship )
colnames(GSE8894_survival)[ (ncol(GSE8894_survival)-1):ncol(GSE8894_survival) ] <- c("c_event","c_event_time")
####KM plot
rownames(LUAD_genes)[11]<-"BCO1"
DataLUADsig<-exprs[rownames(LUAD_genes),]
score<-sapply(colnames(DataLUADsig),function(x){
	y<-LUAD_genes[1,2]*DataLUADsig[rownames(LUAD_genes)[1],x]+
	LUAD_genes[2,2]*DataLUADsig[rownames(LUAD_genes)[2],x]+
	LUAD_genes[3,2]*DataLUADsig[rownames(LUAD_genes)[3],x]
+LUAD_genes[4,2]*DataLUADsig[rownames(LUAD_genes)[4],x]+LUAD_genes[5,2]*DataLUADsig[rownames(LUAD_genes)[5],x]+LUAD_genes[6,2]*DataLUADsig[rownames(LUAD_genes)[6],x]+LUAD_genes[7,2]*DataLUADsig[rownames(LUAD_genes)[7],x]+LUAD_genes[8,2]*DataLUADsig[rownames(LUAD_genes)[8],x]+LUAD_genes[9,2]*DataLUADsig[rownames(LUAD_genes)[9],x]+LUAD_genes[10,2]*DataLUADsig[rownames(LUAD_genes)[10],x]+LUAD_genes[11,2]*DataLUADsig[rownames(LUAD_genes)[11],x]+LUAD_genes[12,2]*DataLUADsig[rownames(LUAD_genes)[12],x]+LUAD_genes[13,2]*DataLUADsig[rownames(LUAD_genes)[13],x]+LUAD_genes[14,2]*DataLUADsig[rownames(LUAD_genes)[14],x]+LUAD_genes[15,2]*DataLUADsig[rownames(LUAD_genes)[15],x]+LUAD_genes[16,2]*DataLUADsig[rownames(LUAD_genes)[16],x]})
median_stratification<-ifelse(score>quantile(score,0.50),"hypoxia","normoxia")
luadsig_median_stratification<-as.data.frame(median_stratification)
####median stratification
survival_luad<-merge(GSE8894_survival,luadsig_median_stratification,by="row.names")
colnames(survival_luad)[ncol(survival_luad)]<-"median_stratification"
survival_luad$median_stratification<-factor(survival_luad$median_stratification,levels=c("normoxia","hypoxia"))
res.cox <- coxph(Surv(c_event_time,c_event) ~ median_stratification, data=survival_luad)
summary(res.cox)#1.83[0.91-3.69] p=0.09

###GSE3141
#HGU133plus2
load( "GSE3141_survival.RData")
load("GSE3141_exprs.RData")

GSE3141_survival$status   <-ifelse(GSE3141_survival$Status=="alive",0,1)
censorship<-60
GSE3141_survival$Overall_survival_time    <- as.numeric( as.character(GSE3141_survival$Overall_survival_time ) )
GSE3141_survival_event_censorship  <- ifelse(GSE3141_survival$Overall_survival_time <= censorship & GSE3141_survival$status == 1 , 1 ,0  )
GSE3141_survival_time_censorship   <- ifelse( GSE3141_survival_event_censorship == 0 & GSE3141_survival$Overall_survival_time >=censorship ,censorship ,GSE3141_survival$Overall_survival_time )   
GSE3141_survival     <- cbind( GSE3141_survival, GSE3141_survival_event_censorship, GSE3141_survival_time_censorship )
colnames(GSE3141_survival)[ (ncol(GSE3141_survival)-1):ncol(GSE3141_survival) ] <- c("c_event","c_event_time")
GSE3141_survival$Cell.type<-as.character(GSE3141_survival$Cell.type)
GSE3141_survival<-GSE3141_survival[GSE3141_survival$Cell.type==" A",]
exprs<-GSE3141_expr[,intersect(colnames(GSE3141_expr),rownames(GSE3141_survival))]
DataLUADsig<-exprs[rownames(LUAD_genes),]
score<-sapply(colnames(DataLUADsig),function(x){
	y<-LUAD_genes[1,2]*DataLUADsig[rownames(LUAD_genes)[1],x]+
	LUAD_genes[2,2]*DataLUADsig[rownames(LUAD_genes)[2],x]+
	LUAD_genes[3,2]*DataLUADsig[rownames(LUAD_genes)[3],x]
+LUAD_genes[4,2]*DataLUADsig[rownames(LUAD_genes)[4],x]+LUAD_genes[5,2]*DataLUADsig[rownames(LUAD_genes)[5],x]+LUAD_genes[6,2]*DataLUADsig[rownames(LUAD_genes)[6],x]+LUAD_genes[7,2]*DataLUADsig[rownames(LUAD_genes)[7],x]+LUAD_genes[8,2]*DataLUADsig[rownames(LUAD_genes)[8],x]+LUAD_genes[9,2]*DataLUADsig[rownames(LUAD_genes)[9],x]+LUAD_genes[10,2]*DataLUADsig[rownames(LUAD_genes)[10],x]+LUAD_genes[11,2]*DataLUADsig[rownames(LUAD_genes)[11],x]+LUAD_genes[12,2]*DataLUADsig[rownames(LUAD_genes)[12],x]+LUAD_genes[13,2]*DataLUADsig[rownames(LUAD_genes)[13],x]+LUAD_genes[14,2]*DataLUADsig[rownames(LUAD_genes)[14],x]+LUAD_genes[15,2]*DataLUADsig[rownames(LUAD_genes)[15],x]+LUAD_genes[16,2]*DataLUADsig[rownames(LUAD_genes)[16],x]})
median_stratification<-ifelse(score>quantile(score,0.50),"hypoxia","normoxia")
luadsig_median_stratification<-as.data.frame(median_stratification)
####median stratification
survival_luad<-merge(GSE3141_survival,luadsig_median_stratification,by="row.names")
colnames(survival_luad)[ncol(survival_luad)]<-"median_stratification"
survival_luad$median_stratification<-factor(survival_luad$median_stratification,levels=c("normoxia","hypoxia"))
res.cox <- coxph(Surv(c_event_time,c_event) ~ median_stratification, data=survival_luad)
summary(res.cox)
fit <- survfit(Surv(c_event_time,c_event)~ median_stratification, data=survival_luad)
par(pty="s")
gplotmeantro1<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantro1$plot <- gplotmeantro1$plot + labs(
  title    = "GSE3141 (Sun signature)"        
)
gplotmeantro1 <- ggpar(
  gplotmeantro1,
  font.title    = c(25, "bold"))
gplotmeantro1$plot<-gplotmeantro1$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantro1$table <- ggpar(gplotmeantro1$table,
                            font.title = list(size = 16))


####Next dataset

load("GSE31210_survival.RData")#GSE31210_survival
load("GSE31210_exprs.RData")


exprs<-exprs[,intersect(colnames(exprs),rownames(GSE31210_survival))]

seed_genes<-rownames(LUAD_genes)
sig_exprs<-madSummary133plus2_Sun(hgu133plus2.db,seed_genes,exprs)

DataLUADsig<-sig_exprs[rownames(LUAD_genes),]
score<-sapply(colnames(DataLUADsig),function(x){
	y<-LUAD_genes[1,2]*DataLUADsig[rownames(LUAD_genes)[1],x]+
	LUAD_genes[2,2]*DataLUADsig[rownames(LUAD_genes)[2],x]+
	LUAD_genes[3,2]*DataLUADsig[rownames(LUAD_genes)[3],x]
+LUAD_genes[4,2]*DataLUADsig[rownames(LUAD_genes)[4],x]+LUAD_genes[5,2]*DataLUADsig[rownames(LUAD_genes)[5],x]+LUAD_genes[6,2]*DataLUADsig[rownames(LUAD_genes)[6],x]+LUAD_genes[7,2]*DataLUADsig[rownames(LUAD_genes)[7],x]+LUAD_genes[8,2]*DataLUADsig[rownames(LUAD_genes)[8],x]+LUAD_genes[9,2]*DataLUADsig[rownames(LUAD_genes)[9],x]+LUAD_genes[10,2]*DataLUADsig[rownames(LUAD_genes)[10],x]+LUAD_genes[11,2]*DataLUADsig[rownames(LUAD_genes)[11],x]+LUAD_genes[12,2]*DataLUADsig[rownames(LUAD_genes)[12],x]+LUAD_genes[13,2]*DataLUADsig[rownames(LUAD_genes)[13],x]+LUAD_genes[14,2]*DataLUADsig[rownames(LUAD_genes)[14],x]+LUAD_genes[15,2]*DataLUADsig[rownames(LUAD_genes)[15],x]+LUAD_genes[16,2]*DataLUADsig[rownames(LUAD_genes)[16],x]})
median_stratification<-ifelse(score>quantile(score,0.50),"hypoxia","normoxia")
luadsig_median_stratification<-as.data.frame(median_stratification)
####median stratification
survival_luad<-merge(GSE31210_survival,luadsig_median_stratification,by="row.names")
colnames(survival_luad)[ncol(survival_luad)]<-"median_stratification"
survival_luad$median_stratification<-factor(survival_luad$median_stratification,levels=c("normoxia","hypoxia"))


GSE31210_fit<-survfit(Surv(censored_os_time,censored_os_status)~median_stratification, data=survival_luad)
fit <- survfit(Surv(censored_os_time,censored_os_status)~median_stratification, data=survival_luad)
par(pty="s")
gplotmeantro7<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantro7$plot <- gplotmeantro7$plot + labs(
  title    = "GSE31210 (Sun signature)"        
)
gplotmeantro7 <- ggpar(
  gplotmeantro7,
  font.title    = c(25, "bold"))
gplotmeantro7$plot<-gplotmeantro7$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantro7$table <- ggpar(gplotmeantro7$table,
                            font.title = list(size = 16))

###MVA

res.cox<-coxph(Surv(censored_os_time,censored_os_status)~age, data=survival_luad)
summary(res.cox)

res.cox<-coxph(Surv(censored_os_time,censored_os_status)~gender, data=survival_luad)
summary(res.cox)

res.cox<-coxph(Surv(censored_os_time,censored_os_status)~smoking_status, data=survival_luad)
summary(res.cox)

survival_luad$Stage<-ifelse(survival_luad$pathological_stage=="II",2,1)

res.cox<-coxph(Surv(censored_os_time,censored_os_status)~Stage, data=survival_luad)
summary(res.cox)


res.cox<-coxph(Surv(censored_os_time,censored_os_status)~Stage+gender+smoking_status+median_stratification, data=survival_luad)
summary(res.cox)



###

load("GSE19188_LUAD_survival.RData")
load("GSE19188_exprs.RData")


exprs<-exprs[,intersect(colnames(exprs),rownames(GSE19188_ad_survival))]
sig_exprs<-madSummary133plus2_Sun(hgu133plus2.db,seed_genes,exprs)
DataLUADsig<-sig_exprs[rownames(LUAD_genes),]
score<-sapply(colnames(DataLUADsig),function(x){
	y<-LUAD_genes[1,2]*DataLUADsig[rownames(LUAD_genes)[1],x]+
	LUAD_genes[2,2]*DataLUADsig[rownames(LUAD_genes)[2],x]+
	LUAD_genes[3,2]*DataLUADsig[rownames(LUAD_genes)[3],x]
+LUAD_genes[4,2]*DataLUADsig[rownames(LUAD_genes)[4],x]+LUAD_genes[5,2]*DataLUADsig[rownames(LUAD_genes)[5],x]+LUAD_genes[6,2]*DataLUADsig[rownames(LUAD_genes)[6],x]+LUAD_genes[7,2]*DataLUADsig[rownames(LUAD_genes)[7],x]+LUAD_genes[8,2]*DataLUADsig[rownames(LUAD_genes)[8],x]+LUAD_genes[9,2]*DataLUADsig[rownames(LUAD_genes)[9],x]+LUAD_genes[10,2]*DataLUADsig[rownames(LUAD_genes)[10],x]+LUAD_genes[11,2]*DataLUADsig[rownames(LUAD_genes)[11],x]+LUAD_genes[12,2]*DataLUADsig[rownames(LUAD_genes)[12],x]+LUAD_genes[13,2]*DataLUADsig[rownames(LUAD_genes)[13],x]+LUAD_genes[14,2]*DataLUADsig[rownames(LUAD_genes)[14],x]+LUAD_genes[15,2]*DataLUADsig[rownames(LUAD_genes)[15],x]+LUAD_genes[16,2]*DataLUADsig[rownames(LUAD_genes)[16],x]})
median_stratification<-ifelse(score>quantile(score,0.50),"hypoxia","normoxia")
luadsig_median_stratification<-as.data.frame(median_stratification)
####median stratification
survival_luad<-merge(GSE19188_ad_survival,luadsig_median_stratification,by="row.names")
colnames(survival_luad)[ncol(survival_luad)]<-"median_stratification"
survival_luad$median_stratification<-factor(survival_luad$median_stratification,levels=c("normoxia","hypoxia"))

fit <- survfit(Surv(censored_os_time,censored_os_status)~median_stratification, data=survival_luad)
par(pty="s")
gplotmeantro9<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantro9$plot <- gplotmeantro9$plot + labs(
  title    = "GSE19188 (Sun signature)"        
)
gplotmeantro9 <- ggpar(
  gplotmeantro9,
  font.title    = c(25, "bold"))
gplotmeantro9$plot<-gplotmeantro9$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantro9$table <- ggpar(gplotmeantro9$table,
                            font.title = list(size = 16))

###GSE41271
GSE41271<-readRDS(file = "GSE41271.rds")
###processing the survival data for GSE42127
#write.table(GSE41271_survival,file="GSE41271_survival.csv",sep=",",row.names=FALSE)
survivalGSE41271<-read.csv(file="GSE41271_survival.csv",sep=",",header=T)
####censoring data
survivalGSE41271$status<-ifelse(survivalGSE41271$vital.statistics.ch1=="A",0,1)
censorship<-60
survivalGSE41271$OS   <- as.numeric( as.character(survivalGSE41271$OS ) )
GSE41271_survival_event_censorship  <- ifelse(survivalGSE41271$OS <= censorship &survivalGSE41271$status == 1 , 1 ,0  )
GSE41271_survival_time_censorship   <- ifelse( GSE41271_survival_event_censorship == 0 & survivalGSE41271$OS>=censorship ,censorship ,
survivalGSE41271$OS )   
survivalGSE41271   <- cbind(survivalGSE41271, GSE41271_survival_event_censorship, GSE41271_survival_time_censorship )
colnames(survivalGSE41271)[ (ncol(survivalGSE41271)-1):ncol(survivalGSE41271) ] <- c("c_event","c_event_time")

survivalGSE41271$rstatus<-ifelse(survivalGSE41271$recurrence.ch1=="N",0,1)
censorship<-60
survivalGSE41271$RFS   <- as.numeric( as.character(survivalGSE41271$RFS ) )
GSE41271_rsurvival_event_censorship  <- ifelse(survivalGSE41271$RFS <= censorship &survivalGSE41271$rstatus == 1 , 1 ,0  )
GSE41271_rsurvival_time_censorship   <- ifelse( GSE41271_rsurvival_event_censorship == 0 & survivalGSE41271$RFS>=censorship ,censorship ,
survivalGSE41271$RFS )   
survivalGSE41271   <- cbind(survivalGSE41271, GSE41271_rsurvival_event_censorship, GSE41271_rsurvival_time_censorship )
colnames(survivalGSE41271)[ (ncol(survivalGSE41271)-1):ncol(survivalGSE41271) ] <- c("c_revent","c_revent_time")
GSE41271_survival<-survivalGSE41271[survivalGSE41271$histology.ch1=="Adenocarcinoma",]
rownames(GSE41271_survival)<-GSE41271_survival$geo_accession
exprs<-GSE41271[,intersect(colnames(GSE41271),rownames(GSE41271_survival))]
exprs<-log2(exprs)
sig_exprs<-madSummaryAgilent_Sun(illuminaHumanv3.db,seed_genes,exprs)
DataLUADsig<-sig_exprs[rownames(LUAD_genes),]
score<-sapply(colnames(DataLUADsig),function(x){
	y<-LUAD_genes[1,2]*DataLUADsig[rownames(LUAD_genes)[1],x]+
	LUAD_genes[2,2]*DataLUADsig[rownames(LUAD_genes)[2],x]+
	LUAD_genes[3,2]*DataLUADsig[rownames(LUAD_genes)[3],x]
+LUAD_genes[4,2]*DataLUADsig[rownames(LUAD_genes)[4],x]+LUAD_genes[5,2]*DataLUADsig[rownames(LUAD_genes)[5],x]+LUAD_genes[6,2]*DataLUADsig[rownames(LUAD_genes)[6],x]+LUAD_genes[7,2]*DataLUADsig[rownames(LUAD_genes)[7],x]+LUAD_genes[8,2]*DataLUADsig[rownames(LUAD_genes)[8],x]+LUAD_genes[9,2]*DataLUADsig[rownames(LUAD_genes)[9],x]+LUAD_genes[10,2]*DataLUADsig[rownames(LUAD_genes)[10],x]+LUAD_genes[11,2]*DataLUADsig[rownames(LUAD_genes)[11],x]+LUAD_genes[12,2]*DataLUADsig[rownames(LUAD_genes)[12],x]+LUAD_genes[13,2]*DataLUADsig[rownames(LUAD_genes)[13],x]+LUAD_genes[14,2]*DataLUADsig[rownames(LUAD_genes)[14],x]+LUAD_genes[15,2]*DataLUADsig[rownames(LUAD_genes)[15],x]+LUAD_genes[16,2]*DataLUADsig[rownames(LUAD_genes)[16],x]})
median_stratification<-ifelse(score>quantile(score,0.50),"hypoxia","normoxia")
luadsig_median_stratification<-as.data.frame(median_stratification)
####median stratification
survival_luad<-merge(GSE41271_survival,luadsig_median_stratification,by="row.names")
colnames(survival_luad)[ncol(survival_luad)]<-"median_stratification"
survival_luad$median_stratification<-factor(survival_luad$median_stratification,levels=c("normoxia","hypoxia"))

res.cox <- coxph(Surv(c_event_time,c_event)~median_stratification, data=survival_luad)
summary(res.cox)
fit <- survfit(Surv(c_event_time,c_event)~median_stratification, data=survival_luad)
par(pty="s")
gplotmeantro3<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantro3$plot <- gplotmeantro3$plot + labs(
  title    = "GSE41271 (Sun signature)"        
)
gplotmeantro3 <- ggpar(
  gplotmeantro3,
  font.title    = c(25, "bold"))
gplotmeantro3$plot<-gplotmeantro3$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantro3$table <- ggpar(gplotmeantro3$table,
                            font.title = list(size = 16))

#####
GSE29013<-readRDS(file = "GSE29013.rds")
GSE29013_survival<-readRDS(file = "GSE29013_survival.rds")
GSE29013_survival$time_to_progression<-as.numeric(GSE29013_survival$time_to_progression) * 12
GSE29013_survival$death_time<-as.numeric(GSE29013_survival$death_time) * 12

####censoring data

censorship<-60
GSE29013_survival_event_censorship  <- ifelse(GSE29013_survival$death_time<= censorship &GSE29013_survival$death_event == 1 , 1 ,0  )
GSE29013_survival_time_censorship   <- ifelse(GSE29013_survival_event_censorship == 0 & GSE29013_survival$death_time>=censorship ,censorship ,
GSE29013_survival$death_time )   
GSE29013_survival  <- cbind(GSE29013_survival, GSE29013_survival_event_censorship, GSE29013_survival_time_censorship )
colnames(GSE29013_survival)[ (ncol(GSE29013_survival)-1):ncol(GSE29013_survival) ] <- c("c_event","c_event_time")

GSE29013_rsurvival_event_censorship  <- ifelse(GSE29013_survival$time_to_progression<= censorship &GSE29013_survival$progression == 1 , 1 ,0  )
GSE29013_rsurvival_time_censorship   <- ifelse(GSE29013_rsurvival_event_censorship == 0 & GSE29013_survival$time_to_progression>=censorship ,censorship ,
GSE29013_survival$time_to_progression)   
GSE29013_survival  <- cbind(GSE29013_survival, GSE29013_rsurvival_event_censorship, GSE29013_rsurvival_time_censorship )
colnames(GSE29013_survival)[ (ncol(GSE29013_survival)-1):ncol(GSE29013_survival) ] <- c("c_revent","c_revent_time")


GSE29013_survival<-GSE29013_survival[GSE29013_survival$histology=="Adenocarcinoma",]

exprs<-GSE29013[,intersect(colnames(GSE29013),rownames(GSE29013_survival))]
sig_exprs<-madSummary133plus2_Sun(hgu133plus2.db,seed_genes,exprs)
DataLUADsig<-sig_exprs[rownames(LUAD_genes),]
score<-sapply(colnames(DataLUADsig),function(x){
	y<-LUAD_genes[1,2]*DataLUADsig[rownames(LUAD_genes)[1],x]+
	LUAD_genes[2,2]*DataLUADsig[rownames(LUAD_genes)[2],x]+
	LUAD_genes[3,2]*DataLUADsig[rownames(LUAD_genes)[3],x]
+LUAD_genes[4,2]*DataLUADsig[rownames(LUAD_genes)[4],x]+LUAD_genes[5,2]*DataLUADsig[rownames(LUAD_genes)[5],x]+LUAD_genes[6,2]*DataLUADsig[rownames(LUAD_genes)[6],x]+LUAD_genes[7,2]*DataLUADsig[rownames(LUAD_genes)[7],x]+LUAD_genes[8,2]*DataLUADsig[rownames(LUAD_genes)[8],x]+LUAD_genes[9,2]*DataLUADsig[rownames(LUAD_genes)[9],x]+LUAD_genes[10,2]*DataLUADsig[rownames(LUAD_genes)[10],x]+LUAD_genes[11,2]*DataLUADsig[rownames(LUAD_genes)[11],x]+LUAD_genes[12,2]*DataLUADsig[rownames(LUAD_genes)[12],x]+LUAD_genes[13,2]*DataLUADsig[rownames(LUAD_genes)[13],x]+LUAD_genes[14,2]*DataLUADsig[rownames(LUAD_genes)[14],x]+LUAD_genes[15,2]*DataLUADsig[rownames(LUAD_genes)[15],x]+LUAD_genes[16,2]*DataLUADsig[rownames(LUAD_genes)[16],x]})
median_stratification<-ifelse(score>quantile(score,0.50),"hypoxia","normoxia")
luadsig_median_stratification<-as.data.frame(median_stratification)
####median stratification
survival_luad<-merge(GSE29013_survival,luadsig_median_stratification,by="row.names")
colnames(survival_luad)[ncol(survival_luad)]<-"median_stratification"
survival_luad$median_stratification<-factor(survival_luad$median_stratification,levels=c("normoxia","hypoxia"))



res.cox <- coxph(Surv(c_event_time,c_event)~median_stratification, data=survival_luad )
summary(res.cox)
fit <- survfit(Surv(c_event_time,c_event)~median_stratification, data=survival_luad )
par(pty="s")
gplotmeantro4<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantro4$plot <- gplotmeantro4$plot + labs(
  title    = "GSE29013 (Sun signature)"        
)
gplotmeantro4 <- ggpar(
  gplotmeantro4,
  font.title    = c(25, "bold"))
gplotmeantro4$plot<-gplotmeantro4$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantro4$table <- ggpar(gplotmeantro4$table,
                            font.title = list(size = 16))



######
GSE50081<-readRDS(file = "GSE50081.rds")
GSE50081_survival<-readRDS(file = "GSE50081_survival.rds")

GSE50081_survival[,54]<-as.numeric(as.character(GSE50081_survival[,54]))*12
GSE50081_survival$Status<-ifelse(GSE50081_survival[,53]=="alive",0,1)
censorship<-60
GSE50081_survival_event_censorship  <- ifelse(GSE50081_survival[,54]<= censorship &GSE50081_survival$Status == 1 , 1 ,0  )
GSE50081_survival_time_censorship   <- ifelse(GSE50081_survival_event_censorship == 0 & GSE50081_survival[,54]>=censorship ,censorship ,
GSE50081_survival[,54])   
GSE50081_survival  <- cbind(GSE50081_survival, GSE50081_survival_event_censorship, GSE50081_survival_time_censorship )
colnames(GSE50081_survival)[ (ncol(GSE50081_survival)-1):ncol(GSE50081_survival) ] <- c("c_event","c_event_time")

GSE50081_survival<-GSE50081_survival[GSE50081_survival$histology=="adenocarcinoma",]
exprs<-GSE50081[,intersect(colnames(GSE50081),rownames(GSE50081_survival))]
sig_exprs<-madSummary133plus2_Sun(hgu133plus2.db,seed_genes,exprs)
DataLUADsig<-sig_exprs[rownames(LUAD_genes),]
score<-sapply(colnames(DataLUADsig),function(x){
	y<-LUAD_genes[1,2]*DataLUADsig[rownames(LUAD_genes)[1],x]+
	LUAD_genes[2,2]*DataLUADsig[rownames(LUAD_genes)[2],x]+
	LUAD_genes[3,2]*DataLUADsig[rownames(LUAD_genes)[3],x]
+LUAD_genes[4,2]*DataLUADsig[rownames(LUAD_genes)[4],x]+LUAD_genes[5,2]*DataLUADsig[rownames(LUAD_genes)[5],x]+LUAD_genes[6,2]*DataLUADsig[rownames(LUAD_genes)[6],x]+LUAD_genes[7,2]*DataLUADsig[rownames(LUAD_genes)[7],x]+LUAD_genes[8,2]*DataLUADsig[rownames(LUAD_genes)[8],x]+LUAD_genes[9,2]*DataLUADsig[rownames(LUAD_genes)[9],x]+LUAD_genes[10,2]*DataLUADsig[rownames(LUAD_genes)[10],x]+LUAD_genes[11,2]*DataLUADsig[rownames(LUAD_genes)[11],x]+LUAD_genes[12,2]*DataLUADsig[rownames(LUAD_genes)[12],x]+LUAD_genes[13,2]*DataLUADsig[rownames(LUAD_genes)[13],x]+LUAD_genes[14,2]*DataLUADsig[rownames(LUAD_genes)[14],x]+LUAD_genes[15,2]*DataLUADsig[rownames(LUAD_genes)[15],x]+LUAD_genes[16,2]*DataLUADsig[rownames(LUAD_genes)[16],x]})
median_stratification<-ifelse(score>quantile(score,0.50),"hypoxia","normoxia")
luadsig_median_stratification<-as.data.frame(median_stratification)
####median stratification
survival_luad<-merge(GSE50081_survival,luadsig_median_stratification,by="row.names")
colnames(survival_luad)[ncol(survival_luad)]<-"median_stratification"
survival_luad$median_stratification<-factor(survival_luad$median_stratification,levels=c("normoxia","hypoxia"))

res.cox <- coxph(Surv(c_event_time,c_event)~median_stratification, data=survival_luad )
summary(res.cox)
fit <- survfit(Surv(c_event_time,c_event)~median_stratification, data=survival_luad )
par(pty="s")
gplotmeantro5<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantro5$plot <- gplotmeantro5$plot + labs(
  title    = "GSE50081 (Sun signature)"        
)
gplotmeantro5 <- ggpar(
  gplotmeantro5,
  font.title    = c(25, "bold"))
gplotmeantro5$plot<-gplotmeantro5$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantro5$table <- ggpar(gplotmeantro5$table,
                            font.title = list(size = 16))

###mva


res.cox <- coxph(Surv(c_event_time,c_event)~as.numeric(survival_luad[,45]), data=survival_luad )
summary(res.cox)

res.cox <- coxph(Surv(c_event_time,c_event)~survival_luad[,51], data=survival_luad )
summary(res.cox)

GSE50081_survival_edited<-survival_luad [!survival_luad[,52]=="Unable to determine",]
GSE50081_survival_edited$smoking<-ifelse(GSE50081_survival_edited$smoking=="Never",0,1)
res.cox <- coxph(Surv(c_event_time,c_event)~smoking, data=GSE50081_survival_edited )
summary(res.cox)

survival_luad$Stage<-sapply(survival_luad[,53],function(x){
	if(x=="1A"|x=="1B"){y=1} else{y=2}})


res.cox <- coxph(Surv(c_event_time,c_event)~Stage, data=survival_luad )
summary(res.cox)

res.cox <- coxph(Surv(c_event_time,c_event)~Stage+median_stratification, data=survival_luad)
summary(res.cox)



#######
GSE30219<-readRDS(file = "GSE30219.rds")
GSE30219_survival<-readRDS(file = "GSE30219_survival.rds")
GSE30219_survival[,43]<-as.numeric(GSE30219_survival[,43]) 
GSE30219_survival<-GSE30219_survival[!is.na(GSE30219_survival[,43]),]
GSE30219_survival$Status<-ifelse(GSE30219_survival[,50]=="ALIVE",0,1)
####censoring data

censorship<-60
GSE30219_survival_event_censorship  <- ifelse(GSE30219_survival[,43]<= censorship &GSE30219_survival$Status == 1 , 1 ,0  )
GSE30219_survival_time_censorship   <- ifelse(GSE30219_survival_event_censorship == 0 & GSE30219_survival[,43]>=censorship ,censorship ,
GSE30219_survival[,43] )   
GSE30219_survival  <- cbind(GSE30219_survival, GSE30219_survival_event_censorship, GSE30219_survival_time_censorship )
colnames(GSE30219_survival)[ (ncol(GSE30219_survival)-1):ncol(GSE30219_survival) ] <- c("c_event","c_event_time")
GSE30219_survival<-GSE30219_survival[GSE30219_survival[,45]=="ADC",]

exprs<-GSE30219[,intersect(colnames(GSE30219),rownames(GSE30219_survival))]

sig_exprs<-madSummary133plus2_Sun(hgu133plus2.db,seed_genes,exprs)
DataLUADsig<-sig_exprs[rownames(LUAD_genes),]
score<-sapply(colnames(DataLUADsig),function(x){
	y<-LUAD_genes[1,2]*DataLUADsig[rownames(LUAD_genes)[1],x]+
	LUAD_genes[2,2]*DataLUADsig[rownames(LUAD_genes)[2],x]+
	LUAD_genes[3,2]*DataLUADsig[rownames(LUAD_genes)[3],x]
+LUAD_genes[4,2]*DataLUADsig[rownames(LUAD_genes)[4],x]+LUAD_genes[5,2]*DataLUADsig[rownames(LUAD_genes)[5],x]+LUAD_genes[6,2]*DataLUADsig[rownames(LUAD_genes)[6],x]+LUAD_genes[7,2]*DataLUADsig[rownames(LUAD_genes)[7],x]+LUAD_genes[8,2]*DataLUADsig[rownames(LUAD_genes)[8],x]+LUAD_genes[9,2]*DataLUADsig[rownames(LUAD_genes)[9],x]+LUAD_genes[10,2]*DataLUADsig[rownames(LUAD_genes)[10],x]+LUAD_genes[11,2]*DataLUADsig[rownames(LUAD_genes)[11],x]+LUAD_genes[12,2]*DataLUADsig[rownames(LUAD_genes)[12],x]+LUAD_genes[13,2]*DataLUADsig[rownames(LUAD_genes)[13],x]+LUAD_genes[14,2]*DataLUADsig[rownames(LUAD_genes)[14],x]+LUAD_genes[15,2]*DataLUADsig[rownames(LUAD_genes)[15],x]+LUAD_genes[16,2]*DataLUADsig[rownames(LUAD_genes)[16],x]})
median_stratification<-ifelse(score>quantile(score,0.50),"hypoxia","normoxia")
luadsig_median_stratification<-as.data.frame(median_stratification)
####median stratification
survival_luad<-merge(GSE30219_survival,luadsig_median_stratification,by="row.names")
colnames(survival_luad)[ncol(survival_luad)]<-"median_stratification"
survival_luad$median_stratification<-factor(survival_luad$median_stratification,levels=c("normoxia","hypoxia"))



res.cox <- coxph(Surv(c_event_time,c_event)~median_stratification, data=survival_luad )
summary(res.cox)
fit <- survfit(Surv(c_event_time,c_event)~ median_stratification, data=survival_luad )
par(pty="s")
gplotmeantro6<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantro6$plot <- gplotmeantro6$plot + labs(
  title    = "GSE30219 (Sun signature)"        
)
gplotmeantro6 <- ggpar(
  gplotmeantro6,
  font.title    = c(25, "bold"))
gplotmeantro6$plot<-gplotmeantro6$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantro6$table <- ggpar(gplotmeantro6$table,
                            font.title = list(size = 16))
                            
                            

load( "GSE42127_LUAD_survival.RData")
load("GSE42127_exprs.RData")


#GSE42127_ad_survival<-GSE42127_ad_survival[which(GSE42127_ad_survival$had_adjuvant_chemo==T),]


exprs<-exprs[,intersect(colnames(exprs),rownames(GSE42127_ad_survival))]
sig_exprs<-madSummaryAgilent_Sun(illuminaHumanv3.db,seed_genes,exprs)
DataLUADsig<-sig_exprs[rownames(LUAD_genes),]
score<-sapply(colnames(DataLUADsig),function(x){
	y<-LUAD_genes[1,2]*DataLUADsig[rownames(LUAD_genes)[1],x]+
	LUAD_genes[2,2]*DataLUADsig[rownames(LUAD_genes)[2],x]+
	LUAD_genes[3,2]*DataLUADsig[rownames(LUAD_genes)[3],x]
+LUAD_genes[4,2]*DataLUADsig[rownames(LUAD_genes)[4],x]+LUAD_genes[5,2]*DataLUADsig[rownames(LUAD_genes)[5],x]+LUAD_genes[6,2]*DataLUADsig[rownames(LUAD_genes)[6],x]+LUAD_genes[7,2]*DataLUADsig[rownames(LUAD_genes)[7],x]+LUAD_genes[8,2]*DataLUADsig[rownames(LUAD_genes)[8],x]+LUAD_genes[9,2]*DataLUADsig[rownames(LUAD_genes)[9],x]+LUAD_genes[10,2]*DataLUADsig[rownames(LUAD_genes)[10],x]+LUAD_genes[11,2]*DataLUADsig[rownames(LUAD_genes)[11],x]+LUAD_genes[12,2]*DataLUADsig[rownames(LUAD_genes)[12],x]+LUAD_genes[13,2]*DataLUADsig[rownames(LUAD_genes)[13],x]+LUAD_genes[14,2]*DataLUADsig[rownames(LUAD_genes)[14],x]+LUAD_genes[15,2]*DataLUADsig[rownames(LUAD_genes)[15],x]+LUAD_genes[16,2]*DataLUADsig[rownames(LUAD_genes)[16],x]})
median_stratification<-ifelse(score>quantile(score,0.50),"hypoxia","normoxia")
luadsig_median_stratification<-as.data.frame(median_stratification)
####median stratification
survival_luad<-merge(GSE42127_ad_survival,luadsig_median_stratification,by="row.names")
colnames(survival_luad)[ncol(survival_luad)]<-"median_stratification"
survival_luad$median_stratification<-factor(survival_luad$median_stratification,levels=c("normoxia","hypoxia"))

GSE42127_fit<-survfit(Surv(censored_os_time,censored_os_status)~median_stratification, data=survival_luad)
fit <- survfit(Surv(censored_os_time,censored_os_status)~ median_stratification, data=survival_luad)
par(pty="s")
gplotmeantro8<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantro8$plot <- gplotmeantro8$plot + labs(
  title    = "GSE42127 (Sun signature)"        
)
gplotmeantro8 <- ggpar(
  gplotmeantro8,
  font.title    = c(25, "bold"))
gplotmeantro8$plot<-gplotmeantro8$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantro8$table <- ggpar(gplotmeantro8$table,
                            font.title = list(size = 16))                            

splots<-list(gplotmeantro1,gplotmeantro3,gplotmeantro4,gplotmeantro5,gplotmeantro6,gplotmeantro7,gplotmeantro8,gplotmeantro9)
s<-arrange_ggsurvplots(splots,ncol=3,nrow=3)
ggsave("/Users/mkhan/Documents/Brian RNA-seq scripts/Figure S_Sun.pdf",plot=s,height=18,width=18)


###### lung Shi signature
###Look at the different LUAD signatures performance in the different cohorts individually
LUAD_genes<-read.table("/Users/mkhan/Documents/PhD work/Brian RNA-seq scripts/luad_stage1signature.csv",sep=",",header=T)
rownames(LUAD_genes)<-LUAD_genes$Gene
load("GSE8894_survival.RData")
load("GSE8894_exprs.RData")
###only work on the lung adenocarcinoma
GSE8894_survival<-GSE8894_survival[GSE8894_survival$cell_type=="Adenocarcinoma",]

exprs<-GSE8894_expr[,intersect(colnames(GSE8894_expr),rownames(GSE8894_survival))]###this line is finding those samples which are common between patient samples clinical and 
####Censor the data
GSE8894_survival$status   <-ifelse(GSE8894_survival$RFS_status=="non_recurrence",0,1)
censorship<-60
GSE8894_survival$rfs_month   <- as.numeric( as.character(GSE8894_survival$rfs_month) )
GSE8894_survival_event_censorship  <- ifelse(GSE8894_survival$rfs_month<= censorship &GSE8894_survival$status == 1 , 1 ,0  )
GSE8894_survival_time_censorship   <- ifelse( GSE8894_survival_event_censorship == 0 & GSE8894_survival$rfs_month>=censorship ,censorship ,GSE8894_survival$rfs_month )   
GSE8894_survival    <- cbind( GSE8894_survival, GSE8894_survival_event_censorship, GSE8894_survival_time_censorship )
colnames(GSE8894_survival)[ (ncol(GSE8894_survival)-1):ncol(GSE8894_survival) ] <- c("c_event","c_event_time")
####KM plot
DataLUADsig<-exprs[rownames(LUAD_genes),]
score<-sapply(colnames(DataLUADsig),function(x){
	y<-LUAD_genes[1,2]*DataLUADsig[rownames(LUAD_genes)[1],x]+
	LUAD_genes[2,2]*DataLUADsig[rownames(LUAD_genes)[2],x]+
	LUAD_genes[3,2]*DataLUADsig[rownames(LUAD_genes)[3],x]
+LUAD_genes[4,2]*DataLUADsig[rownames(LUAD_genes)[4],x]+LUAD_genes[5,2]*DataLUADsig[rownames(LUAD_genes)[5],x]+LUAD_genes[6,2]*DataLUADsig[rownames(LUAD_genes)[6],x]+LUAD_genes[7,2]*DataLUADsig[rownames(LUAD_genes)[7],x]+LUAD_genes[8,2]*DataLUADsig[rownames(LUAD_genes)[8],x]+LUAD_genes[9,2]*DataLUADsig[rownames(LUAD_genes)[9],x]+LUAD_genes[10,2]*DataLUADsig[rownames(LUAD_genes)[10],x]})
###normalise to z-score
score<-scale(score)


stratification<-ifelse(score>0,"hypoxia","normoxia")
stratification<-as.data.frame(stratification)
####median stratification
survival_luad<-merge(GSE8894_survival,stratification,by="row.names")
colnames(survival_luad)[ncol(survival_luad)]<-"stratification"
survival_luad$stratification<-factor(survival_luad$stratification,levels=c("normoxia","hypoxia"))
res.cox <- coxph(Surv(c_event_time,c_event) ~ stratification, data=survival_luad)
summary(res.cox)#1.83[0.91-3.69] p=0.09

###GSE3141
#HGU133plus2
load( "GSE3141_survival.RData")
load("GSE3141_exprs.RData")

GSE3141_survival$status   <-ifelse(GSE3141_survival$Status=="alive",0,1)
censorship<-60
GSE3141_survival$Overall_survival_time    <- as.numeric( as.character(GSE3141_survival$Overall_survival_time ) )
GSE3141_survival_event_censorship  <- ifelse(GSE3141_survival$Overall_survival_time <= censorship & GSE3141_survival$status == 1 , 1 ,0  )
GSE3141_survival_time_censorship   <- ifelse( GSE3141_survival_event_censorship == 0 & GSE3141_survival$Overall_survival_time >=censorship ,censorship ,GSE3141_survival$Overall_survival_time )   
GSE3141_survival     <- cbind( GSE3141_survival, GSE3141_survival_event_censorship, GSE3141_survival_time_censorship )
colnames(GSE3141_survival)[ (ncol(GSE3141_survival)-1):ncol(GSE3141_survival) ] <- c("c_event","c_event_time")
GSE3141_survival$Cell.type<-as.character(GSE3141_survival$Cell.type)
GSE3141_survival<-GSE3141_survival[GSE3141_survival$Cell.type==" A",]
exprs<-GSE3141_expr[,intersect(colnames(GSE3141_expr),rownames(GSE3141_survival))]
DataLUADsig<-exprs[rownames(LUAD_genes),]
score<-sapply(colnames(DataLUADsig),function(x){
	y<-LUAD_genes[1,2]*DataLUADsig[rownames(LUAD_genes)[1],x]+
	LUAD_genes[2,2]*DataLUADsig[rownames(LUAD_genes)[2],x]+
	LUAD_genes[3,2]*DataLUADsig[rownames(LUAD_genes)[3],x]
+LUAD_genes[4,2]*DataLUADsig[rownames(LUAD_genes)[4],x]+LUAD_genes[5,2]*DataLUADsig[rownames(LUAD_genes)[5],x]+LUAD_genes[6,2]*DataLUADsig[rownames(LUAD_genes)[6],x]+LUAD_genes[7,2]*DataLUADsig[rownames(LUAD_genes)[7],x]+LUAD_genes[8,2]*DataLUADsig[rownames(LUAD_genes)[8],x]+LUAD_genes[9,2]*DataLUADsig[rownames(LUAD_genes)[9],x]+LUAD_genes[10,2]*DataLUADsig[rownames(LUAD_genes)[10],x]})
score<-scale(score)
median_stratification<-ifelse(score>0,"hypoxia","normoxia")
luadsig_median_stratification<-as.data.frame(median_stratification)
####median stratification
survival_luad<-merge(GSE3141_survival,luadsig_median_stratification,by="row.names")
colnames(survival_luad)[ncol(survival_luad)]<-"median_stratification"
survival_luad$median_stratification<-factor(survival_luad$median_stratification,levels=c("normoxia","hypoxia"))
res.cox <- coxph(Surv(c_event_time,c_event) ~ median_stratification, data=survival_luad)
summary(res.cox)
fit <- survfit(Surv(c_event_time,c_event)~ median_stratification, data=survival_luad)
par(pty="s")
gplotmeantro1<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantro1$plot <- gplotmeantro1$plot + labs(
  title    = "GSE3141 (Shi signature)"        
)
gplotmeantro1 <- ggpar(
  gplotmeantro1,
  font.title    = c(25, "bold"))
gplotmeantro1$plot<-gplotmeantro1$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantro1$table <- ggpar(gplotmeantro1$table,
                            font.title = list(size = 16))


####Next dataset

load("GSE31210_survival.RData")#GSE31210_survival
load("GSE31210_exprs.RData")


exprs<-exprs[,intersect(colnames(exprs),rownames(GSE31210_survival))]
rownames(LUAD_genes)[10]<-"BCLAF3"
seed_genes<-rownames(LUAD_genes)
sig_exprs<-madSummary133plus2_Shi(hgu133plus2.db,seed_genes,exprs)

DataLUADsig<-sig_exprs[rownames(LUAD_genes),]
score<-sapply(colnames(DataLUADsig),function(x){
	y<-LUAD_genes[1,2]*DataLUADsig[rownames(LUAD_genes)[1],x]+
	LUAD_genes[2,2]*DataLUADsig[rownames(LUAD_genes)[2],x]+
	LUAD_genes[3,2]*DataLUADsig[rownames(LUAD_genes)[3],x]
+LUAD_genes[4,2]*DataLUADsig[rownames(LUAD_genes)[4],x]+LUAD_genes[5,2]*DataLUADsig[rownames(LUAD_genes)[5],x]+LUAD_genes[6,2]*DataLUADsig[rownames(LUAD_genes)[6],x]+LUAD_genes[7,2]*DataLUADsig[rownames(LUAD_genes)[7],x]+LUAD_genes[8,2]*DataLUADsig[rownames(LUAD_genes)[8],x]+LUAD_genes[9,2]*DataLUADsig[rownames(LUAD_genes)[9],x]+LUAD_genes[10,2]*DataLUADsig[rownames(LUAD_genes)[10],x]})
score<-scale(score)
median_stratification<-ifelse(score>0,"hypoxia","normoxia")
luadsig_median_stratification<-as.data.frame(median_stratification)
####median stratification
survival_luad<-merge(GSE31210_survival,luadsig_median_stratification,by="row.names")
colnames(survival_luad)[ncol(survival_luad)]<-"median_stratification"
survival_luad$median_stratification<-factor(survival_luad$median_stratification,levels=c("normoxia","hypoxia"))

res.cox<-coxph(Surv(censored_os_time,censored_os_status)~median_stratification, data=survival_luad)
summary(res.cox)


GSE31210_fit<-survfit(Surv(censored_os_time,censored_os_status)~median_stratification, data=survival_luad)
fit <- survfit(Surv(censored_os_time,censored_os_status)~median_stratification, data=survival_luad)
par(pty="s")
gplotmeantro7<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantro7$plot <- gplotmeantro7$plot + labs(
  title    = "GSE31210 (Shi signature)"        
)
gplotmeantro7 <- ggpar(
  gplotmeantro7,
  font.title    = c(25, "bold"))
gplotmeantro7$plot<-gplotmeantro7$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantro7$table <- ggpar(gplotmeantro7$table,
                            font.title = list(size = 16))

###

load("GSE19188_LUAD_survival.RData")
load("GSE19188_exprs.RData")


exprs<-exprs[,intersect(colnames(exprs),rownames(GSE19188_ad_survival))]
sig_exprs<-madSummary133plus2_Shi(hgu133plus2.db,seed_genes,exprs)
DataLUADsig<-sig_exprs[rownames(LUAD_genes),]
score<-sapply(colnames(DataLUADsig),function(x){
	y<-LUAD_genes[1,2]*DataLUADsig[rownames(LUAD_genes)[1],x]+
	LUAD_genes[2,2]*DataLUADsig[rownames(LUAD_genes)[2],x]+
	LUAD_genes[3,2]*DataLUADsig[rownames(LUAD_genes)[3],x]
+LUAD_genes[4,2]*DataLUADsig[rownames(LUAD_genes)[4],x]+LUAD_genes[5,2]*DataLUADsig[rownames(LUAD_genes)[5],x]+LUAD_genes[6,2]*DataLUADsig[rownames(LUAD_genes)[6],x]+LUAD_genes[7,2]*DataLUADsig[rownames(LUAD_genes)[7],x]+LUAD_genes[8,2]*DataLUADsig[rownames(LUAD_genes)[8],x]+LUAD_genes[9,2]*DataLUADsig[rownames(LUAD_genes)[9],x]+LUAD_genes[10,2]*DataLUADsig[rownames(LUAD_genes)[10],x]})
score<-scale(score)
median_stratification<-ifelse(score>0,"hypoxia","normoxia")
luadsig_median_stratification<-as.data.frame(median_stratification)
####median stratification
survival_luad<-merge(GSE19188_ad_survival,luadsig_median_stratification,by="row.names")
colnames(survival_luad)[ncol(survival_luad)]<-"median_stratification"
survival_luad$median_stratification<-factor(survival_luad$median_stratification,levels=c("normoxia","hypoxia"))

res.cox<- coxph(Surv(censored_os_time,censored_os_status)~median_stratification, data=survival_luad)
summary(res.cox)




fit <- survfit(Surv(censored_os_time,censored_os_status)~median_stratification, data=survival_luad)
par(pty="s")
gplotmeantro9<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantro9$plot <- gplotmeantro9$plot + labs(
  title    = "GSE19188 (Shi signature)"        
)
gplotmeantro9 <- ggpar(
  gplotmeantro9,
  font.title    = c(25, "bold"))
gplotmeantro9$plot<-gplotmeantro9$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantro9$table <- ggpar(gplotmeantro9$table,
                            font.title = list(size = 16))

###GSE41271
GSE41271<-readRDS(file = "GSE41271.rds")
###processing the survival data for GSE42127
#write.table(GSE41271_survival,file="GSE41271_survival.csv",sep=",",row.names=FALSE)
survivalGSE41271<-read.csv(file="GSE41271_survival.csv",sep=",",header=T)
####censoring data
survivalGSE41271$status<-ifelse(survivalGSE41271$vital.statistics.ch1=="A",0,1)
censorship<-60
survivalGSE41271$OS   <- as.numeric( as.character(survivalGSE41271$OS ) )
GSE41271_survival_event_censorship  <- ifelse(survivalGSE41271$OS <= censorship &survivalGSE41271$status == 1 , 1 ,0  )
GSE41271_survival_time_censorship   <- ifelse( GSE41271_survival_event_censorship == 0 & survivalGSE41271$OS>=censorship ,censorship ,
survivalGSE41271$OS )   
survivalGSE41271   <- cbind(survivalGSE41271, GSE41271_survival_event_censorship, GSE41271_survival_time_censorship )
colnames(survivalGSE41271)[ (ncol(survivalGSE41271)-1):ncol(survivalGSE41271) ] <- c("c_event","c_event_time")

survivalGSE41271$rstatus<-ifelse(survivalGSE41271$recurrence.ch1=="N",0,1)
censorship<-60
survivalGSE41271$RFS   <- as.numeric( as.character(survivalGSE41271$RFS ) )
GSE41271_rsurvival_event_censorship  <- ifelse(survivalGSE41271$RFS <= censorship &survivalGSE41271$rstatus == 1 , 1 ,0  )
GSE41271_rsurvival_time_censorship   <- ifelse( GSE41271_rsurvival_event_censorship == 0 & survivalGSE41271$RFS>=censorship ,censorship ,
survivalGSE41271$RFS )   
survivalGSE41271   <- cbind(survivalGSE41271, GSE41271_rsurvival_event_censorship, GSE41271_rsurvival_time_censorship )
colnames(survivalGSE41271)[ (ncol(survivalGSE41271)-1):ncol(survivalGSE41271) ] <- c("c_revent","c_revent_time")
GSE41271_survival<-survivalGSE41271[survivalGSE41271$histology.ch1=="Adenocarcinoma",]
rownames(GSE41271_survival)<-GSE41271_survival$geo_accession
exprs<-GSE41271[,intersect(colnames(GSE41271),rownames(GSE41271_survival))]
exprs<-log2(exprs)
sig_exprs<-madSummaryAgilent_Shi(illuminaHumanv3.db,seed_genes,exprs)
DataLUADsig<-sig_exprs[rownames(LUAD_genes),]
score<-sapply(colnames(DataLUADsig),function(x){
	y<-LUAD_genes[1,2]*DataLUADsig[rownames(LUAD_genes)[1],x]+
	LUAD_genes[2,2]*DataLUADsig[rownames(LUAD_genes)[2],x]+
	LUAD_genes[3,2]*DataLUADsig[rownames(LUAD_genes)[3],x]
+LUAD_genes[4,2]*DataLUADsig[rownames(LUAD_genes)[4],x]+LUAD_genes[5,2]*DataLUADsig[rownames(LUAD_genes)[5],x]+LUAD_genes[6,2]*DataLUADsig[rownames(LUAD_genes)[6],x]+LUAD_genes[7,2]*DataLUADsig[rownames(LUAD_genes)[7],x]+LUAD_genes[8,2]*DataLUADsig[rownames(LUAD_genes)[8],x]+LUAD_genes[9,2]*DataLUADsig[rownames(LUAD_genes)[9],x]+LUAD_genes[10,2]*DataLUADsig[rownames(LUAD_genes)[10],x]})
score<-scale(score)
median_stratification<-ifelse(score>0,"hypoxia","normoxia")
luadsig_median_stratification<-as.data.frame(median_stratification)
####median stratification
survival_luad<-merge(GSE41271_survival,luadsig_median_stratification,by="row.names")
colnames(survival_luad)[ncol(survival_luad)]<-"median_stratification"
survival_luad$median_stratification<-factor(survival_luad$median_stratification,levels=c("normoxia","hypoxia"))

res.cox <- coxph(Surv(c_event_time,c_event)~median_stratification, data=survival_luad)
summary(res.cox)
fit <- survfit(Surv(c_event_time,c_event)~median_stratification, data=survival_luad)
par(pty="s")
gplotmeantro3<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantro3$plot <- gplotmeantro3$plot + labs(
  title    = "GSE41271 (Shi signature)"        
)
gplotmeantro3 <- ggpar(
  gplotmeantro3,
  font.title    = c(25, "bold"))
gplotmeantro3$plot<-gplotmeantro3$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantro3$table <- ggpar(gplotmeantro3$table,
                            font.title = list(size = 16))

#####
GSE29013<-readRDS(file = "GSE29013.rds")
GSE29013_survival<-readRDS(file = "GSE29013_survival.rds")
GSE29013_survival$time_to_progression<-as.numeric(GSE29013_survival$time_to_progression) * 12
GSE29013_survival$death_time<-as.numeric(GSE29013_survival$death_time) * 12

####censoring data

censorship<-60
GSE29013_survival_event_censorship  <- ifelse(GSE29013_survival$death_time<= censorship &GSE29013_survival$death_event == 1 , 1 ,0  )
GSE29013_survival_time_censorship   <- ifelse(GSE29013_survival_event_censorship == 0 & GSE29013_survival$death_time>=censorship ,censorship ,
GSE29013_survival$death_time )   
GSE29013_survival  <- cbind(GSE29013_survival, GSE29013_survival_event_censorship, GSE29013_survival_time_censorship )
colnames(GSE29013_survival)[ (ncol(GSE29013_survival)-1):ncol(GSE29013_survival) ] <- c("c_event","c_event_time")

GSE29013_rsurvival_event_censorship  <- ifelse(GSE29013_survival$time_to_progression<= censorship &GSE29013_survival$progression == 1 , 1 ,0  )
GSE29013_rsurvival_time_censorship   <- ifelse(GSE29013_rsurvival_event_censorship == 0 & GSE29013_survival$time_to_progression>=censorship ,censorship ,
GSE29013_survival$time_to_progression)   
GSE29013_survival  <- cbind(GSE29013_survival, GSE29013_rsurvival_event_censorship, GSE29013_rsurvival_time_censorship )
colnames(GSE29013_survival)[ (ncol(GSE29013_survival)-1):ncol(GSE29013_survival) ] <- c("c_revent","c_revent_time")


GSE29013_survival<-GSE29013_survival[GSE29013_survival$histology=="Adenocarcinoma",]

exprs<-GSE29013[,intersect(colnames(GSE29013),rownames(GSE29013_survival))]
sig_exprs<-madSummary133plus2_Shi(hgu133plus2.db,seed_genes,exprs)
DataLUADsig<-sig_exprs[rownames(LUAD_genes),]
score<-sapply(colnames(DataLUADsig),function(x){
	y<-LUAD_genes[1,2]*DataLUADsig[rownames(LUAD_genes)[1],x]+
	LUAD_genes[2,2]*DataLUADsig[rownames(LUAD_genes)[2],x]+
	LUAD_genes[3,2]*DataLUADsig[rownames(LUAD_genes)[3],x]
+LUAD_genes[4,2]*DataLUADsig[rownames(LUAD_genes)[4],x]+LUAD_genes[5,2]*DataLUADsig[rownames(LUAD_genes)[5],x]+LUAD_genes[6,2]*DataLUADsig[rownames(LUAD_genes)[6],x]+LUAD_genes[7,2]*DataLUADsig[rownames(LUAD_genes)[7],x]+LUAD_genes[8,2]*DataLUADsig[rownames(LUAD_genes)[8],x]+LUAD_genes[9,2]*DataLUADsig[rownames(LUAD_genes)[9],x]+LUAD_genes[10,2]*DataLUADsig[rownames(LUAD_genes)[10],x]})
score<-scale(score)
median_stratification<-ifelse(score>0,"hypoxia","normoxia")
luadsig_median_stratification<-as.data.frame(median_stratification)
####median stratification
survival_luad<-merge(GSE29013_survival,luadsig_median_stratification,by="row.names")
colnames(survival_luad)[ncol(survival_luad)]<-"median_stratification"
survival_luad$median_stratification<-factor(survival_luad$median_stratification,levels=c("normoxia","hypoxia"))



res.cox <- coxph(Surv(c_event_time,c_event)~median_stratification, data=survival_luad )
summary(res.cox)
fit <- survfit(Surv(c_event_time,c_event)~median_stratification, data=survival_luad )
par(pty="s")
gplotmeantro4<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantro4$plot <- gplotmeantro4$plot + labs(
  title    = "GSE29013 (Shi signature)"        
)
gplotmeantro4 <- ggpar(
  gplotmeantro4,
  font.title    = c(25, "bold"))
gplotmeantro4$plot<-gplotmeantro4$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantro4$table <- ggpar(gplotmeantro4$table,
                            font.title = list(size = 16))



######
GSE50081<-readRDS(file = "GSE50081.rds")
GSE50081_survival<-readRDS(file = "GSE50081_survival.rds")

GSE50081_survival[,54]<-as.numeric(as.character(GSE50081_survival[,54]))*12
GSE50081_survival$Status<-ifelse(GSE50081_survival[,53]=="alive",0,1)
censorship<-60
GSE50081_survival_event_censorship  <- ifelse(GSE50081_survival[,54]<= censorship &GSE50081_survival$Status == 1 , 1 ,0  )
GSE50081_survival_time_censorship   <- ifelse(GSE50081_survival_event_censorship == 0 & GSE50081_survival[,54]>=censorship ,censorship ,
GSE50081_survival[,54])   
GSE50081_survival  <- cbind(GSE50081_survival, GSE50081_survival_event_censorship, GSE50081_survival_time_censorship )
colnames(GSE50081_survival)[ (ncol(GSE50081_survival)-1):ncol(GSE50081_survival) ] <- c("c_event","c_event_time")

GSE50081_survival<-GSE50081_survival[GSE50081_survival$histology=="adenocarcinoma",]
exprs<-GSE50081[,intersect(colnames(GSE50081),rownames(GSE50081_survival))]
sig_exprs<-madSummary133plus2_Shi(hgu133plus2.db,seed_genes,exprs)
DataLUADsig<-sig_exprs[rownames(LUAD_genes),]
score<-sapply(colnames(DataLUADsig),function(x){
	y<-LUAD_genes[1,2]*DataLUADsig[rownames(LUAD_genes)[1],x]+
	LUAD_genes[2,2]*DataLUADsig[rownames(LUAD_genes)[2],x]+
	LUAD_genes[3,2]*DataLUADsig[rownames(LUAD_genes)[3],x]
+LUAD_genes[4,2]*DataLUADsig[rownames(LUAD_genes)[4],x]+LUAD_genes[5,2]*DataLUADsig[rownames(LUAD_genes)[5],x]+LUAD_genes[6,2]*DataLUADsig[rownames(LUAD_genes)[6],x]+LUAD_genes[7,2]*DataLUADsig[rownames(LUAD_genes)[7],x]+LUAD_genes[8,2]*DataLUADsig[rownames(LUAD_genes)[8],x]+LUAD_genes[9,2]*DataLUADsig[rownames(LUAD_genes)[9],x]+LUAD_genes[10,2]*DataLUADsig[rownames(LUAD_genes)[10],x]})
score<-scale(score)
median_stratification<-ifelse(score>0,"hypoxia","normoxia")
luadsig_median_stratification<-as.data.frame(median_stratification)
####median stratification
survival_luad<-merge(GSE50081_survival,luadsig_median_stratification,by="row.names")
colnames(survival_luad)[ncol(survival_luad)]<-"median_stratification"
survival_luad$median_stratification<-factor(survival_luad$median_stratification,levels=c("normoxia","hypoxia"))

res.cox <- coxph(Surv(c_event_time,c_event)~median_stratification, data=survival_luad )
summary(res.cox)
fit <- survfit(Surv(c_event_time,c_event)~median_stratification, data=survival_luad )
par(pty="s")
gplotmeantro5<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantro5$plot <- gplotmeantro5$plot + labs(
  title    = "GSE50081 (Shi signature)"        
)
gplotmeantro5 <- ggpar(
  gplotmeantro5,
  font.title    = c(25, "bold"))
gplotmeantro5$plot<-gplotmeantro5$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantro5$table <- ggpar(gplotmeantro5$table,
                            font.title = list(size = 16))



###mva


res.cox <- coxph(Surv(c_event_time,c_event)~as.numeric(survival_luad[,45]), data=survival_luad )
summary(res.cox)

res.cox <- coxph(Surv(c_event_time,c_event)~survival_luad[,51], data=survival_luad )
summary(res.cox)

GSE50081_survival_edited<-survival_luad [!survival_luad[,52]=="Unable to determine",]
GSE50081_survival_edited$smoking<-ifelse(GSE50081_survival_edited$smoking=="Never",0,1)
res.cox <- coxph(Surv(c_event_time,c_event)~smoking, data=GSE50081_survival_edited )
summary(res.cox)

survival_luad$Stage<-sapply(survival_luad[,53],function(x){
	if(x=="1A"|x=="1B"){y=1} else{y=2}})


res.cox <- coxph(Surv(c_event_time,c_event)~Stage, data=survival_luad )
summary(res.cox)

res.cox <- coxph(Surv(c_event_time,c_event)~Stage+median_stratification, data=survival_luad)
summary(res.cox)





















#######
GSE30219<-readRDS(file = "GSE30219.rds")
GSE30219_survival<-readRDS(file = "GSE30219_survival.rds")
GSE30219_survival[,43]<-as.numeric(GSE30219_survival[,43]) 
GSE30219_survival<-GSE30219_survival[!is.na(GSE30219_survival[,43]),]
GSE30219_survival$Status<-ifelse(GSE30219_survival[,50]=="ALIVE",0,1)
####censoring data

censorship<-60
GSE30219_survival_event_censorship  <- ifelse(GSE30219_survival[,43]<= censorship &GSE30219_survival$Status == 1 , 1 ,0  )
GSE30219_survival_time_censorship   <- ifelse(GSE30219_survival_event_censorship == 0 & GSE30219_survival[,43]>=censorship ,censorship ,
GSE30219_survival[,43] )   
GSE30219_survival  <- cbind(GSE30219_survival, GSE30219_survival_event_censorship, GSE30219_survival_time_censorship )
colnames(GSE30219_survival)[ (ncol(GSE30219_survival)-1):ncol(GSE30219_survival) ] <- c("c_event","c_event_time")
GSE30219_survival<-GSE30219_survival[GSE30219_survival[,45]=="ADC",]

exprs<-GSE30219[,intersect(colnames(GSE30219),rownames(GSE30219_survival))]

sig_exprs<-madSummary133plus2_Shi(hgu133plus2.db,seed_genes,exprs)
DataLUADsig<-sig_exprs[rownames(LUAD_genes),]
score<-sapply(colnames(DataLUADsig),function(x){
	y<-LUAD_genes[1,2]*DataLUADsig[rownames(LUAD_genes)[1],x]+
	LUAD_genes[2,2]*DataLUADsig[rownames(LUAD_genes)[2],x]+
	LUAD_genes[3,2]*DataLUADsig[rownames(LUAD_genes)[3],x]
+LUAD_genes[4,2]*DataLUADsig[rownames(LUAD_genes)[4],x]+LUAD_genes[5,2]*DataLUADsig[rownames(LUAD_genes)[5],x]+LUAD_genes[6,2]*DataLUADsig[rownames(LUAD_genes)[6],x]+LUAD_genes[7,2]*DataLUADsig[rownames(LUAD_genes)[7],x]+LUAD_genes[8,2]*DataLUADsig[rownames(LUAD_genes)[8],x]+LUAD_genes[9,2]*DataLUADsig[rownames(LUAD_genes)[9],x]+LUAD_genes[10,2]*DataLUADsig[rownames(LUAD_genes)[10],x]})
score<-scale(score)
median_stratification<-ifelse(score>0,"hypoxia","normoxia")
luadsig_median_stratification<-as.data.frame(median_stratification)
####median stratification
survival_luad<-merge(GSE30219_survival,luadsig_median_stratification,by="row.names")
colnames(survival_luad)[ncol(survival_luad)]<-"median_stratification"
survival_luad$median_stratification<-factor(survival_luad$median_stratification,levels=c("normoxia","hypoxia"))



res.cox <- coxph(Surv(c_event_time,c_event)~median_stratification, data=survival_luad )
summary(res.cox)
fit <- survfit(Surv(c_event_time,c_event)~ median_stratification, data=survival_luad )
par(pty="s")
gplotmeantro6<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantro6$plot <- gplotmeantro6$plot + labs(
  title    = "GSE30219 (Shi signature)"        
)
gplotmeantro6 <- ggpar(
  gplotmeantro6,
  font.title    = c(25, "bold"))
gplotmeantro6$plot<-gplotmeantro6$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantro6$table <- ggpar(gplotmeantro6$table,
                            font.title = list(size = 16))
                            
                            

load( "GSE42127_LUAD_survival.RData")
load("GSE42127_exprs.RData")


#GSE42127_ad_survival<-GSE42127_ad_survival[which(GSE42127_ad_survival$had_adjuvant_chemo==T),]


exprs<-exprs[,intersect(colnames(exprs),rownames(GSE42127_ad_survival))]
sig_exprs<-madSummaryAgilent_Shi(illuminaHumanv3.db,seed_genes,exprs)
DataLUADsig<-sig_exprs[rownames(LUAD_genes),]
score<-sapply(colnames(DataLUADsig),function(x){
	y<-LUAD_genes[1,2]*DataLUADsig[rownames(LUAD_genes)[1],x]+
	LUAD_genes[2,2]*DataLUADsig[rownames(LUAD_genes)[2],x]+
	LUAD_genes[3,2]*DataLUADsig[rownames(LUAD_genes)[3],x]
+LUAD_genes[4,2]*DataLUADsig[rownames(LUAD_genes)[4],x]+LUAD_genes[5,2]*DataLUADsig[rownames(LUAD_genes)[5],x]+LUAD_genes[6,2]*DataLUADsig[rownames(LUAD_genes)[6],x]+LUAD_genes[7,2]*DataLUADsig[rownames(LUAD_genes)[7],x]+LUAD_genes[8,2]*DataLUADsig[rownames(LUAD_genes)[8],x]+LUAD_genes[9,2]*DataLUADsig[rownames(LUAD_genes)[9],x]+LUAD_genes[10,2]*DataLUADsig[rownames(LUAD_genes)[10],x]})
score<-scale(score)
median_stratification<-ifelse(score>0,"hypoxia","normoxia")
luadsig_median_stratification<-as.data.frame(median_stratification)
####median stratification
survival_luad<-merge(GSE42127_ad_survival,luadsig_median_stratification,by="row.names")
colnames(survival_luad)[ncol(survival_luad)]<-"median_stratification"
survival_luad$median_stratification<-factor(survival_luad$median_stratification,levels=c("normoxia","hypoxia"))

res.cox<-coxph(Surv(censored_os_time,censored_os_status)~median_stratification, data=survival_luad)
summary(res.cox)

GSE42127_fit<-survfit(Surv(censored_os_time,censored_os_status)~median_stratification, data=survival_luad)
fit <- survfit(Surv(censored_os_time,censored_os_status)~ median_stratification, data=survival_luad)
par(pty="s")
gplotmeantro8<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantro8$plot <- gplotmeantro8$plot + labs(
  title    = "GSE42127 (Shi signature)"        
)
gplotmeantro8 <- ggpar(
  gplotmeantro8,
  font.title    = c(25, "bold"))
gplotmeantro8$plot<-gplotmeantro8$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantro8$table <- ggpar(gplotmeantro8$table,
                            font.title = list(size = 16)) 


splots<-list(gplotmeantro1,gplotmeantro3,gplotmeantro4,gplotmeantro5,gplotmeantro6,gplotmeantro7,gplotmeantro8,gplotmeantro9)
s<-arrange_ggsurvplots(splots,ncol=3,nrow=3)
ggsave("/Users/mkhan/Documents/Brian RNA-seq scripts/Figure S_Shi.pdf",plot=s,height=18,width=18)


#####Buffa common signature
LUAD_genes<-read.table("/Users/mkhan/Documents/PhD work/Brian RNA-seq scripts/Buffa genes.csv",sep=",")
rownames(LUAD_genes)<-as.character(LUAD_genes[,1])

load("GSE8894_survival.RData")
load("GSE8894_exprs.RData")
###only work on the lung adenocarcinoma
GSE8894_survival<-GSE8894_survival[GSE8894_survival$cell_type=="Adenocarcinoma",]

exprs<-GSE8894_expr[,intersect(colnames(GSE8894_expr),rownames(GSE8894_survival))]###this line is finding those samples which are common between patient samples clinical and 
DataLUADsig<-exprs[rownames(LUAD_genes),]

####Censor the data
GSE8894_survival$status   <-ifelse(GSE8894_survival$RFS_status=="non_recurrence",0,1)
censorship<-60
GSE8894_survival$rfs_month   <- as.numeric( as.character(GSE8894_survival$rfs_month) )
GSE8894_survival_event_censorship  <- ifelse(GSE8894_survival$rfs_month<= censorship &GSE8894_survival$status == 1 , 1 ,0  )
GSE8894_survival_time_censorship   <- ifelse( GSE8894_survival_event_censorship == 0 & GSE8894_survival$rfs_month>=censorship ,censorship ,GSE8894_survival$rfs_month )   
GSE8894_survival    <- cbind( GSE8894_survival, GSE8894_survival_event_censorship, GSE8894_survival_time_censorship )
colnames(GSE8894_survival)[ (ncol(GSE8894_survival)-1):ncol(GSE8894_survival) ] <- c("c_event","c_event_time")
####KM plot
score<-apply(DataLUADsig,2,median)
###normalise to z-score

stratification<-ifelse(score>median(score),"hypoxia","normoxia")
stratification<-as.data.frame(stratification)
####median stratification
survival_luad<-merge(GSE8894_survival,stratification,by="row.names")
colnames(survival_luad)[ncol(survival_luad)]<-"stratification"
survival_luad$stratification<-factor(survival_luad$stratification,levels=c("normoxia","hypoxia"))
res.cox <- coxph(Surv(c_event_time,c_event) ~ stratification, data=survival_luad)
summary(res.cox)
###GSE3141
#HGU133plus2
load( "GSE3141_survival.RData")
load("GSE3141_exprs.RData")

GSE3141_survival$status   <-ifelse(GSE3141_survival$Status=="alive",0,1)
censorship<-60
GSE3141_survival$Overall_survival_time    <- as.numeric( as.character(GSE3141_survival$Overall_survival_time ) )
GSE3141_survival_event_censorship  <- ifelse(GSE3141_survival$Overall_survival_time <= censorship & GSE3141_survival$status == 1 , 1 ,0  )
GSE3141_survival_time_censorship   <- ifelse( GSE3141_survival_event_censorship == 0 & GSE3141_survival$Overall_survival_time >=censorship ,censorship ,GSE3141_survival$Overall_survival_time )   
GSE3141_survival     <- cbind( GSE3141_survival, GSE3141_survival_event_censorship, GSE3141_survival_time_censorship )
colnames(GSE3141_survival)[ (ncol(GSE3141_survival)-1):ncol(GSE3141_survival) ] <- c("c_event","c_event_time")
GSE3141_survival$Cell.type<-as.character(GSE3141_survival$Cell.type)
GSE3141_survival<-GSE3141_survival[GSE3141_survival$Cell.type==" A",]
exprs<-GSE3141_expr[,intersect(colnames(GSE3141_expr),rownames(GSE3141_survival))]
DataLUADsig<-exprs[rownames(LUAD_genes),]
score<-apply(DataLUADsig,2,median)
###normalise to z-score

median_stratification<-ifelse(score>median(score),"hypoxia","normoxia")
luadsig_median_stratification<-as.data.frame(median_stratification)
####median stratification
survival_luad<-merge(GSE3141_survival,luadsig_median_stratification,by="row.names")
colnames(survival_luad)[ncol(survival_luad)]<-"median_stratification"
survival_luad$median_stratification<-factor(survival_luad$median_stratification,levels=c("normoxia","hypoxia"))
res.cox <- coxph(Surv(c_event_time,c_event) ~ median_stratification, data=survival_luad)
summary(res.cox)
fit <- survfit(Surv(c_event_time,c_event)~ median_stratification, data=survival_luad)
par(pty="s")
gplotmeantro1<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantro1$plot <- gplotmeantro1$plot + labs(
  title    = "GSE3141 (Buffa signature)"        
)
gplotmeantro1 <- ggpar(
  gplotmeantro1,
  font.title    = c(25, "bold"))
gplotmeantro1$plot<-gplotmeantro1$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantro1$table <- ggpar(gplotmeantro1$table,
                            font.title = list(size = 16))



###GSE41271
GSE41271<-readRDS(file = "GSE41271.rds")
###processing the survival data for GSE42127
#write.table(GSE41271_survival,file="GSE41271_survival.csv",sep=",",row.names=FALSE)
survivalGSE41271<-read.csv(file="GSE41271_survival.csv",sep=",",header=T)
####censoring data
survivalGSE41271$status<-ifelse(survivalGSE41271$vital.statistics.ch1=="A",0,1)
censorship<-60
survivalGSE41271$OS   <- as.numeric( as.character(survivalGSE41271$OS ) )
GSE41271_survival_event_censorship  <- ifelse(survivalGSE41271$OS <= censorship &survivalGSE41271$status == 1 , 1 ,0  )
GSE41271_survival_time_censorship   <- ifelse( GSE41271_survival_event_censorship == 0 & survivalGSE41271$OS>=censorship ,censorship ,
survivalGSE41271$OS )   
survivalGSE41271   <- cbind(survivalGSE41271, GSE41271_survival_event_censorship, GSE41271_survival_time_censorship )
colnames(survivalGSE41271)[ (ncol(survivalGSE41271)-1):ncol(survivalGSE41271) ] <- c("c_event","c_event_time")

survivalGSE41271$rstatus<-ifelse(survivalGSE41271$recurrence.ch1=="N",0,1)
censorship<-60
survivalGSE41271$RFS   <- as.numeric( as.character(survivalGSE41271$RFS ) )
GSE41271_rsurvival_event_censorship  <- ifelse(survivalGSE41271$RFS <= censorship &survivalGSE41271$rstatus == 1 , 1 ,0  )
GSE41271_rsurvival_time_censorship   <- ifelse( GSE41271_rsurvival_event_censorship == 0 & survivalGSE41271$RFS>=censorship ,censorship ,
survivalGSE41271$RFS )   
survivalGSE41271   <- cbind(survivalGSE41271, GSE41271_rsurvival_event_censorship, GSE41271_rsurvival_time_censorship )
colnames(survivalGSE41271)[ (ncol(survivalGSE41271)-1):ncol(survivalGSE41271) ] <- c("c_revent","c_revent_time")
GSE41271_survival<-survivalGSE41271[survivalGSE41271$histology.ch1=="Adenocarcinoma",]
rownames(GSE41271_survival)<-GSE41271_survival$geo_accession
exprs<-GSE41271[,intersect(colnames(GSE41271),rownames(GSE41271_survival))]
exprs<-log2(exprs)
seed_genes[2]<-"UTP11"
seed_genes[16]<-"C6"
sig_exprs<-madSummaryAgilent_Buffa(illuminaHumanv3.db,seed_genes,exprs)
rownames(LUAD_genes)<-seed_genes
DataLUADsig<-sig_exprs[rownames(LUAD_genes),]
score<-apply(DataLUADsig,2,median)
###normalise to z-score

median_stratification<-ifelse(score>median(score),"hypoxia","normoxia")
luadsig_median_stratification<-as.data.frame(median_stratification)
####median stratification
survival_luad<-merge(GSE41271_survival,luadsig_median_stratification,by="row.names")
colnames(survival_luad)[ncol(survival_luad)]<-"median_stratification"
survival_luad$median_stratification<-factor(survival_luad$median_stratification,levels=c("normoxia","hypoxia"))

res.cox <- coxph(Surv(c_event_time,c_event)~median_stratification, data=survival_luad)
summary(res.cox)
fit <- survfit(Surv(c_event_time,c_event)~median_stratification, data=survival_luad)
par(pty="s")
gplotmeantro3<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantro3$plot <- gplotmeantro3$plot + labs(
  title    = "GSE41271 (Buffa signature)"        
)
gplotmeantro3 <- ggpar(
  gplotmeantro3,
  font.title    = c(25, "bold"))
gplotmeantro3$plot<-gplotmeantro3$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantro3$table <- ggpar(gplotmeantro3$table,
                            font.title = list(size = 16))

#######
                            
load( "GSE42127_LUAD_survival.RData")
load("GSE42127_exprs.RData")
#GSE42127_ad_survival<-GSE42127_ad_survival[which(GSE42127_ad_survival$had_adjuvant_chemo==T),]
exprs<-exprs[,intersect(colnames(exprs),rownames(GSE42127_ad_survival))]
sig_exprs<-madSummaryAgilent_Buffa(illuminaHumanv3.db,seed_genes,exprs)
DataLUADsig<-sig_exprs[rownames(LUAD_genes),]
score<-apply(DataLUADsig,2,median)
###normalise to z-score
median_stratification<-ifelse(score>median(score),"hypoxia","normoxia")
luadsig_median_stratification<-as.data.frame(median_stratification)
####median stratification
survival_luad<-merge(GSE42127_ad_survival,luadsig_median_stratification,by="row.names")
colnames(survival_luad)[ncol(survival_luad)]<-"median_stratification"
survival_luad$median_stratification<-factor(survival_luad$median_stratification,levels=c("normoxia","hypoxia"))

GSE42127_fit<-survfit(Surv(censored_os_time,censored_os_status)~median_stratification, data=survival_luad)
fit <- survfit(Surv(censored_os_time,censored_os_status)~ median_stratification, data=survival_luad)
par(pty="s")
gplotmeantro8<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantro8$plot <- gplotmeantro8$plot + labs(
  title    = "GSE42127 (Buffa signature)"        
)
gplotmeantro8 <- ggpar(
  gplotmeantro8,
  font.title    = c(25, "bold"))
gplotmeantro8$plot<-gplotmeantro8$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantro8$table <- ggpar(gplotmeantro8$table,
                            font.title = list(size = 16)) 

############


load( "GSE72094_LUAD_survival.RData"    )
load( "GSE72094_exprs.RData"    )

exprs<-exprs[,intersect(colnames(exprs),rownames(GSE72094_ad_survival))]
LUAD_genes[,1][10]<-"C20orf20"
LUAD_genes[,1][50]<-"CTSL2"
rownames(LUAD_genes)<-as.character(LUAD_genes[,1])
DataLUADsig<-exprs[rownames(LUAD_genes),]
score<-apply(DataLUADsig,2,median)
median_stratification<-ifelse(score>median(score),"hypoxia","normoxia")
luadsig_median_stratification<-as.data.frame(median_stratification)
####median stratification
survival_luad<-merge(GSE72094_ad_survival,luadsig_median_stratification,by="row.names")
colnames(survival_luad)[ncol(survival_luad)]<-"median_stratification"
survival_luad$median_stratification<-factor(survival_luad$median_stratification,levels=c("normoxia","hypoxia"))
survdiff(Surv(censored_time,censored_status)~median_stratification, data=survival_luad)
GSE72094_fit<-survfit(Surv(censored_time,censored_status)~median_stratification, data=survival_luad)
res.cox <- coxph(Surv(censored_time,censored_status)~median_stratification, data= survival_luad)
summary(res.cox)
fit <- survfit(Surv(censored_time,censored_status)~median_stratification, data= survival_luad)
par(pty="s")
gplotmeantro2<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.01),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantro2$plot <- gplotmeantro2$plot + labs(
  title    = "GSE72094 (Buffa signature)"        
)
gplotmeantro2 <- ggpar(
  gplotmeantro2,
  font.title    = c(25, "bold"))
gplotmeantro2$plot<-gplotmeantro2$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantro2$table <- ggpar(gplotmeantro2$table,
                            font.title = list(size = 16))

###MVA
GSE72094MVA<-readRDS(file="GSE72094.rds")
rownames(survival_luad)<-survival_luad$Row.names
GSE72094_MVA<-merge(GSE72094MVA,survival_luad,by="row.names")

res.cox <- coxph(Surv(censored_time,censored_status)~ GSE72094_MVA[,7], data=GSE72094_MVA)
summary(res.cox)


stage<-sapply(GSE72094_MVA[,6],function(x){
	if(x=="1"|x=="1A"|x=="1B"|x=="2A"|x=="2B"){y=1}
	else{y=2}})
GSE72094_MVA$stage<-as.factor(stage)
res.cox <- coxph(Surv(censored_time,censored_status)~stage, data=GSE72094_MVA)
summary(res.cox)

res.cox <- coxph(Surv(censored_time,censored_status)~GSE72094_MVA[,7]+stage+median_stratification, data=GSE72094_MVA)
summary(res.cox)














splots<-list(gplotmeantro2,gplotmeantro8,gplotmeantro3,gplotmeantro1)
s<-arrange_ggsurvplots(splots,ncol=2,nrow=2)
ggsave("/Users/mkhan/Documents/Brian RNA-seq scripts/Figure S_Buffa.pdf",plot=s,height=12,width=12)

#####



