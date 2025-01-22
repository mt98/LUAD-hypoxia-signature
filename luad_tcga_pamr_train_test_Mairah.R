####This script was used for the development of the TCGA training and test datasets..... 

module load apps/gcc/R/3.6.1 
qrsh -l short -V -cwd -pe smp.pe 8 R --vanilla --interactive
 
scp -r "/Users/mkhan/Documents/LUAD_Train_Exprs.RData" mqbpkmk3@csf3.itservices.manchester.ac.uk:~
scp -r "/Users/mkhan/Documents/Brian RNA-seq scripts/LUAD_Test_Exprs.RData" mqbpkmk3@csf3.itservices.manchester.ac.uk:~
scp -r "/Users/mkhan/Documents/LUAD_Train_Survival.RData" mqbpkmk3@csf3.itservices.manchester.ac.uk:~
scp -r "/Users/mkhan/Documents/Brian RNA-seq scripts/Adeno_Seed_Genes_Expin3CellLines.csv" mqbpkmk3@csf3.itservices.manchester.ac.uk:~


setwd("/mnt/iusers01/cw01/mqbpkmk3")
#### Load libraries
setwd("/Users/mkhan/Documents/PhD work/Brian RNA-seq scripts")
library(pheatmap)
library(hgu133plus2.db)
library( gplots )
library(survival)
library(survminer)
library(pamr)
library(org.Hs.eg.db)
library(limma)


#### Load expression data


 load("LUAD_Train_Exprs.RData")#train_exprs
 load("LUAD_Test_Exprs.RData")#test_exprs
 load("LUAD_Train_Survival.RData")#Train_surv
 load("LUAD_Test_Survival.RData")#Test_surv
 
 

centre_train_exprs<-sweep(train_exprs,1,apply(train_exprs,1,median))
centre_test_exprs<-sweep(test_exprs,1,apply(test_exprs,1,median))

##### Seed genes 

ad_seed_df<-read.csv(file="Adeno_Seed_Genes_Expin3CellLines.csv")#check numerical values are genuine; examine .txt file
ad_seed_genes<-as.character(ad_seed_df$Symbol)
ad_seed_title<-as.character(ad_seed_df$Title)
 
#### Match Expression and Survival data
#common_train<-intersect(colnames(centre_train_exprs),rownames(Train_surv) )
#common_test<-intersect(colnames(centre_test_exprs),rownames(Test_surv) )


train_exprs <- centre_train_exprs[,which( colnames(centre_train_exprs)%in%rownames(Train_surv) )]
train_surv<-Train_surv[colnames(train_exprs),]
select_train_exprs<-train_exprs[which(rownames(train_exprs)%in%ad_seed_genes),]

test_exprs <- centre_test_exprs[,which( colnames(centre_test_exprs)%in%rownames(Test_surv) )]
test_surv<-Test_surv[colnames(test_exprs),]
select_test_exprs<-test_exprs[which(rownames(test_exprs)%in%ad_seed_genes),]


set.seed(321)
kclass<-kmeans(t(select_train_exprs),2,iter.max=1000)
cl_train<-as.integer(kclass$cl)
class_train<-factor(ifelse(kclass$cl==1,"Nor","Hyp"))
hypoxia_state<-factor(ifelse(kclass$cl==1,"Nor","Hyp"),levels=c("Nor","Hyp"))
o_cl_train<-order(cl_train)

class_one_mean<-apply(select_train_exprs[,class_train=="Nor"],1,mean)
class_two_mean<-apply(select_train_exprs[,class_train=="Hyp"],1,mean)
cbind(class_one_mean, class_two_mean)

low_hypoxia_mean<-apply(select_train_exprs[,hypoxia_state=="Nor"],1,mean)
high_hypoxia_mean<-apply(select_train_exprs[,hypoxia_state=="Hyp"],1,mean)
cbind(low_hypoxia_mean, high_hypoxia_mean)

hypoxia_class<-ifelse(cl_train==0,"Hyp","Nor")

###################################################################
#Kmeans Heatmap with pheatmap
annotation_col= data.frame(group = class_train)
rownames(annotation_col)=names(class_train)
ann_colors = list(group= c(Hyp="blue", Nor="red"))
jpeg(file="heatmappamr.jpg")
pheatmap(select_train_exprs[,o_cl_train], fontsize=6, treeheight_row=0,color=colorRampPalette(c("navy", "white", "orange"))(50), 
annotation_colors=ann_colors,cluster_rows=FALSE, cluster_cols=FALSE, clustering_distance_rows="correlation",annotation_col=annotation_col,labels_col=rep("",ncol(select_train_exprs)),scale="row")
dev.off()
scp -r mqbpkmk3@csf3.itservices.manchester.ac.uk:~/heatmappamr.jpg "/Users/mkhan/Documents/Brian RNA-seq scripts"
###EdgeR between the two k-means clusters group and GSEA etc...
#look at colorpanel{gplots} help page for options
set.seed(120)

#train_data<- list(x = select_train_exprs, y = factor(cl_train), genenames = ad_seed_title, geneid = ad_seed_genes)
train_data<- list(x = select_train_exprs, y = hypoxia_state, genenames = ad_seed_title, geneid = ad_seed_genes)
#train_data<- list(x = select_train_exprs, y = as.factor(hypoxia_class), genenames = ad_seed_title, geneid = ad_seed_genes)

train_model<-pamr.train(train_data)
model_cv <- pamr.cv(train_model, train_data,nfold = 10)
model_cv

best_index<-max(which(model_cv$error==min(model_cv$error)))
train_threshold<-model_cv$threshold[best_index]
pamr.confusion(model_cv,train_threshold)

yhat_train<-train_model$yhat[,best_index]
names(yhat_train)<-colnames(select_train_exprs)

survdiff(Surv(censored_time,censored_status)~yhat_train, data=Train_surv)
train_fit<-survfit(Surv(censored_time,censored_status)~yhat_train, data=Train_surv)

fit <- survfit(Surv(censored_time,censored_status)~ yhat_train, data= Train_surv)
par(pty="s")
gplotmeantro<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.0),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantro$plot <- gplotmeantro$plot + labs(

  title    = "TCGA training"        
)
gplotmeantro <- ggpar(
  gplotmeantro,
  font.title    = c(25, "bold"))
gplotmeantro$plot<-gplotmeantro$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantro$table <- ggpar(gplotmeantro$table,
                            font.title = list(size = 16))
###look at the distribution of the clinicopathological variables

train_mva<-merge(Train_surv,yhat_train,by="row.names")

train_mva_hyp<-train_mva[train_mva$y=="Hyp",]
train_mva_nor<-train_mva[train_mva$y=="Nor",]




yhat_test <- pamr.predict(train_model,select_test_exprs, train_threshold)
names(yhat_test)<-colnames(select_test_exprs)
	
survdiff(Surv(censored_time,censored_status)~yhat_test, data=Test_surv)
test_fit<-survfit(Surv(censored_time,censored_status)~yhat_test, data=Test_surv)

res.cox <- coxph(Surv(censored_time,censored_status)~yhat_test, data=Test_surv )
summary(res.cox)
fit <- survfit(Surv(censored_time,censored_status)~ yhat_test, data= Test_surv)
par(pty="s")
gplotmeantro1<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.0),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",22), font.y = c("bold",22), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 7.0,risk.table.height = 0.2,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantro1$plot <- gplotmeantro1$plot + labs(
  title    = "TCGA test"        
)
gplotmeantro1 <- ggpar(
  gplotmeantro1,
  font.title    = c(25, "bold"))
gplotmeantro1$plot<-gplotmeantro1$plot + theme(plot.title = element_text(hjust = 0.5))

gplotmeantro1$table <- ggpar(gplotmeantro1$table,
                            font.title = list(size = 16))
splots<-list(gplotmeantro,gplotmeantro1)
s<-arrange_ggsurvplots(splots,ncol=2,nrow=1)
###look at the distribution of the clinicopathological variables

###########################################
#Is hypoxia_state class list found using KNN###
#associated with enrichment of hypoxia related# 
#functional processes?#########################
#Create limma DEG model fit ###################
#Use to generate log FC for DEG ###############
#Then apply fgsea##############################
###############################################

scp -r "/Users/mkhan/Documents/Brian RNA-seq scripts/human_H_v5p2.Rdata" mqbpkmk3@csf3.itservices.manchester.ac.uk:~
setwd("/Users/brianlane/Projects/Lung Cancer RNAseq/data")
load("human_H_v5p2.Rdata")#Hs.H Broad MolSigDB Hallmark Human


#design<-model.matrix(~class_train)
design<-model.matrix(~hypoxia_state)
head(design)

fit<-lmFit(train_exprs,design)
efit<-eBayes(fit)



deg_table<-topTable(efit,coef=2, adjust.method="BH",number=nrow(train_exprs))
deg_table<-deg_table[deg_table$adj.P.Val<0.05,]
deg_fc<-deg_table$logFC
all_deg_genes<-rownames(deg_table)


entrez_de<-as.character(mget(all_deg_genes,org.Hs.egSYMBOL2EG,ifnotfound=NA))
na_entrez<-which(entrez_de=='NA')
all_entrez<-entrez_de[-na_entrez]
deg_fc<-deg_fc[-na_entrez]
names(deg_fc)<-all_entrez

#range(lapply(pathwaysH,length))

gseaRes<-fgsea(Hs.H,deg_fc,minSize=20,maxSize=200, nperm=1000)

topPathwaysUp <- gseaRes[gseaRes$ES>0&gseaRes$padj<0.05,]
write.table(topPathwaysUp[,1:5],file="/mnt/iusers01/cw01/mqbpkmk3/GSEA.csv",sep=",")


scp -r mqbpkmk3@csf3.itservices.manchester.ac.uk:~/GSEA.csv "/Users/mkhan/Documents/Brian RNA-seq scripts"



topPathwaysUp <- gseaRes[ES > 0][head(order(NES,decreasing=T), n=20), pathway]
topPathwaysDown <- gseaRes[ES < 0][head(order(NES,decreasing=F), n=20), pathway]

topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

#plot top ranked up- and down-regulated processes
#png("Train_GSEPAlot.png", width = 8*300,  height =3* 500, res = 300, pointsize = 8) 
plotGseaTable(Hs.H[topPathways], deg_fc, gseaRes, gseaParam = 0.5)
#dev.off()

#Plot Hallmark Hypoxia enrichment
#fgsea(pathwaysH[["HALLMARK_HYPOXIA"]],all_fc,minSize=20,maxSize=200,nperm=1000)
plotEnrichment(Hs.H[["HALLMARK_HYPOXIA"]],deg_fc)


####################################
deg_table<-topTable(efit,coef=2, adjust.method="BH",number=nrow(train_exprs))
deg_table_up<-deg_table[deg_table$adj.P.Val<0.05&deg_table$logFC>0.58,]
seeds<-deg_table_up[which(ad_seed_genes%in%rownames(deg_table_up)),]
seeds_normal<-deg_table[which(ad_seed_genes%in%rownames(deg_table)),]
write.table(deg_table_up,file="/mnt/iusers01/cw01/mqbpkmk3/degGENES.csv",sep=",")
deg_table_GSEA<-deg_table[deg_table$adj.P.Val<0.01&deg_table$logFC>0.58,]
write.table(deg_table_GSEA,file="/mnt/iusers01/cw01/mqbpkmk3/degGENES_FDR0.01.csv",sep=",")


scp -r mqbpkmk3@csf3.itservices.manchester.ac.uk:~/degGENES.csv "/Users/mkhan/Documents/Brian RNA-seq scripts"

scp -r mqbpkmk3@csf3.itservices.manchester.ac.uk:~/degGENES_FDR0.01.csv "/Users/mkhan/Documents/Brian RNA-seq scripts"


deg_fc<-deg_table$logFC
all_deg_genes<-rownames(deg_table)


entrez_de<-as.character(mget(all_deg_genes,org.Hs.egSYMBOL2EG,ifnotfound=NA))
na_entrez<-which(entrez_de=='NA')
all_entrez<-entrez_de[-na_entrez]
deg_fc<-deg_fc[-na_entrez]
names(deg_fc)<-all_entrez

#range(lapply(pathwaysH,length))

gseaRes<-fgsea(Hs.H,deg_fc,minSize=20,maxSize=200, nperm=1000)


topPathwaysUp <- gseaRes[ES > 0][head(order(NES,decreasing=T), n=20), pathway]
topPathwaysDown <- gseaRes[ES < 0][head(order(NES,decreasing=F), n=20), pathway]

topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

#plot top ranked up- and down-regulated processes
#png("Train_GSEPAlot.png", width = 8*300,  height =3* 500, res = 300, pointsize = 8) 
plotGseaTable(Hs.H[topPathways], deg_fc, gseaRes, gseaParam = 0.5)
#dev.off()

#Plot Hallmark Hypoxia enrichment
#fgsea(pathwaysH[["HALLMARK_HYPOXIA"]],all_fc,minSize=20,maxSize=200,nperm=1000)
plotEnrichment(Hs.H[["HALLMARK_HYPOXIA"]],deg_fc)


### Make a k-means plot in R at my own desktop


library( gplots )
library(survival)
library(survminer)
library(pamr)

#### Load expression data
load("/Users/mkhan/Documents/PhD work/Brian RNA-seq scripts/LUAD_Train_Exprs.RData")#train_exprs
load("/Users/mkhan/Documents/PhD work/Brian RNA-seq scripts/LUAD_Test_Exprs.RData")#test_exprs
load("/Users/mkhan/Documents/PhD work/Brian RNA-seq scripts/LUAD_Train_Survival.RData")#Train_surv
load("/Users/mkhan/Documents/PhD work/Brian RNA-seq scripts/LUAD_Test_Survival.RData")#Test_surv

centre_train_exprs<-sweep(train_exprs,1,apply(train_exprs,1,median))
centre_test_exprs<-sweep(test_exprs,1,apply(test_exprs,1,median))



##### Seed genes 

ad_seed_df<-read.csv(file="/Users/mkhan/Documents/PhD work/Brian RNA-seq scripts/Adeno_Seed_Genes_Expin3CellLines.csv")#check numerical values are genuine; examine .txt file
ad_seed_genes<-as.character(ad_seed_df$Symbol)
ad_seed_title<-as.character(ad_seed_df$Title)
 
#### Match Expression and Survival data
#common_train<-intersect(colnames(centre_train_exprs),rownames(Train_surv) )
#common_test<-intersect(colnames(centre_test_exprs),rownames(Test_surv) )


train_exprs <- centre_train_exprs[,which( colnames(centre_train_exprs)%in%rownames(Train_surv) )]
train_surv<-Train_surv[colnames(train_exprs),]
select_train_exprs<-train_exprs[which(rownames(train_exprs)%in%ad_seed_genes),]

test_exprs <- centre_test_exprs[,which( colnames(centre_test_exprs)%in%rownames(Test_surv) )]
test_surv<-Test_surv[colnames(test_exprs),]
select_test_exprs<-test_exprs[which(rownames(test_exprs)%in%ad_seed_genes),]



set.seed(321)
kclass<-kmeans(t(select_train_exprs),2,iter.max=1000)
cl_train<-as.integer(kclass$cl)
class_train<-factor(ifelse(kclass$cl==1,"Nor","Hyp"))
hypoxia_state<-factor(ifelse(kclass$cl==1,"Nor","Hyp"),levels=c("Nor","Hyp"))
o_cl_train<-order(cl_train)

class_one_mean<-apply(select_train_exprs[,class_train=="Nor"],1,mean)
class_two_mean<-apply(select_train_exprs[,class_train=="Hyp"],1,mean)
cbind(class_one_mean, class_two_mean)

low_hypoxia_mean<-apply(select_train_exprs[,hypoxia_state=="Nor"],1,mean)
high_hypoxia_mean<-apply(select_train_exprs[,hypoxia_state=="Hyp"],1,mean)
cbind(low_hypoxia_mean, high_hypoxia_mean)

hypoxia_class<-ifelse(cl_train==1,"Nor","Hyp")

###################################################################
#Kmeans Heatmap with pheatmap
annotation_col= data.frame(Kmeans_Class = class_train)
rownames(annotation_col)=names(class_train)
ann_colors = list(Kmeans_class = c(Nor="blue", Hyp="yellow"))

pdf(file="/mnt/iusers01/cw01/mqbpkmk3/heatmappamr.pdf")
pheatmap(select_train_exprs[,o_cl_train], fontsize=6, treeheight_row=0,color=(greenred(50)), annotation_colors=ann_colors,cluster_rows=TRUE, cluster_cols=FALSE, clustering_distance_rows="correlation",annotation_col=annotation_col,labels_col=rep("",ncol(select_train_exprs)),scale="row")
#look at colorpanel{gplots} help page for options
dev.off()
scp -r mqbpkmk3@csf3.itservices.manchester.ac.uk:~/heatmappamr.pdf "/Users/mkhan/Documents/Brian RNA-seq scripts"
###EdgeR between the two k-means clusters group and GSEA etc...

set.seed(120)

#train_data<- list(x = select_train_exprs, y = factor(cl_train), genenames = ad_seed_title, geneid = ad_seed_genes)
train_data<- list(x = select_train_exprs, y = hypoxia_state, genenames = ad_seed_title, geneid = ad_seed_genes)
#train_data<- list(x = select_train_exprs, y = as.factor(hypoxia_class), genenames = ad_seed_title, geneid = ad_seed_genes)

train_model<-pamr.train(train_data)
model_cv <- pamr.cv(train_model, train_data,nfold = 10)
best_index<-max(which(model_cv$error==min(model_cv$error)))
train_threshold<-model_cv$threshold[best_index]
pamr.confusion(model_cv,train_threshold)

###work done to get MVA for TCGA-LUAD train and test data
(load("/Users/mkhan/Documents/PhD work/fifth year/LUAD_TCGA_Train_Test_Survival (1).RData"))
clinicaldata<-read.csv("/Users/mkhan/Documents/PhD work/fifth year/luad_tcga_clinical_data.tsv",sep="\t",header=T)
rownames(clinicaldata)<-clinicaldata$Patient.ID
detail_train<-merge(clinicaldata[intersect(rownames(clinicaldata),rownames(Train_surv)),],Train_surv[intersect(rownames(clinicaldata),rownames(Train_surv)),],by="row.names")
rownames(detail_train)<-detail_train$Row.names
clinicaldata_cig<-read.csv( "/Users/mkhan/Documents/PhD work/fifth year/nsclc_tcga.csv",sep=",",header=T)
rownames(clinicaldata_cig)<-clinicaldata_cig$Patient.ID
detail_train_cig<-cbind(detail_train[intersect(rownames(detail_train),rownames(clinicaldata_cig)),],clinicaldata_cig[intersect(rownames(detail_train),rownames(clinicaldata_cig)),])

Training_dataset<-merge(Train_surv,yhat_train,by="row.names")
####make the km curves
res.cox <- coxph(Surv(censored_time,censored_status) ~ y, data =Training_dataset)
summary(res.cox)
fit<-survfit(Surv(censored_time,censored_status) ~ y, data = Training_dataset)
par(pty="s")
gplotmeantro<-ggsurvplot(fit,size=2,censor=TRUE,pval="P=0.0011",legend.title="",linetype=c(3,1),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.1),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",25), font.y = c("bold",25), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 5.0,risk.table.height = 0.3,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmeantro$plot <- gplotmeantro$plot + labs(
  title    = "TCGA training"       
)
gplotmeantro <- ggpar(
  gplotmeantro,
  font.title    = c(25, "bold"))
gplotmeantro$plot<-gplotmeantro$plot + theme(plot.title = element_text(hjust = 0.5))
gplotmeantro$table <- ggpar(gplotmeantro$table,
                            font.title = list(size = 16))
###multivariate analysis for it!!

Training_dataset<-merge(Training_dataset,detail_train[,"Diagnosis.Age"],by="row.names")
colnames(Training_dataset)[ncol(Training_dataset)]<-"Diagnosis.Age"
colnames(Training_dataset)[ncol(Training_dataset)-1]<-"y"

Training_dataset1<-Training_dataset[!is.na(Training_dataset$residual_tumor),]
Training_dataset1$residual_tumorstage<-ifelse(Training_dataset1$residual_tumor=="r0"|Training_dataset1$residual_tumor=="rx","low","high")

res.cox <- coxph(Surv(censored_time,censored_status) ~ residual_tumorstage, data =Training_dataset1)
summary(res.cox)


res.cox <- coxph(Surv(censored_time,censored_status) ~ Diagnosis.Age, data =Training_dataset)
summary(res.cox)

res.cox <- coxph(Surv(censored_time,censored_status) ~ gender, data =Training_dataset)
summary(res.cox)
###Look at differences in age and residual tumor
Training_dataset_hyp<-Training_dataset[Training_dataset$y=="High Hypoxia",]
Training_dataset_nor<-Training_dataset[Training_dataset$y=="Low Hypoxia",]

t.test(Training_dataset_hyp$Diagnosis.Age,Training_dataset_nor$Diagnosis.Age)
chisq.test(rbind(table(Training_dataset_hyp$gender),table(Training_dataset_nor$gender))


Training_dataset1_hyp<-Training_dataset1[Training_dataset1$y=="High Hypoxia",]
Training_dataset1_nor<-Training_dataset1[Training_dataset1$y=="Low Hypoxia",]
chisq.test(rbind(table(Training_dataset1_hyp$residual_tumorstage),table(Training_dataset1_nor$residual_tumorstage)))



Training_dataset2<-Training_dataset[!is.na(Training_dataset$pathologic_stage),]
Training_dataset2$pathologic_stage<-recode(Training_dataset2$pathologic_stage, 'stage i'=0, 'stage ia'=1, 'stage ib'=2, 
                                          'stage ii'=3,'stage iia'=4, 'stage iib'=5,'stage iiia'=6, 'stage iiib'=7,'stage iv'=8)
stagecategory<-sapply(Training_dataset2$pathologic_stage,function(x){
  if(x>5){y=1}
else{y=0}})
Training_dataset2$stagecategory<-stagecategory
res.cox <- coxph(Surv(censored_time,censored_status) ~ stagecategory, data =Training_dataset2)
summary(res.cox)


Training_dataset2_hyp<-Training_dataset2[Training_dataset2$y=="High Hypoxia",]
Training_dataset2_nor<-Training_dataset2[Training_dataset2$y=="Low Hypoxia",]
chisq.test(rbind(table(Training_dataset2_hyp$stagecategory),table(Training_dataset2_nor$stagecategory)))




#####look at the smoking status
Training_dataset3<-merge(Training_dataset,clinicaldata_cig[,"Smoking.History"],by="row.names")
colnames(Training_dataset3)[ncol(Training_dataset3)]<-"Smoking.History"

Training_dataset3<-Training_dataset3[!is.na(Training_dataset3$Smoking.History),]
Training_dataset3$Smokingcategory<-ifelse(Training_dataset3$Smoking.History=="Lifelong Non-Smoker","No","Yes")
res.cox <- coxph(Surv(censored_time,censored_status) ~ Smokingcategory, data =Training_dataset3)
summary(res.cox)

###put in multivariate

Training_dataset4<-Training_dataset1[!is.na(Training_dataset1$pathologic_stage),]
Training_dataset4$pathologic_stage<-recode(Training_dataset4$pathologic_stage, 'stage i'=0, 'stage ia'=1, 'stage ib'=2, 
                                          'stage ii'=3,'stage iia'=4, 'stage iib'=5,'stage iiia'=6, 'stage iiib'=7,'stage iv'=8)
stagecategory<-sapply(Training_dataset4$pathologic_stage,function(x){
  if(x>5){y=1}
else{y=0}})
Training_dataset4$stagecategory<-stagecategory


res.cox <- coxph(Surv(censored_time,censored_status) ~ stagecategory+y+residual_tumorstage, data =Training_dataset4)
summary(res.cox)



####TCGA test 
detail_test<-merge(clinicaldata,Test_surv,by="row.names")
rownames(detail_test)<-detail_test$Row.names
detail_test_cig<-cbind(detail_test[intersect(rownames(detail_test),rownames(clinicaldata_cig)),],clinicaldata_cig[intersect(rownames(detail_test),rownames(clinicaldata_cig)),])

Test_dataset<-merge(Test_surv,yhat_test,by="row.names")
rownames(Test_dataset)<-Test_dataset$Row.names
####make the km curves
res.cox <- coxph(Surv(censored_time,censored_status) ~ y, data =Test_dataset)
summary(res.cox)
fit<-survfit(Surv(censored_time,censored_status) ~ y, data = Test_dataset)
par(pty="s")
gplotmean<-ggsurvplot(fit,size=2,censor=TRUE,pval="P=0.0016",legend.title="",linetype=c(3,1),legend.labs=c("Normoxia","Hypoxia"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=7,pval.coord=c(0,0.1),
                         risk.table = TRUE,palette=c("blue","red"),ylab="Overall survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",25), font.y = c("bold",25), font.tickslab = "bold",
                         font.legend= c( "bold", "black",16),risk.table.fontsize = 5.0,risk.table.height = 0.3,risk.table.col = "strata",risk.table.y.text = FALSE)

gplotmean$plot <- gplotmean$plot + labs(
  title    = "TCGA test"       
)
gplotmean <- ggpar(
  gplotmean,
  font.title    = c(25, "bold"))
gplotmean$plot<-gplotmean$plot + theme(plot.title = element_text(hjust = 0.5))
gplotmean$table <- ggpar(gplotmean$table,
                            font.title = list(size = 16))
splots<-list(gplotmeantro,gplotmean)
surv<-arrange_ggsurvplots(splots, print = TRUE,ncol=2,nrow=1)
###univariate analysis
res.cox <- coxph(Surv(censored_time,censored_status) ~ gender, data =Test_dataset)
summary(res.cox)

####surgical margin status

Test_dataset1<-Test_dataset[!is.na(Test_dataset$pathologic_stage),]
Test_dataset1$pathologic_stage<-recode(Test_dataset1$pathologic_stage, 'stage i'=0, 'stage ia'=1, 'stage ib'=2, 
                                          'stage ii'=3,'stage iia'=4, 'stage iib'=5,'stage iiia'=6, 'stage iiib'=7,'stage iv'=8)
stagecategory_test<-sapply(Test_dataset1$pathologic_stage,function(x){
  if(x>5){y=1}
  else{y=0}})
Test_dataset1$stagecategory_test<-stagecategory_test
res.cox <- coxph(Surv(censored_time,censored_status) ~ stagecategory_test, data =Test_dataset1)
summary(res.cox)
###look at the stage category
Test_dataset1_hyp<-Test_dataset1[Test_dataset1$y=="High Hypoxia",]
Test_dataset1_nor<-Test_dataset1[Test_dataset1$y=="Low Hypoxia",]
chisq.test(rbind(table(Test_dataset1_hyp$stagecategory_test),table(Test_dataset1_nor$stagecategory_test)))





###s
Test_dataset2<-Test_dataset[!is.na(Test_dataset$residual_tumor),]
Test_dataset2$residual_tumorstage<-ifelse(Test_dataset2$residual_tumor=="r0"|Test_dataset2$residual_tumor=="rx","low","high")

res.cox <- coxph(Surv(censored_time,censored_status) ~ residual_tumorstage, data =Test_dataset2)
summary(res.cox)

Test_dataset2_hyp<-Test_dataset2[Test_dataset2$y=="High Hypoxia",]
Test_dataset2_nor<-Test_dataset2[Test_dataset2$y=="Low Hypoxia",]
chisq.test(rbind(table(Test_dataset2_hyp$residual_tumorstage),table(Test_dataset2_nor$residual_tumorstage)))



Test_dataset3<-merge(yhat_test,detail_test,by="row.names")

res.cox <- coxph(Surv(censored_time,censored_status) ~ Diagnosis.Age, data =Test_dataset3)
summary(res.cox)




Test_dataset4<-merge(Test_dataset[intersect(rownames(Test_dataset),rownames(clinicaldata_cig)),],clinicaldata_cig[intersect(rownames(Test_dataset),rownames(clinicaldata_cig)),],by="row.names")
Test_dataset4<-Test_dataset4[!is.na(Test_dataset4$Smoking.History),]
Test_dataset4$Smokingcategory<-ifelse(Test_dataset4$Smoking.History=="Lifelong Non-Smoker","No","Yes")
res.cox <- coxph(Surv(censored_time,censored_status) ~ Smokingcategory, data =Test_dataset4)
summary(res.cox)

###multivariate

Test_dataset5<-Test_dataset1[!is.na(Test_dataset1$residual_tumor),]
Test_dataset5$residual_tumorstage<-ifelse(Test_dataset5$residual_tumor=="r0"|Test_dataset5$residual_tumor=="rx","low","high")

res.cox <- coxph(Surv(censored_time,censored_status) ~ stagecategory_test+y+residual_tumorstage, data =Test_dataset5)
summary(res.cox)

# doing comparisons between train and test set for p values ---------------
###stages
Test_dataset<-merge(Test_surv,yhat_test,by="row.names")
Test_dataset<-Test_dataset[!is.na(Test_dataset$pathologic_stage),]
Test_dataset$pathologic_stage<-recode(Test_dataset$pathologic_stage, 'stage i'=0, 'stage ia'=0, 'stage ib'=0, 
                                      'stage ii'=1,'stage iia'=1, 'stage iib'=1,'stage iiia'=2, 'stage iiib'=2,'stage iv'=3)

Training_dataset<-merge(Train_surv,yhat_train,by="row.names")
Training_dataset<-Training_dataset[!is.na(Training_dataset$pathologic_stage),]
Training_dataset$pathologic_stage<-recode(Training_dataset$pathologic_stage, 'stage i'=0, 'stage ia'=0, 'stage ib'=0, 
                                          'stage ii'=1,'stage iia'=1, 'stage iib'=1,'stage iiia'=2, 'stage iiib'=2,'stage iv'=3)

####smoking category
table(detail_test_cig$Smoking.History)
smoking_category_test<-ifelse(detail_test_cig$Smoking.History=="Current Smoker","Smoking","Non-smoking","")

table(detail_train_cig$Smoking.History)
smoking_category_train<-ifelse(detail_train_cig$Smoking.History=="Current Smoker","Smoking","Non-smoking")










