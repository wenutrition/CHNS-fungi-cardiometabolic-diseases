library("ggplot2")
library("corrplot")
library(ggrepel)
library(permute)
library(vegan)
library(plyr)
library (zCompositions)
library(microbiome)
library(phyloseq)
library(dplyr)
library(NetCoMi)
library(haven)
library(mediation)

#---------------PCOA analysis------------------

df<-read_dta("D:\\真菌分析\\data_pre\\2015_all_fungi_pcoa.dta")
otu<- as.data.frame(df )
phenotype<-read_dta("D:\\真菌分析\\data_pre\\2015_fungi_phenotype_pcoa.dta")
phenotype$veg<-phenotype$dveg+phenotype$lveg
phenotype<- as.data.frame(phenotype)
phenotype$t2d <- as.factor(phenotype$t2d)
phenotype$city <- as.factor(phenotype$city)
phenotype$education <- as.factor(phenotype$education)
phenotype$marrige <- as.factor(phenotype$marrige)
phenotype$smoke <- as.factor(phenotype$smoke)
phenotype$alcohol <- as.factor(phenotype$alcohol)
phenotype$hyper_med <- as.factor(phenotype$hyper_med)
phenotype$xinjigengse <- as.factor(phenotype$xinjigengse)
phenotype$zhongfeng <- as.factor(phenotype$zhongfeng)
phenotype$cancer <- as.factor(phenotype$cancer)
phenotype$sex <- as.factor(phenotype$sex)
phenotype$diabetes_med <- as.factor(phenotype$diabetes_med)
phenotype$hyp <- as.factor(phenotype$hyp)
phenotype$dys <- as.factor(phenotype$dys)
phenotype$pret2d <- as.factor(phenotype$pret2d)
phenotype$famine <- as.factor(phenotype$famine)
phenotype$gut_disease <- as.factor(phenotype$gut_disease)
phenotype$fuxie <- as.factor(phenotype$fuxie)
phenotype$antibiotic_current <- as.factor(phenotype$antibiotic_current)
phenotype$antibiotic_6month <- as.factor(phenotype$antibiotic_6month)
phenotype$probiotics <- as.factor(phenotype$probiotics)
phenotype$kangyan_med <- as.factor(phenotype$kangyan_med)
phenotype$kangsuan_med <- as.factor(phenotype$kangsuan_med)
phenotype$weisuan_med <- as.factor(phenotype$weisuan_med)
phenotype$yogurt_drink <- as.factor(phenotype$yogurt_drink)
phenotype$changdao_shoushu <- as.factor(phenotype$changdao_shoushu)
phenotype$pet <- as.factor(phenotype$pet)
#phenotype<-phenotype[,10:73]

set.seed(666)
fungi_distance = as.matrix(vegdist(otu, method="bray"))

results=NULL
for(i in colnames(phenotype)){
  adonis=adonis(fungi_distance~phenotype[,i], permutations = 1000)
  results=rbind(results,adonis$aov.tab[1,])
}
write.table(results,"D:\\真菌分析\\Final results\\All data\\fungi_permernoval.csv",sep=",")

#---------------------Network analysis-----------

dfm<- readxl::read_xlsx("D:\\真菌分析\\rawdata\\Repeat data.xlsx",sheet = "15") 
dfm<-as.data.frame(dfm)
rownames(dfm1)<-dfm$SampleID
otu<- as.data.frame(subset( dfm, select = -c(SampleID) ))
sparcc <- netConstruct(otu,  
                           measure = "sparcc",alpha = 0.2,
                           verbose = 3,
                           seed = 123456)
#Network analysis

props_single <- netAnalyze(sparcc, 
                           centrLCC = FALSE,
                           avDissIgnoreInf = TRUE,
                           sPathNorm = FALSE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = c("degree", "between", "closeness"),
                           hubQuant = 0.9,
                           lnormFit = TRUE,
                           normDeg = FALSE,
                           normBetw = FALSE,
                           normClose = FALSE,
                           normEigen = FALSE)
summary(props_single)

#Visual network comparison

plot(props_single, repulsion=0.9,rmSingles="all",esize=8,	posCol=c("#F08080"), negCol=c("#1E90FF"),  cexLabels = 2, cexHubLabels = 1.5,
     nodeSizeSpread =8)

#--------------------------DMM analysis---------------
dat <- abundances(physeq)
count <- as.matrix(t(dat))
fit <- lapply(1:5, dmn, count = count, verbose=TRUE)

lplc <- base::sapply(fit, DirichletMultinomial::laplace) # AIC / BIC / Laplace
aic  <- base::sapply(fit, DirichletMultinomial::AIC) # AIC / BIC / Laplace
bic  <- base::sapply(fit, DirichletMultinomial::BIC) # AIC / BIC / Laplace
plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
lines(aic, type="b", lty = 2)
lines(bic, type="b", lty = 3)
best <- fit[[which.min(unlist(aic))]]
best
mixturewt(best)
ass <- apply(mixture(best), 1, which.max)
ass
OUT<-as.data.frame(ass)
write.table(ass,"D:\\真菌分析\\Final results\\All data\\DMM\\All_cluster.csv",sep=",")
for (k in seq(ncol(fitted(best)))) {
  d <- melt(fitted(best))
  colnames(d) <- c("OTU", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    # Arrange OTUs by assignment strength
    arrange(value) %>%
    mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
    # Only show the most important drivers
    filter(abs(value) > quantile(abs(value), 0.99))     
  
  p <- ggplot(d, aes(x = OTU, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top drivers: community type", k))
  print(p)
}

#---------------------------Mediation analysis-------------

#--step1: metabolism module analyses--

df_my<- readxl::read_xlsx("D:\\真菌分析\\Final results\\All data\\metabolism\\mediation_input.xlsx",sheet = "module_diseases") 
df_my<-subset(df_my, outcome!="dys") 
df_my<-subset(df_my, metabolism!="megrey") 
p_adjust = p.adjust(df_my$p,method="BH")
df_my = as.data.frame(cbind(df_my, p_adjust))

df_xm<- readxl::read_xlsx("D:\\真菌分析\\Final results\\All data\\metabolism\\mediation_input.xlsx",sheet = "fungi_module")
df_xm<-subset(df_xm, metabolism!="megrey") 
p_adjust = p.adjust(df_xm$p,method="BH")
df_xm = as.data.frame(cbind(df_xm, p_adjust))

dataset_names <- list('module_disease' = df_my, 'fungi_module' = df_xm)
openxlsx::write.xlsx(dataset_names,rowNames=TRUE, file = 'D:\\真菌分析\\Plot\\data\\代谢模块分析.xlsx') 

df_my<-subset(df_my, p_adjust<0.05 ) 
colnames(df_my)[3:8]=paste(colnames(df_my)[3:8],"_moduledisease",sep = "") # rename the columns for merge
df_xm<-subset(df_xm, p_adjust<0.05 ) 
result=merge(df_xm[,c(1,3:7)],df_my[,c(1,2,4:8)],by="metabolism",all = F)

metabolism<-read_dta("D:\\真菌分析\\data_pre\\metabolism_data.dta")
rownames(metabolism)=metabolism$SampleID
metabolism<-subset(metabolism,select=-c(SampleID))

fungi<-read_dta("D:\\真菌分析\\data_pre\\outcome_data.dta")
rownames(fungi)=fungi$SampleID
fungi<-subset(fungi,select=-c(SampleID))

diseases<-read_dta("D:\\真菌分析\\data_pre\\diseases_data.dta")
rownames(diseases)=diseases$SampleID
diseases<-subset(diseases,select=-c(SampleID))

cov<-read_dta("D:\\真菌分析\\data_pre\\cov_data.dta")
rownames(cov)=cov$SampleID
cov<-subset(cov,select=-c(SampleID))
cov$sex <- as.factor(cov$sex)
cov$district<-as.factor(cov$district)
set.seed(1234)  # 设置随机种子，确保每次运行结果一致

for(i in 1:nrow(result)){
  data=data.frame(cbind(fungi,diseases[,as.character(result$outcome[i])],metabolism[,as.character(result$metabolism[i])],cov))
  data=na.omit(data)
  colnames(data)=c("X","Y","M","age","sex","BMI","district")
  model.m=lmer(M ~ X + age + sex + BMI + (1 | district), data = data)
  model.y=glmer(Y ~ X + M + age + sex + BMI + (1 | district), data = data, family = "binomial")
  summary=summary(mediate(model.m ,model.y,treat = "X", mediator = "M",boot =FALSE))
  result$Pval_mediate[i]=summary$d.avg.p
  result$coef_mediate[i]=summary$d.avg
  result$coefCI_mediate[i]=summary$d.avg.ci
  result$Pval_direct[i]=summary$z.avg.p
  result$coef_direct[i]=summary$z.avg
  result$coefCI_direct[i]=summary$z.avg.ci
  result$Pval_total[i]=summary$tau.p
  result$coef_total[i]=summary$tau.coef
  result$coefCI_total[i]=summary$tau.ci
  result$Pval_ratio[i]=summary$n.avg.p
  result$coef_ratio[i]=summary$n.avg
  result$coefCI_ratio[i]=summary$n.avg.ci
  #inverse mediate
  colnames(data)=c("X","M","Y","age","sex","BMI","district")
  model.m=glmer(M~X+age+sex+BMI+ (1 | district), data = data, family = "binomial")
  model.y=lmer(Y~X+M+age+sex+BMI+ (1 | district), data = data)
  summary=summary(mediate(model.m ,model.y,treat = "X", mediator = "M",boot =FALSE ))
  result$Pval_mediate_inverse[i]=summary$d.avg.p
  result$Pval_direct_inverse[i]=summary$z.avg.p
}

result$Qval_mediate=p.adjust(result$Pval_mediate,method = "BH")
result$Qval_mediate_inverse=p.adjust(result$Pval_mediate_inverse,method = "BH")
write.table(result,file = "D:\\真菌分析\\Final results\\All data\\代谢模块中介分析结果.csv",quote = F,sep = "\t",row.names = F)

#--step2: individual metabolisms analyses--

df_my<- readxl::read_xlsx("D:\\真菌分析\\Final results\\All data\\metabolism\\mediation_input.xlsx",sheet = "metabolism_diseases_mpping") 
df_my<-subset(df_my, outcome!="dys") 
df_my<-subset(df_my, outcome!="hyp") 
df_my<-subset(df_my, color!="brown") # brown完全遮掩，故摒弃
p_adjust = p.adjust(df_my$p,method="BH")
df_my = as.data.frame(cbind(df_my, p_adjust))


df_xm<- readxl::read_xlsx("D:\\真菌分析\\Final results\\All data\\metabolism\\mediation_input.xlsx",sheet = "fungi_metabolism")
p_adjust = p.adjust(df_xm$p,method="BH")
df_xm = as.data.frame(cbind(df_xm, p_adjust))

dataset_names <- list('metabolism_disease' = df_my, 'fungi_metabolism' = df_xm)
openxlsx::write.xlsx(dataset_names,rowNames=TRUE, file = 'D:\\真菌分析\\Plot\\data\\代谢组分析.xlsx') 

df_my<-subset(df_my, p_adjust<0.05 ) 
colnames(df_my)[3:8]=paste(colnames(df_my)[3:8],"_metabolismdisease",sep = "") # rename the columns for merge

df_xm<-subset(df_xm, p_adjust<0.05 ) 
result=merge(df_xm[,c(1,2,4:8)],df_my[,c(1,2,4:8)],by="metabolism",all = F)

metabolism<-read_dta("D:\\真菌分析\\data_pre\\meta_data.dta")
rownames(metabolism)=metabolism$SampleID
metabolism<-subset(metabolism,select=-c(SampleID))

fungi<-read_dta("D:\\真菌分析\\data_pre\\fungi_data.dta")
rownames(fungi)=fungi$SampleID
fungi<-subset(fungi,select=-c(SampleID))

diseases<-read_dta("D:\\真菌分析\\data_pre\\diseases_data.dta")
rownames(diseases)=diseases$SampleID
diseases<-subset(diseases,select=-c(SampleID,hyp))

cov<-read_dta("D:\\真菌分析\\data_pre\\cov_data.dta")
rownames(cov)=cov$SampleID
cov<-subset(cov,select=-c(SampleID))
cov$district<-as.factor(cov$district)
set.seed(123)  # 设置随机种子，确保每次运行结果一致

for(i in 1:nrow(result)){
  data=data.frame(cbind(fungi[,as.character(result$fungi[i])],diseases[,as.character(result$outcome[i])],metabolism[,as.character(result$metabolism[i])],cov))
  data=na.omit(data)
  #diet influence meta through microbiome
  colnames(data)=c("X","Y","M","age","sex","BMI","district")
  model.m=lmer(M ~ X + age + sex + BMI + (1 | district), data = data)
  model.y=glmer(Y ~ X + M + age + sex + BMI + (1 | district), data = data, family = binomial)
  summary=summary(mediate(model.m ,model.y,treat = "X", mediator = "M",boot =FALSE))
  result$Pval_mediate[i]=summary$d.avg.p
  result$coef_mediate[i]=summary$d.avg
  result$coefCI_mediate[i]=summary$d.avg.ci
  result$Pval_direct[i]=summary$z.avg.p
  result$coef_direct[i]=summary$z.avg
  result$coefCI_direct[i]=summary$z.avg.ci
  result$Pval_total[i]=summary$tau.p
  result$coef_total[i]=summary$tau.coef
  result$coefCI_total[i]=summary$tau.ci
  result$Pval_ratio[i]=summary$n.avg.p
  result$coef_ratio[i]=summary$n.avg
  result$coefCI_ratio[i]=summary$n.avg.ci
  #inverse mediate
  colnames(data)=c("X","M","Y","age","sex","BMI","district")
  model.m=glmer(M~X+age+sex+BMI+ (1 | district), data = data, family = binomial)
  model.y=lmer(Y~X+M+age+sex+BMI+ (1 | district), data = data)
  summary=summary(mediate(model.m ,model.y,treat = "X", mediator = "M",boot =FALSE ))
  result$Pval_mediate_inverse[i]=summary$d.avg.p
  result$Pval_direct_inverse[i]=summary$z.avg.p
}

result$Qval_mediate=p.adjust(result$Pval_mediate,method = "BH")
result$Qval_mediate_inverse=p.adjust(result$Pval_mediate_inverse,method = "BH")
write.table(result,file = "D:\\真菌分析\\Final results\\All data\\代谢物中介分析结果1.csv",quote = F,sep = "\t",row.names = F)
