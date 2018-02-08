library(dplyr)
library(splitstackshape)
library(dummies)
library(igraph)
original=read.csv("~/Desktop/课程/3/web/report/ncomms5212s6.csv",header=TRUE,stringsAsFactors = FALSE)
head(original)


disease <- original[,c(1,3,4)]
head(disease)
dis<-cSplit(disease,"Disease.Terms",";",direction = "long")
head(dis)
dis<-cSplit(dis,"Symptom.terms",";",direction = "long")

# one-hot encoding
dis_feature <- dummy.data.frame(dis, names=c("Symptom.terms"), sep="_")
dis_feature <- dis_feature[,-1]

dis_fea <- dis_feature %>%
  group_by(Disease.Terms) %>%
  summarise_all(funs(sum))
str(dis_fea)

# caculate similarity
dis_fea_n <-dis_fea[,-1]
dis_fea_m <- t(dis_fea_n)
dis_sim_per <- cor(dis_fea_m,use="pairwise.complete.obs",method="pearson")

name <- as.character(dis_fea$Disease.Terms)
colnames(dis_sim_per) <- name
rownames(dis_sim_per) <- name


nrow(dis_sim_per)
ncol(dis_sim_per)


dis_fr <- dis_sim_per[1:150,1:150]
name[1]
class(name)

# convert into data frame

mydata1<-data.frame()
dis_sim<-data.frame()
for (i in 1:149){
  for(j in (i+1):150){
    if(dis_fr[i,j]>0 ){
      d1<-name[i]
      d2<-name[j]
      sim<- dis_fr[i,j]
      mydata1<-cbind(d1,d2,sim)
      dis_sim <- rbind(dis_sim,mydata1)}
  }
}

View(dis_sim)

write.csv(dis_sim,"~/Desktop/课程/3/web/report/disease_sim150.csv")






