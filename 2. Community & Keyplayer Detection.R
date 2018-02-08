###################################################################################################
# This code contains five parts:                                                                  #
# 1. (begin from Line 50 ) Pre-liminary Analysis(basic concepts of graphs)                        # 
# 2. (begin from Line 104) Community detection (10 methods & comparison)                          # 
# 3. (begin from Line 335) Analysis on each community                                             #
#                          (based on community detection results obtained by spinglass algorithm) #
# 4. (begin from Line 377) for each community, find the most effective ones to the whole graph    #                        # 
#                          (degree-based/betweenness-based/closeness-based/weight-based)          #
# 5. (begin from Line 536) for each community,find nodes most effective to their own community    #
# 6. (begin from Line 633) find top 5 keyplayers in whole graph(regradless of the communities)    #
###################################################################################################

library(igraph)
library(dplyr)
library(ggplot2)
dis_df=read.csv("C:\\Users\\Yhannah\\Desktop\\CLASS\\web analytics\\disease_sim150.csv",header=TRUE,stringsAsFactors = FALSE)

#plot
str(dis_df)
dis_df<- dis_df[,-1]
dis_df[,1] <- as.character(dis_df[,1])
dis_df[,2] <- as.character(dis_df[,2])
dis=as.matrix(dis_df)

# define function
plotG <- function(g) {
  plot(g, 
       # force-directed layout
       layout=layout.fruchterman.reingold,
       vertex.label=NA,
       vertex.label.font=0.01, vertex.size=4,
       vertex.label.color="dark red",
       vertex.color="orange",
       vertex.frame.color=FALSE, 
       vertex.label.dist=1,
       vertex.label.degree=-pi/2,
       edge.width=E(g)$weight,
       edge.color="gray",
       alpha = 0.3)
}

g=graph_from_edgelist(dis[,1:2], directed = FALSE) 
E(g)$weight=as.numeric(dis[,3])
dg <- decompose.graph(g)
g <- dg[[1]]
plotG(g)


############################
# 1. Pre-liminary Analysis #
############################
summary(g)

#vertex
vcount(g)

#transitivity
transitivity(g)

#Degree of each Node
plot(degree(g),xlab = "Nodes", ylab = "Degree")

#Degree Distribution
plot(degree_distribution(g),xlab = "Nodes", ylab = "Degree Distribution")


#Density
(2*ecount(g))/(gorder(g)*(gorder(g)-1))*100
edge_density(g)

#Diameter
diameter(g)

#Closeness
closeness(g)

#Betweenness
betweenness(g)

#Clique
clique.number(g)

#Vertex with maximum betweenness value
ver = V(g)
ver[which.max(betweenness(g))]

#vertex with maximum closeness centrality
ver[which.max(closeness(g))]


#hub score
a <- as.data.frame(hub.score(g))
name1 <- rownames(a)
hub_score <- as.data.frame(cbind(name1,hub.score(g)$vector))
hub_score <- hub_score[order(-as.numeric(as.character(hub_score$V2))),]

#authority score
b <- as.data.frame(authority.score(g))
name1 <- rownames(b)
authority_score <- as.data.frame(cbind(name1,authority.score(g)$vector))
authority_score <-authority_score[order(-as.numeric(as.character(authority_score$V2))),]

##########################
# 2. community detection #
##########################
# 1. clique-based methods
cliques(g, min=9, max=9)
plot(g)
# clique.number calculates the size of the largest clique(s).
clique.number(g)
#obeying the size limitations given in the min and max arguments. 
cliques(g, min=8)
# find the largest cliques
largest.cliques(g)

# 2. edge.betweenness.community
# Community structure detection based on edge betweenness
eb <- edge.betweenness.community(g)
eb
modularity(eb)

## 3. walktrap 
wc <- walktrap.community(g)
modularity(wc)#the biger the better

membership(wc)
plot(wc, g,layout=layout.fruchterman.reingold,
     vertex.label.font=0.005, vertex.size=3,
     vertex.label.color="dark blue",
     vertex.color="blue",
     vertex.frame.color=FALSE, 
     vertex.label=NA,
     vertex.label.dist=1,
     edge.width=E(g)$weight,
     edge.color="gray",
     alpha = 0.3)

# 4. fastgreedy.community 
#  remove multi-edges and loops
g1 <- simplify(g)
fc <- fastgreedy.community(g1)
membership(fc)
sizes(fc)
dendPlot(fc,cex=.8)
modularity(fc)

# 5.infomap Method
ic <- infomap.community(g1)
membership(ic)
communities(ic)
modularity(ic)

# 6.label.propagation.community Detection
lab <- label.propagation.community(g)
memberships <- list()
memberships$`Label propagation` <- label.propagation.community(g)
membership(lab)
sizes(lab)
plot(lab,g)
modularity(lab)
# get the Number of members in the top Clusters of graph
topClusters_lab <- table(lab$membership) %>% sort(decreasing=TRUE) %>% head(10)
topClusters_lab[1:10]
plot(topClusters_lab, main="Cluster size", ylab="Number of members", type="b", lwd=2)

# 7. spinglass.community
set.seed(2324)
sc <- spinglass.community(g, spins=9)
memberships$`Spinglass` <- sc$membership
modularity(sc)
membership(sc)
sizes(sc)

# 8.K-Core Decomposition
kc <- coreness(g, mode="all")
plot(g, vertex.size=kc, vertex.label=kc)
modularity(g,kc)


# 9. leading.eigenvector.community
lec <- leading.eigenvector.community(g,options=list(maxiter=1000000))
memberships$`Leading eigenvector` <- lec$membership
modularity(lec)
membership(lec)
sizes(lec)


# 10. Grivan-Newman algorithm
# -- 1st we calculate the edge betweenness, merges, etc...
ebc <- edge.betweenness.community(g, directed=F)
ebc
# Now we have the merges/splits and we need to calculate the modularity
# for each merge for this we'll use a function that for each edge
# removed will create a second graph, check for its membership and use
# that membership to calculate the modularity
mods <- sapply(0:ecount(g), function(i){
  g2 <- delete.edges(g, ebc$removed.edges[seq(length=i)])
  cl <- clusters(g2)$membership
  modularity(g,cl)
})
# -- we can now plot all modularities
par(mfrow=c(1,1))
plot(mods, pch=20)

# -- color the nodes according to their membership
g2<-delete.edges(g, ebc$removed.edges[seq(length=which.max(mods)-1)])
V(g)$color=clusters(g2)$membership

# Let's choose a layout for the graph
g$layout <- layout.fruchterman.reingold

# -- plot it
plot(g, vertex.label=NA, vertex.size=3,vertex.label=NA,
     alpha = 0.3)
modularity(g,clusters(g2)$membership)

mods

###########################
# compare all the methods #
###########################
par(mfrow=c(3,3))
plot(wc, g,layout=layout.fruchterman.reingold,
     vertex.label.font=0.005, vertex.size=6,
     vertex.label.color="dark blue",
     vertex.color="blue",
     vertex.frame.color=FALSE, 
     vertex.label=NA,
     vertex.label.dist=1,
     edge.width=E(g)$weight,
     edge.color="gray",
     alpha = 0.1,rotate = FALSE,
     main = "Walktrap",sub = modularity(wc))
plot(fc, g,layout=layout.fruchterman.reingold,
     vertex.label.font=0.005, vertex.size=6,
     vertex.label.color="dark blue",
     vertex.color="blue",
     vertex.frame.color=FALSE, 
     vertex.label=NA,
     vertex.label.dist=1,
     edge.width=E(g)$weight,
     edge.color="gray",
     alpha = 0.1,rotate = FALSE,
     main = "Fastgreedy",sub = modularity(fc))
plot(ic, g,layout=layout.fruchterman.reingold,
     vertex.label.font=0.005, vertex.size=6,
     vertex.label.color="dark blue",
     vertex.color="blue",
     vertex.frame.color=FALSE, 
     vertex.label=NA,
     vertex.label.dist=1,
     edge.width=E(g)$weight,
     edge.color="gray",
     alpha = 0.1,rotate = FALSE,
     main = "Infomap",sub = modularity(ic))
plot(eb, g,layout=layout.fruchterman.reingold,
     vertex.label.font=0.005, vertex.size=6,
     vertex.label.color="dark blue",
     vertex.color="blue",
     vertex.frame.color=FALSE, 
     vertex.label=NA,
     vertex.label.dist=1,
     edge.width=E(g)$weight,
     edge.color="gray",
     alpha = 0.1,rotate = FALSE,
     main = "Edge Betweenness",sub = modularity(eb))
plot(lab, g,layout=layout.fruchterman.reingold,
     vertex.label.font=0.005, vertex.size=6,
     vertex.label.color="dark blue",
     vertex.color="blue",
     vertex.frame.color=FALSE, 
     vertex.label=NA,
     vertex.label.dist=1,
     edge.width=E(g)$weight,
     edge.color="gray",
     alpha = 0.1,rotate = FALSE,
     main = "Label Propagation",sub = modularity(lab))
plot(sc, g,layout=layout.fruchterman.reingold,
     vertex.label.font=0.005, vertex.size=6,
     vertex.label.color="dark blue",
     vertex.color="blue",
     vertex.frame.color=FALSE, 
     vertex.label=NA,
     vertex.label.dist=1,
     edge.width=E(g)$weight,
     edge.color="gray",
     alpha = 0.1,rotate = FALSE,
     main = "spinglass",sub = modularity(sc))
plot(lec, g,layout=layout.fruchterman.reingold,
     vertex.label.font=0.005, vertex.size=6,
     vertex.label.color="dark blue",
     vertex.color="blue",
     vertex.frame.color=FALSE, 
     vertex.label=NA,
     vertex.label.dist=1,
     edge.width=E(g)$weight,
     edge.color="gray",
     alpha = 0.1,rotate = FALSE,
     main = "Leading_eigenvector",sub = modularity(lec))
plot(g, vertex.size=kc,main = "K-Core",
     vertex.label=NA,
     sub = modularity(g,kc),alpha = 0.1)
plot(g, vertex.size=clusters(g2)$membership,main = "Newman",
     vertex.label=NA,
     sub = modularity(g,clusters(g2)$membership),alpha = 0.1)
dev.off()

# compare all the methods
par(mfrow=c(1,1))
method <- c("Walktrap","Fastgreedy", "Infomap","Edge Betweenness", "Label Propagation", "K-Core","spinglass","Leading_eigenvector","Newman")
modularity <- c(modularity(wc),modularity(fc),modularity(ic),modularity(eb),modularity(lab),modularity(g,kc),modularity(sc),modularity(lec),modularity(g,clusters(g2)$membership))
group_size <- c(count(as.data.frame(table(wc$membership)))[[1]],
               count(as.data.frame(table(fc$membership)))[[1]],
               count(as.data.frame(table(ic$membership)))[[1]],
               count(as.data.frame(table(eb$membership)))[[1]],
               count(as.data.frame(table(lab$membership)))[[1]],
               count(as.data.frame(table(kc)))[[1]],
               count(as.data.frame(table(sc$membership)))[[1]],
               count(as.data.frame(table(lec$membership)))[[1]],
               count(as.data.frame(table(clusters(g2)$membership)))[[1]])
mm <- cbind(as.data.frame(method),as.data.frame(modularity),as.data.frame(group_size))
mm$method <- as.character(mm$method)
mm$modularity <- as.numeric(as.character(mm$modularity))


ggplot(data = mm,aes(x=mm$group_size ,y=mm$modularity,fill=factor(mm$method))) +
  geom_point(size=5,shape=24,alpha=0.5,show.legend = TRUE)+
  geom_text(aes(label = mm$method), vjust = 1.4,  hjust ="inward",colour = 'dark red', size = 4)+
  labs(x="group_size",y="Modularity",
       title = "Group size & modularities for each method")



#################################
# 3. Analysis on each community #
#################################
# based on spinglass
sc$membership
sizes(sc)
modularity(sc)
# return graphs of each commuity
ind_com <- function(i){
    c1<- as.data.frame(sc[[i]])
    c2<- as.data.frame(sc[[i]])
    colnames(c1) <- "d1"
    cb <- inner_join(dis_df,c1)
    colnames(c2) <- "d2"
    cb2 <- inner_join(cb,c2)
    cb2[,1] <- as.character(cb2[,1])
    cb2[,2] <- as.character(cb2[,2])
    cb_f=as.matrix(cb2)
    g=graph_from_edgelist(cb_f[,1:2], directed = FALSE) 
    E(g)$weight=as.numeric(cb_f[,3])
    return(g)
}

ind_com(1)
# plot 9 communities
par(mfrow=c(3,3))
for(i in 1:9){
  plot(ind_com(i),vertex.size = degree(ind_com(i)),
       vertex.label.color ="dark red",
       vertex.frame.color=FALSE,vertex.label=NA,
      main=i)
}

modularity(sc)
# plot certain community
par(mfrow=c(1,1))
plot(ind_com(9),vertex.size = degree(ind_com(9)),
     vertex.label.font=1,edge.width=E(g)$weight,
     vertex.label.color="dark red",vertex.label.dist=2)
sc
sc[[4]]
betweenness(g)
########################
# 4. detect keyplayers #
########################
#find keyplayers for each community newman fast algorithm
#cm <- as.data.frame(clusters(g2)$membership)
#name <- rownames(cm)
#cmf <- cbind(name,cm)
#colnames(cmf)[2] <- "community"
#c <- cmf[order(cmf$community ),]

## caculate degrees of every disease
#degree_nf <- as.data.frame(degree(g))
#name <- rownames(degree_nf)
#degree_nf<- cbind(name,degree_nf)
#colnames(degree_nf)[2] <- "degree"
#degree_nf <- degree_nf[order(-degree_nf$degree ),]

# join two tables and to find the keyplayers who have biggest degree in each community
#library(plyr)
#dc <- join(c,degree_nf)
#dclist <- dc %>% group_by(community) %>% ungroup()
#keyplayer <- dclist %>% group_by(community) %>% 
#             top_n(3, degree) %>% ungroup() %>%
#             arrange(community, -degree)
#keyplayer_top1 <- dclist %>% group_by(community) %>% 
#                  top_n(1, degree) %>% ungroup() %>% 
#                  arrange(community, -degree)

###############################################################################
# detect keyplayer based on community results obtained by spinglass algorithm #
###############################################################################

## 1. four measuring standards for detect keyplayer: degree/Betweenness/Closeness/weight

####### - degree_based - #######
# assgin diseases into each community
d_dis <- as.data.frame(cbind(sc$names,sc$membership))
colnames(d_dis) <- c("name","community")
# find number of degrees for every disease 
degree_nf <- as.data.frame(degree(g))
name <- rownames(degree_nf)
degree_nf<- cbind(name,degree_nf)
colnames(degree_nf)[2] <- "degree"
# join two tables
dc <- inner_join(d_dis,degree_nf)
dclist<- dc %>% group_by(community) %>% ungroup()
# get the keyplayers for each community
# (who have most number of degrees in each community)
keyplayer <- dclist %>% group_by(community) %>% 
             top_n(3, degree) %>% ungroup() %>% 
             arrange(community, -degree)

# get the top1 keyplayers for each community
keyplayer_top1 <- dclist %>% group_by(community) %>% 
                  top_n(1, degree) %>% ungroup() %>%
                  arrange(community, -degree)

####### - Betweenness_based - ####### 
betweenness_nf <- as.data.frame(betweenness(g))
name <- rownames(betweenness_nf)
betweenness_nf <- cbind(name,betweenness_nf)
colnames(betweenness_nf)[2] <- "betweenness"
bc <- inner_join(d_dis,betweenness_nf)
bclist <- bc %>% group_by(community) %>% ungroup()
# get the keyplayers for each community(who have largest betweenness in each community)
keyplayer_be <- bclist %>% group_by(community) %>% 
                top_n(3, betweenness) %>% ungroup() %>%
                arrange(community, -betweenness)
# get the top1 keyplayers for each community
keyplayer_top1_be <- bclist %>% group_by(community) %>% 
                     top_n(1, betweenness) %>% ungroup() %>%
                     arrange(community, -betweenness)

####### - Closeness_based - #######
closeness_nf <- as.data.frame(closeness(g))
name <- rownames(closeness_nf)
closeness_nf <- cbind(name,closeness_nf)
colnames(closeness_nf)[2] <- "closeness"
cc <- inner_join(d_dis,closeness_nf)
cclist <- cc%>%group_by(community)%>%ungroup()
# get the keyplayers for each community(who have largest closeness in each community)
keyplayer_cc <- cclist %>% group_by(community) %>% 
                top_n(3, closeness) %>% ungroup() %>%
                arrange(community, -closeness)
# get the top1 keyplayers for each community
keyplayer_top1_cc <- cclist %>% group_by(community) %>% 
                     top_n(1, closeness) %>% ungroup() %>%
                     arrange(community, -closeness)

####### - weight_based - ########
# caculate the weight for each distinct disease in column 1 and column 2 
d1s <- dis_df[,c(1,3)]
d2s <- dis_df[,c(2,3)]
ds_p1 <- d1s %>% group_by(d1) %>% summarise_all(funs(sum)) %>% ungroup()
ds_p2 <- d2s %>% group_by(d2) %>% summarise_all(funs(sum)) %>% ungroup()
colnames(ds_p2)[1] <- "d1"
ds <- rbind(ds_p1,ds_p2)
# sum the weight for each distinct disease 
weight <- ds %>% group_by(d1) %>% summarise_all(funs(sum)) %>% ungroup()
colnames(weight)[1] <- "name"
colnames(weight)[2] <- "w"
cw <- inner_join(d_dis,weight)
cwlist <- cw %>% group_by(community) %>% ungroup()
# get the keyplayers for each community(who have largest weight in each community)
keyplayer_cw <- cwlist %>% group_by(community) %>% 
                top_n(3, w) %>% ungroup() %>%
                arrange(community, -w)
# get the top1 keyplayers for each community
keyplayer_top1_cw <- cwlist %>% group_by(community) %>% 
                     top_n(1, w) %>% ungroup() %>%
                     arrange(community, -w)


## 2 . plot top3 keyplayers for 9 communities respectively
keyplayer %>%
  mutate(name = reorder(name, degree)) %>%
  ggplot(aes(name,degree,fill=factor(community))) +
  geom_col(show.legend = TRUE,alpha = 0.6) +
  labs(title = "Top3 Keyplayers for 9 communities(degree)")+
  #facet_wrap(~ community, scales = "free") +
  coord_flip(expand = TRUE)
keyplayer_be %>%
  mutate(name = reorder(name, betweenness)) %>%
  ggplot(aes(name,betweenness,fill=factor(community))) +
  geom_col(show.legend = TRUE,alpha = 0.6) +
  labs(title = "Top3 Keyplayers for 9 communities(betweenness)")+
  #facet_wrap(~ community, scales = "free") +
  coord_flip(expand = TRUE)
keyplayer_cc %>%
  mutate(name = reorder(name, closeness)) %>%
  ggplot(aes(name,closeness,fill=factor(community))) +
  geom_col(show.legend = TRUE,alpha = 0.6) +
  labs(title = "Top3 Keyplayers for 9 communities(closeness)")+
  #facet_wrap(~ community, scales = "free") +
  coord_flip(expand = TRUE)
keyplayer_cw %>%
  mutate(name = reorder(name, w)) %>%
  ggplot(aes(name,w,fill=factor(community))) +
  geom_col(show.legend = TRUE,alpha = 0.6) +
  labs(title = "Top3 Keyplayers for 9 communities(weight)")+
  #facet_wrap(~ community, scales = "free") +
  coord_flip(expand = TRUE)


## 3. compare the keyplayers based on different measures
adc<- keyplayer_top1[,c(1,2)]
bec <- keyplayer_top1_be[,c(1,2)]
cca <- keyplayer_top1_cc[,c(1,2)]
cww <- keyplayer_top1_cw[,c(1,2)]
df <- inner_join(adc,bec,by="community")
dff <- inner_join(df,cca,by="community")
dff <- inner_join(dff,cww,by="community")
colnames(dff)[1] <- "degree_based"
colnames(dff)[3] <- "betweenness_based"
colnames(dff)[4] <- "closeness_based"
colnames(dff)[5] <- "weight_based"
dff <- dff[order(as.integer(as.character(dff$community))),]
dff <- dff[,c(2,1,3:5)]

#############################################################################
# 5. find keyplayers which are most effective in all 9 groups respectively  #
#############################################################################
#1. degree
keyplayer_in <- data.frame()
for(i in 1:9){
  a <- as.data.frame(degree(ind_com(i)))
  colnames(a)[1] <- "degree"
  name_in <- rownames(a)
  b <- cbind(name_in,a)
  b <- b[order(-b$degree),]
  ba <- cbind(b[b$degree==max(b$degree),],i)
  keyplayer_in <- rbind(keyplayer_in,ba)
}
rownames(keyplayer_in) <- c(1:nrow(keyplayer_in))
colnames(keyplayer_in)[3] <- "community"
keyplayer_in

#2. Closeness
keyplayer_inc <- data.frame()
for(i in 1:9){
  a <- as.data.frame(closeness(ind_com(i)))
  colnames(a)[1] <- "closeness"
  name_in <- rownames(a)
  b <- cbind(name_in,a)
  b <- b[order(-b$closeness),]
  ba <- cbind(b[b$closeness==max(b$closeness),],i)
  keyplayer_inc <- rbind(keyplayer_inc,ba)
}
rownames(keyplayer_inc) <- c(1:nrow(keyplayer_inc))
colnames(keyplayer_inc)[3] <- "community"
keyplayer_inc


#3. Betweenness
keyplayer_inb <- data.frame()
for(i in 1:9){
  a <- as.data.frame(betweenness(ind_com(i)))
  colnames(a)[1] <- "betweenness"
  name_in <- rownames(a)
  b <- cbind(name_in,a)
  b <- b[order(b$betweenness),]
  ba <- cbind(b[b$betweenness==max(b$betweenness),],i)
  keyplayer_inb <- rbind(keyplayer_inb,ba)
}
rownames(keyplayer_inb) <- c(1:nrow(keyplayer_inb))
colnames(keyplayer_inb)[3] <- "community"
keyplayer_inb

#4. weight
keyplayer_inw <- data.frame()
# return dataframe for each commuity
for(i in 1:9){
  c1 <- data.frame()
  c2 <- data.frame()
  c1<- as.data.frame(sc[[i]])
  c2<- as.data.frame(sc[[i]])
  colnames(c1) <- "d1"
  cb <- inner_join(dis_df,c1)
  colnames(c2) <- "d2"
  cb2 <- inner_join(cb,c2)
  cb2[,1] <- as.character(cb2[,1])
  cb2[,2] <- as.character(cb2[,2])
  d1s <- cb2[,c(1,3)]
  d2s <- cb2[,c(2,3)]
  ds_p1 <- d1s %>% group_by(d1) %>% summarise_all(funs(sum)) %>% ungroup()
  ds_p2 <- d2s %>% group_by(d2) %>% summarise_all(funs(sum)) %>% ungroup()
  colnames(ds_p2)[1] <- "d1"
  ds <- rbind(ds_p1,ds_p2)
  # sum the weight for each distinct disease 
  weight <- ds %>% group_by(d1) %>% summarise_all(funs(sum)) %>% ungroup()
  sd <- cbind(weight,i)
  keyplayer_inw <- rbind(keyplayer_inw,sd)
}
colnames(keyplayer_inw)<- c("name","weight","community")
keyplayer_inw_top <- keyplayer_inw %>% group_by(community) %>% 
  top_n(1,weight) %>% ungroup() %>%
  arrange(community, -weight)

###########
# compare #
###########
kin<- keyplayer_in[,c(1,3)]
kib <- keyplayer_inb[,c(1,3)]
kic <- keyplayer_inc[,c(1,3)]
kiw <- keyplayer_inw_top[,c(1,3)]
dfk <- inner_join(kin,kib,by="community")
dfk <- inner_join(dfk,kic,by="community")
dfk <- inner_join(dfk,kiw,by="community")
colnames(dfk)[1] <- "degree_based"
colnames(dfk)[3] <- "betweenness_based"
colnames(dfk)[4] <- "closeness_based"
colnames(dfk)[5] <- "weight_based"
dfk <- dfk[order(as.integer(as.character(dfk$community))),]
dfk <- dfk[,c(2,1,3:5)]


######################################################################
# 6. find top 5 keyplayers in whole graph(ignore the communities)    #
######################################################################
library(influenceR)
keyplayer(g,k=5)
