source("https://bioconductor.org/biocLite.R")
biocLite("graph")
biocLite("Rgraphviz")

library(sna)
library(igraph)
library(spdep)
library(graph)
library(Matrix)
library(Rgraphviz)
library(network)
library(dplyr)
library(geomnet)


setwd("D:/STARNetworkModel/SNR")

nodes<-read.table("NodeXY.txt", header=T, dec=",")

load("datnetwork.rdata")
head(jasa)

load("meandata.rdata")
load("meanintense.rdata")
load("cov2means.rdata")

# extract nodewise means and compute graph related statistics
nodemeans <- unlist(net$mean.intensity)
deg <- unlist(net$n.neighbors)        # degree
g.betweenness <- igraph:::betweenness(net$graph) # betweenness
covar1.means <- as.data.frame(a$datmeans) 
covar2.means <- as.data.frame(a2$datmeans) 

snr.d <- cbind(nodemeans, deg, g.betweenness, covar1.means,covar2.means)
colnames(snr.d) <- c("nodemeans", "deg", "betweenness", "MOVI", "DIS_FAR",
                       "DIS_SALUD","DIS_EDU","DIS_PAR","ESP","EXT","VR_VIVIEND", 
                       "TE", "WIND", "E25","E75","E125","E175","E225","E275",
                       "E325","E375","E425","E475","E525","E575","E625","E675",
                       "E725","E775","E825","E875","E925","E975","E100","RAIN")
snr.dat <- cbind(snr.d, nodes)


adj<-read.table("CastellonAdj.txt", header=F)
graph <- graph_from_adjacency_matrix(as.matrix(adj), "undirected")
bad.vs<-V(graph)[igraph:::degree(graph) == 0] 
graph <-igraph:::delete.vertices(graph, bad.vs)

adj.ge1 <- as_adjacency_matrix(graph, type="both")


A <- graph.adjacency(as.matrix(adj.ge1), mode="directed",diag=FALSE)

g1  <- igraph.to.graphNEL(A) 
sM <- as(g1, "sparseMatrix") 
nb <- mat2listw(sM)$neighbours
class(nb) 
W <- nb2listw(nb, style="W", zero.policy=T) 
W2 <- nb2listw(nb, style="B", zero.policy=TRUE) 

net <- as.network(adj.ge1[20:40,20:40])
red.adj <- adj[20:45, 20:45]

g.red <- graph_from_adjacency_matrix(as.matrix(red.adj), "undirected")

nach.17 <- vector(mode="list")
for(i in 1:3)
nach.17[[i]]<-make_ego_graph(g.red, order=i, 17)[[1]]

sub.g <- fortify(g.red)
sub.g$type <- "Extract from traffic network"

ne1<-na.omit(fortify(nach.17[[1]]))
ne1$type <- "1st order neighbors"


nach2<-na.omit(fortify(nach.17[[2]]))
ne2  <- anti_join(nach2, ne1)
ne2$type <- "2nd order partial neighbors"
nach2$type <- "2nd order neighbors"

nach3<-na.omit(fortify(nach.17[[3]]))
ne3  <- anti_join(nach3, ne2)
ne3$type <- "3rd order partial neighbors"
nach3$type <- "3rd order neighbors"

# joind all neighborhoods and orgiginal graph
subs <- na.omit(rbind(sub.g, ne1, ne2, nach2))



g.data <- graph_from_data_frame(subs, directed = F)
g.data$type <- factor(g.data$type, labels = c("Extract from traffic network", "1st order neighbors",
                                              "2nd order partial neighbors", "2nd order neighbors"))

pdf("figneighbors.pdf",  paper='A4r')
ggraph(g.data, 'igraph',algorithm = 'kk')+ facet_wrap(~type) + 
  theme_bw() +geom_edge_link(edge_alpha = 0.7) +  geom_node_point(size=2, shape=1) + ggforce::theme_no_axes()
dev.off()


p<- ggraph(g.data, 'igraph',algorithm = 'kk')+ facet_wrap(~type) + 
  theme_bw() +geom_edge_link(edge_alpha = 0.7) +  geom_node_point(size=2, shape=1) + ggforce::theme_no_axes()
ggsave("figneighbors.pdf", p,  width = 11.7, height = 8.3)

# local lisa

snr.dat.ge1 <- snr.dat[-c(1095,1218),]
lM <- localmoran(snr.dat.ge1$nodemeans, W, zero.policy=TRUE, na.action = na.omit)
edit(lM)
quadrant <- vector(mode="numeric",length=nrow(lM))
cCMEDV <- snr.dat.ge1$nodemeans- mean(snr.dat.ge1$nodemeans)  
C_mI <- lM[,1] - mean(lM[,1])    
signif <- 0.5 
quadrant[cCMEDV >0 & C_mI>0] <- 4
quadrant[cCMEDV <0 & C_mI<0] <- 1      
quadrant[cCMEDV <0 & C_mI>0] <- 2
quadrant[cCMEDV >0 & C_mI<0] <- 3
quadrant[lM[,5]>signif] <- 0  

oid <- order(V(graph))

printCoefmat(data.frame(lM[oid,],
             check.names=FALSE))

snr.dat.ge1$lM <- quadrant        

# plot local morans I
pdf("mlisa.pdf")
ggplot(snr.dat.ge1, aes(cx,cy)) + geom_point(aes(colour=as.factor(lM)), shape=19, size=1.5) + 
  theme_bw() + geom_tile(aes(fill=as.factor(lM)))+ 
  scale_colour_grey(guide=F) + 
  scale_fill_manual(values=c("#ececec","#c6c6c6", "#939393", "#545454", "#000000"), 
  name="", breaks=c(0,1,2,3,4), labels=c("insignificant","low-low","low-high","high-low","high-high")) +
                scale_y_continuous(name="y-coordinate") + 
                scale_x_continuous(name="x-coordinate") 


# local net G
resG <- localG(snr.dat.ge1$nodemeans, nb2listw(include.self(nb), style = "B"))
snr.dat.ge1$lG <- unlist(as.list(round(resG,1)))

pdf("mG.pdf")
ggplot(snr.dat.ge1, aes(cx,cy)) +  geom_point(alpha=0)+ 
  theme_bw() + geom_tile()+ 
  scale_colour_grey(guide=F) + 
  scale_y_continuous(name="y-coordinate") + 
  scale_x_continuous(name="x-coordinate") + geom_text(aes(label=lG),hjust=0, 
                                                      vjust=0, size=3, check_overlap = T)
dev.off()


gggraph(g.red) + geom_edge_link() + geom_node_point(size=3) + theme_bw()


library(network)
net <- network(adj.red, directed=FALSE)
neigh<-sna:::neighborhood(g,3,neighborhood,type="total",return.all=TRUE)

create_layout(g.red)

plot(g.red, vertex.color="grey", vertex.size=15,
vertex.frame.color="black")

ggraph(g.red, 'igraph', algorithm = 'nidely') +
geom_edge_fan() +ggforce::theme_no_axes() + geom_node_point(size=2)



## global statistics SNA package
moran(snr.dat.ge1$nodemeans, W, zero.policy=TRUE)

L1 <- factor(snr.dat.ge1$nodemeans < mean(snr.dat.ge1$nodemeans), labels=c("L", "H"))
lw <- lag.listw(W, snr.dat.ge1$nodemeans) 
L2 <- factor(lw < mean(lw), labels=c("L", "H"))  
snr.dat.ge1$moran<-paste(L1, L2) 

ggplot(snr.dat.ge1, aes(cx,cy)) + geom_point(aes(colour=as.factor(moran)), shape=19, size=1.5) + 
  theme_bw() + geom_tile(aes(fill=as.factor(moran)))+ 
  scale_colour_grey(guide=F) + 
  scale_fill_manual(values=c("#ececec","#c6c6c6", "#939393", "#545454", "#000000"), 
                    name="", breaks=c(0,1,2,3,4), labels=c("insignificant","low-low","low-high","high-low","high-high")) +
  scale_y_continuous(name="y-coordinate") + 
  scale_x_continuous(name="x-coordinate") 


#### Moran I / Geary C via nacf

y <- snr.dat.ge1$nodemeans
adjm <- as(adj.ge1,"matrix")


mip <- nacf(adjm, y, "moran", partial.neighborhood = TRUE, lag.max = 40)
gcp <- nacf(adjm, y, "geary", partial.neighborhood = TRUE, lag.max = 40)

mic <- nacf(adjm, y, "moran", partial.neighborhood = FALSE, lag.max = 40)
gcc <- nacf(adjm, y, "geary", partial.neighborhood = FALSE, lag.max = 40)

par(mfrow=c(2,2))
plot(mip, type="b")
plot(gcp, type="b")
plot(mic, type="b")
plot(gcc, type="b")
par(mfrow=c(1,1))


moran.plot(y, nb2listw(nb), pch=19, ylab = "spatially lagged nodewise means intensity",
           xlab = "nodewise mean intensity", labels = FALSE)


# global tests
moran.test(y, nb2listw(nb, style="W"))
geary.test(y, nb2listw(nb, style="W"))
globalG.test(y, nb2listw(include.self(nb), style="B"))


### correlogram

lmsp <- lm(y ~ 1)
Mres <- sp.correlogram(nb, residuals(lmsp), order = 10, method = "I",
                     style = "W", zero.policy = TRUE)
plot(Mres, main="", lwd=2)

Cres <- sp.correlogram(nb, residuals(lmsp), order = 10, method = "C",
                       style = "B", zero.policy = TRUE)
plot(Cres, main="", lwd=2)

capture.output(Mres.p <- print(Mres, "bonferroni")[,5]) 
capture.output(Cres.p <- print(Cres, "bonferroni")[,5]) 

# p values of Correlograms
plot(Mres.p, type="l", xlab="lags", ylab="p-values", lwd=3, ylim=c(0,1))
abline(a=0.01, b=0, lty=2, lwd=2, col="gray")
plot(Cres.p, type="l", xlab="lags", ylab="p-values", lwd=3, ylim=c(0,1))
abline(a=0.01, b=0, lty=2, lwd=2, col="gray")

