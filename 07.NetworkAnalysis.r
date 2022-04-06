#  07. Network Analysis
## Figure 5, Figure S8, Figure S3, Table S6, Table S7

## Links to references
### https://github.com/zdk123/SpiecEasi
### https://loimai.github.io/BBobs/Lemonnier_Ushant_front_16S.html

## Package loading
source("https://raw.githubusercontent.com/genomewalker/osd2014_analysis/master/osd2014_16S_asv/lib/graph_lib.R")

## Set Path
spiec <- file.path(paste(path_phy, "/SPIEC-EASI", sep=""))
dir.create(spiec, showWarnings=FALSE)
setwd(spiec)

## remove OV1 and OV2
obj_ps_noOV <- subset_samples(obj_ps, Loc_sec!="Oasis_Valley")

## remove OTUs found majority in Oasis Valley (calculated by: max count in OV / max count in rest of sites > 0.7)
badTaxa <- c('OTU_21', 'OTU_28', 'OTU_33', 'OTU_52', 'OTU_339', 'OTU_97', 'OTU_207', 'OTU_1766', 'OTU_4946', 'OTU_5369', 'OTU_645', 'OTU_5863', 'OTU_194', 'OTU_1164', 'OTU_5465', 'OTU_5306', 'OTU_251', 'OTU_5961')
goodTaxa <- setdiff(taxa_names(obj_ps_noOV), badTaxa)
obj_ps_noOV <- prune_taxa(goodTaxa, obj_ps_noOV)

## Retrieve the 'Common OTUs'
obj_ps_filt3cv = filter_taxa(obj_ps, function(x) sd(x)/mean(x) > 3.0, TRUE) #Filter the taxa using a cutoff of 3.0 for the Coefficient of Variation
obj_ps_filt3cv = filter_taxa(obj_ps_filt3cv, function(x) sum(x > 1) > (0.05*length(x)), TRUE) #Remove taxa not seen more than 1 times in at least 5% of the samples (about n=2).
obj_ps_filt3cv
input <- as.matrix(t(otu_table(obj_ps_filt3cv)))

## Run SPIEC-EASI
print("===== SPIEC-EASI: METHOD MB STARS =====")
Sys.time()
se.mb.stars <- spiec.easi(input, method='mb', sel.criterion = "stars", verbose = TRUE, pulsar.select=TRUE, lambda.min.ratio=0.01, nlambda=20, pulsar.params=list(rep.num=50, ncores=1))

print("===== SPIEC-EASI: METHOD GLASSO STARS =====")
Sys.time()
se.gl.stars <- spiec.easi(input, method='glasso', sel.criterion = "stars", verbose = TRUE, lambda.min.ratio=0.01, nlambda=20, pulsar.params=list(rep.num=50, ncores=1))

print("===== SPIEC-EASI: METHOD MB BSTARS =====")
Sys.time()
se.mb.bstars <- spiec.easi(input, method='mb', sel.criterion = 'bstars', verbose = TRUE, pulsar.select=TRUE, lambda.min.ratio=0.01, nlambda=20, pulsar.params=list(rep.num=50, ncores=1))

print("===== SPIEC-EASI: METHOD GLASSO BSTARS =====")
Sys.time()
se.gl.bstars <- spiec.easi(input, method='glasso', sel.criterion = 'bstars', verbose = TRUE, lambda.min.ratio=0.01, nlambda=20, pulsar.params=list(rep.num=50, ncores=1))

## Evaluate the weights on edges in the network
secor  <- cov2cor(getOptCov(se.gl.stars))
sebeta <- symBeta(getOptBeta(se.mb.stars), mode='maxabs')
elist.gl <- summary(triu(secor*getRefit(se.gl.stars), k=1))
elist.mb <- summary(sebeta)

secor.b  <- cov2cor(getOptCov(se.gl.bstars))
sebeta.b <- symBeta(getOptBeta(se.mb.bstars), mode='maxabs')
elist.gl.b <- summary(triu(secor*getRefit(se.gl.bstars), k=1))
elist.mb.b <- summary(sebeta)

hist(elist.mb[,3], main='', xlab='edge weights', col='forestgreen')
hist(elist.mb.b[,3], main='', xlab='edge weights', col='grey')

hist(elist.gl[,3], main='', xlab='edge weights', col='forestgreen')
hist(elist.gl.b[,3], main='', xlab='edge weights', col='grey')

## Create igraph objects
ig.mb <- adj2igraph(getRefit(se.mb.stars))
ig.mb.b <- adj2igraph(getRefit(se.mb.bstars))
ig.gl <- adj2igraph(forceSymmetric(getRefit(se.gl.stars)))
ig.gl.b <- adj2igraph(forceSymmetric(getRefit(se.gl.bstars)))

## Visualize using igraph plotting
### set size of vertex proportional to clr-mean
vsize    <- rowMeans(clr(input, 1))+6
am.coord <- layout.fruchterman.reingold(ig.mb)

par(mfrow=c(2,2))
plot(ig.mb, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="MB-STARS")
plot(ig.mb.b, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="MB-BSTARS")
plot(ig.gl, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="Glasso-STARS")
plot(ig.gl.b, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="Glasso-BSTARS")

## Look at the degree statistics from the networks inferred by each method (Figure S3)
dd.gl <- degree.distribution(ig.gl)
dd.mb <- degree.distribution(ig.mb)
dd.gl.b <- degree.distribution(ig.gl.b)
dd.mb.b <- degree.distribution(ig.mb.b)

plot(0:(length(dd.gl)-1), dd.gl, ylim=c(0,.1), type='b', ylab="Frequency", xlab="Degree", main="Degree Distributions")
points(0:(length(dd.mb)-1), dd.mb, col="red" , type='b')
legend("topright", c("Glasso-STARS", "MB-STARS"), col=c("black", "red"), pch=1, lty=1)

## Get the weights
se.cor  <- cov2cor(as.matrix(getOptCov(se.gl.stars)))
weighted.adj.mat <- se.cor*getRefit(se.gl.stars)
heatmap(as.matrix(weighted.adj.mat))

## Plot the network
grph.unweighted <- simplify(adj2igraph(forceSymmetric(getRefit(se.gl.bstars)), rmEmptyNodes = TRUE, vertex.attr=list(name=taxa_names(obj_ps_filt3cv))))

se.cor <- cov2cor(as.matrix(getOptCov(se.gl.stars)))
grph.weighted <- simplify(adj2igraph(forceSymmetric(se.cor*getRefit(se.gl.bstars)), rmEmptyNodes = TRUE, vertex.attr=list(name=taxa_names(obj_ps_filt3cv))))

plot(grph.weighted, vertex.size=1, vertex.label=NA)
grph.plot <- grph.weighted
V(grph.plot)
E(grph.plot)

### Custom colors
customvermillion<-rgb(213/255,94/255,0/255)
custombluegreen<-rgb(0/255,158/255,115/255)
customblue<-rgb(0/255,114/255,178/255)
customskyblue<-rgb(86/255,180/255,233/255)
customreddishpurple<-rgb(204/255,121/255,167/255)

### Draft plots
V(grph.plot)$size <- (degree(grph.plot) + 1) # the +1 is to avoid size zero vertices
V(grph.plot)$color <- "black"
plot(grph.plot, vertex.label=NA)

E(grph.plot)$color <- custombluegreen
E(grph.plot)$color[E(grph.plot)$weight<0] <- customreddishpurple
E(grph.plot)$width <- abs(E(grph.plot)$weight)*50
plot(grph.plot, vertex.label=NA)

### Remove edges with very low weight 
weight_threshold <- 0.02
grph.plot <- igraph::delete.edges(grph.plot,which(abs(E(grph.plot)$weight)<weight_threshold)) 
grph.pos <- igraph::delete.edges(grph.plot,which(E(grph.plot)$weight<0)) #Remove negative edges 
plot(grph.pos, vertex.label=NA)

### Clean up the plot size
V(grph.pos)$size <- V(grph.pos)$size/3
plot(grph.pos, vertex.label=NA, edge.curved=0.5)

### Remove unconnected vertices
grph.pos <- igraph::delete.vertices(grph.pos,which(degree(grph.pos)<1)) 
plot(grph.pos, vertex.label=NA)

### Plot with colors for points
grph.pos_deg <- degree(grph.pos, v=V(grph.pos), mode="all")
fine = 500 # this will adjust the resolving power.
graphCol = viridis(fine)[as.numeric(cut(grph.pos_deg, breaks = fine))] #this gives you the colors you want for every point
plot(grph.pos, vertex.color=graphCol, edge.color="black", vertex.label=NA, layout=layout_with_fr(grph.pos))

### Plot betweenness
grph.pos_bw<-betweenness(grph.pos, directed=F)
graphCol = viridis(fine)[as.numeric(cut(grph.pos_bw,breaks = fine))]
plot(grph.pos, vertex.color=graphCol, vertex.label=NA, edge.color="black", layout=layout_with_fr(grph.pos))

### Plot transistivity
grph.pos_tran<-transitivity(grph.pos, type="local")
graphCol = viridis(fine)[as.numeric(cut(grph.pos_tran,breaks = fine))]
plot(grph.pos, vertex.color=graphCol, vertex.label=NA, edge.color="black", layout=layout_with_fr(grph.pos))
grph.pos_tran_gl<-transitivity(grph.pos, type="global")
grph.pos_tran_gl

### Layout with fruchterman-reingold algorithm
plot(grph.pos, vertex.label=NA, layout=layout_with_fr(grph.pos))

## Check degree distributions
dd.grph.pos <- degree.distribution(grph.pos)
plot(0:(length(dd.grph.pos)-1), dd.grph.pos, type='b',
      ylab="Frequency", xlab="Degree", main="Degree Distributions")

## Clustering algorithms
### Greedy algorithm
grph.pos.greedy <- cluster_fast_greedy(grph.pos, weights=E(grph.pos)$weight)
modularity(grph.pos.greedy)
sizes(grph.pos.greedy)
colourCount = length(unique(grph.pos.greedy$membership)) # this will adjust the resolving power.
cluster_col = rainbow(colourCount)[as.numeric(cut(grph.pos.greedy$membership,breaks = colourCount))]
plot(grph.pos, vertex.color=cluster_col, vertex.label=NA, edge.color="black", layout=layout_with_fr(grph.pos))

### Walktrap clustering
grph.pos.walktrap <- cluster_walktrap(grph.pos, weights = E(grph.pos)$weight)
modularity(grph.pos.walktrap)
sizes(grph.pos.walktrap)

### Community structure detecting based on the leading eigenvector of the community matrix
grph.pos.eigen <- cluster_leading_eigen(grph.pos, weights = E(grph.pos)$weight)
modularity(grph.pos.eigen)
sizes(grph.pos.eigen)

### Community structure detection based on edge betweenness
grph.pos.edge <- cluster_edge_betweenness(grph.pos, weights = E(grph.pos)$weight)
modularity(grph.pos.edge)
sizes(grph.pos.edge)

### Communities based on propagating labels
grph.pos.prop <- cluster_label_prop(grph.pos, weights = E(grph.pos)$weight)
modularity(grph.pos.prop)
sizes(grph.pos.prop)

### Louvain clustering method: https://towardsdatascience.com/louvain-algorithm-93fde589f58c
grph.pos.louvain <- cluster_louvain(grph.pos, weights=E(grph.pos)$weight)
modularity(grph.pos.louvain)
sizes(grph.pos.louvain)
colourCount = length(unique(grph.pos.louvain$membership)) # this will adjust the resolving power.
cluster_col = rainbow(colourCount)[as.numeric(cut(grph.pos.louvain$membership,breaks = colourCount))]
plot(grph.pos, vertex.color=cluster_col, vertex.label=NA, edge.color="black", layout=layout_with_fr(grph.pos))

## Go to Cytoscape (Figure 5, Figure S8)
library(RCy3)
cytoscapePing()
createNetworkFromIgraph(grph.pos,"myIgraph")

## Plot distribution of OTUs per cluster (Figure 5, Figure S8)
otu_loc_clust_df <- read.delim("Data_for_donuts.txt")
rownames(otu_loc_clust_df) <- otu_loc_clust_df$OTU_ID
otu_loc_clust_df <- otu_loc_clust_df[,-1]

cv <- function(x) 100*(sd(x)/mean(x))
calc_all <- otu_loc_clust_df %>% group_by(Cluster_ID) %>% summarise_each(funs(mean, sd, cv))
calc_mean <- otu_loc_clust_df %>% group_by(Cluster_ID) %>% summarise_each(funs(mean))

### collapse columns to rows
calc_melt <- reshape2::melt(calc_mean, id.vars=c("Cluster_ID"))

### fix labels
calc_melt$variable <- gsub('FF...Yucca.Flat','Frenchman & Yucca Flat', calc_melt$variable)
calc_melt$variable <- gsub('\\.',' ', calc_melt$variable)
calc_melt$variable <- gsub('Frenchman   Yucca Flat','Frenchman & Yucca Flat', calc_melt$variable)

### fix cluster column
calc_melt$Cluster_ID <- as.factor(calc_melt$Cluster_ID)
unique(calc_melt$variable)
calc_melt$variable <- factor(calc_melt$variable, ordered = TRUE, levels = c("Amargosa Valley", "Ash Meadows", "Death Valley", "Frenchman & Yucca Flat", "Pahute Mesa", "Rainier Mesa", "Spring Mountains", "Oasis Valley"))
head(calc_melt)

### Set up theme
plot_theme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 1, linetype = "solid"),
    panel.border = element_rect(colour="black", size=1, fill=NA),
    panel.spacing.x=unit(2, "lines"),
    panel.grid.major = element_line(size = 0),
    panel.grid.minor = element_line(size = 0),
    axis.text = element_text(size=18, colour="black"),
    axis.title = element_text(face="bold", size=20),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
    legend.position="right",
    legend.key = element_rect(fill = "white"),
    legend.title = element_text(face="bold", size=20),
    legend.text = element_text(size=20))
plot_guides <- guides(colour=FALSE, fill = guide_legend(ncol=1))
plot_nomargins_y <- scale_y_continuous(expand = expansion(mult = c(0, 0)), labels = function(x) paste0(x*100, "%"))
plot_nomargins_x <- scale_x_discrete(expand = expansion(mult = c(0, 0)))

### Plot
ggplot(data=calc_melt, aes(x=2, y=value, fill=variable)) +
    scale_fill_lancet(name = "Location") + 
    geom_bar(aes(), stat="identity", position="fill", color="black", width=0.8) +
    xlim(0.5,2.5) +
    coord_polar("y") + 
    theme_void() +
    facet_wrap(~ Cluster_ID)

save_file_plot <- paste("avg.abund.otus.module.loc.plot.donut.pdf", sep="") #change the file name if need to
ggsave(save_file_plot, scale = 1, width = 13, height = 10, units = c("in"), dpi = 300)
