#  08. Ecological Process
## Figure 6, Figure S2, Table S8

## Links to references
### https://github.com/danczakre/ShaleViralEcology
### https://rdrr.io/cran/NST/man/NST-package.html
### https://github.com/stegen/Stegen_etal_ISME_2013
### https://github.com/FranzKrah/raup_crick
### http://kembellab.ca/r-workshop/biodivR/SK_Biodiversity_R.html
### https://github.com/danczakre/AquiferEcology/blob/master/EPA_BNTI_multithreaded.R
### https://doi.org/10.1093/femsec/fiz160
### https://doi.org/10.1073/pnas.1414261112

## Set path
determ <- file.path(paste(path_phy, "/determ_v_stoch", sep=""))
dir.create(determ, showWarnings=FALSE)
setwd(determ)

## Rarefy
obj_ps_rare <- rarefy_even_depth(obj_ps, rngseed=711, sample.size=min(sample_sums(obj_ps)), replace=F)
obj_ps_rare

## Filter by DVRFS groundwater basin
physeq_PMOV <- subset_samples(obj_ps_rare, Loc_sec=="Pahute_Mesa" | Loc_sec=="Oasis_Valley")
physeq_PMOV <- phyloseq_filter_sample_wise_abund_trim(physeq_PMOV, minabund = 1, rm_zero_OTUs = TRUE)
physeq_PMOV

physeq_PMOV_mix <- subset_samples(obj_ps_rare, Loc_sec=="Pahute_Mesa" | Loc_sec=="Oasis_Valley" | Loc_sec=="Death_Valley" | Loc_sec=="Rainier_Mesa")
physeq_PMOV_mix <- phyloseq_filter_sample_wise_abund_trim(physeq_PMOV_mix, minabund = 1, rm_zero_OTUs = TRUE)
physeq_PMOV_mix

physeq_AFFCR <- subset_samples(obj_ps_rare, Loc_sec=="Rainier_Mesa" | Loc_sec=="Amargosa_Valley" | Loc_sec=="Death_Valley")
physeq_AFFCR <- phyloseq_filter_sample_wise_abund_trim(physeq_AFFCR, minabund = 1, rm_zero_OTUs = TRUE)
physeq_AFFCR

physeq_AFFCR_mix <- subset_samples(obj_ps_rare, Loc_sec=="Rainier_Mesa" | Loc_sec=="Amargosa_Valley" | Loc_sec=="Death_Valley" | Loc_sec=="Ash_Meadows" | Loc_sec=="Spring_Mountains" | Loc_sec=="Frenchman_and_Yucca_Flat")
physeq_AFFCR_mix <- phyloseq_filter_sample_wise_abund_trim(physeq_AFFCR_mix, minabund = 1, rm_zero_OTUs = TRUE)
physeq_AFFCR_mix

physeq_AM <- subset_samples(obj_ps_rare, Loc_sec=="Ash_Meadows" | Loc_sec=="Spring_Mountains" | Loc_sec=="Frenchman_and_Yucca_Flat")
physeq_AM <- phyloseq_filter_sample_wise_abund_trim(physeq_AM, minabund = 1, rm_zero_OTUs = TRUE)
physeq_AM

physeq_AM_mix <- subset_samples(obj_ps_rare, Loc_sec=="Ash_Meadows" | Loc_sec=="Spring_Mountains" | Loc_sec=="Frenchman_and_Yucca_Flat" | Loc_sec=="Death_Valley")
physeq_AM_mix <- phyloseq_filter_sample_wise_abund_trim(physeq_AM_mix, minabund = 1, rm_zero_OTUs = TRUE)
physeq_AM_mix

physeq_AM_lim <- subset_samples(obj_ps_rare, Loc_sec=="Ash_Meadows" | Loc_sec=="Spring_Mountains")
physeq_AM_lim <- phyloseq_filter_sample_wise_abund_trim(physeq_AM_lim, minabund = 1, rm_zero_OTUs = TRUE)
physeq_AM_lim

## Obtain the metadata and phylogenetic tree for each set
metadata_all <- as.data.frame(as.matrix(sample_data(obj_ps_rare)))
metadata_PMOV <- as.data.frame(as.matrix(sample_data(physeq_PMOV)))
metadata_PMOV_mix <- as.data.frame(as.matrix(sample_data(physeq_PMOV_mix)))
metadata_AFFCR <- as.data.frame(as.matrix(sample_data(physeq_AFFCR)))
metadata_AFFCR_mix <- as.data.frame(as.matrix(sample_data(physeq_AFFCR_mix)))
metadata_AM <- as.data.frame(as.matrix(sample_data(physeq_AM)))
metadata_AM_mix <- as.data.frame(as.matrix(sample_data(physeq_AM_mix)))
metadata_AM_lim <- as.data.frame(as.matrix(sample_data(physeq_AM_lim)))

phy_all <- phy_tree(obj_ps_rare)
phy_PMOV <- phy_tree(physeq_PMOV)
phy_PMOV_mix <- phy_tree(physeq_PMOV_mix)
phy_AFFCR <- phy_tree(physeq_AFFCR)
phy_AFFCR_mix <- phy_tree(physeq_AFFCR_mix)
phy_AM <- phy_tree(physeq_AM)
phy_AM_mix <- phy_tree(physeq_AM_mix)
phy_AM_lim <- phy_tree(physeq_AM_lim)

## Transform to relative abundance
comm_all <- decostand(as.data.frame(t(otu_table(obj_ps_rare))), method = "total")
comm_PMOV <- decostand(as.data.frame(t(otu_table(physeq_PMOV))), method = "total")
comm_PMOV_mix <- decostand(as.data.frame(t(otu_table(physeq_PMOV_mix))), method = "total")
comm_AFFCR <- decostand(as.data.frame(t(otu_table(physeq_AFFCR))), method = "total")
comm_AFFCR_mix <- decostand(as.data.frame(t(otu_table(physeq_AFFCR_mix))), method = "total")
comm_AM <- decostand(as.data.frame(t(otu_table(physeq_AM))), method = "total")
comm_AM_mix <- decostand(as.data.frame(t(otu_table(physeq_AM_mix))), method = "total")
comm_AM_lim <- decostand(as.data.frame(t(otu_table(physeq_AM_lim))), method = "total")

## convert phylogenety to a distance matrix
phy.dist_all <- cophenetic(phy_tree(obj_ps_rare))
phy.dist_PMOV <- cophenetic(phy_PMOV)
phy.dist_PMOV_mix <- cophenetic(phy_PMOV_mix)
phy.dist_AFFCR <- cophenetic(phy_AFFCR)
phy.dist_AFFCR_mix <- cophenetic(phy_AFFCR_mix)
phy.dist_AM <- cophenetic(phy_AM)
phy.dist_AM_mix <- cophenetic(phy_AM_mix)
phy.dist_AM_lim <- cophenetic(phy_AM_lim)

## Prepare the data
phy_objs <- c("PMOV", "PMOV_mix", "AFFCR", "AFFCR_mix", "AM", "AM_mix", "AM_lim", "all")
keep_col <- c("Temp_C", "Depth_sampling_m", "Ca_mgL", "SO4_mgL", "Mg_mgL", "Na_mgL") # significant variables for beta diversity RPCA and variables that describe the water chemistry and rock type 

for(i in phy_objs) {
    # Prepare metadata 
        metadata <- get(paste0("metadata_", i))
        metadata_num <- metadata[keep_col] # select the metadata of interest
        metadata_num <- as.data.frame(apply(metadata_num, 2, as.numeric))
        metadata_num[,2] <- metadata_num[,2]+0.001 #add 0.001 to depth to help with the relative-abundance-weighted mean value
        rownames(metadata_num) <- rownames(metadata)
        assign(paste("metadata_mantel_", i, sep=""), metadata_num) #export variable
    # Sum relative abundances for each OTU
        comm <- get(paste0("comm_", i))
        sum_rel_abund <- as.data.frame(mapply(sum,comm))
        colnames(sum_rel_abund) <- c("sum_rel_abund")
        assign(paste("sum_rel_abund_", i, sep=""), sum_rel_abund) #export variable
    # Print out to check
        print(i)
        print(paste("OTU_table ncol:",ncol(comm)," nrow:", nrow(comm)))
        print(paste("sum_rel_abund ncol:",ncol(sum_rel_abund)," nrow:", nrow(sum_rel_abund)))
        print(paste("metadata ncol:",ncol(metadata_num)," nrow:", nrow(metadata_num)))
}

## calculate relative-abundance–weighted mean value
for(i in phy_objs) {
    for(var in colNames){
        # get for loop variables
            metadata <- get(paste0("metadata_mantel_", i))
            OTU_table <- get(paste0("comm_", i))
            sum_rel_abund <- get(paste0("sum_rel_abund_", i))
        # make matrix 
            colNames <- colnames(metadata) # list of variables to calculate environmental optima
            loop_df <- NULL # create an empty dataframe
            mult_abund_par <- OTU_table*metadata[[var]] # multiply rel. abundances with the parameter value
            mult_abund_par <- as.data.frame(mapply(sum,mult_abund_par)) # sum it up for each OTU
            colnames(mult_abund_par) <- paste(var)
            loop_df[[var]] <- mult_abund_par / sum_rel_abund[match(rownames(mult_abund_par), rownames(sum_rel_abund)),] # divide by the total rel. abundance for each OTU
        # Between-OTU differences environmental optima
            weight_par <- scale(data.frame(loop_df)) # normalized as standard normal deviates
            dist_weight_par <- as.matrix(dist(weight_par, method="euclidean", diag = TRUE, upper = TRUE)) # Among-OTU differences
        assign(paste("dist_weight_par_", i, sep=""), dist_weight_par) #export variable
    }
}

## calculate mantel correlogram (Figure S2)
paste("PMOV", Sys.time(), sep=" ")
mantel_PMOV <- mantel.correlog(dist_weight_par_PMOV, phy.dist_PMOV, r.type="pearson", nperm=999, mult="bonferroni", progressive=TRUE, n.class=50, cutoff=FALSE)

paste("PMOV_mix", Sys.time(), sep=" ")
mantel_PMOV_mix <- mantel.correlog(dist_weight_par_PMOV_mix, phy.dist_PMOV_mix, r.type="pearson", nperm=999, mult="bonferroni", progressive=TRUE, n.class=50, cutoff=FALSE)

paste("AFFCR", Sys.time(), sep=" ")
mantel_AFFCR <- mantel.correlog(dist_weight_par_AFFCR, phy.dist_AFFCR, r.type="pearson", nperm=999, mult="bonferroni", progressive=TRUE, n.class=50, cutoff=FALSE)

paste("AFFCR_mix", Sys.time(), sep=" ")
mantel_AFFCR_mix <- mantel.correlog(dist_weight_par_AFFCR_mix, phy.dist_AFFCR_mix, r.type="pearson", nperm=999, mult="bonferroni", progressive=TRUE, n.class=50, cutoff=FALSE)

paste("AM", Sys.time(), sep=" ")
mantel_AM <- mantel.correlog(dist_weight_par_AM, phy.dist_AM, r.type="pearson", nperm=999, mult="bonferroni", progressive=TRUE, n.class=50, cutoff=FALSE)

paste("AM_mix", Sys.time(), sep=" ")
mantel_AM_mix <- mantel.correlog(dist_weight_par_AM_mix, phy.dist_AM_mix, r.type="pearson", nperm=999, mult="bonferroni", progressive=TRUE, n.class=50, cutoff=FALSE)

paste("AM_lim", Sys.time(), sep=" ")
mantel_AM_lim <- mantel.correlog(dist_weight_par_AM_lim, phy.dist_AM_lim, r.type="pearson", nperm=999, mult="bonferroni", progressive=TRUE, n.class=50, cutoff=FALSE)

paste("all", Sys.time(), sep=" ")
mantel_all <- mantel.correlog(dist_weight_par_all, phy.dist_all, r.type="pearson", nperm=999, mult="bonferroni", progressive=TRUE, n.class=50, cutoff=FALSE)

paste("Done", Sys.time(), sep=" ")

### reformat data and plot mantel correlogram
for(i in phy_objs){
    # get variable
        mantel <- get(paste0("mantel_", i))
    # rescale each bin between 0-1
        mantel_scaled <- as.data.frame(mantel$mantel.res)
        mantel_scaled$class.index <- rescale(mantel_scaled$class.index, to=c(0,1))
        colnames(mantel_scaled) <- c("class.index", "n.dist", "Mantel.cor", "P", "P.adj")
    # plot and save figure
        plot_theme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 1, linetype = "solid"),
            panel.border = element_rect(colour="black", size=1, fill=NA),
            panel.grid.major = element_line(size = 0),
            panel.grid.minor = element_line(size = 0),
            axis.text = element_text(size=20, colour="black"),
            axis.title = element_text(face="bold", size=20),
            legend.position="right",
            legend.key = element_rect(fill = "white"),
            legend.title = element_text(face="bold", size=20),
            legend.text = element_text(size=20))
        plot_guide <- guides(fill = guide_legend(override.aes = list(shape = 21, alpha=1)))
        ggplot(mantel_scaled, aes(x=class.index, y=Mantel.cor)) + 
            geom_point(size=4, shape=21, aes(fill=P.adj <0.05)) +
            geom_line(color="black") +
            labs(y = "Mantel correlation", x="Phylogenetic distance") +
            scale_fill_manual(name = 'Significant\n(p-value < 0.05)', values = setNames(c('black','white'), c(T,F))) +
            geom_hline(yintercept=0, size=.2, linetype = "dashed", color="red") + 
            plot_theme + plot_guide
    save_file_plot <- paste("mantel_correlogram_",i,".pdf", sep="")
    ggsave(save_file_plot, path = determ, scale = 1, width = 10, height = 5, units = c("in"), dpi = 300)
    assign(paste("mantel_scaled_", i, sep=""), mantel_scaled) #export variable
}

## Calculate bNTI
### get otu_table
comm_counts_all <- as.data.frame(otu_table(obj_ps_rare))
comm_counts_PMOV <- otu_table(physeq_PMOV)
comm_counts_PMOV_mix <- as.data.frame(otu_table(physeq_PMOV_mix))
comm_counts_AFFCR <- as.data.frame(otu_table(physeq_AFFCR))
comm_counts_AFFCR_mix <- as.data.frame(otu_table(physeq_AFFCR_mix))
comm_counts_AM_mix <- as.data.frame(otu_table(physeq_AM_mix))
comm_counts_AM_lim <- as.data.frame(otu_table(physeq_AM_lim))

### Calculate bMNTD
for(i in phy_objs){
    # get for loop variables
        phy <- get(paste0("phy_",i))
        comm <- get(paste0("comm_counts_",i))
    # Calculate bMNTD
        phylo = match.phylo.data(phy, comm) # Matching the tree to the rarefied OTU dataset
        bMNTD = as.matrix(comdistnt(t(phylo$data), cophenetic(phylo$phy), abundance.weighted = T)) # Calculating bMNTD for my samples
        bMNTD.rand = array(c(-999), dim = c(ncol(phylo$data), ncol(phylo$data), 999)) # Creating 999 'dummy' matrices to put random results into
    #export variables
        assign(paste("phylo_", i, sep=""), phylo) 
        assign(paste("bMNTD_", i, sep=""), bMNTD) 
        assign(paste("bMNTD.rand_", i, sep=""), bMNTD.rand) 
}

### Performing the calculations on using the OTU table but with randomized taxonomic affiliations
paste("PMOV", Sys.time(), sep=" ")
for(i in 1:999){
  bMNTD.rand_PMOV[,,i] = as.matrix(comdistnt(t(phylo_PMOV$data), taxaShuffle(cophenetic(phylo_PMOV$phy)), abundance.weighted = T, exclude.conspecifics = F))} 

paste("PMOV_mix", Sys.time(), sep=" ")
for(i in 1:999){
  bMNTD.rand_PMOV_mix[,,i] = as.matrix(comdistnt(t(phylo_PMOV_mix$data), taxaShuffle(cophenetic(phylo_PMOV_mix$phy)), abundance.weighted = T, exclude.conspecifics = F))} 

paste("AFFCR", Sys.time(), sep=" ")
for(i in 1:999){
  bMNTD.rand_AFFCR[,,i] = as.matrix(comdistnt(t(phylo_AFFCR$data), taxaShuffle(cophenetic(phylo_AFFCR$phy)), abundance.weighted = T, exclude.conspecifics = F))} 

paste("AFFCR_mix", Sys.time(), sep=" ")
for(i in 1:999){
  bMNTD.rand_AFFCR_mix[,,i] = as.matrix(comdistnt(t(phylo_AFFCR_mix$data), taxaShuffle(cophenetic(phylo_AFFCR_mix$phy)), abundance.weighted = T, exclude.conspecifics = F))} 

paste("AM_lim", Sys.time(), sep=" ")
for(i in 1:999){
  bMNTD.rand_AM_lim[,,i] = as.matrix(comdistnt(t(phylo_AM_lim$data), taxaShuffle(cophenetic(phylo_AM_lim$phy)), abundance.weighted = T, exclude.conspecifics = F))} 

paste("AM_mix", Sys.time(), sep=" ")
for(i in 1:999){
  bMNTD.rand_AM_mix[,,i] = as.matrix(comdistnt(t(phylo_AM_mix$data), taxaShuffle(cophenetic(phylo_AM_mix$phy)), abundance.weighted = T, exclude.conspecifics = F))} 

paste("all", Sys.time(), sep=" ")
for(i in 1:999){
  bMNTD.rand_all[,,i] = as.matrix(comdistnt(t(phylo_all$data), taxaShuffle(cophenetic(phylo_all$phy)), abundance.weighted = T, exclude.conspecifics = F))} 

### Calculate bNTI
for(a in phy_objs){
    g <- get(paste0("phylo_",a))
    assign("tree", g)
    bNTI = matrix(c(NA), nrow = ncol(tree$data), ncol = ncol(tree$data))
    assign(paste0("bNTI_",a),bNTI)
    for(i in 1:(ncol(tree$data)-1)){ 
        for(j in (i+1):ncol(tree$data)){
            bMNTD.rand <- get(paste0("bMNTD.rand_",a))
            bNTI <- get(paste0("bNTI_",a))
            bMNTD <- get(paste0("bMNTD_",a))
                m = bMNTD.rand[j,i,] # Just setting all the randomizations for a given comparison to a matrix
                bNTI[j,i] = ((bMNTD[j,i]-mean(m))/sd(m)) # The bNTI calculation
                assign(paste0("bNTI_",a),bNTI)
        }
    }
    bNTI <- get(paste0("bNTI_",a))
    rownames(bNTI) = colnames(tree$data)
    colnames(bNTI) = colnames(tree$data)
    write.csv(bNTI, paste0("bNTI_TotalCounts_",a,".csv"), quote = F) #save csv file
    assign(paste0("bNTI_",a),bNTI)
}

### reformat data and plot 
for(i in phy_objs){
    # get for loop variables
        bNTI <- get(paste0("bNTI_",i))
        metadata <- get(paste0("metadata_",i))
    # melt the bNTI table
        bNTI[upper.tri(bNTI)] = t(bNTI)[upper.tri(bNTI)]
        print(paste(i,dim(bNTI),sep=" "))
        bNTI.melt = reshape2::melt(bNTI)
        bNTI.melt.half <- bNTI.melt %>% distinct(value, .keep_all=TRUE)
        write.csv(bNTI.melt.half, paste0("bNTI_TotalCounts_melt_",i,".csv"), quote = F) #save csv file
    #define each triangle of the plot matric and the diagonal (mi.ids)
        diag(bNTI) = NA # Self comparisons as zero could inflate/deflate differences
        heat = melt(as.matrix(bNTI))
        names(heat) <- c("M1", "M2", "bNTI") 
        mi.ids <- subset(heat, M1 == M2)
        mi.lower <- subset(heat[lower.tri(bNTI),], M1 != M2)
        mi.upper <- subset(heat[upper.tri(bNTI),], M1 != M2)
        plot_data <- rbind(mi.ids, mi.lower) 
    # plot heatmap
        plot_theme <- theme(panel.background = element_blank(),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_blank(),
            axis.title.x = element_blank(),
            axis.text.y = element_text(size=12, colour="black", hjust=1.2),
            axis.text.x = element_text(angle = 45, vjust = 1.2, hjust=1, size=12, colour="black"),
            axis.ticks = element_blank(),
            legend.position="right",
            legend.key = element_rect(fill = "white"),
            legend.title = element_text(face="bold", size=12),
            legend.text = element_text(size=12))
        plot_guide <- guides(fill = guide_colourbar(frame.linewidth = 1, frame.colour = "black", ticks = TRUE, ticks.colour = "black", ticks.linewidth = 1))
        p1 <- ggplot(plot_data, aes(M1, M2, fill=bNTI)) + 
            geom_tile(colour="black", size=0.3) +
            geom_text(data=mi.lower, aes(label=ifelse(abs(bNTI)>2, paste(round(bNTI,1)), "")), size=2)
        meas <- as.character(unique(plot_data$M2))
        p2 <- p1 + scale_colour_identity() + 
            scale_fill_gradient2(low = "#b86902", mid = "#F9F9F9", high = "#050b61", na.value = "#666666") +
            scale_x_discrete(limits=meas[length(meas):1]) + #flip the x axis 
            scale_y_discrete(limits=meas) +
            plot_theme + plot_guide
        save_file_plot <- paste("bNTI_heatmap_",i,".pdf", sep="")
        ggsave(save_file_plot, path = determ, scale = 1, width = 10, height = 7, units = c("in"), dpi = 300)
    #plot box plot
        plot_theme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 1, linetype = "solid"),
            panel.border = element_rect(colour="black", size=1, fill=NA),
            panel.grid.major = element_line(size = 0),
            panel.grid.minor = element_line(size = 0),
            axis.text = element_text(size=20, colour="black"),
            axis.title = element_text(face="bold", size=20),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            legend.position="right",
            legend.key = element_rect(fill = "white"),
            legend.title = element_text(face="bold", size=20),
            legend.text = element_text(size=20))
        bNTI_plots_box <- as.data.frame(bNTI)
        bNTI_plots_box$Sample_abbrev <- metadata$Sample_abbrev[match(rownames(metadata), rownames(bNTI_plots_box))]
        bNTI_plots_box$Loc_sec <- metadata$Loc_sec[match(rownames(metadata), rownames(bNTI_plots_box))]
        bNTI.melt.box = melt(data = bNTI_plots_box, id.vars = c("Sample_abbrev", "Loc_sec"))
        assign(paste0("bNTI.melt.box_",a),bNTI.melt.box)
        bNTI_boxplot <- ggplot(data = bNTI.melt.box, aes(x = Sample_abbrev, y = value))+
            stat_boxplot(geom = "errorbar", width = 0.2) +
            geom_boxplot(outlier.size=0, outlier.colour="white", fill="lightgrey")+
            geom_jitter(size=2, position=position_jitter(0.2), alpha=0.5, fill="lightgrey", shape=21, stroke=1, colour="black") +
            geom_hline(yintercept = c(2,-2), color = "red", linetype = 2, size=1) +
            labs(y="bNTI") +
            plot_theme
        save_file_plot <- paste("bNTI_boxplot_",i,".pdf", sep="")
        ggsave(save_file_plot, path = determ, scale = 1, width = 20, height = 7, units = c("in"), dpi = 300)
}

## Calculate Raup-Crick
### Setting up the data to mesh with the Stegen et al. code
spXsite = comm
spXsite[spXsite>0] = spXsite[spXsite>0]+1 # Adding one to my RPKM value to generate an "RPKM+1" stat. This is necessary to work around Raup-Crick
head(spXsite)

rc.reps = 9999
registerDoMC(cores = 1)

### Count number of sites and total species richness across all plots (gamma)
n_sites = nrow(spXsite)
n_sites
gamma = ncol(spXsite)
gamma

### Build a site by site matrix for the results, with the names of the sites in the row and col names:
results = matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(spXsite), row.names(spXsite)))

### Make the spXsite matrix into a new, pres/abs. matrix:
spXsite.inc = ceiling(spXsite/max(spXsite))

### Create an occurrence vector- used to give more weight to widely distributed species in the null model
occur = apply(spXsite.inc, MARGIN=2, FUN=sum)

### Create an abundance vector- used to give more weight to abundant species in the second step of the null model
abundance = apply(spXsite, MARGIN=2, FUN=sum)

for(null.one in 1:(nrow(spXsite)-1)){
  for(null.two in (null.one+1):nrow(spXsite)){
    null_bray_curtis<-NULL
    null_bray_curtis = foreach(i=1:rc.reps, .packages = c("vegan","picante")) %dopar% {
      # Generates two empty communities of size gamma
          com1<-rep(0,gamma)
          com2<-rep(0,gamma)
      # Add observed number of species to com1, weighting by species occurrence frequencies
          com1[sample(1:gamma, sum(spXsite.inc[null.one,]), replace=FALSE, prob=occur)]<-1
          com1.samp.sp = sample(which(com1>0), (sum(spXsite[null.one,])-sum(com1)), replace=TRUE, prob=abundance[which(com1>0)]);
          com1.samp.sp = cbind(com1.samp.sp,1);
          com1.sp.counts = as.data.frame(tapply(com1.samp.sp[,2], com1.samp.sp[,1], FUN=sum)); colnames(com1.sp.counts) = 'counts'; # head(com1.sp.counts);
          com1.sp.counts$sp = as.numeric(rownames(com1.sp.counts));
          com1[com1.sp.counts$sp] = com1[com1.sp.counts$sp] + com1.sp.counts$counts;
          #sum(com1) - sum(spXsite[null.one,]); # This should be zero if everything worked properly
          rm('com1.samp.sp','com1.sp.counts');
      # Again for com2
          com2[sample(1:gamma, sum(spXsite.inc[null.two,]), replace=FALSE, prob=occur)]<-1
          com2.samp.sp = sample(which(com2>0), (sum(spXsite[null.two,])-sum(com2)), replace=TRUE, prob=abundance[which(com2>0)]);
          com2.samp.sp = cbind(com2.samp.sp,1);
          com2.sp.counts = as.data.frame(tapply(com2.samp.sp[,2], com2.samp.sp[,1], FUN=sum)); colnames(com2.sp.counts) = 'counts'; # head(com2.sp.counts);
          com2.sp.counts$sp = as.numeric(rownames(com2.sp.counts));
          com2[com2.sp.counts$sp] = com2[com2.sp.counts$sp] + com2.sp.counts$counts;
          #sum(com2) - sum(spXsite[null.two,]); # This should be zero if everything worked properly
          rm('com2.samp.sp','com2.sp.counts');
      null.spXsite = rbind(com1,com2); # Null.spXsite
      # Calculates the null Bray-Curtis
        null_bray_curtis[i] = vegdist(null.spXsite, method='bray');
    }; # End of the null loop
    # Unlisting the parallel list
        null_bray_curtis = unlist(null_bray_curtis)
    # Calculates the observed Bray-Curtis
        obs.bray = vegdist(spXsite[c(null.one,null.two),], method='bray');
    # How many null observations is the observed value tied with?
        num_exact_matching_in_null = sum(null_bray_curtis==obs.bray);
    # How many null values are smaller than the observed *dissimilarity*?
        num_less_than_in_null = sum(null_bray_curtis<obs.bray);
    rc = ((num_less_than_in_null + (num_exact_matching_in_null)/2)/rc.reps) # This variation of rc splits ties
    rc = (rc-.5)*2 # Adjusts the range of the  Raup-Crick caclulation to -1 to 1
    results[null.two,null.one] = round(rc,digits=2); # Stores rc into the results matrix
    print(c(null.one,null.two,date())); # Keeps track of position
  }; # End of inner loop
}; # End of outer loop

rc.results = as.dist(results) # Converts results into a distance matrix
write.csv(as.matrix(rc.results), "S3_RCBC_TotalCounts.csv", quote = F)
rm('spXsite.inc')

## Organize the data
rcbc <- read.csv("S3_RCBC_TotalCounts.csv", row.names=1)

bNTI_rcbc <- bNTI_all
rcbc2 <- rcbc

if(identical(x = row.names(rcbc2), y = row.names(bNTI_rcbc)) == FALSE){
    w = which(!(row.names(bNTI_rcbc) %in% row.names(rcbc2)) %in% TRUE)
    rcbc2 = rbind(rcbc2[1:(w-1),], NA, rcbc2[w:length(rcbc2[,1]),])
    rcbc2 = cbind(rcbc2[,1:(w-1)], NA, rcbc2[,w:length(rcbc2[1,])])}

rcbc2[abs(bNTI_rcbc) > 2] = NA # Replacing RCBC values where the bNTI is significant
diag(rcbc2) = NA # Self comparisons as zero could inflate/deflate differences
heat = melt(as.matrix(rcbc2))
names(heat) <- c("M1", "M2", "RCBC") 
mi.ids <- subset(heat, M1 == M2)
mi.lower <- subset(heat[lower.tri(rcbc2),], M1 != M2)
mi.upper <- subset(heat[upper.tri(rcbc2),], M1 != M2)
plot_data <- rbind(mi.ids, mi.lower) 

## Plot heatmap RCBC
plot_theme <- theme(panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size=12, colour="black", hjust=1.2),
    axis.text.x = element_text(angle = 45, vjust = 1.2, hjust=1, size=12, colour="black"),
    axis.ticks = element_blank(),
    legend.position="right",
    legend.key = element_rect(fill = "white"),
    legend.title = element_text(face="bold", size=12),
    legend.text = element_text(size=12))
plot_guide <- guides(fill = guide_colourbar(frame.linewidth = 1, frame.colour = "black", ticks = TRUE, ticks.colour = "black", ticks.linewidth = 1))

p1 <- ggplot(plot_data, aes(M1, M2, fill=RCBC)) + 
    geom_tile(colour="black", size=0.3) +
    geom_text(data=mi.lower, aes(label=ifelse(abs(RCBC)>0.95, paste(round(RCBC,2)), "")), size=2) #+ 
    #geom_text(data=mi.ids, aes(label=M2, colour="black"), position = position_nudge(x = 1.5), size=4)
             
meas <- as.character(unique(plot_data$M2))
p2 <- p1 + scale_colour_identity() + 
    scale_fill_gradient2(low = "#ffb57d", mid = "#F9F9F9", high = "#9984d4", na.value = "#666666") +
    scale_x_discrete(limits=meas[length(meas):1]) + #flip the x axis 
    scale_y_discrete(limits=meas) +
    plot_theme + plot_guide
p2

save_file_plot <- paste("RCBC_heatmap.pdf", sep="")
ggsave(save_file_plot, path = determ, scale = 1, width = 10, height = 7, units = c("in"), dpi = 300)

write.csv(plot_data, "RCBC_bNTIsigOnly.csv", quote = F)

## plot boxplot RCBC
plot_theme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 1, linetype = "solid"),
    panel.border = element_rect(colour="black", size=1, fill=NA),
    strip.background=element_rect(fill='white', colour='white', size = 0),
    strip.text = element_text(face="bold", size=20),
    panel.spacing.x=unit(0.5, "lines"),
    panel.grid.major = element_line(size = 0),
    panel.grid.minor = element_line(size = 0),
    axis.text = element_text(size=20, colour="black"),
    axis.title = element_text(face="bold", size=20),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
    legend.position="right",
    legend.key = element_rect(fill = "white"),
    legend.title = element_text(face="bold", size=20),
    legend.text = element_text(size=20))
plot_guide <- guides(fill = guide_legend(order=2, override.aes = list(shape = 21, alpha=1)),
    size = guide_legend(order=3),
    shape = guide_legend(order=1, override.aes = list(size = 5, color="black", alpha=1)))
plot_nomargins_y <- scale_y_continuous(expand = expansion(mult = c(0, 0)), limits=c(-1.2,1.2), breaks=c(-1, 0, 1))

rcbc2_plots_box <- as.data.frame(rcbc2)
rcbc2_plots_box$Sample_abbrev <- metadata$Sample_abbrev[match(rownames(metadata), rownames(rcbc2_plots_box))]
rcbc2_plots_box$Loc_DV3 <- metadata$Loc_DV3[match(rownames(metadata), rownames(rcbc2_plots_box))]
rcbc2_plots_box$Loc_sec <- metadata$Loc_sec[match(rownames(metadata), rownames(rcbc2_plots_box))]
rcbc2.melt.box = melt(data = rcbc2_plots_box, id.vars = c("Sample_abbrev", "Loc_DV3", "Loc_sec"))

rcbc_boxplot <- ggplot(rcbc2.melt.box, aes(x = Sample_abbrev, y = value))+
    stat_boxplot(geom = "errorbar", width = 0.2) +
    geom_boxplot(outlier.size=0, outlier.color="white", aes(fill=Loc_DV3))+
    scale_fill_manual(name = "Groundwater\nFlow Basin", values = c("AFFCR" = "#74c69d", "AM" = "#f4acb7", "PMOV"="lightgrey")) +
    geom_jitter(size=2, position=position_jitter(0.2), alpha=0.5, fill="lightgrey", shape=21, stroke=1, colour="black") +
    geom_hline(yintercept = c(0.95,-0.95), color = "red", linetype = 2, size=1) +
    labs(y="RCBC") +
    facet_grid(.~Loc_DV3, scales="free")+
    plot_theme + plot_guide + plot_nomargins_y
save_file_plot <- paste("RCBC_boxplot.pdf", sep="")
ggsave(save_file_plot, path = determ, scale = 1, width = 20, height = 7, units = c("in"), dpi = 300)

## plot both RCBC and bNTI
ggarrange(bNTI_boxplot_all + theme(legend.position="none"), rcbc_boxplot + theme(legend.position="none"), ncol = 1, nrow = 2)

save_file_plot <- paste("combo_RCBC_bNTI.pdf", sep="")
ggsave(save_file_plot, path = determ, scale = 1, width = 20, height = 15, units = c("in"), dpi = 300)

## Export data to view on excel and add interpretation column
combo <- cbind(rcbc2.melt.box, bNTI.melt.box)
write.csv(combo, paste(determ,"/combo_melted.csv", sep=""), quote = F)

## Plot combined heatmap
### Create symmetric matrix from long form datatables
#### Get important columns
bNTI_heatmap_data <- bNTI.melt.box[,c("Sample_abbrev", "variable", "value")]
rcbc_heatmap_data <- rcbc2.melt.box[,c("Sample_abbrev", "variable", "value")]

#### Make undirected so that graph matrix will be symmetric
bNTI_heatmap_data <- graph.data.frame(bNTI_heatmap_data, directed=FALSE)
rcbc_heatmap_data <- graph.data.frame(rcbc_heatmap_data, directed=FALSE)

#### add value as a weight attribute
bNTI_heatmap_data <- get.adjacency(bNTI_heatmap_data, attr="value", sparse=FALSE)
rcbc_heatmap_data <- get.adjacency(rcbc_heatmap_data, attr="value", sparse=FALSE)

head(bNTI_heatmap_data)
head(rcbc_heatmap_data)

### Get RCBC values with abs bNTI < 2
combo <- cbind(rcbc2.melt.box, bNTI.melt.box)
colnames(combo) <- c("Sample_abbrev", "Loc_DV3", "Loc_sec", "var_RCBC", "RCBC", "Sample_abbrev_bNTI", "Loc_DV3_bNTI", "Loc_sec_bNTI", "var_bNTI", "bNTI")
combo$RCBC_lim <- combo$RCBC
combo_rcbc_lim <- within(combo, RCBC_lim[abs(bNTI) > 2] <- '2')
head(combo_rcbc_lim)
write.csv(combo_rcbc_lim, paste(determ,"/combo_melted_rcbc_lim.csv", sep=""), quote = F)

#### create the matrix
rcbc_heatmap_data_lim <- combo_rcbc_lim[, c("Sample_abbrev", "var_RCBC", "RCBC_lim")]

#### Make undirected so that graph matrix will be symmetric
rcbc_heatmap_data_lim <- graph.data.frame(rcbc_heatmap_data_lim, directed=FALSE)

#### add value as a weight attribute
rcbc_heatmap_data_lim <- get.adjacency(rcbc_heatmap_data_lim, attr="RCBC_lim", sparse=FALSE)

head(rcbc_heatmap_data_lim)

## Sort rows and columns
bNTI_heatmap_data <- bNTI_heatmap_data[, c('AM1' , 'AM2' , 'AM3' , 'AM4' , 'AM5' , 'AM6' , 'AM7', 'SM1' , 'SM2' , 'SM3' , 'SM4', 'YF1' , 'YF2' , 'YF3', 'FF1' , 'FF2', 'RM1' , 'RM2', 'AV1' , 'AV2' , 'AV3', 'DV1' , 'DV2' , 'DV3', 'OV1', 'OV2', 'PM10', 'PM5' , 'PM16' , 'PM15' , 'PM7' , 'PM4' , 'PM8' , 'PM1' , 'PM2' , 'PM3' , 'PM6' , 'PM9' ,  'PM11' , 'PM12' , 'PM13' , 'PM14')]

bNTI_heatmap_data <- bNTI_heatmap_data[c('AM1' , 'AM2' , 'AM3' , 'AM4' , 'AM5' , 'AM6' , 'AM7', 'SM1' , 'SM2' , 'SM3' , 'SM4', 'YF1' , 'YF2' , 'YF3', 'FF1' , 'FF2', 'RM1' , 'RM2', 'AV1' , 'AV2' , 'AV3', 'DV1' , 'DV2' , 'DV3', 'OV1', 'OV2', 'PM10', 'PM5' , 'PM16' , 'PM15' , 'PM7' , 'PM4' , 'PM8' , 'PM1' , 'PM2' , 'PM3' , 'PM6' , 'PM9' ,  'PM11' , 'PM12' , 'PM13' , 'PM14'), ]

colnames(bNTI_heatmap_data)
rownames(bNTI_heatmap_data)

rcbc_heatmap_data_lim <- rcbc_heatmap_data_lim[, c('AM1' , 'AM2' , 'AM3' , 'AM4' , 'AM5' , 'AM6' , 'AM7', 'SM1' , 'SM2' , 'SM3' , 'SM4', 'YF1' , 'YF2' , 'YF3', 'FF1' , 'FF2', 'RM1' , 'RM2', 'AV1' , 'AV2' , 'AV3', 'DV1' , 'DV2' , 'DV3', 'OV1', 'OV2', 'PM10', 'PM5' , 'PM16' , 'PM15' , 'PM7' , 'PM4' , 'PM8' , 'PM1' , 'PM2' , 'PM3' , 'PM6' , 'PM9' ,  'PM11' , 'PM12' , 'PM13' , 'PM14')]

rcbc_heatmap_data_lim <- rcbc_heatmap_data_lim[c('AM1' , 'AM2' , 'AM3' , 'AM4' , 'AM5' , 'AM6' , 'AM7', 'SM1' , 'SM2' , 'SM3' , 'SM4', 'YF1' , 'YF2' , 'YF3', 'FF1' , 'FF2', 'RM1' , 'RM2', 'AV1' , 'AV2' , 'AV3', 'DV1' , 'DV2' , 'DV3', 'OV1', 'OV2', 'PM10', 'PM5' , 'PM16' , 'PM15' , 'PM7' , 'PM4' , 'PM8' , 'PM1' , 'PM2' , 'PM3' , 'PM6' , 'PM9' ,  'PM11' , 'PM12' , 'PM13' , 'PM14'), ]

colnames(rcbc_heatmap_data_lim)
rownames(rcbc_heatmap_data_lim)

## only lower or upper triangle
### bNTI upper triangle
bNTI_heatmap_data_upper <- bNTI_heatmap_data %>% replace_lower_triangle(by = NA)
rownames(bNTI_heatmap_data_upper) <- bNTI_heatmap_data_upper$rowname
bNTI_heatmap_data_upper <- as.matrix(bNTI_heatmap_data_upper[,-1])
head(bNTI_heatmap_data_upper)

### RCBC lower triangle
rcbc_heatmap_data_lower <- rcbc_heatmap_data_lim %>% replace_upper_triangle(by = NA)
rownames(rcbc_heatmap_data_lower) <- rcbc_heatmap_data_lower$rowname
rcbc_heatmap_data_lower <- as.matrix(rcbc_heatmap_data_lower[,-1])
rcbc_heatmap_data_lower <- suppressWarnings(apply(rcbc_heatmap_data_lower, 2 ,as.numeric))
rownames(rcbc_heatmap_data_lower) <- colnames(rcbc_heatmap_data_lower)
head(rcbc_heatmap_data_lower)

## Plot bNTI heatmap
### Annotate columns
annotation_col = data.frame(
    Location = DATA_PHYLOSEQ_FIXED$Loc_sec,
    Basin = DATA_PHYLOSEQ_FIXED$Loc_DV3)
rownames(annotation_col) = rownames(DATA_PHYLOSEQ_FIXED)

ann_colors = list(
    Location = c(Ash_Meadows = '#e50006', Amargosa_Valley = '#043478', Death_Valley = '#38ab31', Frenchman_and_Yucca_Flat = '#0f87a5', Oasis_Valley = '#9ea8a7', Pahute_Mesa = '#7f468d', Rainier_Mesa = '#fb9e7f', Spring_Mountains = '#9a0020'),
    Basin = c(AFFCR = "#478432", AM = "#a5006b", PMOV = "#5c2d19"))

### Color legend gradient (bNTI)
bk1 <- c(seq(-4.52,-1.9,by=0.1),-1.999)
bk2 <- c(-2.001,seq(-2.1,1.9,by=0.1),1.999)
bk3 <- c(2.000, seq(2.001,8.76,by=0.1))
bk <- c(bk1,bk2,bk3)  #combine the break limits for purpose of graphing

my_palette <- c(colorRampPalette(colors = c("#184e77", "#48cae4"))(n = length(bk1)-1),
              "#283618", "#283618",
              c(colorRampPalette(colors = c("#ffc9b9", "#fefee3", "#f4e285"))(n = length(bk2)-1)),
               "#f4e285", "#aad576",
               c(colorRampPalette(colors = c("#aad576", "#386641"))(n = length(bk3)-1)))

### bNTI heatmap
heatmap_bNTI <- pheatmap(bNTI_heatmap_data_upper,
         color = my_palette, breaks = bk,
         na_col = c("white"),
         border_color = "black",
         show_colnames = TRUE,
         show_rownames = TRUE,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         drop_levels = TRUE,
         fontsize = 14,
         cluster_cols = FALSE, # mat_cluster_cols
         cluster_rows = FALSE, #mat_cluster_rows
         main = "bNTI",
         gaps_col = cumsum(c(11,13,18)),
         gaps_row = cumsum(c(11,13,18)), legend=T)

heatmap_bNTI

save_file_plot <- paste("heatmap_bNTI.svg", sep="")
ggsave(save_file_plot, heatmap_bNTI, path = determ, scale = 1, width = 15, height = 10, units = c("in"), dpi = 300)

## Plot RCBC heatmap
### Annotate columns
annotation_col = data.frame(
    Location = DATA_PHYLOSEQ_FIXED$Loc_sec,
    Basin = DATA_PHYLOSEQ_FIXED$Loc_DV3)
rownames(annotation_col) = rownames(DATA_PHYLOSEQ_FIXED)

ann_colors = list(
    Location = c(Ash_Meadows = '#e50006', Amargosa_Valley = '#043478', Death_Valley = '#38ab31', Frenchman_and_Yucca_Flat = '#0f87a5', Oasis_Valley = '#9ea8a7', Pahute_Mesa = '#7f468d', Rainier_Mesa = '#fb9e7f', Spring_Mountains = '#9a0020'),
    Basin = c(AFFCR = "#478432", AM = "#a5006b", PMOV = "#5c2d19"))

### Color legend gradient (RCBC)
bk1 <- c(seq(-1, -0.94,by=0.1),-0.9499)
bk2 <- c(-0.95,seq(-0.9501,0.94,by=0.3),0.9499)
bk3 <- c(0.95, seq(0.9501,1,by=0.1),1.001)
bk4 <- c(1.002, seq(1.003,2,by=0.1))
bk <- c(bk1,bk2,bk3,bk4)  #combine the break limits for purpose of graphing

my_palette <- c(colorRampPalette(colors = c("#a4133c", "#ff8fa3"))(n = length(bk1)-1),
              "#ff8fa3", "#ffc9b9",
              c(colorRampPalette(colors = c("#ffc9b9", "#fefee3", "#f4e285"))(n = length(bk2)-1)),
               "#f4e285", "#c77dff",
               c(colorRampPalette(colors = c("#c77dff", "#3c096c"))(n = length(bk3)-1)),
               "#d3d3d3","#d3d3d3",
               c(colorRampPalette(colors = c("#d3d3d3", "#d3d3d3"))(n = length(bk4)-1)))

### RCBC heatmap
heatmap_RCBC <- pheatmap(rcbc_heatmap_data_lower,
         color = my_palette, breaks = bk,
         na_col = c("white"),
         border_color = "black",
         show_colnames = TRUE,
         show_rownames = TRUE,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         drop_levels = TRUE,
         fontsize = 14,
         cluster_cols = FALSE, # mat_cluster_cols
         cluster_rows = FALSE, #mat_cluster_rows
         main = "RCBC",
         gaps_col = cumsum(c(11,13,18)),
         gaps_row = cumsum(c(11,13,18)), legend=T)

heatmap_RCBC

save_file_plot <- paste("heatmap_RCBC.svg", sep="")
ggsave(save_file_plot, heatmap_RCBC, path = determ, scale = 1, width = 15, height = 10, units = c("in"), dpi = 300)

# Final plot (Figure 6), manually combine the upper and lower triangles and make the legends nice


## Plot ratios
### Import data
between <- read_excel("combo_melted.xlsx", sheet="between")
head(between)

all1 <- group_by(between, Loc_DV3, Loc_sec, Combo_interp) %>%
                arrange(Combo_interp) %>% 
                summarise(observ = n()) %>% 
                mutate(freq = observ / sum(observ))

### fix cluster column
all1$Loc_sec <- gsub('_',' ', all1$Loc_sec)
all1$Loc_sec <- gsub('and','&', all1$Loc_sec)
all1$Loc_sec <- factor(all1$Loc_sec, ordered = TRUE, levels = c("Oasis Valley", "Pahute Mesa", "Rainier Mesa", "Frenchman & Yucca Flat", "Amargosa Valley", "Ash Meadows", "Spring Mountains", "Death Valley"))
all1$Combo_interp <- factor(all1$Combo_interp, ordered = TRUE, levels = c("Variable Selection", "Homogeneous Selection", "Homogenizing Dispersal", "Dispersal Limitation + Drift", "Undominated"))

### Plot box plot for between communities processes (Figure 6 inset)
plot_theme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 1, linetype = "solid"),
    panel.border = element_rect(colour="black", size=1, fill=NA),
    strip.background=element_rect(fill='white', colour='white', size = 0),
    strip.text = element_text(face="bold", size=20),
    panel.spacing.x=unit(0.5, "lines"),
    panel.grid.major = element_line(size = 0),
    panel.grid.minor = element_line(size = 0),
    axis.text = element_text(size=15, colour="black"),
    axis.title = element_text(face="bold", size=20),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
    legend.position="right",
    legend.key = element_rect(fill = "white"),
    legend.title = element_text(face="bold", size=20),
    legend.text = element_text(size=20))
plot_guides <- guides(colour=FALSE, fill = guide_legend(ncol=1))
plot_nomargins_y <- scale_y_continuous(expand = expansion(mult = c(0, 0)), labels = function(x) paste0(x*100, "%"))
plot_nomargins_x <- scale_x_discrete(expand = expansion(mult = c(0, 0)))

all_between <- ggplot(data=all1, aes(x=Loc_sec, y=freq, fill=Combo_interp)) +
    scale_fill_manual(values=c("Variable Selection"="#80b92c", "Homogeneous Selection"="#2f6f43", "Homogenizing Dispersal"="#f26a8d", "Dispersal Limitation + Drift"="#dd2d4a", "Undominated"="#d6d2d2")) + 
    labs(y = "Ratio", fill = "Community Process", x = "Site") +
    geom_bar(aes(), stat="identity", position="fill", color="black", width=0.9) +
    plot_nomargins_y + plot_nomargins_x + plot_theme + plot_guides
save_file_plot <- paste("ratio_assembly_processes.pdf", sep="") 
ggsave(save_file_plot, path=determ, scale = 1, width = 10, height = 8, units = c("in"), dpi = 300)
