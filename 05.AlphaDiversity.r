#  05. Alpha Diversity
## Figure 3, Figure S6

## Set Path
alpha <- file.path(paste(path_phy, "/alpha_diversity", sep=""))
dir.create(alpha, showWarnings=FALSE)
setwd(alpha)

## obtain number of OTUs per sample
p <- plot_richness(obj_ps, x = "Sample_abbrev", measures = "Observed")
OTUs <- data.frame("OTUs" = p$data$value)
rownames(OTUs) <- rownames(sample_data(obj_ps))
write.csv(OTUs, "OTU_sums_per_sample.csv")

p

## Calculate alpha diversity
alpha.div <- microbiome::alpha(obj_ps, index = "all")
alpha.div.phyloseq <- estimate_richness(obj_ps, split = TRUE, measures = NULL)
write.table(alpha.div,paste("alpha.div.genusGlom.pkgmicrobiome.csv", sep=""), sep=",", row.names = T)
write.table(alpha.div.phyloseq,paste("alpha.div.genusGlom.pkgphyloseq.csv", sep=""), sep=",", row.names = T)

## Calculate Faith's PD
obj_ps.pd <- pd(t(otu_table(obj_ps)), phy_tree(filt_tree), include.root=T) 
head(obj_ps.pd)
write.table(obj_ps.pd,paste("alpha.div.genusGlom.faithsPD.csv", sep=""), sep=",", row.names = T)
alpha.div.phyloseq$PD <- obj_ps.pd$PD[match(rownames(alpha.div.phyloseq), rownames(obj_ps.pd))] #add Faith's PD to alpha.div.phyloseq

## Plot preliminary boxplots
parameters <- c("Loc_sec","Well_spring","rock_type","Piper_group3","Sequence_batch", "Loc_DV3")

func_richness <- llply(as.list(parameters), function(para, x) {
    plot_alpha <- plot_richness(physeq = x, x = para) + 
                    geom_boxplot(notch=FALSE) + 
                    geom_jitter(shape=16, position=position_jitter(0.2)) + 
                    plot_theme
    save_file_plot <- paste("alpha_div.barplot.genusGlom.", para,".pdf", sep="")
    ggsave(save_file_plot, path = alpha, scale = 1, width = 15, height = 6, units = c("in"), dpi = 300)
}, obj_ps)

## Data formatting 
metadata <- as.data.frame(as.matrix(sample_data(obj_ps)))
alpha.div.phyloseq$Read_count <- reads_per_sample$Read_count[match(rownames(alpha.div.phyloseq), rownames(reads_per_sample))] #Add read count per sample
alpha.div.phyloseq$Sample_abbrev <- metadata$Sample_abbrev[match(rownames(alpha.div.phyloseq), rownames(metadata))] # add Sample_abbrev column to 'alpha.div.phyloseq'
alpha.div.metadata <- as.data.frame(merge(metadata, alpha.div.phyloseq, by = "Sample_abbrev")) # merge these two data frames into one

## Check library size effects on alpha diversity
p1 <- ggscatter(alpha.div.metadata, "Shannon", "Read_count", color="Sequence_batch") + stat_cor(method = "pearson")
p2 <- ggscatter(alpha.div.metadata, "Simpson", "Read_count", color="Sequence_batch") + stat_cor(method = "pearson")
p3 <- ggscatter(alpha.div.metadata, "PD", "Read_count", color="Sequence_batch") + stat_cor(method = "pearson")
p4 <- ggscatter(alpha.div.metadata, "Observed", "Read_count", color="Sequence_batch") +
        stat_cor(method = "pearson", label.x = 100, label.y = 50000)

ggarrange(p4, p3, p1, p2, ncol = 4, nrow = 1)
ggsave("alpha_div.library_size.pdf", path = alpha, scale = 1, width = 10, height = 3, units = c("in"), dpi = 300)

## Plot correlation matrix
alpha.div.metadata.corr <- alpha.div.metadata[,13:38]
alpha.div.metadata.corr[,1:26] <- sapply(alpha.div.metadata.corr[,1:26],as.numeric)
head(alpha.div.metadata.corr)
colnames(alpha.div.metadata.corr)
M<-cor(alpha.div.metadata.corr, method="s")

### computing p-values
cor.mtest <- function(mat, ...) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            tmp <- cor.test(mat[, i], mat[, j], ...)
            p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
        }
    }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
p.mat <- cor.mtest(alpha.div.metadata.corr)

### plot
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

pdf(file="alpha_div.corrplot.pdf", width=20, height=20)
corrplot(M, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
         )
dev.off()

alpha.div.metadata.corr.refined <- alpha.div.metadata.corr[c("Ca_mgL", "Temp_C", "Depth_sampling_m", "Shannon", "Simpson", "PD")]
pdf(file="alpha_div.corrplot.refined.pdf", width=10, height=10)
chart.Correlation(alpha.div.metadata.corr.refined, histogram=TRUE, pch=19, method="spearman")
dev.off()

chart.Correlation(alpha.div.metadata.corr.refined, histogram=TRUE, pch=19, method="spearman")

res2 <- rcorr(as.matrix(alpha.div.metadata.corr.refined), type="spearman")
res2$r
res2$P

## Data organizing
alpha.div.metadata2$tritium_BqL_fct <- cut(x = alpha.div.metadata$tritium_BqL, breaks = c(0, 100, Inf))
levels(alpha.div.metadata2$tritium_BqL_fct) <- c("Low", "High")
rownames(alpha.div.metadata2) <- alpha.div.metadata2$Sample_abbrev
head(alpha.div.metadata2)

## Statistical tests (Repeat 1-6 for every variable combinations) (Table S6)
### 1. Anova
aov.obs = aov(Observed ~ Loc_sec, data=alpha.div.metadata2)
summary.lm(aov.obs)

### 2. Extract the residuals & Run Shapiro-Wilk test
aov_residuals <- residuals(object = aov.obs)
shapiro.test(x = aov_residuals)

### 3. Check homogeneity
leveneTest(Observed ~ Loc_sec, data = alpha.div.metadata2)

### 4. Kruskal-Wallis test
kruskal.test(Observed ~ Loc_sec, data=alpha.div.metadata2)

### 5. To do all the pairwise comparisons between groups and correct for multiple comparisons
pairwise.wilcox.test(alpha.div.metadata2$Observed, alpha.div.metadata2$Loc_sec, p.adjust.method="fdr", exact=FALSE, paired=FALSE)

### 6. Plot
par(mfrow=c(2,2))
plot(aov.obs)

### Repeat the above #1-6 without OV1 and OV2
alpha.div.metadata2_noOV <- subset(alpha.div.metadata2, Sample_abbrev!="OV1")
alpha.div.metadata2_noOV <- subset(alpha.div.metadata2_noOV, Sample_abbrev!="OV2")

## More data organizing for making plots
plot_metadata <- alpha.div.metadata2
plot_metadata$Loc_sec <- gsub('_',' ', plot_metadata$Loc_sec)
plot_metadata$Loc_sec <- gsub('and','&', plot_metadata$Loc_sec)
plot_metadata$Loc_sec <- factor(plot_metadata$Loc_sec, ordered=TRUE, level=c("Rainier Mesa", "Spring Mountains", "Pahute Mesa", "Amargosa Valley", "Frenchman & Yucca Flat", "Death Valley", "Ash Meadows", "Oasis Valley"))
plot_metadata$Well_spring <- gsub('Well', 'Groundwater', plot_metadata$Well_spring)
plot_metadata$Well_spring <- gsub('Spring', 'Spring', plot_metadata$Well_spring)
plot_metadata$Well_spring <- factor(plot_metadata$Well_spring, ordered=TRUE, level=c("Groundwater", "Spring", "Tunnel"))
plot_metadata$Piper_group3 <- gsub('Cluster 1', 'Group 1', plot_metadata$Piper_group3)
plot_metadata$Piper_group3 <- gsub('Cluster 2', 'Group 2', plot_metadata$Piper_group3)
plot_metadata$Piper_group3 <- gsub('Cluster 3', 'Group 3', plot_metadata$Piper_group3)
plot_metadata$Piper_group3 <- factor(plot_metadata$Piper_group3, ordered=TRUE, level=c("Group 1", "Group 2", "Group 3"))
plot_metadata$Depth_m <- as.numeric(as.character(plot_metadata$Depth_m))
plot_metadata$Temp_C <- as.numeric(as.character(plot_metadata$Temp_C))
plot_metadata$Depth_sampling_m <- as.numeric(as.character(plot_metadata$Depth_sampling_m))
plot_metadata$pH <- as.numeric(as.character(plot_metadata$pH))

## Check significance between the means (changing variables and also group.by = "Sampling_method")
compare_means(Shannon ~ Loc_sec, data = plot_metadata, group.by = "Sequence_batch", method="wilcox.test", p.adjust.method = "fdr", exact=TRUE)

## Keep variables for plotting
keep_data <- names(plot_metadata) %in% c("Sample_abbrev", "Loc_sec", "Piper_group_ref", "PD", "Shannon", "Simpson", "rock_type", "Well_spring", "tritium_BqL_fct", "Sequence_batch", "Sampling_method")
combo_plot <- plot_metadata[keep_data]
colnames(combo_plot) <- c("Sample_abbrev", "Sampling Method", "Location", "Location Type", "Rock Type", "Overall Chemistry", "Sequence Batch", "Shannon", "Simpson", "PD", "Tritium")

## Melt the dataframe
combo_plot <- reshape2::melt(combo_plot, id=c("Sample_abbrev", "PD", "Shannon", "Simpson"))
combo_plot$value <- factor(combo_plot$value, ordered=TRUE, levels=c("Carbonate", "Volcanic", "Ca-Mg-HCO3", "Na-HCO3", "NaCl", "Groundwater", "Spring", "Tunnel", "Low", "High", "Batch 1", "Batch 2", "Spring_pump", "Pump", "Bailer_Jackpump", "Oasis Valley", "Pahute Mesa", "Rainier Mesa", "Frenchman & Yucca Flat", "Ash Meadows", "Spring Mountains", "Amargosa Valley", "Death Valley"))

## Set up theme
plot_theme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 1, linetype = "solid"),
    panel.border = element_rect(colour="black", size=1, fill=NA),
    strip.background=element_rect(fill='white', colour='white', size = 0),
    strip.text = element_text(face="bold", size=15),
    panel.spacing.x=unit(0.5, "lines"),
    panel.grid.major = element_line(size = 0),
    panel.grid.minor = element_line(size = 0),
    axis.text = element_text(size=15, colour="black"),
    axis.title = element_text(face="bold", size=15),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
    legend.position="right",
    legend.key = element_rect(fill = "white"),
    legend.title = element_text(face="bold", size=15),
    legend.text = element_text(size=15))

plot_guide <- guides(fill = guide_legend(order=2, override.aes = list(shape = 21, alpha=1)),
    shape = guide_legend(order=1, override.aes = list(size = 5, color="black", alpha=1)),
    color = guide_legend(order=3, override.aes = list(size=5, shape = 21, alpha=1)))

## Plot boxplots (Figure S6)
plot_nomargins_PD <- scale_y_continuous(expand = expansion(mult = c(0, 0)), limits=c(0,250), breaks=c(0, 100, 200))
alpha_plot_PD <- ggplot(combo_plot, aes(x=value, y=PD)) +
    stat_boxplot(geom = "errorbar", width = 0.2) +
    geom_boxplot(outlier.size=0, outlier.colour="white", fill="gray88") + 
    geom_jitter(size=2, position=position_jitter(0.2), alpha=0.8, shape=21, fill="grey88") +
    ylab("Faith's PD") +
    facet_grid(. ~ variable, scales="free") +
    plot_theme + plot_nomargins_PD + plot_guide 

plot_nomargins_Shannon <- scale_y_continuous(expand = expansion(mult = c(0, 0)), limits=c(0,7), breaks=c(0, 2, 4, 6))
alpha_plot_Shannon <- ggplot(combo_plot, aes(x=value, y=Shannon)) +
    stat_boxplot(geom = "errorbar", width = 0.2) +
    geom_boxplot(outlier.size=0, outlier.colour="white", fill="gray88") + 
    geom_jitter(size=2, position=position_jitter(0.2), alpha=0.8, shape=21, fill="grey88") +
    ylab("Shannon") +
    facet_grid(. ~ variable, scales="free") +
    plot_theme + plot_nomargins_Shannon + plot_guide 

plot_nomargins_Simpson <- scale_y_continuous(expand = expansion(mult = c(0, 0)), limits=c(0.5,1.2), breaks=c(0.5, 0.75, 1))
alpha_plot_Simpson <- ggplot(combo_plot, aes(x=value, y=Simpson)) +
    stat_boxplot(geom = "errorbar", width = 0.2) +
    geom_boxplot(outlier.size=0, outlier.colour="white", fill="gray88") + 
    geom_jitter(size=2, position=position_jitter(0.2), alpha=0.8, shape=21, fill="grey88") +
    ylab("Simpson") +
    facet_grid(. ~ variable, scales="free") +
    plot_theme + plot_nomargins_Simpson + plot_guide

both <- plot_grid(alpha_plot_PD + theme(axis.text.x = element_blank(), legend.position="none"), 
                  alpha_plot_Shannon + theme(axis.text.x = element_blank(), legend.position="none", strip.text = element_blank()), 
                  alpha_plot_Simpson + theme(legend.position="none", strip.text = element_blank()), 
                  ncol=1, align = "v", axis="b", rel_heights = c(0.5,0.5,1))
both

save_file <- paste("alpha.div.final.shannon.simpson.pd.variables.pdf", sep="")
ggsave(save_file, path = alpha, plot = both, scale = 1, width = 15, height = 7, units = c("in"), dpi = 300)

## Calculate SES.PD, SES.MPD, SES.MNTD
### Check tree
phy <- phy_tree(obj_ps)
phy$tip.label[1:5]
Ntip(phy)
is.rooted(phy)

### convert phylogenety to a distance matrix
phy.dist <- cophenetic(phy)
comm <- as.data.frame(t(otu_table(obj_ps)))
comm <- comm[, phy$tip.label] # make sure counts is in the same order as phy.dist

### Total abundance
comm <- decostand(comm, method = "total") # total abundance
apply(comm, 1, sum) # check total abundance in each sample

### Calculate SES
obj_ps.ses.pd <- ses.pd(comm, phy, null.model="taxa.labels", runs=999, include.root = FALSE) #Faith's PD bias: https://www.biorxiv.org/content/10.1101/579300v1.full
head(obj_ps.ses.pd)
write.csv(obj_ps.ses.pd, "SES.PD.csv")

comm.sesmpd.taxa <- ses.mpd(comm, phy.dist, null.model = "taxa.labels", runs = 999)
head(comm.sesmpd.taxa)
write.csv(comm.sesmpd.taxa, "SES.MPD.csv")

comm.sesmntd.taxa <- ses.mntd(comm, phy.dist, null.model = "taxa.labels", runs = 999)
head(comm.sesmntd.taxa)
write.csv(comm.sesmntd.taxa, "SES.MNTD.csv")

### add to dataframe
plot_metadata <- as.data.frame(plot_metadata)
rownames(plot_metadata) <- plot_metadata$Sample_abbrev
plot_metadata$pd.obs.z <- obj_ps.ses.pd$pd.obs.z[match(rownames(obj_ps.ses.pd), rownames(plot_metadata))]
plot_metadata$mpd.obs.z.taxa <- comm.sesmpd.taxa$mpd.obs.z[match(rownames(comm.sesmpd.taxa), rownames(plot_metadata))]
plot_metadata$mntd.obs.z.taxa <- comm.sesmntd.taxa$mntd.obs.z[match(rownames(comm.sesmntd.taxa), rownames(plot_metadata))]
plot_metadata$pd.obs.p <- obj_ps.ses.pd$pd.obs.p[match(rownames(obj_ps.ses.pd), rownames(plot_metadata))]
plot_metadata$mpd.obs.p.taxa <- comm.sesmpd.taxa$mpd.obs.p[match(rownames(comm.sesmpd.taxa), rownames(plot_metadata))]
plot_metadata$mntd.obs.p.taxa <- comm.sesmntd.taxa$mntd.obs.p[match(rownames(comm.sesmntd.taxa), rownames(plot_metadata))]

### Check significance between the means (changing variables)
compare_means(mntd.obs.z.taxa ~ tritium_BqL_fct, data = plot_metadata, group.by = "Sequence_batch", method="wilcox.test", p.adjust.method = "fdr", exact=TRUE)

### Keep variables for plotting
keep_data <- names(plot_metadata) %in% c("Sample_abbrev", "Loc_sec", "Piper_group_ref", "pd.obs.z", "mpd.obs.z.taxa", "mntd.obs.z.taxa", "pd.obs.p", "mpd.obs.p.taxa", "mntd.obs.p.taxa", "rock_type", "Well_spring", "tritium_BqL_fct", "Sequence_batch", "Sampling_method")
combo_plot <- plot_metadata[keep_data]
colnames(combo_plot) <- c("Sample_abbrev", "Sampling Method", "Location", "Location Type", "Rock Type", "Overall Chemistry", "Sequence Batch", "Tritium", "SES.PD", "SES.MPD", "SES.MNTD", "SES.PD.P", "SES.MPD.P", "SES.MNTD.P")

### Melt the dataframe
combo_plot <- reshape2::melt(combo_plot, id=c("Sample_abbrev","SES.PD", "SES.MPD", "SES.MNTD", "SES.PD.P", "SES.MPD.P", "SES.MNTD.P"))
combo_plot$value <- factor(combo_plot$value, ordered=TRUE, levels=c("Carbonate", "Volcanic", "Ca-Mg-HCO3", "Na-HCO3", "NaCl", "Groundwater", "Spring", "Tunnel", "Low", "High", "Batch 1", "Batch 2", "Spring_pump", "Pump", "Bailer_Jackpump", "Oasis Valley", "Pahute Mesa", "Rainier Mesa", "Frenchman & Yucca Flat", "Ash Meadows", "Spring Mountains", "Amargosa Valley", "Death Valley"))

### Plot boxplots (Figure 3)
plot_nomargins_SESPD <- scale_y_continuous(expand = expansion(mult = c(0, 0)), limits=c(-15,5), breaks=c(-15,-10,-5,0))
alpha_plot_SESPD <- ggplot(combo_plot, aes(x=value, y=SES.PD)) +
    stat_boxplot(geom = "errorbar", width = 0.2) +
    geom_boxplot(outlier.size=0, outlier.colour="white", fill="gray88") + 
    geom_jitter(size=2, position=position_jitter(0.2), alpha=0.8, shape=21, aes(fill=SES.PD.P <0.05)) +
    ylab("SES.PD") +
    facet_grid(. ~ variable, scales="free") +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + 
    scale_fill_manual(name = 'Significance\n(p-value < 0.05)', values = setNames(c('black','#AD1F32'), c(T,F))) +
    plot_theme + plot_guide + plot_nomargins_SESPD

plot_nomargins_SESMPD <- scale_y_continuous(expand = expansion(mult = c(0, 0)), limits=c(-3.5,4), breaks=c(-2,0,2,4))
alpha_plot_SESMPD <- ggplot(combo_plot, aes(x=value, y=SES.MPD)) +
    stat_boxplot(geom = "errorbar", width = 0.2) +
    geom_boxplot(outlier.size=0, outlier.colour="white", fill="gray88") + 
    geom_jitter(size=2, position=position_jitter(0.2), alpha=0.8, shape=21, aes(fill=SES.MPD.P <0.05)) +
    ylab("SES.MPD") +
    facet_grid(. ~ variable, scales="free") +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + 
    scale_fill_manual(name = 'Significance\n(p-value < 0.05)', values = setNames(c('black','#AD1F32'), c(T,F))) +
    plot_theme + plot_guide + plot_nomargins_SESMPD
    
plot_nomargins_SESMNTD <- scale_y_continuous(expand = expansion(mult = c(0, 0)), limits=c(-7,2), breaks=c(-6,-4,-2,0))
alpha_plot_SESMNTD <- ggplot(combo_plot, aes(x=value, y=SES.MNTD)) +
    stat_boxplot(geom = "errorbar", width = 0.2) +
    geom_boxplot(outlier.size=0, outlier.colour="white", fill="gray88") + 
    geom_jitter(size=2, position=position_jitter(0.2), alpha=0.8, shape=21, aes(fill=SES.MNTD.P <0.05)) +
    ylab("SES.MNTD") +
    facet_grid(. ~ variable, scales="free") +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + 
    scale_fill_manual(name = 'Significance\n(p-value < 0.05)', values = setNames(c('black','#AD1F32'), c(T,F))) +
    plot_theme + plot_guide + plot_nomargins_SESMNTD

both <- plot_grid(alpha_plot_SESPD + theme(axis.text.x = element_blank(), legend.position="none"), 
                  alpha_plot_SESMPD + theme(axis.text.x = element_blank(), legend.position="none", strip.text = element_blank()), 
                  alpha_plot_SESMNTD + theme(legend.position="none", strip.text = element_blank()), 
                  ncol=1, align = "v", axis="b", rel_heights = c(0.5,0.5,1))
both

save_file <- paste("alpha.div.final.sespd.mpd.mntd.variables.pdf", sep="")
ggsave(save_file, path = alpha, plot = both, scale = 1, width = 15, height = 7, units = c("in"), dpi = 300)
