# 06. Beta Diversity
## Figure 4, Figure S7, Table 1

## References
### Martino, C.; Morton, J. T.; Marotz, C. A.; Thompson, L. R.; Tripathi, A.; Knight, R.; Zengler, K. A Novel Sparse Compositional Technique Reveals Microbial Perturbations. mSystems 2019, 4 (1). https://doi.org/10.1128/mSystems.00016-19.
### Gloor, G. B.; Macklaim, J. M.; Pawlowsky-Glahn, V.; Egozcue, J. J. Microbiome Datasets Are Compositional: And This Is Not Optional. Front. Microbiol. 2017, 8. https://doi.org/10.3389/fmicb.2017.02224.
### Silverman, J. D.; Washburne, A. D.; Mukherjee, S.; David, L. A. A Phylogenetic Transform Enhances Analysis of Compositional Microbiota Data. eLife 2017, 6, e21887. https://doi.org/10.7554/eLife.21887.
### https://www.davidzeleny.net/anadat-r/doku.php/en:rda_cca_examples

## Set Path
beta <- file.path(paste(path_phy, "/Beta_diversity", sep=""))
dir.create(beta, showWarnings=FALSE)
setwd(beta)

## Phyloseq object
obj_ps

## Phyloseq object of subset samples (remove Batch 2 sequenced and tunnel-collected samples)
obj_ps_sub <- subset_samples(obj_ps_log, Sample_abbrev != "RM1" & Sample_abbrev != "RM2" & Sample_abbrev != "DV1" & Sample_abbrev != "OV1" & Sample_abbrev != "OV2" & Sample_abbrev != "AV1" & Sample_abbrev != "AV2" & Sample_abbrev != "AV3" & Sample_abbrev != "YF2") #For NMDS, also removed DV3 due to stress issues in WUniFrac
obj_ps_sub <- prune_taxa(taxa_sums(obj_ps_sub) > 0, obj_ps_sub) # remove OTUs with zero count
obj_ps_sub

## Set up theme
plot_theme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 1, linetype = "solid"),
    panel.border = element_rect(colour="black", size=1, fill=NA),
    strip.background=element_rect(fill='white', colour='white', size = 0),
    strip.text = element_text(face="bold", size=20),
    panel.spacing.x=unit(0.5, "lines"),
    panel.grid.major = element_line(size = 0),
    panel.grid.minor = element_line(size = 0),
    axis.text = element_text(size=15, colour="black"),
    axis.title = element_text(face="bold", size=20),
    legend.position="right",
    legend.key = element_rect(fill = "transparent"),
    legend.title = element_text(face="bold", size=20),
    legend.text = element_text(size=20),
    legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent", colour = NA))

plot_guide <- guides(fill = guide_legend(order=1, override.aes = list(shape = 21, alpha=1, size = 5)),
    shape = guide_legend(order=2, override.aes = list(size = 5, color="black", alpha=1)),
    color = FALSE)

loc_sec_color = c("Amargosa Valley" = "#043378", "Ash Meadows" = "#E50006", "Death Valley" = "#38AA31", "Frenchman and Yucca Flat" = "#1187A5", "Pahute Mesa" = "#7E478D", "Rainier Mesa" = "#FB9E7F", "Spring Mountains" = "#9B0020", "Oasis Valley" = "#9DA7A7")

## Preliminary ordination plots
### unpound the distance want to use
#dist <- "wunifrac"
#dist <- "bray"
#dist <- "unifrac"
dist <- "jaccard"

ord_meths<- c("DCA", "CCA", "RDA", "MDS", "PCoA", "NMDS")
color <- "Loc_sec" ## make sure to change the geom_point aes below ##

plist = llply(as.list(ord_meths), function(i, physeq, dist){
    ordi = ordinate(physeq, method=i, distance=dist)
    plot_ordination(physeq, ordi, "samples", color=color)
}, obj_ps, dist)

names(plist) <- ord_meths

pdataframe = ldply(plist, function(x){
    df = x$data[, 1:2]
    colnames(df) = c("Axis_1", "Axis_2")
    return(cbind(df, x$data))})

names(pdataframe)[1] = "method"

### Plot multiple ordination methods
plot.ord <- ggplot(pdataframe, aes(Axis_1, Axis_2), color = color, title = "Ordination") +
    geom_point(aes(color = Loc_sec), alpha = 0.7, size = 4) + ## change color ##
    facet_wrap(~method, scales="free") + plot_theme

save_file_plot <- paste(beta, "/multiple.ordination.obj_ps.", color, ".", dist, ".pdf", sep="")
ggsave(save_file_plot, scale = 1, width = 7, height = 5, units = c("in"), dpi = 300)

##==================== NMDS Weighted UniFrac (repeat this for subset of samples obj_ps_sub)
ordi_wunifrac_NMDS = ordinate(obj_ps, method="NMDS", distance="wunifrac")
ordi_wunifrac_NMDS
stressplot(ordi_wunifrac_NMDS)
ordi <- ordi_wunifrac
ordi.scores <- as.data.frame(ordi$points)
ordi.scores$Sample_abbrev <- rownames(ordi.scores)
alpha.div.metadata2$Loc_sec <- gsub('_',' ', alpha.div.metadata2$Loc_sec)
ordi_data <- merge(alpha.div.metadata2, ordi.scores, by = c('Sample_abbrev'))
ordi_data$Loc_sec <- factor(ordi_data$Loc_sec, ordered = TRUE, levels = c("Amargosa Valley", "Ash Meadows", "Death Valley", "Frenchman and Yucca Flat", "Pahute Mesa", "Rainier Mesa", "Spring Mountains", "Oasis Valley"))

### Calculate the hulls for each group
hull <- ordi_data %>%
  group_by(Loc_sec) %>%
  slice(chull(MDS1, MDS2))

### Plot NMDS wUniFrac (Figure 4B)
nmds_whull <- ggplot(ordi_data, aes(x=MDS1, y = MDS2)) +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + geom_vline(xintercept=0, size=.2, linetype = "dashed", color="black") + 
    geom_polygon(data = hull, alpha = 0.1, aes(color=Loc_sec, fill=Loc_sec)) +
    geom_point(size=3, color="black", shape=21, aes(fill=Loc_sec)) +
    xlab(paste("NMDS1")) +
    ylab(paste("NMDS2")) +
    scale_fill_lancet(name = "Location") + 
    scale_color_lancet(name = "") + 
    plot_theme + plot_guide

nmds_whull_label <- nmds_whull + geom_label_repel(aes(label = Sample_abbrev), box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50')

### Plot with rock type
hull <- ordi_data %>%
  group_by(rock_type) %>%
  slice(chull(MDS1, MDS2))

nmds_wuni_rock <- ggplot(ordi_data, aes(x=MDS1, y = MDS2)) +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + geom_vline(xintercept=0, size=.2, linetype = "dashed", color="black") +
    geom_polygon(data = hull, alpha = 0.1, aes(color=rock_type, fill=rock_type)) +
    geom_point(size=3, color="black", shape=21, aes(fill=rock_type)) +
    xlab(paste("NMDS1")) +
    ylab(paste("NMDS2")) +
    scale_fill_lancet(name = "Rock Type") +
    scale_color_lancet(name = "") +
    plot_theme + plot_guide

### Plot with Piper group
hull <- ordi_data %>%
  group_by(Piper_group3) %>%
  slice(chull(MDS1, MDS2))

nmds_wuni_piper <- ggplot(ordi_data, aes(x=MDS1, y = MDS2)) +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + geom_vline(xintercept=0, size=.2, linetype = "dashed", color="black") +
    geom_polygon(data = hull, alpha = 0.1, aes(color=Piper_group_ref, fill=Piper_group_ref)) +
    geom_point(size=3, color="black", shape=21, aes(fill=Piper_group_ref)) +
    xlab(paste("NMDS1")) +
    ylab(paste("NMDS2")) +
    scale_fill_lancet(name = "Overall Chemistry") +
    scale_color_lancet(name = "") +
    plot_theme + plot_guide

### Plot with Sequence Batch
hull <- ordi_data %>%
  group_by(Sequence_batch) %>%
  slice(chull(MDS1, MDS2))

nmds_wuni_seq <- ggplot(ordi_data, aes(x=MDS1, y = MDS2)) +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + geom_vline(xintercept=0, size=.2, linetype = "dashed", color="black") + 
    geom_polygon(data = hull, alpha = 0.1, aes(color=Sequence_batch, fill=Sequence_batch)) +
    geom_point(size=3, color="black", shape=21, aes(fill=Sequence_batch)) +
    xlab(paste("NMDS1")) +
    ylab(paste("NMDS2")) +
    scale_fill_lancet(name = "Sequence batch") + 
    scale_color_lancet(name = "") + 
    plot_theme + plot_guide

### Plot with Sampling method
hull <- ordi_data %>%
  group_by(Sampling_method) %>%
  slice(chull(MDS1, MDS2))

nmds_wuni_method <- ggplot(ordi_data, aes(x=MDS1, y = MDS2)) +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + geom_vline(xintercept=0, size=.2, linetype = "dashed", color="black") +
    geom_polygon(data = hull, alpha = 0.1, aes(color=Sampling_method, fill=Sampling_method)) +
    geom_point(size=3, color="black", shape=21, aes(fill=Sampling_method)) +
    xlab(paste("NMDS1")) +
    ylab(paste("NMDS2")) +
    scale_fill_lancet(name = "Sampling Method") +
    scale_color_lancet(name = "") +
    plot_theme + plot_guide

### Plot with Location type
hull <- ordi_data %>%
  group_by(Well_spring) %>%
  slice(chull(MDS1, MDS2))

nmds_wuni_loctype <- ggplot(ordi_data, aes(x=MDS1, y = MDS2)) +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + geom_vline(xintercept=0, size=.2, linetype = "dashed", color="black") +
    geom_polygon(data = hull, alpha = 0.1, aes(color=Well_spring, fill=Well_spring)) +
    geom_point(size=3, color="black", shape=21, aes(fill=Well_spring)) +
    xlab(paste("NMDS1")) +
    ylab(paste("NMDS2")) +
    scale_fill_lancet(name = "Location Type") +
    scale_color_lancet(name = "") +
    plot_theme + plot_guide

### Plot with TOC
nmds_wuni_toc <- ggplot(ordi_data, aes(x=MDS1, y = MDS2)) +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + geom_vline(xintercept=0, size=.2, linetype = "dashed", color="black") +
    geom_point(shape=21, size = 5, color="black", aes(fill=TOC_mgCL)) +
    xlab(paste("NMDS1")) +
    ylab(paste("NMDS2")) +
    scale_fill_gradient2(name = "TOC (mg-C/L)", breaks = c(0, 20, 40), limits = c(0,40), low="#f72585", mid="#cdb4db", high="#023e8a", midpoint=20) +
    scale_color_lancet(name = "") +
    plot_theme + plot_guide

### Plot with temperature
plot_guide <- guides(fill = guide_colourbar(frame.linewidth = 1, frame.colour = "black", ticks = TRUE, ticks.colour = "black", ticks.linewidth = 1))

nmds_wuni_temp <- ggplot(ordi_data, aes(x=MDS1, y = MDS2)) +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + geom_vline(xintercept=0, size=.2, linetype = "dashed", color="black") + 
    geom_point(shape=21, size = 5, color="black", aes(fill=Temp_C)) +
    xlab(paste("NMDS1")) +
    ylab(paste("NMDS2")) +
    scale_fill_gradient2(name = "Temperature (ºC)", breaks = c(5, 20, 40, 60), limits = c(5,60), low="#f72585", mid="#cdb4db", high="#023e8a", midpoint=30) + 
    scale_color_lancet(name = "") + 
    plot_theme + plot_guide

### Plot with Depth
nmds_wuni_depth <- ggplot(ordi_data, aes(x=MDS1, y = MDS2)) +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + geom_vline(xintercept=0, size=.2, linetype = "dashed", color="black") + 
    geom_point(shape=21, size = 5, color="black", aes(fill=relative_depth)) +
    xlab(paste("NMDS1")) +
    ylab(paste("NMDS2")) +
    scale_fill_gradient2(name = "Depth (m)", breaks = c(0, 400, 800, 1200), limits = c(0,1200), low="#f72585", mid="#cdb4db", high="#023e8a", midpoint=600) + 
    scale_color_lancet(name = "") + 
    plot_theme + plot_guide

both <- plot_grid(nmds_wuni_rock + theme(legend.position="none"),
                    nmds_wuni_seq + theme(legend.position="none"),
                    nmds_wuni_loctype + theme(legend.position="none"),
                    nmds_wuni_temp + theme(legend.position="none"),
                    nmds_wuni_piper + theme(legend.position="none"),
                    nmds_wuni_method + theme(legend.position="none"),
                    nmds_wuni_toc + theme(legend.position="none"),
                    nmds_wuni_depth + theme(legend.position="none"),
                    NULL,
                    ncol=4, align = "v", axis="b")

save_file <- paste("Combo_for_supplementary_wunifrac.svg", sep="")
ggsave(save_file, path = beta, plot = both, scale = 1, width = 15, height = 10, units = c("in"), dpi = 300)

##==================== NMDS UniFrac (repeat this for subset of samples obj_ps_sub)
ordi_unifrac_NMDS = ordinate(obj_ps, method="NMDS", distance="unifrac")
ordi_unifrac_NMDS
stressplot(ordi_unifrac_NMDS)

ordi <- ordi_unifrac_NMDS
ordi.scores <- as.data.frame(ordi$points)
ordi.scores$Sample_abbrev <- rownames(ordi.scores)
head(ordi.scores)
alpha.div.metadata2$Loc_sec <- gsub('_',' ', alpha.div.metadata2$Loc_sec)
ordi_data <- merge(alpha.div.metadata2, ordi.scores, by = c('Sample_abbrev'))
ordi_data$Loc_sec <- factor(ordi_data$Loc_sec, ordered = TRUE, levels = c("Amargosa Valley", "Ash Meadows", "Death Valley", "Frenchman and Yucca Flat", "Pahute Mesa", "Rainier Mesa", "Spring Mountains", "Oasis Valley"))

### Plot NMDS UniFrac (Figure 4A)
### Calculate the hulls for each group
hull <- ordi_data %>%
  group_by(Loc_sec) %>%
  slice(chull(MDS1, MDS2))

nmds_uni_hull <- ggplot(ordi_data, aes(x=MDS1, y = MDS2)) +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + geom_vline(xintercept=0, size=.2, linetype = "dashed", color="black") + 
    geom_polygon(data = hull, alpha = 0.1, aes(color=Loc_sec, fill=Loc_sec)) +
    geom_point(size=3, color="black", shape=21, aes(fill=Loc_sec)) +
    xlab(paste("NMDS1")) +
    ylab(paste("NMDS2")) +
    scale_fill_lancet(name = "Location") + 
    scale_color_lancet(name = "") + 
    plot_theme + plot_guide

nmds_uni_hull_label <- nmds_uni_hull + geom_label_repel(aes(label = Sample_abbrev), box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50')

##==================== Plot both NMDS (Figure 4)
both <- plot_grid(nmds_uni_hull + theme(legend.position="none"),
                 nmds_whull + theme(legend.position="none"), 
                 ncol=2, align = "v", axis="b")
save_file <- paste("wuni.unifrac.nmds.svg", sep="")
ggsave(save_file, path = beta, plot = both, scale = 1, width = 10, height = 5, units = c("in"), dpi = 300)

both_label <- plot_grid(nmds_uni_whull_label + theme(legend.position="none"),
                        nmds_whull_label + theme(legend.position="none"),
                        nmds_uni_whull_sub + theme(legend.position="none"),
                        nmds_whull_sub + theme(legend.position="none"),
                        ncol=2, align = "v", axis="b")
save_file <- paste("wuni.unifrac.nmds.label.svg", sep="")
ggsave(save_file, path = beta, plot = both_label, scale = 1, width = 10, height = 5, units = c("in"), dpi = 300)

##==================== Statistics (repeat this for subset of samples obj_ps_sub)
### https://microbiome.github.io/tutorials/PERMANOVA.html
### https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/multivariate-comparisons-of-microbial-community-composition.html

### Obtain the distance matrix
uwUF.dist = UniFrac(obj_ps, weighted=FALSE, normalized=TRUE)
head(uwUF.dist)

wUF.dist = UniFrac(obj_ps, weighted=TRUE, normalized=TRUE)
head(wUF.dist)

### Organize data
metadata <- as.data.frame(as.matrix(sample_data(obj_ps)))
cols.fct <- c("Loc_sec","Well_spring","rock_type","Piper_group3","Sequence_batch","Loc_DV3", "Sampling_method)
metadata[cols.fct] <- lapply(metadata[cols.fct],as.factor)
metadata[, c(13:32)] <- sapply(metadata[, c(13:32)], as.numeric)
metadata$tritium_BqL_fct <- cut(metadata$tritium_BqL, breaks = c(0, 100, Inf), labels = c("Low", "High"), ordered_result = TRUE)
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$ID_keep
metadata_stats <- metadata

### Adonis
results_uw <- lapply(colnames(metadata_stats), function(x){
  form <- as.formula(paste("uwUF.dist ~ metadata_stats$", x, sep="")) 
  z <- adonis(form, permutations=999)
  return(as.data.frame(z$aov.tab)) #convert anova table to a data frame})

results_w <- lapply(colnames(metadata_stats), function(x){
  form <- as.formula(paste("wUF.dist ~ metadata_stats$", x, sep="")) 
  z <- adonis(form, permutations=999)
  return(as.data.frame(z$aov.tab)) #convert anova table to a data frame})

results_uw[[1]]
names(results_uw) <- colnames(metadata_stats)
results_uw <- do.call(rbind, results_uw)
write.csv(results_uw, "NMDS_unifrac_beta_diversity_adonis_results.csv")

results_w[[1]]
names(results_w) <- colnames(metadata_stats)
results_w <- do.call(rbind, results_w)
write.csv(results_w, "NMDS_wunifrac_beta_diversity_adonis_results.csv")

### Test specific questions (Table 1)
#### location-type variables
adonis2(uwUF.dist~metadata_stats$Loc_sec+metadata_stats$relative_depth+metadata_stats$Well_spring+metadata_stats$Sequence_batch+metadata_stats$Sampling_method, permutations=999, by="margin")
adonis2(wUF.dist~metadata_stats$Loc_sec+metadata_stats$relative_depth+metadata_stats$Well_spring+metadata_stats$Sequence_batch+metadata_stats$Sampling_method, permutations=999, by="margin")

#### geochemical-type variables
adonis2(uwUF.dist~metadata_stats$Temp_C+metadata_stats$rock_type+metadata_stats$pH+metadata_stats$tritium_BqL_fct+metadata_stats$TOC_mgCL+metadata_stats$Sequence_batch+metadata_stats$Sampling_method, permutations=999, by="margin")
adonis2(wUF.dist~metadata_stats$Temp_C+metadata_stats$rock_type+metadata_stats$pH+metadata_stats$tritium_BqL_fct+metadata_stats$TOC_mgCL+metadata_stats$Sequence_batch+metadata_stats$Sampling_method, permutations=999, by="margin")

### Pairwise comparisons
pairwise <- pairwise.adonis2(wUF.dist~Loc_sec, data=metadata_stats, perm=999, method="euclidean")
write.csv(pairwise, "DEICODE_RPCA_beta_diversity_pairwiseadonis_location_results999.csv")

### Check homogeneity (change the variable)
anova(betadisper(wUF.dist, metadata_stats$tritium_BqL_fct))
anova(betadisper(uwUF.dist, metadata_stats$Loc_sec))
permutest(betadisper(wUF.dist, metadata_stats$Loc_sec), pairwise = TRUE)
permutest(betadisper(uwUF.dist, metadata_stats$Loc_sec), pairwise = TRUE)

### ANOSIM for categorical variables (change the variable)
anosim(wUF.dist, metadata_stats$tritium_BqL_fct, permutations = 999)

##==================== CCA
### Prepare metadata
CC_metadata <- alpha.div.metadata2

# Remove Oasis Valley and variables correlated with other variables (spearman's rho > 0.7); only keep geochemical-type variables
CC_metadata <- CC_metadata[!grepl("Oasis_Valley", CC_metadata$Loc_sec),]
CC_metadata$Sr_mgL <- NULL
CC_metadata$Depth_m <- NULL
CC_metadata$Depth_from_water_table_m <- NULL
CC_metadata.num <- select_if(CC_metadata, is.numeric) # only keep numeric variables
CC_metadata.num.std <- decostand(CC_metadata.num, method='standardize') # standardize the variables
str(CC_metadata.num.std)
summary(CC_metadata.num.std)
X <- CC_metadata.num.std
colnames(X)

### Prepare otu table
# remove OV samples
obj_ps_sub <- subset_samples(obj_ps, Sample_abbrev != "OV1" & Sample_abbrev != "OV2")

# remove zero OTU count
obj_ps_sub <- prune_taxa(taxa_sums(obj_ps_sub) > 0, obj_ps_sub)
obj_ps_sub

# Obtain the OTU table
library(QsRutils)
Y.counts <- veganotu(obj_ps_sub)
rownames(Y.counts)

### Test whether species composition is heterogeneous or homogeneous
DCA <- decorana(log1p(Y.counts))
DCA
# Note: axis lengths are > 4: heterogeneous, unimodal method (CCA)

### Conduct CCA
spe.log <- log1p(Y.counts)
spe.hell <- decostand(spe.log, 'hell')

# CCA on all the data
cca <- cca(spe.hell ~ ., data = X)
cca

# CCA on just microbial data to confirm the axis variance without geochemical data
cca_noData <- cca(spe.hell)
cca_noData

# Get the contribution of each axis to the variance
constrained_eig <- cca$CCA$eig/cca$tot.chi*100
unconstrained_eig <- cca$CA$eig/cca$tot.chi*100
expl_var <- c(constrained_eig, unconstrained_eig)
barplot (expl_var[1:20], col = c(rep ('red', length (constrained_eig)), rep ('black', length (unconstrained_eig))),
         las = 2, ylab = '% variation')
  
# Quick plot
ordiplot(cca)

# R2
R2.obs <- RsquareAdj(cca)$r.squared
R2.obs

# Is the CCA worth to interpret? Permutation test of statistical significance
test_all <- anova(cca)
test_all.adj <- test_all
test_all.adj$`Pr(>F)` <- p.adjust(test_all$`Pr(>F)`, method = 'holm')
test_all.adj

# Are the axis significant?
test_axis <- anova(cca, by="axis")
test_axis.adj <- test_axis
test_axis.adj$`Pr(>F)` <- p.adjust(test_axis$`Pr(>F)`, method = 'holm')
test_axis.adj

# Are the geochemical variables significant?
test_margin <- anova(cca, by="margin")
test_margin.adj <- test_margin
test_margin.adj$`Pr(>F)` <- p.adjust(test_margin$`Pr(>F)`, method = 'holm')
test_margin.adj

# Perform forward and backward selection of explanatory variables
cca1 <- cca(spe.hell ~ ., data=X) # full model
cca0 <- cca(spe.hell ~ 1, data=X) # intercept-only (null) model

step.env <- ordistep(cca0, scope=formula(cca1), direction='both', permutations=1000)
step.env

# What are the significant variables after selection?
step.env$anova
step.env_adj <- step.env
step.env_adj$anova$`Pr(>F)` <- p.adjust(step.env$anova$`Pr(>F)`, method = 'holm', n = ncol(X))
step.env_adj$anova

# Perform CCA on the selected variables
cca.lim <- cca(formula = spe.hell ~ Temp_C + Ca_mgL + NO3_mgL + Na_mgL + TOC_mgCL, data = X)
cca.lim

# Get the contribution of each axis to the variance
constrained_eig <- cca.lim$CCA$eig/cca.lim$tot.chi*100
unconstrained_eig <- cca.lim$CA$eig/cca.lim$tot.chi*100
expl_var <- c(constrained_eig, unconstrained_eig)
barplot (expl_var[1:20], col = c(rep ('red', length (constrained_eig)), rep ('black', length (unconstrained_eig))),
         las = 2, ylab = '% variation')

constrained_eig

# Quick plot
ordiplot(cca.lim)

# R2
R2.obs <- RsquareAdj(cca.lim)$r.squared
R2.obs

# Is the CCA worth to interpret? Permutation test of statistical significance
test_all <- anova(cca.lim)
test_all.adj <- test_all
test_all.adj$`Pr(>F)` <- p.adjust(test_all$`Pr(>F)`, method = 'holm')
test_all.adj

# Are the axis significant?
test_axis <- anova(cca.lim, by="axis")
test_axis.adj <- test_axis
test_axis.adj$`Pr(>F)` <- p.adjust(test_axis$`Pr(>F)`, method = 'holm')
test_axis.adj

# Are the geochemical variables significant?
test_margin <- anova(cca.lim, by="margin")
test_margin.adj <- test_margin
test_margin.adj$`Pr(>F)` <- p.adjust(test_margin$`Pr(>F)`, method = 'holm')
test_margin.adj

# obtain site scores
CC_metadata$Label <- rownames(CC_metadata)
sit <- merge(CC_metadata, fortify(cca.lim, axes = 1:2, display = c("wa")), by="Label")
rownames(sit) <- sit$Label
sit$Score <- NULL
sit$Label <- NULL
head(sit)

# obtain biplot vectors
vec <- fortify(cca.lim, axes = 1:2, display = c("bp"))
vec$Score <- NULL
vec$Label <- gsub('Temp_C', 'Temp (ºC)', vec$Label)
vec$Label <- gsub('_mgL', ' (mg/L)', vec$Label)
vec$Label <- gsub('_mgCL', ' (mg-C/L)', vec$Label)
vec

### Plot figure (Figure 4C)
# Set up theme
plot_theme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 1, linetype = "solid"),
    panel.border = element_rect(colour="black", size=1, fill=NA),
    strip.background=element_rect(fill='white', colour='white', size = 0),
    strip.text = element_text(face="bold", size=20),
    panel.spacing.x=unit(0.5, "lines"),
    panel.grid.major = element_line(size = 0),
    panel.grid.minor = element_line(size = 0),
    axis.text = element_text(size=15, colour="black"),
    axis.title = element_text(face="bold", size=20),
    legend.position="right",
    legend.key = element_rect(fill = "transparent"),
    legend.title = element_text(face="bold", size=20),
    legend.text = element_text(size=20),
    legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent", colour = NA))

plot_guide <- guides(fill = guide_legend(order=1, override.aes = list(shape = 21, alpha=1, size = 5)),
    shape = guide_legend(order=2, override.aes = list(size = 5, color="black", alpha=1)),
    color = "none")

# Calculate the hulls for each group
hull <- sit %>%
  group_by(rock_type) %>%
  slice(chull(CCA1, CCA2))

# use these to adjust length of arrows and position of arrow labels
adj.vec <- 2
adj.txt <- 2.1

cca_fig <- ggplot(sit, aes(x=CCA1, y = CCA2)) +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + geom_vline(xintercept=0, size=.2, linetype = "dashed", color="black") +

    # points and polygons
    geom_point(size=3, color="black", aes(fill=Loc_sec, shape=rock_type)) +
    geom_polygon(data = hull, alpha = 0.1, color="black", aes(linetype=rock_type)) +

    # vectors
    geom_segment(data=vec, inherit.aes=F, mapping=aes(x=0, y=0, xend=adj.vec*CCA1, yend=adj.vec*CCA2), arrow=arrow(length=unit(0.2, 'cm'))) +
    geom_text(data=vec, size=6, inherit.aes=F, mapping=aes(x=adj.txt*CCA1, y=adj.txt*CCA2, label=Label)) +

    # make plot pretty
    xlab(paste("CCA1 (4.2%)")) +
    ylab(paste("CCA2 (3.8%)")) +
    scale_shape_manual(name = "Rock Type", values=c(22, 24)) +
    scale_fill_lancet(name = "Location") +
    scale_linetype_manual(name="", values=c(2,3)) +
    #geom_label_repel(aes(label = Sample_abbrev), box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
    plot_theme + plot_guide

cca_fig

save_file <- paste("CCA_draft.pdf", sep="")
ggsave(save_file, path = beta, scale = 1, width = 10, height = 5, units = c("in"), dpi = 300)

##==================== DEICODE RPCA (repeat this for subset of samples obj_ps_sub)
deicode <- file.path(paste(beta, "/deicode_RPCA", sep=""))
dir.create(deicode, showWarnings=FALSE)
setwd(deicode)

### Convert to Qiime2 format
tax <- as(tax_table(obj_ps_filt),"matrix")
tax <- as.data.frame(tax)
tax <- subset(tax, select = -c(FULL_ID))
tax$Kingdom <- paste("k", tax$Kingdom, sep="__")
tax$Phylum <- paste("p", tax$Phylum, sep="__")
tax$Class <- paste("c", tax$Class, sep="__")
tax$Order <- paste("o", tax$Order, sep="__")
tax$Family <- paste("f", tax$Family, sep="__")
tax$Genus <- paste("g", tax$Genus, sep="__")
tax$Species <- paste("s", tax$Species, sep="__")
tax_cols <- c("Kingdom", "Phylum", "Class","Order","Family","Genus","Species")
tax$taxonomy <- do.call(paste, c(tax[tax_cols], sep=";"))
for(co in tax_cols) tax[co]<-NULL
write.table(tax, "tax_for_qiime2.txt", quote=FALSE, col.names=FALSE, sep="\t")

otu <- as(otu_table(obj_ps_filt),"matrix")
otu_biom <- make_biom(data=otu)
write_biom(otu_biom,"otu_biom.biom")
write.table(otu_table(obj_ps_filt), file = "otu_table.txt", sep = "\t", row.names = TRUE, col.names = NA)

write.table(sample_data(obj_ps_filt), file = "metadata_for_qiime2.txt", sep = "\t", row.names = TRUE, col.names = NA)

### Import to Qiime2 (bash)
conda activate qiime2-2020.6

wd=<path to working directory>
cd $wd

### Organize data and convert to biom format
sed 's/"//g' metadata_for_qiime2.txt > metadata_for_qiime2_fixed.txt #also add #SampleID to header
biom convert -i otu_biom.biom -o otu_biom_HDF5.biom --to-hdf5
biom add-metadata -i otu_biom_HDF5.biom -o otu_wTax_metadata.biom --observation-metadata-fp tax_for_qiime2.txt --sc-separated taxonomy --observation-header OTUID,taxonomy --sample-metadata-fp metadata_for_qiime2_fixed.txt

qiime tools import \
--input-path otu_biom_HDF5.biom \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV210Format \
--output-path feature-table.qza

### import tax table to qiime
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path tax_for_qiime2.txt \
--output-path taxonomy.qza

### Summarize data
qiime feature-table summarize \
--i-table feature-table.qza \
--m-sample-metadata-file metadata_for_qiime2_fixed.txt \
--o-visualization summary_vis.qzv

### DEICODE RPCA automatic
qiime deicode auto-rpca \
--i-table feature-table.qza \
--p-min-feature-count 0 \
--p-min-sample-count 0 \
--o-biplot ordination_auto_components.qza \
--p-max-iterations 5 \
--p-min-feature-frequency 0 \
--o-distance-matrix distance_auto_components.qza

qiime emperor biplot \
--i-biplot ordination_auto_components.qza \
--m-sample-metadata-file metadata_for_qiime2_fixed.txt \
--m-feature-metadata-file taxonomy.qza \
--o-visualization biplot_auto_components.qzv

qiime tools view biplot_auto_components.qzv

### Import to RPCA to R
rpca <- read_qza("ordination_auto_components.qza")
metadata <- DATA_PHYLOSEQ_FIXED
rpca_vectors <- as.data.frame(as.matrix(rpca$data$Vectors))
colnames(rpca_vectors)[grep("SampleID", colnames(rpca_vectors))] <- "ID_keep"
metadata$ID_keep <- metadata$Sample_abbrev
rpca_data <- merge(metadata, rpca_vectors, by = c('ID_keep'))
rpca_data %>% convert(num(PC1, PC2)) -> rpca_data_num
rpca$data$ProportionExplained[1]
rpca$data$ProportionExplained[2]
rpca$data$ProportionExplained[3]

xaxis_text <- paste("PC1: 72.5%")
yaxis_text <- paste("PC2: 27.5%")

### Plot ordination (Figure S7)
#### Calculate the hulls for each group
hull <- rpca_data_num %>%
  group_by(Loc_sec) %>%
  slice(chull(PC1, PC2))
  
rpca_hull <- ggplot(rpca_data_num, aes(x=PC1, y = PC2)) +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + geom_vline(xintercept=0, size=.2, linetype = "dashed", color="black") + 
    geom_polygon(data = hull, alpha = 0.1, aes(color=Loc_sec, fill=Loc_sec)) +
    geom_point(size=3, color="black", shape=21, aes(fill=Loc_sec)) +
    xlab(xaxis_text) +
    ylab(yaxis_text) +
    scale_fill_manual(name = "Location", values = loc_sec_color) +
    scale_color_manual(name = "", values = loc_sec_color) +
    plot_theme + plot_guide

save_file <- paste("RPCA_draft.svg", sep="")
ggsave(save_file, path = deicode, scale = 1, width = 9.5, height = 5, units = c("in"), dpi = 300)

### Plot ordination (Figure S7)
#### Calculate the hulls for each group
hull <- rpca_data_num %>%
  group_by(Sequence_batch) %>%
  slice(chull(PC1, PC2))

rpca_hull_Seq <- ggplot(rpca_data_num, aes(x=PC1, y = PC2)) +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + geom_vline(xintercept=0, size=.2, linetype = "dashed", color="black") +
    geom_polygon(data = hull, alpha = 0.1, aes(color=Sequence_batch, fill=Sequence_batch)) +
    geom_point(size=3, color="black", shape=21, aes(fill=Sequence_batch)) +
    xlab(xaxis_text) +
    ylab(yaxis_text) +
    scale_fill_manual(name = "Sequence batch", values = c("black", "red")) +
    scale_color_manual(name = "", values = c("black", "red")) +
    plot_theme + plot_guide
rpca_hull_Seq

save_file <- paste("rpca_draft_sequenceBatch.svg", sep="")
ggsave(save_file, path = beta, scale = 1, width = 9.5, height = 5, units = c("in"), dpi = 300)

##==================== CLR transformation (repeat this for subset of samples obj_ps_sub)
set.seed(200810)

### Obtain the OTU and metadata table
otu_table <- as.data.frame(otu_table(obj_ps))
metadata <- as.data.frame(as.matrix(sample_data(obj_ps)))

### Transform the OTU table and replace 0 values with an estimate
f.n0 <- cmultRepl(t(otu_table), method="CZM", label=0, output="p-counts")

### CLR Transformation
f.clr <- codaSeq.clr(f.n0, IQLR=FALSE)

### Add clr to phyloseq object
ps_clr <- phyloseq(
    otu_table(t(f.clr), taxa_are_rows = TRUE),
    sample_data(sample_data(obj_ps)),
    phy_tree(filt_tree),
    tax_table(tax_table(obj_ps))
)
ps_clr

### Ordinate
ordi = ordinate(ps_clr, "PCoA", "euclidean")

### Get scores and add metadata
ordi.scores <- as.data.frame(ordi$vectors)
ordi.scores$Sample_abbrev <- rownames(ordi.scores)
ordi_data <- merge(alpha.div.metadata2, ordi.scores, by = c('Sample_abbrev'))
ordi_data$Loc_sec <- gsub("_", " ", ordi_data$Loc_sec)
ordi_data$Loc_sec <- factor(ordi_data$Loc_sec, ordered = TRUE, levels = c("Amargosa Valley", "Ash Meadows", "Death Valley", "Frenchman and Yucca Flat", "Pahute Mesa", "Rainier Mesa", "Spring Mountains", "Oasis Valley"))

### Plot ordination (Figure S7)
#### Calculate the hulls for each group
hull <- ordi_data %>%
  group_by(Loc_sec) %>%
  slice(chull(Axis.1, Axis.2))

ordi_clr <- ggplot(ordi_data, aes(x=Axis.1, y = Axis.2)) +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + geom_vline(xintercept=0, size=.2, linetype = "dashed", color="black") +
    geom_polygon(data = hull, alpha = 0.1, aes(color=Loc_sec, fill=Loc_sec)) +
    geom_point(size=5, color="black", shape=21, aes(fill=Loc_sec)) +
    xlab(paste("PC1: 13.6%")) +
    ylab(paste("PC2: 9.6%")) +
    scale_fill_manual(name = "Location", values = loc_sec_color) +
    scale_color_manual(name = "", values = loc_sec_color) +
    plot_theme + plot_guide
ordi_clr

save_file <- paste("PCoA_clr_draft.svg", sep="")
ggsave(save_file, path = beta, scale = 1, width = 9.5, height = 5, units = c("in"), dpi = 300)

### Plot ordination (Figure S7)
# Calculate the hulls for each group
hull <- ordi_data %>%
  group_by(Sequence_batch) %>%
  slice(chull(Axis.1, Axis.2))

# plot whole
ordi_clr_Seq <- ggplot(ordi_data, aes(x=Axis.1, y = Axis.2)) +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + geom_vline(xintercept=0, size=.2, linetype = "dashed", color="black") +
    geom_polygon(data = hull, alpha = 0.1, aes(color=Sequence_batch, fill=Sequence_batch)) +
    geom_point(size=5, color="black", shape=21, aes(fill=Sequence_batch)) +
    xlab(paste("PC1: 13.6%")) +
    ylab(paste("PC2: 9.6%")) +
    scale_fill_manual(name = "Sequence batch", values = c("black", "red")) +
    scale_color_manual(name = "", values = c("black", "red")) +
    plot_theme + plot_guide
ordi_clr_Seq

save_file <- paste("PCoA_clr_draft_sequenceBatch.svg", sep="")
ggsave(save_file, path = beta, scale = 1, width = 9.5, height = 5, units = c("in"), dpi = 300)

##==================== PhILR transformation (repeat this for subset of samples obj_ps_sub)
### Add clr to phyloseq object
ps_ilr <- phyloseq(
    otu_table(t(f.n0), taxa_are_rows = TRUE),
    sample_data(sample_data(obj_ps)),
    phy_tree(filt_tree),
    tax_table(tax_table(obj_ps))
)

### Prepare the data for PhILR
is.rooted(phy_tree(ps_ilr))
is.binary.tree(phy_tree(ps_ilr))
phy_tree(ps_ilr) <- makeNodeLabel(phy_tree(ps_ilr), method="number", prefix='n')
name.balance(phy_tree(ps_ilr), tax_table(ps_ilr), 'n10')
otu.table <- as.data.frame(otu_table(ps_ilr))
otu.table <- as.matrix(t(otu.table))
tree <- phy_tree(ps_ilr)
metadata <- sample_data(ps_ilr)
tax <- tax_table(ps_ilr)

### PhILR
gp.philr <- philr(otu.table, tree, part.weights='enorm.x.gm.counts', ilr.weights='blw.sqrt', return.all=FALSE)

### Ordinate and add metadata
ordi.scores <- as.data.frame(gp.pcoa$vectors)
ordi.scores$Sample_abbrev <- rownames(ordi.scores)
ordi_data <- merge(alpha.div.metadata2, ordi.scores, by = c('Sample_abbrev'))
ordi_data$Loc_sec <- gsub("_", " ", ordi_data$Loc_sec)
ordi_data$Loc_sec <- factor(ordi_data$Loc_sec, ordered = TRUE, levels = c("Amargosa Valley", "Ash Meadows", "Death Valley", "Frenchman and Yucca Flat", "Pahute Mesa", "Rainier Mesa", "Spring Mountains", "Oasis Valley"))

### Plot ordination (Figure S7)
#### Calculate the hulls for each group
hull <- ordi_data %>%
  group_by(Loc_sec) %>%
  slice(chull(Axis.1, Axis.2))

ordi_philr <- ggplot(ordi_data, aes(x=Axis.1, y = Axis.2)) +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + geom_vline(xintercept=0, size=.2, linetype = "dashed", color="black") +
    geom_polygon(data = hull, alpha = 0.1, aes(color=Loc_sec, fill=Loc_sec)) +
    geom_point(size=5, color="black", shape=21, aes(fill=Loc_sec)) +
    xlab(paste("PC1: 15.3%")) +
    ylab(paste("PC2: 11.4%")) +
    scale_fill_manual(name = "Location", values = loc_sec_color) +
    scale_color_manual(name = "", values = loc_sec_color) +
    plot_theme + plot_guide
ordi_philr

save_file <- paste("PCoA_philr_draft.svg", sep="")
ggsave(save_file, path = beta, scale = 1, width = 9.5, height = 5, units = c("in"), dpi = 300)

### Plot ordination (Figure S7)
#### Calculate the hulls for each group
hull <- ordi_data %>%
  group_by(Sequence_batch) %>%
  slice(chull(Axis.1, Axis.2))

ordi_philr_Seq <- ggplot(ordi_data, aes(x=Axis.1, y = Axis.2)) +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + geom_vline(xintercept=0, size=.2, linetype = "dashed", color="black") +
    geom_polygon(data = hull, alpha = 0.1, aes(color=Sequence_batch, fill=Sequence_batch)) +
    geom_point(size=3, color="black", shape=21, aes(fill=Sequence_batch)) +
    xlab(paste("PC1: 15.3%")) +
    ylab(paste("PC2: 11.4%")) +
    scale_fill_manual(name = "Sequence batch", values = c("black", "red")) +
    scale_color_manual(name = "", values = c("black", "red")) +
    plot_theme + plot_guide
ordi_philr_Seq

save_file <- paste("PCoA_philr_draft_sequenceBatch.svg", sep="")
ggsave(save_file, path = beta, scale = 1, width = 9.5, height = 5, units = c("in"), dpi = 300)

## Combine for Figure S7
both <- plot_grid(rpca_hull + theme(legend.position="none"),
                 rpca_hull_Seq + theme(legend.position="none"),
                 rpca_hull_sub + theme(legend.position="none"),
                 ordi_clr + theme(legend.position="none"),
                 ordi_clr_Seq + theme(legend.position="none"),
                 ordi_clr_sub + theme(legend.position="none"),
                 ordi_philr + theme(legend.position="none"),
                 ordi_philr_Seq + theme(legend.position="none"),
                 ordi_philr_sub + theme(legend.position="none"),
                 ncol=3, align = "v", axis="b")

save_file <- paste("Combo_for_supplementary_compositionalTransform.pdf", sep="")
ggsave(save_file, path = beta, plot = both, scale = 1, width = 12, height = 12, units = c("in"), dpi = 300)
