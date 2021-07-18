# 06. Beta Diversity
## Figure 3, Figure S6, Figure S7, Figure S8, Table 1

## Set Path
beta <- file.path(paste(path_phy, "/Beta_diversity", sep=""))
dir.create(beta, showWarnings=FALSE)
setwd(beta)

## Phyloseq object
obj_ps

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

## NMDS Weighted UniFrac
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

### Plot NMDS wUniFrac (Figure 3B, Figure S6B)
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

### Plot with rock type (Figure S7A)
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

### Plot with Piper group (Figure S9B)
hull <- ordi_data %>%
  group_by(Piper_group3) %>%
  slice(chull(MDS1, MDS2))

nmds_wuni_piper <- ggplot(ordi_data, aes(x=MDS1, y = MDS2)) +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + geom_vline(xintercept=0, size=.2, linetype = "dashed", color="black") +
    geom_polygon(data = hull, alpha = 0.1, aes(color=Piper_group3, fill=Piper_group3)) +
    geom_point(size=3, color="black", shape=21, aes(fill=Piper_group3)) +
    xlab(paste("NMDS1")) +
    ylab(paste("NMDS2")) +
    scale_fill_lancet(name = "Piper group") +
    scale_color_lancet(name = "") +
    plot_theme + plot_guide

### Plot with Sequence Batch (Figure S7C)
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

### Plot with Sequence Batch (Figure S7D)
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

### Plot with Location type (Figure S7E)
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

### Plot with TOC (Figure S7F)
nmds_wuni_toc <- ggplot(ordi_data, aes(x=MDS1, y = MDS2)) +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + geom_vline(xintercept=0, size=.2, linetype = "dashed", color="black") +
    geom_point(shape=21, size = 5, color="black", aes(fill=TOC_mgCL)) +
    xlab(paste("NMDS1")) +
    ylab(paste("NMDS2")) +
    scale_fill_gradient2(name = "TOC (mg-C/L)", breaks = c(0, 20, 40), limits = c(0,40), low="#f72585", mid="#cdb4db", high="#023e8a", midpoint=20) +
    scale_color_lancet(name = "") +
    plot_theme + plot_guide

### Plot with temperature (Figure S7G)
plot_guide <- guides(fill = guide_colourbar(frame.linewidth = 1, frame.colour = "black", ticks = TRUE, ticks.colour = "black", ticks.linewidth = 1))

nmds_wuni_temp <- ggplot(ordi_data, aes(x=MDS1, y = MDS2)) +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + geom_vline(xintercept=0, size=.2, linetype = "dashed", color="black") + 
    geom_point(shape=21, size = 5, color="black", aes(fill=Temp_C)) +
    xlab(paste("NMDS1")) +
    ylab(paste("NMDS2")) +
    scale_fill_gradient2(name = "Temperature (ºC)", breaks = c(5, 20, 40, 60), limits = c(5,60), low="#f72585", mid="#cdb4db", high="#023e8a", midpoint=30) + 
    scale_color_lancet(name = "") + 
    plot_theme + plot_guide

### Plot with Depth (Figure S7H)
nmds_wuni_depth <- ggplot(ordi_data, aes(x=MDS1, y = MDS2)) +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + geom_vline(xintercept=0, size=.2, linetype = "dashed", color="black") + 
    geom_point(shape=21, size = 5, color="black", aes(fill=Depth_sampling_m)) +
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

## NMDS UniFrac
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

## Plot both NMDS
both <- plot_grid(nmds_uni_hull + theme(legend.position="none"),
                 nmds_whull + theme(legend.position="none"), 
                 ncol=2, align = "v", axis="b")
save_file <- paste("wuni.unifrac.nmds.svg", sep="")
ggsave(save_file, path = beta, plot = both, scale = 1, width = 10, height = 5, units = c("in"), dpi = 300)

both_label <- plot_grid(nmds_uni_hull_label + theme(legend.position="none"),
                 nmds_whull_label + theme(legend.position="none"), 
                 ncol=2, align = "v", axis="b")
save_file <- paste("wuni.unifrac.nmds.label.svg", sep="")
ggsave(save_file, path = beta, plot = both_label, scale = 1, width = 10, height = 5, units = c("in"), dpi = 300)

## Statistics
### https://microbiome.github.io/tutorials/PERMANOVA.html
### https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/multivariate-comparisons-of-microbial-community-composition.html

### Obtain the distance matrix
uwUF.dist = UniFrac(obj_ps_rel, weighted=FALSE, normalized=TRUE)
head(uwUF.dist)

wUF.dist = UniFrac(obj_ps_rel, weighted=TRUE, normalized=TRUE)
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
adonis2(uwUF.dist~metadata_stats$Loc_sec+metadata_stats$Depth_sampling_m+metadata_stats$Well_spring+metadata_stats$Sequence_batch+metadata_stats$Sampling_method, permutations=999, by="margin")
adonis2(wUF.dist~metadata_stats$Loc_sec+metadata_stats$Depth_sampling_m+metadata_stats$Well_spring+metadata_stats$Sequence_batch+metadata_stats$Sampling_method, permutations=999, by="margin")

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

## DEICODE RPCA
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

#### Calculate the hulls for each group
hull <- rpca_data_num %>%
  group_by(Loc_sec) %>%
  slice(chull(PC1, PC2))

### Plot ordination (Figure S8)
rpca_hull <- ggplot(rpca_data_num, aes(x=PC1, y = PC2)) +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + geom_vline(xintercept=0, size=.2, linetype = "dashed", color="black") + 
    geom_polygon(data = hull, alpha = 0.1, aes(color=Loc_sec, fill=Loc_sec)) +
    geom_point(size=3, color="black", shape=21, aes(fill=Loc_sec)) +
    xlab(xaxis_text) +
    ylab(yaxis_text) +
    scale_fill_lancet(name = "Location") + 
    scale_color_lancet(name = "") + 
    plot_theme + plot_guide

save_file <- paste("RPCA_draft.svg", sep="")
ggsave(save_file, path = deicode, scale = 1, width = 9.5, height = 5, units = c("in"), dpi = 300)

## Statistics on DEICODE RPCA
### Import distance matrix
rpca_dist <- read_qza("distance_auto_components.qza")
rpca_dist.obj <- rpca_dist$data
rpca_dist.obj <- as.dist(rpca_dist.obj, diag = FALSE, upper = FALSE)

### Adonis
#### location-type variables
adonis2(rpca_dist.obj~metadata_stats$Loc_sec+metadata_stats$Depth_sampling_m+metadata_stats$Well_spring+metadata_stats$Sequence_batch+metadata_stats$Sampling_method, permutations=999, by="margin")
#### geochemical-type variables
adonis2(rpca_dist.obj~metadata_stats$Temp_C+metadata_stats$rock_type+metadata_stats$pH+metadata_stats$tritium_BqL_fct+metadata_stats$TOC_mgCL+metadata_stats$Sequence_batch+metadata_stats$Sampling_method, permutations=999, by="margin")