#  02. Filtering OTUs
## Following decontam tutorial: https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
## Table S11; Figure S1

## Phyloseq object
ps_use

## Check library sizes
df <- as.data.frame(sample_data(ps_use)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps_use)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_Control, shape=Sequence_batch)) + geom_point()

## Filtering out OTUs found in DNA extraction control
sample_data(ps_use)$is.neg <- sample_data(ps_use)$Sample_Control == "Control"
contamdf.prev97 <- isContaminant(ps_use, method="prevalence", neg="is.neg", threshold=0.1)

cairo_pdf(file.path(path, "06.phyloseq/threshold_remove_contaminants_decontamV1.4.0_otu97.pdf"), onefile = T)
hist(contamdf.prev97$p, 100, xlim = c(0,1))
dev.off()

table(contamdf.prev97$contaminant)

## Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps_use, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_Control == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_Control == "Sample", ps.pa)

write.table(contamdf.prev97,paste("contamdf.prev.otu97.csv", sep=""), sep=",", row.names = T)

## Manually inspect the table; also assign OTUs as contaminants if are eukaryotic sequences and OTUs without Kingdom assignment 

## Load in new table
fix_cont <- read.csv("contamdf.prev.otu97_manual.csv", row.names=1)
fix_cont$manual_contaminant <- as.logical(fix_cont$manual_contaminant)

head(fix_cont[ which(fix_cont$manual_contaminant=='TRUE'), ])
table(fix_cont$manual_contaminant)

## Make phyloseq object of presence-absence in negative controls and true samples
ps.pa97 <- transform_sample_counts(ps_use, function(abund) 1*(abund>0))
ps.pa.neg97 <- prune_samples(sample_data(ps.pa97)$Sample_Control == "Control", ps.pa97)
ps.pa.pos97 <- prune_samples(sample_data(ps.pa97)$Sample_Control == "Sample", ps.pa97)

## Make data.frame of prevalence in positive and negative samples
df.pa97 <- data.frame(pa.pos=taxa_sums(ps.pa.pos97), pa.neg=taxa_sums(ps.pa.neg97), contaminant=fix_cont$manual_contaminant)

## Prune out the contaminants
ps.noncontam97 <- prune_taxa(!fix_cont$manual_contaminant, ps_use)
ps.noncontam97

## Remove samples with only 0 counts
ps.noncontam97 <- prune_samples(sample_sums(ps.noncontam97) > 0, ps.noncontam97)
ps.noncontam97 <- filter_taxa(ps.noncontam97, function(x) sum(x) > 0, TRUE)
ps.noncontam97

## Remove DNA extraction control
ps.noncontam97.filt <- subset_samples(ps.noncontam97, Sample_Control !="Control")
ps.noncontam97.filt <- filter_taxa(ps.noncontam97.filt, function(x) sum(x) > 0, TRUE)
ps.noncontam97.filt

## Examine the prevalence of each taxa 
## (Following: https://ucdavis-bioinformatics-training.github.io/2017-September-Microbial-Community-Analysis-Workshop/friday/MCA_Workshop_R/phyloseq.html)
prev0 = apply(X = otu_table(ps.noncontam97.filt),
MARGIN = 1,
FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prev0,
TotalAbundance = taxa_sums(ps.noncontam97.filt),
tax_table(ps.noncontam97.filt))

prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps.noncontam97.filt, taxonomic.rank = "Phylum"))

plot_prev <- ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps.noncontam97.filt),color=Phylum)) +
geom_hline(yintercept = 0.02, alpha = 0.5, linetype = 2) +
geom_vline(xintercept = 10, alpha = 0.5, linetype = 2) +
geom_point(size = 2, alpha = 0.7) +
scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
facet_wrap(~Phylum) + theme(legend.position="none")

ggsave("prevalencethreshold.all.noncontam.otu97.pdf", path = path_phy, scale = 1, width = 15, height = 10, units = c("in"), dpi = 300)

### Zoom in to the plot
plot_prev + ylim(0,0.1) + xlim(0,50)
ggsave("prevalencethreshold.zoom.all.noncontam.otu97.pdf", path = path_phy, scale = 1, width = 15, height = 10, units = c("in"), dpi = 300)

### Check the plots. No need for prevalence filtering

## Rarefaction curve (Figure S1)
# Note: in phyloseq_inext script, make global variables out of resl and samplabs
rarefaction_inext <- phyloseq_inext(ps.noncontam97.filt, Q = 0, curve_type = "diversity", correct_singletons = FALSE, endpoint = NULL, knots = 40,multithread = FALSE, show_CI = TRUE, show_sample_labels = TRUE, show_plot = TRUE, justDF = FALSE, add_raw_data = TRUE)

pp <- ggplot(data = rarefaction_inext$data, aes(x = m, y = qD, group = SampleID)) +
    geom_line(data = resl$interpolated, linetype = "solid") +
    geom_line(data = resl$extrapolated, linetype = "dashed") +
    geom_point(data = resl$observed, size = 2) +
    geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL, color = NULL), alpha = 0.2) +
    facet_wrap(~Loc_sec, scales="free") +
    geom_text(data=samplabs, aes(label=SampleID, x=SampSize, y=MaxQD), size = 4, hjust = -0.5) +
    theme(panel.background = element_rect(fill = "white", colour = "black", size = 1, linetype = "solid"),
      panel.border = element_rect(colour="black", size=1, fill=NA),
      strip.background=element_rect(fill='white', colour='white'),
      strip.text = element_text(face="bold", size=10),
      panel.grid.major = element_line(size = 0),
      panel.grid.minor = element_line(size = 0),
      axis.text = element_text(face="bold", size=10, colour="black"),
      axis.title = element_text(face="bold", size=10),
      legend.position="bottom",
      legend.key = element_rect(fill = "white"),
      legend.title = element_blank()) +
    guides(color="none")

pp

ggsave("Rarefraction_curve_otu97.svg", path = path_phy, scale = 1, width = 10, height = 10, units = c("in"), dpi = 300)

## Filter out OTU with less than 10 reads
ps_filt97_10 <- phyloseq_filter_sample_wise_abund_trim(ps.noncontam97.filt, minabund = 10, rm_zero_OTUs = TRUE)
ps_filt97_10

write.table(tax_table(ps_filt97_10), file = "06.phyloseq/TAX_PHYLOSEQ_ps_filt97_great10reads.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(otu_table(ps_filt97_10), file = "06.phyloseq/OTU_PHYLOSEQ_ps_filt97_great10reads.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

after_filter <- as.data.frame(sample_sums(ps_filt97_10))
write.table(after_filter, file = "read_count_after_filtering_otu97_great10abund.csv", append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
