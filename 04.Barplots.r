#  04. Barplots
## Figure S5

## Prepare phyloseq object
filt_tree = read.tree("<path to tree.newick>")
filt_tree

DATA_PHYLOSEQ <- read_excel("<path to metadata>", na = "NA")
DATA_PHYLOSEQ_FIXED <- data.frame(DATA_PHYLOSEQ, row.names = 1)
DATA_PHYLOSEQ_FIXED$Sample_abbrev <- row.names(DATA_PHYLOSEQ_FIXED)
DATA_PHYLOSEQ_FIXED

obj_ps <- phyloseq(
    otu_table(otu_table(ps_filt97_10), taxa_are_rows = TRUE),
    sample_data(DATA_PHYLOSEQ_FIXED),
    phy_tree(filt_tree),
    tax_table(tax_table(ps_filt97_10))
)

obj_ps

## Directory
barplots <- file.path(paste(path_phy, "/barplots", sep=""))
dir.create(barplots, showWarnings = FALSE)

## Set up theme
plot_theme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 1, linetype = "solid"),
    panel.border = element_rect(colour="black", size=1, fill=NA),
    strip.background=element_rect(fill='white', colour='white', size = 0),
    strip.text = element_text(face="bold", size=15),
    panel.spacing.x=unit(2, "lines"),
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
plot_guides <- guides(colour=FALSE, fill = guide_legend(ncol=1))
plot_nomargins_y <- scale_y_continuous(expand = expansion(mult = c(0, 0)), labels = function(x) paste0(x*100, "%"))
plot_nomargins_x <- scale_x_discrete(expand = expansion(mult = c(0, 0)))

## Set up coloring
Phylum_color <- read_excel("<path to list of phylas and corresponding color>", sheet = "Sheet1", na = "NA")
Phylum_color$Color <- paste0("#", Phylum_color$Color)
colors <- distinct(Phylum_color, Phylum, Color)
pal <- colors$Color
names(pal) <- colors$Phylum
pal

## Function to plot all levels >1% at each taxonomic level by relative abundance (Figure S5)
taxrank <- c("Kingdom","Phylum","Class","Order")
func_plotbar <- llply(as.list(taxrank), function(level, x) {
    glom <- tax_glom(x, taxrank=level)
    trans <- transform_sample_counts(glom, function(OTU) OTU/sum(OTU) )
    dat <- psmelt(trans) # create dataframe from phyloseq object
    dat[[level]] <- as.character(dat[[level]]) # convert taxrank to a character vector from a factor
    dat <- dat[!(dat$Abundance == 0),]
    dat[[level]][dat$Abundance <= 0.01] <- 'Less than 1%' # find taxrank whose rel. abund. is less than 1%
    # group the data by taxrank, Sample, and Location
    columns = c("Sample_abbrev", "Loc_sec", "Well_spring", "rock_type", "Piper_group3", "Sequence_batch", "Loc_DV3", "HGU_DV3", "Re_Discharge_DV3")
    dat_grouped <- dat %>%
        group_by_at(vars(one_of(level, columns))) %>%
        summarise(grouped_Abundance = sum(Abundance))
    #set color palette to accommodate the number of genera
        colourCount = length(unique(dat_grouped[[level]]))
        getPalette = colorRampPalette(brewer.pal(8, "Accent"))
    print(colourCount)
    # manipulation of the taxrank list
        dat_grouped[[level]] <- str_replace(dat_grouped[[level]], "NA", "Unclassified") #rename "NA" to "Unclassified"
        dat_grouped[[level]] <- fct_relevel(dat_grouped[[level]], "Unclassified", after = Inf) #make "Unclassified" last
        dat_grouped[[level]] <- fct_relevel(dat_grouped[[level]], "Less than 1%", after = Inf) #make "Less than 1%" last
    # manipulation of location and sample names names
        dat_grouped$Loc_sec <- str_replace_all(dat_grouped$Loc_sec, "_", " ") #replace "_" with space
    # plot
        plot <- ggplot(data=dat_grouped, aes_string(x="Sample_abbrev", y="grouped_Abundance", fill=level)) +
            scale_fill_manual(values=getPalette(colourCount)) + # for plotting only phylum level, change values to object pal
            facet_grid(. ~ Loc_sec, scales="free", space="free") +
            labs(y = "Relative Abundance\n\n", fill = level) +
            geom_bar(aes(), stat="identity", position="fill", color="black", width=0.9) +
            plot_nomargins_y + plot_nomargins_x + plot_theme + plot_guides
    print(plot)
    # save plot
        save_file_plot <- paste("barplot.basic",type, ".", level,"great1perc.pdf", sep="") #change the file name if need to
        ggsave(save_file_plot, path = barplots, scale = 1, width = 40, height = 15, units = c("in"), dpi = 300)
        write.csv(dat_grouped, paste0("grouped_rel_abundances_great1perc",type,".",level,".csv"))
}, obj_ps)
