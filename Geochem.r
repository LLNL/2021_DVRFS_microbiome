#  Geochemistry Analysis of DVRFS
## Figure 2, Table S2, Table S3, Table S9, Table S10

## Load Packages
packages <- c("factoextra", "corrplot", "Hmisc", "funModeling", "PerformanceAnalytics", "naniar", "finalfit", "MissMech", "mice", "rms", "dlookr", "rcompanion", "rgdal", "EnvStats", "pscl", "zCompositions", "robCompositions", "compositions", "rgr", "energy", "missMDA", "ggfortify", "readxl")

### Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

### Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  BiocManager::install(packages[!installed_packages])
}

### Packages loading
invisible(lapply(packages, library, character.only = TRUE))

## Set Paths
path_geochem <- "<path to geochem directory>"
dir.create(quality, showWarnings = FALSE)
setwd(path_geochem)

## Load data
geochem <- read_excel("metadata_table_all_parameters_v2.3.xlsx")
geochem$ID_keep <- geochem$Sample_shortID
geochem <- data.frame(geochem, row.names = 2)
geochem <- geochem[ which(geochem$ID!='EC'), ] # remove Control sample "EC"
colnames(geochem)

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
plot_nomargins_y <- scale_y_continuous(expand = expansion(mult = c(0, 0)))
plot_nomargins_x <- scale_x_discrete(expand = expansion(mult = c(0, 0)))

plot_theme_boxplot <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 1, linetype = "solid"),
    panel.border = element_rect(colour="black", size=1, fill=NA),
    strip.background=element_rect(fill='white', colour='white'),
    strip.text = element_text(size=20, face="bold"),
    panel.grid.major = element_line(size = 0),
    panel.grid.minor = element_line(size = 0),
    axis.text = element_text(size=20, colour="black"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
    axis.title = element_text(face="bold", size=12),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position="right",
    legend.key = element_rect(fill = "white"),
    legend.text = element_text(size=20, colour="black"),
    legend.title = element_blank())
    
## Convert lat/long from decimal to UTM
geochem_latlong <- geochem[c(10:11)]
cord.dec = SpatialPoints(cbind(geochem_latlong$longitude, geochem_latlong$latitude), proj4string = CRS("+proj=longlat")) # Setting existing coordinate as lat-long system
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:32711")) # Transforming coordinate to UTM using EPSG=32711 UTM Zone=11N, nevada east)

### Plot to check if got correct lat/long conversion
par(mfrow = c(1, 2))
plot(cord.dec, axes = TRUE, main = "Lat-Long Coordinates", cex.axis = 0.95)
plot(cord.UTM, axes = TRUE, main = "UTM Coordinates", col = "red", cex.axis = 0.95)

cord.UTM <- as.data.frame(cord.UTM)
colnames(cord.UTM)[1] <- "longitude"
colnames(cord.UTM)[2] <- "latitude"
cord.UTM$ID_keep <- rownames(geochem)

geochem_UTM <- subset(geochem, select = -c(latitude, longitude))
geochem_UTM <- merge(geochem_UTM, cord.UTM, by="ID_keep")

## Summarize the data
describe_geochem <- df_status(geochem_UTM)

## Test Normality
geochem_num <- select_if(geochem_UTM, is.numeric) #obtain only columns with numeric values
geochem_num %>% select_if(~ !any(is.na(.))) -> geochem_num_noNA
rownames(geochem_num) <- geochem_UTM$ID_keep

### sort non-normal distribution variables by p-value
geochem_num %>%
    normality() %>%
    filter(p_value <= 0.01) %>%
    arrange(abs(p_value)) -> geochem_nonnorm

### sort normal distribution variables by p-value
geochem_num %>%
    normality() %>%
    filter(p_value >= 0.01) %>%
    arrange(abs(p_value)) -> geochem_norm

## Skewness
plot_geochem_num <- plot_num(geochem_num) + plot_theme
save_file_plot <- paste("geochem_histogram_summary.pdf", sep="")
ggsave(save_file_plot, path = path_geochem, scale = 1, width = 10, height = 10, units = c("in"), dpi = 300)

geochem_num_skew <- skewness(geochem_num)
#If the skewness of the predictor variable is 0, the data is perfectly symmetrical,
#If the skewness of the predictor variable is less than -1 or greater than +1, the data is highly skewed,
#If the skewness of the predictor variable is between -1 and -0.5 or between +1 and +0.5 then the data is moderately skewed,
#If the skewness of the predictor variable is -0.5 and +0.5, the data is approximately symmetric.

## Initial Correlation analysis
### Only keep variables with count < 50% for BDL and NA
drop_data <- names(geochem_num) %in% c("NH4_mgL","NO2_mgL","Mo_mgL","Zn_mgL","PO4_mgL", "Cu_mgL", "Cr_mgL", "Pb_mgL", "W_mgL")
geochem_num <- geochem_num[!drop_data]
colnames(geochem_num)

geochem_var <- varclus(as.matrix(geochem_num))
summary(geochem_var)
print(round(geochem_var$sim,2))
pdf(file = "geochem_collinearity_varclus_original_data.pdf")
plot(geochem_var)
dev.off()

geochem_var <- varclus(as.matrix(geochem_num3)) # repeat drop data for collinear variables
summary(geochem_var)
print(round(geochem_var$sim,2))
pdf(file = "geochem_collinearity_varclus_before_imputation.pdf")
plot(geochem_var)
dev.off()

## Obtain dataset for imputation
### Test using all variables for imputation. Only keep the numerical variables with ppm unit and <50% count BDL+NA; removed tritium because it causes problems for imputation
drop_data <- names(geochem_num) %in% c("depth_from_watertable_m","Dist_nuc_m","Temp_C","pH","Conductivity_uScm","longitude","latitude","NH4_mgL","NO2_mgL","Mo_mgL","Zn_mgL","PO4_mgL", "Cu_mgL", "Cr_mgL", "Pb_mgL", "W_mgL","tritium_mgL")
geochem_num3 <- geochem_num[!drop_data]
colnames(geochem_num3)

## Impute 0 and missing values with Zcompositions
### get detection limits (Table S4)
dl <- read_excel("metadata_table_all_parameters_v2.3.xlsx", sheet = "DL_R1")
dl <- as.matrix(dl[,-1])
colnames(dl) <- NULL

### Log-ratio EM algorithm (plus)
geochem_num3_lrEMplus_all300 <- lrEMplus(geochem_num3,dl=dl,rob=FALSE,ini.cov="multRepl",delta=0.65,tolerance=0.0001,max.iter=300,closure=10^6)
#No. iterations to converge:  84  (geochem_num3_lrEMplus_all300)

### Convert to a composition
geochem_num3_lrEMplus_acomp <- acomp(geochem_num3_lrEMplus_all300,detectionlimit=dl,total=1,parts=1:15)

## Transform the data clr and ilr
geochem_num3_lrEMplus_clr <- clr(geochem_num3_lrEMplus_acomp, ifclose = FALSE, ifwarn = TRUE)
geochem_num3_lrEMplus_ilr <- ilr(geochem_num3_lrEMplus_acomp)
geochem_num3_lrEMplus_ilrInv <- ilrInv(geochem_num3_lrEMplus_ilr, orig=geochem_num3_lrEMplus_acomp)

## Boxplot comparison of orig data and imputed data
### format the imputed data
geochem_num3_lrEMplus_boxplot <- geochem_num3_lrEMplus_all300
geochem_num3_lrEMplus_boxplot$ID <- geochem_UTM$ID_keep
geochem_num3_lrEMplus_boxplot$Loc_sec <- geochem_UTM$Loc_sec
geochem_num3_lrEMplus_boxplot <- gather(geochem_num3_lrEMplus_boxplot, Geochem, Measurement, -Loc_sec, -ID)
geochem_num3_lrEMplus_boxplot$Type <- c("Imputed")

### format the orig data
geochem_num3_boxplot <- geochem_num3
geochem_num3_boxplot$ID <- geochem_UTM$ID_keep
geochem_num3_boxplot$Loc_sec <- geochem_UTM$Loc_sec
geochem_num3_boxplot <- gather(geochem_num3_boxplot, Geochem, Measurement, -Loc_sec, -ID)
geochem_num3_boxplot$Type <- c("Original")

### combine the data
both <- rbind(geochem_num3_lrEMplus_boxplot, geochem_num3_boxplot)

### box plot samples separated
ggplot(both, aes(x = Loc_sec, y=Measurement, fill=Type)) + geom_boxplot() + facet_wrap(~Geochem, scales="free") +
    scale_fill_brewer(palette="Pastel2") + plot_theme_boxplot
ggsave("geochem_num3_lrEMplus_all300_boxplot_Loc_sec.svg", path = path_geochem, scale = 1, width = 20, height = 20, units = c("in"), dpi = 300)

### box plot samples together
ggplot(both, aes(x = Type, y=Measurement)) + geom_boxplot() + facet_wrap(~Geochem, scales="free") +
scale_fill_brewer(palette="Pastel2") + plot_theme_boxplot
ggsave("geochem_num3_lrEMplus_all300_boxplot.svg", path = path_geochem, scale = 1, width = 20, height = 20, units = c("in"), dpi = 300)

### Collinarity after imputation
geochem_num3_lrEMplus_all_varclus <- geochem_num3_lrEMplus_all300
geochem_num3_lrEMplus_all_varclus$Temp_C <- geochem_num$Temp_C
geochem_num3_lrEMplus_all_varclus$pH <- geochem_num$pH
geochem_num3_lrEMplus_all_varclus$depth_from_watertable_m <- geochem_num$depth_from_watertable_m
geochem_num3_lrEMplus_all_varclus$Conductivity_uScm <- geochem_num$Conductivity_uScm
geochem_num3_lrEMplus_all_varclus$longitude <- geochem_num$longitude
geochem_num3_lrEMplus_all_varclus$latitude <- geochem_num$latitude
geochem_num3_lrEMplus_all_varclus$tritium_mgL <- geochem_num$tritium_mgL

geochem_var <- varclus(as.matrix(geochem_num3_lrEMplus_all_varclus))
summary(geochem_var)
pdf("geochem_num3_lrEMplus_all300_collinearity_varclus.pdf")
plot(geochem_var)
dev.off()
print(round(geochem_var$sim,2))

drop_data <- names(geochem_num) %in% c("Conductivity_uScm", "latitude", "longitude", "Cl_mgL", "tritium_mgL", "depth_from_watertable_m","NH4_mgL","NO2_mgL","Mo_mgL","Zn_mgL","PO4_mgL", "Cu_mgL", "Cr_mgL", "Pb_mgL", "W_mgL") # remove variables that are collinear
geochem_num4 <- geochem_num[!drop_data]
colnames(geochem_num4)

geochem_var <- varclus(as.matrix(geochem_num4))
summary(geochem_var)
print(round(geochem_var$sim,2))
pdf(file = "varclus_geochem_num4.pdf")
plot(geochem_var)
dev.off()

drop_data <- names(geochem_num3_lrEMplus_all_varclus) %in% c("Conductivity_uScm", "latitude", "longitude", "Cl_mgL", "tritium_mgL", "depth_from_watertable_m","NH4_mgL","NO2_mgL","Mo_mgL","Zn_mgL","PO4_mgL", "Cu_mgL", "Cr_mgL", "Pb_mgL", "W_mgL") # remove variables that are collinear
geochem_num3_lrEMplus_all_varclus2 <- geochem_num3_lrEMplus_all_varclus[!drop_data]
colnames(geochem_num3_lrEMplus_all_varclus2)

geochem_var <- varclus(as.matrix(geochem_num3_lrEMplus_all_varclus2))
summary(geochem_var)
print(round(geochem_var$sim,2))
pdf(file = "varclus_geochem_num3_lrEMplus_all_varclus2.pdf")
plot(geochem_var)
dev.off()

## Plot Correlation
pdf(file = "varclus_geochem_num3_lrEMplus_all_varclus_correlation.pdf", width=20, height=20)
chart.Correlation(geochem_num3_lrEMplus_all_varclus2, histogram = TRUE, pch = 19)
dev.off()

## Obtain dataset for phyloseq
metadata_phyloseq <- geochem_num3_lrEMplus_all_varclus2

### match up the samples and values; label the rownames; add other variables categorical; format
geochem2 <- geochem
rownames(metadata_phyloseq) <- geochem2$Sample_abbrev[match(metadata_phyloseq$TIC_mgCL, geochem2$TIC_mgCL)]
metadata_phyloseq$TIC_mgCL <- as.numeric(metadata_phyloseq$TIC_mgCL)
metadata_phyloseq$tritium_BqL <- geochem$tritium_BqL[match(rownames(metadata_phyloseq), geochem$Sample_abbrev)]
metadata_phyloseq$tritium_BqL[is.na(metadata_phyloseq$tritium_BqL)] <- c(0.037) # modify NA to lowest detectable tritium for downstream analyses
metadata_phyloseq$tritium_BqL[metadata_phyloseq$tritium_BqL == 0] <- c(0.037) # modify zeroes to lowest detectable tritium for downstream analyses
metadata_phyloseq$tritium_BqL <- as.numeric(metadata_phyloseq$tritium_BqL)
metadata_phyloseq$Loc_sec <- geochem$Loc_sec[match(rownames(metadata_phyloseq), geochem$Sample_abbrev)]
metadata_phyloseq$Well_spring <- geochem$Well_spring[match(rownames(metadata_phyloseq), geochem$Sample_abbrev)]
metadata_phyloseq$rock_type <- geochem$rock_type[match(rownames(metadata_phyloseq), geochem$Sample_abbrev)]
metadata_phyloseq$Piper_group_ref <- geochem$Piper_group_ref[match(rownames(metadata_phyloseq), geochem$Sample_abbrev)]
metadata_phyloseq$Piper_group_ref <- factor(metadata_phyloseq$Piper_group_ref, ordered = TRUE, levels = c("Ca-Mg-HCO3", "Na-HCO3", "NaCl"))
colnames(metadata_phyloseq)

## PCA Plot
### Keep numeric variables; formatting
drop_data <- names(metadata_phyloseq) %in% c("Loc_sec","Well_spring","rock_type","Piper_group_ref")
metadata_phyloseq_drop <- metadata_phyloseq[!drop_data]
colnames(metadata_phyloseq_drop)
pca_metadata <- metadata_phyloseq_drop
pca_metadata <- pca_metadata %>% rename(TOC=TOC_mgCL, TIC = TIC_mgCL, Ca=Ca_mgL, Fe=Fe_mgL, K=K_mgL, Mg=Mg_mgL, Mn=Mn_mgL, Na=Na_mgL, Nitrate=NO3_mgL, Silica=SiO2_mgL, Sulfate = SO4_mgL, Sr=Sr_mgL, Cs=Cs_mgL, U=U_mgL, Temperature=Temp_C, Tritium=tritium_BqL)

### PCA, variation, and clustering
res.pca = princomp(pca_metadata, cor=TRUE)
var <- get_pca_var(res.pca) # extract results from pca
eig.val <- get_eigenvalue(res.pca)
eig.val

### plot a scree plot
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 25)) + plot_theme + plot_nomargins_y + plot_nomargins_x
ggsave("scree_plot_eigenvalues.TIFF", path = path_geochem, scale = 1, width = 5, height = 5, units = c("in"), dpi = 300)

### total contribution of variables: A reference dashed line is also shown on the barplot. This reference line corresponds to the expected value if the contribution were uniform.
fviz_contrib(res.pca, choice = "var", axes = 1:2, ylim = c(0, 25)) + plot_theme + plot_nomargins_y + plot_nomargins_x
ggsave("contribution_variables_dim1_dim2.TIFF", path = path_geochem, scale = 1, width = 5, height = 5, units = c("in"), dpi = 300)

fviz_contrib(res.pca, choice = "var", axes = 2:3, ylim = c(0, 25)) + plot_theme + plot_nomargins_y + plot_nomargins_x
ggsave("contribution_variables_dim2_dim3.TIFF", path = path_geochem, scale = 1, width = 5, height = 5, units = c("in"), dpi = 300)

fviz_contrib(res.pca, choice = "var", axes = 1, ylim = c(0, 25)) + plot_theme + plot_nomargins_y + plot_nomargins_x
ggsave("contribution_variables_dim1.TIFF", path = path_geochem, scale = 1, width = 5, height = 5, units = c("in"), dpi = 300)

fviz_contrib(res.pca, choice = "var", axes = 2, ylim = c(0, 25)) + plot_theme + plot_nomargins_y + plot_nomargins_x
ggsave("contribution_variables_dim2.TIFF", path = path_geochem, scale = 1, width = 5, height = 5, units = c("in"), dpi = 300)

### Plotting
#### create theme
pca_theme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 1, linetype = "solid"),
    panel.border = element_rect(colour="black", size=1, fill=NA),
    panel.grid.major = element_line(size = 0),
    panel.grid.minor = element_line(size = 0),
    strip.background=element_rect(fill='white', colour='white'),
    strip.text = element_text(size=15),
    axis.text = element_text(size=15, colour="black"),
    axis.title = element_text(face="bold", size=15),
    legend.position="right",
    legend.key = element_rect(fill = "white"),
    legend.title = element_text(size=15, colour="black"))
xlimit <- scale_x_continuous(breaks = c(-4,-2,0,2,4), limits = c(-4,4))
ylimit <- scale_y_continuous(breaks = c(-4,-2,0,2,4,6,8), limits = c(-4,9))

#### Get the data together
pca_sample_loadings <- res.pca$scores
pca_sample_loadings <- as.data.frame(pca_sample_loadings)
rownames(metadata_phyloseq) <- metadata_phyloseq$Sample_ID
metadata_phyloseq <- as.data.frame(metadata_phyloseq)
metadata_phyloseq$PC1 <- pca_sample_loadings$Comp.1[match(rownames(metadata_phyloseq), rownames(pca_sample_loadings))]
metadata_phyloseq$PC2 <- pca_sample_loadings$Comp.2[match(rownames(metadata_phyloseq), rownames(pca_sample_loadings))]

#### Figure 2A (Location)
metadata_phyloseq$Loc_sec <- gsub('_', ' ', metadata_phyloseq$Loc_sec)
metadata_phyloseq$Loc_sec <- factor(metadata_phyloseq$Loc_sec, ordered = TRUE, levels = c("Amargosa Valley", "Ash Meadows", "Death Valley", "Frenchman and Yucca Flat", "Pahute Mesa", "Rainier Mesa", "Spring Mountains", "Oasis Valley"))
loc_sec_color = c("Amargosa Valley" = "#043378", "Ash Meadows" = "#E50006", "Death Valley" = "#38AA31", "Frenchman and Yucca Flat" = "#1187A5", "Pahute Mesa" = "#7E478D", "Rainier Mesa" = "#FB9E7F", "Spring Mountains" = "#9B0020", "Oasis Valley" = "#9DA7A7")
plot_guide <- guides(fill = guide_legend(order=1, override.aes = list(shape = 21, alpha=1, size = 5)),
shape = "none",
color = "none")

figure_2A <- ggplot(metadata_phyloseq, aes(x=PC1, y = PC2)) +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + geom_vline(xintercept=0, size=.2, linetype = "dashed", color="black") +
    geom_point(size=3, color="black", shape=21, aes(fill=Loc_sec)) +
    scale_fill_manual(name = "Location", values = loc_sec_color) +
    scale_color_manual(name = "", values = loc_sec_color) +
    xlab("PC1: 22.1%") + ylab("PC2: 19.2%") +
    pca_theme + xlimit + ylimit + plot_guide +
    geom_text_repel(aes(label = Sample_ID), box.padding = 0.3, point.padding = 0.3, segment.color = 'grey50', max.overlaps = getOption("ggrepel.max.overlaps", default = 20))
figure_2A
ggsave("PCA_location.pdf", scale = 1, width = 7, height = 4, units = c("in"), dpi = 300)

#### Figure 2B (rock type, piper plot)
plot_guide <- guides(fill = guide_legend(order=1, override.aes = list(shape = 21, alpha=1, size = 5, color="black")),
shape = guide_legend(order=1, override.aes = list(alpha=1, size = 5, color="black")))

figure_2B <- ggplot(metadata_phyloseq, aes(x=PC1, y = PC2)) +
    geom_hline(yintercept=0, size=.2, linetype = "dashed", color="black") + geom_vline(xintercept=0, size=.2, linetype = "dashed", color="black") +
    geom_point(size=3, color="black", aes(fill=rock_type, shape=Piper_group_ref)) +
    scale_shape_manual(name = "Major-ion Chemistry", values = c(21, 22, 23)) +
    scale_fill_manual(name = "Rock type", values = c("#ae2012", "#a8dadc")) +
    xlab("PC1: 22.1%") + ylab("PC2: 19.2%") +
    pca_theme + xlimit + ylimit + plot_guide +
    geom_text_repel(aes(label = Sample_ID), box.padding = 0.3, point.padding = 0.3, segment.color = 'grey50', max.overlaps = getOption("ggrepel.max.overlaps", default = 20))
figure_2B
ggsave("PCA_piper_group_ref.pdf", scale = 1, width = 7, height = 4, units = c("in"), dpi = 300)

#### Biplot (Figure 2C)
figure_2C <- fviz_pca_biplot(res.pca,
    # Individuals
        geom.ind = c("point"),
        fill.ind = "white",
        col.ind = "white",
        pointshape = 1, pointsize = 3,
        mean.point = FALSE,
    # Variables
        col.var = "contrib",
        invisible = "quanti.sup",
    # Other
        legend.title = "Contribution (%)",
        repel = TRUE,
        title = NULL,
        labelsize = 5
    ) + pca_theme + xlimit + ylimit + xlab("PC1: 22.1%") + ylab("PC2: 19.2%") +
    scale_fill_gradient2(limits = c(0, 15), breaks=c(0,5,10,15)) + pca_guide
figure_2C
ggsave("PCA_biplot.pdf", scale = 1, width = 7, height = 5, units = c("in"), dpi = 300)

#### Combine Figure 2A,2B,2C
library("cowplot")
both <- plot_grid(figure_2A + theme(legend.position="none"),
                 figure_2B + theme(legend.position="none"),
                 figure_2C + theme(legend.position="none"),
                 ncol=3, align = "v", axis="b")

save_file <- paste("Combo_for_Figure2.svg", sep="")
ggsave(save_file, plot = both, scale = 1, width = 15, height = 5, units = c("in"), dpi = 300)

#### individual only - piper group - convex
fviz_pca_ind(res.pca,
    # Individuals
        geom.ind = c("point","text"),
        fill.ind = metadata_phyloseq$Piper_group3, col.ind = "black",
        pointshape = 21, pointsize = 3, labelsize = 3,
        palette = "jco",
        addEllipses = TRUE,
        ellipse.type = "convex",
        mean.point = FALSE,
    # Other
        legend.title = "Group",
        repel = TRUE
)  + xlimit + ylimit + pca_theme + xlab("PC1 (22.1%)") + ylab("PC2 (19.2%)")
ggsave("PCA_piperGroup_convex.svg", path = path_geochem, scale = 1, width = 7, height = 5, units = c("in"), dpi = 300)

#### individual only - piper group, no ellipse
fviz_pca_ind(res.pca,
    # Individuals
        geom.ind = c("point","text"),
            fill.ind = metadata_phyloseq$Piper_group3, col.ind = "black",
            pointshape = 21, pointsize = 3, labelsize = 3,
            palette = "jco",
            mean.point = FALSE,
        # Other
            legend.title = "Group",
            repel = TRUE
            )  + xlimit + ylimit + pca_theme + xlab("PC1 (22.1%)") + ylab("PC2 (19.2%)")
ggsave("PCA_piperGroup3.svg", path = path_geochem, scale = 1, width = 7, height = 5, units = c("in"), dpi = 300)

#### individual only - location
fviz_pca_ind(res.pca,
    # Individuals
        geom.ind = c("point","text"),
        fill.ind = metadata_phyloseq$Loc_sec, col.ind = "black",
        pointshape = 21, pointsize = 3, labelsize = 3,
        palette = "lancet",
        mean.point = FALSE,
    # Other
        legend.title = "Location",
        repel = TRUE
        )  + xlimit + ylimit + pca_theme + xlab("PC1 (22.1%)") + ylab("PC2 (19.2%)")
ggsave("PCA_locations.svg", path = path_geochem, scale = 1, width = 7, height = 5, units = c("in"), dpi = 300)

#### individual only - location type
fviz_pca_ind(res.pca,
    # Individuals
        geom.ind = c("point","text"),
        fill.ind = metadata_phyloseq$Well_spring, col.ind = "black",
        pointshape = 21, pointsize = 3, labelsize = 3,
        palette = "lancet",
        mean.point = FALSE,
    # Other
        legend.title = "Well or Spring",
        repel = TRUE
        ) + xlimit + ylimit + pca_theme + xlab("PC1 (22.1%)") + ylab("PC2 (19.2%)")
ggsave("PCA_well_spring.svg", path = path_geochem, scale = 1, width = 7, height = 5, units = c("in"), dpi = 300)

## Data Export (Table S2, Table S5)
des_geochem_num3 <- as.data.frame(describe(geochem_num3))
des_geochem_num3_lrEMplus_all300 <- as.data.frame(describe(geochem_num3_lrEMplus_all300))

write.csv(des_geochem_num3, file = "exploratory_description_geochem_num3.csv")
write.csv(des_geochem_num3_lrEMplus_all300, file = "exploratory_description_des_geochem_num3_lrEMplus_all300.csv")
write.csv(metadata_phyloseq, file = "imputed_values_geochem_num3_lrEMplus_all300.csv", row.names=TRUE)
