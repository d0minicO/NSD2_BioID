## script to analyze David's BioID data, plot a volcano plot, do clustering

library(tidyverse)
library(magrittr)
library(limma)
library(edgeR)
library(umap)

#################
#### OPTIONS ####
#################

options(max.print=100)

################
#### inputs ####
################

# specify root folder that contains the data
base = "C:/Users/dowens/OneDrive/Postdoc/Projects/NSD2/David/"

# specify data colnames
cnames =
  c(
    "Gene",
    "Gene_description",
    "DMSO_A1",
    "DMSO_A2",
    "DMSO_B1",
    "DMSO_B2",
    "DMSO_C1",
    "DMSO_C2",
    "DMSO_AVE",
    "DMSO_SD",
    "BFDR_DMSO",
    "UNC8732_A1",
    "UNC8732_A2",
    "UNC8732_B1",
    "UNC8732_B2",
    "UNC8732_C1",
    "UNC8732_C2",
    "UNC8732_AVE",
    "UNC8732_SD",
    "BFDR_UNC8732"
  )

# load the table
dat = read_tsv(paste0(base,"Data_raw.txt"),col_names=cnames)

# SAINT BFDR threshold to use
bfdr_thresh = 0.01

# p threshold for limma
p_thresh = .05

# logFC threshold for limma
lfc_thresh = log2(1.5)



#######################
#### DATA CLEANING ####
#######################

# keep useful cols
dat %<>%
  dplyr::select(-Gene_description,-contains("AVE"),-contains("SD")) %>%
  dplyr::select(Gene,starts_with("DMSO"),starts_with("UNC8732"),contains("BFDR"))

# get into long format
long = gather(dat,sample,SPC,DMSO_A1:UNC8732_C2)


## get a new column just for the lowest (best) BFDR for each protein
## filter based on lower than 0.01
long %<>%
  mutate(BFDR_DMSO_new = as.numeric(BFDR_DMSO),
         BFDR_UNC8732_new = as.numeric(BFDR_UNC8732)) %>%
  dplyr::select(-BFDR_DMSO,-BFDR_UNC8732) %>%
  dplyr::rename(BFDR_DMSO=BFDR_DMSO_new,
                BFDR_UNC8732=BFDR_UNC8732_new) %>%
  mutate(BFDR_min = case_when(
    is.na(BFDR_DMSO) & is.na(BFDR_UNC8732) ~ 10, # dummy value of 10, will be filtered out
    is.na(BFDR_DMSO) & !is.na(BFDR_UNC8732) ~ BFDR_UNC8732,
    !is.na(BFDR_DMSO) & is.na(BFDR_UNC8732) ~ BFDR_DMSO,
    BFDR_DMSO>BFDR_UNC8732 ~ BFDR_UNC8732,
    BFDR_DMSO<BFDR_UNC8732 ~ BFDR_DMSO,
    BFDR_DMSO==BFDR_UNC8732 ~ BFDR_DMSO
  )) %>%
  filter(BFDR_min<=bfdr_thresh)


## get into wide format
wide =
  long %>%
  dplyr::select(Gene,sample,SPC) %>%
  spread(key=sample,value=SPC)

## save the wide format filtered table
write.table(wide,
            paste0(base,"Filtered_wide_format.tsv"),
            col.names = T,
            row.names = F,
            quote=F,
            sep="\t")


#######################
#### NORMALIZATION ####
#######################

## perform voom normalization to prepare data for linear modelling

# Voom (Variance modeling at the observation level) normalization
# works by estimating the mean-variance relationship of the log-counts
# generating a precision weight for each observation
# and modifying the data to stabilize the variance across different levels of intensity.
# This transformation is necessary for analysis with limma
# because it allows limma's linear modeling and
# empirical Bayes moderation methods, originally designed for microarray data,
# to be applied effectively to count data, such as from RNA-seq or BioID spectral counts
# by making the variance of the data more homoscedastic
# and thus fulfilling the assumptions of the linear models used by limma.

# get wide data as numeric matrix
mat = 
  wide %>%
  column_to_rownames("Gene") %>%
  as.matrix()

# Turn your data matrix into a DGEList object (required for calcNormFactors)
y = DGEList(counts=mat)

# Calculate normalization factors to scale the raw library sizes
y = calcNormFactors(y)

# Estimate the common, trended and tagwise dispersion
y = estimateDisp(y)

# Convert to log2-counts per million, with prior.count=2
logCPM = cpm(y, log=TRUE, prior.count=2)

# Perform the voom transformation
v = voom(logCPM, plot=TRUE)

mat = v$E

# Check the first few rows and columns of your transformed data
head(mat)

#mat = log2(mat+1)


##################################################
#### TEST DIFFERENTIAL ENRICHMENT USING LIMMA ####
##################################################

# By including technical replicates as a blocking factor
# the model accounts for this within-group correlation
# reducing the residual variability and potentially
# increasing the power to detect differences between the drug and control treatments

# Design matrix
treatments = rep(c("DMSO", "UNC8732"), each=6)
biological_reps = rep(rep(1:3, each=2), 2)
technical_reps = rep(rep(1:2, 3), 2)
design = model.matrix(~0 + factor(paste(treatments, biological_reps, sep="_")))
colnames(design) = levels(factor(paste(treatments, biological_reps, sep="_")))

# Correlation between technical replicates
corfit = limma::duplicateCorrelation(mat, design, block=technical_reps)

# Update the design matrix
fit = lmFit(mat, design, block=technical_reps, correlation=corfit$consensus)

# Apply contrasts to compare UNC8732 vs DMSO
contrast.matrix = makeContrasts(UNC8732_vs_DMSO = UNC8732_1 + UNC8732_2 + UNC8732_3 - DMSO_1 - DMSO_2 - DMSO_3, levels=design)
fit2 = contrasts.fit(fit, contrast.matrix)

# Perform eBayes moderation
fit2 = eBayes(fit2)

# Get the full table of results
results = topTable(fit2,coef="UNC8732_vs_DMSO",adjust="BH",n=Inf)
head(results)

# tidy table
results %<>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  as_tibble() %>%
  mutate(abs_fc = abs(logFC)) %>%
  arrange(desc(abs_fc))

## save the results
write.table(results,
            paste0(base,"Limma_results.tsv"),
            col.names = T,
            row.names = F,
            quote=F,
            sep="\t")




## filter just significant proteins by p value
sig_ints =
  results %>%
  filter(adj.P.Val<p_thresh)



## filter top n proteins in non-normalized table
mat2 =
  as.data.frame(mat) %>%
  rownames_to_column("Gene") %>%
  filter(Gene %in% sig_ints$Gene) %>%
  dplyr::select(-Gene) %>%
  as.matrix()


########...######.....###...
##.....##.##....##...##.##..
##.....##.##........##...##.
########..##.......##.....##
##........##.......#########
##........##....##.##.....##
##.........######..##.....##


## do the PCA on selected peaks
pc = prcomp(t(mat2), scale. = T, center=T)

## set up values and sample labels to plot using raw ggplot
pc_vals =
  pc$x %>%
  data.frame() %>%
  rownames_to_column("Sample") %>%
  separate(Sample,into=c("Condition","rep"),sep="_")



## get the variance explained for each component
var_explained <- pc$sdev^2/sum(pc$sdev^2)

# plot PC1 and PC2
pc_vals %>%
  ggplot(aes(x=PC1,y=PC2,col=Condition))+
  #geom_point(colour = "black", size = 2,aes(shape=rep)) +
  geom_point(size=1)+
  scale_colour_manual(values=c("#0044aaff","#bf0000ff"))+
  #ggtitle("PCA plot of proteomics","on 455 DEPs")+
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme_bw()+
  theme(panel.grid=element_blank(),
        panel.border = element_rect(size=.1),
        axis.ticks = element_line(size=.1),
        axis.text = element_text(size=6,colour="black",face="bold"),
        axis.title = element_text(size=6,colour="black",face="bold"),
        text=element_text(size=6),
        legend.key.size = unit(.5,"cm"))


ggsave(paste0(base,"PCA.pdf"),
       width=3,
       height=2)








##.....##.##.....##....###....########.
##.....##.###...###...##.##...##.....##
##.....##.####.####..##...##..##.....##
##.....##.##.###.##.##.....##.########.
##.....##.##.....##.#########.##.......
##.....##.##.....##.##.....##.##.......
#######..##.....##.##.....##.##.......

## non transformed data

## change default n_neighbours to allow umap on smaller sample size
custom.settings = umap.defaults
custom.settings$n_neighbors = ncol(mat2)-1


seed = 3
set.seed(seed)

#umap = umap::umap(as.matrix(t(mat)))

# custom config needed
umap = umap::umap(as.matrix(t(mat2)),
                  config = custom.settings)

scores = 
  data.frame(umap$layout) %>%
  rownames_to_column("sample") %>%
  #mutate(sample=gsub("BioID2_GID4","BioID2-GID4",sample)) %>%
  separate(sample, into=c("Condition","rep"),sep="_")



## umap with labelled tissue types
ggplot(data = scores, aes(x = X1, y = X2)) + 
  geom_point(aes(colour = Condition),size=1) + 
  scale_colour_manual(values=c("#0044aaff","#bf0000ff"))+
  labs(x="UMAP1",y="UMAP2")+
  scale_x_reverse()+
  theme_bw() + 
  theme(panel.grid=element_blank(),
        panel.border = element_rect(size=.1),
        axis.ticks = element_line(size=.1),
        axis.text = element_text(size=6,colour="black",face="bold"),
        axis.title = element_text(size=6,colour="black",face="bold"),
        text=element_text(size=6),
        legend.key.size = unit(.5,"cm"))


ggsave(filename=paste0(base,"UMAP.pdf"),
       width=3,
       height=2)