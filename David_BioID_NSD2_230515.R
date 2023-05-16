## script to analyze David's BioID data, and plot a volcano plot

library(tidyverse)
library(magrittr)
library(ggrastr)
library(ggrepel)
library(patchwork)
library(limma)

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
corfit = duplicateCorrelation(mat, design, block=technical_reps)

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




## filter just significant proteins
sig_ints =
  results %>%
  filter(adj.P.Val<p_thresh & abs_fc>lfc_thresh)





######################
#### VOLCANO PLOT ####
######################

to_label = "FBXO22"


# set up universal features of the plots
boxpad = .5
pointpad = .5
minlength = .01
legpos = "bottom"


## set up label column to use in geom_text_repel
volc =
  results %>%
  mutate(label = if_else(
    Gene %in% to_label,
    Gene,
    "nolabel"
  ))

## make minuslogp column
volc %<>%
  mutate(minuslogp=-log(adj.P.Val))

## add significance column
volc %<>%
  mutate(sig=if_else(adj.P.Val<p_thresh & abs_fc>lfc_thresh,"sig","notsig"))


# custom volcano plot
ggplot(volc,aes(logFC, minuslogp, label = label))+
  geom_hline(yintercept = -log10(p_thresh),linetype="dashed",size=.1)+
  geom_vline(xintercept = log2(lfc_thresh),linetype="dashed",size=.1)+
  geom_vline(xintercept = -log2(lfc_thresh),linetype="dashed",size=.1)+
  geom_point(data=subset(volc, sig!="sig"),alpha=.4,shape=".",colour="gray60")+
  geom_point(data=subset(volc, sig=="sig"&label=="nolabel"),alpha=.8,shape=".")+
  geom_point(data=subset(volc, label!="nolabel"),alpha=.95,shape=20,size=.2,colour="red")+
  # decreased proteins number label
  #annotate(geom="text", x=label_x_left, y=label_y, label=prot_nums[1,2],color="black",size=1)+
  # increased proteins number label
  #annotate(geom="text", x=label_x_right, y=label_y, label=prot_nums[2,2],color="black",size=1)+
  ggtitle("NSD2 miniTurbo UNC8732")+
  ylab("-log(p.adj)")+
  xlab("Log2FC")+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(size=.1),
    axis.ticks = element_line(size=.1),
    text=element_text(size=5),
    legend.key.size = unit(5,"mm"),
    title=element_text(size=2.5),
    legend.position = "bottom"
  ) +
  geom_text_repel(data          = subset(volc, label!="nolabel"),
                  colour="black",
                  size          = 1,
                  box.padding   = boxpad,
                  point.padding = pointpad,
                  max.overlaps = 40,
                  force         = 100,
                  segment.size  = 0.1,
                  min.segment.length = minlength,
                  segment.color = "grey50",
                  #direction     = "x")
  )

ggsave(paste0(base,"Volcano_NSD2.pdf"),
       width=1.4,
       height=1.3)
