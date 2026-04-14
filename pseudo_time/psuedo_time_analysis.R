# Pseudo time analysis
# renormalize the protein value to gte the reasonable value

# original meta file--> first file Basal/secretory --> get the corresponding protein matrix 

# Filter meta data
path <- '/Users/christina/Desktop/ProteinResearch/Demo/dat/'
meta <- read.csv(paste0(path,
                        '04_Gene_X_SingleCell_and_annotations/sc_protein_annotations.csv'),
                 row.names=1)
meta_sub <- meta %>%
  filter(Cell_Type %in% c('Basal','Secratory'))
head(meta_sub)
# filter protein matrix to these cells
prot <- read.csv(paste0(path,
                        '04_Gene_X_SingleCell_and_annotations/sc_protein_relative.csv'),
                 row.names=1)
prot <- as.matrix(prot)
head(prot)
prot_sub <- prot[, intersect(colnames(prot), meta_sub$ID)] #  ID is cell ID 
df <- prot_sub
df %>% distinct(colnames(), .keep_all = TRUE)
min_cells <- 200
prot_filt <- prot_sub[rowSums(is.na(prot_sub) == F) >= min_cells, ]

nrow(prot_filt)
nrow( prot_sub)
mean(is.na(prot_filt))

round(nrow(prot_filt)/nrow(prot_sub)*100, 1)

# imputation protein level  use median for imputation

prot_imp <- t(apply(prot_filt, 1, function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
}))

# correlation matrix and do PCA 
cor_mat <- cor(prot_imp, use = "pairwise.complete.obs") # calculate each cell coorelation among all the protein
# eigen decomposition vectors columns=PCs each column one PC
pca_res  <- eigen(cor_mat)
pc_scores <- pca_res$vectors
rownames(pc_scores) <- colnames(prot_imp) 
head(pc_scores)

pc_df <- as.data.frame(pc_scores[,1:5])
head(pc_df)
colnames(pc_df)<- paste0("PC", 1:5)
pc_df$ID <- rownames(pc_df)
pc_df$Cell_Type <- meta_sub$Cell_Type[match(pc_df$ID, meta_sub$ID)]

var_explained <- pca_res$values / sum(pca_res$values) * 100
elbow_df <- data.frame(PC = 1:20, var = var_explained[1:20])

ggplot(elbow_df, aes(x = PC, y = var)) +
  geom_point() + geom_line() +
  theme_classic(base_size = 14) +
  labs(x = "PC", y = "% Variance Explained", title = "Elbow Plot") # selec 5 PC as the final choice


#PCA plot colored by cell type
head(pc_df)

ggplot(pc_df,aes(x=PC1,y=PC2,color=Cell_Type))+
  geom_point(size = 1.5, alpha = 0.7) +
  scale_color_manual(values = c(Basal = "#B50202", Secratory = "#B606C4")) +
  theme_classic(base_size = 14) +
  ggtitle("PCA: Basal vs Secretory (renormalized)")

# 2 protein marker for Basal krt5 and secretory is Arg2

library(seqinr)
library(stringr)

Proc_fasta <- function(path){
  convert_mouse <- read.fasta(path, set.attributes=T, whole.header=T)
  convert_mouse <- names(convert_mouse)
  parse_row <- grep("GN=", convert_mouse, fixed=T)
  split_prot <- str_split(convert_mouse[parse_row], pattern=fixed("GN="))
  gene <- unlist(split_prot)[seq(2, 2*length(split_prot), 2)]
  prot <- unlist(split_prot)[seq(1, 2*length(split_prot), 2)]
  split_gene <- str_split(gene[parse_row], pattern=fixed(" "))
  split_gene <- unlist(split_gene)[seq(1, 3*length(split_gene), 3)]
  split_prot <- str_split(prot[parse_row], pattern=fixed("|"))
  split_prot <- unlist(split_prot)[seq(2, 3*length(split_prot), 3)]
  return(as.data.frame(cbind(split_prot, split_gene)))
}
convert <- Proc_fasta(paste0(path, 'Mouse.fasta'))
krt5_row  <- "Q922U2"                    
agr2_row  <- "O88312" 

df_marker <- data.frame(
  Cell_Type = meta_sub$Cell_Type,
  ID        = meta_sub$ID,
  Krt5      = as.numeric(prot_sub[krt5_row,  meta_sub$ID]),
  Agr2      = as.numeric(prot_sub[agr2_row,  meta_sub$ID])
)
# color PCA by marker abundance 

ggplot(pc_df,aes(x=PC1,y=PC2,color=df_marker$Krt5))+
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_viridis_c(option = "magma", name = "Krt5") +
  theme_classic(base_size = 14) +
  ggtitle("PCA colored by Krt5 (Basal marker)")

ggplot(pc_df,aes(x=PC1,y=PC2,color=df_marker$Agr2))+
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_viridis_c(option = "magma", name = "Agr2") +
  theme_classic(base_size = 14) +
  ggtitle("PCA colored by Agr2 (Secretory marker)")

# Pseudotime Analysis
library(slingshot)

# cluster label
pc_mat <- as.matrix(pc_df[, paste0("PC", 1:5)]) # extract the PC from the overall matrix 
rownames(pc_mat) <- pc_df$ID
stopifnot(length(cl) == nrow(pc_mat))
cl <- factor(pc_df$Cell_Type) 
sce <- slingshot(
  data=pc_mat,
  clusterLabels = cl,
  start.clus    = 'Basal',
  end.clus      = 'Secratory'
)

# extract the psedotime  close to Basal --> the psuedotime small vice versa

pseudotime <- slingPseudotime(sce)[, 1]
pc_df$pseudotime <- pseudotime
head(pc_df)

# pca colored by psedotime
ggplot(pc_df,aes(x=PC1,y=PC2,color=pseudotime))+
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_viridis_c(option = "plasma", name = "pseudotime") +
  theme_classic(base_size = 14) +
  ggtitle("PCA colored by pseudotime")

# marker abundance and pusedotime  scatter
library(tidyr)
marker_pt <- pc_df %>% 
  select(ID, pseudotime, Cell_Type)
marker_pt$Krt5<- df_marker$Krt5
marker_pt$Arg2<- df_marker$Agr2
head(marker_pt)
marker_pt <- marker_pt %>%
  pivot_longer(cols = c(Krt5, Arg2),
               names_to = "protein", values_to = "abundance")

ggplot(marker_pt, aes(x = pseudotime, y = abundance, color = Cell_Type)) +
  geom_point(size = 0.8, alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "black", linewidth = 0.8) +
  facet_wrap(~protein, scales = "free_y") +
  scale_color_manual(values = c(Basal = "#B50202", Secratory = "#B606C4")) +
  theme_classic(base_size = 14) +
  labs(x = "Pseudotime", y = "Protein Abundance",
       title = "Marker proteins along pseudotime trajectory")

# principle curve 
plot(pc_mat[,1], pc_mat[,2],
     col = ifelse(pc_df$Cell_Type == "Basal", "#B50202", "#B606C4"),
     pch = 16, cex = 0.5)
curves <- slingCurves(sce)
for(i in seq_along(curves)){
  curve_coords <- curves[[i]]$s  
  lines(curve_coords[,1], curve_coords[,2], lwd = 2, col = "black")
}



ggplot(marker_pt, aes(x = pseudotime, y = abundance, color = Cell_Type)) +
  geom_point(size = 0.8, alpha = 0.3) +
  geom_smooth(aes(color = Cell_Type, fill = Cell_Type),  
              method = "loess",
              se     = TRUE,
              alpha  = 0.2,
              linewidth = 0.8) +
  facet_wrap(~protein, scales = "free_y") +
  scale_color_manual(values = c(Basal = "#B50202", Secratory = "#B606C4")) +
  scale_fill_manual(values  = c(Basal = "#B50202", Secratory = "#B606C4")) +
  theme_classic(base_size = 14) +
  labs(x = "Pseudotime", y = "Protein Abundance",
       title = "Marker proteins along pseudotime trajectory")


# principle curve with the psuedo time color 
library(RColorBrewer)


pt_vals  <- pc_df$pseudotime
pt_color <- colorRampPalette(c("#440154", "#31688E", 
                               "#35B779", "#FDE725"))(100)

pt_idx   <- cut(pt_vals, breaks = 100, labels = FALSE)

plot(pc_mat[,1], pc_mat[,2],
     col = pt_color[pt_idx],
     pch = 16, cex = 0.6,
     xlab = "PC1", ylab = "PC2",
     main = "Principal curve colored by pseudotime")

# principal curve
curves <- slingCurves(sce)
for(i in seq_along(curves)){
  curve_coords <- curves[[i]]$s
  lines(curve_coords[,1], curve_coords[,2], lwd = 2.5, col = "black")
}

legend("topright",
       legend = c("Early (0)", "Late (0.22)"),
       col    = c("#440154", "#FDE725"),
       pch    = 16, bty = "n")
