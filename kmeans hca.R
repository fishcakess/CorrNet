library(bcdstats)
library(tidyverse)
library(factoextra)
library(paran)
library(clustertend)
library(cluster)
library(dplyr)
library(dendextend)

library(jsonlite)

#load data 

gene_effect_df <- read.csv('median_gene_df.csv', row.names = 1)

pathwayDF <- read.csv('pathwayDF.csv', row.names = 1, check.names = FALSE)

#scale pathwayDF for "z scores"

pathwayDF <- scale(pathwayDF,center = TRUE, scale = TRUE)

#labels

Metabolism <- read.delim('Metabolism.txt')
Metabolism <- as.vector(Metabolism$Metabolism)

Genetic_Information_Processing <- read.delim('Genetic Information Processing.txt')
Genetic_Information_Processing <- as.vector(Genetic_Information_Processing$Genetic.Information.Processing)

Environmental_Information_Processing <- read.delim('Environmental Information Processing.txt')
Environmental_Information_Processing <- as.vector(Environmental_Information_Processing$Environmental.Information.Processing)

Cellular_Processes <- read.delim('Cellular Processes.txt')
Cellular_Processes <- as.vector(Cellular_Processes$Cellular.Processes)

mypathways <- c(Metabolism, Genetic_Information_Processing, Environmental_Information_Processing, Cellular_Processes)
mypathways <- as.factor(mypathways)

#generate my pathwayDF.csvs & corrmatrix
#moved them into blocks of 3 

#all pathways corr matrix

bcdcor <- adjust.corr(pathwayDF,type = c('pearson'),use = c('complete.obs'),adjust = 'BH')
corrmatrix <- as.data.frame(bcdcor$R$r)
qvalues <- as.data.frame(bcdcor$R$P)

squared_euclidean <- 1 - corrmatrix
euclidean <- sqrt(squared_euclidean)

#get corrhistogram

x <- data.frame(lower.triangle(as.matrix(corrmatrix)))

y <- data.frame(y = unlist(x))

y <- y[y$y != 0 & y$y != 1, ]

pdf('all pathways corr histogram.pdf')
hist(y,main = 'Meta, CP, EIP, GIP Correlation Histogram',xlab = 'Pearson r Value')
dev.off()

#HCA

m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

#function to compute agglomerative coefficient
ac <- function(x) {
  agnes(euclidean, method = x)$ac
}

sapply(m, ac)


#wards still best - use euclidean matrix for clustering

#gap stat 

gap_stat <- clusGap(euclidean, FUN=hcut, hc_func="hclust", hc_method="ward.D", K.max = 20, B = 500)

#removed isdiss=FALSE
pdf(file = 'all_gapstat_hca.pdf')
fviz_gap_stat(gap_stat)
dev.off()

#19 for all pathways

hc <- hclust(as.dist(euclidean), method = "ward.D2" )

# Cut tree into 4 groups

sub_grp <- cutree(hc, k = 19)


pdf(file = 'allpathwaysHCAdend.pdf',width = 15,height = 10)
plot(hc, cex = 0.6,main = '',xlab='',sub = '')
rect.hclust(hc, k = 19, border = 2:5)
dev.off()
fviz_cluster(list(data = corrmatrix, cluster = sub_grp),repel = TRUE)

clusters <- cutree(hc,k = 19)

clusters <- tibble(name = names(clusters), value = clusters)

all_hca_json <- jsonlite::toJSON(split(clusters$name, clusters$value))

write(all_hca_json, "all_hca.json")

# k medoids - kmeans

gap_stat <- clusGap(euclidean,FUN = pam,  K.max = 20, B = 500)
pdf(file = 'all_gapstat_pam')
fviz_gap_stat(gap_stat)
dev.off()

pamclusters <- pam(euclidean, diss=TRUE, k = 17)

pamclusters$data <- corrmatrix
pdf(file = 'all_pam_clusters.pdf')
fviz_cluster(pamclusters,repel=TRUE,ellipse = TRUE,geom='point',main = 'Metabolism, CP, EIP, GIP PAM Clusters')
dev.off()

fviz_silhouette(pamclusters)

pamclusterid <- pamclusters$clustering

pamclusterid <- tibble(name = names(pamclusterid), value = pamclusterid)

all_pam_json <- jsonlite::toJSON(split(pamclusterid$name, pamclusterid$value))

write(all_pam_json, "all_pam.json")


#non metabolism

nonmetabolism_pathways <- c(Genetic_Information_Processing,Environmental_Information_Processing,Cellular_Processes)
nonmetabolism_pathwayDF <- pathwayDF[,nonmetabolism_pathways]

nonmetabolism_pathwayDF <- scale(nonmetabolism_pathwayDF,center = TRUE, scale = TRUE)


bcdcor <- adjust.corr(nonmetabolism_pathwayDF,type = c('pearson'),use = c('complete.obs'),adjust = 'BH')
corrmatrix <- as.data.frame(bcdcor$R$r)
qvalues <- as.data.frame(bcdcor$R$P)

squared_euclidean <- 1 - corrmatrix
euclidean <- sqrt(squared_euclidean)



#histogram 

x <- data.frame(lower.triangle(as.matrix(corrmatrix)))

y <- data.frame(y = unlist(x))

y <- y[y$y != 0 & y$y != 1, ]

pdf('no metabolism corr histogram.pdf')
hist(y,main = 'CP, EIP, GIP Correlation Histogram',xlab = 'Pearson r Value')
dev.off()



#write.csv(corrmatrix,'non-metabolism.csv')




m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

#function to compute agglomerative coefficient
ac <- function(x) {
  agnes(euclidean, method = x)$ac
}

sapply(m, ac)


#wards still best - use euclidean matrix for clustering

#gap stat 

gap_stat <- clusGap(euclidean, FUN=hcut, hc_func="hclust", hc_method="ward.D", K.max = 20, B = 500)

#removed isdiss=FALSE
pdf(file='nometa_gapstat_hca.pdf')
fviz_gap_stat(gap_stat)
dev.off()
#19 for all pathways

hc <- hclust(as.dist(euclidean), method = "ward.D2" )

# Cut tree into 4 groups
sub_grp <- cutree(hc, k = 18)
pdf(file='nometabolismHCAdend.pdf',width=15,height=10)
plot(hc, cex = 0.6,main='',xlab='',sub = '')
rect.hclust(hc, k = 18, border = 2:5)
dev.off()

fviz_cluster(list(data = corrmatrix, cluster = sub_grp),repel = TRUE)

clusters <- cutree(hc,k = 18)

clusters <- tibble(name = names(clusters), value = clusters)

no_met_hca_json <- jsonlite::toJSON(split(clusters$name, clusters$value))

write(no_met_hca_json, "no_met_hca.json")

# k medoids


gap_stat <- clusGap(euclidean, FUN=pam, K.max = 20, B = 500)

pdf(file='nomet_gapstat_pam.pdf')
fviz_gap_stat(gap_stat)
dev.off()

pamclusters <- pam(euclidean, diss=TRUE, k = 8)

pamclusters$data <- corrmatrix

fviz_silhouette(pamclusters)

pdf(file='nomet_pam_clusters.pdf')
fviz_cluster(pamclusters,repel=TRUE,ellipse = TRUE,geom='point',main='CP,EIP,GP PAM Clusters')
dev.off()
pamclusterid <- pamclusters$clustering

pamclusterid <- tibble(name = names(pamclusterid), value = pamclusterid)

no_met_pam_json <- jsonlite::toJSON(split(pamclusterid$name, pamclusterid$value))

write(no_met_pam_json, "no_met_pam.json")

#metabolism & EIP
metabolismEIP_pathways <- c(Metabolism,Environmental_Information_Processing)

metabolismEIP_pathwayDF <- pathwayDF[,metabolismEIP_pathways]

metabolismEIP_pathwayDF <- scale(metabolismEIP_pathwayDF,center = TRUE, scale = TRUE)

bcdcor <- adjust.corr(metabolismEIP_pathwayDF,type = c('pearson'),use = c('complete.obs'),adjust = 'BH')
corrmatrix <- as.data.frame(bcdcor$R$r)
qvalues <- as.data.frame(bcdcor$R$P)

#histogram 


x <- data.frame(lower.triangle(as.matrix(corrmatrix)))

y <- data.frame(y = unlist(x))

y <- y[y$y != 0 & y$y != 1, ]

pdf('meta EIP histogram.pdf')
hist(y,main = 'Meta, EIP Correlation Histogram',xlab = 'Pearson r Value')
dev.off()



#write.csv(corrmatrix,'metabolism-EIP.csv')


squared_euclidean <- 1 - corrmatrix
euclidean <- sqrt(squared_euclidean)


#HCA

m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

#function to compute agglomerative coefficient
ac <- function(x) {
  agnes(euclidean, method = x)$ac
}

sapply(m, ac)


#wards still best - use euclidean matrix for clustering

#gap stat 

gap_stat <- clusGap(euclidean, FUN=hcut, hc_func="hclust", hc_method="ward.D", K.max = 20, B = 500)

#removed isdiss=FALSE
pdf(file='metEIP_gapstat_hca')
fviz_gap_stat(gap_stat)
dev.off()

#19 for all pathways

hc <- hclust(as.dist(euclidean), method = "ward.D2" )

# Cut tree into 4 groups

sub_grp <- cutree(hc, k = 19)

pdf('metabolismEIPHCAdend.pdf',width=15,height=10)
plot(hc, cex = 0.6,main='',xlab='',sub='')
rect.hclust(hc, k = 19, border = 2:5)
dev.off()

fviz_cluster(list(data = corrmatrix, cluster = sub_grp),repel = TRUE)

clusters <- cutree(hc,k = 19)

clusters <- tibble(name = names(clusters), value = clusters)

meta_EIP_hca_json <- jsonlite::toJSON(split(clusters$name, clusters$value))

write(meta_EIP_hca_json, "meta_EIP_hca.json")

# k medoids


gap_stat <- clusGap(euclidean, FUN=pam, K.max = 20, B = 500)
pdf(file='metEIP_gapstat_pam.pdf')
fviz_gap_stat(gap_stat)
dev.off()

pamclusters <- pam(euclidean, diss=TRUE, k = 16)

pamclusters$data <- corrmatrix

fviz_silhouette(pamclusters)
pdf('metEIP_pam_clusters.pdf')
fviz_cluster(pamclusters,repel=TRUE,ellipse = TRUE,geom='point',main = 'Metabolism, EIP PAM Clusters')
dev.off()
pamclusterid <- pamclusters$clustering

pamclusterid <- tibble(name = names(pamclusterid), value = pamclusterid)

meta_EIP_pam_json <- jsonlite::toJSON(split(pamclusterid$name, pamclusterid$value))
  
write(meta_EIP_pam_json, "meta_EIP_pam.json")
  