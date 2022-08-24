
library(rjson)
pathwayenrichment <- rjson::fromJSON(file = 'allpathways_pathwayenrichment')
pathwaydistribution <- rjson::fromJSON(file = 'allpathways_pathwaydistribution')

pathwayenrichment <- as.data.frame(pathwayenrichment)
rownames(pathwayenrichment) <-  c('Met','CP','EIP','GIP')
colnames(pathwayenrichment) <- c('1 - 27','2 - 5','3 - 8','4 - 6','5 - 4','6 - 8','7 - 8','8 - 31' ,'9 - 6')
pathwaydistribution <- as.data.frame(pathwaydistribution)


pdf('all pathways network enrichment.pdf', useDingbats = FALSE,height = 10)
par(mar = c(4.1,3.1,3.1,1.1))
dotchart(as.matrix(pathwayenrichment),color=ifelse(t(pathwaydistribution) < 0.05,'red','black'),main='Metabolism, CP, EIP, GIP Network Cluster Enrichment Profiles',xlab = 'Enrichment Score',ylab = 'Cluster & Length of Cluster')
dev.off()

############# remove met 

pathwayenrichment <- rjson::fromJSON(file = 'non-metabolism_pathwayenrichment')
pathwaydistribution <- rjson::fromJSON(file = 'non-metabolism_pathwaydistribution')

pathwayenrichment <- as.data.frame(pathwayenrichment)
rownames(pathwayenrichment) <-  c('CP','EIP','GIP')
colnames(pathwayenrichment) <- c('1 - 8','2 - 8','3 - 14','4 - 6','5 - 13','6 - 4')
pathwaydistribution <- as.data.frame(pathwaydistribution)

pdf('non-metabolism network enrichment.pdf', useDingbats = FALSE,height = 10)
par(mar = c(4.1,3.1,3.1,1.1))
dotchart(as.matrix(pathwayenrichment),color=ifelse(t(pathwaydistribution) < 0.05,'red','black'),main='CP, EIP, GIP Network Cluster Enrichment Profiles',xlab = 'Enrichment Score',ylab = 'Cluster & Length of Cluster')
dev.off()

#met EIP

pathwayenrichment <- rjson::fromJSON(file = 'metabolismEIP_pathwayenrichment')
pathwaydistribution <- rjson::fromJSON(file = 'metabolismEIP_pathwaydistribution')

pathwayenrichment <- as.data.frame(pathwayenrichment)
rownames(pathwayenrichment) <-  c('Met','EIP')
colnames(pathwayenrichment) <- c('1 - 32','2 - 5','3 - 10','4 - 12','5 - 28')
pathwaydistribution <- as.data.frame(pathwaydistribution)


pdf('Metabolism EIP network enrichment.pdf', useDingbats = FALSE,height = 10)
par(mar = c(4.1,3.1,3.1,1.1))
dotchart(as.matrix(pathwayenrichment),color=ifelse(t(pathwaydistribution) < 0.05,'red','black'),main='Metabolism and EIP Network Cluster Enrichment Profiles',xlab = 'Enrichment Score',ylab = 'Cluster & Length of Cluster')
dev.off()


##########CLUSTER ENRICHMENT 


#pam met EIP

pathwayenrichment <- rjson::fromJSON(file = 'meta_EIP_pam_pathwayenrichment')

pathwaydistribution <- rjson::fromJSON(file = 'meta_EIP_pam_pathwaydistribution')

pathwayenrichment <- as.data.frame(pathwayenrichment)
rownames(pathwayenrichment) <-  c('Met','EIP')
colnames(pathwayenrichment) <- c('1 - 8','2 - 8','3 - 7','4 - 3','5 - 16','6 - 7', '7 - 3', '8 - 7', '9 - 5','10 - 4','11 - 3', '12 - 4','13 - 3','14 - 4', '15 - 21', '16 - 2')
pathwaydistribution <- as.data.frame(pathwaydistribution)


pdf('Metabolism EIP pam enrichment.pdf', useDingbats = FALSE,height = 10)
par(mar = c(4.1,3.1,3.1,1.1))
dotchart(as.matrix(pathwayenrichment),color=ifelse(t(pathwaydistribution) < 0.05,'red','black'),main='Metabolism and EIP PAM Cluster Enrichment Profiles',xlab = 'Enrichment Score',ylab = 'Cluster & Length of Cluster')
dev.off()

#hca met EIP


pathwayenrichment <- rjson::fromJSON(file = 'meta_EIP_hca_pathwayenrichment')

pathwaydistribution <- rjson::fromJSON(file = 'meta_EIP_hca_pathwaydistribution')

pathwayenrichment <- as.data.frame(pathwayenrichment)
rownames(pathwayenrichment) <-  c('Met','EIP')
colnames(pathwayenrichment) <- c('1 - 5','2 - 5','3 - 5','4 - 3','5 - 12','6 - 2', '7 - 6', '8 - 8', '9 - 7','10 - 3','11 - 6', '12 - 3','13 - 8','14 - 3', '15 - 2', '16 - 15','17 - 4', '18 -2', '19 - 6')
pathwaydistribution <- as.data.frame(pathwaydistribution)


pdf('Metabolism EIP hca enrichment.pdf', useDingbats = FALSE,height = 10)
par(mar = c(4.1,3.1,3.1,1.1))
dotchart(as.matrix(pathwayenrichment),color=ifelse(t(pathwaydistribution) < 0.05,'red','black'),main='Metabolism and EIP HCA Cluster Enrichment Profiles',xlab = 'Enrichment Score',ylab = 'Cluster & Length of Cluster')
dev.off()

#pam all

pathwayenrichment <- rjson::fromJSON(file = 'all_pam_pathwayenrichment')

pathwaydistribution <- rjson::fromJSON(file = 'all_pam_pathwaydistribution')

pathwayenrichment <- as.data.frame(pathwayenrichment)
rownames(pathwayenrichment) <-  c('Met','CP','EIP','GP')
colnames(pathwayenrichment) <- c('1 - 8','2 - 21','3 - 7','4 - 3','5 - 10','6 - 3', '7 - 8', '8 - 8', '9 - 4','10 - 18','11 - 4', '12 - 8','13 - 5','14 - 8', '15 - 7', '16 - 18', '17 - 4')
pathwaydistribution <- as.data.frame(pathwaydistribution)


pdf('all pam enrichment.pdf', useDingbats = FALSE,height = 15)
par(mar = c(4.1,3.1,3.1,1.1))
dotchart(as.matrix(pathwayenrichment),color=ifelse(t(pathwaydistribution) < 0.05,'red','black'),main='All Pathway Cluster PAM Enrichment Profiles',xlab = 'Enrichment Score',ylab = 'Cluster & Length of Cluster')
dev.off()

#hca all


pathwayenrichment <- rjson::fromJSON(file = 'all_hca_pathwayenrichment')

pathwaydistribution <- rjson::fromJSON(file = 'all_hca_pathwaydistribution')

pathwayenrichment <- as.data.frame(pathwayenrichment)
rownames(pathwayenrichment) <-  c('Met','CP','EIP','GP')
colnames(pathwayenrichment) <- c('1 - 8','2 - 15','3 - 5','4 - 11','5 - 4','6 - 12', '7 - 12', '8 - 3', '9 - 8','10 - 3','11 - 4', '12 - 8','13 - 7','14 - 6', '15 - 4', '16 - 16','17 - 9', '18 - 4', '19 - 5')
pathwaydistribution <- as.data.frame(pathwaydistribution)


pdf('all hca enrichment.pdf', useDingbats = FALSE,height = 15)
par(mar = c(4.1,3.1,3.1,1.1))
dotchart(as.matrix(pathwayenrichment),color=ifelse(t(pathwaydistribution) < 0.05,'red','black'),main='All Pathway Cluster HCA Enrichment Profiles',xlab = 'Enrichment Score',ylab = 'Cluster & Length of Cluster')
dev.off()

#pam no met

pathwayenrichment <- rjson::fromJSON(file = 'no_met_pam_pathwayenrichment')

pathwaydistribution <- rjson::fromJSON(file = 'no_met_pam_pathwaydistribution')

pathwayenrichment <- as.data.frame(pathwayenrichment)
rownames(pathwayenrichment) <-  c('CP','EIP','GP')
colnames(pathwayenrichment) <- c('1 - 5','2 - 10','3 - 8','4 - 3','5 - 21','6 - 6', '7 - 9', '8 - 8')
pathwaydistribution <- as.data.frame(pathwaydistribution)


pdf('no met pam enrichment.pdf', useDingbats = FALSE,height = 10)
par(mar = c(4.1,3.1,3.1,1.1))
dotchart(as.matrix(pathwayenrichment),color=ifelse(t(pathwaydistribution) < 0.05,'red','black'),main='All Pathway Cluster PAM Enrichment Profiles',xlab = 'Enrichment Score',ylab = 'Cluster & Length of Cluster')
dev.off()

#hca no met


pathwayenrichment <- rjson::fromJSON(file = 'no_met_hca_pathwayenrichment')

pathwaydistribution <- rjson::fromJSON(file = 'no_met_hca_pathwaydistribution')

pathwayenrichment <- as.data.frame(pathwayenrichment)
rownames(pathwayenrichment) <-  c('CP','EIP','GP')
colnames(pathwayenrichment) <- c('1 - 2','2 - 3','3 - 4','4 - 2','5 - 3','6 - 2', '7 - 3', '8 - 4', '9 - 3','10 - 5','11 - 11', '12 - 6','13 - 5','14 - 3', '15 - 4', '16 - 3','17 - 5', '18 - 2')
pathwaydistribution <- as.data.frame(pathwaydistribution)


pdf('no met hca enrichment.pdf', useDingbats = FALSE,height = 12)
par(mar = c(4.1,3.1,3.1,1.1))
dotchart(as.matrix(pathwayenrichment),color=ifelse(t(pathwaydistribution) < 0.05,'red','black'),main='All Pathway Cluster HCA Enrichment Profiles',xlab = 'Enrichment Score',ylab = 'Cluster & Length of Cluster')
dev.off()