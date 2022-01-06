### LOADING DATA AND LIBRARIES ################################################
setwd("") #set directory
load("AOP_dataset.Rdata") #can be found in repository


library(tidyverse)
library(CONOR)
library(DescTools)
library("annotate")
library("mouse4302.db") 
library(ggplot2)
library(ggpubr)
library(limma)


### DATA CLEANING #############################################################


## getting probe ids from the studies

study1_ids = row.names(exp_values[[1]])
study2_ids = row.names(exp_values[[2]])
study3_ids = row.names(exp_values[[3]])
study4_ids = row.names(exp_values[[4]])


## getting corresponding gene names and ids for the studies
#### text files downloaded from GEO

study1_geneNames = read.delim("GPL6947-13512.txt", comment.char="#") %>% 
  dplyr::select(c("ID", "Entrez_Gene_ID"))


study2_geneNames = read.delim("GPL1322-26772.txt", comment.char="#") %>% 
  dplyr::select(c("ID", "ENTREZ_GENE_ID"))


study3_geneNames <-select(mouse4302.db, featureNames(data3[[1]]), 
                          c("SYMBOL","ENTREZID", "GENENAME"))

study4_geneNames = read.delim("GPL10333-20415.txt", comment.char="#") %>% 
  dplyr::select(c("ID", "GENE", "GENE_SYMBOL", "GENE_NAME"))



## attaching probe ids/given rownames as a column so the gene ids can be added
## filtering out any genes from GEO that are not in the study
## then filtering out any rows with NAs as the methods are unable to accommodate

study1 = as_tibble(exp_values[[1]]) %>% add_column(ID = study1_ids)
study1_geneNames = study1_geneNames %>% filter(ID %in% study1_ids)
study1_gene = study1 %>% 
  add_column(Entrez_Gene_ID = study1_geneNames$Entrez_Gene_ID) %>% drop_na()

study2 = as_tibble(exp_values[[2]]) %>% add_column(ID = study2_ids)
study2_geneNames = study2_geneNames %>% filter(ID %in% study2_ids)
study2_gene = study2 %>% 
  add_column(ENTREZ_GENE_ID = study2_geneNames$ENTREZ_GENE_ID) %>% drop_na()

study3 = as_tibble(exp_values[[3]]) %>% add_column(PROBEID = study3_ids)
study3_gene = merge(study3, study3_geneNames) %>% 
  dplyr::select(-c(PROBEID, SYMBOL, GENENAME)) %>% drop_na()

study4 = as_tibble(exp_values[[4]]) %>% add_column(ID = study4_ids)
study4_gene = merge(study4, study4_geneNames) %>% 
  dplyr::select(-c(ID, GENE_SYMBOL, GENE_NAME)) %>% drop_na()



## the methods do not work when there are multiple rows of the same gene
## therefore, the mean of each gene will be taken
## this is only on studies 3&4 as they are the only two on the same animal
## studies 1&2 will simply have the gene columns dropped
## the resulting studyi_norm will be ready for analysis

study1_norm = study1_gene %>% dplyr::select(-c(ID, Entrez_Gene_ID))
study2_norm = study2_gene %>% dplyr::select(-c(ID, ENTREZ_GENE_ID))

study3_norm = study3_gene %>% group_by(ENTREZID) %>%
  summarize_all(list(name = mean)) %>% ungroup()
study4_norm = study4_gene  %>% group_by(GENE) %>%
  summarize_all(list(name = mean)) %>% ungroup()

rowNames3 = study3_norm$ENTREZID
rowNames4 = study4_norm$GENE

study3_norm = study3_norm %>% dplyr::select(-c(ENTREZID))
study4_norm = study4_norm %>% dplyr::select(-c(GENE))
row.names(study3_norm) <- rowNames3
row.names(study4_norm) <- rowNames4



## creating two tables combining studies 3&4 and 1&2, respectively,
## pre-normalization to allow for graphical comparison

### studies 3&4 - both on mice

temp_study4 = study4_norm %>% transmute(GSE60541 = rowMeans(study4_norm)) %>% 
  add_column(gene = rowNames4, 
             duplicate = duplicated(c(rowNames3, rowNames4))[22029:44028]) %>%
  filter(duplicate == TRUE) %>% dplyr::select(-c(duplicate))

temp_study3 = study3_norm %>% transmute(GSE85359 = rowMeans(study3_norm)) %>% 
  add_column(gene = rowNames3) %>% filter(gene %in% temp_study4$gene)

mergedMouse = join(temp_study3, temp_study4, by = "gene") %>% 
  dplyr::select(-c(gene))

gene = temp_study3$gene # to keep track of mutual genes

### studies 1&2, on humans and fruit flies respectively

study1_temp = study1_norm %>% transmute(GSE81067 = rowMeans(study1_norm)) %>%
  add_column(ones = rep(1, dim(study1_norm)[1]))
study2_temp = study2_norm %>% transmute(GSE37404 = rowMeans(study2_norm)) %>% 
  add_column(ones = rep(1, dim(study2_norm)[1]))

study1_temp = study1_temp %>% mutate(counter = c(1:35280)) %>% 
  filter(counter <= 18952)

mergedHumanFlies = study2_temp %>% add_column(GSE81067 = study1_temp$GSE81067)



### CROSS-PLATFORM NORMALIZATION ##############################################


set.seed(29)

## studies 3&4

xpn_34 = xpn(platform1.data=study3_norm, platform2.data=study4_norm, 
             p1.names = 0, p2.names = 0)
dwd_34 = dwd(platform1.data = study3_norm, platform2.data = study4_norm, 
             p1.names = 0, p2.names = 0)
eb_34 = eb(platform1.data = study3_norm, platform2.data = study4_norm, 
           p1.names = 0, p2.names = 0)
gq_34 = gq(platform1.data = study3_norm, platform2.data =  study4_norm, 
           p1.names = 0, p2.names = 0)

## studies 1&2
xpn_12 = xpn(platform1.data=study1_norm, platform2.data=study2_norm)
dwd_12 = dwd(platform1.data = study1_norm, platform2.data = study2_norm)
eb_12 = eb(platform1.data = study1_norm, platform2.data = study2_norm)
gq_12 = gq(platform1.data = study1_norm, platform2.data =  study2_norm)




### DATA INTERPRETATION/PLOTIING ##############################################


## mean/mean plots will be used with a x=y line for reference



## creating tibbles with the mean value of each gene across the arrays for each
## successful method

## studies 3&4
xpn_mean34 = xpn_34$x %>% mutate(GSE85359 = as.double(rowMeans(xpn_34$x)), 
                                 GSE60541 = as.double(rowMeans(xpn_34$y))) %>% 
  dplyr::select(c(GSE85359, GSE60541))

dwd_mean34 = dwd_34$x %>% mutate(GSE85359 = as.double(rowMeans(dwd_34$x)), 
                                 GSE60541 = as.double(rowMeans(dwd_34$y))) %>% 
  dplyr::select(c(GSE85359, GSE60541))
gq_mean34 = gq_34$x %>% mutate(GSE85359 = as.double(rowMeans(gq_34$x)), 
                               GSE60541 = as.double(rowMeans(gq_34$y))) %>% 
  dplyr::select(c(GSE85359, GSE60541))
## studies 1&2
xpn_mean12 = xpn_12$x %>% mutate(GSE81067 = as.double(rowMeans(xpn_12$x)), 
                                 GSE37404 = as.double(rowMeans(xpn_12$y))) %>% 
  dplyr::select(c(GSE81067, GSE37404))
dwd_mean12 = dwd_12$x %>% mutate(GSE81067 = as.double(rowMeans(dwd_12$x)), 
                                 GSE37404 = as.double(rowMeans(dwd_12$y))) %>% 
  dplyr::select(c(GSE81067, GSE37404))
gq_mean12 = gq_12$x %>% mutate(GSE81067 = as.double(rowMeans(gq_12$x)), 
                               GSE37404 = as.double(rowMeans(gq_12$y))) %>% 
  dplyr::select(c(GSE81067, GSE37404))


## plotting mean-mean plots, including non-normalized data

## studies 3&4
xpn_mean34_plot = ggplot(data = xpn_mean34, aes(x = GSE85359, y = GSE60541)) +
  geom_point() +
  ggtitle("XPN") +
  xlab("Study 3 (GSE85359)")+
  ylab("Study 4 (GSE60541)")+
  geom_abline(slope=1, intercept = 0, colour = "red")

dwd_mean34_plot = ggplot(data = dwd_mean34, aes(x = GSE85359, y = GSE60541)) +
  geom_point()+
  ggtitle("DWD") +
  xlab("Study 3 (GSE85359)")+
  ylab("Study 4 (GSE60541)")+
  geom_abline(slope=1, intercept = 0, colour = "red")

gq_mean34_plot = ggplot(data = gq_mean34, aes(x = GSE85359, y = GSE60541)) +
  geom_point()+
  ggtitle("GQ") +
  xlab("Study 3 (GSE85359)")+
  ylab("Study 4 (GSE60541)")+
  geom_abline(slope=1, intercept = 0, colour = "red")

nonorm_mean34_plot = ggplot(data = mergedMouse, 
                            aes(x = GSE85359, y = GSE60541)) +
  geom_point()+
  ggtitle("No Normalization") +
  xlab("Study 3 (GSE85359)")+
  ylab("Study 4 (GSE60541)")+
  geom_abline(slope=1, intercept = 0, colour = "red")

## plotting together
figure34 <- ggarrange(nonorm_mean34_plot, xpn_mean34_plot, dwd_mean34_plot, 
                      gq_mean34_plot, ncol = 2, nrow = 2)
figure34 <- annotate_figure(figure34, top = "Comparison of Post-Normlization Log-Transformed Gene Expression Values")
figure34


## studies 1&2
xpn_mean12_plot = ggplot(data = xpn_mean12, aes(x = GSE81067, y = GSE37404)) +
  geom_point() +
  ggtitle("XPN") +
  xlab("Study 1 (GSE81067)")+
  ylab("Study 2 (GSE37404)")+
  geom_abline(slope=1, intercept = 0, colour = "red")

dwd_mean12_plot = ggplot(data = dwd_mean12, aes(x = GSE81067, y = GSE37404)) +
  geom_point()+
  ggtitle("DWD") +
  xlab("Study 1 (GSE81067)")+
  ylab("Study 2 (GSE37404)")+
  geom_abline(slope=1, intercept = 0, colour = "red")

gq_mean12_plot = ggplot(data = gq_mean12, aes(x = GSE81067, y = GSE37404)) +
  geom_point()+
  ggtitle("GQ") +
  xlab("Study 1 (GSE81067)")+
  ylab("Study 2 (GSE37404)")+
  geom_abline(slope=1, intercept = 0, colour = "red")

nonorm_mean12_plot = ggplot(data = mergedHumanFlies, 
                            aes(x = GSE81067, y = GSE37404)) +
  geom_point()+
  ggtitle("No Normalization") +
  xlab("Study 1 (GSE81067)")+
  ylab("Study 2 (GSE37404)")+
  geom_abline(slope=1, intercept = 0, colour = "red")

## plotting together
figure12 <- ggarrange(nonorm_mean12_plot, xpn_mean12_plot, dwd_mean12_plot, 
                      gq_mean12_plot, ncol = 2, nrow = 2)
figure12 <- annotate_figure(figure12, top = "Comparison of Post-Normlization Log-Transformed Gene Expression Values")
figure12



## correlation/concordance results

## studies 3&4
## concordance correlation
xpn_ccc34 = CCC(xpn_mean34$GSE85359, xpn_mean34$GSE60541, 
                ci = "z-transform", conf.level = 0.95, na.rm = FALSE)
dwd_ccc34 = CCC(dwd_mean34$GSE85359, dwd_mean34$GSE60541, 
                ci = "z-transform", conf.level = 0.95, na.rm = FALSE)
gq_ccc34 = CCC(gq_mean34$GSE85359, gq_mean34$GSE60541, 
               ci = "z-transform", conf.level = 0.95, na.rm = FALSE)

xpn_ccc34$rho.c
dwd_ccc34$rho.c
gq_ccc34$rho.c


## studies 1&2
## concordance correlation
xpn_ccc12 = CCC(xpn_mean12$GSE81067, xpn_mean12$GSE37404, 
                ci = "z-transform", conf.level = 0.95, na.rm = FALSE)
dwd_ccc12 = CCC(dwd_mean12$GSE81067, dwd_mean12$GSE37404, 
                ci = "z-transform", conf.level = 0.95, na.rm = FALSE)
gq_ccc12 = CCC(gq_mean12$GSE81067, gq_mean12$GSE37404, 
               ci = "z-transform", conf.level = 0.95, na.rm = FALSE)

xpn_ccc12$rho.c
dwd_ccc12$rho.c
gq_ccc12$rho.c



### MERGING STUDIES 1,3,4
xpn_GSE85359 = xpn_34$x %>% add_column(gene)
xpn_GSE60541 = xpn_34$y %>% add_column(gene)
xpn_mouseData = join(xpn_GSE85359, xpn_GSE60541) %>% dplyr::select(-gene)

dwd_GSE85359 = dwd_34$x %>% add_column(gene)
dwd_GSE60541 = dwd_34$y %>% add_column(gene)
dwd_mouseData = join(dwd_GSE85359, dwd_GSE60541) %>% dplyr::select(-gene)

gq_GSE85359 = gq_34$x %>% add_column(gene)
gq_GSE60541 = gq_34$y %>% add_column(gene)
gq_mouseData = join(gq_GSE85359, gq_GSE60541) %>% dplyr::select(-gene)

## Cross platform normalization
xpn_134 = xpn(xpn_mouseData, study1_norm)
dwd_134 = dwd(dwd_mouseData, study1_norm)
gq_134 = gq(gq_mouseData, study1_norm)


## creating mean-mean tibbles

xpn_mean134 = xpn_134$x %>% mutate(mouse = as.double(rowMeans(xpn_134$x)), 
                                   human = as.double(rowMeans(xpn_134$y))) %>% 
  dplyr::select(c(mouse, human))
dwd_mean134 = dwd_134$x %>% mutate(mouse = as.double(rowMeans(dwd_134$x)), 
                                   human = as.double(rowMeans(dwd_134$y))) %>% 
  dplyr::select(c(mouse, human))
gq_mean134 = gq_134$x %>% mutate(mouse = as.double(rowMeans(gq_134$x)), 
                                 human = as.double(rowMeans(gq_134$y))) %>% 
  dplyr::select(c(mouse, human))

## plotting

xpn_mean134_plot = ggplot(data = xpn_mean134, 
                          aes(x = mouse, y = human)) +
  geom_point() +
  ggtitle("XPN") +
  xlab("Mice")+
  ylab("Human")+
  geom_abline(slope=1, intercept = 0, colour = "red")

dwd_mean134_plot = ggplot(data = dwd_mean134, 
                          aes(x = mouse, y = human)) +
  geom_point() +
  ggtitle("DWD") +
  xlab("Mice")+
  ylab("Human")+
  geom_abline(slope=1, intercept = 0, colour = "red")

gq_mean134_plot = ggplot(data = gq_mean134, 
                         aes(x = mouse, y = human)) +
  geom_point() +
  ggtitle("GQ") +
  xlab("Mice")+
  ylab("Human")+
  geom_abline(slope=1, intercept = 0, colour = "red")

figure134 <- ggarrange(xpn_mean134_plot, dwd_mean134_plot, gq_mean134_plot,
                       ncol = 2, nrow = 2)
figure134 <- annotate_figure(figure134, top = "Comparison of Post-Normlization Log-Transformed Gene Expression Values")
figure134


## concordance
xpn_ccc134 = CCC(xpn_mean134$mouse, xpn_mean134$human, 
                 ci = "z-transform", conf.level = 0.95, na.rm = FALSE)

dwd_ccc134 = CCC(dwd_mean134$mouse, dwd_mean134$human, 
                 ci = "z-transform", conf.level = 0.95, na.rm = FALSE)

gq_ccc134 = CCC(gq_mean134$mouse, gq_mean134$human, 
                ci = "z-transform", conf.level = 0.95, na.rm = FALSE)

xpn_ccc134$rho.c
dwd_ccc134$rho.c
gq_ccc134$rho.c




### MERGING ALL FOUR DATASETS

## creating two datsets from the four normalized ones for each method

xpn_human = xpn_12$x %>% add_column(counter = c(1:18952))
xpn_fruitfly = xpn_12$y %>% add_column(counter = c(1:18952))
xpn_humanFlyData = join(xpn_human, xpn_fruitfly) %>% dplyr::select(-counter)

dwd_human = dwd_12$x %>% add_column(counter = c(1:18952))
dwd_fruitfly = dwd_12$y %>% add_column(counter = c(1:18952))
dwd_humanFlyData = join(dwd_human, dwd_fruitfly) %>% dplyr::select(-counter)

gq_human = gq_12$x %>% add_column(counter = c(1:18952))
gq_fruitfly = gq_12$y %>% add_column(counter = c(1:18952))
gq_humanFlyData = join(gq_human, gq_fruitfly) %>% dplyr::select(-counter)


## Cross platform normalization
xpn_1234 = xpn(xpn_mouseData, xpn_humanFlyData)
dwd_1234 = dwd(dwd_mouseData, dwd_humanFlyData)
gq_1234 = gq(gq_mouseData, gq_humanFlyData)


## creating mean-mean tibbles

xpn_mean1234 = xpn_1234$x %>% mutate(mouse = as.double(rowMeans(xpn_1234$x)), 
                                     human_fruitfly = as.double(rowMeans(xpn_1234$y))) %>% 
  dplyr::select(c(mouse, human_fruitfly))
dwd_mean1234 = dwd_1234$x %>% mutate(mouse = as.double(rowMeans(dwd_1234$x)), 
                                     human_fruitfly = as.double(rowMeans(dwd_1234$y))) %>% 
  dplyr::select(c(mouse, human_fruitfly))
gq_mean1234 = gq_1234$x %>% mutate(mouse = as.double(rowMeans(gq_1234$x)), 
                                   human_fruitfly = as.double(rowMeans(gq_1234$y))) %>% 
  dplyr::select(c(mouse, human_fruitfly))

## plotting
xpn_mean1234_plot = ggplot(data = xpn_mean1234, 
                           aes(x = mouse, y = human_fruitfly)) +
  geom_point() +
  ggtitle("XPN") +
  xlab("Mice")+
  ylab("Human/Fruitfly")+
  geom_abline(slope=1, intercept = 0, colour = "red")

dwd_mean1234_plot = ggplot(data = dwd_mean1234, 
                           aes(x = mouse, y = human_fruitfly)) +
  geom_point() +
  ggtitle("DWD") +
  xlab("Mice")+
  ylab("Human/Fruitfly")+
  geom_abline(slope=1, intercept = 0, colour = "red")

gq_mean1234_plot = ggplot(data = gq_mean1234, 
                          aes(x = mouse, y = human_fruitfly)) +
  geom_point() +
  ggtitle("GQ") +
  xlab("Mice")+
  ylab("Human/Fruitfly")+
  geom_abline(slope=1, intercept = 0, colour = "red")

figure1234 <- ggarrange(xpn_mean1234_plot, dwd_mean1234_plot, gq_mean1234_plot,
                        ncol = 2, nrow = 2)
figure1234 <- annotate_figure(figure1234, top = "Comparison of Post-Normlization Log-Transformed Gene Expression Values")
figure1234


## concordance
xpn_ccc1234 = CCC(xpn_mean1234$mouse, xpn_mean1234$human_fruitfly, 
                  ci = "z-transform", conf.level = 0.95, na.rm = FALSE)

dwd_ccc1234 = CCC(dwd_mean1234$mouse, dwd_mean1234$human_fruitfly, 
                  ci = "z-transform", conf.level = 0.95, na.rm = FALSE)

gq_ccc1234 = CCC(gq_mean1234$mouse, gq_mean1234$human_fruitfly, 
                 ci = "z-transform", conf.level = 0.95, na.rm = FALSE)

xpn_ccc1234$rho.c
dwd_ccc1234$rho.c
gq_ccc1234$rho.c




### DIFFERENTIAL EXPRESSSION ANALYSIS #########################################

## using mice data and then connecting it to noted human genes via orthologs 
## see appendix b for relevant orthologs

## loading in genes of interest
ortho_genes = read.csv("gene conversion.csv")
human_to_mouse = ortho_genes %>% filter(ortho_genes$ENTREZID %in% gene)
missing = ortho_genes %>% filter(!unigene %in% human_to_mouse$unigene)

## study 3 - exposure to 15 Gy

study3_design15 <- pheno_data[[3]][c(10:12,16:18)]
study3_expression15 = as.matrix(xpn_34$x[,c(10:12,16:18)])

study3_design15 <- model.matrix(~ 0+factor(c(1,1,1,2,2,2)))
colnames(study3_design15) <- c("Control","Radiation")
study3_forDE15 <- ExpressionSet(assayData=study3_expression15)
study3_fit15 <- lmFit(study3_forDE15, study3_design15)
study3_contr15 <- makeContrasts(Radiation - Control,
                                levels = colnames(coef(study3_fit15)))
study3_tmp15 <- contrasts.fit(study3_fit15, study3_contr15)
study3_fit15 <- eBayes(study3_tmp15, trend=TRUE, robust=TRUE)
study3_results15 <- decideTests(study3_fit15)
summary(study3_results15)

## filtering to only contain ortholog genes and adding missing required  
## ortholog genes with values of NA
study3_FC15 = as_tibble(study3_fit15$coefficients)
study3_FC15 = as_tibble(study3_FC15) %>% add_column(ENTREZID = gene) %>% 
  filter(ENTREZID %in% human_to_mouse$ENTREZID)
study3_FC15 = join(study3_FC15, human_to_mouse)[,c(1,3)]
colnames(study3_FC15) = c("FC", "unigene")
study3_FC15 = study3_FC15 %>% add_row(unigene = "TRIM22", FC = NA) %>%
  add_row(unigene = "RMI2", FC = NA) %>%
  add_row(unigene = "H4C3", FC = NA)


## study 3 exposure to 12.5 Gy

study3_design12 <- pheno_data[[3]][c(10:15)]
study3_expression12 = as.matrix(xpn_34$x[,c(10:15)])

study3_design12 <- model.matrix(~ 0+factor(c(1,1,1,2,2,2)))
colnames(study3_design12) <- c("Control","Radiation")
study3_forDE12 <- ExpressionSet(assayData=study3_expression12)
study3_fit12 <- lmFit(study3_forDE12, study3_design12)
study3_contr12 <- makeContrasts(Radiation - Control, 
                                levels = colnames(coef(study3_fit12)))
study3_tmp12 <- contrasts.fit(study3_fit12, study3_contr12)
study3_fit12 <- eBayes(study3_tmp12, trend=TRUE, robust=TRUE)
study3_results12 <- decideTests(study3_fit12)
summary(study3_results12)

## filtering out non ortholog genes and adding missing required ortholog genes 
## with values of NA
study3_FC12 = as_tibble(study3_fit12$coefficients)
study3_FC12 = as_tibble(study3_FC12) %>% add_column(ENTREZID = gene) %>% 
  filter(ENTREZID %in% human_to_mouse$ENTREZID)
study3_FC12 = join(study3_FC12, human_to_mouse)[, c(1,3)]
colnames(study3_FC12) = c("FC", "unigene")
study3_FC12 = study3_FC12 %>% add_row(unigene = "TRIM22", FC = NA) %>%
  add_row(unigene = "RMI2", FC = NA) %>%
  add_row(unigene = "H4C3", FC = NA)




## study 4, only expsoure level is 90

study4_design <- pheno_data[[4]]
study4_expression = as.matrix(xpn_34$y)

study4_design <- model.matrix(~ 0+factor(c(1,1,1,1,2,2,2,2,2)))
colnames(study4_design) <- c("Control","Radiation")
study4_forDE <- ExpressionSet(assayData=(study4_expression))
study4_fit <- lmFit(study4_forDE, study4_design)
study4_contr <- makeContrasts(Radiation - Control, 
                              levels = colnames(coef(study4_fit)))
study4_tmp <- contrasts.fit(study4_fit, study4_contr)
study4_fit <- eBayes(study4_tmp, trend=TRUE, robust=TRUE)
study4_results <- decideTests(study4_fit)
summary(study4_results)

## filtering non ortholog genes and adding missing required ortholog genes with 
## values of NA
study4_FC = as_tibble(study4_fit$coefficients)
study4_FC = as_tibble(study4_FC) %>% add_column(ENTREZID = gene) %>% 
  filter(ENTREZID %in% human_to_mouse$ENTREZID)
study4_FC = join(study4_FC, human_to_mouse)[, c(1,3)]
colnames(study4_FC) = c("FC", "unigene")
study4_FC = study4_FC %>% add_row(unigene = "TRIM22", FC = NA) %>%
  add_row(unigene = "RMI2", FC = NA) %>%
  add_row(unigene = "H4C3", FC = NA)



