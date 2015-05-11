setwd("/Users/Michel/Desktop/Research/code/ProteinScavengingScreen/")
options(stringsAsFactors = FALSE)

library(data.table)
library(limma) # for voom
library(dplyr)
library(tidyr)
library(qvalue)
library(org.Mm.eg.db)

source("load_theme_for_figures.R")
theme_set(figure_theme)

source("AnalyzeScreenData/screen_analysis_functions.R") 

## INPUT FILES ##
screen1_file <- "Data/Screen read counts/screen1_results_unfiltered.txt" 
screen2_file <- "Data/Screen read counts/screen2_results_unfiltered.txt"

# filter options
min_guides_per_gene <- 3

#####
#####

# import data
s1_in <- data.table(read.table(screen1_file, header=TRUE))
s2_in <- data.table(read.table(screen2_file, header=TRUE))



# filter out genes with less than 3 guides
# also filter out extreme guides for non-normal gene-comparisons in each replicate (or not, for not.shap.filtered)
shap.filtered <- filter.data(tidy.s)
not.shap.filtered <- filter.data(tidy.s, use.shapiro=FALSE)

# ttests to generate pvalues to test differential essentiality of genes in different conditions 
# using filtered or unfiltered data
shap.ttests <- gen.pvals(shap.filtered)
unfiltered.ttests <- gen.pvals(not.shap.filtered)

save(shap.filtered, not.shap.filtered, shap.ttests, unfiltered.ttests, file="Data/Saved R data/shap_v_noShap.RData")
# save(shap.filtered, not.shap.filtered, file="Data/Saved R data/shap_v_noShap.RData")
# load(file="Data/Saved R data/shap_v_noShap.RData")

# need to run through code from here down!!! 

ttests.prespread <- shap.ttests %>% mutate(comparison = paste(var2, "v", var1, replicate, sep="."))
fold.changes <- ttests.prespread %>% ungroup %>%
  dplyr::select(-c(var1, var2, replicate, mean, n, variance, stderror, statistic, p.value)) %>% spread(comparison, median)
p.values <- ttests.prespread %>% ungroup %>%
  dplyr::select(-c(var1, var2, replicate, mean, n, variance, stderror, statistic, median)) %>% spread(comparison, p.value)

fold.changes <- na.omit(fold.changes)
p.values <- na.omit(p.values)

q.values <- p.values
for (col in colnames(p.values)[-c(1)]) {
  pvals <- p.values[[col]]
  q.values[[col]] = qvalue(pvals)$qvalues
}

q.values.filtered <- q.values %>% filter((NoLeu.v.BSA.1 < 0.25 & NoLeu.v.BSA.2 < 0.25) |
                                           (NoLeu.v.DMEM.1 < 0.25 & NoLeu.v.DMEM.2 < 0.25)) 
sig.genes <- q.values.filtered$gene_symbol

save(fold.changes, p.values, q.values, sig.genes, file="Data/Saved R data/fcs_pvals_qvals.RData")
# load(file="Data/Saved R data/fcs_pvals_qvals.RData")

## calculating fold-changes on a guide-by-guide basis (rather than gene-by-gene)
s1_gather <- s1_in %>% gather(sample, count, init:NoLeu)
s2_gather <- s2_in %>% gather(sample, count, init:NoLeu)
bothScr_gather <- rbind(cbind(s1_gather, screen=1), cbind(s2_gather, screen=2))

bothScr_sums <- bothScr_gather %>% group_by(sample, screen) %>% mutate(totalCounts = sum(count))
bothScr_cpm <- bothScr_sums %>% mutate(cpm = (count+1)/totalCounts*10^6)
bothScr_spread <- bothScr_cpm %>% dplyr::select(-count, -totalCounts) %>% spread(sample, cpm)

bothScr_compCond <- bothScr_spread %>% mutate(DMEM.v.Init = DMEM/init, BSA.v.init = BSA/init, NoLeu.v.Init = NoLeu/init,
                                              NoLeu.v.DMEM = NoLeu/DMEM, NoLeu.v.BSA = NoLeu/BSA, 
                                              NoLeu.v.bothCtrls = NoLeu/(0.5*(DMEM+BSA)))

toJoin <- bothScr_gather %>% filter(sample == "init") %>% dplyr::select(guide, count, screen)
bothScr_compCond <- bothScr_compCond %>% left_join(toJoin, by=c("guide","screen"))
bothScr_compCond <- bothScr_compCond %>% arrange(screen, gene_symbol, NoLeu.v.bothCtrls) %>% dplyr::select(-(init:NoLeu))

foldChanges_summ <- setnames(bothScr_compCond, "count", "init.count") 
save(foldChanges_summ, file="Data/Saved R data/guide_fcs.RData")
# load(file="Data/Saved R data/guide_fcs.RData")


#### library below includes annotations
converter1 = as.character(org.Mm.egSYMBOL2EG)
converter2 = as.character(org.Mm.egGENENAME)

annotations.df <- data.table(gene_symbol = unique(tidy.s$gene_symbol)) %>% 
  mutate(annotation=converter2[converter1[as.character(gene_symbol)]])

save(annotations.df, file="Data/Saved R data/annotations_df.RData")
# load(file="Data/Saved R data/annotations_df.RData")

# for gene-by-gene information: use get_gene_fcs / get_gene_pvals / get_gene_qvals / get_guide_fcs

