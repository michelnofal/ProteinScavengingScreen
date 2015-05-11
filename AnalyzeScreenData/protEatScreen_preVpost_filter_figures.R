# loading shap.filtered and not.shap.filtered
load(file="protein eating screen/protEatScreen_analysis_Rdata/shap_v_noShap.RData")

#####
## making plots of guide fold-changes with and without filtering
diff.no.filter <- ggplot(not.shap.filtered %>% ungroup %>% filter(var1 != "init"), aes(difference)) + geom_histogram() + 
  facet_grid(replicate ~ var1 + var2, scale="free_y")
diff.shap.filter <- ggplot(shap.filtered %>% ungroup %>% filter(var1 != "init"), aes(difference)) + geom_histogram() + 
  facet_grid(replicate ~ var1 + var2, scale="free_y")

ggsave(
  "protein eating screen/figures/screen_summary/foldChanges.noFilter.pdf",
  diff.no.filter,
  width=18,height=18,units="cm",
  dpi=300
)
ggsave(
  "protein eating screen/figures/screen_summary/foldChanges.shapiroFilter.pdf",
  diff.shap.filter,
  width=18,height=18,units="cm",
  dpi=300
)

#####
## making plots of p-values for normality for filtered and unfiltered data
shap.toshow <- shap.filtered %>% summarize(shap.prefilter=mean(shapiro), shap.postfilter=mean(new.shapiro)) %>%
  ungroup %>% filter(var1 != "init" & var2=="NoLeu" & replicate=="2")

normality.noFilter <- ggplot(shap.toshow, aes(shap.prefilter)) + geom_histogram() + 
  facet_grid( ~ var1 + var2, scale="free_y") +
  xlab("\nShapiro normality test p.values without filtering") + ylab("") + scale_y_continuous(breaks=c(), labels=c())
normality.shapiroFilter <- ggplot(shap.toshow, aes(shap.postfilter)) + geom_histogram() + 
  facet_grid( ~ var1 + var2, scale="free_y") +
  xlab("\nShapiro normality test p.values after filtering") + ylab("") + scale_y_continuous(breaks=c(), labels=c())

ggsave(
  "protein eating screen/figures/screen_summary/normality.noFilter.pdf",
  normality.noFilter,
  width=12,height=7,units="cm",
  dpi=300
)
ggsave(
  "protein eating screen/figures/screen_summary/normality.shapiroFilter.pdf",
  normality.shapiroFilter,
  width=12,height=7,units="cm",
  dpi=300
)

#####
# loading shap.ttests and unfiltered.ttests
load(file="protein eating screen/protEatScreen_analysis_Rdata/shap_v_noShap_ttests.RData")

shap.pval.toPlot <- shap.ttests %>% ungroup %>% filter(var1 != "init" & replicate=="2")
unfiltered.pval.toPlot <- unfiltered.ttests %>% ungroup %>% filter(var1 != "init" & replicate=="2")

shap.pval.plot <- ggplot(shap.pval.toPlot, aes(p.value)) + geom_histogram(binwidth=0.02) + facet_grid( ~ var1 + var2, scale="free_y")
unfiltered.pval.plot <- ggplot(unfiltered.pval.toPlot, aes(p.value)) + geom_histogram(binwidth=0.02) + facet_grid( ~ var1 + var2, scale="free_y")

ggsave(
  "protein eating screen/figures/screen_summary/shapiro.pval.plot.pdf",
  shap.pval.plot,
  width=16,height=7,units="cm",
  dpi=300
)
ggsave(
  "protein eating screen/figures/screen_summary/nofilt.pval.plot.pdf",
  unfiltered.pval.plot,
  width=16,height=7,units="cm",
  dpi=300
)

#####
# loading fold.changes , p.values , q.values
load(file="protein eating screen/protEatScreen_analysis_Rdata/fcs_pvals_qvals.RData")

qvals.melt <- melt(q.values, id.vars=c("gene_symbol"), measure.vars=c("DMEM.NoLeu.1","BSA.NoLeu.1","DMEM.NoLeu.2","BSA.NoLeu.2"))
qvals.plot <- ggplot(qvals.melt, aes(value)) + geom_histogram(binwidth=0.02) + facet_wrap( ~ variable, nrow=2)

ggsave(
  "protein eating screen/figures/screen_summary/qvals.plot.pdf",
   qvals.plot,
   width=12,height=12,units="cm",
   dpi=300
)


