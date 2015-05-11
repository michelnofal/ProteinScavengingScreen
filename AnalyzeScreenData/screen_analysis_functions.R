###

# convert input matrices to pairwise comparisons of counts between conditions
convert.pairings = function(id.cols, m) { #}, replication) {    
  # get the left and right sides of the comparison
  dat = as.data.frame(cbind(id.cols, m))
  
  ids = c("guide", "gene_symbol")
  left = cbind(id.cols, m[, pairings[, 1]])
  right = cbind(id.cols, m[, pairings[, 2]])
  leftm = melt(left, id=ids)
  rightm = melt(right, id=ids)
  # now has duplicated factor levels
  leftm$variable = factor(as.character(leftm$variable))
  rightm$variable = factor(as.character(rightm$variable))
  
  # combine
  setnames(leftm, c("variable", "value"), c("var1", "val1"))
  leftm$var2 = rightm$variable
  leftm$val2 = rightm$value
  
  leftm
}

make.tidy = function(dat) {
  # create the matrices to be rearranged
  counts = as.matrix(dat[, c(-1, -2), with=FALSE])
  v = voom(counts, normalize.method = "scale")
  normed = v$E
  w = v$weights
  colnames(w) = colnames(counts)
  
  # rearrange as pairs
  idcols = dat[, 1:2, with=FALSE]
  ret = convert.pairings(idcols, counts)
  normpairs = convert.pairings(idcols, normed)
  
  # not going to use these weights
  # wpairs = convert.pairings(idcols, w)
  
  # ret$w1 = wpairs$val1
  # ret$w2 = wpairs$val2
  ret$norm1 = normpairs$val1
  ret$norm2 = normpairs$val2
  ret
}

filter.data = function(tidy.dat, use.shapiro = TRUE) {   
  # filter out 1) genes with < 3 guides and 2) rows in which counts are 0 for both condition
  filtered.1 <- tidy.dat %>% group_by(gene_symbol, var1, var2, replicate) %>% mutate(difference=norm2-norm1, n=length(guide))
  filtered.1 <- filtered.1 %>% filter((n >= min_guides_per_gene) & !(val1 == 0 & val2 == 0))
  
  # use shapiro test for normality over guide fold-changes for each gene-comparison-replicate group
  # (even if not filtering to compare shapiro p-values with and without filtering)
  filtered.1 <- filtered.1 %>% mutate(shapiro=tryCatch(shapiro.test(difference)$p.value, error=function(x) as.numeric(NA)), 
                                      not.normal=shapiro<0.05)
  
  if(use.shapiro) {
    # filter extreme guides for non-normally distributed gene-comparison-replicate groups
    filtered.1 <- filtered.1 %>% mutate(extreme.guide=(rank(difference) == 1 | rank(-difference) == 1))
    filtered.shapiro <- filtered.1 %>% filter(!(not.normal & extreme.guide))
    
    # re-do shapiro test for normality having removed extreme guides for non-normal groups
    filtered.shapiro <- filtered.shapiro %>% dplyr::select(-not.normal, -extreme.guide)
    filtered.shapiro <- filtered.shapiro %>% mutate(new.n=length(guide)) %>% filter(new.n >= min_guides_per_gene) %>% 
      mutate(new.shapiro=shapiro.test(difference)$p.value)
    
    ret <- filtered.shapiro
  } else {
    ret <- filtered.1
  }
  ret
}

gen.pvals <- function (dat) {
  ttests <- dat %>% summarize(median=median(difference), mean=mean(difference), n=n(), variance=var(difference)) %>% 
    mutate(stderror=sqrt(variance/n), statistic = mean/stderror, p.value= 2*pt(-abs(statistic), n - 1))
  ttests$var1 <- as.character(ttests$var1)
  ttests$var2 <- as.character(ttests$var2)
  ttests
}

add_annot <- function(df) {
  data.table(left_join(annotations.df, df, by="gene_symbol"))
}

get_gene_fcs <- function(gene_name) {
  add_annot(fold.changes) %>% filter(gene_symbol == gene_name)
}

get_gene_pvals <- function(gene_name) {
  add_annot(p.values) %>% filter(gene_symbol == gene_name)
}

get_gene_qvals <- function(gene_name) {
  add_annot(q.values) %>% filter(gene_symbol == gene_name)
}

get_guide_fcs <- function(gene_name) {
  add_annot(foldChanges_summ) %>% filter(gene_symbol == gene_name)
}





