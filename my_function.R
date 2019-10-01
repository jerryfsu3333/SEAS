#############################
a[,c(4,9,14,19,24,29)] <- a[,c(4,9,14,19,24,29)]/sqrt(2)
results <- matrix(0, 2, dim(a)[2], dimnames = list(c('mean', 'sd'), colnames(a)))
for (i in seq_len(dim(a)[2])){
  tmp <- na.omit(a[,i])
  if((i%%5 %in% c(1,2)) & (i <= 30)){
    results[1,i] <- mean(tmp)*100
  }else{
    results[1,i] <- mean(tmp)
  }
  results[2,i] <- sd(tmp)/sqrt(length(tmp))*100
}
results <- round(results,3)

# r_ratio <- sapply(c(3,8,13,18), function(i){
r_ratio <- sapply(c(3,8,13,18,23,28), function(i){
  tmp <- na.omit(a[,i])
  # mean(tmp == 2)*100
  mean(tmp == 1)*100
})
r_ratio <- round(r_ratio, 3)

cat(results[1,4], '(', results[2,4], ')', '&', results[1,31], '(', results[2,31], ')', '&', results[1,3], '(', results[2,3], ')', '&', r_ratio[1], '&', results[1,1], '(', results[2,1], ')', '&', results[1,2], '(', results[2,2], ')', '&', results[1,5], '(', results[2,5], ')\n', sep = '')
cat(results[1,9], '(', results[2,9], ')', '&', results[1,32], '(', results[2,32], ')', '&', results[1,8], '(', results[2,8], ')', '&', r_ratio[2], '&', results[1,6], '(', results[2,6], ')', '&', results[1,7], '(', results[2,7], ')', '&', results[1,10], '(', results[2,10], ')\n',sep = '')
cat(results[1,14], '(', results[2,14], ')', '&', results[1,33], '(', results[2,33], ')', '&', results[1,13], '(', results[2,13], ')', '&', r_ratio[3], '&', results[1,11], '(', results[2,11], ')', '&', results[1,12], '(', results[2,12], ')', '&', results[1,15], '(', results[2,15], ')\n', sep = '')
cat(results[1,19], '(', results[2,19], ')', '&', results[1,34], '(', results[2,34], ')', '&', results[1,18], '(', results[2,18], ')', '&', r_ratio[4], '&', results[1,16], '(', results[2,16], ')', '&', results[1,17], '(', results[2,17], ')', '&', results[1,20], '(', results[2,20], ')\n', sep = '')
cat(results[1,24], '(', results[2,24], ')', '&', results[1,35], '(', results[2,35], ')', '&', results[1,23], '(', results[2,23], ')', '&', r_ratio[5], '&', results[1,21], '(', results[2,21], ')', '&', results[1,22], '(', results[2,22], ')', '&', results[1,25], '(', results[2,25], ')\n', sep = '')
cat(results[1,29], '(', results[2,29], ')', '&', results[1,36], '(', results[2,36], ')', '&', results[1,28], '(', results[2,28], ')', '&', r_ratio[6], '&', results[1,26], '(', results[2,26], ')', '&', results[1,27], '(', results[2,27], ')', '&', results[1,30], '(', results[2,30], ')\n', sep = '')
####################################
# For CovSir
results <- matrix(0, 2, dim(a)[2], dimnames = list(c('mean', 'sd'), colnames(a)))
for (i in seq_len(dim(a)[2])){
  tmp <- na.omit(a[,i])
  if(i %in% c(1,2)){
    results[1,i] <- mean(tmp)*100
  }else{
    results[1,i] <- mean(tmp)
  }
  results[2,i] <- sd(tmp)/sqrt(length(tmp))*100
}
results <- round(results, 3)

r_ratio <- sapply(c(3), function(i){
  tmp <- na.omit(a[,i])
  mean(tmp == 2)*100
  # mean(tmp == 1)*100
})
r_ratio <- round(r_ratio, 3)

cat(results['mean', 'distord_CovSIR'], '(', results['sd', 'distord_CovSIR'], ')', '&', 
    results['mean', 'r_CovSIR'], '(', results['sd', 'r_CovSIR'], ')', '&', 
    results['mean', 'C_CovSIR'], '(', results['sd', 'C_CovSIR'], ')', '&', 
    results['mean', 'IC_CovSIR'], '(', results['sd', 'IC_CovSIR'], ')\n', sep = '')

for (i in 1:ncol(a)){a[,i] <- as.numeric(levels(a[,i]))[a[,i]]}
#############################
# 
results <- matrix(0, 2, dim(a)[2], dimnames = list(c('mean', 'sd'), colnames(a)))

for (i in seq_len(dim(a)[2])){
  tmp <- na.omit(a[,i])
  if((i%%5 %in% c(1,2)) & (i <=20)){
    results[1,i] <- mean(tmp)*100
  }else{
    results[1,i] <- mean(tmp)
  }
  
  results[2,i] <- sd(tmp)/sqrt(length(tmp))*100
}
results <- round(results,2)

r_ratio <- sapply(c(3,8,13,18), function(i){
  # r_ratio <- sapply(c(3,8,13,18,23,28), function(i){
  tmp <- na.omit(a[,i])
  mean(tmp == 2)*100
  # mean(tmp == 1)*100
})
r_ratio <- round(r_ratio, 2)

cat(results[1,4], '(', results[2,4], ')', '&', results[1,21], '(', results[2,21], ')', '&', results[1,3], '(', results[2,3], ')', '&', r_ratio[1], '&', results[1,1], '(', results[2,1], ')',
    '&', results[1,2], '(', results[2,2], ')', '&', results[1,5], '(', results[2,5], ')\n', sep = '')
cat(results[1,9], '(', results[2,9], ')', '&', results[1,22], '(', results[2,22], ')', '&', results[1,8], '(', results[2,8], ')', '&', r_ratio[2], '&', results[1,6], '(', results[2,6], ')', 
    '&', results[1,7], '(', results[2,7], ')', '&', results[1,10], '(', results[2,10], ')\n',sep = '')
cat(results[1,14], '(', results[2,14], ')', '&', results[1,23], '(', results[2,23], ')', '&', results[1,13], '(', results[2,13], ')', '&', r_ratio[3], '&', results[1,11], '(', results[2,11], ')', 
    '&', results[1,12], '(', results[2,12], ')', '&', results[1,15], '(', results[2,15], ')\n', sep = '')
cat(results[1,19], '(', results[2,19], ')', '&', results[1,24], '(', results[2,24], ')', '&', results[1,18], '(', results[2,18], ')', '&', r_ratio[4], '&', results[1,16], '(', results[2,16], ')', 
    '&', results[1,17], '(', results[2,17], ')', '&', results[1,20], '(', results[2,20], ')\n', sep = '')

########################################

a[,c(4,9,14,19)] <- a[,c(4,9,14,19)]/sqrt(2)
results <- matrix(0, 2, dim(a)[2], dimnames = list(c('mean', 'sd'), colnames(a)))

for (i in seq_len(dim(a)[2])){
  tmp <- na.omit(a[,i])
  if((i%%5 %in% c(1,2)) & (i <=20)){
    results[1,i] <- mean(tmp)*100
  }else{
    results[1,i] <- mean(tmp)
  }
  
  results[2,i] <- sd(tmp)/sqrt(length(tmp))*100
}
results <- round(results,3)

r_ratio <- sapply(c(3,8,13,18), function(i){
  # r_ratio <- sapply(c(3,8,13,18,23,28), function(i){
  tmp <- na.omit(a[,i])
  mean(tmp == 2)*100
  # mean(tmp == 1)*100
})
r_ratio <- round(r_ratio, 3)

cat(results[1,4], '(', results[2,4], ')', '&', results[1,21], '(', results[2,21], ')', '&', results[1,3], '(', results[2,3], ')', '&', r_ratio[1], '&', results[1,1], '(', results[2,1], ')',
    '&', results[1,2], '(', results[2,2], ')', '&', results[1,5], '(', results[2,5], ')\n', sep = '')
cat(results[1,9], '(', results[2,9], ')', '&', results[1,22], '(', results[2,22], ')', '&', results[1,8], '(', results[2,8], ')', '&', r_ratio[2], '&', results[1,6], '(', results[2,6], ')', 
    '&', results[1,7], '(', results[2,7], ')', '&', results[1,10], '(', results[2,10], ')\n',sep = '')
cat(results[1,14], '(', results[2,14], ')', '&', results[1,23], '(', results[2,23], ')', '&', results[1,13], '(', results[2,13], ')', '&', r_ratio[3], '&', results[1,11], '(', results[2,11], ')', 
    '&', results[1,12], '(', results[2,12], ')', '&', results[1,15], '(', results[2,15], ')\n', sep = '')
cat(results[1,19], '(', results[2,19], ')', '&', results[1,24], '(', results[2,24], ')', '&', results[1,18], '(', results[2,18], ')', '&', r_ratio[4], '&', results[1,16], '(', results[2,16], ')', 
    '&', results[1,17], '(', results[2,17], ')', '&', results[1,20], '(', results[2,20], ')\n', sep = '')


################################################################################

results <- matrix(0, 2, dim(a)[2], dimnames = list(c('mean', 'sd'), colnames(a)))

for (i in seq_len(dim(a)[2])){
  tmp <- na.omit(a[,i])
  if((i%%6 %in% c(1,2))){
    results[1,i] <- mean(tmp)*100
  }else{
    results[1,i] <- mean(tmp)
  }
  
  results[2,i] <- sd(tmp)/sqrt(length(tmp))*100
}
results <- round(results,3)

r_ratio <- sapply(c(3,9,15,21), function(i){
# r_ratio <- sapply(c(3,9,15,21,27,33), function(i){
  tmp <- na.omit(a[,i])
  # mean(tmp == 2)*100
  mean(tmp == 1)*100
})
r_ratio <- round(r_ratio, 3)

cat(results['mean', 'distord_ssdrsir'], '(', results['sd', 'distord_ssdrsir'], ')', '&', 
    results['mean', 'r_ssdrsir'], '(', results['sd', 'r_ssdrsir'], ')', '&', 
    results['mean', 'C_ssdrsir'], '(', results['sd', 'C_ssdrsir'], ')', '&', 
    results['mean', 'IC_ssdrsir'], '(', results['sd', 'IC_ssdrsir'], ')\n', sep = '')

cat(results['mean', 'distord_ssdrintra'], '(', results['sd', 'distord_ssdrintra'], ')', '&', 
    results['mean', 'r_ssdrintra'], '(', results['sd', 'r_ssdrintra'], ')', '&', 
    results['mean', 'C_ssdrintra'], '(', results['sd', 'C_ssdrintra'], ')', '&', 
    results['mean', 'IC_ssdrintra'], '(', results['sd', 'IC_ssdrintra'], ')\n', sep = '')

cat(results['mean', 'distord_ssdrpfc'], '(', results['sd', 'distord_ssdrpfc'], ')', '&', 
    results['mean', 'r_ssdrpfc'], '(', results['sd', 'r_ssdrpfc'], ')', '&', 
    results['mean', 'C_ssdrpfc'], '(', results['sd', 'C_ssdrpfc'], ')', '&', 
    results['mean', 'IC_ssdrpfc'], '(', results['sd', 'IC_ssdrpfc'], ')\n', sep = '')

cat(results['mean', 'distord_LassoSIR'], '(', results['sd', 'distord_LassoSIR'], ')', '&', 
    results['mean', 'r_LassoSIR'], '(', results['sd', 'r_LassoSIR'], ')', '&', 
    results['mean', 'C_LassoSIR'], '(', results['sd', 'C_LassoSIR'], ')', '&', 
    results['mean', 'IC_LassoSIR'], '(', results['sd', 'IC_LassoSIR'], ')\n', sep = '')

cat(results['mean', 'distord_lasso'], '(', results['sd', 'distord_lasso'], ')', '&',
    results['mean', 'r_lasso'], '(', results['sd', 'r_lasso'], ')', '&',
    results['mean', 'C_lasso'], '(', results['sd', 'C_lasso'], ')', '&',
    results['mean', 'IC_lasso'], '(', results['sd', 'IC_lasso'], ')\n', sep = '')

cat(results['mean', 'distord_rifle'], '(', results['sd', 'distord_rifle'], ')', '&',
    results['mean', 'r_rifle'], '(', results['sd', 'r_rifle'], ')', '&',
    results['mean', 'C_rifle'], '(', results['sd', 'C_rifle'], ')', '&',
    results['mean', 'IC_rifle'], '(', results['sd', 'IC_rifle'], ')\n', sep = '')


######################################
# estimation
library(ggplot2)
library(pracma)
boot_rank <- do.call(rbind, lapply(output, '[[', 1))
boot_s <- do.call(rbind, lapply(output, '[[', 2))
boot_dist <- do.call(rbind, lapply(output, '[[', 3))
boot_nz <- lapply(output, '[[', 4)
ord <- lapply(output, '[[', 5)
ord_est <- lapply(output, '[[', 6)

result <- apply(boot_rank, 2, function(x){c(mean(x, na.rm = TRUE), sd(x, na.rm = TRUE)/sqrt(length(x)) * 100)})
rownames(result) <- c('mean', 's.e.')
colnames(result) <- c('ssdrsir', 'lassosir', 'lasso', 'riflesir')
cat(result['mean', 'ssdrsir'], '&', result['mean', 'lassosir'], '&',result['mean', 'lasso'], '&', result['mean', 'riflesir'], '\n')
cat(result['s.e.', 'ssdrsir'], '&', result['s.e.', 'lassosir'], '&',result['s.e.', 'lasso'], '&', result['s.e.', 'riflesir'], '\n')
# colnames(result) <- c('ssdrsir', 'ssdrintra', 'ssdrpfc', 'lassosir')
# cat(result['mean', 'ssdrsir'], '&', result['mean', 'ssdrintra'], '&', result['mean', 'ssdrpfc'], '&', result['mean', 'lassosir'], '\n')
# cat(result['s.e.', 'ssdrsir'], '&', result['s.e.', 'ssdrintra'], '&', result['s.e.', 'ssdrpfc'], '&', result['s.e.', 'lassosir'], '\n')

result <- apply(boot_dist, 2, function(x){c(mean(x, na.rm = TRUE), sd(x, na.rm = TRUE)/sqrt(length(x)) * 100)})
rownames(result) <- c('mean', 's.e.')
colnames(result) <- c('ssdrsir', 'lassosir', 'lasso', 'riflesir')
cat(result['mean', 'ssdrsir'], '&', result['mean', 'lassosir'], '&',result['mean', 'lasso'], '&', result['mean', 'riflesir'], '\n')
cat(result['s.e.', 'ssdrsir'], '&', result['s.e.', 'lassosir'], '&',result['s.e.', 'lasso'], '&', result['s.e.', 'riflesir'], '\n')


# 1. Top 10 selected variables
# titles = c('SSDR-SIR', 'SSDR-intra', 'SSDR-PFC', 'Lasso-SIR')
titles = c('SSDR-SIR', 'Lasso-SIR', 'Lasso', 'Rifle-SIR')
true_nz <- true_output[[3]]
for (k in 1:4){
  true <- true_nz[[k]]
  nz <- lapply(boot_nz, '[[', k)
  cst <- mean(sapply(nz, function(x){
    if (!is.numeric(x)){
      FALSE
    }else{
      all(true %in% x)
    }
  })) * 100
  cst <- round(cst, digits = 3)
  
  tab <- sort(table(do.call(c, nz)), decreasing = TRUE)
  name <- as.integer(names(tab))
  
  # Select the top 10
  freq <- as.vector(tab[1:min(10,length(tab))]/length(nz) * 100)
  name <- name[1:min(10,length(name))]
  
  g <- ggplot(data.frame(x=1:length(name), y=freq), aes(x=x, y=y, group=1))+
    geom_line()+
    geom_point(size=1)+
    ylim(c(0,100))+
    scale_x_continuous(breaks = 1:length(name), labels = name)+
    xlab('variable index')+
    ylab('Appearance frequency (%)')+
    theme(
      plot.title = element_text(size=16),
      axis.title.x = element_text(size=16),
      axis.title.y = element_text(size=16)
    )+
    labs(title = paste0(titles[k], ' (consistency = ', cst, '%)'))
  print(g)
}

# 2. Cumulative apperance
# titles = c('SSDR-SIR', 'SSDR-intra', 'SSDR-PFC', 'Lasso-SIR')
titles = c('SSDR-SIR', 'Lasso-SIR', 'Lasso', 'Rifle-SIR')
for (k in 1:4){
  nz <- lapply(boot_nz, '[[', k)
  tab <- sort(table(do.call(c, nz)), decreasing = TRUE)
  name <- as.integer(names(tab))
  freq <- sapply(seq_len(length(name)), function(i){
    ind <- name[1:i]
    num <- mean(sapply(nz, function(x){
      if (!is.numeric(x)){
        0
      }else{
        all(ind %in% x)
      }
    })) * 100
  })
  
  # select larger than 10%
  freq <- freq[1:min(10,length(freq))]
  name <- name[1:min(10,length(name))]
  area <- trapz(1:length(freq), freq)/100
  
  g <- ggplot(data.frame(x=1:length(name), y=freq), aes(x=x, y=y, group = 1))+
    geom_line()+
    geom_point(size=1)+
    scale_x_continuous(breaks = 1:length(name), labels = name)+
    ylim(c(0,100))+
    xlab('variable index')+
    ylab('Cumulative appearance frequency (%)')+
    theme(
      plot.title = element_text(size=16),
      axis.title.x = element_text(size=16),
      axis.title.y = element_text(size=16)
    )+
    labs(title = paste0(titles[k], ' (AUC = ', area, ')'))
    # theme(axis.text.x=element_text(angle=45, size = 4, vjust = 0.8))
  print(g)
}


##########################################
## predictions
library(xtable)
err <- sapply(1:4, function(i){
  err <- do.call(rbind, lapply(err_list, '[[', i))
  colMeans(err, na.rm = TRUE)
})

colnames(err) <- c('SSDR-SIR', 'Lasso-SIR', 'Lasso', 'Rifle-SIR')
rownames(err) <- c('Logistic', 'SVM', 'LDA', 'Random Forest')
xtable(err, display = c('d', rep('e',dim(err)[2])),
       caption = 'The classification error rate (\\%) for different methods under different classifiers based on 100 replicates.',
       align = rep('c', dim(err)[2]+1))

# colnames(err) <- c('SSDR-SIR', 'SSDR-intra', 'SSDR-PFC', 'Lasso-SIR')
# rownames(err) <- c('Kernel (l.c.)', 'Random Forest', 'SVM')
# xtable(err, display = c('d', rep('e',dim(err)[2])),
#        caption = 'The mean square error for different methods under different classifiers based on 100 replicates.',
#        align = rep('c', dim(err)[2]+1))

