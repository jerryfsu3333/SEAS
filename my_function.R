########################################################################
# For CovSir
########################################################################
results_covsir <- matrix(0, 2, dim(a)[2], dimnames = list(c('mean', 'sd'), colnames(a)))
for (i in seq_len(dim(a)[2])){
  tmp <- na.omit(a[,i])
  if(i %in% c(1,2)){
    results_covsir[1,i] <- mean(tmp)*100
  }else{
    results_covsir[1,i] <- mean(tmp)
  }
  results_covsir[2,i] <- sd(tmp)/sqrt(length(tmp))*100
}
results_covsir <- round(results_covsir, 3)

r_ratio <- sapply(c(3), function(i){
  tmp <- na.omit(a[,i])
  mean(tmp == 2)*100
  # mean(tmp == 1)*100
})
r_ratio <- round(r_ratio, 3)

cat(results_covsir['mean', 'distord_CovSIR'], '(', results_covsir['sd', 'distord_CovSIR'], ')', '&', 
    results_covsir['mean', 'r_CovSIR'], '(', results_covsir['sd', 'r_CovSIR'], ')', '&', 
    results_covsir['mean', 'C_CovSIR'], '(', results_covsir['sd', 'C_CovSIR'], ')', '&', 
    results_covsir['mean', 'IC_CovSIR'], '(', results_covsir['sd', 'IC_CovSIR'], ')\n', sep = '')

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
  mean(tmp == 2)*100
  # mean(tmp == 1)*100
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


#############################################################################################
# Execution time plot
#############################################################################################

library(reshape2)
library(ggplot2)

url <- "/Users/cengjing/Documents/DIS/Record/2019/Nov2_2019/output"
url_covsir <- "/Users/cengjing/Documents/DIS/Record/2019/Nov2_2019/covsir"
M <- c(1,2,4,6)
m <- c(1,5,10)

for (index1 in M){
  
  times <- c()

  for (index2 in m){
    a1 <- read.table(paste0(url, index1,"_", index2))
    a2 <- read.table(paste0(url_covsir, index1,"_", index2))
    
    results <- matrix(0, 2, dim(a1)[2], dimnames = list(c('mean', 'sd'), colnames(a1)))
    for (i in seq_len(dim(a1)[2])){
      tmp <- na.omit(a1[,i])
      if((i%%6 %in% c(1,2))){
        results[1,i] <- mean(tmp)*100
      }else{
        results[1,i] <- mean(tmp)
      }
      results[2,i] <- sd(tmp)/sqrt(length(tmp))*100
    }
    results <- round(results,3)
    
    results_covsir <- matrix(0, 2, dim(a2)[2], dimnames = list(c('mean', 'sd'), colnames(a2)))
    for (i in seq_len(dim(a2)[2])){
      tmp <- na.omit(a2[,i])
      if(i %in% c(1,2)){
        results_covsir[1,i] <- mean(tmp)*100
      }else{
        results_covsir[1,i] <- mean(tmp)
      }
      results_covsir[2,i] <- sd(tmp)/sqrt(length(tmp))*100
    }
    results_covsir <- round(results_covsir, 3)
    
    times <- cbind(times, c(results[1,seq(6,dim(results)[2], 6)], results_covsir[1,6]))
  }
  # rownames(times) <- c("SEAS-SIR", "SEAS-Intra", "SEAS-PFC", "Lasso-SIR", "Lasso", "Rifle-SIR", "Convex-SIR")
  ## Just record SEAS=SIR
  # times <- times[c("SEAS-SIR", "Lasso-SIR", "Lasso", "Rifle-SIR", "Convex-SIR"),]
  rownames(times) <- c("SEAS-SIR", "SEAS-Intra", "SEAS-PFC", "Lasso-SIR", "Convex-SIR")
  times <- times[c("SEAS-SIR", "Lasso-SIR", "Convex-SIR"),]
  
  times <- melt(times, value.name = "Time")
  times["Time"] <- log(times["Time"])
  
  g <- ggplot(times, aes(x=Var2, y=Time))+
    geom_line(aes(linetype=Var1))+
    # scale_linetype_manual(values = c("solid", "dashed", "twodash", "dotted", "longdash"))+
    scale_linetype_manual(values = c("solid", "twodash", "dotted"))+
    geom_point()+
    xlab("p")+
    ylab("Log(Time)")+
    theme(
      plot.title = element_text(size=15),
      axis.title.x = element_text(size=15),
      axis.title.y = element_text(size=15),
      legend.title = element_text(size=15),
      legend.text = element_text(size=12),
      legend.key.width = unit(1,"cm")
    )+
    scale_x_continuous(breaks = c(1,2,3), labels = c("100", "500", "1000")) +
    # labs(linetype= "Methods", title = paste0("Model (M", index1,")"))
    labs(linetype= "Methods", title = paste0("Model (M", 4,")"))
  g
}

#############################################################################################
# estimation
#############################################################################################

library(xtable)
library(ggplot2)
library(pracma)
boot_rank <- do.call(rbind, lapply(output, '[[', 1))
boot_s <- do.call(rbind, lapply(output, '[[', 2))
boot_nz <- lapply(output, '[[', 3)
boot_dist <- do.call(rbind, lapply(output, '[[', 4))
ord <- lapply(output, '[[', 5)
ord_est <- lapply(output, '[[', 6)

result1 <- apply(boot_dist, 2, function(x){mean(x, na.rm = TRUE)})
std1 <- apply(boot_dist, 2, function(x){sd(x, na.rm = TRUE)/sqrt(nrow(boot_dist))*100})

result2 <- apply(boot_rank, 2, function(x){mean(x, na.rm = TRUE)})
std2 <- apply(boot_rank, 2, function(x){sd(x, na.rm = TRUE)/sqrt(nrow(boot_rank))*100})

result3 <- apply(boot_s, 2, function(x){mean(x, na.rm = TRUE)})
std3 <- apply(boot_s, 2, function(x){sd(x, na.rm = TRUE)/sqrt(nrow(boot_s))*100})

tab <- rbind(paste0(round(result1,3), '(', round(std1,3), ')'), paste0(round(result2,3), '(', round(std2,3), ')'),
      paste0(round(result3,3), '(', round(std3,3), ')'))

colnames(tab) <- c('SEAS-SIR', 'Lasso-SIR', 'Lasso', 'Rifle-SIR')
# colnames(tab) <- c('SEAS-SIR', 'SEAS-Intra', 'SEAS-PFC', 'Lasso-SIR')
xtable(tab, 
       caption = 'The averaged subspace distance $\\bar{\\calD}$, the averaged rank $\\bar{d}$, the averaged sparsity $\\bar{s}$ and the corresponding standard error ($\\times 10^{-2}$) based on original samples and 100 bootstrap samples.',
       align = rep('c', dim(tab)[2]+1))


#############################################################################################
# Appearance frequency plot & cumulative frequency plot
#####################################.########################################################
library(latex2exp)
# 1. Top 10 selected variables
titles = c('SEAS-SIR', 'Lasso-SIR', 'Lasso', 'Rifle-SIR')
# titles = c('SEAS-SIR', 'SEAS-Intra', 'SEAS-PFC', 'Lasso-SIR')
true_nz <- true_output[[3]]
# true_nz <- true_output[[5]]
for (k in 1:4){
  true <- true_nz[[k]]
  nz <- lapply(boot_nz, '[[', k)
  tot <- 103
  # tot <- 102
  # tot <- 100
  # nz <- lapply(ord_est, '[[', k)
  cst <- mean(sapply(nz, function(x){
    if (!is.numeric(x)){
      FALSE
    }else{
      # SMS
      1-(length(union(true, x)) - length(intersect(true, x)))/tot
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
    labs(title = TeX(paste0(titles[k],' ($\\bar{\\mathrm{SMC}}$ = ', signif(cst,3), '%)')))
  print(g)
}

# 2. Cumulative apperance
titles = c('SEAS-SIR', 'SEAS-intra', 'SEAS-PFC', 'Lasso-SIR')
# titles = c('SEAS-SIR', 'Lasso-SIR', 'Lasso', 'Rifle-SIR')
for (k in 1:4){
  # nz <- lapply(ord_est, '[[', k)
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
      plot.title = element_text(size=12),
      axis.title.x = element_text(size=12),
      axis.title.y = element_text(size=12)
    )+
    labs(title = paste0(titles[k], ' (AUC = ', area, ')'))
    # theme(axis.text.x=element_text(angle=45, size = 4, vjust = 0.8))
  print(g)
}

################################################################################################
# Scatterplot
################################################################################################
# Use set.seed(123) for NIR_meat_cont2.R

library(latex2exp)
new_sir <- x %*% directions_sir
new_intra <- x %*% directions_intra
new_pfc <- x %*% directions_pfc
new_lassosir <- x %*% directions_lassosir

i <- 1
data <- data.frame(x = rbind(new_sir[,i,drop=FALSE], new_intra[,i,drop=FALSE], new_pfc[,i,drop=FALSE], new_lassosir[,i,drop=FALSE]),
                   y = rep(y, 4),
                   class = factor(c(rep("SEAS_SIR", nrow(x)), rep("SEAS-intra", nrow(x)), rep("SEAS-PFC", nrow(x)), rep("Lasso-SIR", nrow(x))), levels = c("SEAS_SIR", "SEAS-intra", "SEAS-PFC", "Lasso-SIR")))

g <- ggplot(data, aes(x = x, y = y)) + 
  geom_point() + 
  facet_wrap(~class, nrow = 2, ncol = 2, shrink = FALSE, scales = "free_x") +
  xlab(TeX('$\\beta_1^T\\mathbf{X}$')) +
  ylab(TeX('$\\mathbf{Y}$'))+
  theme_bw()+
  labs(title="Response versus the first component")+
  # labs(title="Response versus the second component")+
  theme(plot.title = element_text(size = 16))
g

# levels = c("SEAS_SIR", "SEAS-intra", "SEAS-PFC", "Lasso-SIR")


################################################################################################
# Histgram
################################################################################################

hist_plot <- function(x, y, title){
  if(!is.factor(y)){y <- factor(y)}
  if(!is.null(dim(x))){x <- drop(x)}
  df <- data.frame(component = x, Type = y)
  means <- sapply(unique(df$Type), function(i){
    mean(df$component[df$Type==i])
  })
  means_df <- data.frame(means = means, Type=unique(df$Type))
  g <- ggplot(df, aes(x=component, color = Type, fill=Type)) +
    # geom_histogram(aes(y=..density..), bins = 50, position = 'identity', alpha=0.5) +
    geom_density(alpha=0.3) +
    geom_vline(data = means_df, aes(xintercept=means, color=Type), linetype='dashed')+
    theme(
      plot.title = element_text(size=16),
      axis.title.x = element_text(size=16),
      axis.title.y = element_text(size=16)
    )+
    scale_fill_grey(end = 0.7)+
    scale_colour_grey(end = 0.7)+
    theme_bw()+
    labs(title = title)
  g
}

x_new <- as.matrix(x) %*% directions_sir
levels(y) <- c("Pork", "Beef")
hist_plot(x_new, y, 'SEAS-SIR')

x_new_lassosir <- as.matrix(x) %*% directions_lassosir
hist_plot(x_new_lassosir, y, 'Lasso-SIR')

x_new_lasso <- as.matrix(x) %*% directions_lasso
hist_plot(x_new_lasso, y, 'Lasso')

x_new_rifle <- as.matrix(x) %*% directions_rifle
hist_plot(x_new_rifle, y, 'Rifle-SIR')


################################################################################################
# Prediction
################################################################################################

library(xtable)
err <- sapply(1:4, function(i){
  tmp <- do.call(rbind, lapply(err_list, '[[', i))
  colMeans(tmp, na.rm = TRUE)*100
})

std <- sapply(1:4, function(i){
  tmp <- do.call(rbind, lapply(err_list, '[[', i))
  apply(tmp, 2, function(x){sd(x, na.rm = TRUE)/sqrt(nrow(tmp))})*100
})

# # Classification
# tab <- matrix(paste0(round(err,2), '(', round(std,2), ')'), 4,4)
# colnames(tab) <- c('SEAS-SIR', 'Lasso-SIR', 'Lasso', 'Rifle-SIR')
# rownames(tab) <- c('Logistic', 'SVM', 'LDA', 'Random Forest')
# xtable(tab, caption = 'The classification error rate (\\%) and standard error ($\\times 10^{-2}$) based on 100 replicates.',
#        align = rep('c', dim(err)[2]+1))

## Regression
tab <- matrix(paste0(round(err,2), '(', round(std,2), ')'), 4,4)
colnames(tab) <- c('SEAS-SIR', 'SEAS-intra', 'SEAS-PFC', 'Lasso-SIR')
rownames(tab) <- c('Linear', 'Kernel (l.c.)', 'Random Forest', 'SVM')
xtable(tab, caption = 'The mean square error ($\\times 10^{-2}$) and standard error ($\\times 10^{-2}$) based on 100 replicates.',
       align = rep('c', dim(err)[2]+1))
