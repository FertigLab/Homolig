#Using matrix factorizations of AADist, describe any group of AA sequences.
ScorePatterns3 <- function(input, metric = 'nmf', error = 'se', pattern_names = NA,
                           interpolate_distance = NA, plotSig = TRUE){
  res = GetPatternScores2(input, metric = metric, pattern_names = pattern_names, 
                          interpolate_distance = interpolate_distance)
  plotres = PlotPatternScores(res, error = error, plotSig = plotSig)
  return(plotres)
}

AverageScoreReport <- function(input, metric = 'ensemble', error = 'se', pattern_names = NA){
  scores = ScorePatterns3(input, metric, error, pattern_names)$data
  fmean = tapply(scores$Score, scores$Factor, mean)
  fsd = tapply(scores$Score, scores$Factor, sd)
  fse = fsd/sqrt(length(input))
  res =data.frame(factor = names(fmean), 
             mean.score = as.numeric(fmean), 
             sd = as.numeric(fsd), 
             se = as.numeric(fse), 
             var = fsd^2)
  rownames(res) = NULL
  return(res)
}
ScoreLogos <- function(input, metric = 'nmf'){
  require(msa)
  require(ggseqlogo)
  pmat = readMatrix2(metric)
  subMat = makeSubstitutionMatrix(pmat)
 x = msa(input, method = 'ClustalW', cluster = 'nj', #substitutionMatrix = 'id', 
         type = 'protein', gapOpening=10, gapExtension=0.2, 
         substitutionMatrix = subMat)
  x = as.character(x)
  gg = ggplot() + geom_logo(x) + theme_logo()
  return(gg)
}
process_logo <- function(sequence_info, metric = 'nmf'){
  require(msa)
  require(ggseqlogo)
  pmat = readMatrix2(metric)
  subMat = makeSubstitutionMatrix(pmat)
  #hdf is a Homolig result object
  if('group' %in% colnames(sequence_info)){
    sequence_info = sequence_info[which(sequence_info$group == 'a'),]
  }
  
  si = sequence_info
  
  seqs = list(
    cdr1 = AAStringSet(si$CDR1), 
    cdr2 = AAStringSet(si$CDR2),
    cdr2.5 = AAStringSet(si$CDR2.5),
    cdr3 = AAStringSet(si$CDR3)
  )
  
  aligned = lapply(seqs, msa, method = 'ClustalOmega',cluster = 'nj', #substitutionMatrix = 'id', 
                   type = 'protein', gapOpening=10, gapExtension=0.2, 
                   substitutionMatrix = subMat)
  aligned_logos = lapply(aligned, FUN = function(x){
    x = as.character(x)
    gg = ggplot() + geom_logo(x) + theme_logo()
    return(gg)
  })
  

  return(aligned_logos)
  
  
}
readMatrix2 <- function(metric){
  load('./factorization-data/2024-05-27_All-MF.rda')
  idx = which(names(allmats) == metric)
  return(allmats[[idx]])
}
PlotPatternScores <- function(pattern_score_data, error = 'se', plotSig = TRUE, ALPHA = 0.10){
  require(ggplot2)
  require(ggprism)
  require(scales) #pretty_breaks()
  out_df = pattern_score_data
  
  sevals = pattern_score_data$se
  sevals[is.na(sevals)] = 0 #for purposes of computing plot limit
  if(max(abs(pattern_score_data$Score) + sevals, na.rm = TRUE ) >= 1.5){
    YLIM = max(abs(pattern_score_data$Score) + pattern_score_data$se, na.rm = TRUE) + 0.1
  }else{
    YLIM=1.5
  }
  if(error == 'sd'){
    ggfacet = ggplot(out_df, aes(x = Position, y = Score, color = Factor)) + 
      geom_path() + 
      theme_prism() + 
      geom_errorbar(aes(ymin = Score - sd, ymax = Score + sd), color = 'grey35') +
      geom_hline(yintercept = 0, size = 1) +
      facet_wrap(~Factor, scales = 'free') + 
      theme(legend.position = 'none') + 
      scale_x_continuous(breaks= pretty_breaks()) + 
      ylim(-YLIM,YLIM)
  }else if(error == 'se'){
    ggfacet = ggplot(out_df, aes(x = Position, y = Score, color = Factor)) + 
      geom_path() + 
      theme_prism() + 
      geom_errorbar(aes(ymin = Score - se, ymax = Score + se), color = 'grey35') +
      geom_hline(yintercept = 0, size = 1) + 
      scale_x_continuous(breaks= pretty_breaks()) + 
      facet_wrap(~Factor, scales = 'free') + 
      theme(legend.position = 'none') + 
      ylim(-YLIM,YLIM)
  }
  
  if(plotSig){
    ggfacet = ggfacet + 
      geom_point(data = out_df[out_df$padj < ALPHA,], 
                 shape = 8, color = 'grey15')
  }
  ggall =  ggplot(out_df, aes(x = Position, y = Score, color = Factor)) + 
    geom_path() + 
    theme_prism() + 
    geom_hline(yintercept = 0, size = 1) 
  return(list(ggfacet= ggfacet, ggall = ggall, data = out_df))
}
PlotWholePatternScores <- function(inputs, metric = 'nmf', error = 'se', GroupFactor = NA, GroupCols = NA){
  require(ggplot2)
  require(ggprism)
  require(scales) #pretty_breaks()
  res = lapply(inputs, GetWholePatternScores, metric = metric)
  data = lapply(res, FUN = function(x) return(x$data))
  data = bind_rows(data, .id = 'group')
  stats = lapply(res, FUN = function(x) return(x$stats))
  stats = bind_rows(stats, .id = 'group')
  
  
  if(max(abs(data$score) , na.rm = TRUE ) >= 1.5){
    YLIM = max(abs(data$score), na.rm = TRUE) + 0.1
  }else{
    YLIM=1.5
  }
  if(!is.na(GroupFactor[1])){
  data$group = factor(data$group, levels = GroupFactor)
  }
  if(is.na(GroupCols[1])){
    GroupCols = rep('grey60', length(GroupFactor))
  }
  spl= split(data, data$variable)
  splstat = split(stats, stats$variable)
  plots = vector(mode = 'list', length = length(spl)); names(plots) = names(spl)
  for(i in 1:length(spl)){
    if(error == 'sd'){
      ggfacet = ggplot(spl[[i]], aes(x = group, y = score, fill = group)) + 
        theme_prism() + 
        geom_hline(yintercept = 0, size = 1) + 
        geom_violin() + 
        geom_errorbar(data = splstat[[i]], aes(ymin = score - sd, ymax = score + sd), color = 'grey35') +
      #  facet_wrap(~variable, scales = 'free') + 
        scale_fill_manual(values = GroupCols) + 
        theme(legend.position = 'none', 
              axis.text.x = element_text(size = 11))  + 
        xlab('') + ylab('') + labs(title = names(spl)[i])
    }else if(error == 'se'){
      ggfacet = ggplot(spl[[i]], aes(x = group, y = score, fill = group)) + 
        theme_prism() + 
        geom_hline(yintercept = 0, size = 1) + 
        geom_violin() + 
        geom_errorbar(data = splstat[[i]], aes(ymin = score - se, ymax = score + se), color = 'grey35') +
       # facet_wrap(~variable, scales = 'free') + 
        scale_fill_manual(values = GroupCols) + 
        theme(legend.position = 'none', 
              axis.text.x = element_text(size = 11)) + 
        xlab('') + ylab('') + labs(title = names(spl)[i])
    }
    
    plots[[i]] = ggfacet
  }
 
  
   
  return(list(plots = plots, data = data, summary = stats))
}
GetPatternScores2 <- function(input, metric = 'ensemble', ALPHA = 0.1, pattern_names = NA, interpolate_distance = NA){
  #Interpolate sequence distance to avoid downstream length-dependent artifacts in msa. 
  #10 October 2023. 
  lentab = table(nchar(input))
  if(is.na(interpolate_distance)){
    #if not specified, select the median of lengths: 
    #interpolate_distance = round(median(nchar(input)))
    #interpolate_distance = as.numeric(names(lentab))[which.max(lentab)[1]] #selects mode
    interpolate_distance = 10 #hard-code 10
  }
  pmat = readMatrix2(metric)
  subMat = makeSubstitutionMatrix(pmat)
  pmat = apply(pmat, 2, FUN = function(x){
    x = (x - mean(x)) / sd(x)
    return(x)
  })
  if(!is.na(pattern_names[1]) & length(pattern_names) == ncol(pmat)){
    colnames(pmat) = pattern_names
  }
  input = toupper(input)
  aln = strsplit(input, split = '')
  nchars = nchar(input)
  
  spl = split(aln, nchars) #split by sequence length. 
  
  spl = spl[order(as.numeric(names(spl)))] #just in case list is not ordered
  lentab = lentab[order(as.numeric(names(lentab)))]
  
  res_by_len = vector(mode = 'list', length = length(spl))
  
  res_by_len = array(data = NA, dim = c(ncol(pmat), interpolate_distance, length(spl)), 
                     dimnames = list(Factor = colnames(pmat), Position = c(1:interpolate_distance), OrigLen = names(lentab)))
  
  sd_by_len = res_by_len
  
  for(s in 1:length(spl)){
    L = as.numeric(names(spl)[s]) #length of this group of sequences. 
    nSeq = length(spl[[s]])               
    w = ComputeInterpolationWeights(orig_len = L,
                                    final_len = interpolate_distance)
    
    m = matrix(data = NA, nrow = nSeq, ncol = L)
    for(r in 1:nSeq){
      m[r,] = spl[[s]][[r]] 
    }
    
    m[m=='-'] = NA
    
    if(nrow(m) > 1){
        out = apply(m, 2, FUN = function(x){
          idx = match(x, rownames(pmat))
          tmp = pmat[idx,]
          tmp = colSums(tmp, na.rm = TRUE) / nrow(tmp) #Reduces pattern score if NAs present. 
        })
        
        out_sd = apply(m, 2, FUN = function(x){
          idx = match(x, rownames(pmat))
          tmp = pmat[idx,]
          sd = apply(tmp, 2, FUN = function(y){
            return(sd(y, na.rm = TRUE))
          })
          return(sd) #colSums(tmp, na.rm = TRUE) / nrow(tmp)
        })
        
    }else{
      idx = match(m, rownames(pmat))
      out = t(pmat[idx,])
      colnames(out) = NULL
      
      out_sd = matrix(data = 0, nrow = nrow(out), ncol = ncol(out), 
                      dimnames = list(rownames(out), colnames(out)))
     
    }
    
    res_by_len[,,s] =  out %*% w #CRITICAL LINE!!!! Performing interpolation to new Length.
    sd_by_len[,,s] = out_sd %*% w

  }
  
  #Now combine results weighted by the number of sequences per 
  #length group. 
  #Thenk you to flodel: 
  #https://stackoverflow.com/questions/9222056/existing-function-to-combine-standard-deviations-in-r
  grand.mean <- function(M, N) {weighted.mean(M, N)}
  grand.sd   <- function(S, M, N) {sqrt(weighted.mean(S^2 + M^2, N) -
                                          weighted.mean(M, N)^2)}
  combined_results = matrix(data = NA, nrow = ncol(pmat), ncol = interpolate_distance, 
                           dimnames = list(Factor = colnames(pmat), Position = c(1:interpolate_distance)))
  combined_sd = combined_results
  
  for(i in 1:nrow(combined_results)){
    for(j in 1:ncol(combined_results)){
      combined_results[i,j] = grand.mean(M = res_by_len[i,j,], N = lentab)
      combined_sd[i,j] = grand.sd(S = sd_by_len[i,j,], M = res_by_len[i,j,], N = lentab)
    }
  }
  
  combined_se = combined_sd / sqrt(length(input)) #Generate SE after combining SDs. 
  
  outsd_df = reshape2::melt(combined_sd); colnames(outsd_df) = c('Factor', 'Position', 'SD')
  outse_df = reshape2::melt(combined_se); colnames(outse_df) = c('Factor', 'Position', 'SE')
  out_df = reshape2::melt(combined_results); colnames(out_df) = c('Factor', 'Position', 'Score')
  out_df$sd = outsd_df$SD
  out_df$se = outse_df$SE
  
  out_df$Factor = factor(out_df$Factor)
  out_df = GetPatternPValues(out_df, m, ALPHA = ALPHA)
  return(out_df)
}
GetWholePatternScores <- function(input, metric = 'nmf', ALPHA = 0.1, pattern_names = NA){
  pmat = readMatrix2(metric)
  subMat = makeSubstitutionMatrix(pmat)
  pmat = apply(pmat, 2, FUN = function(x){
    x = (x - mean(x)) / sd(x)
    return(x)
  })
  if(!is.na(pattern_names[1]) & length(pattern_names) == ncol(pmat)){
    colnames(pmat) = pattern_names
  }
  input = toupper(input)
  aln = strsplit(input, split = '')
  scores = lapply(aln, FUN = function(x){
    idx = match(x, rownames(pmat))
    tmp = pmat[idx,]
    tmp = colSums(tmp, na.rm = TRUE) / nrow(tmp) 
    return(tmp)
  })
  scoremat = as.data.frame(bind_rows(scores))
  stats = data.frame(score = apply(scoremat, 2, mean), 
                     variable = colnames(scoremat),
                     sd = apply(scoremat,2,sd), 
                     se = apply(scoremat,2,FUN = function(x){sd(x)/sqrt(nrow(scoremat))}))
  rownames(stats) = NULL
  scoremat$id = paste0('seq', c(1:nrow(scoremat)))
  df = suppressMessages(reshape2::melt(scoremat, value.name = 'score'))

 return(list(data = df, stats = stats ))
}
ComputeInterpolationWeights <- function(orig_len = 5, final_len = 3){
  #Used internally to interpolate pattern scores of variable length junction sequences. 
  #AAG 10 October 2023 
  bounds = seq(from = 0, to = orig_len, by = orig_len/final_len) + 0.5
  
  new_lower = bounds[1:length(bounds)-1]
  new_upper = bounds[2:length(bounds)]
  
  orig_lower = c(1:orig_len) - 0.5
  orig_upper = c(1:orig_len) + 0.5
  w = matrix(data = 0, nrow = orig_len, ncol = final_len)
  if(final_len == orig_len){
    diag(w) = 1
    return(w)
  }
  for(i in 1:orig_len){
    ilo = orig_lower[i]
    ihi = orig_upper[i]
    occ = 0 #proportion of old interval which has been delegated to new intervals
    j = 0
    while(occ < 1){
      j = j+1
      jlo = new_lower[j]
      jhi = new_upper[j]
      
      if(orig_len > final_len){ #if we must truncate orig sequence.
            if(ilo < jhi & ihi < jhi){
              #if entire range below upper limit
              value = 1 - occ
            }else if(ilo < jhi & ihi >= jhi){
              #if some of range below upper limit
              value =  jhi - ilo
            }else if(ilo >= jhi & ihi >= jhi){
              #if entire range above upper limit 
              value = 0
            }
      }else if(orig_len < final_len){ #if we must extend original sequence.
            if(ilo < jhi & ihi < jhi){
              #if entire range below upper limit
              value = 1 - occ
            }else if(ilo < jhi & ihi >= jhi){
              #if some of range below upper limit
              value =  jhi - ilo - occ
            }else if(ilo >= jhi & ihi >= jhi){
              #if entire range above upper limit 
              value = 0
            }
      }
      w[i,j] = value
      occ = occ + value
    }
  }
  return(w)
}
ComparePatternScores2 <- function(inputs, metric = 'ensemble', ALPHA = 0.1, pattern_names = NA, interpolate_distance = NA){
  #Compare pattern scores between multiple groups.
  #Interpolate sequence distance to avoid downstream length-dependent artifacts in msa. 
  #AAG 7 March 2024
  require(reshape2)
  require(dplyr)
  pmat = readMatrix2(metric)
  subMat = makeSubstitutionMatrix(pmat)
  pmat = apply(pmat, 2, FUN = function(x){
    x = (x - mean(x)) / sd(x)
    return(x)
  })
  if(!is.na(pattern_names[1]) & length(pattern_names) == ncol(pmat)){
    colnames(pmat) = pattern_names
  }
  
  grSizes = unlist(lapply(inputs, length))
  gr_names = names(inputs)
  nGr = length(grSizes)
  group.labels = rep(gr_names, grSizes)
  input = unlist(inputs)
  
  nSeq = length(input)
  input = toupper(input)
  
  tmp = GetPatternScores2(input[1], metric = metric, interpolate_distance = interpolate_distance) #placeholder
  nvar = nrow(tmp) #number of variables 
  
 temp.fun <- function(k){
   GetPatternScores2(k, interpolate_distance = interpolate_distance)$Score
 }
 res =  mclapply(input, FUN = temp.fun)
 resmat = as.matrix(bind_rows(res))
 if(nGr > 2){
 kwp =  apply(resmat, 1, FUN = function(r){
   kruskal.test(r, group.labels)$p.value
 })
 }else{
   #despite continuing to name variable kwp, two-variable version is a wilcox test. 
   kwp = apply(resmat, 1, FUN = function(r){
     wilcox.test(x = r[1:grSizes[1]], y = r[(grSizes[1] + 1): length(r)])$p.value
   })
 }
 
mag = apply(resmat,1, FUN = function(r){
   spl = split(r, group.labels)
   means = unlist(lapply(spl, mean))
   return(var(means)/var(r))
 })
 res = data.frame(
   Factor = tmp$Factor, 
   Position = tmp$Position, 
   Magnitude = mag,
   pval = kwp
 )
 res$padj = BY(res$pval, alpha = ALPHA, silent = TRUE)$adjPValues
 #res = res[order(res$pval, decreasing = FALSE),]
 res = res[order(res$Magnitude, decreasing = TRUE),]
 rownames(res) = NULL
 if(any(res$padj < ALPHA)){
   return(res[res$padj < ALPHA,])
 }else{
   return(res)
 }
}
CompareWholePatternScores <-  function(inputs, metric = 'ensemble', ALPHA = 0.1, pattern_names = NA, interpolate_distance = NA){
  #Compare pattern scores between multiple groups.
  #Interpolate sequence distance to avoid downstream length-dependent artifacts in msa. 
  #AAG 7 March 2024
  require(reshape2)
  require(dplyr)
  pmat = readMatrix2(metric)
  subMat = makeSubstitutionMatrix(pmat)
  pmat = apply(pmat, 2, FUN = function(x){
    x = (x - mean(x)) / sd(x)
    return(x)
  })
  if(!is.na(pattern_names[1]) & length(pattern_names) == ncol(pmat)){
    colnames(pmat) = pattern_names
  }
  
  grSizes = unlist(lapply(inputs, length))
  gr_names = names(inputs)
  nGr = length(grSizes)
  group.labels = rep(gr_names, grSizes)
  input = unlist(inputs)
  
  nSeq = length(input)
  input = toupper(input)
  
  tmp = AverageScoreReport(input[1], metric = metric) #placeholder
  nvar = nrow(tmp) #number of variables 
  
  res =  mclapply(input, FUN = function(k){AverageScoreReport(k)$mean.score})
  resmat = as.matrix(bind_rows(res))
  mag = apply(resmat,1, FUN = function(r){
    spl = split(r, group.labels)
    means = unlist(lapply(spl, mean))
    return(var(means)/var(r))
  })
  if(nGr > 2){
    kwp =  apply(resmat, 1, FUN = function(r){
      kruskal.test(r, group.labels)$p.value
    })
  }else{
    #despite continuing to name variable kwp, two-variable version is a wilcox test. 
    kwp = apply(resmat, 1, FUN = function(r){
      wilcox.test(x = r[1:grSizes[1]], y = r[(grSizes[1] + 1): length(r)])$p.value
    })
  }
  mag = 
  res = data.frame(
    Factor = tmp$factor, 
    Magnitude = mag, 
    pval = kwp
  )
  res$padj = BY(res$pval, alpha = ALPHA, silent = TRUE)$adjPValues
  #res = res[order(res$pval, decreasing = FALSE),]
  res = res[order(res$Magnitude, decreasing = TRUE),]
  rownames(res) = NULL
  if(any(res$padj < ALPHA)){
    return(res[res$padj < ALPHA,])
  }else{
    return(res)
  }
}
PlotPatternScoreComparison2 <- function(inputs, metric = 'ensemble', Factors, Positions, 
                                        pattern_names = NA, GroupFactor = NA, GroupCols = NA, interpolate_distance){
  #Factors is a vector of factor names and msa positions respectively of equal length. 
  #A plot will be generated for each combination of factor + position for each input group. 
  #inputs is a list of sequences by group, with list names corresponding to group names. 
  #12 September 2023 AAG 
  require(msa)
  require(reshape2)
  require(dplyr)
  require(ggplot2)
  require(ggpubr)
  require(ggprism)
  pmat = readMatrix2(metric)
  subMat = makeSubstitutionMatrix(pmat)
  pmat = apply(pmat, 2, FUN = function(x){
    x = (x - mean(x)) / sd(x)
    return(x)
  })
  if(!is.na(pattern_names[1]) & length(pattern_names) == ncol(pmat)){
    colnames(pmat) = pattern_names
  }
  
  grSizes = unlist(lapply(inputs, length))
  gr_names = names(inputs)
  nGr = length(grSizes)
  group.labels = rep(gr_names, grSizes)
  input = unlist(inputs)
  
  nSeq = length(input)
  input = toupper(input)
  
  tmp = GetPatternScores2(input[1], metric = metric) #placeholder
  nvar = nrow(tmp) #number of variables 
  
  temp.fun <- function(k){
    GetPatternScores2(k, interpolate_distance = interpolate_distance)$Score
  }
  res =  mclapply(input, FUN = temp.fun)
  m = t(as.matrix(bind_rows(res)))
  resFactor = tmp$Factor
  resPosition = tmp$Position
  
  plot_args = data.frame(factor = Factors, position = Positions)
  plots = vector(mode = 'list', length = nrow(plot_args))
  
  for(i in 1:nrow(plot_args)){
    Position = plot_args$position[i]
    Factor = plot_args$factor[i]
    x = m[,(resPosition == Position & as.character(resFactor) == Factor)] #subset msa by specific sequence position.
    spl = split(x, rep(gr_names, grSizes))
    spl = lapply(spl, melt)
    df = bind_rows(spl, .id = 'Group')
    df = df[!is.na(df$value),]
    
    if(!is.na(GroupFactor[1])){df$Group = factor(df$Group, levels = GroupFactor)} #order groups for plotting
    
    if(is.na(GroupCols[1])){
    plots[[i]] = ggplot(df, aes(x = Group, y = value, fill = Group)) + 
      geom_hline(yintercept = 0, size = 1) + 
      geom_violin() + 
      geom_boxplot(alpha = 0, color = 'grey45', outlier.shape = NA) + 
      ylab(Factor) + 
      theme_prism() + 
      labs(title = paste0('Position: ', Position)) + 
      theme(legend.position = 'none')
    }else{
      plots[[i]] = ggplot(df, aes(x = Group, y = value, fill = Group)) + 
        geom_hline(yintercept = 0, size = 1) + 
        geom_violin() + 
        geom_boxplot(alpha = 0, color = 'grey45', outlier.shape = NA) + 
        ylab(Factor) + 
        theme_prism() + 
        labs(title = paste0('Position: ', Position)) + 
        scale_fill_manual(values = GroupCols)
        theme(legend.position = 'none')
    }
     # coord_flip() #+ 
      #stat_compare_means(label.x.npc = 'right', label.y.npc = 'center') 
  }
  
  return(plots)
}
GetPatternPValues <- function(pattern_score_data, aa_matrix, ALPHA = 0.10){
  require(mutoss)
  nSeqs = apply(aa_matrix, 2,FUN = function(x){
    return(sum(!is.na(x)))
  })
  
  pattern_score_data$nSeqs = nSeqs[pattern_score_data$Position]
  
  pvals = unlist(apply(pattern_score_data, 1, FUN = function(x){
    y = pnorm(q = as.numeric(x[3]), mean = 0, sd = 1/sqrt(as.numeric(x[6])))
    return(2*min(c(y,1-y)))
  }))
  pattern_score_data$p = pvals
  pattern_score_data$padj =  BY(pattern_score_data$p, silent = TRUE, alpha = ALPHA)$adjPValues
  pattern_score_data$nSeqs <- NULL
  return(pattern_score_data)
}
#Revised function June 2024 
AverageSeqScores <- function(input, metric = 'ensemble'){
  pmat = readMatrix2(metric)
  subMat = makeSubstitutionMatrix(pmat)
  pmat = apply(pmat, 2, FUN = function(x){
    x = (x - mean(x)) / sd(x)
    return(x)
  })
  
  input = toupper(input)
  aln = strsplit(input, split = '')
  scores = mclapply(aln, FUN = function(m){
    idx = match(m, rownames(pmat))
    out = t(pmat[idx,]) #by residue 
    if(nrow(out)==1) return(as.data.frame(out))
    av = rowSums(out)/ncol(out)
    return(av)
  })
scores = as.data.frame(bind_rows(scores))
 
rownames(scores) = make.names(input, unique = TRUE)
  return(scores)
}
#Function to transform pattern matrix into a euclidian distance similarity matrix 
makeSubstitutionMatrix <- function(pmat){
  #Create substitution matrix compatible with msa package 
  #based on pattern matrix. 
  d = as.matrix(dist(pmat))
  reqNames <- c("A", "R", "N", "D", "C", "Q", "E", "G", 
                "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", 
                "Y", "V", "B", "Z", "X", "*")
  missing = reqNames[which(reqNames %in% rownames(d) == FALSE)]
  
  d2 = matrix(data = NA, nrow = length(reqNames), ncol = length(reqNames))
  rownames(d2) = reqNames; colnames(d2) = reqNames
  
  idx = match(rownames(d), rownames(d2))
  d2[idx,idx ] = d
  return(d2)
}

