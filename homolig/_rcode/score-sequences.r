#Using matrix factorizations of AADist, describe any group of AA sequences.
ScorePatterns2 <- function(input, metric = 'nmf10', error = 'se', pattern_names = NA){
  res = GetPatternScores(input, metric = metric, pattern_names = pattern_names)
  plotres = PlotPatternScores(res, error = error)
    return(plotres)
}
ScorePatterns3 <- function(input, metric = 'nmf10', error = 'se', pattern_names = NA,
                           interpolate_distance = NA, plotSig = TRUE){
  res = GetPatternScores2(input, metric = metric, pattern_names = pattern_names, 
                          interpolate_distance = interpolate_distance)
  plotres = PlotPatternScores(res, error = error, plotSig = plotSig)
  return(plotres)
}
ScoreLogos <- function(input, metric = 'nmf10'){
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
process_logo <- function(sequence_info, metric = 'nmf10'){
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
  load('/media/aag7319/WDRed/ZZZ_Homolig/coordinate_objects/2023-06-26MFMats.rda')
  idx = which(names(allmat) == metric)
  return(allmat[[idx]])
}
PlotPatternScores <- function(pattern_score_data, error = 'se', plotSig = TRUE, ALPHA = 0.10){
  require(ggplot2)
  require(ggprism)
  require(scales) #pretty_breaks()
  out_df = pattern_score_data
  
  sevals = pattern_score_data$se
  sevals[is.na(sevals)] = 0 #for purposes of computing plot limit
  if(max(abs(pattern_score_data$Score) + sevals ) >= 1.5){
    YLIM = max(abs(pattern_score_data$Score) + pattern_score_data$se) + 0.1
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
GetPatternScores <- function(input, metric = 'nmf10', ALPHA = 0.1, pattern_names = NA){
  require(msa)
  pmat = readMatrix2(metric)
  subMat = makeSubstitutionMatrix(pmat)
  pmat = apply(pmat, 2, FUN = function(x){
    x = (x - mean(x)) / sd(x)
    return(x)
  })
  if(!is.na(pattern_names[1]) & length(pattern_names) == ncol(pmat)){
    colnames(pmat) = pattern_names
  }

  nSeq = length(input)
  input = toupper(input)
  if(sd(nchar(input)) == 0){ #If all sequences are identical length
    aln = strsplit(input, split = '')
  }else{
    tmp =msa(input, method = 'ClustalW', cluster = 'nj', #substitutionMatrix = 'id', 
             type = 'protein', gapOpening=10, gapExtension=0.2, 
             substitutionMatrix = subMat)
    aln = as.data.frame(tmp@unmasked)
    aln = strsplit(aln$x, split = '')
  }
  
  L = length(aln[[1]])
  m = matrix(data = NA, nrow = nSeq, ncol = L)
  for(r in 1:nSeq){
    m[r,] = aln[[r]] 
  }
  
  m[m=='-'] = NA
  
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
  
  out_se = apply(m, 2, FUN = function(x){
    idx = match(x, rownames(pmat))
    tmp = pmat[idx,]
    se = apply(tmp, 2, FUN = function(y){
      return(sd(y, na.rm = TRUE)/sqrt(sum(!is.na(x))))
    })
    return(se) #colSums(tmp, na.rm = TRUE) / nrow(tmp)
  })
  
  outsd_df = reshape2::melt(out_sd); colnames(outsd_df) = c('Factor', 'Position', 'SD')
  outse_df = reshape2::melt(out_se); colnames(outse_df) = c('Factor', 'Position', 'SE')
  out_df = reshape2::melt(out); colnames(out_df) = c('Factor', 'Position', 'Score')
  out_df$sd = outsd_df$SD
  out_df$se = outse_df$SE
  
  out_df$Factor = factor(out_df$Factor)
  out_df = GetPatternPValues(out_df, m, ALPHA = ALPHA)
  return(out_df)
}
GetPatternScores2 <- function(input, metric = 'nmf10', ALPHA = 0.1, pattern_names = NA, interpolate_distance = NA){
  #Interpolate sequence distance to avoid downstream length-dependent artifacts in msa. 
  #10 October 2023. 
  lentab = table(nchar(input))
  if(is.na(interpolate_distance)){
    #if not specified, select the median of 
    interpolate_distance = round(median(nchar(input)))
    #interpolate_distance = as.numeric(names(lentab))[which.max(lentab)[1]] #selects mode
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
ComparePatternScores <- function(inputs, metric = 'nmf10', ALPHA = 0.1, pattern_names = NA){
  #Compare pattern scores between multiple groups.
  #inputs is a list of sequences by group, with list names corresponding to group names. 
  #11 September 2023 AAG 
  require(msa)
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
   input = unlist(inputs)
 
  nSeq = length(input)
  input = toupper(input)
  if(sd(nchar(input)) == 0){ #If all sequences are identical length
    aln = strsplit(input, split = '')
  }else{
    tmp =msa(input, method = 'ClustalW', cluster = 'nj', #substitutionMatrix = 'id', 
             type = 'protein', gapOpening=10, gapExtension=0.2, 
             substitutionMatrix = subMat)
    aln = as.data.frame(tmp@unmasked)
    aln = strsplit(aln$x, split = '')
  }
  
  L = length(aln[[1]])
  m = matrix(data = NA, nrow = nSeq, ncol = L)
  for(r in 1:nSeq){
    m[r,] = aln[[r]] 
  }
  
  m[m=='-'] = NA
  
  out = apply(m,2,FUN= function(x){
    idx = match(x, rownames(pmat))
    tmp = pmat[idx,] #pattern scores for all composite residues 
    tmp[is.na(tmp)] = 0 #Unsure how to handle this. Missing value problem. Important to consider!
    tmp = as.data.frame.matrix(tmp)
    spl = split(tmp, rep(gr_names, grSizes))
    
    spl = lapply(spl, melt)
    
    df = bind_rows(spl, .id = 'Group')

    df = df[!is.na(df$value),]
    
    kwp =  tapply(df, df$variable, FUN = function(x){kruskal.test(value ~Group, data = x)$p.value})
    #kwp = BY(kwp, alpha = ALPHA, silent = TRUE)$adjPValues #Doing FDR at each Position across all Patterns.
    return(kwp)
  })
 
  
 res= melt(out)
  colnames(res) = c('Factor', 'Position', 'kw_pval')
  res$padj = BY(res$kw_pval, alpha = ALPHA, silent = TRUE)$adjPValues
  res = res[order(res$kw_pval, decreasing = FALSE),]
  rownames(res) = NULL
  if(any(res$padj < ALPHA)){
    return(res[res$padj < ALPHA,])
  }else{
    return(res)
  }

}
ComparePatternScores2 <- function(inputs, metric = 'nmf10', ALPHA = 0.1, pattern_names = NA, interpolate_distance = NA){
  #Compare pattern scores between multiple groups.
  #Interpolate sequence distance to avoid downstream length-dependent artifacts in msa. 
  #AAG 10 October 2023. 
  
  require(reshape2)
  require(dplyr)
  
  #CODE NOT WRITTEN. 
  input = toupper(input)
  aln = strsplit(input, split = '')
  nchars = nchar(input)
  
  spl = split(aln, nchars) #split by sequence length. 
  
  
  out = apply(m,2,FUN= function(x){
    idx = match(x, rownames(pmat))
    tmp = pmat[idx,] #pattern scores for all composite residues 
    tmp[is.na(tmp)] = 0 #Unsure how to handle this. Missing value problem. Important to consider!
    tmp = as.data.frame.matrix(tmp)
    spl = split(tmp, rep(gr_names, grSizes))
    
    spl = lapply(spl, melt)
    
    df = bind_rows(spl, .id = 'Group')
    
    df = df[!is.na(df$value),]
    
    kwp =  tapply(df, df$variable, FUN = function(x){kruskal.test(value ~Group, data = x)$p.value})
    #kwp = BY(kwp, alpha = ALPHA, silent = TRUE)$adjPValues #Doing FDR at each Position across all Patterns.
    return(kwp)
  })
  
  res= melt(out)
  colnames(res) = c('Factor', 'Position', 'kw_pval')
  res$padj = BY(res$kw_pval, alpha = ALPHA, silent = TRUE)$adjPValues
  res = res[order(res$kw_pval, decreasing = FALSE),]
  rownames(res) = NULL
  if(any(res$padj < ALPHA)){
    return(res[res$padj < ALPHA,])
  }else{
    return(res)
  }
  
}
PlotPatternScoreComparison <- function(inputs, metric = 'nmf10', Factors, Positions, pattern_names = NA){
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
  input = unlist(inputs)
  
  nSeq = length(input)
  input = toupper(input)
  if(sd(nchar(input)) == 0){ #If all sequences are identical length
    aln = strsplit(input, split = '')
  }else{
    tmp =msa(input, method = 'ClustalW', cluster = 'nj', #substitutionMatrix = 'id', 
             type = 'protein', gapOpening=10, gapExtension=0.2, 
             substitutionMatrix = subMat)
    aln = as.data.frame(tmp@unmasked)
    aln = strsplit(aln$x, split = '')
  }
  
  L = length(aln[[1]])
  m = matrix(data = NA, nrow = nSeq, ncol = L)
  for(r in 1:nSeq){
    m[r,] = aln[[r]] 
  }
  
  m[m=='-'] = NA
  
  plot_args = data.frame(factor = Factors, position = Positions)
  plots = vector(mode = 'list', length = nrow(plot_args))
  for(i in 1:nrow(plot_args)){
    Position = plot_args$position[i]
    Factor = plot_args$factor[i]
    x = m[,Position] #subset msa by specific sequence position.
    idx = match(x, rownames(pmat))
    tmp = pmat[idx,] #pattern scores for all composite residues 
    tmp[is.na(tmp)] = 0 #Unsure how to handle this. Missing value problem. Important to consider!
    
    tmp = tmp[, Factor] #subset to factor of interest
    spl = split(tmp, rep(gr_names, grSizes))
    spl = lapply(spl, melt)
    df = bind_rows(spl, .id = 'Group')
    df = df[!is.na(df$value),]
    
    plots[[i]] = ggplot(df, aes(x = Group, y = value)) + 
      geom_boxplot() + 
      ylab(Factor) + 
      theme_prism() + 
      coord_flip() + 
      stat_compare_means(label.x.npc = 'right', label.y.npc = 'center') 
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

makeSeqs <- function(L, N, from = NA, probs = NA){
 if(all(is.na(from))){
   from = c("A", "R", "N", "D" ,"C" ,"Q" ,"E" ,"G" ,
            "H" ,"I" ,"L" ,"K" ,"M" ,"F", "P" ,"S" ,
            "T" ,"W" ,"Y" ,"V")
 }
  
  seqs = vector(mode = 'character', length = N)
  tempfun = function(input){
    if(all(is.na(probs))) return(paste0(sample(from, size = L, replace = TRUE), collapse = ''))
    return(paste0(sample(from, size = L, replace = TRUE, prob = probs), collapse = ''))
  }
  
 seqs = unlist( lapply(seqs,tempfun))
 return(seqs)
}

replaceNthPosition <- function(seqs, position, replaceOptions){
  substr(seqs, position, position) <- sample(replaceOptions, length(seqs), replace = TRUE)
  return(seqs)
}
NMF10_NAMES = c('Accessibility', 'Pattern_2', 'Pattern_3', 'Pattern_4',
                'Hydrophobicity', 'Charge', 'Secondary Structure', 'a-Helix Propensity',
                'Composition', 'Steric Hinderance')
#TEST DATA-------------------------------------------
# input  = read.csv('/media/aag7319/IronWolf/ZZZ_H2/test_data/random25k.csv')$sequence[2e3:3e3]
# 
# 
# makeSeqs <- function(l,N){
#   #l = length sequences
#   #N = number sequences
#   mat = readMatrix2('pca02')
#   aas = rownames(mat)
#   seqs = vector(mode = 'character', length = N)
#   seqs = unlist(lapply(seqs, FUN = function(x){
#     return(paste0(sample(aas, l, replace = TRUE), collapse = ''))}))
#   return(seqs)
# }
# replaceNthPosition <- function(seqs, position, replaceOptions){
#   substr(seqs, position, position) <- sample(replaceOptions, length(seqs), replace = TRUE)
#   return(seqs)
# }
# 
# input <- makeSeqs(15, 1e2)
# ns = c(1,3,5,8,12,20,20,12,8,5,3,1)
# lens = c(10:20)
# for(i in 1:length(lens)){
#   temp = paste0('CASS', makeSeqs(lens[i], ns[i]*5), 'QYF')
#   input = c(input,temp)
# }
# #
# dat = GetPatternScores(input, pattern_names = NA)
#  tmp = ScorePatterns3(input, pattern_names = NA)
#
# inputA <- makeSeqs(15, 1e2)
# inputB <- makeSeqs(15,1e2)
#
# inputA <- replaceNthPosition(inputA, 5, c('Q', 'Y', 'F'))
# inputB <- replaceNthPosition(inputB, 10, c('S', 'T'))
#
# tmp = ScorePatterns2(inputA, pattern_names  = nmf10_names)
# tmp = ScorePatterns2(inputB, pattern_names  = nmf10_names)
#
# inputs = list(A = inputA, B = inputB)
# comp =  ComparePatternScores(inputs, pattern_names = nmf10_names)
#

#
# comp_plots = PlotPatternScoreComparison(inputs, Factors = head(comp$Factor),
#                                         Positions = head(comp$Position), pattern_names = nmf10_names)

#msaClustalW(input)