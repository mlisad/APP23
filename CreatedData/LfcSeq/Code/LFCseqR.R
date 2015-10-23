


LFCseq <- 
  function(countsTable, condsAB, norm.method="npseq", nloc=50, pares=T){
    
    #Filter countsTable
    filt <- filter(countsTable)
    countsTable <- filt$countsTable
    keep <- filt$keep
    
    #Estimate sequencing depth 
    if(norm.method=="rpm"){
      depth <- estimateDepth(countsTable, norm="rpm")      
    }else if(norm.method=="deseq"){
      depth <- estimateDepth(countsTable, norm="deseq")      
    }else if(norm.method=="npseq"){
      depth <- npSeq.estimate.depth(countsTable)   
    }
    counts.norm <- t(t(countsTable) / depth)
    
    #Estimate probability of non-DE for each feature
    res <- LFCprob(counts.norm, condsAB=condsAB, nloc=nloc, pares=pares)$prob
    
    prob <- rep(1, length(keep))
    prob[keep] <- res
    
    list(pval=prob)
  }
