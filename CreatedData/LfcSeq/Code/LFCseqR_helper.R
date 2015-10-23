


#-------------------------------------------------------------#
# Compute  log2 fold change
#-------------------------------------------------------------#

LFCmd <-
function (dat = dat, selec = c(1:nrow(dat))) {
  pares <- as.matrix(combn(ncol(dat), 2))

  if (NCOL(pares) > 30) {  
    sub30 <- sample(1:NCOL(pares), size = 30, replace = FALSE)
    pares <- pares[,sub30]
  }
  
  mm <- NULL
  dd <- NULL
  for (i in 1:ncol(pares)) {
    a <- dat[selec,pares[1,i]]
    b <- dat[selec,pares[2,i]]
    mm <- cbind(mm, log(a/b, 2))
  }
  
  list(M = mm)
}

#-------------------------------------------------------------#
# Compute  log2 fold change for random pairs in one condition
#-------------------------------------------------------------#

LFCmd2 <-
  function (dat = dat) {
    
    n <- ncol(dat)
    n2 <- round(n/2)
    
    if(n==2){
      mm <- as.matrix(log(dat[, 1]/dat[, 2], 2))
      dd <- as.matrix(abs(dat[, 1]-dat[, 2]))
    }else{  
        if(n <= 7){
            if(n2==(n/2)){
                firstGroup <- as.matrix(combn(n, n2))
                firstGroup <- firstGroup[, 1:(ncol(firstGroup)/2)]
                firtGroup <- as.matrix(firstGroup)
            }else{
                firstGroup <- as.matrix(combn(n, n2))
            }
        }else{
            firstGroup <- replicate(20, sample(1:n, n2))
        }
        mm <- NULL
        dd <- NULL
        for (i in 1:ncol(firstGroup)){
            a <- apply(as.matrix(dat[, firstGroup[, i]]), 1, mean)
            b <- apply(as.matrix(dat[, setdiff(1:n, firstGroup[, i])]), 1, mean)
            mm <- cbind(mm, log(a/b, 2))
        }
    
    }
    
    list(M = mm)
  }


#-------------------------------------------
#-  Estimate probability of non-DE for each feature
#---------------------------------------------


LFCprob <- 
  function(counts.norm, condsAB, nloc = 50, pares=T){
    
    ABNames <- unique(condsAB)
    factIndex <- list(which(condsAB == ABNames[1]), which(condsAB == ABNames[2]))
    
    if(pares==F){
      MD1 <- LFCmd(dat = counts.norm[, factIndex[[1]]])  
      MD2 <- LFCmd(dat = counts.norm[, factIndex[[2]]])
    }else{
      MD1 <- LFCmd2(dat = counts.norm[, factIndex[[1]]])  
      MD2 <- LFCmd2(dat = counts.norm[, factIndex[[2]]])
    }
    
    resum1 <- rowMeans(counts.norm[, factIndex[[1]]])
    resum2 <- rowMeans(counts.norm[, factIndex[[2]]])   
    MDs <- LFCmd(dat = cbind(resum1, resum2))  
    rowM <- rowMeans(counts.norm)
     
    nloc <- min(nloc, nrow(counts.norm))
    prob <- NULL
    
    for(i in 1:nrow(MDs$M)){  
      dis <- abs(rowM[i] - rowM)
      locGen <- sort(dis, decreasing = F, index = T)$ix[1 : nloc]
      
      noiM <- abs(c(as.numeric(MD1$M[locGen, ]), as.numeric(MD2$M[locGen, ])))
      sigM <- abs(MDs$M[i, ])
      
      num <- sum(sigM < noiM)
      tot <- length(noiM)
      prob <- c(prob, num/tot)
    }
    
    
    list(prob = prob )
    
  }
  
  
  
#----------------------------------
#-  Estimate Depth:  rpm and deseq
#-----------------------------------

estimateDepth <-
function(data, norm){
    if(norm == "rpm"){
        total <- colSums(as.matrix(data))
        depth <- total / 10^6
    }else if(norm == "deseq"){
        logmeans <- rowMeans( log(as.matrix(data)) )
        depth <- apply(data, 2, function(cnts)
            exp(median(log(cnts)-logmeans)))
        
    }
    depth
}

#----------------------------------
#-  Estimate Depth:  npSeq
#-----------------------------------

npSeq.estimate.depth <- function(n, iter=5)
{
  #iter <- 5
  cmeans <- colSums(n) / sum(n)
  for (i in 1:iter)
  {
    n0 <- rowSums(n) %*% t(cmeans)
    prop <- rowSums((n - n0)^2 / (n0 + 1e-08))
    qs <- quantile(prop, c(0.25, 0.75))
    keep <- (prop >= qs[1]) & (prop <= qs[2])
    cmeans <- colMeans(n[keep, ])
    cmeans <- cmeans / sum(cmeans)
  }
  depth <- cmeans / mean(cmeans)
  
  return(depth)
}

#-------------------------------------------
#-   Filter low expressed features
#------------------------------------------------


filter <- 
function(countsTable, ct.sum = 5, ct.mean = 0.5){

  if(is.null(rownames(countsTable))){
    rownames(countsTable) <- 1 : nrow(countsTable)
  }
  
  keep <- (rowMeans(countsTable) > ct.mean) & (rowSums(countsTable) > ct.sum) & 
    apply(countsTable, 1, function(a)all(a>=0))
  
  cat(length(keep) - sum(keep), "gene has been filtered because they contains
        too small number of read across the experiments.", fill = T)
  
  countsTable <- countsTable[keep, ]
  countsTable[countsTable==0] <- 0.5
  
  list(countsTable=countsTable, keep=keep)
}
