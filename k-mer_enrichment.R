


#' Get TF motif enrichment for given gkm-PCs 
#'
#' Applies the minimum hypergeometric test to all PC*motif combinations, for highly-weighted k-mers and lowly-weighted k mers.
#' 
#' @param rotated A numeric matrix containing the k-mer loadings for however many PCs you want analyzed (rows are k-mers, columns are PCs).
#' @param binaryKMerMatchesToTFs A matrix of k-mers by TF motifs, where each row is a k-mer and each column is a TF motif. Each entry in the matrix is 1 if the k-mer matches the TF's motif, and 0 otherwise.
#' @param n_max The maximum number of k-mers to consider for each minimum hypergeometric test. Defaults to 3000. Should be no more than about 10% of the number of k-mers. 
#' @param verbose Print verbose output, describing progression.
#' @return a data.frame containing the following columns: Motif_ID (the motif ID from binaryKMerMatchesToTFs), PC (the PC from rotated), p (minimum hyper geometric ln(p-values)), k (the number of top k-mers that yielded maximal enrichment), log2OR (the log2 odds ratio (observed/expected) of k-mers for the point of maximal enrichment), i (a unique integer for each PC-motif-direction combination).
#' @keywords 
#' @export
#' @examples
#' kmerMat = inputKMerFreqs(sprintf("kMerFiles/%s.freq.gz",sampleDesc$id), IDs = sampleDesc$id)
#' myPCA = doKMerPCA(kmerMat, nPCs = "jackstraw")
#' treatmentPCs = findDistinguishingPCs(myPCA$x[,1:myPCA$nPCs], sampleDesc[c("id","treated")])
#' tfEnrichmentsPBM = getKMerTFEnrichment(myPCA$rotation[,1:myPCA$nPCs], cisbp$binaryPBMZScores);
#' tfEnrichmentsPBM = tfEnrichmentsPBM[order(tfEnrichmentsPBM$p),]
#' tfEnrichmentsPBM = merge(tfEnrichmentsPBM2, cisbp$TFTable[c("Motif_ID","TF_Name")], by="Motif_ID") # add TF names
#' head(tfEnrichmentsPBM[tfEnrichmentsPBM$PC==treatmentPCs$PC[1],], n=20) # show top 20 motifs for the first treatment-distinguishing PC

getKMerTFEnrichment = function(rotated, binaryKMerMatchesToTFs,n_max=3000, verbose=F){
  z=1;
  rotated = rotated[row.names(binaryKMerMatchesToTFs),]; #filter out and sort rows so that they're in the same order
  tfKMerEnrichments = data.frame(Motif_ID = NA, PC=NA, p = NA, k = NA, logOR = NA, direction = NA, i = 1:(ncol(rotated)*ncol(binaryKMerMatchesToTFs)*2), stringsAsFactors = F );
  for(pci in 1:ncol(rotated)){
    curOrder = order(rotated[,pci]); #increasing order
    if (verbose){
      message(sprintf("PC = %i/%i",pci,ncol(rotated)));
    }
    for (tfi in 1:ncol(binaryKMerMatchesToTFs)){
      if (tfi %% 20==0 && verbose){
        message(sprintf("\tTF = %i/%i",tfi,ncol(binaryKMerMatchesToTFs)));
      }
      tfKMerEnrichments$Motif_ID[z:(z+1)] = colnames(binaryKMerMatchesToTFs)[tfi];
      tfKMerEnrichments$PC[z:(z+1)] = colnames(rotated)[pci];
      tfKMerEnrichments$direction[z:(z+1)] = c("low","high")
      curTest = minHG(binaryKMerMatchesToTFs[curOrder,tfi], n_max=n_max); #check for enrichment among low PC weights
      curTestRev = minHG(rev(binaryKMerMatchesToTFs[curOrder,tfi]), n_max=n_max); #check for enrichment among high PC weights
      
      #makeEnrichmentGraph(binaryKMerMatchesToTFs[curOrder,tfi], n_max=3000)
      #makeEnrichmentGraph(rev(binaryKMerMatchesToTFs[curOrder,tfi]), n_max=3000)
      tfKMerEnrichments$p[z] = curTest$lnMinP;
      tfKMerEnrichments$p[z+1] = curTestRev$lnMinP;
      tfKMerEnrichments$k[z] = curTest$k;
      tfKMerEnrichments$k[z+1] = curTestRev$k;
      tfKMerEnrichments$log2OR[z] = curTest$log2OR;
      tfKMerEnrichments$log2OR[z+1] = curTestRev$log2OR;
      z=z+2;
    }
  }
  return(tfKMerEnrichments);
}


#' Do a minimum hyper-geometric test
#'
#' Applies the minimum hypergeometric test to given data, returning the maximal enrichment achieved.
#' 
#' @param x A sorted vector of binary values, where 1 is a "hit" (e.g. cognate k-mer) and 0 is a "miss" (e.g. non-cognate k-mer).
#' @param n_max The maximum number of k-mers to consider for each minimum hypergeometric test. Defaults to testing all possible entries in the vector for maximal enrichment.
#' @return a list containing the lnMinP (natural log(min(p))), k (the number of vector entries included for maximal enrichment), log2OR (log2 odds ratio (observed hits/expected hits)), and n_max (as given to method).
#' @keywords 
#' @export
#' @examples
#' kmerMat = inputKMerFreqs(sprintf("kMerFiles/%s.freq.gz",sampleDesc$id), IDs = sampleDesc$id)
#' myPCA = doKMerPCA(kmerMat, nPCs = "jackstraw")
#' treatmentPCs = findDistinguishingPCs(myPCA$x[,1:myPCA$nPCs], sampleDesc[c("id","treated")])
#' testResultsLow = minHG(as.logical(cisbp$binaryPBMZScores[order(pcs$rotation[row.names(cisbp$binaryPBMZScores),treatmentPCs$PC[1]]),"M0312_1.02.PBM"])); # lowly-weighted k-mers
#' testResultsHigh = minHG(as.logical(cisbp$binaryPBMZScores[order(pcs$rotation[row.names(cisbp$binaryPBMZScores),treatmentPCs$PC[1]],decreasing = T),"M0312_1.02.PBM"])); # highly-weighted k-mers

minHG = function(x, n_max = length(x)-1){
  minP = 0;
  mink = 0;
  logOR = 0;
  n=sum(x); # number of black balls in urn (1s in x)
  m=length(x)-n; #number of white balls in urn (0s in x)
  x=1-x; # make TFBSs=0, so now white=1, black=0
  numX = cumsum(x[1:n_max])
  testI = (1:n_max)[x[1:(n_max-1)]!=x[2:(n_max)]]
  for(i in testI){
    curP = phyper(numX[i], m, n, i, log.p=T) # natural log
    #print(curP)
    if (!is.nan(curP) && curP<minP){
      mink=i;
      minP=curP;
      logOR = log2(((i-numX[i])/i)/(n/length(x)));
    }
  }
  return(list(lnMinP = minP, k=mink, log2OR = logOR, n_max=n_max));
}


#' Makes an enrichment graph for a sorted binary vector
#'
#' Applies the minimum hypergeometric test to given data and plots the enrichment of observed over expected across the data, including the ln(P) value (title), log2(observed/expected) (y-axis; e.g. observed/expected cognate k-mers), and element rank (x-axis, e.g. k-mer). The blue line indicates the point of minimum hypergeometric p-value.
#' 
#' @param x A sorted vector of binary values, where 1 is a "hit" (e.g. cognate k-mer) and 0 is a "miss" (e.g. non-cognate k-mer).
#' @param n_max The maximum number of k-mers to consider for each minimum hypergeometric test. Defaults to testing all possible entries in the vector for maximal enrichment.
#' @param sortedBy An alternate x-axis (instead of vector index) to use for the x-axis of the graph. Generally, what the vector x was sorted by.
#' @return a list containing the plot (plot), the data.frame used to make the plot (rawData), and the minHG test (minHGTest).
#' @keywords 
#' @export
#' @examples
#' kmerMat = inputKMerFreqs(sprintf("kMerFiles/%s.freq.gz",sampleDesc$id), IDs = sampleDesc$id)
#' myPCA = doKMerPCA(kmerMat, nPCs = "jackstraw")
#' treatmentPCs = findDistinguishingPCs(myPCA$x[,1:myPCA$nPCs], sampleDesc[c("id","treated")])
#' plot = makeEnrichmentGraph(as.logical(cisbp$binaryPBMZScores[order(pcs$rotation[row.names(cisbp$binaryPBMZScores),treatmentPCs$PC[1]],decreasing = T),"M0312_1.02.PBM"])); # highly weighted k-mers

makeEnrichmentGraph = function(x, n_max = length(x)-1, sortedBy=NULL){
  testData = data.frame(hits = x, rank =1:length(x), cumu=0, expected = (1:length(x)) * mean(x), OR=1);
  testData$cumu = cumsum(testData$hits)
  testData$cumu[1] = testData$cumu[1]+(testData$expected[1]/2); # add a pseudocount so that it doesn't start at -infinity
  testData$OR= testData$cumu/testData$expected
  testData = rbind(data.frame(hits=0,rank=0,cumu=0,expected=0,OR=1), testData);
  testData$log2OR =log2(testData$OR);
  test = minHG(x, n_max = n_max);
  if (!is.null(sortedBy)){
    testData$sortedBy = sortedBy[1:nrow(testData)];
    p = ggplot(testData, aes(x=sortedBy, y=log2OR)) + geom_line() + theme_bw() + geom_hline(yintercept=0,colour="red")+ggtitle(sprintf("ln(p) = %g",test$lnMinP))+geom_vline(xintercept=testData$sortedBy[testData$rank==test$k], colour="blue") +ylab("log2(O/E)")+xlab("variable"); print(p);
  }else{
    p = ggplot(testData, aes(x=rank, y=log2OR)) + geom_line() + theme_bw() + geom_hline(yintercept=0,colour="red")+ggtitle(sprintf("ln(p) = %g",test$lnMinP))+geom_vline(xintercept=test$k, colour="blue") +ylab("log2(O/E)")+xlab("rank"); print(p);
  }
  
  return (list(plot = p, rawData=testData, minHGTest=test));
}


#' Make a graph showing enrichment of a motif within a PC
#'
#' Applies the minimum hypergeometric test to the current data, and makes a plot where the
#' 
#' @param PC A named numeric vector representing the k-mer loadings of a PC
#' @param binaryData A named logical vector representing k-mer cognate (T) or non-cognate (F) status for a motif.
#' @param n_max The maximum number of k-mers to consider for each minimum hypergeometric test. Defaults to 3000. Should be no more than about 10% of the number of k-mers. 
#' @return a list containing the plot (plot), the data.frame used to make the plot (rawData), and the minHG test (minHGTest).
#' @keywords 
#' @export
#' @examples
#' kmerMat = inputKMerFreqs(sprintf("kMerFiles/%s.freq.gz",sampleDesc$id), IDs = sampleDesc$id)
#' myPCA = doKMerPCA(kmerMat, nPCs = "jackstraw")
#' treatmentPCs = findDistinguishingPCs(myPCA$x[,1:myPCA$nPCs], sampleDesc[c("id","treated")])
#' tfEnrichmentsPBM = getKMerTFEnrichment(myPCA$rotation[,1:myPCA$nPCs], cisbp$binaryPBMZScores);
#' tfEnrichmentsPBM = tfEnrichmentsPBM[order(tfEnrichmentsPBM$p),]
#' tfEnrichmentsPBM = merge(tfEnrichmentsPBM2, cisbp$TFTable[c("Motif_ID","TF_Name")], by="Motif_ID") # add TF names
#' p = makeEnrichmentGraphForPC(pcs$rotation[,treatmentPCs$PC[1]],cisbp$binaryPBMZScores[,head(tfEnrichmentsPBM$Motif_ID[tfEnrichmentsPBM$PC==treatmentPCs$PC[1] & tfEnrichmentsPBM$direction=="low"],n=1)]) #top motif for top treatment-distinguishing PC for lowly-weighted k-mers

makeEnrichmentGraphForPC = function(PC, binaryData, n_max=3000, decreasing=F){
  data = merge(as.data.frame(PC),as.data.frame(binaryData), by=c("row.names"))
  data = data[order(data[[2]],decreasing=decreasing),];
  p = makeEnrichmentGraph(data[[3]], n_max=n_max, sortedBy = data[[2]])
  return(p)
}

