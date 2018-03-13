


#' Find 
#'
#' Identifies which PCs distinguish a feature of interest. If the feature of interest contains two classes, performs a single rank sum test for each PC.
#' If the feature of interest has >2 classes, performs a rank sum test to see if the PC separates the class from all others. 
#' In either case, returns a data frame containing the PC, class, AUROC, and rank sum P-value for how well the PC identifies the class.
#' 
#' @param rotated A numeric matrix containing the k-mer loadings for however many PCs you want analyzed (rows are k-mers, columns are PCs).
#' @param binaryKMerMatchesToTFs A matrix of k-mers by TF motifs, where each row is a k-mer and each column is a TF motif. Each entry in the matrix is 1 if the k-mer matches the TF's motif, and 0 otherwise.
#' @param n_max The maximum number of k-mers to consider for each minimum hypergeometric test. Defaults to 3000. Should be no more than about 10% of the number of k-mers. 
#' @return a data.frame containing the following columns: TF (the TF motif), PC (the PC), pLow (minimum hyper geometric p-value for low-weighted k-mers), kLow (the number of top k-mers that yielded maximal enrichment, for low-weighted k-mers), logORLow (the log odds ratio (observed/expected) of k-mers for the point of maximal enrichment among kowly-weighted k-mers), pHigh (as before, for high-scoring k-mers), kHigh (ibid), logORHigh (ibid), i (a unique integer for each PC-motif combination).
#' @keywords 
#' @export
#' @examples
#' kmerMat = inputKMerFreqs(sprintf("kMerFiles/%s.freq.gz",sampleDesc$id), IDs = sampleDesc$id)
#' myPCA = doKMerPCA(kmerMat, nPCs = "jackstraw")
#` treatmentPCs = findDistinguishingPCs(myPCA$x[,1:myPCA$nPCs], sampleDesc[c("id","treated")])

getKMerTFEnrichment = function(rotated, binaryKMerMatchesToTFs,n_max=3000, verbose=F){
  z=1;
  rotated = rotated[row.names(binaryKMerMatchesToTFs),]; #filter out and sort rows so that they're in the same order
  tfKMerEnrichments = data.frame(TF = NA, PC=NA, pLow = NA, kLow = NA, logORLow = NA, pHigh = NA, kHigh = NA, logORHigh = NA, i = 1:(ncol(rotated)*ncol(binaryKMerMatchesToTFs)), stringsAsFactors = F );
  for(pci in 1:ncol(rotated)){
    curOrder = order(rotated[,pci]); #increasing order
    if (verbose){
      message(sprintf("PC = %i",pci));
    }
    for (tfi in 1:ncol(binaryKMerMatchesToTFs)){
      if (tfi %% 20==0 && verbose){
        message(sprintf("\tTF = %i",tfi));
      }
      tfKMerEnrichments$TF[z] = colnames(binaryKMerMatchesToTFs)[tfi];
      tfKMerEnrichments$PC[z] = colnames(rotated)[pci];
      curTest = minHG(binaryKMerMatchesToTFs[curOrder,tfi], n_max=n_max); #check for enrichment among low PC weights
      curTestRev = minHG(rev(binaryKMerMatchesToTFs[curOrder,tfi]), n_max=n_max); #check for enrichment among high PC weights
      
      #makeEnrichmentGraph(binaryKMerMatchesToTFs[curOrder,tfi], n_max=3000)
      #makeEnrichmentGraph(rev(binaryKMerMatchesToTFs[curOrder,tfi]), n_max=3000)
      tfKMerEnrichments$pLow[z] = curTest$minP;
      tfKMerEnrichments$pHigh[z] = curTestRev$minP;
      tfKMerEnrichments$kLow[z] = curTest$k;
      tfKMerEnrichments$kHigh[z] = curTestRev$k;
      tfKMerEnrichments$logORLow[z] = curTest$logOR;
      tfKMerEnrichments$logORHigh[z] = curTestRev$logOR;
      z=z+1;
    }
  }
  return(tfKMerEnrichments);
}


minHG = function(x, n_max = length(x)-1){
  minP = 0;
  mink = 0;
  logOR = 0;
  n=sum(x); # number of black balls in urn (1s in x)
  m=length(x)-n; #number of white balls in urn (0s in x)
  x=1-x; # make TFBSs=0, so now white=1, black=0
  for(i in 1:n_max){
    if (x[i]!=x[i+1]){
      #print(i)
      curP = phyper(sum(x[1:i]), m, n, i, log.p=T) # natural log
      #print(curP)
      if (!is.nan(curP) && curP<minP){
        mink=i;
        minP=curP;
        logOR = log2((sum(1-x[1:i])/i)/(n/length(x)));
      }
    }
  }
  return(list(minP = minP, k=mink, logOR = logOR));
}


makeEnrichmentGraph = function(x, n_max = length(x)-1, sortedBy=NULL){
  testData = data.frame(hits = x, rank =1:length(x), cumu=0, expected = (1:length(x)) * mean(x), OR=1);
  testData$cumu[1] = testData$hits[1]+(testData$expected[1]/2); # add a pseudocount so that it doesn't start at -infinity
  testData$OR[1] = testData$cumu[1]/testData$expected[1]
  for(i in 2:nrow(testData)){
    testData$cumu[i] = testData$cumu[i-1] + testData$hits[i];
    testData$OR[i] = testData$cumu[i]/testData$expected[i];
  }
  testData = rbind(data.frame(hits=0,rank=0,cumu=0,expected=0,OR=1), testData);
  testData$logOR =log2(testData$OR);
  test = minHG(x, n_max = n_max);
  if (!is.null(sortedBy)){
    testData$sortedBy = sortedBy[1:nrow(testData)];
    p = ggplot(testData, aes(x=sortedBy, y=logOR)) + geom_line() + theme_bw() + geom_hline(yintercept=0,colour="red")+ggtitle(sprintf("log(p) = %g",test$minP))+geom_vline(xintercept=testData$sortedBy[testData$rank==test$k], colour="blue") + coord_cartesian(xlim=sort(c(testData$sortedBy[testData$rank==1],testData$sortedBy[testData$rank==n_max])))+ylab("log2(O/E)")+xlab("variable"); print(p);
  }else{
    p = ggplot(testData, aes(x=rank, y=logOR)) + geom_line() + theme_bw() + geom_hline(yintercept=0,colour="red")+ggtitle(sprintf("log(p) = %g",test$minP))+geom_vline(xintercept=test$k, colour="blue") + coord_cartesian(xlim=c(0,n_max))+ylab("log2(O/E)")+xlab("rank"); print(p);
  }
  
  return (list(plot = p, rawData=testData));
}


makeEnrichmentGraphForPC = function(PC, binaryData, nmax=3000, decreasing=F){
  data = merge(PC,binaryData, by=c("row.names"))
  #str(data)
  data = data[order(data[[2]],decreasing=decreasing),];
  p = makeEnrichmentGraph(data[[3]], n_max=nmax, sortedBy = data[[2]])
  return(p)
}

