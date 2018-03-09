inputKMerFreqs = function(fileNames, IDs=fileNames){
	stopifnot(length(fileNames)==length(IDs))
	curData = read.table(fileNames[i],sep="\t",header=T,stringsAsFactors=FALSE,quote="",row.names=NULL);

	allData=matrix(nrow = nrow(curData),ncol = nrow(sampleDesc));
	row.names(allData)=curData$kMer;
	colnames(allData) = IDs;
	i=1;
	for (i in 1:length(fileNames)){
	  if (i %% 50==1){
	    message(sprintf("Inputting i=%i/%i",i,length(fileNames)))
	  }
	  curData = read.table(sprintf(fileNames[i],sep="\t",header=T,stringsAsFactors=FALSE,quote="",row.names=NULL);
	  allData[,i] = curData$Pct_FG;
	}
	return(allData)
}

doKMerPCA = function(x, nPCs="permutation", scale=T, verbose=F){
	if(scale){
		message("Scaling data")
		scaledData = t(scale(t(x),scale=T, center=T));
  }else{
		scaledData=x;
	}
	rm('x');
	message("Doing PCA")
  pcs = prcomp(t(scaledData), scale.=F, center=F);
	if (nPCs=="jackstraw"){
		message("Finding Significant PCs by Jackstraw")
		require(jackstraw)
  	nPCs = permutationPA(scaledData); 
		pcs$nPCs = nPCs$r
		pcs$nPCMethod = "jackstraw";
	}else if (nPCs=="permutation"){
		message("Finding Significant PCs by partial randomization")
		pcs = findSigPCs(pcs, nfolds=1)
		pcs$nPCMethod = "permutation";
	}else if (is.numeric(nPCs)){
		message(sprintf("Using %i PCs",nPCs))
		pcs$nPCs=nPCs;
		pcs$nPCMethod = "fixed";
	}else{
		stop("Unrecognized 'nPCs' parameter.  Must be numeric or one of 'jackstraw' or 'permutation'");
	}
  rm('scaledData')
  
	message("Doing tSNE")
  ##tSNE only top PCs, including PC1
  pcs$tSNE = tsne(pcs$x[,1:pcs$nPCs]);
  pcs$tSNEProj = data.frame(pcs$tSNE);
  names(pcs$tSNEProj) = c('tSNE1','tSNE2');
  pcs$tSNEProj$ID=row.names(pcs$x);
	return(pcs)
}


findDistinguishingPCs = function(rotatedData, classLabels){
  if (is.data.frame(classLabels)){
    classLabels2 = as.list(classLabels[[2]]);
    names(classLabels2)=classLabels[[1]];
    classLabels=classLabels2;
  }
  stopifnot(is.list(classLabels))
  classes = unique(classLabels);
  allResults = data.frame()
  for(c in classes){
    if (length(classes)>2 || classes[1]==c){
      results = data.frame(PC = colnames(rotatedData), AUROC = NA, P=NA,class=c)
      for (i in 1:ncol(rotatedData)){
        curTest = ranksumROC(rotatedData[names(classLabels)[classLabels==c],i],rotatedData[names(classLabels)[classLabels!=c],i])
        results$AUROC[i] = curTest$AUROC
        results$P[i] = curTest$p.value
      }
      allResults = rbind(allResults, results)
    }
  }
  return(allResults)
}

