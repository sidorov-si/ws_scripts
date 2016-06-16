max2 = function(array) { 
	n = length(array)
	sort(array,partial=n-1)[n-1]
}

dynTT = function(csvfile) { 
	genes = read.csv(file=csvfile,head=T)
	e_m = as.matrix(genes[1:length(genes[,1]),2:length(genes[1,])])
	rownames(e_m) = genes[1:length(genes[,1]),1]
	e_mt = e_m
	cutoff = 100
	count=1
	repeat { 
		maximum2 = apply(e_mt, MARGIN=1,FUN=max2)
		proper.Index = which(maximum2>cutoff) # notice how elegant this shit iz!!! 
		e_mt = as.matrix(e_mt[proper.Index,])
		if ( length(e_mt[,1])>6000 ) {
			cutoff = cutoff + 50
			count = count + 1 
		} else { 
			break
		}
	}
	print(paste("Converged at iteration:",count,",kept genes:",length(e_mt[,1]),",expression cutoff:",cutoff))
}

dynTT2 = function(csv,ngenes) {
	## use on collapsed gene CSV file. New version using less unnessesary butt wiggles. 
	## ngenes should be like 6000, 10000, or 15000 
	
	tt <- read.csv(file=csv,head=T,row.names=1)
	
	cutoff = 10 ## initial guess for NON log-scale 
	count=1
	maxx2 <- apply(tt, MARGIN=1,FUN=max2)
	repeat { 
		ttf <- tt[which(maxx2>cutoff),]
		if (nrow(ttf) > ngenes) {
			cutoff = cutoff + 10
			count = count + 1 
		} else { 
			break
		}
	}
	print(paste("Converged at iteration:",count,",kept genes:",nrow(ttf),",expression cutoff:",cutoff))
	ttfl   <- log2(ttf) 
	ngk  <- ngenes/1000 
	ngk  <- paste(ngk,"K",sep="")
	tag1 <- paste("nolog_cut",ngk,"genes.csv",sep="_")
	tag2 <- paste("log2_cut",ngk,"genes.csv",sep="_")
	out1 <- gsub("collapsed.csv",tag1,csv)
	out2 <- gsub("collapsed.csv",tag2,csv)
	write.table(ttf,quote=F,file=out1,row.names=T,col.names=NA,sep=",")
	write.table(ttfl,quote=F,file=out2,row.names=T,col.names=NA,sep=",")
}

collapse1 = function(pat) { 
## the fucntion implies a certain annotation format, in which expression table 
## comes after 4 columns: probe, gene symbol, Entrez gene ID, and RefSeq trans
## cript ID. Header has all the sample names, too. 
## default pattern is "_final_expression.csv" 
	lf=list.files(pattern=pat)
	print("Processing the following files")
	print(lf) 
	print("******************************")
	for (i in 1:length(lf)) {
	  datANN = read.csv(file=lf[i],head=T)
		e_m = as.matrix(datANN[1:length(datANN[,1]),5:length(datANN[1,])])
		myEntrez = as.vector(datANN[1:length(datANN[,1]),3])
		myProbes = as.vector(datANN[1:length(datANN[,1]),1])
		rownames(e_m) = myProbes
		datCOL = collapseRows(e_m,myEntrez,myProbes,method="MaxMean")
		outfile = gsub(pat,"_entrez_collapsed.csv",lf[i]);
		write.table(datCOL$datETcollapsed,file=outfile,quote=F,row.names=TRUE,col.names=NA,sep=",")
		print(paste("File processed:",lf[i],"Written:",outfile))
	}
}

cluster2 = function(csvfile) {
	genes <- read.csv(file=csvfile,head=T)
	e_m <- as.matrix(genes[1:length(genes[,1]),2:length(genes[1,])])
	rownames(e_m) <- genes[1:length(genes[,1]),1]
	e_mt <- t(e_m)
	names(e_mt) <- rownames(e_m)
	powers <- c(c(1:10), seq(from = 12, to=20, by=2))
	sft <- pickSoftThreshold(e_mt, powerVector = powers, verbose = 5)
	kk <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
	print("Array of R squared values:")
	print(kk)

	for (i in 1:length(kk)) {
		if (kk[i]>0.8) {
			print(paste("R squared first surpasses the threshold (0.8) at power = ",powers[i]))
			break
		}
	}
	print(paste("Clustering is performed using beta (power) of ",powers[i]," and correlation coefficient of ",kk[i]))
	TAG <- gsub("_filtered.csv","",csvfile)
	net_h <- blockwiseModules(e_mt, power = powers[i], minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25,
maxBlockSize=10000,deepSplit=4,numericLabels = TRUE, pamRespectsDendro = FALSE,
saveTOMs = TRUE,saveTOMFileBase = paste(TAG,"_h",sep=""),networkType = "signed hybrid",verbose = 3)
	table(net_h$colors)
	write.table(cbind(rownames(e_m),net_h$colors),file=paste(TAG,"_hybrid_modules.txt",sep=""),quote=F,col.names=F,row.names=F,sep="\t")
	rownames(net_h$MEs) <- colnames(e_m)
	write.table(t(net_h$MEs),file=paste(TAG,"_eigengenes.txt",sep=""),quote=F,col.names=NA,row.names=T,sep="\t")
	print("Clustering is complete. See files *modules.txt for the list of modules and  *eigengenes.txt for averaged module expression")
	print("-------------------------------------------------")
}

listThresholds = function (csvfile) { 
	library(WGCNA)
	
	genes <- read.csv(file=csvfile,head=T)
	e_m <- as.matrix(genes[1:length(genes[,1]),2:length(genes[1,])])
	rownames(e_m) <- genes[1:length(genes[,1]),1]
	e_mt <- t(e_m)
	names(e_mt) <- rownames(e_m)
	
	powers <- c(1:20)
	sft <- pickSoftThreshold(e_mt, powerVector = powers, verbose = 5)
	kk <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
	cat(sprintf("Array of R squared values:\n"))
	cat(sprintf("%s",paste(kk,sep=" ")))
}

wgcna_preprocessed_parasens = function(csvfile,ngenes=6000) {
  ## v 2.0 parameter sensitivity testing
  ## meant to generate clusters for 5 different cutoffs and powers 8-25 (no skipping). 
  
  library (WGCNA)
  stringsAsFactors=F
  enableWGCNAThreads()
	cat(sprintf("Processing file %s, top %d most expressed genes will be considered\n",csvfile,ngenes))

	tt <- read.csv(csvfile,row.names=1)
	ttf <- tt[grep("NONE",tt$Symbol,invert=T),]
	ttf <- ttf[grep("NONE",ttf$Entrez_ID,invert=T),]
	ttc <- collapseRows(ttf[3:ncol(ttf)],ttf$Entrez_ID,rownames(ttf))
	
	cat(sprintf("Initial expression matrix dimensions:\n"))
	cat(sprintf("%d\n",dim(tt)))
	cat(sprintf("Defined probe-only expression matrix dimensions:\n"))
	cat(sprintf("%d\n",dim(ttf)))
	cat(sprintf("Collapsed expression matrix dimensions:\n"))
	cat(sprintf("%d\n",dim(ttc$datETcollapsed)))
	
	exp <- as.data.frame(ttc$datETcollapsed)
	exp$max2 <- apply(exp,1,FUN=max2)
	exp <- exp[order(exp$max2,decreasing=T),]
	exp <- exp[1:ngenes,]
	exp$max2 <- NULL 
	cat(sprintf("Truncated collapsed expression matrix dimensions:\n"))
	cat(sprintf("%d\n",dim(exp)))
	
	texp <- t(exp)
	pwrs <- c(5:25)
	cuts <- c(0.05,0.10,0.15,0.20,0.25)
	cat(sprintf("======================================================================\n"))
	sft <- pickSoftThreshold(texp, powerVector = pwrs, verbose = 0)
	print.data.frame(sft$fitIndices)
	cat(sprintf("======================================================================\n"))
	TAG <- gsub("_preprocessed.csv","",csvfile)
	
for (i in cuts) {
  for (j in pwrs) {
    cat(sprintf("Performing clustering for mergeCutHeight %0.2f and power %d\n",i,j))
    
	  net <- blockwiseModules(texp, power = j, minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = i,
	                          maxBlockSize=10000,deepSplit=4,numericLabels=T,pamRespectsDendro=F,networkType = "signed hybrid",verbose = 0)
	  table(net$colors)
	  write.table(cbind(rownames(exp),net$colors),file=paste(TAG,"_power_",j,"_mch_",i,"_hybrid_modules.txt",sep=""),quote=F,col.names=F,row.names=F,sep="\t")
	  write.table(t(net$MEs),file=paste("power_",j,"_mch_",i,"_eigengenes.txt",sep=""),quote=F,col.names=NA,row.names=T,sep="\t")
	  cat(sprintf("----------------------------------------------------------------------\n"))
  }
}
	
	cat(sprintf("CLUSTERING IS COMPLETE. See files *modules.txt for the list of modules and  *eigengenes.txt for averaged module expression.\n"))
}

vioplot1 = function(file1) {
## violin plots for the file that has a header (states) and columns of different height (expression values). 
	library(vioplot)
	tb <- read.table(file1,header=TRUE,fill=TRUE)
	nstate <- length(tb[1,])
	l<-list()
	
	for (i in 1:nstate) {
		k <- paste("E",i,sep="")
		x[i] <- tb$k[!is.na(tb$k)]
	}
	vioplot(x,col="gold")
}

plotpng = function(pat) { 
## the fucntion to make plots from eigengenes  
## pat should be "_eigengene.txt"  
## have to load library(gplots) before you can use this function 
	lf=list.files(pattern=pat)
	print("Processing the following files")
	print(lf) 
	print("******************************")
	for (i in 1:length(lf)) {
	    outfile = gsub(pat,".png",lf[i]);
		data <- read.table (lf[i],sep="\t",header=T,row.names=1)
		m<-as.matrix(data)
		mt<-t(m)
		png(filename=outfile,height=1200,width=2400)
		heatmap.2(mt,dendrogram="none",Rowv=NA,Colv=NA,key=F,density.info="none",trace="none",col=colorRampPalette(c("blue","white","red")),lmat=rbind(c(0,0,0),c(0,1,2),c(0,3,4)),lhei=c(0.05,4,1),lwid=c(0.05,4,4))
		dev.off()
		print(paste("File processed:",lf[i],"Written:",outfile))
	}
}

violabc = function (tag) {
	library(vioplot)
	Nabc <- paste(tag,".ABC.exp",sep="")
	Nabo <- paste(tag,".ABo.exp",sep="")
	Nbco <- paste(tag,".BCo.exp",sep="")
	Naco <- paste(tag,".ACo.exp",sep="")
	Nao  <- paste(tag,".Ao.exp",sep="")
	Nbo  <- paste(tag,".Bo.exp",sep="")
	Nco  <- paste(tag,".Co.exp",sep="")

	tabc <- read.table(Nabc,header=F)
	tabo <- read.table(Nabo,header=F)
	tbco <- read.table(Nbco,header=F)
	taco <- read.table(Naco,header=F)
	tao  <- read.table(Nao,header=F)
	tbo  <- read.table(Nbo,header=F)
	tco  <- read.table(Nco,header=F)

	Xabc <- log2(tabc$V3)
	Xabo <- log2(tabo$V3)
	Xbco <- log2(tbco$V3)
	Xaco <- log2(taco$V3)
	Xao  <- log2(tao$V3)
	Xbo  <- log2(tbo$V3)
	Xco  <- log2(tco$V3)
	
	Xabcf <- Xabc[!is.infinite(Xabc)]
	Xacof <- Xaco[!is.infinite(Xaco)]
	Xbcof <- Xbco[!is.infinite(Xbco)]
	Xabof <- Xabo[!is.infinite(Xabo)]
	Xaof  <- Xao[!is.infinite(Xao)]
	Xbof  <- Xbo[!is.infinite(Xbo)]
	Xcof  <- Xco[!is.infinite(Xco)]

	vioplot(Xabcf,Xabof,Xbcof,Xacof,Xaof,Xbof,Xcof,names=c("ABC","ABo","BCo","ACo","Ao","Bo","Co"),col="gold",ylim=c(-1,20))
}

box3log = function (tag,column) {
	## build boxplot of 3 things - with log2 transformation of the value 
	ll <- c("3p_R2","5p_R2","fl_R2")
	par(mfrow=c(1,3)) 
	for (i in ll) {
		Nabc <- paste(i,".",tag,".ABC.exp",sep="")
		Nabo <- paste(i,".",tag,".ABo.exp",sep="")
		Nbco <- paste(i,".",tag,".BCo.exp",sep="")
		Naco <- paste(i,".",tag,".ACo.exp",sep="")
		Nao  <- paste(i,".",tag,".Ao.exp",sep="")
		Nbo  <- paste(i,".",tag,".Bo.exp",sep="")
		Nco  <- paste(i,".",tag,".Co.exp",sep="")

		tabc <- read.table(Nabc,header=F)
		tabo <- read.table(Nabo,header=F)
		tbco <- read.table(Nbco,header=F)
		taco <- read.table(Naco,header=F)
		tao  <- read.table(Nao,header=F)
		tbo  <- read.table(Nbo,header=F)
		tco  <- read.table(Nco,header=F)

		Xabc <- log2(abs(tabc[,column]))
		Xabo <- log2(abs(tabo[,column]))
		Xbco <- log2(abs(tbco[,column]))
		Xaco <- log2((taco[,column]))
		Xao  <- log2(abs(tao[,column]))
		Xbo  <- log2(abs(tbo[,column]))
		Xco  <- log2(abs(tco[,column]))
	
		Xabcf <- Xabc[!is.infinite(Xabc)]
		Xacof <- Xaco[!is.infinite(Xaco)]
		Xbcof <- Xbco[!is.infinite(Xbco)]
		Xabof <- Xabo[!is.infinite(Xabo)]
		Xaof  <- Xao[!is.infinite(Xao)]
		Xbof  <- Xbo[!is.infinite(Xbo)]
		Xcof  <- Xco[!is.infinite(Xco)]

		boxplot(Xabcf,Xabof,Xbcof,Xacof,Xaof,Xbof,Xcof,names=c("ABC","ABo","BCo","ACo","Ao","Bo","Co"))
	}
}

box3nolog = function (tag,column) {
	## build boxplot of 3 things 
	ll <- c("3p_R2","5p_R2","fl_R2")
	par(mfrow=c(1,3)) 
	for (i in ll) {
		Nabc <- paste(i,".",tag,".ABC.exp",sep="")
		Nabo <- paste(i,".",tag,".ABo.exp",sep="")
		Nbco <- paste(i,".",tag,".BCo.exp",sep="")
		Naco <- paste(i,".",tag,".ACo.exp",sep="")
		Nao  <- paste(i,".",tag,".Ao.exp",sep="")
		Nbo  <- paste(i,".",tag,".Bo.exp",sep="")
		Nco  <- paste(i,".",tag,".Co.exp",sep="")

		tabc <- read.table(Nabc,header=F)
		tabo <- read.table(Nabo,header=F)
		tbco <- read.table(Nbco,header=F)
		taco <- read.table(Naco,header=F)
		tao  <- read.table(Nao,header=F)
		tbo  <- read.table(Nbo,header=F)
		tco  <- read.table(Nco,header=F)

		Xabc <- abs(tabc[,column])
		Xabo <- abs(tabo[,column])
		Xbco <- abs(tbco[,column])
		Xaco <- abs(taco[,column])
		Xao  <- abs(tao[,column])
		Xbo  <- abs(tbo[,column])
		Xco  <- abs(tco[,column])
	
		Xabcf <- Xabc[!is.infinite(Xabc)]
		Xacof <- Xaco[!is.infinite(Xaco)]
		Xbcof <- Xbco[!is.infinite(Xbco)]
		Xabof <- Xabo[!is.infinite(Xabo)]
		Xaof  <- Xao[!is.infinite(Xao)]
		Xbof  <- Xbo[!is.infinite(Xbo)]
		Xcof  <- Xco[!is.infinite(Xco)]

		boxplot(Xabcf,Xabof,Xbcof,Xacof,Xaof,Xbof,Xcof,names=c("ABC","ABo","BCo","ACo","Ao","Bo","Co"))
	}
}

box3pval = function (tag,column) {
	## build boxplot of 3 things - with log2 transformation of the value 
	ll <- c("3p_R2","5p_R2","fl_R2")
	par(mfrow=c(1,3)) 
	for (i in ll) {
		Nabc <- paste(i,".",tag,".ABC.exp",sep="")
		Nabo <- paste(i,".",tag,".ABo.exp",sep="")
		Nbco <- paste(i,".",tag,".BCo.exp",sep="")
		Naco <- paste(i,".",tag,".ACo.exp",sep="")
		Nao  <- paste(i,".",tag,".Ao.exp",sep="")
		Nbo  <- paste(i,".",tag,".Bo.exp",sep="")
		Nco  <- paste(i,".",tag,".Co.exp",sep="")

		tabc <- read.table(Nabc,header=F)
		tabo <- read.table(Nabo,header=F)
		tbco <- read.table(Nbco,header=F)
		taco <- read.table(Naco,header=F)
		tao  <- read.table(Nao,header=F)
		tbo  <- read.table(Nbo,header=F)
		tco  <- read.table(Nco,header=F)

		Xabc <- -log10(abs(tabc[,column]))
		Xabo <- -log10(abs(tabo[,column]))
		Xbco <- -log10(abs(tbco[,column]))
		Xaco <- -log10((taco[,column]))
		Xao  <- -log10(abs(tao[,column]))
		Xbo  <- -log10(abs(tbo[,column]))
		Xco  <- -log10(abs(tco[,column]))
	
		Xabcf <- Xabc[!is.infinite(Xabc)]
		Xacof <- Xaco[!is.infinite(Xaco)]
		Xbcof <- Xbco[!is.infinite(Xbco)]
		Xabof <- Xabo[!is.infinite(Xabo)]
		Xaof  <- Xao[!is.infinite(Xao)]
		Xbof  <- Xbo[!is.infinite(Xbo)]
		Xcof  <- Xco[!is.infinite(Xco)]

		boxplot(Xabcf,Xabof,Xbcof,Xacof,Xaof,Xbof,Xcof,names=c("ABC","ABo","BCo","ACo","Ao","Bo","Co"))
	}
}

run_deseq = function (pat) {
	## this assumes certain format of names which is not generic. Adjust as needed. 
	## This is purpusefully set up to only work with two conditions. 
	library(DESeq2)
	dir1 <- getwd()
	sampleFiles <- list.files(pattern=pat)
	## this is the best I came up with? 
	sampleCondition <- sub(paste(".*(",bquote(.(pat)),"_[0-9]h).*",sep=""),"\\1",sampleFiles,perl=T)
	##sampleCondition <- sub(".*(3p_R2_[0-9]h).*","\\1",sampleFiles,perl=T)
	sampleTable <- data.frame(sampleName = sampleFiles,fileName = sampleFiles,condition = sampleCondition)
	print("Processing the following files:")
	print(sampleFiles)
	print("Processing the following conditions:")
	print(sampleCondition)
	print("----------------------------------------------------")
	ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,directory = dir1, design= ~ condition)
	ddsHTSeq <- DESeq(ddsHTSeq)
	res <- results(ddsHTSeq)
	resOrdered <- res[order(res$padj),]
	head(resOrdered)
	kk  <- unique(sampleCondition)
	tag <- paste(kk[1],"vs",kk[2],sep="_")
	output <- paste(tag,".DESeq.out",sep="")
	## now let's filter of the annoying extra lines from HTseq output 
	resOrderedFiltered<-resOrdered[grep("__*", rownames(resOrdered),invert=T), ]
	write.table(resOrderedFiltered,file=output,quote=F)
}

run_deseq2_matrix = function (cond1,cond2,condfile,expmatrix) {
	## this one assumes that you have condfile of format <SAMPLE>\t<CONDITION>
	## and expression matrix that has gene IDs in the first column and integer counts in the rest
	## header of the column should contail sample names
  
  stringsAsFactors=F
	library(DESeq2)
  condTable <- read.table(condfile,header=F,col.names=c("sampleID","condition"),stringsAsFactors=F)
	## check.names=F allows to avoid adding "X" to the column names. 
	exp <- read.table(expmatrix,header=T,row.names=1,check.names=F,sep="\t")
	cond12 <- condTable[condTable$condition==cond1 | condTable$condition==cond2,]
	## this is neat
	exp12 <- exp[colnames(exp)%in%cond12$sampleID]
	print("Processing the following samples and conditions:")
	print(cond12)
	print("--------------------------------------------------")
	ddsMatrix<-DESeqDataSetFromMatrix(countData=as.matrix(exp12),colData=cond12,design = ~ condition)
	ddsMatrix <- DESeq(ddsMatrix)
	res <- results(ddsMatrix)
	resOrdered <- res[order(res$padj),]

	tag <- paste(cond1,"vs",cond2,sep="_")
	output <- paste(tag,".DESeq.out",sep="")
	write.table(resOrdered,file=output,quote=F)
}

run_limma_set = function (condfile,expmatrix,log.transformed=F,excluded,cut=T) {
  cond <- read.table(condfile,header=T,row.names=1,col.names=c("Sample_ID","Condition"),sep="\t")	
  excl <- read.table(excluded,header=T,row.names=1,col.names=c("Sample_ID","Condition"),sep="\t")	
  cond <- cond[!row.names(cond) %in% row.names(excl),,drop=F]
  vals <- as.vector(unique(cond$Condition))
  
  for (i in 1:length(vals)) {
    for (j in 1:length(vals)) {
      if (i>j) {
        cat(sprintf("Processing the following pair of conditions: %s vs. %s\n",vals[i],vals[j]))
        run_limma_matrix(vals[i],vals[j],condfile,expmatrix,log.transformed,excluded,cut)
      }
    }
  }
}

run_limma_matrix = function (cond1,cond2,condfile,expmatrix,log.transformed=F,excluded,cut=T) {
	
	library(limma)
	stringsAsFactors=F
	cond <- read.table(condfile,header=T,row.names=1,col.names=c("Sample_ID","Condition"),sep="\t")	
	excl <- read.table(excluded,header=T,row.names=1,col.names=c("Sample_ID","Condition"),sep="\t")	
	exp  <- read.table(expmatrix,header=T,row.names=1,check.names=F,sep="\t")
	k0 <- ncol(exp)
	k1 <- nrow(cond)
	k2 <- nrow(excl)
	cat(sprintf("Initial table: %d columns, condfile: %d samples, to exclude: %d samples\n",k0,k1,k2))
	exp  <- exp[,!colnames(exp) %in% row.names(excl)]
	cond <- cond[!row.names(cond) %in% row.names(excl),,drop=F]
	k0 <- ncol(exp)
	k1 <- nrow(cond)
	k2 <- k0-k1
	cat(sprintf("After removal of the excluded samples: %d columns in expmatrix, %d rows in condfile, %d annotation rows\n",k0,k1,k2))

	exp$max2 <- apply(exp[,colnames(exp) %in% row.names(cond)],1,FUN=max2)
	exp <- exp[order(exp$max2,decreasing=T),]
	
	if (cut) {
	  exp <- exp[1:6000,]
	} 
	
	cond12 <- cond[cond$Condition==cond1 | cond$Condition==cond2,,drop=F]
	exp12 <- exp[colnames(exp) %in% row.names(cond12)]  
	ann <- exp[,!colnames(exp) %in% c(row.names(cond),"max2"),drop=F]
	ann <- ann[order(rownames(ann)),]

	gt <- factor(cond12$Condition, levels=c(cond1,cond2))
	design <- model.matrix(~gt)
	
	cat(sprintf("Processing the following samples and conditions:\n"))
	print(cond12)
	cat(sprintf("--------------------------------------------------\n"))
	
	if (!log.transformed) {
		exp12 <- log2(exp12)
	}
		
	fit <- lmFit(exp12,design)
	fit <- eBayes(fit)
	tt  <- topTable(fit,number=100000)
	tt  <- tt[order(rownames(tt)),]
	tt  <- cbind(tt,ann)
	tt  <- tt[order(tt$adj.P.Val),]
	
	outfile <- paste(cond1,"_vs_",cond2,".limma.out",sep="")
	write.table(tt,outfile,quote=F,sep="\t",row.names=F)
}

run_deseq2_hm = function (cond1,cond2,condfile,expmatrix) {
  ## this one assumes that you have condfile of format <SAMPLE>\t<CONDITION>\t<DONOR>     
  ## version for expression matrix format.        

  
}


run_dss = function (cond1,cond2,condfile,expmatrix) {

	## pairwise differential expression analysis using DSS package. 
	## it is similar to DESeq2 but supposed to have higher sensitivity. 

	library(DSS)
	condTable <- read.table(condfile,header=F,col.names=c("sampleID","condition"),stringsAsFactors=F)
	## check.names=F allows to avoid adding "X" to the column names. 
	exp <- read.table(expmatrix,header=T,row.names=1,check.names=F)
	cond12 <- condTable[condTable$condition==cond1 | condTable$condition==cond2,]
	## this is neat
	exp12 <- exp[colnames(exp)%in%cond12$sampleID]
	print("Package DSS, variant with expression matrix. Processing the following samples and conditions:")
	print(cond12)
	print("--------------------------------------------------")

	## this is package-specific shit now. 
	exp12m <-as.matrix(exp12)
	colnames(exp12m) <- NULL
	seqData <- newSeqCountSet(exp12m,cond12[,2])
	seqData <- estNormFactors(seqData)
	seqData <- estDispersion(seqData)
	res     <- waldTest(seqData,cond1,cond2)

	## write it out, as usual
	tag <- paste(cond1,"vs",cond2,sep="_")
	output <- paste(tag,".DSS.out",sep="")
	write.table(res,file=output,quote=F)
	}
	
	
plotDispShrink <- function(dds) {
	# pick 40 equally spaced genes along the base mean
	bins <- 10^seq(from = 0, to = 5, length = 20)
	pickone <- function(x) {
		if (sum(x) == 0)
		return(NULL)
		if (sum(x) == 1)
		return(which(x))
		sample(which(x), 1)
	}
	up <- sapply(seq_along(bins[-1]), function(i) pickone(mcols(dds)$dispGeneEst > 
	1e-04 & !mcols(dds)$dispOutlier & mcols(dds)$dispGeneEst > mcols(dds)$dispFit &
	mcols(dds)$baseMean > bins[i] & mcols(dds)$baseMean < bins[i + 1]))
	down <- sapply(seq_along(bins[-1]), function(i) pickone(mcols(dds)$dispGeneEst >
	1e-04 & !mcols(dds)$dispOutlier & mcols(dds)$dispGeneEst < mcols(dds)$dispFit &
	mcols(dds)$baseMean > bins[i] & mcols(dds)$baseMean < bins[i + 1]))
	# pick 5 outliers
	bins <- 10^seq(from = 1, to = 4, length = 6)
	outliers <- do.call(c, lapply(seq_along(bins[-1]), function(i) pickone(mcols(dds)$dispGeneEst/mcols(dds)$dispMAP >
	2 & mcols(dds)$dispOutlier & mcols(dds)$baseMean > bins[i] & mcols(dds)$baseMean <
	bins[i + 1])))
	s <- c(up, down, outliers)
	s <- s[!is.na(s)]
	with(mcols(dds[s, ]), plot(baseMean, dispGeneEst, log = "xy", pch = 16,
	xlab = "mean of normalized counts", ylab = "dispersion estimate", yaxt = "n",
	ylim = c(0.001, 100)))
	axis(2, at = 10^(-3:2), label = 10^(-3:2))
	xs <- 10^(-20:50/10)
	lines(xs, dispersionFunction(dds)(xs), col = "red", lwd = 2)
	with(mcols(dds[s, ][!mcols(dds[s, ])$dispOutlier, ]), arrows(baseMean, dispGeneEst,baseMean, dispersion, length = 0.075, col = "dodgerblue", lwd = 2))
	with(mcols(dds[s, ][mcols(dds[s, ])$dispOutlier, ]), segments(baseMean, dispGeneEst, baseMean, dispMAP, col = "dodgerblue", lwd = 2, lty = 3))
	with(mcols(dds[s, ][mcols(dds[s, ])$dispOutlier, ]), points(baseMean, dispersion, cex = 2, col = "dodgerblue", lwd = 2))
	legend("topright", c("MLE", "prior mean", "MAP"), pch = 20, col = c("black","red", "dodgerblue"), bg = "white")
}


doubleNormalFit = function(csv) { 
	library(mixtools)

	## in contrast to dyntt1 and dyntt2, this whole ordeal is happening in log2-transformed space. Deal with it. 
	## initial collapsed files, however, are all NOT log2-transformed. 
	
	tt<-read.csv(csv,header=T,row.names=1)
	cat(sprintf("Processing file with %d rows, %d columns..\n",nrow(tt),ncol(tt)))
	if (max(tt,na.rm=T) > 20) {
		cat(sprintf("The experiment has been identified as NON-log transformed\n"))
		x <- unlist(tt)
		minx <- min(x)
		minr <- max(minx,0)
		xf <- x[x>minr+1]
		xl <- log2(xf) 
		xx <- xl[!is.na(xl)]
		cat(sprintf("Minimum from the experiment: %f, defacto used minimum: %f\n",minx,minr))
		cat(sprintf("Overall values in distribution: %d, after value filtering: %d, after NA filtering: %d\n",length(x),length(xf),length(xx)))
	} else {
		cat(sprintf("The experiment has been identified as LOG2 transformed\n"))
		x <- unlist(tt)
		minx <- min(x)
		minr <- max(minx,0)
		xf <- x[x>minr+0.2]
		xx <- xl[!is.na(xl)]
		cat(sprintf("Minimum from the experiment: %f, defacto used minimum: %f\n",minx,minr))
		cat(sprintf("Overall values in distribution: %d, after value filtering: %d, after NA filtering: %d\n",length(x),length(xf),length(xx)))
	}
	d <- density(xx)
	kk <- normalmixEM(xx,verb=T)
	pngout <- paste(csv,".png",sep="")
	png(filename=pngout,height=900,width=900)
	plot(kk,main2=csv,which=2)
	lines(d)
	dev.off()
	
	## finding the x at which they overlap boils down to 
	l1 <- kk$lambda[1]
	l2 <- kk$lambda[2]
	m1 <- kk$mu[1]
	m2 <- kk$mu[2]
	s1 <- kk$sigma[1]
	s2 <- kk$sigma[2]
	cat(sprintf("The resulting values were found for sum of normal distributions:\n"))
	cat(sprintf("Lambdas: %f,%f, Sigmas: %f,%f, Mu: %f,%f\n",l1,l2,s1,s2,m1,m2))
	a <- s1^2-s2^2
	b <- 2*m1*s2^2-2*m2*s1^2
	c <- s1^2*m2^2-s2^2*m1^2-2*s1^2*s2^2*log((l2*s1)/(l1*s2))
	D <- b^2-4*a*c
	x1 <- (-b+sqrt(D))/(2*a)
	x2 <- (-b-sqrt(D))/(2*a)
	
	ttl <- log2(tt)
	mm <- apply(as.matrix(ttl),MARGIN=1,FUN=max2)
	if (s1<s2) {
		cutx=max(x1,x2)
	} else { 
		cutx=min(x1,x2)
	}
	cat(sprintf("Normal distributions cross at %f or %f; under provided conditions, we choose %f\n",x1,x2,cutx))
	
	if (cutx<8 & !is.na(cutx)) {
		ttxl <- ttl[which(mm>cutx),]
	} else { 
		## in case 2-normal fit failed or produced unreasonable values, let's just take top 10k genes by our max2 criterion
		cat(sprintf("WARNING: The two-normal distribution fit did not converge, using top 10,000 genes."))
		cutoff = 3 ## initial guess for NON log-scale 
		count=1
		maxx2 <- apply(ttl, MARGIN=1,FUN=max2)
		repeat { 
			ttxl <- ttl[which(maxx2>cutoff),]
			if (nrow(ttxl) > 10000) {
				cutoff = cutoff + 0.1
				count = count + 1 
			} else { 
				break
			}
		}
		print(paste("Converged at iteration:",count,",kept genes:",nrow(ttxl),",LOG2 expression cutoff:",cutoff))
	}
	
	## seriously, this works. R is fucking insane. 
	ttnl <- 2^ttxl
	##ttxl <- ttl[rowMeans(ttl)>max(x1,x2),]
	cat(sprintf("Estimated cutoff is %f, number of genes in the resulting filtered file is %d\n",cutx,nrow(ttxl))) 
	
	out1 <- gsub("collapsed.csv","nolog_2normal.csv",csv)
	out2 <- gsub("collapsed.csv","log2_2normal.csv",csv)
	write.table(ttnl,quote=F,file=out1,row.names=T,col.names=NA,sep=",")
	write.table(ttxl,quote=F,file=out2,row.names=T,col.names=NA,sep=",")
	}

	
plotNormFit = function (q) {
	## given a vector of values, should return you a fit of your shitty data. 
	library(MASS)
	qf<-fitdistr(q,"normal")
	mean<-as.numeric(qf$estimate[1])
	sd<-as.numeric(qf$estimate[2])
	low=mean-5*sd
	high=mean+5*sd
	x<-seq(low,high,length=1000)
	y<-dnorm(x,mean=mean, sd=sd)
	dd<-density(q)
	mm<-max(max(y),max(dd$y))
	plot(dd,xlim=c(low,high),ylim=c(0,mm))
	par(new=T)
	plot(x,y, type="l", lwd=1,col="red",xlim=c(low,high),ylim=c(0,mm))
	}
	
make_svg_heatmaps = function (filename) {
	##library(ggplot2)
	##library(reshape)
  ##rm(list=ls(all=TRUE)) ## crucial here!
	tag <- gsub("_eigengenes2.txt","",filename)
	tt <- read.table(filename,header=T,row.names=NULL)
	nsmp <- ncol(tt)-1
	nmod <- nrow(tt)-1
	tt.m <- melt(tt)
	
	names(tt.m) <- c("Module","Sample","value")
	## the following are very important
	tt.m$Module <- as.factor(tt.m$Module)
	sorted_labels <- paste(sort(as.integer(levels(tt.m$Module))))
	tt.m$Module <- factor(tt.m$Module, levels = sorted_labels)
	tt.m$Sample <- with(tt.m, factor(Sample,levels = rev(sort(unique(Sample)))))
	
	## let's evaluate SVG dimensions. We want constant width of 10 in. 
	title_length <- max(nchar(names(tt)))
	svg_height  <- round((10/27.1)*nsmp,digits=3)
	svg_width <- round((10/27.1)*(0.25*title_length+nmod+1),digits=3)
	cat (sprintf("Longest sample title is %d characters, there are %d samples and %d modules, plot width is %f, height is %f\n",title_length,nsmp,nmod,svg_width,svg_height))
	
	base_size <- 9
	for (i in 0:nmod) {
	  ## special type of coloring for easier heatmap reading. Can be adjusted. 
		p <- ggplot(tt.m, aes(Module,Sample)) + geom_tile(aes(fill = value), colour = "black",size=0.5) + scale_fill_gradientn(colours=c("blue","blue","white","red","red"),space="Lab")
		p2 <- p + theme_grey(base_size=base_size)+labs(x="",y="")+scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand =c(0,0))+theme(legend.position ="none",axis.text.y=element_text(size=base_size*1.2,colour="black"),axis.text.x=element_text(size=base_size*1.2,colour="black"))
		xmin <- i+0.5
		xmax <- i+1.5
		ymin <- 0.5
		ymax <- nsmp+0.5
		rect <- data.frame(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
		p3 <- p2 + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),color="black",fill="grey",alpha=0, size=3, inherit.aes = F)
		svgname <- paste(tag,"_module_",i,".svg",sep="")
		##cat(sprintf("Processing module number %d, saving file %s\n",i,svgname))
		ggsave(file=svgname, plot=p3, width=svg_width, height=svg_height,limitsize=F)
		}
}

make_svg_heatmaps2 = function (filename) {
  ## With non-renamed/non-normalized eigengene tables.
  ##library(ggplot2)
  ##library(reshape)
  ##rm(list=ls(all=TRUE)) ## crucial here!
  tag <- gsub("_eigengenes.txt","",filename)
  tt <- read.table(filename,header=T,row.names=1)
  row.names(tt) <- as.numeric(gsub("ME","",row.names(tt)))
  tts <- tt[order(as.numeric(row.names(tt))),]
  tts2 <- tts[,order(colnames(tts))]
  ttn <- as.data.frame(t(apply(tts2,1,FUN=exvector)))
  ttn$Module <- row.names(ttn)
  
  nsmp <- ncol(tt)
  nmod <- nrow(tt)-1 ## number of non-zero modules
  tt.m <- melt(ttn)
  
  names(tt.m) <- c("Module","Sample","value")
  ## the following are very important
  tt.m$Module <- as.factor(tt.m$Module)
  sorted_labels <- paste(sort(as.integer(levels(tt.m$Module))))
  tt.m$Module <- factor(tt.m$Module, levels = sorted_labels)
  tt.m$Sample <- with(tt.m, factor(Sample,levels = rev(sort(unique(Sample)))))
  
  ## let's evaluate SVG dimensions. We want constant width of 10 in. 
  title_length <- max(nchar(names(tt)))
  svg_height  <- round((10/27.1)*nsmp,digits=3)
  svg_width <- round((10/27.1)*(0.25*title_length+nmod+1),digits=3)
  cat (sprintf("Longest sample title is %d characters, there are %d samples and %d modules, plot width is %f, height is %f\n",title_length,nsmp,nmod,svg_width,svg_height))
  
  base_size <- 9
  for (i in 0:nmod) {
    ## special type of coloring for easier heatmap reading. Can be adjusted. 
    p <- ggplot(tt.m, aes(Module,Sample)) + geom_tile(aes(fill = value), colour = "black",size=0.5) + scale_fill_gradientn(colours=c("blue","blue","white","red","red"),space="Lab")
    p2 <- p + theme_grey(base_size=base_size)+labs(x="",y="")+scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand =c(0,0))+theme(legend.position ="none",axis.text.y=element_text(size=base_size*1.2,colour="black"),axis.text.x=element_text(size=base_size*1.2,colour="black"))
    xmin <- i+0.5
    xmax <- i+1.5
    ymin <- 0.5
    ymax <- nsmp+0.5
    rect <- data.frame(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    p3 <- p2 + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),color="black",fill="grey",alpha=0, size=3, inherit.aes = F)
    svgname <- paste(tag,"_module_",i,".svg",sep="")
    ##cat(sprintf("Processing module number %d, saving file %s\n",i,svgname))
    ggsave(file=svgname, plot=p3, width=svg_width, height=svg_height,limitsize=F)
  }
}

eigenorm = function (val,min,max) {
  if (min==0 | max==0) stop ("Line minimum or a maximum equal to 0")
  if (val <= 0) { 
    return(round(val/abs(min),digits=3))
  } else {
    return(round(val/max,digits=3))
  }
}

exvector = function (v) {
  min=min(v)
  max=max(v)
  for (i in 1:length(v)) {
    if (v[i] <= 0) {
      v[i]=round(v[i]/abs(min),digits=3)
    } else { 
      v[i]=round(v[i]/max,digits=3)
    }  
  }
  return(v)
}