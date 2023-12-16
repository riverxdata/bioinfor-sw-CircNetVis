deleteCache=function(){
  system("find ./ -cmin +40 -type f -name 'circRNA_pseudoSeq*' -delete")
  system("find ./ -cmin +40 -type f -name 'RNAhybrid_out*' -delete")
  system("find ./ -cmin +40 -type f -name 'miRanda_out*' -delete")
  system("find ./ -cmin +40 -type f -name 'reactomejson*' -delete")
  return(NULL)
}

getCirmi_RNAhybrid = function(circRNA.in){
  suppressMessages(suppressWarnings(require(Biostrings)))
  load("Annotation_data/genes.exon.all.RData")
  load("Annotation_data/Homo_sapiens.GRCh37.75_geneNames.RData")
  #get library
  #genomeFastaFile ="C:\\Nghiavtr\\Workplaces\\Projects\\common_data\\Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz"
  #fasta_genome = readDNAStringSet(genomeFastaFile)
  #chnames = sapply(names(fasta_genome), function(x) unlist(strsplit(x, " "))[1])
  #names(fasta_genome) = chnames
  #saveRDS(fasta_genome, file = "Annotation_data/Homo_sapiens.GRCh37.75.dna.primary_assembly.rds")
  #fasta_genome=readRDS(file = "Annotation_data/Homo_sapiens.GRCh37.75.dna.primary_assembly.rds")
  
  #txFastaFile ="Annotation_data/Homo_sapiens.GRCh37.75.cdna.all.fa"
  #tx.all.fasta = readDNAStringSet(txFastaFile)
  #tx.all.NAME = sapply(names(tx.all.fasta), function(x) unlist(strsplit(x, " "))[1])
  #save(tx.all.fasta,tx.all.NAME, file="Annotation_data/Homo_sapiens.GRCh37.75.cdna.all.fa.RData")
  #load("Annotation_data/Homo_sapiens.GRCh37.75.cdna.all.fa.RData")
  
    cat("\n CircRNAs not in the RNAhybrid databases: ",circRNA.in)
    x=lapply(circRNA.in,function(z) strsplit(z,"__")[[1]])
    x=do.call(rbind,x)
    chrVec=as.character(x[,1])
    startVec=as.integer(x[,2])
    endVec=as.integer(x[,3])

    cat("\n Generate pseudo-sequence of circRNAs")
    readLen=100
    CPUNUM=1
    registerDoParallel(cores=CPUNUM)
    par.pseudoSeq.fa=foreach(i = 1:length(chrVec)) %dopar% {
      chr=chrVec[i]
      start=startVec[i]
      end=endVec[i]
      
      #### process
      #get information of the precusor gene of the circRNA
      pick=which(genes.exon.all$EXONCHROM==chr & genes.exon.all$EXONSTART>=(start) & genes.exon.all$EXONEND<=(end))[1]
      precursor_gene=genes.exon.all$GENEID[pick]
      myexonInfo=genes.exon.all[genes.exon.all$GENEID %in% precursor_gene,]
      #circRNA_name=paste(chr,start,end,precursor_gene,sep = "__")
      circRNA_name=paste(chr,start,end,precursor_gene,sep = "__")
      ### get pseudo-sequence of circRNA
      #res=getCircPseudoSequence_cDNA(myexonInfo=myexonInfo, tx.all.fasta=tx.all.fasta, tx.all.NAME=tx.all.NAME, circRNA_name=circRNA_name, start=start, end=end, num=-1,readLen=readLen,fasta_genome=NULL)
      res=getCircPseudoSequence_cDNA(myexonInfo=myexonInfo, tx.all.fasta=NULL, tx.all.NAME=NULL, circRNA_name=circRNA_name, start=start, end=end, num=-1,readLen=readLen,fasta_genome=NULL)
      return(res)
    }
    
    pseudoSeq.fa = lapply(par.pseudoSeq.fa, function(p) p$fa)
    pseudoSeq.name = lapply(par.pseudoSeq.fa, function(p) p$name)
    
    pseudoSeq.fa = unlist(pseudoSeq.fa)
    pseudoSeq.name = unlist(pseudoSeq.name)
    
    pseudoSeq.fa = DNAStringSet(pseudoSeq.fa)
    names(pseudoSeq.fa) = pseudoSeq.name
     
    ### generate a random ID
    ranID=floor(runif(1, min=1, max=1000000))
    ### export the pseudo-sequence to file
    pseudoSeq_fn=paste("circRNA_pseudoSeq",ranID,".fa",sep="")
    writeXStringSet(pseudoSeq.fa, file=pseudoSeq_fn)
    
    ##################### run in linux
    cat("\n Run RNAhybrid and collect circRNA-miRNA interactions")
    system('chmod 777 Tools/RNAhybrid')
    #cp /lib/x86_64-linux-gnu/libm.so.6 Tools    
    #system('export LD_LIBRARY_PATH=$PWD/Tools/:$LD_LIBRARY_PATH')
    #Sys.chmod(file, "777", use_umask = FALSE)
    
    circRNA_miRNA_out=paste("RNAhybrid_out",ranID,".txt",sep="")
    miRNA_db="Annotation_data/miRNA_mature_HomoSapiens.fa"
    cmd_RNAhybrid=paste("./Tools/RNAhybrid -m 63780 -n 24 -b 1 -c -s 3utr_human -p 0.05 -q ",miRNA_db," -t ",pseudoSeq_fn, " > ",circRNA_miRNA_out,sep="")
    system(cmd_RNAhybrid)
    
    #####################
    ### collect circRNA-miRNA interactions
    cirmi=read.table(circRNA_miRNA_out,sep=":")
    colnames(cirmi)=c("target","target_len","miRNA","miRNA_len","mfe","pval","target_position","target_5to3","target_binding_regions","miRNA_binding_regions","miRNA_3to5")
    cirmi$miRNA=as.character(cirmi$miRNA)
    cirmi=cirmi[order(cirmi$pval,decreasing = FALSE),]
    #dim(cirmi)
    
    ### delete the files
    system(paste("rm ",pseudoSeq_fn," ",circRNA_miRNA_out,sep=""))
    
    ## build a network
    normID=sapply(as.character(cirmi[,1]), function(z) paste(strsplit(z,"__")[[1]][1:3],collapse ="__"))
    names(normID)=NULL
    cirmi$normID=normID
    return(cirmi)
}


getCirmi_miRanda = function(circRNA.in){
  suppressMessages(suppressWarnings(require(Biostrings)))
  load("Annotation_data/genes.exon.all.RData")
  load("Annotation_data/Homo_sapiens.GRCh37.75_geneNames.RData")

  
  cat("\n CircRNAs not in the miranda databases: ",circRNA.in)
  x=lapply(circRNA.in,function(z) strsplit(z,"__")[[1]])
  x=do.call(rbind,x)
  chrVec=as.character(x[,1])
  startVec=as.integer(x[,2])
  endVec=as.integer(x[,3])
  
  cat("\n Generate pseudo-sequence of circRNAs")
  readLen=100
  CPUNUM=1
  registerDoParallel(cores=CPUNUM)
  par.pseudoSeq.fa=foreach(i = 1:length(chrVec)) %dopar% {
    chr=chrVec[i]
    start=startVec[i]
    end=endVec[i]  
    #### process
    #get information of the precusor gene of the circRNA
    pick=which(genes.exon.all$EXONCHROM==chr & genes.exon.all$EXONSTART>=(start) & genes.exon.all$EXONEND<=(end))[1]
    precursor_gene=genes.exon.all$GENEID[pick]
    myexonInfo=genes.exon.all[genes.exon.all$GENEID %in% precursor_gene,]
    #circRNA_name=paste(chr,start,end,precursor_gene,sep = "__")
    circRNA_name=paste(chr,start,end,precursor_gene,sep = "__")
    ### get pseudo-sequence of circRNA
    res=getCircPseudoSequence_cDNA(myexonInfo=myexonInfo, tx.all.fasta=NULL, tx.all.NAME=NULL, circRNA_name=circRNA_name, start=start, end=end, num=-1,readLen=readLen,fasta_genome=NULL)
    return(res)
  }
  pseudoSeq.fa = lapply(par.pseudoSeq.fa, function(p) p$fa)
  pseudoSeq.name = lapply(par.pseudoSeq.fa, function(p) p$name)
  pseudoSeq.fa = unlist(pseudoSeq.fa)
  pseudoSeq.name = unlist(pseudoSeq.name)
  pseudoSeq.fa = DNAStringSet(pseudoSeq.fa)
  names(pseudoSeq.fa) = pseudoSeq.name
  ### generate a random ID
  ranID=floor(runif(1, min=1, max=1000000))
  ### export the pseudo-sequence to file
  pseudoSeq_fn=paste("circRNA_pseudoSeq",ranID,".fa",sep="")
  writeXStringSet(pseudoSeq.fa, file=pseudoSeq_fn)
  
  ##################### run in linux
  cat("\n Run miRanda and collect circRNA-miRNA interactions")
  system('chmod 777 Tools/miranda')
  circRNA_miRNA_out=paste("miRanda_out",ranID,".txt",sep="")
  miRNA_db="Annotation_data/miRNA_mature_HomoSapiens.fa"
  cmd_miranda=paste("./Tools/miranda ",miRNA_db," ",pseudoSeq_fn, " -quiet | grep '>hsa' > ",circRNA_miRNA_out,sep="")
  system(cmd_miranda)

  cat("\n Run miRanda - done ")
  
  #####################
  ### collect circRNA-miRNA interactions
  cirmi_tmp=readLines(circRNA_miRNA_out)
  
  cirmiSum=NULL
  #hits=grep(">h",cirmi_tmp)
  hits=grep(">>h",cirmi_tmp) #sum-up from hits
  if (length(hits)>0){
    x=strsplit(cirmi_tmp[hits],"\t")
    x1=do.call(rbind,x)
    x1[,1]=gsub(">","",x1[,1])
    cirmiSum=rbind(cirmiSum,x1)
  }

  colnames(cirmiSum)=c("miRNA","circRNA","Tot_Score","Tot_Energy","Max_Score","Max_Energy","Strand","miRNA_Len","circRNA_Len","Positions")
  cirmiSum=as.data.frame(cirmiSum,stringsAsFactors=FALSE)
  
  normID=sapply(cirmiSum$circRNA, function(z) paste(strsplit(z,"__")[[1]][1:3],collapse ="__"))
  names(normID)=NULL
  cirmiSum$normID=normID

  cirmiSum$Max_Score=as.double(cirmiSum$Max_Score)
  cirmiSum$Max_Energy=as.double(cirmiSum$Max_Energy)
  cirmiSum$Tot_Score=as.double(cirmiSum$Tot_Score)
  cirmiSum$Tot_Energy=as.double(cirmiSum$Tot_Energy)

  #keep unique  miRNA-normID, update max_score and max_energy
  x=paste(cirmiSum$miRNA,cirmiSum$normID,sep="_")
  y=table(x)
  x1=tapply(cirmiSum$Max_Score,x,max)
  x2=tapply(cirmiSum$Max_Energy,x,max)
  pick=match(names(x1),x)
  cirmiSum2=cirmiSum[pick,]
  cirmiSum2$Max_Score=x1
  cirmiSum2$Max_Energy=x2

  cirmi=cirmiSum2[,c("miRNA","circRNA","Max_Score","Max_Energy","normID")]
  cat("\n miRanda - done")
  return(cirmi)
}


convertReverseComplement<-function(DNAseq){
  DNAarr=unlist(strsplit(DNAseq,""))
  #reverse
  DNAarr=rev(DNAarr)
  #complement
  Aid=which(DNAarr=="A")
  Tid=which(DNAarr=="T")
  Gid=which(DNAarr=="G")
  Cid=which(DNAarr=="C")
  DNAarr[Aid]="T"
  DNAarr[Tid]="A"
  DNAarr[Gid]="C"
  DNAarr[Cid]="G"
  #result
  DNAseqRc=paste(DNAarr,collapse = "")
  return(DNAseqRc) 
}


getCircPseudoSequence_cDNA <-function (myexonInfo, fasta_genome=NULL,tx.all.fasta=NULL, tx.all.NAME=NULL, circRNA_name, start, end, num=-1,readLen){
  #load tx name
  load("Annotation_data/Homo_sapiens.GRCh37.75.cdna.tx.NAME.RData")
  #output
  circRNA.tx.fa = c() #start with a dumpy DNA sequence - serial mode
  circRNA.tx.name = c() #start with a dumpy DNA sequence - serial mode
  #get all exons completely inside the boundary
  pick1=myexonInfo$EXONSTART == start
  pick2=myexonInfo$EXONEND == end
  mytx=intersect(unique(myexonInfo[pick1,]$TXNAME), unique(myexonInfo[pick2,]$TXNAME))
  #find all transcripts and generate circRNA from those transcripts
  mytx=intersect(mytx,tx.all.NAME) #check if the tx in cdna
  txNum=length(mytx)

  if (txNum > 0){ # if existing transcripts containing both EXONSTART and EXONEND
    #atx=sample(mytx,1)
    if (num > 0){ #if fixed the number of selected tx, normally num=1 for the simulator tool
      num=min(num,length(mytx))
      mytx=sample(mytx,num)
    }
    #load data here
    if (is.null(tx.all.fasta)) load("Annotation_data/Homo_sapiens.GRCh37.75.cdna.all.fa.RData")
    
    for (atx in mytx){ #do for every tx
      mytxex=myexonInfo[myexonInfo$TXNAME==atx,]
      mytxex=mytxex[order(mytxex$EXONSTART),]
      mytxex.len = mytxex$LENGTH

      #covert to transcript location
      mytxex.end.all =cumsum(mytxex$LENGTH)
      mytxex.stat.all = mytxex.end.all - mytxex.len +1
      #select the coding area inside the circRNA region
      pick= which(mytxex$EXONSTART >= start & mytxex$EXONEND <= end) #pick by indexes of exons

      startR= mytxex.stat.all[pick]
      endR= mytxex.end.all[pick]
      ## generate pseudo circularRNAs sequences
      # - strand
      if( mytxex$EXONSTRAND[1] == "-"){
        #pick transcript sequence
        transcript = convertReverseComplement(as.character(tx.all.fasta[which(tx.all.NAME == atx)]))
      }else{
        transcript = as.character(tx.all.fasta[which(tx.all.NAME == atx)])
      }
      
      x=substring(transcript,startR,endR)
      txSeq=paste(x,collapse="")

      # new create a pseudo-circularRNA sequence by adding l-1 nucleotides from the end of the last exon to the begin of the first exon
      l_1=readLen -1 -1
      if (l_1 > nchar(txSeq)) l_1=nchar(txSeq) - 1 #if the circRNA length < l_1
      circRNA.seq=paste(substring(txSeq, nchar(txSeq)-l_1),txSeq,sep="")
      circRNA.as_name=paste(atx,"++",sep="") #do not forget to put the "circRNA", otherwises the format of the output file rawcount.txt of FuSeq is broken  
      #add to result
      circRNA.tx.fa=c(circRNA.tx.fa,circRNA.seq)
      circRNA.tx.name = c(circRNA.tx.name, circRNA.as_name)
    }
  }else{ # if not existing the transcript, build circRNA from the exon contig between two ends
    mytx="NONE"
    #get all exons completely inside the boundary
    pick=myexonInfo$EXONSTART >= start & myexonInfo$EXONEND <= end
    myex=myexonInfo[pick,]
    #select the coding area inside the region. 
    #NOTE: this can not always satisfy the canonical splicing conditions (GT-AG), but we do not care about it this momment
    myex=myex[order(myex$EXONSTART),]
    exstart=myex$EXONSTART/1e6
    exend=myex$EXONEND/1e6
    #get clusters of exons
    exCluster=rep(-1,length(exstart))
    for (j in 1:length(exstart)){
      if (exCluster[j]==-1){
        exCluster[j]=j
        pick=(exstart[j]-exstart)*(exstart[j]-exend)<=0 | (exend[j]-exstart)*(exend[j]-exend)<=0
        #assign to the cluster with min index
        x=exCluster[pick]
        x=x[x>0]
        x=min(x)
        exCluster[pick]=x
      }
    }
    #update exCluster
    exCluster_u=exCluster
    repeat{
      exCluster=exCluster[exCluster]
      if (sum(exCluster_u!=exCluster)==0) break()
      exCluster_u=exCluster
    }

    #create exon regions from the exon clusters
    clusterID=sort(unique(exCluster))
    startR=endR=NULL
    for (j in 1:length(clusterID)){
      cID=clusterID[j]
      cEx=myex[exCluster==cID,]
      startR=c(startR,min(cEx$EXONSTART))
      endR=c(endR,max(cEx$EXONEND))
    }

    #sort again startR
    myorder=order(startR)
    startR=startR[myorder]
    endR=endR[myorder]

    # generate pseudo circularRNAs sequences
    #extract sequences and create a new transcript sequence
#    chrID=as.character(myex$EXONCHROM)[1]
#    x=substring(fasta_genome[chrID],startR,endR)
    if (is.null(fasta_genome)){
      suppressMessages(suppressWarnings(require(BSgenome.Hsapiens.UCSC.hg19)))
      fasta_genome<-BSgenome.Hsapiens.UCSC.hg19    
      chrID=paste0("chr",rep(as.character(myex$EXONCHROM)[1],length(startR)))
      gr = GRanges(seqnames=chrID, ranges=IRanges(start=startR, end=endR))
      x=as.character(getSeq(fasta_genome, gr))
    }else{
      chrID=as.character(myex$EXONCHROM)[1]
      x=substring(fasta_genome[chrID],startR,endR)
    }

    txSeq=paste(x,collapse="")
    # new create a pseudo-circularRNA sequence by adding l-1 nucleotides from the end of the last exon to the begin of the first exon
    l_1=readLen -1 -1
    if (l_1 > nchar(txSeq)) l_1=nchar(txSeq) - 1 #if the circRNA length < l_1
    circRNA.seq=paste(substring(txSeq, nchar(txSeq)-l_1),txSeq,sep="")
    circRNA.as_name=paste(mytx,"++",sep="") #do not forget to put the "circRNA", otherwises the format of the output file rawcount.txt of FuSeq is broken  
    #circRNA.fa=circRNA.seq
    #return(c(circRNA.fa,circRNA.as_name))  
    #add to result
    circRNA.tx.fa=c(circRNA.tx.fa,circRNA.seq)
    circRNA.tx.name = c(circRNA.tx.name, circRNA.as_name)
  }
  
  names(circRNA.tx.fa) = circRNA.tx.name
  x = NULL
  circ.name = NULL
  if(length(circRNA.tx.name) > 0){
    x=tapply(circRNA.tx.name, circRNA.tx.fa,c)

    circ.name = sapply(x, function(n){
      z=paste(n, collapse="=")
      paste(circRNA_name, "__", z, " AScircRNA", sep="")
    })
  }

  names(circ.name) = NULL
  return(list("fa" = names(x), "name" = circ.name))
}



