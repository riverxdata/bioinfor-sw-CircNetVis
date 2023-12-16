##############################
deleteCache=function(){
  system("find ./ -cmin +40 -type f -name 'circRNA_pseudoSeq*' -delete")
  system("find ./ -cmin +40 -type f -name 'RNAhybrid_out*' -delete")
  system("find ./ -cmin +40 -type f -name 'miRanda_out*' -delete")
  system("find ./ -cmin +40 -type f -name 'reactomejson*' -delete")
  return(NULL)
}
##############################
getmiRNAset <- function (){
  circInfo=circInfoData()
  cirmi.net=circInfo$cirmi.net
 
  return (cirmi.net$miRNA)
}
##############################
getGeneSet <- function (){
  circInfo=circInfoData()
  migen.net=circInfo$migen.net
  
  return(migen.net$gene)
}
##############################
parsingInput <- function(x, dataInfo=NULL, intype="any"){
  ##### Cutting inputs by newline, space, comma
  result=unlist(strsplit(x,"\n"))
  result=unique(trimws(result))
  
  result=unlist(strsplit(result," "))
  result=unique(trimws(result))
  
  result=unlist(strsplit(result,","))
  result=unique(trimws(result))
  
  result=result[result!=""]
  

  if (intype=="circRNA" & length(result)>0){
    p=grep("hsa",result)
    if (length(p)>0){
      p1=dataInfo$circRNA.ID %in% result[p]
      r1=dataInfo$normID[p1]
      r2=result[-p]
      result=c(r1,r2)
    }
  }

  return (result)
}
##############################
processcircRNAIDlist <-function(circRNA.in,hsaCirData_small,hsaMirData_small){
    RNAhybrid_cirmi.net = NULL
    miRanda_cirmi.net = NULL
    TargetScan_cirmi.net = NULL

    if (length(circRNA.in) == 0){
      out = list(
      RNAhybrid_cirmi.net = RNAhybrid_cirmi.net,
      miRanda_cirmi.net = miRanda_cirmi.net,
      TargetScan_cirmi.net = TargetScan_cirmi.net
      )    
      return(out)
    }

    x = lapply(circRNA.in, function(z)
    strsplit(z, "__")[[1]])
    x = do.call(rbind, x)
    chrVec = as.character(x[, 1])
    startVec = as.integer(x[, 2])
    endVec = as.integer(x[, 3])
    cat(paste(chrVec, "__", startVec, "__", endVec, sep = ""))
    normID_in = paste(chrVec, startVec, endVec, sep = "__")
    

    useRNAhybrid = TRUE
    normID_remain = normID_in
    RNAhybrid_cirmi.net = NULL
    
    ### 09Aug2022: use index instead of normID in the database
    #normID_in="X__139865340__139866824"
    p = which(hsaCirData_small$normID %in% normID_in)
    if (length(p) > 0) {
      # if existing in the database
      index_in = hsaCirData_small$index[p]
      
      ### collect results from RNAhybrid
      #connect to database
      con <- dbConnect(RSQLite::SQLite(), "RNAhybrid_cirmi_hg19.db")
      
      cirmi = lapply(index_in, function(x) {
        y = dbGetQuery(con,
                       paste(
                         "select * from RNAhybrid_cirmi_hg19 where normID='",
                         x,
                         "'",
                         sep = ""
                       ))
        return(y)
      })
      cirmi = do.call(rbind, cirmi)
      cirmi = cirmi[, c("normID", "miRNA", "mfe", "pval")]
      cirmi$normID = hsaCirData_small$normID[match(cirmi$normID, hsaCirData_small$index)]
      cirmi$miRNA = hsaMirData_small$miRNA[match(cirmi$miRNA, hsaMirData_small$index)]
      
      # Disconnect from the database
      dbDisconnect(con)
      
      if (nrow(cirmi) > 0) {
        #create RNAhybrid_cirmi.net
        RNAhybrid_cirmi.net = cirmi
        RNAhybrid_cirmi.net$pval = as.double(RNAhybrid_cirmi.net$pval)
        #RNAhybrid_cirmi.net$mfe=abs(as.double(RNAhybrid_cirmi.net$mfe))
        RNAhybrid_cirmi.net$mfe = as.double(RNAhybrid_cirmi.net$mfe)
        RNAhybrid_cirmi.net$miRNA = as.character(RNAhybrid_cirmi.net$miRNA)
        RNAhybrid_cirmi.net = as.data.frame(RNAhybrid_cirmi.net, stringsAsFactors =
                                              FALSE)
        
        pick = normID_in %in% cirmi$normID
        normID_remain = normID_in[which(!pick)]
        if (length(normID_remain) == 0)
          useRNAhybrid = FALSE
      }
      
    }
    
    
    if (useRNAhybrid) {
      #run RNAhybrid
      cirmi_remain = getCirmi_RNAhybrid(normID_remain)
      
      if (!is.null(cirmi_remain)) {
        cirmi_remain = cirmi_remain[, c("normID", "miRNA", "mfe", "pval")]
        RNAhybrid_cirmi.net = rbind(RNAhybrid_cirmi.net, cirmi_remain)
      }
    }
    
    RNAhybrid_cirmi.net$circRNA = RNAhybrid_cirmi.net$normID
    RNAhybrid_cirmi.net$pval = as.double(RNAhybrid_cirmi.net$pval)
    RNAhybrid_cirmi.net$mfe = as.double(RNAhybrid_cirmi.net$mfe) #use original value
    RNAhybrid_cirmi.net$weight = RNAhybrid_cirmi.net$mfe
    RNAhybrid_cirmi.net$miRNA = as.character(RNAhybrid_cirmi.net$miRNA)
#    cat("\n RNAhybrid_cirmi.net ", dim(RNAhybrid_cirmi.net))
    
    #####################
    usemiRanda = TRUE
    normID_remain = normID_in
    miRanda_cirmi.net = NULL
    
    ### 09Aug2022: use index instead of normID in the database
    #normID_in="X__139865340__139866824"
    p = which(hsaCirData_small$normID %in% normID_in)
    if (length(p) > 0) {
      # if existing in the database
      index_in = hsaCirData_small$index[p]
      
      ### collect results from miRanda
      #connect to database
      con2 <- dbConnect(RSQLite::SQLite(), "miRanda_cirmi_hg19.db")
      
      miRanda_cirmi = lapply(index_in, function(x) {
        y = dbGetQuery(con2,
                       paste(
                         "select * from miRanda_cirmi_hg19 where normID='",
                         x,
                         "'",
                         sep = ""
                       ))
        return(y)
      })
      miRanda_cirmi = do.call(rbind, miRanda_cirmi)
      #miRanda_cirmi = miRanda_cirmi[, c("normID", "miRNA", "Max_Score", "Max_Energy")]
      miRanda_cirmi = miRanda_cirmi[, c("miRNA","circRNA","normID","Tot_Score","Tot_Energy","Max_Score","Max_Energy","siteNum")]
      
      miRanda_cirmi$normID = hsaCirData_small$normID[match(miRanda_cirmi$normID, hsaCirData_small$index)]
      miRanda_cirmi$miRNA = hsaMirData_small$miRNA[match(miRanda_cirmi$miRNA, hsaMirData_small$index)]
      
      
      # Disconnect from the database
      dbDisconnect(con2)
      
      
      if (nrow(miRanda_cirmi) > 0) {
        miRanda_cirmi.net = miRanda_cirmi
        miRanda_cirmi.net$Max_Score = as.double(miRanda_cirmi.net$Max_Score)
        miRanda_cirmi.net$Max_Energy = as.double(miRanda_cirmi.net$Max_Energy)
        miRanda_cirmi.net$Tot_Score = as.double(miRanda_cirmi.net$Tot_Score)
        miRanda_cirmi.net$Tot_Energy = as.double(miRanda_cirmi.net$Tot_Energy)
        
        miRanda_cirmi.net = as.data.frame(miRanda_cirmi.net, stringsAsFactors =
                                            FALSE)
        pick = normID_in %in% miRanda_cirmi.net$normID
        normID_remain = normID_in[which(!pick)]
        if (length(normID_remain) == 0)
          usemiRanda = FALSE
      }
      
    }
    
    if (usemiRanda) {
      #run miRanda
      miRanda_cirmi_remain_res = getCirmi_miRanda(normID_remain)
      
      if (!is.null(miRanda_cirmi_remain_res)) {
#        miRanda_cirmi_remain = miRanda_cirmi_remain_res[, c("normID", "miRNA", "Max_Score", "Max_Energy")]
        miRanda_cirmi_remain = miRanda_cirmi_remain_res[, c("miRNA","circRNA","normID","Tot_Score","Tot_Energy","Max_Score","Max_Energy","siteNum")]
        miRanda_cirmi.net = rbind(miRanda_cirmi.net, miRanda_cirmi_remain)
      }
    }
    
    miRanda_cirmi.net$circRNA = miRanda_cirmi.net$normID
    miRanda_cirmi.net$Max_Score = as.double(miRanda_cirmi.net$Max_Score)
    miRanda_cirmi.net$Max_Energy = as.double(miRanda_cirmi.net$Max_Energy)
    miRanda_cirmi.net$weight = miRanda_cirmi.net$Max_Energy
    miRanda_cirmi.net$miRNA = as.character(miRanda_cirmi.net$miRNA)
#    cat("\n miRanda_cirmi.net ", dim(miRanda_cirmi.net))
    
    #####################
    useTargetScan = TRUE
    normID_remain = normID_in
    TargetScan_cirmi.net = NULL
    
    #normID_in="X__139865340__139866824"
    p = which(hsaCirData_small$normID %in% normID_in)
    if (length(p) > 0) {
      # if existing in the database
      index_in = hsaCirData_small$index[p]
      
      ### collect results from TargetScan
      #connect to database
      con3 <-
        dbConnect(RSQLite::SQLite(), "TargetScan_cirmi_hg19.db")
      
      TargetScan_cirmi = lapply(index_in, function(x) {
        y = dbGetQuery(
          con3,
          paste(
            "select * from TargetScan_cirmi_hg19 where normID='",
            x,
            "'",
            sep = ""
          )
        )
        return(y)
      })
      TargetScan_cirmi = do.call(rbind, TargetScan_cirmi)
      TargetScan_cirmi = TargetScan_cirmi[, c("normID",
                                              "miRNA",
                                              "st_6mer",
                                              "st_7mer1a",
                                              "st_7merm8",
                                              "st_8mer1a")]
      
      TargetScan_cirmi$normID = hsaCirData_small$normID[match(TargetScan_cirmi$normID, hsaCirData_small$index)]
      TargetScan_cirmi$miRNA = hsaMirData_small$miRNA[match(TargetScan_cirmi$miRNA, hsaMirData_small$index)]
      
      
      # Disconnect from the database
      dbDisconnect(con3)
      
      
      if (nrow(TargetScan_cirmi) > 0) {
        TargetScan_cirmi.net = TargetScan_cirmi
        TargetScan_cirmi.net$st_6mer = as.double(TargetScan_cirmi.net$st_6mer)
        TargetScan_cirmi.net$st_7mer1a = as.double(TargetScan_cirmi.net$st_7mer1a)
        TargetScan_cirmi.net$st_7merm8 = as.double(TargetScan_cirmi.net$st_7merm8)
        TargetScan_cirmi.net$st_8mer1a = as.double(TargetScan_cirmi.net$st_8mer1a)
        TargetScan_cirmi.net$miRNA = as.character(TargetScan_cirmi.net$miRNA)
        TargetScan_cirmi.net$normID = as.character(TargetScan_cirmi.net$normID)
        TargetScan_cirmi.net = as.data.frame(TargetScan_cirmi.net, stringsAsFactors =
                                               FALSE)
        pick = normID_in %in% TargetScan_cirmi.net$normID
        normID_remain = normID_in[which(!pick)]
        if (length(normID_remain) == 0)
          useTargetScan = FALSE
      }
      
    }
    
    if (useTargetScan) {
      #run TargetScan
      TargetScan_cirmi_remain_res = getCirmi_TargetScan(normID_remain)
      
      if (!is.null(TargetScan_cirmi_remain_res)) {
        TargetScan_cirmi_remain = TargetScan_cirmi_remain_res[, c("normID",
                                                                  "miRNA",
                                                                  "st_6mer",
                                                                  "st_7mer1a",
                                                                  "st_7merm8",
                                                                  "st_8mer1a")]
        TargetScan_cirmi.net = rbind(TargetScan_cirmi.net, TargetScan_cirmi_remain)
      }
    }
    
    TargetScan_cirmi.net$st_6mer = as.double(TargetScan_cirmi.net$st_6mer)
    TargetScan_cirmi.net$st_7mer1a = as.double(TargetScan_cirmi.net$st_7mer1a)
    TargetScan_cirmi.net$st_7merm8 = as.double(TargetScan_cirmi.net$st_7merm8)
    TargetScan_cirmi.net$st_8mer1a = as.double(TargetScan_cirmi.net$st_8mer1a)
    TargetScan_cirmi.net$miRNA = as.character(TargetScan_cirmi.net$miRNA)
    TargetScan_cirmi.net$normID = as.character(TargetScan_cirmi.net$normID)
    TargetScan_cirmi.net$circRNA = TargetScan_cirmi.net$normID
#    cat("\n TargetScan_cirmi.net ", dim(TargetScan_cirmi.net))
    
    #add results of TargetScan_cirmi.net, well more complicated now
    RNAhybrid_cirmi.net$int = paste(RNAhybrid_cirmi.net$normID, RNAhybrid_cirmi.net$miRNA, sep = ":")
    miRanda_cirmi.net$int = paste(miRanda_cirmi.net$normID, miRanda_cirmi.net$miRNA, sep = ":")
    TargetScan_cirmi.net$int = paste(TargetScan_cirmi.net$normID, TargetScan_cirmi.net$miRNA, sep = ":")
    out = list(
      #cirmi.net=cirmi.net,migen.net=migen.net, #do not prebuild cirmi.net and migen.net anymore
      RNAhybrid_cirmi.net = RNAhybrid_cirmi.net,
      miRanda_cirmi.net = miRanda_cirmi.net,
      TargetScan_cirmi.net = TargetScan_cirmi.net
    )
    return(out)
}

##############################

processcircRNAIDseqlist <-function(pseudoCircSeq){
  RNAhybrid_cirmi.net = NULL
  miRanda_cirmi.net = NULL
  TargetScan_cirmi.net = NULL

  ##### process RNAhybrid
  RNAhybrid_cirmi_remain_res = runRNAhybrid(pseudoCircSeq)
  if (!is.null(RNAhybrid_cirmi_remain_res)) {
      RNAhybrid_cirmi_remain_res$normID=RNAhybrid_cirmi_remain_res[,1]
      
      cirmiP=data.table(RNAhybrid_cirmi_remain_res)
      x=paste(cirmiP$miRNA,cirmiP$normID,sep="_")
      cirmiP$id=x
      x=NULL
      p=cirmiP[, .I[mfe == min(mfe)], by = id]
      cirmiP=cirmiP[p$V1,]
      p=cirmiP[, .I[pval == min(pval)], by = id]
      cirmiP=cirmiP[p$V1,]
      cirmiP=unique(cirmiP,by="id")
      cirmiP=cirmiP[,id:=NULL]
      RNAhybrid_cirmi_remain_res=as.data.frame(cirmiP)
      
      RNAhybrid_cirmi_remain_res = RNAhybrid_cirmi_remain_res[, c("normID", "miRNA", "mfe", "pval")]

      RNAhybrid_cirmi.net = RNAhybrid_cirmi_remain_res
      RNAhybrid_cirmi.net$circRNA = RNAhybrid_cirmi.net$normID
      RNAhybrid_cirmi.net$pval = as.double(RNAhybrid_cirmi.net$pval)
      RNAhybrid_cirmi.net$mfe = as.double(RNAhybrid_cirmi.net$mfe) #use original value
      RNAhybrid_cirmi.net$weight = RNAhybrid_cirmi.net$mfe
      RNAhybrid_cirmi.net$miRNA = as.character(RNAhybrid_cirmi.net$miRNA)
      RNAhybrid_cirmi.net$int = paste(RNAhybrid_cirmi.net$normID, RNAhybrid_cirmi.net$miRNA, sep = ":")
      #cat("\n RNAhybrid_cirmi.net ", dim(RNAhybrid_cirmi.net))
  }


  ##### process Miranda
  miRanda_cirmi_remain_res = runMiranda(pseudoCircSeq)
  if (!is.null(miRanda_cirmi_remain_res)){
    miRanda_cirmi_remain_res$normID=as.character(miRanda_cirmi_remain_res$circRNA)
    #keep unique  miRNA-normID, update max_score and max_energy
    cirmiSum=miRanda_cirmi_remain_res
    x=paste(cirmiSum$miRNA,cirmiSum$normID,sep="_")
    cirmiSum$id=x
    x=NULL
    cirmiSum=data.table(cirmiSum)
    p=cirmiSum[, .I[Max_Score == max(Max_Score)], by = id]
    cirmiSum=cirmiSum[p$V1,]
    p=cirmiSum[, .I[Tot_Score == max(Tot_Score)], by = id]
    cirmiSum=cirmiSum[p$V1,]
    cirmiSum=unique(cirmiSum,by="id")
    cirmiSum=cirmiSum[,id:=NULL]
    
    cirmiSum=as.data.frame(cirmiSum)
    miRanda_cirmi.net=cirmiSum[,c("miRNA","circRNA","normID","Tot_Score","Tot_Energy","Max_Score","Max_Energy","siteNum")]
    
    # x=paste(miRanda_cirmi_remain_res$miRNA,miRanda_cirmi_remain_res$normID,sep="_")
    # #y=table(x)
    # x1=tapply(miRanda_cirmi_remain_res$Max_Score,x,max)
    # x2=tapply(miRanda_cirmi_remain_res$Max_Energy,x,max)
    # pick=match(names(x1),x)
    # miRanda_cirmi.net=miRanda_cirmi_remain_res[pick,]
    # miRanda_cirmi.net$Max_Score=x1
    # miRanda_cirmi.net$Max_Energy=x2
    # miRanda_cirmi.net=miRanda_cirmi.net[,c("miRNA","circRNA","Max_Score","Max_Energy","normID")]
    
    miRanda_cirmi.net$circRNA = miRanda_cirmi.net$normID
    miRanda_cirmi.net$Max_Score = as.double(miRanda_cirmi.net$Max_Score)
    miRanda_cirmi.net$Max_Energy = as.double(miRanda_cirmi.net$Max_Energy)
    miRanda_cirmi.net$weight = miRanda_cirmi.net$Max_Energy
    miRanda_cirmi.net$miRNA = as.character(miRanda_cirmi.net$miRNA)
    miRanda_cirmi.net$int = paste(miRanda_cirmi.net$normID, miRanda_cirmi.net$miRNA, sep = ":")
  }


  ##### process TargetScan  
  cirmiP = runTargetScan(pseudoCircSeq)
  if (!is.null(cirmiP)){
    #get normID
    normID=cirmiP$a_Gene_ID
    cirmiP$normID=normID
    #make indicies for AS CircRNA
    #asInd=c(1:length(cirmiP$a_Gene_ID))
    asInd=normID
    cirmiP$asInd=asInd
    
    cirmiP=cirmiP[,c("asInd","normID","miRNA_family_ID","Site_type"), with=FALSE]
    dim(cirmiP)
    
    cirmiExport=cirmiP
    combsExport=paste0(cirmiExport$asInd,cirmiExport$miRNA_family_ID,sep="_")
    p=!duplicated(combsExport)
    cirmiExport=cirmiExport[p,]
    combsExport=combsExport[p]
    
    mytabulateFunc <-function(combs){
      # this code is much faster than n=table(combs)
      x=unique(combs)
      y=seq(length(x))
      names(y)=x
      combs1=y[combs]
      n=tabulate(combs1)
      names(n)=x
      return(n)
    }
    
    #siteType=c("6mer","7mer-1a","7mer-m8","8mer-1a")
    # summarize the results, get the total number of binding sites for each type
    # NOTE: we consider the interaction is for an AS circRNA
    d=cirmiP[cirmiP$Site_type=="6mer",] 
    combs=paste0(d$asInd,d$miRNA_family_ID,sep="_")
    n=mytabulateFunc(combs)
    cirmiExport$st_6mer=0
    p=which(combsExport%in%combs)
    cirmiExport$st_6mer[p]=n[combsExport[p]]
    
    d=cirmiP[cirmiP$Site_type=="7mer-1a",]  
    combs=paste0(d$asInd,d$miRNA_family_ID,sep="_")
    n=mytabulateFunc(combs)
    cirmiExport$st_7mer1a=0
    p=which(combsExport%in%combs)
    cirmiExport$st_7mer1a[p]=n[combsExport[p]]
    
    d=cirmiP[cirmiP$Site_type=="7mer-m8",]  
    combs=paste0(d$asInd,d$miRNA_family_ID,sep="_")
    n=mytabulateFunc(combs)
    cirmiExport$st_7merm8=0
    p=which(combsExport%in%combs)
    cirmiExport$st_7merm8[p]=n[combsExport[p]]
    
    d=cirmiP[cirmiP$Site_type=="8mer-1a",]  
    combs=paste0(d$asInd,d$miRNA_family_ID,sep="_")
    n=mytabulateFunc(combs)
    cirmiExport$st_8mer1a=0
    p=which(combsExport%in%combs)
    cirmiExport$st_8mer1a[p]=n[combsExport[p]]
    
    cirmiExport=cirmiExport[,c("asInd","normID","miRNA_family_ID","st_6mer","st_7mer1a","st_7merm8","st_8mer1a"), with=FALSE]
    colnames(cirmiExport)[1]="AScircRNA" #asInd
    colnames(cirmiExport)[3]="miRNA" #miRNA_family_ID
    cirmiExport=as.data.frame(cirmiExport)

    TargetScan_cirmi.net=cirmiExport
    cirmiExport=NULL
    TargetScan_cirmi.net$st_6mer = as.double(TargetScan_cirmi.net$st_6mer)
    TargetScan_cirmi.net$st_7mer1a = as.double(TargetScan_cirmi.net$st_7mer1a)
    TargetScan_cirmi.net$st_7merm8 = as.double(TargetScan_cirmi.net$st_7merm8)
    TargetScan_cirmi.net$st_8mer1a = as.double(TargetScan_cirmi.net$st_8mer1a)
    TargetScan_cirmi.net$miRNA = as.character(TargetScan_cirmi.net$miRNA)
    TargetScan_cirmi.net$normID = as.character(TargetScan_cirmi.net$normID)
    TargetScan_cirmi.net$circRNA = TargetScan_cirmi.net$normID
    TargetScan_cirmi.net$int = paste(TargetScan_cirmi.net$normID, TargetScan_cirmi.net$miRNA, sep = ":")
  }



  out = list(
      #cirmi.net=cirmi.net,migen.net=migen.net, #do not prebuild cirmi.net and migen.net anymore
      RNAhybrid_cirmi.net = RNAhybrid_cirmi.net,
      miRanda_cirmi.net = miRanda_cirmi.net,
      TargetScan_cirmi.net = TargetScan_cirmi.net
  )
  return(out)
}
