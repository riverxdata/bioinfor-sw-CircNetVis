#23Feb2023/Nghia: add targetScan in the filter of cir-mi interaction
library(shiny)
library(shinysky)
library(shinyjs)
library(visNetwork)
#library(igraph)
library(doParallel)
library(DT)
library(rjson)
library(RSQLite)
library(dbplyr)
library(DBI)
library(dplyr)
library(data.table)

#setwd("C:/Nghiavtr/Workplaces/Projects/circularRNAs/mirRNA_interaction/CircRNAanalysis/ShinyRuning/CircNetVis")
###Initializating global variables

#shiny::runApp(system.file("shiny", package = "visNetwork"))
suppressMessages(suppressWarnings(require(Biostrings)))
#suppressMessages(suppressWarnings(require(BSgenome.Hsapiens.UCSC.hg19)))
#load targetscan
load("Annotation_data/targetscan_v72.RData")
load("Annotation_data/hsa_circRNA_mirRNA_annotation_processed_small.RData")
load("Annotation_data/circRNA_RBP.RData") ###  RBP_dat_NET: "circRNA" "RBP"     "SiteNum" "Line"
#load global functions
source("server/functions_global.R")
source("server/functions.R")

getcircinfo = function(input) {
  #  res = reactive({
  res = eventReactive(input$goButton, {
    circRNA.in = parsingInput(input$circRNAinput, hsaCirData_small, intype =
                                "circRNA")
    
    if (length(circRNA.in) > 0){
      circRNAIDlist_res<-processcircRNAIDlist(circRNA.in,hsaCirData_small,hsaMirData_small)
      out = list(
        RNAhybrid_cirmi.net = circRNAIDlist_res$RNAhybrid_cirmi.net,
        miRanda_cirmi.net = circRNAIDlist_res$miRanda_cirmi.net,
        TargetScan_cirmi.net = circRNAIDlist_res$TargetScan_cirmi.net
      )
    }
    
    
    pseudoCircSeq = NULL
    if (input$cbShowcircRNAseqlist && input$circRNAseqlist != "") {
      result = input$circRNAseqlist
      result = unlist(strsplit(result, "\n"))
      result = unique(trimws(result))
      # id = seq(1, length(result), 2)
      # pseudoSeq.name = result[id]
      # pseudoSeq.name = gsub(">","",pseudoSeq.name)
      # pseudoCircSeq = result[-id]
      
      tmp = tempfile()
      write(result,tmp)
      y=readDNAStringSet(tmp)
      system(paste0("rm ",tmp))
      pseudoSeq.name=names(y)
      pseudoSeq.fa=as.character(y)
      
      pseudoCircSeq = DNAStringSet()
      for (i in 1:length(pseudoSeq.fa)){
        rawSeq=pseudoSeq.fa[i]
        #create pseudosequence
        readLen=28 # 28 is the max length of microRNAs
        l_1=readLen -1 -1
        if (l_1 > nchar(rawSeq)) l_1=nchar(rawSeq) - 1 #if the circRNA length < l_1
        rawSeq=paste(substring(rawSeq, nchar(rawSeq)-l_1),rawSeq,sep="")
        rawSeq = DNAStringSet(rawSeq)
        names(rawSeq) = pseudoSeq.name[i]
        pseudoCircSeq=c(pseudoCircSeq,rawSeq)
      }
    }


    if (!is.null(pseudoCircSeq)){
      circRNAIDseq_res<-processcircRNAIDseqlist(pseudoCircSeq)
      out = list(
        RNAhybrid_cirmi.net = circRNAIDseq_res$RNAhybrid_cirmi.net,
        miRanda_cirmi.net = circRNAIDseq_res$miRanda_cirmi.net,
        TargetScan_cirmi.net = circRNAIDseq_res$TargetScan_cirmi.net
      )
      
      #combine results
      # if (length(circRNA.in) > 0){
      #   out = list(
      #     RNAhybrid_cirmi.net = rbind(circRNAIDlist_res$RNAhybrid_cirmi.net,circRNAIDseq_res$RNAhybrid_cirmi.net),
      #     miRanda_cirmi.net = rbind(circRNAIDlist_res$miRanda_cirmi.net,circRNAIDseq_res$miRanda_cirmi.net),
      #     TargetScan_cirmi.net = rbind(circRNAIDlist_res$TargetScan_cirmi.net,circRNAIDseq_res$TargetScan_cirmi.net)
      #       )
      #   }
    }

    return(out)
    
  }) ### The end
  return(res)
}


getCircRBPinfo <- function(input) {
  circRNA.in = parsingInput(input$circRNAinput, hsaCirData_small, intype = "cirRNA")
  #Nghia/27June2023: add the circRNA names from the input with fasta sequence
   if (input$cbShowcircRNAseqlist && input$circRNAseqlist != "") {
      result = input$circRNAseqlist
      result = unlist(strsplit(result, "\n"))
      result = unique(trimws(result))
      tmp = tempfile()
      write(result,tmp)
      y=readDNAStringSet(tmp)
      system(paste0("rm ",tmp))
      pseudoSeq.name=names(y)
      if (length(circRNA.in) > 0)  circRNA.in= unique(c(circRNA.in,pseudoSeq.name)) else circRNA.in= pseudoSeq.name
    }
  #Nghia/27June2023: done 

  CircRBP_int.net = NULL
  
  if (length(circRNA.in) > 0) {
    pick1 = which (
      hsaCirData_small$normID %in% circRNA.in |
        hsaCirData_small$circRNA.ID %in% circRNA.in
    )
    hsaCirData_small = hsaCirData_small[pick1,]
    pick2 = which (RBP_dat_NET$circRNA %in% hsaCirData_small$circRNA.ID)
    CircRBP_int.net1 = RBP_dat_NET[pick2,]
    
    ####obtaining max no of sites
    siteNumLim = as.integer(input$siteNumLim)
    pick3 = which (CircRBP_int.net1$SiteNum >=  siteNumLim)
    
    if (length (pick3) > 0) {
      CircRBP_int.net = CircRBP_int.net1[pick3,]
      CircRBP_int.net$normID = hsaCirData_small$normID[match(CircRBP_int.net$circRNA, hsaCirData_small$circRNA.ID)]
    }
  }
  return(CircRBP_int.net)
}

##############################################
# Define server logic required to draw a histogram
shinyServer(
  function(input, output, session) {
    # mykey = reactive({
    #
    #   if(input$geneSearch == "")
    #     "6__111281617__111283708"
    #   else
    #     input$geneSearch
    # })
    
    #######################
    deleteCacheRes = deleteCache()
    
    ###### common variables in shinyServer()
    circInfoData <- getcircinfo(input)
    
    ############ Displaying cirRNA-miRNA-genes network
    mynetwork <- reactive({
      #cat("\n Plot network")
      
      #obtaining threshold values for filtering miRNA and genes
      pval.thres = as.double(input$RNAhybrid_pvalThres)
      mfe.thres = as.double(input$RNAhybrid_mfeThres)
      weightPercentile.migene.thres = as.double(input$weightPercentilemigeneThres)
      weightScore.migene.thres = as.double(input$weightScoremigeneThres)
      
      
      miRandaScore.thres = as.double(input$miRandaScoreThres)
      
      st_6mer.thres = as.double(input$TargetScan_6mer)
      st_7mer1a.thres = as.double(input$TargetScan_7mer1a)
      st_7merm8.thres = as.double(input$TargetScan_7merm8)
      st_8mer1a.thres = as.double(input$TargetScan_8mer1a)
      
      if (is.na(st_6mer.thres)) st_6mer.thres=0
      if (is.na(st_7mer1a.thres)) st_7mer1a.thres=0
      if (is.na(st_7merm8.thres)) st_7merm8.thres=0
      if (is.na(st_8mer1a.thres)) st_8mer1a.thres=0
      ####obtaining max no of miRNA/genes

      genenumlim = input$genenumlim
      ###########################
      withProgress(message = 'Preparing the circRNA-miRNA-mRNA interaction network, please wait until completed!', value = 0, {
      #    Sys.sleep(0.25)
      circInfo = circInfoData()

      })
      
      if (!is.null(circInfo)) {
        RNAhybrid_cirmi.net = circInfo$RNAhybrid_cirmi.net
        miRanda_cirmi.net = circInfo$miRanda_cirmi.net
        TargetScan_cirmi.net = circInfo$TargetScan_cirmi.net

        #filter cirmi interactions by miRanda score
        if (input$cbUseMiRanda) {
          pick = which(miRanda_cirmi.net$Max_Score > miRandaScore.thres)
          miRanda_cirmi.net = miRanda_cirmi.net[pick, ]
        }
        
        #filter cirmi interactions by p-value - UseRNAhybrid
        if (input$cbUseRNAhybrid) {
          pick = which(
            RNAhybrid_cirmi.net$pval <= pval.thres &
              RNAhybrid_cirmi.net$mfe <= mfe.thres
          )
          #cat("\n condition UseRNAhybrid ",length(pick))
          RNAhybrid_cirmi.net = RNAhybrid_cirmi.net[pick, ]
        }
        
        #filter cirmi interactions by siteType - TargetScan_cirmi.net
        if (input$cbUseTargetScan) {
          pick1 = which(TargetScan_cirmi.net$st_6mer >= st_6mer.thres)
          pick2 = which(TargetScan_cirmi.net$st_7mer1a >= st_7mer1a.thres)
          pick3 = which(TargetScan_cirmi.net$st_7merm8 >= st_7merm8.thres)
          pick4 = which(TargetScan_cirmi.net$st_8mer1a >= st_8mer1a.thres)
          pick = unique(c(pick1, pick2, pick3, pick4))
          TargetScan_cirmi.net = TargetScan_cirmi.net[pick, ]
        }
        
        ### now combine
        TargetScan_cirmi.net$method = "TargetScan"
        miRanda_cirmi.net$method = "miRanda"
        RNAhybrid_cirmi.net$method = "RNAhybrid"
        cirmi.net = TargetScan_cirmi.net[, c("circRNA", "miRNA", "int", "method")]
        cirmi.net = rbind(cirmi.net, miRanda_cirmi.net[, c("circRNA", "miRNA", "int", "method")])
        cirmi.net = rbind(cirmi.net, RNAhybrid_cirmi.net[, c("circRNA", "miRNA", "int", "method")])
        
        p = input$cbUseMiRanda |
          input$cbUseRNAhybrid | input$cbUseTargetScan
        if (!p) {
          if (!input$cbUseMiRanda)
            cirmi.net = cirmi.net[-which(cirmi.net$method == "miRanda"), ]
          if (input$cbUseRNAhybrid)
            cirmi.net = cirmi.net[-which(cirmi.net$method == "RNAhybrid"), ]
          if (input$cbUseTargetScan)
            cirmi.net = cirmi.net[-which(cirmi.net$method == "TargetScan"), ]
        }
        
        #keep only the one found in all active methods
        p = tapply(cirmi.net$method, cirmi.net$int, length)
        p = p[p == length(unique(cirmi.net$method))]
        cirmi.net = cirmi.net[cirmi.net$int %in% names(p), ]
        cirmi.net = cirmi.net[!duplicated(cirmi.net$int), ]
        
        
        ### then create migen.net
        ## get miRNA-mRNA interactions
        miRNAset = unique(as.character(cirmi.net$miRNA))
        pick = targetscan_db$miRNA %in% miRNAset
        migen.net = targetscan_db[pick, ]
        migen.net = migen.net[order(migen.net$weighted.context...score, decreasing = FALSE), ]
        migen.net$weight = migen.net$weighted.context...score
        migen.net$gene = migen.net$Gene.Symbol
        
        #######obtaining interested miRNA list
        #      cat (paste("\n Testing miRNA list 1:",input$cbShowmiRNAList,",", input$miRNAList,sep=""))
        if (input$cbShowmiRNAList && input$miRNAList != "") {
          miRNAlist = parsingInput(input$miRNAList, intype = "miRNA")
          #        cat (paste("\n Testing miRNA list 2:",miRNAlist,sep=""))
          if (length(miRNAlist) > 0) {
            pick = which(cirmi.net$miRNA %in% miRNAlist)
            
            miRNAlist2=paste0(miRNAlist,"-")
            pick1 = sapply(miRNAlist2, function(x) grep(x,cirmi.net$miRNA))
            pick1 = unlist(pick1)
            if (length(pick1) > 0) pick=c(pick,pick1)
            
            cirmi.net = cirmi.net[pick, ]
            
          }
        }
        
        #######obtaining interested gene list
        #cat (paste("\n Testing gene list 3:",input$cbShowGeneList,",", input$miRNAList,sep=""))
        if (input$cbShowGeneList && input$geneList != "") {
          geneList = parsingInput(input$geneList, intype = "mRNA")
          if (length (geneList) > 0) {
            pick = which(migen.net$gene %in% geneList)
            migen.net = migen.net[pick, ]
          }
        }
        
        
        # migen.net filtered by weight threshold - targetscan
        pick = which(
          migen.net$weighted.context...score.percentile >= weightPercentile.migene.thres &
            migen.net$weighted.context...score <= weightScore.migene.thres
        )
        migen.net = migen.net[pick, ]
        
        # cirmi.net=cirmi.net[,c("circRNA","miRNA","weight")]
        # migen.net=migen.net[,c("miRNA","gene","weight")]
        
        cirmi.net = cirmi.net[, c("circRNA", "miRNA")]
        migen.net = migen.net[, c("miRNA", "gene")]
        
        circRNAset = unique(as.character(cirmi.net$circRNA))
        
        # #filtered by max no. of miRNA/genes
        # if (RNAhybrid_numlim=="All") RNAhybrid_numlim=Inf else RNAhybrid_numlim=as.integer(RNAhybrid_numlim)
        # if (RNAhybrid_numlim!=Inf){
        #   #sort by p-value of RNAhybrid
        #   myorder=order(cirmi.net$pval, decreasing = FALSE)
        #   selectID=NULL
        #   for (mycir in circRNAset){
        #     pick=which(cirmi.net$circRNA %in% mycir)
        #     if (length(pick)>RNAhybrid_numlim) pick=pick[1:RNAhybrid_numlim]
        #     selectID=c(selectID,pick)
        #   }
        #   cirmi.net=cirmi.net[selectID,]
        # }
        
        
        
        miRNAset = unique(as.character(cirmi.net$miRNA))
        ######Testing
        #cat(paste("\n nb of miRNA: ",length(miRNAset),sep=""))
        
        if (genenumlim == "All")
          genenumlim = Inf
        else
          genenumlim = as.integer(genenumlim)
        selectID = NULL
        for (mymi in miRNAset) {
          pick = which(migen.net$miRNA %in% mymi)
          if (length(pick) > genenumlim)
            pick = pick[1:genenumlim]
          selectID = c(selectID, pick)
        }
        migen.net = migen.net[selectID, ]
        
        #process nodes
        nodecir = as.character(unique(cirmi.net[, 1]))
        nodemi = as.character(unique(cirmi.net[, 2]))
        nodegen = as.character(unique(migen.net[, 2]))
        
        nodecir.id = c(1:length(nodecir))
        nodemi.id = length(nodecir.id) + c(1:length(nodemi))
        nodegen.id = length(nodecir.id) + length(nodemi.id) + c(1:length(nodegen))
        
        if (length(nodecir) == 0)
          return(NULL)
        
        #    nodecir.shape=rep("square",length(nodecir.id))
        #    nodemi.shape=rep("triangle",length(nodemi.id))
        #    nodegen.shape=rep("dot",length(nodegen.id))
        
        nodecir.shape = rep("square", length(nodecir.id))
        nodemi.shape = rep("crectangle", length(nodemi.id))
        nodegen.shape = rep("oval", length(nodegen.id))#circle
        
        
        nodecir.color = rep("darkred", length(nodecir.id))
        nodemi.color = rep("lightblue", length(nodemi.id))
        nodegen.color = rep("grey", length(nodegen.id))
        
        nodecir.group = rep("circRNA", length(nodecir.id))
        nodemi.group = rep("miRNA", length(nodemi.id))
        nodegen.group = rep("gene", length(nodegen.id))
        
        #process edges
        cirmi.from = nodecir.id[match(cirmi.net[, 1], nodecir)]
        cirmi.to = nodemi.id[match(cirmi.net[, 2], nodemi)]
        
        migen.from = nodemi.id[match(migen.net[, 1], nodemi)]
        migen.to = nodegen.id[match(migen.net[, 2], nodegen)]
        
        # only circRNA-miRNA
        
        # replace normID by circBase ID
        nodecir_show = nodecir
        if (input$cbShowcircbaseID_cirmi) {
          m = match(nodecir_show, hsaCirData_small$normID)
          p1 = which(!is.na(m))
          if (length(p1) > 0)
            nodecir_show[p1] = hsaCirData_small$circRNA.ID[m[p1]]
        }
        
        #nodes.label=c(nodecir,nodemi)
        nodes.label = c(nodecir_show, nodemi)
        
        nodes.id = c(nodecir.id, nodemi.id)
        nodes.shape = c(nodecir.shape, nodemi.shape)
        nodes.color = c(nodecir.color, nodemi.color)
        nodes.group = c(nodecir.group, nodemi.group)
        edges.from = c(cirmi.from)
        edges.to = c(cirmi.to)
        
        # full circRNA-miRNA-genes
        if (length(nodegen) > 0) {
          nodes.id = c(nodes.id, nodegen.id)
          nodes.label = c(nodes.label, nodegen)
          nodes.shape = c(nodes.shape, nodegen.shape)
          nodes.color = c(nodes.color, nodegen.color)
          nodes.group = c(nodes.group, nodegen.group)
          edges.from = c(edges.from, migen.from)
          edges.to = c(edges.to, migen.to)
        }
        
        # customization adding more variables (see visNodes and visEdges)
        nodes <- data.frame(
          id = nodes.id,
          label = nodes.label,
          # labels
          shape = nodes.shape,
          color = nodes.color,
          group = nodes.group
        )
        
        edges.label = NULL
        edges <- data.frame(from = edges.from, to = edges.to)
        
        
        net = visNetwork(nodes, edges, height = "100%", width = "100%") %>% visOptions(highlightNearest = TRUE)  %>% visEdges(smooth = FALSE) %>% visExport(type = "png") %>% visPhysics(stabilization = TRUE)
        return(net)
      }
      
    })
    
    output$circRNAnetwork <- renderVisNetwork({
      mynetwork()
    })
    ################ VisNetwork
  
  
    ############ Displaying cirRNA-RBP network
    output$circRBPnetwork <- renderVisNetwork({
#      suppressMessages(suppressWarnings(require(Biostrings)))
#      load("Annotation_data/circRNA_RBP.RData") ###  RBP_dat_NET: "circRNA" "RBP"     "SiteNum" "Line"
#      load("Annotation_data/hsa_hg19_circRNA.RData") ####### hsaData : circRNA.ID, normID
      
      #cat("\n Plot cirRNA-RBP network")
      ######### obtaining cirRNA normID
      CircRBP_int.net=getCircRBPinfo(input)
      ##########################
      if (!is.null(CircRBP_int.net))  {
        # ####obtaining max no of sites
        # siteNumLim = input$siteNumLim
        # pick3 = which (CircRBP_int.net$SiteNum >=  siteNumLim)
        
        ###############################
#        if (length (pick3) > 0) {
          # CircRBP_int.net = CircRBP_int.net[pick3, ]
          # CircRBP_int.net$normID = hsaCirData_small$normID[match(CircRBP_int.net$circRNA, hsaCirData_small$circRNA.ID)] #### Adding normID
          # #cat (paste("\n RBP_dat_NET$SiteNum=",RBP_dat_NET$SiteNum, sep = ""))
          
          ####Showing only top site RBP
          if (input$cbShowTopSiteRBP) {
            CircRBP_int.net <-
              merge(aggregate(SiteNum ~ normID, data = CircRBP_int.net, max),
                    CircRBP_int.net,
                    all.x = T)
            #cat (paste("\n RBP_dat_NET$SiteNum=",RBP_dat_NET$SiteNum, sep = ""))
          }
          #########Filtering by RBP of interests
          if (input$cbShowRBPList) {
            RBPList = parsingInput(input$RBPList, intype = "RBP")
            
            if (length(RBPList) > 0) {
              pick4 = which (CircRBP_int.net$RBP %in% RBPList)
              CircRBP_int.net = CircRBP_int.net[pick4, ]
            }
          }
          
          #process nodes
          nodecir = as.character(unique(CircRBP_int.net[, "normID"]))
          nodeRBP = as.character(unique(CircRBP_int.net[, "RBP"]))
          
          nodecir.id = c(1:length(nodecir))
          nodeRBP.id = length(nodecir.id) + c(1:length(nodeRBP))
          
          nodecir.shape = rep("square", length(nodecir.id))
          nodeRBP.shape = rep("oval", length(nodeRBP.id))
          
          nodecir.color = rep("darkred", length(nodecir.id))
          nodeRBP.color = rep("orange", length(nodeRBP.id))
          
          nodecir.group = rep("circRNA", length(nodecir.id))
          nodeRBP.group = rep("RBP", length(nodeRBP.id))
          
          #process edges
          cirRBP.from = nodecir.id[match(CircRBP_int.net[, "normID"], nodecir)]
          cirRBP.to = nodeRBP.id[match(CircRBP_int.net[, "RBP"], nodeRBP)]
          
          if (input$cbShowSiteNum) {
            nodeRBPLabel = paste0(nodeRBP, " : ", unique(CircRBP_int.net$SiteNum[match(CircRBP_int.net[, "RBP"], nodeRBP)]), sep =
                                    "")
          }
          else
            nodeRBPLabel = nodeRBP
          # only circRNA-RBP
          # replace normID by circBase ID
          nodecir_show = nodecir
          if (input$cbShowcircbaseID_cirRBP) {
            m = match(nodecir_show, hsaCirData_small$normID)
            p1 = which(!is.na(m))
            if (length(p1) > 0)
              nodecir_show[p1] = hsaCirData_small$circRNA.ID[m[p1]]
          }
          #nodes.label=c(nodecir,nodeRBPLabel)
          nodes.label = c(nodecir_show, nodeRBPLabel)
          
          nodes.id = c(nodecir.id, nodeRBP.id)
          nodes.shape = c(nodecir.shape, nodeRBP.shape)
          nodes.color = c(nodecir.color, nodeRBP.color)
          nodes.group = c(nodecir.group, nodeRBP.group)
          edges.from = c(cirRBP.from)
          edges.to = c(cirRBP.to)
          
          #cat (paste("\n nodes.id=",nodes.id, sep = ""))
          #cat (paste("\n edges.from=",edges.from, sep = ""))
          # customization adding more variables (see visNodes and visEdges)
          nodes <- data.frame(
            id = nodes.id,
            label = nodes.label,
            # labels
            shape = nodes.shape,
            color = nodes.color,
            group = nodes.group
          )
          
          #cat (paste("\n testing point 3, edges.to=",edges.to, sep = ""))
          edges.label = NULL
          edges <- data.frame(from = edges.from, to = edges.to)
          
          #obtaining network height
          #networkHeight = as.integer(input$circRNAnetworkHeight)
          #networkMaxHeight = 1000
          #cat(paste('\n height = "',networkHeight*100/networkMaxHeight,'%"',sep=""))
          
          visNetwork(nodes, edges, height = "800px", width = "100%") %>% visOptions(highlightNearest = TRUE) %>%
            visPhysics(enabled = TRUE) %>% visEdges(smooth = FALSE) %>% visExport()
#        }
      }
      
    })
    
  
  # observe({
  #   input$getSelectedNodes
  #   visNetworkProxy("circRNAnetwork") %>%visGetNodes()
  #   #visNetworkProxy("circRNAnetwork") %>%visGetSelectedNodes()
  # })
  
  
    observeEvent(input$cbdlhtml, {
      if (input$cbdlhtml) {
        shinyjs::show("downloadNetwork")
        #      shinyjs::show("downloadCirmi")
        #      shinyjs::show("downloadmim")
        visNetworkProxy("circRNAnetwork") %>% visGetPositions() %>% visGetNodes() %>% visGetEdges()
      }
      else{
        shinyjs::hide("downloadNetwork")
        #      shinyjs::hide("downloadCirmi")
        #      shinyjs::hide("downloadmim")
      }
    })
  
  
    output$downloadNetwork <- downloadHandler(
      filename = function() {
        paste('network-', Sys.Date(), '.html', sep = '')
      },
      content = function(con) {
        updateCheckboxInput(session, "cbdlhtml", value = FALSE)
        #input$cbdlhtml=FALSE
        #shinyjs::hide("downloadNetwork")
        #########
        cat("\n export network to file ")
        net = mynetwork()
        nodes <- net$x$nodes
        edges <- net$x$edges
        
        network_positions <- input$circRNAnetwork_positions
        ################
        nodes_positions <- NULL
        #      cat(paste("\n Testing exporting network 1, nodes length=",length(nodes),sep=""))
        # format positions
        if (!is.null(network_positions)) {
          nodes_positions <-
            do.call("rbind", lapply(network_positions, function(x) {
              data.frame(x = x$x, y = x$y)
            }))
          nodes_positions$id <- names(network_positions)
          #        cat(paste("\n Testing exporting network 2 - id, length=",length (nodes_positions$id),sep=""))
        }
        
        cat(paste(
          "\n Testing exporting network 3, length=",
          length(nodes_positions),
          sep = ""
        ))
        ###### Saving the network
        if (!is.null(nodes_positions)) {
          nodes_save <- merge(nodes, nodes_positions, by = "id", all = T)
        } else  {
          nodes_save <- nodes
        }#else
        
        netExport = visNetwork(nodes = nodes_save,
                               edges = edges,
                               height = "800px") %>% visOptions(highlightNearest = TRUE) %>% visPhysics(enabled = FALSE) %>% visEdges(smooth = FALSE)
        
        visSave(netExport, con)
        
        
      }
    )
  
    getIntTables <- function() {
      circInfo = circInfoData()
      
      net = mynetwork()
      nodes <- net$x$nodes
      edges <- net$x$edges
      edges_from = nodes$label[match(edges$from, nodes$id)] #circRNA
      edges_to = nodes$label[match(edges$to, nodes$id)] #miRNA
      
      #if using circRNA.ID, then convert to normID
      m = match(edges_from, hsaCirData_small$circRNA.ID)
      p1 = which(!is.na(m))
      if (length(p1) > 0)
        edges_from[p1] = hsaCirData_small$normID[m[p1]]
      
      
      ### collect the current results
      RNAhybrid_cirmi.net = circInfo$RNAhybrid_cirmi.net
      miRanda_cirmi.net = circInfo$miRanda_cirmi.net
      TargetScan_cirmi.net = circInfo$TargetScan_cirmi.net
      
      RNAhybrid_cirmi.net = RNAhybrid_cirmi.net[which(
        RNAhybrid_cirmi.net$miRNA %in% c(edges_from, edges_to) &
          RNAhybrid_cirmi.net$normID %in% c(edges_from, edges_to)
      ), ]
      miRanda_cirmi.net = miRanda_cirmi.net[which(
        miRanda_cirmi.net$miRNA %in% c(edges_from, edges_to) &
          miRanda_cirmi.net$normID %in% c(edges_from, edges_to)
      ), ]
      TargetScan_cirmi.net = TargetScan_cirmi.net[which(
        TargetScan_cirmi.net$miRNA %in% c(edges_from, edges_to) &
          TargetScan_cirmi.net$normID %in% c(edges_from, edges_to)
      ), ]
      
      ### now combine
      cirmi.net.short = RNAhybrid_cirmi.net[, c("circRNA", "miRNA", "int")]
      cirmi.net.short = cirmi.net.short[cirmi.net.short$int %in% miRanda_cirmi.net$int, ]
      cirmi.net.short = cirmi.net.short[cirmi.net.short$int %in% TargetScan_cirmi.net$int, ]
      
      ### then create migen.net
      ## get miRNA-mRNA interactions
      miRNAset = unique(as.character(cirmi.net.short$miRNA))
      pick = targetscan_db$miRNA %in% miRNAset
      migen.net = targetscan_db[pick, ]
      migen.net = migen.net[order(migen.net$weighted.context...score, decreasing = FALSE), ]
      #      migen.net$weight=migen.net$weighted.context...score
      #      migen.net$gene=migen.net$Gene.Symbol
      
      weightPercentile.migene.thres = as.double(input$weightPercentilemigeneThres)
      weightScore.migene.thres = as.double(input$weightScoremigeneThres)
      pick = which(
        migen.net$weighted.context...score.percentile >= weightPercentile.migene.thres &
          migen.net$weighted.context...score <= weightScore.migene.thres
      )
      migen.net = migen.net[pick, ]
      
      pick = migen.net$Gene.Symbol %in% c(edges_from, edges_to) &
        migen.net$miRNA %in% c(edges_from, edges_to)
      migen.net = migen.net[pick, ]
      rownames(migen.net) = NULL
      
      #get table of cirmi or cirmi.net
      RNAhybrid_cirmi.net[c("circRNA", "weight")] <- NULL
      #RNAhybrid_cirmi.net[c("normID","circRNA","miRNA")]<-NULL
      miRanda_cirmi.net[c("normID", "circRNA", "miRNA", "weight")] <-
        NULL
      TargetScan_cirmi.net[c("normID", "circRNA", "miRNA")] <- NULL
      
      cirmi.net = RNAhybrid_cirmi.net
      cirmi.net = merge(cirmi.net, miRanda_cirmi.net, by = c("int", "int"))
      cirmi.net = merge(cirmi.net, TargetScan_cirmi.net, by = c("int", "int"))
      #ordering by Max_Score
      cirmi.net=cirmi.net[order(cirmi.net$Max_Score, decreasing=TRUE),]
      
      m = match(cirmi.net$normID, hsaCirData_small$normID)
      p1 = which(!is.na(m))
      if (length(p1) > 0)
        cirmi.net$int[p1] = hsaCirData_small$circRNA.ID[m[p1]]
      colnames(cirmi.net)[1] = "circRNA"
      x=colnames(cirmi.net)
      x=c(setdiff(x,"normID"),"normID") #put normID to the end
      cirmi.net=cirmi.net[,x]
      
      
      rownames(cirmi.net)=NULL
      rownames(migen.net)=NULL
      
      return(list(migen.net = migen.net, cirmi.net = cirmi.net))
      #return(cirmi.net)
    }

  
    output$cirmiTableOut <-
      renderDataTable(server = FALSE, {
        #add server = FALSE to renderDataTable to allow it to fetch the whole table instead of just the part displayed
        data = getIntTables()$cirmi.net
        datatable(
          data,
          extensions = 'Buttons',
          options = list(
            scrollX = TRUE,
            lengthMenu = c(5, 10, 15),
            paging = TRUE,
            searching = TRUE,
            fixedColumns = TRUE,
            autoWidth = FALSE,
            pageLength = 5,
            ordering = TRUE,
            dom = 'Bfrtip',
            columnDefs = list(list(width = '10%', targets = c(2))),
            #dom = 'tB', #disable search and paging
            #buttons = c('copy', 'csv', 'excel','pdf')
            buttons =
              list(
                'copy',
                #'print',
                list(
                  extend = 'csv',
                  filename = paste0('cir-mir_', Sys.Date())
                ),
                list(
                  extend = 'excel',
                  filename = paste0('cir-mir_', Sys.Date())
                ),
                list(
                  extend = 'pdf',
                  filename = paste0('cir-mir_', Sys.Date())
                )
              )
          )
        )# %>% formatStyle(columns = c(2), width='20px')

      })
    
    
    output$migenTableOut <-
      renderDT(server = FALSE, {
        #add server = FALSE to renderDataTable to allow it to fetch the whole table instead of just the part displayed
        data = getIntTables()$migen.net
        datatable(
          data,
          extensions = 'Buttons',
          options = list(
            scrollX = TRUE,
            lengthMenu = c(5, 10, 15),
            paging = TRUE,
            searching = TRUE,
            fixedColumns = TRUE,
            autoWidth = TRUE,
            pageLength = 5,
            ordering = TRUE,
            dom = 'Bfrtip',
            #dom = 'frtBip',
            #dom = 'tB', #disable search and paging
            #buttons = c('copy', 'csv', 'excel','pdf')
            buttons =
              list(
                'copy',
                #'print',
                list(
                  extend = 'csv',
                  filename = paste0('mir-mRNA_', Sys.Date())
                ),
                list(
                  extend = 'excel',
                  filename = paste0('mir-mRNA_', Sys.Date())
                ),
                list(
                  extend = 'pdf',
                  filename = paste0('mir-mRNA_', Sys.Date())
                )
              )
          )
        )
      })
    
  
    ##########
    # observeEvent(input$goButton,{
    #   showNotification("Building the circRNA-miRNA-mRNA interaction network, please wait until completed!")
    # })
    
    observeEvent(input$cbShowcircRNAseqlist, {
      if (input$cbShowcircRNAseqlist)
        shinyjs::show("circRNAseqlist")
      else
        shinyjs::hide("circRNAseqlist")
    })
    
    observeEvent(input$cbShowmiRNAList, {
      if (input$cbShowmiRNAList)
        shinyjs::show("miRNAList")
      else
        shinyjs::hide("miRNAList")
    })
    
    observeEvent(input$cbShowGeneList, {
      if (input$cbShowGeneList)
        shinyjs::show("geneList")
      else
        shinyjs::hide("geneList")
    })
    
    observeEvent(input$cbShowRBPList, {
      if (input$cbShowRBPList)
        shinyjs::show("RBPList")
      else
        shinyjs::hide("RBPList")
    })
    
    
    observeEvent(input$cbUseRNAhybrid, {
      if (input$cbUseRNAhybrid) {
        shinyjs::show("RNAhybrid_pvalThres")
        shinyjs::show("RNAhybrid_mfeThres")
      }
      else{
        shinyjs::hide("RNAhybrid_pvalThres")
        shinyjs::hide("RNAhybrid_mfeThres")
      }
    })
    
    observeEvent(input$cbUseMiRanda, {
      if (input$cbUseMiRanda)
        shinyjs::show("miRandaScoreThres")
      else
        shinyjs::hide("miRandaScoreThres")
    })
    
    observeEvent(input$cbUseTargetScan, {
      if (input$cbUseTargetScan) {
        shinyjs::show("TargetScan_6mer")
        shinyjs::show("TargetScan_7mer1a")
        shinyjs::show("TargetScan_7merm8")
        shinyjs::show("TargetScan_8mer1a")
      }
      else{
        shinyjs::hide("TargetScan_6mer")
        shinyjs::hide("TargetScan_7mer1a")
        shinyjs::hide("TargetScan_7merm8")
        shinyjs::hide("TargetScan_8mer1a")
      }
    })
    

    ##########
    output$selectedMiRNAs <- renderText({
      if (!is.null(input$circRNAnetwork_nodes)) {
        mynodes = input$circRNAnetwork_nodes
        #     if(!is.null(input$circRNAnetwork_selectedNodes)){
        #       mynodes=input$circRNAnetwork_selectedNodes
        res_group = sapply(mynodes, function(x)
          x["group"])
        res_label = sapply(mynodes, function(x)
          x["label"])
        pick = which(res_group == "miRNA")
        res = res_label[pick]
        res = paste(res, collapse = " ")
        return(res)
      }
    })
    
    output$selectedGenes <- renderText({
      if (!is.null(input$circRNAnetwork_nodes)) {
        mynodes = input$circRNAnetwork_nodes
        #     if(!is.null(input$circRNAnetwork_selectedNodes)){
        #       mynodes=input$circRNAnetwork_selectedNodes
        res_group = sapply(mynodes, function(x)
          x["group"])
        res_label = sapply(mynodes, function(x)
          x["label"])
        pick = which(res_group == "gene")
        res = res_label[pick]
        res = paste(res, collapse = " ")
        return(res)
      }
    })
    
    #output$showMiRNA <- renderDT({
    output$showMiRNA <- renderText({
      # #wrong, should from network
      # circInfo=circInfoData()
      # cirmi.net=circInfo$cirmi.net
      # migen.net=circInfo$migen.net
      # circRNAset=unique(as.character(cirmi.net$circRNA))
      # miRNAset=unique(as.character(cirmi.net$miRNA))
      
      net = mynetwork()
      nodes <- net$x$nodes
      #    edges <- net$x$edges
      #    edges_from=nodes$label[match(edges$from,nodes$id)] #circRNA
      #    edges_to=nodes$label[match(edges$to,nodes$id)] #miRNA
      #    node_group=nodes$group
      miRNAset = unique(nodes$label[nodes$group %in% "miRNA"])
      return(miRNAset)
      
    })
    
    #output$showGene <- renderDT({
    output$showGene <- renderText({
      # #wrong, should from network
      # circInfo=circInfoData()
      # cirmi.net=circInfo$cirmi.net
      # migen.net=circInfo$migen.net
      # geneset=unique(as.character(migen.net$gene))
      net = mynetwork()
      nodes <- net$x$nodes
      #    edges <- net$x$edges
      #    edges_from=nodes$label[match(edges$from,nodes$id)] #circRNA
      #    edges_to=nodes$label[match(edges$to,nodes$id)] #miRNA
      #    node_group=nodes$group
      geneset = unique(nodes$label[nodes$group %in% "gene"])
      return(geneset)
    })
    
    observe({
      input$goPWButton
      visNetworkProxy("circRNAnetwork") %>% visGetNodes()
      #visNetworkProxy("circRNAnetwork") %>%visGetSelectedNodes()
    })
    

    ############################
    output$reactomeAPI <- renderDataTable(server = FALSE,{
      geneset = NULL
      if (input$rbgset == "net") {
        #if (!input$cbgset){
        #get only genes from the network
        #visNetworkProxy("circRNAnetwork") %>%visGetNodes()
        if (!is.null(input$circRNAnetwork_nodes)) {
          mynodes = input$circRNAnetwork_nodes
          #     if(!is.null(input$circRNAnetwork_selectedNodes)){
          #       mynodes=input$circRNAnetwork_selectedNodes
          res_group = sapply(mynodes, function(x)
            x["group"])
          res_label = sapply(mynodes, function(x)
            x["label"])
          pick = which(res_group == "gene")
          geneset = res_label[pick]
          geneset = paste(geneset, collapse = " ")
        }
      } else{
        #NO, should from the network, not all genes
        #get all genes
        circInfo = circInfoData()
        cirmi.net = circInfo$cirmi.net
        migen.net = circInfo$migen.net
        geneset = unique(as.character(migen.net$gene))
        geneset = trimws(geneset)
      }
      
      if (is.null (geneset)) {
        return(DT::datatable(matrix(0, 1, 1)))
      } else {
        geneset_str = paste(geneset, collapse = " ")
        if (trimws(geneset_str) == "")
          return(DT::datatable(matrix(0, 1, 1)))
        ### generate a random ID
        ranID = floor(runif(1, min = 1, max = 1000000))
        ### export the pseudo-sequence to Json file
        json_fn = paste("reactomejson_", ranID, ".js", sep = "")
        cmd_curl = paste(
          'curl -X POST "https://reactome.org/AnalysisService/identifiers/projection?interactors=false&sortBy=ENTITIES_PVALUE&order=ASC&resource=TOTAL&pValue=1&includeDisease=true" -H "accept: application/json" -H "content-type: text/plain" -d "',
          geneset_str,
          '" > ',
          json_fn,
          sep = ""
        )
        
        system(cmd_curl)
        res <- fromJSON(file = json_fn
        )
   
        system(paste("rm ", json_fn, sep = ""))
        out = lapply(res$pathways, unlist)
        out = do.call(rbind, out)
        out = out[, -c(2, 4, 5, 7, 8, 14)]
        
        # Rename column "stId"
        colnames(out)[1] <- "Link"
        #cat (colnames (out))
        
        #creating hyperlink of stID to reactome page
        out[, c(1)] <-
          sapply(out[, c(1)], function(x)
            toString(tags$a(
              href = paste0('https://reactome.org/PathwayBrowser/#/', x, sep = ""), x
            )))
        
        #creating datatable
        out = DT::datatable(
          out,
          extensions = c('Buttons'),
          options = list(dom = 'frtBip',
                         scrollX = TRUE,
                         lengthMenu = c(5, 10, 15),
                         paging = TRUE,
                         searching = TRUE,
                         fixedColumns = TRUE,
                         autoWidth = TRUE,
                         pageLength = 5,
                         ordering = TRUE,
                         buttons = list(
                           list(
                             extend = 'csv',
                             filename = 'geneset_page',
                             text = 'Download this page',
                             exportOptions = list(modifier = list(page = "current"))
                           ),
                           list(
                             extend = 'csv',
                             filename = 'geneset_all',
                             text = 'Download all',
                             exportOptions = list(modifier = list(page = "all"))
                           )
                         ))
          ,
          escape = FALSE
        )
        return(out)
        datatable
      }
    })
    
    output$cirmiTargetScanPlots <- renderPlot({
      cirmi.net = getIntTables()$cirmi.net
        if (nrow(cirmi.net) > 0) {
          par(mfrow = c(2, 2), mar = c(2, 1, 2, 2) + 0.1)
          hist(cirmi.net$st_8mer1a, main = "TargetScan - 8mer1a", xlab = "")
          hist(cirmi.net$st_7merm8, main = "TargetScan - 7merm8", xlab = "")
          hist(cirmi.net$st_7mer1a, main = "TargetScan - 7mer1a", xlab = "")
          hist(cirmi.net$st_6mer, main = "TargetScan - 6mer", xlab = "")
        }
    })
    
    output$cirmiRNAhybridPlots <- renderPlot({
      cirmi.net = getIntTables()$cirmi.net
      if (nrow(cirmi.net) > 0) {
        par(mfrow = c(1, 2), mar = c(2, 1, 2, 2) + 0.1)
        hist(cirmi.net$pval, main = "RNAhybrid - pvalue", xlab = "")
        hist(cirmi.net$mfe, main = "RNAhybrid - mfe", xlab = "")
      }
    })

    output$cirmiMiRandaPlots <- renderPlot({
      cirmi.net = getIntTables()$cirmi.net
      if (nrow(cirmi.net) > 0) {
        par(mfrow = c(3, 2), mar = c(2, 1, 2, 2) + 0.1)
        hist(cirmi.net$Max_Score, main = "miRanda - Max_Score", xlab = "")
        hist(cirmi.net$Max_Energy, main = "miRanda - Max_Energy", xlab = "")
        hist(cirmi.net$Tot_Score, main = "miRanda - Tot_Score", xlab = "")
        hist(cirmi.net$Tot_Energy, main = "miRanda - Tot_Energy", xlab = "")
        hist(cirmi.net$siteNum, main = "miRanda - siteNum", xlab = "")

      }
    })
   
    output$mimRNAPlots <- renderPlot({
      migen.net = getIntTables()$migen.net

      if (nrow(migen.net) > 0) {
        migen_geneNum = table(migen.net$miRNA)
        migen_wcs = migen.net$weighted.context...score
        migen_wcsp = migen.net$weighted.context...score.percentile
        
        par(mfrow = c(2, 2), mar = c(2, 1, 2, 2) + 0.1)
        #plot(0,0,yaxt = "n",xaxt = "n",xlab = "",ylab = "",bty = "n",pch = "")
        hist(migen_geneNum, main = "Gene # of a miRNA", xlab = "")
        hist(migen_wcs, main = "Weighted context score", xlab = "")
        hist(migen_wcsp, main = "Score percentile", xlab = "")
      }
      
    })
    
    
    output$circRBPStatsPlots <- renderPlot({
      CircRBP_int.net=getCircRBPinfo(input)
      if (!is.null(CircRBP_int.net)){
        plot(0,0,yaxt = "n",xaxt = "n",xlab = "",ylab = "",bty = "n",pch = "")
        hist(CircRBP_int.net$SiteNum, main = "", xlab = "Site number")
      }
    })
    
    output$circRBPTableOut <-
      renderDT(server = FALSE, {
        #add server = FALSE to renderDataTable to allow it to fetch the whole table instead of just the part displayed
        data = getCircRBPinfo(input)
        if (!is.null(data)){
          data=data[,c("normID","circRNA", "RBP", "SiteNum")]
        }
        datatable(
          data,
          extensions = 'Buttons',
          options = list(
            scrollX = TRUE,
            lengthMenu = c(5, 10, 15),
            paging = TRUE,
            searching = TRUE,
            fixedColumns = TRUE,
            autoWidth = TRUE,
            pageLength = 5,
            ordering = TRUE,
            dom = 'Bfrtip',
            #dom = 'frtBip',
            #dom = 'tB', #disable search and paging
            #buttons = c('copy', 'csv', 'excel','pdf')
            buttons =
              list(
                'copy',
                #'print',
                list(
                  extend = 'csv',
                  filename = paste0('circRBP_', Sys.Date())
                ),
                list(
                  extend = 'excel',
                  filename = paste0('circRBP_', Sys.Date())
                ),
                list(
                  extend = 'pdf',
                  filename = paste0('circRBP_', Sys.Date())
                )
              )
          )
        )
      })
    
  }) #end of shinyServer



