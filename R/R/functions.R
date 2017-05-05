joinIhmMd <- function(x,
                      myClasses=c("4_likelyPath","5_Pathogenic")){
  x %>% ungroup %>%
    filter(consensusClass %in% myClasses) %>%
    left_join(interactions,by="aaPos") %>%
    left_join(mdFunctional,by="aaPos") %>%
    separate(consensusClass,into=c("classLabel","class"),sep="_") %>%
    arrange(desc(class)) %>%
    mutate(ihmLabel=paste(interaction," (",state,")",sep="")) %>%
    mutate(mdLabel=paste(activity," (",label.y,")",sep="")) %>%
    group_by(aaPos) %>%
    summarise(
      # dnaVariant=paste(unique(Variant_HGVS),collapse=","),
      aaVariant=paste(unique(Variant_protein),collapse=","),
      cases=sum(cases),
      # class,
      ihm=paste(unique(ihmLabel),collapse=","),
      ihmType=paste(unique(type),collapse=","),
      md=paste(unique(mdLabel),collapse=",")
    ) %>%
    ungroup %>%
    select (-aaPos) %>%
    (function(x){x[x=="NA (NA)"] <- ""; return(x)}) %>%
    (function(x){x[x=="NA"] <- ""; return(x)}) %>%
    return
}

#####

tidyMyInput <- function(x,sourceLabel="sourceLabel"){
  # colnames are written to first row. Move to names, and remove spaces/brackets
  names(x) <- x[1,] %>% 
    sub(" ","_",.) %>%
    sub("/",":",.) %>%
    gsub("[()]","",.)
  x <- x[-1,1:8] %>%
    mutate(source=sourceLabel) %>%
    tidyr::separate(col="Cases:Patients_tested",
                    sep="/",
                    into=c("cases","sequenced"))
  
  x$Clinical_classification[grep("VUS",x$Clinical_classification)] <- "VUS"
  x$Clinical_classification[grep("Class 3",x$Clinical_classification)] <- "VUS"
  x$Clinical_classification[x$Clinical_classification %in% c("Likely Pathogenic","Class 4-Likely pathogenic")]  <- "likelyPath"
  x$Clinical_classification[x$Clinical_classification %in% c("Pathogenic","Class 5-Certainly pathogenic")]  <- "Pathogenic"
  return(x)
}

#####

myMerge <- function(x){
  x <- x %>%
    group_by(Variant_HGVS,aaPos,Variant_protein) %>%
    summarise(
      cases=sum(as.numeric(cases)),
      records=n(),
      classification=paste(unique(Clinical_classification),collapse=","),
      source=paste(source,collapse=","),
      consensusClass=NA
    ) %>%
    mutate(source=gsub("OMGL,PLMM","both",source)) %>%
    mutate(source=gsub("PLMM,OMGL","both",source)) %>%
    arrange(desc(records),desc(classification)) %>%
    select(-records)
  x$consensusClass[x$classification %in% c("likelyPath")] <- "4_likelyPath"
  x$consensusClass[x$classification %in% c("VUS","VUS,likelyPath","likelyPath,VUS")] <- "3_VUS"
  x$consensusClass[x$classification %in% c("Pathogenic","Pathogenic,likelyPath","likelyPath,Pathogenic")] <- "5_Pathogenic"
  return(x)
}

#####

loadExAC <- function(){
  a <- read.delim("../data-raw/exac/Mutation_inExAC_sarcomericGenesOnly.txt") %>% tbl_df
  b <- read.delim("../data-raw/exac/GeneInfo.txt") %>% tbl_df
  c <- read.delim("../data-raw/exac/ExACPASSData.txt.gz") %>% tbl_df
  d <- read.delim("../data-raw/exac/ExACPopulations.txt") %>% tbl_df
  
    #retrieve all ExAC variants
    right_join(
      select(b,gene_id,gene_symbol,gene_canon_transcript),  #b = gene info
      a,  #b = variant info
      by = c("gene_id" = "mut_gene_id")
    ) %>%
    #filter for rare MYH7 missense
    filter(gene_symbol=="MYH7",mut_effect=="missense",mut_exac_frequency<0.0001) %>%
    select(mut_id,
           #gene_canon_transcript,
           #mut_chromosome,mut_chr_start_source,mut_chr_end_source,
           mut_coding,mut_protein,
           mut_position,
           mut_exac_count,mut_exac_frequency) %>%
    # annotate with NFE specific counts
    left_join(
      select(d,ep_id,ep_code) %>% # d = ExAC population codes
        right_join(c,by = c("ep_id" = "ed_ep_id")) %>%  #c = ExAC population-specific AC/AFs
        filter(ep_code=="NFE") %>%
        select(ed_mut_id,NFE_AN=ed_total,NFE_AC=ed_allele_count,NFE_AF=ed_frequency),
      by=c("mut_id" = "ed_mut_id")
    ) %>%
      return
}

#####

myCaseControlSummaryTable <- function(caseData=hcmMis,
                                      controlData=exacVars,
                                      nCases=6112, #MYH7/HCM
                                      #exacVars %>% select(NFE_AN) %>% unlist %>% as.numeric() %>% median()
                                      nControls=66728/2, #exac median AN / 2
                                      counts="cases" 
){
  allMyosin <- caseControl(myosinTotal,caseData,controlData,nCases,nControls,counts)
  rownames(allMyosin) <- "allMyosin"
  
  myosinHead <- caseControl(myosinHeadCluster,caseData,controlData,nCases,nControls,counts)
  rownames(myosinHead) <- "myosinHead"
  
  myosinNotHead <- myosinTotal[!myosinTotal %in% myosinHeadCluster]
  myosinNotHead <- caseControl(myosinNotHead,caseData,controlData,nCases,nControls,counts)
  rownames(myosinNotHead) <- "myosinNotHead"
  
  allIhmInteractions <- caseControl(group0,caseData,controlData,nCases,nControls,counts)
  rownames(allIhmInteractions) <- "allIhmInteractions"
  
  priming <- caseControl(group1,caseData,controlData,nCases,nControls,counts)
  rownames(priming) <- "priming"
  
  anchoring <- caseControl(group2,caseData,controlData,nCases,nControls,counts)
  rownames(anchoring) <- "anchoring"
  
  stabilising <- caseControl(group3,caseData,controlData,nCases,nControls,counts)
  rownames(stabilising) <- "stabilising"
  
  scaffolding <- caseControl(group4,caseData,controlData,nCases,nControls,counts)
  rownames(scaffolding) <- "scaffolding"
  
  converterSphere <- caseControl(homburgerSphere,caseData,controlData,nCases,nControls,counts)
  rownames(converterSphere) <- "converterSphere_preStroke"
  
  mesaSurface <- caseControl(homburgerSurface,caseData,controlData,nCases,nControls,counts)
  rownames(mesaSurface) <- "mesaSurface_preStroke"
  
  mdFunctionalSites <- caseControl(motorFunctional,caseData,controlData,nCases,nControls,counts)
  rownames(mdFunctionalSites) <- "mdFunctionalSites"
  
  
  rbind(
    allMyosin,myosinHead,myosinNotHead,
    allIhmInteractions,
    priming,anchoring,stabilising,scaffolding,
    # myosinHeadAndIhm,myosinHeadNotIhm,
    converterSphere,mesaSurface,
    mdFunctionalSites
  ) %>%
    select(-nCases,-nControls,-pValueFisher) %>%
    rename(or=orFisher,ci=ciFisher) %>%
    return
}

#####

compareHcmVsDcm <- function(myInteraction=interactions$aaPos,interactionLabel="allIhm"){
  data.frame(
    interaction=interactionLabel,
    dcmVarsOnIhm=dcmMis %>%
      filter(consensusClass %in% c("4_likelyPath","5_Pathogenic")) %>% 
      filter(aaPos %in% myInteraction) %>% 
      ungroup %>% select(Variant_HGVS) %>% nrow,
    dcmVarsTotal=dcmMis %>%
      filter(consensusClass %in% c("4_likelyPath","5_Pathogenic")) %>% 
      ungroup %>% select(Variant_HGVS) %>% nrow,
    hcmVarsOnIhm=hcmMis %>%
      filter(consensusClass %in% c("4_likelyPath","5_Pathogenic")) %>% 
      filter(aaPos %in% myInteraction) %>% 
      ungroup %>% select(Variant_HGVS) %>% nrow,
    hcmVarsTotal=hcmMis %>%
      filter(consensusClass %in% c("4_likelyPath","5_Pathogenic")) %>% 
      ungroup %>% select(Variant_HGVS) %>% nrow
  ) %>%
    mutate(dcmRateOnIhm=dcmVarsOnIhm/dcmVarsTotal,
           hcmRateOnIhm=hcmVarsOnIhm/hcmVarsTotal,
           rateRatio=dcmRateOnIhm/hcmRateOnIhm,
           binomP=binom.test(
             x=dcmVarsOnIhm,
             n=dcmVarsTotal,
             p=hcmRateOnIhm
           )$p.value
    ) %>%
    mutate(dcmRateOnIhm=signif(dcmRateOnIhm,2),
           hcmRateOnIhm=signif(hcmRateOnIhm,2),
           rateRatio=signif(rateRatio,2),
           binomP=signif(binomP,2)
    ) %>%
    return
}

#####

wrapForCC <- function(x){
  # wrap for case control
  # a wrapper that renames columns from hcmMis or dcmMis so that they can be passed as control to the caseControl() function.  Original function was written using column names for ExAC data.
  x %>% 
    rename(mut_position=aaPos,
           NFE_AC=cases) %>%
    return
}

#####

countOverlap <- function(x,
                         myClasses=c("4_likelyPath","5_Pathogenic")){
  x %>% ungroup %>%
    filter(consensusClass %in% myClasses) %>%
    left_join(interactions,by="aaPos") %>%
    left_join(mdFunctional,by="aaPos") %>%
    mutate(ihm=!is.na(interaction),
           md=!is.na(activity),
           label=ifelse(ihm+md==2,"both",
                        ifelse(ihm,"ihm",
                               ifelse(md,"md","none")))
    ) %>%
    group_by(label) %>%
    summarise(nVars=length(unique(Variant_HGVS))) %>%
    return
}

#####

myTabulate <- function(x,label="label"){
  x %>%
    group_by(Variant_protein,aaPos) %>%
    summarise(variants=n(),
              cases=sum(cases)) %>%
    group_by(aaPos) %>%
    summarise(variants=sum(variants),
              cases=sum(cases),
              substitutions=n()) %>%
    ungroup %>%
    summarise(label=label,
              variants=sum(variants),
              substitutions=sum(substitutions),
              residues=n(),
              cases=sum(cases)) %>%
    return()
}

myShaper <- function(x=hcmMis,interaction=myosinTotal){
  filter(x,aaPos %in% interaction)
}

myTabulator <- function(x=hcmMis){
  output <- rbind(
    x %>% myShaper(myosinTotal) %>% myTabulate("allMyosin"),
    x %>% myShaper(myosinHeadCluster) %>% myTabulate("myosinHead"),
    x %>% myShaper(homburgerSurface) %>% myTabulate("mesaSurface"),
    x %>% myShaper(motorFunctional) %>% myTabulate("mdFunctional"),
    x %>% myShaper(group0) %>% myTabulate("allIhmInteractions"),
    x %>% myShaper(group1) %>% myTabulate("priming"),
    x %>% myShaper(group2) %>% myTabulate("anchoring"),
    x %>% myShaper(group3) %>% myTabulate("stabilising"),
    x %>% myShaper(group4) %>% myTabulate("scaffolding")
  )
  for (i in unique(interactions$interaction)){
    nextTarget <- interactions %>%
      filter(interaction==i) %>%
      select(aaPos) %>%
      unlist
    nextRow <- x %>% myShaper(nextTarget) %>% myTabulate(i)
    output <- bind_rows(output,nextRow)
  }
  return(output)
}