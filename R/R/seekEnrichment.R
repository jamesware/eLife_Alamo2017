toggleCounts <- function(x,counts="cases"){
  #a helper function for seekEnrichment
  #toggles between counting variant alleles & variant alleles
  if(nrow(x)==0){return(0)}  
  if(counts=="cases"){
    x %>% select(cases) %>% sum(na.rm=T) %>% return
  } else if(counts=="sites"){
    x %>% select(Variant_HGVS) %>% unique %>% nrow
  }
}

seekEnrichment <- function(interactions=myosinHeadCluster,
                           myClass="5_Pathogenic",
                           caseData=hcmMis,
                           controlData=NULL,
                           aaLength=1935,  # based on length of ENST00000355349, equivalent to P12883
                           counts="cases", 
                           classLabel=NULL){
  # inputs
  # a list of aa residues of interest (e.g. in a subdomain or interaction zone)
  # a table of variants in cases
  # determines whether region of interest is over-represented in cases
  caseVarsTarget = caseData %>%
    filter(consensusClass %in% myClass) %>%
    filter(aaPos %in% unlist(interactions)) %>%
    ungroup() %>%
    toggleCounts(counts)
  caseVarsTotal = caseData %>%
    filter(consensusClass %in% myClass) %>%
    ungroup() %>%
    toggleCounts(counts)
  aaLengthTarget = length(unique(unlist(interactions)))
  
  if(is.null(classLabel)){
    if(length(myClass)==1){myClassLabel=myClass} else {myClassLabel="unspecified"}
  } else {myClassLabel=classLabel}
  
  data.frame(
    class=myClassLabel,
    caseVarsTarget,
    caseVarsTotal,
    varsRateOnTarget=signif(caseVarsTarget/caseVarsTotal,3),
    aaLengthTarget,
    aaLengthTotal = aaLength,
    aaRateOnTarget=signif(aaLengthTarget/aaLength,3)
  ) %>%
    mutate(pValueBinom=binom.test(x=caseVarsTarget,
                                  n=caseVarsTotal,
                                  p=aaLengthTarget/aaLengthTotal,
                                  alternative="two")$p.value %>% signif(3),
           rateRatio=(varsRateOnTarget/aaRateOnTarget) %>% signif(3),
           pValueFisher=fisher.test(matrix(c(caseVarsTarget,
                                             caseVarsTotal-caseVarsTarget,
                                             aaLengthTarget,
                                             aaLengthTotal-aaLengthTarget),2,2)
           )$p.value  %>% signif(3),
           orFisher=fisher.test(matrix(c(caseVarsTarget,
                                         caseVarsTotal-caseVarsTarget,
                                         aaLengthTarget,
                                         aaLengthTotal-aaLengthTarget),2,2)
           )$estimate  %>% signif(3),
           ciFisher=fisher.test(matrix(c(caseVarsTarget,
                                         caseVarsTotal-caseVarsTarget,
                                         aaLengthTarget,
                                         aaLengthTotal-aaLengthTarget),2,2)
           )$conf.int %>%
             signif(3) %>%
             paste(collapse="-")
    ) %>%
    return()
}


seekEnrichmentTable <- function(interactions=bh_IHM,
                                caseData=hcmMis,
                                controlData=NULL,
                                aaLength=1935,
                                counts="cases"){
  rbind(
    seekEnrichment(interactions,myClass="5_Pathogenic",caseData,controlData,aaLength,counts),
    seekEnrichment(interactions,myClass="4_likelyPath",caseData,controlData,aaLength,counts),
    seekEnrichment(interactions,myClass="3_VUS",caseData,controlData,aaLength,counts)
  ) %>% return()
} 


formatRes <- function(x){
  x %>%
    select(aaPos) %>%
    unique %>%
    unlist %>%
    return
}