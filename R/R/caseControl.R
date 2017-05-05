toggleCounts2 <- function(x,counts="cases"){
  #a helper function for caseControl
  #toggles between counting variant alleles & variant alleles for controls
  if(nrow(x)==0){return(0)}
  
  if(counts=="cases"){
    x %>% select(NFE_AC) %>% sum(na.rm=T) %>% return
  } else if(counts=="sites"){
    x %>% select(mut_coding) %>% unique %>% nrow
  }
}

caseControl <- function(target,
                        caseData=hcmMis,
                        controlData=exacVars,
                        nCases=6112, #MYH7/HCM
                        #exacVars %>% select(NFE_AN) %>% unlist %>% as.numeric() %>% median()
                        nControls=66728/2, #exac median AN / 2
                        counts="cases" ){# or "sites"){
  # inputs
  # a list of aa residues of interest (e.g. in a subdomain or interaction zone)
  # a table of variants in cases
  # determines whether region of interest is over-represented in cases
  caseVarsTarget = caseData %>%
    filter(aaPos %in% unlist(target)) %>%
    ungroup() %>%
    toggleCounts(counts)
  contVarsTarget = controlData %>%
    filter(mut_position %in% unlist(target)) %>%
    ungroup() %>%
    toggleCounts2(counts)
  
  data.frame(
    caseVarsTarget,
    # select(cases) %>%
    # sum,
    nCases,
    casePrev=(caseVarsTarget/nCases)  %>% signif(3),
    contVarsTarget,
    nControls,
    controlPrev=(contVarsTarget/nControls)  %>% signif(3)
  ) %>%
    mutate(pValueBinom=binom.test(x=caseVarsTarget,
                                  n=nCases,
                                  p=contVarsTarget/nControls,
                                  alternative="two")$p.value %>% signif(3),
           pValueFisher=fisher.test(matrix(c(caseVarsTarget,
                                             nCases-caseVarsTarget,
                                             contVarsTarget,
                                             nControls-contVarsTarget),2,2)
           )$p.value %>% signif(3),
           orFisher=fisher.test(matrix(c(caseVarsTarget,
                                         nCases-caseVarsTarget,
                                         contVarsTarget,
                                         nControls-contVarsTarget),2,2)
           )$estimate %>% signif(3),
           ciFisher=fisher.test(matrix(c(caseVarsTarget,
                                         nCases-caseVarsTarget,
                                         contVarsTarget,
                                         nControls-contVarsTarget),2,2)
           )$conf.int %>%
             signif(3) %>%
             paste(collapse="-"),
           ef=((orFisher-1)/orFisher)  %>% signif(3)
    ) %>%
    return()
}