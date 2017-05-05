defineAATable <- function(){
  charged <- data.frame(aa=c("R","K","H","D","E"),
                        charge=c("+2","+1","+1","-1","-1"))
  nonCharged <- data.frame(aa=c("A","N","C","Q","G","I","L","M","F","P","S","T","W","Y","V"),
                           charge="0")
  aaTable <- bind_rows(charged,nonCharged)
  return(aaTable)
}


annotateCharges <- function(x,
                              # refLabel="refAA",
                              # altLabel="altAA",
                              myTable=aaTable){
    # x = table to be annotated
    # refLabel = name of column in x that contains the original (ref) AA
    # altLabel = name of column in x that contains the substitute (alt) AA
    # aaTable = a data.frame that includes a column of AAs (named "aa") and a column of charges (named "charge")
    #####
    # returns 4 extra columns:
    # originalCharge = charge on refAA
    # substituteCharge = charge on altAA
    # deltaQ = charge change (quantitative)
    # chargeChange = charge change (logical)
    refCharges <- select(myTable,aa,refAACharge=charge)
    altCharges <- select(myTable,aa,altAACharge=charge)
  x %>%
      left_join(refCharges,by=c("refAA"="aa")) %>%
      left_join(altCharges,by=c("altAA"="aa")) %>%
      mutate(deltaQ=as.numeric(altAACharge)-as.numeric(refAACharge),
             chargeChange=(deltaQ!=0))  %>%
      return
  }  


# aaTable <- defineAATable()
# 
# hcmMis %>%
#   mutate(
#     refAA = Variant_protein %>% gsub("p.","",.) %>% stringr::str_extract("^[[:alpha:]]"),
#     altAA = Variant_protein %>% stringr::str_extract("[[:alpha:]]$")
#   ) %>%
#   annotateCharges
