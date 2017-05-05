### Define motor domain activityal sites

defineMD <- function(myColours=c("darkblue","cadetblue","grey","orange")){
  rbind(
    #Nucleotide binding pocket
    data.frame(aaPos=126:134,
               activity="nucleotideBinding",
               label="atpBindingI",
               colour=myColours[1]),
    data.frame(aaPos=178:187,
               activity="nucleotideBinding",
               label="atpBindingII",
               colour=myColours[1]),
    data.frame(aaPos=238:245,
               activity="nucleotideBinding",
               label="atpBindingIII",
               colour=myColours[1]),
    data.frame(aaPos=269,
               activity="nucleotideBinding",
               label="atpBindingIV",
               colour=myColours[1]),
    data.frame(aaPos=461:471,
               activity="nucleotideBinding",
               label="atpBindingV",
               colour=myColours[1]),
    #Actin-myosin interface
    data.frame(aaPos=401:416,
               activity="actinMyosin",
               label="actinInterfaceI",
               colour=myColours[2]),
    data.frame(aaPos=527:556,
               activity="actinMyosin",
               label="actinInterfaceII",
               colour=myColours[2]),
    data.frame(aaPos=564:577,
               activity="actinMyosin",
               label="actinInterfaceIII",
               colour=myColours[2]),
    #Relay
    data.frame(aaPos=490:516,
               activity="relay",
               label="relay",
               colour=myColours[3]),
    #Converter
    data.frame(aaPos=710:777,
               activity="converter",
               label="converter",
               colour=myColours[4])
  )
}