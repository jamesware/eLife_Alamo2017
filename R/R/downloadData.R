downloadWalsh <- function(){
  system2("wget",
          args=c("--continue",
                 "--no-directories",
                 "--timestamping",
                 "--accept zip",
                 "http://www.nature.com/gim/journal/v19/n2/extref/gim201690x1.zip")
  )
  system2("mv",
          args=c("gim201690x1.zip",
                 "../data-raw/gim201690x1.zip"))
  system2("unzip",
          args=c("../data-raw/gim201690x1",
                 "*.xlsx",
                 "-d ../data-raw/"))
}
downloadHomburger <- function(){
  system2("wget",
          args=c("--continue",
                 "--no-directories",
                 "--timestamping",
                 "--accept xlsx",
                 "http://www.pnas.org/content/suppl/2016/05/27/1606950113.DCSupplemental/pnas.1606950113.sd02.xlsx")
  )
  system2("mv",
          args=c("pnas.1606950113.sd02.xlsx",
                 "../data-raw/pnas.1606950113.sd02.xlsx"))
}

