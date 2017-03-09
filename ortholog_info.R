# Libraries ---------------------------------------------------------------
library(magrittr)
library(stringr)
library(dtplyr)
library(dplyr)
library(data.table)
load('Human.info.Rdata')
load('Mouse.info.Rdata')
# Prepare data from NCBI and MGI------------------------------------------------------
homo.ncbi <-
  read.csv(
    'homologene.data',
    header = F,
    sep = '\t',
    stringsAsFactors = F
  ) %>% select(V1, V2, V3, V4) %>% as.data.table()
hum.homo.ncbi <-
  homo.ncbi[V2 == 9606, .(V1, V3)] %>% `names<-`(c('homoid', 'Entrez.hum'))
mus.homo.ncbi <-
  homo.ncbi[V2 == 10090, .(V1, V3)] %>% `names<-`(c('homoid', 'Entrez.mus'))
homo.pairs.ncbi <- full_join(hum.homo.ncbi, mus.homo.ncbi) %>% as.data.table() #NCBI

homo.mgi <- read.csv('HOM_MouseHumanSequence.rpt',
                     sep = '\t',
                     stringsAsFactors = F) %>% as.data.table()
hum.homo.mgi <- homo.mgi[NCBI.Taxon.ID==9606,.(HomoloGene.ID,
                                               EntrezGene.ID,
                                               SWISS_PROT.IDs)] %>% `names<-`(c('homoid',
                                                                                'Entrez.hum',
                                                                                'Accession.hum'))
mus.homo.mgi <- homo.mgi[NCBI.Taxon.ID==10090,.(HomoloGene.ID,
                                               EntrezGene.ID,
                                               SWISS_PROT.IDs)] %>% `names<-`(c('homoid',
                                                                                'Entrez.mus',
                                                                                'Accession.mus'))
homo.pairs.mgi <- full_join(hum.homo.mgi,mus.homo.mgi) %>% as.data.table() #MGI


orth.pairs.ens <- read.csv('ensembl_orthlogpair.txt',
                     sep = '\t',
                     stringsAsFactors = F) %>% `names<-`(c('Ensembl.hum','Ensembl.mus')) %>% as.data.table()

# Merge -------------------------------------------------------------------
homo.pairs.ncbi[,homoid:=NULL]
homo.pairs.ncbi.pair <- homo.pairs.ncbi[!is.na(Entrez.hum) & !is.na(Entrez.mus)]

homo.pairs.mgi[,homoid:=NULL]
homo.pairs.mgi.pair <- homo.pairs.mgi[!is.na(Entrez.hum) & !is.na(Entrez.mus)]

homo.pairs <- full_join(homo.pairs.ncbi.pair, homo.pairs.mgi.pair) #Homolog pairs combined form MGI and NCBI.

pair.map <- homo.pairs %>%
  full_join(hum.info, by = 'Entrez.hum') %>% full_join(orth.pairs.ens) %>%  filter(!duplicated(.)) %>% as.data.table()

blank2na <- function(vec){
  vec[!str_detect(vec,'\\w')] <- NA
  vec
}

pair.map[str_detect(Accession.hum.x,','), Accession.hum.x := NA]
pair.map[str_detect(Accession.mus,','), Accession.mus:= NA]
pair.map <- pair.map[,lapply(.SD,blank2na)]
pair.map[,`:=`( Accession.hum = coalesce(Accession.hum.y,Accession.hum.x),
                   Accession.hum.x = NULL,
                   Accession.hum.y = NULL)] # Fix Accession.hum.
setkey(pair.map,NULL)
pair.map <- unique(pair.map)



pair.map <- pair.map %>% 
  left_join(mus.info[!is.na(Entrez.mus)],
            by = 'Entrez.mus') %>% as.data.table()
pair.map[,`:=`( Accession.mus = coalesce(Accession.mus.y,Accession.mus.x),
                Ensembl.mus = coalesce(Ensembl.mus.x,Ensembl.mus.y),
                   Accession.mus.x = NULL,
                   Accession.mus.y = NULL,
                   Ensembl.mus.x = NULL,
                   Ensembl.mus.y = NULL)] # Merge by Entrez.

pair.map <- pair.map %>% 
  left_join(mus.info[!is.na(Ensembl.mus)],
            by = 'Ensembl.mus') %>% as.data.table()
pair.map[,`:=`( Accession.mus = coalesce(Accession.mus.y,Accession.mus.x),
                Symbol.mus = coalesce(Symbol.mus.y,Symbol.mus.x),
                Entrez.mus = coalesce(Entrez.mus.y,Entrez.mus.x),
                Accession.mus.x = NULL,
                Accession.mus.y = NULL,
                Symbol.mus.x = NULL,
                Symbol.mus.y = NULL,
                Entrez.mus.x = NULL,
                Entrez.mus.y = NULL)]
# Merge by Ensembl.



setkey(pair.map,NULL)
pair.map <- unique(pair.map)
pair.map.pair <- pair.map[!is.na(Accession.hum) & !is.na(Accession.mus)]
# orthlog.map.hum <- orthlog.map[!is.na(Accession.hum) & is.na(Accession.mus)]
# orthlog.map.mus <- orthlog.map[is.na(Accession.hum) & !is.na(Accession.mus)]

# Pickup the unpaired -----------------------------------------------------

unpaired.hum <-
  hum.info[!is.na(Accession.hum) & !Accession.hum %in% pair.map.pair$Accession.hum][!is.na(Symbol.hum)]
unpaired.mus <-
  mus.info[!is.na(Accession.mus) & !Accession.mus %in% pair.map.pair$Accession.mus][!is.na(Symbol.mus)]

unpaired.mus[,Symbol.hum := toupper(Symbol.mus)]
manual.paired <- inner_join(unpaired.hum,unpaired.mus)

ortholog.pair <- bind_rows(pair.map.pair,manual.paired) %>% filter(!duplicated(.)) %>% select(Ensembl.hum,
                                                                                              Entrez.hum,
                                                                                              Symbol.hum,
                                                                                              Accession.hum,
                                                                                              Ensembl.mus,
                                                                                              Entrez.mus,
                                                                                              Symbol.mus,
                                                                                              Accession.mus)
ortholog.hum <- pair.map[!Accession.hum %in% pair.map.pair$Accession.hum & is.na(Accession.mus) & !is.na(Accession.hum)][!is.na(
  Entrez.mus)|!is.na(Ensembl.mus)] %>% select(Ensembl.hum,
                                  Entrez.hum,
                                  Symbol.hum,
                                  Accession.hum,
                                  Ensembl.mus,
                                  Entrez.mus,
                                  Symbol.mus,
                                  Accession.mus)
ortholog.mus <- pair.map[!Accession.mus %in% pair.map.pair$Accession.mus & is.na(Accession.hum) & !is.na(Accession.mus)][!is.na(
  Entrez.hum)|!is.na(Ensembl.hum)] %>% select(Ensembl.hum,
                                  Entrez.hum,
                                  Symbol.hum,
                                  Accession.hum,
                                  Ensembl.mus,
                                  Entrez.mus,
                                  Symbol.mus,
                                  Accession.mus)


save(ortholog.pair,file = 'Ortholog.pair.Rdata')
save(ortholog.hum,file = 'Ortholog.hum.Rdata')
save(ortholog.mus,file = 'Ortholog.mus.Rdata')

# End ---------------------------------------------------------------------

