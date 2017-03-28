# 
# 
# 
# x <- outobj$data$validated
# 
# 
# hyperlinks <- function(.table, org) {
#     db <- unique(.table$database)
#     if (db > 1) stop("invalid table; # of databases > 1")
# 
#     org <- translate_org_for_url(org)
# 
#     mir     <- as.character(.table$mature_mirna_id)
#     symbol  <- as.character(.table$target_symbol)
#     if (db == "diana_microt") {
#         ensembl <- as.character(.table$target_ensembl) 
#     } else ensembl <- NULL
# 
# }
# 
# 
# translate_org_for_url <- function(.table, org) {
# 
#     if (.table == "mirecords") {
#         switch(org, hsa = "Homo+sapiens", mmu = "Mus+musculus",
#                rno = "Rattus+norvegicus")
#     } else if (.table == "elmmo") {
#         switch(org, hsa = "hg", mmu = "mm", rno = "rn")
#     } else if (.table == "microcosm") {
#         switch(org, hsa = "2964", mmu = "3876", rno = "5171")
#     } else if (.table == "miranda") { 
#         switch(org, hsa = "9606", mmu = "10090", rno = "10116")
#     } else if (.table == "mirdb") { 
#         switch(org, hsa = "Human", mmu = "Mouse", rno = "Rat")
#     } else if (.table == "pita") {
#         switch(org, hsa = "Human", mmu = "Mouse")
#     } else if (.table == "targetscan") {
#         switch(org, 
#                hsa = list(abbrev = "vert", cap = "Human"), 
#                mmu = list(abbrev = "mmu", cap = "Mouse"))
#     } else {
#         org 
#     }
# }
# 
# trim_mir <- function(.table, x) {
#     x <- as.character(x)
#     ifelse(.table %in% c("mirecords", "elmmo", "microcosm", "pita", "targetscan"),
#            sub("-5p|-3p", "", x), 
#            x)
# }
# 
# # mirecords 
# #  s      <- translated org
# #  mir    <- translate x$mature_mirna_id
# #  symbol <- translate x$target_symbol
# 
#     function(.table, org) {}
# 
# ## VALIDATED -- org, mir, symbol
# url_validated <- function(.url, org, mir, symbol) {
#     org     <- translate_org_for_url(.table, org)
#     mir     <- as.character(x$mature_mirna_id[m])
#     symbol  <- as.character(x$target_symbol[m])
#     ensembl <- as.character(x$target_ensembl[m]) # only diana_microt
# 
#     sprintf(.url, org, as.character(mir), as.character(symbol))
# 
# }
# 
# 
# url_mirecords <- function(.table, org) {
#     s <- translate_org_for_url(org)
#     mir <- trim
#     paste0("http://mirecords.biolead.org/interactions.php?species=", 
#            s, "&mirna_acc=", 
#            trim_mir(.table, .table$mature_mirna_id), 
#            "&targetgene_type=symbol&targetgene_info=", 
#            as.character(.table$target_symbol), "&v=yes&search_int=Search")
# }
#     mirtarbase   <- sprintf(paste0("http://mirtarbase.mbc.nctu.edu.tw/php/search.php",
#                                    "?org=%s&mirnas=%s&targets=%s&opt=adv"), 
#                             org, mir, symbol)
#     tarbase      <- sprintf(paste0("http://diana.imis.athena-innovation.gr/DianaTools/index.php?",
#                                    "r=tarbase/index&mirnas=%s&genes=%s"), 
#                             mir, symbol)
# ## DISEASEDRUG -- symbol
#     mir2disease  <- sprintf(paste0("http://watson.compbio.iupui.edu:8080/miR2Disease/searchTarget.jsp",
#                                    "?SearchUnit=target&SearchText=%s",
#                                    "&checkbox2=Causal&checkbox2=Unspecified"),
#                             symbol)
#     # pharmaco_mir -- links don't work
#     phenomir     <- sprintf(paste0("http://mips.helmholtz-muenchen.de/phenomir/main/list/searchform2?query=%s",
#                                    "&selectedview=mirs&searchtype=fuzzy"), 
#                             symbol)
# ## PREDICTED -- org, mir, symbol, ensemble
#     diana_microt <- sprintf(paste0("http://diana.imis.athena-innovation.gr/DianaTools/index.php",
#                                    "?r=microT_CDS/results&genes=%s&mirnas=%s&threshold=0"), 
#                             ensembl, mir)
#     elmmo <- sprintf(paste0("http://www.mirz.unibas.ch/ElMMo3/?organism=%s",
#                             "&cellType=all&miRNAs[]=%s",
#                             "&predict=Predict+miRNAs+targets+!"), 
#                      s_elmmo, mir)
#     microcosm <- sprintf(paste0("http://www.ebi.ac.uk/enright-srv/microcosm/",
#                                 "cgi-bin/targets/v5/hit_list.pl?genome_id=%s",
#                                 "&mirna_id=%s&external_name=%s"), 
#                          s_microcosm, mir, symbol)
#     miranda <- c(symbol, s)
#     miranda <- paste0("http://www.microrna.org/microrna/searchGenes.do?gene=%s",
#                       symbol, "&organism=", symbol, s_miranda)
#     mirdb <- c(s, symbol)
#     mirdb <- paste0("http://mirdb.org/cgi-bin/search.cgi?species=", s,
#                    "&searchType=gene&geneChoice=symbol&searchBox=", symbol)
#     # pictar -- links don't work
#     pita <- c(s, mir, symbol)
#     pita <- paste0("http://genie.weizmann.ac.il/cgi-bin/search_mir07_prediction.pl?Organism=",
#                   s, "&microRNAs=", mir, "&Genes=", symbol,
#                   "&MinimumSeed=7&AllowSingleGU=1&AllowSingleMismatch=1&MinConservation=0&",
#                   "FlankOption=0_0")
# 
# 
# url_targetscan <- function(.table, org) {
#     org <- translate_org_for_url(.table, org)
#     paste0("http://www.targetscan.org/cgi-bin/targetscan/", org$abbrev,
#            "_61/targetscan.cgi?species=", org$cap, "Human&gid=",
#            as.character(.table$target_symbol), "&mirg=", 
#            trim_mir(.table, .table$mature_mirna_id))
# }
# 
# 
# 
