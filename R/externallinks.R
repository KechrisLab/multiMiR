
#' Add External Database Link for Each of the multiMiR Result Entry
#' 
#' This is an internal multiMiR function that is not intended to be used
#' directly.  Please use \code{get.multimir}.
#' 
#' @param x table/dataset returned by multimir db
#' @param org Organism (see \code{get.multimir})
#' @export add.multimir.links
add.multimir.links <- function(x, org) {
    # To add external database link for each of the multiMiR result entry
    if (nrow(x) == 0) return(x)

    links <- rep(NA, nrow(x))
    db    <- as.character(unique(x$database))
    for (d in db) {
        m   <- which(x$database == d)
        mir <- as.character(x$mature_mirna_id[m])

        if (d == "mirecords") {
            # NOTE: need to resolve miRNA IDs with '*' in mirecords
            mir <- sub("-5p", "", mir)
            mir <- sub("-3p", "", mir)
            symbol <- as.character(x$target_symbol[m])
            if (org == "hsa") {
                s <- "species=Homo+sapiens"
            } else if (org == "mmu") {
                s <- "species=Mus+musculus"
            } else if (org == "rno") {
                s <- "species=Rattus+norvegicus"
            }
            links[m] <- paste("http://mirecords.biolead.org/interactions.php?",
                              s, "&mirna_acc=", mir,
                              "&targetgene_type=symbol&targetgene_info=", 
                              symbol, "&v=yes&search_int=Search", sep = "")
        } else if (d == "mirtarbase") {
            symbol   <- as.character(x$target_symbol[m])
            links[m] <- 
                paste("http://mirtarbase.mbc.nctu.edu.tw/php/search.php?org=",
                      org, "&mirnas=", mir, "&targets=", symbol, "&opt=adv", 
                      sep = "")
        } else if (d == "tarbase") {
            symbol   <- as.character(x$target_symbol[m])
            links[m] <- 
                paste("http://diana.imis.athena-innovation.gr/DianaTools/index.php?r=tarbase/index&mirnas=",
                      mir, "&genes=", symbol, sep = "")
        } else if (d == "mir2disease") {
            # NOTE: Can only search by miRNA, gene or disease alone - here use gene
            symbol   <- as.character(x$target_symbol[m])
            links[m] <- 
                paste("http://watson.compbio.iupui.edu:8080/miR2Disease/searchTarget.jsp?SearchUnit=target&SearchText=",
                      symbol, "&checkbox2=Causal&checkbox2=Unspecified", 
                      sep = "")
        } else if (d == "pharmaco_mir") {
            # NOTE: Links don't work

        } else if (d == "phenomir") {
            # NOTE: search by gene
            symbol   <- as.character(x$target_symbol[m])
            links[m] <- 
                paste("http://mips.helmholtz-muenchen.de/phenomir/main/list/searchform2?query=",
                      symbol, "&selectedview=mirs&searchtype=fuzzy", 
                      sep = "")
        } else if (d == "diana_microt") {
            ensembl <- as.character(x$target_ensembl[m])
            links[m] <- 
                paste("http://diana.imis.athena-innovation.gr/DianaTools/index.php?r=microT_CDS/results&genes=",
                      ensembl, "&mirnas=", mir, "&threshold=0", 
                      sep = "")
        } else if (d == "elmmo") {
            # NOTE: Need RefSeq accession for the gene - use miRNA only
            mir <- sub("-5p", "", mir)
            mir <- sub("-3p", "", mir)
            if (org == "hsa") {
                s <- "organism=hg"
            } else if (org == "mmu") {
                s <- "organism=mm"
            } else if (org == "rno") {
                s <- "organism=rn"
            }
            links[m] <- paste("http://www.mirz.unibas.ch/ElMMo3/?", s,
                              "&cellType=all&miRNAs[]=", mir,
                              "&predict=Predict+miRNAs+targets+!", 
                              sep = "")
        } else if (d == "microcosm") {
            mir <- sub("-5p", "", mir)
            mir <- sub("-3p", "", mir)
            if (org == "hsa") {
                s <- "genome_id=2964"
            } else if (org == "mmu") {
                s <- "genome_id=3876"
            } else if (org == "rno") {
                s <- "genome_id=5171"
            }
            symbol   <- as.character(x$target_symbol[m])
            links[m] <-
                paste("http://www.ebi.ac.uk/enright-srv/microcosm/cgi-bin/targets/v5/hit_list.pl?",
                      s, "&mirna_id=", mir, "&external_name=", symbol, 
                      sep = "")
        } else if (d == "miranda") {
            # NOTE: Could only search by gene or miRNA - use gene here
            if (org == "hsa") {
                s <- "organism=9606"
            } else if (org == "mmu") {
                s <- "organism=10090"
            } else if (org == "rno") {
                s <- "organism=10116"
            }
            symbol   <- as.character(x$target_symbol[m])
            links[m] <-
                paste("http://www.microrna.org/microrna/searchGenes.do?gene=",
                      symbol, "&", s, sep = "")
        } else if (d == "mirdb") {
            # NOTE: Could only search by gene or miRNA - use gene here
            if (org == "hsa") {
                s <- "species=Human"
            } else if (org == "mmu") {
                s <- "species=Mouse"
            } else if (org == "rno") {
                s <- "species=Rat"
            }
            symbol <- as.character(x$target_symbol[m])
            links[m] <- paste("http://mirdb.org/cgi-bin/search.cgi?", s,
                              "&searchType=gene&geneChoice=symbol&searchBox=", 
                              symbol, sep = "")
        } else if (d == "pictar") {
            # NOTE: Links don't work

        } else if (d == "pita") {
            mir <- sub("-5p", "", mir)
            mir <- sub("-3p", "", mir)
            if (org == "hsa") {
                s <- "Organism=Human"
            } else if (org == "mmu") {
                s <- "Organism=Mouse"
            }
            symbol <- as.character(x$target_symbol[m])
            links[m] <-
                paste("http://genie.weizmann.ac.il/cgi-bin/search_mir07_prediction.pl?",
                      s, "&microRNAs=", mir, "&Genes=", symbol,
                      "&MinimumSeed=7&AllowSingleGU=1&AllowSingleMismatch=1&MinConservation=0&FlankOption=0_0",
                      sep = "")
        } else if (d == "targetscan") {
            mir <- sub("-5p", "", mir)
            mir <- sub("-3p", "", mir)
            symbol <- as.character(x$target_symbol[m])
            if (org == "hsa") {
                links[m] <-
                    paste("http://www.targetscan.org/cgi-bin/targetscan/vert_61/targetscan.cgi?species=Human&gid=",
                          symbol, "&mirg=", mir, sep = "")
            } else if (org == "mmu") {
                links[m] <-
                    paste("http://www.targetscan.org/cgi-bin/targetscan/mmu_61/targetscan.cgi?species=Mouse&gid=",
                          symbol, "&mirg=", mir, sep = "")
            }
        }
    }

    x = data.frame(x, DB.link = links)
    return(x)
}

