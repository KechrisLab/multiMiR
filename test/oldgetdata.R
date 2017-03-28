
# #' Get Validated microRNA-target Interactions from the multiMiR Package
# #' 
# #' This is an internal multiMiR function that is not intended to be used
# #' directly.  Please use \code{get.multimir}.
# #' 
# # @examples
# # 	 get.multimir.validated(mirna = "hsa-miR-18a-3p")
# # 
# get.multimir.validated <- function(url = NULL, org = "hsa", 
#                                    mirna = NULL, target = NULL) {
# 
#     if (!is.null(url)) deprecate_arg("url")
#     if (is.null(mirna) & is.null(target)) return(NULL)
# 
#     # get.multimir.by.table for each table in c() 
# 	rtn <- lapply(c("mirecords", "mirtarbase", "tarbase"), get.multimir.by.table, 
# 				  url = url, org = org, mirna = mirna, target = target)
# 	do.call(rbind, rtn)
# 
# }
# 
# #' Get Predicted microRNA-target Interactions from the multiMiR Package
# #' 
# #' This is an internal multiMiR function that is not intended to be used
# #' directly.  Please use \code{get.multimir}.
# #' 
# get.multimir.predicted <- function(url = NULL, org = "hsa", mirna = NULL, 
#                                    target = NULL, cutoff = NULL, 
# 								   cutoff.type = "p", site = "conserved") {
# 
#     if (!is.null(url)) deprecate_arg("url")
#     if (is.null(mirna) & is.null(target)) return(NULL)
# 
#     # get.multimir.by.table for each table in c() 
# 	rtn <- lapply(c("diana_microt", "elmmo", "microcosm", "miranda", "mirdb",
# 					"pictar", "pita", "targetscan"),
# 				  get.multimir.by.table, org = org, mirna = mirna, 
# 				  target = target, predicted.cutoff = cutoff,
# 				  predicted.cutoff.type = cutoff.type, predicted.site = site)
# 	do.call(rbind, rtn)
# 
# }
# 
# #' Get microRNA-disease/drug interactions from the multiMiR Package
# #' 
# #' This is an internal multiMiR function that is not intended to be used
# #' directly.  Please use \code{get.multimir}.
# #' 
# get.multimir.disease <- function(url = NULL,
#                                  org = "hsa",
#                                  mirna = NULL,
#                                  target = NULL,
#                                  disease.drug = NULL) {
# 
#     if (!is.null(url)) deprecate_arg("url")
# 	if (is.null(mirna) & is.null(target) & is.null(disease.drug)) return(NULL)
# 
#     # get.multimir.by.table for each table in c() 
# 	rtn <- lapply(c("mir2disease", "pharmaco_mir", "phenomir"),
# 				  get.multimir.by.table, org = org, mirna = mirna, 
# 				  target = target, disease.drug = disease.drug)
# 	do.call(rbind, rtn)
# 
# }
 
 
 
# #' Get microRNA/target Information from a Given Table in the multiMiR Package
# #' 
# #' This is an internal multiMiR function that is not intended to be used
# #' directly.  Please use \code{get.multimir}.
# #'
# # @param url                   PLACEHOLDER
# # @param org                   PLACEHOLDER
# # @param mirna                 PLACEHOLDER
# # @param target                PLACEHOLDER
# # @param table                 PLACEHOLDER
# # @param disease.drug          PLACEHOLDER
# # @param predicted.cutoff      PLACEHOLDER
# # @param predicted.cutoff.type PLACEHOLDER
# # @param predicted.site        PLACEHOLDER
# # @param mirna.table           PLACEHOLDER
# # @param target.table          PLACEHOLDER
# get.multimir.by.table <- function(table                 = NULL,
#                                   url                   = NULL,
#                                   org                   = c("hsa", "mmu", "rno"),
#                                   mirna                 = NULL,
#                                   target                = NULL,
#                                   disease.drug          = NULL,
#                                   predicted.cutoff      = NULL,
#                                   predicted.cutoff.type = "p",
#                                   predicted.site        = "conserved",
#                                   mirna.table           = "mirna",
#                                   target.table          = "target") {
#     # To get miRNA-target interactions in a given table
#     
#     # prepare query for validated target table
# 	if (table %in% c("mirecords", "mirtarbase", "tarbase")) {
# 		qry <- query_validated(table        = table,
#                                mirna        = mirna,
#                                target       = target,
# 							   org          = org,
#                                mirna.table  = mirna.table,
# 							   target.table = target.table)
# 	}
# 
# 	get.multimir.by.table(table                 = table,
#                           org                   = org,
# 						  mirna                 = mirna,
# 						  target                = target,
# 						  disease.drug          = disease.drug,
# 						  predicted.cutoff      = predicted.cutoff,
# 						  predicted.cutoff.type = predicted.cutoff.type,
# 						  predicted.site        = predicted.site)
#     
#     # prepare query for predicted target table
# 	if (table %in% c("diana_microt", "elmmo", "microcosm", "miranda", "mirdb",
# 					 "pictar", "pita", "targetscan")) {
# 		qry <- query_predicted(table                 = table,
#                                mirna                 = mirna,
#                                target                = target,
# 							   org                   = org,
#                                mirna.table           = mirna.table,
# 							   target.table          = target.table,
#                                predicted.cutoff      = predicted.cutoff,
#                                predicted.cutoff.type = predicted.cutoff.type,
#                                predicted.site        = predicted.site)
# 	}
# 
#     
#     # prepare query for miRNA-disease/drug table
# 	if (table %in% c("mir2disease", "pharmaco_mir", "phenomir")) {
# 
# 		qry <- query_disease(table        = table,
#                              mirna        = mirna,
#                              target       = target,
#                              org          = org,
# 							 mirna.table  = mirna.table,
#                              target.table = target.table,
# 							 disease.drug = disease.drug)
# 	}
# 
#     
#     # query the database
#     result <- search.multimir(query = qry)
#     if (!is.null(result)) 
#         result <- cbind(database = table, result)
#     if (table %in% c("mir2disease", "pharmaco_mir", "phenomir")) 
#         result <- unique(result)
# 
#     return(result)
# }
# 
# 
