
# get.multimirold(
example3 <- get.multimir(org = "mmu", 
                         target  = "Gnb1",
                         table   = "miranda",
                         summary = TRUE,
                         predicted.cutoff      = 35,
                         predicted.cutoff.type = "p",
                         predicted.site        = "all")

example3 <- get.multimir(org = "mmu", 
                         target  = "Gnb1",
                         table   = "predicted",
                         summary = TRUE,
                         predicted.cutoff      = 35,
                         predicted.cutoff.type = "p",
                         predicted.site        = "all")

example2 = get.multimir(disease.drug='cisplatin', table='disease.drug')
example1 = get.multimir(mirna='hsa-miR-18a-3p', summary=TRUE)


assertWarning()
assertWarning()

list(org = "mmu", 
     target  = "Gnb1",
     table   = "predicted",
     summary = TRUE,
     predicted.cutoff      = 35,
     predicted.cutoff.type = "p",
     predicted.site        = "all",
     mirna = NULL,
     disease.drug = NULL) %>% list2env(env=globalenv())

