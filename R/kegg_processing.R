#' KEGG processing
#' 
#' Wrapper function for the \code{\link{retrieve_data}} function, building the `KEGG` and `DRUG` data frames if needed (`DrugMineR` comes pre-packaged with these data frames.)
#' 
#' @param save_dir User defined save directory for parsed data
#' @return KEGG and FRUG data frames.
# a wrapper to process all the data from kegg and write the data into somewhere
kegg_processing <- function(save_dir) {
  
  message("Processing data from KEGG. This will take a while.")
  message("")
  message("Retrieve all data...")
  KEGG <- retrieve_data()
  
  message("Annotate drugs...")
  prep <- vector(mode = "list", length = nrow(KEGG$drug))
  dg <- vector(mode = "list", length = nrow(KEGG$drug))
  uses <- vector(mode = "list", length = nrow(KEGG$drug))
  tgts <- vector(mode = "list", length = nrow(KEGG$drug))
  indic <- vector(mode = "list", length = nrow(KEGG$drug))
  pths <- vector(mode = "list", length = nrow(KEGG$drug))
  for (i in seq(1, nrow(KEGG$drug))) {
    message(paste("Drug:", KEGG$drug$name[i]))
    y1 <- kegg_query(query_string = KEGG$drug$name[i])
    d1 <- drug_info(y1)
    prep[[i]] <- d1$mixture
    dg[[i]] <- d1$group
    uses[[i]] <- d1$efficacy
    tgts[[i]] <- d1$targets
    indic[[i]] <- d1$indication
    pths[[i]] <- d1$pathway
    mps[[i]] <- d1$properties
  }
  
  message("Building data frames...")
  prep <- dplyr::bind_rows(prep)
  dg <- dplyr::bind_rows(dg)
  uses <- dplyr::bind_rows(uses)
  tgts <- dplyr::bind_rows(tgts)
  indic <- dplyr::bind_rows(indic)
  pths <- dplyr::bind_rows(pths)
  mps <- dplyr::bind_rows(mps)
  
  prep <- prep %>% 
    dplyr::mutate(., drug_name = KEGG$drug$name[match(prep$id, KEGG$drug$id)])
  
  dg <- dg %>% 
    dplyr::mutate(., group_name = KEGG$drug_group$group_name[match(dg$drug_group, KEGG$drug_group$group_id)]) %>% 
    dplyr::mutate(., drug_name = KEGG$drug$name[match(id, KEGG$drug$id)])
  
  uses <- uses %>% 
    dplyr::mutate(., drug_name = KEGG$drug$name[match(uses$id, KEGG$drug$id)])
  
  indic <- indic %>% 
    dplyr::mutate(., indication = KEGG$disease$disease_name[match(indic$drug_disease, KEGG$disease$disease_id)]) %>% 
    dplyr::mutate(., drug_name = KEGG$drug$name[match(id, KEGG$drug$id)])
  
  pths <- pths %>% 
    dplyr::mutate(., pathway_name = KEGG$pathway$pathway_name[match(pths$pathway_id, KEGG$pathway$pathway_id)]) %>% 
    dplyr::mutate(., drug_name = KEGG$drug$name[match(id, KEGG$drug$id)])

  tgts <- tgts %>% 
    dplyr::mutate(., gene_symbol = AnnotationDbi::mapIds(x = org.Hs.eg.db, keys = as.character(tgts$target_entrez), column = "SYMBOL", keytype = "ENTREZID")) %>% 
    dplyr::mutate(., drug_name = KEGG$drug$name[match(id, KEGG$drug$id)])
  
  DRUG <- list(drug_group = dg,
               drug_uses = uses,
               disease_indication = indic,
               targets = tgts,
               pathways = pths,
               preparation = prep)
  
    
  message("Saving...")
  save(file = paste0(save_dir, "/", "KEGG_data.RData"), KEGG, DRUG)
  
  message("Done.")
}