#' Compound query
#' 
#' This function will perform a query of a given compound/drug in KEGG or PubChem and return available data based on queries to the REST APIs for either database. 
#' 
#' @param query_string A character string with a compound/drug name
#' @return A list of tibbles that contain general information, target information, pathways where compound is involved, diseases associated with compound/drug, common uses for a drug of interest, and KEGG drug groups for compound of interest.
#' 
#' @examples 
#' compound_query("aspirin")
compound_query <- function(query_string) {
  # make vectors to hold all results
  general <- vector(mode = "list", length = length(query_string))
  targets <- vector(mode = "list", length = length(query_string))
  pathways <- vector(mode = "list", length = length(query_string))
  diseases <- vector(mode = "list", length = length(query_string))
  uses <- vector(mode = "list", length = length(query_string))
  dgroups <- vector(mode = "list", length = length(query_string))
  
  for (i in seq(1, length(query_string))) {
    
    # first search kegg
    a1 <- kegg_query(query_string[i])
    if (!is.na(a1)) {
      a2 <- extract_info(a1)
      a3 <- property_extractor(a2$external_ids$pubchem_cid)
      
      # general properties
      general[[i]] <- tibble::tibble(name = query_string[i],
                                     kegg_id = a2$properties$drug_id,
                                     pubchem_cid = a2$external_ids$pubchem_cid,
                                     drugbank_id = a2$external_ids$drugbank_id,
                                     cas_rn = a2$external_ids$cas_id,
                                     formula = a2$properties$formula,
                                     exact_mass = a2$properties$exact_mass,
                                     molecular_weight = a2$properties$molecular_weight,
                                     chebi_id = a2$external_ids$chebi_id,
                                     canonical_smiles = a3$canonical_smiles,
                                     inchikey = a3$inchikey)
      
      # targets
      if (length(a2$targets$target_entrez) != 0) { 
        tmp <- AnnotationDbi::mapIds(x = org.Hs.eg.db, keys = a2$targets$target_entrez, 
                                     keytype = "ENTREZID", 
                                     column = "SYMBOL") 
        targets[[i]] <- tibble::tibble(name = query_string[i],
                                       kegg_id = a2$properties$drug_id,
                                       pubchem_cid = a2$external_ids$pubchem_cid,
                                       drugbank_id = a2$external_ids$drugbank_id,
                                       target_entrez = a2$targets$target_entrez,
                                       target_symbol = tmp)
        
      } else { 
        targets[[i]] <- tibble::tibble(name = query_string[i],
                                       kegg_id = a2$properties$drug_id,
                                       pubchem_cid = a2$external_ids$pubchem_cid,
                                       drugbank_id = a2$external_ids$drugbank_id,
                                       target_entrez = NA,
                                       target_symbol = NA)
        
      }
      
      
      # pathways
      if (length(a2$pathway$pathway_id) != 0) { 
        pathways[[i]] <- tibble::tibble(name = query_string[i],
                                        kegg_id = a2$properties$drug_id,
                                        pubchem_cid = a2$external_ids$pubchem_cid,
                                        drugbank_id = a2$external_ids$drugbank_id,
                                        pathway_id = a2$pathway$pathway_id,
                                        pathway_name = KEGG$pathway$pathway_name[match(pathway_id, KEGG$pathway$pathway_id)])  
      } else {
        pathways[[i]] <- tibble::tibble(name = query_string[i],
                                        kegg_id = a2$properties$drug_id,
                                        pubchem_cid = a2$external_ids$pubchem_cid,
                                        drugbank_id = a2$external_ids$drugbank_id,
                                        pathway_id = NA,
                                        pathway_name = NA)
      }
      
      
      
      # diseases
      if (length(a2$indication$drug_disease) != 0) {
        diseases[[i]] <- tibble::tibble(name = query_string[i],
                                        kegg_id = a2$properties$drug_id,
                                        pubchem_cid = a2$external_ids$pubchem_cid,
                                        drugbank_id = a2$external_ids$drugbank_id,
                                        disease_id = a2$indication$drug_disease,
                                        disease_name = KEGG$disease$disease_name[match(disease_id, KEGG$disease$disease_id)])  
      } else {
        diseases[[i]] <- tibble::tibble(name = query_string[i],
                                        kegg_id = a2$properties$drug_id,
                                        pubchem_cid = a2$external_ids$pubchem_cid,
                                        drugbank_id = a2$external_ids$drugbank_id,
                                        disease_id = NA,
                                        disease_name = NA)
      }
      
      
      # uses
      if (length(a2$efficacy$drug_efficacy) != 0) {
        uses[[i]] <- tibble::tibble(name = query_string[i],
                                    kegg_id = a2$properties$drug_id,
                                    pubchem_cid = a2$external_ids$pubchem_cid,
                                    drugbank_id = a2$external_ids$drugbank_id,
                                    drug_uses = a2$efficacy$drug_efficacy)  
      } else {
        uses[[i]] <- tibble::tibble(name = query_string[i],
                                    kegg_id = a2$properties$drug_id,
                                    pubchem_cid = a2$external_ids$pubchem_cid,
                                    drugbank_id = a2$external_ids$drugbank_id,
                                    drug_uses = NA)
      }
      
      
      # drug groups
      if (length(a2$group$drug_group) != 0) {
        dgroups[[i]] <- tibble::tibble(name = query_string[i],
                                       kegg_id = a2$properties$drug_id,
                                       pubchem_cid = a2$external_ids$pubchem_cid,
                                       drugbank_id = a2$external_ids$drugbank_id,
                                       group_id = a2$group$drug_group,
                                       group_name = KEGG$drug_group$group_name[match(a2$group$drug_group, KEGG$drug_group$group_id)])
        
      } else {
        dgroups[[i]] <- tibble::tibble(name = query_string[i],
                                       kegg_id = a2$properties$drug_id,
                                       pubchem_cid = a2$external_ids$pubchem_cid,
                                       drugbank_id = a2$external_ids$drugbank_id,
                                       group_id = NA,
                                       group_name = NA)
        
      }
    } else {
      # if you don't find anything, search pubchem
      b1 <- name2cid(query_string[i])
      if (!is.na(b1$pubchem_cid)) {
        b2 <- property_extractor(b1$pubchem_cid)
        b3 <- pubchem2cas(b1$pubchem_cid)
        if (!is.na(b3$cas_rn)) {
          b4 <- which(DRUGBANK$general$cas_number == b3$cas_rn)
          if (length(b4) != 0) {
            b5 <- DRUGBANK$general$primary_key[b4]
            # general properties
            general[[i]] <- tibble::tibble(name = query_string[i],
                                           kegg_id = NA,
                                           pubchem_cid = b3$pubchem_cid,
                                           drugbank_id = b5,
                                           cas_rn = b3$cas_rn,
                                           formula = b2$molecular_formula,
                                           exact_mass = NA,
                                           molecular_weight = b2$molecular_weight,
                                           chebi_id = DRUGBANK$external_ids[DRUGBANK$external_ids$parent_key == b5, ] %>% dplyr::filter(., resource == "ChEBI") %>% dplyr::select(., identifier) %>% as.numeric,
                                           canonical_smiles = b2$canonical_smiles,
                                           inchikey = b2$inchikey)
            
            # targets
            targets[[i]] <- tibble::tibble(name = query_string[i],
                                           kegg_id = NA,
                                           pubchem_cid = b3$pubchem_cid,
                                           drugbank_id = b5,
                                           target_entrez = NA,
                                           target_symbol = NA)
            
            # pathways
            pathways[[i]] <- tibble::tibble(name = query_string[i],
                                            kegg_id = NA,
                                            pubchem_cid = b3$pubchem_cid,
                                            drugbank_id = b5,
                                            pathway_id = NA,
                                            pathway_name = if(length(which(DRUGBANK$pathway$parent_key == b5)) != 0) {
                                              DRUGBANK$pathway$name[DRUGBANK$pathway$parent_key == b5] } else { NA })
            
            # diseases
            diseases[[i]] <- tibble::tibble(name = query_string[i],
                                            kegg_id = NA,
                                            pubchem_cid = b3$pubchem_cid,
                                            drugbank_id = b5,
                                            disease_id = NA,
                                            disease_name = NA)
            
            # uses
            uses[[i]] <- tibble::tibble(name = query_string[i],
                                        kegg_id = NA,
                                        pubchem_cid = b3$pubchem_cid,
                                        drugbank_id = b5,
                                        drug_uses = NA)
            
            # drug groups
            dgroups[[i]] <- tibble::tibble(name = query_string[i],
                                           kegg_id = NA,
                                           pubchem_cid = b3$pubchem_cid,
                                           drugbank_id = b5,
                                           group_id = NA,
                                           group_name = if (length(which(DRUGBANK$category$parent_key == b5)) != 0) { 
                                             tolower(DRUGBANK$category$category[DRUGBANK$category$parent_key == b5])
                                             } else { NA })
            
            
          } else { # <-- if it has no drugbank id but has cas
            # general properties
            general[[i]] <- tibble::tibble(name = query_string[i],
                                           kegg_id = NA,
                                           pubchem_cid = b3$pubchem_cid,
                                           drugbank_id = NA,
                                           cas_rn = b3$cas_rn,
                                           formula = b2$molecular_formula,
                                           exact_mass = NA,
                                           molecular_weight = b2$molecular_weight,
                                           chebi_id = NA,
                                           canonical_smiles = b2$canonical_smiles,
                                           inchikey = b2$inchikey)
            
            # targets
            targets[[i]] <- tibble::tibble(name = query_string[i],
                                           kegg_id = NA,
                                           pubchem_cid = b3$pubchem_cid,
                                           drugbank_id = NA,
                                           target_entrez = NA,
                                           target_symbol = NA)
            
            # pathways
            pathways[[i]] <- tibble::tibble(name = query_string[i],
                                            kegg_id = NA,
                                            pubchem_cid = b3$pubchem_cid,
                                            drugbank_id = NA,
                                            pathway_id = NA,
                                            pathway_name = NA)
            
            # diseases
            diseases[[i]] <- tibble::tibble(name = query_string[i],
                                            kegg_id = NA,
                                            pubchem_cid = b3$pubchem_cid,
                                            drugbank_id = NA,
                                            disease_id = NA,
                                            disease_name = NA)
            
            # uses
            uses[[i]] <- tibble::tibble(name = query_string[i],
                                        kegg_id = NA,
                                        pubchem_cid = b3$pubchem_cid,
                                        drugbank_id = NA,
                                        drug_uses = NA)
            
            # drug groups
            dgroups[[i]] <- tibble::tibble(name = query_string[i],
                                           kegg_id = NA,
                                           pubchem_cid = b3$pubchem_cid,
                                           drugbank_id = NA,
                                           group_id = NA,
                                           group_name = NA)
          }
        } else { # if it has no cas but has pubchem id
          # general properties
          general[[i]] <- tibble::tibble(name = query_string[i],
                                         kegg_id = NA,
                                         pubchem_cid = b3$pubchem_cid,
                                         drugbank_id = NA,
                                         cas_rn = NA,
                                         formula = b2$molecular_formula,
                                         exact_mass = NA,
                                         molecular_weight = b2$molecular_weight,
                                         chebi_id = NA,
                                         canonical_smiles = b2$canonical_smiles,
                                         inchikey = b2$inchikey)
          
          # targets
          targets[[i]] <- tibble::tibble(name = query_string[i],
                                         kegg_id = NA,
                                         pubchem_cid = b3$pubchem_cid,
                                         drugbank_id = NA,
                                         target_entrez = NA,
                                         target_symbol = NA)
          
          # pathways
          pathways[[i]] <- tibble::tibble(name = query_string[i],
                                          kegg_id = NA,
                                          pubchem_cid = b3$pubchem_cid,
                                          drugbank_id = NA,
                                          pathway_id = NA,
                                          pathway_name = NA)
          
          # diseases
          diseases[[i]] <- tibble::tibble(name = query_string[i],
                                          kegg_id = NA,
                                          pubchem_cid = b3$pubchem_cid,
                                          drugbank_id = NA,
                                          disease_id = NA,
                                          disease_name = NA)
          
          # uses
          uses[[i]] <- tibble::tibble(name = query_string[i],
                                      kegg_id = NA,
                                      pubchem_cid = b3$pubchem_cid,
                                      drugbank_id = NA,
                                      drug_uses = NA)
          
          # drug groups
          dgroups[[i]] <- tibble::tibble(name = query_string[i],
                                         kegg_id = NA,
                                         pubchem_cid = b3$pubchem_cid,
                                         drugbank_id = NA,
                                         group_id = NA,
                                         group_name = NA)
          }
      } else { # <-- has no pubchem cid
        # general properties
        general[[i]] <- tibble::tibble(name = query_string[i],
                                       kegg_id = NA,
                                       pubchem_cid = NA,
                                       drugbank_id = NA,
                                       cas_rn = NA,
                                       formula = NA,
                                       exact_mass = NA,
                                       molecular_weight = NA,
                                       chebi_id = NA,
                                       canonical_smiles = NA,
                                       inchikey = NA)
        
        # targets
        targets[[i]] <- tibble::tibble(name = query_string[i],
                                       kegg_id = NA,
                                       pubchem_cid = NA,
                                       drugbank_id = NA,
                                       target_entrez = NA,
                                       target_symbol = NA)
        
        # pathways
        pathways[[i]] <- tibble::tibble(name = query_string[i],
                                        kegg_id = NA,
                                        pubchem_cid = NA,
                                        drugbank_id = NA,
                                        pathway_id = NA,
                                        pathway_name = NA)
        
        # diseases
        diseases[[i]] <- tibble::tibble(name = query_string[i],
                                        kegg_id = NA,
                                        pubchem_cid = NA,
                                        drugbank_id = NA,
                                        disease_id = NA,
                                        disease_name = NA)
        
        # uses
        uses[[i]] <- tibble::tibble(name = query_string[i],
                                    kegg_id = NA,
                                    pubchem_cid = NA,
                                    drugbank_id = NA,
                                    drug_uses = NA)
        
        # drug groups
        dgroups[[i]] <- tibble::tibble(name = query_string[i],
                                       kegg_id = NA,
                                       pubchem_cid = NA,
                                       drugbank_id = NA,
                                       group_id = NA,
                                       group_name = NA)
      }
    } 
  }
  
  # pull everything together
  # i'll make a list of things
  general <- dplyr::bind_rows(general)
  targets <- dplyr::bind_rows(targets)
  pathways <- dplyr::bind_rows(pathways)
  diseases <- dplyr::bind_rows(diseases)
  uses <- dplyr::bind_rows(uses)
  dgroups <- dplyr::bind_rows(dgroups)
  
  # done
  return(list(general = general,
              targets = targets,
              pathways = pathways,
              diseases = diseases,
              drug_uses = uses,
              drug_groups = dgroups))
}