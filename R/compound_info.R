#' Compound info
#' 
#' Given a KEGG query obtained using the \link{\code(kegg_query)} function, returns information on a chemical compound. Returned elements are drug ids for compound (when available), pathways the compound is involved in, chemical formula, and external ids (CAS, PubChem, ChEBI, and DrugBank)  for compound.
#' 
#' @param kegg_content The content of a KEGG query
#' @return A tibble.
compound_info <- function(kegg_content) {
  
  # i'll need to carry this over
  compound_id <- str_extract(string = str_split(string = kegg_content, pattern = "NAME")[[1]][1], pattern = "(C\\d{5})")
  
  drug_map <- tibble::tibble(id = compound_id,
                            drug_id = unique(unlist(str_extract_all(string = kegg_content, pattern = "(D\\d{5})"))))
  
  form <- tibble::tibble(id = compound_id,
                         formula = str_remove_all(string = str_split(str_split(kegg_content, "FORMULA")[[1]][2], "\n")[[1]][1], pattern = " "))
  
  pathways <- tibble::tibble(id = compound_id,
                             pathway = str_extract_all(string = kegg_content, pattern = "(map\\d{5})")[[1]])
  
  # molecular properties
  drug_formula <- str_remove_all(str_split(str_split(drug_content, "FORMULA")[[1]][2], "\n")[[1]][1], " ")
  drug_mass <- as.numeric(str_remove_all(str_split(str_split(drug_content, "EXACT_MASS")[[1]][2], "\n")[[1]][1], " "))
  drug_weight <- as.numeric(str_remove_all(str_split(str_split(drug_content, "MOL_WEIGHT")[[1]][2], "\n")[[1]][1], " "))
  
  mprops <- tibble::tibble(formula = drug_formula, 
                           exact_mass = drug_mass, 
                           molecular_weight = drug_weight)
  
  # external ids
  cas_id <- str_split(str_split(kegg_content, "CAS: ")[[1]][2], "\n")[[1]][1]
  pubchem_sid <- str_split(str_split(kegg_content, "PubChem: ")[[1]][2], "\n")[[1]][1]
  chebi_id <- str_split(str_split(kegg_content, "ChEBI: ")[[1]][2], "\n")[[1]][1]
  drugbank_id <- str_split(str_split(kegg_content, "DrugBank: ")[[1]][2], "\n")[[1]][1]
  ext_ids <- tibble::tibble(id = compound_id,
                            cas_id = cas_id, 
                            pubchem_sid = pubchem_sid, 
                            chebi_id = chebi_id, 
                            drugbank_id = drugbank_id)
 
  return(list(drug_compound = drug_map,
              formula = form,
              pathways = pathways,
              external_ids = ext_ids, 
              properties = mprops))
}