#' Parse DrugBank
#' 
#' Wrapper function to parse DrugBank XML file and generate a tibble for usage in DrugMineR. Built on the `dbparser` package. The data that is provided with `DrugMineR` was built using DrugBank v5.1.4.
#' 
#' @param drugbank_xml Full XML database file for DrugBank ([from here](https://www.drugbank.ca/releases/latest))
#' @return A tibble for DrugBank information
#' 
#' @example 
#' drugbank_parsed <- parse_drugbank("path/to/drugbank.xml")
parse_drugbank <- function(drugbank_xml) {
  
  message("Loading DrugBank data to memory...")
  dbparser::get_xml_db_rows(xml_db_name = drugbank_xml)
  
  message("Extracting general drug data...")
  drugs <- parse_drug()
  
  message("Extracting drug status...")
  drug_groups <- parse_drug_groups()
  
  # message("Extracting drug target actions...")
  # drug_targets_actions <- parse_drug_targets_actions()
  
  message("Extracting drug categories...")
  drug_categories <- dbparser::parse_drug_categories()
  # drug_classification <- dbparser::parse_drug_classification()
  
  message("Extracting drug pathways...")
  drug_pathway <- dbparser::parse_drug_pathway()
  # drug_targets <- dbparser::parse_drug_targets()
  
  message("Extracting drug external identifiers...")
  drug_extids <- dbparser::parse_drug_external_identifiers()
  # drug_synonyms <- dbparser::parse_drug_synonyms()
  
  # message("Extracting drug groups...")
  # drug_expprop <- dbparser::parse_drug_experimental_properties()
  # drug_doses <- dbparser::parse_drug_dosages()
  # drug_all <- dbparser::parse_drug_all()
  
  
  general <- drugs %>% dplyr::select(., primary_key,
                                     secondary_key, 
                                     type, 
                                     name, 
                                     description, 
                                     cas_number, 
                                     indication, 
                                     pharmacodynamics, 
                                     mechanism_of_action, 
                                     protein_binding, 
                                     toxicity)
  
  message("Putting it all together...")
  DRUGBANK <- list(general = general, 
                   category = drug_categories %>% dplyr::select(., parent_key, category),
                   pathway = drug_pathway %>% dplyr::select(., parent_key, name, smpdb_id, category),
                   status = drug_groups %>% dplyr::select(drugbank_id, group),
                   external_ids = drug_extids %>% as_tibble %>% dplyr::select(parent_key, resource, identifier))
  
  message("Done")
  return(DRUGBANK)
  
}


