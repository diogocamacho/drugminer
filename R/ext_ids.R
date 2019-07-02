ext_ids <- function(drug_content) {
  cas_id <- str_split(str_split(drug_content, "CAS: ")[[1]][2], "\n")[[1]][1]
  pubchem_cid <- str_split(str_split(drug_content, "PubChem: ")[[1]][2], "\n")[[1]][1]
  chebi_id <- str_split(str_split(drug_content, "ChEBI: ")[[1]][2], "\n")[[1]][1]
  drugbank_id <- str_split(str_split(drug_content, "DrugBank: ")[[1]][2], "\n")[[1]][1]
  return(tibble::tibble(cas_id = cas_id, pubchem_cid = pubchem_cid, chebi_id = chebi_id, drugbank_id = drugbank_id))
}