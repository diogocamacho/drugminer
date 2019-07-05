#' KEGG 
#' 
#' Data from KEGG. Data was parsed with the \link{\code(retrieve_data)} function.
#' 
#' @format A list of tibbles:
#' \describe{
#'   \item{drug}{A tibble with KEGG drug IDs and drug names}
#'   \item{drug_group}{A tibble with KEGG drug group IDs and drug group names}
#'   \item{disease}{A tibble with KEGG disease IDs and disease names}
#'   \item{pathway}{A tibble with KEGG pathway IDs and pathway names}
#'   \item{synonyms}{A tibble with KEGG drug IDs, drug names, and drug synonyms}
#'   \item{compound}{A tibble withe KEGG compound IDs and compound name}
#' }
#' 
#' @source \url{https://www.kegg.jp}
"KEGG"


#' DRUG
#' 
#' Data from KEGG, collected for drugs. Data was parsed with the \link{\code(retrieve_data)} function.
#' 
#' @format A list of tibbles:
#' \describe{
#'   \item{drug_group}{A tibble mapping KEGG drug IDs to KEGG drug group IDs, group name, and drug name}
#'   \item{drug_uses}{A tibble mapping KEGG drug IDs to known uses for the drugs (efficacy) and drug name}
#'   \item{disease_indication}{A tibble mapping KEGG drug IDs to KEGG disease IDs, indication name, and drug name}
#'   \item{targets}{A tibble mapping KEGG drug IDs to gene Entrez IDs, gene symbols, and drug name}
#'   \item{pathways}{A tibble mapping KEGG drug IDs to KEGG pathway IDs, pathway name, and drug name}
#'   \item{preparation}{A tibble mapping KEGG drug IDs to drug name and information if record is a known drug mixture (mixture == 1) or not (mixture == 0).}
#' }
#' 
#' @source \url{https://www.kegg.jp}
"DRUG"


#' DrugBank
#' 
#' Data from DrugBank.
#' 
#' @format A list of tibbles:
#' \describe{
#'   \item{general}{}
#'   \item{category}{}
#'   \item{pathway}{}
#'   \item{status}{}
#'   \item{external_ids}{}
#' }
#' 
#' @source \url{https://www.drugbank.ca}
"DRUGBANK"