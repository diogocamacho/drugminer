#' Retrieve data
#' 
#' Retrieves all data from KEGG using the REST API. This is a no-parameter function.
#' 
#' @return A list of tibbles: compounds, drugs, drug groups, pathways, and disease indications.
#' 
#' @examples 
#' kegg_data <- retrieve_data()
retrieve_data <- function(...) {
  
  ## all compounds in kegg
  message("Downloading KEGG compounds...")
  query_url <- "http://rest.kegg.jp/list/compound"
  compound_content <- content(POST(query_url), as = "text")
  kegg_compounds <- str_split(compound_content, "cpd:")[[1]]
  kegg_compounds <- kegg_compounds[-which(kegg_compounds == "")]
  tmp <- vector(mode = "list", length = length(kegg_compounds))
  for (i in seq(1, length(tmp))) {
    a1 <- str_split(kegg_compounds[i], "\t")[[1]] 
    a2 <- str_split(a1[2], ";")[[1]]
    a2 <- tolower(gsub("\n", "", sapply(a2, function(y) gsub("^ ", "", y))))
    tmp[[i]] <- tibble::tibble(id = a1[1], 
                               name = tolower(a2[1]),
                               synonyms = a2[1:length(a2)][-1])
  }
  compound <- dplyr::bind_rows(tmp)
  rm(tmp)
  
  
  ## all drugs in kegg
  message("Downloading KEGG drugs...")
  query_url <- "http://rest.kegg.jp/list/drug"
  compound_content <- content(POST(query_url), as = "text")
  kegg_drugs <- str_split(compound_content, "dr:")[[1]]
  kegg_drugs <- kegg_drugs[-which(kegg_drugs == "")]
  tmp <- vector(mode = "list", length = length(kegg_drugs))
  for (i in seq(1, length(tmp))) {
    a1 <- str_split(str_split(gsub(pattern = "\t", replacement = ":", kegg_drugs[i]), ";")[[1]][1], ":")[[1]]
    a2 <- gsub("\n", "", a1[2])
    a3 <- str_split(string = a2, pattern = " \\(")[[1]][1]
    tmp[[i]] <- tibble::tibble(id = a1[1], 
                               name = tolower(a3))
  }
  drug <- dplyr::bind_rows(tmp)
  rm(tmp)
  
  
  ## all drug groups
  message("Downloading KEGG drug groups...")
  query3 <- paste0("http://rest.kegg.jp/list/dgroup/")
  group_content <- httr::content(POST(query3), type = "text", encoding = "UTF-8")
  a1 <- str_extract_all(group_content, "(DG\\d{5}\\t[:print:]{1,}\\n)")[[1]]
  a2 <- str_split(a1, "\t", simplify = TRUE)
  a3 <- tolower(str_split(a2[, 2], "\n", simplify = TRUE)[, 1])
  drug_group <- tibble::tibble(group_id = a2[, 1],
                               group_name = a3)
  
  
  ## all disease -> drugs
  message("Downloading KEGG drug-disease indications...")
  query3 <- paste0("http://rest.kegg.jp/list/disease/")
  disease_content <- httr::content(POST(query3), type = "text", encoding = "UTF-8")
  a1 <- str_extract_all(disease_content, "(ds:H[:digit:]{5}\\t[:print:]{1,}\\n)")[[1]]
  a2 <- str_split(a1, "\t", simplify = TRUE)
  a3 <- str_split(a2[, 2], "\n", simplify = TRUE)[, 1]
  a4 <- str_remove(a2[, 1], pattern = "ds:")
  disease <- tibble::tibble(disease_id = a4,
                            disease_name = a3)
  
  ## all disease -> drugs
  message("Downloading KEGG pathways...")
  query3 <- paste0("http://rest.kegg.jp/list/pathway/")
  pathway_content <- httr::content(POST(query3), type = "text", encoding = "UTF-8")
  y2 <- str_split(string = gsub("\n", "", str_split(string = pathway_content, "path:")[[1]]), pattern = "\t", simplify = TRUE)
  y2 <- y2[-which(y2[, 1] == ""), ]
  pathway <- tibble::tibble(pathway_id = y2[, 1],
                            pathway_name = y2[, 2])
  
  return(list(compound = compound,
              drug = drug, 
              drug_group = drug_group, 
              disease = disease, 
              pathway = pathway))
}

