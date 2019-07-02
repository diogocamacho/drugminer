#' KEGG drug query
#' 
#' Extracts information from KEGG for a given drug of interest
kegg_query <- function(query_string) {
  
  id <- KEGG$drug$id[which(KEGG$drug$name == query_string)]
  
  if (length(id) > 1) id <- id[1]
  
  qurl <- paste0("http://rest.kegg.jp/get/", id)
  cont <- httr::content(POST(qurl), type = "text", encoding = "UTF-8")
  return(cont)
}