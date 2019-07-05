#' KEGG drug query
#' 
#' Extracts information from KEGG for a given drug or compound of interest using the KEGG REST API. The POST request is done using the `httr` package. Given that this is a character vector that is unformated, it is likely unreadable and it is recommended that you run \link{\code(extract_info)} on the output of this function.
#' 
#' @param query_string A character string with compound or drug name.
#' @return The content of the POST request.
#' 
#' @examples 
#' Query a drug:
#' aspirin <- kegg_query("aspirin") 
#' 
#' Query a compound:
#' nad_syn <- kegg_query("dpnh")
kegg_query <- function(query_string) {
  
  x <- which(KEGG$drug$name == query_string)
  if (length(x) != 0) {
    id <- KEGG$drug$id[x]
  } else {
    x <- which(KEGG$synonyms$synonyms == query_string)
    if (length(x) != 0) {
      id <- KEGG$synonyms$id[x]
    } else {
      x <- which(KEGG$compound$name == query_string)
      if (length(x) != 0) {
        id <- KEGG$compound$compound_id[x]
      } else {
        x <- which(KEGG$compound$synonyms == query_string)
        if (length(x) != 0) {
          id <- KEGG$compound$compound_id[x]
        } else {
          stop("Can't find that drug in KEGG.")
        }
      }
    }
  }
  
  if (length(id) > 1) id <- id[1]
  qurl <- paste0("http://rest.kegg.jp/get/", id)
  cont <- httr::content(POST(qurl), type = "text", encoding = "UTF-8")
  return(cont)
}


