#' Extract info
#' 
#' Given the output of a POST request to KEGG REST API, returns information on either drug or compound queried using the \code{\link{drug_info}} or \code{\link{compound_info}}, respectively. The input to the function comes from the \code{\link{kegg_query}} output.
#' 
#' @param kegg_content REST API content for a KEGG query on a drug or compound.
#' @return Information on compound or drug, including chemical formula, external ids, disease indications, or pathways involved. Please refer to the \code{\link{drug_info}} or \code{\link{compound_info}} for details.
#' 
#' @example 
#' q <- kegg_query(query_string = "aspirin")
#' info <- extract_info(kegg_content = q)
extract_info <- function(kegg_content) {
  a1 <- str_split(kegg_content, "\n")[[1]][1] # <-- check if it's compound or drug
  if (is.na(str_extract(str_split(a1, "\n")[[1]][1], "(D\\d{5})"))) {
    message("Compound query. Extracting info.")
    res <- compound_info(kegg_content)
  } else {
    message("Drug query. Extracting info.")
    res <- drug_info(kegg_content)
  }
  return(res)
}