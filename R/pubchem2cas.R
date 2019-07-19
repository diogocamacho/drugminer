#' PubChem to CAS number 
#' 
#' Extracts a CAS number for a given PubChem CID.
#' 
#' @param query A PubChem CID
#' @return A tibble with PubChem CID and corresponding CAS number when found
#' 
#' @examples 
#' pubchem2cas("16038120")
pubchem2cas <- function(query) {
  message(paste("Query has", length(query), "elements."))
  res <- vector(mode = "list", length = length(query))
  for(i in seq(1, length(res))) {
    message(paste("Converting element", i))
    if (is.na(query[i])) {
      res[[i]] <- tibble::tibble(pubchem_cid = query[i],
                                 cas_rn = NA)
    } else {
      qurl <- paste0("http://cts.fiehnlab.ucdavis.edu/service/convert/PubChem CID/CAS/", query[i])
      qurl <- URLencode(qurl)
      qcont <- httr::content(GET(qurl))
      res[[i]] <- tibble::tibble(pubchem_cid = query[i],
                                 cas_rn = unlist(qcont)[4])
    }
  }
  res <- dplyr::bind_rows(res)
  return(res)
}