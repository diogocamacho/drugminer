#' Convert PubChem SID to PubChem CID
#' 
#' Function converts a PubChem SID into a PubChem CID using the NCBI's REST API for PubChem. 
#' 
#' @param query A PubChem SID
#' @return A tibble with PubChem SID and corresponding PubChem CID
#' 
#' @example 
#' sid2cid(query = "7847177")
sid2cid <- function(query) {
  
  # query is an sid. should have a check for this.
  
  cid <- vector(mode = "list", length = length(query))
  
  for (i in seq(1, length(cid))) {
    prolog <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    input <- paste0("/substance/sid")
    output <- "/cids/JSON?"
    qurl <- paste0(prolog, input, output, query[i])
    
    Sys.sleep(0.2)
    
    pubchem_content <- try(httr::content(POST(qurl,
                                              body = paste0("sid=", query[i]), 
                                              type = "text", 
                                              encoding = "UTF-8")),
                           silent = TRUE)
    # pubchem_content <- fromJSON(pubchem_content)
    
    if (names(pubchem_content) == "Fault") {
      cid[[i]] <- tibble::tibble(pubchem_sid = query[i], pubchem_cid = NA)
    } else {
      x1 <- unlist(pubchem_content)
      cid[[i]] <- tibble::tibble(pubchem_sid = query[i], pubchem_cid = unlist(x1)[2])      
    }
  }
  cid <- dplyr::bind_rows(cid)
  return(cid)
}