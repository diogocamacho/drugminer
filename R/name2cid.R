#' Find PubChem CID from compound name
#' 
#' Function converts a name query into a PubChem CID using the NCBI's REST API for PubChem. 
#' 
#' @param query A compound name (as a string)
#' @return A tibble with name and corresponding PubChem CID
#' 
#' @example 
#' name2cid(query = "vioxx")
name2cid <- function(query) {
  
  if (length(query) > 1) {
    cid <- vector(mode = "list", length = length(query))
    wait_time <- 0.2
    for (i in seq(1, length(cid))) {
      prolog <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
      input <- "/compound/name"
      output <- "/cids/JSON?"
      qurl <- paste0(prolog, input, output, query[i])
      
      pubchem_content <- try(httr::content(POST(qurl,
                                                body = paste0("name=", query[i]), 
                                                type = "text", 
                                                encoding = "UTF-8")),
                             silent = TRUE)
      # pubchem_content <- fromJSON(pubchem_content)
      
      if (names(pubchem_content) == "Fault") {
        cid[[i]] <- tibble::tibble(name = tolower(query[i]), pubchem_cid = NA)
      } else {
        x1 <- unlist(pubchem_content)[[1]][1]
        cid[[i]] <- tibble::tibble(name = tolower(query[i]), pubchem_cid = x1)      
      }
      
      Sys.sleep(wait_time)
      
    }
    cid <- dplyr::bind_rows(cid)
    return(cid)
  } else {
    prolog <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    input <- "/compound/name"
    output <- "/cids/JSON?"
    qurl <- paste0(prolog, input, output, query)
    
    pubchem_content <- try(httr::content(POST(qurl,
                                              body = paste0("name=", query), 
                                              type = "text", 
                                              encoding = "UTF-8")),
                           silent = TRUE)
    
    if (names(pubchem_content) == "Fault") {
      cid <- tibble::tibble(name = tolower(query), pubchem_cid = NA)
    } else {
      x1 <- unlist(pubchem_content)[[1]][1]
      cid <- tibble::tibble(name = tolower(query), pubchem_cid = x1)      
    }
    return(cid)
  }
}