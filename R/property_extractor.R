#' Property extractor
#' 
#' Extracts molecular properties for a compound of interest based on its PubChem CID using the PubChem REST API.
#' 
#' @param pubchem_cid A PubChem CID
#' @return A tibble with CID, molecular formula, molecular weight, canonical SMILES, and InChIKey.
#' 
#' @examples 
#' property_extractor("5090")
property_extractor <- function(pubchem_cid) {
   
   props <- "MolecularFormula,MolecularWeight,CanonicalSMILES,InChIKey"
   
   if (length(pubchem_cid) > 1) {
     props_list <- vector(mode = "list", length = length(pubchem_cid))
     wait_time <- 0.2
     for (i in seq(1, length(props_list))) {
       prolog <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
       input <- "/compound/cid"
       output <- paste0("/property/", props, "/JSON?")
       qurl <- paste0(prolog, input, output, pubchem_cid[i])
       
       pubchem_content <- try(httr::content(POST(qurl,
                                                 body = paste0("cid=", pubchem_cid[i]), 
                                                 type = "text", 
                                                 encoding = "UTF-8")),
                              silent = TRUE)
       # pubchem_content <- fromJSON(pubchem_content)
       
       if (names(pubchem_content) == "Fault") {
         res[[i]] <- tibble::tibble(cid = pubchem_cid, 
                               molecular_formula = NA,
                               molecular_weight = NA,
                               canonical_smiles = NA,
                               inchikey = NA)
       } else {
         x1 <- unlist(pubchem_content)
         res[[i]] <- tibble::tibble(cid = pubchem_cid, 
                               molecular_formula = x1[2],
                               molecular_weight = x1[3],
                               canonical_smiles = x1[4],
                               inchikey = x1[5])
       }
       
       Sys.sleep(wait_time)
       
     }
     res <- dplyr::bind_rows(res)
   } else {
     prolog <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
     input <- "/compound/cid"
     output <- paste0("/property/", props, "/JSON?")
     qurl <- paste0(prolog, input, output, pubchem_cid)
     
     pubchem_content <- try(httr::content(POST(qurl,
                                               body = paste0("cid=", pubchem_cid), 
                                               type = "text", 
                                               encoding = "UTF-8")),
                            silent = TRUE)
     
     if (names(pubchem_content) == "Fault") {
       res <- tibble::tibble(cid = pubchem_cid, 
                             molecular_formula = NA,
                             molecular_weight = NA,
                             canonical_smiles = NA,
                             inchikey = NA)
     } else {
       x1 <- unlist(pubchem_content)
       res <- tibble::tibble(cid = pubchem_cid, 
                             molecular_formula = x1[2],
                             molecular_weight = x1[3],
                             canonical_smiles = x1[4],
                             inchikey = x1[5])
     }
   }
   return(res)
 }




