drug_synonyms <- function(drug_content) {
  a1 <- stringr::str_split(string = drug_content, pattern = "NAME")[[1]][2]
  a1 <- stringr::str_split(string = a1, pattern = "FORMULA")[[1]][1]
  a1 <- stringr::str_split(string = a1, pattern = "PRODUCT")[[1]][1]
  a1 <- stringr::str_split(string = a1, pattern = "EFFICACY")[[1]][1]
  a1 <- stringr::str_split(string = a1, pattern = "COMPONENT")[[1]][1]
  a1 <- stringr::str_split(string = a1, pattern = "REMARK")[[1]][1]
  a1 <- stringr::str_split(string = a1, pattern = "SOURCE")[[1]][1]
  a1 <- stringr::str_split(string = a1, pattern = "CLASS")[[1]][1]
  a1 <- stringr::str_split(string = a1, pattern = "SEQUENCE")[[1]][1]
  a1 <- stringr::str_split(string = a1, pattern = "\n")
  a1 <- unlist(a1)
  a1 <- stringr::str_remove_all(string = a1, pattern = "^ *")
  a1 <- gsub(pattern = ";", replacement = "", x = a1)
  a1 <- str_remove_all(string = a1, pattern = "\\([:print:]{1,}\\)")
  a1 <- gsub(pattern = " $", replacement = "", a1)
  a1 <- a1[a1 != ""]
  a1 <- unique(a1)
  a1 <- tolower(a1)
  return(a1)  
}