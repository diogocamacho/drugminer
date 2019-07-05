drug_info <- function(drug_content) {
  
  # i'll need to carry this over
  drug_id <- str_extract(string = str_split(string = drug_content, pattern = "NAME")[[1]][1], pattern = "(D\\d{5})")
  
  
  # is it a mixture?
  y1 <- str_split(string = drug_content, pattern = "NAME")[[1]][1]
  y2 <- grep("mixture", y1, ignore.case = TRUE)
  if (length(y2) == 0) {
    drug_mixture <- tibble::tibble(id = drug_id,
                                   mixture = 0)  
  } else {
    drug_mixture <- tibble::tibble(id = drug_id,
                                   mixture = 1)
  }
  
  
  # drug group
  y1 <- str_extract_all(string = drug_content, pattern = "(DG\\d{5})")[[1]]
  if (length(y1) != 0) {
    drug_group <- tibble::tibble(id = drug_id,
                                 drug_group = unique(y1))
  } else {
    drug_group <- tibble::tibble(id = drug_id,
                                 drug_group = NA)
  }
  
  
  # drug efficacy
  y1 <- str_split(str_split(drug_content, "EFFICACY")[[1]][2], "\n")[[1]][1]
  if (!is.na(y1)) {
    drug_efficacy <- str_remove(string = y1, pattern = "^ *")
    drug_efficacy <- unique(tolower(str_remove_all(string = str_split(drug_efficacy, ",")[[1]], pattern = "^ ")))
    drug_efficacy <- tibble::tibble(id = drug_id,
                                    drug_efficacy = drug_efficacy)
  } else {
    drug_efficacy <- tibble::tibble(id = drug_id,
                                    drug_efficacy = NA)
  }
  
  
  # drug targets
  # y1 <- str_split(str_split(drug_content, "TARGET")[[1]][2], "\n")[[1]][1]
  y1 <- str_extract_all(string = drug_content, pattern = "(HSA:\\d{1,})")[[1]]
  y1 <- gsub("HSA:", "", y1)
  drug_targets <- tibble::tibble(id = drug_id,
                                 target_entrez = y1)
  
  # if (!is.na(y1)) {
  #   drug_targets <- unique(str_split(str_remove(string = y1, pattern = "^ *"), " ")[[1]][1])
  #   drug_targets <- tibble::tibble(id = drug_id,
  #                                  drug_targets = drug_targets)
  # } else {
  #   drug_targets <- tibble::tibble(id = drug_id,
  #                                  drug_targets = NA)
  # }
  
  
  # disease indication
  y1 <- str_split(drug_content, "DISEASE")[[1]][2]
  if (!is.na(y1)) {
    drug_disease <- str_extract_all(string = y1, pattern = "(DS:[A-Z]\\d{5})")[[1]]
    drug_disease <- unique(gsub(pattern = "DS:", replacement = "", x = drug_disease))
    drug_disease <- tibble::tibble(id = drug_id,
                                   drug_disease = drug_disease)
  } else {
    drug_disease <- tibble::tibble(id = drug_id,
                                   drug_disease = NA)
  }
  
  # affected pathways
  y1 <- str_extract_all(string = drug_content, pattern = "(hsa\\d{5})")[[1]]
  y1 <- gsub("hsa", "map", y1)
  drug_pathways <- tibble::tibble(id = drug_id,
                                 pathway_id = y1)
  
  
  return(list(mixture = drug_mixture,
              group = drug_group,
              efficacy = drug_efficacy,
              targets = drug_targets,
              indication = drug_disease,
              pathway = drug_pathways))
}