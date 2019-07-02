biophys_properties <- function(drug_content) {
  # formula
  drug_formula <- str_remove_all(str_split(str_split(drug_content, "FORMULA")[[1]][2], "\n")[[1]][1], " ")
  
  # exact mass
  drug_mass <- as.numeric(str_remove_all(str_split(str_split(drug_content, "EXACT_MASS")[[1]][2], "\n")[[1]][1], " "))
  
  # molecular weight
  drug_weight <- as.numeric(str_remove_all(str_split(str_split(drug_content, "MOL_WEIGHT")[[1]][2], "\n")[[1]][1], " "))
  
  return(tibble::tibble(formula = drug_formula, exact_mass = drug_mass, molecular_weight = drug_weight))
}