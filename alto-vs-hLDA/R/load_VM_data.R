
load_VM_data <- function() {
  
  data <- 
    alto::vm_data$counts
  
  tax <- alto::vm_data$taxonomy %>% as_tibble()
  tax <- 
    tax %>% 
    mutate(
      ASV_name = 
        str_c(
          "asv", row_number(),
          str_replace_na(Genus, ""),
          str_replace_na(Species,"")
        ),
      ASV_name = 
        ASV_name %>% str_replace_all(" ","") %>% str_replace_all("/","")
    )
  
  colnames(data) <- tax$ASV_name
  
  data
  
}
