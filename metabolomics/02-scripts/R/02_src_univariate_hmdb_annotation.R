library(readr)

metanr_packages <- function(){
  metr_pkgs <- c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "SSPA", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest", "RBGL", "edgeR", "fgsea", "devtools", "crmn")
  list_installed <- installed.packages()
  new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
  if(length(new_pkgs)!=0){if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install(new_pkgs)
    print(c(new_pkgs, " packages added..."))
  }
  
  if((length(new_pkgs)<1)){
    print("No new packages added...")
  }
}

library(httr)

##check hmdb ids of our list of compounds
check_metaboanalyst_db <- function(name_vector){
  name.vec <- name_vector
  toSend = list(queryList = name.vec, inputType = "name")
  
  call <- "https://www.xialab.ca/api/mapcompounds"
  
  query_results <<- httr::POST(call, body = toSend, encode = "json")
  
  query_results$status_code==200
  
  query_results_text <<- content(query_results, "text", encoding = "UTF-8")
  query_results_json <- rjson::fromJSON(query_results_text)
  query_results_table <- t(rbind.data.frame(query_results_json))
  rownames(query_results_table) <- query_results_table[,1]
  
  return(query_results_table)
}

##reading the annotation table

gcms <- read_tsv("03-data/raw_data/Metabolomics_IL10_primary_metabolism_id_To_taxonomy.txt")
lipids <- read_tsv("03-data/raw_data/Metabolomics_IL10_lipidomics_id_To_taxonomy.txt")
amines <- read_tsv("03-data/raw_data/Metabolomics_IL10_biogenic_amines._metabolism_d_To_taxonomy.txt")

#========= reading the json table of HMDB ===========#
library(fuzzyjoin)
hmdb <-  rjson::fromJSON(file = "03-data/databases/hmdb_full.json")

chem_compound <- unique(amines.df.tl1a.Named$name)

lookup_hmdb <- function(compound) {
  # Check if compound matches any 'name' or 'synonyms' in HMDB entries
  
  compound <<- compound
  
  if (nchar(compound) > 55){
    compound <- substr(compound, nchar(compound)-55, nchar(compound))
  }
  
  matching_entries <- lapply(hmdb, function(entry) {
    name_match <- grepl(compound, entry$name, ignore.case = TRUE)
    synonyms_match <- any(grepl(compound, entry$synonyms, ignore.case = TRUE))
    name_match | synonyms_match
  })
  
  # Extract the names of the matching entries
  matching_entry_names <- names(hmdb)[unlist(matching_entries)]
  
  # Return a list containing entry names and taxonomy information
  
  taxonomy<- lapply(matching_entry_names, function(entry_name) {
    list(input_name = compound,
      hmdb_id = entry_name, 
         name = c(hmdb[[entry_name]]$name, hmdb[[entry_name]]$synonyms),
         tax_kingdom = hmdb[[entry_name]]$tax_kingdom,
         tax_superclass = hmdb[[entry_name]]$tax_superclass,
         tax_class = hmdb[[entry_name]]$tax_class,
         tax_subclass = hmdb[[entry_name]]$tax_subclass,
         tax_direct_parent = hmdb[[entry_name]]$tax_direct_parent)
  })
  
  taxonomy <<- taxonomy
  
  if (length(taxonomy) == 0){
        taxonomy[[1]] <- list(input_name = compound,
                   hmdb_id = NA, 
                   name = NA,
                   tax_kingdom = NA,
                   tax_superclass = NA,
                   tax_class = NA,
                   tax_subclass = NA,
                   tax_direct_parent = NA)
        
        output <- data.table::rbindlist(taxonomy, fill = TRUE) %>% 
          mutate(string_dist = adist(name, compound))
  }
  
  else{
    output <- data.table::rbindlist(taxonomy, fill = TRUE) %>% 
      mutate(string_dist = adist(name, compound)) %>%
      filter(string_dist == min(string_dist)) %>%
      distinct()
  }
  
  return(output)
}

# test <- lookup_hmdb("(2-Methylpropyl)-3,6-dioxopiperazin-2-yl]propanoic acid")

# Apply the lookup function to each chemical compound
result <- lapply(chem_compound, lookup_hmdb)
