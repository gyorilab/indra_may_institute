install.packages("tidyverse")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("MSstats")
# The following initializes usage of Bioc devel
# BiocManager::install(version='devel')
# BiocManager::install("MSstatsBioNet")
install.packages('devtools')
BiocManager::install('BiocStyle')
devtools::install_github("Vitek-Lab/MSstatsBioNet", build_vignettes = TRUE)
library(tidyverse)
library(MSstatsBioNet)
library(MSstats)

# Talus 2025 Dataset - Starter Code
input = data.table::fread("Talus-2025/model.csv")

## Doxorubicin
filtered_doxo = input %>% filter(Label == "Doxo-DMSO")
annotated_df_doxo = annotateProteinInfoFromIndra(filtered_doxo, "Uniprot_Mnemonic")
filtered_df_doxo = annotated_df_doxo %>% filter(is.na(issue))
edges_doxo = getPathwaysFromIndra(filtered_df_doxo, main_target = "doxorubicin", target_type = "Drug")

## Aclarubicin or aclacinomycin A 
filtered_acla = input %>% filter(Label == "Acla-DMSO")
annotated_df_acla = annotateProteinInfoFromIndra(filtered_acla, "Uniprot_Mnemonic")
filtered_df_acla = annotated_df_acla %>% filter(is.na(issue))
edges_acla = getPathwaysFromIndra(filtered_df_acla, main_target = "74619@CHEBI", target_type = "Drug")



# CPTAC analysis starter code - CDK7 to CDK1 matches with INDRA on phosphorylation of T161

## Import and filter data
input = data.table::fread("CPTAC/phospho_adjusted_PTM_model.csv")
adjusted_input = input %>% filter(Adjusted == TRUE)

## Clean up protein identifiers and PTM annotations
adjusted_input$GlobalProtein <- sub(".*\\|", "", adjusted_input$GlobalProtein)
adjusted_input$PTM = adjusted_input$Protein
adjusted_input$PTM = sub("^[^_]*_[^_]*_", "", adjusted_input$PTM)
adjusted_input$Protein = adjusted_input$GlobalProtein

## Filter and summarize PTM data
filtered = adjusted_input %>% filter(adj.pvalue < 0.001)
grouped_df = filtered %>% group_by(Protein) %>% summarise(adj.pvalue = median(adj.pvalue), log2FC = median(log2FC), issue = NA, ptms = paste(PTM, collapse="; "))

## Annotate protein information
annotated_df = annotateProteinInfoFromIndra(grouped_df, "Uniprot_Mnemonic")

## Filter by protein function
functional_df = annotated_df %>% filter(IsKinase == TRUE)

## Query subnetwork from INDRA
subnetwork = getSubnetworkFromIndra(
  functional_df, 
  statement_types = c("Phosphorylation"),
  evidence_count_cutoff = 10,
)

## Add PTM information to nodes
subnetwork$nodes <- subnetwork$nodes %>%
  left_join(functional_df %>% select(Protein, ptms), by = c("id" = "Protein"))

## Visualize networks
visualizeNetworks(
  subnetwork$nodes, 
  subnetwork$edges, 
  node_label_column = "hgncName", 
  logfcCutoff = 0.01,
  pvalueCutoff = 0.05
)


## Get pathways for clear drivers
main_driver = "VHL_HUMAN" # Feel free to replace with a different protein
filtered_pathways = adjusted_input %>% filter(adj.pvalue < 0.001 | str_detect(GlobalProtein, main_driver))
grouped_df_pathways = filtered_pathways %>% group_by(Protein) %>% summarise(adj.pvalue = median(adj.pvalue), log2FC = median(log2FC), issue = NA, ptms = paste(PTM, collapse="; "))
annotated_df_pathways = annotateProteinInfoFromIndra(grouped_df_pathways, "Uniprot_Mnemonic")
edges = getPathwaysFromIndra(annotated_df_pathways, main_target = main_driver)
edges <- edges %>%
  left_join(annotated_df_pathways %>% filter(adj.pvalue < 0.05) %>% select(Protein, adj.pvalue, ptms), by = c("target" = "Protein"))

