# -*- coding: UTF-8 -*-
#' Find diseases from the VCF file
#'
#' This function finds diseases from a VCF file by comparing the variants in the VCF file with those in various databases.
#' @param file the name of the VCF file.
#' @return A data frame containing variant information from various databases
#' @importFrom data.table fread setDT
#' @import sqldf
#' @import dplyr
#' @importFrom stats setNames
#' @importFrom utils data install.packages write.table
#' @export
#'

DisVar <- function(file = "file_name.vcf", GWASdb = TRUE, GRASP = TRUE, GWASCat = TRUE, GAD = TRUE, JohnO = TRUE, p_value = 1e-7) {
  # Load required libraries
  required_packages <- c("sqldf", "data.table", "dplyr")
  for (package in required_packages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package)
      suppressPackageStartupMessages(library(package, character.only = TRUE))
    }
  }

  # Validate inputs
  if (!file.exists(file)) {
    stop("File does not exist")
  }
  if (!grepl(".vcf$", file, ignore.case = TRUE)) {
    stop("Input file must be a VCF file")
  }

  #Check database parameter
  databases <- list(
    "GWASdb" = GWASdb,
    "GRASP" = GRASP,
    "GWASCat" = GWASCat,
    "GAD" = GAD,
    "JohnO" = JohnO)

  if (!any(unlist(databases))) {
    stop("Please select at least one database")
  }


  #Read files
  cat("\nReading files...\n")
  GWASdb_GRCh38 <- V1 <- P_value <- Rsid <- Chr <- V4 <- V5 <- Ref <- Alt <- Gwas_trait <- Gene <- Variant_type <- V6 <- V7 <- V8 <- V2 <- GRASP_GRCh38 <-GWAS_catalog_GRCh38 <- JnO_GRCh38 <- Disease <- Confident <- NULL
  GWASdb_GRCh38 <- DisVar::GWASdb_GRCh38
  GRASP_GRCh38 <- DisVar::GRASP_GRCh38
  GWAS_catalog_GRCh38 <- DisVar::GWAS_catalog_GRCh38
  GAD_GRCh38 <- DisVar::GAD_GRCh38
  JnO_GRCh38 <- DisVar::JnO_GRCh38
  suppressWarnings(variant_data <- fread(file, sep = "\t", stringsAsFactors=FALSE, showProgress=TRUE, header =TRUE))
  colnames(variant_data) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8")

  . <- NULL
  chr <- c()
  allele_db_list <- c()
  allele_sample_list <- c()
  variant_data$V1<-gsub("chr","",as.character(variant_data$V1))
  variant_data$V2 <- as.integer(variant_data$V2)
  cat("Reading files...DONE\n")

  #Search databases
  cat("Searching...\n")

  result_df <- NULL

  for (db in names(databases)) {
    if (databases[[db]]) {
      op_df <- NULL
      switch(db,
             "GWASdb" = {
               op_df <- setDT(GWASdb_GRCh38)[setDT(variant_data), on = .(Position = V2, Chr = V1)][P_value < p_value, .(Rsid, Chr, Position, V4, V5, Ref, Alt, P_value, Gwas_trait, Gene, Variant_type, V6, V7, V8)]
               op_df$DB <- "GWASdb"
             },
             "GRASP" = {
               op_df <- setDT(GRASP_GRCh38)[setDT(variant_data), on = .(Position = V2, Chr = V1)][P_value < p_value, .(Rsid, Chr, Position, V4, V5, Ref, Alt, P_value, Gwas_trait, Gene, Variant_type, V6, V7, V8)]
               op_df$DB <- "GRASP"
             },
             "GWASCat" = {
               op_df <- setDT(GWAS_catalog_GRCh38)[setDT(variant_data), on = .(Position = V2, Chr = V1)][P_value < p_value, .(Rsid, Chr, Position, V4, V5, Ref, Alt, P_value, Gwas_trait, Gene, Variant_type, V6, V7, V8)]
               op_df$DB <- "GWAS Catalog"
             },
             "GAD" = {
               GAD_GRCh38$P_value <- 0
               op_df <- setDT(GAD_GRCh38)[setDT(variant_data), on = .(Position = V2, Chr = V1)][P_value < p_value, .(Rsid , Chr , Position, V4, V5, Ref, Alt, P_value, Gwas_trait , Gene , Variant_type, V6, V7, V8)]
               op_df$P_value <- NA
               op_df$DB <- "GADCDC"
             },
             "JohnO" = {
               op_df <- setDT(JnO_GRCh38)[setDT(variant_data), on = .(Position = V2, Chr = V1)][P_value < p_value, .(Rsid, Chr, Position, V4, V5, Ref, Alt, P_value, Gwas_trait, Gene, Variant_type, V6, V7, V8)]
               op_df$DB <- "Johnson and O'Donnell"
             }
      )
      if (nrow(op_df) > 0) {
        if (is.null(result_df)) {
          result_df <- op_df
        } else {
          result_df <- rbind(result_df, op_df, fill = TRUE)
        }
      }
    }
  }

  cat("Searching...DONE\n")

  cat("Processing results...\n")

  # If no variant is found, Exit program
  if (is.null(result_df) || nrow(result_df) == 0) {
    stop("No variant was found in the selected databases")
    # EXIT command
  }


  #create results in a table
  align_df <- setNames(data.frame(matrix(ncol = 13, nrow = nrow(result_df))), c("Disease", "Chrom", "Position", "Gene", "Variant ID", "Variant Type", "Allele Sample", "Allele DB", "Confident", "DB", "Qual", "Filter", "Info"))

  align_df["Disease"] <- sqldf('SELECT Gwas_trait FROM result_df')
  align_df["Chrom"] <- sqldf('SELECT Chr FROM result_df')
  align_df["Position"] <- sqldf('SELECT Position FROM result_df')
  align_df["Gene"] <- sqldf('SELECT Gene FROM result_df')
  align_df["Variant ID"] <- sqldf('SELECT Rsid FROM result_df')
  align_df["Variant Type"] <- sqldf('SELECT Variant_type FROM result_df')

  result_df$Ref <- as.character(result_df$Ref)
  result_df$Alt <- as.character(result_df$Alt)
  for (i in 1:nrow(result_df))
  {
    ref_db <- result_df[i,6]
    alt_db <- result_df[i,7]
    if (is.na(ref_db)) {allele_db_list <- append(allele_db_list,NA)}
    else
    {
      allele_db = paste(ref_db, alt_db, sep = ">")
      allele_db_list <- append(allele_db_list,allele_db)
    }
  }
  align_df["Allele DB"] <- allele_db_list

  for (i in 1:nrow(result_df))
  {
    ref_sample <- result_df[i,4]
    alt_sample <- result_df[i,5]
    allele_sample = paste(ref_sample, alt_sample, sep = ">")
    allele_sample_list <- append(allele_sample_list,allele_sample)
  }
  align_df["Allele Sample"] <- allele_sample_list

  align_df["Confident"] <- sqldf('SELECT P_value FROM result_df')
  align_df["DB"] <- sqldf('SELECT DB FROM result_df')
  align_df["Qual"] <- sqldf('SELECT V6 FROM result_df')
  align_df["Filter"] <- sqldf('SELECT V7 FROM result_df')
  align_df["Info"] <- sqldf('SELECT V8 FROM result_df')

  aligned_df <- align_df %>% arrange(Disease, Confident)
  aligned_df$Disease <- as.character(aligned_df$Disease)
  aligned_df$Disease[duplicated(aligned_df$Disease)] <- ''
  colnames(aligned_df)[9] <- "P-value"
  cat("Processing results...DONE\n")
  cat("Generating result file...\n")
  output_file <- sub(".vcf", "_DisVar.txt", file)
  write.table(aligned_df, file = output_file, quote = FALSE, sep = '\t', row.names = FALSE)
  cat("Generating result file...DONE\n")
  cat("The output file is saved as:", output_file, "in the directory:", getwd(), "\n")
}

#' GAD_GRCh38 Dataset
#'
#' This dataset contains GWAS information from the GAD database.
#'
#' @name GAD_GRCh38
#' @format A data frame with columns Rsid, Gene, Gwas_trait, Disease_class, PubmedID, Chr, Position, Ref, Alt, P_value, Variant_type
GAD_GRCh38_man <- data(GAD_GRCh38, envir = environment())

#' GRASP_GRCh38 Dataset
#'
#' This dataset contains GWAS information from the GRASP database.
#'
#' @name GRASP_GRCh38
#' @format A data frame with columns Rsid, Chr_37, Position_37, P_value, Pubmedid, Gwas_trait, Phenotype_escription, Chr, Position, Gene, Ref, Alt, Variant_type
GRASP_GRCh38_man <- data(GRASP_GRCh38, envir = environment())

#' GWAS_catalog_GRCh38 Dataset
#'
#' This dataset contains GWAS information from the GWAS catalog database.
#'
#' @name GWAS_catalog_GRCh38
#' @format A data frame with columns Rsid, Chr_37, Position_37, Gwas_trait, Gene, Variant_type, P_value, Pubmedid, Chr, Position, Ref, Alt
GWAS_catalog_GRCh38_man <- data(GWAS_catalog_GRCh38, envir = environment())

#' GWASdb_GRCh38 Dataset
#'
#' This dataset contains GWAS information from the GWASdb database.
#'
#' @name GWASdb_GRCh38
#' @format A data frame with columns Rsid, Chr_37, Position_37, Ref, Alt, P_value, Gwas_trait, Gene, Variant_type, Pubmedid, Chr, Position
GWASdb_GRCh38_man <- data(GWASdb_GRCh38, envir = environment())

#' JnO_GRCh38 Dataset
#'
#' This dataset contains GWAS information from the Johnson and O'donnell database.
#'
#' @name JnO_GRCh38
#' @format A data frame with columns Rsid, Gwas_trait, P_value Gene, Validation, Pubmedid, Chr, Position, Ref, Alt, Variant_type
JnO_GRCh38_man <- data(JnO_GRCh38, envir = environment())
