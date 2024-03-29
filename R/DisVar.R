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

DisVar <- function(file = "file_name.vcf"){
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

  #Read files
  cat("\nReading files...\n")
  GWASdb_GRCh38 <- V1 <- P_value <- Rsid <- Chr <- V4 <- V5 <- Ref <- Alt <- Gwas_trait <- Gene <- Variant_type <- V6 <- V7 <- V8 <- V2 <- GRASP_GRCh38 <-GWAS_catalog_GRCh38 <- JnO_GRCh38 <- Disease <- Confident <- NULL
  GWASdb_GRCh38 <- DisVar::GWASdb_GRCh38
  GRASP_GRCh38 <- DisVar::GRASP_GRCh38
  GWAS_catalog_GRCh38 <- DisVar::GWAS_catalog_GRCh38
  GAD_GRCh38 <- DisVar::GAD_GRCh38
  JnO_GRCh38 <- DisVar::JnO_GRCh38
  grep <- "grep -v '^#'"
  suppressWarnings(variant_data <- fread(cmd = paste(grep, file), sep = "\t", stringsAsFactors=FALSE, showProgress=TRUE, header =FALSE))

  . <- NULL
  chr <- c()
  allele_db_list <- c()
  allele_sample_list <- c()
  variant_data$V1<-gsub("chr","",as.character(variant_data$V1))
  variant_data$V2 <- as.integer(variant_data$V2)
  cat("Reading files...DONE\n")

  #Search databases
  cat("Searching...\n")


  #find variant that in the GWASdb database

  op_df_GWASdb <- setDT(GWASdb_GRCh38)[setDT(variant_data),
                                       on = .(Position = V2, Chr = V1)][P_value < 0.0000001, .(Rsid, Chr, Position, V4, V5, Ref, Alt, P_value, Gwas_trait, Gene, Variant_type, V6, V7, V8)]

  if (nrow(op_df_GWASdb) > 0)
  {
    op_df_GWASdb$DB <- "GWASdb"
  }

  #find variant that in the GRASP database
  op_df_GRASP <- setDT(GRASP_GRCh38)[setDT(variant_data),
                                     on = .(Position = V2, Chr = V1)][P_value < 0.0000001, .(Rsid, Chr, Position, V4, V5, Ref, Alt, P_value, Gwas_trait, Gene, Variant_type, V6, V7, V8)]

  if (nrow(op_df_GRASP) > 0)
  {
    op_df_GRASP$DB <- "GRASP"
  }

  #find variant that in the GWAS catalog database
  op_df_GWAS_catalog <- setDT(GWAS_catalog_GRCh38)[setDT(variant_data),
                                                   on = .(Position = V2, Chr = V1)][P_value < 0.0000001, .(Rsid, Chr, Position, V4, V5, Ref, Alt, P_value, Gwas_trait, Gene, Variant_type, V6, V7, V8)]

  if (nrow(op_df_GWAS_catalog) > 0)
  {
    op_df_GWAS_catalog$DB <- "GWAS Catalog"
  }

  #find variant that in the GAD database
  GAD_GRCh38$P_value <- 0
  op_df_GAD <- setDT(GAD_GRCh38)[setDT(variant_data),
                                 on = .(Position = V2, Chr = V1)][P_value < 0.0000001, .(Rsid , Chr , Position, V4, V5, Ref, Alt, P_value, Gwas_trait , Gene , Variant_type, V6, V7, V8)]
  op_df_GAD$P_value <- NA

  if (nrow(op_df_GAD) > 0)
  {
    op_df_GAD$DB <- "GADCDC"
  }

  #find variant that in the Johnson and O'Donnell database
  op_df_JnO <- setDT(JnO_GRCh38)[setDT(variant_data),
                                 on = .(Position = V2, Chr = V1)][P_value < 0.0000001,
                                                                  .(Rsid, Chr, Position, V4, V5, Ref, Alt, P_value, Gwas_trait, Gene, Variant_type, V6, V7, V8)]

  if (nrow(op_df_JnO) > 0)
  {
    op_df_JnO$DB <- "Johnson and O'Donnell"
  }
  op_df <- rbind(op_df_GWASdb,op_df_GRASP,op_df_GWAS_catalog,op_df_GAD,op_df_JnO, fill=TRUE)

  cat("Searching...DONE\n")

  cat("Processing results...\n")
  #If no variant is found, Exit program
  if (nrow(op_df) == 0)
  {
    stop(paste('No variant was found in the database'))
    #EXIT command
  }

  #create results in a table
  align_df <- setNames(data.frame(matrix(ncol = 13, nrow = nrow(op_df))), c("Disease", "Chrom", "Position", "Gene", "Variant ID", "Variant Type", "Allele Sample", "Allele DB", "Confident", "DB", "Qual", "Filter", "Info"))

  align_df["Disease"] <- sqldf('SELECT Gwas_trait FROM op_df')
  align_df["Chrom"] <- sqldf('SELECT Chr FROM op_df')
  align_df["Position"] <- sqldf('SELECT Position FROM op_df')
  align_df["Gene"] <- sqldf('SELECT Gene FROM op_df')
  align_df["Variant ID"] <- sqldf('SELECT Rsid FROM op_df')
  align_df["Variant Type"] <- sqldf('SELECT Variant_type FROM op_df')

  op_df$Ref <- as.character(op_df$Ref)
  op_df$Alt <- as.character(op_df$Alt)
  for (i in 1:nrow(op_df))
  {
    ref_db <- op_df[i,6]
    alt_db <- op_df[i,7]
    if (is.na(ref_db)) {allele_db_list <- append(allele_db_list,NA)}
    else
    {
      allele_db = paste(ref_db, alt_db, sep = ">")
      allele_db_list <- append(allele_db_list,allele_db)
    }
  }
  align_df["Allele DB"] <- allele_db_list

  for (i in 1:nrow(op_df))
  {
    ref_sample <- op_df[i,4]
    alt_sample <- op_df[i,5]
    allele_sample = paste(ref_sample, alt_sample, sep = ">")
    allele_sample_list <- append(allele_sample_list,allele_sample)
  }
  align_df["Allele Sample"] <- allele_sample_list

  align_df["Confident"] <- sqldf('SELECT P_value FROM op_df')
  align_df["DB"] <- sqldf('SELECT DB FROM op_df')
  align_df["Qual"] <- sqldf('SELECT V6 FROM op_df')
  align_df["Filter"] <- sqldf('SELECT V7 FROM op_df')
  align_df["Info"] <- sqldf('SELECT V8 FROM op_df')

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
