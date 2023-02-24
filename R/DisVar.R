#' Find diseases from the VCF file
#'
#' This function finds diseases from a VCF file by comparing the variants in the VCF file with those in various databases.
#' @param file the name of the VCF file.
#' @param skip the number of lines to skip at the beginning of the VCF file. (default: 0)
#' @return A data frame containing variant information from various databases
#' @export
#'
#' @examples DisVar("file_name.vcf")
#'
DisVar <- function(file = "vcf_file_name.vcf"){
  if (!require("sqldf")) install.packages("sqldf")
  if (!require("data.table")) install.packages("data.table")

  # Validate inputs
  if (!file.exists(file)) {
    stop("File does not exist")
  }
  if (!grepl(".vcf", file, ignore.case = TRUE)) {
    stop("Input file must be a VCF file")
  }

  #Read files
  print("Reading files...")
  data(GWASdb_GRCh38)
  data(GRASP_GRCh38)
  data(GWAS_catalog_GRCh38)
  data(GAD_GRCh38)
  data(JnO_GRCh38)
  grep <- "grep -v '^#'"
  #variant_data <- fread(file = paste0(file), select = 1:5, header=TRUE, stringsAsFactors=TRUE, skip=paste0(skip))
  #variant_data <- fread(file = paste0(file), header=TRUE, stringsAsFactors=TRUE, skip=paste0(skip))
  variant_data <- fread(paste(grep, file), select = 1:5, stringsAsFactors=TRUE)

  col_list <- paste0("V", 1:5)
  colnames(variant_data) <- col_list

  chr <- c()
  allele_db_list <- c()
  allele_sample_list <- c()
  variant_data$V1<-gsub("chr","",as.character(variant_data$V1))
  print("Reading files...DONE")


  #Search databases
  print("Searching...")

  #find variant that in the GWASdb database
  op_df_GWASdb <- sqldf("
                        SELECT
                        GWASdb_GRCh38.Rsid,
                        GWASdb_GRCh38.Chr,
                        GWASdb_GRCh38.Position,
                        variant_data.V4,
                        variant_data.V5,
                        GWASdb_GRCh38.Ref,
                        GWASdb_GRCh38.Alt,
                        GWASdb_GRCh38.P_value,
                        GWASdb_GRCh38.Gwas_trait
                        FROM GWASdb_GRCh38, variant_data
                        WHERE variant_data.V2 = GWASdb_GRCh38.Position AND variant_data.V1 = GWASdb_GRCh38.Chr AND GWASdb_GRCh38.P_value < 0.0000001
                        ")

  if (nrow(op_df_GWASdb) > 0)
  {
    op_df_GWASdb$DB <- "GWASdb"
  }

  #find variant that in the GRASP database
  GRASP_GRCh38$Ref = NA
  GRASP_GRCh38$Alt = NA
  op_df_GRASP <- sqldf("
                        SELECT
                        GRASP_GRCh38.Rsid,
                        GRASP_GRCh38.Chr,
                        GRASP_GRCh38.Position,
                        variant_data.V4,
                        variant_data.V5,
                        GRASP_GRCh38.Ref,
                        GRASP_GRCh38.Alt,
                        GRASP_GRCh38.P_value,
                        GRASP_GRCh38.Phenotype
                        FROM GRASP_GRCh38, variant_data
                        WHERE variant_data.V2 = GRASP_GRCh38.Position AND variant_data.V1 = GRASP_GRCh38.Chr AND GRASP_GRCh38.P_value < 0.0000001
                       ")
  colnames(op_df_GRASP)[9] <- "Gwas_trait"
  if (nrow(op_df_GRASP) > 0)
  {
    op_df_GRASP$DB <- "GRASP"
  }

  #find variant that in the GWAS catalog database
  GWAS_catalog_GRCh38$Ref = NA
  GWAS_catalog_GRCh38$Alt = NA
  op_df_GWAS_catalog <- sqldf("
                        SELECT
                        GWAS_catalog_GRCh38.Rsid,
                        GWAS_catalog_GRCh38.Chr,
                        GWAS_catalog_GRCh38.Position,
                        variant_data.V4,
                        variant_data.V5,
                        GWAS_catalog_GRCh38.Ref,
                        GWAS_catalog_GRCh38.Alt,
                        GWAS_catalog_GRCh38.P_value,
                        GWAS_catalog_GRCh38.Gwas_trait
                        FROM GWAS_catalog_GRCh38, variant_data
                        WHERE variant_data.V2 = GWAS_catalog_GRCh38.Position AND variant_data.V1 = GWAS_catalog_GRCh38.Chr AND GWAS_catalog_GRCh38.P_value < 0.0000001
                        ")

  if (nrow(op_df_GWAS_catalog) > 0)
  {
    op_df_GWAS_catalog$DB <- "GWAS Catalog"
  }

  #find variant that in the GAD database
  GAD_GRCh38$Ref = NA
  GAD_GRCh38$Alt = NA
  GAD_GRCh38$P_value = NA
  op_df_GAD <- sqldf("
                      SELECT
                      GAD_GRCh38.Rsid,
                      GAD_GRCh38.Chr,
                      GAD_GRCh38.Position,
                      variant_data.V4,
                      variant_data.V5,
                      GAD_GRCh38.Ref,
                      GAD_GRCh38.Alt,
                      GAD_GRCh38.P_value,
                      GAD_GRCh38.Gwas_trait
                      FROM GAD_GRCh38, variant_data
                      WHERE variant_data.V2 = GAD_GRCh38.Position AND variant_data.V1 = GAD_GRCh38.Chr")

  if (nrow(op_df_GAD) > 0)
  {
    op_df_GAD$DB <- "GADCDC"
  }

  #find variant that in the Johnson and O'Donnell database
  JnO_GRCh38$Ref = NA
  JnO_GRCh38$Alt = NA
  op_df_JnO <- sqldf("
                      SELECT
                      JnO_GRCh38.Rsid,
                      JnO_GRCh38.Chr,
                      JnO_GRCh38.Position,
                      variant_data.V4,
                      variant_data.V5,
                      JnO_GRCh38.Ref,
                      JnO_GRCh38.Alt,
                      JnO_GRCh38.P_value,
                      JnO_GRCh38.Gwas_trait
                      FROM JnO_GRCh38, variant_data
                      WHERE variant_data.V2 = JnO_GRCh38.Position AND variant_data.V1 = JnO_GRCh38.Chr AND JnO_GRCh38.P_value < 0.0000001
                        ")

  if (nrow(op_df_JnO) > 0)
  {
    op_df_JnO$DB <- "Johnson and O'Donnell"
  }
  op_df <- rbind(op_df_GWASdb,op_df_GRASP,op_df_GWAS_catalog,op_df_GAD,op_df_JnO)

  print("Searching...DONE")

  print("Processing results...")
  #If no variant is found, Exit programe
  if (nrow(op_df) == 0)
  {
    stop(paste('No variant was found in the database'))
    #EXIT command
  }

  #create results in a table
  align_df <- setNames(data.frame(matrix(ncol = 8, nrow = nrow(op_df))), c("Disease", "Chr", "Position", "Variant_id", "Allele sample", "Allele DB", "Confident", "DB"))


  align_df["Disease"] <- sqldf('SELECT Gwas_trait FROM op_df')
  align_df["Chr"] <- sqldf('SELECT Chr FROM op_df')
  align_df["Position"] <- sqldf('SELECT Position FROM op_df')
  align_df["Variant_id"] <- sqldf('SELECT Rsid FROM op_df')

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
  align_df["Allele sample"] <- allele_sample_list

  align_df["Confident"] <- sqldf('SELECT P_value FROM op_df')
  align_df["DB"] <- sqldf('SELECT DB FROM op_df')

  #aligned_df <- align_df[order(align_df$Disease),]
  aligned_df <- align_df %>% arrange(Disease, Confident)
  aligned_df$Disease <- as.character(aligned_df$Disease)
  aligned_df$Disease[duplicated(aligned_df$Disease)] <- ''
  colnames(aligned_df)[7] <- "Confident/P-value"
  print("Processing results...DONE")
  print("Generating result file...")
  write.table(aligned_df,file= paste0(file,"_diseases_output.txt"), quote = FALSE, sep = '\t', row.names = FALSE)
  print("Generating result file...DONE")
}
