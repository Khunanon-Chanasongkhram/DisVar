% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DisVar.R
\name{DisVar}
\alias{DisVar}
\title{Find diseases from the VCF file}
\usage{
DisVar(
  file,
  GWASdb = TRUE,
  GRASP = TRUE,
  GWASCat = TRUE,
  GAD = TRUE,
  JohnO = TRUE,
  ClinVar = TRUE,
  p_value = 5e-08,
  merge_result = FALSE,
  runOnShiny = FALSE,
  output_path = NULL
)
}
\arguments{
\item{file}{A character string specifying the path to the input VCF file (Variant Call Format) for analysis.}

\item{GWASdb}{Logical. If \code{TRUE}, use the GWAS database for variant analysis. Defaults to \code{TRUE}.}

\item{GRASP}{Logical. If \code{TRUE}, use the GRASP database for variant analysis. Defaults to \code{TRUE}.}

\item{GWASCat}{Logical. If \code{TRUE}, use the GWAS Catalog database for variant analysis. Defaults to \code{TRUE}.}

\item{GAD}{Logical. If \code{TRUE}, use the Genetic Association Database (GAD) for variant analysis. Defaults to \code{TRUE}.}

\item{JohnO}{Logical. If \code{TRUE}, use the JohnO database for variant analysis. Defaults to \code{TRUE}.}

\item{ClinVar}{Logical. If \code{TRUE}, use the ClinVar database for variant analysis. Defaults to \code{TRUE}.}

\item{p_value}{Numeric. The p-value threshold for filtering significant associations. Variants with p-values greater than this threshold are excluded. Defaults to \code{1e-7}.}

\item{runOnShiny}{Logical. If \code{TRUE}, the function is optimized for running within a Shiny application. Defaults to \code{FALSE}.}
}
\value{
A data frame containing variant information from various databases
}
\description{
This function finds diseases from a VCF file by comparing the variants in the VCF file with those in various databases.
}
\note{
The GRASP database is now loaded from four parts (\code{GRASP_GRCh38_part1}, \code{GRASP_GRCh38_part2}, \code{GRASP_GRCh38_part3}, and \code{GRASP_GRCh38_part4}) and merged internally.
}
