Khunanon Chanasongkhram
24/01/2023

<!-- README.md is generated from README.Rmd. Please edit that file -->

# DisVar: An R library for identifying variants associated with diseases using personal genetic information

<!-- badges: start -->
<!-- badges: end -->

DisVar is an R package that finds diseases from a VCF file by comparing
the variants in the VCF file with those in various databases. The
package currently supports five databases: GWASdb, GRASP, GWAS Catalog,
GAD, and Johnson and O’Donnell.

## Features

- Finds diseases from a VCF file by comparing the variants in the VCF
  file with those in various databases.

- Supports five databases: GWASdb, GRASP, GWAS Catalog, GAD, and Johnson
  and O’Donnell.

## Installation

You can install the development version of DisVar from
[GitHub](https://github.com/) to install the package, you can use the
`devtools` package and run the following command:

``` r
# install.packages("devtools")
devtools::install_github("Khunanon-Chanasongkhram/DisVar")
```

## Dependencies

The package depends on the following packages:

- `sqldf`
- `data.table`

## Usage

The main function in the package is `DisVar()`, which takes two
arguments:

- `file`: the name of the VCF file.

- `skip`: the number of lines to skip at the beginning of the VCF file.
  (default: 0)

The function returns a data frame containing variant information from
various databases.

### Examples

Here’s an example of how to use the function:

``` r
library(DisVar)
DisVar("file_name.vcf", skip = 28)
#> Loading required package: sqldf
#> Loading required package: gsubfn
#> Loading required package: proto
#> Loading required package: RSQLite
#> Loading required package: data.table
#> Warning: package 'data.table' was built under R version 4.2.2
#> [1] "Reading files..."
#> [1] "Reading files...DONE"
#> [1] "Searching..."
#> [1] "Searching...DONE"
#> [1] "Processing results..."
#> [1] "Processing results...DONE"
#> [1] "Generating result file..."
#> [1] "Generating result file...DONE"
```

You will get this output file: `file_name`\_diseases_output.txt

When you test the program with the test file (`file_name.vcf`), the
results are shown in the following table.

| Disease                                                  | Chr | Position | Variant_id | Allele sample | Allele DB | Confident/P-value | DB           |
|----------------------------------------------------------|-----|----------|------------|---------------|-----------|-------------------|--------------|
| Breast cancer                                            | 3   | 4700592  | rs6762644  | A\>G          | A\>G      | 2.00E-12          | GWASdb       |
|                                                          | 3   | 4700592  | rs6762644  | A\>G          | NA        | 2.20E-12          | GRASP        |
| Breast cancer (Estrogen receptor positive breast cancer) | 3   | 4700592  | rs6762644  | A\>G          | NA        | 1.40E-08          | GRASP        |
| Breast cancer (invasive breast cancer)                   | 3   | 4700592  | rs6762644  | A\>G          | NA        | 1.20E-09          | GRASP        |
|                                                          | 3   | 4700592  | rs6762644  | A\>G          | NA        | 2.00E-12          | GWAS Catalog |
|                                                          | 3   | 4700592  | rs6762644  | A\>G          | NA        | 4.00E-18          | GWAS Catalog |
| Breast_cancer                                            | 3   | 4700592  | rs6762644  | A\>G          | NA        | 9.00E-12          | GWAS Catalog |

## Additional Resources

- The VCF format is a standard file format for storing genetic variation
  data, and it is widely used in bioinformatics research. More
  information on the VCF format can be found at the [VCF specification
  page](http://samtools.github.io/hts-specs/VCFv4.3.pdf)

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/DisVar)](https://cran.r-project.org/package=DisVar)
