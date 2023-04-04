Khunanon Chanasongkhram
04/04/2023

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

## Installing the Package in R Studio

You can install the development version of DisVar from
[GitHub](https://github.com/) to install the package, you can use the
`devtools` package and run the following command:

``` r
# install.packages("devtools")
devtools::install_github("Khunanon-Chanasongkhram/DisVar")
```

## Installing the Package on Linux servers

### Step 1: Installing Minoconda 3

To install Miniconda via the terminal, you can follow these steps:
1.Create a directory to install Miniconda into:

    mkdir -p ~/miniconda3

2.Download the latest Python 3 based install script for Linux 64 bit:

    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh

3.Run the install script:

    bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3

After then, you can restart your shell, and conda will be ready for use.
For further details regarding conda installation please see [Miniconda
linux
installers](https://docs.conda.io/en/latest/miniconda.html#linux-installers)

### Step 2: Installing R using conda

Create a new conda environment for R by running the following command:

    conda create -n r-environment r-base

Please answer “y” for the question “Proceed (\[y\]/n)?”. This will
create a new environment named ‘r-environment’ with R installed. Now you
can activate the new environment by running:

    conda activate r-environment

and deactivate an active environment by running:

    conda deactivate

### Step 3: Installing the DisVar Package

After activate r-environment via conda. Run R in Linux by following
command:

    R

Then install the package

    remotes::install_github(repo = "Khunanon-Chanasongkhram/DisVar",dependencies = TRUE,upgrade = TRUE)

Now the DisVar is ready to use!

## Dependencies

The package depends on the following packages:

- `sqldf`
- `data.table`

## Usage

The main function in the package is `DisVar()`, which takes one
arguments:

- `file`: the name of the VCF file.

The function returns a data frame containing variant information from
various databases.

### Examples

Here’s an example of how to use the function:

``` r
library(DisVar)
DisVar("file_name.vcf")
#> Loading required package: sqldf
#> Warning: package 'sqldf' was built under R version 4.2.2
#> Loading required package: gsubfn
#> Loading required package: proto
#> Loading required package: RSQLite
#> Loading required package: data.table
#> Warning: package 'data.table' was built under R version 4.2.2
#> data.table 1.14.8 using 6 threads (see ?getDTthreads).  Latest news: r-datatable.com
#> Loading required package: dplyr
#> Warning: package 'dplyr' was built under R version 4.2.2
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:data.table':
#> 
#>     between, first, last
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
#> 
#> Reading files...
#> Reading files...DONE
#> Searching...
#> Searching...DONE
#> Processing results...
#> Processing results...DONE
#> Generating result file...
#> Generating result file...DONE
#> The output file is saved as: file_name_diseases_output.txt in the directory: C:/DisVar
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
