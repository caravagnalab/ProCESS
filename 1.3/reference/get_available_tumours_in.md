# Getting the tumour types available in a setup

This method returns the tumour types available for a set-up code.

## Usage

``` r
get_available_tumours_in(setup_code)
```

## Arguments

- setup_code:

  The set-up code whose available tumour types are requested.

## Value

A data frame reporting the types available for a set-up code.

## See also

[`MutationEngine`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine_class.md)
[`MutationEngine()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine.md)

## Examples

``` r
# get the types available for the "demo" set-up code
get_available_tumours_in("demo")
#> Downloading reference genome...
#> Reference genome downloaded
#> Decompressing reference genome...done
#> Downloading signature files...
#> Signature file downloaded
#> Downloading driver mutation file...
#> Driver mutation file downloaded
#> Decompressing driver mutation file...done
#> Downloading passenger CNAs file...
#> Passenger CNAs file downloaded
#> Decompressing passenger CNAs file...done
#> Downloading germline...
#> Germline downloaded
#> Decompressing mutations...
#> done
#>           type
#> 1          ACC
#> 2          AML
#> 3         ANGS
#> 4         ANSC
#> 5          BCC
#> 6         BLCA
#> 7         BRCA
#> 8        CCRCC
#> 9         CESC
#> 10        CHOL
#> 11      CLLSLL
#> 12         CML
#> 13        COAD
#> 14    COADREAD
#> 15        CSCC
#> 16          ES
#> 17        ESCA
#> 18         GBM
#> 19        GIST
#> 20         HCC
#> 21      HGGNOS
#> 22        HNSC
#> 23        KICH
#> 24      LGGNOS
#> 25        LICA
#> 26        LINC
#> 27        LIRI
#> 28         LMS
#> 29         LNM
#> 30        LUAD
#> 31        LUSC
#> 32         MBL
#> 33         MEL
#> 34        MLYM
#> 35         NBL
#> 36      NETNOS
#> 37         NHL
#> 38         NPC
#> 39       NSCLC
#> 40          OS
#> 41         OVT
#> 42        PAAD
#> 43       PANET
#> 44         PCM
#> 45        PGNG
#> 46      PLMESO
#> 47        PRAD
#> 48        PRCC
#> 49         RCC
#> 50        READ
#> 51         RMS
#> 52        SACA
#> 53     SARCNOS
#> 54        SCLC
#> 55        SKCM
#> 56 SOFT_TISSUE
#> 57        STAD
#> 58        THCA
#> 59        THYM
#> 60        UCEC
#> 61         UCS
#> 62       VULVA
#> 63        WDTC
```
