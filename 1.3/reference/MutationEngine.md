# Creating a mutation engine

This function downloads and sets up the data requires by a mutation
engine. Finally, it builds an object of the class
[`MutationEngine`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine_class.md)
itself.

## Usage

``` r
MutationEngine(directory, reference_src,
               SBS_signatures_src, indel_signatures_src,
               drivers_src, passenger_CNAs_src,
               germline_src)

MutationEngine(setup_code="demo")
```

## Arguments

- setup_code:

  The set-up code (optional).

- directory:

  The set-up directory (mandatory when `setup_code` is *not* provided).

- reference_src:

  The reference genome path or URL (mandatory when `setup_code` is *not*
  provided).

- SBS_signatures_src:

  The SBS signature file path or URL (mandatory when `setup_code` is
  *not* provided).

- indel_signatures_src:

  The indel signature file path or URL (mandatory when `setup_code` is
  *not* provided).

- drivers_src:

  The driver mutation file path or URL (mandatory when `setup_code` is
  *not* provided).

- passenger_CNAs_src:

  The passenger CNAs file path or URL (mandatory when `setup_code` is
  *not* provided).

- germline_src:

  The germline directory path or URL (mandatory when `setup_code` is
  *not* provided).

- germline_subject:

  The germline subject (optional).

- context_sampling:

  The number of reference contexts per context in the index (optional:
  default value is 30).

- COSMIC_account:

  A named list containing "email" and "password" of a valid COSMIC
  account (required to download mutational signatures from COSMIC site).

- max_index_size:

  The maximum size of an admitted indel and, as a consequence, the
  maximum size of a motif stored in the repeated sequence index
  (optional: default value is 50).

- max_repetition_storage:

  The maximum number of repetitions per type stored in the repeated
  sequence index (optional: default value is 500000).

- driver_CNA_min_distance:

  The minimum distance between a driver mutation and a passenger CNA.

- tumour_type:

  The type of tumour. This is currently used to select the admissible
  passenger CNAs. If any passenger CNA in the dataset is admissible, use
  the the empty string `""` (optional: default value is `""`).

- avoid_homozygous_losses:

  An optional Boolean flag to avoid homozygous losses. When set to
  `TRUE`, passenger CNAs will be exclusively applied to regions covered
  by two alleles at least. (default: TRUE).

- quiet:

  An optional Boolean flag to avoid the progress bar (default: FALSE).

## Value

A mutation engine object.

## Details

There are two building modalities: the first one is more general, but it
requires to specify all the data sources; the second one adopts some
pre-set configurations, but it is sufficient in many cases.

The first building modality requires to specify the directory in which
the data must be saved, the path or URL of the reference sequence, the
mutational signatures, the driver SNVs file, the passenger CNAs file,
and the germline data directory through the `directory`,
`reference_src`, `SBS_src`, `drivers_src`, `passenger_CNAs_src`, and
`germline_src`, respectively.

The second building modality exclusively requires a set-up code
(parameter `setup_code`). The list of supported set-up codes can be
obtained by using the function
[`get_mutation_engine_codes()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_mutation_engine_codes.md).

Whenever the mutational signatures are meant to be downloaded from the
COSMIC site, a valid COSMIC account is needed and can be provided by the
parameter `COSMIC_account`.

The number of context sampling is an optional parameter that allows
sampling the reference contexts while building the context index. This
parameter, which is set to 100 by default, specifies how many
occurrences of the same context must be identified before adding one of
them to the context index. The larger the number of context sampling,
the larger the context index. On the other side, the lower the number of
context sampling, the lower the number of sites in the reference genome
that can be affected by simulated mutations.

If the parameters of a mutation engine construction match those of a
previous construction, then the corresponding reference sequence, the
SBS file, and the previously built context index are loaded from the
set-up directory avoiding further computations.

## See also

[`get_mutation_engine_codes()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_mutation_engine_codes.md)
provides a list of the supported set-up codes;
[`MutationEngine$get_germline_subjects()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-get_germline_subjects.md)
to get the available germline subjects;
[`MutationEngine$set_germline_subject()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-set_germline_subject.md)
to set the active germline subject;
[`MutationEngine$get_active_germline()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-get_active_germline.md)
to get the active germline subject.

[`MutationEngine`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine_class.md)

## Examples

``` r
# set the reference and SBS URLs
reference_url <- paste0("https://ftp.ensembl.org/pub/grch37/release-111/",
  "fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.",
  "dna.chromosome.22.fa.gz")
sbs_url <- paste0("https://zenodo.org/records/15875185/files/",
  "SBS_demo_signatures.txt")
indel_url <- paste0("https://zenodo.org/records/15875185/files/",
  "indel_demo_signatures.txt")
drivers_url <- paste0("https://zenodo.org/records/15875185/files/",
  "driver_mutations_hg19.csv.bz2")
passenger_CNAs_url <- paste0("https://zenodo.org/records/15875185/",
  "files/passenger_CNAs_hg19.csv.bz2")
germline_url <- paste0("https://zenodo.org/records/15875185/files/",
  "germline_data_demo.tar.bz2")

# build a mutation engine and save the required files into the
# directory "Test". The `drivers_url` parameter is optional, but
# it is suggested to avoid passenger mutations on driver loci.
m_engine <- MutationEngine(directory = "Test",
  reference_src = reference_url,
  SBS_signatures_src = sbs_url,
  indel_signatures_src = indel_url,
  drivers_src = drivers_url,
  passenger_CNAs_src = passenger_CNAs_url,
  germline_src = germline_url)
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
#> Building context index...
#> 
 [█---------------------------------------] 0% [00m:00s] Processing chr. 22                                                                   

 [█████████████████-----------------------] 40% [00m:00s] Processing chr. 22                                                                  

 [█████████████████████████████████-------] 81% [00m:02s] Processing chr. 22                                                                  

 [████████████████████████████████████████] 100% [00m:02s] Context index built                                                                

#> 
 [█---------------------------------------] 0% [00m:00s] Saving context index                                                                 

 [████████████████████████████████████████] 100% [00m:00s] Context index saved                                                                

#> done
#> Building repeated sequence index...
#> 
 [█---------------------------------------] 0% [00m:00s] Reading 22                                                                           

 [████████████████████████████████████████] 100% [00m:00s] Reading 22                                                                         

#> 
 [████████████████████████████████████████] 100% [00m:00s] Reading 22                                                                         

#> 
 [████████████████████████████████████████] 100% [00m:00s] Reading 22                                                                         

#> 
 [████████████████████████████████████████] 100% [00m:00s] Reading 22                                                                         

#> 
 [████████████████████████████████████████] 100% [00m:00s] Reading 22                                                                         

#> 
 [████████████████████████████████████████] 100% [00m:00s] Reading 22                                                                         

#> 
 [████████████████████████████████████████] 100% [00m:00s] Reading 22                                                                         

#> 
 [████████████████████████████████████████] 100% [00m:00s] RS index built                                                                     

#> 
 [█---------------------------------------] 0% [00m:00s] Saving RS index                                                                      

 [█---------------------------------------] 0% [00m:00s] Saving RS index                                                                      

 [████████████████████--------------------] 48% [00m:02s] Saving RS index                                                                     
done
#> 
 [████████████████████████████████████████] 100% [00m:02s] RS index saved                                                                     

#> 
 [█---------------------------------------] 0% [00m:00s] Loading germline                                                                     

 [████████████████████████████████████████] 100% [00m:00s] Germline loaded                                                                    

#> 
 [█---------------------------------------] 0% [00m:00s] Saving germline                                                                      

 [████████████████████████████████████████] 100% [00m:00s] Germline saved                                                                     


# if the parameters of a mutation engine construction match those of a
# previous construction, then the corresponding reference sequence,
# the SBS file, and the previously built context index are loaded from
# the set-up directory avoiding further computations.
m_engine <- MutationEngine("Test", reference_url, sbs_url, indel_url,
  drivers_url, passenger_CNAs_url, germline_url)
#> 
 [█---------------------------------------] 0% [00m:00s] Loading context index                                                                

 [████████████████████████████████████████] 100% [00m:00s] Context index loaded                                                               

#> 
 [█---------------------------------------] 0% [00m:00s] Loading RS index                                                                     

 [█████████████---------------------------] 31% [00m:01s] Loading RS index                                                                    

 [██████████████████████████--------------] 63% [00m:02s] Loading RS index                                                                    

 [██████████████████████████████████████--] 93% [00m:03s] Loading RS index                                                                    

 [████████████████████████████████████████] 100% [00m:03s] RS index loaded                                                                    

#> 
 [█---------------------------------------] 0% [00m:00s] Loading germline                                                                     

 [████████████████████████████████████████] 100% [00m:00s] Germline loaded                                                                    


# if the `context_sampling` parameter changes, a new context index is
# built, while neither the reference sequence nor the SBS file are
# downloaded again.
m_engine <- MutationEngine("Test", reference_url, sbs_url, indel_url,
  drivers_url, passenger_CNAs_url, germline_url,
  context_sampling = 50)
#> Building context index...
#> 
 [█---------------------------------------] 0% [00m:00s] Processing chr. 22                                                                   

 [████████████████████--------------------] 49% [00m:01s] Processing chr. 22                                                                  

 [████████████████████████████████████----] 89% [00m:02s] Processing chr. 22                                                                  

 [████████████████████████████████████████] 100% [00m:02s] Context index built                                                                

#> 
 [█---------------------------------------] 0% [00m:00s] Saving context index                                                                 

 [████████████████████████████████████████] 100% [00m:00s] Context index saved                                                                

#> done
#> 
 [█---------------------------------------] 0% [00m:00s] Loading RS index                                                                     

 [█████████████---------------------------] 32% [00m:01s] Loading RS index                                                                    

 [██████████████████████████--------------] 64% [00m:02s] Loading RS index                                                                    

 [██████████████████████████████████████--] 94% [00m:03s] Loading RS index                                                                    

 [████████████████████████████████████████] 100% [00m:03s] RS index loaded                                                                    

#> 
 [█---------------------------------------] 0% [00m:00s] Loading germline                                                                     

 [████████████████████████████████████████] 100% [00m:00s] Germline loaded                                                                    


# a further construction with the same parameters avoids both
# downloads and context index construction.
m_engine <- MutationEngine("Test", reference_url, sbs_url, indel_url,
  drivers_url, passenger_CNAs_url, germline_url,
  context_sampling = 50)
#> 
 [█---------------------------------------] 0% [00m:00s] Loading context index                                                                

 [████████████████████████████████████████] 100% [00m:00s] Context index loaded                                                               

#> 
 [█---------------------------------------] 0% [00m:00s] Loading RS index                                                                     

 [█████████████---------------------------] 31% [00m:01s] Loading RS index                                                                    

 [█████████████████████████---------------] 62% [00m:02s] Loading RS index                                                                    

 [█████████████████████████████████████---] 90% [00m:03s] Loading RS index                                                                    

 [████████████████████████████████████████] 100% [00m:03s] RS index loaded                                                                    

#> 
 [█---------------------------------------] 0% [00m:00s] Loading germline                                                                     

 [████████████████████████████████████████] 100% [00m:00s] Germline loaded                                                                    


m_engine
#> MutationEngine
#>  Passenger rates
#> 
#>  Driver mutations
#> 
#>  Timed Exposure
#>    SBS Timed Exposures
#> 
#>    indel Timed Exposures
#> 

# the parameters `directory`, `reference_src`, `SBS_src`, `drivers_src`,
# `passenger_CNAs_src`, and `germline_src` can be avoided by providing
# the `setup_code` parameter. The set-up code `demo` is provided among
# those available for testing purpose.
m_engine <- MutationEngine(setup_code = "demo")
#> 
 [█---------------------------------------] 0% [00m:00s] Loading context index                                                                

 [████████████████████████████████████████] 100% [00m:00s] Context index loaded                                                               

#> 
 [█---------------------------------------] 0% [00m:00s] Loading RS index                                                                     

 [█████████████---------------------------] 31% [00m:01s] Loading RS index                                                                    

 [██████████████████████████--------------] 63% [00m:02s] Loading RS index                                                                    

 [█████████████████████████████████████---] 92% [00m:03s] Loading RS index                                                                    

 [████████████████████████████████████████] 100% [00m:03s] RS index loaded                                                                    

#> 
 [█---------------------------------------] 0% [00m:00s] Loading germline                                                                     

 [████████████████████████████████████████] 100% [00m:00s] Germline loaded                                                                    


# the `context_sampling` can be used also when a pre-defined set-up
# configuration is adopted.
m_engine <- MutationEngine(setup_code = "demo", context_sampling = 50)
#> Building context index...
#> 
 [█---------------------------------------] 0% [00m:00s] Processing chr. 22                                                                   

 [████████████████████--------------------] 49% [00m:01s] Processing chr. 22                                                                  

 [████████████████████████████████████----] 89% [00m:02s] Processing chr. 22                                                                  

 [████████████████████████████████████████] 100% [00m:02s] Context index built                                                                

#> 
 [█---------------------------------------] 0% [00m:00s] Saving context index                                                                 

 [████████████████████████████████████████] 100% [00m:00s] Context index saved                                                                

#> done
#> 
 [█---------------------------------------] 0% [00m:00s] Loading RS index                                                                     

 [█████████████---------------------------] 30% [00m:01s] Loading RS index                                                                    

 [█████████████████████████---------------] 60% [00m:02s] Loading RS index                                                                    

 [████████████████████████████████████----] 89% [00m:03s] Loading RS index                                                                    

 [████████████████████████████████████████] 100% [00m:03s] RS index loaded                                                                    

#> 
 [█---------------------------------------] 0% [00m:00s] Loading germline                                                                     

 [████████████████████████████████████████] 100% [00m:00s] Germline loaded                                                                    


m_engine
#> MutationEngine
#>  Passenger rates
#> 
#>  Driver mutations
#> 
#>  Timed Exposure
#>    SBS Timed Exposures
#> 
#>    indel Timed Exposures
#> 

# remove the "demo" directory
unlink("demo", recursive = TRUE)

# Some of the pre-defined configurations requires to download the mutational
# signatures from the COSMIC site which requires an account (e.g., "GRCh37"
# and "GRCh38"). The COSMIC account can be passed to `MutationEngine()` as
# follows
m_engine <- MutationEngine(setup_code = "demo",
  COSMIC_account = list(email = "foo@bar.org",
  password = "********"))
#> Downloading reference genome...
#> Reference genome downloaded
#> Decompressing reference genome...done
#> Downloading signature files...
#> Signature file downloaded
#> Downloading driver mutation file...
#> Driver mutation file downloaded
#> Decompressing driver mutation file...
#> Error in getWritablePathname.Arguments(static, ...): File already exists: /Users/alberto/Library/CloudStorage/Dropbox/Lavoro/Code/ProCESS/ProCESS/docs/reference/demo/drivers.txt
m_engine
#> MutationEngine
#>  Passenger rates
#> 
#>  Driver mutations
#> 
#>  Timed Exposure
#>    SBS Timed Exposures
#> 
#>    indel Timed Exposures
#> 

# remove the "demo" directory
unlink("demo", recursive = TRUE)

# In alternative, pre-download the mutational signatures and pass their
# paths to `MutationEngine()` as parameters.
m_engine <- MutationEngine(setup_code = "demo",
  SBS_signatures_src = "Test/SBS_signatures.txt",
  indel_signatures_src = "Test/indel_signatures.txt")
#> Downloading reference genome...
#> Reference genome downloaded
#> Decompressing reference genome...done
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
#> Building context index...
#> 
 [█---------------------------------------] 0% [00m:00s] Processing chr. 22                                                                   

 [█████████████████-----------------------] 40% [00m:00s] Processing chr. 22                                                                  

 [█████████████████████████████████-------] 81% [00m:02s] Processing chr. 22                                                                  

 [████████████████████████████████████████] 100% [00m:02s] Context index built                                                                

#> 
 [█---------------------------------------] 0% [00m:00s] Saving context index                                                                 

 [█---------------------------------------] 0% [00m:00s] Saving context index                                                                 

 [████████████████████████████████████████] 100% [00m:00s] Context index saved                                                                

#> done
#> 
 [█---------------------------------------] 0% [00m:00s] Loading RS index                                                                     

 [█---------------------------------------] 1% [00m:06s] Loading RS index                                                                     

 [███████████-----------------------------] 27% [00m:07s] Loading RS index                                                                    

 [███████████████████████-----------------] 55% [00m:08s] Loading RS index                                                                    

 [█████████████████████████████████-------] 80% [00m:09s] Loading RS index                                                                    

 [████████████████████████████████████████] 100% [00m:10s] RS index loaded                                                                    

#> 
 [█---------------------------------------] 0% [00m:00s] Loading germline                                                                     

 [█---------------------------------------] 1% [00m:01s] Loading germline                                                                     

 [████████████████████████████████████████] 100% [00m:01s] Germline loaded                                                                    

m_engine
#> MutationEngine
#>  Passenger rates
#> 
#>  Driver mutations
#> 
#>  Timed Exposure
#>    SBS Timed Exposures
#> 
#>    indel Timed Exposures
#> 

# remove the "Test" and "demo" directories
unlink("Test", recursive = TRUE)
unlink("demo", recursive = TRUE)
```
