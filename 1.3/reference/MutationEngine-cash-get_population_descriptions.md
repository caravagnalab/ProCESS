# Getting the population descriptions

This method returns the population descriptions.

## Value

A data frame containing the population descriptions.

## Details

The population descriptions are stored in a data frame describing the
populations. The column `code` contains the population codes; the
columns `description` and `long description` report a brief and a long
description for the populations, respectively.

## See also

[`MutationEngine`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine_class.md)

## Examples

``` r
# build a mutation engine
m_engine <- MutationEngine(setup_code = "demo")
#> 
 [█---------------------------------------] 0% [00m:00s] Loading context index                                                                

 [████████████████████████████████████████] 100% [00m:00s] Context index loaded                                                               

#> 
 [█---------------------------------------] 0% [00m:00s] Loading RS index                                                                     

 [█████████████---------------------------] 32% [00m:01s] Loading RS index                                                                    

 [███████████████████████████-------------] 65% [00m:02s] Loading RS index                                                                    

 [███████████████████████████████████████-] 96% [00m:03s] Loading RS index                                                                    

 [████████████████████████████████████████] 100% [00m:03s] RS index loaded                                                                    

#> 
 [█---------------------------------------] 0% [00m:00s] Loading germline                                                                     

 [████████████████████████████████████████] 100% [00m:00s] Germline loaded                                                                    


# get the active germline subject data frame
head(m_engine$get_population_descriptions(), 5)
#>   code          description                    long.description
#> 1  CHB          Han Chinese       Han Chinese in Beijing, China
#> 2  JPT             Japanese            Japanese in Tokyo, Japan
#> 3  CHS Southern Han Chinese                   Han Chinese South
#> 4  CDX          Dai Chinese Chinese Dai in Xishuangbanna, China
#> 5  KHV      Kinh Vietnamese   Kinh in Ho Chi Minh City, Vietnam
```
