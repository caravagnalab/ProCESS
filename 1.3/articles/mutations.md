# Simulating Mutations

> *Disclaimer:* ProCESS/CLONES implements probability distributions
> using the C++11 random number distribution classes. Since the standard
> does not specify the underlying algorithms, their implementations are
> left to the compiler. Consequently, the simulation output depends on
> the compiler used to build
> [CLONES](https://github.com/albertocasagrande/CLONES), and the results
> reported in this article may differ from those obtained by the reader.

ProCESS/CLONES can simulate genomic mutations in the cells represented
in a `SampleForest` according to the specified SBS and indel mutational
signatures (see (Alexandrov et al. 2020)). This process is performed by
the class `MutationEngine`, which also takes into account the mutation
rate of the simulated species and gives the chance to dynamically change
the signatures.

### Setting Up Mutation Engine

The creation of an object of the type `MutationEngine` requires
downloading the reference sequence and the signature files and building
the corresponding two signature indices. The function
[`MutationEngine()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine.md)
performs all these steps in a single call.

This function is quite flexible and allows customising the mutation
engine in many ways. However, the vast majority of ProCESS users would
aim for a standard setup, for instance, involving the human genome.

Because of this, ProCESS provides some predefined setups. A list
containing them can be obtained by invoking the function
[`get_mutation_engine_codes()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_mutation_engine_codes.md).

[`library`](https://rdrr.io/r/base/library.html)`(`[`"ProCESS"`](https://caravagnalab.github.io/ProCESS/1.3)`)`` ``#> `` ``#> Attaching package: 'ProCESS'`` ``#> The following object is masked from 'package:utils':`` ``#> `` ``#> example`` `` `[`get_mutation_engine_codes`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_mutation_engine_codes.md)`(``)`` ``#> name description`` ``#> 1 GRCh37 Homo sapiens (GRCh37)`` ``#> 2 GRCh38 Homo sapiens (GRCh38)`` ``#> 3 demo A demonstrative set-up`

Let us build a `MutationEngine` by using the `"demo"` setup.

`# building a mutation engine by using the "demo" setup`` ``m_engine`` ``<-`` `[`MutationEngine`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine.md)`(``setup_code ``=`` ``"demo"``)`` ``#> Downloading reference genome...`` ``#> Reference genome downloaded`` ``#> Decompressing reference genome...done`` ``#> Downloading signature files...`` ``#> Signature file downloaded`` ``#> Downloading driver mutation file...`` ``#> Driver mutation file downloaded`` ``#> Decompressing driver mutation file...done`` ``#> Downloading passenger CNAs file...`` ``#> Passenger CNAs file downloaded`` ``#> Decompressing passenger CNAs file...done`` ``#> Downloading germline...`` ``#> Germline downloaded`` ``#> Decompressing mutations...`` ``#> done`` ``#> Building context index...`` ``#> [█---------------------------------------] 0% [00m:00s] Processing chr. 22 [█████████████████-----------------------] 40% [00m:00s] Processing chr. 22 [█████████████████████████████████-------] 81% [00m:02s] Processing chr. 22 [████████████████████████████████████████] 100% [00m:02s] Context index built`` ``#> [█---------------------------------------] 0% [00m:00s] Saving context index [████████████████████████████████████████] 100% [00m:00s] Context index saved`` ``#> done`` ``#> Building repeated sequence index...`` ``#> [█---------------------------------------] 0% [00m:00s] Reading 22 [████████████████████████████████████████] 100% [00m:01s] Reading 22`` ``#> [████████████████████████████████████████] 100% [00m:01s] Reading 22`` ``#> [████████████████████████████████████████] 100% [00m:01s] Reading 22`` ``#> [████████████████████████████████████████] 100% [00m:01s] Reading 22`` ``#> [████████████████████████████████████████] 100% [00m:01s] Reading 22`` ``#> [████████████████████████████████████████] 100% [00m:01s] Reading 22`` ``#> [████████████████████████████████████████] 100% [00m:01s] Reading 22`` ``#> [████████████████████████████████████████] 100% [00m:01s] RS index built`` ``#> [█---------------------------------------] 0% [00m:00s] Saving RS index [█---------------------------------------] 0% [00m:00s] Saving RS index [█████████████████████-------------------] 51% [00m:02s] Saving RS indexdone`` ``#> [████████████████████████████████████████] 100% [00m:02s] RS index saved`` ``#> [█---------------------------------------] 0% [00m:00s] Loading germline [████████████████████████████████████████] 100% [00m:00s] Germline loaded`` ``#> [█---------------------------------------] 0% [00m:00s] Saving germline [████████████████████████████████████████] 100% [00m:00s] Germline saved`` `` ``m_engine`` ``#> MutationEngine`` ``#> Passenger rates`` ``#> `` ``#> Driver mutations`` ``#> `` ``#> Timed Exposure`` ``#> SBS Timed Exposures`` ``#> `` ``#> indel Timed Exposures`

The above call creates the directory `demo`, downloads all the data
required by the mutation engine, and builds it.

[`dir.exists`](https://rdrr.io/r/base/files2.html)`(``"demo"``)`` ``#> [1] TRUE`` `` `[`list.files`](https://rdrr.io/r/base/list.files.html)`(``"demo"``)`` ``#> [1] "context_index_30.cif" "drivers.txt" `` ``#> [3] "germline_data" "indel_signatures.txt" `` ``#> [5] "passenger_CNAs.txt" "reference.fasta" `` ``#> [7] "rs_index_50_500000.rsif" "SBS_signatures.txt" `` ``#> [9] "sources.csv"`

The execution of the function
[`MutationEngine()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine.md)
may take some time, but it is meant to be performed one for all and, as
long as the user does not need to change the reference genome or the
signature files, it is no longer required. In this spirit, any call to
[`MutationEngine()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine.md)
exclusively performs the building process sub-tasks that have not been
already fulfilled.

`# building a mutation engine by using the "demo" set-up configuration`` ``m_engine`` ``<-`` `[`MutationEngine`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine.md)`(``setup_code ``=`` ``"demo"``)`` ``#> [█---------------------------------------] 0% [00m:00s] Loading context index [████████████████████████████████████████] 100% [00m:00s] Context index loaded`` ``#> [█---------------------------------------] 0% [00m:00s] Loading RS index [█████████████---------------------------] 32% [00m:01s] Loading RS index [██████████████████████████--------------] 64% [00m:02s] Loading RS index [██████████████████████████████████████--] 94% [00m:03s] Loading RS index [████████████████████████████████████████] 100% [00m:03s] RS index loaded`` ``#> [█---------------------------------------] 0% [00m:00s] Loading germline [████████████████████████████████████████] 100% [00m:00s] Germline loaded`

#### Mutant Genetic Characterization

Once the mutation engine has been built, we can define a mutant genotype
and declare species mutation rates.

Let us consider the simulation of Section “*Two populations with
epigenetic state*” in [this
article](https://caravagnalab.github.io/ProCESS/1.3/articles/sampling.md).
It involves the two mutants `A` and `B`. Both of them have two possible
epigenetic states, `E1` and `E2`, leading to the four species `A[E1]`,
`A[E2]`, `B[E1]`, and `B[E2]`, respectively. Each of these species has a
passenger mutation rate that must be declared to the mutation engine
before the sample forest is labelled. Let 1e-9, 3e-8, 8e-7, and 5e-8 be
the passenger SNV rates for the species `A[E1]`, `A[E2]`, `B[E1]`, and
`B[E2]`, respectively. Let the indel rates be 0 for all the species, but
`A[E1]`, which has an indel rate of 1e-8. Furthermore, let 0, 1e-11, 0,
and 0 be the passenger CNA rates of the same species, respectively.

The two mutants may also be genetically characterised by some driver
mutations. The driver mutations associated with each mutant must occur
in any cell belonging to that mutant. Hence, they must be declared to
the mutant engine before the labelling.

The method
[`MutationEngine$add_mutant()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-add_mutant.md)
takes care of all these declarations.

For the sake of example, let us assume that `A` is characterised by one
driver mutation on chromosome 22, while `B` has three driver mutations
on the same chromosome.

The genetic specification of the mutant `A` can be declared as it
follows.

`# add the mutant "A" characterized by one driver SNV on the allele 1 of`` ``# chromosome 22, one indel deletion on the same chromosome, two CNAs (a`` ``# deletion on the allele 1 of chromosome 22 and an amplification on a random`` ``# allele of the same chromosome) and, finally, a whole genome doubling event`` ``# (WGD). The mutant has two epigenetic states and its species "A[E1]" and`` ``# "A[E2]" have passenger SNV rates 1e-9 and 3e-9, respectively, passenger indel`` ``# rates 1e-9 and 0, respectively, and passenger CNA rates 0 and 1e-11,`` ``# respectively.`` ``m_engine``$``add_mutant``(``mutant_name ``=`` ``"A"``,`` `` passenger_rates ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``"E1"`` ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``SNV ``=`` ``1e-9``, indel ``=`` ``1e-9``)``,`` `` ``"E2"`` ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``SNV ``=`` ``3e-9``, CNA ``=`` ``1e-11``)``)``,`` `` drivers ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(`[`SNV`](https://caravagnalab.github.io/ProCESS/1.3/reference/SNV.md)`(``"22"``, ``46510210``, ``"C"``, allele ``=`` ``1``)``,`` `` `[`Mutation`](https://caravagnalab.github.io/ProCESS/1.3/reference/Mutation.md)`(``"22"``, ``16085675``, ``"GCCTCCCGA"``, ``"G"``)``,`` `` `[`CNA`](https://caravagnalab.github.io/ProCESS/1.3/reference/CNA.md)`(``type ``=`` ``"A"``, chr ``=`` ``"22"``,`` `` from ``=`` ``10303470``, len ``=`` ``200000``)``,`` `` `[`CNA`](https://caravagnalab.github.io/ProCESS/1.3/reference/CNA.md)`(``"D"``, ``"22"``, ``5010000``, ``200000``,`` `` allele ``=`` ``1``)``,`` `` ``WGD``)``)`` ``#> [█---------------------------------------] 0% [00m:00s] Retrieving "A" SIDs [█████████████████████-------------------] 50% [00m:03s] Retrieving "A" SIDs [████████████████████████████████████████] 100% [00m:03s] "A"'s SIDs validated`` `` ``m_engine`` ``#> MutationEngine`` ``#> Passenger rates`` ``#> "A[E1]":`` ``#> [0,inf): {SNV: 1e-09, indel: 1e-09}`` ``#> "A[E2]":`` ``#> [0,inf): {SNV: 3e-09, CNA: 1e-11}`` ``#> `` ``#> Driver mutations`` ``#> "A":`` ``#> (chr22(46510210)[A>C]) on allele 1`` ``#> (chr22(16085675)[GCCTCCCGA>G]) on random allele`` ``#> CNA("A",chr22(10303470), len: 200000)`` ``#> CNA("D",chr22(5010000), len: 200000, allele: 1)`` ``#> Whole genome duplication`` ``#> `` ``#> Timed Exposure`` ``#> SBS Timed Exposures`` ``#> `` ``#> indel Timed Exposures`

In the above code, the driver CNAs, SNVs, and indels are declared by
calling the functions
[`CNA()`](https://caravagnalab.github.io/ProCESS/1.3/reference/CNA.md),
[`SNV()`](https://caravagnalab.github.io/ProCESS/1.3/reference/SNV.md),
and
[`Mutation()`](https://caravagnalab.github.io/ProCESS/1.3/reference/Mutation.md),
respectively. These functions allow specifying the allele in which the
mutation must lie. The mutant `A` is also characterised by *whole genome
doubling* (`WGD`): a genomic event that simultaneously duplicates all
the chromosome alleles.

The driver mutations are applied in the order specified. For instance,
the following two pieces of code specify two different genomic
characterisations for the same mutant `E`.

`m_engine``$``add_mutant``(``"E"``, ``passenger_rates``,`` `` drivers ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(`[`SNV`](https://caravagnalab.github.io/ProCESS/1.3/reference/SNV.md)`(``"22"``, ``10510210``, ``"C"``, allele ``=`` ``1``)``,`` `` ``WGD``)``)`

`m_engine``$``add_mutant``(``"E"``, ``passenger_rates``,`` `` drivers ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``WGD``,`` `` `[`SNV`](https://caravagnalab.github.io/ProCESS/1.3/reference/SNV.md)`(``"22"``, ``10510210``, ``"C"``, allele ``=`` ``1``)``)``)`

The former snippet places an SNV and, afterwards, duplicates all the
alleles producing a genome in which two alleles contain the placed SNV.
Instead, the latter snippet requires a whole genome doubling event and,
then, places the SNV yielding a single occurrence of the SNV in the
final genome.

As far as SNVs and indels are concerned, ProCESS provides users with a
more compact and, at times, convenient approach. The mutation engine
stores a list of known driver mutations and labels each with a code that
can be used during mutant declaration. The corresponding data frame can
be obtained by using the method
[`MutationEngine$get_known_drivers()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-get_known_drivers.md).

[`library`](https://rdrr.io/r/base/library.html)`(`[`dplyr`](https://dplyr.tidyverse.org)`)`` ``#> `` ``#> Attaching package: 'dplyr'`` ``#> The following objects are masked from 'package:stats':`` ``#> `` ``#> filter, lag`` ``#> The following objects are masked from 'package:base':`` ``#> `` ``#> intersect, setdiff, setequal, union`` `` ``m_engine``$``get_known_drivers``(``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``chr`` ``==`` ``"22"``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `[`head`](https://rdrr.io/r/utils/head.html)`(``)`` ``#> chr from to ref alt mutation_type driver_gene driver_code`` ``#> 1 22 20073503 20073503 G A SNV DGCR8 DGCR8 S6N`` ``#> 2 22 20073503 20073503 G A SNV DGCR8 DGCR8 S6N`` ``#> 3 22 20073503 20073503 G A SNV DGCR8 DGCR8 S6N`` ``#> 4 22 20073503 20073503 G A SNV DGCR8 DGCR8 S6N`` ``#> 5 22 20073511 20073511 C T SNV DGCR8 DGCR8 P9S`` ``#> 6 22 20073511 20073511 C T SNV DGCR8 DGCR8 P9S`` ``#> driver_CDS tumour_type`` ``#> 1 c.17G>A DLBCLNOS`` ``#> 2 c.17G>A PLMESO`` ``#> 3 c.17G>A WDTC`` ``#> 4 c.17G>A WT`` ``#> 5 c.25C>T DLBCLNOS`` ``#> 6 c.25C>T PLMESO`

The code of the known driver mutations can be used in place of the full
specification as follows.

`# add the mutant "B" characterized by two driver SNVs on chromosome 22 (no`` ``# CNA) and two epigenetic states. The first SNV is "DGCR8 P26L" and it must`` ``# lay in the allele 1. The remaining SNV is specified by using the SNV function`` ``# as done above. The species "B[E1]" and "B[E2]" have passenger SNV rates 8e-9`` ``# and 2e-8, respectively, and CNA rates 0 for both species.`` ``m_engine``$``add_mutant``(``"B"``, `[`list`](https://rdrr.io/r/base/list.html)`(``"E1"`` ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``SNV ``=`` ``8e-9``)``, ``"E2"`` ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``SNV ``=`` ``2e-8``)``)``,`` `` `[`list`](https://rdrr.io/r/base/list.html)`(`[`list`](https://rdrr.io/r/base/list.html)`(``"DGCR8 P26L"``, allele ``=`` ``1``)``, ``# the first SNV`` `` `[`SNV`](https://caravagnalab.github.io/ProCESS/1.3/reference/SNV.md)`(``"22"``, ``12028576``, ``"G"``)``)``)`` ``# the second SNV`` ``#> [█---------------------------------------] 0% [00m:00s] Retrieving "B" SIDs [████████████████████████████████████████] 100% [00m:00s] "B"'s SIDs validated`` `` ``m_engine`` ``#> MutationEngine`` ``#> Passenger rates`` ``#> "A[E1]":`` ``#> [0,inf): {SNV: 1e-09, indel: 1e-09}`` ``#> "A[E2]":`` ``#> [0,inf): {SNV: 3e-09, CNA: 1e-11}`` ``#> "B[E1]":`` ``#> [0,inf): {SNV: 8e-09}`` ``#> "B[E2]":`` ``#> [0,inf): {SNV: 2e-08}`` ``#> `` ``#> Driver mutations`` ``#> "A":`` ``#> (chr22(46510210)[A>C]) on allele 1`` ``#> (chr22(16085675)[GCCTCCCGA>G]) on random allele`` ``#> CNA("A",chr22(10303470), len: 200000)`` ``#> CNA("D",chr22(5010000), len: 200000, allele: 1)`` ``#> Whole genome duplication`` ``#> "B":`` ``#> DGCR8 P26L (chr22(20073563)[C>T]) on allele 1`` ``#> (chr22(12028576)[N>G]) on random allele`` ``#> `` ``#> Timed Exposure`` ``#> SBS Timed Exposures`` ``#> `` ``#> indel Timed Exposures`

#### Mutational Exposures

The probability of a mutation occurring depends on both its genomic and
environmental context.

A *signature* is a mutation probability distribution over mutation
contexts or mutation structure. *SBS (single-base substitution)
signatures* provide for any genomic context (i.e., a triplet of bases)
the probability that an SNV occurs on that context. On the contrary, *ID
(indel) signatures* associate the probability of an indel to its length
and structure (see (Alexandrov et al. 2020)).

The signature depends on the environmental context and, as a result,
more than one signature may be active at the same time with different
probabilities. A *mutational exposure* (or *exposure*) is a discrete
probability distribution among signatures.

To simulate passenger mutations of a given type, we need to specify a
default exposure for that type. This can be achieved as follows.

`# add SNV and indel default exposures. This will be used from simulated time 0`` ``# up to the successive exposure change.`` ``m_engine``$``add_exposure``(``coefficients ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``SBS13 ``=`` ``0.2``, SBS1 ``=`` ``0.8``)``)`` ``m_engine``$``add_exposure``(`[`c`](https://rdrr.io/r/base/c.html)`(``ID2 ``=`` ``0.6``, ID13 ``=`` ``0.2``, ID21 ``=`` ``0.2``)``)`

Further exposures can also be defined by specifying an activation time
for each of them, i.e., the new exposures will be used from that until
the next exposure change.

`# add a new SNV exposure that will be used from simulated`` ``# time 100 up to the successive exposure change.`` ``m_engine``$``add_exposure``(``time ``=`` ``100``, `[`c`](https://rdrr.io/r/base/c.html)`(``SBS5 ``=`` ``0.3``, SBS2 ``=`` ``0.2``, SBS3 ``=`` ``0.5``)``)`` `` ``m_engine`` ``#> MutationEngine`` ``#> Passenger rates`` ``#> "A[E1]":`` ``#> [0,inf): {SNV: 1e-09, indel: 1e-09}`` ``#> "A[E2]":`` ``#> [0,inf): {SNV: 3e-09, CNA: 1e-11}`` ``#> "B[E1]":`` ``#> [0,inf): {SNV: 8e-09}`` ``#> "B[E2]":`` ``#> [0,inf): {SNV: 2e-08}`` ``#> `` ``#> Driver mutations`` ``#> "A":`` ``#> (chr22(46510210)[A>C]) on allele 1`` ``#> (chr22(16085675)[GCCTCCCGA>G]) on random allele`` ``#> CNA("A",chr22(10303470), len: 200000)`` ``#> CNA("D",chr22(5010000), len: 200000, allele: 1)`` ``#> Whole genome duplication`` ``#> "B":`` ``#> DGCR8 P26L (chr22(20073563)[C>T]) on allele 1`` ``#> (chr22(12028576)[N>G]) on random allele`` ``#> `` ``#> Timed Exposure`` ``#> SBS Timed Exposures`` ``#> [0, 100[: {"SBS1": 0.8, "SBS13": 0.2}`` ``#> [100, ∞[: {"SBS2": 0.2, "SBS3": 0.5, "SBS5": 0.3}`` ``#> `` ``#> indel Timed Exposures`` ``#> [0, ∞[: {"ID13": 0.2, "ID2": 0.6, "ID21": 0.2}`

The default exposures and the simultaneous changes in SNV and ID
exposures can be specified by a single call to the function
[`MutationEngine$add_exposure()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-add_exposure.md)
as follows.

`m_engine``$``add_exposure``(``time ``=`` ``120``, `[`c`](https://rdrr.io/r/base/c.html)`(``SBS5 ``=`` ``0.3``, SBS2 ``=`` ``0.2``, SBS3 ``=`` ``0.5``,`` `` ID1 ``=`` ``0.8``, ID9 ``=`` ``0.2``)``)`` `` ``m_engine`` ``#> MutationEngine`` ``#> Passenger rates`` ``#> "A[E1]":`` ``#> [0,inf): {SNV: 1e-09, indel: 1e-09}`` ``#> "A[E2]":`` ``#> [0,inf): {SNV: 3e-09, CNA: 1e-11}`` ``#> "B[E1]":`` ``#> [0,inf): {SNV: 8e-09}`` ``#> "B[E2]":`` ``#> [0,inf): {SNV: 2e-08}`` ``#> `` ``#> Driver mutations`` ``#> "A":`` ``#> (chr22(46510210)[A>C]) on allele 1`` ``#> (chr22(16085675)[GCCTCCCGA>G]) on random allele`` ``#> CNA("A",chr22(10303470), len: 200000)`` ``#> CNA("D",chr22(5010000), len: 200000, allele: 1)`` ``#> Whole genome duplication`` ``#> "B":`` ``#> DGCR8 P26L (chr22(20073563)[C>T]) on allele 1`` ``#> (chr22(12028576)[N>G]) on random allele`` ``#> `` ``#> Timed Exposure`` ``#> SBS Timed Exposures`` ``#> [0, 100[: {"SBS1": 0.8, "SBS13": 0.2}`` ``#> [100, 120[: {"SBS2": 0.2, "SBS3": 0.5, "SBS5": 0.3}`` ``#> [120, ∞[: {"SBS2": 0.2, "SBS3": 0.5, "SBS5": 0.3}`` ``#> `` ``#> indel Timed Exposures`` ``#> [0, 120[: {"ID13": 0.2, "ID2": 0.6, "ID21": 0.2}`` ``#> [120, ∞[: {"ID1": 0.8, "ID9": 0.2}`

#### Passenger CNAs and Tumour Type

The passenger CNAs applied during the simulation depend on tumour type.
The type of tumour can be specified when building the mutation engine by
using the parameters `tumour_type` or `tumour_study`. For instance, if
the passenger CNA file used to build the mutation engine contains some
of the data identified in breast carcinoma in the UK, the specific CNA
can be applied during the simulation as follows.

`m_engine`` ``<-`` `[`MutationEngine`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine.md)`(``setup_code ``=`` ``"demo"``,`` `` tumour_type ``=`` ``"BRCA"``,`` `` tumour_study ``=`` ``"UK"``)`

The complete list of available tumour types in a set-up can be obtained
by using the function
[`get_available_tumours_in()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_available_tumours_in.md).

[`get_available_tumours_in`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_available_tumours_in.md)`(``setup_code ``=`` ``"demo"``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `[`head`](https://rdrr.io/r/utils/head.html)`(``)`` ``#> type`` ``#> 1 ACC`` ``#> 2 AML`` ``#> 3 ANGS`` ``#> 4 ANSC`` ``#> 5 BCC`` ``#> 6 BLCA`

#### Germline Mutations

ProCESS allows users to apply germline mutations from one of the
subjects in the germline data used to build the mutation engine. This
feature will enable users to simulate a specific cancer type’s evolution
on an individual with the desired gender and ethnicity.

The available subjects, together with their sex and ethnicity, can be
obtained by using the method
[`MutationEngine$get_germline_subjects()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-get_germline_subjects.md).

`subjects`` ``<-`` ``m_engine``$``get_germline_subjects``(``)`` `` ``subjects`` ``#> sample pop super_pop gender`` ``#> 1 NA18941 JPT EAS female`` ``#> 2 NA20513 TSI EUR male`

The column `sample` contains the names of the available subjects. The
columns `pop` and `super_pop` report the subjects’ population and
super-population codes. The last column, `gender`, includes the subject
gender.

The method
[`MutationEngine$get_population_descriptions()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-get_population_descriptions.md)
clarifies the meaning of the codes reported in the pop columns.

`m_engine``$``get_population_descriptions``(``)`` ``#> code description`` ``#> 1 CHB Han Chinese`` ``#> 2 JPT Japanese`` ``#> 3 CHS Southern Han Chinese`` ``#> 4 CDX Dai Chinese`` ``#> 5 KHV Kinh Vietnamese`` ``#> 6 CHD Denver Chinese`` ``#> 7 CEU CEPH`` ``#> 8 TSI Tuscan`` ``#> 9 GBR British`` ``#> 10 FIN Finnish`` ``#> 11 IBS Spanish`` ``#> 12 YRI Yoruba`` ``#> 13 LWK Luhya`` ``#> 14 GWD Gambian`` ``#> 15 MSL Mende`` ``#> 16 ESN Esan`` ``#> 17 ASW African-American SW`` ``#> 18 ACB African-Caribbean`` ``#> 19 MXL Mexican-American`` ``#> 20 PUR Puerto Rican`` ``#> 21 CLM Colombian`` ``#> 22 PEL Peruvian`` ``#> 23 GIH Gujarati`` ``#> 24 PJL Punjabi`` ``#> 25 BEB Bengali`` ``#> 26 STU Sri Lankan`` ``#> 27 ITU Indian`` ``#> long.description`` ``#> 1 Han Chinese in Beijing, China`` ``#> 2 Japanese in Tokyo, Japan`` ``#> 3 Han Chinese South`` ``#> 4 Chinese Dai in Xishuangbanna, China`` ``#> 5 Kinh in Ho Chi Minh City, Vietnam`` ``#> 6 Chinese in Denver, Colorado (pilot 3 only)`` ``#> 7 Utah residents (CEPH) with Northern and Western European ancestry `` ``#> 8 Toscani in Italia `` ``#> 9 British in England and Scotland `` ``#> 10 Finnish in Finland `` ``#> 11 Iberian populations in Spain `` ``#> 12 Yoruba in Ibadan, Nigeria`` ``#> 13 Luhya in Webuye, Kenya`` ``#> 14 Gambian in Western Division, The Gambia `` ``#> 15 Mende in Sierra Leone`` ``#> 16 Esan in Nigeria`` ``#> 17 African Ancestry in Southwest US `` ``#> 18 African Caribbean in Barbados`` ``#> 19 Mexican Ancestry in Los Angeles, California`` ``#> 20 Puerto Rican in Puerto Rico`` ``#> 21 Colombian in Medellin, Colombia`` ``#> 22 Peruvian in Lima, Peru`` ``#> 23 Gujarati Indian in Houston, TX`` ``#> 24 Punjabi in Lahore, Pakistan`` ``#> 25 Bengali in Bangladesh`` ``#> 26 Sri Lankan Tamil in the UK`` ``#> 27 Indian Telugu in the UK`

The method
[`MutationEngine$get_active_germline()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-get_active_germline.md)
returns the the active germline subject.

`m_engine``$``get_active_germline``(``)`` ``#> sample pop super_pop gender`` ``#> 1 NA18941 JPT EAS female`

Users can change the germline subject using the method
[`MutationEngine$set_germline_subject()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-set_germline_subject.md).

When a subject is selected for the first time, ProCESS builds a binary
representation of the subject’s genome and save it for future use. This
step may take a few minutes. However, all the successive selections of
the same subject directly load the binary file.

`m_engine``$``set_germline_subject``(``subjects``[``2``, ``"sample"``]``)`` ``#> [█---------------------------------------] 0% [00m:00s] Loading germline [████████████████████████████████████████] 100% [00m:00s] Germline loaded`` `` ``m_engine``$``get_active_germline``(``)`` ``#> sample pop super_pop gender`` ``#> 1 NA20513 TSI EUR male`

### Building Phylogenetic Forests

The configured mutation engine can be used to label each node in a
sample forest with mutations.

Since the mutation engine has been configured to deal with the
simulation performed in Section “*Two populations with epigenetic
state*” in [this
article](https://caravagnalab.github.io/ProCESS/1.3/articles/sampling.md),
we can load its sample forest from the file `"sample_forest.sff"` as
saved there.

`sample_forest`` ``<-`` `[`load_forest`](https://caravagnalab.github.io/ProCESS/1.3/reference/load_forest.md)`(``"sample_forest.sff"``)`` ``#> [█---------------------------------------] 0% [00m:00s] Loading forest [████████████████████████████████████████] 100% [00m:00s] Forest loaded`` `` ``# place mutations on the sample forest assuming 1000 pre-neoplastic SNVs and`` ``# 500 indels`` ``phylo_forest`` ``<-`` ``m_engine``$``place_mutations``(``sample_forest``, ``1000``, ``500``)`` ``#> [█---------------------------------------] 0% [00m:00s] Placing mutations [████████████████████████████████████████] 100% [00m:00s] Mutations placed`` `` ``phylo_forest`` ``#> PhylogeneticForest`` ``#> # of trees: 1`` ``#> # of nodes: 21182`` ``#> # of leaves: 5425`` ``#> samples: {"S_1_1", "S_1_2", "S_2_1", "S_2_2"}`` ``#> `` ``#> # of emerged SNVs and indels: 16731`` ``#> # of emerged CNAs: 37`

The phylogenetic forest stores all mutations, labelling the sampled
cells represented by its leaves. Users can retrieve such data by using
the methods
[`PhylogeneticForest``get_sampled_cell_mutations()](../reference/PhylogeneticForest-cash-get_sampled_cell_mutations.md)</code> and <code>[PhylogeneticForest``get_sampled_cell_CNAs()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_sampled_cell_CNAs.md).

[`library`](https://rdrr.io/r/base/library.html)`(`[`dplyr`](https://dplyr.tidyverse.org)`)`` `` ``# select the first mutations among all the mutations occurring in`` ``# the genomes of the sampled cells`` ``phylo_forest``$``get_sampled_cell_mutations``(``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `[`head`](https://rdrr.io/r/utils/head.html)`(``)`` ``#> chr from allele ref alt cause nature cell_id`` ``#> 1 22 16133982 0 T A SBS1 pre-neoplastic 23281`` ``#> 2 22 16270440 0 T A SBS1 pre-neoplastic 23281`` ``#> 3 22 16277531 0 C CT ID1 pre-neoplastic 23281`` ``#> 4 22 16391790 0 A C SBS1 pre-neoplastic 23281`` ``#> 5 22 16440524 0 T C SBS1 pre-neoplastic 23281`` ``#> 6 22 16507752 0 A T SBS1 pre-neoplastic 23281`` `` ``# select the first CNAs among all the mutations occurring in`` ``# the genomes of the sampled cells`` ``phylo_forest``$``get_sampled_cell_CNAs``(``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `[`head`](https://rdrr.io/r/utils/head.html)`(``)`` ``#> chr begin end type allele src.allele cause nature cell_id`` ``#> 1 22 10303470 10503469 A 2 0 driver 23281`` ``#> 2 22 5010000 5209999 D 1 NA driver 23281`` ``#> 3 22 10303470 10503469 A 2 0 driver 56990`` ``#> 4 22 5010000 5209999 D 1 NA driver 56990`` ``#> 5 22 10303470 10503469 A 2 0 driver 92307`` ``#> 6 22 5010000 5209999 D 1 NA driver 92307`` `` ``# get the sampled cells`` ``sampled_cells`` ``<-`` ``phylo_forest``$``get_nodes``(``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `` `[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``!`[`is.na`](https://rdrr.io/r/base/NA.html)`(``.data``$``sample``)``)`` `` ``# show the first of them`` ``sampled_cells`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `[`head`](https://rdrr.io/r/utils/head.html)`(``)`` ``#> cell_id ancestor depth mutant epistate sample birth_time`` ``#> 1 4263 1988 14 A E2 S_1_1 132.3400`` ``#> 2 9258 9136 34 A E2 S_1_1 168.2328`` ``#> 3 10018 9148 19 A E2 S_1_1 172.3136`` ``#> 4 11446 10053 30 A E2 S_1_1 179.5106`` ``#> 5 11720 10386 36 A E2 S_1_1 180.6825`` ``#> 6 11764 10253 43 A E2 S_1_1 180.9052`` `` ``` # get the identifier of the 3rd cell in `sampled_cells` ``` ``cell_id`` ``<-`` ``sampled_cells``[``3``, ``1``]`` `` ``` # get the node whose associated cell have `cell_id` as identifier ``` ``node`` ``<-`` ``phylo_forest``$``get_node``(``cell_id``)`` ``node`` ``#> PhylogeneticForestNode(cell_id = 10018, species = "A[E2]")`` `` ``` # compute the genome of the cell having `cell_id` as identifier ``` ``genome`` ``<-`` ``node``$``get_genome``(``)`` ``genome`` ``#> GenomeMutations: 1 chrs 6 alleles`` `` ``` # get the SNVs and the indels of the 3rd cell in `sampled_cells` ``` ``genome``$``get_mutations``(``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `[`head`](https://rdrr.io/r/utils/head.html)`(``)`` ``#> chr from allele ref alt cause nature`` ``#> 1 22 16133982 0 T A SBS1 pre-neoplastic`` ``#> 2 22 16270440 0 T A SBS1 pre-neoplastic`` ``#> 3 22 16277531 0 C CT ID1 pre-neoplastic`` ``#> 4 22 16391790 0 A C SBS1 pre-neoplastic`` ``#> 5 22 16440524 0 T C SBS1 pre-neoplastic`` ``#> 6 22 16507752 0 A T SBS1 pre-neoplastic`` `` ``` # get the CNAs of the 3rd cell in `sampled_cells` ``` ``genome``$``get_CNAs``(``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `[`head`](https://rdrr.io/r/utils/head.html)`(``)`` ``#> chr begin end type allele src.allele cause nature`` ``#> 1 22 10303470 10503469 A 2 0 driver`` ``#> 2 22 5010000 5209999 D 1 NA driver`

The method
[`PhylogeneticForest``get_samples_info()](../reference/PhylogeneticForest-cash-get_samples_info.md)</code> is analogous to <code>[SampleForest``get_samples_info()`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForest-cash-get_samples_info.md):
it returns a data frame containing information about the forest samples.
However, the former method adds two columns to the result of the latter:
`DNA_quantity` and `equivalent_normal_cells`. The former column reports
the total quantity of tumoral DNA in the sample, i.e., the sum of the
lengths of all alleles in the sample. This quantity is a natural number,
nevertheless, it is stored as a real number as it usually exceeds the
largest natural number that can be natively represented by R. The column
“`equivalent_normal_cells`” instead contains the number of normal cells
that contain as much DNA as the sample tumour cells.

`phylo_forest``$``get_samples_info``(``)`` ``#> name id xmin ymin xmax ymax tumour_cells tumour_cells_in_bbox time`` ``#> 1 S_1_1 0 480 480 530 530 2576 2576 445.2994`` ``#> 2 S_1_2 1 500 500 550 550 1609 1609 445.2994`` ``#> 3 S_2_1 2 399 540 423 564 620 620 506.0268`` ``#> 4 S_2_2 3 549 365 573 389 620 620 506.0268`` ``#> DNA_quantity equivalent_normal_cells`` ``#> 1 531196062192 5176.889`` ``#> 2 330425908811 3220.239`` ``#> 3 127222173306 1239.872`` ``#> 4 127235323680 1240.000`

The method
[`PhylogeneticForest$get_germline_mutations()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_germline_mutations.md)
returns the SNVs and the indels in the germline.

`# extract the germline mutation`` ``phylo_forest``$``get_germline_mutations``(``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `[`head`](https://rdrr.io/r/utils/head.html)`(``)`` ``#> chr from allele ref alt cause nature`` ``#> 1 22 16051493 0 G A germinal`` ``#> 2 22 16052167 0 A AAAAC germinal`` ``#> 3 22 16053659 0 A C germinal`` ``#> 4 22 16054740 0 A G germinal`` ``#> 5 22 16055942 0 C T germinal`` ``#> 6 22 16058070 0 A G germinal`

Users can also identify the cell in which a mutation emerged even when
the cell was not sampled.

`` # get the genome of the cell having as the identifier `cell_id` ``` ``genome`` ``<-`` ``phylo_forest``$``get_node``(``cell_id``)``$``get_genome``(``)`` `` ``# select one of the mutations`` ``mutation_row`` ``<-`` ``genome``$``get_mutations``(``with_germline``=``FALSE``)``[``2``, ``]`` `` ``# rebuild the corresponding mutation`` ``mutation`` ``<-`` `[`Mutation`](https://caravagnalab.github.io/ProCESS/1.3/reference/Mutation.md)`(``mutation_row``[``"chr"``]``[``1``, ``]``,`` `` ``mutation_row``[``"from"``]``[``1``, ``]``,`` `` ``mutation_row``[``"ref"``]``[``1``, ``]``,`` `` ``mutation_row``[``"alt"``]``[``1``, ``]``)`` `` ``# get the identifier of the oldest cells in which the mutation occurs`` ``phylo_forest``$``get_first_occurrences``(``mutation``)`` ``#> [[1]]`` ``#> [1] 1`

The exposures used in placing the mutations on the cells in the
phylogenetic forest can be obtained by using the method
[`PhylogeneticForest$get_exposures()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_exposures.md).

`# get the exposures used in placing the mutations`` ``phylo_forest``$``get_exposures``(``)`` ``#> time signature exposure type`` ``#> 1 0 SBS1 0.8 SNV`` ``#> 2 0 SBS13 0.2 SNV`` ``#> 3 100 SBS2 0.2 SNV`` ``#> 4 100 SBS3 0.5 SNV`` ``#> 5 100 SBS5 0.3 SNV`` ``#> 6 120 SBS2 0.2 SNV`` ``#> 7 120 SBS3 0.5 SNV`` ``#> 8 120 SBS5 0.3 SNV`` ``#> 9 0 ID13 0.2 indel`` ``#> 10 0 ID2 0.6 indel`` ``#> 11 0 ID21 0.2 indel`` ``#> 12 120 ID1 0.8 indel`` ``#> 13 120 ID9 0.2 indel`

The method
[`MutationEngine$get_bulk_allelic_fragmentation()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-get_bulk_allelic_fragmentation.md)
returns a data frame reporting the allelic type per genome fragment.

`# get the name of the first sample`` ``sample_name`` ``<-`` ``phylo_forest``$``get_samples_info``(``)``[[``"name"``]``]``[``1``]`` `` ``# print the bulk allelic fragmentation`` ``phylo_forest``$``get_bulk_allelic_fragmentation``(``sample_name``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `[`head`](https://rdrr.io/r/utils/head.html)`(``)`` ``#> chr begin end major minor ratio`` ``#> 1 22 1 5009999 2 2 1.0000000`` ``#> 2 22 5010000 5209999 2 0 1.0000000`` ``#> 3 22 5210000 10303469 2 2 1.0000000`` ``#> 4 22 10303470 10503469 4 2 1.0000000`` ``#> 5 22 10503470 17098700 2 2 1.0000000`` ``#> 6 22 17098701 17899999 2 2 0.9996118`

Instead, the method
[`MutationEngine$get_cell_allelic_fragmentation()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-get_cell_allelic_fragmentation.md)
returns the allelic fragmentation per cell.

`# print the cell allelic fragmentation`` ``phylo_forest``$``get_cell_allelic_fragmentation``(``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `[`head`](https://rdrr.io/r/utils/head.html)`(``)`` ``#> cell_id chr begin end major minor`` ``#> 1 23281 22 1 5009999 2 2`` ``#> 2 23281 22 5010000 5209999 2 0`` ``#> 3 23281 22 5210000 10303469 2 2`` ``#> 4 23281 22 10303470 10503469 4 2`` ``#> 5 23281 22 10503470 51304566 2 2`` ``#> 6 56990 22 1 5009999 2 2`

The details about the SNV and indel signatures adopted during the
evolution are available in the mutation engine and they can be retrieved
by using the methods
[`MutationEngine``get_SNV_signatures()](../reference/MutationEngine-cash-get_SNV_signatures.md)</code> and <code>[MutationEngine``get_indel_signatures()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-get_indel_signatures.md).

`# get the SNV signatures used in placing the mutations`` ``m_engine``$``get_SNV_signatures``(``)``[``1``:``6``, ``1``:``5``]`` ``#> Type SBS1 SBS2 SBS3 SBS4`` ``#> 1 A[C>A]A 0.01041667 0.01041667 0.01041667 0.01041667`` ``#> 2 A[C>A]C 0.01041667 0.01041667 0.01041667 0.01041667`` ``#> 3 A[C>A]G 0.01041667 0.01041667 0.01041667 0.01041667`` ``#> 4 A[C>A]T 0.01041667 0.01041667 0.01041667 0.01041667`` ``#> 5 A[C>G]A 0.01041667 0.01041667 0.01041667 0.01041667`` ``#> 6 A[C>G]C 0.01041667 0.01041667 0.01041667 0.01041667`` `` ``# get the indel signatures used in placing the mutations`` ``m_engine``$``get_indel_signatures``(``)``[``1``:``6``, ``1``:``5``]`` ``#> Type ID1 ID2 ID3 ID4`` ``#> 1 1:Del:C:0 0.01204819 0.01204819 0.01204819 0.01204819`` ``#> 2 1:Del:C:1 0.01204819 0.01204819 0.01204819 0.01204819`` ``#> 3 1:Del:C:2 0.01204819 0.01204819 0.01204819 0.01204819`` ``#> 4 1:Del:C:3 0.01204819 0.01204819 0.01204819 0.01204819`` ``#> 5 1:Del:C:4 0.01204819 0.01204819 0.01204819 0.01204819`` ``#> 6 1:Del:C:5 0.01204819 0.01204819 0.01204819 0.01204819`

Finally, the data of the subject whose germline corresponds to wild-type
genome in the phylogenetic forest can be obtained by the method
[`PhylogeneticForest$get_germline_subject()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_germline_subject.md).

`phylo_forest``$``get_germline_subject``(``)`` ``#> sample pop super_pop gender`` ``#> 1 NA20513 TSI EUR male`

### Storing Phylogenetic Forests

As in the case of the sample forests, the phylogenetic forests can be
saved by using the method
[`PhylogeneticForest$save()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-save.md)
and load by the function
[`load_forest()`](https://caravagnalab.github.io/ProCESS/1.3/reference/load_forest.md).

`# save the phylogenetic forest in the file "phylo_forest.pff"`` ``phylo_forest``$``save``(``"phylo_forest.pff"``)`` ``#> [█---------------------------------------] 0% [00m:00s] Saving forest [████████████████████████████████████████] 100% [00m:00s] Forest saved`` `` ``# loading the saved forest`` ``loaded_phylo_forest`` ``<-`` `[`load_forest`](https://caravagnalab.github.io/ProCESS/1.3/reference/load_forest.md)`(``"phylo_forest.pff"``)`` ``#> [█---------------------------------------] 0% [00m:00s] Loading forest [████████████████████████████████████████] 100% [00m:00s] Forest loaded`` `` ``loaded_phylo_forest`` ``#> PhylogeneticForest`` ``#> # of trees: 1`` ``#> # of nodes: 21182`` ``#> # of leaves: 5425`` ``#> samples: {"S_1_1", "S_1_2", "S_2_1", "S_2_2"}`` ``#> `` ``#> # of emerged SNVs and indels: 16731`` ``#> # of emerged CNAs: 37`

### Getting and Setting the Reference Genome Path

The phylogenetic forest object contains the reference genome FASTA file
path. The methods
[`PhylogeneticForest``get_reference_path()](../reference/PhylogeneticForest-cash-get_reference_path.md)</code> and <code>[PhylogeneticForest``set_reference_path()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-set_reference_path.md)
can be used to get the path and set it, respectively.

`phylo_forest``$``get_reference_path``(``)`` ``#> [1] "/Users/alberto/Library/CloudStorage/Dropbox/Lavoro/Code/ProCESS/ProCESS/vignettes/demo/reference.fasta"`` `` ``phylo_forest``$``set_reference_path``(``"demo/reference.fasta"``)`` `` ``phylo_forest``$``get_reference_path``(``)`` ``#> [1] "/Users/alberto/Library/CloudStorage/Dropbox/Lavoro/Code/ProCESS/ProCESS/vignettes/demo/reference.fasta"`

Alexandrov, Ludmil B, Jaegil Kim, Nicholas J Haradhvala, et al. 2020.
“The Repertoire of Mutational Signatures in Human Cancer.” *Nature* 578
(7793): 94–101.
