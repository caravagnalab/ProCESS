# Customizing mutation engine

> *Disclaimer:* ProCESS/CLONES implements probability distributions
> using the C++11 random number distribution classes. Since the standard
> does not specify the underlying algorithms, their implementations are
> left to the compiler. Consequently, the simulation output depends on
> the compiler used to build
> [CLONES](https://github.com/albertocasagrande/CLONES), and the results
> reported in this article may differ from those obtained by the reader.

ProCESS/CLONES can simulate genomic mutations in the cells represented
in a `SampleForest` according to the specified SBS and indel mutational
signatures (see the [article about
mutations](https://caravagnalab.github.io/ProCESS/1.3/articles/mutations.md)).
This task is achieved by the mutation engine, which uses two indices on
the reference genome to place the most probable mutations according to
signatures: the *context index* and the *repeated sequence index*.

The function
[`MutationEngine()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine.md)
builds the mutation engine and these two indices.

## Customizing the mutation engine

While the most convenient way to invoke
[`MutationEngine()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine.md)
is by using the parameter `setup_code` (see the [article about
mutations](https://caravagnalab.github.io/ProCESS/1.3/articles/mutations.md)),
there are cases in which a custom setup is required. In these cases, the
function
[`MutationEngine()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine.md)
can be called by specifying the name of the setup directory, the path or
URL of the reference sequence, the signature files, the driver mutation
file, the passenger CNAs file, and the germline data directory.

[`library`](https://rdrr.io/r/base/library.html)`(`[`"ProCESS"`](https://caravagnalab.github.io/ProCESS/1.3)`)`` ``#> `` ``#> Attaching package: 'ProCESS'`` ``#> The following object is masked from 'package:utils':`` ``#> `` ``#> example`` `` ``reference_url`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"https://ftp.ensembl.org/pub/grch37/release-111/"``,`` `` ``"fasta/homo_sapiens/dna/Homo_sapiens.GRCh37."``,`` `` ``"dna.chromosome.22.fa.gz"``)`` `` ``sbs_url`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"https://zenodo.org/records/15875185/files/"``,`` `` ``"SBS_demo_signatures.txt"``)`` `` ``indel_url`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"https://zenodo.org/records/15875185/files/"``,`` `` ``"indel_demo_signatures.txt"``)`` `` ``drivers_url`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"https://zenodo.org/records/15875185/files/"``,`` `` ``"driver_mutations_hg19.csv.bz2"``)`` `` ``passenger_cnas_url`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"https://zenodo.org/records/15875185/files/"``,`` `` ``"passenger_CNAs_hg19.csv.bz2"``)`` `` ``germline_url`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"https://zenodo.org/records/15875185/files/"``,`` `` ``"germline_data_demo.tar.bz2"``)`` `` ``# build a mutation engine and place all the files in the directory "Test"`` ``m_engine`` ``<-`` `[`MutationEngine`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine.md)`(``directory ``=`` ``"Test"``,`` `` reference_src ``=`` ``reference_url``,`` `` SBS_signatures_src ``=`` ``sbs_url``,`` `` indel_signatures_src ``=`` ``indel_url``,`` `` drivers_src ``=`` ``drivers_url``,`` `` passenger_CNAs_src ``=`` ``passenger_cnas_url``,`` `` germline_src ``=`` ``germline_url``)`` ``#> Downloading reference genome...`` ``#> Reference genome downloaded`` ``#> Decompressing reference genome...done`` ``#> Downloading signature files...`` ``#> Signature file downloaded`` ``#> Downloading driver mutation file...`` ``#> Driver mutation file downloaded`` ``#> Decompressing driver mutation file...done`` ``#> Downloading passenger CNAs file...`` ``#> Passenger CNAs file downloaded`` ``#> Decompressing passenger CNAs file...done`` ``#> Downloading germline...`` ``#> Germline downloaded`` ``#> Decompressing mutations...`` ``#> done`` ``#> Building context index...`` ``#> [█---------------------------------------] 0% [00m:00s] Processing chr. 22 [█████████████████-----------------------] 40% [00m:00s] Processing chr. 22 [█████████████████████████████████-------] 81% [00m:02s] Processing chr. 22 [████████████████████████████████████████] 100% [00m:02s] Context index built`` ``#> [█---------------------------------------] 0% [00m:00s] Saving context index [████████████████████████████████████████] 100% [00m:00s] Context index saved`` ``#> done`` ``#> Building repeated sequence index...`` ``#> [█---------------------------------------] 0% [00m:00s] Reading 22 [████████████████████████████████████████] 100% [00m:01s] Reading 22`` ``#> [████████████████████████████████████████] 100% [00m:01s] Reading 22`` ``#> [████████████████████████████████████████] 100% [00m:01s] Reading 22`` ``#> [████████████████████████████████████████] 100% [00m:01s] Reading 22`` ``#> [████████████████████████████████████████] 100% [00m:01s] Reading 22`` ``#> [████████████████████████████████████████] 100% [00m:01s] Reading 22`` ``#> [████████████████████████████████████████] 100% [00m:01s] Reading 22`` ``#> [████████████████████████████████████████] 100% [00m:01s] RS index built`` ``#> [█---------------------------------------] 0% [00m:00s] Saving RS index [█---------------------------------------] 0% [00m:00s] Saving RS index [█████████████████████-------------------] 51% [00m:02s] Saving RS indexdone`` ``#> [████████████████████████████████████████] 100% [00m:02s] RS index saved`` ``#> [█---------------------------------------] 0% [00m:00s] Loading germline [████████████████████████████████████████] 100% [00m:00s] Germline loaded`` ``#> [█---------------------------------------] 0% [00m:00s] Saving germline [████████████████████████████████████████] 100% [00m:00s] Germline saved`` `` `[`dir.exists`](https://rdrr.io/r/base/files2.html)`(``"Test"``)`` ``#> [1] TRUE`

### The Context Index

The SBS signatures set, for every possible context (i.e., a triplet of
consecutive nucleotides), the probability for an SNV to occur in such a
context (see (Alexandrov et al. 2020)). Hence, the function
[`MutationEngine()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine.md)
builds a context index of the reference sequence to place SNVs on the
forest cells. The complete context index of a sequence whose size is up
to 4Gbps takes 4 times the length of the sequence itself, both in memory
and on the disk, because 4 bytes per nucleotide are required to store
the context position. For instance, the complete context index of a
human genome takes about 12 GB. In order to avoid such a requirement in
memory and disk space,
[`MutationEngine()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine.md)
allows sampling the reference genome contexts and storing only some of
them in the context index. This is achieved by the optional parameter
`context_sampling`, which specifies how many occurrences of the same
context must be identified before adding one of them to the context
index. The larger the number of context sampling, the larger the context
index. On the other side, the lower the number of context sampling, the
lower the number of sites in the reference genome that can be affected
by simulated mutations. The `context_sampling` is set to 100 by default,
but it can be specified during the
[`MutationEngine()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine.md)
call.

`m_engine`` ``<-`` `[`MutationEngine`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine.md)`(``setup_code ``=`` ``"demo"``)`` ``#> Downloading reference genome...`` ``#> Reference genome downloaded`` ``#> Decompressing reference genome...done`` ``#> Downloading signature files...`` ``#> Signature file downloaded`` ``#> Downloading driver mutation file...`` ``#> Driver mutation file downloaded`` ``#> Decompressing driver mutation file...done`` ``#> Downloading passenger CNAs file...`` ``#> Passenger CNAs file downloaded`` ``#> Decompressing passenger CNAs file...done`` ``#> Downloading germline...`` ``#> Germline downloaded`` ``#> Decompressing mutations...`` ``#> done`` ``#> Building context index...`` ``#> [█---------------------------------------] 0% [00m:00s] Processing chr. 22 [█████████████████-----------------------] 40% [00m:00s] Processing chr. 22 [█████████████████████████████████-------] 81% [00m:02s] Processing chr. 22 [████████████████████████████████████████] 100% [00m:02s] Context index built`` ``#> [█---------------------------------------] 0% [00m:00s] Saving context index [████████████████████████████████████████] 100% [00m:00s] Context index saved`` ``#> done`` ``#> Building repeated sequence index...`` ``#> [█---------------------------------------] 0% [00m:00s] Reading 22 [████████████████████████████████████████] 100% [00m:01s] Reading 22`` ``#> [████████████████████████████████████████] 100% [00m:01s] Reading 22`` ``#> [████████████████████████████████████████] 100% [00m:01s] Reading 22`` ``#> [████████████████████████████████████████] 100% [00m:01s] Reading 22`` ``#> [████████████████████████████████████████] 100% [00m:01s] Reading 22`` ``#> [████████████████████████████████████████] 100% [00m:01s] Reading 22`` ``#> [████████████████████████████████████████] 100% [00m:01s] Reading 22`` ``#> [████████████████████████████████████████] 100% [00m:01s] RS index built`` ``#> [█---------------------------------------] 0% [00m:00s] Saving RS index [█---------------------------------------] 0% [00m:00s] Saving RS index [███████████████████---------------------] 45% [00m:02s] Saving RS indexdone`` ``#> [████████████████████████████████████████] 100% [00m:02s] RS index saved`` ``#> [█---------------------------------------] 0% [00m:00s] Loading germline [████████████████████████████████████████] 100% [00m:00s] Germline loaded`` ``#> [█---------------------------------------] 0% [00m:00s] Saving germline [████████████████████████████████████████] 100% [00m:00s] Germline saved`` `` ``` # get the size of the context index when `context_sampling` is 100 ``` ``utils``:::`[`format.object_size`](https://rdrr.io/r/utils/object.size.html)`(`[`file.size`](https://rdrr.io/r/base/file.info.html)`(``"demo/context_index_100.cif"``)``,`` `` standard ``=`` ``"auto"``)`` ``#> [1] "NA bytes"`

Let us rebuild the context index, sampling one context every 50.

`# building a mutation engine by using the "demo" set-up configuration`` ``# and setting context_sampling to 50`` ``m_engine`` ``<-`` `[`MutationEngine`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine.md)`(``setup_code ``=`` ``"demo"``, context_sampling ``=`` ``50``)`` ``#> Building context index...`` ``#> [█---------------------------------------] 0% [00m:00s] Processing chr. 22 [█████████████████-----------------------] 40% [00m:00s] Processing chr. 22 [█████████████████████████████████-------] 81% [00m:02s] Processing chr. 22 [████████████████████████████████████████] 100% [00m:02s] Context index built`` ``#> [█---------------------------------------] 0% [00m:00s] Saving context index [████████████████████████████████████████] 100% [00m:00s] Context index saved`` ``#> done`` ``#> [█---------------------------------------] 0% [00m:00s] Loading RS index [█████████████---------------------------] 31% [00m:01s] Loading RS index [█████████████████████████---------------] 62% [00m:02s] Loading RS index [█████████████████████████████████████---] 92% [00m:03s] Loading RS index [████████████████████████████████████████] 100% [00m:03s] RS index loaded`` ``#> [█---------------------------------------] 0% [00m:00s] Loading germline [████████████████████████████████████████] 100% [00m:00s] Germline loaded`

The above call uses the already downloaded data in the directory
`"demo"` to build the context and produces the file
`"demo/context_index_50.cif"`.

`` # get the size of the context index when `context_sampling` is 50 ``` ``utils``:::`[`format.object_size`](https://rdrr.io/r/utils/object.size.html)`(`[`file.size`](https://rdrr.io/r/base/file.info.html)`(``"demo/context_index_50.cif"``)``, ``"auto"``)`` ``#> [1] "2.7 Mb"`

### The Repeated Sequence Index

Analogously to the SBS signatures, the ID (indel) signatures set, for
some kinds of repeated sequences, the probability for an indel to occur
in them (see (Alexandrov et al. 2020)). Thus, the function
[`MutationEngine()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine.md)
also builds an index for the repeated sequences in the reference. As in
the case of the SBS signatures, the ID signatures can be memory eager.
To mitigate the memory requirement,
[`MutationEngine()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine.md)
accepts as parameters the maximum size of the considered indel (the
parameter `max_index_size`) and the maximum number of reference
positions stored per repetition type (the parameter
`max_repetition_storage`). The larger the two parameters, the larger the
repetition index. As for building time, the larger the max_index_size,
the longer it takes to construct the index. On the contrary,
`max_repetition_storage` does not directly affect the computation time.
By default, the two parameters are set to 50 and 500000, respectively.
However, their values can be specified during the
[`MutationEngine()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine.md)
call.

`# get the size of the repeated index when the default parameter values`` ``# are used`` ``utils``:::`[`format.object_size`](https://rdrr.io/r/utils/object.size.html)`(`[`file.size`](https://rdrr.io/r/base/file.info.html)`(``"demo/rs_index_50_500000.rsif"``)``, ``"auto"``)`` ``#> [1] "212.5 Mb"`

Let us rebuild the context index, sampling one context every 50.

`# building a mutation engine by using the "demo" set-up configuration`` ``# and setting max_repetition_storage to 10M`` ``m_engine`` ``<-`` `[`MutationEngine`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine.md)`(``setup_code ``=`` ``"demo"``,`` `` max_repetition_storage ``=`` ``1e7``)`` ``#> [█---------------------------------------] 0% [00m:00s] Loading context index [████████████████████████████████████████] 100% [00m:00s] Context index loaded`` ``#> Building repeated sequence index...`` ``#> [█---------------------------------------] 0% [00m:00s] Reading 22 [████████████████████████████████████████] 100% [00m:01s] Reading 22`` ``#> [████████████████████████████████████████] 100% [00m:01s] Reading 22`` ``#> [████████████████████████████████████████] 100% [00m:01s] Reading 22`` ``#> [████████████████████████████████████████] 100% [00m:01s] Reading 22`` ``#> [████████████████████████████████████████] 100% [00m:01s] Reading 22`` ``#> [████████████████████████████████████████] 100% [00m:01s] Reading 22`` ``#> [████████████████████████████████████████] 100% [00m:01s] Reading 22`` ``#> [████████████████████████████████████████] 100% [00m:01s] RS index built`` ``#> [█---------------------------------------] 0% [00m:00s] Saving RS index [█---------------------------------------] 0% [00m:00s] Saving RS index [█---------------------------------------] 0% [00m:01s] Saving RS index [█---------------------------------------] 0% [00m:02s] Saving RS index [█---------------------------------------] 0% [00m:03s] Saving RS index [█---------------------------------------] 0% [00m:04s] Saving RS index [█████-----------------------------------] 12% [00m:06s] Saving RS index [████████████████------------------------] 38% [00m:07s] Saving RS index [██████████████████████████--------------] 63% [00m:08s] Saving RS index [███████████████████████████████████-----] 87% [00m:09s] Saving RS indexdone`` ``#> [████████████████████████████████████████] 100% [00m:09s] RS index saved`` ``#> [█---------------------------------------] 0% [00m:00s] Loading germline [████████████████████████████████████████] 100% [00m:00s] Germline loaded`

The above call uses the already downloaded data in the directory
`"demo"` to build the repeated sequence index and produces the file
`"demo/rs_index_50_10000000.rsif"`.

`` # get the size of the context index when `context_sampling` is 50 ``` ``utils``:::`[`format.object_size`](https://rdrr.io/r/utils/object.size.html)`(`[`file.size`](https://rdrr.io/r/base/file.info.html)`(``"demo/rs_index_50_10000000.rsif"``)``,`` `` ``"auto"``)`` ``#> [1] "732.1 Mb"`

Alexandrov, Ludmil B, Jaegil Kim, Nicholas J Haradhvala, et al. 2020.
“The Repertoire of Mutational Signatures in Human Cancer.” *Nature* 578
(7793): 94–101.
