# nf-core/methylseq

## v1.7dev

### Pipeline Updates

* Added Picard CollectInsertSizeMetrics and Picard CollectGcBiasMetrics
* Improved qulimap and preseq by adding `samtools sort` and `samtools index` step in the Bismark aligner
* Added BISCUIT aligner as an optional aligner, with all relative steps (alignment, mark duplicates with [samblaster](https://github.com/GregoryFaust/samblaster), methylation extraction, QC for biscuit, and optional [Epi-read](https://huishenlab.github.io/biscuit/epiread_format/) file creation with SNP information).

### New software

* samblaster `0.1.26`
* bedtools `2.30.0`
* biscuit `0.3.16.20200420`
* bcftools`1.10`
* parallel `20201122`
* gawk `5.1.0`

## [v1.6](https://github.com/nf-core/methylseq/releases/tag/1.6) - 2021-03-26

**:warning: Breaking change!**

In line with a standardisation change across all of nf-core, we have changed the main parameter name for supplying files to the pipeline.
In this release, please use `--input` instead of `--reads`.
The parameter still works in the same way as before.

### Pipeline Updates

* Increased resources for `fastqc` process ([#143](https://github.com/nf-core/methylseq/issues/143))
* Raised Nextflow version requirement to `20.07.1`
* Updated template to nf-core/tools 1.13.3
* Renamed `--reads` to `--input`
* Added new `--maxins` and `--minins` parameters to pass on to Bismark
* New `--em_seq` preset
  * Sets `bismark_maxins = 1000`, `clip_r1 = 8`, `clip_r2 = 8`, `three_prime_clip_r1 = 8`, `three_prime_clip_r2 = 8`
* New `--publish_dir_mode` parameter to customise results folder behaviour
* Fix bug on AWS for `bismark_hisat` known splice file ([#177](https://github.com/nf-core/methylseq/issues/177))
* Moved parameter documentation into new `nextflow_schema.json` file
  * This improves web documentation and enables `nf-core launch` functionality. See <https://nf-co.re/launch?pipeline=methylseq>
* Added a `-profile test_full` config for running the pipeline with a full-size test dataset
  * See [the config file](https://github.com/nf-core/methylseq/blob/dev/conf/test_full.config) for details
  * This will be used for automated release tests on AWS, results browsable on the website

### Software updates

* Python base `3.7.3` > `3.8.8`
* markdown `3.1.1` > `3.3.4`
* pymdown-extensions `6.0` > `8.1.1`
* pygments `2.6.1` > `2.8.1`
* pigz `2.3.4` > `2.6`
* samtools `1.9` > `1.11`
* TrimGalore! `0.6.5` > `0.6.6`
* Bowtie2 `2.3.5` > `2.4.2`
* Hisat2 `2.2.0` > `2.2.1`
* Bismark `0.22.3` > `0.23.0`
* Picard `2.22.2` > `2.25.1`
* MethylDackel `0.5.0` > `0.5.2`
* MultiQC `1.8` > `1.10`

## [v1.5](https://github.com/nf-core/methylseq/releases/tag/1.5) - 2020-04-09

### New features

* Added multicore support for `TrimGalore!`
* Improved the multicore support for Bismark methXtract for more parallelisation ([#121](https://github.com/nf-core/methylseq/issues/121))
* Added `--cytosine_report` option to tell Bismark to give reports for all cytosines in the genome.
* Added options `--bismark_align_cpu_per_multicore` and `--bismark_align_cpu_per_multicore` to customise how Bismark align `--multicore` is decided ([#124](https://github.com/nf-core/methylseq/issues/124))
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

### Software updates

* _New_: pigz `2.3.4`
* Python base `2.7` > `3.7.3`
* FastQC `0.11.8` > `0.11.9`
* TrimGalore! `0.6.4` > `0.6.5`
* HiSAT2 `2.1.0` > `2.2.0`
* Bismark `0.22.2` > `0.22.3`
* Qualimap `2.2.2c` > `2.2.2d`
* Picard `2.21.3` > `2.22.2`
* MethylDackel `0.4.0` > `0.5.0`
* MultiQC `1.7` > `1.8`

### Pipeline Updates

* Fixed bug where the iGenomes config was loaded after the references were set ([#121](https://github.com/nf-core/methylseq/issues/121))
* Removed awsbatch config profile because it is now served by [nf-core/configs](https://github.com/nf-core/configs)
* Tidied up the summary log messages when starting the pipeline
  * Fewer messages saying what you're _not_ doing, sanitised the order of some logs and removed a few things
* Slightly refactored the code for trimming parameters
* Updated template to tools 1.9

### Bug fixes

* Fixed error where MethylDackel would consume the Nextflow channels and not work with more than one sample [#140](https://github.com/nf-core/methylseq/issues/140)

## [v1.4](https://github.com/nf-core/methylseq/releases/tag/1.4) - 2019-11-19

### New features

* Changed all parameter names to `snake_case`
* Added `--local_alignment` option to run Bismark with the `--local` flag to allow soft-clipping of reads.
* Added support for bismark's [SLAM-seq mode](https://github.com/FelixKrueger/Bismark/blob/master/CHANGELOG.md#slam-seq-mode)
* Added support for running bismark with HISAT2 as an aligner option [#85](https://github.com/nf-core/methylseq/issues/85)
* Added support for centralized configuration profiles [nf-core/configs](https://github.com/nf-core/configs)
* Add `--meth_cutoff` parameter to change default for `bismark_methylation_extractor`
  * eg. use `--meth_cutoff 5` on the command line or `params.meth_cutoff = 5` to require 5 overlapping reads to call a methylation site.
* Added `--methyl_kit` option to run MethylDackel with the `--methylKit` flag, producing output suitable for the methylKit R package.

### Software updates

* _new dependency_: hisat2 `2.1.0`
* _new dependency_: r-markdown `1.1`
* TrimGalore! `0.5.0` > `0.6.4`
* Bismark `0.20.0` > `0.22.2`
* Bowtie2 `2.3.4.3` > `2.3.5`
* Picard `2.18.21` > `2.21.3`
* Qualimap `2.2.2b` > `2.2.2c`
* MethylDackel `0.3.0` > `0.4.0`

### Pipeline updates

* Keep memory in GB for samtools, to avoid problems with unit conversion ([#99](https://github.com/nf-core/methylseq/issues/99))
* Changed `params.container` for `process.container`
* Synchronised with version 1.7 of the nf-core/tools template

### Bug fixes

* Fixed a bug that caused conda dependencies to be resolved very slowly
* Allowed some spare memory in the samtools sort steps, avoiding crashes for some users ([#81](https://github.com/nf-core/methylseq/issues/81))

## [v1.3](https://github.com/nf-core/methylseq/releases/tag/1.3) - 2019-02-01

### New features

* Added [preseq](http://smithlabresearch.org/software/preseq/) analysis to calculate sample complexity.
  * This new step can help decide sufficient sequencing depth has been reached.

### Bug fixes

* Fixed new bug that meant pipeline only worked with one sample at a time [#66](https://github.com/nf-core/methylseq/issues/66)
  * Introduced in previous release. TrimGalore onwards would only process one sample.

## [v1.2](https://github.com/nf-core/methylseq/releases/tag/1.2) - 2019-01-02

### New features

* Trim 9bp from both ends of both reads for PBAT mode.
* Save `where_are_my_files.txt` to results directory to inform the user about missing intermediate files [#42](https://github.com/nf-core/methylseq/issues/42)

### Software updates

* Fastqc `0.11.7` > `0.11.8`
* Bowtie2 `2.3.4.2` > `2.3.4.3`
* Bismark `0.19.1` > `0.20.0`
* Qualimap `2.2.2a` > `2.2.2b`
* Picard `2.18.11` > `2.18.21`
* MultiQC `1.6` > `1.7`

### Bug fixes

* Fixed error when running the pipeline with `--unmapped`
  * Previously, could result in error `Error ~ No such variable: bismark_unmapped`
* Fixed error where single-sample reports could mix up log files [#48](https://github.com/nf-core/methylseq/issues/48)
* Fixed bug in MultiQC process that skipped results from some tools
* Supply available memory as argument to Picard MarkDuplicates

## [v1.1](https://github.com/nf-core/methylseq/releases/tag/1.1) - 2018-08-09

* Tests simplified - now work by simply using the `test` config profile
  * eg: `nextflow run nf-core/methylseq -profile test,docker`
  * Removed previous `run_test.sh` script and data
* New `Singularity` build script for direct compatibility with [singularity-hub](https://singularity-hub.org/)
* Minor improvements to the docs
* A number of boilerplate nf-core code updates
* Updated `process$name` nextflow syntax to avoid warnings in new versions of nextflow
* Updated software tools
  * `trim-galore` `v0.4.5` update to `0.5.0`
  * `samtools` `v1.8` update to `1.9`
  * `bowtie2` `v2.3.4.1` update to `2.3.4.2`
  * `multiqc` `v1.5` update to `1.6`
  * `picard` `v2.18.2` update to `2.18.11`
  * `bwameth` `v0.2.0` update to `0.2.2`

## [v1.0](https://github.com/nf-core/methylseq/releases/tag/1.0) - 2018-04-17

Version 1.0 marks the first release of this pipeline under the nf-core flag. It also marks a significant step up in the maturity of the workflow, with everything now in a single script and both aligner workflows fully supported.

* Renamed and moved [SciLifeLab/NGI-MethylSeq](https://github.com/SciLifeLab/NGI-MethylSeq/) to [nf-core/methylseq](https://github.com/nf-core/methylseq/)
* Merged bwa-meth and bismark pipeline scripts, now chosen with `--aligner` flag
* Refactored multi-core parameters for Bismark alignment and methylation extraction
* Rewrote most of the documentation
* Changed the Docker container to use Bioconda installations

---

Previous to these releases, this pipeline was called [SciLifeLab/NGI-MethylSeq](https://github.com/SciLifeLab/NGI-MethylSeq):

## v0.4dev

* Fixed MultiQC channel bug
* Integrated config for QBiC Tuebingen
* Numerous small container bugfixes
* Refactored how the config is loaded
* Fix for resource limit function, improved resource request defaults
* Fix for iGenomes base path in configs

## [v0.3.1](https://github.com/SciLifeLab/NGI-MethylSeq/releases/tag/0.3.1) - 2017-09-05

* Include base profile name and documentation about Singularity.
* Testing automated docker hub image tagging for releases.

## [v0.3](https://github.com/SciLifeLab/NGI-MethylSeq/releases/tag/0.3) - 2017-09-01

* Fix `--rrbs` mode ([#24](https://github.com/SciLifeLab/NGI-MethylSeq/issues/24))
* Fixed fairly major bug where only a single sample would run past alignment
* Merged test scripts and rewrote to use command line flags / new travis script.
* Refactored software version collection code to be more resilient and cleaner / easier to maintain.
* Dropped support for environment modules and added support for use of Singularity on UPPMAX

## [v0.2](https://github.com/SciLifeLab/NGI-MethylSeq/releases/tag/0.2) - 2017-07-17

First (semi-) stable release of the new NGI-MethylSeq pipeline, as we head towards deployment in production.
