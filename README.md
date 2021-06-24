# ![nf-core/methylseq](docs/images/nf-core-methylseq_logo.png)

[![GitHub Actions CI Status](https://github.com/nf-core/methylseq/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/methylseq/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/nf-core/methylseq/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/methylseq/actions?query=workflow%3A%22nf-core+linting%22)
[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/methylseq/results)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.1343417-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.1343417)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.03.0--edge-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23methylseq-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/methylseq)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)
[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/methylseq** is a bioinformatics analysis pipeline used for Methylation (Bisulfite) sequencing data. It pre-processes raw data from FastQ inputs, aligns the reads and performs extensive quality-control on the results.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker / Singularity containers making installation trivial and results highly reproducible.

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/methylseq/results).

## Pipeline Summary

The pipeline allows you to choose between running either [Bismark](https://github.com/FelixKrueger/Bismark) or [bwa-meth](https://github.com/brentp/bwa-meth) / [MethylDackel](https://github.com/dpryan79/methyldackel).
Choose between workflows by using `--aligner bismark` (default, uses bowtie2 for alignment), `--aligner bismark_hisat` or `--aligner bwameth`.

| Step                                         | Bismark workflow | bwa-meth workflow     |
|----------------------------------------------|------------------|-----------------------|
| Generate Reference Genome Index _(optional)_ | Bismark          | bwa-meth              |
| Raw data QC                                  | FastQC           | FastQC                |
| Adapter sequence trimming                    | Trim Galore!     | Trim Galore!          |
| Align Reads                                  | Bismark          | bwa-meth              |
| Deduplicate Alignments                       | Bismark          | Picard MarkDuplicates |
| Extract methylation calls                    | Bismark          | MethylDackel          |
| Sample report                                | Bismark          | -                     |
| Summary Report                               | Bismark          | -                     |
| Alignment QC                                 | Qualimap         | Qualimap              |
| Sample complexity                            | Preseq           | Preseq                |
| Project Report                               | MultiQC          | MultiQC               |

## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation) (`>=20.04.0`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```bash
    nextflow run nf-core/methylseq -profile test,<docker/singularity/podman/shifter/charliecloud/conda/institute>
    ```

    * Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
    * If you are using `singularity` then the pipeline will auto-detect this and attempt to download the Singularity images directly as opposed to performing a conversion from Docker images. If you are persistently observing issues downloading Singularity images directly due to timeout or network issues then please use the `--singularity_pull_docker_container` parameter to pull and convert the Docker image instead. It is also highly recommended to use the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) settings to store the images in a central location for future pipeline runs.
    * If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

    ```bash
    nextflow run nf-core/methylseq -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> --input samplesheet.csv --genome GRCh37
    ```

See [usage docs](https://nf-co.re/methylseq/usage) for all of the available options when running the pipeline.

## Documentation

The nf-core/methylseq pipeline comes with documentation about the pipeline: [usage](https://nf-co.re/methylseq/usage) and [output](https://nf-co.re/methylseq/output).

## Credits

These scripts were originally written for use at the [National Genomics Infrastructure](https://portal.scilifelab.se/genomics/) at [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden.

* Main author:
  * Phil Ewels ([@ewels](https://github.com/ewels/))
* Contributors:
  * Rickard Hammarén ([@Hammarn](https://github.com/Hammarn/))
  * Alexander Peltzer ([@apeltzer](https://github.com/apeltzer/))
  * Patrick Hüther ([@phue](https://github.com/phue/))

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#methylseq` channel](https://nfcore.slack.com/channels/methylseq) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use  nf-core/methylseq for your analysis, please cite it using the following doi: [10.5281/zenodo.1343417](https://doi.org/10.5281/zenodo.1343417)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
