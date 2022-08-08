# Introduction

Imputation refers to the statistical inference of unobserved genotypes. It is achieved by using known haplotypes in a population, for instance from the HapMap or the 1000 Genomes Project in humans, thereby allowing to test for association between a trait of interest (e.g. a disease) and experimentally untyped genetic variants, but whose genotypes have been statistically inferred ("imputed"). Genotype imputation is usually performed on SNPs. More information can be found [here][imputation].

The benefit of performing imputation on GWAS data is increasing SNP density, allowing for more accurate computations.

This repository dows not accept data from GWAS in its current state. It accepts only .vcf[.gz] files.

It uses a local instance of the [Michigan Imputation Server][MIS] (MIS) to perform the imputation.

Note: Phasing (done with Eagle) can also be done through the MIS API, but in the current implementation the process in done in the workflow. Might change in the near future.

# Pre-requisites

## Java

More info can he found [here][nxt-docs]

## Nextflow

You can install it with this:
```
curl -s https://get.nextflow.io | bash
chmod +x nextflow
mv nextflow /bin/ (or your PATH)
```

## Bcftools

You can install it with this:
```bash
wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2
tar -xvf bcftools-1.9.tar.bz2
cd bcftools-1.9
./configure --prefix=/usr/local
make
make install
cd ..
rm -rf bcftools-1.9*
```

## Eagle

You can install it with this:
```bash
wget https://data.broadinstitute.org/alkesgroup/Eagle/downloads/Eagle_v2.4.1.tar.gz
gunzip Eagle_v2.4.1.tar.gz
tar xvf Eagle_v2.4.1.tar
mv Eagle_v2.4.1/eagle /usr/local/bin/
rm -r Eagle_v2.4.1
rm -r Eagle_v2.4.1.tar
```

Eagle reference panels can be found [here][EAGLE]. 

## 7z

```bash
sudo apt install p7zip-full p7zip-rar
```

## Michigan Imputation Server

Get the image:
```bash
docker pull genepi/imputationserver:v1.4.1
```

Create in your prefered directory a persistent data folder:
```bash
mkdir -p MICHIGAN_IMPUTATION_SERVER/data
```
Start it with the default Hapmap2 Reference Panel. To install aditional references, check this [documentation][GENEAPI].

To start using:
```
docker run -d -p 8080:80 -v MICHIGAN_IMPUTATION_SERVER/data/:/data/ genepi/imputationserver:v1.4.1
```

Follow the [instructions][GENEAPI] on how to connect to the web interface in `localhost:8080` to see your imputation jobs after initiating them through Nextflow.

You will have to login into the web server with **admin** and default password **admin1978** mentioned in the MIS instructions. You can create other profiles or use this one. In your user, go to "Profile"->"API Access" and press "Create API Token". Save this token in an environment variable or in a secure place. It will be needed in the run given to `params.token`.

In the "Admin Panel" go to "Applications" and install "1000 Genomes Phase 3" in version 2.0.0 (take a coffee since it will take some time to get setup).

## Fasta 

The fasta file for human is located [here][FASTA].

# Build Docker image

All dependencies (except MIS) can be build with the `Dockerfile`. If using **buildkit**, you can use with `DOCKER_BUILDKIT=1`. The name of the image can be the same as the one in the following example or another one.

```console
DOCKER_BUILDKIT=1 docker build -t genoimputation:v1.0.0 .
```

To use it in your commands:
```console
nextflow run geno-imputation/main.nf -with-docker genoimputation:v1.0.0
```

Or define it in `nextflow.config` either in a new `profile` or in the `docker` scope:

```groovy
process.container = "genoimputation:v1.0.0"
docker {
    enabled = true
}
```

# Run

## Command line options

The parameter description is as follows:
```
USAGE

Mandatory arguments:
-------------------
--vcf                 FILE    Input .vcf file (it need to have .csi index from bcftools)
--ref_assembly        FILE    FASTA reference file                                      
--outdir              PATH    where output should be stored                             
--eagle_ref_panel     PATH    where chromosome reference files are located              
--token               STR     Michigan Imputation Server (MIS) local token              
--imputation_job_dir  PATH    where MIS 'data' folder is located on disk                

Optional arguments:
------------------
--max_memory          STR     specifies memory to allocate                              
--max_cpus            INT     specifies number of logical cpus required                 
--max_time            FILE    how long a process is allowed to run                      
--monochrome_logs     FILE    whether to display logs using just one color (white)      

```

## Mock example

First export token or save it in ENV:
```bash
export TOKEN=*********************************
```
If you have every pre-requisite installed the following command should work when adapted to your user paths/data.

Note: although in this example all mandatory arguments are given, if running on the same machine, they can be set in `nextflow.config` and only given the `--vcf`. Also, new config profiles can be created and parameters changed for other connected machines. The profile can be activated with `-profile`.
```
nextflow run geno-imputation/main.nf \
    --vcf /HDD/example.vcf.gz \
    --ref_assembly /HDD/ref/Homo_sapiens_assembly38.fasta \
    --eagle_ref_panel /HDD/eagle/reference/panel/ \
    --token $TOKEN \
    --imputation_job_dir MICHIGAN_IMPUTATION_SERVER/data/jobs \
    --outdir /HDD/imputation_analysis/results
```

Make sure `--imputation_job_dir` is corresponding to a real path. This is not
refering to where would you like to save the output, this parameter needs the
correct imputation server path to where the `localhost:8080` directs.

[GENEAPI]: https://github.com/genepi/imputationserver-docker
[EAGLE]: https://alkesgroup.broadinstitute.org/Eagle/#x1-320005.3.2
[FASTA]: https://www.ncbi.nlm.nih.gov/genome?term=human&cmd=DetailsSearch
[imputation]: https://en.wikipedia.org/wiki/Imputation_(genetics)
[nxt-docs]: https://www.nextflow.io/docs/latest/getstarted.html
[MIS]: https://github.com/genepi/imputationserver-docker