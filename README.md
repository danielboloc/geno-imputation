# Introduction

# Pre-requisites

## Java
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

Follow the [instructions][GENEAPI] on how to connect to the web interface in localhost:8080 to see your imputation jobs after initiating them through Nextflow.

You will have to login into the web server with **admin** and default password **admin1978** mentioned in the documentation. You can create other profiles or use this one. In your user, go to "Profile"->"API Access" and press "Create API Token". Save this token in an environment variable or in a secure place. It will be needed in the run given to `params.token`.

In the "Admin Panel" go to "Applications" and install "1000 Genomes Phase 3" in version 2.0.0 (take a coffee since it will take some time to get setup).

# Fasta 

The fasta file for human is located [here][FASTA].

# Run

## Command line options

## Mock example

First export token or save it in ENV:
```bash
export TOKEN=*********************************
```
If you have every pre-requisite installed the following command should work when adapted to user paths/data:
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
correct imputation server path to where the localhost:8080 directs.

[GENEAPI]: https://github.com/genepi/imputationserver-docker
[EAGLE]: https://alkesgroup.broadinstitute.org/Eagle/#x1-320005.3.2
[FASTA]: https://www.ncbi.nlm.nih.gov/genome?term=human&cmd=DetailsSearch