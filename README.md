# rnaseq-staj

clone this repo and then run the commands in `initialize.sh` script manually. Unfortunately, running the mamba installation in batch failes since you cannot source the .bashrc file. So, please run the steps manually.

then you can run snakemake. please change the `annotation` file for your own analysis

for large fastq files, you might need to increase ulimit for sorted BAM files

```
ulimit -S -n 4096
```

## using aegea for aws launch

the code below launches a spot instance at AWS (instance c6a.8xlarge, 32cpu 64GB ram) with spot price of 0.65 and with 250GB root disk (please refer to [this page](https://instances.vantage.sh/) for pricing details.)

```bash
aegea launch  staj --ubuntu-linux-ami --spot-price 0.65 --instance-type c6a.8xlarge --storage /=250GB
```

then you can ssh into the machine simply with `aegea ssh staj`
