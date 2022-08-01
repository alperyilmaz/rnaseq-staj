# rnaseq-staj

clone this repo and then run the commands in `initialize.sh` script manually. Unfortunately, running the mamba installation in batch failes since you cannot source the .bashrc file. So, please run the steps manually.

then you can run snakemake. please change the `annotation` file for your own analysis

## using aegea for aws launch

the code below launches a spot instance at AWS (instance c6a.8xlarge, 32cpu 64GB ram) with spot price of 0.65 and with 250GB root disk

```bash
aegea launch  staj --ubuntu-linux-ami --spot-price 0.65 --instance-type c6a.8xlarge --storage /=250GB
```

then you can ssh into the machine simply with `aegea ssh staj`
