# vcf_filtering
workflow and documentation of filtering of VCFs

The workflow moves in stages

First, with a VCF in an unfiltered folder in the Josephs Lab shared HPCC space

```bash
bcftools query -l <VCF>
```

Then, to evaluate which individuals made it into the VCF, I use ==checkoutVCF.R==

| Item      | Description |
| ----------- | ----------- |
| Individuals      | contains lists of specific individuals in VCF       |
| *.R   | Rscripts that are called within the VCF filtering script (aside from checkout file)         |
