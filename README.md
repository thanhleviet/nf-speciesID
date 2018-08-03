A small nextflow script for screening contaminants from Illumina short reads and estimate bacterial genome size

It uses centrifuge for reads classification; some options (kmc, bbtools, mash) for estimation of genome size.

## Dependencies:

All the tools can be installed via [bioconda](http://bioconda.github.io/)

    - centrifuge
    - mash
    - krona

- **Centrifuge database and Krona taxonomy**: https://genomics.sschmeier.com/downloads.html#downloads

- **Mash RefSketch database**: https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh

## How to run:

- Edit location to the centrifuge database and Mash
- Edit pattern of reads if needed

`nextflow run speciation.nf --input /path/to/folder/of/reads --output /path/to/folder/of/results`

