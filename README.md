# PanGenomeGrowth
The pangenome growth calculator

First running the `./dist/PanGenomeGrowth.jar`

```
java -jar ./dist/PanGenomeGrowth.jar -f vcfFile.gz -n sampleNum -idx index -chr chrlist -core 0.95 -common 0.05 -superLT 3000000 -r region.bed -N T -h
*********The parameters*******************
-f              [required] The VCF File in gzipped format
-n              [required] The number of samples
-idx            [optional] The column index (start from 0) of the first sample in vcf, default: 9
-chr            [optional] The file of list of chrs to be kept in results. All chrs in one line and seperating by tab, default: all chromosomes
-core           [optional] The percentage of all haplotyes represents the core, default: 0.95
-common         [optional] The percentage of all haplotypes represents the common, default: 0.05
-superLT                [optional] The length threshold of super larger segment that need to remove, default: 3000000
-r              [optional] The regions stored in bed file format that need to remove, default: no
-N              [optional] True(T) or False(F) remove regions with N or not, default: False
-h               print this help message
*******************************************
```

Then using the `Rscript --vanilla depthProcessor.R` to get the plot

```
Rscript --vanilla depthProcessor.R theVCFFile[required] sampleNum[required] sampleOrderFile[optional]
```
