# PWAS
This tutorial is created for Proteome-wide association study (PWAS) and further performing a bivariate conditional analysis to interpret PWAS findings in the context of those from Transcriptome-wide association studies (TWAS). 



Introduction
============

This tutorial is created for Proteome-wide association study (PWAS) and
further performing a bivariate conditional analysis to interpret PWAS
findings in the context of those from Transcriptome-wide association
studies (TWAS).

The pipeline of PWAS is mainly from
[TWAS/FUSION](http://gusevlab.org/projects/fusion/) with some
modifications. The bivariate conditional analysis was described in our
recent paper. Please cite both the manuscript for primary TWAS methods
and the recent pQTL and PWAS paper:

> - Gusev, et al. “Integrative approaches for large-scale transcriptome-wide association studies” 2016 Nature Genetics
> - Zhang, Chaterjee, et al. “Large Bi-Ethnic Study of Plasma Proteome Leads to Comprehensive Mapping of cis-pQTL and Models for Proteome-wide Association Studies” **doi**:
https://doi.org/10.1101/2021.03.15.435533



Note that we stored the trained models and required data for two ethnic groups separately (<span style="background-color: #d9d9d9">PWAS\_[EA/AA]</span>). [EA/AA] in file or directory names represents the ancestry where we trained these models (<span style="background-color: #d9d9d9">EA</span> represents for <span style="background-color: #d9d9d9">European Americans</span> and <span style="background-color: #d9d9d9">AA</span> represents for <span style="background-color: #d9d9d9">African Americans</span>). We recommend matching the ancestry to your GWAS summary data. 


Installation
============

Launch R and install required libraries: 
``` {.r language="R"}
install.packages(c('optparse','RColorBrewer','readr','stringr','dplyr'))
install.packages('plink2R-master/plink2R/',repos=NULL)
```

Note: We noticed some people may have issues with the installation of
plink2R. Pleas try the command below in R if the above one does not
work.

    devtools::install_github("carbocation/plink2R/plink2R", ref="carbocation-permit-r361")

Download and unzip the scripts:

    mkdir PWAS
    cd PWAS
	
    wget https://jhupwas.s3.amazonaws.com/scripts.zip
    unzip scripts.zip && rm scripts.zip

    
Download and unzip the required data:


    wget https://jhupwas.s3.amazonaws.com/LDref.zip
    unzip LDref.zip && rm LDref.zip

    ## choose to download models for EA / AA
    wget https://jhupwas.s3.amazonaws.com/PWAS_EA.zip
    unzip PWAS_EA.zip && rm PWAS_EA.zip
    wget https://jhupwas.s3.amazonaws.com/PWAS_AA.zip
    unzip PWAS_AA.zip && rm PWAS_AA.zip
    
    ## GTEx V7 list for performing conditional analysis
    wget https://jhupwas.s3.amazonaws.com/GTEx_V7_list.zip
    unzip GTEx_V7_list.zip && rm GTEx_V7_list.zip
    
    
We also provides some example output results:

    wget https://jhupwas.s3.amazonaws.com/Results.zip
    unzip Results.zip && rm Results.zip
    



PWAS – analysis and output
==========================


Required data to perform PWAS
-----------------------------

LD reference data (LDref in our data directory), weights
(e.g. Plasma\_Protein\_weights\_[EA/AA] in our data directory), and
weight-loading files (e.g. Plasma\_Protein\_[EA/AA]\_hg19.pos in our data
directory).



Input 1: GWAS summary statistics
--------------------------------

Same format of input of the summary statistics: A white-space separated
table with a header row containing <span style="background-color: #d9d9d9">SNP (rsid)</span>, <span style="background-color: #d9d9d9">A1 (effect allele)</span>, <span style="background-color: #d9d9d9">A2 (other allele)</span>, and <span style="background-color: #d9d9d9">Z (z-score)</span>.

Input 2: Weigths
----------------

Imputation models were available for 1,318 and 1,368 significant cis-heritable plasma proteins which have significant non-zero cis-heritability estimated by GCTA (p<0.01) in EA and AA, respectively. Currently, only enet (elastic-net) and top1 (best single SNP) models are available. We recommend using enet models because they have better prediction accuracy. 

Inside <span style="background-color: #d9d9d9">PWAS\_[EA/AA]</span>, weights are stored in <span style="background-color: #d9d9d9">Plasma\_Protein\_weights\_[EA/AA]</span>. Weights are loaded from <span style="background-color: #d9d9d9">Plasma\_Protein\_[EA/AA]\_[hg19/hg38].pos</span>. [hg19/hg38] is the genome build which is used to locate those proteins. We recommend using the <span style="background-color: #d9d9d9">[hg19]</span> to match TWAS models based on GTEx V7 from FUSION.


Input 3: LD reference data 
--------------------------

LD reference data for European (EUR) and African (AFR) ancestry from 1000 Genomes Project (1000G) are stored in <span style="background-color: #d9d9d9">LDref</span>.

Example for performing association test
---------------------------------------

The association test was performed by chromosome. Below is an example
using data on chromosome 22:
    
    CHR=22
    #until [ $CHR -lt 1 ]
    #do
    Rscript ./scripts/PWAS.assoc_test.R \
    --sumstats ./sumstats.txt \
    --weights ./PWAS_EA/Plasma_Protein_EA_hg19.pos \
    --weights_dir ./PWAS_EA/Plasma_Protein_weights_EA/ \
    --ref_ld_chr ./LDref/EUR/chr \
    --force_model enet \
    --chr ${CHR} \
    --out ./Results/Example_PWAS_by_chr/chr${CHR}.out
    #let CHR-=1
    #done

Output: Gene-disease association
--------------------------------

The output file has the same format as
[TWAS/FUSION](http://gusevlab.org/projects/fusion/#typical-analysis-and-output).
(The only differences are eQTL <span>&#10230;</span> pQTL; TWAS <span>&#10230;</span> PWAS)

Bivariate conditional analysis for TWAS and PWAS
======================================

We provide a framework for interpreting results from PWAS in conjunction
with analogous results from TWAS. As demonstrated in our paper, genetic correlation of gene-expression in tissues and protein level in plasma due to *cis* genetic regulation is moderate. As a result, it is possible to explore whether the PWAS signals could be explained by cis-genetic regulation of the expression of nearby (1Mb region around) genes and vice versa using a bivariate conditional analysis framework.

Required data to perform conditional analysis
---------------------------------------------

TWAS output for only **[significant](http://gusevlab.org/projects/fusion/weights/GTEX7.txt)** cis-heritable gene expressions performed using **enet** models, PWAS output from our pipeline performed using **enet** models, and pre-imputed *cis*-regulated plasma proteins and gene expressions for 1000G reference individuals which are stored in <span style="background-color:#d9d9d9">1000G\_imputed\_[EA/AA]</span> inside <span style="background-color: #d9d9d9">PWAS\_[EA/AA]</span>).

Input 1: PWAS and TWAS output tables
------------------------------------

The output tables in the FUSION format are required. Here we also require **enet** models for both PWAS and TWAS (trained with the option --force\_model enet), because the imputed plasma protein and gene expression levels for reference individuals were pre-computed using enet models (stored in <span style="background-color:#d9d9d9">1000G\_imputed\_[EA/AA]</span> inside <span style="background-color: #d9d9d9">PWAS\_[EA/AA]</span>). 

Output tables for all chromosomes from TWAS/PWAS are needed to be merged together. Example R code:

``` {.r language="R"}
library(dplyr)
library(readr)
results <- tibble()
for (chr in 1:22) {
    results <- rbind(results, read_tsv(paste0("./Results/Example_PWAS_by_chr/chr", chr, ".out")))
    if(chr==6){
    	results <- rbind(results, read_tsv(paste0("./Results/Example_PWAS_by_chr/chr6.out.MHC")))
    }
}
write_tsv(results, "./Results/PWAS.out")
```

TWAS output tables for different tissues should be stored separately with the name of <span style="background-color: #d9d9d9">*tissue*.out</span> in a same folder. For example, the TWAS in whole blood should be stored with name of <span style="background-color: #d9d9d9">Whole\_Blood.out</span>.

Input 2: list of TWAS tissues to be analyzed
--------------------------------------------

The full list of TWAS tissues to be analyzed for conditional analysis
should be stored in a text file. One tissue per line. An example is the
full list of all GTEx V7 tissues, <span style="background-color: #d9d9d9">GTEx\_V7\_tissue\_list.txt</span>.

Input 3: imputed plasma protein and gene expression reference data
--------------------------------------------------------

*Cis*-genetic regulated plasma proteins and gene expressions were
pre-imputed on publicly available 1000G individuals. This data is stored in <span style="background-color: #d9d9d9">1000G\_imputed\_[EA/AA]</span> inside <span style="background-color: #d9d9d9">PWAS\_[EA/AA]</span>.

Example for performing conditional analysis of TWAS and PWAS
-------------------------------------------------------------

Here we give an example for performing the conditional analysis of TWAS and
PWAS. In this example, we stored the PWAS output table in <span style="background-color: #d9d9d9">./Results/PWAS.out</span>, and
all-tissue TWAS output tables inside the folder of <span style="background-color: #d9d9d9">./Results/TWAS.out/</span>.
    
    Rscript ./scripts/PWAS.conditional.R \
    --PWAS ./Results/PWAS.out \
    --TWAS ./Results/TWAS.out/ \
    --tissue_list ./GTEx_V7_tissue_list.txt \
    --tissue_n_gene ./GTEx_V7_n_gene.rds \
    --imputed_P ./PWAS_EA/1000G_imputed_EA/1000G_imputed_Plasma_Protein.txt \
    --imputed_T ./PWAS_EA/1000G_imputed_EA/1000G_imputed_FUSION/ \
    --out ./Results/ConditionalAnalysis/

Output: bivariate conditional analysis of PWAS and TWAS in all tissues
-----------------------------------------------------------------------

This code will first generate a <span style="background-color: #d9d9d9">tissue.RDat</span> file for each tissue. These RDat files contain
the flowing variables: (below we use T to represent gene expressions,
and P to represent plasma proteins)

-   <span style="background-color: #d9d9d9">dat.sentinel.pwas</span>: The regional sentinel PWAS genes (+/- 500Kb).

-   <span style="background-color: #d9d9d9">PcT.z</span>: z-score of P conditional on T.

-   <span style="background-color: #d9d9d9">PcT.p</span>: p-value of P conditional on T.

-   <span style="background-color: #d9d9d9">TcP.z</span>: z-score of T conditional on P.

-   <span style="background-color: #d9d9d9">TcP.p</span>: p-value of T conditional on P.

-   <span style="background-color: #d9d9d9">twas.hit</span>: the most significant nearby TWAS gene (+/- 500Kb) for each
    regional sentinel PWAS gene.

-   <span style="background-color: #d9d9d9">twas.p</span>: p-value of the most significant nearby TWAS gene regional
    sentinel PWAS gene.

-   <span style="background-color: #d9d9d9">dist</span>: distance (bp) between each sentinel PWAS gene and its most
    significant nearby TWAS gene.

-   <span style="background-color: #d9d9d9">corr</span>: cis-regulated genetic correlation between each sentinel PWAS gene and its most significant nearby TWAS gene.
    

It will then generate summary tables for tissue-specific conditional analysis and an all-tissue conditional analysis. The columns are:

-   <span style="background-color: #d9d9d9"> PWAS\_hit </span>: PWAS significant gene

-   <span style="background-color: #d9d9d9"> PWAS\_p </span>: p-value o PWAS significant gene

-   <span style="background-color: #d9d9d9"> TWAS\_hit </span>: the most significant nearby TWAS gene (+/- 500Kb) for each regional sentinel PWAS gene.

-   <span style="background-color: #d9d9d9"> TWAS\_p </span>: p-value of the most significant nearby TWAS gene regional sentinel PWAS gene.

-   <span style="background-color: #d9d9d9"> Dist\_of\_hits </span>: distance (bp) between each sentinel PWAS gene and its most significant nearby TWAS gene.

-   <span style="background-color: #d9d9d9"> Corr\_of\_hits </span>: cis-regulated genetic correlation between each sentinel PWAS gene and its most significant nearby TWAS gene.

-   <span style="background-color: #d9d9d9"> PcT\_p </span>: p-value of P conditional on T.

-   <span style="background-color: #d9d9d9"> TcP\_p </span>: p-value of T conditional on P.

-   <span style="background-color: #d9d9d9"> min\_TWAS\_hit </span>: the most significant nearby TWAS gene (+/- 500Kb) for each regional sentinel PWAS gene in all-tissue analysis.

-   <span style="background-color: #d9d9d9"> min\_TWAS\_Tissue </span>: the tissue with the most significant nearby TWAS gene (+/- 500Kb) for each regional sentinel PWAS gene in all-tissue analysis.

-   <span style="background-color: #d9d9d9"> min\_TWAS\_p </span>: the p-value of the most significant nearby TWAS gene (+/- 500Kb) for each regional sentinel PWAS gene in all-tissue analysis.

-   <span style="background-color: #d9d9d9"> N\_significant\_tissues\_in\_TWAS </span>: the number of significant TWAS tissues in all-tissue analysis.

-   <span style="background-color: #d9d9d9"> all\_significant\_tissues\_in\_TWAS </span>: all significant TWAS tissues in all-tissue analysis.

