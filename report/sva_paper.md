---
title: "Assessing the best practice for surrogate variable analysis in epigenome-wide association studies"
author: "Thomas Battram1, Matthew Suderman1, Paul Yousefi1, James Staley1"
output:
  fig_caption: no
---

1MRC Integrative Epidemiology Unit, University of Bristol

















## Introduction

Epigenome-wide association studies (EWAS) examine the association between differences in DNA methylation and a trait of interest. High-throughput experiments are required to measure DNA methylation across thousands or millions of CpG sites in large cohorts. These experiments generate unwanted variation in DNA methylation that is related not to the trait of interest being measured, but instead some batch variable. In order to counteract these batch effects and potentially other unwanted variation, surrogate variables (SVs) can be used as covariates in an EWAS. 

Surrogate variable analysis (SVA) was introduced by Leek et al. in 2007 as a method to deal with unwanted variation in gene expression measurement with expression arrays. Since then it has been applied to various different array-based 'omic' measurements including the popular Illumina blah arrays. 

__Add something about why computational efficiency is important__ -e.g. sva is often the most time consuming step + it must be performed for each phenotype of interest

The computational efficiency has been improved upon by Chen et al. 2017. According to Chen et al., SmartSVA gives the same SVs as the original SVA, but generates these SVs ten to fifty times faster. A further simple method to improve computational efficiency whilst not compromising the SVs produced could be performing SVs on a representative subset of the array.

One major problem that researchers experience with SVA is the number of SVs suggested by the SVA or SmartSVA package. Often this can be so high that including them in a linear model would result in major over-fitting. It may also be very low, suggesting few SVs are needed to account for unwanted variation, when in reality more are needed. This has lead to many researchers simply adjusting for an arbitrary number of SVs such as ten or twenty. 

__Add something about the fact there number of SVs as well as subsetting the array may differ massively between datasets__

In this paper we aimed to test the length of time taken to generate SVs using various different parameters, compare the SVs generated and find out how many SVs are needed to account for variation in typical covariates.

## Methods

### Data

#### Simulated data
Continuous data was simulated using the rnorm() function and binary data was simulated using sample() in R. 

#### Avon Longitudinal Study of Parents and Children (ALSPAC)
The real data in the study mainly came from one cohort: the Avon Longitudinal Study of Parents and Children (ALSPAC). ALSPAC recruited pregnant women in the Bristol and Avon area, United Kingdom, with an expected delivery date between April 1991 and December 1992 (http://www.bris.ac.uk/alspac/). Over 14,000 pregnancies have been followed up (both children and parents) throughout the life-course. Full details of the cohort has been published previously (Fraser et al. 2013). This study uses phenotypic and DNA methylation data from the mothers (N = 900).

10 continuous traits were extracted from the same timepoint blood was drawn for DNA methylation measurements. A summary of the phenotypes are present in the __Supplementary Material__ and full details of all the data is available through a fully searchable data dictionary: http://www.bris.ac.uk/alspac/researchers/data-access/data-dictionary/

These traits do not neccessarily represent independent phenotypes and as such we wanted to prevent correlated traits skewing results. The absolute Pearson's correlation coefficient between each trait was subtracted from one (1 - absolute r). Then traits were greedily selected having 1 - aboslute r < 0.2 with any other trait.

Ethical approval for the study was obtained from the ALSPAC Ethics and Law Committee and from the UK National Health Service Local Research Ethics Committees. Written informed consent was obtained from both the parent/guardian and, after the age of 16, children provided written assent.

#### DNA methylation data
DNA methylation were measured using the Illumina Infinium HumanMethylation450 (HM450) BeadChip. The data went through quality control and was normalised before use. Full details can be found in the __Supplementary Material__.

DNA methylation data generated from blood collected at a single clinic visit was used for each of the participants. 

Cell proportions (CD8+ and CD4+ T cells, B cells, monocytes, natural killer cells, and granulocytes) were estimated using an algorithm proposed by Houseman et al. (Houseman et al. 2012).

### Surrogate variable analyses

#### Parameters tested
* SVA package used - between sva and smartsva
* Number of CpGs - increasing from 20,000 to 300,000 by increments of 20,000 + the full set of CpGs (4.83103 &times; 10<sup>5</sup>)  
* Number of samples - increasing from 100 to 900 by increments of 100
* Number of SVs - comparison of 5, 10, 15, 20, 60 SVs

#### SVA package
Time taken to run the analyses was compared as were the SVs generated by assessing the Pearson's correlation between SVs and how much of the variance of each SV generated using the sva package was explained by the same SV generated using the smartsva package.

We set the number of iterations used by each package to be the same and kept the other parameters as their defaults. 
 
Where unspecified, smartsva was chosen to run analyses.

#### Effectiveness of subsetting CpGs
Time taken to run the analyses and the SVs generated by subsets of CpGs were compared to the same factors when running the analyses using all CpGs on the 450k array. To compare the SVs we took each of the SVs generated using all CpGs on the 450k array and calculated the variance of each of these SVs explained by all the SVs generated using the subset of CpGs. i.e:

$$SV_i \sim SV_1' + SV_2' + ... + SV_n'$$

where _SV~i~_ is the ith SV generated using all CpGs from the 450k array. 

For the number of SVs and the number of samples we just assessed run time.

#### Number of SVs required

The variance explained of often used covariates in EWAS by SVs was examined. The covariates chosen were: cell proportions derived using the Housemann method, two batch variables (BCD_id and MSA4Plate_id), age, and the top 10 genomic PCs. 

Normally distributed, random continuous and binary variables were used as well as 10 randomly selected, independent (-0.2 < r < 0.2) variables from the FOM1 timepoint.

Firstly, the number of SVs required was estimated using the num.sv() function from the sva package, with the "leek" method and by using random matrix theory as suggested by the smartsva package authors. The highest number of SVs estimated from either method was used as the maximum number of SVs required. 

Then, these SVs were generated for each trait using all of the probes on the 450k array.

Finally, the cumulative variance explained by SVs for each covariate was examined. 

#### Replication and whole genome estimates

Number of CpGs required to generate 20 SVs for an EPIC array dataset and X other 450k datasets from GEO was also examined.

Finally we designed a function that suggests how many CpGs should be used given the number of SVs needed and the number of CpGs in total.

#### Analysis
All analyses were done using R (XX), with the SmartSVA and SVA packages.

## Results
The correlation between each of 20 SVs produced by the SmartSVA and SVA package was 1 or -1 when the parameters used were the same (__Supplementary Table 1__) and the time taken by SmartSVA was faster when the sample size and number of CpGs were at the maximum, median = 4.3 times faster. When using the default alpha and B values in the SmartSVA package there was little change in time taken to using SmartSVA with the default parameters from the SVA package: __Supplementary Figure 2__.

The time taken to produce 20 SVs increased with number of samples and CpGs used to create the SVs (__Figure 1__). However, the number of SVs created and whether the phenotype distribution was binary or continuous didn't effect the time very much (__Supplementary Figure 5__).

![plot of chunk time_plot](figure/time_plot-1.png)

__Figure 1__: Differences in time taken to finish running surrogate variable analysis with the smartsva.rcpp function in the SmartSVA package when varying the number of CpGs and samples.

When subsetting the 450k array, the number of CpGs required to accurately estimate SVs was examined. When producing 20 SVs, over 90% of the variance of SV 1- produced using the entire array was explained by all 20 SVs produced using a random subset of 20,000 CpGs. As the SV number increased, the variance explained by all SVs estimated from a subset of CpGs decreased (__Figure 2__). To explain over 90% of variance of all SVs,  CpGs were required. Subsetting to the most variable CpGs produced SVs that explained less of the variance than when subsetting to random CpGs (__Supplementary Figure 3__).

Similar findings were observed across X other 450k datasets extracted from GEO __either a table in main manuscript or figures in supplementary__

Using a dataset with DNA methylation measured with the EPIC array we observed X.

![plot of chunk ncpg_plot](figure/ncpg_plot-1.png)

__Figure 2__: Variance captured of each of the first 20 SVs made using 450k CpGs by all 20 SVs made using random subsets of varying numbers of CpGs. 

The num.sv() function from the SVA package estimated fewer SVs for each of the X phenotypes tested, with a median of X fewer, compared to estimates using random matrix theory (__Supplementary Figure 4__). 

The variance explained of various common EWAS covariates by the number of SVs estimated by random matrix theory varied substantially (__Figure 3__). As more SVs were added into the model, less additional variance was explained. For the 10 genomic PCs, all the SVs explained nearly none of the variance, but for the other covariates, a varying amount of variance was explained, which was fairly consistent across all the real and simulated traits assessed (__Supplementary Figure 5__). <!-- The variance explained by 10, 20, 30 and the total number of SVs for each covariate for a simulated continuous trait and a real trait are shown in __Table 1__. -->

![plot of chunk covs_var_exp](figure/covs_var_exp-1.png)

__Figure 3__: Variance of often used covariates in EWAS explained by the number of surrogate variables suggested by random matrix theory.

<!-- ```{r nsv_summary, results = "asis"}
```
 -->
## Discussion
Surrogate variable analyses has been shown to adequately control for unwanted variation in EWAS. However, a method to conduct SVA that increases computational efficiency, whilst maintaining the desired SV function is needed as datasets increase in size. Here we demonstrate that by randomly subsetting the number of CpGs to X% of the number of CpGs used for the analyses to calculate 20 SVs using the SmartSVA package, the time taken decreases by X-fold, compared to using the full set of CpGs and the SVA package, whilst maintaing SV functionality.



How does time taken to perform SVA vary with sample number, SVA package and CpG number?


# Supplementary Material

## Results



What is the correlation between SVs produced by smartsva and sva?

![plot of chunk model_param_time_plot](figure/model_param_time_plot-1.png)

How does time taken to perform SVA vary with alpha and B values?

![plot of chunk time_plot2](figure/time_plot2-1.png)

How does time taken to perform SVA vary with number of SVs and distribution of phenotype (continuous or binary)?

![plot of chunk est_num_plot](figure/est_num_plot-1.png)

Number of SVs required estimated using num.sv() and random matrix theory
