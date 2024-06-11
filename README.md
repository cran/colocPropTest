# colocPropTest

An R package for proportionality testing for colocalisation.  

Colocalisation is the occurence of two traits sharing a causal variant in a given genetic region. One popular method for inferring whether two traits have a colocalising signal is "coloc" [1].  This is based on fine mapping the summary statistics for two traits, and performing some joint analysis. However, fine mapping can be inaccurate when information varies between SNPs (eg due to varying sample size), particularly in large samples, and this inaccuracy can propogate to coloc.

This motivated us to consider an older method of testing colocalisation, based on testing proportionality of regression coefficients between two traits [3,4]. There were a number of issues with this approach which motivated the development of coloc:

1. it tests the null hypothesis that coefficients are proportional, and thus is difficult to interpret because failure to reject can either correspond to colocalisation, or lack of power

2. the test is based on the coefficients from a joint multi-SNP model, which can be hard to reconstruct accurately from the marginal test statistics typically available

2. the degrees of freedom of the test is n-1, where n is the number of SNPs, which can lead to lack of power

Here we solve those issues by

1. only advising to conduct the test when standard evidence points to convincing association with both traits (eg minimum p values $< 10^{-8}$)

2. using marginal coefficients directly, with the LD between the variants used to infer the covariance of those test statistics

3. performing many tests of pairs of variants rather than one test of n variants, and using false discovery rates to infer whether any of those tests were true rejections of the null of colocalisation

## References

1. Giambartolomei, C. et al. Bayesian Test for Colocalisation between Pairs of Genetic Association Studies Using Summary Statistics. PLOS Genet. 10, e1004383 (2014).

2. Kanai, M. et al. Meta-analysis fine-mapping is often miscalibrated at single-variant resolution. Cell Genomics 2, 100210 (2022).

3.	Plagnol, V., Smyth, D. J., Todd, J. A. & Clayton, D. G. Statistical independence of the colocalized association signals for type 1 diabetes and RPS26 gene expression on chromosome 12q13. Biostatistics 10, 327–334 (2009).

4.	Wallace, C. Statistical Testing of Shared Genetic Control for Potentially Related Traits. Genet. Epidemiol. 37, 802–813 (2013).

