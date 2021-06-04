# whole_genome_analysis_pipeline



* [June 2017: WES QC Process (last updated June 2017)](#wes-qc-process)
   1. [Files](#files)
   2. [QC comparisons](#qc-comparisons)
   3. [Sample QC Parameters](#sample-qc-parameters)
   4. [Variant QC Parameters](#variant-qc-parameters)
   5. [QC pipeline steps](#qc-pipeline-steps)
* [QC pipeline examples](#qc-pipeline-examples)


## QC comparisons

**Categorical Groupings**
  * Cohorts/waves
  * C-Project
  * Capture platform
  * Sex
  * Affection status
  * Project specific categories

**Quantitative parameters**
  * Top Principal components
  * % Coverage (callRate)
  * Mean Depth (dpMean)
  * Singleton Synonymous rate (nSingleton when restricting to synonymous variants)
  * Non-ExAC / non-discovEHR singleton rate
        * 92k ExAC/GnomAD sites with Psych-exomes and MiGEN removed: `gs://exome-qc/resources/gnomad.r2.0.1/gnomad.r2.0.1.nonpsych.variants.tsv`
	* discovEHR sites `/humgen/atgu1/fs03/shared_resources/discovEHR/GHS_Freeze_50.L3DP10.pVCF.frq.hail.vcf`
	* ExAC sites `/humgen/atgu1/fs03/shared_resources/ExAC/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz`
	* Non-Psych ExAC sites `ftp://ftp.broadinstitute.org/distribution/ExAC_release/release1/subsets/ExAC_nonpsych.r1.sites.vep.vcf.gz`


## Sample QC Parameters

**From GATK/Picard metadata**
  * Freemix Contamination
  * Chimeric read percentage
  * PCT_TARGET_BASES_10X
  * PCT_TARGET_BASES_20X

**From Hail Sample QC Annotation Table**
**Primary QC parameters** 
  * callRate
  * rTiTv
  * rHetHomVar
  * rInsertionDeletion
  * nSingleton
  * dpMean

**Secondary QC parameters**
  * nCalled
  * nNotCalled
  * nHomRef
  * nHet
  * nHomVar
  * nSNP
  * nInsertion
  * nDeletion
  * nTransition
  * nTransversion
  * dpStDev
  * gqMean
  * gqStDev
  * nNonRef


## Variant QC Parameters

**Variant Quality**
  * VQSLOD - Variant Quality Score in LOD space (or new RF posterior probabilities)
  * pHWE - Hard-Weinberg Equilibrium
  * AC - allele count
  * Median depth
  * QD - quality by depth

**Genotype Quality**
  * Depth
  * PHRED likelihood (PL)
  * Genotype Quality (GQ - same as PL in joint-called VCF)
  * Allele Depth (AD)


## QC pipeline steps

1. **Dataset QC prior to menipulating VCF**
	* **GOAL: Understanding the project/phenotype data you are working with**
	* List and understand all sample phenotypes provided
	* Match up phenotype file IDs with VCF IDs
		* Resolve any inconsistencies before moving forward
	* Start spreadsheet/table of datasets, capture platforms
		* Split data up into category with the most designations 
	* Write up paragraph of sample collection and descriptives
  * Know what genome reference your sequence data is mapped
  * Look over the VCF meta-data at the top of the VCF file
	* Read VCF and phenotype/sample file into Hail
  * Generate sample QC metrics from raw VCF

2. **Variant QC pre-filtering**
	* **GOAL: Remove variants that are highly unlikely to be analyzed**
    * Remove variants that fail VQSR
    * Remove variants in telomeres / centromeres
    * Remove variants in low-complexity regions
    * Generate sample QC metrics from pre-filtered VCF
	  * VEP annotation    

3. Creating a high quality hard call subset of the data
  * **GOAL: Use a smaller subset of variants to analyze relatedness and population ancestry**
  
  ........work from here

2. **Outlier sample QC: part 1**
	* **GOAL: Remove samples that are contaminated or have poor sequencing levels**
	* Use pre-filtered VCF
	* Plot values below before using DEFAULT filters to ensure you are not throwing away large amounts of samples
    * freemix contamination filtering (DEFAULT > 5%)
    * Chimeric read filtering (DEFAULT > 1.4% )
    * Call rate filtering (DEFAULT < 90%)
    * Mean Depth coverage filtering (DEFAULT < 20)

3. **Sex check**
	* **GOAL: remove samples where genotype sex does not equal reported sex**
	* Filter out variants within PAR coordinates (https://en.wikipedia.org/wiki/Pseudoautosomal_region)
	* Reported males shoud have X chromosome F-statistic from 0.8 to 1
	* Reported females shoud have X chromosome F-statistic from -0.2 to 0.4
    * Remove samples with ambiguous imputed sex (X chromosome inbreeding coefficient between 0.4 and 0.8)
	* Large-scale sex check errors are indicative of ID mismatching or upstream data mishandling

4. **Principal components analysis**
	* **GOAL: Determine general ancestry of cohort**
	* Use raw VCF (so as not to exclude any common variants)
	* Subset to common variant VCF (use 9k set: pruned_9k_common_variants_t.bim or generate your own)
    * Run Hail PCA with 1K Genomes samples
    * Create plots of PCs
	    * Case / control coloring
	    * Cohort coloring
	* Assigning samples to a particular ancestry
		* Random forest approach (used by TJ Singh)
		* SNPWEIGHTS: https://www.hsph.harvard.edu/alkes-price/software/ (used by Chia-Yen Chen)

5. **Outlier sample QC: part 2**
 	* **GOAL: remove samples within cohort that have unusual variant properties**
	* Use pre-filtered VCF
	* Examine per-cohort variation in:
		* TiTv ratio
		* Het/Hom ratio 
		* Insertion/Deletion ratio
	* Plots with colors defined by assigned ancestry
		* Different ancestries have significant mean differences in these ratios
	* Filter out within cohort outliers (DEFAULT > 4 Std. deviations within a particular ancestry)

6. **Principal components filtering**
    * **GOAL: match case and controls within a common genetic ancestry**
	* Run Hail PCA without 1K Genomes samples 
	* If retaining multiple ancestries, make sure to define ancestry groups in phenotype file
	* PCA filtering (no DEFAULT filtering parameters)
		* 2-dimensional centroid approach (using distance from 1KG ancestry mean PC; used by Chia-Yen Chen)
		* pair-matching cases to controls (R package: optmatch; R function: pairmatch(); used by Mitja Kurki)
	* Within each ancestry group, re-run PCA
		* Re-evaluate PC dispersion  
		* Add these PCs as covariates to phenotype file

7. **Identity-by-descent (IBD) filtering**
    * **GOAL: remove 1st and 2nd degree relatives from population based sample**
    * Within each ancestry group, calculate IBD using common variant VDS
    * Plot proportion of 0 shared and 1 shared alleles
    * IBD filtering on PI-HAT value (DEFAULT > 0.2)

8. **Variant QC**
	* **GOAL: Remove low quality/somatic variants**
	* Use pre-filtered VCF with outlier samples removed
    * Filter variants with low call rate (DEFAULT < 95%)
	    * Split variant call rate by capture, case/control status, or other category
    * Remove monoallelic variants: no alt or no ref calls
	* Genotype filtering (set to filtered variants to missing)
		* Filter by Depth (DEFAULT < 10)
		* Filter by GQ (DEFAULT < 25)
		* HET calls: Filter by allele depth/balance (AB/AD; DEFAULT < 0.25)
			* Additional pAB filter on higher depth (binomial test - DEFAULT p < 1e-9)
		* HOMREF calls: (AB/AD; DEFAULT > 0.1)
		* HOMALT calls: (AB/AD; DEFAULT < 0.9)
    * Remove variants not in HWE (DEFAULT p < 1e-6)
    * Generate sample QC metrics from variant-QC'ed VCF

9. **Assessing variant QC**
    * **GOAL: Determine if more stringent variant QC is needed**
    * Examine QC parameters across 3 filtering steps:
	    * Raw VCF
	    * Pre-filtered VCF
	    * Variant QC'ed VCF
	* QC parameters:
		* Number of SNPs / Indels
		* TiTv ratio
		* Het/Hom ratio
		* Ins/Del ratio
		* Singleton synonymous rate
	* Primary categories:
		* Case/control (should be equal between groups)
		* Cohort (should vary predictably)
	* Determine whether additional variant filtering needs to be done

10. **Final sample QC**
	* **GOAL: See if any samples are outliers after variant QC**
    * Call rate
    * Median Depth
    * TiTv
    * Singleton synonymous rate



## QC pipeline examples

**Andrea Ganna's steps for running QC and analysis in WGS data:**

[Analysis Plan](https://storage.googleapis.com/wgspd-wgs-v1-noshared/Analysis_plan_METSIM.md.html)

Includes:
 - Normalizing phenotypes
 - pyhail quality control
 - VEP annotation
 - Score generation
 - Association analysis

[WGSPD](https://storage.googleapis.com/wgspd-wgs-v1/wgspd_guide_shared.html)

Includes:
 - Data structure on the cloud
 - Downloading files/data from the cloud
 - Running Hail on the cloud
 - QC filtering
 - Generating PCs
 - Annotation 











# Case/control Whole Exome Sequencing - Burden and Association pipeline

Basic Requirements:
 
 * A QC-ready annotated VCF/VDS or set of VCF/VDS files
 * Samples annotated with PCs

## Whole Exome Burden

Primary Hypotheses:

 * Do cases have a higher overall burden of rare deleterious variants than controls?
	 - Does this burden decrease as allele frequency increase?
 * Is the burden concentrated in genes with evidence of selective constraint?
	 - Does the burden effect attenuate outside of these constrained genes?
 * Sanity check: Is the rate of rare synonymous variants similar in cases and controls?
	 - Highly significant differences suggest QC issues still persist

General Method:

 * Aggregate per-individual counts on selected annotation/allele frequency
 * Testing exome-wide variant count using logistic regression
 * Include first 10 PCs and sex as covariates
	
Comparisons:

 * Stratify by cohort
	 - Require cases and controls within any cohort designation 
 * Stratify by primary coding annotations
	 - synonymous
	 - non-damaging missense
	 - damaging missense (PolyPhen damaging and SIFT deleterious)
	 - Protein truncating variants (frameshift_variant, splice_acceptor_variant, splice_donor_variant, stop_gained)
 * Stratify by Allele frequency
	 - URV: non-ExAC singletons
	 - other rare variant categories: ExAC AC=1, ExAC AC=2-10, ExAC AC > 10

Infomative Graphics:

 * Forest plots of Odds ratio, 95% CI, and p-value


## Gene-based association

Primary Hypotheses:

 * Do any genes associate with our phenotype after correcting for multiple testing?
 * Is there inflation of the median test statistic (lambda)
	 - Inflation suggests that there is underlying population stratification/QC not accounted for
	 - Can also mean a polygenic signal in well-powered cohorts

General Method:

 * Aggregate per-individual counts on selected annotation/allele frequency
 * Testing per-gene variant count using logistic regression
 * Include first 10 PCs and sex as covariates

Infomative Graphics:

 * QQ plots of expected/actual p-values [example .ipynb](https://github.com/Nealelab/epi25/blob/master/notebooks/variant_association.ipynb)
 * Gene-level Manhattan plot [example .R function](https://github.com/Nealelab/array_CNV_analysis/blob/master/Manhattan.R)
 * Forest plots of Odds ratio, 95% CI, and p-value


## Portal information

 * Example Case/Control WES Portal: [ibd.broadinstitute.org](https://ibd.broadinstitute.org)
 * See the **WES_meta-analysis_results_file.tsv** for a guideline for the per-variant / per-gene results information required for the disease web browser interface


## Basic association tests:

 * Fisher's Exact test
	 - Separate individuals into discrete binary categories (carrier / non-carrier)
	 - fisher.test() in R
	 - in Hail: http://discuss.hail.is/t/simple-rare-variant-burden-testing-with-fisher-exact-test/210

 * Poisson rate test
	 - Compares variant counts when the mean and variance are equal (generally for ultra-rare variants only)
	 - poisson.test() in R

 * Logistic (or Firth regression) in Hail
	 - Predicts case/control status by allele frequency and allows for covariates
	 - Requires some asymptotic assumptions, which can be unmet at low allele counts (< 20)
	 
 * Kernel association - SKAT test
	 - R implementation: [SKAT](https://cran.r-project.org/web/packages/SKAT/index.html)
	 - Random effects model with covariates
	 - Handles alleles at various frequencies

 * Example of logistic and fishers-exact test using Hail (courtesy of Mitja Kurki) [Hail Topic](http://discuss.hail.is/t/simple-rare-variant-burden-testing-with-fisher-exact-test/210)
 * R packages for rare variant testing
	 - [GENESIS](https://bioconductor.org/packages/release/bioc/html/GENESIS.html)
	 - [SeqMeta](https://cran.r-project.org/web/packages/seqMeta/seqMeta.pdf)

## burden_testing.Rmd

 * A walkthrough of various gene-based tests, strategies, and Hail commands for testing rare variant burden


