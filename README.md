# Whole genome analysis pipeline


## Sections

   * [Repository Description](#repository-description)
   * [WGS considerations](#wgs-considerations)
   * [QC pipeline steps](#qc-pipeline-steps)
   1. [Dataset QC prior to manipulating VCF](#dataset-qc-prior-to-manipulating-vcf)
   2. [Dataset pre-filtering](#dataset-pre-filtering)
   3. [High quality hard call subset of the data](#high-quality-hard-call-subset-of-the-data)
   4. [Outlier sample QC part 1](#outlier-sample-qc-part-1)
   5. [Sex check sample QC](#sex-check-sample-qc)
   6. [Identity-by-descent filtering](#identity-by-descent-filtering)
   7. [Principal components analysis](#principal-components-analysis)
   8. [Principal components filtering](#principal-components-filtering)
   9. [Outlier sample QC part 2](#outlier-sample-qc-part-2)
   10. [Variant QC](#variant-qc)
   11. [Assessing variant QC](#assessing-variant-qc)
   * [Sample QC Parameters](#sample-qc-parameters)
   * [Variant QC Parameters](#variant-qc-parameters)


<br/><br/><br/>

## Repository description

__This repo provides a template for processing WGS data in Hail__
  * __README__: Considerations of WGS data and QC pipeline steps
  * __Resources__: Links and other things relevant to WGS pipelines
  * __Modules__: Hail script templates for each pipeline step

Inspired by gnomAD blog posts:
  * [gnomAD v2.1](https://macarthurlab.org/2018/10/17/gnomad-v2-1/)
  * [gnomAD v3.0](https://macarthurlab.org/2019/10/16/gnomad-v3-0/)

<br/><br/><br/>

## WGS considerations

__File size:__
  * WGS data is orders of magnitdues larger than GWAS array (fixed variant size) and exome capture sequencing (~1% of genome)
  * File size will scale with sample size
  * Strategies to deal with large files:
  	* Use [Hail](https://hail.is/)
  	* Run everything on the [cloud](https://sites.google.com/a/broadinstitute.org/atgu/google-cloud-platform-start-up)
  	* [Re-partitioning](https://hail.is/docs/0.2/hail.MatrixTable.html#hail.MatrixTable.repartition) the Hail matrix table
  	* Use a [Sparse matrix table](https://hail.is/docs/0.2/experimental/vcf_combiner.html#working-with-sparse-matrix-tables) 

__File formats:__
  * Per-sample gVCFs: VCF format data for each sample, where all variant sites and reference blocks are listed for all contigs
  * Ref-blocks: Summary depth/genotype quality information in the gVCF across a genomic interval where the sample has no variant sites

__Common and rare variation:__
  * WGS data contains SNP/indel calls across the allele frequency spectrum, whereas array (MAF > 0.1%)
  * Low coverage WGS (< 10x) aided by imputation
  * High coverage WGS (> 20x) robust across all SNP/indels
  * Multiple levels of analysis
  	* Pulling out common variants to include in GWAS meta-analysis (GWAMA)
  	* Pulling out rare coding variants to include in Exome meta-analysis
  	* Looking at rare non-coding variation (good luck!)

__Repetitive regions:__
  * WGS contains challenging genomic regions often ignored in array and exome capture sequencing
  * Regions with highly repetitive sequence that is difficult for short read sequening to align properly to the reference genome library
  	* Telomeres/centromeres of each chromosome
  	* Segmental duplications - large repetitive chumks
  	* Low complexity regions (LCRs)
  * These regions are often excluded early on in an analysis

<br/><br/><br/>

## QC pipeline steps

### Index of steps:

   1. [Dataset QC prior to manipulating VCF](#dataset-qc-prior-to-manipulating-vcf)
   2. [Dataset pre-filtering](#dataset-pre-filtering)
   3. [High quality hard call subset of the data](#high-quality-hard-call-subset-of-the-data)
   4. [Outlier sample QC part 1](#outlier-sample-qc-part-1)
   5. [Sex check sample QC](#sex-check-sample-qc)
   6. [Identity-by-descent filtering](#identity-by-descent-filtering)
   7. [Principal components analysis](#principal-components-analysis)
   8. [Principal components filtering](#principal-components-filtering)
   9. [Outlier sample QC part 2](#outlier-sample-qc-part-2)
   10. [Variant QC](#variant-qc)
   11. [Assessing variant QC](#assessing-variant-qc)

<br/><br/>

### Dataset QC prior to manipulating VCF

  * **GOAL: Understanding the project/phenotype data you are working with**
    * List and understand all sample phenotypes provided
    * Match up phenotype file IDs with genetic data IDs
	* Resolve any inconsistencies before moving forward
    * Start spreadsheet/table of datasets, sequence platforms 
    * Know what genome reference your sequence data is mapped (most WGS data is hg38 / GRCh38)
    * Determine whether you are using a dense or sparse matrix table
    * Look over the VCF meta-data at the top of the VCF file
	* Read VCF and phenotype/sample file into Hail
    * Generate sample QC metrics from raw VCF
    * Write up paragraph of sample collection and descriptives

**Categorical variables often analyzed in datasets**
  * Cohort
  * Sequencing wave
  * C-Project (Broad specific)
  * Sequencing Plate
  * Sex
  * Affection status
  * Continental ancestry groupings / reported ancestry
 
**Quantitative parameters often analyzed in datasets**
  * Age
  * Top genetic principal components
  * Number of variant sites (n_snps)
  * Singleton rate (n_singleton)
  * Het / hom variant ratio (r_het_hom_var)

<br/><br/>

### Dataset pre-filtering
  * **GOAL: Remove variants that are highly unlikely to be analyzed**
    * Remove variants that fail VQSR (or AS-VQSR)
    * Remove variants in telomeres / centromeres
    * Remove variants in low-complexity regions
    * Remove variants in segmental duplication regions
    * Generate sample QC metrics from pre-filtered VCF
    * VEP annotate remaining sites    

<br/><br/>

### High quality hard call subset of the data
  * **GOAL: Use a smaller subset of variants (i.e. 50-500k variants) to analyze relatedness and population ancestry**
    * Bi-allelic SNVs only
    * Common variants (MAF > 0.1%)
    * Call rate > 99%
    * LD-pruned with a cutoff of r2 = 0.1  
  * If running from raw VCF calls
    * Run split.multi to maximize bi-allelic sites
    * Filter to PASS sites and SNV 
    * Run variant.qc to get MAF

<br/><br/>


### Outlier sample QC part 1
  * **GOAL: Remove samples that are contaminated or have poor sequencing levels**
    * Use pre-filtered dataset
    * Plot values below before using DEFAULT filters to ensure you are not throwing away large amounts of samples
      * freemix contamination filtering (DEFAULT > 0.05)
      * Chimeric read filtering (DEFAULT > 0.05 )
      * Call rate filtering (DEFAULT < 0.85)
      * Mean Depth coverage filtering (DEFAULT < 15)
      * Small insert size: Median insert size < 250bp

<br/><br/>

### Sex check sample QC
  * **GOAL: remove samples where genotype sex does not equal reported sex**
    * Filter out variants within PAR coordinates (https://en.wikipedia.org/wiki/Pseudoautosomal_region)
      * Reported males shoud have X chromosome F-statistic from 0.8 to 1
      * Reported females shoud have X chromosome F-statistic from -0.2 to 0.4
      * Remove samples with ambiguous imputed sex (X chromosome F-statistic between 0.4 and 0.8)
    * Large-scale sex check errors are indicative of ID mismatching or upstream data mishandling

<br/><br/>

### Identity-by-descent filtering
  * **GOAL: remove 1st and 2nd degree relatives from population based sample**
    * Use HQ hardcall dataset
    * Consider which IBD algorithm to use ([KING](https://pubmed.ncbi.nlm.nih.gov/20926424/), [PC-AiR](https://pubmed.ncbi.nlm.nih.gov/25810074/), [PC-relate](https://www.rdocumentation.org/packages/GENESIS/versions/2.2.2/topics/pcrelate))
    * Which algorithm you use will depend on sample size, ancestry composition, and knowledge of family pedigrees in data
    * Unsure? Start by trying genetic relatedness using Hail [pc-relate](https://hail.is/docs/0.2/methods/relatedness.html#hail.methods.pc_relate)
    * Plot proportion of 0 shared vs 1 shared alleles
    * IBD filtering on PI-HAT value (DEFAULT > 0.2)

<br/><br/>

### Principal components analysis
  * **GOAL: Determine general ancestry of cohort**
    * Use HQ hardcall dataset
    * Consider which PCA strategy to use
      * PCA comparing against self-identified ancestry information (PROS: easy if you have this info CONS: dependent on quality of self-id info) 
      * PCA including reference datasets with known ancestry (PROS: good truth dataset CONS: can be hard to harmonize variants)
      * PCA using gnomAD PCs (PROS: no need to merge variants CONS: shrinkage in PCs due to variant mismatch
    * Unsure? Start with Hail [hwe_normalized_pca](https://hail.is/docs/0.2/methods/genetics.html#hail.methods.hwe_normalized_pca)
    * Run with and without reference panel data (1KG or 1KG+HGDP genome reference panel)
    * Create plots of PCs
      * Case / control coloring
      * Cohort coloring
    * Assigning samples to a particular ancestry
      * [Random forest model](https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html)

<br/><br/>

### Principal components filtering
  * **GOAL: match case and controls within a common genetic ancestry**
    * If retaining multiple ancestries, make sure to define ancestry groups in phenotype file
    * PCA filtering (no DEFAULT filtering parameters)
      * 2-dimensional centroid approach (using distance from 1KG ancestry mean PC; used by Chia-Yen Chen)
      * pair-matching cases to controls (R package: optmatch; R function: pairmatch(); used by Mitja Kurki)
      * Within each ancestry group, re-run PCA
		* Re-evaluate PC dispersion  
		* Add these PCs as covariates to phenotype file

<br/><br/>

### Outlier sample QC part 2
 * **GOAL: remove samples within continental ancestry groupings that have unusual variant properties**
   * Examine variation in:
	* TiTv ratio
	* Het/Hom ratio 
	* Insertion/Deletion ratio
	* Plots with colors defined by assigned ancestry
		* Different ancestries have significant mean differences in these ratios
	* Filter out within cohort outliers (DEFAULT > 4 Std. deviations within a particular ancestry)

<br/><br/>

### Variant QC
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
    * Generate sample QC metrics

<br/><br/>

### Assessing variant QC
  * **GOAL: Determine if more stringent variant QC is needed**
    * Examine QC parameters across 3 filtering steps:
	    * Pre-filtered VCF
	    * Sample QC'ed VCF
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

<br/><br/>







## Sample QC Parameters

**From GATK/Picard metadata**
  * Freemix Contamination
  * Chimeric read percentage

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

<br/><br/><br/>

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

<br/><br/><br/>






# Case/control Whole Genome Sequencing - Burden and Association pipeline

Basic Requirements:
 
 * A QC-ready annotated matrix table
 * Samples annotated with PCs

## Whole Genome Burden

Primary Hypotheses:

 * Do cases have a higher overall burden of rare deleterious variants than controls?
	 - Does this burden decrease as allele frequency increase?
 * Is the burden concentrated in genes with evidence of selective constraint?
	 - Does the burden effect attenuate outside of these constrained genes?
 * Sanity check: Is the rate of rare variants similar in cases and controls?
	 - Highly significant differences suggest QC issues still persist

General Method:

 * Aggregate per-individual counts on selected annotation/allele frequency
 * Testing deleterious coding variant count using logistic regression
 * Include first 10 PCs and sex as covariates
	
Comparisons:

 * Stratify by cohort
	 - Require cases and controls within any cohort designation 
 * Stratify by commonly assessed variant annotations
	 - Genic / non-genic
	 - Conservation/constraint
	 - CADD severity prediction
	 - Damaging missense (PolyPhen damaging, SIFT deleterious, MPC > 2)
	 - Protein truncating variants (frameshift_variant, splice_acceptor_variant, splice_donor_variant, stop_gained)
 * Stratify by Allele frequency
	 - Ultra-rare variation: non-gnomAD singletons
	 - Dataset + gnomAD AC < 5
	 - Doubleton distributions

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

 * Example Case/Control WES Portal: [SCHEMA](https://schema.broadinstitute.org)

## Basic association tests:

 * Fisher's Exact test
	 - Separate individuals into discrete binary categories (carrier / non-carrier)
	 - fisher.test() in R
	 - in Hail: http://discuss.hail.is/t/simple-rare-variant-burden-testing-with-fisher-exact-test/210

 * Poisson rate test
	 - Compares variant counts when the mean and variance are equal (generally for de-novo / ultra-rare variants only)
	 - poisson.test() in R

 * Logistic (or Firth regression) in Hail
	 - Predicts case/control status by allele frequency and allows for covariates
	 - Requires some asymptotic assumptions, which can be unmet at low allele counts (< 20)

## Mixed model association tests

 * [SAIGE](https://github.com/weizhouUMICH/SAIGE)
 	 - Mixed-model assocation with saddle-point approximation (SPA) to control for imbalanced case/control ratio
 	 - Both individal variant and gene-based tests (using SKAT) available

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


