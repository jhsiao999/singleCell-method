---
layout: page
title: "Single-cell sequencing data analysis and methods"
tagline: 
---

### Home
  * [Admixture clustering](#admixture-clustering)
  * [Cell-cycle phase modeling](#assign-cell-cycle)
    # [Methods](#methods)
    * [Exploratory analysis](#explore)
    * [Sinusoidal fitting](#sinusoidal-cellcycler)
    * [Nonparametric smoothing](#nonparametric-cellcycler)
  * [Batch effect](#batch)
  * [Differential testing](#testing)


---

### Admixture clustering <a id = 'admixture-clustering'></a>

* Human iPSC data (Tung et al., 2015)
  * Endogeneous genes
     * [Structure plot of iPSC by batch, individual, cell cycle](project/analysis/cell_phase_analysis.html)
     * [Structure plot per individual and per individual & batch](project/analysis/structure_per_individual.html)
     * [Investigating the properties of genes - low counts](project/analysis/low_counts_genes.html)
     * [Investigating ERCC spike ins - RUV normalization](project/analysis/RUV_normalization.html)
     * [iPSC Structure plots before and after batch correction](project/analysis/batch_effect_all_genes.html)
  * Cell-cycle genes
     * [Cell-cycle scores](project/analysis/cell_cycle_score_analysis.html)
     * [Structure plots using before batch, individual correction](project/analysis/clustering_cell_cycle_genes.html)
     * [Structure plots after batch effect correction](project/analysis/batch_effect_cell_cycle_genes.html)

* LCL data (Tung et al., 2015)
	*  [Endogeneous genes](project/analysis/lcl_structure.html)
	*  [Cell-cycle genes](project/analysis/lcl_structure_cell_cycle_genes.html)
* iPSC +LCL pooled data analysis (Tung et al., 2015)
	* [Admixture analysis (All genes + Cell cycle genes)](project/analysis/ipsc_lcl_structure.html)
	* [Gene annotations of iPSC and LCL data](project/analysis/gene_annotations_ipsc_lcl.html)
	* [Patterns of counts for cluster driving genes](project/analysis/gene_patterns_iPSC_LCL.html)

* Mouse ESCs staged for cell-cycle phase (Buettner et al., 2015)
	* [Structure analysis](project/analysis/marioni_structure_all_genes.html)
	* [Gene annotations](project/analysis/gene_annotations_marioni.html)

* Mouse embryos from oocyte to blastocyst stages (Deng et al., 2015)
	* [Structure analysis](project/analysis/deng_structure_all_genes.html)
	* [Gene annotations (K=7)](project/analysis/gene_annotations_deng.html)


---

### Cell-cycle phase modeling (Macosko et al., 2015) <a id = 'assign-cell-cycle'></a>

#### Methods and Materials <a id="methods"></a>
* [Our model](project/docs/cell_reorder.pdf) 

#### Exploratory analysis <a id="explore"></a>

* [Human iPSC data] - ( [Macosko method](project/analysis/cell_ordering_iPSC.html), [explore cell order](project/analysis/cell_cycle_score_analysis.html) )

#### Sinusoidal fitting <a id="sinusoidal-cellcycler"></a>

* Human iPSC data (Tung et al. 2015)
  * Annotated cell cycle genes (2015-11-04 - [Estimation](project/analysis/yoav_cellcycleR_cellcycle_genes.html), [Gene annotations](project/analysis/yoav_cellcycleR_postprocessing_cellcycle_genes.html) ) ( 2016-01-28, NA19239 - [Estimation & annotations](project/analysis/yoav_cellcycleR_cellcycle_genes-2016-01-28.html), [Compared with PCA](project/analysis/pca-sinu-ipsc-19239.html) ) 
  * All genes (2015-11-11 - [Estimation](project/analysis/yoav_cellcycleR_all_genes.html), [Gene annotations]((project/analysis/yoav_cellcycleR_postprocessing_all_genes.html)) ) ( 2016-01-31, NA19239 - [Estimation & annotations](project/analysis/gilad-ipsc-all-genes-2016-01-31.html))

	* [On non ribosomal genes + gene annotations](project/analysis/yoav_cellycleR_non_ribosomal.html)  
	* [On CDC, cyclin and cell cycle genes](project/analysis/yoav_cellcycleR_cdc_cyclin.html)
* [Human LCL data, Tung et al., 2015](project/analysis/lcl_cellcycleR.html)
* [Mouse ESCs, Buettner et al., 2015](project/analysis/marioni_cellcycleR.html)
* [Oscope data, Leng et al., 2015](project/analysis/oscope_cellcycleR.html)
* [Monocle data, Trapnell et al., 2014](project/analysis/monocle_cellcycleR.html)
* [Botstein yeast data (cdc, elu and alpha)](project/analysis/yeast_cellcycleR.html)


* [Runtime comparisons: parallel vs. lmFit](project/analysis/sin_cell_order_iter-runtime.html)


#### Nonparametric smoothing <a id="nonparametric-cellcycler"></a>

* [Wavelets](project/analysis/wavelet_validation_check.html)
* [SMASH (smoothing wavelets)](project/analysis/smash_validation_check.html)
* [LOESS (locally weighted scatterplot smoothing) ](project/analysis/loess_validation_check.html)
* [Smoothing splines](project/analysis/splines_validation_check.html)
* [Testing effectiveness of nonparametric cellcycleR](project/analysis/nonparametric_cellcycleR_tests.html)
* [Nonparametric vs Sinusoidal cellcyclerR]
  * [Case Study 1](project/analysis/cellcycleR_compare1.html)
  * [Case Study 2](project/analysis/cellcycler_compare2.html)
  * [Case Study 3](project/analysis/cellcycleR_compare3.html)
  * [Case Study 4](project/analysis/cellcycleR_compare4.html)
* [Comparison of different nonparameteric smoothers in cellcycleR](project/analysis/nonparametric_cellcycler_methods_compare.html)

### Cell classification using classtpx <a id="classtpx"></a>

* [classtpx package](https://github.com/kkdey/classtpx)
* Simulation Runs (K=2) to validate classtpx - ([Simulation Run 1](project/analysis/classtpx_simulation_run_1.html),
[Simulation Run 2](project/analysis/classtpx_simulation_run_2.html))
* Simulation Runs (K=3) to validate classtpx - ([Simulation Run 1](project/analysis/classtpx_simulation_run_3.html),
[Simulation Run 2](project/analysis/classtpx_simulation_run_4.html))


#### Comparing to the other methods <a id="comparisons"></a>
* PCA comparison with cellcycleR - ([Simulation](project/analysis/pca_snr_compare.html), [Gilad-2015-NA19101](project/analysis/pca-sinu-ipsc-19239.html)
* [Nonparameteric smoothing - make ends meet](project/analysis/np_smoother_constraint.html)
* [Simulation study to check how cell states and ribosomal genes may affect cellcycleR](project/analysis/cellcycler_with_ribosomal_sim.html)
* [Testing Penalized Matrix Decomposition on cellcycle genes](project/analysis/pmd_cellcycler_test_1.html)


#### Correlated gene expression patterns <a id = "correlated-expression"></a>
* [Homegeneous correlation](project/analysis/gene-correlation-sinusoidal.html)

---

### Batch effect <a id = 'batch'></a>
* [ERCC genes](project/analysis/ercc-pca.html)
  * [RUVg](project/analysis/ercc-ruvg.html)
  * [SVA package](project/analysis/ercc-sva.html)
* [Bulk RNA]
  * [RUVg zibrafish](project/analysis/ercc-ruvg-paper-data.html)

---

### Differential testing <a id = 'testing'></a>
* [Exploring ECDF](project/analysis/count-cumulative.html)


### Differential testing <a id = 'differential-testing'></a>

* [Exploring ECDF](project/analysis/count-cumulative.html)
