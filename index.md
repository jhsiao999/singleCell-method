---
layout: page
title: "Single-cell sequencing data analysis and methods"
tagline: 
---

### Home
  * [Admixture clustering](#admixture-clustering)
  * [Cell-cycle phase modeling](#assign-cell-cycle)
  * [Batch effect](#batch)

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

#### [Methods and Materials](project/docs/cell_reorder.pdf)

#### Applications

* [Human iPSC data] - ( [Macosko method](project/analysis/cell_ordering_iPSC.html), [explore cell order](project/analysis/cell_cycle_score_analysis.html) )
* Human iPSC data (Tung et al. 2015)
  * [On annotated cell cycle genes](project/analysis/yoav_cellcycleR_cellcycle_genes.html)
     * [Gene annotations](project/analysis/yoav_cellcycleR_postprocessing_cellcycle_genes.html)
  * [On all genes](project/analysis/yoav_cellcycleR_all_genes.html)
     * [Gene annotations](project/analysis/yoav_cellcycleR_postprocessing_all_genes.html)
	* [On non ribosomal genes + gene annotations](project/analysis/yoav_cellycleR_non_ribosomal.html)  
	* [On CDC, cyclin and cell cycle genes](project/analysis/yoav_cellcycleR_cdc_cyclin.html)
* [Human LCL data, Tung et al., 2015](project/analysis/lcl_cellcycleR.html)
* [Mouse ESCs, Buettner et al., 2015](project/analysis/marioni_cellcycleR.html)
* [Oscope data, Leng et al., 2015](project/analysis/oscope_cellcycleR.html)
* [Monocle data, Trapnell et al., 2014](project/analysis/monocle_cellcycleR.html)
* [Botstein yeast data (cdc, elu and alpha)](project/analysis/yeast_cellcycleR.html)

#### Nonparametric smoothing <a id="nonparametric-cellcycler"></a>
* [Wavelets](project/analysis/wavelet_validation_check.html)
* [SMASH (smoothing wavelets)](project/analysis/smash_validation_check.html)
* [LOESS (locally weighted scatterplot smoothing) ](project/analysis/loess_validation_check.html)
* [Smoothing splines](project/analysis/splines_validation_check.html)
* [Testing effectiveness of nonparametric cellcycleR](project/analysis/nonparametric_cellcycleR_tests.html)
* [Nonparametric vs Sinusoidal cellcyclerR]
  * [Case Study 1](project/analysis/cellcycleR_compare1.html)
  * [Case Study 2](project/analysis/cellcycleR_compare2.html)
  * [Case Study 3](project/analysis/cellcycleR_compare3.html)
  * [Case Study 4](project/analysis/cellcycleR_compare4.html)


---

### Batch effect <a id = 'batch'></a>
* [ERCC genes](project/analysis/ercc-pca.html)
  * [RUVg](project/analysis/ercc-ruvg.html)
  * [SVA package](project/analysis/ercc-sva.html)
* [Bulk RNA]
  * [RUVg zibrafish](project/analysis/ercc-ruvg-paper-data.html)



