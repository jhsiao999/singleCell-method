---
layout: page
title: "Single-cell sequencing data analysis and methods"
tagline: 
---

### Home
  * [Admixture clustering](#admixture-clustering)
  * [Cell-cycle phase assignment](#assign-cell-cycle)
  * [cellcycleR](#cellcycleR)

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

### Cell-cycle phase assignment (Macosko et al., 2015) <a id = 'assign-cell-cycle'></a>
* [Cell ordering of the iPSCs - Macosko method](project/analysis/cell_ordering_iPSC.html)
* [Cell cycle genes scores - recovered cell order](project/analysis/cell_cycle_score_analysis.html)

---

### cellcycleR applications <a id = 'cellcycleR'></a>
* [Oscope data, Leng et al., 2015](project/analysis/oscope_cellcycleR.html)
* [Mouse ESCs, Buettner et al., 2015](project/analysis/marioni_cellcycleR.html)
* [Human iPSC data, Tung et al., 2015](project/analysis/yoav_cellcycleR.html)
	* [Gene annotations](project/analysis/yoav_cellcycleR_postprocessing.html)
* [Human LCL data, Tung et al., 2015](project/analysis/lcl_cellcycleR.html)
* [Monocle data, Trapnell et al., 2014](project/analysis/monocle_cellcycleR.html)
* [Botstein yeast data (cdc, elu and alpha)](project/analysis/yeast_cellcycleR.html)


