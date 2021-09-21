# **RNA-Seq Worflow**

## Quality analysis
### Software & packages
* FastQC
* FastQ Screen
* Trimmomatic
* BUSCO
### Scripts
* fastqAnalysis.py
* runTrimmomatic.py
* clean_fasta.py
<p>&nbsp;</p>

## Transcriptome assembly
### Software & packages
* Trinity
### Scripts
* trinityMultiple.py
<p>&nbsp;</p>
  
## Mapping and gene count
### Software & packages
* STAR + HTSeq
* RSEM
### Scripts
* STARmultiple.py
<p>&nbsp;</p>
  
## Differential Expression analysis
### Software & packages
* DESeq/DESeq2
* edgeR
* NOISeq
### Scripts
* DE_analysis_DESeq2.R
* DE_analysis_edgeR.R
* DE_analysis_NOISeq.R
* DE_analysis_overlap.R
* get_DEGs_tables.py
* get_unique_DEGs.py
* write_annotations_DEGs.py
* get_longer_transcript.py
* get_OppositeRegulationDegs.py
<p>&nbsp;</p>
  
## Functional analysis
### Software & packages
* Blastx
* gProfiler
* REVIGO
* MapMan
### Scripts
* get_Carabica_annotations.py
* write_annotations_genes.py
* runBlastx.py
* get_unique_homologs.py
* get_gProfiler_input.py
* clean_gProfiler.py
* generate_GMT.py
* get_GO_category.py
* get_specific_GO_terms.py
* get_TairsTFfromiTAK.py
* get_KOgenes.py
* get_MapMan_input.py
<p>&nbsp;</p>
  
## Plotting
### Software & packages
* PCA
* Venn diagrams
* MA and volcano
* Enrichment dotplots
### Scripts
* heatmap.R
* pca_DESeq2.R
* plot_enrichment_gProfiler.R
