# lncRNA expression and regulation across different brain regions in primates.
## **Abstract:** <br>
Human and non-human primates have strikingly similar genomes, but they strongly differ in many  brain-based processes (e.g., behaviour and cognition). While the functions of protein-coding genes have been extensively studied, rather little is known about the role of non-coding RNAs such as long noncoding RNAs (lncRNAs). Here, we predicted lncRNAs and analysed their expression pattern across different brain regions of human and non-human primates (chimpanzee, gorilla, and gibbon). Our analysis identified shared orthologous and non-orthologous lncRNAs, showing striking differences in the genomic features. Differential expression analysis of the shared orthologous lncRNAs from humans and chimpanzees revealed distinct expression patterns in subcortical regions (striatum, hippocampus) and neocortical areas while retaining a homogeneous expression in the cerebellum. Co-expression analysis of lncRNAs and protein-coding genes revealed massive proportions of co-expressed pairs in neocortical regions of humans compared to chimpanzees. Network analysis of co-expressed pairs revealed the distinctive role of the hub-acting orthologous lncRNAs in a region- and species-specific manner. Overall, our study provides novel insight into lncRNA driven gene regulatory landscape, neural regulation, brain evolution, and constitutes a resource for primate’s brain lncRNAs. (http://primbrainlnc.bio.uni-mainz.de/)


## **Script Information:**
### Command_for_the_lncRNA_predictions.txt:
- All the commands for tools used for predicting lncRNAs.

### Downstream_data_analysis_of_expression.R
- Rscripts used to get the expression pattern of differentially expressed orthologous lncRNAs from both and either species.

### Human_DEGs_visualization.R
- Provides information about the protein-coding genes differentially expressed across any of the brain regions.
- Extracts the information of variance stabilized read counts, normalized read counts of respective differentially expressed genes.
- Heatmap and boxplot visualization of the expression of differentially expressed orthologs lncRNAs from both and either of the species.

### Co-expression_calculation.R
- Calculates the correlation of expression between Human lncRNAs from each cluster's lncRNA and protein-coding genes.

### Correlation_input_preparation_and_plotting.R
- Prepares the dataframe for getting circular plot showing correlations.

### Plot_statistics_Pos_vs_Neg_co-expressed_Pairs.R
- Plots the statistics of positively and negatively co-expressed lncRNAs and protein-coding genes.


