A pathway score summarizes the expression of a set of functionally related genes . 
A Gene Ontology set of 989 GO Biological Process terms was used to create a curated 
selection of pathways capturing liver parenchymal and non-parenchymal cellular functions 
and pathways (Supplementary Table S9). Building on the methodology described in PCxN, 
we used a rank-based approach to define the pathway scores, 
where the pathway score is the sum of the adjusted ranks of the genes in the pathway annotation 
scaled by the square root of the number of genes in the pathway. 
First, the ranks based on the UMI counts are calculated per gene for each nucleus solving ties by selecting the minimum. 
Then, we scale and center the ranks across each nucleus. 
In order to account for the effect of rank sparsity for each gene we split the scaled and centered ranks by their sign 
(positive or negative) and regress out with a linear model the effect of the number of genes detected and the log of the total number of UMIs. 
Finally, we use the removeBatchEffect function from limma to adjust the pathway scores for batch effects. 
The same approach was used to estimate a score for the curated signatures described by Sánchez-Taltavull et al. (proliferating Kupffer cells), 
and by Niethamer et al. (influenza-injury signature).
