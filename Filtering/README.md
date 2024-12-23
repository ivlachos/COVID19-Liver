We filtered out nuclei with fewer than 400 UMIs, 200 genes, or greater than 20% of UMIs mapped to mitochondrial genes. 
Furthermore, we discarded samples with less than 100 nuclei. We retained all nuclei that pass the quality metrics described above. 
Subsequently, snRNA-seq data from individual samples were combined into a single expression matrix and analyzed using Seurat v.3.2.3. 
The UMI counts for each nuclei were divided by the total counts for that nuclei, and multiplied by a scale factor of 10,000. 
Then, values are log-transformed using log1p resulting in log(1+10,000*UMIs/Total UMIs) for each nucleus. 
Subsequently, highly variable genes were identified using Seurat’s FindVariableFeatures function. 
Then, data dimensionality was reduced to the top 15 principal components by PCA using the top 2000 highly variable genes. 
The lower dimensional embedding was then corrected for technical noise using each sample as a separate batch with Harmony. 
Neighbors were computed using the Harmony-corrected embedding. 
The UMAP and Leiden clusters were computed using the resulting nearest neighbor graph. 
