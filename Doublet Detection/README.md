We used a two-step procedure to identify doublets. First, we identified doublets in each sample with the re-implementation of the Scrublet algorithm in Pegasus. 
Second, we integrated and clustered all samples and identified clusters significantly enriched for doublets. All nuclei in the enriched clusters were flagged as potential doublets. 


In brief, we integrated the nuclei that passed the quality control, normalized each nuclei to feature counts per 100K counts (FP100K) and log transformed the expression values (log(FP100k + 1)), 
selected highly variable genes, computed the first 30 principal components (PCs), corrected the PCs for batch effects using Harmony, 
and clustered the cells using the Harmony corrected embedding with the Leiden algorithm. 
Then, we tested if each cluster is significantly enriched for doublets using Fisher extract test controlling at a False Discovery Rate of 5%. 
Among the significantly enriched clusters, we selected those with more than 60% of nuclei identified as potential doublets and marked all nuclei in these clusters as doublets.
