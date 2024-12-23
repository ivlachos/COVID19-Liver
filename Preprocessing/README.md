The raw sequencing reads were demultiplexed using Cell Ranger mkfastq (10x Genomics). 
We trimmed the reads from the BIDMC liver samples for polyA tails and the template switching oligo 5’- AAGCAGTGGTATCAACGCAGAGTACATrGrGrG -3’ with cutadapt v.2.7.
The reads were aligned to generate the count matrix using Cell Ranger count (10x Genomics) on the Terra Cloud platform (https://app.terra.bio/) with the cellranger_workflow in Cumulus (https://github.com/klarman-cell-observatory/cumulus). 
