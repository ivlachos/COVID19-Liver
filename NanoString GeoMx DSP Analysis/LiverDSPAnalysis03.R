rm(list=ls())
options(stringsAsFactors = FALSE)

library(data.table)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(limma)

library(parallel)
library(tidyverse)
library(rafalib)
library(pheatmap)
library(patchwork)
library(googledrive)
library(RColorBrewer)
library(DT)
Sys.setenv("RSTUDIO_PANDOC" = "/data/resources/tools/pandoc/")


# ====> Liver DSP Data <====
# ==== Raw Data ====
# **SegmentProperties.txt - The current annotations that we have available for each segment within each ROI. 
# This will be updated once the BIDMC team provides complete annotations for the slides which will then be used for analysis
SegmentProperties = fread(
  "/data/work/Projects/BrScRNAseq/NanoString/Liver_normalization/20210614/preprocessing_data/BIDMC_liver_only_20210614_SegmentProperties.txt",
  data.table = FALSE
)

rownames(SegmentProperties) = SegmentProperties$Sample_ID
SegmentProperties$Patient = factor(SegmentProperties$Patient)
SegmentProperties$lobule = factor(SegmentProperties$lobule)


# **TargetProperties.txt - Annotation of the probes & pools from which they are derived - 1:1 with the rows of the BioProbeCountMatrix
TargetProperties = fread("/data/work/Projects/BrScRNAseq/NanoString/Liver_normalization/20210614/preprocessing_data/BIDMC_liver_only_20210614_TargetProperties.txt")

# Count Matrices (targets/probes as rows, samples as columns; column 1 is target name):

# *BIDMC_Covid_BioProbeCountMatrix.txt - uncollapsed probe count matrix. 
# This is effectively the deduplicated count matrix after processing FASTQ files. 
# No QC has been performed on this beyond standard FASTQ processing. Pipeline information can be provided.
BioProbeCountMatrix = fread(
  "/data/work/Projects/BrScRNAseq/NanoString/Liver_normalization/20210614/preprocessing_data/BIDMC_liver_only_20210614_BioProbeCountMatrix.txt",
  data.table = FALSE
)


# *BIDMC_Covid_TargetCountMatrix.txt - collapsed target count matrix. 
# This is effectively the raw data file if you do not want to work at a per probe level 
# (only impacts negative controls and COVID spike-ins which have more than 1 probe per target). 
# Targets are reported as the geometric mean of the probes after QC'ing probes for outliers which 
# may add noise to the target expression in certain AOIs. More information can be provided about outlier calling methods.
TargetCountMatrix = fread(
  "/data/work/Projects/BrScRNAseq/NanoString/M-317 BIDMC_COVID WTA/preliminary_data_20200714/BIDMC_Covid_TargetCountMatrix.txt",
  data.table = FALSE
)

# **NegNorm* - Count matrix normalized against the geometric mean of the negative control probes. 
NegNormCountMatrix = fread(
  "/data/work/Projects/BrScRNAseq/NanoString/Liver_normalization/20210614/preprocessing_data/BIDMC_liver_only_20210614_NegNorm_TargetCountMatrix.txt",
  data.table = FALSE
)

SampleIDs = grep("DSP",colnames(BioProbeCountMatrix),value=T)


# ==== WTA Negative Probes ====
NegProbeIDs = BioProbeCountMatrix$ProbeDisplayName[BioProbeCountMatrix$TargetName == "Neg Probe"]
TargetsWTA = unique(TargetProperties$TargetName[TargetProperties$Pooling == "1"])

NegProbeCounts = BioProbeCountMatrix[BioProbeCountMatrix$TargetName == "Neg Probe",SampleIDs] 
rownames(NegProbeCounts) = NegProbeIDs

TotalCountsWTA = colSums(BioProbeCountMatrix[BioProbeCountMatrix$TargetName %in% TargetsWTA,SampleIDs])

# ===== COVID Negative Probes =====
COVIDNegProbeIDs = BioProbeCountMatrix$ProbeDisplayName[BioProbeCountMatrix$TargetName == "SARS-CoV-2 Neg"]
TargetsCOVID = unique(TargetProperties$TargetName[TargetProperties$Pooling == "2"])

COVIDNegProbeCounts = BioProbeCountMatrix[BioProbeCountMatrix$TargetName == "SARS-CoV-2 Neg",SampleIDs] 
rownames(COVIDNegProbeCounts) = COVIDNegProbeIDs

TotalCountsCOVID = colSums(BioProbeCountMatrix[BioProbeCountMatrix$TargetName %in% TargetsCOVID,SampleIDs])

# ==== Liver annotation ===== 
zone_annot = readRDS(
  "/data/work/Projects/BrScRNAseq/NanoString/NK/Liver/ZoneAnnotation.RDS"
)
rownames(zone_annot) = gsub(".","-",rownames(zone_annot),fixed=T)

liver_annot = readRDS(
  "/data/work/Projects/BrScRNAseq/NanoString/LiverAnnotation_CD45Extended.RDS"
)
rownames(liver_annot) = gsub(".","-",rownames(liver_annot),fixed=T)

table(liver_annot[,c("SlideAnnotation","Patient")])


liver_annot$log_aoi_size = log(liver_annot$aoi_size)

# ===== Select ROIs ====
liver_annot = liver_annot[!grepl("CD45",liver_annot$SlideAnnotation),]
liver_annot$nuclei_counts = SegmentProperties[rownames(liver_annot),"nuclei_counts"]
liver_annot$log_nuclei_counts = log(as.numeric(liver_annot$nuclei_counts))



# ==== WTA Raw Counts ====
TargetsWTA = setdiff(TargetsWTA,"Neg Probe")
RawCountsWTA = BioProbeCountMatrix[BioProbeCountMatrix$TargetName %in% TargetsWTA,rownames(liver_annot)]
rownames(RawCountsWTA) = BioProbeCountMatrix$TargetName[BioProbeCountMatrix$TargetName %in% TargetsWTA]

# ==== COVID Raw Counts ====
TargetsCOVID = setdiff(TargetsCOVID,"SARS-CoV-2 Neg")
covid_genome = c("ORF1ab","ORF1ab_REV","S")

RawCountsCOVID = TargetCountMatrix[TargetCountMatrix$TargetName %in% TargetsCOVID,rownames(liver_annot)]
rownames(RawCountsCOVID) = TargetCountMatrix$TargetName[TargetCountMatrix$TargetName %in% TargetsCOVID]



# ==== COVID Score =====
covid_test = lapply(
  1:nrow(liver_annot),
  function(ic){
    
    # subset data for current ROI
    tmpdat = c(
      RawCountsCOVID[covid_genome,ic],
      COVIDNegProbeCounts[,rownames(liver_annot)[ic]]
    )
    names(tmpdat) = c(covid_genome,rownames(COVIDNegProbeCounts))
    
    tmplabels = rep("NegProbe",length(tmpdat))
    tmplabels[names(tmpdat) %in% covid_genome] = "COVID"
    
    # estimate gene set score
    zrnk_obs = c(scale(rank(tmpdat)))
    covid_score = sqrt(length(covid_genome))*mean(zrnk_obs[tmplabels == "COVID"])
    
    # permute labels and estimate null distribution of gene set score
    bootLst = mclapply(
      1:1E5,
      function(ic){
        bootlabels = sample(x=tmplabels,size=length(tmplabels))
        
        sqrt(length(covid_genome))*mean(zrnk_obs[bootlabels == "COVID"])
      },
      mc.cores = 8
    )
    bootLst = unlist(bootLst)
    
    return(
      c(
        covid_score = covid_score,
        scovid_score = (covid_score - mean(bootLst))/sd(bootLst),
        p.boot = mean(bootLst > covid_score)
      )
    )
    
    
  }
)

names(covid_test) = rownames(liver_annot)
covid_test = do.call(rbind,covid_test)
covid_test = as.data.frame(covid_test)


# adjust for ROI size and number of nuclei
currentLst = split(x=rownames(liver_annot),f=list(liver_annot$Patient,liver_annot$region))

scovid_scoreLst = lapply(
  currentLst,
  function(current_rois){
    removeBatchEffect(
      x= t(covid_test[current_rois,"scovid_score",drop=F]),
      covariates = liver_annot[current_rois,c("log_nuclei_counts","log_aoi_size")]
    )
  }
)

adjregion_scovid_score = unlist(scovid_scoreLst)
names(adjregion_scovid_score) = unlist(currentLst)
adjregion_scovid_score = adjregion_scovid_score[rownames(liver_annot)]

# ===== Test: Patients ======
# fit linear model
patient_lfit01 = lmFit(
  object = matrix(adjregion_scovid_score,nrow=1),
  design = model.matrix( ~ 0 + Patient, liver_annot)
)

AllPatients = as.character(1:4)
PatientComparisonLst = lapply(
  AllPatients,
  function(sck){
    eval(
      parse(
        text = paste0(
          "patientContrast01 = makeContrasts(\n",
          "Patient",sck,"vsAll = ",
          "Patient",sck," - (",
          paste(
            paste0(
              "Patient",setdiff(AllPatients,sck)
            ),
            collapse = " + "
          ),
          ")/",length(AllPatients) -1,
          ",\nlevels = colnames(patient_lfit01$coefficients)\n)"
        )
      )
    )
    
    eval(
      parse(
        text = paste0(
          "patientContrast02 = makeContrasts(\n",
          paste(
            sapply(
              setdiff(AllPatients,sck),
              function(x){
                paste0(
                  "Patient",sck,"vs",x," = ",
                  "Patient",sck," - ","Patient",x
                )
              }
            ),
            collapse=",\n"
          ),
          ",\nlevels = colnames(patient_lfit01$coefficients)\n)"
        )
      )
    )
    
    patientContrast_lfit01 <- contrasts.fit(patient_lfit01,patientContrast01)
    patientContrast_lfit01 <- eBayes(patientContrast_lfit01)
    patientContrast_toptab01 <- topTable(patientContrast_lfit01)
    
    patientContrast_lfit02 <- contrasts.fit(patient_lfit01,patientContrast02)
    patientContrast_lfit02 <- eBayes(patientContrast_lfit02)
    patientContrast_toptab02 <- topTable(patientContrast_lfit02)
    
    pairwise_decide = decideTests(patientContrast_lfit02)
    pairwise_toptab = do.call(
      rbind,lapply(
        colnames(patientContrast_lfit02$coefficients), 
        function(coefname){
          topTable(patientContrast_lfit02,coef=coefname)
        }
      )
    )
    rownames(pairwise_toptab) = colnames(patientContrast_lfit02$coefficients)
    pairwise_toptab$decideTest = c(pairwise_decide[,rownames(pairwise_toptab)])
    
    return(
      list(
        ContrastVsAll = patientContrast_toptab01,
        ContrastPairwise = patientContrast_toptab02,
        ContrastPairwise2 = pairwise_toptab
      )
    )
    
    
  }
)

names(PatientComparisonLst) = AllPatients



# ===== Regions ====
AllRegions = unique(liver_annot$region)


RegionComparisonLst = lapply(
  AllPatients,
  function(sck){
    current_rois = rownames(liver_annot)[liver_annot$Patient == sck]
    # region lfit
    region_lfit01 = lmFit(
      object = matrix(adjregion_scovid_score[current_rois],nrow=1),
      design = model.matrix( ~ 0 + region, liver_annot[current_rois,])
    )
    # loop thru regions
    region_res = lapply(
      AllRegions,
      function(current_region){
        eval(
          parse(
            text = paste0(
              "regionContrast01 = makeContrasts(\n",
              "region",current_region,"vsAll = ",
              "region",current_region," - (",
              paste(
                paste0(
                  "region",setdiff(AllRegions,current_region)
                ),
                collapse = " + "
              ),
              ")/",length(AllRegions) -1,
              ",\nlevels = colnames(region_lfit01$coefficients)\n)"
            )
          )
        )
        
        eval(
          parse(
            text = paste0(
              "regionContrast02 = makeContrasts(\n",
              paste(
                sapply(
                  setdiff(AllRegions,current_region),
                  function(x){
                    paste0(
                      "region",current_region,"vs",x," = ",
                      "region",current_region," - ","region",x
                    )
                  }
                ),
                collapse=",\n"
              ),
              ",\nlevels = colnames(region_lfit01$coefficients)\n)"
            )
          )
        )
        
        
        regionContrast_lfit01 <- contrasts.fit(region_lfit01,regionContrast01)
        regionContrast_lfit01 <- eBayes(regionContrast_lfit01)
        regionContrast_toptab01 <- topTable(regionContrast_lfit01)
        
        regionContrast_lfit02 <- contrasts.fit(region_lfit01,regionContrast02)
        regionContrast_lfit02 <- eBayes(regionContrast_lfit02)
        regionContrast_toptab02 <- topTable(regionContrast_lfit02)
        
        pairwise_decide = decideTests(regionContrast_lfit02)
        pairwise_toptab = do.call(
          rbind,lapply(
            colnames(regionContrast_lfit02$coefficients), 
            function(coefname){
              topTable(regionContrast_lfit02,coef=coefname)
            }
          )
        )
        rownames(pairwise_toptab) = colnames(regionContrast_lfit02$coefficients)
        pairwise_toptab$decideTest = c(pairwise_decide[,rownames(pairwise_toptab)])
        
        
        
        return(
          list(
            ContrastVsAll = regionContrast_toptab01,
            ContrastPairwise = regionContrast_toptab02,
            ContrastPairwise2 = pairwise_toptab
          )
        )
        
        
      }
    )
    names(region_res) = AllRegions
    return(region_res)
    
  }
)
names(RegionComparisonLst) = AllPatients


# ===== Save Results ====
res_covid = data.frame(
  liver_annot[,c("Patient","region")],
  COVIDScore = c(scale(adjregion_scovid_score[rownames(liver_annot)]))
)
res_covid$Patient = paste0("L",res_covid$Patient)


saveRDS(
  res_covid,
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/",
    "LiverDSP/COVIDEnrichment/",
    "LiverDSPCovidEnrichmentScores.RDS"
  )
)

saveRDS(
  RegionComparisonLst,
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/",
    "LiverDSP/COVIDEnrichment/",
    "RegionComparisonLst.RDS"
  )
)

saveRDS(
  PatientComparisonLst,
  paste0(
    "/data/work/Projects/BrScRNAseq/output/BIBroadLiver/SelectLiversNoDoublets/",
    "LiverDSP/COVIDEnrichment/",
    "PatientComparisonLst.RDS"
  )
)


# ==== Results Overview ====
# Patient Comparison
for(ip in AllPatients){
  cat(paste0("Patient ",ip,"\n"))
  print(PatientComparisonLst[[ip]]$ContrastVsAll )# %>% select(-c(adj.P.Val,B)))
}

for(ip in AllPatients){
  cat(paste0("Patient ",ip,"\n"))
  print(PatientComparisonLst[[ip]]$ContrastPairwise ) # %>% select(-c(adj.P.Val,`F`)))
}

for(ip in AllPatients){
  print(PatientComparisonLst[[ip]]$ContrastPairwise2 ) # %>% select(-c(adj.P.Val,`F`)))
  cat("\n\n\n")
}

# Region Comparison (Patient 1)
RegionComparisonLst[["1"]]




# ==== Patient Plots =====
patient_bp = ggplot(
  res_covid,
  aes(x=Patient,y=COVIDScore)
) + 
  geom_boxplot() +
  theme_cowplot() +
  theme(
    axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,size=11),
    legend.position = "none"
  ) +
  ylab("") + xlab("") 



# ==== Patient1: Region Results =====
liver_annot$x_coords = as.numeric(liver_annot$x_coords )
liver_annot$y_coords = as.numeric(liver_annot$y_coords )

# group ROIs
liver_annot$SpatialGroup = zone_annot[rownames(liver_annot),c("SpatialGroup")]
liver_annot["DSP-1012340031703-H-B06","SpatialGroup"] = "G1_2"
liver_annot["DSP-1012340031703-H-B07","SpatialGroup"] = "G1_1"
liver_annot["DSP-1012340031703-H-C09","SpatialGroup"] = "G1_3"
liver_annot["DSP-1012340031703-H-C05","SpatialGroup"] = "G1_4"

# plot using spatial coordinates
plotdf = data.frame(
  liver_annot,
  Exprs = adjregion_scovid_score[rownames(liver_annot)]
) %>% filter(Patient == "1")

plotdf$Exprs = c(scale(plotdf$Exprs))

ggplot(
  data = plotdf,
  aes(x=x_coords,y=y_coords,shape=region,fill=Exprs)
)+ 
  geom_point(size=3) +
  scale_shape_manual(
    values = c(
      "Portal" = 22,"Zone1" = 23,"Zone2" = 24,"Zone3" = 25
    )
  ) +
  scale_fill_gradient2(low="blue",mid="lightgray",high="red",name="COVID Score") +
  scale_x_reverse() +
  scale_y_reverse() +
  xlab("") + ylab("") + 
  #facet_wrap(vars(SpatialGroup), scales = "free")
  guides(shape = guide_legend(title="Region"))






