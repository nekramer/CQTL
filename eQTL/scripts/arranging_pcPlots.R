library(ggplot2)
library(gridExtra)


loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


age <- loadRData("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/pca/ALL_Age.RData")
cause_of_death <- loadRData("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/pca/ALL_CauseOfDeath.RData")
condition <- loadRData("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/pca/ALL_Condition.RData")
dna_reagent_batch <- loadRData("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/pca/ALL_DNAReagentBatch.RData")
fragment_batch <- loadRData("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/pca/ALL_Fragment_Batch.RData")
genotyping_batch <- loadRData("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/pca/ALL_GenotypingBatch.RData")
oa_grade_avg <- loadRData("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/pca/ALL_OAGradeAvg.RData")
race <- loadRData("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/pca/ALL_Race.RData")
rna_extractionkit_batch <- loadRData("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/pca/ALL_RNAextractionKitBatch.RData")
sequencing_batch <- loadRData("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/pca/ALL_SequencingBatch.RData")
sex <- loadRData("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/pca/ALL_Sex.RData")


age <- age + theme(plot.subtitle = element_blank())
cause_of_death <- cause_of_death + theme(plot.subtitle = element_blank())
condition <- condition + theme(plot.subtitle = element_blank())
dna_reagent_batch <- dna_reagent_batch + theme(plot.subtitle = element_blank())
fragment_batch <- fragment_batch + theme(plot.subtitle = element_blank())
genotyping_batch <- genotyping_batch + theme(plot.subtitle = element_blank())
oa_grade_avg <- oa_grade_avg + theme(plot.subtitle = element_blank())
race <- race + theme(plot.subtitle = element_blank())
rna_extractionkit_batch <- rna_extractionkit_batch + theme(plot.subtitle = element_blank())
sequencing_batch <- sequencing_batch + theme(plot.subtitle = element_blank())
sex <- sex + theme(plot.subtitle = element_blank())


pdf(file = "/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/ALL_pcs.pdf",
    width = 14, height = 30)

grid.arrange(age, cause_of_death, condition, dna_reagent_batch, fragment_batch, genotyping_batch,
          oa_grade_avg, race, rna_extractionkit_batch, sequencing_batch, sex, ncol = 2,
          padding =2)
dev.off()


CTL_age <- loadRData("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/pca/CTL_Age.RData")
FNF_age <- loadRData("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/pca/FNF_Age.RData")
CTL_cause_of_death <- loadRData("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/pca/CTL_CauseOfDeath.RData")
FNF_cause_of_death <- loadRData("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/pca/FNF_CauseOfDeath.RData")
CTL_dna_reagent_batch <- loadRData("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/pca/CTL_DNAReagentBatch.RData")
FNF_dna_reagent_batch <- loadRData("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/pca/FNF_DNAReagentBatch.RData")

CTL_fragment_batch <- loadRData("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/pca/CTL_Fragment_Batch.RData")
FNF_fragment_batch <- loadRData("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/pca/FNF_Fragment_Batch.RData")

CTL_genotyping_batch <- loadRData("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/pca/CTL_GenotypingBatch.RData")
FNF_genotyping_batch <- loadRData("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/pca/FNF_GenotypingBatch.RData")
CTL_oa_grade_avg <- loadRData("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/pca/CTL_OAGradeAvg.RData")
FNF_oa_grade_avg <- loadRData("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/pca/FNF_OAGradeAvg.RData")
CTL_race <- loadRData("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/pca/CTL_Race.RData")
FNF_race <- loadRData("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/pca/FNF_Race.RData")
CTL_rna_extractionkit_batch <- loadRData("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/pca/CTL_RNAextractionKitBatch.RData")
FNF_rna_extractionkit_batch <- loadRData("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/pca/FNF_RNAextractionKitBatch.RData")
CTL_sequencing_batch <- loadRData("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/pca/CTL_SequencingBatch.RData")
FNF_sequencing_batch <- loadRData("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/pca/FNF_SequencingBatch.RData")
CTL_sex <- loadRData("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/pca/CTL_Sex.RData")
FNF_sex <- loadRData("/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/output/pca/FNF_Sex.RData")


pdf(file = "/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch/Condition_PCs.pdf",
    width = 14, height = 45)
grid.arrange(CTL_age, FNF_age,
             CTL_cause_of_death, FNF_cause_of_death, 
             CTL_dna_reagent_batch, FNF_dna_reagent_batch,
             CTL_fragment_batch, FNF_fragment_batch,
             CTL_genotyping_batch, FNF_genotyping_batch,
             CTL_oa_grade_avg, FNF_oa_grade_avg,
             CTL_race, FNF_race, 
             CTL_rna_extractionkit_batch, FNF_rna_extractionkit_batch, 
             CTL_sequencing_batch, FNF_sequencing_batch, 
             CTL_sex, FNF_sex, 
             ncol = 2,
             padding =2)
dev.off()
