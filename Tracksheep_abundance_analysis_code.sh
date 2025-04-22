####main codes for 微生物驱动性情变异（Tracksheep study）
###for questions, please contact ShuaiYang (shuai_yang8@stu.njmu.edu.cn)


###############################################part01_metagenomic raw data processing###############################################
### step_1 remove host-contaminated sequences
export PATH=~/soft/Trimmomatic-0.33:$PATH
sample=F1
threads=20
data_dir=~/tracksheep/base_clean_data
result_dir=~/results
KNEADDATA_DB_HUMAN_GENOME=/~/database/kneaddata/hg37
METAPHLAN_DATABASE=~/database/metaphlan4
LOG_DIR=~/results/log
TRIMMOMATIC_PATH=~/Trimmomatic-0.39
CHECKM2_DB=~/database/checkm2_db/uniref100.KO.1.dmnd

echo "# step_1.1 kneadata...."
mkdir -p ${result_dir}/kneaddata/${sample}/filtering_data/
mkdir -p ${LOG_DIR}/kneaddata
kneaddata --input ${data_dir}/${sample}_1.fq.gz --input ${data_dir}/${sample}_2.fq.gz -t ${threads} --trimmomatic ${TRIMMOMATIC_PATH} --output-prefix ${sample} --reference-db ${KNEADDATA_DB_HUMAN_GENOME} --serial --reorder --output ${result_dir}/kneaddata/${sample}/filtering_data/   --log ${LOG_DIR}/kneaddata/${sample}.log

echo "# step_1.2 move files...."
mkdir -p ${result_dir}/kneaddata/pair_reads
mv ${result_dir}/kneaddata/${sample}/filtering_data/${sample}_paired_1.fastq ${result_dir}/kneaddata/${sample}/filtering_data/${sample}_paired_2.fastq ${result_dir}/kneaddata/pair_reads
rm -r ${result_dir}/kneaddata/${sample}/filtering_data/
gzip ${result_dir}/kneaddata/pair_reads/${sample}_paired_1.fastq
gzip ${result_dir}/kneaddata/pair_reads/${sample}_paired_2.fastq

### step_2 taxonomy abundance profile
metaphlan --bowtie2db ${METAPHLAN_DATABASE}  ${result_dir}/kneaddata/pair_reads/${sample}_paired_1.fastq.gz,${result_dir}/kneaddata/pair_reads/${sample}_paired_2.fastq.gz --index mpa_v202403_GMR --input_type fastq --nproc ${threads} --output_file ${result_dir}/metaphlan/single_sample/${sample}.profile -s ${result_dir}/metaphlan/single_sample/${sample}.sam.bz2 --bowtie2out ${result_dir}/metaphlan/single_sample/${sample}.bowtie2_out.bz2  --unclassified_estimation 2> ${LOG_DIR}/metaphlan/${sample}.log 

### step_3 functional pathways profile
mkdir -p ${LOG_DIR}/humann
mkdir -p ${result_dir}/humann/${sample}
time humann -i ${result_dir}/kneaddata/pair_reads/${sample}_paired_1.fastq.gz \
      -i ${result_dir}/kneaddata/pair_reads/${sample}_paired_2.fastq.gz \
      --taxonomic-profile ${result_dir}/metaphlan/single_sample/${sample}.profile \
      --threads ${threads} \
      --output ${result_dir}/humann/${sample}/ \
      2> ${LOG_DIR}/humann/${sample}.log

### step_4 use clean reads to perform de novo assembly
echo "running megahit==================================="
echo "#step_4 single sample assembly...."
mkdir -p ${result_dir}/megahit/assembly/${sample}
mkdir -p ${LOG_DIR}/megahit
time megahit -1 ${result_dir}/kneaddata/pair_reads/${sample}_paired_1.fastq.gz -2 ${result_dir}/kneaddata/pair_reads/${sample}_paired_2.fastq.gz -o ${result_dir}/megahit/assembly/${sample} -f -t ${threads} --out-prefix ${sample} 2> ${LOG_DIR}/megahit/${sample}_assembly.log

### step_5 metagenomic binning
echo "running binning===================================="
echo "#step5.1 build index"
mkdir -p  ${result_dir}/metabat/${sample}/index
mkdir -p ${LOG_DIR}/metabat
time bowtie2-build --threads ${threads} -f  ${result_dir}/megahit/assembly/${sample}/${sample}.contigs.fa ${result_dir}/metabat/${sample}/index/${sample} 2> ${LOG_DIR}/metabat/${sample}_index.log
echo "#step5.2 get sam"
mkdir -p ${result_dir}/metabat/${sample}/bam
time bowtie2 --threads  ${threads} -x ${result_dir}/metabat/${sample}/index/${sample} -1 ${result_dir}/kneaddata/pair_reads/${sample}_paired_1.fastq.gz -2 ${result_dir}/kneaddata/pair_reads/${sample}_paired_2.fastq.gz | samtools sort -@ ${threads} -O bam -o ${result_dir}/metabat/${sample}/bam/${sample}_sort.bam 2> ${LOG_DIR}/metabat/${sample}_bam.log
echo "#step5.3 stat contig depth"
mkdir -p ${result_dir}/metabat/${sample}/depth
time jgi_summarize_bam_contig_depths --outputDepth  ${result_dir}/metabat/${sample}/depth/${sample}_depth.txt ${result_dir}/metabat/${sample}/bam/${sample}_sort.bam 
echo "#step5.4 binning"
mkdir -p ${result_dir}/metabat/${sample}/bin
time metabat -i ${result_dir}/megahit/assembly/${sample}/${sample}.contigs.fa -m 1500 -a ${result_dir}/metabat/${sample}/depth/${sample}_depth.txt -t ${threads} -o ${result_dir}/metabat/${sample}/bin/${sample} 2> ${LOG_DIR}/metabat/${sample}_bin.log
rm ${result_dir}/metabat/${sample}/bam/${sample}_sort.bam

### step_6 genome quality assessing
mkdir -p ${result_dir}/checkM2/${sample}
mkdir -p ${LOG_DIR}/checkM2
time checkm2 predict --threads ${threads} --input ${result_dir}/metabat/${sample}/bin -x fa --database_path ${CHECKM2_DB} --output-directory ${result_dir}/checkM2/${sample}/  --force 2> ${LOG_DIR}/checkM2/${sample}.log

### step_7 get high quality bins for prokka and GTDB-tk
mkdir -p ${result_dir}/metabat/${sample}/high_quality_bins
mkdir -p ${LOG_DIR}/metabat/MAG
awk -F"\t" 'NR>1 && $2 >= 50 && $3 < 10 && ($2 - 5 * $3) > 50 {print $1}' ${result_dir}/checkM2/${sample}/quality_report.tsv | while read bin; do
mv ${result_dir}/metabat/${sample}/bin/$bin.fa ${result_dir}/metabat/${sample}/high_quality_bins/
done
if [ $? -ne 0 ]; then
2> ${LOG_DIR}/metabat/MAG/${sample}.log

### step_8 GTDB-tk Taxonomic classification
gtdbtk classify_wf --genome_dir input  --extension fa --out_dir output --cpus 12 --force --debug --skip_ani_screen






###############################################part02_Microbiome explains temperament variation among merino sheep###############################################
### step_1 Calculate the inter-individual temperament variation distance
library(vegan)
library(tibble)
data <- read.csv("197phenotypes.csv", row.names = 1)
data_matrix <- as.matrix(data)
data_matrix_no_na <- na.omit(data_matrix)
manhattan_distance <- vegdist(data_matrix_no_na, method = "manhattan")
manhattan_matrix <- as.matrix(manhattan_distance)
write.csv(manhattan_matrix, file = "temperament_manhattan_distance_matrix.csv", row.names = TRUE)

### step_2 Calculate the explanatory power of the microbiome on temperament variation
abundance_data <- read.csv("197filtered_relative_species_abundance.csv", row.names = 1)
results_df <- data.frame(
  Variable = character(),
  R_squared = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)
for (i in 1:ncol(abundance_data)) {
  set.seed(66)
  single_abundance <- abundance_data[, i]
  bacteria_name <- colnames(abundance_data)[i]
  adonis_result <- adonis2(manhattan_matrix ~ single_abundance, permutations = 999)
  Df <- adonis_result$Df[1]
  SumOfSqs <- adonis_result$SumOfSqs[1]
  R2 <- adonis_result$R2[1]
  F_stat <- adonis_result$`F`[1]
  p_value <- adonis_result$`Pr(>F)`[1]
  results_df <- rbind(results_df, data.frame(
    species = bacteria_name,
    Df = Df,
    SumOfSqs = SumOfSqs,
    R2 = R2,
    F_stat = F_stat,
    p_value = p_value,
    stringsAsFactors = FALSE
  ))
}
write.csv(results_df, file = "adonis_results_all_bacteria66.csv", row.names = FALSE)
### step_3 Remove multicollinearity
species_data <- read.csv("significant_species_abundance_data.csv", row.names = 1)
cor_matrix <- cor(species_data, use = "pairwise.complete.obs")
write.csv(cor_matrix, file = "significant_bacteria_cor_matrix.csv", row.names = TRUE)
cutoff <- 0.9
highly_correlated_pairs <- which(abs(cor_matrix) > cutoff, arr.ind = TRUE)
highly_correlated_pairs <- highly_correlated_pairs[highly_correlated_pairs[, 1] != highly_correlated_pairs[, 2], ]
result <- data.frame(
  Species_1 = rownames(cor_matrix)[highly_correlated_pairs[, 1]],
  Species_2 = colnames(cor_matrix)[highly_correlated_pairs[, 2]],
  Correlation = cor_matrix[highly_correlated_pairs]
)






###############################################part03_Association between rumen microbiome and temperament phenotypes###############################################
### step_1 Association between species abundance and temperament
raw_phenotype <- read.csv("197phenotypes.csv", row.names = 1)
data <- read.csv("197species_abundance.csv", row.names = 1)
phenotype_columns <- colnames(raw_phenotype)[!colnames(raw_phenotype) %in% c("sample_ID")]
microbe_columns <- colnames(data)[!colnames(data) %in% c("sample_ID")]
all_p_values <- data.frame(Phenotype = character(), Microbe = character(), P_value = numeric(), stringsAsFactors = FALSE)
analyze_ass <- function(pheno_df, micro_df, phenotype, microbe, output_dir) {
  merged_data <- merge(pheno_df, micro_df, by = "row.names", all = TRUE)
  rownames(merged_data) <- merged_data$Row.names
  merged_data$Row.names <- NULL
  valid_data <- merged_data[complete.cases(merged_data[, c(phenotype, microbe)]), c(phenotype, microbe)]
  valid_data <- valid_data[is.finite(valid_data[[microbe]]), ]
  correlation_test <- cor.test(valid_data[[phenotype]], valid_data[[microbe]], method = "spearman")
  lm_model <- lm(as.formula(paste(microbe, "~", phenotype)), data = valid_data)
  beta_coefficient <- coef(lm_model)[2]
  p_value <- correlation_test$p.value
  correlation_coefficient <- correlation_test$estimate
  global_p_values <- get("all_p_values", envir = .GlobalEnv)
  global_p_values <- rbind(global_p_values, data.frame(Phenotype = phenotype, Microbe = microbe, P_value = p_value, Correlation = correlation_coefficient, Beta = beta_coefficient, stringsAsFactors = FALSE))
  assign("all_p_values", global_p_values, envir = .GlobalEnv)
}
for (phenotype in phenotype_columns) {
  phenotype_output_dir <- file.path("clr_corr_phenotype_2", phenotype)
  if (!dir.exists(phenotype_output_dir)) {
    dir.create(phenotype_output_dir, recursive = TRUE)
  }
  for (microbe in microbe_columns) {
    analyze_ass(raw_phenotype, data, phenotype, microbe, phenotype_output_dir)
  }    
}
global_p_values <- get("all_p_values", envir = .GlobalEnv)
global_p_values$Adjusted_P_value <- p.adjust(global_p_values$P_value, method = "BH")
write.csv(global_p_values, file = "all_adj_p_values.csv", row.names = FALSE)

### step_2 Association between species and pathway
raw_pathway <- read.csv("197pathway_abundance.csv", row.names = 1)
data <- read.csv("197species_abundance.csv", row.names = 1)
pathway_columns <- colnames(raw_pathway)[!colnames(raw_pathway) %in% c("sample_ID")]
microbe_columns <- colnames(data)[!colnames(data) %in% c("sample_ID")]
all_p_values <- data.frame(pathway = character(), Microbe = character(), P_value = numeric(), stringsAsFactors = FALSE)
analyze_ass <- function(pheno_df, micro_df, pathway, microbe, output_dir) {
  merged_data <- merge(pheno_df, micro_df, by = "row.names", all = TRUE)
  rownames(merged_data) <- merged_data$Row.names
  merged_data$Row.names <- NULL
  valid_data <- merged_data[complete.cases(merged_data[, c(pathway, microbe)]), c(pathway, microbe)]
  valid_data <- valid_data[is.finite(valid_data[[microbe]]), ]
  correlation_test <- cor.test(valid_data[[pathway]], valid_data[[microbe]], method = "spearman")
  lm_model <- lm(as.formula(paste(microbe, "~", pathway)), data = valid_data)
  beta_coefficient <- coef(lm_model)[2]
  p_value <- correlation_test$p.value
  correlation_coefficient <- correlation_test$estimate
  global_p_values <- get("all_p_values", envir = .GlobalEnv)
  global_p_values <- rbind(global_p_values, data.frame(pathway = pathway, Microbe = microbe, P_value = p_value, Correlation = correlation_coefficient, Beta = beta_coefficient, stringsAsFactors = FALSE))
  assign("all_p_values", global_p_values, envir = .GlobalEnv)
}
for (pathway in pathway_columns) {
  pathway_output_dir <- file.path("clr_corr_pathway_2", pathway)
  if (!dir.exists(pathway_output_dir)) {
    dir.create(pathway_output_dir, recursive = TRUE)
  }
  for (microbe in microbe_columns) {
    analyze_ass(raw_pathway, data, pathway, microbe, pathway_output_dir)
  }    
}
global_p_values <- get("all_p_values", envir = .GlobalEnv)
global_p_values$Adjusted_P_value <- p.adjust(global_p_values$P_value, method = "BH")
write.csv(global_p_values, file = "all_adj_p_values.csv", row.names = FALSE)






#########################part04_Mediation linkages among rumen MetaCyc pathway, baseline plasma metabolites and temperament phenotypes#########################
### step_1 mediation analysis: microbial impacts on host temperament through metabolites
all=cbind(pheno,microbe,diet)
ass=read.delim("ass_path_pheno.txt",header = T,sep = "\t")
ass=ass[which(ass$Pval<0.05),]
colnames(ass)[3:4]=paste("di_pheno",colnames(ass)[3:4],sep = "_")
tmp=read.delim("ass_path_meta.txt",header = T,sep = "\t")
tmp=tmp[which(tmp$Pval<0.05),]
colnames(tmp)[3:4]=paste("di_microbe",colnames(tmp)[3:4],sep = "_")
tmp1=read.delim("ass_meta_pheno.txt",header = T,sep = "\t")
tmp1=tmp1[which(tmp1$Pval<0.05),]
colnames(tmp1)[3:5]=paste("pheno_path",colnames(tmp1)[3:5],sep = "_")
tmp1$name=paste(tmp1$path,ass$pheno,sep = "_with_")
ass=merge(ass[,1:4],tmp[,1:4],by="di",all = F)
ass$name=paste(ass$microbe,ass$pheno,sep = "_with_")
ass=merge(ass,tmp1[,3:6],by="name",all = F)
ass=ass[,c(2,3,6,4:5,7:11)]
##mediation
ass$Pval_mediate=NA
ass$Pval_direct=NA
ass$med_prop=NA
ass$Pval_mediate_inverse=NA
ass$Pval_direct_inverse=NA
library(mediation)
for(i in 1:nrow(ass)){
  data=all[,c(as.character(ass$di[i]),as.character(ass$pheno[i]),as.character(ass$microbe[i]))]
  data=na.omit(data)
  colnames(data)=c("X","Y","M")
  #MetaCyc pathway influence temperament through metabolites
    model.m=lm(M~X,data)
    model.y=lm(Y~X+M,data)
    summary=summary(mediate(model.m ,model.y,treat = "X", mediator = "M",boot = T,sims = 1000))
    ass$Pval_mediate[i]=summary$d.avg.p
    ass$Pval_direct[i]=summary$z.avg.p
    ass$med_prop[i]=summary$n0
  #inverse mediate
  colnames(data)=c("Y","X","M")
  model.m=lm(M~X,data)
  model.y=lm(Y~X+M,data)
  summary=summary(mediate(model.m ,model.y,treat = "X", mediator = "M",boot = T,sims = 1000))
  ass$Pval_mediate_inverse[i]=summary$d.avg.p
  ass$Pval_direct_inverse[i]=summary$z.avg.p
}