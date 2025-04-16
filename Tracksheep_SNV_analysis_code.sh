####main codes for gut microbial genetic variations（Tracksheep study）
###for questions, please contact ShuaiYang (shuai_yang8@stu.njmu.edu.cn)


###############################################part01_metagenome_asembly###############################################
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

### step_2 use clean reads to perform de novo assembly
echo "running megahit==================================="
echo "#step_2 single sample assembly...."
mkdir -p ${result_dir}/megahit/assembly/${sample}
mkdir -p ${LOG_DIR}/megahit
time megahit -1 ${result_dir}/kneaddata/pair_reads/${sample}_paired_1.fastq.gz -2 ${result_dir}/kneaddata/pair_reads/${sample}_paired_2.fastq.gz -o ${result_dir}/megahit/assembly/${sample} -f -t ${threads} --out-prefix ${sample} 2> ${LOG_DIR}/megahit/${sample}_assembly.log

### step_3 metagenomic binning
echo "running binning===================================="
echo "#step3.1 build index"
mkdir -p  ${result_dir}/metabat/${sample}/index
mkdir -p ${LOG_DIR}/metabat
time bowtie2-build --threads ${threads} -f  ${result_dir}/megahit/assembly/${sample}/${sample}.contigs.fa ${result_dir}/metabat/${sample}/index/${sample} 2> ${LOG_DIR}/metabat/${sample}_index.log
echo "#step3.2 get sam"
mkdir -p ${result_dir}/metabat/${sample}/bam
time bowtie2 --threads  ${threads} -x ${result_dir}/metabat/${sample}/index/${sample} -1 ${result_dir}/kneaddata/pair_reads/${sample}_paired_1.fastq.gz -2 ${result_dir}/kneaddata/pair_reads/${sample}_paired_2.fastq.gz | samtools sort -@ ${threads} -O bam -o ${result_dir}/metabat/${sample}/bam/${sample}_sort.bam 2> ${LOG_DIR}/metabat/${sample}_bam.log
echo "#step3.3 stat contig depth"
mkdir -p ${result_dir}/metabat/${sample}/depth
time jgi_summarize_bam_contig_depths --outputDepth  ${result_dir}/metabat/${sample}/depth/${sample}_depth.txt ${result_dir}/metabat/${sample}/bam/${sample}_sort.bam 
echo "#step3.4. binning"
mkdir -p ${result_dir}/metabat/${sample}/bin
time metabat -i ${result_dir}/megahit/assembly/${sample}/${sample}.contigs.fa -m 1500 -a ${result_dir}/metabat/${sample}/depth/${sample}_depth.txt -t ${threads} -o ${result_dir}/metabat/${sample}/bin/${sample} 2> ${LOG_DIR}/metabat/${sample}_bin.log
rm ${result_dir}/metabat/${sample}/bam/${sample}_sort.bam

### step_4 genome quality assessing
mkdir -p ${result_dir}/checkM2/${sample}_test
mkdir -p ${LOG_DIR}/checkM2
time checkm2 predict --threads ${threads} --input ${result_dir}/metabat/${sample}/bin -x fa --database_path ${CHECKM2_DB} --output-directory ${result_dir}/checkM2/${sample}_test/  --force 2> ${LOG_DIR}/checkM2/${sample}_test.log

### step_5 Taxonomic classification
gtdbtk classify_wf --genome_dir input  --extension fa --out_dir output --cpus 12 --force --debug --skip_ani_screen






###############################################part02_Microbial_genetic_diversity&phylogenetic_analysis###############################################
### 1. Calculate Kimura distance matrix for species
distmat -sequence "$input_file" -nucmethod 2 -outfile "$output_file"

### 2. Calculate Manhattan distance for each phenotype/metabolite
for (phenotype in colnames(pheno_data)) {
  phenotype_data <- pheno_data[, phenotype, drop = FALSE]
  distance_matrix <- as.matrix(dist(phenotype_data, method = "manhattan"))
  rownames(distance_matrix) <- sample_ID
  colnames(distance_matrix) <- sample_ID
  }

### 3. Association between species genetic distance matrix and phenotype/metabolite variation matrices
genetic_files <- list.files(genetic_path, full.names = TRUE, pattern = "\\.csv$")
phenotypic_files <- list.files(phenotypic_path, full.names = TRUE, pattern = "\\.csv$")
get_upper_tri <- function(mat) {
  return(mat[upper.tri(mat)])
}
results <- data.frame(Species = character(), Phenotype = character(), Correlation = numeric(), P_value = numeric(), stringsAsFactors = FALSE)
for (g_file in genetic_files) {
  genetic_distance <- as.matrix(read.csv(g_file, row.names = 1))
  species_name <- tools::file_path_sans_ext(basename(g_file))
  for (p_file in phenotypic_files) {
    phenotypic_distance <- as.matrix(read.csv(p_file, row.names = 1))
    phenotype_name <- tools::file_path_sans_ext(basename(p_file))
    common_samples <- intersect(rownames(genetic_distance), rownames(phenotypic_distance))
    if (length(common_samples) < 3) {  
      cat("Note: Omitting", species_name, "and", phenotype_name, "due to inadequate sample count\n")
      next
    }
    genetic_distance <- genetic_distance[common_samples, common_samples]
    phenotypic_distance <- phenotypic_distance[common_samples, common_samples]
    genetic_vector <- get_upper_tri(genetic_distance)
    phenotypic_vector <- get_upper_tri(phenotypic_distance)
    cor_test_result <- cor.test(genetic_vector, phenotypic_vector, method = "spearman")
    results <- rbind(results, data.frame(Species = species_name, 
                                         Phenotype = phenotype_name,
                                         Correlation = cor_test_result$estimate, 
                                         P_value = cor_test_result$p.value))
  }
}}

### 4. Phylogenetic tree clustering
library(ape)
library(ggtree)
tree <- hclust(dist(dis))
phylo_tree <- as.phylo(tree)
cluster_assignments <- cutree(tree, k = 2)  # k = 2 indicates partitioning into two clusters






###############################################part03_SNV_identification&subsequent_analysis###############################################
### 1. Filter bins of medium quality based on the criteria according to the results obtained from checkm: completeness ≥50%, contamination <10%.
awk -F '\t' 'NR>1 {if ($6 >= 50 && $7 < 10) print $1}' checkm_filtered_results.txt > medium_quality_bins.txt

### 2. Cluster the medium quality genomes at species level
dRep dereplicate output \
    -g medium_quality_bins/*.fa \
    --S_algorithm fastANI \
    --multiround_primary_clustering \
    --greedy_secondary_clustering \
    -ms 10000 \
    -pa 0.9 \
    -sa 0.95 \
    -nc 0.30 \
    -cm larger \
    -p 24
    
### 3. SNV_identification
### 3.1 make a scaffold-to-bin file
parse_stb.py --reverse -f input/* -o tracksheep.stb
### 3.2 build index
bowtie2-build tracksheep.fa tracksheep --large-index --threads ${threads}
### 3.3 get bam file
bowtie2 -p ${threads} -x tracksheep -1 ${dir}/${sample}_1.fastq.gz -2 ${dir}/${sample}_2.fastq.gz | samtools sort -@ ${threads} -O bam -o - > bam/${sample}.bam
### 3.4 get SNV file
inStrain profile bam/${sample}.bam tracksheep.fa -o instrain/${sample} -p ${threads} -g tracksheep_genes.fna -s tracksheep.stb --database_mode

### 4 Association between SNVs and metabolites
results <- data.frame(
  scaffold_position_ref_base = character(),
  metabolites = character(),
  p_value = numeric(),
  group1_mean = numeric(),
  group2_mean = numeric(),
  group1_median = numeric(),
  group2_median = numeric(),
  stringsAsFactors = FALSE)
for (pos in unique(snp_data$scaffold_position_ref_base)) {
  snp_subset <- snp_data[snp_data$scaffold_position_ref_base == pos, ]
  snp_subset_clean <- snp_subset[, colSums(is.na(snp_subset)) == 0] 
  snp_subset_clean <- snp_subset_clean[, snp_subset_clean[ ,2:ncol(snp_subset_clean)] != "NA"]
  valid_bases <- na.omit(as.vector(snp_subset_clean[, 2:ncol(snp_subset_clean)]))
  unique_bases <- unique(valid_bases)
  snp_subset_clean1 <- snp_subset_clean[, -1]
  group1_samples <- colnames(snp_subset_clean1)[snp_subset_clean1[ ,1:ncol(snp_subset_clean1)] == unique_bases[1]]
  group2_samples <- colnames(snp_subset_clean1)[snp_subset_clean1[ ,1:ncol(snp_subset_clean1)] == unique_bases[2]]
  for (meta_name in colnames(meta_data)) {
    group1_meta <- na.omit(meta_data[group1_samples, meta_name])
    group2_meta <- na.omit(meta_data[group2_samples, meta_name])
    test_result <- wilcox.test(group1_meta, group2_meta)
    results <- rbind(results, data.frame(
      scaffold_position_ref_base = pos,
      metabolites = meta_name,
      p_value = test_result$p.value,
      group1_mean = mean(group1_meta, na.rm = TRUE),
      group2_mean = mean(group2_meta, na.rm = TRUE),
      group1_median = median(group1_meta, na.rm = TRUE),
      group2_median = median(group2_meta, na.rm = TRUE),
      group1_basetype = as.character(unique_bases[1]),
      group2_basetype = as.character(unique_bases[2])
    ))
  }
}






###############################################part04_Association between SNV-related metabolites and phenotypes###############################################
### 1 Association between neurophenotype and metabolites
for (phenotype in phenotype_columns) {
  for (meta_name in metabolites_columns) {
    message("Processing: ", phenotype, " vs ", meta_name)
    valid_data <- merge(raw_phenotype, data, by = "row.names") %>% 
      tibble::column_to_rownames("Row.names") %>%
      .[complete.cases(.[, c(phenotype, meta_name)]), ] %>%
      .[is.finite(.[[microbe]]), ]
    cor_test <- cor.test(valid_data[[phenotype]], valid_data[[meta_name]], 
                        method = "spearman", exact = FALSE)
    lm_model <- lm(as.formula(paste(meta_name, "~", phenotype)), data = valid_data)
    all_p_values <- rbind(all_p_values, data.frame(
      Phenotype = phenotype,
      Metabolites = meta_name,
      P_value = cor_test$p.value,
      Correlation = cor_test$estimate,
      Beta = coef(lm_model)[2],
      stringsAsFactors = FALSE
    ))
  }
}