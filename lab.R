#' ---
#' title: "Lab 3, homework"
#' author: "Urszula Sawczuk, Michal Wietecki,"
#' ---

setwd("/Users/michalwietecki/Desktop/pw/9_sem_barca/bsg/bsgLAB3")

install.packages("data.table", type="binary")
library(data.table)

install.packages("HardyWeinberg", type = "binary")
library(genetics)
library(HardyWeinberg)

# Linkage Disequilibrium

# 1.

data <- read.table("./data/FOXP2/FOXP2.dat", header = TRUE)

n_ind <- nrow(data)
n_snps <- ncol(data) - 1 # minus the id column

table(unlist(data[, -1]), useNA = "always")

# there are 104 individuals, 543 SNPs and NO data is missing 

# 2.

g1 <- genotype(data$rs34684677)
g2 <- genotype(data$rs2894715)

ld_res <- LD(g1, g2)
print(ld_res)

# p-val = 5.77645e-06 thus the association of
# alleles of those two SNPs is significant

# 3. 
D <- ld_res$D
p1 <- summary(g1)$allele.freq[, 2]
p2 <- summary(g2)$allele.freq[, 2]

# get names of the alleles
A <- names(p1)[1]
a <- names(p1)[2]
B <- names(p2)[1]
b <- names(p2)[2]

p_A <- p1[1]
p_a <- p1[2]
p_B <- p2[1]
p_b <- p2[2]

hap_freqs <- c(
  p_A * p_B + D,
  p_A * p_b - D,
  p_a * p_B - D,
  p_a * p_b + D
)

names(hap_freqs) <- c(
  paste0(A, B),
  paste0(A, b),
  paste0(a, B),
  paste0(a, b)
)

print(hap_freqs)
cat("most common haplotype is: ", names(which.max(hap_freqs)))

# 4. 

bim <- read.table("./data/FOXP2/FOXP2.bim", header = FALSE, stringsAsFactors = FALSE)
n_snps <- nrow(bim)

counts_mat <- matrix(0, nrow = n_snps, ncol = 3)

for (i in 1:n_snps) {
  snp_id <- trimws(bim[i, 2])
  
  #reference
  ref_alleles <- sort(c(trimws(bim[i, 5]), trimws(bim[i, 6])))
  homo1 <- paste0(ref_alleles[1], ref_alleles[1])
  hetero <- paste0(ref_alleles[1], ref_alleles[2])
  homo2 <- paste0(ref_alleles[2], ref_alleles[2])
  
  # for the
  clean_geno <- gsub("/", "", as.character(data[[snp_id]]))
  
  # alphabetical sorting
  clean_geno <- sapply(clean_geno, function(x) {
    paste(sort(unlist(strsplit(x, ""))), collapse = "")
  })
  
  # genotypes counting
  counts_mat[i, 1] <- sum(clean_geno == homo1)
  counts_mat[i, 2] <- sum(clean_geno == hetero)
  counts_mat[i, 3] <- sum(clean_geno == homo2)
}

valid_rows <- rowSums(counts_mat) > 0
counts_mat_clean <- counts_mat[valid_rows, , drop = FALSE]
length(counts_mat_clean)/3

hw_results <- HWChisqMat(counts_mat_clean, cc = 0)
p_values <- hw_results$pvalvec

rejected_count <- sum(p_values < 0.05, na.rm = TRUE)
expected_chance <- 0.05 * nrow(counts_mat_clean)

cat("count of rejected HWE: ", rejected_count, "\n")
cat("expected from chance: ", expected_chance, "\n")


