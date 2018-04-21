library("gridExtra")
library("ggplot2")

input_dir <- dirname(snakemake@input[[1]])
output_dir <- dirname(snakemake@output[[1]])

vcf_snps <- read.csv(
    snakemake@input[["snps"]],
    header = TRUE,
    na.strings=c("", "NA"),
    sep = "\t")
vcf_indels <- read.csv(
    snakemake@input[["indels"]],
    header = TRUE,
    na.strings=c("","NA"),
    sep = "\t")
vcf <- rbind(vcf_snps, vcf_indels)
vcf$variant <- factor(
    c(rep("SNPs", dim(vcf_snps)[1]),
      rep("Indels", dim(vcf_indels)[1])))

snp_color <- "#A9E2E4"
indel_color <- "#F4CCCA"

QUAL <- ggplot(vcf, aes(x = QUAL, fill = variant)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 100, size = 0.7) +
  xlim(0, 10000)

DP <- ggplot(vcf, aes(x = DP, fill = variant)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = c(10, 100)) +
  xlim(0, 1000)

QD <- ggplot(vcf, aes(x = QD, fill = variant)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 5, size = 0.7)

FS <- ggplot(vcf, aes(x = FS, fill = variant)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = c(60, 200), size = 0.7) +
  ylim(0,0.1) +
  xlim(0, 250)

MQ <- ggplot(vcf, aes(x = MQ, fill = variant)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 35, size = 0.7)

MQRankSum <- ggplot(vcf, aes(x = MQRankSum, fill = variant)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = -12.5, size = 0.7)

SOR <- ggplot(vcf, aes(x=SOR, fill=variant)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = c(4, 10), size = 1, color = c(snp_color, indel_color))

ReadPosRankSum <- ggplot(vcf, aes(x = ReadPosRankSum, fill = variant)) +
  geom_density(alpha = 0.3) +
  geom_vline(
      xintercept = c(-10, 10),
      size = 0.7) +
  xlim(-15, 15)

svg(snakemake@output[[1]], height = 20, width = 15)
theme_set(theme_gray(base_size = 18))
grid.arrange(QUAL, QD, DP, FS, MQ, MQRankSum, SOR, ReadPosRankSum, nrow = 4)
dev.off()
