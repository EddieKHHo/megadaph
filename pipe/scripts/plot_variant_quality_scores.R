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

DP <- ggplot(vcf, aes(x = DP, fill = variant)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = c(10, 6200))

QD <- ggplot(vcf, aes(x = QD, fill = variant)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 2, size = 0.7)

FS <- ggplot(vcf, aes(x = FS, fill = variant)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = c(60, 200), size = 0.7) + ylim(0,0.1)

MQ <- ggplot(vcf, aes(x = MQ, fill = variant)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 40, size = 0.7)

MQRankSum <- ggplot(vcf, aes(x = MQRankSum, fill = variant)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = -20, size = 0.7)

SOR <- ggplot(vcf, aes(x=SOR, fill=variant)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = c(4, 10), size = 1, color = c(snp_color, indel_color))

ReadPosRankSum <- ggplot(vcf, aes(x = ReadPosRankSum, fill = variant)) +
  geom_density(alpha = 0.3) +
  geom_vline(
      xintercept = c(-10, 10, -20, 20),
      size = 1,
      colour = c(snp_color, snp_color, indel_color, indel_color)) +
  xlim(-30, 30)

svg(snakemake@output[[1]], height = 20, width = 15)
theme_set(theme_gray(base_size = 18))
grid.arrange(QD, DP, FS, MQ, MQRankSum, SOR, ReadPosRankSum, nrow = 4)
dev.off()
