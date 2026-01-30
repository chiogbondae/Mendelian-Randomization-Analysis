# MENDELIAN RANDOMIZATION ANALYSIS

library(TwoSampleMR)
library(MendelianRandomization)
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(grid)

tsh_gwas_path <- paste0(base_path, "TSH_GWAS-SNPs.txt")
af_tsh_exposure_path <- paste0(base_path, "AF_TSH-exposure.txt")
tsh_stroke_outcome_path <- paste0(base_path, "TSH_Stroke-SNPs.txt")
af_stroke_outcome_path <- paste0(base_path, "AF_TSH-stroke_outcome.txt")
exclusion_path <- paste0(base_path, "exclusion.txt")

#  TSH exposure data
tsh_gwas_data <- read_tsv(tsh_gwas_path, col_types = cols())
tsh_exposure <- format_data(
  tsh_gwas_data,
  type = "exposure",
  snp_col = "RSID",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1",
  pval_col = "P"
)
tsh_exposure$exposure <- "TSH"
tsh_exposure <- tsh_exposure %>%
  filter(pval.exposure < 5e-8) %>%
  distinct(SNP, .keep_all = TRUE) %>%
  arrange(pval.exposure) %>%
  head(20)
tsh_exposure$f_statistic <- (tsh_exposure$beta.exposure^2) / (tsh_exposure$se.exposure^2)

# stroke outcome data
stroke_raw <- read.table(tsh_stroke_outcome_path, header = TRUE, sep = "", stringsAsFactors = FALSE)
colnames(stroke_raw) <- c("MarkerName", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P.value")
stroke_outcome <- format_data(
  stroke_raw,
  type = "outcome",
  snp_col = "MarkerName",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1",
  pval_col = "P.value"
)
stroke_outcome$outcome <- "Stroke"

# Harmonised TSH–stroke dataset
tsh_stroke_harmonised <- harmonise_data(
  exposure_dat = tsh_exposure,
  outcome_dat = stroke_outcome,
  action = 2
)

# Univariable MR
univariable_results <- mr(tsh_stroke_harmonised, method_list = "mr_ivw")
tsh_ivw <- univariable_results %>% filter(method == "Inverse variance weighted")

#  AF and TSH exposure effects
af_tsh_data <- read_tsv(af_tsh_exposure_path, col_types = cols())
af_data <- af_tsh_data %>%
  filter(Phenotype == "AF") %>%
  select(SNP, beta, se) %>%
  rename(beta_af = beta, se_af = se)
tsh_data <- af_tsh_data %>%
  filter(Phenotype == "TSH") %>%
  select(SNP, beta, se) %>%
  rename(beta_tsh = beta, se_tsh = se)

#  stroke outcome for AF SNPs
stroke_af <- read.table(af_stroke_outcome_path, header = TRUE, sep = "", stringsAsFactors = FALSE)
colnames(stroke_af) <- c("SNP", "Allele1", "Allele2", "Freq1", "beta_stroke", "se_stroke", "P")

# Multivariable MR
excluded_snps <- read_lines(exclusion_path)
mvmr_data <- af_data %>%
  inner_join(tsh_data, by = "SNP") %>%
  inner_join(stroke_af, by = "SNP") %>%
  filter(!SNP %in% excluded_snps)

bx <- cbind(TSH = mvmr_data$beta_tsh, AF  = mvmr_data$beta_af)
bxse <- cbind(TSH = mvmr_data$se_tsh, AF  = mvmr_data$se_af)
mvmr_input <- mr_mvinput(
  bx = bx, bxse = bxse,
  by = mvmr_data$beta_stroke, byse = mvmr_data$se_stroke
)
mvmr_results <- mr_mvivw(mvmr_input)
tsh_direct <- mvmr_results@Estimate["BxTSH"]

# Effect comparison table
tsh_direct_se <- mvmr_results@StdError["BxTSH"]
effect_comparison <- data.frame(
  Model = c("Univariable MR", "Multivariable MR"),
  Beta = c(tsh_ivw$b, tsh_direct),
  SE = c(tsh_ivw$se, tsh_direct_se)
) %>%
  mutate(Lower_CI = Beta - 1.96 * SE, Upper_CI = Beta + 1.96 * SE)

# Conditional F-statistics
t_tsh <- mvmr_data$beta_tsh / mvmr_data$se_tsh
t_af <- mvmr_data$beta_af / mvmr_data$se_af
instrument_correlation <- cor(mvmr_data$beta_tsh, mvmr_data$beta_af)

Q <- cbind(t_tsh, t_af)
Q_squared <- Q^2
f_tsh_sw <- mean(Q_squared[,1]) - (mean(Q[,1]*Q[,2])^2 / mean(Q_squared[,2]))
f_af_sw <- mean(Q_squared[,2]) - (mean(Q[,1]*Q[,2])^2 / mean(Q_squared[,1]))

proper_f_stats <- data.frame(
  Exposure = c("TSH", "AF"),
  Conditional_F_SW = c(f_tsh_sw, f_af_sw),
  Unconditional_F = c(mean(t_tsh^2), mean(t_af^2)),
  Instrument_Correlation = c(instrument_correlation, instrument_correlation)
)
write_csv(proper_f_stats, paste0(base_path, "proper_conditional_f_statistics.csv"))

# Scatter plot
scatter_data <- tsh_stroke_harmonised %>%
  mutate(
    weight = 1/(se.outcome^2),
    f_stat_category = cut(f_statistic, breaks = c(0,10,30,100,Inf),
                          labels = c("Weak (<10)","Moderate (10-30)","Strong (30-100)","Very Strong (>100)"),
                          include.lowest = TRUE)
  )

scatter_plot <- ggplot(scatter_data, aes(x = beta.exposure, y = beta.outcome)) +
  geom_errorbar(aes(ymin = beta.outcome - 1.96*se.outcome, ymax = beta.outcome + 1.96*se.outcome),
                color = "gray70", width = 0) +
  geom_errorbarh(aes(xmin = beta.exposure - 1.96*se.exposure, xmax = beta.exposure + 1.96*se.exposure),
                 color = "gray70", height = 0) +
  geom_point(aes(size = weight, fill = f_stat_category), shape = 21, color = "black", alpha = 0.8) +
  geom_abline(slope = tsh_ivw$b, intercept = 0, color = "#E41A1C", linewidth = 1.2) +
  scale_fill_manual(values = c("Weak (<10)"="#FF9999","Moderate (10-30)"="#FFCC00",
                               "Strong (30-100)"="#66CC99","Very Strong (>100)"="#3399CC")) +
  labs(title="TSH → Stroke", x="Effect of SNP on TSH", y="Effect of SNP on Stroke",
       fill="F-statistic", size="Weight") +
  theme_minimal(base_size=14)

ggsave(paste0(base_path,"tsh_stroke_scatter_enhanced.png"), scatter_plot, width=10, height=8, dpi=300)

# Forest plot
forest_data <- mr_singlesnp(tsh_stroke_harmonised) %>%
  filter(!is.na(b) & !is.na(se)) %>%
  mutate(
    SNP_short = gsub("rs", "", SNP),
    lower_ci = b - 1.96 * se,
    upper_ci = b + 1.96 * se,
    pval_category = cut(p, breaks = c(0,0.001,0.01,0.05,1),
                        labels = c("P < 0.001","P < 0.01","P < 0.05","P ≥ 0.05"),
                        include.lowest = TRUE)
  ) %>%
  arrange(desc(abs(b)))

ivw_estimate <- data.frame(
  SNP_short = "IVW",
  b = tsh_ivw$b,
  se = tsh_ivw$se,
  lower_ci = tsh_ivw$b - 1.96*tsh_ivw$se,
  upper_ci = tsh_ivw$b + 1.96*tsh_ivw$se,
  pval_category = "IVW"
)

plot_data <- bind_rows(forest_data, ivw_estimate) %>%
  mutate(SNP_short = factor(SNP_short, levels = c(forest_data$SNP_short,"IVW")))

forest_plot <- ggplot(plot_data, aes(x=b, y=SNP_short)) +
  geom_vline(xintercept = 0, linetype="dashed", color="gray50") +
  geom_errorbarh(aes(xmin=lower_ci, xmax=upper_ci, color=pval_category), height=0.3, size=0.8) +
  geom_point(aes(fill=pval_category), shape=21, size=3, color="black") +
  scale_color_manual(values=c("P < 0.001"="#7F0000","P < 0.01"="#FF6600",
                              "P < 0.05"="#8A2BE2","P ≥ 0.05"="#228B22","IVW"="#E41A1C")) +
  scale_fill_manual(values=c("P < 0.001"="#7F0000","P < 0.01"="#FF6600",
                             "P < 0.05"="#8A2BE2","P ≥ 0.05"="#228B22","IVW"="#E41A1C")) +
  labs(title="Individual SNP Effects: TSH → Stroke",
       subtitle="Forest plot of Wald ratio estimates",
       x="Causal Effect (β) with 95% CI",
       y="SNP (rsID)",
       color="P-value", fill="P-value") +
  theme_minimal(base_size=14) +
  theme(axis.text.y=element_text(family="mono"), panel.grid.major.y=element_blank())

ggsave(paste0(base_path,"tsh_stroke_forest_robust.png"), forest_plot, width=10, height=12, dpi=300)

# Mediation plot
pie_data <- data.frame(
  Component=c("Direct Effect","AF-Mediated"),
  Percentage=c(abs(tsh_direct/tsh_ivw$b)*100, abs(1 - tsh_direct/tsh_ivw$b)*100)
)

pie_chart <- ggplot(pie_data, aes(x="", y=Percentage, fill=Component)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  geom_text(aes(label=paste0(round(Percentage,1),"%")), position=position_stack(vjust=0.5), size=5, fontface="bold") +
  scale_fill_manual(values=c("Direct Effect"="#4DAF4A","AF-Mediated"="#377EB8")) +
  labs(title="TSH Effect Decomposition via AF", fill="Pathway") +
  theme_minimal(base_size=14) +
  theme(axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(), panel.grid=element_blank())

pathway_diagram <- ggplot() +
  annotate("text", x=1,y=3,label="TSH",size=8,fontface="bold",color="#E41A1C") +
  annotate("text", x=2,y=3,label="AF",size=8,fontface="bold",color="#377EB8") +
  annotate("text", x=3,y=3,label="Stroke",size=8,fontface="bold",color="#4DAF4A") +
  geom_segment(aes(x=1.2,y=3,xend=1.8,yend=3), arrow=arrow(length=unit(0.3,"cm")), linewidth=1.5, color="#E41A1C") +
  geom_segment(aes(x=2.2,y=3,xend=2.8,yend=3), arrow=arrow(length=unit(0.3,"cm")), linewidth=1.5, color="#377EB8") +
  annotate("text", x=1.5,y=3.2,label=paste0("β=",round(tsh_direct,3)), size=4,fontface="bold") +
  annotate("text", x=2.5,y=3.2,label=paste0("β=",round(mvmr_results@Estimate["BxAF"],3)), size=4,fontface="bold") +
  annotate("text", x=2,y=2.6,label=paste0("β=",round(tsh_ivw$b - tsh_direct,3)), size=4,fontface="bold", color="#E41A1C") +
  xlim(0.5,3.5) + ylim(2,3.5) + labs(title="Causal Pathway Diagram") +
  theme_void(base_size=14)

mediation_plot <- ggarrange(pathway_diagram, pie_chart, ncol=2, widths=c(1.2,1))
ggsave(paste0(base_path,"mediation_pathway_diagram.png"), mediation_plot, width=14, height=6, dpi=300)

# 1. Exact P-value for Univariable IVW (TSH -> Stroke)
cat("Univariable IVW P-value:", tsh_ivw$pval, "\n")

# 2. Exact P-values for Multivariable MR (TSH and AF)
cat("Multivariable MR P-values (TSH, AF):", mvmr_results@Pvalue, "\n")

summary_table <- data.frame(
  Exposure = c("TSH (Total)", "TSH (Direct)", "AF (Direct)"),
  Beta = c(tsh_ivw$b, mvmr_results@Estimate["BxTSH"], mvmr_results@Estimate["BxAF"]),
  P_Value = c(tsh_ivw$pval, mvmr_results@Pvalue[1], mvmr_results@Pvalue[2])
)

print(summary_table)

