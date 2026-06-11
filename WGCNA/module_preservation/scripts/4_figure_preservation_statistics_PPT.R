rm(list = ls())
library(ggplot2)
library(patchwork)

out_dir <- "~/Documents/projects/ISF/ISF_fede/code/WGCNA/res/figures"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

set.seed(42)

# ─────────────────────────────────────────────────────────────────────────────
# Simulate small protein sets to illustrate the statistics
# 6 proteins in a module; show IF (preserved) vs SERUM (not preserved)
# ─────────────────────────────────────────────────────────────────────────────
prots <- c("P1","P2","P3","P4","P5","P6")
n     <- 6

# IF: tight positive correlations within module
make_cor <- function(mean_r, noise) {
  m <- matrix(mean_r + rnorm(n*n, 0, noise), n, n)
  m <- (m + t(m)) / 2
  diag(m) <- 1
  m <- pmax(pmin(m, 1), -1)
  dimnames(m) <- list(prots, prots)
  m
}
cor_IF   <- make_cor(0.72, 0.08)   # tight
cor_SERUM <- make_cor(0.05, 0.30)  # no structure

# kIM (intramodular connectivity) = sum of abs correlations with other module members
kim_IF    <- (rowSums(abs(cor_IF))    - 1) / (n - 1)
kim_SERUM <- (rowSums(abs(cor_SERUM)) - 1) / (n - 1)
kim_IF    <- sort(kim_IF,    decreasing = TRUE)
kim_SERUM <- kim_SERUM[names(kim_IF)]   # same protein order as IF

# ─────────────────────────────────────────────────────────────────────────────
# Panel A — DENSITY: correlation matrices (IF vs SERUM)
# ─────────────────────────────────────────────────────────────────────────────
mat_to_long <- function(mat, fluid) {
  df <- as.data.frame(as.table(mat))
  colnames(df) <- c("x","y","r")
  df$fluid <- fluid
  df
}
df_cor <- rbind(mat_to_long(cor_IF, "IF (reference)"),
                mat_to_long(cor_SERUM, "SERUM (test)"))
df_cor$x <- factor(df_cor$x, levels = prots)
df_cor$y <- factor(df_cor$y, levels = rev(prots))
df_cor$r_plot <- ifelse(df_cor$x == df_cor$y, NA, df_cor$r)

pal <- colorRampPalette(c("steelblue","white","firebrick"))(100)

pA <- ggplot(df_cor, aes(x=x, y=y, fill=r_plot)) +
  geom_tile(color="white", linewidth=0.5) +
  facet_wrap(~fluid, ncol=2) +
  scale_fill_gradientn(colors=pal, limits=c(-1,1), na.value="grey90",
                       name="r", guide=guide_colorbar(barheight=5)) +
  labs(title="DENSITY statistics",
       subtitle="Are proteins still highly correlated within the module in SERUM?",
       x=NULL, y=NULL) +
  theme_bw(base_size=10) +
  theme(
    plot.title    = element_text(face="bold", size=12),
    plot.subtitle = element_text(size=8.5, color="grey30"),
    strip.text    = element_text(face="bold", size=10),
    axis.text     = element_text(size=8),
    panel.grid    = element_blank(),
    legend.position = "right"
  ) +
  # Annotation
  annotate("text", x=3.5, y=0.3,
           label="High: module is dense", color="firebrick",
           size=3, hjust=0.5, data=data.frame(fluid="IF (reference)")) +
  annotate("text", x=3.5, y=0.3,
           label="Low: density lost", color="grey50",
           size=3, hjust=0.5, data=data.frame(fluid="SERUM (test)"))

# ─────────────────────────────────────────────────────────────────────────────
# Panel B — CONNECTIVITY: kIM ranking (IF vs SERUM)
# ─────────────────────────────────────────────────────────────────────────────
df_kim <- data.frame(
  protein   = names(kim_IF),
  kIM_IF    = as.numeric(kim_IF),
  kIM_SERUM = as.numeric(kim_SERUM),
  rank_IF   = seq_along(kim_IF)
)
df_kim$rank_SERUM <- rank(-df_kim$kIM_SERUM, ties.method="first")
df_kim$protein    <- factor(df_kim$protein, levels=rev(names(kim_IF)))

# Long format for dots
df_long_k <- rbind(
  data.frame(protein=df_kim$protein, kIM=df_kim$kIM_IF,    fluid="IF (reference)"),
  data.frame(protein=df_kim$protein, kIM=df_kim$kIM_SERUM, fluid="SERUM (test)")
)
df_long_k$fluid <- factor(df_long_k$fluid, levels=c("IF (reference)","SERUM (test)"))

pB <- ggplot(df_long_k, aes(x=fluid, y=kIM, group=protein, color=protein)) +
  geom_line(linewidth=0.8, alpha=0.7) +
  geom_point(size=3) +
  scale_color_brewer(palette="Set2", name="Protein") +
  labs(title="CONNECTIVITY statistics",
       subtitle="Does the hub-protein hierarchy survive in SERUM?\n(cor.kIM: Spearman r between kIM in IF and kIM in SERUM)",
       x=NULL, y="Intramodular connectivity (kIM)") +
  theme_bw(base_size=10) +
  theme(
    plot.title    = element_text(face="bold", size=12),
    plot.subtitle = element_text(size=8.5, color="grey30"),
    axis.text.x   = element_text(face="bold", size=10),
    legend.position = "right"
  )

# ─────────────────────────────────────────────────────────────────────────────
# Panel C — ZSUMMARY: how statistics combine
# ─────────────────────────────────────────────────────────────────────────────
stats_df <- data.frame(
  stat     = c("meanCor\n(avg pairwise r)",
               "propVarExpl\n(variance by PC1)",
               "cor.kIM\n(connectivity rank)",
               "cor.kME\n(eigengene membership)",
               "cor.cor\n(full corr. matrix)"),
  family   = c("Density","Density","Connectivity","Connectivity","Connectivity"),
  Z_example = c(1.1, 0.8, 1.6, 1.3, 1.9)
)
stats_df$stat   <- factor(stats_df$stat, levels=rev(stats_df$stat))
stats_df$family <- factor(stats_df$family, levels=c("Density","Connectivity"))

pC <- ggplot(stats_df, aes(x=Z_example, y=stat, fill=family)) +
  geom_col(width=0.6) +
  geom_vline(xintercept=0, color="black", linewidth=0.3) +
  geom_text(aes(label=paste0("Z = ", Z_example)), hjust=-0.15, size=3.5) +
  scale_fill_manual(values=c("Density"="steelblue3","Connectivity"="darkorange2"),
                    name="Statistic family") +
  scale_x_continuous(limits=c(0, 3.2), expand=c(0,0)) +
  annotate("segment", x=mean(stats_df$Z_example), xend=mean(stats_df$Z_example),
           y=0.3, yend=5.7, linetype="dashed", color="grey40", linewidth=0.7) +
  annotate("text", x=mean(stats_df$Z_example), y=5.9,
           label=paste0("Zsummary = ", round(mean(stats_df$Z_example),1)),
           size=3.5, fontface="bold", color="grey20", hjust=0.5) +
  labs(title="ZSUMMARY: combining all statistics",
       subtitle="Average of individual Z-scores across density + connectivity families\n(each Z = (observed - mean_permuted) / sd_permuted  over 200 permutations)",
       x="Z-score", y=NULL) +
  theme_bw(base_size=10) +
  theme(
    plot.title    = element_text(face="bold", size=12),
    plot.subtitle = element_text(size=8.5, color="grey30"),
    axis.text.y   = element_text(size=9),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom"
  )

# ─────────────────────────────────────────────────────────────────────────────
# Assemble
# ─────────────────────────────────────────────────────────────────────────────
fig <- (pA | pB) / pC +
  plot_annotation(
    title    = "How modulePreservation() evaluates co-expression structure",
    subtitle = "Langfelder et al., PLoS Comp. Biol. 2011 | WGCNA::modulePreservation()",
    theme = theme(
      plot.title    = element_text(face="bold", size=14, hjust=0.5),
      plot.subtitle = element_text(size=9, color="grey40", hjust=0.5)
    )
  )

ggsave(file.path(out_dir, "module_preservation_statistics_PPT.pdf"),
       plot=fig, width=14, height=10)
ggsave(file.path(out_dir, "module_preservation_statistics_PPT.png"),
       plot=fig, width=14, height=10, dpi=300)

message("Saved in: ", out_dir)
