#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

##To run: system2("Rscript", "02_stats_analysis.R")
# KEEP FOR PAPER

if (file.exists("relpos_run_meta.rds")) {
  meta <- readRDS("relpos_run_meta.rds")
  message("[CACHE META] created_at=", meta$created_at,
          " | gw_sample_n=", meta$gw_sample_n,
          " | gw_relpos_n=", meta$gw_relpos_n)
} else {
  message("[CACHE META] relpos_run_meta.rds not found")
}
disc_rel  <- readRDS("disc_rel.rds")
gw_relpos <- readRDS("gw_relpos.rds")

# Clean rel_pos in each (do NOT assume upstream filtering)
disc_rel2 <- disc_rel  %>% filter(is.finite(rel_pos), rel_pos >= 0, rel_pos <= 1)
gw_rel2   <- gw_relpos %>% filter(is.finite(rel_pos), rel_pos >= 0, rel_pos <= 1)

combo <- bind_rows(disc_rel2, gw_rel2)

# Wilcoxon
combo$source <- factor(combo$source, levels = c("Discordant", "Genome-wide"))
wilc <- wilcox.test(rel_pos ~ source, data = combo, exact = FALSE)

# Bin helper
bin_relpos <- function(r) case_when(
  is.na(r) ~ NA_character_,
  r < 1/3 ~ "early",
  r < 2/3 ~ "middle",
  TRUE ~ "late"
)

disc_bins <- table(bin_relpos(disc_rel2$rel_pos))
gw_bins   <- table(bin_relpos(gw_rel2$rel_pos))

bin_tab <- rbind(
  Discordant   = c(early = as.integer(disc_bins["early"]),
                   middle = as.integer(disc_bins["middle"]),
                   late = as.integer(disc_bins["late"])),
  `Genome-wide` = c(early = as.integer(gw_bins["early"]),
                    middle = as.integer(gw_bins["middle"]),
                    late = as.integer(gw_bins["late"]))
)
bin_tab[is.na(bin_tab)] <- 0L

bin_test <- if (any(bin_tab < 5)) fisher.test(bin_tab) else chisq.test(bin_tab)

late_mat <- matrix(
  c(bin_tab["Discordant","late"],
    bin_tab["Discordant","early"] + bin_tab["Discordant","middle"],
    bin_tab["Genome-wide","late"],
    bin_tab["Genome-wide","early"] + bin_tab["Genome-wide","middle"]),
  nrow = 2, byrow = TRUE,
  dimnames = list(source = c("Discordant","Genome-wide"), bin = c("late","non_late"))
)
late_fisher <- fisher.test(late_mat)

# Summary table
optA_summary <- tibble(
  source = c("Discordant","Genome-wide"),
  n_relpos = c(nrow(disc_rel2), nrow(gw_rel2)),
  relpos_median = c(median(disc_rel2$rel_pos), median(gw_rel2$rel_pos)),
  relpos_IQR = c(IQR(disc_rel2$rel_pos), IQR(gw_rel2$rel_pos)),
  early = c(bin_tab["Discordant","early"], bin_tab["Genome-wide","early"]),
  middle = c(bin_tab["Discordant","middle"], bin_tab["Genome-wide","middle"]),
  late = c(bin_tab["Discordant","late"], bin_tab["Genome-wide","late"])
) %>%
  mutate(prop_late = late / pmax(1, early + middle + late))

write_csv(optA_summary, "trunc_position_domain_nmd_OPTION_A_summary.csv")

# Tests txt
test_lines <- c(
  "---- OPTION A: Discordant vs genome-wide truncations ----",
  paste0("Discordant usable rel_pos n = ", nrow(disc_rel2)),
  paste0("Genome-wide usable rel_pos n = ", nrow(gw_rel2)),
  "",
  "Wilcoxon test (rel_pos Discordant vs Genome-wide):",
  paste0("W = ", unname(wilc$statistic), ", p = ", signif(wilc$p.value, 4)),
  "",
  "Early/middle/late bin table:",
  capture.output(print(bin_tab)),
  "",
  paste0(bin_test$method, ": p = ", signif(bin_test$p.value, 4)),
  "",
  "Late vs non-late (Discordant vs Genome-wide):",
  paste0(
    "OR = ", signif(unname(late_fisher$estimate), 4),
    " (95% CI ", signif(late_fisher$conf.int[1], 4), "–", signif(late_fisher$conf.int[2], 4), ")",
    ", Fisher p = ", signif(late_fisher$p.value, 4)
  )
)

writeLines(test_lines, "trunc_position_domain_nmd_OPTION_A_tests.txt")

message("[DONE] 02_stats_analysis.R wrote summary + tests")
