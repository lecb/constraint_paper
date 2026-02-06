suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})
#KEEP FOR PAPER

GNOMAD_LOF_BGZ <- "gnomad.v2.1.1.all_lofs.txt.bgz"
stopifnot(file.exists(GNOMAD_LOF_BGZ))

OUT_RDS <- "cache/trunc_variants_universe.rds"
TRUNC_CSQ <- c("stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant")

v <- read_tsv(GNOMAD_LOF_BGZ, comment = "#", show_col_types = FALSE, progress = TRUE)

v2 <- v %>%
  transmute(
    chr  = if_else(str_detect(chrom, "^chr"), chrom, paste0("chr", chrom)),
    pos  = as.integer(pos),
    gene = vapply(strsplit(gene_symbols, ","), `[`, character(1), 1),
    consequence = most_severe_consequence
  ) %>%
  filter(!is.na(chr), !is.na(pos), !is.na(gene)) %>%
  filter(gene %in% UNIVERSE$Gene) %>%     # <-- FIX HERE
  filter(consequence %in% TRUNC_CSQ)

saveRDS(v2, OUT_RDS)
message("Wrote: ", OUT_RDS, " | n=", nrow(v2), " | genes=", length(unique(v2$gene)))


suppressPackageStartupMessages({
  library(dplyr)
  library(GenomicRanges)
  library(rtracklayer)
  library(readr)
  library(ggplot2)
})

OUT_CSV <- "variant_qc_summary.csv"
OUT_PNG <- "SuppFig_variant_qc_overlap.png"

v <- readRDS("cache/trunc_variants_universe.rds")

# handle DISCORDANT as tibble vs vector
disc_vec <- if (is.data.frame(DISCORDANT)) DISCORDANT$Gene else DISCORDANT

v <- v %>%
  mutate(group = if_else(gene %in% disc_vec, "discordant", "background"))

stopifnot(nrow(v) > 0)

gr_v <- GRanges(seqnames = v$chr, ranges = IRanges(start = v$pos, end = v$pos))

gr_seg    <- import(BED_SEG)
gr_segdup <- import(BED_SEGDUP)
gr_lowmp  <- import(BED_LOWMP)  # your blacklist BED from Option A

v$is_low_complexity <- countOverlaps(gr_v, gr_seg)    > 0
v$is_segdup         <- countOverlaps(gr_v, gr_segdup) > 0
v$is_problematic    <- countOverlaps(gr_v, gr_lowmp)  > 0

summ_one <- function(flag, feature_name) {
  tmp <- v %>%
    group_by(group) %>%
    summarise(
      feature   = feature_name,
      n         = n(),
      n_overlap = sum(flag, na.rm = TRUE),
      prop      = n_overlap / n,
      .groups   = "drop"
    )
  tab <- table(v$group, flag)
  p <- if (all(dim(tab) == c(2,2))) fisher.test(tab)$p.value else NA_real_
  tmp %>% mutate(fisher_p = p)
}

summ <- bind_rows(
  summ_one(v$is_low_complexity, "low_complexity"),
  summ_one(v$is_segdup,         "segdup"),
  summ_one(v$is_problematic,    "problematic_blacklist")
) %>% arrange(feature, group)

write_csv(summ, OUT_CSV)

p <- ggplot(summ, aes(x = feature, y = prop, shape = group)) +
  geom_point(size = 3) +
  labs(x = NULL, y = "Proportion of truncating variants overlapping feature") +
  theme_bw()

ggsave(OUT_PNG, p, width = 7, height = 4, dpi = 300)

summ_clean <- summ %>%
  filter(feature %in% c("low_complexity", "segdup"))

write_csv(summ_clean, "./tables/SuppTable_variant_qc_lowcomplexity_segdup.csv")

library(dplyr)
library(readr)
library(stringr)
library(tidyr)

v <- read_csv("gnomad_variants_all.csv", show_col_types = FALSE)

# pick a "variant id" column if present; otherwise fall back to chr:pos:ref:alt
if (!("variant_id" %in% names(v))) {
  v <- v %>%
    mutate(variant_id = paste(chrom, pos, ref, alt, sep=":"))
}

tab_splice_free <- v %>%
  mutate(cons_l = tolower(consequence)) %>%
  filter(str_detect(cons_l, "stop_gained|frameshift_variant")) %>%   # splice-free
  mutate(group = if_else(set == "Discordant", "Discordant", "Non-discordant")) %>%
  group_by(gene_symbol, group) %>%
  summarise(
    n_distinct_sf = n_distinct(variant_id),
    n_stop        = n_distinct(variant_id[str_detect(cons_l, "stop_gained")]),
    n_frameshift  = n_distinct(variant_id[str_detect(cons_l, "frameshift_variant")]),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = group,
    values_from = c(n_distinct_sf, n_stop, n_frameshift),
    values_fill = 0
  ) %>%
  arrange(desc(n_distinct_sf_Discordant))

write_csv(tab_splice_free, "./tables/SuppTable_splice_free_stop_frameshift_counts.csv")
tab_splice_free

library(readr)
library(dplyr)

plof  <- read_csv("ancestry_plof_summary_by_gene.csv",  show_col_types = FALSE) %>%
  mutate(analysis = "pLoF_including_splice")

trunc <- read_csv("ancestry_trunc_summary_by_gene.csv", show_col_types = FALSE) %>%
  mutate(analysis = "truncating_only")

# Sanity check: same genes
stopifnot(setequal(plof$Gene, trunc$Gene))

supp_ancestry <- bind_rows(plof, trunc) %>%
  arrange(analysis, Gene)

OUT <- "./tables/SuppTable_ancestry_distribution_top15.csv"
write_csv(supp_ancestry, OUT, na = "")

supp_ancestry
