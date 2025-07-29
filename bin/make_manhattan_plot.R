#!/usr/bin/env Rscript

#Usage:
#Rscript make_manhattan_plot.R [assoc_tool(regenie,plink,saige)] [full_path_to_assoc_sumstat] [additional title text] [output name]

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ggrepel)
})

args <- commandArgs(trailingOnly=TRUE)

#Significance Threshold:
sig <- 5*10^-8
#Candidate Threshold:
candi <- 5e-5

if(sig>=candi){stop("Significance threshold must be smaller than candidate threshold!")}


if (length(args) < 1 || !(tolower(args[1]) %in% c("saige", "regenie", "plink", "plink2"))) {
  cat("Error: First argument must be 'saige', 'regenie', or 'plink2'\n")
  quit(status = 1)
}

assoc_tool <- tolower(args[1])

if (length(args) >= 3) {
  added_text <- args[3]
}else{ 
  added_text<- NA_character_
}

if (length(args) >= 4) {
  output_name <- args[4]
}else{
  output_name <- paste0(sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(args[2])))
}
if(!is.na(added_text)){
  title_text <- paste0("Exome, Primary sclerosing cholangits, single-variant test, ",added_text)
}else{
  title_text <- paste0("Exome, Primary sclerosing cholangits, single-variant test")
}

if(assoc_tool=="saige"){
  sumstat <- fread(args[2]) %>%
    rename(CHROM="CHR",ID=MarkerID, GENPOS=POS, ALLELE1=Allele1,ALLELE2=Allele2,"p.value"="p.value") %>%
    mutate(LOG10P=-log10(p.value),
           CHROM=CHROM %>% gsub("X","23",.) %>% as.integer())
  
  print_tool <- "SAIGE"
}else if(assoc_tool=="regenie"){
  sumstat <- fread(args[2])%>%
    mutate(p.value=10**-LOG10P)
  
  print_tool <- "REGENIE"
}else if(assoc_tool %in% c("plink","plink2")){
  sumstat <- fread(args[2]) %>%
    rename(CHROM="#CHROM", GENPOS=POS, ALLELE1=REF,ALLELE2=A1,"p.value"=P) %>%
    mutate(LOG10P=-log10(p.value),
           ID=gsub("unk_","",ID),
           CHROM=CHROM %>% gsub("X","23",.) %>% as.integer()) 
  
  print_tool <- "PLINK2"
}else{
  cat("\nError: First argument must be 'saige', 'regenie', or 'plink2'\n")
  quit(status = 1)
}
cat(paste0("\nReading sumstats ",args[2]))


exome_only_sumstat <- sumstat %>%
  dplyr::filter(!is.na(LOG10P)) %>%
  dplyr::filter(!is.na(p.value)) %>%
  mutate(is_annotate=case_when(p.value<=sig ~ "sig",
                               p.value>=sig & p.value<=candi ~ "candi",
                               p.value>=candi ~ "none")
  ) %>%
  mutate(pointcolor=case_when(is_annotate=="none" ~ "black",
                              is_annotate=="sig" ~ RColorBrewer::brewer.pal(8,"Accent")[c(6)],
                              is_annotate=="candi" ~ RColorBrewer::brewer.pal(8,"Accent")[c(7)])
  )  %>%
  # Variant thinning for faster plotting:
  mutate(
    bp_bin = cut(GENPOS, breaks = 100),                   # Bin by basepair position
    p_bin = cut(LOG10P, breaks = seq(0, max(LOG10P, na.rm=T), by = .1))  # Bin by p-value (log scale)
  ) %>%
  group_by(CHROM, bp_bin, p_bin) %>%
  #slice_sample(n = 10) %>%  # Keep up to 10 variants per bin
  group_modify(~ {
    # Check if this p_bin is "close" to significance
    p_bin_center <- mean((.x$LOG10P))  # Rough estimate of bin's p-values
    if (p_bin_center < candi) {
      # If p-value is low (e.g. < 1e-5), keep all
      .x
    } else {
      # Else, downsample to max 10
      slice_sample(.x, n = 10)
    }
  }) %>% 
  ungroup()


exome_only_cum <- exome_only_sumstat %>% 
  group_by(CHROM) %>% 
  summarise(max_bp = max(GENPOS)) %>% 
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
  select(CHROM, bp_add)


exome_only_sumstat <- exome_only_sumstat  %>% 
  inner_join(exome_only_cum, by = "CHROM")  %>% 
  mutate(bp_cum = GENPOS + bp_add)

exome_only_axis_set <- exome_only_sumstat %>% 
  group_by(CHROM) %>% 
  summarize(center = mean(bp_cum))

exome_only_ylim <- exome_only_sumstat %>% 
  filter(p.value %in% min(p.value)) %>% last() %>%
  mutate(ylim = abs(floor(log10(p.value))) + 2) %>% 
  pull(ylim)




singlemarker_plot <- ggplot(exome_only_sumstat, aes(y=LOG10P, x=bp_cum)) +
  annotate(y=-0.2,  yend=-0.2,  x=exome_only_axis_set[exome_only_axis_set$CHROM == max(exome_only_axis_set$CHROM),]$center-2500000, xend=exome_only_axis_set[exome_only_axis_set$CHROM == min(exome_only_axis_set$CHROM),]$center+2500000, colour="black", lwd=0.75, geom="segment") +
  annotate(y=-0.03, yend=max(exome_only_sumstat$LOG10P) %>% ceiling(), x=min(exome_only_sumstat$bp_cum), xend=min(exome_only_sumstat$bp_cum), colour="black", lwd=0.75, geom="segment") +#30.03,
  annotate(y=-0.03, yend=max(exome_only_sumstat$LOG10P) %>% ceiling(), x=max(exome_only_sumstat$bp_cum), xend=max(exome_only_sumstat$bp_cum), colour="black", lwd=0.75, geom="segment") +#30.03,
  # Candidate threshold line
  annotate(y=-log10(candi)+0.03, yend=-log10(candi), x=min(exome_only_sumstat$bp_cum), xend=max(exome_only_sumstat$bp_cum), colour="lightgrey", lwd=0.75, geom="segment") +
  # Significance threshold line
  annotate(y=-log10(sig)+0.03, yend=-log10(sig), x=min(exome_only_sumstat$bp_cum), xend=max(exome_only_sumstat$bp_cum), colour="darkgrey", lwd=0.75, geom="segment") +#30.03,
  # Plotting non-significant variants
  geom_point(aes(color=as.factor(CHROM)), alpha=0.8, size=1.3) +
  # Color non-significant variants by chr color
  scale_color_manual(values = rep(RColorBrewer::brewer.pal(8,"Accent")[c(1,3)], 23)) +
  
  scale_x_continuous(label = factor(exome_only_axis_set$CHROM, levels = seq(1:23)), breaks= exome_only_axis_set$center, expand = expansion(mult = c(0, .005))) +
  scale_y_continuous(label = c(0,10,20,30,40), breaks= c(0,10,20,30,40), expand = c(0, .05)) + 
  # Plotting significant and candidate Variants
  geom_point(data=exome_only_sumstat%>%filter(is_annotate!="none"), color=exome_only_sumstat%>% filter(is_annotate!="none") %>% pull(pointcolor), alpha=0.8, size=1.3) +
  
  geom_label_repel(aes(label=ID), data=exome_only_sumstat %>%
                     filter(is_annotate!="none") %>%
                     group_by(CHROM) %>%
                     mutate(bp_block=round(bp_cum/5e6)) %>%
                     group_by(CHROM,bp_block) %>%
                     arrange(-LOG10P,.by_group = TRUE ) %>%
                     slice(1),
                   size=3,
                   fontface = 'bold', segment.color = 'grey80',
                   #nudge_x = -100,
                   nudge_y = 15,
                   direction = "both",
                   #hjust = 0.5,
                   segment.size=0.1,
                   #force=1,
                   max.overlaps = 30)+
  labs(title=title_text, x="Chromosome", y=expression(paste(-log[10], "(", italic(P), "-value)")) )+
  theme_bw() +
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
# singlemarker_plot


ggsave(filename=paste0(output_name,"_manhattan.jpg"),plot=singlemarker_plot, height=7, width = 13)
#ggsave(filename=paste0(output_name,"_manhattan.svg"),plot=singlemarker_plot, height=7, width = 13)
cat(paste0("\nManhattan plot saved as ",paste0(output_name,"_manhattan.jpg")))
ps_assoc <- sumstat %>% filter(!is.na(p.value)) %>% pull(p.value)

inflation <- function(ps) {
  chisq <- qchisq(1 - ps, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  lambda
}

gg_qqplot <- function(ps, ci = 0.95, thin=T, thin.obs.places=2) {
  
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  
  if(thin){
    df <- df %>% mutate(observed=observed %>% round(thin.obs.places),
                        expected=expected %>% round(thin.obs.places)) %>% unique()
  }
  
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed), shape = 1, size = 2.5) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po)
}

qqplotfull <-  gg_qqplot(ps_assoc) +
  theme_bw(base_size = 24) +
  annotate(
    geom = "text",
    x = -Inf,
    y = Inf,
    hjust = -0.15,
    vjust = 1 + 0.15 * 3,
    label = paste0(sprintf("lambda = %.3f", inflation(ps_assoc)),"; N Variants: ",length(ps_assoc)),
    size = 8
  ) +
  theme(
    axis.ticks = element_line(linewidth = 0.5),
    panel.grid = element_blank()
    # panel.grid = element_line(size = 0.5, color = "grey80")
  )

ggsave(filename=paste0(output_name,"_qqplot.jpg"),plot=qqplotfull, height=7, width = 13)
cat(paste0("\nQQplot saved as ",output_name,"_qqplot.jpg\n"))