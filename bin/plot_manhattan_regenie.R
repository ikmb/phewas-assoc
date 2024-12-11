library(data.table)
library(tidyverse)
library(ggrepel)

args = commandArgs(trailingOnly=TRUE)

regenie_sum <- fread(args[1]) %>% mutate(p.value=10**-LOG10P)
output_name <- args[2]

data_cum <- regenie_sum %>% 
  group_by(CHROM) %>% 
  summarise(max_bp = max(GENPOS)) %>% 
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
  select(CHROM, bp_add)


regenie_sum <- regenie_sum  %>% 
  inner_join(data_cum, by = "CHROM")  %>% 
  mutate(bp_cum = GENPOS + bp_add)

axis_set <- regenie_sum %>% 
  group_by(CHROM) %>% 
  summarize(center = mean(bp_cum))

ylim <- regenie_sum %>% 
  filter(p.value %in% min(p.value)) %>% last() %>%
  mutate(ylim = abs(floor(log10(p.value))) + 2) %>% 
  pull(ylim)

#sig <- 0.05 / nrow(regenie_sum)
sig <- 5*10^-8

regenie_sum <- regenie_sum %>% 
  mutate(is_annotate=ifelse(LOG10P>=-log10(sig), "yes","no"))


singlemarker_plot <- ggplot(regenie_sum, aes(y=LOG10P, x=bp_cum)) +
  annotate(y=-0.2,  yend=-0.2,  x=axis_set[axis_set$CHROM == max(axis_set$CHROM),]$center-2500000, xend=axis_set[axis_set$CHROM == min(axis_set$CHROM),]$center+2500000, colour="black", lwd=0.75, geom="segment") +
  annotate(y=-0.03, yend=max(regenie_sum$LOG10P), x=min(regenie_sum$bp_cum), xend=min(regenie_sum$bp_cum), colour="black", lwd=0.75, geom="segment") +#30.03,
  annotate(y=-0.03, yend=max(regenie_sum$LOG10P), x=max(regenie_sum$bp_cum), xend=max(regenie_sum$bp_cum), colour="black", lwd=0.75, geom="segment") +#30.03,
  annotate(y=-log10(sig)+0.03, yend=-log10(sig), x=min(regenie_sum$bp_cum), xend=max(regenie_sum$bp_cum), colour="grey", lwd=0.75, geom="segment") +#30.03,
  geom_point(aes(color=as.factor(CHROM)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(RColorBrewer::brewer.pal(8,"Accent")[c(1,3)], 23)) +
  scale_x_continuous(label = factor(axis_set$CHROM, levels = seq(1:23)), breaks= axis_set$center, expand = expand_scale(mult = c(0, .005))) +
  scale_y_continuous(label = c(0,10,20,30), breaks= c(0,10,20,30), expand = c(0,.5)) + 
  geom_point(data=regenie_sum[(is_annotate=="yes"),], color=RColorBrewer::brewer.pal(8,"Accent")[6], alpha=0.8, size=1.3) +
  geom_label_repel(aes(label=ID), data=regenie_sum[(LOG10P >= 5),], 
                   size=3, 
                   fontface = 'bold', segment.color = 'grey80', 
                   nudge_x = -100, nudge_y = 0, direction = "x", hjust = 0.5, 
                   segment.size=0.1, force=1, max.overlaps = 10)+
  labs(title="", x="Chromosome", y=expression(paste(-log[10], "(", italic(P), "-value)")) )+
  theme_bw() +
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(filename=paste0(output_name,"_manhattan.jpg"),plot=singlemarker_plot, height=7, width = 13)