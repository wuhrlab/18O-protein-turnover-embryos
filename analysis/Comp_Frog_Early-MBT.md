Comparing frog 2C and gastrulation data
================
Edward Cruz
2026-02-19

``` r
rm(list=ls(all=T))
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(tidyr)
library(ggplot2)
library(patchwork)

knitr::opts_chunk$set(fig.path = "figures/Comp_Frog_Early-MBT/")
```

``` r
Frog_T5_all.min <- c(161, 221, 281, 401, 521, 641,
                     881, 1605, 1961, 2918, 4405) - 161
Frog_T6_all.min <- c(150, 210, 270, 390, 510, 635,
                     930, 1591, 1950, 2910, 4394) - 150
Frog_Early.min <- c(Frog_T5_all.min + Frog_T6_all.min)/2

Frog_MBT.min <- c(1472, 1502, 1535, 1592, 1714, 1832,
                  2412, 2934, 3270, 3749, 4260, 4411) - 1472

colorblind_palette <- c("#CC79A7", "#D55E00", "#0072B2", "#F0E442",
                        "#009E73", "#56B4E9", "#E69F00", "#999999")
```

``` r
early_theo_peps <- read.csv("Files/Fits/XLA_O18-Early-TS_Theo-Peptide-Decay.csv")
mbt_theo_peps <-read.csv("Files/Fits/XLA_O18-MBT-TS_Theo-Peptide-Decay.csv")
```

``` r
shared_theo_peps <- merge(early_theo_peps, mbt_theo_peps, by="Peptide",
                          suffixes=c("_2C", "_MBT"))

head(shared_theo_peps)
```

    ##                     Peptide Pep_k1_2C  Pep_k2_2C Pep_k1_MBT Pep_k2_MBT
    ## 1          AAAADGMEQMEMDESR         1 0.09458704          1 0.22193346
    ## 2       AAAAGNEAASLFLATHGAK         1 0.06110211          1 0.13301161
    ## 3             AAAANQGPDVLNK         1 0.03039152          1 0.06693866
    ## 4 AAAASAAEAGITSTTGDDSDEALLK         1 0.11717971          1 0.22609981
    ## 5                   AAAASVR         1 0.03409022          1 0.06933514
    ## 6                 AAADIAENK         1 0.02711900          1 0.06110995

``` r
median_2cell_residues <- read.csv("Files/Fits/AA_examples/Median_2cell_residue-prob.csv")
median_gast_residues <- read.csv("Files/Fits/AA_examples/Median_Gast_residue-prob.csv")


ex.p <- ggplot() +
  geom_point(data=median_2cell_residues, 
             aes(x=Time, y=Fit_values), color="#0072B2", size=2,
             shape=1) +
  geom_point(data=median_gast_residues, 
             aes(x=Time, y=Fit_values), color="#CC79A7", size=2,
             shape=1) +

  # # geom_line(data=ex.fit_data[ex.fit_data$Type=="2C_Min",],
  # #           aes(x=Time, y=Fit_values), color="#0072B2", size=1.5,
  # #           alpha=0.25, linetype="dashed") +
  # geom_line(data=ex.fit_data[ex.fit_data$Type=="2C_Median",],
  #           aes(x=Time, y=Fit_values), color="#0072B2", size=1.5) +  
  # 
  # # geom_line(data=ex.fit_data[ex.fit_data$Type=="MBT_Min",],
  # #           aes(x=Time, y=Fit_values), color="#CC79A7", size=1.5,
  # #           alpha=0.25, linetype="dashed") +
  # geom_line(data=ex.fit_data[ex.fit_data$Type=="MBT_Median",],
  #           aes(x=Time, y=Fit_values), color="#CC79A7", size=1.5) +  
  theme_bw() +
  labs(x="Hours after Labeling", y="Prob. of light peptide") +
  theme_bw() +
  theme(axis.text=element_text(size=21,colour="black"),
        plot.title = element_text(size=24, hjust=0.5),
        axis.title = element_text(size=22),
        panel.border = element_rect(linewidth=2),
        panel.grid.major = element_line(color = "grey70", linewidth = 0.25),
        panel.grid.minor = element_line(color = "grey70", linewidth = 0.25),
        aspect.ratio = 1,
        legend.position = "none") +
  coord_cartesian(ylim=c(0,1), xlim=c(0,6))


# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=4, height=4, res=300)

ex.p
```

![](figures/Comp_Frog_Early-MBT/Comparing%202C%20and%20gastrulation%20residue%20probabilties-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
median_2cell_residues[median_2cell_residues$Fit_values<0.5,][1,]
```

    ##    Time Fit_values
    ## 19  0.3  0.4833439

``` r
median_gast_residues[median_gast_residues$Fit_values<0.5,][1,]
```

    ##    Time Fit_values
    ## 10 0.15  0.4979259

``` r
early_fits.df <- read.csv("Files/Fits/Dyn-Model/XLA_O18_T5-6_FinalFits_ALL.csv")
mbt_fits.df <- read.csv("Files/Fits/Dyn-Model/XLA_O18_T9-10_FinalFits_ALL.csv")

length(unique(c(early_fits.df$Protein_ID, mbt_fits.df$Protein_ID)))
```

    ## [1] 9606

``` r
early_fits.df <- early_fits.df[early_fits.df$Estimate=="interval",]
mbt_fits.df <- mbt_fits.df[mbt_fits.df$Estimate=="interval",]

early_min_kd <- log(2)/(4244*3)
early_max_hl <- log(2)/early_min_kd/60

mbt_min_kd <- log(2)/(2939*3)
mbt_max_hl <- log(2)/mbt_min_kd/60

combined_ts <- merge(early_fits.df[c("Protein_ID", "Human_Gene", "kD",
                                     "T5_kD", "T6_kD",
                                     "HL_Hrs", "mClass", "Number_Pep_Fits")],
                     mbt_fits.df[c("Protein_ID", "kD",
                                   "T9_kD", "T10_kD",                                   
                                   "HL_Hrs", "mClass", "Number_Pep_Fits")],
                     by="Protein_ID", suffixes=c("_Early", "_MBT"))

combined_ts[(combined_ts$mClass_Early=="Deg") &
              (combined_ts$kD_Early < early_min_kd), "mClass_Early"] <- "Deg-Long"
combined_ts[(combined_ts$mClass_MBT=="Deg") &
              (combined_ts$kD_MBT < mbt_min_kd), "mClass_MBT"] <- "Deg-Long"

combined_ts["Comb_kD"] <- (combined_ts$kD_Early + combined_ts$kD_MBT)/2
combined_ts["Combined_Class"] <- paste0(combined_ts$mClass_Early, "_", combined_ts$mClass_MBT)
# combined_ts <- combined_ts[!combined_ts$Comb_kD==1e-6,]

combined_ts <- combined_ts %>%
  arrange(desc(kD_Early))
```

``` r
ss_bp.all <- combined_ts
ss_bp.all <- ss_bp.all[ss_bp.all$mClass_Early=="Deg" | ss_bp.all$mClass_Early=="Deg-Long" |
                         ss_bp.all$mClass_MBT=="Deg" | ss_bp.all$mClass_MBT=="Deg-Long",]

ss_bp.all["Log10_Early"] <- abs(ss_bp.all$kD_Early)
ss_bp.all["Log10_Late"] <- abs(ss_bp.all$kD_MBT)
ss_bp.all["Log10_Early"] <- log10(ss_bp.all$Log10_Early)
ss_bp.all["Log10_Late"] <- log10(ss_bp.all$Log10_Late)
ss_bp.all[ss_bp.all$Log10_Early<(-6), "Log10_Early"] <- (-6)
ss_bp.all[ss_bp.all$Log10_Late<(-6), "Log10_Late"] <- (-6)


p1 <- ggplot() +
  geom_point(data=ss_bp.all[ss_bp.all$Human_Gene %in%
                              c("CCNB1", "KIF22", "GMNN"),],
             aes(x=Log10_Late, y=Log10_Early),
             size=3, stroke=1,
             color="#56B4E9") +

  geom_point(data=ss_bp.all[ss_bp.all$Human_Gene %in%
                              c("CDT1", "CDC6"),],
             aes(x=Log10_Late, y=Log10_Early),
             size=3, stroke=1,
             color="#CC79A7") +

  geom_point(data=ss_bp.all[ss_bp.all$Human_Gene %in%
                              c("TDRKH"),],
             aes(x=Log10_Late, y=Log10_Early),
             size=3, stroke=1,
             color="#F0E442") +
  geom_point(data=ss_bp.all[ss_bp.all$Human_Gene %in%
                              c("GMNN", "CCNB1", "KIF22",
                                "CDT1", "CDC6", "TDRKH"),],
             aes(x=Log10_Late, y=Log10_Early),
             size=3, stroke=1, shape=1,
             color="black") +
    
  geom_point(data=ss_bp.all,
             aes(x=Log10_Late, y=Log10_Early),
             shape=1, size=3, stroke=1, alpha=0.25,
             color="black") +
  geom_hline(yintercept=log10(early_min_kd), linetype="dashed") +
  geom_vline(xintercept=log10(mbt_min_kd), linetype="dashed") +  

  geom_abline(slope=1, intercept=0, color="black", linetype="dotted",
              linewidth=1) +
  labs(x=expression("Gastrulation log"[10]*"(min"^-1*")"),
       y=expression("2-cell log"[10]*"(min"^-1*")")) +
  theme_bw() +
  theme(axis.text=element_text(size=22,colour="black"),
        plot.title = element_text(size=24, hjust=0.5),
        axis.title = element_text(size=20),
        axis.line.x.bottom = element_line(color = "black", linewidth = 1),
        axis.line.y.left   = element_line(color = "black", linewidth = 1),       
        # panel.border = element_rect(linewidth=1.5),
        panel.border = element_blank(),        
        panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
        aspect.ratio = 1,
        legend.position = "none",
        plot.margin = margin(t = 10, r = 25, b = 10,
                             l = 10, unit = "pt")) +
  coord_cartesian(xlim = c(-6,-1),
                  ylim = c(-6,-1))

# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=5, height=5, res=300)

p1
```

![](figures/Comp_Frog_Early-MBT/Biplot%20of%20kDs%20between%20experiments-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
Frog_T5_all.min <- c(161, 221, 281, 401, 521, 641,
                     881, 1605, 1961, 2918, 4405)
Frog_T5_light.min <- Frog_T5_all.min[c(1:2, 4:5, 7:8, 10:11)]

Frog_T6_all.min <- c(150, 210, 270, 390, 510, 635,
                     930, 1591, 1950, 2910, 4394)
Frog_T6_light.min <- Frog_T6_all.min[c(1:2, 4:5, 7:8, 10:11)]

#--------------------------------

early_heavy_time <- (Frog_T5_all.min + Frog_T6_all.min)/2
early_light_time <- (Frog_T5_light.min + Frog_T6_light.min)/2

mbt_heavy_time <- c(1472, 1502, 1535, 1592, 1714, 1832,
                    2412, 2934, 3270, 3749, 4260, 4411)
mbt_light_time <- mbt_heavy_time[c(1,3,5,7,9,11:12)]
```

``` r
plot_protein <- function(protein_id, title_color="black") {
  
  sub_early <- early_fits.df[early_fits.df$Protein_ID == protein_id,]
  sub_mbt <- mbt_fits.df[mbt_fits.df$Protein_ID == protein_id,]

  early_light <- as.numeric(sub_early[10:17])
  early_light <- early_light / early_light[1]
  
  early_label <- as.numeric(sub_early[18:28])
  early_label <- early_label / early_label[1]
      
  mbt_light <- as.numeric(sub_mbt[10:16])
  mbt_light <- mbt_light / mbt_light[1]  
  
  mbt_label <- as.numeric(sub_mbt[17:28])  
  mbt_label <- mbt_label / mbt_label[1]
    
  all_data <- data.frame(Time=c(early_light_time, early_heavy_time,
                                mbt_light_time, mbt_heavy_time),
                         Values=c(early_light, early_label,
                                  mbt_light, mbt_label),
                         Group=c(rep("2C_Light", length(early_light)),
                                 rep("2C_Label", length(early_label)),
                                 rep("Gast_Light", length(mbt_light)),
                                 rep("Gast_Label", length(mbt_label))))
  all_data$Time <- all_data$Time/60
  all_data <- all_data[!all_data$Group %in% c("2C_Light", "Gast_Light"),]
  
  
  p1 <- ggplot() +
    geom_line(data=all_data, aes(x=Time, y=Values, color=Group),
              size=2) +
    geom_hline(yintercept = 1, linetype="dashed", linewidth=1,
               color="black") +
    geom_vline(xintercept=early_heavy_time[1]/60,
               linetype="dotted", linewidth=1) +
    geom_vline(xintercept=mbt_heavy_time[1]/60,
               linetype="dotted", linewidth=1) + 
    annotate("text", x=30, y=1.15, label="2-cell", color="#A63D00",
             hjust=0, size=10) +
    annotate("text", x=50, y=1.15, label="Gast.", color="#E69F00",
             hjust=0, size=10) +
    
    annotate("text", x=45, y=0.7,
             label=paste(round(as.numeric(sub_early["HL_Hrs"]),1), "hrs"),
             color="#A63D00", hjust=0, size=10) +
    annotate("text", x=45, y=0.5,
             label=paste(round(as.numeric(sub_mbt["HL_Hrs"]),1), "hrs"),
             color="#E69F00", hjust=0, size=10) +
    # annotate("text", x=45, y=0.5,
    #          label="Stable",
    #          color="#E69F00", hjust=0, size=10) +
    
    labs(title=as.character(sub_early["Human_Gene"])) +
    theme_bw() +
    theme(axis.text=element_text(size=21,colour="black"),
          plot.title = element_text(size=26, hjust=0.5, color=title_color),
          # axis.title = element_text(size=20),
          axis.title = element_blank(),          
          panel.border = element_rect(linewidth=2),
          panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
          panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
          legend.position = "none",
          plot.margin = margin(t = 10, r = 25, b = 10,
                               l = 10, unit = "pt")) +
    coord_cartesian(ylim=c(0,1.25)) +
    scale_color_manual(values=c("#A63D00", "#E69F00"))

  return(p1) }

# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=6, height=3.5, res=300)

plot_protein("XBmRNA10569|XBXL10_1g5638")
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once per session.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![](figures/Comp_Frog_Early-MBT/Example%20proteins-1.png)<!-- -->

``` r
plot_protein("XBmRNA51997|XBXL10_1g27704")
```

![](figures/Comp_Frog_Early-MBT/Example%20proteins-2.png)<!-- -->

``` r
plot_protein("XBmRNA79044|XBXL10_1g42043")
```

![](figures/Comp_Frog_Early-MBT/Example%20proteins-3.png)<!-- -->

``` r
plot_protein("XBmRNA34319|XBXL10_1g18654")
```

![](figures/Comp_Frog_Early-MBT/Example%20proteins-4.png)<!-- -->

``` r
plot_protein("XBmRNA79355|XBXL10_1g42215")
```

![](figures/Comp_Frog_Early-MBT/Example%20proteins-5.png)<!-- -->

``` r
plot_protein("XBmRNA69786|XBXL10_1g37002")
```

![](figures/Comp_Frog_Early-MBT/Example%20proteins-6.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
combined_ts["Class_Shift"] <- combined_ts$Combined_Class
combined_ts[combined_ts$Class_Shift ==
              "Deg-Long_Deg-Long",
            "Class_Shift"] <- "Always Long"
combined_ts[combined_ts$Class_Shift ==
              "Deg-Long_Deg",
            "Class_Shift"] <- "Long to Deg"
combined_ts[combined_ts$Class_Shift ==
              "Deg-Long_Flat",
            "Class_Shift"] <- "Long to Stable"

combined_ts[combined_ts$Class_Shift ==
              "Deg_Deg",
            "Class_Shift"] <- "Always Deg"
combined_ts[combined_ts$Class_Shift ==
              "Deg_Deg-Long",
            "Class_Shift"] <- "Deg to Long"
combined_ts[combined_ts$Class_Shift ==
              "Deg_Flat",
            "Class_Shift"] <- "Deg to Stable"

combined_ts[combined_ts$Class_Shift ==
              "Flat_Flat",
            "Class_Shift"] <- "Always Stable"
combined_ts[combined_ts$Class_Shift ==
              "Flat_Deg-Long",
            "Class_Shift"] <- "Stable to Long"
combined_ts[combined_ts$Class_Shift ==
              "Flat_Deg",
            "Class_Shift"] <- "Stable to Deg"

table(combined_ts$Class_Shift)
```

    ## 
    ##     Always Deg    Always Long  Always Stable    Deg to Long  Deg to Stable 
    ##            807             76           4608            266            590 
    ##    Long to Deg Long to Stable  Stable to Deg Stable to Long 
    ##              3            280             61            113

``` r
early_classes <- as.data.frame(table(combined_ts$mClass_MBT))
colnames(early_classes) <- c("Category", "Count")

early_classes$Category <- factor(early_classes$Category,
                                 levels=c("Deg", "Deg-Long", "Flat"))

early_classes$Count / sum(early_classes$Count) * 100
```

    ## [1] 12.801293  6.687243 80.511464

``` r
p1 <- ggplot(early_classes, aes(x = "", y = Count, fill = Category)) +
  geom_bar(stat = "identity", color="black", alpha=0.8) +
  labs(x = NULL,        # Removes the blank x-axis label
       # y = "Number of proteins",
       fill = "Class") +
  theme_bw() +
  theme(axis.text=element_text(size=16,colour="black"),
        plot.title = element_text(size=24, hjust=0.5),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(linewidth=1.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
        legend.position = "none",
        plot.margin = margin(t = 10, r = 25, b = 10,
                             l = 10, unit = "pt")) +
  scale_fill_manual(values=c("#5E4FA2", "#56B4E9", "#999999"))


# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=2, height=5, res=300)

p1
```

![](figures/Comp_Frog_Early-MBT/Example%20of%20how%20bar%20plots%20for%20supplemental%20figures%20were%20generated-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
table(combined_ts$Class_Shift)
```

    ## 
    ##     Always Deg    Always Long  Always Stable    Deg to Long  Deg to Stable 
    ##            807             76           4608            266            590 
    ##    Long to Deg Long to Stable  Stable to Deg Stable to Long 
    ##              3            280             61            113

``` r
summary.data <- data.frame(
  category = c("Long", "Deg", "Stable"),
  count = c(265, 807, 590))

# summary.data <- data.frame(
#   category = c("Long", "Deg", "Stable"),
#   count = c(76, 3, 281))

# summary.data <- data.frame(
#   category = c("Long", "Deg", "Stable"),
#   count = c(112, 61, 4609))

x <- round(summary.data$count / sum(summary.data$count) * 100,1)
names(x) <- summary.data$category
x
```

    ##   Long    Deg Stable 
    ##   15.9   48.6   35.5

``` r
summary.data["category"] <- factor(summary.data$category,
                                   levels=summary.data$category)

# Create a pie chart
deg.pie <- ggplot(summary.data, aes(x = 1, y = count, fill = category)) +
  geom_bar(width = 1, stat = "identity", color = "black", alpha=0.8) +
  coord_polar(theta = "y") +
  xlim(0.5, 1.5) +
  scale_fill_manual(values = c("#56B4E9", "#5E4FA2", "#767676")) +
  theme_void()

deg.pie
```

![](figures/Comp_Frog_Early-MBT/Example%20of%20how%20pie%20charts%20for%20supplemental%20figures%20were%20generated-1.png)<!-- -->

``` r
ggsave('test.png', bg="transparent", plot=deg.pie,width=4,height=4,units="in")
```

``` r
plot_deg_hist <- function(group) {
  
  sub_data <- combined_ts[combined_ts$Class_Shift==group,]$HL_Hrs_Early
  print(median(sub_data))
  sub_data <- data.frame(x=sub_data)
  
  p1 <- ggplot(sub_data, aes(x = x)) +
    geom_histogram(aes(y = after_stat(density)), 
                   binwidth = 5,
                   fill = "grey70", 
                   color = "black", 
                   alpha = 0.7) +
    geom_density(color = "darkred", linewidth = 1.5) +
    geom_vline(xintercept = median(sub_data$x), linetype="dashed",
               linewidth=1) +
    labs(title = group,
         x = "2-cell half-life (Hours)",
         y = "Density") +
    theme_bw() +
    theme(axis.text=element_text(size=21,colour="black"),
          plot.title = element_text(size=22, hjust=0.5),
          axis.title = element_text(size=22),
          panel.border = element_rect(linewidth=2),
          panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
          panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
          legend.position = "none",
          plot.margin = margin(t = 10, r = 25, b = 10,
                               l = 10, unit = "pt")) +
    coord_cartesian(ylim=c(0,0.015))
  
  return(p1) }


# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=7, height=3.5, res=300)

plot_deg_hist("Always Deg")
```

    ## [1] 46.99289

![](figures/Comp_Frog_Early-MBT/Plotting%20HLs%20for%20each%20combined%20group-1.png)<!-- -->

``` r
plot_deg_hist("Deg to Long")
```

    ## [1] 121.9124

![](figures/Comp_Frog_Early-MBT/Plotting%20HLs%20for%20each%20combined%20group-2.png)<!-- -->

``` r
plot_deg_hist("Deg to Stable")
```

    ## [1] 127.5671

![](figures/Comp_Frog_Early-MBT/Plotting%20HLs%20for%20each%20combined%20group-3.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
median(combined_ts[combined_ts$Class_Shift %in% c("Deg to Long", "Deg to Stable"),]$HL_Hrs_Early)
```

    ## [1] 125.0095

``` r
sessionInfo()
```

    ## R version 4.5.3 (2026-03-11 ucrt)
    ## Platform: x86_64-w64-mingw32/x64
    ## Running under: Windows 11 x64 (build 26200)
    ## 
    ## Matrix products: default
    ##   LAPACK version 3.12.1
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.utf8 
    ## [2] LC_CTYPE=English_United States.utf8   
    ## [3] LC_MONETARY=English_United States.utf8
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.utf8    
    ## 
    ## time zone: America/Los_Angeles
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] patchwork_1.3.2 ggplot2_4.0.3   tidyr_1.3.2     dplyr_1.2.0    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.6       compiler_4.5.3     tidyselect_1.2.1   systemfonts_1.3.2 
    ##  [5] scales_1.4.0       textshaping_1.0.5  yaml_2.3.12        fastmap_1.2.0     
    ##  [9] R6_2.6.1           labeling_0.4.3     generics_0.1.4     knitr_1.51        
    ## [13] tibble_3.3.1       pillar_1.11.1      RColorBrewer_1.1-3 rlang_1.1.7       
    ## [17] xfun_0.58          S7_0.2.1           otel_0.2.0         cli_3.6.5         
    ## [21] withr_3.0.3        magrittr_2.0.4     digest_0.6.39      grid_4.5.3        
    ## [25] rstudioapi_0.19.0  lifecycle_1.0.5    vctrs_0.7.1        evaluate_1.0.5    
    ## [29] glue_1.8.0         farver_2.1.2       ragg_1.5.2         rmarkdown_2.31    
    ## [33] purrr_1.2.2        tools_4.5.3        pkgconfig_2.0.3    htmltools_0.5.9
