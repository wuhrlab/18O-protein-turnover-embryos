Fly Turnover Analysis
================
Edward Cruz
2026-05-12

``` r
rm(list=ls(all=T))
library(lsa)
```

    ## Loading required package: SnowballC

``` r
library(fgsea)
library(stats)
library(GO.db)
```

    ## Loading required package: AnnotationDbi

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Loading required package: generics

    ## 
    ## Attaching package: 'generics'

    ## The following objects are masked from 'package:base':
    ## 
    ##     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
    ##     setequal, union

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
    ##     unsplit, which.max, which.min

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: IRanges

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:grDevices':
    ## 
    ##     windows

    ## 

``` r
library(tidyr)
```

    ## 
    ## Attaching package: 'tidyr'

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     expand

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:AnnotationDbi':
    ## 
    ##     select

    ## The following objects are masked from 'package:IRanges':
    ## 
    ##     collapse, desc, intersect, setdiff, slice, union

    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     first, intersect, rename, setdiff, setequal, union

    ## The following object is masked from 'package:Biobase':
    ## 
    ##     combine

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, setequal, union

    ## The following object is masked from 'package:generics':
    ## 
    ##     explain

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(ggplot2)
library(stringr)
library(KEGGREST)
library(rrvgo)
library(org.Hs.eg.db)
```

    ## 

``` r
library(clusterProfiler)
```

    ## clusterProfiler v4.18.4 Learn more at https://yulab-smu.top/contribution-knowledge-mining/
    ## 
    ## Please cite:
    ## 
    ## T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan,
    ## X Fu, S Liu, X Bo, and G Yu. clusterProfiler 4.0: A universal
    ## enrichment tool for interpreting omics data. The Innovation. 2021,
    ## 2(3):100141

    ## 
    ## Attaching package: 'clusterProfiler'

    ## The following object is masked from 'package:AnnotationDbi':
    ## 
    ##     select

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     slice

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     rename

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

``` r
library(minpack.lm)
library(Biostrings)
```

    ## Loading required package: XVector

    ## Loading required package: Seqinfo

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

``` r
library(Peptides)
library(ggpubr)
library(forcats)
library(patchwork)

knitr::opts_chunk$set(fig.path = "figures/fly_turnover_analysis/")
```

``` r
#Ultimately - collected these minutes apart so used one vector
Fly_T6.time.min <- c(257, 287, 320, 351, 386, 441, 500,
                     565, 682, 868, 995, 1218, 1514) - 257
light_time <- c(0, Fly_T6.time.min[c(4, 8, 11, 13)])

Fly_T6_hpf_time <- ((Fly_T6.time.min + 257) + (Fly_T6.time.min + 117))/2
Fly_T6_hpf_time <- Fly_T6_hpf_time / 60
Fly_T6_hpf_light_time <- c(Fly_T6_hpf_time[1], Fly_T6_hpf_time[c(4, 8, 11, 13)])

#------------------------------------------

fly_blast_fits <- 
  read.csv("Files/Fits/Dyn-Model/Fly_O18_T6-7_FinalFits_ALL.csv")

head(fly_blast_fits)
```

    ##   Protein_ID Gene_Symbol         kD    HL_Hrs Estimate      T6_kD      T7_kD
    ## 1     Q9VFD5      CG6966 0.01539471 0.7504170    point 0.01539471         NA
    ## 2     Q9VMA3         cup 0.01486216 0.7773063 interval 0.01538270 0.01434163
    ## 3     Q8IPM1         srl 0.01408818 0.8200105    point         NA 0.01408818
    ## 4     Q9VED4        bard 0.01264941 0.9132802    point         NA 0.01264941
    ## 5     P09085         cad 0.01217044 0.9492226    point         NA 0.01217044
    ## 6     Q9VNG1      CG2182 0.01191980 0.9691817 interval 0.01182441 0.01201520
    ##    deltaBIC     T0_L       N13       N14       N15       N16     T0_H       O1
    ## 1  7.832369 1.430443 0.3653358 0.6719770 1.3041889 1.2280553 3.365396 2.688617
    ## 2 10.631890 3.149509 0.6416273 0.4164807 0.3952220 0.3971613 3.778549 2.663961
    ## 3 11.704248 1.706303 0.5332616 1.4637543 0.9769972 0.3196836 3.141349 2.952435
    ## 4 23.447389 2.822755 1.7116034 0.2114468 0.1230994 0.1310951 4.380005 2.508002
    ## 5 19.866379 2.215896 1.0709636 0.5444190 0.5127122 0.6560089 3.627582 2.951559
    ## 6 18.057844 1.141686 1.0725223 1.0495665 0.9212109 0.8150147 3.075938 2.791074
    ##         O2        O3        O4        O5        O6        O7        O8
    ## 1 1.144322 0.8776167 0.7173222 0.4833720 0.6023776 0.4890450 0.5963080
    ## 2 1.149432 0.9644864 0.6122611 0.5266244 0.4758346 0.5489199 0.4854475
    ## 3 1.405408 1.1108608 0.6885891 0.7144997 0.3689344 0.5798764 0.5538892
    ## 4 1.952674 1.7729541 0.7960269 0.4390828 0.2705379 0.2474848 0.1318265
    ## 5 2.063225 1.0804270 0.5648673 0.5653617 0.3289622 0.3461073 0.3201948
    ## 6 1.871093 1.2094890 0.7646547 0.6387289 0.4852980 0.4144466 0.4046215
    ##          O9       O10       O11       O12 mClass
    ## 1 0.5719438 0.5362716 0.3889168 0.5384916    Deg
    ## 2 0.4306688 0.4159514 0.4550532 0.4928102    Deg
    ## 3 0.3652615 0.4694481 0.2848861 0.3645625    Deg
    ## 4 0.1570708 0.1219197 0.1053457 0.1170699    Deg
    ## 5 0.2565727 0.2858838 0.3429187 0.2663379    Deg
    ## 6 0.3697731 0.3911873 0.2936905 0.2900057    Deg

``` r
max_fit_kd <- log(2)/(tail(Fly_T6.time.min, 1)*3)
max_fit_hl <- log(2)/max_fit_kd/60
cat("Max Half-life (hrs):\t", max_fit_hl) 
```

    ## Max Half-life (hrs):  62.85

``` r
fly_GO_map <- read.csv("Systems/Fly_GO-Terms.csv")
fly_GO_map <- fly_GO_map[fly_GO_map$Protein_ID %in%
                           fly_blast_fits$Protein_ID,]
fly_GO_map <- apply(fly_GO_map, 1, function(x){

  BP_terms <- trimws(strsplit(x["GO_BP"], ";")[[1]])
  BP_id <- sub(".*\\[(GO:\\d+)\\].*", "\\1", BP_terms)
  BP_desc <- sub("\\s*\\[GO:.*\\]", "", BP_terms)

  MF_terms <- trimws(strsplit(x["GO_MF"], ";")[[1]])
  MF_id <- sub(".*\\[(GO:\\d+)\\].*", "\\1", MF_terms)
  MF_desc <- sub("\\s*\\[GO:.*\\]", "", MF_terms)

  CC_terms <- trimws(strsplit(x["GO_CC"], ";")[[1]])
  CC_id <- sub(".*\\[(GO:\\d+)\\].*", "\\1", CC_terms)
  CC_desc <- sub("\\s*\\[GO:.*\\]", "", CC_terms)

  out_df <- data.frame(Protein_ID=rep(x["Protein_ID"],
                                      length(c(BP_id, MF_id, CC_id))),
                       ID=c(BP_id, MF_id, CC_id),
                       Desc=c(BP_desc, MF_desc, CC_desc),
                       Type=c(rep("BP", length(BP_id)),
                              rep("MF", length(MF_id)),
                              rep("CC", length(CC_id))))

  return(out_df) }) %>% bind_rows

#This dataframe is solely to count terms and remove <10 annotations
fly_GO_counts <- data.frame(table(fly_GO_map$ID))
fly_GO_counts <- fly_GO_counts[fly_GO_counts$Freq > 20,]

#df of GO terms for identified and mapped proteins
fly_GO_annot <- unique(fly_GO_map[c("ID", "Type", "Desc")])
fly_GO_annot <- fly_GO_annot[fly_GO_annot$ID %in% fly_GO_counts$Var1,]
row.names(fly_GO_annot) <- NULL

head(fly_GO_annot)
```

    ##           ID Type                                  Desc
    ## 1 GO:0007143   BP       female meiotic nuclear division
    ## 2 GO:0008104   BP    intracellular protein localization
    ## 3 GO:0000226   BP microtubule cytoskeleton organization
    ## 4 GO:0000278   BP                    mitotic cell cycle
    ## 5 GO:0007052   BP          mitotic spindle organization
    ## 6 GO:0008017   MF                   microtubule binding

``` r
enrich_kd <- fly_blast_fits[fly_blast_fits$Protein_ID %in% fly_GO_map$Protein_ID,
                            c("Protein_ID", "kD", "O12")]
enrich_kd <- enrich_kd %>% arrange(kD, desc(O12))

enrich_kd["Rank"] <- seq(1, nrow(enrich_kd))
enrich_kd["Signed_Rank"] <- enrich_kd$Rank - mean(enrich_kd$Rank)
enrich_kd["Signed_Rank"] <- enrich_kd$Signed_Rank / max(enrich_kd$Signed_Rank)
enrich_kd <- enrich_kd %>%
  select(Signed_Rank, Protein_ID, kD) %>%
  arrange(desc(Signed_Rank))

protein_stats <- setNames(enrich_kd$Signed_Rank, enrich_kd$Protein_ID)

head(protein_stats)
```

    ##    Q9VFD5    Q9VMA3    Q8IPM1    Q9VED4    P09085    Q9VNG1 
    ## 1.0000000 0.9996429 0.9992858 0.9989288 0.9985717 0.9982146

``` r
cat("\n\n")
```

``` r
tail(protein_stats)
```

    ##     Q9VTC4     A1Z9S1     Q7KBL8     Q7K2W3     A1ZB61     Q9VPN3 
    ## -0.9982146 -0.9985717 -0.9989288 -0.9992858 -0.9996429 -1.0000000

``` r
filtered_fly_GO <- fly_GO_map
filtered_fly_GO <- filtered_fly_GO[filtered_fly_GO$ID %in% fly_GO_annot$ID,]

fly_GO_ref <- lapply(unique(filtered_fly_GO$ID),
       function(x){

         term <- fly_GO_map[fly_GO_map$ID == x,]$Protein_ID

         return(term) })

names(fly_GO_ref) <- unique(filtered_fly_GO$ID)

fly_GO_ref[1]
```

    ## $`GO:0007143`
    ##  [1] "A0A0B4JD97" "A0A0B4KG66" "A0A1Z1CH00" "P18431"     "P25992"    
    ##  [6] "P34739"     "P47938"     "P52304"     "Q24087"     "Q24141"    
    ## [11] "Q24152"     "Q27297"     "Q27889"     "Q7KNM2"     "Q7KRY6"    
    ## [16] "Q95TJ9"     "Q95TN8"     "Q9NI63"     "Q9U405"     "Q9VEZ3"    
    ## [21] "Q9VMA3"     "Q9VVN4"     "Q9W0S8"     "Q9XZL8"     "Q9VU45"    
    ## [26] "Q0KHR4"

``` r
# Run enrichment
set.seed(123)   # reproducible fgsea permutations
res <- fgseaMultilevel(
  pathways = fly_GO_ref,
  stats    = protein_stats,
  minSize  = 10,
  maxSize  = 500)

#Order your results by significance
res <- res[order(res$padj), ]

# fgsea's built-in collapsing function
# Compares the leading edge genes of nested pathways and drops redundant ones
collapsed_pathways <- collapsePathways(res, fly_GO_ref, protein_stats)

#Filter main results to keep only the 'parent' representative terms
res <- res[res$pathway %in% collapsed_pathways$mainPathways, ]
res <- res[res$padj < 0.05,] %>% as.data.frame
res <- merge(res, fly_GO_annot, by.x="pathway",  by.y="ID")
res <- res %>% arrange(desc(NES))
res <- res[c("pathway", "NES", "padj", "size", "Type", "Desc")]

head(res[c("NES", "Type", "Desc")])
```

    ##        NES Type
    ## 1 2.969762   MF
    ## 2 2.823140   CC
    ## 3 2.820233   MF
    ## 4 2.777002   MF
    ## 5 2.690267   BP
    ## 6 2.672947   BP
    ##                                                                    Desc
    ## 1 RNA polymerase II cis-regulatory region sequence-specific DNA binding
    ## 2                                                apical plasma membrane
    ## 3 DNA-binding transcription factor activity, RNA polymerase II-specific
    ## 4                             DNA-binding transcription factor activity
    ## 5                                                         axon guidance
    ## 6                                            motor neuron axon guidance

``` r
custom_GO_names <- c("GO:0051301", "GO:0006260", "GO:0004674", "GO:0000981",
                     "GO:0061630", "GO:0007155",
                     "GO:0006457", "GO:0006635", "GO:0005759", "GO:0000502",
                     "GO:0022626", "GO:0008017")
names(custom_GO_names) <- c("Cell division", "DNA replication", "Protein kinase activity",
                            "Transcription factor activity", "Ubiquitin protein ligase activity",
                            "Cell adhesion", "Protein folding", "Fatty acid β-oxidation",
                            "Mitochondrial matrix", "Proteasome complex", "Cytosolic ribosome",
                            "Microtubule binding")

custom_GO_names <- data.frame(GO=custom_GO_names,
                              Custom_Name=names(custom_GO_names))
custom_GO_terms <- merge(res, custom_GO_names,
                         by.x="pathway", by.y="GO") %>% arrange(desc(NES))
custom_GO_terms
```

    ##       pathway       NES         padj size Type
    ## 1  GO:0000981  2.820233 1.424762e-11   95   MF
    ## 2  GO:0004674  2.611589 1.149820e-09  113   MF
    ## 3  GO:0061630  2.361176 3.011684e-07   91   MF
    ## 4  GO:0007155  2.219842 1.935209e-04   34   BP
    ## 5  GO:0008017  2.217395 3.804691e-05   75   MF
    ## 6  GO:0051301  2.148162 5.124795e-06  112   BP
    ## 7  GO:0006635 -2.117250 6.924532e-04   25   BP
    ## 8  GO:0006457 -2.306698 5.003317e-06   79   BP
    ## 9  GO:0005759 -2.797669 2.483875e-16  183   CC
    ## 10 GO:0000502 -2.920408 9.820419e-10   31   CC
    ## 11 GO:0022626 -3.667790 6.666764e-26   72   CC
    ##                                                                     Desc
    ## 1  DNA-binding transcription factor activity, RNA polymerase II-specific
    ## 2                               protein serine/threonine kinase activity
    ## 3                                      ubiquitin protein ligase activity
    ## 4                                                          cell adhesion
    ## 5                                                    microtubule binding
    ## 6                                                          cell division
    ## 7                                              fatty acid beta-oxidation
    ## 8                                                        protein folding
    ## 9                                                   mitochondrial matrix
    ## 10                                                    proteasome complex
    ## 11                                                    cytosolic ribosome
    ##                          Custom_Name
    ## 1      Transcription factor activity
    ## 2            Protein kinase activity
    ## 3  Ubiquitin protein ligase activity
    ## 4                      Cell adhesion
    ## 5                Microtubule binding
    ## 6                      Cell division
    ## 7             Fatty acid β-oxidation
    ## 8                    Protein folding
    ## 9               Mitochondrial matrix
    ## 10                Proteasome complex
    ## 11                Cytosolic ribosome

``` r
p1 <- ggplot(custom_GO_terms, aes(x = NES, size = size,
                                  y = reorder(Custom_Name, NES))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(alpha = 1) +
  geom_point(shape = 1, color = "black", stroke=0.8) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 50)) +  # wrap for display only
  scale_color_continuous(low = "grey80", high = "grey30", trans = "reverse") +
  labs(
    x = "Normalized enrichment score",
    y = "GO term",
    size = "Set size"
  ) +
  theme_bw() +
  theme(
    legend.text = element_text(size = 12, colour = "black"),
    legend.title = element_text(size = 12, colour = "black"),    
    axis.text = element_text(size = 13, colour = "black"),    
    plot.title = element_text(size = 18, hjust = 0.5),
    axis.title.x = element_text(size=16),
    axis.title.y = element_blank(),    
    panel.border = element_rect(linewidth = 1),
    plot.margin = margin(t = 5, r = 10, b = 15, l = 7, unit = "pt"),
    legend.position = c(0.25,0.75),
    legend.background = element_rect(colour = "black", linewidth = 0.5,
                                     fill = "white")
  ) +
  coord_cartesian(xlim=c(-4,3))

# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=6, height=4, res=300)

p1
```

![](figures/fly_turnover_analysis/GO%20term%20plot-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
target_bins <- data.frame(
  Macro_Compartment = c("Nucleus", "Cytoplasm", "Mitochondrion", "Endoplasmic reticulum", 
                        "Plasma membrane", "Golgi apparatus", "Extracellular region", 
                        "Chromatin", "Vesicle"),
  Macro_GO_ID = c("GO:0005634", "GO:0005737", "GO:0005739", "GO:0005783", 
                  "GO:0005886", "GO:0005794", "GO:0005576", 
                  "GO:0000785", "GO:0031982"))

#---------------------------

offspring_map <- as.list(GOCCOFFSPRING[target_bins$Macro_GO_ID])

mapping_df <- lapply(names(offspring_map),
                               function(parent_id) {
  data.frame(Macro_GO_ID = parent_id,
             Specific_GO_ID = offspring_map[[parent_id]],
             stringsAsFactors = FALSE) }) %>% bind_rows

mapping_df <- rbind(data.frame(Macro_GO_ID = target_bins$Macro_GO_ID,
                               Specific_GO_ID = target_bins$Macro_GO_ID),
                    mapping_df)

mapping_df <- merge(mapping_df, target_bins, by="Macro_GO_ID")
head(mapping_df)
```

    ##   Macro_GO_ID Specific_GO_ID Macro_Compartment
    ## 1  GO:0000785     GO:0033503         Chromatin
    ## 2  GO:0000785     GO:0000785         Chromatin
    ## 3  GO:0000785     GO:0000123         Chromatin
    ## 4  GO:0000785     GO:0000124         Chromatin
    ## 5  GO:0000785     GO:0000786         Chromatin
    ## 6  GO:0000785     GO:0000791         Chromatin

``` r
fly_macro_CC <- fly_blast_fits[c("Protein_ID", "Gene_Symbol",
                                 "kD", "HL_Hrs", "mClass")]
# fly_macro_CC <- fly_macro_CC[!fly_macro_CC$mClass=="Noise // Weak Evidence",]
# fly_macro_CC[!fly_macro_CC$mClass=="Quantifiable degradation",
#              "HL_Hrs"] <- max_fit_hl

fly_macro_CC <- merge(fly_macro_CC,
                      unique(fly_GO_map[c("Protein_ID", "ID")]),
                      by="Protein_ID")
fly_macro_CC <- fly_macro_CC[fly_macro_CC$ID %in%
                               mapping_df$Specific_GO_ID,]

fly_macro_CC <- merge(fly_macro_CC, mapping_df[2:3],
                      by.x="ID", by.y="Specific_GO_ID")
fly_macro_CC["ID"] <- NULL
fly_macro_CC <- unique(fly_macro_CC) %>% arrange(Protein_ID)

CC_prot_counts <- data.frame(table(fly_macro_CC$Protein_ID))
fly_macro_CC <- fly_macro_CC[fly_macro_CC$Protein_ID %in%
                               CC_prot_counts[CC_prot_counts$Freq<=2,]$Var1,]

all_prot_CC <- unique(fly_macro_CC[-7])
all_prot_CC["Macro_Compartment"] <- "All proteins"

fly_macro_CC <- rbind(fly_macro_CC, all_prot_CC)
fly_macro_CC["log_kD"] <- log10(abs(fly_macro_CC$kD))
```

``` r
CC_ordering <- lapply(unique(fly_macro_CC$Macro_Compartment),
       function(x){
         
         sub_df <- fly_macro_CC[fly_macro_CC$Macro_Compartment == x,]
         
         max_hl_rows <- nrow(sub_df[sub_df$HL_Hrs == max_fit_hl,])

         percent_collapsed <- max_hl_rows/nrow(sub_df)
         median_hl <- median(sub_df$HL_Hrs)
          
         out_line <- data.frame(CC=x, HL=median_hl,
                                Collapse_Percent = round(percent_collapsed*100,0),
                                Deg_Percent = 100-round(percent_collapsed*100,0))

         return(out_line) }) %>% bind_rows %>%
  arrange(CC != "All proteins", -HL, -Collapse_Percent)

fly_macro_CC$Macro_Compartment <- factor(fly_macro_CC$Macro_Compartment,
                                         levels=CC_ordering$CC)

head(CC_ordering)
```

    ##                      CC       HL Collapse_Percent Deg_Percent
    ## 1          All proteins 19.74685                2          98
    ## 2         Mitochondrion 26.58172                2          98
    ## 3 Endoplasmic reticulum 25.62422                1          99
    ## 4             Cytoplasm 22.25863                2          98
    ## 5       Golgi apparatus 18.90098                0         100
    ## 6  Extracellular region 18.74924                8          92

``` r
p_kd_CC <- ggplot(data=fly_macro_CC,
                  aes(x = Macro_Compartment, y = HL_Hrs)) +

  geom_boxplot(fill="#5E4FA2", alpha=0.7,
               outliers=FALSE, width = 0.7, color = "black") +
  
  geom_vline(xintercept = 1.5, linetype="dashed",
             color = "black", linewidth = 1) +

  stat_summary(fun = median, geom = "point", 
               shape = 21, size = 3, fill = "red3", color = "black") +

  # geom_text(data = CC_ordering,
  #           aes(x = CC, y = max_fit_hl, label = paste0(Deg_Percent, "%")),
  #           vjust = -0.6, size = 5, inherit.aes = FALSE) +
    
  theme_bw() +
  labs(x = "Cellular Compartment",
       y = "Half-life (hours)") +
  
  theme(axis.text.x = element_text(angle = 35, size=15, hjust = 1, colour="black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=16,colour="black"),
        axis.title.y = element_text(size=16, hjust=0.2),        
        plot.title = element_blank(),
        panel.border = element_rect(linewidth=1.2),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
        plot.margin = margin(t = 10, r = 25, b = 10,
                             l = 10, unit = "pt")) +
  coord_cartesian(ylim = c(0,65)) +
  scale_y_continuous(breaks=seq(0,60,20))
  
p_kd_CC
```

![](figures/fly_turnover_analysis/Bottom%20of%20Compartment%20plot-1.png)<!-- -->

``` r
prop_df <- CC_ordering %>%
  pivot_longer(cols = -c(CC, HL),
               names_to = "Type",
               values_to = "Percent")
prop_df$CC <- factor(prop_df$CC, levels=CC_ordering$CC)

# --- top track: 100% stacked proportion bar ---
p_bar <- ggplot(prop_df, aes(x = CC, y = Percent, fill = Type)) +
  geom_col(width = 0.7, colour = "black", alpha=0.7) +
  geom_vline(xintercept = 1.5, linetype = "dashed", colour = "black", linewidth = 1) +
  scale_fill_manual(values = c(Deg_Percent = "#5E4FA2",
                               Collapse_Percent = "#767676"),
                    labels = c("No measurable degradation", "Degrading")) +
  labs(y = "% of prot.") +
  theme_bw() +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y  = element_text(size = 16, colour = "black"),
        axis.title.y = element_text(size = 16),
        panel.border = element_rect(linewidth = 1.2),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        plot.margin = margin(t = 15, r = 25,
                             b = 0, l = 0, unit = "pt")) +
  scale_y_continuous(breaks = c(0,50,100), expand = expansion(mult = c(0, 0)))

# --- bottom: your existing boxplot, x-title/labels live here ---
p_box <- p_kd_CC + theme(plot.margin = margin(t = 0, r = 25, b = 0, l = 10, unit = "pt"))

# --- stack, aligned x-axis ---
p_combined <- p_bar / p_box + plot_layout(heights = c(1, 4.3))

# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=8, height=5, res=300)

p_combined
```

![](figures/fly_turnover_analysis/Complete%20compartment%20plot-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
fly_deg_proteins <- fly_blast_fits$Protein_ID

interproscan_fly <- 
  read.csv("Systems/Fly_Interproscan/Dmel_Uniprot_Proteome_Interproscan.csv")

interproscan_fly <- interproscan_fly[interproscan_fly$Protein_ID %in%
                                       fly_deg_proteins,]
interproscan_fly <- unique(interproscan_fly[c("Protein_ID", "Signature_Accession",
                                              "Signature_Description")])

sub_IP_fly <- data.frame(table(interproscan_fly$Signature_Accession)) %>%
  arrange(Freq)
sub_IP_fly <- sub_IP_fly[sub_IP_fly$Freq>=4,]$Var1
sub_IP_fly <- sapply(sub_IP_fly, as.character)

length(sub_IP_fly)
```

    ## [1] 455

``` r
sub_IP_fly[1:10]
```

    ##  [1] "PF00002" "PF00046" "PF00082" "PF00104" "PF00105" "PF00128" "PF00136"
    ##  [8] "PF00175" "PF00180" "PF00241"

``` r
domain_sub <- fly_blast_fits[c("Protein_ID", "Gene_Symbol", "kD", "HL_Hrs", "mClass")]
domain_sub <- merge(domain_sub,
                    interproscan_fly[c("Protein_ID", "Signature_Accession",
                                       "Signature_Description")],
                    by="Protein_ID") %>% arrange(HL_Hrs)
# domain_sub <-  domain_sub[domain_sub$Signature_Accession %in%
#                             sub_IP_fly,] %>% arrange(desc(Avg_kD))

head(domain_sub)
```

    ##   Protein_ID Gene_Symbol         kD    HL_Hrs mClass Signature_Accession
    ## 1     Q9VFD5      CG6966 0.01539471 0.7504170    Deg             PF00023
    ## 2     Q9VFD5      CG6966 0.01539471 0.7504170    Deg             PF12796
    ## 3     Q9VMA3         cup 0.01486216 0.7773063    Deg             PF10477
    ## 4     Q8IPM1         srl 0.01408818 0.8200105    Deg             PF00076
    ## 5     P09085         cad 0.01217044 0.9492226    Deg             PF00046
    ## 6     A1ZAP7       Mov10 0.01184692 0.9751438    Deg             PF13087
    ##                                            Signature_Description
    ## 1                                                 Ankyrin repeat
    ## 2                                     Ankyrin repeats (3 copies)
    ## 3 Nucleocytoplasmic shuttling protein for mRNA cap-binding EIF4E
    ## 4                                          RNA recognition motif
    ## 5                                                    Homeodomain
    ## 6                                                     AAA domain

``` r
pfam_domains <- c("PF00225", "PF01437", "PF00096", "PF00595", "PF00567",
                  "PF02984", "PF00028", "PF00176", "PF02207", "PF00439",
                  "PF13445", "PF00632", "PF00069", "PF07525", "PF00046", 
                  "PF00022", "PF00472", "PF00191", "PF01399", "PF00012",
                  "PF00227", "PF00118", "PF00106", "PF00023")

names(pfam_domains) <-
  c("Kinesin motor", "Plexin repeat", "Znf-C2H2", "PDZ", "Tudor",
    "Cyclin-Cterm", "Cadherin", "SNF2-related", "UBR box", "Bromo",
    "RING-Znf", "HECT", "Protein kinase", "SOCS box", "Homeo",
    "Actin", "RF-1", "Annexin", "PCI", "Hsp70",
    "Proteasome subunit", "TCP1-cnp60 chaperonin", "Short chain dehydrogenase",
    "Ankyrin repeat")
```

``` r
domain_all_data <- lapply(pfam_domains, function(pfam){

  sub_df <- domain_sub[domain_sub$Signature_Accession == pfam,]

  sub_df[sub_df$kD<0, "HL_Hrs"] <- max_fit_hl
  sub_df[sub_df$kD<max_fit_kd, "HL_Hrs"] <- max_fit_hl
  sub_df$Signature_Description <- rep(names(pfam_domains[pfam_domains==pfam]),
                                      nrow(sub_df))

  return(sub_df)}) %>% bind_rows

domain_all_data <- domain_all_data[domain_all_data$Signature_Accession %in% sub_IP_fly,]

domain_all_data <- domain_all_data[!domain_all_data$Signature_Description %in%
    c("Actin", "RF-1", "Annexin", "Hsp70", "Proteasome subunit",
      "TCP1-cnp60 chaperonin", "Short chain dehydrogenase"),]

head(domain_all_data)
```

    ##   Protein_ID Gene_Symbol          kD   HL_Hrs mClass Signature_Accession
    ## 1     Q9VKH9        cana 0.002741630 4.213717    Deg             PF00225
    ## 2     Q9VKI0        cmet 0.002211184 5.224555    Deg             PF00225
    ## 3     Q9VIP4         neb 0.001876834 6.155289    Deg             PF00225
    ## 4     Q9VSW5      Klp67A 0.001839068 6.281690    Deg             PF00225
    ## 5     Q9VB25      Klp98A 0.001723566 6.702645    Deg             PF00225
    ## 6     Q9VRK9      Klp64D 0.001722376 6.707277    Deg             PF00225
    ##   Signature_Description
    ## 1         Kinesin motor
    ## 2         Kinesin motor
    ## 3         Kinesin motor
    ## 4         Kinesin motor
    ## 5         Kinesin motor
    ## 6         Kinesin motor

``` r
median(fly_blast_fits[fly_blast_fits$Protein_ID %in% domain_all_data[domain_all_data$Signature_Accession == "PF01399",]$Protein_ID,]$HL_Hrs)
```

    ## [1] 30.1634

``` r
deg_domains <- ggplot(domain_all_data, aes(x = reorder(Signature_Description, -HL_Hrs,
                                                       FUN = median),
                                           y = HL_Hrs)) +

  # geom_hline(yintercept = max_fit_hl, linetype = "dashed",
  #            color = "black", linewidth=1) +
  geom_boxplot(outliers = FALSE, fill = "grey90", color = "black", alpha = 0.7,
               linewidth=0.8) +
  stat_summary(fun = median, geom = "point", 
               shape = 21, size = 3, fill = "red3", color = "black") +
  labs(title = "Protein Half-Lives by InterPro Domain",
       x = "Protein domain", y = "Half-life (hours)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 24,colour="black",
                                   angle=35, hjust=1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=28,colour="black"),
        axis.title.y = element_text(size=28),
        plot.title = element_blank(),
        panel.border = element_rect(linewidth=2),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
        legend.position = "none",
        plot.margin = margin(t = 10, r = 25, b = 10,
                             l = 10, unit = "pt"))

# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=16, height=5.5, res=300)

deg_domains
```

![](figures/fly_turnover_analysis/Domain%20plot-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
fly_iupred <- read.csv("Systems/iupred2a/Dmel_Uniprot_IUPred2A-results.csv")
fly_iupred$Protein_ID <- sapply(strsplit(fly_iupred$Protein_ID,
                                         split = "\\|"),
                                "[", 2)

fly_idr_res <- fly_blast_fits[c("Protein_ID", "Gene_Symbol",
                                "kD", "HL_Hrs", "mClass")]
fly_idr_res <- merge(fly_idr_res, fly_iupred, by="Protein_ID") %>%
  arrange(desc(kD))

fly_idr_res["IDR_Class"] <- "N/A"
fly_idr_res[fly_idr_res$Mean_IUPred>=0.5,"IDR_Class"] <- "Disordered"
fly_idr_res[fly_idr_res$Mean_IUPred<0.5,"IDR_Class"] <- "Ordered"

# --- Are disordered proteins more likely to be degrading AT ALL? ------------
fly_idr_res$IDR_Class <- factor(fly_idr_res$IDR_Class,
                                levels = c("Ordered", "Disordered"))

prop_tab <- table(IDR_Class = fly_idr_res$IDR_Class,
                  Degrading = fly_idr_res$mClass)

# chi-square (Fisher is needlessly slow at this n; result is identical in spirit)
prop_test <- chisq.test(prop_tab)

print(prop_test)
```

    ## 
    ##  Pearson's Chi-squared test with Yates' continuity correction
    ## 
    ## data:  prop_tab
    ## X-squared = 0.46509, df = 1, p-value = 0.4953

``` r
cat("\n\n\n")
```

``` r
# --- AMONG degraders, do disordered ones degrade faster? ------------

deg <- fly_idr_res %>% filter(mClass=="Deg") 

rate_test <- t.test(log10(kD) ~ IDR_Class, data = deg)
print(rate_test)
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  log10(kD) by IDR_Class
    ## t = -13.816, df = 1278.6, p-value < 2.2e-16
    ## alternative hypothesis: true difference in means between group Ordered and group Disordered is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.1756745 -0.1319889
    ## sample estimates:
    ##    mean in group Ordered mean in group Disordered 
    ##                -3.151676                -2.997844

``` r
IDR_boxplot_df <- cbind(fly_idr_res[c("IDR_Class", "HL_Hrs")],
                        data.frame(Set=rep("All", nrow(fly_idr_res))))
IDR_boxplot_df <- rbind(IDR_boxplot_df,
                        cbind(deg[c("IDR_Class", "HL_Hrs")],
                              data.frame(Set=rep("Degrading", nrow(deg)))))

p1 <- ggplot(IDR_boxplot_df,
             aes(x = IDR_Class, y = HL_Hrs, fill = Set,
                 group=interaction(IDR_Class, Set))) +
  geom_boxplot(width = 0.5, color = "black",
               outliers = FALSE, alpha=0.7,
               position=position_dodge(width = 0.6)) +
  stat_summary(fun = median, geom = "point",
               shape = 21, size = 3, fill = "red3", color = "black",
               position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = c(All = "#D55E00", Degrading = "#5E4FA2")) +
  theme_bw() +
  theme(axis.text.y = element_text(size=16, colour="black"),
        axis.text.x = element_text(size = 15, colour="black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_blank(),
        plot.title = element_blank(),
        panel.border = element_rect(linewidth=1.2),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
        legend.position = "none",
        plot.margin = margin(t = 0, r = 25, b = 0, l = 10, unit = "pt")) +
  labs(y=expression("Half-life (hours)")) +
  coord_cartesian(ylim=c(0,60)) +
  scale_y_continuous(breaks=c(0, 10, 30, 50))
```

``` r
dist_idr_df <- 
  rbind(data.frame(table(fly_idr_res[fly_idr_res$IDR_Class=="Ordered",]$mClass)),
        data.frame(table(fly_idr_res[fly_idr_res$IDR_Class=="Disordered",]$mClass)))
dist_idr_df["IDR"] <- c("Ordered", "Ordered", "Disordered", "Disordered")
dist_idr_df[1:2, "Freq"] <- dist_idr_df$Freq[1:2] / 4872
dist_idr_df[3:4, "Freq"] <- dist_idr_df$Freq[3:4] / 993
dist_idr_df$IDR <- factor(dist_idr_df$IDR, levels=c("Ordered", "Disordered"))
dist_idr_df$Var1 <- factor(dist_idr_df$Var1, levels=c("Flat", "Deg"))

p2 <- ggplot(dist_idr_df, aes(x = IDR, y = Freq*100, fill = Var1)) +
  geom_col(width = 0.6, colour = "black", alpha=0.7) +
  scale_fill_manual(values = c(Deg = "#5E4FA2",
                               Flat = "#767676"),
                    labels = c("No measurable degradation", "Degrading")) +
  labs(y = "% of prot.") +
  theme_bw() +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y  = element_text(size = 16, colour = "black"),
        axis.title.y = element_text(size = 18),
        panel.border = element_rect(linewidth = 1.2),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        plot.margin = margin(t = 0, r = 25,
                             b = 0, l = 0, unit = "pt")) +
  scale_y_continuous(breaks = c(0,50,100), expand = expansion(mult = c(0, 0)))


# --- stack, aligned x-axis ---
p_combined <- p2 / p1 + plot_layout(heights = c(1, 2.5))

# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=3.7, height=3.5, res=300)

p_combined
```

![](figures/fly_turnover_analysis/Final%20disorder%20plot-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
all_tracks_df <- fly_blast_fits
all_tracks_df["Light_FC"] <- apply(all_tracks_df, 1, function(row){
  
  light_fc <- as.numeric(row["N16"]) / as.numeric(row["T0_L"])
  
  return(light_fc) })

inc_min_fc <- 1.439371 
dec_min_fc <- 0.6218052 

all_tracks_df["Light_Type"] <- rep("Unchanging", nrow(all_tracks_df))
all_tracks_df[all_tracks_df$Light_FC > inc_min_fc, "Light_Type"] <- "Increasing"
all_tracks_df[all_tracks_df$Light_FC < dec_min_fc, "Light_Type"] <- "Decreasing"

all_tracks_df["Merged_Class"] <- paste0(all_tracks_df$Light_Type, "_", all_tracks_df$mClass)

table(all_tracks_df$Merged_Class)
```

    ## 
    ##  Decreasing_Deg Decreasing_Flat  Increasing_Deg Increasing_Flat  Unchanging_Deg 
    ##             724              44            1157             186            3477 
    ## Unchanging_Flat 
    ##             279

``` r
summary.data <- data.frame(
  category = c("No measurable degradation", "Balanced Protein Turnover",
               "Protein Degradation Only", "Net Protein Increase with Degradation",
               "Protein Synthesis Only"),
  count = c(279+44, 3477, 724, 1157, 186))

sum(summary.data$count)
```

    ## [1] 5867

``` r
summary.data$category
```

    ## [1] "No measurable degradation"            
    ## [2] "Balanced Protein Turnover"            
    ## [3] "Protein Degradation Only"             
    ## [4] "Net Protein Increase with Degradation"
    ## [5] "Protein Synthesis Only"

``` r
cat("\n")
```

``` r
round(summary.data$count / sum(summary.data$count) * 100,1)
```

    ## [1]  5.5 59.3 12.3 19.7  3.2

``` r
summary.data["category"] <- factor(summary.data$category,
                                   levels=c("No measurable degradation",
                                            "Protein Synthesis Only",
                                            "Net Protein Increase with Degradation",
                                            "Balanced Protein Turnover",
                                            "Protein Degradation Only"))

# Create a pie chart
deg.pie <- ggplot(summary.data, aes(x = 1, y = count, fill = category)) +
  geom_bar(width = 1, stat = "identity", color = "black", alpha=0.8,
           linewidth=0.25) +
  coord_polar(theta = "y") +
  xlim(0.5, 1.5) +
  scale_fill_manual(values = c("grey100", "grey80", "grey55", "grey20", "black")) +
  theme_void()

deg.pie
```

![](figures/fly_turnover_analysis/Pie%20chart%20of%20kinetic%20classes-1.png)<!-- -->

``` r
ggsave('test.png', bg="transparent", plot=deg.pie,width=4,height=4,units="in")
```

``` r
plot_protein <- function(protein_id){

  protein_sub <- fly_blast_fits[fly_blast_fits$Protein_ID==protein_id,]  

  print(protein_sub$Gene_Symbol)
  
  light_data <- as.numeric(protein_sub[9:13])
  light_data <- light_data / light_data[1]
  
  heavy_data <- as.numeric(protein_sub[14:26])
  heavy_data <- heavy_data / heavy_data[1]

  plot_data <- data.frame(Time=c(Fly_T6_hpf_light_time,
                                 Fly_T6_hpf_time),
                          Values=c(light_data, heavy_data),
                          Group=c(rep("Control",
                                      length(Fly_T6_hpf_light_time)),
                                  rep("O18",
                                      length(Fly_T6_hpf_time))))  

  p1 <- ggplot() +
    geom_line(data=plot_data[plot_data$Group=="O18",],
              aes(x=Time, y=Values, group=Group),
              color="#36A893", size=3) +
    geom_line(data=plot_data[plot_data$Group=="Control",],
              aes(x=Time, y=Values, group=Group),
              color="#767676", size=3) +    
    theme_bw() +
    theme(axis.text=element_text(size=18,colour="black"),
          plot.title = element_text(size=22, hjust=0.5),
          axis.title = element_blank(),
          panel.border = element_rect(linewidth=2),
          panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
          panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
          legend.position = "none",
          aspect.ratio = 1,
          plot.margin = margin(t = 5, r = 5, b = 5,
                               l = 5, unit = "pt")) +    
    geom_hline(yintercept = 1, linetype="dashed", linewidth=1) +
    coord_cartesian(ylim=c(0,3))
  
  return(p1) }

# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=10, height=5, res=300)

plot_protein("Q94529") |
plot_protein("Q9VEK7") |
plot_protein("Q9W2E9") |
plot_protein("A0A4D6K881") |
plot_protein("Q9VUI0")
```

    ## [1] "Gs1l"

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once per session.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## [1] "CREG"
    ## [1] "Ip6k"
    ## [1] "fwd"
    ## [1] "gnu"

![](figures/fly_turnover_analysis/Example%20proteins-1.png)<!-- -->

``` r
# i <- 1
# i <- i + 1
# print(i)
# all_tracks_df[all_tracks_df$Merged_Class == "Decreasing_Deg","Gene_Symbol"][i]
# 
# plot_protein(all_tracks_df[all_tracks_df$Merged_Class == "Decreasing_Deg","Protein_ID"][i])

# dev.off() #Uncomment for new image!
```

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
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] patchwork_1.3.2        forcats_1.0.1          ggpubr_0.6.3          
    ##  [4] Peptides_2.4.6         Biostrings_2.78.0      Seqinfo_1.0.0         
    ##  [7] XVector_0.50.0         minpack.lm_1.2-4       clusterProfiler_4.18.4
    ## [10] org.Hs.eg.db_3.22.0    rrvgo_1.22.0           KEGGREST_1.50.0       
    ## [13] stringr_1.6.0          ggplot2_4.0.3          dplyr_1.2.0           
    ## [16] tidyr_1.3.2            GO.db_3.22.0           AnnotationDbi_1.72.0  
    ## [19] IRanges_2.44.0         S4Vectors_0.48.1       Biobase_2.70.0        
    ## [22] BiocGenerics_0.56.0    generics_0.1.4         fgsea_1.36.2          
    ## [25] lsa_0.73.4             SnowballC_0.7.1       
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RColorBrewer_1.1-3      rstudioapi_0.19.0       jsonlite_2.0.0         
    ##   [4] tidydr_0.0.6            umap_0.2.10.0           magrittr_2.0.4         
    ##   [7] ggtangle_0.1.2          farver_2.1.2            rmarkdown_2.31         
    ##  [10] ragg_1.5.2              fs_2.1.0                vctrs_0.7.1            
    ##  [13] memoise_2.0.1           ggtree_4.0.5            askpass_1.2.1          
    ##  [16] rstatix_0.7.3           htmltools_0.5.9         broom_1.0.13           
    ##  [19] Formula_1.2-5           gridGraphics_0.5-1      htmlwidgets_1.6.4      
    ##  [22] plyr_1.8.9              cachem_1.1.0            igraph_2.3.2           
    ##  [25] mime_0.13               lifecycle_1.0.5         pkgconfig_2.0.3        
    ##  [28] gson_0.1.0              Matrix_1.7-5            R6_2.6.1               
    ##  [31] fastmap_1.2.0           shiny_1.13.0            digest_0.6.39          
    ##  [34] aplot_0.2.9             enrichplot_1.30.5       colorspace_2.1-2       
    ##  [37] ggnewscale_0.5.2        RSpectra_0.16-2         textshaping_1.0.5      
    ##  [40] RSQLite_3.53.2          labeling_0.4.3          abind_1.4-8            
    ##  [43] polyclip_1.10-7         httr_1.4.8              compiler_4.5.3         
    ##  [46] bit64_4.8.2             fontquiver_0.2.1        withr_3.0.3            
    ##  [49] backports_1.5.1         S7_0.2.1                BiocParallel_1.44.0    
    ##  [52] carData_3.0-6           DBI_1.3.0               ggforce_0.5.0          
    ##  [55] R.utils_2.13.0          ggsignif_0.6.4          MASS_7.3-65            
    ##  [58] openssl_2.4.2           rappdirs_0.3.4          tools_4.5.3            
    ##  [61] otel_0.2.0              ape_5.8-1               scatterpie_0.2.6       
    ##  [64] httpuv_1.6.17           R.oo_1.27.1             glue_1.8.0             
    ##  [67] nlme_3.1-169            GOSemSim_2.36.0         promises_1.5.0         
    ##  [70] grid_4.5.3              gridBase_0.4-7          cluster_2.1.8.2        
    ##  [73] reshape2_1.4.5          snow_0.4-4              gtable_0.3.6           
    ##  [76] R.methodsS3_1.8.2       data.table_1.18.4       car_3.1-5              
    ##  [79] xml2_1.5.2              ggrepel_0.9.8           pillar_1.11.1          
    ##  [82] yulab.utils_0.2.4       later_1.4.8             splines_4.5.3          
    ##  [85] tweenr_2.0.3            treeio_1.34.0           lattice_0.22-9         
    ##  [88] bit_4.6.0               tidyselect_1.2.1        fontLiberation_0.1.0   
    ##  [91] tm_0.7-18               knitr_1.51              fontBitstreamVera_0.1.1
    ##  [94] NLP_0.3-2               xfun_0.58               pheatmap_1.0.13        
    ##  [97] stringi_1.8.7           lazyeval_0.2.3          ggfun_0.2.0            
    ## [100] yaml_2.3.12             evaluate_1.0.5          codetools_0.2-20       
    ## [103] wordcloud_2.6           gdtools_0.5.1           tibble_3.3.1           
    ## [106] qvalue_2.42.0           ggplotify_0.1.3         cli_3.6.5              
    ## [109] xtable_1.8-8            reticulate_1.46.0       systemfonts_1.3.2      
    ## [112] treemap_2.4-4           Rcpp_1.1.1-1.1          png_0.1-9              
    ## [115] parallel_4.5.3          blob_1.3.0              DOSE_4.4.0             
    ## [118] slam_0.1-55             tidytree_0.4.7          ggiraph_0.9.6          
    ## [121] scales_1.4.0            purrr_1.2.2             crayon_1.5.3           
    ## [124] rlang_1.1.7             cowplot_1.2.0           fastmatch_1.1-8
