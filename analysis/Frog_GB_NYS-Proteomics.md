Frog No Spinout Proteomics
================
Gloria Bao
2026-05-04

``` r
rm(list=ls(all=T))
library(tidyr)
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
library(stringr)
library(data.table)
```

    ## 
    ## Attaching package: 'data.table'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     between, first, last

``` r
library(ggplot2)
library(patchwork)
library(TOSTER)
library(minpack.lm)
library(enviPat)
```

    ## 
    ##  
    ##  Welcome to enviPat version 2.8 
    ##  
    ##  Check www.envipat.eawag.ch for an interactive online version
    ##  
    ##  Check www.envimass.ch for a full workflow integration
    ##  

``` r
library(ggtext)

knitr::opts_chunk$set(fig.path = "figures/Frog_GB_NYS-Proteomics/")
```

``` r
human_map <- read.csv("Files/Reference/Xen10_to_Human_Name_Assign.csv")
human_map["Protein_ID"] <- sapply(human_map$Protein_ID, function(x){
  split_id <- strsplit(x, "\\|")[[1]]
  reformat <- paste0(split_id[3], "|", split_id[2])

  return(reformat) })

human_map <- human_map[c("Protein_ID", "XLA_Gene", "Human_Gene")]

frog_yolk_set <- c("XBmRNA35014|XBXL10_1g18990", "XBmRNA39801|XBXL10_1g21498",
                   "XBmRNA35013|XBXL10_1g18989", "XBmRNA35012|XBXL10_1g18988",
                   "XBmRNA67639|XBXL10_1g35803")
names(frog_yolk_set) <- c("vtgb1.L", "vtgb1.S", "vtga1.L", "vtga2.L", "serpina")

keratins <- human_map[grepl("KRT", human_map$Human_Gene),]$Protein_ID
```

``` r
sn.cols <- c("126 Sn", "127n Sn", "127c Sn", "128n Sn", "128c Sn", "129n Sn",
             "129c Sn", "130n Sn", "130c Sn", "131n Sn", "131c Sn", "132n Sn",
             "132c Sn", "133n Sn", "133c Sn", "134n Sn", "134c Sn", "135n Sn")

GB_timepoints <- c(176, 286, 420, 807, 1488, 1920, 2930, 3641, 4406) - 176
```

“XBgroup” IDs are a custom FASTA artifact resulting from shared protein
sequences between multiple Xenbase annotated genes. Thus, we cannot
uniquely identify them, and they are a large source of noise in frog
data.

-\> Not removed in other scripts, but here, we are looking for a stable
and consistent protein set for normalization.

``` r
WI_Rep1 <- read.csv("Data/Gloria/GB_Frog-WideIsolation_NoSpinout_Rep1.csv")
WI_Rep1 <- WI_Rep1[paste0("X", gsub(" ", ".", sn.cols))]
WI_Rep1 <- apply(WI_Rep1, 2, sum)

R1_Corr_Factor <- WI_Rep1[1:9]
R1_Corr_Factor <- R1_Corr_Factor / sum(R1_Corr_Factor)

R2_Corr_Factor <- WI_Rep1[10:18]
R2_Corr_Factor <- R2_Corr_Factor / sum(R2_Corr_Factor)

R1_Corr_Factor
```

    ##   X126.Sn  X127n.Sn  X127c.Sn  X128n.Sn  X128c.Sn  X129n.Sn  X129c.Sn  X130n.Sn 
    ## 0.1008305 0.1146148 0.1147921 0.1158107 0.1129955 0.1167795 0.1079033 0.1092616 
    ##  X130c.Sn 
    ## 0.1070119

``` r
cat("\n")
```

``` r
R2_Corr_Factor
```

    ##   X131n.Sn   X131c.Sn   X132n.Sn   X132c.Sn   X133n.Sn   X133c.Sn   X134n.Sn 
    ## 0.11738052 0.11686939 0.11833482 0.11219488 0.11414736 0.10855309 0.09795259 
    ##   X134c.Sn   X135n.Sn 
    ## 0.10409315 0.11047420

``` r
WI_Rep2 <- read.csv("Data/Gloria/GB_Frog-WideIsolation_NoSpinout_Rep2.csv")
WI_Rep2 <- WI_Rep2[paste0("X", gsub(" ", ".", sn.cols))]
WI_Rep2 <- apply(WI_Rep2, 2, sum)

R3_Corr_Factor <- WI_Rep2[1:9]
R3_Corr_Factor <- R3_Corr_Factor / sum(R3_Corr_Factor)

R4_Corr_Factor <- WI_Rep2[10:18]
R4_Corr_Factor <- R4_Corr_Factor / sum(R4_Corr_Factor)

R3_Corr_Factor
```

    ##   X126.Sn  X127n.Sn  X127c.Sn  X128n.Sn  X128c.Sn  X129n.Sn  X129c.Sn  X130n.Sn 
    ## 0.1086108 0.1173898 0.1122544 0.1152877 0.1071504 0.1195903 0.1057918 0.1065790 
    ##  X130c.Sn 
    ## 0.1073458

``` r
cat("\n")
```

``` r
R4_Corr_Factor
```

    ##  X131n.Sn  X131c.Sn  X132n.Sn  X132c.Sn  X133n.Sn  X133c.Sn  X134n.Sn  X134c.Sn 
    ## 0.1261513 0.1122925 0.1166380 0.1059871 0.1123987 0.1066832 0.1075474 0.1019724 
    ##  X135n.Sn 
    ## 0.1103293

``` r
filter_raw <- function(csv.file, timepoints) {

  #Selecting columns needed for filtering and data analysis
  raw.columns <- c('Protein ID', 'Parsimony', 'Peptide', 'z',  sn.cols)

  data.df <- fread(input = csv.file, select = raw.columns) %>%
    as.data.frame
  data.df["Peptide"] <- sapply(data.df$Peptide, function(x){
    return(str_split(x, "\\.")[[1]][2]) })
  
  #'*Parsimony, REV sequences, contaminants, and oxM*
  data.df <- data.df[data.df$Parsimony %in% c("U", "R"),] #Drop NA Parsimony
  data.df <- data.df[!grepl("#", data.df$`Protein ID`),] #Remove reverse sequences
  data.df <- data.df[!grepl('contaminant', data.df$`Protein ID`),] #Remove contaminants
  data.df <- data.df[!grepl("\\*", data.df$`Peptide`),] #Remove oxM
  data.df <- data.df[data.df$z==2 | data.df$z==3, ]

  #'*Missed Cleavages*
  #Remove peptides with missed cleavages
  # ... Contains KR after removing last AA
  data.df <- data.df[!(grepl("K",
                        str_sub(data.df$Peptide, end=-2))),]
  data.df <- data.df[!(grepl("R",
                        str_sub(data.df$Peptide, end=-2))),]

  #Removing unnecessary columns and renaming data
  data.df["Parsimony"] <- NULL

  return(data.df) }

RTS_Set1 <- filter_raw("Data/Gloria/ORC_07253_RTS_Frog_NoSpinout_Yolk-Proteomics_Shallow.csv")
RTS_Set2 <- filter_raw("Data/Gloria/ORC_07366_RTS_Frog_NoSpinout_Yolk-Proteomics_Shallow.csv")

RTS_Set1 <- RTS_Set1[!grepl("group", RTS_Set1$`Protein ID`),]
RTS_Set1 <- RTS_Set1[!RTS_Set1$`Protein ID` %in% keratins,]

RTS_Set2 <- RTS_Set2[!grepl("group", RTS_Set2$`Protein ID`),]
RTS_Set2 <- RTS_Set2[!RTS_Set2$`Protein ID` %in% keratins,]

head(RTS_Set1)
```

    ##                    Protein ID      Peptide z  126 Sn 127n Sn 127c Sn 128n Sn
    ## 1  XBmRNA35013|XBXL10_1g18989      QHENEQR 3 11.4551 14.1228 13.9901 14.3483
    ## 3  XBmRNA55087|XBXL10_1g29309     EYAESQLK 2 53.8398 55.4299 47.9111 56.9241
    ## 5    XBmRNA7340|XBXL10_1g4020      WGTDEEK 2 79.7770 81.9105 68.5892 70.4796
    ## 7  XBmRNA33552|XBXL10_1g18211     ELGNEAYK 2 49.0822 51.2799 49.8422 52.2617
    ## 9    XBmRNA3457|XBXL10_1g1945 EGDVMMGSQVAR 2 35.2513 37.9659 40.1183 38.9938
    ## 11 XBmRNA40540|XBXL10_1g21843   DMQGMPVTAR 2 64.5773 77.9467 61.4015 71.6063
    ##    128c Sn 129n Sn 129c Sn 130n Sn 130c Sn 131n Sn 131c Sn 132n Sn 132c Sn
    ## 1  12.9634 12.4230 12.2325  8.4804 10.8187  9.6381 10.3571  8.7963 12.7986
    ## 3  56.4190 57.9171 47.8865 47.4001 44.5853 56.0758 56.7916 51.7926 52.6630
    ## 5  70.2971 79.2369 67.6992 67.9748 67.6311 78.7870 79.0480 85.4225 83.8757
    ## 7  45.1650 54.2248 49.0750 49.7434 51.6407 55.8149 52.1540 53.1469 53.1415
    ## 9  36.6759 32.2224 36.6611 31.1191 44.8024 41.8267 38.2939 45.5681 27.6578
    ## 11 77.9972 74.5432 60.8453 71.2854 67.2345 62.1005 65.9374 70.9993 73.6229
    ##    133n Sn 133c Sn 134n Sn 134c Sn 135n Sn
    ## 1   9.8531 11.7858 15.7111  8.3097  8.5812
    ## 3  48.8719 44.3038 41.4107 47.1527 40.3759
    ## 5  79.7862 83.0429 64.1957 76.2393 73.9838
    ## 7  51.3289 49.7720 41.7834 44.8603 50.9481
    ## 9  37.7091 26.2234 38.1023 43.0783 44.5492
    ## 11 60.7431 62.8484 56.3753 61.1369 66.2509

``` r
head(RTS_Set2)
```

    ##                   Protein ID       Peptide z  126 Sn 127n Sn 127c Sn 128n Sn
    ## 1 XBmRNA39801|XBXL10_1g21498      ETHGQTGR 3 59.1885 62.7746 61.7420 61.1299
    ## 3 XBmRNA35014|XBXL10_1g18990      DTHGQTGR 3 74.0743 78.3361 70.6345 60.9218
    ## 4 XBmRNA35014|XBXL10_1g18990      DTHGQTGR 3 34.6047 47.6696 53.3988 39.2941
    ## 6 XBmRNA75253|XBXL10_1g40003    GPSPQEAIQK 2 35.4484 40.3419 40.1487 40.2010
    ## 7 XBmRNA46194|XBXL10_1g24744 SASDEELETPTDK 2  4.2791  7.0040  6.8924  9.6190
    ## 9   XBmRNA2076|XBXL10_1g1265       APDFLHR 3  0.0000  0.0000  4.7702  3.7602
    ##   128c Sn 129n Sn 129c Sn 130n Sn 130c Sn 131n Sn 131c Sn 132n Sn 132c Sn
    ## 1 53.7413 56.5912 48.7337 45.0465 46.3677 62.8648 49.1632 47.6139 49.5930
    ## 3 79.3649 70.9781 57.2185 68.1106 61.1510 71.6556 71.6051 71.3251 66.5476
    ## 4 43.5120 48.2118 37.8944 37.0275 34.7585 54.3180 46.7800 42.9241 39.5627
    ## 6 45.9213 50.0405 40.8960 40.9375 49.5518 50.5885 46.6387 42.0563 32.3125
    ## 7  6.4542 12.1578  7.1485  5.2025  7.5457  5.3150  7.3307 11.6026  7.1319
    ## 9  3.7663  5.1045  3.6110  5.4177  6.4097  5.8688  5.1437  3.4200  2.4084
    ##   133n Sn 133c Sn 134n Sn 134c Sn 135n Sn
    ## 1 57.0508 41.9397 57.6904 35.9604 37.1509
    ## 3 76.1059 67.9077 59.6317 51.4821 57.5756
    ## 4 37.2237 40.6668 36.4775 31.2753 31.9395
    ## 6 39.5679 37.8025 41.7287 41.1385 43.6880
    ## 7  5.6776  9.4893  7.3389  8.9699  9.3376
    ## 9  4.7931  4.9433  5.0523  2.7270  5.4441

``` r
split_experiments <- function(df, cols) {

  new_df <- df[c("Protein ID", "Peptide", "z", cols)]
  colnames(new_df)[1] <- "Protein_ID"
  
  #Apply sn filter for 9 cols
  new_df["sum_sn"] <- rowSums(new_df[cols])
  new_df <- new_df[(new_df$z==2 & new_df$sum_sn>129)|
                     (new_df$z==3 & new_df$sum_sn>265),]

  #Remove charge and sum channels
  new_df <- new_df %>% select(-z) %>%
    group_by(Protein_ID, Peptide) %>%
    dplyr::summarise(across(where(is.numeric), \(x) sum(x)),
                     .groups = "drop") %>%
    as.data.frame()

  #Mean center resulting data and change column names
  new_df[cols] <- new_df[cols] / rowMeans(new_df[cols])
  colnames(new_df)[3:11] <- paste0("T", seq(0,8))
  
  return(new_df) }

GB_R1_Raw <- split_experiments(RTS_Set1, sn.cols[1:9])
GB_R2_Raw <- split_experiments(RTS_Set1, sn.cols[10:18])
GB_R3_Raw <- split_experiments(RTS_Set2, sn.cols[1:9])
GB_R4_Raw <- split_experiments(RTS_Set2, sn.cols[10:18])

head(GB_R1_Raw)
```

    ##                  Protein_ID          Peptide        T0        T1        T2
    ## 1 XBmRNA10020|XBXL10_1g5370         LTLGIIPK 0.9601936 1.0926405 1.0051322
    ## 2 XBmRNA10020|XBXL10_1g5370 WPEVDDDSIEDLGEVK 0.9073059 1.0435138 0.9894413
    ## 3 XBmRNA10074|XBXL10_1g5397  VDQSILTGESVSVIK 1.0075397 0.9862539 0.9282156
    ## 4 XBmRNA10112|XBXL10_1g5418      FVDSQLLPLIK 0.9364854 1.1171375 0.9283796
    ## 5 XBmRNA10112|XBXL10_1g5418         FVIATSTK 0.9425833 1.1145501 0.9269990
    ## 6 XBmRNA10112|XBXL10_1g5418      HQEGEIFDTEK 0.9377173 1.0610844 1.0034841
    ##          T3        T4        T5        T6        T7        T8    sum_sn
    ## 1 1.0968792 0.9908442 1.0181969 1.0537681 0.9022256 0.8801198 1401.5930
    ## 2 1.1434640 0.9154860 0.8663576 1.0427388 0.9581167 1.1335759  245.0210
    ## 3 0.7426104 0.7108594 0.9010631 1.1367185 0.9887605 1.5979789  151.8756
    ## 4 0.9644777 1.0117961 1.1186051 1.0162725 0.9303821 0.9764640 1751.9780
    ## 5 1.0049079 1.0343342 1.0516852 0.9868918 0.9768655 0.9611829 1290.8965
    ## 6 1.0263888 1.0870015 1.1037166 0.9490354 0.9276968 0.9038751 1541.8692

``` r
normalize_exp <- function(df, corr_factor) {

  #Applying wide isolation correction factor to channels
  ratio.df <- as.matrix(df[paste0("T", seq(0,8))])
  corr.ratios <- sweep(ratio.df, 2, corr_factor, `/`)
  corr.ratios <- sweep(corr.ratios, 1, rowMeans(corr.ratios), `/`)  
    
  df[paste0("T", seq(0,8))] <- corr.ratios
  
  return(df) }

GB_R1_Norm <- normalize_exp(GB_R1_Raw, R1_Corr_Factor)
GB_R2_Norm <- normalize_exp(GB_R2_Raw, R2_Corr_Factor)
GB_R3_Norm <- normalize_exp(GB_R3_Raw, R3_Corr_Factor)
GB_R4_Norm <- normalize_exp(GB_R4_Raw, R4_Corr_Factor)
```

``` r
#Determining protein quant using weighted peptides based on signal
protein_rollup <- function(df) {

  roll_df <- lapply(unique(df$Protein_ID), function(x){

    sub_df <- df[df$Protein_ID == x,]
    sub_df["Fraction"] <- sub_df$sum_sn / sum(sub_df$sum_sn)
    sub_df[paste0("T", seq(0,8))] <- sub_df[paste0("T", seq(0,8))] * sub_df$Fraction
    
    out_df <- apply(sub_df[paste0("T", seq(0,8))], 2, sum)
    out_df <- out_df / mean(out_df)
    
    out_df <- data.frame(t(out_df))
    out_df <- cbind(data.frame(Protein_ID = sub_df[1,1]),
                    out_df)

    return(out_df) }) %>% bind_rows

  return(roll_df) }

#Gloria: PQ stands for Protein Quant
GB_R1_PQ <- protein_rollup(GB_R1_Norm)
GB_R2_PQ <- protein_rollup(GB_R2_Norm)
GB_R3_PQ <- protein_rollup(GB_R3_Norm)
GB_R4_PQ <- protein_rollup(GB_R4_Norm)

head(GB_R1_PQ)
```

    ##                    Protein_ID        T0        T1       T2        T3        T4
    ## 1   XBmRNA10020|XBXL10_1g5370 1.0489485 1.0516766 0.970196 1.0585214 0.9628635
    ## 2   XBmRNA10074|XBXL10_1g5397 1.1018408 0.9488487 0.891632 0.7070680 0.6936995
    ## 3   XBmRNA10112|XBXL10_1g5418 1.0341460 1.0810057 0.912364 0.9475512 1.0438143
    ## 4 XBmRNA101811|XBXL10_1g44445 0.9402422 0.8875174 1.059509 0.9954347 1.1450018
    ## 5   XBmRNA10181|XBXL10_1g5453 1.1229956 0.9890111 1.069758 0.9691299 0.9363885
    ## 6 XBmRNA101826|XBXL10_1g44449 1.3337979 1.0625183 0.876152 0.9744851 0.9406689
    ##          T5        T6        T7        T8
    ## 1 0.9468672 1.0829071 0.9255168 0.9525030
    ## 2 0.8508191 1.1616275 0.9978659 1.6465984
    ## 3 1.0525122 1.0142491 0.9293799 0.9849776
    ## 4 0.9846637 0.9319271 1.0341768 1.0215275
    ## 5 1.0856018 1.0027712 0.8665818 0.9577621
    ## 6 0.8731844 0.8475007 0.9890466 1.1026460

``` r
#Calculating median data for plot
GB_AllReps_WI <- rbind(GB_R1_PQ, GB_R2_PQ, GB_R3_PQ, GB_R4_PQ) %>%
  group_by(Protein_ID) %>%
  dplyr::summarise(across(where(is.numeric), \(x) median(x)),
                   .groups = "drop")

#Recentering just in case
GB_AllReps_WI[paste0("T", seq(0,8))] <- 
  sweep(GB_AllReps_WI[paste0("T", seq(0,8))], 1,
        rowMeans(GB_AllReps_WI[paste0("T", seq(0,8))]), `/`)  

consistent_proteins <- intersect(intersect(GB_R1_PQ$Protein_ID, GB_R2_PQ$Protein_ID),
                                 intersect(GB_R3_PQ$Protein_ID, GB_R4_PQ$Protein_ID))

GB_AllReps_WI <- GB_AllReps_WI[GB_AllReps_WI$Protein_ID %in%
                                 consistent_proteins,]
head(GB_AllReps_WI)
```

    ## # A tibble: 6 × 10
    ##   Protein_ID                  T0    T1    T2    T3    T4    T5    T6    T7    T8
    ##   <chr>                    <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
    ## 1 XBmRNA10020|XBXL10_1g53… 1.03  1.06  0.967 1.08  0.960 1.02  0.981 0.917 0.988
    ## 2 XBmRNA10074|XBXL10_1g53… 0.946 0.817 0.802 0.774 0.864 0.836 1.10  1.26  1.60 
    ## 3 XBmRNA10112|XBXL10_1g54… 1.03  0.993 0.982 0.998 1.02  1.04  0.993 0.956 0.991
    ## 4 XBmRNA101811|XBXL10_1g4… 1.04  0.977 1.01  0.998 0.983 1.08  0.971 0.988 0.949
    ## 5 XBmRNA10181|XBXL10_1g54… 1.12  1.01  1.03  0.969 1.02  0.988 0.918 0.910 1.03 
    ## 6 XBmRNA101826|XBXL10_1g4… 0.997 1.04  0.983 1.02  1.05  0.936 0.915 1.01  1.05

``` r
#Gloria: DNU stands for Do Not Use (DNU)
DNU_Raw1_PQ <- protein_rollup(GB_R1_Raw)
DNU_Raw2_PQ <- protein_rollup(GB_R2_Raw)
DNU_Raw3_PQ <- protein_rollup(GB_R3_Raw)
DNU_Raw4_PQ <- protein_rollup(GB_R4_Raw)

GB_AllReps_DNU <- rbind(DNU_Raw1_PQ, DNU_Raw2_PQ, DNU_Raw3_PQ, DNU_Raw4_PQ) %>%
  group_by(Protein_ID) %>%
  dplyr::summarise(across(where(is.numeric), \(x) median(x)),
                   .groups = "drop")

#Recentering just in case
GB_AllReps_DNU[paste0("T", seq(0,8))] <- 
  sweep(GB_AllReps_DNU[paste0("T", seq(0,8))], 1,
        rowMeans(GB_AllReps_DNU[paste0("T", seq(0,8))]), `/`)  

GB_AllReps_DNU <- GB_AllReps_DNU[GB_AllReps_DNU$Protein_ID %in%
                                   consistent_proteins,]
head(GB_AllReps_DNU)
```

    ## # A tibble: 6 × 10
    ##   Protein_ID                  T0    T1    T2    T3    T4    T5    T6    T7    T8
    ##   <chr>                    <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
    ## 1 XBmRNA10020|XBXL10_1g53… 1.01  1.10  1.00  1.12  0.970 1.02  0.938 0.876 0.957
    ## 2 XBmRNA10074|XBXL10_1g53… 0.958 0.850 0.849 0.797 0.883 0.857 1.04  1.18  1.58 
    ## 3 XBmRNA10112|XBXL10_1g54… 1.04  1.03  1.02  0.995 1.04  1.05  0.938 0.908 0.982
    ## 4 XBmRNA101811|XBXL10_1g4… 1.08  1.00  1.06  1.03  0.966 1.07  0.914 0.945 0.937
    ## 5 XBmRNA10181|XBXL10_1g54… 1.12  1.05  1.07  0.990 1.04  0.998 0.841 0.870 1.02 
    ## 6 XBmRNA101826|XBXL10_1g4… 1.05  1.08  1.03  1.01  1.03  0.941 0.857 0.976 1.02

``` r
exp_decay <- function(t_hrs, a, b) {
  return(a * exp(-b * t_hrs)) }

fit_yolk_trend <- function(df1, df2, df3, df4, med_df){

  calc_yolk_vec <- function(df) {
    yolk_med <- apply(df[df$Protein_ID %in% frog_yolk_set,2:10], 2, median)
    yolk_med <- yolk_med / mean(yolk_med)    
    return(yolk_med) }
    
  df1_yolk_med <- calc_yolk_vec(df1)
  df2_yolk_med <- calc_yolk_vec(df2)
  df3_yolk_med <- calc_yolk_vec(df3)
  df4_yolk_med <- calc_yolk_vec(df4)  

  # df1_yolk_med <- df1_yolk_med / df1_yolk_med[1]
  # df2_yolk_med <- df2_yolk_med / df2_yolk_med[1]
  # df3_yolk_med <- df3_yolk_med / df3_yolk_med[1]
  # df4_yolk_med <- df4_yolk_med / df4_yolk_med[1]
  
  #---------------------------------------
  
  yolk_df <- data.frame(Time=rep(c((GB_timepoints + 176)/60), 4),
                        Value=c(df1_yolk_med, df2_yolk_med,
                                df3_yolk_med, df4_yolk_med),
                        Exp=c(rep("T1", 9), rep("T2", 9),
                              rep("T3", 9), rep("T4", 9)))

  lin <- lm(log(Value) ~ Time, data = yolk_df, weights = yolk_df$Value^2)
  a0 <- unname(exp(coef(lin)[1]))
  b0 <- unname(-coef(lin)[2])

  yolk_fit_nls <- nlsLM(Value ~ a * exp(-b * Time),
                        data = yolk_df,
                        start = list(a = a0, b = b0),
                        lower = c(a=0, b=0))

  smooth_fit <-
    data.frame(Time = (GB_timepoints+176)/60)
  smooth_fit["Predict"] <- exp_decay(smooth_fit$Time,
                                     coef(yolk_fit_nls)['a'],
                                     coef(yolk_fit_nls)['b'])

  yolk_meds <- med_df[med_df$Protein_ID %in% frog_yolk_set,]
  yolk_meds <- apply(yolk_meds[paste0("T", seq(0,8))], 2, median)

  yolk_meds <- yolk_meds / mean(yolk_meds)
  # yolk_meds <- yolk_meds / yolk_meds[1]
  yolk_meds <- data.frame(Time = (GB_timepoints+176)/60,
                          Values=yolk_meds)

  p1 <- ggplot() +
    geom_line(data=yolk_df, aes(x=Time, y=Value, group=Exp),
              color="grey70", size=2, alpha=0.7) +
    
    # geom_line(data=yolk_meds,
    #           aes(x=Time, y=Values, group=1),
    #           color="#CC3333", size=3, alpha=0.7) +
    
    geom_line(data=yolk_meds,
              aes(x=Time, y=Values, group=1),
              color="#E69F00", size=3, alpha=0.7) +
    geom_line(data=smooth_fit,
              aes(x=Time, y=Predict, group=1),
              color="purple", size=3, alpha=0.7) +
    
    geom_hline(yintercept = 1, linetype="dotted", linewidth=1.2, color="black") +
    labs(x = "Hours post-fertilization", y = "Rel. abundance") +
    theme_bw() +
    theme(axis.text = element_text(size = 20, colour = "black"),
          plot.title = element_text(size = 20, hjust = 0.5),
          axis.title = element_text(size = 24),
          aspect.ratio = 1,
          legend.text = element_text(size=16),
          legend.title = element_text(size=16),
          panel.border = element_rect(linewidth = 2),
          panel.grid.major = element_line(color = "grey70", linewidth = 0.25),
          panel.grid.minor = element_line(color = "grey70", linewidth = 0.25),
          plot.margin = margin(10, 15, 10, 15, unit = "pt")) +
    coord_cartesian(ylim=c(0,1.5))
   
  print(p1)
  
  return(coef(yolk_fit_nls)) }

# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=4.5, height=4.5, res=300)

# fit_yolk_trend(DNU_Raw1_PQ, DNU_Raw2_PQ, DNU_Raw3_PQ, DNU_Raw4_PQ,
#                GB_AllReps_DNU)
yolk_fit_values <- fit_yolk_trend(GB_R1_PQ, GB_R2_PQ, GB_R3_PQ, GB_R4_PQ,
                                  GB_AllReps_WI)
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once per session.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![](figures/Frog_GB_NYS-Proteomics/Exponential%20decay%20plots-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
generate_yolk_corr <- function(DNU_df, raw_df){
  theory_median <- exp_decay((GB_timepoints+176)/60, yolk_fit_values['a'],
                             yolk_fit_values['b'])
  theory_median <- theory_median / mean(theory_median)
  
  raw_average <- DNU_df[DNU_df$Protein_ID %in% frog_yolk_set,]
  raw_average <- apply(raw_average[2:10], 2, median)
  raw_average <- raw_average / mean(raw_average)
  
  corr_factor <- raw_average / theory_median

  yolk_norm <- raw_df
  yolk_norm[3:11] <- sweep(raw_df[3:11], 2, corr_factor, `/`)
  yolk_norm[3:11] <- sweep(yolk_norm[3:11], 1,
                           rowMeans(yolk_norm[3:11]), `/`)

  yolk_PQ <- protein_rollup(yolk_norm)
  
  return(yolk_PQ) }

GB_R1_YolkPQ <- generate_yolk_corr(DNU_Raw1_PQ, GB_R1_Raw)
GB_R2_YolkPQ <- generate_yolk_corr(DNU_Raw2_PQ, GB_R2_Raw)
GB_R3_YolkPQ <- generate_yolk_corr(DNU_Raw3_PQ, GB_R3_Raw)
GB_R4_YolkPQ <- generate_yolk_corr(DNU_Raw4_PQ, GB_R4_Raw)

head(GB_R1_YolkPQ)
```

    ##                    Protein_ID        T0        T1        T2        T3        T4
    ## 1   XBmRNA10020|XBXL10_1g5370 1.0732356 1.0139600 1.0388836 1.0885937 0.9401716
    ## 2   XBmRNA10074|XBXL10_1g5397 1.1297724 0.9167832 0.9568055 0.7287141 0.6788054
    ## 3   XBmRNA10112|XBXL10_1g5418 1.0594794 1.0436556 0.9782364 0.9757405 1.0205791
    ## 4 XBmRNA101811|XBXL10_1g44445 0.9619463 0.8556291 1.1344402 1.0236412 1.1179415
    ## 5   XBmRNA10181|XBXL10_1g5453 1.1477955 0.9525443 1.1442949 0.9956170 0.9133649
    ## 6 XBmRNA101826|XBXL10_1g44449 1.3652363 1.0248299 0.9385625 1.0025751 0.9188749
    ##          T5        T6        T7        T8
    ## 1 0.9352053 1.0726554 0.9035697 0.9337251
    ## 2 0.8421466 1.1530986 0.9762904 1.6175838
    ## 3 1.0409291 1.0060018 0.9085346 0.9668436
    ## 4 0.9724727 0.9230400 1.0095798 1.0013091
    ## 5 1.0711132 0.9922379 0.8451440 0.9378883
    ## 6 0.8627841 0.8398183 0.9659825 1.0813365

``` r
#Calculating median data for plot
GB_YolkNorm_Merge <- rbind(GB_R1_YolkPQ, GB_R2_YolkPQ,
                           GB_R3_YolkPQ, GB_R4_YolkPQ) %>%
  group_by(Protein_ID) %>%
  dplyr::summarise(across(where(is.numeric), \(x) median(x)),
                   .groups = "drop")

GB_YolkNorm_Merge[paste0("T", seq(0,8))] <- 
  sweep(GB_YolkNorm_Merge[paste0("T", seq(0,8))], 1,
        rowMeans(GB_YolkNorm_Merge[paste0("T", seq(0,8))]), `/`)  

GB_YolkNorm_Merge <- GB_YolkNorm_Merge[GB_YolkNorm_Merge$Protein_ID %in%
                                          consistent_proteins,]

head(GB_YolkNorm_Merge)
```

    ## # A tibble: 6 × 10
    ##   Protein_ID                  T0    T1    T2    T3    T4    T5    T6    T7    T8
    ##   <chr>                    <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
    ## 1 XBmRNA10020|XBXL10_1g53… 1.05  1.03  1.03  1.11  0.923 0.988 0.970 0.913 0.991
    ## 2 XBmRNA10074|XBXL10_1g53… 0.943 0.809 0.845 0.790 0.830 0.820 1.09  1.24  1.63 
    ## 3 XBmRNA10112|XBXL10_1g54… 1.05  0.980 1.01  1.01  0.974 1.02  0.987 0.942 1.02 
    ## 4 XBmRNA101811|XBXL10_1g4… 1.05  0.987 1.05  1.02  0.936 1.04  0.962 0.979 0.971
    ## 5 XBmRNA10181|XBXL10_1g54… 1.11  1.01  1.08  0.988 0.970 0.965 0.903 0.911 1.06 
    ## 6 XBmRNA101826|XBXL10_1g4… 1.01  1.01  1.02  1.04  1.02  0.920 0.909 0.982 1.09

``` r
kmeans_plot <- function(data, time_lst, k=5, seed=123, order=seq(1,k,1)) {
  
  set.seed(seed)

  data <- as.data.frame(data)
  row.names(data) <- data$Protein_ID
  data[c("Protein_ID")] <- NULL

  # K-means clustering
  km.all <- kmeans(data, centers = k)
  
  all.kresult <- as.data.frame(km.all$centers)
  all.kresult["Cluster"] <- row.names(all.kresult)
  all.kresult["Size"] <- km.all$size
  row.names(all.kresult) <- NULL
  
  all.kresult <- all.kresult %>%
    pivot_longer(cols = -c(Cluster, Size),
                 names_to = "Time", values_to = "Value")

  all.kresult["Time_Hrs"] <- 
    rep(time_lst, nrow(all.kresult) / length(time_lst))  

  plots <- lapply(seq_len(k), function(i){
    sub_cl <- all.kresult[all.kresult$Cluster == i,]
    # sub_cl["Value"] <- sub_cl$Value / sub_cl$Value[1]
    
    p1 <- ggplot() +
      geom_line(data = sub_cl, aes(x = Time_Hrs, y = Value, group = 1),
                linewidth = 3, color = "black") +
      geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1) +
      theme_bw() +
      labs(x = "", y = "Rel. Abundance",
           title = unique(sub_cl$Size)) +
      theme(axis.text = element_text(size = 18, colour = "black"),
            plot.title = element_text(size = 18, hjust = 0.5),
            axis.title = element_text(size = 22),
            axis.title.y = element_blank(),
            aspect.ratio = 1,
            panel.border = element_rect(linewidth = 2),
            legend.position = "none",
            panel.grid.major = element_line(color = "grey70", linewidth = 0.25),
            panel.grid.minor = element_line(color = "grey70", linewidth = 0.25)) +
      coord_cartesian(ylim = c(0, 3))
    
    return(p1) })
  
  # Combine all plots with patchwork
  combined_plot <- Reduce(`|`, plots[order])
  
  print(combined_plot)
  
  out.clusters <- data.frame(Protein_ID = names(km.all$cluster),
                             Cluster = km.all$cluster)
  row.names(out.clusters) <- NULL

  return(out.clusters) }  
  

# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=10, height=4, res=300)

yolk_kmeans <- kmeans_plot(GB_YolkNorm_Merge,
                           time_lst=(GB_timepoints+176)/60,
                           order=c(4,3,1,5,2))
```

![](figures/Frog_GB_NYS-Proteomics/kmeans%20plots%20for%20each%20dataset-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
fk_fit_proteins <- yolk_kmeans[yolk_kmeans$Cluster %in% c(4,2),]$Protein_ID

cat("FK Fit Proteins: ", length(fk_fit_proteins), "\n\n")
```

    ## FK Fit Proteins:  911

``` r
fk_fit_proteins[1:5]
```

    ## [1] "XBmRNA10020|XBXL10_1g5370"   "XBmRNA10112|XBXL10_1g5418"  
    ## [3] "XBmRNA101811|XBXL10_1g44445" "XBmRNA10181|XBXL10_1g5453"  
    ## [5] "XBmRNA101826|XBXL10_1g44449"

``` r
external_norm_set <- yolk_kmeans[yolk_kmeans$Cluster == 4,]$Protein_ID

external_norm_set <- sapply(external_norm_set, function(x){
  
  all_values <- as.numeric(GB_YolkNorm_Merge[GB_YolkNorm_Merge$Protein_ID == x,2:10])
  max_change <- max(all_values)/min(all_values)
  
  if (max_change < 1.2 & max_change > 0.8) {
    return(x) }
  else { return(NA) }}, USE.NAMES = FALSE)

external_norm_set <- external_norm_set[!is.na(external_norm_set)]

cat("Potential Stable Proteins: ", length(external_norm_set), "\n\n")
```

    ## Potential Stable Proteins:  656

``` r
external_norm_set[1:5]
```

    ## [1] "XBmRNA10112|XBXL10_1g5418"   "XBmRNA101811|XBXL10_1g44445"
    ## [3] "XBmRNA101837|XBXL10_1g44455" "XBmRNA10209|XBXL10_1g5470"  
    ## [5] "XBmRNA1039|XBXL10_1g760"

``` r
filter_lblf_data <- function(data.df) {

  data.df["Protein.ID"] <- #Reorganizing order of IDs to match T8 search
    paste0(sapply(strsplit(data.df$`Protein.ID`, split = "\\|"), "[", 3), "|",
           sapply(strsplit(data.df$`Protein.ID`, split = "\\|"), "[", 2)) 

  data.df <- data.df[data.df$Parsimony %in% c("U"),] #Drop NA Parsimony
  data.df <- data.df[!grepl("#", data.df$`Protein.ID`),] #Remove reverse sequences
  data.df <- data.df[!grepl('contaminant', data.df$`Protein.ID`),] #Remove contaminants
  data.df <- data.df[!grepl("\\*", data.df$`Trimmed.Peptide`),] #Remove oxM

  #Remove peptides with missed cleavages
  # ... Contains KR after removing last AA
  data.df <- data.df[!(grepl("K",
                        str_sub(data.df$`Trimmed.Peptide`, end=-2))),]
  data.df <- data.df[!(grepl("R",
                        str_sub(data.df$`Trimmed.Peptide`, end=-2))),]

  #Max signal within peptide group
  data.df <- data.df %>% group_by(`Trimmed.Peptide`) %>%
    slice_max(`Sum.S.N`, with_ties = FALSE) %>%
    ungroup

  data.df <- data.df[c('Protein.ID', 'Trimmed.Peptide', "z",
                       'ScanF', "Parent.Scan", "Isolation.m.z", "Theo.m.z",
                       paste("Max.Sn.M", seq(0,90,1), sep=""))]
                       # paste("Area.M", seq(0,90,1), sep=""))]

  data.df <- data.df[data.df$`Max.Sn.M0` > 20 &
                       data.df$`Max.Sn.M1` > 20 , ]
  
  return(data.df) }

T5_O12_lblfree <- read.csv("Data/XLA_O18/LabelFree/ORC_05461_XLA-O18-T5_LabelFree_O12.csv")
T5_O12_lblfree <- filter_lblf_data(T5_O12_lblfree)

head(T5_O12_lblfree)
```

    ## # A tibble: 6 × 98
    ##   Protein.ID     Trimmed.Peptide     z  ScanF Parent.Scan Isolation.m.z Theo.m.z
    ##   <chr>          <chr>           <int>  <int>       <int>         <dbl>    <dbl>
    ## 1 XBmRNA1930|XB… AAAIAFAK            2  42249       42181          382.     382.
    ## 2 XBmRNA33552|X… AAALEFLNR           2 134175      134116          503.     503.
    ## 3 XBmRNA26944|X… AADELNAFLAEEADK     2 198919      198887          804.     804.
    ## 4 XBmRNA7505|XB… AAEADQIIEYLK        2 210260      210203          683.     682.
    ## 5 XBmRNA80431|X… AAEAVDDIPFGITS…     2 266050      265980         1057.    1057.
    ## 6 XBmRNA60082|X… AAEDGMGEYLFDK       2 154388      154331          724.     723.
    ## # ℹ 91 more variables: Max.Sn.M0 <dbl>, Max.Sn.M1 <dbl>, Max.Sn.M2 <dbl>,
    ## #   Max.Sn.M3 <dbl>, Max.Sn.M4 <dbl>, Max.Sn.M5 <dbl>, Max.Sn.M6 <dbl>,
    ## #   Max.Sn.M7 <dbl>, Max.Sn.M8 <dbl>, Max.Sn.M9 <dbl>, Max.Sn.M10 <dbl>,
    ## #   Max.Sn.M11 <dbl>, Max.Sn.M12 <dbl>, Max.Sn.M13 <dbl>, Max.Sn.M14 <dbl>,
    ## #   Max.Sn.M15 <dbl>, Max.Sn.M16 <dbl>, Max.Sn.M17 <dbl>, Max.Sn.M18 <dbl>,
    ## #   Max.Sn.M19 <dbl>, Max.Sn.M20 <dbl>, Max.Sn.M21 <dbl>, Max.Sn.M22 <dbl>,
    ## #   Max.Sn.M23 <dbl>, Max.Sn.M24 <dbl>, Max.Sn.M25 <dbl>, Max.Sn.M26 <dbl>, …

``` r
# Amino acid residue compositions (no terminal groups)
aa_composition <- list(
  A = c(C=3, H=5, N=1, O=1, S=0),
  V = c(C=5, H=9, N=1, O=1, S=0),
  I = c(C=6, H=11, N=1, O=1, S=0),
  L = c(C=6, H=11, N=1, O=1, S=0),
  M = c(C=5, H=9, N=1, O=1, S=1),  
  F = c(C=9, H=9, N=1, O=1, S=0),
  Y = c(C=9, H=9, N=1, O=2, S=0),  
  W = c(C=11, H=10, N=2, O=1, S=0),
  C = c(C=3, H=5, N=1, O=1, S=1),  
  G = c(C=2, H=3, N=1, O=1, S=0),  
  P = c(C=5, H=7, N=1, O=1, S=0),
  S = c(C=3, H=5, N=1, O=2, S=0),
  T = c(C=4, H=7, N=1, O=2, S=0),
  N = c(C=4, H=6, N=2, O=2, S=0),
  Q = c(C=5, H=8, N=2, O=2, S=0),  
  H = c(C=6, H=7, N=3, O=1, S=0),  
  #Can consider adding +1H if problems arise
  R = c(C=6, H=12, N=4, O=1, S=0),
  K = c(C=6, H=12, N=2, O=1, S=0),  
  #Can consider subtracting -1H if problems arise  
  D = c(C=4, H=5, N=1, O=3, S=0),
  E = c(C=5, H=7, N=1, O=3, S=0))


get_light_formula <- function(sequence, charge) {
  total <- c(C=0, H=0, N=0, O=0, S=0)
  total_O <- 0
  
  for (aa in strsplit(sequence, "")[[1]]) {
    total <- total + aa_composition[[aa]] 
    
    if (aa == "C") {
      # Add NEM modification per cysteine
      total["C"] <- total["C"] + 6
      total["H"] <- total["H"] + 7
      total["N"] <- total["N"] + 1
      total["O"] <- total["O"] + 2 }    
    
    }
  
  # Add terminal H2O
  total["H"] <- total["H"] + 2
  total["O"] <- total["O"] + 1
  
  # Add proton(s) for charged ion
  total["H"] <- total["H"] + charge

  return(total) }

data("isotopes")
```

``` r
convolve_spectra <- function(iso_df, theo_mz, mz_max) {
  
  #kwargs - Orbitrap @120k MS1
  resolution_at_200 <- 120000
  step=0.0005
  
  #Orbitrap resolution decreases with sqrt(m/z)
  Res <- resolution_at_200 * sqrt(200 / theo_mz) #Resolution
  fwhm <- theo_mz / Res #full width at half maximum
  sigma <- fwhm / (2*sqrt(2*log(2))) #Deviation for Gaussian dist.
    
  #Calulate range and creating empty vector for convolving
  mz_min <- theo_mz - 1
  mz_max <- mz_max + 1
  mz_grid <- seq(mz_min, mz_max, by = step)
  smooth_abundance <- numeric(length(mz_grid))
  
  for (i in seq_len(nrow(iso_df))) {
    #Loop each m/z and intensity
    mz <- iso_df$`m/z`[i]
    intensity <- iso_df$abundance[i]
    
    #Generate gaussian distribution assuming:
    # ... a mean at m/z with sigma from resolution
    peak_dist <- intensity * dnorm(mz_grid, mean = mz, sd = sigma)
    smooth_abundance <- smooth_abundance +
      intensity * dnorm(mz_grid, mean = mz, sd = sigma) #Combine for all m/z
  }
  
  #Normalize maximum to 100 and convert to dataframe
  smooth_abundance <- 100 * smooth_abundance / max(smooth_abundance)
  gaus_convolved <- data.frame(mz = mz_grid, abundance = smooth_abundance)  

  return(gaus_convolved) }
```

``` r
# Class A: no side-chain O (m=2, s=1): R in {0,1}, k in {0,1,2}
# rows = k, cols = R
M_A <- matrix(c(
  1,   0,    # k=0
  1/2, 1/2,  # k=1
  0,   1     # k=2
), nrow=3, byrow=TRUE)

# Class B: one side-chain O (S,N,Q) (m=3, s=2): R in {0,1,2}, k in {0..3}
M_B <- matrix(c(
  1,   0,   0,   # k=0
  1/3, 2/3, 0,   # k=1
  0,   2/3, 1/3, # k=2
  0,   0,   1    # k=3
), nrow=4, byrow=TRUE)

# Class C: two side-chain O (D,E) (m=4, s=3): R in {0,1,2,3}, k in {0..4}
M_C <- matrix(c(
  1,   0,   0,   0,   # k=0
  1/4, 3/4, 0,   0,   # k=1
  0,   1/2, 1/2, 0,   # k=2
  0,   0,   3/4, 1/4, # k=3
  0,   0,   0,   1    # k=4
), nrow=5, byrow=TRUE)
```

``` r
XLA_AA_df <- read.csv("Data/XLA_O18/Isotopic_Envelopes/Frog_Early_AA_Matrix.csv")

fk_map <- list(
  X = as.numeric(XLA_AA_df[XLA_AA_df$AA=="X" &
                             XLA_AA_df$Sample=="col009b_T8",
                c("C12.PARENT", "O18.label.1", "O18.label.2")]),
  U = as.numeric(XLA_AA_df[XLA_AA_df$AA=="U" &
                        XLA_AA_df$Sample=="col009b_T8",
                c("C12.PARENT", "O18.label.1", "O18.label.2", "O18.label.3")]),
  Z = as.numeric(XLA_AA_df[XLA_AA_df$AA=="Z" &
                        XLA_AA_df$Sample=="col009b_T8",
                c("C12.PARENT", "O18.label.1", "O18.label.2",
                  "O18.label.3", "O18.label.4")])
)

fk_map
```

    ## $X
    ## [1] 0.1580079 0.4007321 0.4412600
    ## 
    ## $U
    ## [1] 0.03150128 0.18970297 0.39823788 0.38055787
    ## 
    ## $Z
    ## [1] 0.02409882 0.03165132 0.09513571 0.33629094 0.51282321

``` r
ms_for_residue <- function(aa) {
  classA <- c("R","H","K","C","G","P","A","V",
              "I","L","M","F","W","T","Y") # per essential rule
  classB <- c("S","N","Q")
  classC <- c("D","E")
  if (aa %in% classA) return("X")
  if (aa %in% classB) return("U")
  if (aa %in% classC) return("Z")
  stop(sprintf("Unknown residue %s", aa))
}

matrix_for_class <- function(residue){
  switch(ms_for_residue(residue), X=M_A, U=M_B, Z=M_C)
} 


PR_of_residue <- function(aa, aa_map) {
  M  <- matrix_for_class(aa)                 # M_A / M_B / M_C
  fk <- aa_map[[ ms_for_residue(aa) ]]       # class-median f_k
  pr <- as.numeric(fk %*% M)                 # row-vector times matrix
  pr / sum(pr)                               # safety normalize
}


peptide_PR <- function(seq, aa_map) {
  aas <- strsplit(seq, "")[[1]] #Split peptide sequence
  #Assign first probabilities to K to begin convolution
  K <- PR_of_residue(aas[1], aa_map)

  for (aa in aas[-1]) {
    # base::convolve does correlation by default, so we reverse b
    K <- convolve(K, rev( PR_of_residue(aa, aa_map) ), type = "open")
  }

  K / sum(K) }

peptide_PR("PEPTIDE", fk_map)
```

    ##  [1] 5.410905e-07 7.490890e-06 5.955788e-05 3.524864e-04 1.633084e-03
    ##  [6] 6.130325e-03 1.908604e-02 4.923393e-02 1.043528e-01 1.793653e-01
    ## [11] 2.399118e-01 2.289898e-01 1.348337e-01 3.604325e-02

``` r
mix_labeled_envelope <- function(env_nat, Kdist, z, delta_neutral = 2.004245) {
  delta_mz <- delta_neutral / z #m/z delta for each oxygen (dependent on charge)
  out <- NULL #Initialize empty object
  for (K in 0:(length(Kdist)-1)) { #kdist has length +1 because it includes 0
    shifted <- env_nat #Copy of starting natural positions
    
    #Shift m/z axis by K * (mass shift per label)
    shifted$mz <- shifted$mz + K * delta_mz 
    
    #Kdist[K+1] is P(K), the probability of having K heavy oxygens
    # ... We shifted k back for math, but need +1 for R indexing
    shifted$abundance <- shifted$abundance * Kdist[K+1]
    
    #Probabilities are rbinded to be able to aggregate
    out <- if (is.null(out)) shifted else rbind(out, shifted)
  }
  
  #Adds identical m/z and rescale so tallest peak is 100
  agg <- aggregate(abundance ~ mz, data = out, sum)
  agg$abundance <- 100 * agg$abundance / max(agg$abundance)
  agg
}


combine_envelopes <- function(start_envelope, label_envelope,
                              start_fraction, end_fraction) {
  
  start_pro <- start_envelope
  start_pro["abundance"] <- start_pro$abundance * start_fraction

  label_pro <- label_envelope
  label_pro["abundance"] <- label_pro$abundance * end_fraction  
    
  final_envelope <- rbind(start_pro, label_pro)
  final_envelope <- aggregate(abundance ~ mz, data = final_envelope, sum)
  
  final_envelope$abundance <- 100 * final_envelope$abundance / max(final_envelope$abundance)
  final_envelope  
}
```

``` r
peak_picker <- function(conv_df, start_mz, charge,
                        step_Da = 1.002123, win_ppm = 10) {
  step_mz = step_Da / charge
  
  # Ensure ascending by m/z
  conv_df <- conv_df[order(conv_df$mz), ]
  max_mz  <- max(conv_df$mz, na.rm = TRUE)

  picked_mz <- numeric(0)
  picked_I  <- numeric(0)

  seed <- start_mz

  # Step until run past the last m/z in the data
  while (seed <= max_mz) {
    # ±ppm window around the current target m/z ("seed")
    tol <- seed * win_ppm * 1e-6

    # Indices of points inside the window
    idx <- which(conv_df$mz >= (seed - tol) & conv_df$mz <= (seed + tol))

    if (length(idx) > 0) {
      # Take the apex (highest intensity) within the window
      sub_mz <- conv_df$mz[idx]
      sub_I  <- conv_df$abundance[idx]
      j      <- which.max(sub_I) 
  
      picked_mz <- c(picked_mz, sub_mz[j])
      picked_I  <- c(picked_I,  sub_I[j])
    }

    # Advance the target by the fixed step
    seed <- seed + step_mz
  }
  
  # Return exactly what’s in the data, no normalization
  data.frame(
    idx = paste0("M", seq_along(picked_mz) - 1),
    mz  = picked_mz,
    intensity = picked_I
  )
}
```

``` r
shift_right_accumulate <- function(v, p = 0.05) {
  n <- length(v)

  out <- v * (1 - p)         # everyone keeps (1-p) of themselves...
  if (n > 1) out[2:n] <- out[2:n] + p * v[1:(n-1)]  # ...and gets p from left
  out[n] <- v[n] + p * v[n-1]  # last bin keeps 100% of itself + p% from left

  out <- out / sum(out)
  out
}

shift_left_accumulate <- function(v, p = 0.05) {
  n <- length(v)
  
  out <- v * (1 - p)   # everyone keeps (1 - p) of themselves
  if (n > 1) {
    out[1:(n-1)] <- out[1:(n-1)] + p * v[2:n] }  # gets p from right
  out[1] <- v[1] + p * v[2] # first bin keeps 100% of itself + p from right

  out <- out / sum(out)
  out }

fk_optim_shift <- lapply(fk_map, shift_left_accumulate, p = 0.26)
```

``` r
library(parallel)
num_cores <- detectCores() - 1

#You should parallelize this function -> extremely slow and I'm done writing this script
envelope_fit_optim <- function(l_frac, ss_start_env, ss_lbl_env,
                               theo_mz, pep_z, obs_peaks) {
  heavy_synthesized <- 1 - l_frac

  #Final merged envelope based on optimizer l_frac
  merged_envelope <- combine_envelopes(ss_start_env, ss_lbl_env,
                                       l_frac, heavy_synthesized)
  
  #Simulating GFY's peak picking method
  late_shift_peaks <- peak_picker(merged_envelope, theo_mz, pep_z)
  late_shift_peaks["intensity"] <- late_shift_peaks$intensity / 
    max(late_shift_peaks$intensity) * 100 #Normalize to max 100

  #Combining predicted and observed data  
  l_peak_comparison <- merge(late_shift_peaks, obs_peaks,
                           by="idx", all.x=TRUE,
                           suffixes=c("_Pred", "_Obs")) %>% arrange(mz)
  
  theo_data <- as.numeric(l_peak_comparison$intensity_Pred)
  obs_data <- as.numeric(l_peak_comparison$intensity_Obs)

  rss <- (theo_data - obs_data)**2
  rss <- sum(rss)
  
  return(rss) }
```

``` r
estimate_synthesis <- function(row, pred_fk_map) {
  
  #------------ Start Prediction ------------------
  
  #Setting peptide sequence, charge, and theoretical m/z
  protein_id <- as.character(row["Protein.ID"])
  peptide_sequence <- row["Trimmed.Peptide"] %>% as.character
  pep_z <- row["z"] %>% as.numeric
  theo_mz <- row["Theo.m.z"] %>% as.numeric

  plot_title <- paste0("Peptide: ", peptide_sequence)

  #Calculating atomic formula and removing unnecessary items    
  named_atomic_vec <- get_light_formula(peptide_sequence, pep_z)
  named_atomic_vec <- named_atomic_vec[named_atomic_vec > 0] #Removing unneeded

  #Collapse into formula for enviPat::isopattern() from
  formula <- paste0(names(named_atomic_vec), named_atomic_vec,
                    collapse = "")

  #Calculate formula assuming infinitely possible resolution at ~ 1% threshold
  start_prediction <- suppressMessages( #Suppressing annoying message
    isopattern(isotopes, chemforms = formula,
               threshold = 0.01, charge = pep_z)[[1]] )
  start_prediction <- as.data.frame(start_prediction)
  
  #Convolving spectra assuming gaussian peaks
  start_envelope <- convolve_spectra(start_prediction, theo_mz,
                                     max(start_prediction$`m/z`))
  
  #------------ Label Prediction ------------------  

  #Calculating k probability distirbutions based on peptide sequence  
  K_label_shifted <- peptide_PR(peptide_sequence, pred_fk_map)
  
  #Shifting starting envelope based on weighted probabilities
  labeled_envelope <- mix_labeled_envelope(start_envelope, K_label_shifted,
                                           z = pep_z)

  #------------ Synthesis Prediction ------------------          

  #Data from file
  measured_DDA <- row[grepl("Max.Sn", names(row))]
  measured_DDA <- data.frame(idx=gsub("Max.Sn.", "", names(measured_DDA)),
                             intensity = as.numeric(measured_DDA))

  #Normalizing observed data such that maximum is 100 similar to predictions
  measured_DDA$intensity <- 100 * measured_DDA$intensity / max(measured_DDA$intensity)
  
  #This returns the fraction of light based on optimizing the observed data  
  synthesis_pred <- optimize(f=envelope_fit_optim, interval = c(0, 1),
                             ss_start_env = start_envelope,
                             ss_lbl_env = labeled_envelope,
                             theo_mz = theo_mz, pep_z = pep_z,
                             obs_peaks = measured_DDA)
  obj_error <- synthesis_pred$objective
  
  #------------ Visualization ------------------      
    
  light_remaining <- synthesis_pred$minimum
  heavy_synthesized <- 1 - light_remaining

  #Final merged envelope based on fitting
  merged_envelope <- combine_envelopes(start_envelope, labeled_envelope,
                                       light_remaining, heavy_synthesized)
  late_shift_peaks <- peak_picker(merged_envelope, theo_mz, pep_z)
  late_shift_peaks["intensity"] <- late_shift_peaks$intensity /
    max(late_shift_peaks$intensity) * 100

  #Combining predicted and observed data  
  l_peak_comparison <- merge(late_shift_peaks, measured_DDA,
                             by="idx", all.x=TRUE,
                             suffixes=c("_Pred", "_Obs")) %>% arrange(mz)  
  l_peak_comparison["idx"] <- as.numeric(gsub("M", "", l_peak_comparison$idx))
  
  peaks_signal <- sum(l_peak_comparison$intensity_Obs)
  
  outline <- data.frame(Protein_ID=protein_id, Peptide=peptide_sequence,
                        Signal=peaks_signal, Fit_Error = obj_error,
                        LightF = light_remaining, HeavyF = heavy_synthesized)  
  
  return(outline) }

#Setup cluster
cl <- makeCluster(num_cores)

#Export necessary variables/functions to the cluster nodes
clusterExport(cl,
              varlist=c("fk_optim_shift", "get_light_formula", "aa_composition",
                        "isotopes", "convolve_spectra", "peptide_PR",
                        "PR_of_residue", "matrix_for_class", "M_A", "M_B",
                        "M_C", "ms_for_residue", "mix_labeled_envelope",
                        "combine_envelopes", "envelope_fit_optim",
                        "peak_picker", "estimate_synthesis"))
invisible(clusterEvalQ(cl, library(enviPat)))
invisible(clusterEvalQ(cl, library(dplyr)))
invisible(clusterEvalQ(cl, library(ggplot2)))

norm_set_synthesis <- T5_O12_lblfree[T5_O12_lblfree$Protein.ID %in% external_norm_set,]
norm_set_synthesis <- parApply(cl, norm_set_synthesis, 1,
                               function(x){ estimate_synthesis(x, fk_optim_shift) })
norm_set_synthesis <- norm_set_synthesis %>% bind_rows
norm_set_synthesis <- norm_set_synthesis[!is.na(norm_set_synthesis$Signal),] %>%
  arrange(HeavyF)

write.csv(norm_set_synthesis, "Data/XLA_Norm/final_frog_synthesis_calc.csv",
          row.names = FALSE)

# Importing when not trying to rerun section above
# norm_set_synthesis <- read.csv("Data/XLA_Norm/050526_final_frog_synthesis_calc.csv")

stopCluster(cl) #Stop Cluster
```

``` r
p1 <- ggplot(norm_set_synthesis, aes(x = HeavyF)) +

  geom_histogram(
    aes(fill = after_stat(xmax <= .01)),
    binwidth = .01,
    boundary = 0,
    color = "black"
  ) +

  scale_fill_manual(
    values = c("TRUE" = "grey40",
               "FALSE" = "grey90"),
    guide = "none"
  ) +
  labs(x="Mixing coefficient of two-envelope model (f)", y="Peptides") +
  theme_bw() +
  theme(axis.text=element_text(size=21,colour="black"),
        plot.title = element_text(size=24, hjust=0.5),
        axis.title = element_text(size=22),
        panel.border = element_rect(linewidth=1.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
        legend.position = "none",
        plot.margin = margin(t = 10, r = 25, b = 10,
                             l = 10, unit = "pt"))
  
length(unique(norm_set_synthesis[norm_set_synthesis$HeavyF < 0.01,]$Protein_ID))
```

    ## [1] 179

``` r
# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=5, height=3, res=300)

p1
```

![](figures/Frog_GB_NYS-Proteomics/Synthesis%20fit%20distribution-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
visualize_predictions <- function(row, pred_fk_map) {
  
  #------------ Start Prediction ------------------
  
  #Setting peptide sequence, charge, and theoretical m/z
  protein_id <- as.character(row["Protein.ID"])
  peptide_sequence <- row["Trimmed.Peptide"] %>% as.character
  pep_z <- row["z"] %>% as.numeric
  theo_mz <- row["Theo.m.z"] %>% as.numeric

  plot_title <- paste0("Pep: ", peptide_sequence)

  #Calculating atomic formula and removing unnecessary items    
  named_atomic_vec <- get_light_formula(peptide_sequence, pep_z)
  named_atomic_vec <- named_atomic_vec[named_atomic_vec > 0] #Removing unneeded

  #Collapse into formula for enviPat::isopattern() from
  formula <- paste0(names(named_atomic_vec), named_atomic_vec,
                    collapse = "")

  #Calculate formula assuming infinitely possible resolution at ~ 1% threshold
  start_prediction <- suppressMessages( #Suppressing annoying message
    isopattern(isotopes, chemforms = formula,
               threshold = 0.01, charge = pep_z)[[1]] )
  start_prediction <- as.data.frame(start_prediction)
  
  #Convolving spectra assuming gaussian peaks
  start_envelope <- convolve_spectra(start_prediction, theo_mz,
                                     max(start_prediction$`m/z`))
  
  #------------ Label Prediction ------------------  

  #Calculating k probability distirbutions based on peptide sequence  
  K_label_shifted <- peptide_PR(peptide_sequence, pred_fk_map)
  
  #Shifting starting envelope based on weighted probabilities
  labeled_envelope <- mix_labeled_envelope(start_envelope, K_label_shifted,
                                           z = pep_z)

  #------------ Synthesis Prediction ------------------          

  #Data from file
  measured_DDA <- row[grepl("Max.Sn", names(row))]
  measured_DDA <- data.frame(idx=gsub("Max.Sn.", "", names(measured_DDA)),
                             intensity = as.numeric(measured_DDA))

  #Normalizing observed data such that maximum is 100 similar to predictions
  measured_DDA$intensity <- 100 * measured_DDA$intensity / max(measured_DDA$intensity)
  
  #This returns the fraction of light based on optimizing the observed data  
  synthesis_pred <- optimize(f=envelope_fit_optim, interval = c(0, 1),
                             ss_start_env = start_envelope,
                             ss_lbl_env = labeled_envelope,
                             theo_mz = theo_mz, pep_z = pep_z,
                             obs_peaks = measured_DDA)

  #------------ Visualization ------------------      
    
  light_remaining <- synthesis_pred$minimum
  heavy_synthesized <- 1 - light_remaining

  #Final merged envelope based on fitting
  merged_envelope <- combine_envelopes(start_envelope, labeled_envelope,
                                       light_remaining, heavy_synthesized)
  late_shift_peaks <- peak_picker(merged_envelope, theo_mz, pep_z)
  late_shift_peaks["intensity"] <- late_shift_peaks$intensity /
    max(late_shift_peaks$intensity) * 100

  #Combining predicted and observed data  
  l_peak_comparison <- merge(late_shift_peaks, measured_DDA,
                             by="idx", all.x=TRUE,
                             suffixes=c("_Pred", "_Obs")) %>% arrange(mz)  
  l_peak_comparison["idx"] <- as.numeric(gsub("M", "", l_peak_comparison$idx))
  
  plot_title <- paste0(plot_title, "<br>*f* = ", round(heavy_synthesized, 2))
  
  p1 <- ggplot() +
    geom_bar(data = l_peak_comparison, aes(x=idx, y=intensity_Obs),
             stat="identity", fill="#009E73", alpha=0.6, color="black") +
    geom_bar(data = l_peak_comparison, aes(x=idx, y=intensity_Pred),
             fill="#F0E442", alpha=0.6, stat="identity", color="black") +
    labs(y="Norm. abundance", x="MS1 peak", title=plot_title) +
    theme_bw() +
    theme(axis.text= element_text(size=20,colour="black"),
          plot.title = element_markdown(size = 18, hjust = 0.5),
          axis.title= element_text(size=20,colour="black"),          
          panel.border = element_rect(linewidth=1),
          panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
          panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
          legend.position = "none",
          plot.margin = margin(t = 10, r = 25, b = 10,
                               l = 10, unit = "pt")) +
    coord_cartesian(xlim=c(-1,51)) +
    scale_x_continuous(
      breaks = c(0, 10, 20, 30, 40, 50),
      labels = c("M0", "M10", "M20", "M30", "M40", "M50"))
    
  return(p1) }

example_df <- T5_O12_lblfree[T5_O12_lblfree$Trimmed.Peptide %in%
                                  c("GFGFVTFENVDDAK", "STNPGISIGDIAK",
                                    "IITMLPSSANAVEAYSGSNGILK",
                                    "VVEIAPAAQLDPQLR"),]

example_plots <- apply(example_df, 1,
                       function(x){ visualize_predictions(x, fk_optim_shift) })
example_plots
```

    ## [[1]]

![](figures/Frog_GB_NYS-Proteomics/Prediction%20visualization-1.png)<!-- -->

    ## 
    ## [[2]]

![](figures/Frog_GB_NYS-Proteomics/Prediction%20visualization-2.png)<!-- -->

    ## 
    ## [[3]]

![](figures/Frog_GB_NYS-Proteomics/Prediction%20visualization-3.png)<!-- -->

    ## 
    ## [[4]]

![](figures/Frog_GB_NYS-Proteomics/Prediction%20visualization-4.png)<!-- -->

``` r
plot_env_proteins <- function(protein_id) {
  
  plot_df <- GB_YolkNorm_Merge[GB_YolkNorm_Merge$Protein_ID == protein_id,]
  # plot_df[2:10] <- plot_df[2:10] / plot_df[[2]]

  plot_df <- plot_df %>% pivot_longer(cols = -Protein_ID)
  plot_df["Time"] <- rep((GB_timepoints+176)/60, length(set))
  
  p1 <- ggplot() +
    geom_line(data=plot_df, aes(x=Time, y=value, group=1),
              linewidth=3) +
    geom_hline(yintercept = 1, linetype = "dashed", size = 1) + 
    labs(y="Rel. abundance", x="Hours post-fertilization",
         title=protein_id) +
    theme_bw() +
    theme(axis.text= element_text(size=20,colour="black"),
          plot.title = element_text(size=18, vjust=0.5, hjust=0.5),
          axis.title= element_text(size=20,colour="black"),
          panel.border = element_rect(linewidth=1),
          aspect.ratio=1,
          panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
          panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
          legend.position = "none",
          plot.margin = margin(t = 10, r = 25, b = 10,
                               l = 10, unit = "pt")) +
    coord_cartesian(ylim=c(0,3))

  return(p1) }


# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=9, height=4.5, res=300)

plot_env_proteins(as.character(example_df[1,"Protein.ID"])) | example_plots[[1]]
```

![](figures/Frog_GB_NYS-Proteomics/Example%20Plots%20of%20some%20fitted%20envelopes%20for%20Supp.-1.png)<!-- -->

``` r
plot_env_proteins(as.character(example_df[3,"Protein.ID"])) | example_plots[[3]]
```

![](figures/Frog_GB_NYS-Proteomics/Example%20Plots%20of%20some%20fitted%20envelopes%20for%20Supp.-2.png)<!-- -->

``` r
plot_env_proteins(as.character(example_df[2,"Protein.ID"])) | example_plots[[2]]
```

![](figures/Frog_GB_NYS-Proteomics/Example%20Plots%20of%20some%20fitted%20envelopes%20for%20Supp.-3.png)<!-- -->

``` r
plot_env_proteins(as.character(example_df[4,"Protein.ID"])) | example_plots[[4]]
```

![](figures/Frog_GB_NYS-Proteomics/Example%20Plots%20of%20some%20fitted%20envelopes%20for%20Supp.-4.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
stable_prot_corr <- function(DNU_df, raw_df, ref_proteins) {

  flat_median <- DNU_df[DNU_df$Protein_ID %in% ref_proteins,]  
  flat_median <- apply(flat_median[2:10], 2, median)
  flat_median <- flat_median / mean(flat_median)

  # flat_median <- raw_df[raw_df$Protein_ID %in% ref_proteins,]  
  # flat_median <- apply(flat_median[3:11], 2, median)
  # flat_median <- flat_median / mean(flat_median)
  
  final_norm <- raw_df

  final_norm[3:11] <- sweep(final_norm[3:11], 2, flat_median, `/`)
  final_norm[3:11] <- sweep(final_norm[3:11], 1,
                           rowMeans(final_norm[3:11]), `/`)

  final_PQ <- protein_rollup(final_norm)

  return(final_PQ) }

stable_proteins <- 
  unique(norm_set_synthesis[norm_set_synthesis$HeavyF < 0.01,]$Protein_ID)

GB_R1_Final_PQ <- stable_prot_corr(DNU_Raw1_PQ, GB_R1_Raw, stable_proteins)
GB_R2_Final_PQ <- stable_prot_corr(DNU_Raw2_PQ, GB_R2_Raw, stable_proteins)
GB_R3_Final_PQ <- stable_prot_corr(DNU_Raw3_PQ, GB_R3_Raw, stable_proteins)
GB_R4_Final_PQ <- stable_prot_corr(DNU_Raw4_PQ, GB_R4_Raw, stable_proteins)

head(GB_R1_Final_PQ)
```

    ##                    Protein_ID        T0        T1        T2        T3        T4
    ## 1   XBmRNA10020|XBXL10_1g5370 0.9921728 1.0202522 1.0281373 1.0764107 0.9569899
    ## 2   XBmRNA10074|XBXL10_1g5397 1.0417323 0.9200835 0.9444619 0.7187058 0.6891544
    ## 3   XBmRNA10112|XBXL10_1g5418 0.9789732 1.0495653 0.9676166 0.9643504 1.0382850
    ## 4 XBmRNA101811|XBXL10_1g44445 0.8876647 0.8593678 1.1206637 1.0103566 1.1358560
    ## 5   XBmRNA10181|XBXL10_1g5453 1.0623872 0.9596192 1.1338402 0.9856879 0.9308265
    ## 6 XBmRNA101826|XBXL10_1g44449 1.2640597 1.0327780 0.9302904 0.9928999 0.9367468
    ##          T5        T6        T7        T8
    ## 1 0.9290970 1.0825569 0.9494448 0.9649383
    ## 2 0.8344608 1.1607373 1.0232237 1.6674403
    ## 3 1.0335962 1.0147447 0.9541592 0.9987094
    ## 4 0.9643364 0.9298672 1.0589249 1.0329626
    ## 5 1.0653854 1.0026202 0.8891509 0.9704826
    ## 6 0.8584499 0.8488822 1.0166125 1.1192806

``` r
#Calculating median data for plot
GB_FinalNorm_Merge <- rbind(GB_R1_Final_PQ, GB_R2_Final_PQ,
                            GB_R3_Final_PQ, GB_R4_Final_PQ) %>%
  group_by(Protein_ID) %>%
  dplyr::summarise(across(where(is.numeric), \(x) median(x)),
                   .groups = "drop")

GB_FinalNorm_Merge[paste0("T", seq(0,8))] <- 
  sweep(GB_FinalNorm_Merge[paste0("T", seq(0,8))], 1,
        rowMeans(GB_FinalNorm_Merge[paste0("T", seq(0,8))]), `/`)  

head(GB_FinalNorm_Merge)
```

    ## # A tibble: 6 × 10
    ##   Protein_ID                  T0    T1    T2    T3    T4    T5    T6    T7    T8
    ##   <chr>                    <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
    ## 1 XBmRNA10020|XBXL10_1g53… 1.00  1.03  0.975 1.08  0.950 1.01  1.01  0.948 0.995
    ## 2 XBmRNA10074|XBXL10_1g53… 0.942 0.796 0.821 0.762 0.848 0.812 1.11  1.29  1.62 
    ## 3 XBmRNA10112|XBXL10_1g54… 0.996 0.975 0.989 0.998 1.01  1.02  1.02  0.985 1.01 
    ## 4 XBmRNA10135|XBXL10_1g54… 1.05  1.18  0.971 1.07  0.940 0.815 0.981 0.973 1.02 
    ## 5 XBmRNA1017|XBXL10_1g750  0.887 0.997 0.945 0.949 1.05  0.941 1.01  1.11  1.11 
    ## 6 XBmRNA101811|XBXL10_1g4… 1.02  0.977 1.02  0.997 0.954 1.06  1.00  1.01  0.961

``` r
plot_proteins <- function(set, data) {
  
  name_df <- data.frame(Protein_ID = set, XLA_Gene=names(set))
  
  plot_df <- data[data$Protein_ID %in% set,]
  plot_df[2:10] <- plot_df[2:10] / plot_df[[2]]

  plot_df <- plot_df %>% pivot_longer(cols = -Protein_ID)
  plot_df["Time"] <- rep((GB_timepoints+176)/60, length(set))
  
  plot_df <- merge(name_df, plot_df, by="Protein_ID")
  
  p1 <- ggplot() +
    geom_line(data=plot_df, aes(x=Time, y=value, group=XLA_Gene,
                                color=XLA_Gene),
              linewidth=1.5) +
    geom_point(data=plot_df, aes(x=Time, y=value, color=XLA_Gene),
               size=3) +
    geom_hline(yintercept = 1, linetype = "dashed", size = 1) +    
    theme_bw() +
    labs(x = "Hours post-fertilization", y = "Rel. Abundance") +
    theme(axis.text = element_text(size = 18, colour = "black"),
          plot.title = element_text(size = 18, hjust = 0.5),
          axis.title = element_text(size = 22),
          aspect.ratio = 1,
          legend.text = element_text(size=16),
          legend.title = element_text(size=16),          
          panel.border = element_rect(linewidth = 2),
          panel.grid.major = element_line(color = "grey70", linewidth = 0.25),
          panel.grid.minor = element_line(color = "grey70", linewidth = 0.25)) +
    coord_cartesian(ylim=c(0,1.25))

  return(p1) }

plot_proteins(frog_yolk_set, GB_AllReps_DNU)
```

![](figures/Frog_GB_NYS-Proteomics/Visualizing%20effects%20of%20different%20normalizations-1.png)<!-- -->

``` r
plot_proteins(frog_yolk_set, GB_AllReps_WI)
```

![](figures/Frog_GB_NYS-Proteomics/Visualizing%20effects%20of%20different%20normalizations-2.png)<!-- -->

``` r
plot_proteins(frog_yolk_set, GB_YolkNorm_Merge)
```

![](figures/Frog_GB_NYS-Proteomics/Visualizing%20effects%20of%20different%20normalizations-3.png)<!-- -->

``` r
plot_proteins(frog_yolk_set, GB_FinalNorm_Merge)
```

![](figures/Frog_GB_NYS-Proteomics/Visualizing%20effects%20of%20different%20normalizations-4.png)<!-- -->

``` r
# The normalization reference set read by the model-fitting documents is provided as
# Data/XLA_Norm/XLA-O18_YolkNormSet_T8-NYS-Decay.csv. Re-running this chunk writes the
# candidate list to a separate file.
write.csv(data.frame(Protein_ID=stable_proteins),
          "Data/XLA_Norm/XLA-O18_YolkNormSet_T8-NYS-Decay_candidates.csv",
          row.names=FALSE)

write.csv(GB_YolkNorm_Merge,
          "Data/XLA_Norm/XLA-O18_GB-FinalDataset.csv",
          row.names=FALSE)
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
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] ggtext_0.1.2      enviPat_2.8       minpack.lm_1.2-4  TOSTER_0.8.6     
    ##  [5] patchwork_1.3.2   ggplot2_4.0.3     data.table_1.18.4 stringr_1.6.0    
    ##  [9] dplyr_1.2.0       tidyr_1.3.2      
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] ggdist_3.3.3         utf8_1.2.6           generics_0.1.4      
    ##  [4] xml2_1.5.2           stringi_1.8.7        digest_0.6.39       
    ##  [7] magrittr_2.0.4       evaluate_1.0.5       grid_4.5.3          
    ## [10] RColorBrewer_1.1-3   fastmap_1.2.0        purrr_1.2.2         
    ## [13] scales_1.4.0         cli_3.6.5            rlang_1.1.7         
    ## [16] litedown_0.10        commonmark_2.0.0     cowplot_1.2.0       
    ## [19] withr_3.0.3          yaml_2.3.12          otel_0.2.0          
    ## [22] tools_4.5.3          vctrs_0.7.1          R6_2.6.1            
    ## [25] lifecycle_1.0.5      pkgconfig_2.0.3      pillar_1.11.1       
    ## [28] gtable_0.3.6         glue_1.8.0           Rcpp_1.1.1-1.1      
    ## [31] xfun_0.58            tibble_3.3.1         tidyselect_1.2.1    
    ## [34] rstudioapi_0.19.0    knitr_1.51           farver_2.1.2        
    ## [37] htmltools_0.5.9      rmarkdown_2.31       labeling_0.4.3      
    ## [40] compiler_4.5.3       S7_0.2.1             markdown_2.0        
    ## [43] distributional_0.7.1 gridtext_0.1.6
