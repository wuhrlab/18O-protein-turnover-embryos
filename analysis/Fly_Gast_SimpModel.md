Fly Gast kD Fits - Simplified
================
Edward Cruz
2026-06-09

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
library(stringr)
library(ggplot2)
library(purrr)
library(patchwork)
library(minpack.lm)
library(parallel)

knitr::opts_chunk$set(fig.path = "figures/Fly_Gast_SimpModel/")
num_cores <- detectCores() - 1  # Save one core for system stability
```

Minimum HPF: 117 min Max HPF: 257 min

Time after labeling start is the same -\> Average start is 3 hrs 7 min

``` r
sn.cols <- c("126 Sn", "127n Sn", "127c Sn", "128n Sn", "128c Sn", "129n Sn",
             "129c Sn", "130n Sn", "130c Sn", "131n Sn", "131c Sn", "132n Sn",
             "132c Sn", "133n Sn", "133c Sn", "134n Sn", "134c Sn", "135n Sn")

label_order <- c("T0A", "T0B", "O1", "O2", "O3", "O4", "O5", "O6", "O7",
                 "O8", "O9", "O10", "O11", "O12", "N13", "N14", "N15", "N16")

raw_heavy_labels <- c("T0A", "T0B", label_order[c(3:14)])
raw_light_labels <- c("T0A", "T0B", label_order[c(15:18)])

#Files to normalize
Fly_O18_T6  <- read.csv("Data/Dmel_O18/Fly_O18-T6_Filtered_Data.csv")
Fly_O18_T7  <- read.csv("Data/Dmel_O18/Fly_O18-T7_Filtered_Data.csv")

#Ultimately - collected these minutes apart so used one vector
Fly_T6.time.min <- c(257, 287, 320, 351, 386, 441, 500,
                     565, 682, 868, 995, 1218, 1514) - 257
light_time <- c(0, Fly_T6.time.min[c(4, 8, 11, 13)])

Fly_T6_hpf_time <- ((Fly_T6.time.min + 257) + (Fly_T6.time.min + 117))/2
Fly_T6_hpf_time <- Fly_T6_hpf_time / 60
Fly_T6_hpf_light_time <- c(Fly_T6_hpf_time[1], Fly_T6_hpf_time[c(4, 8, 11, 13)])

colorblind_palette <- c("#CC79A7", "#D55E00", "#0072B2", "#F0E442",
                        "#009E73", "#56B4E9", "#E69F00", "#999999",
                        "#5E4FA2")
```

``` r
Dmel.genes <- read.csv("Files/Reference/Dmel_GeneTable_Descriptions.csv")
protein_annotation <- Dmel.genes[c("Uniprot_ID", "Name", "Description")]
colnames(protein_annotation)[1] <- "Protein_ID"

Dmel.genes <- Dmel.genes[c("Uniprot_ID", "Gene_Symbol")]
colnames(Dmel.genes)[1] <- "Protein_ID"
```

``` r
Fly_O18_T6 <- Fly_O18_T6[!rowSums(Fly_O18_T6[raw_heavy_labels]) == 0,]
Fly_O18_T7 <- Fly_O18_T7[!rowSums(Fly_O18_T7[raw_heavy_labels]) == 0,]

Fly_O18_T6["StartPool"] <- apply(Fly_O18_T6[raw_heavy_labels], 1,
      function(row){
        starting_pool <- sum(row[c("T0A", "T0B", "O1")]) / sum(row)
        return(starting_pool) })

Fly_O18_T6["HeavyWeight"] <- apply(Fly_O18_T6[c(raw_light_labels, raw_heavy_labels)], 1,
      function(row){

        light_sn <- row[1:length(raw_light_labels)]
        heavy_sn <- row[(length(raw_light_labels)+1):length(row)]
        
        return(mean(light_sn)/(mean(heavy_sn) + 1e-6)) })

Fly_O18_T7["StartPool"] <- apply(Fly_O18_T7[raw_heavy_labels], 1,
      function(row){
        starting_pool <- sum(row[c("T0A", "T0B", "O1")]) / sum(row)
        return(starting_pool) })

Fly_O18_T7["HeavyWeight"] <- apply(Fly_O18_T7[c(raw_light_labels, raw_heavy_labels)], 1,
      function(row){

        light_sn <- row[1:length(raw_light_labels)]
        heavy_sn <- row[(length(raw_light_labels)+1):length(row)]
        
        return(mean(light_sn)/(mean(heavy_sn) + 1e-6)) })


#This sets a floor to quantify based on signal
# -> Requires (T0 + O1)/14 > 0.15
# -> Requires heavy signal to not outweigh light signal
Fly_O18_T6 <- Fly_O18_T6[(Fly_O18_T6$StartPool > 0.15) &
                           (Fly_O18_T6$HeavyWeight > 0.8),]
Fly_O18_T7 <- Fly_O18_T7[(Fly_O18_T7$StartPool > 0.15) &
                           (Fly_O18_T7$HeavyWeight > 0.8),]
```

``` r
#Splitting light vs heavy experiments
split_experiments <- function(df){
  base_df <- df[c("Protein_ID", "Peptide", "Pep_k1", "Pep_k2")]
  
  heavy_df <- df[raw_heavy_labels]
  heavy_df <- heavy_df/rowSums(heavy_df)
  colnames(heavy_df)[1:2] <- c("T0A_H", "T0B_H")
  
  light_df <- df[raw_light_labels]
  light_df <- light_df/rowSums(light_df)
  colnames(light_df)[1:2] <- c("T0A_L", "T0B_L")
    
  split_df <- cbind(base_df, heavy_df, light_df,
                    df[c("sum_sn", "StartPool", "HeavyWeight")])  
  
  return(split_df) }

Fly_O18_T6 <- split_experiments(Fly_O18_T6)
Fly_O18_T7 <- split_experiments(Fly_O18_T7)

raw_heavy_labels <- c("T0A_H", "T0B_H", label_order[c(3:14)])
raw_light_labels <- c("T0A_L", "T0B_L", label_order[c(15:18)])

heavy_labels <- c("T0_H", label_order[c(3:14)])
light_labels <- c("T0_L", label_order[c(15:18)])
```

``` r
#'*Wide Isolation Norm Dmel O18 T6*
WideIso_T6_vec <- read.csv(
  "Data/Dmel_O18/Norm_WideIso/ORC_04724_T6-EdNew_WideIso-Norm.csv")[
    paste("X", gsub(" ", ".", sn.cols), sep="")]
WideIso_T6_vec <- apply(WideIso_T6_vec, 2, median)

WideIso_T6_H <- WideIso_T6_vec[c(1:14)]
WideIso_T6_H <- WideIso_T6_H / sum(WideIso_T6_H)

WideIso_T6_L <- WideIso_T6_vec[c(1:2, 15:18)]
WideIso_T6_L <- WideIso_T6_L / sum(WideIso_T6_L)

#'*Wide Isolation Norm Dmel O18 T7*
WideIso_T7_vec <- read.csv(
  "Data/Dmel_O18/Norm_WideIso/ORC_04730_T7-EdNew_WideIso-Norm.csv")[
    paste("X", gsub(" ", ".", sn.cols), sep="")]
WideIso_T7_vec <- apply(WideIso_T7_vec, 2, median)

WideIso_T7_H <- WideIso_T7_vec[c(1:14)]
WideIso_T7_H <- WideIso_T7_H / sum(WideIso_T7_H)

WideIso_T7_L <- WideIso_T7_vec[c(1:2, 15:18)]
WideIso_T7_L <- WideIso_T7_L / sum(WideIso_T7_L)
```

``` r
norm_Fly_split <- function(df, H_error, L_error) {
  
  #Divide each original ratio by the pipet error
  # ... Divide each column by the row sum
  H_ratio.df <- as.matrix(df[raw_heavy_labels])
  H_corr.ratios <- sweep(H_ratio.df, 2, H_error, `/`)
  H_corr.ratios <- sweep(H_corr.ratios, 1, rowMeans(H_corr.ratios), `/`)
  
  #Divide each original ratio by the pipet error
  # ... Divide each column by the row sum
  L_ratio.df <- as.matrix(df[raw_light_labels])
  L_corr.ratios <- sweep(L_ratio.df, 2, L_error, `/`)
  L_corr.ratios <- sweep(L_corr.ratios, 1, rowMeans(L_corr.ratios), `/`)
  
  #Normalized dataframe then shifting to 1 -> (1/18 Ratio) = 1    
  norm.df <- df
  norm.df[raw_heavy_labels] <- H_corr.ratios
  norm.df[raw_light_labels] <- L_corr.ratios
  
  return(as.data.frame(norm.df)) }

heavy_labels <- c("T0_H", label_order[c(3:14)])
light_labels <- c("T0_L", label_order[c(15:18)])

WI_T6_df <- norm_Fly_split(Fly_O18_T6, WideIso_T6_H, WideIso_T6_L)
WI_T7_df <- norm_Fly_split(Fly_O18_T7, WideIso_T7_H, WideIso_T7_L)
```

``` r
fit_yolk_trend <- function(df1, df2) {
  
  quant_protein <- function(x, target_df){

    pq <- target_df[target_df$Protein_ID == x,]
    pq$Fraction <- pq$sum_sn / sum(pq$sum_sn)
    pq[c(raw_heavy_labels, raw_light_labels)] <- 
      pq[c(raw_heavy_labels, raw_light_labels)] * pq$Fraction
    
    rollup <- apply(pq[c(raw_heavy_labels, raw_light_labels)], 2, sum)
    
    light_ch <- rollup[raw_light_labels] / mean(rollup[raw_light_labels])
    heavy_ch <- rollup[raw_heavy_labels] / mean(rollup[raw_heavy_labels])

    return(data.frame(t(c(light_ch, heavy_ch)))) }

  yolk_sub1 <- df1[df1$Protein_ID %in% yolk.proteins,]  
  yolk_sub1 <- lapply(unique(yolk_sub1$Protein_ID), function(x) {
    quant_protein(x, df1)}) %>% bind_rows
  yolk_sub1 <- apply(yolk_sub1, 2, median)
  
  yolk_sub2 <- df2[df2$Protein_ID %in% yolk.proteins,]  
  yolk_sub2 <- lapply(unique(yolk_sub2$Protein_ID), function(x) {
    quant_protein(x, df2)}) %>% bind_rows
  yolk_sub2 <- apply(yolk_sub2, 2, median)
  
  lch_1 <- yolk_sub1[raw_light_labels]
  hch_1 <- yolk_sub1[raw_heavy_labels]

  lch_2 <- yolk_sub2[raw_light_labels]
  hch_2 <- yolk_sub2[raw_heavy_labels]

  norm_vector <- function(vec) {
    norm_t <- vec / mean(vec)
    return(norm_t) }

  all_data <- data.frame(Time=c(rep(c(0, light_time), 2),
                                rep(c(0, Fly_T6.time.min), 2)),
                         Values=c(norm_vector(lch_1),
                                  norm_vector(lch_2),
                                  norm_vector(hch_1),
                                  norm_vector(hch_2)))  

  lin <- lm(log(Values) ~ Time, data = all_data, weights = all_data$Values^2)
  a0 <- unname(exp(coef(lin)[1]))
  b0 <- unname(-coef(lin)[2])  

  yolk_fit_nls <- nlsLM(Values ~ a * exp(-b * Time),
                        data = all_data,
                        start = list(a = a0, b = b0),
                        control = list(maxiter = 1000))

  #-----------------------------------------------------------------  
  # Plots are mostly shown with fixed start for simplification
  
  plot_vector <- function(vec) {

    norm_t <- c(mean(vec[1:2]), vec[3:length(vec)])
    norm_t <- norm_t / mean(norm_t)
    norm_t <- norm_t / norm_t[1]
    
    return(norm_t) }
  
  plot_data <- data.frame(Time=c(rep(light_time/60, 2),
                                 rep(Fly_T6.time.min/60, 2)),
                          Values=c(plot_vector(lch_1), plot_vector(lch_2),
                                   plot_vector(hch_1), plot_vector(hch_2)),
                          Set=c(rep("LCH1", length(light_time)),
                                rep("LCH2", length(light_time)),
                                rep("HCH1", length(Fly_T6.time.min)),
                                rep("HCH2", length(Fly_T6.time.min))))  

  median_WI <- aggregate(Values ~ Time, data = plot_data, FUN = median)
  median_WI$Values <- median_WI$Values / mean(median_WI$Values)
  
  smooth_fit <-
    data.frame(Time = Fly_T6.time.min)
  smooth_fit["Predict"] <- exp_decay(smooth_fit$Time,
                                     coef(yolk_fit_nls)['a'],
                                     coef(yolk_fit_nls)['b'])
  smooth_fit$Time <- smooth_fit$Time/60
  smooth_fit$Predict <- smooth_fit$Predict / smooth_fit$Predict[1]
  
  p1 <- ggplot() +
    geom_line(data=plot_data, aes(x=Time+(187/60), y=Values, group=Set),
              color="grey60", alpha=0.8, linewidth=2) +
    geom_line(data=smooth_fit, aes(x=Time+(187/60), y=Predict), linewidth=3, alpha=0.6,
              color="#D55E00") +
    geom_hline(yintercept=1, linewidth=1, linetype="dashed", color="black") +    
    labs(x = "Hours post-fertilization", y = "Rel. Abundance") +
    theme_bw() +
    theme(axis.text = element_text(size = 22, colour = "black"),
          plot.title = element_text(size = 20, hjust = 0.5),
          axis.title = element_text(size = 22),
          aspect.ratio = 1,
          legend.text = element_text(size=16),
          legend.title = element_text(size=16),
          panel.border = element_rect(linewidth = 2),
          panel.grid.major = element_line(color = "grey70", linewidth = 0.25),
          panel.grid.minor = element_line(color = "grey70", linewidth = 0.25),
          plot.margin = margin(t = 10, r = 25, b = 10,
                               l = 10, unit = "pt")) +
    coord_cartesian(ylim=c(0,1.1)) +
    scale_y_continuous(breaks=seq(0,1.25,0.25))

  print(p1)
  
  return(coef(yolk_fit_nls)) }


exp_decay <- function(t_hrs, a, b) {
  return(a * exp(-b * t_hrs)) }

yolk.proteins <- c("P02843", "P02844", "P06607") #Yp1, Yp2, Yp3
names(yolk.proteins) <- c("Yp1", "Yp2", "Yp3")

# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=4.5, height=4.5, res=300)

yolk_fit <- fit_yolk_trend(WI_T6_df, WI_T7_df)
```

![](figures/Fly_Gast_SimpModel/Fitting%20yolk%20trend%20to%20YPs-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
theory_heavy <- exp_decay(c(0, Fly_T6.time.min), yolk_fit['a'],
                          yolk_fit['b'])
theory_heavy <-  theory_heavy / mean(theory_heavy)

theory_light <- exp_decay(c(0, light_time), yolk_fit['a'],
                          yolk_fit['b'])
theory_light <- theory_light / mean(theory_light)
```

``` r
corr_by_theory <- function(df, ch, theory) {
  
  sub_df <- df[df$Protein_ID %in% yolk.proteins, c(ch, "sum_sn")]
  sub_df[ch] <- sub_df[ch] / rowMeans(sub_df[ch])

  sub_df["Fraction"] <- sub_df$sum_sn / sum(sub_df$sum_sn)
  sub_df <- sub_df[ch] * sub_df$Fraction
  
  yolk_sum <- apply(sub_df, 2, sum)
  yolk_sum <- yolk_sum / mean(yolk_sum)  

  corr_factor <- yolk_sum / theory  

  return(corr_factor) }
```

``` r
T6_light_corr <- corr_by_theory(Fly_O18_T6, raw_light_labels, theory_light)
T7_light_corr <- corr_by_theory(Fly_O18_T7, raw_light_labels, theory_light)

T6_heavy_corr <- corr_by_theory(Fly_O18_T6, raw_heavy_labels, theory_heavy)
T7_heavy_corr <- corr_by_theory(Fly_O18_T7, raw_heavy_labels, theory_heavy)

#--------------------------------------

norm_T6.df <- norm_Fly_split(Fly_O18_T6, T6_heavy_corr, T6_light_corr)
norm_T6.df <- norm_T6.df[(norm_T6.df$T0A_H > 0.8) & (norm_T6.df$T0B_H > 0.8),]
norm_T6.df <- norm_T6.df[(norm_T6.df$T0A_L > 0.15) & (norm_T6.df$T0B_L > 0.15),]

norm_T7.df <- norm_Fly_split(Fly_O18_T7, T7_heavy_corr, T7_light_corr)
norm_T7.df <- norm_T7.df[(norm_T7.df$T0A_H > 0.8) & (norm_T7.df$T0B_H > 0.8),]
norm_T7.df <- norm_T7.df[(norm_T7.df$T0A_L > 0.15) & (norm_T7.df$T0B_L > 0.15),]

head(norm_T6.df)
```

    ##   Protein_ID                Peptide Pep_k1      Pep_k2    T0A_H    T0B_H
    ## 1 A0A023GRW3 DYTAVATPQPILSDTEAATTSR      1 0.019371857 1.303900 1.270988
    ## 2 A0A023GRW3     GLDALTNEEGVLTALQQR      1 0.019832737 1.265698 1.317057
    ## 3 A0A023GRW3        LNDYVPEAGPPAISK      1 0.014837079 1.270216 1.362397
    ## 4 A0A023GRW3              VDLEPTSIR      1 0.012141327 1.343749 1.349311
    ## 5 A0A0A1EI90               IGLLDNAK      1 0.009538950 1.064771 1.778539
    ## 6 A0A0B4JCU3               FVHPMLTR      1 0.009296024 1.212077 1.397410
    ##         O1       O2       O3       O4       O5        O6        O7        O8
    ## 1 1.510216 1.201291 1.302594 1.220548 1.084274 1.0166657 1.1029814 0.8208088
    ## 2 1.128195 1.295572 1.180233 1.307920 1.217075 1.3838368 1.0263448 0.8139898
    ## 3 1.499694 1.388829 1.100761 1.171800 1.191816 0.9716466 1.0379431 0.8687891
    ## 4 1.318210 1.522413 1.174517 1.293097 1.216326 1.0853129 0.9292052 0.8010225
    ## 5 1.211027 1.074438 1.437891 1.703253 1.544567 1.0642907 0.5367249 0.4558719
    ## 6 1.027079 1.289248 1.378472 1.436120 1.205077 1.1999782 0.9092770 0.7111652
    ##          O9       O10       O11       O12     T0A_L    T0B_L       N13
    ## 1 0.6877578 0.6191191 0.4683843 0.3904710 1.0268422 1.001644 1.0326742
    ## 2 0.7498945 0.5868968 0.4747499 0.2525365 0.9892746 1.030158 1.1507062
    ## 3 0.7539540 0.5418534 0.4733439 0.3669569 0.9944799 1.067418 0.9333578
    ## 4 0.7250322 0.5284877 0.3709620 0.3423541 0.9954990 1.000340 0.9251612
    ## 5 0.6574010 0.5849428 0.5620504 0.3242326 0.8031171 1.342452 0.5979955
    ## 6 0.6507393 0.5807409 0.5659068 0.4367095 1.0673453 1.231435 1.0189032
    ##         N14       N15       N16    sum_sn StartPool HeavyWeight
    ## 1 1.1463651 0.9318308 0.8606436  927.5658 0.2853168    1.212107
    ## 2 0.9240675 0.9290608 0.9767327  702.6622 0.2612201    1.224303
    ## 3 1.0856087 0.9315683 0.9875669 1742.4237 0.2882710    1.209673
    ## 4 0.9538238 1.1014163 1.0237600 1638.2761 0.2819314    1.273946
    ## 5 0.5469141 1.4274735 1.2820479 1920.9276 0.2867982    1.223964
    ## 6 1.2116081 0.7321607 0.7385477 6752.5860 0.2559441    1.100363

``` r
head(norm_T7.df)
```

    ##   Protein_ID                Peptide Pep_k1      Pep_k2     T0A_H    T0B_H
    ## 1 A0A021WW64             TIFAVGSFLR      1 0.010841243 0.9212154 1.130516
    ## 2 A0A023GRW3              ALEASVGTK      1 0.010926241 1.1505976 1.467493
    ## 3 A0A023GRW3 DYTAVATPQPILSDTEAATTSR      1 0.019371857 1.1732288 1.189437
    ## 4 A0A023GRW3        LNDYVPEAGPPAISK      1 0.014837079 1.1001393 1.262364
    ## 5 A0A023GRW3               YDTPDVSK      1 0.008920047 1.3776451 1.354358
    ## 6 A0A0B4JCU3               FVHPMLTR      1 0.009296024 1.4831847 1.282635
    ##         O1       O2       O3       O4       O5       O6        O7        O8
    ## 1 1.202098 1.377231 1.358160 1.265195 1.063152 1.137159 0.7150183 0.8309828
    ## 2 1.275488 1.306297 1.383087 1.208066 1.174868 1.017244 0.8076463 0.6919215
    ## 3 1.365236 1.221310 1.305683 1.197688 1.147293 1.194384 1.1683207 0.9329174
    ## 4 1.337477 1.385796 1.303549 1.195905 1.095942 1.147555 0.9785424 1.0077188
    ## 5 1.395159 1.559796 1.496317 1.045942 1.230764 0.972694 0.6359049 0.4849814
    ## 6 1.053881 1.469697 1.378617 1.351864 1.115459 1.065355 0.9221624 0.6395322
    ##          O9       O10       O11       O12     T0A_L    T0B_L       N13      N14
    ## 1 0.7084855 0.8545672 0.7434766 0.6927439 0.8274319 1.009792 0.8893860 1.015585
    ## 2 0.7548576 0.6267256 0.5205432 0.6151657 1.0062061 1.276215 0.8725061 1.004478
    ## 3 0.6939246 0.6107283 0.4487766 0.3510713 1.1031616 1.112198 1.0332460 1.105487
    ## 4 0.7378945 0.6353862 0.4728578 0.3388724 1.0571271 1.206280 1.0373294 1.056652
    ## 5 0.9429039 0.4278131 0.6208550 0.4548663 1.2690023 1.240631 0.8242670 0.924551
    ## 6 0.6520420 0.5706565 0.5252368 0.4896774 1.3391236 1.151629 1.1260500 1.008496
    ##         N15       N16    sum_sn StartPool HeavyWeight
    ## 1 1.1531000 1.1047047  285.8936 0.2152183    1.180488
    ## 2 0.8376930 1.0029019 1799.4349 0.2613783    1.212682
    ## 3 0.8516663 0.7942409 1169.1872 0.2492197    1.128819
    ## 4 0.9097850 0.7328262 2723.7827 0.2474700    1.104694
    ## 5 0.8883993 0.8531492  281.9851 0.2809648    1.162342
    ## 6 0.7181354 0.6565656 6182.6616 0.2589954    1.173099

``` r
start_check <- lapply(intersect(norm_T6.df$Protein_ID, norm_T7.df$Protein_ID),
                      function(protein_id){

  determine_start <- function(x, exp_df) {
    
    sub_df <- exp_df[exp_df$Protein_ID == x,]
    sub_df["T0_H"] <- (sub_df$T0A_H + sub_df$T0B_H)/2    
    sub_df["T0_L"] <- (sub_df$T0A_L + sub_df$T0B_L)/2        
    sub_df["Fraction"] <- sub_df$sum_sn / sum(sub_df$sum_sn)  
    
    sub_df <- sub_df[c(heavy_labels, light_labels, "Fraction")]
    sub_df[heavy_labels] <- sweep(sub_df[heavy_labels], 1,
                                  rowMeans(sub_df[heavy_labels]), `/`)
    sub_df[light_labels] <- sweep(sub_df[light_labels], 1,
                                  rowMeans(sub_df[light_labels]), `/`)
    
    h_rollup <- sub_df[heavy_labels] * sub_df$Fraction
    h_rollup <- apply(h_rollup, 2, sum)
    h_rollup <- h_rollup / mean(h_rollup)

    l_rollup <- sub_df[light_labels] * sub_df$Fraction
    l_rollup <- apply(l_rollup, 2, sum)
    l_rollup <- l_rollup / mean(l_rollup)

    return(c(h_rollup, l_rollup)) }
                        
  T6_med <- determine_start(protein_id, norm_T6.df)
  T7_med <- determine_start(protein_id, norm_T7.df)

  outline <- data.frame(Protein_ID = protein_id,
                        T6_T0_H = T6_med["T0_H"],
                        T7_T0_H = T7_med["T0_H"],
                        T0_H_FC = T6_med["T0_H"] / T7_med["T0_H"],                        
                        T6_T0_L = T6_med["T0_L"],
                        T7_T0_L = T7_med["T0_L"],
                        T0_L_FC = T6_med["T0_L"] / T7_med["T0_L"])

  return(outline) }) %>% bind_rows

start_check["Keep"] <- rep(0)
start_check[abs(log(start_check$T0_H_FC)) < log(1.5) &
            abs(log(start_check$T0_L_FC)) < log(1.5), "Keep"] <- 1
remove_proteins <- start_check[start_check$Keep==0,]$Protein_ID

norm_T6.df <- norm_T6.df[!norm_T6.df$Protein_ID %in% remove_proteins,]
norm_T7.df <- norm_T7.df[!norm_T7.df$Protein_ID %in% remove_proteins,]
```

``` r
# Determing where synthesis isn't relevant
pep_decay_rate <- function(fraction, k1, k2) {
  (-1 * log(fraction / k1)) / k2 }

#-----------------------------------------------------------------------

AA_peptide_cutoff <- function(row, l_time) {

  pep_k1 <- as.numeric(row["Pep_k1"])
  pep_k2 <- as.numeric(row["Pep_k2"])

  #Trimming data to where below 10%
  time_cutoff <- pep_decay_rate(0.1, pep_k1, pep_k2)
  light_trim <- l_time[l_time < time_cutoff]

  #If depletes too fast, then floor to minimum of 2
  if (length(light_trim) < 2) { light_trim <- l_time[1:2] }

  #Calculating real target decrease of AA probabilities
  f_synth_calc <- exp(-pep_k2 * tail(light_trim, 1))
  
  # Extend by one timepoint while f_synth is still above 10% 
  while (f_synth_calc > 0.1 && length(light_trim) < length(l_time)) {
    next_idx <- length(light_trim) + 1
    light_trim <- l_time[1:next_idx]
    f_synth_calc <- exp(-pep_k2 * tail(light_trim, 1)) }  

  n_timepoints <- length(light_trim)

  out_df <- data.frame(Peptide=row[1],
                       Final_Prob=f_synth_calc,
                       N_timepoints=n_timepoints)

  return(out_df) }

AA_pep_T6 <- norm_T6.df[c("Peptide", "Pep_k1", "Pep_k2")]
AA_pep_T6 <- apply(AA_pep_T6, 1, function(x) {
  AA_peptide_cutoff(x, light_time)}) %>% bind_rows
row.names(AA_pep_T6) <- NULL

AA_pep_T7 <- norm_T7.df[c("Peptide", "Pep_k1", "Pep_k2")]
AA_pep_T7 <- apply(AA_pep_T7, 1, function(x) {
  AA_peptide_cutoff(x, light_time)}) %>% bind_rows
row.names(AA_pep_T7) <- NULL

norm_T6.df <- merge(norm_T6.df, AA_pep_T6, by="Peptide")
norm_T7.df <- merge(norm_T7.df, AA_pep_T7, by="Peptide")
```

``` r
# Control channels: M0 approaches kt/kd asymptotically
simple_normal_model <- function(t, m0_n, kd, kt) {
  (kt / kd) + (m0_n - (kt / kd)) * exp(-kd * t) }

# Heavy channels: original pool decays, synthesis contributes light peptide
simple_heavy_model <- function(t, m0_h, kd, kt, k1, k2) {
  m0_h * exp(-kd * t) + #Simplified model - no baseline
    (kt * k1 / (kd - k2)) * (exp(-k2 * t) - exp(-kd * t)) }

# Combined kt/kd fits for simplified model
simple_pred_model <- function(t, heavy_ch, ctrl_ch, m0_heavy, m0_ctrl,
                              kd, kt, pep_k1, pep_k2) {
  
  light_time <- tail(t, ctrl_ch)
  heavy_time <- t[1:heavy_ch]

  ctrl_pred <- simple_normal_model(light_time, m0_ctrl, kd, kt)
  heavy_pred <- simple_heavy_model(heavy_time, m0_heavy, kd, kt,
                                   pep_k1, pep_k2)

  pred <- c(heavy_pred, ctrl_pred)
  
  return(pred) }
```

``` r
# Control channel: two-pool model
# Dynamic pool approaches kt/kd; stable pool sits at phi*m0_n
twopool_normal_model <- function(t, m0_n, kd, kt, phi) {
  dynamic <- (kt / kd) + ((1 - phi) * m0_n - (kt / kd)) * exp(-kd * t)
  stable <- phi * m0_n
  
  dynamic + stable }

# Heavy channel: two-pool model
# Dynamic pool decays + synthesis adds light protein; stable pool sits at phi*m0_h
twopool_heavy_model <- function(t, m0_h, kd, kt, k1, k2, phi) {
  dynamic_decay <- (1 - phi) * m0_h * exp(-kd * t)
  synthesis <- (kt * k1 / (kd - k2)) * (exp(-k2 * t) - exp(-kd * t))
  stable <- phi * m0_h
  
  dynamic_decay + synthesis + stable }

# Combined kt/kd fits for simplified model
twopool_pred_model <- function(t, heavy_ch, ctrl_ch, m0_heavy, m0_ctrl,
                               kd, kt, pep_k1, pep_k2, phi) {
  
  light_time <- tail(t, ctrl_ch)
  heavy_time <- t[1:heavy_ch]

  ctrl_pred <- twopool_normal_model(light_time, m0_ctrl, kd, kt, phi)
  heavy_pred <- twopool_heavy_model(heavy_time, m0_heavy, kd, kt,
                                    pep_k1, pep_k2, phi)

  pred <- c(heavy_pred, ctrl_pred)
  
  return(pred) }
```

``` r
peptide_fit <- function(row, l_time, h_time) {

  #--- Model input ----------------------------  
  
  pep_k1 <- as.numeric(row["Pep_k1"])
  pep_k2 <- as.numeric(row["Pep_k2"])

  heavy_data <- as.numeric(row[raw_heavy_labels])

  light_cutoff <- as.numeric(row["N_timepoints"])
  light_trim <- l_time[1:light_cutoff]
  sub_light <- as.numeric(row[raw_light_labels])[1:(light_cutoff + 1)]

  #Input for nlsLM
  all_data <- c(heavy_data, sub_light)
  all_time <- c(c(0, h_time), c(0, light_trim))

  light_n <- light_cutoff + 1
  heavy_n <- length(h_time) + 1

  #--- Starting Guess ----------------------------  

  #Starting channel parameters
  m0_heavy_start <- mean(heavy_data[1:2])
  m0_ctrl_start <- mean(sub_light[1:2])

  #kD starting parameter 
  se_heavy <- pmax(heavy_data[1:(3+1)], 1e-6) #Prevent 0 after log
  se_h_time <- c(0, h_time[1:3])
  
  se_fit <- lm(log(se_heavy) ~ se_h_time) #Approximating linearly

  kd_start <- -coef(se_fit)[2]
  kd_start <- max(kd_start, 1e-6)

  #kT starting parameter
  se_light <- c(m0_ctrl_start, pmax(tail(sub_light, 1), 1e-6))
  se_l_time <- c(light_trim[1], tail(light_trim, 1))
  
  initial_slope <- (se_light[2] - se_light[1]) / (se_l_time[2] - se_l_time[1])
  kt_start <- initial_slope + kd_start * se_light[1]
  kt_start <- max(kt_start, 1e-6)  # keep non-negative

  #Baseline starting parameter
  phi_start <- max(tail(heavy_data, 1), 0.01)  
  
  #--- Bounds ----------------------------  
  
  m0_heavy_lower <- max(m0_heavy_start * 0.5, 0.1)
  m0_ctrl_lower <- max(m0_ctrl_start * 0.5, 0.1)
  C_upper <- mean(tail(heavy_data, 3)) * 1.5
  kd_upper <- log(2) / 5    # Shortest half-life ~ 5 min
  kt_upper <- kd_upper * m0_heavy_start * 2  
  
  #--- No Baseline Model ----------------------------    

  simple_fit <- tryCatch( #No baseline model
    nlsLM(all_data ~ simple_pred_model(all_time, heavy_n, light_n,
                                       m0_heavy, m0_ctrl, kd, kt,
                                       pep_k1, pep_k2),
          start = list(m0_heavy = m0_heavy_start,
                       m0_ctrl = m0_ctrl_start,
                       kd = kd_start,
                       kt = kt_start),
          lower = c(m0_heavy = m0_heavy_lower,
                    m0_ctrl = m0_ctrl_lower,
                    kd = 1e-6,
                    kt = 0),
          upper = c(m0_heavy = heavy_data[1] * 1.2,
                    m0_ctrl = sub_light[1] * 1.2,
                    kd = kd_upper,
                    kt = kt_upper),
          control = nls.lm.control(maxiter = 1000)),
    error = function(e) NULL )

  #Extracting fit features from object  
  if (is.null(simple_fit)) {
    simple_pred_values <- c(m0_heavy = NA, m0_ctrl = NA, kd = NA, kt = NA)
    simple_predict <- rep(NA, length(all_data))
    simple_RSS <- NA }
  else {
    simple_pred_values <- coef(simple_fit) 
    simple_predict <- predict(simple_fit)
    simple_RSS <- deviance(simple_fit) }

  #--- No Baseline Model ----------------------------    

  max_phi <- 0.25
  if (phi_start >= max_phi) { phi_start <- max_phi - 0.01 }
  
  baseline_fit <- tryCatch( #No baseline model
    nlsLM(all_data ~ twopool_pred_model(all_time, heavy_n, light_n,
                                        m0_heavy, m0_ctrl, kd, kt,
                                        pep_k1, pep_k2, phi),
          start = list(m0_heavy = m0_heavy_start,
                       m0_ctrl = m0_ctrl_start,
                       kd = kd_start,
                       kt = kt_start,
                       phi = phi_start),
          lower = c(m0_heavy = m0_heavy_lower,
                    m0_ctrl = m0_ctrl_lower,
                    kd = 1e-6,
                    kt = 0,
                    phi = 0),
          upper = c(m0_heavy = heavy_data[1] * 1.2,
                    m0_ctrl = sub_light[1] * 1.2,
                    kd = kd_upper,
                    kt = kt_upper,
                    phi = max_phi),
          control = nls.lm.control(maxiter = 1000)),
    error = function(e) NULL )    

  #Extracting fit features from object
  if (is.null(baseline_fit)) {
    baseline_pred_values <- c(m0_heavy = NA, m0_ctrl = NA, kd = NA, kt = NA,
                              phi = NA)
    baseline_predict <- rep(NA, length(all_data))
    baseline_RSS <- NA }
  else {
    baseline_pred_values <- coef(baseline_fit)
    baseline_predict <- predict(baseline_fit)
    baseline_RSS <- deviance(baseline_fit) }
  
  #--- Output ----------------------------

  if ( is.na(simple_RSS) ) { deltaBIC <- NA }
  else if ( is.na(baseline_RSS) ) { deltaBIC <- NA }
  else {
    total_n <- length(all_data)
    simple_BIC <- total_n * log(simple_RSS/total_n) + 4 * log(total_n)
    baseline_BIC <- total_n * log(baseline_RSS/total_n) + 5 * log(total_n)

    deltaBIC <- ifelse(baseline_BIC < simple_BIC, abs(baseline_BIC - simple_BIC),
                       -1 * abs(baseline_BIC - simple_BIC)) }

  # plot(c(0, l_time), as.numeric(row[raw_light_labels]), ylim=c(0, 3))
  # lines(c(0, light_trim), simple_predict[15:length(all_data)],
  #       col="red", lty = "dashed")
  # lines(c(0, light_trim), predict(baseline_fit)[15:length(all_data)],
  #       col="blue", lty = "dashed")
  # 
  # plot(c(0, h_time), heavy_data, ylim=c(0, 3))
  # lines(c(0, h_time), simple_predict[1:14], col="red")
  # lines(c(0, h_time), predict(baseline_fit)[1:14], col="blue")
  
  col_names <- c("Peptide", "Pep_k1", "Pep_k2", "N_Light", "deltaBIC",
                 paste0("Simple_", names(simple_pred_values)),
                 paste0("Baseline_", names(baseline_pred_values)))
  out_df <- data.frame(matrix(nrow=0, ncol = length(col_names)))
  colnames(out_df) <- col_names

  out_df[1,"Peptide"] <- row["Peptide"]
  out_df[1, 2:14] <- c(pep_k1, pep_k2, light_cutoff, deltaBIC,
                      simple_pred_values, baseline_pred_values)
  
  return(out_df) }


#---------------------------------------------------------------------------------------

fit_k_rates <- function(protein_id, pep_df, l_time, h_time) {

  pep_df <- pep_df[pep_df$Protein_ID == protein_id,]
  pepfit <- apply(pep_df, 1, function(x) { peptide_fit(x, l_time, h_time) })
  pepfit <- do.call(rbind, pepfit)
  pepfit <- cbind(data.frame(Protein_ID = rep(protein_id, nrow(pepfit))),
                  pepfit) 

  return(pepfit) }

# fit_k_rates("P14785", norm_T6.df, light_time, Fly_T6.time.min)
# fit_k_rates("P14785", norm_T7.df, light_time, Fly_T6.time.min)

#Setup cluster
cl <- makeCluster(num_cores)

#Export necessary variables/functions to the cluster nodes
clusterExport(cl,
              varlist=c("raw_light_labels", "raw_heavy_labels",
                        "light_time", "Fly_T6.time.min",
                        "norm_T6.df", "norm_T7.df", "fit_k_rates",
                        "peptide_fit", "pep_decay_rate",
                        "simple_normal_model",
                        "simple_heavy_model", "simple_pred_model",
                        "twopool_pred_model",
                        "twopool_normal_model", "twopool_heavy_model"))
invisible(clusterEvalQ(cl, library(minpack.lm)))

T6_pep_fits <- parLapply(cl, unique(norm_T6.df$Protein_ID), function(x){
  fit_k_rates(x, norm_T6.df, light_time, Fly_T6.time.min) }) %>%
  bind_rows

T7_pep_fits <- parLapply(cl, unique(norm_T7.df$Protein_ID), function(x){
  fit_k_rates(x, norm_T7.df, light_time, Fly_T6.time.min) }) %>%
  bind_rows

stopCluster(cl) #Stop Cluster

head(T6_pep_fits)
```

    ##       Protein_ID           Peptide Pep_k1      Pep_k2 N_Light  deltaBIC
    ## 1         Q9I7R7 AAAAAAAATNQVQLAGK      1 0.015478680       3 -3.685094
    ## 41242     Q9I7R7        TPNESDVDPK      1 0.013313740       3 -7.033032
    ## 2         Q9W2E9       AAAAAAGLAGK      1 0.009145265       3 10.895908
    ## 9200      Q9W2E9    ELDFYQNIPQDILK      1 0.017339778       3  6.829749
    ## 9573      Q9W2E9          ELNEGGFK      1 0.011942307       3 16.775457
    ## 24213     Q9W2E9       LLLLNQSTVIK      1 0.013586221       3  4.184563
    ##       Simple_m0_heavy Simple_m0_ctrl   Simple_kd   Simple_kt Baseline_m0_heavy
    ## 1            1.469821      0.5756223 0.002148805 0.004894240          1.493237
    ## 41242        1.402184      0.4975947 0.001955186 0.004749963          1.419537
    ## 2            2.142342      0.5043628 0.005756680 0.004751768          2.196892
    ## 9200         2.349692      0.4452723 0.006584075 0.006440900          2.386514
    ## 9573         1.694507      0.5573716 0.002858347 0.002451342          1.799429
    ## 24213        1.465574      0.4374989 0.001884116 0.002683137          1.465574
    ##       Baseline_m0_ctrl Baseline_kd Baseline_kt Baseline_phi
    ## 1            0.6452374 0.003827689 0.006047281   0.25000000
    ## 41242        0.6157427 0.002205476 0.004578262   0.06981496
    ## 2            0.4988638 0.008102381 0.006224474   0.10475177
    ## 9200         0.4432271 0.008259270 0.007570390   0.07938473
    ## 9573         0.5269289 0.007084930 0.005045704   0.25000000
    ## 24213        0.4374989 0.004634846 0.004858480   0.25000000

``` r
head(T7_pep_fits)
```

    ##       Protein_ID         Peptide Pep_k1      Pep_k2 N_Light    deltaBIC
    ## 1         Q9W2E9     AAAAAAGLAGK      1 0.009145265       3 -10.4960788
    ## 22794     Q9W2E9    LIDFAHTAFVPR      1 0.012553276       3 -27.3070332
    ## 2         Q9VM36      AAAAAALAAK      1 0.009505175       3  -5.1385539
    ## 10973     Q9VM36  EVSPAPVFNIFTPR      1 0.014888111       3   0.2742286
    ## 26419     Q9VM36 LVVDSDGGGSPLAQK      1 0.014545172       3  -3.6686305
    ## 37599     Q9VM36        TAATFLNR      1 0.009289932       3  -8.2453304
    ##       Simple_m0_heavy Simple_m0_ctrl   Simple_kd   Simple_kt Baseline_m0_heavy
    ## 1            2.016617      0.6154977 0.005051726 0.004688959          1.978626
    ## 22794        2.276933      0.6473962 0.006016795 0.005577055          2.129363
    ## 2            1.410669      0.9639478 0.002159530 0.004441200          1.414909
    ## 10973        1.561992      1.2231511 0.002133573 0.002530200          1.617399
    ## 26419        1.633016      1.1319403 0.002644045 0.004720329          1.618898
    ## 37599        1.481828      1.1209672 0.002441609 0.004557408          1.463767
    ##       Baseline_m0_ctrl Baseline_kd Baseline_kt Baseline_phi
    ## 1            0.5994839 0.009940814 0.007659983         0.25
    ## 22794        0.6265199 0.010232755 0.008134466         0.25
    ## 2            1.0735891 0.003815949 0.005691816         0.25
    ## 10973        1.1935039 0.004625811 0.004645933         0.25
    ## 26419        1.1368141 0.002361357 0.004294229         0.00
    ## 37599        1.2731003 0.004025814 0.005511714         0.25

``` r
max_fit_kd <- log(2)/(tail(Fly_T6.time.min, 1)*3)
max_fit_hl <- log(2)/max_fit_kd/60
cat("Max Half-life (hrs):\t", max_fit_hl) 
```

    ## Max Half-life (hrs):  62.85

``` r
T6_pep_fits <- T6_pep_fits[!is.na(T6_pep_fits$Simple_kd),]
norm_T6.df <- norm_T6.df[norm_T6.df$Peptide %in% T6_pep_fits$Peptide,]

T7_pep_fits <- T7_pep_fits[!is.na(T7_pep_fits$Simple_kd),]
norm_T7.df <- norm_T7.df[norm_T7.df$Peptide %in% T7_pep_fits$Peptide,]
```

``` r
pep_fit_comp <- T6_pep_fits
pep_fit_comp <- pep_fit_comp[!is.na(pep_fit_comp$deltaBIC),]

pep_fit_comp["Model"] <- rep("Unselected", nrow(pep_fit_comp))
pep_fit_comp[pep_fit_comp$deltaBIC < 0, "Model"] <- "Simple"
pep_fit_comp[pep_fit_comp$deltaBIC > 0, "Model"] <- "Two-pool"
# pep_fit_comp[pep_fit_comp$deltaBIC > 10, "deltaBIC"] <- 10

table(pep_fit_comp$Model)
```

    ## 
    ##   Simple Two-pool 
    ##    44867     3722

``` r
p1 <- ggplot() +
  geom_point(data=pep_fit_comp[pep_fit_comp$Model == "Simple",],
             aes(x=log10(Simple_kd), y=log10(Baseline_kd)),
             size=2, alpha=0.5) +
  geom_abline(slope=1, intercept=0, linetype="dashed", linewidth=1,
              color="red") +
  geom_vline(xintercept = log10(max_fit_kd),
             linewidth=1, color="hotpink") +
  geom_hline(yintercept = log10(max_fit_kd),
             linewidth=1, color="hotpink") +  
  theme_bw() +
  theme(aspect.ratio=1) +
  labs(x=expression("Simple log"[10]*"(min"^-1*")"),
       y=expression("Baseline log"[10]*"(min"^-1*")"),
       title="Prefers Simple Model") +
  coord_cartesian(xlim=c(-6, -1), ylim=c(-6, -1))

p2 <- ggplot() +
  geom_point(data=pep_fit_comp[pep_fit_comp$Model == "Two-pool",],
             aes(x=log10(Simple_kd), y=log10(Baseline_kd),
                 color=deltaBIC),
             size=2, alpha=0.5) +
  geom_abline(slope=1, intercept=0, linetype="dashed", linewidth=1,
              color="red") +
  geom_vline(xintercept = log10(max_fit_kd),
             linewidth=1, color="hotpink") +
  geom_hline(yintercept = log10(max_fit_kd),
             linewidth=1, color="hotpink") +  
  theme_bw() +
  theme(aspect.ratio=1) +
  labs(x=expression("Simple log"[10]*"(min"^-1*")"),
       y=expression("Baseline log"[10]*"(min"^-1*")"),
       title="Prefers Baseline Model") +
  coord_cartesian(xlim=c(-6, -1), ylim=c(-6, -1))


# #Uncomment for new image!
# tiff("graph1.tiff", units="in",
#      width=9, height=4.5, res=300)

p1 | p2
```

![](figures/Fly_Gast_SimpModel/Exploratory%20analysis%20not%20included%20in%20manuscript-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
protein_rollup <- function(protein_id, all_reps){

  sub_data <-  all_reps[all_reps$Protein_ID==protein_id,]

  final_calc_lst <- lapply(unique(sub_data$Exp), function(x){
    
    sub_exp <- sub_data[sub_data$Exp == x,]
    sub_exp["Fraction"] <- sub_exp$sum_sn / sum(sub_exp$sum_sn)
  
    sub_exp[light_labels] <- sub_exp[light_labels] * sub_exp$Fraction
    light_avg <- apply(sub_exp[light_labels], 2, sum)
    light_avg <- light_avg / mean(light_avg)
  
    sub_exp[heavy_labels] <- sub_exp[heavy_labels] * sub_exp$Fraction
    heavy_avg <- apply(sub_exp[heavy_labels], 2, sum)
    heavy_avg <- heavy_avg / mean(heavy_avg)

    return(c(light_avg, heavy_avg)) })
  
  merged_matrix <- do.call(rbind, final_calc_lst)
  merged_avg <- apply(merged_matrix, 2, mean)

  merged_avg[light_labels] <-
    merged_avg[light_labels] / mean(merged_avg[light_labels])
  merged_avg[heavy_labels] <-
    merged_avg[heavy_labels] / mean(merged_avg[heavy_labels])

  protein_avg <- cbind(data.frame(Protein_ID = protein_id),
                       data.frame(t(merged_avg)))

  return(protein_avg) }


protein_med_merge <- rbind(cbind(norm_T6.df,
                                 data.frame(Exp=rep("T6", nrow(norm_T6.df)))),
                           cbind(norm_T7.df,
                                 data.frame(Exp=rep("T7", nrow(norm_T7.df)))))

protein_med_merge["T0_H"] <- (protein_med_merge$T0A_H + protein_med_merge$T0B_H)/2
protein_med_merge["T0_L"] <- (protein_med_merge$T0A_L + protein_med_merge$T0B_L)/2
protein_med_merge <- protein_med_merge[c("Protein_ID", "Peptide", heavy_labels,
                                         light_labels, "sum_sn", "Exp")]

protein_med_merge[light_labels] <-
    sweep(protein_med_merge[light_labels], 1,
          rowMeans(protein_med_merge[light_labels]), `/`)

protein_med_merge[heavy_labels] <-
    sweep(protein_med_merge[heavy_labels], 1,
          rowMeans(protein_med_merge[heavy_labels]), `/`)

protein_med_merge[light_labels] <-
    sweep(protein_med_merge[light_labels], 1,
          rowMeans(protein_med_merge[light_labels]), `/`)

protein_med_merge[heavy_labels] <-
    sweep(protein_med_merge[heavy_labels], 1,
          rowMeans(protein_med_merge[heavy_labels]), `/`)

#Setup cluster
cl <- makeCluster(num_cores)

#Export necessary variables/functions to the cluster nodes
clusterExport(cl,
              varlist=c("protein_med_merge", "protein_rollup",
                        "light_labels", "heavy_labels"))

protein_med_merge <- parLapply(cl, unique(protein_med_merge$Protein_ID),
       function(x){ protein_rollup(x, protein_med_merge) }) %>% bind_rows

stopCluster(cl) #Stop Cluster

# protein_med_merge <- merge(Dmel.genes, protein_med_merge,
#                            by="Protein_ID", all.y=TRUE)
# protein_med_merge[is.na(protein_med_merge$Gene_Symbol),] <- "N/A"

head(protein_med_merge)
```

    ##   Protein_ID      T0_L       N13       N14       N15       N16     T0_H
    ## 1     Q9I7R7 0.5313668 0.8710916 1.2123800 1.2782350 1.1069266 1.275749
    ## 2     Q9W2E9 0.5438511 0.5535339 0.8671569 1.2662170 1.7692411 2.344955
    ## 3     Q9VM36 0.8648921 1.4350585 1.2537853 0.8199373 0.6263269 1.249177
    ## 4     Q05825 1.0639623 1.0323401 1.0233273 0.8879462 0.9924240 1.090494
    ## 5     A1Z7P3 1.0703280 1.1402480 1.1130509 0.8983857 0.7779873 1.134712
    ## 6     Q9VQF7 0.4160977 0.6779988 1.2934199 1.3296812 1.2828024 1.009165
    ##         O1       O2       O3       O4        O5        O6        O7        O8
    ## 1 1.578229 1.625240 1.424104 1.365290 1.2052000 1.0161671 1.0197758 0.8734550
    ## 2 2.200645 2.034438 1.580315 1.260032 0.8567635 0.6714147 0.6192902 0.4871164
    ## 3 1.350394 1.647791 1.513066 1.346929 1.2428441 1.1681994 1.0762592 0.8131878
    ## 4 1.096226 1.112457 1.143692 1.135890 1.1165950 1.1283786 0.9811863 0.9274079
    ## 5 1.139799 1.262229 1.204924 1.242674 1.1942737 1.1158906 1.0763806 0.9550416
    ## 6 1.219946 1.196597 1.248324 1.356657 1.3328551 1.2254250 1.0562172 0.9800123
    ##          O9       O10       O11       O12
    ## 1 0.6005151 0.4568678 0.2940004 0.2654074
    ## 2 0.3049155 0.2556798 0.1941214 0.1903133
    ## 3 0.5976583 0.4383453 0.3076702 0.2484774
    ## 4 0.8898213 0.8421844 0.7990755 0.7365918
    ## 5 0.8373459 0.7497583 0.5855531 0.5014185
    ## 6 0.7565328 0.6919483 0.4712030 0.4551165

# BIC Calculation

``` r
merge_pep_fits <- function(protein_id, exp1, exp2) {

  exp1_sub <- exp1[exp1$Protein_ID == protein_id,]
  exp2_sub <- exp2[exp2$Protein_ID == protein_id,]

  pep_data <- rbind(exp1_sub[c("Peptide", "Pep_k1", "Pep_k2")],
                    exp2_sub[c("Peptide", "Pep_k1", "Pep_k2")])
  pep_data <- unique(pep_data)

  if (nrow(exp1_sub) > 0 & nrow(exp2_sub) > 0) {

    fitted_values_1 <- exp1_sub[c("Simple_m0_heavy", "Simple_kd",
                                  "Simple_kt")]
    fitted_values_1 <- apply(fitted_values_1, 2, median)
    
    fitted_values_2 <- exp2_sub[c("Simple_m0_heavy",
                                  "Simple_kd", "Simple_kt")]
    fitted_values_2 <- apply(fitted_values_2, 2, median)
    
    estimate_desc <- "interval"
    exp1_kD <- fitted_values_1["Simple_kd"]
    exp2_kD <- fitted_values_2["Simple_kd"]    
    merged_values <- (fitted_values_1 + fitted_values_2)/2 }

  else if (nrow(exp1_sub) > 0 & nrow(exp2_sub) == 0) {

    estimate_desc <- "point"            
    merged_values <- exp1_sub[c("Simple_m0_heavy", "Simple_kd", "Simple_kt")]
    merged_values <- apply(merged_values, 2, median) 
    
    exp1_kD <- merged_values["Simple_kd"]
    exp2_kD <- NA }  

  else if (nrow(exp1_sub) == 0 & nrow(exp2_sub) > 0) {

    estimate_desc <- "point"    
    merged_values <- exp2_sub[c("Simple_m0_heavy", "Simple_kd", "Simple_kt")]
    merged_values <- apply(merged_values, 2, median)
    
    exp1_kD <- NA
    exp2_kD <- merged_values["Simple_kd"] }
  
  else { cat("Logic Error"); invokeRestart("abort") }

  names(merged_values) <- c("m0_heavy", "kd", "kt")
  
  pep_predict <- apply(pep_data[2:3], 1, function(sub_pep){

    predict_fit <- simple_heavy_model(Fly_T6.time.min, merged_values["m0_heavy"],
                                      merged_values["kd"], merged_values["kt"],
                                      sub_pep["Pep_k1"], sub_pep["Pep_k2"])
    
    return(data.frame(t(predict_fit))) })
  
  pep_predict <- do.call(rbind, pep_predict)
  pep_predict <- apply(pep_predict, 2, median)

  out_df <- data.frame(t(pep_predict))
  colnames(out_df) <- paste0("Fit_", heavy_labels)
  
  out_df <- cbind(data.frame(Protein_ID = protein_id,
                             kD = merged_values["kd"],
                             HL_Hrs = log(2)/merged_values["kd"]/60),                             
                             Estimate = estimate_desc,
                             T6_kD = exp1_kD, T7_kD = exp2_kD,
                  out_df)
  
  return(out_df) }

Deg_fits.df <- lapply(protein_med_merge$Protein_ID,
                      function(x) { merge_pep_fits(x, T6_pep_fits, T7_pep_fits) }) %>%
  bind_rows() %>% arrange(desc(kD))

head(Deg_fits.df)
```

    ##        Protein_ID         kD    HL_Hrs Estimate      T6_kD      T7_kD Fit_T0_H
    ## kd...1     Q9VFD5 0.01539471 0.7504170    point 0.01539471         NA 2.929405
    ## kd...2     Q9VMA3 0.01486216 0.7773063 interval 0.01538270 0.01434163 3.041118
    ## kd...3     Q8IPM1 0.01408818 0.8200105    point         NA 0.01408818 2.803578
    ## kd...4     Q9VED4 0.01264941 0.9132802    point         NA 0.01264941 3.400825
    ## kd...5     P09085 0.01217044 0.9492226    point         NA 0.01217044 3.114322
    ## kd...6     Q9VNG1 0.01191980 0.9691817 interval 0.01182441 0.01201520 2.728303
    ##          Fit_O1   Fit_O2    Fit_O3    Fit_O4    Fit_O5     Fit_O6     Fit_O7
    ## kd...1 1.983756 1.294051 0.8675530 0.5532929 0.2738329 0.12934374 0.05689680
    ## kd...2 1.993018 1.250774 0.8066763 0.4911285 0.2247297 0.09690259 0.03825395
    ## kd...3 2.093416 1.500632 1.0877691 0.7499252 0.4118304 0.21308312 0.10159221
    ## kd...4 2.413816 1.655138 1.1609221 0.7776667 0.4141334 0.21053224 0.09984398
    ## kd...5 2.263455 1.579918 1.1197868 0.7545693 0.4016675 0.20211359 0.09399044
    ## kd...6 2.160221 1.646179 1.2608573 0.9231171 0.5553132 0.31577473 0.16656747
    ##             Fit_O8       Fit_O9      Fit_O10      Fit_O11      Fit_O12
    ## kd...1 0.013137672 0.0013160458 2.784652e-04 1.875479e-05 5.466843e-07
    ## kd...2 0.007134431 0.0004874736 7.740818e-05 3.021407e-06 4.003099e-08
    ## kd...3 0.026019385 0.0028346834 6.090654e-04 3.970124e-05 1.019169e-06
    ## kd...4 0.026027105 0.0030590669 7.075558e-04 5.392853e-05 1.760474e-06
    ## kd...5 0.023318065 0.0024833684 5.336222e-04 3.559632e-05 9.729899e-07
    ## kd...6 0.050885804 0.0072600417 1.864837e-03 1.648609e-04 6.250107e-06

``` r
Deg_fits.df["deltaBIC"] <- apply(Deg_fits.df, 1, function(row){

  protein_id <- row["Protein_ID"]
  dyn_fit <- as.numeric(row[grepl("Fit", names(row))])
  
  med_data <- protein_med_merge[protein_med_merge$Protein_ID == protein_id,]
  med_data <- as.numeric(med_data[heavy_labels])
  
  total_n <- length(med_data)
  
  flat_line <- rep(1, total_n)

  model_RSS <- sum((med_data - dyn_fit)^2)
  flat_RSS <- sum((med_data - flat_line)^2)
    
  model_BIC <- total_n * log(model_RSS/total_n) + 4 * log(total_n)
  flat_BIC <- total_n * log(flat_RSS/total_n) + 0 * log(total_n)
  
  deltaBIC <- ifelse(model_BIC < flat_BIC, abs(model_BIC - flat_BIC),
                     -1 * abs(model_BIC - flat_BIC))

  return(deltaBIC) })

head(Deg_fits.df)
```

    ##        Protein_ID         kD    HL_Hrs Estimate      T6_kD      T7_kD Fit_T0_H
    ## kd...1     Q9VFD5 0.01539471 0.7504170    point 0.01539471         NA 2.929405
    ## kd...2     Q9VMA3 0.01486216 0.7773063 interval 0.01538270 0.01434163 3.041118
    ## kd...3     Q8IPM1 0.01408818 0.8200105    point         NA 0.01408818 2.803578
    ## kd...4     Q9VED4 0.01264941 0.9132802    point         NA 0.01264941 3.400825
    ## kd...5     P09085 0.01217044 0.9492226    point         NA 0.01217044 3.114322
    ## kd...6     Q9VNG1 0.01191980 0.9691817 interval 0.01182441 0.01201520 2.728303
    ##          Fit_O1   Fit_O2    Fit_O3    Fit_O4    Fit_O5     Fit_O6     Fit_O7
    ## kd...1 1.983756 1.294051 0.8675530 0.5532929 0.2738329 0.12934374 0.05689680
    ## kd...2 1.993018 1.250774 0.8066763 0.4911285 0.2247297 0.09690259 0.03825395
    ## kd...3 2.093416 1.500632 1.0877691 0.7499252 0.4118304 0.21308312 0.10159221
    ## kd...4 2.413816 1.655138 1.1609221 0.7776667 0.4141334 0.21053224 0.09984398
    ## kd...5 2.263455 1.579918 1.1197868 0.7545693 0.4016675 0.20211359 0.09399044
    ## kd...6 2.160221 1.646179 1.2608573 0.9231171 0.5553132 0.31577473 0.16656747
    ##             Fit_O8       Fit_O9      Fit_O10      Fit_O11      Fit_O12
    ## kd...1 0.013137672 0.0013160458 2.784652e-04 1.875479e-05 5.466843e-07
    ## kd...2 0.007134431 0.0004874736 7.740818e-05 3.021407e-06 4.003099e-08
    ## kd...3 0.026019385 0.0028346834 6.090654e-04 3.970124e-05 1.019169e-06
    ## kd...4 0.026027105 0.0030590669 7.075558e-04 5.392853e-05 1.760474e-06
    ## kd...5 0.023318065 0.0024833684 5.336222e-04 3.559632e-05 9.729899e-07
    ## kd...6 0.050885804 0.0072600417 1.864837e-03 1.648609e-04 6.250107e-06
    ##         deltaBIC
    ## kd...1  7.832369
    ## kd...2 10.631890
    ## kd...3 11.704248
    ## kd...4 23.447389
    ## kd...5 19.866379
    ## kd...6 18.057844

``` r
true_shuffle <- function(x) {
  
  n <- length(x)
  original <- x
  shuffled <- sample(x)

  #Keep re-rolling until NO elements match their starting position
  # -> This prevents correlation for fastest degraders
  #Also, stop from last timepoint ending up with the first one
  # -> This prevents anti-correlation for fastest degraders
  while(any(shuffled == original) || shuffled[2] == original[1] ||
        shuffled[n] == original[1] || shuffled[n-1] == original[1]) {
    shuffled <- sample(x) }

  T0_spot <- which(shuffled == "T0_H")
  shuffle_expanded <- append(shuffled[-T0_spot], sample(c("T0A_H", "T0B_H")),
                             after = T0_spot - 1)
  
  return(shuffle_expanded) }
```

``` r
final_fly_gast <- merge(Dmel.genes, Deg_fits.df, by="Protein_ID", all.y=TRUE)
final_fly_gast <- final_fly_gast[c(colnames(final_fly_gast)[1:7],
                                   "deltaBIC")]
final_fly_gast <- merge(final_fly_gast, protein_med_merge, by="Protein_ID")
final_fly_gast <- final_fly_gast %>% arrange(desc(kD))

head(final_fly_gast)
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
    ##          O9       O10       O11       O12
    ## 1 0.5719438 0.5362716 0.3889168 0.5384916
    ## 2 0.4306688 0.4159514 0.4550532 0.4928102
    ## 3 0.3652615 0.4694481 0.2848861 0.3645625
    ## 4 0.1570708 0.1219197 0.1053457 0.1170699
    ## 5 0.2565727 0.2858838 0.3429187 0.2663379
    ## 6 0.3697731 0.3911873 0.2936905 0.2900057

``` r
dyn_model_pep_fit <- function(row, l_time, h_time, heavy_header) {

  #--- Model input ----------------------------  
  
  pep_k1 <- as.numeric(row["Pep_k1"])
  pep_k2 <- as.numeric(row["Pep_k2"])

  heavy_data <- as.numeric(row[heavy_header])

  light_cutoff <- as.numeric(row["N_timepoints"])
  light_trim <- l_time[1:light_cutoff]
  sub_light <- as.numeric(row[raw_light_labels])[1:(light_cutoff + 1)]

  #Input for nlsLM
  all_data <- c(heavy_data, sub_light)
  all_time <- c(c(0, h_time), c(0, light_trim))

  light_n <- light_cutoff + 1
  heavy_n <- length(h_time) + 1

  #--- Starting Guess ----------------------------  

  #Starting channel parameters
  m0_heavy_start <- mean(heavy_data[1:2])
  m0_ctrl_start <- mean(sub_light[1:2])

  #kD starting parameter 
  se_heavy <- pmax(heavy_data[1:(3+1)], 1e-6) #Prevent 0 after log
  se_h_time <- c(0, h_time[1:3])
  
  se_fit <- lm(log(se_heavy) ~ se_h_time) #Approximating linearly

  kd_start <- -coef(se_fit)[2]
  kd_start <- max(kd_start, 1e-6)

  #kT starting parameter
  se_light <- c(m0_ctrl_start, pmax(tail(sub_light, 1), 1e-6))
  se_l_time <- c(light_trim[1], tail(light_trim, 1))
  
  initial_slope <- (se_light[2] - se_light[1]) / (se_l_time[2] - se_l_time[1])
  kt_start <- initial_slope + kd_start * se_light[1]
  kt_start <- max(kt_start, 1e-6)  # keep non-negative

  #--- Bounds ----------------------------  
  
  m0_heavy_lower <- max(m0_heavy_start * 0.5, 0.1)
  m0_ctrl_lower <- max(m0_ctrl_start * 0.5, 0.1)
  kd_upper <- log(2) / 5    # Shortest half-life ~ 5 min
  kt_upper <- kd_upper * m0_heavy_start * 2  
  
  #--- No Baseline Model ----------------------------    

  simple_fit <- tryCatch( #No baseline model
    nlsLM(all_data ~ simple_pred_model(all_time, heavy_n, light_n,
                                       m0_heavy, m0_ctrl, kd, kt,
                                       pep_k1, pep_k2),
          start = list(m0_heavy = m0_heavy_start,
                       m0_ctrl = m0_ctrl_start,
                       kd = kd_start,
                       kt = kt_start),
          lower = c(m0_heavy = m0_heavy_lower,
                    m0_ctrl = m0_ctrl_lower,
                    kd = 1e-6,
                    kt = 0),
          upper = c(m0_heavy = heavy_data[1] * 1.2,
                    m0_ctrl = sub_light[1] * 1.2,
                    kd = kd_upper,
                    kt = kt_upper),
          control = nls.lm.control(maxiter = 1000)),
    error = function(e) NULL )

  #Extracting fit features from object  
  if (is.null(simple_fit)) {
    simple_pred_values <- c(m0_heavy = NA, m0_ctrl = NA,
                            kd = NA, kt = NA) }
  else {
    simple_pred_values <- coef(simple_fit) }

  return(data.frame(t(simple_pred_values))) }
```

``` r
#Setup cluster
cl <- makeCluster(num_cores)
clusterSetRNGStream(cl, 123)   # reproducible permutation null

#Export necessary variables/functions to the cluster nodes
clusterExport(cl,
              varlist=c("light_labels", "heavy_labels",
                        "raw_light_labels", "raw_heavy_labels",
                        "Fly_T6.time.min", "light_time",
                        "norm_T6.df", "norm_T7.df", "protein_med_merge",
                        "dyn_model_pep_fit", "true_shuffle",
                        "simple_normal_model",
                        "simple_heavy_model", "simple_pred_model"))
invisible(clusterEvalQ(cl, library(minpack.lm)))

NULL_BIC_values <- parApply(cl, final_fly_gast, 1, function(row) {

  protein_id <- row["Protein_ID"]
  total_n <- length(heavy_labels)
  flat_predict <- rep(1, total_n)

  median_data <- as.numeric(row[heavy_labels])
  T6_check <- is.na(row["T6_kD"])
  T7_check <- is.na(row["T7_kD"])

  real_deltaBIC <- as.numeric(row["deltaBIC"])

  if (T6_check) {
    T7_peps <- norm_T7.df[norm_T7.df$Protein_ID == protein_id,]
    pep_data <- T7_peps[c("Peptide", "Pep_k1", "Pep_k2")] }
  else if (T7_check) {
    T6_peps <- norm_T6.df[norm_T6.df$Protein_ID == protein_id,]
    pep_data <- T6_peps[c("Peptide", "Pep_k1", "Pep_k2")] }
  else {
    T6_peps <- norm_T6.df[norm_T6.df$Protein_ID == protein_id,]
    T7_peps <- norm_T7.df[norm_T7.df$Protein_ID == protein_id,]

    pep_data <- rbind(T6_peps[c("Peptide", "Pep_k1", "Pep_k2")],
                      T7_peps[c("Peptide", "Pep_k1", "Pep_k2")])
    pep_data <- unique(pep_data) }

  NULL_BIC_deltas <- sapply(1:10, function(dummy) {

    #-- Fit ------------------------------------------

    new_header <- true_shuffle(heavy_labels)
    # new_header <- raw_heavy_labels
    prot_header <- new_header[new_header != "T0B_H"]
    prot_header[prot_header == "T0A_H"] <- "T0_H"

    scrambled_median <- as.numeric(row[prot_header])

    if (T6_check) {

      T7_fits <- apply(T7_peps, 1, function(x) {
        dyn_model_pep_fit(x, light_time, Fly_T6.time.min, new_header) } )
      T7_fits <- do.call(rbind, T7_fits)
      merged_values <- apply(T7_fits, 2, median, na.rm = TRUE) }

    else if (T7_check) {

      T6_fits <- apply(T6_peps, 1, function(x) {
        dyn_model_pep_fit(x, light_time, Fly_T6.time.min, new_header) } )
      T6_fits <- do.call(rbind, T6_fits)
      merged_values <- apply(T6_fits, 2, median, na.rm = TRUE) }

    else {

      T6_fits <- apply(T6_peps, 1, function(x) {
        dyn_model_pep_fit(x, light_time, Fly_T6.time.min, new_header) } )
      T6_fits <- do.call(rbind, T6_fits)
      T6_fits <- apply(T6_fits, 2, median, na.rm = TRUE)

      T7_fits <- apply(T7_peps, 1, function(x) {
        dyn_model_pep_fit(x, light_time, Fly_T6.time.min, new_header) } )
      T7_fits <- do.call(rbind, T7_fits)
      T7_fits <- apply(T7_fits, 2, median, na.rm = TRUE)

      merged_values <- (T6_fits + T7_fits)/2 }

    #-- Prediction ------------------------------------------
    if (any(is.na(merged_values))) {
      deltaBIC <- (-100)
      return(deltaBIC) }
    else { } # ... normal BIC computation ...

    pep_predict <- apply(pep_data[2:3], 1, function(sub_pep){

      predict_fit <- simple_heavy_model(Fly_T6.time.min, merged_values["m0_heavy"],
                                        merged_values["kd"], merged_values["kt"],
                                        sub_pep["Pep_k1"], sub_pep["Pep_k2"])

      return(data.frame(t(predict_fit))) })

    pep_predict <- do.call(rbind, pep_predict)
    pep_predict <- apply(pep_predict, 2, median)

    #-- BIC Calculation ------------------------------------------

    model_RSS <- sum((scrambled_median - pep_predict)^2)
    flat_RSS <- sum((scrambled_median - flat_predict)^2)

    model_BIC <- total_n * log(model_RSS/total_n) + 4 * log(total_n)
    flat_BIC <- total_n * log(flat_RSS/total_n) + 0 * log(total_n)

    deltaBIC <- ifelse(model_BIC < flat_BIC, abs(model_BIC - flat_BIC),
                       -1 * abs(model_BIC - flat_BIC))

    return(deltaBIC) })

  return(NULL_BIC_deltas) })

NULL_BIC_values <- as.vector(NULL_BIC_values)
NULL_BIC_values[1:20]
```

    ##  [1] -11.516785 -11.150376 -11.770823 -10.611328 -11.362982 -10.449153
    ##  [7] -12.713022 -11.495564 -12.652571  -8.837024 -10.269309 -10.290269
    ## [13] -12.306836 -10.243774 -13.213287 -12.221833 -12.469832 -15.065883
    ## [19] -13.511406 -13.296494

``` r
stopCluster(cl) #Stop Cluster
```

``` r
real_scores <- final_fly_gast$deltaBIC

current_cutoff <- 20.0
current_FDR <- 0

#Walk down by 0.1 until FDR hits 5%
while(current_FDR < 0.05) {

  #How many real proteins pass this cutoff?
  real_hits <- sum(real_scores > current_cutoff)

  # How many noise proteins pass this cutoff?
  # - Divide by dataset multiplier size
  expected_false_positives <- sum(NULL_BIC_values > current_cutoff) / 10

  #Calculate FDR
  current_FDR <- expected_false_positives / real_hits

  #Step the cutoff down
  current_cutoff <- current_cutoff - 0.1 }

#The loop breaks when FDR hits 5%.
final_BIC_threshold <- current_cutoff + 0.1

paste("The 5% FDR Threshold is Delta BIC >", final_BIC_threshold)
```

    ## [1] "The 5% FDR Threshold is Delta BIC > -7.8"

``` r
plot_data <- data.frame(
  BIC_Score = c(NULL_BIC_values, real_scores),
  Distribution = c(rep("Null (Noise)", length(NULL_BIC_values)),
                   rep("Real Data", length(real_scores))))
plot_data["Distribution"] <- factor(plot_data$Distribution,
                                    levels=c("Null (Noise)", "Real Data"))
plot_data$Weight <- ifelse(plot_data$Distribution == "Real Data", 1, 1/10)

hist_deltaBIC <- ggplot(plot_data, aes(x = BIC_Score, fill = Distribution, weight = Weight)) +
  geom_histogram(position = "identity", binwidth = 1, alpha=0.5) +
  geom_vline(xintercept = final_BIC_threshold, linetype="dashed", color="black",
             linewidth=1) +
  labs(x = expression(Delta * " BIC Score"),
       y = "Expected Frequency") +
  theme_bw() +
  theme(axis.text=element_text(size=21,colour="black"),
        plot.title = element_text(size=24, hjust=0.5),
        axis.title = element_text(size=22),
        panel.border = element_rect(linewidth=2),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
        legend.position = "none",
        plot.margin = margin(t = 10, r = 25, b = 10,
                             l = 10, unit = "pt")) +
  scale_fill_manual(values = c("Null (Noise)" = "black", "Real Data" = "#56B4E9"))

# #Uncomment for new image!
# tiff("fly_gast_deltaBIC_hist_10variations.tiff", units="in",
#      width=6, height=4.5, res=300)

hist_deltaBIC
```

![](figures/Fly_Gast_SimpModel/BIC%20cutoff%20plot-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
final_BIC_threshold <- 2

final_fly_gast["mClass"] <- rep("Flat", nrow(final_fly_gast))
final_fly_gast[(final_fly_gast$deltaBIC > final_BIC_threshold),
               "mClass"] <- "Deg"
final_fly_gast[final_fly_gast$kD < max_fit_kd, "HL_Hrs"] <- max_fit_hl

table(final_fly_gast$mClass)
```

    ## 
    ##  Deg Flat 
    ## 5358  509

``` r
head(final_fly_gast)
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
ss_bp.all <- final_fly_gast

ss_bp.all <- ss_bp.all[!is.na(ss_bp.all$T6_kD),]
ss_bp.all <- ss_bp.all[!is.na(ss_bp.all$T7_kD),]

ss_bp.all["Avg_kD"] <- (ss_bp.all$T6_kD + ss_bp.all$T7_kD)/2

ss_bp.all["Log10_T6"] <- log10(ss_bp.all$T6_kD)
ss_bp.all["Log10_T7"] <- log10(ss_bp.all$T7_kD)
ss_bp.all[ss_bp.all$Log10_T6<(-5), "Log10_T6"] <- (-6)
ss_bp.all[ss_bp.all$Log10_T7<(-5), "Log10_T7"] <- (-6)

cor(ss_bp.all$Log10_T6, ss_bp.all$Log10_T7)^2
```

    ## [1] 0.5254217

``` r
cor(ss_bp.all[ss_bp.all$mClass=="Deg",]$Log10_T6,
    ss_bp.all[ss_bp.all$mClass=="Deg",]$Log10_T7)^2
```

    ## [1] 0.917243

``` r
p1 <- ggplot() +
  geom_point(data=ss_bp.all[ss_bp.all$mClass == "Deg",],
             aes(x=Log10_T6, y=Log10_T7),
             shape=1, size=3, stroke=1, alpha=0.5) +  

  geom_hline(yintercept=log10(max_fit_kd), size=1) +
  geom_vline(xintercept=log10(max_fit_kd), size=1) +

  geom_abline(slope=1, intercept=0, color="black", linetype="dotted",
              linewidth=1) +
  labs(x=expression("R1 log"[10]*"(min"^-1*")"),
       y=expression("R2 log"[10]*"(min"^-1*")")) +
  theme_bw() +
  theme(axis.text=element_text(size=24,colour="black"),
        plot.title = element_text(size=24, hjust=0.5),
        axis.title = element_text(size=22),
        panel.border = element_rect(linewidth=2),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
        aspect.ratio = 1,
        legend.position = "none",
        plot.margin = margin(t = 10, r = 25, b = 10,
                             l = 10, unit = "pt")) +
  coord_cartesian(xlim = c(-6,0),
                  ylim = c(-6,0)) +
  scale_color_manual(values=c("#5E4FA2", "#767676"))
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once per session.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
# #Uncomment for new image!
# tiff("fly_with_badproteins.tiff", units="in",
#      width=5, height=5, res=300)

p1
```

![](figures/Fly_Gast_SimpModel/Degradation%20rate%20biplot%20for%20degrading%20proteins-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
kmeans_plot <- function(data, k = 5, seed = 123,
                        order=seq(1,k,1)) {
  
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
  
  all.kresult[light_labels] <- apply(all.kresult[light_labels], 1,
                                       function(row){
    row <- row / mean(row)
    row <- row / row[1]
    return(data.frame(t(row))) }) %>% bind_rows

  all.kresult[heavy_labels] <- apply(all.kresult[heavy_labels], 1,
                                       function(row){
    row <- row / mean(row)
    row <- row / row[1]
    return(data.frame(t(row))) }) %>% bind_rows

  
  all.kresult <- all.kresult %>%
    pivot_longer(cols = -c(Cluster, Size),
                 names_to = "Time", values_to = "Value")
  
  whole_t_vec <- c(Fly_T6_hpf_light_time, Fly_T6_hpf_time)

  all.kresult["Time_Hrs"] <- rep(whole_t_vec,
                                 nrow(all.kresult) /
                                   length(whole_t_vec))


  plots <- lapply(seq_len(k), function(i){
    sub_cl <- all.kresult[all.kresult$Cluster == i,]

    light.df <- sub_cl[sub_cl$Time %in% light_labels,]
    heavy.df <- sub_cl[sub_cl$Time %in% heavy_labels,]
    heavy.df <- heavy.df[!heavy.df$Time == "O1",]
    
    ggplot() +
      geom_line(data = heavy.df, aes(x = Time_Hrs, y = Value, group = 1),
                size = 3, color = "#36A893") +
      geom_line(data = light.df, aes(x = Time_Hrs, y = Value, group = 1),
                size = 3, color = "#767676") +
      geom_hline(yintercept = 1, linetype = "dashed", size = 1) +
      theme_bw() +
      labs(x = "", y = "Rel. Abundance",
           title = paste0("n=", unique(sub_cl$Size))) +
      theme(axis.text = element_text(size = 18, colour = "black"),
            plot.title = element_text(size = 18, hjust = 0.5),
            axis.title = element_text(size = 22),
            axis.title.y = element_blank(),
            aspect.ratio = 1,
            panel.border = element_rect(linewidth = 1.25),
            legend.position = "none",
            panel.grid.major = element_line(color = "grey70", linewidth = 0.25),
            panel.grid.minor = element_line(color = "grey70", linewidth = 0.25)) +
      coord_cartesian(ylim = c(0, 3.5))
  })
  
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

export_kmeans <- kmeans_plot(protein_med_merge, order = c(1,5,4,3,2))
```

![](figures/Fly_Gast_SimpModel/kmeans%20plot%20of%20protein%20averages-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
fly_kd_hist <- final_fly_gast[c("Protein_ID", "kD", "mClass")]
fly_kd_hist["Plot_kD"] <- fly_kd_hist$kD
fly_kd_hist[fly_kd_hist$kD > 0.001, "Plot_kD"] <- 0.001

fly_kd_hist$mClass <- factor(fly_kd_hist$mClass,
                             levels = c("Flat", "Deg"))

#Define exact bin edges to guarantee the bars and the outline align perfectly
bw <- 0.000025
brks <- seq(1e-6, 0.001 + bw, by = bw)

#Calculate half-life x-intercepts (kd = ln(2) / minutes)
kd_1day <- log(2) / (24 * 60)
kd_12hr <- log(2) / (12 * 60)
kd_6hr  <- log(2) / (6 * 60)

fly_kd_pdf <- ggplot(fly_kd_hist, aes(x = Plot_kD, fill = mClass,
                                      color = mClass)) +
  # after_stat(density) normalizes the area to 1
  # position = "identity" prevents stacking so they overlap naturally
  geom_histogram(aes(y = after_stat(count / sum(count)), alpha = mClass), 
                 breaks = brks, position = "identity", 
                 color=NA) +
  
  # ONLY draws an outline for "Deg"
  geom_step(aes(y = after_stat(count / sum(count))),
            stat = "bin", breaks = brks, direction = "mid",
            linewidth = 1, show.legend = FALSE) +
  
  geom_vline(xintercept = c(kd_1day, kd_12hr, kd_6hr),
             linetype = "dashed", color = "black", 
             linewidth = 1) +

  geom_vline(xintercept = 0.0005514297,
             linetype = "dashed", color = "red3",
             linewidth = 1) +

  scale_fill_manual(values = c("Flat" = "grey55", "Deg" = "#5E4FA2")) + 
  scale_color_manual(values = c("Flat" = NA, "Deg" = "#5E4FA2")) +
  scale_alpha_manual(values = c("Flat" = 0.9, "Deg" = 0.7)) +
  
  # LABELS & THEME: Matching your reference image
  labs(x = expression("k"["d"]~"(min"^{-1}*") | Proteins = 5,867"), y = "Density") +
  theme_bw() +
  theme(axis.text=element_text(size=22,colour="black"),
        plot.title = element_text(size=24, hjust=0.5),
        axis.title = element_text(size=23),
        panel.border = element_rect(linewidth=2),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
        legend.position = "none",
        plot.margin = margin(t = 10, r = 25, b = 10,
                             l = 10, unit = "pt")) +
  coord_cartesian(xlim=c(1e-6, 0.001),
                  ylim=c(0, 0.5)) +
  scale_x_continuous(breaks=c(0, 0.00025, 0.0005, 0.00075,  0.001),
                     labels=c("0", "0.00025", "0.0005",
                              "0.00075",  "0.001"))

# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=9, height=3.5, res=300)

fly_kd_pdf
```

    ## Warning: Removed 40 rows containing missing values or values outside the scale range
    ## (`geom_step()`).

![](figures/Fly_Gast_SimpModel/Fly%20kd%20distribution-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
plot_data <- final_fly_gast[final_fly_gast$mClass == "Deg",]
plot_data["Label"] <- rep("Base", nrow(plot_data))
plot_data[plot_data$kD < max_fit_kd, "Label"] <- "Collapsed"
plot_data[plot_data$kD < max_fit_kd, "HL_Hrs"] <- max_fit_hl
plot_data["Label"] <- factor(plot_data$Label, levels = c("Collapsed", "Base"))

median_HL <- median(plot_data$HL_Hrs, na.rm = TRUE)

hl_hist <- ggplot(data = plot_data, aes(x = HL_Hrs, fill=Label, alpha=Label)) +
  geom_histogram(bins = 30, color = "black", linewidth = 0.5) +
  geom_vline(xintercept = median_HL, color = "red3", linetype = "dashed", linewidth = 1.2) +
  labs(x = "Half-life (hours)",
       y = "No. of proteins",
       title = " ") +
  theme_bw() +
  theme(axis.text=element_text(size=18,colour="black"),
        plot.title = element_text(size=22, hjust=0.5),
        axis.title = element_text(size=20),
        panel.border = element_rect(linewidth=1.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
        legend.position = "none",
        plot.margin = margin(t = 10, r = 25, b = 10,
                             l = 10, unit = "pt")) +
  coord_cartesian(ylim=c(0,500)) +
  scale_fill_manual(values=c("Collapsed" = "#9ECAE1", "Base" = "#5E4FA2")) +
  scale_alpha_manual(values=c("Collapsed" = 1, "Base" = 0.7))

# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=8, height=3.5, res=300)

hl_hist
```

![](figures/Fly_Gast_SimpModel/Fly%20HL%20distribution%20for%20degrading%20proteins-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
plot_example_protein <- function(protein_id) {

  sub_df <- final_fly_gast[final_fly_gast$Protein_ID == protein_id,]
  print(sub_df$HL_Hrs)
  
  sub_df <- sub_df[c("Protein_ID", light_labels, heavy_labels)]

  sub_df[light_labels] <- sub_df[light_labels] / sub_df[[light_labels[1]]]
  sub_df[heavy_labels] <- sub_df[heavy_labels] / sub_df[[heavy_labels[1]]]
    
  sub_df <- sub_df %>%
    pivot_longer(cols=-Protein_ID)
  sub_df["Label"] <- c(rep("Control", length(light_labels)),
                       rep("Label", length(heavy_labels)))
  sub_df["Label"] <- factor(sub_df$Label, levels=c("Label", "Control"))
  
  sub_df["Time"] <- c(Fly_T6_hpf_light_time, Fly_T6_hpf_time) 

  ggplot() +
    geom_line(data=sub_df[sub_df$Label=="Label",],
              aes(x=Time, y=value, group=1),
              size=3, color="#36A893") +
    geom_line(data=sub_df[sub_df$Label=="Control",],
              aes(x=Time, y=value, group=1),
              size=3, color="#767676") +
    geom_hline(yintercept=1, linewidth=1, linetype="dashed") +
    theme_bw() +
    labs(x="Hours post-fertilization", y="Rel. Abundance") +
    theme(axis.text=element_text(size=21,colour="black"),
          plot.title = element_text(size=24, hjust=0.5),
          axis.title = element_text(size=22),
          panel.border = element_rect(linewidth=2),
          panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
          panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
          aspect.ratio = 1,
          legend.position = "none",
          plot.margin = margin(t = 10, r = 25, b = 10,
                               l = 10, unit = "pt")) +
    coord_cartesian(ylim=c(0,2.5)) }

# # #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=4, height=4, res=300)

plot_example_protein("P18824")
```

    ## [1] 4.590503

![](figures/Fly_Gast_SimpModel/Plot%20of%20example%20proteins-1.png)<!-- -->

``` r
plot_example_protein("P07487")
```

    ## [1] 47.97983

![](figures/Fly_Gast_SimpModel/Plot%20of%20example%20proteins-2.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
write.csv(final_fly_gast,
          "Files/Fits/Dyn-Model/Fly_O18_T6-7_FinalFits_ALL.csv",
          row.names = FALSE)
```

``` r
peptide_count <- data.frame(table(rbind(norm_T6.df[1:2],
                                        norm_T7.df[1:2])$Protein_ID))
colnames(peptide_count) <- c("Protein_ID", "Number_Pep_Fits")

head(peptide_count)
```

    ##   Protein_ID Number_Pep_Fits
    ## 1 A0A021WW64               1
    ## 2 A0A023GRW3               8
    ## 3 A0A0A1EI90               1
    ## 4 A0A0B4JCU3               6
    ## 5 A0A0B4JCZ0              13
    ## 6 A0A0B4JCZ1               3

``` r
supp_df <- merge(protein_annotation, final_fly_gast,
                 by="Protein_ID", all.y=TRUE) %>%
  dplyr::select(Protein_ID, Gene_Symbol, Name, Description, everything())
supp_df <- merge(supp_df, peptide_count, by="Protein_ID") %>%
  arrange(desc(kD))

#--------------------------------------------

supp_df["Light_FC"] <- apply(supp_df, 1, function(row){
  
  light_fc <- as.numeric(row["N16"]) / as.numeric(row["T0_L"])
  
  return(light_fc) })

inc_min_fc <- 1.439371 
dec_min_fc <- 0.6218052 

#--------------------------------------------

supp_df["Control_Type"] <- rep("Unchanging", nrow(supp_df))
supp_df[supp_df$Light_FC > inc_min_fc, "Control_Type"] <- "Increasing"
supp_df[supp_df$Light_FC < dec_min_fc, "Control_Type"] <- "Decreasing"

supp_df["Is_Degrading"] <- rep("Unknown", nrow(supp_df))
supp_df[supp_df$mClass == "Deg", "Is_Degrading"] <- "Yes"
supp_df[supp_df$mClass == "Flat", "Is_Degrading"] <- "No"

#--------------------------------------------

supp_df["Combined_Class"] <- paste0(supp_df$Control_Type, "_",
                                    supp_df$Is_Degrading)
supp_df[supp_df$Combined_Class == "Decreasing_No",
        "Control_Type"] <- "Unchanging"
supp_df["Combined_Class"] <- paste0(supp_df$Control_Type, "_",
                                    supp_df$Is_Degrading)

supp_df[supp_df$Combined_Class == "Unchanging_No",
        "Combined_Class"] <- "No measurable change"
supp_df[supp_df$Combined_Class == "Unchanging_Yes",
        "Combined_Class"] <- "Protein synthesis & deg."
supp_df[supp_df$Combined_Class == "Increasing_No",
        "Combined_Class"] <- "Primarily protein synthesis"
supp_df[supp_df$Combined_Class == "Increasing_Yes",
        "Combined_Class"] <- "Protein increase with deg."
supp_df[supp_df$Combined_Class == "Decreasing_Yes",
        "Combined_Class"] <- "Primarily protein degradation"

#--------------------------------------------

supp_df[c("mClass", "Light_FC")] <- NULL

colnames(supp_df) <-
  c("Protein ID", "FlyBase Symbol", "Name", "Description",
    "Fitted kD", "Reported Half-life (Hrs)", "Estimate Confidence",
    "Replicate 1 kD", "Replicate 2 kD", "BIC against NULL",
    paste0("Control: ", round(Fly_T6_hpf_light_time, 1), " hpf"),
    paste0("18O: ", round(Fly_T6_hpf_time, 1), " hpf"),
    "Peptides Quantified", "Control Type", "Is Degrading?",
    "Grouped Class")

supp_df[(supp_df$`Is Degrading?`=="No") &
          (supp_df$`Reported Half-life (Hrs)` < max_fit_hl),
        "Reported Half-life (Hrs)"] <- NA
supp_df["Grouped Class"] <- NULL
supp_df <- rbind(supp_df[!is.na(supp_df$`Reported Half-life (Hrs)`),],
                 supp_df[is.na(supp_df$`Reported Half-life (Hrs)`),])

dir.create("Files/Supp_Tables", showWarnings = FALSE, recursive = TRUE)
write.csv(supp_df, "Files/Supp_Tables/Fly_Gastrulation_Supplementary_Table.csv",
          row.names = FALSE)
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
    ## [1] minpack.lm_1.2-4 patchwork_1.3.2  purrr_1.2.2      ggplot2_4.0.3   
    ## [5] stringr_1.6.0    tidyr_1.3.2      dplyr_1.2.0     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.6       compiler_4.5.3     tidyselect_1.2.1   scales_1.4.0      
    ##  [5] yaml_2.3.12        fastmap_1.2.0      R6_2.6.1           labeling_0.4.3    
    ##  [9] generics_0.1.4     knitr_1.51         tibble_3.3.1       pillar_1.11.1     
    ## [13] RColorBrewer_1.1-3 rlang_1.1.7        stringi_1.8.7      xfun_0.58         
    ## [17] S7_0.2.1           otel_0.2.0         cli_3.6.5          withr_3.0.3       
    ## [21] magrittr_2.0.4     digest_0.6.39      grid_4.5.3         rstudioapi_0.19.0 
    ## [25] lifecycle_1.0.5    vctrs_0.7.1        evaluate_1.0.5     glue_1.8.0        
    ## [29] farver_2.1.2       rmarkdown_2.31     tools_4.5.3        pkgconfig_2.0.3   
    ## [33] htmltools_0.5.9
