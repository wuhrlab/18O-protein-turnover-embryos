Frog 2C kD Fits - Simplified
================
Edward Cruz
2026-06-05

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
library(minpack.lm)
library(parallel)

knitr::opts_chunk$set(fig.path = "figures/Frog_2C_SimpModel/")
num_cores <- detectCores() - 1  # Save one core for system stability
```

``` r
human_map <- read.csv("Files/Reference/Xen10_to_Human_Name_Assign.csv")
human_map["Protein_ID"] <- sapply(human_map$Protein_ID, function(x){
  split_id <- strsplit(x, "\\|")[[1]]
  reformat <- paste0(split_id[3], "|", split_id[2])

  return(reformat) })

protein_annotation <- human_map[c("Protein_ID", "Description")]
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

light.channels <- c("T0", "N1", "N3", "N4",
                    "N7", "N9",  "N11", "N12")
heavy.channels <- c("T0", "O1", "O2", "O3", "O4", "O5",
                    "O7", "O9", "O10", "O11", "O12")

Frog_T5_all.min <- c(161, 221, 281, 401, 521, 641,
                     881, 1605, 1961, 2918, 4405) - 161
Frog_T5_light.min <- Frog_T5_all.min[c(1:2, 4:5, 7:8, 10:11)]

Frog_T6_all.min <- c(150, 210, 270, 390, 510, 635,
                     930, 1591, 1950, 2910, 4394) - 150
Frog_T6_light.min <- Frog_T6_all.min[c(1:2, 4:5, 7:8, 10:11)]

avg_light_time <- (Frog_T5_light.min + Frog_T6_light.min)/2
avg_heavy_time <- (Frog_T5_all.min + Frog_T6_all.min)/2

#Files to normalize
XLA_O18_T5  <- read.csv("Data/XLA_O18/XLA_O18-T5_Filtered_Data.csv")
XLA_O18_T6  <- read.csv("Data/XLA_O18/XLA_O18-T6_Filtered_Data.csv")

norm_set <- read.csv("Data/XLA_Norm/XLA-O18_YolkNormSet_T8-NYS-Decay.csv")$Protein_ID
```

``` r
#Removing things we can't say anything about -> human keratin contam. + frog spinout
XLA_O18_T5 <- XLA_O18_T5[!XLA_O18_T5$Protein_ID %in% c(frog_yolk_set, keratins),]
XLA_O18_T6 <- XLA_O18_T6[!XLA_O18_T6$Protein_ID %in% c(frog_yolk_set, keratins),]

XLA_O18_T5 <- XLA_O18_T5[!rowSums(XLA_O18_T5[heavy.channels]) == 0,]
XLA_O18_T6 <- XLA_O18_T6[!rowSums(XLA_O18_T6[heavy.channels]) == 0,]

XLA_O18_T5["StartPool"] <- apply(XLA_O18_T5[heavy.channels], 1,
      function(row){
        starting_pool <- sum(row[c("T0", "O1")]) / sum(row)
        return(starting_pool) })

XLA_O18_T5["HeavyWeight"] <- apply(XLA_O18_T5[c(light.channels, heavy.channels)], 1,
      function(row){

        light_sn <- row[1:length(light.channels)]
        heavy_sn <- row[(length(light.channels)+1):length(row)]
        
        return(mean(light_sn)/(mean(heavy_sn) + 1e-6)) })

XLA_O18_T6["StartPool"] <- apply(XLA_O18_T6[heavy.channels], 1,
      function(row){
        starting_pool <- sum(row[c("T0", "O1")]) / sum(row)
        return(starting_pool) })

XLA_O18_T6["HeavyWeight"] <- apply(XLA_O18_T6[c(light.channels, heavy.channels)], 1,
      function(row){

        light_sn <- row[1:length(light.channels)]
        heavy_sn <- row[(length(light.channels)+1):length(row)]
        
        return(mean(light_sn)/(mean(heavy_sn) + 1e-6)) })


#This filter sets a floor to quantify based on signal
# -> Requires (T0 + O1)/11 > 0.15
# -> Requires heavy signal to not outweigh light signal
XLA_O18_T5 <- XLA_O18_T5[(XLA_O18_T5$StartPool > 0.15) &
                           (XLA_O18_T5$HeavyWeight > 0.8),]
XLA_O18_T6 <- XLA_O18_T6[(XLA_O18_T6$StartPool > 0.15) &
                           (XLA_O18_T6$HeavyWeight > 0.8),]
```

``` r
split_experiments <- function(df){
  base_df <- df[c("Protein_ID", "Peptide", "Pep_k1", "Pep_k2")]

  heavy_df <- df[heavy.channels]
  heavy_df <- heavy_df/rowSums(heavy_df)
  heavy_df <- heavy_df/rowMeans(heavy_df)
  colnames(heavy_df)[1] <- c("T0_H")

  light_df <- df[light.channels]
  light_df <- light_df/rowSums(light_df)
  light_df <- light_df/rowMeans(light_df)  
  colnames(light_df)[1] <- c("T0_L")  
  
  split_df <- cbind(base_df, heavy_df, light_df,
                    df[c("sum_sn", "StartPool", "HeavyWeight")])

  return(split_df) }
  
XLA_O18_T5 <- split_experiments(XLA_O18_T5)
XLA_O18_T6 <- split_experiments(XLA_O18_T6)

light.channels <- c("T0_L", "N1", "N3", "N4",
                    "N7", "N9",  "N11", "N12")
heavy.channels <- c("T0_H", "O1", "O2", "O3", "O4", "O5",
                    "O7", "O9", "O10", "O11", "O12")
```

``` r
norm_Frog_split <- function(df, L_error, H_error) {
  
  L_error <- apply(df[df$Protein_ID %in% norm_set, light.channels],
                   2, median)
  L_error <- L_error / mean(L_error)
  H_error <- apply(df[df$Protein_ID %in% norm_set, heavy.channels],
                   2, median)
  H_error <- H_error / mean(H_error)

  print(c(L_error, H_error))
  cat("\n")
  #Divide each original ratio by the pipet error
  # ... Divide each column by the row sum
  H_ratio.df <- as.matrix(df[heavy.channels])
  H_corr.ratios <- sweep(H_ratio.df, 2, H_error, `/`)
  H_corr.ratios <- sweep(H_corr.ratios, 1, rowMeans(H_corr.ratios), `/`)

  #Divide each original ratio by the pipet error
  # ... Divide each column by the row sum
  L_ratio.df <- as.matrix(df[light.channels])
  L_corr.ratios <- sweep(L_ratio.df, 2, L_error, `/`)
  L_corr.ratios <- sweep(L_corr.ratios, 1, rowMeans(L_corr.ratios), `/`)

  #Normalized dataframe then shifting to 1 -> (1/18 Ratio) = 1
  norm.df <- df
  norm.df[heavy.channels] <- H_corr.ratios
  norm.df[light.channels] <- L_corr.ratios

  return(as.data.frame(norm.df)) }

norm_T5.df <- norm_Frog_split(XLA_O18_T5, L_error_T5, H_error_T5)
```

    ##      T0_L        N1        N3        N4        N7        N9       N11       N12 
    ## 1.2141990 1.0721352 1.1271691 1.0485090 1.0844280 0.9295866 0.7808328 0.7431402 
    ##      T0_H        O1        O2        O3        O4        O5        O7        O9 
    ## 1.1375238 1.1562778 1.0045707 1.0664445 1.0708354 1.0851515 0.9560543 0.9962837 
    ##       O10       O11       O12 
    ## 0.8601071 0.9227582 0.7439929

``` r
norm_T5.df <- norm_T5.df[norm_T5.df$T0_H > 0.8,]
norm_T5.df <- norm_T5.df[norm_T5.df$T0_L > 0.15,]

norm_T6.df <- norm_Frog_split(XLA_O18_T6, L_error_T6, H_error_T6)
```

    ##      T0_L        N1        N3        N4        N7        N9       N11       N12 
    ## 1.1356730 1.0747581 1.0427457 1.1003610 1.0840537 1.0710511 0.8201571 0.6712004 
    ##      T0_H        O1        O2        O3        O4        O5        O7        O9 
    ## 1.1070418 1.1673554 1.0091527 1.0031861 1.0489369 1.1628553 0.9664357 1.0410583 
    ##       O10       O11       O12 
    ## 0.9134799 0.8778372 0.7026608

``` r
norm_T6.df <- norm_T6.df[norm_T6.df$T0_H > 0.8,]
norm_T6.df <- norm_T6.df[norm_T6.df$T0_L > 0.15,]
```

``` r
# Determing where synthesis isn't relevant
pep_decay_rate <- function(fraction, k1, k2) {
  (-1 * log(fraction / k1)) / k2 }

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

AA_pep_T5 <- norm_T5.df[c("Peptide", "Pep_k1", "Pep_k2")]
AA_pep_T5 <- apply(AA_pep_T5, 1, function(x) {
  AA_peptide_cutoff(x, Frog_T5_light.min)}) %>% bind_rows
row.names(AA_pep_T5) <- NULL

AA_pep_T6 <- norm_T6.df[c("Peptide", "Pep_k1", "Pep_k2")]
AA_pep_T6 <- apply(AA_pep_T6, 1, function(x) {
  AA_peptide_cutoff(x, Frog_T6_light.min)}) %>% bind_rows
row.names(AA_pep_T6) <- NULL

norm_T5.df <- merge(norm_T5.df, AA_pep_T5, by="Peptide")
norm_T6.df <- merge(norm_T6.df, AA_pep_T6, by="Peptide")
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

  heavy_data <- as.numeric(row[heavy.channels])

  light_cutoff <- as.numeric(row["N_timepoints"])
  light_trim <- l_time[1:light_cutoff]
  sub_light <- as.numeric(row[light.channels])[1:light_cutoff]
  
  #Input for nlsLM
  all_data <- c(heavy_data, sub_light)  
  all_time <- c(h_time, light_trim)

  heavy_n <- length(all_time) - length(light_trim)
  light_n <- length(light_trim)
  
  #--- Starting Guess ----------------------------  

  #Starting channel parameters
  m0_heavy_start <- heavy_data[1]
  m0_ctrl_start <- sub_light[1]

  #kD starting parameter 
  se_heavy <- pmax(heavy_data[1:3], 1e-6) #Prevent 0 after log
  se_h_time <- h_time[1:3]
  
  se_fit <- lm(log(se_heavy) ~ se_h_time) #Approximating linearly

  kd_start <- -coef(se_fit)[2]
  kd_start <- max(kd_start, 1e-6)

  #kT starting parameter
  se_light <- c(sub_light[1], pmax(tail(sub_light, 1), 1e-6))
  se_l_time <- c(light_trim[1], tail(light_trim, 1))

  initial_slope <- (se_light[2] - se_light[1]) / (se_l_time[2] - se_l_time[1])
  kt_start <- initial_slope + kd_start * se_light[1]
  kt_start <- max(kt_start, 1e-6)  # keep non-negative
  
  #Baseline starting parameter
  phi_start <- max(mean(tail(heavy_data, 3)), 0.01)

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

  # print(deltaBIC)
  # plot(h_time, all_data[1:11])
  # lines(h_time, simple_predict[1:11], col="red")
  # lines(h_time, baseline_predict[1:11], col="blue")
  # 
  # plot(l_time, as.numeric(row[light.channels]))
  # lines(light_trim, simple_predict[12:length(all_data)], col="red", lty="dashed")
  # lines(light_trim, baseline_predict[12:length(all_data)], col="blue", lty="dashed")

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

# fit_k_rates("XBmRNA10569|XBXL10_1g5638", norm_T5.df, Frog_T5_light.min, Frog_T5_all.min)
# fit_k_rates("XBmRNA10569|XBXL10_1g5638", norm_T6.df, Frog_T6_light.min, Frog_T6_all.min)

#Setup cluster
cl <- makeCluster(num_cores)

#Export necessary variables/functions to the cluster nodes
clusterExport(cl,
              varlist=c("light.channels", "heavy.channels",
                        "Frog_T5_light.min", "Frog_T5_all.min",
                        "Frog_T6_light.min", "Frog_T6_all.min",
                        "norm_T5.df", "norm_T6.df", "fit_k_rates",
                        "peptide_fit", "pep_decay_rate", "simple_normal_model",
                        "simple_heavy_model", "simple_pred_model",
                        "twopool_pred_model",
                        "twopool_normal_model", "twopool_heavy_model"))
invisible(clusterEvalQ(cl, library(minpack.lm)))

T5_pep_fits <- parLapply(cl, unique(norm_T5.df$Protein_ID), function(x){
  fit_k_rates(x, norm_T5.df, Frog_T5_light.min, Frog_T5_all.min) }) %>%
  bind_rows

T6_pep_fits <- parLapply(cl, unique(norm_T6.df$Protein_ID), function(x){
  fit_k_rates(x, norm_T6.df, Frog_T6_light.min, Frog_T6_all.min) }) %>%
  bind_rows

stopCluster(cl) #Stop Cluster

head(T5_pep_fits)
```

    ##                      Protein_ID         Peptide Pep_k1     Pep_k2 N_Light
    ## 1     XBmRNA11342|XBXL10_1g6076   AAAAAAGASGSQK      1 0.06037187       2
    ## 3677  XBmRNA11342|XBXL10_1g6076         AVLNTLK      1 0.01421940       3
    ## 8351  XBmRNA11342|XBXL10_1g6076       EEDMMELLK      1 0.04945150       2
    ## 11396 XBmRNA11342|XBXL10_1g6076 ESTVQGAFSDSAELK      1 0.09112961       2
    ## 12127 XBmRNA11342|XBXL10_1g6076       EVYIAMGQK      1 0.02857930       3
    ## 13562 XBmRNA11342|XBXL10_1g6076   FNEQIEPIYALLR      1 0.04326544       2
    ##        deltaBIC Simple_m0_heavy Simple_m0_ctrl    Simple_kd    Simple_kt
    ## 1     -2.553404       1.0137200      0.3636715 3.389666e-05 0.0015819719
    ## 3677  -2.639058       0.9784899      0.4145035 1.531724e-06 0.0003882838
    ## 8351  -2.565712       0.9664322      0.2596390 2.459337e-05 0.0024252763
    ## 11396 -2.557181       1.1069536      0.4329860 9.947948e-05 0.0001281779
    ## 12127 -2.639030       0.9782866      0.3616827 1.628690e-06 0.0004563400
    ## 13562 -2.315780       1.0627316      0.4414723 5.447875e-05 0.0000000000
    ##       Baseline_m0_heavy Baseline_m0_ctrl  Baseline_kd  Baseline_kt Baseline_phi
    ## 1             1.0141785        0.3635604 4.614119e-05 0.0015865195         0.25
    ## 3677          0.9784790        0.4145077 1.516455e-06 0.0003882339         0.00
    ## 8351          0.9664322        0.2599947 2.406264e-05 0.0024113407         0.00
    ## 11396         1.1096032        0.4325153 1.427217e-04 0.0001472977         0.25
    ## 12127         0.9782866        0.3616737 2.172210e-06 0.0004564568         0.25
    ## 13562         1.0637406        0.4407579 7.641822e-05 0.0000000000         0.25

``` r
head(T6_pep_fits)
```

    ##                      Protein_ID                 Peptide Pep_k1     Pep_k2
    ## 1     XBmRNA11342|XBXL10_1g6076           AAAAAAGASGSQK      1 0.06037187
    ## 3792  XBmRNA11342|XBXL10_1g6076                 AVLNTLK      1 0.01421940
    ## 8818  XBmRNA11342|XBXL10_1g6076               EEDMMELLK      1 0.04945150
    ## 11904 XBmRNA11342|XBXL10_1g6076         ESTVQGAFSDSAELK      1 0.09112961
    ## 11997 XBmRNA11342|XBXL10_1g6076 ETEDINLAYEILEPEPTDIPALR      1 0.08589515
    ## 12619 XBmRNA11342|XBXL10_1g6076               EVYIAMGQK      1 0.02857930
    ##       N_Light  deltaBIC Simple_m0_heavy Simple_m0_ctrl    Simple_kd
    ## 1           2 -2.561797       0.9921021      0.4216439 1.108492e-05
    ## 3792        3 -2.451570       1.0509472      0.4967873 6.616293e-05
    ## 8818        2 -2.530117       1.0516884      0.4330745 8.081394e-05
    ## 11904       2 -2.555062       1.0445287      0.4280146 3.928233e-05
    ## 11997       2 -2.501104       0.9860225      0.4864391 1.000000e-06
    ## 12619       3 -2.639145       1.0051651      0.4062885 1.869287e-05
    ##          Simple_kt Baseline_m0_heavy Baseline_m0_ctrl  Baseline_kd  Baseline_kt
    ## 1     0.0013501033         0.9923337        0.4215875 1.506083e-05 0.0013522318
    ## 3792  0.0003697222         1.0526375        0.4963413 9.289572e-05 0.0003768129
    ## 8818  0.0020428531         1.0535932        0.4325331 1.132653e-04 0.0020646972
    ## 11904 0.0000000000         1.0445458        0.4278729 5.341090e-05 0.0000000000
    ## 11997 0.0000000000         0.9982913        0.4889384 1.000000e-06 0.0000000000
    ## 12619 0.0005038895         1.0049627        0.4062673 1.844505e-05 0.0005037976
    ##       Baseline_phi
    ## 1             0.25
    ## 3792          0.25
    ## 8818          0.25
    ## 11904         0.25
    ## 11997         0.25
    ## 12619         0.00

``` r
max_fit_kd <- log(2)/(tail(avg_heavy_time, 1)*3)
max_fit_hl <- log(2)/max_fit_kd/60
cat("Max Half-life (hrs):\t", max_fit_hl) 
```

    ## Max Half-life (hrs):  212.2

``` r
T5_pep_fits <- T5_pep_fits[!is.na(T5_pep_fits$Simple_kd),]
norm_T5.df <- norm_T5.df[norm_T5.df$Peptide %in% T5_pep_fits$Peptide,]

T6_pep_fits <- T6_pep_fits[!is.na(T6_pep_fits$Simple_kd),]
norm_T6.df <- norm_T6.df[norm_T6.df$Peptide %in% T6_pep_fits$Peptide,]
```

``` r
pep_fit_comp <- T6_pep_fits
pep_fit_comp <- pep_fit_comp[!is.na(pep_fit_comp$deltaBIC),]

pep_fit_comp["Model"] <- rep("Unselected", nrow(pep_fit_comp))
pep_fit_comp[pep_fit_comp$deltaBIC < 0, "Model"] <- "Simple"
pep_fit_comp[pep_fit_comp$deltaBIC > 0, "Model"] <- "Two-pool"
pep_fit_comp[pep_fit_comp$deltaBIC > 10, "deltaBIC"] <- 10

table(pep_fit_comp$Model)
```

    ## 
    ##   Simple Two-pool 
    ##    45986      602

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

![](figures/Frog_2C_SimpModel/Exploratory%20analysis%20not%20included%20in%20manuscript-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
protein_rollup <- function(protein_id, all_reps){

  sub_data <-  all_reps[all_reps$Protein_ID==protein_id,]

  final_calc_lst <- lapply(unique(sub_data$Exp), function(x){
    
    sub_exp <- sub_data[sub_data$Exp == x,]
    sub_exp["Fraction"] <- sub_exp$sum_sn / sum(sub_exp$sum_sn)
  
    sub_exp[light.channels] <- sub_exp[light.channels] * sub_exp$Fraction
    light_avg <- apply(sub_exp[light.channels], 2, sum)
    light_avg <- light_avg / mean(light_avg)
  
    sub_exp[heavy.channels] <- sub_exp[heavy.channels] * sub_exp$Fraction
    heavy_avg <- apply(sub_exp[heavy.channels], 2, sum)
    heavy_avg <- heavy_avg / mean(heavy_avg)
    
    return(c(light_avg, heavy_avg)) })
  
  merged_matrix <- do.call(rbind, final_calc_lst)
  merged_avg <- apply(merged_matrix, 2, mean)

  merged_avg[light.channels] <-
    merged_avg[light.channels] / mean(merged_avg[light.channels])
  merged_avg[heavy.channels] <-
    merged_avg[heavy.channels] / mean(merged_avg[heavy.channels])

  protein_avg <- cbind(data.frame(Protein_ID = protein_id),
                       data.frame(t(merged_avg)))

  return(protein_avg) }


protein_med_merge <- rbind(cbind(norm_T5.df,
                                 data.frame(Exp=rep("T5", nrow(norm_T5.df)))),
                           cbind(norm_T6.df,
                                 data.frame(Exp=rep("T6", nrow(norm_T6.df)))))

#Setup cluster
cl <- makeCluster(num_cores)

#Export necessary variables/functions to the cluster nodes
clusterExport(cl,
              varlist=c("protein_med_merge", "protein_rollup",
                        "light.channels", "heavy.channels"))

protein_med_merge <- parLapply(cl, unique(protein_med_merge$Protein_ID),
       function(x){ protein_rollup(x, protein_med_merge) }) %>% bind_rows

stopCluster(cl) #Stop Cluster

# protein_med_merge <- merge(human_map, protein_med_merge, by="Protein_ID", all.y=TRUE)
# protein_med_merge[is.na(protein_med_merge$Human_Gene), "XLA_Gene"] <- "N/A"
# protein_med_merge[is.na(protein_med_merge$Human_Gene), "Human_Gene"] <- "N/A"

head(protein_med_merge)
```

    ##                   Protein_ID      T0_L        N1        N3        N4        N7
    ## 1  XBmRNA11342|XBXL10_1g6076 0.4208845 0.4785994 0.5448444 0.5970521 0.7701754
    ## 2 XBmRNA73944|XBXL10_1g39228 0.9733531 0.9328750 0.9601760 0.9760373 0.9566142
    ## 3 XBmRNA19289|XBXL10_1g10324 0.9022768 0.9658400 0.9292869 0.9316272 0.9121772
    ## 4 XBmRNA75459|XBXL10_1g40108 0.9874260 1.0954007 0.8733955 0.9715309 0.9146859
    ## 5 XBmRNA57343|XBXL10_1g30445 0.8780641 0.9208206 0.8533880 0.8856882 0.8438882
    ## 6 XBmRNA69844|XBXL10_1g37022 0.9824866 0.9443270 0.9410123 0.9800028 0.9327650
    ##          N9       N11      N12      T0_H        O1        O2       O3        O4
    ## 1 0.9908260 1.7004397 2.497179 1.0131095 1.0265790 1.1017480 1.062352 0.9458617
    ## 2 0.9777540 1.0052733 1.217917 0.9992422 1.0188498 0.9762403 1.007476 1.0169969
    ## 3 0.9267686 1.0580668 1.373957 0.9877894 1.0330499 1.0422463 1.000547 0.9998502
    ## 4 0.9008309 1.1480111 1.108719 1.1208384 1.2488393 1.2581402 1.140366 1.0283165
    ## 5 0.8196696 1.0448356 1.753646 1.0293644 0.9813028 1.0834973 1.104164 1.0408364
    ## 6 0.9867207 0.9905895 1.242096 1.0257464 1.0362610 0.9706270 1.014064 1.0258966
    ##          O5        O7        O9       O10       O11       O12
    ## 1 1.0340756 1.0120790 0.9775329 0.9897608 0.9204845 0.9164171
    ## 2 1.0048174 0.9967572 1.0353408 0.9872515 1.0001491 0.9568787
    ## 3 1.0194296 1.0118611 0.9781498 1.0104561 0.9620053 0.9546154
    ## 4 0.8938632 1.0070439 0.8561329 0.8957045 0.7840280 0.7667273
    ## 5 0.9895954 1.0576162 0.9605185 0.9580319 0.9235430 0.8715296
    ## 6 1.0365603 1.0086087 1.0496626 0.9627951 0.9988585 0.8709193

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

    predict_fit <- simple_heavy_model(avg_heavy_time, merged_values["m0_heavy"],
                                      merged_values["kd"], merged_values["kt"],
                                      sub_pep["Pep_k1"], sub_pep["Pep_k2"])
    
    return(data.frame(t(predict_fit))) })

  pep_predict <- do.call(rbind, pep_predict)
  pep_predict <- apply(pep_predict, 2, median)
  
  out_df <- data.frame(t(pep_predict))
  colnames(out_df) <- paste0("Fit_", heavy.channels)
  
  out_df <- cbind(data.frame(Protein_ID = protein_id,
                             kD = merged_values["kd"],
                             HL_Hrs = log(2)/merged_values["kd"]/60),                             
                             Estimate = estimate_desc,
                             T5_kD = exp1_kD, T6_kD = exp2_kD,
                  out_df)
  
  return(out_df) }

Deg_fits.df <- lapply(protein_med_merge$Protein_ID,
                      function(x) { merge_pep_fits(x, T5_pep_fits, T6_pep_fits) }) %>%
  bind_rows() %>% arrange(desc(kD))

head(Deg_fits.df)
```

    ##                        Protein_ID         kD    HL_Hrs Estimate      T5_kD
    ## kd...1  XBmRNA10569|XBXL10_1g5638 0.05817981 0.1985646 interval 0.05196905
    ## kd...2 XBmRNA35507|XBXL10_1g19214 0.05236107 0.2206306    point         NA
    ## kd...3 XBmRNA34126|XBXL10_1g18541 0.05055962 0.2284917    point         NA
    ## kd...4 XBmRNA40199|XBXL10_1g21682 0.04860139 0.2376980    point 0.04860139
    ## kd...5   XBmRNA1720|XBXL10_1g1100 0.03471313 0.3327977    point         NA
    ## kd...6 XBmRNA20588|XBXL10_1g10943 0.03404127 0.3393661 interval 0.02838204
    ##             T6_kD Fit_T0_H    Fit_O1     Fit_O2       Fit_O3       Fit_O4
    ## kd...1 0.06439057 9.479480 0.5856523 0.06478722 1.460836e-03 3.640090e-05
    ## kd...2 0.05236107 6.415616 0.3198320 0.01576685 3.734582e-05 8.623444e-08
    ## kd...3 0.05055962 9.303371 0.5244588 0.03524892 3.151491e-04 4.712569e-06
    ## kd...4         NA 6.060752 0.5522329 0.09547474 6.934862e-03 5.902776e-04
    ## kd...5 0.03471313 8.540926 1.1795156 0.17803734 6.060285e-03 3.330955e-04
    ## kd...6 0.03970050 8.330511 1.3240351 0.18197992 3.136513e-03 5.290678e-05
    ##              Fit_O5       Fit_O7       Fit_O9      Fit_O10      Fit_O11
    ## kd...1 8.432540e-07 2.267036e-10 1.294815e-19 2.189885e-24 3.521232e-37
    ## kd...2 1.721461e-10 2.091208e-16 8.649746e-32 9.442580e-40 3.965072e-61
    ## kd...3 7.242768e-08 8.266022e-12 5.175940e-22 2.796312e-27 2.108872e-41
    ## kd...4 4.801949e-05 2.005395e-07 1.389713e-13 9.188132e-17 2.742074e-25
    ## kd...5 2.128960e-05 6.018231e-08 1.593733e-14 6.414810e-18 5.060155e-27
    ## kd...6 8.177263e-07 9.076559e-11 5.248636e-21 2.721383e-26 1.838180e-40
    ##             Fit_O12
    ## kd...1 5.233040e-57
    ## kd...2 2.711809e-94
    ## kd...3 2.733373e-63
    ## kd...4 1.681343e-38
    ## kd...5 3.945890e-41
    ## kd...6 2.008412e-62

``` r
Deg_fits.df["deltaBIC"] <- apply(Deg_fits.df, 1, function(row){

  protein_id <- row["Protein_ID"]
  dyn_fit <- as.numeric(row[grepl("Fit", names(row))])
  
  med_data <- protein_med_merge[protein_med_merge$Protein_ID == protein_id,]
  med_data <- as.numeric(med_data[heavy.channels])
  
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

    ##                        Protein_ID         kD    HL_Hrs Estimate      T5_kD
    ## kd...1  XBmRNA10569|XBXL10_1g5638 0.05817981 0.1985646 interval 0.05196905
    ## kd...2 XBmRNA35507|XBXL10_1g19214 0.05236107 0.2206306    point         NA
    ## kd...3 XBmRNA34126|XBXL10_1g18541 0.05055962 0.2284917    point         NA
    ## kd...4 XBmRNA40199|XBXL10_1g21682 0.04860139 0.2376980    point 0.04860139
    ## kd...5   XBmRNA1720|XBXL10_1g1100 0.03471313 0.3327977    point         NA
    ## kd...6 XBmRNA20588|XBXL10_1g10943 0.03404127 0.3393661 interval 0.02838204
    ##             T6_kD Fit_T0_H    Fit_O1     Fit_O2       Fit_O3       Fit_O4
    ## kd...1 0.06439057 9.479480 0.5856523 0.06478722 1.460836e-03 3.640090e-05
    ## kd...2 0.05236107 6.415616 0.3198320 0.01576685 3.734582e-05 8.623444e-08
    ## kd...3 0.05055962 9.303371 0.5244588 0.03524892 3.151491e-04 4.712569e-06
    ## kd...4         NA 6.060752 0.5522329 0.09547474 6.934862e-03 5.902776e-04
    ## kd...5 0.03471313 8.540926 1.1795156 0.17803734 6.060285e-03 3.330955e-04
    ## kd...6 0.03970050 8.330511 1.3240351 0.18197992 3.136513e-03 5.290678e-05
    ##              Fit_O5       Fit_O7       Fit_O9      Fit_O10      Fit_O11
    ## kd...1 8.432540e-07 2.267036e-10 1.294815e-19 2.189885e-24 3.521232e-37
    ## kd...2 1.721461e-10 2.091208e-16 8.649746e-32 9.442580e-40 3.965072e-61
    ## kd...3 7.242768e-08 8.266022e-12 5.175940e-22 2.796312e-27 2.108872e-41
    ## kd...4 4.801949e-05 2.005395e-07 1.389713e-13 9.188132e-17 2.742074e-25
    ## kd...5 2.128960e-05 6.018231e-08 1.593733e-14 6.414810e-18 5.060155e-27
    ## kd...6 8.177263e-07 9.076559e-11 5.248636e-21 2.721383e-26 1.838180e-40
    ##             Fit_O12 deltaBIC
    ## kd...1 5.233040e-57 59.78006
    ## kd...2 2.711809e-94 19.02400
    ## kd...3 2.733373e-63 53.52402
    ## kd...4 1.681343e-38 16.69609
    ## kd...5 3.945890e-41 50.29648
    ## kd...6 2.008412e-62 51.62111

``` r
true_shuffle <- function(x) {
  
  n <- length(x)
  original <- x
  shuffled <- sample(x)

  #Keep re-rolling until NO elements match their starting position
  # -> This prevents correlation for fastest degraders
  #Also, stop from last timepoint ending up with the first one
  # -> This prevents anti-correlation for fastest degraders
  while(any(shuffled == original) || shuffled[n] == original[1]) {
    shuffled <- sample(x) }
  
  return(shuffled) }
```

``` r
final_xla_2C <- merge(human_map, Deg_fits.df, by="Protein_ID", all.y=TRUE)
final_xla_2C[is.na(final_xla_2C$XLA_Gene), "XLA_Gene"] <- "N/A"
final_xla_2C[is.na(final_xla_2C$Human_Gene), "Human_Gene"] <- "N/A"
final_xla_2C <- final_xla_2C[c(colnames(final_xla_2C)[1:8],
                               "deltaBIC")]
final_xla_2C <- merge(final_xla_2C, protein_med_merge, by="Protein_ID")
final_xla_2C <- final_xla_2C %>% arrange(desc(kD))

head(final_xla_2C)
```

    ##                   Protein_ID   XLA_Gene Human_Gene         kD    HL_Hrs
    ## 1  XBmRNA10569|XBXL10_1g5638    ccnb1.S      CCNB1 0.05817981 0.1985646
    ## 2 XBmRNA35507|XBXL10_1g19214  ccnb1.2.L      CCNB1 0.05236107 0.2206306
    ## 3 XBmRNA34126|XBXL10_1g18541  LOC398156      PTTG1 0.05055962 0.2284917
    ## 4 XBmRNA40199|XBXL10_1g21682  ccnb1.2.S      CCNB1 0.04860139 0.2376980
    ## 5   XBmRNA1720|XBXL10_1g1100    ccna2.L      CCNA2 0.03471313 0.3327977
    ## 6 XBmRNA20588|XBXL10_1g10943    ccna1.S      CCNA1 0.03404127 0.3393661
    ##   Estimate      T5_kD      T6_kD deltaBIC      T0_L        N1        N3
    ## 1 interval 0.05196905 0.06439057 59.78006 1.7147274 1.0893752 1.1090663
    ## 2    point         NA 0.05236107 19.02400 0.5715232 0.3167322 0.4122321
    ## 3    point         NA 0.05055962 53.52402 0.7593631 0.2741591 0.3558655
    ## 4    point 0.04860139         NA 16.69609 0.8122257 0.5816428 0.5083152
    ## 5    point         NA 0.03471313 50.29648 0.3719585 0.2490067 0.3498724
    ## 6 interval 0.02838204 0.03970050 51.62111 2.1104195 1.5523841 1.7135697
    ##          N4        N7         N9        N11          N12     T0_H        O1
    ## 1 1.2411665 1.8280703 0.92944429 0.06468336 2.346668e-02 9.479349 0.6087783
    ## 2 0.4270137 1.8856945 1.67675998 1.13389467 1.576150e+00 6.415712 0.3157785
    ## 3 0.4436930 1.9537817 1.95445907 0.69835989 1.560319e+00 9.303156 0.5306161
    ## 4 0.3789010 1.0024334 1.50784145 1.53953923 1.669101e+00 6.061506 0.5321921
    ## 5 0.3315833 1.4536714 1.20658237 1.33815842 2.699167e+00 8.538092 1.2239487
    ## 6 1.7205897 0.8069156 0.07777375 0.01829269 5.503508e-05 8.399209 1.3292753
    ##           O2        O3        O4         O5        O7         O9       O10
    ## 1 0.00000000 0.2046160 0.1690556 0.15700874 0.1380140 0.06531466 0.1457002
    ## 2 0.05760827 0.4069183 0.6042256 0.44871175 0.5156761 0.49982006 0.2410755
    ## 3 0.00000000 0.1234161 0.2365997 0.18184335 0.1569108 0.00000000 0.0000000
    ## 4 0.20740221 0.4029011 0.3367262 0.26181593 0.7716691 0.50081417 0.4388957
    ## 5 0.00000000 0.2604006 0.2245385 0.15593400 0.2258017 0.17048824 0.1304815
    ## 6 0.24473731 0.3145507 0.2879517 0.08938614 0.1435368 0.00000000 0.0000000
    ##          O11       O12
    ## 1 0.03216316 0.0000000
    ## 2 0.74512284 0.7493515
    ## 3 0.18012809 0.2873295
    ## 4 0.48312743 1.0029505
    ## 5 0.07031518 0.0000000
    ## 6 0.11168100 0.0796721

``` r
dyn_model_pep_fit <- function(row, l_time, h_time, heavy_header) {
  
  #--- Model input ----------------------------  
  
  pep_k1 <- as.numeric(row["Pep_k1"])
  pep_k2 <- as.numeric(row["Pep_k2"])

  heavy_data <- as.numeric(row[heavy_header])

  light_cutoff <- as.numeric(row["N_timepoints"])
  light_trim <- l_time[1:light_cutoff]
  sub_light <- as.numeric(row[light.channels])[1:light_cutoff]

  #Input for nlsLM
  all_data <- c(heavy_data, sub_light)
  all_time <- c(h_time, light_trim)

  heavy_n <- length(all_time) - length(light_trim)
  light_n <- length(light_trim)
  
  #--- Starting Guess ----------------------------  
  
  #Starting channel parameters
  m0_heavy_start <- heavy_data[1]
  m0_ctrl_start <- sub_light[1]

  #kD starting parameter 
  se_heavy <- pmax(heavy_data[1:3], 1e-6) #Prevent 0 after log
  se_h_time <- h_time[1:3]
  
  se_fit <- lm(log(se_heavy) ~ se_h_time) #Approximating linearly

  kd_start <- -coef(se_fit)[2]
  kd_start <- max(kd_start, 1e-6)

  #kT starting parameter
  se_light <- c(sub_light[1], pmax(tail(sub_light, 1), 1e-6))
  se_l_time <- c(light_trim[1], tail(light_trim, 1))

  initial_slope <- (se_light[2] - se_light[1]) / (se_l_time[2] - se_l_time[1])
  kt_start <- initial_slope + kd_start * se_light[1]
  kt_start <- max(kt_start, 1e-6)  # keep non-negative
  
  #--- Bounds ----------------------------  
  
  m0_heavy_lower <- max(heavy_data[1] * 0.5, 0.1)
  m0_ctrl_lower <- max(sub_light[1] * 0.5, 0.1)
  kd_upper <- log(2) / 5    # Shortest half-life ~ 5 min
  kt_upper <- kd_upper * heavy_data[1] * 2    
  
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
    simple_pred_values <- c(m0_heavy = NA, m0_ctrl = NA, kd = NA, kt = NA) }
  else {
    simple_pred_values <- coef(simple_fit) }
  
  return(data.frame(t(simple_pred_values))) }
```

``` r
#Setup cluster
cl <- makeCluster(num_cores)

#Export necessary variables/functions to the cluster nodes
clusterExport(cl,
              varlist=c("light.channels", "heavy.channels",
                        "Frog_T5_light.min", "Frog_T5_all.min",
                        "Frog_T6_light.min", "Frog_T6_all.min",
                        "avg_heavy_time",
                        "norm_T5.df", "norm_T6.df", "protein_med_merge",
                        "dyn_model_pep_fit", "true_shuffle",
                        "simple_normal_model",
                        "simple_heavy_model", "simple_pred_model"))
invisible(clusterEvalQ(cl, library(minpack.lm)))

NULL_BIC_values <- parApply(cl, Deg_fits.df, 1, function(row) {

  protein_id <- row["Protein_ID"]
  total_n <- length(heavy.channels)
  flat_predict <- rep(1, total_n)

  median_data <- protein_med_merge[protein_med_merge$Protein_ID == protein_id,]
  T5_check <- is.na(row["T5_kD"])
  T6_check <- is.na(row["T6_kD"])

  real_deltaBIC <- as.numeric(row["deltaBIC"])

  if (T5_check) {
    T6_peps <- norm_T6.df[norm_T6.df$Protein_ID == protein_id,]
    pep_data <- T6_peps[c("Peptide", "Pep_k1", "Pep_k2")] }
  else if (T6_check) {
    T5_peps <- norm_T5.df[norm_T5.df$Protein_ID == protein_id,]
    pep_data <- T5_peps[c("Peptide", "Pep_k1", "Pep_k2")] }
  else {
    T5_peps <- norm_T5.df[norm_T5.df$Protein_ID == protein_id,]
    T6_peps <- norm_T6.df[norm_T6.df$Protein_ID == protein_id,]

    pep_data <- rbind(T5_peps[c("Peptide", "Pep_k1", "Pep_k2")],
                      T6_peps[c("Peptide", "Pep_k1", "Pep_k2")])
    pep_data <- unique(pep_data) }

  NULL_BIC_deltas <- sapply(1:10, function(dummy) { #No difference between 10 and 100

    #-- Fit ------------------------------------------

    # new_header <- heavy.channels
    new_header <- true_shuffle(heavy.channels)
    scrambled_median <- as.numeric(median_data[new_header])

    if (T5_check) {

      T6_fits <- apply(T6_peps, 1, function(x) {
        dyn_model_pep_fit(x, Frog_T6_light.min, Frog_T6_all.min, new_header) } )
      T6_fits <- do.call(rbind, T6_fits)
      merged_values <- apply(T6_fits, 2, median, na.rm = TRUE) }

    else if (T6_check) {

      T5_fits <- apply(T5_peps, 1, function(x) {
        dyn_model_pep_fit(x, Frog_T5_light.min, Frog_T5_all.min, new_header) } )
      T5_fits <- do.call(rbind, T5_fits)
      merged_values <- apply(T5_fits, 2, median, na.rm = TRUE) }

    else {

      T5_fits <- apply(T5_peps, 1, function(x) {
        dyn_model_pep_fit(x, Frog_T5_light.min, Frog_T5_all.min, new_header) } )
      T5_fits <- do.call(rbind, T5_fits)
      T5_fits <- apply(T5_fits, 2, median, na.rm = TRUE)

      T6_fits <- apply(T6_peps, 1, function(x) {
        dyn_model_pep_fit(x, Frog_T6_light.min, Frog_T6_all.min, new_header) } )
      T6_fits <- do.call(rbind, T6_fits)
      T6_fits <- apply(T6_fits, 2, median, na.rm = TRUE)

      merged_values <- (T5_fits + T6_fits)/2 }

    #-- Prediction ------------------------------------------
    if (any(is.na(merged_values))) {
      deltaBIC <- (-100)
      return(deltaBIC) }
    else { } # ... normal BIC computation ...

    pep_predict <- apply(pep_data[2:3], 1, function(sub_pep){

      predict_fit <- simple_heavy_model(avg_heavy_time, merged_values["m0_heavy"],
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

    ##  [1] -11.019303 -10.180922 -10.271292 -10.495663 -10.965415 -10.979382
    ##  [7] -11.019303 -10.385405 -11.019303 -10.972231  -9.623352  -8.593963
    ## [13]  -8.308638 -10.051262  -8.777796 -10.899971  -9.898571  -9.497135
    ## [19]  -8.579501 -10.136944

``` r
stopCluster(cl) #Stop Cluster
```

``` r
real_scores <- Deg_fits.df$deltaBIC

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

    ## [1] "The 5% FDR Threshold is Delta BIC > -5.70000000000001"

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
# tiff("frog_2C_deltaBIC_hist_10variations.tiff", units="in",
#      width=6, height=4.5, res=300)

hist_deltaBIC
```

![](figures/Frog_2C_SimpModel/BIC%20cutoff%20plot-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
final_BIC_threshold <- 2 #Negative doesn't make sense

final_xla_2C["mClass"] <- rep("Flat", nrow(final_xla_2C))
final_xla_2C[(final_xla_2C$deltaBIC > final_BIC_threshold),
             "mClass"] <- "Deg"
final_xla_2C[final_xla_2C$kD < max_fit_kd, "HL_Hrs"] <- max_fit_hl

head(final_xla_2C)
```

    ##                   Protein_ID   XLA_Gene Human_Gene         kD    HL_Hrs
    ## 1  XBmRNA10569|XBXL10_1g5638    ccnb1.S      CCNB1 0.05817981 0.1985646
    ## 2 XBmRNA35507|XBXL10_1g19214  ccnb1.2.L      CCNB1 0.05236107 0.2206306
    ## 3 XBmRNA34126|XBXL10_1g18541  LOC398156      PTTG1 0.05055962 0.2284917
    ## 4 XBmRNA40199|XBXL10_1g21682  ccnb1.2.S      CCNB1 0.04860139 0.2376980
    ## 5   XBmRNA1720|XBXL10_1g1100    ccna2.L      CCNA2 0.03471313 0.3327977
    ## 6 XBmRNA20588|XBXL10_1g10943    ccna1.S      CCNA1 0.03404127 0.3393661
    ##   Estimate      T5_kD      T6_kD deltaBIC      T0_L        N1        N3
    ## 1 interval 0.05196905 0.06439057 59.78006 1.7147274 1.0893752 1.1090663
    ## 2    point         NA 0.05236107 19.02400 0.5715232 0.3167322 0.4122321
    ## 3    point         NA 0.05055962 53.52402 0.7593631 0.2741591 0.3558655
    ## 4    point 0.04860139         NA 16.69609 0.8122257 0.5816428 0.5083152
    ## 5    point         NA 0.03471313 50.29648 0.3719585 0.2490067 0.3498724
    ## 6 interval 0.02838204 0.03970050 51.62111 2.1104195 1.5523841 1.7135697
    ##          N4        N7         N9        N11          N12     T0_H        O1
    ## 1 1.2411665 1.8280703 0.92944429 0.06468336 2.346668e-02 9.479349 0.6087783
    ## 2 0.4270137 1.8856945 1.67675998 1.13389467 1.576150e+00 6.415712 0.3157785
    ## 3 0.4436930 1.9537817 1.95445907 0.69835989 1.560319e+00 9.303156 0.5306161
    ## 4 0.3789010 1.0024334 1.50784145 1.53953923 1.669101e+00 6.061506 0.5321921
    ## 5 0.3315833 1.4536714 1.20658237 1.33815842 2.699167e+00 8.538092 1.2239487
    ## 6 1.7205897 0.8069156 0.07777375 0.01829269 5.503508e-05 8.399209 1.3292753
    ##           O2        O3        O4         O5        O7         O9       O10
    ## 1 0.00000000 0.2046160 0.1690556 0.15700874 0.1380140 0.06531466 0.1457002
    ## 2 0.05760827 0.4069183 0.6042256 0.44871175 0.5156761 0.49982006 0.2410755
    ## 3 0.00000000 0.1234161 0.2365997 0.18184335 0.1569108 0.00000000 0.0000000
    ## 4 0.20740221 0.4029011 0.3367262 0.26181593 0.7716691 0.50081417 0.4388957
    ## 5 0.00000000 0.2604006 0.2245385 0.15593400 0.2258017 0.17048824 0.1304815
    ## 6 0.24473731 0.3145507 0.2879517 0.08938614 0.1435368 0.00000000 0.0000000
    ##          O11       O12 mClass
    ## 1 0.03216316 0.0000000    Deg
    ## 2 0.74512284 0.7493515    Deg
    ## 3 0.18012809 0.2873295    Deg
    ## 4 0.48312743 1.0029505    Deg
    ## 5 0.07031518 0.0000000    Deg
    ## 6 0.11168100 0.0796721    Deg

``` r
ss_bp.all <- final_xla_2C

ss_bp.all <- ss_bp.all[!is.na(ss_bp.all$T5_kD),]
ss_bp.all <- ss_bp.all[!is.na(ss_bp.all$T6_kD),]

ss_bp.all["Avg_kD"] <- (ss_bp.all$T5_kD + ss_bp.all$T6_kD)/2

ss_bp.all["Log10_T5"] <- log10(ss_bp.all$T5_kD)
ss_bp.all["Log10_T6"] <- log10(ss_bp.all$T6_kD)
ss_bp.all[ss_bp.all$Log10_T5<(-8), "Log10_T5"] <- (-8)
ss_bp.all[ss_bp.all$Log10_T6<(-8), "Log10_T6"] <- (-8)

cor(ss_bp.all[ss_bp.all$mClass=="Deg",]$Log10_T5,
    ss_bp.all[ss_bp.all$mClass=="Deg",]$Log10_T6)^2
```

    ## [1] 0.8213983

``` r
cor(ss_bp.all[ss_bp.all$mClass=="Deg" &
                ss_bp.all$Avg_kD > max_fit_kd,]$Log10_T5,
    ss_bp.all[ss_bp.all$mClass=="Deg" &
                ss_bp.all$Avg_kD > max_fit_kd,]$Log10_T6)^2
```

    ## [1] 0.8205654

``` r
p1 <- ggplot() +
  # geom_point(data=ss_bp.all,
  #            aes(x=Log10_T5, y=Log10_T6, color=mClass),
  #            shape=1, size=3, stroke=1, alpha=0.5) +
  geom_point(data=ss_bp.all[ss_bp.all$mClass == "Deg",],
             aes(x=Log10_T5, y=Log10_T6),
             shape=1, size=3, stroke=1, alpha=0.5) +  
  
  geom_hline(yintercept=log10(max_fit_kd), size=1, linetype="dashed") +
  geom_vline(xintercept=log10(max_fit_kd), size=1, linetype="dashed") +

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
# tiff("2C_kD_biplot.tiff", units="in",
#      width=5, height=5, res=300)

p1
```

![](figures/Frog_2C_SimpModel/Degradation%20rate%20biplot%20for%20degrading%20proteins-1.png)<!-- -->

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
  
  all.kresult[light.channels] <- apply(all.kresult[light.channels], 1,
                                       function(row){
    row <- row / mean(row)
    row <- row / row[1]
    return(data.frame(t(row))) }) %>% bind_rows

  all.kresult[heavy.channels] <- apply(all.kresult[heavy.channels], 1,
                                       function(row){
    row <- row / mean(row)
    row <- row / row[1]
    return(data.frame(t(row))) }) %>% bind_rows
  
  all.kresult <- all.kresult %>%
    pivot_longer(cols = -c(Cluster, Size),
                 names_to = "Time", values_to = "Value")
  
  whole_t_vec <- c((Frog_T5_light.min + Frog_T6_light.min)/2,
                   (Frog_T5_all.min + Frog_T6_all.min)/2)
  whole_t_vec <- whole_t_vec + ((161+150)/2)
  
  all.kresult["Time_Hrs"] <- rep(whole_t_vec/60,
                                 nrow(all.kresult) /
                                   length(whole_t_vec))


  plots <- lapply(seq_len(k), function(i){
    sub_cl <- all.kresult[all.kresult$Cluster == i,]

    light.df <- sub_cl[sub_cl$Time %in% light.channels,]
    heavy.df <- sub_cl[sub_cl$Time %in% heavy.channels,]

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

export_kmeans <- kmeans_plot(protein_med_merge, order=c(4,1,5,3,2))
```

![](figures/Frog_2C_SimpModel/kmeans%20plot%20of%20protein%20averages-1.png)<!-- -->

``` r
# kmeans_plot(final_2C_proteins[final_2C_proteins$Protein_ID %in% norm_set,], k=3)

# dev.off() #Uncomment for new image!
```

``` r
plot_example_protein <- function(protein_id) {

  sub_df <- final_xla_2C[final_xla_2C$Protein_ID == protein_id,]
  print(sub_df$HL_Hrs)
  print(sub_df$Human_Gene)
  
  sub_df <- sub_df[c("Protein_ID", light.channels, heavy.channels)]

  sub_df[light.channels] <- sub_df[light.channels] / sub_df[[light.channels[1]]]
  sub_df[heavy.channels] <- sub_df[heavy.channels] / sub_df[[heavy.channels[1]]]
    
  sub_df <- sub_df %>%
    pivot_longer(cols=-Protein_ID)
  sub_df["Label"] <- c(rep("Control", length(light.channels)),
                       rep("Label", length(heavy.channels)))
  sub_df["Label"] <- factor(sub_df$Label, levels=c("Label", "Control"))
  
  sub_df["Time"] <- c(avg_light_time, avg_heavy_time) 
  sub_df["Time"] <- sub_df$Time + (161+150)/2
  sub_df["Time"] <- sub_df$Time/60

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

plot_example_protein("XBmRNA51997|XBXL10_1g27704")
```

    ## [1] 1.233804
    ## [1] "GMNN"

![](figures/Frog_2C_SimpModel/Example%20proteins%20in%20manuscript-1.png)<!-- -->

``` r
plot_example_protein("XBmRNA56795|XBXL10_1g30178")
```

    ## [1] 212.2
    ## [1] "GAPDH"

![](figures/Frog_2C_SimpModel/Example%20proteins%20in%20manuscript-2.png)<!-- -->

``` r
plot_example_protein("XBmRNA79044|XBXL10_1g42043")
```

    ## [1] 0.8130045
    ## [1] "KIF22"

![](figures/Frog_2C_SimpModel/Example%20proteins%20in%20manuscript-3.png)<!-- -->

``` r
plot_example_protein("XBmRNA50813|XBXL10_1g27137")
```

    ## [1] 4.916948
    ## [1] "SGO1"

![](figures/Frog_2C_SimpModel/Example%20proteins%20in%20manuscript-4.png)<!-- -->

``` r
plot_example_protein("XBmRNA20427|XBXL10_1g10860")
```

    ## [1] 6.093857
    ## [1] "ESPL1"

![](figures/Frog_2C_SimpModel/Example%20proteins%20in%20manuscript-5.png)<!-- -->

``` r
plot_example_protein("XBmRNA71859|XBXL10_1g38013")
```

    ## [1] 8.316108
    ## [1] "SOX3"

![](figures/Frog_2C_SimpModel/Example%20proteins%20in%20manuscript-6.png)<!-- -->

``` r
plot_example_protein("XBmRNA46716|XBXL10_1g25041")
```

    ## [1] 21.69022
    ## [1] "FBXO11"

![](figures/Frog_2C_SimpModel/Example%20proteins%20in%20manuscript-7.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
frog_kd_hist <- final_xla_2C[c("Protein_ID", "kD", "mClass")]
frog_kd_hist["Plot_kD"] <- frog_kd_hist$kD
frog_kd_hist[frog_kd_hist$kD > 0.001, "Plot_kD"] <- 0.001

frog_kd_hist$mClass <- factor(frog_kd_hist$mClass,
                              levels = c("Flat", "Deg"))

#Define exact bin edges to guarantee the bars and the outline align perfectly
bw <- 0.000025
brks <- seq(1e-6, 0.001 + bw, by = bw)

#Calculate half-life x-intercepts (kd = ln(2) / minutes)
kd_1day <- log(2) / (24 * 60)
kd_12hr <- log(2) / (12 * 60)
kd_6hr  <- log(2) / (6 * 60)

frog_kd_pdf <- ggplot(frog_kd_hist, aes(x = Plot_kD, fill = mClass,
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

  geom_vline(xintercept = log(2) / (70 * 60),
             linetype = "dashed", color = "red3",
             linewidth = 1) +

  scale_fill_manual(values = c("Flat" = "grey55", "Deg" = "#5E4FA2")) + 
  scale_color_manual(values = c("Flat" = NA, "Deg" = "#5E4FA2")) +
  scale_alpha_manual(values = c("Flat" = 0.9, "Deg" = 0.7)) +
  
  # LABELS & THEME: Matching your reference image
  labs(x = expression("k"["d"]~"(min"^{-1}*") | Proteins = 8,806"), y = "Density") +
  theme_bw() +
  theme(axis.text=element_text(size=21,colour="black"),
        plot.title = element_text(size=24, hjust=0.5),
        axis.title = element_text(size=22),
        panel.border = element_rect(linewidth=1.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
        legend.position = "none",
        plot.margin = margin(t = 10, r = 25, b = 10,
                             l = 10, unit = "pt")) +
  coord_cartesian(xlim=c(1e-6, 0.001)) +
  scale_x_continuous(breaks=c(0, 0.00025, 0.0005, 0.00075, 0.001),
                     labels=c("0", "0.00025", "0.0005",
                              "0.00075", "0.001"))

# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=9, height=4, res=300)

frog_kd_pdf
```

    ## Warning: Removed 40 rows containing missing values or values outside the scale range
    ## (`geom_step()`).

![](figures/Frog_2C_SimpModel/Frog%20kd%20distribution-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
plot_data <- final_xla_2C[final_xla_2C$mClass == "Deg",]
plot_data["Label"] <- rep("Base", nrow(plot_data))
plot_data[plot_data$kD < max_fit_kd, "Label"] <- "Collapsed"
plot_data[plot_data$kD < max_fit_kd, "HL_Hrs"] <- max_fit_hl
plot_data["Label"] <- factor(plot_data$Label, levels = c("Collapsed", "Base"))

median_HL <- median(plot_data$HL_Hrs, na.rm = TRUE)

hl_hist <- ggplot(data = plot_data, aes(x = HL_Hrs, fill=Label, alpha=Label)) +
  geom_histogram(bins = 20, color = "black", linewidth = 0.5) +
  geom_vline(xintercept = median_HL, color = "red3", linetype = "dashed", linewidth = 1.2) +
  labs(x = "Half-life (Hours)",
       y = "Number of proteins") +
  theme_bw() +
  theme(axis.text=element_text(size=21,colour="black"),
        plot.title = element_text(size=22, hjust=0.5),
        axis.title = element_text(size=22),
        panel.border = element_rect(linewidth=1.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
        legend.position = "none",
        plot.margin = margin(t = 10, r = 25, b = 10,
                             l = 10, unit = "pt")) +
  scale_fill_manual(values=c("Collapsed" = "#9ECAE1", "Base" = "#5E4FA2")) +
  scale_alpha_manual(values=c("Collapsed" = 1, "Base" = 0.7))


# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=9, height=4, res=300)

hl_hist
```

![](figures/Frog_2C_SimpModel/Frog%20HL%20distribution%20for%20degrading%20proteins-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
peptide_count <- data.frame(table(rbind(norm_T5.df[1:2],
                                        norm_T6.df[1:2])$Protein_ID))
colnames(peptide_count) <- c("Protein_ID", "Number_Pep_Fits")

final_xla_2C <- merge(final_xla_2C, peptide_count, by="Protein_ID") %>% 
  arrange(desc(kD))

head(final_xla_2C)
```

    ##                   Protein_ID   XLA_Gene Human_Gene         kD    HL_Hrs
    ## 1  XBmRNA10569|XBXL10_1g5638    ccnb1.S      CCNB1 0.05817981 0.1985646
    ## 2 XBmRNA35507|XBXL10_1g19214  ccnb1.2.L      CCNB1 0.05236107 0.2206306
    ## 3 XBmRNA34126|XBXL10_1g18541  LOC398156      PTTG1 0.05055962 0.2284917
    ## 4 XBmRNA40199|XBXL10_1g21682  ccnb1.2.S      CCNB1 0.04860139 0.2376980
    ## 5   XBmRNA1720|XBXL10_1g1100    ccna2.L      CCNA2 0.03471313 0.3327977
    ## 6 XBmRNA20588|XBXL10_1g10943    ccna1.S      CCNA1 0.03404127 0.3393661
    ##   Estimate      T5_kD      T6_kD deltaBIC      T0_L        N1        N3
    ## 1 interval 0.05196905 0.06439057 59.78006 1.7147274 1.0893752 1.1090663
    ## 2    point         NA 0.05236107 19.02400 0.5715232 0.3167322 0.4122321
    ## 3    point         NA 0.05055962 53.52402 0.7593631 0.2741591 0.3558655
    ## 4    point 0.04860139         NA 16.69609 0.8122257 0.5816428 0.5083152
    ## 5    point         NA 0.03471313 50.29648 0.3719585 0.2490067 0.3498724
    ## 6 interval 0.02838204 0.03970050 51.62111 2.1104195 1.5523841 1.7135697
    ##          N4        N7         N9        N11          N12     T0_H        O1
    ## 1 1.2411665 1.8280703 0.92944429 0.06468336 2.346668e-02 9.479349 0.6087783
    ## 2 0.4270137 1.8856945 1.67675998 1.13389467 1.576150e+00 6.415712 0.3157785
    ## 3 0.4436930 1.9537817 1.95445907 0.69835989 1.560319e+00 9.303156 0.5306161
    ## 4 0.3789010 1.0024334 1.50784145 1.53953923 1.669101e+00 6.061506 0.5321921
    ## 5 0.3315833 1.4536714 1.20658237 1.33815842 2.699167e+00 8.538092 1.2239487
    ## 6 1.7205897 0.8069156 0.07777375 0.01829269 5.503508e-05 8.399209 1.3292753
    ##           O2        O3        O4         O5        O7         O9       O10
    ## 1 0.00000000 0.2046160 0.1690556 0.15700874 0.1380140 0.06531466 0.1457002
    ## 2 0.05760827 0.4069183 0.6042256 0.44871175 0.5156761 0.49982006 0.2410755
    ## 3 0.00000000 0.1234161 0.2365997 0.18184335 0.1569108 0.00000000 0.0000000
    ## 4 0.20740221 0.4029011 0.3367262 0.26181593 0.7716691 0.50081417 0.4388957
    ## 5 0.00000000 0.2604006 0.2245385 0.15593400 0.2258017 0.17048824 0.1304815
    ## 6 0.24473731 0.3145507 0.2879517 0.08938614 0.1435368 0.00000000 0.0000000
    ##          O11       O12 mClass Number_Pep_Fits
    ## 1 0.03216316 0.0000000    Deg               2
    ## 2 0.74512284 0.7493515    Deg               1
    ## 3 0.18012809 0.2873295    Deg               1
    ## 4 0.48312743 1.0029505    Deg               1
    ## 5 0.07031518 0.0000000    Deg               1
    ## 6 0.11168100 0.0796721    Deg               3

``` r
supp_df <- merge(protein_annotation, final_xla_2C,
                 by="Protein_ID", all.y=TRUE) %>%
  arrange(desc(kD)) %>%
  dplyr::select(Protein_ID, XLA_Gene, Human_Gene,
                Description, everything())

supp_df["Light_FC"] <- apply(supp_df, 1, function(row){

  light_fc <- as.numeric(row["N12"]) / as.numeric(row["T0_L"])

  return(light_fc) })

inc_min_fc <- quantile(supp_df[supp_df$Protein_ID %in% norm_set,]$Light_FC, 0.99)
dec_min_fc <- quantile(supp_df[supp_df$Protein_ID %in% norm_set,]$Light_FC, 0.01)

#--------------------------------------------

supp_df["Control_Type"] <- rep("Modify", nrow(supp_df))
supp_df[supp_df$Light_FC > inc_min_fc, "Control_Type"] <- "Increasing"
supp_df[supp_df$Light_FC < dec_min_fc, "Control_Type"] <- "Decreasing"
supp_df[(supp_df$Light_FC < inc_min_fc) &
          (supp_df$Light_FC > dec_min_fc), "Control_Type"] <- "Unchanging"

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
  c("Protein ID", "XenBase Annotation", "Human Protein", "Description",
    "Fitted kD", "Reported Half-life (Hrs)", "Estimate Confidence",
    "Replicate 1 kD", "Replicate 2 kD", "BIC against NULL",
    paste0("Control: ", round((avg_light_time+ ((161+150)/2))/60, 1), " hpf"),
    paste0("18O: ", round((avg_heavy_time + ((161+150)/2))/60, 1), " hpf"),
    "Peptides Quantified", "Control Type", "Is Degrading?",
    "Grouped Class")

supp_df[(supp_df$`Is Degrading?`=="No") &
          (supp_df$`Reported Half-life (Hrs)` < max_fit_hl),
        "Reported Half-life (Hrs)"] <- NA
supp_df["Grouped Class"] <- NULL
supp_df <- rbind(supp_df[!is.na(supp_df$`Reported Half-life (Hrs)`),],
                 supp_df[is.na(supp_df$`Reported Half-life (Hrs)`),])

dir.create("Files/Supp_Tables", showWarnings = FALSE, recursive = TRUE)
write.csv(supp_df, "Files/Supp_Tables/Frog_2-cell_Supplementary_Table.csv",
          row.names = FALSE)
```

``` r
write.csv(final_xla_2C,
          "Files/Fits/Dyn-Model/XLA_O18_T5-6_FinalFits_ALL.csv",
          row.names = FALSE)
```

``` r
norm_export <- supp_df[supp_df$`Protein ID` %in% norm_set,]

set.seed(123)
norm_export <- norm_export[sample(nrow(norm_export)),][1:3]
norm_export <- merge(norm_export, protein_annotation,
                     by.x="Protein ID", by.y="Protein_ID")

dir.create("Files/Supp_Tables", showWarnings = FALSE, recursive = TRUE)
write.csv(norm_export, "Files/Supp_Tables/Norm-Proteins_Supplementary_Table.csv",
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
    ## [1] minpack.lm_1.2-4 patchwork_1.3.2  ggplot2_4.0.3    tidyr_1.3.2     
    ## [5] dplyr_1.2.0     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] vctrs_0.7.1        cli_3.6.5          knitr_1.51         rlang_1.1.7       
    ##  [5] xfun_0.58          otel_0.2.0         purrr_1.2.2        generics_0.1.4    
    ##  [9] S7_0.2.1           labeling_0.4.3     glue_1.8.0         htmltools_0.5.9   
    ## [13] scales_1.4.0       rmarkdown_2.31     grid_4.5.3         evaluate_1.0.5    
    ## [17] tibble_3.3.1       fastmap_1.2.0      yaml_2.3.12        lifecycle_1.0.5   
    ## [21] compiler_4.5.3     RColorBrewer_1.1-3 pkgconfig_2.0.3    rstudioapi_0.19.0 
    ## [25] farver_2.1.2       digest_0.6.39      R6_2.6.1           tidyselect_1.2.1  
    ## [29] pillar_1.11.1      magrittr_2.0.4     withr_3.0.3        tools_4.5.3       
    ## [33] gtable_0.3.6
