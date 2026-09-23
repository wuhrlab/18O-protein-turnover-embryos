Frog Gast SimpModel
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
library(ggplot2)
library(patchwork)
library(minpack.lm)
library(parallel)

num_cores <- detectCores() - 1  # Save one core for system stability
knitr::opts_chunk$set(fig.path = "figures/Frog_Gast_SimpModel/")
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

light.channels <- c("T0", "N2", "N4", "N6", "N8", "N10", "N11")
heavy.channels <- c("T0", "O1", "O2", "O3", "O4", "O5",
                    "O6", "O7", "O8", "O9", "O10", "O11")

Frog_T9_all.min <- c(1472, 1502, 1535, 1592, 1714, 1832,
                     2412, 2934, 3270, 3749, 4260, 4411) - 1472
Frog_T9_light.min <- Frog_T9_all.min[c(1,3,5,7,9,11:12)]

#We're run in parallel from different mothers
avg_light_time <- Frog_T9_light.min
avg_heavy_time <- Frog_T9_all.min

#Files to normalize
XLA_O18_T9  <- read.csv("Data/XLA_O18/XLA_O18-T9_Filtered_Data.csv")
XLA_O18_T10  <- read.csv("Data/XLA_O18/XLA_O18-T10_Filtered_Data.csv")

norm_set <- read.csv("Data/XLA_Norm/XLA-O18_YolkNormSet_T8-NYS-Decay.csv")$Protein_ID
```

``` r
#Removing things we can't say anything about -> human keratin contam. + frog spinout
XLA_O18_T9 <- XLA_O18_T9[!XLA_O18_T9$Protein_ID %in% c(frog_yolk_set, keratins),]
XLA_O18_T10 <- XLA_O18_T10[!XLA_O18_T10$Protein_ID %in% c(frog_yolk_set, keratins),]

XLA_O18_T9 <- XLA_O18_T9[!rowSums(XLA_O18_T9[heavy.channels]) == 0,]
XLA_O18_T10 <- XLA_O18_T10[!rowSums(XLA_O18_T10[heavy.channels]) == 0,]

XLA_O18_T9["StartPool"] <- apply(XLA_O18_T9[heavy.channels], 1,
      function(row){
        starting_pool <- sum(row[c("T0", "O1")]) / sum(row)
        return(starting_pool) })

XLA_O18_T9["HeavyWeight"] <- apply(XLA_O18_T9[c(light.channels, heavy.channels)], 1,
      function(row){

        light_sn <- row[1:length(light.channels)]
        heavy_sn <- row[(length(light.channels)+1):length(row)]
        
        return(mean(light_sn)/(mean(heavy_sn) + 1e-6)) })

XLA_O18_T10["StartPool"] <- apply(XLA_O18_T10[heavy.channels], 1,
      function(row){
        starting_pool <- sum(row[c("T0", "O1")]) / sum(row)
        return(starting_pool) })

XLA_O18_T10["HeavyWeight"] <- apply(XLA_O18_T10[c(light.channels, heavy.channels)], 1,
      function(row){

        light_sn <- row[1:length(light.channels)]
        heavy_sn <- row[(length(light.channels)+1):length(row)]
        
        return(mean(light_sn)/(mean(heavy_sn) + 1e-6)) })


#This filter sets a floor to quantify based on signal
# -> Requires (T0 + O1)/11 > 0.15
# -> Requires heavy signal to not outweigh light signal
XLA_O18_T9 <- XLA_O18_T9[(XLA_O18_T9$StartPool > 0.15) &
                           (XLA_O18_T9$HeavyWeight > 0.8),]
XLA_O18_T10 <- XLA_O18_T10[(XLA_O18_T10$StartPool > 0.15) &
                             (XLA_O18_T10$HeavyWeight > 0.8),]
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
  
XLA_O18_T9 <- split_experiments(XLA_O18_T9)
XLA_O18_T10 <- split_experiments(XLA_O18_T10)

light.channels <- c("T0_L", "N2", "N4", "N6", "N8", "N10", "N11")
heavy.channels <- c("T0_H", "O1", "O2", "O3", "O4", "O5",
                    "O6", "O7", "O8", "O9", "O10", "O11")
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

norm_T9.df <- norm_Frog_split(XLA_O18_T9, L_error_T9, H_error_T9)
```

    ##      T0_L        N2        N4        N6        N8       N10       N11      T0_H 
    ## 1.0462468 1.0893495 1.1479329 1.0626476 0.9987278 0.8361145 0.8189810 1.0326380 
    ##        O1        O2        O3        O4        O5        O6        O7        O8 
    ## 1.1393531 1.0555069 1.0392099 1.0484202 1.1497304 1.0283104 1.0536524 0.8324753 
    ##        O9       O10       O11 
    ## 0.9215045 0.8222788 0.8769202

``` r
norm_T9.df <- norm_T9.df[norm_T9.df$T0_H > 0.8,]
norm_T9.df <- norm_T9.df[norm_T9.df$T0_L > 0.15,]

norm_T10.df <- norm_Frog_split(XLA_O18_T10, L_error_T10, H_error_T10)
```

    ##      T0_L        N2        N4        N6        N8       N10       N11      T0_H 
    ## 1.1361189 1.0385251 1.0580857 1.0471988 0.9227267 0.8977043 0.8996406 1.0867109 
    ##        O1        O2        O3        O4        O5        O6        O7        O8 
    ## 1.2540582 1.0677327 0.9889145 1.0533855 1.1013085 0.9949280 0.9969123 0.8512694 
    ##        O9       O10       O11 
    ## 0.8965379 0.8585746 0.8496674

``` r
norm_T10.df <- norm_T10.df[norm_T10.df$T0_H > 0.8,]
norm_T10.df <- norm_T10.df[norm_T10.df$T0_L > 0.15,]
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

AA_pep_T9 <- norm_T9.df[c("Peptide", "Pep_k1", "Pep_k2")]
AA_pep_T9 <- apply(AA_pep_T9, 1, function(x) {
  AA_peptide_cutoff(x, Frog_T9_light.min)}) %>% bind_rows
row.names(AA_pep_T9) <- NULL

AA_pep_T10 <- norm_T10.df[c("Peptide", "Pep_k1", "Pep_k2")]
AA_pep_T10 <- apply(AA_pep_T10, 1, function(x) {
  AA_peptide_cutoff(x, Frog_T9_light.min)}) %>% bind_rows
row.names(AA_pep_T10) <- NULL

norm_T9.df <- merge(norm_T9.df, AA_pep_T9, by="Peptide")
norm_T10.df <- merge(norm_T10.df, AA_pep_T10, by="Peptide")
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

  print(row)
  print(m0_heavy_lower)
  cat("\n\n")
  
  plot(h_time, all_data[1:12], ylim=c(0,3))
  lines(h_time, simple_predict[1:12], col="red")
  lines(h_time, baseline_predict[1:12], col="blue")

  plot(l_time, as.numeric(row[light.channels]), ylim=c(0,3))
  lines(light_trim, simple_predict[13:length(all_data)], col="red", lty="dashed")
  lines(light_trim, baseline_predict[13:length(all_data)], col="blue", lty="dashed")

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

# fit_k_rates("XBmRNA1670|XBXL10_1g1076", norm_T9.df, Frog_T9_light.min, Frog_T9_all.min)
# fit_k_rates("XBmRNA1670|XBXL10_1g1076", norm_T10.df, Frog_T9_light.min, Frog_T9_all.min)

# fit_k_rates("XBmRNA4619|XBXL10_1g2646", norm_T9.df, Frog_T9_light.min, Frog_T9_all.min)
# fit_k_rates("XBmRNA4619|XBXL10_1g2646", norm_T10.df, Frog_T9_light.min, Frog_T9_all.min)

#Setup cluster
cl <- makeCluster(num_cores)

#Export necessary variables/functions to the cluster nodes
clusterExport(cl,
              varlist=c("light.channels", "heavy.channels",
                        "Frog_T9_light.min", "Frog_T9_all.min",
                        "norm_T9.df", "norm_T10.df", "fit_k_rates",
                        "peptide_fit", "pep_decay_rate", "simple_normal_model",
                        "simple_heavy_model", "simple_pred_model",
                        "twopool_pred_model",
                        "twopool_normal_model", "twopool_heavy_model"))
invisible(clusterEvalQ(cl, library(minpack.lm)))

T9_pep_fits <- parLapply(cl, unique(norm_T9.df$Protein_ID), function(x){
  fit_k_rates(x, norm_T9.df, Frog_T9_light.min, Frog_T9_all.min) }) %>%
  bind_rows

T10_pep_fits <- parLapply(cl, unique(norm_T10.df$Protein_ID), function(x){
  fit_k_rates(x, norm_T10.df, Frog_T9_light.min, Frog_T9_all.min) }) %>%
  bind_rows

stopCluster(cl) #Stop Cluster

head(T9_pep_fits)
```

    ##                       Protein_ID             Peptide Pep_k1     Pep_k2 N_Light
    ## 1     XBmRNA73944|XBXL10_1g39228    AAAADGMEQMEMDESR      1 0.22193346       2
    ## 8408  XBmRNA73944|XBXL10_1g39228  EELQLLQEQGSYVGEVVR      1 0.14049298       2
    ## 9120  XBmRNA73944|XBXL10_1g39228 EHAPSIIFMDEIDSIGSSR      1 0.24647935       2
    ## 9760  XBmRNA73944|XBXL10_1g39228             ELFVMAR      1 0.06383364       2
    ## 11683 XBmRNA73944|XBXL10_1g39228            EVIELPVK      1 0.05804911       2
    ## 12877 XBmRNA73944|XBXL10_1g39228             FIGEGAR      1 0.04925093       2
    ##        deltaBIC Simple_m0_heavy Simple_m0_ctrl    Simple_kd   Simple_kt
    ## 1     -2.575996       0.9870011      0.7426055 1.000000e-06 0.002199183
    ## 8408  -2.602970       0.9919972      0.9417978 1.000000e-06 0.000994895
    ## 9120  -2.647707       1.0027575      0.9714300 1.310254e-06 0.000000000
    ## 9760  -2.639208       0.9649619      0.8520336 2.573130e-06 0.002675317
    ## 11683 -2.676742       1.0218185      1.0410635 1.653023e-05 0.000000000
    ## 12877 -2.639097       1.0181603      1.0298985 1.653631e-05 0.000000000
    ##       Baseline_m0_heavy Baseline_m0_ctrl  Baseline_kd  Baseline_kt Baseline_phi
    ## 1             0.9882686        0.7425166 1.000000e-06 0.0021988328         0.25
    ## 8408          0.9918077        0.9418406 1.000000e-06 0.0009930859         0.25
    ## 9120          1.0057776        0.9699627 4.357253e-06 0.0000000000         0.25
    ## 9760          0.9649509        0.8520359 3.373275e-06 0.0026752705         0.25
    ## 11683         1.0214417        1.0422838 2.118309e-05 0.0000000000         0.25
    ## 12877         1.0178553        1.0297961 1.627868e-05 0.0000000000         0.00

``` r
head(T10_pep_fits)
```

    ##                       Protein_ID             Peptide Pep_k1     Pep_k2 N_Light
    ## 1     XBmRNA31787|XBXL10_1g17150           AAAAAGAAK      1 0.05286154       2
    ## 16860 XBmRNA31787|XBXL10_1g17150           GSGASGSFK      1 0.13279410       2
    ## 31702 XBmRNA31787|XBXL10_1g17150 PSSGPSVSEQIVTAVSASK      1 0.24919666       2
    ## 35117 XBmRNA31787|XBXL10_1g17150           SGVSLAALK      1 0.10753959       2
    ## 38995 XBmRNA31787|XBXL10_1g17150         TLAAAGYDVDK      1 0.06434523       2
    ## 2     XBmRNA73944|XBXL10_1g39228    AAAADGMEQMEMDESR      1 0.22193346       2
    ##        deltaBIC Simple_m0_heavy Simple_m0_ctrl    Simple_kd   Simple_kt
    ## 1     -2.638934       0.9684102      0.2657964 1.000000e-06 0.001916433
    ## 16860 -2.762997       1.0224265      0.2464643 3.680279e-05 0.002551913
    ## 31702 -2.759841       0.9945131      0.2709698 1.000000e-06 0.001160128
    ## 35117 -2.639098       1.0017751      0.2880129 1.878687e-05 0.002195705
    ## 38995 -2.654148       1.0132545      0.2190363 4.142279e-05 0.002269697
    ## 2     -2.605956       0.9826187      0.9324408 1.000000e-06 0.003292750
    ##       Baseline_m0_heavy Baseline_m0_ctrl  Baseline_kd  Baseline_kt Baseline_phi
    ## 1             0.9682724        0.2658564 1.000000e-06 0.0019136111         0.25
    ## 16860         1.0283180        0.3116009 3.561803e-05 0.0015291284         0.00
    ## 31702         0.9984860        0.3336643 1.000000e-06 0.0001700572         0.25
    ## 35117         1.0015894        0.2880477 1.854698e-05 0.0021941988         0.00
    ## 38995         1.0127828        0.2191004 5.480712e-05 0.0022685328         0.25
    ## 2             0.9841242        0.9321306 1.000000e-06 0.0032984323         0.25

``` r
max_fit_kd <- log(2)/(tail(avg_heavy_time, 1)*3)
max_fit_hl <- log(2)/max_fit_kd/60
cat("Max Half-life (hrs):\t", max_fit_hl) 
```

    ## Max Half-life (hrs):  146.95

``` r
T9_pep_fits <- T9_pep_fits[!is.na(T9_pep_fits$Simple_kd),]
norm_T9.df <- norm_T9.df[norm_T9.df$Peptide %in% T9_pep_fits$Peptide,]

T10_pep_fits <- T10_pep_fits[!is.na(T10_pep_fits$Simple_kd),]
norm_T10.df <- norm_T10.df[norm_T10.df$Peptide %in% T10_pep_fits$Peptide,]
```

``` r
pep_fit_comp <- T9_pep_fits
pep_fit_comp <- pep_fit_comp[!is.na(pep_fit_comp$deltaBIC),]

pep_fit_comp["Model"] <- rep("Unselected", nrow(pep_fit_comp))
pep_fit_comp[pep_fit_comp$deltaBIC < 0, "Model"] <- "Simple"
pep_fit_comp[pep_fit_comp$deltaBIC > 0, "Model"] <- "Two-pool"
pep_fit_comp[pep_fit_comp$deltaBIC > 10, "deltaBIC"] <- 10

table(pep_fit_comp$Model)
```

    ## 
    ##   Simple Two-pool 
    ##    44487      188

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

![](figures/Frog_Gast_SimpModel/Exploratory%20analysis%20not%20included%20in%20manuscript-1.png)<!-- -->

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


protein_med_merge <- rbind(cbind(norm_T9.df,
                                 data.frame(Exp=rep("T9", nrow(norm_T9.df)))),
                           cbind(norm_T10.df,
                                 data.frame(Exp=rep("T10", nrow(norm_T10.df)))))

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

    ##                   Protein_ID      T0_L        N2        N4        N6        N8
    ## 1 XBmRNA73944|XBXL10_1g39228 0.9902634 0.9757891 0.9875302 0.9704393 0.9930499
    ## 2 XBmRNA80099|XBXL10_1g42597 0.9605735 0.9598869 1.0242157 0.9051549 1.0958439
    ## 3 XBmRNA52689|XBXL10_1g28053 0.9423872 0.9154595 1.0454842 0.9754304 1.0210515
    ## 4 XBmRNA67792|XBXL10_1g35883 1.0168566 1.0028495 1.0743979 1.0203253 0.9825850
    ## 5  XBmRNA18430|XBXL10_1g9862 1.0016065 0.9433747 0.9690117 0.9556214 1.0710597
    ## 6 XBmRNA56795|XBXL10_1g30178 0.9923938 0.9735855 0.9987633 0.9920880 0.9884279
    ##         N10       N11      T0_H        O1        O2        O3        O4
    ## 1 1.0545842 1.0283438 1.0188206 1.0053499 1.0003034 0.9906056 1.0330509
    ## 2 0.9999542 1.0543710 0.9826704 1.0635280 0.9878824 0.9780995 0.9837637
    ## 3 0.9656949 1.1344924 0.9944406 0.9743798 1.0046140 1.0401318 1.0060431
    ## 4 0.9318439 0.9711416 1.0409522 1.0200658 1.0573793 1.0322114 1.0252243
    ## 5 1.0263005 1.0330256 1.0288425 0.9886572 1.0634323 1.0288698 0.8412957
    ## 6 1.0136293 1.0411122 1.0010623 0.9752478 1.0063893 1.0053325 0.9602227
    ##          O5        O6        O7        O8        O9       O10       O11
    ## 1 1.0046370 1.0055938 0.9878375 1.0031570 1.0019442 0.9836920 0.9650080
    ## 2 1.0492517 0.9862489 0.9983364 0.9784702 1.0343958 0.9403012 1.0170517
    ## 3 0.9721912 0.9673939 0.9819609 1.0403005 1.0622400 1.0364744 0.9198298
    ## 4 1.0341145 1.0269136 1.0065007 1.0070726 0.9582988 0.9172194 0.8740476
    ## 5 1.0046205 0.9749319 0.9629220 1.0400987 0.9950444 0.9855695 1.0857156
    ## 6 0.9802172 1.0053298 0.9670408 1.0219768 1.0032779 1.0431828 1.0307201

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
                             T9_kD = exp1_kD, T10_kD = exp2_kD,
                  out_df)
  
  return(out_df) }

Deg_fits.df <- lapply(protein_med_merge$Protein_ID,
                      function(x) { merge_pep_fits(x, T9_pep_fits, T10_pep_fits) }) %>%
  bind_rows() %>% arrange(desc(kD))

head(Deg_fits.df)
```

    ##                        Protein_ID          kD    HL_Hrs Estimate       T9_kD
    ## kd...1 XBmRNA48969|XBXL10_1g26175 0.016365203 0.7059156    point          NA
    ## kd...2  XBmRNA18031|XBXL10_1g9676 0.006690799 1.7266179    point 0.006690799
    ## kd...3 XBmRNA63717|XBXL10_1g33737 0.005468247 2.1126428    point 0.005468247
    ## kd...4 XBmRNA64557|XBXL10_1g34224 0.005440849 2.1232813    point 0.005440849
    ## kd...5 XBmRNA34126|XBXL10_1g18541 0.005127610 2.2529897 interval 0.005334910
    ## kd...6 XBmRNA23143|XBXL10_1g12383 0.004715197 2.4500468    point 0.004715197
    ##             T10_kD Fit_T0_H   Fit_O1   Fit_O2    Fit_O3    Fit_O4    Fit_O5
    ## kd...1 0.016365203 5.055340 3.285542 1.931567 0.7608722 0.1033316 0.0149820
    ## kd...2          NA 3.351341 2.773260 2.226041 1.5203195 0.6720945 0.3051748
    ## kd...3          NA 3.124033 2.683564 2.248694 1.6483977 0.8460290 0.4437677
    ## kd...4          NA 3.234032 2.746986 2.295513 1.6834205 0.8667924 0.4561306
    ## kd...5 0.004920311 3.186846 2.869813 2.441134 1.8243675 0.9759892 0.5329318
    ## kd...6          NA 3.059682 2.716938 2.330192 1.7812893 1.0020872 0.5744695
    ##              Fit_O6       Fit_O7       Fit_O8       Fit_O9      Fit_O10
    ## kd...1 1.130643e-06 2.204478e-10 9.020841e-13 3.555019e-16 8.298559e-20
    ## kd...2 6.297934e-03 1.915943e-04 2.023212e-05 8.206809e-07 2.687338e-08
    ## kd...3 1.861035e-02 1.071750e-03 1.706675e-04 1.243390e-05 7.604479e-07
    ## kd...4 1.943522e-02 1.135376e-03 1.824714e-04 1.346949e-05 8.353981e-07
    ## kd...5 2.723164e-02 1.873423e-03 3.345028e-04 2.868918e-05 2.088218e-06
    ## kd...6 3.728656e-02 3.181332e-03 6.524596e-04 6.818126e-05 6.126994e-06
    ##             Fit_O11
    ## kd...1 7.011287e-21
    ## kd...2 9.784756e-09
    ## kd...3 3.330199e-07
    ## kd...4 3.673592e-07
    ## kd...5 9.627534e-07
    ## kd...6 3.006298e-06

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

    ##                        Protein_ID          kD    HL_Hrs Estimate       T9_kD
    ## kd...1 XBmRNA48969|XBXL10_1g26175 0.016365203 0.7059156    point          NA
    ## kd...2  XBmRNA18031|XBXL10_1g9676 0.006690799 1.7266179    point 0.006690799
    ## kd...3 XBmRNA63717|XBXL10_1g33737 0.005468247 2.1126428    point 0.005468247
    ## kd...4 XBmRNA64557|XBXL10_1g34224 0.005440849 2.1232813    point 0.005440849
    ## kd...5 XBmRNA34126|XBXL10_1g18541 0.005127610 2.2529897 interval 0.005334910
    ## kd...6 XBmRNA23143|XBXL10_1g12383 0.004715197 2.4500468    point 0.004715197
    ##             T10_kD Fit_T0_H   Fit_O1   Fit_O2    Fit_O3    Fit_O4    Fit_O5
    ## kd...1 0.016365203 5.055340 3.285542 1.931567 0.7608722 0.1033316 0.0149820
    ## kd...2          NA 3.351341 2.773260 2.226041 1.5203195 0.6720945 0.3051748
    ## kd...3          NA 3.124033 2.683564 2.248694 1.6483977 0.8460290 0.4437677
    ## kd...4          NA 3.234032 2.746986 2.295513 1.6834205 0.8667924 0.4561306
    ## kd...5 0.004920311 3.186846 2.869813 2.441134 1.8243675 0.9759892 0.5329318
    ## kd...6          NA 3.059682 2.716938 2.330192 1.7812893 1.0020872 0.5744695
    ##              Fit_O6       Fit_O7       Fit_O8       Fit_O9      Fit_O10
    ## kd...1 1.130643e-06 2.204478e-10 9.020841e-13 3.555019e-16 8.298559e-20
    ## kd...2 6.297934e-03 1.915943e-04 2.023212e-05 8.206809e-07 2.687338e-08
    ## kd...3 1.861035e-02 1.071750e-03 1.706675e-04 1.243390e-05 7.604479e-07
    ## kd...4 1.943522e-02 1.135376e-03 1.824714e-04 1.346949e-05 8.353981e-07
    ## kd...5 2.723164e-02 1.873423e-03 3.345028e-04 2.868918e-05 2.088218e-06
    ## kd...6 3.728656e-02 3.181332e-03 6.524596e-04 6.818126e-05 6.126994e-06
    ##             Fit_O11 deltaBIC
    ## kd...1 7.011287e-21 30.18289
    ## kd...2 9.784756e-09 30.07042
    ## kd...3 3.330199e-07 36.28110
    ## kd...4 3.673592e-07 35.48077
    ## kd...5 9.627534e-07 58.05644
    ## kd...6 3.006298e-06 31.52686

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
final_xla_gast <- merge(human_map, Deg_fits.df, by="Protein_ID", all.y=TRUE)
final_xla_gast[is.na(final_xla_gast$XLA_Gene), "XLA_Gene"] <- "N/A"
final_xla_gast[is.na(final_xla_gast$Human_Gene), "Human_Gene"] <- "N/A"
final_xla_gast <- final_xla_gast[c(colnames(final_xla_gast)[1:8],
                                   "deltaBIC")]
final_xla_gast <- merge(final_xla_gast, protein_med_merge, by="Protein_ID")
final_xla_gast <- final_xla_gast %>% arrange(desc(kD))

head(final_xla_gast)
```

    ##                   Protein_ID      XLA_Gene Human_Gene          kD    HL_Hrs
    ## 1 XBmRNA48969|XBXL10_1g26175        odc1.S       ODC1 0.016365203 0.7059156
    ## 2  XBmRNA18031|XBXL10_1g9676       mkrn3.S      MKRN1 0.006690799 1.7266179
    ## 3 XBmRNA63717|XBXL10_1g33737  LOC108697796      OTOGL 0.005468247 2.1126428
    ## 4 XBmRNA64557|XBXL10_1g34224    pou5f3.3.L     POU3F2 0.005440849 2.1232813
    ## 5 XBmRNA34126|XBXL10_1g18541     LOC398156      PTTG1 0.005127610 2.2529897
    ## 6 XBmRNA23143|XBXL10_1g12383       fgfr4.L      FGFR2 0.004715197 2.4500468
    ##   Estimate       T9_kD      T10_kD deltaBIC      T0_L        N2        N4
    ## 1    point          NA 0.016365203 30.18289 0.6059934 1.0607010 1.0597994
    ## 2    point 0.006690799          NA 30.07042 1.7788198 1.4251402 1.8423530
    ## 3    point 0.005468247          NA 36.28110 0.7108096 0.5875566 0.7258997
    ## 4    point 0.005440849          NA 35.48077 3.1445249 2.0883277 1.0526202
    ## 5 interval 0.005334910 0.004920311 58.05644 1.7174742 1.8785722 1.6392702
    ## 6    point 0.004715197          NA 31.52686 0.7613206 0.8685615 0.7168250
    ##          N6        N8       N10       N11     T0_H       O1       O2        O3
    ## 1 1.0606275 0.5936988 1.1214435 1.4977365 4.837079 3.892699 1.604626 0.4356522
    ## 2 0.9779518 0.3061656 0.3386701 0.3308995 3.623325 2.312575 2.404421 1.4603548
    ## 3 1.0068287 1.3715976 1.3318800 1.2654278 3.022044 2.923162 2.143615 1.5629446
    ## 4 0.1852186 0.1454455 0.2019847 0.1818784 3.371389 2.377706 2.571615 1.6405517
    ## 5 0.5693240 0.3303836 0.4365402 0.4284358 3.222400 2.846712 2.513728 1.6351682
    ## 6 1.6847620 0.6625976 1.0728950 1.2330384 2.976464 2.502071 2.816976 1.6903271
    ##          O4        O5         O6          O7         O8         O9        O10
    ## 1 0.5031902 0.4761279 0.00000000 0.000000000 0.02265585 0.00000000 0.00000000
    ## 2 0.6604300 0.4915625 0.16823169 0.108968454 0.18779537 0.21790087 0.17368580
    ## 3 0.9227063 0.4191432 0.08576604 0.135173815 0.32385485 0.21384630 0.12129625
    ## 4 0.7671788 0.5582766 0.11847637 0.205411033 0.14970946 0.07423544 0.05804739
    ## 5 1.0852578 0.5494354 0.04253118 0.003653468 0.03180038 0.03617200 0.00000000
    ## 6 0.8928577 0.5474962 0.00000000 0.000000000 0.00000000 0.00000000 0.31930356
    ##          O11
    ## 1 0.22796976
    ## 2 0.19075001
    ## 3 0.12644722
    ## 4 0.10740307
    ## 5 0.03314245
    ## 6 0.25450494

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
                        "Frog_T9_light.min", "Frog_T9_all.min",
                        "avg_heavy_time",
                        "norm_T9.df", "norm_T10.df", "protein_med_merge",
                        "dyn_model_pep_fit", "true_shuffle",
                        "simple_normal_model",
                        "simple_heavy_model", "simple_pred_model"))
invisible(clusterEvalQ(cl, library(minpack.lm)))

NULL_BIC_values <- parApply(cl, Deg_fits.df, 1, function(row) {

  protein_id <- row["Protein_ID"]
  total_n <- length(heavy.channels)
  flat_predict <- rep(1, total_n)

  median_data <- protein_med_merge[protein_med_merge$Protein_ID == protein_id,]
  T9_check <- is.na(row["T9_kD"])
  T10_check <- is.na(row["T10_kD"])

  real_deltaBIC <- as.numeric(row["deltaBIC"])

  if (T9_check) {
    T10_peps <- norm_T10.df[norm_T10.df$Protein_ID == protein_id,]
    pep_data <- T10_peps[c("Peptide", "Pep_k1", "Pep_k2")] }
  else if (T10_check) {
    T9_peps <- norm_T9.df[norm_T9.df$Protein_ID == protein_id,]
    pep_data <- T9_peps[c("Peptide", "Pep_k1", "Pep_k2")] }
  else {
    T9_peps <- norm_T9.df[norm_T9.df$Protein_ID == protein_id,]
    T10_peps <- norm_T10.df[norm_T10.df$Protein_ID == protein_id,]

    pep_data <- rbind(T9_peps[c("Peptide", "Pep_k1", "Pep_k2")],
                      T10_peps[c("Peptide", "Pep_k1", "Pep_k2")])
    pep_data <- unique(pep_data) }

  NULL_BIC_deltas <- sapply(1:10, function(dummy) {

    #-- Fit ------------------------------------------

    # new_header <- heavy.channels
    new_header <- true_shuffle(heavy.channels)
    scrambled_median <- as.numeric(median_data[new_header])

    if (T9_check) {

      T10_fits <- apply(T10_peps, 1, function(x) {
        dyn_model_pep_fit(x, Frog_T9_light.min, Frog_T9_all.min, new_header) } )
      T10_fits <- do.call(rbind, T10_fits)
      merged_values <- apply(T10_fits, 2, median, na.rm = TRUE) }

    else if (T10_check) {

      T9_fits <- apply(T9_peps, 1, function(x) {
        dyn_model_pep_fit(x, Frog_T9_light.min, Frog_T9_all.min, new_header) } )
      T9_fits <- do.call(rbind, T9_fits)
      merged_values <- apply(T9_fits, 2, median, na.rm = TRUE) }

    else {

      T9_fits <- apply(T9_peps, 1, function(x) {
        dyn_model_pep_fit(x, Frog_T9_light.min, Frog_T9_all.min, new_header) } )
      T9_fits <- do.call(rbind, T9_fits)
      T9_fits <- apply(T9_fits, 2, median, na.rm = TRUE)

      T10_fits <- apply(T10_peps, 1, function(x) {
        dyn_model_pep_fit(x, Frog_T9_light.min, Frog_T9_all.min, new_header) } )
      T10_fits <- do.call(rbind, T10_fits)
      T10_fits <- apply(T10_fits, 2, median, na.rm = TRUE)

      merged_values <- (T9_fits + T10_fits)/2 }

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

    ##  [1] -13.994805 -13.994805 -13.283519 -10.087538 -13.287135 -13.994805
    ##  [7]  -7.334864  -8.619607  -9.140929 -13.994805 -10.261999  -9.870672
    ## [13] -12.779932  -4.450976 -13.896132 -12.862949 -12.744673 -11.965842
    ## [19]  -8.489175 -13.221589

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

    ## [1] "The 5% FDR Threshold is Delta BIC > -4.90000000000001"

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
# tiff("frog_gast_deltaBIC_hist_10variations.tiff", units="in",
#      width=6, height=4.5, res=300)

hist_deltaBIC
```

![](figures/Frog_Gast_SimpModel/BIC%20cutoff%20plot-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
final_BIC_threshold <- 2

final_xla_gast["mClass"] <- rep("Flat", nrow(final_xla_gast))
final_xla_gast[(final_xla_gast$deltaBIC > final_BIC_threshold),
               "mClass"] <- "Deg"
final_xla_gast[final_xla_gast$kD < max_fit_kd, "HL_Hrs"] <- max_fit_hl

head(final_xla_gast)
```

    ##                   Protein_ID      XLA_Gene Human_Gene          kD    HL_Hrs
    ## 1 XBmRNA48969|XBXL10_1g26175        odc1.S       ODC1 0.016365203 0.7059156
    ## 2  XBmRNA18031|XBXL10_1g9676       mkrn3.S      MKRN1 0.006690799 1.7266179
    ## 3 XBmRNA63717|XBXL10_1g33737  LOC108697796      OTOGL 0.005468247 2.1126428
    ## 4 XBmRNA64557|XBXL10_1g34224    pou5f3.3.L     POU3F2 0.005440849 2.1232813
    ## 5 XBmRNA34126|XBXL10_1g18541     LOC398156      PTTG1 0.005127610 2.2529897
    ## 6 XBmRNA23143|XBXL10_1g12383       fgfr4.L      FGFR2 0.004715197 2.4500468
    ##   Estimate       T9_kD      T10_kD deltaBIC      T0_L        N2        N4
    ## 1    point          NA 0.016365203 30.18289 0.6059934 1.0607010 1.0597994
    ## 2    point 0.006690799          NA 30.07042 1.7788198 1.4251402 1.8423530
    ## 3    point 0.005468247          NA 36.28110 0.7108096 0.5875566 0.7258997
    ## 4    point 0.005440849          NA 35.48077 3.1445249 2.0883277 1.0526202
    ## 5 interval 0.005334910 0.004920311 58.05644 1.7174742 1.8785722 1.6392702
    ## 6    point 0.004715197          NA 31.52686 0.7613206 0.8685615 0.7168250
    ##          N6        N8       N10       N11     T0_H       O1       O2        O3
    ## 1 1.0606275 0.5936988 1.1214435 1.4977365 4.837079 3.892699 1.604626 0.4356522
    ## 2 0.9779518 0.3061656 0.3386701 0.3308995 3.623325 2.312575 2.404421 1.4603548
    ## 3 1.0068287 1.3715976 1.3318800 1.2654278 3.022044 2.923162 2.143615 1.5629446
    ## 4 0.1852186 0.1454455 0.2019847 0.1818784 3.371389 2.377706 2.571615 1.6405517
    ## 5 0.5693240 0.3303836 0.4365402 0.4284358 3.222400 2.846712 2.513728 1.6351682
    ## 6 1.6847620 0.6625976 1.0728950 1.2330384 2.976464 2.502071 2.816976 1.6903271
    ##          O4        O5         O6          O7         O8         O9        O10
    ## 1 0.5031902 0.4761279 0.00000000 0.000000000 0.02265585 0.00000000 0.00000000
    ## 2 0.6604300 0.4915625 0.16823169 0.108968454 0.18779537 0.21790087 0.17368580
    ## 3 0.9227063 0.4191432 0.08576604 0.135173815 0.32385485 0.21384630 0.12129625
    ## 4 0.7671788 0.5582766 0.11847637 0.205411033 0.14970946 0.07423544 0.05804739
    ## 5 1.0852578 0.5494354 0.04253118 0.003653468 0.03180038 0.03617200 0.00000000
    ## 6 0.8928577 0.5474962 0.00000000 0.000000000 0.00000000 0.00000000 0.31930356
    ##          O11 mClass
    ## 1 0.22796976    Deg
    ## 2 0.19075001    Deg
    ## 3 0.12644722    Deg
    ## 4 0.10740307    Deg
    ## 5 0.03314245    Deg
    ## 6 0.25450494    Deg

``` r
ss_bp.all <- final_xla_gast

ss_bp.all <- ss_bp.all[!is.na(ss_bp.all$T9_kD),]
ss_bp.all <- ss_bp.all[!is.na(ss_bp.all$T10_kD),]

ss_bp.all["Avg_kD"] <- (ss_bp.all$T9_kD + ss_bp.all$T10_kD)/2

ss_bp.all["Log10_T9"] <- log10(ss_bp.all$T9_kD)
ss_bp.all["Log10_T10"] <- log10(ss_bp.all$T10_kD)
ss_bp.all[ss_bp.all$Log10_T9<(-8), "Log10_T9"] <- (-8)
ss_bp.all[ss_bp.all$Log10_T10<(-8), "Log10_T10"] <- (-8)

cor(ss_bp.all$Log10_T9, ss_bp.all$Log10_T10)^2
```

    ## [1] 0.4242768

``` r
cor(ss_bp.all[ss_bp.all$mClass=="Deg",]$Log10_T9,
    ss_bp.all[ss_bp.all$mClass=="Deg",]$Log10_T10)^2
```

    ## [1] 0.8168571

``` r
p1 <- ggplot() +
  geom_point(data=ss_bp.all[ss_bp.all$mClass == "Deg",],
             aes(x=Log10_T9, y=Log10_T10),
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

![](figures/Frog_Gast_SimpModel/Degradation%20rate%20biplot%20for%20degrading%20proteins-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
peptide_count <- data.frame(table(rbind(norm_T9.df[1:2],
                                        norm_T10.df[1:2])$Protein_ID))
colnames(peptide_count) <- c("Protein_ID", "Number_Pep_Fits")

final_xla_gast <- merge(final_xla_gast, peptide_count, by="Protein_ID") %>%
  arrange(desc(kD))

head(final_xla_gast)
```

    ##                   Protein_ID      XLA_Gene Human_Gene          kD    HL_Hrs
    ## 1 XBmRNA48969|XBXL10_1g26175        odc1.S       ODC1 0.016365203 0.7059156
    ## 2  XBmRNA18031|XBXL10_1g9676       mkrn3.S      MKRN1 0.006690799 1.7266179
    ## 3 XBmRNA63717|XBXL10_1g33737  LOC108697796      OTOGL 0.005468247 2.1126428
    ## 4 XBmRNA64557|XBXL10_1g34224    pou5f3.3.L     POU3F2 0.005440849 2.1232813
    ## 5 XBmRNA34126|XBXL10_1g18541     LOC398156      PTTG1 0.005127610 2.2529897
    ## 6 XBmRNA23143|XBXL10_1g12383       fgfr4.L      FGFR2 0.004715197 2.4500468
    ##   Estimate       T9_kD      T10_kD deltaBIC      T0_L        N2        N4
    ## 1    point          NA 0.016365203 30.18289 0.6059934 1.0607010 1.0597994
    ## 2    point 0.006690799          NA 30.07042 1.7788198 1.4251402 1.8423530
    ## 3    point 0.005468247          NA 36.28110 0.7108096 0.5875566 0.7258997
    ## 4    point 0.005440849          NA 35.48077 3.1445249 2.0883277 1.0526202
    ## 5 interval 0.005334910 0.004920311 58.05644 1.7174742 1.8785722 1.6392702
    ## 6    point 0.004715197          NA 31.52686 0.7613206 0.8685615 0.7168250
    ##          N6        N8       N10       N11     T0_H       O1       O2        O3
    ## 1 1.0606275 0.5936988 1.1214435 1.4977365 4.837079 3.892699 1.604626 0.4356522
    ## 2 0.9779518 0.3061656 0.3386701 0.3308995 3.623325 2.312575 2.404421 1.4603548
    ## 3 1.0068287 1.3715976 1.3318800 1.2654278 3.022044 2.923162 2.143615 1.5629446
    ## 4 0.1852186 0.1454455 0.2019847 0.1818784 3.371389 2.377706 2.571615 1.6405517
    ## 5 0.5693240 0.3303836 0.4365402 0.4284358 3.222400 2.846712 2.513728 1.6351682
    ## 6 1.6847620 0.6625976 1.0728950 1.2330384 2.976464 2.502071 2.816976 1.6903271
    ##          O4        O5         O6          O7         O8         O9        O10
    ## 1 0.5031902 0.4761279 0.00000000 0.000000000 0.02265585 0.00000000 0.00000000
    ## 2 0.6604300 0.4915625 0.16823169 0.108968454 0.18779537 0.21790087 0.17368580
    ## 3 0.9227063 0.4191432 0.08576604 0.135173815 0.32385485 0.21384630 0.12129625
    ## 4 0.7671788 0.5582766 0.11847637 0.205411033 0.14970946 0.07423544 0.05804739
    ## 5 1.0852578 0.5494354 0.04253118 0.003653468 0.03180038 0.03617200 0.00000000
    ## 6 0.8928577 0.5474962 0.00000000 0.000000000 0.00000000 0.00000000 0.31930356
    ##          O11 mClass Number_Pep_Fits
    ## 1 0.22796976    Deg               1
    ## 2 0.19075001    Deg               1
    ## 3 0.12644722    Deg               1
    ## 4 0.10740307    Deg               1
    ## 5 0.03314245    Deg               3
    ## 6 0.25450494    Deg               1

``` r
supp_df <- merge(protein_annotation, final_xla_gast,
                 by="Protein_ID", all.y=TRUE) %>%
  arrange(desc(kD)) %>%
  dplyr::select(Protein_ID, XLA_Gene, Human_Gene,
                Description, everything())

supp_df["Is_Degrading"] <- rep("Unknown", nrow(supp_df))
supp_df[supp_df$mClass == "Deg", "Is_Degrading"] <- "Yes"
supp_df[supp_df$mClass == "Flat", "Is_Degrading"] <- "No"

supp_df[c("mClass")] <- NULL
colnames(supp_df) <-
  c("Protein ID", "XenBase Annotation", "Human Protein", "Description",
    "Fitted kD", "Reported Half-life (Hrs)", "Estimate Confidence",
    "Replicate 1 kD", "Replicate 2 kD", "BIC against NULL",
    paste0("Control: ", round((avg_light_time+1472)/60, 1), " hpf"),
    paste0("18O: ", round((avg_heavy_time+1472)/60, 1), " hpf"),
    "Peptides Quantified", "Is Degrading?")

supp_df[(supp_df$`Is Degrading?`=="No") &
          (supp_df$`Reported Half-life (Hrs)` < max_fit_hl),
        "Reported Half-life (Hrs)"] <- NA
supp_df <- rbind(supp_df[!is.na(supp_df$`Reported Half-life (Hrs)`),],
                 supp_df[is.na(supp_df$`Reported Half-life (Hrs)`),])

dir.create("Files/Supp_Tables", showWarnings = FALSE, recursive = TRUE)
write.csv(supp_df, "Files/Supp_Tables/Frog_Gastrulation_Supplementary_Table.csv",
          row.names = FALSE)
```

``` r
write.csv(final_xla_gast,
          "Files/Fits/Dyn-Model/XLA_O18_T9-10_FinalFits_ALL.csv",
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
