Frog AA Incorporation
================
Edward Cruz
2025-04-25

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
library(splines)
library(stringr)
library(ggplot2)
library(patchwork)
library(data.table)
```

    ## 
    ## Attaching package: 'data.table'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     between, first, last

``` r
library(minpack.lm)
library(parallel)

num_cores <- detectCores() - 1  # Save one core for system stability
knitr::opts_chunk$set(fig.path = "figures/Frog_AA_Incorporation_MBT-T2/")
```

``` r
aa.profiles <- read.csv("Files/Reference/AA_Profiles.csv")

head(aa.profiles)
```

    ##   AA       Name   Formula     MW Nitrogen Carbon Oxygen Hydrogen Sulfur
    ## 1  A    alanine   C3H7NO2  89.10        1      3      2        7      0
    ## 2  R   arginine C6H14N4O2 174.20        4      6      2       14      0
    ## 3  N asparagine  C4H8N2O3 132.12        2      4      3        8      0
    ## 4  D  aspartate   C4H7NO4 133.11        1      4      4        7      0
    ## 5  C   cysteine  C3H7NO2S 121.16        1      3      2        7      1
    ## 6  E  glutamate   C5H9NO4 147.13        1      5      4        9      0

``` r
#Time in minutes of timepoints collected, spline time is coarser
#Note: These are the exact timepoint matches for proteomics samples
AA_T2.time <- c(0, 30, 60, 120, 180, 300, 485,
                725, 1440, 1800, 2886)
```

``` r
merged_frog_table <- read.csv("Data/XLA_AAA/XLA_O18-T1-T2_results-merge_AA.csv")
merged_frog_table <- merged_frog_table[!grepl("blank", merged_frog_table$Sample),]

early_exp_id <- "col009b"
mbt_exp_id <- "col009b2"

MBT_AA_table <- merged_frog_table[grepl(paste(mbt_exp_id, "_", sep=""),
                                        merged_frog_table$Sample),]

cat("Missing AAs: ", 
    aa.profiles$Name[!aa.profiles$Name %in% unique(MBT_AA_table$compound)],
    "\n")
```

    ## Missing AAs:  cysteine leucine

``` r
#Removing inaccurate timepoint
AA_T2.time <- c(0, 30, 60, 120, 180, 300, 485,
                725, 1440, 2886)

MBT_AA_table <- MBT_AA_table[!MBT_AA_table$Sample == "col009b2_T9",]
```

``` r
#'*Looping through unique amino acids and reformatting table*
MBT_AA_matrix <- sapply(unique(MBT_AA_table$compound), function(x){
  aa.data <- MBT_AA_table[MBT_AA_table$compound==x,] #selecting AA

  aa_wide <- aa.data %>% #Reformating columns long-wise
    pivot_wider(
      id_cols = c(compound, Sample),
      names_from = isotopeLabel,
      values_from = fraction.all ) %>%
    arrange(compound, Sample)  

  #Binding total number of O18-labels
  aa_wide <- cbind(data.frame(Total=rep(sum(grepl("O18", colnames(aa_wide))),
                                        nrow(aa_wide))),
                   aa_wide)
    
  return(aa_wide) })

MBT_AA_matrix <- bind_rows(MBT_AA_matrix)

#'*Calculating median values for missing AAs*
MBT_med_values <- lapply(paste("T", c(seq(0,8,1),10), sep=""), function(time){
  #Calculating median for timepoints collected
  t_data <- MBT_AA_matrix[MBT_AA_matrix$Sample == paste("col009b2_", time, sep=""),]

    
  sub_2 <- t_data[t_data$compound %in% c("arginine", "cysteine",
                                         "glycine", "valine", "histidine",
                                         "isoleucine", "leucine", "lysine",
                                         "methionine", "phenylalanine",
                                         "proline", "tryptophan"),]
  med_2 <- apply(sub_2[c("C12 PARENT", "O18-label-1",
                         "O18-label-2")], 2, median)
  med_2 <- med_2 / sum(med_2) #Necessary to remake row sum = 1

  
  sub_3 <- t_data[t_data$compound %in% c("glutamine"),]  
  med_3 <- apply(sub_3[c("C12 PARENT", "O18-label-1",
                         "O18-label-2", "O18-label-3")], 2, median)
  med_3 <- med_3 / sum(med_3) #Necessary to remake row sum = 1

    
  sub_4 <- t_data[t_data$compound %in% c("glutamate", "aspartate"),]
  med_4 <- apply(sub_4[c("C12 PARENT", "O18-label-1",
                         "O18-label-2", "O18-label-3",
                         "O18-label-4")], 2, median)
  med_4 <- med_4 / sum(med_4) #Necessary to remake row sum = 1
  

  row_out <- data.frame(Name=rep("median", 3),
                        AA=c("X", "U", "Z"),
                        Total=c(2,3,4),
                        Sample=paste("col009b2_", time, sep=""))
  
  return(cbind(row_out, bind_rows(med_2, med_3, med_4))) })

#Merging and appending to table of all amino acids
MBT_med_values <- do.call(rbind, MBT_med_values)
MBT_AA_matrix <- merge(aa.profiles[c("Name", "AA")],
                       MBT_AA_matrix, by.x="Name",
                       by.y="compound")

MBT_AA_matrix <- bind_rows(MBT_AA_matrix, MBT_med_values)
MBT_AA_matrix$Name <- sub("^(\\w)", "\\U\\1", MBT_AA_matrix$Name, perl = TRUE)

head(MBT_AA_matrix)
```

    ##      Name AA Total       Sample C12 PARENT O18-label-1 O18-label-2 O18-label-3
    ## 1 Alanine  A     2  col009b2_T0  0.9941659  0.00583406  0.00000000          NA
    ## 2 Alanine  A     2  col009b2_T1  0.9034068  0.07065839  0.02593479          NA
    ## 3 Alanine  A     2 col009b2_T10  0.2810308  0.49171793  0.22725124          NA
    ## 4 Alanine  A     2  col009b2_T2  0.8849653  0.07462412  0.04041059          NA
    ## 5 Alanine  A     2  col009b2_T3  0.8065851  0.10372388  0.08969100          NA
    ## 6 Alanine  A     2  col009b2_T4  0.7314560  0.15039951  0.11814450          NA
    ##   O18-label-4
    ## 1          NA
    ## 2          NA
    ## 3          NA
    ## 4          NA
    ## 5          NA
    ## 6          NA

``` r
write.csv(MBT_AA_matrix, "Data/XLA_O18/Isotopic_Envelopes/Frog_MBT_AA_Matrix.csv",
          row.names = FALSE)
```

``` r
#Frog Essential AAs
measured_AAs <- unique(MBT_AA_matrix$AA)

no_Oxy_AA <- c("A", "R", "C", "G", "V",
               "H", "I", "L", "K",
               "M", "F", "P", "W",
               "T", "Y", "X") 
no_Oxy_AA <- no_Oxy_AA[no_Oxy_AA %in% setdiff(measured_AAs, c("A"))]

one_Oxy_AA <- c("Q", "S", "N", "U") 
one_Oxy_AA <- one_Oxy_AA[one_Oxy_AA %in% setdiff(measured_AAs, c("N"))]

two_Oxy_AA <- c("D", "E", "Z")
two_Oxy_AA <- two_Oxy_AA[two_Oxy_AA %in% measured_AAs]
```

``` r
no_Oxy.splines <- do.call(rbind, lapply(no_Oxy_AA, function(AA) {

  # Example setup
  interpolated_times <- seq(0, 2886, by = 1)
  sub.data <- MBT_AA_matrix[MBT_AA_matrix$AA == AA,]

  sub.data["time"] <- gsub("col009b2_T", "", sub.data$Sample)  
    
  sub.data["Sample"] <- gsub("col009b2_T", "", sub.data$Sample) %>% as.numeric
  sub.data <- sub.data %>% arrange(Sample)

  if (AA=="Y") { knot_positions <- c(50, 200) }   
  else if (AA %in% c("M")) { knot_positions <- c(20, 200) }     
  else if (AA %in% c("P", "H")) { knot_positions <- c(50, 500) }     
  else if (AA %in% c("A")) { knot_positions <- c(50, 00) }       
  else { knot_positions <- c(120, 300)   }
  

  # Fit spline models
  c12_fit <- lm(`C12 PARENT` ~ bs(AA_T2.time, degree = 3,
                                  knots = knot_positions), data = sub.data)
  o18_1_fit <- lm(`O18-label-1` ~ bs(AA_T2.time, degree = 3,
                                     knots = knot_positions), data = sub.data)
  o18_2_fit <- lm(`O18-label-2` ~ bs(AA_T2.time, degree = 3,
                                     knots = knot_positions), data = sub.data)

  # Predict spline values
  c12_pred <- predict(c12_fit, newdata = data.frame(AA_T2.time = interpolated_times))
  o18_1_pred <- predict(o18_1_fit, newdata = data.frame(AA_T2.time = interpolated_times))
  o18_2_pred <- predict(o18_2_fit, newdata = data.frame(AA_T2.time = interpolated_times))

  # Combine into one dataframe for ggplot
  spline_df <- data.frame(AA=rep(AA, length(interpolated_times)),
                          Labels=rep(2, length(interpolated_times)),
                          Time=rep(interpolated_times, 3),
                          Type=c(rep("C12", length(interpolated_times)),
                                 rep("O18_1", length(interpolated_times)),
                                 rep("O18_2", length(interpolated_times))),
                          Value=c(c12_pred, o18_1_pred, o18_2_pred))
  spline_df["AA"] <- AA
  spline_df[spline_df$Value > 1, "Value"] <- 1
  spline_df[spline_df$Value < 0, "Value"] <- 0

  # x <- ggplot() +
  #   geom_point(data=data.frame(x=AA_T2.time,
  #                              y=sub.data$`C12 PARENT`),
  #              aes(x=x,y=y), size=3) +
  #   geom_line(data=spline_df[spline_df$Type=="C12",],
  #              aes(x=Time,y=Value), size=1) +
  #   theme_bw() +
  #   coord_cartesian(ylim=c(0,1))
  # print(x)
  # invokeRestart("abort")
      
  return(spline_df)
    
}))
```

``` r
one_Oxy.splines <- do.call(rbind, lapply(one_Oxy_AA,
                                         function(AA) {

  # Example setup
  interpolated_times <- seq(0, 2886, by = 1)
  sub.data <- MBT_AA_matrix[MBT_AA_matrix$AA == AA,]
  sub.data["Sample"] <- gsub("col009b2_T", "", sub.data$Sample) %>% as.numeric
  sub.data <- sub.data %>% arrange(Sample)


  if (AA=="U") { knot_positions <- c(120,500) }
  
  else if (AA=="S") { knot_positions <- c(60,300) }
  else { knot_positions <- c(180,500) }
    
  # Fit spline models
  c12_fit <- lm(`C12 PARENT` ~ bs(AA_T2.time, degree = 3,
                                  knots = knot_positions), data = sub.data)
  o18_1_fit <- lm(`O18-label-1` ~ bs(AA_T2.time, degree = 3,
                                     knots = knot_positions), data = sub.data)
  o18_2_fit <- lm(`O18-label-2` ~ bs(AA_T2.time, degree = 3,
                                     knots = knot_positions), data = sub.data)
  o18_3_fit <- lm(`O18-label-3` ~ bs(AA_T2.time, degree = 3,
                                     knots = knot_positions), data = sub.data)
  
  # Predict spline values
  c12_pred <- predict(c12_fit, newdata = data.frame(AA_T2.time = interpolated_times))
  o18_1_pred <- predict(o18_1_fit, newdata = data.frame(AA_T2.time = interpolated_times))
  o18_2_pred <- predict(o18_2_fit, newdata = data.frame(AA_T2.time = interpolated_times))
  o18_3_pred <- predict(o18_3_fit, newdata = data.frame(AA_T2.time = interpolated_times))

  # Combine into one dataframe for ggplot
  spline_df <- data.frame(AA=rep(AA, length(interpolated_times)),
                          Labels=rep(3, length(interpolated_times)),
                          Time=rep(interpolated_times, 4),
                          Type=c(rep("C12", length(interpolated_times)),
                                 rep("O18_1", length(interpolated_times)),
                                 rep("O18_2", length(interpolated_times)),
                                 rep("O18_3", length(interpolated_times))),
                          Value=c(c12_pred, o18_1_pred, o18_2_pred, o18_3_pred))
  spline_df["AA"] <- AA
  spline_df[spline_df$Value > 1, "Value"] <- 1
  spline_df[spline_df$Value < 0, "Value"] <- 0

  # x <- ggplot() +
  #   geom_point(data=data.frame(x=AA_T2.time,
  #                              y=sub.data$`C12 PARENT`),
  #              aes(x=x,y=y), size=3) +
  #   geom_line(data=spline_df[spline_df$Type=="C12",],
  #              aes(x=Time,y=Value), size=1) +
  #   theme_bw() +
  #   coord_cartesian(ylim=c(0,1))
  # print(x)
  
  return(spline_df)

}))
```

``` r
two_Oxy.splines <- do.call(rbind, lapply(two_Oxy_AA,
                                         function(AA) {

  # Example setup
  interpolated_times <- seq(0, 2886, by = 1)
  sub.data <- MBT_AA_matrix[MBT_AA_matrix$AA == AA,]
  sub.data["Sample"] <- gsub("col009b2_T", "", sub.data$Sample) %>% as.numeric
  sub.data <- sub.data %>% arrange(Sample)

  knot_positions <- c(120,200)    

  # Fit spline models
  c12_fit <- lm(`C12 PARENT` ~ bs(AA_T2.time, degree = 3,
                                  knots = knot_positions), data = sub.data)
  o18_1_fit <- lm(`O18-label-1` ~ bs(AA_T2.time, degree = 3,
                                     knots = knot_positions), data = sub.data)
  o18_2_fit <- lm(`O18-label-2` ~ bs(AA_T2.time, degree = 3,
                                     knots = knot_positions), data = sub.data)
  o18_3_fit <- lm(`O18-label-3` ~ bs(AA_T2.time, degree = 3,
                                     knots = knot_positions), data = sub.data)
  o18_4_fit <- lm(`O18-label-4` ~ bs(AA_T2.time, degree = 3,
                                     knots = knot_positions), data = sub.data)
    
  # Predict spline values
  c12_pred <- predict(c12_fit, newdata = data.frame(AA_T2.time = interpolated_times))
  o18_1_pred <- predict(o18_1_fit, newdata = data.frame(AA_T2.time = interpolated_times))
  o18_2_pred <- predict(o18_2_fit, newdata = data.frame(AA_T2.time = interpolated_times))
  o18_3_pred <- predict(o18_3_fit, newdata = data.frame(AA_T2.time = interpolated_times))
  o18_4_pred <- predict(o18_4_fit, newdata = data.frame(AA_T2.time = interpolated_times))
  
  # Combine into one dataframe for ggplot
  spline_df <- data.frame(AA=rep(AA, length(interpolated_times)*5),
                          Labels=rep(4, length(interpolated_times)*5),
                          Time=rep(interpolated_times, 5),
                          Type=c(rep("C12", length(interpolated_times)),
                                 rep("O18_1", length(interpolated_times)),
                                 rep("O18_2", length(interpolated_times)),
                                 rep("O18_3", length(interpolated_times)),
                                 rep("O18_4", length(interpolated_times))),
                          Value=c(c12_pred, o18_1_pred, o18_2_pred, o18_3_pred,
                                  o18_4_pred))
  spline_df["AA"] <- AA
  spline_df[spline_df$Value > 1, "Value"] <- 1
  spline_df[spline_df$Value < 0, "Value"] <- 0

  # x <- ggplot() +
  #   geom_point(data=data.frame(x=AA_T2.time,
  #                              y=sub.data$`C12 PARENT`),
  #              aes(x=x,y=y), size=3) +
  #   geom_line(data=spline_df[spline_df$Type=="C12",],
  #              aes(x=Time,y=Value), size=1) +
  #   theme_bw() +
  #   coord_cartesian(ylim=c(0,1))
  # print(x)
  
  return(spline_df)
    
}))
```

``` r
all_splines.merge <- rbind(no_Oxy.splines, one_Oxy.splines, two_Oxy.splines)
all_splines.merge <- do.call(rbind,
              lapply(unique(all_splines.merge$AA), function(AA) {

  #Looping through each time to verify rowsum=1 for each AA
  c_df <- lapply(seq(0, 2886, by = 1), function(t) {
    sub_t <- all_splines.merge[all_splines.merge$AA==AA &
                              all_splines.merge$Time==t,]
    sub_t["Value"] <- sub_t$Value / sum(sub_t$Value)
    
    return(sub_t) })
  
  c_df <- do.call(rbind, c_df)

  return(c_df) }))

head(all_splines.merge)
```

    ##      AA Labels Time  Type       Value
    ## 1     R      2    0   C12 0.996867571
    ## 2888  R      2    0 O18_1 0.003132429
    ## 5775  R      2    0 O18_2 0.000000000
    ## 2     R      2    1   C12 0.984802319
    ## 2889  R      2    1 O18_1 0.015197681
    ## 5776  R      2    1 O18_2 0.000000000

``` r
all_spline.p <- lapply(no_Oxy_AA, function(ind_AA){
  
  sub.data <- MBT_AA_matrix[MBT_AA_matrix$AA == ind_AA,]
  sub.data["Sample"] <- gsub("col009b2_T", "", sub.data$Sample) %>% 
    as.numeric
  sub.data <- sub.data %>% arrange(Sample)
  sub.data["Time"] <- (AA_T2.time+1470)/60
  
  spline.data <- no_Oxy.splines[no_Oxy.splines$AA==ind_AA,]
  spline.data["Time"] <- (spline.data$Time+1470)/60
  
  AA_spline.p <- ggplot() +
    # geom_point(data=sub.data, aes(x=Time, y=`O18-label-1`),
    #            color="#009E73", size=3) +
    geom_line(data=spline.data[spline.data$Type=="O18_1",],
              aes(x=Time, y=`Value`),
              color="#009E73", size=1) +
    # geom_point(data=sub.data, aes(x=Time, y=`O18-label-2`),
    #            color="#B45100", size=3) +
    geom_line(data=spline.data[spline.data$Type=="O18_2",],
              aes(x=Time, y=`Value`),
              color="#B45100", size=1) +    
    geom_point(data=sub.data, aes(x=Time, y=`C12 PARENT`),
               color="black", size=3) +
    geom_line(data=spline.data[spline.data$Type=="C12",],
              aes(x=Time, y=`Value`),
              color="black", size=3, alpha=0.4) +
    labs(x="Hours after Labeling", y="Relative Abundance",
         title=unique(sub.data$Name)) +
    theme_bw() +
    theme(axis.text=element_text(size=21,colour="black"),
          # plot.title = element_text(size=24, hjust=0.5, colour="#5E4FA2"),
          plot.title = element_text(size=24, hjust=0.5, colour="black"),          
          # axis.title = element_text(size=22),
          axis.title = element_blank(),                    
          panel.border = element_rect(linewidth=2),
          panel.grid.major = element_line(color = "grey70", linewidth = 0.25),
          panel.grid.minor = element_line(color = "grey70", linewidth = 0.25),
          aspect.ratio = 1,
          legend.position = "none") +
    coord_cartesian(xlim=c(24, 50), ylim=c(0,1)) +
    scale_x_continuous(breaks=seq(20,50,5))
  
  return(AA_spline.p) })
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once per session.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
all_spline.p <- c(all_spline.p,
                  lapply(one_Oxy_AA, function(ind_AA){

  sub.data <- MBT_AA_matrix[MBT_AA_matrix$AA == ind_AA,]
  sub.data["Sample"] <- gsub("col009b2_T", "", sub.data$Sample) %>% 
    as.numeric
  sub.data <- sub.data %>% arrange(Sample)
  sub.data["Time"] <- (AA_T2.time+1470)/60

  spline.data <- one_Oxy.splines[one_Oxy.splines$AA==ind_AA,]
  spline.data["Time"] <- (spline.data$Time+1470)/60

  AA_spline.p <- ggplot() +
    # geom_point(data=sub.data, aes(x=Time, y=`O18-label-1`),
    #            color="#009E73", size=3) +
    geom_line(data=spline.data[spline.data$Type=="O18_1",],
              aes(x=Time, y=`Value`),
              color="#009E73", size=1) +
    # geom_point(data=sub.data, aes(x=Time, y=`O18-label-2`),
    #            color="#D55E00", size=3) +
    geom_line(data=spline.data[spline.data$Type=="O18_2",],
              aes(x=Time, y=`Value`),
              color="#D55E00", size=1) +
    # geom_point(data=sub.data, aes(x=Time, y=`O18-label-3`),
    #            color="#0072B2", size=3) +
    geom_line(data=spline.data[spline.data$Type=="O18_3",],
              aes(x=Time, y=`Value`),
              color="#0072B2", size=1) +
    geom_point(data=sub.data, aes(x=Time, y=`C12 PARENT`),
               color="black", size=3) +
    geom_line(data=spline.data[spline.data$Type=="C12",],
              aes(x=Time, y=`Value`),
              color="black", size=3, alpha=0.4) +
    labs(x="Hours after Labeling", y="Relative Abundance",
         title=unique(sub.data$Name)) +
    theme_bw() +
    theme(axis.text=element_text(size=21,colour="black"),
          # plot.title = element_text(size=24, hjust=0.5, colour="#E7298A"),
          plot.title = element_text(size=24, hjust=0.5, colour="black"),          
          # axis.title = element_text(size=22),
          axis.title = element_blank(),                    
          panel.border = element_rect(linewidth=2),
          panel.grid.major = element_line(color = "grey70", linewidth = 0.25),
          panel.grid.minor = element_line(color = "grey70", linewidth = 0.25),
          aspect.ratio = 1,
          legend.position = "none") +
    coord_cartesian(xlim=c(24, 50), ylim=c(0,1)) +
    scale_x_continuous(breaks=seq(20,50,5))

  return(AA_spline.p) }) )


all_spline.p <- c(all_spline.p,
                  lapply(two_Oxy_AA, function(ind_AA){

  sub.data <- MBT_AA_matrix[MBT_AA_matrix$AA == ind_AA,]
  sub.data["Sample"] <- gsub("col009b2_T", "", sub.data$Sample) %>% 
    as.numeric
  sub.data <- sub.data %>% arrange(Sample)
  sub.data["Time"] <- (AA_T2.time+1470)/60

  spline.data <- two_Oxy.splines[two_Oxy.splines$AA==ind_AA,]
  spline.data["Time"] <- (spline.data$Time+1470)/60

  AA_spline.p <- ggplot() +
    # geom_point(data=sub.data, aes(x=Time, y=`O18-label-1`),
    #            color="#009E73", size=3) +
    geom_line(data=spline.data[spline.data$Type=="O18_1",],
              aes(x=Time, y=`Value`),
              color="#009E73", size=1) +
    # geom_point(data=sub.data, aes(x=Time, y=`O18-label-2`),
    #            color="#D55E00", size=3) +
    geom_line(data=spline.data[spline.data$Type=="O18_2",],
              aes(x=Time, y=`Value`),
              color="#D55E00", size=1) +
    # geom_point(data=sub.data, aes(x=Time, y=`O18-label-3`),
    #            color="#0072B2", size=3) +
    geom_line(data=spline.data[spline.data$Type=="O18_3",],
              aes(x=Time, y=`Value`),
              color="#0072B2", size=1) +
    # geom_point(data=sub.data, aes(x=Time, y=`O18-label-4`),
    #            color="#CC79A7", size=3) +
    geom_line(data=spline.data[spline.data$Type=="O18_4",],
              aes(x=Time, y=`Value`),
              color="#CC79A7", size=1) +
    geom_point(data=sub.data, aes(x=Time, y=`C12 PARENT`),
               color="black", size=3) +
    geom_line(data=spline.data[spline.data$Type=="C12",],
              aes(x=Time, y=`Value`),
              color="black", size=3, alpha=0.4) +
    labs(x="Hours after Labeling", y="Relative Abundance",
         title=unique(sub.data$Name)) +
    theme_bw() +
    theme(axis.text=element_text(size=21,colour="black"),
          # plot.title = element_text(size=24, hjust=0.5, colour="#F0E442"),
          plot.title = element_text(size=24, hjust=0.5, colour="black"),          
          # axis.title = element_text(size=22),
          axis.title = element_blank(),          
          panel.border = element_rect(linewidth=2),
          panel.grid.major = element_line(color = "grey70", linewidth = 0.25),
          panel.grid.minor = element_line(color = "grey70", linewidth = 0.25),
          aspect.ratio = 1,
          legend.position = "none") +
    coord_cartesian(xlim=c(24, 50), ylim=c(0,1)) +
    scale_x_continuous(breaks=seq(20,50,5))

  return(AA_spline.p) }) )

names(all_spline.p) <- c(no_Oxy_AA, one_Oxy_AA, two_Oxy_AA)
```

``` r
# #Uncomment for new image!
# tiff("graphs.tiff", units="in",
#      width=9, height=5, res=300)

all_spline.p[["X"]] | all_spline.p[["U"]] | all_spline.p[["Z"]]
```

![](figures/Frog_AA_Incorporation_MBT-T2/Median%20plots%20of%20AAs-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=14, height=13, res=300)

(all_spline.p[["R"]] | all_spline.p[["G"]] | all_spline.p[["V"]] | all_spline.p[["H"]]) /
  (all_spline.p[["I"]] | all_spline.p[["K"]] | all_spline.p[["M"]] | all_spline.p[["F"]]) /
  (all_spline.p[["P"]] | all_spline.p[["W"]] | all_spline.p[["T"]] | all_spline.p[["Y"]]) /
  (all_spline.p[["S"]] | all_spline.p[["Q"]] | all_spline.p[["D"]] | all_spline.p[["E"]])
```

![](figures/Frog_AA_Incorporation_MBT-T2/Supp%20Figure%20of%20all%20amino%20acids-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
AA_matrix_fit <- data.frame(Time=seq(0, 2886, by = 1))

AA_matrix_fit <- cbind(AA_matrix_fit,
                       do.call(cbind, lapply(unique(all_splines.merge$AA),
                                             function(AA){
                                               
  sub_data <- all_splines.merge[all_splines.merge$AA==AA&
                                  all_splines.merge$Type=="C12",]
  AA_labels <- unique(sub_data$Labels)

  if (AA %in% no_Oxy_AA) { p_ao <- sub_data$Value^(1/2)  }
  else if (AA %in% one_Oxy_AA) { p_ao <- sub_data$Value^(2/3)  }
  else if (AA %in% two_Oxy_AA) { p_ao <- sub_data$Value^(3/4)  }  
  
  out_col <- data.frame(AA_Name=p_ao)
  colnames(out_col) <- AA
  
  return(out_col) })))

matrix_header <- colnames(AA_matrix_fit)[2:length(colnames(AA_matrix_fit))]
AA_matrix_fit <- AA_matrix_fit %>% as.matrix

head(AA_matrix_fit)
```

    ##      Time         R         G         V         H         I         K         M
    ## [1,]    0 0.9984326 1.0000000 1.0000000 0.9968675 1.0000000 1.0000000 1.0000000
    ## [2,]    1 0.9923721 0.9962097 0.9959525 0.9934185 0.9943906 0.9963840 0.9623039
    ## [3,]    2 0.9857971 0.9893887 0.9903153 0.9899853 0.9888188 0.9881214 0.9267310
    ## [4,]    3 0.9779393 0.9826164 0.9844516 0.9865901 0.9826666 0.9799149 0.8933373
    ## [5,]    4 0.9701377 0.9758931 0.9785521 0.9832503 0.9764724 0.9717650 0.8621711
    ## [6,]    5 0.9623929 0.9692191 0.9726959 0.9799652 0.9703278 0.9636719 0.8332706
    ##              F         P         W         T         Y         X         Q
    ## [1,] 0.9998505 0.9989682 0.9910346 0.9985729 0.9962528 0.9971213 0.9979930
    ## [2,] 0.9964953 0.9975390 0.9898931 0.9956706 0.9996855 0.9921195 0.9978888
    ## [3,] 0.9929146 0.9960867 0.9887517 0.9927584 1.0000000 0.9871556 0.9955871
    ## [4,] 0.9886707 0.9946529 0.9876106 0.9898361 1.0000000 0.9814970 0.9932916
    ## [5,] 0.9844501 0.9932369 0.9859981 0.9869039 1.0000000 0.9758578 0.9909314
    ## [6,] 0.9802529 0.9918379 0.9843299 0.9839619 1.0000000 0.9702628 0.9880730
    ##              S         U         D         E         Z
    ## [1,] 0.9960640 0.9985739 0.9965708 0.9970807 0.9972562
    ## [2,] 0.9694649 0.9971973 0.9951050 0.9880594 0.9916547
    ## [3,] 0.9399824 0.9945680 0.9886465 0.9784743 0.9836214
    ## [4,] 0.9110366 0.9919118 0.9790672 0.9660765 0.9727788
    ## [5,] 0.8826282 0.9892536 0.9695393 0.9534496 0.9620563
    ## [6,] 0.8547583 0.9865934 0.9600628 0.9409275 0.9514497

``` r
filter_raw <- function(csv.file) {

  sn.cols <- c("126 Sn", "127n Sn", "127c Sn", "128n Sn", "128c Sn", "129n Sn",
               "129c Sn", "130n Sn", "130c Sn", "131n Sn", "131c Sn", "132n Sn",
               "132c Sn", "133n Sn", "133c Sn", "134n Sn", "134c Sn", "135n Sn")

  #Selecting columns needed for filtering and data analysis
  raw.columns <- c('Protein ID', 'Parsimony', 'Trimmed Peptide', sn.cols,
                   'Theo m/z', 'z', 'Isolation m/z', "Isolation Specificity")

  
  cat("Reading",  tail(unlist(strsplit(csv.file, "/")), n = 1), "...\n\n")  
    
  data.df <- fread(input = csv.file, select = raw.columns) %>%
    as.data.frame
  
  cat("Initial rows:\t\t\t\t\t", nrow(data.df), "\n")

  #'*Parsimony, REV sequences, contaminants, and oxM*
  data.df <- data.df[data.df$Parsimony %in% c("R", "U"),] #Drop NA Parsimony
  data.df <- data.df[!grepl("#", data.df$`Protein ID`),] #Remove reverse sequences
  data.df <- data.df[!grepl('contaminant', data.df$`Protein ID`),] #Remove contaminants
  data.df <- data.df[!grepl("\\*", data.df$`Trimmed Peptide`),] #Remove oxM
  
  cat("After Initial Parsimony/REV/Contam/oxM Filter:\t", nrow(data.df), "\n")

  data.df <- data.df[data.df$z==2 | data.df$z==3, ]

  cat("After Charge State Filter:\t\t\t", nrow(data.df), "\n")  
      
  # data.df <- data.df[data.df$`Isolation Specificity` > 0.75,]

  cat("After Isolation Specificity Filter:\t\t", nrow(data.df), "\n")
  
  #'*M0 Filter*  
  #Calculate theoretical/observed masses + "error"
  # ... Then select rows within error range
  data.df <- data.df[((data.df['z']*data.df['Isolation m/z']) >
                        ((data.df['z']*data.df['Theo m/z']) - 0.1)) &
                       ((data.df['z']*data.df['Isolation m/z']) <
                          ((data.df['z']*data.df['Theo m/z']) + 0.1)), ]
  
  cat("After M0 Filter:\t\t\t\t", nrow(data.df), "\n")  

  #'*Missed Cleavages*
  #Remove peptides with missed cleavages
  # ... Contains KR after removing last AA
  data.df <- data.df[!(grepl("K",
                        str_sub(data.df$`Trimmed Peptide`, end=-2))),]
  data.df <- data.df[!(grepl("R",
                        str_sub(data.df$`Trimmed Peptide`, end=-2))),]
  cat("After removing Missed Cleavages:\t\t", nrow(data.df), "\n")
  
  #'*TMTpro Signal*
  #Signal per channel
  # data.df["Sn_per_channel"] <- apply(data.df, 1, function(x) {
  #   return(all(as.numeric(x[sn.cols]) >= 4)) })
  # data.df <- data.df[data.df$Sn_per_channel==TRUE,]

  #Total signal
  data.df["sum_sn"] <- rowSums(data.df[sn.cols])
  data.df <- data.df[(data.df$z==2 & data.df$sum_sn>257)|
                       (data.df$z==3 & data.df$sum_sn>529),]  
    
  cat("After TMTpro Signal Filter:\t\t\t", nrow(data.df), "\n")  
  
  #Removing unnecessary columns and renaming data
  data.df[c("Parsimony", "Sn_per_channel", "z",
            "Theo m/z", "Isolation m/z", "Isolation Specificity")] <- NULL

  #Summing signal across peptides
  data.df <- data.df %>%
    group_by(`Protein ID`, `Trimmed Peptide`) %>%
    dplyr::summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)),
                     .groups = "drop") %>%
    as.data.frame()  
  
  # #Max signal within peptide group
  # data.df <- data.df %>% group_by(`Trimmed Peptide`) %>%
  #   filter(sum_sn == max(sum_sn)) %>% # ... Select row with max signal
  #   ungroup

  cat("Final Unique Peptides:\t\t\t\t", nrow(data.df), "\n\n\n")

  return(data.df) }
  
# filter_raw("Data/XLA_O18/XLA-O18-T5_RTS-18plex_LongestGFY.csv")
```

``` r
#'*MBT Timeseries Replicate 1*
XLA_T9.df <- filter_raw("Data/XLA_O18/XLA-O18-T9_RTS-18plex_Fractionated.csv")
```

    ## Reading XLA-O18-T9_RTS-18plex_Fractionated.csv ...
    ## 
    ## Initial rows:                     330099 
    ## After Initial Parsimony/REV/Contam/oxM Filter:    156871 
    ## After Charge State Filter:            156871 
    ## After Isolation Specificity Filter:       156871 
    ## After M0 Filter:              136312 
    ## After removing Missed Cleavages:      122048 
    ## After TMTpro Signal Filter:           83335 
    ## Final Unique Peptides:                48946

``` r
colnames(XLA_T9.df) <- c("Protein_ID", "Peptide",
                         "T0", "O1", "O2", "O3", "O4", "O5", "O6", "O7", "O8", "O9",
                         "O10", "O11", "N2", "N4", "N6", "N8", "N10", "N11", "sum_sn")

#'*MBT Timeseries Replicate 2*
XLA_T10.df <- filter_raw("Data/XLA_O18/XLA-O18-T10_RTS-18plex_Fractionated.csv")
```

    ## Reading XLA-O18-T10_RTS-18plex_Fractionated.csv ...
    ## 
    ## Initial rows:                     338537 
    ## After Initial Parsimony/REV/Contam/oxM Filter:    157751 
    ## After Charge State Filter:            157751 
    ## After Isolation Specificity Filter:       157751 
    ## After M0 Filter:              136926 
    ## After removing Missed Cleavages:      123614 
    ## After TMTpro Signal Filter:           87511 
    ## Final Unique Peptides:                49688

``` r
colnames(XLA_T10.df) <- c("Protein_ID", "Peptide",
                          "T0", "O1", "O2", "O3", "O4", "O5", "O6", "O7", "O8", "O9",
                          "O10", "O11", "N2", "N4", "N6", "N8", "N10", "N11", "sum_sn")
```

``` r
#Exponential decay function with two parameters
decay_fcn <- function(t, k1, k2) { k1* exp(-k2*t) }

#Function to fit all peptide decay
peptide_AA_fit <- function(peptide) {
  # Z not used
  # CHLA -> X
  # YN -> U
  mod_pep <- gsub("[CLA]", "X", peptide)
  mod_pep <- gsub("[N]", "U", mod_pep)  
  
  pep_table <- table(strsplit(mod_pep,"")) #Splitting to structure matrix

  #Creating empty matrix of the entire amino acid table
  peptide_matrix <- matrix(0, nrow = nrow(AA_matrix_fit),
                           ncol = ncol(AA_matrix_fit) - 1)  
  colnames(peptide_matrix) <- colnames(AA_matrix_fit)[-1]

  #Replacing AA column with number of amino acids
  for (i in names(pep_table)) { peptide_matrix[,i] <- as.numeric(pep_table[i]) }

  calc_per_AA <- AA_matrix_fit[,-1] ^ peptide_matrix
  pep_deg.tot <- apply(calc_per_AA, 1, prod)
  pep_deg.tot <- pep_deg.tot / pep_deg.tot[1]  

  # Fitting using nlsLM which is less limited for starting values
  #   decay_fcn is an exponential decay function with two parameters defined above
  decay.fit <- nlsLM(pep_deg.tot ~ decay_fcn(seq(0, 2886, by = 1), 1, k2),
                     start = list(k2=0.04),
                     lower = c(0),
                     upper = c(Inf),
                     control = list(maxiter = 1000))
  fit.values <- as.numeric(coef(decay.fit)) #Extracting fit  

  return(data.frame(t(c(1, fit.values)))) }    


#All unique identified peptides
mbt_xla_peps <- unique(c(XLA_T9.df$Peptide, XLA_T10.df$Peptide))

#Applying function to fit each peptide
mbt_xla_pep_fits <- data.frame(Peptide = mbt_xla_peps)

#------------------------------------------------------------

#Setup cluster
cl <- makeCluster(num_cores)

#Export necessary variables/functions to the cluster nodes
clusterExport(cl,
              varlist=c("mbt_xla_peps", "peptide_AA_fit", "decay_fcn",
                        "AA_matrix_fit"))
invisible(clusterEvalQ(cl, library(minpack.lm)))

mbt_xla_pep_fits <- cbind(mbt_xla_pep_fits,
  do.call(rbind, lapply(mbt_xla_peps, peptide_AA_fit)))
colnames(mbt_xla_pep_fits)[2:3] <- c("Pep_k1", "Pep_k2")

stopCluster(cl) #Stop Cluster

#------------------------------------------------------------

head(mbt_xla_pep_fits)
```

    ##                  Peptide Pep_k1     Pep_k2
    ## 1              AQQNNVEHK      1 0.04624886
    ## 2            DVVFEFPEFQL      1 0.06807309
    ## 3               EIEVGAGR      1 0.06360606
    ## 4               ELNITAAK      1 0.04928538
    ## 5               HVVFIAQR      1 0.03990100
    ## 6 TLTAVHDAILEDLVYPSEIVGR      1 0.15254588

``` r
filtered_XLA_T9.df <- merge(XLA_T9.df,
                            mbt_xla_pep_fits[c("Peptide", "Pep_k1", "Pep_k2")],
                            by="Peptide") %>%
  arrange(Protein_ID) %>% select(Protein_ID, everything())
filtered_XLA_T10.df <- merge(XLA_T10.df,
                             mbt_xla_pep_fits[c("Peptide", "Pep_k1", "Pep_k2")],
                             by="Peptide") %>%
  arrange(Protein_ID) %>% select(Protein_ID, everything())

write.csv(mbt_xla_pep_fits,
          "Files/Fits/XLA_O18-MBT-TS_Theo-Peptide-Decay.csv",
          row.names = FALSE)
write.csv(filtered_XLA_T9.df,
          "Data/XLA_O18/XLA_O18-T9_Filtered_Data.csv",
          row.names = FALSE)
write.csv(filtered_XLA_T10.df,
          "Data/XLA_O18/XLA_O18-T10_Filtered_Data.csv",
          row.names = FALSE)
```

``` r
min_pep_k2 <- mbt_xla_pep_fits %>%
  filter(min(Pep_k2)==Pep_k2)

max_pep_k2 <- mbt_xla_pep_fits %>%
  filter(max(Pep_k2)==Pep_k2)

med_pep_k2 <- mbt_xla_pep_fits %>%
  mutate(diff = abs(Pep_k2 - median(Pep_k2))) %>%
  filter(diff == min(diff)) %>%
  select(-diff)
med_pep_k2 <- med_pep_k2[1,]

example_fits <- rbind(min_pep_k2, med_pep_k2, max_pep_k2)
example_fits["Type"] <- c("Min", "Median", "Max")

ex.fit_data <- lapply(example_fits$Peptide, function(peptide) {
  
  ex.fit <- example_fits[example_fits$Peptide==peptide,]

  curve_values <- decay_fcn(seq(0,2886,1), ex.fit$Pep_k1, ex.fit$Pep_k2)  

  example.df <- data.frame(Time=seq(0,2886,1)/60,
                           Fit_values=curve_values,
                           Type=rep(ex.fit$Type, length(curve_values)))
  
  return(example.df) })

ex.fit_data <- do.call(rbind, ex.fit_data) 

ex.p <- ggplot() +
  geom_line(data=ex.fit_data[ex.fit_data$Type=="Median",], 
             aes(x=Time, y=Fit_values), color="#009E73", size=1.5) +
  geom_line(data=ex.fit_data[ex.fit_data$Type=="Min",], 
             aes(x=Time, y=Fit_values), color="#CC79A7", size=1.5) +
  
  theme_bw() +
  labs(x="Hours after Labeling", y="Probability of light peptide") +
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
```

``` r
round(log(2)/example_fits$Pep_k2,0)
```

    ## [1] 37  9  2

``` r
# #Uncomment for new image!
# tiff("Graphs/Frog_AA_Incorporation/XLA_T2_Example_Fits.tiff", units="in",
#      width=4.5, height=4.5, res=300)

ex.p
```

![](figures/Frog_AA_Incorporation_MBT-T2/Example%20fits%20plots-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
example_fits
```

    ##               Peptide Pep_k1     Pep_k2   Type
    ## 1             HPYWPHK      1 0.01893787    Min
    ## 2          LQLFSTQDGR      1 0.08032852 Median
    ## 3 DSSSSPASTASSGSSASLK      1 0.34424193    Max

``` r
#Function to fit all peptide decay
peptide_AA_fit <- function(peptide) {
  # Z not used
  # CHLA -> X
  # YN -> U
  mod_pep <- gsub("[CLA]", "X", peptide)
  mod_pep <- gsub("[N]", "U", mod_pep)  
  
  pep_table <- table(strsplit(mod_pep,"")) #Splitting to structure matrix

  #Creating empty matrix of the entire amino acid table
  peptide_matrix <- matrix(0, nrow = nrow(AA_matrix_fit),
                           ncol = ncol(AA_matrix_fit) - 1)  
  colnames(peptide_matrix) <- colnames(AA_matrix_fit)[-1]

  #Replacing AA column with number of amino acids
  for (i in names(pep_table)) { peptide_matrix[,i] <- as.numeric(pep_table[i]) }

  calc_per_AA <- AA_matrix_fit[,-1] ^ peptide_matrix
  pep_deg.tot <- apply(calc_per_AA, 1, prod)
  pep_deg.tot <- pep_deg.tot / pep_deg.tot[1]  

  # # Fitting using nlsLM which is less limited for starting values
  # #   decay_fcn is an exponential decay function with two parameters defined above
  # decay.fit <- nlsLM(pep_deg.tot ~ decay_fcn(seq(0, 2886, by = 1), 1, k2),
  #                    start = list(k2=0.04),
  #                    lower = c(0),
  #                    upper = c(Inf),
  #                    control = list(maxiter = 1000))
  # fit.values <- as.numeric(coef(decay.fit)) #Extracting fit  

  # plot(seq(0, 2886, by = 1)/60, pep_deg.tot,
  #      xlim=c(0,8), ylim=c(0,1), pch=16)
  # lines(seq(0,8*60,1)/60, decay_fcn(seq(0,8*60,1), 1, fit.values[1]),
  #       col="red")

  return(pep_deg.tot) }

median_pep_residues <- data.frame(Time = seq(0, 2886, by = 1)/60,
                                  Fit_values = peptide_AA_fit("LQLFSTQDGR"))
longest_pep_residues <- data.frame(Time = seq(0, 2886, by = 1)/60,
                                   Fit_values = peptide_AA_fit("HPYWPHK"))

res.p <- ggplot() +
  geom_point(data=longest_pep_residues, 
             aes(x=Time, y=Fit_values), color="#0072B2", size=2,
             shape=1) +
  geom_point(data=median_pep_residues, 
             aes(x=Time, y=Fit_values), color="#E7298A", size=2,
             shape=1) +

  theme_bw() +
  labs(x="Hours after Labeling", y="Probability of light peptide") +
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

write.csv(median_pep_residues, "Files/Fits/AA_examples/Median_Gast_residue-prob.csv",
          row.names = FALSE)

# #Uncomment for new image!
# tiff("Graphs/Frog_AA_Incorporation/XLA_T2_Example_Fits.tiff", units="in",
#      width=4, height=4, res=300)

res.p
```

![](figures/Frog_AA_Incorporation_MBT-T2/Plots%20of%20example%20fits-1.png)<!-- -->

``` r
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
    ## [1] parallel  splines   stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ## [1] minpack.lm_1.2-4  data.table_1.18.4 patchwork_1.3.2   ggplot2_4.0.3    
    ## [5] stringr_1.6.0     dplyr_1.2.0       tidyr_1.3.2      
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.6       compiler_4.5.3     tidyselect_1.2.1   scales_1.4.0      
    ##  [5] yaml_2.3.12        fastmap_1.2.0      R6_2.6.1           labeling_0.4.3    
    ##  [9] generics_0.1.4     knitr_1.51         tibble_3.3.1       pillar_1.11.1     
    ## [13] RColorBrewer_1.1-3 rlang_1.1.7        stringi_1.8.7      xfun_0.58         
    ## [17] S7_0.2.1           otel_0.2.0         cli_3.6.5          withr_3.0.3       
    ## [21] magrittr_2.0.4     digest_0.6.39      grid_4.5.3         rstudioapi_0.19.0 
    ## [25] lifecycle_1.0.5    vctrs_0.7.1        evaluate_1.0.5     glue_1.8.0        
    ## [29] farver_2.1.2       rmarkdown_2.31     purrr_1.2.2        tools_4.5.3       
    ## [33] pkgconfig_2.0.3    htmltools_0.5.9
