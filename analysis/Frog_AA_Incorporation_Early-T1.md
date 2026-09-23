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
library(scam)
```

    ## This is scam 1.2-22.

``` r
library(parallel)

knitr::opts_chunk$set(fig.path = "figures/Frog_AA_Incorporation_Early_T1/")
num_cores <- detectCores() - 1  # Save one core for system stability
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
AA_T1.time <- c(0, 57, 120, 182, 296, 475, 657, 1080, 1440) #Add 149 to get hpf
```

``` r
merged_frog_table <- read.csv("Data/XLA_AAA/XLA_O18-T1-T2_results-merge_AA.csv")
merged_frog_table <- merged_frog_table[!grepl("blank", merged_frog_table$Sample),]

early_exp_id <- "col009b"
mbt_exp_id <- "col009b2"

early_AA_table <- merged_frog_table[grepl(paste(early_exp_id, "_", sep=""),
                                          merged_frog_table$Sample),]

cat("Missing AAs: ", 
    aa.profiles$Name[!aa.profiles$Name %in% unique(early_AA_table$compound)],
    "\n")
```

    ## Missing AAs:  cysteine histidine leucine tyrosine

``` r
#'*Looping through unique amino acids and reformatting table*
early_AA_matrix <- sapply(unique(early_AA_table$compound), function(x){
  aa.data <- early_AA_table[early_AA_table$compound==x,] #selecting AA

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

early_AA_matrix <- bind_rows(early_AA_matrix)

#'*Calculating median values for missing AAs*
early_med_values <- lapply(paste("T", seq(0,8,1), sep=""), function(time){
  #Calculating median for timepoints collected
  t_data <- early_AA_matrix[grepl(time, early_AA_matrix$Sample),]
  
  t_data <- t_data[!t_data$compound %in% c("alanine", "threonine",
                                           "tyrosine", "serine",
                                           "asparagine"),]
  
  sub_2 <- t_data[t_data$Total==2,]      
  med_2 <- apply(sub_2[c("C12 PARENT", "O18-label-1",
                         "O18-label-2")], 2, median)
  med_2 <- med_2 / sum(med_2) #Necessary to remake row sum = 1

  
  sub_3 <- t_data[t_data$Total==3,]      
  med_3 <- apply(sub_3[c("C12 PARENT", "O18-label-1",
                         "O18-label-2", "O18-label-3")], 2, median)
  med_3 <- med_3 / sum(med_3) #Necessary to remake row sum = 1
  
  sub_4 <- t_data[t_data$Total==4,]      
  med_4 <- apply(sub_4[c("C12 PARENT", "O18-label-1",
                         "O18-label-2", "O18-label-3",
                         "O18-label-4")], 2, median)
  med_4 <- med_4 / sum(med_4) #Necessary to remake row sum = 1
  

  row_out <- data.frame(Name=rep("median", 3),
                        AA=c("X", "U", "Z"),
                        Total=c(2,3,4),
                        Sample=paste("col009b_", time, sep=""))
  
  return(cbind(row_out, bind_rows(med_2, med_3, med_4))) })

#Merging and appending to table of all amino acids
early_med_values <- do.call(rbind, early_med_values)
early_AA_matrix <- merge(aa.profiles[c("Name", "AA")],
                         early_AA_matrix, by.x="Name",
                         by.y="compound")

early_AA_matrix <- bind_rows(early_AA_matrix, early_med_values)
early_AA_matrix$Name <- sub("^(\\w)", "\\U\\1", early_AA_matrix$Name, perl = TRUE)

head(early_AA_matrix)
```

    ##      Name AA Total     Sample C12 PARENT O18-label-1 O18-label-2 O18-label-3
    ## 1 Alanine  A     2 col009b_T0  0.9981057 0.001894258   0.0000000          NA
    ## 2 Alanine  A     2 col009b_T1  0.5469703 0.137989921   0.3150398          NA
    ## 3 Alanine  A     2 col009b_T2  0.5186425 0.121910796   0.3594467          NA
    ## 4 Alanine  A     2 col009b_T3  0.6286285 0.102806275   0.2685652          NA
    ## 5 Alanine  A     2 col009b_T4  0.6308999 0.101144723   0.2679554          NA
    ## 6 Alanine  A     2 col009b_T5  0.7079944 0.090237951   0.2017676          NA
    ##   O18-label-4
    ## 1          NA
    ## 2          NA
    ## 3          NA
    ## 4          NA
    ## 5          NA
    ## 6          NA

``` r
write.csv(early_AA_matrix, "Data/XLA_O18/Isotopic_Envelopes/Frog_Early_AA_Matrix.csv",
          row.names = FALSE)
```

``` r
#Frog Essential AAs
measured_AAs <- unique(early_AA_matrix$AA)
#Removed asparagine -> no measured O18-label-3
#Removed alanine -> interference on run

no_Oxy_AA <- c("A", "R", "C", "G", "V",
               "H", "I", "L", "K",
               "M", "F", "P", "W",
               "T", "Y", "X") #Threonine and Tyrosine included here because of model
no_Oxy_AA <- no_Oxy_AA[no_Oxy_AA %in% setdiff(measured_AAs, c("A"))]

one_Oxy_AA <- c("Q", "S", "N", "U") 
one_Oxy_AA <- one_Oxy_AA[one_Oxy_AA %in% setdiff(measured_AAs, c("N"))]

two_Oxy_AA <- c("D", "E", "Z")
two_Oxy_AA <- two_Oxy_AA[two_Oxy_AA %in% measured_AAs]
```

``` r
no_Oxy.splines <- do.call(rbind, lapply(no_Oxy_AA, function(AA) {

  # Example setup
  interpolated_times <- seq(0, 1440, by = 1)
  sub.data <- early_AA_matrix[early_AA_matrix$AA == AA,]
  sub.data <- sub.data %>% arrange(Sample)
  
  # Define knot positions manually
  knot_positions <- c(90, 300)

  # Fit spline models
  c12_fit <- lm(`C12 PARENT` ~ bs(AA_T1.time, degree = 3,
                                  knots = knot_positions), data = sub.data)
  o18_1_fit <- lm(`O18-label-1` ~ bs(AA_T1.time, degree = 3,
                                     knots = knot_positions), data = sub.data)
  o18_2_fit <- lm(`O18-label-2` ~ bs(AA_T1.time, degree = 3,
                                     knots = knot_positions), data = sub.data)

  # Predict spline values
  c12_pred <- predict(c12_fit, newdata = data.frame(AA_T1.time = interpolated_times))
  o18_1_pred <- predict(o18_1_fit, newdata = data.frame(AA_T1.time = interpolated_times))
  o18_2_pred <- predict(o18_2_fit, newdata = data.frame(AA_T1.time = interpolated_times))


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
  #   geom_point(data=data.frame(x=AA_T1.time,
  #                              y=sub.data$`C12 PARENT`),
  #              aes(x=x,y=y), size=3) +
  #   geom_line(data=spline_df[spline_df$Type=="C12",],
  #              aes(x=Time,y=Value), size=1) +
  #   theme_bw() +
  #   coord_cartesian(ylim=c(0,1))
  # 
  # print(x)
  
  return(spline_df)
  
}))
```

``` r
one_Oxy.splines <- do.call(rbind, lapply(one_Oxy_AA,
                                         function(AA) {

  # Example setup
  interpolated_times <- seq(0, 1440, by = 1)
  sub.data <- early_AA_matrix[early_AA_matrix$AA == AA,]
  sub.data <- sub.data %>% arrange(Sample)

  knot_positions <- c(90,300)    

  # Fit spline models
  c12_fit <- lm(`C12 PARENT` ~ bs(AA_T1.time, degree = 3,
                                  knots = knot_positions), data = sub.data)
  o18_1_fit <- lm(`O18-label-1` ~ bs(AA_T1.time, degree = 3,
                                     knots = knot_positions), data = sub.data)
  o18_2_fit <- lm(`O18-label-2` ~ bs(AA_T1.time, degree = 3,
                                     knots = knot_positions), data = sub.data)
  o18_3_fit <- lm(`O18-label-3` ~ bs(AA_T1.time, degree = 3,
                                     knots = knot_positions), data = sub.data)
  
  # Predict spline values
  c12_pred <- predict(c12_fit, newdata = data.frame(AA_T1.time = interpolated_times))
  o18_1_pred <- predict(o18_1_fit, newdata = data.frame(AA_T1.time = interpolated_times))
  o18_2_pred <- predict(o18_2_fit, newdata = data.frame(AA_T1.time = interpolated_times))
  o18_3_pred <- predict(o18_3_fit, newdata = data.frame(AA_T1.time = interpolated_times))

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
  #   geom_point(data=data.frame(x=AA_T1.time,
  #                              y=sub.data$`C12 PARENT`),
  #              aes(x=x,y=y), size=3) +
  #   geom_line(data=spline_df[spline_df$Type=="C12",],
  #              aes(x=Time,y=Value), size=1) +
  #   theme_bw() +
  #   coord_cartesian(ylim=c(0,1))
  # 
  # print(x)
    
  return(spline_df)
  
}))
```

``` r
two_Oxy.splines <- do.call(rbind, lapply(two_Oxy_AA,
                                         function(AA) {

  # Example setup
  interpolated_times <- seq(0, 1440, by = 1)
  sub.data <- early_AA_matrix[early_AA_matrix$AA == AA,]
  sub.data <- sub.data %>% arrange(Sample)

  knot_positions <- c(300,500)    

  # Fit spline models
  c12_fit <- lm(`C12 PARENT` ~ bs(AA_T1.time, degree = 3,
                                  knots = knot_positions), data = sub.data)
  o18_1_fit <- lm(`O18-label-1` ~ bs(AA_T1.time, degree = 3,
                                     knots = knot_positions), data = sub.data)
  o18_2_fit <- lm(`O18-label-2` ~ bs(AA_T1.time, degree = 3,
                                     knots = knot_positions), data = sub.data)
  o18_3_fit <- lm(`O18-label-3` ~ bs(AA_T1.time, degree = 3,
                                     knots = knot_positions), data = sub.data)
  o18_4_fit <- lm(`O18-label-4` ~ bs(AA_T1.time, degree = 3,
                                     knots = knot_positions), data = sub.data)
    
  # Predict spline values
  c12_pred <- predict(c12_fit, newdata = data.frame(AA_T1.time = interpolated_times))
  o18_1_pred <- predict(o18_1_fit, newdata = data.frame(AA_T1.time = interpolated_times))
  o18_2_pred <- predict(o18_2_fit, newdata = data.frame(AA_T1.time = interpolated_times))
  o18_3_pred <- predict(o18_3_fit, newdata = data.frame(AA_T1.time = interpolated_times))
  o18_4_pred <- predict(o18_4_fit, newdata = data.frame(AA_T1.time = interpolated_times))
  
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
  #   geom_point(data=data.frame(x=AA_T1.time,
  #                              y=sub.data$`C12 PARENT`),
  #              aes(x=x,y=y), size=3) +
  #   geom_line(data=spline_df[spline_df$Type=="C12",],
  #              aes(x=Time,y=Value), size=1) +
  #   theme_bw() +
  #   coord_cartesian(ylim=c(0,1))
  # 
  # print(x)
    
  return(spline_df)
  
}))
```

``` r
all_splines.merge <- rbind(no_Oxy.splines, one_Oxy.splines, two_Oxy.splines)
all_splines.merge <- do.call(rbind,
              lapply(unique(all_splines.merge$AA), function(AA) {

  #Looping through each time to verify rowsum=1 for each AA
  c_df <- lapply(seq(0, 1440, by = 1), function(t) {
    sub_t <- all_splines.merge[all_splines.merge$AA==AA &
                              all_splines.merge$Time==t,]
    sub_t["Value"] <- sub_t$Value / sum(sub_t$Value)
    
    return(sub_t) })
  
  c_df <- do.call(rbind, c_df)

  return(c_df) }))

head(all_splines.merge)
```

    ##      AA Labels Time  Type       Value
    ## 1     R      2    0   C12 0.996670681
    ## 1442  R      2    0 O18_1 0.001947345
    ## 2883  R      2    0 O18_2 0.001381974
    ## 2     R      2    1   C12 0.989931000
    ## 1443  R      2    1 O18_1 0.008474179
    ## 2884  R      2    1 O18_2 0.001594821

``` r
all_spline.p <- lapply(no_Oxy_AA, function(ind_AA){
  
  sub.data <- early_AA_matrix[early_AA_matrix$AA == ind_AA,]
  sub.data <- sub.data %>% arrange(Sample)
  sub.data["Time"] <- (AA_T1.time+149)/60
  
  spline.data <- no_Oxy.splines[no_Oxy.splines$AA==ind_AA,]
  spline.data["Time"] <- (spline.data$Time+149)/60
  
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
          plot.title = element_text(size=24, hjust=0.5, colour="black"),          
          # plot.title = element_text(size=24, hjust=0.5, colour="#5E4FA2"),
          # axis.title = element_text(size=22),
          axis.title = element_blank(),                    
          panel.border = element_rect(linewidth=2),
          panel.grid.major = element_line(color = "grey70", linewidth = 0.25),
          panel.grid.minor = element_line(color = "grey70", linewidth = 0.25),
          aspect.ratio = 1,
          legend.position = "none") +
    coord_cartesian(xlim=c(2,(1440+149)/60), ylim=c(0,1)) +
    scale_x_continuous(breaks=seq(5,25,5))

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

  sub.data <- early_AA_matrix[early_AA_matrix$AA == ind_AA,]
  sub.data <- sub.data %>% arrange(Sample)
  sub.data["Time"] <- (AA_T1.time+149)/60

  spline.data <- one_Oxy.splines[one_Oxy.splines$AA==ind_AA,]
  spline.data["Time"] <- (spline.data$Time+149)/60

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
          plot.title = element_text(size=24, hjust=0.5, colour="black"),
          # plot.title = element_text(size=24, hjust=0.5, colour="#E7298A"),
          # axis.title = element_text(size=22),
          axis.title = element_blank(),
          panel.border = element_rect(linewidth=2),
          panel.grid.major = element_line(color = "grey70", linewidth = 0.25),
          panel.grid.minor = element_line(color = "grey70", linewidth = 0.25),
          aspect.ratio = 1,
          legend.position = "none") +
    coord_cartesian(xlim=c(2,(1440+149)/60), ylim=c(0,1)) +
    scale_x_continuous(breaks=seq(5,25,5))

  return(AA_spline.p) }) )


all_spline.p <- c(all_spline.p,
                  lapply(two_Oxy_AA, function(ind_AA){

  sub.data <- early_AA_matrix[early_AA_matrix$AA == ind_AA,]
  sub.data <- sub.data %>% arrange(Sample)
  sub.data["Time"] <- (AA_T1.time+149)/60

  spline.data <- two_Oxy.splines[two_Oxy.splines$AA==ind_AA,]
  spline.data["Time"] <- (spline.data$Time+149)/60

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
          plot.title = element_text(size=24, hjust=0.5, colour="black"),
          # plot.title = element_text(size=24, hjust=0.5, colour="#F0E442"),
          # axis.title = element_text(size=22),
          axis.title = element_blank(),
          panel.border = element_rect(linewidth=2),
          panel.grid.major = element_line(color = "grey70", linewidth = 0.25),
          panel.grid.minor = element_line(color = "grey70", linewidth = 0.25),
          aspect.ratio = 1,
          legend.position = "none") +
    coord_cartesian(xlim=c(2,(1440+149)/60), ylim=c(0,1)) +
    scale_x_continuous(breaks=seq(5,25,5))

  return(AA_spline.p) }) )

names(all_spline.p) <- c(no_Oxy_AA, one_Oxy_AA,
                         two_Oxy_AA)

# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=9, height=5, res=300)

# all_spline.p[["R"]] | all_spline.p[["Q"]] | all_spline.p[["E"]]

# dev.off() #Uncomment for new image!
```

``` r
# #Uncomment for new image!
# tiff("Graphs/Frog_AA_Incorporation/Median_XLA_T1_AA-Decay.tiff", units="in",
#      width=9, height=5, res=300)

all_spline.p[["X"]] | all_spline.p[["U"]] | all_spline.p[["Z"]]
```

![](figures/Frog_AA_Incorporation_Early_T1/Median%20plots%20of%20AAs-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=13, height=12, res=300)

(all_spline.p[["R"]] | all_spline.p[["G"]] | all_spline.p[["V"]] | all_spline.p[["I"]] ) /
  (all_spline.p[["K"]] | all_spline.p[["M"]] | all_spline.p[["F"]] | all_spline.p[["P"]]) /
  (all_spline.p[["W"]] | all_spline.p[["T"]] | all_spline.p[["T"]] | all_spline.p[["T"]]) /
  (all_spline.p[["S"]] | all_spline.p[["Q"]] | all_spline.p[["D"]] | all_spline.p[["E"]])
```

![](figures/Frog_AA_Incorporation_Early_T1/Supp%20Figure%20of%20all%20amino%20acids-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
AA_matrix_fit <- data.frame(Time=seq(0, 1440, by = 1))

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

    ##      Time         R         G         V         I         K         M         F
    ## [1,]    0 0.9983340 0.9998720 0.9989307 0.9982563 0.9980076 0.9986071 0.9987260
    ## [2,]    1 0.9949528 0.9996709 0.9974590 0.9943551 0.9983918 0.9875786 0.9955145
    ## [3,]    2 0.9915838 0.9993930 0.9959728 0.9904905 0.9987214 0.9764687 0.9923298
    ## [4,]    3 0.9882270 0.9990971 0.9944723 0.9866625 0.9989971 0.9654851 0.9891718
    ## [5,]    4 0.9848826 0.9987832 0.9929578 0.9828711 0.9992097 0.9546294 0.9860404
    ## [6,]    5 0.9815503 0.9984516 0.9914295 0.9791162 0.9991486 0.9439031 0.9829356
    ##              P         W         T         X         Q         S         U
    ## [1,] 0.9994393 0.9915790 0.9984928 0.9987204 0.9974083 0.9971706 0.9974083
    ## [2,] 0.9942146 0.9897164 0.9974798 0.9961691 0.9961993 0.9778029 0.9961993
    ## [3,] 0.9888408 0.9878959 0.9964074 0.9935745 0.9949542 0.9586455 0.9949542
    ## [4,] 0.9835350 0.9861172 0.9953032 0.9909920 0.9936735 0.9398432 0.9936735
    ## [5,] 0.9782974 0.9843801 0.9941674 0.9884215 0.9923578 0.9213888 0.9923578
    ## [6,] 0.9731278 0.9826843 0.9930006 0.9858630 0.9910074 0.9032754 0.9910074
    ##              D         E         Z
    ## [1,] 0.9884870 0.9881641 0.9891702
    ## [2,] 0.9866716 0.9854335 0.9868107
    ## [3,] 0.9831089 0.9820241 0.9832349
    ## [4,] 0.9786308 0.9759312 0.9778546
    ## [5,] 0.9741746 0.9698748 0.9725025
    ## [6,] 0.9697399 0.9638545 0.9671783

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
#'*Early Timeseries Replicate 1*
XLA_T5.df <- filter_raw("Data/XLA_O18/XLA-O18-T5_RTS-18plex_LongestGFY.csv")
```

    ## Reading XLA-O18-T5_RTS-18plex_LongestGFY.csv ...
    ## 
    ## Initial rows:                     325659 
    ## After Initial Parsimony/REV/Contam/oxM Filter:    155215 
    ## After Charge State Filter:            155215 
    ## After Isolation Specificity Filter:       155215 
    ## After M0 Filter:              132285 
    ## After removing Missed Cleavages:      117340 
    ## After TMTpro Signal Filter:           88585 
    ## Final Unique Peptides:                49730

``` r
colnames(XLA_T5.df) <- c("Protein_ID", "Peptide",
                         "T0", "O1", "O2", "O3", "O4", "O5", "O7", "O9", "O10",
                         "O11", "O12", "N1", "N3", "N4", "N7", "N9",
                         "N11", "N12", "sum_sn")

#'*Early Timeseries Replicate 2*
XLA_T6.df <- filter_raw("Data/XLA_O18/XLA-O18-T6_RTS-18plex_LongestGFY.csv")
```

    ## Reading XLA-O18-T6_RTS-18plex_LongestGFY.csv ...
    ## 
    ## Initial rows:                     331914 
    ## After Initial Parsimony/REV/Contam/oxM Filter:    160072 
    ## After Charge State Filter:            160072 
    ## After Isolation Specificity Filter:       160072 
    ## After M0 Filter:              136323 
    ## After removing Missed Cleavages:      118811 
    ## After TMTpro Signal Filter:           88389 
    ## Final Unique Peptides:                50417

``` r
colnames(XLA_T6.df) <- c("Protein_ID", "Peptide",
                         "T0", "O1", "O2", "O3", "O4", "O5", "O7", "O9", "O10",
                         "O11", "O12", "N1", "N3", "N4", "N7", "N9",
                         "N11", "N12", "sum_sn")

head(XLA_T5.df)
```

    ##   Protein_ID                Peptide       T0       O1       O2       O3
    ## 1 XBgroup104              AQQNNVEHK  65.5842  87.7493  53.0013  75.5488
    ## 2 XBgroup104            DVVFEFPEFQL 173.6741 163.5202 161.0787 168.0296
    ## 3 XBgroup104               EIEVGAGR  24.1504  20.5260  29.0353  22.6757
    ## 4 XBgroup104               ELNITAAK 202.1141 212.6018 155.3938 200.8254
    ## 5 XBgroup104               HVVFIAQR 898.3401 984.1120 732.6715 937.9363
    ## 6 XBgroup104 TLTAVHDAILEDLVYPSEIVGR 374.1061 344.2254 314.8003 370.9311
    ##         O4       O5       O7       O9      O10      O11      O12       N1
    ## 1  69.2270  65.6506  41.1070  63.2307  39.0110  53.6571  36.0788  46.5791
    ## 2 142.9165 132.3541 127.6454 140.0324 106.9388 126.3283 122.1724 148.2365
    ## 3  14.6748  14.4232  14.4842  20.1995  16.9965  12.5243  19.0686  17.2215
    ## 4 178.0862 178.2336 132.2850 187.8982 135.9825 154.4422 120.8135 159.9677
    ## 5 783.9324 798.9270 680.1434 825.4029 584.2716 775.6975 597.4344 716.6382
    ## 6 325.9760 303.9890 278.5779 325.4775 243.8627 299.5034 278.4892 274.8966
    ##         N3       N4       N7       N9      N11      N12     sum_sn
    ## 1  64.2125  57.7312  60.4782  52.5800  32.1251  34.5394   998.0913
    ## 2 163.5156 136.2404 144.0379 133.9678 110.1057 115.1145  2515.9089
    ## 3  20.2093  16.1441  21.8170  14.6610  16.2793  19.4522   334.5429
    ## 4 178.7240 174.3278 167.1650 145.8722 101.6695 123.2108  2909.6133
    ## 5 845.3745 730.5324 747.1786 706.4852 496.8501 599.9997 13441.9278
    ## 6 358.5770 253.5004 316.0018 289.1671 226.9415 252.7842  5431.8072

``` r
head(XLA_T6.df)
```

    ##   Protein_ID     Peptide       T0       O1       O2       O3       O4       O5
    ## 1 XBgroup104   AQQNNVEHK 156.4828 135.7434 135.3264 159.1960 143.2329 150.0653
    ## 2 XBgroup104 DVVFEFPEFQL 354.8156 320.8339 343.5246 388.3506 365.9719 363.9255
    ## 3 XBgroup104    EIEVGAGR  43.2024  30.7374  41.1918  31.8476  42.2891  30.9209
    ## 4 XBgroup104    ELNITAAK  41.7534  39.9121  33.7158  44.0203  44.0107  40.7731
    ## 5 XBgroup104    HVVFIAQR 476.8073 459.9247 465.7442 490.8083 510.5731 541.1917
    ## 6 XBgroup104     MFSTSAK  35.9688  24.8003  36.8372  39.6326  35.2698  28.1325
    ##         O7       O9      O10      O11      O12       N1       N3       N4
    ## 1 141.5782 146.6952 116.2993 117.6814  91.7409 148.5631 143.8358 118.6437
    ## 2 327.8006 367.3011 282.5061 298.0320 227.2869 363.4394 327.1348 331.1368
    ## 3  29.0257  31.8337  32.1546  29.5228  24.8822  31.5533  28.2626  28.8684
    ## 4  40.9190  37.0442  38.4914  31.4106  29.2499  44.5847  34.0522  31.4433
    ## 5 396.4272 478.1729 402.6733 433.1485 320.6805 516.7107 456.4655 411.9480
    ## 6  28.8639  31.2527  15.5212  18.8204  16.0954  44.0017  25.8528  35.2672
    ##         N7       N9      N11      N12    sum_sn
    ## 1 135.2696 149.5379 116.0691 155.1960 2461.1570
    ## 2 328.7638 356.8906 290.3327 406.1308 6044.1777
    ## 3  34.0982  37.5765  21.9854  37.1742  587.1268
    ## 4  43.1018  39.6662  36.7613  42.8346  693.7446
    ## 5 511.2939 475.9609 372.2627 459.0928 8179.8862
    ## 6  36.3043  22.1903  44.0437 208.3953  727.2501

``` r
#Exporting for S-L Comp
dir.create("Data/S-L_Comp", showWarnings = FALSE, recursive = TRUE)
write.csv(AA_matrix_fit, "Data/S-L_Comp/AA_T1_Decay_Matrix.csv",
          row.names = FALSE)
```

``` r
#Exponential decay function with two parameters
decay_fcn <- function(t, k1, k2) { k1* exp(-k2*t) }

#Function to fit all peptide decay
peptide_AA_fit <- function(peptide) {
  # Z not used
  # CHLA -> X
  # YN -> U
  mod_pep <- gsub("[CHLA]", "X", peptide)
  mod_pep <- gsub("[YN]", "U", mod_pep)  

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
  decay.fit <- nlsLM(pep_deg.tot ~ decay_fcn(seq(0, 1440, by = 1), 1, k2),
                     start = list(k2=0.04),
                     lower = c(0),
                     upper = c(Inf),
                     control = list(maxiter = 1000))
  fit.values <- as.numeric(coef(decay.fit)) #Extracting fit  

  # plot(seq(0, 1440, by = 1)/60, pep_deg.tot,
  #      xlim=c(0,6), ylim=c(0,1), pch=16)
  # lines(seq(0,6*60,1)/60, decay_fcn(seq(0,6*60,1), 1, fit.values[1]),
  #       col="red")
  # print(decay.fit)
  # invokeRestart("abort")
  
  return(data.frame(t(c(1, fit.values)))) }

#All unique identified peptides
early_xla_peps <- unique(c(XLA_T5.df$Peptide, XLA_T6.df$Peptide))

#Applying function to fit each peptide
early_xla_pep_fits <- data.frame(Peptide = early_xla_peps)

#------------------------------------------------------------

#Setup cluster
cl <- makeCluster(num_cores)

#Export necessary variables/functions to the cluster nodes
clusterExport(cl,
              varlist=c("early_xla_peps", "peptide_AA_fit", "decay_fcn",
                        "AA_matrix_fit"))
invisible(clusterEvalQ(cl, library(minpack.lm)))

early_xla_pep_fits <- cbind(early_xla_pep_fits,
                            do.call(rbind, parLapply(cl, early_xla_peps, peptide_AA_fit)))
colnames(early_xla_pep_fits)[2:3] <- c("Pep_k1", "Pep_k2")

stopCluster(cl) #Stop Cluster

#------------------------------------------------------------

head(early_xla_pep_fits)
```

    ##                  Peptide Pep_k1     Pep_k2
    ## 1              AQQNNVEHK      1 0.02249180
    ## 2            DVVFEFPEFQL      1 0.03833003
    ## 3               EIEVGAGR      1 0.02556036
    ## 4               ELNITAAK      1 0.02181584
    ## 5               HVVFIAQR      1 0.01952115
    ## 6 TLTAVHDAILEDLVYPSEIVGR      1 0.08094138

``` r
filtered_XLA_T5.df <- merge(XLA_T5.df,
                            early_xla_pep_fits,
                            by="Peptide") %>%
  arrange(Protein_ID) %>% select(Protein_ID, everything())

filtered_XLA_T6.df <- merge(XLA_T6.df,
                            early_xla_pep_fits,
                            by="Peptide") %>%
  arrange(Protein_ID) %>% select(Protein_ID, everything())


write.csv(early_xla_pep_fits,
          "Files/Fits/XLA_O18-Early-TS_Theo-Peptide-Decay.csv",
          row.names = FALSE)
write.csv(filtered_XLA_T5.df,
          "Data/XLA_O18/XLA_O18-T5_Filtered_Data.csv",
          row.names = FALSE)
write.csv(filtered_XLA_T6.df,
          "Data/XLA_O18/XLA_O18-T6_Filtered_Data.csv",
          row.names = FALSE)
```

``` r
min_pep_k2 <- early_xla_pep_fits %>%
  filter(min(Pep_k2)==Pep_k2)

max_pep_k2 <- early_xla_pep_fits %>%
  filter(max(Pep_k2)==Pep_k2)

med_pep_k2 <- early_xla_pep_fits %>%
  mutate(diff = abs(Pep_k2 - median(Pep_k2))) %>%
  filter(diff == min(diff)) %>%
  select(-diff)
med_pep_k2 <- med_pep_k2[1,]

example_fits <- rbind(min_pep_k2, med_pep_k2, max_pep_k2)
example_fits["Type"] <- c("Min", "Median", "Max")

ex.fit_data <- lapply(example_fits$Peptide, function(peptide) {
  
  ex.fit <- example_fits[example_fits$Peptide==peptide,]

  curve_values <- decay_fcn(seq(0,1440,1), ex.fit$Pep_k1, ex.fit$Pep_k2)  

  example.df <- data.frame(Time=seq(0,1440,1)/60,
                           Fit_values=curve_values,
                           Type=rep(ex.fit$Type, length(curve_values)))
  
  return(example.df) })

ex.fit_data <- do.call(rbind, ex.fit_data) 

ex.p <- ggplot() +
  geom_line(data=ex.fit_data[ex.fit_data$Type=="Median",], 
             aes(x=Time, y=Fit_values), color="#E7298A", size=1.5) +
  geom_line(data=ex.fit_data[ex.fit_data$Type=="Min",], 
             aes(x=Time, y=Fit_values), color="#0072B2", size=1.5) +

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

    ## [1] 71 17  2

``` r
# #Uncomment for new image!
# tiff("Graphs/Frog_AA_Incorporation/XLA_T1_Example_Fits.tiff", units="in",
#      width=4, height=4, res=300)

ex.p
```

![](figures/Frog_AA_Incorporation_Early_T1/Example%20fits%20plots-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
example_fits
```

    ##                Peptide Pep_k1      Pep_k2   Type
    ## 1              WNGFGGK      1 0.009806907    Min
    ## 2       MLGNEICLCIVGNK      1 0.040398972 Median
    ## 3 SSESSSSSSSESSSSSESSR      1 0.332779069    Max

``` r
#Function to fit all peptide decay
peptide_AA_fit <- function(peptide) {
  # Z not used
  # CHLA -> X
  # YN -> U
  mod_pep <- gsub("[CHLA]", "X", peptide)
  mod_pep <- gsub("[YN]", "U", mod_pep)  

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
  # decay.fit <- nlsLM(pep_deg.tot ~ decay_fcn(seq(0, 1440, by = 1), 1, k2),
  #                    start = list(k2=0.04),
  #                    lower = c(0),
  #                    upper = c(Inf),
  #                    control = list(maxiter = 1000))
  # fit.values <- as.numeric(coef(decay.fit)) #Extracting fit

  # plot(seq(0, 1440, by = 1)/60, pep_deg.tot,
  #      xlim=c(0,8), ylim=c(0,1), pch=16)
  # lines(seq(0,8*60,1)/60, decay_fcn(seq(0,8*60,1), 1, fit.values[1]),
  #       col="red")

  return(pep_deg.tot) }

median_pep_residues <- data.frame(Time = seq(0, 1440, by = 1)/60,
                                  Fit_values = peptide_AA_fit("MLGNEICLCIVGNK"))
longest_pep_residues <- data.frame(Time = seq(0, 1440, by = 1)/60,
                                   Fit_values = peptide_AA_fit("WNGFGGK"))

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

write.csv(median_pep_residues, "Files/Fits/AA_examples/Median_2cell_residue-prob.csv",
          row.names = FALSE)

# #Uncomment for new image!
# tiff("Graphs/Frog_AA_Incorporation/XLA_T1_Example_Fits.tiff", units="in",
#      width=4, height=4, res=300)

res.p
```

![](figures/Frog_AA_Incorporation_Early_T1/Plots%20of%20example%20fits-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
median_pep_residues[median_pep_residues$Fit_values < 0.5,][1,]
```

    ##    Time Fit_values
    ## 19  0.3  0.4833439

``` r
longest_pep_residues[longest_pep_residues$Fit_values < 0.5,][1,]
```

    ##        Time Fit_values
    ## 74 1.216667  0.4951807

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
    ## [1] scam_1.2-22       minpack.lm_1.2-4  data.table_1.18.4 patchwork_1.3.2  
    ## [5] ggplot2_4.0.3     stringr_1.6.0     dplyr_1.2.0       tidyr_1.3.2      
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Matrix_1.7-5       gtable_0.3.6       compiler_4.5.3     tidyselect_1.2.1  
    ##  [5] scales_1.4.0       yaml_2.3.12        fastmap_1.2.0      lattice_0.22-9    
    ##  [9] R6_2.6.1           labeling_0.4.3     generics_0.1.4     knitr_1.51        
    ## [13] tibble_3.3.1       pillar_1.11.1      RColorBrewer_1.1-3 rlang_1.1.7       
    ## [17] stringi_1.8.7      xfun_0.58          S7_0.2.1           otel_0.2.0        
    ## [21] cli_3.6.5          mgcv_1.9-4         withr_3.0.3        magrittr_2.0.4    
    ## [25] digest_0.6.39      grid_4.5.3         rstudioapi_0.19.0  nlme_3.1-169      
    ## [29] lifecycle_1.0.5    vctrs_0.7.1        evaluate_1.0.5     glue_1.8.0        
    ## [33] farver_2.1.2       rmarkdown_2.31     purrr_1.2.2        tools_4.5.3       
    ## [37] pkgconfig_2.0.3    htmltools_0.5.9
