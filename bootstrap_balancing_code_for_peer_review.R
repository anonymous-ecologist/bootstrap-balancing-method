#Install required packages----

#List of packages
pkg_list <- c("ggplot2", "dplyr", "readxl", "smatr", "ggpmisc", "ggExtra", "grid", "gridExtra", "cowplot")

#Function to check and install missing packages
check_and_install <- function(pkg){
  if (!require(pkg, character.only = TRUE)){
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

#Apply the function to the list of packages
sapply(pkg_list, check_and_install)

#Load the required packages

library(ggplot2)
library(dplyr)
library(readxl)
library(smatr)
library(ggpmisc)
library(ggExtra)
library(grid)
library(gridExtra)
library(cowplot)

#Validating the proposed bootstrapping method (Section 2)----

#Import Tallo database

path_to_tallo<- "" #write path to file on machine

Tallo <- read_excel(path_to_tallo)

#Options

options(scipen = 1)

##Cleaning Tallo database----

#Remove any rows that are missing data for stem diameter or height

Tallo_complete <- subset(Tallo, !is.na(height_m) & !is.na(stem_diameter_cm))

#Remove height outliers, as specified in Tallo publication

Tallo_complete_no_ht_outliers<- Tallo_complete[which(Tallo_complete$height_outlier=='N'),]

##Fitting raw Tallo database----

#Make the "height" column numeric

Tallo_complete_no_ht_outliers$height_m<- as.numeric(Tallo_complete_no_ht_outliers$height_m) 

#Fit SMA model to full dataset

sma(data=Tallo_complete_no_ht_outliers, height_m~stem_diameter_cm, method = "SMA", log = "XY", alpha = 0.05)

##Create balanced dataset and fit with SMA model----

#Define the breaks for the diameter classes

breaks_log <- seq(log10(1), log10(max(Tallo_complete_no_ht_outliers$stem_diameter_cm)), by = log10(2))

breaks<- 10^breaks_log

#Classify diameters into logarithmic bins and summarize

diameter_summary <- Tallo_complete_no_ht_outliers %>%
  mutate(diameter_class = cut(stem_diameter_cm, breaks = breaks, right = FALSE, include.lowest = TRUE, labels = paste(breaks[-length(breaks)], "-", breaks[-1]))) %>%
  group_by(diameter_class) %>%
  summarise(count = n(), .groups = 'drop')

#Generate a new dataset which only contains stems <=128 cm DBH

Tallo_filtered <- subset(Tallo_complete_no_ht_outliers, Tallo_complete_no_ht_outliers$stem_diameter_cm < 128)

#Add diameter classes to filtered dataset

Tallo_filtered <- Tallo_filtered %>%
  mutate(diameter_class = cut(stem_diameter_cm, breaks = breaks, right = FALSE, include.lowest = TRUE, labels = paste(breaks[-length(breaks)], "-", breaks[-1])))

#Subset filtered dataset by randomly sampling 100 diameters for every 1 cm increment
set.seed(999)

Tallo_balanced <- Tallo_filtered %>%
  group_by(diameter_class) %>%
  slice_sample(n = 10000, replace = FALSE) %>%
  ungroup()

#Fit subset data

sma(data=Tallo_balanced, formula = height_m~stem_diameter_cm, method = "SMA", log = "XY", slope.test = 2/3)

##Create imbalanced dataset and fit with SMA model----

#Which diameter "class" has most observations?

Tallo_summary_filtered <- table(Tallo_filtered$diameter_class)

Tallo_summary_filtered_df<- as.data.frame(Tallo_summary_filtered)

max_obs<- max(Tallo_summary_filtered_df$Freq) #10 - 11 cm 

#Rename two columns in filtered summary

Tallo_summary_filtered_df = rename(Tallo_summary_filtered_df, diameter_class = Var1, freq = Freq)

#Add a new column representing the percent observations for each diameter class as a function of max observations

Tallo_summary_filtered_df<- Tallo_summary_filtered_df %>% 
  mutate(perc_max = (freq/max_obs)*100)

Tallo_summary_filtered_df$freq_dist_real<- (ceiling(Tallo_summary_filtered_df$perc_max))*100

#Remove rows which contain less than 10000 observations

Tallo_summary_filtered_df <- subset(Tallo_summary_filtered_df, freq >= 100)

#Randomly sample x number of points (x corresponding to freq_dist_real) of each diameter class

#First, merge subsetted and summary df's

Tallo_subset_merged <- merge(Tallo_subset, Tallo_summary_filtered_df, by = "diameter_class")

#Generate a vector of levels corresponding to diameter classes

diameter_classes<- as.vector(Tallo_summary_filtered_df$diameter_class)

#Make the diameter class column a factor with specific levels

Tallo_subset_merged$diameter_class<- factor(Tallo_subset_merged$diameter_class, levels = diameter_classes)

#Split the merged dataframe by stem diameter

Tallo_subset_split <- split(Tallo_subset_merged, Tallo_subset_merged$diameter_class)

#Generate a vector of random seeds

set.seed(123)  # Set a seed for reproducibility
seeds <- sample.int(n = length(Tallo_subset_split), size = length(Tallo_subset_split))

#Function to sample rows

sample_rows <- function(df, seed) {
  sample_size <- df$freq_dist_real[1]  # Get the frequency for this diameter
  set.seed(seed)  # Set the seed for this sampling iteration
  slice_sample(df, n = sample_size)  # Sample rows without replacement
}

#Apply the function to each dataframe in the list

sampled_dfs <- mapply(sample_rows, Tallo_subset_split, seeds, SIMPLIFY = FALSE)

#Combine all the sampled dataframes back into one dataframe

Tallo_imbalanced <- do.call(rbind, sampled_dfs)

#Fit SMA model to imbalanced dataset

sma(data=Tallo_imbalanced, formula = height_m~stem_diameter_cm, method = "SMA", log = "XY")

##Use bootstrap-balancing to balance imbalanced dataset and fit with SMA model----

#Determine the number of observations to bootstrap for each bin
bootstrap_summary <- Tallo_imbalanced %>%
  group_by(diameter_class) %>%
  summarize(original_count = n(),
            bootstrap_count = if_else(n() < 10000, 10000 - n(), 0)) %>%
  ungroup() 

#Use bootstrapping to balance imbalanced dataset
set.seed(999) #Set seed to "999" for reproducibility

Tallo_bootstrap_balanced <- Tallo_imbalanced %>%
  group_by(diameter_class) %>%
  do({
    # Get the bootstrap count for the current group
    current_diameter_class <- unique(.$diameter_class)
    current_bootstrap_count <- bootstrap_summary$bootstrap_count[bootstrap_summary$diameter_class == current_diameter_class]
    
    # Create a flag for the original data as FALSE
    original_data <- mutate(., Resampled = FALSE)
    
    # Bind the original data with the bootstrapped samples
    # and flag the bootstrapped data as TRUE
    bind_rows(
      original_data,
      mutate(slice_sample(., n = current_bootstrap_count, replace = TRUE), Resampled = TRUE)
    )
  }) %>%
  ungroup()

#Fit SMA model to bootstrapped dataset

sma(data=Tallo_bootstrap_balanced, height_m~stem_diameter_cm, log = "XY", method = "SMA", alpha = 0.05)

##Testing for common slope and elevation among three datasets----

#Generating three datasets

#Balanced dataset

Tallo_balanced_clean <- dplyr::select(Tallo_balanced, height_m, stem_diameter_cm)

Tallo_balanced_clean<- Tallo_balanced_clean %>% mutate(Dataset = "Balanced")

#Imbalanced dataset

Tallo_imbalanced_clean <- dplyr::select(Tallo_imbalanced, height_m, stem_diameter_cm)

Tallo_imbalanced_clean<- Tallo_imbalanced_clean %>% mutate(Dataset = "Unbalanced")

#Bootstrap balanced dataset

Tallo_bootstrap_balanced_clean<- dplyr::select(Tallo_bootstrap_balanced, height_m, stem_diameter_cm)

Tallo_bootstrap_balanced_clean<- Tallo_bootstrap_balanced_clean %>% mutate(Dataset = "Bootstrapped")

#Merge "clean" datasets for comparison

Tallo_datasets_clean<- rbind(Tallo_balanced_clean, Tallo_imbalanced_clean, Tallo_bootstrap_balanced_clean)

Tallo_datasets_clean$Dataset<- as.factor(Tallo_datasets_clean$Dataset)

Tallo_datasets_clean<- as.data.frame(Tallo_datasets_clean)

#Compare slopes 

Tallo_SMA_model_slope_test<- sma(height_m~stem_diameter_cm*Dataset, data=Tallo_datasets_clean, log = "XY", method = "SMA", multcomp = TRUE, multcompmethod = "adjusted")

summary(Tallo_SMA_model_slope_test)

#Compare Y-intercepts

Tallo_SMA_model_elevation_test<- sma(height_m~stem_diameter_cm+Dataset, data=Tallo_datasets_clean, log = "XY", method = "SMA", multcomp = TRUE, multcompmethod = "adjusted")

summary(Tallo_SMA_model_elevation_test)

##Figure 1----

#Plot balanced dataset

Tallo_balanced_subfigure<- ggplot(data=Tallo_balanced, aes(x=stem_diameter_cm, y=height_m))+
  geom_point(alpha=0.1, color="navy")+
  stat_ma_line(method = "SMA", color="darkorange", se=TRUE)+
  scale_x_log10()+
  scale_y_log10(breaks= c(3,10,30,100), limits = c(1,100))+
  theme_classic()+
  annotate("text", x = 1, y = 100, label = "b = 0.67 (95% CI: 0.67 – 0.68)", hjust = 0, vjust = 1)+
  labs(tag = 'a')+
  theme(aspect.ratio = 1,
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.tag = element_text(size = 14, face = "bold"))

Tallo_balanced_subfigure <- ggMarginal(Tallo_balanced_subfigure, type="histogram", binwidth= 0.3, color = "navy", fill="navy", alpha=0.1)

#Plot imbalanced dataset

Tallo_imbalanced_subfigure<- ggplot(data=Tallo_imbalanced, aes(x=stem_diameter_cm, y=height_m))+
  geom_point(alpha=0.1, color="navy")+
  stat_ma_line(method = "SMA", color="darkorange", se=TRUE)+
  scale_x_log10()+
  scale_y_log10(breaks= c(3,10,30,100), limits = c(1,100))+
  theme_classic()+
  annotate("text", x = 1, y = 100, label = "b = 0.78 (95% CI: 0.78 – 0.79)", hjust = 0, vjust = 1)+
  labs(tag = 'b')+
  theme(aspect.ratio = 1,
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.tag = element_text(size = 14, face = "bold"))

Tallo_imbalanced_subfigure<- ggMarginal(Tallo_imbalanced_subfigure, type="histogram", binwidth= 0.3, color = "navy", fill="navy", alpha=0.1)

#Plot bootstrap-balanced dataset

Tallo_bootstrap_balanced_subfigure<- ggplot(data=Tallo_bootstrap_balanced, aes(x=stem_diameter_cm, y=height_m))+
  geom_point(aes(colour = Resampled, alpha = Resampled))+
  scale_color_manual(values = c("navy", "cyan2"))+
  scale_alpha_manual(values = c(0.1, 0.03))+
  stat_ma_line(method = "SMA", color="darkorange", se=TRUE)+
  scale_x_log10()+
  scale_y_log10(breaks= c(3,10,30,100), limits = c(1,100))+
  theme_classic()+
  annotate("text", x = 1, y = 100, label = "b = 0.67 (95% CI: 0.67 – 0.68)", hjust = 0, vjust = 1)+
  labs(tag = 'c') +
  theme(aspect.ratio = 1,
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.tag = element_text(size = 14, face = "bold"),
        legend.position = "none")

Tallo_bootstrap_balanced_subfigure<- ggMarginal(Tallo_bootstrap_balanced_subfigure, type="histogram", binwidth= 0.3, groupColour = TRUE, groupFill = TRUE, alpha=0.1)

#Arrange subfigures into a single three panel figure and save

Fig_1_Ytitle <- textGrob("Plant height (m)", gp = gpar(fontsize = 14))

Fig_1_Xtitle <- textGrob("Stem diameter (cm)", gp = gpar(fontsize = 14))

Fig_1<- grid.arrange(Tallo_balanced_subfigure, Tallo_imbalanced_subfigure, Tallo_bootstrap_balanced_subfigure, nrow=1, 
                     left=textGrob("Plant height (m)",rot = 90, gp=gpar(fontsize=14)),
                     bottom=textGrob("Stem diameter (cm)", gp=gpar(fontsize=14)))

ggsave(filename = "Fig_1.jpg", plot=Fig_1, width = 30, height = 10, units = "cm", dpi = 600)

#Applying the proposed bootstrapping method (Section 3)----

#Import AnimalTraits database

path_to_animaltraits<- "" #write path to file on machine

AnimalTraits <- read_excel(path_to_animaltraits, 
                           col_types = c("text", "text", "text", 
                                         "text", "text", "text", "text", "text", 
                                         "numeric", "text", "text", "text", 
                                         "numeric", "text", "numeric", "numeric", 
                                         "numeric", "text", "numeric", "text", 
                                         "text", "numeric", "numeric", "numeric", 
                                         "numeric", "text", "numeric", 
                                         "numeric", "numeric", "numeric", 
                                         "numeric", "numeric", "numeric", 
                                         "numeric", "numeric", "text", "numeric", 
                                         "numeric", "text", "text", "text", 
                                         "numeric", "text"))

#Options

options(scipen = 1)

##Cleaning and filtering AnimalTraits database----

#First, retain only paired observations

AnimalTraits_complete <- subset(AnimalTraits, !is.na(metabolic_rate) & !is.na(body_mass))

#Generate dataset of just ectotherms

AnimalTraits_ectotherms<- AnimalTraits_complete %>%
  filter(class %in% c("Amphibia", "Arachnida", "Insecta", "Malacostraca", "Clitellata", "Gastropoda", "Reptilia", "Chilopoda"))

#Filter only basal metabolic rate (also called resting or standard metabolic rate)

AnimalTraits_complete_ecto_BMR<- AnimalTraits_ectotherms %>%
  filter(`metabolic rate - method` %in% c("basal metabolic rate", "resting metabolic rate", "standard metabolic rate"))

##Fit raw ectotherm BMR data with an SMA model----

#Fit an SMA model to BMR data

sma(metabolic_rate~body_mass, data=AnimalTraits_complete_ecto_BMR, method = "SMA", log = "XY", slope.test = 3/4)

##Balance the raw ectotherm BRM dataset and fit with SMA model----

# Define the breaks for the body mass classes

min_ecto<- min(AnimalTraits_complete_ecto_BMR$body_mass)

max_ecto<- max(AnimalTraits_complete_ecto_BMR$body_mass)

breaks_ecto<- c(1e-7, 1e-6, 1e-5, 1e-4,1e-3, 1e-2,1e-1, 1, 10)

# Determine maximum sample size across bins

summary <- AnimalTraits_complete_ecto_BMR %>%
  mutate(body_mass_class = cut(body_mass, breaks = breaks_ecto, right = FALSE, include.lowest = TRUE, labels = paste(breaks_ecto[-length(breaks_ecto)], "-", breaks_ecto[-1]))) %>%
  group_by(body_mass_class) %>%
  summarise(original_count = n())

max_obs<- max(summary$original_count)

# Classify body masses into logarithmic bins and determine bootstrap sample size

body_mass_summary_ecto_BMR <- AnimalTraits_complete_ecto_BMR %>%
  mutate(body_mass_class = cut(body_mass, breaks = breaks_ecto, right = FALSE, include.lowest = TRUE, labels = paste(breaks_ecto[-length(breaks_ecto)], "-", breaks_ecto[-1]))) %>%
  group_by(body_mass_class) %>%
  summarise(original_count = n(),
            bootstrap_count = if_else(n() < max_obs, max_obs - n(), 0)) %>%
  ungroup() 

#Add body mass classes to full dataset

AnimalTraits_complete_ecto_BMR <- AnimalTraits_complete_ecto_BMR %>%
  mutate(body_mass_class = cut(body_mass, breaks = breaks_ecto, right = FALSE, include.lowest = TRUE, labels = paste(breaks_ecto[-length(breaks_ecto)], "-", breaks_ecto[-1])))

#Generate sample dataset for plotting

set.seed(999) #Set seed to "999" for reproducibility

bootstrapped_data_sample_ecto_BMR <- AnimalTraits_complete_ecto_BMR %>%
  group_by(body_mass_class) %>%
  do({
    # Get the bootstrap count for the current group
    current_body_mass_class <- unique(.$body_mass_class)
    current_bootstrap_count <- body_mass_summary_ecto_BMR$bootstrap_count[body_mass_summary_ecto_BMR$body_mass_class == current_body_mass_class]
    
    # Create a flag for the original data as FALSE
    original_data <- mutate(., Resampled = FALSE)
    
    # Bind the original data with the bootstrapped samples
    # and flag the bootstrapped data as TRUE
    bind_rows(
      original_data,
      mutate(slice_sample(., n = current_bootstrap_count, replace = TRUE), Resampled = TRUE)
    )
  }) %>%
  ungroup()

#Fit model to bootstrapped dataset

sma(metabolic_rate~body_mass, data=bootstrapped_data_sample_ecto_BMR, method = "SMA", log = "XY", slope.test = 3/4)

##Figure 2----

log_format <- function(x) {
  parse(text = paste0("10^{", round(log10(x)), "}")) #create unique log format for displaying units
}

#Plot raw ectotherm BMR dataset

ecto_BMR_raw_subfigure<- ggplot(data=AnimalTraits_complete_ecto_BMR, aes(x=body_mass, y=metabolic_rate))+
  geom_point(alpha=0.5, color="navy")+
  stat_ma_line(method = "SMA", color="darkorange", se=TRUE)+
  scale_x_log10(labels = log_format)+
  scale_y_log10(labels = log_format)+
  theme_classic()+
  annotate("text", x = 10^-7, y = 1, label = "b = 0.80 (95% CI: 0.76 - 0.83)", hjust = 0, vjust = 1)+
  labs(tag = 'a') +
  theme(aspect.ratio = 1,
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.tag = element_text(size = 14, face = "bold"))

ecto_BMR_raw_subfigure<- ggMarginal(ecto_BMR_raw_subfigure, type="histogram", binwidth=1, color = "navy", fill="navy", alpha=0.1)

#Plot balanced ectotherm BMR data

ecto_BMR_bootstrap_balanced_subfigure<- ggplot(data=bootstrapped_data_sample_ecto_BMR, aes(x=body_mass, y=metabolic_rate))+
  geom_point(aes(colour = Resampled, alpha = Resampled))+
  scale_color_manual(values = c("navy", "cyan2"))+
  scale_alpha_manual(values = c(0.5, 0.3))+
  stat_ma_line(method = "SMA", color="darkorange", se=TRUE)+
  scale_x_log10(labels = log_format)+
  scale_y_log10(labels = log_format)+
  theme_classic()+
  annotate("text", x = 10^-7, y = 1, label = "b = 0.75 (95% CI: 0.73 - 0.77)", hjust = 0, vjust = 1)+
  labs(tag = 'b') +
  theme(aspect.ratio = 1,
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.tag = element_text(size = 14, face = "bold"),
        legend.position = "none")

ecto_BMR_bootstrap_balanced_subfigure<- ggMarginal(ecto_BMR_bootstrap_balanced_subfigure, type="histogram", binwidth=1, groupColour = TRUE, groupFill = TRUE, alpha=0.1)

#Arrange subfigures into a single two panel figure and save

Fig_2<- grid.arrange(ecto_BMR_raw_subfigure, ecto_BMR_bootstrap_balanced_subfigure, nrow=1, 
                     left=textGrob("Basal metabolic rate (W)",rot = 90, gp=gpar(fontsize=14)),
                     bottom=textGrob("Body mass (kg)", gp=gpar(fontsize=14)))

ggsave(filename = "Fig_2.jpg", plot=Fig_2, width = 20, height = 10, units = "cm", dpi=600)

#Common data processing methods reduce statistical power and widen confidence intervals (Section 4)----

##Import and clean dataset----

path_to_conduit_dataset<- "" #write path to file on machine

df_full <- read_excel(path_to_conduit_dataset)

#Subset aboveground data for full dataset

df_full_aboveground <- df_full[which(df_full$Organ=='Leaf'|df_full$Organ=='Twig'| df_full$Organ=='Branch'| df_full$Organ=='Stem'), ]
df_full_aboveground$Organ<- factor(df_full_aboveground$Organ, levels=(c("Leaf", "Twig", "Branch", "Stem"))) 

##Simulating the effects of data preprocessing decisions----

###Bootstrapping data at different sample sizes, averaging the data, and fitting a model----

## Function to perform random sampling, fit SMA model, and extract coefficients

run_iteration <- function(df_full_aboveground, sample_size) {
  # Random sampling of D across each and every L
  bootstrapped_data <- df_full_aboveground %>%
    group_by(L) %>%
    slice_sample(n = sample_size, replace = TRUE) %>%
    ungroup()
  
  # Calculate the average DAVG for each L
  averaged_data <- bootstrapped_data %>%
    group_by(L) %>%
    summarize(DAVG = mean(DAVG))
  
  # Fit a SMA model to averaged_data
  bootstrapped_scaling <- sma(data = averaged_data, DAVG ~ L, log = "XY", method = "SMA")
  
  # Extract coefficients
  bootstrapped_scaling_summary <- bootstrapped_scaling$groupsummary
  bootstrapped_results <- c(
    slope = bootstrapped_scaling_summary$Slope,
    y_intercept = bootstrapped_scaling_summary$Int,
    slope_lowCI = bootstrapped_scaling_summary$Slope_lowCI,
    slope_highCI = bootstrapped_scaling_summary$Slope_highCI
  )
  
  # Return a list of coefficients
  return(bootstrapped_results)
}

## Specify the sample sizes and number of iterations
sample_sizes <- c(1, 5, 10, 25, 50, 75, 100,1000)
num_iterations <- 1000  

## Create an empty dataframe to store the results
bootstrapped_results_df <- data.frame(
  iteration = integer(),
  mean_slope = double(),
  mean_y_intercept = double(),
  mean_slope_lowCI = double(),
  mean_slope_highCI = double(),
  true_lowCI = double(),
  true_highCI = double(),
  se_slope = double(),
  se_y_intercept = double(),
  se_slope_lowCI = double(),
  se_slope_highCI = double(),
  cv_slope= double(),
  cv_slope_lowCI = double(),
  cv_slope_highCI = double(),
  sample_size = integer()
)

## Loop through different sample sizes
for (sample_size in sample_sizes) {
  for (num_iteration in num_iterations) {
    slopes <- numeric(num_iteration)
    y_intercepts <- numeric(num_iteration)
    slope_lowCIs <- numeric(num_iteration)
    slope_highCIs <- numeric(num_iteration)
    
    # Run multiple iterations
    for (iteration in 1:num_iteration) {
      # Set seed for this iteration based on the sample size and iteration
      set.seed(iteration + sample_size)
      
      # Run the iteration and extract the coefficients
      coefficients_iteration <- run_iteration(df_full_aboveground, sample_size)
      
      # Store the results
      slopes[iteration] <- coefficients_iteration["slope"]
      y_intercepts[iteration] <- coefficients_iteration["y_intercept"]
      slope_lowCIs[iteration] <- coefficients_iteration["slope_lowCI"]
      slope_highCIs[iteration] <- coefficients_iteration["slope_highCI"]
    }
    
    # Calculate mean and standard error for each coefficient
    mean_slope <- mean(slopes)
    se_slope <- sd(slopes) / sqrt(num_iteration)
    
    mean_y_intercept <- mean(y_intercepts)
    se_y_intercept <- sd(y_intercepts) / sqrt(num_iteration)
    
    mean_slope_lowCI <- mean(slope_lowCIs)
    se_slope_lowCI <- sd(slope_lowCIs) / sqrt(num_iteration)
    
    mean_slope_highCI <- mean(slope_highCIs)
    se_slope_highCI <- sd(slope_highCIs) / sqrt(num_iteration)
    
    true_lowCI<- mean_slope - 1.96 * ((sd(slopes))/(sqrt(sample_size*169)))
    true_highCI<- mean_slope + 1.96 * ((sd(slopes))/(sqrt(sample_size*169)))
    
    cv_slope<- sd(slopes) / mean(slopes)
    cv_slope_lowCI<- sd(slope_lowCIs) / mean(slope_lowCIs)
    cv_slope_highCI<- sd(slope_highCIs) / mean(slope_highCIs)
    
    
    
    # Store the results in the dataframe
    bootstrapped_results_df <- rbind(bootstrapped_results_df, c(iteration, mean_slope, mean_y_intercept, mean_slope_lowCI, mean_slope_highCI, se_slope, se_y_intercept, se_slope_lowCI, se_slope_highCI, true_lowCI, true_highCI, cv_slope, cv_slope_lowCI, cv_slope_highCI, sample_size))
  }
}

## Rename the columns in bootstrapped_results_df
colnames(bootstrapped_results_df) <- c("iteration", "mean_slope", "mean_y_intercept", "mean_slope_lowCI", "mean_slope_highCI", "se_slope", "se_y_intercept", "se_slope_lowCI", "se_slope_highCI", "true_lowCI", "true_highCI", "cv_slope", "cv_slope_lowCI", "cv_slope_highCI",  "sample_size")

##Store the results

bootstrap_simulation_arith_mean_summary<- bootstrapped_results_df

###Bootstrapping data at different sample sizes, averaging the data and weighing by conductance (Dh), and fitting a model----

## Function to perform random sampling, fit SMA model, and extract coefficients

run_iteration <- function(df_full_aboveground, sample_size) {
  # Random sampling of D across each and every L
  bootstrapped_data <- df_full_aboveground %>%
    group_by(L) %>%
    slice_sample(n = sample_size, replace = TRUE) %>%
    ungroup()
  
  # Calculate the average DAVG for each L
  averaged_data <- bootstrapped_data %>%
    group_by(L) %>%
    summarize(DhAVG = sum(DAVG^5)/sum(DAVG^4))
  
  # Fit a SMA model to averaged_data
  bootstrapped_scaling <- sma(data = averaged_data, DhAVG ~ L, log = "XY", method = "SMA")
  
  # Extract coefficients
  bootstrapped_scaling_summary <- bootstrapped_scaling$groupsummary
  bootstrapped_results <- c(
    slope = bootstrapped_scaling_summary$Slope,
    y_intercept = bootstrapped_scaling_summary$Int,
    slope_lowCI = bootstrapped_scaling_summary$Slope_lowCI,
    slope_highCI = bootstrapped_scaling_summary$Slope_highCI
  )
  
  # Return a list of coefficients
  return(bootstrapped_results)
}

## Specify the sample sizes and number of iterations
sample_sizes <- c(1, 5, 10, 25, 50, 75, 100,1000)
num_iterations <- 1000  

## Create an empty dataframe to store the results
bootstrapped_results_df <- data.frame(
  iteration = integer(),
  mean_slope = double(),
  mean_y_intercept = double(),
  mean_slope_lowCI = double(),
  mean_slope_highCI = double(),
  se_slope = double(),
  se_y_intercept = double(),
  se_slope_lowCI = double(),
  se_slope_highCI = double(),
  cv_slope= double(),
  cv_slope_lowCI = double(),
  cv_slope_highCI = double(),
  sample_size = integer()
)

## Loop through different sample sizes
for (sample_size in sample_sizes) {
  for (num_iteration in num_iterations) {
    slopes <- numeric(num_iteration)
    y_intercepts <- numeric(num_iteration)
    slope_lowCIs <- numeric(num_iteration)
    slope_highCIs <- numeric(num_iteration)
    
    # Run multiple iterations
    for (iteration in 1:num_iteration) {
      # Set seed for this iteration based on the sample size and iteration
      set.seed(iteration + sample_size)
      
      # Run the iteration and extract the coefficients
      coefficients_iteration <- run_iteration(df_full_aboveground, sample_size)
      
      # Store the results
      slopes[iteration] <- coefficients_iteration["slope"]
      y_intercepts[iteration] <- coefficients_iteration["y_intercept"]
      slope_lowCIs[iteration] <- coefficients_iteration["slope_lowCI"]
      slope_highCIs[iteration] <- coefficients_iteration["slope_highCI"]
    }
    
    # Calculate mean and standard error for each coefficient
    mean_slope <- mean(slopes)
    se_slope <- sd(slopes) / sqrt(num_iteration)
    
    mean_y_intercept <- mean(y_intercepts)
    se_y_intercept <- sd(y_intercepts) / sqrt(num_iteration)
    
    mean_slope_lowCI <- mean(slope_lowCIs)
    se_slope_lowCI <- sd(slope_lowCIs) / sqrt(num_iteration)
    
    mean_slope_highCI <- mean(slope_highCIs)
    se_slope_highCI <- sd(slope_highCIs) / sqrt(num_iteration)
    
    cv_slope<- sd(slopes) / mean(slopes)
    cv_slope_lowCI<- sd(slope_lowCIs) / mean(slope_lowCIs)
    cv_slope_highCI<- sd(slope_highCIs) / mean(slope_highCIs)
    
    
    
    # Store the results in the dataframe
    bootstrapped_results_df <- rbind(bootstrapped_results_df, c(iteration, mean_slope, mean_y_intercept, mean_slope_lowCI, mean_slope_highCI, se_slope, se_y_intercept, se_slope_lowCI, se_slope_highCI, cv_slope, cv_slope_lowCI, cv_slope_highCI, sample_size))
  }
}

## Rename the columns in bootstrapped_results_df
colnames(bootstrapped_results_df) <- c("iteration", "mean_slope", "mean_y_intercept", "mean_slope_lowCI", "mean_slope_highCI", "se_slope", "se_y_intercept", "se_slope_lowCI", "se_slope_highCI","cv_slope", "cv_slope_lowCI", "cv_slope_highCI",  "sample_size")

##Store the results

bootstrap_simulation_hydr_mean_summary<- bootstrapped_results_df

###Bootstrapping data at different sample sizes and fitting a model without averaging----

## Function to perform random sampling, fit SMA model, and extract coefficients
run_iteration <- function(df_full_aboveground, sample_size) {
  # Random sampling of D across each and every L
  bootstrapped_data <- df_full_aboveground %>%
    group_by(L) %>%
    slice_sample(n = sample_size, replace = TRUE) %>%
    ungroup()
  
  # Fit a SMA model to bootstrapped_data
  bootstrapped_scaling <- sma(data = bootstrapped_data, DAVG ~ L, log = "XY", method = "SMA")
  
  # Extract coefficients
  bootstrapped_scaling_summary <- bootstrapped_scaling$groupsummary
  bootstrapped_results <- c(
    slope = bootstrapped_scaling_summary$Slope,
    y_intercept = bootstrapped_scaling_summary$Int,
    slope_lowCI = bootstrapped_scaling_summary$Slope_lowCI,
    slope_highCI = bootstrapped_scaling_summary$Slope_highCI
  )
  
  # Return a list of coefficients
  return(bootstrapped_results)
}

## Specify the sample sizes and number of iterations
sample_sizes <- c(1, 5, 10, 25, 50, 75, 100,1000)
num_iterations <- 1000  # You can adjust the number of iterations

## Create an empty dataframe to store the results
bootstrapped_results_df <- data.frame(
  iteration = integer(),
  mean_slope = double(),
  mean_y_intercept = double(),
  mean_slope_lowCI = double(),
  mean_slope_highCI = double(),
  se_slope = double(),
  se_y_intercept = double(),
  se_slope_lowCI = double(),
  se_slope_highCI = double(),
  sample_size = integer()
)

## Loop through different sample sizes
for (sample_size in sample_sizes) {
  for (num_iteration in num_iterations) {
    slopes <- numeric(num_iteration)
    y_intercepts <- numeric(num_iteration)
    slope_lowCIs <- numeric(num_iteration)
    slope_highCIs <- numeric(num_iteration)
    
    # Run multiple iterations
    for (iteration in 1:num_iteration) {
      # Set seed for this iteration based on the sample size and iteration
      set.seed(iteration + sample_size)
      
      # Run the iteration and extract the coefficients
      coefficients_iteration <- run_iteration(df_full_aboveground, sample_size)
      
      # Store the results
      slopes[iteration] <- coefficients_iteration["slope"]
      y_intercepts[iteration] <- coefficients_iteration["y_intercept"]
      slope_lowCIs[iteration] <- coefficients_iteration["slope_lowCI"]
      slope_highCIs[iteration] <- coefficients_iteration["slope_highCI"]
    }
    
    # Calculate mean and standard error for each coefficient
    mean_slope <- mean(slopes)
    se_slope <- sd(slopes) / sqrt(num_iteration)
    
    mean_y_intercept <- mean(y_intercepts)
    se_y_intercept <- sd(y_intercepts) / sqrt(num_iteration)
    
    mean_slope_lowCI <- mean(slope_lowCIs)
    se_slope_lowCI <- sd(slope_lowCIs) / sqrt(num_iteration)
    
    mean_slope_highCI <- mean(slope_highCIs)
    se_slope_highCI <- sd(slope_highCIs) / sqrt(num_iteration)
    
    cv_slope<- sd(slopes) / mean(slopes)
    cv_slope_lowCI<- sd(slope_lowCIs) / mean(slope_lowCIs)
    cv_slope_highCI<- sd(slope_highCIs) / mean(slope_highCIs)
    
    # Store the results in the dataframe
    bootstrapped_results_df <- rbind(bootstrapped_results_df, c(iteration, mean_slope, mean_y_intercept, mean_slope_lowCI, mean_slope_highCI, se_slope, se_y_intercept, se_slope_lowCI, se_slope_highCI, cv_slope, cv_slope_lowCI, cv_slope_highCI, sample_size))
  }
}

## Rename the columns in bootstrapped_results_df
colnames(bootstrapped_results_df) <- c("iteration", "mean_slope", "mean_y_intercept", "mean_slope_lowCI", "mean_slope_highCI", "se_slope", "se_y_intercept", "se_slope_lowCI", "se_slope_highCI", "cv_slope", "cv_slope_lowCI", "cv_slope_highCI", "sample_size")

## Store the summary in separate dataframe

bootstrap_simulation_raw_summary<- bootstrapped_results_df

##Figure 3----

##Make fake plot to extract legend from

Model_data<- c("Mean conduit diameter (arithmetic)", "Mean conduit diameter (hydraulic)", "Raw data")
Value<- c(1,2,3)

df_data<- data.frame(Model_data, Value)

df_data$Model_data<- factor(df_data$Model_data, levels = c("Mean conduit diameter (arithmetic)", "Mean conduit diameter (hydraulic)", "Raw data"))

data_grp<- levels(df_data$Model_data)
data_color<- c("red", "yellow2", "blue3")
names(data_color) <- data_grp

data_fake_plot<- ggplot(data=df_data, aes(x=Model_data, y=Value))+
  geom_bar(aes(fill= Model_data), stat = "identity")+
  scale_fill_manual(values = data_color)+
  labs(fill = "Data preprocessing")+
  theme(legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

data_fake_plot

data_fake_plot_legend<- cowplot::get_legend(data_fake_plot)

##Build actual plot

bootstrap_simulation_V2_plot<- ggplot()+
  geom_path(data=bootstrap_simulation_arith_mean_summary, aes(y=mean_slope, x=sample_size), color="red", linewidth=1.5, linetype = "solid")+
  geom_ribbon(data=bootstrap_simulation_arith_mean_summary, aes(x=sample_size, ymin=mean_slope_lowCI, ymax=mean_slope), fill="red", alpha=0.15)+
  geom_ribbon(data=bootstrap_simulation_arith_mean_summary, aes(x=sample_size, ymin=mean_slope, ymax=mean_slope_highCI), fill="red", alpha=0.15)+
  geom_path(data=bootstrap_simulation_hydr_mean_summary, aes(y=mean_slope, x=sample_size), color="yellow2", linewidth=1.5, linetype = "solid")+
  geom_ribbon(data=bootstrap_simulation_hydr_mean_summary, aes(x=sample_size, ymin=mean_slope_lowCI, ymax=mean_slope), fill="yellow2", alpha=0.15)+
  geom_ribbon(data=bootstrap_simulation_hydr_mean_summary, aes(x=sample_size, ymin=mean_slope, ymax=mean_slope_highCI), fill="yellow2", alpha=0.15)+
  geom_path(data=bootstrap_simulation_raw_summary, aes(y=mean_slope, x=sample_size), color="blue3", linewidth=1.5, linetype = "solid")+
  geom_ribbon(data=bootstrap_simulation_raw_summary, aes(x=sample_size, ymin=mean_slope_lowCI, ymax=mean_slope), fill="blue3", alpha=0.15)+
  geom_ribbon(data=bootstrap_simulation_raw_summary, aes(x=sample_size, ymin=mean_slope, ymax=mean_slope_highCI), fill="blue3", alpha=0.15)+
  xlab("Number of conduits measured per sample (dimensionless)")+
  ylab(expression(paste("Fitted scaling exponent, ", italic("b"), " (dimensionless)")))+
  scale_x_log10(breaks=c(1, 5, 10, 25, 50, 100,1000))+
  theme_classic()+
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

bootstrap_simulation_V2_plot

bootstrap_simulation_V2_plot_combined<- ggdraw() +
  draw_plot(bootstrap_simulation_V2_plot, x = 0, y = 0, width = 1, height = 1)+
  draw_grob(data_fake_plot_legend, x = 0.6, y = 0.85, width = 0.1, height = 0.1)

bootstrap_simulation_V2_plot_combined

#Save plot!

ggsave("Fig_3.jpg", plot = bootstrap_simulation_V2_plot_combined, width=6, height = 6, dpi = 600)