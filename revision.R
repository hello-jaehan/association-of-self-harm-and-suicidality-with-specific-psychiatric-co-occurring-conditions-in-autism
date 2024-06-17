# Load required libraries
library(meta)
library(ggplot2)
library(dmetar)

# Read the data
data_path <- "/Users/jayhaankim/Desktop/Research/1_First_Ongoing/11_ASD self-harm by comorbidity/2_Data extraction/data_for_R_revision/suicidality.csv"
data <- read.csv(data_path)
data$N <- as.numeric(data$N)
data$per_boy <- as.numeric(data$per_boy)
data$mean_age <- as.numeric(data$mean_age)
data$NOS <- as.numeric(data$NOS)
data$estimate <- as.numeric(data$estimate)
data$SE <- as.numeric(data$SE)

# Ensure the plots directory exists
plots_dir <- "/Users/jayhaankim/Desktop/Research/1_First_Ongoing/11_ASD self-harm by comorbidity/2_Data extraction/data_for_R_revision/plots"
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

# Function to perform meta-analysis for each comorbidity
meta_analysis_by_comorbidity <- function(data) {
  comorbidities <- unique(data$comorbid_uni)
  
  for (comorbidity in comorbidities) {
    subset_data <- subset(data, comorbid_uni == comorbidity)
    
    # Remove rows with missing estimate or SE
    subset_data <- subset_data[!is.na(subset_data$estimate) & !is.na(subset_data$SE), ]
    
    if (nrow(subset_data) > 0) {
      # Check the data before meta-analysis
      print(paste("Processing comorbidity:", comorbidity))
      
      # Perform random-effects meta-analysis
      meta_res <- metagen(TE = log(subset_data$estimate), seTE = subset_data$SE, data = subset_data, studlab = subset_data$authorYear, sm = "OR", method.tau = "REML", hakn = TRUE, comb.random = TRUE, comb.fixed = FALSE)
      print(meta_res)
      
      print(">>>>>>>>>>>>>>>>>>Subgroup analysis")
      print(update(meta_res, subgroup = adjusted, tau.common = TRUE))
      
      # Save forest plot
      forest_plot_path <- file.path(plots_dir, paste0("forest_plot_", comorbidity, ".png"))
      png(forest_plot_path, width = 800, height = 600)
      forest(meta_res)
      dev.off()
      
      if (nrow(subset_data) > 2) {
        # Assess publication bias
        print(">>>>>>>>>>>>>>>>>>Publication bias")
        metabias_res <- metabias(meta_res, k.min = 3)
        print(metabias_res)
        
        # Save funnel plot
        funnel_plot_path <- file.path(plots_dir, paste0("funnel_plot_", comorbidity, ".png"))
        png(funnel_plot_path, width = 800, height = 600)
        funnel(meta_res)
        dev.off()
        
        # Trim-and-Fill method if Egger's test p-value < 0.05
        if (metabias_res$p.value < 0.05) {
          print("Performing Trim-and-Fill method due to significant Egger's test")
          trimfill_res <- trimfill(meta_res)
          print(trimfill_res)
          
          # Save trim-and-fill funnel plot
          trimfill_funnel_plot_path <- file.path(plots_dir, paste0("trimfill_funnel_plot_", comorbidity, ".png"))
          png(trimfill_funnel_plot_path, width = 800, height = 600)
          funnel(trimfill_res)
          dev.off()
        }
        
        # Pcurve analysis
        pcurve_result <- tryCatch({
          pcurve(meta_res)
        }, error = function(e) {
          message("Error in pcurve: ", e$message)
          NULL
        })
        if (!is.null(pcurve_result)) {
          print(pcurve_result)
        }
      }
    }
  }
}

# Function to perform meta-regression and save results as a CSV file
meta_regression <- function(data, output_path) {
  comorbidities <- unique(data$comorbid_uni)
  
  # Initialize an empty data frame to store results
  results <- data.frame(
    Comorbidity = character(),
    Moderator = character(),
    k = integer(),
    Coefficient = character(),
    P = numeric(),
    NA_reason = character(),
    stringsAsFactors = FALSE
  )
  
  for (comorbidity in comorbidities) {
    subset_data <- subset(data, comorbid_uni == comorbidity)
    
    for (i in c("mean_age", "per_boy", "N", "NOS", "adjusted")) {
      subset_data_i <- subset_data[!is.na(subset_data[[i]]), ]  # Remove rows with NA in the specific covariate
      k <- nrow(subset_data_i)
      
      if (k > 3) {
        meta_res <- metagen(TE = log(estimate), seTE = SE, data = subset_data_i, sm = "OR", method.tau = "REML")
        if (k > 1) {  # Ensure at least two observations
          tryCatch({
            reg_res <- metareg(meta_res, as.formula(paste("~", i)))
            coef <- coef(summary(reg_res))
            estimate <- round(coef[2, "estimate"], 4)
            ci_lb <- round(coef[2, "ci.lb"], 4)
            ci_ub <- round(coef[2, "ci.ub"], 4)
            pval <- round(coef[2, "pval"], 4)
            results <- rbind(results, data.frame(
              Comorbidity = comorbidity,
              Moderator = i,
              k = k,
              Coefficient = paste0(estimate, " (", ci_lb, " to ", ci_ub, ")"),
              P = pval,
              NA_reason = ""
            ))
          }, error = function(e) {
            results <- rbind(results, data.frame(
              Comorbidity = comorbidity,
              Moderator = i,
              k = k,
              Coefficient = NA,
              P = NA,
              NA_reason = paste("Error:", e$message)
            ))
          })
        } else {
          results <- rbind(results, data.frame(
            Comorbidity = comorbidity,
            Moderator = i,
            k = k,
            Coefficient = NA,
            P = NA,
            NA_reason = "k<2"
          ))
        }
      } else {
        results <- rbind(results, data.frame(
          Comorbidity = comorbidity,
          Moderator = i,
          k = k,
          Coefficient = NA,
          P = NA,
          NA_reason = "k<4"
        ))
      }
    }
  }
  
  # Save results to CSV file
  write.csv(results, file = output_path, row.names = FALSE)
}

# Perform meta-analysis by comorbidity and subgroup analysis by adjusted
meta_analysis_by_comorbidity(data)

# Perform meta-regression
output_path <- "/Users/jayhaankim/Downloads/meta_regression_results.csv"
meta_regression(data, output_path)



