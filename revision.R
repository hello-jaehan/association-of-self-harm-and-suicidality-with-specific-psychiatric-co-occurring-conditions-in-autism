# Load required libraries
library(meta)
library(ggplot2)
library(dmetar)

# Read the data
data_path <- "/Users/jayhaankim/Desktop/Research/1_First_Ongoing/11_ASD self-harm by comorbidity/2_Data extraction/data_for_R_revision/suicidality.csv"
data <- read.csv(data_path)
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
        }
        }
      }
    }
  }
}

  # Save results to CSV file
  write.csv(results, file = output_path, row.names = FALSE)
}

# Perform meta-analysis by comorbidity and subgroup analysis by adjusted
meta_analysis_by_comorbidity(data)
