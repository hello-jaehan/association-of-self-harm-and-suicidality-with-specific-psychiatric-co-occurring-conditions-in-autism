library(meta)
library(ggplot2)
library(dmetar)

data_path <- "/Users/datapath"
data <- read.csv(data_path)
data$N <- as.numeric(data$N)
data$per_boy <- as.numeric(data$per_boy)
data$mean_age <- as.numeric(data$mean_age)
data$NOS <- as.numeric(data$NOS)
data$estimate <- as.numeric(data$estimate)
data$SE <- as.numeric(data$SE)

plots_dir <- "/Users/datapath"
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

meta_analysis_by_comorbidity <- function(data) {
  comorbidities <- unique(data$comorbid_uni)
  
  for (comorbidity in comorbidities) {
    subset_data <- subset(data, comorbid_uni == comorbidity)
    subset_data <- subset_data[!is.na(subset_data$estimate) & !is.na(subset_data$SE), ]
    
    if (nrow(subset_data) > 0) {
      print(paste("Processing comorbidity:", comorbidity))
      meta_res <- metagen(TE = log(subset_data$estimate), seTE = subset_data$SE, data = subset_data, studlab = subset_data$authorYear, sm = "OR", method.tau = "REML", hakn = TRUE, comb.random = TRUE, comb.fixed = FALSE)
      print(meta_res)
      
      if (nrow(subset_data) > 1) {
        print(">>>>>>>>>>>>>>>>>>Publication bias")
        metabias_res <- metabias(meta_res, k.min = 3)
        print(metabias_res)
        
        pcurve_result <- tryCatch({
          pcurve(meta_res)
        }, error = function(e) {
          message("Error in pcurve: ", e$message)
          NULL
        })
        
        if (!is.null(pcurve_result)) {
          print(pcurve_result)
        }
        
        print(">>>>>>>>>>>>Subgroup analysis")
        subgroup_res <- update.meta(meta_res, byvar = subset_data$adjusted, tau.common = TRUE)
        print(subgroup_res)

        png(filename = file.path(plots_dir, paste0("funnel_plot_", comorbidity, ".png")), width = 1000, height = 1200, res = 300)
        funnel(meta_res)
        dev.off()
        
        png(filename = file.path(plots_dir, paste0("forest_plot_", comorbidity, ".png")), width = 3000, height = 1000, res = 300)
        forest(meta_res)
        dev.off()
        
        png(filename = file.path(plots_dir, paste0("pcurve_plot_", comorbidity, ".png")), width = 2000, height = 1500, res = 300)
        pcurve_result <- tryCatch({
          pcurve(meta_res)
        }, error = function(e) {
          message("Error in pcurve: ", e$message)
          NULL
        })
        dev.off()
      }
    }
  }
}

meta_regression <- function(data, output_path) {
  comorbidities <- unique(data$comorbid_uni)
  
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
    write.csv(results, file = output_path, row.names = FALSE)
}

meta_analysis_by_comorbidity(data)

output_path <- "/Users/datapath/meta_regression_results.csv"
meta_regression(data, output_path)



