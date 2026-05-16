library(MuMIn)
library(viridis)

setup_tumour <- function(border_growth_model, growth_rate, death_rate=0.01)
{
  width <- 7000
  height <- 7000
  
  sim <- TissueSimulation("Test", width, height)
  
  sim$border_growth_model <- border_growth_model
  
  # Clone 1
  sim$add_mutant(name = "Clone 1", growth_rates = growth_rate,
                 death_rates = death_rate)
  
  # place the first tumour cell
  sim$place_cell("Clone 1", width/2, height/2)
  
  return(sim)
}

grow_tumour_up_to <- function(sim, n)
{
  # let the tumour evolve until the clone 1 consists of 30K cells
  sim$run_up_to_size("Clone 1", n)
  
  return(sim)
}

compare_complexity <- function(data, x_col, y_col, title=NULL) {
  clean_df <- data[data[[x_col]] > 0, ]
  clean_df <- clean_df[order(clean_df[[x_col]]), ]
  
  x <- as.numeric(clean_df[[x_col]])
  y <- as.numeric(clean_df[[y_col]])
  
  mod_const   <- lm(y ~ 1)
  mod_lin   <- lm(y ~ I(x))
  mod_logn <- lm(y ~ I(log(x)))
  mod_inv <- lm(y ~ 1/I(x))
  #mod_nlogn <- lm(y ~ x:log(x)+x)
  #mod_nnlogn <- lm(y ~ x^2:log(x)+x)
  #mod_nsqrtnlogn <- lm(y ~ x:sqrt(x):log(x)+x)
  #mod_nsqrtn <- lm(y ~ x:sqrt(x)+x)
  #mod_quad  <- lm(y ~ poly(x, 2))
  
  mod_nlogn <- lm(y ~ I(x*log(x)))
  mod_nnlogn <- lm(y ~ I(x^2*log(x)))
  mod_nsqrtnlogn <- lm(y ~ I(x*sqrt(x)*log(x)))
  mod_nsqrtnlogn_nlogn <- lm(y ~ I(x*sqrt(x)*log(x))+I(x*log(x)))
  mod_nsqrtn_nlogn <- lm(y ~ I(x*sqrt(x))+I(x*log(x)))
  mod_nsqrtn <- lm(y ~ I(x*sqrt(x)))
  mod_quad  <- lm(y ~ I(x^2))
  
  fit <- try(lm(log(y + 0.0001) ~ x), silent = TRUE)
  if(inherits(fit, "lm")) {
    start_values <- list(a = exp(coef(fit)[1]), b = coef(fit)[2])
  } else {
    start_values <- list(a = 0.01, b = 0.01) # Sensible defaults for time
  }
  
  exp_AIC <- NA
  exp_AICc <- NA

  tryCatch({
    mod_exp <- nls(y ~ a * exp(b * x),
                   start = start_values,
                   control = nls.control(minFactor = 1/8192, maxiter = 200))
    exp_AIC <- AIC(mod_exp)
    exp_AICc <- AICc(mod_exp)
  }, error = function(e) {})
  
  # 3. Create Comparison Table
  results <- data.frame(
    Model = c("O(1)", "O(1/n)", "O(n)", "O(n log(n))",
              "O(n^{3/2})", "O(n^{3/2}+nlog(n))", "O(n^{3/2} log(n))",
              "O(n^{3/2} log(n)+nlog(n))", "O(n^2)", "O(n^2 log(n))"),
    AIC = c(AIC(mod_const), AIC(mod_inv), AIC(mod_lin), AIC(mod_nlogn),
            AIC(mod_nsqrtn), AIC(mod_nsqrtn_nlogn), AIC(mod_nsqrtnlogn),
            AIC(mod_nsqrtnlogn_nlogn),
            AIC(mod_quad), AIC(mod_nnlogn)),
    AICc = c(AICc(mod_const), AIC(mod_inv), AICc(mod_lin), AICc(mod_nlogn),
             AICc(mod_nsqrtn), AICc(mod_nsqrtn_nlogn),  AICc(mod_nsqrtnlogn),
             AICc(mod_nsqrtnlogn_nlogn),
             AICc(mod_quad), AICc(mod_nnlogn)),
    R_Squared = c(summary(mod_const)$r.squared,
                  summary(mod_inv)$r.squared,
                  summary(mod_lin)$r.squared,
                  summary(mod_nlogn)$r.squared,
                  summary(mod_nsqrtn)$r.squared,
                  summary(mod_nsqrtn_nlogn)$r.squared,
                  summary(mod_nsqrtnlogn)$r.squared,
                  summary(mod_nsqrtnlogn_nlogn)$r.squared,
                  summary(mod_quad)$r.squared,
                  summary(mod_nnlogn)$r.squared)
  )

  if (!is.null(title)) {
    cat(paste(title, "\n"))
    print(results)
  }
  
  return(results)
}

latex_AICc_row <- function(AICcs, AICc_delta_threshold, first_columns=c()) {
  AICcs <- round(AICcs, digits = 2)
  
  min_AICc <- min(AICcs, na.rm = TRUE)
  
  cat(paste(first_columns, collapse = "&"))
  
  for (AICc in AICcs) {
    text_AICc <- paste0("\\num{", AICc, "}")
    text_AICc <- AICc

    if (!is.na(AICc)) {
      AICc_delta <- 100*(AICc-min_AICc)
      
      if (AICc_delta/abs(min_AICc)<=AICc_delta_threshold) {
        text_AICc <- paste0("\\textbf{", text_AICc,"}")
        
        if (AICc_delta==0) {
          text_AICc <- paste0("\\underline{", text_AICc, "}")
        }
      }
    }
    cat(paste0("&", text_AICc))
  }
  cat("\\\\\n")
}

compare_complexity <- function(data, x_col, y_col, title=NULL) {
  clean_df <- data[data[[x_col]] > 0, ]
  clean_df <- clean_df[order(clean_df[[x_col]]), ]
  
  x <- as.numeric(clean_df[[x_col]])
  y <- as.numeric(clean_df[[y_col]])
  
  mod_const   <- lm(y ~ 1)
  mod_inv <- lm(y ~ I(1/x))
  mod_lin   <- lm(y ~ I(x))
  mod_logn <- lm(y ~ I(log(x))+x)
  
  mod_nlogn <- lm(y ~ I(x*log(x)))
  mod_nsqrtn <- lm(y ~ I(x*sqrt(x)))
  mod_quad  <- lm(y ~ I(x^2))
  mod_cube  <- lm(y ~ I(x^3))
  
  fit <- try(lm(log(y + 0.0001) ~ x), silent = TRUE)
  if(inherits(fit, "lm")) {
    start_values <- list(a = exp(coef(fit)[1]), b = coef(fit)[2])
  } else {
    start_values <- list(a = 0.01, b = 0.01) # Sensible defaults for time
  }
  
  exp_AIC <- NA
  exp_AICc <- NA
  
  tryCatch({
    mod_exp <- nls(y ~ a * exp(b * x),
                   start = start_values,
                   control = nls.control(minFactor = 1/8192, maxiter = 200))
    exp_AIC <- AIC(mod_exp)
    exp_AICc <- AICc(mod_exp)
  }, error = function(e) {})
  
  # 3. Create Comparison Table
  results <- data.frame(
    Model = c("O(1/n)", "O(n)", "O(n log(n))",
              "O(n^{3/2})", "O(n^2)", #"O(n^3)",
              "2^n"),
    AIC = c(AIC(mod_inv), AIC(mod_lin), AIC(mod_nlogn),
            AIC(mod_nsqrtn), AIC(mod_quad), #AIC(mod_cube),
            exp_AIC),
    AICc = c(AICc(mod_inv), AICc(mod_lin), AICc(mod_nlogn),
             AICc(mod_nsqrtn), AICc(mod_quad), #AICc(mod_cube),
             exp_AICc),
    R_Squared = c(summary(mod_inv)$r.squared,
                  summary(mod_lin)$r.squared,
                  summary(mod_nlogn)$r.squared,
                  summary(mod_nsqrtn)$r.squared,
                  summary(mod_quad)$r.squared,
                  #summary(mod_cube)$r.squared,
                  NA)
  )
  
  if (!is.null(title)) {
    cat(paste(title, "\n"))
    print(results)
  }
  
  return(results)
}

get_sample_forest_filename <- function(i, k, j, border_growth_v, sample_on_edge,
                                       duplication_rate, death_rate)
{
  forest_filename <-paste0("sample_forest_", i, "_",
                           k, "_", j, "_")
  if (border_growth_v) {
    forest_filename <-paste0(forest_filename, "border")
  } else {
    
    forest_filename <-paste0(forest_filename, "uniform")
  }
  
  if (sample_on_edge) {
    forest_filename <-paste0(forest_filename, "_edge")
  } else {
    forest_filename <- paste0(forest_filename, "_center")
  }
  
  forest_filename <- paste0(forest_filename, "_", duplication_rate,
                            "_", death_rate, ".sff")
  
  return(forest_filename)
}

palette_six <- function()
{
  p <- magma(6)
  
  return(c(p[1], p[2], p[4], p[6], p[5], p[3]))
  
  my_colors <- c("#00C3C6", "#FF6A68", "#D478FF", "#6BB000")

  palette_gen <- colorRampPalette()

  p <- palette_gen(6)
  
  return(c(p[1], p[2], p[4], p[6], p[5], p[3]))
}

save_image <- function(section, name, extension="png", directory="images",
                       width = 6, height = 4, dpi = 600)
{
  if (!dir.exists(directory)) {
    dir.create(directory)
  }

  filename <- paste0(section, "-", name, ".", extension)
  
  filepath <- file.path(directory, filename)
  ggsave(filepath, width = width, height = height, dpi = dpi)
  
  return(filepath)
}