#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
rm(list = ls())

# TODO problème pdflatex sur pdf ordi spectro -> à installer!

library(openxlsx)
library(smoothr)
library(shiny)
library(DT)
library(rmarkdown)
library(ggplot2)
# library(pdftools)
# install.packages("pdftools")

# constants
CSV_SKIP_ROWS = 16
SMOOTH_GAUSSIAN_SIGMA = 3
INF_POINT_LEFT_LLOW = 353 # nm
INF_POINT_LEFT_LHIGH = 473 # nm  29.07.2025 => moved from 400 to 473 according to Collet et al. 2012
INF_POINT_LEFT_DEFAULT = round((INF_POINT_LEFT_LLOW + INF_POINT_LEFT_LHIGH) / 2)
INF_POINT_RIGHT_LLOW = 479 # nm  29.07.2025 => moved from 450 to 479 according to Collet et al. 2012
INF_POINT_RIGHT_LHIGH = 550 # nm
INF_POINT_RIGHT_DEFAULT = round((INF_POINT_RIGHT_LLOW + INF_POINT_RIGHT_LHIGH) / 2)
CHECK_LEFT = 350
CHECK_RIGHT = 750
HB_LAMBDA = 415
BL_LAMBDA = 476
WARNING_LEVEL_OK = 0
WARNING_LEVEL_WARN = 1
WARNING_LEVEL_STOP = 2
RESIDUALS_THRESHOLD_WARNING = 1e-4
RESIDUALS_THRESHOLD_ERROR = 1e-3

# functions
load_data_xl <- function(file_path) {
  # dummy path for debugging
  # file_path = "C:/Users/flori/OneDrive - univ-angers.fr/Documents/Home/Hospital/Spectre LCR/R_LCR/data/test 01-01-2026/Spectre LCR + bili HEUVELINE Maxime 24.01.1989.xlsx"
  # load and format accordingly
  absorbance_df <- read.xlsx(file_path, startRow=CSV_SKIP_ROWS)
  # keep only desired columns
  absorbance_df <- absorbance_df[ ,c(1, ncol(absorbance_df))]
  # rename columns
  colnames(absorbance_df) <- c("lambda", "raw_abs")
  # return
  absorbance_df
}

load_data_tst <- function(file_path) {
  # dummy path for debugging
  # file_path = "C:/Users/flori/OneDrive - univ-angers.fr/Documents/Home/Hospital/Spectre LCR/R_LCR/data/test 01-01-2026/Spectre LCR + bili HEUVELINE Maxime 24.01.1989.tst"
  # load raw text data
  rawfiledata <- readLines(file_path, skipNul=T)
  # reencode
  rawfiledata <- iconv(rawfiledata, from = "ISO-8859-1", to = "UTF-8")
  # decode data
  file_values <- strsplit(rawfiledata, " ")[[1]]
  warning("Ignore warning right below: converting raw values read in .tst file to numeric, str parts will produce NAs.")
  file_values <- as.numeric(file_values)
  read_index = 1
  n_traces = 0
  results_df <- data.frame(lambda=seq(350, 750, 2))
  while(TRUE) {
    # check if there are still values to read
    if (read_index + 1 > length(file_values)) {
      break
    }
    # read item at current position
    cur_value <- file_values[read_index]
    next_value <- file_values[read_index + 1]
    # check no NA (numeric values only)
    if ((!is.na(cur_value)) & (!is.na(next_value))) {
      if ((cur_value == "3500") & (next_value != "7500")) {
        trace_l <- numeric(0)
        trace_t <- numeric(0)
        trace_valid = TRUE
        for(trace_index in seq(350, 750, 2)) {
          if (read_index + 1 > length(file_values)) {
            print("File ended in the middle of trace data, stopping trace read")
            trace_valid = F
            break
          }
          cur_value <- file_values[read_index] / 10
          next_value <- file_values[read_index + 1]
          if (is.na(cur_value)) {
            print(paste0("NA lambda found at position: ", read_index, ", expected: ", trace_index, ", stopping trace read"))
            trace_valid = F
            break
          }
          if (is.na(next_value)) {
            print(paste0("NA transmittance found at position: ", read_index, ", expected: ", trace_index, ", stopping trace read"))
            trace_valid = F
            break
          }
          if (cur_value != trace_index) {
            print(paste0("Inconsistent lambda found at position: ", read_index, ", expected: ", trace_index, ", got: ", cur_value, ", stopping trace read"))
            trace_valid = F
            break
          }
          trace_l <- c(trace_l, cur_value)
          trace_t <- c(trace_t, next_value)
          read_index = read_index + 2
        }
        # QC trace
        if (trace_valid) {
          if (nrow(results_df) != length(trace_l)) {
            print("Trace data length is inconsistent, discarding trace")
            trace_valid = F
          }
          if (any(results_df$lambda != trace_l)) {
            print("Trace data lambdas list is inconsistent, discarding trace")
            trace_valid = F
          }
          if (any(is.na(trace_l))) {
            print("Trace data lambdas list contains missing values, discarding trace")
            trace_valid = F
          }
          if (any(is.na(trace_t))) {
            print("Trace data transmisttances list contains missing values, discarding trace")
            trace_valid = F
          }
        }
        # store trace if valid
        if (trace_valid) {
          n_traces = n_traces + 1
          results_df[, paste0("t_", n_traces)] <- trace_t
        }
        # next steps
      }
    }
    read_index = read_index + 1
  }
  # check that everything worked well
  if (nrow(results_df) != 201) {
    stop(paste0("Inconsistent number of rows after reading tst file data [found nrow=", nrow(results_df), "], discarding loaded data"))
  }
  # check that everything worked well
  if (ncol(results_df) < 2) {
    stop(paste0("No valid data found [found ncol=", ncol(results_df), " < 2], discarding loaded data"))
  }
  if (ncol(results_df) != 3) {
    warning(paste0("Inconsistent number of columns [found ncol=", ncol(results_df), "] after reading tst file data!"))
  }
  # post-process
  for(trace_index in 1:n_traces) {
    results_df[, paste0("abs_", trace_index)] <- log10(1 / results_df[, paste0("t_", n_traces)])
  }
  # encapsulate into a final df with correct column names
  absorbance_df <- data.frame(lambda=results_df$lambda, raw_abs=results_df[, paste0("abs_", n_traces)])
  # return
  absorbance_df
}

load_data <- function(file_path) {
  ext <- tools::file_ext(file_path)
  if (ext == "xlsx") {
    return(load_data_xl(file_path=file_path))
  }
  if (ext == "tst") {
    return(load_data_tst(file_path=file_path))
  }
  stop(paste0("Unhandled extension: ", ext, " for file ", file_path))
}

# called when file is loaded, to make sure the absorbance_df contains all required info
intialize_spectrum_data <- function(absorbance_df) {
  # interpolate so we make sure we have data between 350 and 750
  interpolated_data <- approx(x=absorbance_df$lambda, y=absorbance_df$raw_abs, xout=c(min(absorbance_df$lambda):max(absorbance_df$lambda)))
  absorbance_df <- data.frame(lambda=interpolated_data$x,
                              raw_abs=interpolated_data$y)

  # add smoothed version of the trace
  # new version -> packages "smoothr"
  absorbance_df$smooth_abs = ksmooth(x=absorbance_df$lambda, absorbance_df$raw_abs, bandwidth=SMOOTH_GAUSSIAN_SIGMA)$y
  # old version -> package "smoother"
  # absorbance_df$smooth_abs = smth.gaussian(absorbance_df$raw_abs, window=SMOOTH_GAUSSIAN_SIGMA, tails=T)
  
  # by default tangent is missing (no value), but still initialized
  absorbance_df$tangent_auto = NA
  absorbance_df$residuals_auto = NA
  absorbance_df$tangent_manual = NA
  absorbance_df$residuals_manual = NA
  
  # return
  absorbance_df
}

# also do a manual tangent method
create_tangent_data <- function(absorbance_df, left_lambda, right_lambda) {
  lloc = which(absorbance_df$lambda == left_lambda)
  rloc = which(absorbance_df$lambda == right_lambda)
  
  # define candidate tangent function
  x2 <- right_lambda
  x1 <- left_lambda
  y2 <- absorbance_df$smooth_abs[rloc]
  y1 <- absorbance_df$smooth_abs[lloc]
  a=(y2-y1)/(x2-x1)
  b=y2-a*x2
  # compute tangent values
  tangent_values <- absorbance_df$lambda * a + b
  # subtract tangent from abs curve
  spectrum_residuals <- absorbance_df$smooth_abs - tangent_values

  # return
  list(tangent_values=tangent_values, residuals=spectrum_residuals)
}

# called to automatically create a tangent and return its parameters
get_auto_best_tangent_parameters <- function(absorbance_df) {
  # create a copy of the spectrum_data (tangent data will be stored in the same)
  absorbance_df <- cbind(absorbance_df)

  # define inflexion points candidates windows
  check_filter <- (absorbance_df$lambda >= CHECK_LEFT) & (absorbance_df$lambda <= CHECK_RIGHT)
  left_infpoint_filter <- (absorbance_df$lambda >= INF_POINT_LEFT_LLOW) & (absorbance_df$lambda <= INF_POINT_LEFT_LHIGH)
  right_infpoint_filter <- (absorbance_df$lambda >= INF_POINT_RIGHT_LLOW) & (absorbance_df$lambda <= INF_POINT_RIGHT_LHIGH)
  
  # we will also have to remember the locations (two points) from which the tangent is created
  best_left_lambda = NA
  best_right_lambda = NA

  # test each possible combination
  min_residuals <- 999
  for (lloc in which(left_infpoint_filter)) {
    for (rloc in which(right_infpoint_filter)) {
      # define candidate tangent function
      x2 <- absorbance_df$lambda[rloc]
      x1 <- absorbance_df$lambda[lloc]
      y2 <- absorbance_df$smooth_abs[rloc]
      y1 <- absorbance_df$smooth_abs[lloc]
      a=(y2-y1)/(x2-x1)
      b=y2-a*x2
      # compute tangent values
      tangent_values <- absorbance_df$lambda * a + b
      # check if tangent crosses the line
      # subtract tangent from abs curve
      spectrum_residuals <- absorbance_df$smooth_abs - tangent_values
      # n_crossing_points <- sum((spectrum_residuals[1:(length(spectrum_residuals) - 1)] < 0) == (spectrum_residuals[2:length(spectrum_residuals)] >= 0))
      sum_residuals <- sum(-spectrum_residuals[check_filter][spectrum_residuals[check_filter] < 0])
      if (sum_residuals < min_residuals) {
        # keep this tangent as the best
        min_residuals <- sum_residuals
        best_left_lambda = x1
        best_right_lambda = x2
      }
    }
  }
  
  # return
  list(left_lambda=best_left_lambda, right_lambda=best_right_lambda)
}

compute_quantitative_results <- function(absorbance_df, method) {
  computation_results <- list()
  
  warning_level = WARNING_LEVEL_OK
  warning_items = list()
  if (all(absorbance_df$raw_abs == -1)) {
    warning_items[[length(warning_items) + 1]] <- "Données invalides: impossible de charger le fichier spectre."
    warning_level <- max(warning_level, WARNING_LEVEL_STOP)
  }
  if (any(absorbance_df$raw_abs < 0)) {
    warning_items[[length(warning_items) + 1]] <- "Valeurs d'absorbance < 0!"
    warning_level <- max(warning_level, WARNING_LEVEL_WARN)
  }
  if (any(absorbance_df$raw_abs[(absorbance_df$lambda < 550)] < 0)) {
    warning_items[[length(warning_items) + 1]] <- "Valeurs d'absorbance < 0 entre 350 et 550nm!"
    warning_level <- max(warning_level, WARNING_LEVEL_STOP)
  }
  if (method!="auto") {
    warning_items[[length(warning_items) + 1]] <- "Méthode manuelle choisie: vérifier la cohérence des résultats!"
    warning_level <- max(warning_level, WARNING_LEVEL_WARN)
  }
  if (method=="auto") {
    tg_values <- absorbance_df$tangent_auto
  } else {
    tg_values <- absorbance_df$tangent_manual
  }
  tangent_crosses_0_at_lambda <- max(absorbance_df$lambda[tg_values > 0])
  if (tangent_crosses_0_at_lambda < BL_LAMBDA) {
    warning_items[[length(warning_items) + 1]] <- paste0("La tangeante passe par zéro avant ", BL_LAMBDA, "nm!")
    warning_level <- max(warning_level, WARNING_LEVEL_WARN)
  }
  spectrum_residuals <- (absorbance_df$smooth_abs - tg_values)[(absorbance_df$lambda >= INF_POINT_LEFT_LLOW) & (absorbance_df$lambda <= INF_POINT_RIGHT_LHIGH)]
  # n_crossing_points <- sum((spectrum_residuals[1:(length(spectrum_residuals) - 1)] < 0) == (spectrum_residuals[2:length(spectrum_residuals)] >= 0))
  sum_total_residuals <- sum(-spectrum_residuals[spectrum_residuals < 0])
  if (sum_total_residuals > RESIDUALS_THRESHOLD_ERROR) {
    if (method=="auto") {
      warning_items[[length(warning_items) + 1]] <- paste0("Impossible de trouver une tangente sans couper la courbe entre ", INF_POINT_LEFT_LLOW, " et ", INF_POINT_RIGHT_LHIGH, "nm!")
      warning_level <- max(warning_level, WARNING_LEVEL_STOP)
    } else {
      warning_items[[length(warning_items) + 1]] <- paste0("Les positions choisies sont invalides: la tangente passe au-dessus de la courbe entre ", INF_POINT_LEFT_LLOW, " et ", INF_POINT_RIGHT_LHIGH, "nm!")
      warning_level <- max(warning_level, WARNING_LEVEL_STOP)
    }
  } else if (sum_total_residuals > RESIDUALS_THRESHOLD_WARNING) {
    if (method=="auto") {
      warning_items[[length(warning_items) + 1]] <- paste0("Impossible de trouver une tangente sans couper la courbe entre ", INF_POINT_LEFT_LLOW, " et ", INF_POINT_RIGHT_LHIGH, "nm!")
      warning_level <- max(warning_level, WARNING_LEVEL_WARN)
    } else {
      warning_items[[length(warning_items) + 1]] <- paste0("Les positions choisies sont douteuses: la tangente passe au-dessus de la courbe entre ", INF_POINT_LEFT_LLOW, " et ", INF_POINT_RIGHT_LHIGH, "nm!")
      warning_level <- max(warning_level, WARNING_LEVEL_WARN)
    }
  }
  
  hb_trace_do <- absorbance_df$smooth_abs[absorbance_df$lambda == HB_LAMBDA]
  bl_trace_do <- absorbance_df$smooth_abs[absorbance_df$lambda == BL_LAMBDA]
  
  if (method=="auto") {
    hb_tangent_do <- absorbance_df$tangent_auto[absorbance_df$lambda == HB_LAMBDA]
    bl_tangent_do <- absorbance_df$tangent_auto[absorbance_df$lambda == BL_LAMBDA]
  } else {
    hb_tangent_do <- absorbance_df$tangent_manual[absorbance_df$lambda == HB_LAMBDA]
    bl_tangent_do <- absorbance_df$tangent_manual[absorbance_df$lambda == BL_LAMBDA]
  }
  
  hb_delta_do <- hb_trace_do - hb_tangent_do
  bl_delta_do <- bl_trace_do - bl_tangent_do
  
  hb_concl <- ifelse(hb_delta_do >= 0.1, ">0.1", ifelse(hb_delta_do >= 0.02, ">0.02", "<0.02"))
  bl_concl <- ifelse(bl_delta_do >= 0.007, ">0.07", "<0.007")

  com_concl <- "?"
  if (bl_concl == "<0.007") {
    if (hb_concl == ">0.1") {
      com_concl = "C1 : Présence d’hémoglobine en quantité trop importante pouvant masquer la présence de bilirubine. Hémorragie sous-arachnoïdienne non exclue."
    } else if (hb_concl == ">0.02") {
      com_concl = "C3 : Présence d’hémoglobine en faible quantité mais absence de bilirubine. Pas d’argument en faveur d’une hémorragie sous-arachnoïdienne ou compatible avec une hémorragie sous-arachnoïdienne débutante (LCR prélevé moins de 12 h après les premiers signes) ou une ponction lombaire traumatique."
    } else if (hb_concl == "<0.02") {
      com_concl = "C2 : Examen spectrophotométrique négatif. Absence d’hémoglobine et de bilirubine. Pas d’argument en faveur d’une hémorragie sous-arachnoïdienne."
    }
  } else {
    if ((hb_concl == ">0.1") || (hb_concl == ">0.02")) {
      com_concl = "C4 : Présence de bilirubine et d’hémoglobine. Résultat compatible avec une hémorragie sous-arachnoïdienne. Attention, la présence de bilirubine dans le LCR n’est pas spécifique d’une hémorragie sous-arachnoïdienne. D’autres étiologies peuvent être responsables : méningites, hémorragie sous-durale, malignité..."
    } else {
      com_concl = "Delta DO bilirubine positif! Regarder la protéinorachie et la bilirubinémie plasmatique pour conclure."
    }
  }
  
  computation_results <- list(warning_level=warning_level,
                              warning_items=warning_items,
                              hb_trace_do=hb_trace_do, bl_trace_do=bl_trace_do,
                              hb_tangent_do=hb_tangent_do, bl_tangent_do=bl_tangent_do,
                              hb_delta_do=hb_delta_do, bl_delta_do=bl_delta_do,
                              hb_concl=hb_concl, bl_concl=bl_concl,
                              com_concl=com_concl
                              )
  
  computation_results
}

plot_spectrum <- function(lambda, raw_abs, smooth_abs, tangent, tangent_parameters, adjust_scale) {
  absorbance_df <- data.frame(lambda=lambda,
                              raw_abs=raw_abs,
                              smooth_abs=smooth_abs,
                              tangent=tangent)
  
  ggp <- ggplot(data=absorbance_df) +
    geom_line(mapping=aes(x=lambda, y=raw_abs), color="gray", linewidth=.75) +
    geom_line(mapping=aes(x=lambda, y=smooth_abs), color="black", linewidth=.75)
  
  if (adjust_scale) {
    ylim_min <- min(0, absorbance_df$raw_abs)
    ylim_max <- max(0.1, max(absorbance_df$raw_abs))
  } else {
    ylim_min <- min(0, absorbance_df$raw_abs)
    ylim_max <- max(absorbance_df$raw_abs)
  }
  
  if (all(!is.na(absorbance_df$tangent))) {
    ggp <- ggp +
      geom_line(mapping=aes(x=lambda, y=tangent), color="red", linewidth=.75)
    
    ggp <- ggp +
      geom_segment(data=data.frame(x=HB_LAMBDA,
                                   xend=HB_LAMBDA,
                                   y=0,
                                   yend=absorbance_df$raw_abs[absorbance_df$lambda == HB_LAMBDA]),
                   mapping=aes(x=x, y=y, xend=xend, yend=yend), color="black", linetype='dotted', linewidth=1) +
      geom_segment(data=data.frame(x=BL_LAMBDA,
                                   xend=BL_LAMBDA,
                                   y=0,
                                   yend=absorbance_df$raw_abs[absorbance_df$lambda == BL_LAMBDA]),
                   mapping=aes(x=x, y=y, xend=xend, yend=yend), color="black", linetype='dotted', linewidth=1)
    
    ggp <- ggp +
      geom_point(data=data.frame(x=c(tangent_parameters$left_lambda,
                                     tangent_parameters$right_lambda),
                                 y=c(absorbance_df$raw_abs[absorbance_df$lambda == tangent_parameters$left_lambda],
                                     absorbance_df$raw_abs[absorbance_df$lambda == tangent_parameters$right_lambda])),
                 mapping=aes(x=x, y=y), color="red", size=2)
  }

  ggp <- ggp +
    ylim(ylim_min, ylim_max) +
    xlab("Longueur d'onde (nm)") +
    ylab("Absorbance") +
    guides() +
    theme_minimal()
  
  ggp
}


################################################################################
################################################################################


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Recherche de pigments biliaires dans le LCR"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          fileInput("file_input_tst", "Choisir un fichier spectre", accept = ".tst"),
          fileInput("file_input_xl", "Choisir un fichier .xlsx", accept = ".xlsx"),
          uiOutput("spectrum_parameter_inputs"),
          uiOutput("spectrum_parameter_manual_inputs"),
          uiOutput("spectrum_inputs_outputs"),
        ),

        # Show a plot of the generated distribution
        mainPanel(
          fluidPage(
            plotOutput("main_plot"),
            # plotlyOutput("main_plot_plotly"),
          )
        )
    )
)


################################################################################
################################################################################


# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
    reactive_sample_data <- reactiveValues(data_loaded = NULL,
                                           original_fn = NULL,
                                           absorbance_df = NULL,
                                           auto_tangent_parameters = NULL,
                                           auto_computation_results = NULL,
                                           manual_computation_results = NULL,
                                           temp_ggp = NULL,
                                           temp_pltly = NULL,
                                           )
    
    # file input (XL version)
    observeEvent(input$file_input_xl, {
      file <- input$file_input_xl
      ext <- tools::file_ext(file$datapath)
      
      req(file)
      validate(need(ext == "xlsx", "Choisir un fichier .xlsx!"))
      
      # load and format accordingly
      file_path <- file$datapath
      # file_path = "C:/Users/flori/OneDrive - univ-angers.fr/Documents/Home/Hospital/Spectre LCR/R_LCR/data/test 01-01-2026/Spectre LCR + bili BOURRIGAULT Clément 18-12-2003.tst"
      absorbance_df <- NULL
      tryCatch({
        absorbance_df <- load_data(file_path)
      }, error = function(e) {
        showModal(modalDialog(
          title = "Impossible de charger le fichier spectre",
          paste0("Erreur: ", e$message),
          easyClose = TRUE,
          footer = NULL
        ))
      })
      
      req(absorbance_df)
      
      # post process data
      absorbance_df <- intialize_spectrum_data(absorbance_df=absorbance_df)

      # auto tangent
      auto_tangent_parameters <- get_auto_best_tangent_parameters(absorbance_df)
      
      # create tangent from auto input
      auto_tangent_output <- create_tangent_data(absorbance_df=absorbance_df,
                                                 left_lambda=auto_tangent_parameters$left_lambda,
                                                 right_lambda=auto_tangent_parameters$right_lambda)

      # store
      absorbance_df$tangent_auto <- auto_tangent_output$tangent_values
      absorbance_df$residuals_auto <- auto_tangent_output$residuals
      
      # compute quantification results
      computation_results <- compute_quantitative_results(absorbance_df, method="auto")
      
      # store
      reactive_sample_data$absorbance_df = absorbance_df
      reactive_sample_data$auto_tangent_parameters = auto_tangent_parameters
      reactive_sample_data$auto_computation_results = computation_results
      reactive_sample_data$original_fn <- tools::file_path_sans_ext(basename(file$name))
      reactive_sample_data$data_loaded <- TRUE
    })
    
    # file input (.tst version)
    observeEvent(input$file_input_tst, {
      file <- input$file_input_tst
      ext <- tools::file_ext(file$datapath)
      
      req(file)
      validate(need(ext == "tst", "Choisir un fichier .tst!"))
      
      # load and format accordingly
      file_path <- file$datapath
      # load data
      absorbance_df <- NULL
      tryCatch({
        absorbance_df <- load_data(file_path)
      }, error = function(e) {
        showModal(modalDialog(
          title = "Impossible de charger le fichier spectre",
          paste0("Erreur: ", e$message),
          easyClose = TRUE,
          footer = NULL
        ))
      })
      
      req(absorbance_df)
      
      # post process data
      absorbance_df <- intialize_spectrum_data(absorbance_df=absorbance_df)
      
      # auto tangent
      auto_tangent_parameters <- get_auto_best_tangent_parameters(absorbance_df)

      # create tangent from auto input
      auto_tangent_output <- create_tangent_data(absorbance_df=absorbance_df,
                                                 left_lambda=auto_tangent_parameters$left_lambda,
                                                 right_lambda=auto_tangent_parameters$right_lambda)

      # store
      absorbance_df$tangent_auto <- auto_tangent_output$tangent_values
      absorbance_df$residuals_auto <- auto_tangent_output$residuals
      
      # compute quantification results
      computation_results <- compute_quantitative_results(absorbance_df, method="auto")

      # store
      reactive_sample_data$absorbance_df = absorbance_df
      reactive_sample_data$auto_tangent_parameters = auto_tangent_parameters
      reactive_sample_data$auto_computation_results = computation_results
      reactive_sample_data$original_fn <- tools::file_path_sans_ext(basename(file$name))
      reactive_sample_data$data_loaded <- TRUE
    })
    
    # ggplot2 output
    output$main_plot <- renderPlot({
      req(reactive_sample_data$data_loaded)
      
      if (is.null(input$auto_tangent_checkbox)) {
        return(NULL)
      }
      
      graphics_tangent_parameters = NULL
      if (input$auto_tangent_checkbox) {
        graphics_tangent_parameters <- reactive_sample_data$auto_tangent_parameters
        tangent_values <- reactive_sample_data$absorbance_df$tangent_auto
      } else {
        graphics_tangent_parameters <- list(left_lambda=input$left_tangent_lambda,
                                            right_lambda=input$right_tangent_lambda)
        tangent_values <- reactive_sample_data$absorbance_df$tangent_manual
      }
      
      # make figure
      ggp <- plot_spectrum(lambda=reactive_sample_data$absorbance_df$lambda,
                           raw_abs=reactive_sample_data$absorbance_df$raw_abs,
                           smooth_abs=reactive_sample_data$absorbance_df$smooth_abs,
                           tangent=tangent_values,
                           tangent_parameters=graphics_tangent_parameters,
                           adjust_scale=ifelse(is.null(input$adjust_scale), F, input$adjust_scale))
      
      reactive_sample_data$temp_ggp <- ggp
      
      # return
      reactive_sample_data$temp_ggp
    },
    width = "auto",
    height = 800,
    res=100)
    
    output$spectrum_parameter_inputs <- renderUI({
      req(reactive_sample_data$data_loaded)
      
      tagList(
        htmlOutput("warning_output"),
        HTML("<strong>Paramètres</strong>"),
        checkboxInput("adjust_scale", "Ajuster l'échelle", TRUE),
        checkboxInput("auto_tangent_checkbox", "Déterminer automatiquement la tangente", TRUE),
      )
    })
    
    output$spectrum_parameter_manual_inputs <- renderUI({
      req(reactive_sample_data$data_loaded)
      
      if (is.null(input$auto_tangent_checkbox)) {
        return(NULL)
      }
      
      if (input$auto_tangent_checkbox) {
        return (NULL)
      }
      
      tagList(
        sliderInput("left_tangent_lambda", "Point gauche", INF_POINT_LEFT_LLOW, INF_POINT_LEFT_LHIGH, INF_POINT_LEFT_DEFAULT),
        sliderInput("right_tangent_lambda", "Point droit", INF_POINT_RIGHT_LLOW, INF_POINT_RIGHT_LHIGH, INF_POINT_RIGHT_DEFAULT),
      )
    })
    
    observeEvent(input$left_tangent_lambda, {
      # create tangent from auto input
      manual_tangent_output <- create_tangent_data(absorbance_df=reactive_sample_data$absorbance_df,
                                                   left_lambda=input$left_tangent_lambda,
                                                   right_lambda=input$right_tangent_lambda)
      
      # store
      reactive_sample_data$absorbance_df$tangent_manual <- manual_tangent_output$tangent_values
      reactive_sample_data$absorbance_df$residuals_manual <- manual_tangent_output$residuals
      
      # compute quantification results
      computation_results <- compute_quantitative_results(reactive_sample_data$absorbance_df, method="manual")
      
      # store
      reactive_sample_data$manual_computation_results = computation_results
    })
    
    observeEvent(input$right_tangent_lambda, {
      # create tangent from auto input
      manual_tangent_output <- create_tangent_data(absorbance_df=reactive_sample_data$absorbance_df,
                                                   left_lambda=input$left_tangent_lambda,
                                                   right_lambda=input$right_tangent_lambda)
      
      # store
      reactive_sample_data$absorbance_df$tangent_manual <- manual_tangent_output$tangent_values
      reactive_sample_data$absorbance_df$residuals_manual <- manual_tangent_output$residuals
      
      # compute quantification results
      computation_results <- compute_quantitative_results(reactive_sample_data$absorbance_df, method="manual")
      
      # store
      reactive_sample_data$manual_computation_results = computation_results
    })

    output$spectrum_inputs_outputs <- renderUI({
      req(reactive_sample_data$data_loaded)
      # if (is.null(input$auto_tangent_checkbox)) {
      #   req(FALSE)
      # }
      if (is.null(input$auto_tangent_checkbox)) {
        return(NULL)
      }
      if (input$auto_tangent_checkbox) {
        computation_results <- reactive_sample_data$auto_computation_results
      } else {
        computation_results <- reactive_sample_data$manual_computation_results
      }
      req(computation_results$warning_level < WARNING_LEVEL_STOP)
      
      tagList(
        HTML("<strong>Positions de la tangente</strong>"),
        DTOutput("tangent_dt"),
        HTML("<br>"),
        HTML("<strong>Calculs</strong>"),
        DTOutput("computs_dt"),
        HTML("<strong>Résultats</strong>"),
        DTOutput("results_dt"),
        HTML("<br>"),
        HTML("<strong>Conclusion</strong>"),
        HTML("<br>"),
        htmlOutput("comm_output"),
        HTML("<br>"),
        fluidRow(
          div(
            style = "display: flex; justify-content: center; align-items: center;",
            downloadButton("report", "Générer PDF")
          )
        )
      )
    })
    
    # dt output
    output$tangent_dt <- renderDT({
      req(reactive_sample_data$data_loaded)
      
      graphics_tangent_parameters = NULL
      tangent_method = "?"
      if (input$auto_tangent_checkbox) {
        tangent_method = "automatique"
        graphics_tangent_parameters <- reactive_sample_data$auto_tangent_parameters
        computation_results <- reactive_sample_data$auto_computation_results
      } else {
        tangent_method = "manuelle"
        computation_results <- reactive_sample_data$manual_computation_results
        graphics_tangent_parameters <- list(left_lambda=input$left_tangent_lambda,
                                            right_lambda=input$right_tangent_lambda)
      }
      
      req(graphics_tangent_parameters)
      req(computation_results)
      req(computation_results$warning_level < WARNING_LEVEL_STOP)
      
      tangent_parameters <- reactive_sample_data$tangent_parameters
      output_dt <- data.frame(`Tangente`=c("Méthode",
                                           "λ point gauche (nm)",
                                           "λ point droit (nm)"),
                              `Valeur`=c(tangent_method,
                                         paste0(graphics_tangent_parameters$left_lambda, " nm"),
                                         paste0(graphics_tangent_parameters$right_lambda, " nm"))
      )
      
      # return
      datatable(output_dt, options = list(dom = 't'))
    })
    
    # dt output
    output$computs_dt <- renderDT({
      req(reactive_sample_data$data_loaded)
      
      if (input$auto_tangent_checkbox) {
        computation_results <- reactive_sample_data$auto_computation_results
      } else {
        computation_results <- reactive_sample_data$manual_computation_results
      }
      
      req(computation_results)
      req(computation_results$warning_level < WARNING_LEVEL_STOP)
      
      output_dt <- data.frame(`Paramètre`=c("DO hémoglobine (courbe)",
                                            "DO hémoglobine (tangente)",
                                            "DO bilirubine (courbe)",
                                            "DO bilirubine (tangente)"),
                              `Mesure`=c(round(computation_results$hb_trace_do, 4),
                                         round(computation_results$hb_tangent_do, 4),
                                         round(computation_results$bl_trace_do, 4),
                                         round(computation_results$bl_tangent_do, 4)))
      
      # return
      datatable(output_dt, options = list(dom = 't'))
    })
    
    # dt output
    output$results_dt <- renderDT({
      req(reactive_sample_data$data_loaded)
      
      if (input$auto_tangent_checkbox) {
        computation_results <- reactive_sample_data$auto_computation_results
      } else {
        computation_results <- reactive_sample_data$manual_computation_results
      }
      
      req(computation_results)
      req(computation_results$warning_level < WARNING_LEVEL_STOP)
      
      output_dt <- data.frame(`Paramètre`=c("Delta DO hémoglobine",
                                            "Delta DO bilirubine"),
                              `Résultat`=c(round(computation_results$hb_delta_do, 4),
                                           round(computation_results$bl_delta_do, 4)),
                              `Interprétation`=c(computation_results$hb_concl,
                                                 computation_results$bl_concl))

      # return
      datatable(output_dt, options = list(dom = 't'))
    })
    
    # comm output
    output$comm_output <- renderUI({
      req(reactive_sample_data$data_loaded)
      
      if (input$auto_tangent_checkbox) {
        computation_results <- reactive_sample_data$auto_computation_results
      } else {
        computation_results <- reactive_sample_data$manual_computation_results
      }
      
      req(computation_results)
      req(computation_results$warning_level < WARNING_LEVEL_STOP)
      
      HTML(computation_results$com_concl)
    })
    
    # warning output
    output$warning_output <- renderUI({
      req(reactive_sample_data$data_loaded)
      
      if (input$auto_tangent_checkbox) {
        computation_results <- reactive_sample_data$auto_computation_results
      } else {
        computation_results <- reactive_sample_data$manual_computation_results
      }
      
      req(computation_results)
      
      if (computation_results$warning_level == WARNING_LEVEL_OK) {
        return(HTML(""))
      } else if (computation_results$warning_level == WARNING_LEVEL_WARN) {
        warn_color = "orange"
      } else {
        warn_color = "red"
      }
      
      warning_items <- computation_results$warning_items
      warning_text <- paste0("Avertissements:<ul>", paste(paste0("<li>", warning_items, "</li>"), collapse=""), "</ul>")
      warning_text <- paste0("<b style='color:", warn_color, ";'>", warning_text, "</b>")
      HTML(warning_text)
    })
    
    # generate report
    output$report <- downloadHandler(
      filename = function() {
        paste(reactive_sample_data$original_fn, sep = '.', "pdf")
      },
      
      content = function(file) {
        src <- normalizePath('report.Rmd')
        
        # temporarily switch to the temp dir, in case you do not have write
        # permission to the current working directory
        owd <- setwd(tempdir())
        on.exit(setwd(owd))
        file.copy(src, 'report.Rmd', overwrite = TRUE)
        
        out <- render('report.Rmd', pdf_document())
        file.rename(out, file)
      }
    )
    
    session$onSessionEnded(function() {
      stopApp()
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
