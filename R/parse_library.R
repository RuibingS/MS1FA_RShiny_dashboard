# Improt CBIO
#' @import MetaboCoreUtils
#' @import doParallel
#' @import enviPat
#'
#' Parse a single compound from a MoNA MSP file
#'
#' This function is to parses a NIST format library.
#' It can filter the compound based on the spectrum type, ion mode, instrument_type, Instrument
#'
#' @param file_path A character vector of lines.
#' @param spectrum_type A string specifying the spectrum type to filter (e.g., "MS1", "MS2"). Set to NULL to skip filtering.
#' @param IonPolarity A string specifying the ion mode to filter (e.g., "pos" for positive, "neg" for negative). Set to NULL to skip filtering.
#' @param InstType A string specifying the instrument_type to filter (e.g., "ESI-TOFF","ESI"). Set to NULL to skip filtering.
#' @param ColEnergy  A string specifying the collision_energy to filter (e.g., "10 V", "20 V"). Set to NULL to skip filtering.
#' @param InstName A string specifying the instrument to filter (e.g., "maXis HD"). Set to NULL to skip filtering.

#' @return A list containing the parsed compound information. If the compound doesn't match the filter, it returns NULL.
#' @export
#'
library(MetaboCoreUtils)
library(doParallel)
library(enviPat)

data(isotopes)



parse_single_compound_NIST <- function(compound_lines,
                                       spectrum_type = NULL,
                                       ionPolarity = NULL,
                                       ioniMethod = NULL,
                                       colEnergy = NULL,
                                       instName = NULL) {
  compound <- list()
  peaks <- NULL

  # skip if "Num Peaks" is missing
  if (!any(grepl("^Num Peaks: ", compound_lines))) {
    return(NULL)
  }

  is_ms2 <- FALSE

  for (line in compound_lines) {


    # parse each line
    if (grepl(": ", line)) {
      infor_value <- unlist(strsplit(line, ": ", fixed = TRUE))
      infor <- trimws(infor_value[1])
      value <- trimws(infor_value[2])
      compound[[infor]] <- value
    }

    # MS level check
    if (grepl("^MSMS:", line)) {
      is_ms2 <- TRUE
    } else {
      is_ms2 <- FALSE
    }

    if (!is.null(spectrum_type)) {
      if ((spectrum_type == "MS1" && is_ms2) || (spectrum_type == "MS2" && !is_ms2)) {
        return(NULL)
      }
    }

    # Apply filter conditions for ionPolarity, ioniMethod, colEnergy, and instName
    if (!is.null(ionPolarity) && grepl("^IonPolarity: ", line)) {
      value <- sub("^IonPolarity: ", "", line)
      if (value != ionPolarity) {
        return(NULL)
      }
    }

    if (!is.null(ioniMethod) && grepl("^IoniMethod: ", line)) {
      value <- sub("^IoniMethod: ", "", line)
      if (value != ioniMethod) {
        return(NULL)
      }
    }

    if (!is.null(colEnergy) && grepl("^ColEnergy: ", line)) {
      value <- sub("^ColEnergy: ", "", line)
      if (value != colEnergy) {
        return(NULL)
      }
    }

    if (!is.null(instName) && grepl("^InstName: ", line)) {
      value <- sub("^InstName: ", "", line)
      if (value != instName) {
        return(NULL)
      }
    }

    # Process peak data
    if (grepl("^Num Peaks: ", line)) {
      peaks_start <- which(compound_lines == line) + 1
      peaks_line <- paste(compound_lines[peaks_start:length(compound_lines)], collapse = " ")
      peaks_values <- as.numeric(unlist(strsplit(peaks_line, "\\s+")))

      peaks_values <- peaks_values[!is.na(peaks_values)]

      if (length(peaks_values) %% 2 != 0) {
        stop("Uneven number of m/z and intensity values in peak data.")
      }

      peaks_df <- data.frame(
        mz = peaks_values[seq(1, length(peaks_values), 2)],
        rel_intensity = peaks_values[seq(2, length(peaks_values), 2)]
      )

      compound[["Peaks"]] <- peaks_df
    }

    # calculate exact mass using enviPat package
    if (grepl("^Formula: ", line)) {
      infor_value <- unlist(strsplit(line, ": ", fixed = TRUE))
      formula <- trimws(infor_value[2])

      if (!exists("isotopes", where = "package:enviPat")) {
        data("isotopes", package = "enviPat", envir = environment())
      }

      exactMass_value <- tryCatch({
        enviPat::isopattern(
          isotopes,
          formula,
          threshold = 0.1,
          plotit = FALSE,
          charge = FALSE,
          emass = 0.00054858,
          algo = 1
        )[[1]][1]
      }, warning = function(w) {
        message("Warning encountered in isopattern: ", conditionMessage(w))
        NA
      }, error = function(e) {
        message("Error in isopattern: ", conditionMessage(e))
        NA
      })

      compound[["ExactMass"]] <- exactMass_value
    }
  }

  return(compound)
}




parse_library_file_parallel_NIST <- function(file_path,
                                             spectrum_type = NULL,
                                             ionPolarity = NULL,
                                             ioniMethod = NULL,
                                             colEnergy = NULL,
                                             instName = NULL) {

  if (!file.exists(file_path)) {
    stop("Error: The file does not exist. Please check the file path.")
  }


  lines <- readLines(file_path, warn = FALSE)


  name_positions <- grep("^Name: ", lines)

  num_compounds <- length(name_positions)

  required_packages <- c("MetaboCoreUtils", "doParallel", "enviPat","stringr")
  n_cores <- detectCores() - 1
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)


  chunk_size <- ceiling(num_compounds / n_cores)
  chunk_indices <- split(seq_along(name_positions), ceiling(seq_along(name_positions) / chunk_size))


  process_chunk <- function(indices) {
    compounds <- list()
    for (i in indices) {
      start <- name_positions[i]
      end <- ifelse(i < length(name_positions), name_positions[i+1]-1 , length(lines))

      compound_lines <- lines[start:end]
      compound <- parse_single_compound_NIST(compound_lines, spectrum_type,
                                             ionPolarity,
                                             ioniMethod,
                                             colEnergy,
                                             instName)


      if (!is.null(compound)) {
        # check x$ExactMass and other numeric fields
        if (!is.null(compound$ExactMass) && is.numeric(compound$ExactMass)) {
          print(compound$ExactMass)
        } else {
          message("ExactMass is NULL or not numeric for this compound.")
        }

        compounds <- append(compounds, list(compound))
      }
    }
    return(compounds)
  }


  clusterExport(cl, varlist = c("parse_single_compound_NIST"))


  compounds <- foreach(i = chunk_indices, .combine = c, .packages = required_packages
  ) %dopar% {
    process_chunk(i)
  }


  stopCluster(cl)

  return(compounds)
}








