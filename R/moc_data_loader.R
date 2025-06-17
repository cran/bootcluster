#' Load Multiple-Objective Clustering (MOC) Datasets
#'
#' @description Loads and processes datasets for multiple-objective clustering analysis.
#' The function loads CSV files from a specified directory and processes them by removing
#' NA columns.
#'
#' @param data_dir Directory containing the CSV datasets (default: current working directory)
#'
#' @return A list containing:
#' \describe{
#'   \item{datasets}{Named list of processed datasets}
#' }
#'
#' @importFrom stats prcomp
#' @importFrom utils read.csv
#'
#' @examples
#' \dontrun{
#' # Load datasets
#' result <- load_moc_datasets("path/to/MOC_Data")
#' 
#' # Access a specific dataset
#' spiral <- result$datasets$Spiral
#' }
#'
#' @export
load_moc_datasets <- function(data_dir = getwd()) {
    
    # Save current working directory
    original_dir <- getwd()
    on.exit(setwd(original_dir), add = TRUE)
    
    # Change to data directory
    setwd(data_dir)
    
    # Get the dataset files
    csv_files <- list.files(pattern = "[.]csv$")
    if (length(csv_files) == 0) {
        stop("No CSV files found in directory: ", data_dir)
    }
    
    datasets <- lapply(csv_files, read.csv)
    
    # Get the dataset files names for renaming
    filenames_extensions <- csv_files
    
    # Remove the extension
    filenames <- sub(pattern = "(.*)\\..*$", replacement = "\\1", filenames_extensions)
    
    # Set names
    names(datasets) <- filenames
    
    # Remove all NA's columns
    for (i in 1:length(datasets)) {
        temp <- datasets[[i]]
        # remove all NA columns
        temp_new <- temp[, colSums(is.na(temp)) < nrow(temp)]
        datasets[[i]] <- temp_new
    }
    
    
    # Return results
    return(list(
        datasets = datasets
    ))
}
