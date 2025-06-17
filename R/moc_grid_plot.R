#' Create a Grid Plot of MOC Results
#'
#' @description Creates a grid plot with datasets as rows and clustering methods as columns.
#' This function is designed to visualize multiple datasets and methods in a single plot.
#'
#' @param results Results from analyze_moc_datasets
#' @param dataset_names Names of datasets to plot (default: all datasets in results)
#' @param methods Methods to plot (default: all available methods)
#' @param plot_file Output file for plots (default: NULL)
#' @param format Output format, either "pdf" or "eps" (default: "pdf")
#' @param mar Margins for plots (default: c(2, 2, 2, 1))
#' @param cex Text size multiplier (default: 0.7)
#' @param point_size Point size for scatter plots (default: 0.8)
#' @param family Font family (default: "serif" for Times New Roman)
#' @param label_style Whether to add row/column labels (default: TRUE)
#' @param maintain_aspect_ratio Whether to maintain aspect ratio in PDF (default: TRUE)
#'
#' @return Invisibly returns the layout information
#'
#' @importFrom grDevices rainbow dev.off pdf postscript
#' @importFrom graphics par plot title
#'
#' @export
plot_moc_grid <- function(results, 
                        dataset_names = NULL, 
                        methods = c("kmeans", "hierarchical", "spectral", 
                                  "fastgreedy", "metis", "hmetis"),
                        plot_file = NULL,
                        format = "pdf",
                        mar = c(2, 2, 2, 1),
                        cex = 0.7,
                        point_size = 0.8,
                        family = "serif",
                        label_style = TRUE,
                        maintain_aspect_ratio = TRUE) {
    
    # Get all dataset names if not specified
    if (is.null(dataset_names)) {
        dataset_names <- names(results$results)
    }
    
    # Filter out any dataset names that don't exist in results
    valid_datasets <- names(results$results)
    dataset_names <- dataset_names[dataset_names %in% valid_datasets]
    
    if (length(dataset_names) == 0) {
        stop("No valid dataset names provided")
    }
    
    # Filter out datasets with errors
    valid_datasets <- sapply(dataset_names, function(name) {
        result <- results$results[[name]]
        !is.null(result) && is.null(result$error)
    })
    dataset_names <- dataset_names[valid_datasets]
    
    if (length(dataset_names) == 0) {
        stop("No valid datasets without errors")
    }
    
    # Calculate dimensions to fit on one page while maintaining aspect ratio
    n_datasets <- length(dataset_names)
    n_methods <- length(methods)
    
    if (maintain_aspect_ratio) {
        # Standard paper size (letter) is 8.5 x 11 inches
        # Calculate width and height to fit on one page
        max_width <- 8
        max_height <- 10.5
        
        # Calculate aspect ratio based on number of rows and columns
        aspect_ratio <- n_datasets / n_methods
        
        # Calculate dimensions to fit on one page
        if (aspect_ratio > max_height / max_width) {
            # Height limited
            height <- max_height
            width <- height / aspect_ratio
        } else {
            # Width limited
            width <- max_width
            height <- width * aspect_ratio
        }
    } else {
        # Default dimensions
        width <- 15
        height <- 3 * n_datasets
    }
    
    # Open graphics device if plot_file is provided
    if (!is.null(plot_file)) {
        # Check if file extension is provided
        if (!grepl("\\.[^\\.]+$", plot_file)) {
            # Add extension based on format
            plot_file <- paste0(plot_file, ".", format)
        }
        
        if (tolower(format) == "eps") {
            postscript(file = plot_file, width = width, height = height, 
                      horizontal = FALSE, onefile = FALSE, paper = "special",
                      family = family)
        } else {
            # Default to PDF
            pdf(file = plot_file, width = width, height = height, family = family)
        }
    }
    
    # Save old par settings and restore on exit
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)
    
    # Set up the layout
    par(mfrow = c(n_datasets, n_methods), 
        mar = mar, 
        oma = c(2, 2, 1, 1),  # Reduced top outer margin since no title
        cex = cex,
        cex.main = 0.9,
        cex.axis = 0.7,
        cex.lab = 0.7,
        family = family)  # Set font family to Times New Roman
    
    # Plot each dataset as a row
    for (i in 1:n_datasets) {
        dataset_name <- dataset_names[i]
        result <- results$results[[dataset_name]]
        
        # Skip if result has error
        if (is.null(result) || !is.null(result$error)) {
            next
        }
        
        pca_plot <- result$pca_plot
        k <- length(unique(result$true_clusters))
        
        # Plot each method as a column
        for (j in 1:n_methods) {
            method <- methods[j]
            
            # Get membership data
            if (method %in% c("kmeans", "hierarchical", "spectral")) {
                membership <- result$ensemble_results$product$individual_results[[method]]$membership
            } else {
                membership <- result$ensemble_results$product[[method]]$membership
            }
            
            # Format method name properly (k-means, METIS and hMETIS)
            display_method <- method
            if (method == "kmeans") {
                display_method <- "k-means"
            } else if (method == "metis") {
                display_method <- "METIS"
            } else if (method == "hmetis") {
                display_method <- "hMETIS"
            }
            
            # Create plot title with row/column labels if requested
            if (label_style) {
                # Create labels like A1, A2, B1, B2, etc.
                row_label <- LETTERS[i]
                col_label <- j
                # Simplified title with just the label and method name
                plot_title <- paste0(row_label, col_label, ") ", display_method)
            } else {
                plot_title <- paste(dataset_name, "k =", k, display_method)
            }
            
            # Create the plot
            plot(pca_plot,
                 col = rainbow(k)[membership],
                 pch = 20,  # Use small filled circles
                 cex = point_size,
                 main = plot_title,
                 xlab = "",
                 ylab = "")
            
            # Add axis labels only for the leftmost column and bottom row
            # Use even larger line values to move labels further from axis
            if (j == 1) {
                title(ylab = "y", line = 2.2)  # Further increased space for y-axis label
            }
            if (i == n_datasets) {
                title(xlab = "x", line = 2.2)  # Further increased space for x-axis label
            }
        }
    }
    
    # No overall title as requested
    
    # Close graphics device if it was opened
    if (!is.null(plot_file)) {
        dev.off()
    }
    
    # Return layout information
    invisible(list(
        datasets = dataset_names,
        methods = methods,
        layout = c(n_datasets, n_methods)
    ))
}
