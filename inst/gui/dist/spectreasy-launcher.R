# Spectreasy portable setup helper
# This script installs nothing automatically. Run it in the R version that
# should host Spectreasy and the local cockpit backend.

spectreasy_launcher <- function() {
    minimum_r <- numeric_version("4.5.0")
    gui_dependencies <- c("httpuv", "later", "plumber", "jsonlite")
    print_install_code <- function() {
        cat('setup_packages <- c("remotes", "httpuv", "later", "plumber", "jsonlite")\n')
        cat('missing_packages <- setup_packages[!vapply(setup_packages, requireNamespace, logical(1), quietly = TRUE)]\n')
        cat('if (length(missing_packages)) {\n')
        cat('  install.packages(missing_packages, repos = "https://cloud.r-project.org")\n')
        cat('}\n')
        cat('remotes::install_github("pkheisig/spectreasy", force = TRUE)\n')
    }
    r_executable <- file.path(
        R.home("bin"),
        if (.Platform$OS.type == "windows") "R.exe" else "R"
    )
    r_executable <- normalizePath(r_executable, mustWork = FALSE)

    cat("\nSpectreasy local setup\n")
    cat(strrep("-", 52), "\n", sep = "")
    cat("R version   : ", R.version.string, "\n", sep = "")
    cat("R executable: ", r_executable, "\n", sep = "")
    cat("Library paths:\n")
    cat(paste0("  - ", .libPaths(), collapse = "\n"), "\n")
    cat(strrep("-", 52), "\n", sep = "")

    if (getRversion() < minimum_r) {
        cat("Spectreasy requires R 4.5.0 or newer.\n")
        cat("Download a current R release from https://cran.r-project.org/\n")
        return(invisible(FALSE))
    }

    package_path <- suppressWarnings(find.package("spectreasy", quiet = TRUE))
    if (!length(package_path) || !nzchar(package_path)) {
        cat("Spectreasy is not installed in this R library.\n\n")
        cat("Copy and run:\n\n")
        print_install_code()
        return(invisible(FALSE))
    }

    load_error <- tryCatch({
        loadNamespace("spectreasy")
        NULL
    }, error = function(error) conditionMessage(error))
    if (!is.null(load_error)) {
        cat("Spectreasy is installed but could not be loaded.\n")
        cat("Reason: ", load_error, "\n\n", sep = "")
        cat("Re-run the installation command above to repair dependencies.\n")
        return(invisible(FALSE))
    }

    missing_gui_dependencies <- gui_dependencies[
        !vapply(gui_dependencies, requireNamespace, logical(1), quietly = TRUE)
    ]
    if (length(missing_gui_dependencies)) {
        cat("Spectreasy is installed, but cockpit dependencies are missing: ")
        cat(paste(missing_gui_dependencies, collapse = ", "), "\n\n", sep = "")
        cat("Copy and run:\n\n")
        cat('install.packages(c(')
        cat(paste(sprintf('"%s"', missing_gui_dependencies), collapse = ", "))
        cat('), repos = "https://cloud.r-project.org")\n')
        return(invisible(FALSE))
    }

    cat("Spectreasy : ", as.character(utils::packageVersion("spectreasy")), "\n", sep = "")
    cat("Status      : ready\n")
    if (!interactive()) {
        cat("\nRun spectreasy::spectreasy_gui() in an interactive R session.\n")
        return(invisible(TRUE))
    }

    answer <- trimws(tolower(readline("Launch the Spectreasy cockpit now? [Y/n]: ")))
    if (!answer %in% c("n", "no")) {
        spectreasy::spectreasy_gui()
    }
    invisible(TRUE)
}

spectreasy_launcher()
