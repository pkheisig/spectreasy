export const launchCommand = 'spectreasy::spectreasy_gui()'
export const installAnalysisCommand = 'spectreasy::install_analysis_dependencies()'

export const installPackageCommand = `setup_packages <- c(
  "remotes",
  "httpuv",
  "later",
  "plumber",
  "jsonlite"
)

missing_packages <- setup_packages[
  !vapply(setup_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages)) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

remotes::install_github(
  "pkheisig/spectreasy",
  force = TRUE
)`

export const installAndLaunchCommand = `${installPackageCommand}

${launchCommand}`
