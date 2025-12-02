#generate report html file

# 1. Check and install required packages
required_packages <- c("rmarkdown", "knitr", "rprojroot")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(rprojroot)

# 2. Locate project root directory
root_result <- try(find_root(has_file("Metadoon_Report.Rmd")), silent = TRUE)

if (inherits(root_result, "try-error")) {
  message("Warning: Could not identify project root using rprojroot. Using current working directory.")
  root_dir <- getwd()
} else {
  root_dir <- root_result
}

# 3. Set Working Directory
setwd(root_dir)
message("📂 Working directory set to: ", root_dir)

# 4. Define filenames
rmd_file_name <- "Metadoon_Report.Rmd"
# We define the output name here, but we will also use the string directly later
# to avoid errors if the Rmd script wipes the environment variables.
html_file_name <- "Metadoon_Report.html"

# 5. Render the Report
message("⏳ Starting report rendering...")

rmarkdown::render(
  input = rmd_file_name, 
  output_file = html_file_name,
  output_dir = root_dir, 
  output_format = "html_document",
  quiet = FALSE
)

# 6. Finalize
# FIX: We use the string "Metadoon_Report.html" directly here.
# This prevents the "object not found" error if the Rmd file contains 'rm(list=ls())'
# which wipes variables from memory during execution.
final_output_name <- "Metadoon_Report.html"
full_output_path <- file.path(getwd(), final_output_name)

if (file.exists(full_output_path)) {
  cat("\n✅ HTML report successfully generated:\n", full_output_path, "\n")
  
  # Attempt to open the default browser
  tryCatch({
    utils::browseURL(full_output_path)
  }, error = function(e) {
    message("Report generated! (Could not open browser automatically)")
  })
  
} else {
  stop("❌ Error: The process finished, but the HTML file was not found.")
}