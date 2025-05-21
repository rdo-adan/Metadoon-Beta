# generate_report.R

if (!requireNamespace("rmarkdown", quietly = TRUE)) install.packages("rmarkdown")
if (!requireNamespace("knitr", quietly = TRUE)) install.packages("knitr")

output_file <- "Metadoon_Report.html"
rmd_file <- "Metadoon_Report.Rmd"

# Ensure the working directory is the same as the location of this script
setwd(dirname(normalizePath(rmd_file)))


# Render HTML from RMarkdown
rmarkdown::render(rmd_file, output_file = output_file, output_format = "html_document")

# Notify user
cat("\n✅ HTML report successfully generated: ", output_file, "\n")

# Open in default browser (cross-platform)
utils::browseURL(output_file)