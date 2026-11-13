# data-raw/dental.R
#
# One-time script to build data/dental.rda from data-raw/dentex.csv.
# Run this locally (from the package root, e.g. with devtools::load_all()
# or just Rscript data-raw/dental.R) whenever the source CSV changes.
# It requires the 'usethis' package (install.packages("usethis") if needed).

dental <- read.csv("data-raw/dentex.csv", header = TRUE, stringsAsFactors = FALSE)

# ID and momid are identifiers / cluster labels, not modeling factors -
# leave as character/integer. All other 0/1 variables are stored as
# integer 0/1 (matching the coding used throughout the gmediate/gmedsens
# examples and vignette, where these are treated as 0/1 numeric or
# binomial-family responses/predictors).
dental$ID    <- as.character(dental$ID)
dental$momid <- as.character(dental$momid)

usethis::use_data(dental, overwrite = TRUE)
