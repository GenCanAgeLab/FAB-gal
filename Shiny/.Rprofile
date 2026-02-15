source("renv/activate.R")

if (interactive() && requireNamespace("renv", quietly = TRUE)) {
  status <- renv::status()
  if (!status$synchronized) {
    message(
      "\nThis project uses renv.\n",
      "Required packages are missing.\n",
      "Run `renv::restore()` to install them.\n"
    )
  }
}
