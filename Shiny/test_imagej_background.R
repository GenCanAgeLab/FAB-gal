# Test ImageJ Subtract Background via rJava + EBImage
# Usage (in R): source('Shiny/test_imagej_background.R')

suppressPackageStartupMessages({
  library(EBImage)
  library(RBioFormats)
  library(rJava)
})

# ---- CONFIG ----
# Set this to the path of ij.jar bundled in the project
ij_jar_path <- 'ij/ij.jar'

if (!file.exists(ij_jar_path)) {
  stop("ij.jar not found at: ", ij_jar_path, "\nPlace ij.jar in the Shiny/ folder.")
}

# ---- INIT JAVA ----
.jinit(parameters = "-Xmx1g")
.jaddClassPath(ij_jar_path)

# ---- LOAD IMAGE (example) ----
# Use RBioFormats to read .czi and other formats
img_path <- 'Examples/Sample_cells_1.czi'

if (!file.exists(img_path)) {
  stop("Image not found at: ", img_path)
}

# Read via RBioFormats (returns array)
img_array <- read.image(img_path,normalize = F)

colorMode(img_array) <- 'Grayscale'

img_array <- img_array[,,1]
display(normalize(img_array))



# ---- EBImage -> ImageJ ByteProcessor ----
mat <- imageData(img_array)      # matrix (y, x) in [0-255]
h <- nrow(mat)
w <- ncol(mat)

message("Original matrix dims (h x w): ", h, " x ", w)
message("Original mat range: ", min(mat), " - ", max(mat))

# Save original to compare
orig_img <- Image(mat / 255, colormode = "Grayscale")
message("Original image stored in 'orig_img'")

# Try WITHOUT transpose - ByteProcessor may handle this differently
byte_array <- as.raw(as.vector(mat))
fp <- .jnew("ij.process.ByteProcessor", w, h, byte_array)

# ---- RUN SUBTRACT BACKGROUND ----
bs <- .jnew("ij.plugin.filter.BackgroundSubtracter")

radius <- 50.0
createBackground <- FALSE
lightBackground <- FALSE
useParaboloid <- FALSE
doPresmooth <- TRUE
correctCorners <- TRUE

bs$rollingBallBackground(fp, radius, createBackground,
                         lightBackground, useParaboloid,
                         doPresmooth, correctCorners)

# ---- EXTRACT RESULT ----
pixels <- fp$getPixels()
message("Output pixels class: ", class(pixels))
message("Output pixels length: ", length(pixels))

# Convert raw bytes to numeric
pixels_num <- as.numeric(pixels)
message("Output pixels range: ", min(pixels_num, na.rm=TRUE), " - ", max(pixels_num, na.rm=TRUE))

# Reconstruct: ByteProcessor stores as column-major (Fortran order)
out_mat <- matrix(pixels_num, nrow = h, ncol = w, byrow = FALSE)  # Column-major
bg_result <- Image(out_mat / 255, colormode = "Grayscale")
message("Background subtraction result stored in 'bg_result'")
message("Use display(bg_result) to view or display(orig_img) to compare")

