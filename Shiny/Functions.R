
# Read image ----
  myreadimg <- function(imgpath){
    img_obj <- read.image(imgpath,
                      proprietary.metadata = F,
                      read.metadata = T,
                      normalize = F)
    colorMode(img_obj) <- 'Grayscale'
    img_obj
  }

# Safely extract a channel from an image object----
  get_channel <- function(imgobj, channel){
    ch <- as.numeric(channel)
    if (length(ch) == 0) {
      stop("channel cannot be NULL")
    }
    if (is.na(ch)) {
      stop("A channel must be provided.")
    }
    ndim <- length(dim(imgobj))
    if (ndim == 2){
      if (ch == 1){
        return(imgobj)
      } else {
        stop("Trying to get a channel from an image with less dimensions")
      }
    }
    if (ndim > 2){
      if (dim(imgobj)[3] >= ch){
        return(imgobj[,,ch])
      }
      else {
        stop("Channel not present in image")
      }
    }
  }

# Extract bit depth from image metadata ----
  get_bitdepth <- function(img_obj){
    tryCatch({
      2^coreMetadata(img_obj)$bitsPerPixel 
    },
    error = function(e) {
      print("image bit depth not available")
      NULL}
    )
  }
  
# Get pixel resolution----
  get_pxarea <- function(img_obj){
    pxsize <- c(globalMetadata(img_obj)$XResolution, globalMetadata(img_obj)$YResolution)
    oriunit <- tolower(globalMetadata(img_obj)$ResolutionUnit)
    resunit <- NA
    if (oriunit %in% c("cm","centimeter")){
      pxsize <- pxsize / 10000
      resunit <- "micron^2"
    }
    if (oriunit %in% c("um","micron","micrometer")){
      resunit <- "micron^2"
    }
    if (!is.na(resunit)){
      return( 1/prod(pxsize))
    } else {
      return(NULL)
    }
    }

# Subtract background using ImageJ rolling ball algorithm----
  subtract_background <- function(ebimg, radius = 50.0) {
    
    # Ensure radius is double for Java (not integer)
    radius <- as.double(radius)
    
    # ---- VALIDATE INPUT ----
    if (!inherits(ebimg, c("Image", "array", "matrix"))) {
      stop("Input must be an EBImage Image object or numeric array/matrix")
    }
    
    # Extract data if EBImage object
    if (inherits(ebimg, "Image")) {
      mat <- EBImage::imageData(ebimg)
    } else {
      mat <- ebimg
    }
    
    # Ensure 2D
    if (!is.matrix(mat) && !is.array(mat)) {
      stop("Input must be 2D (matrix or 2D array)")
    }
    if (length(dim(mat)) != 2) {
      stop("Input must be 2D")
    }
    
    # Get dimensions BEFORE any type conversion
    h <- nrow(mat)
    w <- ncol(mat)
    
    # ---- PREPARE DATA ----
    # Ensure 0-255 range for ByteProcessor
    if (max(mat) <= 1.0) {
      # Normalize from [0, 1] to [0, 255]
      mat <- mat * 255
    }
    # Convert to integer while preserving matrix shape
    mat <- matrix(as.integer(mat), nrow = h, ncol = w)
    
    # ---- CREATE JAVA OBJECTS (exactly as in test_imagej_background.R) ----
    byte_array <- as.raw(as.vector(mat))
    fp <- .jnew("ij.process.ByteProcessor", w, h, byte_array)
    
    # ---- RUN SUBTRACT BACKGROUND ----
    bs <- .jnew("ij.plugin.filter.BackgroundSubtracter")
    
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
    pixels_num <- as.numeric(pixels)
    
    # Reconstruct: column-major (byrow=FALSE)
    out_mat <- matrix(pixels_num, nrow = h, ncol = w, byrow = FALSE)
    
    # ---- RETURN AS EBIMAGE ----
    result <- EBImage::Image(out_mat, colormode = "Grayscale")
    return(result)
  }

# Helper: Check if Java initialized ----
  .jniInitialized <- function() {
    tryCatch({
      .jcall("java/lang/System", "S", "getProperty", "java.version")
      TRUE
    }, error = function(e) {
      FALSE
    })
  }
  
  


# Function to report pixel intensity ----
  calc_bgalstats <- function(img_obj,imgth,coords) {
    perA <- sum(imgth == 1) * 100 / length(imgth)
    mfi <- mean(img_obj)
    if (!is.null(coords)){
      cint <- tryCatch({
        img_obj[coords[1],coords[2]]
      }, error= function(e) {return(NA)}
      )
    } else {cint <- NA}
    sprintf("%d  MFI = %.2f  Sel_Area = %.1f%%",cint,mfi,perA)
  }    
  
  
#  Functions for measure and run.all ----
  
  # Emtpy results to store batch results
  empty_res <- function(){
    data.frame(
      File = character(0), 
      NImagePx = integer(0),
      Area = double(0),
      NposPx = integer(0),
      RawIntDen = integer(0),
      Bckg_MFI = integer(0),
      CTFperAreaPx = double(0),
      CTFperAreaUnits = double(0),
      Nuclei_Count = integer(0),
      CTFperNuclei = double(0),
      status = character(0)
    )
  }
  
# Function to process a single image
process_single_image <- function(img_path, appsets,outdir=NULL) {
  errors <- NULL
  # Load image
  img_obj <- tryCatch(
    {myreadimg(img_path)},
    error = function(e) {
      msg = "Error reading file"
      showNotification(paste(msg,":",basename(img_path)),type='error')
      errors <<- c(errors,msg)
      return(NULL)
    }
  )
  # stop if reading error
  if (is.null(img_obj)){
    return(c(as.list(rep(NA,9)),as.list(errors)))
  }
  # Process nuclei
  nres <- tryCatch({
    if (appsets$nchan == "") {
      list("Count"=NA,"Maskn"=NULL)
    } else {
      process_nuclei(get_channel(img_obj, appsets$nchan), appsets)
    }
    }, error = function(e) {
      msg <- paste("Error processing nuclei",e)
      showNotification(paste(msg, ":",basename(img_path)),type='error')
      errors <<- c(errors, msg)
      list("Count"=NA,"Maskn"=NULL)
    })
  # Save nuclei mask if requested and present
  outfname <- file.path(outdir,gsub("\\..*$",".png",basename(img_path)))
  if (!is.null(nres$Maskn)){
    tryCatch({
      if (file.exists(outfname)){stop("File already exists")}
      writeImage(nres$Maskn,outfname)
      }, error = function(e){
      msg <- paste("Error saving nuclei mask. File already exits?")
      showNotification(paste(msg, ":",basename(img_path)),type='error')
      errors <<- c(errors, msg)
    })
  }

  # Process sabgal
  sres <- tryCatch({
    process_sabgal(get_channel(img_obj, appsets$schan), appsets)
    }, error = function(e) {
      msg <- paste("Error processing sabgal", e$message)
      showNotification(paste(msg, ":",basename(img_path)),type='error')
      errors <<- c(errors, msg)
      return(as.list(rep(NA,7)))
    })
  # Collect all errors
  status <- if (!is.null(errors)) {
    paste(errors, collapse = "; ")
  } else {
    "Success"
  }
  res <- c(sres,
           nres$Count,
           sres[[4]]/nres$Count,
           as.list(status))
  return(res)
}
  
  
  # Processing nuclei
  process_nuclei <- function(img_obj, appsets){
    # Run Thresholding
    img_obj <- run_thresh(img_obj, appsets)
    # Check if any objects detected
    if (max(img_obj) == 0) {return(list("Count"=0,"Maskn"=NULL))}
    # Apply segmentation
    # img_obj <- bwlabel(img_obj)
    img_obj <- watershed(distmap(img_obj))
    # Remove objects on edges if activated
    if (appsets$rmborder == TRUE) {
      onedge <- unique(c(
        img_obj[1,],
        img_obj[nrow(img_obj),],
        img_obj[,1],
        img_obj[,ncol(img_obj)])
      )
      img_obj <- rmObjects(img_obj,onedge , reenumerate = F)
    }
    # Compute features and filter
    # feat <- computeFeatures.shape(img_obj)[,'s.area']
    feat <- c(table(c(img_obj),exclude = 0)) # Way faster
    # Filter by area using appsets$nsize
    feat.f <- feat[feat >= appsets$nsize[1] & feat <= appsets$nsize[2]]
    # Return nuclei mask if enabled
    maskn <- if (appsets$saven == TRUE){
      rmlabs <- setdiff(names(feat),names(feat.f))
      colorLabels(rmObjects(img_obj,rmlabs))
    } else {NULL}
    return(list("Count"=length(feat.f),"Maskn"=maskn))
  }

  # Run thresholding for nuclei
  run_thresh <- function(img_obj, input){
    # Gaussian blur
    if (input$sigma > 0){
      img_obj <- gblur(img_obj, input$sigma)
    }
    # Background correction
    if (input$sback == TRUE){
      img_obj <- subtract_background(img_obj, input$radius)
    }
    # Apply Threshold
    img_obj <- img_obj >= input$thres.n[1] & img_obj <= input$thres.n[2]
    # Fill holes
    result <- fillHull(img_obj)
    return(result)
  }


# Process bgal
  process_sabgal <- function(img_obj, appsets){
    img_th <- img_obj >= appsets$thres.s[1] & img_obj <= appsets$thres.s[2]
    NimagePx <- length(img_obj)
    Area <- NimagePx*as.numeric(appsets$pxarea)
    NpxPos <- sum(img_th == 1) 
    RawIntDen <- sum(img_obj[img_th == 1])
    CTF <- RawIntDen - (as.numeric(appsets$bmfi)*NpxPos)
    CTFpx <- CTF / NpxPos
    CTFpa <- CTF / Area
    return(list(
      NimagePx,
      Area,
      NpxPos,
      RawIntDen,
      as.numeric(appsets$bmfi),
      CTFpx,
      CTFpa)
      )
  }
  

# Function to process all images in batch
  process_all_images <- function(image_files, input, outdir) {
    n_files <- length(image_files)
    resdf <- empty_res()
    withProgress(message = "Processing images", value = 0, {
      # Process each image
      for (i in seq_along(image_files)) {
        res <- process_single_image(image_files[i], input, outdir)
        resdf[i,] <- c(list(basename(image_files[i])),res)
        incProgress(1/n_files, detail = paste("Doing image", i, "of", n_files))
      }
    resdf
    })
  }
  
