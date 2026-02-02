
# Read image ----
  myreadimg <- function(imgpath){
    img <- read.image(imgpath,
                      proprietary.metadata = F,
                      read.metadata = T,
                      normalize = F)
    colorMode(img) <- 'Grayscale'
    img
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
  get_bitdepth <- function(img){
    tryCatch({
      2^coreMetadata(img)$bitsPerPixel 
    },
    error = function(e) {
      print("image bit depth not available")
      NULL}
    )
  }
  
# Get pixel resolution----
  get_pxarea <- function(img){
    pxsize <- c(globalMetadata(img)$XResolution, globalMetadata(img)$YResolution)
    oriunit <- tolower(globalMetadata(img)$ResolutionUnit)
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

# Subtract background----
  subtract_background <- function(img, radius){
    kernel <- ball_kernel(radius, normalize = TRUE)
    background <- filter2(img, kernel)
    img <- img - background
    img[img < 0] <- 0
    return(img)
  }
  
  # Function to create a 2D ball kernel
  # (matching scikit-image implementation)
  ball_kernel <- function(radius, normalize = TRUE) {
    # Create coordinate matrices - using ceiling of radius like Python version
    size <- 2 * ceiling(radius) + 1
    center <- ceiling(radius) + 1
    
    # Create meshgrid coordinates
    coords <- expand.grid(
      x = seq(-ceiling(radius), ceiling(radius)),
      y = seq(-ceiling(radius), ceiling(radius))
    )
    
    # Calculate sum of squares and distance
    sum_of_squares <- coords$x^2 + coords$y^2
    distance_from_center <- sqrt(sum_of_squares)
    
    # Create kernel: z = sqrt(r^2 - (x^2 + y^2)) where inside sphere
    kernel_values <- sqrt(pmax(radius^2 - sum_of_squares, 0))
    
    # Set values outside radius to 0 (more practical than Inf for our use case)
    kernel_values[distance_from_center > radius] <- 0
    
    # Convert to matrix
    kernel <- matrix(kernel_values, nrow = size, ncol = size, byrow = TRUE)
    
    # Normalize kernel if requested
    if (normalize) {
      kernel <- kernel / sum(kernel)
    }
    
    return(kernel)
  }
  


# Function to report pixel intensity ----
  calc_bgalstats <- function(img,imgth,coords) {
    perA <- sum(imgth == 1) * 100 / length(imgth)
    mfi <- mean(img)
    if (!is.null(coords)){
      cint <- tryCatch({
        img[coords[1],coords[2]]
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
process_single_image <- function(img_path, appsets) {
  errors <- NULL
  # Load image
  img <- tryCatch(
    {myreadimg(img_path)},
    error = function(e) {
      msg = "Error reading file"
      showNotification(paste(msg,":",basename(img_path)),type='error')
      errors <<- c(errors,msg)
    }
  )
  # stop if reading error
  if (is.null(img)){
    return(c(as.list(rep(NA,10)),as.list(errors)))
  }
  # Process nuclei
  nres <- tryCatch({
    if (appsets$nchan == "") {
      NA }
    else {
      process_nuclei(get_channel(img, appsets$nchan), appsets)
      }
    }, error = function(e) {
      msg <- paste("Error processing nuclei")
      showNotification(paste(msg, ":",basename(img_path)),type='error')
      errors <<- c(errors, msg)
      return(NA)
    })
  # Process sabgal
  sres <- tryCatch({
    process_sabgal(get_channel(img, appsets$schan), appsets)
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
           nres,
           sres[[4]]/nres,
           as.list(status))
  print(res)
  return(res)
}
  
  
  # Processing nuclei
  process_nuclei <- function(img, appsets){
    # Run Thresholding
    img <- run_thresh(img, appsets)
    # Check if any objects detected
    if (max(img) == 0) {return(0)}
    # Apply segmentation
    # img <- bwlabel(img)
    img <- watershed(distmap(img))
    # Remove objects on edges if activated
    if (appsets$rmborder == TRUE) {
      onedge <- unique(c(
        img[1,],
        img[nrow(img),],
        img[,1],
        img[,ncol(img)])
      )
      img <- rmObjects(img,onedge , reenumerate = F)
    }
    # Compute features and filter
    # feat <- computeFeatures.shape(img)[,'s.area']
    feat <- c(table(c(img),exclude = 0)) # Way faster
    # Filter by area using appsets$nsize
    feat <- feat[feat >= appsets$nsize[1] & feat <= appsets$nsize[2]]
    return(length(feat))
  }

  # Run thresholding for nuclei
  run_thresh <- function(img, input){
    # Gaussian blur
    if (input$sigma > 0){
      img <- gblur(img, input$sigma)
    }
    # Background correction
    if (input$sback == TRUE){
      img <- subtract_background(img, input$radius)
    }
    # Apply Threshold
    img <- img >= input$thres.n[1] & img <= input$thres.n[2]
    # Fill holes
    result <- fillHull(img)
    return(result)
  }


# Process bgal
  process_sabgal <- function(img, appsets){
    img_th <- img >= appsets$thres.s[1] & img <= appsets$thres.s[2]
    NimagePx <- length(img)
    Area <- NimagePx*as.numeric(appsets$pxarea)
    NpxPos <- sum(img_th == 1) 
    RawIntDen <- sum(img[img_th == 1])
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
  process_all_images <- function(image_files, input) {
    n_files <- length(image_files)
    resdf <- empty_res()
    withProgress(message = "Processing images", value = 0, {
      # Process each image
      for (i in seq_along(image_files)) {
        res <- process_single_image(image_files[i], input)
        resdf[i,] <- c(list(basename(image_files[i])),res)
        incProgress(1/n_files, detail = paste("Doing image", i, "of", n_files))
      }
    resdf
    })
  }
  
