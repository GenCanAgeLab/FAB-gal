# Load libraries and functions ----
library(shiny)
library(shinyFiles)
library(bslib)
require(RBioFormats)
require(EBImage)
require(rJava)

# Initialize Java for ImageJ background subtraction
.jinit(parameters = "-Xmx1g")
ij_jar_path <- 'ij/ij.jar'
if (file.exists(ij_jar_path)) {
  .jaddClassPath(ij_jar_path)
} else {
  warning("ij.jar not found at:", ij_jar_path, "\nBackground subtraction may not work.")
}

# Source functions
source("Functions.R")

# Define UI ----
ui <- fluidPage(
  ## Application title ----
  titlePanel(
    windowTitle = "FAβGal - SA-β-Galactosidase Analysis",
    fluidRow(
      style = "text-align: center;",
      column(
        width = 3,
        a(
          img(
            src = "Logo.png",
            align = "rigth",
            height = 80
          ),
          href="https://github.com/antotartier/FAB-gal",
          target="_blank"
        )
      ),
      column(
        width = 6,
        h3("FAβGal - SA-β-Galactosidase Analysis"),
        h4(
          a(
            "Genomes, Cancer and Aging (GenCanAge) Group",
            href = "https://github.com/GenCanAgeLab",
            target = "_blank"
          )
        )
      ),
      column(
        width = 3,
        a(
          img(
          src = "Logo_gencanage.jpeg",
          align = "rigth",
          height = 80
          ),
          href="https://github.com/GenCanAgeLab",
          target="_blank"
        )
        )
    )
  ),
  ## SidebarPanel ----
  sidebarPanel(
    width = 3,
    style = "background: #ffffff;text-align: center;",
    ### Channel selection ----
    wellPanel(
      style = "background: #cce6d0",
      div(strong("Channel Selection")),
      fluidRow(
        column(
          width = 6,
          align = "center",
          strong("SAbGal")
        ),
        column(
          width = 6,
          selectizeInput(
            "schan",
            label = NULL,
            choices = NULL,
          )
        )
      ),
      fluidRow(
        column(
          width = 6,
          strong("Nuclei")
        ),
        column(
          width = 6,
          selectizeInput(
            "nchan",
            label = NULL,
            choices = NULL,
            options = list(
              allowEmptyOption = TRUE,
              showEmptyOptionInDropdown = TRUE,
              emptyOptionLabel = "None"
            )
          )
        )
      )
    ),
    ### Nuclei options wellpanel ----
    wellPanel(
      style = "background: #eddcd3",
      align = "center",
      div(strong("Nuclei segmentation"), style = "text-align:center"),
      br(),
      # sliderInput(
      #   "sigma",
      #   label = "Gaussian blur",
      #   min = 0,
      #   max = 10,
      #   value = 0,
      #   step = 0.1
      # ),
      checkboxInput(
        'sback',
        label = 'Substract background',
        value = FALSE
      ),
      sliderInput(
        "radius",
        label = "Rolling ball radius",
        min = 21,
        max = 1001,
        value = 51,
        step = 1
      ),
      actionButton('ApplyTh', "Count and filter", class = "btn-primary"),
      checkboxInput('rmborder', "Remove on edges", FALSE),
      sliderInput('nsize', 'Area (in pixels)', 0, 10000, c(0, 10000))
    ),
    ### Pixel size, MFI and output dir panel ----
    wellPanel(
      style = "background: #dccce6",
      strong('Pixel area (microns^2)'),
      fluidRow(
        column(
          width = 7,
          align = 'center',
          textInput('pxarea', NULL, "", placeholder = "Ej. 0.45")
        ),
        column(
          width = 5,
          align = 'center',
          actionButton('getpxa', 'Get')
        )
      ),
      strong('Background fluorescence'),
      fluidRow(
        column(
          width = 7,
          align = 'center',
          textInput('bmfi', NULL, "", placeholder = "Ej. 23")
        ),
        column(
          width = 5,
          align = 'center',
          actionButton('getbmfi', 'Get')
        )
      ),
      checkboxInput("saven","Save nuclei segmentation",value = F),
      textOutput("outdirtext"),
      shinyDirButton(
        id = 'outdir',
        label = 'Browse',
        title = 'output dir'
      )
    )
  ),
  ## Main panel ----
  mainPanel(
    width = 9,
    style = "text-align: center;",
    ### File selection section ----
    wellPanel(
      style = "background: #eddcd3;",
      fluidRow(
        column(
          width = 3,
          div(strong("Choose a folder")),
          div(shinyDirButton(
            id = 'dir',
            label = 'Browse',
            title = 'dire'
          ))
        ),
        column(
          width = 6,
          strong("Choose a file"),
          selectizeInput(
            'imgpath',
            label = NULL,
            selected = NULL,
            choices = NULL,
            width = "100%"
          )
        ),
        column(
          width = 3,
          strong("Filter files (Regex)"),
          textInput(
            'file_filter',
            label = NULL,
            placeholder = "e.g., .czi or .tif",
            width = "100%"
          )
        )
      ),
      fluidRow(
        actionButton('run', "Measure current", class = "btn-primary"),
        actionButton('run.all', "Run for all images", class = "btn-primary"),
        actionButton("help", "Click Here To Learn More", class = "btn btn-info")
      )
    ),
    ### Images and sliders ----
    fluidRow(
      #### Nuclei ----
      column(
        width = 6,
        style = "padding-right: 2px",
        wellPanel(
          sliderInput(
            inputId = "thres.n",
            label = "Nuclei Threshold",
            min = 0,
            max = 255,
            value = c(0, 255),
            round = 2
          ),
          actionButton('runotsu.n', 'Auto'),
          actionButton('reset.thn', 'Reset'),
          plotOutput(
            outputId = "nplot",
            brush = brushOpts(
              id = "nplot_brush",
              direction = "xy",
              delayType = "debounce",
              delay = 500,
              clip = TRUE,
              resetOnNew = TRUE
            )
          ),
          actionButton('reset.zoom.n', 'Reset Zoom')
        ),
        verbatimTextOutput('stats.n')
      ),
      column(
        ### SABGAL ----
        width = 6,
        style = "padding-left: 2px",
        wellPanel(
          sliderInput(
            inputId = "thres.s",
            label = "SA-β-Gal Threshold",
            min = 0,
            max = 255,
            value = c(0, 255),
            round = 2
          ),
          actionButton('runotsu.s', 'Auto'),
          actionButton('reset.ths', 'Reset'),
          plotOutput(
            outputId = "splot",
            brush = brushOpts(
              id = "splot_brush",
              direction = "xy",
              delayType = "debounce",
              delay = 500,
              clip = TRUE,
              resetOnNew = TRUE
            )
          ),
          actionButton('reset.zoom.s', 'Reset Zoom')
        ),
        verbatimTextOutput('stats.s')
      )
    ),
    ## Results table ----
    fluidRow(
      style = "align-text:center",
      h4("Batch Analysis Results"),
      downloadButton("download_results", "Download Results CSV"),
      actionButton('reset', 'Reset results')
    ),
    fluidRow(
      tableOutput("results_table"),
      br()
    )
  ),
  ## Footer ----
  fluidRow(
    style = "text-align: center;",
    p(
      "Developed by Antonio Garcia-Bernardo Tartiere, Alejandro P. Ugalde and David Roiz del Valle"
    ),
    a(
      "github.com/GenCanAgeLab/SAbGal_quant",
      href = "https://github.com/antotartier/FAB-gal",
      target = "_blank"
    )
  ),
  br(),
  fluidRow(
    style = "text-align: center;",
    a(
      img(
      src = "https://www.uniovi.es/documents/39158/11ff14cf-90c5-892e-473f-829945ed1733",
      align = "rigth",
      height = 80
      ),
      href="https://www.uniovi.es/"
    ),
    a(
      img(
        src = "https://ispa-finba.es/wp-content/uploads/2019/12/ISPA_marca_principal_color-150x150.png",
        align = "rigth",
        height = 80
      ),
      href="https://ispa-finba.es/"
    ),
    a(
      img(
        src = "Logo_iuopa.jpg",
        align = "rigth",
        height = 80
      ),
      href="https://www.unioviedo.es/IUOPA/"
    ),
    a(
      img(
        src = "  https://www.aei.gob.es/sites/default/files/inline-images/miciu-cofinanciadoUE-aei.png",
        align = "rigth",
        height = 80
      ),
      href="https://www.aei.gob.es/"
    )
  ),
  
  ## CSS Styling ----
  tags$head(
    tags$style(
      '
      .well {padding:2px; margin-bottom: 2px;}
      .form-group {margin-bottom: 2px;}
      .checkbox label {font-weight: bold;}
      .shiny-notification {overflow-wrap: break-word;}
      '
    )
  )
)


# Define Server ----
server <- function(input, output, session) {
  ## Directory handling ----

  roots <- c("Examples" = './Examples', 'Home' = fs::path_home(), getVolumes()())
  # roots <- getVolumes()()
  
  ### Input directory ----
  ## Client connection for file system
  shinyDirChoose(input, 'dir', roots = roots, session = session,defaultRoot = "Home")
  
  ## Handeling selected dir
  ## It also prevent updating if the browse dialog is cancelled
  current_mydir <- reactiveVal(NULL)
  observeEvent(input$dir, {
    parsed <- tryCatch(
      parseDirPath(roots = roots, input$dir),
      error = function(e) character(0)
    )
    if (is.character(parsed) && length(parsed) > 0 && nzchar(parsed)) {
      current_mydir(parsed)
    }
  })
  mydir <- reactive({
    validate(need(current_mydir(),"Please select a folder"))
    current_mydir()
  })

  ### File handeling ----
  ## Reactive for filtered files
  filtered_files <- reactive({
    validate(need(mydir(), "open a folder"))
    all_files <- list.files(mydir(), pattern = "\\.", full.names = FALSE)
    # Apply file filter if provided
    if (!is.null(input$file_filter) && nzchar(input$file_filter)) {
      all_files <- all_files[grepl(
        input$file_filter,
        all_files,
        ignore.case = TRUE
      )]
    }
    all_files
  })

  ## Observer to update file selectizeinput
  observe({
    updateSelectizeInput(
      session,
      'imgpath',
      choices = filtered_files(),
      selected = "" # Don't pre-select any file
    )
  })

  ## Manage selected image 
  # a reactiveVal is used to allow reset it when open a folder)
  imgpath <- reactiveVal(NULL)
  
  observeEvent(input$imgpath, {
    imgpath(input$imgpath)
  })

  ### Actions when changing directory ----
  observeEvent(mydir(), {
    updatechan(TRUE)
    imgpath(NULL)
    updateSliderInput(session, 'thres.n', value = c(0, 255), min = 0, max = 255)
    updateSliderInput(session, 'thres.s', value = c(0, 255), min = 0, max = 255)
    updateTextInput(session, 'bmfi', value = "")
    updateTextInput(session, 'pxarea', value = "")
    outdir(mydir()) # Set output dir to current dir
  })
  
  ### Output directory ----
  shinyDirChoose(input, 'outdir', roots = roots, session = session, 
                 defaultRoot = "Home")
  
  # A reactiveVal is used to allow setting it by different means
  outdir <- reactiveVal(NULL)
  observe({
    outdir(parseDirPath(roots, input$outdir))
  })
  # Notify output dir when save nuclei is enabled
  observeEvent(input$saven,{
    req(outdir())
    if (input$saven == TRUE){
      showNotification("Output directory will be:",outdir(),type="message")
    }
  })
  
  observeEvent(outdir(),{
    if (input$saven == TRUE){
      showNotification("Output directory will be:",outdir(),type="message")
    }
  })
  
  
  
  ## Image loading ----

  ### Read image ----
  valid_img <- reactive({
    validate(need(mydir(),"Please select a folder"))
    validate(need(imgpath(),"Please select an image"))
    
    # Read image
    imgobj <- tryCatch(
      {
        myreadimg(file.path(mydir(), imgpath()))
      },
      error = function(e) {
        showNotification(
          paste0("Unknown file format,cannot read image file: ", imgpath()),
          type = 'error'
        )
        return(NULL)
      }
    )
    if (is.null(imgobj)){return(NULL)}
    if (length(dim(imgobj)) > 3) {
      showNotification(
        "Image with more than 3 dimesions. Only images with XYC dimensions are supported",
        type = 'error'
      )
      return(NULL)
    }
    if (!grepl('^XYC', coreMetadata(imgobj)$dimensionOrder)) {
      showNotification(
        paste0(
          "Image dimesion order should be XYC. Your image is: coreMetadata(imgobj)$dimensionOrder"
        ),
        type = 'error'
      )
      return(NULL)
    }
    imgobj
  })


  ### Number of channels ----
  nc <- reactive({
    req(valid_img())
    coreMetadata(valid_img())$sizeC
  })
  
  ### Image bitdepth ----
  img_bitdepth <- reactive({
    req(valid_img())
    get_bitdepth(valid_img())
  })

  ### Valid image actions ----
  
  # Logic to update channels only when the first image is open
  # or when the image has less channels that currently set
  updatechan <- reactiveVal(TRUE)
  
  observeEvent(valid_img(),{
    req(valid_img(), nc())
    # Actual channel selection (it there is one)
    n <- input$nchan
    s <- input$schan
    maxchan <- max(as.numeric(c(n,s, 1)), na.rm = T)
    # update selectize inputs if required
    if (updatechan() == TRUE | nc() < maxchan) {
      if (nc() >= 2) {
        n <- 1
        s <- 2
      } else {
        n <- ""
        s <- 1
      }
      updatechan(FALSE)
    }
    # If image has only when channel, set it to sabgal
    if (nc() == 1){
      n <- ""
      s <- 1
    }
    updateSelectizeInput(session,'nchan',choices = 1:nc(),selected = n)
    updateSelectizeInput(session,'schan',choices = 1:nc(),selected = s)
    # Reset zoom when loading a new image
    crop_region(NULL)
  })

  
  ### Pixel physical size ----

  px_area <- reactiveVal(NA)

  #### Observer for text input
  observeEvent(input$getpxa, {
    req(valid_img())
    pxa <- tryCatch(
      {
        get_pxarea(valid_img())
      },
      error = function(e) {
        showNotification(
          paste0('Error: ', e),
          type = 'error'
        )
      }
    )
    if (is.null(pxa)) {
      px_area(NA)
      return(NULL)
      }
    if (pxa$pxa == ""){
      showNotification(pxa$error,type='error')
      px_area(NA)
    } else {
      if (!is.null(pxa$error)){
        showNotification(pxa$error,type='warning')
      }
      px_area(pxa$pxa)
      updateTextInput(session, 'pxarea', value = pxa$pxa)
    }
  })
  
  #### Observer for get button
  observe({
    if (grepl('[^0-9\\.]', input$pxarea)) {
      showNotification(
        "pixel area must be a numeric value. Ej 0.02",
        type = 'error'
      )
      updateTextInput(session, 'pxarea', value = "")
      px_area(NA)
    } else {
      px_area(as.numeric(input$pxarea))
    }
  })

  ### Background fluorescence ----

  bmfi <- reactiveVal(NA)
  ### Observer for text input
  observeEvent(input$getbmfi, {
    req(sabgal(),sabgal_stats())
    bmfi <- sabgal_stats()$mfi
    if (!is.null(bmfi)) {
      updateTextInput(session, 'bmfi', value = bmfi)
    }
  })
  ### Observer for get button
  observe({
    if (grepl('[^0-9\\.]', input$bmfi)) {
      showNotification(
        "Background fluorescence must be a positive numeric value. Ej 23.2",
        type = 'error'
      )
      updateTextInput(session, 'bmfi', value = "")
      bmfi(NA)
    } else {
      bmfi(as.numeric(input$bmfi))
    }
  })

  
  observe({print("pxarea");print(input$pxarea);print("bmfi");print(input$bmfi)})
  
  ## Nuclei  processing ----

  ### Nuclei Preprocessing ----
  nuclei <- reactive({
    req(valid_img())
    validate(need(input$nchan, "Select a channel"))
    tryCatch(
      get_channel(valid_img(), input$nchan),
      error = function(e) {
        return(NULL)
      }
    )
  })

  # Reactive to store nuclei blurred
  # nuclei_b <- reactive({
  #   req(nuclei())
  #   if (input$sigma > 0) {
  #     gblur(nuclei(), input$sigma)
  #   } else {
  #     nuclei()
  #   }
  # })

  # Reactive to store nuclei background subtracted
  nuclei_bgs <- reactive({
    # req(nuclei_b())
    req(nuclei())
    if (input$sback == TRUE) {
      tryCatch({
        # subtract_background(nuclei_b(), input$radius)
        subtract_background(nuclei(), input$radius)
      }, error = function(e) {
        showNotification(paste("Error in background subtraction:", e$message), type = "error")
        # nuclei_b()
        nuclei()
      })
    } else {
      # nuclei_b()
      nuclei()
    }
  })

  # Reactive to store nuclei thresholded
  nuclei_th <- reactive({
    req(nuclei_bgs())
    imgobj <- nuclei_bgs()
    imgobj <- imgobj >= input$thres.n[1] & imgobj <= input$thres.n[2]
    fillHull(imgobj)
  })

  # ObserveEvent for otsu
  observeEvent(input$runotsu.n, {
    deft <- round(otsu(nuclei_bgs(), c(0, img_bitdepth() - 1), img_bitdepth()))
    updateSliderInput(session, 'thres.n', value = c(deft, img_bitdepth() - 2))
  })

  # ObserveEvent for reset
  observeEvent(input$reset.thn, {
    updateSliderInput(session, 'thres.n', value = c(0, img_bitdepth() - 1))
  })


  ### Nuclei Segmentatiobn ----
 
  # Handle count and filter status
  apply_th <- reactiveVal(FALSE)

  # Reactive expression to generate the final nuclei image
  nuclei_seg <- reactiveVal(NULL)
  nuclei_objects <- reactiveVal(NULL)

  # Segment when clicking apply
  observeEvent(input$ApplyTh, {
    req(nuclei_th())
    if (all(input$thres.n == c(0, img_bitdepth() - 1))) {
      showNotification(
        "Nuclei threshold has not been set. Please set a threshold before counting",
        type = 'warning'
      )
      return(NULL)
    }
    if (length(unique(c(nuclei_th()))) == 1) {
      showNotification(
        "Empty or full image after thresholding, modify the threshold and try againg"
      )
      return(NULL)
    }
    n_seg <- watershed(distmap(nuclei_th()))
    # n_objects <- computeFeatures.shape(n_seg)[,'s.area']
    n_objects <- c(table(c(n_seg), exclude = 0)) # Way faster
    apply_th(TRUE)
    # Set slider area
    nmax <- as.integer(max(n_objects))
    updateSliderInput(session, 'nsize', value=c(50,nmax),max = nmax)
    nuclei_seg(n_seg)
    nuclei_objects(n_objects)
  })

  filtoutByArea <- reactive({
    nobjs <- nuclei_objects()
    names(nobjs[nobjs < input$nsize[1] | nobjs > input$nsize[2]])
  })

  # Update nuclei on edges values
  filtoutByEdge <- reactive({
    req(nuclei_seg())
    if (input$rmborder == TRUE) {
      unique(c(
        nuclei_seg()[1, ],
        nuclei_seg()[nrow(nuclei_seg()), ],
        nuclei_seg()[, 1],
        nuclei_seg()[, ncol(nuclei_seg())]
      ))
    } else {
      NULL
    }
  })

  # Reactive to compute nuclei filtered
  nuclei_filtered <- reactive({
    req(nuclei_seg())
    labs2Rm <- c(filtoutByArea(), filtoutByEdge())
    if (!is.null(labs2Rm)) {
      rmObjects(nuclei_seg(), labs2Rm, reenumerate = F)
    } else {
      nuclei_seg()
    }
  })

  # Reactive to compute nuclei filtered
  nuclei_count <- reactive({
    req(nuclei_filtered())
    length(unique(c(nuclei_filtered()))) - 1
  })

  # Reset apply_th when loading a new image
  observeEvent(imgpath(), {
    apply_th(FALSE)
  })

  # Reset apply_th when changing threshold
  observeEvent(input$thres.n, {
    apply_th(FALSE)
  })

  ## SABGAL Processing ----

  # Create sabgal image when selecting a channel
  sabgal <- reactive({
    req(valid_img())
    validate(need(input$schan,"Please select a channel"))
    tryCatch(
      get_channel(valid_img(), input$schan),
      error = function(e) {
        return(NULL)
      }
    )
  })

  # Reactive storing thresholded sabgal
  sabgal_th <- reactive({
    req(sabgal())
    sabgal() >= input$thres.s[1] & sabgal() <= input$thres.s[2]
  })

  # ObserveEvent for otsu
  observeEvent(input$runotsu.s, {
    req(sabgal())
    deft <- round(otsu(sabgal(), c(0, img_bitdepth() - 1), img_bitdepth()), 1)
    updateSliderInput(session, 'thres.s', value = c(deft, img_bitdepth() - 2))
  })

  # ObserveEvent for reset
  observeEvent(input$reset.ths, {
    updateSliderInput(session, 'thres.s', value = c(0, img_bitdepth() - 1))
  })

  # Reactive to get sabgal coordinates
  # hover_coords <- reactive({
  #   # req(input$splot_hover) # Renmove or bgalstats will be invalidated
  #   coords <- unlist(input$splot_hover[c('x', 'y')])
  #   if (!is.null(coords)) {
  #     if (any(coords < 0)) {
  #        NULL
  #     } else {
  #       round(coords)
  #     }
  #   }
  # })

  # Bgal stats
  sabgal_stats <- reactive({
    req(sabgal(), sabgal_th())
    calc_bgalstats(sabgal(), sabgal_th(), crop_region())
  })

  ### Measure and runAll ----

  batch_results <- reactiveVal(empty_res())

  # Reset button
  observeEvent(input$reset, {
    batch_results(empty_res())
  })

  # Observer for "Measure" button
  observeEvent(input$run, {
    req(imgpath(), input$schan)
    checkrun(input,img_bitdepth())
    withProgress(
      {
        # Process image
        imgpath <- file.path(mydir(), imgpath())
        res <- process_single_image(imgpath, input, outdir())
        res <- c(list(imgpath()), res)
        # Update batch_results
        cbres <- batch_results()
        cbres[nrow(cbres) + 1, ] <- res
        batch_results(cbres)
      },
      message = "Calculating current image"
    )
  })

  # Observer for "Run all images" button
  observeEvent(input$run.all, {
    req(mydir())
    req(input$schan)
    # Get filtered files using the reactive value
    if (length(filtered_files()) == 0) {
      showNotification(
        "No files found matching the current filter",
        type = "message"
      )
      return()
    }
    # Process all images
    checkrun(input,img_bitdepth())
    image_files <- file.path(mydir(), filtered_files())
    # Process all images
    batch_results(process_all_images(image_files, input,outdir()))
  })

  # Observer for "Click Here To Learn More"
  observeEvent(input$help, {
    browseURL("https://github.com/antotartier/FAB-gal/wiki")
  })

  ## Outputs ----

  ### Nuclei Plots ----
  
  # Reactive holding the type of nuclei display
  nuclei2disp <- reactive({
    req(nuclei())
    if (apply_th() == TRUE) {
      req(nuclei_filtered())
      return(colorLabels(nuclei_filtered()))}
    req(img_bitdepth())
    if (all(input$thres.n == c(0, img_bitdepth() - 1))) {
      req(nuclei_bgs())
      return(normalize(
        nuclei_bgs(), 
        inputRange = c(0, img_bitdepth() - 1)
        ))
    }
    req(nuclei_th)
    nuclei_th()
  })
  
  # Reactive to store the crop region (coordinates only)
  crop_region <- reactiveVal(NULL)
  
  # Reactive expressions to compute displayed images based on crop region
  nplot_displayed <- reactive({
    req(nuclei2disp())
    img <- nuclei2disp()
    region <- crop_region()
    if (!is.null(region)) {
      # Handle both 2D and 3D images (grayscale vs RGB)
      if (length(dim(img)) == 2) {
        img <- img[region$xmin:region$xmax, region$ymin:region$ymax]
      } else {
        img <- img[region$xmin:region$xmax, region$ymin:region$ymax, ]
      }
    }
    img
  })
  
  # Apply zoom when brush is drawn on nuclei plot
  observeEvent(input$nplot_brush, {
    req(nplot_displayed())
    region <- get_crop_region_from_brush(nplot_displayed(), input$nplot_brush, crop_region())
    if (!is.null(region)) {
      crop_region(region)
    }
  })
  
  # Reset zoom for nuclei
  observeEvent(input$reset.zoom.n, {
    crop_region(NULL)
  })
  
  
  output$nplot <- renderPlot({
    display(nplot_displayed(), method = 'raster')
  })

  ### SABGAL Plot ----
  
  # Reactive holding the sabgal image to display
  sabgal2disp <- reactive({
    req(img_bitdepth())
    if (all(input$thres.s == c(0, img_bitdepth() - 1))) {
      req(sabgal())
      normalize(sabgal(), inputRange = c(0, img_bitdepth() - 1))
    } else {
      req(sabgal_th())
      sabgal_th()
    }
  })
  
  # Reactive holding the cropped sabgal image to display
  splot_displayed <- reactive({
    img <- sabgal2disp()
    region <- crop_region()
    if (!is.null(region)) {
      img <- img[region$xmin:region$xmax, region$ymin:region$ymax]
      }
    img
  })
  
  # Apply zoom when brush is drawn on sabgal plot
  observeEvent(input$splot_brush, {
    req(splot_displayed())
    region <- get_crop_region_from_brush(splot_displayed(), input$splot_brush, crop_region())
    if (!is.null(region)) {
      crop_region(region)
    }
  })
  
  # Reset zoom for sabgal
  observeEvent(input$reset.zoom.s, {
    crop_region(NULL)
  })
  
  output$splot <- renderPlot({
    display(splot_displayed(), method = 'raster')
  })

  ### Nuclei stats  ----
  output$stats.n <- renderPrint({
    validate(need(nuclei(), "Select an image and channel"))
    if (!apply_th()) {
      outstr <- "NA (press count)"
    } else {
      outstr <- as.character(nuclei_count())
    }
    cat('Nuclei count: ', outstr)
  })

  ### SABGAL stats ----
  output$stats.s <- renderPrint({
    req(sabgal_stats())
    sabst <- sabgal_stats()
    cat(sprintf("MFI = %.2f  Sel_Area = %.2f%%",sabst$mfi,sabst$perA))
    })

  ### Results table ----
  output$results_table <- renderTable({
    results <- batch_results()
    if (nrow(results) > 0) {
      results
    }
  })

  ### Download results ----
  output$download_results <- downloadHandler(
    filename = function() {
      paste("FluoroSABGAL_results_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(batch_results(), file, row.names = FALSE)
    }
  )
}


# Run the application  ----
shinyApp(ui = ui, server = server)
