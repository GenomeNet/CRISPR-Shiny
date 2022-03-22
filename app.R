## app.R ##
library(shinydashboard)
library(shiny)
library(keras)
library(deepG)
library(ggplot2)
library(dplyr)
library(circlize)
library(shinycssloaders)
library(DT)
library(hdf5r)
library(ComplexHeatmap) # sudo apt-get install libcairo2-dev
#library(plotly)impo
library(zoo)

batch_size <- 30000

# crispr_gap: After what numbber of negative predictions to start new CRISPR region.
# conf_cutoff: Confidence decision threshold.
# pos_rate: What percentage of predictions must be above conf_cutoff.
# min_seq_len: Minimum sequence length.
filter_crispr <- function(states_df,
                          crispr_gap = 10,
                          conf_cutoff = 0.5,
                          pos_rate = 0.8,
                          min_seq_len = 120,
                          maxlen = 200) {
  
  stopifnot(all(c("conf_CRISPR", "conf_non_CRISPR", "seq_end") %in% colnames(states_df)))
  step <- states_df$seq_end[2] - states_df$seq_end[1]
  crispr_list <- list()
  states_df <- states_df %>% dplyr::filter(conf_CRISPR > conf_cutoff)
  states_df <- states_df[order(states_df$seq_end), ]
  row_num <- nrow(states_df)
  crispr_index <- 1
  
  crispr_start <- states_df$seq_end[1]
  for (i in 1:(row_num-1)) {
    
    l_name <- paste0("CRISPR_region_", crispr_index)
    current_pos <- states_df$seq_end[i]
    next_pos <- states_df$seq_end[i+1]
    
    if ((abs(current_pos - next_pos) > crispr_gap) & (i != (row_num-1))) {
      
      index <- (states_df$seq_end >= crispr_start) & (states_df$seq_end <= current_pos)
      crispr_list[[l_name]] <- states_df[index, ]
      crispr_start <- next_pos
      crispr_index <- crispr_index + 1
    }
    
    if (i == (row_num-1)) {
      if (abs(current_pos - next_pos) <= crispr_gap) {
        index <- states_df$seq_end >= crispr_start
        crispr_list[[l_name]] <- states_df[index, ]
      } else {
        index <- (states_df$seq_end >= crispr_start) & (states_df$seq_end <= current_pos)
        crispr_list[[l_name]] <- states_df[index, ]
        
        # single sample at end 
        crispr_index <- crispr_index + 1
        l_name <- paste0("CRISPR_region_", crispr_index)
        crispr_list[[l_name]] <- states_df[nrow(states_df), ]
      }
    }
  } 
  
  # filter by positivity rate
  for (i in names(crispr_list)) {
    df <- crispr_list[[i]]
    seq_len <- df$seq_end[nrow(df)] - df$seq_end[1] + 1
    ## consider step size
    num_possible_pos_pred <- ((seq_len - 1)/step) + 1
    cov_rate <- nrow(df)/num_possible_pos_pred 
    if (cov_rate < pos_rate) {
      crispr_list[[i]] <- NULL
    } 
  }
  
  # filter by size 
  for (i in names(crispr_list)) {
    df <- crispr_list[[i]]
    seq_len <- df$seq_end[nrow(df)] - df$seq_end[1] + 1
    if (seq_len < min_seq_len) {
      crispr_list[[i]] <- NULL
    } 
  }
  
  for (i in names(crispr_list)) {
    df <- crispr_list[[i]]
    df$seq_middle <- df$seq_end - (maxlen/2)
    crispr_list[[i]] <- df
  }
  crispr_list
}

model <- keras::load_model_hdf5("models/crispr_model_el4.hdf5", compile = FALSE)

num_layers <- length(model$get_config()$layers)
layer_name <- model$get_config()$layers[[num_layers]]$name
layer2_name <- model$get_config()$layers[[num_layers-1]]$name

dummy <- readChar("dummy.txt",
                  file.info("dummy.txt")$size)

ui <- dashboardPage(#skin = "black",
  dashboardHeader(title = "CRISPR prediction", titleWidth = 350, disable = T),
  
  dashboardSidebar(
width = 350, disable = F,
                   box(width = 12, solidHeader = TRUE,
                       p(
                         class = "text-muted",
                         paste("Deep-CRISPR finder. This service will predict the CRISPR propability of a input sequence.")
                         )),
                   box(width = 12, solidHeader = TRUE,
                       uiOutput('resetable_input'),
                       
                       
                       actionButton("run", "Predict"),
                       actionButton("reset_input", "Reset inputs")
                   ),
                   
                   box(width = 12, solidHeader = TRUE,
                       p(class = "text-muted",
                         paste("Funded by Priority Programme 2141. Much more than defence: the multiple functions and facets of CRISPR-Cas. Funded in part by BMBF Project GenomeNet (031L0199A/031L0199B)")
                       ))
                   
  ),
  dashboardBody(
    tags$head(tags$style(HTML('
    
                                #text {height: 200px; word-wrap: break-word; font-family: monospace; overflow-y:scroll;}

                                /* logo */
                                .skin-blue .main-header .logo {
                                background-color: #ffffff;
                                }
                                
                                /* logo when hovered */
                                .skin-blue .main-header .logo:hover {
                                background-color: #ffffff;
                                }
                                
                                /* navbar (rest of the header) */
                                .skin-blue .main-header .navbar {
                                background-color: #ffffff;
                                }
                              
                                .box-title {
                                  color:#000000;
                                  font-family: "Source Sans Pro","Helvetica Neue","Helvetica","Arial","sans-serif";
                                  background:#e8e8e8
                                }
                            
                                .box.box-solid.box-primary>.box-header {
                                  color:#000000;
                                  background:#e8e8e8
                                                    }
                                
                                .box.box-solid.box-primary{
                                border-bottom-color:#e8e8e8;
                                border-left-color:#e8e8e8;
                                border-right-color:#e8e8e8;
                                border-top-color:#e8e8e8;
                                }
                                                                
                                /* main sidebar */
                                .skin-blue .main-sidebar {
                                background-color: #e8e8e8;
                                }
                                
                                /* control lable */
                                .skin-blue .control-label {
                                background-color: #55555;
                                color:#000000;
                                }
                                
                                  /* control lable */
                                .checkbox {
                                background-color: #55555;
                                color:#000000;
                                }
                                
                                /* active selected tab in the sidebarmenu */
                                .skin-blue .main-sidebar .sidebar .sidebar-menu .active a{
                                background-color: #ff0000;
                                }
                                
                                /* other links in the sidebarmenu */
                                .skin-blue .main-sidebar .sidebar .sidebar-menu a{
                                background-color: #00ff00;
                                color: #000000;
                                }
                                
                                /* other links in the sidebarmenu when hovered */
                                .skin-blue .main-sidebar .sidebar .sidebar-menu a:hover{
                                background-color: #ff69b4;
                                }
                                /* toggle button when hovered  */
                                .skin-blue .main-header .navbar .sidebar-toggle:hover{
                                background-color: #ff69b4;
                                }

                                /* body */
                                .content-wrapper, .right-side {
                                background-color: #ffffff;
                                }
                                
                                '))),
    
  
    
    fluidRow(
      box(
        id = "box2", title = "Model Confidence Plot", width = 12, solidHeader = TRUE, status = "primary",
        plotOutput("plot1", height = 250) %>% withSpinner(color="#000000", type = 8),
      ),
    ),
      fluidRow(
        box(title = "Matches", width = 8, solidHeader = TRUE, status = "primary",
            dataTableOutput('table'), 
            uiOutput("downloadDataUI"),
            

        ),
        box(title = "Post-processing", width = 4, solidHeader = TRUE, status = "primary",
                 numericInput("threshold_conf", "Min average confidence:", 0.8, step = .1, min = 0, max = 1),
                 numericInput("threshold_pos_rate", "Positive rate:", 0.8, step = .1, min = 0, max = 0.99),
                 numericInput("threshold_min_seq_len", "Min array length:", 100, step = 10,  min = 1, max = 100),
                 numericInput("threshold_gap", "Min gap between loci:", 50, step = 10, min = 1, max = 1000),
                 numericInput("adjacent", "Number of nt to plot adjacent to loci:", 50, step = 10, min = 1, max = 1000))
    
      ),
      
      fluidRow(
      
      box(
        title = "Detailed View", width = 12, solidHeader = TRUE, status = "primary",
        plotOutput("plot3", height = 250)) %>% withSpinner(color="#000000", type = 8),
      ),
    fluidRow(
      
      box(
        title = "States View", width = 12, solidHeader = TRUE, status = "primary",
        plotOutput("plot4", height = 450) %>% withSpinner(color="#000000", type = 8),
     
      )#,
      #  box(
      #    title = "Importance scores", width = 12, solidHeader = TRUE, status = "primary",
      #    plotOutput("plot5", height = 250) %>% withSpinner(color="#000000", type = 8),
      #  )
    ),
    fluidRow(
      box( title = "Sequence", width = 12, solidHeader = TRUE, status = "primary",
           textOutput("text")
           )
      
   
      # box(width = 12, status = "warning",
      #    p(class = "text-muted",
      #       paste("No predictions have been computed."))),
      
    )
   
  )
)

server <- function(input, output) {
  
  
  # download button
  output$downloadDataUI <- renderUI({
    req(dataSummary())
    downloadButton('downloadData', label = 'Download Table')

    })
  

  
  output$resetable_input <- renderUI({
    times <- input$reset_input
    div(id=letters[(times %% length(letters)) + 1],
        
        helpText("Specify the input sequence (nt only)"),
        
        textAreaInput("text", "Sequence input", 
                      height = 150,
                      width = 350,
                      value = dummy, 
                      placeholder = NULL),
        helpText("Upload a file in the .FASTA format (only first sequence will be considered)"),
        
        fileInput("fasta_path", width = 350, label = "or Input a FASTA file"),
        helpText("Predictions will be performed within a sliding window over the whole sequence."),
        
    )
  })
  
  
  #https://f000.backblazeb2.com/file/bioinf/crispr/logo.png
  output$png <- renderUI({
    tags$a(img(src = "https://f000.backblazeb2.com/file/bioinf/crispr/logo.png", width = "50px"), href = "http://deepg.de", target = "_blank")
  })
  
  fastaData <- reactive({
    req(input$fasta_path$datapath)
    fasta_file <- microseq::readFasta(input$fasta_path$datapath) #tibble
    fasta_file$Sequence[1]
  })
  
  # the nt string
  inputSequenceFull <- reactive({
    if(is.null(input$fasta_path$datapath)) {
      as.character(input$text) # use textarea data
    } else {
      fastaData() # use uploaded file
    }
  })
    
  # the nt string subsetted after clicking on table
  inputSequenceSubset <- reactive({
   req(inputSequenceFull())
   req(dataSummary())
   req(input$table_row_last_clicked)
   subset <- substr(inputSequenceFull(), 
                    dataSummary()[input$table_row_last_clicked,]$Start,
                    dataSummary()[input$table_row_last_clicked,]$End)
   subset
  })
  
  
  dataSummary <- reactive({
    req(dataInput())
    states_df = dataInput()
    crispr_list <- filter_crispr(states_df = dataInput(),
                                 crispr_gap = input$threshold_gap,
                                 conf_cutoff = input$threshold_conf,
                                 pos_rate = input$threshold_pos_rate,
                                 min_seq_len = input$threshold_min_seq_len)
    
    step <- states_df$seq_end[2] - states_df$seq_end[1]
    placeholder <- rep(0, length(crispr_list))
    
    crispr_summary_df <- data.frame(seq_len = placeholder, cov_rate = placeholder, 
                                    average_confidence = placeholder,
                                    first_pred = placeholder, last_pred = placeholder)
    count <- 1
    for (i in names(crispr_list)) {
      df <- crispr_list[[i]]
      seq_len <- df$seq_end[nrow(df)] - df$seq_end[1] + 1
      crispr_summary_df$seq_len[count] <- seq_len
      num_possible_pos_pred <- ((seq_len - 1)/step) + 1
      crispr_summary_df$cov_rate[count] <- round(nrow(df)/num_possible_pos_pred, digits = 3)
      crispr_summary_df$first_pred[count] <- min(df$seq_middle)
      crispr_summary_df$last_pred[count] <- max(df$seq_middle)
      # get average crispr confidence
      subset_index <- states_df$seq_end >= min(df$seq_end) & states_df$seq_end <= max(df$seq_end)
      states_df_subset <- states_df[subset_index, ]
      crispr_summary_df$average_confidence[count] <- round(mean(states_df_subset$conf_CRISPR), digits = 3)
      count <- count + 1
    }
    
    colnames(crispr_summary_df) <- c("Length", "Coverage", "Confidence", "Start", "End")
    crispr_summary_df <- crispr_summary_df[order(-crispr_summary_df$Length),]
    crispr_summary_df
  })
  

  
  
  dataInput <- eventReactive(input$run, {
    withProgress(value = 0, message = 'Running predictions',
                 detail = 'Please wait' ,{
                   incProgress(amount = .25, message = "Loading sequence") 
                   
                   if(is.null(input$fasta_path$datapath)) {
                     input_sequence <- as.character(input$text) # use textarea data
                   } else {
                     input_sequence <- fastaData()
                   }
                   
                   if (file.exists("tmp_lvl1.h5")) file.remove("tmp_lvl1.h5")
                   incProgress(amount = .25, message = "Inference") 
                   deepG::writeStates(model = model,
                                      layer_name = layer_name,
                                      sequence = input_sequence,
                                      round_digits = 4,
                                      filename = "tmp_lvl1.h5",
                                      batch.size = batch_size,
                                      step = 1)
                   incProgress(amount = .5, message = "Loading results") 
                   states <- readRowsFromH5(h5_path = "tmp_lvl1.h5", complete = TRUE, getTargetPositions = TRUE)
                   incProgress(amount = 1, message = "Processing results") 
                   pred <- states[[1]]
                   position <- states[[2]] - 1
                   df <- cbind(pred, position) %>% as.data.frame()
                   colnames(df) <- c("conf_CRISPR", "conf_non_CRISPR", "seq_end")
                   df$seq_middle <- df$seq_end - 100
                   df
                 })
  })
  
  
  
#  eventReactive(input$reset_input, {
#    if(is.null(input$fasta_path$datapath)) {
#      input_sequence <- as.character(input$text) # use textarea data
#    } else {
#      input_sequence <- fastaData()
#    }  })
  
  dataInput2 <- eventReactive(input$run, {
    withProgress(value = 0, message = 'Layer2',
                 detail = 'Please wait' ,{
                   incProgress(amount = .25, message = "Loading sequence") 
                   
                   if(is.null(input$fasta_path$datapath)) {
                     input_sequence <- as.character(input$text) # use textarea data
                   } else {
                     input_sequence <- fastaData()
                   }
                   
                   if (file.exists("tmp_lvl2.h5")) file.remove("tmp_lvl2.h5")
                   incProgress(amount = .25, message = "Inference intermediate layer") 
                   deepG::writeStates(model = model,
                                      layer_name = layer2_name,
                                      sequence = input_sequence,
                                      round_digits = 4,
                                      filename = "tmp_lvl2.h5",
                                      batch.size = batch_size,
                                      step = 1)
                   incProgress(amount = .5, message = "Loading results") 
                   states <- readRowsFromH5(h5_path = "tmp_lvl2.h5", complete = TRUE, getTargetPositions = TRUE)
                   pred <- states[[1]]
                   position <- states[[2]] - 1
                   
                   pred <- data.matrix(pred)
                   rownames(pred) <- paste0("pos ", position)
                   colnames(pred) <- paste0("neuron ", 1:ncol(pred))
                   pred
                 })
  })
  
  output$plot1 <- renderPlot({
    req(dataInput2())
    req(dataInput())
    req(input$threshold_min_seq_len)
    p <- ggplot(data = dataInput(), aes(seq_middle, conf_CRISPR)) #+ geom_point(size = .1, alpha = .2) 
    p <- p + xlab(paste0("Position in sequence")) + ylab("CRISPR confidence")
    #p <- p + geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") 
    p <- p + geom_hline(yintercept = input$threshold_conf, linetype = "dashed", color = "darkgreen") 
    p <- p + theme_classic()
    p <- p + geom_line(aes(y=rollmean(conf_CRISPR, input$threshold_min_seq_len,
                                      na.pad = TRUE)), color = "black", alpha = 1, size = 1) 
    
  #  if(!is.null(input$table_row_last_clicked)){
  #    start <- dataSummary()[input$table_row_last_clicked,]$Start - 10
  #    end <- dataSummary()[input$table_row_last_clicked,]$End + 10
  #    p <- p + geom_vline(xintercept = start, color = "red", linetype = "dashed")
  #    p <- p + geom_vline(xintercept = end, color = "red", linetype = "dashed")
  #  }
    p
  })
  
  output$plot2 <- renderPlot({
    p <- ggplot(data = dataInput(), aes(conf_CRISPR)) 
    p <- p + geom_histogram()
    #p <- p + geom_vline(xintercept = 0.5, color = "black", linetype = "dashed")
    p <- p + theme_classic()
    p
  })
  
  output$table <- renderDataTable(dataSummary(), rownames = FALSE, selection = 'single')
  
  output$plot3 <- renderPlot({
    req(input$table_row_last_clicked)
    req(dataSummary())
    req(dataInput())
    start <- dataSummary()[input$table_row_last_clicked,]$Start
    end <- dataSummary()[input$table_row_last_clicked,]$End 
    p <- ggplot(data = dataInput(), aes(seq_middle, conf_CRISPR))
    p <- p + xlab("Sequence position") + ylab("CRISPR confidence")
    p <- p + geom_vline(xintercept = start, color = "red", linetype = "dashed")
    p <- p + geom_vline(xintercept = end, color = "red", linetype = "dashed")
    p <- p + geom_point(size = .1) 
    p <- p + xlim(start - (input$adjacent), end + (input$adjacent))
    p <- p + geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") 
    p <- p + theme_classic() 
    p
  })
  
  
  output$plot4 <- renderPlot({
    req(inputSequenceFull())
    req(input$table_row_last_clicked)
    req(dataInput())
    req(dataInput2())
    req(dataSummary())

    # get start and end based on fasta length
    start_pos <- ifelse(dataSummary()[input$table_row_last_clicked,]$Start - (input$adjacent) < 0,
                        0,
                        dataSummary()[input$table_row_last_clicked,]$Start - (input$adjacent))
    
    
    end_pos <- ifelse(dataSummary()[input$table_row_last_clicked,]$End - (input$adjacent) > nchar(inputSequenceFull()),
                        nchar(inputSequenceFull()),
                      dataSummary()[input$table_row_last_clicked,]$End  + (input$adjacent))
    

    subset <- dataInput2()[start_pos:end_pos,]
    
    max_value <- max(subset)
    
    loci_length <- dataSummary()[input$table_row_last_clicked,]$End - dataSummary()[input$table_row_last_clicked,]$Start
    
    col_fun = colorRamp2(c(0, max_value/2, max_value), c("white", "orange", "black"))
    ht <- Heatmap(data.matrix(t(subset)),
                  name = "activation",
            #      column_split = c(rep("pre", input$adjacent),
             #                                     rep("main", loci_length),
              #                                    rep("post", input$adjacent)), 
                  col = col_fun,
                  show_row_names = F,
                  border = TRUE,
                  column_title = "Sequence position in loci",
                  row_title = "Neuron", 
                  show_heatmap_legend = FALSE,
                  show_column_names = F,
                  use_raster = F,
                  cluster_rows = F,
                  row_names_gp = gpar(fontsize = 2),
                  column_names_gp = gpar(fontsize = 2),
                  cluster_columns = F)
 
    
    
    ht
  })
  
    
  output$plot5 <- renderPlot({
    req(input$table_row_last_clicked)
    req(dataInput())
    
    if(is.null(input$fasta_path$datapath)) {
      input_sequence <- as.character(input$text) # use textarea data
    } else {
      input_sequence <- fastaData()
    }
    
    maxlen <- 200
    vocabulary = c("A", "C", "G", "T")
    amb_nuc_token <- "B"
    tokenizer_pred <- keras::fit_text_tokenizer(keras::text_tokenizer(char_level = TRUE, lower = TRUE, oov_token = "0"), c(vocabulary, amb_nuc_token))
    nucSeq <- paste(input_sequence, collapse = "")
    nucSeq <- keras::texts_to_sequences(tokenizer_pred, nucSeq)[[1]] - 1
    input <- sequenceToArrayLabel(sequence = nucSeq, maxlen = maxlen, vocabulary = vocabulary, startInd = 1,
                                  ambiguous_nuc = "equal", nuc_dist = NULL, use_quality = FALSE, quality_vector = NULL)
    
    ig <- integrated_gradients(m_steps = 300,
                               baseline_type = "zero",
                               input_seq = input,
                               target_class_idx = 1,
                               model = model,
                               pred_stepwise = FALSE)
    
    ig_array <- as.array(ig)
    hm <- heatmaps_integrated_grad(integrated_grads = ig_array, input_seq = input)[[1]]
    hm
    })
  
  
  output$downloadData <- downloadHandler(
    filename = function() {"results.csv"},
    content = function(file) {
      write.csv(dataSummary(), file, quote = F, row.names = F, sep = ";")
    })
  
  
  # sequence box
  output$text <- renderPrint({ 
    cat(paste0(inputSequenceSubset(), collapse = " "))
  })
    
}

shinyApp(ui, server)