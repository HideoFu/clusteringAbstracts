#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(readr)
library(stringr)
library(dplyr)
library(XML)

library(tm)
library(slam)

library(dendextend)


n_term <- 10  # Number of words to be shown


# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Clustering Abstracts"),
   

   sidebarLayout(
      sidebarPanel(
        fileInput("file", "Select File:"),
        fileInput("sw", "Stop Words File:"),  # Stop words
        
        # Clustering
        wellPanel(
          sliderInput("n_cl", "Cluster into:", value = 1, min = 1, max = 10),  # Cluster No.
          actionButton("draw", "Cluster"),
          downloadButton("dl_sum", "Top_words")
        ),
        
        # Labeling
        wellPanel(
          textAreaInput("kw", "Keywords:", 
                        placeholder = "Separate words with ','. Regex can be used."),  # Key word
          actionButton("label", "Label"),
          downloadButton("dl_lab", "KW count")
        ),
        
        
        conditionalPanel(condition = "input.n_cl > 1",
                         uiOutput("ui_sel"),
                         downloadButton("dl", "Subset")
                         )
        
      ),
      
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("dendrogram"),
         tableOutput("freq_wd"),
         tableOutput("kw_sum")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  # Make UI
  output$ui_sel <- renderUI({
    text <- 1
    
    if(input$n_cl > 1){
      for(i in 2:input$n_cl){
        text <- c(text, i)
      }
    }
    
    selectInput("cl", "Clusters to DL:", choices = text, multiple = T)
  })
  
  
  # Read Data
  raw_data <- reactive({
     ext <- str_sub(input$file$name, -3,-1)
     
     if (ext == "csv"){
       data <- read_csv(input$file$datapath)
       if("clust" %in% names(data)){
         data <- data[,-which(names(data) == "clust")]
       }
     } else if(ext == "xml"){
       xml <-xmlTreeParse(input$file$datapath, useInternalNodes = T)
       
       id <- xpathSApply(xml, "//PMID", xmlValue)
       title <- xpathSApply(xml, "//ArticleTitle", xmlValue)
       abst <- xpathSApply(xml, "//Abstract", xmlValue)
       
       data <- data_frame(doc_id = id, text = abst, title = title)
     }
     
     return(data)
  })
   
  # Make Clean Corpus
   clean_corp <- reactive({
     ## clean special chars
     data <- raw_data()
     
     # data$Abstract <- gsub("\\r*\\n*", "", data$Abstract)
     # data$PubDate <- data$PubDate %>% 
     #   str_replace_all("\\r\\n", "/") %>%
     #   str_replace_all("^/|/$", "")
     # data$AuthorList <- gsub("\\r\\n", ", ", data$AuthorList)
     # data$KeywordList <- data$KeywordList %>%
     #   str_replace_all("\\r\\n", ", ") %>%
     #   str_replace("^, ", "")
     # 
     # 
     # ## rename / sort columns
     # data <- data %>%
     #   mutate(doc_id = seq(1, nrow(data))) %>%
     #   select(doc_id, Abstract)
     
     colnames(data)[1:2] <- c("doc_id","text")
     
     
     ## convert to corpus
     source <- DataframeSource(data)
     corp <- VCorpus(source)
     
     
     
     ### Stemming function
     stemming <- function(text){
       char_vec <- unlist(strsplit(text, split = " "))
       stem_doc <- stemDocument(char_vec)
       new_text <- paste(stem_doc, collapse = " ")
       return(new_text)
     }
     
     ### Function for remove "[number]%", "[num]/[num]"
     rm_res_num <- function(text){
       # new_text <- gsub("([[:digit:]]|\\.)+%", "", text)
       # new_text <- gsub("([[:digit:]]|/)", "", new_text)
       new_text <- text %>% str_remove_all("([[:digit:]]|\\.)+%") %>% 
         str_remove_all("[[:digit:]]/[[:digit:]]")
       
       return(new_text)
     }
     
     
     clean_corp <- tm_map(corp, content_transformer(tolower)) %>%
       tm_map(content_transformer(rm_res_num)) %>%
       tm_map(removePunctuation, preserve_intra_word_dashes = T) %>%
       tm_map(stripWhitespace) %>%
       tm_map(content_transformer(stemming))  # must place after remove punctuations
     
     ## clean corpus
     
     if(!is.null(input$sw)){
       sw <- scan(input$sw$datapath, what = character())
       sw <- unique(stemDocument(sw))
       clean_corp <- tm_map(clean_corp, removeWords, sw)
       
       # sw <- paste0(" ", sw, " ")
       # sw <- paste(sw, collapse = "|")
       # rm_sw <- function(text){
       #   #new_text <- gsub(sw, "", text)
       #   new_text <- str_replace_all(text, sw, " ")
       #   return(new_text)
       # }
       # clean_corp <- tm_map(clean_corp, content_transformer(rm_sw))
     }
     
     return(clean_corp)
     })
   
   # Make TDM
   tdm <- reactive({
     TermDocumentMatrix(clean_corp(), control = list(weighting = weightBin, wordLengths = c(2, Inf)))
   })
   
   
   # Make hc cluster
   hc <- reactive({
     # Calculate Cosine Distance
     cos_dist_mat <- 1 - crossprod_simple_triplet_matrix(tdm())/(sqrt(col_sums(tdm()^2) %*% t(col_sums(tdm()^2))))
     
     # Make Cluster
     hc <- hclust(as.dist(cos_dist_mat), method = "ward.D")
   })
   
   
   # Make dendrogram obj
   dend <- reactive({
     as.dendrogram(hc())
   })
   
   
   # Cut into n_cl cluster
   clust <- reactive({
     input$draw
     
     isolate({
       if(input$n_cl > 1){
         clust <- cutree(dend(), input$n_cl)
       }
     })
   })
   
   # Make key words into vector
   kw_vec <- reactive({
     input$label
     
     isolate({
       kws <- unlist(str_split(input$kw, pattern = "[,.;][[:blank:]]*"))
       kws <- kws %>%
         tolower() %>%
         stemDocument()
       
       kws <- paste0("^", kws, "$")
       
       terms <- tdm()$dimnames$Terms
       
       term <- terms[grepl(paste(kws, collapse = "|"), terms)]
       
       return(term)
     })
   })
   
   # Draw Dendrogram (w/ Cluster)
   ## select doc_id which has KWs
   docs <- reactive({
     input$label
     
     isolate({
       if(input$kw != ""){
         kws <- paste0("^", kw_vec(), "$") %>%
           paste(collapse = "|")
         
         tdm <- tdm()
         
         row_no <- grep(kws, tdm$dimnames$Terms)
         docs <- tdm$dimnames$Docs[tdm$j[tdm$i %in% row_no]]
       }
     })
   })
   
   ## Draw
   output$dendrogram <- renderPlot({
     input$draw
     input$label
     
     isolate({
       if(is.null(input$file)){
         return(NULL)
       }
       
       dend <- dend()
       
       ## Emphasize document including key word
       if(input$kw != ""){
         dend <- branches_attr_by_labels(dend, docs(), "red")
       }
       
       plot(dend)
       
       # Cluster on dendrogram
       if(input$n_cl > 1){
         # Order clusters as in dendrogram
         cl <- unique(clust()[order.dendrogram(dend())])
         
         rect.dendrogram(dend(), k = input$n_cl, border = "red", 
                         text = cl, text_cex = 2, text_col = "red")
       }
     })
   })
   
   
   # Make DTM
   dtm <- reactive({
     DocumentTermMatrix(clean_corp(), control = list(weighting = weightBin, wordLengths = c(2, Inf)))
   })
   
   # Make freq word table
   freq_tbl <- reactive({
     input$draw
     
     isolate({
       if(input$n_cl > 1){
         
         dtm_mat <- as.matrix(dtm())
         
         p_words <- colSums(dtm_mat)/sum(dtm_mat)
         
         # calulate word possibility deviations of each cluster over whole documents
         cl_words <- lapply(unique(clust()), function(x){
           # filter rows
           rows <- as.matrix(dtm_mat[clust() == x, ])
           # select cols
           rows <- as.matrix(rows[, colSums(rows) >0])
           # deviation
           colSums(rows) / sum(rows) - p_words[colnames(rows)]
         })
         
         
         # Summary Table
         data.frame(cluster = unique(clust()),
                     size = as.integer(table(clust())),
                     top_words = sapply(cl_words, function(d){
                       paste(
                         names(d)[order(d, decreasing = TRUE)][1:n_term],
                         collapse = ", "
                       )
                     }),
                     stringsAsFactors = F)
       }
     })
   })
   
   output$freq_wd <- renderTable({
     freq_tbl()
   })
   
   #Key words count table
   ## Count documents including kw for each cluster
   kw_sum_df <- reactive({
     input$label
     
     isolate({
       # dtm: Document Temr Matrix in Simple triplet matrix
       # kw_vec():  key words vector
       # clust(): cutree obj
       
       
       dtm_df <- as.matrix(dtm()) %>% as_data_frame()
       
       # make table: keyword count per cluster
       kw_sum <- dtm_df %>%
         select(kw_vec()) %>%
         bind_cols(cluster = clust()) %>%
         group_by(cluster) %>%
         summarise_all(funs(sum)) %>%
         left_join(freq_tbl(), by = "cluster") %>%
         select(cluster, size, kw_vec())
       
       
       # mutate counts to percentage
       size <- kw_sum$size
       total <- sum(size)
       
       in_each <- round(kw_sum[,kw_vec()] / size * 100, 1)
       in_total <- round(kw_sum[,"size"] / total * 100, 0)
       
       ## all KWs percentage
       kws_in_each <- vector()
       
       for (i in 1:length(kw_sum$cluster)){
         clust <- clust()
         docs <- docs()
         
         id_cl <- clust[clust == i]  # named vector of doc_id in cluster i
         kws_counts <- sum(unique(docs) %in% names(id_cl))
         
         kws_in_each[i] <- round(kws_counts / size[i] * 100, 0)
       }
       
       
       ## percentage into character
       in_perc <- function(df){
         as_data_frame(apply(apply(df, 2, as.character), 2, paste0, "%"))
       }
       
       in_each_perc <- in_perc(in_each)
       in_total_perc <- in_perc(in_total)
       
       kws_perc <- in_perc(data.frame(KWs = kws_in_each))
       
       
       ## make into data_frame
       kw_sum2 <- kw_sum %>%
         select(cluster) %>%
         bind_cols(in_total_perc, kws_perc, in_each_perc)
         
       
       return(kw_sum2)
     })
   })
   
   ## render Table
   output$kw_sum <- renderTable({
     input$label
     
     isolate({
       if(input$kw == ""){
         return(NULL)
       }
       
       return(kw_sum_df())
     })
   })
   
   
   # DL
   ## subset
   output$dl <- downloadHandler(
                                filename = function(){
                                  root <- str_remove(input$file$name, "\\.csv$")
                                  cls <- paste(as.vector(as.integer(input$cl)), collapse = "_")
                                  paste0(root, "-", cls, ".csv")
                                  },
                                content = function(file){
                                  
                                  dl_data <- raw_data() %>%
                                    bind_cols(data_frame(clust = clust())) %>%
                                    filter(clust %in% as.vector(as.integer(input$cl)))
                                  
                                  write.csv(dl_data, file, row.names = F)
                                }
                                )
   
   # top words list
   output$dl_sum <- downloadHandler(
                                    filename = function(){
                                      root <- str_remove(input$file$name, "\\.csv$")
                                      paste0(root, "_top_words.csv")
                                    },
                                    content = function(file){
                                      write.csv(freq_tbl(), file, row.names = F)
                                    }
                                    )

   # key words counts
   output$dl_lab <- downloadHandler(
                                    filename = function(){
                                      root <- str_remove(input$file$name, "\\.csv$")
                                      paste0(root, "_kw_counts.csv")
                                      },
                                    content = function(file){
                                      write.csv(kw_sum_df(), file, row.names = F)
                                    }
                                    )
     
   
}

# Run the application 
shinyApp(ui = ui, server = server)

