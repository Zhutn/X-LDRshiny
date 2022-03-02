# Load packages ----
conf=read.table("X-LD.conf", as.is = T)

# install.packages(c("shiny","bsplus","RColorBrewer","corrplot"))
library(shiny)
library(bsplus)
library(RColorBrewer)
library(corrplot)

unzip = "unzip -o -q plink.zip"
 
if(length(grep("linux",sessionInfo()$platform, ignore.case = TRUE))>0) {
  print("linux")
  system(paste0(unzip))
  system("chmod a+x ./plink/plink_linux")
  plink2 = "./plink/plink_linux"
} else if(length(grep("apple",sessionInfo()$platform, ignore.case = TRUE))>0) {
  print("apple")
  system(paste0(unzip))
  system("chmod a+x ./plink/plink_mac")
  plink2 = "./plink/plink_mac"
  #system("git rev-list head --max-count 1 > gitTag.txt")
} else {
  print("windows")
  system(paste0(unzip))
  system("chmod a+x ./plink/plink_win.exe")
  plink2 = "./plink/plink_win.exe"
}


options(shiny.maxRequestSize=conf[1,2]*1024^2, shiny.launch.browser=T)
gTag=read.table("gitTag.txt")

# Define UI for X-LD Application
ui <- fluidPage(
        theme = "style.css",
        div(style = "padding: 1px 0px; width: '100%'",
        titlePanel(
          title = "",
          windowTitle = "X-LD"
        )
      ),
      navbarPage(
        title = div(
          span(
            HTML("<input type=button style='font-size:30px;border:0;height:35px' value='X-LD' onclick=\"window.history.go(-1)\">"),
            style = "position: relative; top: 30%; transform: translateY(-50%);"
          )
        ),
        id = "inNavbar",
        tabPanel(
          title = "Data Input",
          value = "datainput",
          fluidRow(
            column(
              4,
              fileInput(
                "file_input",
                paste0('Source files (.bim, .bed, .fam) [< ', conf[1,2],' MB]'),
                multiple = TRUE,
                accept = c("bed", "fam", "bim")
              ) %>%
              shinyInput_label_embed(
                icon("question-circle") %>%
                bs_embed_popover(
                  title = "The chromosome index (the first column of the .bim file) must be numeric", content = "", placement = "right"
                )
              )
            ),
            column(
              4,
              numericInput('autosome', "Autosome number", value=22)
            )
          ),
          hr(),
          fluidRow(
            column(
              4,
              radioButtons(
                'bred',
                'Population type',
                choices = list('Outbred' = 'outbred', 'Inbred' = 'inbred'),
                selected = 'outbred'
                #inline = T
              ) %>%
              shinyInput_label_embed(
                icon("question-circle") %>%
                bs_embed_popover(
                  title = "INBRED is chosen if your sample has homogenous genome, otherwise choose OUTBRED", content = "", placement = "right"
                )
              )
            ),
            column(
              4,
              sliderInput(
                'maf_cut',
                'MAF threshold',
                 value = 0.05, min = 0.01, max = 0.1, step = 0.01
              ) %>%
              shinyInput_label_embed(
                icon("question-circle") %>%
                bs_embed_popover(
                  title = "Marker with allele frequency lower than the given MAF threshold will be filtered out", content = "", placement = "right"
                )
              )
            ),
            column(
              4,
              sliderInput(
                'geno_cut',
                'Missing genotype rates',
                value = 0.2, min = 0, max = 1, step = 0.1
              ) %>%
              shinyInput_label_embed(
                icon("question-circle") %>%
                bs_embed_popover(
                  title = "Marker with missing call rates exceeding the provided value will be filtered out", content = "", placement = "right"
                )
              )
            )
          ),
          hr(),
          fluidRow(
            column(
              4,
              numericInput('ld_windows', "LD pruning (windows size) [optional]", value="") %>%
              shinyInput_label_embed(
                icon("question-circle") %>%
                bs_embed_popover(
                  title = "Window size in variant count or kilobase", content = "", placement = "right"
                )
              )
            ),
            column(
              4,
              numericInput('step_size', "LD pruning (step size) [optional]", value="") %>%
              shinyInput_label_embed(
                icon("question-circle") %>%
                bs_embed_popover(
                  title = "Variant count to shift the window at the end of each step", content = "", placement = "right"
                )
              )
            ),
            column(
              4,
              sliderInput(
                'r_square',
                'LD pruning (R square) [optional]',
                value = '0', min = 0, max = 1, step = 0.1
              ) %>%
              shinyInput_label_embed(
                icon("question-circle") %>%
                bs_embed_popover(
                  title = "Pairs of variants in the current window with squared correlation greater than the threshold will be filtered out", content = "", placement = "right"
                )
              )
            )
          ),      
          fluidRow(
            column(
              4,
              actionButton(
              'run',
              'Run (X-LD)!'
              )
            )
          )
        ),
        tabPanel(
          title = "Visualization",
          value = "visualization",
          fluidRow(
            column(12, 
              mainPanel(
                tabsetPanel(
                  id = 'X-LDFunctions',
                  type = 'tabs',
                  tabPanel(
                    'chromosome level LD',
                    plotOutput('me'),
                    htmlOutput('me_note')
                  ),
                  tabPanel(
                    'chromosome level LD (scaled)',
                    plotOutput('me_scale'),
                    htmlOutput('me_scale_note')
                  )
                 )
              )
            )
          ),
          hr(),
          fluidRow(
            column(
              12,
              mainPanel(
                column(3,
                  downloadButton(
                    'figure', 
                    'Figure'
                  )
                ),
                column(3,
                  downloadButton(
                    'LD', 
                    'LD matrix'
                  )
                )
              )
            )
          )
        ),
        tabPanel(
          title = "About",
          value = "about",
          tags$h3("Source Code"),
          tags$p(HTML("For X-LD core algorithm implementation and R Shiny code in this web tool please refer to")),
          tags$p(HTML("<a href=\"https://github.com/huangxin0221/X-LDRshiny\" target=\"_blank\">GitHub repository: X-LDRshiny.</a>")),
          tags$br(),
          tags$h3("Citation"),
          tags$p(HTML("<a>Xin Huang et al, X-LD: a fast and effective algorithm for estimating inter-chromosomal linkage disequilibrium (Under review).</a>")),
          tags$br(),
          tags$p(HTML(paste("Git version:", gTag[1,1])))
      
        )
      )
)

# Define server logic required to draw a histogram

server <- function(input, output, session) {
  # Plot on the web
  currentFile <- reactive({
    withProgress(message="X-LD:", value=0, {
      incProgress(1/3, detail = paste0(" check filesets ..."))
      FileLoad=0
      str=""
      if(length(which(grepl("*.bed", input$file_input$name)))  != 1) {
        str=paste(str, "No bed file found.") 
      } else {
        FileLoad=FileLoad+1
      }
      
      if(length(which(grepl("*.bim", input$file_input$name)))  != 1) {
        str=paste(str, "\nNo bim file found.")
      } else {
        FileLoad=FileLoad+1
      }
      
      if(length(which(grepl("*.fam", input$file_input$name)))  != 1) {
        str=paste(str, "\nNo fam file found.")
      } else {
        FileLoad=FileLoad+1
      }
      
      if (FileLoad < 3) {
        showNotification(str, duration = 5, type = "error")
        return()
      } else if (FileLoad > 3) {
        showNotification("More than 3 files selected", duration = 5, type="error")
      }
      
      idx=grep(".bed$", input$file_input$datapath)
      if (length(idx)==1) {
        rt=substr(input$file_input$datapath[idx], 1, nchar(input$file_input$datapath[idx])-4)
      }
      for (i in 1:3) {
        if (i != idx) {
          f1 = input$file_input$datapath[i]
          tl = substr(f1, nchar(f1)-2, nchar(f1))
          file.symlink(f1, paste0(rt, ".", tl))
        }
      }
            
      incProgress(1/3, detail = paste0(" check chromosome ..."))
      froot = substr(input$file_input$datapath[idx], 1, nchar(input$file_input$datapath[idx])-4)
      get_chr = read.table(paste0(froot,'.bim'),header=F,colClasses = c("character","NULL","NULL","NULL","NULL","NULL"))
      if (length(which(is.na(as.numeric(get_chr[,1]))))>0){
        showNotification("The chromosome index in the .bim file must be numeric!", duration = 5, type="error")
        stop("The chromosome index in the .bim file must be numeric! Refresh to continue.")
      }
      
      return (froot)
    })
  })
  
  mark <- gsub('[-: ]','',as.character(Sys.time()))
  
  observeEvent(input$run, {
    updateNavbarPage(session, "inNavbar", selected = "visualization")
    updateTabsetPanel(session, "X-LDFunctions","X-LD")
    
    autosome=as.numeric(input$autosome)
    
    froot = currentFile()
    withProgress(message="X-LD:", value=0, {
      n=4  
      Chr_Me_Matr <- data.frame(matrix(NA,nrow=autosome,ncol=autosome))
      rownames(Chr_Me_Matr)[1:autosome] <- c(paste0("chr",seq(1,autosome,1)))
      colnames(Chr_Me_Matr)[1:autosome] <- c(paste0("chr",seq(1,autosome,1)))
      # Chr marker
      Chr_Mark_Matr <- data.frame(matrix(NA,nrow=autosome,ncol=autosome))
      rownames(Chr_Mark_Matr)[1:autosome] <- c(paste0("chr",seq(1,autosome,1)))
      colnames(Chr_Mark_Matr)[1:autosome] <- c(paste0("chr",seq(1,autosome,1)))
      # MAF+Geno+LD
      time1 = proc.time()
      incProgress(1/n, detail = paste0(" QC ..."))
      
      if (input$ld_windows=="" | input$step_size=="" | input$r_square==0){
        QC1 = paste0(plink2, " --bfile ", froot, " --chr-set ", autosome, " --allow-extra-chr --set-missing-var-ids @:# --autosome --snps-only --make-bed --out ", froot,".1.autosome.snp")
        QC2 = paste0(plink2, " --bfile ", froot, ".1.autosome.snp --chr-set ", autosome, " --maf ", input$maf_cut, " --geno ", input$geno_cut, " --make-bed --out ", froot,".3.core")
        cat("QC...\nExtract autosome SNP variants...\n\n")
        system(QC1)
        cat("\nQC...\nMAF and missing rate...\n\n")
        system(QC2)
        cat("\nQC...\nFinished...\n\n")
      }else{
        QC1 = paste0(plink2, " --bfile ", froot, " --chr-set ", autosome, " --allow-extra-chr --set-missing-var-ids @:# --autosome --snps-only --make-bed --out ", froot,".1.autosome.snp")
        QC2 = paste0(plink2, " --bfile ", froot, ".1.autosome.snp --chr-set ", autosome, " --maf ", input$maf_cut, " --geno ", input$geno_cut, " --make-bed --out ", froot,".2.maf.geno")
        cat("QC...\nExtract autosome SNP variants...\n\n")
        system(QC1)
        cat("\nQC...\nMAF and missing rate...\n\n")
        system(QC2)
        QC3 = paste0(plink2, " --bfile ", froot, ".2.maf.geno --indep-pairwise ",input$ld_windows," ",input$step_size," ",input$r_square," --out ", froot ,".2.maf.geno.LD", " --chr-set ",autosome)
        QC4 = paste0(plink2, " --bfile ", froot, ".2.maf.geno --extract ", froot,".2.maf.geno.LD.prune.in"," --make-bed --out ", froot,".3.core ", " --chr-set ",autosome)
        system(QC3)
        system(QC4)
        cat("\nQC...\nFinished...\n\n")
      }

      #sc=ifelse(input$bred == 'inbred', 2, 1)
      #nn<-nrow(read.table(paste0(froot, ".3.core.fam"), as.is = T, header = F, colClasses = c("character","NULL","NULL","NULL","NULL","NULL")))
      time2 = proc.time()
      time = (time2-time1)[3][[1]]
      print(paste0(froot,' takes ',time,' seconds to finish the QC step.'))
      
      incProgress(2/n, detail = paste0(" making grm and calculate chromosome level LD ..."))
      # GRM construct + Me calculate
      time3 = proc.time()
      nn<-nrow(read.table(paste0(froot, ".3.core.fam"), as.is = T, header = F, colClasses = c("character","NULL","NULL","NULL","NULL","NULL")))
      sc=ifelse(input$bred == 'inbred', 2, 1)
      offDiag_Matr <- data.frame(matrix(NA,nrow=nn*(nn-1)/2,ncol=autosome))
    
      for(i in 1:autosome){
        Chr_Mark_Num_cmd = paste0("cat ",froot,".3.core.bim | grep -w '^",i,"' | wc -l")
        Chr_Mark_Num = system(Chr_Mark_Num_cmd,intern = TRUE)
        Chr_Mark_Num <- as.numeric(Chr_Mark_Num)
        # Determine whether chromosomes exist
        if(Chr_Mark_Num==0){
          next
        }else{
          Chr_Mark_Matr[i,i] <- Chr_Mark_Num
          incProgress(2/n, detail = paste0(" making grm for chromosome ",i ," ..."))
          GRM_cmd = paste0(plink2, " --bfile ",froot,".3.core --chr ",i," --chr-set 90 --allow-extra-chr --allow-no-sex --make-grm-gz --out ",froot,".chr",i)
          system(GRM_cmd)
          gz=gzfile(paste0(froot,".chr",i, ".grm.gz"))
          grm=read.table(gz, as.is = T)
          offDiag_Matr[,i] = grm[grm[,1]!=grm[,2], 4]/sc
          offDiag_2 = (grm[grm[,1]!=grm[,2], 4]/sc)^2
          Me=1/mean(offDiag_2,na.rm=TRUE)
          Chr_Me_Matr[i,i] <- as.numeric(Me)
        }
      }
      time4 = proc.time()
      time = (time4-time3)[3][[1]]
      print(paste0(froot,' takes ',time,' seconds to finish the decomposition of me.'))
      for(i in 1:autosome){
        if(Chr_Mark_Matr[i,i]==0){
          next
        }else{
          SNP1=as.numeric(Chr_Mark_Matr[i,i])
          offDiag1 = offDiag_Matr[,i]
          for(j in 1:autosome){
            if(Chr_Mark_Matr[j,j]==0){
              next
            }else{
              if(i<j){
                SNP2=as.numeric(Chr_Mark_Matr[j,j])
                Chr_Mark_Matr[i,j]=SNP1+SNP2
                offDiag2 = offDiag_Matr[,j] 
                Weight_offDiag_2 = ((SNP1/(SNP1+SNP2))*offDiag1+(SNP2/(SNP1+SNP2))*offDiag2)^2
                Me=1/mean(Weight_offDiag_2,na.rm=TRUE)
                Chr_Me_Matr[i,j] <- as.numeric(Me)
              }
            }
          }
          
        }
      }      
      
      # Fill lower triangle
      Chr_Me_Matr[lower.tri(Chr_Me_Matr)] <- t(Chr_Me_Matr)[lower.tri(Chr_Me_Matr)]
      Chr_Mark_Matr[lower.tri(Chr_Mark_Matr)] <- t(Chr_Mark_Matr)[lower.tri(Chr_Mark_Matr)]
      # Remove empty chromosome
      Chr_Me_Matr <- Chr_Me_Matr[apply(Chr_Me_Matr,1,function(y) any(!is.na(y))),]
      Chr_Me_Matr <- Chr_Me_Matr[,apply(Chr_Me_Matr,2,function(y) any(!is.na(y)))]
      
      Chr_Mark_Matr <- Chr_Mark_Matr[apply(Chr_Mark_Matr,1,function(y) any(!is.na(y))),]
      Chr_Mark_Matr <- Chr_Mark_Matr[,apply(Chr_Mark_Matr,2,function(y) any(!is.na(y)))]
      
      adjustdata <- function(data) {
        data<-cbind(rownames(data),data)
      }      
      LD <- as.data.frame(matrix(NA,nrow=nrow(Chr_Me_Matr),ncol=nrow(Chr_Me_Matr)))
      colnames(LD) <- rownames(Chr_Me_Matr)
      rownames(LD) <- rownames(Chr_Me_Matr)
      for(i in 1:nrow(LD)){
        for(j in 1:ncol(LD)){
          # intra-chromosomal LD
          if(i==j){
            LD[i,i] <- (Chr_Mark_Matr[i,i]^2/Chr_Me_Matr[i,i]-Chr_Mark_Matr[i,i])/(Chr_Mark_Matr[i,i]*(Chr_Mark_Matr[i,i]-1))
          }else{
            # inter-chromosomal LD
            LD[i,j] <- (Chr_Mark_Matr[i,j]^2/Chr_Me_Matr[i,j]-Chr_Mark_Matr[i,i]^2/Chr_Me_Matr[i,i]-Chr_Mark_Matr[j,j]^2/Chr_Me_Matr[j,j])/(2*Chr_Mark_Matr[i,i]*Chr_Mark_Matr[j,j])
          }
        }
      }

      LD <- as.matrix(LD)
      LD[LD < 0] <- 1e-300
      LD_Final <- adjustdata(LD)
      colnames(LD_Final)[1] <- ""
      Output3=paste0(froot,"_Chromosome_Level_LD.txt")
      write.table(LD_Final,file=Output3,quote = FALSE,sep="\t",row.names = FALSE)
      # LD scale calculate
      LD_Scale <- as.data.frame(matrix(NA,nrow=nrow(Chr_Me_Matr),ncol=nrow(Chr_Me_Matr)))
      colnames(LD_Scale) <- rownames(Chr_Me_Matr)
      rownames(LD_Scale) <- rownames(Chr_Me_Matr)
      for(i in 1:nrow(LD_Scale)){
        for(j in 1:nrow(LD_Scale)){
          # within chr
          if(i==j){
            LD_Scale[i,j] <- 1
          }else{
            LD_Scale[i,j] <- LD[i,j]/(sqrt(LD[i,i])*sqrt(LD[j,j]))
          }
        }
      }
      

      LD_Scale <- as.matrix(LD_Scale)
      LD_Scale_Final <- adjustdata(LD_Scale)
      colnames(LD_Scale_Final)[1] <- ""
      Output4=paste0(froot,"_Chromosome_Level_LD_Scale.txt")
      write.table(LD_Scale_Final,file=Output4,quote = FALSE,sep="\t",row.names = FALSE)
    

    withProgress(message="X-LD complete! Visualizing:", value=0, {    
      Output5=paste0(froot,"_Chromosome_Level_LD.png")
      Output6=paste0(froot,"_Chromosome_Level_LD_Scale.png")
      # Log conversion
      Log_LD <- -log10(LD)          
      col1 <- colorRampPalette(c("steelblue","white"))
      col2 <- colorRampPalette(c("white","steelblue"))               
    
      incProgress(2/n, detail = paste0(" LD plot ... "))
      output$me <- renderPlot({
        corrplot(
          Log_LD,
          is.corr=FALSE,
          method = "color",
          type = "up",
          col = col1(20),
          order = "original",
          cl.length=5,
          addgrid.col="white",
          cl.pos = "r",
          mar=c(0, 0, 1, 0),
          tl.col="black",tl.cex = 0.5,cl.ratio = 0.2,
          title = expression(paste("Decomposition of ",m[e]," [",-log[10],"LD]"))
        )
      })
      
      png(file=Output5,width = 600, height = 600)  
      corrplot(
        Log_LD,
        is.corr=FALSE,
        method = "color",
        type = "up",
        col = col1(20),
        order = "original",
        cl.length=5,
        addgrid.col="white",
        cl.pos = "r",
        mar=c(0, 0, 1, 0),
        tl.col="black",tl.cex = 0.5,cl.ratio = 0.2,
        title = expression(paste("Decomposition of ",m[e]," [",-log[10],"LD]"))
      )        
      dev.off()
      
      incProgress(2/n, detail = paste0(" LD (scaled) plot ... "))
      output$me_scale <- renderPlot({
        corrplot(
          LD_Scale,
          is.corr=FALSE,
          method = "color",
          type = "up",
          col = col2(20),
          order = "original",
          cl.length=5,
          addgrid.col="white",
          cl.pos = "r",
          mar=c(0, 0, 1, 0),
          tl.col="black",tl.cex = 1,cl.ratio = 0.2,
          title = expression(paste("Decomposition of ",m[e]," [",-log[10],"LD (scaled)]"))
        )
      })
      
      png(file=Output6,width = 600, height = 600)
      corrplot(
        LD_Scale,
        is.corr=FALSE,
        method = "color",
        type = "up",
        col = col2(20),
        order = "original",
        cl.length=5,
        addgrid.col="white",
        cl.pos = "r",
        mar=c(0, 0, 1, 0),
        tl.col="black",tl.cex = 1,cl.ratio = 0.2,
        title = expression(paste("Decomposition of ",m[e]," [",-log[10],"LD (scaled)]"))
      )
      dev.off()  
      
      output$figure <- downloadHandler( 
        filename = function(){
          paste0('X-LD_Fig.zip')
        },
        content = function(file) {
          fs <- c(paste0(froot,"_Chromosome_Level_LD.png"),paste0(froot,"_Chromosome_Level_LD_Scale.png"))
          zip(zipfile = file, files = fs)
        }
      )
      output$LD <- downloadHandler(          
        filename = function(){
          paste0("X-LD_Table.zip")
        },
        content = function(file) {
          fs <- c(paste0(froot,"_Chromosome_Level_LD.txt"),paste0(froot,"_Chromosome_Level_LD_Scale.txt"))
          zip(zipfile = file, files = fs)
        }
      )                  
        
    })
              
  })
  })
}
# Run the application
shinyApp(ui = ui, server = server)
