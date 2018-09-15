# Load packages
library("shiny")
library("shinythemes")
library("dplyr")
library("readr")
library("ape")
library("RColorBrewer")
library("phytools")
library("paleotree")
library("ggplot2")
library("cowplot")
library("viridis")
library("rmarkdown")

## Load functions
makeTransparent<-function(someColor,alpha=10){
    newColor<-col2rgb(someColor)
    apply(newColor,2,function(curcoldata){
        rgb(red=curcoldata[1],green=curcoldata[2],blue=curcoldata[3],
            alpha=alpha,maxColorValue=255)
    })
}

####    BRANCHING.TIMES COMPATIVEL COM EXTINÇÃO     ####
branching.times.with.extinct<-function(phy)
{
    n <- length(phy$tip.label)
    N <- dim(phy$edge)[1]
    xx <- numeric(phy$Nnode)
    interns <- which(phy$edge[, 2] > n)
    for (i in interns)
        xx[phy$edge[i, 2] - n] <- xx[phy$edge[i,1] - n] + phy$edge.length[i]
    ##parte nova, novo algoritimo para calcular qual o tempo maximo da arvore
    ntip<-Ntip(phy)
    sum<-numeric(ntip)
    index<-numeric()
    x<-numeric()
    for(i in 1:ntip)
    {
        node<-i
        while(node!=ntip+1)
        {
            index<-which(phy$edge[,2]==node)
            sum[i]<-sum[i] + phy$edge.length[index]
            node<-phy$edge[index,1]
        }
    }
    depth <- max(sum)
    ##volta a ser usado o codigo da função original
    xx <- depth - xx
    names(xx) <- if (is.null(phy$node.label)) 
        (n + 1):(n + phy$Nnode)
    else phy$node.label
    xx    
}

## Organizing fossil data for plot

fossil.for.plot <- function(x){
    first.last <- plyr::ldply(x, function(y){y$taxa.data})
    sampling.occs <- setNames(plyr::llply(x, function(y){y$sampling.times}), NULL)
    occs <- data.frame(samp.times = 0, taxon.id = 0) #matrix(NA, ncol = length(sampling.occs[which.max(unlist(lapply(sampling.occs, length)))][[1]]))
    for(i in 1:dim(first.last)[1]){
        #print(i)
        if(length(sampling.occs[[i]]) == 0){
            occs <- rbind(occs, data.frame(samp.times = c(first.last$orig.time[i], first.last$ext.time[i]), taxon.id = rep(i, 2)))
        } else {
            temp <- sampling.occs[[i]]
            occs <- rbind(occs, data.frame(samp.times = temp, taxon.id = rep(i, length(temp))))
        }
    }
    return(list(first.last, occs[-1, ]))
}

# Define UI
ui <- navbarPage(title = "Exploração de Dados", theme = shinytheme("flatly"),
                 tabPanel(title = "Filogenias Moleculares",
                          sidebarLayout(
                              sidebarPanel = sidebarPanel(
                                  sliderInput("lambda", label = h3("Taxa de Especiação"), min = 0, 
                                              max = 0.5, value = 0.2, step = 0.05),
                                  sliderInput("mu", label = h3("Taxa de Extinção"), min = 0, 
                                              max = 0.5, value = 0, step = 0.05),
                                  sliderInput("sim.time", label = h3("Tempo"), min = 0, 
                                              max = 40, value = 10),
                                  #checkboxInput("ltt", label = "Mostrar LTT Plot", value = FALSE),
                                  actionButton(inputId = "run.sim", label = "Rodar Simulação")#,
                                  # radioButtons('format', 'Formato do Documento', c('PDF', 'HTML', 'Word'),
                                  #              inline = TRUE, selected = 'PDF'),
                                  # downloadButton("report.phylo", "Gerar relatório")
                              ),
                              # Output: Description, lineplot, and reference
                              mainPanel = mainPanel(
                                  plotOutput(outputId = "phylo", height = "800px")
                              )
                          )
                 ),
                 tabPanel(title = "Registro Fóssil",
                          sidebarLayout(
                              sidebarPanel = sidebarPanel(
                                  sliderInput("lambda.fossil", label = h3("Taxa de Especiação"), min = 0, 
                                              max = 0.5, value = 0.2, step = 0.1),
                                  sliderInput("mu.fossil", label = h3("Taxa de Extinção"), min = 0, 
                                              max = 0.5, value = 0, step = 0.1),
                                  sliderInput("sim.time.fossil", label = h3("Número Total de Espécies"), min = 30, 
                                              max = 50, value = 30),
                                  sliderInput("completeness", label = h3("Preservação (%)"), min = 0.001, 
                                              max = 100, value = 100),
                                  radioButtons("firstLast", label = "Tipo de Gráfico", choices = list("Duração" = 1, "Ocorrências" = 2), selected = 1, inline = TRUE),
                                  actionButton(inputId = "run.fossil.sim", label = "Rodar Simulação")#,
                                  # radioButtons('format', 'Formato do Documento', c('PDF', 'HTML', 'Word'),
                                  #              inline = TRUE, selected = 'PDF'),
                                  # downloadButton("report.fossil", "Gerar relatório")
                              ),
                              # Output: Description, lineplot, and reference
                              mainPanel = mainPanel(
                                  fluidRow(
                                      column(12,
                                             fluidRow(
                                                 column(6, plotOutput(outputId = "fossil.plot", height = "800px")),
                                                 column(6,
                                                        plotOutput(outputId = "duration.hist", height = "400px"),
                                                        plotOutput(outputId = "div.tt", height = "400px"))
                                             )
                                      )
                                  )
                              )
                          )
                 )
)



server <- function(input, output) {    
    ## Phylo simulation
    lambda.sim <- reactive({
        input$lambda
    })
    mu.sim <- reactive({
        input$mu
    })
    sim.time <- reactive({
        input$sim.time
    })
    # lttPlot <- reactive({
    #     input$ltt
    # })
    sim.phylo.res <- reactiveValues(phy = NULL)
    observeEvent(
        input$run.sim, {
            sim.phylo.res$phy <- rlineage(lambda.sim(), mu.sim(), sim.time())
        })
    output$phylo <- renderPlot({
        if(is.null(sim.phylo.res$phy)){
            return()
        }
        par(mfrow = c(2, 2))
        plot(sim.phylo.res$phy, edge.col = makeTransparent("blue", alpha = 90), show.tip.label = FALSE, edge.width = 3);axisPhylo()
        ltt.plot(sim.phylo.res$phy, log = "y", ylab = "Número de Espécies", xlab = "Tempo", main = "Filogenia Completa")
        plot(drop.fossil(sim.phylo.res$phy), edge.col = makeTransparent("blue", alpha = 90), show.tip.label = FALSE, edge.width = 3);axisPhylo()
        ltt.plot(drop.fossil(sim.phylo.res$phy), log = "y", ylab = "Número de Espécies", xlab = "Tempo", main = "Filogenia Molecular")
        })
    output$report.phylo <- downloadHandler(
        filename = function(){
            paste("report", sep = ".", switch(input$format, PDF = 'pdf', HTML = 'html', Word = 'docx')
            )
        },
        content = function(file){
            out <- rmarkdown::render("report.Rmd", output_format = switch(input$format, PDF = pdf_document(), HTML = html_document(), Word = word_document()), output_file = filename())
            file.copy(out, file)
        },
        contentType = "application/pdf"
    )
    ## Fossil simulation
    lambda.fossil <- reactive({
        input$lambda.fossil
    })
    mu.fossil <- reactive({
        input$mu.fossil
    })
    sim.time.fossil <- reactive({
        input$sim.time.fossil
    })
    fossil.comp <- reactive({
        input$completeness/100 - 0.00000001
    })
    firstLast <- reactive({
        input$firstLast
    })
    fossil.sim <- reactiveValues(lin = NULL)
    observeEvent(
        input$run.fossil.sim, {
            sim.record <- simFossilRecord(p = lambda.fossil(), q = mu.fossil(), r = sProb2sRate(R = fossil.comp()), nTotalTaxa = c(sim.time.fossil() -5, sim.time.fossil() + 5))
            fossil.sim$lin <- fossil.for.plot(sim.record)
        })
    output$fossil.plot <- renderPlot({
        if(is.null(fossil.sim$lin)){
            return()
        }
        if(firstLast() == 1){
            ggplot(data = fossil.sim$lin[[1]]) +
                geom_segment(aes(x = -orig.time, xend = -ext.time, y = taxon.id, yend = taxon.id, colour = taxon.id), size = 2) +
                scale_colour_viridis(direction = -1) +
                theme(legend.position = "none") +
                labs(x = "Tempo (Milhões de Anos)", y = "Linhagem")
        } else {
            ggplot(data = fossil.sim$lin[[2]]) +
                geom_point(mapping = aes(x = -samp.times, y = taxon.id, colour = taxon.id), size = 1.5, alpha = 0.5) +
                scale_colour_viridis(direction = -1) +
                theme(legend.position = "none") +
                labs(x = "Tempo (Milhões de Anos)", y = "Linhagem") +
                xlim(-max(fossil.sim$lin[[1]]$orig.time), 0)
        }
    })
    output$duration.hist <- renderPlot({
        if(is.null(fossil.sim$lin)){
            return()
        }
        ggplot(data = fossil.sim$lin[[1]]) +
            geom_histogram(aes(x = orig.time - ext.time)) +
            theme(legend.position = "none") +
            labs(x = "Duração da Linhagem (Milhões de Anos)")
    })
    output$div.tt <- renderPlot({
        if(is.null(fossil.sim$lin)){
            return()
        }
        ggplot(data = aggregate(fossil.sim$lin[[2]]$samp.times, by = list(fossil.sim$lin[[2]]$taxon.id), FUN = function(x){diff(range(x))})) +
            geom_histogram(aes(x = x)) +
            theme(legend.position = "none") +
            labs(x = "Duração da Linhagem (Milhões de Anos)")
    })
    output$report.fossil <- downloadHandler(
        filename = function(){
            paste("report_fossil", sep = ".", switch(input$format, PDF = 'pdf', HTML = 'html', Word = 'docx')
            )
        },
        content = function(file){
            out <- render("report_fossil.Rmd", output_format = switch(input$format, PDF = pdf_document(), HTML = html_document(), Word = word_document()))
            file.copy(out, file)
        },
        contentType = "application/pdf"
    )
}

# Create Shiny object
app <- shinyApp(ui = ui, server = server)

runApp(app)


### gerar relatório
### summary da duração das espécies
### simular árvores com declínio para evidenciar padrões do ltt plot
