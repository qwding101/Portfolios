library(shiny); library(shinythemes); library(plotly); library(DT)
shinyUI(
  navbarPage(theme = shinytheme("simplex"), # 採用的模板
  "DSR Project", # 標題
  
  tabPanel("Abstract",
    h1("Facebook同溫層初探", align = "center"),
    h3("計劃說明"),
    p("　　過去已有研究層探討政黨立場的同溫層現象 (Colleoni et al., 2014)，但未曾以議題為單位但討同溫層現象（如：是否挺同婚的人，往往也支持廢死）。因此，",
      span("本研究分析社群網站Facebook的資料，利用以下四個議題（統獨、廢死、同志婚姻、核能發電）的八個面向（挺廢死、反廢死、挺同、反同、擁核、反核、統一、獨立），將臺灣的臉書用戶利用群聚分析 (cluster analysis) 予以分群，觀察臺灣是否存在同溫層現象", style = "color:brown"),
      "。我們預期透過人們支持或反對各議題的立場，可以透過演算法區分出不同的「議題同溫層」，進而檢視各同溫層的特質。",
      style = "font-size: 14pt"),
    h3("問題摘要"),
    withTags({
      ul(
      li("台灣FB用戶真的有「議題」同溫層嗎？"),li("哪一個同溫層的人最多？"),
      li("同溫層真的很「同溫」嗎？"),li("各同溫層的特性為何？"), style = "font-size: 14pt"
      )
    })
    
           ),
  
  
  navbarMenu("Results",
    
    tabPanel("EDA",
      h1("Exploratory Data Analysis"),
        mainPanel(
          tabsetPanel(# position = "above",
                      tabPanel("同溫層人數分配",
                               p("已排除Cluster 2（無外顯意見，共336432人）之資料"),
                               plotlyOutput("distPlot")),
                      tabPanel("議題正反立場人數",
                               p("已排除Cluster 2（無外顯意見，共336432人）之資料"),
                               plotlyOutput("IssuePlot"))
                      #tabPanel("Summary", verbatimTextOutput("summary")), 
                      #tabPanel("Table", tableOutput("table"))
                    )
                  )
            ),
    
    tabPanel("3D Vis",
             h1("k-means Clustering"),
             p("以下分析皆已排除Cluster 2 (無外顯意見)之資料。"),
             sidebarLayout(
               # sidebarPanel(
               # 
               #   conditionalPanel(condition = "input.CA3D == '支持：同婚、臺獨、廢死'||input.CA3D == '支持：恐同、統一、挺死'",
               #        htmlOutput("CA3Dsidertext")
               #                  ),
               #   conditionalPanel(condition = "input.CA3D == '立場加總'",
               #        htmlOutput("CA3Dsidertext"),
               #        hr(),
               #        radioButtons("dimchange", label = "您想自由調整哪一種圖的維度？", 
               #                           c("支持","立場加總"), selected = "立場加總", inline = T),
               #        uiOutput("dimchoice")
               #                  ),
                 # conditionalPanel(condition = "input.CA3D == '支持：自由繪製'",
                 #      htmlOutput("CA3Dsidertext"),
                 #      hr(),
                 #      radioButtons("dimchange", label = "您想自由調整哪一種圖的維度？", 
                 #                         c("支持","立場加總"), selected = "支持", inline = T),
                 #      uiOutput("dimchoice")
                 #                 )
                 #            ),
               sidebarPanel(
                 htmlOutput("CA3Dsidertext"),
                 shiny::hr(),
                 radioButtons("dimchange", label = "您想自由調整哪一種圖的維度？（結果請點選「自由繪製」的欄位）", c("支持","立場加總"),
                              inline = T),
                 uiOutput("dimchoice")
                            ),
               mainPanel(
                 tabsetPanel(id = "CA3D", # position = "above",
                             tabPanel("支持：同婚、臺獨、廢死", value = "3DS", plotlyOutput("CA3D_S")),
                             tabPanel("支持：恐同、統一、挺死", value = "3DO", plotlyOutput("CA3D_O")),
                             tabPanel("支持（自由繪製）", value = "3DSOF", plotlyOutput("CA3D_SOF")),
                             tabPanel("立場加總（自由繪製）", value = "3DT", htmlOutput("CA3D_TF_Error"), plotlyOutput("CA3D_TF"))
                            )
                        )
               
                          )
             
            ),
    
    tabPanel("Features on Issues",
      h1("Feature of Each Cluster on Issues"),
      h2("哪一個同溫層最團結？"),
      sidebarLayout( # Generate a row with a sidebar
        sidebarPanel( # Define the sidebar with one input
          selectInput(inputId = "Cluster", label = "同溫層（在各議題的分數）",
                      choices=c("Cluster 1"="1", "Cluster 2"="2", "Cluster 3"="3",
                                "Cluster 4"="4", "Cluster 5"="5", "Cluster 6"="6",
                                "Cluster 7"="7", "Cluster 8"="8", "Cluster 9"="9",
                                "Cluster 10"="10"),
                      selected="Cluster1"),
          selectInput(inputId = "Issue", label = "該同溫層在 _____ 立場的Z分數分配",
                      choices=c("廢死 (EndPenalty)"="EndPenalty", "臺獨 (Independence)"="Independence", "同婚 (GayMarriage)"="GayMarriage",
                                "反核 (Antinuclear)"="Antinuclear", "挺死 (Penalty)"="Penalty", "統一 (Reunification)"="Reunification",
                                "恐同 (Homophobia)"="Homophobia", "擁核 (Nuclear)"="Nuclear"),
                      selected="反核 (Antinuclear)"),
          shiny::hr(),
          tableOutput("Cltable"),
          helpText("Facebook Data from 2016, March.")
                      ),
        
        mainPanel(
          tabsetPanel(id = "FeatIssue", # position = "above",
            tabPanel("各議題分數", value = "FeatIssue1",
                     p("紅線表95%信賴區間。", style = "color:red"),
                     p("分數 (Z-score)越高表示越支持該立場"),
                     plotlyOutput("Features", height = "450px")),
            tabPanel("分數分配", value = "FeatIssue2", plotlyOutput("CluFreq")),
            tabPanel("所有同溫層", value = "FeatIssue3", plotOutput("AllCluFreq")),
            tabPanel("所有同溫層（去除Cluster 2）", value = "FeatIssue4", plotOutput("CluFreqno2"))
            #tabPanel("Summary", verbatimTextOutput("summary")), 
            #tabPanel("Table", tableOutput("table"))
          )
               )
              )
            ),
    
     tabPanel("Top FanPages",
       h1("各同溫層按讚粉專排名"),
       htmlOutput("numtext"),
       mainPanel(
         tabsetPanel(id = "Fan", # position = "above",
                     tabPanel("反核", value = "1", DT::dataTableOutput("Fan_table1")),
                     tabPanel("無外顯意見", value = "2", DT::dataTableOutput("Fan_table2")),
                     tabPanel("同婚", value = "3", DT::dataTableOutput("Fan_table3")),
                     tabPanel("恐同", value = "4", DT::dataTableOutput("Fan_table4")),
                     tabPanel("積極挺死", value = "5", DT::dataTableOutput("Fan_table5")),
                     tabPanel("挺死", value = "6", DT::dataTableOutput("Fan_table6")),
                     tabPanel("臺獨", value = "7", DT::dataTableOutput("Fan_table7")),
                     tabPanel("統一", value = "8", DT::dataTableOutput("Fan_table8")),
                     tabPanel("擁核", value = "9", DT::dataTableOutput("Fan_table9")),
                     tabPanel("廢死", value = "10", DT::dataTableOutput("Fan_table10"))
                     # shiny::dataTableOutput suggests please use DT::renderDataTable() instead
                    )
              )
     
           ),
    
    tabPanel("Other Features",
             h1("各同溫層的其他特徵"),
             mainPanel(
               tabsetPanel(# position = "above",
                           tabPanel("政治",
                                    plotOutput("PotPlot")),
                           tabPanel("宗教",
                                    plotOutput("RelPlot")),
                           tabPanel("媒體",
                                    plotOutput("MedPlot")),
                           tabPanel("有機食物",
                                    plotOutput("OrgPlot"))
                           #tabPanel("Summary", verbatimTextOutput("summary")), 
                           #tabPanel("Table", tableOutput("table"))
               )
             )
    )
  )#, 
  #tabPanel("Team",
  #         h1("國立臺灣大學", align = "center"), h2("R與資料科學導論", align = "center"),
  #         h3("第八組", align = "center")#,
            # withTags(
            #   table(tr(td("心理四"),td("B02207004"),td("　莊昀軒"), align = "center"),
            #         tr(td("心理四"),td("B02207030"),td("　林庭羽"), align = "center"),
            #         tr(td("心理四"),td("B02207069"),td("　丁麒文"), align = "center"),
            #         tr(td("歷史四"),td("B02203057"),td("　簡韻真"), align = "center"),
            #         tr(td("經濟碩二"),td("R04323033"),td("　尤筱筑"), align = "center"),
            #         tr(td("經濟碩二"),td("R04323017"),td("　黃詩媛"), align = "center"),
            #         align = "center",style = "font-size: 14pt")
            # )
   #       )
  
)
)