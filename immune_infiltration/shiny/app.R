# RShiny App for Immune Infiltration Analysis Dashboard
# 免疫浸润分析 - 交互式展示工具

library(shiny)
library(DT)
library(glue)
library(bslib)

# 项目路径配置
PROJECT_DIR <- "/media/desk16/share/secure/immune_infiltration"
RESULT_DIR <- file.path(PROJECT_DIR, "result")
LOG_DIR <- file.path(PROJECT_DIR, "logs")

# 方法信息
methods_info <- list(
  "1" = list(name = "CIBERSORT", desc = "CIBERSORT免疫细胞反卷积", output_dir = "cibersort_output"),
  "2" = list(name = "EPIC", desc = "EPIC免疫浸润分析", output_dir = "epic_output"),
  "3" = list(name = "ESTIMATE", desc = "ESTIMATE免疫评分", output_dir = "estimate_output"),
  "4" = list(name = "IPS", desc = "IPS免疫表性评分", output_dir = "ips_output"),
  "5" = list(name = "MCPcounter", desc = "MCP-counter免疫细胞计数", output_dir = "mcpcounter_output"),
  "6" = list(name = "ssGSEA", desc = "ssGSEA基因集富集分析", output_dir = "ssgsea_output"),
  "7" = list(name = "TIMER", desc = "TIMER免疫浸润", output_dir = "timer_output"),
  "8" = list(name = "xCell", desc = "xCell细胞富集评分", output_dir = "xcell_output")
)

# UI
ui <- fluidPage(
  theme = bs_theme(bootswatch = "cosmo"),

  titlePanel("🧬 免疫浸润分析 Dashboard"),

  sidebarLayout(
    sidebarPanel(
      h4("📊 项目状态"),
      verbatimTextOutput("project_status"),

      hr(),

      h4("⚙️ 命令生成器"),
      selectInput("method_select", "选择分析方法:",
                  choices = names(methods_info),
                  selected = "1",
                  selectize = FALSE),

      textInput("input_file", "表达矩阵文件:", value = "expression.csv"),

      textInput("output_dir", "输出目录:", value = "result"),

      actionButton("gen_cmd", "生成运行命令", class = "btn-primary"),

      hr(),

      h4("📋 运行命令"),
      verbatimTextOutput("generated_cmd"),
      actionButton("copy_cmd", "复制命令", class = "btn-sm"),

      hr(),

      downloadButton("download_log", "下载最新日志")
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("📁 结果文件",
                 h4("输出文件结构"),
                 DTOutput("file_table")),

        tabPanel("📈 日志查看",
                 h4("最新运行日志"),
                 verbatimTextOutput("log_content")),

        tabPanel("🖼️ 图片预览",
                 selectInput("img_select", "选择图片:",
                             choices = NULL),
                 imageOutput("img_preview", height = "500px")),

        tabPanel("📊 统计信息",
                 fluidRow(
                   column(6, card(card_header("运行次数"), textOutput("run_count"))),
                   column(6, card(card_header("结果文件数"), textOutput("file_count")))
                 ),
                 h4("各方法结果统计"),
                 tableOutput("method_stats"))
      )
    )
  )
)

# Server
server <- function(input, output, session) {

  # 项目状态
  output$project_status <- renderText({
    status <- "✅ 项目已完成\n"
    status <- paste0(status, "📁 结果目录: ", RESULT_DIR, "\n")
    status <- paste0(status, "📝 日志目录: ", LOG_DIR)
    status
  })

  # 生成命令
  observeEvent(input$gen_cmd, {
    method_name <- methods_info[[input$method_select]]$name
    cmd <- glue('../run_immune_infiltration.sh \\
  --config ../config/config.ini \\
  -o {input$output_dir} \\
  --methods {input$method_select}')

    output$generated_cmd <- renderText(cmd)
  })

  # 文件列表
  output$file_table <- renderDT({
    files <- list.files(RESULT_DIR, recursive = TRUE, full.names = FALSE)
    df <- data.frame(
      路径 = files,
      类型 = ifelse(grepl("\\.(png|pdf|jpg)$", files), "图片", "数据")
    )
    datatable(df, options = list(pageLength = 20))
  })

  # 日志内容
  output$log_content <- renderText({
    log_files <- list.files(LOG_DIR, pattern = "\\.log$", full.names = TRUE)
    if (length(log_files) > 0) {
      latest_log <- tail(sort(log_files), 1)
      readLines(latest_log, n = 100)
    } else {
      "暂无日志"
    }
  })

  # 图片选择
  observe({
    img_files <- list.files(RESULT_DIR, pattern = "\\.(png|pdf|jpg)$", recursive = TRUE)
    updateSelectInput(session, "img_select", choices = img_files)
  })

  # 图片预览
  output$img_preview <- renderImage({
    req(input$img_select)
    img_path <- file.path(RESULT_DIR, input$img_select)

    # 根据文件类型返回
    if (grepl("\\.png$|\\.jpg$", img_path)) {
      list(src = img_path, contentType = "image/png", alt = "图片预览")
    } else if (grepl("\\.pdf$", img_path)) {
      # PDF 需要特殊处理
      list(src = img_path, contentType = "application/pdf", alt = "PDF预览")
    }
  }, deleteFile = FALSE)

  # 统计信息
  output$run_count <- renderText({
    length(list.files(LOG_DIR, pattern = "\\.log$"))
  })

  output$file_count <- renderText({
    length(list.files(RESULT_DIR, recursive = TRUE))
  })

  output$method_stats <- renderTable({
    data.frame(
      方法 = sapply(methods_info, function(x) x$name),
      描述 = sapply(methods_info, function(x) x$desc),
      状态 = "✅ 完成"
    )
  })

  # 下载日志
  output$download_log <- downloadHandler(
    filename = function() {
      paste0("immune_log_", Sys.Date(), ".txt")
    },
    content = function(file) {
      log_files <- list.files(LOG_DIR, pattern = "\\.log$", full.names = TRUE)
      if (length(log_files) > 0) {
        latest_log <- tail(sort(log_files), 1)
        file.copy(latest_log, file)
      }
    }
  )

  # 复制命令
  observeEvent(input$copy_cmd, {
    cmd <- isolate(output$generated_cmd())
    if (!is.null(cmd) && cmd != "") {
      writeLines(cmd, "temp_cmd.sh")
      system("xclip -selection clipboard < temp_cmd.sh")
      showNotification("命令已复制到剪贴板!", type = "message")
    }
  })
}

# Run
shinyApp(ui = ui, server = server)
