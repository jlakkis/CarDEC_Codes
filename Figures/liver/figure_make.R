setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

files = list.files()
files = files[tools::file_ext(files) == 'csv']

data = matrix(nrow = 29, ncol = 4)

start = 1

for(file in files) {
  tmp = read.csv(file)[,-1]
  end = nrow(tmp) + start - 1
  
  colnames(data) = colnames(tmp)
  data[c(start:end),] = as.matrix(tmp)
  
  start = end + 1
}

colnames(data) = c("Time", "Memory", "Method", "Percent")
data = as.data.frame(data)
data$Time = as.numeric(data$Time)
data$Time = data$Time/3600
data$Percent = as.numeric(data$Percent)

argmax = function(x) which(x == max(x))

data = data[-argmax(data$Time),]

myplot = ggplot2::ggplot(data,ggplot2::aes_string(x = 'Percent', y = 'Time', color = "Method")) +
  ggplot2::geom_point() +
  ggplot2::geom_line()

myplot = myplot + ggplot2::xlab("Percent of Liver Data") + ggplot2::ylab("Time (Hours)") + ggplot2::theme(axis.title = ggplot2::element_text(size=18))

ggsave("scaleability_plot.tiff", height = 7.4, width = 10, units = "in", dpi = 300, myplot)