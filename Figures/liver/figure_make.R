setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

packages = c("ggplot2")

for (package in packages) {
  if(package %in% rownames(installed.packages()) == F) {
    install.packages(package)
  }
}

library(ggplot2)

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

myplot = ggplot(data,ggplot2::aes_string(x = 'Percent', y = 'Time', color = "Method")) +
  geom_point() +
  geom_line()

myplot = myplot + xlab("Percent of Liver Data") + ylab("Time (Hours)") + theme(axis.title = element_text(size=18))

ggsave("scaleability_plot.tiff", height = 7.4, width = 10, units = "in", dpi = 300, myplot)