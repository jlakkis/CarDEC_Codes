setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

files = list.files()
files = files[tools::file_ext(files) == 'csv']

data = matrix(nrow = 27, ncol = 4)

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
data$`Time` = as.numeric(levels(data$`Time`))[data$`Time`]
data$Time = data$Time/3600
data$`Percent` = as.numeric(levels(data$`Percent`))[data$`Percent`]

argmax = function(x) which(x == max(x))

data

data = data[-argmax(data$Time),]

myplot = ggplot(data,aes_string(x = 'Percent', y = 'Time', color = "Method")) +
  geom_point() +
  geom_line()

myplot = myplot + xlab("Percent of Liver Data") + ylab("Time (Hours)") + theme(axis.title = element_text(size=18))

ggsave("scaleability_plot.tiff", height = 7.4, width = 10, units = "in", dpi = 300, myplot)