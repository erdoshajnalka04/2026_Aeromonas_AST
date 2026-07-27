
library(gplots)
data <- read.csv("C:/Users/erdos/Desktop/AST heatmap.csv",
                 sep = ";", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

colnames(data) <- trimws(colnames(data))
rownames(data) <- trimws(data[,1])
data <- data[,-1]
data <- as.data.frame(lapply(data, function(x) {
  x <- gsub(",", ".", x)
  as.numeric(x)
}), row.names = rownames(data))


breakpoints <- list(
  ERY = c(S = 10, I = NA, R = 9),
  AMP = c(S = 10, I = NA, R = 9),
  ATM = c(S = 10, I = NA, R = 9),
  COL = c(S = 10, I = NA, R = 9),
  VAN = c(S = 10, I = NA, R = 9),
  OXC = c(S = 10, I = NA, R = 9),
  MEM = c(S = 23, I = 20, R = 19),
  CHL = c(S = 18, I = 13, R = 12),
  CIP = c(S = 21, I = 16, R = 15),
  TET = c(S = 15, I = 12, R = 11),
  GMN = c(S = 15, I = 13, R = 12),
  SXT = c(S = 16, I = 11, R = 10)
)

classify_SIR <- function(value, bp) {
  if (is.na(value)) return("I")   
  
  
  if (is.na(bp["I"])) {
    if (value >= bp["S"]) return("S")
    return("R")
  }
  
  
  if (value >= bp["S"]) return("S")
  if (value <= bp["R"]) return("R")
  return("I")
}

sir_matrix <- data
for (ab in colnames(data)) {
  sir_matrix[, ab] <- vapply(
    data[, ab],
    classify_SIR,
    bp = breakpoints[[ab]],
    FUN.VALUE = character(1)
  )
}


numeric_matrix <- sir_matrix
numeric_matrix[numeric_matrix == "S"] <- 0
numeric_matrix[numeric_matrix == "I"] <- 1
numeric_matrix[numeric_matrix == "R"] <- 2

numeric_matrix <- apply(numeric_matrix, 2, as.numeric)
rownames(numeric_matrix) <- rownames(data)

numeric_matrix[is.na(numeric_matrix)] <- 1



my_colors <- c("lightgreen", "lightyellow", "pink")  # S, I, R


x11(width = 14, height = 14)

heatmap.2(numeric_matrix,
          trace = "none",
          dendrogram = "both",
          Rowv = TRUE,
          Colv = TRUE,
          col = my_colors,
          key = FALSE,
          density.info = "none",
          margins = c(10, 10),
          main = "Antibiotic resistance of water-borne Aeromonas samples",
          rowsep = 1:nrow(numeric_matrix),
          colsep = 1:ncol(numeric_matrix),
          sepcolor = "grey",
          sepwidth = c(0.003, 0.003)
)

legend("topleft",
       legend = c("S = Sensitive", "I = Intermedier", "R = Resistant"),
       fill = my_colors,
       border = "black",
       box.lwd = 0.7,
       cex = 1.0,
       title = "Legend")

