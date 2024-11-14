# given a sample name, x1, x2, y1, y2, output the selection formatted in xenium explorer
# NOTE: x and y are flipped in seurat
library(dplyr)
helper_selection <- function(name_selection, x1,x2,y1,y2){
  if(x1 > x2){
    tmp <- x1
    x1 <- x2
    x2 <- tmp
  }
  if(y1 > y2){
    tmp <- y1
    y1 <- y2
    y2 <- tmp
  }
  
  line1 <- paste0('#Selection name: ', name_selection)
  line2 <- paste0('#Area (µm^2): ', abs((x1 - x2)*(y1 - y2)) %>% format(., nsmall=2))
  line3 <- paste0('X,Y')
  line4 <- paste0(x1, ',', y1)
  line5 <- paste0(x2, ',', y1)
  line6 <- paste0(x2, ',', y2)
  line7 <- paste0(x1, ',', y2)
  line8 <- paste0(x1, ',', y1)
  cat(line1, line2, line3, line4, line5, line6, line7, line8, sep = '\n')
}


name_selection <- 'pouR1'
x1 <- 400
x2 <- 1000
y1 <- 1400
y2 <- 2000
helper_selection(name_selection = name_selection, x1 = x1, x2 = x2, y1 = y1, y2 = y2)




















