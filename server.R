library(shiny) 
library(shinyhelper) 
library(data.table) 
library(Matrix) 
library(DT) 
library(magrittr) 
library(ggplot2) 
library(ggrepel) 
library(hdf5r) 
library(ggdendro) 
library(gridExtra) 
RNAconf = readRDS("RNAconf.rds")
RNAdef  = readRDS("RNAdef.rds")
RNAgene = readRDS("RNAgene.rds")
RNAmeta = readRDS("RNAmeta.rds")



ATACconf = readRDS("ATACconf.rds")
ATACdef  = readRDS("ATACdef.rds")
ATACgene = readRDS("ATACgene.rds")
ATACmeta = readRDS("ATACmeta.rds")



Peaksconf = readRDS("Peaksconf.rds")
Peaksdef  = readRDS("Peaksdef.rds")
Peaksgene = readRDS("Peaksgene.rds")
Peaksmeta = readRDS("Peaksmeta.rds")



### Useful stuff 
# Colour palette 
cList = list(c("grey85","#FFF7EC","#FEE8C8","#FDD49E","#FDBB84", 
               "#FC8D59","#EF6548","#D7301F","#B30000","#7F0000"), 
             c("#4575B4","#74ADD1","#ABD9E9","#E0F3F8","#FFFFBF", 
               "#FEE090","#FDAE61","#F46D43","#D73027")[c(1,1:9,9)], 
             c("#FDE725","#AADC32","#5DC863","#27AD81","#21908C", 
               "#2C728E","#3B528B","#472D7B","#440154")) 
names(cList) = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple") 
 
# Panel sizes 
pList = c("400px", "600px", "800px") 
names(pList) = c("Small", "Medium", "Large") 
pList2 = c("500px", "700px", "900px") 
names(pList2) = c("Small", "Medium", "Large") 
pList3 = c("600px", "800px", "1000px") 
names(pList3) = c("Small", "Medium", "Large") 
sList = c(18,24,30) 
names(sList) = c("Small", "Medium", "Large") 
lList = c(5,6,7) 
names(lList) = c("Small", "Medium", "Large") 
 
# Function to extract legend 
g_legend <- function(a.gplot){  
  tmp <- ggplot_gtable(ggplot_build(a.gplot))  
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")  
  legend <- tmp$grobs[[leg]]  
  legend 
}  
 
# Plot theme 
sctheme <- function(base_size = 24, XYval = TRUE, Xang = 0, XjusH = 0.5){ 
  oupTheme = theme( 
    text =             element_text(size = base_size, family = "Helvetica"), 
    panel.background = element_rect(fill = "white", colour = NA), 
    axis.line =   element_line(colour = "black"), 
    axis.ticks =  element_line(colour = "black", size = base_size / 20), 
    axis.title =  element_text(face = "bold"), 
    axis.text =   element_text(size = base_size), 
    axis.text.x = element_text(angle = Xang, hjust = XjusH), 
    legend.position = "bottom", 
    legend.key =      element_rect(colour = NA, fill = NA) 
  ) 
  if(!XYval){ 
    oupTheme = oupTheme + theme( 
      axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
      axis.text.y = element_blank(), axis.ticks.y = element_blank()) 
  } 
  return(oupTheme) 
} 
 
### Common plotting functions 
# Plot cell information on dimred 
scDRcell <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inpsub1, inpsub2, 
                     inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt, inplab){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID, 
                       inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID),  
                   with = FALSE] 
  colnames(ggData) = c("X", "Y", "val", "sub") 
  rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y)) 
  bgCells = FALSE 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    bgCells = TRUE 
    ggData2 = ggData[!sub %in% inpsub2] 
    ggData = ggData[sub %in% inpsub2] 
  } 
  if(inpord == "Max-1st"){ 
    ggData = ggData[order(val)] 
  } else if(inpord == "Min-1st"){ 
    ggData = ggData[order(-val)] 
  } else if(inpord == "Random"){ 
    ggData = ggData[sample(nrow(ggData))] 
  } 
  
  # Do factoring if required 
  if(!is.na(inpConf[UI == inp1]$fCL)){ 
    ggCol = strsplit(inpConf[UI == inp1]$fCL, "\\|")[[1]] 
    names(ggCol) = levels(ggData$val) 
    ggLvl = levels(ggData$val)[levels(ggData$val) %in% unique(ggData$val)] 
    ggData$val = factor(ggData$val, levels = ggLvl) 
    ggCol = ggCol[ggLvl] 
  } 
 
  # Actual ggplot 
  ggOut = ggplot(ggData, aes(X, Y, color = val)) 
  if(bgCells){ 
    ggOut = ggOut + 
      geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 16) 
  } 
  ggOut = ggOut + 
    geom_point(size = inpsiz, shape = 16) + xlab(inpdrX) + ylab(inpdrY) + 
    sctheme(base_size = sList[inpfsz], XYval = inptxt) 
  if(is.na(inpConf[UI == inp1]$fCL)){ 
    ggOut = ggOut + scale_color_gradientn("", colours = cList[[inpcol]]) + 
      guides(color = guide_colorbar(barwidth = 15)) 
  } else { 
    sListX = min(nchar(paste0(levels(ggData$val), collapse = "")), 200) 
    sListX = 0.75 * (sList - (1.5 * floor(sListX/50))) 
    ggOut = ggOut + scale_color_manual("", values = ggCol) + 
      guides(color = guide_legend(override.aes = list(size = 5),  
                                  nrow = inpConf[UI == inp1]$fRow)) + 
      theme(legend.text = element_text(size = sListX[inpfsz])) 
    if(inplab){ 
      ggData3 = ggData[, .(X = mean(X), Y = mean(Y)), by = "val"] 
      lListX = min(nchar(paste0(ggData3$val, collapse = "")), 200) 
      lListX = lList - (0.25 * floor(lListX/50)) 
      ggOut = ggOut + 
        geom_text_repel(data = ggData3, aes(X, Y, label = val), 
                        color = "grey10", bg.color = "grey95", bg.r = 0.15, 
                        size = lListX[inpfsz], seed = 42) 
    } 
  } 
  if(inpasp == "Square") { 
    ggOut = ggOut + coord_fixed(ratio = rat) 
  } else if(inpasp == "Fixed") { 
    ggOut = ggOut + coord_fixed() 
  } 
  return(ggOut) 
} 
 
scDRnum <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2, 
                    inpH5, inpGene, inpsplt){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID), 
                   with = FALSE] 
  colnames(ggData) = c("group", "sub") 
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData$val2 = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
  ggData[val2 < 0]$val2 = 0 
  h5file$close_all() 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  
  # Split inp1 if necessary 
  if(is.na(inpConf[UI == inp1]$fCL)){ 
    if(inpsplt == "Quartile"){nBk = 4} 
    if(inpsplt == "Decile"){nBk = 10} 
    ggData$group = cut(ggData$group, breaks = nBk) 
  } 
  
  # Actual data.table 
  ggData$express = FALSE 
  ggData[val2 > 0]$express = TRUE 
  ggData1 = ggData[express == TRUE, .(nExpress = .N), by = "group"] 
  ggData = ggData[, .(nCells = .N), by = "group"] 
  ggData = ggData1[ggData, on = "group"] 
  ggData = ggData[, c("group", "nCells", "nExpress"), with = FALSE] 
  ggData[is.na(nExpress)]$nExpress = 0 
  ggData$pctExpress = 100 * ggData$nExpress / ggData$nCells 
  ggData = ggData[order(group)] 
  colnames(ggData)[3] = paste0(colnames(ggData)[3], "_", inp2) 
  return(ggData) 
} 
# Plot gene expression on dimred 
scDRgene <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inpsub1, inpsub2, 
                     inpH5, inpGene, 
                     inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID, 
                       inpConf[UI == inpsub1]$ID),  
                   with = FALSE] 
  colnames(ggData) = c("X", "Y", "sub") 
  rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y)) 
  
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData$val = h5data$read(args = list(inpGene[inp1], quote(expr=))) 
  ggData[val < 0]$val = 0 
  h5file$close_all() 
  bgCells = FALSE 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    bgCells = TRUE 
    ggData2 = ggData[!sub %in% inpsub2] 
    ggData = ggData[sub %in% inpsub2] 
  } 
  if(inpord == "Max-1st"){ 
    ggData = ggData[order(val)] 
  } else if(inpord == "Min-1st"){ 
    ggData = ggData[order(-val)] 
  } else if(inpord == "Random"){ 
    ggData = ggData[sample(nrow(ggData))] 
  } 
   
  # Actual ggplot 
  ggOut = ggplot(ggData, aes(X, Y, color = val)) 
  if(bgCells){ 
    ggOut = ggOut + 
      geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 16) 
  } 
  ggOut = ggOut + 
    geom_point(size = inpsiz, shape = 16) + xlab(inpdrX) + ylab(inpdrY) + 
    sctheme(base_size = sList[inpfsz], XYval = inptxt) +  
    scale_color_gradientn(inp1, colours = cList[[inpcol]]) + 
      guides(color = guide_colorbar(barwidth = 15)) 
  if(inpasp == "Square") { 
    ggOut = ggOut + coord_fixed(ratio = rat) 
  } else if(inpasp == "Fixed") { 
    ggOut = ggOut + coord_fixed() 
  } 
  return(ggOut) 
} 
 
# Plot gene coexpression on dimred 
bilinear <- function(x,y,xy,Q11,Q21,Q12,Q22){ 
  oup = (xy-x)*(xy-y)*Q11 + x*(xy-y)*Q21 + (xy-x)*y*Q12 + x*y*Q22 
  oup = oup / (xy*xy) 
  return(oup) 
} 
scDRcoex <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inp2, 
                     inpsub1, inpsub2, inpH5, inpGene, 
                     inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID, 
                       inpConf[UI == inpsub1]$ID),  
                   with = FALSE] 
  colnames(ggData) = c("X", "Y", "sub") 
  rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y)) 
  
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData$val1 = h5data$read(args = list(inpGene[inp1], quote(expr=))) 
  ggData[val1 < 0]$val1 = 0 
  ggData$val2 = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
  ggData[val2 < 0]$val2 = 0 
  h5file$close_all() 
  bgCells = FALSE 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    bgCells = TRUE 
    ggData2 = ggData[!sub %in% inpsub2] 
    ggData = ggData[sub %in% inpsub2] 
  } 
  
  # Generate coex color palette 
  cInp = strsplit(inpcol, "; ")[[1]] 
  if(cInp[1] == "Red (Gene1)"){ 
    c10 = c(255,0,0) 
  } else if(cInp[1] == "Orange (Gene1)"){ 
    c10 = c(255,140,0) 
  } else { 
    c10 = c(0,255,0) 
  } 
  if(cInp[2] == "Green (Gene2)"){ 
    c01 = c(0,255,0) 
  } else { 
    c01 = c(0,0,255) 
  } 
  c00 = c(217,217,217) ; c11 = c10 + c01 
  nGrid = 16; nPad = 2; nTot = nGrid + nPad * 2 
  gg = data.table(v1 = rep(0:nTot,nTot+1), v2 = sort(rep(0:nTot,nTot+1))) 
  gg$vv1 = gg$v1 - nPad ; gg[vv1 < 0]$vv1 = 0; gg[vv1 > nGrid]$vv1 = nGrid 
  gg$vv2 = gg$v2 - nPad ; gg[vv2 < 0]$vv2 = 0; gg[vv2 > nGrid]$vv2 = nGrid 
  gg$cR = bilinear(gg$vv1, gg$vv2, nGrid, c00[1], c10[1], c01[1], c11[1]) 
  gg$cG = bilinear(gg$vv1, gg$vv2, nGrid, c00[2], c10[2], c01[2], c11[2]) 
  gg$cB = bilinear(gg$vv1, gg$vv2, nGrid, c00[3], c10[3], c01[3], c11[3]) 
  gg$cMix = rgb(gg$cR, gg$cG, gg$cB, maxColorValue = 255) 
  gg = gg[, c("v1", "v2", "cMix")] 
  
  # Map colours 
  ggData$v1 = round(nTot * ggData$val1 / max(ggData$val1)) 
  ggData$v2 = round(nTot * ggData$val2 / max(ggData$val2)) 
  ggData$v0 = ggData$v1 + ggData$v2 
  ggData = gg[ggData, on = c("v1", "v2")] 
  if(inpord == "Max-1st"){ 
    ggData = ggData[order(v0)] 
  } else if(inpord == "Min-1st"){ 
    ggData = ggData[order(-v0)] 
  } else if(inpord == "Random"){ 
    ggData = ggData[sample(nrow(ggData))] 
  } 
  
  # Actual ggplot 
  ggOut = ggplot(ggData, aes(X, Y)) 
  if(bgCells){ 
    ggOut = ggOut + 
      geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 16) 
  } 
  ggOut = ggOut + 
    geom_point(size = inpsiz, shape = 16, color = ggData$cMix) + 
    xlab(inpdrX) + ylab(inpdrY) + 
    sctheme(base_size = sList[inpfsz], XYval = inptxt) + 
    scale_color_gradientn(inp1, colours = cList[[1]]) + 
    guides(color = guide_colorbar(barwidth = 15)) 
  if(inpasp == "Square") { 
    ggOut = ggOut + coord_fixed(ratio = rat) 
  } else if(inpasp == "Fixed") { 
    ggOut = ggOut + coord_fixed() 
  } 
  return(ggOut) 
} 
 
scDRcoexLeg <- function(inp1, inp2, inpcol, inpfsz){ 
  # Generate coex color palette 
  cInp = strsplit(inpcol, "; ")[[1]] 
  if(cInp[1] == "Red (Gene1)"){ 
    c10 = c(255,0,0) 
  } else if(cInp[1] == "Orange (Gene1)"){ 
    c10 = c(255,140,0) 
  } else { 
    c10 = c(0,255,0) 
  } 
  if(cInp[2] == "Green (Gene2)"){ 
    c01 = c(0,255,0) 
  } else { 
    c01 = c(0,0,255) 
  } 
  c00 = c(217,217,217) ; c11 = c10 + c01 
  nGrid = 16; nPad = 2; nTot = nGrid + nPad * 2 
  gg = data.table(v1 = rep(0:nTot,nTot+1), v2 = sort(rep(0:nTot,nTot+1))) 
  gg$vv1 = gg$v1 - nPad ; gg[vv1 < 0]$vv1 = 0; gg[vv1 > nGrid]$vv1 = nGrid 
  gg$vv2 = gg$v2 - nPad ; gg[vv2 < 0]$vv2 = 0; gg[vv2 > nGrid]$vv2 = nGrid 
  gg$cR = bilinear(gg$vv1, gg$vv2, nGrid, c00[1], c10[1], c01[1], c11[1]) 
  gg$cG = bilinear(gg$vv1, gg$vv2, nGrid, c00[2], c10[2], c01[2], c11[2]) 
  gg$cB = bilinear(gg$vv1, gg$vv2, nGrid, c00[3], c10[3], c01[3], c11[3]) 
  gg$cMix = rgb(gg$cR, gg$cG, gg$cB, maxColorValue = 255) 
  gg = gg[, c("v1", "v2", "cMix")] 
  
  # Actual ggplot 
  ggOut = ggplot(gg, aes(v1, v2)) + 
    geom_tile(fill = gg$cMix) + 
    xlab(inp1) + ylab(inp2) + coord_fixed(ratio = 1) + 
    scale_x_continuous(breaks = c(0, nTot), label = c("low", "high")) + 
    scale_y_continuous(breaks = c(0, nTot), label = c("low", "high")) + 
    sctheme(base_size = sList[inpfsz], XYval = TRUE) 
  return(ggOut) 
} 
 
scDRcoexNum <- function(inpConf, inpMeta, inp1, inp2, 
                        inpsub1, inpsub2, inpH5, inpGene){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inpsub1]$ID), with = FALSE] 
  colnames(ggData) = c("sub") 
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData$val1 = h5data$read(args = list(inpGene[inp1], quote(expr=))) 
  ggData[val1 < 0]$val1 = 0 
  ggData$val2 = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
  ggData[val2 < 0]$val2 = 0 
  h5file$close_all() 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  
  # Actual data.table 
  ggData$express = "none" 
  ggData[val1 > 0]$express = inp1 
  ggData[val2 > 0]$express = inp2 
  ggData[val1 > 0 & val2 > 0]$express = "both" 
  ggData$express = factor(ggData$express, levels = unique(c("both", inp1, inp2, "none"))) 
  ggData = ggData[, .(nCells = .N), by = "express"] 
  ggData$percent = 100 * ggData$nCells / sum(ggData$nCells) 
  ggData = ggData[order(express)] 
  colnames(ggData)[1] = "expression > 0" 
  return(ggData) 
} 
 
# Plot violin / boxplot 
scVioBox <- function(inpConf, inpMeta, inp1, inp2, 
                     inpsub1, inpsub2, inpH5, inpGene, 
                     inptyp, inppts, inpsiz, inpfsz){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID), 
                   with = FALSE] 
  colnames(ggData) = c("X", "sub") 
  
  # Load in either cell meta or gene expr
  if(inp2 %in% inpConf$UI){ 
    ggData$val = inpMeta[[inpConf[UI == inp2]$ID]] 
  } else { 
    h5file <- H5File$new(inpH5, mode = "r") 
    h5data <- h5file[["grp"]][["data"]] 
    ggData$val = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
    ggData[val < 0]$val = 0 
    set.seed(42) 
    tmpNoise = rnorm(length(ggData$val)) * diff(range(ggData$val)) / 1000 
    ggData$val = ggData$val + tmpNoise 
    h5file$close_all() 
  } 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  
  # Do factoring 
  ggCol = strsplit(inpConf[UI == inp1]$fCL, "\\|")[[1]] 
  names(ggCol) = levels(ggData$X) 
  ggLvl = levels(ggData$X)[levels(ggData$X) %in% unique(ggData$X)] 
  ggData$X = factor(ggData$X, levels = ggLvl) 
  ggCol = ggCol[ggLvl] 
  
  # Actual ggplot 
  if(inptyp == "violin"){ 
    ggOut = ggplot(ggData, aes(X, val, fill = X)) + geom_violin(scale = "width") 
  } else { 
    ggOut = ggplot(ggData, aes(X, val, fill = X)) + geom_boxplot() 
  } 
  if(inppts){ 
    ggOut = ggOut + geom_jitter(size = inpsiz, shape = 16) 
  } 
  ggOut = ggOut + xlab(inp1) + ylab(inp2) + 
    sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +  
    scale_fill_manual("", values = ggCol) +
    theme(legend.position = "none")
  return(ggOut) 
} 
 
# Plot proportion plot 
scProp <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2, 
                   inptyp, inpflp, inpfsz){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inp2]$ID, 
                       inpConf[UI == inpsub1]$ID),  
                   with = FALSE] 
  colnames(ggData) = c("X", "grp", "sub") 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  ggData = ggData[, .(nCells = .N), by = c("X", "grp")] 
  ggData = ggData[, {tot = sum(nCells) 
                      .SD[,.(pctCells = 100 * sum(nCells) / tot, 
                             nCells = nCells), by = "grp"]}, by = "X"] 
  
  # Do factoring 
  ggCol = strsplit(inpConf[UI == inp2]$fCL, "\\|")[[1]] 
  names(ggCol) = levels(ggData$grp) 
  ggLvl = levels(ggData$grp)[levels(ggData$grp) %in% unique(ggData$grp)] 
  ggData$grp = factor(ggData$grp, levels = ggLvl) 
  ggCol = ggCol[ggLvl] 
  
  # Actual ggplot 
  if(inptyp == "Proportion"){ 
    ggOut = ggplot(ggData, aes(X, pctCells, fill = grp)) + 
      geom_col() + ylab("Cell Proportion (%)") 
  } else { 
    ggOut = ggplot(ggData, aes(X, nCells, fill = grp)) + 
      geom_col() + ylab("Number of Cells") 
  } 
  if(inpflp){ 
    ggOut = ggOut + coord_flip() 
  } 
  ggOut = ggOut + xlab(inp1) + 
    sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +  
    scale_fill_manual("", values = ggCol) + 
    theme(legend.position = "right") 
  return(ggOut) 
} 
 
# Get gene list 
scGeneList <- function(inp, inpGene){ 
  geneList = data.table(gene = unique(trimws(strsplit(inp, ",|;|
")[[1]])), 
                        present = TRUE) 
  geneList[!gene %in% names(inpGene)]$present = FALSE 
  return(geneList) 
} 
 
# Plot gene expression bubbleplot / heatmap 
scBubbHeat <- function(inpConf, inpMeta, inp, inpGrp, inpPlt, 
                       inpsub1, inpsub2, inpH5, inpGene, inpScl, inpRow, inpCol, 
                       inpcols, inpfsz, save = FALSE){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Identify genes that are in our dataset 
  geneList = scGeneList(inp, inpGene) 
  geneList = geneList[present == TRUE] 
  shiny::validate(need(nrow(geneList) <= 50, "More than 50 genes to plot! Please reduce the gene list!")) 
  shiny::validate(need(nrow(geneList) > 1, "Please input at least 2 genes to plot!")) 
   
  # Prepare ggData 
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData = data.table() 
  for(iGene in geneList$gene){ 
    tmp = inpMeta[, c("sampleID", inpConf[UI == inpsub1]$ID), with = FALSE] 
    colnames(tmp) = c("sampleID", "sub") 
    tmp$grpBy = inpMeta[[inpConf[UI == inpGrp]$ID]] 
    tmp$geneName = iGene 
    tmp$val = h5data$read(args = list(inpGene[iGene], quote(expr=))) 
    ggData = rbindlist(list(ggData, tmp)) 
  } 
  h5file$close_all() 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  shiny::validate(need(uniqueN(ggData$grpBy) > 1, "Only 1 group present, unable to plot!")) 
   
  # Aggregate 
  ggData$val = expm1(ggData$val) 
  ggData = ggData[, .(val = mean(val), prop = sum(val>0) / length(sampleID)), 
                  by = c("geneName", "grpBy")] 
  ggData$val = log1p(ggData$val) 
   
  # Scale if required 
  colRange = range(ggData$val) 
  if(inpScl){ 
    ggData[, val:= scale(val), keyby = "geneName"] 
    colRange = c(-max(abs(range(ggData$val))), max(abs(range(ggData$val)))) 
  } 
   
  # hclust row/col if necessary 
  ggMat = dcast.data.table(ggData, geneName~grpBy, value.var = "val") 
  tmp = ggMat$geneName 
  ggMat = as.matrix(ggMat[, -1]) 
  rownames(ggMat) = tmp 
  if(inpRow){ 
    hcRow = dendro_data(as.dendrogram(hclust(dist(ggMat)))) 
    ggRow = ggplot() + coord_flip() + 
      geom_segment(data = hcRow$segments, aes(x=x,y=y,xend=xend,yend=yend)) + 
      scale_y_continuous(breaks = rep(0, uniqueN(ggData$grpBy)), 
                         labels = unique(ggData$grpBy), expand = c(0, 0)) + 
      scale_x_continuous(breaks = seq_along(hcRow$labels$label), 
                         labels = hcRow$labels$label, expand = c(0, 0.5)) + 
      sctheme(base_size = sList[inpfsz]) + 
      theme(axis.title = element_blank(), axis.line = element_blank(), 
            axis.ticks = element_blank(), axis.text.y = element_blank(), 
            axis.text.x = element_text(color="white", angle = 45, hjust = 1)) 
    ggData$geneName = factor(ggData$geneName, levels = hcRow$labels$label) 
  } else { 
    ggData$geneName = factor(ggData$geneName, levels = rev(geneList$gene)) 
  } 
  if(inpCol){ 
    hcCol = dendro_data(as.dendrogram(hclust(dist(t(ggMat))))) 
    ggCol = ggplot() + 
      geom_segment(data = hcCol$segments, aes(x=x,y=y,xend=xend,yend=yend)) + 
      scale_x_continuous(breaks = seq_along(hcCol$labels$label), 
                         labels = hcCol$labels$label, expand = c(0.05, 0)) + 
      scale_y_continuous(breaks = rep(0, uniqueN(ggData$geneName)), 
                         labels = unique(ggData$geneName), expand=c(0,0)) + 
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) + 
      theme(axis.title = element_blank(), axis.line = element_blank(), 
            axis.ticks = element_blank(), axis.text.x = element_blank(), 
            axis.text.y = element_text(color = "white")) 
    ggData$grpBy = factor(ggData$grpBy, levels = hcCol$labels$label) 
  } 
   
  # Actual plot according to plottype 
  if(inpPlt == "Bubbleplot"){ 
    # Bubbleplot 
    ggOut = ggplot(ggData, aes(grpBy, geneName, color = val, size = prop)) + 
      geom_point() +  
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +  
      scale_x_discrete(expand = c(0.05, 0)) +  
      scale_y_discrete(expand = c(0, 0.5)) + 
      scale_size_continuous("proportion", range = c(0, 8), 
                            limits = c(0, 1), breaks = c(0.00,0.25,0.50,0.75,1.00)) + 
      scale_color_gradientn("expression", limits = colRange, colours = cList[[inpcols]]) + 
      guides(color = guide_colorbar(barwidth = 15)) + 
      theme(axis.title = element_blank(), legend.box = "vertical") 
  } else { 
    # Heatmap 
    ggOut = ggplot(ggData, aes(grpBy, geneName, fill = val)) + 
      geom_tile() +  
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) + 
      scale_x_discrete(expand = c(0.05, 0)) +  
      scale_y_discrete(expand = c(0, 0.5)) + 
      scale_fill_gradientn("expression", limits = colRange, colours = cList[[inpcols]]) + 
      guides(fill = guide_colorbar(barwidth = 15)) + 
      theme(axis.title = element_blank()) 
  } 
     
  # Final tidy 
  ggLeg = g_legend(ggOut) 
  ggOut = ggOut + theme(legend.position = "none") 
  if(!save){ 
    if(inpRow & inpCol){ggOut =  
      grid.arrange(ggOut, ggLeg, ggCol, ggRow, widths = c(7,1), heights = c(1,7,2),  
                   layout_matrix = rbind(c(3,NA),c(1,4),c(2,NA)))  
    } else if(inpRow){ggOut =  
      grid.arrange(ggOut, ggLeg, ggRow, widths = c(7,1), heights = c(7,2),  
                   layout_matrix = rbind(c(1,3),c(2,NA)))  
    } else if(inpCol){ggOut =  
      grid.arrange(ggOut, ggLeg, ggCol, heights = c(1,7,2),  
                   layout_matrix = rbind(c(3),c(1),c(2)))  
    } else {ggOut =  
      grid.arrange(ggOut, ggLeg, heights = c(7,2),  
                   layout_matrix = rbind(c(1),c(2)))  
    }  
  } else { 
    if(inpRow & inpCol){ggOut =  
      arrangeGrob(ggOut, ggLeg, ggCol, ggRow, widths = c(7,1), heights = c(1,7,2),  
                  layout_matrix = rbind(c(3,NA),c(1,4),c(2,NA)))  
    } else if(inpRow){ggOut =  
      arrangeGrob(ggOut, ggLeg, ggRow, widths = c(7,1), heights = c(7,2),  
                  layout_matrix = rbind(c(1,3),c(2,NA)))  
    } else if(inpCol){ggOut =  
      arrangeGrob(ggOut, ggLeg, ggCol, heights = c(1,7,2),  
                  layout_matrix = rbind(c(3),c(1),c(2)))  
    } else {ggOut =  
      arrangeGrob(ggOut, ggLeg, heights = c(7,2),  
                  layout_matrix = rbind(c(1),c(2)))  
    }  
  } 
  return(ggOut) 
} 
 
 
 
 
 
### Start server code 
shinyServer(function(input, output, session) { 
  ### For all tags and Server-side selectize 
  observe_helpers() 
 optCrt="{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }" 
  updateSelectizeInput(session, "RNAa1inp2", choices = names(RNAgene), server = TRUE, 
                       selected = RNAdef$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "RNAa3inp1", choices = names(RNAgene), server = TRUE, 
                       selected = RNAdef$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "RNAa3inp2", choices = names(RNAgene), server = TRUE, 
                       selected = RNAdef$gene2, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "RNAb2inp1", choices = names(RNAgene), server = TRUE, 
                       selected = RNAdef$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "RNAb2inp2", choices = names(RNAgene), server = TRUE, 
                       selected = RNAdef$gene2, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "RNAc1inp2", server = TRUE, 
                       choices = c(RNAconf[is.na(fID)]$UI,names(RNAgene)), 
                       selected = RNAconf[is.na(fID)]$UI[1], options = list( 
                         maxOptions = length(RNAconf[is.na(fID)]$UI) + 3, 
                         create = TRUE, persist = TRUE, render = I(optCrt))) 
 
 
  ### Plots for tab a1 
  output$RNAa1sub1.ui <- renderUI({ 
    sub = strsplit(RNAconf[UI == input$RNAa1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("RNAa1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$RNAa1sub1non, { 
    sub = strsplit(RNAconf[UI == input$RNAa1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "RNAa1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$RNAa1sub1all, { 
    sub = strsplit(RNAconf[UI == input$RNAa1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "RNAa1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$RNAa1oup1 <- renderPlot({ 
    scDRcell(RNAconf, RNAmeta, input$RNAa1drX, input$RNAa1drY, input$RNAa1inp1,  
             input$RNAa1sub1, input$RNAa1sub2, 
             input$RNAa1siz, input$RNAa1col1, input$RNAa1ord1, 
             input$RNAa1fsz, input$RNAa1asp, input$RNAa1txt, input$RNAa1lab1) 
  }) 
  output$RNAa1oup1.ui <- renderUI({ 
    plotOutput("RNAa1oup1", height = pList[input$RNAa1psz]) 
  }) 
  output$RNAa1oup1.pdf <- downloadHandler( 
    filename = function() { paste0("RNA",input$RNAa1drX,"_",input$RNAa1drY,"_",  
                                   input$RNAa1inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$RNAa1oup1.h, width = input$RNAa1oup1.w, useDingbats = FALSE, 
      plot = scDRcell(RNAconf, RNAmeta, input$RNAa1drX, input$RNAa1drY, input$RNAa1inp1,   
                      input$RNAa1sub1, input$RNAa1sub2, 
                      input$RNAa1siz, input$RNAa1col1, input$RNAa1ord1,  
                      input$RNAa1fsz, input$RNAa1asp, input$RNAa1txt, input$RNAa1lab1) ) 
  }) 
  output$RNAa1oup1.png <- downloadHandler( 
    filename = function() { paste0("RNA",input$RNAa1drX,"_",input$RNAa1drY,"_",  
                                   input$RNAa1inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$RNAa1oup1.h, width = input$RNAa1oup1.w, 
      plot = scDRcell(RNAconf, RNAmeta, input$RNAa1drX, input$RNAa1drY, input$RNAa1inp1,   
                      input$RNAa1sub1, input$RNAa1sub2, 
                      input$RNAa1siz, input$RNAa1col1, input$RNAa1ord1,  
                      input$RNAa1fsz, input$RNAa1asp, input$RNAa1txt, input$RNAa1lab1) ) 
  }) 
  output$RNAa1.dt <- renderDataTable({ 
    ggData = scDRnum(RNAconf, RNAmeta, input$RNAa1inp1, input$RNAa1inp2, 
                     input$RNAa1sub1, input$RNAa1sub2, 
                     "RNAgexpr.h5", RNAgene, input$RNAa1splt) 
    datatable(ggData, rownames = FALSE, extensions = "Buttons", 
              options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>% 
      formatRound(columns = c("pctExpress"), digits = 2) 
  }) 
   
  output$RNAa1oup2 <- renderPlot({ 
    scDRgene(RNAconf, RNAmeta, input$RNAa1drX, input$RNAa1drY, input$RNAa1inp2,  
             input$RNAa1sub1, input$RNAa1sub2, 
             "RNAgexpr.h5", RNAgene, 
             input$RNAa1siz, input$RNAa1col2, input$RNAa1ord2, 
             input$RNAa1fsz, input$RNAa1asp, input$RNAa1txt) 
  }) 
  output$RNAa1oup2.ui <- renderUI({ 
    plotOutput("RNAa1oup2", height = pList[input$RNAa1psz]) 
  }) 
  output$RNAa1oup2.pdf <- downloadHandler( 
    filename = function() { paste0("RNA",input$RNAa1drX,"_",input$RNAa1drY,"_",  
                                   input$RNAa1inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$RNAa1oup2.h, width = input$RNAa1oup2.w, useDingbats = FALSE, 
      plot = scDRgene(RNAconf, RNAmeta, input$RNAa1drX, input$RNAa1drY, input$RNAa1inp2,  
                      input$RNAa1sub1, input$RNAa1sub2, 
                      "RNAgexpr.h5", RNAgene, 
                      input$RNAa1siz, input$RNAa1col2, input$RNAa1ord2, 
                      input$RNAa1fsz, input$RNAa1asp, input$RNAa1txt) ) 
  }) 
  output$RNAa1oup2.png <- downloadHandler( 
    filename = function() { paste0("RNA",input$RNAa1drX,"_",input$RNAa1drY,"_",  
                                   input$RNAa1inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$RNAa1oup2.h, width = input$RNAa1oup2.w, 
      plot = scDRgene(RNAconf, RNAmeta, input$RNAa1drX, input$RNAa1drY, input$RNAa1inp2,  
                      input$RNAa1sub1, input$RNAa1sub2, 
                      "RNAgexpr.h5", RNAgene, 
                      input$RNAa1siz, input$RNAa1col2, input$RNAa1ord2, 
                      input$RNAa1fsz, input$RNAa1asp, input$RNAa1txt) ) 
  }) 
   
   
  ### Plots for tab a2 
  output$RNAa2sub1.ui <- renderUI({ 
    sub = strsplit(RNAconf[UI == input$RNAa2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("RNAa2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$RNAa2sub1non, { 
    sub = strsplit(RNAconf[UI == input$RNAa2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "RNAa2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$RNAa2sub1all, { 
    sub = strsplit(RNAconf[UI == input$RNAa2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "RNAa2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$RNAa2oup1 <- renderPlot({ 
    scDRcell(RNAconf, RNAmeta, input$RNAa2drX, input$RNAa2drY, input$RNAa2inp1,  
             input$RNAa2sub1, input$RNAa2sub2, 
             input$RNAa2siz, input$RNAa2col1, input$RNAa2ord1, 
             input$RNAa2fsz, input$RNAa2asp, input$RNAa2txt, input$RNAa2lab1) 
  }) 
  output$RNAa2oup1.ui <- renderUI({ 
    plotOutput("RNAa2oup1", height = pList[input$RNAa2psz]) 
  }) 
  output$RNAa2oup1.pdf <- downloadHandler( 
    filename = function() { paste0("RNA",input$RNAa2drX,"_",input$RNAa2drY,"_",  
                                   input$RNAa2inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$RNAa2oup1.h, width = input$RNAa2oup1.w, useDingbats = FALSE, 
      plot = scDRcell(RNAconf, RNAmeta, input$RNAa2drX, input$RNAa2drY, input$RNAa2inp1,   
                      input$RNAa2sub1, input$RNAa2sub2, 
                      input$RNAa2siz, input$RNAa2col1, input$RNAa2ord1,  
                      input$RNAa2fsz, input$RNAa2asp, input$RNAa2txt, input$RNAa2lab1) ) 
  }) 
  output$RNAa2oup1.png <- downloadHandler( 
    filename = function() { paste0("RNA",input$RNAa2drX,"_",input$RNAa2drY,"_",  
                                   input$RNAa2inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$RNAa2oup1.h, width = input$RNAa2oup1.w, 
      plot = scDRcell(RNAconf, RNAmeta, input$RNAa2drX, input$RNAa2drY, input$RNAa2inp1,   
                      input$RNAa2sub1, input$RNAa2sub2, 
                      input$RNAa2siz, input$RNAa2col1, input$RNAa2ord1,  
                      input$RNAa2fsz, input$RNAa2asp, input$RNAa2txt, input$RNAa2lab1) ) 
  }) 
   
  output$RNAa2oup2 <- renderPlot({ 
    scDRcell(RNAconf, RNAmeta, input$RNAa2drX, input$RNAa2drY, input$RNAa2inp2,  
             input$RNAa2sub1, input$RNAa2sub2, 
             input$RNAa2siz, input$RNAa2col2, input$RNAa2ord2, 
             input$RNAa2fsz, input$RNAa2asp, input$RNAa2txt, input$RNAa2lab2) 
  }) 
  output$RNAa2oup2.ui <- renderUI({ 
    plotOutput("RNAa2oup2", height = pList[input$RNAa2psz]) 
  }) 
  output$RNAa2oup2.pdf <- downloadHandler( 
    filename = function() { paste0("RNA",input$RNAa2drX,"_",input$RNAa2drY,"_",  
                                   input$RNAa2inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$RNAa2oup2.h, width = input$RNAa2oup2.w, useDingbats = FALSE, 
      plot = scDRcell(RNAconf, RNAmeta, input$RNAa2drX, input$RNAa2drY, input$RNAa2inp2,   
                      input$RNAa2sub1, input$RNAa2sub2, 
                      input$RNAa2siz, input$RNAa2col2, input$RNAa2ord2,  
                      input$RNAa2fsz, input$RNAa2asp, input$RNAa2txt, input$RNAa2lab2) ) 
  }) 
  output$RNAa2oup2.png <- downloadHandler( 
    filename = function() { paste0("RNA",input$RNAa2drX,"_",input$RNAa2drY,"_",  
                                   input$RNAa2inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$RNAa2oup2.h, width = input$RNAa2oup2.w, 
      plot = scDRcell(RNAconf, RNAmeta, input$RNAa2drX, input$RNAa2drY, input$RNAa2inp2,   
                      input$RNAa2sub1, input$RNAa2sub2, 
                      input$RNAa2siz, input$RNAa2col2, input$RNAa2ord2,  
                      input$RNAa2fsz, input$RNAa2asp, input$RNAa2txt, input$RNAa2lab2) ) 
  }) 
   
   
  ### Plots for tab a3 
  output$RNAa3sub1.ui <- renderUI({ 
    sub = strsplit(RNAconf[UI == input$RNAa3sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("RNAa3sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$RNAa3sub1non, { 
    sub = strsplit(RNAconf[UI == input$RNAa3sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "RNAa3sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$RNAa3sub1all, { 
    sub = strsplit(RNAconf[UI == input$RNAa3sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "RNAa3sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$RNAa3oup1 <- renderPlot({ 
    scDRgene(RNAconf, RNAmeta, input$RNAa3drX, input$RNAa3drY, input$RNAa3inp1,  
             input$RNAa3sub1, input$RNAa3sub2, 
             "RNAgexpr.h5", RNAgene, 
             input$RNAa3siz, input$RNAa3col1, input$RNAa3ord1, 
             input$RNAa3fsz, input$RNAa3asp, input$RNAa3txt) 
  }) 
  output$RNAa3oup1.ui <- renderUI({ 
    plotOutput("RNAa3oup1", height = pList[input$RNAa3psz]) 
  }) 
  output$RNAa3oup1.pdf <- downloadHandler( 
    filename = function() { paste0("RNA",input$RNAa3drX,"_",input$RNAa3drY,"_",  
                                   input$RNAa3inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$RNAa3oup1.h, width = input$RNAa3oup1.w, useDingbats = FALSE, 
      plot = scDRgene(RNAconf, RNAmeta, input$RNAa3drX, input$RNAa3drY, input$RNAa3inp1,  
                      input$RNAa3sub1, input$RNAa3sub2, 
                      "RNAgexpr.h5", RNAgene, 
                      input$RNAa3siz, input$RNAa3col1, input$RNAa3ord1, 
                      input$RNAa3fsz, input$RNAa3asp, input$RNAa3txt) ) 
  }) 
  output$RNAa3oup1.png <- downloadHandler( 
    filename = function() { paste0("RNA",input$RNAa3drX,"_",input$RNAa3drY,"_",  
                                   input$RNAa3inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$RNAa3oup1.h, width = input$RNAa3oup1.w, 
      plot = scDRgene(RNAconf, RNAmeta, input$RNAa3drX, input$RNAa3drY, input$RNAa3inp1,  
                      input$RNAa3sub1, input$RNAa3sub2, 
                      "RNAgexpr.h5", RNAgene, 
                      input$RNAa3siz, input$RNAa3col1, input$RNAa3ord1, 
                      input$RNAa3fsz, input$RNAa3asp, input$RNAa3txt) ) 
  }) 
   
  output$RNAa3oup2 <- renderPlot({ 
    scDRgene(RNAconf, RNAmeta, input$RNAa3drX, input$RNAa3drY, input$RNAa3inp2,  
             input$RNAa3sub1, input$RNAa3sub2, 
             "RNAgexpr.h5", RNAgene, 
             input$RNAa3siz, input$RNAa3col2, input$RNAa3ord2, 
             input$RNAa3fsz, input$RNAa3asp, input$RNAa3txt) 
  }) 
  output$RNAa3oup2.ui <- renderUI({ 
    plotOutput("RNAa3oup2", height = pList[input$RNAa3psz]) 
  }) 
  output$RNAa3oup2.pdf <- downloadHandler( 
    filename = function() { paste0("RNA",input$RNAa3drX,"_",input$RNAa3drY,"_",  
                                   input$RNAa3inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$RNAa3oup2.h, width = input$RNAa3oup2.w, useDingbats = FALSE, 
      plot = scDRgene(RNAconf, RNAmeta, input$RNAa3drX, input$RNAa3drY, input$RNAa3inp2,  
                      input$RNAa3sub1, input$RNAa3sub2, 
                      "RNAgexpr.h5", RNAgene, 
                      input$RNAa3siz, input$RNAa3col2, input$RNAa3ord2, 
                      input$RNAa3fsz, input$RNAa3asp, input$RNAa3txt) ) 
  }) 
  output$RNAa3oup2.png <- downloadHandler( 
    filename = function() { paste0("RNA",input$RNAa3drX,"_",input$RNAa3drY,"_",  
                                   input$RNAa3inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$RNAa3oup2.h, width = input$RNAa3oup2.w, 
      plot = scDRgene(RNAconf, RNAmeta, input$RNAa3drX, input$RNAa3drY, input$RNAa3inp2,  
                      input$RNAa3sub1, input$RNAa3sub2, 
                      "RNAgexpr.h5", RNAgene, 
                      input$RNAa3siz, input$RNAa3col2, input$RNAa3ord2, 
                      input$RNAa3fsz, input$RNAa3asp, input$RNAa3txt) ) 
  }) 
     
   
  ### Plots for tab b2 
  output$RNAb2sub1.ui <- renderUI({ 
    sub = strsplit(RNAconf[UI == input$RNAb2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("RNAb2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$RNAb2sub1non, { 
    sub = strsplit(RNAconf[UI == input$RNAb2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "RNAb2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$RNAb2sub1all, { 
    sub = strsplit(RNAconf[UI == input$RNAb2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "RNAb2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$RNAb2oup1 <- renderPlot({ 
    scDRcoex(RNAconf, RNAmeta, input$RNAb2drX, input$RNAb2drY,   
             input$RNAb2inp1, input$RNAb2inp2, input$RNAb2sub1, input$RNAb2sub2, 
             "RNAgexpr.h5", RNAgene, 
             input$RNAb2siz, input$RNAb2col1, input$RNAb2ord1, 
             input$RNAb2fsz, input$RNAb2asp, input$RNAb2txt) 
  }) 
  output$RNAb2oup1.ui <- renderUI({ 
    plotOutput("RNAb2oup1", height = pList2[input$RNAb2psz]) 
  }) 
  output$RNAb2oup1.pdf <- downloadHandler( 
    filename = function() { paste0("RNA",input$RNAb2drX,"_",input$RNAb2drY,"_",  
                                    input$RNAb2inp1,"_",input$RNAb2inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$RNAb2oup1.h, width = input$RNAb2oup1.w, useDingbats = FALSE, 
      plot = scDRcoex(RNAconf, RNAmeta, input$RNAb2drX, input$RNAb2drY,  
                      input$RNAb2inp1, input$RNAb2inp2, input$RNAb2sub1, input$RNAb2sub2, 
                      "RNAgexpr.h5", RNAgene, 
                      input$RNAb2siz, input$RNAb2col1, input$RNAb2ord1, 
                      input$RNAb2fsz, input$RNAb2asp, input$RNAb2txt) ) 
  }) 
  output$RNAb2oup1.png <- downloadHandler( 
    filename = function() { paste0("RNA",input$RNAb2drX,"_",input$RNAb2drY,"_",  
                                    input$RNAb2inp1,"_",input$RNAb2inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$RNAb2oup1.h, width = input$RNAb2oup1.w, 
      plot = scDRcoex(RNAconf, RNAmeta, input$RNAb2drX, input$RNAb2drY,  
                      input$RNAb2inp1, input$RNAb2inp2, input$RNAb2sub1, input$RNAb2sub2, 
                      "RNAgexpr.h5", RNAgene, 
                      input$RNAb2siz, input$RNAb2col1, input$RNAb2ord1, 
                      input$RNAb2fsz, input$RNAb2asp, input$RNAb2txt) ) 
  }) 
  output$RNAb2oup2 <- renderPlot({ 
    scDRcoexLeg(input$RNAb2inp1, input$RNAb2inp2, input$RNAb2col1, input$RNAb2fsz) 
  }) 
  output$RNAb2oup2.ui <- renderUI({ 
    plotOutput("RNAb2oup2", height = "300px") 
  }) 
  output$RNAb2oup2.pdf <- downloadHandler( 
    filename = function() { paste0("RNA",input$RNAb2drX,"_",input$RNAb2drY,"_",  
                                    input$RNAb2inp1,"_",input$RNAb2inp2,"_leg.pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = 3, width = 4, useDingbats = FALSE, 
      plot = scDRcoexLeg(input$RNAb2inp1, input$RNAb2inp2, input$RNAb2col1, input$RNAb2fsz) ) 
  }) 
  output$RNAb2oup2.png <- downloadHandler( 
    filename = function() { paste0("RNA",input$RNAb2drX,"_",input$RNAb2drY,"_",  
                                    input$RNAb2inp1,"_",input$RNAb2inp2,"_leg.png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = 3, width = 4, 
      plot = scDRcoexLeg(input$RNAb2inp1, input$RNAb2inp2, input$RNAb2col1, input$RNAb2fsz) ) 
  }) 
  output$RNAb2.dt <- renderDataTable({ 
    ggData = scDRcoexNum(RNAconf, RNAmeta, input$RNAb2inp1, input$RNAb2inp2, 
                         input$RNAb2sub1, input$RNAb2sub2, "RNAgexpr.h5", RNAgene) 
    datatable(ggData, rownames = FALSE, extensions = "Buttons", 
              options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>% 
      formatRound(columns = c("percent"), digits = 2) 
  }) 
     
   
  ### Plots for tab c1 
  output$RNAc1sub1.ui <- renderUI({ 
    sub = strsplit(RNAconf[UI == input$RNAc1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("RNAc1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$RNAc1sub1non, { 
    sub = strsplit(RNAconf[UI == input$RNAc1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "RNAc1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$RNAc1sub1all, { 
    sub = strsplit(RNAconf[UI == input$RNAc1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "RNAc1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$RNAc1oup <- renderPlot({ 
    scVioBox(RNAconf, RNAmeta, input$RNAc1inp1, input$RNAc1inp2, 
             input$RNAc1sub1, input$RNAc1sub2, 
             "RNAgexpr.h5", RNAgene, input$RNAc1typ, input$RNAc1pts, 
             input$RNAc1siz, input$RNAc1fsz) 
  }) 
  output$RNAc1oup.ui <- renderUI({ 
    plotOutput("RNAc1oup", height = pList2[input$RNAc1psz]) 
  }) 
  output$RNAc1oup.pdf <- downloadHandler( 
    filename = function() { paste0("RNA",input$RNAc1typ,"_",input$RNAc1inp1,"_",  
                                   input$RNAc1inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$RNAc1oup.h, width = input$RNAc1oup.w, useDingbats = FALSE, 
      plot = scVioBox(RNAconf, RNAmeta, input$RNAc1inp1, input$RNAc1inp2, 
                      input$RNAc1sub1, input$RNAc1sub2, 
                      "RNAgexpr.h5", RNAgene, input$RNAc1typ, input$RNAc1pts, 
                      input$RNAc1siz, input$RNAc1fsz) ) 
  }) 
  output$RNAc1oup.png <- downloadHandler( 
    filename = function() { paste0("RNA",input$RNAc1typ,"_",input$RNAc1inp1,"_",  
                                   input$RNAc1inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$RNAc1oup.h, width = input$RNAc1oup.w, 
      plot = scVioBox(RNAconf, RNAmeta, input$RNAc1inp1, input$RNAc1inp2, 
                      input$RNAc1sub1, input$RNAc1sub2, 
                      "RNAgexpr.h5", RNAgene, input$RNAc1typ, input$RNAc1pts, 
                      input$RNAc1siz, input$RNAc1fsz) ) 
  }) 
     
   
### Plots for tab c2 
  output$RNAc2sub1.ui <- renderUI({ 
    sub = strsplit(RNAconf[UI == input$RNAc2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("RNAc2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$RNAc2sub1non, { 
    sub = strsplit(RNAconf[UI == input$RNAc2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "RNAc2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$RNAc2sub1all, { 
    sub = strsplit(RNAconf[UI == input$RNAc2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "RNAc2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
output$RNAc2oup <- renderPlot({ 
  scProp(RNAconf, RNAmeta, input$RNAc2inp1, input$RNAc2inp2,  
         input$RNAc2sub1, input$RNAc2sub2, 
         input$RNAc2typ, input$RNAc2flp, input$RNAc2fsz) 
}) 
output$RNAc2oup.ui <- renderUI({ 
  plotOutput("RNAc2oup", height = pList2[input$RNAc2psz]) 
}) 
output$RNAc2oup.pdf <- downloadHandler( 
  filename = function() { paste0("RNA",input$RNAc2typ,"_",input$RNAc2inp1,"_",  
                                 input$RNAc2inp2,".pdf") }, 
  content = function(file) { ggsave( 
    file, device = "pdf", height = input$RNAc2oup.h, width = input$RNAc2oup.w, useDingbats = FALSE, 
    plot = scProp(RNAconf, RNAmeta, input$RNAc2inp1, input$RNAc2inp2,  
                  input$RNAc2sub1, input$RNAc2sub2, 
                  input$RNAc2typ, input$RNAc2flp, input$RNAc2fsz) ) 
  }) 
output$RNAc2oup.png <- downloadHandler( 
  filename = function() { paste0("RNA",input$RNAc2typ,"_",input$RNAc2inp1,"_",  
                                 input$RNAc2inp2,".png") }, 
  content = function(file) { ggsave( 
    file, device = "png", height = input$RNAc2oup.h, width = input$RNAc2oup.w, 
    plot = scProp(RNAconf, RNAmeta, input$RNAc2inp1, input$RNAc2inp2,  
                  input$RNAc2sub1, input$RNAc2sub2, 
                  input$RNAc2typ, input$RNAc2flp, input$RNAc2fsz) ) 
  }) 
     
   
  ### Plots for tab d1 
  output$RNAd1sub1.ui <- renderUI({ 
    sub = strsplit(RNAconf[UI == input$RNAd1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("RNAd1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$RNAd1sub1non, { 
    sub = strsplit(RNAconf[UI == input$RNAd1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "RNAd1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$RNAd1sub1all, { 
    sub = strsplit(RNAconf[UI == input$RNAd1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "RNAd1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$RNAd1oupTxt <- renderUI({ 
    geneList = scGeneList(input$RNAd1inp, RNAgene) 
    if(nrow(geneList) > 50){ 
      HTML("More than 50 input genes! Please reduce the gene list!") 
    } else { 
      oup = paste0(nrow(geneList[present == TRUE]), " genes OK and will be plotted") 
      if(nrow(geneList[present == FALSE]) > 0){ 
        oup = paste0(oup, "<br/>", 
                     nrow(geneList[present == FALSE]), " genes not found (", 
                     paste0(geneList[present == FALSE]$gene, collapse = ", "), ")") 
      } 
      HTML(oup) 
    } 
  }) 
  output$RNAd1oup <- renderPlot({ 
    scBubbHeat(RNAconf, RNAmeta, input$RNAd1inp, input$RNAd1grp, input$RNAd1plt, 
               input$RNAd1sub1, input$RNAd1sub2, "RNAgexpr.h5", RNAgene, 
               input$RNAd1scl, input$RNAd1row, input$RNAd1col, 
               input$RNAd1cols, input$RNAd1fsz) 
  }) 
  output$RNAd1oup.ui <- renderUI({ 
    plotOutput("RNAd1oup", height = pList3[input$RNAd1psz]) 
  }) 
  output$RNAd1oup.pdf <- downloadHandler( 
    filename = function() { paste0("RNA",input$RNAd1plt,"_",input$RNAd1grp,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$RNAd1oup.h, width = input$RNAd1oup.w, 
      plot = scBubbHeat(RNAconf, RNAmeta, input$RNAd1inp, input$RNAd1grp, input$RNAd1plt, 
                        input$RNAd1sub1, input$RNAd1sub2, "RNAgexpr.h5", RNAgene, 
                        input$RNAd1scl, input$RNAd1row, input$RNAd1col, 
                        input$RNAd1cols, input$RNAd1fsz, save = TRUE) ) 
  }) 
  output$RNAd1oup.png <- downloadHandler( 
    filename = function() { paste0("RNA",input$RNAd1plt,"_",input$RNAd1grp,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$RNAd1oup.h, width = input$RNAd1oup.w, 
      plot = scBubbHeat(RNAconf, RNAmeta, input$RNAd1inp, input$RNAd1grp, input$RNAd1plt, 
                        input$RNAd1sub1, input$RNAd1sub2, "RNAgexpr.h5", RNAgene, 
                        input$RNAd1scl, input$RNAd1row, input$RNAd1col, 
                        input$RNAd1cols, input$RNAd1fsz, save = TRUE) ) 
  }) 
   
   
   optCrt="{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }" 
  updateSelectizeInput(session, "ATACa1inp2", choices = names(ATACgene), server = TRUE, 
                       selected = ATACdef$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "ATACa3inp1", choices = names(ATACgene), server = TRUE, 
                       selected = ATACdef$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "ATACa3inp2", choices = names(ATACgene), server = TRUE, 
                       selected = ATACdef$gene2, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "ATACb2inp1", choices = names(ATACgene), server = TRUE, 
                       selected = ATACdef$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "ATACb2inp2", choices = names(ATACgene), server = TRUE, 
                       selected = ATACdef$gene2, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "ATACc1inp2", server = TRUE, 
                       choices = c(ATACconf[is.na(fID)]$UI,names(ATACgene)), 
                       selected = ATACconf[is.na(fID)]$UI[1], options = list( 
                         maxOptions = length(ATACconf[is.na(fID)]$UI) + 3, 
                         create = TRUE, persist = TRUE, render = I(optCrt))) 
 
 
  ### Plots for tab a1 
  output$ATACa1sub1.ui <- renderUI({ 
    sub = strsplit(ATACconf[UI == input$ATACa1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("ATACa1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$ATACa1sub1non, { 
    sub = strsplit(ATACconf[UI == input$ATACa1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "ATACa1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$ATACa1sub1all, { 
    sub = strsplit(ATACconf[UI == input$ATACa1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "ATACa1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$ATACa1oup1 <- renderPlot({ 
    scDRcell(ATACconf, ATACmeta, input$ATACa1drX, input$ATACa1drY, input$ATACa1inp1,  
             input$ATACa1sub1, input$ATACa1sub2, 
             input$ATACa1siz, input$ATACa1col1, input$ATACa1ord1, 
             input$ATACa1fsz, input$ATACa1asp, input$ATACa1txt, input$ATACa1lab1) 
  }) 
  output$ATACa1oup1.ui <- renderUI({ 
    plotOutput("ATACa1oup1", height = pList[input$ATACa1psz]) 
  }) 
  output$ATACa1oup1.pdf <- downloadHandler( 
    filename = function() { paste0("ATAC",input$ATACa1drX,"_",input$ATACa1drY,"_",  
                                   input$ATACa1inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$ATACa1oup1.h, width = input$ATACa1oup1.w, useDingbats = FALSE, 
      plot = scDRcell(ATACconf, ATACmeta, input$ATACa1drX, input$ATACa1drY, input$ATACa1inp1,   
                      input$ATACa1sub1, input$ATACa1sub2, 
                      input$ATACa1siz, input$ATACa1col1, input$ATACa1ord1,  
                      input$ATACa1fsz, input$ATACa1asp, input$ATACa1txt, input$ATACa1lab1) ) 
  }) 
  output$ATACa1oup1.png <- downloadHandler( 
    filename = function() { paste0("ATAC",input$ATACa1drX,"_",input$ATACa1drY,"_",  
                                   input$ATACa1inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$ATACa1oup1.h, width = input$ATACa1oup1.w, 
      plot = scDRcell(ATACconf, ATACmeta, input$ATACa1drX, input$ATACa1drY, input$ATACa1inp1,   
                      input$ATACa1sub1, input$ATACa1sub2, 
                      input$ATACa1siz, input$ATACa1col1, input$ATACa1ord1,  
                      input$ATACa1fsz, input$ATACa1asp, input$ATACa1txt, input$ATACa1lab1) ) 
  }) 
  output$ATACa1.dt <- renderDataTable({ 
    ggData = scDRnum(ATACconf, ATACmeta, input$ATACa1inp1, input$ATACa1inp2, 
                     input$ATACa1sub1, input$ATACa1sub2, 
                     "ATACgexpr.h5", ATACgene, input$ATACa1splt) 
    datatable(ggData, rownames = FALSE, extensions = "Buttons", 
              options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>% 
      formatRound(columns = c("pctExpress"), digits = 2) 
  }) 
   
  output$ATACa1oup2 <- renderPlot({ 
    scDRgene(ATACconf, ATACmeta, input$ATACa1drX, input$ATACa1drY, input$ATACa1inp2,  
             input$ATACa1sub1, input$ATACa1sub2, 
             "ATACgexpr.h5", ATACgene, 
             input$ATACa1siz, input$ATACa1col2, input$ATACa1ord2, 
             input$ATACa1fsz, input$ATACa1asp, input$ATACa1txt) 
  }) 
  output$ATACa1oup2.ui <- renderUI({ 
    plotOutput("ATACa1oup2", height = pList[input$ATACa1psz]) 
  }) 
  output$ATACa1oup2.pdf <- downloadHandler( 
    filename = function() { paste0("ATAC",input$ATACa1drX,"_",input$ATACa1drY,"_",  
                                   input$ATACa1inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$ATACa1oup2.h, width = input$ATACa1oup2.w, useDingbats = FALSE, 
      plot = scDRgene(ATACconf, ATACmeta, input$ATACa1drX, input$ATACa1drY, input$ATACa1inp2,  
                      input$ATACa1sub1, input$ATACa1sub2, 
                      "ATACgexpr.h5", ATACgene, 
                      input$ATACa1siz, input$ATACa1col2, input$ATACa1ord2, 
                      input$ATACa1fsz, input$ATACa1asp, input$ATACa1txt) ) 
  }) 
  output$ATACa1oup2.png <- downloadHandler( 
    filename = function() { paste0("ATAC",input$ATACa1drX,"_",input$ATACa1drY,"_",  
                                   input$ATACa1inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$ATACa1oup2.h, width = input$ATACa1oup2.w, 
      plot = scDRgene(ATACconf, ATACmeta, input$ATACa1drX, input$ATACa1drY, input$ATACa1inp2,  
                      input$ATACa1sub1, input$ATACa1sub2, 
                      "ATACgexpr.h5", ATACgene, 
                      input$ATACa1siz, input$ATACa1col2, input$ATACa1ord2, 
                      input$ATACa1fsz, input$ATACa1asp, input$ATACa1txt) ) 
  }) 
   
   
  ### Plots for tab a2 
  output$ATACa2sub1.ui <- renderUI({ 
    sub = strsplit(ATACconf[UI == input$ATACa2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("ATACa2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$ATACa2sub1non, { 
    sub = strsplit(ATACconf[UI == input$ATACa2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "ATACa2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$ATACa2sub1all, { 
    sub = strsplit(ATACconf[UI == input$ATACa2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "ATACa2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$ATACa2oup1 <- renderPlot({ 
    scDRcell(ATACconf, ATACmeta, input$ATACa2drX, input$ATACa2drY, input$ATACa2inp1,  
             input$ATACa2sub1, input$ATACa2sub2, 
             input$ATACa2siz, input$ATACa2col1, input$ATACa2ord1, 
             input$ATACa2fsz, input$ATACa2asp, input$ATACa2txt, input$ATACa2lab1) 
  }) 
  output$ATACa2oup1.ui <- renderUI({ 
    plotOutput("ATACa2oup1", height = pList[input$ATACa2psz]) 
  }) 
  output$ATACa2oup1.pdf <- downloadHandler( 
    filename = function() { paste0("ATAC",input$ATACa2drX,"_",input$ATACa2drY,"_",  
                                   input$ATACa2inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$ATACa2oup1.h, width = input$ATACa2oup1.w, useDingbats = FALSE, 
      plot = scDRcell(ATACconf, ATACmeta, input$ATACa2drX, input$ATACa2drY, input$ATACa2inp1,   
                      input$ATACa2sub1, input$ATACa2sub2, 
                      input$ATACa2siz, input$ATACa2col1, input$ATACa2ord1,  
                      input$ATACa2fsz, input$ATACa2asp, input$ATACa2txt, input$ATACa2lab1) ) 
  }) 
  output$ATACa2oup1.png <- downloadHandler( 
    filename = function() { paste0("ATAC",input$ATACa2drX,"_",input$ATACa2drY,"_",  
                                   input$ATACa2inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$ATACa2oup1.h, width = input$ATACa2oup1.w, 
      plot = scDRcell(ATACconf, ATACmeta, input$ATACa2drX, input$ATACa2drY, input$ATACa2inp1,   
                      input$ATACa2sub1, input$ATACa2sub2, 
                      input$ATACa2siz, input$ATACa2col1, input$ATACa2ord1,  
                      input$ATACa2fsz, input$ATACa2asp, input$ATACa2txt, input$ATACa2lab1) ) 
  }) 
   
  output$ATACa2oup2 <- renderPlot({ 
    scDRcell(ATACconf, ATACmeta, input$ATACa2drX, input$ATACa2drY, input$ATACa2inp2,  
             input$ATACa2sub1, input$ATACa2sub2, 
             input$ATACa2siz, input$ATACa2col2, input$ATACa2ord2, 
             input$ATACa2fsz, input$ATACa2asp, input$ATACa2txt, input$ATACa2lab2) 
  }) 
  output$ATACa2oup2.ui <- renderUI({ 
    plotOutput("ATACa2oup2", height = pList[input$ATACa2psz]) 
  }) 
  output$ATACa2oup2.pdf <- downloadHandler( 
    filename = function() { paste0("ATAC",input$ATACa2drX,"_",input$ATACa2drY,"_",  
                                   input$ATACa2inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$ATACa2oup2.h, width = input$ATACa2oup2.w, useDingbats = FALSE, 
      plot = scDRcell(ATACconf, ATACmeta, input$ATACa2drX, input$ATACa2drY, input$ATACa2inp2,   
                      input$ATACa2sub1, input$ATACa2sub2, 
                      input$ATACa2siz, input$ATACa2col2, input$ATACa2ord2,  
                      input$ATACa2fsz, input$ATACa2asp, input$ATACa2txt, input$ATACa2lab2) ) 
  }) 
  output$ATACa2oup2.png <- downloadHandler( 
    filename = function() { paste0("ATAC",input$ATACa2drX,"_",input$ATACa2drY,"_",  
                                   input$ATACa2inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$ATACa2oup2.h, width = input$ATACa2oup2.w, 
      plot = scDRcell(ATACconf, ATACmeta, input$ATACa2drX, input$ATACa2drY, input$ATACa2inp2,   
                      input$ATACa2sub1, input$ATACa2sub2, 
                      input$ATACa2siz, input$ATACa2col2, input$ATACa2ord2,  
                      input$ATACa2fsz, input$ATACa2asp, input$ATACa2txt, input$ATACa2lab2) ) 
  }) 
   
   
  ### Plots for tab a3 
  output$ATACa3sub1.ui <- renderUI({ 
    sub = strsplit(ATACconf[UI == input$ATACa3sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("ATACa3sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$ATACa3sub1non, { 
    sub = strsplit(ATACconf[UI == input$ATACa3sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "ATACa3sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$ATACa3sub1all, { 
    sub = strsplit(ATACconf[UI == input$ATACa3sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "ATACa3sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$ATACa3oup1 <- renderPlot({ 
    scDRgene(ATACconf, ATACmeta, input$ATACa3drX, input$ATACa3drY, input$ATACa3inp1,  
             input$ATACa3sub1, input$ATACa3sub2, 
             "ATACgexpr.h5", ATACgene, 
             input$ATACa3siz, input$ATACa3col1, input$ATACa3ord1, 
             input$ATACa3fsz, input$ATACa3asp, input$ATACa3txt) 
  }) 
  output$ATACa3oup1.ui <- renderUI({ 
    plotOutput("ATACa3oup1", height = pList[input$ATACa3psz]) 
  }) 
  output$ATACa3oup1.pdf <- downloadHandler( 
    filename = function() { paste0("ATAC",input$ATACa3drX,"_",input$ATACa3drY,"_",  
                                   input$ATACa3inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$ATACa3oup1.h, width = input$ATACa3oup1.w, useDingbats = FALSE, 
      plot = scDRgene(ATACconf, ATACmeta, input$ATACa3drX, input$ATACa3drY, input$ATACa3inp1,  
                      input$ATACa3sub1, input$ATACa3sub2, 
                      "ATACgexpr.h5", ATACgene, 
                      input$ATACa3siz, input$ATACa3col1, input$ATACa3ord1, 
                      input$ATACa3fsz, input$ATACa3asp, input$ATACa3txt) ) 
  }) 
  output$ATACa3oup1.png <- downloadHandler( 
    filename = function() { paste0("ATAC",input$ATACa3drX,"_",input$ATACa3drY,"_",  
                                   input$ATACa3inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$ATACa3oup1.h, width = input$ATACa3oup1.w, 
      plot = scDRgene(ATACconf, ATACmeta, input$ATACa3drX, input$ATACa3drY, input$ATACa3inp1,  
                      input$ATACa3sub1, input$ATACa3sub2, 
                      "ATACgexpr.h5", ATACgene, 
                      input$ATACa3siz, input$ATACa3col1, input$ATACa3ord1, 
                      input$ATACa3fsz, input$ATACa3asp, input$ATACa3txt) ) 
  }) 
   
  output$ATACa3oup2 <- renderPlot({ 
    scDRgene(ATACconf, ATACmeta, input$ATACa3drX, input$ATACa3drY, input$ATACa3inp2,  
             input$ATACa3sub1, input$ATACa3sub2, 
             "ATACgexpr.h5", ATACgene, 
             input$ATACa3siz, input$ATACa3col2, input$ATACa3ord2, 
             input$ATACa3fsz, input$ATACa3asp, input$ATACa3txt) 
  }) 
  output$ATACa3oup2.ui <- renderUI({ 
    plotOutput("ATACa3oup2", height = pList[input$ATACa3psz]) 
  }) 
  output$ATACa3oup2.pdf <- downloadHandler( 
    filename = function() { paste0("ATAC",input$ATACa3drX,"_",input$ATACa3drY,"_",  
                                   input$ATACa3inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$ATACa3oup2.h, width = input$ATACa3oup2.w, useDingbats = FALSE, 
      plot = scDRgene(ATACconf, ATACmeta, input$ATACa3drX, input$ATACa3drY, input$ATACa3inp2,  
                      input$ATACa3sub1, input$ATACa3sub2, 
                      "ATACgexpr.h5", ATACgene, 
                      input$ATACa3siz, input$ATACa3col2, input$ATACa3ord2, 
                      input$ATACa3fsz, input$ATACa3asp, input$ATACa3txt) ) 
  }) 
  output$ATACa3oup2.png <- downloadHandler( 
    filename = function() { paste0("ATAC",input$ATACa3drX,"_",input$ATACa3drY,"_",  
                                   input$ATACa3inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$ATACa3oup2.h, width = input$ATACa3oup2.w, 
      plot = scDRgene(ATACconf, ATACmeta, input$ATACa3drX, input$ATACa3drY, input$ATACa3inp2,  
                      input$ATACa3sub1, input$ATACa3sub2, 
                      "ATACgexpr.h5", ATACgene, 
                      input$ATACa3siz, input$ATACa3col2, input$ATACa3ord2, 
                      input$ATACa3fsz, input$ATACa3asp, input$ATACa3txt) ) 
  }) 
     
   
  ### Plots for tab b2 
  output$ATACb2sub1.ui <- renderUI({ 
    sub = strsplit(ATACconf[UI == input$ATACb2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("ATACb2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$ATACb2sub1non, { 
    sub = strsplit(ATACconf[UI == input$ATACb2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "ATACb2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$ATACb2sub1all, { 
    sub = strsplit(ATACconf[UI == input$ATACb2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "ATACb2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$ATACb2oup1 <- renderPlot({ 
    scDRcoex(ATACconf, ATACmeta, input$ATACb2drX, input$ATACb2drY,   
             input$ATACb2inp1, input$ATACb2inp2, input$ATACb2sub1, input$ATACb2sub2, 
             "ATACgexpr.h5", ATACgene, 
             input$ATACb2siz, input$ATACb2col1, input$ATACb2ord1, 
             input$ATACb2fsz, input$ATACb2asp, input$ATACb2txt) 
  }) 
  output$ATACb2oup1.ui <- renderUI({ 
    plotOutput("ATACb2oup1", height = pList2[input$ATACb2psz]) 
  }) 
  output$ATACb2oup1.pdf <- downloadHandler( 
    filename = function() { paste0("ATAC",input$ATACb2drX,"_",input$ATACb2drY,"_",  
                                    input$ATACb2inp1,"_",input$ATACb2inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$ATACb2oup1.h, width = input$ATACb2oup1.w, useDingbats = FALSE, 
      plot = scDRcoex(ATACconf, ATACmeta, input$ATACb2drX, input$ATACb2drY,  
                      input$ATACb2inp1, input$ATACb2inp2, input$ATACb2sub1, input$ATACb2sub2, 
                      "ATACgexpr.h5", ATACgene, 
                      input$ATACb2siz, input$ATACb2col1, input$ATACb2ord1, 
                      input$ATACb2fsz, input$ATACb2asp, input$ATACb2txt) ) 
  }) 
  output$ATACb2oup1.png <- downloadHandler( 
    filename = function() { paste0("ATAC",input$ATACb2drX,"_",input$ATACb2drY,"_",  
                                    input$ATACb2inp1,"_",input$ATACb2inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$ATACb2oup1.h, width = input$ATACb2oup1.w, 
      plot = scDRcoex(ATACconf, ATACmeta, input$ATACb2drX, input$ATACb2drY,  
                      input$ATACb2inp1, input$ATACb2inp2, input$ATACb2sub1, input$ATACb2sub2, 
                      "ATACgexpr.h5", ATACgene, 
                      input$ATACb2siz, input$ATACb2col1, input$ATACb2ord1, 
                      input$ATACb2fsz, input$ATACb2asp, input$ATACb2txt) ) 
  }) 
  output$ATACb2oup2 <- renderPlot({ 
    scDRcoexLeg(input$ATACb2inp1, input$ATACb2inp2, input$ATACb2col1, input$ATACb2fsz) 
  }) 
  output$ATACb2oup2.ui <- renderUI({ 
    plotOutput("ATACb2oup2", height = "300px") 
  }) 
  output$ATACb2oup2.pdf <- downloadHandler( 
    filename = function() { paste0("ATAC",input$ATACb2drX,"_",input$ATACb2drY,"_",  
                                    input$ATACb2inp1,"_",input$ATACb2inp2,"_leg.pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = 3, width = 4, useDingbats = FALSE, 
      plot = scDRcoexLeg(input$ATACb2inp1, input$ATACb2inp2, input$ATACb2col1, input$ATACb2fsz) ) 
  }) 
  output$ATACb2oup2.png <- downloadHandler( 
    filename = function() { paste0("ATAC",input$ATACb2drX,"_",input$ATACb2drY,"_",  
                                    input$ATACb2inp1,"_",input$ATACb2inp2,"_leg.png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = 3, width = 4, 
      plot = scDRcoexLeg(input$ATACb2inp1, input$ATACb2inp2, input$ATACb2col1, input$ATACb2fsz) ) 
  }) 
  output$ATACb2.dt <- renderDataTable({ 
    ggData = scDRcoexNum(ATACconf, ATACmeta, input$ATACb2inp1, input$ATACb2inp2, 
                         input$ATACb2sub1, input$ATACb2sub2, "ATACgexpr.h5", ATACgene) 
    datatable(ggData, rownames = FALSE, extensions = "Buttons", 
              options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>% 
      formatRound(columns = c("percent"), digits = 2) 
  }) 
     
   
  ### Plots for tab c1 
  output$ATACc1sub1.ui <- renderUI({ 
    sub = strsplit(ATACconf[UI == input$ATACc1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("ATACc1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$ATACc1sub1non, { 
    sub = strsplit(ATACconf[UI == input$ATACc1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "ATACc1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$ATACc1sub1all, { 
    sub = strsplit(ATACconf[UI == input$ATACc1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "ATACc1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$ATACc1oup <- renderPlot({ 
    scVioBox(ATACconf, ATACmeta, input$ATACc1inp1, input$ATACc1inp2, 
             input$ATACc1sub1, input$ATACc1sub2, 
             "ATACgexpr.h5", ATACgene, input$ATACc1typ, input$ATACc1pts, 
             input$ATACc1siz, input$ATACc1fsz) 
  }) 
  output$ATACc1oup.ui <- renderUI({ 
    plotOutput("ATACc1oup", height = pList2[input$ATACc1psz]) 
  }) 
  output$ATACc1oup.pdf <- downloadHandler( 
    filename = function() { paste0("ATAC",input$ATACc1typ,"_",input$ATACc1inp1,"_",  
                                   input$ATACc1inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$ATACc1oup.h, width = input$ATACc1oup.w, useDingbats = FALSE, 
      plot = scVioBox(ATACconf, ATACmeta, input$ATACc1inp1, input$ATACc1inp2, 
                      input$ATACc1sub1, input$ATACc1sub2, 
                      "ATACgexpr.h5", ATACgene, input$ATACc1typ, input$ATACc1pts, 
                      input$ATACc1siz, input$ATACc1fsz) ) 
  }) 
  output$ATACc1oup.png <- downloadHandler( 
    filename = function() { paste0("ATAC",input$ATACc1typ,"_",input$ATACc1inp1,"_",  
                                   input$ATACc1inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$ATACc1oup.h, width = input$ATACc1oup.w, 
      plot = scVioBox(ATACconf, ATACmeta, input$ATACc1inp1, input$ATACc1inp2, 
                      input$ATACc1sub1, input$ATACc1sub2, 
                      "ATACgexpr.h5", ATACgene, input$ATACc1typ, input$ATACc1pts, 
                      input$ATACc1siz, input$ATACc1fsz) ) 
  }) 
     
   
### Plots for tab c2 
  output$ATACc2sub1.ui <- renderUI({ 
    sub = strsplit(ATACconf[UI == input$ATACc2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("ATACc2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$ATACc2sub1non, { 
    sub = strsplit(ATACconf[UI == input$ATACc2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "ATACc2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$ATACc2sub1all, { 
    sub = strsplit(ATACconf[UI == input$ATACc2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "ATACc2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
output$ATACc2oup <- renderPlot({ 
  scProp(ATACconf, ATACmeta, input$ATACc2inp1, input$ATACc2inp2,  
         input$ATACc2sub1, input$ATACc2sub2, 
         input$ATACc2typ, input$ATACc2flp, input$ATACc2fsz) 
}) 
output$ATACc2oup.ui <- renderUI({ 
  plotOutput("ATACc2oup", height = pList2[input$ATACc2psz]) 
}) 
output$ATACc2oup.pdf <- downloadHandler( 
  filename = function() { paste0("ATAC",input$ATACc2typ,"_",input$ATACc2inp1,"_",  
                                 input$ATACc2inp2,".pdf") }, 
  content = function(file) { ggsave( 
    file, device = "pdf", height = input$ATACc2oup.h, width = input$ATACc2oup.w, useDingbats = FALSE, 
    plot = scProp(ATACconf, ATACmeta, input$ATACc2inp1, input$ATACc2inp2,  
                  input$ATACc2sub1, input$ATACc2sub2, 
                  input$ATACc2typ, input$ATACc2flp, input$ATACc2fsz) ) 
  }) 
output$ATACc2oup.png <- downloadHandler( 
  filename = function() { paste0("ATAC",input$ATACc2typ,"_",input$ATACc2inp1,"_",  
                                 input$ATACc2inp2,".png") }, 
  content = function(file) { ggsave( 
    file, device = "png", height = input$ATACc2oup.h, width = input$ATACc2oup.w, 
    plot = scProp(ATACconf, ATACmeta, input$ATACc2inp1, input$ATACc2inp2,  
                  input$ATACc2sub1, input$ATACc2sub2, 
                  input$ATACc2typ, input$ATACc2flp, input$ATACc2fsz) ) 
  }) 
     
   
  ### Plots for tab d1 
  output$ATACd1sub1.ui <- renderUI({ 
    sub = strsplit(ATACconf[UI == input$ATACd1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("ATACd1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$ATACd1sub1non, { 
    sub = strsplit(ATACconf[UI == input$ATACd1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "ATACd1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$ATACd1sub1all, { 
    sub = strsplit(ATACconf[UI == input$ATACd1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "ATACd1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$ATACd1oupTxt <- renderUI({ 
    geneList = scGeneList(input$ATACd1inp, ATACgene) 
    if(nrow(geneList) > 50){ 
      HTML("More than 50 input genes! Please reduce the gene list!") 
    } else { 
      oup = paste0(nrow(geneList[present == TRUE]), " genes OK and will be plotted") 
      if(nrow(geneList[present == FALSE]) > 0){ 
        oup = paste0(oup, "<br/>", 
                     nrow(geneList[present == FALSE]), " genes not found (", 
                     paste0(geneList[present == FALSE]$gene, collapse = ", "), ")") 
      } 
      HTML(oup) 
    } 
  }) 
  output$ATACd1oup <- renderPlot({ 
    scBubbHeat(ATACconf, ATACmeta, input$ATACd1inp, input$ATACd1grp, input$ATACd1plt, 
               input$ATACd1sub1, input$ATACd1sub2, "ATACgexpr.h5", ATACgene, 
               input$ATACd1scl, input$ATACd1row, input$ATACd1col, 
               input$ATACd1cols, input$ATACd1fsz) 
  }) 
  output$ATACd1oup.ui <- renderUI({ 
    plotOutput("ATACd1oup", height = pList3[input$ATACd1psz]) 
  }) 
  output$ATACd1oup.pdf <- downloadHandler( 
    filename = function() { paste0("ATAC",input$ATACd1plt,"_",input$ATACd1grp,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$ATACd1oup.h, width = input$ATACd1oup.w, 
      plot = scBubbHeat(ATACconf, ATACmeta, input$ATACd1inp, input$ATACd1grp, input$ATACd1plt, 
                        input$ATACd1sub1, input$ATACd1sub2, "ATACgexpr.h5", ATACgene, 
                        input$ATACd1scl, input$ATACd1row, input$ATACd1col, 
                        input$ATACd1cols, input$ATACd1fsz, save = TRUE) ) 
  }) 
  output$ATACd1oup.png <- downloadHandler( 
    filename = function() { paste0("ATAC",input$ATACd1plt,"_",input$ATACd1grp,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$ATACd1oup.h, width = input$ATACd1oup.w, 
      plot = scBubbHeat(ATACconf, ATACmeta, input$ATACd1inp, input$ATACd1grp, input$ATACd1plt, 
                        input$ATACd1sub1, input$ATACd1sub2, "ATACgexpr.h5", ATACgene, 
                        input$ATACd1scl, input$ATACd1row, input$ATACd1col, 
                        input$ATACd1cols, input$ATACd1fsz, save = TRUE) ) 
  }) 
   
   
   optCrt="{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }" 
  updateSelectizeInput(session, "Peaksa1inp2", choices = names(Peaksgene), server = TRUE, 
                       selected = Peaksdef$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "Peaksa3inp1", choices = names(Peaksgene), server = TRUE, 
                       selected = Peaksdef$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "Peaksa3inp2", choices = names(Peaksgene), server = TRUE, 
                       selected = Peaksdef$gene2, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "Peaksb2inp1", choices = names(Peaksgene), server = TRUE, 
                       selected = Peaksdef$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "Peaksb2inp2", choices = names(Peaksgene), server = TRUE, 
                       selected = Peaksdef$gene2, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "Peaksc1inp2", server = TRUE, 
                       choices = c(Peaksconf[is.na(fID)]$UI,names(Peaksgene)), 
                       selected = Peaksconf[is.na(fID)]$UI[1], options = list( 
                         maxOptions = length(Peaksconf[is.na(fID)]$UI) + 3, 
                         create = TRUE, persist = TRUE, render = I(optCrt))) 
 
 
  ### Plots for tab a1 
  output$Peaksa1sub1.ui <- renderUI({ 
    sub = strsplit(Peaksconf[UI == input$Peaksa1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("Peaksa1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$Peaksa1sub1non, { 
    sub = strsplit(Peaksconf[UI == input$Peaksa1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "Peaksa1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$Peaksa1sub1all, { 
    sub = strsplit(Peaksconf[UI == input$Peaksa1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "Peaksa1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$Peaksa1oup1 <- renderPlot({ 
    scDRcell(Peaksconf, Peaksmeta, input$Peaksa1drX, input$Peaksa1drY, input$Peaksa1inp1,  
             input$Peaksa1sub1, input$Peaksa1sub2, 
             input$Peaksa1siz, input$Peaksa1col1, input$Peaksa1ord1, 
             input$Peaksa1fsz, input$Peaksa1asp, input$Peaksa1txt, input$Peaksa1lab1) 
  }) 
  output$Peaksa1oup1.ui <- renderUI({ 
    plotOutput("Peaksa1oup1", height = pList[input$Peaksa1psz]) 
  }) 
  output$Peaksa1oup1.pdf <- downloadHandler( 
    filename = function() { paste0("Peaks",input$Peaksa1drX,"_",input$Peaksa1drY,"_",  
                                   input$Peaksa1inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$Peaksa1oup1.h, width = input$Peaksa1oup1.w, useDingbats = FALSE, 
      plot = scDRcell(Peaksconf, Peaksmeta, input$Peaksa1drX, input$Peaksa1drY, input$Peaksa1inp1,   
                      input$Peaksa1sub1, input$Peaksa1sub2, 
                      input$Peaksa1siz, input$Peaksa1col1, input$Peaksa1ord1,  
                      input$Peaksa1fsz, input$Peaksa1asp, input$Peaksa1txt, input$Peaksa1lab1) ) 
  }) 
  output$Peaksa1oup1.png <- downloadHandler( 
    filename = function() { paste0("Peaks",input$Peaksa1drX,"_",input$Peaksa1drY,"_",  
                                   input$Peaksa1inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$Peaksa1oup1.h, width = input$Peaksa1oup1.w, 
      plot = scDRcell(Peaksconf, Peaksmeta, input$Peaksa1drX, input$Peaksa1drY, input$Peaksa1inp1,   
                      input$Peaksa1sub1, input$Peaksa1sub2, 
                      input$Peaksa1siz, input$Peaksa1col1, input$Peaksa1ord1,  
                      input$Peaksa1fsz, input$Peaksa1asp, input$Peaksa1txt, input$Peaksa1lab1) ) 
  }) 
  output$Peaksa1.dt <- renderDataTable({ 
    ggData = scDRnum(Peaksconf, Peaksmeta, input$Peaksa1inp1, input$Peaksa1inp2, 
                     input$Peaksa1sub1, input$Peaksa1sub2, 
                     "Peaksgexpr.h5", Peaksgene, input$Peaksa1splt) 
    datatable(ggData, rownames = FALSE, extensions = "Buttons", 
              options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>% 
      formatRound(columns = c("pctExpress"), digits = 2) 
  }) 
   
  output$Peaksa1oup2 <- renderPlot({ 
    scDRgene(Peaksconf, Peaksmeta, input$Peaksa1drX, input$Peaksa1drY, input$Peaksa1inp2,  
             input$Peaksa1sub1, input$Peaksa1sub2, 
             "Peaksgexpr.h5", Peaksgene, 
             input$Peaksa1siz, input$Peaksa1col2, input$Peaksa1ord2, 
             input$Peaksa1fsz, input$Peaksa1asp, input$Peaksa1txt) 
  }) 
  output$Peaksa1oup2.ui <- renderUI({ 
    plotOutput("Peaksa1oup2", height = pList[input$Peaksa1psz]) 
  }) 
  output$Peaksa1oup2.pdf <- downloadHandler( 
    filename = function() { paste0("Peaks",input$Peaksa1drX,"_",input$Peaksa1drY,"_",  
                                   input$Peaksa1inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$Peaksa1oup2.h, width = input$Peaksa1oup2.w, useDingbats = FALSE, 
      plot = scDRgene(Peaksconf, Peaksmeta, input$Peaksa1drX, input$Peaksa1drY, input$Peaksa1inp2,  
                      input$Peaksa1sub1, input$Peaksa1sub2, 
                      "Peaksgexpr.h5", Peaksgene, 
                      input$Peaksa1siz, input$Peaksa1col2, input$Peaksa1ord2, 
                      input$Peaksa1fsz, input$Peaksa1asp, input$Peaksa1txt) ) 
  }) 
  output$Peaksa1oup2.png <- downloadHandler( 
    filename = function() { paste0("Peaks",input$Peaksa1drX,"_",input$Peaksa1drY,"_",  
                                   input$Peaksa1inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$Peaksa1oup2.h, width = input$Peaksa1oup2.w, 
      plot = scDRgene(Peaksconf, Peaksmeta, input$Peaksa1drX, input$Peaksa1drY, input$Peaksa1inp2,  
                      input$Peaksa1sub1, input$Peaksa1sub2, 
                      "Peaksgexpr.h5", Peaksgene, 
                      input$Peaksa1siz, input$Peaksa1col2, input$Peaksa1ord2, 
                      input$Peaksa1fsz, input$Peaksa1asp, input$Peaksa1txt) ) 
  }) 
   
   
  ### Plots for tab a2 
  output$Peaksa2sub1.ui <- renderUI({ 
    sub = strsplit(Peaksconf[UI == input$Peaksa2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("Peaksa2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$Peaksa2sub1non, { 
    sub = strsplit(Peaksconf[UI == input$Peaksa2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "Peaksa2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$Peaksa2sub1all, { 
    sub = strsplit(Peaksconf[UI == input$Peaksa2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "Peaksa2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$Peaksa2oup1 <- renderPlot({ 
    scDRcell(Peaksconf, Peaksmeta, input$Peaksa2drX, input$Peaksa2drY, input$Peaksa2inp1,  
             input$Peaksa2sub1, input$Peaksa2sub2, 
             input$Peaksa2siz, input$Peaksa2col1, input$Peaksa2ord1, 
             input$Peaksa2fsz, input$Peaksa2asp, input$Peaksa2txt, input$Peaksa2lab1) 
  }) 
  output$Peaksa2oup1.ui <- renderUI({ 
    plotOutput("Peaksa2oup1", height = pList[input$Peaksa2psz]) 
  }) 
  output$Peaksa2oup1.pdf <- downloadHandler( 
    filename = function() { paste0("Peaks",input$Peaksa2drX,"_",input$Peaksa2drY,"_",  
                                   input$Peaksa2inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$Peaksa2oup1.h, width = input$Peaksa2oup1.w, useDingbats = FALSE, 
      plot = scDRcell(Peaksconf, Peaksmeta, input$Peaksa2drX, input$Peaksa2drY, input$Peaksa2inp1,   
                      input$Peaksa2sub1, input$Peaksa2sub2, 
                      input$Peaksa2siz, input$Peaksa2col1, input$Peaksa2ord1,  
                      input$Peaksa2fsz, input$Peaksa2asp, input$Peaksa2txt, input$Peaksa2lab1) ) 
  }) 
  output$Peaksa2oup1.png <- downloadHandler( 
    filename = function() { paste0("Peaks",input$Peaksa2drX,"_",input$Peaksa2drY,"_",  
                                   input$Peaksa2inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$Peaksa2oup1.h, width = input$Peaksa2oup1.w, 
      plot = scDRcell(Peaksconf, Peaksmeta, input$Peaksa2drX, input$Peaksa2drY, input$Peaksa2inp1,   
                      input$Peaksa2sub1, input$Peaksa2sub2, 
                      input$Peaksa2siz, input$Peaksa2col1, input$Peaksa2ord1,  
                      input$Peaksa2fsz, input$Peaksa2asp, input$Peaksa2txt, input$Peaksa2lab1) ) 
  }) 
   
  output$Peaksa2oup2 <- renderPlot({ 
    scDRcell(Peaksconf, Peaksmeta, input$Peaksa2drX, input$Peaksa2drY, input$Peaksa2inp2,  
             input$Peaksa2sub1, input$Peaksa2sub2, 
             input$Peaksa2siz, input$Peaksa2col2, input$Peaksa2ord2, 
             input$Peaksa2fsz, input$Peaksa2asp, input$Peaksa2txt, input$Peaksa2lab2) 
  }) 
  output$Peaksa2oup2.ui <- renderUI({ 
    plotOutput("Peaksa2oup2", height = pList[input$Peaksa2psz]) 
  }) 
  output$Peaksa2oup2.pdf <- downloadHandler( 
    filename = function() { paste0("Peaks",input$Peaksa2drX,"_",input$Peaksa2drY,"_",  
                                   input$Peaksa2inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$Peaksa2oup2.h, width = input$Peaksa2oup2.w, useDingbats = FALSE, 
      plot = scDRcell(Peaksconf, Peaksmeta, input$Peaksa2drX, input$Peaksa2drY, input$Peaksa2inp2,   
                      input$Peaksa2sub1, input$Peaksa2sub2, 
                      input$Peaksa2siz, input$Peaksa2col2, input$Peaksa2ord2,  
                      input$Peaksa2fsz, input$Peaksa2asp, input$Peaksa2txt, input$Peaksa2lab2) ) 
  }) 
  output$Peaksa2oup2.png <- downloadHandler( 
    filename = function() { paste0("Peaks",input$Peaksa2drX,"_",input$Peaksa2drY,"_",  
                                   input$Peaksa2inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$Peaksa2oup2.h, width = input$Peaksa2oup2.w, 
      plot = scDRcell(Peaksconf, Peaksmeta, input$Peaksa2drX, input$Peaksa2drY, input$Peaksa2inp2,   
                      input$Peaksa2sub1, input$Peaksa2sub2, 
                      input$Peaksa2siz, input$Peaksa2col2, input$Peaksa2ord2,  
                      input$Peaksa2fsz, input$Peaksa2asp, input$Peaksa2txt, input$Peaksa2lab2) ) 
  }) 
   
   
  ### Plots for tab a3 
  output$Peaksa3sub1.ui <- renderUI({ 
    sub = strsplit(Peaksconf[UI == input$Peaksa3sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("Peaksa3sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$Peaksa3sub1non, { 
    sub = strsplit(Peaksconf[UI == input$Peaksa3sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "Peaksa3sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$Peaksa3sub1all, { 
    sub = strsplit(Peaksconf[UI == input$Peaksa3sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "Peaksa3sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$Peaksa3oup1 <- renderPlot({ 
    scDRgene(Peaksconf, Peaksmeta, input$Peaksa3drX, input$Peaksa3drY, input$Peaksa3inp1,  
             input$Peaksa3sub1, input$Peaksa3sub2, 
             "Peaksgexpr.h5", Peaksgene, 
             input$Peaksa3siz, input$Peaksa3col1, input$Peaksa3ord1, 
             input$Peaksa3fsz, input$Peaksa3asp, input$Peaksa3txt) 
  }) 
  output$Peaksa3oup1.ui <- renderUI({ 
    plotOutput("Peaksa3oup1", height = pList[input$Peaksa3psz]) 
  }) 
  output$Peaksa3oup1.pdf <- downloadHandler( 
    filename = function() { paste0("Peaks",input$Peaksa3drX,"_",input$Peaksa3drY,"_",  
                                   input$Peaksa3inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$Peaksa3oup1.h, width = input$Peaksa3oup1.w, useDingbats = FALSE, 
      plot = scDRgene(Peaksconf, Peaksmeta, input$Peaksa3drX, input$Peaksa3drY, input$Peaksa3inp1,  
                      input$Peaksa3sub1, input$Peaksa3sub2, 
                      "Peaksgexpr.h5", Peaksgene, 
                      input$Peaksa3siz, input$Peaksa3col1, input$Peaksa3ord1, 
                      input$Peaksa3fsz, input$Peaksa3asp, input$Peaksa3txt) ) 
  }) 
  output$Peaksa3oup1.png <- downloadHandler( 
    filename = function() { paste0("Peaks",input$Peaksa3drX,"_",input$Peaksa3drY,"_",  
                                   input$Peaksa3inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$Peaksa3oup1.h, width = input$Peaksa3oup1.w, 
      plot = scDRgene(Peaksconf, Peaksmeta, input$Peaksa3drX, input$Peaksa3drY, input$Peaksa3inp1,  
                      input$Peaksa3sub1, input$Peaksa3sub2, 
                      "Peaksgexpr.h5", Peaksgene, 
                      input$Peaksa3siz, input$Peaksa3col1, input$Peaksa3ord1, 
                      input$Peaksa3fsz, input$Peaksa3asp, input$Peaksa3txt) ) 
  }) 
   
  output$Peaksa3oup2 <- renderPlot({ 
    scDRgene(Peaksconf, Peaksmeta, input$Peaksa3drX, input$Peaksa3drY, input$Peaksa3inp2,  
             input$Peaksa3sub1, input$Peaksa3sub2, 
             "Peaksgexpr.h5", Peaksgene, 
             input$Peaksa3siz, input$Peaksa3col2, input$Peaksa3ord2, 
             input$Peaksa3fsz, input$Peaksa3asp, input$Peaksa3txt) 
  }) 
  output$Peaksa3oup2.ui <- renderUI({ 
    plotOutput("Peaksa3oup2", height = pList[input$Peaksa3psz]) 
  }) 
  output$Peaksa3oup2.pdf <- downloadHandler( 
    filename = function() { paste0("Peaks",input$Peaksa3drX,"_",input$Peaksa3drY,"_",  
                                   input$Peaksa3inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$Peaksa3oup2.h, width = input$Peaksa3oup2.w, useDingbats = FALSE, 
      plot = scDRgene(Peaksconf, Peaksmeta, input$Peaksa3drX, input$Peaksa3drY, input$Peaksa3inp2,  
                      input$Peaksa3sub1, input$Peaksa3sub2, 
                      "Peaksgexpr.h5", Peaksgene, 
                      input$Peaksa3siz, input$Peaksa3col2, input$Peaksa3ord2, 
                      input$Peaksa3fsz, input$Peaksa3asp, input$Peaksa3txt) ) 
  }) 
  output$Peaksa3oup2.png <- downloadHandler( 
    filename = function() { paste0("Peaks",input$Peaksa3drX,"_",input$Peaksa3drY,"_",  
                                   input$Peaksa3inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$Peaksa3oup2.h, width = input$Peaksa3oup2.w, 
      plot = scDRgene(Peaksconf, Peaksmeta, input$Peaksa3drX, input$Peaksa3drY, input$Peaksa3inp2,  
                      input$Peaksa3sub1, input$Peaksa3sub2, 
                      "Peaksgexpr.h5", Peaksgene, 
                      input$Peaksa3siz, input$Peaksa3col2, input$Peaksa3ord2, 
                      input$Peaksa3fsz, input$Peaksa3asp, input$Peaksa3txt) ) 
  }) 
     
   
  ### Plots for tab b2 
  output$Peaksb2sub1.ui <- renderUI({ 
    sub = strsplit(Peaksconf[UI == input$Peaksb2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("Peaksb2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$Peaksb2sub1non, { 
    sub = strsplit(Peaksconf[UI == input$Peaksb2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "Peaksb2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$Peaksb2sub1all, { 
    sub = strsplit(Peaksconf[UI == input$Peaksb2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "Peaksb2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$Peaksb2oup1 <- renderPlot({ 
    scDRcoex(Peaksconf, Peaksmeta, input$Peaksb2drX, input$Peaksb2drY,   
             input$Peaksb2inp1, input$Peaksb2inp2, input$Peaksb2sub1, input$Peaksb2sub2, 
             "Peaksgexpr.h5", Peaksgene, 
             input$Peaksb2siz, input$Peaksb2col1, input$Peaksb2ord1, 
             input$Peaksb2fsz, input$Peaksb2asp, input$Peaksb2txt) 
  }) 
  output$Peaksb2oup1.ui <- renderUI({ 
    plotOutput("Peaksb2oup1", height = pList2[input$Peaksb2psz]) 
  }) 
  output$Peaksb2oup1.pdf <- downloadHandler( 
    filename = function() { paste0("Peaks",input$Peaksb2drX,"_",input$Peaksb2drY,"_",  
                                    input$Peaksb2inp1,"_",input$Peaksb2inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$Peaksb2oup1.h, width = input$Peaksb2oup1.w, useDingbats = FALSE, 
      plot = scDRcoex(Peaksconf, Peaksmeta, input$Peaksb2drX, input$Peaksb2drY,  
                      input$Peaksb2inp1, input$Peaksb2inp2, input$Peaksb2sub1, input$Peaksb2sub2, 
                      "Peaksgexpr.h5", Peaksgene, 
                      input$Peaksb2siz, input$Peaksb2col1, input$Peaksb2ord1, 
                      input$Peaksb2fsz, input$Peaksb2asp, input$Peaksb2txt) ) 
  }) 
  output$Peaksb2oup1.png <- downloadHandler( 
    filename = function() { paste0("Peaks",input$Peaksb2drX,"_",input$Peaksb2drY,"_",  
                                    input$Peaksb2inp1,"_",input$Peaksb2inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$Peaksb2oup1.h, width = input$Peaksb2oup1.w, 
      plot = scDRcoex(Peaksconf, Peaksmeta, input$Peaksb2drX, input$Peaksb2drY,  
                      input$Peaksb2inp1, input$Peaksb2inp2, input$Peaksb2sub1, input$Peaksb2sub2, 
                      "Peaksgexpr.h5", Peaksgene, 
                      input$Peaksb2siz, input$Peaksb2col1, input$Peaksb2ord1, 
                      input$Peaksb2fsz, input$Peaksb2asp, input$Peaksb2txt) ) 
  }) 
  output$Peaksb2oup2 <- renderPlot({ 
    scDRcoexLeg(input$Peaksb2inp1, input$Peaksb2inp2, input$Peaksb2col1, input$Peaksb2fsz) 
  }) 
  output$Peaksb2oup2.ui <- renderUI({ 
    plotOutput("Peaksb2oup2", height = "300px") 
  }) 
  output$Peaksb2oup2.pdf <- downloadHandler( 
    filename = function() { paste0("Peaks",input$Peaksb2drX,"_",input$Peaksb2drY,"_",  
                                    input$Peaksb2inp1,"_",input$Peaksb2inp2,"_leg.pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = 3, width = 4, useDingbats = FALSE, 
      plot = scDRcoexLeg(input$Peaksb2inp1, input$Peaksb2inp2, input$Peaksb2col1, input$Peaksb2fsz) ) 
  }) 
  output$Peaksb2oup2.png <- downloadHandler( 
    filename = function() { paste0("Peaks",input$Peaksb2drX,"_",input$Peaksb2drY,"_",  
                                    input$Peaksb2inp1,"_",input$Peaksb2inp2,"_leg.png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = 3, width = 4, 
      plot = scDRcoexLeg(input$Peaksb2inp1, input$Peaksb2inp2, input$Peaksb2col1, input$Peaksb2fsz) ) 
  }) 
  output$Peaksb2.dt <- renderDataTable({ 
    ggData = scDRcoexNum(Peaksconf, Peaksmeta, input$Peaksb2inp1, input$Peaksb2inp2, 
                         input$Peaksb2sub1, input$Peaksb2sub2, "Peaksgexpr.h5", Peaksgene) 
    datatable(ggData, rownames = FALSE, extensions = "Buttons", 
              options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>% 
      formatRound(columns = c("percent"), digits = 2) 
  }) 
     
   
  ### Plots for tab c1 
  output$Peaksc1sub1.ui <- renderUI({ 
    sub = strsplit(Peaksconf[UI == input$Peaksc1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("Peaksc1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$Peaksc1sub1non, { 
    sub = strsplit(Peaksconf[UI == input$Peaksc1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "Peaksc1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$Peaksc1sub1all, { 
    sub = strsplit(Peaksconf[UI == input$Peaksc1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "Peaksc1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$Peaksc1oup <- renderPlot({ 
    scVioBox(Peaksconf, Peaksmeta, input$Peaksc1inp1, input$Peaksc1inp2, 
             input$Peaksc1sub1, input$Peaksc1sub2, 
             "Peaksgexpr.h5", Peaksgene, input$Peaksc1typ, input$Peaksc1pts, 
             input$Peaksc1siz, input$Peaksc1fsz) 
  }) 
  output$Peaksc1oup.ui <- renderUI({ 
    plotOutput("Peaksc1oup", height = pList2[input$Peaksc1psz]) 
  }) 
  output$Peaksc1oup.pdf <- downloadHandler( 
    filename = function() { paste0("Peaks",input$Peaksc1typ,"_",input$Peaksc1inp1,"_",  
                                   input$Peaksc1inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$Peaksc1oup.h, width = input$Peaksc1oup.w, useDingbats = FALSE, 
      plot = scVioBox(Peaksconf, Peaksmeta, input$Peaksc1inp1, input$Peaksc1inp2, 
                      input$Peaksc1sub1, input$Peaksc1sub2, 
                      "Peaksgexpr.h5", Peaksgene, input$Peaksc1typ, input$Peaksc1pts, 
                      input$Peaksc1siz, input$Peaksc1fsz) ) 
  }) 
  output$Peaksc1oup.png <- downloadHandler( 
    filename = function() { paste0("Peaks",input$Peaksc1typ,"_",input$Peaksc1inp1,"_",  
                                   input$Peaksc1inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$Peaksc1oup.h, width = input$Peaksc1oup.w, 
      plot = scVioBox(Peaksconf, Peaksmeta, input$Peaksc1inp1, input$Peaksc1inp2, 
                      input$Peaksc1sub1, input$Peaksc1sub2, 
                      "Peaksgexpr.h5", Peaksgene, input$Peaksc1typ, input$Peaksc1pts, 
                      input$Peaksc1siz, input$Peaksc1fsz) ) 
  }) 
     
   
### Plots for tab c2 
  output$Peaksc2sub1.ui <- renderUI({ 
    sub = strsplit(Peaksconf[UI == input$Peaksc2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("Peaksc2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$Peaksc2sub1non, { 
    sub = strsplit(Peaksconf[UI == input$Peaksc2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "Peaksc2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$Peaksc2sub1all, { 
    sub = strsplit(Peaksconf[UI == input$Peaksc2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "Peaksc2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
output$Peaksc2oup <- renderPlot({ 
  scProp(Peaksconf, Peaksmeta, input$Peaksc2inp1, input$Peaksc2inp2,  
         input$Peaksc2sub1, input$Peaksc2sub2, 
         input$Peaksc2typ, input$Peaksc2flp, input$Peaksc2fsz) 
}) 
output$Peaksc2oup.ui <- renderUI({ 
  plotOutput("Peaksc2oup", height = pList2[input$Peaksc2psz]) 
}) 
output$Peaksc2oup.pdf <- downloadHandler( 
  filename = function() { paste0("Peaks",input$Peaksc2typ,"_",input$Peaksc2inp1,"_",  
                                 input$Peaksc2inp2,".pdf") }, 
  content = function(file) { ggsave( 
    file, device = "pdf", height = input$Peaksc2oup.h, width = input$Peaksc2oup.w, useDingbats = FALSE, 
    plot = scProp(Peaksconf, Peaksmeta, input$Peaksc2inp1, input$Peaksc2inp2,  
                  input$Peaksc2sub1, input$Peaksc2sub2, 
                  input$Peaksc2typ, input$Peaksc2flp, input$Peaksc2fsz) ) 
  }) 
output$Peaksc2oup.png <- downloadHandler( 
  filename = function() { paste0("Peaks",input$Peaksc2typ,"_",input$Peaksc2inp1,"_",  
                                 input$Peaksc2inp2,".png") }, 
  content = function(file) { ggsave( 
    file, device = "png", height = input$Peaksc2oup.h, width = input$Peaksc2oup.w, 
    plot = scProp(Peaksconf, Peaksmeta, input$Peaksc2inp1, input$Peaksc2inp2,  
                  input$Peaksc2sub1, input$Peaksc2sub2, 
                  input$Peaksc2typ, input$Peaksc2flp, input$Peaksc2fsz) ) 
  }) 
     
   
  ### Plots for tab d1 
  output$Peaksd1sub1.ui <- renderUI({ 
    sub = strsplit(Peaksconf[UI == input$Peaksd1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("Peaksd1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$Peaksd1sub1non, { 
    sub = strsplit(Peaksconf[UI == input$Peaksd1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "Peaksd1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$Peaksd1sub1all, { 
    sub = strsplit(Peaksconf[UI == input$Peaksd1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "Peaksd1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$Peaksd1oupTxt <- renderUI({ 
    geneList = scGeneList(input$Peaksd1inp, Peaksgene) 
    if(nrow(geneList) > 50){ 
      HTML("More than 50 input genes! Please reduce the gene list!") 
    } else { 
      oup = paste0(nrow(geneList[present == TRUE]), " genes OK and will be plotted") 
      if(nrow(geneList[present == FALSE]) > 0){ 
        oup = paste0(oup, "<br/>", 
                     nrow(geneList[present == FALSE]), " genes not found (", 
                     paste0(geneList[present == FALSE]$gene, collapse = ", "), ")") 
      } 
      HTML(oup) 
    } 
  }) 
  output$Peaksd1oup <- renderPlot({ 
    scBubbHeat(Peaksconf, Peaksmeta, input$Peaksd1inp, input$Peaksd1grp, input$Peaksd1plt, 
               input$Peaksd1sub1, input$Peaksd1sub2, "Peaksgexpr.h5", Peaksgene, 
               input$Peaksd1scl, input$Peaksd1row, input$Peaksd1col, 
               input$Peaksd1cols, input$Peaksd1fsz) 
  }) 
  output$Peaksd1oup.ui <- renderUI({ 
    plotOutput("Peaksd1oup", height = pList3[input$Peaksd1psz]) 
  }) 
  output$Peaksd1oup.pdf <- downloadHandler( 
    filename = function() { paste0("Peaks",input$Peaksd1plt,"_",input$Peaksd1grp,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$Peaksd1oup.h, width = input$Peaksd1oup.w, 
      plot = scBubbHeat(Peaksconf, Peaksmeta, input$Peaksd1inp, input$Peaksd1grp, input$Peaksd1plt, 
                        input$Peaksd1sub1, input$Peaksd1sub2, "Peaksgexpr.h5", Peaksgene, 
                        input$Peaksd1scl, input$Peaksd1row, input$Peaksd1col, 
                        input$Peaksd1cols, input$Peaksd1fsz, save = TRUE) ) 
  }) 
  output$Peaksd1oup.png <- downloadHandler( 
    filename = function() { paste0("Peaks",input$Peaksd1plt,"_",input$Peaksd1grp,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$Peaksd1oup.h, width = input$Peaksd1oup.w, 
      plot = scBubbHeat(Peaksconf, Peaksmeta, input$Peaksd1inp, input$Peaksd1grp, input$Peaksd1plt, 
                        input$Peaksd1sub1, input$Peaksd1sub2, "Peaksgexpr.h5", Peaksgene, 
                        input$Peaksd1scl, input$Peaksd1row, input$Peaksd1col, 
                        input$Peaksd1cols, input$Peaksd1fsz, save = TRUE) ) 
  }) 
   
   
      
}) 
 
 
 
 