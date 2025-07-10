
# Instruction: The function of this script is to creat QuSAGE polt based on rds file and genesets

library(qusage)
library(Seurat)

#Fuctions: Load input files online
LoadInputFile <- function () {
  crtwd <<- getwd()
  print(paste("current work directory:", crtwd, sep = " "))
  args <- commandArgs(T)
  print(args)
  species <<- args[1]
  file_rds <<- args[2]
  setwd(args[3])
  crtwd <<- getwd()
  print(paste("current work directory:", crtwd, sep = " "))
  input = c()
  if ( species == "human" ) {
    if (grepl("MSigDB$", args[4], perl = T)) {
      gmtDB <- unlist(strsplit(args[5], "[,]"))
      gmtDB <- paste(args[4], "_human/", gmtDB, sep = "")
    } else {
        input <- unlist(strsplit(args[4], "[,]"))
        if (grepl("MSigDB$", args[5], perl = T)) {
          gmtDB <- unlist(strsplit(args[6], "[,]"))
          gmtDB <- paste(args[5], "_human/", gmtDB, sep = "")
        }
      }
  }
  if ( species == "rat" ) {
    if (grepl("MSigDB$", args[4], perl = T)) {
      gmtDB <- unlist(strsplit(args[5], "[,]"))
      gmtDB <- paste(args[4], "_rat/", gmtDB, sep = "")
    } else {
        input <- unlist(strsplit(args[4], "[,]"))
        if (grepl("MSigDB$", args[5], perl = T)) {
          gmtDB <- unlist(strsplit(args[6], "[,]"))
          gmtDB <- paste(args[5], "_rat/", gmtDB, sep = "")
        }
      }
  }
  #print(class(gmtDB))
  print(class(input))
  #list_gmt <<- c(input, gmtDB)
  list_gmt <<- c(input)
  print(list_gmt)
}

#Fuctions: Load input files offline
LoadInputFile.offline <- function () {
  crtwd <<- getwd()
  print(paste("current work directory:", crtwd, sep = " "))
  args <- commandArgs(T)
  print(args)
  file_rds <<- paste(crtwd, args[1], sep = "/")
  list_gmt <<- args[2:length(args)]
  list_gmt <<- paste(crtwd, list_gmt, sep = "/")
  print(list_gmt)
}

# Function: Get expression matrix and cluster list from rds file
# Input:    file_rds
# Output:   eset (expression matrix)
#           ident (list of clusters for each cell)
GetInfoFromRDS <- function (file_rds) {
  seuset <- readRDS(file_rds)
  date()
  eset <<- as.matrix(seuset@data)
  ident <<- as.character(seuset@ident)
  date()
  print("loading rds file finish")
}

# Function: Write gene set info into MSIG.genesets from gmt file or txt file with GeneSets info
# Input:    file_gmt
# Output:   MSIG.genesets
GetGeneSets <- function (file_gmt) {
  #print(paste('file_gmt:', file_gmt))
  if (grepl("txt$", file_gmt, perl = T, ignore.case = T)) {
    info <- read.table(file_gmt)
    colnames(info) <- c("type", "gene")
    MSIG.genesets <<- list()
    for (i in unique(info$type)) {
      set <- list(as.character(info[info$type==i,2]))
      MSIG.genesets <<- c(MSIG.genesets, set)
    }
    names(MSIG.genesets) <<- unique(info$type)
    #print(MSIG.genesets)
  }
  if ( grepl("gmt$", file_gmt, perl = T, ignore.case = T)) {
    MSIG.genesets <<- read.gmt(file_gmt)
    length(MSIG.genesets)
  }
  #MSIG.genesets
}


# Function: Create color schemes for plot
SetPlotArg <- function () {
  colors40 <<- c("darkblue", "orange", "green3", "red", "purple", "saddlebrown",
		"pink1", "darkgrey", "yellow", "deepskyblue", "darkgreen",
		"blue", "darkmagenta", "red4", "yellow4", "aquamarine2",
		"magenta1", "lemonchiffon1", "grey21", "tomato1", "hotpink4",
		"deeppink1", "deepskyblue4", "slategrey", "greenyellow",
		"mediumpurple", "darkorange1", "orchid3", "peachpuff1",
		"lightgreen", "cornflowerblue", "brown3", "tan", "darkorchid4",
		"hotpink", "blue3", "chartreuse4", "khaki1", "cyan", "salmon4")
  colors2  <<- c("black", rainbow(60))
  colors80 <<- c(colors40, rep("grey", 40))
  lwds50 <<- c(rep(5,49), 10)
}

# Instruction: This function is modified code based on plotDensityCurve from qusage library
plotDC <- function (QSarray, path.index = 1:numPathways(QSarray), zeroLine = TRUE,
		    addVIF = !is.null(QSarray$vif), col = NULL, plot = TRUE,
		    add = FALSE, xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL,
		    type = "l", ...) {
    if (is.character(path.index)) {
        path.index = match(path.index, names(QSarray$pathways))
    }
    if (all(is.na(QSarray$path.mean[path.index]))) {
        stop("no non-missing pathways in path.index")
    }
    scaleFactor = pdfScaleFactor(QSarray, addVIF = addVIF)
    if (is.null(xlim)) {
        xlim = range(calcBayesCI(QSarray, addVIF = addVIF)[, path.index], na.rm = T)
        xlim = xlim + (xlim - mean(xlim)) * 2
    }
    if (is.null(ylim)) {
        ylim = range(t(QSarray$path.PDF[, path.index])/scaleFactor[path.index], na.rm = T)
    }
    if (is.null(xlab))  xlab = "Gene Set Activity"
    if (is.null(ylab))  ylab = "Density"
    if (is.null(col))   col = par("col")
    if (!add & plot) {
        plot(0, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, ...)
    }
    if (length(col) != length(path.index)) {
        col = rep(col, length.out = length(path.index))
    }
    if (zeroLine & plot) abline(v = 0, lty = 2)
    retVal = list()
    for (i in 1:length(path.index)) {
        path = path.index[i]
        x = getXcoords(QSarray, path, addVIF = addVIF)
        y = QSarray$path.PDF[, path]/scaleFactor[path]
        if (plot) {
            lines(x, y, col = col[i], type = type, ...)
        }
        retVal[[i]] = data.frame(x, y)
    }
    names(retVal) = colnames(QSarray$path.PDF)[path.index]
    invisible(retVal)
}

# Instruction: This function is the original code from qusage library
pdfScaleFactor <- function(QSarray, addVIF=!is.null(QSarray$vif)){
  if(is.null(QSarray$vif) && addVIF){stop("vif is undefined for QSarray object. addVIF can not be set to true.")}
  sif = sapply(1:numPathways(QSarray),function(i){ifelse(addVIF,sqrt(QSarray$vif[i]),1)})
  sif[is.na(sif)] = 1
  pdfSum = colSums(QSarray$path.PDF)

  ##get the pdf range
  ranges = QSarray$ranges * 2

  ##the scale factor is essentially the distance between points in the x coordinates times the (current) sum of the pdf
  scaleFactor = (ranges*sif) / (QSarray$n.points-1) * pdfSum
  scaleFactor
}

# Function: Create QuSAGE results include DensityCurvesPlots, ConfidenceIntervalPlots and DataTable based on qs.results
# Input:    qs.results (output of qusage)
#           table_name (file name of DataTable)
#           plotDC_name(file name of DensityCurvesPlot)
#           plotCI_name(file name of ConfidenceIntervaPlots)
#           cluster    (cluster ID)
GetResultFiles <- function (
    qs.results,
    table_name,
    plotDC_name,
    plotCI_name,
    cluster
  ) {
  clusID <- paste("Cluster", cluster)
  p.vals <- pdf.pVal(qs.results)
  #class(p.vals)
  #print(p.vals)
  q.vals <- p.adjust(p.vals, method="fdr")
  #class(q.vals)
  topSets <- order(p.vals)
  #class(topSets)
  #print(topSets[1:25])

  #print(MSIG.genesets)
  table <- qsTable(qs.results, number = (length(MSIG.genesets)))
  #print(table[(1:5),(1:3)])
  write.table(table, file=table_name, sep="\t", row.names = F)
  table.order <- as.numeric(row.names(table))
  #print(table.order)


  ordered.legend = character()
  for (i in 1:length(table.order)) {
    ordered.legend[i] <- names(MSIG.genesets)[(table.order[i])]
    if (i > 5) next
    print(paste(table.order[i],names(MSIG.genesets)[table.order[i]],ordered.legend[i],sep=" "))
  }
  #print(names(MSIG.genesets))
  #print(ordered.legend)

  curve.num = length(table.order)
  top.num = curve.num

  png(file=plotDC_name, width = 1800, height = 1000)
  par(mgp = c(3,1,0), mar=c(6,6,3,2))
  plotDC(qs.results, path.index=table.order[top.num:1], lwd=3, col=colors80[top.num:1],
		    cex.axis=2, cex.xaxis=2, cex.lab=2, main=clusID, cex.main=2)
  plotDC(qs.results, path.index=table.order[1], lwd=9, col=colors80[1],
		    cex.axis=2, cex.xaxis=2, cex.lab=2,
		    add=T)
  legend("topleft", legend = ordered.legend[1:top.num], text.col = colors80, bty="n")#xpd=T)
  dev.off()

  pdf(file=sub("png", "pdf", plotDC_name), width = 27, height = 15)
  par(mgp = c(3,1,0), mar=c(6,6,3,3))
  plotDensityCurves(qs.results, path.index=table.order[top.num:1], lwd=5, col=colors80[top.num:1],
		    cex.axis=2, cex.xaxis=2, cex.lab=2,
		    xlab="Gene Set Activity", ylab="Density", main=clusID, cex.main=2)
  plotDensityCurves(qs.results, path.index=table.order[1],
		    lwd=10, col=colors80[1], cex.axis=2, cex.xaxis=2, cex.lab=2,
		    add=T)
  legend("topleft", legend = ordered.legend[1:top.num], text.col = colors80,
	 bty="o", bg="transparent") #xpd = T)
  dev.off()


  png(file=plotCI_name, width = 1800, height = 1000)
  par(mgp = c(6,1,0))
  plotCIs(qs.results, path.index=topSets[1:top.num], mar=c(30,15,3,1), lwd=3,
	  col=colors2[1:top.num], cex.axis=2, cex.xaxis=3, cex.lab=2,
	  ylab="Gene Set Activity", main=clusID, cex.main=2)
  dev.off()
}

# Function: Set output file names inclue DensityCurvesPlot, ConfidenceIntervalPlot and DataTable
# Input:    cluster (cluster name)
SetResultFileName <- function (cluster) {
  label_rds <<- "TestData"
  label_gmt <<- sub(".gmt$", "", file_gmt, perl = T, ignore.case = T)
  label_gmt <<- sub(".symbols$", "", label_gmt, perl = T, ignore.case = T)
  label_gmt <<- sub(".v6.2$", "", label_gmt, perl = T, ignore.case = T)
  label_gmt <<- sub("^.+/", "", label_gmt, perl = T, ignore.case = T)
  #print(paste('label_gmt:', label_gmt))
  label <- paste(label_gmt, "Cluster", cluster, sep = "_")
  outwd <- paste(crtwd, "qusage_result", label_gmt, sep = "/")
  if (!file.exists(outwd)) dir.create(outwd, recursive = T)
  table_name  <<- paste(outwd, "/", "table_", label, ".txt", sep = "")
  plotDC_name <<- paste(outwd, "/", "plotDC_", label, ".png", sep = "")
  #DC_legend   <<- paste(outwd, "/", "plotDC_legend_", label, ".png", sep = "")
  plotCI_name <<- paste(outwd, "/", "plotCI_", label, ".png", sep = "")
}

main <- function () {
  print(date())
  date()
  LoadInputFile()
  GetInfoFromRDS (file_rds)
  SetPlotArg()
  for (f in list_gmt) {
    file_gmt <<- f
    print(file_gmt)
    GetGeneSets(file_gmt)
    for (i in unique(ident)) {
      print(paste("Cluster", i, date(), sep = " "))
      labels <- ident
      labels[labels!=i] = "B"
      labels[labels==i] = "C"
      contrast <- "C-B"
      qs.results <- qusage(eset, labels, contrast, MSIG.genesets)
      #class(qs.results)
      #head(qs.results)
      SetResultFileName(i)
      GetResultFiles(qs.results, table_name, plotDC_name, plotCI_name, i)
    }
    date()
  }
}

main()
