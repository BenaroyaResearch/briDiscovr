# Private function to read in data from FCS files and associate the data with
# subject IDs. Accepts a data frame with at least two columns:'filename',
# which should contain paths to .fcs files, and the field indicated by
# 'indexField' which will be used as the index in the returned list
buildFcsList <- function(fileFrame, indexField = "subject", ...) {
  stopifnot(all(c(indexField, "filename") %in% colnames(fileFrame)))
  stopifnot(is.data.frame(fileFrame))
  tmpList <- list()
  for (currIdx in unique(fileFrame[[indexField]])) {
    currFcsName = as.character(fileFrame[which(fileFrame[[indexField]]==currIdx), "filename"])
    if(file.info(currFcsName)$size > 0){
      tmpList[[currIdx]] <- flowCore::read.FCS(currFcsName, transformation = FALSE, ...)
    }
  }
  return(tmpList)
}

# Private function to generate a list of colors of arbitrary length.
getColorList <- function(n, seed = 42){
  set.seed(seed)
  qualColorPals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  allColors = unlist(mapply(RColorBrewer::brewer.pal, qualColorPals$maxcolors, rownames(qualColorPals)))

  retList = c()

  if (n > length(allColors)){
    message("Note: ", n, " colors requested, ", length(allColors), " unique colors available. Some colors will be repeated.")
    while(n > length(allColors)){
      retList = c(retList, sample(allColors))
      n = n - length(allColors)
    }
  }
  retList = c(retList, sample(allColors)[1:n])
  return(retList)
}
