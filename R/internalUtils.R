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
