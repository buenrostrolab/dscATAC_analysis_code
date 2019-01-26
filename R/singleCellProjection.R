library(SummarizedExperiment)

projectCells <- function(scSE, # SE object of scATAC-seq cell profiles to score/project
                 refPCs, # PC loadings derived using bulk ATAC-seq counts (Peaks x PCs). Peaks should match with those in single cell SE object
                 doInChunks=TRUE # Split up the work in peak-wise chunks?
         
){
  
  if(nrow(scSE)!=nrow(refPCs))
    stop("sc-ATAC-seq counts object and PC reference do not have matching dimensions..\\n")
  
  if(doInChunks){
    # Run sequentially and merge results at the end
    chunksize <- 1000
    cat("Computing scores for cells sequentially in groups of size ",chunksize," ..\\n")
    starts <- seq(1,ncol(scSE),chunksize)
    
  } else {
    starts <- 1
  }
  
  
  scores.l <- list()
  
  for(i in 1:length(starts)){
    beginning <- starts[i]
    if(i==length(starts)) # If it's the last one
    {
      ending <- ncol(scSE)
    } else {
      ending <- starts[i]+chunksize-1
    }
    
    cat("Computing scores for cells: ",beginning," to ", ending,"..\\n")
    
    m <- SummarizedExperiment::assay(scSE[,beginning:ending])
    
    cellMeans <- Matrix::colMeans(m)
    
    # Center cell counts based on its mean FRIP count
    
    cCounts <- Matrix::t(Matrix::t(m)/cellMeans)
    
    # Score cells based on the ATAC reference
    # Matrix multiplication
    # Notation is (Rows x Columns)
    # (Cells x Peaks) * (Peaks x PC loadings) = (Cells X PC scores)
    
    scores.l[[i]] <- as.matrix(Matrix::t(cCounts) %*% refPCs)
    
    gc()
  }
  
  scATAC.scores <- plyr::ldply(scores.l,rbind)
  cat("Done!\\n\\n")
  
  rownames(scATAC.scores) <- colnames(scSE)
  
  return(as.matrix(scATAC.scores))
}


# Single cell ATAC-seq counts SummarizedExperimen object
countsFileUrl <- "https://github.com/buenrostrolab/dscATAC_analysis_code/blob/master/data/jdb-cell_CD34_C1sorted_matchedpeaks.rds"
# PC Loadings from bulk ATAC-seq counts (reference cell types)
pcFileUrl <- "https://github.com/buenrostrolab/dscATAC_analysis_code/blob/master/data/pcLoadings.rds"

# Load single cell counts object
SE <- readRDS(url(countsFileUrl))
SE

# Load bulk reference PC loadings
pcLoadings <- readRDS(url(pcFileUrl))
print(dim(pcLoadings))

# Project single cell data in PC space to get single cell scores
CD34JDBscScores <- projectCells(scSE = SE,refPCs = pcLoadings,doInChunks = TRUE)
