###############################################################################
#' Precisely identify TSSs from bam files, paired end bam files, bed files,
#' BigWig files, tss files, or tss tables.
#'
#' @description getTSS function is used to precisely identify TSSs from multiple
#' input file formats. The files include users' home-made alignment files (bam format)
#' or downloaded files from public databases. See inputFilesType for details on
#' the supported input file formats.
#'
#' @usage getTSS(object, sequencingQualityThreshold = 10, 
#' mappingQualityThreshold = 20, softclippingAllowed = FALSE)
#'
#' @param object A TSSr object.
#' @param sequencingQualityThreshold Used only if inputFilesType == "bam" or
#' "bamPairedEnd", otherwise ignored.
#' @param mappingQualityThreshold Used only if inputFilesType == "bam" or
#' "bamPairedEnd", otherwise ignored.
#' @param softclippingAllowed Used only if inputFilesType == "bam" or
#' "bamPairedEnd". Default is FALSE.
#' @return Large List of elements - one element for each sample
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(exampleTSSr)
#' #getTSS(exampleTSSr)
#' }
setGeneric("getTSS",function(object
                             ,sequencingQualityThreshold = 10
                             ,mappingQualityThreshold = 20
                             ,softclippingAllowed = FALSE
                             ,useMultiCore=FALSE
                             ,numCores=NULL)standardGeneric("getTSS"))
#' @rdname getTSS
#' @export
setMethod("getTSS",signature(object = "TSSr"), function(object
                                    ,sequencingQualityThreshold
                                    ,mappingQualityThreshold
                                    ,softclippingAllowed
                                    ,useMultiCore
                                    ,numCores){
  ##initialize values
  Genome <- .getGenome(object@genomeName)
  sampleLabels <- object@sampleLabels
  inputFilesType <- object@inputFilesType
  if (length(object@sampleLabelsMerged) == 0) {
    object@sampleLabelsMerged <- sampleLabels
  }
  objName <- deparse(substitute(object))
  if(inputFilesType == "bam" | inputFilesType == "bamPairedEnd"){

#### MULTICORE CODE STARTS HERE (only for BAM files for now) ####

    if(useMultiCore) {
      if (is.null(numCores)) {
         numCores <- detectCores()
      }
      print(paste("process is running on", numCores, "cores..."))      
    
      inputFilesID_forMultiCore <- as.list(c(1:length(object@inputFiles)))

      results <- mclapply(inputFilesID_forMultiCore
                          ,function(x) {
                            tssMC <- .getTSS_from_bam(object@inputFiles[x]
                            ,Genome
                            ,sampleLabels[x]
                            ,inputFilesType
                            ,sequencingQualityThreshold
                            ,mappingQualityThreshold
                            ,softclippingAllowed)
                          }
                          ,mc.cores = numCores)

#### MERGE THE TABLES FROM SEPARATE CORES ####

    tss <- NULL
    
      for (i in results) {
        if (is.null(tss)) {
          tss = i
        } else {
          tss = merge(tss, i, by = c("chr", "pos", "strand"),all=TRUE)
        }
      }
    
    tss=tss[order(tss$chr, tss$pos),]

    tss[is.na(tss)] = 0

#### END OF MERGE ####

#### MULTICORE CODE ENDS HERE ####

    } else {

#### original code if useMultiCore is FALSE ####

    tss <- .getTSS_from_bam(object@inputFiles
                     ,Genome
                     ,sampleLabels
                     ,inputFilesType
                     ,sequencingQualityThreshold
                     ,mappingQualityThreshold
                     ,softclippingAllowed)

    } #### end of else

  }else if(inputFilesType == "bed"){
    tss <- .getTSS_from_bed(object@inputFiles, Genome, sampleLabels)
  }else if(inputFilesType == "BigWig"){
    tss <- .getTSS_from_BigWig(object@inputFiles,Genome, sampleLabels)
  }else if(inputFilesType == "tss"){
    tss <- .getTSS_from_tss(object@inputFiles, sampleLabels)
  }else if(inputFilesType == "TSStable"){
    tss <- .getTSS_from_TSStable(object@inputFiles, sampleLabels)
  }

#################################################################################

  setorder(tss, "strand","chr","pos")
  # get library sizes
  object@librarySizes <- colSums(tss[,4:ncol(tss), drop = FALSE], na.rm = TRUE)

  object@TSSrawMatrix <- tss
  object@TSSprocessedMatrix <- tss
  assign(objName, object, envir = parent.frame())
})


