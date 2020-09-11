#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

##########################
#Functions include breakpoint_chrom_position, unknown_chrom_position, breakpoint_counts, cell_counts -> stored in R script recombination_rates_functions.R
#Global variable -> vcf file that is stored as "df"
#Functions take df as an argument, which specifies the indiviudal and chromosome it orginated from
#Calculates the coordinates of breakpoints and unknown blocks, and then determines the breakpoints and cell counts per 1Mb bin


#Created by Charlotte Michelle Nguyen
#November, 2017

##########################

##########################
###
###
###    Functions
###
###
#########################
#Computes coordinates of breakpoints in every cell of file
breakpoint_chrom_position <- function(v) { #v - vector with the parental states within a cell; it takes the df as global variable.
  
  #Test function with one cell in vcf file
  #v <- df[,3]
  
  #Stores first and last coordinate of each M/P block
  bp <- c()
  #Neutral parental state of the chromosome (vs. M/P)
  state <- "I"
  
  #Scans cell from the first to last position of chromosome to detect parental origin (M or P)
  for (i in 1:length(v)) {
    
    #Scans for first coordinate of M block
    if (v[i] != "." & v[i] == "M" & v[i] != state) {
      state <- "M"
      #Stores first coordinate of M block
      bp <- c(bp, i)
      #Scans for previous P block end coordinate up until new M block
      bp <- c(bp, suppressWarnings(max(which(v[1:i] == "P"))))
      
      #Scans for first coordinate of new P block
    } else if (v[i] != "." & v[i] == "P" & v[i] != state) {
      state <- "P"
      #Stores the first coordinate of P block 
      bp <- c(bp, i)
      #Scans for previous M block end coordinate up until new P block
      bp <- c(bp, suppressWarnings(max(which(v[1:i] == "M"))))
    }
  } #end of going along every element of the vector with the parental states
  
  #Removes all NA/null values of bp
  bp <- bp[is.finite(bp)]
  
  #Breakpoint computation requires min. 2 values to compute coordinates of breakpoint between M and P
  if(length(bp) > 1) {
    
    #Sorts coordinates in ascending order
    bp <- sort(c(bp)) 
    
    #Breakpoint requires more than 1 parental state in cell 
  } else {
    bp <- 0
  }
  
  #List of chromosome positions for breakpoint coordinates
  bpCoord <- c()
  
  #Removes first coordinate of bp if length >2 because first coordinate is the beginning of an M or P block
  #Stores cooresponding bp coordinate as chromosome position
  if (length(bp) > 2) {
    bpCoord <- c(bpCoord, rownames(df)[bp][2:length(bp)]) 
    #if bp coordinate = 2, store cooresponding chromosome positions of both coordinates as it represents a true breakpoint
  } else if (length(bp) == 2) {
    bpCoord <- c(bpCoord, rownames(df)[bp])
  } else {
    bpCoord <- 0
  }
  
  #Returns the size of different breakpoints 
  #print(bp)
  return(as.numeric(bpCoord)) 
}
##----
#--


#Function to compute size of unknown blocks
unknown_chrom_position <- function(v) { #v - vector with the parental states within a cell; it takes the df as global variable.
  
  #Test function with one cell in vcf file
  #v <- df[,109]
  
  #uCoord stores first and last coordinates of unknown blocks
  uCoord <- c()
  #Neutral parental state of the chromosome
  state <- "I"
  
  #Scans cell from  first to last position of chromosome to detect unknown positions
  for (i in 1:length(v)) {
    
    #Scans for first coordinate of new U block
    if (v[i] != "." & v[i] == "U" & v[i] != state) {
      state <- "U"
      #Stores first cooridnate of U block 
      uCoord <- c(uCoord, i) 
      
      #Scans for new M block and finds previous end point of U block
    } else if (v[i] != "." & v[i] == "M" & v[i] != state) {
      state <- "M"
      #Stores U block end coordinate
      uCoord <- c(uCoord, suppressWarnings(max(which(v[1:i] == "U"))))
      
      #Scans for new P block and finds previous end point of U block
    } else if (v[i] != "." & v[i] == "P" & v[i] != state) {
      state <- "P"
      #Stores U block end cooridnate in uCoord
      uCoord <- c(uCoord, suppressWarnings(max(which(v[1:i] == "U"))))
      
      #If U is the last state of chromosome, finds end point of U block
    } else if (i == length(v) & state == "U") {
      uCoord <- c(uCoord, suppressWarnings(max(which(v == "U"))))
    }
    
  } #end of going along every element of the vector with the parental states
  
  
  #Removes all NA/null values of uCoord and duplicate values
  uCoord <- unique(uCoord[is.finite(uCoord)])
  uChrom <- c()
  
  #size_unknown computation requires min. 2 values to compute distance between first and last U coordinate
  if(length(uCoord) > 1) {
    
    #Sorts coordinates of uCoord in ascending order
    #uCoord <- sort(c(uCoord))
    #print(ordereduCoord)
    
    #Stores size of unknown blocks, computed by difference of chromosome positions (df[,2] = chromosome positions) from 
    #the first and last coordinates of U blocks  
    #uCoordBlocks <- c()
    #Compute difference between last and first point of U blocks
    #for (i in seq(1, (length(uCoord)-1), 2)) {
    #uCoordBlocks <- c(uCoordBlocks, (df[,2][uCoord[i+1]] - df[,2][uCoord[i]]))
    uChrom <- c(uChrom, rownames(df)[uCoord]) 
    
    #No size_unknown computation if there are no U blocks
  } else {
    
    uChrom <- 0
  }
  
  #Returns the size of different unknown blocks 
  #print(uChrom) 
  return(as.numeric(uChrom))
}
##---------
#-----


#Test individual cell
#cell <- list_breakpoint_coordinates[[5]]
#for (cell in list_breakpoint_coordinates) {

breakpoint_counts <- function(cell) { #cell - vector with the breakpoint coordinates within a cell; 
  #it takes the list_breakpoint_cooridnate as global variable
  #cell <- list_breakpoint_coordinates[[54]]
  #FUnction only applies to cells with breakpoints
  if (length(cell) > 1) {
    #Check each pair of breakpoint coordinates (first and last)
    for (i in seq(1, length(cell)-1, 2)) {
      #Determine bin of the first coordinate
      (max(which(bin_start < cell[i])))
      #Determine bin of the last coordinate
      (min(which(bin_end > cell [i+1])))
      #If first and last cooridnate of breakpoint are in the same bin, increment count by 1
      if ((max(which(bin_start < cell[i]))) ==  (min(which(bin_end > cell [i+1])))) {
        bin_counts[(max(which(bin_start < cell[i])))] <- bin_counts[(max(which(bin_start < cell[i])))] + 1
        #If the difference between bin_end and bin_start > 1
      } else if ( (min(which(bin_end > cell [i+1]))) - (max(which(bin_start < cell[i]))) > 1 ) {
        #Store the value of the difference between bin_end and bin_start
        x <- (min(which(bin_end > cell [i+1]))) - (max(which(bin_start < cell[i])))
        #If the difference is even, increment count of median bin
        if (x %% 2 == 0) {
          bin_counts[(max(which(bin_start < cell[i]))) + x/2] <- bin_counts[(max(which(bin_start < cell[i]))) + x/2] + 1
          #If the difference is odd, increment count of smaller ofmedian bin
        } else {
          bin_counts[(max(which(bin_start < cell[i]))) + (x-1)/2] <- bin_counts[(max(which(bin_start < cell[i]))) + x/2] + 1
        }
        #If the difference between bin_end and bin_start = 1
      } else {
        #If 50% of breakpoint is covered by first bin, increment count of first bin
        if ((cell[i+1]-cell[i])/2 + bin_start[(max(which(bin_start < cell[i])))] < bin_end[(max(which(bin_start < cell[i])))]) {
          bin_counts[(max(which(bin_start < cell[i])))] <-  bin_counts[(max(which(bin_start < cell[i])))] + 1
          #If 50% of breakpoint not covered by first bin, increment count of second bin
        } else {
          bin_counts[(min(which(bin_end > cell [i+1])))] <- bin_counts[(min(which(bin_end > cell [i+1])))] + 1
        }
      }
    }
  } else {bin_counts <- 0}
  
  
  return(bin_counts)
}

##---------
#-----

#list_tmp <- list_breakpoint_counts

#list_tmp <- list_breakpoint_counts[which(lengths(list_breakpoint_counts) > 1)]
#m_tmp <- matrix(unlist(list_tmp), ncol = length(bin_counts), byrow = T)
#apply(m_tmp, 2, sum)


#############################
###
### End of functions
###
#############################

##------------------------------
#############################
###
###         Main
###
#############################
#this is the first script required to analyse recombination maps. 
# it requires the output of HMM to compute recombination rates in windows of XMb size. 
#it output files with number of breakpoints per XMb windows.
# it requires a folder specified in "globalPath" to contain the HMM inferred states for the different individuals
# together with a directory containing all the HapMap and deCODE recombination maps. This folder name is specified
# in "hapMapDeCodeRecMapDir" 
#source("/home/ialves/Dropbox/singleCellProject/phasing_donors/wgs_SC_phasing/GT.fullMatrices.Examples/scripts_Method2/scripts_resolution/breakpoint_functions.R")
#globalPath <- "/home/ialves/Dropbox/singleCellProject/phasing_donors/wgs_SC_phasing/GT.fullMatrices.Examples/inferenceMatrices"

source("/.mounts/labs/awadallalab/scratch/ialves/sc_recombination/scripts/breakpoint_functions.R")
globalPath <- "/.mounts/labs/awadallalab/scratch/ialves/sc_recombination/inferenceMatrices"


hapMapDeCodeRecMapDir <- "recombinationMaps_HapMap_deCODE"
.libPaths( c("/.mounts/labs/awadallalab/private/flamaze/R_3.6.0_packages", .libPaths() ) )
#import library to execute functions in parallel 
require(parallel)
#library(randomcoloR) # not installed in the cluster
require(RColorBrewer)
no_cores <- 4
#below is a vector of 30 distinct colors got from: https://medialab.github.io/iwanthue/
distColourVector <- c("#aa485b","#57cc63","#c952be","#54a126","#7862cf","#91bf3c","#da3e7f","#4aaa5a","#a04990","#beb236",
                      "#637fc5","#dc9130","#4ab3d2","#d8572c","#5bcea8","#d6424d","#50a37d","#c88dd8","#498132","#e283a3",
                      "#2b7650","#96547c","#9eb366","#a44d2d","#646e19","#e3896e","#646f36","#d1a661","#9b6931","#8e7e2f")


chrNb <- as.numeric(args[1])
#chrNb <- 21

print(paste0("Computing recombination break points for chrm: ", chrNb))
chrTag <- paste0("chr", chrNb)
hapMapRecombFile <- paste0("genetic_map_GRCh37_", chrTag, ".txt")
deCodeRecombFile <- paste0("male.map.", chrTag, ".txt")
#open file centrome region across chrms
centromerePos <- read.table(paste0(globalPath, "/centromereOnlyReg.txt"), header = F)
centromerePos <- centromerePos[match(paste0("chr", 1:22), centromerePos$V2), c(2,3,4)]

if (chrNb < 8) {
  
  pathToFile_tmp  <- paste0(globalPath, "/", chrTag, "/")
  dir.create(pathToFile_tmp)
  lapply(list.files(path=paste0(globalPath, "/", chrTag, "p/"), 
                    pattern = "*.stateSeqHMM_clean.txt"), function(x) { file.copy(paste0(globalPath, "/", chrTag, "p/",x),
                                                                                  pathToFile_tmp, recursive = F, copy.mode = T) }) 
  lapply(list.files(path=paste0(globalPath, "/", chrTag, "q/"), 
                    pattern = "*.stateSeqHMM_clean.txt"), function(x) { file.copy(paste0(globalPath, "/", chrTag, "q/",x),
                                                                                  pathToFile_tmp, recursive = F, copy.mode = T) }) 
  inferredHap_tmp <- list.files(path=paste0(globalPath, "/", chrTag), 
                                pattern = paste0("*", chrTag, "p.stateSeqHMM_clean.txt"))
  setwd(pathToFile_tmp)
  cat(paste0("Sample", "\t", "StartPos","\t","EndPos"), file = paste0(pathToFile_tmp, "chr", chrNb, "_centromereRegions.txt"), sep = "\n")
  
  for (i in 1:length(inferredHap_tmp)) {
    #i <- 1
    indID <- unlist(strsplit(inferredHap_tmp[i], split = "\\."))[1]  
    headMat_tmp1 <- scan(inferredHap_tmp[i], what = character(), nlines = 1)
    
    if(file.exists(paste0(indID, ".", chrTag, "q.stateSeqHMM_clean.txt"))) { 
      mat_tmp2_name <- paste0(indID, ".", chrTag, "q.stateSeqHMM_clean.txt")
      headMat_tmp2 <- scan(mat_tmp2_name, what = character(), nlines = 1)
      length(headMat_tmp1) == length(headMat_tmp2)
      intersectCells <- intersect(headMat_tmp1, headMat_tmp2)
      if (sum(is.na(match(intersectCells, headMat_tmp2))) > 0) {
        print("Chromosome arms of chr: ", chrNb, " do not have info for the same cells")
        break;
      } else {
        m1 <- read.table(file = inferredHap_tmp[i], header = TRUE, stringsAsFactors = FALSE)
        m1 <- m1[,match(intersectCells, colnames(m1))]
        dim(m1)
        m2 <- read.table(file = mat_tmp2_name, header = TRUE, stringsAsFactors = FALSE)
        m2 <- m2[,match(intersectCells, colnames(m2))]
        dim(m2)
        
        full_m_name <- paste0(indID, ".", chrTag, ".stateSeqHMM_clean.txt")
        write.table(rbind(m1,m2), file=full_m_name, col.names = T, row.names = T, quote = F)
        file.copy(paste0(globalPath, "/", hapMapDeCodeRecMapDir, "/", hapMapRecombFile),
                  pathToFile_tmp, recursive = F, copy.mode = T)
        file.copy(paste0(globalPath, "/", hapMapDeCodeRecMapDir, "/", deCodeRecombFile),
                  pathToFile_tmp, recursive = F, copy.mode = T)
        print(paste0("Full matrix of inferred chromosomes printed for chr: ", chrNb, " and ind: ", indID))
        varPos <- as.numeric(rownames(rbind(m1,m2)))
        maxPosIndShortArm <- max(which(varPos < centromerePos[which(centromerePos$V2 == chrTag),2]))
        minPosIndLongArm <- min(which(varPos > centromerePos[which(centromerePos$V2 == chrTag),3]))
        cat(paste0(indID, "\t", varPos[maxPosIndShortArm],"\t", varPos[minPosIndLongArm]), file = paste0(pathToFile_tmp, "chr", chrNb, "_centromereRegions.txt"),sep = "\n", append = T)
            
      }
      rm(headMat_tmp1, headMat_tmp2, m1, m2, mat_tmp2_name)
    }
}
  lapply(list.files(path=pathToFile_tmp, 
                    pattern = "*p.stateSeqHMM_clean.txt"), function(x) { file.remove(x) }) 
  lapply(list.files(path=pathToFile_tmp, 
                    pattern = "*q.stateSeqHMM_clean.txt"), function(x) { file.remove(x) }) 

} else {
  pathToFile_tmp  <- paste0(globalPath, "/", chrTag, "/")
  setwd(pathToFile_tmp)
  cat(paste0("Sample", "\t", "StartPos","\t","EndPos"), file = paste0(pathToFile_tmp, "chr", chrNb, "_centromereRegions.txt"), sep = "\n")
  inferredHap_tmp <- list.files(path=paste0(globalPath, "/", chrTag), 
                                pattern = paste0("*", chrTag, ".stateSeqHMM_clean.txt"))
  for (i in 1:length(inferredHap_tmp)) {
    indID <- unlist(strsplit(inferredHap_tmp[i], split = "\\."))[1]  
    #i <- 1
    m1 <- read.table(file = inferredHap_tmp[i], header = TRUE, stringsAsFactors = FALSE)
    varPos <- as.numeric(rownames(m1))
    indPosBeforeCentromere <- which(varPos < centromerePos[which(centromerePos$V2 == chrTag),2])
    if(length(indPosBeforeCentromere) > 0) {
      maxPosIndShortArm <- max(which(varPos < centromerePos[which(centromerePos$V2 == chrTag),2]))
      minPosIndLongArm <- min(which(varPos > centromerePos[which(centromerePos$V2 == chrTag),3])) 
      cat(paste0(indID, "\t", varPos[maxPosIndShortArm],"\t", varPos[minPosIndLongArm]), file = paste0(pathToFile_tmp, "chr", chrNb, "_centromereRegions.txt"),sep = "\n", append = T)
      
    } else {
      #maxPosIndShortArm <- 1
      minPosIndLongArm <- min(which(varPos > centromerePos[which(centromerePos$V2 == chrTag),3])) 
      cat(paste0(indID, "\t", 1,"\t", varPos[minPosIndLongArm]), file = paste0(pathToFile_tmp, "chr", chrNb, "_centromereRegions.txt"),sep = "\n", append = T)
    }
    rm(m1,varPos,indID)
  }
}


pathToFile <- paste0(globalPath, "/", chrTag, "/resolution")
dir.create(pathToFile)
lapply(list.files(path=paste0(globalPath, "/", chrTag), 
                  pattern = "*.stateSeqHMM_clean.txt"), function(x) { file.copy(paste0(globalPath, "/", chrTag, "/", x),
                                                                                pathToFile, recursive = F, copy.mode = T) }) 
lapply(list.files(path=paste0(globalPath, "/", chrTag), 
                  pattern = "*.stateSeqHMM_clean.txt"), function(x) { file.remove(x) }) 
file.copy(paste0(globalPath, "/", hapMapDeCodeRecMapDir, "/", hapMapRecombFile),
          pathToFile, recursive = F, copy.mode = T)
file.copy(paste0(globalPath, "/", hapMapDeCodeRecMapDir, "/", deCodeRecombFile),
          pathToFile, recursive = F, copy.mode = T)

setwd(pathToFile)

####---------------
###-------------
#Section applies for finding breakpoints/bin 

chrSizes <- c(249250621, 243199373,198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 
              141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 
              81195210, 78077248, 59128983, 63025520, 48129895, 51304566)

#Chromosome 15
chrmLength <- chrSizes[chrNb]
print(chrmLength)
windowSize <- 1000000 #1Mb
bin_start<- c(seq(1, chrmLength, by=windowSize))
bin_end <- c(seq(windowSize, chrmLength, by=windowSize), chrmLength)

bin_counts <- vector(mode = "numeric" , length = length(bin_start))
nb_of_breakpoints <- list()
medSizeBreakpoint <- list()
bp_size_list_clean <- list()
nb_of_cells <- list()
q.lim <- .99 #no longer used. Based now on the 3*sd(mean)

#ADD commands defining the data frame file name
# open the data frame and save is as df
dataframes <- list.files(pattern = paste0("*", chrTag, ".stateSeqHMM_clean.txt$"))

cat(paste0("Ind", "\t", "MeanNbOfBps", "\t", "NbOfCells", "\t", "MedBpSize(kb)", "\n"), 
        file=paste0("summaryBp_", chrTag, ".txt"), sep = "")

for (file in 1:length(dataframes)) {
  
  #Takes vcf file specifying individual and chromosome # to detect breakpoints within single cells 
  #file <- 11
  indID <- unlist(strsplit(dataframes[file], split = "\\."))[1]  
  df <- read.table(file = dataframes[file], header = TRUE, stringsAsFactors = FALSE)
  
  #List of breakpoints, unknown, M, and P blocks
  list_breakpoint_coordinates <- list()
  list_unknown_coordinates <- list()
  list_breakpoint_counts <- list()
  #List of breakpoints, unknown, M, and P blocks
  list_breakpoint_size <- list()
  list_unknown_size <- list()
  list_M_size <- list()
  list_P_size <- list()
  
  #parLapply -> speeds up functions by splitting the computation in paralell 
  # no_cores <- detectCores()-1
  # cl <- makeCluster(no_cores)
  # clusterExport(cl, c("df", "breakpoint_chrom_position", "unknown_chrom_position"))
  list_breakpoint_coordinates <- mclapply(1:ncol(df), function(x) { breakpoint_chrom_position(df[,x]) }, mc.cores = no_cores)
  list_unknown_coordinates <- mclapply(1:ncol(df), function(x) { unknown_chrom_position(df[,x]) }, mc.cores = no_cores)
  #stopCluster(cl)
  
  listCellsWithRecomb <- list_breakpoint_coordinates[lengths(list_breakpoint_coordinates) > 1]
  length(list_breakpoint_coordinates) 
  names(listCellsWithRecomb) <- colnames(df)[which(lengths(list_breakpoint_coordinates) > 1)]
  cat(paste0("Cell_ID", "\t","BreakpointCoordinates"), file=sub(".stateSeqHMM_clean.txt", "_bpCoordinates_perCell.txt", dataframes[file]), sep = "\n")
  for(cell in 1:length(listCellsWithRecomb)) {
    #cell <- 1 
    cat(paste0(names(listCellsWithRecomb)[cell], "\t", paste(listCellsWithRecomb[[cell]], collapse = "\t")),
        file=sub(".stateSeqHMM_clean.txt", "_bpCoordinates_perCell.txt", dataframes[file]), sep = "\n", append = T)
  }
  #}
  #parLapply -> speeds up functions by splitting the computation in paralell 
  # no_cores <- detectCores()-1
  # cl <- makeCluster(no_cores)
  # clusterExport(cl, c("df", "breakpoint_computation", "size_unknown", "size_M", "size_P"))
  list_breakpoint_size <- mclapply(1:ncol(df), function(x) { breakpoint_computation(df[,x]) }, mc.cores = no_cores)
  length(list_breakpoint_size)
  list_unknown_size <- mclapply(1:ncol(df), function(x) { size_unknown(df[,x]) }, mc.cores = no_cores)
  list_M_size <- mclapply(1:ncol(df), function(x) { size_M(df[,x]) }, mc.cores = no_cores)
  list_P_size <- mclapply(1:ncol(df), function(x) { size_P(df[,x]) }, mc.cores = no_cores)
  #stopCluster(cl)

  nb_of_breakpoints[[file]] <- lengths(list_breakpoint_size)
  #extremeNb_bp_cells <- which(nb_of_breakpoints[[file]] > mean(nb_of_breakpoints[[file]])+3*sd(nb_of_breakpoints[[file]]))
  extremeNb_bp_cells <- which(nb_of_breakpoints[[file]] > 5)
  
  #removing cells with too many bp from the list of BP coord and unkown coord
  #added by Feb 4 2020
  if(length(extremeNb_bp_cells) > 0) {  
    list_breakpoint_coordinates <- list_breakpoint_coordinates[-c(extremeNb_bp_cells)]
    list_unknown_coordinates <- list_unknown_coordinates[-c(extremeNb_bp_cells)] #added to remove cells with too many bp
    bp_size_list_clean[[file]] <- unlist(list_breakpoint_size[-c(extremeNb_bp_cells)])
    nb_of_cells[[file]] <- length(list_breakpoint_size[-c(extremeNb_bp_cells)])
    nb_of_breakpoints[[file]] <- length(bp_size_list_clean[[file]])
    
  } else {
    bp_size_list_clean[[file]] <- unlist(list_breakpoint_size)
    nb_of_cells[[file]] <- length(list_breakpoint_size)
    
  }
  ##-- end of the part added on Feb 4 2020
  
  list_breakpoint_counts <- lapply(1:length(list_breakpoint_coordinates), function(x) { breakpoint_counts(list_breakpoint_coordinates[[x]]) }) #changed to remove cells with too many bp
  #lapply(1:54, function(x) { breakpoint_counts(list_breakpoint_coordinates[[x]]) }) #changed to remove cells with too many bp
  #Sum list_breakpoint_counts
  bp <- Reduce('+', list_breakpoint_counts[which(lengths(list_breakpoint_counts) > 1)])
  
  #Export the list_breakpoint_counts to corresponding file
  write.table(bp, file = paste0(sub(".txt", "_bp_counts.txt", dataframes[file])), row.names = FALSE, col.names = FALSE)
  
  
  cell_counts <- vector(mode = "numeric", length = length(bin_start))
  for (i in 1:length(cell_counts)) { #setting up the nb of informative cells to the number of surveyed cells after filtering. 
    cell_counts[i] <- length(list_unknown_coordinates)}
  
  #Export the cell_counts to corresponding file
  cc <- cell_counts
  write.table(cc, file = paste0(sub(".txt", "_cell_counts.txt", dataframes[file])), row.names = FALSE, col.names = FALSE)
  
  #RESOLUTION plot
  pdf(file = paste0(sub("stateSeqHMM_clean.txt", "resolutionWGS.95.pdf", dataframes[file])), height = 5, width = 8)
  par(mfrow=c(1,1), mar=c(5,4,2,4))
  medSizeBreakpoint_clean <- median(bp_size_list_clean[[file]])/1000
  histDist <- hist(bp_size_list_clean[[file]][bp_size_list_clean[[file]] > 0]/1000, 
                   breaks = c(seq(0, max(bp_size_list_clean[[file]])/1000, by = 100), max(bp_size_list_clean[[file]])+10/1000), plot = F)
  cdist <- ecdf(bp_size_list_clean[[file]][bp_size_list_clean[[file]] > 0]/1000)
  barplot(histDist$counts,  
          xlab="Size of interval (kb)",
          ylab = "Number of breakpoints", xaxt='n', space = 0, main = paste0(indID, "_", chrTag))
  lines(x = c(0:length(histDist$counts)), y=c(0,cdist(histDist$mids)*max(histDist$counts)), col ='red')
  axis(4, at=seq(from = 0, to = max(histDist$counts), length.out = 11), labels=seq(0, 1, 0.1), col = 'red', col.axis = 'red')
  axis(1, at=seq(0,length(histDist$counts), length.out = length(c(seq(0, max(bp_size_list_clean[[file]])/1000, by = 500)))+1), 
       labels=c(seq(0, max(bp_size_list_clean[[file]])/1000, by = 500), (max(bp_size_list_clean[[file]])+10)/1000), col = 'black', col.axis = 'black')
  mtext(side = 4, line = 3, 'Cumulative Density', col = 'red')
  legend(length(histDist$counts)-length(histDist$counts)/2, max(histDist$counts), legend=paste0("Mean nb of breakpoints: ", round(sum(nb_of_breakpoints[[file]])/nb_of_cells[[file]], digits = 2)), bty = "n")  
  legend(length(histDist$counts)-length(histDist$counts)/2, max(histDist$counts)-(max(histDist$counts)/6), legend=paste0("Nb of cells: ", nb_of_cells[[file]]), bty = "n")  
  legend(length(histDist$counts)-length(histDist$counts)/2, max(histDist$counts)-(max(histDist$counts)/4), legend=paste0("Median size kb: ", round(medSizeBreakpoint_clean, digits = 2)), bty = "n")  
  dev.off()
  
  #printing in a file summary statistics
  cat(paste0(indID, "\t", round(sum(nb_of_breakpoints[[file]])/nb_of_cells[[file]], digits = 2),
	                "\t", nb_of_cells[[file]], "\t", round(medSizeBreakpoint_clean, digits = 2), "\n"), 
          file=paste0("summaryBp_", chrTag, ".txt"), sep = "", append = TRUE) 

} 
####################
##
## Plotting recombination rates along the chroms
##
####################

#computing recombination rates from deCode
deCode_chr <- data.frame(read.table(file=deCodeRecombFile, header = F, na.strings = "NA"))
deCodecM <- c()

for (binbNb in 1:length(bin_counts)) {
  
  #binNb <- 1
  index_tmp <- which(deCode_chr[,3] >= bin_start[binbNb] & deCode_chr[,3] < bin_end[binbNb])
  
  if (length(index_tmp) > 0) {
    if (sum(is.na(deCode_chr[index_tmp, 4])) != length(deCode_chr[index_tmp, 4])) {
      deCodecM[binbNb] <- sum(deCode_chr[index_tmp, 4][!is.na(deCode_chr[index_tmp, 4])])
    } else {
      deCodecM[binbNb] <- 0
    }
  } else {
    deCodecM[binbNb] <- 0
  }
}
##-----------

#computing HapMap recombination rates
hapMap_chr <- data.frame(read.table(file=hapMapRecombFile, header = T, na.strings = "NA"))
hapMapcM <- c()
startRecomb <- 0

for (binNb in 1:length(bin_counts)) {
  
  #binNb <- 1
  index_tmp <- which(hapMap_chr[,2] >= bin_start[binNb] & hapMap_chr[,2] < bin_end[binNb])
  if (length(index_tmp) > 0) {
    
    if (binNb == 1) {
      hapMapcM[binNb] <- max(hapMap_chr[index_tmp,4])  
    } else {
      hapMapcM[binNb] <- max(hapMap_chr[index_tmp,4])-startRecomb
    }
    
    startRecomb <- max(hapMap_chr[index_tmp,4])
  } else {
    
    hapMapcM[binNb] <- 0
  }
  
}

####-----------------
# computing recombination rates from single cell sperm cells
empiricalRecombRates <- list()
QC_tag <- list()
filesToAnalyse <- list.files(pattern = "*bp_counts.txt$")
openIndCenLocation  <- read.table(paste0(pathToFile_tmp,  chrTag,"_centromereRegions.txt"), header = T)
for (files in 1:length(filesToAnalyse)) {
  #files <- 2
  indTag <- unlist(strsplit(filesToAnalyse[files], split="\\."))[1]
  QC_tag[[files]] <- strsplit(filesToAnalyse[files], split = "\\.")[[1]][1]
  empiricalCountsF <- scan(file=filesToAnalyse[files], what = numeric())
  empiricalCellCountsF <- scan(file=sub("bp", "cell", filesToAnalyse[files]), what = numeric())
  
  cenLocation <- openIndCenLocation[which(openIndCenLocation[,1] == indTag),2:3]
  binCentromere <- c(max(which(bin_start <= as.numeric(as.vector(cenLocation$StartPos)))), min(which(bin_end > as.numeric(as.vector(cenLocation$EndPos)))))
  empiricalCountsF[binCentromere[1]:binCentromere[2]] <- 0  
  empiricalRecombRates[[indTag]] <- (empiricalCountsF/empiricalCellCountsF)*100
  
}
##------------
colors <- sample(distColourVector,length(empiricalRecombRates))

#########################
##
## computing the cumulative proportions
##
#########################
pdf(file=paste0("genetic_vs_physicalDist_", chrTag, ".pdf"), width=6, height=6)
par(mar=c(5,4,2,2))
plot(1:length(hapMapcM), y=cumsum(hapMapcM), type="l", xaxt='n', xlim=c(0,max(cumsum(hapMapcM))/(windowSize/1e6)), ylim=c(0,max(max(cumsum(hapMapcM)), max(unlist(lapply(empiricalRecombRates, cumsum))))), xlab="Physical Distance (Mb)", 
     ylab="Genetic Distance (cM)", col="grey80", lwd=2, lty=4)
axis(1, at = seq(0, max(cumsum(hapMapcM))/(windowSize/1e6), by = 20/(windowSize/1e6)), labels = seq(0, max(cumsum(hapMapcM)), by =20))

lines(1:length(hapMapcM), y=cumsum(deCodecM), type="l", xlim=c(0,max(cumsum(deCodecM))), ylim=c(0,max(cumsum(deCodecM))), xlab="Physical Distance (Mb)", 
      ylab="Genetic Distance (cM)", col="black",lwd=2, lty=3)
rect((centromerePos[which(centromerePos[,1] == chrTag),2]/1e6)/(windowSize/1e6), 0, (centromerePos[which(centromerePos[,1] == chrTag),3]/1e6)/(windowSize/1e6), max(max(cumsum(hapMapcM)), max(unlist(lapply(empiricalRecombRates, cumsum)))),
     col = "grey88", border = NA)

for (files in 1:length(empiricalRecombRates)) {
  
  lines(1:length(hapMapcM), y=cumsum(empiricalRecombRates[[files]]), type="l", xlim=c(0,max(cumsum(deCodecM))), ylim=c(0,max(cumsum(deCodecM))), xlab="Physical Distance (Mb)", 
        ylab="Genetic Distance (cM)", col=colors[files],lwd=2)
  
}
legend("topright", legend = c("deCODE", "HapMap", unlist(QC_tag)), lty = c(3,4,rep(1,length(empiricalRecombRates))), lwd = 2, col = c("black", "grey80", colors), bty = "n", cex = 0.60)
legend("topleft", legend = chrTag, bty = "n", cex = 1.5)
dev.off()
##------------
#--------
pdf(file = paste0("recombRates_SNPDensity_chr", chrTag, "_2.pdf"), width = 11, height = 4)
par(mfrow=c(1,1), mar=c(5,4,2,2))
##----------
# set up the plot 
y_range <- range(unlist(empiricalRecombRates))
x_range <- range(-6:length(bin_counts))
tickMarks <- seq(1,length(bin_counts), by = 20)

plot(x_range, y_range, type="n", xlab=paste0("Chromosome position (" , chrTag,")"),
     ylab="cM/Mb", xaxt='n')
axis(1, at = tickMarks, labels = format(apply(rbind(bin_start, bin_end), 2, mean), scientific = T, digits = 2)[tickMarks])
legend("topleft", legend = c("deCODE", "HapMap", unlist(QC_tag)), lty = c(3,4,rep(1,length(empiricalRecombRates))), lwd = 2, col = c("black", "grey80", colors), bty = "n", cex = 0.60)
rect((centromerePos[which(centromerePos[,1] == chrTag),2]/1e6)/(windowSize/1e6), 0, (centromerePos[which(centromerePos[,1] == chrTag),3]/1e6)/(windowSize/1e6), y_range[2],
     col = "grey88", border = NA)

for (files in 1:length(empiricalRecombRates)) {
  
  lines(unlist(empiricalRecombRates[files]), type = "l", col = colors[files], lwd=1.5)  
  
  
}
lines(deCodecM, type = "l", lty= 3, lwd = 2, col = "black")
lines(hapMapcM, type = "l", lty= 4, lwd = 2, col = "grey80")
dev.off()
##-------------

##########################
##
##
## END
##
##
##########################
