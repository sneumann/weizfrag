## devtools::install_github("rformassspectrometry/MsBackendMsp")
## devtools::install_github("ipb-halle/MetFragRelaunched/MetFragR/rpackage/MetFragR")

library(Spectra)
library(MsBackendMsp)

options(java.parameters = c("-Xms512m", "-Xmx1024m"))
library(rcdk)
library(metfRag)

#fls <- dir(system.file("extdata", package = "MsBackendMsp"),
#           full.names = TRUE, pattern = "msp$")
#fls <- "/home/sneumann/src/MetFragRelaunched/MetFragLib/src/main/resources/MoNA-export-LC-MS.mb"
fls <- "/tmp/MoNA-export-LC-MS-MS_Spectra-20241014.msp"
#fls <- "/media/MoNA-export-LC-MS-MS_Spectra.msp"
#fls <- "/media/MoNA-export-LC-MS-MS_Spectra-half.msp"
#fls <- "/media/MoNA-export-LC-MS-MS_Spectra-200.mb"
#fls <- "/tmp/head-MoNA-export-LC-MS-MS_Spectra.msp"

if (TRUE) {
  ## Import a single file.
  # register(SerialParam())
  register(MulticoreParam(4))
  sps <- Spectra(fls[1], source = MsBackendMsp())
  save(sps, file="/tmp/sps.rdata")
} else {
  #load("/tmp/sps.rdata")
  load("/tmp/MoNA-export-LC-MS-MS_Spectra-20221217.rdata")
}

## Filter all Peaks relative intensity < 0.005
## https://rformassspectrometry.github.io/Spectra/articles/Spectra.html#data-manipulations
## and help("applyProcessing")
##
keep_peaks_relative <- function(x, thresh_relative = 0.05) {
  x >= max(x, na.rm = TRUE) * thresh_relative
}
sps <- filterIntensity(sps, intensity = keep_peaks_relative, thresh_relative = 0.005)

## Filter all apectra with >100 peaks after intensity filtering
##
np <- as.numeric(spectraData(sps, "Num.Peaks")[,1])
sps <- sps[np<1000]

inchikey <- spectraData(sps, columns = "InChIKey")[,1]
ion_mode <- spectraData(sps, columns = "Ion_mode")[,1]
comments <- spectraData(sps, columns = "Comments")[,1]

## Filter any spectra without InChIKey, ion_mode or comments
hasNA <- is.na(inchikey) | is.na(ion_mode) | is.na(comments)
sps <- sps[!hasNA]

inchikey <- spectraData(sps, columns = "InChIKey")[,1]
ion_mode <- spectraData(sps, columns = "Ion_mode")[,1]
comments <- spectraData(sps, columns = "Comments")[,1]
 
computed <- as.data.frame(t(sapply(comments, function(com) {
  entries <- scan(text=com, what="character", quiet = TRUE)
  compsmiles <- sub("^computed SMILES=(.*)$", "\\1", grep("^computed SMILES=", entries, value=TRUE))
  if (length(compsmiles) != 1) 
    compsmiles <- NA
  compinchi <- sub("^InChI=(.*)$", "\\1", grep("^InChI=", entries, value=TRUE))
  if (length(compinchi) != 1) 
    compinchi <- NA
  return(cbind(compsmiles=compsmiles, compinchi=compinchi))
}, USE.NAMES = FALSE)), stringsAsFactors = FALSE)
colnames(computed) <- c("compsmiles", "compinchi")

stopifnot(   length(inchikey)   == length(ion_mode) 
             && length(ion_mode)   == length(comments) 
             && length(comments)   == nrow(computed) )

rm (comments)
gc()

sp <- get.smiles.parser()

fps <- sapply(computed[,"compsmiles"], function(s) {
  m <- parse.smiles(s, smiles.parser = sp)[[1]]
  fp <- ""
  if (!is.null(m) && class(m) == "jobjRef" && m@jclass == "org/openscience/cdk/interfaces/IAtomContainer") {
    fp <- calculateFingerprintFromSmiles(s)
  } else {
    cat("Broken SMILES: ", s, "\n")
    fp <- ""
  }
})

smilesOK <- sapply (fps, function(f) f!="")
table(smilesOK)

spd <- spectraData(sps, c("Name", "PrecursorMZ"))
hasPrecursor <- !is.na(spd[, "PrecursorMZ"])
table(hasPrecursor)

## Fix comma/dot issue in PrecursorMZ
spd[, "PrecursorMZ"] <- gsub(",", ".", spd[, "PrecursorMZ"])

peakLists <- peaksData(sps)


stopifnot(   length(smilesOK)   == length(hasPrecursor) 
             && length(smilesOK)   == length(sps)
             && length(smilesOK)   == length(peakLists)
             && length(smilesOK)   == nrow(spd) )

len <-  length(which(smilesOK==TRUE & hasPrecursor==TRUE))
               


sink("/tmp/MoNA-export-LC-MS-MS_Spectra-20241014-0.005.mb")
for (i in seq(1, length(sps))[smilesOK & hasPrecursor ]) {
  if ( (i/len*100)%% 10 == 0) {
    message(i/len*100, "%")
  }
  spec <- sps[i]
  
  cat(
    "# SampleName = ", spd[i,"Name"], "\n",
    "# InChI = ", computed[i,"compinchi"], "\n",
    "# InChIKey = ", inchikey[i], "\n",
    "# IsPositiveIonMode = ", ifelse(ion_mode[i]=="P", "True", "False"), "\n",
    "# IonizedPrecursorMass = ", spd[i, "PrecursorMZ"], "\n",
    "# NumPeaks = ", nrow(peakLists[[i]]), "\n",
    "# MolecularFingerPrint = ", fps[i], "\n", sep="")
  
  write.table(peakLists[[i]], row.names = F, col.names = F)
  cat("\n")
}
sink(NULL)

## Create CSV with InChiKey and Fingerprint
## egrep 'InChIKey|MolecularFingerPrint' MoNA-export-LC-MS-MS_Spectra.mb | cut -d"=" -f 2 | paste -s -d' \n' | sort | uniq  >new.mb.fp 

## Join two of such CSVs to detect (mis)matches
## join new.mb.fp old.mb.fp | while read I F1 F2 ; do if [ "$F1" == "$F2" ] ; then echo "OK" ; else echo "mismatch: $I" ; echo $F1 ; echo $F2  ; fi ; done 




