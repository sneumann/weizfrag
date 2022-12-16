## devtools::install_github("rformassspectrometry/MsBackendMsp")
## devtools::install_github("rformassspectrometry/MsBackendMassBank")
## devtools::install_github("ipb-halle/MetFragRelaunched/MetFragR/rpackage/MetFragR")

library(Spectra)
library(MsBackendMassbank)

options(java.parameters = c("-Xms512m", "-Xmx1024m"))
library(rcdk)
library(metfRag)

fls <- list.files("/vol/massbank/data/MassBank-data/",
           full.names = TRUE, pattern = "txt$", recursive = TRUE)

## During development: use less compounds
subsample <- sample(seq(1,length(fls)), size=0.1*length(fls))
fls <- fls[subsample]

# create data frame to indicate with metadata blocks shall be read.
metaDataBlocks <- data.frame(metadata = c("ac", "ch", "sp", "ms",
                                          "record", "pk", "comment"),
                             read = rep(TRUE, 7))

if (TRUE) {
  ## Import a single file.
  sps <- Spectra(fls,
                 source = MsBackendMassbank(),
                 backeend = MsBackendDataFrame(),
                 metaBlock = metaDataBlocks)
  save(sps, file="/tmp/sps.rdata")
} else {
  load("/tmp/sps.rdata")
}

## Filter all Peaks relative intensity < 0.005
## https://rformassspectrometry.github.io/Spectra/articles/Spectra.html#data-manipulations
## and help("applyProcessing")
##
keep_peaks_relative <- function(x, thresh_relative = 0.05) {
  x >= max(x, na.rm = TRUE) * thresh_relative
}
sps <- filterIntensity(sps, intensity = keep_peaks_relative, thresh_relative = 0.005)

# all spectraVariables possible in MassBank are read
spectraVariables(sps)

## Filter all apectra with >100 peaks after intensity filtering
##
np <- as.numeric(spectraData(sps, "pknum")[,1])
sps <- sps[np<1000]

inchikey <- spectraData(sps, columns = "inchikey")[,1]
ion_mode <- sps$polarity

# Sanity check
stopifnot(   length(inchikey)   == length(ion_mode))

gc()

sp <- get.smiles.parser()

fps <- sapply(sps$smiles, function(s) {
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
len <-  length(which((smilesOK==TRUE)))
table(smilesOK)

spd <- spectraData(sps, c("accession", "precursorMz", "inchi", "inchikey"))

peakLists <- peaksData(sps)

sink("/tmp/MoNA-export-LC-MS-MS_Spectra-0.005.mb")
for (i in seq(1, length(sps))[smilesOK & !is.na(spd[, "precursorMz"]) ]) {
  if ( (i/len*100)%% 10 == 0) {
    message(i/len*100, "%")
  }

  spec <- sps[i]
  
  cat(
    "# SampleName = ", spd[i,"accession"], "\n",
    "# InChI = ", spd[i,"inchi"], "\n",
    "# InChIKey = ", spd[i, "inchikey"], "\n",
    ## MetFrag wants ion mode TRUE=positive, FALSE=negative
    "# IsPositiveIonMode = ", ifelse(ion_mode[i]==1, "True", "False"), "\n",
    "# IonizedPrecursorMass = ", spd[i, "precursorMz"], "\n",
    "# NumPeaks = ", nrow(peakLists[[i]]), "\n",
    "# MolecularFingerPrint = ", fps[i], "\n", sep="")
  
  write.table(peakLists[[i]], row.names = F, col.names = F)
  cat("\n")
}
sink(NULL)



