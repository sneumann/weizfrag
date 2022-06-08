## devtools::install_github("rformassspectrometry/MsBackendMsp")
## devtools::install_github("sneumann/MsBackendMS")
library(Spectra)
library(MsCoreUtils)
library(MsBackendMS)
library(MsBackendMsp)
library(metfRag)
library(IRanges)

setwd("/home/sneumann/src/weizfrag")
be <- MsBackendMS()

## Read the CHemistry Metadata
metadata <- read.delim("msExports/WeizMassV2.tsv", stringsAsFactors=F)
ion_mode <- metadata[, "ionization"]
cpdid    <- metadata[, "compound_id"]
filename <- paste("msExports/", cpdid, ".", tolower(ion_mode), ".ms", sep="")

## Remove the ones where I deleted the *.ms lacking MS2 
fileexists <- file.exists(filename)
metadata <- metadata[fileexists, ]
ion_mode <- metadata[, "ionization"]
cpdid    <- metadata[, "compound_id"]
filename <- paste("msExports/", cpdid, ".", tolower(ion_mode), ".ms", sep="")

## Extract metadata
inchi <- metadata[, "InChI"]
inchikey <- metadata[, "InChiKey"]
smiles <- metadata[, "Smiles"]

## Read actual SIRIUS input files
ms <- backendInitialize(be, filename)

## Calculate Fingerprints for MetFrag
options(java.parameters = c("-Xms512m", "-Xmx1024m"))
library(rcdk)
sp <- get.smiles.parser()

fps <- sapply(smiles, function(s) {
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

inchikeyOK <- inchikey!=""

sink("/tmp/weizfragV2.mb")
for (i in seq(1, length(ms))[smilesOK & inchikeyOK]) {
  if ( (i/len*100)%% 10 == 0) {
    message(i/len*100, "%")
  }
  spec <- ms[i]

  cat(
    "# SampleName = ", spectraData(spec)[,"compound"], "\n",
    "# InChI = ", inchi[i], "\n",
    "# InChIKey = ", inchikey[i], "\n",
    "# IsPositiveIonMode = ", ifelse(ion_mode[i]=="POSITIVE", "True", "False"), "\n",
    "# IonizedPrecursorMass = ", spectraData(spec)[,"parentmass"], "\n",
    "# NumPeaks = ", nrow(peaksData(spec)[[1]]), "\n",
    "# MolecularFingerPrint = ", fps[i], "\n", sep="")
   
    write.table(peaksData(spec), row.names = F, col.names = F)
    cat("\n")
}
sink(NULL)

## Try to write *.MSP files
export(Spectra(ms), backend = MsBackendMsp(), file = "/tmp/WeizMassV2.msp")

## Create CSV with InChiKey and Fingerprint
## egrep 'InChIKey|MolecularFingerPrint' MoNA-export-LC-MS-MS_Spectra.mb | cut -d"=" -f 2 | paste -s -d' \n' | sort | uniq  >new.mb.fp 

## Join two of such CSVs to detect (mis)matches
## join new.mb.fp old.mb.fp | while read I F1 F2 ; do if [ "$F1" == "$F2" ] ; then echo "OK" ; else echo "mismatch: $I" ; echo $F1 ; echo $F2  ; fi ; done 




