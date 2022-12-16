# Creating spectral libraries for MetFrag

## MoNA

First get a recent MoNA download in MSP format from https://mona.fiehnlab.ucdavis.edu/downloads
https://mona.fiehnlab.ucdavis.edu/rest/downloads/retrieve/8e7a9605-728e-4758-b090-477ff1e5be73

The `convert-mona.R` is an example to read the MSP file,
apply a small bit of filtering that improves runtime afterwards,
and calculates the Fingerprints with RCDK,
and finally outputs the library as an `*.mb`
