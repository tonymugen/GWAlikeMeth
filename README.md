This is an R package to perform genome-wide association studies on replicated data. It implements the [EMMAX](http://genetics.cs.ucla.edu/emmax/) method, including cases where each genotype has more than one observation. In addition, if there are several traits it pays to analyze them together. While the GWA itself is done separately for each trait, some components of the mixed model and the SNP regression are common among traits and thus some time savings can be achieved by considering traits from the same data set together.

The SNP regression portion is multi-threaded. This is implemented so as not to interfere with multi-threading of linear algebra libraries used in R. More details can be found in the documentation.

To install, make sure you have the `devtools` package on your system, and then run `install_github("tonymugen/GWAlikeMeth")`. It should work under any Unix-like system, but has not been tested under Windows.

