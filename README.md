# polyG

## Installation



### Install prerequsite R packages
```r
# Install missing R packages from CRAN:
required.R.packages <- c('here','gplots','RColorBrewer','colorspace')
new.R.packages <- required.R.packages[!(required.R.packages %in% installed.packages()[,"Package"])]
if(length(new.R.packages)>0) install.packages(new.R.packages)

# Install missing R packages from Bioconductor:
required.BC.packages <- c('phyloseq','dendextend','ape','phangorn','adephylo')
new.BC.packages <- required.BC.packages[!(required.BC.packages %in% installed.packages()[,"Package"])]
if(length(new.BC.packages)>0) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(new.BC.packages)
}
```

### Install the polyG package:
```r
# Clone this GitHub repo. Then in R, run:
install.packages('[PATH TO CLONED REPO]/polyG/pkg', type='src', repos=NULL)
```

## Running poly-G pipeline
