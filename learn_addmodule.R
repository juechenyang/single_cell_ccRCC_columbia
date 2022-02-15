# Function arguments
DefaultAssay(test) = "RNA"
object = test
features = hla_markers
pool = rownames(object)
nbin = 24
ctrl = 100
k = FALSE
name = "hla_enriched"
seed = 1


# Find how many gene lists were provided. In this case just one.
cluster.length <- length(x = features)

# Pull the expression data from the provided Seurat object
assay.data <- GetAssayData(object = object)
# For all genes, get the average expression across all cells (named vector)
data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
# Order genes from lowest average expression to highest average expression
data.avg <- data.avg[order(data.avg)]

# Use ggplot2's cut_number function to make n groups with (approximately) equal numbers of observations. 
#The 'rnorm(n = length(data.avg))/1e+30' part adds a tiny bit of noise to the data, presumably to break ties.
data.cut <- ggplot2::cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30,
                                n = nbin,
                                labels = FALSE,
                                right = FALSE)

# Set the names of the cuts as the gene names
names(x = data.cut) <- names(x = data.avg)

# Create an empty list the same length as the number of input gene sets. This will contain the names of the control genes
ctrl.use <- vector(mode = "list", length = cluster.length)

# For each of the input gene lists:
for (i in 1:cluster.length) {
  # Get the gene names from the input gene set as a character vector  
  features.use <- features[[i]]
  
  # Loop through the provided genes (1:num_genes) and for each gene, find ctrl (default=100) genes from the same expression bin (by looking in data.cut):
  for (j in 1:length(x = features.use)) {
    # Within this loop, 'data.cut[features.use[j]]' gives us the expression bin number. We then sample `ctrl` genes from that bin without replacement and add the gene names to ctrl.use.
    ctrl.use[[i]] <- c(ctrl.use[[i]],
                       names(x = sample(x = data.cut[which(x = data.cut == data.cut[features.use[j]])],
                                        size = ctrl,
                                        replace = FALSE)))
  }
}

# Have a quick look at what's in ctrl.use:
class(ctrl.use)
## [1] "list"
length(ctrl.use)
## [1] 1
class(ctrl.use[[1]])
## [1] "character"
# There should be length(features.use)*ctrl genes (i.e. 20*100):
length(ctrl.use[[1]])


# Remove any repeated gene names - even though we set replace=FALSE when we sampled genes from the same expression bin, there may be more than two genes in our input gene list that fall in the same expression bin, so we can end up sampling the same gene more than once.
ctrl.use <- lapply(X = ctrl.use, FUN = unique)


## Get control gene scores

# Create an empty matrix with dimensions;
# number of rows equal to the number of gene sets (just one here)
# number of columns equal to number of cells in input Seurat object
ctrl.scores <- matrix(data = numeric(length = 1L),
                      nrow = length(x = ctrl.use),
                      ncol = ncol(x = object))

# Loop through each provided gene set and add to the empty matrix the mean expression of the control genes in each cell
for (i in 1:length(ctrl.use)) {
  # Get control gene names as a vector  
  features.use <- ctrl.use[[i]]
  # For each cell, calculate the mean expression of *all* of the control genes 
  ctrl.scores[i, ] <- Matrix::colMeans(x = assay.data[features.use,])
}


## Get scores for input gene sets

# Similar to the above, create an empty matrix
features.scores <- matrix(data = numeric(length = 1L),
                          nrow = cluster.length,
                          ncol = ncol(x = object))

# Loop through input gene sets and calculate the mean expression of these genes for each cell
for (i in 1:cluster.length) {
  features.use <- features[[i]]
  data.use <- assay.data[features.use, , drop = FALSE]
  features.scores[i, ] <- Matrix::colMeans(x = data.use)
}