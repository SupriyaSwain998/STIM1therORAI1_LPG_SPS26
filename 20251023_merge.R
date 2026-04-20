getwd ()
setwd ("/shared/mendel/teams/laporte/supriya/Trancriptomics/01_Analysis/")
ext <- "-extractedReadCounts.tab"
files <- list.files(path = "04_Merged/", pattern = ext)

# Read in all files into a list of data frames
listOfFiles <- lapply(files, function(x) {
  read.table(paste0("04_Merged/", x), header = TRUE)
})

# Reduce the list of files by merging them by 'gene_ID' column
tbl <- Reduce(function(x, y) merge(x, y, all = TRUE, by = "gene_ID"), listOfFiles)

# Replace any NA values with 0
tbl[is.na(tbl)] <- 0

# Write the merged table to a CSV file
write.csv(tbl, file = "04_Merged/mergedReadCounts.csv", quote = FALSE, row.names = FALSE)
