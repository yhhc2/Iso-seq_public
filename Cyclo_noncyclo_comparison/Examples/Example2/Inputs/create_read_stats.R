# Read the matrix from a CSV file
counts_matrix <- read.csv("matrix.csv", row.names = 1, check.names = FALSE)

# Create the output file
output_file <- "dummy_read_stats.txt"

# Open the file for writing
file_conn <- file(output_file, "w")

# Write the header
writeLines("id pbid", file_conn)

# Write the read IDs based on the counts matrix
for (pbid in rownames(counts_matrix)) {
  for (ps in colnames(counts_matrix)) {
    count <- counts_matrix[pbid, ps]
    for (i in 1:count) {
      writeLines(sprintf("%s_readID\t%s", ps, pbid), file_conn)
    }
  }
}

# Close the file connection
close(file_conn)
