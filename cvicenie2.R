# Installation of packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

# Installation of biostrings
BiocManager::install("Biostrings")

# Load library
library('Biostrings')

###############################################
# Loading DNA file
# Seqinr method - result is a list
install.packages('seqinr')
library('seqinr')
seq_s <- read.fasta('fishes.fna.gz')

# Biostrings method - result is DNA stringset (Different Data structure)
seq_b <- readDNAStringSet('fishes.fna.gz')
  
# Explore fasta files
length(seq_b)
width(seq_b[1])  
names(seq_b)  
seq1 <- seq_b[1]  
seq1_sequence <- seq_b[[1]]  
seq1_string <- toString(seq_b[1])  
help(XStringSet)

# Global alignment of 2 sequences
seq6 <- seq_b[6]  
seq6_sequence <- seq_b[[6]] 

seq24 <- seq_b[24]  
seq24_sequence <- seq_b[[24]] 

align_seqs <- pairwiseAlignment(seq6_sequence, seq24_sequence, 
                                substitutionMatrix= 'BLOSUM62', 
                                gapOpening = -1, 
                                gapExtension = 1, 
                                scoreOnly = FALSE)

################################################################################
# Regular expressions are case sensitive! 
names_list <-  c("anna", "jana", "kamil", "norbert", "pavel", "petr", "stanislav", "zuzana")

# GREP returns position from the list satysfying the regex
grep("jana", names_list, perl = TRUE)

# Search for all names containing letter n at least once
grep("n+", names_list, perl = TRUE)

# Search for all names containing letters nn
grep("n{2}", names_list, perl = TRUE)

# Search for all names starting with n
grep("^n", names_list, perl = TRUE)

# Search for names Anna or Jana#
grep("Anna|Jana", names_list, perl = TRUE)
grep("anna|jana", names_list, perl = TRUE)
grep("anna|Jana", names_list, perl = TRUE)

# Search for names starting with z and ending with a:
grep("^z.*a$", names_list, perl=TRUE)

################################################################################
# Demultiplexing of sequencing data

fMID <-           'ACGAGTGCGT' 
rev_compl_rMID <- 'ACGCACTCGT'


RMID <-           'ACGAGTGCGT' 
rev_compl_rMID <- 'TGCTCACGCA'

paste0('^', fMID, '.*', rev_compl_rMID, '$')
is_in <- grep("^ACGAGTGCGT.+ACGCACTCGT$|^ACGAGTGCGT.+TGCTCACGCA$", seq_b, perl=TRUE)
print(length(is_in))

################################################################################
# Function 
desc_data <- read.csv('fishes_MIDs.csv', sep=';')

path <- 'fishes.fna.gz'
fMids <-desc_data[2] 
rMids <- desc_data[3]
sample_labels <- desc_data[1] 

demultiplexer <- function(path, fMids, rMids, sample_labels) {
  seq <- readDNAStringSet(path)
  
  # Create dataframe for generating txt file 
  labels <- sample_labels[1]
  counts <- c()
  for (i in 1:18)
  {
    # Create regular expression
    revComplRMid <- reverseComplement(DNAString(rMids[i, 1]))
    revComplFMid <- reverseComplement(DNAString(fMids[i, 1]))
    regex <- paste0('^', fMids[i,1],'.*', revComplRMid, '$',
                    '|^', rMids[i, 1], '.*', revComplFMid, '$' )

    # Find positions
    positions <- grep(regex, seq, perl=TRUE)
    count <- length(positions)
    counts <- c(counts, count)
    
    
    found_positions <- seq[positions]
    fMids_len <- nchar(fMids[i, 1])
    rMids_len <- nchar(rMids[i, 1])
    trimmed_sequences <- c()
    
    for (x in 1:length(found_positions))
    { sequence <- found_positions[[x]]
      last_pos <- length(sequence) - rMids_len
      trimmed_sequence <- sequence[fMids_len:last_pos]
      trimmed_sequences <- c(trimmed_sequences, sequence)
    }
    
    filename <- paste0('Results/Organisms/', sample_labels[i,1], '.fasta')
    write.fasta(trimmed_sequences, names(found_positions), filename, open = "w", nbchar = 60, as.string = FALSE)
  }
  
  
  results_df <- data.frame(labels, counts)
  write.table(results_df, "Results/results.txt", sep=",", row.names=FALSE)
}


demultiplexer(path, fMids, rMids, sample_labels)


