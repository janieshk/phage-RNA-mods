# In silico search for RluA-modified stem-loop substrates

# To run script: Rscript UURHA_search.R <fasta> [--outdir=DIR] [--outstem=NAME] [--nodraw]
## <fasta> --> fasta file containing sequence(s) of interest
## [--outdir=DIR] --> specify desired output directory
## [--outstem=NAME] --> file naming tag
## [--nodraw] --> specify "nodraw"  if do not want secondary structures of hits drawn and outputted as .svg files

## Note: Change below bin locations to match user's locations
## Note: Cannot have duplicate identical FASTA IDs


library (Biostrings)

RNALfold_bin <- "/media/data0/janiek/miniforge3/envs/rna_struct/bin/RNALfold"
RNAfold_bin <- "/media/data0/janiek/miniforge3/envs/rna_struct/bin/RNAfold"
RNAeval_bin <- "/media/data0/janiek/miniforge3/envs/rna_struct/bin/RNAeval"
RNAPLOT_BIN <- "/media/data0/janiek/miniforge3/envs/rna_struct/bin/RNAplot"

# Defining parameters
stem_min <- 3
stem_max <- 6
loop_min <-4
loop_max <- 16
sequence_motif <- "UU[AG][ACU]A"
RNALfold_L <- 50


# Defining restrictions on exact orientation of UURHA in the loop
# 'off = 0' means the first U is the 5' base of the closing pair.
ok_offset <- function(loop_len, off) {
    if (loop_len == 4) return(off == 0)
    if (loop_len == 5) return(off == 0)
    if (loop_len == 6) return(off == 1)
    if (loop_len == 7) return(off == 1)
    if (loop_len >= 8  && loop_len <= 10) return(off == 2)
    if (loop_len >= 11 && loop_len <= 14) return(off >= 2 && off <= 4)
    if (loop_len == 15) return(off >= 2 && off <= 5)
    if (loop_len == 16) return(off >= 2 && off <= 5)
    FALSE
}


# 1. Function to call RNALfold on one input sequence and give the dot-bracket notation
run_rnalfold <- function(seq_str, L = RNALfold_L) {
    out <- system2(RNALfold_bin, args = c(paste0("-L", L)), input = seq_str, stdout = TRUE)

    rx <- "^([.()]+)\\s+\\(\\s*(-?[0-9.]+)\\)\\s+([0-9]+)\\s*$" # regex to grab the dot-bracket notation line from output
    m  <- regmatches(out, regexec(rx, out)) # find where the ^regex matches
    keep <- lengths(m) == 4
    if (!any(keep)) return(data.frame()) # return empty df if no regex matches

    # make a df for every matched line that contains columns structure, energy, start
    do.call(rbind, lapply(m[keep], function(x) {
        data.frame(structure = x[[2]], energy = as.numeric(x[[3]]), start = as.integer(x[[4]]), stringsAsFactors = FALSE)
    }))
}



# 2. Function to check for the hairpin + RluA motif
check_hairpin <- function(db, motif_start, motif_end) {
    chars <- strsplit(db, "")[[1]] # break the dot-bracket string into vector of single characters to be able to grab individual positions!!
    n <- length(chars)


    # 2a. Make a "map" of the base-pairing:
    ## function that will describe which bases pair with which bases:
    partner <- rep(NA_integer_, n) # initially just sticking an NA into every position (will be overwritten by basepaired bases only)
    stack <- integer(0)
    # "walk" from left to right, and the position of every ( is saved in 'stack'
    for (i in seq_len(n)) {
        if (chars[i] == "(") {
            stack <- c(stack, i)
        # every occurrence of ) pops the most recent unmatched ( off 'stack' and their positions are saved
        } else if (chars[i] == ")") {
            j <- stack[length(stack)]; stack <- stack[-length(stack)]
            partner[i] <- j; partner[j] <- i
        }
    }

    # 2b. Scan for a clean hairpin (aka ^ paired i,j that encloses only unpaired bases aka a loop region)
    ## loop through and skip over anything that isn't an opening (
    ## for each (, get its partner j
    ## if it has no partner, just skip
    for (i in seq_len(n)) {
        if (chars[i] != "(") next
        j <- partner[i]
        if (is.na(j) || j <= i) next


    # this is the candidate loop (aka the positions between the basepair i and j)
    loop_id <- (i + 1):(j - 1)
    # check that loop length is right
    if (length(loop_id) < loop_min || length(loop_id) > loop_max) next 
    # check that all those loop bases are indeed unpaired (should all be ".")
    if (!all(chars[loop_id] == ".")) next 


    # 2c. Measure the length of contiguous stem and check if it meets requirements
    ## starting with innermost closing pair (i,j) walk through the stem with a and b, keeping length saved in 'stem_len'
    stem_len <- 0; a <- i; b <- j
    while (a >= 1 && b <= n && chars[a] == "(" && chars[b] == ")" && !is.na(partner[a]) && partner[a] == b) {  # !is.na(partner[a]) is checking if partner exists; partner[a] == b is checking if a pairs with b
        stem_len <- stem_len + 1; a <- a - 1; b <- b + 1
    }
    if (stem_len < stem_min || stem_len > stem_max) next


    # 2d. Now enforce the RluA sequence motif constraints:

    ## defining constraints:
    loop_len <- j - i - 1
    off <- motif_start - i
    if (off < 0) next
    if (!ok_offset(loop_len, off)) next
    if (motif_end >= j) next   # (the full UURHA shouldn't extend past the loop)
    loop_start <- i + 1

    # if all these criteria are met, the hairpin's info is returned. If no valid hairpin found, then returns NULL
    return(list(loop_start = i + 1, loop_end = j - 1, stem_len = stem_len, close_i = i, close_j = j))
    }
    NULL
}




# 3. Extract the free energy of just the stemloop substructure
eval_structure <- function(seq, db) {
    out <- system2(RNAeval_bin, input = paste(seq, db, sep = "\n"), stdout = TRUE)
    as.numeric(sub(".*\\(\\s*(-?[0-9.]+)\\)\\s*$", "\\1", out[length(out)]))
}

# get stem-only (stacking) free energy (loop penalty excluded)
stem_only_energy <- function(seq, db) {
    out <- system2(RNAeval_bin, args = "-v", input = paste(seq, db, sep = "\n"), stdout = TRUE)
    sl <- grep("Interior loop", out, value = TRUE)
    if (!length(sl)) return(NA_real_)
    sum(as.numeric(sub(".*:\\s*(-?[0-9]+)\\s*$", "\\1", sl))) / 100   # dekacal -> kcal
}



# 3. Need to test every single outputted RNALfold structure outputted by RNALfold for all of above defined criteria. 
# Should return any substructures within those RNALfold structures as hits if meet criteria

## make a hits list to fill (if RNALfold didn't find any hits, then just returns empty df)
evaluate_sequence <- function(seq_id, full_seq, substructs) {
    hits <- list()
    if (nrow(substructs) == 0) return(data.frame())

## grab values from every RNALfold structure outputted:
for (r in seq_len(nrow(substructs))) {
    db <- substructs$structure[r] #the structure's dot-bracket
    coord <- substructs$start[r] #where the structure starts
    sublen <- nchar(db)
    subseq <- substr(full_seq, coord, coord + sublen - 1) #the actual sequence corresponding to the structure

    ## find all desired motif (UURHA) and their start positions
    ## if gregexpr returns -1, means none were found, so skip
    motif_pos <- gregexpr(sequence_motif, subseq)[[1]]
    if (motif_pos[1] == -1) next

    ## for each motif occurrence at position m (m + 4 is its end since UURHA is 5 nt long), run the check_hairpin function from above
    for (m in motif_pos) {
        v <- check_hairpin(db, m, m + 4)
        if (!is.null(v)) {

            # hairpin fragment = stem + loop + stem
            hp_start <- max(1, v$close_i - v$stem_len + 1)
            hp_end <- min(sublen, v$close_j + v$stem_len - 1)
            hp_seq <- substr(subseq, hp_start, hp_end)
            hp_db <- substr(db, hp_start, hp_end)

            genome_motif_start <- coord + m - 1 # position-within-RNALfold-struture converted to genomic coordinate
            hits[[length(hits) + 1]] <- data.frame(
                seq_id = seq_id,
                genome_motif_start = coord + m - 1,
                genome_motif_end = coord + m + 3,
                motif = substr(subseq, m, m + 4),
                stem_len = v$stem_len,
                loop_len = v$close_j - v$close_i - 1,
                offset = m - v$close_i,   # 0 = Psi is 5' base of closing pair
                loop_seq = substr(subseq, v$loop_start, v$loop_end),
                structure_start = coord,
                structure_energy = substructs$energy[r],
                structure = db,
                hairpin_seq = hp_seq,
                hairpin_db = hp_db,
                stemloop_energy = eval_structure(hp_seq, hp_db),
                stem_stack_energy = stem_only_energy(hp_seq, hp_db),
                stringsAsFactors = FALSE
            )
        }
    }
}
if (length(hits) == 0) return(data.frame())
df <- do.call(rbind, hits)

# if the same motif shows up in multiple overlapping substructures, then keep the lowest-energy/most stable structure per motif site
df <- df[order(df$genome_motif_start, df$structure_energy), ]
    df[!duplicated(df[, c("seq_id", "genome_motif_start")]), ]
}





# 4. Now running all the above:

# Either Option 1: Uncomment the below, replace "MS2_DNA.fasta" with the desired fasta FILE PATH, then source() this script. (It must be a path; the script reads it itself.)

#infile <- "MS2_DNA.fasta"

# Or Option 2: Default as is, command-line arguments read in
args <- commandArgs(trailingOnly = TRUE)

usage <- function(msg = NULL) {
    if (!is.null(msg)) cat("Error:", msg, "\n\n")
    cat("Usage: Rscript UURHA_search.R <fasta> [--outdir=DIR] [--outstem=NAME] [--nodraw]\n")
    quit(status = 1)
}

# get the value out of a --flag=value argument (default if the flag is absent)
get_flag <- function(flag, default = NULL) {
    hit <- grep(paste0("^", flag, "="), args, value = TRUE)
    if (!length(hit)) return(default)
    sub(paste0("^", flag, "="), "", hit[length(hit)])
}

DRAW <- !any(args == "--nodraw")   # --nodraw works if written in any position
positional <- args[!grepl("^--", args)]

# fasta comes from argument 1, or from an `infile` set by hand above (option 1)
if (length(positional) >= 1) {
    infile <- positional[1]
} else if (!exists("infile")) {
    usage("no input fasta given")
}
if (!file.exists(infile)) usage(paste0("input fasta not found: ", infile))

outstem <- get_flag("--outstem", tools::file_path_sans_ext(basename(infile)))

# all results go into <outdir>/<date>__<outstem>__RluA_hits/
OUTDIR <- file.path(get_flag("--outdir", "."), paste0(Sys.Date(), "__", outstem, "__RluA_hits"))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(OUTDIR)) usage(paste0("could not create output directory: ", OUTDIR))
OUTDIR <- normalizePath(OUTDIR)
cat("Writing results to", OUTDIR, "\n")



# (Can take DNA or RNA input, uppercase or lowercase)
seqs <- readBStringSet(infile)
chars <- gsub("T", "U", toupper(as.character(seqs)))
bad <- grepl("[^ACGU]", chars)
if (any(bad)) {
    usage(paste0("non-ACGT/U character(s) in ", sum(bad), " sequence(s); first offender: ",
                 names(seqs)[which(bad)[1]]))
}
seqs_rna <- RNAStringSet(chars)
output <- file.path(OUTDIR, paste0("RNALfold_stemloop_hits_", outstem, "_", Sys.Date(), ".csv"))

all_hits <- list()
first_write <- TRUE

for (k in seq_along(seqs_rna)) {
  sname <- names(seqs_rna)[k]
  full <- as.character(seqs_rna[[k]])
  cat(sprintf("Folding %s (%d nt)...\n", sname, nchar(full)))

  substructs <- run_rnalfold(full, L = RNALfold_L)
  res <- evaluate_sequence(sname, full, substructs)

  if (nrow(res) > 0) {
    write.table(res, output, sep = ",", row.names = FALSE, col.names = first_write, append = !first_write)
    first_write <- FALSE
    all_hits[[length(all_hits) + 1]] <- res
    cat(sprintf("  -> %d hit(s)\n", nrow(res)))
  } else {
    cat("  -> no hits\n")
  }
}

final <- if (length(all_hits)) do.call(rbind, all_hits) else data.frame()
cat(sprintf("\nDONE. %d total hit(s) written to %s\n", nrow(final), output))
if (nrow(final)) print(final[, c("seq_id","genome_motif_start","motif","stem_len","loop_seq","structure_energy")])


if (nrow(final)) {
    print(final[, c("seq_id","genome_motif_start","motif","stem_len","loop_len", "offset","loop_seq","stemloop_energy","stem_stack_energy")])
}

file_basename <- outstem


write.table(final, file = paste0(OUTDIR, "/", Sys.Date(), "_", file_basename, ".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)



# print out hits per sequence:
all_ids <- names(seqs_rna)
per_cds <- table(factor(final$seq_id, levels = all_ids))

cat(sprintf("\n%d hits across %d of %d sequences (%.2f per sequence, %.3f per kb)\n",
            nrow(final), sum(per_cds > 0), length(all_ids),
            nrow(final) / length(all_ids),
            1000 * nrow(final) / sum(width(seqs))))

cat("\ndistribution -- how many sequences have N hits:\n")
print(table(hits = as.integer(per_cds)))

cds_tab <- data.frame(
    seq_id = all_ids,
    length = width(seqs),
    hits = as.integer(per_cds),
    per_kb = round(1000 * as.integer(per_cds) / width(seqs), 3),
    stringsAsFactors = FALSE
)
cds_tab <- cds_tab[order(-cds_tab$hits, -cds_tab$per_kb), ]

write.table(cds_tab, file.path(OUTDIR, paste0(Sys.Date(), "_", file_basename, "_hits_per_cds.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)

cat("\ntop 20:\n"); print(head(cds_tab, 20), row.names = FALSE)


# 5. Draw the secondary structure (draw w ViennaRNA / RNAplot)
## Draws unless the --nodraw flag is passed in

if (DRAW) {

    draw_hit <- function(final, seqs_rna, r, outdir = OUTDIR) {

        sid <- final$seq_id[r]
        db <- final$structure[r]
        sub_st <- final$structure_start[r]
        sublen <- nchar(db)
        full <- as.character(seqs_rna[[sid]])
        subseq <- substr(full, sub_st, sub_st + sublen - 1)

        # line1 = sequence, line2 = dot-bracket
        input_txt <- paste(subseq, db, sep = "\n")

        old <- setwd(outdir)
        on.exit(setwd(old), add = TRUE)

        system2(RNAPLOT_BIN, args = c("--output-format=svg"),input = input_txt, stdout = TRUE)

        # renaming:
        fname <- sprintf("%s_pos%d_%s.svg", sid, final$genome_motif_start[r], final$motif[r])
        if (file.exists("rna.svg")) {
            file.rename("rna.svg", fname)
            cat("Wrote", file.path(outdir, fname), "\n")
        } else {
            cat("WARNING: RNAplot produced no rna.svg for", sid,
                "pos", final$genome_motif_start[r], "\n")
        }
    }

    # draw all hits:
    for (r in seq_len(nrow(final))) draw_hit(final, seqs_rna, r)
}
