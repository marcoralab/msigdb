library(tidyverse)
library(dbplyr)
library(plyranges)
library(here)

dir.create("input", showWarnings = FALSE)
dir.create("output", showWarnings = FALSE)

species = "human"
gencode_release <- "43"
msigdb_release <- "2023.2"

species <- str_to_lower(species)

if (species == "human") {
  gencode_species <- species
  msigdb_species <- "Hs"
} else if (species == "mouse") {
  gencode_species <- species
  msigdb_species <- "Mm"
} else {
  stop(str_glue("Species {species} not supported"))
}

gencode_url_base <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode"
msigdb_url_base <- "https://data.broadinstitute.org/gsea-msigdb/msigdb"

gencode_gtf <- str_glue("gencode.v{gencode_release}.basic.annotation.gtf")
msigdb_db <- str_glue("msigdb_v{msigdb_release}.{msigdb_species}.db")
msigdb_log <- str_glue("msigdb_v{msigdb_release}.{msigdb_species}.gencode.v{gencode_release}.basic.annotation.log")

if (species == "human") {
  msigdb_chip <- str_glue("Human_Ensembl_Gene_ID_MSigDB.v{msigdb_release}.{msigdb_species}.chip")
} else if (species == "mouse") {
  msigdb_chip <- str_glue("Human_Ensembl_Gene_ID_Mouse_Orthologs_MSigDB.v{msigdb_release}.{msigdb_species}.chip")
}

gencode_gtf_gz <- str_glue("{gencode_gtf}.gz")
msigdb_db_zip <- str_glue("{msigdb_db}.zip")

gencode_gtf_gz_url <- str_glue("{gencode_url_base}_{gencode_species}/release_{gencode_release}/{gencode_gtf_gz}")
msigdb_db_zip_url <- str_glue("{msigdb_url_base}/release/{msigdb_release}.{msigdb_species}/{msigdb_db_zip}")
msigdb_chip_url <- str_glue("{msigdb_url_base}/annotations/{species}/{msigdb_chip}")

gencode_gtf_gz_destfile <- here("input", gencode_gtf_gz)
msigdb_db_zip_destfile <- here("input", msigdb_db_zip)
msigdb_db_destfile <- here("input", msigdb_db)
msigdb_chip_destfile <- here("input", msigdb_chip)

# download (g)zipped files
if(!file.exists(gencode_gtf_gz_destfile)) download.file(url = gencode_gtf_gz_url, destfile = gencode_gtf_gz_destfile)
if(!file.exists(msigdb_db_destfile) && !file.exists(msigdb_db_zip_destfile)) download.file(url = msigdb_db_zip_url, destfile = msigdb_db_zip_destfile)
if(!file.exists(msigdb_chip_destfile)) download.file(url = msigdb_chip_url, destfile = msigdb_chip_destfile)

if(file.exists(msigdb_db_zip_destfile)) {
  # unzip zipped file
  unzip(msigdb_db_zip_destfile, exdir = here("input"))

  # remove zip file
  file.remove(msigdb_db_zip_destfile)
}

# process gencode gtf
gtf <- rtracklayer::import(gzfile(gencode_gtf_gz_destfile))

gencode.genes <- gtf |> filter(type == "gene") |> as_tibble()

stopifnot(nrow(gencode.genes) == nrow(distinct(gencode.genes, gene_id)))

gencode.genes <- gencode.genes |> filter(str_ends(gene_id, "_PAR_Y", negate = TRUE)) # seqnames %in% paste0("chr", 1:22)

gencode.genes |>
  group_by(gene_type) |>
  summarize(number = n()) |>
  print(n = Inf)

ensgenes <- gencode.genes |>
  select(gene_id, gene_name) |>
  separate(gene_id, into = c("ensgene", NA), sep = "\\.", remove = FALSE)

ensgenes

stopifnot(nrow(ensgenes) == nrow(distinct(ensgenes, gene_id)))

# process msigdb chip 
chip <- read_tsv(msigdb_chip_destfile, progress = FALSE, show_col_types = FALSE) |>
  select(ensgene = `Probe Set ID`, symbol = `Gene Symbol`)

nrow(chip)

nrow(ensgenes)

chip |> filter(ensgene %in% ensgenes$ensgene) |> nrow()

ensgenes |> inner_join(chip, by = "ensgene") |> filter(gene_name != symbol) |> print(n = Inf)

# process msigdb db
db <- DBI::dbConnect(RSQLite::SQLite(), dbname = msigdb_db_destfile, flags = RSQLite::SQLITE_RO)

RSQLite::dbListTables(db)

gene_set <- tbl(db, "gene_set") |> as_tibble()
gene_set_gene_symbol <- tbl(db, "gene_set_gene_symbol") |> as_tibble()
gene_symbol <- tbl(db, "gene_symbol") |> as_tibble()
gene_set_details <- tbl(db, "gene_set_details") |> as_tibble()

DBI::dbDisconnect(db)

# size of genesets in input
gene_set_size_input <- gene_set |>
  inner_join(gene_set_gene_symbol, join_by(id == gene_set_id)) |> 
  count(standard_name)

# join tables
gene_set.genes <- gene_set |>
  inner_join(gene_set_gene_symbol, join_by(id == gene_set_id)) |> 
  inner_join(gene_symbol, join_by(gene_symbol_id == id)) |> 
  inner_join(ensgenes, join_by(symbol == gene_name), relationship = "many-to-many") |>
  inner_join(gene_set_details, join_by(id == gene_set_id)) |>
  select(id, standard_name, systematic_name, collection_name, description_brief, symbol, NCBI_id, gene_id, ensgene) |>
  separate(collection_name, into = c("cat", "subcat"), sep = ":", remove = FALSE, extra = "merge")

gene_set.genes

# size of genesets in output
gene_set_size_output <- gene_set.genes |>
  count(standard_name)

# size of genesets in input and output, check for large mismatch
gene_set_size <- gene_set_size_input |> 
  left_join(gene_set_size_output, by = "standard_name", suffix = c("_input", "_output")) |>
  filter(abs(n_input - n_output) >= (n_input * 0.1)) |>
  write_tsv(msigdb_log) |> 
  print(n = Inf)

# write GMT files to output dir
for (col in unique(gene_set.genes$collection_name)) {
  suffix <- str_replace_all(col, "[[:punct:]]", "\\.")
  gene_set.genes |> filter(collection_name == col) |>
    select(standard_name, description_brief, gene_id) |>
    arrange(standard_name, gene_id) |>
    nest(data = gene_id) |>
    mutate(gene_ids = map_chr(data, ~ str_flatten(unique(.x$gene_id), collapse = "\t"))) |> 
    select(standard_name, description_brief, gene_ids) |>
    write_tsv(str_glue("output/msigdb_v{msigdb_release}.{msigdb_species}.gencode.v{gencode_release}.basic.annotation.{suffix}.gmt"), col_names = FALSE)
}

for (col in unique(gene_set.genes$cat)) {
  suffix <- col
  gene_set.genes |> filter(cat == col) |>
    select(standard_name, description_brief, gene_id) |>
    arrange(standard_name, gene_id) |>
    nest(data = gene_id) |>
    mutate(gene_ids = map_chr(data, ~ str_flatten(unique(.x$gene_id), collapse = "\t"))) |> 
    select(standard_name, description_brief, gene_ids) |>
    write_tsv(str_glue("output/msigdb_v{msigdb_release}.{msigdb_species}.gencode.v{gencode_release}.basic.annotation.{suffix}.gmt"), col_names = FALSE)
}

################################################################################

read_gmt <- function (file, start = 1, end = -1) 
{
  if (!grepl("\\.gmt$", file)[1]) {
    stop("Pathway information must be a .gmt file")
  }
  geneSetDB = scan(file, what = "character", n = end - start + 
                     1, skip = start - 1, sep = "\n", quiet = TRUE)
  geneSetDB = suppressWarnings(strsplit(geneSetDB, "\t"))
  names(geneSetDB) = sapply(geneSetDB, "[", 1)
  geneSetDB = lapply(geneSetDB, "[", -1:-2)
  geneSetDB = lapply(geneSetDB, function(x) {
    x[which(x != "")]
  })
  geneSetDB <- geneSetDB[sapply(geneSetDB, length) > 0 & !is.na(names(geneSetDB))]
  return(geneSetDB)
}