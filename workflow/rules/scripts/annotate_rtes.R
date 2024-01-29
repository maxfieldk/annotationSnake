source("~/data/common/myDefaults.r")
library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(ORFik)

conf <- c(
    conf <- configr::read.config(file = "conf/config.yaml")
)
tryCatch(
    {
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("inputs", list(
            r_annotation_fragmentsjoined = "annotations/repeatmasker.gtf.rformatted.fragmentsjoined.csv"
        ), env = globalenv())
        assign("outputs", list(
            repeatmasker_annotation = "annotations/repeatmasker_annotation.csv"
        ), env = globalenv())
    }
)

rmfragments <- read_csv(inputs$r_annotation_fragmentsjoined, col_names = TRUE)

rmfamilies <- rmfragments %>%
    dplyr::select(gene_id, family) %>%
    mutate(repeat_superfamily = case_when(
        str_detect(family, "^LINE") ~ "LINE",
        str_detect(family, "^SINE") ~ "SINE",
        str_detect(family, "^LTR") ~ "LTR",
        str_detect(family, "^DNA") ~ "DNA",
        str_detect(family, "^Satellite") ~ "SAT",
    )) %>%
    mutate(repeat_superfamily = replace_na(repeat_superfamily, "Other")) %>%
    mutate(rte_superfamily = case_when(
        str_detect(family, "^LINE") ~ "LINE",
        str_detect(family, "^SINE") ~ "SINE",
        str_detect(family, "^LTR") ~ "LTR",
    )) %>%
    mutate(rte_superfamily = replace_na(rte_superfamily, "Other")) %>%
    mutate(rte_family = case_when(
        str_detect(family, "^LINE/L1") ~ "L1",
        str_detect(family, "^SINE/Alu") ~ "Alu",
        str_detect(family, "^LTR/ERV") ~ "ERV",
    )) %>%
    mutate(rte_family = replace_na(rte_family, "Other")) %>%
    mutate(rte_subfamily = case_when(
        grepl("LTR(.)*LTR/ERVK$", family, perl = TRUE) ~ "HERVK",
        grepl("HERVK(.)*int$", family, perl = TRUE) ~ "HERVK",
        grepl("LTR(.)*LTR/ERVL$", family, perl = TRUE) ~ "HERVL",
        grepl("HERVL(.)*int$", family, perl = TRUE) ~ "HERVL",
        grepl("LTR5#LTR/ERVK$", family, perl = TRUE) ~ "HERVK",
        grepl("LTR5A.*ERVK$", family, perl = TRUE) ~ "HERVK",
        grepl("LTR5B.*ERVK$", family, perl = TRUE) ~ "HERVK",
        grepl("LTR5_Hs.*ERVK$", family, perl = TRUE) ~ "HERVK",
        grepl("^SINE/Alu/AluY$", family, perl = TRUE) ~ "AluY",
        grepl("^LINE/L1/L1HS$", family, perl = TRUE) ~ "L1HS",
        grepl("^LINE/L1/L1PA2$", family, perl = TRUE) ~ "L1PA2",
        grepl("^LINE/L1/L1PA3$", family, perl = TRUE) ~ "L1PA3",
        grepl("^LINE/L1/L1PA4$", family, perl = TRUE) ~ "L1PA4",
        grepl("^LINE/L1/L1PA5$", family, perl = TRUE) ~ "L1PA5",
        grepl("^LINE/L1/L1PA6$", family, perl = TRUE) ~ "L1PA6"
    )) %>%
    mutate(rte_subfamily = replace_na(rte_subfamily, "Other")) %>%
    mutate(rte_subfamily_limited = case_when(
        grepl("LTR(.)*LTR/ERVK$", family, perl = TRUE) ~ "HERVK",
        grepl("HERVK(.)*int$", family, perl = TRUE) ~ "HERVK",
        grepl("LTR5#LTR/ERVK$", family, perl = TRUE) ~ "HERVK",
        grepl("LTR5A.*ERVK$", family, perl = TRUE) ~ "HERVK",
        grepl("LTR5B.*ERVK$", family, perl = TRUE) ~ "HERVK",
        grepl("LTR5_Hs.*ERVK$", family, perl = TRUE) ~ "HERVK",
        grepl("^SINE/Alu/AluY$", family, perl = TRUE) ~ "AluY",
        grepl("^LINE/L1/L1HS$", family, perl = TRUE) ~ "L1HS"
    )) %>%
    mutate(rte_subfamily_limited = replace_na(rte_subfamily_limited, "Other")) %>%
    mutate(l1_subfamily_limited = case_when(
        grepl("^LINE/L1/L1HS$", family, perl = TRUE) ~ "L1HS",
        grepl("^LINE/L1/L1PA2$", family, perl = TRUE) ~ "L1PA2",
        grepl("^LINE/L1/L1PA3$", family, perl = TRUE) ~ "L1PA3",
        grepl("^LINE/L1/L1PA4$", family, perl = TRUE) ~ "L1PA4",
        grepl("^LINE/L1/L1PA5$", family, perl = TRUE) ~ "L1PA5",
        grepl("^LINE/L1/L1PA6$", family, perl = TRUE) ~ "L1PA6"
    )) %>%
    mutate(l1_subfamily_limited = replace_na(l1_subfamily_limited, "Other"))


rmlengthreq <- rmfragments %>%
    full_join(rmfamilies) %>%
    dplyr::select(gene_id, rte_subfamily, length) %>%
    mutate(rte_length = case_when(
        str_detect(rte_subfamily, "L1HS") & length > 5999 ~ "L1HS >6kb",
        str_detect(rte_subfamily, "L1HS") & length < 5999 ~ "L1HS <6kb",
        str_detect(rte_subfamily, "L1PA2") & length > 5999 ~ "L1PA2 >6kb",
        str_detect(rte_subfamily, "L1PA2") & length < 5999 ~ "L1PA2 <6kb",
        str_detect(rte_subfamily, "L1PA3") & length > 5999 ~ "L1PA3 >6kb",
        str_detect(rte_subfamily, "L1PA3") & length < 5999 ~ "L1PA3 <6kb",
        str_detect(rte_subfamily, "L1PA4") & length > 5999 ~ "L1PA4 >6kb",
        str_detect(rte_subfamily, "L1PA4") & length < 5999 ~ "L1PA4 <6kb",
        str_detect(rte_subfamily, "L1PA5") & length > 5999 ~ "L1PA5 >6kb",
        str_detect(rte_subfamily, "L1PA5") & length < 5999 ~ "L1PA5 <6kb",
        str_detect(rte_subfamily, "L1PA6") & length > 5999 ~ "L1PA6 >6kb",
        str_detect(rte_subfamily, "L1PA6") & length < 5999 ~ "L1PA6 <6kb",
        str_detect(rte_subfamily, "HERVL") & length > 6999 ~ "HERVL >7kb",
        str_detect(rte_subfamily, "HERVL") & length < 6999 ~ "HERVL <7kb",
        str_detect(rte_subfamily, "HERVK") & length > 6999 ~ "HERVK >7kb",
        str_detect(rte_subfamily, "HERVK") & length < 6999 ~ "HERVK <7kb",
        str_detect(rte_subfamily, "AluY") & length > 299 ~ "AluY >300bp",
        str_detect(rte_subfamily, "AluY") & length < 299 ~ "AluY <300bp",
        str_detect(rte_subfamily, "Other") ~ "Other",
    )) %>%
    dplyr::select(gene_id, rte_length)


# Annotate Intactness
fa <- Rsamtools::FaFile("ref.fa")

ranges <- GRanges(rmfragments)
l1hsranges <- ranges[grepl("L1HS", ranges$gene_id)]
l1pa2ranges <- ranges[grepl("L1PA2", ranges$gene_id)]
l1hspa2ranges <- c(l1hsranges, l1pa2ranges)
gene_ids <- l1hspa2ranges$gene_id
l1hspa2ss <- getSeq(fa, l1hspa2ranges)
mcols(l1hspa2ss) <- mcols(l1hspa2ranges)
names(l1hspa2ss) <- mcols(l1hspa2ss)$gene_id
flss <- l1hspa2ss[width(l1hspa2ss) > 5999]

# orf analysis
pos <- ORFik::findORFs(flss, startCodon = "ATG", longestORF = TRUE, minimumLength = 333)
names(pos) <- names(flss[as.integer(names(pos))])
# FILTERING ON THE BASIS OF ORF LENGTHS
pass <- c()
for (element in names(pos)) {
    if (3828 %in% width(pos[[element]]) & 1017 %in% width(pos[[element]])) {
        pass <- c(pass, element)
    }
}
l1hspass <- pass[pass %>% str_detect("L1HS")]
l1pa2pass <- pass[pass %>% str_detect("L1PA2")]

intactness <- rmfragments %>%
    dplyr::select(gene_id) %>%
    mutate(l1_intactness = case_when(
        gene_id %in% l1hspass ~ "L1HS ORFs Intact",
        gene_id %in% l1pa2pass ~ "L1PA2 ORFs Intact",
        TRUE ~ "Other"
    ))



# rmfamilies <- read_csv(outputs$r_annotation_families, col_names = TRUE)

# for annotation purposes, I will have to have the location of nonreference inserts be their insertion site
rmfragments_ref <- rmfragments %>% filter(!str_detect(seqnames, "nonref"))
rmfragments_nonref <- rmfragments %>% filter(str_detect(seqnames, "nonref"))
seqnamesNonref <- rmfragments_nonref %$% seqnames
insert_seqnames <- c()
insert_start <- c()
insert_end <- c()
for (seqnames in seqnamesNonref) {
    split <- str_split(seqnames, "_") %>% unlist()
    splitlen <- length(split)
    insert_seqnames <- c(insert_seqnames, split[splitlen - 2])
    insert_start <- c(insert_start, split[splitlen - 1])
    insert_end <- c(insert_end, split[splitlen])
}
rmfragments_nonref$insert_seqnames <- insert_seqnames
rmfragments_nonref$insert_start <- insert_start
rmfragments_nonref$insert_end <- insert_end

rmfragments_nonrefgr <- GRanges(rmfragments_nonref %>% dplyr::select(-seqnames, -start, -end) %>% dplyr::relocate(gene_id, insert_seqnames, source, type, insert_start, insert_end, strand) %>% dplyr::rename(seqnames = insert_seqnames, start = insert_start, end = insert_end))
rmfragments_refgr <- GRanges(rmfragments_ref)
rmfragmentsgr_properinsertloc <- c(rmfragments_refgr, rmfragments_nonrefgr)


# genes, centromere and telomere
refseq <- import(conf$ref_refseq)
refseqdf <- as.data.frame(refseq) %>% tibble()
coding_transcripts <- refseq[(mcols(refseq)$type == "transcript" & grepl("^NM", mcols(refseq)$transcript_id))]

cytobandsdf <- read_delim(conf$ref_cytobands, col_names = FALSE, delim = "\t")
cytobands <- cytobandsdf %>%
    dplyr::rename(seqnames = X1, start = X2, end = X3) %>%
    GRanges()
centromere <- cytobandsdf %>%
    filter(X5 == "acen") %>%
    dplyr::rename(seqnames = X1, start = X2, end = X3) %>%
    GRanges()
telomeredf <- read_delim(conf$ref_telomere, col_names = FALSE, delim = "\t")
telomere <- telomeredf %>%
    dplyr::rename(seqnames = X1, start = X2, end = X3) %>%
    GRanges()

getannotation <- function(to_be_annotated, regions_of_interest, property, name_in, name_out) {
    inregions <- to_be_annotated %>% subsetByOverlaps(regions_of_interest, invert = FALSE)
    inregions$prop <- name_in
    outregions <- to_be_annotated %>% subsetByOverlaps(regions_of_interest, invert = TRUE)
    outregions$prop <- name_out
    merged <- c(inregions, outregions)
    annot <- merged %>%
        as.data.frame() %>%
        tibble() %>%
        dplyr::select(gene_id, prop)
    annot[[property]] <- annot$prop
    annot <- annot %>% dplyr::select(-prop)
    return(annot)
}

genic_annot <- getannotation(rmfragmentsgr_properinsertloc, coding_transcripts, "genic", "Genic", "NonGenic")
cent_annot <- getannotation(rmfragmentsgr_properinsertloc, centromere, "centromeric", "Centromeric", "NonCentromeric")
telo_annot <- getannotation(rmfragmentsgr_properinsertloc, telomere, "telomeric", "Telomeric", "NonTelomeric")
genic_annot %$% genic %>% table()
cent_annot %$% centromeric %>% table()
telo_annot %$% telomeric %>% table()

region_annot <- full_join(genic_annot, cent_annot) %>% full_join(telo_annot)

annots %>% filter(rte_subfamily == "L1HS")

annots <- rmfamilies %>%
    full_join(rmlengthreq %>% rename_at(vars(-gene_id), ~ paste0(., "_req"))) %>%
    full_join(intactness %>% rename_at(vars(-gene_id), ~ paste0(., "_req"))) %>%
    full_join(region_annot %>% rename_at(vars(-gene_id), ~ paste0(., "_loc")))

write_csv(annots, outputs$repeatmasker_annotation)
