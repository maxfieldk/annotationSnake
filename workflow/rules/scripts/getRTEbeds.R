source("~/data/common/myDefaults.r")
library(GenomicRanges)
library(rtracklayer)

conf <- c(
    conf <- configr::read.config(file = "conf/config.yaml")
)

tryCatch(
    {
        inputs <- snakemake@input
        outputs <- snakemake@output
        print("sourced snake variables")
    },
    error = function(e) {
        print("not sourced snake variables")
        assign("inputs", list(
            "r_annotation_fragmentsjoined" = "annotations/repeatmasker.gtf.rformatted.fragmentsjoined.csv",
            "r_repeatmasker_annotation" = "annotations/repeatmasker_annotation.csv"
        ), env = globalenv())
        assign("outputs", list(
            "outfile" = "annotations/rte_beds/outfile.txt"
        ), env = globalenv())
    }
)

outputdir <- dirname(outputs$outfile)
dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)

r_annotation_fragmentsjoined <- read_csv(inputs$r_annotation_fragmentsjoined)
r_repeatmasker_annotation <- read_csv(inputs$r_repeatmasker_annotation)
ann <- r_annotation_fragmentsjoined %>%
    left_join(r_repeatmasker_annotation)

### ONTOLOGY DEFINITION
{
    annot_colnames <- colnames(r_repeatmasker_annotation)
    annot_colnames_good <- annot_colnames[!(annot_colnames %in% c("gene_id", "family"))]
    ontologies <- annot_colnames_good[str_detect(annot_colnames_good, "family")]
    small_ontologies <- ontologies[grepl("subfamily", ontologies)]

    big_ontologies <- ontologies[!grepl("subfamily", ontologies)]
    big_ontology_groups <- c()
    for (ontology in big_ontologies) {
        big_ontology_groups <- c(big_ontology_groups, ann %>%
            pull(!!sym(ontology)) %>%
            unique())
    }
    big_ontology_groups <- big_ontology_groups %>% unique()

    modifiers <- annot_colnames_good[!str_detect(annot_colnames_good, "family")]
    region_modifiers <- modifiers[str_detect(modifiers, "_loc$")]
    element_req_modifiers <- modifiers[str_detect(modifiers, "_req$")]
}

groups_that_have_been_run <- c()
groups_not_to_run <- c()
for (ontology in ontologies) {
    ontology_groups <- ann %>%
        pull(!!sym(ontology)) %>%
        unique()
    ontology_groups <- ontology_groups[ontology_groups != "Other"]
    for (group in ontology_groups) {
        if (!(group %in% groups_that_have_been_run | group %in% groups_not_to_run | group %in% big_ontology_groups)) {
            groups_that_have_been_run <- c(groups_that_have_been_run, group)
            groupframe <- ann %>% filter(!!sym(ontology) == group)
            eligible_modifiers <- c()
            for (modifier in modifiers) {
                values_present <- ann %>%
                    filter(!!sym(ontology) == group) %>%
                    pull(!!sym(modifier)) %>%
                    unique()
                if ((length(values_present) > 1) | !("Other" %in% values_present)) {
                    eligible_modifiers <- c(eligible_modifiers, modifier)
                }
            }
            for (modifier in eligible_modifiers) {
                values_present <- ann %>%
                    filter(!!sym(ontology) == group) %>%
                    pull(!!sym(modifier)) %>%
                    unique()
                if ((length(values_present) > 1) | !("Other" %in% values_present)) {
                    eligible_modifiers <- c(eligible_modifiers, modifier)
                }
                eligible_filter_modifiers <- c(eligible_modifiers[grepl("_req$", eligible_modifiers)], "ALL")
                eligible_facet_modifiers <- c(eligible_modifiers[grepl("_loc$", eligible_modifiers)], "ALL")
                eligible_modifier_combinations <- expand.grid(filter_var = eligible_filter_modifiers, facet_var = eligible_facet_modifiers, stringsAsFactors = FALSE)
            }
            for (i in seq(1, length(rownames(eligible_modifier_combinations)))) {
                filter_var <- eligible_modifier_combinations[i, ]$filter_var
                facet_var <- eligible_modifier_combinations[i, ]$facet_var
                groupframe <- ann %>% filter(!!sym(ontology) == group)
                if (filter_var != "ALL") {
                    groupframe <- groupframe %>% filter(str_detect(!!sym(filter_var), ">|Intact"))
                }
                if (facet_var != "ALL") {
                    facet_values <- groupframe %>%
                        pull(!!sym(facet_var)) %>%
                        unique()
                    for (facet_value in facet_values) {
                        ranges <- groupframe %>%
                            filter(!!sym(facet_var) == facet_value) %>%
                            dplyr::select(seqnames, start, end, gene_id, pctdiv, strand) %>%
                            dplyr::rename(name = gene_id, score = pctdiv) %>%
                            GRanges()
                        export(ranges, paste0(outputdir, "/", group, "_", filter_var, "_", facet_value, ".bed"))
                    }
                } else {
                    ranges <- groupframe %>%
                        dplyr::select(seqnames, start, end, gene_id, pctdiv, strand) %>%
                        dplyr::rename(name = gene_id, score = pctdiv) %>%
                        GRanges()
                    export(ranges, paste0(outputdir, "/", group, "_", filter_var, "_ALL.bed"))
                }
            }
        }
    }
}

x <- data.frame()
write.table(x, file = outputs$outfile, col.names = FALSE)
