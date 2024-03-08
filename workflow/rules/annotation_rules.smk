def inputToMergeGTFs(wildcards):
      checkpoint_output = checkpoints.split_genome_fasta.get(**wildcards).output[0]
      return expand("repeatmasker/ref/{chromosome}/ref.part_{chromosome}.fa.gtf",
                  chromosome=glob_wildcards(os.path.join(checkpoint_output, "ref.part_{chromosome}.fa")).chromosome)

checkpoint split_genome_fasta:
    input:
        ref = config["ref"],
        ref_index = config["ref_index"]

    output:
        directory("ref/split")
    conda:
        "evo"
    shell:
        """
cp {input.ref} ref.fa
cp {input.ref_index} ref.fa.fai
mkdir -p $(dirname {output})
mkdir -p ref/split
seqkit split --by-id ref.fa -O ref/split/
        """

rule repeatmasker:
    input:
        chr_fasta = "ref/split/ref.part_{chromosome}.fa"
    output:
        rmout = "repeatmasker/ref/{chromosome}/ref.part_{chromosome}.fa.out"
    resources:
        cpus_per_task =20,
        runtime = 1200,
        mem_mb = 50000
    conda:
        "omics"
    shell:
        """
mkdir -p $(dirname {output})
/oscar/data/jsedivy/mkelsey/tools/RepeatMasker/RepeatMasker -species human -pa {resources.cpus_per_task} -gff {input.chr_fasta} -dir repeatmasker/ref/{wildcards.chromosome}
        """

rule getGtfs:
    input:
        rmout = "repeatmasker/ref/{chromosome}/ref.part_{chromosome}.fa.out"
    output:
        "repeatmasker/ref/{chromosome}/ref.part_{chromosome}.fa.gtf"
    conda:
        "evo"
    shell:
        """
workflow/scripts/outToGtf.sh {input.rmout} {output}
        """

rule mergeGtfsandCleanupRM:
    input:
        inputToMergeGTFs
    output:
        gtf = "repeatmasker/repeatmasker.gtf"
    conda:
        "evo"
    shell:
        """
cat {input} > {output.gtf}
sort -k1,1V -k4,4n -k5,5n {output.gtf} > tmp.gtf
ls -d RM_* | xargs rm -r
mv tmp.gtf {output.gtf}
        """


rule process_gtf:
    input:
        gtf = "repeatmasker/repeatmasker.gtf"
    output:
        repmask_gff2 = "annotations/repeatmasker.gff2",
        repmask_gff3 = "annotations/repeatmasker.gff3",
        r_annotation = "annotations/repeatmasker.gtf.rformatted.csv",
        r_annotation_fragmentsjoined = "annotations/repeatmasker.gtf.rformatted.fragmentsjoined.csv"
    conda:
        "evo2"
    script:
        "scripts/process_gtf.R"

rule move_refseq:
    input:
        ref_refseq_gff3 = config["ref_refseq_gff3"]
    output:
        refseq_gff3 = "annotations/refseq.gff3"
    shell:
        """
cp {input.ref_refseq_gff3} {output.refseq_gff3}
        """


rule complete_gff3:
    input:
        gff3 = "annotations/{annotation}.gff3",
    output:
        complete_gff3 = "annotations/{annotation}.complete.gff3",
    resources:
        mem_mb = 128000,
        runtime = 300
    conda:
        "omics"
    shell:
        """
agat_convert_sp_gxf2gxf.pl --gff {input.gff3} -o {output.complete_gff3}.temp
awk '!/#/ {{print}}' {output.complete_gff3}.temp | awk '{{FS="\t";OFS="\t"}} $4 < 900000000000 {{print}}' | sort -k1,1V -k4,4n -k5,5n > {output.complete_gff3}
rm  {output.complete_gff3}.temp
        """


rule gff_to_gtf:
    input:
        complete_gff3 = "annotations/{annotation}.complete.gff3",
    output:
        gtf = "annotations/{annotation}.complete.gtf",
    conda:
        "omics"
    resources:
        mem_mb = 128000,
        runtime = 300
    shell:
        """
agat_convert_sp_gff2gtf.pl --gff {input.complete_gff3} -o {output.gtf}.temp
awk '!/#/ {{print}}' {output.gtf}.temp | awk '{{FS="\t";OFS="\t"}} $4 < 900000000000 {{print}}' | sort -k1,1V -k4,4n -k5,5n > {output.gtf}
rm {output.gtf}.temp
        """

rule gff_to_bed12:
    input:
        complete_gff3 = "annotations/{annotation}.complete.gff3",
    output:
        bed12 = "annotations/{annotation}.complete.bed"
    resources:
        mem_mb = 128000,
        runtime = 300
    conda:
        "omics"
    shell:
        """
agat_convert_sp_gff2bed.pl --gff {input.complete_gff3} -o {output.bed12}.temp
awk '!/#/ {{print}}' {output.bed12}.temp | awk '{{FS="\t";OFS="\t"}} $4 < 900000000000 {{print}}' | sort -k1,1V -k2,2n -k3,3n > {output.bed12}
rm {output.bed12}.temp
        """

rule merge_genes_and_repeats_gff:
    input:
        complete_repeat_gff3 = "annotations/repeatmasker.complete.gff3",
        complete_refseq_gff3 = "annotations/refseq.complete.gff3"
    output:
        merged_gff3 = "annotations/repeatmasker_refseq.complete.gff3"
    resources:
        mem_mb = 128000,
        runtime = 300
    conda:
        "omics"
    shell:
        """
agat_sp_merge_annotations.pl -f {input.complete_refseq_gff3} -f {input.complete_repeat_gff3} -o {output.merged_gff3}.temp
awk '!/#/ {{print}}' {output.merged_gff3}.temp | awk '{{FS="\t";OFS="\t"}} $4 < 900000000000 {{print}}' | sort -k1,1V -k4,4n -k5,5n >  {output.merged_gff3}
rm {output.merged_gff3}.temp
        """

rule merge_OG_genes_and_repeats_gff3:
    input:
        repeatmasker = "annotations/repeatmasker.gff3",
        refseq = "annotations/refseq.gff3"
    output:
        merged_gff3 = "annotations/repeatmasker_refseq.gff3"
    resources:
        mem_mb = 128000,
        runtime = 300
    conda:
        "omics"
    shell:
        """
awk '!/#/ {{print}}' {input.repeatmasker} | awk '{{FS="\t";OFS="\t"}} $4 < 900000000000 {{print}}' > {output.merged_gff3}.temp
awk '!/#/ {{print}}' {input.refseq} | awk '{{FS="\t";OFS="\t"}} $4 < 900000000000 {{print}}' >> {output.merged_gff3}.temp
cat {output.merged_gff3}.temp | sort -k1,1V -k4,4n -k5,5n > {output.merged_gff3}
rm {output.merged_gff3}.temp
        """


rule merge_genes_and_repeats_gtf:
    input:
        complete_repeat_gtf = "annotations/repeatmasker.complete.gtf",
        complete_refseq_gtf = "annotations/refseq.complete.gtf"
    output:
        merged_gtf = "annotations/repeatmasker_refseq.complete.gtf"
    resources:
        mem_mb = 128000,
        runtime = 300
    conda:
        "omics"
    shell:
        """
agat_sp_merge_annotations.pl -f {input.complete_refseq_gtf} -f {input.complete_repeat_gtf} -o {output.merged_gtf}.temp
awk '!/#/ {{print}}' {output.merged_gtf}.temp | awk '{{FS="\t";OFS="\t"}} $4 < 900000000000 {{print}}' | sort -k1,1V -k4,4n -k5,5n > {output.merged_gtf}
rm {output.merged_gtf}.temp
        """

rule tabixindex:
    input:
        annot = "annotations/{annot}"
    output:
        gz = "annotations/{annot}.gz",
        index = "annotations/{annot}.gz.tbi"
    resources:
        mem_mb = 60000,
        runtime = 300
    conda:
        "omics"
    shell:
        """

awk '!/#/ {{print}}' {input.annot} | sort -k1,1V -k4,4n -k5,5n -t '\t'| bgzip > {input.annot}.gz
tabix -p gff {input.annot}.gz
        """

rule makeTxDB:
    input:
        refseq = "annotations/refseq.gff3",
        repeatmasker = "annotations/repeatmasker.complete.gff3"
    output:
        txdb = "annotations/repeatmasker_refseq.complete.sqlite",
        txdbrefseq = "annotations/refseq.sqlite",
        txdbrepeatmasker = "annotations/repeatmasker.complete.sqlite"
    resources:
        mem_mb = 40000
    conda:
        "repeatanalysis"
    script:
        "scripts/txdb.R"

rule get_transcriptome:
    input:
        gtf = "annotations/{annot}.gtf",
    output:
        fa = "annotations/{annot}.fa"
    params:
        genome = "ref.fa"
    resources:
        mem_mb = 128000,
        runtime = 300
    conda:
        "omics"
    shell:
        """
agat_sp_extract_sequences.pl -g {input.gtf} -f {params.genome} -t exon --merge -o {output.fa}
        """

rule annotate_rtes:
    input:
        r_annotation_fragmentsjoined = "annotations/repeatmasker.gtf.rformatted.fragmentsjoined.csv",
    output:
        r_repeatmasker_annotation = "annotations/repeatmasker_annotation.csv",
    conda:
        "evo2"
    script:
        "scripts/annotate_rtes.R"

rule getRTEbeds:
    input:
        r_annotation_fragmentsjoined = "annotations/repeatmasker.gtf.rformatted.fragmentsjoined.csv",
        r_repeatmasker_annotation = "annotations/repeatmasker_annotation.csv",
    output:
        outfile = "annotations/rte_beds/outfile.txt"
    conda:
        "evo2"
    script:
        "scripts/getRTEbeds.R"

rule element_analysis:
    input:
        r_annotation_fragmentsjoined = "annotations/repeatmasker.gtf.rformatted.fragmentsjoined.csv"
    params:
        l13 = config["l13fasta"]
    output:
        outfile = "RefAnalysis/l1element_evo.outfile"
    conda:
        "evo2"
    script:
        "scripts/element_analysis.R"


rule make_star_index:
    input:
        reference = "ref.fa"
    output:
        outfile = "ref_indeces/make_star_index.out"
    resources:
        cpus_per_task = 32,
        mem_mb = 128000,
        runtime = 300
    conda:
        "star"
    shell:
        """
mkdir -p $(dirname {output.outfile})
STAR --runThreadN  30 --runMode genomeGenerate --genomeDir $(dirname {output.outfile})/star_index --genomeFastaFiles {input.reference}
touch {output.outfile}
        """

rule cpgIslandFun:
    params:
        ref_cpgislands = config["ref_cpgislands"]
    output:
        cpg_islands_fullinfo = "annotations/cpg_islands.tsv",
        cpg_islands = "annotations/cpg_islands.bed",
        cpgi_shores = "annotations/cpgi_shores.bed",
        cpgi_shelves = "annotations/cpgi_shelves.bed"
    resources:
        cpus_per_task = 2,
        mem_mb = 20000,
        runtime = 60
    conda:
        "ds"
    script:
        "scripts/cpgIslandFun.R"

rule copySelectAnnotations:
    params:
        ref_cytobands = config["ref_cytobands"],
        ref_telomere = config["ref_telomere"]
    output:
        ref_cytobands = "annotations/cytobands.bed",
        ref_telomere = "annotations/telomeres.bed"
    shell:
        """
cp {params.ref_cytobands} {output.ref_cytobands}
cp {params.ref_telomere} {output.ref_telomere}
        """
