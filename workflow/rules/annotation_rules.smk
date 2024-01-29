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


rule complete_gff3:
    input:
        gff3 = "annotations/repeatmasker.gff3",
    output:
        complete_gff3 = "annotations/repeatmasker.complete.gff3",
    resources:
        mem_mb = 128000,
        runtime = 300
    conda:
        "omics"
    shell:
        """
agat_convert_sp_gxf2gxf.pl --gff {input.gff3} --out {output.complete_gff3}
        """


rule gff_to_gtf:
    input:
        complete_gff3 = "annotations/repeatmasker.complete.gff3",
    output:
        gtf = "annotations/repeatmasker.complete.gtf",
    conda:
        "omics"
    resources:
        mem_mb = 128000,
        runtime = 300
    shell:
        """
agat_convert_sp_gff2gtf.pl --gff {input.complete_gff3} -o {output.gtf}
        """

rule gff_to_bed12:
    input:
        complete_gff3 = "annotations/repeatmasker.complete.gff3",
    output:
        bed12 = "annotations/repeatmasker.complete.bed"
    resources:
        mem_mb = 128000,
        runtime = 300
    conda:
        "omics"
    shell:
        """
agat_convert_sp_gff2bed.pl --gff {input.complete_gff3} -o {output.bed12}
        """


rule annotate_rtes:
    input:
        r_annotation_fragmentsjoined = "annotations/repeatmasker.gtf.rformatted.fragmentsjoined.csv",
    output:
        repeatmasker_annotation = "annotations/repeatmasker_annotation.csv",
    conda:
        "evo2"
    script:
        "scripts/annotate_rtes.R"

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