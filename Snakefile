REFERENCE = "config/reference.fasta"
GENOME_ANNOTATION = "config/genome_annotation.gff3"
PATHOGEN_JSON = "config/pathogen.json"
TAXON_ID = 11089
DROPPED_STRAINS = "config/dropped_strains.txt"
MEDIOCRE_STRAINS = "config/mediocre_strains.txt"

COLORS = ("config/colors.tsv",)
LAT_LONGS = ("config/lat_longs.tsv",)
CLADES = "config/clades_new.tsv"
AUSPICE_CONFIG = "config/auspice_config.json"

rule all:
    input:
        "auspice/yellow-fever.json",
        "dataset.zip",

rule fetch_ncbi_dataset_package:
    output:
        dataset_package="data/ncbi_dataset.zip",
    shell:
        """
        datasets download virus genome taxon {TAXON_ID} \
            --no-progressbar \
            --filename {output.dataset_package} \
        """


rule extract_ncbi_dataset_sequences:
    input:
        dataset_package=rules.fetch_ncbi_dataset_package.output.dataset_package,
    output:
        ncbi_dataset_sequences="data/sequences.fasta",
    shell:
        """
        unzip -jp {input.dataset_package} ncbi_dataset/data/genomic.fna > {output.ncbi_dataset_sequences}
        """


rule format_ncbi_dataset_report:
    input:
        dataset_package=rules.fetch_ncbi_dataset_package.output.dataset_package,
    output:
        ncbi_dataset_tsv="data/metadata_post_extract.tsv",
    shell:
        """
        dataformat tsv virus-genome \
            --package {input.dataset_package} \
            > {output.ncbi_dataset_tsv}
        """


rule prealign:
    input:
        sequences="data/sequences.fasta",
        reference=REFERENCE,
    output:
        prealigned_fasta="data/prealigned_sequences.fasta",
        prealigned_tsv="data/prealigned_nextclade.tsv",
    shell:
        """
        nextclade run \
        --retry-reverse-complement \
        --input-ref={input.reference} \
        --output-fasta={output.prealigned_fasta} \
        --output-tsv={output.prealigned_tsv} \
        {input.sequences}
        """


rule filter:
    message:
        """
        Filter out dropped strains, sequences with coverage 
        under {params.min_coverage} and sequences that failed prealignment
        """
    input:
        sequences="data/sequences.fasta",
        metadata=rules.format_ncbi_dataset_report.output.ncbi_dataset_tsv,
        exclude=DROPPED_STRAINS,
        prealigned_fasta="data/prealigned_sequences.fasta",
        prealigned_tsv="data/prealigned_nextclade.tsv",
    output:
        filtered_sequences="data/filtered_sequences.fasta",
        filtered_metadata="data/filtered_metadata.tsv",
    params:
        min_coverage=0.95,
    shell:
        """
        python scripts/filter_reformat.py \
            --all-sequences {input.sequences} \
            --all-metadata {input.metadata} \
            --dropped-strains {input.exclude} \
            --prealigned-sequences {input.prealigned_fasta} \
            --prealigned-tsv {input.prealigned_tsv} \
            --output-sequences {output.filtered_sequences} \
            --output-metadata {output.filtered_metadata} \
            --min-coverage {params.min_coverage}
        """

rule curate:
    message:
        """
        Curate date fields in the metadata into consistent YYYY-MM-DD format
        """
    input:
        metadata=rules.filter.output.filtered_metadata,
    output:
        curated_metadata="data/curated_metadata.tsv",
    params:
        id_column="Accession",
        date_column="date",
    shell:
        """
        augur curate format-dates \
            --metadata {input.metadata} \
            --id-column {params.id_column} \
            --date-fields {params.date_column} \
            --expected-date-formats %Y %Y-%m \
            --output-metadata {output.curated_metadata} \
	    
        """


rule align:
    message:
        """
        Aligning sequences to {input.reference}
        """
    input:
        sequences=rules.filter.output.filtered_sequences,
        reference=REFERENCE,
    output:
        alignment="data/aligned_sequence.fasta",
    shell:
        """
        nextclade run \
        --retry-reverse-complement \
        --min-seed-cover 0.2 \
        --input-ref={input.reference} \
        --output-fasta={output.alignment} \
        --include-reference=false \
        {input.sequences}
        """


rule tree:
    message:
        "Building tree"
    input:
        alignment=rules.align.output.alignment,
    output:
        tree="data/tree_raw.nwk",
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --tree-builder-args="-czb"\
            --output {output.tree}
        """


rule refine:
    message:
        """
        Refining tree. Clock rate is taken from https://github.com/nextstrain/yellow-fever/blob/main/phylogenetic/defaults/config.yaml
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
          - clock rate set to {params.clock_rate} ± {params.clock_std}
          - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation
        """
    input:
        tree=rules.tree.output.tree,
        alignment=rules.align.output.alignment,
        metadata=rules.curate.output.curated_metadata
    output:
        tree="results/tree_refined.nwk",
        node_data="results/branch_lengths.json",
    params:
        coalescent="opt",
        date_inference="marginal",
        id_column="Accession",
        root="mid_point",
        clock_rate=0.0002,
        clock_std=0.00001,
        clock_filter_iqd=4,
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --clock-rate {params.clock_rate} \
            --clock-std-dev {params.clock_std} \
            --clock-filter-iqd {params.clock_filter_iqd} \
            --year-bounds 1900 2025 \
            --root {params.root} \
            --metadata-id-columns {params.id_column}
        """


rule ancestral:
    message:
        "Reconstructing ancestral sequences and mutations"
    input:
        tree=rules.refine.output.tree,
        alignment=rules.align.output,
        reference=REFERENCE,
    output:
        node_data="results/nt_muts.json",
    params:
        inference="joint",
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --root-sequence {input.reference} \
            --inference {params.inference}
        """


rule translate:
    message:
        "Translating amino acid sequences"
    input:
        tree=rules.refine.output.tree,
        node_data=rules.ancestral.output.node_data,
        annotation=GENOME_ANNOTATION,
    output:
        node_data="results/aa_muts.json",
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.annotation} \
            --output-node-data {output.node_data} \
        """


rule clades:
    message:
        "Adding internal clade labels"
    input:
        tree=rules.refine.output.tree,
        aa_muts=rules.translate.output.node_data,
        nuc_muts=rules.ancestral.output.node_data,
        clades=CLADES,
    output:
        node_data="data/clades_raw.json",
    shell:
        """
        augur clades \
            --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data {output.node_data} 2>&1 | tee {log}
        """


rule traits:
    message:
        "Inferring ancestral traits for {params.columns!s}"
    input:
        tree=rules.refine.output.tree,
        metadata=rules.curate.output.curated_metadata,
    output:
        node_data="results/traits.json",
    params:
        columns="'Region'",
        id_column="Accession",
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output-node-data {output.node_data} \
            --columns {params.columns} \
            --confidence \
            --metadata-id-columns {params.id_column}
        """


rule export:
    message:
        "Exporting data files for for auspice"
    input:
        tree=rules.refine.output.tree,
        metadata=rules.curate.output.curated_metadata,
        clades=rules.clades.output.node_data,
        branch_lengths=rules.refine.output.node_data,
        traits=rules.traits.output.node_data,
        nt_muts=rules.ancestral.output.node_data,
        aa_muts=rules.translate.output.node_data,
        colors=COLORS,
        lat_longs=LAT_LONGS,
        auspice_config=AUSPICE_CONFIG,
    output:
        auspice_json="auspice/yellow-fever.json",
    params:
        id_column="Accession",
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} {input.clades} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_json} \
            --metadata-id-columns {params.id_column}
        """

rule generate_example_data:
    message:
        "Generating example dataset"
    input:
        sequences=rules.extract_ncbi_dataset_sequences.output.ncbi_dataset_sequences,
        clades=rules.clades.output.node_data,
        dropped_strains=DROPPED_STRAINS,
        mediocre_strains=MEDIOCRE_STRAINS,
    output:
        output_fasta="data/example_sequences.fasta",
    params:
        n_sequences=20,
        balance=True,
        random_seed=42,
    shell:
        """
        python scripts/example_data.py \
            --input-fasta {input.sequences} \
            -n {params.n_sequences} \
            --clade-membership {input.clades} \
            --seed {params.random_seed} \
            --output {output.output_fasta} \
            --bad-strains {input.dropped_strains} \
            --mediocre-strains {input.mediocre_strains}
        """

rule assemble_dataset:
    message:
        "Assembling NextClade dataset"
    input:
        tree = rules.export.output.auspice_json,
        reference = REFERENCE,
        annotation = GENOME_ANNOTATION,
        sequences = rules.generate_example_data.output.output_fasta,
        pathogen = PATHOGEN_JSON,
        readme = "README.md",
        changelog = "CHANGELOG.md",
    params:
        pathogen = "out-dataset/pathogen.json",
    output:
        tree = "out-dataset/tree.json",
        reference = "out-dataset/reference.fasta",
        annotation = "out-dataset/genome_annotation.gff3",
        sequences = "out-dataset/sequences.fasta",
        readme = "out-dataset/README.md",
        changelog = "out-dataset/CHANGELOG.md",
        dataset_zip = "dataset.zip",
    shell:
        """
        cp {input.tree} {output.tree}
        cp {input.reference} {output.reference}
        cp {input.annotation} {output.annotation}
        cp {input.sequences} {output.sequences}
        cp {input.pathogen} {params.pathogen}
        cp {input.readme} {output.readme}
        cp {input.changelog} {output.changelog}
        zip -rj dataset.zip  out-dataset/*
        """

