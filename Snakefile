include: "scripts/workflow_messages.snkm"

# -----------------------------------------------------------------------------
# Parameters and Settings
# -----------------------------------------------------------------------------
# Define general parameters, filtering thresholds, and workflow options.

REFERENCE_ACCESSION =   "<accession>"   # define the reference sequence
TAXON_ID =              <...>           # define the taxon id of your virus 
GENES =                 ["VP4", "VP2", "VP3", "VP1", "2A", "2B", "2C", "3A", "3B", "3C", "3D"] # specify the genes that will be included. The names must match the annotation & reference gene/product names
ALLOWED_DIVERGENCE =    "3000"          # TODO: lower this threshold to exclude outliers
MIN_DATE =              "1960-01-01"    # all sequences collected beforehand will be excluded
MIN_LENGTH =            "6000"          # is 6000 bp for whole genome build on Nextstrain
MAX_SEQS =              "10000"         # TODO: set lower to subsample the tree
ROOTING =               "ancestral_sequence"  # mid_point, outgroup, reference, ancestral sequence
ID_FIELD=               "accession"     # either accession or strain, used for meta-id-column in augur

FETCH_SEQUENCES = True              # whether to fetch sequences from NCBI Virus via ingest workflow
STATIC_ANCESTRAL_INFERRENCE = True  # whether to use the static inferred ancestral sequence
INFERRENCE_RERUN = False            # whether to rerun the inference of the ancestral sequence worfkflow (inferred-root)

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
# Define all input and resource file locations.

SEQUENCES =             "data/sequences.fasta"              # Input nucleotide sequences in FASTA format.
METADATA =              "data/metadata.tsv"                 # Associated sequence metadata in TSV format.

INFERRED_SEQ_PATH = "results/sequences_with_ancestral.fasta" if STATIC_ANCESTRAL_INFERRENCE else SEQUENCES
INFERRED_META_PATH = "results/metadata_with_ancestral.tsv" if STATIC_ANCESTRAL_INFERRENCE else "results/metadata.tsv"

GFF_PATH =              "dataset/genome_annotation.gff3"    # Reference genome annotation in GFF3 format.
PATHOGEN_JSON =         "dataset/pathogen.json"             # Pathogen definition file for Nextclade QC and alignment.
README_PATH =           "dataset/README.md"                 # Dataset description and usage notes.
CHANGELOG_PATH =        "dataset/CHANGELOG.md"              # Log of dataset releases and version history.
REFERENCE_PATH =        "dataset/reference.fasta"           # Reference genome sequence in FASTA format.

AUSPICE_CONFIG =        "resources/auspice_config.json"     # Configuration for Auspice visualization.
EXCLUDE =               "resources/exclude.txt"             # List of sequences to exclude from the build.
CLADES =                "resources/clades.tsv"              # Clade definitions for Nextclade/Nextstrain builds.
ACCESSION_STRAIN =      "resources/accession_strain.tsv"    # Optional mapping of accession to corrected strain names
INCLUDE_EXAMPLES =      "resources/include_examples.txt"    # List of sequences to include in example sequences (rule subsample_example_sequences)
COLORS =                "resources/colors.tsv"              # Color assignments for metadata fields.
COLORS_SCHEMES =        "resources/color_schemes.tsv"       # Preset color schemes for Auspice.
GENBANK_PATH =          "resources/reference.gbk"           # Reference genome in GenBank format (for gene mapping)
INFERRED_ANCESTOR =     "resources/inferred-root.fasta"     # Inferred ancestral sequence (used as alternative reference)
# -----------------------------------------------------------------------------

configfile: PATHOGEN_JSON

# -----------------------------------------------------------------------------
# Workflow Rules
# -----------------------------------------------------------------------------
# Define the main workflow steps and dependencies.
# Each rule represents a key stage in the dataset build.

rule all:
    input:
        auspice = "results/auspice.json",
        augur_jsons = "test_out/",
        data = "dataset.zip",
        seqs = "results/example_sequences.fasta",
        **({"root": INFERRED_ANCESTOR} if STATIC_ANCESTRAL_INFERRENCE else {})


if FETCH_SEQUENCES == True:
    rule fetch:
        input:
            dir = "ingest"
        output:
            sequences=SEQUENCES,
            metadata=METADATA
        shell:
            """
            cd {input.dir} 
            snakemake --cores 9 all
            cd ../
            """

rule curate:
    message:
        """
        Cleaning up metadata with augur merge & augur curate
        """
    input:
        meta = METADATA,  # Path to input metadata file
        strains = ACCESSION_STRAIN  # Strain - accession lookup table
    params:
        strain_id_field = ID_FIELD,
    output:
        metadata = "results/metadata.tsv",  # Final output file for publications metadata
    shell:
        """
        augur merge --metadata metadata={input.meta} strains={input.strains}\
            --metadata-id-columns {params.strain_id_field} \
            --output-metadata metadata.tmp
        augur curate normalize-strings \
            --metadata metadata.tmp \
            --id-column {params.strain_id_field} \
            --output-metadata {output.metadata}

        rm metadata.tmp
        """

rule add_reference_to_include:
    """
    Create an include file for augur filter
    """
    input:
        "resources/include.txt",
    output:
        "results/include.txt",
    shell:
        """
        cat {input} >> {output}
        echo "{REFERENCE_ACCESSION}" >> {output}
        echo ancestral_sequence >> {output}
        """

if STATIC_ANCESTRAL_INFERRENCE and INFERRENCE_RERUN:
    rule static_inferrence:
        message:
            """
            Running "inferred-root" snakefile for inference of the ancestral root. 
            This reference will be included in the Nextclade reference tree.
            WARNING: This will overwrite your {output.seq} & {output.meta} files!
            """
        input:
            dir = "inferred-root",
            dataset_path = "dataset",
            meta = rules.curate.output.metadata,
            seq = SEQUENCES,
            meta_ancestral = "resources/static_inferred_metadata.tsv",
            include = "results/include.txt"
        params:
            strain_id_field = ID_FIELD,
        output:
            inref = INFERRED_ANCESTOR,
            seq = INFERRED_SEQ_PATH,
            meta = INFERRED_META_PATH,
        threads: workflow.cores
        shell:
            r"""
            set -euo pipefail

            echo "Cleaning previous results..."
            rm -rf {input.dir}/results/* {input.dir}/resources/inferred-root.fasta

            echo "Running inferred-root workflow..."
            cd {input.dir}
            snakemake --cores {threads} all_sub
            cd - > /dev/null

            echo "Combining sequences with ancestral root..."
            cat {input.seq} {output.inref} > {output.seq}

            echo "Merging metadata..."
            augur merge \
                --metadata metadata={input.meta} ancestral={input.meta_ancestral} \
                --metadata-id-columns {params.strain_id_field} \
                --output-metadata {output.meta}

            echo "Static ancestral inference completed successfully!"
            """

rule index_sequences:
    message:
        """
        Creating an index of sequence composition for filtering
        """
    input:
        sequences = INFERRED_SEQ_PATH,
    output:
        sequence_index = "results/sequence_index.tsv"
    shell:
        """
        augur index \
            --sequences {input.sequences} \
            --output {output.sequence_index}
        """

rule filter:
    """
    Exclude sequences from before {MIN_DATE} and subsample to {MAX_SEQS} sequences.
    Only take sequences longer than {MIN_LENGTH}
    """
    input:
        sequences = INFERRED_SEQ_PATH,
        sequence_index = rules.index_sequences.output.sequence_index,
        metadata = INFERRED_META_PATH,
        include = rules.add_reference_to_include.output,
    output:
        filtered_sequences = "results/filtered_sequences_raw.fasta",
        filtered_metadata = "results/filtered_metadata_raw.tsv",
    params: 
        min_date="" if MIN_DATE == "" else "--min-date " + MIN_DATE,
        min_length="" if MIN_LENGTH == "" else "--min-length " + MIN_LENGTH,
        max_seqs=MAX_SEQS,
        # categories = "country year", #TODO: add subsampling per category?
        strain_id_field = ID_FIELD,
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            {params.min_length} \
            {params.min_date} \
            --include {input.include} \
            --subsample-max-sequences {params.max_seqs} \
            --output-sequences {output.filtered_sequences} \
            --output-metadata {output.filtered_metadata}
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference} using Nextclade3.
        """
    input:
        sequences = rules.filter.output.filtered_sequences,
        reference = REFERENCE_PATH,
        annotation = GFF_PATH,
    output:
        alignment = "results/aligned.fasta",
        tsv = "results/nextclade.tsv",
    params:
        translation_template = lambda w: "results/translations/cds_{cds}.translation.fasta",
        penalty_gap_extend = config["alignmentParams"]["penalityGapExtend"],
        penalty_gap_open = config["alignmentParams"]["penaltyGapOpen"],
        penalty_gap_open_in_frame = config["alignmentParams"]["penaltyGapOpenInFrame"],
        penalty_gap_open_out_of_frame = config["alignmentParams"]["penaltyGapOpenOutOfFrame"],
        kmer_length = config["alignmentParams"]["kmerLength"],
        kmer_distance = config["alignmentParams"]["kmerDistance"],
        min_match_length = config["alignmentParams"]["minMatchLength"],
        allowed_mismatches = config["alignmentParams"]["allowedMismatches"],
        min_length = config["alignmentParams"]["minLength"],
        gap_alignment_side = config["alignmentParams"]["gapAlignmentSide"],  
        min_seed_cover = config["alignmentParams"]["minSeedCover"],
    shell:
        """
        nextclade3 run \
        -j {threads} \
        {input.sequences} \
        --input-ref {input.reference} \
        --input-annotation {input.annotation} \
        --alignment-preset high-diversity \
        --penalty-gap-open {params.penalty_gap_open} \
        --penalty-gap-extend {params.penalty_gap_extend} \
        --penalty-gap-open-in-frame {params.penalty_gap_open_in_frame} \
        --penalty-gap-open-out-of-frame {params.penalty_gap_open_out_of_frame} \
        --gap-alignment-side {params.gap_alignment_side} \
        --kmer-length {params.kmer_length} \
        --kmer-distance {params.kmer_distance} \
        --min-match-length {params.min_match_length} \
        --allowed-mismatches {params.allowed_mismatches} \
        --min-seed-cover {params.min_seed_cover} \
        --min-length {params.min_length} \
        --max-alignment-attempts 5 \
        --include-reference false \
        --output-tsv {output.tsv} \
        --output-translations {params.translation_template} \
        --output-fasta {output.alignment} 
        """


rule get_outliers:
    """
    Automatically identify sequences with >{ALLOWED_DIVERGENCE} substitutions
    (likely to be sequencing errors or low quality/misannotated sequences) and put them in outliers.txt
    """
    input:
        nextclade = rules.align.output.tsv,
    output:
        outliers = "results/outliers.txt",
        tmp = "tmp/outliers.txt",
    params:
        allowed_divergence = lambda w: ALLOWED_DIVERGENCE,
    shell:
        """
        tsv-filter -H -v --is-numeric totalSubstitutions {input.nextclade} \
        > {output.tmp}
        tsv-filter -H \
            --is-numeric totalSubstitutions \
            --gt totalSubstitutions:{params.allowed_divergence} \
            {input.nextclade} \
        | tail -n +2 >> {output.tmp}
        cat {output.tmp} \
        | tsv-select -H -f seqName \
        | tail -n +2 > {output.outliers}
        """


rule exclude:
    """
    Rule to allow for manual and automatic exclusion of sequences
    without triggering a new subsampling that could
    surface new bad sequences resulting in an infinite loop
    """
    input:
        sequences = rules.align.output.alignment,
        sequence_index = rules.index_sequences.output.sequence_index,
        metadata = INFERRED_META_PATH,
        exclude = EXCLUDE,
        outliers = rules.get_outliers.output.outliers,
        example = INCLUDE_EXAMPLES,

    params:
        strain_id_field = ID_FIELD,
    output:
        filtered_sequences = "results/filtered_aligned.fasta",
        filtered_metadata = "results/filtered_metadata.tsv",
        strains = "results/tree_strains.txt",
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --exclude {input.exclude} {input.outliers} {input.example} \
            --output-sequences {output.filtered_sequences} \
            --output-metadata {output.filtered_metadata} \
            --output-strains {output.strains}
        """

rule tree:
    message:
        """
        Creating a maximum likelihood tree
        """
    input:
        alignment = rules.exclude.output.filtered_sequences,
    output:
        tree = "results/tree_raw.nwk",
    threads: workflow.cores
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --nthreads {threads}\
            --output {output.tree} \
        """

rule refine:
    input:
        tree=rules.tree.output.tree,
        alignment=rules.exclude.output.filtered_sequences,
    output:
        tree="results/tree.nwk",
        node_data="results/branch_lengths.json",
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --root {ROOTING} \
            --keep-polytomies \
            --divergence-unit mutations-per-site \
            --output-node-data {output.node_data} \
            --output-tree {output.tree}
        """


rule ancestral:
    input:
        tree=rules.refine.output.tree,
        alignment=rules.exclude.output.filtered_sequences,
        annotation=GENBANK_PATH,
    output:
        node_data="results/muts.json",
        ancestral_sequences="results/ancestral_sequences.fasta",
    params:
        translation_template=r"results/translations/cds_%GENE.translation.fasta",
        output_translation_template=r"results/translations/cds_%GENE.ancestral.fasta",
        genes=" ".join(GENES),
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --annotation {input.annotation} \
            --root-sequence {input.annotation} \
            --genes {params.genes} \
            --translations {params.translation_template} \
            --output-node-data {output.node_data} \
            --output-translations {params.output_translation_template}\
            --output-sequences {output.ancestral_sequences}
        """

rule clades:
    input:
        tree = rules.refine.output.tree,
        mutations = rules.ancestral.output.node_data,
        clades = CLADES
    output:
        json = "results/clades.json",
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.mutations} \
            --clades {input.clades} \
            --output-node-data {output.json}
        """

rule get_dates:
    """Create ordering for color assignment"""
    input:
        metadata = rules.exclude.output.filtered_metadata
    output:
        ordering = "results/color_ordering.tsv"
    run:
        import pandas as pd
        column = "date"
        meta = pd.read_csv(input.metadata, delimiter='\t')

        if column not in meta.columns:
            print(f"The column '{column}' does not exist in the file.")
            sys.exit(1)

        deflist = meta[column].dropna().tolist()
        # Store unique values (ordered)
        deflist = sorted(set(deflist))
        if "XXXX-XX-XX" in deflist:
            deflist.remove("XXXX-XX-XX")

        result_df = pd.DataFrame({
            'column': ['date'] * len(deflist),
            'value': deflist
        })

        result_df.to_csv(output.ordering, sep='\t', index=False, header=False)
        

rule colors:
    """Assign colors based on ordering"""
    input:
        ordering=rules.get_dates.output.ordering,
        color_schemes=COLORS_SCHEMES,
        colors=COLORS,
    output:
        colors="results/colors_dates.tsv",
        final_colors="results/final_colors.tsv"
    shell:
        """
        python3 scripts/assign-colors.py \
            --ordering {input.ordering} \
            --color-schemes {input.color_schemes} \
            --output {output.colors}

        echo -e '\ndate\tXXXX-XX-XX\t#a6acaf' >> {output.colors}

        cat {output.colors} {input.colors} >> {output.final_colors}
        """

rule export: 
    input:
        tree = rules.refine.output.tree,
        metadata = rules.exclude.output.filtered_metadata,
        mutations = rules.ancestral.output.node_data,
        branch_lengths = rules.refine.output.node_data,
        clades = rules.clades.output.json, # dummy_clades if not set yet
        auspice_config = AUSPICE_CONFIG,
        colors = rules.colors.output.final_colors
    params:
        strain_id_field = ID_FIELD,
        fields="region country date",
    output:
        auspice = "results/auspice.json",
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --auspice-config {input.auspice_config} \
            --color-by-metadata {params.fields} \
            --colors {input.colors} \
            --include-root-sequence \
            --node-data {input.mutations} {input.branch_lengths} {input.clades} \
            --output {output.auspice}
        """


rule subsample_example_sequences:
    input:
        all_sequences = INFERRED_SEQ_PATH,
        metadata = INFERRED_META_PATH,
        incl_examples = INCLUDE_EXAMPLES,
        exclude = EXCLUDE,
        outliers = rules.get_outliers.output.outliers,
    output:
        example_sequences = "results/example_sequences.fasta",
    params:
        strain_id_field = ID_FIELD,
    shell:
        """
        augur filter \
            --sequences {input.all_sequences} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --min-date 2015 --group-by year --subsample-max-sequences 30  \
            --include {input.incl_examples} \
            --exclude {input.exclude} {input.outliers} \
            --exclude-ambiguous-dates-by year \
            --probabilistic-sampling \
            --output-sequences {output.example_sequences}
        """

rule assemble_dataset:
    input:
        tree = rules.export.output.auspice,
        reference = REFERENCE_PATH,
        annotation = GFF_PATH,
        sequences = rules.subsample_example_sequences.output.example_sequences,
        pathogen = PATHOGEN_JSON,
        readme = README_PATH,
        changelog = CHANGELOG_PATH,
    output:
        tree = "out-dataset/tree.json",
        reference = "out-dataset/reference.fasta",
        annotation = "out-dataset/genome_annotation.gff3",
        sequences = "out-dataset/sequences.fasta",
        pathogen = "out-dataset/pathogen.json",
        readme = "out-dataset/README.md",
        changelog = "out-dataset/CHANGELOG.md",
        dataset_zip = "dataset.zip",
    shell:
        """
        cp {input.tree} {output.tree}
        cp {input.reference} {output.reference}
        cp {input.annotation} {output.annotation}
        cp {input.sequences} {output.sequences}
        cp {input.pathogen} {output.pathogen}
        cp {input.readme} {output.readme}
        cp {input.changelog} {output.changelog}
        zip -rj dataset.zip  out-dataset/*
        """


rule test:
    input:
        dataset = rules.assemble_dataset.output.dataset_zip,
        sequences = rules.assemble_dataset.output.sequences,
    output:
        output = directory("test_out"),
    shell:
        """
        nextclade3 run \
            --input-dataset {input.dataset} \
            --output-all {output.output} \
            {input.sequences}
        """

rule mutLabels:
    input:
        table = "results/nextclade.tsv",
        clade = "results/clades_metadata.tsv",
        json = PATHOGEN_JSON,
    params:
        min_proportion = 0.2,
        high_threshold_proportion = 0.70,
        clades_high_threshold = ["clade1","clade2","..."],
        clades_to_drop = ["unassigned"],
    output:
        clade_meta = "results/clades_mut_metadata.tsv",
        properties = "results/virus_properties.json",
        json = "out-dataset/pathogen.json"
    shell:
        """
        augur merge \
            --metadata meta={input.table} clade={input.clade} \
            --metadata-id-columns meta=seqName clade=accession \
            --output-metadata {output.clade_meta}

        python3 scripts/generate_virus_properties.py \
            --clade_meta {output.clade_meta} \
            --properties {output.properties} \
            --min-prop {params.min_proportion} \
            --high-min-prop {params.high_threshold_proportion} \
            --high-prop-clades "{params.clades_high_threshold}" \
            --exclude-clades "{params.clades_to_drop}"


        jq --slurpfile v {output.properties} \
           '.mutLabels.nucMutLabelMap = $v[0].nucMutLabelMap |
            .mutLabels.nucMutLabelMapReverse = $v[0].nucMutLabelMapReverse' \
           {input.json} > {output.json}

        zip -rj dataset.zip  out-dataset/*
        """


rule clean:
    shell:
        """
        rm -r results out-dataset test_out dataset.zip tmp
        rm ingest/data/* data/*
        rm resources/inferred-root.fasta
        rm -r inferred-root/results/* inferred-root/resources/*
        """
