# Enterovirus <...> dataset with reference <...>

| Key                  | Value                                                                 |
|----------------------|-----------------------------------------------------------------------|
| authors              | <...>, [ENPEN](https://escv.eu/european-non-polio-enterovirus-network-enpen/)                                                 |
| name                 | Enterovirus <...>                                                       |
| reference            | [<...>](https://www.ncbi.nlm.nih.gov/nuccore/<...>)         |
| workflow             | https://github.com/...                          |
| path                 | `TBD`                                                                 |
| clade definitions    |  e.g., A–C                                                                  |

## Scope of this dataset

Based on the full genome sequence, this dataset uses the <...> reference sequence (<GenBank Link>), originally isolated in <...>. It provides a framework for quality control, clade assignment, and mutation calling.

***Note:** The official RefSeq or reference sequence is substantially diverged from currently circulating strains. This is common for many enterovirus datasets, in contrast to some other virus datasets (e.g., seasonal influenza) where the reference is updated more frequently to match recent sequences.*

To address this, the dataset is *rooted* on a Static Inferred Ancestor — a phylogenetically reconstructed ancestral sequence near the tree root. This provides a stable reference point that can be used, optionally, as an alternative for mutation calling.

## Features

This dataset supports:

- Assignment of subgenotypes
- Phylogenetic placement
- Sequence quality control (QC)


## Subgenogroups of Enterovirus <...>

<...> is divided into subgenogroups A, B, C, ... 

These designations are based on the phylogenetic structure and mutations, and are widely used in molecular epidemiology, similar to subgenotype systems for other enteroviruses. Unlike influenza (H1N1, H3N2) or SARS-CoV-2, there is no universal, standardized global lineage nomenclature for enteroviruses. Naming follows conventions from published studies and surveillance practices.

## Reference types

This dataset includes several reference points used in analyses:
- *Reference:*: RefSeq or similarly established reference sequence.

- *Parent:* The nearest ancestral node of a sample in the tree, used to infer branch-specific mutations.

- *Clade founder:* The inferred ancestral node defining a clade (e.g., A2, B3). Mutations "since clade founder" describe changes that define that clade.

- *Static Inferred Ancestor:* Reconstructed ancestral sequence inferred with an outgroup, representing the likely founder of <your EV>. Serves as a stable reference.

- *Tree root:* Corresponds to the root of the tree, it may change in future updates as more data become available.

All references use the coordinate system of the Fermon sequence.

## Issues & Contact
- For questions or suggestions, please [open an issue](https://github.com/enterovirus-phylo/dataset-template-inferred-root/issues).

## What is a Nextclade dataset?

A Nextclade dataset includes the reference sequence, genome annotations, tree, clade definitions, and QC rules. Learn more in the [Nextclade documentation](https://docs.nextstrain.org/projects/nextclade/en/stable/user/datasets.html).