#### **What are the requirements of the input?**

CircNetVis can get input from

- List of circRNAs using the circRNA ID format from circbase.org (Glažar et al., 2014) (e.g., hsa\_circ\_0001946)

- List of circRNAs using the normID format as chr\_\_start\_\_end from Circall (Nguyen et al., 2021) which concatenates chromosome name and start/stop coordinates, e.g., X\_\_139865340\_\_139866824.

- List of circRNA sequences in the fasta format for example:
``` sh
>hsa_circ_0001946
GGTTTCCGATGGCACCTGTGTCAAGGTCTTCCAACAACTCCGGGTCTTCCAGCGACTTCAAGTCTTCC
AATAATCTCAAGGTCTTCCAGATAATCCTGAGCTTCCAGAAAATCCACATCTTCCAGACAATCCATGT
CTTCCGGACAATCCATGTCTTCCAAGAAGCTCCAAGTCTTCCAGTAAATCAAGTCTTCCAGCAAATCC
AGTCTTCCAGCAATTACTGGTCTTCCACCAAATCCAGATCTTCCAGGAAAATCCACGTCTTCCAGGAA
ATCCATGTCTTCCAATAATTTCAAGGTCTTCCATCAAATACAGATCTTCCAGCTAATCCATGTCTTCC
AGAAAAATCTGTGTCTTCCACCAAATCCAAGTCTTCCAGTAAATCTAGTTCTTCCAGAAAAATCTAGA
TCTTCCAGTCAATCAGTGTCTTCCAGAAAGAAATCCAGGTCTTCCAGTCAATCAGTGTCTTCCAGAAA
GAAATCCAGGTCTTCCAGTCAGTCAGTGTCTTCCAGAAAAATCTACGTCTTCCACCAAATCCAGGTCT
TCCAGTCAATCCACATCTTCCGGAAAAAATCCAGGTCTTCCAGCCAATATATGTCTTCCTGAAGATCC
ACGTCTTCCAGAAAATCCATGTCTTCCAGAAAATCCATGTCTTCCAGTAACCTCCCAGTCTTCCAGAA
AATCCACGTCTTCCCAACAATCCAAGTCTTCCGGATAATTTGGGTCTTCCTGAAAATCTACGTCTTCC
AAAAAAGCCATGTCTTCCAGAAAATCCACATCTTCCAATGGCCTCCAGGTCTTCCAGACTATCCATGT
CTTCCAGAAAATCCTTGTCTTCCCTTAAATCTATAGCTTCCAAAAAATCCGGGTCTTCCAGGAAATCC
GTGTCTTCCAGCAAGTCCACGTCTTCCAACAAAGCCATGTCTTCCAGACTATCCATGTCTTCCAGAAA
ATCCTTGTCTTCCCTCAAATCCATAGCTTCCGAAAAATCCAGGTCTTCCAGGAAATCCGTGTCTTCCA
GCAAATCCACGTCTTCCAACAAAGCCATGTCTTCCATCAAATTAATGTCTTCCAGCCTACTTGTGTCT
TCCAACAAAGGTACGTCTTCCAACAAAGGTACGTCTTCCAACAAAGGTATGTCTTCCAACAAAGGTAC
GTCTTCCAGAAAATCCACGTCTTCCAACCAAGCCATGTCTTCCAGAAAATCCACGTCTTCCAGAAAAT
ATATGTCTTCCAACTAAGCTACGTCTTCCAACAAATCCATGTCTTCCTATATCTCCAGGTCTTCCAGC
ATCTCCAGGGCTTCCAGCATCTGCTCGTCTTCCAACATCTCCACGTCTTCCAGCATCTCTGTGTCTTC
CAGCATCTTCATGTCTTCCAACAACTACCCAGTCTTCCATCAACTGGCTCAATATCCATGTCTTCCAA
CGTCTCCAGTGTGCTGATCTTCTGACATTCAGGTCTTCCAGTGTCTGCAATATCCAG
```

#### **How the circRNA sequence is generated as the input for the circRNA-miRNA prediction tools**

For each circRNA, CircNetVis generates pseudo-sequence by using Circall-simulator (Nguyen et al., 2021). Given a sequence of a circRNA, last 27 bases (max miRNA length minus one) are added to the beginning of the sequence to capture the back-splicing junction region of the circRNA. The added bases can sometimes produce duplidated binding sites in the prediction of circRNA-miRNA interactions, however this rarely happens.

#### **Settings for circRNA-miRNA interactions**

Three commonly used prediction tools can be selected to discover circRNA-miRNA interactions: TargetScan, RNAhybrid and miRanda. Users can utilize the results from a single tool or intersection of results of different tools.

-  Filtering for TargetScan

Users can set the minimum number of binding sites of the interactions for each site type including 6mer, 7mer-1a, 7mer-m8, 8mer-1a. The interactions which were sucessful validated by wetlab work usually have higher numbers of biding sites.

-  Filtering for RNAhybrid

For RNAhybrid, users can filer the its results by p-value and/or minimum of free energy (mfe, negative value). In many studies, mfe no greater than -18 is preferably used as the filter.

-  Filtering for miRanda

For miRanda, the value of Max_Score from the output of the tool is usually used for the filtering. From the literature, Max_Score no less than 140 or more is often applied.

**Details of the command lines to run the prediction tools**

-   For TargetScan:

``` sh
perl Tools/targetscan.pl $miRNA_input $circRNA_input targetscan_out.txt
```

where the miRNA_input file contain data of all mircRNAs from the miRbase.org database (Kozomara et al., 2019) and the circRNA_input file contain the data of the circRNAs of interest. The results provided in targetscan_out.txt are input into CircNetVis.

-   For RNAhybrid:

``` sh
Tools/RNAhybrid -m $maxCircRNALen -n $maxmiRNALen -b 1 -c -s 3utr_human -p 1.0 -q $miRNA_fa -t $circRNA_fa > RNAhybrid_out.txt
```

where maxCircRNALen is the maximum sequence length of the circRNAs of interest (plus 1); maxmiRNALen is the maximum sequence length of the miRNAs (plus 1); miRNA_fa file contains the fasta sequence of all miRNAs from the miRbase database; circRNA_fa file consists of the fasta sequence of the circRNAs. The results provided in RNAhybrid_out.txt are input into CircNetVis. The value of p ("-p 1.0") indicates that we all output with p-value \<= 1.0; the value of b ("-b 1") means we keep only one hit per target.

-   For miRanda:

``` sh
Tools/miranda $miRNA_fa $circRNA_fa -sc 100 -quiet | grep '>>hsa' > miranda_out.txt
```

where miRNA_fa file contains the fasta sequence of all miRNAs from the miRbase database; circRNA_fa file consists of the fasta sequence of the circRNAs. The results provided in miranda_out.txt are input into CircNetVis. The value of sc ("-sc 100") indicates that we keep only interactions with the score of at least 100.

#### **Settings for miRNA-mRNA interactions**

The interactions of miRNA-mRNA are collected directly from the TargetScan database version 72 for visualization. Two common parameters are used for filtering: weighted-context score (the lower is better) and weighted-context score percentile (the higher is better); please refer to the website of TargetScan for the details of the parameters. Users also can set the number of genes displaying in the interaction networks.

#### **Settings for circRNA-RBP interactions**

The interactions of circRNA-RBP are collected from the CircInteractome database (Dudekula et al., 2016). Users filter the interactions by the number of binding sites and limit the sets of RBP displaying in the interaction networks.

###### **References**
-   Dudekula,D.B. et al. (2016) CircInteractome: A web tool for exploring circular RNAs and their interacting proteins and microRNAs. RNA Biol., 13, 34--42.
-   Glažar,P. et al. (2014) circBase: a database for circular RNAs. RNA, 20, 1666--1670.
-   Kozomara,A. et al. (2019) miRBase: from microRNA sequences to function. Nucleic Acids Res., 47, D155--D162.
-   Nguyen,D. et al. (2021) Circall: fast and accurate methodology for discovery of circular RNAs from paired-end RNA-sequencing data. BMC Bioinformatics.
-   Jassal,B. et al. (2020) The reactome pathway knowledgebase. Nucleic Acids Res., 48, D498--D503.


