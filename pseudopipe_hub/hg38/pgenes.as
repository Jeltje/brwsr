table pseudopipePgenes
"Bed 12+6 file with pseudogenes and their (gene, transcript, and protein) parents."
    (
    string chrom;      "Chromosome (or contig, scaffold, etc.)"
    uint   chromStart; "Start position in chromosome"
    uint   chromEnd;   "End position in chromosome"
    string name;       "Name of item"
    uint   score;      "Score from 0-1000"
    char[1] strand;    "+ or -"
    uint thickStart;   "Start of where display should be thick (start codon)"
    uint thickEnd;     "End of where display should be thick (stop codon)"
    uint reserved;     "Used as itemRgb as of 2004-11-22"
    int blockCount;    "Number of blocks"
    int[blockCount] blockSizes;  "Comma separated list of block sizes"
    int[blockCount] chromStarts; "Start positions relative to chromStart"
    string hugo;                "Pseudogene of"
    string pgeneType;           "pseudogene, processed_pseudogene, unprocessed_pseudogene"
    string ppgene;              "Pseudopipe pseudogene ID"
    string gene;                "Ensembl pseudogene ID"
    string geneParent;          "Ensembl gene parent ID (url)"
    string proteinParent;       "Ensembl protein parent ID"
    )
