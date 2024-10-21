table spliceAI
"Bed 9+4 file with Ensembl Gene IDs and ABsplice values per tissue."
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
    float AIscore;       "spliceAI score"
    string spliceType; "donor_gain, donor_loss, acceptor_gain, or acceptor_loss"
    string relativePos;   "Relative location of donor or acceptor affected by this variant"
    lstring gene;      "Gene ID" 
    )
