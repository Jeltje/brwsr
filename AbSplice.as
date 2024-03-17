table abSplice
"Bed 9+5 file with Ensembl Gene IDs and ABsplice values per tissue."
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
    string ENSGid;     "Ensembl Gene ID"
    string hugoId;     "hugo Gene ID"
    float spliceABscore; "AbSplice highest score for this position"
    lstring maxScore;   "All tissues containing the highest score"
    lstring tissues;   "All 49 GTEX tissues with ABSplice value (empty if none were provided)"
    )
