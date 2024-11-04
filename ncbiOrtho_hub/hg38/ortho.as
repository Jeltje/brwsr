table orthologs
"Bed 6+3 file for NCBI orthologs"
    (
    string chrom;      "Reference sequence chromosome or scaffold"
    uint   chromStart; "Start position in chromosome"
    uint   chromEnd;   "End position in chromosome"
    string name;       "Short Name of item"
    uint   score;      "Score from 0-1000"
    char[1] strand;    "+ or -"
    string hugo;       "gene symbol"
    lstring url;       "urls to orthologs"
    string ortho;      "ortholog symbol"
    )
