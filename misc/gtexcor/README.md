# gtexcor
Correct some of the inconsistencies introduced by t2g/agar due to premature evaluation of GTEx. 

Current list of corrections:
1. wrong XS tag
2. multimappers sharing the same genomic coordinates due to overlapping transcripts/loci on opposite strands
3. overabundance of multimappers in low complexity regions
