# ABE-APOBEC
The progress of RBP and Base editor


## bactools Vcf files progress
Run vcf_process.py on wsl or windows to process vcf files.
It need four parameters:
   s, # strand: 1=FWD/neg, 2=REV/pos
   input, # Input VCF directory (including ProQ-rABE_FWD/REV.vcf[.gz])
   output, # Output Directory
   samples # Optional: list of sample names; if not provided, they will be automatically inferred from the VCF header
 python vcf_process.py [S] \
   --input \
   --output
 
