# SRA submission files  
The two main files in this folder (`SRA_metadata-hybseqv2.txt` and `Plant.1.0-Hyb-Seqv2.xlsx`) were
used to submit fastq data to NCBI SRA (Bioproject #1223965; http://www.ncbi.nlm.nih.gov/bioproject/1223965.  

`checkEquip.R` added 3/21/2025 to check whether sequencing platform could bias the latitude results (AH); it does not. Files output:  
* `latitudeVinstrument.pdf` : two boxplots  
* `t.tests.txt` : two-sample t-tests of above question  
