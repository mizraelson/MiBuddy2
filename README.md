# MiBuddy2
The script joins migec,mixcr and vdjtools softwares in a single pipeline, creating a report file with basic repertoire descriptors. 

usage:

Mibuddy2.py barcodes_file.txt -s [hsa/mmu/spalax]

+ Barcodes file should be created according to MiGec software needs (https://migec.readthedocs.io/en/latest/checkout.html)

+ In order to force MiXCR export specific chains, a sample_id in the barcodes_file.txt should include one or several of the following values ['TRA', 'TRB', 'TRG', 'TRD', 'TCR', 'IGH', 'IGK', 'IGL', 'IG'].

+ To include isotype (mixcr assemble -OseparateByC=true ...) pass an optional -ig parameter to the main command or include 'IG' in the sample_id.
