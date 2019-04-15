# MiBuddy2
The script joins migec,mixcr and vdjtools softwares in a single pipeline, creating a report file with basic repertoire descriptors. 

usage:

Mibuddy2.py barcodes_file.txt -s [hsa/mmu/spalax]

+ Make sure that MiGec,MiXCR and VDJtools sowtawe is installed and the paths to the executable files are set so the above sowtware is run by simply typing migec, mixcr or vdjtools. 

+ Barcodes file should be created according to MiGec software needs (https://migec.readthedocs.io/en/latest/checkout.html)

+ In order to force MiXCR export specific chains, a sample_id in the barcodes_file.txt should include one or several of the following values ['TRA', 'TRB', 'TRG', 'TRD', 'TCR', 'IGH', 'IGK', 'IGL', 'IG'].

+ To include isotype (mixcr assemble -OseparateByC=true ...) pass an optional -ig parameter to the main command or include 'IG' in the sample_id.

+ Using '\_' in the sample_id will result in separating values in between '\_' into different metadata colomns in the metadata.txt file as well as in the report.txt output file.
