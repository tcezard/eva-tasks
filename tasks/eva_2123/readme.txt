
The purpose of these scripts is to help check wether a project was accessioned using chromosome
names, instead of chromosome genbank accessions.

For that, you provide a study, and the script looks in our FTP (/nfs/ftp/pub/databases/eva/PRJ*),
takes the first variant of the accessioned files, and queries it in the identifiers webservice. The
VCFs in the FTP don't reflect the chr naming used in the DB, but the webservices do.

I used the next command to loop over all the studies in https://docs.google.com/spreadsheets/d/1iNgvCaZNy7QOWid8XPXBzOptCN-WRaPAYDILfYXA6vM/edit#gid=0

cat accessioned_studies.txt | while read line ; do echo "getting $line" >&2 ; ./get_accession_study_assembly_contig_from_study_in_ftp.sh $line; done > studies_and_example_variants.txt 2>studies_and_example_variants.err

