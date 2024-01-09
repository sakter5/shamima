cd #!/bin/bash
keyp=/home/shamima/.aspera/connect/etc/asperaweb_id_dsa.openssh
ncbiftp=anonftp@ftp.ncbi.nlm.nih.gov
cd /groups/songli_lab
for f1 in `cat SRR_files_NCBI_ftp.txt | awk '{ print $1 }'`

do
~/.aspera/connect/bin/ascp -T -l10000M -i $keyp $ncbiftp:$f1 /groups/songli_lab/RegNetRNAseq/data/SRA


done





