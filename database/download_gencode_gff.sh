wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh37_mapping/gencode.v37lift37.annotation.gff3.gz -O gencode.GRCh37lift37.annotation.gff3.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gff3.gz -O gencode.GRCh38.annotation.gff3.gz

awk '$1~/#/ || $3=="gene"' gencode.GRCh37lift37.annotation.gff3 > gencode.GRCh37lift37.annotation.only_genes.gff3
awk '$1~/#/ || $3=="gene"' gencode.GRCh38.annotation.gff3 > gencode.GRCh38.annotation.only_genes.gff3

