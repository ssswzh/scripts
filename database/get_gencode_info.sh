#!/usr/bin/bash

if [ $1 == "-h" -o $1 == "--help" -o $1 == "-help" ] || [ $# == 0 ]
then printf "\nUSAGE: \n    sh $0 annotation.only_genes.gff gene.list\n"
printf "    echo \"CYP21A1P\"|sh $0 annotation.only_genes.gff - \n\n"
printf "gene.list format: \n"
printf "gene1; gene2; gene3\n"
printf "gene4; gene5\n\n"
exit 1
fi

gff=$1
genelist=$2

awk \
    'BEGIN {
        FS="\t"; OFS=" ";
        print "Gene Info\tTotal Bases"
    }
    
    NR==FNR {
        split($9, info, ";");
        split(info[4], gene_name, "=");
        split(info[6], hgnc, "=");
        glist[gene_name[2]] = gene_name[2]":"$1":"$4"-"$5":"$7":"hgnc[2]"; ";
        dist[gene_name[2]] = $5-$4+1
    }
    
    NR>FNR {
        total_dist = 0;
        split($0, genes, ";");
        for(i=1; i<=length(genes); i++){
            gsub(/ /, "", genes[i]);
            gene = toupper(genes[i]);
            if(gene in glist){
                printf glist[gene];
                total_dist += dist[gene];
            } else {
                printf "?("gene"); "
            }
        }
        print "\t"total_dist
    }' $gff $genelist 
