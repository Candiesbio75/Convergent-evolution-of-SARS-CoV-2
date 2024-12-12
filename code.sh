#!/bin/bash

####################
Virus_genome="./db/NC_045512.fasta"
NGSdata_path="./CleanData/"
lib_list=$(ls ${NGSdata_path})

# Step1 =============================
## mapping
bwa index -p ./db/NC_045512 ${Virus_genome} 
Virus_index="./db/NC_045512"

mkdir s1_bwa_mapping

for i in $(cat ../lib_list)
do

bwa_options=(
 -t 10 
 -M 
 -R \'@RG\tID:$i\tSM:$i\tLB:$i\' 
 "${Virus_index}" 
 "${NGSdata_path}""$i"/*
)

echo \
    "bwa mem "${bwa_options[@]}" | 
        samtools view --threads 10 -bF 4 -q30 | 
        samtools sort --threads 10 \
        > s1_bwa_mapping/"$i".bam" | 
tr '\n' ' ' | 
    >> s1_bwa_cmdlist
done

ParaFly -c s1_bwa_cmdlist -CPU 4 -failed_cmds s1_bwa_cmdlist.failed


# Step2 =============================
## call SNP

for i in $(ls s1_bwa_align/*.bam)
do
samtools index $i
done
ls ../step1_BWA/*.bam > listbam

mkdir step2_freebayes_CallSNP

# considering the RAM limitation, call 1000bp each time
for i in {0..28}
do

freebayes_options=(
 --haplotype-length 0 
 -m 20 
 -q 20 
 --ploidy 1 
 --pooled-continuous 
 -F 0.01 
 -r NC_045512.2:"$i"001-`expr "$i" + 1`000 
 -f "${Virus_genome}" 
 -L listbam_IFN 
)

echo \
    "freebayes "${freebayes_options[@]}" \
	> step2_freebayes_CallSNP/"$i".pooled.vcf" | 
tr '\n' ' ' | 
    >> s2_freebayes_cmdlist
done

freebayes_options=(
 --haplotype-length 0 
 -m 20 
 -q 20 
 --ploidy 1 
 --pooled-continuous 
 -F 0.01 
 -r NC_045512.2:29001-29903 
 -f "${Virus_genome}" 
 -L listbam_IFN 
)

echo \
    "freebayes "${freebayes_options[@]}" \
	> step2_freebayes_CallSNP/29.pooled.vcf" | 
tr '\n' ' ' | 
    >> s2_freebayes_cmdlist
	
ParaFly -c s2_freebayes_cmdlist -CPU 4 -failed_cmds s2_freebayes_cmdlist.failed

cat step2_freebayes_CallSNP/0.pooled..vcf | grep "^#" > step2_freebayes_CallSNP/all.pooled.vcf
cat step2_freebayes_CallSNP/*.pooled.vcf | grep -v "^#" | sort -t $'\t' -k 2 -n >> step2_freebayes_CallSNP/all.pooled.vcf


############################################################################################################################################
#########call indel
Virus_genome="./db/NC_045512.fasta"
NGSdata_path="./CleanData/"
lib_list=$(ls ${NGSdata_path})

mkdir s3_blat_CallIndel

# extract variation data of indel sites 
grep "^#" step2_freebayes_CallSNP/all.pooled.vcf > s3_blat_CallIndel/indel.pooled.vcf
grep -v "^#" step2_freebayes_CallSNP/all.pooled.vcf | awk '{if (length($4) > 1 || length($5) > 1) print $0}' \
    >> s3_blat_CallIndel/indel.pooled.vcf

# normalize indel sites
bcftools norm -f ${Virus_genome} -m - s3_blat_CallIndel/indel.pooled.vcf \
    > s3_blat_CallIndel/indel.pooled.norm.vcf
vcfallelicprimitives s3_blat_CallIndel/indel.pooled.norm.vcf --keep-info --keep-geno \
    > s3_blat_CallIndel/indel.pooled.norm.vcfallelicprimitives.vcf

# make blat reference sequences
cat s3_blat_CallIndel/indel.pooled.norm.vcfallelicprimitives.vcf | 
    grep -v "^#" | 
	cut -f2,4,5 | 
	sort -k1,1n -k2,2 -k3,3 | 
	uniq | 
	perl Find_Deletion.pl \
	> s3_blat_CallIndel/indel_blat.database
## mannually drop identical complete reference sequence duplicates of the same position
## mannually write multiallelic indels grouping table based on indel site and longest indel coverage of that site : groupToHeader.tsv

# blat
for i in $(cat lib_list)
do
seqtk seq -A ${NGSdata_path}${i}/* > s3_blat_CallIndel/${i}.reads.fasta
done
for i in $(cat lib_list)
do
blat indel_blat.database ${i}.reads.fasta -out=blast8 ${i}.blat
done

# summerize blat result
for i in $(cat lib_list)
do
echo -e ${i}"\t"100 >> s3_blat_CallIndel/library_info.txt
done

for i in $(cat lib_list)
do
grep -H -v "#" ${i}_IFN.blat | 
    awk '$3 >= 80' | 
	sed 's/\.blat:/\t/g' | 
	perl ../SpanFilter.pl | 
	awk '
	BEGIN {
        FS = "\t";
        OFS = "\t";
    }

    NR == FNR {
        group[$2] = $1;         
        next;                   
    }

    {
        key = $3;               
        if (key in group) {     
            $0 = $0 OFS group[key]; 
    }
        print;                  
    }' IFNgroupToHeader.tsv - | 
	sort | 
	awk '
    BEGIN {
        FS = "\t";               
        OFS = "\t";              
        max_val = -1;          
    }

    {
        key = $1 "|" $2 "|" $14;  
        if (!(key in seen) || $13 > values[key]) {
            seen[key] = 1;         
            values[key] = $13;     
            rows[key] = $0;        
        }
    }

    END {
        for (key in rows) {   
            print rows[key];
        }
    }
	' | 
	sort \
	> ${i}_bestfilter_IFN.blat

perl BlatSummary.pl ${i}_bestfilter_IFN.blat | \
    sort | \
	awk '
    BEGIN {
        FS="\t"; 
        OFS="\t"}

    NR==FNR { 
        group[$2]=$1;
        next} 
    
	{
	    key=$2; 
        if (key in group) {
            $0=$0 OFS group[key]; 
    }
    print;
    }' IFNgroupToHeader.tsv - \
	> ${i}_bestfilter_IFN.summary
done

cat *_bestfilter_IFN.summary > all_bestfilter_IFN.summary