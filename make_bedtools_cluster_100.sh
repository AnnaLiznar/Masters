#/bin/bash 

cd /data1/Sebastian_sense_piRNAs/


mkdir -p ./Cluster_100

#for i in *piRNAs.bed; do 

#sort-bed  $i > ${i%%_mapped_weighted_filtered_genes_sense_piRNAs.bed}_sorted.bed
#done 

#for i in *_sorted.bed; do

#bedtools cluster -i $i -d 100 > /data1/Sebastian_sense_piRNAs/Cluster_100/${i%%sorted.bed}Cluster_100.bed

#done 


cd ./Cluster_100/ 
#mkdir -p ./Uniq
#for uniquely mapped piRNAs 

#for f in *_Cluster_100.bed; do
#cat $f | awk '{if($4 ==1)print}' > ./Uniq/${f%%.bed}_unique.bed
#done 


path=/data1/Sebastian_sense_piRNAs/Cluster_100/Uniq/
cd $path
##start python program 

#thresh
for k in *_unique.bed; do

python /data1/Sebastian_sense_piRNAs/Cluster_Sebastian.py --infile $path/$k --min_piRNAs 200 --outfile $path/${k%%.bed}_200.bed 
done 
#keep only unique -> if gene id occurs more than obce it is discarded
for k in *_200.bed; do

python /data1/Sebastian_sense_piRNAs/Cluster_Sebastian2.py --infile $path/$k --outfile $path/${k%%.bed}_unique.bed
done
#merge reps into one df
mkdir -p /Merged/

python '/data1/Sebastian_sense_piRNAs/Sebastian-3-merger.py' --infile_path '/data1/Sebastian_sense_piRNAs/Cluster_100/Uniq/' --motif 200_unique.bed --outfile '/data1/Sebastian_sense_piRNAs/Cluster_100/Uniq/Merged/small_sense_stretches.bed'


#threshold
cd /data1/Sebastian_sense_piRNAs/Cluster_100/Uniq/Merged/Mapped/
mkdir -p ./Thresh

for file in *.bed; do 
cat $file | awk '{if($13>15)print}' | awk '{ print $1,"\t"$2,"\t"$3,"\t"$4,"\t""250","\t""-" }' > /data1/Sebastian_sense_piRNAs/Cluster_100/Uniq/Merged/Mapped/Thresh/${file%%.bed}_15.bed
done

cd ./Thresh 
for file in *_15.bed; do 

tail -n +2 $file > ${file%%.bed}_noheader.bed
done 
for file in *_noheader.bed; do 

sort -u -k1,2 $file > ${file%%.bed}_sorted.bed
done 


