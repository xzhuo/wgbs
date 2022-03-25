# modified from https://www.danielecook.com/aggregate-fastqc-reports/
# Run this script in a directory containing zip files from fastqc. It aggregates images of each type in individual folders
# So looking across data is quick.

dir="$1"
zips=`ls $dir/*.zip`

for i in $zips; do
    unzip -d $dir -o $i &>/dev/null;
done

fastq_folders=${zips/.zip/}

rm -rf $dir/fq_aggregated # Remove aggregate folder if present
mkdir $dir/fq_aggregated

# Rename Files within each using folder name.
for folder in $fastq_folders; do
    folder=${folder%.zip}
    img_files=`ls ${folder}/Images/*png`;
    for img in $img_files; do
        img_name=$(basename "$img");
        img_name=${img_name%.*}
        new_name=${folder};
        mkdir -p $dir/fq_aggregated/${img_name};
        mv $img $dir/fq_aggregated/${img_name}/$(basename "${folder%_fastqc}").png;
    done;
done;


# Concatenate Summaries
for folder in $fastq_folders; do
    folder=${folder%.zip}
    cat ${folder}/summary.txt >> $dir/fq_aggregated/summary.txt
done;

# Concatenate Statistics
for folder in $fastq_folders; do
    folder=${folder%.zip}
    head -n 10 ${folder}/fastqc_data.txt | tail -n 7 | awk -v f=$(basename "${folder%_fastqc}") '{ print $0 "\t" f }' >> $dir/fq_aggregated/statistics.txt
    rm -rf ${folder}
done
