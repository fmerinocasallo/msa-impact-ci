#!/bin/bash

# gene_ids=("COX1" "COX2" "COX3")
# gene_ids=("D-loop" "tRNA-Phe" "tRNA-Pro")

tools=("clustalo"
       "mafft_parttree")
       # "mafft_parttree_retree-1"
       # "mafft_parttree_retree-2_partsize-1000")

fms=("unweighted" "weighted")
cms=("entropy" "variance")

alignment_path="/home/franxesk/Documents/alignments/2016-07-12-bis/data/dbs/"
report_path="/home/franxesk/Documents/alignments/2016-07-12-bis/data/reports/detailed/"
gen_report_path="/home/franxesk/Documents/alignments/2016-07-12-bis/conservation-index/"

# basename=("hmtDNA_2015_11_09_processables_aligned_with_"
basename=("hmtDNA_2015_11_09_without_D-loop_processables_aligned_with_"
          "less_1.00_detailed.txt")

base_command_report="python3.4 "$gen_report_path"generate_report.py" 
base_command_diff="python diff_alignments.py"

# Let's see how many columns have different conservation score regarding on the
# alignment tool used
# echo "gene,1st tool,2nd tool,frequencies method,conservation method,haplogroup,#diffs"
echo "1st tool,2nd tool,frequencies method,conservation method,min_diff,diff gap,diff residue,diff frequency"
# for gene_id in "${gene_ids[@]}"
# do
for i in $(seq 0 $((${#tools[@]}-1)))
do
    for j in $(seq $((i+1)) $((${#tools[@]}-1)))
    do
        for fm in "${fms[@]}"
        do
            for cm in "${cms[@]}"
            do
                alignment1=$alignment_path${basename[0]}${tools[$i]}".fasta"
                # file1=${basename[0]}$gene_id"_"${basename[1]}"_"${tools[$i]}"_"$fm"_"$cm"_"${basename[2]}
                report1=$report_path${basename[0]}${tools[$i]}"_"$fm"_"$cm"_"${basename[1]}
                if [ ! -f $report1 ]
                then
                    command_report=$base_command_report" -if "$alignment1
                    command_report=$command_report" -st dna -rt detailed -fm "$fm
                    command_report=$command_report" -cm "$cm" -co less -th 1.00"
                    command_report=$command_report" -od "$report_path
                    eval $command_report
                fi

                alignment2=$alignment_path${basename[0]}${tools[$j]}".fasta"
                # file2=${basename[0]}$gene_id"_"${basename[1]}"_"${tools[$j]}"_"$fm"_"$cm"_"${basename[2]}
                report2=$report_path${basename[0]}${tools[$j]}"_"$fm"_"$cm"_"${basename[1]}
                if [ ! -f $report2 ]
                then
                    command_report=$base_command_report" -if "$alignment2
                    command_report=$command_report" -st dna -rt detailed -fm "$fm
                    command_report=$command_report" -cm "$cm" -co less -th 1.00"
                    command_report=$command_report" -od "$report_path
                    eval $command_report
                fi

                for min_diff_pow in $(seq -1 1)
                do
                    min_diff=`bc -l <<< "scale=4; 10^${min_diff_pow}"`
                    # diff_="data/diffs/complete_"$fm"_"$cm"_"$min_diff".txt"
                    diff_="data/diffs/without_D-loop_"$fm"_"$cm"_"$min_diff".txt"

                    diff_command=$base_command_diff" -f1 "$report1" -f2 "$report2
                    diff_command=$diff_command" -df "$min_diff" > "$diff_
                    eval ${diff_command}

                    # echo -n $gene_id","${tools[$i]}","${tools[$j]}","$fm","$cm","$haplogroup","
                    echo -n ${tools[$i]}","${tools[$j]}","$fm","$cm","$min_diff","
                    gaps=`cat $diff_ | grep "gap" | wc -l`
                    echo -n $gaps","
                    residues=`cat $diff_ | grep "res" | wc -l`
                    echo -n $residues","
                    frequencies=`cat $diff_ | grep "freq" | wc -l`
                    echo $frequencies
                done
            done
        done
    echo
    done
done

