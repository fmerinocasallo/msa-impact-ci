#/bin/bash

gene_ids=( "tRNA-Phe"  "rRNA-12S"  "tRNA-Val"  "rRNA-16S" "tRNA-Leu1"
                "ND1"  "tRNA-Ile"  "tRNA-Gln"  "tRNA-Met"       "ND2"
           "tRNA-Trp"  "tRNA-Ala"  "tRNA-Asn"  "tRNA-Cys"  "tRNA-Tyr"
               "COX1" "tRNA-Ser1"  "tRNA-Asp"      "COX2"  "tRNA-Lys"
               "ATP8"      "ATP6"      "COX3"  "tRNA-Gly"       "ND3"
           "tRNA-Arg"      "ND4L"       "ND4"  "tRNA-His" "tRNA-Ser2"
          "tRNA-Leu2"       "ND5"       "ND6"  "tRNA-Glu"      "CYTB"
           "tRNA-Thr"  "tRNA-Pro");
            # "D-loop")

# tools=("clustalo"
#        "mafft_auto" "mafft_linsi" "mafft_parttree"
#        "muscle_auto" "muscle_fastdna" "muscle_largein")
tools=("mafft_auto" "mafft_parttree")

fms=("unweighted" "weighted")
cms=("entropy" "variance")

basename=("data/reports/global/hmtDNA_2014_12_04_aligned"
          "less_1.0_detailed.txt")

# Let's see how many columns have different conservation score regarding on the
# alignment tool used
echo "gene,1st tool,2nd tool,frequencies method,conservation method,#diffs"
base_command=("diff -U 0 " "grep -c ^@")
for gene_id in "${gene_ids[@]}"
do
    for i in $(seq 0 $((${#tools[@]}-1)))
    do
        for j in $(seq $((i+1)) $((${#tools[@]}-1)))
        do
            for fm in "${fms[@]}"
            do
                for cm in "${cms[@]}"
                do
                    file1=${basename[0]}"_"${tools[$i]}"_"$gene_id"_"$fm"_"$cm"_"${basename[1]}
                    sed "/\(^>\|^There\)/d" $file1 > tmp1.txt
                    file2=${basename[0]}"_"${tools[$j]}"_"$gene_id"_"$fm"_"$cm"_"${basename[1]}
                    sed "/\(^>\|^There\)/d" $file2 > tmp2.txt
                    command_=${base_command[0]}" tmp1.txt tmp2.txt | "${base_command[1]}
                    echo -n $gene_id","${tools[$i]}","${tools[$j]}","$fm","$cm","
                    eval ${command_}
                done
            done
        done
    done
    echo
done

# We don't need them anymore
rm tmp1.txt tmp2.txt

