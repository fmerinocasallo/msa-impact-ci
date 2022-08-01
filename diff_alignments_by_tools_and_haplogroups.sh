#/bin/bash

gene_ids=("COX1" "COX2" "COX3")
tools=("clustalo"
       "mafft_auto" "mafft_linsi" "mafft_parttree"
       "muscle_auto" "muscle_fastdna" "muscle_largein")

fms=("unweighted" "weighted")
cms=("entropy" "variance")

haplogroups=("L0" "L1" "L2" "L3" "L4" "L5" "L6"  "M"  "C"  "E"  "G"  "Q"  "Z"
              "D"  "N"  "A"  "I"  "O"  "S"  "W"  "X"  "Y" "R0" "B6"  "F"  "J"
              "P"  "K"  "B" "HV"  "T"  "R"  "H"  "V"  "U")
              # "&"  "mt-MRCA")

basename=("data/reports/hmtDNA_2015_11_09_with_haplogroups" "aligned_with"
          "less_1.0_detailed.txt")

# Let's see how many columns have different conservation score regarding on the
# alignment tool used
echo "gene,1st tool,2nd tool,frequencies method,conservation method,haplogroup,#diffs"
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
                    for haplogroup in "${haplogroups[@]}"
                    do
                        file1=${basename[0]}"_"$gene_id"_"$haplogroup"_"${basename[1]}"_"${tools[$i]}"_"$fm"_"$cm"_"${basename[2]}
                        sed "/\(^>\|^There\)/d" $file1 > tmp1.txt
                        file2=${basename[0]}"_"$gene_id"_"$haplogroup"_"${basename[1]}"_"${tools[$j]}"_"$fm"_"$cm"_"${basename[2]}
                        sed "/\(^>\|^There\)/d" $file2 > tmp2.txt
                        command_=${base_command[0]}" tmp1.txt tmp2.txt | "${base_command[1]}
                        echo -n $gene_id","${tools[$i]}","${tools[$j]}","$fm","$cm","$haplogroup","
                        eval ${command_}
                    done
                done
            done
        done
    done
    echo
done

# We don't need them anymore
rm tmp1.txt tmp2.txt

