#/bin/bash

gene_ids=("COX1" "COX2" "COX3")
clustalo=("clustalo -i" "--auto --output-order=input-order -o")
mafft_auto=("mafft --auto --quiet --thread -1" ">")
mafft_parttree=("mafft --parttree --retree 2 --partsize  1000 --quiet --thread -1" ">")
mafft_linsi=("mafft --localpair --maxiterate 1000 --quiet --thread -1" ">")
muscle_auto=("muscle -quiet -in" ">")
muscle_fastdna=("muscle -maxiters 1 -diags -quiet -in" ">")
muscle_largein=("muscle -maxiters 2 -quiet -in" ">")

haplogroups=("L0" "L1" "L2" "L3" "L4" "L5" "L6"  "M"  "C"  "E"  "G"  "Q"  "Z"
              "D"  "N"  "A"  "I"  "O"  "S"  "W"  "X"  "Y" "R0" "B6"  "F"  "J"
              "P"  "K"  "B" "HV"  "T"  "R"  "H"  "V"  "U")
              # "&"  "mt-MRCA")

basename=("data/dbs/hmtDNA_2015_11_09_with_haplogroups")

for gene_id in "${gene_ids[@]}"
do
    for haplogroup in "${haplogroups[@]}"
    do
        input_filename=${basename[0]}"_"$gene_id"_"$haplogroup".fasta"

        echo
        echo $input_filename
        echo

        # ClustalO
        echo -n "Aligning with clustalo..."
        output_filename=${basename[0]}"_"$gene_id"_"$haplogroup"_aligned_with_clustalo.fasta"

        command_=${clustalo[0]}" "$input_filename" "${clustalo[1]}" "$output_filename
        eval ${command_}
        echo "done"

        # Mafft
        output_filename=${basename[0]}"_"$gene_id"_"$haplogroup"_aligned_with_mafft_auto.fasta"

        echo -n "Aligning with mafft auto..."
        command_=${mafft_auto[0]}" "$input_filename" "${mafft_auto[1]}" "$output_filename
        eval ${command_}
        echo "done"

        output_filename=${basename[0]}"_"$gene_id"_"$haplogroup"_aligned_with_mafft_linsi.fasta"

        echo -n "Aligning with mafft linsi..."
        command_=${mafft_linsi[0]}" "$input_filename" "${mafft_linsi[1]}" "$output_filename
        eval ${command_}
        echo "done"

        output_filename=${basename[0]}"_"$gene_id"_"$haplogroup"_aligned_with_mafft_parttree.fasta"

        echo -n "Aligning with mafft parttree..."
        command_=${mafft_parttree[0]}" "$input_filename" "${mafft_parttree[1]}" "$output_filename
        eval ${command_}
        echo "done"

        # Muscle
        output_filename=${basename[0]}"_"$gene_id"_"$haplogroup"_aligned_with_muscle_auto.fasta"

        echo -n "Aligning with muscle auto..."
        command_=${muscle_auto[0]}" "$input_filename" "${muscle_auto[1]}" "$output_filename
        eval ${command_}
        echo "done"

        output_filename=${basename[0]}"_"$gene_id"_"$haplogroup"_aligned_with_muscle_fastdna.fasta"

        echo -n "Aligning with muscle fastdna..."
        command_=${muscle_fastdna[0]}" "$input_filename" "${muscle_fastdna[1]}" "$output_filename
        eval ${command_}
        echo "done"

        echo -n "Aligning with muscle largein..."
        output_filename=${basename[0]}"_"$gene_id"_"$haplogroup"_aligned_with_muscle_largein.fasta"

        command_=${muscle_largein[0]}" "$input_filename" "${muscle_largein[1]}" "$output_filename
        eval ${command_}
        echo "done"
    done
done

