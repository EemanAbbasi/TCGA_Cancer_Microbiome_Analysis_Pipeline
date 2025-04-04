###!/bin/bash 
#!/bin/env sh
###!/usr/bin/env python

#--media M9 --mediadb M9_LB.tsv
#cd /Users/eemanabbasi/Desktop/Thesis_Project_Two/

if [ -d species_combinations ]; then
 rm -rf species_combinations
 mkdir species_combinations

fi

#if [ -d species_detailed ]; then
  #rm -rf species_detailed
  #mkdir species_detailed

#fi


current_dir=$(pwd)
#--media 'LB' --mediadb /Users/eemanabbasi/Desktop/Thesis_Project_Two/smetana_code/M9_LB.tsv
#cd /Users/eemanabbasi/Desktop/Thesis_Project_Two/nematode_models

#-m M9 --mediadb M9_LB.tsv
declare array2=();
declare df=();
declare counter=0; 

species_names=$1
path_to_save=$2
media=$3
media_comp_file=$4
echo "$media"


names=$(python3 generate.combs.py $species_names)
#names=$(python3 /Users/eemanabbasi/Desktop/Thesis_Project_Two/smetana_code/generate.combs.py $species_names)

SAVEIFS=$IFS   # Save current IFS (Internal Field Separator)
IFS=$'\n'      # Change IFS to newline char
names=($names) # split the `names` string into an array by the same name
IFS=$SAVEIFS   # Restore original IFS

for (( i=0; i<${#names[@]}; i++ ))
do
    my=`echo "${names[$i]}"   | tr ',' ' '`
    #model_name=$(printf '/Users/eemanabbasi/Desktop/Thesis_Project_Two/nematode_models/%s ' $my)
    #df=$(python3 smetana_main.py $model_name "${names[$i]}" --media 'M9' --mediadb M9_LB.tsv --output /Users/eemanabbasi/Desktop/Thesis_Project_Two/species_combinations/)
    #df=$(python3 /Users/eemanabbasi/Desktop/Thesis_Project_Two/smetana_code/smetana_main.py $model_name --media 'M9' --mediadb /Users/eemanabbasi/Desktop/Thesis_Project_Two/smetana_code/M9_LB.tsv  --output /Users/eemanabbasi/Desktop/Thesis_Project_Two/species_combinations/)

    model_name=$(printf '%s ' $my)
    echo $model_name

    if [  $media == 'complete' ]; then
        echo "in  complete"
        df=$(python3 smetana_main.py $model_name --output ${current_dir}/\${species_combinations})
        db=$(python3 smetana_main.py -d $model_name --output ${current_dir}/\${species_combinations})

    #else
        #df=$(python3 smetana_main.py $model_name --media $media --mediadb $media_comp_file --output ${current_dir}/\${species_combinations})
    fi


done


#python3 parwise_info.py ${current_dir}/\${species_combinations $path_to_save
