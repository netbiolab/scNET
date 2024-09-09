#!/usr/bin/bash
#scNET.sh [input] [filename] [outdir] [BScore] [SVcore] [BSreuse] [SVreuse] [package path]

## part1
input=$1
filename=$2
outdir=$3
BScore=$4
SVcore=$5
BSreuse=$6
SVreuse=$7
packagepath=$8


cd $packagepath

start=$(date)
start_sec=$(date +%s)
/bin/time -vv ./functions/scFGN_1_BS.R -i $input -o $filename -d $outdir -n $BScore -r $BSreuse -gc example/input/HNv3_coverage_genes.rds -p $packagepath 2> BS_mem_tmp &
/bin/time -vv ./functions/scFGN_1_SV.R -i $input -o $filename -d $outdir -n $SVcore -r $SVreuse -c 99 -gc example/input/HNv3_coverage_genes.rds -p $packagepath 2> SV_mem_tmp &
/bin/time -vv ./functions/scFGN_1_SC.R -i $input -o $filename -d $outdir -gc example/input/HNv3_coverage_genes.rds 2> SC_mem_tmp &
wait
./functions/scFGN_1.5_Bin_automation.R -o $filename -d $outdir -bs ${filename}_BS_PCCnet_possorted.binlls -sv ${filename}_SV_PCCnet_possorted.binlls -sc ${filename}_SC_PCCnet_possorted.binlls
BS_bin_suggest=$(cat BS_bin_tmp)
SV_bin_suggest=$(cat SV_bin_tmp)
SC_bin_suggest=$(cat SC_bin_tmp)
./functions/scFGN_2_fit_curve_script.R -bl ${filename}_BS_PCCnet_possorted.binlls -b ${filename}_BS_PCCnet_possorted -sl ${filename}_SV_PCCnet_possorted.binlls -s ${filename}_SV_PCCnet_possorted -scl ${filename}_SC_PCCnet_possorted.binlls -sc ${filename}_SC_PCCnet_possorted -o $filename -d $outdir || echo "error for scFNG_2" >> ${filename}_log.txt 

BS_memory=$(echo "$(cat BS_mem_tmp | grep -o 'Maximum resident set size (kbytes):\s[0-9]*' | awk -F '[ ]' '{print $6}')*0.00000095367432"|bc)
SV_memory=$(echo "$(cat SV_mem_tmp | grep -o 'Maximum resident set size (kbytes):\s[0-9]*' | awk -F '[ ]' '{print $6}')*0.00000095367432"|bc)

rm BS_bin_tmp
rm SV_bin_tmp
rm SC_bin_tmp
rm BS_mem_tmp
rm SV_mem_tmp
rm SC_mem_tmp
echo "completed Bigscale, SAVER and SuperCell"


./functions/scFGN_3_overlap.R -b ${filename}_BS_fit -s ${filename}_SV_fit -c ${filename}_SC_fit -d $outdir -o $filename || echo "error for scFNG_3" >> ${filename}_log.txt  
./functions/scFGN_4_integration_1_jointable.py ${filename}_joined_net ${filename}_BS_fit ${filename}_SV_fit ${filename}_SC_fit || echo "error for scFNG_4-1" >> ${filename}_log.txt 
./functions/scFGN_4_integration_2_optimize_weight_inwsum.py "./example/input/gold_standard_symbol_HNv3" ${filename}_joined_net ${filename}_wsum_report || echo "error for scFNG_4-2" >> ${filename}_log.txt 
report=$(./functions/scFGN_4_integration_3_report_reader.py ${filename}_wsum_report)
weight=$(echo $report | grep -o -e 'Weight\s[0-9]*' -e 'Weight\sPositive\sinfinity' | awk -F '[ ]' '{print $2}')
cutoff=$(echo $report | grep -o 'Threshold\s[0-9]*\.[0-9]*' | awk -F '[ ]' '{print $2}')
coverage=$(cat "final_coverage.txt")
rm "final_coverage.txt"
if [ ${weight} = "Positive" ]
then
    ./functions/scFGN_4_integration_4_summation.py ${filename}_joined_net ${filename}_WS_net max $cutoff || echo "error for scFGN_4-4" >> ${filename}_log.txt
elif [ ${weight} -eq 0 ]
then
    ./functions/scFGN_4_integration_4_summation.py ${filename}_joined_net ${filename}_WS_net naivesum $cutoff || echo "error for scFGN_4-4" >> ${filename}_log.txt
else
    ./functions/scFGN_4_integration_4_summation.py ${filename}_joined_net ${filename}_WS_net wsum $cutoff -w $weight || echo "error for scFNG_4-4" >> ${filename}_log.txt 
fi
./functions/scFGN_5_WS_final_fit.R -i ${filename}_WS_net -d $outdir -o $filename -p $packagepath || echo "error for scFNG_5" >> ${filename}_log.txt 
end=$(date)
end_sec=$(date +%s)
time_elap=$(($end_sec - $start_sec))
time_elap_min=$(echo "scale=2;$time_elap/60"|bc)

if test -f "${filename}_WS_net_final" ; then
    echo -e "\n\nAUTO2_log\nStart Time : ${start}\nEnd Time : ${end}\nTime Elapsed : ${time_elap_min} mins\nBigScale Bin : ${BSbin}\nSAVER Bin : ${SVbin}\nCutoff : ${cutoff}\nWeight : ${weight}\nFinalCoverage : ${coverage}" >> ${filename}_log.txt
    echo "completed combined network creation"
else
    echo -e "error occured\nStart Time : ${start}\nError Time : ${end}"
    echo -e "\n\nAUTO2_ERROR_log\nStart Time : ${start}\nERROR Time : ${end}\nTime Elapsed : ${time_elap_min} mins\nBigScale Bin : ${BSbin}\nSAVER Bin : ${SVbin}\nCutoff : ${cutoff}\nWeight : ${weight}\nFinalCoverage : ${coverage}" >> ${filename}_log.txt
fi