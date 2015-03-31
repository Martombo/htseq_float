for k in {0..29000}; do awk -v k=$k 'BEGIN{OFS="\t"}{$4=$4+k; print}' prova_wholeLen.gtf >prova_tmp.gtf; htseq-count -a 0 -s no prova_uniqueMapped.sam prova_tmp.gtf >>output; done
