cat <(head -n1 *predictions.tsv | grep -v "^==>" | head -n1) <(tail -n+2 *predictions.tsv | egrep -v '^==>|^$') > all_predictions.tsv

rm *snps_predictions.tsv
