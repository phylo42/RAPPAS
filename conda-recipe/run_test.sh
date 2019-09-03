#basic execution test
java -jar $(which RAPPAS.jar) -h
t0=$(echo $?)

#verify that dabatase can be computed
rappas -w test -s nucl -p b --omega 2.0 -b $(which phyml) -t test/150.newick -r test/150.fasta --use_unrooted
test -f test/DB_session_k8_a2.0_f1_t-2.40824.union
t1=$(echo $?)

#verify that jplace can be computed 
rappas -w test -p p -d test/DB_session_k8_a2.0_f1_t-2.40824.union -q test/q.fasta
test -f test/placements_q.fasta_union.jplace
t2=$(echo $?)

#remove tests outputs
rm -f test/align.reduced
rm -f test/DB_session_k8_a2.0_f1_t-2.40824.union
rm -f test/placements_q.fasta_union.jplace
rm -rf test/logs
rm -rf test/AR
rm -rf test/extended_trees

#both test successful ?
if [[ $t0 == 0 && $t1 == 0 && $t2 == 0 ]]
  then
    echo "All tests successful !" ;
    exit 0
  else
    echo "Test failed"
    if [[ $t0 != 0 ]]; then echo "launch = failed"; fi
    if [[ $t1 != 0 ]]; then echo "dbbluild = failed"; fi
    if [[ $t2 != 0 ]]; then echo "placement = failed"; fi
    exit 1
fi
