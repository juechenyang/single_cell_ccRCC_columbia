while getopts f:n:m: flag
do
    case "${flag}" in
        f) rfile=${OPTARG};;
        n) nThread=${OPTARG};;
        m) memory=${OPTARG};;
    esac
done
echo "number of thread: $nThread";
echo "memory: $memory"
echo "command is: $rfile"

bsub -M $memory -n $nThread -W 120:00 -R "span[hosts=1]" \
-J jason -o /users/yanv5j/logs/%J.out -e /users/yanv5j/logs/%J.err \
"Rscript $rfile"