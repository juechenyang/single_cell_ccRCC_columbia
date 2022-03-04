while getopts c:n:m: flag
do
    case "${flag}" in
        c) command=${OPTARG};;
        n) nThread=${OPTARG};;
        m) memory=${OPTARG};;
    esac
done
echo "number of thread: $nThread";
echo "memory: $memory"
echo "command is: $command"

bsub -M $memory -n $nThread -W 120:00 -R "span[hosts=1]" \
-J jason -o /users/yanv5j/logs/%J.out -e /users/yanv5j/logs/%J.err \
"$command"