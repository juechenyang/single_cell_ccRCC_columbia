#BSUB -W 120:00
#BSUB -M 80000
#BSUB -n 8
#BSUB -R "span[hosts=1]"
#BSUB -J jason_infercnv
#BSUB -o /users/yanv5j/logs/%J.out
#BSUB -e /users/yanv5j/logs/%J.err

Rscript run_infercnv.R
