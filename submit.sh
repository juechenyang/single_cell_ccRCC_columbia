#BSUB -W 120:00
#BSUB -M 96000
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -J jason_infercnv
#BSUB -o /users/yanv5j/logs/%J.out
#BSUB -e /users/yanv5j/logs/%J.err

Rscript get_filtered_integrated_cancer_cells.R
