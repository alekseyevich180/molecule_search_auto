#!/bin/bash
#PJM -L rscgrp=a-pj24001724
#PJM -L node=1
#PJM --mpi proc=120
#PJM -L elapse=6:00:00
#PJM -j

module load intel
module load impi

source /home/pj24001724/ku40000345/venv/TPE_search/.venv/bin/activate
python optuna_ad_search.py
deactivate