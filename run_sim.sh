#!/bin/bash

export PATH=$PATH:"C:\Program Files\R\R-4.2.1\bin"
nohup Rscript --vanilla R/server_simulations_workflow.R > simulations.out &
