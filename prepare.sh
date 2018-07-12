#!/bin/sh
source /public/software/intel/ifort2015/composer_xe_2015.0.090/bin/compilervars.sh intel64

source /public/software/intel/ifort2015/impi_latest/intel64/bin/mpivars.sh

make
make run
