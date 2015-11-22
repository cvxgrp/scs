#!/bin/bash

# clean up any remaining files from previous runs
./cleanup
cd .. && make purge && cd r

# run autoconf to generate config
# usually not necessary unless configure file missing or out of date
#autoconf configure.ac > configure
#chmod +x configure

# install SCS
R CMD INSTALL .

# clean up again
./cleanup

# run basic example
R --no-save < demo/randomLp.R

