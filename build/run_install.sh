#!/bin/bash

pwd=`pwd`

flibinstall=${pwd}

cd ${pwd}/../gsl-2.5;
./configure --prefix=${flibinstall};
make;
make install;
cd ${pwd};

cd ${pwd}/../pyWaveformGenerator;
./mkconf.sh;
./configure;
C_INCLUDE_PATH=${pwd}/include:${C_INCLUDE_PATH} LD_LIBRARY_PATH=${pwd}/lib:${pwd}/lib64:${LD_LIBRARY_PATH} make;
cd ${pwd}
