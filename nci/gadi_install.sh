#!/bin/bash

module purge

module load R/4.2.1
module load gdal/3.5.0
module load proj/8.1.1
module load geos/3.8.0 udunits/2.2.26 unixODBC/2.3.7
module load gsl/2.6  # for the Rfast package
module load intel-compiler
module load openmpi/4.1.0

exit_on_failure() {
    if [[ $1 != 0 ]]
    then
        echo $2
        exit 1
    fi
}

if [[ -z ${R_LIBS+x} ]]
then
    echo "Installing R packages to `R --slave -e '.libPaths()[1]'`. If this is not correct, change the R_LIBS environment variable."
else
    echo "Installing R packages to '$R_LIBS'. If this is not correct, change the R_LIBS environment variable."
fi

read -p "Continue? (y|Y)" -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]
then
  R --slave -e "install.packages(c('rgdal', 'R.utils', 'raster', 'rasterVis', 'entropy', 'rstudioapi', 'moments', 'rgeos', 'BalancedSampling', 'ggplot2', 'reshape2', 'plyr', 'scales', 'snow', 'Rfast'), repos='https://mirror.aarnet.edu.au/pub/CRAN/', lib='$R_LIBS', depdencies=TRUE)"
  R --slave -e "install.packages('Rmpi', configure.args=c('--with-Rmpi-type=OPENMPI --with-Rmpi-include=/apps/openmpi/4.1.0/include/ --with-Rmpi-libpath=/apps/openmpi/4.1.0/lib/'), repos='https://mirror.aarnet.edu.au/pub/CRAN/', lib='$R_LIBS', depdencies=TRUE)"
  exit_on_failure $? "Error: could not install R packages."
  CURDIR=`pwd`
  cd $HOME
  # 'clhs' is not available from CRAN at time of writing (14-04-2020), so need to download from archive and install.
  curl -O https://cran.r-project.org/src/contrib/Archive/clhs/clhs_0.7-2.tar.gz
  tar -xzf clhs_0.7-2.tar.gz
  R --slave -e "install.packages('$HOME/clhs', repos=NULL, type='source', lib='$R_LIBS')"
    curl -O https://cran.r-project.org/src/contrib/Archive/SDMTools/SDMTools_1.1-221.2.tar.gz
  tar -xzvf SDMTools_1.1-221.2.tar.gz
  R --slave -e "install.packages('$HOME/SDMTools', repos=NULL, type='source', lib='$R_LIBS')"
  RET=$?
  rm clhs_0.7-2.tar.gz SDMTools_1.1-221.2.tar.gz
  rm -r clhs SDMTools
  cd $CURDIR
  exit_on_failure $RET "Error: could not install 'clhs'."

  echo "'geo-sampling' install complete. Ensure you set 'R_LIBS' environment variable to the location where you installed these R packages when submitting a job."
fi
