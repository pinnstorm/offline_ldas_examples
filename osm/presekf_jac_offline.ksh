#!/bin/ksh
#
#  Job to retrieve SEKF EDA Jacobians offline LDAS
#
#  Author: original: D. Fairbairn and K. Ochi and E. Pinnington, June 2023
#
#  - EDA retrieved and converted to netCDF4
# 
# set -e
#
# Setup for extraction
#-----------------
#  Set up default values 
# Script to perform monthly run
set -exu


##=================================================
## --- 1 --- setup default variables
##===================================

## default
EDATELABEL=${DATELABEL}

ZDTFORC=$(( $FCFREQ * 3600 ))   # Forcing frequency in seconds 

## functions
gridT(){
resol=$1
gtype=$2  
case ${gtype} in
  "l_2"|"_3"|"_2") GT=N;;
  "_full") GT=F;;
  "_4") GT=O;;
  *) print "GTYPE not available!"; exit -1;;
esac
NLAT=$( gaussgr -r $1 -g $2 )
echo ${GT}${NLAT}
}
gridTF(){

NLAT=$(grib_ls -x -p N -w count=1 $1 | head -n3 | tail -n1 | awk '{print $1}' )
isO=$(grib_ls -x -p isOctahedral -w count=1 $1 | head -n3 | tail -n1 | awk '{print $1}')
GT=N
if [[ $isO == 1 ]] ; then GT=O ; fi
echo ${GT}${NLAT}
}

## check and create DATADIR
if [ ! -d $DATADIR ] ; then mkdir -p $DATADIR ; fi
if [ ! -d $DATADIR/jac ] ; then mkdir -p $DATADIR/jac ; fi
cd $DATADIR/jac

## Check if interpolation is required : 

MRESOL=$(gridT $RESOL $GTYPE )

##====================================================================================
## --- 2 --- extract ASCAT, IMS and SLVs and convert assimilation grib files to netCDF 
##====================================================================================

extract_sekf_jacobians.ksh -e $EDATELABEL -i $FCEXPVER -d $DATADIR/jac -f $FCFREQ -c $FCCLASS -s $FCSTREAM -r ${MRESOL}

##========================================
## netcdf file template
function create_cdl {

var=$1
long_name=$2
units=$3
ddate=$4

yyyy=$(echo $ddate | cut -c1-4)
mm=$(echo $ddate | cut -c5-6)
dd=$(echo $ddate | cut -c7-8)
if [[ $dd = '00' ]]; then
  # for LERA runs dd == 1
  dd='01'
fi

cat > ${var}.cdl << EOF

netcdf ${var} {
dimensions:
        time = UNLIMITED ; 
        x  = $NPOINTS;
variables:
        int    x(x);
                x:long_name ="grid points ";
                x:units="-";
        float  lat(x) ;
                lat:long_name = "latitude" ;
                lat:units = "degrees_north" ;
        float  lon(x) ;
                lon:long_name = "longidude" ;
                lon:units = "degrees_east" ;
        double time(time) ;
                time:long_name = "time " ;
                time:units = "seconds since $yyyy-$mm-$dd 00:00:00" ;
        float ${var}(time,x) ;
                ${var}:long_name = "${long_name}" ;
                ${var}:units = "${units}" ;

// global attributes:
                :SOURCE = "ECMWF" ;
                :GRID_POINTS = "gaussian grid " ;
                :CONVERTED = " from grib files with grib_api"; 
}
EOF


rm -f ${var}.nc
ncgen -b ${var}.cdl
rm ${var}.cdl
}

CVARS=""

  ## create empty nc forcing files 
  NPOINTS=$(ncdump -h $FIXDIR/surfclim | grep "x = " | awk '{print $3}' )


##=================================================
## --- 3 --- Convert grib files to netCDF 
##=================================================

  #SEKF SSM Jacobians
  ##===================================

  for HH in 00 06 12 18
  do

    CVARS_SEKF_SSM_jac="covar_ssm_${HH}_swvl1 covar_ssm_${HH}_swvl2 covar_ssm_${HH}_swvl3  
    var_${HH}_swvl1 var_${HH}_swvl2 var_${HH}_swvl3"
    
    for ff in $CVARS_SEKF_SSM_jac
    do
      if [[ -r ${ff}.grb ]]; then
        CVARS="$CVARS ${ff}"
      else
        exit 1
      fi
    done


    create_cdl 'covar_ssm_'$HH'_swvl1' 'covar (swvl1,swvl1)' '-' $DATELABEL 
    create_cdl 'covar_ssm_'$HH'_swvl2' 'covar (swvl1,swvl2)' '-' $DATELABEL 
    create_cdl 'covar_ssm_'$HH'_swvl3' 'covar (swvl1,swvl3)' '-' $DATELABEL 
    create_cdl 'var_'$HH'_swvl1' 'var (swvl1)' '-' $DATELABEL 
    create_cdl 'var_'$HH'_swvl2' 'var (swvl2)' '-' $DATELABEL 
    create_cdl 'var_'$HH'_swvl3' 'var (swvl3)' '-' $DATELABEL 
  done



  #SEKF T2m/RH2m Jacobians
  ##===================================

  for HH in 00 06 12 18
  do
    CVARS_SEKF_SLV="covar_RH2m_${HH}_swvl1 covar_RH2m_${HH}_swvl2 covar_RH2m_${HH}_swvl3 covar_T2m_${HH}_swvl1 covar_T2m_${HH}_swvl2 covar_T2m_${HH}_swvl3"
    for ff in $CVARS_SEKF_SLV
    do
      if [[ -r ${ff}.grb ]]; then
        CVARS="$CVARS ${ff}"
      else
        exit 1
      fi
    done

    create_cdl 'covar_T2m_'${HH}'_swvl1' 'covar (T2m,swvl1)' '-' $DATELABEL 
    create_cdl 'covar_T2m_'${HH}'_swvl2' 'covar (T2m,swvl2)' '-' $DATELABEL 
    create_cdl 'covar_T2m_'${HH}'_swvl3' 'covar (T2m,swvl3)' '-' $DATELABEL 
    create_cdl 'covar_RH2m_'${HH}'_swvl1' 'covar (RH2m,swvl1)' '-' $DATELABEL 
    create_cdl 'covar_RH2m_'${HH}'_swvl2' 'covar (RH2m,swvl2)' '-' $DATELABEL 
    create_cdl 'covar_RH2m_'${HH}'_swvl3' 'covar (RH2m,swvl3)' '-' $DATELABEL 

  done



for var in $CVARS 
do
  echo $var 
  ## 3.2 create input namelist 
cat > input.nam <<EOF
&INPUT
  GRIB_FILE='${var}.grb'
  NETCDF_FILE='${var}.nc'
  VAR_NAME='${var}'
  INFO_FILE='$FIXDIR/surfclim'
  ZDTFORC=${ZDTFORC}
/
EOF

  ## 3.3 run conv_forcing to transform grib files to netCDF4
  ${BINS}/conv_forcing 
  if [ $? -ne 0 ]; then
    echo "some probelm in conv_forcing"
    exit -1
  fi 
  rm -f input.nam
  rm ${var}.grb

done

mv $DATADIR/jac/*.nc $DATADIR

exit 0
