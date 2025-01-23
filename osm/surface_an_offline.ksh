#!/bin/ksh
#   Job to run the simplified extended Kalman filter over 12-hour windows
#   1. Link all netcdf files of model and observations to run directory
#   2. For each 12-hour window, the SEKF is called from the python script, with
#      the increment at the end of the assimilation window
#   3. netCDF output files are produced and then analyses saved in the lag family
#  Author: original: D. Fairbairn (10/02/2021)
#

set -eux
LOAD_MODULE netcdf4
LOAD_MODULE nco
LOAD_MODULE python3

##=================================================
## 1. DEFAULT

INIDATE=$(echo $INIBEGINDATE | cut -c1-8)

PastCycle=$( dateincr -d ${DATELABEL} -1 ) #Last cycle date
NDAYS=0.5

todate=$(date -d ${DATELABEL} '+%%s')
ascat_start=$( date -d 19920101 '+%%s' )

if [ $todate -lt $ascat_start ]; then
  LUSE_ASCAT="false"
fi  

ims_start=$( date -d 20100901 '+%%s' )
if [ $todate -lt $ims_start ]; then
  scf_threshold=0.3
  qc_threshold=0.8
else
  scf_threshold=0.5
  qc_threshold=1.0
fi

echo $RUNDIR
##================================================
## 1. Link model and observations to run directory

#Setup results and run directories
if [[ $FAMILY == */lw00/* ]] ; then
  RESDIR=${DATA}/osm/outdata/${DATELABEL}/lw00/$SUBFSFAMILY 
  RUNDIR=${DATA}/osm/outdata/${DATELABEL}/lw00/an_dir
  if [ -d $RUNDIR ] ; then rm -rf $RUNDIR ; fi
  mkdir -p $RUNDIR
elif [[ $FAMILY == */lw12/* ]] ; then
  RESDIR=${DATA}/osm/outdata/${DATELABEL}/lw12/$SUBFSFAMILY
  RUNDIR=${DATA}/osm/outdata/${DATELABEL}/lw12/an_dir
  if [ -d $RUNDIR ] ; then rm -rf $RUNDIR ; fi
  mkdir -p $RUNDIR
fi

cd $RUNDIR

# Link forcing files and observations to RUNDIR location
ln -s -f ${RESDIR}/*.nc .
ln -s -f ${RESDIR}/*.grib .
ln -s -f ${FIXDIR}/surfclim .

# Link AWE to run directory

ln -s -f $RUNDIR/../control/AWC.nc .
LSLV_ASSIM=true
if [[ $LSLV_ASSIM = true ]]; then
  ln -s -f $RUNDIR/../control/o_d2m.nc .
fi

echo $PWD

## functions
gridT(){
resol=$1
gtype=$2  
case ${gtype} in
  "l_2"|"_3") GT=N;;
  "_4") GT=O;;
  *) print "GTYPE not available!"; exit -1;;
esac
NLAT=$( gaussgr -r $1 -g $2 )
echo ${GT}${NLAT}
}

MRESOL=$(gridT $RESOL $GTYPE )

## SEKF ################################
LESNML=${LESNML:-false}
NCSNEC=${NCSNEC:-1}

ncks -v SoilMoist restartout.nc SoilMoist_control.nc  # Make copy of soil moisture for SEKF control
ncks -v SWE       restartout.nc SWE_control.nc        # Make copy of SWE for SEKF control
ncks -v snowdens  restartout.nc snowdens_control.nc   # Make copy of snow density for SEKF control
if [[ $LESNML = true ]]; then
  ncks -v SWEML      restartout.nc SWEML_control.nc      # Make copy of multi-layered SWE for SEKF control
  ncks -v snowdensML restartout.nc snowdensML_control.nc # Make copy of multi-layered snow density for SEKF control
fi

if [[ $DATELABEL = $INIDATE ]]; then # first month of run
  initime=true 
else
  initime=false
fi

INI_BIAS="false"
if [[ $LUSE_ASCAT = true ]] && [[ $LASCAT_ADAPTIVE_BC = true ]]; then 

  if [[ $DATELABEL = $INIDATE ]] && [[ $FAMILY == */lw00/* ]]; then 
    INI_BIAS="true"
    if [[ $LASCAT_ADAPTIVE_BC == "true" ]] && [[ $ADAPTIVE_BC_PATH == coldstart ]] || [[ $ADAPTIVE_BC_PATH == standard ]]; then 
      #do -cat -apply,-selname,SoilMoist,time [ restartout.nc ] sekf_ascat_bias_in.nc
      ncks -O -d nlevs,1 -v SoilMoist,time restartout.nc sekf_ascat_bias_in.nc
      ncrename -h -O -v SoilMoist,sekf_ascat_bias_in sekf_ascat_bias_in.nc
      INI_BIAS=$ADAPTIVE_BC_PATH
    fi
  elif [[ $FAMILY == */lw00/* ]] && [[ -f "${RUNDIR}/../../../${PastCycle}/lw12/sekf_ascat_bias_in.nc" ]]; then
    cp ${RUNDIR}/../../../${PastCycle}/lw12/sekf_ascat_bias_in.nc .
  elif [[ $FAMILY == */lw12/* ]] && [[ -f "${RUNDIR}/../../lw00/sekf_ascat_bias_in.nc" ]]; then
    cp ${RUNDIR}/../../lw00/sekf_ascat_bias_in.nc . #retrieve previous T2m
  fi
fi

ncks -O -v SWE o_gg.nc swe.nc

echo "here........ "

cp T2m_forecast.nc SEKF_save.nc
cp T2m_forecast.nc SEKF_snow_save.nc

#   Analysis script
if [[ $FAMILY == */lw00/* ]]; then
  if [[ $LSLV_ASSIM = true ]]; then
    if [[ -f "${RUNDIR}/../../../${PastCycle}/lw12/T2m_forecast.nc" ]]; then
      ln -s ${RUNDIR}/../../../${PastCycle}/lw12/T2m_forecast.nc T2m_prev.nc #retrieve previous T2m
      ln -s ${RUNDIR}/../../../${PastCycle}/lw12/RH2m_forecast.nc RH2m_prev.nc #retrieve previous RH2m
    else
      cp T2m_forecast.nc T2m_prev.nc
      cp RH2m_forecast.nc RH2m_prev.nc
    fi
  fi

  PASTDIR=${RUNDIR}/../../../${PastCycle}/lw12

  if [[ $LSEKF_SM = true ]]; then
    $TMPDIR/compute_SEKF_Jacobians.py -t 00 -e 06 -r $RUNDIR -p ${PASTDIR} -j $LUSE_EDA_JACOB -i $initime -f ${FREQ_PP_OSM}    #Run python script to
 
    $TMPDIR/sekf_ascat_slv_assim.py -t 00 -e 06 -s $LSAVE_SEKF -a $LUSE_ASCAT -j $LUSE_EDA_JACOB -f ${FREQ_PP_OSM} -b ${LASCAT_ADAPTIVE_BC} -i $INI_BIAS #Run python script to assimilate ASCAT and SLV obs

    ncks -h -a -A SoilMoist_control.nc restartout.nc 

    if [[ $LUSE_EDA_JACOB = false ]]; then
      cp -f ssm_forecast.nc ${DATA}/osm/outdata/%YMD%/lw00/ssm_forecast.nc 
    fi

    if [[ $LASCAT_ADAPTIVE_BC = true ]]; then
      cp -f sekf_ascat_bias_in.nc ${DATA}/osm/outdata/%YMD%/lw00/sekf_ascat_bias_in.nc
    fi
  fi

  if [[ $LSEKF_SD = true ]]; then
    $TMPDIR/sekf_snow_assim.py -t 00 -e 06 -s $LSAVE_SEKF -m $LESNML -n $NCSNEC -f ${FREQ_PP_OSM} -b ${scf_threshold} -q ${qc_threshold}

    ncks -h -a -A SWE_control.nc restartout.nc
    ncks -h -a -A snowdens_control.nc restartout.nc
    if [[ $LESNML = true ]]; then
      ncks -h -a -A SWEML_control.nc restartout.nc
      ncks -h -a -A snowdensML_control.nc restartout.nc
    fi
  fi

  cp -f restartout.nc ${DATA}/osm/outdata/%YMD%/lw00/restartout.nc
  cp -f T2m_forecast.nc ${DATA}/osm/outdata/%YMD%/lw00/T2m_forecast.nc
  cp -f RH2m_forecast.nc ${DATA}/osm/outdata/%YMD%/lw00/RH2m_forecast.nc

else
  if [[ $LSLV_ASSIM = true ]]; then
    if [[ -f "${RUNDIR}/../../lw00/T2m_forecast.nc" ]]; then
      ln -s ${RUNDIR}/../../lw00/T2m_forecast.nc T2m_prev.nc #retrieve previous T2m
      ln -s ${RUNDIR}/../../lw00/RH2m_forecast.nc RH2m_prev.nc #retrieve previous RH2m
    else
      cp T2m_forecast.nc T2m_prev.nc
      cp RH2m_forecast.nc RH2m_prev.nc
    fi
  fi

  PASTDIR=${RUNDIR}/../../lw00

  if [[ $LSEKF_SM = true ]]; then
    $TMPDIR/compute_SEKF_Jacobians.py -t 12 -e 18 -r $RUNDIR -p ${PASTDIR} -j $LUSE_EDA_JACOB -i $initime -f ${FREQ_PP_OSM}     #Run python script to
 
    $TMPDIR/sekf_ascat_slv_assim.py -t 12 -e 18 -s $LSAVE_SEKF -a $LUSE_ASCAT -j $LUSE_EDA_JACOB -f ${FREQ_PP_OSM} -b ${LASCAT_ADAPTIVE_BC} -i $INI_BIAS  #Run python script to assimilate ASCAT and SLV obs

    ncks -h -a -A SoilMoist_control.nc restartout.nc 

    if [[ $LUSE_EDA_JACOB = false ]]; then
      cp -f ssm_forecast.nc ${DATA}/osm/outdata/%YMD%/lw12/ssm_forecast.nc 
    fi

    if [[ $LASCAT_ADAPTIVE_BC = true ]]; then
      cp -f sekf_ascat_bias_in.nc ${DATA}/osm/outdata/%YMD%/lw12/sekf_ascat_bias_in.nc 
    fi
  fi

  if [[ $LSEKF_SD = true ]]; then
    $TMPDIR/sekf_snow_assim.py -t 12 -e 18 -s $LSAVE_SEKF -m $LESNML -n $NCSNEC -f ${FREQ_PP_OSM} -b ${scf_threshold} -q ${qc_threshold}

    ncks -h -a -A SWE_control.nc restartout.nc
    ncks -h -a -A snowdens_control.nc restartout.nc
    if [[ $LESNML = true ]]; then
      ncks -h -a -A SWEML_control.nc restartout.nc
      ncks -h -a -A snowdensML_control.nc restartout.nc
    fi
  fi

  cp -f restartout.nc ${DATA}/osm/outdata/%YMD%/lw12/restartout.nc
  cp -f T2m_forecast.nc ${DATA}/osm/outdata/%YMD%/lw12/T2m_forecast.nc
  cp -f RH2m_forecast.nc ${DATA}/osm/outdata/%YMD%/lw12/RH2m_forecast.nc
fi


#clean obs.
rm -f Ascat*.nc ascat*.nc RH*.nc U*.nc V*.nc H_Ascat*.nc SCF*.nc #T2*.nc

exit 0 
