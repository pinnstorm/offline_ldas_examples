#!/usr/bin/env python3
#   Python code to compute the Jacobians, either through the EDA (LUSE_EDA_JACOB=="true") or the finite difference approach 
#   (LUSE_EDA_JACOB=="false")
#
#  Author: original: D. Fairbairn (10/02/2022)
#
#  1. Reads in either (i) variances/covariances for EDA Jacobians or (ii) control and perturbed runs for finite differences
#  2. Constructs the H matrix based on the Jacobians

import os
import sys, getopt
import numpy as np
from netCDF4 import Dataset
import sekf_functions

if __name__ == "__main__":

    nargs = len(sys.argv)

    resol = []
    Time00 = []
    Time06 = []

    scriptname = sys.argv[0].split("/")[-1]
    helptext = (
        scriptname
        + """
     -t <Time at start of window>    REQUIRED either 00 or 12 depending on assimilation window                                                                                                                                                  
     -e <Time after 6 hours>    REQUIRED either 06 or 18 depending on assimilation window                                                                                                                                                  
    """
    )

    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "t:e:r:p:j:i:f:",
            ["oTime00=", "oTime06=", "oRUNDIR=", "oPASTCYCLE=", "oLUSE_EDA_JACOB=","oInitime=", "ofreq="],
        )

    except getopt.GetoptError:
        print(helptext)
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-h", "-help", "--help"):
            print(helptext)
            sys.exit()
        elif opt in ("-t", "--oTime00"):
            Time00 = str(arg)
        elif opt in ("-e", "--oTime06"):
            Time06 = str(arg)
        elif opt in ("-r", "--oRUNDIR"):
            RUNDIR = str(arg)
        elif opt in ("-p", "--oPASTCYCLE"):
            PASTCYCLE = str(arg)
        elif opt in ("-j", "--oLUSE_EDA_JACOB"):
            LUSE_EDA_JACOB = str(arg)
        elif opt in ("-i", "--oInitime"):
            initime = str(arg)
        elif opt in ("-f", "--ofreq"):
            freq = float(arg)

    Times = [Time00, Time06]

    #1. File reading: If LUSE_EDA_JACOB=True, then read in variance and covariance data from netCDF files,
    #   else read in perturbed and unperturbed runs for finite differences
    #    --------------------------------------------------------
    SLV_assim = True

    Control_vars = ["swvl1", "swvl2", "swvl3"]
    Obs = ["ssm", "T2m", "RH2m"]
    Model = dict()
    Data = dict()

    nsteps=int(12.0/freq)
    s_freq=int(1.0/freq)

    if (LUSE_EDA_JACOB=="true"):
        variance = dict()
        covariance = dict()
    else:
        directories=dict()
        perturbation=dict()
        forecast=dict()
        finite_difference = dict()
        ini_pert=dict()

        directories["controldir_00H"]=PASTCYCLE+"/"
        directories["pertdir_00H"]=PASTCYCLE+"/"

        directories["controldir_06H"]=RUNDIR+"/../"+"control/"
        directories["pertdir_06H"]=RUNDIR+"/../"

        ini_pert["pertl1"]=0.01*70.0
        ini_pert["pertl2"]=0.01*210.0
        ini_pert["pertl3"]=0.01*720.0

    H = dict()

    for tt in Times:

        for xx, var in enumerate(Control_vars):
            print("LUSE_EDA_JACOB", LUSE_EDA_JACOB)
            if (LUSE_EDA_JACOB=="true"):
                variance[var + "_" + tt + "_f"] = Dataset(
                    "var_" + tt + "_" + var + ".nc", "r"
                )
                variance[var + "_" + tt] = variance[var + "_" + tt + "_f"].variables[
                    "var_" + tt + "_" + var
                ][
                    0
                ]  # Replace WG1 with swvl1 etc...
                variance[var + "_" + tt + "_f"].close()

                for oo in Obs:
                    covariance[oo + "_" + var + "_" + tt + "_f"] = Dataset(
                        "covar_" + oo + "_" + tt + "_" + var + ".nc", "r"
                    )
                    covariance[oo + "_" + var + "_" + tt] = covariance[
                        oo + "_" + var + "_" + tt + "_f"
                    ].variables["covar_" + oo + "_" + tt + "_" + var][
                        0
                    ]  # Replace WG1 with swvl1 etc...
                    covariance[oo + "_" + var + "_" + tt + "_f"].close()

            else:
                pert_num=xx+1
                for oo in Obs:

                    if (tt==Times[1]):
                        perturbation[oo + "_" + var + "_f"] = Dataset(directories["pertdir_06H"]+oo+"_pertl" + str(pert_num) + ".nc", "r")
                        forecast[oo + "_" + var + "_f"] = Dataset(directories["controldir_06H"]+oo+"_forecast.nc", "r")
                        
                        if oo=="ssm":

                            print(perturbation[oo + "_" + var + "_f"].variables["SoilMoist"][int(nsteps/2),0,:])
                            print(forecast[oo + "_" + var + "_f"].variables["SoilMoist"][int(nsteps/2),0,:])

                            finite_difference[oo + "_" + var + "_" + tt] =  (perturbation[oo + "_" + var + "_f"].variables["SoilMoist"][int(nsteps/2),0,:] - forecast[oo + "_" + var + "_f"].variables["SoilMoist"][int(nsteps/2),0,:])/ini_pert["pertl"+str(pert_num)]
                        else:
                            finite_difference[oo + "_" + var + "_" + tt] =  (perturbation[oo + "_" + var + "_f"].variables[oo][int(nsteps/2)-1,:] - forecast[oo + "_" + var + "_f"].variables[oo][int(nsteps/2)-1,:])/0.01

                        perturbation[oo + "_" + var + "_f"].close()
                        forecast[oo + "_" + var + "_f"].close()

                    elif (tt==Times[0]) and (initime=="false"):

                        perturbation[oo + "_" + var + "_f"] = Dataset(directories["pertdir_00H"]+oo+"_pertl" + str(pert_num) + ".nc", "r")
                        forecast[oo + "_" + var + "_f"] = Dataset(directories["controldir_00H"]+oo+"_forecast.nc", "r")
                        
                        if oo=="ssm":
                            finite_difference[oo + "_" + var + "_" + tt] =  (perturbation[oo + "_" + var + "_f"].variables["SoilMoist"][-1,0,:] - forecast[oo + "_" + var + "_f"].variables["SoilMoist"][-1,0,:])/ini_pert["pertl"+str(pert_num)]
                            
                        else:
                            finite_difference[oo + "_" + var + "_" + tt] =  (perturbation[oo + "_" + var + "_f"].variables[oo][-1,:] - forecast[oo + "_" + var + "_f"].variables[oo][-1,:])/0.01

                        perturbation[oo + "_" + var + "_f"].close()
                        forecast[oo + "_" + var + "_f"].close()

        if (LUSE_EDA_JACOB=="true"):
            EDAH_TAPER = 0.6

            H[tt] = sekf_functions.EDA_Jacobian_calc(
                covariance["T2m_swvl1_" + tt][:],
                covariance["T2m_swvl2_" + tt][:],
                covariance["T2m_swvl3_" + tt][:],
                covariance["RH2m_swvl1_" + tt][:],
                covariance["RH2m_swvl2_" + tt][:],
                covariance["RH2m_swvl3_" + tt][:],
                covariance["ssm_swvl1_" + tt][:],
                covariance["ssm_swvl2_" + tt][:],
                covariance["ssm_swvl3_" + tt][:],
                variance["swvl1_" + tt][:],
                variance["swvl2_" + tt][:],
                variance["swvl3_" + tt][:],
                EDAH_TAPER,
            )
            np.save("variance_swvl3_"+tt, np.array(variance["swvl3_" + tt][:]))

        else:
            if ((tt=="06") or (tt=="18")) or (initime=="false"):
                print(tt)    
                H[tt] = np.array([[finite_difference["T2m_swvl1_" + tt],\
                                   finite_difference["T2m_swvl2_" + tt],\
                                   finite_difference["T2m_swvl3_" + tt]],\
                                  [finite_difference["RH2m_swvl1_" + tt],\
                                   finite_difference["RH2m_swvl2_" + tt],\
                                   finite_difference["RH2m_swvl3_" + tt]],\
                                  [finite_difference["ssm_swvl1_" + tt],\
                                   finite_difference["ssm_swvl2_" + tt],\
                                   finite_difference["ssm_swvl3_" + tt]]],dtype="float64")

    if (initime=="true"):
      H[Times[0]]=H[Times[1]]

#    print("Times 1:")

#    print("T2m_swvl1_", np.mean(abs(finite_difference["T2m_swvl1_" + Times[0]])), len(np.where(finite_difference["T2m_swvl1_" + Times[0]][:]>0)[0]))
#    print("T2m_swvl2_", np.mean(abs(finite_difference["T2m_swvl2_" + Times[0]])), len(np.where(finite_difference["T2m_swvl2_" + Times[0]][:]>0)[0]))
#    print("T2m_swvl3_", np.mean(abs(finite_difference["T2m_swvl3_" + Times[0]])), len(np.where(finite_difference["T2m_swvl3_" + Times[0]][:]>0)[0]))
#    print("RH2m_swvl1_", np.mean(abs(finite_difference["RH2m_swvl1_" + Times[0]])), len(np.where(finite_difference["RH2m_swvl1_" + Times[0]][:]>0)[0]))
#    print("RH2m_swvl2_", np.mean(abs(finite_difference["RH2m_swvl2_" + Times[0]])), len(np.where(finite_difference["RH2m_swvl2_" + Times[0]][:]>0)[0]))
#    print("RH2m_swvl3_", np.mean(abs(finite_difference["RH2m_swvl3_" + Times[0]])), len(np.where(finite_difference["RH2m_swvl3_" + Times[0]][:]>0)[0]))
#    print("ssm_swvl1_", np.mean(abs(finite_difference["ssm_swvl1_" + Times[0]])), len(np.where(finite_difference["ssm_swvl1_" + Times[0]][:]>0)[0]))
#    print("ssm_swvl2_", np.mean(abs(finite_difference["ssm_swvl2_" + Times[0]])), len(np.where(finite_difference["ssm_swvl2_" + Times[0]][:]>0)[0]))
#    print("ssm_swvl3_", np.mean(abs(finite_difference["ssm_swvl3_" + Times[0]])), len(np.where(finite_difference["ssm_swvl3_" + Times[0]][:]>0)[0]))

#    print("Times 1:")

#    print("T2m_swvl1_", np.mean(abs(finite_difference["T2m_swvl1_" + Times[1]])), len(np.where(finite_difference["T2m_swvl1_" + Times[1]][:]>0)[0]))
#    print("T2m_swvl2_", np.mean(abs(finite_difference["T2m_swvl2_" + Times[1]])), len(np.where(finite_difference["T2m_swvl2_" + Times[1]][:]>0)[0]))
#    print("T2m_swvl3_", np.mean(abs(finite_difference["T2m_swvl3_" + Times[1]])), len(np.where(finite_difference["T2m_swvl3_" + Times[1]][:]>0)[0]))
#    print("RH2m_swvl1_", np.mean(abs(finite_difference["RH2m_swvl1_" + Times[1]])), len(np.where(finite_difference["RH2m_swvl1_" + Times[1]][:]>0)[0]))
#    print("RH2m_swvl2_", np.mean(abs(finite_difference["RH2m_swvl2_" + Times[1]])), len(np.where(finite_difference["RH2m_swvl2_" + Times[1]][:]>0)[0]))
#    print("RH2m_swvl3_", np.mean(abs(finite_difference["RH2m_swvl3_" + Times[1]])), len(np.where(finite_difference["RH2m_swvl3_" + Times[1]][:]>0)[0]))
#    print("ssm_swvl1_", np.mean(abs(finite_difference["ssm_swvl1_" + Times[1]])), len(np.where(finite_difference["ssm_swvl1_" + Times[1]][:]>0)[0]))
#    print("ssm_swvl2_", np.mean(abs(finite_difference["ssm_swvl2_" + Times[1]])), len(np.where(finite_difference["ssm_swvl2_" + Times[1]][:]>0)[0]))
#    print("ssm_swvl3_", np.mean(abs(finite_difference["ssm_swvl3_" + Times[1]])), len(np.where(finite_difference["ssm_swvl3_" + Times[1]][:]>0)[0]))

    Hmat = np.concatenate((H[Times[0]], H[Times[1]]), axis=0)
    np.save("HMAT", Hmat)


