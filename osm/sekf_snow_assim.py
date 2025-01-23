#!/usr/bin/env python3
#   Python code to run the simplified extended Kalman filter
#
#  Author: original: D. Fairbairn (31/01/2020)
#
#  1. Reads in model output and SCF obs
#    (The observations are pre-extracted in netCDF4 files over the land regions)
#  2. The SEKF assimilates SCF pseudo-obs at T+0 and the analysis increment is added at the
#    end of the assimilation window (T+12)

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
            "t:e:s:m:n:f:b:q:",
            ["oTime00=", "oTime06=", "oLSAVE_SEKF=", "oLESNML=", "oNCSNEC=", "ofreq=", "oscf_threshold=", "oqc_threshold="],
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
        elif opt in ("-s", "--oLSAVE_SEKF"):
            LSAVE_SEKF = str(arg)
        elif opt in ("-m", "--oLESNML"):
            LESNML = str(arg)
        elif opt in ("-n", "--oNCSNEC"):
            NCSNEC = int(arg)
        elif opt in ("-f", "--ofreq"):
            freq = float(arg)
        elif opt in ("-b", "--oscf_threshold"):
            scf_threshold = float(arg)
        elif opt in ("-q", "--oqc_threshold"):
            qc_threshold = float(arg)

    nsteps=int(12.0/freq)

    s_freq=int(1.0/freq)

    Times = [Time00, Time06]

    Model = dict()
    Data = dict()

    snml = ""
    if LESNML == "true":
        snml = "ML"
    else:
        snml = ""

    Model["SWE_control_f"] = Dataset("SWE_control.nc", "r+")
    Model["SWE_control"] = Model["SWE_control_f"].variables["SWE"]
    Model["snowdens_control_f"] = Dataset("snowdens_control.nc", "r+")
    Model["snowdens_control"] = Model["snowdens_control_f"].variables["snowdens"]

    if LESNML == "true":
        Model["SWEML_control_f"] = Dataset("SWEML_control.nc", "r+")
        Model["SWEML_control"] = Model["SWEML_control_f"].variables["SWEML"]
        Model["snowdensML_control_f"] = Dataset("snowdensML_control.nc", "r+")
        Model["snowdensML_control"] = Model["snowdensML_control_f"].variables["snowdensML"]
        Model["SD_control"] = Model["SWE_control"][:] * 0.0
        for k in range(NCSNEC):
            Model["SD_control"][:] += Model["SWEML_control"][k, :] / Model["snowdensML_control"][k, :]
    else:
        Model["SD_control"] = Model["SWE_control"][:] / Model["snowdens_control"][:]

    Model["SWE_background_f"] = Dataset("SWE_forecast.nc", "r")
    Model["SWE_background"] = Model["SWE_background_f"].variables["SWE" + snml][0:nsteps:s_freq]
    Model["SWE_background_f"].close()

    Model["snowdens_background_f"] = Dataset("snowdens_forecast.nc", "r")
    Model["snowdens_background"] = Model["snowdens_background_f"].variables["snowdens" + snml][0:nsteps:s_freq]
    Model["snowdens_background_f"].close()

    if LESNML == "true":
        Model["SD_background"] = Model["SWE_background"][:, 0, :] * 0.0
        for k in range(NCSNEC):
            Model["SD_background"][:, :] += Model["SWE_background"][:, k, :] / Model["snowdens_background"][:, k, :]
    else:
        Model["SD_background"] = Model["SWE_background"][:, :] / Model["snowdens_background"][:, :]

    if LSAVE_SEKF == "true":
        Model["SEKF_save_f"] = Dataset("SEKF_snow_save.nc", "r+")
        Model["SEKF_save"] = Model["SEKF_save_f"].variables["T2m"]
        Model["SEKF_save"][:, :] = 0.0

    for tt in Times:

        Data["T2m" + tt + "_f"] = Dataset("T2m" + tt + ".nc", "r")
        Data["T2m" + tt] = Data["T2m" + tt + "_f"].variables["T2m" + tt][:].flatten()
        Data["T2m" + tt + "_f"].close()
        Data["T2m" + tt] = Data["T2m" + tt][:] - 273.16

        Data["U10m" + tt + "_f"] = Dataset("U10m" + tt + ".nc", "r")
        Data["U10m" + tt] = Data["U10m" + tt + "_f"].variables["U10m" + tt][0]
        Data["U10m" + tt + "_f"].close()

        Data["V10m" + tt + "_f"] = Dataset("V10m" + tt + ".nc", "r")
        Data["V10m" + tt] = Data["V10m" + tt + "_f"].variables["V10m" + tt][0]
        Data["V10m" + tt + "_f"].close()

    ImsChk = os.path.exists("SCF" + Time00 + ".nc")
    if ImsChk:
        Data["SCF" + Time00 + "_f"] = Dataset("SCF" + Time00 + ".nc", "r")
        Data["SCF" + Time00] = Data["SCF" + Time00 + "_f"].variables["SCF" + Time00][0]
        Data["SCF" + Time00 + "_f"].close()
        Data["SCF_clim_f"] = Dataset("SCF_clim.nc", "r")
        Data["SCF_clim"] = Data["SCF_clim_f"].variables["SCF_clim"][0]
        Data["SCF_clim_f"].close()
    else:
        Data["SCF" + Time00] = Data["T2m" + Time00][:] * 0.0
        Data["SCF" + Time00][:] = -1.0
        Data["SCF_clim"] = Data["T2m" + Time00][:] * 0.0
        Data["SCF_clim"][:] = -1.0

    Model["sdfor_f"] = Dataset("surfclim", "r")
    Model["sdfor"] = Model["sdfor_f"].variables["sdfor"][:]
    Model["sdfor_f"].close()

    Model["glacier_mask_f"] = Dataset("surfclim", "r")
    Model["glacier_mask"] = Model["glacier_mask_f"].variables["glm"][:]
    Model["glacier_mask_f"].close()

    increment = Model["SD_control"][:] * 0.0

    if LSAVE_SEKF == "true":
        Model["SEKF_save"][0, :] = Model["SD_control"][:]
        Model["SEKF_save"][1, :] = Data["SCF" + Time00][:]
        Model["SEKF_save"][2, :] = Data["SCF_clim"][:]

    p_screen = np.where(
        (Data["SCF" + Time00][:] >= 0.0)
        & (Model["SD_control"][:] >= 0.0)
        & (Model["glacier_mask"][:] <= 0.5)
        & (Model["sdfor"][:] <= 250.0)
        & (Data["SCF_clim"][:] >= 0.0)
        & ((Data["SCF" + Time00][:] >= scf_threshold) | (Data["SCF_clim"][:] <= qc_threshold))
    )[0]

    if np.size(p_screen) > 0:

        qc_monit = Model["SD_control"][:] * 0.0

        # Background and observations errors
        # --------------------------------------------------------

        # Background errors
        B_ERR = 0.03

        Bmat = np.array(
            [
                [B_ERR ** 2],
            ],
            dtype="float32",
        )

        # Observation errors
        R_ERR_SD = 0.04

        # Observation-error covariance matrix (diagonal)
        Rmat = np.array(
            [
                [R_ERR_SD ** 2],
            ],
            dtype="float32",
        )

        # Observation operator
        # --------------------------------------------------------

        # Model values in observation space
        hx = np.array(
            [
                [Model["SD_background"][:, :]],
            ],
            dtype="float32",
        )

        # Observations themselves
        Data["SD" + Time00] = Model["SD_background"][0, :]
        p_ims = np.where( (Model["SD_background"][0, :] < 0.01) & (Data["SCF" + Time00][:] >= scf_threshold) )[0]
        Data["SD" + Time00][p_ims] = 0.03
        p_ims = np.where( Data["SCF" + Time00][:] < scf_threshold )[0]
        Data["SD" + Time00][p_ims] = 0.0

        Y00 = np.array(
            [
                [Data["SD" + Time00][:]],
            ],
            dtype="float32",
        )

        # Observation operator Jacobians
        # --------------------------------------------------------

        H = dict()
        H[Time00] = np.array(
            [
                [Model["SD_control"][:] * 0.0 + 1.0],
            ],
            dtype="float32",
        )

        Hmat = H[Time00]
        Htmat = Hmat.transpose()

        # Setup matrices for points
        # --------------------------------------------------------

        innov_mat = dict()
        innov_mat[Time00] = np.array(Y00[:, 0, p_screen] - hx[:, 0, 0, p_screen]) * 0.0

        p_dep = dict()
        p_dep[Time00] = np.where(
            (abs(Y00[0, 0, p_screen] - hx[0, 0, 0, p_screen]) <= 5.0*np.sqrt(B_ERR**2 + R_ERR_SD**2))
        )[0]

        innov_mat[Time00][:, p_dep[Time00]] = np.array(
            Y00[:, 0, p_screen[p_dep[Time00]]] - hx[:, 0, 0, p_screen[p_dep[Time00]]]
        )
        innov_all = innov_mat[Time00]

        INC, K_gain = sekf_functions.Inc_calc_ecmwf(
            p_screen,
            Bmat[:, :],
            Hmat[:, :, p_screen],
            Htmat[p_screen, :, :],
            Rmat[:, :],
            innov_all,
        )

        print("Increments min max average:")
        print(min(INC[:, 0]), max(INC[:, 0]), np.mean(abs(INC[:, 0])))

        # Save increments to full global array before writing out:
        increment[p_screen] = INC[:, 0]

        # QC by background screen-level temperature
        p_qc = np.where(
            (increment[:] > 0.0)
            & (increment[:] > 0.16 - Data["T2m" + Time00][:]*0.016)
        )[0]
        increment[p_qc] = 0.0
        qc_monit[p_qc] = 1.0

        if LSAVE_SEKF == "true":
            Model["SEKF_save"][3, :] = increment[:]
            Model["SEKF_save"][4, :] = qc_monit[:]

        # Modify snow density
        p_inc = np.where(
            (increment[:] > 0.0)
            & (Model["SD_control"][:] == 0.0)
        )[0]
        if np.size(p_inc) > 0:
            freshdens = Model["SD_control"][:] * 0.0
            freshdens[p_inc] = 109 + 6.0*Data["T2m" + Time00][p_inc] \
                             + 26.0*np.power(Data["U10m" + Time00][p_inc]**2 + Data["V10m" + Time00][p_inc]**2, 0.25)
            Model["snowdens_control"][p_inc] = np.clip(freshdens[p_inc],100.0,450.0)
            if LESNML == "true":
                Model["snowdensML_control"][0, p_inc] = np.clip(freshdens[p_inc],100.0,450.0)

        # Add increments to each layer
        if LESNML == "true":
            p_inc = np.where(
                (increment[:] > 0.0)
            )[0]
            if np.size(p_inc) > 0:
                Model["SWEML_control"][0, p_inc] = Model["SWEML_control"][0, p_inc] + increment[p_inc]*Model["snowdensML_control"][0, p_inc]

            kmax = Model["SD_control"][:] * 0.0
            for k in range(NCSNEC):
                p_max = np.where(
                    (Model["SWEML_control"][k, :] <= 0.0)
                )[0]
                kmax[p_max] = 1.0

                p_inc = np.where(
                    (increment[:] < 0.0)
                    & (kmax[:] == 0.0)
                )[0]
                if np.size(p_inc) > 0:
                    wk = Model["SWEML_control"][k, p_inc] / Model["SD_control"][p_inc]
                    Model["SWEML_control"][k, p_inc] = np.clip(Model["SWEML_control"][k, p_inc]+wk*increment[p_inc],0.0,None)
                    p_inc = np.where(
                        (increment[:] < 0.0)
                        & (kmax[:] == 0.0)
                        & (Model["SWEML_control"][k, :] < 0.000001)
                    )[0]
                    if np.size(p_inc) > 0:
                        Model["SWEML_control"][k, p_inc] = 0.0
                        Model["snowdensML_control"][k, p_inc] = 100.0

        p_inc = np.where(
            (increment[:] != 0.0)
        )[0]
        if np.size(p_inc) > 0:
            Model["SD_control"][p_inc] = np.clip(Model["SD_control"][p_inc]+increment[p_inc],0.0,None)
            Model["SWE_control"][p_inc] = np.clip(Model["SWE_control"][p_inc]+increment[p_inc]*Model["snowdens_control"][p_inc],0.0,None)

    # Lower limit for SWE
    p_inc = np.where(
        (Model["glacier_mask"][:] <= 0.5)
        & (Model["SD_control"][:] > 0.0)
        & (Model["SD_control"][:] < 0.001)
    )[0]
    if np.size(p_inc) > 0:
        if LESNML == "true":
            Model["SWEML_control"][0, p_inc] = 0.0
            Model["snowdensML_control"][0, p_inc] = 100.0
        Model["SD_control"][p_inc] = 0.0
        Model["SWE_control"][p_inc] = 0.0
        Model["snowdens_control"][p_inc] = 100.0

    # Upper limit for SWE
    SD_CAP = 3.0

    p_inc = np.where(
        (Model["glacier_mask"][:] <= 0.5)
        & (Model["SD_control"][:] > SD_CAP)
        & (Model["SD_control"][:] < 30.0)
    )[0]
    if np.size(p_inc) > 0:
        if LESNML == "true":
            Model["SWEML_control"][:, p_inc] = Model["SWEML_control"][:, p_inc] * SD_CAP / Model["SD_control"][p_inc]
        Model["SD_control"][p_inc] = SD_CAP
        Model["SWE_control"][p_inc] = SD_CAP * Model["snowdens_control"][p_inc]

    # Close netCDF file
    # --------------------------------------------------------
    Model["SWE_control_f"].close()
    Model["snowdens_control_f"].close()

    if LESNML == "true":
        Model["SWEML_control_f"].close()
        Model["snowdensML_control_f"].close()

    if LSAVE_SEKF == "true":
        Model["SEKF_save"][5, :] = Model["SD_control"][:]
        Model["SEKF_save_f"].close()

