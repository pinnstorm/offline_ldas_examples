#!/usr/bin/env python3
#   Python code to run the simplified extended Kalman filter
#
#  Author: original: D. Fairbairn (31/01/2020)
#
#  1. Reads in model output and Jacobians/SLV obs/ascat obs
#    (The observations and SEKF Jacobians are pre-extracted in netCDF4 files over the land regions)
#  2. The SEKF assimilates SLV pseudo-obs and ASCAT at T+0 and T+6 and the analysis increment is added at the
#    end of the assimilation window (T+12)

import os
import sys, getopt
import numpy as np
from netCDF4 import Dataset
import sekf_functions
import copy

VAR_TOL = 1e-5

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
            "t:e:s:a:j:f:b:i:",
            ["oTime00=", "oTime06=", "oLSAVE_SEKF=", "oLUSE_ASCAT=", "oLUSE_EDA_JACOB=", "ofreq=", "oLASCAT_ADAPTIVE_BC", "oINI_BIAS"],
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
        elif opt in ("-a", "--oLUSE_ASCAT"):
            LUSE_ASCAT = str(arg)
        elif opt in ("-j", "--oLUSE_EDA_JACOB"):
            LUSE_EDA_JACOB = str(arg)
        elif opt in ("-f", "--ofreq"):
            freq = float(arg)
        elif opt in ("-b", "--oLASCAT_ADAPTIVE_BC"):
            LASCAT_ADAPTIVE_BC = str(arg)
        elif opt in ("-i", "--oINI_BIAS"):
            INI_BIAS = str(arg)

    nsteps=int(12.0/freq)

    s_freq=int(1.0/freq)

    Times = [Time00, Time06]

    #   If operational config, set ascat obs and background errors to operational values
    #   Otherwise they will use ASCAT noise and AWC respectively
    #    --------------------------------------------------------
    oper_config = True

    #    Read in model and observations from netCDF files
    #    --------------------------------------------------------
    SLV_assim = True

    Control_vars = ["swvl1", "swvl2", "swvl3"]
    Obs = ["ssm", "T2m", "RH2m"]
    Model = dict()
    Data = dict()

    Model["SM_control_f"] = Dataset("SoilMoist_control.nc", "r+")
    Model["SM_control"] = Model["SM_control_f"].variables["SoilMoist"]

    if LSAVE_SEKF == "true":
        Model["SEKF_save_f"] = Dataset("SEKF_save.nc", "r+")
        Model["SEKF_save"] = Model["SEKF_save_f"].variables["T2m"]

    Model["SM_background_f"] = Dataset("ssm_forecast.nc", "r")
    Model["SM_background"] = Model["SM_background_f"].variables["SoilMoist"][0:nsteps:s_freq]
    Model["SM_background_f"].close()

    #    for SLV_assim (if applicable):
    Model["T2m_background_f"] = Dataset("T2m_forecast.nc", "r")
    Model["T2m_background"] = Model["T2m_background_f"].variables["T2m"][(s_freq-1):nsteps+1:s_freq]
    Model["T2m_background_f"].close()

    Model["RH2m_background_f"] = Dataset("RH2m_forecast.nc", "r")
    lat_vals = Model["RH2m_background_f"]["lat"][:]
    lon_vals = Model["RH2m_background_f"]["lon"][:]
    Model["RH2m_background"] = Model["RH2m_background_f"].variables["RH2m"][(s_freq-1):nsteps+1:s_freq]
    Model["RH2m_background_f"].close()

    Model["T2m_prev_f"] = Dataset("T2m_prev.nc", "r")
    Model["T2m_prev"] = Model["T2m_prev_f"].variables["T2m"][(s_freq-1):nsteps+1:s_freq]
    Model["T2m_prev_f"].close()

    Model["RH2m_prev_f"] = Dataset("RH2m_prev.nc", "r")
    Model["RH2m_prev"] = Model["RH2m_prev_f"].variables["RH2m"][(s_freq-1):nsteps+1:s_freq]
    Model["RH2m_prev_f"].close()

    # Set T2m and RH2m to start from T+0
    Model["T2m_background"][1:, :] = Model["T2m_background"][:-1, :]
    Model["RH2m_background"][1:, :] = Model["RH2m_background"][:-1, :]
    Model["T2m_background"][0, :] = Model["T2m_prev"][-1, :]
    Model["RH2m_background"][0, :] = Model["RH2m_prev"][-1, :]


    if LUSE_ASCAT == "true":
        Model["AWC_f"] = Dataset("AWC.nc", "r")
        Model["AWC"] = Model["AWC_f"].variables["AWC"][0]
        Model["AWC_f"].close()

    Model["swe_f"] = Dataset("swe.nc", "r")
    Model["swe"] = Model["swe_f"].variables["SWE"][0]
    Model["swe_f"].close()

    for tt in Times:

        if LUSE_ASCAT == "true":
            Data["Asc" + tt + "_f"] = Dataset("Ascat" + tt + ".nc", "r")
            Data["Asc" + tt] = Data["Asc" + tt + "_f"].variables["Ascat" + tt][0]
            Data["Asc" + tt + "NN_f"] = Dataset("N_Ascat" + tt + ".nc", "r")
            Data["Asc" + tt + "NN"] = Data["Asc" + tt + "NN_f"].variables[
                "N_Ascat" + tt
            ][0]
            Data["Asc" + tt + "HH_f"] = Dataset("H_Ascat" + tt + ".nc", "r")
            Data["Asc" + tt + "HH"] = Data["Asc" + tt + "HH_f"].variables[
                "H_Ascat" + tt
            ][0]

            Data["Asc" + tt + "_f"].close()
            Data["Asc" + tt + "NN_f"].close()
            Data["Asc" + tt + "HH_f"].close()

        #     if SLV_assim:
        Data["T2m" + tt + "_f"] = Dataset("T2m" + tt + ".nc", "r")
        Data["T2m" + tt] = Data["T2m" + tt + "_f"].variables["T2m" + tt][:].flatten()
        Data["RH2m" + tt + "_f"] = Dataset("RH2m" + tt + ".nc", "r")
        Data["RH2m" + tt] = Data["RH2m" + tt + "_f"].variables["RH2m" + tt][:].flatten()

        #      if SLV_assim:
        Data["T2m" + tt + "_f"].close()
        Data["RH2m" + tt + "_f"].close()

        #     Conversion of dew-point temperature obs to relative humidity:
        #     *** Moved this to sekf_extract_RH2m_T2m.ksh in Mars commands ***
        # r2es = 611.21
        # r3les = 17.502
        # r4les = 32.19
        # T0 = 273.16

        # Data["RH2m" + tt] = (
        #     (
        #         r2es
        #         * np.exp(r3les * (Data["RH2m" + tt] - T0) / (Data["RH2m" + tt] - r4les))
        #     )
        #     / (
        #         r2es
        #         * np.exp(r3les * (Data["T2m" + tt] - T0) / (Data["T2m" + tt] - r4les))
        #     )
        # ) * 100.0


    Bqc = 0.1  # Reject ASCAT SM departures exceeding 0.1 m3/m3
    Model["SEKF_save"][:, :] = 0.0

    # Increment calculation
    # -----------------------------------------------------------------------------------------------------------------

    # Model values in observation space
    hx = np.array(
        [
            [Model["T2m_background"][:, :]],
            [Model["RH2m_background"][:, :]],
            [Model["SM_background"][:, 0, :] / 70.0],
        ],
        dtype="float32",
    )

    # Observations themselves

    if LUSE_ASCAT == "true":
        Y00 = np.array(
            [
                [Data["T2m" + Times[0]]],
                [Data["RH2m" + Times[0]]],
                [Data["Asc" + Times[0]]],
            ],
            dtype="float32",
        )
        Y06 = np.array(
            [
                [Data["T2m" + Times[1]]],
                [Data["RH2m" + Times[1]]],
                [Data["Asc" + Times[1]]],
            ],
            dtype="float32",
        )
    else:
        Y00 = np.array(
            [
                [Data["T2m" + Times[0]]],
                [Data["RH2m" + Times[0]]],
                [Data["RH2m" + Times[0]] * 0.0],
            ],
            dtype="float32",
        )
        Y06 = np.array(
            [
                [Data["T2m" + Times[1]]],
                [Data["RH2m" + Times[1]]],
                [Data["RH2m" + Times[1]] * 0.0],
            ],
            dtype="float32",
        )

    #    Background and observations errors
    #    --------------------------------------------------------

    #   As in ECMWF operations,
    B_ERR = 0.02

    Bmat = np.array(
        [[B_ERR ** 2, 0.0, 0.0], [0.0, B_ERR ** 2, 0.0], [0.0, 0.0, B_ERR ** 2]],
        dtype="float32",
    )

    # Test a depth varying clim-B derived from EDA SM spread but using B_ERR
    # Bmat = (B_ERR**2 * np.diag([1., 0.7, 0.6])).astype("float32")

    # Ascat observation errors
    R_ERR_asc00 = 0.0025
    R_ERR_asc06 = R_ERR_asc00

    # T2m/RH2m Observation errors
    R_ERR_T2m = 1.0
    R_ERR_RH2m = 16.0

    # Observation-error covariance matrix (diagonal)
    Rmat = np.array(
        [
            [R_ERR_T2m, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, R_ERR_RH2m, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, R_ERR_asc00, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, R_ERR_T2m, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, R_ERR_RH2m, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, R_ERR_asc06],
        ],
        dtype="float32",
    )

    #    Observation operator Jacobians
    #    --------------------------------------------------------
    Hmat = np.load("HMAT.npy")
    Htmat = Hmat.transpose()  # Hmat.transpose()

    #    Screening and increment calculation
    # -----------------------------------------------------------------------------------------------------------------   

    if LUSE_EDA_JACOB == "true":

        var_swvl1_time00 = Dataset(f"var_{Time00}_swvl1.nc", "r")[f"var_{Time00}_swvl1"][0]
        var_swvl1_time06 = Dataset(f"var_{Time06}_swvl1.nc", "r")[f"var_{Time06}_swvl1"][0]
        var_swvl2_time00 = Dataset(f"var_{Time00}_swvl2.nc", "r")[f"var_{Time00}_swvl2"][0]
        var_swvl2_time06 = Dataset(f"var_{Time06}_swvl2.nc", "r")[f"var_{Time06}_swvl2"][0]
        var_swvl3_time00=np.load("variance_swvl3_"+Time00+".npy")
        var_swvl3_time06=np.load("variance_swvl3_"+Time06+".npy")

        p_screen = np.where(
            (Y00[0, 0, :] > 274.0)
            & (Y06[0, 0, :] > 274.0)
            & (Model["swe"][:] < 0.01)
            & (Model["SM_background"][0, 0, :] > 0.0)
            & (abs(Hmat[2, 0, :]) < 2.0)
            & (abs(Hmat[2, 1, :]) < 2.0)
            & (abs(Hmat[2, 2, :]) < 2.0)
            & (abs(Hmat[0, 0, :]) < 50.0)
            & (abs(Hmat[0, 1, :]) < 50.0)
            & (abs(Hmat[0, 2, :]) < 50.0)
            & (abs(Hmat[1, 0, :]) < 500.0)
            & (abs(Hmat[1, 1, :]) < 500.0)
            & (abs(Hmat[1, 2, :]) < 500.0)
            & (abs(Hmat[5, 0, :]) < 2.0)
            & (abs(Hmat[5, 1, :]) < 2.0)
            & (abs(Hmat[5, 2, :]) < 2.0)
            & (abs(Hmat[3, 0, :]) < 50.0)
            & (abs(Hmat[3, 1, :]) < 50.0)
            & (abs(Hmat[3, 2, :]) < 50.0)
            & (abs(Hmat[4, 0, :]) < 500.0)
            & (abs(Hmat[4, 1, :]) < 500.0)
            & (abs(Hmat[4, 2, :]) < 500.0)
            & (var_swvl1_time00 > VAR_TOL)
            & (var_swvl1_time06 > VAR_TOL)
            & (var_swvl2_time00 > VAR_TOL)
            & (var_swvl2_time06 > VAR_TOL)
            & (var_swvl3_time00 > VAR_TOL)
            & (var_swvl3_time06 > VAR_TOL)
        )[0]
    else:
        p_screen = np.where(
            (Y00[0, 0, :] > 274.0)
            & (Y06[0, 0, :] > 274.0)
            & (Model["swe"][:] < 0.01)
            & (Model["SM_background"][0, 0, :] > 0.0)
            & (abs(Hmat[2, 0, :]) < 2.0)
            & (abs(Hmat[2, 1, :]) < 2.0)
            & (abs(Hmat[2, 2, :]) < 2.0)
            & (abs(Hmat[0, 0, :]) < 50.0)
            & (abs(Hmat[0, 1, :]) < 50.0)
            & (abs(Hmat[0, 2, :]) < 50.0)
            & (abs(Hmat[1, 0, :]) < 500.0)
            & (abs(Hmat[1, 1, :]) < 500.0)
            & (abs(Hmat[1, 2, :]) < 500.0)
            & (abs(Hmat[5, 0, :]) < 2.0)
            & (abs(Hmat[5, 1, :]) < 2.0)
            & (abs(Hmat[5, 2, :]) < 2.0)
            & (abs(Hmat[3, 0, :]) < 50.0)
            & (abs(Hmat[3, 1, :]) < 50.0)
            & (abs(Hmat[3, 2, :]) < 50.0)
            & (abs(Hmat[4, 0, :]) < 500.0)
            & (abs(Hmat[4, 1, :]) < 500.0)
            & (abs(Hmat[4, 2, :]) < 500.0)
        )[0]


    #    Setup matrices for points
    # -----------------------------------------------------------------------------------------------------------------

    innov_mat = dict()

    innov_mat[Time00] = np.array(Y00[:, 0, p_screen] - hx[:, 0, 0, p_screen]) * 0.0
    innov_mat[Time06] = np.array(Y06[:, 0, p_screen] - hx[:, 0, 6, p_screen]) * 0.0

    p_slv = dict()
    if SLV_assim:

        # Points where SLV passes quality control:
        p_slv[Time00] = np.where(
            (abs(Y00[0, 0, p_screen] - hx[0, 0, 0, p_screen]) < 5.0)
            & (abs(Y00[1, 0, p_screen] - hx[1, 0, 0, p_screen]) < 20.0)
        )[0]
        innov_mat[Time00][0:2, p_slv[Time00]] = np.array(
            Y00[0:2, 0, p_screen[p_slv[Time00]]]
            - hx[0:2, 0, 0, p_screen[p_slv[Time00]]]
        )

        p_slv[Time06] = np.where(
            (abs(Y06[0, 0, p_screen] - hx[0, 0, 6, p_screen]) < 5.0)
            & (abs(Y06[1, 0, p_screen] - hx[1, 0, 6, p_screen]) < 20.0)
        )[0]
        innov_mat[Time06][0:2, p_slv[Time06]] = np.array(
            Y06[0:2, 0, p_screen[p_slv[Time06]]]
            - hx[0:2, 0, 6, p_screen[p_slv[Time06]]]
        )

    p_ascat = dict()
    if LUSE_ASCAT == "true":

        # Points where ascat passes quality control:
        Bqc = 0.1
        p_ascat[Time00] = np.where(
            (abs(Y00[2, 0, p_screen] - hx[2, 0, 0, p_screen]) < Bqc)
        )[0]
        innov_mat[Time00][2, p_ascat[Time00]] = np.array(
            Y00[2, 0, p_screen[p_ascat[Time00]]]
            - hx[2, 0, 0, p_screen[p_ascat[Time00]]]
        )

        p_ascat[Time06] = np.where(
            (abs(Y06[2, 0, p_screen] - hx[2, 0, 6, p_screen]) < Bqc)
        )[0]
        innov_mat[Time06][2, p_ascat[Time06]] = np.array(
            Y06[2, 0, p_screen[p_ascat[Time06]]]
            - hx[2, 0, 6, p_screen[p_ascat[Time06]]]
        )

    innov_all = np.concatenate((innov_mat[Time00], innov_mat[Time06]), axis=0)

    #    Setup ascat bias correction
    # -----------------------------------------------------------------------------------------------------------------
    if (LASCAT_ADAPTIVE_BC == "true") and (LUSE_ASCAT == "true"):

      Bias_file = Dataset("sekf_ascat_bias_in.nc", "r+")
      Model["SM_bias"] = Bias_file.variables["sekf_ascat_bias_in"]

      if INI_BIAS == "coldstart":
        Model["SM_bias"][0, :] = Model["SM_bias"][0, :] * 0.0
      elif INI_BIAS == "standard":
        Model["SM_bias"][0, :] = 0.0
        Model["SM_bias"][0, p_screen[p_ascat[Time00]]] = (innov_mat[Time00][2, p_ascat[Time00]])
        Model["SM_bias"][0, p_screen[p_ascat[Time06]]] = (innov_mat[Time06][2, p_ascat[Time06]])

      bias_innov = dict()

      SM1 = np.array(Model["SM_bias"][0, :])      
      SM1[SM1 > 0.1] = 0.1
      SM1[SM1 < -0.1] = -0.1
      Model["SM_bias"][0, :] = SM1

      bias_innov[Time00] = copy.deepcopy(innov_mat[Time00])
      bias_innov[Time06] = copy.deepcopy(innov_mat[Time06])

      bias_innov[Time00][2, p_ascat[Time00]] -= Model["SM_bias"][
      0, p_screen[p_ascat[Time00]]
      ]
      bias_innov[Time06][2, p_ascat[Time06]] -= Model["SM_bias"][
      0, p_screen[p_ascat[Time06]]
      ]

      bias_innov_all = np.concatenate((bias_innov[Time00], bias_innov[Time06]), axis=0)

      INC, K_gain = sekf_functions.Inc_calc_ecmwf(
        p_screen,
        Bmat[:, :],
        Hmat[:, :, p_screen],
        Htmat[p_screen, :, :],
        Rmat[:, :],
        bias_innov_all,
      )

      bias_innov_ascat = np.concatenate((bias_innov[Time00][2:3, :], bias_innov[Time06][2:3, :]), axis=0)

      gamma = 0.25

      Hmat_ascat=np.concatenate((Hmat[2:3, :, p_screen], Hmat[5:6, :, p_screen]),axis=0)
      Htmat_ascat=Hmat_ascat.transpose()


      Rmat_ascat=np.array(
        [
            [R_ERR_asc00, 0.0],
            [0.0, R_ERR_asc06],
        ],
        dtype="float32",
    )

      print(np.shape(bias_innov_ascat),np.shape(bias_innov),np.shape(Hmat_ascat),np.shape(Htmat_ascat),np.shape(Hmat),np.shape(Htmat))

      # Bias calculation:
      INC_beta = sekf_functions.Inc_calc_bias(
        p_screen,
        Bmat[:, :],
        Hmat_ascat,
        Htmat_ascat,
        Rmat_ascat[:, :],
        bias_innov_ascat,
        gamma,
      )

      increments_b = Model["SM_bias"][0, :] * 0.0
      increments_b[p_screen] = INC_beta[:]


      Model["SM_bias"][0, :] += increments_b[:]

      print("Obs bias min max average:")
      print(
        min(Model["SM_bias"][0, :]),
        max(Model["SM_bias"][0, :]),
        np.mean(Model["SM_bias"][0, :]),
      )

    else:

      bias_innov=None
      INC, K_gain = sekf_functions.Inc_calc_ecmwf(
        p_screen,
        Bmat[:, :],
        Hmat[:, :, p_screen],
        Htmat[p_screen, :, :],
        Rmat[:, :],
        innov_all,
      )

    print("Increments min max average (layers 1-3):")
    print(min(INC[:, 0]), max(INC[:, 0]), np.mean(abs(INC[:, 0])))
    print(min(INC[:, 1]), max(INC[:, 1]), np.mean(abs(INC[:, 1])))
    print(min(INC[:, 2]), max(INC[:, 2]), np.mean(abs(INC[:, 2])))

    # Save increments to full global array before writing out:
    increments = Model["SM_control"][:] * 0.0
    increments[0, p_screen] = INC[:, 0]
    increments[1, p_screen] = INC[:, 1]
    increments[2, p_screen] = INC[:, 2]

    if LSAVE_SEKF == "true":

        print(len(INC[:, 0][abs(INC[:, 0]) > 0.1]))
        INC[:, 0][abs(INC[:, 0]) > 0.1] = 0.0
        print(len(INC[:, 1][abs(INC[:, 1]) > 0.1]))
        INC[:, 1][abs(INC[:, 1]) > 0.1] = 0.0
        print(len(INC[:, 2][abs(INC[:, 2]) > 0.1]))
        INC[:, 2][abs(INC[:, 2]) > 0.1] = 0.0

        # Save first guess departures
        innovations = Model["SEKF_save"][0:6, :] * 0.0
        innovations[0:6, p_screen] = innov_all[0:6, :]
        Model["SEKF_save"][0:6, :] = innovations[0:6, :]

        # Save analysis increments
        Model["SEKF_save"][6:9, :] = increments[0:3, :]

    if (LASCAT_ADAPTIVE_BC == "true"):
        # Save bias correction terms
        Model["SEKF_save"][9, :] = Model["SM_bias"][0, :]
        Model["SEKF_save"][10, :] = increments_b[:]
    else:

        #Save Kalman gain G0101 to G0103
        kalman_gain = Model["SEKF_save"][9:12, :] * 0.0
        kalman_gain[:, p_screen] = K_gain[:, 0, 0:3].swapaxes(0, 1)
        Model["SEKF_save"][9:12, :] = kalman_gain[:, :]

        #       Close netCDF file
        #       --------------------------------------------------------
        Model["SEKF_save_f"].close()

    #  Soil moisture increment added to control variable final timestep for initializing following assimilation window

    if (LUSE_ASCAT == "true"):

      innov_mat[Time00+'_an_depar'] = np.array(Y00[2, 0, p_screen] - hx[2, 0, 0, p_screen]) * 0.0
      innov_mat[Time06+'_an_depar'] = np.array(Y06[2, 0, p_screen] - hx[2, 0, 6, p_screen]) * 0.0

      if (LASCAT_ADAPTIVE_BC == "true"):
        innov_mat[Time00+'_an_depar'][p_ascat[Time00]] = bias_innov[Time00][2, p_ascat[Time00]] - increments[0,p_screen[p_ascat[Time00]]]
        innov_mat[Time06+'_an_depar'][p_ascat[Time06]] = bias_innov[Time06][2, p_ascat[Time06]] - increments[0,p_screen[p_ascat[Time06]]]
      else:
        innov_mat[Time00+'_an_depar'][p_ascat[Time00]] = innov_mat[Time00][2, p_ascat[Time00]] - increments[0,p_screen[p_ascat[Time00]]]
        innov_mat[Time06+'_an_depar'][p_ascat[Time06]] = innov_mat[Time06][2, p_ascat[Time06]] - increments[0,p_screen[p_ascat[Time06]]]

    # Define wilting point array for screening 
    soil_type_arr = Dataset("surfclim", "r").variables["sotype"][:]
    soil_dic = {0: np.nan,
                1: 0.059,
                2: 0.151,
                3: 0.133,
                4: 0.279,
                5: 0.335,
                6: 0.267,
                7: 0.151}
    zwilt_arr = np.array([soil_dic[x] for x in soil_type_arr[:]])
    # Find indices of swvl2/3 < z_wilt
    swvl2_idx  = np.where(Model["SM_background"][0, 1, :] < zwilt_arr)[0]
    swvl3_idx  = np.where(Model["SM_background"][0, 2, :] < zwilt_arr)[0]
    # Turn off SLV assimilation where swvl2/3 < z_wilt by setting increment elements to zero
    increments[1, swvl2_idx] = 0
    increments[2, swvl3_idx] = 0

    increments[0, :] = increments[0, :] * 70.0
    increments[1, :] = increments[1, :] * 210.0
    increments[2, :] = increments[2, :] * 720.0

    Model["SM_control"][0:3, :] += increments[0:3, :]

    SM1 = np.array(Model["SM_control"][0, :])
    SM2 = np.array(Model["SM_control"][1, :])
    SM3 = np.array(Model["SM_control"][2, :])

    #  Ensure soil moisture is higher than residual value and lower than saturated value

    SM1 = [
        0.01 if x < 0.01 else (0.7660 * 70.0) if x > (0.7660 * 70.0) else x for x in SM1
    ]
    SM2 = [
        0.01 if x < 0.01 else (0.7660 * 210.0) if x > (0.7660 * 210.0) else x
        for x in SM2
    ]
    SM3 = [
        0.01 if x < 0.01 else (0.7660 * 720.0) if x > (0.7660 * 720.0) else x
        for x in SM3
    ]

    Model["SM_control"][0, :] = SM1
    Model["SM_control"][1, :] = SM2
    Model["SM_control"][2, :] = SM3

    sekf_functions.summarise(
        lat_vals,
        lon_vals,
        p_screen,
        p_slv,
        p_ascat,
        innov_mat,
        Data,
        Times,
        LUSE_ASCAT,
        SLV_assim,
        LASCAT_ADAPTIVE_BC,
        bias_innov,
        Model
    )

    #  Close open files
    Model["SM_control_f"].close()
    if (LASCAT_ADAPTIVE_BC == "true"):
      Bias_file.close()

