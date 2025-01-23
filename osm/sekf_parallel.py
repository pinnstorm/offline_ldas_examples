#!/usr/bin/env python3
import argparse
import logging

import dask
import netCDF4 as nc
import numpy as np
import xarray as xr

logging.basicConfig(level=logging.INFO)


CTRLVAR = ["swvl1", "swvl2", "swvl3"]
BACK_ERR = 0.01
WORK_DIR = "/ec/res4/scratch/daep/offline_hres_surf_ana/"
# WORK_DIR = ""

OBS_convdic = {"T2m": "T2m", "RH2m": "RH2m", "ssm": "Ascat"}
OBS_ERR = {"ssm": 0.0025, "RH2m": 16.0, "T2m": 1.0}

EDAH_TAPER = 0.6
HMAT_CHUNK_SIZE = -1  # 3800000
SEKF_CHUNK_SIZE = -1  # 100000

freq = 1
nsteps = int(12.0 / freq)
s_freq = int(1.0 / freq)


def conv_td2rh(t2m_da: xr.DataArray, rh2m_da: xr.DataArray) -> xr.DataArray:
    """Function to convert dew point temperature to relative humidity

    :param t2m_da: 2m temperature dataarray
    :param rh2m_da: 2m dew point temperature dataarray
    :return: 2m relative humidtiy dataarray
    """
    r2es = 611.21
    r3les = 17.502
    r4les = 32.19
    T0 = 273.16
    num = r2es * np.exp(r3les * (rh2m_da - T0) / (rh2m_da - r4les))
    den = r2es * np.exp(r3les * (t2m_da - T0) / (t2m_da - r4les))
    return 100 * (num / den)


@dask.delayed
def open_nc_da(fname: str, var_name: str, t_idx: int = None) -> xr.DataArray:
    """Dask delayed function that opens a netcdf file using xarray

    :param fname: path to file
    :param var_name: name of variable to open
    :param t_idx: index along the time deimension to open, defaults to None
    :return: dataarray for give path and variable
    """
    da = xr.open_dataset(fname, engine="netcdf4")[var_name]
    if "x" in da.coords.keys():
        da = da.drop_vars("x")
    if "time" in da.coords.keys():
        da = da.drop_vars("time")
    if t_idx is not None:
        da = da.isel(time=t_idx)
    return da


def eda_jacobian_calc(
    var_ds: xr.Dataset,
    covar_ds: xr.Dataset,
    obs: list = ["ssm", "RH2m", "T2m"],
    tt: str = "00",
    EDAH_TAPER: float = 0.6,
) -> xr.DataArray:
    """Function that calculates the jacobian of the observation operator using variances and
    covariances from the Ensemble of Data Assimilations (EDA)

    :param var_ds: Dataset containing soil moisture variances at three depths
    :param covar_ds: Dataset containing covariances between soil moisture at three depths and the assimilated observations
    :param tt: time step at which to calculated jacobian, defaults to "00"
    :param EDAH_TAPER: taper coefficient to impact of observations at depth, defaults to 0.6
    :return: dataarray of observation operator jacobian
    """
    hmat = np.array(
        [
            [
                covar_ds[f"covar_{ob}_{tt}_swvl1"] / var_ds[f"var_{tt}_swvl1"],
                (1.0 / (1.0 + EDAH_TAPER))
                * covar_ds[f"covar_{ob}_{tt}_swvl2"]
                / var_ds[f"var_{tt}_swvl2"],
                (1.0 / (1.0 + (2.0 * EDAH_TAPER)))
                * covar_ds[f"covar_{ob}_{tt}_swvl3"]
                / var_ds[f"var_{tt}_swvl3"],
            ]
            for ob in obs
        ],
        dtype="float32",
    )
    hmat_da = xr.DataArray(
        np.array(hmat),
        dims=["obs", "ctrlvec", "x"],
        name="hmat",
    )
    return hmat_da


def calc_inc_sekf(
    h_mat: xr.DataArray, innov: xr.DataArray, bmat_inv: np.ndarray, rmat_inv: np.ndarray
) -> xr.DataArray:
    """Function calculating incrememnts for the simplified extended Kalman filter (SEKF)

    :param h_mat: Linearised observation operator
    :param innov: innovation dataarray of (y - hxb)
    :param bmat_inv: inverse of the background error covariance matrix
    :param rmat_inv: inverse of the observation error covariance matrix
    :return: increments for the control vector at t0
    """
    h_tmp = h_mat.transpose(
        "obs",
        "ctrlvec",
        "x",
    ).values
    #  R^{-1}H
    invR_H = np.einsum("ab,bdi->iad", rmat_inv, h_tmp)
    #  H^TR^{-1}H
    Ht_invR_H = np.einsum("iab,ibd->iad", h_tmp.T, invR_H)
    #  B^{-1}+H^TR^{-1}H
    invB_Ht_invR_H = Ht_invR_H + bmat_inv
    #  (B^{-1}+H^TR^{-1}H)^{-1}
    inv_invB_Ht_invR_H = np.linalg.inv(invB_Ht_invR_H)
    #  H^TR^{-1}
    Ht_invR = np.einsum("iab,bd->iad", h_tmp.T, rmat_inv)
    #  K1 = (B^{-1}+H^TR^{-1}H)^{-1}H^TR^{-1}
    K1 = np.einsum("iab,ibd->iad", inv_invB_Ht_invR_H, Ht_invR)
    #  dx = K1(yo-hx)
    Zinc = np.einsum("iab,bi->ia", K1, innov.values)
    increments = xr.DataArray(
        np.array(Zinc).T.reshape(h_mat.isel(obs=0, drop=True).shape),
        coords=h_mat.isel(obs=0, drop=True).coords,
        dims=h_mat.isel(obs=0, drop=True).dims,
        name="incs",
    )
    return increments


@dask.delayed
def make_hxb(data: dict) -> xr.Dataset:
    """From dictionary of xr.DataArrays creates model estimate to observations hxb

    :param data: dictionary of model data in obesrvations space
    :return: model in observation space hxb
    """
    hxb = {
        f"ssm{TIMES[0]}": (data["ssm_fc"].isel(time=0, nlevs=0, drop=True) / 70.0),
        f"RH2m{TIMES[0]}": data["rh2m_pfc"],
        f"T2m{TIMES[0]}": data["t2m_pfc"],
        f"ssm{TIMES[1]}": (data["ssm_fc"].isel(time=1, nlevs=0, drop=True) / 70.0),
        f"RH2m{TIMES[1]}": data["rh2m_fc"],
        f"T2m{TIMES[1]}": data["t2m_fc"],
    }
    hxb = xr.Dataset(hxb)
    hxb = hxb[[f"{ob}{tt}" for tt in TIMES for ob in OBS]]
    return hxb


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="SEKF assimilation script for offline LDAS experiments"
    )
    parser.add_argument(
        "-t",
        "--TIME00",
        help="start time of assimlation window (e.g. 00 or 12)",
        required=True,
        default="00",
    )
    parser.add_argument(
        "-e",
        "--TIME01",
        help="end time of assimlation window (e.g. 06 or 18)",
        required=True,
        default="06",
    )
    parser.add_argument(
        "-s",
        "--LSAVE_SEKF",
        help="switch to save SEKF increments",
        default="false",
    )
    parser.add_argument(
        "-a",
        "--LUSE_ASCAT",
        help="switch to use ASCAT SM observations",
        default="true",
    )
    parser.add_argument(
        "-j",
        "--LUSE_EDA_JACOB",
        help="switch to use EDA in Hmat Jacobian creation",
        default="true",
    )
    parser.add_argument(
        "-i",
        "--INITIME",
        help="switch for initial time of experiment",
        default="false",
    )
    args = parser.parse_args()
    logging.info(f"Parsed args = {args}")

    TIMES = [args.TIME00, args.TIME01]

    if args.LUSE_ASCAT == "true":
        logging.info("ASCAT used in assimilation")
        OBS = ["ssm", "RH2m", "T2m"]
        ASC_OBS = [OBS_convdic[ob] for ob in OBS]
    else:
        logging.info("SLV only assimilation")
        OBS = ["RH2m", "T2m"]
        ASC_OBS = [OBS_convdic[ob] for ob in OBS]

    logging.info("Opening necessary data for assimilation cycle...")
    data = dict(
        sm_control=open_nc_da(
            f"{WORK_DIR}SoilMoist_control.nc",
            "SoilMoist",
        ),
        ssm_fc=open_nc_da(
            f"{WORK_DIR}ssm_forecast.nc",
            "SoilMoist",
            [0, 5],
        ),
        t2m_fc=open_nc_da(f"{WORK_DIR}T2m_forecast.nc", "T2m", 5),
        rh2m_fc=open_nc_da(f"{WORK_DIR}RH2m_forecast.nc", "RH2m", 5),
        t2m_pfc=open_nc_da(f"{WORK_DIR}T2m_prev.nc", "T2m", -1),
        rh2m_pfc=open_nc_da(f"{WORK_DIR}RH2m_prev.nc", "RH2m", -1),
        swe=open_nc_da(f"{WORK_DIR}swe.nc", "SWE", 0),
    )
    logging.info("Opening observations to be used in assimilation cycle...")
    ob_dic = {
        f"{ob}{tt}": open_nc_da(f"{WORK_DIR}{ob}{tt}.nc", f"{ob}{tt}", 0)
        for tt in TIMES
        for ob in ASC_OBS
    }
    logging.info(
        "Opening variances/covariances to be used in creation of linearised obs operator Hmat..."
    )
    var_dic = {
        f"var_{tt}_{cv}": open_nc_da(
            f"{WORK_DIR}var_{tt}_{cv}.nc",
            f"var_{tt}_{cv}",
            0,
        )
        for tt in TIMES
        for cv in CTRLVAR
    }
    covar_dic = {
        f"covar_{ob}_{tt}_{cv}": open_nc_da(
            f"{WORK_DIR}covar_{ob}_{tt}_{cv}.nc",
            f"covar_{ob}_{tt}_{cv}",
            0,
        )
        for tt in TIMES
        for ob in OBS
        for cv in CTRLVAR
    }
    logging.info("Opening files in parallel using Dask and load into memory...")
    data, ob_dic, var_dic, covar_dic = dask.compute(data, ob_dic, var_dic, covar_dic)

    logging.info("Creating xarray Datasets from the dictionaries of variables...")
    var_ds = xr.Dataset(var_dic)
    covar_ds = xr.Dataset(covar_dic)
    # rename ascat to ssm to be consistent with hxb and ensure our dimensions are in order
    obs_ds = xr.Dataset(ob_dic)
    if args.LUSE_ASCAT == "true":
        obs_ds = obs_ds.rename({f"Ascat{tt}": f"ssm{tt}" for tt in TIMES})
    logging.info("Converting dew point temperature to relative humidity...")
    for tt in TIMES:
        obs_ds[f"RH2m{tt}"] = conv_td2rh(obs_ds[f"T2m{tt}"], obs_ds[f"RH2m{tt}"])

    # create hxb
    hxb_ds = make_hxb(data)
    # create linearised hxb using Dask map_blocks function to apply eda_jacobian_cal over chunks
    # hmat = var_ds.chunk({"x": HMAT_CHUNK_SIZE}).map_blocks(
    #     eda_jacobian_calc, [covar_ds.chunk({"x": HMAT_CHUNK_SIZE}), "06"]
    # )
    hmat_dic = {
        tt: var_ds.chunk({"x": HMAT_CHUNK_SIZE}).map_blocks(
            eda_jacobian_calc, [covar_ds.chunk({"x": HMAT_CHUNK_SIZE}), OBS, tt]
        )
        for tt in TIMES
    }
    logging.info("Triggering computation of hxb and hmat in parallel with Dask...")
    hxb_ds, hmat = dask.compute(hxb_ds, hmat_dic)
    # Convert Datasets to DataArrays for use in assimilation functions
    obs_da = obs_ds.to_array(dim="obs")
    hxb_da = hxb_ds.to_array(dim="obs")
    # Create Hmat of correct shape as a DataArray (in previous scripts Hmat at time t1 instead of t0 was used for both time, corresponding to an ini_time arg)
    if args.INITIME == "true":
        hmat_da = xr.concat([hmat_dic[TIMES[-1]], hmat_dic[TIMES[-1]]], dim="obs")
    else:
        hmat_da = xr.concat(hmat_dic.values(), dim="obs")
    # Slight discrepencies with previous scripts removed by replacing NaN and inf values with zero, these broadly correspond to coastal or poorly conditioned points, so may want to remove from assimilation anyway?
    hmat_da = hmat_da.where((hmat_da.notnull()) & (~np.isinf(hmat_da)), 999)
    # Set obs dimension of DataArray to have same values as observation DataArray obs_da
    hmat_da["obs"] = obs_da.obs

    logging.info("Finding indicies to screen from computation...")
    if args.LUSE_ASCAT == "true":
        indices = np.where(
            (obs_da.sel(obs=[f"T2m{TIMES[0]}", f"T2m{TIMES[1]}"]).min("obs") > 274.0)
            & (
                abs(hmat_da.sel(obs=[f"T2m{TIMES[0]}", f"T2m{TIMES[1]}"])).max(
                    ("obs", "ctrlvec")
                )
                < 50.0
            )
            & (
                abs(hmat_da.sel(obs=[f"RH2m{TIMES[0]}", f"RH2m{TIMES[1]}"])).max(
                    ("obs", "ctrlvec")
                )
                < 500.0
            )
            & (
                abs(hmat_da.sel(obs=[f"ssm{TIMES[0]}", f"ssm{TIMES[1]}"])).max(
                    ("obs", "ctrlvec")
                )
                < 2.0
            )
            & (data["swe"] < 0.01)
            & (data["ssm_fc"].isel(time=0, nlevs=0, drop=True) > 0.0)
            & (var_dic[f"var_{TIMES[0]}_swvl3"] > 0.0001)
            & (var_dic[f"var_{TIMES[1]}_swvl3"] > 0.0001)
        )[0]
        logging.info(f"Number of indicies found = {len(indices)}")
    else:
        indices = np.where(
            (obs_da.sel(obs=[f"T2m{TIMES[0]}", f"T2m{TIMES[1]}"]).min("obs") > 274.0)
            & (
                abs(hmat_da.sel(obs=[f"T2m{TIMES[0]}", f"T2m{TIMES[1]}"])).max(
                    ("obs", "ctrlvec")
                )
                < 50.0
            )
            & (
                abs(hmat_da.sel(obs=[f"RH2m{TIMES[0]}", f"RH2m{TIMES[1]}"])).max(
                    ("obs", "ctrlvec")
                )
                < 500.0
            )
            & (data["swe"] < 0.01)
            & (data["ssm_fc"].isel(time=0, nlevs=0, drop=True) > 0.0)
            & (var_dic[f"var_{TIMES[0]}_swvl3"] > 0.0001)
            & (var_dic[f"var_{TIMES[1]}_swvl3"] > 0.0001)
        )[0]
        logging.info(f"Number of indicies found = {len(indices)}")

    logging.info("Calcualting innovations (y - hxb)...")
    innov_ds = obs_ds.isel(x=indices) - hxb_ds.isel(x=indices)
    logging.info("Setting any innovations greater than specified tolerance to zero...")
    dep_tol = {"T2m": 5, "RH2m": 20, "ssm": 0.1}
    for key in OBS:
        for hr in [TIMES[0], TIMES[1]]:
            innov_ds[key + hr] = innov_ds[key + hr].where(
                abs(innov_ds[key + hr]) < dep_tol[key], 0
            )

    # Define inverse of Bmat and Rmat
    Bmat_inv = np.linalg.inv(np.diag([BACK_ERR ** 2] * 3))
    Rmat_inv = np.linalg.inv(np.diag([OBS_ERR[ob] for ob in OBS] * 2))

    logging.info(
        "Calculating increments from SEKF using Dask map_blocks function to parallelize computation..."
    )
    inc_template = (
        hmat_da.isel(x=indices).isel(obs=0, drop=True).chunk({"x": SEKF_CHUNK_SIZE})
    )
    innov_da = innov_ds.to_array(dim="obs", name="innov").chunk(
        {"obs": -1, "x": SEKF_CHUNK_SIZE}
    )
    incs = (
        hmat_da.isel(x=indices)
        .chunk({"obs": -1, "x": SEKF_CHUNK_SIZE})
        .map_blocks(
            calc_inc_sekf,
            [
                innov_da,
                Bmat_inv,
                Rmat_inv,
            ],
            template=inc_template,
        )
        .compute()
    )

    # Print some diagnostics
    INC = incs.values
    logging.info("Printing increments min max average for layers 1-3:")
    for x in range(3):
        logging.info(
            f"min: {np.min(INC[x]):.6f}, max: {np.max(INC[x]):.6f}, mu: {np.mean(abs(INC[x])):.6f}"
        )

    if args.LSAVE_SEKF == "true":
        logging.info("Saving SEKF incremenets to file...")
        da = xr.open_dataset(f"{WORK_DIR}T2m_forecast.nc").isel(time=slice(0, 3))
        da["T2m"][:] = 0
        da["T2m"][:, indices] = INC
        logging.info(f"SEKF increments saved under ini time: {da.time.values[0]}")
        da.to_netcdf("SEKF_save.nc", mode="w")

    logging.info("Converting values for different soil levels...")
    for x, nlev_conv in enumerate([70, 210, 720]):
        incs[x] *= nlev_conv
    logging.info("Adding increments to control vector...")
    sm_tmp = data["sm_control"].values
    sm_tmp[:3, indices] += incs

    logging.info(
        "Opening netcdf file to restart model with and updating Soil State with new values..."
    )
    # this was previously done using NCKS but appears to be much faster using netCDF4 while results in memory
    da = nc.Dataset(
        # f"{WORK_DIR}test/restartout.nc",
        f"{WORK_DIR}restartout.nc",
        "r+",
        clobber=True,
    )
    da["SoilMoist"][:] = sm_tmp
    da.close()

    logging.info("Job complete!")
