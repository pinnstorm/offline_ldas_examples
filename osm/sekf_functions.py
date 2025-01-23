#   Python code of the SEKF functions
#
#  Author: original: D. Fairbairn (31/01/2020)
#
#  1. Function to compute observation operator Jacobian matrix
#  2. SEKF increment function
#  3. Summary of USED observations stats

import os
import sys
import numpy as np
from sys import argv
from tabulate import tabulate

# 1. Observation operator Jacobians
#    --------------------------------------------------------
def EDA_Jacobian_calc(
    CovT2m_swvl1,
    CovT2m_swvl2,
    CovT2m_swvl3,
    CovRH2m_swvl1,
    CovRH2m_swvl2,
    CovRH2m_swvl3,
    CovSSM_swvl1,
    CovSSM_swvl2,
    CovSSM_swvl3,
    Var_swvl1,
    Var_swvl2,
    Var_swvl3,
    EDAH_TAPER,
):

    # Observation operator Jacobians
    Hjac = np.array(
        [
            [
                CovT2m_swvl1 / Var_swvl1,
                (1.0 / (1.0 + EDAH_TAPER)) * CovT2m_swvl2 / Var_swvl2,
                (1.0 / (1.0 + (2.0 * EDAH_TAPER))) * CovT2m_swvl3 / Var_swvl3,
            ],
            [
                CovRH2m_swvl1 / Var_swvl1,
                (1.0 / (1.0 + EDAH_TAPER)) * CovRH2m_swvl2 / Var_swvl2,
                (1.0 / (1.0 + (2.0 * EDAH_TAPER))) * CovRH2m_swvl3 / Var_swvl3,
            ],
            [
                CovSSM_swvl1 / Var_swvl1,
                (1.0 / (1.0 + EDAH_TAPER)) * CovSSM_swvl2 / Var_swvl2,
                (1.0 / (1.0 + (2.0 * EDAH_TAPER))) * CovSSM_swvl3 / Var_swvl3,
            ],
        ],
        dtype="float32",
    )

    return Hjac

    # Increment calculation


# 2. SEKF increment calculation
# -----------------------------------------------------------------------------------------------------------------
def Inc_calc_ecmwf(points, B_mat, H_mat, H_matT, R_mat, innov):

    invB = np.linalg.inv(B_mat)
    invR = np.linalg.inv(R_mat)

    #       R^{-1}H
    invR_H = np.einsum("ab,bdi->iad", invR, H_mat)
    #       H^TR^{-1}H
    Ht_invR_H = np.einsum("iab,ibd->iad", H_matT, invR_H)
    #       B^{-1}+H^TR^{-1}H
    invB_Ht_invR_H = Ht_invR_H + invB

    #       (B^{-1}+H^TR^{-1}H)^{-1}
    inv_invB_Ht_invR_H = np.linalg.inv(invB_Ht_invR_H)
    #       H^TR^{-1}
    Ht_invR = np.einsum("iab,bd->iad", H_matT, invR)
    #       K1=(B^{-1}+H^TR^{-1}H)^{-1}H^TR^{-1}
    K1 = np.einsum("iab,ibd->iad", inv_invB_Ht_invR_H, Ht_invR)
    #       dx=K1(yo-hx)
    Zinc = np.einsum("iab,bi->ia", K1, innov)

    return Zinc, K1


# 2. Alternative SEKF increment calculation
# -----------------------------------------------------------------------------------------------------------------
def Inc_calc(points, B_mat, Hmat, H_matT, Rmat, innov):

    BHtmat = np.einsum("ab,ibd->iad", B_mat, H_matT)

    HBHtmat = np.einsum("abi,ibd->iad", Hmat, BHtmat)

    print(Hmat, np.shape(Hmat), BHtmat, np.shape(BHtmat))

    HBHT_R = np.array(HBHtmat + np.swapaxes(Rmat, 0, 1))  # checked

    print("HBHT_R", np.shape(HBHT_R))

    inv_HBHT_R = np.linalg.inv(HBHT_R)
    print("inv", np.shape(inv_HBHT_R))

    K1 = np.einsum("iab,ibd->iad", BHtmat, inv_HBHT_R)

    print("K1", np.shape(K1))

    print("innov", np.shape(innov))

    Zinc = np.einsum("iab,bi->ia", K1, innov)

    return Zinc, K1


# 3. Friedmann bias calculation
# -----------------------------------------------------------------------------------------------------------------
def Inc_calc_bias(points, B_mat, Hmat, H_matT, Rmat, innov, gamma):

    BHtmat = np.einsum("ab,ibd->iad", B_mat, H_matT)

    HBHtmat = np.einsum("abi,ibd->iad", Hmat, BHtmat)

    HBHT_R = np.array(HBHtmat + (1.0 - gamma) * np.swapaxes(Rmat, 0, 1))  # checked

    inv_HBHT_R = np.linalg.inv(HBHT_R)

    K1 = gamma * np.einsum("iab,ibd->iad", BHtmat, inv_HBHT_R)

    Bias = np.einsum("iab,bi->ia", K1, innov)[:,0]

    return Bias


# 3. USED observations stastistics
# -----------------------------------------------------------------------------------------------------------------
def summarise(
    lat,
    lon,
    p_s,
    p_slvs,
    p_asc,
    innovs,
    Obs,
    obs_times,
    LUSE_ASCAT,
    SLV_assim,
    LASCAT_ADAPTIVE_BC, 
    innovs_bias,   
    model,
):

    headers = [
        "time",
        "obstype",
        "number",
        "min lat",
        "max lat",
        "min lon",
        "max lon",
        "min obsvalue",
        "max obsvalue",
    ]
    headers.extend(
        [
            "min state",
            "max state",
            "min depar",
            "mean depar",
            "max depar",
        ]
    )
    table = []

    #   print ASCAT statistics:
    for tt in obs_times:
        if LUSE_ASCAT == "true":
            line = [
                tt,
                "ASCAT fg_depar" + tt,
                len(p_asc[tt]),
                lat[p_s[p_asc[tt]]].min(),
                lat[p_s[p_asc[tt]]].max(),
                lon[p_s[p_asc[tt]]].min(),
                lon[p_s[p_asc[tt]]].max(),
                Obs["Asc" + tt][p_s[p_asc[tt]]].min(),
                Obs["Asc" + tt][p_s[p_asc[tt]]].max(),
            ]
            table.append(line)

            line.extend(
                [
                    (model["SM_background"][0, 0, p_s[p_asc[tt]]] / 70.0).min(),
                    (model["SM_background"][0, 0, p_s[p_asc[tt]]] / 70.0).max(),
                    innovs[tt][2, p_asc[tt]].min(),
                    innovs[tt][2, p_asc[tt]].mean(),
                    innovs[tt][2, p_asc[tt]].max(),
                ]
            )

            if (LASCAT_ADAPTIVE_BC=="true"):

              line = [
                tt,
                "ASCAT BC_depar" + tt,
                len(p_asc[tt]),
                lat[p_s[p_asc[tt]]].min(),
                lat[p_s[p_asc[tt]]].max(),
                lon[p_s[p_asc[tt]]].min(),
                lon[p_s[p_asc[tt]]].max(),
                Obs["Asc" + tt][p_s[p_asc[tt]]].min(),
                Obs["Asc" + tt][p_s[p_asc[tt]]].max(),
              ]
              table.append(line)

              line.extend(
                [
                    (model["SM_background"][0, 0, p_s[p_asc[tt]]] / 70.0).min(),
                    (model["SM_background"][0, 0, p_s[p_asc[tt]]] / 70.0).max(),
                    innovs_bias[tt][2, p_asc[tt]].min(),
                    innovs_bias[tt][2, p_asc[tt]].mean(),
                    innovs_bias[tt][2, p_asc[tt]].max(),
                ]
              )

            line = [
                tt,
                "ASCAT an_depar" + tt,
                len(p_asc[tt]),
                lat[p_s[p_asc[tt]]].min(),
                lat[p_s[p_asc[tt]]].max(),
                lon[p_s[p_asc[tt]]].min(),
                lon[p_s[p_asc[tt]]].max(),
                Obs["Asc" + tt][p_s[p_asc[tt]]].min(),
                Obs["Asc" + tt][p_s[p_asc[tt]]].max(),
            ]
            table.append(line)

            line.extend(
                [
                    (model["SM_control"][0, p_s[p_asc[tt]]] / 70.0).min(),
                    (model["SM_control"][0, p_s[p_asc[tt]]] / 70.0).max(),
                    innovs[tt+'_an_depar'][p_asc[tt]].min(),
                    innovs[tt+'_an_depar'][p_asc[tt]].mean(),
                    innovs[tt+'_an_depar'][p_asc[tt]].max(),
                ]
            )

    #   print T2m statistics:
    for tt in obs_times:
        if SLV_assim:
            line = [
                tt,
                "T2m" + tt,
                len(p_slvs[tt]),
                lat[p_s[p_slvs[tt]]].min(),
                lat[p_s[p_slvs[tt]]].max(),
                lon[p_s[p_slvs[tt]]].min(),
                lon[p_s[p_slvs[tt]]].max(),
                Obs["T2m" + tt][p_s[p_slvs[tt]]].min(),
                Obs["T2m" + tt][p_s[p_slvs[tt]]].max(),
            ]
            table.append(line)

            line.extend(
                [
                    model["T2m_background"][0, p_s[p_slvs[tt]]].min(),
                    model["T2m_background"][0, p_s[p_slvs[tt]]].max(),
                    innovs[tt][0, p_slvs[tt]].min(),
                    innovs[tt][0, p_slvs[tt]].mean(),
                    innovs[tt][0, p_slvs[tt]].max(),
                ]
            )

    #   print RH2m statistics:
    for tt in obs_times:
        if SLV_assim:
            line = [
                tt,
                "RH2m" + tt,
                len(p_slvs[tt]),
                lat[p_s[p_slvs[tt]]].min(),
                lat[p_s[p_slvs[tt]]].max(),
                lon[p_s[p_slvs[tt]]].min(),
                lon[p_s[p_slvs[tt]]].max(),
                Obs["RH2m" + tt][p_s[p_slvs[tt]]].min(),
                Obs["RH2m" + tt][p_s[p_slvs[tt]]].max(),
            ]
            table.append(line)

            line.extend(
                [
                    model["RH2m_background"][0, p_s[p_slvs[tt]]].min(),
                    model["RH2m_background"][0, p_s[p_slvs[tt]]].max(),
                    innovs[tt][1, p_slvs[tt]].min(),
                    innovs[tt][1, p_slvs[tt]].mean(),
                    innovs[tt][1, p_slvs[tt]].max(),
                ]
            )

    print(
        "######################### Summary of USED observations stats ###########################"
    )
    print(tabulate(table, headers=headers))

    return
