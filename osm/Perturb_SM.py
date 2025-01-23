#!/usr/bin/env python3
#   Python code to perturb the initial SM values for the finite difference approximation of the SEKF EDA
#   Jacobians
#  Author: original: D. Fairbairn (01/02/2022)
#
import os
import sys, getopt
from netCDF4 import Dataset

if __name__ == "__main__":

    nargs = len(sys.argv)

    scriptname = sys.argv[0].split("/")[-1]
    helptext = (
        scriptname
        + """
     -f <filename>           REQUIRED. eg: restartin.nc                                                                                                                                                                
     -v <Control variable>    REQUIRED either SM layer 1,2,3 window                                                                                                                                                  
    """
    )

    try:
      opts, args = getopt.getopt(
      sys.argv[1:],
      "f:v:i:",
      ["ofilename=", "ovariable=", "oinidate="],
      )
    except getopt.GetoptError:
      print(helptext)
      sys.exit(2)

    for opt, arg in opts:
        if opt in ("-h", "-help", "--help"):
            print(helptext)
            sys.exit()
        elif opt in ("-f", "--ofilename"):
            filename = arg
        elif opt in ("-v", "--ovariable"):
            variable = str(arg)
        elif opt in ("-i", "--oinidate"):
            inidate = str(arg)


    ncobj = Dataset(filename, 'r+')
    field = ncobj.variables['SoilMoist']

    if (inidate==".TRUE."):
        field[:,:] += 0.01
    else:
        if (variable=="1"):
          field[0,:] += 0.01*70.0    
        elif (variable=="2"):
          field[1,:] += 0.01*210.0    
        elif (variable=="3"):
          field[2,:] += 0.01*720.0    
        ncobj.close()


