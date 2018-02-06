# MahonFitting v1.0
<img src="https://raw.githubusercontent.com/trappitsch/MahonFitting/master/icon/fitico.png?token=AKKyQI8F3m2egXde2bbREh4UHeIYhwH7ks5affc6wA%3D%3D" width="128" height="128">

<!--toc-->

## Introduction

This routine performs linear regressions to a set of measurement that have uncertainties in both axes, which we subsequently call *x* and *y*. The regression is based on the approach by Mahon (1996) and is the same regression as used by York (1969), however, with the corrected calculation for the uncertainties. The code presented here accompanies the publication of Trappitsch et al. (2018). If you want, please cite Trappitsch et al. (2018) as a reference for pointing at this tool.

## Usage
### Preparing your data set
The fitting program can read text files and csv files. The general format of your file has to be:

Column 1 | Column 2                  | Column 3 | Column 4                 |
---------|---------------------------|----------|--------------------------|
x values | uncertainties of x values	| y values | uncertainties of y values |

For text files you can seperate the individual values with a space or a tab, for csv values with a comma. There has to be a line break to separate between different datapoints.

**Please to not use a header in your datafile!**

The testfiles that are provided are the files that come from the dataset published in Mahon (1996) and can be used for validation.

### Running a regression
If you want to calculate a regression, please follow this procedure. We assume you have prepared a file *data.txt* that is prepared as described above.

1. Click *Load txt* and select *data.txt*
2. Select if your uncertainties are 1&sigma; or 2&sigma;. The uncertainty in the output will be the same as for the input.
2. Press *Calculate*
3. Your results will be given below, these results can be directly copied. The program also writes an output file into the folder that you have selected.

If you would like to see your data along with the regression in a quick plot hit the *Plot* button. This gives you a rough overview. if you would like you can save the plot using the button at the bottom.

### Using a fixed intercept
If you want fix the y-axis intercept, click the check box and type in the intercept value that you want to fix. The software then artificially adds another datapoint to the regression with the *(x, y)* values *(0, f0)*, where *f0* is the typed in intercept. The uncertainties for the intercept are chosen such that they are 10<sup>12</sup> times lower than the smallest uncertainty in your dataset. 


## Pre-compiled packages
You can find full standalone programs that include the necessary python libraries in the released folder or by clicking [here](https://github.com/LLNL/MahonFitting/releases). All precompiled packages have been generated using the [pyinstaller](http://www.pyinstaller.org) package.

## Instalation from source
There is no real installation necessary, however, if you would like to run the source using your own python distriution, you can download the the the [mahon.py](https://github.com/LLNL/MahonFitting/blob/master/mahon.py) file and run it from your python distribution. The program was written in python2.7 and the required packages are:

* tkinter
* numpy
* matplotlib

Of course you can always run the script from a terminal by calling

    python mahon.py

where `python` is the shortcut for your python distribution.

### OSX
Here are step by step instructions for creating an exectuable for running your program. 

1. Put the files [mahon.py](https://github.com/LLNL/MahonFitting/blob/master/mahon.py) and [mahon.command](https://github.com/LLNL/MahonFitting/blob/master/mahon.command) into the same folder (which we will call `$APPFOLDER` here)
2. If you would like to use your own python distribution or your system's python distribution is not at `/usr/bin/python`, modify the first line of the command file to point to your python distribution (for most OSX systems, this step can be skipped)
3. Now we need to make the command file an executable:
 31. Open a terminal and navigate to the folder in which you put the script by typing  
 	  `cd $APPFOLDER`
 32. Make the *mahon.command* file executable by typing  
 	  `chmod 755 mahon.py`
4. Now you can run the command script by either double clicking it in the finder or by running  
`./mahon.command`  
from the terminal.
5. You can create a shortcut of the mahon.command in the Dock and download the [icon](https://github.com/LLNL/MahonFitting/blob/master/icon/fitico.icns) to associate it with the command file in the dock.

### Windows
Windows does not come natively with python. If you have your own python installation, chances are you already know how to run the script. If you don't have your own python and would like to still run from source (or if the executable does not work for some reason), please follow these instructions.

1. Install a python interpreter on your system. One example to use would be [Anaconda](https://www.anaconda.com/download/), which can be downloaded by clicking on the link. Make sure to download **python2.7**
2. Put the file [mahon.py](https://github.com/LLNL/MahonFitting/blob/master/mahon.py) into a a folder (we will name this folder `$APPFOLDER` here)
2. Anaconda comes with a python editor called *Spyder*. Open *mahon.py* in *Spyder* and click the run button, the program should start now.
3. If you would like to create an icon on the Desktop follow these simple steps:
 31. Create a shortcut of *Spyder* (which is actually an editor written in python) on your Desktop
 32. Rename this shortcut to *MahonFitter* or any other name you like
 33. Open the properties of this file
 34. In the entrybox that is named *Target:*, there are two programs called, first is the call to the executable python distribution, the second entry is the call to the spyder python script. Change this path such that is says:  
     `PATH_TO_PYTHON\pythonw.exe $APPFOLDER\mahon.py`
 35. Hit okay to apply the changes. Double clicking this icon should now start the program.
 36. If you want to use the [icon](https://github.com/LLNL/MahonFitting/blob/master/icon/fitico.ico), download it and point to the downloaded location after you clicked *Change Icon...* in the properties manager

## Contact
If you find a bug or problem with the software or if you can't get it to run on your machine (and have read this readme file), feel free to contact me (Reto Trappitsch) at <trappitsch1@llnl.gov>.

## References
Mahon K. I. (1996) The New "York" Regression: Application of an Improved Statistical Method to Geochemistry, *International Geology Review*, 38:293-303.

York D. (1969) Least squares fitting of a straight line with correlated errors, *Earth & Planetary Science Letters*, 5:320-324.

Trappitsch R., Boehnke P., Stephan T., Telus M., Savina M. R., Pardo O., Davis A. M., Dauphas N., Pellin M. J., adn Huss G. (2018) New Constraints for the Low Abundance of <sup>60</sup>Fe in the Early Solar System, *The Astrophysical Journal Letters*, in prep.

## Release
Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at the Lawrence Livermore National Laboratory. Written by Reto Trappitsch <trappitsch1@llnl.gov>LLNL-CODE-745740 All rights reserved. This file is part of MahonFitting v1.0Please also read this link – Our [Disclaimer](https://github.com/LLNL/MahonFitting/blob/master/DISCLAIMER) and [GNU General Public License](https://github.com/LLNL/MahonFitting/blob/master/LICENSE).This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License (as published by the Free Software Foundation) version 2, dated June 1991.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA