# AlphaTraj

AlphaTraj is software designed for analyzing pocket dynamics. 

## Features

Compared to other pocket dynamics analysis software, AlphaTraj has the following features:

1. Divides protein pockets into several sub-pockets (implemented through modified [AlphaSpace2.0](https://github.com/RedesignScience/AlphaSpace2/tree/master) code).
2. Analyzes various dynamic properties of sub-pockets, such as lifetime, fluctuation of pocket volume, average nonpolarizability of the pocket, etc.
3. Examines relationships between sub-pockets, such as coexistence rate and correlation.
4. Implements a scoring function tailored for sub-pockets to characterize their quality.
5. Develops a scoring function for the entire pocket to automatically find the optimal conformations suitable for inhibitor design.

## Usage

Please refer to the [tutorial](./tutorial/AlphaTrajTutorial_ENG.md).

## Citation

If AlphaTraj is used in your research, please cite the following paper in your manuscript:  
Lv, N.; Cao, Z. Subpocket-Based Analysis Approach for the Protein Pocket Dynamics. *J. Chem. Theory Comput.* 2024, 20 (11), 4909–4920.  
Literature link: https://pubs.acs.org/doi/10.1021/acs.jctc.4c00476

## Additional Information:
If you have any questions or suggestions, please feel free to contact me at lvnan@htu.edu.cn or lvnan@stu.xmu.edu.cn 


# <font color=red>**Notice!!!**</font>  
Due to an update in the `mdtraj` version, the `_sasa()` function interface has been modified (version 1.10.0 uses the old interface, while the current version 1.10.2 uses the new interface), which may cause runtime errors.  

The last few lines of the error message are as follows:  
```
File "mdtraj/geometry/src/_geometry.pyx", line 227, in mdtraj.geometry._geometry._sasa
TypeError: _sasa() takes exactly 6 positional arguments (1 given)
```

This issue has now been resolved, and AlphaTraj can run normally on both the new and old versions of `mdtraj`.  
Users experiencing this error are advised to download the latest version of AlphaTraj.