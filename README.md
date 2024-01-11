We study the effect of periodic glucose injection on glycolytic oscillations by using a allosteric model of phosphofructokinase proposed by Goldbeter. (A. Goldbeter and R. Lefever, Biophys. J. 12, 1302 (1972))

- deterministic_figure.py (A Python script to generate a trajectory of glycolytic oscillations)

- deterministic_findT.py (A Python script to calculate the period of glycolytic oscillations)

- lyapunov_Benettin_realsys.py (A Python script to calculate the Lyapunov exponent averaged over the N x driving period of glycolytic oscillations)

User Guide

All codes need two input parameters: Amplitude, Input Period
- deterministic_figure.py gives trajectory figure named "Trajectory_A..._Tin....png"
- deterministic_findT.py gives the mean period of oscillations and standard deviation of period. If standard deviation of period is large, the period of oscillations is irregular,i.e., quasi-periodic or chaotic.
- lyapunov_Benettin_realsys.py gives the Lyapunov exponent.
