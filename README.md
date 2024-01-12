We study the effect of periodic glucose injection on glycolytic oscillations by using a allosteric model of phosphofructokinase proposed by Goldbeter. (A. Goldbeter and R. Lefever, Biophys. J. 12, 1302 (1972))

- deterministic_figure.py (A Python script to generate a trajectory of glycolytic oscillations)

- deterministic_findT.py (A Python script to calculate the period of glycolytic oscillations)

- lyapunov_Benettin_realsys.py (A Python script to calculate the Lyapunov exponent averaged over the N x driving period of glycolytic oscillations)


All python scripts need two input parameters: (i) amplitude and (2) period of external driving. 


In order to generate a dynamic trajectory, say, at input amplitude of 0.0025 mM/s (=0.5 υ_0) and period of T_ext=100 s, and calculate the resulting oscillatory period (T_obs) and two Lyapunov exponents, execute the script by typing the followings: 

“python deterministic_figure.py 0.0025 100” --  generates a figure of trajectory,  “Trajectory_A0.0025_Tin100.png”

“python deteministic_findT.py 0.0025 100” -- calculates the mean period of oscillations (T_obs) and its standard deviation. If the standard deviation is large, it means that the trajectory is irregular, either quasi-periodic or chaotic.

“python lyapunov_Benettin_realsys.py 0.0025 100” -- calculates the two Lyapunov exponents except 0. Our system subject to an external periodic driving can have three Lyapunov exponents, but the Lyapunov exponent corresponding external period is always 0.

