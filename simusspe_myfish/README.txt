SPECIFIC TO SIMULATION SETTINGS AND SCENARIOS
=============================================

each scenario (called in the c++ command line)
requires a specific input .dat e.g. baseline.dat,
high_profit_grounds.dat, etc.
gathering the high level parameters of the simulation

config.dat
is used for fine tuning of some criticial parameters
for calibration purpose
This file also lists the index of the implicit_stocks
which are the stocks that will take several shortcuts
at various steps of the simulation because being too data-poor 


OTHERS
================================


FOr convenience and easy detection
for the discrete events occuring at various time scale
during the simulation, some data will just gather in STL vectors
the index of the relevant time steps (in hours spent):
 
tstep_years_2012_2017.dat
tstep_semesters_2012_2017.dat
tstep_quarters_2012_2017.dat
tstep_months_2012_2017.dat