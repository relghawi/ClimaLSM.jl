"""This file contains simulation variables for running CliMA LSM on the US-NR1
fluxtower site. This includes both the domain variables and timestepping 
variables for running the simulation."""

# DOMAIN SETUP:

# Column dimensions - separation of layers at the top and bottom of the column:
dz_bottom = FT(1.25)
dz_top = FT(0.05)

# Stem and leaf compartments and their heights:
n_stem = Int64(1)
n_leaf = Int64(1)
h_leaf = FT(6.5) # m
h_stem = FT(7.5) # m

# TIME STEPPING:

# Starting time:
t0 = Float64(120 * 3600 * 24)# start mid year to avoid snow

# Time step size:
dt = Float64(40)

# Number of timesteps between saving output:
n = 45
