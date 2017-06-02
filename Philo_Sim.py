#!/usr/bin/python

# Philo's Free Body Simulation
# ----------------------------
# Author: Ara Kourchians
# ----------------------------

# TODO List:
# TODO: UPDATE Ma, Ve, and Area calculations with correct equations!

import matplotlib
import numpy as np
# import matplotlib.pyplot as plt
import sys
import os.path

from vehicle import *


def run_sim(vehicle):
    # Simulation
    print("\nPhilo Sim")
    print("-----------------------")

    # Sim Constants
    n = 1000000
    dt = 0.001

    # Flight Variables
    flight_time = 0  # sec

    # Flight Sim
    # ----------------------------------------
    print("Initial Vehicle Dry Mass (kg): \t\t\t %.6f" % vehicle.dry_mass)
    print("Initial Propellant Mass (kg): \t\t\t %.6f" % vehicle.propellant_mass)
    print("Initial Vehicle Wet Mass (kg): \t\t\t %.6f" % vehicle.wet_mass)
    print("Exhaust Mass Flow (kg/s): \t\t\t %.6f" % vehicle.mass_flow)

    # Using Rocket Equation
    delta_v = vehicle.engine.V_e * np.log(vehicle.wet_mass / vehicle.dry_mass)
    print("Delta V (m/s): \t\t\t\t\t %.6f" % delta_v)

    # Print Engine Isp
    print("Engine Specifc Impulse (sec): \t\t\t %.6f" % vehicle.engine.Isp)

    # Using Mass Flow
    flight_time = vehicle.propellant_mass / vehicle.mass_flow
    print("Constant Thrust - Flight Time (sec): \t\t %.6f" % flight_time)

    flight_time = 0

    # Using Numerical Approx, mass loss thrust compensation
    # 	This function compensates for the change in thrust due to
    #	mass loss during flight. More propellant you use, the less thrust
    #	you need to null-out m*g.
    while (vehicle.propellant_mass > 0):
        vehicle.veh_mass = vehicle.dry_mass + vehicle.propellant_mass
        vehicle.Fnull = vehicle.veh_mass * env.g
        vehicle.mass_flow = vehicle.Fnull / vehicle.engine.calc_Ve()
        vehicle.propellant_mass -= vehicle.mass_flow * dt
        #print ("%s" % vehicle.Fnull)
        flight_time += dt

    print("Constant Acceleration - Flight Time (sec): \t %.6f" % flight_time)
    print("-----------------------\n")


# ----------------------------------------

# x = np.zeros ( (n,3) )

# Plot Data
# plt.figure()
# plt.plot(x[:,0])
# plt.show()

if __name__ == "__main__":
    vehicle = config.load_vehicle(sys.argv)
    run_sim(vehicle)
