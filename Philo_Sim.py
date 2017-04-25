# Philo's Free Body Simulation
# ----------------------------
# Author: Ara Kourchians
# ----------------------------

import matplotlib
import numpy as np
import matplotlib.pyplot as plt

# Natural Constants
R = 8.314	# J/(K*mol) or ((Kg*m^2)/s^2)/(K*mol)
g = 9.81	# m/s^2
P_amb = 101325	# Pa
T_amb = 298.15	# K


# Rocket Classes 
class Tank:
	def __init__(self, volume, pressure, tank_mass, reg_min, reg_max, temp):
		self.volume = volume
		self.pressure = pressure
		self.tank_mass = tank_mass
		self.reg_min = reg_min
		self.reg_max = reg_max
		self.temp = temp
class Regulator:
	def __init__(self, in_min, in_max, out_min, out_max, reg_set, reg_in, reg_out):
		self.in_min = in_min
		self.in_max = in_max
		self.out_min = out_min
		self.out_max = out_max
		self.reg_set = reg_set
		self.reg_in = reg_in
		self.reg_out = reg_out

class Fuel:
	def __init__(self, fuel_type, cp, cv, molar_mass):
		self.fuel_type = fuel_type
		self.cp = cp				#@ 0.0C
		self.cv = cv				#@ 0.0C
		self.molar_mass = molar_mass
		self.specific_heat_ratio = cp/cv
		self.fuel_mass = 0.0

class Engine:
	def __init__(self,Isp):
		self.Isp = Isp
		
class Vehicle:
	def __init__(self, avionics_mass, mech_mass, tank, engine, fuel):
		self.avionics_mass = avionics_mass
		self.mech_mass = mech_mass
		self.tank = tank
		self.engine = engine
		self.fuel = fuel

# Rocket Dictionary
tanks = {'NIN-90-4500SO-GG':Tank(0.00147484, 3.1026e+7, 1.49, 3.103e+6, 5.516e+6, 298.15), \
	 'NIN-NINJAGREY68D':Tank(0.00111432, 3.1026e+7, 1.36, 3.103e+6, 5.516e+6, 298.15), \
	 'NIN-NINJAGREY50D':Tank(0.00081935, 3.1026e+7, 1.04, 3.103e+6, 5.516e+6, 298.15)}

fuel = {'Air':Fuel("Compressed Dry Air", 1.005, 0.718, 0.02897), \
	'CO2':Fuel("Carbon Dioxide", 0.846, 0.657, 0.04401)}

#regulator = {'':

# Simulation
print "Philo Sim"
print "-----------------------"

# Sim Constants
n = 1000000
dt = 0.001

# Tank Variables
tank_vol = 0.00147484		# m^3
tank_press = 3.1026e+7		# Pa
tank_temp = 298.15		# K 
tank_mass = 1.49		# kg

# Fuel Variables
fuel_type = "Compressed Dry Air"
cp = 1.005
cv = 0.718
fuel_molar_mass = 0.02897   	# kg/mol
specific_heat_ratio = cp/cv
fuel_mass = ((tank_press*tank_vol)/(R*tank_temp))*fuel_molar_mass

# Vehicle Variables
avionics_mass = 0.2		# kg
mech_mass = 0.2			# kg
wet_mass = tank_mass + fuel_mass + avionics_mass + mech_mass
dry_mass = tank_mass + avionics_mass + mech_mass
veh_mass = wet_mass

# Engine Variables
Isp = 100			# sec (TO BE DETERMINED)
V_e = Isp*g			# m/s
F_null = veh_mass*g		# N (kg*m/s^2)
mass_flow = F_null/V_e		# kg/s

# Flight Variables
flight_time = 0			# sec


# Flight Sim
# ----------------------------------------
print "Initial Fuel Mass (kg): \t\t\t %.8f" % fuel_mass
print "Initial Vehicle Wet Mass (kg): \t\t\t %.8f" % wet_mass
print "Exhaust Mass Flow (kg/s): \t\t\t %.8f" % mass_flow

# Using Rocket Equation
#delta_v = V_e*np.log(wet_mass/dry_mass)
#print "Delta V (m/s): \t\t\t %.8f" % delta_v

# Using Mass Flow
flight_time = fuel_mass/mass_flow
print "Fixed Mass Flow - Flight Time (sec): \t\t %.8f" % flight_time
flight_time = 0.0

# Using Numerical Approx, mass loss thrust compensation
# 	This function compensates for the change in thrust due to
#	mass loss during flight. More fuel you use, the less thrust 
#	you need to null-out m*g.
while (fuel_mass > 0):
	veh_mass = tank_mass + fuel_mass + avionics_mass + mech_mass
	F_null = veh_mass*g
	mass_flow = F_null/V_e
	fuel_mass -= mass_flow*dt
	flight_time += dt 
	#print veh_mass
	
print "Fixed Position - Flight Time (sec): \t\t %.8f" % flight_time
# ----------------------------------------

#x = np.zeros ( (n,3) )

# Plot Data
#plt.figure()
#plt.plot(x[:,0])
#plt.show()

