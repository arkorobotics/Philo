#!/usr/bin/python

# Philo's Free Body Simulation
# ----------------------------
# Author: Ara Kourchians
# ----------------------------

import matplotlib
import numpy as np
import math
#import matplotlib.pyplot as plt
import json
import sys

# Natural Constants
R = 8.314	# J/(K*mol) or ((Kg*m^2)/s^2)/(K*mol)
g = 9.81	# m/s^2
P_amb = 101325	# Pa
T_amb = 298.15	# K


# Rocket Classes 
class Tank:
	def __init__(self, volume, pressure, tank_mass, reg_min, reg_max, temp):
		self.volume = volume		# m^3
		self.pressure = pressure	# Pa
		self.tank_mass = tank_mass	# kg
		self.reg_min = reg_min
		self.reg_max = reg_max
		self.temp = temp 			# K
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
		self.k = cp/cv # specific heat ratio
		print("vehicle.fuel.k: %f" %self.k)

class Engine:
	def __init__(self, fuel, A_throat, A_exit, Isp, Pc_max, Tc):
		self.fuel = fuel
		self.Isp = Isp
		self.A_throat = A_throat
		self.A_exit = A_exit
		self.Pc_max = Pc_max
		self.Tc = Tc
		self.ve_avg = self.Isp * g

		self.mass_flow = 0.0
		self.Pc = 0.0

	def ideal_thrust_coefficient(self, Pc):
		"""
			Calculate the ideal thrust coefficient of the engine from the chamber pressure.
			Equation 3-30 of Sutton and Biblarz, assuming p2 = p3.
		"""

		return math.sqrt( (2*self.fuel.k)/(self.fuel.k-1) * pow(2/(self.fuel.k + 1), (self.fuel.k + 1)/(self.fuel.k-1)) * ( 1 - pow(P_amb/self.Pc, (self.fuel.k - 1)/self.fuel.k) ) )

	def thrust_coefficient_viscous_effect(self):
		"""
			Calculates the viscuous effect on the thrust coefficient.
			Equation 5 of Hashem.
			Note - Requires the Reynolds number of the throat...
		"""
		return 17.6 * math.exp((0.0032 * self.A_exit) / self.A_throat)

	def thrust_coefficient(self, Pc):
		"""
			Returns the adjusted thrust coefficient of the engine
		"""
		return self.ideal_thrust_coefficient(Pc) # just the ideal coefficient for now. Should be C_ideal - C_viscuous_effect

	def calc_mass_flow(self, F):
		"""
			Returns the mass flow for a given thrust setting. Assumes Pexit = Pambient
		"""
		return F/self.v_e(F)

	def calc_Pc(self, F):
		"""
			Returns the required chamber pressure for a given thrust setting and mass flow
		"""
		return F / (self.thrust_coefficient() * self.A_throat)

	def calc_F(self, Pc):
		"""
			Returns the thrust produced at a given chamber pressure
		"""
		return self.thrust_coefficient() * self.A_throat * Pc

	def v_e(self, F):
		"""
			Returns the nozzle exit velocity given the current chamber pressure and temperature (P_chamber and T_chamber).
			This ignores chamber velocity V_1 and assumes ideal thermodynamic expansion such that nozzle exit pressure == ambient pressure
		"""
		return math.sqrt( (2*self.fuel.k)/(self.fuel.k-1) * R * self.Tc * (1 - pow(P_amb/self.calc_Pc(F), (self.fuel.k-1)/self.fuel.k) ) )

class Vehicle:
	def __init__(self, avionics_mass, mech_mass, tank, engine, fuel):
		self.avionics_mass = avionics_mass
		self.mech_mass = mech_mass
		self.tank = tank
		self.engine = engine
		self.fuel = fuel
		self.fuel_mass = ((self.tank.pressure * self.tank.volume)/(R*self.tank.temp)) * self.fuel.molar_mass

		self.dry_mass = self.tank.tank_mass + self.avionics_mass + self.mech_mass
		self.wet_mass = self.dry_mass + self.fuel_mass
		self.veh_mass = self.wet_mass # Initial (wet) condition

		self.Fnull = self.wet_mass * g # Initial
		self.mass_flow = self.Fnull/(self.engine.Isp * g) # Initial

	def calc_Fnull(self):
		return self.wet_mass * g


#regulator = {'':
def run_sim(vehicle):
	# Simulation
	print "Philo Sim"
	print "-----------------------"

	# Sim Constants
	n = 1000000
	dt = 0.001

	# Flight Variables
	flight_time = 0			# sec


	# Flight Sim
	# ----------------------------------------
	print("Initial Vehicle Dry Mass (kg): %f" %vehicle.dry_mass)
	print("Initial Fuel Mass (kg): %f" %vehicle.fuel_mass)
	print("Initial Vehicle Wet Mass (kg): %f" %vehicle.wet_mass)
	print("Exhaust Mass Flow (kg/s): %f" %vehicle.mass_flow)

	# Using Rocket Equation
	delta_v = vehicle.engine.ve_avg * np.log(vehicle.wet_mass/vehicle.dry_mass)
	print("Delta V (m/s): %f" %delta_v)

	# Using Mass Flow
	flight_time = vehicle.fuel_mass/vehicle.mass_flow
	print("Flight Time (sec): %f" %flight_time)

	flight_time = 0

	# Using Numerical Approx, mass loss thrust compensation
	# 	This function compensates for the change in thrust due to
	#	mass loss during flight. More fuel you use, the less thrust 
	#	you need to null-out m*g.
	while (vehicle.fuel_mass > 0):
		vehicle.Fnull = vehicle.calc_Fnull()
		vehicle.engine.mass_flow = vehicle.engine.calc_mass_flow(vehicle.Fnull)
		vehicle.fuel_mass -= vehicle.mass_flow*dt
		flight_time += dt 
	#print veh_mass
		
	print("\nEND OF FLIGHT")
	print("Flight Time (sec): %f" %flight_time)
	# ----------------------------------------

	#x = np.zeros ( (n,3) )

	# Plot Data
	#plt.figure()
	#plt.plot(x[:,0])
	#plt.show()

def load_config(cfgfile):
	cfg = None
	tanks = None
	fuels = None

	with open(cfgfile) as cfgson:
		cfg = json.load(cfgson)

	with open("tanks.json") as tankson:
		tanks = json.load(tankson)

	with open("fuels.json") as fuelson:
		fuels = json.load(fuelson)

	cfgdat = {
				"tank": None,
				"avionics_mass": None,
				"mech_mass": None,
				"fuel": None,
				"A_throat": None,
				"A_exit": None,
				"Isp": None,
				"Pc_max": None,
				"Tc": None
	}

	# get a list of the sections in the config
	sections = [section for section in cfg]
	
	# parse the config section by section and load values
	# into cfgdat
	if "Vehicle" in sections:
		sec = cfg["Vehicle"]

		if "tank" in sec:
			cfgdat["tank"] = tanks[sec["tank"]]
		if "avionics_mass" in sec:
			cfgdat["avionics_mass"] = sec["avionics_mass"]
		if "mech_mass" in sec:
			cfgdat["mech_mass"] = sec["mech_mass"]
	
	if "Engine" in sections:
		sec = cfg["Engine"]

		if "fuel" in sec:
			cfgdat["fuel"] = fuels[sec["fuel"]]
		if "A_throat" in sec:
			cfgdat["A_throat"] = sec["A_throat"]
		if "A_exit" in sec:
			cfgdat["A_exit"] = sec["A_exit"]
		if "Isp" in sec:
			cfgdat["Isp"] = sec["Isp"]
		if "Pc_max" in sec:
			cfgdat["Pc_max"] = sec["Pc_max"]
		if "Tc" in sec:
			cfgdat["Tc"] = sec["Tc"]
		

	# Initialize a Tank object
	tankdat = cfgdat["tank"]
	tankobj = Tank(tankdat["volume"], tankdat["pressure"], tankdat["tank_mass"], tankdat["reg_min"], tankdat["reg_max"], tankdat["temp"])
	
	# Initialize a Fuel Object
	fueldat = cfgdat["fuel"]
	fuelobj = Fuel(fueldat["fuel_type"], fueldat["cp"], fueldat["cv"], fueldat["molar_mass"])

	# Initialize an Engine Object
	engineobj = Engine(fuelobj, cfgdat["A_throat"], cfgdat["A_exit"], cfgdat["Isp"], cfgdat["Pc_max"], cfgdat["Tc"])

	# Create the vehicle object
	vehicleobj = Vehicle(cfgdat["avionics_mass"], cfgdat["mech_mass"], tankobj, engineobj, fuelobj)

	return vehicleobj
if __name__ == "__main__":
	if len(sys.argv) > 1:
		vehicle = load_config(sys.argv[1])

		run_sim(vehicle)



