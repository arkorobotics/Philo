# Load vehicle config

import json
import sys
import os.path

from .env import *
from .vehicle import *


def load_vehicle(cfgfile):
    # Default Config File
    philo_cfg_file = "philo_cfg.json"

    if len(cfgfile) > 1:
        if os.path.isfile(cfgfile[1]):
            philo_cfg_file = cfgfile[1]
            # vehicle = config.load_config(cfgfile[1])
        else:
            print("\n\033[91m{}\033[00m".format("ERROR: "))
            print("Config file not found.")

            print("\033[93m{}\033[00m".format("WARNING: "))
            print("Loading default config: %s" % philo_cfg_file)
    else:
        print("\n\033[93m{}\033[00m".format("WARNING:"))
        print(" Config file undefined... loading default config: %s" % philo_cfg_file)

    cfg = None
    tanks = None
    propellants = None
    regulators = None
    heaters = None
    engines = None

    with open(philo_cfg_file) as cfgson:
        cfg = json.load(cfgson)

    with open("vehicle/tanks.json") as tankson:
        tanks = json.load(tankson)

    with open("vehicle/propellants.json") as propellantson:
        propellants = json.load(propellantson)

    with open("vehicle/regulators.json") as regulatorson:
        regulators = json.load(regulatorson)

    with open("vehicle/heaters.json") as heaterson:
        heaters = json.load(heaterson)

    with open("vehicle/engines.json") as engineson:
        engines = json.load(engineson)

    cfgdat = {
        "avionics_mass": None,
        "mech_mass": None,
        "tank": None,
        "propellant": None,
        "regulator": None,
        "heaters": None,
        "engine": None
    }

    # get a list of the sections in the config
    sections = [section for section in cfg]

    # parse the config section by section and load values
    # into cfgdat
    if "Vehicle" in sections:
        sec = cfg["Vehicle"]

        if "avionics_mass" in sec:
            cfgdat["avionics_mass"] = sec["avionics_mass"]
        if "mech_mass" in sec:
            cfgdat["mech_mass"] = sec["mech_mass"]
        if "tank" in sec:
            cfgdat["tank"] = tanks[sec["tank"]]
        if "propellant" in sec:
            cfgdat["propellant"] = propellants[sec["propellant"]]
        if "regulator" in sec:
            cfgdat["regulator"] = regulators[sec["regulator"]]
        if "heater" in sec:
            cfgdat["heater"] = heaters[sec["heater"]]
        if "engine" in sec:
            cfgdat["engine"] = engines[sec["engine"]]

    if "Environment" in sections:
        sec = cfg["Environment"]

        if "P_ambient" in sec:
            cfgdat["P_ambient"] = sec["P_ambient"]

    # Initialize a Tank object
    tankdat = cfgdat["tank"]
    tankobj = Tank(tankdat["volume"], tankdat["pressure"], tankdat["tank_mass"], tankdat["reg_min"], tankdat["reg_max"],
                   tankdat["reg_out"], tankdat["temp"])

    # Initialize a Propellant Object
    propellantdat = cfgdat["propellant"]
    propellantobj = Propellant(propellantdat["propellant_type"], propellantdat["cp"], propellantdat["cv"],
                               propellantdat["molar_mass"])

    # Initialize a Regulator Object
    regulatordat = cfgdat["regulator"]
    regulatorobj = Regulator(regulatordat["in_min"], regulatordat["in_max"], regulatordat["out_min"],
                             regulatordat["out_max"], regulatordat["reg_in"], regulatordat["reg_out"])

    # Initialize a Heater Object
    heaterdat = cfgdat["heater"]
    heaterobj = Heater(heaterdat["T_chamber"])

    # Initialize an Engine Object
    enginedat = cfgdat["engine"]
    engineobj = Engine(enginedat["Ae"], enginedat["At"], propellantobj, regulatordat["reg_out"], cfgdat["P_ambient"], heaterobj)

    # Create the vehicle object
    vehicleobj = Vehicle(cfgdat["avionics_mass"], cfgdat["mech_mass"], tankobj, propellantobj, regulatorobj, heaterobj,
                         engineobj)

    return vehicleobj
