import matplotlib
import numpy as np
# import matplotlib.pyplot as plt
from .env import *


# Rocket Classes
class Propellant:
    def __init__(self, propellant_type, cp, cv, molar_mass):
        self.propellant_type = propellant_type
        self.cp = cp  # @ 0.0C
        self.cv = cv  # @ 0.0C
        self.molar_mass = molar_mass
        self.gamma = cp / cv


class Tank:
    def __init__(self, volume, pressure, tank_mass, reg_min, reg_max, reg_out, temp):
        self.volume = volume  # m^3
        self.pressure = pressure  # Pa
        self.tank_mass = tank_mass  # kg
        self.reg_min = reg_min
        self.reg_max = reg_max
        self.reg_out = reg_out
        self.temp = temp  # K


class Regulator:
    def __init__(self, in_min, in_max, out_min, out_max, reg_in, reg_out):
        self.in_min = in_min
        self.in_max = in_max
        self.out_min = out_min
        self.out_max = out_max
        self.reg_in = reg_in
        self.reg_out = reg_out


class Heater:
    def __init__(self, T_chamber):
        self.T_chamber = T_chamber


class Engine:
    def __init__(self, Ae, At, propellant, P_o, P_b, heater):
        self.Ae = Ae
        self.At = At
        self.propellant = propellant
        self.heater = heater
        self.P_b = P_b
        self.P_o = P_o #self.calc_Po()
        self.V_e = self.calc_Ve()
        self.Isp = self.V_e / g
        self.mass_flow = self.calc_mass_flow()
        self.mass_flow_max = self.calc_mass_flow_max()


    def calc_Po(self, margin=1.0):
        P_o = self.P_b/((2/(self.propellant.gamma+1))**(self.propellant.gamma/(self.propellant.gamma-1)))
        P_o *= margin
        
        return P_o

    def calc_Ve(self):
        """
            Updates V_e and Isp based on current pressure setting P_o.
            Returns V_e
        """
        V_e = np.sqrt(((self.heater.T_chamber * R) / self.propellant.molar_mass) \
                           * ((2 * self.propellant.gamma) / (self.propellant.gamma - 1)) \
                           * (1 - (self.P_b / self.P_o) ** (
                                   (self.propellant.gamma - 1) / self.propellant.gamma)))

        self.Isp = V_e / g
        return V_e

    def calc_mass_flow(self):
        mass_flow = (self.P_b*self.propellant.gamma*self.At) \
                         * np.sqrt((1/(self.propellant.gamma*R*T_amb)) \
                         *((2/(self.propellant.gamma+1))**((self.propellant.gamma+1)/(self.propellant.gamma-1))))

        return mass_flow

    def calc_mass_flow_max(self):
        # Calculated at nozzle
        self.Ae_optimal = 6.35e-7
        self.Ma = 1
        mass_flow_max = 0

        mass_flow_max = (self.Ae_optimal*self.Ma*self.P_o \
                            *np.sqrt(self.propellant.gamma/(R*self.heater.T_chamber))) \
                            / ((1+((self.propellant.gamma-1)*(self.Ma**2))/2)** \
                              ((self.propellant.gamma+1)/(2*(self.propellant.gamma-1))))

        return mass_flow_max
        
    '''
    TODO: Make functions later
    Ma1 = 1
    P_1 = P_b
    P_0 = P_1*(((1+((gamma-1)/2)*(Ma_t**2))**(gamma/(gamma-1))))
    print ("P_o: %0.2f Pa" % P_0)

    T_1 = T_b
    T_0 = T_1*(1+(((gamma-1)/2)*(Ma_t)))
    print ("T_o: %0.2f K" % T_0)

    rho_1 = 1
    rho_0 = rho_1*(((1+((gamma-1)/2)*(Ma_t**2))**(1/(gamma-1))))
    '''

class Vehicle:
    def __init__(self, avionics_mass, mech_mass, tank, propellant, regulator, heater, engine):
        self.avionics_mass = avionics_mass
        self.mech_mass = mech_mass
        self.tank = tank
        self.engine = engine
        self.propellant = propellant
        self.propellant_mass = ((self.tank.pressure * self.tank.volume) / (R * self.tank.temp)) \
                               * self.propellant.molar_mass

        self.dry_mass = self.tank.tank_mass + self.avionics_mass + self.mech_mass
        self.wet_mass = self.dry_mass + self.propellant_mass
        self.veh_mass = self.wet_mass  # Initial (wet) condition

        # self.Fnull_max = self.wet_mass * g
        # self.mass_flow_max = self.Fnull_max/self.engine.V_e_max

        self.Fnull = self.wet_mass * g
        self.mass_flow = self.Fnull / self.engine.V_e
