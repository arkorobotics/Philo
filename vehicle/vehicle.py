import numpy as np
from .env import *


# Rocket Classes
class Propellant(object):
    def __init__(self, propellant_type, cp, cv, molar_mass):
        self.propellant_type = propellant_type
        self.cp = cp  # @ 0.0C
        self.cv = cv  # @ 0.0C
        self.molar_mass = molar_mass
        self.gamma = cp / cv
        # TODO: ADD ENTHALPY


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
    def __init__(self, chamber_temp):
        self.chamber_temp = chamber_temp


class Engine:
    def __init__(self, a_e, a_t, propellant, p_o, p_b, heater):
        self.a_e = a_e
        self.a_t = a_t
        self.propellant = propellant
        self.heater = heater
        self.p_b = p_b
        self.p_o = p_o
        self.v_e = self.calc_v_e()
        self.isp = self.v_e / g
        self.mass_flow = self.calc_mass_flow()
        self.mass_flow_max = self.calc_mass_flow_max()

    @staticmethod
    def calc_isp(v_e):
        isp = v_e / g
        return isp  # sec

    @staticmethod
    def calc_area_ratio(ma, gamma):
        a_ratio = (1 / ma) * \
                  (((1 + ((gamma - 1) / 2) * (ma ** 2)) / (1 + ((gamma - 1) / 2))) ** ((gamma + 1) / (2 * (gamma - 1))))
        return a_ratio

    @staticmethod
    def calc_a_e(p_o, mass_flow, molar_mass, t_o, ma_e,  gamma):
        a_e = (mass_flow * (1 + (gamma - 1) * ((ma_e ** 2) / 2)) ** ((gamma + 1) / (2 * (gamma - 1)))) / (
               ma_e * p_o * np.sqrt((molar_mass * gamma) / (R * t_o)))
        return a_e  # m^2

    @staticmethod
    def calc_a_t(a_e, ma_e, gamma):
        a_t = (a_e * ma_e) * (
            ((2 / (gamma + 1)) * (1 + ((gamma - 1) / 2) * (ma_e ** 2))) ** ((-(gamma + 1)) / (2 * (gamma - 1))))
        return a_t

    def calc_mass_flow(self):
        mass_flow = (self.p_b*self.propellant.gamma*self.a_t) \
                    * np.sqrt((1/(self.propellant.gamma*(R/self.propellant.molar_mass)*ambient_temp)) \
                    * ((2/(self.propellant.gamma+1))**((self.propellant.gamma+1)/(self.propellant.gamma-1))))
        return mass_flow


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

        self.thrust = self.wet_mass * g
        self.mass_flow = self.Fnull / self.engine.V_e

    def update_veh_mass(self):
        self.propellant_mass = ((self.tank.pressure * self.tank.volume) / (R * self.tank.temp)) \
                               * self.propellant.molar_mass
        self.veh_mass = self.dry_mass + self.propellant_mass

    def calc_f_null(self):
        f_null = self.veh_mass * g
        return f_null                        # N

    def calc_mass_flow(self, v_e):
        mass_flow = self.thrust / v_e
        return mass_flow                    # kg/s

    @staticmethod
    def calc_delta_v(self, v_e, m_o, m_f):
        delta_v = v_e * np.log(m_o / m_f)
        return delta_v
