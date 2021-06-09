from amuse.units import units, constants
import numpy as np 

class ExternalStellarEnergy:
    def __init__(self, star,planet, planet_lag_time=0.66 | units.s, planet_apsidal_motion_constant = 0.25):
        self.star = star
        self.planet = planet
        self.planet_lag_time = planet_lag_time
        self.planet_apsidal_motion_constant = planet_apsidal_motion_constant
    
    def n(self,separation): #mean motion
        #if separation is None:
        #    separation = self.separation
            #print('q')

        #check_n.append(np.sqrt(constants.G * (self.star.mass+self.planet.mass)/(separation**3)).value_in(units.s**(-1)))
        return np.sqrt(constants.G * (self.star.mass+self.planet.mass)/(separation**3))
 
    
    def period(self,separation):
        #if separation is None:בכל
        #    separation = self.separation
        #    print('qp')

        return 2.0*np.pi/self.n(self,separation)


    def tides_energy_on_planet(self,time_step, separation, eccentricity, tide_model='weak'):
       # weak_tide_timescale = self.planet.radius**3/(constants.G*self.planet.mass*self.planet_lag_time)
        if tide_model == 'weak':  # Hut 1981 equation A10
            return (time_step/self.period(separation))*Hut81Tides(separation, eccentricity, self.planet.radius, self.planet.mass, self.star.mass, self.planet_lag_time,
                              self.planet_apsidal_motion_constant)
        elif tide_model =='dynamical':
            if eccentricity > 0.8:
                return (time_step/self.period(separation))*MoeKratter18Tides(separation, eccentricity, self.planet.radius, self.planet.mass, self.star.mass)
            else:
                return (time_step/self.period(separation))*Hut81Tides(separation, eccentricity, self.planet.radius, self.planet.mass, self.star.mass,
                                  self.planet_lag_time,
                                  self.planet_apsidal_motion_constant)

        return 0.0 | units.erg / units.s

def Hut81Tides(separation, eccentricity,radius_planet, m_planet, M_star, tau, k):
#TODO: match between the input parameters and the 'real ones', the name of the function
    #this is A12 in Hut81
    return (9.0/2.0) * (constants.G**2) * (m_planet + M_star) * (m_planet**2) * \
           (radius_planet**5) * k * tau * (separation**(-9.0)) * (1-eccentricity**2)**(-15.0/2.0) * \
           (eccentricity**2) *\
            (1.0 + (15.0/4.0)*(eccentricity**2) + (15.0/8.0)*(eccentricity**4) + (5.0/64.0)*(eccentricity**6))

def MoeKratter18Tides(separation, eccentricity,radius_planet, m_planet, M_star):
    f_dyn = 0.1  # Pre-MS primaries, below eq. 9 MoeKratter2018

    return f_dyn * (constants.G**1.5) * (M_star**2.5) * (separation**7.5) * ((1 - eccentricity)**9) * \
           (M_star + m_planet) / (
                2 * constants.pi * m_planet * (radius_planet**10))

#f_dyn*((M_star+m_planet)/m_planet)*(constants.G*M_star**2/radius_planet)*\
#        (separation*(1-eccentricity)/radius_planet)**(-9)
#f_dyn * (constants.G**1.5) * (M_star**2.5) * (separation**7.5) * ((1 - eccentricity)**9) * \
 #          (M_star + m_planet) / (
 #               2 * constants.pi * m_planet * (radius_planet**10))
