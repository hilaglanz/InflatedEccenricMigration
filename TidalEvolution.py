from  OrbitalEvolution import *
import numpy as np
from amuse.units import units, constants
from scipy.integrate import odeint

separation_units = units.cm
time_units = units.s
energy_units = units.erg

class TidalEvolution(OrbitalEvolution):
    def __init__(self, planet, star, initial_separation, initial_eccentricity,
                     planet_lag_time=0.66 | units.s, planet_apsidal_motion_constant=0.25, eccentricity_cuttoff=1e-5,
                     tide_model="weak", dynamical_limit = 0.3):
        OrbitalEvolution.__init__(self, planet, star, initial_separation, initial_eccentricity,
                                  eccentricity_cuttoff=eccentricity_cuttoff)
        self.tide_model = tide_model
        self.planet_lag_time = planet_lag_time
        self.planet_apsidal_motion_constant = planet_apsidal_motion_constant
        self.dynamical_limit = dynamical_limit

    def f(self, eccentricity=None):  # Hamers & Tremaine 2017 eq. 18
        if eccentricity is None:
            eccentricity = self.eccentricity

        return (1 + (45.0 / 14.0) * (eccentricity ** 2) + 8.0 * (eccentricity ** 4) + (685.0 / 224.0) * (
                    eccentricity ** 6) +
                (255.0 / 448.0) * (eccentricity ** 8) + (25.0 / 1792.0) *
                (eccentricity ** 10)) / (1.0 + 3.0 * (eccentricity ** 2) + (3.0 / 8.0) * (eccentricity ** 4))
   
    def transition_dynamical_to_weak(self,f_dyn=0.01):
        return 2.*f_dyn*self.planet.radius**3*(1-self.eccentricity**2)**7.5/\
            (21*constants.G*self.planet.mass*self.planet_apsidal_motion_constant*self.planet_lag_time*self.period()*(1-self.eccentricity)**9*self.f()*self.eccentricity**2)
  
    def evolve(self, time, number_of_time_steps=1000):
        if self.separation * (1-self.eccentricity) <= self.Roche_radius():
            print('planet is ripped apart because the pericenter is smaller than its Roche radius')
            raise Exception("planet is ripped apart")

        if self.eccentricity <= self.eccentricity_cutoff:
            print('we reached e=', self.eccentricity_cutoff, ' couldn\'t evolve further')
            raise Exception("we reached e=0 or a=0 AU, couldn\'t evolve further")

        time_array = np.linspace(0, time.value_in(time_units), number_of_time_steps)
        initial_conditions = np.array([self.separation.value_in(separation_units), self.eccentricity])
        if self.tide_model == 'weak':  # Hut 1981 equation A10
            y = odeint(weak_tides_ODE, initial_conditions, time_array, args=(self,))
        elif self.tide_model == 'dynamical':
            if self.eccentricity > self.dynamical_limit:
                transition = self.transition_dynamical_to_weak()
                if transition>1.:
                    y = odeint(dynamical_tides_ODE, initial_conditions, time_array, args=(self,))
                    e = y[:, 1][-1]
                    if e < self.dynamical_limit * 0.6 or e >= 1.0:
                        print("transition during evolution with dynamical model, decreasing timesteps to", time/2.0)
                        old_separation = self.separation
                        self.evolve(time/2.0, number_of_time_steps)
                        print("middle e: ",self.eccentricity)
                        self.evolve(time/2.0, number_of_time_steps)
                        y = np.array([[self.separation.value_in(separation_units),self.eccentricity]])
                        self.separation = old_separation
                else:
                    y = odeint(weak_tides_ODE, initial_conditions, time_array, args=(self,))
            else:
                y = odeint(weak_tides_ODE, initial_conditions, time_array, args=(self,))
       
        else:
            print('This tide model is not allowed')
            raise Exception("tide model " + self.tide_model + " is not allowed")

        a = y[:, 0]
        e = y[:, 1]
        self.previous_separation = self.separation
        self.separation = a[-1] | separation_units
        self.eccentricity = e[-1]

        print("orbit evolved to separation: ", self.separation.as_quantity_in(units.AU), " eccentricity: ", self.eccentricity)
  
    def tides_energy_on_planet(self):
        if self.separation <= self.Roche_radius():
            print('planet is ripped apart because the separation is smaller than its Roche radius')
            raise Exception("planet is ripped apart")

        if self.eccentricity <= self.eccentricity_cutoff:
            print('we reached e=', self.eccentricity_cutoff, ' couldn\'t evolve further')
            raise Exception("we reached e=0 or a=0 AU, couldn\'t evolve further")

        if self.tide_model == 'weak':  # Hut 1981 equation A10
            tidal_heat = self.weak_tides_energy_on_planet()

        elif self.tide_model == 'dynamical':
            if self.eccentricity > self.dynamical_limit:
                transition = self.transition_dynamical_to_weak()
                if transition > 1.:
                    tidal_heat = MoeKratter18Tides(self.separation, self.eccentricity, self.planet.radius,
                                                   self.planet.mass,
                                                   self.star.mass)/self.period()
                else: 
                    tidal_heat = self.weak_tides_energy_on_planet()
            else:
                tidal_heat = self.weak_tides_energy_on_planet()
        else:
            print('This tide model is not allowed')
            raise Exception("tide model " + self.tide_model + " is not allowed")

        return tidal_heat

    def weak_tides_energy_on_planet(self):
        return (self.orbital_energy()/self.separation)*dadt_weak(separation=self.separation,
                                                                                   eccentricity=self.eccentricity,
                                                                                   orbital_evolution_model=self)

    def planet_bulge_from_tides(self):
        if self.tide_model == 'weak':
            return self.planet.radius - self.planet_bulge_height_from_weak_tides()
        else:
            return self.planet.radius

    def planet_bulge_height_from_weak_tides(self):
        return (self.star.mass/self.planet.mass)*self.planet.radius*(self.planet.radius/self.separation)**3


def weak_tides_ODE(y, time, orbital_evolution_model):  # Hamers & Tremaine 2017 eqns. 17
    separation, eccentricity = y
    separation = separation | separation_units
    dadt = dadt_weak(separation,eccentricity,orbital_evolution_model)
    dedt = dedt_weak(separation,eccentricity,orbital_evolution_model)

    return [dadt.value_in(separation_units / time_units), dedt.value_in(time_units ** (-1))]

def dadt_weak(separation,eccentricity,orbital_evolution_model):
    return -21.0 * orbital_evolution_model.planet_apsidal_motion_constant * (
        orbital_evolution_model.n(separation)) ** 2 * \
           orbital_evolution_model.planet_lag_time * (
                       orbital_evolution_model.star.mass / orbital_evolution_model.planet.mass) * \
           (orbital_evolution_model.planet.radius / separation) ** 5 * separation * eccentricity ** 2 * \
           orbital_evolution_model.f(eccentricity) / ((1 - eccentricity ** 2) ** (15.0 / 2.0))

def dedt_weak(separation,eccentricity,orbital_evolution_model):
    return -(21.0 / 2.0) * orbital_evolution_model.planet_apsidal_motion_constant * (
        orbital_evolution_model.n(separation)) ** 2 * \
           orbital_evolution_model.planet_lag_time * (
                       orbital_evolution_model.star.mass / orbital_evolution_model.planet.mass) * \
           (orbital_evolution_model.planet.radius / separation) ** 5 * eccentricity * \
           orbital_evolution_model.f(eccentricity) / ((1 - eccentricity ** 2) ** (13.0 / 2.0))


def dynamical_tides_ODE(y, time, orbital_evolution_model):  # MoeKratter2018
    separation, eccentricity = y
    separation = separation | separation_units
    E_orb = constants.G * orbital_evolution_model.planet.mass * \
        orbital_evolution_model.star.mass / (2 * separation)

    Delta_E = MoeKratter18Tides(separation=separation, eccentricity=eccentricity,
                                radius_planet=orbital_evolution_model.planet.radius,
                                m_planet=orbital_evolution_model.planet.mass,M_star=orbital_evolution_model.star.mass)
    dadt = -(separation / orbital_evolution_model.period(separation)) * Delta_E / E_orb
    dedt = ((1 - eccentricity) / separation) * dadt

    return [dadt.value_in(separation_units / time_units), dedt.value_in(time_units ** (-1))]


def Hut81Tides(separation, eccentricity, radius_planet, m_planet, M_star, tau, k):
    # this is A12 in Hut81
    return (9.0 / 2.0) * (constants.G ** 2) * (m_planet + M_star) * (M_star ** 2) * \
           (radius_planet ** 5) * k * tau * (separation ** (-9.0)) * (1 - eccentricity ** 2) ** (-15.0 / 2.0) * \
           (eccentricity ** 2) * \
           (1.0 + (15.0 / 4.0) * (eccentricity ** 2) + (15.0 / 8.0) * (eccentricity ** 4) + (5.0 / 64.0) * (
                       eccentricity ** 6))

def MoeKratter18Tides(separation, eccentricity, radius_planet, m_planet, M_star):
    f_dyn = 0.01  # Pre-MS primaries, below eq. 9 MoeKratter2018
    Delta_E = f_dyn * ((M_star + m_planet) / m_planet) * \
              (constants.G * M_star ** 2 / radius_planet) * \
              (separation * (1 - eccentricity) / radius_planet) ** (-9)

    return Delta_E
