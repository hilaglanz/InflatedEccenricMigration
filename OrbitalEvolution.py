import numpy as np
from amuse.units import units,constants

separation_units = units.cm
time_units = units.s
energy_units = units.erg

check_rate1 = []
check_time = []
check_n = []

class OrbitalEvolution:
    def __init__(self, planet, star, initial_separation,initial_eccentricity, eccentricity_cuttoff=1e-5,
                 tide_model="weak"):
        self.planet = planet
        self.star = star 
        self.separation = initial_separation
        self.previous_separation = initial_separation
        self.eccentricity = initial_eccentricity
        self.eccentricity_cutoff = eccentricity_cuttoff
        self.tide_model = tide_model
        
    def n(self, separation=None):  # mean motion
        if separation is None:
            separation = self.separation

        return np.sqrt(constants.G * (self.star.mass + self.planet.mass) / (separation ** 3))

    def period(self, separation=None):
        if separation is None:
            separation = self.separation

        return 2.0 * np.pi / self.n(separation)

    def orbital_energy(self):

        return -constants.G*self.planet.mass*self.star.mass/(2*self.separation)

    def evolve(self, time):
        if self.eccentricity <= self.eccentricity_cutoff:
            print('we reached e=',self.eccentricity_cutoff,' couldn\'t evolve further')
            raise Exception("we reached e=0 or a=0 AU, couldn\'t evolve further")

        if self.separation * (1-self.eccentricity) <= self.Roche_radius():
            print('planet is ripped apart because the pericenter is smaller than its Roche radius')
            raise Exception("planet is ripped apart")

        raise Exception("not implemented")

    
    def tides_energy_on_planet(self):

        return 0.0 | units.erg / units.s
    
    def irradiation_flux(self):
        try:
            luminosity = self.star.luminosity
        except:
            print("couldnt find luminoisty of the star, assuming solar value")
            luminosity = 1.0 | units.LSun
        
        return luminosity * (self.planet.radius/self.separation)**2/(1-self.eccentricity**2)**0.5

    def Roche_radius(self,eta=2.7):#eta from Guillochon et al 2011
        return eta*self.planet.radius*(self.star.mass/self.planet.mass)**(1/3)
