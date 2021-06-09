# EccentricPlanets

This is a script to evolve planets affected by external sources, both in their dynamics and internal evolution.

In order to use it, you need to install the following:
- Python 3.6 and above
- AMUSE 13.2 and above - <href="https://amusecode.github.io/"> <b> including MESA </b>
- numpy
- scipy
- shutil


In order to run your configuration run the following in corresponding to the location of amuse.sh and EvolvePlanet.sh:

<i> ./amuse.sh EvolvePlanet.py  --output_path=<OUTPUT>  --evolution_time=<TIME> --time_step=<TSTEP> --initial_separation=<SEPARATION> --initial_planet_mass=<MASS> --initial_planet_radius=<RADII> --initial_eccentricity=<ECCENTRICITY> --eccentricity_cutoff=<LOW_E> --star_luminosity=<STELLAR_L> --saved_planet=<MODEL> --initial_model=<MOD_FILE> --extension_opacity_file_path=<KAPPA_LOW> --tidal_efficiency=<T_Eff> --irradiation_efficiency=<L_Eff> --tides_model=<T_model> </i>

or run <i>./amuse.sh EvolvePlanet.py --help </i> to see relevant information.


<b><Large> Use your own code for dynamical evolution </Large></b>

You can make your own class of orbital evolution that inherites from OrbitalEvolution.py,
you should read to OrbitalEvolution.__init__(self, planet, star, initial_separation, initial_eccentricity,
                                  eccentricity_cuttoff=eccentricity_cuttoff) in the class __init__, and implement an evolve method that evolves self.separation and self.eccentricity


<b><Large> Use pre built MESA model for the planet </Large></b>

simply add <i> --initial_model=<MOD_FILE> </i> to your running script


<b><Large> Use MESA model created by this code </Large></b>

add <i> --saved_planet=<MODEL_FILE> </i> to your running script
  
  


  For any additional help please contact Hila Glanz- glanz@tx.technion.ac.il
