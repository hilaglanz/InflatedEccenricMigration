import os, sys, time
import argparse
from OrbitalEvolution import *
from amuse.units import units
from amuse.datamodel import Particle
from InternalPlanetEvolution import *
from TidalEvolution import *

def main(initial_mod_file, opacity_extension ,saved_planet, initial_planet_mass, initial_planet_radius, star_mass,
         star_luminosity, initial_planet_metalicity, output_path,
         initial_separation, initial_eccentricity, eccentricity_cutoff, radius_cutoff, total_time,
         time_step, tide_model, tidal_efficiency, radiation_efficiency):
    #initializes the planet to be evolved
    plotting = "a,e,R,Tc,Teff,Pc,Psurf,Ltides,Lrad,t\n" #header for plotting
    #initializes planet with MESA
    planet_model = InternalPlanetEvolution(initial_mod_file ,initial_planet_mass, initial_planet_radius,
                                           initial_planet_metalicity, output_path, saved_model=saved_planet ,new_opacity_file=opacity_extension)
    initial_planet_radius = planet_model.planet.radius
    #initializes the calculator of the external energy source
    star = Particle(mass= star_mass) # could be evolved with MESA too but is constant for current purpose
    star.luminosity = star_luminosity # must be set if the star model was not evolved with a stellar evolution code
    star.radius = 1.0 | units.RSun
    # if one wants to evolve the star during the simulation its values should be updated in orbitalEvolution at each step

    #initializes the orbital evolution simulator
    orbit_model = TidalEvolution(planet=planet_model.planet, star=star,
                               initial_separation=initial_separation, initial_eccentricity=initial_eccentricity,
                               eccentricity_cuttoff=eccentricity_cutoff, tide_model=tide_model)


    print("All models are ready to evolve")
    thermal_problematic_times = []
    last_tidal_heat = 0.0 | (units.erg / units.s)
    last_radiation_heat = 0.0 | (units.erg / units.s)
    try:
        thermal_time_step = time_step
        leftover_time = 0.0 | units.yr
        evolution_time = 0.0 | units.yr
        while evolution_time < total_time and planet_model.planet.radius > radius_cutoff:
        #for evolution_time in range(0, int(total_time.value_in(units.yr)), int(time_step.value_in(units.yr))):
            plotting += str(orbit_model.separation.value_in(units.cm)) + "," + str(orbit_model.eccentricity) + \
                        "," + str(planet_model.planet.radius.value_in(units.cm)) + "," + \
                        str(planet_model.planet.central_temperature.value_in(units.K)) + "," + \
                        str(planet_model.planet.temperature.value_in(units.K)) + "," + \
                        str(planet_model.planet_pressure_profile[0].value_in(units.erg / units.cm**3)) + "," + \
                        str(planet_model.planet_pressure_profile[-1].value_in(units.erg / units.cm**3)) + "," + \
                        str(last_tidal_heat.value_in(units.erg / units.s)) + "," + \
                        str(last_radiation_heat.value_in(units.erg / units.s)) + "," + \
                        str((evolution_time).value_in(units.s))+"\n"
            orbit_model.planet = planet_model.planet # updating planet properties before another orbital evolution
            current_time_step = time_step#TODO: min(time_step, planet_model.time_scale())
            orbit_model.evolve(current_time_step) # performing orbital evolution for time_step, using tide_model
            # now ready to inject the external energy from the current parameters of the system
            last_tidal_heat = tidal_efficiency * orbit_model.tides_energy_on_planet()
            planet_model.inject_energy(last_tidal_heat,
                                       shell=planet_model.find_layer_of_radius(orbit_model.planet_bulge_from_tides()))
            last_radiation_heat = radiation_efficiency * orbit_model.irradiation_flux()
            planet_model.inject_energy(last_radiation_heat, shell=-1, add_to_current=True)
            old_planet_age = planet_model.planet.age
            thermal_time_step = current_time_step + leftover_time # next step evolve additional time as what is left, if no leftover will be just time_step again
            planet_model.evolve_planet(thermal_time_step) # evolving the internal planet
            leftover_time = (old_planet_age + thermal_time_step - planet_model.planet.age)
            if thermal_time_step != current_time_step:
                print ("different time steps for thermal and orbital evolutions")
                thermal_problematic_times.append(str((evolution_time).value_in(units.s)))

            print ("time: ", (evolution_time) + current_time_step, "eccentricity: ", orbit_model.eccentricity, \
                "semimajor: ", orbit_model.separation.as_quantity_in(units.AU), \
                "tidal heat: ", last_tidal_heat.as_quantity_in(units.erg / units.s),  \
                "radiation heat: ", last_radiation_heat.as_quantity_in(units.erg / units.s), \
                "central luminosity:",planet_model.planet.get_luminosity_profile()[0],"planet: ", planet_model.planet)

            evolution_time += current_time_step


    except Exception as ex:
        print ("couldnt reach end of simulation, ", ex)
        print ("you can check MESA output for more information about internal errors")

    if planet_model.planet.radius <= radius_cutoff:
        print ("Planet has a Jupiter radius. Simulation terminated")

    tides_str = "_tides"
    rad_str = "_radiation"
    if tidal_efficiency == 0:
        tides_str = "no" + tides_str
    else:
        tides_str = "with" + tides_str
    if radiation_efficiency == 0:
        rad_str = "no" + rad_str
    else:
        rad_str = "with" + rad_str
    plotting_path = output_path + "/resultCSV{0}.csv".format(tides_str + "_" + rad_str + "_" +
                                                             str(round(initial_planet_mass.value_in(units.MJupiter),1)) + "MJ-" +
                                                             str(round(initial_planet_radius.value_in(units.RJupiter),3)) + "RJ-" +
                                                             str(initial_separation.value_in(units.AU)) + "AU-" +
                                                             str(initial_eccentricity) +
                                                             "-time" + str(time.localtime().tm_year) + "-" +
                            str(time.localtime().tm_mon) + "-" + str(time.localtime().tm_mday) + "-" +
                            str(time.localtime().tm_hour) + ":" + str(time.localtime().tm_min) + ":" +
                            str(time.localtime().tm_sec))
    file = open(plotting_path, mode="w")
    file.write(plotting)
    try:
        if len(thermal_problematic_times) > 0:
            file.write("\n problems with time steps in the following points: ")
            file.write('\n'.join(thermal_problematic_times))
    except:
        print("couldnt write problematic timings: ", thermal_problematic_times)

    file.close()
    print ("plotted a,e,R,Teff,Tc,t to ", plotting_path)


def InitParser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--initial_model', type=str, default=None ,help='path to the initial planet model created '
                                                                        'externally with MESA. for example- .mod file.')
    parser.add_argument('--extension_opacity_file_path', type=str, default=None ,help='path to the extention of opacity '
                                                                                      'table gs98 from new MESA version, '
                                                                                      'especially for low temperature'
                                                                                      'for example- lowT_Freedman11_z0.02')
    parser.add_argument('--saved_planet', type=str, default=None ,help='path to a planet model which was created with this code')
    parser.add_argument('--initial_planet_mass', type=float, default=1, help='if no initial model provided- create a new preMS '
                                                           'star with this mass in MJupiter')
    parser.add_argument('--initial_planet_radius', type=float, default=3, help='if no initial model provided- '
                                                             'evolve this preMS to this radius in RJupiter')
    parser.add_argument('--star_mass', type=float, default=1.0, help='stellar mass in MSun')
    parser.add_argument('--star_luminosity', type=float, default=1.0, help='stellar luminosity in LSun')
    parser.add_argument('--initial_planet_metalicity', type=float,default=0.02, help='if no initial model provided- take this metallicity')
    parser.add_argument('--evolution_time', type=float, required=True, help='total evolution time in Myrs')
    parser.add_argument('--time_step', type=float, required=True, help='time of each step in yrs')
    parser.add_argument('--initial_separation', type=float, required=True, help='initial separation between star and '
                                                                                'planet in AU')
    parser.add_argument('--initial_eccentricity', type=float, required=True, help='initial eccentricity between star and '
                                                                                'planet')

    parser.add_argument('--eccentricity_cutoff', type=float, default=10**-3, required=False, help='lower limit of the eccentricity to stop the simulation')
    parser.add_argument('--radius_cutoff', type=float, default=1.09582, required=False,
                        help='lower limit of the planet radius in RJupiter')
    parser.add_argument('--tidal_efficiency', type=float, default=1, required=False, help='efficiency for tidal energy')
    parser.add_argument('--irradiation_efficiency', type=float, default=1, required=False, help='efficiency for radiation energy')
    parser.add_argument('--save_snapshots', type=bool, default = False, help='whether we should save the snapshots or not, '
                                                           'will be saved in the output dir')
    parser.add_argument('--output_path', type=str, required=True, help='path to simulations output directory')
    parser.add_argument('--tides_model', type=str, required=False, default="weak", help='tides model for the calculation, default is weak tides using Hut81')
    return parser


if __name__ == "__main__":
    parser = InitParser()
    args = parser.parse_args()
    main(args.initial_model, args.extension_opacity_file_path, args.saved_planet, args.initial_planet_mass | units.MJupiter,
         args.initial_planet_radius | units.RJupiter,
         args.star_mass | units.MSun, args.star_luminosity | units.LSun, args.initial_planet_metalicity,
         args.output_path, args.initial_separation | units.AU, args.initial_eccentricity, args.eccentricity_cutoff,
         args.radius_cutoff | units.RJupiter,
         (args.evolution_time * (10**6)) | units.yr, args.time_step | units.yr, args.tides_model,
         args.tidal_efficiency, args.irradiation_efficiency)
