import os
import time
import shutil
import math
import numpy as np
import re
import pickle
from amuse.community.mesa.interface import MESA
from amuse.units import units, constants
from amuse.support.exceptions import AmuseException
from amuse.datamodel import Particle, Particles

find_nums_in_row = lambda row: np.array([float(x) for x in re.findall('[\-\.0-9]+', row)])

class InternalPlanetEvolution:
    def __init__(self,initial_mod_file,initial_mass, initial_radius, initial_metalicity, output_path,
                 saved_model = None, new_opacity_file=None):
        mesa_output_file = os.path.join(output_path,"mesa_output_{0}.log".format
                            (str(time.localtime().tm_year) + "-" +
                            str(time.localtime().tm_mon) + "-" + str(time.localtime().tm_mday) + "-" +
                            str(time.localtime().tm_hour) + ":" + str(time.localtime().tm_min) + ":" +
                            str(time.localtime().tm_sec)))
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        try:
            self.evolutionCode = MESA(redirection='file',redirect_file=mesa_output_file)
            print("plotting MESA logs to ", mesa_output_file)
        except:
            self.evolutionCode = MESA()

        new_mesa_path = output_path + '/MESA'
        self.mesa_parh = new_mesa_path
        self.copy_data_dir_to_new_mesa_path(new_mesa_path)

        if initial_mod_file is not None:
            self.copy_model_into_zams(initial_mod_file,
                                      zams_to_extend= new_mesa_path + '/data/star_data/starting_models/zams_z' +
                                                     str(int(initial_metalicity*100))+'m2.data',
                                      new_name= new_mesa_path + '/data/star_data/starting_models/zams_z' +
                                                str(int(initial_metalicity*1000))+'m3.data')

        self.minimal_Teff = 10 ** 2.1 | units.K # for OPAL EOS, will be changed by concatanation
        self.kappa_prefix = 'gs98'
        if new_opacity_file is not None:
            self.concat_new_opacities('gs98',new_opacity_file,'GRPG20')
            self.kappa_prefix = 'GRPG20'

        print("creating new inlist file")
        self.copy_and_change_inlist(initial_mod_file, self.evolutionCode, new_mesa_path)
        self.evolutionCode.set_MESA_paths(new_mesa_path + '/inlist',new_mesa_path+'/data', new_mesa_path+'/data')
        self.set_zero_winds()
        self.evolutionCode.set_min_timestep_stop_condition(0.0 | units.s)
        self.evolutionCode.initialize_code()
        print("code has been initialized")
        self.planet = None
        if saved_model is not None:
            try:
                with open(saved_model,'rb') as saved_planet_file:
                    initial_planet = pickle.load(saved_planet_file)
                    self.evolutionCode.parameters.stabilize_new_stellar_model_flag = False
                    self.planet = self.evolutionCode.new_particle_from_model(initial_planet,initial_planet.age)
                    print("saved planet model loaded")
            except Exception as ex:
                print("couldn't load saved model")
        if self.planet is None:
            print("creating planet")
            #will still create a pre main sequence particle if the provided profile is of a very low mass star
            self.planet = self.evolutionCode.particles.add_particle(Particle(mass=initial_mass))
        self.set_zero_winds()
        self.evolutionCode.set_min_timestep_stop_condition(0.0 | units.s)

        print(self.planet)
        print("now evolving to the initial radius")
        while self.planet.radius > initial_radius:
            self.planet.evolve_one_step()
            print(self.planet.radius)
        model_file_name = str(self.planet.mass.value_in(units.MJupiter)) + "MJupiter_" + \
                          str(round(self.planet.radius.value_in(units.RJupiter),3))+"RJupiter"
        if not os.path.exists(os.path.join(output_path,model_file_name)):
            with open(os.path.join(output_path,model_file_name), 'wb') as openedFile:
                pickle.dump(StellarModel(self.planet), openedFile, pickle.HIGHEST_PROTOCOL)
                print("model saved to"+os.path.join(output_path,model_file_name))
        print("planet is ready")
        self.planet_pressure_profile = self.planet.get_pressure_profile()
        print(self.evolutionCode.particles)

    def set_zero_winds(self):
        self.evolutionCode.set_RGB_wind_scheme(0)
        self.evolutionCode.set_AGB_wind_scheme(0)
        self.evolutionCode.set_reimers_wind_efficiency(0.0)
        self.evolutionCode.set_blocker_wind_efficiency(0.0)
        self.evolutionCode.set_de_jager_wind_efficiency(0.0)
        self.evolutionCode.set_dutch_wind_efficiency(0.0)

    def copy_data_dir_to_new_mesa_path(self, new_mesa_path):
        print("copying mesa data dir into new working dir")
        if os.path.exists(new_mesa_path+'/data'):
            print("path "+new_mesa_path+"/data already exists, no need to copy")
            return
        try:
            os.makedirs(new_mesa_path)
        except(OSError):
            pass
        try:
            shutil.copytree(self.evolutionCode.default_path_to_MESA_data, new_mesa_path+'/data')
            print("succesfully coppied")
        except OSError as error:
            print(error)
            print(error.message)

    def copy_and_change_inlist(self, initial_model, evolution_code_instance, new_path):
        new_inlist = new_path + '/inlist'
        if os.path.exists(new_inlist):
            os.remove(new_inlist)
        if initial_model is not None:
            try:
                shutil.copyfile(initial_model, os.path.join(new_path,os.path.basename(initial_model)))
            except(OSError):
                print("there was a problem copying the saved model")
                print (OSError.message)
            line_found = False
            done_changing = False
        else:
            line_found = True
            done_changing = True
        with open(evolution_code_instance.default_path_to_inlist, 'r') as old_inlist_file:
            with open(new_inlist, 'w') as new_inlist_file:
                for line in old_inlist_file:
                    if not line_found:
                        if "load_saved_model = " in line:
                            line = line.replace("false", "true")
                            line_found = True
                    elif not done_changing:
                        splited = line.split("saved_model_name = ") # not really needed, amuse doesnt work like that,
                                                                    ## inserted this mod data into zams file instead
                        if len(splited) > 1:
                            line = ''.join(line.split(splited[1]) + [ "'" + os.path.basename(initial_model) + "'", '\n'])
                            done_changing = True
                    #if 'report' in line:
                    #    line = line.replace('false','true')
                    if 'timestep_limit' in line:
                        line = line.replace('1d-6','1d-9')
                    if 'kappa_file_prefix' in line:
                        line = line.replace('gn93',self.kappa_prefix)
                        #new_inlist_file.write("kappa_lowT_prefix = 'lowT_Freedman11'\n")
                    if 'RGB_wind_scheme' in line:
                        line = line.replace('1', '0')
                        line = line.replace('2', '0')
                    if 'AGB_wind_scheme' in line:
                        line = line.replace('1', '0')
                        line = line.replace('2', '0')
                    if 'T_mix_limit' in line:
                        line = line.replace('1d4', '0')
                    if 'wind_envelope_limit' in line:
                        line = line.replace('-1', '2')
                    if 'logQ_limit' in line:
                        line = line.replace('3.9','4.3')
                    if 'change_net' in line:
                        line = line.replace('false','true')
                    if 'new_net_name' in line:
                        line = line.replace("new_net_name = ''","new_net_name = 'basic.net'")
                    new_inlist_file.write(line)


                new_lines = [
                "change_initial_net =.true.",
                "eos_file_prefix = 'mesa'",
                "change_lnPgas_flag = .true.",
                "new_lnPgas_flag = .true."
                ]
                new_inlist_file.writelines(new_lines)

    def copy_model_into_zams(self, model_file=None, zams_file=None,zams_to_extend = None,new_name=None):
        '''
        This is a patch to import a new MESA model into the old version used by AMUSE.
        The model is converted into the old version of .data zams file and is inserted in the right place in the file
        zams_to_extend.
        If a ready zams_file is supplied it will just copy it to the working by AMUSE directory.

        :param model_file: a model file created with old or new MESA version
        :param zams_file: ready zmas file to be used
        :param zams_to_extend: a zams_file of the old version, to be extened with the supplied model
        '''
        if zams_file is not None:
            print("coping ready zams file")
            shutil.copyfile(zams_file, os.path.join(self.evolutionCode.default_path_to_MESA_data, os.path.basename(zams_file)))
            print("successfully coppied")
        elif model_file is not None and zams_to_extend is not None:
            print("extending old zmas file with supplied model")
            mass = 0
            n_shells = 0
            matrix_began = False
            matrix_end = False
            matrix_lines = []
            with open(model_file, 'r') as original_model:
                for line in original_model.readlines():
                    if matrix_end: # no need for the extra lines at the end of mod file
                        if len(matrix_lines) > 1:
                            matrix_lines.append("\n")
                        break
                    if mass == 0 and "M/Msun" not in line:
                        continue
                    if "M/Msun" in line:
                        mass_str = [val for val in line.split(" ") if val != ""][-1].replace('D','E')
                        mass = mass_str.replace('\r','').replace('\n', '')
                        matrix_lines.append("                          M/Msun      " + mass_str)
                        continue
                    if n_shells == 0 and "n_shells" not in line:
                        continue
                    if "n_shells" in line:
                        n_shells_str = [val for val in line.split(" ") if val != ""][-1]
                        n_shells = n_shells_str.replace('\r','').replace('\n', '')
                        matrix_lines.append("                        n_shells        " + n_shells_str)
                        matrix_lines.append("\n")
                        continue
                    if not matrix_began:
                        if "lnd" not in line:
                            continue
                        else:
                            matrix_began = True
                            vals_indices_to_remove = []
                            headers = [header for header in line.split(' ') if header != ""]
                            zams_headers = ["lnd","lnT","lnR","L","dq","h1","he3","he4","c12","n14","o16","ne20","mg24"]
                            for i in range(len(headers)):
                                if headers[i] not in zams_headers:
                                    vals_indices_to_remove.append(i+2) # the plus 2 is for the extra column of the shell which doesnt have a header and the count which starts with 1
                            matrix_lines.append("  ".join(zams_headers))
                            continue
                    if int([val for val in line.split(' ') if val != ""][0]) == int(n_shells): #this is the last line of the table
                        matrix_end = True
                    line = line.replace('D','E')
                    still_working = False
                    count_vals = 0
                    new_line = ""
                    prev_i = -3 # remain same spaces for the first val in line, i,e- the shell index
                    for internal_i in range(len(line)):
                        if line[internal_i] == " ":
                            still_working = False
                            continue

                        if not still_working:
                            count_vals += 1

                        if count_vals in vals_indices_to_remove:
                            still_working = True
                            prev_i = internal_i
                            continue

                        if still_working:  # working on previous wanted data
                            new_line += line[internal_i]
                            prev_i = internal_i
                            continue

                        # new wanted value, adding previoud spaces before
                        new_line += ''.join([" " for r in range(internal_i - prev_i - 3)])
                        prev_i = internal_i
                        new_line += line[internal_i]
                        still_working = True

                    if new_line[-1] != '\n':
                        new_line += '\n'
                    matrix_lines.append(new_line)

            if zams_to_extend is not None:
                to_write = []
                with open(zams_to_extend,'r') as opened_zams:
                    lines = opened_zams.readlines()
                    to_write.extend(lines[:8])
                    for i in range(8,len(lines)):
                        curr_mass = [val for val in lines[i].split(' ') if val != ""][0].replace('\r', '').replace('\n','')
                        if float(curr_mass) == float(mass):
                            print("zams file already contains the model")
                            return
                        if float(curr_mass) > float(mass):
                            mass_to_write = str(round(float(mass),8))
                            if len(mass_to_write.split('.')[-1])<8:
                                mass_to_write = str(round(float(mass),9))[:-1]
                            to_write.append("		"+mass_to_write+"			  "+n_shells_str)
                            break
                    pasted = False
                    while "M/Msun" not in lines[i]:
                        to_write.append(lines[i])
                        i += 1

                    matrix_lines[3] = lines[i+3]
                    #we now need to paste our model
                    while i <=len(lines):
                        curr_mass = [val for val in lines[i].split(' ') if val != ""][-1].replace('\n','').replace('\r','')
                        if float(curr_mass) == float(mass):
                            print("zams file already contains the model")
                            return
                        elif float(curr_mass) > float(mass):
                            to_write.extend(matrix_lines)
                            to_write.extend(lines[i:])
                            break
                        else:
                            to_write.append(lines[i])
                            i += 1
                            to_write.append(lines[i])
                            lines_count = int([val for val in lines[i].split(' ')
                                               if val != ""][-1].replace('\r','').replace('\n','')) + 3
                            to_write.extend(lines[i:lines_count])
                            i += lines_count

                with open(new_name,'w') as opened_zams:
                    opened_zams.writelines(to_write)

                print("new zams file was created")


    def concat_new_opacities(self, kappa_old_prefix, kappa_file_2,kappa_new_prefix,logT_bnd=3.75):
        '''
        Forming new opacity files for each file with the kappa_old prefix, extended by kappa_file_2
        :param kappa_old_prefix: prefix of original kappa files
        :param kappa_file_2: extension to the old kappa files
        :param kappa_new_prefix: the new prefix to be used with MESA via AMUSE
        :param logT_bnd: temperature boundary between two opacity tables
        :return:
        '''
        kap_data = os.path.join(self.mesa_parh,"data","kap_data")
        z = float(kappa_file_2.split('_')[-1].split('.data')[0].replace('z',''))
        m = int(math.ceil(abs(math.log10(z))))
        z = int(z*(10**m))
        kappa_files_to_copy_asis = [os.path.join(kap_data,f_name) for f_name in os.listdir(kap_data) if kappa_old_prefix in f_name
                       and ("_co_" in f_name or "z"+str(z)+"m"+str(m) not in f_name)]
        for file_to_copy in kappa_files_to_copy_asis:
            try:
                shutil.copyfile(file_to_copy, file_to_copy.replace(kappa_old_prefix,kappa_new_prefix))
            except OSError as ex:
                print("there was a problem copying file")
                print (ex)

        kappa_files = [os.path.join(kap_data,f_name) for f_name in os.listdir(kap_data) if kappa_old_prefix in f_name
                       and "_co_" not in f_name and "z"+str(z)+"m"+str(m) in f_name]
        extra_kappa = KappaFileData(kappa_file_2)
        print ("extending kappa files")
        for kappa_file in kappa_files:
            current_kappa = KappaFileData(kappa_file, prefix=kappa_old_prefix)
            new_kappa_file = KappaFileData(kappa_file, prefix=kappa_old_prefix)
            new_kappa_file.path = kappa_file.replace(kappa_old_prefix,kappa_new_prefix)
            new_kappa_file.prefix = kappa_new_prefix
            new_kappa_file.logT_min = min(current_kappa.logT_min,extra_kappa.logT_min)
            self.minimal_Teff = max(self.minimal_Teff, 10 ** new_kappa_file.logT_min | units.K)
            new_kappa_file.logT_max = max(current_kappa.logT_max,extra_kappa.logT_max)
            table_lines = []
            #now choose which one comes before the other
            if new_kappa_file.logT_min == extra_kappa.logT_min:
                first = extra_kappa
                second = current_kappa
            else:
                first = current_kappa
                second = extra_kappa
            #find the intersection
            first_end_i = 0
            for i in range(7,len(first.lines)):
                if find_nums_in_row(first.lines[i])[0] >= logT_bnd:
                    first_end_i= i
                    break

            second_start_i = 0
            for i in range(7,len(second.lines)):
                if find_nums_in_row(second.lines[i])[0] >= logT_bnd:
                    second_start_i = i
                    break
            new_kappa_file.logTs = int((first_end_i -1 - 7) + second.logTs - (second_start_i-7))
            #put in the data
            table_lines.extend(first.read_partial_lines(new_kappa_file.logR_line)[:first_end_i-1])
            table_lines.extend(second.read_partial_lines(new_kappa_file.logR_line)[second_start_i:])
            #replace first lines with new values
            table_lines[0] = table_lines[0].replace('\n','').replace('\r','') + " and " + second.lines[0]
            table_lines[2] = table_lines[2].replace(first.X,new_kappa_file.X)
            table_lines[2] = table_lines[2].replace(first.Z,new_kappa_file.Z)
            table_lines[2] = table_lines[2].replace(str(first.logR_min),str(new_kappa_file.logR_min))
            table_lines[2] = table_lines[2].replace(str(first.logR_max),str(new_kappa_file.logR_max))
            table_lines[2] = table_lines[2].replace(str(first.logT_min),str(new_kappa_file.logT_min))
            table_lines[2] = table_lines[2].replace(str(first.logT_max),str(new_kappa_file.logT_max))
            table_lines[2] = table_lines[2].replace(str(first.logRs),str(new_kappa_file.logRs))
            table_lines[2] = table_lines[2].replace(str(first.logTs),str(new_kappa_file.logTs))
            if new_kappa_file.logR_max == second.logR_max:
                table_lines[5] = second.lines[5]
            #save the new file
            with open(new_kappa_file.path, 'w') as new_kappa_file_opened:
                new_kappa_file_opened.writelines(table_lines)

        print ("all kappa files saved")



    def evolve_planet(self, time_step, backups=100):
        if backups == 0:
            print("number of MESA backups reached, cannot evolve further")
            return #raise Exception("number of MESA backups reached, cannot evolve further")

        if time_step <= 0.0 | units.yr:
            print ("evolved an extra ", (-1.0 * time_step))
            return

        if self.planet.central_temperature == 1.0 | units.K:
            print ("evolution with MESA terminated")
            raise Exception("evolution with MESA terminated")

        if self.planet.temperature <= self.minimal_Teff:
            print ("evolution with MESA terminatied due to lower limit of temperature as the minimal in the opacity table")
            return

        print("evolving lonely planet for " + str(time_step))
        self.planet.reset_number_of_backups_in_a_row()
        evolved_age = self.planet.age + time_step
        self.planet.time_step = time_step
        old_planet = self.planet.copy()
        try:
            self.planet.evolve_for(time_step)
        except AmuseException as ex:
            if "backups" not in str(ex):
                print("cannot evolve further with MESA ", ex)
                if self.planet.age == old_planet.age:
                    if self.planet.central_temperature != old_planet.central_temperature:
                        # TODO: recover old model somehow, copy doesnt work
                        print ("central temperature changed evethough the age didnt check MESA log file")
                return

            print("increasing backups")
            self.evolve_planet(evolved_age - self.planet.age, backups - 1)

        #self.planet.time_step = time_step

    def find_layer_of_radius(self, radius):
        radius_profile = self.planet.get_radius_profile()
        for i in range(len(radius_profile)):
            if radius_profile[i] > radius:
                if i > 0 and (radius - radius_profile[i-1]) < (radius_profile[i] - radius):
                    return i-1
                else:
                    return i

    def inject_energy(self, extra_energy, shell=0, add_to_current=False):
        radius_profile = self.planet.get_radius_profile()
        density_profile = self.planet.get_density_profile()
        self.planet_pressure_profile = self.planet.get_pressure_profile()
        dmass_profile = self.planet.get_mass_profile() * self.planet.mass
        gravity_acceleration = constants.G * dmass_profile.cumsum() / self.planet.radius**2
        small_scale_height_profile = np.sqrt(self.planet_pressure_profile / (constants.G * density_profile**2))
        big_scale_height_profile = self.planet_pressure_profile/(density_profile*gravity_acceleration)
        scale_height_profile = small_scale_height_profile
        scale_height_profile[0] = self.planet.radius
        scale_height_profile[1:] = np.minimum(small_scale_height_profile[1:], big_scale_height_profile[1:])
        sigma = 0.5 * scale_height_profile[shell]
        print(sigma.as_quantity_in(units.RSun) , "around shell=" , shell)
        extra_profile = extra_energy * (1.0 / \
                        (np.sqrt(2. * np.pi * sigma ** 2) * 4.0 * np.pi * density_profile * radius_profile**2)) * \
                        np.exp(-0.5 * ((radius_profile - radius_profile[shell]) / sigma) ** 2 )
        print ((extra_profile * dmass_profile).sum().as_quantity_in(units.erg/units.s), extra_energy.as_quantity_in(units.erg/units.s))
        if add_to_current:
            current_extra = self.planet.get_extra_heat_profile()
            extra_profile = extra_profile + current_extra

        self.planet.set_extra_heat_profile(extra_profile)

    def time_scale(self):
        '''
        :return: Thermal timescale
        '''
        return (constants.G * self.planet.mass**2)/(2*self.planet.radius*self.planet.luminosity)


class KappaFileData:
    def __init__(self, path, prefix = None):
        self.prefix = prefix
        self.path = path
        if not os.path.exists(path):
            print (path + " doesnt exist")
            return
        with open(path, 'r') as kappa_file:
            self.lines = kappa_file.readlines()
            self.indices_line = self.lines[2]
            self.logR_line = self.lines[5]
            self.retrieve_special_line_params()

    def retrieve_special_line_params(self):
        splitted = [obj for obj in self.indices_line.split(' ') if obj != ""]
        self.X = splitted[2]
        self.Z = splitted[3]
        self.logRs = int(splitted[4])
        self.logR_min = float(splitted[5])
        self.logR_max = float(splitted[6])
        self.logTs = int(splitted[7])
        self.logT_min = float(splitted[8])
        self.logT_max = float(splitted[9])
        self.logR_intervals =  (self.logR_max - self.logR_min)/(self.logRs-1)

    def read_lines(self):
        with open(self.path,'r') as opened_file:
            self.lines = opened_file.readlines()


    def read_partial_lines(self, logR_line):
        logR_values_for_int = find_nums_in_row(logR_line)
        logR_values = find_nums_in_row(self.lines[5])
        if len(logR_values) == len(logR_values_for_int) and \
                len([logR_values[i] for i in range(len(logR_values)) if logR_values[i] == logR_values_for_int[i]]) > 0:
            return self.lines #no need to resample

        cutted_lines = self.lines[:7]
        printing_format = '{:-11.6f}' # format of old MESA kappa files

        for i in range(7,len(self.lines)):
            if not self.lines[i].strip():
                cutted_lines.append(self.lines[i])
                continue
            line_nums = find_nums_in_row(self.lines[i])
            new_line_values = np.interp(logR_values_for_int, logR_values,line_nums[1:]) #interpolating for the wanted sampling values
            new_line = printing_format.format(line_nums[0])+' '*4+''.join(printing_format.format(n) for n in new_line_values)
            if new_line[-1] != '\n':
                new_line += '\n'
            cutted_lines.append(new_line)

        return cutted_lines



class StellarModel:
    def __init__(self, star):
        self.radius = star.get_radius_profile()
        self.rho = star.get_density_profile()
        self.temperature = star.get_temperature_profile()
        self.luminosity = star.get_luminosity_profile()
        composition = star.get_chemical_abundance_profiles()
        self.composition = composition
        self.X_H = composition[0]
        self.X_He = composition[1] + composition[2]
        self.X_C = composition[3]
        self.X_N = composition[4]
        self.X_O = composition[5]
        self.X_Ne = composition[6]
        self.X_Mg = composition[7]
        self.X_Si = composition[7]*0.0
        self.X_Fe = composition[7]*0.0
        self.dmass = star.get_mass_profile() * star.mass
        self.age = star.age


