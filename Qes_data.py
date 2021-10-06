#!/usr/bin/env python3
import os
import glob
import re
import numpy as np
import math
import pandas as pd
import h5py
import subprocess
from shutil import copyfile, copy2
import autocorr_pydiag
import autocorr_johannes
import sim_launcher

def get_num_species():
    return int(os.environ.get('ADAPTATION_NUMBER_OF_SPECIES'))

def read_nrg_file(path):
    if get_num_species() == 2:
        times = pd.read_table(path, sep="\s+",
                            header=None, skiprows=lambda x: x % 3 != 0)[0]
        names=["n", "u", "T_\parallel", "T_\perp", "Gamma_es", "Gamma_em", "Q_es", "Q_em"]
        frame0 = pd.read_table(path, sep="\s+",
                            header=None, names=names,
                            skiprows=lambda x: x % 3 != 1)
        frame1 = pd.read_table(path, sep="\s+",
                            header=None, names=names,
                            skiprows=lambda x: x % 3 != 2)
        frame = frame0.join(frame1, lsuffix='0', rsuffix='1')
    elif get_num_species() == 1:
        times = pd.read_table(path, sep="\s+",
                            header=None, skiprows=lambda x: x % 2 != 0)[0]
        names=["n0", "u0", "T_\parallel0", "T_\perp0", "Gamma_es0", "Gamma_em0", "Q_es0", "Q_em0"]
        frame = pd.read_table(path, sep="\s+",
                            header=None, names=names,
                            skiprows=lambda x: x % 2 != 1)
    else:
        raise NotImplementedError
    frame["time"] = times
    return frame

# to append nrg.dat slices from consecutive runs, put them in a list


def file_list_to_data_frame(filelist):
    frame = pd.DataFrame()
    list_ = []
    for file_ in filelist:
        assert os.path.isfile(file_)
        csv = read_nrg_file(file_)
#         csv.set_index("time",drop=False)
        list_.append(csv)
    frame = pd.concat(list_, ignore_index=True).drop_duplicates(
        subset="time").reset_index(drop=True)
#     frame = frame.set_index("time",drop=False)
    return frame


def get_in_sim_time_per_timestep(frame):
    # in the nrg file, every entry stands for 100 time steps #TODO read the actual nrg spacing from parameters template
    num_timesteps = 100 * len(frame['time'])
    in_sim_time = frame['time'].iloc[-1] - frame['time'].iloc[0]
    return in_sim_time / num_timesteps


def get_time_per_timestep(gene_last_output_file):
    try:
        with open(gene_last_output_file, 'r') as infile:
            text = infile.read()
            pattern_tasks = r'We have \s* (.*)\s*MPI tasks.'
            match_tasks = re.search(pattern_tasks, text)
            num_tasks = int(match_tasks.group(1))

            pattern_time = r'Time per time step:\s*(.*)\s*sec'
            match_time = re.search(pattern_time, text)
            time = float(match_time.group(1))
        return time*num_tasks
    except AttributeError:
        print("AttributeError")
        return None


def get_cost_per_time(gene_last_output_file, frame):
    try:
        in_sim_time_per_timestep = get_in_sim_time_per_timestep(frame)
        time_per_timestep = get_time_per_timestep(gene_last_output_file)
        return time_per_timestep/in_sim_time_per_timestep
    except (AttributeError, TypeError):
        return None


def get_autocorrelation_time(frame):
    assert frame['time'].is_monotonic
    t_p = [0.] * get_num_species()
    t_j = [0.] * get_num_species()
    for species in range(get_num_species()):
        t_j[species], _, _, _ = autocorr_johannes.get_autocorr_time(
            frame['time'].copy().values, frame['Q_es'+str(species)].copy().values)
        t_p[species] = autocorr_pydiag.autocorrtime_1d(
            frame["Q_es"+str(species)].copy(), frame['time'].copy())
        relative = t_j[species]/t_p[species]
        if relative > 1.11 or relative < 0.9:
            print("autocorrelations :" + str(t_j) + " and " + str(t_p) + \
                " -- there is something weird going on potentially")
    return t_p


def get_time_error(frame, printOutput=False, autocorr_time=None):
    # if not autocorr_time:
    #    autocorr_time = get_autocorrelation_time(frame)
    serror = [0.] * get_num_species()
    for species in range(get_num_species()):
        _, serror[species], nwin = autocorr_johannes.get_mean_error(frame['time'].copy().values, 
            frame['Q_es'+str(species)].copy().values, printOutput=printOutput)
    # Q_em is always 0 in our examples
    # Alejandro says: does not need to be > 0assert frame['Q_es'].ge(frame["Q_em"]).all()
    return serror


class Qes_data:
    def __init__(self, level_vector, \
            prob_prepath=os.environ.get('ADAPTATION_PROB_PATH'),\
            gene_path=os.environ.get('ADAPTATION_GENE3D_PATH'), skiplevels=[]):
        self.gene_path = gene_path
        if ~(gene_path == None):
            assert os.path.isdir(gene_path)
        self.level_vector=level_vector
        self.prob_directory = sim_launcher.get_prob_path(prob_prepath, level_vector)
        self.frame = None
        self.crop_time = float(os.environ.get('ADAPTATION_PARAM_CROP_TIME'))
        self.results_csv = os.environ.get('ADAPTATION_RESULTS_CSV')
        self.skiplevels=skiplevels

    def set_csv(self, delta=None):
        return self.set_result_to_csv(self.get_result(), self.frame.time.iloc[0],
                       self.how_long_run(),
                       cost_per_time=get_cost_per_time(self.get_last_output_file(), self.frame),
                       autocorrelation_time=get_autocorrelation_time(self.frame),
                       time_error=self.get_time_error(), delta=delta)

    def set_result_to_csv(self, qes, qes_from, qes_to, cost_per_time=None, autocorrelation_time=None, time_error=None, delta=None, geometry=''):
        current_results = pd.read_csv(self.results_csv, index_col=0)
        d = {}
        lkeys = [str('l_'+str(dim)) for dim in range(len(self.level_vector))]
        for dim in range(len(self.level_vector)):
            d[lkeys[dim]] = int(self.level_vector[dim])
        d['qes0'] = qes[0]
        d['qes_from'] = qes_from
        d['qes_to'] = qes_to
        d['cost_per_time'] = cost_per_time
        d['time_error0'] = np.nan if time_error is None else time_error[0]
        d['autocorrelation_time0'] = np.nan if autocorrelation_time is None else autocorrelation_time[0]
        d['delta0'] = np.nan if delta is None else float(delta[0])
        if get_num_species() == 2:
            d['qes1'] = qes[1]
            d['time_error1'] = np.nan if time_error is None else time_error[1]
            d['autocorrelation_time1'] = np.nan if autocorrelation_time is None else autocorrelation_time[1]
            d['delta1'] = np.nan if delta is None else float(delta[1])

        new_index = geometry + sim_launcher.l_vec_to_string(self.level_vector)
        #print(current_results.loc[new_index])
        current_results.loc[new_index]=d
        try:
            current_results[lkeys] = current_results[lkeys].astype(int)
        except ValueError:
            pass
        current_results.to_csv(self.results_csv)

    def set_delta_to_csv(self, delta, geometry=''):
        return self.set_csv(delta)

    def get_from_csv(self, level_vector):
        current_results = pd.read_csv(self.results_csv, index_col=0)
        ew_index = sim_launcher.l_vec_to_string(level_vector)
        if new_index in current_results.index:
            # print(new_index)
            return current_results.loc[new_index, :]
        else:
            print(new_index + " not in csv!")
            return None

    def get_cost_per_time_csv(self, level_vector):
        return self.get_from_csv(level_vector)['cost_per_time']

    def get_autocorrelation_time_csv(self, level_vector):
        return self.get_from_csv(level_vector)['autocorrelation_time']

    def get_time_error_csv(self):
        return self.get_from_csv(self.level_vector)['time_error']

    def get_last_output_file(self):
        outfiles = glob.glob(os.path.join(self.prob_directory, '*.out'))
        outfiles.sort()
        return outfiles[-1]

    def get_frame(self, nrg_append=None):
        if self.frame is not None:
            return self.frame
        if (self.gene_path is None):
            return None
        # get all nrg_ files, filter out the ones that don't exist
        if nrg_append:
            nrg_files = [self.prob_directory + '/out/nrg_' +
                         nrg_append + str(y) for y in range(1, 30)]
        else:
            nrg_files = [self.prob_directory +
                         '/out/nrg_' + str(y) for y in range(1, 30)]
        # assume that nrg.dat is the latest
        nrg_files.append(self.prob_directory + '/out/nrg.dat')
        nrg_files = list(filter(lambda i: os.path.isfile(i), nrg_files))

        if not nrg_files:
            self.frame = None
        else:
            print(nrg_files)
            self.frame = file_list_to_data_frame(nrg_files)
            # crop away the peak
            self.frame = self.frame.loc[self.frame['time'] >= self.crop_time]
            self.frame.reset_index(inplace=True, drop=True)
            if len(self.frame) < 0:
                self.frame = None
        return self.frame

    def get_flux(self, level_vector, interpolate=False, interpolate_onto_number_of_points=1000):
        raise NotImplementedError("This was not reworked for em simulations!")
        probname = sim_launcher.l_vec_to_string(level_vector) 
        flux_filename = './flux_diags/flux_profile_ions_' + str(probname) + '.h5'

        with h5py.File(flux_filename,  "r") as f:
            for item in f.attrs.keys():
                print(item + ":", f.attrs[item])
            Q_es_ions = f['Q_es_ions']
            x_a = f['x_a']
            d = {'Q_es_ions': np.array(Q_es_ions), 'x_a': np.array(x_a)}
            flux = pd.DataFrame(data=d)
        if interpolate:
            x_a_range = [flux['x_a'].min(), flux['x_a'].max()]
            x_a = np.linspace(
                x_a_range[0], x_a_range[1], interpolate_onto_number_of_points)

            # cf. https://stackoverflow.com/questions/10464738/interpolation-on-dataframe-in-pandas
            flux.set_index('x_a', inplace=True, drop=False)

            # Interpolation technique to use. One of:

            # 'linear': Ignore the index and treat the values as equally spaced. This is the only method supported on MultiIndexes.
            # 'time': Works on daily and higher resolution data to interpolate given length of interval.
            # 'index', 'values': use the actual numerical values of the index.
            # 'pad': Fill in NaNs using existing values.
            # 'nearest', 'zero', 'slinear', 'quadratic', 'cubic', 'spline', 'barycentric', 'polynomial': Passed to scipy.interpolate.interp1d. These methods use the numerical values of the index. Both 'polynomial' and 'spline' require that you also specify an order (int), e.g. df.interpolate(method='polynomial', order=5).
            # 'krogh', 'piecewise_polynomial', 'spline', 'pchip', 'akima': Wrappers around the SciPy interpolation methods of similar names. See Notes.
            # 'from_derivatives': Refers to scipy.interpolate.BPoly.from_derivatives which replaces 'piecewise_polynomial' interpolation method in scipy 0.18.

    #         flux = flux.reindex(flux.index.union(Xresampled)).interpolate('linear').loc[Xresampled]
            flux = flux.reindex(flux.index.union(x_a)).interpolate(
                'polynomial', order=5).loc[x_a]
            flux.reset_index(inplace=True, drop=True)

        return flux, x_a

    def get_bin_averaged_flux(self, level_vector, number_of_bins=2, interpolate=False):
        raise NotImplementedError("This was not reworked for em simulations!")
        try:
            flux, x_a = self.get_flux(level_vector, interpolate)
        except OSError:
            return None

        x_a_range = [flux['x_a'].min(), flux['x_a'].max()]
        bin_bounds = np.linspace(x_a_range[0], x_a_range[1], number_of_bins+1)
        bin_frames = [flux[(flux['x_a'] > bin_bounds[i]) & (flux['x_a']
                           <= bin_bounds[i+1])] for i in range(number_of_bins)]

        if interpolate:
            bin_averaged_flux = [np.mean(bin_frames[i]['Q_es_ions'])
                                 for i in range(number_of_bins)]
        else:
            bin_averaged_flux = [autocorr_johannes.get_mean(bin_frames[i]['x_a'].copy().values,
                                                            bin_frames[i]['Q_es_ions'].copy().values) for i in range(number_of_bins)]

        #print(bin_averaged_flux)
        return bin_averaged_flux

    def how_long_run(self):
        frame = self.get_frame()
        if frame is not None and len(frame) > 0:
            return frame.time.iloc[-1]
        else:
            return 0

    def get_result(self, up_to_time=None):
        if self.level_vector in self.skiplevels:
            if get_num_species() == 2:
                result = (math.inf, math.inf)
            elif get_num_species() == 1:
                result = (math.inf)

        if self.gene_path is None:
            result = self.get_from_csv()
            if result is not None:
                if get_num_species() == 2:
                    result = (result['qes0'], result['qes1'])
                elif get_num_species() == 1:
                    result = (result['qes0'])
        else:
            self.frame = self.get_frame()
            if up_to_time:
                # crop away later values
                self.frame = self.frame.loc[self.frame['time'] <= up_to_time]
            if self.frame is None or len(self.frame) == 0:
                result = None
                frame = None
            else:
                # have a minimum length to base predicitons on
                if self.frame.loc[self.frame['time'] > 2*self.crop_time].empty:
                    result = None
                else:
                    # how much do the times vary? raise error if more than 20%
                    # assert self.frame.time.diff().std() < 0.2 * self.frame.time.diff().mean(),\
                    #        "mean is " + str(self.frame.time.diff().mean()) +\
                    #        " and deviation " + str(self.frame.time.diff().std())
                    pass
                #result = self.frame['Q_es'].mean()
                result = [0.] * get_num_species()
                for species in range(get_num_species()):
                    result[species] = autocorr_johannes.get_mean(
                        self.frame['time'].copy().values, self.frame['Q_es'+str(species)].copy().values)
                #cost_per_time=get_cost_per_time(self.get_last_output_file(), self.frame)
                #self.set_result_to_csv(result, self.frame.time.iloc[0],
                #       self.how_long_run(),
                #       cost_per_time=cost_per_time,#TODO
                #       autocorrelation_time=get_autocorrelation_time(self.frame),
                #       time_error=self.get_time_error())
        return result

    def get_result_gplot(self):
        gplot = os.path.join(self.gene_path, "tools/gplot")
        out = os.path.join(self.prob_directory, "out")
        #output = subprocess.Popen(gplot +" "+ out +"/nrg* -m " +str(self.crop_time), \
        #                            stdout=subprocess.PIPE, shell=True)
        # the output of the gplot script contains something like
        # run 0, spec 0: mean(<Q_{es}>)=0.782097 +/- 0.139969 (stddev)
        # run 0, spec 1: mean(<Q_{es}>)=0.118065 +/- 0.0188793 (stddev)
        #output = output.stdout.read().decode("utf-8") #TODO this also returns the gnuplot window
        #print(output)
        output = "run 0, spec 0: mean(<Q_{es}>)=0.782097 +/- 0.139969 (stddev)"
        for spec in range(get_num_species()):
            expression = r'spec\s'+str(spec)+r':\smean\(<Q_\{es\}>\)=(.*)\s\+'
            pattern = re.compile(expression)
            match = pattern.search(output)
            if match is None:
                qes = None
            else:
                qes = match.group(1)
            print("Qes gplot " + str(qes))

    def set_delta(self, delta):
        if (self.gene_path is None):
            #self.set_surplus_to_csv(level_vector, self.prob_prefix, surplus)
            pass
        else:
            self.frame = self.get_frame()
            result = self.get_result()

            self.set_result_to_csv(self.level_vector, result, frame.time.iloc[0],
                                   self.how_long_run(level_vector),
                                   cost_per_time=get_cost_per_time(
                                       self.get_last_output_file(), self.frame),  # TODO
                                   autocorrelation_time=get_autocorrelation_time(
                                       self.frame),
                                   time_error=self.get_time_error(), geometry=self.prob_prefix, surplus=surplus)

    def get_time_error(self):
        if self.gene_path is None:
            time_error = self.get_time_error_csv(self.level_vector)
        else:
            time_error = get_time_error(self.frame)
        for species in range(get_num_species()):
            nwin = autocorr_johannes.get_number_of_autocorrelation_windows(self.frame['time'].copy().values, self.frame['Q_es'+str(species)].copy().values)
            print("NWIN " +str(nwin) + " time/autocorrelation_time " + str(self.how_long_run()/ self.get_autocorrelation_time(species)))
        return time_error

    def how_many_nwin_run(self):
        nwin = [0] * get_num_species()
        for species in range(get_num_species()):
            nwin[species] = autocorr_johannes.get_number_of_autocorrelation_windows(self.frame['time'].copy().values, self.frame['Q_es'+str(species)].copy().values)
            #print("NWIN " +str(nwin) + " time/autocorrelation_time " + str(self.how_long_run()/ self.get_autocorrelation_time(species)))
        return min(nwin)

    def get_autocorrelation_time(self, species):
        if self.gene_path is None:
            autocorrelation_time = self.get_autocorrelation_time_csv(
                species)
        else:
            autocorrelation_time = get_autocorrelation_time(self.frame)[species]
        return autocorrelation_time

    def get_total_cost(self):
        return self.get_cost_per_time() * (self.frame['time'].iloc[-1])

    def get_cost_per_time(self):
        if self.gene_path is None:
            cost = self.get_cost_per_time_csv(self.level_vector)
        else:
            cost = get_cost_per_time(self.get_last_output_file(), self.frame)
        return cost



# if called directly, run tests
if __name__ == "__main__":
    #lv = [5, 6, 6, 5, 3]
    lv = [5,5,5,6,3]
    #prob_prepath="/hppfs/scratch/02/di39qun2/gene3d-flw-simulations/"
    #gene_path="/hppfs/work/pn34mi/di68xux2/myGene3d/"
    qes_data = Qes_data(lv)
    #qes_gplot = qes_data.get_result_gplot()
    qes = qes_data.get_result()
    frame = qes_data.get_frame()
    print(frame)
    print(qes)

    autocorr = get_autocorrelation_time(frame)
    sem = get_time_error(frame, printOutput=True)
    print(autocorr, sem)

    #frame_after_pause = qes_grid.get_frame(lv, '_')
    #autocorr = get_autocorrelation_time(frame_after_pause)
    #sem = get_time_error(frame_after_pause, printOutput=True)
    #print(autocorr, sem)

    print(qes_data.get_last_output_file())
    last_output_file = qes_data.get_last_output_file()
    in_sim_time_per_timestep = get_in_sim_time_per_timestep(frame)
    time_per_timestep = get_time_per_timestep(last_output_file)
    cost = get_cost_per_time(last_output_file, frame)
    #cost_grid = qes_grid.get_total_cost()
    print(cost, in_sim_time_per_timestep, time_per_timestep)
