#!/usr/bin/env python3

import os
import sys
import subprocess
import numpy as np
import pandas as pd
import h5py

from Qes_data import get_num_species

diagnostics_dir = './flux_diags/'

diagnostics_df = pd.DataFrame(data={
    'QoI': ["Q_es", "Qes_ky", "Qem_ky", "Ges_ky", "Gem_ky", "T", "T_b", "n", "n_b"],
    'diagnostics_filename': ["flux_profile_", "flux_spectra_Qes_", "flux_spectra_Qem_", "flux_spectra_Ges_", "flux_spectra_Gem_", "profile_", "profile_", "profile_", "profile_"],
    'x_axis_name': ["x_a", "ky", "ky", "ky", "ky", "x_a", "x_a", "x_a", "x_a"]})

relativeRescale = os.environ.get(
    'ADAPTATION_POSTPROCESSING_RELATIVE_COMBINATION', 'True').lower() in ['true', '1']
rollingAvgNumPoints = int(os.environ.get(
    'ADAPTATION_POSTPROCESSING_ROLLING_AVG_NUM_POINTS', '1'))

# (relative rescale currently only works for qes diagnostics)
for diagnostics_index in [0]:  # range(len(diagnostics_df)):
    QoI = diagnostics_df['QoI'][diagnostics_index]
    diagnostics_filename = diagnostics_dir + "profile_"
    diagnostics_filename = diagnostics_dir + \
        diagnostics_df['diagnostics_filename'][diagnostics_index]

    # In[4]:

    prob_prefix = ''
    combiSchemeMode = 'oldSet.csv'
    #combiSchemeMode='allSet.csv'
    startTimeForAverage = ""
    endTimeForAverage = ""

    if (len(sys.argv) > 1):
        combiSchemeMode = sys.argv[1]
    if (len(sys.argv) > 2):
        startTimeForAverage = sys.argv[2]
    if (len(sys.argv) > 3):
        endTimeForAverage = sys.argv[3]

    def get_combiScheme(prefix, mode=combiSchemeMode, dropzeros=True):
        combiScheme = pd.read_csv(prefix+combiSchemeMode, index_col=0)
        if dropzeros:
            combiScheme = combiScheme[combiScheme['coefficient'] != 0].reset_index(
                drop=True)
    #     combiScheme = combiScheme[1:].reset_index(drop=True)
        assert(abs(combiScheme['coefficient'].sum() - 1) < 1e-6)
        return combiScheme

    combiScheme = get_combiScheme(
        prob_prefix, dropzeros=True if combiSchemeMode == 'oldSet.csv' else False)

    def printnan(results):
        for probname in results.index:
            cost_per_time = results['cost_per_time'][probname]
            if cost_per_time == 'cost_per_time':
                cost_per_time = 'nan'
            cost_per_time = float(cost_per_time)
            if (np.isnan(cost_per_time)):
                print("cost_per_time is not known for " + probname)

    def get_cost(results, probname, time_simulated=None):
        """returns the time needed for simulation, the unit is total core-seconds"""
        cost_per_time = results['cost_per_time'][probname]
        if cost_per_time == 'cost_per_time':
            cost_per_time = 'nan'
        cost_per_time = float(cost_per_time)
        if (np.isnan(cost_per_time)):
            print(probname)
        if not time_simulated:
            time_simulated = float(results['qes_to'][probname])
        return cost_per_time * time_simulated

    def get_total_cost(results, combischeme, time_simulated=None):
        costs = [get_cost(results, probname, time_simulated)
                 for probname in combischeme['probname']]
    #     print(costs)
        return np.sum(costs)

    qes_results = pd.read_csv(os.environ.get(
        'ADAPTATION_RESULTS_CSV'), index_col=0)
    printnan(qes_results)

    # try printing the computational costs for running the full scheme (including the zero-coefficient grids)
    try:
        # to test:
        # get_cost(qes_results, "prob_5_5_5_5_4")
        combiSchemeCost = get_combiScheme(prob_prefix, dropzeros=False)
        print("Running the scheme " + combiSchemeMode + " took approximately " +
              str(get_total_cost(qes_results, combiSchemeCost)/3600) + " core-h")
    except Exception as e:
        print("Could not print cost: ")
        print(e)

    # In[6]:

    # extract the flux profiles
    subprocess.run("./extract_flux_profile.sh " +
                   str(startTimeForAverage) + " " + str(endTimeForAverage), shell=True)

    def get_filenames(probname):
        if get_num_species() == 2:
            fnames = [[diagnostics_filename+'ions_' + str(
                x) + '.h5', diagnostics_filename+'electrons_' + str(x) + '.h5'] for x in probname]
        elif get_num_species() == 1:
            fnames = [[diagnostics_filename+'ions_' +
                       str(x) + '.h5'] for x in probname]
        return fnames
    filenames = get_filenames(combiScheme['probname'])
    print(filenames)

    probname = combiScheme['probname'][0]

    def filenames_to_fluxes(flux_filenames, probnames, resample=True):
        fluxes = {}

        for i in range(len(probnames)):
            probname = probnames[i]
            fluxes[probname] = {}
            for species in range(get_num_species()):
                try:
                    with h5py.File(flux_filenames[i][species],  "r") as f:
                        #                 print(flux_filenames[i][species])
                        for item in f.attrs.keys():
                            print(item + ":", f.attrs[item])
                        try:
                            Q_es = f[QoI + '_ions' if species ==
                                     0 else QoI + '_electrons']
                        except KeyError:
                            Q_es = f[QoI]
                        Q_es = np.array(Q_es)
                        x_a = f[diagnostics_df['x_axis_name']
                                [diagnostics_index]]
                        SI_conv = f['SI_conv']
                        d = {QoI: np.array(
                            Q_es) * np.array(SI_conv), 'x_a': np.array(x_a)}
                        fluxes[probname][species] = pd.DataFrame(data=d)
                        # print(fluxes[probname][species][QoI].rolling(window=rollingAvgNumPoints, center=True).sum(), fluxes[probname][species][QoI])
                        fluxes[probname][species][QoI] = fluxes[probname][species][QoI].rolling(
                            window=rollingAvgNumPoints, center=True).sum().div(rollingAvgNumPoints)
                except Exception as e:
                    print("error reading filename: " + flux_filenames[i][species] + " " + str(
                        flux_filenames) + " " + str(i) + " " + str(species))
                    raise e
        if resample:
            Xresampled = set()
            for probname in fluxes:
                x_coords_prob = fluxes[probname][0]['x_a']
                Xresampled = Xresampled | set(x_coords_prob)
            # if in spectral space, remove zero value:
            if (diagnostics_df['x_axis_name'][diagnostics_index] == "ky"):
                Xresampled.discard(0.)
            Xresampled = sorted(Xresampled)

            # cf. https://stackoverflow.com/questions/10464738/interpolation-on-dataframe-in-pandas
            for i in range(len(fluxes)):
                probname = probnames[i]
                for species in range(get_num_species()):
                    modflux = fluxes[probname][species]
                    modflux.set_index('x_a', inplace=True, drop=False)
                    # to avoid interpolate error for newer pandas versions:
                    # cf. https://github.com/pandas-dev/pandas/issues/33956
                    commonIndex = pd.Float64Index(modflux.index.union(
                        Xresampled), dtype=np.float64, name=200.)
                    modflux = modflux.reindex(commonIndex)

                    #Interpolation technique to use. One of:

                    #'linear': Ignore the index and treat the values as equally spaced. This is the only method supported on MultiIndexes.
                    #'time': Works on daily and higher resolution data to interpolate given length of interval.
                    #'index', 'values': use the actual numerical values of the index.
                    #'pad': Fill in NaNs using existing values.
                    #'nearest', 'zero', 'slinear', 'quadratic', 'cubic', 'spline', 'barycentric', 'polynomial': Passed to scipy.interpolate.interp1d. These methods use the numerical values of the index. Both 'polynomial' and 'spline' require that you also specify an order (int), e.g. df.interpolate(method='polynomial', order=5).
                    #'krogh', 'piecewise_polynomial', 'spline', 'pchip', 'akima': Wrappers around the SciPy interpolation methods of similar names. See Notes.
                    #'from_derivatives': Refers to scipy.interpolate.BPoly.from_derivatives which replaces 'piecewise_polynomial' interpolation method in scipy 0.18.
                    # modflux = modflux.interpolate('linear').loc[Xresampled]
                    modflux = modflux.interpolate(
                        'polynomial', order=1).loc[Xresampled]

                    modflux.reset_index(inplace=True, drop=True)
                    # fill values at the ends with zeroes
                    modflux.fillna(0., inplace=True)
                    fluxes[probname][species] = modflux

            return fluxes, Xresampled
        else:
            return fluxes

    fluxes, Xresampled = filenames_to_fluxes(
        filenames, combiScheme['probname'])

    def get_qes_trapezoidal(fluxes, probname):
        qes = [0.]*get_num_species()
        for species in range(get_num_species()):
            # print(fluxes[probname][species])
            qes[species] = np.trapz(
                fluxes[probname][species]["Q_es"], x=fluxes[probname][species]["x_a"], axis=0)
        return qes

    def get_qes_csv(results, probname):
        qes0 = results['qes0'][probname]
        if qes0 == 'qes0':
            return None
        qes0 = float(qes0)
        if get_num_species() == 2:
            qes = [qes0, float(results['qes1'][probname])]
        else:
            qes = [qes0]
        return qes

    if relativeRescale:
        # get combined average QoI
        qesCombined = [0.]*get_num_species()
        qesCombinedTrap = [0.]*get_num_species()
        if QoI is "Q_es":  # or QoI is "Qes_ky":
            # read from qes_results.csv
            for component in combiScheme.itertuples(index=False):
                probname = component.probname
                coefficient = float(component.coefficient)
                qesProb = get_qes_csv(qes_results, probname)
                qesProbTrap = get_qes_trapezoidal(fluxes, probname)
                print("Qes " + probname + " ", qesProb, qesProbTrap)
                for species in range(get_num_species()):
                    qesCombined[species] += coefficient * qesProb[species]
                    qesCombinedTrap[species] += coefficient * \
                        qesProbTrap[species]
            print("qes combined: " + str(qesCombined) + str(qesCombinedTrap))
        else:
            # get from curves by trapezoidal rule
            raise NotImplementedError
        # rescale all the component fluxes
        for probname in combiScheme['probname']:
            for species in range(get_num_species()):
                print("rescaling ", probname, species)
                csvRescaleFactor = qesCombined[species] / \
                    get_qes_csv(qes_results, probname)[species]
                trapRescaleFactor = qesCombinedTrap[species] / \
                    get_qes_trapezoidal(fluxes, probname)[species]
                if abs(csvRescaleFactor - trapRescaleFactor) / csvRescaleFactor > 0.1:
                    print("different rescaling relations! ",
                          csvRescaleFactor, trapRescaleFactor)
                for q in range(len(fluxes[probname][species][QoI])):
                    # fluxes[probname][species][QoI][q] *= csvRescaleFactor
                    fluxes[probname][species][QoI][q] *= trapRescaleFactor

    # In[7]:

    def get_combi_flux(fluxes, combiScheme):
        # print(len(fluxes))
        combi_flux = []
        for species in range(get_num_species()):
            combi_flux_species = [fluxes[combiScheme.loc[x]['probname']][species]
                                  [QoI] * combiScheme.loc[x]['coefficient'] for x in range(len(fluxes))]
            combi_flux_species = pd.DataFrame(data=combi_flux_species).sum()
            combi_flux_species = pd.DataFrame(
                data={QoI: combi_flux_species, 'x_a': Xresampled})
            combi_flux.append(combi_flux_species)
        return combi_flux
    combi_flux = get_combi_flux(fluxes, combiScheme)

    # clip away negative values
    for c in combi_flux:
        c[QoI].clip(lower=0., inplace=True)

    # write out combined flux as h5 and .txt data
    combi_filename_begin = "./" + QoI + "_" + combiSchemeMode[:-4] + "_Combi_"
    with h5py.File(combi_filename_begin + "ions.h5",  "w") as f:
        f.create_dataset('x_a', data=combi_flux[0]['x_a'])
        f.create_dataset(QoI, data=combi_flux[0][QoI])
    np.savetxt(combi_filename_begin + "ions.txt", combi_flux[0])
    if get_num_species() > 1:
        with h5py.File(combi_filename_begin + "electrons.h5",  "w") as f:
            f.create_dataset('x_a', data=combi_flux[1]['x_a'])
            f.create_dataset(QoI, data=combi_flux[1][QoI])
        np.savetxt(combi_filename_begin + "electrons.txt", combi_flux[1])