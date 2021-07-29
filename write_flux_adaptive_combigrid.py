#!/usr/bin/env python3

import os
import sys
import subprocess
import numpy as np
import pandas as pd
import h5py
import scipy.interpolate

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


def interpolate_1d_qty(xGene, xTango, qtyGene):
    """Interpolate a 1D quantity from Gene radial grid to Tango radial grid.

    Inputs:
      xGene             rho_tor grid as provided by Gene (1D array)
      xTango            rho_tor grid for Tango on which to interpolate the quantity (1D array)
      qtyGene         quantity evaluated on Genegrid points xGene (1D array)

    Outputs:
      qtyTango          quantity interpolated onto xTango grid points (1D array)
    """
    interpolator = scipy.interpolate.InterpolatedUnivariateSpline(
        xGene, qtyGene, k=1)
    qtyTango = interpolator(xTango)
    return qtyTango


def get_combiScheme(prefix, mode, dropzeros=True):
    combiScheme = pd.read_csv(prefix+mode, index_col=0)
    if dropzeros:
        combiScheme = combiScheme[combiScheme['coefficient'] != 0].reset_index(
            drop=True)
#   combiScheme = combiScheme[1:].reset_index(drop=True)
    assert(abs(combiScheme['coefficient'].sum() - 1) < 1e-6)
    return combiScheme


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


def get_filenames(probname, diagnostics_filename):
    if get_num_species() == 2:
        fnames = [[diagnostics_filename+'ions_' + str(
            x) + '.h5', diagnostics_filename+'electrons_' + str(x) + '.h5'] for x in probname]
    elif get_num_species() == 1:
        fnames = [[diagnostics_filename+'ions_' +
                   str(x) + '.h5'] for x in probname]
    return fnames


def csv_to_flux(flux_filename, column_names):
    flux = pd.read_table(flux_filename, sep="\s+",
                         header=None, names=column_names)
    return flux


def filenames_to_fluxes(flux_filenames, probnames, QoI, x_name, resample=True):
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
                    x_a = f[x_name]
                    # make sure that SI_conv is an array qty
                    try:
                        SI_conv = f['SIA_conv']
                    except KeyError as ke:
                        # fall back to 'SI_conf' if SIA_conf does not exist
                        SI_conv = f['SI_conv']
                    SI_conv = np.array(SI_conv)
                    assert(len(SI_conv) > 1)
                    assert(len(SI_conv) == len(Q_es))
                    d = {QoI: np.array(Q_es), 'x_a': np.array(
                        x_a), 'SI_conv': np.array(SI_conv)}
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
        if (x_name == "ky"):
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
                # Interpolation technique to use. One of:

                # 'linear': Ignore the index and treat the values as equally spaced. This is the only method supported on MultiIndexes.
                # 'time': Works on daily and higher resolution data to interpolate given length of interval.
                # 'index', 'values': use the actual numerical values of the index.
                # 'pad': Fill in NaNs using existing values.
                # 'nearest', 'zero', 'slinear', 'quadratic', 'cubic', 'spline', 'barycentric', 'polynomial': Passed to scipy.interpolate.interp1d. These methods use the numerical values of the index. Both 'polynomial' and 'spline' require that you also specify an order (int), e.g. df.interpolate(method='polynomial', order=5).
                # 'krogh', 'piecewise_polynomial', 'spline', 'pchip', 'akima': Wrappers around the SciPy interpolation methods of similar names. See Notes.
                # 'from_derivatives': Refers to scipy.interpolate.BPoly.from_derivatives which replaces 'piecewise_polynomial' interpolation method in scipy 0.18.
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


def get_qes_trapezoidal(fluxes, probname, withSIconv=True):
    qes = [0.]*get_num_species()
    for species in range(get_num_species()):
        # print(fluxes[probname][species])
        if withSIconv:
            flux = fluxes[probname][species]["Q_es"] * \
                fluxes[probname][species]["SI_conv"]
        else:
            flux = fluxes[probname][species]["Q_es"]
        qes[species] = np.trapz(
            flux, x=fluxes[probname][species]["x_a"], axis=0)
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


def get_combi_flux(fluxes, combiScheme, QoI, x):
    # print(len(fluxes))
    combi_flux = []
    for species in range(get_num_species()):
        combi_flux_species_list = [fluxes[combiScheme.loc[x]['probname']][species]
                                   [QoI] * combiScheme.loc[x]['coefficient'] for x in range(len(fluxes))]
        combi_si_conv_list = [fluxes[combiScheme.loc[x]['probname']][species]
                              ["SI_conv"] * combiScheme.loc[x]['coefficient'] for x in range(len(fluxes))]
        combi_flux_species = pd.DataFrame(
            data=combi_flux_species_list).sum()
        combi_si_conv = pd.DataFrame(
            data=combi_si_conv_list).sum()
        combi_flux_species = pd.DataFrame(
            data={QoI: combi_flux_species, 'x_a': x, 'SI_conv': combi_si_conv})
        combi_flux.append(combi_flux_species)
    return combi_flux


def write_flux(combiSchemeMode, startTimeForAverage, endTimeForAverage):
    # (relative rescale currently only works for qes diagnostics)
    for diagnostics_index in [0]:  # range(len(diagnostics_df)):
        QoI = diagnostics_df['QoI'][diagnostics_index]
        diagnostics_filename = diagnostics_dir + "profile_"
        diagnostics_filename = diagnostics_dir + \
            diagnostics_df['diagnostics_filename'][diagnostics_index]

        # In[4]:

        prob_prefix = ''

        combiScheme = get_combiScheme(
            prob_prefix, combiSchemeMode, dropzeros=True if combiSchemeMode == 'oldSet.csv' else False)

        csv_name = os.environ.get('ADAPTATION_RESULTS_CSV')
        qes_results = pd.read_csv(csv_name, index_col=0)
        printnan(qes_results)

        # try printing the computational costs for running the full scheme (including the zero-coefficient grids)
        try:
            # to test:
            # get_cost(qes_results, "prob_5_5_5_5_4")
            combiSchemeCost = get_combiScheme(
                prob_prefix, combiSchemeMode, dropzeros=False)
            print("Running the scheme " + combiSchemeMode + " took approximately " +
                  str(get_total_cost(qes_results, combiSchemeCost)/3600) + " core-h")
        except Exception as e:
            print("Could not print cost: ")
            print(e)

        # In[6]:

        # extract the flux profiles
        subprocess.run("./extract_flux_profile.sh " +
                       str(startTimeForAverage) + " " + str(endTimeForAverage), shell=True)

        filenames = get_filenames(
            combiScheme['probname'], diagnostics_filename)
        print(filenames)

        probname = combiScheme['probname'][0]

        fluxes, Xresampled = filenames_to_fluxes(
            filenames, combiScheme['probname'], QoI, diagnostics_df['x_axis_name']
            [diagnostics_index])

        if relativeRescale:
            powerNormalization = True
            # get combined average QoI
            qesCombined = [0.]*get_num_species()
            qesCombinedTrap = [0.]*get_num_species()
            if QoI is "Q_es":  # or QoI is "Qes_ky":
                # read from qes_results.csv
                for component in combiScheme.itertuples(index=False):
                    probname = component.probname
                    coefficient = float(component.coefficient)
                    qesProb = get_qes_csv(qes_results, probname)
                    qesProbTrap = get_qes_trapezoidal(
                        fluxes, probname, withSIconv=powerNormalization)
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
                        get_qes_trapezoidal(
                            fluxes, probname, withSIconv=powerNormalization)[species]
                    if abs(csvRescaleFactor - trapRescaleFactor) / csvRescaleFactor > 0.1:
                        print("different rescaling relations! ",
                              csvRescaleFactor, trapRescaleFactor)
                    for q in range(len(fluxes[probname][species][QoI])):
                        # fluxes[probname][species][QoI][q] *= csvRescaleFactor
                        fluxes[probname][species][QoI][q] *= trapRescaleFactor

        # In[7]:
        combi_flux = get_combi_flux(fluxes, combiScheme, QoI, Xresampled)

        # clip away negative values
        for c in combi_flux:
            c[QoI].clip(lower=0., inplace=True)

        # write out combined flux as h5 and .txt data
        # combi_filename_begin = "./" + QoI + \
        #    "_" + combiSchemeMode[:-4] + "_Combi_"
        # with h5py.File(combi_filename_begin + "ions.h5",  "w") as f:
        #    f.create_dataset('x_a', data=combi_flux[0]['x_a'])
        #    f.create_dataset(QoI, data=combi_flux[0][QoI])
        # #Do the same for electrons
        # if get_num_species() > 1:
        #     with h5py.File(combi_filename_begin + "electrons.h5",  "w") as f:
        #         f.create_dataset('x_a', data=combi_flux[1]['x_a'])
        #         f.create_dataset(QoI, data=combi_flux[1][QoI])
        #     np.savetxt(combi_filename_begin + "electrons.txt", combi_flux[1])

        # ABN
        with h5py.File('./flux_diags/flux_profile_ions_prob_5_5_5_5_3.h5', 'r') as hf:
            try:
                SIA_conv_file = np.array(hf.get('/SIA_conv')).T
                SI_conv = np.array(hf.get('/SI_conv')).T
            except KeyError as ke:
                # fall back to 'SI_conf' if SIA_conf does not exist
                print(hf.get('/SI_conv'))
                SIA_conv_file = np.array(hf.get('/SI_conv')).T
                SI_conv = np.array(hf.get('/SI_noA_conv')).T
            assert(not isinstance(SI_conv, list))
            x_a_file = np.array(hf.get('/x_a')).T

        # make sure ABN's way of getting SI_conv is the same as re-using it
        # print(SIA_conv_file, combi_flux[0]["SI_conv"])
        SIA_conv = interpolate_1d_qty(x_a_file, np.array(
            combi_flux[0]['x_a']), SIA_conv_file)
        print(SIA_conv, len(SIA_conv), combi_flux[0]["SI_conv"])

        np.savetxt("heat_flux_ions.txt", np.c_[
                   combi_flux[0]['x_a'], combi_flux[0][QoI]])
        np.savetxt("heat_flux_ions_" + str(startTimeForAverage) + "_" + str(endTimeForAverage) + ".txt",
                   np.c_[combi_flux[0]['x_a'], combi_flux[0][QoI], combi_flux[0][QoI]*combi_flux[0]['SI_conv']])
        np.savetxt("heat_flux_seed_ions", np.c_[
                   combi_flux[0]['x_a'], combi_flux[0][QoI]*SI_conv])


if __name__ == "__main__":
    combiSchemeMode = 'oldSet.csv'
    # combiSchemeMode='allSet.csv'
    startTimeForAverage = ""
    endTimeForAverage = ""

    if (len(sys.argv) > 1):
        combiSchemeMode = sys.argv[1]
    if (len(sys.argv) > 2):
        startTimeForAverage = sys.argv[2]
    if (len(sys.argv) > 3):
        endTimeForAverage = sys.argv[3]
    write_flux(combiSchemeMode, startTimeForAverage, endTimeForAverage)
