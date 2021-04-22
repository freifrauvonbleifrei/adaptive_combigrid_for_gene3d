#!/usr/bin/env python3

import os
import sys
import subprocess
import pkg_resources

import math
import numpy as np
import pandas as pd

from bokeh.layouts import gridplot, column, row
from bokeh.models import Legend, Band, ColumnDataSource, Span, RangeSlider, Range, Panel, Plot, Range1d, Circle, LinearColorMapper, ColorBar, HoverTool
from bokeh.plotting import figure, show, output_notebook
from bokeh.io import export_svgs, export_png
bokeh_version = pkg_resources.get_distribution('bokeh').version

from Qes_data import get_num_species

import h5py


# In[2]:


# #have wider cells in notebook
# from IPython.core.display import display, HTML
# display(HTML("<style>.container { width:100% !important; }</style>"))


# In[3]:


diagnostics_dir = './flux_diags/'

diagnostics_df = pd.DataFrame(data={
    'QoI': ["Q_es","Qes_ky","Qem_ky","Ges_ky","Gem_ky", "T", "T_b", "n", "n_b"],
    'diagnostics_filename': ["flux_profile_", "flux_spectra_Qes_", "flux_spectra_Qem_", "flux_spectra_Ges_", "flux_spectra_Gem_", "profile_", "profile_", "profile_", "profile_"],
    'x_axis_name' : ["x_a", "ky", "ky", "ky", "ky", "x_a", "x_a", "x_a", "x_a"]})

relativeRescale = os.environ.get('ADAPTATION_POSTPROCESSING_RELATIVE_COMBINATION', 'True').lower() in ['true', '1']
rollingAvgNumPoints = int(os.environ.get('ADAPTATION_POSTPROCESSING_ROLLING_AVG_NUM_POINTS', '1'))


# (relative rescale currently only works for qes diagnostics)
for diagnostics_index in [0,1]: #range(len(diagnostics_df)):
    QoI = diagnostics_df['QoI'][diagnostics_index]
    diagnostics_filename = diagnostics_dir + "profile_"
    diagnostics_filename = diagnostics_dir + diagnostics_df['diagnostics_filename'][diagnostics_index]

    # In[4]:

    prob_prefix = ''
    combiSchemeMode='oldSet.csv'
    #combiSchemeMode='allSet.csv'
    if (len(sys.argv) > 1):
        combiSchemeMode = sys.argv[1]

    def get_combiScheme(prefix, mode=combiSchemeMode, dropzeros=True):
        combiScheme = pd.read_csv(prefix+combiSchemeMode, index_col=0)
        if dropzeros:
            combiScheme = combiScheme[combiScheme['coefficient'] != 0].reset_index(drop=True)
    #     combiScheme = combiScheme[1:].reset_index(drop=True)
        assert(abs( combiScheme['coefficient'].sum() - 1) < 1e-6)
        return combiScheme

    combiScheme = get_combiScheme(prob_prefix, dropzeros=True)


    # In[5]:


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
        costs = [get_cost(results, probname, time_simulated) for probname in combischeme['probname']]
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
    subprocess.run("./extract_flux_profile.sh", shell=True)

    def get_filenames(probname):
        if get_num_species() == 2:
            fnames = [[diagnostics_filename+'ions_'+ str(x) +'.h5', diagnostics_filename+'electrons_'+ str(x) +'.h5'] for x in probname]
        elif get_num_species() == 1:
            fnames = [[diagnostics_filename+'ions_'+ str(x) +'.h5'] for x in probname]
        return fnames
    filenames = get_filenames(combiScheme['probname'])
    print(filenames)

    probname = combiScheme['probname'][0]

    def filenames_to_fluxes(flux_filenames, probnames, resample=True):
        fluxes = {}

        for i in range(len(probnames)):
            probname = probnames[i]
            fluxes[probname]={}
            for species in range(get_num_species()):
                try:
                    with h5py.File(flux_filenames[i][species],  "r") as f:
        #                 print(flux_filenames[i][species])
                        for item in f.attrs.keys():
                            print(item + ":", f.attrs[item])
                        try:
                            Q_es = f[QoI + '_ions'  if species == 0 else QoI + '_electrons']
                        except KeyError:
                            Q_es = f[QoI]
                        Q_es = np.array(Q_es)
                        x_a = f[diagnostics_df['x_axis_name'][diagnostics_index]]
                        d = {QoI: np.array(Q_es), 'x_a': np.array(x_a)}
                        fluxes[probname][species] = pd.DataFrame(data=d)
                        # print(fluxes[probname][species][QoI].rolling(window=rollingAvgNumPoints, center=True).sum(), fluxes[probname][species][QoI])
                        fluxes[probname][species][QoI] = fluxes[probname][species][QoI].rolling(window=rollingAvgNumPoints, center=True).sum().div(rollingAvgNumPoints)
                except Exception as e:
                    print("error reading filename: " + flux_filenames[i][species] + " " + str(flux_filenames) + " " + str(i) + " " + str(species))
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

    fluxes, Xresampled = filenames_to_fluxes(filenames, combiScheme['probname'])

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
            # print("Qes " + probname + " ", qesProb, qesProbTrap)
            for species in range(get_num_species()):
                qesCombined[species] += coefficient * qesProb[species]
                qesCombinedTrap[species] += coefficient * qesProbTrap[species]
        print("qes combined: " + str(qesCombined) + str(qesCombinedTrap))
    else:
        # get from curves by trapezoidal rule
        raise NotImplementedError

    if relativeRescale:
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
                    fluxes[probname][species][QoI][q] *= csvRescaleFactor




    # In[7]:


    def get_combi_flux(fluxes, combiScheme):
        # print(len(fluxes))
        combi_flux= []
        for species in range(get_num_species()):
            combi_flux_species = [fluxes[combiScheme.loc[x]['probname']][species][QoI] * combiScheme.loc[x]['coefficient'] for x in range(len(fluxes))]
            combi_flux_species = pd.DataFrame(data=combi_flux_species).sum()
            combi_flux_species = pd.DataFrame(data={QoI: combi_flux_species, 'x_a': Xresampled})
            combi_flux.append(combi_flux_species)
        return combi_flux
    combi_flux = get_combi_flux(fluxes, combiScheme)


    # In[8]:



    from bokeh.palettes import Category10_10 as palette
    import itertools

    # create a color iterator
    global colors
    colors = itertools.cycle(palette)

    # get profile plot
    def get_figure(x, width=None, height=None):
        x_range = Range1d(start=x.min(), end=x.max())
        if (diagnostics_df['x_axis_name'][diagnostics_index] == "ky"):
            return figure(plot_width=width if width else 1000, plot_height=height if height else 500, x_range=x_range, y_axis_type="log")#, x_axis_type="log") # why won't this work???
        else:
            return figure(plot_width=width if width else 1000, plot_height=height if height else 500, x_range=x_range)

    def add_to_plot(plot, df, label=None, color=None, display_legend=False):
        if not color:
            color = next(colors)
        source = ColumnDataSource(df)
        if bokeh_version > "2.":
            q = plot.line("x_a", QoI, color=color, source=source, legend_label=label)
        else:
            q = plot.line("x_a", QoI, color=color, source=source)
        plot.add_tools(HoverTool(renderers=[q], tooltips=[("",label)]))
        if display_legend and bokeh_version < "2.":
            legend = Legend(items=[
                (label,   [q]),
            ], location=(0, -30))
            plot.add_layout(legend, 'right')
        return plot

    def get_plot(df, label=None, color=None, display_legend=False, width=None, height=None):
        global colors 
        colors = itertools.cycle(palette)
        plot = get_figure(df["x_a"], width=width, height=height)
        plot.xaxis.axis_label = 'radial coordinate ρ_tor' if (diagnostics_df['x_axis_name'][diagnostics_index] == "x_a") else 'spectral coordinate k_y'
        plot.yaxis.axis_label = QoI #'ion heat flux Q_es'
        return add_to_plot(plot, df, label, color, display_legend)


    # $\rho_{tor}$, $Q_j^{es}$

    # In[9]:

    ## plot reference results (here diagnostics are appended with "flw")
    # filenames_ref = get_filenames(['flw'])
    # fluxes_ref, _ = filenames_to_fluxes(filenames_ref, ["ref"])
    # plot_ref = get_plot(fluxes_ref["ref"][0], label="reference ions", display_legend=True, width=1400)
    # if get_num_species() > 1:
    #     plot_ref = add_to_plot(plot_ref, fluxes_ref["ref"][1], label="reference electrons", display_legend=True)
    # plot_ref.output_backend = "svg"
    # export_svgs([plot_ref], filename=QoI + "_ref.svg")
    # show(plot_ref)


    # In[10]:


    plot_combi = get_plot(combi_flux[0], label="combi ions" ,display_legend=True, width=1400)
    if get_num_species() > 1:
        plot_combi = add_to_plot(plot_combi, combi_flux[1], label="combi electrons", display_legend=True)
    plot_combi.output_backend = "svg"
    export_svgs([plot_combi], filename=QoI + "_" + combiSchemeMode[:-4] + "Combi.svg")
    # export_png(plot_combi, filename=QoI + "_" + combiSchemeMode[:-4] + "Combi.png")
    output_notebook()
    show(plot_combi)

    ## have reference and combi plot in one
    # plot_combi = add_to_plot(plot_combi, fluxes_ref["ref"][0], label="reference ions", display_legend=True)
    # if get_num_species() > 1:
    #     plot_combi = add_to_plot(plot_combi, fluxes_ref["ref"][1], label="reference electrons", display_legend=True)
    # export_svgs([plot_combi], filename=QoI + "_" + combiSchemeMode[:-4] + "Combi_and_ref.svg")


    # In[11]:


    fluxes_multiplied_with_coefficient = {}
    multiplyWithCoefficient = True

    for i in range(len(fluxes)):
        probname = combiScheme.iloc[i]['probname']
        fluxes_multiplied_with_coefficient[probname] = {}
        for species in range(get_num_species()):
            fluxes_multiplied_with_coefficient[probname][species] = fluxes[probname][species].copy()
            if multiplyWithCoefficient:
                fluxes_multiplied_with_coefficient[probname][species][QoI] = fluxes[probname][species][QoI] * combiScheme.iloc[i]['coefficient']

    plot_scheme_mult_in_one = get_plot( fluxes_multiplied_with_coefficient[combiScheme.iloc[0]['probname']][0], label=combiScheme.iloc[0]['probname']+" ions")
    if get_num_species() > 1:
        plot_scheme_mult_in_one = add_to_plot(plot_scheme_mult_in_one, fluxes_multiplied_with_coefficient[combiScheme.iloc[0]['probname']][1],label=combiScheme.iloc[0]['probname']+" electrons")
    for i in range(len(combiScheme[1:])):
        plot_scheme_mult_in_one = add_to_plot(plot_scheme_mult_in_one, fluxes_multiplied_with_coefficient[combiScheme.iloc[i+1]['probname']][0], label=combiScheme.iloc[i+1]['probname']+" ions")
        if get_num_species() > 1:
            plot_scheme_mult_in_one = add_to_plot(plot_scheme_mult_in_one, fluxes_multiplied_with_coefficient[combiScheme.iloc[i+1]['probname']][1], label=combiScheme.iloc[i+1]['probname']+" electrons")


    # In[12]:


    # grid = gridplot((plots)+[plot_combi,], ncols=2)
    #don übereinander plotten, auch mit Koeffizienten multipliziert
    #DONE alle verfügbaren Gitter einbeziehen
    #TODO adaption durch l2-unterschied im Flux (oder höhe oder lage des Maximums)
    plot_scheme_mult_in_one.output_backend = "svg"
    export_svgs([plot_scheme_mult_in_one], filename=QoI + "_" + combiSchemeMode[:-4] + "AllInOneMultipliedWithCoefficients.svg")
    # export_png(plot_scheme_mult_in_one, filename=QoI + "_" + combiSchemeMode[:-4] + "AllInOneMultipliedWithCoefficients.png")
    output_notebook()
    show(plot_scheme_mult_in_one)


    # In[13]:


    plot_scheme_in_one = get_plot(fluxes[combiScheme.iloc[0]['probname']][0],label=combiScheme.iloc[0]['probname']+" ions", display_legend=True)
    if get_num_species() > 1:
        plot_scheme_in_one = add_to_plot(plot_scheme_in_one, fluxes[combiScheme.iloc[0]['probname']][1], label=combiScheme.iloc[0]['probname']+" electrons", display_legend=True)
    for i in range(len(combiScheme[1:])):
        plot_scheme_in_one = add_to_plot(plot_scheme_in_one, fluxes[combiScheme.iloc[i+1]['probname']][0], label=combiScheme.iloc[i+1]['probname']+" ions")
        if get_num_species() > 1:
            plot_scheme_in_one = add_to_plot(plot_scheme_in_one, fluxes[combiScheme.iloc[i+1]['probname']][1], label=combiScheme.iloc[i+1]['probname']+" electrons")
    plot_scheme_in_one.output_backend = "svg"
    export_svgs([plot_scheme_in_one], filename=QoI + "_" + combiSchemeMode[:-4] + "AllInOne.svg")
    # export_png(plot_scheme_in_one, filename=QoI + "_" + combiSchemeMode[:-4] + "AllInOne.png")
    output_notebook()
    show(plot_scheme_in_one)


    # In[14]:
    # plots = [get_plot(fluxes[combiScheme.loc[i]['probname']][s], label=combiScheme.iloc[i]['probname']+str(s), color = 'orange' if combiScheme.iloc[i]['coefficient'] > 0 else 'blue', display_legend=True)for s in [0,1] for i in range(len(fluxes))]
    # output_notebook()
    # show(gridplot([plots]))
