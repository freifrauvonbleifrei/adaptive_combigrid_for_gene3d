#!/usr/bin/env python3

from write_flux_adaptive_combigrid import get_combiScheme, get_filenames, \
    get_combi_flux, get_qes_csv, get_qes_trapezoidal, printnan, filenames_to_fluxes, \
    get_cost, get_total_cost
from Qes_data import get_num_species
import os
import sys
import subprocess
import pkg_resources

import math
import numpy as np
import pandas as pd
import h5py

from bokeh.layouts import gridplot, column, row
from bokeh.models import Legend, Band, ColumnDataSource, Span, RangeSlider, Range, Panel, Plot, Range1d, Circle, LinearColorMapper, ColorBar, HoverTool
from bokeh.plotting import figure, show, output_notebook
from bokeh.io import export_svgs, export_png
bokeh_version = pkg_resources.get_distribution('bokeh').version


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
for diagnostics_index in [0]: #range(len(diagnostics_df)):
    QoI = diagnostics_df['QoI'][diagnostics_index]
    diagnostics_filename = diagnostics_dir + "profile_"
    diagnostics_filename = diagnostics_dir + diagnostics_df['diagnostics_filename'][diagnostics_index]

    # In[4]:

    prob_prefix = ''
    combiSchemeMode='oldSet.csv'
    #combiSchemeMode='allSet.csv'
    if (len(sys.argv) > 1):
        combiSchemeMode = sys.argv[1]


    combiScheme = get_combiScheme(prob_prefix, combiSchemeMode, dropzeros=True if combiSchemeMode == 'oldSet.csv' else False)

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

    filenames = get_filenames(combiScheme['probname'], diagnostics_filename)
    print(filenames)

    probname = combiScheme['probname'][0]

    fluxes, Xresampled = filenames_to_fluxes(
        filenames, combiScheme['probname'], QoI,
        diagnostics_df['x_axis_name'][diagnostics_index])

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
                    qesCombinedTrap[species] += coefficient * qesProbTrap[species]
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
                    fluxes[probname][species][QoI][q] *= csvRescaleFactor
                    # fluxes[probname][species][QoI][q] *= trapRescaleFactor




    # In[7]:
    combi_flux = get_combi_flux(fluxes, combiScheme, QoI, Xresampled)

    for c in combi_flux:
        c[QoI].clip(lower=0., inplace=True)

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

    svgSuffix = "_" + ("relative" if relativeRescale else "absolute") + "_" + str(rollingAvgNumPoints) + ".svg"

    # # plot reference results (here diagnostics are appended with "2mw")
    # filenames_ref = get_filenames(['2mw'], diagnostics_filename)
    # fluxes_ref, _ = filenames_to_fluxes(filenames_ref, ["ref"], QoI,
    #     diagnostics_df['x_axis_name'][diagnostics_index])
    # plot_ref = get_plot(fluxes_ref["ref"][0], label="reference ions", display_legend=True, width=1400)
    # if get_num_species() > 1:
    #     plot_ref = add_to_plot(plot_ref, fluxes_ref["ref"][1], label="reference electrons", display_legend=True)
    # plot_ref.output_backend = "svg"
    # export_svgs([plot_ref], filename=QoI + "_ref" + svgSuffix)


    # In[10]:

    plot_combi = get_plot(combi_flux[0], label="combi ions" ,display_legend=True, width=1400, color="black")
    if get_num_species() > 1:
        plot_combi = add_to_plot(plot_combi, combi_flux[1], label="combi electrons", display_legend=True, color="black")
    plot_combi.output_backend = "svg"
    export_svgs([plot_combi], filename=QoI + "_" + combiSchemeMode[:-4] + "Combi" + svgSuffix)
    # export_png(plot_combi, filename=QoI + "_" + combiSchemeMode[:-4] + "Combi.png")
    # output_notebook()
    # show(plot_combi)

    # # add "best" component to plot
    # best_component="prob_8_8_5_5_3"
    # filenames_best = get_filenames([best_component], diagnostics_filename)
    # fluxes_best, _ = filenames_to_fluxes(filenames_best, [best_component], QoI,
    #     diagnostics_df['x_axis_name'][diagnostics_index])
    # plot_combi = add_to_plot(plot_combi, fluxes_best[best_component][0], label=best_component, display_legend=True, color="grey")
    # if get_num_species() > 1:
    #     plot_combi = add_to_plot(plot_combi, fluxes_best[best_component][1], label=best_component, display_legend=True, color="grey")

    # # have reference and combi plot in one
    # plot_combi = add_to_plot(plot_combi, fluxes_ref["ref"][0], label="reference ions", display_legend=True, color="turquoise")
    # if get_num_species() > 1:
    #     plot_combi = add_to_plot(plot_combi, fluxes_ref["ref"][1], label="reference electrons", display_legend=True, color="turquoise")
    # export_svgs([plot_combi], filename=QoI + "_" + combiSchemeMode[:-4] + "Combi_and_ref" + svgSuffix)


    # In[11]:


    # fluxes_multiplied_with_coefficient = {}
    # multiplyWithCoefficient = True

    # for i in range(len(fluxes)):
    #     probname = combiScheme.iloc[i]['probname']
    #     fluxes_multiplied_with_coefficient[probname] = {}
    #     for species in range(get_num_species()):
    #         fluxes_multiplied_with_coefficient[probname][species] = fluxes[probname][species].copy()
    #         if multiplyWithCoefficient:
    #             fluxes_multiplied_with_coefficient[probname][species][QoI] = fluxes[probname][species][QoI] * combiScheme.iloc[i]['coefficient']

    # plot_scheme_mult_in_one = get_plot( fluxes_multiplied_with_coefficient[combiScheme.iloc[0]['probname']][0], label=combiScheme.iloc[0]['probname']+" ions")
    # if get_num_species() > 1:
    #     plot_scheme_mult_in_one = add_to_plot(plot_scheme_mult_in_one, fluxes_multiplied_with_coefficient[combiScheme.iloc[0]['probname']][1],label=combiScheme.iloc[0]['probname']+" electrons")
    # for i in range(len(combiScheme[1:])):
    #     plot_scheme_mult_in_one = add_to_plot(plot_scheme_mult_in_one, fluxes_multiplied_with_coefficient[combiScheme.iloc[i+1]['probname']][0], label=combiScheme.iloc[i+1]['probname']+" ions")
    #     if get_num_species() > 1:
    #         plot_scheme_mult_in_one = add_to_plot(plot_scheme_mult_in_one, fluxes_multiplied_with_coefficient[combiScheme.iloc[i+1]['probname']][1], label=combiScheme.iloc[i+1]['probname']+" electrons")


    # # In[12]:


    # # grid = gridplot((plots)+[plot_combi,], ncols=2)
    # #don übereinander plotten, auch mit Koeffizienten multipliziert
    # #DONE alle verfügbaren Gitter einbeziehen
    # #TODO adaption durch l2-unterschied im Flux (oder höhe oder lage des Maximums)
    # plot_scheme_mult_in_one.output_backend = "svg"
    # export_svgs([plot_scheme_mult_in_one], filename=QoI + "_" + combiSchemeMode[:-4] + "AllInOneMultipliedWithCoefficients" + svgSuffix)
    # # export_png(plot_scheme_mult_in_one, filename=QoI + "_" + combiSchemeMode[:-4] + "AllInOneMultipliedWithCoefficients.png")
    # output_notebook()
    # show(plot_scheme_mult_in_one)


    # In[13]:


    plot_scheme_in_one = get_plot(fluxes[combiScheme.iloc[0]['probname']][0],label=combiScheme.iloc[0]['probname']+" ions", display_legend=True)
    if get_num_species() > 1:
        plot_scheme_in_one = add_to_plot(plot_scheme_in_one, fluxes[combiScheme.iloc[0]['probname']][1], label=combiScheme.iloc[0]['probname']+" electrons", display_legend=True)
    for i in range(len(combiScheme[1:])):
        plot_scheme_in_one = add_to_plot(plot_scheme_in_one, fluxes[combiScheme.iloc[i+1]['probname']][0], label=combiScheme.iloc[i+1]['probname']+" ions")
        if get_num_species() > 1:
            plot_scheme_in_one = add_to_plot(plot_scheme_in_one, fluxes[combiScheme.iloc[i+1]['probname']][1], label=combiScheme.iloc[i+1]['probname']+" electrons")
    plot_scheme_in_one = add_to_plot(plot_scheme_in_one, combi_flux[0], label="combi ions", color="black")
    if get_num_species() > 1:
        plot_scheme_in_one = add_to_plot(plot_scheme_in_one, combi_flux[1], label="combi electrons", color="black")

    # # add reference plot
    # plot_scheme_in_one = add_to_plot(plot_scheme_in_one, fluxes_ref["ref"][0], label="reference ions", display_legend=True, color="turquoise")
    # if get_num_species() > 1:
    #     plot_scheme_in_one = add_to_plot(plot_scheme_in_one, fluxes_ref["ref"][1], label="reference electrons", display_legend=True, color="pink")

    plot_scheme_in_one.output_backend = "svg"
    export_svgs([plot_scheme_in_one], filename=QoI + "_" + combiSchemeMode[:-4] + "AllInOne" + svgSuffix)
    # export_png(plot_scheme_in_one, filename=QoI + "_" + combiSchemeMode[:-4] + "AllInOne.png")
    output_notebook()
    show(plot_scheme_in_one)


    # In[14]:
    # plots = [get_plot(fluxes[combiScheme.loc[i]['probname']][s], label=combiScheme.iloc[i]['probname']+str(s), color = 'orange' if combiScheme.iloc[i]['coefficient'] > 0 else 'blue', display_legend=True)for s in [0,1] for i in range(len(fluxes))]
    # output_notebook()
    # show(gridplot([plots]))
