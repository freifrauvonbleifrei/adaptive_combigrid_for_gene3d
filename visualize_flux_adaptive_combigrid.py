#!/usr/bin/env python3

import os
import sys
import subprocess

import math
import numpy as np
import pandas as pd

from bokeh.layouts import gridplot, column, row
from bokeh.models import Legend, Band, ColumnDataSource, Span, RangeSlider, Range, Panel, Plot, Range1d, Circle, LinearColorMapper, ColorBar, HoverTool
from bokeh.plotting import figure, show, output_notebook
from bokeh.io import export_svgs

from Qes_data import get_num_species

import h5py


# In[2]:


# #have wider cells in notebook
# from IPython.core.display import display, HTML
# display(HTML("<style>.container { width:100% !important; }</style>"))


# In[3]:


QoI = "Q_es"


# In[4]:

prob_prefix = ''
combiSchemeMode='oldSet.csv'
combiSchemeMode='allSet.csv'
if (len(sys.argv) > 1):
    combiSchemeMode = sys.argv[1]

def get_combiScheme(prefix, mode=combiSchemeMode, dropzeros=True):
    combiScheme = pd.read_csv(prefix+combiSchemeMode, index_col=0)
    if dropzeros:
        combiScheme = combiScheme[combiScheme['coefficient'] != 0].reset_index(drop=True)
#     combiScheme = combiScheme[1:].reset_index(drop=True)
    assert(abs( combiScheme['coefficient'].sum() - 1) < 1e-6)
    return combiScheme

combiScheme = get_combiScheme(prob_prefix, dropzeros=False)

lx_range = combiScheme['l_0'].copy().drop_duplicates().tolist()
spacings = [2**lx -1 for lx in lx_range]
denominator_spacings = 2047
combiScheme
# filenames
# lx_range


# In[5]:


def printnan(results):
    for probname in results.index:
        cost_per_time = results['cost_per_time'][probname]
        if (np.isnan(cost_per_time)):
            print("cost_per_time is not known for " + probname)

def get_cost(results, probname, time_simulated=None):
    """returns the time needed for simulation, the unit is total core-seconds"""
    cost_per_time = results['cost_per_time'][probname]
    if (np.isnan(cost_per_time)):
        print(probname)
    if not time_simulated:
        time_simulated = results['qes_to'][probname]
    return cost_per_time * time_simulated

def get_total_cost(results, combischeme, time_simulated=None):
    costs = [get_cost(results, probname, time_simulated) for probname in combischeme['probname']]
#     print(costs)
    return np.sum(costs)


qes_results = pd.read_csv(os.environ.get('ADAPTATION_RESULTS_CSV'), index_col=0)
printnan(qes_results)

get_cost(qes_results, "prob_5_5_5_5_4")

combiSchemeCost = get_combiScheme(prob_prefix, dropzeros=False)
print("Running the scheme " + combiSchemeMode + " took approximately " + str(get_total_cost(qes_results, combiSchemeCost)/3600) + " core-h")


# In[6]:

# extract the flux profiles
subprocess.run("./extract_flux_profiles.sh", shell=True)

def get_filenames(probname):
    fnames = [('./flux_diags/flux_profile_ions_'+ str(x) +'.h5', './flux_diags/flux_profile_electrons_'+ str(x) +'.h5') for x in probname]
    return fnames
filenames = get_filenames(combiScheme['probname'])

probname = combiScheme['probname'][0]

def filenames_to_fluxes(flux_filenames, resample=True):
    fluxes = {}

    for i in range(len(combiScheme['probname'])):
        probname = combiScheme['probname'][i]
        fluxes[probname]={}
        for species in [0,1]:
            with h5py.File(flux_filenames[i][species],  "r") as f:
#                 print(flux_filenames[i][species])
                for item in f.attrs.keys():
                    print(item + ":", f.attrs[item])    
                Q_es = f['Q_es_ions'  if species == 0 else 'Q_es_electrons'] 
                x_a = f['x_a']
                d = {'Q_es': np.array(Q_es), 'x_a': np.array(x_a)}
                fluxes[probname][species] = pd.DataFrame(data=d)
    if resample:
        x_a_range = [fluxes[probname][0]['x_a'].min(), fluxes[probname][0]['x_a'].max()]
        Xresampled = np.linspace(x_a_range[0],x_a_range[1],denominator_spacings)

        # cf. https://stackoverflow.com/questions/10464738/interpolation-on-dataframe-in-pandas
        for i in range(len(fluxes)):
            probname = combiScheme['probname'][i]
            for species in [0,1]:
                modflux = fluxes[probname][species]
                modflux.set_index('x_a',inplace =True, drop=False)

                #Interpolation technique to use. One of:

                #'linear': Ignore the index and treat the values as equally spaced. This is the only method supported on MultiIndexes.
                #'time': Works on daily and higher resolution data to interpolate given length of interval.
                #'index', 'values': use the actual numerical values of the index.
                #'pad': Fill in NaNs using existing values.
                #'nearest', 'zero', 'slinear', 'quadratic', 'cubic', 'spline', 'barycentric', 'polynomial': Passed to scipy.interpolate.interp1d. These methods use the numerical values of the index. Both 'polynomial' and 'spline' require that you also specify an order (int), e.g. df.interpolate(method='polynomial', order=5).
                #'krogh', 'piecewise_polynomial', 'spline', 'pchip', 'akima': Wrappers around the SciPy interpolation methods of similar names. See Notes.
                #'from_derivatives': Refers to scipy.interpolate.BPoly.from_derivatives which replaces 'piecewise_polynomial' interpolation method in scipy 0.18.

        #         modflux = modflux.reindex(modflux.index.union(Xresampled)).interpolate('linear').loc[Xresampled]
                modflux = modflux.reindex(modflux.index.union(Xresampled)).interpolate('polynomial', order=5).loc[Xresampled]
                modflux.reset_index(inplace=True, drop=True)
                fluxes[probname][species] = modflux
        return fluxes, Xresampled
    else:
        return fluxes
    
fluxes, Xresampled = filenames_to_fluxes(filenames)
# fluxes = filenames_to_fluxes(filenames, False)

# fluxes
# len(Xresampled), len(fluxes[probname][0]['Q_es'])
# fluxes[probname][0], fluxes[probname][1]


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
    return figure(plot_width=width if width else 1000, plot_height=height if height else 500, x_range=x_range)

def get_plot_xy(x, y, color=None):
    if not color:
        color="navy"
    plot = get_figure(x)
    q = plot.line(x, y, color=color, alpha=0.5)
    return plot

def add_to_plot(plot, df, label=None, color=None, display_legend=False):
    if not color:
        color = next(colors)
    source = ColumnDataSource(df)
    q = plot.line("x_a", QoI, color=color, source=source)
    plot.add_tools(HoverTool(renderers=[q], tooltips=[("",label)]))
    if display_legend:
        legend = Legend(items=[
            (label,   [q]),
        ], location=(0, -30))
        plot.add_layout(legend, 'right')
    return plot

def get_plot(df, label=None, color=None, display_legend=False, width=None, height=None):
    global colors 
    colors = itertools.cycle(palette)
    plot = get_figure(df["x_a"], width=width, height=height)
    plot.xaxis.axis_label = 'radial coordinate ρ_tor'
    plot.yaxis.axis_label = 'ion heat flux Q_es'
    return add_to_plot(plot, df, label, color, display_legend)


# $\rho_{tor}$, $Q_j^{es}$

# In[9]:


# fluxes_alejandro = []
# prefixes_alejandro = ['aim', 'dbm', 'eim00', 'ftm', 'kjm', 'eim00_hr']
# for p in prefixes_alejandro:
#     with h5py.File("flux_profiles_alejandro.h5",  "r") as f:
#         group = f['/']
#         Qes_ions = group['Qes_' + p]
#         xarray = group['xarray_' + p]
#         d = {'Q_es_ions': np.array(Qes_ions), 'x_a': np.array(xarray)}
#         fluxes_alejandro.append(pd.DataFrame(data=d))
# # fluxes_alejandro[2]
# plot_ref = get_plot(fluxes_alejandro[0], label=prefixes_alejandro[0], display_legend=True, width=1500)
# for i in range(len(prefixes_alejandro[1:])):
#     plot_ref = add_to_plot(plot_ref, fluxes_alejandro[i+1], label=prefixes_alejandro[i+1], display_legend=True)
# output_notebook()
# show(plot_ref)


# In[10]:


plot_combi = get_plot(combi_flux[0], label="combi ions" ,display_legend=True, width=1400)
if get_num_species() > 1:
    plot_combi = add_to_plot(plot_combi, combi_flux[1], label="combi electrons", display_legend=True)
plot_combi.output_backend = "svg"
export_svgs([plot_combi], filename=combiSchemeMode[:-4] + "Combi.svg")
output_notebook()
show(plot_combi)


# In[11]:


fluxes_multiplied_with_coefficient = {}
multiplyWithCoefficient = True

for i in range(len(fluxes)):
    probname = combiScheme.iloc[i]['probname']
    fluxes_multiplied_with_coefficient[probname] = {}
    for species in [0,1]:
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
export_svgs([plot_scheme_mult_in_one], filename=combiSchemeMode[:-4] + "AllInOneMultipliedWithCoefficients.svg")
output_notebook()
show(plot_scheme_mult_in_one)


# In[13]:


plot_scheme_in_one = get_plot(fluxes[combiScheme.iloc[0]['probname']][0],label=combiScheme.iloc[0]['probname']+" ions")
if get_num_species() > 1:
    plot_scheme_in_one = add_to_plot(plot_scheme_in_one, fluxes[combiScheme.iloc[0]['probname']][1], label=combiScheme.iloc[0]['probname']+" electrons")
for i in range(len(combiScheme[1:])):
    plot_scheme_in_one = add_to_plot(plot_scheme_in_one, fluxes[combiScheme.iloc[i+1]['probname']][0], label=combiScheme.iloc[i+1]['probname']+" ions")
    if get_num_species() > 1:
        plot_scheme_in_one = add_to_plot(plot_scheme_in_one, fluxes[combiScheme.iloc[i+1]['probname']][1], label=combiScheme.iloc[i+1]['probname']+" electrons")
plot_scheme_in_one.output_backend = "svg"
export_svgs([plot_scheme_in_one], filename=combiSchemeMode[:-4] + "AllInOne.svg")
output_notebook()
show(plot_scheme_in_one)


# In[14]:



plots = [get_plot(fluxes[combiScheme.loc[i]['probname']][s], label=combiScheme.iloc[i]['probname']+str(s), color = 'orange' if combiScheme.iloc[i]['coefficient'] > 0 else 'blue', display_legend=True)for s in [0,1] for i in range(len(fluxes))]
output_notebook()
show(gridplot([plots]))
