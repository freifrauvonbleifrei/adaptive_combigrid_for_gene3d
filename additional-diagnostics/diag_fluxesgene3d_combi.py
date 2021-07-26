# -*- coding: utf-8 -*-

from tkinter import END

import h5py
import matplotlib.pyplot as plt
import numpy as np
import utils.averages as averages
from diagnostics.baseplot import Plotting
from diagnostics.diagnostic import Diagnostic
from pathlib import Path
import utils.unitconversion as unitcon

class DiagFluxesgene3dC(Diagnostic):

    # pylint: disable=invalid-name
    def __init__(self, avail_vars=None, specnames=None):
        super().__init__()
        self.avail_vars = avail_vars
        self.specnames = specnames
        self.opts = {"spec": {'tag': "Spec", 'values': self.list_specnames()}}
        self.set_defaults()
        self.options_gui = Diagnostic.OptionsGUI()

    def set_options(self, run_data, species):

        class FluxStep:
            def __init__(self, specnames, electromagnetic):
                self.__em__ = electromagnetic
                for spec in specnames:
                    setattr(self, spec, self.__SpecSpectra(self.__em__))

            class __SpecSpectra:
                def __init__(self, em):
                    self.Qes = []
                    self.Ges = []

                    if em:
                        self.Qem = []
                        self.Gem = []

        self.geom = run_data.geometry
        self.run_data = run_data
        self.specnames = self.run_data.specnames

        self.flux_step = FluxStep(self.specnames, self.run_data.electromagnetic)

        self.get_needed_vars(['mom'])

        self.get_spec_from_opts()

        self.area = averages.get_area(self.geom)*self.run_data.pars["Lref"]**2  #this normalizes the area
        Qturb = self.run_data.gyrobohm_SI(quantity="Qturb")
        self.norm_q = Qturb

        return self.needed_vars

    def execute(self, data, step, last_step):

        for spec in self.specnames:

            flux_step = getattr(self.flux_step, spec)

            Qes_yz_av = np.average(getattr(getattr(data, 'mom_' + spec), "Q_es")(step.time, step.mom), weights=self.geom.jacobian.T, axis=(1, 2))*self.run_data.pars["temp" + spec]*self.run_data.pars["dens" + spec]
            Ges_yz_av = np.average(getattr(getattr(data, 'mom_' + spec), "Gamma_es")(step.time, step.mom), weights=self.geom.jacobian.T, axis=(1, 2))*self.run_data.pars["dens" + spec]

            flux_step.Qes.append(Qes_yz_av)
            flux_step.Ges.append(Ges_yz_av)

            if self.run_data.electromagnetic:

                Qem_yz_av = np.average(getattr(getattr(data, 'mom_' + spec), "Q_em")(step.time, step.mom), weights=self.geom.jacobian.T, axis=(1, 2))
                Gem_yz_av = np.average(getattr(getattr(data, 'mom_' + spec), "Gamma_em")(step.time, step.mom), weights=self.geom.jacobian.T, axis=(1, 2))

                flux_step.Qem.append(Qem_yz_av)
                flux_step.Gem.append(Gem_yz_av)


    def plot(self, time_requested, output=None, out_folder=None):

        for spec in self.specnames:
            spec_flux = getattr(self.flux_step, spec)
            name_hdf5 = Path('flux_profile_{}'.format(spec) + '.h5')
            file = h5py.File(out_folder / name_hdf5, 'w')
            file["/x_a"] = self.run_data.spatialgrid.x_a
            file["/time"] = time_requested
            file["/SI_conv"] = self.area*self.norm_q*1E-03
            file["/SI_noA_conv"] = self.norm_q*1E+03
            file["/Q_es_" + spec] = np.array(averages.mytrapz(getattr(spec_flux, "Qes"), time_requested))*np.array(self.run_data.pars["dens" + spec])*np.array(self.run_data.pars["temp" + spec])
            print(np.array(self.run_data.pars["dens" + spec]))
            print(np.array(self.run_data.pars["temp" + spec]))
            file["/G_es_" + spec] = averages.mytrapz(getattr(spec_flux, "Ges"), time_requested)
            if self.run_data.electromagnetic:
                file["/Q_em_" + spec] = averages.mytrapz(getattr(spec_flux, "Qem"), time_requested)
                file["/G_em_" + spec] = averages.mytrapz(getattr(spec_flux, "Gem"), time_requested)
