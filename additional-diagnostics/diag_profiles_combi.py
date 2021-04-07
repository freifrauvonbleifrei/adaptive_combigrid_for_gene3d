#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import h5py
import matplotlib.pyplot as plt
import numpy as np

import utils.averages as averages
from diagnostics.baseplot import Plotting
from diagnostics.diagnostic import Diagnostic

from pathlib import Path

class DiagProfilesC(Diagnostic):
    def __init__(self, avail_vars=None, specnames=None):
        super().__init__()

        self.avail_vars = avail_vars
        self.specnames = specnames
        self.opts = {"spec": {'tag': "Spec", 'values': self.list_specnames()}}
        self.set_defaults()
        self.options_gui = Diagnostic.OptionsGUI()

    def set_options(self, run_data, species):

        class ProfilesStep:
            def __init__(self, specnames):
                for spec in specnames:
                    setattr(self, spec, self.__Specprofiles())

            class __Specprofiles:
                def __init__(self):
                    self.T = []
                    self.n = []
                    self.u = []

        self.run_data = run_data

        self.geom = run_data.geometry

        self.rhostarref = self.run_data.pars["rhostar"] * self.run_data.pars["minor_r"]

        self.profiles_step = ProfilesStep(self.specnames)

        self.get_needed_vars(['mom'])

        self.get_spec_from_opts()

        return self.needed_vars

    def execute(self, data, step, last_step):
        for spec in self.specnames:

            profiles_step = getattr(self.profiles_step, spec)

            # 3D
            dens = getattr(getattr(data, 'mom_' + spec), "n")(step.time, step.mom)
            T_par = getattr(getattr(data, 'mom_' + spec), "T_par")(step.time, step.mom)
            T_perp = getattr(getattr(data, 'mom_' + spec), "T_per")(step.time, step.mom)
            u_par = getattr(getattr(data, 'mom_' + spec), "u_par")(step.time, step.mom)

            profiles_step.T.append(np.average(np.squeeze(1/3 * T_par + 2/3 * T_perp),
                                                  weights=self.geom.jacobian.T, axis=(1, 2)))

            profiles_step.n.append(np.average(np.squeeze(dens),
                                                  weights=self.geom.jacobian.T, axis=(1, 2)))

            profiles_step.u.append(np.average(np.squeeze(u_par),
                                                  weights=self.geom.jacobian.T, axis=(1, 2)))

    def plot(self, time_requested, output=None, out_folder=None):

        if self.run_data.x_global or self.run_data.is3d:
            y_lbl = r'$x/a$'
            y_ax = self.run_data.spatialgrid.x_a
        else:
            y_lbl = r'$x/\rho_{ref}$'
            y_ax = self.run_data.spatialgrid.x

        for i_s, spec in enumerate(self.specnames):

            omn_b = self.run_data.profilesdata.omn0s[:, i_s]
            omt_b = self.run_data.profilesdata.omt0s[:, i_s]

            temp = np.array(getattr(getattr(self.profiles_step, spec), 'T'))*self.rhostarref*self.run_data.pars["temp" + spec]*self.run_data.pars["Tref"]
            n = np.array(getattr(getattr(self.profiles_step, spec), 'n'))*self.rhostarref*self.run_data.pars["dens" + spec]*self.run_data.pars["nref"]
            u = np.array(getattr(getattr(self.profiles_step, spec), 'u'))*self.rhostarref

            if (self.run_data.pars["radial_dependence"]): 

                T_b = self.run_data.profilesdata.T0s[:, i_s]
                n_b = self.run_data.profilesdata.n0s[:, i_s]

            else:

                T_b = np.exp(-omt_b*(self.run_data.spatialgrid.x_a - self.run_data.pars["x0"]))*self.run_data.pars["temp" + spec]*self.run_data.pars["Tref"]
                n_b = np.exp(-omn_b*(self.run_data.spatialgrid.x_a - self.run_data.pars["x0"]))*self.run_data.pars["dens" + spec]*self.run_data.pars["nref"]

            omt = -np.gradient(np.log(T_b + temp), self.run_data.spatialgrid.x_a, axis=1) / \
                self.run_data.pars["minor_r"]
            omn = -np.gradient(np.log(n_b + n), self.run_data.spatialgrid.x_a, axis=1) / \
                self.run_data.pars["minor_r"]

            name_hdf5 = Path('profile_{}'.format(spec) + '.h5')
            file = h5py.File(out_folder / name_hdf5, 'w')
            file["/x_a"] = y_ax
            file["/n_b_" + spec] = n_b
            file["/T_b_" + spec] = T_b
            file["/omn_b_" + spec] = omn_b
            file["/omt_b_" + spec] = omt_b
            file["/n_" + spec] = n_b + averages.mytrapz(n, time_requested)
            file["/T_" + spec] = T_b + averages.mytrapz(temp, time_requested)
            file["/omn_" + spec] = averages.mytrapz(omn, time_requested)
            file["/omT_" + spec] = averages.mytrapz(omt, time_requested)
            file.close()
