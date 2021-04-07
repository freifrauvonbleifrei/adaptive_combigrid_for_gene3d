from tkinter import END
import matplotlib.pyplot as plt
import numpy as np
import h5py
import utils.aux_func as aux_func
import utils.averages as averages
from diagnostics.diagnostic import Diagnostic
from diagnostics.baseplot import Plotting
from pathlib import Path


class DiagFluxSpectragene3dC(Diagnostic):
    # pylint: disable=invalid-name
    def __init__(self, avail_vars=None, specnames=None):
        super().__init__()
     
        self.avail_vars = avail_vars
        self.specnames = specnames
   
        self.opts = {"spec": {'tag': "Spec", 'values': self.list_specnames()}}

        self.set_defaults()
        self.options_gui = Diagnostic.OptionsGUI()

    def set_options(self, run_data, specnames=None):

        class FluxStep:
            def __init__(self, specnames, electromagnetic):
                self.__em__ = electromagnetic
                for spec in specnames:
                    setattr(self, spec, self.__SpecSpectra(self.__em__))

            class __SpecSpectra:
                def __init__(self, em):
                    self.Qes = self.__FluxSpectra()
                    self.Ges = self.__FluxSpectra()

                    if em:
                        self.Qem = self.__FluxSpectra()
                        self.Gem = self.__FluxSpectra()

                class __FluxSpectra:
                    def __init__(self):
                        self.ky = []
                        self.t = []

        self.geom = run_data.geometry
        self.run_data = run_data
        self.specnames = self.run_data.specnames

        self.flux_step = FluxStep(self.specnames, self.run_data.electromagnetic)

        self.get_needed_vars(['field', 'mom'])

        # This is special since we dont save stuff by default, so just do what we need
        self.get_spec_from_opts()

        return self.needed_vars

    def execute(self, data, step, last_step):

        def flux_spectra_xz_global_av(var, min_x, max_x):
            if self.run_data.x_global:
                jaco3d = np.broadcast_to(self.geom.jacobian[:, np.newaxis, :], (
                    self.run_data.pnt.nz0, self.run_data.pnt.nky0, self.run_data.pnt.nx0))
                return 2*np.average(var[: , :, :],
                                    weights=jaco3d.T[: , :, :], axis=(0, 2))
            elif self.run_data.is3d:
                jac_xz = np.mean(self.geom.jacobian.T, axis=1)
                jaco3d_xz = np.broadcast_to(jac_xz[:, np.newaxis, :], (
                    self.run_data.pnt.nx0, self.run_data.pnt.ny0, self.run_data.pnt.nz0))
                return np.average(var[: , :, :],
                                  weights=jaco3d_xz[: , :, :], axis=(0, 2))

        def compute_ky(flux_spec, var, geometry):
            temp = flux_spectra_xz_global_av(var, 0, self.run_data.pnt.nx0-1)
            flux_spec.ky.append(temp)
            flux_spec.t.append(np.sum(temp))
            return flux_spec

        phi = getattr(getattr(data, 'field'), "phi")(step.time, step.field)

        if self.run_data.x_global:
            vE_x = -1j*self.run_data.spatialgrid.ky[np.newaxis, :, np.newaxis]*phi/ \
                   self.run_data.pars["Bref"]/self.geom.Cxy[:, np.newaxis, np.newaxis]

        elif self.run_data.is3d:
            phi_c = np.fft.fft(phi, axis=1)/self.run_data.pnt.ny0
            vE_x = -1j*self.run_data.spatialgrid.ky[np.newaxis, :, np.newaxis]*phi_c/ \
                   self.run_data.pars["Bref"]/self.geom.Cxy[:, np.newaxis, np.newaxis]

        if self.run_data.electromagnetic:
            A_par = getattr(getattr(data, 'field'), "A_par")(step.time, step.field)

            if self.run_data.x_global:
                B_x = 1j*self.run_data.spatialgrid.ky[np.newaxis, :,
                         np.newaxis]*A_par/self.geom.Cxy[:, np.newaxis, np.newaxis]
            elif self.run_data.is3d:
                A_par_c = np.fft.fft(A_par, axis=1)/self.run_data.pnt.ny0
                B_x = -1j*self.run_data.spatialgrid.ky[np.newaxis, :,
                          np.newaxis]*A_par_c/self.geom.Cxy[:, np.newaxis, np.newaxis]

        for i_s, spec in enumerate(self.specnames):

            n0 = self.run_data.profilesdata.n0s[:, i_s]
            T0 = self.run_data.profilesdata.T0s[:, i_s]
            T_b = T0/self.run_data.pars["temp" + spec]/self.run_data.pars["Tref"]
            n_b = n0/self.run_data.pars["dens" + spec]/self.run_data.pars["nref"]

            if self.run_data.x_global:
                dens = getattr(getattr(data, 'mom_' + spec), "dens")(step.time, step.mom)
                T_par = getattr(getattr(data, 'mom_' + spec), "T_par")(step.time, step.mom)
                T_perp = getattr(getattr(data, 'mom_' + spec), "T_perp")(step.time, step.mom)
                u_par = getattr(getattr(data, 'mom_' + spec), "u_par")(step.time, step.mom)

            elif self.run_data.is3d:
                dens_real = getattr(getattr(data, 'mom_' + spec), "n")(step.time, step.mom)
                T_par_real = getattr(getattr(data, 'mom_' + spec), "T_par")(step.time, step.mom)
                T_perp_real = getattr(getattr(data, 'mom_' + spec), "T_per")(step.time, step.mom)
                u_par_real = getattr(getattr(data, 'mom_' + spec), "u_par")(step.time, step.mom)

                dens = np.fft.fft(dens_real, axis=1)/self.run_data.pnt.ny0
                T_par = np.fft.fft(T_par_real, axis=1)/self.run_data.pnt.ny0
                T_perp = np.fft.fft(T_perp_real, axis=1)/self.run_data.pnt.ny0
                u_par = np.fft.fft(u_par_real, axis=1)/self.run_data.pnt.ny0

            flux_step = getattr(self.flux_step, spec)

            G_es = np.conj(vE_x)*dens
            flux_step.Ges = compute_ky(flux_step.Ges, G_es, self.geom)

            var1 = (0.5*T_par + T_perp)*n_b[:, np.newaxis, np.newaxis] + 1.5*dens*T_b[:, np.newaxis,
                                                                                  np.newaxis]
            Q_es = np.conj(vE_x)*var1
            flux_step.Qes = compute_ky(flux_step.Qes, Q_es, self.geom)

            if self.run_data.is3d:

                jac_xz = np.mean(self.geom.jacobian.T, axis=1)
                jaco3d_xz = np.broadcast_to(jac_xz[:, np.newaxis, :], (
                    self.run_data.pnt.nx0, self.run_data.pnt.ny0, self.run_data.pnt.nz0))
                G_es_complex_total = np.average(G_es, weights=jaco3d_xz, axis=(0, 2))
                Q_es_complex_total = np.average(Q_es, weights=jaco3d_xz, axis=(0, 2))

            if self.run_data.electromagnetic:
                if self.run_data.x_global:

                    Q_par = getattr(getattr(data, 'mom_' + spec), "q_par")(step.time, step.mom)
                    Q_perp = getattr(getattr(data, 'mom_' + spec), "q_perp")(step.time, step.mom)

                    G_em = B_x*np.conj(u_par)
                    flux_step.Gem = compute_ky(flux_step.Gem, G_em, self.geom)
                    Q_em = B_x*np.conj(Q_par + Q_perp)
                    flux_step.Qem = compute_ky(flux_step.Qem, Q_em, self.geom)

                elif self.run_data.is3d:
                    G_em = B_x*np.conj(u_par)
                    flux_step.Gem = compute_ky(flux_step.Gem, G_em, self.geom)
                    Q_em=G_em*0
                    flux_step.Qem = compute_ky(flux_step.Qem, Q_em, self.geom)


    def plot(self, time_requested, output=None, out_folder=None):
        """ For each selected species we have one figure with six subplots.
            Left is vs. kx, right vs. ky; columnwise we plot, log-log, log-lin, lin-in
            Dashed lines are negative values in log-log plot
            Dashed lines are k multiplied values in log-lin plot"""

        if output:
            output.info_txt.insert(END, "Flux specta:\n")

        ky = self.run_data.spatialgrid.ky
        if self.run_data.x_global:
            kymax = self.run_data.pnt.nky0
        elif self.run_data.is3d:
            kymax = int(self.run_data.pnt.ny0/2)

        for spec in self.specnames:

            spec_flux = getattr(self.flux_step, spec)

            for flux in vars(spec_flux).keys():

                flux_ky = np.real(averages.mytrapz(getattr(getattr(spec_flux, flux), 'ky'), time_requested))

                name_hdf5 = Path('flux_spectra_{}_{}'.format(flux, spec) +  '.h5')
                file = h5py.File(out_folder / name_hdf5, 'w')
                file["/ky"] = ky[0:kymax]
                file["/" + flux + '_ky'] = flux_ky[0:kymax]
                file.close()
