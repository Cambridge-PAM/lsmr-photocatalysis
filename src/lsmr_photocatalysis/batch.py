from .experiment import _Experiment
from .settings import FIGSIZE, DPI, LIGHTBOXBLUE
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm as cm
from matplotlib import colors as colors

class Batch(_Experiment):
    """
    Batch photocatalysis reaction in the lightbox.

    Args:
        data: Data folder name.
        start: Reaction start timestamp in `'hh-mm-ss-msmsms'`.
        c0: Initial concentration/M of substrate of interest.
        details: Details of other reagents and concentrations.
        x: Dictionary key for `SUBSTRATE` properties, see `settings.py`.
        film: Polymer film name. Defaults to no film.
        light: LED colour. Defaults to blue light.
        wavelength: LED wavelength/m. Defaults to 455 nm.
        power: LED power/W. Defaults to `LIGHTBOXBLUE` 105 W m-2.
        correction: Default choice for plots to be baseline corrected or not. Defaults to `True`.
        duration: How long to plot for `concvstime()` and `spectrumvstime()`. Defaults to `1e10`.
    
    Attributes:
        points: Time/s points to extract data, `None` for `Batch`
        x: LaTeX text formatting for substrate concentration.
        x0: LaTeX text formatting for substrate initial concentration.
        peak: Substrate absorption peak wavelength/m.
        e: Substrate extinction coefficient/M m-1.
        l: Path length/m of reactor chamber.
        bregions: Wavelength/nm regions to define as baseline.
        lmin1: Wavelength/nm minimum for `spectrumvstimeplot()`.
        lmax1: Wavelength/nm maximum for `spectrumvstimeplot()`.
        lmin2: Wavelength/nm minimum for `spectrumvsconcplot()`.
        lmax2: Wavelength/nm maximum for `spectrumvsconcplot()`.
    """

    def __init__(self, data: str, start: str, c0: float, details: str, x: str,
                 film: str = 'no', light: str = 'blue', wavelength: str = 455e-9, power: float = LIGHTBOXBLUE,
                 correction: bool = True, duration: float = 1e10):
        super().__init__(data, c0, details, x, film, light, wavelength, power, correction=correction, duration=duration)
        self.start: str = start # override from _Experiment

    def concvstimeplot(self, ln: bool = False, linreg: bool = True, correction: bool | None = None, duration: float | None = None, 
                       exclude: list[int] = [], istart: int = 0, iend: int | None = None, trial: int = 1, 
                       r2threshold: float | None = None, interceptthreshold: float | None = None) -> None:
        """
        Plots concentration/μM vs time/s OR ln(c0-concentration/M) vs time/s for a first order reaction with optional linear regression.
        
        Args:
            ln: Whether to plot ln(c0-concentration/M) vs time/s instead of concentration/μM vs time/s. Defaults to `False`.
            linreg: Whether to perform linear regression at the initial region. Defaults to `True`.
            correction: Overrides instance's default attribute for plots to be baseline corrected or not. Defaults to `None` then points to the instance's attribute.
            duration: Overrides instance's default attribute for how long to plot. Defaults to `None` then points to the instance's attribute.
        """

        t, c = self._concvstime(correction=correction, duration=duration)
        if ln:
            creg = np.log(self.c0-c) # regression with ln concentrations
            cplot = creg # plot in ln concentrations also, no scaling
        else:
            creg = c # regression with concentrations in M
            cplot = c*1e6 # plot concentrations in μM
        n = self._ratephotons() # M s-1 of photons

        fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI, constrained_layout=True)
        ax.plot(t, cplot, color='black', alpha=0.2) # plot data
        if linreg:
            fit, imin, imax = self._linreg(t, creg, exclude=exclude, istart=istart, iend=iend, trial=trial, 
                                           r2threshold=r2threshold, interceptthreshold=interceptthreshold)
            x1 = np.linspace(t[imin], t[imax], 100)
            x2 = np.linspace(t[imax], t[imax]*2, 100) # extrapolate linear regression line
            y1 = (fit.slope*x1 + fit.intercept)
            y2 = (fit.slope*x2 + fit.intercept)
            if ln == False:
                y1 *= 1e6 # plot in μM
                y2 *= 1e6 # plot in μM
            ax.plot(x1, y1, color='black') # plot linear regression line over the initial linear region
            ax.plot(x2, y2, color='black', linestyle='--') # extrapolate linear regression line slightly
            a1, b1, ha1, va1, a2, b2, ha2, va2 = self._textalign(fit.slope)
            if ln: # linear regression results for ln concentrations
                ax.text(a1, b1,
                        rf'$\Phi$ = $\mathrm{{\frac{{rate_{{reaction}}}}{{rate_{{photon}}}}}}$ = {self._numerr(-fit.slope*self.c0/n, fit.stderr*self.c0/n, prefix=1e-2)} %''\n'
                        rf'$k$ = $-$slope = {self._numerr(-fit.slope, fit.stderr, prefix=1e-6)} $\times$ 10$^{{-6}}$ s$^{{-1}}$''\n'
                        f'intercept = {self._numerr(fit.intercept, fit.intercept_stderr, sci=False)}\n'
                        f'$r^2$ = {fit.rvalue**2:.3f}',
                        transform=ax.transAxes, ha=ha1, va=va1)
            else: # linear regression results for concentrations
                ax.text(a1, b1,
                        rf'$\Phi$ = $\mathrm{{\frac{{rate_{{reaction}}}}{{rate_{{photon}}}}}}$ = {self._numerr(fit.slope/n, fit.stderr/n, prefix=1e-2)} %''\n'
                        rf'$k$ = {self._numerr(fit.slope/self.c0, fit.stderr/self.c0, prefix=1e-6)} $\times$ 10$^{{-6}}$ s$^{{-1}}$''\n'
                        f'rate = slope = {self._numerr(fit.slope, fit.stderr, prefix=1e-9)} nM s$^{{-1}}$\n'
                        f'intercept = {self._numerr(fit.intercept, fit.intercept_stderr, prefix=1e-6)} μM\n'
                        f'$r^2$ = {fit.rvalue**2:.3f}',
                        transform=ax.transAxes, ha=ha1, va=va1)
        else:
            a1, b1, ha1, va1, a2, b2, ha2, va2 = self._textalign(-1 if ln else 1)

        ax.text(a2, b2, f'{self.x0} = {self.c0*1e6:.0f} μM\n{self.details}', transform=ax.transAxes, ha=ha2, va=va2) # reaction conditions

        fig.canvas.draw() # calculate ticks preliminarily
        xscale = np.diff(ax.get_xticks())[0]
        yscale = np.diff(ax.get_yticks())[0]
        xmin = np.floor(min(t)/xscale)*xscale
        xmax = np.ceil(max(t)/xscale)*xscale
        if abs(min(cplot)) < yscale: # if the smallest c value is approximately zero
            ymin = 0

        fig.canvas.draw() # calculate ticks preliminarily
        xscale = np.diff(ax.get_xticks())[0]
        yscale = np.diff(ax.get_yticks())[0]
        xmin = np.floor(min(t)/xscale)*xscale
        xmax = np.ceil(max(t)/xscale)*xscale
        ymin = np.floor((min(cplot)-0.3*yscale)/yscale)*yscale
        ymax = np.ceil((max(cplot)+0.3*yscale)/yscale)*yscale
        if abs(min(cplot)) <= 0.3*yscale: # if the smallest c value is approximately zero
            ymin = 0

        ax.set_xlim(xmin, xmax)
        ax.set_xticks(np.arange(xmin, xmax+abs(1e-10*xmax), xscale))
        ax.set_ylim(ymin, ymax)
        ax.set_yticks(np.arange(ymin, ymax+abs(1e-10*ymax), yscale))
        ax.set_xlabel(f'{self.ttext} / s')
        if ln:
            ax.set_ylabel(f'ln({self.x0}$-${self.x} / M)')
            ax.set_title(f'ln({self.x0}$-${self.x} / M) against {self.ttext} / s {self._desc()}')
        else:
            ax.set_ylabel(f'{self.x} / μM')
            ax.set_title(f'{self.x} / μM against {self.ttext} / s {self._desc()}')
        plt.show()

    def spectrumvstimeplot(self, correction: bool | None = None, duration: float | None = None,
                           lmin: float | None = None, lmax: float | None = None, amin: float | None = None, amax: float | None = None) -> None:
        """
        Plots the UV-Vis absorption spectrum as a heat map over time/s

        Args:
            correction: Overrides instance's default attribute for plots to be baseline corrected or not. Defaults to `None` then points to the instance's attribute.
            duration: Overrides instance's default attribute for how long to plot. Defaults to `None` then points to the instance's attribute.
            lmin: Wavelength/nm minimum. Defaults to `None` then points to the instance's attribute `lmin2`.
            lmax: Wavelength/nm maximum. Defaults to `None` then points to the instance's attribute `lmin2`.
            amin: Absorbance minimum. Defaults to `None` then points to the minimum value in the data.
            amax: Absorbance maximum. Defaults to `None` then points to the maximum value in the data.
        """

        if lmin is None:
            lmin = self.lmin2
        if lmax is None:
            lmax = self.lmax2

        t, l, a = self._spectrumvstime(correction=correction, duration=duration)
        mask = (lmin <= l) & (l <= lmax) 
        a = a[mask, :]
        l = l[mask]

        if amin == None:
            amin = np.floor(a.min()/0.1)*0.1 # round down to nearest 0.1
        if amax == None:
            amax = np.ceil(a.max()/0.1)*0.1 # round up to nearest 0.1

        fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI, constrained_layout=True)
        cmap = plt.get_cmap('viridis_r')
        ax.contourf(t, l, a, levels=100, vmin=amin, vmax=amax, cmap=cmap)
        ax.axhline(y=self.peak, color='black', linestyle='--') # plot vertical line at the peak
        ax.text(0.98, 0.98, f'{self.x0} = {self.c0*1e6:.0f} μM\n{self.details}', transform=ax.transAxes, ha='right', va='top') # reaction conditions

        norm = colors.Normalize(vmin=amin, vmax=amax)
        sm = cm.ScalarMappable(norm=norm, cmap=cmap)
        cbar = plt.colorbar(sm, ax=ax)
        cbar.set_ticks(np.arange(amin, amax+abs(1e-10*amax), 0.1))
        cbar.set_label("absorbance")

        fig.canvas.draw() # calculate ticks preliminarily
        xscale = np.diff(ax.get_xticks())[0]
        yscale = np.diff(ax.get_yticks())[0]
        xmin = np.floor(min(t)/xscale)*xscale

        ax.set_xticks(np.arange(xmin, np.ceil(max(t))+abs(1e-10*max(t)), xscale))
        ax.set_ylim(lmin, lmax)
        ax.set_yticks(np.arange(lmin, lmax+abs(1e-10*lmax), yscale))
        ax.set_xlabel(f'{self.ttext} / s')
        ax.set_ylabel(r'$\lambda$ / nm')
        ax.set_title(rf'UV-Vis spectra for {self.x[1:-1]}{self._textcorr(correction)[0]} against {self.ttext} / s {self._desc()}')
        plt.show()

    # uv-vis absorption spectrum at time t for each run
    @staticmethod
    def spectrumvsconcplot(expts: list[Batch], t: float, correction: bool | None = None, inspect: int = 0,
                           details: str | None = None, cbarlabel: str | None = None, cbarticks: list[float] | None = None, 
                           lmin: float | None = None, lmax: float | None = None) -> None:
        """
        Plots the UV-Vis absorption spectrum at a common time/s point for a series of experiments on the same axis.

        Args:
            expts: List of `Batch` experiment objects to plot.
            t: Common time/s point where the spectra are taken.
            correction: Overrides instance's default attribute for plots to be baseline corrected or not. Defaults to `None` then points to the instance's attribute.
            inspect: Index of the time/s point for which to show baseline correction. Defaults to `0`, i.e. the first time/s point.
            details: Overrides instance's default attribute for the text description of the experimental details. Defaults to `None` then points to the instance's attribute.
            cbarlabel: Label for the independent variable that differs between the `Batch` experiments. Defaults to `None` then points to inital concentration `c0`.
            cbarticks: Tick lables for the independent variable that differs between the `Batch` experiments. Defaults to `None` then points to inital concentration `c0`.
            lmin: Wavelength/nm minimum. Defaults to `None` then points to the instance's attribute `lmin1`.
            lmax: Wavelength/nm maximum. Defaults to `None` then points to the instance's attribute `lmin1`.
        """

        if lmin is None:
            lmin = expts[0].lmin1
        if lmax is None:
            lmax = expts[0].lmax1

        if np.isclose(t % 15, 0): # if values of t are divisible by 15 s, quote t in min, otherwise in s
            tunit = 'min'
            tfactor = 60
        else:
            tunit = 's'
            tfactor = 1

        if correction == None:
            correction = expts[0].correction # use the default correction choice if not specified
        if correction:
            plots = [True, False] # plot both corrected and non-corrected for comparison if correction is specified
        else:
            plots = [False] # plot only non-corrected for comparison if no correction is specified

        for corr in plots:
            fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI, constrained_layout=True)
            cmap = plt.get_cmap('viridis_r')
            col = cmap(np.linspace(0, 1, len(expts)))

            amax = 0
            amin = 10
            for i, expt in enumerate(expts):
                l, a, b = expt._spectrumvstimepoints(points=[t], correction=corr)
                mask = (lmin <= l) & (l <= lmax)
                l = l[mask]
                a = a[:, mask]
                b = b[:, mask]
                ax.plot(l, a[0], color=col[i])
                    
                if np.max(a) > amax:
                    amax = np.max(a)
                if np.min(a) < amin: # if found a new minimum
                    amin = np.min(a)
                    if l[a[0] == amin][0] >= (lmin+lmax)/2: # if the lowest point in the plot is on the right
                        pos, ha = 0.02, 'left' # align baseline label to the left
                        start = 0
                        end = int(len(l)/2)
                    else:
                        pos, ha = 0.98, 'right' # align baseline label to the right
                        start = int(len(l)/2)
                        end = len(l)
                if i == inspect:
                    b0 = b[0] # save the baseline for plotting later

            ax.axvline(x=expts[0].peak, color='black', linestyle='--', alpha=0.2) # plot vertical line at the peak

            label, params = expts[0]._textcorr(corr)
            if details is None:
                details = expts[0].details
            ax.text(0.98, 0.98, f'{expts[0].ttext} = {t/tfactor} {tunit}\n\n{details}\n\n{params}', 
                    transform=ax.transAxes, ha='right', va='top') # reaction conditions and baseline correction parameters
            
            norm = colors.Normalize(vmin=0, vmax=len(expts)-1)
            sm = cm.ScalarMappable(norm=norm, cmap=cmap)
            cbar = fig.colorbar(sm, ax=ax)
            if cbarlabel is None:
                cbarticks = [f'{expt.c0*1e6:.0f}' for expt in expts]
                cbarlabel = f'{expts[0].x0} / μM'
            cbar.set_ticks(range(len(expts)), labels=cbarticks)
            cbar.set_label(cbarlabel)

            fig.canvas.draw() # calculate ticks preliminarily
            yscale = np.diff(ax.get_yticks())[0]
            ymax = np.ceil((amax+0.3*yscale)/yscale)*yscale # automatic y axis maximum
            if corr: 
                ymin = 0 # set minimum to 0 for baseline corrected spectra
            else:
                ymin = np.floor((amin-0.3*yscale)/yscale)*yscale # automatic y axis minimum
                if correction:
                    ax.plot(l, b0, color='black', linestyle = '--') # plot the first baseline for non-corrected spectra
                    ax.text(pos, min(b0[start:end])-0.1*yscale, 
                            f'interpolated baseline for\n{cbarlabel.split('/')[0].strip()} = {cbarticks[inspect]} {cbarlabel.split('/')[-1].strip()}', 
                            transform=ax.get_yaxis_transform(), ha=ha, va='top')

            ax.set_xlim(lmin, lmax)
            ax.set_ylim(ymin, ymax)
            ax.set_yticks(np.arange(ymin, ymax+abs(1e-10*ymax), yscale))
            ax.set_xlabel(r'$\lambda$ / nm')
            ax.set_ylabel('absorbance')
            ax.set_title(f'UV-Vis spectra for {expts[0].x[1:-1]}{label} against {cbarlabel} {expts[0]._desc()}')
            plt.show()