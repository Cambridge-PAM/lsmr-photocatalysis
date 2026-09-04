from .experiment import _Experiment
from .settings import FIGSIZE, DPI
import numpy as np
from numpy.typing import NDArray
import matplotlib.pyplot as plt
from matplotlib import cm as cm
from matplotlib import colors as colors

class StoppedFlow(_Experiment):
    """
    Stopped-flow with reaction/LED/film in 50 ul reactor chamber.
    Irradiation for each reaction time, then reactor chamber immediately moved to UV-Vis.
    Flushed with fresh starting reaction mixture and repeated for a new reaction time.

    Args:
        data: Data folder name.
        t: Reaction time/s points.
        c0: Initial concentration/M of substrate of interest.
        details: Details of other reagents and concentrations.
        x: Dictionary key for `SUBSTRATE` properties, see `settings.py`.
        film: Polymer film name. Defaults to no film.
        light: LED colour. Defaults to blue light.
        wavelength: LED wavelength/m.
        power: LED power/W.
        correction: Default choice for plots to be baseline corrected or not. Defaults to `True`.
        reverse: Denotes that the reaction time/s points were recorded from largest to smallest.
    
    Attributes:
        points: Time/s points to extract data, `None` for `StoppedFlow`.
        start: Meaningless `''` for `StoppedFlow`.
        x: LaTeX text formatting for substrate concentration.
        x0: LaTeX text formatting for substrate initial concentration.
        peak: Substrate absorption peak wavelength/m.
        e: Substrate extinction coefficient/M m-1.
        l: Path length/m of reactor chamber.
        bregions: Wavelength/nm regions to define as baseline.
        lmin1: Wavelength/nm minimum for `spectrumvstimeplot()`.
        lmax1: Wavelength/nm maximum for `spectrumvstimeplot()`.
        lmin2: Meaningless for `StoppedFlow`.
        lmax2: Meaningless for `StoppedFlow`.
        duration: Meaningless for `StoppedFlow`. Defaults to `1e10`.
    """

    def __init__(self, data: str, t: NDArray[np.float64], c0: float, details: str, x: str, film: str, light: str, wavelength: float, power: float,
                 correction: bool = True, reverse: bool = False):
        super().__init__(data, c0, details, x, film, light, wavelength, power, correction=correction)
        if reverse:
            self.t: NDArray[np.float64] = t[::-1]
        else:
            self.t: NDArray[np.float64] = t # reaction times in s, new attribute not inherited from _Experiment
        self.reverse: bool = reverse # whether the reaction times were recorded in descending order, new attribute not inherited from _Experiment

    def concvstimeplot(self, ln: bool = False, correction: bool | None = None, 
                       exclude: list[int] = [], istart: int = 0, iend: int | None = None, trial: int = 1, 
                       r2threshold: float | None = None, interceptthreshold: float | None = None) -> None:
        """
        Plots concentration/mM vs time/s OR ln(c0-concentration/M) vs time/s for a first order reaction with optional linear regression.
        
        Args:
            ln: Whether to plot ln(c0-concentration/M) vs time/s instead of concentration/mM vs time/s. Defaults to `False`.
            linreg: Whether to perform linear regression at the initial region. Defaults to `True`.
            correction: Overrides instance's default attribute for plots to be baseline corrected or not. Defaults to `None` then points to the instance's attribute.
            duration: Overrides instance's default attribute for how long to plot. Defaults to `None` then points to the instance's attribute.
        """
        
        c = self._concvstimepoints(correction=correction)
        if self.reverse:
            c = c[::-1]
        if ln:
            creg = np.log(self.c0-c) # regression with ln concentrations
            cplot = creg # plot in ln concentrations also, no scaling
        else:
            creg = c # regression with concentrations in M
            cplot = c*1e3 # plot concentrations in mM
        n = self._ratephotons() # M s-1 of photons
            
        fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI, constrained_layout=True)
        ax.plot(self.t, cplot, 'x', color='black') # plot data points
        fit, imin, imax = self._linreg(self.t, creg, exclude=exclude, istart=istart, iend=iend, trial=trial, 
                                       r2threshold=r2threshold, interceptthreshold=interceptthreshold)
        x = np.linspace(self.t[imin], self.t[imax], 100)
        y = (fit.slope*x + fit.intercept)
        if ln == False:
            y = y*1e3 # plot in mM
        ax.plot(x, y, color='black') # plot linear regression line over the initial linear region
        ax.plot(self.t[exclude], cplot[exclude], 'o', markersize=8, markerfacecolor='none', markeredgecolor='red') # highlight excluded points

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
                    rf'rate = slope = {self._numerr(fit.slope, fit.stderr, prefix=1e-6)} $\mathrm{{\mu}}$M s$^{{-1}}$''\n'
                    f'intercept = {self._numerr(fit.intercept, fit.intercept_stderr, prefix=1e-3)} mM\n'
                    f'$r^2$ = {fit.rvalue**2:.3f}',
                    transform=ax.transAxes, ha=ha1, va=va1)
        ax.text(a2, b2, f'{self.x0} = {self.c0*1e3:.0f} mM\n{self.details}', transform=ax.transAxes, ha=ha2, va=va2) # reaction conditions

        fig.canvas.draw() # calculate ticks preliminarily
        xscale = np.diff(ax.get_xticks())[0]
        yscale = np.diff(ax.get_yticks())[0]
        xmin = max(np.floor((min(self.t)-0.3*xscale)/xscale)*xscale, 0) # dont let x axis t < 0
        xmax = np.ceil((max(self.t)+0.3*xscale)/xscale)*xscale
        ymin = np.floor((min(cplot)-0.3*yscale)/yscale)*yscale
        ymax = np.ceil((max(cplot)+0.3*yscale)/yscale)*yscale
        if fit.slope == 0: # give more y axis ticks
            ymin -= yscale
            ymax += yscale

        ax.set_xlim(xmin, xmax)
        ax.set_xticks(np.arange(xmin, xmax+abs(1e-10*xmax), xscale))
        ax.set_ylim(ymin, ymax)
        ax.set_yticks(np.arange(ymin, ymax+abs(1e-10*ymax), yscale))
        ax.set_xlabel(f'{self.ttext} / s')
        if ln:
            ax.set_ylabel(rf'ln({self.x0}$-${self.x} / M)')
            ax.set_title(rf'ln({self.x0}$-${self.x} / M)'f' against {self.ttext} / s {self._desc()}')
        else:
            ax.set_ylabel(f'{self.x} / mM')
            ax.set_title(f'{self.x} / mM against {self.ttext} / s {self._desc()}')
        plt.show()

    def spectrumvstimeplot(self, correction: bool | None = None, inspect: int = 0, lmin: float | None = None, lmax: float | None = None) -> None:
        """
        Plots the UV-Vis absorption spectrum at each time/s point on the same axis.

        Args:
            correction: Overrides instance's default attribute for plots to be baseline corrected or not. Defaults to `None` then points to the instance's attribute.
            inspect: Index of the time/s point for which to show baseline correction. Defaults to `0`, i.e. the first time/s point.
            lmin: Wavelength/nm minimum. Defaults to `None` then points to the instance's attribute `lmin1`.
            lmax: Wavelength/nm maximum. Defaults to `None` then points to the instance's attribute `lmin1`.
        """

        if lmin is None:
            lmin = self.lmin1
        if lmax is None:
            lmax = self.lmax1

        if np.isclose(self.t[0] % 15, 0): # if values of t are divisible by 15 s, quote t in min, otherwise in s
            tunit = 'min'
            tfactor = 60
        else:
            tunit = 's'
            tfactor = 1

        if correction == None:
            correction = self.correction # use the default correction choice if not specified
        if correction:
            plots = [True, False] # plot both corrected and non-corrected for comparison if correction is specified
        else:
            plots = [False] # plot only non-corrected for comparison if no correction is specified
        
        for corr in plots:
            l, a, b = self._spectrumvstimepoints(correction=corr)
            mask = (lmin <= l) & (l <= lmax)
            l = l[mask]
            a = a[:, mask]
            b = b[:, mask]

            if self.reverse:
                a = a[::-1]
                b = b[::-1]

            fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI, constrained_layout=True)
            cmap = plt.get_cmap('viridis_r')
            col = cmap(np.linspace(0, 1, len(self.t)))

            for i in range(len(self.t)):
                ax.plot(l, a[i], color=col[i])
            ax.axvline(x=self.peak, color='black', linestyle='--', alpha=0.2) # plot vertical line at the peak

            label, params = self._textcorr(corr)
            ax.text(0.98, 0.98, f'{self.x0} = {self.c0*1e3:.0f} mM\n{self.details}\n\n{params}', 
                    transform=ax.transAxes, ha='right', va='top') # reaction conditions and baseline correction parameters

            norm = colors.Normalize(vmin=0, vmax=len(self.t)-1)
            sm = cm.ScalarMappable(norm=norm, cmap=cmap)
            cbar = fig.colorbar(sm, ax=ax)
            cbar.set_ticks(range(len(self.t)), labels=self.t/tfactor)
            cbar.set_label(f'{self.ttext} / {tunit}')

            fig.canvas.draw() # calculate ticks preliminarily
            yscale = np.diff(ax.get_yticks())[0]
            ymax = np.ceil((np.max(a)+0.3*yscale)/yscale)*yscale # automatic y axis maximum
            if corr: 
                ymin = 0 # set minimum to 0 for baseline corrected spectra
            else:
                ymin = np.floor((np.min(a)-0.3*yscale)/yscale)*yscale # automatic y axis minimum
                if correction:
                    if l[np.unravel_index(np.argmin(a), a.shape)[1]] >= (lmin+lmax)/2: # if the lowest point in the plot is on the right
                        pos, ha = 0.02, 'left' # align baseline label to the left
                        start = 0
                        end = int(len(l)/2)
                    else:
                        pos, ha = 0.98, 'right' # align baseline label to the right
                        start = int(len(l)/2)
                        end = len(l)
                    ax.plot(l, b[inspect], color='black', linestyle = '--') # plot the first baseline for non-corrected spectra
                    ax.text(pos, np.min(b[:, start:end])-0.1*yscale,
                            f'interpolated baseline for\n{self.ttext} = {self.t[inspect]/tfactor} {tunit}', 
                            transform=ax.get_yaxis_transform(), ha=ha, va='top')

            ax.set_xlim(lmin, lmax)
            ax.set_ylim(ymin, ymax)
            ax.set_yticks(np.arange(ymin, ymax+abs(1e-10*ymax), yscale))
            ax.set_xlabel(r'$\lambda$ / nm')
            ax.set_ylabel('absorbance')
            ax.set_title(f'UV-Vis spectra for {self.x[1:-1]}{label} against {self.ttext} / {tunit} {self._desc()}')
            plt.show()