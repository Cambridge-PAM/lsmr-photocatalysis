from .stoppedflow import StoppedFlow
from .settings import FIGSIZE, DPI
import numpy as np
from numpy.typing import NDArray
import matplotlib.pyplot as plt
from matplotlib import cm as cm
from matplotlib import colors as colors

class ContinuousFlow(StoppedFlow):
    """
    Continuous flow with reaction/LED/film and UV-Vis at opposite ends of a single 250 ul reactor chamber.
    Flow rate changes according to flow scheme at different time points to change residence time tR.

    Args:
        data: Data folder name.
        tRs: Residence time/s points.
        points: Time/s points to extract data.
        start: Reaction start timestamp in `'hh-mm-ss-msmsms'`.
        c0: Initial concentration/M of substrate of interest.
        details: Details of other reagents and concentrations.
        x: Dictionary key for `SUBSTRATE` properties, see `settings.py`.
        film: Polymer film name. Defaults to no film.
        light: LED colour. Defaults to blue light.
        wavelength: LED wavelength/m.
        power: LED power/W.
        correction: Default choice for plots to be baseline corrected or not. Defaults to `True`.
        reverse: Denotes that the reaction time/s points were recorded from largest to smallest.
        duration: How long to plot for `detectorconcplot()` and `detectorspectrumplot()`. Defaults to `1e10`.
    
    Attributes:
        x: LaTeX text formatting for substrate concentration.
        x0: LaTeX text formatting for substrate initial concentration.
        peak: Substrate absorption peak wavelength/m.
        e: Substrate extinction coefficient/M m-1.
        l: Path length/m of reactor chamber.
        bregions: Wavelength/nm regions to define as baseline.
        lmin1: Wavelength/nm minimum for `spectrumvstimeplot()`.
        lmax1: Wavelength/nm maximum for `spectrumvstimeplot()`.
        lmin2: Wavelength/nm minimum for `detectorspectrumplot()`
        lmax2: Wavelength/nm maximum for `detectorspectrumplot()`
    """

    ttext = r'$t_\mathrm{R}$' # residence time

    def __init__(self, data: str, tRs: NDArray[np.float64], points: NDArray[np.float64], start: str, 
                 c0: float, details: str, x: str, film: str, light: str, wavelength: float, power: float,
                 correction: bool = True, reverse: bool = False, duration: float = 1e10):
        super().__init__(data, tRs, c0, details, x, film, light, wavelength, power, correction=correction, reverse=reverse)
        self.points: NDArray[np.float64] = points # override from _Experiment
        self.start: str = start
        self.duration: float = duration

    def detectorconcplot(self, correction: bool | None = None, duration: float | None = None) -> None:
        """
        Plots the concentration/M of the substrate passing through the UV-Vis detector over course of the entire experiment.

        Args:
            correction: Overrides instance's default attribute for plots to be baseline corrected or not. Defaults to `None` then points to the instance's attribute.
            duration: Overrides instance's default attribute for how long to plot. Defaults to `None` then points to the instance's attribute.
        """
        
        t, c = self._concvstime(correction=correction, duration=duration)
        crec = self._concvstimepoints(correction=correction)
        cplot = c*1e3 # plot concentrations in mM
        crecplot = crec*1e3

        fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI, constrained_layout=True)
        ax.plot(t, cplot, color='black', alpha=0.2) # plot the overall detector reading
        ax.plot(self.points, crecplot, 'x', color='black') # plot each recorded point

        fig.canvas.draw() # calculate ticks preliminarily
        xscale = np.diff(ax.get_xticks())[0]
        yscale = np.diff(ax.get_yticks())[0]
        xmin = np.floor(min(t)/xscale)*xscale
        xmax = np.ceil(max(t)/xscale)*xscale
        ymin = np.floor((min(cplot)-0.3*yscale)/yscale)*yscale
        ymax = np.ceil((max(cplot)+0.3*yscale)/yscale)*yscale
        if abs(min(cplot)) <= 0.3*yscale: # if the smallest c value is approximately zero
            ymin = 0

        mask = ((xmax - 2*xscale) <= t) & (t <= xmax) # mask near the end of the plot
        if np.average(cplot[mask], weights=np.linspace(0, 1, len(cplot[mask]))) <= (ymin+ymax)/2: # average of c values weighted towards indices on the right 
            b, ha, va = 0.97, 'right', 'top' # put the text on the top right if weighted average is on the bottom
        else:
            b, ha, va = 0.03, 'right', 'bottom' # put the text on the bottom right if weighted average is on the bottom
        ax.text(0.98, b, f'{self.x0} = {self.c0*1e3:.0f} mM\n{self.details}', transform=ax.transAxes, ha=ha, va=va) # reaction conditions

        ax.set_xlim(xmin, xmax)
        ax.set_xticks(np.arange(xmin, xmax+abs(1e-10*xmax), xscale))
        ax.set_ylim(ymin, ymax)
        ax.set_yticks(np.arange(ymin, ymax+abs(1e-10*ymax), yscale))
        ax.set_xlabel('time / s')
        ax.set_ylabel(f'{self.x} / mM')
        ax.set_title(f'Detected {self.x} / mM over time / s {self._desc()}')
        plt.show()

    def detectorspectrumplot(self, correction: bool | None = None, duration: float | None = None,
                             lmin: float | None = None, lmax: float | None = None, amin: float | None = None, amax: float | None = None) -> None:
        """
        Plots what the UV-Vis detector sees as a heat map over the course of the entire experiment.

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
        for size, edge, color in [(7, 3, 'black'), (6, 1.5, 'white')]:
            ax.plot(self.points, np.full(len(self.points), self.peak), 'x', markersize=size, markeredgewidth=edge, color=color) # plot each recorded point
        ax.text(0.98, 0.98, f'{self.x0} = {self.c0*1e3:.0f} mM\n{self.details}', transform=ax.transAxes, ha='right', va='top') # reaction conditions

        norm = colors.Normalize(vmin=amin, vmax=amax)
        sm = cm.ScalarMappable(norm=norm, cmap=cmap)
        cbar = plt.colorbar(sm, ax=ax)
        cbar.set_ticks(np.arange(amin, amax+abs(1e-10*amax), 0.1))
        cbar.set_label('absorbance')

        fig.canvas.draw() # calculate ticks preliminarily
        xscale = np.diff(ax.get_xticks())[0]
        yscale = np.diff(ax.get_yticks())[0]
        xmin = np.floor(min(t)/xscale)*xscale

        ax.set_xticks(np.arange(xmin, np.ceil(max(t))+abs(1e-10*max(t)), xscale))
        ax.set_ylim(lmin, lmax)
        ax.set_yticks(np.arange(lmin, lmax+abs(1e-10*lmax), yscale))
        ax.set_xlabel('time / s')
        ax.set_ylabel(r'$\lambda$ / nm')
        ax.set_title(f'UV-Vis spectra for {self.x[1:-1]}{self._textcorr(correction)[0]} over time / s {self._desc()}')
        plt.show()

    @staticmethod
    def flowpoints(tRs: NDArray[np.float64], tbuffer: float, vrxn: float, vfull: float) -> NDArray[np.float64]:
        """
        Calculates the time/s point intervals between flow rate switches according to the flow scheme.

        Args:
            tR: Residence times/s
            tbuffer: Time/s after the solution fully fills the reactor chamber before recording
            vrxn: Reaction volume/ul under LED/film
            vfull: Full reactor chamber volume/ul
        
        Returns:
            points: time/s point intervals between flow rate switches.
        """
        flowrates = vrxn/tRs # flow rates for each tR
        points = np.round(vfull/flowrates + 2*tbuffer)
        return points

    @staticmethod
    def recpoints(tRs: NDArray[np.float64], tbuffer: float, vrxn: float, vfull: float) -> NDArray[np.float64]:
        """
        Calculates the time/s points where an absorbance measurement should be recorded/taken according to the flow scheme,
        cumulative from the start of the experiment.

        Args:
            tR: Residence times/s
            tbuffer: Time/s after the solution fully fills the reactor chamber before recording
            vrxn: Reaction volume/ul under LED/film
            vfull: Full reactor chamber volume/ul
        
        Returns:
            points: time/s points for measurement recording/taking.
        """
        flowrates = vrxn/tRs # flow rates for each tR
        points = []
        prevtime = 0
        for i in range(len(tRs)):
            points.append(round(prevtime + vfull/flowrates[i] + tbuffer))
            prevtime = points[-1] + tbuffer
        return np.array(points)