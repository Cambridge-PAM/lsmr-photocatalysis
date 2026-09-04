from .settings import DATA_DIR, SUBSTRATE, S120VC
import numpy as np
import scipy.stats as stats
import scipy.constants as constants
from pybaselines import Baseline
from typing import TextIO
from numpy.typing import NDArray

class _Experiment:
    """
    Private class outlining attributes for a generic photocatalysis experiment with several private utility methods.

    Args:
        data: Data folder name.
        c0: Initial concentration/M of substrate of interest.
        details: Details of other reagents and concentrations.
        x: Dictionary key for `SUBSTRATE` properties, see `settings.py`.
        film: Polymer film name.
        light: LED colour.
        wavelength: LED wavelength/m.
        power: LED power/W.
        correction: Default choice for plots to be baseline corrected or not. Defaults to `True`.
        duration: How long to plot for `detectorconcplot()` and `detectorspectrumplot()` in `ContinuousFlow` 
                  and `concvstime()` and `spectrumvstime()` in `Batch`, meaningless for `StoppedFlow`. Defaults to `1e10`.
    
    Attributes:
        points: Time/s points to extract data for `ContinuousFlow`, `None` for `StoppedFlow` and `Batch`.
        start: Reaction start timestamp in `'hh-mm-ss-msmsms'` for `ContinuousFlow` and `Batch`, meaningless `''` for `StoppedFlow`.
        x: LaTeX text formatting for substrate concentration.
        x0: LaTeX text formatting for substrate initial concentration.
        peak: Substrate absorption peak wavelength/m.
        e: Substrate extinction coefficient/M m-1.
        l: Path length/m of reactor chamber.
        bregions: Wavelength/nm regions to define as baseline.
        lmin1: Wavelength/nm minimum for `spectrumvstimeplot()`.
        lmax1: Wavelength/nm maximum for `spectrumvstimeplot()`.
        lmin2: Wavelength/nm minimum for `detectorspectrumplot()` in `ContinuousFlow` 
               or `spectrumvsconcplot()` in `Batch`, meaningless `''` for `StoppedFlow`.
        lmax2: Wavelength/nm maximum for `detectorspectrumplot()` in `ContinuousFlow` 
               or `spectrumvsconcplot()` in `Batch`, meaningless `''` for `StoppedFlow`.
    """

    ttext = '$t$' # reaction time

    def __init__(self, data: str, c0: float, details: str, x: str, film: str, light: str, wavelength: float, power: float,
                 correction: bool = True, duration: float = 1e10):
    
        self.data: str = data 
        self.points: np.typing.NDArray[np.float64] | None = None
        self.start: str = ''
        self.c0: str = c0 
        self.details: str = details 
        self.x, self.x0, self.peak, self.e, self.l, self.bregions, self.lmin1, self.lmax1, self.lmin2, self.lmax2 = SUBSTRATE[x] # substrate properties
        self.film: str = film
        self.light: str = light
        self.wavelength: float = wavelength
        self.intensity: float = power/S120VC 
        self.correction: bool = correction 
        self.duration: float = duration

    def _desc(self):
        """Returns a formatted description of the photocatalysis light setup to be put in the plot title."""
        if self.film == 'no':
            article = ''
        else:
            article = 'a'
        return f'with {article}\n{self.film} film under {self.light} light ({int(self.wavelength*1e9)} nm, {round(self.intensity)} W m$^{{-2}}$)'

    def _ratephotons(self):
        """Return the rate of photons (M s-1) for a given LED intensity/W m-2, wavelength/m, path length/m."""
        h = constants.physical_constants['Planck constant'][0]
        c = constants.physical_constants['speed of light in vacuum'][0]
        n = constants.physical_constants['Avogadro constant'][0]
        return self.intensity/(self.l*1e3*n*h*c/self.wavelength)

    @staticmethod
    def _numerr(num: float, err: float, prefix: int = 1, sci: bool = True) -> str:
        """Returns a formatted string of `num` ± `err` with the correct dp in scientific notation with optional unit prefix specifications."""

        num = num*(1/prefix) if np.isfinite(num) else 0 # scale by the unit prefix
        err = err*(1/prefix) if np.isfinite(err) else 0

        i = int(f'{err:.0e}'.split('e')[-1])+1 # find the exponent such that the error has 1sf in the first decimal place

        if sci and prefix==1: # scientific notation to 1dp by default
            return rf'${num*(10**-i):.1f}$ ± ${err*(10**-i):.1f}$ $\times$ 10$^{{{i}}}$'
        else: # use more dp (implied default if a unit prefix is specified)
            return f'${num:.{max(-i+1,0)}f}$ ± ${err:.{max(-i+1,0)}f}$'

    @staticmethod
    def _textalign(slope: float) -> tuple[float, float, str, str, float, float, str, str]:
        """Returns parameters for aligning textboxes according to the slope direction for `concvstimeplot()`."""

        if slope > 2e-10:
            a1, b1, ha1, va1 = 0.98, 0.03, 'right', 'bottom'
            a2, b2, ha2, va2 = 0.02, 0.97, 'left', 'top'
        else:
            a1, b1, ha1, va1 = 0.98, 0.97, 'right', 'top'
            a2, b2, ha2, va2 = 0.02, 0.03, 'left', 'bottom'
        return a1, b1, ha1, va1, a2, b2, ha2, va2

    def _textcorr(self, correction: bool) -> tuple[str, str]:
        """Returns details of baseline parameters as a text description to be put in the title and textbox of a plot."""

        if correction is None:
            correction = self.correction
        if correction:
            params = 'baseline set at\n' # baseline interpolated between set regions
            for start, end in self.bregions:
                params += f'{start}\u2013{end} nm,\n'
            params = params[:-1] + '\nthen interpolated'
            return ' (baseline corrected)', params
        else:
            return ' (raw data)', ''

    def _elapsed(self, timestamp: str) -> float:
        """Returns the elapsed time/s from the start of experiment for a timestamp in the format `'hh-mm-ss-msmsms'`."""

        def parse(ts):
            parts = ts.split('-')
            if len(parts) != 4:
                raise ValueError('Timestamps must use the format "hh-mm-ss-msmsms".')

            h, m, s, ms = (int(part) for part in parts)
            if not (0 <= m < 60 and 0 <= s < 60 and 0 <= ms < 1000):
                raise ValueError('Minutes, seconds, and milliseconds must be valid time units.')

            return h * 3600 + m * 60 + s + ms * 0.001

        return parse(timestamp) - parse(self.start)

    def _readspectrum(self, file: TextIO, correction: bool = True) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
        """
        Reads OceanView data file, extracts the UV-Vis absorption spectrum and an produces an interpolated baseline.
        
        Args:
            file: `TextIO` items corresponding to OceanView data files after a `f.open('r') as file:`.
            correction: Applies baseline correction to the absorbances if `True`.
        
        Returns:
            l: Wavelengths/nm.
            a: Absorbances.
            b: Baseline.
        """
        
        l = [] # wavelengths
        a = [] # absorbances

        for line in file:
            try: # line contains spectrum data
                l.append(float(line.split('	')[0]))
                a.append(float(line.split('	')[1]))
            except: # line does not contain spectrum data
                continue

        l = np.array(l)
        a = np.array(a)

        # baseline correction
        if self.bregions is not None:
            mask = np.zeros(len(a), dtype=bool) # new mask for baseline regions
            ba = np.array([]) # baseline absorbances
            for i, (start, end) in enumerate(self.bregions):

                tempmask = np.zeros(len(a), dtype=bool) # new temporary mask for averaging regions of a
                tempmask |= (start <= l) & (l <= end) # finding the region
                mask |= tempmask # apply to overall mask also
                ave = np.mean(a[tempmask]) # average absorbance in the region
                length = len(a[tempmask])

                if i == 0:
                    tempmask = np.zeros(len(a), dtype=bool) # new temporary mask for averaging regions of a
                    tempmask |= (l < start) # extend the first baseline all the way before the start
                    mask |= tempmask # apply to overall mask also
                    ba = np.concatenate((ba, ave*np.ones(len(a[tempmask])))) # appending ave absorbances to ba

                ba = np.concatenate((ba, ave*np.ones(length))) # appending ave absorbances to ba
                
                if i == len(self.bregions)-1:
                    tempmask = np.zeros(len(a), dtype=bool) # new temporary mask for averaging regions of a
                    tempmask |= (end < l) # extend the last baseline all the way after the end
                    mask |= tempmask # apply to overall mask also
                    ba = np.concatenate((ba, ave*np.ones(len(a[tempmask])))) # appending ave absorbances to ba
            baseline = Baseline(l)
            b, _ = baseline.interp_pts(l, baseline_points=np.column_stack((l[mask], ba)), interp_method='quadratic') # baseline from interpolating between bregions
        else:
            b = np.zeros(len(a))

        if correction:
            a = a - b # apply baseline correction to the absorbances
        
        return l, a, b

    def _readpeak(self, wavelengths: NDArray[np.float64], absorbances: NDArray[np.float64], pwidth: int) -> float:
        """Extracts average absorbance at peak wavelength and within `pwidth`"""
        mask = (self.peak-pwidth/2 <= wavelengths) & (wavelengths <= self.peak+pwidth/2) # mask for wavelengths within peak width
        return np.mean(absorbances[mask]) # average absorbance within peak width

    def _concvstime(self, correction: bool | None = None, duration: float | None = None, pwidth: int = 10) -> tuple[np.typing.NDArray[np.float64], np.typing.NDArray[np.float64]]:
        """
        Return the concentration/M of the substrate over time/s from OceanView data.

        Args:
            correction: Whether to apply baseline correction. Defaults to `None` which then points to the instance's attribute.
            duration: Limits the data to a certain time duration/s. Defaults to `None` which then points to the instance's attribute.
            pwidth: Wavelength/nm window around the substrate absorbance peak for absorbance value averaging. Defaults to `10`.

        Returns:
            t: Times/s.
            c: Concentrations/M.
        """

        if correction is None:
            correction = self.correction
        if duration is None:
            duration = self.duration

        a = [] # array of peak absorbances
        t = [] # array of time in seconds
        read = False

        for f in sorted((DATA_DIR / self.data).glob('*.txt'), key=lambda f: f.stem[-12:]): # sort files by timestamp
            timestamp = f.stem[-12:]

            if self.start == '': # if no start time defined, start with the first file
                self.start = timestamp

            time = self._elapsed(timestamp) # convert timestamp to time in s

            if time >= duration:
                read = False # stop reading files
                break
            elif time >= 0:
                read = True # start reading files
            
            if read:
                t.append(time) # list of time points
                with f.open('r') as file:
                    wavelengths, absorbances, _ = self._readspectrum(file, correction=correction)
                    a.append(self._readpeak(wavelengths, absorbances, pwidth))

        c = np.array(a)/(self.e*self.l) # convert absorbance to concentration using beer-lambert law
        return np.array(t), c

    def _concvstimepoints(self, correction: bool | None = None, pwidth: int = 10, twidth: int = 5) -> NDArray[np.float64]:
        """
        Return the concentration/M of the substrate at given time/s points from OceanView data.

        Args:
            correction: Whether to apply baseline correction. Defaults to `None` which then points to the instance's attribute.
            pwidth: Wavelength window/nm around the substrate absorbance peak for absorbance value averaging. Defaults to `10`.
            twidth: Time/s window around each time points to average the concentrations. Defaults to `5`.

        Returns:
            c: Concentrations/M.
        """

        if correction is None:
            correction = self.correction
        
        a = [] # array of peak absorbances
        files = sorted((DATA_DIR / self.data).glob('*.txt'), key=lambda f: f.stem[-12:]) # sort files by timestamp
        
        for i in range(len(files)):

            if self.points is not None: # read files from certain time points (in seconds)
                timestamp = files[i].stem[-12:]
                
                if self.start == '':  # if no start time defined, start with the first file
                    self.start = timestamp

                time = round(self._elapsed(timestamp)) # convert timestamp to time in s, round to integer

                if len(a) >= len(self.points): # if all points have been read, break
                    break

                if time >= self.points[len(a)]: # if time is more than or equal to the next point to be searched for
                    atemp = [] # temporary list to store absorbances within twidth seconds of timestamp
                    for j in range(-twidth//2 + 1, twidth//2 + 1): # search for absorbance peak in files within twidth seconds of timestamp
                        with files[i+j].open('r') as file:
                            wavelengths, absorbances, _ = self._readspectrum(file, correction=correction)
                            atemp.append(self._readpeak(wavelengths, absorbances, pwidth))
                    a.append(np.mean(np.array(atemp))) # add the averaged absorbances

            else: # read all files
                with files[i].open('r') as file:
                    wavelengths, absorbances, _ = self._readspectrum(file, correction=correction)
                    a.append(self._readpeak(wavelengths, absorbances, pwidth))

        c = np.array(a)/(self.e*self.l)  # convert absorbance to concentration using beer-lambert law
        return c

    # uv-vis absorbance spectrum over time from oceanoptics data
    def _spectrumvstime(self, correction: bool | None = None, duration: float | None = None) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
        """
        Return the UV-Vis absorption spectrum over time/s from OceanView data

        Args:
            correction: Whether to apply baseline correction. Defaults to `None` which then points to the instance's attribute.
            duration: Limits the data to a certain time/s duration. Defaults to `None` which then points to the instance's attribute.

        Returns:
            t: Time/s.
            l: Wavelengths/nm.
            a: Absorbances for each time in t (2D array).
        """
        if correction is None:
            correction = self.correction
        if duration is None:
            duration = self.duration

        t = [] # time in seconds
        a = [] # absorbances
        read = False

        for f in sorted((DATA_DIR / self.data).glob('*.txt'), key=lambda f: f.stem[-12:]): # sort files by timestamp
            timestamp = f.stem[-12:]

            if self.start == '': # if no start time defined, start with the first file
                self.start = timestamp
            
            time = self._elapsed(timestamp) # convert timestamp to time in s

            if time >= duration:
                read = False # stop reading files
                break
            elif time >= 0:
                read = True # start reading files

            if read:
                t.append(time) # list of time points
                with f.open('r') as file:
                    wavelengths, absorbances, _ = self._readspectrum(file, correction=correction)
                    if len(a) == 0: # if first file
                        l = wavelengths # create list of wavelengths
                    a.append(absorbances) # append absorbances for this time point

        return np.array(t), l, np.array(a).T # convert lists to numpy array and transpose a to have wavelengths as rows and time points as columns

    # uv-vis absorbance spectrum at given time points (in seconds) from oceanoptics data
    def _spectrumvstimepoints(self, points: NDArray[np.float64] | None = None, correction: bool | None = None) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
        """
        Return the UV-Vis absorption spectrum given time/s points from OceanView data

        Args:
            correction: Whether to apply baseline correction. Defaults to `None` which then points to the instance's attribute.
            points: Time/s points to obtain the spectrum. Defaults to `None` which then points to the instance's attribute.

        Returns:
            l: Wavelengths/nm.
            a: Absorbances.
            b: Baseline.
        """
        
        if correction is None:
            correction = self.correction
        if points is None:
            points = self.points

        a = [] # absorbances
        b = [] # baselines
        
        for f in sorted((DATA_DIR / self.data).glob('*.txt'), key=lambda f: f.stem[-12:]): # sort files by timestamp
            
            if points is not None: # read files from certain time points (in seconds)
                timestamp = f.stem[-12:]

                if self.start == '': # if no start time defined, start with the first file
                    self.start = timestamp

                time = round(self._elapsed(timestamp)) # convert timestamp to time in s, round to integer
                
                if len(a) >= len(points): # if all points have been read, break
                    break

                if time >= points[len(a)]: # if time is more than or equal to the next point to be searched for
                    pass # continue on to read the file
                else:
                    continue # skip this file

            with f.open('r') as file: # arrive here to read the file if time matches the specified points or read all files
                wavelengths, absorbances, baseline = self._readspectrum(file, correction=correction)
                if len(a) == 0: # if first file
                    l = wavelengths # create list of wavelengths
                a.append(absorbances) # append absorbances for this time point
                b.append(baseline)

        return l, np.asarray(a), np.asarray(b)

    @staticmethod
    def _linreg(t: NDArray[np.float64], c: NDArray[np.float64], 
                exclude: list[int] = [], istart: int = 0, iend: int | None = None, trial: int = 1, 
                r2threshold: float | None = None, interceptthreshold: float | None = None) -> tuple[None, int, int]:
        """
        Linear regression analysis, two possible modes:
        1. Regression of points from `istart` to `iend`, excluding points in specified in `exclude` (e.g. outliers).
        2. Iterative regression trials from `istart` by `trial` no. points, then `trial+1`, ..., until `iend` or `r2threshold` or `interceptthreshold` is exceeded.

        Args:
            t: Time points/s
            c: Concentrations/M or ln(concentrations/M)
            exclude: List of point indexes to exclude. Defaults to `[]`.
            istart: Point index to start the regression from. Defaults to `0`.
            iend: Point index to end the regression at. Defaults to `None` which then points to the last point in the data.
            trial: Number of points after istart to begin iterative regression trials. Defaults to `1`.
            r2threshold: Threshold for coefficient of determination r2. Defaults to `None`.
            interceptthreshold: Threshold for y-axis intercept. Defaults to `None`.

        Returns:
            fit: SciPy `LingressResult` object with attributes `slope`, `stderr`, `intercept`, `intercept_stderr`, `rvalue`, etc.
            imin: Index of first point involved in the regression, for plotting use later.
            imax: Index of last point involved in the regression, for plotting use later.
        """

        if iend == None: # by default iend is the last point available
            iend = len(t)-1

        # create mask to remove points outside the range of interest
        mask = np.ones(len(t), dtype=bool)
        mask[:istart] = False
        mask[iend+1:] = False
        if exclude: # linear regression excluding some outlier points specified in exclude
            mask[exclude] = False

        if r2threshold is not None or interceptthreshold is not None: # run regression until r2 or intercept exceeds threshold
            stop = False
            for i in range(trial+istart, len(t)):
                
                if i+1 < len(mask):
                    if mask[i] == False or mask[i+1] == False: 
                        continue # skip if the current point and/or the next point is excluded

                tempmask = mask.copy() # create independent copy of mask
                tempmask[i+2:] = False # only include points up to i+1
                fit = stats.linregress(t[tempmask], c[tempmask]) # linear regression on the points within the mask

                if r2threshold is not None:
                    if fit.rvalue**2 < r2threshold:
                        stop = True
                if interceptthreshold is not None:
                    if abs(fit.intercept) > abs(interceptthreshold):
                        stop = True

                if stop: # stop when r2 or intercept exceeds thresholds
                    tempmask = mask.copy() # create independent copy of mask
                    tempmask[i+1:] = False # only include points up to i, i.e. the previous fit
                    mask = tempmask # update mask
                    break

        fit = stats.linregress(t[mask], c[mask]) # linear regression on the points within the mask
        imin = np.argmax(mask) # first point included in the regression
        imax = len(mask) - np.argmax(mask[::-1]) - 1 # last point included in the regression
        
        return fit, imin, imax