import numpy as np
from pathlib import Path

"""
Path to the directory where the experimental data is stored
e.g. RFS: Path("Z:/Ocean Optics/yhn24") or relative directory: Path(__file__).parent.parent.parent / 'data'
"""

DATA_DIR = Path(__file__).parent.parent.parent / 'data'

"""
Dictionary of LaTeX text formatting for Matplotlib rendering, add new entries if necessary
"""
X = {'XV': r'[XV$^{+\!\cdot\!}$]',
     'XV0': r'[XV$^{2\!+\!}$]$_0$',
     'EV': r'[EV$^{+\!\cdot\!}$]',
     'EV0': r'[EV$^{2\!+\!}$]$_0$',
     'BV': r'[BV$^{+\!\cdot\!}$]',
     'BV0': r'[BV$^{2\!+\!}$]$_0$',
     'R': r'[Ru(bpy)$_3]^{2\!+\!}$',
     'RU': r'[Ru(bpy)$_3^{2\!+\!}$]',
     'RU0': r'[Ru(bpy)$_3^{2\!+\!}$]$_0$',
     }

"""
Constants
"""

FIGSIZE = (8, 5) # figure size for all plots
DPI = 400 # figure dpi resolution

S120VC = np.pi*(9.5e-3/2)**2 # area/m2 of the thorlabs s120vc sensor of the power meter 
LIGHTBOXBLUE = 105*S120VC # power/W for blue light in the light box (105 W m-2 x s120vc sensor area)
LIGHTBOXUV = 241*S120VC # power/W for uv light in the light box (241 W m-2 x s120vc sensor area)

"""
Flow scheme and flow reactor chamber parameters (tbuffer, vrxn, vfull), 
More detailed definitions and calculations in the report/spreadsheet formulas

tbuffer:  Time/s after the solution fully fills the reactor chamber before recording
vrxn:     Reaction volume/ul under LED/film
vfull:    Full reactor chamber volume/ul
"""

FLOW_PARAMS = {'2x50ul': (30, 50, (2*50)+(2*3.3)+8.4), # (2 x chamber vol) + (2 x luer vol) + tube vol
               '250ul':  (30, 12*8.75*0.8, 250), # film length/mm x chamber width/mm x chamber depth/mm
               }

"""
Substrate parameters to choose from, add new entries for future experiments

x:          LaTeX text formatting for substrate concentration
x0:         LaTeX text formatting for substrate initial concentration
peak:       Substrate absorption peak wavelength/m
e:          Substrate extinction coefficient/M m-1
l:          Path length/m of reactor chamber
bregions:   Wavelength/nm regions to define as baseline
lmin1:      Wavelength/nm minimum for spectrumvstimeplot()
lmax1:      Wavelength/nm maximum for spectrumvstimeplot()
lmin2:      Wavelength/nm minimum for detectorspectrumplot() or spectrumvsconcplot() 
lmax2:      Wavelength/nm maximum for detectorspectrumplot() or spectrumvsconcplot()
"""
#                                x        x0        peak e       l       bregions                           lmin1 lmax1 lmin2 lmax2
SUBSTRATE = {'EVrubpybatch':    (X['EV'], X['EV0'], 600, 1.22e6, 1e-2,   [(300,330), (430,450), (770,780)], 300,  800,  450,  700),
             'BVrubpybatch':    (X['BV'], X['BV0'], 550, 0.78e6, 1e-2,   [(300,330), (430,450), (770,780)], 300,  800,  450,  700),
             'EVfluorbatch':    (X['EV'], X['EV0'], 600, 1.22e6, 1e-2,   None,                              300,  800,  450,  700),
             'BVfluorbatch':    (X['BV'], X['BV0'], 550, 0.78e6, 1e-2,   None,                              300,  800,  450,  700),
             'EVrubpystopflow': (X['EV'], X['EV0'], 600, 1.22e6, 700e-6, [(300,330), (430,450), (770,780)], 300,  800,  450,  700),
             'BVrubpystopflow': (X['BV'], X['BV0'], 540, 0.78e6, 700e-6, [(300,330), (430,450), (770,780)], 300,  800,  450,  700),
             'EVrubpy2x50ul':   (X['EV'], X['EV0'], 600, 1.22e6, 700e-6, [(300,330), (430,450), (770,780)], 300,  800,  450,  700),
             'Rubpy':           (X['RU'], X['RU0'], 453, 1.48e6, 1e-2,   [(150,200), (600,1000)],           200,  1000, 340,  500),
             'EVrubpy':         (X['EV'], X['EV0'], 600, 1.22e6, 800e-6, [(200,250), (770,780)],            200,  1000, 500,  660),
             'BVrubpy':         (X['BV'], X['BV0'], 540, 0.78e6, 800e-6, [(200,250), (770,780)],            200,  1000, 500,  660),
             'EVxanthenegreen': (X['EV'], X['EV0'], 600, 1.22e6, 800e-6, [(200,250), (750,1000)],           200,  1000, 580,  660),
             'EVxantheneblue':  (X['EV'], X['EV0'], 600, 1.22e6, 800e-6, [(200,250), (770,780)],            200,  1000, 580,  660),
             'EVxantheneuv':    (X['EV'], X['EV0'], 600, 1.22e6, 800e-6, [(200,250), (750,1000)],           200,  1000, 580,  660),
             }