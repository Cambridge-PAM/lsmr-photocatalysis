# LSMR Photocatalysis

Python package and Jupyter notebooks for analysing Luminescent Solar Microreactor (LSMR) photocatalysis results and an automated pump control script.

## Features
- Data extraction from OceanOptics/OceanView files
- UV-Vis spectra baseline correction
- Linear regression to find apparent quantum yield and reaction rate
- UV-Vis spectra visualisation over experiment time

## File structure

```
lsmr-photocatalysis/
|
├── data/ # manually create this folder in your local repo (optional, ignored by git) or host the data files elsewhere
|   ├── date_PC_XV_SED_Film_LED_Batch/
|   ├── date_PC_XV_SED_Film_LED_StoppedFlow/
|   └── date_PC_XV_SED_Film_LED_Flow/
|       ├── data1.txt
|       ├── data2.txt
|       └── ...
|
├── notebooks/ # create new folders for new experiments and store the Jupyter notebooks here for data analysis
|   ├── bargraphs/
|   ├── eosiny/
|   ├── fluorescein/
|   └── rubpy/
|       ├── date_flow_led_film_pc_xv_sed_1.ipynb
|       ├── date_flow_led_film_pc_xv_sed_2.ipynb
|       └── ...
|
├── scripts/
|   └── pump.py # pump control script for automatic flow rate changes in continuous flow experiments
|
├── src/
|   └── lsmr_photocatalysis/ # package of classes and functions to import into Jupyter notebooks during data analysis
|       ├── __init__.py
|       ├── batch.py
|       ├── continuousflow.py
|       ├── experiment.py
|       ├── settings.py # change directory to data files and add other experimental parameters/constants here 
|       └── stoppedflow.py
|
├── .gitignore
├── environment.yml
├── pyproject.toml
└── README.md
```

## Initial setup

Install Anaconda/Miniconda (Python), Git and VS Code if not already done so.

Clone the repository to a local folder on your computer using Git Bash.

```bash
git clone https://github.com/Cambridge-PAM/lsmr-photocatalysis.git
cd lsmr-photocatalysis
```
Create and activate the Python environment.

```sh
conda env create -f environment.yml
conda activate lsmr-photocatalysis
```

In `./src/lsmr_photocatalysis/settings.py`, change `DATA_DIR` to the local directory where the data files are stored, such as on RFS (performance will be slower on a network drive),

```python
DATA_DIR = Path("Z:/Ocean Optics/yhn24")
```

or as a relative directory by creating a new folder `./data` with the data files in the local repository itself (this is the default behaviour).

```python
DATA_DIR = Path(__file__).parent.parent.parent / 'data'
```
## General workflow

- Copy the raw data files from RFS into the directory chosen during setup (e.g. `./data`).
- Create a Jupyter notebook in `./notebooks` for a new experiment.
- Refer to the contents of an exisiting notebook as an example. But generally,
    - Import the appropropriate class from the `lsmr_photocatalysis` package, i.e.
        - `Batch`
        - `StoppedFlow`
        - `ContinuousFlow`
    - Create an instance of the class and specify its attributes (refer to the files in `./src/lsmr_photocatalysis` for more information), e.g.
        - `data` Data folder name
        - `c0`: Initial concentration/M of substrate of interest
        - `details`: Details of other reagents and concentrations
        - `x`: Dictionary key for `SUBSTRATE` properties, see `settings.py`
        - `film`: Polymer film name
        - `light`: LED colour
        - `wavelength`: LED wavelength/m
        - `power`: LED power/W
        - `correction`: Default choice for plots to be baseline corrected or not. Defaults to `True`
        - etc.
    - Call methods on that instance to plot the desired graphs, e.g.
        - `concvstimeplot()`
        - `detectorconcplot()`
        - `detectorspectrumplot()`
        - `spectrumvstimeplot()`
        - `spectrumvsconcplot()`
- Add additional constants/parameters in `./src/lsmr_photocatalysis/settings.py` for new experiments when required.

## Pump control script
- Install and setup the Python environment and Git repository on the lab computer (refer to [Initial setup](#initial-setup)).
- Navigate to `./scripts` and open `pump.py` in VS Code.
- Check "Device Manager" on the lab computer for the COM port that the pump is connected to and edit the `com` variable accordingly.
- Edit the `tRs` numpy array of residence times for the current experiment.
- Add a new dictonary entry in `FLOW_PARAMS` in `./src/lsmr_photocatalysis/settings.py` if a new continuous flow setup has been created.
- Run the `pump.py` script in VS Code or in the terminal when ready.

## Contact

Yu Hung Ng  
St Catharine's College, Cambridge  
yhn24@cam.ac.uk  
brandon.ngyuhung@gmail.com