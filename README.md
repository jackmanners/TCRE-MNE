# EDF Processing and Topo Map Plotting
This is a basic script designed to assist with the processing of EDF (European Data Format) files and the generation of topographic maps from power spectral density data. The script is specifically set up to work with TCRE (Tipolar Concentric Ring Electrodes) and eEEG (Emulated Electroencephalography) data captured via the Compumedics [Grael 4K PSG:EEG](https://www.compumedics.com.au/en/products/grael-4k-psg-eeg/), but it should also be functional with and `.edf` files.

The main processing is performed using the [MEG + EEG Analysis & Visualisation](https://mne.tools/stable/index.html) Python package, which provides powerful tools for EEG and MEG (Magnetoencephalography) data analysis.

### Requirements
This script uses Python 3.10.5 and relies on the `mne` package, which can be installed using `pip install mne`.

# Usage

The functions can be found in the `process_plot_edf.py` file. Alternatively, if the required package is installed, the functions can be directly called from the command line. Details are provided below, and additional help can be obtained by using the command `python process_plot_edf.py --help`.

**Process EDF Files and Plot Topo Maps**
- Command: `topo_plot_edfs`
- This command processes EDF files and generates topographic maps from the data.
- Required Arguments:
    - `--edf-path <EDF-PATH>`: Path to the EDF file or folder. Subfolders will be searched.
- Optional Arguments:
    - `--output-path <OUTPUT-PATH>`: Path to save the plot(s).
        - Default: The current directory.
    - `--channel-prefix <CHANNEL-PREFIX>`: EDF channel prefix ('e': EEG, 't': TCRE). Multiple prefixes can be provided.
        - Default: None.
    - `--duration <DURATION>`: Duration in seconds.
        - Default: The entire recording.
    - `-h, --help`: Help

Example Usage:
```
python process_plot_edf.py process_edfs --edf-path 'edfs/PX01 - N1.edf' --channel-prefix e t --duration 180
python process_plot_TCRE_edf.py topo_plot_edfs --edf-path 'C:\TCRE\EDF - Sleep' --channel-prefix e t
```

**Combine Topo Maps**
- Command: `combine_topo_plots`
- This command combines multiple topographic maps into a single plot.
- Required Arguments:
    - `--figure-path <FIGURE-PATH>`: Path to the folder containing the individual (.png) topo maps.
- Optional Arguments:
    - `--output-path <OUTPUT-PATH>`: Path to save the combined plot(s).
        - Default: The current directory.
    - `--pxs <PXS>`: Single ID or list of IDs to filter figures by (e.g., `<PX01>` or `<PX01 PX02 PX03>`).
        - Default: Combine all files and save as `stacked_plot.png`.
    - `-h, --help`: Help

Example Usage:
```
python process_plot_edf.py combine_topo_plots --figure-path figures --output-path figures/combined --pxs PX01 PX02 PX03
```