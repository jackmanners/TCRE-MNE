import argparse
import glob
import os
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
import mne
from matplotlib import gridspec

# plt.switch_backend('agg')

def plot_topo_from_edf(edf_path, plt_savepath=None, duration=None, channel_prefix=None):
    """
    Plots topographic maps from EEG data stored in an EDF file.

    Parameters:
        edf_path (str): Path to the EDF file.
        channel_prefix (str, optional): Prefix of the EEG channels to select. Default is None.
        plt_savepath (str, optional): Directory path to save the generated plot. Default is None (CWD).
        duration (float, optional): Duration of the data to plot in seconds. Default is None (plots the entire duration).

    Returns:
        None
    """

    # Read the EDF file and load the data
    raw = mne.io.read_raw_edf(
        edf_path,
        eog=['E1', 'E2'],
        misc=['EMG', 'EMG2', 'SpO2', 'Pulse', 'Derived HR', 'ManPosition', 'ManLights', 'ECG'],
        infer_types=True,
        preload=True)

    # Select EEG channels only
    raw.pick(['eeg'])

    # Crop the data if a specific duration is provided
    if duration:
        raw.crop(tmax=duration)

    # Select EEG channels based on the provided prefix and rename them
    if channel_prefix:
        eeg_channels = [ch_name for ch_name in raw.ch_names if ch_name.startswith(channel_prefix)]
        raw.pick_channels(eeg_channels)
        raw.rename_channels(lambda x: x[1:])

    # Rename specific channels
    channel_rename_dict = {'FpZ': 'Fpz', 'CpZ': 'CPz'}
    raw.rename_channels(channel_rename_dict)

    # Set the montage to standard 10-20
    montage = mne.channels.make_standard_montage('standard_1020')
    raw.set_montage(montage)

    # Compute power spectral density and plot topographic map
    spectrum = raw.compute_psd()
    cmap = mpl.colormaps['Reds']
    fig = spectrum.plot_topomap(bands = {
        'Delta (0.5-4.5 Hz)': (0.5, 4.5),
        'Theta (4.5-8 Hz)': (4.5, 8),
        'Alpha (8-12 Hz)': (8, 12),
        'Sigma (12-15 Hz)': (12, 15),
        'Beta (15-32 Hz)': (12, 30)
        }, res=300, size=2, sphere=0.1, 
        sensors=True, dB=True, extrapolate='head',
        ch_type='eeg', cmap=(cmap, False), show=False)
    # Adjust figure margins
    fig.subplots_adjust(left=0.05)

    # Generate save name based on channel prefix
    plt_savename = edf_path.split('/')[-1].split('.')[0]
    if channel_prefix == 'e':
        savename = f"{plt_savename} (EEG)"
    elif channel_prefix == 't':
        savename = f"{plt_savename} (TCRE)"
    else:
        savename = plt_savename

    # Save the figure and close it
    fig.tight_layout()
    if plt_savepath: savename = f"{plt_savepath}/{savename}"
    fig.savefig(savename, bbox_inches='tight', dpi=300)
    plt.close(fig)


def combine_topo_plots(filenames, save_path=None, dpi=3000):
    """
    Combines multiple topographic plot images into a single figure for EEG and TCRE plots.

    Parameters:
        filenames (list): List of file names of the topographic plot images.
        save_path (str): File path to save the combined figure.
        dpi (int, optional): Dots per inch for the saved figure. Default is 500.

    Returns:
        None
    """

    # Separate EEG and TCRE filenames
    eeg_filenames = [filename for filename in filenames if '(EEG)' in filename]
    tcre_filenames = [filename for filename in filenames if '(TCRE)' in filename]

    # Extract EEG and TCRE stages from filenames
    eeg_stages = [re.findall(r' - (\w+)', filename)[0] for filename in eeg_filenames]
    tcre_stages = [re.findall(r' - (\w+)', filename)[0] for filename in tcre_filenames]

    # Finds unique stages and orders these appropriately
    stage_order = ['N1', 'N2', 'N3', 'REM']
    unique_stages = list(set(eeg_stages + tcre_stages))
    unique_stages = sorted(unique_stages, key=lambda x: stage_order.index(x) if x in stage_order else float('inf'))
    num_rows = max(len(unique_stages), 1)

    # Create a GridSpec object
    fig = plt.figure(figsize=(4, num_rows/2 * 0.5))
    gs = gridspec.GridSpec(num_rows, 3, width_ratios=[20, 0.5, 20], figure=fig, wspace=0, hspace=0.01)

    # Plot EEG images in the first column
    for i, (filename, stage) in enumerate(zip(eeg_filenames, eeg_stages)):
        ax = plt.subplot(gs[i, 0])
        image = plt.imread(filename)
        ax.imshow(image)
        ax.axis('off')
        if i == 0: 
            title = ax.set_title('EEG')
            title.set_fontsize(6)

    # Add sleep-stage headers in the second column
    for i, stage in enumerate(unique_stages):
        ax = plt.subplot(gs[i, 1])
        ax.text(0.5, 0.5, stage, ha='center', va='center', fontsize=4)
        ax.axis('off')

    # Plot TCRE images in the third column
    for i, (filename, stage) in enumerate(zip(tcre_filenames, tcre_stages)):
        ax = plt.subplot(gs[i, 2])
        image = plt.imread(filename)
        ax.imshow(image)
        ax.axis('off')
        if i == 0: 
            title = ax.set_title('TCRE')
            title.set_fontsize(6)

    if not save_path: save_path = 'stacked_plot'

    # Adjust spacing between subplots, save, and close figure
    plt.subplots_adjust(wspace=0.001)
    plt.tight_layout()
    plt.savefig(save_path, bbox_inches='tight', dpi=dpi)
    plt.close()
    
    
def plot_spect_from_edf(edf_paths, plt_savepath=None, duration=None, channel_prefix=None):
    """
    Plots spectrum density plots from EEG data stored in an EDF file.

    Parameters:
        edf_path (str): Path to the EDF file.
        channel_prefix (str, optional): Prefix of the EEG channels to select. Default is None.
        plt_savepath (str, optional): Directory path to save the generated plot. Default is None (CWD).
        duration (float, optional): Duration of the data to plot in seconds. Default is None (plots the entire duration).

    Returns:
        None
    """

    fig, ax = plt.subplots()

    for edf_path in edf_paths:
        # Read the EDF file and load the data
        raw = mne.io.read_raw_edf(
            edf_path,
            eog=['E1', 'E2'],
            misc=['EMG', 'EMG2', 'SpO2', 'Pulse', 'Derived HR', 'ManPosition', 'ManLights', 'ECG'],
            infer_types=True,
            preload=True)

        # Select EEG channels only
        raw.pick(['eeg'])

        # Crop the data if a specific duration is provided
        if duration:
            raw.crop(tmax=duration)

        # Select EEG channels based on the provided prefix and rename them
        if channel_prefix:
            eeg_channels = [ch_name for ch_name in raw.ch_names if ch_name.startswith(channel_prefix)]
            raw.pick_channels(eeg_channels)
            raw.rename_channels(lambda x: x[1:])

        # Rename specific channels
        channel_rename_dict = {'FpZ': 'Fpz', 'CpZ': 'CPz'}
        raw.rename_channels(channel_rename_dict)

        # Set the montage to standard 10-20
        montage = mne.channels.make_standard_montage('standard_1020')
        raw.set_montage(montage)

        # Compute power spectral density
        spectrum = raw.compute_psd()

        # Plot the spectrum
        spectrum.plot(
            dB=True,
            average=True,
            axes=ax,
            show=False
        )

    # Generate save name based on channel prefix
    plt_savename = edf_paths[0].split('/')[-1].split('.')[0]
    plt_savename = edf_paths[0].split('\\')[-1].split('.')[0]
    if channel_prefix == 'e':
        savename = f"{plt_savename} (EEG)"
    elif channel_prefix == 't':
        savename = f"{plt_savename} (TCRE)"
    else:
        savename = plt_savename

    # Save the figure and close it
    fig.tight_layout()
    if plt_savepath:
        savename = os.path.join(plt_savepath, savename)
    fig.savefig(savename, bbox_inches='tight', dpi=300)
    plt.close(fig)


def parse_cli_args():
    parser = argparse.ArgumentParser(description='Process EDF files and plot topographic PSD brain-maps',
                                     formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers(dest='command')

    # Process EDFs sub-command
    process_parser = subparsers.add_parser(
        'topo_plot_edfs',
        help='Process EDF files and plot topo maps\n'
             'Run with "python process_plot_edf.py topo_plot_edfs --edf-path <EDF-PATH>"\n'
             'Use command "python process_plot_edf.py topo_plot_edfs --help" for more information & options\n'
    )
    process_parser.add_argument('--edf-path', metavar='<edf_path>', required=True,
                                help='Path to the EDF file or folder\n'
                                '')
    process_parser.add_argument('--output-path', metavar='<output_path>', default=None, help='Directory to save plots in')
    process_parser.add_argument('--channel-prefix', metavar='<channel_prefix>', default=None, nargs='+', choices=['e', 't', 'e t'], help='Channel prefix ("e", "t", or "e t" for both)')
    process_parser.add_argument('--duration', metavar='<duration>', type=int, default=None, help='Duration in seconds')

    # Combine topo maps sub-command
    combine_parser = subparsers.add_parser(
        'combine_topo_plots',
        help='Combine multiple topo maps into a single plot\n'
             'Run with "python process_plot_edf.py combine_topo_plots --figure-folder <FIGURE-FOLDER>"\n'
             'Use command "python process_plot_edf.py combine_topo_plots --help" for more information & options\n'
    )
    combine_parser.add_argument('--figure-folder', metavar='<figure_folder>', required=True, help='Initial figure directory')
    combine_parser.add_argument('--output-path', metavar='<output_path>', default=None, help='Directory to save the combined plot(s)')
    combine_parser.add_argument('--pxs', metavar='<pxs>', default=None, nargs='+', help='Single or list of px IDs (e.g., <PX01> or <PX01 PX02 PX03>)')
    
    # Spectral plot EDFs sub-command
    process_parser = subparsers.add_parser(
        'spect_plot_edfs',
        help='Process EDF files and plot spectral density maps\n'
             'Run with "python process_plot_edf.py spect_plot_edfs --edf-path <EDF-PATH>"\n'
             'Use command "python process_plot_edf.py spect_plot_edfs --help" for more information & options\n'
    )
    process_parser.add_argument('--edf-path', metavar='<edf_path>', required=True,
                                help='Path to the EDF file or folder\n'
                                '')
    process_parser.add_argument('--output-path', metavar='<output_path>', default=None, help='Directory to save plots in')
    process_parser.add_argument('--channel-prefix', metavar='<channel_prefix>', default=None, nargs='+', choices=['e', 't'], help='Channel prefix ("e", "t", or "e t" for both)')
    process_parser.add_argument('--duration', metavar='<duration>', type=int, default=None, help='Duration in seconds')

    args = parser.parse_args()

    return args

if __name__ =='__main__':
    args = parse_cli_args()

    if args.command == 'topo_plot_edfs':
        if os.path.isfile(args.edf_path) and args.edf_path.endswith('.edf'):
            for channel_prefix in args.channel_prefix:
                plot_topo_from_edf(edf_path=args.edf_path,
                                   plt_savepath=args.output_path,
                                   duration=args.duration,
                                   channel_prefix=channel_prefix)
        else:
            files = glob.glob(f"{args.edf_path}/**/*.edf", recursive=True)
            for file in files:
                for channel_prefix in args.channel_prefix:
                    plot_topo_from_edf(edf_path=file,
                                       plt_savepath=args.output_path,
                                       duration=args.duration,
                                       channel_prefix=channel_prefix)
                
    elif args.command == 'combine_topo_plots':
        if args.pxs:
            for px in args.pxs:
                filenames = glob.glob(f'{args.figure_folder}/*{px}*.png')
                save_path = f"{args.output_path}/{px}" if args.output_path else px
                combine_topo_plots(filenames, save_path)
        else: 
            filenames = glob.glob(f'{args.figure_folder}/**.png')
            save_path = f"{args.output_path}/stacked_plot" if args.output_path else 'stacked_plot'
            combine_topo_plots(filenames, save_path)
            
    elif args.command == 'spect_plot_edfs':
        if os.path.isfile(args.edf_path) and not args.edf_path.endswith('.edf'):
            print("Only works with multiple EDFs (Names contining sleep stage)")
        else:
            files = glob.glob(f"{args.edf_path}/**/*.edf", recursive=True)
            for channel_prefix in args.channel_prefix:
                plot_spect_from_edf(edf_paths=files,
                                    plt_savepath=args.output_path,
                                    duration=args.duration,
                                    channel_prefix=channel_prefix)
    else:
        print("Invalid command. Available commands: topo_plot_edfs, combine_topo_plots")
