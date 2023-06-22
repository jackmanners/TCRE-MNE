import glob
import pickle

from tqdm import tqdm

from psga.psga import io

main_folder = r'R:\CMPH-Undersea Decision Superiority Defence Project\2020-86-20 (LED lighting project)\16. LAB Study Data\PSG Files\Sleep Opportunities\LED*'
folders = glob.glob(main_folder)
for folder in tqdm(folders):
    filename = folder.split('\\')[-1]

    raw, hypnogram, events = io.read_psg_compumedics(folder=folder)
    
    hypnogram.to_csv(f'{filename}_hypnogram.csv')
    events.to_csv(f'{filename}_events.csv')
    with open(f'{filename}_raw.pickle', 'wb') as f:
        pickle.dump(raw, f)