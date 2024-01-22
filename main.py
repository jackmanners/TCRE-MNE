import glob
import pickle
import pandas as pd


from tqdm import tqdm

# from psga.psga import io

# main_folder = r'R:\CMPH-Undersea Decision Superiority Defence Project\2020-86-20 (LED lighting project)\16. LAB Study Data\PSG Files\Sleep Opportunities\LED*'
# folders = glob.glob(main_folder)
# for folder in tqdm(folders):
#     filename = folder.split('\\')[-1]

#     raw, hypnogram, events = io.read_psg_compumedics(folder=folder)
    
#     hypnogram.to_csv(f'{filename}_hypnogram.csv')
#     events.to_csv(f'{filename}_events.csv')
#     with open(f'{filename}_raw.pickle', 'wb') as f:
#         pickle.dump(raw, f)

columns_to_read = ['EEG=0, TCRE=1', 'absolute_delta', 'absolute_theta', 'absolute_alpha', 
                   'absolute_sigma', 'absolute_beta', 'deltaR', 'thetaR', 
                   'alphaR', 'sigmaR', 'betaR', 'Noisy Epochs Sleep (more than 400)',
                   'SleepStage', 'SleepStageOnset', 'SleepStageDuration', 'Epoch Start Time']

files = glob.glob(r'C:\TCRE\Data\**.xlsx')

# Initialize an empty DataFrame to hold concatenated data
concatenated_df = pd.DataFrame()

for file in tqdm(files):
    # Read the specified columns from the first sheet
    sheet_name = pd.ExcelFile(file).sheet_names[0]
    df = pd.read_excel(file, sheet_name=sheet_name, usecols=columns_to_read)

    # Add a 'source' column to identify data from different files if necessary
    df['electrode'] = sheet_name

    # Concatenate the DataFrame to the main DataFrame
    concatenated_df = pd.concat([concatenated_df, df], ignore_index=True)

# Save the concatenated DataFrame to a CSV file
concatenated_df.to_csv('concatenated_data.csv', index=False)
print("All files have been concatenated into one CSV file.")