##import mne
##
##
##
### Load the standard_1020 montage
##montage_standard_1020 = mne.channels.make_standard_montage('standard_1020')
##
##
##### Choose what channels you want to keep
##### Here, the chosen channels correspond to EasyCapM3 33-channel layout 
####kept_channels = ['Fp1',
####'F7',
####'Fp2',
####'F8',
####'F3',
####'Fz',
####'F4',
####'AFz',
####'FT7',
####'FC3',
####'FCz',
####'FC4',
####'FT8',
####'T8',
####'C4',
####'Cz',
####'C3',
####'T7',
####'TP9',
####'TP7',
####'CP3',
####'CPz',
####'CP4',
####'TP8',
####'TP10',
####'P8',
####'P4',
####'Pz',
####'P3',
####'P7',
####'O1',
####'Oz',
####'O2']
####
####
##### get indices of kept_channels in full 10-20 montage 
####ind = []
####for channel in kept_channels:
####    if channel in mont1020.ch_names:
####        ind.append(mont1020.ch_names.index(channel))
####
####mont1020_new = mont1020.copy()
####
##### Keep only the desired channels
####mont1020_new.ch_names = [mont1020.ch_names[x] for x in ind]
####kept_channel_info = [mont1020.dig[x+3] for x in ind]
####
##### Keep the first three rows as they are the fiducial points information
####mont1020_new.dig = mont1020.dig[0:3] + kept_channel_info
####
##### Plot original montage
####mont1020.plot()
####
##### Plot new montage 
####mont1020_new.plot()
##
##
##
##
### Get the positions of electrodes (X,Y,Z spatial coordinates; MNI space)
##
### Load the 'standard_1020' montage
##montage_standard_1005 = mne.channels.make_standard_montage('standard_1005')
##
### Get the positions of electrodes
##positions = montage_standard_1005.get_positions()
##
### Extract x, y, z coordinates
##x_coordinates = []
##y_coordinates = []
##z_coordinates = []
##
##for pos in positions.values():
##    if isinstance(pos, dict) and 'r' in pos:
##        x_coordinates.append(pos['r'][0])
##        y_coordinates.append(pos['r'][1])
##        z_coordinates.append(pos['r'][2])
##
### Print or use these coordinates as needed
##print("X coordinates:", x_coordinates)
##print("Y coordinates:", y_coordinates)
##print("Z coordinates:", z_coordinates)
##
##
##
##
##



import pandas as pd

# Specify the file path
file_path = '/Users/jane_chesley/Documents/Github/Dataset_dynamic/EEG-analysis/input/10-20toMNI.xlsx'

# Load the Excel file
data = pd.read_excel(file_path)

# Display the loaded data
print(data)
