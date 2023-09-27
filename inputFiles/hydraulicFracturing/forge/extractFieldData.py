import pandas as pd

# Define the CSV file path
csv_file = "cycle4.csv"

# Read the CSV into a DataFrame
df = pd.read_csv(csv_file, header=None, names=["Timestamp", "Tubing pressure", "Casing pressure"])

# Convert the 'Timestamp' column to datetime objects
df['Timestamp'] = pd.to_datetime(df['Timestamp'], format="%m/%d/%Y %H:%M:%S")

# Calculate the time difference in seconds from the first timestamp
df['TimeInSeconds'] = (df['Timestamp'] - df['Timestamp'].iloc[0]).dt.total_seconds()

# Define cut off time
cutoff_timestamp = pd.Timestamp('2019-04-26 11:08:10')

# Create a mask to filter rows before the cutoff timestamp
mask = df['Timestamp'] <= cutoff_timestamp

# Use the mask to extract data until the cutoff timestamp
filtered_df = df[mask]

# Define the path to the HDF5 file
hdf5_file = "cycle4_dfit.h5"

# Create a new DataFrame with the selected columns
selected_columns = filtered_df[['TimeInSeconds', 'Tubing pressure', 'Casing pressure']]

# Export the selected DataFrame to the HDF5 file
selected_columns.to_hdf(hdf5_file, key='data', mode='w')

