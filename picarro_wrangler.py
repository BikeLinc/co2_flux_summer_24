import os
import pandas as pd
from datetime import datetime
import pytz
from tqdm import tqdm
import tkinter as tk
from tkinter import filedialog, simpledialog

def import_picarro_file(filename, selected_columns=None, resample_freq=None):
    """
    Import a single Picarro .dat file into a Pandas DataFrame.

    Parameters:
    - filename: Path to the .dat file.
    - selected_columns: List of column names to keep (optional)
    - resample_freq: Frequency string for resampling (e.g., '1min', '5min', '1H')

    Returns:
    - df: Pandas DataFrame containing the imported data.
    """
    # Define column names and data types
    column_names = [
        "DATE", "TIME", "FRAC_DAYS_SINCE_JAN1", "FRAC_HRS_SINCE_JAN1", "JULIAN_DAYS",
        "EPOCH_TIME", "ALARM_STATUS", "INST_STATUS", "CavityPressure", "CavityTemp",
        "DasTemp", "EtalonTemp", "WarmBoxTemp", "species", "MPVPosition", "OutletValve",
        "solenoid_valves", "CO", "CO2", "CO2_dry", "CH4", "CH4_dry", "H2O",
        "h2o_reported", "b_h2o_pct", "peak_14", "peak84_raw"
    ]
    
    # Define data types for each column
    column_types = {
        "DATE": "string",
        "TIME": "string",
        "FRAC_DAYS_SINCE_JAN1": float,
        "FRAC_HRS_SINCE_JAN1": float,
        "JULIAN_DAYS": float,
        "EPOCH_TIME": float,
        "ALARM_STATUS": float,
        "INST_STATUS": float,
        "CavityPressure": float,
        "CavityTemp": float,
        "DasTemp": float,
        "EtalonTemp": float,
        "WarmBoxTemp": float,
        "species": float,
        "MPVPosition": float,
        "OutletValve": float,
        "solenoid_valves": float,
        "CO": float,
        "CO2": float,
        "CO2_dry": float,
        "CH4": float,
        "CH4_dry": float,
        "H2O": float,
        "h2o_reported": float,
        "b_h2o_pct": float,
        "peak_14": float,
        "peak84_raw": float
    }

    try:
        # Read the file
        df = pd.read_csv(
            filename,
            sep=r'\s+',  # Fixed: Use raw string to avoid escape sequence warning
            names=column_names,
            dtype=column_types,
            skiprows=1,
            engine='c',
            na_values=["NA", "NaN", "nan", "---"]
        )

        # Combine DATE and TIME into a single datetime column called "T"
        df["T"] = pd.to_datetime(
            df["DATE"] + " " + df["TIME"], format="%Y-%m-%d %H:%M:%S.%f", errors="coerce"
        )

        # Drop invalid rows and unnecessary columns
        df = df.drop(columns=["DATE", "TIME"])
        df = df.dropna(subset=["T"])
        
        # Simply subtract 7 hours to convert from UTC to MST
        df["T"] = df["T"] - pd.Timedelta(hours=7)
        
        # Set T as the index
        df = df.set_index("T")
        
        # Keep only selected columns if specified
        if selected_columns and all(col in df.columns for col in selected_columns):
            df = df[selected_columns]
            
        # Resample data if frequency is provided
        if resample_freq:
            df = df.resample(resample_freq).mean()

        return df

    except pd.errors.ParserError as pe:
        print(f"[PARSER ERROR] Failed to parse {filename}: {pe}")
    except Exception as e:
        print(f"[ERROR] An error occurred while processing {filename}: {e}")

    return pd.DataFrame()  # Return empty DataFrame on failure

def retime_data(df, freq='1min'):
    """
    Resample time series data to a specified frequency.
    
    Parameters:
    - df: DataFrame with datetime index
    - freq: Frequency string for resampling (e.g., '1min', '5min', '1H')
    
    Returns:
    - Resampled DataFrame
    """
    return df.resample(freq).mean()

def import_picarro_subfolders(folder, selected_columns=None, resample_freq=None):
    """
    Concatenate Picarro .dat datasets from subfolders into a single Pandas DataFrame.

    Parameters:
    - folder: Path to the root folder containing subfolders with .dat files.
    - selected_columns: List of column names to keep (optional)
    - resample_freq: Frequency string for resampling (e.g., '1min', '5min', '1H')

    Returns:
    - combined_df: Pandas DataFrame containing all the concatenated data.
    """
    # Find all .dat files in the folder and its subfolders
    dat_files = [
        os.path.join(dirpath, file)
        for dirpath, _, files in os.walk(folder)
        for file in files if file.endswith(".dat")
    ]

    combined_data = []

    print(f"[INFO] Found {len(dat_files)} .dat files. Processing...")
    for file in tqdm(dat_files, desc="Importing files"):
        df = import_picarro_file(file, selected_columns, resample_freq)
        if not df.empty:
            combined_data.append(df)

    # Concatenate all data into a single DataFrame
    if combined_data:
        combined_df = pd.concat(combined_data)
        combined_df = combined_df.sort_index()  # Sort by time
    else:
        combined_df = pd.DataFrame()
        print("[WARNING] No data was loaded. Please check the input files.")

    return combined_df

def select_and_process():
    """
    Open a file dialog to select a file or folder and process accordingly.
    """
    root = tk.Tk()
    root.withdraw()  # Hide the main window
    
    # Ask the user if they want to select a file or folder
    file_or_folder = simpledialog.askstring("Select Type", 
                                            "Enter 'file' to select a .dat file or 'folder' to select a directory:",
                                            initialvalue="file")
    
    path = None
    if file_or_folder and file_or_folder.lower() == 'file':
        path = filedialog.askopenfilename(
            title="Select a .dat file",
            filetypes=[("DAT files", "*.dat"), ("All files", "*.*")]
        )
    elif file_or_folder:
        path = filedialog.askdirectory(
            title="Select a folder containing .dat files"
        )
    
    if not path:  # User canceled
        print("No file or folder selected.")
        return None
    
    # Ask for resampling frequency before processing data
    root = tk.Tk()
    root.withdraw()
    resample_freq = simpledialog.askstring("Input", 
                                         "Enter resampling frequency (e.g., '1min', '5min', '1H'):", 
                                         initialvalue='1min')
    
    if resample_freq:
        print(f"[INFO] Data will be resampled to {resample_freq} frequency")
    
    selected_columns = ["CO", "CO2", "CO2_dry", "CH4", "CH4_dry", "H2O"]
    
    if os.path.isfile(path) and path.endswith('.dat'):
        # Process single file
        print(f"Processing file: {path}")
        df = import_picarro_file(path, selected_columns, resample_freq)
        return df
    
    elif os.path.isdir(path):
        # Process folder
        print(f"Processing folder: {path}")
        return import_picarro_subfolders(path, selected_columns, resample_freq)
    
    else:
        print("Selected path is not a .dat file or folder.")
        return None

if __name__ == "__main__":
    # Use the file dialog to select file or folder
    data = select_and_process()
    
    if data is not None and not data.empty:
        # Save to CSV
        root = tk.Tk()
        root.withdraw()
        save_path = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
            title="Save processed data as"
        )
        
        if save_path:
            # No timezone conversion needed - datetimes are already naive
            data.to_csv(save_path)
            print(f"[INFO] Data saved to {save_path}")
        
        # Display summary
        print("\n[INFO] Data Summary:")
        print(data.describe())
        print("\n[INFO] First few rows of the data:")
        print(data.head())
    else:
        print("[INFO] No data was loaded.")



