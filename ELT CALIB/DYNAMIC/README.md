# CO2 Sensor Calibration Tool - User Guide

## Overview

This tool generates calibration models for CO2 sensors using both linear models and neural networks. It processes data from specified DAQ and reference files to create calibration equations based on CO2 concentration, temperature, and humidity values.

## Requirements

- MATLAB (R2019b or newer recommended)
- Neural Network Toolbox
- Statistics and Machine Learning Toolbox
- Input data files (see Data Format section)

## Getting Started

1. Make sure your input data is properly formatted (see Data Format section)
2. Place your data files in the proper directories
3. Run `calibration_v2.m` script
4. Configure parameters through the displayed dialog
5. Review the generated models and output files

## Data Format

### Required Input Files

1. **CSV Mapping File**: A file named `overlapping_periods.csv` containing columns:
   - `DAQFile`: Path to DAQ sensor data file
   - `LicorFile`: Path to LICOR reference data file (if applicable)
   - `PICARROFile`: Path to PICARRO reference data file (if not using LICOR)
   - `OverlapStart`: Start time of calibration period (MM/dd/yyyy HH:mm:ss)
   - `OverlapEnd`: End time of calibration period (MM/dd/yyyy HH:mm:ss)
   - `SensorA`: ID of the first sensor in the DAQ file
   - `SensorB`: ID of the second sensor in the DAQ file
   - `Type`: Type of data ("AMB" or "CL")

### Data Organization

#### CSV Structure Example
Here's an example of how your `overlapping_periods.csv` should be structured:

```
DAQFile,LicorFile,PICARROFile,OverlapStart,OverlapEnd,SensorA,SensorB,Type
sensor_data_20230601.daq,licor_data_20230601.licor,,06/01/2023 10:00:00,06/01/2023 14:30:00,K30_1,K30_2,AMB
sensor_data_20230602.daq,,picarro_data_20230602.picarro,06/02/2023 09:15:00,06/02/2023 16:45:00,SBA_5,X,CL
```

#### File Naming Conventions
- **DAQ files**: Should use a `.daq` extension (contains sensor measurements)
- **LICOR files**: Should use a `.licor` extension (reference CO2 data)
- **PICARRO files**: Should use a `.picarro` extension (alternative reference data)

You can use any naming convention for the actual filenames, but the extensions help the script identify the file types.

#### Using "X" in the CSV
- In the `SensorA` or `SensorB` columns, an "X" means that channel doesn't contain a valid sensor.
  - Example: If `SensorB` contains "X", only data from `SensorA` will be processed.
- In the `LicorFile` or `PICARROFile` columns, leave one empty if you're using the other.
  - The script expects either a LICOR file or a PICARRO file, but not both.

#### File Requirements
1. **DAQ files**: Must contain columns:
   - `T` (timestamp)
   - `CA`, `TA`, `HA` (Channel A: CO2, temperature, humidity)
   - `CB`, `TB`, `HB` (Channel B: CO2, temperature, humidity)
   
2. **LICOR files**: Must contain columns:
   - `T` (timestamp)
   - `C` (reference CO2 concentration)
   
3. **PICARRO files**: Must contain columns:
   - `T` (timestamp)
   - Either `CO2` or `CO2_sync` (reference CO2 concentration)

### Directory Structure

The script expects the following directory structure:
```
co2_flux_summer_24_old/
├── DATA/
│   └── CALIB/
│       ├── [Your DAQ data files]
│       └── [Your LICOR/PICARRO reference files]
├── ELT CALIB/
│   └── DYNAMIC/
│       ├── calibration_v2.m
│       ├── getCalibrationParameters.m
│       └── overlapping_periods.csv
└── UTILS/
    └── [Utility functions]
```

## Parameter Configuration

When you run the script, a parameter dialog will appear. Here's what each parameter means:

### Data Selection
- **Data Type**: Choose "AMB" for ambient, "CL" for closed-loop (enhanced), or "Both" for all data sets

### Preprocessing Parameters
- **Smoothing Duration**: Window size for smoothing data (in minutes)
- **Retiming Duration**: Interval for resampling data (in minutes)
- **Outlier Percentiles**: Range for identifying outliers (e.g., [2, 98])
- **Outlier Removal**: Enable/disable removal of outlier data points
- **Reference Min/Max**: Valid range for reference CO2 values (ppm)

### Data Partitioning
- **Number of Bins**: Number of bins for partitioning data
- **Target Bin Count**: Target number of points per bin
- **Train/Validation/Eval Fractions**: Data split proportions (must sum to 1.0)
- **Replicate Data**: Enable/disable data replication for underrepresented bins
- **Random Split**: Enable/disable random splitting (vs. sequential)

### Model Parameters
- **Predictors**: Variables to use as predictors (e.g., ["X_C", "X_T", "X_H"])
- **Target Variable**: Variable to predict (usually "Y_C")
- **Reference Ranges**: CO2 ranges for separate calibrations
- **Range Labels**: Labels for each calibration range
- **Neural Network Layers**: Number of neurons in each hidden layer
- **Max Epochs**: Maximum training epochs for neural networks
- **Max Training Tries**: Number of attempts to train neural networks
- **Target RMSE**: Target error for stopping neural network training

## Execution Flow

1. Parameter configuration dialog appears
2. Script loads input files based on settings
3. Data preprocessing and alignment occurs
4. Sensor data is collected and cleaned
5. Data is partitioned into training/validation/evaluation sets
6. Linear and neural network models are trained
7. Model performance is evaluated and compared
8. Models and figures are saved to output directories

## Output Files

### Models
- Saved in the "models" directory with filename format:
  - `[SensorID]-linear-[Predictors]-[Date].mat`: Linear models
  - `[SensorID]-net-[Predictors]-[Date].mat`: Neural network models

### Figures
- Saved in the "figs_training" directory:
  - Raw data plots
  - Data partitioning visualizations
  - Residual plots
  - Time series predictions
  - Cross-correlation corrections

### Reports
- `model_comparison_all_sensors.csv`: Table comparing all models' performance

## Evaluating Saved Models

The script includes functionality to evaluate previously saved models on new data:
1. Run the script through the model training
2. When prompted, select model files to evaluate
3. Review the evaluation plots and performance metrics

## Troubleshooting

### Common Issues

1. **Missing Data Files**: 
   - Ensure files exist in the correct directories
   - Check CSV has correct file paths

2. **Parameter Errors**:
   - Ensure fraction values sum to 1.0
   - Verify all numeric parameters are valid numbers
   - Check array formats match expected formats

3. **Not Enough Data**:
   - Increase date range in overlapping_periods.csv
   - Reduce bin count or target bin count parameters

4. **Poor Model Performance**:
   - Try different combinations of predictors
   - Adjust neural network architecture
   - Check reference data quality

### Additional Common Issues

5. **CSV Format Issues**:
   - Ensure date formats match MM/dd/yyyy HH:mm:ss exactly
   - Verify that file extensions are correctly specified (.daq, .licor, .picarro)
   - Check that every DAQ file has at least one reference file (either LICOR or PICARRO)
   - Make sure SensorA and SensorB values match the expected sensor IDs or "X"

6. **Channel Assignment**:
   - Channel A data is mapped to SensorA in the CSV
   - Channel B data is mapped to SensorB in the CSV
   - Incorrect channel assignment will result in data mismatch