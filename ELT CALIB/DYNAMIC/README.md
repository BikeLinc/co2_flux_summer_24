# CO2 Flux Sensor Calibration Tool

This tool generates calibration models for all sensors listed in the overlapping_periods.csv file. The calibration procedure has been designed to produce both linear regression and neural network models that can be later used to convert raw sensor readings into accurate CO2 measurements.

## Getting Started

1. Clone the repository by excecuting the following command

`git clone -b gui-version https://github.com/BikeLinc/co2_flux_summer_24.git`

2. Ensure your data files are properly organized in the `../../DATA/CALIB/` directory
3. Make sure `overlapping_periods.csv` file is present and properly formatted
4. Launch MATLAB and navigate to `calibration_v2.m` script.
5. Run the `calibration_v2.m` script, select `Change Folder` if prompted.

## CSV File Format Requirements

The `overlapping_periods.csv` file must contain the following columns:

| Column Name | Description |
|-------------|-------------|
| Type | Data type identifier ('AMB' or 'CL') |
| DAQFile | Path to the DAQ data file (relative to DATA/CALIB/) |
| LicorFile | Path to the Licor reference data file (relative to DATA/CALIB/) |
| PICARROFile | Path to the PICARRO reference data file (relative to DATA/CALIB/) |
| OverlapStart | Start time of the overlap period (format: MM/dd/yyyy HH:mm:ss) |
| OverlapEnd | End time of the overlap period (format: MM/dd/yyyy HH:mm:ss) |
| SensorA | ID of the first sensor in the DAQ file (e.g., 'A1', 'B3', or 'X' if not used) |
| SensorB | ID of the second sensor in the DAQ file (e.g., 'A1', 'B3', or 'X' if not used) |

Example CSV format:
```
Type,DAQFile,LicorFile,PICARROFile,OverlapStart,OverlapEnd,SensorA,SensorB
AMB,daq_file1.csv,licor_file1.csv,,05/10/2024 10:30:00,05/10/2024 15:45:00,A1,B2
CL,,picarro_file1.csv,daq_file2.csv,05/12/2024 08:00:00,05/12/2024 11:30:00,A2,B3
```

## Parameter Configuration UI

The script uses a two-panel approach for configuring calibration parameters:

### Panel 1: Data Processing Parameters

This panel collects parameters related to data preprocessing and filtering:

| Parameter | Description |
|-----------|-------------|
| Smooth duration | Window size (in minutes) for data smoothing |
| Retime duration | Interval (in minutes) for resampling data to regular time intervals |
| Outlier percentiles | Array defining the percentile range for outlier removal, e.g., `[2, 98]` |
| Enable outlier removal | Toggle outlier removal (true/false) |
| Reference minimum value | Minimum acceptable reference CO2 value (ppm) |
| Reference maximum value | Maximum acceptable reference CO2 value (ppm) |
| Data type to include | Select data subset to process ('AMB', 'CL', or 'Both') |

### Panel 2: Model Training Parameters

This panel collects parameters related to model training and evaluation:

| Parameter | Description |
|-----------|-------------|
| Number of bins | Number of bins to divide the data into for stratified sampling |
| Target bin count | Minimum number of samples per bin for data replication |
| Training fraction | Fraction of data used for model training (e.g., 0.5 = 50%) |
| Validation fraction | Fraction of data used for validation (e.g., 0.1 = 10%) |
| Evaluation fraction | Fraction of data used for final evaluation (e.g., 0.4 = 40%) |
| Enable data replication | Toggle data replication for balancing bins (true/false) |
| Use random split | Toggle between random or time-ordered data splitting (true/false) |
| Neuron layers | Array defining neural network hidden layer sizes, e.g., `[16, 16]` |
| Maximum epochs | Maximum training iterations for neural networks |
| Predictors | Comma-separated list of predictor variables (e.g., 'X_C,X_T,X_H') |
| Reference ranges | Comma-separated CO2 concentration ranges (format: '390-450,450-1200,390-1200') |
| Range labels | Comma-separated labels for each range (e.g., '390_450,450_1200,390_1200') |

**Note**: The target variable is fixed to 'Y_C' (reference CO2 concentration) and doesn't require user input.

## Output Files

The calibration tool produces the following outputs:

1. Trained models saved in `models/` directory with the naming pattern:
   - Linear models: `<sensorID>-linear-<predictorSet>-<rangeLabel>-<date>.mat`
   - Neural Network models: `<sensorID>-net-<predictorSet>-<rangeLabel>-<date>.mat`

2. Figures saved in `figs_training/` directory, showing:
   - Data partitions 
   - Raw data overlaps
   - Residual plots
   - Time series comparisons between models and reference data

3. Model performance metrics saved in `model_comparison_all_sensors.csv`

## Model Evaluation

The tool includes a built-in model evaluation function that allows you to:
1. Select previously saved models
2. Apply them to evaluation datasets
3. Compare their performance with visual plots and metrics
