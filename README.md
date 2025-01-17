# CO₂ Flux Project - Summer 2024

This repository contains firmware, data processing scripts, and design files for the CO₂ Flux Chamber project (2022-2024) conducted by the Casey Air Quality Lab at Fort Lewis College. The system enables spatial monitoring of soil CO₂ fluxes in wildfire-affected areas using a dynamic flux chamber system with dual ELT S300-3V CO₂ sensors calibrated against a LICOR LI-7810 and PICARRO CRDS.

**Note: This repository is inactive as of 01-10-2025 and is no longer maintained.**

---

# Project Overview

## Purpose
The CO₂ Flux Project aims to measure soil CO₂ emissions in wildfire-affected areas using dynamic flux chambers and cost-effective $70 NDIR sensors as alternatives to $5000 commercial IR gas analyzers. The project includes prototype development, lab validation, and calibration studies.

## Key Contributions
- **Modeling:** Mass balance models and literature reviews informed the design of optimal flux chambers and theoretical models to interpret data.
- **Data Analysis:** Processed datasets to calculate CO₂ fluxes and assess sensor performance.
- **Fabrication:** Designed and built hardware, including PCBs and data acquisition systems.

## Progress Highlights
- **Lab Validation:** Successfully validated the flux chamber system in controlled conditions.
- **Prototype Development:** Built and tested the initial device.
- **Calibration Studies:** Conducted calibrations for ELT S300-3V sensors using LICOR LI-7810 and PICARRO CRDS as benchmarks.
- **Manuscript:** Completed the first draft of the manuscript.

## Challenges
- **microSD Card Data Logging Issues:** Overloaded SPI interface resolved by implementing caching methods and creating multiple data files.
- **ELT S300 Baseline Shift:** No significant shifts detected; older calibrations exhibit inconsistent predictions over time.
- **Unlabeled Datasets:** Inadequate labeling and documentation practices caused data mismanagement at deployment sites.

---

# Repository Structure

| Folder          | Description                                                                 |
|-----------------|-----------------------------------------------------------------------------|
| **DATA**        | Contains datasets and metadata. File extensions (.picarro, .daq, .licor) indicate the instrument used. Naming convention: `MM_DD_YY_A_B.<EXT>`. Flux datasets primarily use .daq or .csv files. |
| **ELT CALIBRATION** | Calibration scripts for ELT S300-3V sensors. Key script: `ELT CALIB/DYNAMIC/calibration_v2.m`. |
| **LICOR CALIBRATION** | Calibration scripts for LI-7810 IR gas analyzer. |
| **MFC CALIBRATION** | Calibration scripts for ALICAT MFCs. |
| **MODELING**     | Scripts for theoretical modeling of dynamic flux chambers. Key scripts: `MODELING/measurement_range.mlx`, `MODELING/flow_rate_sweet_spot_new.m`. |
| **PARAM_SPACE**  | Scripts for evaluating calibration parameter sensitivity. Key script: `PARAM_SPACE/parameter_space.m`. |
| **PCB**          | Hardware design files and firmware. Key files: `PCB/firmware/firmware.ino`, `PCB/pcb_rev4/` for PCB design files. |
| **FLUXES**       | Flux calculation and data processing scripts. Note: Structural issues persist. Key script: `ARCHIVED_CODE/ARCHIVED_CODE/co2_flux_summer_24/FLUXES/flux_measurement.m`. |
| **UTILS**        | Shared utility functions and configuration files. Key script: `CONFIG.m`. |

---

# Contact

For historical questions about this project, please contact the Casey Air Quality Lab at Fort Lewis College.

