# Pipeline

1. Import Data
    - `importdaqfile`
    - `importlicorfile`
2. Correct Timestamps
    - `detectoffset`
    - `applyoffset`
3. Smooth Dataset
    1. Windowing Function
    `movmean`
    2. High Pass
    `hipass`
    3. Retiming
    `retime`
4. Apply Specific Correction

    - `load correction`
    - `apply correction`
5. Merge Datasets w/ Reference if Necessary

    - `merge`
6. Determine Steady State Periods

    - `movmean < threshold`

7. Calculate Flux

    - Use gradient method on non-steady state values.
    - Use steady-state model on assumed steady state periods.

8. Format Data

    - Organize and document all actions performed on the datasets.
    - Organize all of the data to one .csv file showing important numbers.

9. Display Data

    - Display important information

10. Generate Plots

    - Generate plots
    - Display Important Plots

11. Save Data

## Information
- What ELTs map to which DAQ column.
- Calibration or Colocation Dataset
- The flowrate of the 