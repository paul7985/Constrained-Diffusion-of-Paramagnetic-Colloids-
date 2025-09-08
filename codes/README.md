# Colloidal Particles Motion Analysis

This repository contains MATLAB scripts developed for the analysis of the motion of **paramagnetic colloids** subjected to a vertical magnetic field.  
The project was conducted during a research internship at **LCP-A2MC**, and is based on the processing of particle trajectories extracted by microscopy.

The main objective is to compute and interpret the **Mean Squared Displacement (MSD)** under various conditions, taking into account possible system drift.  
Several approaches have been explored: global vs. local treatments, with or without correction, individual vs. collective analyses — in order to highlight a potential **caging effect** (particles confined by their neighbors).

Each script is described in context, with a clear explanation of its role, inputs/outputs, and use cases.

---

## 📂 Repository Organization

Scripts are organized according to their main function:

1. **Standard MSD calculations**  
2. **Drift correction** (global or local, including moving average methods)  
3. **Restricted analyses** (quiet zones, individual particles)  
4. **Parametric studies** (dependence on correction radius, time, etc.)  
5. **Comparisons across experiments**  
6. **Visualization and utilities**  

---

## 1️⃣ Standard MSD Calculations

- **`calcul_MSD.m`**  
  Loads multiple CSV files containing particle trajectories. Applies time/position conversion and computes the average MSD. Useful for a quick first look at dynamics without drift correction.

- **`calcul_msd_r_david_OK.m`**  
  A more complete version of the above. Includes pixel → μm conversion, file selection dialog, and computes MSD for each file before averaging. Reliable for producing final curves with physical units.

- **`msd_individuel.m`**  
  Computes MSD for each particle individually (no averaging). Ideal to visualize disparities between particles, detect outliers, and assess dispersion in dynamics.

- **`MSD_1_intervalle_STD.m`**  
  Computes MSD at a single lag time across all files. Useful for quick comparisons between experiments (e.g. at τ = 10s).

---

## 2️⃣ Drift Correction (Local or Global Methods)

- **`MSD_moyenne_mobile_traj.m`**, **`test_MSD_moyenne_mobile.m`**, **`test_MSD_traj_moyenne_mobile.m`**  
  Apply local drift correction using a moving average: the mean velocity of neighboring particles within a given radius is subtracted from each trajectory. Corrects convection or global flow while preserving local dynamics. Essential to reveal diffusive motion.

- **`MSD_rayon_MA.m`**  
  Explores the influence of the neighborhood radius on MSD results. Provides comparative plots to select the best trade-off between drift removal and signal preservation.

---

## 3️⃣ Restricted Analyses (Quiet Zones, Individual Particles)

- **`msd_zone_restreinte.m`**, **`msd_zone_restreinte_et_traj.m`**  
  Filter trajectories to keep only those within a defined spatial zone (X,Y). Allows isolation of regions with weak drift or more homogeneous behavior. The `_et_traj` version also plots the selected trajectories.

- **`msd_solo_comprehension.m`**  
  Studies MSD for a manually chosen single particle. Useful for error diagnosis or testing whether unusual behavior is generalizable.

---

## 4️⃣ Parametric Studies

- **`msd_fonction_r.m`**, **`msd_fonction_r_2.m`**, **`msd_fonction_r_max_test.m`**  
  Test different neighborhood radii in local drift correction and study their effect on MSD. Produce plots like “MSD vs radius” (e.g. Figure 11 in the report).

- **`msd_solo_fonction_r.m`**, **`msd_solo_fonction_r_2.m`**  
  Same objective as above, but applied to a single particle. Evaluates whether radius dependence is uniform or trajectory-specific.

---

## 5️⃣ Comparisons Across Experiments

- **`code_comparaison_msd.m`**, **`code_comparaison_msd_sans_drift.m`**,  
  **`code_comparaison_msd_1intervalle_sans_code_comparaison_msd_sans_drift_v2_test.m`**  
  Compare multiple experiments (or regions of a single experiment), with or without drift correction. Some average over multiple files, others handle a single file. Useful for testing reproducibility and assessing parameter effects (magnetic field, concentration, particle size, etc.).

---

## 6️⃣ Visualization and Utilities

- **`Plot_superposition_image.m`**  
  Overlays particle trajectories onto a background image (e.g. from ImageJ). Helps validate that detected trajectories match raw images.

- **`plot_trajectoires.m`**  
  Displays all detected trajectories in one plot. Useful for spotting strong drift regions or artifacts.

- **`Track_Particles.m`**, **`track.m`**  
  Particle tracking functions. Generate trajectories by linking positions across frames, using a maximum displacement parameter (`maxdisp`).

---

✨ This README serves as a **reference guide** to navigate the scripts and understand their role in the colloid motion analysis project.
