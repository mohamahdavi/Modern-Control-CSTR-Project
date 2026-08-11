# Modern Control Project – CSTR Reactor

This repository contains a multi‑phase project for the *Modern Control* course, 
focusing on the modelling, analysis, and control of a non‑isothermal continuous 
stirred‑tank reactor (CSTR).

## Project Phases

| Phase | Status        | Description |
|-------|---------------|-------------|
| 1     | ✅ Completed  | Nonlinear modelling, linearisation, stability analysis, sensitivity analysis, canonical forms, controllability & observability checks, model validation (RMSE). |
| 2     | ✅ Completed  | State‑feedback control, domain of attraction mapping, invariant transfer‑function realisations, servo control with integral action, Luenberger and PI observer design, observer‑based control, actuator saturation and anti‑windup, robustness to 30 % inlet flow increase. |

## Repository Content

### 📁 Phase 1

All files for Phase 1 are in the [`Phase1/`](Phase1/) folder.

- **Report**: `Phase1/cstr_phase1.pdf` (English) and LaTeX source.
- **Figures**: Pole‑zero map, step response, sensitivity plots, validation curves, RMSE analysis.
- **MATLAB Code**: Scripts for nonlinear simulation, linearisation, stability, canonical forms, controllability/observability, and RMSE validation.

**Highlights of Phase 1:**
- Derived nonlinear dimensionless state‑space model from the reference paper.
- Linearised around the operating point (poles at –1.1182, –1.8354, asymptotically stable).
- Sensitivity analysis on heat transfer coefficient \(U\).
- Derived controller and observer canonical forms; verified invariance of transfer function.
- Proved controllability and observability via five methods (rank, Jordan, Gramians, PBH).
- Validated linear model against nonlinear plant with step/sinusoidal inputs; maximum RMSE ≈ 2×10⁻⁴.

---

### 📁 Phase 2

All files for Phase 2 are in the [`Phase2/`](Phase2/) folder.

- **Report**: `Phase2/cstr_phase2.pdf` (English) and LaTeX source.
- **Figures**: Time responses, domain of attraction, realisation comparisons, servo control performance, integrator windup, observer speed vs. noise, observer‑based control, saturation effects, anti‑windup, robustness tests.
- **MATLAB Code**: Scripts for controller design, observer design, simulation with nonlinear plant, saturation handling, and robustness analysis.

**Highlights of Phase 2:**
- Designed state‑feedback gain (poles at –3, –4.5) and mapped its limited **domain of attraction** on the nonlinear plant.
- Constructed six different realisations (controller, observer, Jordan forms) and confirmed that the **transfer function is invariant** under similarity transformations.
- Implemented a **servo controller with integral action** to eliminate steady‑state error; observed **integrator windup** under large disturbances.
- Designed both a **Luenberger observer** (poles at –9, –12) and a **PI observer** (augmented poles at –1, –9, –12) – the PI observer removes steady‑state estimation error.
- Studied the **speed‑vs‑noise trade‑off** by tuning observer poles; found the “mid” setting (poles –1, –9, –12) offers the best compromise.
- Combined observer with controller and compared against ideal state feedback – transient performance is slightly degraded but acceptable.
- Addressed **actuator saturation** (coolant flow bounded between 0 and 2) and implemented a simple **anti‑windup** mechanism (freezing integrator) to maintain stability for gains up to factor 4.
- Performed a **robustness test** with a +30 % increase in inlet flow rate \(F_{os}\); the system remained stable and tracked the desired setpoint even with saturation.


## Course Details

- **Course:** Modern Control (25792)
- **Instructor:** Prof. Hossein Pourshamsaei
- **Date:** August 2026

## Authors

- **Mohammad Reza Mahdavi**
- **Mohammad Mehdi Barzegar**

---

*The project is complete. All code, figures, and reports are provided in their respective phase folders.*
