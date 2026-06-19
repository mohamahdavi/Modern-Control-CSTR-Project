# Modern Control Project – CSTR Reactor

This repository contains a multi‑phase project for the *Modern Control* course, 
focusing on the modelling, analysis, and control of a non‑isothermal continuous 
stirred‑tank reactor (CSTR).

## Project Phases

| Phase | Status        | Description |
|-------|---------------|-------------|
| 1     | ✅ Completed  | Nonlinear modelling, linearisation, stability analysis, sensitivity analysis, canonical forms, controllability & observability checks, model validation (RMSE). |
| 2     | 🚧 Planned    | Linear controller design (LQR, PID), observer design, and simulation on the nonlinear plant. |
| 3     | 🚧 Planned    | Advanced control techniques (to be defined). |

## Current Content – Phase 1

All files for Phase 1 are in the [`Phase1/`](Phase1/) folder.

### 📄 Report
- `Phase1/cstr_phase1.pdf` – Full English report with all derivations, tables, and figures.
- `Phase1/cstr_phase1.tex` – LaTeX source.

### 📊 Figures
- `Phase1/figures/` – Pole‑zero map, step response, sensitivity plots, validation curves, RMSE analysis.

### 💻 MATLAB Code
- `Phase1/matlab/` – Scripts for nonlinear simulation, linearisation, stability analysis, canonical transformations, and RMSE validation.

### 🔍 Summary of Phase 1
- Extracted the nonlinear dimensionless state‑space model from the reference paper.
- Linearised around the operating point and computed system matrices.
- Confirmed asymptotic stability (poles at $-1.2113 \pm j0.3594$).
- Performed sensitivity analysis on the heat transfer coefficient $U$.
- Derived controller and observer canonical forms.
- Proved controllability and observability via five methods (matrix rank, Jordan form, Gramians, Popov–Belevitch–Hautus).
- Validated the linear model against the nonlinear plant with step and sinusoidal inputs; maximum RMSE $\approx 2\times 10^{-4}$.

## Future Phases
- Phase 2 will design linear controllers (LQR, PID) and observers, then test them on the nonlinear model.
- Phase 3 may explore robust or nonlinear control approaches.

## Course Details
- **Course:** Modern Control (25792)
- **Instructor:** Prof. Hossein Pourshamsaei
- **Date:** June 2026

## Authors
- **Mohammad Reza Mahdavi**
- **Mohammad Mahdi Barzegar**

*Status: In Progress – remaining sections will be added soon.*
