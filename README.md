# Modern Control Project – CSTR Reactor

Nonlinear modeling, linearization, and control design for a non‑isothermal continuous stirred‑tank reactor (CSTR).

**Electrical Engineering Department, Sharif University of Technology**

## Phases Overview

| Phase | Status | Summary |
|-------|--------|---------|
| 1     | ✅      | Nonlinear modeling, stability & sensitivity, canonical forms, controllability/observability (5 methods), linear model validation (RMSE). |
| 2     | ✅      | State‑feedback control, domain of attraction, invariant realizations, servo + integral action, Luenberger & PI observers, saturation/anti‑windup, robustness to +30% inlet flow. |

---

### 📁 Phase 1 – [`Phase1/`](Phase1/)

- Report (PDF + LaTeX), figures, MATLAB scripts.
- **Highlights:** Derived nonlinear model; linearized at equilibrium (poles –1.1182, –1.8354); sensitivity analysis on \(U\); verified controllability/observability via rank, Jordan, Gramians, PBH; validated against nonlinear plant (RMSE ≈ 2×10⁻⁴).

### 📁 Phase 2 – [`Phase2/`](Phase2/)

- Report (PDF + LaTeX), figures, MATLAB scripts.
- **Highlights:** Designed state feedback (poles –3, –4.5); mapped domain of attraction; showed transfer‑function invariance across 6 realizations; added integral action (servo); built Luenberger and PI observers; handled actuator saturation with anti‑windup; confirmed robustness under +30% \(F_{os}\).

---

## Course Details

- **Course:** Modern Control (25792)
- **Instructor:** Prof. Hossein Pourshamsaei
- **Department:** Electrical Engineering, Sharif University of Technology
- **Date:** August 2026

## Authors

- Mohammad Reza Mahdavi
- Mohammad Mehdi Barzegar

---

*Project complete – all code, figures, and reports are included.*
