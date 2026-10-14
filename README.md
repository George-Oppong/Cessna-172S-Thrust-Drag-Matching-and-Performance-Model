# Cessna 172S Thrust–Drag Matching and Performance Model

A Python performance model of the Cessna 172S that computes drag and available thrust across speed and altitude, solves for the minimum-drag speed three different ways, and checks the model against the aircraft's Pilot's Operating Handbook (POH).

**Author:** [George Oppong Wiafe] · **Contact:** [georgekyekyeku@gmail.com] · **Status:** complete (undergraduate project)

![Thrust-drag matching](thrust_drag_matching.png)

## What this project does

1. **Models the aircraft.** ISA atmosphere, parabolic drag polar, and a fixed-pitch propeller power model, with every parameter traced to a source (see [Sources](#sources)).
2. **Finds the minimum-drag speed three ways:** an analytical closed form, gradient descent, and differential evolution.
3. **Validates against the POH:** cruise power (37 table rows), maximum speed, rate of climb at seven altitudes, service ceiling, and best glide.
4. **Tests sensitivity** to the parameters that are not well known.

## Method

**Atmosphere.** US Standard Atmosphere 1976, troposphere (density from ISA temperature and pressure).

**Drag.** `D = ½ρV²S (CD0 + k CL²)`, with `CL = W / (½ρV²S)` and `k = 1/(π·AR·e)`.

**Minimum drag (analytical).** Setting `dD/dV = 0` gives `CD0 = k·CL²`, so

- `CL_md = √(CD0/k)`
- `V_md = √(2W / (ρ·S·CL_md))`
- `D_min = 2W·√(CD0·k)`, which does not depend on altitude. Only the speed at which it occurs changes with density.

**Propulsion.** `P_avail = η(J) · BHP_WOT(h) · RPM/2700`, thrust = `P_avail / V`.

- Wide-open-throttle power falls linearly from 180 hp at sea level to 80 hp at 20,000 ft (Lycoming chart, as reproduced by Harada).
- Propeller efficiency η(J) is Harada's Eq. 84, a regression to X-Plane simulator output (no measured curve for this propeller was found).
- Full-throttle RPM is assumed to rise linearly from 2350 (static) to 2700 at the POH maximum speed.

**Aircraft data (POH):** weight 2550 lb, wing area 174 ft², span 36 ft 1 in (AR = 7.48), 180 hp at 2700 RPM, 76 in fixed-pitch propeller, stall speed 53 KCAS (gives CLmax = 1.54).

**Drag parameters:** CD0 = 0.0344 and e = 0.747 (Min et al., 2015).

## Results

### 1. Minimum-drag speed: three methods agree

At 2000 m: V_md = 42.36 m/s, D_min = 1004 N.

| Method | Error in V_md | Function evaluations |
|---|---|---|
| Analytical | exact | 0 |
| Gradient descent (4 starting points) | < 3×10⁻⁷ m/s | 37–51 |
| Differential evolution | 5×10⁻⁵ m/s | 225 |

The problem is smooth with a single minimum, so the closed form is exact and needs no search. The numerical methods serve as a check on the algebra. A genetic or evolutionary method would only be justified for a non-smooth or multi-objective version of the problem.

### 2. Validation against the POH

| Check | Model | POH | Comment |
|---|---|---|---|
| Cruise %BHP (37 rows) | RMS error 1.8 points (bias −1.2) | Fig 5-8 | No parameter tuned to these rows. All rows have J = 0.66–0.74, so this tests essentially one propeller operating point. |
| Max speed, sea level | 131 kt | 126 KTAS | +4%. The RPM assumption uses the POH V_max, so this check is only semi-independent. |
| Service ceiling | 16,400 ft | 14,000 ft | Overpredicted (see below). |
| Best-glide speed | V_md = 74.6 kt | 68 KIAS | Model L/D at 68 kt is only 2% below its maximum (11.1 vs 11.3); the optimum is flat. |

**Rate of climb (ISA):**

| Altitude (ft) | POH (fpm) | Model (fpm) | Excess-power error |
|---|---|---|---|
| 0 | 729 | 830 | +5.8 kW |
| 4,000 | 597 | 669 | +4.2 kW |
| 8,000 | 408 | 499 | +5.3 kW |
| 12,000 | 221 | 316 | +5.5 kW |

![Validation](validation.png)

### 3. Why climb is overpredicted

The excess-power error is a nearly constant +4 to +6 kW (about 7 hp, 4% of rated power) at every altitude. The percentage error grows with altitude only because the climb rate itself shrinks.

It is not a drag-polar error. Matching sea-level climb by raising CD0 alone would require CD0 ≈ 0.045, which pushes cruise power error to 12.6 points and drops V_max to 118 kt. The gap most likely lies on the propulsion side, where the propeller efficiency at low advance ratio is the least certain input, but this was not isolated further.

### 4. Sensitivity

| Case | Cruise RMS (pts) | V_max (kt) | SL climb (fpm) | Ceiling (ft) |
|---|---|---|---|---|
| Baseline (Min et al.) | 1.8 | 131 | 830 | 16,419 |
| CD0/e from Harada build-up (0.0279 / 0.699) | 9.2 | 141 | 868 | 17,171 |
| CD0/e from Harada simulator (0.0384 / 0.673) | 5.5 | 125 | 754 | 14,645 |
| η(J) × 0.75 (peak ≈ 0.65) | 18.4 | 115 | 456 | 10,581 |
| Static RPM 2300 / 2400 | 1.8 | 131 | 824 / 835 | 16,373 / 16,464 |

Propeller efficiency is the dominant uncertainty. No single drag parameter set matches cruise, top speed and climb together, which is the same CD0–η coupling the cruise data cannot resolve on its own.

## Limitations

- Propeller efficiency η(J) comes from a simulator regression, not measured data for this propeller.
- Constant-torque and linear-RPM assumptions for full-throttle power (A1, A2 in the script header).
- Climb speeds are treated as equivalent airspeed (no instrument or compressibility correction).
- The drag parameters are published estimates for the type, not measurements on a specific aircraft. Published values for CD0 range from about 0.021 to 0.046 depending on method.
- The Gagg–Farrar power-lapse formula used as a cross-check was written from recall and has not been verified against its source.

## Run it

```bash
pip install numpy scipy matplotlib
python thrust_drag_matching.py
```

This prints the three tables above and writes `thrust_drag_matching.png` and `validation.png` in the folder you run it from.

## Sources

- Cessna Aircraft Company, *Cessna 172S Pilot's Operating Handbook*, Rev 5, 2004: geometry, engine, stall speed (Fig 5-3), climb (Fig 5-6), cruise (Fig 5-8).
- S. Min, E. Harrison, H. Jimenez, D. Mavris, "Development of Aerodynamic Modeling and Calibration Methods for General Aviation Aircraft Performance Analysis – A Survey and Comparison of Models," AIAA 2015-2853: CD0 = 0.0344 and e = 0.7469 for the C172S.
- Harada, *Validation of Classical Aircraft Performance Methods Using Flight Test Data from a Cessna 172S*, M.S. thesis, The Ohio State University [year]: build-up and simulator drag parameters, Lycoming WOT power line (Fig 15), propeller efficiency (Eq. 84).
- *U.S. Standard Atmosphere, 1976*, NASA-TM-X-74335.
- S. Gudmundsson, *General Aviation Aircraft Design: Applied Methods and Procedures*, 2013: Gagg–Farrar power lapse (cross-check only).
