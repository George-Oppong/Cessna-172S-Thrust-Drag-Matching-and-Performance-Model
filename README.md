# Cessna-172S-Thrust-Drag-Matching-and-Performance-Model
A Python performance model of the Cessna 172S that computes drag and available thrust across speed and altitude, solves for the minimum-drag speed three different ways, and checks the model against the aircraft's Pilot's Operating Handbook (POH).
Author: George Oppong Wiafe· 
Contact:georgekyekyeku@gmail.com 
Status: Undergraduate project

**What this project does**

Models the aircraft. ISA atmosphere, parabolic drag polar, and a fixed-pitch propeller power model, with every parameter traced to a source (see Sources).
Finds the minimum-drag speed three ways: an analytical closed form, gradient descent, and differential evolution.
Validates against the POH: cruise power (37 table rows), maximum speed, rate of climb at seven altitudes, service ceiling, and best glide.
Tests sensitivity to the parameters that are not well known.

Method
Atmosphere. US Standard Atmosphere 1976, troposphere (density from ISA temperature and pressure).

Drag. D = ½ρV²S (CD0 + k CL²), with CL = W / (½ρV²S) and k = 1/(π·AR·e).

Minimum drag (analytical). Setting dD/dV = 0 gives CD0 = k·CL², so

CL_md = √(CD0/k)
V_md = √(2W / (ρ·S·CL_md))
D_min = 2W·√(CD0·k), which does not depend on altitude. Only the speed at which it occurs changes with density.

Propulsion. P_avail = η(J) · BHP_WOT(h) · RPM/2700, thrust = P_avail / V.

Wide-open-throttle power falls linearly from 180 hp at sea level to 80 hp at 20,000 ft (Lycoming chart, as reproduced by Harada).
Propeller efficiency η(J) is Harada's Eq. 84, a regression to X-Plane simulator output (no measured curve for this propeller was found).
Full-throttle RPM is assumed to rise linearly from 2350 (static) to 2700 at the POH maximum speed.

Aircraft data (POH): weight 2550 lb, wing area 174 ft², span 36 ft 1 in (AR = 7.48), 180 hp at 2700 RPM, 76 in fixed-pitch propeller, stall speed 53 KCAS (gives CLmax = 1.54).

Drag parameters: CD0 = 0.0344 and e = 0.747 (Min et al., 2015).

Results

1. Minimum-drag speed: three methods agree

At 2000 m: V_md = 42.36 m/s, D_min = 1004 N.

Method	Error in V_md	Function evaluations
Analytical	exact	0
Gradient descent (4 starting points)	< 3×10⁻⁷ m/s	37–51
Differential evolution	5×10⁻⁵ m/s	225


Sources
Cessna Aircraft Company, Cessna 172S Pilot's Operating Handbook, Rev 5, 2004: geometry, engine, stall speed (Fig 5-3), climb (Fig 5-6), cruise (Fig 5-8).
S. Min, E. Harrison, H. Jimenez, D. Mavris, "Development of Aerodynamic Modeling and Calibration Methods for General Aviation Aircraft Performance Analysis – A Survey and Comparison of Models," AIAA 2015-2853: CD0 = 0.0344 and e = 0.7469 for the C172S.
Harada, Validation of Classical Aircraft Performance Methods Using Flight Test Data from a Cessna 172S, M.S. thesis, The Ohio State University [year]: build-up and simulator drag parameters, Lycoming WOT power line (Fig 15), propeller efficiency (Eq. 84).
U.S. Standard Atmosphere, 1976, NASA-TM-X-74335.
S. Gudmundsson, General Aviation Aircraft Design: Applied Methods and Procedures, 2013: Gagg–Farrar power lapse (cross-check only).
