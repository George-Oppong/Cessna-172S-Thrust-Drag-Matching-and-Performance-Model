"""
Cessna 172S Thrust-Drag Matching, Minimum-Drag Optimization and POH Validation
===============================================================================
1. Drag/thrust model with a sourced parameter set (see SOURCES below).
2. Minimum-drag speed solved three ways: analytical, gradient descent, differential evolution.
3. Validation against the POH (Cessna 172S, Rev 5): cruise %BHP, V_max, climb, service ceiling.
4. Sensitivity analysis on the parameters that are NOT known well.

SOURCES
  [POH]    Cessna 172S Pilot's Operating Handbook, Rev 5 (2004): weight, wing area, span, engine,
           propeller diameter, stall speed (Fig 5-3), climb (Fig 5-6), cruise (Fig 5-8).
  [MIN]    Min, Harrison, Jimenez, Mavris, AIAA 2015-2853: CD0 = 0.0344 (POH best-glide method,
           Table IX), K = 0.0568 / e = 0.7469 (Kroo method, Sec V.C.2).
  [HAR]    Harada, MS thesis, Ohio State: build-up CD0 = 0.0279, e = 0.699 (Sec 2.4.7); simulator
           CD0 = 0.0384, e = 0.673 (Sec 6.7); WOT power vs altitude from the Lycoming chart
           (Fig 15); propeller efficiency eta(J) (Eq. 84, a regression to X-PLANE output).
  [GF]     Gagg-Farrar power lapse as recalled (Gudmundsson 2013): cross-check only, unverified.

ASSUMPTIONS (not from a source; sensitivity-tested below)
  A1  Full-throttle RPM rises linearly from the static value to 2700 at the POH max speed.
  A2  WOT power is proportional to RPM (constant torque).
  A3  Climb IAS ~ EAS (no instrument/compressibility correction).
"""
import numpy as np
from dataclasses import dataclass, replace
from scipy.optimize import differential_evolution, brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import os
OUT = os.path.dirname(os.path.abspath(__file__))
LBF, KT, FT, HP, FT2, IN = 4.4482216, 0.514444, 0.3048, 745.7, 0.09290304, 0.0254

# ----------------------------- ISA (US Standard Atmosphere 1976, troposphere) -----------------------------
T0, P0, LAPSE, G, R_AIR = 288.15, 101325.0, 0.0065, 9.80665, 287.05
RHO0 = P0 / (R_AIR * T0)

def rho(h_m):
    T = T0 - LAPSE * h_m
    return P0 * (T / T0) ** (G / (R_AIR * LAPSE)) / (R_AIR * T)

# ----------------------------- AIRCRAFT (POH) -----------------------------
W = 2550 * LBF                       # [POH] max weight
S = 174 * FT2                        # [POH] Fig 1-1 note 4
SPAN = (36 + 1 / 12) * FT            # [POH] 36'-1"
AR = SPAN**2 / S                     # = 7.48
D_PROP = 76 * IN                     # [POH] fixed-pitch propeller diameter
VS_KCAS = 53                         # [POH] Fig 5-3, flaps up, power off
CLMAX = 2 * W / (RHO0 * S * (VS_KCAS * KT) ** 2)   # ~1.54 (UIUC value in [MIN] Table III: 1.5404)
BHP_RATED, RPM_RATED = 180.0, 2700.0                # [POH]
POH_VMAX_KT, POH_CEILING_FT, POH_BEST_GLIDE_KIAS = 126, 14000, 68

def v_stall(h_m):
    return np.sqrt(2 * W / (rho(h_m) * S * CLMAX))

# ----------------------------- MODEL -----------------------------
@dataclass(frozen=True)
class Model:
    cd0: float = 0.0344              # [MIN]
    e: float = 0.747                 # [MIN]
    eta_scale: float = 1.0           # scales [HAR] Eq. 84 (sensitivity)
    rpm_static: float = 2350.0       # [POH] static full-throttle RPM is 2300-2400
    v_rpm_rated: float = POH_VMAX_KT * KT   # A1

    @property
    def k(self):
        return 1.0 / (np.pi * AR * self.e)

    def drag(self, V, h):
        r = rho(h)
        CL = W / (0.5 * r * V**2 * S)
        return 0.5 * r * V**2 * S * (self.cd0 + self.k * CL**2)

    # --- propulsion ---
    def rpm_wot(self, V):                                   # A1
        return self.rpm_static + (RPM_RATED - self.rpm_static) * np.minimum(V / self.v_rpm_rated, 1.0)

    def eta(self, V, rpm):                                  # [HAR] Eq. 84 (X-PLANE regression)
        J = np.clip(V / (rpm / 60.0 * D_PROP), 0.0, 1.0)
        return self.eta_scale * np.clip(-3.50*J**4 + 6.71*J**3 - 5.31*J**2 + 2.83*J, 0.0, None)

    @staticmethod
    def bhp_wot(h):                                         # [HAR] Fig 15: 180 hp @ SL -> 80 hp @ 20,000 ft
        return BHP_RATED - 5.0 * (h / FT) / 1000.0

    def p_avail(self, V, h):                                # W; A2
        rpm = self.rpm_wot(V)
        return self.eta(V, rpm) * self.bhp_wot(h) * HP * rpm / RPM_RATED

    def t_avail(self, V, h):
        return self.p_avail(V, h) / V

BASE = Model()

def bhp_gagg_farrar(h):                                    # [GF] cross-check, as recalled
    s = rho(h) / RHO0
    return BHP_RATED * (s - 0.117) / (1 - 0.117)

# =============================================================================================
# MINIMUM-DRAG SPEED: THREE METHODS
# =============================================================================================
def v_md_analytical(m, h):
    """dD/dV = 0  ->  CD0 = k CL^2 : CL_md = sqrt(CD0/k), V_md = sqrt(2W/(rho S CL_md))."""
    return np.sqrt(2 * W / (rho(h) * S * np.sqrt(m.cd0 / m.k)))

def d_min_analytical(m):
    """D_min = 2 W sqrt(CD0 k): independent of altitude."""
    return 2 * W * np.sqrt(m.cd0 * m.k)

class Counted:
    def __init__(self, m, h):
        self.m, self.h, self.n = m, h, 0
    def __call__(self, V):
        self.n += 1
        return float(self.m.drag(V, self.h))

def gradient_descent(m, h, V0, tol=1e-6, max_iter=500):
    """1-D gradient descent: central-difference gradient, Armijo backtracking line search."""
    f, V, vmin = Counted(m, h), V0, 1.05 * v_stall(h)
    for it in range(max_iter):
        dv = 1e-5
        g = (f(V + dv) - f(V - dv)) / (2 * dv)
        if abs(g) < tol:
            break
        fV, a = f(V), 1.0
        while a > 1e-12:
            Vn = V - a * g
            if Vn > vmin and f(Vn) <= fV - 1e-4 * a * g**2:
                break
            a *= 0.5
        V = Vn
    return V, float(m.drag(V, h)), f.n, it + 1

def de_solve(m, h, seed=42):
    """polish=False so scipy's gradient-based L-BFGS-B polish does not contaminate the result."""
    f = Counted(m, h)
    r = differential_evolution(lambda x: f(x[0]), [(1.05 * v_stall(h), 90.0)], seed=seed, tol=1e-10, polish=False)
    return r.x[0], r.fun, f.n, r.nit

def compare_methods(m=BASE, h=2000.0):
    Vr, Dr = v_md_analytical(m, h), d_min_analytical(m)
    print(f"\n[1] MINIMUM-DRAG SPEED at h = {h:.0f} m   (analytical: V_md = {Vr:.4f} m/s, D_min = {Dr:.2f} N)")
    print(f"{'Method':<30}{'V_opt m/s':>11}{'D_opt N':>10}{'|V err|':>11}{'f-evals':>9}{'iters':>7}")
    print(f"{'Analytical':<30}{Vr:11.4f}{Dr:10.2f}{0:11.1e}{0:9d}{0:7d}")
    vs = v_stall(h)
    for mult in (1.1, 1.4, 1.8, 2.2):
        V, D, n, it = gradient_descent(m, h, mult * vs)
        print(f"{'Gradient descent (V0=%.2f Vs)' % mult:<30}{V:11.4f}{D:10.2f}{abs(V-Vr):11.1e}{n:9d}{it:7d}")
    V, D, n, it = de_solve(m, h)
    print(f"{'Differential evolution':<30}{V:11.4f}{D:10.2f}{abs(V-Vr):11.1e}{n:9d}{it:7d}")
    worst = max(abs(de_solve(m, hh)[0] - v_md_analytical(m, hh)) for hh in (0, 1000, 2000, 3000, 4000))
    print(f"  DE vs analytical, worst |V error| over 0-4000 m: {worst:.1e} m/s;  D_min is altitude-independent: {Dr:.2f} N")

# =============================================================================================
# POH DATA (transcribed from the uploaded POH, Section 5)
# =============================================================================================
CRUISE_STD = [  # (pressure alt ft, RPM, %BHP, KTAS)  Fig 5-8, standard-temperature column, 2550 lb
 (2000,2550,77,118),(2000,2500,73,115),(2000,2400,64,110),(2000,2300,57,104),(2000,2200,50,97),(2000,2100,44,90),
 (4000,2600,77,120),(4000,2550,73,117),(4000,2500,69,115),(4000,2400,61,109),(4000,2300,54,102),(4000,2200,48,96),(4000,2100,42,89),
 (6000,2650,77,122),(6000,2600,73,119),(6000,2500,65,114),(6000,2400,57,108),(6000,2300,51,101),(6000,2200,45,94),
 (8000,2700,77,124),(8000,2650,72,122),(8000,2600,68,119),(8000,2500,61,112),(8000,2400,54,106),(8000,2300,48,99),(8000,2200,43,92),
 (10000,2700,72,123),(10000,2650,68,120),(10000,2600,64,117),(10000,2500,57,111),(10000,2400,51,104),(10000,2300,46,97),
 (12000,2650,64,119),(12000,2600,61,116),(12000,2500,54,109),(12000,2400,49,102),(12000,2300,44,95)]

ROC_TEMPS = [-20, 0, 20, 40]
ROC_POH = [  # (alt ft, climb KIAS, ROC fpm at -20, 0, 20, 40 degC)  Fig 5-6, flaps up, full throttle, 2550 lb
 (0,74,[855,785,710,645]),(2000,73,[760,695,625,560]),(4000,73,[685,620,555,495]),
 (6000,73,[575,515,450,390]),(8000,72,[465,405,345,285]),(10000,72,[360,300,240,180]),(12000,72,[255,195,135,None])]

def poh_roc_isa(alt_ft, vals):
    t_isa = 15.0 - 1.98 * alt_ft / 1000.0
    pts = [(t, v) for t, v in zip(ROC_TEMPS, vals) if v is not None]
    return float(np.interp(t_isa, [p[0] for p in pts], [p[1] for p in pts]))

# =============================================================================================
# VALIDATION
# =============================================================================================
def cruise_check(m):
    """Zero-fit prediction of POH %BHP: P_shaft = D V / eta(J) (J from the tabulated RPM)."""
    pred, poh = [], []
    for h_ft, rpm, pct, kta in CRUISE_STD:
        h, V = h_ft * FT, kta * KT
        pred.append(m.drag(V, h) * V / m.eta(V, rpm) / (BHP_RATED * HP) * 100)
        poh.append(pct)
    pred, poh = np.array(pred), np.array(poh)
    err = pred - poh
    return pred, poh, err.mean(), np.sqrt((err**2).mean())

def v_max(m, h=0.0):
    V = np.linspace(35, 95, 600)
    f = m.p_avail(V, h) - m.drag(V, h) * V
    idx = np.where((f[:-1] > 0) & (f[1:] <= 0))[0]
    if len(idx) == 0:
        return np.nan
    i = idx[-1]
    return brentq(lambda v: m.p_avail(v, h) - m.drag(v, h) * v, V[i], V[i + 1])

def roc_fpm(m, V, h):
    return (m.p_avail(V, h) - m.drag(V, h) * V) / W / FT * 60.0

def climb_check(m):
    rows = []
    for alt_ft, kias, vals in ROC_POH:
        h = alt_ft * FT
        V = kias * KT / np.sqrt(rho(h) / RHO0)                # A3
        rows.append((alt_ft, poh_roc_isa(alt_ft, vals), float(roc_fpm(m, V, h))))
    return rows

def service_ceiling_ft(m):
    Vg = np.linspace(28, 65, 150)
    f = lambda h: max(roc_fpm(m, Vg[Vg > 1.1 * v_stall(h)], h)) - 100.0
    try:
        return brentq(f, 3000.0, 7000.0) / FT
    except ValueError:
        return np.nan

def validate(m=BASE):
    pred, poh, bias, rms = cruise_check(m)
    Js = [r[3] * KT / (r[1] / 60 * D_PROP) for r in CRUISE_STD]
    print(f"\n[2] VALIDATION AGAINST THE POH (CD0/e from [MIN], eta(J) from [HAR]; nothing tuned to these rows)")
    print(f"  Cruise %BHP (37 rows, Fig 5-8): mean error {bias:+.1f} points, RMS {rms:.1f} points")
    print(f"      NOTE: all rows have J = {min(Js):.2f}-{max(Js):.2f}, so this tests ONE propeller operating point")
    for h_ft in (2000, 4000, 6000, 8000, 10000, 12000):
        sel = [i for i, r in enumerate(CRUISE_STD) if r[0] == h_ft]
        print(f"      {h_ft:>6} ft: mean error {np.mean(pred[sel]-poh[sel]):+5.1f} points")
    vm = v_max(m) / KT
    print(f"  V_max at sea level: model {vm:.0f} kt vs POH {POH_VMAX_KT} KTAS ({(vm/POH_VMAX_KT-1)*100:+.1f}%)   (semi-independent: A1 uses the POH V_max as the speed where RPM reaches 2700)")
    rows = climb_check(m)
    print(f"  Max rate of climb, ISA (Fig 5-6):")
    print(f"      {'alt ft':>7}{'POH fpm':>9}{'model fpm':>11}{'error':>8}{'excess-power error':>20}")
    for a, p, mo in rows:
        print(f"      {a:7d}{p:9.0f}{mo:11.0f}{(mo/p-1)*100:+7.0f}%{(mo-p)*FT/60*W/1000:+15.1f} kW")
    mae = np.mean([abs(mo / p - 1) * 100 for _, p, mo in rows])
    print(f"      mean absolute error {mae:.0f}%")
    sc = service_ceiling_ft(m)
    print(f"  Service ceiling (ROC = 100 fpm): model {sc:,.0f} ft vs POH {POH_CEILING_FT:,} ft")
    h0 = 0.0
    LD_max = W / d_min_analytical(m)
    Vb = POH_BEST_GLIDE_KIAS * KT
    LD_b = W / m.drag(Vb, h0)
    print(f"  Best glide: model V_md = {v_md_analytical(m, h0)/KT:.1f} kt (L/D max {LD_max:.1f}); at POH {POH_BEST_GLIDE_KIAS} KIAS model L/D = {LD_b:.1f} ({(LD_b/LD_max-1)*100:+.0f}% vs max)")
    print(f"  Power lapse cross-check: Lycoming-chart line vs Gagg-Farrar (as recalled):")
    for h_ft in (8000, 14000, 20000):
        print(f"      {h_ft:>6} ft: {m.bhp_wot(h_ft*FT):6.1f} hp  vs  {bhp_gagg_farrar(h_ft*FT):6.1f} hp")
    return pred, poh, rows

def sensitivity():
    print(f"\n[3] SENSITIVITY (what the conclusions depend on)")
    cases = [("Baseline [MIN] 0.0344 / 0.747", BASE),
             ("CD0/e [HAR] build-up 0.0279/0.699", replace(BASE, cd0=0.0279, e=0.699)),
             ("CD0/e [HAR] simulator 0.0384/0.673", replace(BASE, cd0=0.0384, e=0.673)),
             ("eta(J) x 0.75 (~0.65 peak, HAR Ch.3)", replace(BASE, eta_scale=0.75)),
             ("Static RPM 2300", replace(BASE, rpm_static=2300.0)),
             ("Static RPM 2400", replace(BASE, rpm_static=2400.0))]
    print(f"  {'case':<40}{'cruise RMS':>11}{'Vmax kt':>9}{'SL ROC':>8}{'ceiling ft':>12}{'V_md kt':>9}{'D_min N':>9}")
    for name, m in cases:
        _, _, _, rms = cruise_check(m)
        print(f"  {name:<40}{rms:11.1f}{v_max(m)/KT:9.0f}{climb_check(m)[0][2]:8.0f}{service_ceiling_ft(m):12,.0f}"
              f"{v_md_analytical(m,0)/KT:9.1f}{d_min_analytical(m):9.0f}")
    tgt = climb_check(BASE)[0][1]
    c = brentq(lambda x: climb_check(replace(BASE, cd0=x))[0][2] - tgt, 0.03, 0.08)
    mc = replace(BASE, cd0=c)
    print(f"  Diagnostic: CD0 needed to match sea-level climb alone = {c:.4f}; that gives cruise RMS "
          f"{cruise_check(mc)[3]:.1f} points and V_max {v_max(mc)/KT:.0f} kt -> the climb gap is not a drag-polar error")
    print(f"  Closed-form +/-20% on the baseline (V_md at 2000 m / D_min):")
    for lab, m in (("CD0 -20%", replace(BASE, cd0=BASE.cd0*0.8)), ("CD0 +20%", replace(BASE, cd0=BASE.cd0*1.2)),
                   ("e   -20%", replace(BASE, e=BASE.e*0.8)), ("e   +20%", replace(BASE, e=BASE.e*1.2))):
        print(f"      {lab}: V_md = {v_md_analytical(m,2000.0):.2f} m/s, D_min = {d_min_analytical(m):.0f} N")

# =============================================================================================
# PLOTS
# =============================================================================================
def plot_matching(m=BASE):
    V = np.linspace(28, 80, 300)
    fig, ax = plt.subplots(figsize=(8, 6))
    for h_ft, c in ((0, "C0"), (4000, "C1"), (8000, "C2"), (12000, "C3")):
        h = h_ft * FT
        ax.plot(V, m.drag(V, h), c, label=f"Drag, {h_ft} ft")
        ax.plot(V, m.t_avail(V, h), c + "--", label=f"Thrust avail (WOT), {h_ft} ft")
        ax.plot(v_md_analytical(m, h), d_min_analytical(m), "ko", ms=4)
    ax.axvline(POH_VMAX_KT * KT, color="gray", ls=":", lw=1)
    ax.text(POH_VMAX_KT * KT + 0.3, 2300, "POH V_max\n(SL)", fontsize=8, color="gray")
    ax.set_xlabel("True airspeed (m/s)"); ax.set_ylabel("Force (N)"); ax.set_ylim(0, 2600)
    ax.set_title("Cessna 172S thrust-drag matching (dots: analytical min-drag points)")
    ax.legend(fontsize=8, ncol=2); ax.grid(True)
    fig.savefig(OUT + "thrust_drag_matching.png", dpi=150, bbox_inches="tight"); plt.close(fig)

def plot_validation(pred, poh, rows):
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    ax[0].plot([35, 85], [35, 85], "k-", lw=1)
    ax[0].plot(poh, pred, "o", ms=5)
    ax[0].set_xlabel("POH %BHP (Fig 5-8)"); ax[0].set_ylabel("Model %BHP"); ax[0].grid(True)
    ax[0].set_title("Cruise power: model vs POH (37 rows)")
    a = [r[0] for r in rows]
    ax[1].plot(a, [r[1] for r in rows], "ko-", label="POH (ISA, interpolated)")
    ax[1].plot(a, [r[2] for r in rows], "s--", label="Model")
    ax[1].set_xlabel("Pressure altitude (ft)"); ax[1].set_ylabel("Max rate of climb (fpm)"); ax[1].grid(True)
    ax[1].set_title("Climb: model vs POH"); ax[1].legend()
    fig.savefig(OUT + "validation.png", dpi=150, bbox_inches="tight"); plt.close(fig)

if __name__ == "__main__":
    print(f"Geometry: W = {W:.0f} N, S = {S:.2f} m^2, span = {SPAN:.2f} m, AR = {AR:.3f}, CLmax (from POH stall) = {CLMAX:.3f}")
    compare_methods()
    pred, poh, rows = validate()
    sensitivity()
    plot_matching(); plot_validation(pred, poh, rows)
