"""
Magnetopause Pressure Balance & Geomagnetic Perturbation Simulation
====================================================================

Author: Murad Tulunay (February 2026)
Part of: Tulunay Geomagnetic External Modulation Hypothesis

This script models the interaction between Solar Wind dynamic pressure
and Earth's magnetic field to estimate:
  1. Magnetopause standoff distance (Chapman-Ferraro)
  2. Surface magnetic perturbation (ΔB) from external compression
  3. Order-of-magnitude scale calculation (Appendix A)

Developed with AI assistance (Claude, ChatGPT, Gemini).
"""

import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Physical Constants
# =============================================================================
MU_0 = 4 * np.pi * 1e-7        # Vacuum permeability (T·m/A)
M_P = 1.6726e-27                # Proton mass (kg)
R_EARTH = 6.371e6               # Earth radius (m)
B_EQUATOR = 3.0e-5              # Earth's equatorial surface field (T) ~ 30,000 nT

# Earth's magnetic dipole moment
# From B_eq = μ₀ M / (4π R³)  →  M = 4π R³ B_eq / μ₀
EARTH_DIPOLE_MOMENT = 4 * np.pi * R_EARTH**3 * B_EQUATOR / MU_0  # ~7.8e22 A·m²


# =============================================================================
# Appendix A: Order-of-Magnitude Current Requirement
# =============================================================================
def appendix_a_scale_calculation():
    """
    Calculate the current required to generate Earth's full magnetic field
    from a single ionospheric current loop (Appendix A of the paper).

    For a circular loop of radius R carrying current I,
    the field at the center is: B = μ₀ I / (2R)
    Solving for I: I = 2RB / μ₀
    """
    R_iono = R_EARTH + 300e3     # Ionosphere F-layer altitude (~300 km)
    B_target = 50e-6             # Target field: 50 μT (full Earth field)

    I_required = (2 * R_iono * B_target) / MU_0

    print("=" * 65)
    print("APPENDIX A: Order-of-Magnitude Scale Calculation")
    print("=" * 65)
    print(f"  Ionosphere radius:     {R_iono/1e6:.2f} x 10^6 m")
    print(f"  Target field:          {B_target*1e6:.0f} uT")
    print(f"  Required current:      {I_required:.2e} A")
    print(f"                       = {I_required/1e9:.2f} GigaAmperes")
    print(f"  Observed Sq currents:  ~10^5 - 10^6 A (100 kA - 1 MA)")
    print(f"  Gap factor:            ~{I_required/1e6:.0f}x too large")
    print()
    print(f"  Observed perturbation: ~20-80 nT (Sq variation)")
    print(f"  Fraction of main field: {80/50000*100:.2f}%")
    print(f"  -> Ionosphere CANNOT generate the main field")
    print(f"  -> But CAN modulate it (~0.04-2% perturbation)")
    print("=" * 65)
    print()

    return I_required


# =============================================================================
# Magnetopause Standoff Distance (Chapman-Ferraro)
# =============================================================================
def magnetopause_distance(v_sw, n_sw):
    """
    Calculate magnetopause standoff distance using pressure balance.

    Chapman-Ferraro model:
        Solar wind dynamic pressure = Magnetic pressure at magnetopause
        (1/2) rho v^2 = B^2 / (2 mu_0)

    For a dipole, B at distance R on the equator:
        B(R) = mu_0 M / (4 pi R^3)

    Solving for R_mp:
        R_mp = [ mu_0 M^2 / (8 pi^2 P_dyn) ]^(1/6)

    Parameters:
        v_sw : float - Solar wind velocity (m/s)
        n_sw : float - Solar wind proton density (particles/m^3)

    Returns:
        R_mp : float - Magnetopause distance in Earth radii
    """
    # Dynamic pressure: P_dyn = (1/2) n m_p v^2
    P_dyn = 0.5 * n_sw * M_P * v_sw**2

    # Standoff distance from pressure balance
    R_mp = (MU_0 * EARTH_DIPOLE_MOMENT**2 / (8 * np.pi**2 * P_dyn))**(1.0/6.0)

    return R_mp / R_EARTH  # Return in Earth radii


def surface_perturbation_dst(R_mp_RE):
    """
    Estimate Dst-equivalent surface perturbation (nT) from magnetopause compression.

    Uses empirical scaling calibrated to the May 2024 superstorm:
      May 2024: R_mp ~ 5.04 R_E, Dst ~ -412 nT
      Quiet:    R_mp ~ 10 R_E,   Dst ~ 0 nT

    The perturbation scales linearly with compression ratio (R_quiet/R_mp - 1),
    calibrated so that compression ratio 2 (10/5) yields ~412 nT.
    """
    R_quiet = 10.0  # Quiet-time magnetopause ~10 R_E

    if R_mp_RE >= R_quiet:
        return 0.0

    compression = R_quiet / R_mp_RE
    delta_dst = 412.0 * (compression - 1)  # Calibrated to May 2024

    return delta_dst


# =============================================================================
# Simulation Scenarios
# =============================================================================
def run_scenarios():
    """Run the main simulation scenarios from the paper."""

    scenarios = [
        ("Quiet Day",           400e3,  5e6,   "Normal conditions"),
        ("Moderate Storm (G2)", 600e3,  15e6,  "Moderate activity"),
        ("Severe Storm (G4)",   800e3,  30e6,  "Strong CME impact"),
        ("May 2024 (G5)",       900e3,  50e6,  "Actual event estimate"),
        ("Carrington-class",    1500e3, 80e6,  "Historical extreme"),
    ]

    print("=" * 95)
    print("SIMULATION RESULTS: Magnetopause Distance & Surface Perturbation")
    print("=" * 95)
    print(f"  {'Scenario':<23} {'V_sw':>9} {'n':>10} {'R_mp':>10} {'Est. Dst':>10}   {'Note'}")
    print(f"  {'':23} {'(km/s)':>9} {'(cm^-3)':>10} {'(R_E)':>10} {'(nT)':>10}")
    print("-" * 95)

    results = []
    for name, v, n, note in scenarios:
        R_mp = magnetopause_distance(v, n)
        dst = surface_perturbation_dst(R_mp)
        n_cm3 = n / 1e6
        v_kms = v / 1e3

        print(f"  {name:<23} {v_kms:>8.0f}  {n_cm3:>9.0f}  {R_mp:>9.2f}  {dst:>9.0f}    {note}")
        results.append((name, v_kms, n_cm3, R_mp, dst))

    print("-" * 95)
    print()
    print("  Validation against observations:")
    print(f"    May 2024 observed:  Magnetopause = 5.04 R_E, Dst = -412 nT (Hayakawa et al., 2024)")
    print(f"    Carrington 1859:    Estimated Dst = -850 to -1050 nT")
    print()

    return results


# =============================================================================
# Visualization
# =============================================================================
def plot_results(results):
    """Generate publication-quality plots."""

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    names = [r[0] for r in results]
    R_mp_values = [r[3] for r in results]
    dst_values = [r[4] for r in results]

    colors = ['#2ecc71', '#f1c40f', '#e74c3c', '#c0392b', '#8e44ad']

    # --- Plot 1: Magnetopause distance ---
    ax1 = axes[0]
    bars1 = ax1.bar(names, R_mp_values, color=colors, edgecolor='black', linewidth=0.5)
    ax1.axhline(y=6.6, color='orange', linestyle=':', linewidth=1.5,
                label='Geosynchronous orbit (6.6 R_E)')
    ax1.axhline(y=5.04, color='red', linestyle='--', linewidth=1.5,
                label='May 2024 observed (5.04 R_E)')
    ax1.set_ylabel('Magnetopause Distance (Earth Radii)', fontsize=12)
    ax1.set_title('Magnetopause Compression\nUnder Solar Wind Pressure', fontsize=13, fontweight='bold')
    ax1.legend(fontsize=9, loc='upper right')
    ax1.set_ylim(0, 14)
    ax1.tick_params(axis='x', rotation=30, labelsize=9)

    for bar, val in zip(bars1, R_mp_values):
        ax1.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.3,
                f'{val:.1f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

    # --- Plot 2: Estimated Dst perturbation ---
    ax2 = axes[1]
    bars2 = ax2.bar(names, dst_values, color=colors, edgecolor='black', linewidth=0.5)
    ax2.axhline(y=412, color='red', linestyle='--', linewidth=1.5,
                label='May 2024 actual: Dst = -412 nT')
    ax2.set_ylabel('Estimated |Dst| Perturbation (nT)', fontsize=12)
    ax2.set_title('Surface Magnetic Perturbation\nFrom External Compression', fontsize=13, fontweight='bold')
    ax2.legend(fontsize=9)
    ax2.tick_params(axis='x', rotation=30, labelsize=9)

    for bar, val in zip(bars2, dst_values):
        ax2.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 15,
                f'{val:.0f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

    plt.tight_layout()
    plt.savefig('simulation_results.png', dpi=150, bbox_inches='tight')
    print("  -> Plot saved: simulation_results.png")
    plt.close()

    # --- Plot 3: Continuous magnetopause curve ---
    fig2, ax3 = plt.subplots(figsize=(10, 6))

    v_range = np.linspace(300e3, 1600e3, 300)
    n_range = 5e6 + (80e6 - 5e6) * ((v_range - 300e3) / (1600e3 - 300e3))**1.3

    R_mp_curve = [magnetopause_distance(v, n) for v, n in zip(v_range, n_range)]

    ax3.plot(v_range/1e3, R_mp_curve, 'b-', linewidth=2.5, label='Modeled R_mp')
    ax3.axhline(y=10.0, color='green', linestyle=':', alpha=0.6, label='Quiet-time (~10 R_E)')
    ax3.axhline(y=6.6, color='orange', linestyle=':', linewidth=1.5, label='Geosynchronous (6.6 R_E)')
    ax3.axhline(y=5.04, color='red', linestyle='--', linewidth=1.5, label='May 2024: 5.04 R_E')
    ax3.axhspan(0, 6.6, alpha=0.08, color='red')
    ax3.set_xlabel('Solar Wind Velocity (km/s)', fontsize=12)
    ax3.set_ylabel('Magnetopause Standoff Distance (Earth Radii)', fontsize=12)
    ax3.set_title('Magnetopause Distance vs Solar Wind Speed\n'
                  '(Tulunay External Modulation Hypothesis)', fontsize=13, fontweight='bold')
    ax3.legend(fontsize=10, loc='upper right')
    ax3.set_ylim(0, 14)
    ax3.set_xlim(300, 1600)
    ax3.grid(True, alpha=0.3)
    ax3.annotate('Satellite danger zone', xy=(1400, 4), fontsize=10,
                color='red', alpha=0.6, ha='center')

    plt.tight_layout()
    plt.savefig('magnetopause_curve.png', dpi=150, bbox_inches='tight')
    print("  -> Plot saved: magnetopause_curve.png")
    plt.close()

    # --- Plot 4: Perturbation percentage of main field ---
    fig3, ax4 = plt.subplots(figsize=(10, 5))

    dst_curve = [surface_perturbation_dst(r) for r in R_mp_curve]
    pct_curve = [d / 50000 * 100 for d in dst_curve]

    ax4.fill_between(v_range/1e3, pct_curve, alpha=0.3, color='purple')
    ax4.plot(v_range/1e3, pct_curve, 'purple', linewidth=2, label='External perturbation (% of main field)')
    ax4.axhline(y=0.82, color='red', linestyle='--', linewidth=1,
                label='May 2024: ~0.82% (412 nT / 50,000 nT)')
    ax4.set_xlabel('Solar Wind Velocity (km/s)', fontsize=12)
    ax4.set_ylabel('Perturbation as % of Main Field', fontsize=12)
    ax4.set_title('External Contribution to Total Geomagnetic Field\n'
                  'Key Question: Is this enough to steer secular variation?',
                  fontsize=13, fontweight='bold')
    ax4.legend(fontsize=10)
    ax4.grid(True, alpha=0.3)
    ax4.set_xlim(300, 1600)

    plt.tight_layout()
    plt.savefig('perturbation_percentage.png', dpi=150, bbox_inches='tight')
    print("  -> Plot saved: perturbation_percentage.png")
    plt.close()


# =============================================================================
# Summary Table for Paper
# =============================================================================
def print_paper_table():
    """Print the results table formatted for the whitepaper."""

    print("=" * 75)
    print("TABLE FOR PAPER (Section 7: Perturbation Analysis)")
    print("=" * 75)
    print()
    print("  Solar Wind      | V_sw    | Standoff  | Est. dB  | % of Main")
    print("  Regime           | (km/s)  | (R_E)     | (nT)     | Field    ")
    print("  -----------------+---------+-----------+----------+----------")

    table_data = [
        ("Quiet",       400e3,  5e6),
        ("Moderate G2",  600e3,  15e6),
        ("Severe G4",    800e3,  30e6),
        ("May 2024 G5",  900e3,  50e6),
        ("Carrington",   1500e3, 80e6),
    ]

    for name, v, n in table_data:
        R = magnetopause_distance(v, n)
        dst = surface_perturbation_dst(R)
        pct = dst / 50000 * 100
        print(f"  {name:<17}| {v/1e3:>6.0f}  | {R:>8.2f}  | {dst:>7.0f}  | {pct:>6.2f}%")

    print()
    print("  May 2024 observed: R_mp = 5.04 R_E, Dst = -412 nT (0.82%)")
    print("=" * 75)
    print()


# =============================================================================
# Main
# =============================================================================
if __name__ == "__main__":
    print()
    print("=" * 66)
    print("  Tulunay Geomagnetic External Modulation Hypothesis")
    print("  Magnetopause Pressure Balance Simulation")
    print("  Author: Murad Tulunay (February 2026)")
    print("=" * 66)
    print()

    print(f"  Earth dipole moment: {EARTH_DIPOLE_MOMENT:.2e} A.m^2")
    print(f"  (Literature value:   ~7.94 x 10^22 A.m^2)")
    print()

    # 1. Scale calculation (Appendix A)
    appendix_a_scale_calculation()

    # 2. Main simulation scenarios
    results = run_scenarios()

    # 3. Paper-ready table
    print_paper_table()

    # 4. Generate plots
    try:
        plot_results(results)
        print("\n  All simulations complete.")
        print("  Output files:")
        print("    - simulation_results.png      : Bar charts (R_mp and Dst)")
        print("    - magnetopause_curve.png      : Continuous R_mp vs V_sw")
        print("    - perturbation_percentage.png  : External % of main field")
    except Exception as e:
        print(f"\n  Plot generation failed ({e}), but numerical results above are valid.")

    print()
