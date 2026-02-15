# Tulunay Geomagnetic External Modulation Hypothesis

**External Modulation of Geomagnetic Secular Variation: A Hypothesis on Solar Wind Coupling, The May 2024 Storm, and the Tesla Circuit**

**Author:** Murad Tulunay (B.Sc. Computer Engineering, Istanbul Technical University)  
**Date:** February 2026  
**Contact:** murad.tulunay@gmail.com

---

## Summary

This repository contains a hypothesis paper and supporting simulation code proposing that the Solar Wind's dynamic pressure and ionospheric day/night conductivity asymmetry act as a significant **external steering mechanism** for Earth's magnetic pole position and rapid geomagnetic impulses ("jerk-like" events).

**Core idea:** Earth's magnetosphere is not a passive shield — it is an active, externally-driven current system that modulates the geomagnetic vector observed at the surface.

### Key Findings (Simulation)

| Solar Scenario | Wind Speed | Surface ΔB | Equivalent Effect |
|---|---|---|---|
| Quiet | 400 km/s | ~20–30 nT | Minor diurnal variation |
| Storm (G4/G5) | 800 km/s | ~100–300 nT | Significant jerk-like signature |
| Extreme (Carrington) | 1500 km/s | >1000 nT | Major vector deviation |
| **May 2024 actual** | ~800–1000 km/s | **Dst ≈ −412 nT** | Magnetopause at 5.04 Rₑ |

### Scale Honesty (Appendix A)

Generating Earth's full ~50 μT field from a single ionospheric current loop would require ~0.5 GigaAmperes. Observed ionospheric currents are ~10⁵–10⁶ A. **The ionosphere cannot generate the main field, but it is perfectly scaled to modulate it** (~0.02–2% perturbation depending on space weather).

---

## Repository Contents

```
├── README.md                          # This file
├── paper/
│   └── External_Modulation_Hypothesis_v2.pdf   # Full paper (PDF)
│   
├── simulation/
│   └── magnetopause_pressure_balance.py        # Python simulation code
└── LICENSE
```

## Call for Falsification

This is a **falsifiable proposition**. We invite teams with access to ESA Swarm data and/or ground magnetometer arrays to test:

- **Test A:** Correlate OMNIWeb Solar Wind dynamic pressure with historical geomagnetic jerk events
- **Test B:** Separate internal/external field harmonics using Swarm A/B/C vector data  
- **Test C:** Compare Earth's external field response to Mars's induced magnetic environment

Even a negative result would be valuable — it would bound external contributions and clarify what cannot be explained without core dynamics.

## How to Run the Simulation

```bash
pip install numpy matplotlib
cd simulation/
python magnetopause_pressure_balance.py
```

## Acknowledgements

Developed with the assistance of AI research tools (Anthropic Claude, OpenAI ChatGPT, and Google Gemini) for data synthesis, simulation logic, literature review, and document drafting. Responsibility for all claims and interpretation remains with the author.

Dedicated to the engineering intuition of Nikola Tesla, and to every curious mind that asks "what if?"

## License

This work is licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). You are free to share and adapt with attribution.

## Citation

If you reference this work:

```
Tulunay, M. (2026). External Modulation of Geomagnetic Secular Variation: 
A Hypothesis on Solar Wind Coupling, The May 2024 Storm, and the Tesla Circuit. 
GitHub: https://github.com/MuradTulunay/tulunay-geomagnetic-external-modulation-hypothesis
```
