# Ethics Statement for TFM-Chemistry

## Synthetic Data Disclaimer
All figures, datasets, and results in this repository are **theoretical mock-ups** generated from analytical equations (Eqs. 1–7 of the paper). They serve to:
- Illustrate the Time Field Model (TFM) framework  
- Guide future experimental/computational validation  
- Provide parameterized benchmarks for HPC implementations  

No actual high-performance computing (HPC) simulations or laboratory experiments were conducted to generate these results.

---

## Reproducibility & Transparency
### Mock Data Generation
- Scripts in [`mock_data/`](mock_data/) produce synthetic datasets via simplified analytical models (e.g., exponential decay, cosine oscillations).  
- Parameters (Γ_chem, λβ², etc.) are defined in [`src/utils/tfm_parameters.py`](src/utils/tfm_parameters.py).

### Validation Limitations
- Bond energy comparisons in [`validation/`](validation/) use idealized quantum-chemical references, not real experimental data.  
- Rydberg shifts in `vs_rydberg_data.csv` are extrapolated from Eq. 2, not spectroscopic measurements.

---

## Call for Independent Validation
We urge researchers to:
1. **Test TFM predictions experimentally**:  
   - Compare high-\(n\) Rydberg spectra (e.g., hydrogen/cesium) against Eq. 2 deviations.  
   - Measure ultra-cold reaction rates for Eq. 6 oscillatory signatures.  
2. **Implement HPC-scale PDE solvers**:  
   - Extend [`src/utils/theoretical_models.py`](src/utils/theoretical_models.py) with FEniCS/PETSc solvers.  
   - Validate multi-atom coherence (Eq. 7) on clusters (guide: [`docs/future_hpc_guide.md`](docs/future_hpc_guide.md)).  

---

## Licensing & Attribution
- **Code**: MIT License ([LICENSE](LICENSE))  
- **Data**: CC-BY-4.0 (mock data)  
- **Citation**: Required if using this repository for derivative work (see [README.md#citation](README.md#citation)).  

---
## This Repository Contains Theoretical Frameworks, Not Real Simulations
- **All data is synthetic**: Generated from simplified equations (no HPC/experiments).  
- **Figures are illustrative**: Use `mock_data/` scripts to regenerate them.  
- **Validation needed**: Parameters in `src/params/` must be tested in real simulations.  

★ **Your Role**: Implement `templates/` on HPC clusters and compare to experiments!
---

## Contact  
For questions or collaboration proposals:  
📧 [alifayyaz@live.com](mailto:alifayyaz@live.com)  
📄 [Full Paper Preprint](INSERT_LINK_HERE)  
---
## Core Disclaimer  
**This repository does NOT contain**:  
- Results from real HPC simulations  
- Experimentally validated data  
- Actual PDE solutions for multi-atom systems  

**What it DOES provide**:  
- Parameterized equations (Eqs. 1–7) for others to implement  
- Mock data scripts to illustrate theoretical predictions  
- Templates for future HPC/experimental validation  

**Figures 1–5 are illustrative** and generated from simplified analytical models.  
