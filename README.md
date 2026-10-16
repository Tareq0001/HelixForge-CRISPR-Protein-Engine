# HelixForge 3D | Computational Genomics, CRISPR-Cas9 Studio & Protein Folding Engine

[![Live Workstation](https://img.shields.io/badge/Live_Workstation-HelixForge_3D-10b981?style=for-the-badge&logo=dna)](https://tareq0001.github.io/HelixForge-CRISPR-Protein-Engine/)
[![Architecture](https://img.shields.io/badge/Engine-Molecular_Dynamics_%26_CRISPR_Cas9-06b6d4?style=for-the-badge)](https://github.com/Tareq0001/HelixForge-CRISPR-Protein-Engine)

**HelixForge 3D** is an interactive computational biology, CRISPR-Cas9 genome editing, and 3D protein folding molecular dynamics workstation engineered by **Tareq Aboushi**. It unites real-time 3D polypeptide conformational folding (Lennard-Jones 12-6 & hydrogen-bond energy minimization), live Ramachandran $(\phi, \psi)$ dihedral phase-space plotting, CRISPR-Cas9 sgRNA PAM scanning & thermodynamic scoring, and Smith-Waterman dynamic programming sequence alignment in a single zero-dependency platform.

---

## Key Architectural Pillars

1. **Ab-Initio 3D Protein Folding & Molecular Dynamics (`src/molecular_dynamics.py`)**
   - Real-time 3D Canvas/WebGL polypeptide renderer supporting **Ribbon + Secondary Structure**, **CPK Space-Filling Van der Waals Spheres**, **Ball-and-Stick All-Atom**, and **Electrostatic Charge Mapping**.
   - **Lennard-Jones 12-6 Potential & Langevin Thermostat**:
     $$V_{\text{LJ}}(r_{ij}) = 4\varepsilon \left[ \left(\frac{\sigma}{r_{ij}}\right)^{12} - \left(\frac{\sigma}{r_{ij}}\right)^6 \right] + \sum_{\text{H-bonds}} E_{\text{HB}}(\theta, r)$$
   - Trigger **Denature (Thermal Unfold)** and watch the polypeptide collapse down its thermodynamic folding funnel ($\Delta G$ in kcal/mol) into its native minimum-RMSD conformation.
2. **Real-Time Ramachandran $(\phi, \psi)$ Dihedral Phase Plot**
   - Tracks backbone torsion angles $\phi \in [-180^\circ, 180^\circ]$ and $\psi \in [-180^\circ, 180^\circ]$ across right-handed $\alpha$-helices, $\beta$-sheets, and left-handed $\alpha_L$ regions.
3. **CRISPR-Cas9 Guide RNA (sgRNA) Scanner & Base Editor (`src/crispr_engine.py`)**
   - Scans genomic loci (*BRCA1*, *HBB Sickle-Cell $\beta$-Globin*, *TP53*, *PCSK9*) for `5'-NGG-3'` SpCas9 and `5'-TTTV-3'` Cas12a PAM motifs.
   - Computes **SantaLucia Nearest-Neighbor Duplex Melting Temperature ($T_m$)**, **Doench Rule Set 2 On-Target Cleavage Efficiency**, and **CFD Off-Target Specificity**.
   - Interactive **Single-Base Mutator & Cas9 DSB Cut Simulator** with live mRNA transcription, 64-codon ribosome translation, and Kyte-Doolittle hydropathy profile.
4. **Smith-Waterman Dynamic Programming Alignment Matrix (`src/sequence_alignment.py`)**
   - Interactive local/global sequence alignment matrix with traceback visualization.

---

## Live Interactive Workstation
Launch directly in your browser: **[https://tareq0001.github.io/HelixForge-CRISPR-Protein-Engine/](https://tareq0001.github.io/HelixForge-CRISPR-Protein-Engine/)**
