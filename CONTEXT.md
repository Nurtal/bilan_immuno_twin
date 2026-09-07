# Bilan Immuno Twin — Digital Twin of the Immune System

Mechanistic digital twin initialized from a patient's immunological profile (cytometry), simulating immune population dynamics via ODEs to predict differential responses to immunotherapies.

## Language

**Bilan immunologique**:
Patient's immune status snapshot at a given time, measured exclusively by flow cytometry (cell populations and activation markers).
_Avoid_: bilan sérologique, dosages cytokiniques

**Jumeau (numérique immunitaire)**:
Computational model representing a specific patient's immune system state and its dynamic response to perturbations.
_Avoid_: modèle statistique, prédicteur

**Population immunitaire**:
A quantified cell subset measured by cytometry (CD4+ T, CD8+ T, B, NK, Treg, monocytes, etc.).
_Avoid_: sous-population, lignée cellulaire

**Panel cytometrique**:
The set of flow cytometry markers and populations measured in the immunological profile.
_Avoid_: panel, marker set

**Graphe de régulations**:
Directed graph whose nodes are immune populations and edges encode known regulatory interactions (activation/inhibition). Topology is fixed from immunological literature; edge weights vary per patient.
_Avoid_: réseau, matrice de corrélation

**Perturbation**:
An exogenous intervention applied to the model (immunotherapy, infection, cytokine) that modifies the ODE dynamics. Encoded as a modification of kinetic parameters (e.g., anti-PD1 increases CD8+ activation rate).
_Avoid_: traitement, stimulus

**Calibration bayésienne**:
Process of estimating patient-specific kinetic parameters from cytometry data using Bayesian inference with literature-derived priors. MVP uses MAP estimation with bootstrap confidence intervals.
_Ajustement paramétrique_: calibration

**Score de réponse**:
Quantitative prediction of patient's likely response to a given immunotherapy, derived from the twin's simulation.
_Avoid_: prédiction binaire, diagnostic

**Fonction de Hill**:
Mathematical function used to model saturating interactions in the ODEs: `f(x) = x^n / (K^n + x^n)`. Captures biological thresholds and saturation effects without excessive complexity.
_Avoid_: cinétique de Michaelis-Menten (for regulatory interactions), mass-action bilinéaire

**Stiffness ( raideur)**:
Property of ODE systems where different variables evolve on vastly different timescales, requiring implicit solvers (Radau, BDF) rather than explicit ones (RK45).
_Avoid_: non-stiff

**MAP (Maximum A Posteriori)**:
Point estimate of model parameters that maximizes the posterior probability given patient data and priors. Used as the calibration method for the MVP.
_Avoid_: MCMC complet (for the MVP)

**Variable quasi-statique**:
A model variable computed algebraically from other state variables at each timestep, without its own ODE. Used for fast-evolving species (cytokines) whose dynamics are negligible on the simulation timescale.
_Avoid_: variable d'état (for cytokines), paramètre fixe

**Sous-population Th**:
CD4+ T helper subsets classified by function: Th1 (pro-inflammatory, IFN-γ), Th2 (anti-inflammatory, IL-4/IL-10), Th17 (pro-inflammatory, IL-17). Each is a separate node in the regulation graph.
_Avoid_: Th global, CD4+ non-spécifié
