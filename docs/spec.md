# Specification: Jumeau numérique immunitaire (Bilan Immuno Twin)

## Problem Statement

From the user's perspective: a clinician has a patient's bilan immunologique — a flow cytometry snapshot of their immune system at a single point in time. They want to know what will happen to that patient's immune system if a specific immunotherapy is administered. Currently there is no way to answer this question: each patient is different, and the only way to know how they'll respond is to treat them and observe, which is costly, slow, and sometimes harmful.

The problem is building a computational model that takes a single cytometric snapshot and predicts — per patient — how the immune system will evolve in response to a perturbation, and to flag patients whose response is differential or unexpected.

## Solution

A mechanistic digital twin of the immune system (the **jumeau numérique immunitaire**), built as a Python CLI package. It:

1. Reads a patient's bilan immunologique (cytometry cell populations) as the initial state.
2. Represents the immune interactions between populations as a **graphe de régulations** (fixed topology from immunology literature, patient-specific edge weights).
3. Simulates the temporal evolution of the populations via a system of ODEs using **fonctions de Hill** for saturating interactions and **variables quasi-statiques** for cytokines.
4. Applies a **perturbation** (immunotherapy) modeled as a modification of kinetic parameters.
5. **Calibre** patient-specific kinetic parameters via **MAP** estimation with bootstrap confidence intervals.
6. Outputs time trajectories of populations plus a **score de réponse**, consumed through the CLI.

The user selects the use-case: predicting pre-treatment response, or diagnosing why a patient does not respond.

## User Stories

1. As a clinician, I want to provide a patient's bilan immunologique (cytométrie), so that the jumeau can be initialized from that specific patient's immune state.

2. As a clinician, I want the jumeau to accept a single cytometric snapshot (no longitudinal history required), so that the model works even when only one bilan is available.

3. As a clinician, I want the jumeau to represent the immune populations from my panel as distinct nodes (CD8+ T, CD4+ Th subsets, B, NK, Treg, monocytes), so that the simulation reflects the populations I can actually measure.

4. As a researcher, I want the graphe de régulations to encode known immunology (activation/inhibition edges) with a fixed topology from the literature, so that the model's structure is biologically grounded rather than learned from scratch.

5. As a researcher, I want the edge weights of the graphe de régulations to be patient-specific, so that each patient's twin captures their individual immune reactivity.

6. As a researcher, I want the Th subsets (Th1, Th2, Th17) to be separate nodes in the graphe, so that distinct pro- and anti-inflammatory behaviors are modeled rather than lumped into a single CD4+ population.

7. As a clinician, I want to specify an immunotherapy perturbation (e.g., anti-PD1, anti-TNF, corticoïde), so that the jumeau simulates the immune response to that specific treatment.

8. As a clinician, I want the perturbation to be expressed as a documented modification of kinetic parameters, so that the treatment's effect is interpretable ("anti-PD1 increases CD8+ activation").

9. As a researcher, I want the ODE system to use fonctions de Hill for regulatory interactions, so that biological thresholds and saturation effects are represented.

10. As a researcher, I want the ODE solver to handle stiffness (raideur) via an implicit method (Radau/BDF), so that the system remains numerically stable despite differing timescales.

11. As a researcher, I want cytokines (IFN-γ, IL-10, TNF-α, IL-6, IL-4, IL-17) modeled as variables quasi-statiques, so that fast-evolving signals are computed algebraically without their own ODEs.

12. As a clinician, I want the simulation to run from the bilan time to a configurable horizon (with an adaptive option that detects reaching an steady state), so that I can see trajectories over a clinically relevant window.

13. As a clinician, I want the jumeau to output time-trajectories of each population with confidence intervals, so that I can judge both the trend and the uncertainty of the prediction.

14. As a clinician, I want a score de réponse derived from the simulation (not a hard binary), so that I can interpret the predicted response with nuance rather than a false-precision yes/no.

15. As a clinician, I want to use the jumeau to predict whether a patient will respond to a treatment before it is prescribed, so that I can make an informed treatment decision.

16. As a clinician, I want to use the jumeau to diagnose why a patient does not respond to a treatment, so that I can identify the mechanism of non-response.

17. As a clinician, I want the jumeau to flag differential or unexpected response trajectories, so that I am alerted to patients whose immune reaction diverges from typical expectations.

18. As a researcher, I want to calibrate patient-specific kinetic parameters from the cytometric data using MAP estimation with literature-derived priors, so that the twin reflects the individual patient without requiring extensive longitudinal data.

19. As a researcher, I want the calibration to produce bootstrap confidence intervals, so that the reliability of the parameter estimates is quantified.

20. As a developer, I want all interactions with the jumeau through a single Python CLI, so that the tool is scriptable, reproducible, and CI-friendly.

21. As a developer, I want the CLI to read a bilan immunologique and a perturbation specification, and emit structured output (JSON trajectories + score), so that results are machine-parseable.

22. As a developer, I want the CLI to emit CSV/JSON, so that results can be analyzed in downstream tooling.

23. As a researcher, I want a self-consistency validation mode (fixed parameters → simulate with noise → recalibrate → recover parameters), so that the framework can be validated on synthetic data before clinical use.

24. As a researcher, I want unit tests that verify mathematical correctness of the fonctions de Hill and the ODE system, so that numerical errors are caught at the lowest seam.

25. As a developer, I want integration tests that drive the full pipeline (bilan → calibration → simulation → score) through the CLI, so that the end-to-end flow is verified.

26. As a researcher, I want regression tests that pin simulation outputs across versions, so that model changes that silently alter results are caught.

27. As a clinician, I want the ability to simulate "what-if" scenarios to explore how different perturbations affect a patient, so that I can reason about treatment alternatives.

## Implementation Decisions

The MVP is a **CLI-only Python package** (ADR-0004). No notebooks, no REST API.

**Modules** (interfaces — not file paths):
- **Model module**: the ODE right-hand side. Owns the state vector (populations as variables d'état), uses fonctions de Hill for regulatory interactions (ADR-0001), computes cytokines as variables quasi-statiques (ADR-0005). Exposes a function that takes current state + kinetic parameters and returns the time derivative.
- **Graph module**: defines the graphe de régulations — nodes (populations including separate Th1/Th2/Th17, CD8+, B, NK, Treg, monocytes), fixed topology from literature, patient-specific edge weights.
- **Perturbation module**: maps an immunotherapy to modifications of kinetic parameters (ADR-0003). Contract: `perturbation → {parameter: relative_change}` and a documented catalogue (anti-PD1 raises CD8+ activation, anti-TNF reduces TNF-α effect, corticoïdes depress multiple pro-inflammatory signals).
- **Calibration module**: MAP estimation with literature priors, bootstrap confidence intervals (ADR-0002). Produces patient-specific kinetic parameters from a bilan.
- **Simulation module**: integrates the ODE system with an implicit solver (`solve_ivp` Radau) over a configured horizon, with an adaptive option that stops at steady state. Emits trajectories.
- **Response module**: derives the score de réponse (continuous) from the trajectories; not a binary classifier.
- **CLI module**: the single entry point. Subcommands for calibration, simulation, and scoring, reading a bilan + perturbation and emitting JSON/CSV.

**Architectural decisions carried from ADRs** (all respected in this spec):
- ADR-0001: Hill functions (not mass-action, not Michaelis-Menten) for regulatory interactions.
- ADR-0002: MAP + bootstrap (not full MCMC) for calibration bayésienne.
- ADR-0003: perturbations modify kinetic parameters (not forcing terms, not topology changes).
- ADR-0004: CLI-only interface.
- ADR-0005: cytokines as variables quasi-statiques.

**Solver**: `scipy.integrate.solve_ivp` with Radau, chosen for stiffness handling.

**Data for MVP**: no real patient data required. Validation starts with synthetic data via self-consistency.

## Testing Decisions

A good test verifies **external behavior** (what the jumeau does), not implementation details (how it does it). For a scientific model, the highest-value tests assert mathematical and numerical correctness and end-to-end pipeline fidelity.

**Seams** (confirmed with the user):
- **Primary seam — the CLI**: integration tests drive the full pipeline (bilan → calibration → simulation → score) through the CLI and assert on the emitted JSON/CSV output (trajectories, score) and process exit behavior. This is the highest seam and the one a user actually interacts with.
- **Secondary seam — pure-math internals**: unit tests assert mathematical correctness of the fonctions de Hill (limits, saturation, thresholds) and the ODE system right-hand side (e.g., conservation, sign of derivatives, zero-stability) below the CLI.

**Modules tested**:
- Model (Hill functions, derivative correctness).
- Calibration (MAP recovers known parameters on synthetic data; bootstrap intervals contain the truth).
- Simulation (solver stability, steady-state detection, horizon handling).
- CLI end-to-end (calibration → simulation → score pipeline, output schema).

**Prior art**: none in this repo yet — greenfield. The testing approach should follow the repo's conventions once the project scaffolding (test runner choice) is established by the implementing agent. Synthetic validation follows the self-consistency pattern (fixed parameters → simulate → add noise → recalibrate → recover), which becomes the reference regression fixture.

## Out of Scope

- Real clinical patient data / hospital cohort integration (starts with synthetic data only).
- Full MCMC / variational Bayesian inference (MVP uses MAP + bootstrap only).
- Cytokines as full dynamic state variables (kept as variables quasi-statiques per ADR-0005).
- Non-immunotherapy perturbations (infections, exogeneous cytokines) in the MVP.
- Multi-context / dashboard / API surface (CLI only, ADR-0004).
- Molecular-level signaling cascades (NF-κB, receptor dynamics) — out of scope by design (population-level granularity).
- Automatic perturbation catalogue expansion — only the opening therapies are implemented in the MVP.

## Further Notes

- This is the very first piece of work in the repository; the implementing agent must first scaffold the Python package and choose the test runner, following the CLI-only constraint in ADR-0004.
- The bilan is a single cytometric snapshot; the model intentionally works without longitudinal data, which is the key constraint driving the calibration approach (MAP with priors rather than data-hungry fitting).
- The score de réponse is deliberately continuous, not a binary response classification, to avoid false precision for a clinical user.
- Validation sequencing (self-consistency → public cohort → clinical data) is a later phase and only the first step (synthetic self-consistency) is in scope here.
