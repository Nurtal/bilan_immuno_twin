# Calibration bayésienne par MAP + bootstrap

L'estimation des paramètres cinétiques patient-specific se fait par MAP (Maximum A Posteriori) avec des intervalles de confiance par bootstrap, plutôt que par MCMC complet ou inférence variationnelle.

## Contexte

Le jumeau nécessite d'estimer des paramètres cinétiques (taux de prolifération, de mort, d'activation) pour chaque patient à partir de données cytometriques. Les données cliniques sont souvent limitées en nombre de points temporels, rendant le MCMC complet coûteux et potentiellement instable.

## Considérées

- **MCMC complet (Stan/PyMC)** : gold standard bayésien, mais nécessite des centaines d'échantillons et un temps de calcul important.
- **Inférence variationnelle (VI)** : plus rapide, mais approximation de la-posterior qui peut être biaisée.
- **MAP + bootstrap** : estimation ponctuelle + intervalles de confiance par rééchantillonnage. Suffisant pour un MVP, interprétable, rapide.

## Décision

MAP + bootstrap pour l'MVP. Les priors sont tirés de la littérature immunologique. Les intervalles de confiance par bootstrap (100-500 rééchantillons) donnent une estimation de l'incertitude. Migration vers MCMC complet envisagée pour la phase clinique.
