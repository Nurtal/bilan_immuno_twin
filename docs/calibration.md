# Calibration bayésienne (MAP + bootstrap)

Le jumeau est rendu patient-spécifique en estimant les **taux de croissance
basaux de chaque population** (`growth:<population>`) à partir d'un seul bilan
(une seule snapshot cytométrique — pas d'historique longitudinal requis).

## Hypothèse de modélisation

Le bilan observé `y` est (approximativement) un **état stationnaire** du graphe
de régulations. À l'équilibre (`dx_i/dt = 0`), la condition se résout en forme
fermée pour chaque population non nulle :

```
g*_i = death_i / ( (basal_i + act_i(y)) · damp_i(y) · (1 − y_i / M_i) )
```

où `act_i(y)` et `damp_i(y)` (activation et inhibition saturantes par fonctions
de Hill) ne dépendent que de l'état observé. `g*_i` est l'estimateur du
maximum de vraisemblance local.

## MAP avec priors de littérature

Le prior log-normal (moyenne 0.35/jour, écart-type log 0.5) provient de la
littérature immunologique. Le MAP `g_map` rétrécit `log g*` vers la moyenne du
prior :

```
log g_map = ( w_data · log g* + w_prior · log μ ) / ( w_data + w_prior )
```

Le poids de la donnée est atténué pour les populations proches de leur capacité
de charge (le niveau observé y contraint alors faiblement le taux de croissance,
et le prior doit dominer pour éviter des estimations extrêmes dues au bruit de
mesure).

## Intervalles de confiance par bootstrap

Le bilan observé est rééchantillonné en log-normal multiplicatif
(`CV` ≈ 15 %, typique de la cytométrie en flux) ; le MAP est recalculé sur
chaque réplica. Les intervalles de confiance à 95 % sont les percentiles 2.5 et
97.5 de la distribution bootstrap de chaque paramètre.

## Validation

Validé par auto-consistance (mode `validate` sur données synthétiques) : des
paramètres fixés → simulation → bruit → recalibration → récupération des
paramètres à ±25 % relatifs, intervalles contenant la vérité pour la plupart
des populations.