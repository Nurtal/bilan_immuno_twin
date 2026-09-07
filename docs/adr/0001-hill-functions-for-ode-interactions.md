# Fonction de Hill comme formalisme des interactions ODE

Le modèle utilise des fonctions de Hill (`f(x) = x^n / (K^n + x^n)`) pour représenter les interactions régulatrices dans les ODEs plutôt que des cinétiques mass-action bilinéaires.

## Contexte

Les interactions immunitaires présentent des effets de seuil et de saturation biologiques : un signal doit atteindre un certain niveau pour activer une réponse, et au-delà d'un certain point, l'effet est plafonné. La mass-action simple (`β·x·y`) modélise mal ces phénomènes.

## Considérées

- **Mass-action simple** : `dx/dt = β·x·y - δ·x`. Plus simple, mais ne capture ni seuils ni saturations.
- **Michaelis-Menten** : standard en biochimie enzymatique, mais moins adapté aux régulations inter-populations.
- **Fonctions de Hill** : `f(x) = x^n / (K^n + x^n)`. Standard en modélisation des systèmes biologiques, paramètres interprétables (K = seuil, n = coopérativité).

## Décision

Hill functions. Le paramètre n (coopérativité) est fixé à des valeurs raisonnables (généralement 1-4) et K est estimé par calibration patient-specific.
