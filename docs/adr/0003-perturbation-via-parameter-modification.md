# Encodage des perturbations par modification paramétrique

Les immunothérapies modifient les paramètres cinétiques du modèle plutôt que d'ajouter des termes de forçage externe ou de modifier la topologie du graphe.

## Contexte

Il faut un moyen de représenter l'effet d'un traitement (anti-PD1, anti-TNF, corticoïdes) dans les ODEs. Le choix a un impact sur l'interprétabilité clinique et la maintenance du modèle.

## Considérées

- **Termes de forçage externe** : `dx/dt = f(x) + u(t)`. Simple mais difficile à interpréter cliniquement.
- **Modification de topologie** : un traitement "découpe" ou "ajoute" un arc du graphe. Puissant mais complexe à maintenir.
- **Modification paramétrique** : un traitement modifie la valeur d'un ou plusieurs paramètres cinétiques (ex: anti-PD1 → β_CDL8↑ de 40%). Naturel et interprétable.

## Décision

Modification paramétrique. Chaque perturbation est définie par un mapping : `traitement → {paramètre: delta_relatif}`. Les deltas sont documentés et calibrés sur des données de littérature. Exemple : anti-PD1 augmente le taux d'activation des CD8+ de 30-50% selon les études.
