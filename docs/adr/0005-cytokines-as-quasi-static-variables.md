# Cytokines comme variables quasi-statiques

Les cytokines (IFN-γ, IL-10, TNF-α, IL-6, IL-4, IL-17…) sont modélisées comme des variables quasi-statiques calculées algébriquement à partir des populations cellulaires, plutôt que comme des variables d'état avec leurs propres ODEs.

## Contexte

Les cytokines évoluent sur des échelles de temps très differentes des populations cellulaires : heures pour les cytokines vs jours/semaines pour les cellules. Sur l'horizon de simulation typique (2-4 semaines), les cytokines atteignent quasi-instantanément leur équilibre par rapport aux populations.

## Considérées

- **Variables d'état (ODEs séparées)** : plus réaliste, mais double le nombre d'EDOs, nécessite des paramètres cinétiques cytokine supplémentaires difficiles à calibrer.
- **Paramètres fixes** : trop rigide, ne permet pas de capturer les feedbacks cytokine-cellule.
- **Variables quasi-statiques** : calculées algébriquement à chaque pas de temps. `[IFN-γ] = f(Th1, NK, CD8+)`. Réaliste, simple, interpretable.

## Décision

Quasi-statiques. Les équations algébriques de cytokines sont définies dans le graphe de régulations et documentées. Chaque cytokine est une fonction linéaire ou de Hill des populations qui la produisent.
