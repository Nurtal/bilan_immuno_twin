# CLI Python comme seule interface d'interaction

Le jumeau s'interface exclusivement via une CLI Python, sans notebooks Jupyter ni API REST.

## Contexte

La majorité des outils de recherche en bioinformatique utilisent des notebooks Jupyter pour l'exploration. Cependant, pour un jumeau clinique destiné à la reproductibilité et potentiellement à l'intégration dans un pipeline hospitalier, la CLI offre des avantages clés.

## Considérées

- **Notebooks Jupyter** : excellent pour l'exploration interactive, mais mauvaise reproductibilité, difficultés de versioning, output non-structuré.
- **API REST** : intégration facile dans des dashboards cliniques, mais over-engineering pour un MVP, nécessite un serveur.
- **CLI Python** : scriptable, reproductible, output structuré (JSON/CSV), compatible avec des pipelines CI/CD, pas de dépendance serveur.

## Décision

CLI Python uniquement. Les résultats sont exportés en JSON/CSV pour analyse ultérieure. Si un dashboard clinique est requis plus tard, l'API REST sera ajoutée comme couche au-dessus de la CLI.
