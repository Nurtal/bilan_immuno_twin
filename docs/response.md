# Score de réponse

Le score de réponse est **continu**, délibérément différent d'un classifieur
binaire, pour éviter une fausse précision clinique.

## Définition

Pour un traitement donné, chaque population possède une **direction cible**
(+1 = une réponse favorable accroît la population, −1 = une réponse favorable la
réduit, 0 = neutre) :

| Traitement | CD8 | Th1 | Th2 | Th17 | B | NK | Treg | Monocytes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| anti-PD1 | +1 | +1 | 0 | 0 | 0 | +1 | 0 | 0 |
| anti-TNF | 0 | −1 | 0 | −1 | 0 | 0 | 0 | −1 |
| corticoïde | 0 | −1 | 0 | −1 | 0 | 0 | 0 | −1 |

Le score est la somme pondérée des log-2 des changements de pli (fold-change)
de fin de simulation, perturbé vs non perturbé :

```
score = Σ_i direction_i · log2(fold_i)
```

- `score > 0` : réponse **favorable**
- `score ≈ 0` : réponse **neutre** (pas de réponse nette)
- `score < 0` : réponse **défavorable**

## Interprétation

Le score est émis tel quel dans la sortie CLI (pleine précision continue) avec
une bande qualitative (favorable / neutre / défavorable) donnée comme aide de
lecture — ce n'est **pas** une classification binaire.

## Réponse différentielle

Un patient est marqué **differential** quand son score s'écarte du score de
référence — obtenu sur un bilan immunologique typique avec le même traitement
et le même horizon — de plus de `DIFFERENTIAL_SCORE_TOL` (0,5). Les populations
**unexpected** sont celles dont le fold contredit la direction cible du
traitement ou sort de la plage attendue (`[0.5, 2.0]`).