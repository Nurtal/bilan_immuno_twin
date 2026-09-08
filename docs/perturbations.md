# Perturbation catalogue (immunothérapies)

Les perturbations sont encodées par **modification de paramètres cinétiques**
(ADR-0003). Chaque traitement est un mapping `{identifiant_de_paramètre:
delta_relatif}` où un delta de `0.4` = multiplication par `1.4` (+40 %) et un
delta de `-0.6` = multiplication par `0.4` (−60 %).

## Identifiants de paramètres

| Identifiant | Paramètre |
| --- | --- |
| `growth:<population>` | taux de prolifération de base (par jour) |
| `death:<population>` | taux de mortalité de base (par jour) |
| `weight:<source>-><cible>:act` | poids d'un arc d'activation du graphe |
| `weight:<source>-><cible>:sup` | poids d'un arc d'inhibition du graphe |
| `capacity:<population>` | capacité de charge d'une population |
| `basal:<population>` | activation basale d'une population |

## Catalogue

### anti-PD1 — blocage de PD-1

Lève l'épuisement des lymphocytes T CD8+ : augmentation de la prolifération
CD8+, réduction de leur mortalité, et amplification de l'activation CD8+
dépendante de l'IFN-γ.

| Paramètre | Delta |
| --- | --- |
| `growth:CD8` | +40 % |
| `death:CD8` | −10 % |
| `weight:IFNg->CD8:act` | +30 % |

### anti-TNF — blocage du TNF-α

Réduit le signal pro-inflammatoire TNF-α qui active les monocytes/macrophages.

| Paramètre | Delta |
| --- | --- |
| `weight:TNFa->Monocytes:act` | −60 % |

### corticoïde — corticoïdes

Déprime plusieurs signaux pro-inflammatoires : croissance Th1, Th17 et
monocytes, ainsi que la signalisation TNF-α/IL-6.

| Paramètre | Delta |
| --- | --- |
| `growth:Th1` | −35 % |
| `growth:Th17` | −40 % |
| `growth:Monocytes` | −30 % |
| `weight:TNFa->Monocytes:act` | −50 % |
| `weight:IL6->Th17:act` | −40 % |
| `weight:IL6->Monocytes:act` | −30 % |