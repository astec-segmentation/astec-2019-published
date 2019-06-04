********************************************************************************
***************************** Format for Morpheme_lineage **********************
********************************************************************************

I - INTRODUCTION
----------------

Fichier texte (format .txt) permettant la sauvegarde de correspondances d'objets
(cliques) et plus spécifiquement de cellules segmentées provenant de séquences 
d'images ND+t. 


À chaque instant de la séquence, une cellule est décrite par une clé selon des 
règles de nomenclature spécifiées dans ce document.

Le format proposé permet l'écriture de correspondances de groupes de cellules
vers des groupes de cellules (cliques) et peut s'étendre au delà de 
l'établissement de correspondances intra-séquence.

Les lignes débutant par le caractère '#' sont ignorées (commentaires).

II - CLÉ D'UN OBJET
-------------------

La clé d'un objet débute par le caractère '[', se termine par le caractère ']',
et est constituée d'une juxtaposition de descripteurs de la forme :
<CODE_DESCRIPTEUR> <valeur_descripteur>
avec un caractère espace séparant le code de la valeur du descripteur.

Les descripteurs d'une clé sont séparés d'une virgule ','.

Ainsi, la clé associée à un objet suit la syntaxe suivante :
[<CODE_1> <valeur_1>, ... , <CODE_N> <valeur_N>]

Note 1 : le code associé à un descripteur ne doit contenir aucun espace.
Note 2 : tous les objets contenus dans un même fichier de lignage doivent être 
identifiés à l'aide du même ensemble de descripteurs.

III - GROUPE D'OBJETS
---------------------


Un groupe d'objets débute par le caractère '{', se termine par le caractère '}' 
et contient la juxtaposition de la liste des clés associées aux objets qui 
composent le groupe où les éléments sont séparés par une virgule ',' :
{} : groupe vide
{[<C_1> <v_1>, ... , <C_N> <v_N>]} : singleton
{[<C_1> <v_1>, <C_2> <v_2>], [<C_1> <w_1>, <C_2> <w_2>]} : groupe de 2 objets
{<clé_1>, ..., <clé_N>} : groupe de N objets  
  
  
IV - CLIQUE
-----------

Une clique décrit une correspondance entre deux groupes d'objets.
Le format impose :
 - une clique par ligne
 - deux groupes par clique
 - le séparateur entre les deux groupes de la clique est le caractère '-'

{[<C> <v>]} - {[<C> <w>]} : correspondance entre deux singletons

V - CAS PARTICULIER DES LIGNÉES DE CELLULES INTRA-SÉQUENCES
-----------------------------------------------------------

  V.1 - Clé d'une cellule
  -----------------------
  Une séquence temporelle se compose de N images segmentées identifiées selon 
  l'instant de l'acquisition (pour les embryons d'ascidies, une image par 
  minute).

  L'instant d'acquisition compte parmi les descripteurs avec le code 'T'.
  Ainsi, une cellule (c'est-à-dire un objet) d'une image prise à l'instant 192
  aura une clé contenant le descripteur 'T 192'.

  Pour une image donnée, une étiquette (label) est associée à chaque cellule.
  Le label compte parmi les descripteurs avec le code 'L'.
  Ainsi, une cellule étiquetée 84 dans une image aura une clé contenant le
  descripteur 'L 84'.

  De façon optionnelle, il est aussi possible d'enrichir l'information relative 
  à une cellule en ajoutant un descripteur permettant d'identifier l'embryon 
  auquel appartient la cellule (c'est-à-dire permettant de savoir à quelle 
  séquence la cellule se réfère). Ce descripteur est associé au code 'E'.
  Ainsi, une cellule appartenant à l'embryon nommé '161213_lucie' aura une clé
  contenant le descripteur 'E 161213_lucie'.

  Exemple : cellule décrite par temps=192, label=84, embryon='161213_lucie'
            clé associée : [T 192, L 84, E 161213_lucie]
 
  V.2 - Groupe de cellules
  ------------------------
  Dans une même séquence, on construit des groupes de cellules qui possèdent 
  chacune le même temps :
  {[T <x>, L <y>], [T <w>, L <z>]} : non-accepté si <x> différent de <w>
  {[T <x>, L <y>], [T <x>, L <z>]} : accepté (<x> est alors le temps associé au
                                     groupe)
  
  
  V.3 - Clique de cellules
  ------------------------
  Une clique de cellules met en correspondance deux groupes de cellules G et H
  ayant respectivement g et h comme temps associés avec la contrainte 
  h = next(g) où next décrit l'instant "suivant" de la séquence temporelle. 
  Si on considère une séquence à la résolution d'une image par minute, alors 
  la contrainte devient h = g + 1.
  
  Une cellule n'ayant pas de correspondance s'associe au groupe vide.
  
  Exemple :
  {[T 8, L 55]} - {[T 9, L 55]}
  {[T 9, L 55]} - {[T 10, L 68],[T 10, L 69]}
  {[T 10, L 68]} - {[T 11, L 68]}
  {[T 10, L 69]} - {[T 11, L 69]}

  {[T 191, L 1441]} - {[T 192, L 1458],[T 192, L 1459]}
  {[T 192, L 1458]} - {}
  {[T 192, L 1459]} - {}
  
  
