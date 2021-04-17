# collision--root-

dans une carte, des boules s'entrechoquent de maniere élastique : sans perte d'énergie.
Chaque boule est considérée creuse, et dedans des billes s'entrechoquent de maniere élastique.

On calcule les temps de collision. On calcule pour une distance petite (delta)
On prend t_min et on avance tous les temps (et on fait avancer les billes+boules).

on regarde les objets en situation de collision : si la distance est plus petite que 2 delta.

Si en plus elles VONT se collisionner (si temps de collision >0 pour billes et billes ou pour boules et boules)
(et si distance de collison < 3 *delta pour billes et boules)

alors on calcule la collision : changement de vitesse.
Et on recalcule les temps de collisions pour les objets modifiés.
Par exemple si une boule change, les temps de collisions avec ses billes et avec les autres boules changent.

les temps de collision sont calculés : (x+t*v)^2 = r^2 avec x la différence de position entre les objets. Idem pour v.
les changements de vitesse : on regarde le vecteur e de norme 1 dirigé d'un centre à l'autre.
La partie orthogonale à e ne change pas. La partie parallèle change :
https://fr.wikipedia.org/wiki/Choc_%C3%A9lastique#Cas_en_une_dimension

On remplit un histogramme (télécharger et installer root.cern) : 
on regarde les vitesses des billes à l'intérieur de chaque boule.
On ajoute v_bille - <v_boule> : l'écart de vitesse entre la bille et la vitesse moyenne. C'est l'énergie interne de la boule.
On ajoute la vitesse moyenne dans la boule : <v_boule>.

Cela donne deux histogrammes. Il faut explorer la situation en changeant les parametres. Ceux-ci sont :

Carte(int n, int m, double temperature, double fraction, double rapport_masse, double rayon);
n : nombre de boules
m : nombre de billes par boule
temperature : de la gaussienne des vitesses initiales
fraction : dans la carte, on prend le cube inscrit, qu'on découpe en par exemple en 3*3*3. Donc le rayon de la boule c'est le rayon du cube / (3 * fraction). Fraction >1 !
rapport_masse : la masse de la carte c'est 1, la masse des boules c'est 1/rapport_masse, la masse des billes c'est 1/rapport_masse ^2.
Le rayon c'est le rayon de la carte ...

Apres je calcule pour t variant de (par exemple) 5 à 200 (de 0 à 5 je n'enregistre pas les données).
et aussi pour les calculs, je dis qu'il y a une distance de collision (10^-5 ici) d'écart entre les boules/billes.
Parceque l'ordinateur gere mal les calculs réels = 0. C'est comme si d remplaçait 0 ...

notez les types : 0,1,2 : carte boule bille.
Les liste de temps : carte avec boules et boules avec boules : listes dans la carte.
boule avec billes et billes avec billes : listes dans la boule.

les conditions initiales changent peu le résultat. Seules comptent l'énergie moyenne, et le temps (court) au bout duquel le systeme "thermalise" : la répartition en gaussienne thermique.

Enjoy reading !
