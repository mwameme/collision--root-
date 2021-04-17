#include "carte.h"
#include <vector>
#include <TRandom.h>
#include <TMath.h>
#include <iostream>

using namespace std;


void Bille::collision(Bille& autre) {
	//coordonnées sphériques.
	double r;
	if (this->m_type == autre.m_type)
		r = 2 * this->m_r;
	else
		r = abs(this->m_r - autre.m_r);
	double e[3] = { 0,0,0 };
	for (int i(0); i < 3; ++i) 
		e[i] = (this->m_x[i] - autre.m_x[i]) / r; //le vecteur de norme 1 reliant les centres des billes

	double m = ((this->m_v[0] - autre.m_v[0]) * e[0] + (this->m_v[1] - autre.m_v[1]) * e[1] + (this->m_v[2] - autre.m_v[2]) * e[2]); //la vitesse relative entre les centres
	double vPar[3] = { m * e[0],m * e[1],m * e[2] }; //la vitesse parallèle
	//double vOrth[3] = { (m_v[0] - autre.m_v[0]) - vPar[0],(m_v[1] - autre.m_v[1]) - vPar[1] ,(m_v[2] - autre.m_v[2]) - vPar[2] };//la vitesse orthogonale : se conserve !
	for (int i(0); i < 3; ++i) {
		this->m_v[i] -= vPar[i]*(2* autre.m_m /(this->m_m + autre.m_m));
		autre.m_v[i] += + vPar[i]*2*this->m_m/(this->m_m+autre.m_m);
	}
	return;
}


double carre(double nombre) {
	return nombre * nombre;
}

double Bille::v2(Bille& autre) {
	return carre(this->m_v[0] - autre.m_v[0]) + carre(this->m_v[1] - autre.m_v[1]) + carre(this->m_v[2] - autre.m_v[2]);
}



double Bille::calculerTempsBille(Bille& autre, double d) {
	double r2 = 0;
	if (this->m_type == autre.m_type)
		r2 = carre((m_r + autre.m_r + d) );
	else
		r2 = carre(abs(this->m_r - autre.m_r) - d);

	double x[3] = { this->m_x[0] - autre.m_x[0], this->m_x[1] - autre.m_x[1], this->m_x[2] - autre.m_x[2] };
	double v[3] = { this->m_v[0] - autre.m_v[0] , this->m_v[1] - autre.m_v[1] , this->m_v[2] - autre.m_v[2] };
	double delta = carre(x[0] * v[0] + x[1] * v[1] + x[2] * v[2]) - (carre(v[0]) + carre(v[1]) + carre(v[2])) * ( carre(x[0]) + carre(x[1]) + carre(x[2]) - (r2));
	if (delta < 0)
		return -1;
	if ((v[0] == 0) && (v[1] == 0) && (v[2] == 0))
		return -1;
	double t1 = (-(x[0] * v[0] + x[1] * v[1] + x[2] * v[2]) + sqrt(delta)) / (carre(v[0] ) + carre(v[1]) + carre(v[2]));
	double t2 = (-(x[0] * v[0] + x[1] * v[1] + x[2] * v[2]) - sqrt(delta)) / (carre(v[0]) + carre(v[1]) + carre(v[2]));
	if ((t1 < 0) && (t2 < 0))
		return -1;
	if (t1 < 0)
		return t2;
	if (t2 < 0)
		return t1;
	if (t1 < t2)
		return t1;
	return t2;
}//le plus petit temps positif ... ou négatif == impossible


void Boule::calculerTempsBoule(double d) {
	for(int i(1);i< m_M;++i)
		for (int j(0); j < i; ++j){
			double res = this->m_liste_billes[i].calculerTempsBille(this->m_liste_billes[j],d);
			this->m_liste_temps_billes[i][j] = res;
			this->m_liste_temps_billes[j][i] = res;
		}

		for (int i(0); i < m_M; ++i)
			this->m_liste_temps_boule[i] = this->m_liste_billes[i].calculerTempsBille(*this,d);
}

void Carte::calculerTempsCarte(double d) {
	for (int i(0); i < m_N; ++i) //calculer dans les boules : entre les billes et billes avec boule
		this->m_liste_boules[i].calculerTempsBoule(d);

	for (int i(1); i < m_N; ++i) //calculer entre les boules
		for (int j(0); j < i; ++j){
			double a = this->m_liste_boules[i].calculerTempsBille(this->m_liste_boules[j],d);
			this->m_liste_temps_boules[i][j] = a;
			this->m_liste_temps_boules[j][i] = a;
		}

		for (int i(0); i < m_N; ++i)
			this->m_liste_temps_carte[i] = this->m_liste_boules[i].calculerTempsBille(*this,d); //calculer boule + carte

		return;
}

void Carte::avancerTempsCarte(double temps) {
	for (int i(0); i < m_N; ++i)
		this->m_liste_boules[i].avancerTempsBoule(temps); //dans chaque boule et entre billes
	this->avancerTempsBille(temps); // temps de collison entre boules + carte
	for (int i(0); i < m_N; ++i)
		this->m_liste_temps_carte[i] -= temps;
	for (int i(1); i < m_N; ++i)
		for (int j(0); j < i; ++j) {
			this->m_liste_temps_boules[i][j] -= temps; //temps de collision entre boules
			this->m_liste_temps_boules[j][i] -= temps;
		}

	return;
}

void Boule::avancerTempsBoule(double temps) { //avance dans la boule, entre les billes, et leurs positions (boule + billes)

	for (int i(0); i < m_M; ++i)
		this->m_liste_billes[i].avancerTempsBille(temps);
	this->avancerTempsBille(temps);
	for (int i(0); i < m_M; ++i)
		for (int j(0); j < i; ++j) {
			this->m_liste_temps_billes[i][j] -= temps;
			this->m_liste_temps_billes[j][i] -= temps;
		}
	for (int i(0); i < m_M; ++i)
		this->m_liste_temps_boule[i] -= temps;

	return;
}


void Bille::avancerTempsBille(double temps) {
	for (int i(0); i < 3; ++i)
		this->m_x[i] = this->m_x[i] + temps * m_v[i];
	return;
}

double Bille::distance(Bille& autre) {
	double somme = 0;
	for (int i(0); i < 3; ++i)
		somme += carre(this->m_x[i] - autre.m_x[i]);
	return sqrt(somme);
}

void Carte::calculCollision(double d) {
	bool modifie;
	int compteur = 0;
label:

label4:
	compteur = 0;
	modifie = false;
	for (int i(0); i < this->m_N; ++i)
		for (int j(1); j < this->m_M; ++j)
			for (int k(0); k < j; ++k) //dans la boule i entre les billes j et k ...
				if (this->m_liste_boules[i].m_liste_billes[j].distance(this->m_liste_boules[i].m_liste_billes[k]) < (this->m_liste_boules[i].m_liste_billes[j].m_r * 2 + d * 2))
					if (this->m_liste_boules[i].m_liste_billes[j].calculerTempsBille(this->m_liste_boules[i].m_liste_billes[k], 0) > 0) {
						this->m_liste_boules[i].m_liste_billes[j].collision(this->m_liste_boules[i].m_liste_billes[k]); //on fait la collision
						this->m_liste_boules[i].m_liste_temps_billes[j][k] = -1; //les billes ne se retouchent plus ...
						this->m_liste_boules[i].m_liste_temps_billes[k][j] = -1;
						//puis on recalcule le temps : les deux billes avec la boule i , les deux billes, les deux billes avec les autres billes ...
						this->m_liste_boules[i].m_liste_temps_boule[j] = m_liste_boules[i].calculerTempsBille(this->m_liste_boules[i].m_liste_billes[j], d);
						this->m_liste_boules[i].m_liste_temps_boule[k] = m_liste_boules[i].calculerTempsBille(this->m_liste_boules[i].m_liste_billes[k], d);
						//m_liste_boules[i].m_liste_temps_billes[j][k] = m_liste_boules[i].m_liste_billes[j].calculerTempsBille(m_liste_boules[i].m_liste_billes[k], d);
						//m_liste_boules[i].m_liste_temps_billes[k][j] = m_liste_boules[i].m_liste_temps_billes[j][k];
						for (int l(0); l < m_M; ++l) {
							//if ((l == j) || (l == k))
							//	continue;
							if (j != l) {
								this->m_liste_boules[i].m_liste_temps_billes[j][l] = this->m_liste_boules[i].m_liste_billes[j].calculerTempsBille(this->m_liste_boules[i].m_liste_billes[l], d);
								this->m_liste_boules[i].m_liste_temps_billes[l][j] = this->m_liste_boules[i].m_liste_temps_billes[j][l];
							}
							if (k != l) {
								this->m_liste_boules[i].m_liste_temps_billes[k][l] = this->m_liste_boules[i].m_liste_billes[k].calculerTempsBille(this->m_liste_boules[i].m_liste_billes[l], d);
								this->m_liste_boules[i].m_liste_temps_billes[l][k] = this->m_liste_boules[i].m_liste_temps_billes[k][l];
							}

						}
						//goto label4;
						cout << "label4" << endl;
						//modifie = true;
						compteur++;
					}
	cout << "compteur4 : " << compteur << endl;
	if (modifie)
		goto label;

label3:
	modifie = false;
	compteur = 0;

	//dans la boule i, avec la bille j ...
	for (int i(0); i < this->m_N; ++i)
		for (int j(0); j < this->m_M; ++j)
			if (this->m_liste_boules[i].distance(this->m_liste_boules[i].m_liste_billes[j]) > (m_liste_boules[i].m_r - m_liste_boules[i].m_liste_billes[j].m_r - d * 2))
				if (this->m_liste_boules[i].calculerTempsBille(this->m_liste_boules[i].m_liste_billes[j], 0) *
					sqrt(this->m_liste_boules[i].m_liste_billes[j].v2(this->m_liste_boules[i])) < 3 * d) { //si il y a collision ... entre boule i et bille j
					this->m_liste_boules[i].collision(this->m_liste_boules[i].m_liste_billes[j]); //calcul de la collision
					//this->m_liste_boules[i].m_liste_temps_boule[j] = this->m_liste_boules[i].calculerTempsBille(this->m_liste_boules[i].m_liste_billes[j], d); //on recalculle la collision suivante ...
							//on met à jour la boule i avec : la carte, les autres boules, + (la boule et les billes, les billes). 
					this->m_liste_temps_carte[i] = this->calculerTempsBille(this->m_liste_boules[i], d);
					for (int k(0); k < this->m_N; ++k) { //boule i avec boule k ...
						if (k == i)
							continue;
						this->m_liste_temps_boules[i][k] = this->m_liste_boules[i].calculerTempsBille(this->m_liste_boules[k], d);
						this->m_liste_temps_boules[k][i] = this->m_liste_temps_boules[i][k];
					}
					for (int k(0); k < this->m_M; ++k) {//la boule i avec ses billes
						//if (k == j)
						//	continue;
						this->m_liste_boules[i].m_liste_temps_boule[k] = this->m_liste_boules[i].calculerTempsBille(this->m_liste_boules[i].m_liste_billes[k], d);
					}
					//les billes k avec la bille j
					for (int k(0); k < m_M; ++k) {
						if (k == j)
							continue;
						this->m_liste_boules[i].m_liste_temps_billes[j][k] = this->m_liste_boules[i].m_liste_billes[j].calculerTempsBille(this->m_liste_boules[i].m_liste_billes[k], d);
						this->m_liste_boules[i].m_liste_temps_billes[k][j] = this->m_liste_boules[i].m_liste_temps_billes[j][k];
					}
					//goto label3;
					//cout << "label3" << endl;
					//modifie = true;
					compteur++;
				}
	cout << "compteur3 : " << compteur << endl;
	if (modifie)
		goto label;




label2:
	modifie = false;
	compteur = 0;
	for (int i(1); i < this->m_N; ++i) //entre boule i et boule j :
		for (int j(0); j < i; ++j)
			if (this->m_liste_boules[i].distance(this->m_liste_boules[j]) < (this->m_liste_boules[i].m_r * 2 + d * 2)) //si  dans la zone de collision
				if (this->m_liste_boules[i].calculerTempsBille(this->m_liste_boules[j], 0) > 0) {
					this->m_liste_boules[i].collision(this->m_liste_boules[j]); //collision entre i et j;
					this->m_liste_temps_boules[i][j] = -1; //normalement elles ne se retouchent plus ...
					this->m_liste_temps_boules[j][i] = -1;

					//Puis recalculer temps pour : boule k avec i et k avec j, et i avec ses billes, et j avec ses billes, et i et j avec la carte...
					for (int k(0); k < m_N; ++k) {
						//if ((i == k) || (j == k))
						//	continue;
						if (k != i) {
							this->m_liste_temps_boules[i][k] = this->m_liste_boules[i].calculerTempsBille(this->m_liste_boules[k], d);
							this->m_liste_temps_boules[k][i] = this->m_liste_temps_boules[i][k];
						}
						if (k != j) {
							this->m_liste_temps_boules[j][k] = this->m_liste_boules[j].calculerTempsBille(this->m_liste_boules[k], d);
							this->m_liste_temps_boules[k][j] = this->m_liste_temps_boules[j][k];
						}
					}
					this->m_liste_temps_carte[i] = this->calculerTempsBille(this->m_liste_boules[i], d);
					this->m_liste_temps_carte[j] = this->calculerTempsBille(this->m_liste_boules[j], d);
					for (int k(0); k < m_M; ++k) {
						this->m_liste_boules[i].m_liste_temps_boule[k] = this->m_liste_boules[i].m_liste_billes[k].calculerTempsBille(this->m_liste_boules[i], d);
					}
					for (int k(0); k < m_M; ++k) {
						this->m_liste_boules[j].m_liste_temps_boule[k] = this->m_liste_boules[j].m_liste_billes[k].calculerTempsBille(this->m_liste_boules[j], d);
					}
					//m_liste_temps_boules[i][j] = m_liste_boules[i].calculerTempsBille(m_liste_boules[j], d);
					//m_liste_temps_boules[j][i] = m_liste_temps_boules[i][j];
					//goto label2;
					//cout << "label2" << endl;
					//modifie = true;
					compteur++;
				}
	cout << "compteur2 : " << compteur << endl;
	if (modifie)
		goto label;


label1:
	modifie = false;
	compteur = 0;
	//on commence par la collision de la carte avec une boule
	for (int i(0); i < this->m_N; ++i)
		if (this->distance(this->m_liste_boules[i]) > (this->m_r - this->m_liste_boules[i].m_r) - d * 2) //si dans la distance de collision (d*2) : carte + boule
			if (this->calculerTempsBille(m_liste_boules[i], 0) *
				sqrt(this->v2(this->m_liste_boules[i])) < 3 * d) //si il y aura une collision !
			{
				this->collision(this->m_liste_boules[i]); //on change les vitesses, puis on recalcule les temps
				for (int j(0); j < m_N; ++j) { //on calcule le temps pour la boule i avec les autres boules (j)
					if (j == i)
						continue;
					this->m_liste_temps_boules[i][j] = this->m_liste_boules[i].calculerTempsBille(this->m_liste_boules[j], d);
					this->m_liste_temps_boules[j][i] = this->m_liste_temps_boules[i][j];
				}
				for (int j(0); j < m_N; ++j) //dont la carte avec la boule j ...
					this->m_liste_temps_carte[j] = this->calculerTempsBille(this->m_liste_boules[j], d);//on calcule le temps pour la carte avec les autres boules
				for (int j(0); j < m_M; ++j) //on calcule la boule i avec ses billes (j)
					this->m_liste_boules[i].m_liste_temps_boule[j] = this->m_liste_boules[i].m_liste_billes[j].calculerTempsBille(this->m_liste_boules[i], d); //entre la boule et ses billes. On ne calcule pas entre les billes !
				//fin du recalcul de temps ! pour la collsion de la carte avec la boule i.
				//goto label1;
				//cout << "label1" << endl;
				//modifie = true;
				compteur++;
			}
	cout << "compteur1 : " << compteur << endl;
	if (modifie)
		goto label;


	return;
}

double Carte::tempsMin() {
	double t = 0;
	for (int i(0); i < this->m_N; ++i)
		if (this->m_liste_temps_carte[i] > 0)
			t = this->m_liste_temps_carte[i]; //on prend une valeur de t >0, quelconque

	for (int i(0); i < this->m_N; ++i)
		if ((this->m_liste_temps_carte[i] < t) && (this->m_liste_temps_carte[i]>0))
			t = this->m_liste_temps_carte[i];

	for (int i(1); i < this->m_N; ++i)
		for (int j(0); j < i; ++j)
			if ((this->m_liste_temps_boules[i][j] < t) && (this->m_liste_temps_boules[i][j]>0))
				t = this->m_liste_temps_boules[i][j];

	for(int i(0);i< this->m_N;++i)
		for(int j(0);j< this->m_M;++j)
			if ((this->m_liste_boules[i].m_liste_temps_boule[j] <t) && (this->m_liste_boules[i].m_liste_temps_boule[j] >0))
				t = this->m_liste_boules[i].m_liste_temps_boule[j];

	for(int i(0);i< this->m_N;++i)
		for(int j(1);j< this->m_M;++j)
			for(int k(0);k<j;++k)
				if ((this->m_liste_boules[i].m_liste_temps_billes[j][k] >0) && (this->m_liste_boules[i].m_liste_temps_billes[j][k] <t))
					t = this->m_liste_boules[i].m_liste_temps_billes[j][k];

	return t; //le plus petit temps positif
}


void Carte::iterer(double tempsMin, double tempsMax,double d) {
	double t = 0;
	this->calculerTempsCarte(d);

	while(t < tempsMax) {
		double deltaT = this->tempsMin();
		this->avancerTempsCarte(deltaT);
		if (t > tempsMin)
			remplirHisto(deltaT);
		this->calculCollision(d);
		t += deltaT;
		cout << "temps cumule : " << t << endl;
	}
}

void Carte::remplirHisto(double t) {
	for (int i(0); i < m_N; ++i) {
		double vitesse[3] = { 0,0,0 };
		for (int s(0); s < 3; ++s)
			vitesse[s] += m_liste_boules[i].m_v[s] * masse_Boule;
		for (int j(0); j < m_M; ++j)
			for (int s(0); s < 3; ++s)
				vitesse[s] += m_liste_boules[i].m_liste_billes[j].m_v[s] * masse_Bille;
		for (int s(0); s < 3; ++s)
			vitesse[s] = vitesse[s] / (masse_Boule + m_M * masse_Bille);
		m_histo_Boules->Fill(vitesse[0], t);
		m_histo_Boules->Fill(vitesse[1], t);
		m_histo_Boules->Fill(vitesse[2], t);
		//cout << "vitesse : " << vitesse[0] << endl;
		for (int j(0); j < m_M; ++j) {
			m_histo_Billes->Fill(m_liste_boules[i].m_liste_billes[j].m_v[0] - vitesse[0],t);
			m_histo_Billes->Fill(m_liste_boules[i].m_liste_billes[j].m_v[1] - vitesse[1],t);
			m_histo_Billes->Fill(m_liste_boules[i].m_liste_billes[j].m_v[2] - vitesse[2],t);
		}
	}
}

Carte::Carte() {

};

Carte::Carte(int n, int m, double temperature, double fraction, double rapport_masse,double rayon) { //n boules, n*n billes (n=15), fraction de l'espace (r), rapport des masses (=3). Ne pas oublier le type.
	this->m_N = n;
	this->m_M = m;
	this->masse_Carte = 1;
	this->masse_Boule = 1 / rapport_masse;
	this->masse_Bille = 1 / carre(rapport_masse);

	for (int j(0); j < this->m_N; ++j)
		this->m_liste_temps_boules.push_back(vector<double>(m_N));

	cout << "allocation liste temps boules" << endl;


	this->m_liste_temps_carte.resize(m_N);


	//TRandom essai();
	m_liste_boules.reserve(sizeof(vector<int>) + m_N * (sizeof(int) + sizeof(Bille) + m_M * sizeof(Bille) + sizeof(vector<int>) * 3 + m_M * (sizeof(vector<int>) + m_M * sizeof(double)) + m_M * sizeof(double)));
	cout << "resize liste boules" << endl;

	double a = rayon * 2 / sqrt(3);
	int nCoord = ceil(pow(m_N, 1. / 3.));
	double rayon_boule = a / (nCoord * fraction);

	int k[3] = { 0,0,0 };

	for (int i(0); i < m_N; ++i) {
		Boule boule(temperature, fraction, rapport_masse, rayon_boule,&essai,m);
		boule.m_x[0] = boule.m_x[0] - a / 2 + a / nCoord * (.5 + k[0]);
		boule.m_x[1] = boule.m_x[1] - a / 2 + a / nCoord * (.5 + k[1]);
		boule.m_x[2] = boule.m_x[2] - a / 2 + a / nCoord * (.5 + k[2]);
		for (int j(0); j < m_M; ++j) {
			boule.m_liste_billes[j].m_x[0] += -a / 2 + a / nCoord * (.5 + k[0]);
			boule.m_liste_billes[j].m_x[1] += -a / 2 + a / nCoord * (.5 + k[1]);
			boule.m_liste_billes[j].m_x[2] += -a / 2 + a / nCoord * (.5 + k[2]);
		}
		k[0] += 1;
		if (k[0] == nCoord) {
			k[0] = 0;
			k[1] += 1;
			if (k[1] == nCoord) {
				k[1] = 0;
				k[2] += 1;
			}
		}
		this->m_liste_boules.push_back(boule);
		cout << "ajouter Boule " << i << endl;
	}

	m_histo_Billes = new TH1D("billes", "energie interne", 150,  -3 *temperature, 3 * temperature);
	m_histo_Boules = new TH1D("boules", "energie des boules", 150, -3 * temperature, 3 * temperature);

	cout << "allouer histogrammes" << endl;


	this->masse_Carte = 1;
	this->masse_Boule = 1 / rapport_masse;
	this->masse_Bille = 1 / carre(rapport_masse);

	double m_x[3] = {0.,0.,0.};
	this->m_m = 1;
	this->m_type = 0;
	this->m_r = rayon;

	double impulsion[3] = { 0,0,0 };
	cout << "avant calcul vitesse" << endl;
	cout << "nombre de boules" << this->m_liste_boules.size() << endl;

	for (int i(0); i < this->m_N; ++i) {
		impulsion[0] += m_liste_boules[i].m_v[0] * (this->masse_Boule);
		impulsion[1] += m_liste_boules[i].m_v[1] * (this->masse_Boule);
		impulsion[2] += m_liste_boules[i].m_v[2] * (this->masse_Boule);
		cout << "vitesse boule " <<i<< endl;
		cout << "taille liste " << m_liste_boules[i].m_liste_billes.size() << endl;
		for (int j(0); j < this->m_M; ++j) {
			impulsion[0] += m_liste_boules[i].m_liste_billes[j].m_v[0] * (this->masse_Bille);
			impulsion[1] += m_liste_boules[i].m_liste_billes[j].m_v[1] * (this->masse_Bille);
			impulsion[2] += m_liste_boules[i].m_liste_billes[j].m_v[2] * (this->masse_Bille);
		}
	}
	cout << "calcul vitesse fin" << endl;

	for (int k(0); k < 3; ++k) {
		this->m_v[k] = - impulsion[k];
		this->m_x[k] = 0;
	}



}

Boule::Boule() {

}

Boule::Boule(double temperature, double fraction, double rapport_masse, double rayon_boule,TRandom* essai,int m) {
	double vitesse[3];
	vitesse[0] = essai->Gaus(0, temperature);
	vitesse[1] = essai->Gaus(0, temperature);
	vitesse[2] = essai->Gaus(0, temperature);

	this->m_M = m;
	this->m_liste_billes.resize(0);


	int nCoord = ceil(pow(m_M, 1. / 3.));
	double a = 2 * rayon_boule / sqrt(3);

	int k[3] = { 0,0,0 };
	cout << "allocation boule ..." << endl;

	for (int i(0); i < this->m_M; ++i) {
		Bille bille;
		bille.m_x[0] = a / nCoord * (.5 + k[0]) - a / 2;
		bille.m_x[1] = a / nCoord * (.5 + k[1]) - a / 2;
		bille.m_x[2] = a / nCoord * (.5 + k[2]) - a / 2;

		bille.m_v[0] = essai->Gaus(vitesse[0], temperature);
		bille.m_v[1] = essai->Gaus(vitesse[1], temperature);
		bille.m_v[2] = essai->Gaus(vitesse[2], temperature);

		k[0] += 1;
		if (k[0] == nCoord) {
			k[0] = 0;
			k[1] += 1;
			if (k[1] == nCoord) {
				k[1] = 0;
				k[2] += 1;
			}
		}

		bille.m_r = (a / nCoord) * sqrt(3) / 2 / fraction;
		bille.m_type = 2;
		bille.m_m = 1 / carre(rapport_masse);
		this->m_liste_billes.push_back(bille);
		cout << "bille : " << i << endl;
	}

	double d_vitesse[3] = { 0,0,0 };
	for (int i(0); i < this->m_M; ++i) {
		d_vitesse[0] += m_liste_billes[i].m_v[0] - vitesse[0];
		d_vitesse[1] += m_liste_billes[i].m_v[1] - vitesse[1];
		d_vitesse[2] += m_liste_billes[i].m_v[2] - vitesse[2];
	}
	for (int k(0); k < 3; ++k) {

		this->m_x[k] = 0;
		this->m_v[k] = -d_vitesse[0] / rapport_masse;
	}
	this->m_r = rayon_boule;
	this->m_m = 1 / rapport_masse;
	this->m_type = 1;

	for(int j(0);j< this->m_M;++j)
		this->m_liste_temps_billes.push_back(vector<double> (this->m_M));

	this->m_liste_temps_boule.resize(this->m_M);
}

Bille::Bille() {

}

Carte::~Carte() {
	delete m_histo_Billes;
	delete m_histo_Boules;

}