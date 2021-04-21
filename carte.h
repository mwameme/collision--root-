#pragma once
#include <TApplication.h>

#include <vector>
#include <TH1D.h>
#include <TRandom.h>
#include <TMath.h>



class Bille {
public:
	double m_r;
	double m_x[3];
	double m_v[3];
	double m_m;
	int m_type; //0 si carte, 1 si boule, 2 si bille

	void collision(Bille& autre); //dépend du type
	double calculerTempsBille(Bille& autre,double d); //calcule le temps de collision entre deux billes
		//calcule le temps de collision à r(1+0,0001),
		//et inverse les vitesses en dessous de r(1+0,0002), tant que collision existe à r*1;
		//e vaut (1+ ...)
		
	void avancerTempsBille(double temps); //x + v*t
	double distance(Bille& autre);
	double v2(Bille& autre);
	bool enCollision(Bille& autre,double d);
	Bille(); //générer la température gaussienne, la  position quelconque
};

class Boule : public Bille {
public:
	std::vector<Bille> m_liste_billes;
	std::vector<std::vector<double>> m_liste_temps_billes; //collision entre billes et la boule
	std::vector<double> m_liste_temps_boule;

	//std::vector< std::vector<bool>> m_collision; //vrai si dR < r(1+0,0001) et collision existe (r*1)
	//double m_vitesse[3]; //vitesse moyenne dans la boule
	int m_M; //nombre de billes par boule

	void calculerTempsBoule(double d); //calcule le temps entre les billes Remplit le tablea
				//d = 0 + 10^-4
	void avancerTempsBoule(double temps); //avance pour toutes les billes, et la boule
	Boule(double temperature, double fraction, double rapport_masse, double rayon_boule,TRandom* essai,int m);//remplir n billes
	Boule();
};

class Carte : public Bille{ // d'impulsion totale nulle. De masse a  choisir.
public:
	std::vector<Boule> m_liste_boules;
	std::vector<std::vector<double>> m_liste_temps_boules; //collision entre boules
	std::vector<double> m_liste_temps_carte; //collision avec le contour ... qui est une bille !

	TH1D* m_histo_Billes;
	TH1D* m_histo_Boules;

	double masse_Carte;
	double masse_Boule;
	double masse_Bille;
	int m_N; //nombre de boules par carte
	int m_M;

	TRandom essai;

	void calculerTempsCarte(double d); //remplit le tableau Temps
	void avancerTempsCarte(double temps); //avance pour toutes les boules, et pour la carte
	void calculCollision(double d);

	void iterer(double tempsMin, double tempsMax,double d);
	double tempsMin();
	void remplirHisto(double t);

	double energie();

	Carte(int n,int m,double temperature,double fraction, double rapport_masse,double rayon); //remplir n boules de n billes, avec une température ! fraction du volume occupé
	Carte();
	~Carte();
};

