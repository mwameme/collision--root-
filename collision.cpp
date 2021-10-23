#include <TApplication.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH1.h>
#include <TF1.h>
#include <thread>
#include <chrono>
#include "carte.h"

#include <iostream>
using namespace std;


//TCanvas canvas ("fCanvas", "fCanvas", 600, 400);

int main(int argc, char* argv[])
{
    double fitGaussBoule[20];
    double fitGaussBille[20];
    double tabEnergie[20];

    TApplication app("app", &argc, argv);

    /*
    for (int i = 9; i <= 14; ++i) {


        //TCanvas* canvas = new TCanvas("fCanvas", "fCanvas", 600, 400);
        //c1.Divide(2, 1);

        //TH1D h("h", "h", 10, 0, 10);
        //h.Fill(1);
        //h.Draw();

        //Carte(int n, int m, double temperature, double fraction, double rapport_masse, double rayon);
        cout << "carte numero " << i << endl;
        Carte carte(20, 5*i, 0.5, 5., 1., 1);
        carte.iterer(2, 40 , 0.000001);

        TCanvas c1("c1", "Candle Decay", 800, 600);
        c1.cd(1);
        carte.m_histo_Billes->Draw();
        c1.cd(2);
        carte.m_histo_Boules->Draw();

        c1.Update();
        c1.Draw();

        TF1* f1 = new TF1("f1", "gaus", -.5, .5);
        // set initial parameters (not really needed for gaus)
        f1->SetParameters(carte.m_histo_Boules->GetMaximum(), carte.m_histo_Boules->GetMean(), carte.m_histo_Boules->GetRMS());
        carte.m_histo_Boules->Fit("f1");
        fitGaussBoule[i - 1] = f1->GetParameter(2);

        f1->SetParameters(carte.m_histo_Billes->GetMaximum(), carte.m_histo_Billes->GetMean(), carte.m_histo_Billes->GetRMS());
        carte.m_histo_Billes->Fit("f1");
        fitGaussBille[i - 1] = f1->GetParameter(2);

        tabEnergie[i - 1] = sqrt(carte.energie()/(carte.masse_Carte * 1 + carte.masse_Boule * carte.m_N + carte.masse_Bille * carte.m_M * carte.m_N));

        delete f1;
        //std::this_thread::sleep_for(std::chrono::seconds(3));
    }
    */

     //Un essai unique
    int nbBoule, nbBille;
    double temperature, fracVolume, fracMasse1, fracMasse2;
    cout << "nombre de boules ?" << endl;
    cin >> nbBoule;
    cout << "nombre de billes dans chaque boule ?" << endl;
    cin >> nbBille;
    cout << "temperature : vitesses initiales ?" << endl;
    cin >> temperature;
    cout << "fraction de volume : tailles relatives boules/billes ?" << endl;
    cin >> fracVolume;
    cout << "fraction masse sphere/boule ?"  << endl;
    cin >> fracMasse1;
    cout << "fraction masse boule/bille ?" << endl;
    cin >> fracMasse2;
    
    int nbCaseBoule, nbCaseBille;
    cout << "nombre case histogramme Boule" << endl;
    cin >> nbCaseBoule;

    cout << "nombre case histogramme Bille" << endl;
    cin >> nbCaseBille;


    Carte carte(nbBoule, nbBille, temperature, fracVolume, fracMasse1, fracMasse2,nbCaseBoule,nbCaseBille); // nbr boules, nbr billes, température (répartition de vitesses) , fraction du volume, fraction masse sphere/boule, fraction masse boule/bille
   
    double tempsMin, tempsMax, epsilon;
    cout << "temps de début d'enregistrement ?" << endl;
    cin >> tempsMin;
    cout << "temps fin d'enregistrement ?" << endl;
    cin >> tempsMax;
    cout << "epsilon : 0.000001 ?" << endl;
    cin >> epsilon;

    
    carte.iterer(tempsMin,tempsMax,epsilon);  // temps de commencement d'enregistrement, temps de simulation, erreur pour calculer la collision

    TCanvas c1("c1", "Candle Decay", 800, 600); //creer la fenetre de l'histograme
    c1.Divide(2, 1);
    c1.cd(1);

    carte.m_histo_Billes->Draw();

    TF1* f2 = new TF1("f2", "gaus", -.5, .5); //pour fit le resultat avec une gaussienne
    f2->SetParameters(carte.m_histo_Billes->GetMaximum(), carte.m_histo_Billes->GetMean(), carte.m_histo_Billes->GetRMS()); 
    carte.m_histo_Billes->Fit("f2"); //fit

    c1.cd(2);
    carte.m_histo_Boules->Draw(); //dessiner l'histo 1

    TF1* f1 = new TF1("f1", "gaus", -.5, .5); //idem ...
    // set initial parameters (not really needed for gaus)
    f1->SetParameters(carte.m_histo_Boules->GetMaximum(), carte.m_histo_Boules->GetMean(), carte.m_histo_Boules->GetRMS());
    carte.m_histo_Boules->Fit("f1");

    c1.Update();
    c1.Draw(); 

    delete f1;
    delete f2;
    
    //résultat dans la console ...


    /*
    cout << endl << endl;
    cout << "resultat" << endl;
    for (int i(8); i < 14; ++i) {
        cout << "i = " << 5 * i + 5 << endl;
        cout << "deviation boule = " << fitGaussBoule[i] << endl;
        cout << "deviation bille = " << fitGaussBille[i] << endl;
        cout << "deviation vitesse = " << tabEnergie[i] << endl;
        cout << endl;
    }
    */


    app.Run();


    return 0;
}
