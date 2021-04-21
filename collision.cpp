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
    double fitGaussBoule[12];
    double fitGaussBille[12];
    double tabEnergie[12];

    TApplication app("app", &argc, argv);
    for (int i = 1; i <= 12; ++i) {


        //TCanvas* canvas = new TCanvas("fCanvas", "fCanvas", 600, 400);
        //TCanvas c1("c1", "Candle Decay", 800, 600);
        //c1.Divide(2, 1);

        //TH1D h("h", "h", 10, 0, 10);
        //h.Fill(1);
        //h.Draw();

        //Carte(int n, int m, double temperature, double fraction, double rapport_masse, double rayon);
        cout << "carte numero " << i << endl;
        Carte carte(15, 5*i, 0.5, 3., 3., 1);
        carte.iterer(5, 200 , 0.00001);

        /*c1.cd(1);
        carte.m_histo_Billes->Draw();
        c1.cd(2);
        carte.m_histo_Boules->Draw();

        c1.Update();
        c1.Draw();*/

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
    cout << endl << endl;
    cout << "resultat" << endl;
    for (int i(0); i < 12; ++i) {
        cout << "i = " << 5 * i + 5 << endl;
        cout << "deviation boule = " << fitGaussBoule[i] << endl;
        cout << "deviation bille = " << fitGaussBille[i] << endl;
        cout << "deviation vitesse = " << tabEnergie[i] << endl;
        cout << endl;
    }

    app.Run();


    return 0;
}
