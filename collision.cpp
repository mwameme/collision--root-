#include <TApplication.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH1.h>
#include <thread>
#include <chrono>
#include "carte.h"

#include <iostream>
using namespace std;


//TCanvas canvas ("fCanvas", "fCanvas", 600, 400);

int main(int argc, char* argv[])
{
    TApplication app("app", &argc, argv);

    //TCanvas* canvas = new TCanvas("fCanvas", "fCanvas", 600, 400);
    TCanvas c1("c1", "Candle Decay", 800, 600);
    c1.Divide(2, 1);

    //TH1D h("h", "h", 10, 0, 10);
    //h.Fill(1);
    //h.Draw();

    //Carte(int n, int m, double temperature, double fraction, double rapport_masse, double rayon);
    Carte carte(15, 15, 0.5, 3, 4, 1);
    carte.iterer(3,200,0.000001);

    c1.cd(1);
    carte.m_histo_Billes->Draw();
    c1.cd(2);
    carte.m_histo_Boules->Draw();

    c1.Update();
    c1.Draw();

    //std::this_thread::sleep_for(std::chrono::seconds(3));
    app.Run();


    return 0;
}
