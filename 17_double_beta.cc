// Tue Oct 3, 2024
// Rayhaneh Dehghani; Master's student 
// Queen's University 
// SNO+ Experiment
//
// This script is for the purpose of making a stack plot of the response function convolved with the 2-neutrino 2-beta spectrum of 
// Tellurium-130. This signal is then overlayed with the unconvolved 2-neutrino 2-beta spectrum along with other bbackground signals
// in the energy region of 0.1-3 MeV.
//

#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TLatex.h"
#include "time.h"
#include "TTree.h"
#include "TPad.h"
#include "TStyle.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

using namespace std;

int tellurium()
{
	
	// defining the histogram to fill the spectra with
	TH1D *Energy = new TH1D("Stack_Energy", "2n2b_Energy", 500, 0.0, 5);
        Energy -> GetXaxis() -> SetTitle("2n2b_MC_Energy_spectrum");
        Energy -> GetYaxis() -> SetTitle("event_count");

	TH1D *TeConv_hist = new TH1D("convolved tellurium", "Total", 500, 0.0, 5);
	TeConv_hist -> GetXaxis() -> SetTitle("Energy [MeV]");
	TeConv_hist -> GetYaxis() -> SetTitle("event_count");
	// The Path for all the MC background and all other files
	//-------------------------------------------------------------------------------------------
        TString Directory = "/data/snoplus/production/scint/7.0.9/ntuple/ScintFit_2p2Te130_2n2bRun/";
        TString filename = "ScintFit_2p2Te130_2n2bRun_r";
        TString filename_2 = "_s0_p0.ntuple.root";
	// Include the path for other background as well

	TFile *file1 = TFile::Open("/home/rdehghani/Thesis/BiPo_studies/Te1302n2b_spectrum_2p2_optics.root");//MC 2n2b energy spectrum
	TFile *file2 = TFile::Open("/home/rdehghani/Thesis/BiPo_studies/Te130_0nu_spectrum_2p2_optics.root");//MC 0nu2b energy spectrum
	TFile *file3 = TFile::Open("/home/rdehghani/Thesis/BiPo_studies/B8_nu_spectrum_2p2_optics.root");// MC B8 nu spectrum
	TFile *file4 = TFile::Open("/home/rdehghani/Thesis/BiPo_studies/Tl208_nu_spectrum_2p2_optics.root");
	TFile *file5 = TFile::Open("/home/rdehghani/Thesis/BiPo_studies/BiPo212_Leakage_energy_spectrum_2p2_optics.root");
	//-------------------------------------------------------------------------------------------
	TH1D *Te130_2n2b_Energy;
	TH1D *Te130_0nu_Energy;
	TH1D *B8_nu_energy;
	TH1D *Tl208_Energy;
	TH1D *EnergyBiPo212;

	//Get the histograms from files
	Te130_2n2b_Energy = (TH1D*)file1 -> Get("Te130_2n2b_Energy");
	Te130_0nu_Energy =  (TH1D*)file2 -> Get("Te130_0nu_Energy");
        B8_nu_energy = (TH1D*)file3 -> Get("B8_Energy");
	Tl208_Energy = (TH1D*)file4 -> Get("Tl208_Energy");
	EnergyBiPo212 = (TH1D*)file5 -> Get("BiPo212_Energy");	
	//-------------------------------------------------------------------------------------------
	int i;
	int nbins = 500;
	float Te[nbins];
	float Te_Conv[nbins];
	float ResponseF[nbins];
	
	for(i = 0; i<nbins; i++) Te[i] = Te130_2n2b_Energy->GetBinContent(i+1);
        //-------------------------------------------------------------------------------------------
	// defining the response function
	double A = 0.0001;
        double B = 0.007;
        for( int i=0; i<nbins; i++){
                if (i <60){
        ResponseF[i] = A *exp((-1*i) * B);
                }
                else ResponseF[i]=0;
        }
        ResponseF[0] = 1-A;

        cout << "A =" << A << endl;
        cout << "B =" << B << endl;
        cout << "You passed the second barrier; taken bin cntents into arrays" << endl;
        //-------------------------------------------------------------------------------------------
	int m = 500;
	int n = 500;

	 for (int i=0 ; i < m; i++)//convolved histogram loop
        {
                double C = 0;
                for (int j= 0 ; j <= i ; j++)//MC loop, since "j" is the MC index
                {
                if (i-j >= 0) //check to omly get the possitive side of the response function
                {
                        C += Te[j] * ResponseF[i - j];//getting the parameters and adding to the previous one
                }
                }
                TeConv_hist->SetBinContent(i+1 , C);
        }//the end of the Perform convolutoin loop

	
	Long64_t numEvents_0nu = Te130_0nu_Energy -> Integral();
        cout << "Number of events in histogram " << Te130_0nu_Energy << ": " << numEvents_0nu << endl;

        Long64_t numEvents_B8 = B8_nu_energy -> Integral();
        cout << "Number of events in B8 histogram = " << numEvents_B8 << endl;

	Long64_t numEvents_Tl208 = Tl208_Energy -> Integral();
	cout << "Number of events in Tl208 histogram = " << numEvents_Tl208 << endl;

        Long64_t numEvents_bipo212Leakage = EnergyBiPo212 -> Integral();
	cout << "Number of events in BiPo212_Leakage histogram = " << numEvents_bipo212Leakage << endl; 

	Long64_t numEvents_twonu = Te130_2n2b_Energy -> Integral();
	cout << "Number of Events in Te130_2n2b = " << numEvents_twonu << endl;
	 
	// Current background rates
	TeConv_hist -> Scale(4344789.56/numEvents_twonu);//
        Te130_2n2b_Energy -> Scale(4344789.56/numEvents_twonu);//6.05254
        Te130_0nu_Energy -> Scale(76.4506/numEvents_0nu);//46.8727
        B8_nu_energy -> Scale(1021.4251/numEvents_B8);//
        Tl208_Energy -> Scale(16137.5569/numEvents_Tl208);
        EnergyBiPo212 -> Scale(2282.9266/numEvents_bipo212Leakage);
	
	//----------Set fill style to transparent
	gStyle->SetFillStyle(1001);
	//-------Setline color of each spectrum----
	TeConv_hist->SetLineColor(kRed);
	TeConv_hist->SetLineWidth(2);

	Te130_2n2b_Energy->SetLineColor(kBlue);
	Te130_2n2b_Energy->SetLineWidth(2);

	Te130_0nu_Energy->SetLineColor(kGreen);
	Te130_0nu_Energy->SetLineWidth(2);

	B8_nu_energy->SetLineColor(kCyan);
	B8_nu_energy->SetLineWidth(2);

	Tl208_Energy-> SetLineColor(kBlack);
	Tl208_Energy->SetLineWidth(2);

	EnergyBiPo212-> SetLineColor(kMagenta);
	EnergyBiPo212->SetLineWidth(2);

	TCanvas *c4 = new TCanvas("c4", "c4", 800, 600);
        c4 -> SetTickx();
        c4 -> SetTicky();
	
	// add each physocal proccess to the same histogram
	Te130_2n2b_Energy -> Draw("HIST");
	EnergyBiPo212 -> Draw("HIST same");
	Te130_0nu_Energy -> Draw("HIST same");
	TeConv_hist -> Draw("HIST same");
	B8_nu_energy-> Draw("HIST same");
	Tl208_Energy -> Draw("HIST same");
	

	TLegend *legend = new TLegend(0.6, 0.6, 0.9, 0.9);
        legend -> AddEntry(TeConv_hist, "Convolved Te1302n2b spectrum");
        legend -> AddEntry(Te130_2n2b_Energy, "Te1302n2b spectrum");
        legend -> AddEntry(Te130_0nu_Energy, "0nu energy spectrum");
	legend -> AddEntry(Tl208_Energy, "Tl208 energy spectrum");
	legend -> AddEntry(EnergyBiPo212, "BiPo212_Leakage");
	legend -> AddEntry(B8_nu_energy, "B8_solar");
	legend -> Draw();
        gStyle -> SetOptStat(0);
        gStyle -> SetOptTitle(0);
/*
	TFile* TeConv_file = new TFile("TeConvolved_2n2b_spectrum_2p2_optics.root", "Recreate");//done
	TeConv_hist-> Write();
*/
return 0;
}
	




