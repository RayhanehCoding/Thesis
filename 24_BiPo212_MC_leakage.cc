// Author : Rayhaneh Dehghani
// Compile the code in Root; using ".L 18_zeronu_Te.cc " 
// For the purpose of the Master's degree thises

// This script is for finding the BiPo212 MC energy spectrum that leak past the cuts and appear as single energentic events
// 
// Selection cuts:  will be defined ; energy, Nhits_low, Nhits_high, FV, FitValid, event_index.
// Classifiers : alphaBeta212 




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

int bipo212()
{

	TH1D *EnergyBiPo212 = new TH1D("BiPo212_Energy", "BiPo212_Energy", 500, 0.0, 5);
        EnergyBiPo212 -> GetXaxis() -> SetTitle("BiPo212_MC_energy_spectrum");
        EnergyBiPo212 -> GetYaxis() -> SetTitle("event_count");

        TH2D *BiPo212_Energy_posr = new TH2D("BiPo212_Energy_posr","BiPo212_Energy_posr", 500, 0, 5, 600, 0, 6000);
        BiPo212_Energy_posr -> GetXaxis() -> SetTitle("Energy");
        BiPo212_Energy_posr -> GetYaxis() -> SetTitle("posr");
	
	TH2D *nhitscleaned_BiPo212 = new TH2D("nhitscleaned_BiPo212", "nhitscleaned_BiPo212", 500, 0, 1000, 600, 0, 6000);
        nhitscleaned_BiPo212 -> GetXaxis() -> SetTitle("NhitsC_BiPo212");
        nhitscleaned_BiPo212 -> GetYaxis() -> SetTitle("Event_position_r");

	//-------------------------------------------------------------------------------------------
	TString Directory = "/data/snoplus/production/scint/7.0.8/ntuple/ScintFit_2p2Bipo212Run/";
	TString filename = "ScintFit_2p2Bipo212Run_r";
	TString filename_2 = "_s0_p0.ntuple.root";
	//-------------------------------------------------------------------------------------------
        int NhL=20;
       	int Nhits_low;
	int NhH = 700;
        int Nhits_high ;
	//int F = 5500;//[mm]
        int FV_cut = 3300;//[mm]
	
        int nhitsCleaned;
        double EL = 0;//[Mev]
	double energy_low = EL;
	double EH = 2.75;//[Mev]
        double energy_high = EH;

	Double_t alphaBeta212;
	Double_t mcPosr;
        double posr, posx, posy, posz;
        bool fitValid;
	double energy;

	int passed_evIndex = 0;
        int passed_fitValid = 0;
        int passed_mcPosr = 0;
        int passed_posr = 0;
	int passed_alphaBeta212 = 0;
        //-------------------------------------------------------------------------------------------
	for (int Q=300000 ; Q<306498 ; Q=Q+1)
	//for (int Q=300000 ; Q<301000 ; Q=Q+1)
	{
	char numstr[100];
	sprintf(numstr, "%d", Q);
	TString MC_BiPo212_path = Directory + filename + numstr + filename_2 ;

	ifstream input(MC_BiPo212_path);
        if (input.is_open())
        {
	cout << MC_BiPo212_path << endl;

	TChain *chain = new TChain ("output");
        chain -> Add(MC_BiPo212_path);

	int evIndex;
	chain -> SetBranchAddress("evIndex",&evIndex);//Only for MC files to control re-trigger plotting
	chain -> SetBranchAddress("energy",&energy);
        chain -> SetBranchAddress("posr",&posr);
        chain -> SetBranchAddress("posx",&posx);
        chain -> SetBranchAddress("posy",&posy);
        chain -> SetBranchAddress("posz",&posz);
        chain -> SetBranchAddress("fitValid",&fitValid);
        chain -> SetBranchAddress("nhitsCleaned",&nhitsCleaned);
        chain -> SetBranchAddress("alphaBeta212",&alphaBeta212);
	chain -> SetBranchAddress("mcPosr",&mcPosr);
	int entries=chain->GetEntries();
	int Events = chain->GetEntries("evIndex <= 0");
	//-------------------------------------------------------------------------------------------
		for (int i=0; i <= entries; i=i+1)
		{
			chain->GetEntry(i);
			//use the right order of cuts to keep track of what happens

			// Condition #1
                        if (evIndex <= 0) {
                                passed_evIndex++; // Increment counter if condition is true

				// Condition #2
                                if  (fitValid ==1) {
                                        passed_fitValid++; // Increment counter if condition is true

					// Condition #3
					if (alphaBeta212 > -5) {
						passed_alphaBeta212++;

						// Condition #4
                                        	if (mcPosr <= 6000) {
                                                	passed_mcPosr++; // Increment counter if condition is true

							// Condition #5
                                                	if (posr <= FV_cut) {
                                                        	passed_posr++; // Increment counter if condition is true

								EnergyBiPo212 -> Fill(energy);//fill with whatever I want to fill the histogram; i.e energy, posr
								BiPo212_Energy_posr -> Fill(energy, posr);
								nhitscleaned_BiPo212 -> Fill(nhitsCleaned, posr);
							}
						}
					}
				}
			}// end of condition #1
		}
	input.close();
        delete chain;
	}//the end of the input file
	//sum all entriies from each chain	
	else cout <<"file does not exist" << endl;
	}	
	cout << "Number of events passing [evIndex] cut = " << passed_evIndex << endl;
        cout << "Number of events passing [evIndex and fitValid] cuts = " << passed_fitValid << endl;
	cout << "Number of events passing [evIndex and fitValid and alphaBeta212] cuts = " << passed_alphaBeta212 << endl;
        cout << "Number of events passing [evIndex and fitValid and alphaBeta212 and mcPosr] cuts = " << passed_mcPosr << endl;
        cout << "Number of events passing [evIndex and fitValid and alphaBeta212 and mcPosr and posr] cuts = " << passed_posr << endl;
	//print the counted numbers
        //----------------------------------------------------------------------------------------
	TFile* BiPo212_Leakage_energy_file = new TFile("BiPo212_Leakage_energy_spectrum_2p2_optics.root", "Recreate");//done

	EnergyBiPo212 -> Sumw2();

	//Ploting histograms;
	//-------------------------
	TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
	c1 -> SetTickx();
	c1 -> SetTicky();
	c1 -> SetGridx();
	c1 -> SetGridy();
	EnergyBiPo212 -> Draw();
	//Po_energy -> SetFillColor(kViolet);
	//Po_energy -> Fit ("Fit_gaus");
	//EnergyBiPo212->Scale(1.0/EnergyBiPo212->Integral());
	EnergyBiPo212 -> Write();
	//-------------------------
	TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
        c2 -> SetTickx();
        c2 -> SetTicky();
        c2 -> SetGridx();
        c2 -> SetGridy();	
	BiPo212_Energy_posr -> Draw("colz");
	BiPo212_Energy_posr-> Write();
	//-------------------------
	TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
        c3 -> SetTickx();
        c3 -> SetTicky();
        c3 -> SetGridx();
        c3 -> SetGridy();
	nhitscleaned_BiPo212 -> Draw("colz");
	nhitscleaned_BiPo212 -> Write(); 


return 0;


}
