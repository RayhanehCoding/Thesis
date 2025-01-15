// Author : Rayhaneh Dehghani
// Compile the code in Root; using ".L 18_zeronu_Te.cc " 
// For the purpose of the Master's degree thesis

// This script is for finding the 2-neutrino 2-beta signal in Monte Carlo for selected run ranges
// (that describe the detector state in terms of reconstruction and physics)
// Selection cuts will be defined ; energy, Nhits_low, Nhits_high, FV, FitValid, event_index, etc. 




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

int zeronu()
{

	TH1D *Energy0nu = new TH1D("Te130_0nu_Energy", "Te130_0nu_Energy", 500, 0.0, 5);
        Energy0nu -> GetXaxis() -> SetTitle("0nu_MC_Energy_spectrum");
        Energy0nu -> GetYaxis() -> SetTitle("event_count");

        TH2D *Energy_posr = new TH2D("0nEnergy_posr","0nuEnergy_posr", 500, 0, 5, 600, 0, 6000);
        Energy_posr -> GetXaxis() -> SetTitle("Energy");
        Energy_posr -> GetYaxis() -> SetTitle("posr");
	
	TH2D *nhitscleaned_0nuTe130 = new TH2D("nhitscleaned_0nuTe130", "nhitscleaned_0nuTe130", 500, 0, 1000, 600, 0, 6000);
        nhitscleaned_0nuTe130 -> GetXaxis() -> SetTitle("NhitsC_Te130_0nu");
        nhitscleaned_0nuTe130 -> GetYaxis() -> SetTitle("Event_position_r");
	//-------------------------------------------------------------------------------------------
	TString Directory = "/data/snoplus/production/scint/7.0.8/ntuple/ScintFit_2p2Te130_0nuRun/";
	TString filename = "ScintFit_2p2Te130_0nuRun_r";
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

	double alphaBeta212;
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
	vector<double> Te0nu_energy_vec, Te0nu_nhits_vec;	
		
	for (int Q=300000 ; Q<306498 ; Q=Q+1)
	{
	char numstr[100];
	sprintf(numstr, "%d", Q);
	TString MC_0nuTe130_path = Directory + filename + numstr + filename_2 ;
	//TString File_name = filename + numstr + filename_2 ;
		
	ifstream input(MC_0nuTe130_path);
        if (input.is_open())
        {
	cout << MC_0nuTe130_path << endl;

	TChain *chain = new TChain ("output");
        chain -> Add(MC_0nuTe130_path);

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
	//Double_t output = chain->GetEntries("evIndex <= 0");	
	//-------------------------------------------------------------------------------------------
		//cout << "The data file is open" << endl;
		for (int i=0; i <= entries; i=i+1)
		{
			chain->GetEntry(i);
			
			// Condition #1
			if (evIndex <= 0) {
				passed_evIndex++; // Increment counter if condition is true
	
				// Condition #2
				if  (fitValid ==1) {
					passed_fitValid++; // Increment counter if condition is true

					// Condition #3
                    if (alphaBeta212 > -5) {
                        passed_alphaBeta212++;

                        // Condition #3
                        if (mcPosr <= 6000) {
                            passed_mcPosr++; // Increment counter if condition is true

                            // Condition #4
                            if (posr <= FV_cut) {
                                passed_posr++; // Increment counter if condition is true
					

					Energy0nu -> Fill(energy);//fill with what ever I want to fill the histogram; i.e energy, posr
					Energy_posr -> Fill(energy, posr);
					nhitscleaned_0nuTe130 -> Fill(nhitsCleaned, posr);
						}
					}
				}			

			}//the end of first condition
		}
		}
		input.close();
		delete chain;
		}//the end of the for loop for entries
		else cout <<"file does not exist" << endl;
		}//the end of the input file 
	//else cout <<"file does not exist" << endl;
	cout << "Number of events passing [evIndex] cut = " << passed_evIndex << endl;
        cout << "Number of events passing [evIndex & fitValid] cuts = " << passed_fitValid << endl;
	cout << "Number of events passing [evIndex & fitValid & alphaBeta212] cuts = "<< passed_alphaBeta212 << endl;
	cout << "Number of events passing [evIndex & fitValid & alphaBeta212 & mcPosr] cuts = " << passed_mcPosr << endl;
	cout << "Number of events passing [evIndex & fitValid & alphaBeta212 & mcPosr & posr] cuts = " << passed_posr << endl;	
        //----------------------------------------------------------------------------------------
	TFile* Te130_0nu_energy_file = new TFile("Te130_0nu_spectrum_2p2_optics.root", "Recreate");//done

	Energy0nu -> Sumw2();

	//Ploting histograms;
	//-------------------------
	TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
	c1 -> SetTickx();
	c1 -> SetTicky();
	c1 -> SetGridx();
	c1 -> SetGridy();
	Energy0nu->Draw();
	//Po_energy -> SetFillColor(kViolet);
	//Po_energy -> Fit ("Fit_gaus");
	//Energy0nu->Scale(1.0/Energy0nu->Integral());
	Energy0nu -> Write();
	//-------------------------
	TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
        c2 -> SetTickx();
        c2 -> SetTicky();
        c2 -> SetGridx();
        c2 -> SetGridy();	
	Energy_posr->Draw("colz");
	Energy_posr-> Write();
	//-------------------------
	TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
        c3 -> SetTickx();
        c3 -> SetTicky();
        c3 -> SetGridx();
        c3 -> SetGridy();
	nhitscleaned_0nuTe130-> Draw("colz");
	nhitscleaned_0nuTe130 -> Write();
return 0;


}
