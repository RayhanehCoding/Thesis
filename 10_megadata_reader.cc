// Author: Rayhaneh Dehghani
// Description: This script is for reading data files and plotting the rlevant histograms of the BiPo214 finder MegaData files.
// Run: with .L 10_megadata_reader in root
// Last update: 8/Jan/2023 

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

void HistoMaker()
{
	//defining Histograms:
	
	//NhitsC vs. position_r
	//2D histogram for Pollonium 214 nhitsC vs postion r; How many PMT hits per single event.
	TH2D *Po_nhitscleaned = new TH2D("Po214_NhitsC", "Po214_NhitsC", 200, 0, 450, 200, 0, 6000);
	Po_nhitscleaned -> GetXaxis() -> SetTitle("NhitsC_Po214");
	Po_nhitscleaned -> GetYaxis() -> SetTitle("Event_position_r");
	Po_nhitscleaned -> SetTitle("Po214_NhitsC_May_2022_#300032_301122_deltaR_2500");
	

	//2D histogram for Bismuth 214 nhitsC vs position r; How many PMT hits per single event.
	TH2D *Bi_nhitscleaned = new TH2D("Bi214_NhitsC", "Bi214_NhitsC", 200, 300, 1000, 200, 0., 6000);
	Bi_nhitscleaned -> GetXaxis() -> SetTitle("NhitsC_Bi214");
	Bi_nhitscleaned -> GetYaxis() -> SetTitle("Event_position_r");
	Bi_nhitscleaned -> SetTitle("Bi214_NhitsC_May_2022_#300032_301122_deltaR_2500");
	// --------------------------------------------------------------------------
	// NHitsC vs. position x
	TH2D *Po_nhitscleaned_xz = new TH2D("Po_nhitXZ","Po_nhitXZ",200, -6000,6000, 200, -6000, 6000);
	Po_nhitscleaned_xz -> GetXaxis() -> SetTitle("Nhits_x");
	Po_nhitscleaned_xz -> GetYaxis() -> SetTitle("Nhits_z");
	Po_nhitscleaned_xz -> SetTitle("Po_nhitXZ_May_2022_#300032_301122_deltaR_2500");
	// --------------------------------------------------------------------------
	
	//Energy; 1D histogram for Pollonium 214 Energy 
	TH1D *Po_energy = new TH1D("Po214_Energy", "Po214_Energy", 500, 0, 5);
	Po_energy -> GetXaxis() -> SetTitle("Po214_Energy_spectrum");
	Po_energy -> GetYaxis() -> SetTitle("Event_Count");
	Po_energy -> SetTitle("Po214_Energy_May_2022_#300032_301122_deltaR_2500");
	//1D histogram for Bismuth 214 Energy
	
	TH1D *Bi_energy = new TH1D("Bi214_Energy", "Bi214_Energy", 500, 0, 5);
	Bi_energy -> GetXaxis() -> SetTitle("Bi214_Energy_spectrum");
	Bi_energy -> GetYaxis() -> SetTitle("Event_Count");
	Bi_energy -> SetTitle("Bi214_Energy_May_2022_#300032_301122_deltaR_2500");

	// --------------------------------------------------------------------------
	//Decay Time of each tags; 1D line graph of the deacy tie of all the tagged BiPos(214)
	TH1D *BiPo_Decay_Time = new TH1D("BiPo_decay_time", "BiPo_decay_time", 100, 0, 5000000);
	BiPo_Decay_Time -> GetXaxis() -> SetTitle("Decay_Time");
	BiPo_Decay_Time -> GetYaxis() -> SetTitle("Event_Count");
	BiPo_Decay_Time -> SetTitle("BiPo_decay_time_May_2022_#300032_301122_deltaR_2500");
	// --------------------------------------------------------------------------
	// Delta_r of tagged events
	TH1D *BiPo_Delta_r = new TH1D("BiPo_Delta_r","BiPo_Delta_r", 100, 0, 4000);
	BiPo_Delta_r -> GetXaxis() -> SetTitle("Delta_r");
	BiPo_Delta_r -> GetYaxis() -> SetTitle("Event_Count");
	BiPo_Delta_r -> SetTitle("BiPo_Delta_r_May_2022_#300032_301122_deltaR_2500");
	//Events per month; 2D line graph; Plots events in each day vs. the date 
	TH2D *Monthly_events = new TH2D("Monthly_events", "Monthly_events", 1000, 0,2000, 500, 0,600);
	Monthly_events -> GetXaxis() -> SetTitle("Run Number");
	Monthly_events -> GetYaxis() -> SetTitle("event count");
	//Monthly_events -> SetTitle("Monthly_events_June_2022_#301123_301306_deltaR_1000");

	// defining the path to the MegaData files:
	
	TString Directory = "/home/rdehghani/Thesis/BiPo_studies/BiPo_Files/Megadata_2022/";
	TString Filename = "Megadata_BiPo214_2022_File_";
	TString Filetype = "_long_halflife.txt";
	TString Filepath;

	vector<double> Po_posr_vec, Po_posx_vec, Po_posy_vec, Po_posz_vec, Po_energy_vec, Po_nhitsC_vec, Po_time_vec;
	vector<double> Bi_posr_vec, Bi_posx_vec, Bi_posy_vec, Bi_posz_vec, Bi_energy_vec, Bi_nhitsC_vec, Bi_time_vec;
	vector<double> decay_time_vec, delta_r_vec;
	//for (int Q=300468 ; Q<301122 ; Q=Q+1)
	for (int Q=300032 ; Q<301123 ; Q=Q+1) //uncomment for runs in May2022 2p2gpl optics
	//for (int Q=301123 ; Q<301923 ; Q=Q+1) //uncomment for runs in June2022 2p2gpl optics
	{
		char numstr[20];
		sprintf(numstr, "%d", Q);
		Filepath = Directory + Filename + numstr + Filetype ;
		cout << Filepath << endl ;

		ifstream input(Filepath);
		if (input.is_open())
		{
			cout << "The data file is open" << endl;
			double a ,b, c, d, e, f ,g, h, k, l, m, n, o, p, s, t;
			//reading the datafile into vectors;
			while (input >> a >> b >> c >> d >> e >> f >> g >> h >> k >> l >> m >> n >> o >> p >> s >> t)
			{
				//if (s<2500 && t<2624000)//16_half_life
				if (s<1000 && t<2624000)//16_half_life
				//if (s<700 && t<1312000)//8_Half_life
				{
				Po_energy_vec.push_back(a);
				Bi_energy_vec.push_back(b);
				Po_posr_vec.push_back(c);
				Bi_posr_vec.push_back(d);
				Po_posx_vec.push_back(e);
				Bi_posx_vec.push_back(f);
				Po_posy_vec.push_back(g);
				Bi_posy_vec.push_back(h);
				Po_posz_vec.push_back(k);
				Bi_posz_vec.push_back(l);
				Po_nhitsC_vec.push_back(m);
				Bi_nhitsC_vec.push_back(n);
				Po_time_vec.push_back(o);
				Bi_time_vec.push_back(p);
				delta_r_vec.push_back(s);
				decay_time_vec.push_back(t);
				}
			}
		//	cout << "Number of events is =" <<  Po_energy_vec.size() << endl;
			input.close();
		}
		else cout << "Error : Unable to input the data file" << endl;
	//	int event_count = Megadata_vec.size();
	//	cout << "event_count = " << event_count << endl;
		cout << "Number of events is =" <<  Po_energy_vec.size() << endl;
		//int EventCount = Po_energy_vec.size();	
	}
	//filling histograms;
	for (int i=0; i < Po_energy_vec.size(); i=i+1)
	{
		Po_energy -> Fill (Po_energy_vec[i]);
		Bi_energy -> Fill (Bi_energy_vec[i]);
		Po_nhitscleaned -> Fill (Po_nhitsC_vec[i] , Po_posr_vec[i]);
		Bi_nhitscleaned -> Fill (Bi_nhitsC_vec[i], Bi_posr_vec[i]);
		BiPo_Decay_Time -> Fill (decay_time_vec[i]);
		BiPo_Delta_r -> Fill (delta_r_vec[i]);
		Po_nhitscleaned_xz -> Fill (Po_posx_vec[i], Po_posz_vec[i]);
		//Monthly_events -> Fill (EventCount , Q ) ;
	}
//-------------------------
//Defining fit functions//
TF1 *Fit_gaus = new TF1("Fit_gaus", "gaus", 0.6, 1);
Fit_gaus -> SetLineColor(kBlack);

TFile* Po_energy_file = new TFile("Rebin_BiPo214_May_2p2_optics_1000mm_16HL.root", "Recreate");//done


Po_energy -> Sumw2();
Bi_energy -> Sumw2();

//Ploting histograms;	
//-------------------------
TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
c1 -> SetTickx();
c1 -> SetTicky();
c1 -> SetGridx();
c1 -> SetGridy();
Po_energy -> Draw();
//Po_energy -> SetFillColor(kViolet);
//Po_energy -> Fit ("Fit_gaus");
Po_energy -> Scale(1.0/Po_energy->Integral());
Po_energy -> Write();
//-------------------------
TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
c2 -> SetTickx();
c2 -> SetTicky();
c2 -> SetGridx();
c2 -> SetGridy();
Bi_energy -> Draw();
//Bi_energy -> SetFillColor(kBlue);
Bi_energy -> Scale(1.0/Bi_energy ->Integral());
Bi_energy -> Write();
//-------------------------
TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
c3 -> SetTickx();
c3 -> SetTicky();
c3 -> SetGridx();
c3 -> SetGridy();
Po_nhitscleaned -> Draw("colz");
Po_nhitscleaned -> Scale(1.0/Po_nhitscleaned ->Integral());
Po_nhitscleaned -> Write();
//-------------------------
TCanvas *c4 = new TCanvas("c4", "c4", 800, 600);
c4 -> SetTickx();
c4 -> SetTicky();
c4 -> SetGridx();
c4 -> SetGridy();
Bi_nhitscleaned -> Draw("colz");
Bi_nhitscleaned -> Scale(1.0/Bi_nhitscleaned ->Integral());
Bi_nhitscleaned -> Write();
//-------------------------
TCanvas *c5 = new TCanvas("c5", "c5", 800, 600);
c5 -> SetTickx();
c5 -> SetTicky();
c5 -> SetGridx();
c5 -> SetGridy();
BiPo_Decay_Time -> Draw();
BiPo_Decay_Time -> Scale(1.0/BiPo_Decay_Time ->Integral());
//BiPo_Decay_Time -> SetFillColor(kViolet-9);
BiPo_Decay_Time -> Write();
//-------------------------
TCanvas *c6 = new TCanvas("c6","c6",800, 600);
c6 -> SetTickx();
c6 -> SetTicky();
c6 -> SetGridx();
c6 -> SetGridy();
BiPo_Delta_r -> Draw();
BiPo_Delta_r -> Scale(1.0/BiPo_Delta_r ->Integral());
//BiPo_Delta_r -> SetFillColor ( kRed-9 );
BiPo_Delta_r -> Write();
//-------------------------
TCanvas *c7 = new TCanvas("c7", "c7", 800, 600);
c7 -> SetTickx();
c7 -> SetTicky();
c7 -> SetGridx();
c7 -> SetGridy();
Po_nhitscleaned_xz -> Draw("colz");
Po_nhitscleaned_xz -> Scale(1.0/Po_nhitscleaned_xz ->Integral());
Po_nhitscleaned_xz -> Write();
//-------------------------
TCanvas *c8 = new TCanvas("c8", "c8", 800, 600);
Po_energy -> Draw();
Bi_energy -> Draw("SAME");

return 0; 
} 




