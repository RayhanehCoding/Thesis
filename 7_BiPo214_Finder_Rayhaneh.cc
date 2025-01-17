// \Authur : Rayhaneh Dehghani --- Master's program --- Queen's University --- SNO+ Experiment
//
// \This code will look for 214_Bismuth_Polonuim tagged coincidences 
// \and saves them into .txt files
//
// \Important: This code only looks at the high energy tail of the prompt and late eventrs, (energy >= 1 M eV)
//
//Run the code using:
//root
//.L 7_BiPo214_Finder_Rayhaneh.cc+ 


#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TLatex.h"
#include <time.h>
#include <iostream>  // cout ..
#include <fstream>  // ofstream
#include <math.h>
#include <string>
#include <TString.h>

#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TChain.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <map>
#include <ostream>
#include "TApplication.h"
#include <stdlib.h>
#include <stdio.h>
#include <TMath.h>
#include <vector>
#include <iterator>
#include <TEllipse.h>
#include <TGaxis.h>

///////////////////////////////////////////////////////////////////////////////////////////
                                                                                      
                                                                                      
//                        last update : Oct 18  2022                                   
                                                                                       
                                                                                       
///////////////////////////////////////////////////////////////////////////////////////////

// First we need to find the Bismuth 214 events (prompt event) and then look for events with
// a specific energy that occurs some specific time after the
// Bismuth 214 is found (late event), aka, Polonium 214.

// 1;
// Defining the data path:
// 2; 
// Specify your preferred cuts on "position" and "energy" and "nhits" and the "fiducial volume";
// 3; 
// Define the main loop
// 4;
// Define all the necessary variables and vectors for the output file
// 5;
// Define all te variables to be used
// 6;
// Define all of the branches you need to use
// 7;
// check to see if a run is more than 30 minutes
// 8;
// begin the serach for BiPo tags

	using namespace std;

	int mother_loop()

	{

	clock_t tStart = clock();
	//start = clock();
//---------------------------- Plotting for Po_214 -------------------------------------- 

	// This Plot tells us where in the detctor had how many PMT hits for a single Bi214 event
	TH2D *Po_nhitsCleaned = new TH2D("Po214_nhitsCleaned", "Po214_nhitsCleaned", 200, 0, 450, 200, 0., 6000);
	Po_nhitsCleaned -> GetXaxis() -> SetTitle("nhitsCleaned Po-214");
	Po_nhitsCleaned -> GetYaxis() -> SetTitle("Pos-r Po-214");

	// Plotting a 3D histogram for the purpose of seeing what event occured where; 
	TH3D *Po_Total_3D  = new TH3D("3D_Po-214","Po214_overall_deacyT_nhits_deltaR", 200, 0, 400, 200, 0, 5000, 200, 0, 500000);
	Po_Total_3D -> GetXaxis() -> SetTitle("nhitsCleaned Po-214");
	Po_Total_3D -> GetYaxis() -> SetTitle("Delta_r");
	Po_Total_3D -> GetZaxis() -> SetTitle("Decay_Time");

	// Plotting a 3D histogram as a projection of the distribution of events in the detector:
	TH3D *Po214_Pos_distribution_3D = new TH3D("Po214_Pos_distribution", "3D Po-214", 200, -6000, 6000, 200, -6000, 6000, 200, -6000, 6000);
	Po214_Pos_distribution_3D  -> GetXaxis() -> SetTitle("posx");	
	Po214_Pos_distribution_3D  -> GetYaxis() -> SetTitle("posy");	
	Po214_Pos_distribution_3D  -> GetZaxis() -> SetTitle("posz");

	// Plotting for the Polonium 214 energy spectrum:
	TH1D *Po_energy_plot  = new TH1D("Po214_energy_distibution", "Po214_energy_distibution", 500, 0, 10);
	Po_energy_plot-> GetXaxis() -> SetTitle("Po214_energy");
	//Po_energy_plot -> GetYaxis() -> SetTitle("Po_nhitsCleaned");

	// Plotting for the Bismuth 214 energy spectrum:	
	TH1D *Bi_energy_plot  = new TH1D("Bi214_energy", "Bi214_energy_distibution", 500, 0, 10);
	Bi_energy_plot -> GetXaxis() -> SetTitle("Bi214_energy");
	//Bi_energy_plot -> GetYaxis() -> SetTitle("Bi214_nhitsCleaned");

	// The posx vs. posr projection of the nhitsCleaned distribution:
	TH2D *Po_nhits_distribution_xr = new TH2D("Po214_xr","Po214_nhitsCleaned_distribution_xr", 300, 0., 6000, 300, 0., 6000);
	Po_nhits_distribution_xr -> GetXaxis() -> SetTitle("pos-x");
	Po_nhits_distribution_xr -> GetYaxis() -> SetTitle("pos-r");

	// The posx vs. posz projection of the nhitsCleaned distribution:
	TH2D *Po_nhits_distribution_xz = new TH2D("Po214_xz", "Po214_nhitsCleaned_distribution_xz", 300, -6000, 6000, 300, -6000, 6000);
	Po_nhits_distribution_xz -> GetXaxis() -> SetTitle("pos-x");
	Po_nhits_distribution_xz -> GetYaxis() -> SetTitle("pos-z");

	// The Decay time shows how long it takes a single Bi214(Beta) to decay into a Po214 alpha:
	TH1D *Po_214_Decay_Time  = new TH1D("Po214_Decay_Time", "Po214_Decay_Time", 200, 0., 5000000);	
	Po_214_Decay_Time -> GetXaxis() -> SetTitle("Po214_Deacy_time");

	// Delta_r is the distance from the first appearance of the Beta signal to when the corresponding Po214 is detected: 
	TH1D *BiPo214_Delta_r = new TH1D("BiPo214_Delta_r", "BiPo214", 200, 0., 4000);
	BiPo214_Delta_r -> GetXaxis() -> SetTitle("Delta_r");

	// Plotting the number of PMT hits for each Po214 event:
	TH1D *Po_nhits_plot = new TH1D("Po214_nhits_distribution_per_event", "Po214_nhits_distribution_", 200, 0., 450);
	Po_nhits_plot -> GetXaxis() -> SetTitle("Po214_nhitsCleaned");

//---------------------------- Plotting for Bi_214 --------------------------------------

	// This Plot tells us where in the detctor had how many PMT hits for a single Bi214 event;
	TH2D *Bi_nhitsCleaned = new TH2D("Bi214_nhitsCleaned_per_event", "Bi214", 200, 300, 1000, 200, 0., 6000);
	Bi_nhitsCleaned -> GetXaxis() -> SetTitle("Bi214_nhitsCleaned");	
	Bi_nhitsCleaned -> GetYaxis() -> SetTitle("Bi214_posr");

	// The nhitsC distribution for Bismuth 214 :
	TH1D *Bi_nhits_plot = new TH1D("Bi_nhits_plot", "Bi214", 200, 300, 1000);
	Bi_nhits_plot -> GetXaxis() -> SetTitle("1D nhitsCleaned");

//---------------------------------------------------------------------------------------
	//------------Defining a path to store the energy values inot.------------------------
        TString Data_Storage_path = "/home/rdehghani/Thesis/Data_files/";
        TString Energy_Output_Path = "/home/rdehghani/Thesis/Data_files/Energy_Data/" ;
	
	//----------- Define a path to get the desired run numbers from ----------------
	
//	TString run_files = "runlist_Rayhaneh_2.txt"; //the file with just the run number list in it
//	TString run_files = "run_list.txt"; // this is a runlist only for small tests.
	TString run_files = "run_list_2.txt"; 
	TString run_list_1 = "/home/rdehghani/R2/";

        // Defining the path to get each ntuple data file for each run:
      //  TString ntuple_path_1 = "/data/snoplus/processed_data/Scintillator/official_process/7.0.0/ntuple/";
	TString ntuple_path_1 = "/data/snoplus/processed_data/scint/7.0.4/ntuple/";
       	TString ntuple_path_2 = "Analysis30_r0000"; 
       	TString ntuple_path_3 = "_s0*";	

	//Defining a path to store the found data into:
//	TString   = Data_Storage_path + [run-number] + 

	// This is the main loop. 
	// Each loop will go through a single run
//	for (int Q = 271672; Q < 271673 ; Q = Q + 1) // if I'm using only a range of run numbers 

        TString Run_List_path = run_list_1 + run_files;
	vector <int> run_vec;
	ifstream input(Run_List_path);
        // open the .txt file with the run range in it:
	if (input.is_open())
        {
               // cout << Data_directory_path << endl;
                cout << " Run List path is = " << Run_List_path << endl ;  
                int index;
                while (input >> index)
                        {
                        run_vec.push_back(index) ;
                        }
		cout << " vector size is: " << run_vec.size() << endl; 
                input.close();
        }

	else cout << " Ntuple file does not exist " << endl; 

//	char numstr[20];
//	sprintf (numstr, "%d", index);
	//Defining a path to store the found data into:
	
	TString Po_energy_output = "_Po214_energy.txt";
        TString Po_eventID_output = "_Po214_eventID.txt";
        TString Po_posx_output  = "_Po214_posx.txt";
        TString Po_posy_output  = "_Po214_posy.txt";
        TString Po_posz_output  = "_Po214_posz.txt";
        TString Po_posr_output  = "_Po214_posr.txt";
        TString Po_nhitsC_output = "_Po214_nhitsCleaned.txt";

        TString Bi_nhitsC_output = "_Bi214_nhitsCleaned.txt";
        TString Bi_posx_output  = "_Bi214_posx.txt";
        TString Bi_posy_output  = "_Bi214_posy.txt";
        TString Bi_posz_output  = "_Bi214_posz.txt";
        TString Bi_posr_output  = "_Bi214_posr.txt";
        TString Bi_eventID_output = "_Bi214_eventID.txt";
        TString Bi_energy_output = "_Bi214_energy.txt";
	
	int runs = (int)run_vec.size(); // geting the number of runs that we want to proccess
	
	//**************************************************************************************
	
	/* The main loop that loops through run numbers in the .txt file one by one
	 for each run, the loop will check and look for BiPo's
	 and fills in the histograms 	
	*/

//	ofstream outfile;
//        outfile.open("BiPo214_Data_May.txt");	
	for (int run_number=0; run_number < runs; run_number++ )
//	for (int run_number=0; run_number < 2; run_number++ )

	{
		char nameString [20];
		sprintf (nameString,  " Output %d .txt", index) ;
		ofstream outfile;
	        outfile.open(nameString);

		// define different names
		//int Q; // the run nuumber
		char numstr[20];
		sprintf(numstr, "%d", run_vec[run_number]);
		TString Data_directory_path = ntuple_path_1 + ntuple_path_2 + numstr + ntuple_path_3 ;
		TString Po_energy_path =  Data_Storage_path + numstr + Po_energy_output ;

		//Run_List_path = run_list_1 + run_files ;
	        
		cout << Data_directory_path << endl; // ********************* 5 ****************
		cout << Po_energy_path << endl;

		cout << " Processing run " << index << endl; //***************** 6 ***************

		TChain *chain = new TChain ("output");
		chain -> Add (Data_directory_path);	
		chain -> Add (Po_energy_path);
		// --------------------------- ----------------------------------------
		//--------------2--Specify your preferred nhit cuts: ----------------------

		int BL = 400 ; 		// the Lower cut on Bismuth_214 nhits 
		int BH = 1000; 		// the Higher cut on Bismuth_214 nhits
		int PL = 160; 		// the Lower cut on Polonium_214 nhits
		int PH = 400; 		// the Higher cut on Polonium_214 nhits
		int HL = 8 ;            // the multiple of half-life (Half Life)

		int Bi_Low = BL ;
		int Bi_High = BH ;
		int Po_Low = PL ;
		int Po_High = PH;
		
		int FV_cut = 5500 ;       // [milimeter]  // Specify the Fiducial Volume cut:
		double half_life = HL * 164000;  // [nano seconds]   // variacble for how many half lives we are looking into

		//------------------------3-- Create vectors of the output parameters -------------
		//each vector element corresponds to an individual event.
			
		// Polonuim 214	
		vector<double>  	    Po214_posX_vec;
		vector<double> 		    Po214_posY_vec;
		vector<double>              Po214_posZ_vec;
		vector<double> 		    Po214_posr_vec;
		vector<double> 		  Po214_energy_vec;
		vector<double> 		Po214_momentum_vec;
		vector<double> 		   Po214_nhits_vec;
       		vector<double> 	    Po214_nhitsCleaned_vec;
		vector<int>    		 Po214_eventID_vec;
		
		// Bismuth 214
		vector<int>    		 Bi214_eventID_vec;
		vector<double> 		    Bi214_posX_vec;
		vector<double> 		    Bi214_posY_vec;
		vector<double> 		    Bi214_posZ_vec;
       		vector<double> 		    Bi214_posr_vec;
       		vector<double> 		  Bi214_energy_vec;
       		vector<double> 		Bi214_momentum_vec; // I want to look at the momentum of the outgoing particle.
		vector<double> 		   Bi214_nhits_vec;
		 //----------------------------------------5 & 6---------------------------------------
		 //-------------------- Define the variables and branches to be used- ----------------- 		

		 //TChain *chain2 = new TChain("output_data");

		double energy;
		chain -> SetBranchAddress("energy",             &energy);
		
		double skyShine;
		chain -> SetBranchAddress("skyShine",         &skyShine);
			
		double posr;
		chain -> SetBranchAddress("posr",                 &posr);	
		double posx;
		chain -> SetBranchAddress("posx",                 &posx);
		double posy;
		chain -> SetBranchAddress("posy",                 &posy);
		double posz;
		chain -> SetBranchAddress("posz",                 &posz);
				
		double dirx;
		chain -> SetBranchAddress("dirx",                 &dirx);
		double diry;
		chain -> SetBranchAddress("diry",                 &diry);
		double dirz;
		chain -> SetBranchAddress("dirz",                 &dirz);
				
		int nhits;
		chain -> SetBranchAddress("nhits",               &nhits);
		int nhitsCleaned;
		chain -> SetBranchAddress("nhitsCleaned", &nhitsCleaned);
			
		int eventID;
		chain -> SetBranchAddress("eventID",           &eventID);
		   	int runID;	//the signiture of each event and run in the dataset
		chain -> SetBranchAddress("runID",               &runID);
			
		int uTDays;
		chain -> SetBranchAddress("uTDays",             &uTDays);	
		int uTSecs;
		chain -> SetBranchAddress("uTSecs",             &uTSecs);
		int  uTNSecs; 	// gives you the time of each event or each run duration
		chain -> SetBranchAddress("uTNSecs",           &uTNSecs);
			
		bool fitValid;
		chain -> SetBranchAddress("fitValid",         &fitValid);
			
		bool partialFit;//if returns 1, shows the data passed the data cleaning cuts
		chain -> SetBranchAddress("partialFit",     &partialFit);
			
		bool scintFit;
		chain -> SetBranchAddress("scintFit",         &scintFit);
			
		ULong64_t clockCount50;	// the time of events
		chain -> SetBranchAddress("clockCount50", &clockCount50);
   	    		
		ULong64_t dcFlagged;
		chain -> SetBranchAddress("dcFlagged",       &dcFlagged);
			
		ULong64_t dcApplied;//data cleaning cuts
		chain -> SetBranchAddress("dcApplied",       &dcApplied);

		int entries = chain->GetEntries();
		//-------------------------------------------------------------------------------------
		//-----------------------------------------7;------------------------------------------
		// ---------------- Now check and see if a run is more than 30 minutes in duration-----
		// ------------------------------------------------------------------------------------
	
		//Check to see id a run is more than 30 minutes
		//	if (runtime > )

		//	}
		//-------------------------------------- 8 --------------------------------------------
		//------------------Now the search for Bi_Po_214 pairs begins -------------------------

		//--------------- Finding the Bismuth_214 first------------------
		int Bi_Found = 0; // this will tell you how many events you found
		int Po_Found = 0; // this will tell you how many Po214s you found
		ULong64_t Bi_Time;
		ULong64_t Po_Time;
		double offset; // to get the offset of the AV for every event 
		// this piece of code gets the value of the z-offset of the AV for every run, because it is constantly changing due to
                // AV recuiculation and adding PPO and stuff, (we are not in a stable data taking period.)	
				
		RAT::DB *db = RAT::DB::Get();
                RAT::DS::Run rundb;
                rundb.SetRunID(index);
                db->BeginOfRun(rundb);
                RAT::DBLinkPtr dblink2=db->GetLink("AV_OFFSET_RUNTIME");
                vector<double> off=dblink2->GetDArray("position");
                offset = off[2];
                cout << "offset is =" << offset << endl;


//		ofstream outfile;
//		outfile.open("Po_Data_December.txt");
		for(int BiFinder = 0; BiFinder < entries; BiFinder++)
			{
			chain->GetEntry(BiFinder);
			ULong64_t Bi_Time = clockCount50;
			double Bi_posx = posx;
			double Bi_posy = posy;
		 	double Bi_posz = posz;
			double Bi_posr = sqrt ( (Bi_posx**2) + (Bi_posy**2) + (Bi_posz-offset)**2 );
			int    Bi_nhits= nhitsCleaned;
			double Bi_energy = energy;	
		        double me_0 = 511000; //[eV]
                      	//ouble Bi_momentum;
			if (nhitsCleaned < BH 
				&& nhitsCleaned > BL
			//	&& partialFit == 1
			//	&& scintFit == 1
				&& fitValid == 1
	  		//	&& skyShine != -99999
	  			&& Bi_posr <= FV_cut && Bi_posr > 0
				&& Bi_energy != -99999
				// && dcFlagged ?
	  			    )
				{
					Bi_Found++;
					for (int PoFinder = BiFinder+1 ; PoFinder < entries; PoFinder++)
					{
						chain->GetEntry(PoFinder);
			      			ULong64_t Po_Time = clockCount50;
						double Po_posx = posx;
						double Po_posy = posy;
						double Po_posz = posz;
						double Po_posr = sqrt ( (Po_posx**2) + (Po_posy**2) + (Po_posz-offset)**2 );
						int    Po_nhits= nhitsCleaned;
					        double Po_energy = energy;
						double m_alpha = 3727300000; ///[eV]	

						if (nhitsCleaned < PH
                                               		&& nhitsCleaned > PL
                                             		&& fitValid == 1
						//	&& skyShine != -99999
							&& Po_posr <= FV_cut && Po_posr > 0
							&& Po_energy != -99999
							)
						{

						//	double Po_momentum = sqrt (Po_energy**2 - m_alpha**2);
							ULong64_t Decay_Time = (Po_Time - Bi_Time) * 20;
							//Decay time is the time window that we should expect a Po_214 event
							//after Bi_214 event is observed.
							//The time we find is in units of "nano seconds".
							//
							//Delta_r is the position cut on the Bi-Po coincidences, meaning how much apart
							//the two events are.
							double Delta_r = sqrt ( (Po_posx - Bi_posx)**2 + (Po_posy - Bi_posy)**2 + (Po_posz-offset - Bi_posz-offset)**2 ) ;
								
							// I need to loop through all of the tags and not only the ones that satisfy the 
							// if statement conditions for being plotted. 
							// Plotting and `saving data are 2 different things :)
							if(Decay_Time > 400  && Decay_Time <= half_life 
									     && Delta_r <= 2500
									     ) 
							{
								cout << " found one " << endl;
								Po_Found++;
								
								if (outfile.is_open())
								{
								outfile << Bi_energy <<"\t"<< Po_energy <<"\t"<< Bi_posr <<"\t"<< Po_posr <<"\t"<< Bi_nhits <<"\t"<< Po_nhits <<"\t"<<  endl;
								outfile << Bi_posx <<"\t"<< Po_posx <<"\t"<<  Bi_posy <<"\t"<< Po_Posy  <<"\t"<< Bi_posz <<"\t"<< Po_posz <<"\t"<< endl;
								outfile << Bi_Time  <<"\t"<< Po_Time <<"\t"<< Delta_r <<"\t"<< Decay_Time <<  endl;
									
								}
								else cout << "ERROR" << endl ;

								Po_214_Decay_Time 	  -> Fill (Decay_Time) ;
								BiPo214_Delta_r           -> Fill (Delta_r) ;
								Po_nhitsCleaned           -> Fill (Po_nhits, Po_posr) ;
								Po214_Pos_distribution_3D -> Fill (Po_posx, Po_posy, Po_posz) ;
								Po_energy_plot            -> Fill (Po_energy) ;
								Po_nhits_distribution_xz  -> Fill (Po_posx, Po_posz) ;
								Po_nhits_distribution_xr  -> Fill (Po_posx, Po_posr) ;


								Bi_energy_plot            -> Fill (Bi_energy) ;
								Bi_nhitsCleaned           -> Fill (Bi_nhits, Bi_posr) ;

								cout << " Delta_r = "     << Delta_r    << endl ;


							}//end of decay time condition
					

						else	PoFinder = entries;
						}
						//	} //if condition on Po214
					} //the end of PoFinder loop
				} //if condition on Bi214		
		//	delete chain;
		//	utput.close();	
				} //the end of BiFinder loop
				//Let's see how many events we found:
	//		 delete chain;
	//        }
//	outfile.close();
	cout << " You found " << Bi_Found << " Bismuth 214 events " << endl;
	cout << " You found " << Po_Found << " Polonium 214 events " << endl;

	if (Bi_Found == Po_Found)
	
	{	
		
		cout << " The numbers match and you have found " << Bi_Found << " pairs " << endl;
	}

	else {
		cout << " the numbers DO NOT match " << endl;
	}

delete chain;

outfile.close();

	}

//outfile.close();

//--------------------- Define your Fit functions ------------------

	TF1 *Fit_time = new TF1("Fit_time", "expo", 0, 10);

	TF1 *Fit_energy = new TF1("Fit_energy", "expo", 1, 10);


//----------------------Plotting for Polonium 214:------------------------------

	TFile *BiPo_214_tagged_pairs_high_energy_March2022  = new TFile("/home/rdehghani/Thesis/Data_files/BiPo_214_tagged_pairs_high_energy.root", "Recreate");
	
	BiPo_214_tagged_pairs_high_energy_March2022 -> cd();

	//TCanvas *c1 = new TCanvas("c1", "multipad", 1200, 1000);
	//c1 -> Divide(2,2,0,0);

	// c1 -> cd(1);
	// Plotting the decay time of 214_Po
	
	
	TCanvas *c2 = new TCanvas("c2", "multipad", 800, 600);
	c2 -> cd();
	c2 -> SetTickx();
	c2 -> SetTicky();
	c2 -> SetGridx();
	c2 -> SetGridy();
	Po_214_Decay_Time  -> Draw();
	Po_214_Decay_Time -> SetFillColor(kViolet-9);
	Po_214_Decay_Time  -> Write ();
	Po_214_Decay_Time  -> Fit("Fit_time");
	

	
//	Fit_time -> SetLineColor (  );
	TLegend *leg2 = new TLegend(0.7, 0.7, .9, .9);
	leg2 -> AddEntry (Po_214_Decay_Time, " Measured Data ", "l" );
	leg2 -> AddEntry (Fit_time, "Fit Function", "l");
	leg2 -> Draw();
	

	// c1 -> cd(2);
	// Plotting the Delta-r (the distace from 214-Bi to 214-Po)
	TCanvas *c3 = new TCanvas("c3", "multipad", 800, 600);
	c3 -> cd();
	c3 -> SetTickx();
	c3 -> SetTicky();
	c3 -> SetGridx();
        c3 -> SetGridy();
	BiPo214_Delta_r  -> SetFillColor ( kRed-9 );
	BiPo214_Delta_r  -> Draw();
	BiPo214_Delta_r  -> Write ();
	// c1 -> cd(3);

	// Plotting nhitsCleaned for 214-Polonium
	TCanvas *c4 = new TCanvas("c4", "multipad", 800, 600);
	c4 -> cd();
	c4 -> SetTickx();
	c4 -> SetTicky();
	c4 -> SetGridx();
        c4 -> SetGridy();
	Po_nhitsCleaned -> Draw ( "colz" );
	Po_nhitsCleaned -> Write ();

// Plotting nhitsC for xr plane;
	TCanvas *c5 = new TCanvas("c5", "multipad", 800, 600);
	c5 -> cd();
	c5 -> SetTickx();
	c5 -> SetTicky();
	c5 -> SetGridx();
        c5 -> SetGridy();
	Po_nhits_distribution_xr -> Draw ( "colz") ;
	Po_nhits_distribution_xr -> Write ();

// Plotting nhitsC for xz plane;
	TCanvas *c8 = new TCanvas("c8", "multipad", 800, 600);
	c8 -> cd();
	c8 -> SetTickx();
	c8 -> SetTicky();
	c8 -> SetGridx();
        c8 -> SetGridy();
	Po_nhits_distribution_xz ->  Draw ( "colz" );
	Po_nhits_distribution_xz -> Write();

// Plotting the the total info of Polonuim in 3D;
	TCanvas *c7 = new TCanvas("c7", "multipad", 800, 600);
	c7 -> cd();
	Po_Total_3D -> Draw ( "LEGO1" );
	Po_Total_3D -> Write ();

	/*

TCanvas *c11 = new TCanvas("c11", "multipad", 800, 600);
c11 -> cd();
Po_energy_plot -> Draw ("colz");
Po_energy_plot -> Write();
*/

	TCanvas *c9 = new TCanvas("c9", "multipad", 800, 600);
	c9 -> cd();
	c9 -> SetTickx();
	c9 -> SetTicky();
	gStyle -> SetPalette(1);
	Po214_Pos_distribution_3D -> Draw ( "LEGO1" );
	Po214_Pos_distribution_3D -> Write();

//----------------------Plotting for Bi-Po 214 energy:---------------------------------
	TCanvas *c20 = new TCanvas("energy spectrum of Po214", "multipad", 800, 600);
	//c20 -> Divide(1,2);
	c20 -> cd();
	c20 -> SetTickx();
        c20 -> SetTicky();
        c20 -> SetGridx();
        c20 -> SetGridy();
	//c20 -> cd();
	Po_energy_plot -> Draw ( );
	Po_energy_plot -> SetFillColor ( kViolet );
	Po_energy_plot -> Write ();
//	Po_energy_plot -> Fit ( "Fit_energy" );

//	Po_energy_plot -> Fit("gaus");

// ------------------	
	TCanvas *c21 = new TCanvas("energy spectrum of Bi214", "multipad", 800, 600);
	c21 -> cd();	
	c21 -> SetTickx();
	c21 -> SetTicky();
	c21 -> SetGridx();
	c21 -> SetGridy();
	Bi_energy_plot -> Draw ();
	Bi_energy_plot -> SetFillColor ( kBlue );
	Bi_energy_plot -> Write ();

// -----------------------------

	TCanvas *c11 = new TCanvas("c11", "multipad", 800, 600);
	c11 -> cd();
	c11 -> SetTickx();
	c11 -> SetTicky();
        c11 -> SetGridy();
        c11 -> SetGridx();
	Bi_nhitsCleaned -> Draw ( "colz" );
	Bi_nhitsCleaned -> Write ( );
	

	
//	TLegend *leg11 = new TLegend (0.7, 0.7, 0.9, 0.9);
	
//	leg11 -> Draw ();
/*

// Plotting 214_Bismuth nhitsCleaned distrubution
TCanvas *c10 = new TCanvas("c10", "multipad", 800, 600);
c10 -> cd();
Bi_nhitsCleaned -> Draw("colz");
Bi_nhitsCleaned -> Write();

//-------------------------------------------------------------------------------
TCanvas *c11 = new TCanvas("c11", "multipad", 800, 600);
c11 -> cd();
Po_energy_plot -> Draw ("colz");
Po_energy_plot -> Write();


TCanvas *c12 = new TCanvas("c12", "multipad", 800, 600);
c12 -> cd();
Bi_energy_plot -> Draw ("colz");
Bi_energy_plot -> Write();

*/

//	TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);

//	leg -> Draw();



// Priinting the execution time of my code:
//clock_t tStart = clock();

	//end = clock();
//	double time_taken = double (clock() - tStart) / (double) CLOCKS_PER_SEC ;
//	cout << "Execution time : " << time_taken << endl;
	//`cout << " sec " << endl;

// ---------------- saving the data ----------------
// I will save the data into root files and .txt files for furthur work



	
return 0;

}//close the mother_loop function


