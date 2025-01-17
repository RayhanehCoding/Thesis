#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TLatex.h"
#include "time.h"
#include "TTree.h"
#include "TPad.h"
#include "TSpectrum.h"
#include "TSystem.h"

using namespace std;

void RF_Convolution()
{

	// load the histograms from seperate .root files
        // data histogram file
        TFile *file1 = TFile::Open("/home/rdehghani/Thesis/BiPo_studies/Rebin_BiPo214_May_2p2_optics_1000mm_16HL.root");
        //MC histogram file
        TFile *file2 = TFile::Open("/home/rdehghani/Thesis/BiPo_studies/Rebin_MC_BiPo214_May_2p2_optics_1000mm_16HL.root");
        //Response function histogram file(the delta-function gotten from the deconvolution script)
        TFile *file3 = TFile::Open("/home/rdehghani/Thesis/BiPo_studies/Decovoluted_MC_Data_BiPo214_May_2p2_optics_1000mm_16HL.root");
  	// Path to the custom function (test function)
	//----------------------------------------------------------------------------------------
	//declare the type of these histograms
	TH1D *datahist, *MChist;
        TH1F *RF_signal, *RF_test;
        cout << "You passed the first barrier; files are open" << endl;

        datahist = (TH1D*)file1 -> Get("Po214_Energy");
        MChist = (TH1D*)file2 -> Get("Po214_Energy"); 
	
	TH1D *convolved_hist = new TH1D("Response_MC", "Convolution_hist", 500, 0., 5);
	convolved_hist ->GetXaxis()->SetTitle("Energy [MeV]");
        //----------------------------------------------------------------------------------------
	Int_t i;
	Int_t nbins = 500;
	float MC[nbins];
	float data[nbins];//the raw po214 data
	float ResponseF[nbins];// response function
	float C[i];//the convolved reasulting signal

	for(i = 0; i<nbins; i++) data[i] = datahist->GetBinContent(i+1);
	for(i = 0; i<nbins; i++) MC[i] = MChist->GetBinContent(i+1);
        //----------------------------------------------------------------------------------------
	//define the delta function and the exponential function:
	cout << "You're on track" << endl;	
	double A = 0.0002;
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
        //----------------------------------------------------------------------------------------
	//Ploting the convolved signals
        TH1D* Raw_data_signal = new TH1D("Po214_raw_data", "Po214_raw_data", 500, 0, 5);

	//Perform the convolution:
	int m = MChist->GetNbinsX();
	cout << m << endl;
    	int n = nbins;
	cout << n << endl;

	//for (int i=0 ; i < n+m-1; i++)
	for (int i=0 ; i < m; i++)//convolved histogram loop
       	{
		double C = 0;
		for (int j= 0 ; j <= i ; j++)//MC loop, since "j" is the MC index
		{
		if (i-j >= 0) //check to omly get the possitive side of the response function
		{
			C += MC[j] * ResponseF[i - j];//getting the parameters and adding to the previous one 
		}
		}
		convolved_hist->SetBinContent(i+1 , C);
	}//the end of the Perform convolutoin loop
	// Fill raw data signal histogram
   	 for (int i = 0; i < 500; i++) {
       		Raw_data_signal->SetBinContent(i+1, data[i]);
	 }
        //----------------------------------------------------------------------------------------	
	// In this part, we normalize both our signals by the area under the curves so the number of data point don't affect our plots
	Raw_data_signal -> Scale(1/Raw_data_signal->Integral());
	convolved_hist-> Scale(1/convolved_hist->Integral());
	
	//Fitting the convolved_hist (which we hope to look like a data signal) to compare with the actual data-gauss-fit parameters	

	Raw_data_signal -> SetLineColor(kRed);
	convolved_hist-> SetLineColor(kBlack);
	
	TCanvas *c4 = new TCanvas("c4", "c4", 800, 600);
	c4 -> SetTickx();
	c4 -> SetTicky();
	c4 -> SetGridx();
	c4 -> SetGridy();
	convolved_hist -> Draw("HIST");
	Raw_data_signal -> Draw("HIST same");
		
}


