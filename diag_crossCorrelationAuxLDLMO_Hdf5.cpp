/*
This is a compilable program to run cross-correlation analysis on multiple channels using data from QRaw_RDCF files.
Make sure that you have a working version of diana sourced before running this program.
When compiling, compile with g++ -std=c++0x -D_HAVE_FFTW3_ `root-config --cflags --glibs` -lRooFit `diana-config --cflags --libs`


*/

#include "H5Cpp.h"
#include <hdf5.h>
#include <hdf5_hl.h>
#include <H5Cpp.h>
//#include "/usr/include/hdf5_hl.h"


#include <cstdlib>
#include <complex>
#include <cmath>
#include <iostream>
#include <string>
#include <getopt.h>
#include <sstream>

#include <TPaveText.h>
#include <TString.h>
#include <TMath.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TAxis.h>
#include <TRandom.h>
#include <TMarker.h>
#include <TH2F.h>
#include <TSystem.h>
#include <TStyle.h>

#include "QCuore.hh"
#include "QComplex.hh"
#include "QVector.hh"
#include "QVectorC.hh"
#include "QMatrix.hh"
#include "QMatrixC.hh"
#include "QRdcfRootFileReader.hh"
#include "QRdcfRootFileReaderHandler.hh"
#include "QRealComplexFFTW3.hh"
#include "QFFT.hh"
#include "QGlobalHandle.hh"
#include "QChain.hh"
#include "QRawWaveformGetter.hh"
#include "QDetTrgParamsDerivativeHandle.hh"
#include "QDetTrgParamsDerivative.hh"
#include "QGlobalDataManager.hh"
#include "QDetChannelHandle.hh"

#include <gsl/gsl_fit.h>
using namespace Cuore;


int SearchForTrigger(QVector timestream,double Threshold, double AvgWin, double Debounce, double adc2mV) {

	// threshold in derParams is in mV per ms
  // threshold is Multiplied by AvgWin to speed up calculation
  // Getting the SAME th normalization as in QApolloDerivativeTrigger.cc
  Threshold *= (AvgWin/adc2mV); //Threshold avg slope for determining a pulse (in units of ADC * Sampling Frequency)
	Debounce *= 2; 	
	AvgWin *= 2; 	
	
  // Number of samples above threshold
  int fCount = 0;

  // Loop from start to end in derivative
  unsigned int bLow, bUp;
  int k=0;
	int TrigNumber = 0;
	int DeadTime = 100*2;
	while(k < (int) timestream.Size()-AvgWin-DeadTime){
    bLow = k;
    bUp = k + AvgWin;

    if(bLow >= timestream.Size() || bUp >= timestream.Size()) return 0;
    // Evaluate derivative jump in a wider window (width=AvgWin)
    int low = timestream[bLow];
    int up = timestream[bUp];
    int delta = up - low;

    // If delta exceeds threshold, check for triggers
    if (delta > Threshold) {
      fCount++;
      // If above threshold for enough points (Debounce), fill trigger
      if (fCount == Debounce) {
        TrigNumber++;
        fCount = 0;
        k += DeadTime;
      }
    } else {
      // Not above threshold, reset count
      fCount = 0;
    }
		k++;
  }
  return TrigNumber;
}




bool hasPulse(QVector timestream,double thresh, double average, double debounce, double adc2mV,double MaxValue_MaxInRMS = 4.5){
    //Call the trigger func
    int ntrigs = SearchForTrigger(timestream, thresh,  average,  debounce, adc2mV);;
    
    //first: maxmin in window in RMS cut
    if(ntrigs == 0){
      if(timestream.GetMax()/timestream.GetRMS()>MaxValue_MaxInRMS){ //4.5
        ntrigs+=1;
        std::cout<<"\tMAXMININWININRMS"<<std::endl;
      }
    }
    

    return ntrigs>0;
}


void GetSamples(const QRdcfRootFileReader* reader, const int ch, const Long64_t &startns, const Long64_t &endns, QVector& wave){
    QError err = QERR_SUCCESS;
    err = reader->GetSamples(wave, ch, startns, endns);

    if (err == QERR_OUT_OF_RANGE){
        return;
    }
    return;
}

void GetTimestream(const QRdcfRootFileReader* reader, 
                   const int &ch, const double &tstart, const double &tend, QVector& wave,
                   int run,
                   int DownsamplingFactor=1){
  Long64_t startns = (Long64_t) (std::round(tstart*1.0e9));
  Long64_t endns = (Long64_t) (std::round(tend*1.0e9));

	QVector temp;
  GetSamples(reader, ch, startns, endns, temp);

  //removing slope
  double fBaselineIntercept,fBaselineSlope;
  double cov00, cov01, cov11,fBaselineRMS;

  // interpolate the baseline to remove it
  size_t npoints = temp.Size();
  Double_t* TimeArray = new Double_t[npoints];
  for(size_t i = 0; i < npoints; i++){
    TimeArray[i] = (Double_t)i;
  }
  gsl_fit_linear (TimeArray, 1, temp.GetArray(), 1, npoints, 
          &fBaselineIntercept, &fBaselineSlope, &cov00, &cov01, &cov11, &fBaselineRMS);
    

  //downsampling to the expected size
  wave.Resize(int(npoints/DownsamplingFactor));
  wave.Initialize(0);
  for(size_t i = 0; i < npoints; i+=DownsamplingFactor){
    wave[uint(i/DownsamplingFactor)] = (temp[i]-TimeArray[i]*fBaselineSlope-fBaselineIntercept);///GainByChannel[ch];
  }


  //downsampling to 2kHz for LD2LD
  // if(run<500936){
  //   wave.Resize(int(npoints/2));
  //   wave.Initialize(0);
  //   for(size_t i = 0; i < npoints; i+=2){
  //     wave[uint(i/2)] = (temp[i]-TimeArray[i]*fBaselineSlope-fBaselineIntercept);///GainByChannel[ch];
  //   }
  // }
  // else{
  //   wave = temp;
  //   }
  // }

  

}

void GetFFT(const QVector& wave, QVectorC&fft, QRealComplexFFTW3& transformer, const bool trunc=false) {
    transformer.TransformToFreq(wave, fft, trunc);
}

void plotTimestream(QVector wave, int run, int ch, int t0, double fSampleFreq, TString output, TString prefix){
	std::cout<<"Saving the wave"<<std::endl;
    TFile* newfile = new TFile(Form("%s/timestreams/%s_run%d_ch%04d_t%d.root",output.Data(),prefix.Data(), run, ch,t0),"RECREATE");
    TCanvas *can_ts = new TCanvas("can_ts", "timestream", 1400, 1000);

    double tlength = wave.Size()/fSampleFreq;
    can_ts->cd();
    TGraph *gWave = new TGraph();
    gWave= wave.GetGraph(fSampleFreq);

	//adding the differential one too
	TGraph *gWaveDiff = new TGraph(wave.Size()-1);
	for(size_t kk = 0; kk<wave.Size()-1; kk++){
		gWaveDiff->SetPoint(kk,kk/fSampleFreq,abs(wave[kk+1]-wave[kk]));
	}
		gWave->SetTitle("Timestream; time [s]; Voltage [a.u.]");
    gWave->SetLineColor(kBlack);
    gWave->SetLineWidth(2);
    gWave->Draw("AL");
    gWave->GetXaxis()->SetRangeUser(0.0,tlength);

    gWaveDiff->SetLineColor(kBlue);
    gWaveDiff->SetLineWidth(2);
    gWaveDiff->Draw("LSAME");
    gWaveDiff->GetXaxis()->SetRangeUser(0.0,tlength);

    gPad->SetGridx();
    gPad->SetGridy();
    can_ts->Update();

    TPaveText *pt = new TPaveText(0.8, 0.8, 0.99 ,0.99,"brNDC");
    pt->AddText(Form("Run = %d", run));
    pt->AddText(Form("ch = %d", ch));
    pt->AddText(Form ("t_start = %d", t0));
    pt->Draw();

    newfile->cd();
    can_ts->Write("can_ts");
    can_ts->SaveAs(Form("%s/timestreamspng/%s_run%d_ch%04d_t%d.png",output.Data(),prefix.Data(),run,ch,t0),"RECREATE");
    newfile->Close();
    delete can_ts;
    delete pt;
    delete newfile;
    delete gWave;
    delete gWaveDiff;
    return;
}

void plotSinglePSD(QVector psd, int run, int ch, double fNyqFreq, TString output, TString prefix){
    // This plots a single power spectral denisty

    TCanvas *can_psd = new TCanvas("can_psd", "psd", 1000, 1000);
    can_psd->cd();
    //this is the "sampling time," like the sampling frequency but in Fourier space
    //it is (1.0/df), where df is the size of each frequency bin
    double df = fNyqFreq/psd.Size();
    double tsamp = 1.0/df;
    TGraph *gpsd = psd.GetGraph(tsamp);
    gpsd->SetTitle("Power Spectral Density; Frequency [Hz]; PSD [mV^2/Hz]");
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetLogy();
    gPad->SetLogx();
    gpsd->GetXaxis()->SetRangeUser(5e-2,fNyqFreq);
    //gpsd->GetYaxis()->SetRangeUser(1e-7,1e3);
    gpsd->SetLineColor(kBlack);
    gpsd->SetLineWidth(2);
    gpsd->Draw("AL");

    TPaveText *pt = new TPaveText(0.8, 0.8, 0.99 ,0.99,"brNDC");
    pt->AddText(Form("Run = %d", run));
    pt->AddText(Form("ch = %d", ch));
    pt->Draw();
    can_psd->Update();

    TFile* newfile = new TFile(Form("%s/psd/%s_run%d_ch%04d.root",output.Data(),prefix.Data(),run, ch),"RECREATE");
    newfile->cd();
    can_psd->Write("can_psd");
    can_psd->SaveAs(Form("%s/psd/%s_run%d_ch%04d.png",output.Data(), prefix.Data(), run, ch),"RECREATE");
    newfile->Close();
    delete can_psd;
    delete pt;
    delete newfile;
    return;
}
//BEGIN NEW FUNCTION FROM TRISTAN
//I wanted to make a function that just outputs the diagonal h5 file(1MB/RunHr), not the full matrix h5 file(5GB/RunHr)
//copying the structure of the WriteMats2H5 function directly below, but not needing the full matrix
void WriteDiagMats2H5(H5::Group* G_D, H5::DataSpace* D_D, 
          const QMatrixC &mat, QVector boloANPS, QVector auxANPSvec,double SamplingFrequency){
  int nxbins = mat.GetNRow();
  int nybins = mat.GetNCol();
  if(nxbins != nybins){
    std::cout << "ERROR: MATRIX NOT SQUARE" << std::endl;
    return;
  }
	//from matrix to array
	std::cout<<"Transcribing to arrays"<<std::endl;
	std::cout<<"Matrix size = "<<nxbins<<"*"<<nybins<<" = "<<nxbins*nybins<<std::endl;
	//size_t dimtot = size_t(nxbins*nybins);
	double diagM[nxbins];
	double diagP[nxbins];

	hsize_t  cdims[1];
	//cdims[0] = 1;
  cdims[0] = nybins;//1000;
	H5::DSetCreatPropList ds_creatplist;  // create dataset creation prop list
  ds_creatplist.setChunk( 1, cdims );  // then modify it for compression
  ds_creatplist.setDeflate( 9 );

  for(int i = 0; i<nxbins;i++){
		diagM[i] = mat(i,i).GetMagnitude()/std::sqrt(boloANPS(i)*auxANPSvec(i));
		diagP[i] = mat(i,i).GetPhase();    
  }

  
	std::cout<<"Write Diags"<<std::endl;
	G_D->createDataSet("Magnitude",H5::PredType::NATIVE_DOUBLE,*D_D).write(diagM,H5::PredType::NATIVE_DOUBLE);
	G_D->createDataSet("Phase",H5::PredType::NATIVE_DOUBLE,*D_D).write(diagP,H5::PredType::NATIVE_DOUBLE);
  
  H5::DataSpace dspace(H5S_SCALAR);
  G_D->createAttribute("SamplingFrequency",H5::PredType::NATIVE_DOUBLE,dspace).write(H5::PredType::NATIVE_DOUBLE,&SamplingFrequency);	

}
//END NEW FUNCTION FROM TRISTAN

void WriteMats2H5(H5::Group* G_M,H5::Group* G_D,
				  H5::DataSpace* D_M, H5::DataSpace* D_D,
				  QMatrixC mat,QVector boloANPS,QVector auxANPSvec,double SamplingFrequency)//, int run, int ch, int chaux, double fNyqFreq, TString output)
{
    int nxbins = mat.GetNRow();
    int nybins = mat.GetNCol();
    if(nxbins != nybins){
        std::cout << "ERROR: MATRIX NOT SQUARE" << std::endl;
        return;
    }
	//from matrix to array
	std::cout<<"Transcribing to arrays"<<std::endl;
	std::cout<<"Matrix size = "<<nxbins<<"*"<<nybins<<" = "<<nxbins*nybins<<std::endl;
	size_t dimtot = size_t(nxbins*nybins);
	double* MagVersRow = new double[dimtot];
	double* PhaseVersRow = new double[dimtot];
	double diagM[nxbins];
	double diagP[nxbins];

	hsize_t  cdims[1];
	//cdims[0] = 1;
  cdims[0] = nybins;//1000;
	H5::DSetCreatPropList ds_creatplist;  // create dataset creation prop list
  ds_creatplist.setChunk( 1, cdims );  // then modify it for compression
  ds_creatplist.setDeflate( 9 );
	
	for (int i = 0; i < nxbins; i++){
        for (int bin = 0; bin < nybins; bin++)
		{
//			std::cout<<"nx = "<<i<<" , ny = "<<bin<<" --> index = "<<i*(nxbins)+bin<<std::endl;
			MagVersRow[i*(nxbins)+bin] = double(mat(i,bin).GetMagnitude()/std::sqrt(boloANPS(i)*auxANPSvec(bin))); 
			PhaseVersRow[i*(nxbins)+bin] = double(mat(i,bin).GetPhase()); 
			//MagVersRow[i][bin] = double(mat(i,bin).GetMagnitude()/std::pow(boloANPS(i)*auxANPSvec(bin),0.5)); 
			//PhaseVersRow[i][bin] = double(mat(i,bin).GetPhase()); 
		}
	
		
		diagM[i] = mat(i,i).GetMagnitude()/std::sqrt(boloANPS(i)*auxANPSvec(i));
		diagP[i] = mat(i,i).GetPhase();
	
	}
	std::ostringstream oss;
	std::cout<<"Write Magn"<<std::endl;
	oss<<"Magnitude";//<<i;
	G_M->createDataSet(oss.str().c_str(),H5::PredType::NATIVE_DOUBLE,*D_M,ds_creatplist).write(MagVersRow,H5::PredType::NATIVE_DOUBLE);	
	oss.str("");
	oss.clear();

	std::cout<<"Write Phase"<<std::endl;
	oss<<"Phase";//<<i;
	G_M->createDataSet(oss.str().c_str(),H5::PredType::NATIVE_DOUBLE,*D_M,ds_creatplist).write(PhaseVersRow,H5::PredType::NATIVE_DOUBLE);	
	oss.str("");
	oss.clear();

	delete[] MagVersRow;
	delete[] PhaseVersRow;
	
	//to HDF5
	
	//std::cout<<"Transcribing matrix to HDf5"<<std::endl;
//	for (int i = 0; i < nxbins; i++){
//		std::cout<<"\trow "<<i<<std::endl;
	//	G_M->createDataSet("Magnitude",H5::PredType::NATIVE_DOUBLE,*D_M).write(MagVersRow,H5::PredType::NATIVE_DOUBLE);
	//	G_M->createDataSet("Phase",H5::PredType::NATIVE_DOUBLE,*D_M).write(PhaseVersRow,H5::PredType::NATIVE_DOUBLE);	
//	}
	
	std::cout<<"Write Diags"<<std::endl;
	G_D->createDataSet("Magnitude",H5::PredType::NATIVE_DOUBLE,*D_D).write(diagM,H5::PredType::NATIVE_DOUBLE);
	G_D->createDataSet("Phase",H5::PredType::NATIVE_DOUBLE,*D_D).write(diagP,H5::PredType::NATIVE_DOUBLE);
  
  H5::DataSpace dspace(H5S_SCALAR);
  G_D->createAttribute("SamplingFrequency",H5::PredType::NATIVE_DOUBLE,dspace).write(H5::PredType::NATIVE_DOUBLE,&SamplingFrequency);	

}

void plotCorrelationMatrix(QMatrixC mat,QVector boloANPS,QVector auxANPSvec, int run, int ch, int chaux, double fNyqFreq, TString output){

    int nxbins = mat.GetNRow();
    int nybins = mat.GetNCol();
    if(nxbins != nybins){
        std::cout << "ERROR: MATRIX NOT SQUARE" << std::endl;
        return;
    }
    //int nfreqbins = nxbins;

    //The amphist contains the amplitudes of the cross-correlation matrix
    //TH2F* amphist  = new TH2F(Form("CorrMag_ch%04d_aux%04d",ch,chaux),Form("Frequency Correlation: ch%04d vs. aux%04d; Freq ch%04d; Freq aux%04d",ch, chaux, ch, chaux),nfreqbins,0,fNyqFreq,nfreqbins,0,fNyqFreq);
    
    //The phasehist contains all of the phases at each frequency
    //TH2F* phasehist = new TH2F(Form("CorrPhase_ch%04d_aux%04d",ch,chaux),Form("Frequency Correlation Phase: ch%04d vs. aux%04d; Freq ch%04d; Freq aux%04d",ch, chaux, ch, chaux),nfreqbins,0,fNyqFreq,nfreqbins,0,fNyqFreq);
	
	//berettam: putting also a txt dump for python
	std::ofstream outM(Form("%s/txtfiles/CrossCorrelationMag_run%d_ch%04d_aux%04d.txt",output.Data(),run,ch,chaux),std::ios::out);
	std::ofstream outP(Form("%s/txtfiles/CrossCorrelationPhase_run%d_ch%04d_aux%04d.txt",output.Data(),run,ch,chaux),std::ios::out);
	std::ofstream outMDiag(Form("%s/txtfiles/CrossCorrelationDiag_run%d_ch%04d_aux%04d.txt",output.Data(),run,ch,chaux),std::ios::out);
	std::ofstream outMDiagP(Form("%s/txtfiles/CrossCorrelationDiagPhase_run%d_ch%04d_aux%04d.txt",output.Data(),run,ch,chaux),std::ios::out);
	
    for (int i = 0; i < nxbins; i++){
        QVectorC row0 = mat.GetRow(i);
        for (int bin = 0; bin < nybins; bin++){ 
         //   amphist->SetBinContent(i,bin, row0(bin).GetMagnitude());   
			outM << row0(bin).GetMagnitude()/std::pow(boloANPS(i)*auxANPSvec(bin),0.5) << ","; 
			outP << row0(bin).GetPhase() << ","; 
      //      phasehist->SetBinContent(i, bin, row0(bin).GetPhase());
        }
		outMDiag << row0(i).GetMagnitude()/std::pow(boloANPS(i)*auxANPSvec(i),0.5)<<std::endl;
		outMDiagP << row0(i).GetPhase()<<std::endl;
		outM << std::endl;
		outP << std::endl;
    }
	outMDiag.close();
	outMDiagP.close();
	outM.close();
	outP.close();

//    TCanvas* c1 = new TCanvas("c1","c1", 3200,2000);

 //   gStyle->SetOptStat(0);
  //  amphist->Draw("COLZ"); //phasehist->Draw("COLZ");

   // c1->SetLogz();
   // gPad->SetGridx();
   // gPad->SetGridy();
   // c1->Update();
//
//    TFile* newfile = new TFile(Form("%s/CrossCorrelation_run%d_ch%04d_aux%04d.root",output.Data(),run,ch,chaux),"RECREATE");
//    newfile->cd();
//    c1->Write("c1");
    //c1->SaveAs(Form("%s/CrossCorrelation_run%d_ch%04d_aux%04d.png",output.Data(),run,ch,chaux),"RECREATE");
//    newfile->Close();
   // delete c1;
    return;
}

void parseListToValues(const std::string& listOfValues, std::vector<std::string>& parsedList) {
    if (listOfValues != "") {
        std::stringstream ss(listOfValues);
        while (ss.good()) {
            std::string substr;
            getline(ss, substr, ' ');
            parsedList.push_back(substr);
        }
    }
}



void Usage()
{

  std::cout << std::endl << std::endl;
  std::cout << "---------------------------------------------------------------------------------" << std::endl;
  if( std::getenv("USER") != NULL )
    std::cout << "Hi " << std::getenv("USER") <<"! The usage of this wonderful program is: " << std::endl;
  else
    std::cout << "Hi there! The usage of this wonderful program is:" << std::endl;
  std::cout<<"./CrossCorrLDLMOH5"<<std::endl;
  std::cout<<"Run with the following options "<<std::endl;
  std::cout<<"options "<<std::endl;
  std::cout<<"------------------------------------------------------------"<<std::endl;
  std::cout<<"-r (or --run)             run to be processed                                 [REQUIRED] "<<std::endl;
  std::cout<<"-c (or --cahnnel)         channel to be analyzed                              [REQUIRED] "<<std::endl;
  std::cout<<"-t (or --tstart)          start time for scan                                 [REQUIRED]"<<std::endl;
  std::cout<<"-T (or --tstop)           stop time for the scan                              [REQUIRED]"<<std::endl;
  std::cout<<"-w (or --twindow)         window length to do the scan                        [REQUIRED]"<<std::endl;
  std::cout<<"-o (or --outputpath)      path where to save the output                       [REQUIRED]"<<std::endl;
  std::cout<<"-f (or --trigparamsfile)  file for the trg params. If not provided, using DB  [empty]"<<std::endl; 
  std::cout<<"-g (or --gainfile)        file for the CH gains. If not provided, using 1     [empty]"<<std::endl; 
  std::cout<<"-a (or --adc2mv)          window length to do the scan                        [8.02e-2 mV/ADC]"<<std::endl;
  std::cout<<"-s (or --samplingfreq)    path where to save the output                       [1000 Hz]"<<std::endl;
  std::cout<<"-A (or --auxCH)           comma-separated list of aux channels                [1,2,4,5,6,11,12,14,15,16,17,18]"<<std::endl;
  std::cout<<"-C (or --allCH)           comma-separated list of all channels                [1,2,4,5,6,11,12,14,15,16,17,18]"<<std::endl;
  std::cout<<"-h (or --help)            prints this help                                    [no opt]"<<std::endl;
  std::cout<<std::endl;  
  std::cout<<"Compile with:"<<std::endl;
  std::cout<<"g++ -o CrossCorrLDLMOH5 -std=c++0x -I/mnt/disk1/software/hdf5/install/hdf5-1.12.0/include -L/mnt/disk1/software/hdf5/install/hdf5-1.12.0/lib -lhdf5_cpp -lhdf5_hl_cpp -lhdf5 -lhdf5_hl -D_HAVE_FFTW3_ `root-config --cflags --glibs` -lRooFit `diana-config --cflags --libs` -O3 crossCorrelationAuxLDLMO_Hdf5.cpp"<<std::endl;
}


std::vector<std::string> getWords(std::string s,std::string delimiter){
    std::vector<std::string> res;
    int pos = 0;
    
    while(pos != -1){
      pos = s.find(",");
      res.push_back(s.substr(0,pos));
      s.erase(0,pos+delimiter.length());
    }
    return res;
}


#ifndef __CINT__
int main(int argc, char** argv){


  
  if(argc<2 ){
    std::cout<< "No arguments provided!!!" << std::endl;
    Usage();
    return 1;
  }

  int run = 0;
  int ch = 0;
  int tstart = 0;
  int tend = 0;
  int twindow = 0;
  // TString output = "";
  bool InputTrigFile = false;
  bool InputGainFile = false;
  
  double adc2mV = 8.0108642578125e-02;
  double fSampleFreq = 1000; //berettam: putting a fixed value, should change according to channel FIXME
	
  std::string outpath = "";
  std::string trigfile = "";
  
  std::string gainfile = "";
  std::string test = "";
  
  std::string auxchannelsList;// = {1,2,4,5,6,11,12,14,15,16,17,18};
  std::string AllchannelsList;// = {1,2,4,5,6,11,12,14,15,16,17,18};
  
  // Parse the options
  static struct option long_options[] = {
                                        { "run",            required_argument,  0,  'r'},
                                        { "channel",        required_argument,  0,  'c'},
                                        { "tstart",         required_argument,  0,  't'},
                                        { "tstop",          required_argument,  0,  'T'},
                                        { "twindow",        required_argument,  0,  'w'},
                                        { "outpath",        required_argument,  0,  'o'},
                                        { "trigparamsfile", optional_argument,  0,  'f'},
                                        { "gainfile",       optional_argument,  0,  'g'},
                                        { "adc2mv",         optional_argument,  0,  'a'},
                                        { "samplingfreq",   optional_argument,  0,  's'},
                                        { "auxCH",          optional_argument,  0,  'A'},
                                        { "allCH",          optional_argument,  0,  'C'},
                                        { "help",           no_argument,        0,  'h'},
                                        {0, 0, 0, 0}
  };
  const char* const short_options = "r:c:t:T:w:o:f:g:a::s::A:C:h";
  int c;

  
  while ((c = getopt_long(argc, argv, short_options, long_options, nullptr)) != -1 ) {
    
    switch (c) {
      case 'r': {
        run = atoi(optarg);
        break;
      }
      case 'c':{
        ch=atoi(optarg);
        break;
      }
      case 't':{
        tstart= atoi(optarg);
        break;
      }
      case 'T':{
        tend= atoi(optarg);
        break;
      }
      case 'w':{
        twindow= atoi(optarg);
        break;
      }
      case 'o':{
        outpath=optarg;
        std::cout<<"\n\n\noutdir = "<<outpath<<"\n\n\n";
        break;
      }
      case 'f':{
        if(optarg == NULL){
          InputTrigFile = false;
        }
        else{
          std::cout<<"QUA!!"<<std::endl;
          InputTrigFile = true;
          trigfile = optarg;
        }
        break;
      }    
      case 'g':{
        if(optarg == NULL){
          InputGainFile = false;
          std::cout<<"QUAN!!"<<std::endl;
        }
        else{
          std::cout<<"QUA!!"<<std::endl;
          InputGainFile = true;
          gainfile = optarg;
        }
        break;
      }    
      case 'a':{
        adc2mV= (optarg==NULL) ?  8.0108642578125e-02 : atof(optarg);
        break;
      }
      case 's':{
        if(optarg) {fSampleFreq = atoi(optarg);}
        //fSampleFreq= (optarg==NULL) ? 1000 : atof(optarg);
        //std::cout << "fSampleFreq is : " << fSampleFreq << std::endl;
        break;
      }      
      case 'A':{
        auxchannelsList = optarg;
        break;
      }
      case 'C':{
        AllchannelsList = optarg;
        break;
      }
      case 'h':{
        Usage();
        exit(1);
      }
      default: {
        std::cout<<"Unknown argument!!!!"<<std::endl;
        exit(1);
      }
    }
  }
  //Creating the output filename
  TString output = Form("%s/run%d", outpath.c_str(),run);




    /********************************THESE MAY NEED TO BE RE-DEFINED FOR CUPID BDPT TESTS********************************/
    //double fSampleFreq_LD = 2000; 

    // fSampleFreq_LMO = 1000;
    
    //berettam: building the channel list
     std::vector<int> auxchannels;
     std::vector<int> Allchannels;

    //fill the list from the input
    if (auxchannelsList.empty()){
      auxchannels = {1,2,4,5,6,11,12,14,15,16,17,18};
    }
    else {
      std::string delimiter = ",";
      std::vector<std::string> parsing = getWords(auxchannelsList,delimiter);
      for(size_t j = 0;j<parsing.size();j++){
        auxchannels.push_back(atoi(parsing[j].c_str()));
      }
    }

    
    if (AllchannelsList.empty()){
      Allchannels = {1,2,4,5,6,11,12,14,15,16,17,18};
    }
    else {
      std::string delimiter = ",";
      std::vector<std::string> parsing = getWords(AllchannelsList,delimiter);
      for(size_t j = 0;j<parsing.size();j++){
        Allchannels.push_back(atoi(parsing[j].c_str()));
      }
    }

    std::cout<<std::endl;
    std::cout<<"AuxChannels Vector"<<std::endl;
    for(size_t j = 0;j<auxchannels.size();j++) std::cout<<auxchannels[j]<<std::endl;
    std::cout<<std::endl;

    std::cout<<std::endl;
    std::cout<<"AllChannels Vector"<<std::endl;
    for(size_t j = 0;j<Allchannels.size();j++) std::cout<<Allchannels[j]<<std::endl;
    std::cout<<std::endl;

    //check that the main ch is in both lists, if not add it!
    if(std::find(auxchannels.begin(), auxchannels.end(), ch) == auxchannels.end()) {
      /* v does not contain x */
      auxchannels.push_back(ch);
    }
    if(std::find(Allchannels.begin(), Allchannels.end(), ch) == Allchannels.end()) {
      /* v does not contain x */
      Allchannels.push_back(ch);
    }

    std::map< int,std::vector<int> > SideMap;
    SideMap[18].push_back(6);
    SideMap[17].push_back(6);
    SideMap[17].push_back(5);
    SideMap[16].push_back(5);
    SideMap[16].push_back(4);
    SideMap[15].push_back(4);
    SideMap[12].push_back(2);
    SideMap[12].push_back(1);
    SideMap[11].push_back(1);
	
    SideMap[1].push_back(11);
    SideMap[1].push_back(12);
    SideMap[2].push_back(12);
    SideMap[4].push_back(15);
    SideMap[4].push_back(16);
    SideMap[5].push_back(17);
    SideMap[5].push_back(16);
    SideMap[6].push_back(17);
    SideMap[6].push_back(18);
	

	/********************************************************************************************************************/


  /********************************************************************************************************************/
  //READING the DAQ parameters form DB	
  std::cout<<std::endl<<"Getting DAQ parameters from the DB"<<std::endl;
  QGlobalDataManager dm;
  std::map<int,double> adc2mV_byCH;
  std::map<int,double> SF_byCH;  
  std::map<int,int> DownsamplingFactor_byCH;

  for(size_t i =0; i<Allchannels.size(); i++){
    
    QDetChannelHandle handle_DAQ(Allchannels[i]);
    handle_DAQ.SetRun(run);
    dm.Get(&handle_DAQ,"DB");
    QDetChannel detChannel = handle_DAQ.Get();
    adc2mV_byCH[Allchannels[i]] = 1000.*( detChannel.fDaqSet.fVAdcLimitMax - detChannel.fDaqSet.fVAdcLimitMin ) / pow(2,(double)( detChannel.fDig.fDaqNBits ));
    SF_byCH[Allchannels[i]] = (double) detChannel.fDaqSet.fSamplingFrequency; 
    DownsamplingFactor_byCH[Allchannels[i]] = int(SF_byCH[Allchannels[i]]/fSampleFreq);
  
    std::cout<<Form("CH %02d \t SF = %.0f \t adc2mV = %.2e \t DownSampling = %d",Allchannels[i],SF_byCH[Allchannels[i]],adc2mV_byCH[Allchannels[i]],DownsamplingFactor_byCH[Allchannels[i]])<<std::endl;
	}
  
  std::cout<<"Finished Getting DAQ parameters from the DB"<<std::endl<<std::endl<<std::endl;
  /********************************************************************************************************************/



	//declaring the reading structure
  QRdcfRootFileReaderHandler& handler  = QRdcfRootFileReaderHandler::GetInstance();
  QError err = QERR_SUCCESS;
	//map to keep the readers
	std::map< int , const QRdcfRootFileReader* > readersCH;
	for(size_t i =0; i<Allchannels.size(); i++){
		readersCH[Allchannels[i]] = handler.GetReader(run, Allchannels[i], err);
		if(!readersCH[Allchannels[i]] || err != QERR_SUCCESS){
			std::cout << err << std::endl;
			std::cout << "No file found for ch: "<< ch << std::endl;
			return 231192;
		}
	}
	//check if the original channel is in the vector
	if (!std::count(Allchannels.begin(), Allchannels.end(), ch)) {
    std::cout << "Original channel is not in list, adding the reader";
		readersCH[ch] = handler.GetReader(run, ch, err);
		if(!Allchannels[ch] || err != QERR_SUCCESS){
			std::cout << err << std::endl;
			std::cout << "No file found for ch: "<< ch << std::endl;
			return 231192;
		}
	}


  //handy frequency-based parameters
  double fNyqFreq = fSampleFreq / 2.0;
  double freqBinWidth = (double) 1.0/twindow;
  int nfreqbins = fNyqFreq/freqBinWidth+1;

  //Create the QVector(C) instances we will need
  QVector ts(0);
  QVector tside(0);
  QVector ts2(0);
  QVector tsaux(0);
  QVectorC fft1(0);
  QVectorC auxfft(0);

  //Get the trigger thresholds for this channel-runs from the database
  std::map< int,double > threshmVCH;
  std::map< int,double > averageCH;
  std::map< int,double > debounceCH;
  if(InputTrigFile){
    //open the file
    std::cout<<"Reading trigger from file: "<<trigfile<<std::endl;
    std::cout<<"assuming this format for each line: "<<std::endl;
    std::cout<<"ch average threshold debounce"<<std::endl;
    std::string line;
    ifstream trigfilename(trigfile) ;
    while (getline(trigfilename,line))
    {
      std::vector<std::string> vec;
      //std::cout<<line<<std::endl;
      parseListToValues(line,vec);
      averageCH[stoi(vec[0])] = stof(vec[1]);
      threshmVCH[stoi(vec[0])] = stof(vec[2]);
      debounceCH[stoi(vec[0])] = stof(vec[3]);
    }
    trigfilename.close();
	}
	else{
		//trigger params from the DB	
		std::cout<<"Getting the trigger parameters from the DB"<<std::endl;
		// QGlobalDataManager dm;
		QDetTrgParamsDerivative *result = NULL;
		for(size_t i =0; i<auxchannels.size(); i++){
			
			QDetChannelHandle handle(auxchannels[i]); 
			handle.SetRun(run);
			dm.Get(&handle,"DB");
			QDetChannel detChannel = handle.Get();
			result = (QDetTrgParamsDerivative *)(detChannel.fTrgSet.fTriggers[0].GetTrgParams());
			threshmVCH[auxchannels[i]] =result->fThreshold;
			averageCH[auxchannels[i]] =result->fAverage;
			debounceCH[auxchannels[i]] =result->fDebounce;
		}
		if (!std::count(auxchannels.begin(), auxchannels.end(), ch)) {
			std::cout << "Original channel is not in list, adding the trigger params";
			QDetChannelHandle handle(ch);
			handle.SetRun(run);
			dm.Get(&handle,"DB");
			QDetChannel detChannel = handle.Get();
			result = (QDetTrgParamsDerivative *)(detChannel.fTrgSet.fTriggers[0].GetTrgParams());
			threshmVCH[ch] =result->fThreshold;
			averageCH[ch] =result->fAverage;
			debounceCH[ch] =result->fDebounce;
		}	
	}	
	//printing the channels parameters
	std::cout<<"Using the following parameters for triggering:"<<std::endl;
	std::cout<<"CH\tAvg\tTH\tDeb"<<std::endl;
	//if the params were read from the DB, write them in a txt file
	// so it is easier to customize them in case
	std::ofstream outtrigp ;
	if(!InputTrigFile){
		outtrigp.open(Form("TrigParams_Run%d.txt",run),std::ios::out);
	}
	for(size_t i =0; i<auxchannels.size(); i++){
		std::cout<<auxchannels[i]<<" "<<averageCH[auxchannels[i]]<<" "<<threshmVCH[auxchannels[i]]<<" "<<debounceCH[auxchannels[i]]<<std::endl;
		if(!InputTrigFile){
			outtrigp<<auxchannels[i]<<" "<<averageCH[auxchannels[i]]<<" "<<threshmVCH[auxchannels[i]]<<" "<<debounceCH[auxchannels[i]]<<std::endl;
		}
	}
	if(!InputTrigFile) outtrigp.close();
	std::cout<<"--> Read trigger Done"<<std::endl<<std::endl;

  /*
  * Get the gain for each channel
  */
	std::cout<<"--> Reading Gain"<<std::endl<<std::endl;
  
  std::map< int, double> GainByChannel;
  if(InputGainFile){
    //open the file
    std::cout<<"Reading Gain from file: "<<gainfile<<std::endl;
    std::cout<<"assuming this format for each line: "<<std::endl;
    std::cout<<"ch Gain"<<std::endl;
    std::string line;
    ifstream gainfilename(gainfile) ;
    while (getline(gainfilename,line))
    {
      std::vector<std::string> vec;
      //std::cout<<line<<std::endl;
      parseListToValues(line,vec);
      GainByChannel[stoi(vec[0])] = stof(vec[1]);
    }
    gainfilename.close();
  }
  else{
    for(size_t i =0; i<auxchannels.size(); i++){
      GainByChannel[auxchannels[i]] = 1.0;
    }
  }
  std::cout<<"Using the following gains:"<<std::endl;
  for(size_t i =0; i<auxchannels.size(); i++){
      std::cout<<auxchannels[i]<<" --> "<<GainByChannel[auxchannels[i]]<<std::endl;
  }

	std::cout<<"--> Reading Gain done"<<std::endl<<std::endl;


    //now create an instance of a QRealComplexFFTW3 which we can re-use when building all of the FFTs we need
    int N = twindow*fSampleFreq; // take the number of points of the LMO: will be the ultimate limit. LD will be downsampled
    QRealComplexFFTW3 fTransformer = QRealComplexFFTW3(N);
	  fTransformer.SetWindowType(QFFT::WT_Welch);


    //finally, some counters we need
    double t0 = tstart;
    

	//loop over the noise events
	//for the channel of interest
	//saving the ok windows w and their t0 in a map
	std::cout<<"Looping over time for ch "<<ch<<std::endl;
	std::map< double, QVector> VecTCH;
    while(t0 < tend){ 
        if((int)t0 % 3600 == 0){
            std::cout << "Reached time " << t0/3600 << " hours" << std::endl;
        }   
        std::cout << "Reached time " << t0 << " s" << std::endl;

        //Get timestream from the channel of interest
        GetTimestream(readersCH[ch], ch,t0,t0+twindow, ts,run,DownsamplingFactor_byCH[ch]);
//        std::cout<<ch<<") "<<t0<<" - "<<twindow<<" - "<<ts.Size()<<" - "<<N<<std::endl; 
		    if((int) ts.Size() != N){
            std::cout << "CH = " << ch << std::endl;
            std::cout << "Size did not match at t0 = " << t0 << std::endl;
            std::cout << Form("Expected Size = %d - Got size = %d",N,ts.Size()) << std::endl;
            std::cout << "Exiting the program " << std::endl;
            break; //EXIT ONCE THE TIMESTREAM IS NOT SUFFICIENTLY LONG. DO NOT DELETE
        }
		//plotTimestream(ts, run, ch, t0, fSampleFreq_LMO, output, "MainEvent");
              
		//look for pulse with derivative trigger
		//this uses a home-made derivative trigger to throw away pulse events (see the hasPulse function)
		//It is slow, but it works without needing a list of triggers as an input
    bool HasATrigger = hasPulse(ts,threshmVCH[ch],averageCH[ch],debounceCH[ch],adc2mV_byCH[ch]);
    if(HasATrigger){
      std::cout<<"Has Trigger"<<std::endl;
      // plotTimestream(ts, run, ch, t0, fSampleFreq_LMO, output, "triggeredEvent");
			t0 += twindow;
			continue;
		}
    else{ //check that the eventual sides do not have a trigger!
      std::cout<<"Side search "<<SideMap[ch].size()<<std::endl;
      if(SideMap[ch].size()>0){
        bool DoesSideTrigger = false;
        for (std::vector<int>::const_iterator it = SideMap[ch].begin(); it!=SideMap[ch].end() && !DoesSideTrigger; ++it) {
          std::cout<<" --> side = "<<*it<<std::endl;
          GetTimestream(readersCH[*it], *it,t0,t0+twindow, tside,run,DownsamplingFactor_byCH[*it]);
          if(*it<10){
            DoesSideTrigger = hasPulse(tside,threshmVCH[*it]*10,averageCH[*it],debounceCH[*it],adc2mV_byCH[*it],8);               
          }
          else {
            DoesSideTrigger = hasPulse(tside,threshmVCH[*it],averageCH[*it],debounceCH[*it],adc2mV_byCH[*it]);    
          }
          if(DoesSideTrigger) std::cout<<" --> TRIGGERED"<<std::endl;
        }     
        if(DoesSideTrigger) {
          t0 += twindow;
          continue;
        }   
      }
    }
		
		//saving in the map the good couple
		VecTCH[t0] = ts;
    std::cout<<"\nAdding pulse, total = "<<VecTCH.size()<<std::endl;
    //increasing the time
		t0+= twindow;
	}

	std::cout<<"-->Done"<<std::endl<<std::endl;

	//now doing the same for the aux channels. 
	//for each case I build the matrix to avoid the initial big vector of matrices
	//cycle on the aux detectors
	double adc2psd = adc2mV*adc2mV/((nfreqbins-1)*2)/fSampleFreq;
	double adc2fft = std::pow(adc2psd, 0.5);
	double count = 0; //to count the windows in the avgs, has to be double since it is used to rescale 
	//ch by ch correlation --> it is a vector now, I have a single channel
	std::vector<double> CHCHCorr;
	//tfile/ttree to save the vectors
	std::cout<<"Saving TTree here:" <<Form("%s/PulsesFileLDLMOs_BaseCH%d.root",output.Data(),ch)<<std::endl;
	TFile* fpulse = new TFile(Form("%s/rootfiles/PulsesFileLDLMOs_BaseCH%d.root",output.Data(),ch),"RECREATE");
	TTree* tvec = new TTree("tvec","Tree with the noise windows");
	TFile* fnoise = new TFile(Form("%s/rootfiles/NoiseFileLDLMOs_BaseCH%d.root",output.Data(),ch),"RECREATE");
	TTree* tnoise = new TTree("tnoise","Tree with the noise power spectra");
	std::vector<double> pulseref;
	std::vector<double> pulseaux ;
	int chref = ch;
	int chauxsave =0;//= 38;
	double startt =0;
	tvec->Branch("t0",&startt);
	tvec->Branch("chref",&chref);
	tvec->Branch("vref",&pulseref);
	tvec->Branch("vaux",&pulseaux);
	tvec->Branch("chaux",&chauxsave);

	//define hdf5 file
	TString FileNameH5_M = Form("%s/h5files/Matrices_BaseCH%d.h5",output.Data(),ch);
	TString FileNameH5_D = Form("%s/h5files/Diag_BaseCH%d.h5",output.Data(),ch);
	std::cout<<"Create HDF5 Files"<<std::endl;
	H5::H5File fileM (FileNameH5_M.Data(), H5F_ACC_TRUNC);
	H5::H5File fileD (FileNameH5_D.Data(), H5F_ACC_TRUNC);
	//H5::H5File *fileM = new H5::H5File(FileNameH5_M.Data(), H5F_ACC_TRUNC);
	//H5::H5File *fileD = new H5::H5File(FileNameH5_D.Data(), H5F_ACC_TRUNC);
	//save a dataset about the shape of the matrix
	hsize_t dimmatV[1];
	dimmatV[0] = 2;
	std::cout<<"\n\n\nCreate HDF5 Dataspacer"<<std::endl;
	H5::DataSpace dataspace_Shape(1, dimmatV);
	int shape[2] = {nfreqbins,nfreqbins};
	std::cout<<"\n\n\nwrite HDF5 DataSet"<<std::endl;
	fileM.createDataSet("Shape",H5::PredType::NATIVE_INT,dataspace_Shape).write(shape,H5::PredType::NATIVE_INT);	
	//fileM->createDataSet("Shape",H5::PredType::NATIVE_INT,dataspace_Shape).write(shape,H5::PredType::NATIVE_INT);	
	std::cout<<"\n\n\n"<<std::endl;
	
	//makign one group for ech aux detector
	std::cout<<"Create HDF5 Groups"<<std::endl;
	std::vector<H5::Group*> group_M;
	std::vector<H5::Group*> group_D;
	std::ostringstream oss;
	for (int i = 0; i < (int) auxchannels.size(); i++){
        oss<<"/AuxCH_"<<auxchannels[i];
		group_M.push_back(new H5::Group(fileM.createGroup(oss.str().c_str())));
		group_D.push_back(new H5::Group(fileD.createGroup(oss.str().c_str())));
		//group_M.push_back(new H5::Group(fileM->createGroup(oss.str().c_str())));
		//group_D.push_back(new H5::Group(fileD->createGroup(oss.str().c_str())));
        oss.str("");
        oss.clear();	
	}
	//Create the dataspace for the twi kind of data I will save
	std::cout<<"Create HDF5 DataSpaces"<<std::endl;
	hsize_t dimmat[1];
	dimmat[0] = nfreqbins*nfreqbins,//sizeof(double[nfreqbins]); //size_t(nfreqbins);
	std::cout<<"dim 1 = "<<1<<" -- dim 2 = "<<dimmat[0]<<std::endl;
	H5::DataSpace *dataspace_M = new H5::DataSpace(1, dimmat);	
	
	hsize_t dimmatD[1];
	dimmatD[0] = nfreqbins,//sizeof(double[nfreqbins]); //size_t(nfreqbins);
	std::cout<<"dim 1 = "<<1<<" -- dim 2 = "<<dimmatD[0]<<std::endl;
	H5::DataSpace *dataspace_D = new H5::DataSpace(1, dimmatD);	
	
	
	
	int chnoiseps = 0;
	std::vector<double> npsvec ;
	tnoise->Branch("ch",&chnoiseps);
	tnoise->Branch("nps",&npsvec);
	//cycle
	for (int i = 0; i < (int) auxchannels.size(); i++){
    int chaux = auxchannels[i];
    //declaring the matrix
		QMatrixC auxmatrix(nfreqbins,nfreqbins);
		//declaring the anps for the aux channel
    QVector auxANPSvec((uint)nfreqbins);    
    QVector boloANPS((uint)nfreqbins);    
		//initialize to zero	
		for(size_t jj = 0; jj<(uint)nfreqbins; jj++){
			auxANPSvec(jj)=0;
			boloANPS(jj)=0;
		}
	
		//initialize the counter
		count = 0;
		std::cout<<"CrossC with CH "<<chaux<<std::endl;
		//cycling the ok t0s

    QMatrixC temp(nfreqbins,nfreqbins);
		for(std::map<double,QVector>::iterator iter = VecTCH.begin(); iter != VecTCH.end(); ++iter)
		{
			double t0i =  iter->first;
			QVector ts = iter->second;
      //int n = ts.Size();
			std::cout<<"t0 = "<<t0i<<") "<<std::endl;

      //get the time stream
			GetTimestream(readersCH[chaux],chaux,t0i, t0i+twindow, tsaux,run,DownsamplingFactor_byCH[chaux]);
			//check if there is a pulse 	
      if(hasPulse(tsaux,threshmVCH[chaux],averageCH[chaux],debounceCH[chaux], adc2mV_byCH[chaux])){
        std::cout<<"Found pulse in the aux"<<std::endl;
        continue;
      }
      else{ //check that the eventual sides do not have a trigger!
        std::cout<<"CHannel "<<chaux<<" - Side search"<<SideMap[chaux].size()<<std::endl;
        if(SideMap[chaux].size()>0){
          bool DoesSideTrigger = false;
          for (std::vector<int>::const_iterator it = SideMap[chaux].begin(); it!=SideMap[chaux].end() && !DoesSideTrigger; ++it) {
            GetTimestream(readersCH[*it], *it,t0i,t0i+twindow, tside,run,DownsamplingFactor_byCH[*it]);
            DoesSideTrigger = hasPulse(tside,threshmVCH[*it],averageCH[*it],debounceCH[*it],adc2mV_byCH[*it]);
            
          }  
          if(DoesSideTrigger) continue;      
        }
      }
			//arriving here if no pulse is found --> double good
			//now doing the actual FFT-based analysis
			
			//Applying a rectangular windowing
      //TODO: add more windowing options
      ts[0] = 0;
      ts[ts.Size()-1] = 0;
      tsaux[0] = 0;
      tsaux[tsaux.Size()-1] = 0;
			
      //Then normalizing by the gain
      ts/=GainByChannel[ch];
      tsaux/=GainByChannel[chaux];

			//Gget the fft of the base channel
      std::cout<<"FFt main"<<std::endl;
      fTransformer.TransformToFreq(ts, fft1, true);//trunc = true -- returns a truncated fft
      fft1 *= adc2fft;
      QVector psd1 = fft1.GetMagnitudesSquare();
		
			//doing the same for the aux
      std::cout<<"FFt aux"<<std::endl;
      fTransformer.TransformToFreq(tsaux, auxfft, true);
			auxfft *= adc2fft;
      QVector psdaux = auxfft.GetMagnitudesSquare();

      
      


			//dimensional checks
			int nbins = auxfft.Size();
      std::cout << "nbins var is = " << nbins << std::endl;

      std::cout << "nfreqbins var is = " << nfreqbins << std::endl;

      std::cout << "fSampleFreq var is = " << fSampleFreq << std::endl;
      //return 1;

      if (nbins != (int) fft1.Size()){ 
          std::cout << "Error. Two channels have different PSD lengths. Cannot build a square matrix!" << std::endl;
          std::cout << Form("CH %d and Ch %d",ch,chaux) << std::endl;
          std::cout << "Exiting..." << std::endl;
          return 1;
      }
      if(nbins != nfreqbins){
          std::cout << "ERROR: nbins = " << nbins << " while nfreqbins = " << nfreqbins << std::endl;
          std::cout << "Exiting..." << std::endl;
          return 1;
      }   

			//in case a print is needed
			//plotTimestream(ts, run, ch, t0i, fSampleFreq, output, "MainEvent");
			//plotTimestream(tsaux, run, chaux, t0i, fSampleFreq, output, "AuxEvent");
			
			//saving the pulses to the ttree

      /* Commenting out the Root File creation for the pulses to see if this is the memory leak

      std::cout<<"Save the pulses to the ttree"<<std::endl;
      startt = t0i;
      chauxsave = chaux;
      pulseref.clear();   
      pulseaux.clear();        
      for(size_t iv = 0; iv<ts.Size(); iv++){
        pulseref.push_back(ts(iv));
        pulseaux.push_back(tsaux(iv));
      }
      tvec->Fill();	

      */

			/*
      * Where the magic happens
      */ 				
			
			//outer-product the two ffts row-by-row to fill a matrix
      std::cout<<"Outer prod"<<std::endl;
      std::cout << "Initializing temp" << nbins << std::endl;
      //QMatrixC temp(nbins,nbins);
      //QMatrixC temp = ::new QMatrixC(nbins,nbins);
      std::cout<<"Adding Setting rows of Temp"<<std::endl;             
      
      for(int j = 0; j < nbins; j++){
          QComplex rowval = fft1(j).Conj();
          temp.SetRow(j,rowval*auxfft);    
      }
      std::cout<<"Adding temp to auxmatrix"<<std::endl;             
      auxmatrix += temp;
			//build the anps  
      std::cout<<"Building ANPS aux"<<std::endl;             
      auxANPSvec += psdaux; 
      std::cout<<"Building ANPS target"<<std::endl;       
			boloANPS += psd1;
			//counting
      
      std::cout<<"Counting"<<std::endl;
			count++;
		
		}//cycle over times is over --> from now working on the avgs
		std::cout<<"-->Done"<<std::endl<<std::endl;

		

		std::cout << "Final count: " << count << " noise windows" << std::endl;
		double scale = (double) 1.0/count;
		std::cout << scale << std::endl;
		//normalize the ANPSs
		boloANPS *= scale;
		auxANPSvec *= scale;
    //save anps to ttree
		chnoiseps = chaux;
		npsvec.clear();           
		for(size_t iv = 0; iv<auxANPSvec.Size(); iv++){
			npsvec.push_back(auxANPSvec(iv));
		}
		tnoise->Fill();	

	
		//normalize the matrix
		std::cout << "Normalizing the matrix" << std::endl;
		auxmatrix *= scale;
		//Get the actual corrmatrix
		QMatrix CorrMat = auxmatrix.Magnitude();
		//Normalize the correlation matrix from 0 to 1

    //Because for this, I'm only interested in diagonals, replacing this with only one for loop
    //Commented out previous for loop for normalizing matrix
    /*
		for(int j = 0; j < (int) CorrMat.GetNRow(); j++){     
			for(int k = 0; k < (int) CorrMat.GetNCol(); k++){
				CorrMat(j,k) /= std::sqrt(boloANPS(j)*auxANPSvec(k));
				//if(j == k && j % 50 == 0) std::cout << "Rho[j,j] for j = " << j << ": " << CorrMat(j,k) << std::endl; 
			}
		}
    */

    for(int j = 0; j < (int) CorrMat.GetNRow(); j++){
      CorrMat(j,j) /= std::sqrt(boloANPS(j)*auxANPSvec(j));
    }

		std::cout << Form("Plotting aux %04d",chaux) << std::endl;

    /*
		WriteMats2H5(group_M[i],group_D[i],
						  dataspace_M, dataspace_D,
						  auxmatrix,boloANPS,auxANPSvec,fSampleFreq);
    */
    WriteDiagMats2H5(group_D[i],dataspace_D,auxmatrix,boloANPS,auxANPSvec,fSampleFreq);
    
		//berettam:
		//add the method to get the CH-CH correlation
		double tempchchcorr = 0;
		double normsum = 0;
		for(size_t kk = 0; kk<boloANPS.Size();kk++){
			normsum+=boloANPS(kk);
		}
		//just lookign at the diagonal and summing over that
		for(int j = 0; j < (int) CorrMat.GetNRow(); j++){     
			tempchchcorr += std::pow(CorrMat(j,j),2)*boloANPS(j)/normsum;
		}
		CHCHCorr.push_back(std::sqrt(tempchchcorr));
	}//-->cycle over aux devices is over	
		//closing the ttree and saing into the tfile
		fpulse->cd();
		tvec->Write();
		fpulse->Close();	
		fnoise->cd();
		tnoise->Write();
		fnoise->Close();	

	//cl;ose the Hdf5
	fileM.close();	
	fileD.close();
	delete dataspace_M;	
	delete dataspace_D;	
	
	//berettam: write the chch corr matrix
	std::ofstream outchch (Form("%s/chchcorr/CH%04d_toAuxCorr.txt",output.Data(),ch),std::ios::out);
    outchch<<"CH,aux,corr"<<std::endl;
	for (int i = 0; i < (int) auxchannels.size(); i++){
		outchch<<ch<<","<<auxchannels[i]<<","<<CHCHCorr[i]<<std::endl;		
	}

	outchch.close();
    std::cout << "Correlation matrix saved" << std::endl;
    return 0;
}
#endif
