

//
//  Following the example of tutorials/unfold/testUnfold1.c
//


#include "histio.c"

TH2F* get_hist( const char* hname ) {
   TH2F* hp = (TH2F*) gDirectory -> FindObject( hname ) ;
   if ( hp == 0x0 ) { printf("\n\n *** can't find %s\n\n", hname ) ; gSystem->Exit(-1) ; }
   return hp ;
}

//----------

   void unfold1( const char* hist_name = "h_log10_tau1b_gen_vs_obs_dnn", //"h_log10_x_gen_vs_obs_dnn",
                 int ngen = 1e5,
                 const char* input_file = "unfold-hists-input-nbins_gen020_obs035.root" ) {

      gDirectory -> Delete( "h*" ) ;

      loadHist( input_file ) ;

      printf("\n\n") ;
      gDirectory -> ls() ;
      printf("\n\n") ;

      TH2F* h_in_gen_vs_obs = get_hist( hist_name ) ;

      TH1* h_obs_source = h_in_gen_vs_obs -> ProjectionX( "h_obs_source" ) ;
      TH1* h_gen_source = h_in_gen_vs_obs -> ProjectionY( "h_gen_source" ) ;

      TH1* h_obs_random = (TH1*) h_obs_source -> Clone( "h_obs_random" ) ;
      h_obs_random->Reset() ;
      h_obs_random->FillRandom( h_obs_source, ngen ) ;


      TUnfoldDensity unfold( h_in_gen_vs_obs, TUnfold::kHistMapOutputVert ) ;

      int return_status = unfold.SetInput( h_obs_random ) ;
      printf("  Return status for SetInput: %d\n", return_status ) ;

      //========================================================================
      // the unfolding is done here
      //
      // scan L curve and find best point
      Int_t nScan=30;
      // use automatic L-curve scan: start with taumin=taumax=0.0
      Double_t tauMin=0.0;
      Double_t tauMax=0.0;
      Int_t iBest;
      TSpline *logTauX,*logTauY;
      TGraph *lCurve;

      // if required, report Info messages (for debugging the L-curve scan)
      Int_t oldinfo=gErrorIgnoreLevel;
      gErrorIgnoreLevel=kInfo;

      // this method scans the parameter tau and finds the kink in the L curve
      // finally, the unfolding is done for the best choice of tau
      iBest=unfold.ScanLcurve(nScan,tauMin,tauMax,&lCurve,&logTauX,&logTauY);

      // if required, switch to previous log-level
      gErrorIgnoreLevel=oldinfo;


      //==========================================================================
      // print some results
      //
      std::cout<<"tau="<<unfold.GetTau()<<"\n";
      std::cout<<"chi**2="<<unfold.GetChi2A()<<"+"<<unfold.GetChi2L()
               <<" / "<<unfold.GetNdf()<<"\n";

      //==========================================================================
      // create graphs with one point to visualize the best choice of tau
      //
      Double_t t[1],x[1],y[1];
      logTauX->GetKnot(iBest,t[0],x[0]);
      logTauY->GetKnot(iBest,t[0],y[0]);
      TGraph *bestLcurve=new TGraph(1,x,y);
      TGraph *bestLogTauLogChi2=new TGraph(1,t,x);

 //// TCanvas* can = new TCanvas( "can", "", 50, 50, 1500, 500 ) ;

 //// can -> Divide(3,1) ;


 //// can -> cd(1) ;
 //// lCurve->Draw("alp") ;
 //// bestLcurve->SetMarkerStyle(20) ;
 //// bestLcurve->Draw("p") ;

 //// can -> cd(2) ;
 //// logTauX->Draw() ;

 //// can -> cd(3) ;
 //// logTauY->Draw() ;


      //==========================================================================
      // retrieve results into histograms

      // get unfolded distribution
      TH1 *histMunfold=unfold.GetOutput("Unfolded");

      // get unfolding result, folded back
      TH1 *histMdetFold=unfold.GetFoldedOutput("FoldedBack");

      // get error matrix (input distribution [stat] errors only)
      // TH2D *histEmatData=unfold.GetEmatrix("EmatData");

      // get total error matrix:
      //   migration matrix uncorrelated and correlated systematic errors
      //   added in quadrature to the data statistical errors
      TH2 *histEmatTotal=unfold.GetEmatrixTotal("EmatTotal");

      // create data histogram with the total errors
      int nGen = histMunfold -> GetNbinsX() ;
      float xminGen = histMunfold -> GetXaxis() -> GetXmin() ;
      float xmaxGen = histMunfold -> GetXaxis() -> GetXmax() ;

      TH1D *histTotalError=
         new TH1D("TotalError","TotalError",nGen,xminGen,xmaxGen);
      for(Int_t bin=1;bin<=nGen;bin++) {
        histTotalError->SetBinContent(bin,histMunfold->GetBinContent(bin));
        histTotalError->SetBinError
           (bin,TMath::Sqrt(histEmatTotal->GetBinContent(bin,bin)));
      }


      // get global correlation coefficients
      // for this calculation one has to specify whether the
      // underflow/overflow bins are included or not
      // default: include all bins
      // here: exclude underflow and overflow bins
      TH2 *gHistInvEMatrix;
      TH1 *histRhoi=unfold.GetRhoItotal("rho_I",
                                        0, // use default title
                                        0, // all distributions
                                        "*[UO]", // discard underflow and overflow bins on all axes
                                        kTRUE, // use original binning
                                        &gHistInvEMatrix // store inverse of error matrix
                                        );



      gStyle -> SetOptStat(0) ;

      h_obs_source -> SetLineColor(4) ;
      h_gen_source -> SetLineColor(2) ;

      TH1* h_gen_compare = (TH1*) h_gen_source -> Clone( "h_gen_compare" ) ;
      h_gen_compare -> Scale( ( histMunfold -> Integral() )/( h_gen_source -> Integral() ) ) ;

      TCanvas* can1 = (TCanvas*) gDirectory -> FindObject( "can1" ) ;
      if ( can1 == 0x0 ) can1 = new TCanvas( "can1", "", 50, 50, 1000, 1000 ) ;
      can1 -> Clear() ;
      can1 -> cd() ;

      can1 -> Divide(2,2) ;


      can1 -> cd(1) ;
      h_in_gen_vs_obs -> Draw("colz") ;

      can1 -> cd(2) ;
      histMunfold -> Draw() ;
      h_gen_compare -> Draw("same hist") ;

      can1 -> cd(3) ;
      h_obs_source -> Draw("hist") ;

      can1 -> cd(4) ;
      h_gen_source -> Draw("hist") ;


      can1->SaveAs("test.png");
      can1->SaveAs("test.pdf");

   }













