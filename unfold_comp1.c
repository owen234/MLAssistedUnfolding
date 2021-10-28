

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

   void unfold_comp1( const char* hist_name_a = "h_log10_x_gen_vs_obs_dnn",
                      const char* hist_name_b = "h_log10_x_gen_vs_obs_e",
                      int ngen = 1e5,
                      const char* meth_a_label = "DNN",
                      const char* meth_b_label = "e",
                      const char* var_label = "log10(x)",
                      const char* input_file = "unfold-hists-input-nbins_gen020_obs050.root" ) {

      gDirectory -> Delete( "h*" ) ;

      loadHist( input_file ) ;

      printf("\n\n") ;
      gDirectory -> ls() ;
      printf("\n\n") ;

      TH2F* h_in_gen_vs_obs_a = get_hist( hist_name_a ) ;
      TH2F* h_in_gen_vs_obs_b = get_hist( hist_name_b ) ;

      TH1* h_obs_source_a = h_in_gen_vs_obs_a -> ProjectionX( "h_obs_source_a" ) ;
      TH1* h_gen_source_a = h_in_gen_vs_obs_a -> ProjectionY( "h_gen_source_a" ) ;

      TH1* h_obs_source_b = h_in_gen_vs_obs_b -> ProjectionX( "h_obs_source_b" ) ;
      TH1* h_gen_source_b = h_in_gen_vs_obs_b -> ProjectionY( "h_gen_source_b" ) ;


      TH1* h_obs_random_a = (TH1*) h_obs_source_a -> Clone( "h_obs_random_a" ) ;
      h_obs_random_a->Reset() ;
      h_obs_random_a->FillRandom( h_obs_source_a, ngen ) ;

      TH1* h_obs_random_b = (TH1*) h_obs_source_b -> Clone( "h_obs_random_b" ) ;
      h_obs_random_b->Reset() ;
      h_obs_random_b->FillRandom( h_obs_source_b, ngen ) ;


      TUnfoldDensity unfold_a( h_in_gen_vs_obs_a, TUnfold::kHistMapOutputVert ) ;
      TUnfoldDensity unfold_b( h_in_gen_vs_obs_b, TUnfold::kHistMapOutputVert ) ;

      int return_status_a = unfold_a.SetInput( h_obs_random_a ) ;
      printf("  Return status for SetInput A: %d\n", return_status_a ) ;

      int return_status_b = unfold_b.SetInput( h_obs_random_b ) ;
      printf("  Return status for SetInput B: %d\n", return_status_b ) ;



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
      iBest=unfold_a.ScanLcurve(nScan,tauMin,tauMax,&lCurve,&logTauX,&logTauY);
      iBest=unfold_b.ScanLcurve(nScan,tauMin,tauMax,&lCurve,&logTauX,&logTauY);

      // if required, switch to previous log-level
      gErrorIgnoreLevel=oldinfo;


      //==========================================================================
      // print some results
      //
      std::cout<<"A tau="<<unfold_a.GetTau()<<"\n";
      std::cout<<"A chi**2="<<unfold_a.GetChi2A()<<"+"<<unfold_a.GetChi2L() <<" / "<<unfold_a.GetNdf()<<"\n";

      std::cout<<"B tau="<<unfold_b.GetTau()<<"\n";
      std::cout<<"B chi**2="<<unfold_b.GetChi2A()<<"+"<<unfold_b.GetChi2L() <<" / "<<unfold_b.GetNdf()<<"\n";

      //==========================================================================
      // retrieve results into histograms

      // get unfolded distribution
      TH1 *histMunfold_a = unfold_a.GetOutput("Unfolded_a");
      TH1 *histMunfold_b = unfold_b.GetOutput("Unfolded_b");

      // get unfolding result, folded back
      TH1 *histMdetFold_a = unfold_a.GetFoldedOutput("FoldedBack_a");
      TH1 *histMdetFold_b = unfold_b.GetFoldedOutput("FoldedBack_b");

      // get error matrix (input distribution [stat] errors only)
      // TH2D *histEmatData=unfold.GetEmatrix("EmatData");

      // get total error matrix:
      //   migration matrix uncorrelated and correlated systematic errors
      //   added in quadrature to the data statistical errors
      TH2 *histEmatTotal_a=unfold_a.GetEmatrixTotal("EmatTotal_a");
      TH2 *histEmatTotal_b=unfold_b.GetEmatrixTotal("EmatTotal_b");

      // create data histogram with the total errors
      int nGen = histMunfold_a -> GetNbinsX() ;
      float xminGen = histMunfold_a -> GetXaxis() -> GetXmin() ;
      float xmaxGen = histMunfold_a -> GetXaxis() -> GetXmax() ;

      TH1D *histTotalError_a = new TH1D("TotalError_a","TotalError_a",nGen,xminGen,xmaxGen);
      TH1D *histTotalError_b = new TH1D("TotalError_b","TotalError_b",nGen,xminGen,xmaxGen);
      for(Int_t bin=1;bin<=nGen;bin++) {
        histTotalError_a->SetBinContent(bin,histMunfold_a->GetBinContent(bin));
        histTotalError_a->SetBinError(bin,TMath::Sqrt(histEmatTotal_a->GetBinContent(bin,bin)));
        histTotalError_b->SetBinContent(bin,histMunfold_b->GetBinContent(bin));
        histTotalError_b->SetBinError(bin,TMath::Sqrt(histEmatTotal_b->GetBinContent(bin,bin)));
      }


      // get global correlation coefficients
      // for this calculation one has to specify whether the
      // underflow/overflow bins are included or not
      // default: include all bins
      // here: exclude underflow and overflow bins
      TH2 *gHistInvEMatrix;
      TH1 *histRhoi_a = unfold_a.GetRhoItotal("rho_I_a",
                                        0, // use default title
                                        0, // all distributions
                                        "*[UO]", // discard underflow and overflow bins on all axes
                                        kTRUE, // use original binning
                                        &gHistInvEMatrix // store inverse of error matrix
                                        );
      TH1 *histRhoi_b = unfold_b.GetRhoItotal("rho_I_b",
                                        0, // use default title
                                        0, // all distributions
                                        "*[UO]", // discard underflow and overflow bins on all axes
                                        kTRUE, // use original binning
                                        &gHistInvEMatrix // store inverse of error matrix
                                        );



      gStyle -> SetOptStat(0) ;

      h_obs_source_a -> SetLineColor(4) ;
      h_gen_source_a -> SetLineColor(2) ;

      TH1* h_gen_compare_a = (TH1*) h_gen_source_a -> Clone( "h_gen_compare_a" ) ;
      h_gen_compare_a -> Scale( ( histMunfold_a -> Integral() )/( h_gen_source_a -> Integral() ) ) ;


      h_obs_source_b -> SetLineColor(4) ;
      h_gen_source_b -> SetLineColor(2) ;

      TH1* h_gen_compare_b = (TH1*) h_gen_source_b -> Clone( "h_gen_compare_b" ) ;
      h_gen_compare_b -> Scale( ( histMunfold_b -> Integral() )/( h_gen_source_b -> Integral() ) ) ;


      TH1* h_unfold_err_a = (TH1*) histMunfold_a -> Clone( "h_unfold_err_a" ) ;
      TH1* h_unfold_err_b = (TH1*) histMunfold_b -> Clone( "h_unfold_err_b" ) ;

      TH1* h_unfold_err_ratio = (TH1*) histMunfold_b -> Clone( "h_unfold_err_ratio" ) ;

      h_unfold_err_ratio -> SetYTitle( "Unfolded error ratio" ) ;
      h_unfold_err_ratio -> SetTitle( "" ) ;

      h_unfold_err_a -> SetFillColor( kBlue-9 ) ;
      h_unfold_err_b -> SetFillColor( kRed-9 ) ;
      h_unfold_err_a -> SetFillStyle( 1001 ) ;
      h_unfold_err_b -> SetFillStyle( 1001 ) ;

      h_unfold_err_b -> SetTitle( "Unfolded error" ) ;

      for ( int bi=1; bi<=histMunfold_a->GetNbinsX(); bi++ ) {
         h_unfold_err_a -> SetBinContent( bi, 0. ) ;
         h_unfold_err_b -> SetBinContent( bi, 0. ) ;
         float err_a = h_unfold_err_a -> GetBinError( bi ) ;
         float err_b = h_unfold_err_b -> GetBinError( bi ) ;
         float ratio = 0  ;
         if ( err_b > 0 ) ratio = err_a / err_b ;
         h_unfold_err_ratio -> SetBinContent( bi, ratio ) ;
         h_unfold_err_ratio -> SetBinError( bi, 0. ) ;
      }

      h_unfold_err_ratio -> SetMaximum(1.1) ;


      TCanvas* can1 = (TCanvas*) gDirectory -> FindObject( "can1" ) ;
      if ( can1 == 0x0 ) can1 = new TCanvas( "can1", "", 50, 50, 800, 1300 ) ;
      can1 -> Clear() ;
      can1 -> cd() ;

      can1 -> Divide(2,3) ;


      can1 -> cd(1) ;
      h_in_gen_vs_obs_a -> Draw("colz") ;

      can1 -> cd(2) ;
      histMunfold_a -> Draw() ;
      h_gen_compare_a -> Draw("same hist") ;


      can1 -> cd(3) ;
      h_in_gen_vs_obs_b -> Draw("colz") ;

      can1 -> cd(4) ;
      histMunfold_b -> Draw() ;
      h_gen_compare_b -> Draw("same hist") ;

      can1 -> cd(5) ;
      h_unfold_err_ratio -> Draw( "hist" ) ;
      gPad -> SetGridy(1) ;

      can1 -> cd(6) ;
      h_unfold_err_b -> Draw( "E2" ) ;
      h_unfold_err_a -> Draw("E2 same") ;
      gPad -> SetGridy(1) ;
      h_unfold_err_a -> Draw("axig same") ;


   }













