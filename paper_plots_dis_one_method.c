

//
//  Following the example of tutorials/unfold/testUnfold1.c
//


//
//  This makes one plot per canvas and saves each plot to its own pdf file.
//
//

#include "histio.c"

#ifndef helpers
#define helpers

TH2F* get_hist( const char* hname ) {
   TH2F* hp = (TH2F*) gDirectory -> FindObject( hname ) ;
   if ( hp == 0x0 ) { printf("\n\n *** can't find %s\n\n", hname ) ; gSystem->Exit(-1) ; }
   return hp ;
}

void SetupCorrelationPalette() {

      static Bool_t initialized = kFALSE ;
      static Int_t colors[90] ;

      const Int_t Number = 3 ;
      Double_t Length[Number] = {0.,0.5, 1.} ;
      Double_t Red[Number] = {0.,1.,1.} ;
      Double_t Green[Number] = {0.,1.,0.} ;
      Double_t Blue[Number] = {1.,1.,0.} ;

      if ( !initialized ) {
         Int_t fi = TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,90);
         for ( int i=0; i<90; i++ ) colors[i] = fi+i ;
         initialized = kTRUE ;
         return ;
      }
      gStyle -> SetPalette( 90, colors ) ;

}

void Setup2DhistPalette() {
      gStyle -> SetPalette( kBird ) ;
}

#endif

//----------

   void paper_plots_dis_one_method(
                      const char* hist_name_a = "h_log10_x_gen_vs_obs_dnn",
                      const char* method_name = "DNN",
                      const char* var_name = "x",
                      int ngen = 1e5,
                      const char* input_file = "example-input-nbins_gen010_obs020.root"
                      ) {

      char htitle[1000] ;
      char fname[1000] ;

      gSystem -> Exec( "mkdir -p paper-plots" ) ;

      TRandom3* tran = new TRandom3(1249) ;

      gDirectory -> Delete( "h*" ) ;
      gStyle -> SetPalette( kBird ) ;

      loadHist( input_file ) ;

      printf("\n\n") ;
      gDirectory -> ls() ;
      printf("\n\n") ;

      TH2F* h_in_gen_vs_obs_a = get_hist( hist_name_a ) ;

      TH1* h_obs_source_a = h_in_gen_vs_obs_a -> ProjectionX( "h_obs_source_a" ) ;
      TH1* h_gen_source_a = h_in_gen_vs_obs_a -> ProjectionY( "h_gen_source_a" ) ;



      TH1* h_obs_random_a = (TH1*) h_obs_source_a -> Clone( "h_obs_random_a" ) ;
      h_obs_random_a->Reset() ;
      h_obs_random_a->FillRandom( h_obs_source_a, ngen, tran ) ;




     //-- Add underflows and overflows to fake data!
      float underflow_mean ;
      float overflow_mean ;
      int nbins_obs = h_obs_source_a -> GetNbinsX() ;
      float nobs_source ;
      float nobs_integral ;
      float nobs_underflow ;
      float nobs_overflow ;



      nobs_integral = h_obs_source_a -> Integral() ;
      nobs_underflow = h_obs_source_a -> GetBinContent(0) ;
      nobs_overflow = h_obs_source_a -> GetBinContent( nbins_obs+1 ) ;

      printf("  h_obs_source_a :  integral %9.5f, underflow %9.5f, overflow %9.5f\n", nobs_integral, nobs_underflow, nobs_overflow ) ;

      nobs_source = nobs_integral ;
      nobs_source += nobs_underflow ;
      nobs_source += nobs_overflow ;

      underflow_mean = ngen * ( nobs_underflow / nobs_integral ) ; // ngen only fills inside of histogram (integral)
      overflow_mean = ngen * ( nobs_overflow / nobs_integral ) ;

      h_obs_random_a -> SetBinContent( 0, tran->Poisson( underflow_mean ) ) ;
      h_obs_random_a -> SetBinContent( nbins_obs+1, tran->Poisson( overflow_mean ) ) ;






      TUnfoldDensity unfold_a( h_in_gen_vs_obs_a, TUnfold::kHistMapOutputVert ) ;

      int return_status_a = unfold_a.SetInput( h_obs_random_a ) ;
      printf("  Return status for SetInput A: %d\n", return_status_a ) ;




      //========================================================================
      // the unfolding is done here
      //
      // scan L curve and find best point
      Int_t nScan=30;
      // use automatic L-curve scan: start with taumin=taumax=0.0
      Double_t tauMin=0.0;
      Double_t tauMax=0.0;
      Int_t iBest;
      TSpline *logTauX_a,*logTauY_a;
      TSpline *logTauX_b,*logTauY_b;
      TGraph *lCurve_a;
      TGraph *lCurve_b;

      // if required, report Info messages (for debugging the L-curve scan)
      Int_t oldinfo=gErrorIgnoreLevel;
      gErrorIgnoreLevel=kInfo;

      // this method scans the parameter tau and finds the kink in the L curve
      // finally, the unfolding is done for the best choice of tau
      iBest=unfold_a.ScanLcurve( nScan, tauMin, tauMax, &lCurve_a, &logTauX_a, &logTauY_a );

      // if required, switch to previous log-level
      gErrorIgnoreLevel=oldinfo;


      //==========================================================================
      // print some results
      //
      std::cout<<"A tau="<<unfold_a.GetTau()<<"\n";
      std::cout<<"A chi**2="<<unfold_a.GetChi2A()<<"+"<<unfold_a.GetChi2L() <<" / "<<unfold_a.GetNdf()<<"\n";


      //==========================================================================
      // retrieve results into histograms

      // get unfolded distribution
      TH1 *histMunfold_a = unfold_a.GetOutput("Unfolded_a");
      sprintf( htitle, "Unfolded distribution, %s", method_name ) ;
      histMunfold_a -> SetTitle( htitle ) ;

      // get unfolding result, folded back
      TH1 *histMdetFold_a = unfold_a.GetFoldedOutput("FoldedBack_a");


      // get total error matrix:
      //   migration matrix uncorrelated and correlated systematic errors
      //   added in quadrature to the data statistical errors
      TH2 *histEmatTotal_a=unfold_a.GetEmatrixTotal("EmatTotal_a");
      histEmatTotal_a -> SetTitle( "Covariance matrix, a" ) ;


      //-- calculate the correlation coefficients matrix
      TH2* correlation_matrix_a = (TH2*) histEmatTotal_a -> Clone( "correlation_matrix_a" ) ;
      sprintf( htitle, "Correlation coefficients, %s", method_name ) ;
      correlation_matrix_a -> SetTitle( htitle ) ;

      for ( int xbi=1; xbi<=histEmatTotal_a->GetNbinsX(); xbi++ ) {
         for ( int ybi=1; ybi<=histEmatTotal_a->GetNbinsY(); ybi++ ) {

            float rho = 1. ;
            float sigma_x2, sigma_y2 ;

            sigma_x2 = histEmatTotal_a->GetBinContent( xbi, xbi ) ;
            sigma_y2 = histEmatTotal_a->GetBinContent( ybi, ybi ) ;
            if ( sigma_x2 > 0 && sigma_y2 > 0 ) {
               rho = histEmatTotal_a -> GetBinContent( xbi, ybi )  / sqrt( sigma_x2 * sigma_y2 ) ;
            }
            correlation_matrix_a -> SetBinContent( xbi, ybi, rho ) ;


         } // ybi
      } // xbi

      correlation_matrix_a -> SetMinimum( -1. ) ;

      correlation_matrix_a -> SetMaximum(  1. ) ;



      // create data histogram with the total errors
      int nGen = histMunfold_a -> GetNbinsX() ;
      float xminGen = histMunfold_a -> GetXaxis() -> GetXmin() ;
      float xmaxGen = histMunfold_a -> GetXaxis() -> GetXmax() ;

      TH1D *histTotalError_a = new TH1D("TotalError_a","TotalError_a",nGen,xminGen,xmaxGen);
      for(Int_t bin=1;bin<=nGen;bin++) {
        histTotalError_a->SetBinContent(bin,histMunfold_a->GetBinContent(bin));
        histTotalError_a->SetBinError(bin,TMath::Sqrt(histEmatTotal_a->GetBinContent(bin,bin)));
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

      histRhoi_a -> SetTitle( "Global correlation, a" ) ;


      gStyle -> SetOptStat(0) ;

      h_obs_source_a -> SetLineColor(4) ;
      h_gen_source_a -> SetLineColor(2) ;

      TH1* h_gen_compare_a = (TH1*) h_gen_source_a -> Clone( "h_gen_compare_a" ) ;
      h_gen_compare_a -> Scale( ( histMunfold_a -> Integral() )/( h_gen_source_a -> Integral() ) ) ;




      TH1* h_unfold_err_a = (TH1*) histMunfold_a -> Clone( "h_unfold_err_a" ) ;


      h_unfold_err_a -> SetFillColor( kBlue-9 ) ;
      h_unfold_err_a -> SetFillStyle( 1001 ) ;


      for ( int bi=1; bi<=histMunfold_a->GetNbinsX(); bi++ ) {
         h_unfold_err_a -> SetBinContent( bi, 0. ) ;
      }








      histRhoi_a->SetMaximum(1.1) ;


      gStyle -> SetPadRightMargin(0.18) ;
      gStyle -> SetPadLeftMargin(0.15) ;
      gStyle -> SetPadBottomMargin(0.15) ;
      gStyle -> SetTitleBorderSize(0) ;
      gStyle -> SetTitleY(0.975) ;



      int can_width = 600 ;
      int can_height = 600 ;

      int can_spacing = 20 ;


     //-----

      TCanvas* can1 = (TCanvas*) gDirectory -> FindObject( "can1" ) ;
      if ( can1 == 0x0 ) { printf( "Making can1\n") ; can1 = new TCanvas( "can1", "", can_spacing, can_spacing, can_width, can_height ) ; }
      can1 -> Clear() ;
      can1 -> cd() ;

      h_in_gen_vs_obs_a -> SetTitleOffset( 1.2, "x" ) ;
      h_in_gen_vs_obs_a -> SetTitleOffset( 1.6, "y" ) ;
      sprintf( htitle, "Response matrix, %s", method_name ) ;
      h_in_gen_vs_obs_a -> SetTitle( htitle ) ;

      h_in_gen_vs_obs_a -> Draw("colz") ;
      TExec* change_hist_palette = new TExec( "change_hist_palette", "Setup2DhistPalette();" );
      change_hist_palette->Draw() ;
      h_in_gen_vs_obs_a -> Draw("colz same") ;
      h_in_gen_vs_obs_a -> Draw("axis same") ;

      can1 -> Update() ; can1 -> Draw() ; gSystem -> ProcessEvents() ;

      sprintf( fname, "paper-plots/dis-response-%s-%s.pdf", var_name, method_name ) ;
      can1 -> SaveAs( fname ) ;
      sprintf( fname, "paper-plots/dis-response-%s-%s.png", var_name, method_name ) ;
      can1 -> SaveAs( fname ) ;


     //-----

      gStyle -> SetPadRightMargin(0.05) ;

      TCanvas* can2 = (TCanvas*) gDirectory -> FindObject( "can2" ) ;
      if ( can2 == 0x0 ) { printf( "Making can2\n") ; can2 = new TCanvas( "can2", "", 2*can_spacing + can_width , can_spacing, can_width, can_height ) ; }
      can2 -> Clear() ;
      can2 -> cd() ;

      histMunfold_a -> SetTitleOffset( 1.2, "x" ) ;
      histMunfold_a -> SetTitleOffset( 1.8, "y" ) ;
      histMunfold_a -> SetYTitle( "Events" ) ;

      if ( strcmp( var_name, "x" ) == 0 ) { histMunfold_a -> SetXTitle( "log10(x)" ) ; }
      if ( strcmp( var_name, "y" ) == 0 ) { histMunfold_a -> SetXTitle( "log10(y)" ) ; }

      histMunfold_a -> SetLineWidth(3) ;
      h_gen_compare_a -> SetLineWidth(3) ;

      histMunfold_a -> Draw() ;
      h_gen_compare_a -> Draw("same hist") ;
      histMunfold_a -> Draw( "same" ) ;

      TLegend* legend = new TLegend( 0.40, 0.20, 0.70, 0.35 ) ;
      legend -> AddEntry( histMunfold_a, "Unfolded" ) ;
      legend -> AddEntry( h_gen_compare_a, "Gen" ) ;
      legend -> Draw() ;

      can2 -> Update() ; can2 -> Draw() ; gSystem -> ProcessEvents() ;

      sprintf( fname, "paper-plots/dis-unfolded-%s-%s.pdf", var_name, method_name ) ;
      can2 -> SaveAs( fname ) ;
      sprintf( fname, "paper-plots/dis-unfolded-%s-%s.png", var_name, method_name ) ;
      can2 -> SaveAs( fname ) ;

     //-----

      gStyle -> SetPadRightMargin(0.18) ;

      TCanvas* can3 = (TCanvas*) gDirectory -> FindObject( "can3" ) ;
      if ( can3 == 0x0 ) { printf( "Making can3\n") ; can3 = new TCanvas( "can3", "", 3*can_spacing + 2*can_width , can_spacing, can_width, can_height ) ; }
      can3 -> Clear() ;
      can3 -> cd() ;

      correlation_matrix_a -> SetTitleOffset( 1.2, "x" ) ;
      correlation_matrix_a -> SetTitleOffset( 1.6, "y" ) ;

      if ( strcmp( var_name, "x" ) == 0 ) { correlation_matrix_a -> SetXTitle( "log10(x)" ) ; correlation_matrix_a -> SetYTitle( "log10(x)" ) ; }
      if ( strcmp( var_name, "y" ) == 0 ) { correlation_matrix_a -> SetXTitle( "log10(y)" ) ; correlation_matrix_a -> SetYTitle( "log10(y)" ) ; }

      correlation_matrix_a -> Draw( "colz" ) ;
      TExec* change_cor_palette = new TExec( "change_cor_palette", "SetupCorrelationPalette();" );
      change_cor_palette->Draw() ;
      correlation_matrix_a -> Draw( "colz same" ) ;
      correlation_matrix_a -> Draw( "axis same" ) ;

      can3 -> Update() ; can3 -> Draw() ; gSystem -> ProcessEvents() ;

      sprintf( fname, "paper-plots/dis-cormat-%s-%s.pdf", var_name, method_name ) ;
      can3 -> SaveAs( fname ) ;
      sprintf( fname, "paper-plots/dis-cormat-%s-%s.png", var_name, method_name ) ;
      can3 -> SaveAs( fname ) ;

     //-----


      printf("\n\n\n") ;

      printf(" cut and paste for this:\n\n") ;

      printf("     paper_plots_dis_one_method(\"%s\",\"%s\",\"%s\",%d,\"%s\")\n", hist_name_a, method_name, var_name, ngen, input_file ) ;


      printf("\n\n\n") ;

   }













