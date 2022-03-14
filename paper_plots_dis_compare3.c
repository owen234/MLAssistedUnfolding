

//
//  Following the example of tutorials/unfold/testUnfold1.c
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

   void paper_plots_dis_compare3(
                      const char* hist_name_a = "h_log10_x_gen_vs_obs_dnn",
                      const char* hist_name_b = "h_log10_x_gen_vs_obs_esigma",
                      const char* hist_name_c = "h_log10_x_gen_vs_obs_e",
                      const char* var_name = "x",
                      int ngen = 1e5,
                      const char* input_file = "paper-plots-input-1D-nbins_gen010_obs020.root",
                      float max_error = -1.
                      ) {

      char fname[1000] ;

      TRandom3* tran_a = new TRandom3(1249) ;
      TRandom3* tran_b = new TRandom3(1249) ;
      TRandom3* tran_c = new TRandom3(1249) ;


      gDirectory -> Delete( "h*" ) ;
      gStyle -> SetPalette( kBird ) ;

      loadHist( input_file ) ;

      printf("\n\n") ;
      gDirectory -> ls() ;
      printf("\n\n") ;

      TH2F* h_in_gen_vs_obs_a = get_hist( hist_name_a ) ;
      TH2F* h_in_gen_vs_obs_b = get_hist( hist_name_b ) ;
      TH2F* h_in_gen_vs_obs_c = get_hist( hist_name_c ) ;

      TH1* h_obs_source_a = h_in_gen_vs_obs_a -> ProjectionX( "h_obs_source_a" ) ;
      TH1* h_gen_source_a = h_in_gen_vs_obs_a -> ProjectionY( "h_gen_source_a" ) ;

      TH1* h_obs_source_b = h_in_gen_vs_obs_b -> ProjectionX( "h_obs_source_b" ) ;
      TH1* h_gen_source_b = h_in_gen_vs_obs_b -> ProjectionY( "h_gen_source_b" ) ;

      TH1* h_obs_source_c = h_in_gen_vs_obs_c -> ProjectionX( "h_obs_source_c" ) ;
      TH1* h_gen_source_c = h_in_gen_vs_obs_c -> ProjectionY( "h_gen_source_c" ) ;


      TH1* h_obs_random_a = (TH1*) h_obs_source_a -> Clone( "h_obs_random_a" ) ;
      h_obs_random_a->Reset() ;
      h_obs_random_a->FillRandom( h_obs_source_a, ngen, tran_a ) ;

      TH1* h_obs_random_b = (TH1*) h_obs_source_b -> Clone( "h_obs_random_b" ) ;
      h_obs_random_b->Reset() ;
      h_obs_random_b->FillRandom( h_obs_source_b, ngen, tran_b ) ;

      TH1* h_obs_random_c = (TH1*) h_obs_source_c -> Clone( "h_obs_random_c" ) ;
      h_obs_random_c->Reset() ;
      h_obs_random_c->FillRandom( h_obs_source_c, ngen, tran_c ) ;




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

      h_obs_random_a -> SetBinContent( 0, tran_a->Poisson( underflow_mean ) ) ;
      h_obs_random_a -> SetBinContent( nbins_obs+1, tran_a->Poisson( overflow_mean ) ) ;






      nobs_integral = h_obs_source_b -> Integral() ;
      nobs_underflow = h_obs_source_b -> GetBinContent(0) ;
      nobs_overflow = h_obs_source_b -> GetBinContent( nbins_obs+1 ) ;

      printf("  h_obs_source_b :  integral %9.5f, underflow %9.5f, overflow %9.5f\n", nobs_integral, nobs_underflow, nobs_overflow ) ;

      nobs_source = nobs_integral ;
      nobs_source += nobs_underflow ;
      nobs_source += nobs_overflow ;


      underflow_mean = ngen * ( nobs_underflow / nobs_integral ) ; // ngen only fills inside of histogram (integral)
      overflow_mean = ngen * ( nobs_overflow / nobs_integral ) ;


      h_obs_random_b -> SetBinContent( 0, tran_b->Poisson( underflow_mean ) ) ;
      h_obs_random_b -> SetBinContent( nbins_obs+1, tran_b->Poisson( overflow_mean ) ) ;






      nobs_integral = h_obs_source_c -> Integral() ;
      nobs_underflow = h_obs_source_c -> GetBinContent(0) ;
      nobs_overflow = h_obs_source_c -> GetBinContent( nbins_obs+1 ) ;

      printf("  h_obs_source_c :  integral %9.5f, underflow %9.5f, overflow %9.5f\n", nobs_integral, nobs_underflow, nobs_overflow ) ;

      nobs_source = nobs_integral ;
      nobs_source += nobs_underflow ;
      nobs_source += nobs_overflow ;


      underflow_mean = ngen * ( nobs_underflow / nobs_integral ) ; // ngen only fills inside of histogram (integral)
      overflow_mean = ngen * ( nobs_overflow / nobs_integral ) ;


      h_obs_random_c -> SetBinContent( 0, tran_c->Poisson( underflow_mean ) ) ;
      h_obs_random_c -> SetBinContent( nbins_obs+1, tran_c->Poisson( overflow_mean ) ) ;








      TUnfoldDensity unfold_a( h_in_gen_vs_obs_a, TUnfold::kHistMapOutputVert ) ;
      TUnfoldDensity unfold_b( h_in_gen_vs_obs_b, TUnfold::kHistMapOutputVert ) ;
      TUnfoldDensity unfold_c( h_in_gen_vs_obs_c, TUnfold::kHistMapOutputVert ) ;

      int return_status_a = unfold_a.SetInput( h_obs_random_a ) ;
      printf("  Return status for SetInput A: %d\n", return_status_a ) ;

      int return_status_b = unfold_b.SetInput( h_obs_random_b ) ;
      printf("  Return status for SetInput B: %d\n", return_status_b ) ;

      int return_status_c = unfold_c.SetInput( h_obs_random_c ) ;
      printf("  Return status for SetInput B: %d\n", return_status_c ) ;



      //========================================================================
      // the unfolding is done here
      //
      // scan L curve and find best point
      Int_t nScan=30;
      // use automatic L-curve scan: start with taumin=taumax=0.0
      Double_t tauMin=0.0;
      Double_t tauMax=0.0;
      //Double_t tauMin=0.00001;
      //Double_t tauMax=0.1;
      Int_t iBest;
      TSpline *logTauX_a,*logTauY_a;
      TSpline *logTauX_b,*logTauY_b;
      TSpline *logTauX_c,*logTauY_c;
      TGraph *lCurve_a;
      TGraph *lCurve_b;
      TGraph *lCurve_c;

      // if required, report Info messages (for debugging the L-curve scan)
      Int_t oldinfo=gErrorIgnoreLevel;
      gErrorIgnoreLevel=kInfo;

      // this method scans the parameter tau and finds the kink in the L curve
      // finally, the unfolding is done for the best choice of tau
      iBest=unfold_a.ScanLcurve( nScan, tauMin, tauMax, &lCurve_a, &logTauX_a, &logTauY_a );
      iBest=unfold_b.ScanLcurve( nScan, tauMin, tauMax, &lCurve_b, &logTauX_b, &logTauY_b );
      iBest=unfold_c.ScanLcurve( nScan, tauMin, tauMax, &lCurve_c, &logTauX_c, &logTauY_c );

      // if required, switch to previous log-level
      gErrorIgnoreLevel=oldinfo;


      //==========================================================================
      // print some results
      //
      std::cout<<"A tau="<<unfold_a.GetTau()<<"\n";
      std::cout<<"A chi**2="<<unfold_a.GetChi2A()<<"+"<<unfold_a.GetChi2L() <<" / "<<unfold_a.GetNdf()<<"\n";

      std::cout<<"B tau="<<unfold_b.GetTau()<<"\n";
      std::cout<<"B chi**2="<<unfold_b.GetChi2A()<<"+"<<unfold_b.GetChi2L() <<" / "<<unfold_b.GetNdf()<<"\n";

      std::cout<<"C tau="<<unfold_c.GetTau()<<"\n";
      std::cout<<"C chi**2="<<unfold_c.GetChi2A()<<"+"<<unfold_c.GetChi2L() <<" / "<<unfold_c.GetNdf()<<"\n";

      //==========================================================================
      // retrieve results into histograms

      // get unfolded distribution
      TH1 *histMunfold_a = unfold_a.GetOutput("Unfolded_a");
      TH1 *histMunfold_b = unfold_b.GetOutput("Unfolded_b");
      TH1 *histMunfold_c = unfold_c.GetOutput("Unfolded_c");
      histMunfold_a -> SetTitle( "Unfolded, a" ) ;
      histMunfold_b -> SetTitle( "Unfolded, b" ) ;
      histMunfold_c -> SetTitle( "Unfolded, c" ) ;

      // get unfolding result, folded back
      TH1 *histMdetFold_a = unfold_a.GetFoldedOutput("FoldedBack_a");
      TH1 *histMdetFold_b = unfold_b.GetFoldedOutput("FoldedBack_b");
      TH1 *histMdetFold_c = unfold_c.GetFoldedOutput("FoldedBack_c");

      // get error matrix (input distribution [stat] errors only)
      // TH2D *histEmatData=unfold.GetEmatrix("EmatData");

      // get total error matrix:
      //   migration matrix uncorrelated and correlated systematic errors
      //   added in quadrature to the data statistical errors
      TH2 *histEmatTotal_a=unfold_a.GetEmatrixTotal("EmatTotal_a");
      TH2 *histEmatTotal_b=unfold_b.GetEmatrixTotal("EmatTotal_b");
      TH2 *histEmatTotal_c=unfold_c.GetEmatrixTotal("EmatTotal_c");
      histEmatTotal_a -> SetTitle( "Covariance matrix, a" ) ;
      histEmatTotal_b -> SetTitle( "Covariance matrix, b" ) ;
      histEmatTotal_c -> SetTitle( "Covariance matrix, c" ) ;


      //-- calculate the correlation coefficients matrix
      TH2* correlation_matrix_a = (TH2*) histEmatTotal_a -> Clone( "correlation_matrix_a" ) ;
      TH2* correlation_matrix_b = (TH2*) histEmatTotal_b -> Clone( "correlation_matrix_b" ) ;
      TH2* correlation_matrix_c = (TH2*) histEmatTotal_c -> Clone( "correlation_matrix_c" ) ;
      correlation_matrix_a -> SetTitle( "Correlation coefficients, a" ) ;
      correlation_matrix_b -> SetTitle( "Correlation coefficients, b" ) ;
      correlation_matrix_c -> SetTitle( "Correlation coefficients, c" ) ;

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

            sigma_x2 = histEmatTotal_b->GetBinContent( xbi, xbi ) ;
            sigma_y2 = histEmatTotal_b->GetBinContent( ybi, ybi ) ;
            if ( sigma_x2 > 0 && sigma_y2 > 0 ) {
               rho = histEmatTotal_b -> GetBinContent( xbi, ybi )  / sqrt( sigma_x2 * sigma_y2 ) ;
            }
            correlation_matrix_b -> SetBinContent( xbi, ybi, rho ) ;

            sigma_x2 = histEmatTotal_c->GetBinContent( xbi, xbi ) ;
            sigma_y2 = histEmatTotal_c->GetBinContent( ybi, ybi ) ;
            if ( sigma_x2 > 0 && sigma_y2 > 0 ) {
               rho = histEmatTotal_c -> GetBinContent( xbi, ybi )  / sqrt( sigma_x2 * sigma_y2 ) ;
            }
            correlation_matrix_c -> SetBinContent( xbi, ybi, rho ) ;

         } // ybi
      } // xbi

      correlation_matrix_a -> SetMinimum( -1. ) ;
      correlation_matrix_b -> SetMinimum( -1. ) ;
      correlation_matrix_c -> SetMinimum( -1. ) ;

      correlation_matrix_a -> SetMaximum(  1. ) ;
      correlation_matrix_b -> SetMaximum(  1. ) ;
      correlation_matrix_c -> SetMaximum(  1. ) ;



      // create data histogram with the total errors
      int nGen = histMunfold_a -> GetNbinsX() ;
      float xminGen = histMunfold_a -> GetXaxis() -> GetXmin() ;
      float xmaxGen = histMunfold_a -> GetXaxis() -> GetXmax() ;

      TH1D *histTotalError_a = new TH1D("TotalError_a","TotalError_a",nGen,xminGen,xmaxGen);
      TH1D *histTotalError_b = new TH1D("TotalError_b","TotalError_b",nGen,xminGen,xmaxGen);
      TH1D *histTotalError_c = new TH1D("TotalError_c","TotalError_c",nGen,xminGen,xmaxGen);
      for(Int_t bin=1;bin<=nGen;bin++) {
        histTotalError_a->SetBinContent(bin,histMunfold_a->GetBinContent(bin));
        histTotalError_a->SetBinError(bin,TMath::Sqrt(histEmatTotal_a->GetBinContent(bin,bin)));
        histTotalError_b->SetBinContent(bin,histMunfold_b->GetBinContent(bin));
        histTotalError_b->SetBinError(bin,TMath::Sqrt(histEmatTotal_b->GetBinContent(bin,bin)));
        histTotalError_c->SetBinContent(bin,histMunfold_c->GetBinContent(bin));
        histTotalError_c->SetBinError(bin,TMath::Sqrt(histEmatTotal_c->GetBinContent(bin,bin)));
      }


      // get global correlation coefficients
      // for this calculation one has to specify whether the
      // underflow/overflow bins are included or not
      // default: include all bins
      // here: exclude underflow and overflow bins
      //
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

      TH1 *histRhoi_c = unfold_c.GetRhoItotal("rho_I_c",
                                        0, // use default title
                                        0, // all distributions
                                        "*[UO]", // discard underflow and overflow bins on all axes
                                        kTRUE, // use original binning
                                        &gHistInvEMatrix // store inverse of error matrix
                                        );

      histRhoi_a -> SetTitle( "Global correlation" ) ;
      histRhoi_b -> SetTitle( "Global correlation" ) ;
      histRhoi_c -> SetTitle( "Global correlation" ) ;


      gStyle -> SetOptStat(0) ;



      h_obs_source_a -> SetLineColor(4) ;
      h_gen_source_a -> SetLineColor(2) ;

      TH1* h_gen_compare_a = (TH1*) h_gen_source_a -> Clone( "h_gen_compare_a" ) ;
      h_gen_compare_a -> Scale( ( histMunfold_a -> Integral() )/( h_gen_source_a -> Integral() ) ) ;



      h_obs_source_b -> SetLineColor(4) ;
      h_gen_source_b -> SetLineColor(2) ;

      TH1* h_gen_compare_b = (TH1*) h_gen_source_b -> Clone( "h_gen_compare_b" ) ;
      h_gen_compare_b -> Scale( ( histMunfold_b -> Integral() )/( h_gen_source_b -> Integral() ) ) ;



      h_obs_source_c -> SetLineColor(4) ;
      h_gen_source_c -> SetLineColor(2) ;

      TH1* h_gen_compare_c = (TH1*) h_gen_source_c -> Clone( "h_gen_compare_c" ) ;
      h_gen_compare_c -> Scale( ( histMunfold_c -> Integral() )/( h_gen_source_c -> Integral() ) ) ;







      TH1* h_unfold_err_a = (TH1*) histMunfold_a -> Clone( "h_unfold_err_a" ) ;
      TH1* h_unfold_err_b = (TH1*) histMunfold_b -> Clone( "h_unfold_err_b" ) ;
      TH1* h_unfold_err_c = (TH1*) histMunfold_c -> Clone( "h_unfold_err_c" ) ;

      h_unfold_err_a -> SetTitle( "Unfolded error" ) ;
      h_unfold_err_b -> SetTitle( "Unfolded error" ) ;
      h_unfold_err_c -> SetTitle( "Unfolded error" ) ;


      TH1* h_unfold_err_ratio_b = (TH1*) histMunfold_b -> Clone( "h_unfold_err_ratio_b" ) ;
      TH1* h_unfold_err_ratio_c = (TH1*) histMunfold_c -> Clone( "h_unfold_err_ratio_c" ) ;

      h_unfold_err_ratio_b -> SetYTitle( "Unfolded error ratio" ) ;
      h_unfold_err_ratio_b -> SetTitle( "Unfolded error ratio" ) ;

      h_unfold_err_ratio_c -> SetYTitle( "Unfolded error ratio" ) ;
      h_unfold_err_ratio_c -> SetTitle( "Unfolded error ratio" ) ;

      h_unfold_err_b -> SetTitle( "Unfolded error" ) ;

      for ( int bi=1; bi<=histMunfold_a->GetNbinsX(); bi++ ) {
         float err_a = h_unfold_err_a -> GetBinError( bi ) ;
         float err_b = h_unfold_err_b -> GetBinError( bi ) ;
         float err_c = h_unfold_err_c -> GetBinError( bi ) ;
         h_unfold_err_a -> SetBinContent( bi, err_a ) ;
         h_unfold_err_b -> SetBinContent( bi, err_b ) ;
         h_unfold_err_c -> SetBinContent( bi, err_c ) ;
         h_unfold_err_a -> SetBinError( bi, 0. ) ;
         h_unfold_err_b -> SetBinError( bi, 0. ) ;
         h_unfold_err_c -> SetBinError( bi, 0. ) ;
         float ratio = 0  ;
         if ( err_b > 0 ) ratio = err_a / err_b ;
         h_unfold_err_ratio_b -> SetBinContent( bi, ratio ) ;
         h_unfold_err_ratio_b -> SetBinError( bi, 0. ) ;
         if ( err_c > 0 ) ratio = err_a / err_c ;
         h_unfold_err_ratio_c -> SetBinContent( bi, ratio ) ;
         h_unfold_err_ratio_c -> SetBinError( bi, 0. ) ;
      }

      h_unfold_err_ratio_b -> SetMaximum(1.1) ;
      h_unfold_err_ratio_c -> SetMaximum(1.1) ;

      histRhoi_a->SetMaximum(1.1) ;
      histRhoi_b->SetMaximum(1.1) ;
      histRhoi_c->SetMaximum(1.1) ;

      TExec* change_hist_palette = new TExec( "change_hist_palette", "Setup2DhistPalette();" );
      TExec* change_cor_palette = new TExec( "change_cor_palette", "SetupCorrelationPalette();" );


      gStyle -> SetPadRightMargin(0.08) ;
      gStyle -> SetPadLeftMargin(0.15) ;
      gStyle -> SetPadBottomMargin(0.15) ;
      gStyle -> SetTitleBorderSize(0) ;
      gStyle -> SetTitleY(0.975) ;



      int can_width = 600 ;
      int can_height = 600 ;

      int can_spacing = 20 ;


      float lx, ly, lh, lw ;

     //-----

      TCanvas* can1 = (TCanvas*) gDirectory -> FindObject( "can1" ) ;
      if ( can1 == 0x0 ) { printf( "Making can1\n") ; can1 = new TCanvas( "can1", "", can_spacing, can_spacing, can_width, can_height ) ; }
      can1 -> Clear() ;
      can1 -> cd() ;

      histRhoi_a -> SetLineColor(1) ;
      histRhoi_b -> SetLineColor(4) ;
      histRhoi_c -> SetLineColor(2) ;

      histRhoi_a -> SetLineWidth(3) ;
      histRhoi_b -> SetLineWidth(3) ;
      histRhoi_c -> SetLineWidth(3) ;

      histRhoi_a -> SetLineStyle(1) ;
      histRhoi_b -> SetLineStyle(9) ;
      histRhoi_c -> SetLineStyle(7) ;

      histRhoi_c -> SetYTitle("rho") ;

      histRhoi_c -> SetTitleOffset( 1.4, "y" ) ;
      histRhoi_c -> SetTitleOffset( 1.2, "x" ) ;

      histRhoi_c -> SetTitleSize( 0.045, "x" ) ;
      histRhoi_c -> SetTitleSize( 0.045, "y" ) ;

      if ( strcmp( var_name, "x" ) == 0 ) { histRhoi_c -> SetXTitle( "log10(x)" ) ; }
      if ( strcmp( var_name, "y" ) == 0 ) { histRhoi_c -> SetXTitle( "log10(y)" ) ; }
      if ( strcmp( var_name, "y" ) == 0 ) { histRhoi_c -> SetNdivisions( 605 ) ; }

      histRhoi_c -> Draw("hist") ;
      histRhoi_b -> Draw("hist same") ;
      histRhoi_a -> Draw("hist same") ;

      gPad -> SetGridy(1) ;

      if ( strcmp( var_name, "x" ) == 0 ) { lx = 0.20 ; ly = 0.75 ; lw = 0.25 ; lh = 0.13 ; }
      if ( strcmp( var_name, "y" ) == 0 ) { lx = 0.65 ; ly = 0.75 ; lw = 0.25 ; lh = 0.13 ; }

      TLegend* legend1 = new TLegend( lx, ly, lx+lw, ly+lh ) ;
      legend1 -> AddEntry( histRhoi_c, "electron" ) ;
      legend1 -> AddEntry( histRhoi_b, "Sigma" ) ;
      legend1 -> AddEntry( histRhoi_a, "DNN" ) ;

      legend1 -> Draw() ;

      can1 -> Update() ; can1 -> Draw() ; gSystem -> ProcessEvents() ;

      sprintf( fname, "paper-plots/dis-global-correlation-comparison-%s.pdf", var_name ) ;
      can1 -> SaveAs( fname ) ;
      sprintf( fname, "paper-plots/dis-global-correlation-comparison-%s.png", var_name ) ;
      can1 -> SaveAs( fname ) ;





     //-----

      TCanvas* can2 = (TCanvas*) gDirectory -> FindObject( "can2" ) ;
      if ( can2 == 0x0 ) { printf( "Making can2\n") ; can2 = new TCanvas( "can2", "", 2*can_spacing + can_width , can_spacing, can_width, can_height ) ; }
      can2 -> Clear() ;
      can2 -> cd() ;

      h_unfold_err_ratio_b -> SetLineWidth(3) ;
      h_unfold_err_ratio_c -> SetLineWidth(3) ;

      h_unfold_err_ratio_b -> SetLineColor(4) ;
      h_unfold_err_ratio_c -> SetLineColor(2) ;

      h_unfold_err_ratio_b -> SetLineStyle(9) ;
      h_unfold_err_ratio_c -> SetLineStyle(7) ;


      h_unfold_err_ratio_b -> SetTitleOffset( 1.4, "y" ) ;
      h_unfold_err_ratio_b -> SetTitleOffset( 1.2, "x" ) ;
      h_unfold_err_ratio_b -> SetTitleSize( 0.045, "x" ) ;
      h_unfold_err_ratio_b -> SetTitleSize( 0.045, "y" ) ;

      h_unfold_err_ratio_b -> SetYTitle( "Error ratio" ) ;

      if ( strcmp( var_name, "x" ) == 0 ) { h_unfold_err_ratio_b -> SetXTitle( "log10(x)" ) ; }
      if ( strcmp( var_name, "y" ) == 0 ) { h_unfold_err_ratio_b -> SetXTitle( "log10(y)" ) ; }
      if ( strcmp( var_name, "y" ) == 0 ) { h_unfold_err_ratio_b -> SetNdivisions( 605 ) ; }

      h_unfold_err_ratio_b -> Draw("hist") ;
      h_unfold_err_ratio_c -> Draw("hist same") ;

      gPad -> SetGridy(1) ;

      if ( strcmp( var_name, "x" ) == 0 ) { lx = 0.55 ; ly = 0.75 ; lw = 0.35 ; lh = 0.13 ; }
      if ( strcmp( var_name, "y" ) == 0 ) { lx = 0.20 ; ly = 0.75 ; lw = 0.35 ; lh = 0.13 ; }

      TLegend* legend2 = new TLegend( lx, ly, lx+lw, ly+lh ) ;
      legend2 -> AddEntry( h_unfold_err_ratio_c, "DNN / electron" ) ;
      legend2 -> AddEntry( h_unfold_err_ratio_b, "DNN / Sigma" ) ;

      legend2 -> Draw() ;

      can2 -> Update() ; can2 -> Draw() ; gSystem -> ProcessEvents() ;

      sprintf( fname, "paper-plots/dis-error-ratio-comparison-%s.pdf", var_name ) ;
      can2 -> SaveAs( fname ) ;
      sprintf( fname, "paper-plots/dis-error-ratio-comparison-%s.png", var_name ) ;
      can2 -> SaveAs( fname ) ;





     //-----

      TCanvas* can3 = (TCanvas*) gDirectory -> FindObject( "can3" ) ;
      if ( can3 == 0x0 ) { printf( "Making can3\n") ; can3 = new TCanvas( "can3", "", 3*can_spacing + 2*can_width , can_spacing, can_width, can_height ) ; }
      can3 -> Clear() ;
      can3 -> cd() ;


      h_unfold_err_a -> SetLineColor(1) ;
      h_unfold_err_b -> SetLineColor(4) ;
      h_unfold_err_c -> SetLineColor(2) ;

      h_unfold_err_a -> SetLineWidth(3) ;
      h_unfold_err_b -> SetLineWidth(3) ;
      h_unfold_err_c -> SetLineWidth(3) ;

      h_unfold_err_a -> SetLineStyle(1) ;
      h_unfold_err_b -> SetLineStyle(9) ;
      h_unfold_err_c -> SetLineStyle(7) ;


      h_unfold_err_c -> SetTitle( "Unfolded error" ) ;
      h_unfold_err_c -> SetYTitle( "Unfolded error (events)" ) ;
      h_unfold_err_c -> SetTitleOffset( 1.4, "y" ) ;
      h_unfold_err_c -> SetTitleOffset( 1.2, "x" ) ;
      h_unfold_err_c -> SetTitleSize( 0.045, "x" ) ;
      h_unfold_err_c -> SetTitleSize( 0.045, "y" ) ;

      if ( strcmp( var_name, "x" ) == 0 ) { h_unfold_err_c -> SetXTitle( "log10(x)" ) ; }
      if ( strcmp( var_name, "y" ) == 0 ) { h_unfold_err_c -> SetXTitle( "log10(y)" ) ; }
      if ( strcmp( var_name, "y" ) == 0 ) { h_unfold_err_c -> SetNdivisions( 605 ) ; }

      if ( max_error > 0 ) h_unfold_err_c -> SetMaximum( max_error ) ;


      h_unfold_err_c -> Draw("hist") ;
      h_unfold_err_b -> Draw("hist same") ;
      h_unfold_err_a -> Draw("hist same") ;

      gPad -> SetGridy(1) ;

      if ( strcmp( var_name, "x" ) == 0 ) { lx = 0.20 ; ly = 0.75 ; lw = 0.25 ; lh = 0.13 ; }
      if ( strcmp( var_name, "y" ) == 0 ) { lx = 0.65 ; ly = 0.75 ; lw = 0.25 ; lh = 0.13 ; }

      TLegend* legend3 = new TLegend( lx, ly, lx+lw, ly+lh ) ;
      legend3 -> AddEntry( h_unfold_err_c, "electron" ) ;
      legend3 -> AddEntry( h_unfold_err_b, "Sigma" ) ;
      legend3 -> AddEntry( h_unfold_err_a, "DNN" ) ;

      legend3 -> Draw() ;

      can3 -> Update() ; can3 -> Draw() ; gSystem -> ProcessEvents() ;

      sprintf( fname, "paper-plots/dis-error-comparison-%s.pdf", var_name ) ;
      can3 -> SaveAs( fname ) ;
      sprintf( fname, "paper-plots/dis-error-comparison-%s.png", var_name ) ;
      can3 -> SaveAs( fname ) ;



      printf("\n\n\n") ;

      printf(" cut and paste for this:\n\n") ;

      printf("     paper_plots_dis_compare3(\"%s\",\"%s\",\"%s\",\"%s\",%d,\"%s\")\n",
          hist_name_a, hist_name_b, hist_name_c, var_name, ngen, input_file ) ;


      printf("\n\n\n") ;

   }













