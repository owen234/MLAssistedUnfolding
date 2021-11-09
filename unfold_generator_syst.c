

//
//  Following the example of tutorials/unfold/testUnfold1.c
//


#include "histio.c"
#include "draw_obs_in_gen_slices.c"

TH2F* get_hist( const char* hname ) {
   TH2F* hp = (TH2F*) gDirectory -> FindObject( hname ) ;
   if ( hp == 0x0 ) { printf("\n\n *** can't find %s\n\n", hname ) ; gSystem->Exit(-1) ; }
   return hp ;
}

//----------

   void unfold_generator_syst( const char* hist_name = "h_log10_x_gen_vs_obs_dnn",
                      const char* input_file_a = "h1-input-nbins_gen020_obs050-b2d-rapgap.root",
                      const char* input_file_b = "h1-input-nbins_gen020_obs050-b2d-django.root",
                      const char* input_name_a = "Rapgap",
                      const char* input_name_b = "Djangoh",
                      const char* method_name = "DNN",
                      int ngen = 1e7,
                      float diff_max = -1
                      ) {

      gSystem -> Exec( "mkdir -p plots" ) ;

      TRandom3 tran(12345) ;

      gDirectory -> Delete( "h*" ) ;
      gStyle -> SetPalette( kBird ) ;

      loadHist( input_file_a, input_name_a ) ;
      loadHist( input_file_b, input_name_b ) ;


      printf("\n\n") ;
      gDirectory -> ls() ;
      printf("\n\n") ;

      char hname[1000] ;
      sprintf( hname, "%s_%s", hist_name, input_name_a ) ;
      TH2F* h_in_gen_vs_obs_a = get_hist( hname ) ;

      sprintf( hname, "%s_%s", hist_name, input_name_b ) ;
      TH2F* h_in_gen_vs_obs_b = get_hist( hname ) ;


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

      /////underflow_mean = ngen * ( nobs_underflow / nobs_source ) ;
      /////overflow_mean = ngen * ( nobs_overflow / nobs_source ) ;

      underflow_mean = ngen * ( nobs_underflow / nobs_integral ) ; // ngen only fills inside of histogram (integral)
      overflow_mean = ngen * ( nobs_overflow / nobs_integral ) ;

      h_obs_random_a -> SetBinContent( 0, tran.Poisson( underflow_mean ) ) ;
      h_obs_random_a -> SetBinContent( nbins_obs+1, tran.Poisson( overflow_mean ) ) ;



      nobs_integral = h_obs_source_b -> Integral() ;
      nobs_underflow = h_obs_source_b -> GetBinContent(0) ;
      nobs_overflow = h_obs_source_b -> GetBinContent( nbins_obs+1 ) ;

      printf("  h_obs_source_b :  integral %9.5f, underflow %9.5f, overflow %9.5f\n", nobs_integral, nobs_underflow, nobs_overflow ) ;

      nobs_source = nobs_integral ;
      nobs_source += nobs_underflow ;
      nobs_source += nobs_overflow ;


      //////underflow_mean = ngen * ( nobs_underflow / nobs_source ) ;
      //////overflow_mean = ngen * ( nobs_overflow / nobs_source ) ;

      underflow_mean = ngen * ( nobs_underflow / nobs_integral ) ; // ngen only fills inside of histogram (integral)
      overflow_mean = ngen * ( nobs_overflow / nobs_integral ) ;


      h_obs_random_b -> SetBinContent( 0, tran.Poisson( underflow_mean ) ) ;
      h_obs_random_b -> SetBinContent( nbins_obs+1, tran.Poisson( overflow_mean ) ) ;








      TUnfoldDensity unfold_a_with_a( h_in_gen_vs_obs_a, TUnfold::kHistMapOutputVert ) ;
      TUnfoldDensity unfold_a_with_b( h_in_gen_vs_obs_b, TUnfold::kHistMapOutputVert ) ;
      TUnfoldDensity unfold_b_with_a( h_in_gen_vs_obs_a, TUnfold::kHistMapOutputVert ) ;
      TUnfoldDensity unfold_b_with_b( h_in_gen_vs_obs_b, TUnfold::kHistMapOutputVert ) ;

      int return_status ;

      return_status = unfold_a_with_a.SetInput( h_obs_random_a ) ;
      printf("  Return status for SetInput A with A: %d\n", return_status ) ;

      return_status = unfold_a_with_b.SetInput( h_obs_random_a ) ;
      printf("  Return status for SetInput A with B: %d\n", return_status ) ;

      return_status = unfold_b_with_a.SetInput( h_obs_random_b ) ;
      printf("  Return status for SetInput B with A: %d\n", return_status ) ;

      return_status = unfold_b_with_b.SetInput( h_obs_random_b ) ;
      printf("  Return status for SetInput B with B: %d\n", return_status ) ;





      //========================================================================
      // the unfolding is done here
      //
      // scan L curve and find best point
      Int_t nScan=30;
      // use automatic L-curve scan: start with taumin=taumax=0.0
      Double_t tauMin=0.0;
      Double_t tauMax=0.0;
      Int_t iBest;
      TSpline *logTauX_a_with_a,*logTauY_a_with_a;
      TSpline *logTauX_a_with_b,*logTauY_a_with_b;
      TSpline *logTauX_b_with_a,*logTauY_b_with_a;
      TSpline *logTauX_b_with_b,*logTauY_b_with_b;
      TGraph *lCurve_a_with_a;
      TGraph *lCurve_a_with_b;
      TGraph *lCurve_b_with_a;
      TGraph *lCurve_b_with_b;

      // if required, report Info messages (for debugging the L-curve scan)
      Int_t oldinfo=gErrorIgnoreLevel;
      gErrorIgnoreLevel=kInfo;

      // this method scans the parameter tau and finds the kink in the L curve
      // finally, the unfolding is done for the best choice of tau

      printf("\n\n ===== Begin unfolding %s with %s\n", input_name_a, input_name_a ) ;
      iBest=unfold_a_with_a.ScanLcurve( nScan, tauMin, tauMax, &lCurve_a_with_a, &logTauX_a_with_a, &logTauY_a_with_a );

      printf("\n\n ===== Begin unfolding %s with %s\n", input_name_a, input_name_b ) ;
      iBest=unfold_a_with_b.ScanLcurve( nScan, tauMin, tauMax, &lCurve_a_with_b, &logTauX_a_with_b, &logTauY_a_with_b );


      printf("\n\n ===== Begin unfolding %s with %s\n", input_name_b, input_name_a ) ;
      iBest=unfold_b_with_a.ScanLcurve( nScan, tauMin, tauMax, &lCurve_b_with_a, &logTauX_b_with_a, &logTauY_b_with_a );

      printf("\n\n ===== Begin unfolding %s with %s\n", input_name_b, input_name_b ) ;
      iBest=unfold_b_with_b.ScanLcurve( nScan, tauMin, tauMax, &lCurve_b_with_b, &logTauX_b_with_b, &logTauY_b_with_b );


      // if required, switch to previous log-level
      gErrorIgnoreLevel=oldinfo;


      //==========================================================================
      // print some results
      //
      printf("\n\n %11s with %11s : tau= %9.7f ,  chi2= (%9.5f + %9.5f)/%d\n",
         input_name_a, input_name_a, unfold_a_with_a.GetTau(), unfold_a_with_a.GetChi2A(), unfold_a_with_a.GetChi2L(), unfold_a_with_a.GetNdf() ) ;

      printf("\n\n %11s with %11s : tau= %9.7f ,  chi2= (%9.5f + %9.5f)/%d\n",
         input_name_a, input_name_b, unfold_a_with_b.GetTau(), unfold_a_with_b.GetChi2A(), unfold_a_with_b.GetChi2L(), unfold_a_with_b.GetNdf() ) ;

      printf("\n\n %11s with %11s : tau= %9.7f ,  chi2= (%9.5f + %9.5f)/%d\n",
         input_name_b, input_name_a, unfold_b_with_a.GetTau(), unfold_b_with_a.GetChi2A(), unfold_b_with_a.GetChi2L(), unfold_b_with_a.GetNdf() ) ;

      printf("\n\n %11s with %11s : tau= %9.7f ,  chi2= (%9.5f + %9.5f)/%d\n",
         input_name_b, input_name_b, unfold_b_with_b.GetTau(), unfold_b_with_b.GetChi2A(), unfold_b_with_b.GetChi2L(), unfold_b_with_b.GetNdf() ) ;

      printf("\n\n") ;

      //==========================================================================
      // retrieve results into histograms

      // get unfolded distribution
      TH1 *histMunfold_a_with_a = unfold_a_with_a.GetOutput("Unfolded_a_with_a");
      TH1 *histMunfold_a_with_b = unfold_a_with_b.GetOutput("Unfolded_a_with_b");
      TH1 *histMunfold_b_with_a = unfold_b_with_a.GetOutput("Unfolded_b_with_a");
      TH1 *histMunfold_b_with_b = unfold_b_with_b.GetOutput("Unfolded_b_with_b");

      char htitle[1000] ;

      sprintf( htitle, "Unfolded %s with %s", input_name_a, input_name_a ) ;
      histMunfold_a_with_a -> SetTitle( htitle ) ;
      sprintf( htitle, "Unfolded %s with %s", input_name_a, input_name_b ) ;
      histMunfold_a_with_b -> SetTitle( htitle ) ;
      sprintf( htitle, "Unfolded %s with %s", input_name_b, input_name_a ) ;
      histMunfold_b_with_a -> SetTitle( htitle ) ;
      sprintf( htitle, "Unfolded %s with %s", input_name_b, input_name_b ) ;
      histMunfold_b_with_b -> SetTitle( htitle ) ;




      TH1* h_gen_compare_a_with_a = (TH1*) h_gen_source_a -> Clone( "h_gen_compare_a_with_a" ) ;
      h_gen_compare_a_with_a -> Scale( ( histMunfold_a_with_a -> Integral() )/( h_gen_source_a -> Integral() ) ) ;

      TH1* h_gen_compare_a_with_b = (TH1*) h_gen_source_a -> Clone( "h_gen_compare_a_with_b" ) ;
      h_gen_compare_a_with_b -> Scale( ( histMunfold_a_with_b -> Integral() )/( h_gen_source_a -> Integral() ) ) ;


      TH1* h_gen_compare_b_with_a = (TH1*) h_gen_source_b -> Clone( "h_gen_compare_b_with_a" ) ;
      h_gen_compare_b_with_a -> Scale( ( histMunfold_b_with_a -> Integral() )/( h_gen_source_b -> Integral() ) ) ;

      TH1* h_gen_compare_b_with_b = (TH1*) h_gen_source_b -> Clone( "h_gen_compare_b_with_b" ) ;
      h_gen_compare_b_with_b -> Scale( ( histMunfold_b_with_b -> Integral() )/( h_gen_source_b -> Integral() ) ) ;


      h_gen_compare_a_with_a -> SetLineColor(2) ;
      h_gen_compare_a_with_b -> SetLineColor(2) ;
      h_gen_compare_b_with_a -> SetLineColor(2) ;
      h_gen_compare_b_with_b -> SetLineColor(2) ;


      TH1* h_diff_a_with_a = (TH1*) histMunfold_a_with_a -> Clone( "h_diff_a_with_a" ) ;
      TH1* h_diff_a_with_b = (TH1*) histMunfold_a_with_b -> Clone( "h_diff_a_with_b" ) ;
      TH1* h_diff_b_with_a = (TH1*) histMunfold_b_with_a -> Clone( "h_diff_b_with_a" ) ;
      TH1* h_diff_b_with_b = (TH1*) histMunfold_b_with_b -> Clone( "h_diff_b_with_b" ) ;

      h_diff_a_with_a -> Add( h_gen_compare_a_with_a, -1. ) ;
      h_diff_a_with_b -> Add( h_gen_compare_a_with_a, -1. ) ;
      h_diff_b_with_a -> Add( h_gen_compare_a_with_a, -1. ) ;
      h_diff_b_with_b -> Add( h_gen_compare_a_with_a, -1. ) ;

      float max_abs_diff = 0. ;
      if ( h_diff_a_with_a -> GetMaximum() > max_abs_diff ) max_abs_diff = h_diff_a_with_a -> GetMaximum() ;
      if ( h_diff_a_with_b -> GetMaximum() > max_abs_diff ) max_abs_diff = h_diff_a_with_b -> GetMaximum() ;
      if ( h_diff_b_with_a -> GetMaximum() > max_abs_diff ) max_abs_diff = h_diff_b_with_a -> GetMaximum() ;
      if ( h_diff_b_with_b -> GetMaximum() > max_abs_diff ) max_abs_diff = h_diff_b_with_b -> GetMaximum() ;
      if ( fabs( h_diff_a_with_a -> GetMinimum() ) > max_abs_diff ) max_abs_diff = fabs( h_diff_a_with_a -> GetMinimum() ) ;
      if ( fabs( h_diff_a_with_b -> GetMinimum() ) > max_abs_diff ) max_abs_diff = fabs( h_diff_a_with_b -> GetMinimum() ) ;
      if ( fabs( h_diff_b_with_a -> GetMinimum() ) > max_abs_diff ) max_abs_diff = fabs( h_diff_b_with_a -> GetMinimum() ) ;
      if ( fabs( h_diff_b_with_b -> GetMinimum() ) > max_abs_diff ) max_abs_diff = fabs( h_diff_b_with_b -> GetMinimum() ) ;

      if ( diff_max > 0 && diff_max > max_abs_diff ) max_abs_diff = diff_max ;
      //max_abs_diff = diff_max ;

      h_diff_a_with_a -> SetMaximum( 1.1 * max_abs_diff ) ;
      h_diff_a_with_b -> SetMaximum( 1.1 * max_abs_diff ) ;
      h_diff_b_with_a -> SetMaximum( 1.1 * max_abs_diff ) ;
      h_diff_b_with_b -> SetMaximum( 1.1 * max_abs_diff ) ;

      h_diff_a_with_a -> SetMinimum( -1.1 * max_abs_diff ) ;
      h_diff_a_with_b -> SetMinimum( -1.1 * max_abs_diff ) ;
      h_diff_b_with_a -> SetMinimum( -1.1 * max_abs_diff ) ;
      h_diff_b_with_b -> SetMinimum( -1.1 * max_abs_diff ) ;



      h_diff_a_with_a -> SetLineColor(2) ;
      h_diff_a_with_b -> SetLineColor(2) ;
      h_diff_b_with_a -> SetLineColor(4) ;
      h_diff_b_with_b -> SetLineColor(4) ;

      h_diff_a_with_a -> SetLineWidth(2) ;
      h_diff_a_with_b -> SetLineWidth(2) ;
      h_diff_b_with_a -> SetLineWidth(2) ;
      h_diff_b_with_b -> SetLineWidth(2) ;




      TLine* tl = new TLine() ;
      tl -> SetLineColor(2) ;

      float lx = 0.25 ;
      float ly = 0.15 ;
      float lw = 0.50 ;
      float lh = 0.15 ;

      float lw2 = 0.30 ;
      float lh2 = 0.15 ;

      char label[100] ;

      TLegend* legend ;


      gStyle -> SetOptStat(0) ;

      TCanvas* can1 = (TCanvas*) gDirectory -> FindObject( "can1" ) ;
      if ( can1 == 0x0 ) can1 = new TCanvas( "can1", "", 50, 50, 1300, 1100 ) ;
      can1 -> Clear() ;
      can1 -> cd() ;

      can1 -> Divide(2,3) ;

      int ci(1) ;



      can1 -> cd( ci++ ) ;
      histMunfold_a_with_a -> SetTitle( method_name ) ;
      histMunfold_a_with_a -> Draw() ;
      h_gen_compare_a_with_a -> Draw("same hist") ;

      legend = new TLegend( lx, ly, lx+lw, ly+lh ) ;
      sprintf( label, "%s unfolded with %s", input_name_a, input_name_a ) ;
      legend -> AddEntry( histMunfold_a_with_a, label ) ;
      sprintf( label, "Gen %s", input_name_a ) ;
      legend -> AddEntry( h_gen_compare_a_with_a, label ) ;
      legend -> Draw() ;



      can1 -> cd( ci++ ) ;
      histMunfold_a_with_b -> SetTitle( method_name ) ;
      histMunfold_a_with_b -> Draw() ;
      h_gen_compare_a_with_b -> Draw("same hist") ;

      legend = new TLegend( lx, ly, lx+lw, ly+lh ) ;
      sprintf( label, "%s unfolded with %s", input_name_a, input_name_b ) ;
      legend -> AddEntry( histMunfold_a_with_b, label ) ;
      sprintf( label, "Gen %s", input_name_a ) ;
      legend -> AddEntry( h_gen_compare_a_with_a, label ) ;
      legend -> Draw() ;



      can1 -> cd( ci++ ) ;
      histMunfold_b_with_b -> SetTitle( method_name ) ;
      histMunfold_b_with_b -> Draw() ;
      h_gen_compare_b_with_b -> Draw("same hist") ;

      legend = new TLegend( lx, ly, lx+lw, ly+lh ) ;
      sprintf( label, "%s unfolded with %s", input_name_b, input_name_b ) ;
      legend -> AddEntry( histMunfold_b_with_b, label ) ;
      sprintf( label, "Gen %s", input_name_b ) ;
      legend -> AddEntry( h_gen_compare_b_with_b, label ) ;
      legend -> Draw() ;



      can1 -> cd( ci++ ) ;
      histMunfold_b_with_a -> SetTitle( method_name ) ;
      histMunfold_b_with_a -> Draw() ;
      h_gen_compare_b_with_a -> Draw("same hist") ;

      legend = new TLegend( lx, ly, lx+lw, ly+lh ) ;
      sprintf( label, "%s unfolded with %s", input_name_b, input_name_a ) ;
      legend -> AddEntry( histMunfold_b_with_a, label ) ;
      sprintf( label, "Gen %s", input_name_b ) ;
      legend -> AddEntry( h_gen_compare_b_with_a, label ) ;
      legend -> Draw() ;





      can1 -> cd( ci++ ) ;
      h_diff_a_with_a -> SetTitle("Closure tests") ;
      h_diff_a_with_a -> Draw("hist") ;
      h_diff_b_with_b -> Draw("hist same") ;
      gPad -> SetGridy(1) ;
      //tl -> DrawLine( h_diff_a_with_a -> GetXaxis() -> GetXmin(), 0.,   h_diff_a_with_a -> GetXaxis() -> GetXmax(), 0. ) ;

      legend = new TLegend( lx, ly, lx+lw2, ly+lh2 ) ;
      sprintf( label, "%s", input_name_a ) ;
      legend -> AddEntry( h_diff_a_with_a, label ) ;
      sprintf( label, "%s", input_name_b ) ;
      legend -> AddEntry( h_diff_b_with_b, label ) ;
      legend -> Draw() ;





      can1 -> cd( ci++ ) ;
      h_diff_b_with_a -> SetTitle("Generator systematic") ;
      h_diff_b_with_a -> Draw("hist") ;
      h_diff_a_with_b -> Draw("hist same") ;
      gPad -> SetGridy(1) ;
      //tl -> DrawLine( h_diff_b_with_a -> GetXaxis() -> GetXmin(), 0.,   h_diff_b_with_a -> GetXaxis() -> GetXmax(), 0. ) ;

      legend = new TLegend( lx, ly, lx+lw, ly+lh ) ;
      sprintf( label, "%s unfolded with %s", input_name_b, input_name_a ) ;
      legend -> AddEntry( h_diff_b_with_a, label ) ;
      sprintf( label, "%s unfolded with %s", input_name_a, input_name_b ) ;
      legend -> AddEntry( h_diff_a_with_b, label ) ;
      legend -> Draw() ;



      char fname[1000] ;
      sprintf( fname, "plots/generator-syst-%s.png", method_name ) ;
      can1 -> SaveAs( fname ) ;
      sprintf( fname, "plots/generator-syst-%s.pdf", method_name ) ;
      can1 -> SaveAs( fname ) ;


      printf("\n\n\n") ;
      printf("  Cut and paste for this:\n") ;
      printf(" unfold_generator_syst(\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%d,%9.2f)",
         hist_name, input_file_a, input_file_b, input_name_a, input_name_b, method_name, ngen, diff_max ) ;
      printf("\n\n\n") ;

      return ;



   }













