

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

#endif

//----------

   void paper_plots_dis_gen_syst(  const char* var_name = "x",
                      int seed = 12345,
                      float diff_max = -1,
                      const char* input_file_a = "paper-plots-input-1D-nbins_gen010_obs020.root",
                      const char* input_file_b = "paper-plots-input-1D-nbins_gen010_obs020-django.root",
                      const char* input_name_a = "Rapgap",
                      const char* input_name_b = "Djangoh",
                      int ngen = 1e7
                      ) {

      char htitle[100] ;
      char fname[1000] ;

      gSystem -> Exec( "mkdir -p paper-plots" ) ;

      char response_hname_base[100] ;
      if ( strcmp( var_name, "x" ) == 0 ) { sprintf( response_hname_base, "h_log10_x_gen_vs_obs" ) ; }
      if ( strcmp( var_name, "y" ) == 0 ) { sprintf( response_hname_base, "h_log10_y_gen_vs_obs" ) ; }


      TRandom3* tran = new TRandom3( seed ) ;

      gDirectory -> Delete( "h*" ) ;
      gStyle -> SetPalette( kBird ) ;

      loadHist( input_file_a, input_name_a ) ;
      loadHist( input_file_b, input_name_b ) ;

      char methods[3][100] ;
      sprintf( methods[0], "e" ) ;
      sprintf( methods[1], "esigma" ) ;
      sprintf( methods[2], "dnn" ) ;

      char methods_title[3][100] ;
      sprintf( methods_title[0], "electron" ) ;
      sprintf( methods_title[1], "Sigma" ) ;
      sprintf( methods_title[2], "DNN" ) ;

      int method_lc[3] ;
      method_lc[0] = 2 ;
      method_lc[1] = 4 ;
      method_lc[2] = 1 ;

      int method_ls[3] ;
      method_ls[0] = 7 ;
      method_ls[1] = 9 ;
      method_ls[2] = 1 ;

      printf("\n\n") ;
      gDirectory -> ls() ;
      printf("\n\n") ;

      char hname[1000] ;

      TH2F* h_response_a[3] ;
      TH2F* h_response_b[3] ;

      TH1*  h_obs_source_a[3] ;
      TH1*  h_obs_source_b[3] ;

      TH1*  h_gen_source_a[3] ;
      TH1*  h_gen_source_b[3] ;

      TH1*  h_obs_random_a[3] ;
      TH1*  h_obs_random_b[3] ;

      TUnfoldDensity* tud_a_with_a[3] ;
      TUnfoldDensity* tud_a_with_b[3] ;
      TUnfoldDensity* tud_b_with_a[3] ;
      TUnfoldDensity* tud_b_with_b[3] ;

      TH1*  h_unfolded_a_with_a[3] ;
      TH1*  h_unfolded_a_with_b[3] ;
      TH1*  h_unfolded_b_with_a[3] ;
      TH1*  h_unfolded_b_with_b[3] ;

      TH1* h_gen_compare_a_with_a[3] ;
      TH1* h_gen_compare_a_with_b[3] ;
      TH1* h_gen_compare_b_with_a[3] ;
      TH1* h_gen_compare_b_with_b[3] ;

      TH1* h_diff_a_with_a[3] ;
      TH1* h_diff_a_with_b[3] ;
      TH1* h_diff_b_with_a[3] ;
      TH1* h_diff_b_with_b[3] ;

      TH1* h_abs_diff_a_with_a[3] ;
      TH1* h_abs_diff_a_with_b[3] ;
      TH1* h_abs_diff_b_with_a[3] ;
      TH1* h_abs_diff_b_with_b[3] ;

      TH1* h_ave_abs_diff_not_same[3] ;


      for ( int mi=0; mi<3; mi++ ) {

         sprintf( hname, "%s_%s_%s", response_hname_base, methods[mi], input_name_a ) ;
         h_response_a[mi] = get_hist( hname ) ;
         sprintf( hname, "%s_%s_%s", response_hname_base, methods[mi], input_name_b ) ;
         h_response_b[mi] = get_hist( hname ) ;

         sprintf( hname, "h_obs_source_%s_a", methods[mi] ) ;
         h_obs_source_a[mi] = h_response_a[mi] -> ProjectionX( hname ) ;
         sprintf( hname, "h_obs_source_%s_b", methods[mi] ) ;
         h_obs_source_b[mi] = h_response_b[mi] -> ProjectionX( hname ) ;

         sprintf( hname, "h_gen_source_%s_a", methods[mi] ) ;
         h_gen_source_a[mi] = h_response_a[mi] -> ProjectionY( hname ) ;
         sprintf( hname, "h_gen_source_%s_b", methods[mi] ) ;
         h_gen_source_b[mi] = h_response_b[mi] -> ProjectionY( hname ) ;

         sprintf( hname, "h_obs_random_%s_a", methods[mi] ) ;
         h_obs_random_a[mi] = (TH1*) h_obs_source_a[mi] -> Clone( hname ) ;
         h_obs_random_a[mi] -> Reset() ;
         h_obs_random_a[mi] -> FillRandom( h_obs_source_a[mi], ngen, tran ) ;

         sprintf( hname, "h_obs_random_%s_b", methods[mi] ) ;
         h_obs_random_b[mi] = (TH1*) h_obs_source_b[mi] -> Clone( hname ) ;
         h_obs_random_b[mi] -> Reset() ;
         h_obs_random_b[mi] -> FillRandom( h_obs_source_b[mi], ngen, tran ) ;





        //-- Add underflows and overflows to fake data!
         float underflow_mean ; float overflow_mean ;
         int nbins_obs = h_obs_source_a[mi] -> GetNbinsX() ;
         float nobs_source ; float nobs_integral ; float nobs_underflow ; float nobs_overflow ;



         nobs_integral = h_obs_source_a[mi] -> Integral() ;
         nobs_underflow = h_obs_source_a[mi] -> GetBinContent(0) ;
         nobs_overflow = h_obs_source_a[mi] -> GetBinContent( nbins_obs+1 ) ;

         printf("  h_obs_source_a :  integral %9.5f, underflow %9.5f, overflow %9.5f\n", nobs_integral, nobs_underflow, nobs_overflow ) ;

         nobs_source = nobs_integral ;
         nobs_source += nobs_underflow ;
         nobs_source += nobs_overflow ;

         underflow_mean = ngen * ( nobs_underflow / nobs_integral ) ; // ngen only fills inside of histogram (integral)
         overflow_mean = ngen * ( nobs_overflow / nobs_integral ) ;

         h_obs_random_a[mi] -> SetBinContent( 0, tran -> Poisson( underflow_mean ) ) ;
         h_obs_random_a[mi] -> SetBinContent( nbins_obs+1, tran -> Poisson( overflow_mean ) ) ;



         nobs_integral = h_obs_source_b[mi] -> Integral() ;
         nobs_underflow = h_obs_source_b[mi] -> GetBinContent(0) ;
         nobs_overflow = h_obs_source_b[mi] -> GetBinContent( nbins_obs+1 ) ;

         printf("  h_obs_source_b :  integral %9.5f, underflow %9.5f, overflow %9.5f\n", nobs_integral, nobs_underflow, nobs_overflow ) ;

         nobs_source = nobs_integral ;
         nobs_source += nobs_underflow ;
         nobs_source += nobs_overflow ;

         underflow_mean = ngen * ( nobs_underflow / nobs_integral ) ; // ngen only fills inside of histogram (integral)
         overflow_mean = ngen * ( nobs_overflow / nobs_integral ) ;

         h_obs_random_b[mi] -> SetBinContent( 0, tran -> Poisson( underflow_mean ) ) ;
         h_obs_random_b[mi] -> SetBinContent( nbins_obs+1, tran -> Poisson( overflow_mean ) ) ;



        //--- set up and run unfolding

         tud_a_with_a[mi] = new TUnfoldDensity( h_response_a[mi] , TUnfold::kHistMapOutputVert ) ;
         tud_a_with_b[mi] = new TUnfoldDensity( h_response_b[mi] , TUnfold::kHistMapOutputVert ) ;
         tud_b_with_a[mi] = new TUnfoldDensity( h_response_a[mi] , TUnfold::kHistMapOutputVert ) ;
         tud_b_with_b[mi] = new TUnfoldDensity( h_response_b[mi] , TUnfold::kHistMapOutputVert ) ;

         int return_status ;

         return_status = tud_a_with_a[mi] -> SetInput( h_obs_random_a[mi] ) ;
         printf(" %s :  Return status for SetInput A with A: %d\n", methods[mi], return_status ) ;

         return_status = tud_a_with_b[mi] -> SetInput( h_obs_random_a[mi] ) ;
         printf(" %s :  Return status for SetInput A with B: %d\n", methods[mi], return_status ) ;

         return_status = tud_b_with_a[mi] -> SetInput( h_obs_random_b[mi] ) ;
         printf(" %s :  Return status for SetInput B with A: %d\n", methods[mi], return_status ) ;

         return_status = tud_b_with_b[mi] -> SetInput( h_obs_random_b[mi] ) ;
         printf(" %s :  Return status for SetInput B with B: %d\n", methods[mi], return_status ) ;


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

         printf("\n\n ===== %s : Begin unfolding %s with %s\n", methods[mi], input_name_a, input_name_a ) ;
         iBest = tud_a_with_a[mi] -> ScanLcurve( nScan, tauMin, tauMax, &lCurve_a_with_a, &logTauX_a_with_a, &logTauY_a_with_a );

         printf("\n\n ===== %s : Begin unfolding %s with %s\n", methods[mi], input_name_a, input_name_b ) ;
         iBest = tud_a_with_b[mi] -> ScanLcurve( nScan, tauMin, tauMax, &lCurve_a_with_b, &logTauX_a_with_b, &logTauY_a_with_b );


         printf("\n\n ===== %s : Begin unfolding %s with %s\n", methods[mi], input_name_b, input_name_a ) ;
         iBest = tud_b_with_a[mi] -> ScanLcurve( nScan, tauMin, tauMax, &lCurve_b_with_a, &logTauX_b_with_a, &logTauY_b_with_a );

         printf("\n\n ===== %s : Begin unfolding %s with %s\n", methods[mi], input_name_b, input_name_b ) ;
         iBest = tud_b_with_b[mi] -> ScanLcurve( nScan, tauMin, tauMax, &lCurve_b_with_b, &logTauX_b_with_b, &logTauY_b_with_b );


         // if required, switch to previous log-level
         gErrorIgnoreLevel=oldinfo;


         //==========================================================================
         // print some results
         //
         printf("\n\n %s : %11s with %11s : tau= %9.7f ,  chi2= (%9.5f + %9.5f)/%d\n",
            methods[mi], input_name_a, input_name_a, tud_a_with_a[mi] -> GetTau(), tud_a_with_a[mi] -> GetChi2A(), tud_a_with_a[mi] -> GetChi2L(), tud_a_with_a[mi] -> GetNdf() ) ;

         printf("\n\n %s : %11s with %11s : tau= %9.7f ,  chi2= (%9.5f + %9.5f)/%d\n",
            methods[mi], input_name_a, input_name_b, tud_a_with_b[mi] -> GetTau(), tud_a_with_b[mi] -> GetChi2A(), tud_a_with_b[mi] -> GetChi2L(), tud_a_with_b[mi] -> GetNdf() ) ;

         printf("\n\n %s : %11s with %11s : tau= %9.7f ,  chi2= (%9.5f + %9.5f)/%d\n",
            methods[mi], input_name_b, input_name_a, tud_b_with_a[mi] -> GetTau(), tud_b_with_a[mi] -> GetChi2A(), tud_b_with_a[mi] -> GetChi2L(), tud_b_with_a[mi] -> GetNdf() ) ;

         printf("\n\n %s : %11s with %11s : tau= %9.7f ,  chi2= (%9.5f + %9.5f)/%d\n",
            methods[mi], input_name_b, input_name_b, tud_b_with_b[mi] -> GetTau(), tud_b_with_b[mi] -> GetChi2A(), tud_b_with_b[mi] -> GetChi2L(), tud_b_with_b[mi] -> GetNdf() ) ;

         printf("\n\n") ;



         //==========================================================================
         // retrieve results into histograms

         // get unfolded distribution

         sprintf( hname, "h_unfolded_a_with_a_%s", methods[mi] ) ;
         h_unfolded_a_with_a[mi] = tud_a_with_a[mi] -> GetOutput( hname ) ;

         sprintf( hname, "h_unfolded_a_with_b_%s", methods[mi] ) ;
         h_unfolded_a_with_b[mi] = tud_a_with_b[mi] -> GetOutput( hname ) ;

         sprintf( hname, "h_unfolded_b_with_a_%s", methods[mi] ) ;
         h_unfolded_b_with_a[mi] = tud_b_with_a[mi] -> GetOutput( hname ) ;

         sprintf( hname, "h_unfolded_b_with_b_%s", methods[mi] ) ;
         h_unfolded_b_with_b[mi] = tud_b_with_b[mi] -> GetOutput( hname ) ;



         sprintf( hname, "h_gen_compare_a_with_a_%s", methods[mi] ) ;
         h_gen_compare_a_with_a[mi] = (TH1*) h_gen_source_a[mi] -> Clone( hname ) ;
         h_gen_compare_a_with_a[mi] -> Scale( ( h_unfolded_a_with_a[mi] -> Integral() )/( h_gen_source_a[mi] -> Integral() ) ) ;

         sprintf( hname, "h_gen_compare_a_with_b_%s", methods[mi] ) ;
         h_gen_compare_a_with_b[mi] = (TH1*) h_gen_source_a[mi] -> Clone( hname ) ;
         h_gen_compare_a_with_b[mi] -> Scale( ( h_unfolded_a_with_b[mi] -> Integral() )/( h_gen_source_a[mi] -> Integral() ) ) ;

         sprintf( hname, "h_gen_compare_b_with_a_%s", methods[mi] ) ;
         h_gen_compare_b_with_a[mi] = (TH1*) h_gen_source_b[mi] -> Clone( hname ) ;
         h_gen_compare_b_with_a[mi] -> Scale( ( h_unfolded_b_with_a[mi] -> Integral() )/( h_gen_source_b[mi] -> Integral() ) ) ;

         sprintf( hname, "h_gen_compare_b_with_b_%s", methods[mi] ) ;
         h_gen_compare_b_with_b[mi] = (TH1*) h_gen_source_b[mi] -> Clone( hname ) ;
         h_gen_compare_b_with_b[mi] -> Scale( ( h_unfolded_b_with_b[mi] -> Integral() )/( h_gen_source_b[mi] -> Integral() ) ) ;

         h_gen_compare_a_with_a[mi] -> SetLineColor( method_lc[mi] ) ;
         h_gen_compare_a_with_b[mi] -> SetLineColor( method_lc[mi] ) ;
         h_gen_compare_b_with_a[mi] -> SetLineColor( method_lc[mi] ) ;
         h_gen_compare_b_with_b[mi] -> SetLineColor( method_lc[mi] ) ;

         h_gen_compare_a_with_a[mi] -> SetLineStyle( method_ls[mi] ) ;
         h_gen_compare_a_with_b[mi] -> SetLineStyle( method_ls[mi] ) ;
         h_gen_compare_b_with_a[mi] -> SetLineStyle( method_ls[mi] ) ;
         h_gen_compare_b_with_b[mi] -> SetLineStyle( method_ls[mi] ) ;


         sprintf( hname, "h_diff_a_with_a_%s", methods[mi] ) ;
         h_diff_a_with_a[mi] = (TH1*) h_unfolded_a_with_a[mi] -> Clone( hname ) ;

         sprintf( hname, "h_diff_a_with_b_%s", methods[mi] ) ;
         h_diff_a_with_b[mi] = (TH1*) h_unfolded_a_with_b[mi] -> Clone( hname ) ;

         sprintf( hname, "h_diff_b_with_a_%s", methods[mi] ) ;
         h_diff_b_with_a[mi] = (TH1*) h_unfolded_b_with_a[mi] -> Clone( hname ) ;

         sprintf( hname, "h_diff_b_with_b_%s", methods[mi] ) ;
         h_diff_b_with_b[mi] = (TH1*) h_unfolded_b_with_b[mi] -> Clone( hname ) ;

         h_diff_a_with_a[mi] -> Add( h_gen_compare_a_with_a[mi], -1. ) ;
         h_diff_a_with_b[mi] -> Add( h_gen_compare_a_with_b[mi], -1. ) ;
         h_diff_b_with_a[mi] -> Add( h_gen_compare_b_with_a[mi], -1. ) ;
         h_diff_b_with_b[mi] -> Add( h_gen_compare_b_with_b[mi], -1. ) ;

       //--- 2022-03-14 : Change to relative (normalized) difference instead of absolute difference.
         h_diff_a_with_a[mi] -> Divide( h_gen_compare_a_with_a[mi] ) ;
         h_diff_a_with_b[mi] -> Divide( h_gen_compare_a_with_a[mi] ) ;
         h_diff_b_with_a[mi] -> Divide( h_gen_compare_a_with_b[mi] ) ;
         h_diff_b_with_b[mi] -> Divide( h_gen_compare_a_with_b[mi] ) ;

         printf("  h_diff_a_with_a        : %d in %9.5f to %9.5f\n", h_diff_a_with_a[mi]->GetNbinsX(), h_diff_a_with_a[mi]->GetBinLowEdge(1), h_diff_a_with_a[mi]->GetBinLowEdge(h_diff_a_with_a[mi]->GetNbinsX()+1) ) ;
         printf("  h_gen_compare_a_with_a : %d in %9.5f to %9.5f\n", h_gen_compare_a_with_a[mi]->GetNbinsX(), h_gen_compare_a_with_a[mi]->GetBinLowEdge(1), h_gen_compare_a_with_a[mi]->GetBinLowEdge(h_gen_compare_a_with_a[mi]->GetNbinsX()+1) ) ;

      } // mi


      for ( int mi=0; mi<3; mi++ ) {

         sprintf( hname, "h_abs_diff_a_with_a_%s", methods[mi] ) ;
         h_abs_diff_a_with_a[mi] = (TH1*) h_diff_a_with_a[mi] -> Clone( hname ) ;

         sprintf( hname, "h_abs_diff_a_with_b_%s", methods[mi] ) ;
         h_abs_diff_a_with_b[mi] = (TH1*) h_diff_a_with_b[mi] -> Clone( hname ) ;

         sprintf( hname, "h_abs_diff_b_with_a_%s", methods[mi] ) ;
         h_abs_diff_b_with_a[mi] = (TH1*) h_diff_b_with_a[mi] -> Clone( hname ) ;

         sprintf( hname, "h_abs_diff_b_with_b_%s", methods[mi] ) ;
         h_abs_diff_b_with_b[mi] = (TH1*) h_diff_b_with_b[mi] -> Clone( hname ) ;

         for ( int bi=1; bi <= h_diff_a_with_a[mi] -> GetNbinsX() ; bi++ ) {

            h_abs_diff_a_with_a[mi] -> SetBinContent( bi, fabs( h_diff_a_with_a[mi] -> GetBinContent(bi) ) ) ;
            h_abs_diff_a_with_b[mi] -> SetBinContent( bi, fabs( h_diff_a_with_b[mi] -> GetBinContent(bi) ) ) ;
            h_abs_diff_b_with_a[mi] -> SetBinContent( bi, fabs( h_diff_b_with_a[mi] -> GetBinContent(bi) ) ) ;
            h_abs_diff_b_with_b[mi] -> SetBinContent( bi, fabs( h_diff_b_with_b[mi] -> GetBinContent(bi) ) ) ;

         } // bi

         sprintf( hname, "h_ave_abs_diff_not_same_%s", methods[mi] ) ;
         h_ave_abs_diff_not_same[mi] = (TH1*) h_abs_diff_a_with_b[mi] -> Clone( hname ) ;
         h_ave_abs_diff_not_same[mi] -> Add( h_abs_diff_b_with_a[mi] ) ;
         h_ave_abs_diff_not_same[mi] -> Scale( 0.5 ) ;

      } // mi












    //--- Below here is all plotting stuff.


      float max_abs_diff = 0. ;
      for ( int mi=0; mi<3; mi++ ) {
         if ( h_diff_a_with_a[mi] -> GetMaximum() > max_abs_diff ) max_abs_diff = h_diff_a_with_a[mi] -> GetMaximum() ;
         if ( h_diff_a_with_b[mi] -> GetMaximum() > max_abs_diff ) max_abs_diff = h_diff_a_with_b[mi] -> GetMaximum() ;
         if ( h_diff_b_with_a[mi] -> GetMaximum() > max_abs_diff ) max_abs_diff = h_diff_b_with_a[mi] -> GetMaximum() ;
         if ( h_diff_b_with_b[mi] -> GetMaximum() > max_abs_diff ) max_abs_diff = h_diff_b_with_b[mi] -> GetMaximum() ;
         if ( fabs( h_diff_a_with_a[mi] -> GetMinimum() ) > max_abs_diff ) max_abs_diff = fabs( h_diff_a_with_a[mi] -> GetMinimum() ) ;
         if ( fabs( h_diff_a_with_b[mi] -> GetMinimum() ) > max_abs_diff ) max_abs_diff = fabs( h_diff_a_with_b[mi] -> GetMinimum() ) ;
         if ( fabs( h_diff_b_with_a[mi] -> GetMinimum() ) > max_abs_diff ) max_abs_diff = fabs( h_diff_b_with_a[mi] -> GetMinimum() ) ;
         if ( fabs( h_diff_b_with_b[mi] -> GetMinimum() ) > max_abs_diff ) max_abs_diff = fabs( h_diff_b_with_b[mi] -> GetMinimum() ) ;
      } // mi

      if ( diff_max > 0 ) max_abs_diff = diff_max ;

     //--- 2022-03-14 : Set by hand
      if ( strcmp( var_name, "x" ) == 0 ) { max_abs_diff = 0.11 ; }
      if ( strcmp( var_name, "y" ) == 0 ) { max_abs_diff = 0.11 ; }

      for ( int mi=0; mi<3; mi++ ) {

         h_diff_a_with_a[mi] -> SetMaximum( 1.1 * max_abs_diff ) ;
         h_diff_a_with_b[mi] -> SetMaximum( 1.1 * max_abs_diff ) ;
         h_diff_b_with_a[mi] -> SetMaximum( 1.1 * max_abs_diff ) ;
         h_diff_b_with_b[mi] -> SetMaximum( 1.1 * max_abs_diff ) ;

         h_diff_a_with_a[mi] -> SetMinimum( -1.1 * max_abs_diff ) ;
         h_diff_a_with_b[mi] -> SetMinimum( -1.1 * max_abs_diff ) ;
         h_diff_b_with_a[mi] -> SetMinimum( -1.1 * max_abs_diff ) ;
         h_diff_b_with_b[mi] -> SetMinimum( -1.1 * max_abs_diff ) ;

         h_diff_a_with_a[mi] -> SetLineColor( method_lc[mi] ) ;
         h_diff_a_with_b[mi] -> SetLineColor( method_lc[mi] ) ;
         h_diff_b_with_a[mi] -> SetLineColor( method_lc[mi] ) ;
         h_diff_b_with_b[mi] -> SetLineColor( method_lc[mi] ) ;

         h_diff_a_with_a[mi] -> SetLineStyle( method_ls[mi] ) ;
         h_diff_a_with_b[mi] -> SetLineStyle( method_ls[mi] ) ;
         h_diff_b_with_a[mi] -> SetLineStyle( method_ls[mi] ) ;
         h_diff_b_with_b[mi] -> SetLineStyle( method_ls[mi] ) ;

         h_diff_a_with_a[mi] -> SetLineWidth( 3 ) ;
         h_diff_a_with_b[mi] -> SetLineWidth( 3 ) ;
         h_diff_b_with_a[mi] -> SetLineWidth( 3 ) ;
         h_diff_b_with_b[mi] -> SetLineWidth( 3 ) ;

         h_abs_diff_a_with_a[mi] -> SetMaximum( 1.1 * max_abs_diff ) ;
         h_abs_diff_a_with_b[mi] -> SetMaximum( 1.1 * max_abs_diff ) ;
         h_abs_diff_b_with_a[mi] -> SetMaximum( 1.1 * max_abs_diff ) ;
         h_abs_diff_b_with_b[mi] -> SetMaximum( 1.1 * max_abs_diff ) ;

         h_abs_diff_a_with_a[mi] -> SetMinimum(-0.1 * max_abs_diff ) ;
         h_abs_diff_a_with_b[mi] -> SetMinimum(-0.1 * max_abs_diff ) ;
         h_abs_diff_b_with_a[mi] -> SetMinimum(-0.1 * max_abs_diff ) ;
         h_abs_diff_b_with_b[mi] -> SetMinimum(-0.1 * max_abs_diff ) ;

         h_abs_diff_a_with_a[mi] -> SetLineColor( method_lc[mi] ) ;
         h_abs_diff_a_with_b[mi] -> SetLineColor( method_lc[mi] ) ;
         h_abs_diff_b_with_a[mi] -> SetLineColor( method_lc[mi] ) ;
         h_abs_diff_b_with_b[mi] -> SetLineColor( method_lc[mi] ) ;

         h_abs_diff_a_with_a[mi] -> SetLineStyle( method_ls[mi] ) ;
         h_abs_diff_a_with_b[mi] -> SetLineStyle( method_ls[mi] ) ;
         h_abs_diff_b_with_a[mi] -> SetLineStyle( method_ls[mi] ) ;
         h_abs_diff_b_with_b[mi] -> SetLineStyle( method_ls[mi] ) ;

         h_abs_diff_a_with_a[mi] -> SetLineWidth( 3 ) ;
         h_abs_diff_a_with_b[mi] -> SetLineWidth( 3 ) ;
         h_abs_diff_b_with_a[mi] -> SetLineWidth( 3 ) ;
         h_abs_diff_b_with_b[mi] -> SetLineWidth( 3 ) ;

         h_ave_abs_diff_not_same[mi] -> SetMaximum(  1.1 * max_abs_diff ) ;
         h_ave_abs_diff_not_same[mi] -> SetMinimum( -0.1 * max_abs_diff ) ;

         h_ave_abs_diff_not_same[mi] -> SetLineColor( method_lc[mi] ) ;

         h_ave_abs_diff_not_same[mi] -> SetLineStyle( method_ls[mi] ) ;

         h_ave_abs_diff_not_same[mi] -> SetLineWidth( 3 ) ;

         if ( strcmp( var_name, "x" ) == 0 ) { h_diff_a_with_a[mi] -> SetXTitle( "log10(x)" ) ; }
         if ( strcmp( var_name, "x" ) == 0 ) { h_diff_a_with_b[mi] -> SetXTitle( "log10(x)" ) ; }
         if ( strcmp( var_name, "x" ) == 0 ) { h_diff_b_with_a[mi] -> SetXTitle( "log10(x)" ) ; }
         if ( strcmp( var_name, "x" ) == 0 ) { h_diff_b_with_b[mi] -> SetXTitle( "log10(x)" ) ; }


         if ( strcmp( var_name, "y" ) == 0 ) { h_diff_a_with_a[mi] -> SetXTitle( "log10(y)" ) ; }
         if ( strcmp( var_name, "y" ) == 0 ) { h_diff_a_with_b[mi] -> SetXTitle( "log10(y)" ) ; }
         if ( strcmp( var_name, "y" ) == 0 ) { h_diff_b_with_a[mi] -> SetXTitle( "log10(y)" ) ; }
         if ( strcmp( var_name, "y" ) == 0 ) { h_diff_b_with_b[mi] -> SetXTitle( "log10(y)" ) ; }

         h_diff_a_with_a[mi] -> SetTitleOffset( 1.4, "x" ) ;
         h_diff_a_with_b[mi] -> SetTitleOffset( 1.4, "x" ) ;
         h_diff_b_with_a[mi] -> SetTitleOffset( 1.4, "x" ) ;
         h_diff_b_with_b[mi] -> SetTitleOffset( 1.4, "x" ) ;

         ///////h_diff_a_with_a[mi] -> SetYTitle( "Unfolded - Gen (events)" ) ;
         ///////h_diff_a_with_b[mi] -> SetYTitle( "Unfolded - Gen (events)" ) ;
         ///////h_diff_b_with_a[mi] -> SetYTitle( "Unfolded - Gen (events)" ) ;
         ///////h_diff_b_with_b[mi] -> SetYTitle( "Unfolded - Gen (events)" ) ;

         h_diff_a_with_a[mi] -> SetYTitle( "(Unfolded - Gen)/Gen" ) ;
         h_diff_a_with_b[mi] -> SetYTitle( "(Unfolded - Gen)/Gen" ) ;
         h_diff_b_with_a[mi] -> SetYTitle( "(Unfolded - Gen)/Gen" ) ;
         h_diff_b_with_b[mi] -> SetYTitle( "(Unfolded - Gen)/Gen" ) ;

         h_diff_a_with_a[mi] -> SetTitleOffset( 1.9, "y" ) ;
         h_diff_a_with_b[mi] -> SetTitleOffset( 1.9, "y" ) ;
         h_diff_b_with_a[mi] -> SetTitleOffset( 1.9, "y" ) ;
         h_diff_b_with_b[mi] -> SetTitleOffset( 1.9, "y" ) ;

         h_diff_a_with_a[mi] -> SetTitle( "Closure test" ) ;

         sprintf( htitle, "Generator model systematic, %s", methods_title[mi] ) ;
         h_diff_a_with_b[mi] -> SetTitle( htitle ) ;
         h_diff_b_with_a[mi] -> SetTitle( htitle ) ;
         h_diff_b_with_b[mi] -> SetTitle( htitle ) ;

         h_diff_a_with_a[mi] -> SetTitleSize( 0.045, "x" ) ;
         h_diff_a_with_b[mi] -> SetTitleSize( 0.045, "x" ) ;
         h_diff_b_with_a[mi] -> SetTitleSize( 0.045, "x" ) ;
         h_diff_b_with_b[mi] -> SetTitleSize( 0.045, "x" ) ;

         h_diff_a_with_a[mi] -> SetTitleSize( 0.045, "y" ) ;
         h_diff_a_with_b[mi] -> SetTitleSize( 0.045, "y" ) ;
         h_diff_b_with_a[mi] -> SetTitleSize( 0.045, "y" ) ;
         h_diff_b_with_b[mi] -> SetTitleSize( 0.045, "y" ) ;

         h_diff_a_with_a[mi] -> SetNdivisions( 605 ) ;
         h_diff_a_with_b[mi] -> SetNdivisions( 605 ) ;
         h_diff_b_with_a[mi] -> SetNdivisions( 605 ) ;
         h_diff_b_with_b[mi] -> SetNdivisions( 605 ) ;
      }








      gStyle -> SetPadRightMargin(0.08) ;
      gStyle -> SetPadLeftMargin(0.17) ;
      gStyle -> SetPadBottomMargin(0.15) ;
      gStyle -> SetTitleBorderSize(0) ;
      gStyle -> SetTitleY(0.975) ;
      gStyle -> SetOptStat(0) ;

      gStyle -> SetTitleX(0.2) ;




      int can_width = 600 ;
      int can_height = 600 ;

      int can_spacing = 20 ;

      float lx, ly, lw, lh ;
      char ltitle[100] ;


     //-----

      TCanvas* can1 = (TCanvas*) gDirectory -> FindObject( "can1" ) ;
      if ( can1 == 0x0 ) { printf( "Making can1\n") ; can1 = new TCanvas( "can1", "", can_spacing, can_spacing, can_width, can_height ) ; }
      can1 -> Clear() ;
      can1 -> cd() ;

      h_diff_a_with_a[0] -> Draw("") ;
      h_diff_a_with_a[1] -> Draw("same") ;
      h_diff_a_with_a[2] -> Draw("same") ;

      gPad -> SetGridy(1) ;

      lx = 0.65 ; ly = 0.73 ; lw = 0.24; lh = 0.15 ;

      TLegend* legend1 = new TLegend( lx, ly, lx+lw, ly+lh ) ;
      legend1 -> AddEntry( h_diff_a_with_a[0], methods_title[0] ) ;
      legend1 -> AddEntry( h_diff_a_with_a[1], methods_title[1] ) ;
      legend1 -> AddEntry( h_diff_a_with_a[2], methods_title[2] ) ;

      legend1 -> Draw() ;

      sprintf( fname, "paper-plots/dis-closure-%s.png", var_name ) ;
      can1 -> SaveAs( fname ) ;
      sprintf( fname, "paper-plots/dis-closure-%s.pdf", var_name ) ;
      can1 -> SaveAs( fname ) ;

      can1 -> Update() ; can1 -> Draw() ; gSystem -> ProcessEvents() ;



     //-----

      TCanvas* can2 = (TCanvas*) gDirectory -> FindObject( "can2" ) ;
      if ( can2 == 0x0 ) { printf( "Making can2\n") ; can2 = new TCanvas( "can2", "", 2*can_spacing + can_width , can_spacing, can_width, can_height ) ; }
      can2 -> Clear() ;
      can2 -> cd() ;

      h_diff_a_with_b[0] -> SetLineStyle(1) ;
      h_diff_a_with_b[0] -> SetMarkerStyle(20) ;
      h_diff_a_with_b[0] -> SetMarkerSize(1.5) ;
      h_diff_a_with_b[0] -> SetMarkerColor( method_lc[0] ) ;
      h_diff_b_with_a[0] -> SetLineStyle(9) ;
      h_diff_b_with_a[0] -> SetMarkerStyle(26) ;
      h_diff_b_with_a[0] -> SetMarkerSize(1.5) ;
      h_diff_b_with_a[0] -> SetMarkerColor( method_lc[0] ) ;

      h_diff_a_with_b[0] -> Draw("") ;
      h_diff_b_with_a[0] -> Draw("same") ;
      h_diff_a_with_b[0] -> Draw("hist same") ;
      h_diff_b_with_a[0] -> Draw("hist same") ;

      gPad -> SetGridy(1) ;

      if ( strcmp( var_name, "x" ) == 0 ) {
         lx = 0.28 ; ly = 0.73 ; lw = 0.44; lh = 0.15 ;
      } else {
         lx = 0.52 ; ly = 0.73 ; lw = 0.44; lh = 0.15 ;
      }

      TLegend* legend2 = new TLegend( lx, ly, lx+lw, ly+lh ) ;

      sprintf( ltitle, "%s with %s", input_name_a, input_name_b ) ;
      legend2 -> AddEntry( h_diff_a_with_b[0], ltitle  ) ;
      sprintf( ltitle, "%s with %s", input_name_b, input_name_a ) ;
      legend2 -> AddEntry( h_diff_b_with_a[0], ltitle ) ;

      legend2 -> Draw() ;

      sprintf( fname, "paper-plots/dis-gen-syst-%s-%s.png", var_name, methods[0] ) ;
      can2 -> SaveAs( fname ) ;
      sprintf( fname, "paper-plots/dis-gen-syst-%s-%s.pdf", var_name, methods[0] ) ;
      can2 -> SaveAs( fname ) ;

      can2 -> Update() ; can2 -> Draw() ; gSystem -> ProcessEvents() ;





     //-----

      TCanvas* can3 = (TCanvas*) gDirectory -> FindObject( "can3" ) ;
      if ( can3 == 0x0 ) { printf( "Making can3\n") ; can3 = new TCanvas( "can3", "", 3*can_spacing + 2*can_width , can_spacing, can_width, can_height ) ; }
      can3 -> Clear() ;
      can3 -> cd() ;

      h_diff_a_with_b[1] -> SetLineStyle(1) ;
      h_diff_a_with_b[1] -> SetMarkerStyle(20) ;
      h_diff_a_with_b[1] -> SetMarkerSize(1.5) ;
      h_diff_a_with_b[1] -> SetMarkerColor( method_lc[1] ) ;
      h_diff_b_with_a[1] -> SetLineStyle(9) ;
      h_diff_b_with_a[1] -> SetMarkerStyle(26) ;
      h_diff_b_with_a[1] -> SetMarkerSize(1.5) ;
      h_diff_b_with_a[1] -> SetMarkerColor( method_lc[1] ) ;

      h_diff_a_with_b[1] -> Draw("") ;
      h_diff_b_with_a[1] -> Draw("same") ;
      h_diff_a_with_b[1] -> Draw("hist same") ;
      h_diff_b_with_a[1] -> Draw("hist same") ;


      gPad -> SetGridy(1) ;


      TLegend* legend3 = new TLegend( lx, ly, lx+lw, ly+lh ) ;

      sprintf( ltitle, "%s with %s", input_name_a, input_name_b ) ;
      legend3 -> AddEntry( h_diff_a_with_b[1], ltitle  ) ;
      sprintf( ltitle, "%s with %s", input_name_b, input_name_a ) ;
      legend3 -> AddEntry( h_diff_b_with_a[1], ltitle ) ;

      legend3 -> Draw() ;

      sprintf( fname, "paper-plots/dis-gen-syst-%s-%s.png", var_name, methods[1] ) ;
      can3 -> SaveAs( fname ) ;
      sprintf( fname, "paper-plots/dis-gen-syst-%s-%s.pdf", var_name, methods[1] ) ;
      can3 -> SaveAs( fname ) ;

      can3 -> Update() ; can2 -> Draw() ; gSystem -> ProcessEvents() ;





     //-----

      TCanvas* can4 = (TCanvas*) gDirectory -> FindObject( "can4" ) ;
      if ( can4 == 0x0 ) { printf( "Making can4\n") ; can4 = new TCanvas( "can4", "", 4*can_spacing + 3*can_width , can_spacing, can_width, can_height ) ; }
      can4 -> Clear() ;
      can4 -> cd() ;

      h_diff_a_with_b[2] -> SetLineStyle(1) ;
      h_diff_a_with_b[2] -> SetMarkerStyle(20) ;
      h_diff_a_with_b[2] -> SetMarkerSize(1.5) ;
      h_diff_a_with_b[2] -> SetMarkerColor( method_lc[2] ) ;
      h_diff_b_with_a[2] -> SetLineStyle(9) ;
      h_diff_b_with_a[2] -> SetMarkerStyle(26) ;
      h_diff_b_with_a[2] -> SetMarkerSize(1.5) ;
      h_diff_b_with_a[2] -> SetMarkerColor( method_lc[2] ) ;

      h_diff_a_with_b[2] -> Draw("") ;
      h_diff_b_with_a[2] -> Draw("same") ;
      h_diff_a_with_b[2] -> Draw("hist same") ;
      h_diff_b_with_a[2] -> Draw("hist same") ;


      gPad -> SetGridy(1) ;


      TLegend* legend4 = new TLegend( lx, ly, lx+lw, ly+lh ) ;

      sprintf( ltitle, "%s with %s", input_name_a, input_name_b ) ;
      legend4 -> AddEntry( h_diff_a_with_b[2], ltitle  ) ;
      sprintf( ltitle, "%s with %s", input_name_b, input_name_a ) ;
      legend4 -> AddEntry( h_diff_b_with_a[2], ltitle ) ;

      legend4 -> Draw() ;


      sprintf( fname, "paper-plots/dis-gen-syst-%s-%s.png", var_name, methods[2] ) ;
      can4 -> SaveAs( fname ) ;
      sprintf( fname, "paper-plots/dis-gen-syst-%s-%s.pdf", var_name, methods[2] ) ;
      can4 -> SaveAs( fname ) ;

      can4 -> Update() ; can2 -> Draw() ; gSystem -> ProcessEvents() ;









   }













