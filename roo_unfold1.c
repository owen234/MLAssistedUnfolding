
#include "histio.c"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldTUnfold.h"

TH2F* get_hist( const char* hname ) {
   TH2F* hp = (TH2F*) gDirectory -> FindObject( hname ) ;
   if ( hp == 0x0 ) { printf("\n\n *** can't find %s\n\n", hname ) ; gSystem->Exit(-1) ; }
   return hp ;
}

//
//
//  Before doing
//
//    .L roo_unfold.c
//
//  do
//
//   gSystem -> Load( "../RooUnfold/build/libRooUnfold.dylib" ) ;
//
//

   void roo_unfold1( const char* hist_name = "h_log10_x_gen_vs_obs_dnn",
                     int ngen = 1e6,
                     const char* input_file = "h1-input-nbins_gen020_obs100-b2c.root" ) {

      gStyle -> SetOptStat(0) ;

      TRandom3 tran(12345) ;

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


     //-- Add underflows and overflows to fake data!
      float underflow_mean ;
      float overflow_mean ;
      int nbins_obs = h_obs_source -> GetNbinsX() ;
      float nobs_source ;
      float nobs_integral ;
      float nobs_underflow ;
      float nobs_overflow ;



      nobs_integral = h_obs_source -> Integral() ;
      nobs_underflow = h_obs_source -> GetBinContent(0) ;
      nobs_overflow = h_obs_source -> GetBinContent( nbins_obs+1 ) ;

      printf("  h_obs_source :  integral %9.5f, underflow %9.5f, overflow %9.5f\n", nobs_integral, nobs_underflow, nobs_overflow ) ;

      nobs_source = nobs_integral ;
      nobs_source += nobs_underflow ;
      nobs_source += nobs_overflow ;

      ////////underflow_mean = ngen * ( nobs_underflow / nobs_source ) ;
      ////////overflow_mean = ngen * ( nobs_overflow / nobs_source ) ;

      underflow_mean = ngen * ( nobs_underflow / nobs_integral ) ; // FillRandom only fills inside of histogram (integral)
      overflow_mean = ngen * ( nobs_overflow / nobs_integral ) ;

      h_obs_random -> SetBinContent( 0, tran.Poisson( underflow_mean ) ) ;
      h_obs_random -> SetBinContent( nbins_obs+1, tran.Poisson( overflow_mean ) ) ;

      h_obs_random -> SetBinError( 0, sqrt( underflow_mean ) ) ;
      h_obs_random -> SetBinError( nbins_obs+1, sqrt( overflow_mean ) ) ;




      //////////RooUnfoldResponse* rur = new RooUnfoldResponse( h_obs_source, h_gen_source, h_in_gen_vs_obs ) ;
      RooUnfoldResponse* rur = new RooUnfoldResponse( 0, 0, h_in_gen_vs_obs ) ;

      rur -> UseOverflow(true) ;

      RooUnfoldBayes  unfold( rur, h_obs_random, 4 ) ;
      //RooUnfoldTUnfold  unfold( rur, h_obs_random ) ;


      TH1D* hReco = (TH1D*) unfold.Hreco() ;
      hReco -> SetName( "h_unfolded" ) ;
      hReco -> Draw() ;

      TH1* h_gen_compare = (TH1*) h_gen_source -> Clone( "h_gen_compare" ) ;
      h_gen_compare -> Scale( ( hReco -> Integral() )/( h_gen_source -> Integral() ) ) ;

      h_gen_compare -> SetLineColor(2) ;
      h_gen_compare -> Draw("hist same") ;


   }












