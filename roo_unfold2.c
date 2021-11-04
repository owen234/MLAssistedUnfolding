
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
//    .L roo_unfold2.c
//
//  do
//
//   gSystem -> Load( "../RooUnfold/build/libRooUnfold.dylib" ) ;
//
//

   void roo_unfold2( const char* rur_name = "rur_log10_x_gen_vs_obs_e",
                     int ngen = 1e6,
                     int method_index = 1,
                     int n_iter = 1000,  // this is huge because I hacked RooUnfoldBayes to start with a flat prior.
                     const char* input_file = "example-input-nbins_gen010_obs020.root" ) {

      gStyle -> SetPadRightMargin(0.15) ;

      TCanvas* can_ru = (TCanvas*) gDirectory -> FindObject( "can_ru" ) ;
      if ( can_ru == 0x0 ) { can_ru = new TCanvas( "can_ru", "", 50, 50, 1500, 1500 ) ; }
      can_ru -> Clear() ;
      can_ru -> Divide(3,2) ;

      gStyle -> SetOptStat(0) ;

      TRandom3 tran(12345) ;

      gDirectory -> Delete( "*" ) ;

      TFile f( input_file, "read" ) ;
      f.ls() ;

      RooUnfoldResponse* rur = (RooUnfoldResponse*) f.Get( rur_name ) ;
      if ( rur == 0x0 ) { printf("\n\n *** can't find %s\n\n", rur_name ) ; return ; }


      TH2* h_in_gen_vs_obs = rur -> Hresponse() ;

      TH1* h_obs_source = rur -> Hmeasured() ;
      TH1* h_gen_source = rur -> Htruth() ;

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



      rur -> UseOverflow(true) ;

      RooUnfold* unfold ;
      if ( method_index == 1 ) {
         unfold = new RooUnfoldBayes( rur, h_obs_random, n_iter ) ;
      } else if ( method_index == 2 ) {
         unfold = new RooUnfoldTUnfold( rur, h_obs_random ) ;
      } else {
         printf("\n\n *** I don't know method_index = %d\n\n", method_index ) ;
         return ;
      }


      //TH1D* hReco = (TH1D*) unfold -> Hreco() ;
      TH1D* hReco = (TH1D*) unfold -> Hreco(RooUnfold::kCovariance) ;
      //TH1D* hReco = (TH1D*) unfold -> Hreco(RooUnfold::kCovToy) ;


      unfold -> PrintTable( cout ) ;

      TMatrixD unfolding_cov_mat = unfold -> Ereco() ;

      TMatrixD unfolding_inverse_cov_mat = unfold -> Wreco() ;


      unfolding_cov_mat.Print() ;

      TH2F* h_unfold_cov_mat = new TH2F( "h_unfold_cov_mat", "Unfolding cov. mat.",
          unfolding_cov_mat.GetNcols(), -0.5, unfolding_cov_mat.GetNcols()-0.5,
          unfolding_cov_mat.GetNcols(), -0.5, unfolding_cov_mat.GetNcols()-0.5 ) ;

      TH2F* h_unfold_cor_mat = new TH2F( "h_unfold_cor_mat", "Unfolding cov. mat.",
          unfolding_cov_mat.GetNcols(), -0.5, unfolding_cov_mat.GetNcols()-0.5,
          unfolding_cov_mat.GetNcols(), -0.5, unfolding_cov_mat.GetNcols()-0.5 ) ;

      TH1F* h_global_correlation_coeff = new TH1F( "h_global_correlation_coeff", "Global cor. coef.",
         unfolding_cov_mat.GetNcols(), -0.5, unfolding_cov_mat.GetNcols()-0.5 ) ;

      for ( int ri=0; ri<unfolding_cov_mat.GetNcols(); ri++ ) {
         int rbi = ri+1 ;
         for ( int ci=0; ci<unfolding_cov_mat.GetNcols(); ci++ ) {
            int cbi = ci+1 ;
            h_unfold_cov_mat -> SetBinContent( rbi, cbi, unfolding_cov_mat[ri][ci] ) ;
            float rho = 1. ;
            float err_i = sqrt( unfolding_cov_mat[ri][ri] ) ;
            float err_j = sqrt( unfolding_cov_mat[ci][ci] ) ;
            if ( err_i > 0 && err_j > 0 ) rho = unfolding_cov_mat[ri][ci] / ( err_i * err_j ) ;
            h_unfold_cor_mat -> SetBinContent( rbi, cbi, rho ) ;
         } // ci
         float global_rho = 1. ;
         if ( unfolding_cov_mat[ri][ri] != 0 && unfolding_inverse_cov_mat[ri][ri] != 0 ) {
            global_rho = sqrt( 1. - 1. / ( unfolding_cov_mat[ri][ri] * unfolding_inverse_cov_mat[ri][ri] ) ) ;
         }
         h_global_correlation_coeff -> SetBinContent( rbi, global_rho ) ;
      } // ri

      /////h_unfold_cov_mat -> Print("all") ;



      const TMatrixD& response_matrix = rur -> Mresponse() ;
      response_matrix.Print() ;

      TH2F* h_response_matrix_normalized = new TH2F( "h_response_matrix_normalized", "Response matrix, normalized",
          response_matrix.GetNrows(), -0.5, response_matrix.GetNrows()-0.5,
          response_matrix.GetNcols(), -0.5, response_matrix.GetNcols()-0.5
          ) ;
      for ( int i=0; i<response_matrix.GetNrows(); i++ ) {
         for ( int j=0; j<response_matrix.GetNcols(); j++ ) {
            h_response_matrix_normalized -> SetBinContent( i+1, j+1,  response_matrix[i][j] ) ;
         } // j
      } // i



      TH1* h_gen_compare = (TH1*) h_gen_source -> Clone( "h_gen_compare" ) ;
      h_gen_compare -> Scale( ( hReco -> Integral() )/( h_gen_source -> Integral() ) ) ;




      int ci = 1;

      can_ru -> cd( ci++ ) ;

      h_in_gen_vs_obs -> DrawCopy("colz") ;

      can_ru -> cd( ci++ ) ;
      h_response_matrix_normalized -> DrawCopy( "colz") ;

      can_ru -> cd( ci++ ) ;

      hReco -> SetName( "h_unfolded" ) ;
      hReco -> DrawCopy() ;

      h_gen_compare -> SetLineColor(2) ;
      h_gen_compare -> DrawCopy("hist same") ;




      can_ru -> cd( ci++ ) ;

      h_unfold_cov_mat -> DrawCopy("colz") ;


      can_ru -> cd( ci++ ) ;

      h_unfold_cor_mat -> DrawCopy("colz") ;


      can_ru -> cd( ci++ ) ;
      h_global_correlation_coeff -> DrawCopy() ;




 ///  RooUnfoldErrors* rue = new RooUnfoldErrors( 100, unfold, h_gen_compare ) ;
 ///  TH1* unfolding_error = rue -> UnfoldingError() ;
 ///  TH1* unfolding_toy_rms = rue -> RMSResiduals() ;


 ///  can_ru -> cd( ci++ ) ;
 ///  unfolding_error -> DrawCopy() ;
 ///  unfolding_toy_rms -> DrawCopy("same") ;


      f.Close() ;


      gDirectory -> ls() ;


      printf("\n\n\n") ;

      printf("   Cut and paste for this:\n\n") ;

      printf("     roo_unfold2(\"%s\",%d,%d,%d,\"%s\")\n", rur_name, ngen, method_index, n_iter, input_file ) ;

      printf("\n\n\n") ;

   }












