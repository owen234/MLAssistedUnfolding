
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

   void roo_unfold_comp2( const char* rur_name = "rur_log10_x_gen_vs_obs_e",
                     int ngen = 1e5,
                     int n_iter = 10000, // this is huge because I hacked RooUnfoldBayes to have a flat prior instead of the truth distribution.
                     const char* input_file = "example-input-nbins_gen010_obs020.root" ) {

      gStyle -> SetPadRightMargin(0.15) ;

      TCanvas* can_ru = (TCanvas*) gDirectory -> FindObject( "can_ru" ) ;
      if ( can_ru == 0x0 ) { can_ru = new TCanvas( "can_ru", "", 50, 50, 1700, 800 ) ; }
      can_ru -> Clear() ;

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


      RooUnfold* unfold_ibu = new RooUnfoldBayes( rur, h_obs_random, n_iter ) ;
      RooUnfold* unfold_tu  = new RooUnfoldTUnfold( rur, h_obs_random ) ;



      TH1D* hReco_ibu = (TH1D*) unfold_ibu -> Hreco(RooUnfold::kCovariance) ;
      TMatrixD unfolding_cov_mat_ibu = unfold_ibu -> Ereco() ;
      TVectorD unfolding_err_ibu = unfold_ibu -> ErecoV( RooUnfold::kCovariance ) ;
      ////TVectorD unfolding_err_ibu = unfold_ibu -> ErecoV( RooUnfold::kErrors ) ;
      TMatrixD unfolding_inverse_cov_mat_ibu = unfold_ibu -> Wreco() ;

      TH1D* hReco_tu = (TH1D*) unfold_tu -> Hreco(RooUnfold::kCovariance) ;
      TMatrixD unfolding_cov_mat_tu = unfold_tu -> Ereco() ;
      TVectorD unfolding_err_tu = unfold_tu -> ErecoV( RooUnfold::kCovariance ) ;
      ////TVectorD unfolding_err_tu = unfold_tu -> ErecoV( RooUnfold::kErrors ) ;
      TMatrixD unfolding_inverse_cov_mat_tu = unfold_tu -> Wreco() ;


      /////unfold_ibu -> PrintTable( cout ) ;
      //////unfolding_cov_mat_ibu.Print() ;



      TH2F* h_unfold_cov_mat_ibu = new TH2F( "h_unfold_cov_mat_ibu", "IBU: Unfolding cov. mat.",
          unfolding_cov_mat_ibu.GetNcols(), -0.5, unfolding_cov_mat_ibu.GetNcols()-0.5,
          unfolding_cov_mat_ibu.GetNcols(), -0.5, unfolding_cov_mat_ibu.GetNcols()-0.5 ) ;

      TH2F* h_unfold_cor_mat_ibu = new TH2F( "h_unfold_cor_mat_ibu", "IBU: Unfolding cov. mat.",
          unfolding_cov_mat_ibu.GetNcols(), -0.5, unfolding_cov_mat_ibu.GetNcols()-0.5,
          unfolding_cov_mat_ibu.GetNcols(), -0.5, unfolding_cov_mat_ibu.GetNcols()-0.5 ) ;

      TH1F* h_global_correlation_coeff_ibu = new TH1F( "h_global_correlation_coeff_ibu", "IBU: Global cor. coef.",
         unfolding_cov_mat_ibu.GetNcols(), -0.5, unfolding_cov_mat_ibu.GetNcols()-0.5 ) ;




      TH2F* h_unfold_cov_mat_tu = new TH2F( "h_unfold_cov_mat_tu", "TUnfold: Unfolding cov. mat.",
          unfolding_cov_mat_tu.GetNcols(), -0.5, unfolding_cov_mat_tu.GetNcols()-0.5,
          unfolding_cov_mat_tu.GetNcols(), -0.5, unfolding_cov_mat_tu.GetNcols()-0.5 ) ;

      TH2F* h_unfold_cor_mat_tu = new TH2F( "h_unfold_cor_mat_tu", "TUnfold: Unfolding cov. mat.",
          unfolding_cov_mat_tu.GetNcols(), -0.5, unfolding_cov_mat_tu.GetNcols()-0.5,
          unfolding_cov_mat_tu.GetNcols(), -0.5, unfolding_cov_mat_tu.GetNcols()-0.5 ) ;

      TH1F* h_global_correlation_coeff_tu = new TH1F( "h_global_correlation_coeff_tu", "TUnfold: Global cor. coef.",
         unfolding_cov_mat_tu.GetNcols(), -0.5, unfolding_cov_mat_tu.GetNcols()-0.5 ) ;



      for ( int ri=0; ri<unfolding_cov_mat_ibu.GetNcols(); ri++ ) {
         int rbi = ri+1 ;
         for ( int ci=0; ci<unfolding_cov_mat_ibu.GetNcols(); ci++ ) {

            int cbi = ci+1 ;

            h_unfold_cov_mat_ibu -> SetBinContent( rbi, cbi, unfolding_cov_mat_ibu[ri][ci] ) ;
            h_unfold_cov_mat_tu  -> SetBinContent( rbi, cbi, unfolding_cov_mat_tu[ri][ci] ) ;

            float rho = 1. ;
            float err_i = 1 ;
            float err_j = 1 ;

            err_i = sqrt( unfolding_cov_mat_ibu[ri][ri] ) ;
            err_j = sqrt( unfolding_cov_mat_ibu[ci][ci] ) ;
            if ( err_i > 0 && err_j > 0 ) rho = unfolding_cov_mat_ibu[ri][ci] / ( err_i * err_j ) ;
            h_unfold_cor_mat_ibu -> SetBinContent( rbi, cbi, rho ) ;
            if ( ri == ci ) {
               printf("   bin %3d     IBU: sqrt(cov_ii) = %9.2f, hist err = %9.2f ,   error = %9.2f  \n", ri, err_i, hReco_ibu->GetBinError( ri ), unfolding_err_ibu[ri] ) ;
            }

            err_i = sqrt( unfolding_cov_mat_tu [ri][ri] ) ;
            err_j = sqrt( unfolding_cov_mat_tu [ci][ci] ) ;
            if ( err_i > 0 && err_j > 0 ) rho = unfolding_cov_mat_tu[ri][ci] / ( err_i * err_j ) ;
            h_unfold_cor_mat_tu  -> SetBinContent( rbi, cbi, rho ) ;
            if ( ri == ci ) {
               printf("   bin %3d     T U: sqrt(cov_ii) = %9.2f, hist err = %9.2f ,   error = %9.2f  \n", ri, err_i, hReco_tu->GetBinError( ri ), unfolding_err_tu[ri] ) ;
            }




         } // ci

         float global_rho = 1. ;

         if ( unfolding_cov_mat_ibu[ri][ri] != 0 && unfolding_inverse_cov_mat_ibu[ri][ri] != 0 ) {
            global_rho = sqrt( 1. - 1. / ( unfolding_cov_mat_ibu[ri][ri] * unfolding_inverse_cov_mat_ibu[ri][ri] ) ) ;
         }
         h_global_correlation_coeff_ibu -> SetBinContent( rbi, global_rho ) ;

         if ( unfolding_cov_mat_tu [ri][ri] != 0 && unfolding_inverse_cov_mat_tu [ri][ri] != 0 ) {
            global_rho = sqrt( 1. - 1. / ( unfolding_cov_mat_tu [ri][ri] * unfolding_inverse_cov_mat_tu [ri][ri] ) ) ;
         }
         h_global_correlation_coeff_tu  -> SetBinContent( rbi, global_rho ) ;

      } // ri

      /////h_unfold_cov_mat_ibu -> Print("all") ;



////  const TMatrixD& response_matrix = rur -> Mresponse() ;
////  response_matrix.Print() ;

////  TH2F* h_response_matrix_normalized = new TH2F( "h_response_matrix_normalized", "Response matrix, normalized",
////      response_matrix.GetNrows(), -0.5, response_matrix.GetNrows()-0.5,
////      response_matrix.GetNcols(), -0.5, response_matrix.GetNcols()-0.5
////      ) ;
////  for ( int i=0; i<response_matrix.GetNrows(); i++ ) {
////     for ( int j=0; j<response_matrix.GetNcols(); j++ ) {
////        h_response_matrix_normalized -> SetBinContent( i+1, j+1,  response_matrix[i][j] ) ;
////     } // j
////  } // i



      TH1* h_gen_compare_ibu = (TH1*) h_gen_source -> Clone( "h_gen_compare_ibu" ) ;
      h_gen_compare_ibu -> Scale( ( hReco_ibu -> Integral() )/( h_gen_source -> Integral() ) ) ;

      TH1* h_gen_compare_tu  = (TH1*) h_gen_source -> Clone( "h_gen_compare_tu " ) ;
      h_gen_compare_tu  -> Scale( ( hReco_tu  -> Integral() )/( h_gen_source -> Integral() ) ) ;





      TH1* h_unfolding_err_ibu = (TH1*) hReco_ibu -> Clone( "h_unfolding_err_ibu" ) ;
      TH1* h_unfolding_err_tu  = (TH1*) hReco_tu  -> Clone( "h_unfolding_err_tu"  ) ;

      for ( int bi=1; bi <= hReco_ibu -> GetNbinsX(); bi++ ) {
         h_unfolding_err_ibu -> SetBinContent( bi , hReco_ibu -> GetBinError( bi ) ) ;
         h_unfolding_err_tu -> SetBinContent( bi , hReco_tu -> GetBinError( bi ) ) ;
         h_unfolding_err_ibu -> SetBinError( bi , 0. ) ;
         h_unfolding_err_tu -> SetBinError( bi , 0. ) ;
      }



      TH1* h_unfolding_result_diff = (TH1*) hReco_ibu -> Clone( "h_unfolding_result_diff" ) ;
      h_unfolding_result_diff -> Add( hReco_tu, -1. ) ;




      char htitle[1000] ;


      can_ru -> Divide(4,2) ;

      int ci = 1;




      can_ru -> cd( ci++ ) ;
      hReco_ibu -> SetName( "h_unfolded_ibu" ) ;
      sprintf( htitle, "IBU: %s", hReco_ibu -> GetTitle() ) ;
      hReco_ibu -> SetTitle( htitle ) ;
      hReco_ibu -> DrawCopy() ;

      h_gen_compare_ibu -> SetLineColor(2) ;
      h_gen_compare_ibu -> DrawCopy("hist same") ;




      can_ru -> cd( ci++ ) ;
      hReco_tu  -> SetName( "h_unfolded_tu " ) ;
      sprintf( htitle, "TUnfold: %s", hReco_tu -> GetTitle() ) ;
      hReco_tu -> SetTitle( htitle ) ;
      hReco_tu  -> DrawCopy() ;

      h_gen_compare_tu  -> SetLineColor(2) ;
      h_gen_compare_tu  -> DrawCopy("hist same") ;



      can_ru -> cd( ci++ ) ;
      h_global_correlation_coeff_tu -> SetLineColor(2) ;
      h_global_correlation_coeff_ibu -> DrawCopy() ;
      h_global_correlation_coeff_tu -> DrawCopy("same") ;


      can_ru -> cd( ci++ ) ;
      h_unfolding_err_tu -> SetLineColor(2) ;
      h_unfolding_err_ibu -> DrawCopy( "hist" ) ;
      h_unfolding_err_tu -> DrawCopy( "hist same" ) ;






      can_ru -> cd( ci++ ) ;
      h_unfold_cor_mat_ibu -> DrawCopy("colz") ;

      can_ru -> cd( ci++ ) ;
      h_unfold_cor_mat_tu  -> DrawCopy("colz") ;





      can_ru -> cd( ci++ ) ;
      h_unfolding_result_diff -> DrawCopy( "hist" ) ;


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

      printf("     roo_unfold_comp2(\"%s\",%d,%d,\"%s\")\n", rur_name, ngen, n_iter, input_file ) ;

      printf("\n\n\n") ;

   }












