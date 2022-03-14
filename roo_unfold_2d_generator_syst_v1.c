
#include "histio.c"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldTUnfold.h"

//
//
//  Before doing
//
//    .L roo_unfold_2d_generator_syst_v1.c
//
//  do
//
//   gSystem -> Load( "../RooUnfold/build/libRooUnfold.dylib" ) ;
//
//

//-------

TH2F* get_hist( const char* hname ) {
   TH2F* hp = (TH2F*) gDirectory -> FindObject( hname ) ;
   if ( hp == 0x0 ) { printf("\n\n *** can't find %s\n\n", hname ) ; gSystem->Exit(-1) ; }
   return hp ;
}

std::set<int> unused_global_bins ;

//-------

TH2F* trim_unused_bins( TH2F* hp, RooUnfoldResponse* rur ) {

   if ( unused_global_bins.size() <= 0 ) {
      TH1* h_gen_source = rur -> Htruth() ;
      for ( int i=1; i<=h_gen_source->GetNbinsX(); i++ ) {
         double x = h_gen_source->GetXaxis()->GetBinCenter( i ) ;
         for ( int j=1; j<=h_gen_source->GetNbinsY(); j++ ) {
            double y = h_gen_source->GetYaxis()->GetBinCenter( j ) ;
            float entries = h_gen_source -> GetBinContent( i, j ) ;
            if ( entries <= 0 ) {
               int global_bin = rur -> FindBin( (TH1*) h_gen_source, x, y ) ;
               printf("  trim_unused_bins:  %2d, %2d  (%9.5f, %9.5f)  global_bin %3d unused.\n", i, j, x, y, global_bin ) ;
               unused_global_bins.insert( global_bin ) ;
            }
         } // j
      } // i
   }

   int nb_input = hp -> GetNbinsX() ;
   int nb_trimmed = nb_input - unused_global_bins.size() ;

   char hname[100] ;
   sprintf( hname, "%s_trimmed", hp -> GetName() ) ;
   TH2F* hpr = new TH2F( hname, hp->GetTitle(), nb_trimmed, -0.5, nb_trimmed-0.5,  nb_trimmed, -0.5, nb_trimmed-0.5 ) ;
   int out_i = 1 ;
   for ( int in_i=1; in_i <= nb_input; in_i++ ) {
      //printf(" trim_unused_bins: in_i %3d\n", in_i ) ;
      ////if ( unused_global_bins.contains( in_i ) ) continue ; // NFG
      if ( unused_global_bins.find( in_i-1 ) != unused_global_bins.end() ) {
         //printf(" trim_unused_bins: skipping in_i %2d\n", in_i ) ;
         continue ;
      }
      int out_j = 1 ;
      for ( int in_j=1; in_j <= nb_input; in_j++ ) {
         //printf(" trim_unused_bins: in_j %3d\n", in_j ) ;
         ////if ( unused_global_bins.contains( in_j ) ) continue ; // NFG
         if ( unused_global_bins.find( in_j-1 ) != unused_global_bins.end() ) {
            //printf(" trim_unused_bins: skipping in_j %2d\n", in_j ) ;
            continue ;
         }
         //printf(" trim_unused_bins: Setting content for %2d, %2d to %9.5f\n", out_i, out_j, hp -> GetBinContent( in_i, in_j ) ) ;
         hpr -> SetBinContent( out_i, out_j, hp -> GetBinContent( in_i, in_j ) ) ;
         hpr -> SetBinError( out_i, out_j, hp -> GetBinError( in_i, in_j ) ) ;
         out_j ++ ;
      } // in_j
      out_i ++ ;
   } // in_i

   return hpr ;

}


//-------

TH1F* trim_unused_bins( TH1F* hp, RooUnfoldResponse* rur ) {

   if ( unused_global_bins.size() <= 0 ) {
      TH1* h_gen_source = rur -> Htruth() ;
      for ( int i=1; i<=h_gen_source->GetNbinsX(); i++ ) {
         double x = h_gen_source->GetXaxis()->GetBinCenter( i ) ;
         for ( int j=1; j<=h_gen_source->GetNbinsX(); j++ ) {
            double y = h_gen_source->GetYaxis()->GetBinCenter( j ) ;
            float entries = h_gen_source -> GetBinContent( i, j ) ;
            if ( entries <= 0 ) {
               int global_bin = rur -> FindBin( (TH1*) h_gen_source, x, y ) ;
               printf("  trim_unused_bins:  %2d, %2d  (%9.5f, %9.5f)  global_bin %3d unused.\n", i, j, x, y, global_bin ) ;
               unused_global_bins.insert( global_bin ) ;
            }
         } // j
      } // i
   }

   int nb_input = hp -> GetNbinsX() ;
   int nb_trimmed = nb_input - unused_global_bins.size() ;

   char hname[100] ;
   sprintf( hname, "%s_trimmed", hp -> GetName() ) ;
   TH1F* hpr = new TH1F( hname, hp->GetTitle(), nb_trimmed, -0.5, nb_trimmed-0.5  ) ;
   int out_i = 1 ;
   for ( int in_i=1; in_i <= nb_input; in_i++ ) {
      if ( unused_global_bins.find( in_i-1 ) != unused_global_bins.end() ) {
         continue ;
      }
      hpr -> SetBinContent( out_i, hp -> GetBinContent( in_i ) ) ;
      hpr -> SetBinError( out_i, hp -> GetBinError( in_i ) ) ;
      out_i ++ ;
   } // in_i

   return hpr ;

}

//-------

   void roo_unfold_2d_generator_syst_v1(
                     const char* rur_name = "rur_2D_log10_q2_vs_log10_x_dnn",
                     const char* input_file_a = "example-input-nbins_gen008_obs016.root",
                     const char* input_file_b = "django-input-nbins_gen008_obs016.root",
                     const char* input_name_a = "Rapgap",
                     const char* input_name_b = "Djangoh",
                     const char* method_name = "DNN",
                     int ngen = 1e7,
                     int method_index = 2,
                     int n_iter = 1000  // this is huge because I hacked RooUnfoldBayes to start with a flat prior.
                     ) {

      gStyle -> SetPalette( kBird ) ;
      gStyle -> SetPadRightMargin(0.15) ;

      TCanvas* can1 = (TCanvas*) gDirectory -> FindObject( "can1" ) ;
      if ( can1 == 0x0 ) { can1 = new TCanvas( "can1", "", 50, 50, 1200, 1200 ) ; }

      /////////TCanvas* can2 = (TCanvas*) gDirectory -> FindObject( "can2" ) ;
      /////////if ( can2 == 0x0 ) { can2 = new TCanvas( "can2", "", 1250, 50, 400, 800 ) ; }


      gStyle -> SetOptStat(0) ;

      TRandom3 tran(12345) ;

      gDirectory -> Delete( "*" ) ;

      TFile f_a( input_file_a, "read" ) ;
      printf("\n\n File A: %s  %s\n", input_file_a, input_name_a ) ;
      f_a.ls() ;

      TFile f_b( input_file_b, "read" ) ;
      printf("\n\n File B: %s  %s\n", input_file_b, input_name_b ) ;
      f_b.ls() ;


      char htitle[1000] ;


      RooUnfoldResponse* rur_a = (RooUnfoldResponse*) f_a.Get( rur_name ) ;
      if ( rur_a == 0x0 ) { printf("\n\n *** can't find %s in A\n\n", rur_name ) ; return ; }

      RooUnfoldResponse* rur_b = (RooUnfoldResponse*) f_b.Get( rur_name ) ;
      if ( rur_b == 0x0 ) { printf("\n\n *** can't find %s in B\n\n", rur_name ) ; return ; }





      TH2* h_in_gen_vs_obs_a = rur_a -> Hresponse() ;
      TH1* h_obs_source_a = rur_a -> Hmeasured() ;
      TH1* h_gen_source_a = rur_a -> Htruth() ;

      TH1* h_obs_random_a = (TH1*) h_obs_source_a -> Clone( "h_obs_random_a" ) ;
      h_obs_random_a->Reset() ;
      h_obs_random_a->FillRandom( h_obs_source_a, ngen ) ;

      sprintf( htitle, "Response matrix for log10 Q2 vs log10 x, %s", method_name ) ;
      h_in_gen_vs_obs_a -> SetTitle( htitle ) ;
      h_in_gen_vs_obs_a -> SetXTitle( "Reconstructed bin number" ) ;
      h_in_gen_vs_obs_a -> SetYTitle( "Gen bin number" ) ;







      TH2* h_in_gen_vs_obs_b = rur_b -> Hresponse() ;
      TH1* h_obs_source_b = rur_b -> Hmeasured() ;
      TH1* h_gen_source_b = rur_b -> Htruth() ;

      TH1* h_obs_random_b = (TH1*) h_obs_source_b -> Clone( "h_obs_random_b" ) ;
      h_obs_random_b->Reset() ;
      h_obs_random_b->FillRandom( h_obs_source_b, ngen ) ;

      sprintf( htitle, "Response matrix for log10 Q2 vs log10 x, %s", method_name ) ;
      h_in_gen_vs_obs_b -> SetTitle( htitle ) ;
      h_in_gen_vs_obs_b -> SetXTitle( "Reconstructed bin number" ) ;
      h_in_gen_vs_obs_b -> SetYTitle( "Gen bin number" ) ;





      rur_a -> UseOverflow(false) ;
      rur_b -> UseOverflow(false) ;




      RooUnfold* unfold_a_with_a ;
      RooUnfold* unfold_a_with_b ;
      RooUnfold* unfold_b_with_a ;
      RooUnfold* unfold_b_with_b ;

      if ( method_index == 1 ) {

         unfold_a_with_a = new RooUnfoldBayes( rur_a, h_obs_random_a, n_iter ) ;
         unfold_a_with_b = new RooUnfoldBayes( rur_b, h_obs_random_a, n_iter ) ;
         unfold_b_with_a = new RooUnfoldBayes( rur_a, h_obs_random_b, n_iter ) ;
         unfold_b_with_b = new RooUnfoldBayes( rur_b, h_obs_random_b, n_iter ) ;

      } else if ( method_index == 2 ) {

         unfold_a_with_a = new RooUnfoldTUnfold( rur_a, h_obs_random_a ) ;
         unfold_a_with_b = new RooUnfoldTUnfold( rur_b, h_obs_random_a ) ;
         unfold_b_with_a = new RooUnfoldTUnfold( rur_a, h_obs_random_b ) ;
         unfold_b_with_b = new RooUnfoldTUnfold( rur_b, h_obs_random_b ) ;

      } else {

         printf("\n\n *** I don't know method_index = %d\n\n", method_index ) ;
         return ;

      }


      TH1* hReco_a_with_a = (TH1*) unfold_a_with_a -> Hreco(RooUnfold::kCovariance) ;
      TH1* hReco_a_with_b = (TH1*) unfold_a_with_b -> Hreco(RooUnfold::kCovariance) ;
      TH1* hReco_b_with_a = (TH1*) unfold_b_with_a -> Hreco(RooUnfold::kCovariance) ;
      TH1* hReco_b_with_b = (TH1*) unfold_b_with_b -> Hreco(RooUnfold::kCovariance) ;


      sprintf( htitle, "Unfolded, log10 Q2 vs log10 x, %s, %s with %s", method_name, input_name_a, input_name_a ) ;
      hReco_a_with_a -> SetTitle( htitle ) ;
      hReco_a_with_a -> SetXTitle( "Unfolded log10 x" ) ;
      hReco_a_with_a -> SetYTitle( "Unfolded log10 Q2" ) ;

      sprintf( htitle, "Unfolded, log10 Q2 vs log10 x, %s, %s with %s", method_name, input_name_a, input_name_b ) ;
      hReco_a_with_b -> SetTitle( htitle ) ;
      hReco_a_with_b -> SetXTitle( "Unfolded log10 x" ) ;
      hReco_a_with_b -> SetYTitle( "Unfolded log10 Q2" ) ;

      sprintf( htitle, "Unfolded, log10 Q2 vs log10 x, %s, %s with %s", method_name, input_name_b, input_name_a ) ;
      hReco_b_with_a -> SetTitle( htitle ) ;
      hReco_b_with_a -> SetXTitle( "Unfolded log10 x" ) ;
      hReco_b_with_a -> SetYTitle( "Unfolded log10 Q2" ) ;

      sprintf( htitle, "Unfolded, log10 Q2 vs log10 x, %s, %s with %s", method_name, input_name_b, input_name_b ) ;
      hReco_b_with_b -> SetTitle( htitle ) ;
      hReco_b_with_b -> SetXTitle( "Unfolded log10 x" ) ;
      hReco_b_with_b -> SetYTitle( "Unfolded log10 Q2" ) ;






   ///printf("\n\n Unfolding results for A:\n") ;
   ///unfold_a -> PrintTable( cout ) ;
   ///printf("\n\n") ;

      TMatrixD unfolding_cov_mat_a_with_a = unfold_a_with_a -> Ereco() ;
      TMatrixD unfolding_inverse_cov_mat_a_with_a = unfold_a_with_a -> Wreco() ;
      TVectorD unfolding_err_a_with_a = unfold_a_with_a -> ErecoV( RooUnfold::kCovariance ) ;
      TVectorD unfolding_val_a_with_a = unfold_a_with_a -> Vreco() ;
      TVectorD gen_val_a_with_a( unfolding_val_a_with_a.GetNrows() ) ;

      TMatrixD unfolding_cov_mat_a_with_b = unfold_a_with_b -> Ereco() ;
      TMatrixD unfolding_inverse_cov_mat_a_with_b = unfold_a_with_b -> Wreco() ;
      TVectorD unfolding_err_a_with_b = unfold_a_with_b -> ErecoV( RooUnfold::kCovariance ) ;
      TVectorD unfolding_val_a_with_b = unfold_a_with_b -> Vreco() ;
      TVectorD gen_val_a_with_b( unfolding_val_a_with_b.GetNrows() ) ;

      TMatrixD unfolding_cov_mat_b_with_a = unfold_b_with_a -> Ereco() ;
      TMatrixD unfolding_inverse_cov_mat_b_with_a = unfold_b_with_a -> Wreco() ;
      TVectorD unfolding_err_b_with_a = unfold_b_with_a -> ErecoV( RooUnfold::kCovariance ) ;
      TVectorD unfolding_val_b_with_a = unfold_b_with_a -> Vreco() ;
      TVectorD gen_val_b_with_a( unfolding_val_b_with_a.GetNrows() ) ;

      TMatrixD unfolding_cov_mat_b_with_b = unfold_b_with_b -> Ereco() ;
      TMatrixD unfolding_inverse_cov_mat_b_with_b = unfold_b_with_b -> Wreco() ;
      TVectorD unfolding_err_b_with_b = unfold_b_with_b -> ErecoV( RooUnfold::kCovariance ) ;
      TVectorD unfolding_val_b_with_b = unfold_b_with_b -> Vreco() ;
      TVectorD gen_val_b_with_b( unfolding_val_b_with_b.GetNrows() ) ;









      sprintf( htitle, "Unfolding result error" ) ;
      TH1F* h_unfolding_result_err_a_with_a = new TH1F( "h_unfolding_result_err_a_with_a", htitle, unfolding_err_a_with_a.GetNrows(), -0.5, unfolding_err_a_with_a.GetNrows()-0.5 ) ;
      TH1F* h_unfolding_result_err_a_with_b = new TH1F( "h_unfolding_result_err_a_with_b", htitle, unfolding_err_a_with_b.GetNrows(), -0.5, unfolding_err_a_with_b.GetNrows()-0.5 ) ;
      TH1F* h_unfolding_result_err_b_with_a = new TH1F( "h_unfolding_result_err_b_with_a", htitle, unfolding_err_b_with_a.GetNrows(), -0.5, unfolding_err_b_with_a.GetNrows()-0.5 ) ;
      TH1F* h_unfolding_result_err_b_with_b = new TH1F( "h_unfolding_result_err_b_with_b", htitle, unfolding_err_b_with_b.GetNrows(), -0.5, unfolding_err_b_with_b.GetNrows()-0.5 ) ;
      for ( int i=0; i<unfolding_err_a_with_a.GetNrows(); i++ ) {
         h_unfolding_result_err_a_with_a -> SetBinContent( i+1, unfolding_err_a_with_a[i] ) ;
         h_unfolding_result_err_a_with_b -> SetBinContent( i+1, unfolding_err_a_with_b[i] ) ;
         h_unfolding_result_err_b_with_a -> SetBinContent( i+1, unfolding_err_b_with_a[i] ) ;
         h_unfolding_result_err_b_with_b -> SetBinContent( i+1, unfolding_err_b_with_b[i] ) ;
      }







      sprintf( htitle, "Unfolding cov. mat.  %s", method_name ) ;
      TH2F* h_unfold_cov_mat_a_with_a = new TH2F( "h_unfold_cov_mat_a_with_a", htitle,
          unfolding_cov_mat_a_with_a.GetNcols(), -0.5, unfolding_cov_mat_a_with_a.GetNcols()-0.5,
          unfolding_cov_mat_a_with_a.GetNcols(), -0.5, unfolding_cov_mat_a_with_a.GetNcols()-0.5 ) ;

      sprintf( htitle, "Unfolding cov. mat.  %s", method_name ) ;
      TH2F* h_unfold_cov_mat_a_with_b = new TH2F( "h_unfold_cov_mat_a_with_b", htitle,
          unfolding_cov_mat_a_with_b.GetNcols(), -0.5, unfolding_cov_mat_a_with_b.GetNcols()-0.5,
          unfolding_cov_mat_a_with_b.GetNcols(), -0.5, unfolding_cov_mat_a_with_b.GetNcols()-0.5 ) ;

      sprintf( htitle, "Unfolding cov. mat.  %s", method_name ) ;
      TH2F* h_unfold_cov_mat_b_with_a = new TH2F( "h_unfold_cov_mat_b_with_a", htitle,
          unfolding_cov_mat_b_with_a.GetNcols(), -0.5, unfolding_cov_mat_b_with_a.GetNcols()-0.5,
          unfolding_cov_mat_b_with_a.GetNcols(), -0.5, unfolding_cov_mat_b_with_a.GetNcols()-0.5 ) ;

      sprintf( htitle, "Unfolding cov. mat.  %s", method_name ) ;
      TH2F* h_unfold_cov_mat_b_with_b = new TH2F( "h_unfold_cov_mat_b_with_b", htitle,
          unfolding_cov_mat_b_with_b.GetNcols(), -0.5, unfolding_cov_mat_b_with_b.GetNcols()-0.5,
          unfolding_cov_mat_b_with_b.GetNcols(), -0.5, unfolding_cov_mat_b_with_b.GetNcols()-0.5 ) ;





      sprintf( htitle, "Unfolding cor. mat.  %s", method_name ) ;
      TH2F* h_unfold_cor_mat_a_with_a = new TH2F( "h_unfold_cor_mat_a_with_a", htitle,
          unfolding_cov_mat_a_with_a.GetNcols(), -0.5, unfolding_cov_mat_a_with_a.GetNcols()-0.5,
          unfolding_cov_mat_a_with_a.GetNcols(), -0.5, unfolding_cov_mat_a_with_a.GetNcols()-0.5 ) ;

      sprintf( htitle, "Unfolding cor. mat.  %s", method_name ) ;
      TH2F* h_unfold_cor_mat_a_with_b = new TH2F( "h_unfold_cor_mat_a_with_b", htitle,
          unfolding_cov_mat_a_with_b.GetNcols(), -0.5, unfolding_cov_mat_a_with_b.GetNcols()-0.5,
          unfolding_cov_mat_a_with_b.GetNcols(), -0.5, unfolding_cov_mat_a_with_b.GetNcols()-0.5 ) ;

      sprintf( htitle, "Unfolding cor. mat.  %s", method_name ) ;
      TH2F* h_unfold_cor_mat_b_with_a = new TH2F( "h_unfold_cor_mat_b_with_a", htitle,
          unfolding_cov_mat_b_with_a.GetNcols(), -0.5, unfolding_cov_mat_b_with_a.GetNcols()-0.5,
          unfolding_cov_mat_b_with_a.GetNcols(), -0.5, unfolding_cov_mat_b_with_a.GetNcols()-0.5 ) ;

      sprintf( htitle, "Unfolding cor. mat.  %s", method_name ) ;
      TH2F* h_unfold_cor_mat_b_with_b = new TH2F( "h_unfold_cor_mat_b_with_b", htitle,
          unfolding_cov_mat_b_with_b.GetNcols(), -0.5, unfolding_cov_mat_b_with_b.GetNcols()-0.5,
          unfolding_cov_mat_b_with_b.GetNcols(), -0.5, unfolding_cov_mat_b_with_b.GetNcols()-0.5 ) ;








      TH1F* h_global_correlation_coeff_a_with_a = new TH1F( "h_global_correlation_coeff_a_with_a", "Global cor. coef.",
         unfolding_cov_mat_a_with_a.GetNcols(), -0.5, unfolding_cov_mat_a_with_a.GetNcols()-0.5 ) ;
      TH1F* h_global_correlation_coeff_a_with_b = new TH1F( "h_global_correlation_coeff_a_with_b", "Global cor. coef.",
         unfolding_cov_mat_a_with_b.GetNcols(), -0.5, unfolding_cov_mat_a_with_b.GetNcols()-0.5 ) ;
      TH1F* h_global_correlation_coeff_b_with_a = new TH1F( "h_global_correlation_coeff_b_with_a", "Global cor. coef.",
         unfolding_cov_mat_b_with_a.GetNcols(), -0.5, unfolding_cov_mat_b_with_a.GetNcols()-0.5 ) ;
      TH1F* h_global_correlation_coeff_b_with_b = new TH1F( "h_global_correlation_coeff_b_with_b", "Global cor. coef.",
         unfolding_cov_mat_b_with_b.GetNcols(), -0.5, unfolding_cov_mat_b_with_b.GetNcols()-0.5 ) ;

      for ( int ri=0; ri<unfolding_cov_mat_a_with_a.GetNcols(); ri++ ) {

         int rbi = ri+1 ;

         for ( int ci=0; ci<unfolding_cov_mat_a_with_a.GetNcols(); ci++ ) {

            int cbi = ci+1 ;

            h_unfold_cov_mat_a_with_a -> SetBinContent( rbi, cbi, unfolding_cov_mat_a_with_a[ri][ci] ) ;
            h_unfold_cov_mat_a_with_b -> SetBinContent( rbi, cbi, unfolding_cov_mat_a_with_b[ri][ci] ) ;
            h_unfold_cov_mat_b_with_a -> SetBinContent( rbi, cbi, unfolding_cov_mat_b_with_a[ri][ci] ) ;
            h_unfold_cov_mat_b_with_b -> SetBinContent( rbi, cbi, unfolding_cov_mat_b_with_b[ri][ci] ) ;

            float rho = 1. ;
            float err_i, err_j ;

            err_i = 1. ;
            err_j = 1. ;
            if ( unfolding_cov_mat_a_with_a[ri][ri] > 0 ) err_i = sqrt( unfolding_cov_mat_a_with_a[ri][ri] ) ;
            if ( unfolding_cov_mat_a_with_a[ci][ci] > 0 ) err_j = sqrt( unfolding_cov_mat_a_with_a[ci][ci] ) ;
            if ( err_i > 0 && err_j > 0 ) rho = unfolding_cov_mat_a_with_a[ri][ci] / ( err_i * err_j ) ;
            h_unfold_cor_mat_a_with_a -> SetBinContent( rbi, cbi, rho ) ;

            err_i = 1. ;
            err_j = 1. ;
            if ( unfolding_cov_mat_a_with_b[ri][ri] > 0 ) err_i = sqrt( unfolding_cov_mat_a_with_b[ri][ri] ) ;
            if ( unfolding_cov_mat_a_with_b[ci][ci] > 0 ) err_j = sqrt( unfolding_cov_mat_a_with_b[ci][ci] ) ;
            if ( err_i > 0 && err_j > 0 ) rho = unfolding_cov_mat_a_with_b[ri][ci] / ( err_i * err_j ) ;
            h_unfold_cor_mat_a_with_b -> SetBinContent( rbi, cbi, rho ) ;

            err_i = 1. ;
            err_j = 1. ;
            if ( unfolding_cov_mat_b_with_a[ri][ri] > 0 ) err_i = sqrt( unfolding_cov_mat_b_with_a[ri][ri] ) ;
            if ( unfolding_cov_mat_b_with_a[ci][ci] > 0 ) err_j = sqrt( unfolding_cov_mat_b_with_a[ci][ci] ) ;
            if ( err_i > 0 && err_j > 0 ) rho = unfolding_cov_mat_b_with_a[ri][ci] / ( err_i * err_j ) ;
            h_unfold_cor_mat_b_with_a -> SetBinContent( rbi, cbi, rho ) ;

            err_i = 1. ;
            err_j = 1. ;
            if ( unfolding_cov_mat_b_with_b[ri][ri] > 0 ) err_i = sqrt( unfolding_cov_mat_b_with_b[ri][ri] ) ;
            if ( unfolding_cov_mat_b_with_b[ci][ci] > 0 ) err_j = sqrt( unfolding_cov_mat_b_with_b[ci][ci] ) ;
            if ( err_i > 0 && err_j > 0 ) rho = unfolding_cov_mat_b_with_b[ri][ci] / ( err_i * err_j ) ;
            h_unfold_cor_mat_b_with_b -> SetBinContent( rbi, cbi, rho ) ;


         } // ci

         float global_rho_a_with_a = 1. ;
         if ( unfolding_cov_mat_a_with_a[ri][ri] != 0 && unfolding_inverse_cov_mat_a_with_a[ri][ri] != 0 ) {
            global_rho_a_with_a = 0. ;
            float sqrt_arg = 1. - 1. / ( unfolding_cov_mat_a_with_a[ri][ri] * unfolding_inverse_cov_mat_a_with_a[ri][ri] ) ;
            if ( sqrt_arg > 0 ) global_rho_a_with_a = sqrt( sqrt_arg ) ;
         }
         h_global_correlation_coeff_a_with_a -> SetBinContent( rbi, global_rho_a_with_a ) ;

         float global_rho_a_with_b = 1. ;
         if ( unfolding_cov_mat_a_with_b[ri][ri] != 0 && unfolding_inverse_cov_mat_a_with_b[ri][ri] != 0 ) {
            global_rho_a_with_b = 0. ;
            float sqrt_arg = 1. - 1. / ( unfolding_cov_mat_a_with_b[ri][ri] * unfolding_inverse_cov_mat_a_with_b[ri][ri] ) ;
            if ( sqrt_arg > 0 ) global_rho_a_with_b = sqrt( sqrt_arg ) ;
         }
         h_global_correlation_coeff_a_with_b -> SetBinContent( rbi, global_rho_a_with_b ) ;

         float global_rho_b_with_a = 1. ;
         if ( unfolding_cov_mat_b_with_a[ri][ri] != 0 && unfolding_inverse_cov_mat_b_with_a[ri][ri] != 0 ) {
            global_rho_b_with_a = 0. ;
            float sqrt_arg = 1. - 1. / ( unfolding_cov_mat_b_with_a[ri][ri] * unfolding_inverse_cov_mat_b_with_a[ri][ri] ) ;
            if ( sqrt_arg > 0 ) global_rho_b_with_a = sqrt( sqrt_arg ) ;
         }
         h_global_correlation_coeff_b_with_a -> SetBinContent( rbi, global_rho_b_with_a ) ;

         float global_rho_b_with_b = 1. ;
         if ( unfolding_cov_mat_b_with_b[ri][ri] != 0 && unfolding_inverse_cov_mat_b_with_b[ri][ri] != 0 ) {
            global_rho_b_with_b = 0. ;
            float sqrt_arg = 1. - 1. / ( unfolding_cov_mat_b_with_b[ri][ri] * unfolding_inverse_cov_mat_b_with_b[ri][ri] ) ;
            if ( sqrt_arg > 0 ) global_rho_b_with_b = sqrt( sqrt_arg ) ;
         }
         h_global_correlation_coeff_b_with_b -> SetBinContent( rbi, global_rho_b_with_b ) ;


      } // ri

      /////h_unfold_cov_mat -> Print("all") ;





      TH1* h_gen_compare_a_with_a = (TH1*) h_gen_source_a -> Clone( "h_gen_compare_a_with_a" ) ;
      h_gen_compare_a_with_a -> Scale( ( hReco_a_with_a -> Integral() )/( h_gen_source_a -> Integral() ) ) ;

      TH1* h_gen_compare_a_with_b = (TH1*) h_gen_source_a -> Clone( "h_gen_compare_a_with_b" ) ;
      h_gen_compare_a_with_b -> Scale( ( hReco_a_with_b -> Integral() )/( h_gen_source_a -> Integral() ) ) ;

      TH1* h_gen_compare_b_with_a = (TH1*) h_gen_source_b -> Clone( "h_gen_compare_b_with_a" ) ;
      h_gen_compare_b_with_a -> Scale( ( hReco_b_with_a -> Integral() )/( h_gen_source_b -> Integral() ) ) ;

      TH1* h_gen_compare_b_with_b = (TH1*) h_gen_source_b -> Clone( "h_gen_compare_b_with_b" ) ;
      h_gen_compare_b_with_b -> Scale( ( hReco_b_with_b -> Integral() )/( h_gen_source_b -> Integral() ) ) ;



      {
         int vi = 0 ;
         for ( int ybi = 1; ybi <= h_gen_compare_a_with_a -> GetNbinsY(); ybi++ ) {
            for ( int xbi = 1; xbi <= h_gen_compare_a_with_a -> GetNbinsX(); xbi++ ) {

               gen_val_a_with_a[vi] = h_gen_compare_a_with_a -> GetBinContent( xbi, ybi ) ;
               gen_val_a_with_b[vi] = h_gen_compare_a_with_b -> GetBinContent( xbi, ybi ) ;
               gen_val_b_with_a[vi] = h_gen_compare_b_with_a -> GetBinContent( xbi, ybi ) ;
               gen_val_b_with_b[vi] = h_gen_compare_b_with_b -> GetBinContent( xbi, ybi ) ;

               vi ++ ;

            } // xbi
         } // ybi
      }


      printf("\n\n A with A: Unfolding results:\n") ;
      for ( int i=0; i<unfolding_val_a_with_a.GetNrows(); i++ ) {
         if ( unfolding_err_a_with_a[i] == 0 ) continue ;
         printf("A with A: global bin %3d :  gen = %12.1f, unfold =  %12.1f +/- %12.1f, diff = %12.1f,  diff/err = %12.3f  global correlation %12.3f\n",
            i, gen_val_a_with_a[i], unfolding_val_a_with_a[i], unfolding_err_a_with_a[i], (unfolding_val_a_with_a[i]-gen_val_a_with_a[i]),
            (unfolding_val_a_with_a[i]-gen_val_a_with_a[i])/unfolding_err_a_with_a[i],
            h_global_correlation_coeff_a_with_a->GetBinContent( i+1 ) ) ;
      }
      printf("\n\n\n") ;

      printf("\n\n A with B: Unfolding results:\n") ;
      for ( int i=0; i<unfolding_val_a_with_b.GetNrows(); i++ ) {
         if ( unfolding_err_a_with_b[i] == 0 ) continue ;
         printf("A with B: global bin %3d :  gen = %12.1f, unfold =  %12.1f +/- %12.1f, diff = %12.1f,  diff/err = %12.3f  global correlation %12.3f\n",
            i, gen_val_a_with_b[i], unfolding_val_a_with_b[i], unfolding_err_a_with_b[i], (unfolding_val_a_with_b[i]-gen_val_a_with_b[i]),
            (unfolding_val_a_with_b[i]-gen_val_a_with_b[i])/unfolding_err_a_with_b[i],
            h_global_correlation_coeff_a_with_b->GetBinContent( i+1 ) ) ;
      }
      printf("\n\n\n") ;

      printf("\n\n B with A: Unfolding results:\n") ;
      for ( int i=0; i<unfolding_val_b_with_a.GetNrows(); i++ ) {
         if ( unfolding_err_b_with_a[i] == 0 ) continue ;
         printf("B with A: global bin %3d :  gen = %12.1f, unfold =  %12.1f +/- %12.1f, diff = %12.1f,  diff/err = %12.3f  global correlation %12.3f\n",
            i, gen_val_b_with_a[i], unfolding_val_b_with_a[i], unfolding_err_b_with_a[i], (unfolding_val_b_with_a[i]-gen_val_b_with_a[i]),
            (unfolding_val_b_with_a[i]-gen_val_b_with_a[i])/unfolding_err_b_with_a[i],
            h_global_correlation_coeff_b_with_a->GetBinContent( i+1 ) ) ;
      }
      printf("\n\n\n") ;

      printf("\n\n B with B: Unfolding results:\n") ;
      for ( int i=0; i<unfolding_val_b_with_b.GetNrows(); i++ ) {
         if ( unfolding_err_b_with_b[i] == 0 ) continue ;
         printf("B with B: global bin %3d :  gen = %12.1f, unfold =  %12.1f +/- %12.1f, diff = %12.1f,  diff/err = %12.3f  global correlation %12.3f\n",
            i, gen_val_b_with_b[i], unfolding_val_b_with_b[i], unfolding_err_b_with_b[i], (unfolding_val_b_with_b[i]-gen_val_b_with_b[i]),
            (unfolding_val_b_with_b[i]-gen_val_b_with_b[i])/unfolding_err_b_with_b[i],
            h_global_correlation_coeff_b_with_b->GetBinContent( i+1 ) ) ;
      }
      printf("\n\n\n") ;







      sprintf( htitle, "1D view: unfolded  %s, %s with %s", method_name, input_name_a, input_name_a ) ;
      TH1F* h_1d_unfolded_val_a_with_a = new TH1F( "h_1d_unfolded_val_a_with_a", htitle, unfolding_val_a_with_a.GetNrows(), -0.5, unfolding_val_a_with_a.GetNrows()-0.5 ) ;
      TH1F* h_1d_gen_val_a_with_a = new TH1F( "h_1d_gen_val_a_with_a", "1D view: gen", unfolding_val_a_with_a.GetNrows(), -0.5, unfolding_val_a_with_a.GetNrows()-0.5 ) ;

      sprintf( htitle, "1D view: unfolded  %s, %s with %s", method_name, input_name_a, input_name_b ) ;
      TH1F* h_1d_unfolded_val_a_with_b = new TH1F( "h_1d_unfolded_val_a_with_b", htitle, unfolding_val_a_with_b.GetNrows(), -0.5, unfolding_val_a_with_b.GetNrows()-0.5 ) ;
      TH1F* h_1d_gen_val_a_with_b = new TH1F( "h_1d_gen_val_a_with_b", "1D view: gen", unfolding_val_a_with_a.GetNrows(), -0.5, unfolding_val_a_with_a.GetNrows()-0.5 ) ;

      sprintf( htitle, "1D view: unfolded  %s, %s with %s", method_name, input_name_b, input_name_a ) ;
      TH1F* h_1d_unfolded_val_b_with_a = new TH1F( "h_1d_unfolded_val_b_with_a", htitle, unfolding_val_b_with_a.GetNrows(), -0.5, unfolding_val_b_with_a.GetNrows()-0.5 ) ;
      TH1F* h_1d_gen_val_b_with_a = new TH1F( "h_1d_gen_val_b_with_a", "1D view: gen", unfolding_val_b_with_a.GetNrows(), -0.5, unfolding_val_b_with_a.GetNrows()-0.5 ) ;

      sprintf( htitle, "1D view: unfolded  %s, %s with %s", method_name, input_name_b, input_name_b ) ;
      TH1F* h_1d_unfolded_val_b_with_b = new TH1F( "h_1d_unfolded_val_b_with_b", htitle, unfolding_val_b_with_b.GetNrows(), -0.5, unfolding_val_b_with_b.GetNrows()-0.5 ) ;
      TH1F* h_1d_gen_val_b_with_b = new TH1F( "h_1d_gen_val_b_with_b", "1D view: gen", unfolding_val_b_with_b.GetNrows(), -0.5, unfolding_val_b_with_b.GetNrows()-0.5 ) ;





      for ( int i=0; i<unfolding_val_a_with_a.GetNrows(); i++ ) {

         h_1d_unfolded_val_a_with_a -> SetBinContent( i+1, unfolding_val_a_with_a[i] ) ;
         h_1d_unfolded_val_a_with_a -> SetBinError( i+1, unfolding_err_a_with_a[i] ) ;
         h_1d_gen_val_a_with_a -> SetBinContent( i+1, gen_val_a_with_a[i] ) ;

         h_1d_unfolded_val_a_with_b -> SetBinContent( i+1, unfolding_val_a_with_b[i] ) ;
         h_1d_unfolded_val_a_with_b -> SetBinError( i+1, unfolding_err_a_with_b[i] ) ;
         h_1d_gen_val_a_with_b -> SetBinContent( i+1, gen_val_a_with_b[i] ) ;

         h_1d_unfolded_val_b_with_a -> SetBinContent( i+1, unfolding_val_b_with_a[i] ) ;
         h_1d_unfolded_val_b_with_a -> SetBinError( i+1, unfolding_err_b_with_a[i] ) ;
         h_1d_gen_val_b_with_a -> SetBinContent( i+1, gen_val_b_with_a[i] ) ;

         h_1d_unfolded_val_b_with_b -> SetBinContent( i+1, unfolding_val_b_with_b[i] ) ;
         h_1d_unfolded_val_b_with_b -> SetBinError( i+1, unfolding_err_b_with_b[i] ) ;
         h_1d_gen_val_b_with_b -> SetBinContent( i+1, gen_val_b_with_b[i] ) ;

      }







      TH2F* h_unfold_cov_mat_trimmed_a_with_a = trim_unused_bins( h_unfold_cov_mat_a_with_a, rur_a ) ;
      TH2F* h_unfold_cor_mat_trimmed_a_with_a = trim_unused_bins( h_unfold_cor_mat_a_with_a, rur_a ) ;
      TH1F* h_global_correlation_coeff_trimmed_a_with_a = trim_unused_bins( h_global_correlation_coeff_a_with_a, rur_a ) ;
      TH1F* h_unfolding_result_err_trimmed_a_with_a = trim_unused_bins( h_unfolding_result_err_a_with_a, rur_a ) ;
      TH1F* h_1d_unfolded_val_trimmed_a_with_a = trim_unused_bins( h_1d_unfolded_val_a_with_a, rur_a ) ;
      TH1F* h_1d_gen_val_trimmed_a_with_a = trim_unused_bins( h_1d_gen_val_a_with_a, rur_a ) ;

      TH2F* h_unfold_cov_mat_trimmed_a_with_b = trim_unused_bins( h_unfold_cov_mat_a_with_b, rur_a ) ;
      TH2F* h_unfold_cor_mat_trimmed_a_with_b = trim_unused_bins( h_unfold_cor_mat_a_with_b, rur_a ) ;
      TH1F* h_global_correlation_coeff_trimmed_a_with_b = trim_unused_bins( h_global_correlation_coeff_a_with_b, rur_a ) ;
      TH1F* h_unfolding_result_err_trimmed_a_with_b = trim_unused_bins( h_unfolding_result_err_a_with_b, rur_a ) ;
      TH1F* h_1d_unfolded_val_trimmed_a_with_b = trim_unused_bins( h_1d_unfolded_val_a_with_b, rur_a ) ;
      TH1F* h_1d_gen_val_trimmed_a_with_b = trim_unused_bins( h_1d_gen_val_a_with_b, rur_a ) ;

      TH2F* h_unfold_cov_mat_trimmed_b_with_a = trim_unused_bins( h_unfold_cov_mat_b_with_a, rur_b ) ;
      TH2F* h_unfold_cor_mat_trimmed_b_with_a = trim_unused_bins( h_unfold_cor_mat_b_with_a, rur_b ) ;
      TH1F* h_global_correlation_coeff_trimmed_b_with_a = trim_unused_bins( h_global_correlation_coeff_b_with_a, rur_b ) ;
      TH1F* h_unfolding_result_err_trimmed_b_with_a = trim_unused_bins( h_unfolding_result_err_b_with_a, rur_b ) ;
      TH1F* h_1d_unfolded_val_trimmed_b_with_a = trim_unused_bins( h_1d_unfolded_val_b_with_a, rur_b ) ;
      TH1F* h_1d_gen_val_trimmed_b_with_a = trim_unused_bins( h_1d_gen_val_b_with_a, rur_b ) ;

      TH2F* h_unfold_cov_mat_trimmed_b_with_b = trim_unused_bins( h_unfold_cov_mat_b_with_b, rur_b ) ;
      TH2F* h_unfold_cor_mat_trimmed_b_with_b = trim_unused_bins( h_unfold_cor_mat_b_with_b, rur_b ) ;
      TH1F* h_global_correlation_coeff_trimmed_b_with_b = trim_unused_bins( h_global_correlation_coeff_b_with_b, rur_b ) ;
      TH1F* h_unfolding_result_err_trimmed_b_with_b = trim_unused_bins( h_unfolding_result_err_b_with_b, rur_b ) ;
      TH1F* h_1d_unfolded_val_trimmed_b_with_b = trim_unused_bins( h_1d_unfolded_val_b_with_b, rur_b ) ;
      TH1F* h_1d_gen_val_trimmed_b_with_b = trim_unused_bins( h_1d_gen_val_b_with_b, rur_b ) ;





      h_1d_unfolded_val_trimmed_a_with_a -> SetXTitle( "Global unfolded bin number") ;

      h_unfold_cov_mat_trimmed_a_with_a -> SetXTitle( "Global unfolded bin number") ;
      h_unfold_cov_mat_trimmed_a_with_a -> SetYTitle( "Global unfolded bin number") ;

      h_unfold_cor_mat_trimmed_a_with_a -> SetXTitle( "Global unfolded bin number") ;
      h_unfold_cor_mat_trimmed_a_with_a -> SetYTitle( "Global unfolded bin number") ;

      h_global_correlation_coeff_trimmed_a_with_a -> SetXTitle( "Global unfolded bin number") ;

      h_unfolding_result_err_trimmed_a_with_a -> SetXTitle( "Global unfolded bin number") ;

      h_global_correlation_coeff_trimmed_a_with_a -> SetLineWidth(2) ;

      h_unfolding_result_err_trimmed_a_with_a -> SetLineWidth(2) ;




      h_1d_unfolded_val_trimmed_a_with_b -> SetXTitle( "Global unfolded bin number") ;

      h_unfold_cov_mat_trimmed_a_with_b -> SetXTitle( "Global unfolded bin number") ;
      h_unfold_cov_mat_trimmed_a_with_b -> SetYTitle( "Global unfolded bin number") ;

      h_unfold_cor_mat_trimmed_a_with_b -> SetXTitle( "Global unfolded bin number") ;
      h_unfold_cor_mat_trimmed_a_with_b -> SetYTitle( "Global unfolded bin number") ;

      h_global_correlation_coeff_trimmed_a_with_b -> SetXTitle( "Global unfolded bin number") ;

      h_unfolding_result_err_trimmed_a_with_b -> SetXTitle( "Global unfolded bin number") ;

      h_global_correlation_coeff_trimmed_a_with_b -> SetLineWidth(2) ;

      h_unfolding_result_err_trimmed_a_with_b -> SetLineWidth(2) ;




      h_1d_unfolded_val_trimmed_b_with_a -> SetXTitle( "Global unfolded bin number") ;

      h_unfold_cov_mat_trimmed_b_with_a -> SetXTitle( "Global unfolded bin number") ;
      h_unfold_cov_mat_trimmed_b_with_a -> SetYTitle( "Global unfolded bin number") ;

      h_unfold_cor_mat_trimmed_b_with_a -> SetXTitle( "Global unfolded bin number") ;
      h_unfold_cor_mat_trimmed_b_with_a -> SetYTitle( "Global unfolded bin number") ;

      h_global_correlation_coeff_trimmed_b_with_a -> SetXTitle( "Global unfolded bin number") ;

      h_unfolding_result_err_trimmed_b_with_a -> SetXTitle( "Global unfolded bin number") ;

      h_global_correlation_coeff_trimmed_b_with_a -> SetLineWidth(2) ;

      h_unfolding_result_err_trimmed_b_with_a -> SetLineWidth(2) ;





      h_1d_unfolded_val_trimmed_b_with_b -> SetXTitle( "Global unfolded bin number") ;

      h_unfold_cov_mat_trimmed_b_with_b -> SetXTitle( "Global unfolded bin number") ;
      h_unfold_cov_mat_trimmed_b_with_b -> SetYTitle( "Global unfolded bin number") ;

      h_unfold_cor_mat_trimmed_b_with_b -> SetXTitle( "Global unfolded bin number") ;
      h_unfold_cor_mat_trimmed_b_with_b -> SetYTitle( "Global unfolded bin number") ;

      h_global_correlation_coeff_trimmed_b_with_b -> SetXTitle( "Global unfolded bin number") ;

      h_unfolding_result_err_trimmed_b_with_b -> SetXTitle( "Global unfolded bin number") ;

      h_global_correlation_coeff_trimmed_b_with_b -> SetLineWidth(2) ;

      h_unfolding_result_err_trimmed_b_with_b -> SetLineWidth(2) ;




      TH1F* h_unfolded_minus_gen_a_with_a = (TH1F*) h_1d_unfolded_val_trimmed_a_with_a -> Clone( "h_unfolded_minus_gen_a_with_a" ) ;
      h_unfolded_minus_gen_a_with_a -> Add( h_1d_gen_val_trimmed_a_with_a, -1. ) ;

      TH1F* h_unfolded_minus_gen_a_with_b = (TH1F*) h_1d_unfolded_val_trimmed_a_with_b -> Clone( "h_unfolded_minus_gen_a_with_b" ) ;
      h_unfolded_minus_gen_a_with_b -> Add( h_1d_gen_val_trimmed_a_with_b, -1. ) ;

      TH1F* h_unfolded_minus_gen_b_with_a = (TH1F*) h_1d_unfolded_val_trimmed_b_with_a -> Clone( "h_unfolded_minus_gen_b_with_a" ) ;
      h_unfolded_minus_gen_b_with_a -> Add( h_1d_gen_val_trimmed_b_with_a, -1. ) ;

      TH1F* h_unfolded_minus_gen_b_with_b = (TH1F*) h_1d_unfolded_val_trimmed_b_with_b -> Clone( "h_unfolded_minus_gen_b_with_b" ) ;
      h_unfolded_minus_gen_b_with_b -> Add( h_1d_gen_val_trimmed_b_with_b, -1. ) ;



      h_unfolded_minus_gen_a_with_a -> SetLineColor(2) ;
      h_unfolded_minus_gen_b_with_b -> SetLineColor(4) ;

      h_unfolded_minus_gen_a_with_b -> SetLineColor(2) ;
      h_unfolded_minus_gen_b_with_a -> SetLineColor(4) ;


      float hist_diff_max = 0 ;

      if ( h_unfolded_minus_gen_a_with_a -> GetMaximum() > hist_diff_max ) hist_diff_max = h_unfolded_minus_gen_a_with_a -> GetMaximum() ;
      if ( h_unfolded_minus_gen_a_with_b -> GetMaximum() > hist_diff_max ) hist_diff_max = h_unfolded_minus_gen_a_with_b -> GetMaximum() ;
      if ( h_unfolded_minus_gen_b_with_a -> GetMaximum() > hist_diff_max ) hist_diff_max = h_unfolded_minus_gen_b_with_a -> GetMaximum() ;
      if ( h_unfolded_minus_gen_b_with_b -> GetMaximum() > hist_diff_max ) hist_diff_max = h_unfolded_minus_gen_b_with_b -> GetMaximum() ;

      if ( fabs(h_unfolded_minus_gen_a_with_a -> GetMinimum()) > hist_diff_max ) hist_diff_max = fabs( h_unfolded_minus_gen_a_with_a -> GetMinimum() ) ;
      if ( fabs(h_unfolded_minus_gen_a_with_b -> GetMinimum()) > hist_diff_max ) hist_diff_max = fabs( h_unfolded_minus_gen_a_with_b -> GetMinimum() ) ;
      if ( fabs(h_unfolded_minus_gen_b_with_a -> GetMinimum()) > hist_diff_max ) hist_diff_max = fabs( h_unfolded_minus_gen_b_with_a -> GetMinimum() ) ;
      if ( fabs(h_unfolded_minus_gen_b_with_b -> GetMinimum()) > hist_diff_max ) hist_diff_max = fabs( h_unfolded_minus_gen_b_with_b -> GetMinimum() ) ;

      h_unfolded_minus_gen_a_with_a -> SetMaximum( 1.1 * hist_diff_max ) ;
      h_unfolded_minus_gen_a_with_b -> SetMaximum( 1.1 * hist_diff_max ) ;
      h_unfolded_minus_gen_b_with_a -> SetMaximum( 1.1 * hist_diff_max ) ;
      h_unfolded_minus_gen_b_with_b -> SetMaximum( 1.1 * hist_diff_max ) ;

      h_unfolded_minus_gen_a_with_a -> SetMinimum( -1.1 * hist_diff_max ) ;
      h_unfolded_minus_gen_a_with_b -> SetMinimum( -1.1 * hist_diff_max ) ;
      h_unfolded_minus_gen_b_with_a -> SetMinimum( -1.1 * hist_diff_max ) ;
      h_unfolded_minus_gen_b_with_b -> SetMinimum( -1.1 * hist_diff_max ) ;

      int ci ;

      can1 -> cd() ;
      can1 -> Clear() ;
      can1 -> Divide(2,3) ;

      ci = 1 ;

      can1 -> cd( ci++ ) ;
      h_1d_unfolded_val_trimmed_a_with_a -> DrawCopy() ;
      h_1d_gen_val_trimmed_a_with_a -> SetLineColor(2) ;
      h_1d_gen_val_trimmed_a_with_a -> DrawCopy("same hist" ) ;

      can1 -> cd( ci++ ) ;
      h_1d_unfolded_val_trimmed_a_with_b -> DrawCopy() ;
      h_1d_gen_val_trimmed_a_with_b -> SetLineColor(2) ;
      h_1d_gen_val_trimmed_a_with_b -> DrawCopy("same hist" ) ;

      can1 -> cd( ci++ ) ;
      h_1d_unfolded_val_trimmed_b_with_b -> DrawCopy() ;
      h_1d_gen_val_trimmed_b_with_b -> SetLineColor(2) ;
      h_1d_gen_val_trimmed_b_with_b -> DrawCopy("same hist" ) ;

      can1 -> cd( ci++ ) ;
      h_1d_unfolded_val_trimmed_b_with_a -> DrawCopy() ;
      h_1d_gen_val_trimmed_b_with_a -> SetLineColor(2) ;
      h_1d_gen_val_trimmed_b_with_a -> DrawCopy("same hist" ) ;



      can1 -> cd( ci++ ) ;
      h_unfolded_minus_gen_a_with_a -> DrawCopy("hist") ;
      h_unfolded_minus_gen_b_with_b -> DrawCopy("hist same") ;

      can1 -> cd( ci++ ) ;
      h_unfolded_minus_gen_a_with_b -> DrawCopy("hist") ;
      h_unfolded_minus_gen_b_with_a -> DrawCopy("hist same") ;



      printf("\n\n\n") ;

///   printf("   Cut and paste for this:\n\n") ;

///   printf("     roo_unfold_2d_compare_v1(\"%s\",\"%s\",%d,%d,%d,\"%s\",\"%s\",\"%s\")\n", rur_name_a, rur_name_b, ngen, method_index, n_iter, input_file, method_name_a, method_name_b ) ;

///   printf("\n\n\n") ;

      unused_global_bins.clear() ;


      saveHist("foo.root","*") ;

   }












