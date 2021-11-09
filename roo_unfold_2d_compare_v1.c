
#include "histio.c"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldTUnfold.h"

//
//
//  Before doing
//
//    .L roo_unfold_2d_v1.c
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

   void roo_unfold_2d_compare_v1(
                     const char* rur_name_a = "rur_2D_log10_q2_vs_log10_x_dnn",
                     const char* rur_name_b = "rur_2D_log10_q2_vs_log10_x_e",
                     int ngen = 1e5,
                     int method_index = 2,
                     int n_iter = 1000,  // this is huge because I hacked RooUnfoldBayes to start with a flat prior.
                     const char* input_file = "example-input-nbins_gen008_obs016.root",
                     const char* method_name_a = "DNN",
                     const char* method_name_b = "electron"
                     ) {

      gStyle -> SetPalette( kBird ) ;
      gStyle -> SetPadRightMargin(0.15) ;

      TCanvas* can1 = (TCanvas*) gDirectory -> FindObject( "can1" ) ;
      if ( can1 == 0x0 ) { can1 = new TCanvas( "can1", "", 50, 50, 1200, 1200 ) ; }

      TCanvas* can2 = (TCanvas*) gDirectory -> FindObject( "can2" ) ;
      if ( can2 == 0x0 ) { can2 = new TCanvas( "can2", "", 1250, 50, 400, 800 ) ; }


      gStyle -> SetOptStat(0) ;

      TRandom3 tran(12345) ;

      gDirectory -> Delete( "*" ) ;

      TFile f( input_file, "read" ) ;
      f.ls() ;


      char htitle[1000] ;


      RooUnfoldResponse* rur_a = (RooUnfoldResponse*) f.Get( rur_name_a ) ;
      if ( rur_a == 0x0 ) { printf("\n\n *** can't find %s\n\n", rur_name_a ) ; return ; }

      RooUnfoldResponse* rur_b = (RooUnfoldResponse*) f.Get( rur_name_b ) ;
      if ( rur_b == 0x0 ) { printf("\n\n *** can't find %s\n\n", rur_name_b ) ; return ; }





      TH2* h_in_gen_vs_obs_a = rur_a -> Hresponse() ;

      TH1* h_obs_source_a = rur_a -> Hmeasured() ;
      TH1* h_gen_source_a = rur_a -> Htruth() ;

      TH1* h_obs_random_a = (TH1*) h_obs_source_a -> Clone( "h_obs_random_a" ) ;
      h_obs_random_a->Reset() ;
      h_obs_random_a->FillRandom( h_obs_source_a, ngen ) ;

      sprintf( htitle, "Response matrix for log10 Q2 vs log10 x, %s", method_name_a ) ;
      h_in_gen_vs_obs_a -> SetTitle( htitle ) ;
      h_in_gen_vs_obs_a -> SetXTitle( "Reconstructed bin number" ) ;
      h_in_gen_vs_obs_a -> SetYTitle( "Gen bin number" ) ;




      TH2* h_in_gen_vs_obs_b = rur_b -> Hresponse() ;

      TH1* h_obs_source_b = rur_b -> Hmeasured() ;
      TH1* h_gen_source_b = rur_b -> Htruth() ;

      TH1* h_obs_random_b = (TH1*) h_obs_source_b -> Clone( "h_obs_random_b" ) ;
      h_obs_random_b->Reset() ;
      h_obs_random_b->FillRandom( h_obs_source_b, ngen ) ;

      sprintf( htitle, "Response matrix for log10 Q2 vs log10 x, %s", method_name_b ) ;
      h_in_gen_vs_obs_b -> SetTitle( htitle ) ;
      h_in_gen_vs_obs_b -> SetXTitle( "Reconstructed bin number" ) ;
      h_in_gen_vs_obs_b -> SetYTitle( "Gen bin number" ) ;





      rur_a -> UseOverflow(false) ;
      rur_b -> UseOverflow(false) ;




      RooUnfold* unfold_a ;
      RooUnfold* unfold_b ;
      if ( method_index == 1 ) {
         unfold_a = new RooUnfoldBayes( rur_a, h_obs_random_a, n_iter ) ;
         unfold_b = new RooUnfoldBayes( rur_b, h_obs_random_b, n_iter ) ;
      } else if ( method_index == 2 ) {
         unfold_a = new RooUnfoldTUnfold( rur_a, h_obs_random_a ) ;
         unfold_b = new RooUnfoldTUnfold( rur_b, h_obs_random_b ) ;
      } else {
         printf("\n\n *** I don't know method_index = %d\n\n", method_index ) ;
         return ;
      }


      TH1* hReco_a = (TH1*) unfold_a -> Hreco(RooUnfold::kCovariance) ;
      TH1* hReco_b = (TH1*) unfold_b -> Hreco(RooUnfold::kCovariance) ;

      sprintf( htitle, "Unfolded, log10 Q2 vs log10 x, %s", method_name_a ) ;
      hReco_a -> SetTitle( htitle ) ;
      hReco_a -> SetXTitle( "Unfolded log10 x" ) ;
      hReco_a -> SetYTitle( "Unfolded log10 Q2" ) ;


      sprintf( htitle, "Unfolded, log10 Q2 vs log10 x, %s", method_name_b ) ;
      hReco_b -> SetTitle( htitle ) ;
      hReco_b -> SetXTitle( "Unfolded log10 x" ) ;
      hReco_b -> SetYTitle( "Unfolded log10 Q2" ) ;




   ///printf("\n\n Unfolding results for A:\n") ;
   ///unfold_a -> PrintTable( cout ) ;
   ///printf("\n\n") ;

      TMatrixD unfolding_cov_mat_a = unfold_a -> Ereco() ;

      TMatrixD unfolding_inverse_cov_mat_a = unfold_a -> Wreco() ;

      TVectorD unfolding_err_a = unfold_a -> ErecoV( RooUnfold::kCovariance ) ;

      TVectorD unfolding_val_a = unfold_a -> Vreco() ;

      TVectorD gen_val_a( unfolding_val_a.GetNrows() ) ;





   ///printf("\n\n Unfolding results for B:\n") ;
   ///unfold_b -> PrintTable( cout ) ;
   ///printf("\n\n") ;

      TMatrixD unfolding_cov_mat_b = unfold_b -> Ereco() ;

      TMatrixD unfolding_inverse_cov_mat_b = unfold_b -> Wreco() ;

      TVectorD unfolding_err_b = unfold_b -> ErecoV( RooUnfold::kCovariance ) ;

      TVectorD unfolding_val_b = unfold_b -> Vreco() ;

      TVectorD gen_val_b( unfolding_val_b.GetNrows() ) ;








      sprintf( htitle, "Unfolding result error" ) ;
      TH1F* h_unfolding_result_err_a = new TH1F( "h_unfolding_result_err_a", htitle, unfolding_err_a.GetNrows(), -0.5, unfolding_err_a.GetNrows()-0.5 ) ;
      TH1F* h_unfolding_result_err_b = new TH1F( "h_unfolding_result_err_b", htitle, unfolding_err_b.GetNrows(), -0.5, unfolding_err_b.GetNrows()-0.5 ) ;
      for ( int i=0; i<unfolding_err_a.GetNrows(); i++ ) {
         h_unfolding_result_err_a -> SetBinContent( i+1, unfolding_err_a[i] ) ;
         h_unfolding_result_err_b -> SetBinContent( i+1, unfolding_err_b[i] ) ;
      }







      sprintf( htitle, "Unfolding cov. mat.  %s", method_name_a ) ;
      TH2F* h_unfold_cov_mat_a = new TH2F( "h_unfold_cov_mat_a", htitle,
          unfolding_cov_mat_a.GetNcols(), -0.5, unfolding_cov_mat_a.GetNcols()-0.5,
          unfolding_cov_mat_a.GetNcols(), -0.5, unfolding_cov_mat_a.GetNcols()-0.5 ) ;

      sprintf( htitle, "Unfolding cov. mat.  %s", method_name_b ) ;
      TH2F* h_unfold_cov_mat_b = new TH2F( "h_unfold_cov_mat_b", htitle,
          unfolding_cov_mat_b.GetNcols(), -0.5, unfolding_cov_mat_b.GetNcols()-0.5,
          unfolding_cov_mat_b.GetNcols(), -0.5, unfolding_cov_mat_b.GetNcols()-0.5 ) ;




      sprintf( htitle, "Unfolding cor. mat.  %s", method_name_a ) ;
      TH2F* h_unfold_cor_mat_a = new TH2F( "h_unfold_cor_mat_a", htitle,
          unfolding_cov_mat_a.GetNcols(), -0.5, unfolding_cov_mat_a.GetNcols()-0.5,
          unfolding_cov_mat_a.GetNcols(), -0.5, unfolding_cov_mat_a.GetNcols()-0.5 ) ;

      sprintf( htitle, "Unfolding cor. mat.  %s", method_name_b ) ;
      TH2F* h_unfold_cor_mat_b = new TH2F( "h_unfold_cor_mat_b", htitle,
          unfolding_cov_mat_b.GetNcols(), -0.5, unfolding_cov_mat_b.GetNcols()-0.5,
          unfolding_cov_mat_b.GetNcols(), -0.5, unfolding_cov_mat_b.GetNcols()-0.5 ) ;







      TH1F* h_global_correlation_coeff_a = new TH1F( "h_global_correlation_coeff_a", "Global cor. coef.",
         unfolding_cov_mat_a.GetNcols(), -0.5, unfolding_cov_mat_a.GetNcols()-0.5 ) ;

      TH1F* h_global_correlation_coeff_b = new TH1F( "h_global_correlation_coeff_b", "Global cor. coef.",
         unfolding_cov_mat_b.GetNcols(), -0.5, unfolding_cov_mat_b.GetNcols()-0.5 ) ;

      for ( int ri=0; ri<unfolding_cov_mat_a.GetNcols(); ri++ ) {
         int rbi = ri+1 ;
         for ( int ci=0; ci<unfolding_cov_mat_a.GetNcols(); ci++ ) {
            int cbi = ci+1 ;
            h_unfold_cov_mat_a -> SetBinContent( rbi, cbi, unfolding_cov_mat_a[ri][ci] ) ;
            h_unfold_cov_mat_b -> SetBinContent( rbi, cbi, unfolding_cov_mat_b[ri][ci] ) ;

            float rho = 1. ;
            float err_i, err_j ;

            err_i = 1. ;
            err_j = 1. ;
            if ( unfolding_cov_mat_a[ri][ri] > 0 ) err_i = sqrt( unfolding_cov_mat_a[ri][ri] ) ;
            if ( unfolding_cov_mat_a[ci][ci] > 0 ) err_j = sqrt( unfolding_cov_mat_a[ci][ci] ) ;
            if ( err_i > 0 && err_j > 0 ) rho = unfolding_cov_mat_a[ri][ci] / ( err_i * err_j ) ;
            h_unfold_cor_mat_a -> SetBinContent( rbi, cbi, rho ) ;

            err_i = 1. ;
            err_j = 1. ;
            if ( unfolding_cov_mat_b[ri][ri] > 0 ) err_i = sqrt( unfolding_cov_mat_b[ri][ri] ) ;
            if ( unfolding_cov_mat_b[ci][ci] > 0 ) err_j = sqrt( unfolding_cov_mat_b[ci][ci] ) ;
            if ( err_i > 0 && err_j > 0 ) rho = unfolding_cov_mat_b[ri][ci] / ( err_i * err_j ) ;
            h_unfold_cor_mat_b -> SetBinContent( rbi, cbi, rho ) ;

         } // ci

         float global_rho_a = 1. ;
         if ( unfolding_cov_mat_a[ri][ri] != 0 && unfolding_inverse_cov_mat_a[ri][ri] != 0 ) {
            global_rho_a = 0. ;
            float sqrt_arg = 1. - 1. / ( unfolding_cov_mat_a[ri][ri] * unfolding_inverse_cov_mat_a[ri][ri] ) ;
            if ( sqrt_arg > 0 ) global_rho_a = sqrt( sqrt_arg ) ;
         }
         h_global_correlation_coeff_a -> SetBinContent( rbi, global_rho_a ) ;

         float global_rho_b = 1. ;
         if ( unfolding_cov_mat_b[ri][ri] != 0 && unfolding_inverse_cov_mat_b[ri][ri] != 0 ) {
            global_rho_b = 0. ;
            float sqrt_arg = 1. - 1. / ( unfolding_cov_mat_b[ri][ri] * unfolding_inverse_cov_mat_b[ri][ri] ) ;
            if ( sqrt_arg > 0 ) global_rho_b = sqrt( sqrt_arg ) ;
         }
         h_global_correlation_coeff_b -> SetBinContent( rbi, global_rho_b ) ;

      } // ri

      /////h_unfold_cov_mat -> Print("all") ;





      TH1* h_gen_compare_a = (TH1*) h_gen_source_a -> Clone( "h_gen_compare_a" ) ;
      //////h_gen_compare_a -> Scale( ngen /( h_gen_source_a -> Integral() ) ) ;
      h_gen_compare_a -> Scale( ( hReco_a -> Integral() )/( h_gen_source_a -> Integral() ) ) ;

      TH1* h_gen_compare_b = (TH1*) h_gen_source_b -> Clone( "h_gen_compare_b" ) ;
      /////h_gen_compare_b -> Scale( ngen /( h_gen_source_b -> Integral() ) ) ;
      h_gen_compare_b -> Scale( ( hReco_b -> Integral() )/( h_gen_source_b -> Integral() ) ) ;


      {
         int vi = 0 ;
         for ( int ybi = 1; ybi <= h_gen_compare_a -> GetNbinsY(); ybi++ ) {
            for ( int xbi = 1; xbi <= h_gen_compare_a -> GetNbinsX(); xbi++ ) {
               gen_val_a[vi] = h_gen_compare_a -> GetBinContent( xbi, ybi ) ;
               gen_val_b[vi] = h_gen_compare_b -> GetBinContent( xbi, ybi ) ;
               vi ++ ;
            } // xbi
         } // ybi
      }


      printf("\n\n A: Unfolding results:\n") ;
      for ( int i=0; i<unfolding_val_a.GetNrows(); i++ ) {
         if ( unfolding_err_a[i] == 0 ) continue ;
         printf("A: global bin %3d :  gen = %9.3f, unfold =  %9.3f +/- %9.3f, diff = %9.3f,  diff/err = %9.3f  global correlation %9.3f\n",
            i, gen_val_a[i], unfolding_val_a[i], unfolding_err_a[i], (unfolding_val_a[i]-gen_val_a[i]),
            (unfolding_val_a[i]-gen_val_a[i])/unfolding_err_a[i],
            h_global_correlation_coeff_a->GetBinContent( i+1 ) ) ;
      }
      printf("\n\n\n") ;

      printf("\n\n B: Unfolding results:\n") ;
      for ( int i=0; i<unfolding_val_b.GetNrows(); i++ ) {
         if ( unfolding_err_b[i] == 0 ) continue ;
         printf("A: global bin %3d :  gen = %9.3f, unfold =  %9.3f +/- %9.3f, diff = %9.3f,  diff/err = %9.3f  global correlation %9.3f\n",
            i, gen_val_b[i], unfolding_val_b[i], unfolding_err_b[i], (unfolding_val_b[i]-gen_val_b[i]),
            (unfolding_val_b[i]-gen_val_b[i])/unfolding_err_b[i],
            h_global_correlation_coeff_b->GetBinContent( i+1 ) ) ;
      }
      printf("\n\n\n") ;







      sprintf( htitle, "1D view: unfolded  %s", method_name_a ) ;
      TH1F* h_1d_unfolded_val_a = new TH1F( "h_1d_unfolded_val_a", htitle, unfolding_val_a.GetNrows(), -0.5, unfolding_val_a.GetNrows()-0.5 ) ;
      TH1F* h_1d_gen_val_a = new TH1F( "h_1d_gen_val_a", "1D view: gen", unfolding_val_a.GetNrows(), -0.5, unfolding_val_a.GetNrows()-0.5 ) ;


      sprintf( htitle, "1D view: unfolded  %s", method_name_b ) ;
      TH1F* h_1d_unfolded_val_b = new TH1F( "h_1d_unfolded_val_b", htitle, unfolding_val_b.GetNrows(), -0.5, unfolding_val_b.GetNrows()-0.5 ) ;
      TH1F* h_1d_gen_val_b = new TH1F( "h_1d_gen_val_b", "1D view: gen", unfolding_val_b.GetNrows(), -0.5, unfolding_val_b.GetNrows()-0.5 ) ;




      for ( int i=0; i<unfolding_val_a.GetNrows(); i++ ) {

         h_1d_unfolded_val_a -> SetBinContent( i+1, unfolding_val_a[i] ) ;
         h_1d_unfolded_val_a -> SetBinError( i+1, unfolding_err_a[i] ) ;
         h_1d_gen_val_a -> SetBinContent( i+1, gen_val_a[i] ) ;

         h_1d_unfolded_val_b -> SetBinContent( i+1, unfolding_val_b[i] ) ;
         h_1d_unfolded_val_b -> SetBinError( i+1, unfolding_err_b[i] ) ;
         h_1d_gen_val_b -> SetBinContent( i+1, gen_val_b[i] ) ;

      }



//    TH1* h_diff = (TH1*) hReco -> Clone( "h_diff" ) ;
//    h_diff -> SetTitle("Unfolded - Gen") ;
//    h_diff -> Add( h_gen_compare, -1. ) ;


//    TH1* h_diff_over_gen = (TH1*) h_diff -> Clone( "h_diff_over_gen" ) ;
//    h_diff_over_gen -> SetTitle("(Unfolded-Gen)/Gen") ;
//    h_diff_over_gen -> Divide( h_gen_compare ) ;


//    TH1* h_diff_over_err = (TH1*) h_diff -> Clone( "h_diff_over_err" ) ;
//    h_diff_over_err -> SetTitle("(Unfolded-Gen)/Error") ;
//    for ( int xbi=1; xbi<=hReco->GetNbinsX(); xbi++ ) {
//       for ( int ybi=1; ybi<=hReco->GetNbinsY(); ybi++ ) {
//          float err = hReco->GetBinError( xbi, ybi ) ;
//          float diff = h_diff -> GetBinContent( xbi, ybi ) ;
//          float diff_over_err = 0. ;
//          if ( err > 0 ) {
//             diff_over_err = diff / err ;
//          }
//          h_diff_over_err -> SetBinContent( xbi, ybi, diff_over_err ) ;
//          h_diff_over_err -> SetBinError( xbi, ybi, 0. ) ;
//          /////printf("  h_diff_over_err :   %2d, %2d   (%9.3f,%9.3f)   %9.3f / %9.3f = %9.3f\n",
//          /////  xbi, ybi, h_diff_over_gen->GetXaxis()->GetBinCenter(xbi), h_diff_over_gen->GetYaxis()->GetBinCenter(ybi),
//          /////  diff, err, diff_over_err ) ;
//       } // ybi
//    } // xbi




      TH2F* h_unfold_cov_mat_trimmed_a = trim_unused_bins( h_unfold_cov_mat_a, rur_a ) ;
      TH2F* h_unfold_cor_mat_trimmed_a = trim_unused_bins( h_unfold_cor_mat_a, rur_a ) ;
      TH1F* h_global_correlation_coeff_trimmed_a = trim_unused_bins( h_global_correlation_coeff_a, rur_a ) ;
      TH1F* h_unfolding_result_err_trimmed_a = trim_unused_bins( h_unfolding_result_err_a, rur_a ) ;
      TH1F* h_1d_unfolded_val_trimmed_a = trim_unused_bins( h_1d_unfolded_val_a, rur_a ) ;
      TH1F* h_1d_gen_val_trimmed_a = trim_unused_bins( h_1d_gen_val_a, rur_a ) ;

      TH2F* h_unfold_cov_mat_trimmed_b = trim_unused_bins( h_unfold_cov_mat_b, rur_b ) ;
      TH2F* h_unfold_cor_mat_trimmed_b = trim_unused_bins( h_unfold_cor_mat_b, rur_b ) ;
      TH1F* h_global_correlation_coeff_trimmed_b = trim_unused_bins( h_global_correlation_coeff_b, rur_b ) ;
      TH1F* h_unfolding_result_err_trimmed_b = trim_unused_bins( h_unfolding_result_err_b, rur_b ) ;
      TH1F* h_1d_unfolded_val_trimmed_b = trim_unused_bins( h_1d_unfolded_val_b, rur_b ) ;
      TH1F* h_1d_gen_val_trimmed_b = trim_unused_bins( h_1d_gen_val_b, rur_b ) ;


      h_1d_unfolded_val_trimmed_a -> SetXTitle( "Global unfolded bin number") ;
      h_1d_unfolded_val_trimmed_b -> SetXTitle( "Global unfolded bin number") ;

      h_unfold_cov_mat_trimmed_a -> SetXTitle( "Global unfolded bin number") ;
      h_unfold_cov_mat_trimmed_a -> SetYTitle( "Global unfolded bin number") ;
      h_unfold_cov_mat_trimmed_b -> SetXTitle( "Global unfolded bin number") ;
      h_unfold_cov_mat_trimmed_b -> SetYTitle( "Global unfolded bin number") ;

      h_unfold_cor_mat_trimmed_a -> SetXTitle( "Global unfolded bin number") ;
      h_unfold_cor_mat_trimmed_a -> SetYTitle( "Global unfolded bin number") ;
      h_unfold_cor_mat_trimmed_b -> SetXTitle( "Global unfolded bin number") ;
      h_unfold_cor_mat_trimmed_b -> SetYTitle( "Global unfolded bin number") ;

      h_global_correlation_coeff_trimmed_a -> SetXTitle( "Global unfolded bin number") ;
      h_global_correlation_coeff_trimmed_b -> SetXTitle( "Global unfolded bin number") ;

      h_unfolding_result_err_trimmed_a -> SetXTitle( "Global unfolded bin number") ;
      h_unfolding_result_err_trimmed_b -> SetXTitle( "Global unfolded bin number") ;



      h_global_correlation_coeff_trimmed_a -> SetLineWidth(2) ;
      h_global_correlation_coeff_trimmed_b -> SetLineWidth(2) ;

      h_unfolding_result_err_trimmed_a -> SetLineWidth(2) ;
      h_unfolding_result_err_trimmed_b -> SetLineWidth(2) ;


      TH1F* h_unfolding_result_error_ratio = (TH1F*) h_unfolding_result_err_trimmed_a -> Clone( "h_unfolding_result_error_ratio" ) ;
      h_unfolding_result_error_ratio -> Divide( h_unfolding_result_err_trimmed_b ) ;
      sprintf( htitle, "Unfolding error ratio %s / %s", method_name_a, method_name_b ) ;
      h_unfolding_result_error_ratio -> SetTitle( htitle ) ;
      h_unfolding_result_error_ratio -> SetYTitle( htitle ) ;
      h_unfolding_result_error_ratio -> SetXTitle( "Global unfolded bin number") ;




      int ci ;

      can1 -> cd() ;
      can1 -> Clear() ;
      can1 -> Divide(3,3) ;

      ci = 1 ;


    //--- main plots for A

      can1 -> cd( ci++ ) ;
      h_in_gen_vs_obs_a -> DrawCopy("colz") ;

      can1 -> cd( ci++ ) ;
      hReco_a -> DrawCopy( "colz" ) ;

      can1 -> cd( ci++ ) ;
      h_1d_unfolded_val_trimmed_a -> DrawCopy() ;
      h_1d_gen_val_trimmed_a -> SetLineColor(2) ;
      h_1d_gen_val_trimmed_a -> DrawCopy("same hist") ;



    //--- main plots for B

      can1 -> cd( ci++ ) ;
      h_in_gen_vs_obs_b -> DrawCopy("colz") ;

      can1 -> cd( ci++ ) ;
      hReco_b -> DrawCopy( "colz" ) ;

      can1 -> cd( ci++ ) ;
      h_1d_unfolded_val_trimmed_b -> DrawCopy() ;
      h_1d_gen_val_trimmed_b -> SetLineColor(2) ;
      h_1d_gen_val_trimmed_b -> DrawCopy("same hist") ;



    //--- plots with both superimposed

      can1 -> cd( ci++ ) ;
      h_global_correlation_coeff_trimmed_a -> SetMaximum(1.1) ;
      h_global_correlation_coeff_trimmed_a -> SetLineColor(4) ;
      h_global_correlation_coeff_trimmed_b -> SetLineColor(2) ;
      h_global_correlation_coeff_trimmed_a -> DrawCopy("hist") ;
      h_global_correlation_coeff_trimmed_b -> DrawCopy("same hist") ;
      gPad -> SetGridy(1) ;



      can1 -> cd( ci++ ) ;
      h_unfolding_result_err_trimmed_a -> SetLineColor(4) ;
      h_unfolding_result_err_trimmed_b -> SetLineColor(2) ;
      h_unfolding_result_err_trimmed_b -> DrawCopy("hist") ;
      h_unfolding_result_err_trimmed_a -> DrawCopy("hist same") ;
      gPad -> SetGridy(1) ;


      can1 -> cd( ci++ ) ;
      h_unfolding_result_error_ratio -> SetMaximum( 1.1 ) ;
      h_unfolding_result_error_ratio -> DrawCopy( "hist" ) ;
      gPad -> SetGridy(1) ;



      can1 -> Update() ;
      can1 -> Draw() ;
      gSystem -> ProcessEvents() ;


    //===============

      can2 -> cd() ;
      can2 -> Clear() ;
      can2 -> Divide(1,2) ;

        Int_t nb = 90 ;

        const Int_t Number = 3 ;
        Double_t Length[Number] = {0.,0.5, 1.} ;
        Double_t Red[Number] = {0.,1.,1.} ;
        Double_t Green[Number] = {0.,1.,0.} ;
        Double_t Blue[Number] = {1.,1.,0.} ;
        TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);

      h_unfold_cor_mat_trimmed_a -> SetContour( nb ) ;
      h_unfold_cor_mat_trimmed_b -> SetContour( nb ) ;


      ci = 1 ;

      can2 -> cd( ci++ ) ;
      h_unfold_cor_mat_trimmed_a -> SetMinimum(-1.) ;
      h_unfold_cor_mat_trimmed_a -> SetMaximum( 1.) ;
      h_unfold_cor_mat_trimmed_a -> DrawCopy("colz") ;

      can2 -> cd( ci++ ) ;
      h_unfold_cor_mat_trimmed_b -> SetMinimum(-1.) ;
      h_unfold_cor_mat_trimmed_b -> SetMaximum( 1.) ;
      h_unfold_cor_mat_trimmed_b -> DrawCopy("colz") ;



//    can_ru -> Divide(6,2) ;

//    int ci = 1;

//    can_ru -> cd( ci++ ) ;

//    h_in_gen_vs_obs -> DrawCopy("colz") ;

//    ////////////can_ru -> cd( ci++ ) ;
//    ////////////h_response_matrix_normalized -> DrawCopy( "colz") ;

//    can_ru -> cd( ci++ ) ;

//    hReco -> SetName( "h_unfolded" ) ;
//    hReco -> DrawCopy("colz") ;



//    can_ru -> cd( ci++ ) ;
//    h_gen_compare -> DrawCopy("colz") ;

//    can_ru -> cd( ci++ ) ;
//    h_diff -> DrawCopy("colz") ;

//    can_ru -> cd( ci++ ) ;
//    h_diff_over_gen -> SetMaximum(1.) ;
//    h_diff_over_gen -> SetMinimum(-1.) ;
//    h_diff_over_gen -> DrawCopy("colz") ;


//    can_ru -> cd( ci++ ) ;
//    h_diff_over_err -> DrawCopy("colz") ;

//    can_ru -> cd( ci++ ) ;
//    h_1d_unfolded_val_trimmed -> DrawCopy() ;
//    h_1d_gen_val_trimmed->SetLineColor(2) ;
//    h_1d_gen_val_trimmed->DrawCopy("same hist" ) ;


//    can_ru -> cd( ci++ ) ;
//    ////h_unfold_cov_mat -> DrawCopy("colz") ;
//    h_unfold_cov_mat_trimmed -> DrawCopy("colz") ;


//    can_ru -> cd( ci++ ) ;
//    ///h_unfold_cor_mat -> SetMinimum(-1.) ;
//    ///h_unfold_cor_mat -> SetMaximum( 1.) ;
//    ///h_unfold_cor_mat -> DrawCopy("colz") ;
//    h_unfold_cor_mat_trimmed -> SetMinimum(-1.) ;
//    h_unfold_cor_mat_trimmed -> SetMaximum( 1.) ;
//    h_unfold_cor_mat_trimmed -> DrawCopy("colz") ;


//    can_ru -> cd( ci++ ) ;
//    ////h_global_correlation_coeff -> DrawCopy() ;
//    h_global_correlation_coeff_trimmed -> SetMaximum(1.1) ;
//    h_global_correlation_coeff_trimmed -> DrawCopy("hist") ;


//    can_ru -> cd( ci++ ) ;
//    ////h_unfolding_result_err -> DrawCopy() ;
//    h_unfolding_result_err_trimmed -> DrawCopy("hist") ;




//    f.Close() ;


//    gDirectory -> ls() ;


      printf("\n\n\n") ;

      printf("   Cut and paste for this:\n\n") ;

      printf("     roo_unfold_2d_compare_v1(\"%s\",\"%s\",%d,%d,%d,\"%s\",\"%s\",\"%s\")\n", rur_name_a, rur_name_b, ngen, method_index, n_iter, input_file, method_name_a, method_name_b ) ;

      printf("\n\n\n") ;

      unused_global_bins.clear() ;

   }












