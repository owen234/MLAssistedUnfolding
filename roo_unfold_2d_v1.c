
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

   void roo_unfold_2d_v1( const char* rur_name = "rur_2D_log10_q2_vs_log10_x_dnn",
                     int ngen = 1e6,
                     int method_index = 2,
                     int n_iter = 1000,  // this is huge because I hacked RooUnfoldBayes to start with a flat prior.
                     const char* input_file = "example-input-nbins_gen008_obs016.root" ) {

      gStyle -> SetPadRightMargin(0.15) ;

      TCanvas* can_ru = (TCanvas*) gDirectory -> FindObject( "can_ru" ) ;
      if ( can_ru == 0x0 ) { can_ru = new TCanvas( "can_ru", "", 50, 50, 1800, 800 ) ; }
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


      rur -> UseOverflow(false) ;

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
      TH1* hReco = (TH1*) unfold -> Hreco(RooUnfold::kCovariance) ;
      //TH1D* hReco = (TH1D*) unfold -> Hreco(RooUnfold::kCovToy) ;


      unfold -> PrintTable( cout ) ;

      TMatrixD unfolding_cov_mat = unfold -> Ereco() ;

      TMatrixD unfolding_inverse_cov_mat = unfold -> Wreco() ;

      TVectorD unfolding_err = unfold -> ErecoV( RooUnfold::kCovariance ) ;

      TVectorD unfolding_val = unfold -> Vreco() ;


      TVectorD gen_val( unfolding_val.GetNrows() ) ;



      TH1F* h_unfolding_result_err = new TH1F( "h_unfolding_result_err", "Unfolding result error", unfolding_err.GetNrows(), -0.5, unfolding_err.GetNrows()-0.5 ) ;
      for ( int i=0; i<unfolding_err.GetNrows(); i++ ) {
         h_unfolding_result_err -> SetBinContent( i+1, unfolding_err[i] ) ;
      }



      //////unfolding_cov_mat.Print() ;

      TH2F* h_unfold_cov_mat = new TH2F( "h_unfold_cov_mat", "Unfolding cov. mat.",
          unfolding_cov_mat.GetNcols(), -0.5, unfolding_cov_mat.GetNcols()-0.5,
          unfolding_cov_mat.GetNcols(), -0.5, unfolding_cov_mat.GetNcols()-0.5 ) ;

      TH2F* h_unfold_cor_mat = new TH2F( "h_unfold_cor_mat", "Unfolding cor. mat.",
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
      //////response_matrix.Print() ;

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


      {
         int vi = 0 ;
         for ( int ybi = 1; ybi <= h_gen_compare -> GetNbinsY(); ybi++ ) {
            for ( int xbi = 1; xbi <= h_gen_compare -> GetNbinsX(); xbi++ ) {
               gen_val[vi] = h_gen_compare -> GetBinContent( xbi, ybi ) ;
               vi ++ ;
            } // xbi
         } // ybi
      }

      printf("\n\n Unfolding results:\n") ;
      for ( int i=0; i<unfolding_val.GetNrows(); i++ ) {
         if ( unfolding_err[i] == 0 ) continue ;
         printf("  global bin %3d :  gen = %9.3f, unfold =  %9.3f +/- %9.3f, diff = %9.3f,  diff/err = %9.3f  global correlation %9.3f\n",
            i, gen_val[i], unfolding_val[i], unfolding_err[i], (unfolding_val[i]-gen_val[i]),
            (unfolding_val[i]-gen_val[i])/unfolding_err[i],
            h_global_correlation_coeff->GetBinContent( i+1 ) ) ;
      }
      printf("\n\n\n") ;



      TH1F* h_1d_unfolded_val = new TH1F( "h_1d_unfolded_val", "1D view: unfolded", unfolding_val.GetNrows(), -0.5, unfolding_val.GetNrows()-0.5 ) ;
      TH1F* h_1d_gen_val = new TH1F( "h_1d_gen_val", "1D view: gen", unfolding_val.GetNrows(), -0.5, unfolding_val.GetNrows()-0.5 ) ;
      for ( int i=0; i<unfolding_val.GetNrows(); i++ ) {
         h_1d_unfolded_val -> SetBinContent( i+1, unfolding_val[i] ) ;
         h_1d_unfolded_val -> SetBinError( i+1, unfolding_err[i] ) ;
         h_1d_gen_val -> SetBinContent( i+1, gen_val[i] ) ;
      }



      TH1* h_diff = (TH1*) hReco -> Clone( "h_diff" ) ;
      h_diff -> SetTitle("Unfolded - Gen") ;
      h_diff -> Add( h_gen_compare, -1. ) ;


      TH1* h_diff_over_gen = (TH1*) h_diff -> Clone( "h_diff_over_gen" ) ;
      h_diff_over_gen -> SetTitle("(Unfolded-Gen)/Gen") ;
      h_diff_over_gen -> Divide( h_gen_compare ) ;


      TH1* h_diff_over_err = (TH1*) h_diff -> Clone( "h_diff_over_err" ) ;
      h_diff_over_err -> SetTitle("(Unfolded-Gen)/Error") ;
      for ( int xbi=1; xbi<=hReco->GetNbinsX(); xbi++ ) {
         for ( int ybi=1; ybi<=hReco->GetNbinsY(); ybi++ ) {
            float err = hReco->GetBinError( xbi, ybi ) ;
            float diff = h_diff -> GetBinContent( xbi, ybi ) ;
            float diff_over_err = 0. ;
            if ( err > 0 ) {
               diff_over_err = diff / err ;
            }
            h_diff_over_err -> SetBinContent( xbi, ybi, diff_over_err ) ;
            h_diff_over_err -> SetBinError( xbi, ybi, 0. ) ;
            /////printf("  h_diff_over_err :   %2d, %2d   (%9.3f,%9.3f)   %9.3f / %9.3f = %9.3f\n",
            /////  xbi, ybi, h_diff_over_gen->GetXaxis()->GetBinCenter(xbi), h_diff_over_gen->GetYaxis()->GetBinCenter(ybi),
            /////  diff, err, diff_over_err ) ;
         } // ybi
      } // xbi




      TH2F* h_unfold_cov_mat_trimmed = trim_unused_bins( h_unfold_cov_mat, rur ) ;
      TH2F* h_unfold_cor_mat_trimmed = trim_unused_bins( h_unfold_cor_mat, rur ) ;
      TH1F* h_global_correlation_coeff_trimmed = trim_unused_bins( h_global_correlation_coeff, rur ) ;
      TH1F* h_unfolding_result_err_trimmed = trim_unused_bins( h_unfolding_result_err, rur ) ;
      TH1F* h_1d_unfolded_val_trimmed = trim_unused_bins( h_1d_unfolded_val, rur ) ;
      TH1F* h_1d_gen_val_trimmed = trim_unused_bins( h_1d_gen_val, rur ) ;



      can_ru -> Divide(6,2) ;

      int ci = 1;

      can_ru -> cd( ci++ ) ;

      h_in_gen_vs_obs -> DrawCopy("colz") ;

      ////////////can_ru -> cd( ci++ ) ;
      ////////////h_response_matrix_normalized -> DrawCopy( "colz") ;

      can_ru -> cd( ci++ ) ;

      hReco -> SetName( "h_unfolded" ) ;
      hReco -> DrawCopy("colz") ;



      can_ru -> cd( ci++ ) ;
      h_gen_compare -> DrawCopy("colz") ;

      can_ru -> cd( ci++ ) ;
      h_diff -> DrawCopy("colz") ;

      can_ru -> cd( ci++ ) ;
      h_diff_over_gen -> SetMaximum(1.) ;
      h_diff_over_gen -> SetMinimum(-1.) ;
      h_diff_over_gen -> DrawCopy("colz") ;


      can_ru -> cd( ci++ ) ;
      h_diff_over_err -> DrawCopy("colz") ;

      can_ru -> cd( ci++ ) ;
      h_1d_unfolded_val_trimmed -> DrawCopy() ;
      h_1d_gen_val_trimmed->SetLineColor(2) ;
      h_1d_gen_val_trimmed->DrawCopy("same hist" ) ;


      can_ru -> cd( ci++ ) ;
      ////h_unfold_cov_mat -> DrawCopy("colz") ;
      h_unfold_cov_mat_trimmed -> DrawCopy("colz") ;


      can_ru -> cd( ci++ ) ;
      ///h_unfold_cor_mat -> SetMinimum(-1.) ;
      ///h_unfold_cor_mat -> SetMaximum( 1.) ;
      ///h_unfold_cor_mat -> DrawCopy("colz") ;
      h_unfold_cor_mat_trimmed -> SetMinimum(-1.) ;
      h_unfold_cor_mat_trimmed -> SetMaximum( 1.) ;
      h_unfold_cor_mat_trimmed -> DrawCopy("colz") ;


      can_ru -> cd( ci++ ) ;
      ////h_global_correlation_coeff -> DrawCopy() ;
      h_global_correlation_coeff_trimmed -> SetMaximum(1.1) ;
      h_global_correlation_coeff_trimmed -> DrawCopy("hist") ;


      can_ru -> cd( ci++ ) ;
      ////h_unfolding_result_err -> DrawCopy() ;
      h_unfolding_result_err_trimmed -> DrawCopy("hist") ;




      f.Close() ;


      gDirectory -> ls() ;


      printf("\n\n\n") ;

      printf("   Cut and paste for this:\n\n") ;

      printf("     roo_unfold_2d_v1(\"%s\",%d,%d,%d,\"%s\")\n", rur_name, ngen, method_index, n_iter, input_file ) ;

      printf("\n\n\n") ;

      unused_global_bins.clear() ;

   }












