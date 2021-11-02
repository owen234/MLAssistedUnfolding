

   void draw_obs_in_gen_slices( TH2F* hp, bool make_canvas = true, int last_slice_bin = -1, int first_slice_bin = -1  ) {

      gStyle -> SetOptStat(0) ;

      if ( hp == 0x0 ) { printf("\n\n *** Null pointer!\n\n") ; return ; }

      int nbinsy = hp -> GetNbinsY() ;

      TH1* h_slices[100] ;

      float ymax = 0. ;
      for ( int ybi=1; ybi<=nbinsy; ybi++ ) {

         int hi = ybi-1;

         char hname[100] ;
         sprintf( hname, "%s_yb%02d", hp->GetName(), ybi ) ;
         h_slices[hi] = hp -> ProjectionX( hname, ybi, ybi ) ;

         float area = h_slices[hi] -> Integral() ;
         if ( area > 0 ) h_slices[hi] -> Scale( 1./area ) ;

         float hmax = h_slices[hi] -> GetMaximum() ;
         printf("   ybin %3d : hmax = %9.5f\n", ybi, hmax ) ;
         if ( h_slices[hi] -> GetMaximum() > ymax ) { ymax = hmax ; }

      } // ybi

      if ( make_canvas ) {
         TCanvas* can_doigs = (TCanvas*) gDirectory->FindObject("can_doigs") ;
         if ( can_doigs == 0x0 ) can_doigs = new TCanvas( "can_doigs", "", 50, 50, 800, 800 ) ;
         can_doigs -> Clear() ;
         can_doigs -> cd() ;
      }

      int last_slice = nbinsy ;
      if ( last_slice_bin > 0 ) {
         last_slice = last_slice_bin ;
      }

      int first_slice = 0 ;
      if ( first_slice_bin > 0 ) {
         first_slice = first_slice_bin ;
      }

      for ( int hi=first_slice; hi<last_slice; hi++ ) {
         h_slices[hi] -> SetLineColor( 12 + hi+1 ) ;
         h_slices[hi] -> SetLineWidth( 2 ) ;
         if ( hi == first_slice ) {
            h_slices[hi] -> SetMaximum( 1.1 * ymax ) ;
            h_slices[hi] -> Draw( "histc" ) ;
         } else {
            h_slices[hi] -> Draw( "histc same" ) ;
         }
      }




   }


