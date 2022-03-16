

   void draw_2d_slices( TH1* hp, bool do_logy = false, const char* htitle = "" ) {

      char hname[100] ;

      int nbx = hp -> GetNbinsX() ;
      int nby = hp -> GetNbinsY() ;

      int first_xbin = 2 ;
      int last_xbin = nbx-1 ;

      int first_ybin = 2 ;
      int last_ybin = nby-1 ;

      float xaxis_min = hp -> GetXaxis() -> GetBinLowEdge( first_xbin ) ;
      float xaxis_max = hp -> GetXaxis() -> GetBinLowEdge( last_xbin + 1 ) ;
      int nbins_slice = last_xbin - first_xbin + 1 ;


      float hist_max = hp -> GetMaximum() ;
      printf("\n\n Hist max : %9.1f\n\n", hist_max ) ;

      gStyle -> SetLegendTextSize(0.035) ;
      TLegend* legend = new TLegend( 0.20, 0.63, 0.80, 0.87 ) ;


      for ( int ybi = first_ybin; ybi <= last_ybin; ybi ++ ) {

         sprintf( hname, "%s_yb%02d", hp -> GetName(), ybi ) ;
         printf("  Making %s\n", hname ) ;

         TH1F* h1dslice = new TH1F( hname, hname, nbins_slice, xaxis_min, xaxis_max ) ;

         for ( int xbi = first_xbin; xbi <= last_xbin; xbi ++ ) {
            float xaxis_val = hp -> GetXaxis() -> GetBinCenter( xbi ) ;
            int slice_bin = h1dslice -> FindBin( xaxis_val ) ;
            h1dslice -> SetBinContent( slice_bin, hp -> GetBinContent( xbi, ybi ) ) ;
            h1dslice -> SetBinError( slice_bin, hp -> GetBinError( xbi, ybi ) ) ;
            printf("  bin %2d : %9.1f +/- %5.1f | ", slice_bin, hp -> GetBinContent( xbi, ybi ), hp -> GetBinError( xbi, ybi ) ) ;
         } // xbi
         printf("\n") ;

         h1dslice -> SetLineColor( 11 + ybi ) ;
         h1dslice -> SetLineWidth( 2 ) ;
         h1dslice -> SetXTitle( hp -> GetXaxis() -> GetTitle() ) ;
         h1dslice -> SetYTitle( "Unfolded events" ) ;
         h1dslice -> SetMarkerStyle(0) ;
         h1dslice -> SetTitle( htitle ) ;
         h1dslice -> SetTitleOffset( 1.4, "y" ) ;
         if ( ybi == first_ybin ) {
            if ( do_logy ) {
               h1dslice -> SetMaximum( 340.2 * hist_max ) ;
               h1dslice -> SetMinimum( 1. ) ;
            } else {
               h1dslice -> SetMaximum( 1.7 * hist_max ) ;
            }
            h1dslice -> DrawCopy() ;
            h1dslice -> DrawCopy( "hist l same" ) ;
            if ( do_logy ) { gPad -> SetLogy(1) ; }
         } else {
            h1dslice -> DrawCopy("same") ;
            h1dslice -> DrawCopy("hist l same") ;
         }

         char label[100] ;
         float bin_log10q2min = hp -> GetYaxis() -> GetBinLowEdge( ybi ) ;
         float bin_log10q2max = hp -> GetYaxis() -> GetBinLowEdge( ybi+1 ) ;
         sprintf( label, "bin %2d : Q2 range [%4.0f, %4.0f]", ybi, pow(10.,bin_log10q2min), pow(10,bin_log10q2max)  ) ;
         legend -> AddEntry( h1dslice, label ) ;

      } // ybi

      legend -> Draw() ;




   }



