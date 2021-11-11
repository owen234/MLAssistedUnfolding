


#include "histio.c"

#ifndef helpers
#define helpers

TH2F* get_hist( const char* hname ) {
   TH2F* hp = (TH2F*) gDirectory -> FindObject( hname ) ;
   if ( hp == 0x0 ) { printf("\n\n *** can't find %s\n\n", hname ) ; gSystem->Exit(-1) ; }
   return hp ;
}

void SetupCorrelationPalette() {

      static Bool_t initialized = kFALSE ;
      static Int_t colors[90] ;

      const Int_t Number = 3 ;
      Double_t Length[Number] = {0.,0.5, 1.} ;
      Double_t Red[Number] = {0.,1.,1.} ;
      Double_t Green[Number] = {0.,1.,0.} ;
      Double_t Blue[Number] = {1.,1.,0.} ;

      if ( !initialized ) {
         Int_t fi = TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,90);
         for ( int i=0; i<90; i++ ) colors[i] = fi+i ;
         initialized = kTRUE ;
         return ;
      }
      gStyle -> SetPalette( 90, colors ) ;

}

void Setup2DhistPalette() {
      gStyle -> SetPalette( kBird ) ;
}

#endif

//----------

   void paper_plots_dis_gen_response_comp(
                      const char* hist_name = "h_log10_x_gen_vs_obs_dnn",
                      const char* method_name = "DNN",
                      const char* var_name = "x",
                      float max_diff = -1.,
                      const char* input_file_a = "paper-plots-input-1D-nbins_gen010_obs020.root",
                      const char* input_file_b = "paper-plots-input-1D-nbins_gen010_obs020-django.root",
                      const char* input_name_a = "Rapgap",
                      const char* input_name_b = "Djangoh"
                      ) {

      char htitle[1000] ;
      char hname[1000] ;
      char fname[1000] ;

      gSystem -> Exec( "mkdir -p paper-plots" ) ;

      gDirectory -> Delete( "h*" ) ;
      gStyle -> SetPalette( kBird ) ;

      gStyle -> SetTitleX(0.05) ;

      loadHist( input_file_a, input_name_a ) ;
      loadHist( input_file_b, input_name_b ) ;

      printf("\n\n") ;
      gDirectory -> ls() ;
      printf("\n\n") ;

      sprintf( hname, "%s_%s", hist_name, input_name_a ) ;
      TH2F* h_in_gen_vs_obs_a = get_hist( hname ) ;

      sprintf( hname, "%s_%s", hist_name, input_name_b ) ;
      TH2F* h_in_gen_vs_obs_b = get_hist( hname ) ;



      TH2F* h_normalized_response_a = (TH2F*) h_in_gen_vs_obs_a -> Clone( "h_normalized_response_a" ) ;
      for ( int ybi=1; ybi<=h_in_gen_vs_obs_a->GetNbinsX(); ybi++ ) {
         float row_sum = 0. ;
         for ( int xbi=1; xbi<=h_in_gen_vs_obs_a->GetNbinsX(); xbi++ ) {
            row_sum += h_in_gen_vs_obs_a -> GetBinContent( xbi, ybi ) ;
         } // xbi
         for ( int xbi=1; xbi<=h_in_gen_vs_obs_a->GetNbinsX(); xbi++ ) {
            if ( row_sum > 0 ) {
               h_normalized_response_a -> SetBinContent( xbi, ybi, ( h_in_gen_vs_obs_a -> GetBinContent( xbi, ybi ) ) / row_sum ) ;
               h_normalized_response_a -> SetBinError( xbi, ybi, 0. ) ;
            } else {
               h_normalized_response_a -> SetBinContent( xbi, ybi, 0. ) ;
               h_normalized_response_a -> SetBinError( xbi, ybi, 0. ) ;
            }
         } // xbi
      } // ybi



      TH2F* h_normalized_response_b = (TH2F*) h_in_gen_vs_obs_b -> Clone( "h_normalized_response_b" ) ;
      for ( int ybi=1; ybi<=h_in_gen_vs_obs_b->GetNbinsX(); ybi++ ) {
         float row_sum = 0. ;
         for ( int xbi=1; xbi<=h_in_gen_vs_obs_b->GetNbinsX(); xbi++ ) {
            row_sum += h_in_gen_vs_obs_b -> GetBinContent( xbi, ybi ) ;
         } // xbi
         for ( int xbi=1; xbi<=h_in_gen_vs_obs_b->GetNbinsX(); xbi++ ) {
            if ( row_sum > 0 ) {
               h_normalized_response_b -> SetBinContent( xbi, ybi, ( h_in_gen_vs_obs_b -> GetBinContent( xbi, ybi ) ) / row_sum ) ;
               h_normalized_response_b -> SetBinError( xbi, ybi, 0. ) ;
            } else {
               h_normalized_response_b -> SetBinContent( xbi, ybi, 0. ) ;
               h_normalized_response_b -> SetBinError( xbi, ybi, 0. ) ;
            }
         } // xbi
      } // ybi


      TH2F* h_normalized_response_diff = (TH2F*) h_normalized_response_a -> Clone( "h_normalized_response_diff" ) ;
      h_normalized_response_diff -> Add( h_normalized_response_b, -1. ) ;








      gStyle -> SetOptStat(0) ;
      gStyle -> SetPadRightMargin(0.18) ;
      gStyle -> SetPadLeftMargin(0.15) ;
      gStyle -> SetPadBottomMargin(0.15) ;
      gStyle -> SetTitleBorderSize(0) ;
      gStyle -> SetTitleY(0.975) ;
      gStyle -> SetTitleW(0.95) ;



      int can_width = 600 ;
      int can_height = 600 ;

      int can_spacing = 20 ;


     //-----

      TCanvas* can1 = (TCanvas*) gDirectory -> FindObject( "can1" ) ;
      if ( can1 == 0x0 ) { printf( "Making can1\n") ; can1 = new TCanvas( "can1", "", can_spacing, can_spacing, can_width, can_height ) ; }
      can1 -> Clear() ;
      can1 -> cd() ;

      h_normalized_response_a -> SetTitleSize( 0.045, "x" ) ;
      h_normalized_response_a -> SetTitleSize( 0.045, "y" ) ;

      h_normalized_response_a -> SetTitleOffset( 1.2, "x" ) ;
      h_normalized_response_a -> SetTitleOffset( 1.6, "y" ) ;
      sprintf( htitle, "Normalized response matrix, %s", input_name_a ) ;
      h_normalized_response_a -> SetTitle( htitle ) ;
      if ( strcmp( var_name, "y" ) == 0 ) { h_normalized_response_a -> SetNdivisions( 605 ) ; }

      h_normalized_response_a -> Draw("colz") ;
      TExec* change_hist_palette2 = new TExec( "change_hist_palette2", "Setup2DhistPalette();" );
      change_hist_palette2->Draw() ;
      h_normalized_response_a -> Draw("colz same") ;
      h_normalized_response_a -> Draw("axis same") ;

      can1 -> Update() ; can1 -> Draw() ; gSystem -> ProcessEvents() ;

      sprintf( fname, "paper-plots/dis-response-%s-%s-%s.pdf", var_name, method_name, input_name_a ) ;
      can1 -> SaveAs( fname ) ;
      sprintf( fname, "paper-plots/dis-response-%s-%s-%s.png", var_name, method_name, input_name_a ) ;
      can1 -> SaveAs( fname ) ;


     //-----

      TCanvas* can2 = (TCanvas*) gDirectory -> FindObject( "can2" ) ;
      if ( can2 == 0x0 ) { printf( "Making can2\n") ; can2 = new TCanvas( "can2", "", 2*can_spacing + can_width, can_spacing, can_width, can_height ) ; }
      can2 -> Clear() ;
      can2 -> cd() ;

      h_normalized_response_b -> SetTitleSize( 0.045, "x" ) ;
      h_normalized_response_b -> SetTitleSize( 0.045, "y" ) ;

      h_normalized_response_b -> SetTitleOffset( 1.2, "x" ) ;
      h_normalized_response_b -> SetTitleOffset( 1.6, "y" ) ;
      sprintf( htitle, "Normalized response matrix, %s", input_name_b ) ;
      h_normalized_response_b -> SetTitle( htitle ) ;
      if ( strcmp( var_name, "y" ) == 0 ) { h_normalized_response_b -> SetNdivisions( 605 ) ; }

      h_normalized_response_b -> Draw("colz") ;
      change_hist_palette2->Draw() ;
      h_normalized_response_b -> Draw("colz same") ;
      h_normalized_response_b -> Draw("axis same") ;

      can2 -> Update() ; can2 -> Draw() ; gSystem -> ProcessEvents() ;

      sprintf( fname, "paper-plots/dis-response-%s-%s-%s.pdf", var_name, method_name, input_name_b ) ;
      can2 -> SaveAs( fname ) ;
      sprintf( fname, "paper-plots/dis-response-%s-%s-%s.png", var_name, method_name, input_name_b ) ;
      can2 -> SaveAs( fname ) ;


     //-----

      TCanvas* can3 = (TCanvas*) gDirectory -> FindObject( "can3" ) ;
      if ( can3 == 0x0 ) { printf( "Making can3\n") ; can3 = new TCanvas( "can3", "", 3*can_spacing + 2*can_width , can_spacing, can_width, can_height ) ; }
      can3 -> Clear() ;
      can3 -> cd() ;

      h_normalized_response_diff -> SetMaximum( 0.02 ) ;
      h_normalized_response_diff -> SetMinimum(-0.02 ) ;

      h_normalized_response_diff -> SetTitleSize( 0.045, "x" ) ;
      h_normalized_response_diff -> SetTitleSize( 0.045, "y" ) ;

      h_normalized_response_diff -> SetTitleOffset( 1.2, "x" ) ;
      h_normalized_response_diff -> SetTitleOffset( 1.6, "y" ) ;
      sprintf( htitle, "Response matrix difference, %s, %s - %s", method_name, input_name_a, input_name_b ) ;
      h_normalized_response_diff -> SetTitle( htitle ) ;
      if ( strcmp( var_name, "y" ) == 0 ) { h_normalized_response_diff -> SetNdivisions( 605 ) ; }

      h_normalized_response_diff -> Draw("colz") ;
      TExec* change_cor_palette = new TExec( "change_cor_palette", "SetupCorrelationPalette();" );
      change_cor_palette->Draw() ;
      h_normalized_response_diff -> Draw("colz same") ;
      h_normalized_response_diff -> Draw("axis same") ;

      can3 -> Update() ; can3 -> Draw() ; gSystem -> ProcessEvents() ;

      sprintf( fname, "paper-plots/dis-response-%s-%s-diff.pdf", var_name, method_name ) ;
      can3 -> SaveAs( fname ) ;
      sprintf( fname, "paper-plots/dis-response-%s-%s-diff.png", var_name, method_name ) ;
      can3 -> SaveAs( fname ) ;











     //-----


      printf("\n\n\n") ;

      printf(" cut and paste for this:\n\n") ;

      printf("     paper_plots_dis_gen_response_comp(\"%s\",\"%s\",\"%s\",%9.5f,\"%s\",\"%s\",\"%s\",\"%s\")\n",
                      hist_name, method_name, var_name, max_diff, input_file_a, input_file_b, input_name_a, input_name_b ) ;


      printf("\n\n\n") ;

   }













