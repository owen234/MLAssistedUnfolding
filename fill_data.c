#define fill_data_cxx
#include "fill_data.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "histio.c"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldTUnfold.h"

void fill_data::Loop( const char* rur_file_name )
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();

   char hname[100] ;

   gDirectory -> Delete( "h*" ) ;

   printf("\n\n Ntuple has %lld entries\n\n", nentries ) ;


  //--- Take the RooUnfoldResponse object used in the unfolding to define the observed histogram
  //     to ensure consistency.

   TFile* tf_in = new TFile( rur_file_name, "read" ) ;
   if ( tf_in == 0x0 ) { printf("\n\n *** Bad input file for RooUnfoldResponse objects: %s\n\n", rur_file_name ) ; gSystem -> Exit(-1) ; }

   tf_in -> ls() ;

   char rur_name_a[100] ;
   RooUnfoldResponse* rur_a ;
   TH1* h_template ;

   sprintf( rur_name_a, "rur_2D_log10_q2_vs_log10_y_dnn" ) ;
   rur_a = (RooUnfoldResponse*) tf_in -> Get( rur_name_a ) ;
   if ( rur_a == 0x0 ) { printf("\n\n *** can't find %s\n\n", rur_name_a ) ; return ; }

   h_template = rur_a -> Hmeasured() ;

   printf("\n\n") ;
   printf("   Template histogram:  %s\n", h_template -> GetName() ) ;
   printf("     x axis:  nbins %2d,  [%9.5f, %9.5f]\n",
     h_template -> GetNbinsX(),
     h_template -> GetXaxis() -> GetBinLowEdge(1), 
     h_template -> GetXaxis() -> GetBinLowEdge( h_template -> GetNbinsX() + 1 ) ) ;
   printf("     y axis:  nbins %2d,  [%9.5f, %9.5f]\n",
     h_template -> GetNbinsY(),
     h_template -> GetYaxis() -> GetBinLowEdge(1), 
     h_template -> GetYaxis() -> GetBinLowEdge( h_template -> GetNbinsX() + 1 ) ) ;
   printf("\n\n") ;


   TString tsname_dnn( h_template -> GetName() ) ;
   tsname_dnn.ReplaceAll( "dummy", "data" ) ;
   tsname_dnn.Append("_dnn") ;
   TH1* h_data_dnn = (TH1*) h_template -> Clone( tsname_dnn.Data() ) ;
   h_data_dnn -> Reset() ;
   h_data_dnn -> SetTitle( tsname_dnn.Data() ) ;

   TString tsname ;

   tsname = tsname_dnn ;
   tsname.ReplaceAll( "dnn", "e" ) ;
   TH1* h_data_e = (TH1*) h_template -> Clone( tsname.Data() ) ;
   h_data_e -> Reset() ;
   h_data_e -> SetTitle( tsname.Data() ) ;

   tsname = tsname_dnn ;
   tsname.ReplaceAll( "dnn", "esigma" ) ;
   TH1* h_data_esigma = (TH1*) h_template -> Clone( tsname.Data() ) ;
   h_data_esigma -> Reset() ;
   h_data_esigma -> SetTitle( tsname.Data() ) ;

   tsname = tsname_dnn ;
   tsname.ReplaceAll( "dnn", "isigma" ) ;
   TH1* h_data_isigma = (TH1*) h_template -> Clone( tsname.Data() ) ;
   h_data_isigma -> Reset() ;
   h_data_isigma -> SetTitle( tsname.Data() ) ;

   tsname = tsname_dnn ;
   tsname.ReplaceAll( "dnn", "da" ) ;
   TH1* h_data_da = (TH1*) h_template -> Clone( tsname.Data() ) ;
   h_data_da -> Reset() ;
   h_data_da -> SetTitle( tsname.Data() ) ;




   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      int ei = jentry ;

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if ( ei%100 == 0 ) {
         printf(" --- Event: %7d / %lld    %6.3f\r", ei, nentries, (1.*ei)/(1.*nentries) ) ;
         fflush(stdout) ;
      }

      if ( dnn_Q2 < 240 ) continue ;


      if ( ei < 50 ) {
         printf("  %3d : dnn  %9.5f, %9.5f  |  e  %9.5f, %9.5f  |  esigma %9.5f, %9.5f  |  isigma %9.5f, %9.5f  |  da  %9.5f, %9.5f |\n",
           ei,
           log10( dnn_y ), log10( dnn_Q2 ),
           log10( obs_y_e ), log10( obs_Q2_e ),
           log10( obs_y_eSigma ), log10( obs_Q2_eSigma ),
           log10( obs_y_ISigma ), log10( obs_Q2_ISigma ),
           log10( obs_y_DA ), log10( obs_Q2_DA ) ) ;
      }

      h_data_dnn -> Fill( log10( dnn_y ), log10( dnn_Q2 ) ) ;
      h_data_e -> Fill( log10( obs_y_e ), log10( obs_Q2_e ) ) ;
      h_data_esigma -> Fill( log10( obs_y_eSigma ), log10( obs_Q2_eSigma ) ) ;
      h_data_isigma -> Fill( log10( obs_y_ISigma ), log10( obs_Q2_ISigma ) ) ;
      h_data_da -> Fill( log10( obs_y_DA ), log10( obs_Q2_DA ) ) ;


   } // jentry

   h_data_dnn -> DrawCopy("colz") ;

   saveHist( "data-hists.root", "*" ) ;

   loadHist( "data-hists.root" ) ;
   gDirectory -> ls() ;

   printf("\n\n\n") ;

} // Loop










