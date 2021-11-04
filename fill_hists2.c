#define fill_hists2_cxx
#include "fill_hists2.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "histio.c"


//
//
//  Before doing
//
//    .L fill_hists2.c
//
//  do
//
//   gSystem -> Load( "../RooUnfold/build/libRooUnfold.dylib" ) ;
//
//

void fill_hists2::Loop( int nbins_gen, int nbins_obs, bool verbose, int last_event, int first_event, const char* out_file )
{

   //////float wgt = 1. ; // for athena only

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();

   printf("\n\n Ntuple has %lld entries\n\n", nentries ) ;

   gDirectory -> Delete( "h*" ) ;

   float xmin = pow( 10, -2.5 ) ;
   float xmax = pow( 10, -0.3 ) ;

   float ymin = pow( 10, -2.3 ) ;
   float ymax = pow( 10, -0.2 ) ;

   float q2min = 1e2 ;
   float q2max = 1e4 ;






   TH2F* h_dummy_log10_q2_vs_log10_x_gen = new TH2F( "h_dummy_log10_q2_vs_log10_x_gen", "Gen Template for 2D unfolding, log10 Q2 vs log10 x",
        nbins_gen, log10(xmin), log10(xmax),    nbins_gen, log10(q2min), log10(q2max) ) ;

   TH2F* h_dummy_log10_q2_vs_log10_x_obs = new TH2F( "h_dummy_log10_q2_vs_log10_x_obs", "Obs Template for 2D unfolding, log10 Q2 vs log10 x",
        nbins_obs, log10(xmin), log10(xmax),    nbins_obs, log10(q2min), log10(q2max) ) ;

   RooUnfoldResponse* rur_2D_log10_q2_vs_log10_x_dnn = new RooUnfoldResponse( h_dummy_log10_q2_vs_log10_x_obs, h_dummy_log10_q2_vs_log10_x_gen,
         "rur_2D_log10_q2_vs_log10_x_dnn", "rur_2D_log10_q2_vs_log10_x_dnn" ) ;
   RooUnfoldResponse* rur_2D_log10_q2_vs_log10_x_e   = new RooUnfoldResponse( h_dummy_log10_q2_vs_log10_x_obs, h_dummy_log10_q2_vs_log10_x_gen,
         "rur_2D_log10_q2_vs_log10_x_e", "rur_2D_log10_q2_vs_log10_x_e" ) ;
   RooUnfoldResponse* rur_2D_log10_q2_vs_log10_x_isigma = new RooUnfoldResponse( h_dummy_log10_q2_vs_log10_x_obs, h_dummy_log10_q2_vs_log10_x_gen,
         "rur_2D_log10_q2_vs_log10_x_isigma", "rur_2D_log10_q2_vs_log10_x_isigma" ) ;
   RooUnfoldResponse* rur_2D_log10_q2_vs_log10_x_da = new RooUnfoldResponse( h_dummy_log10_q2_vs_log10_x_obs, h_dummy_log10_q2_vs_log10_x_gen,
         "rur_2D_log10_q2_vs_log10_x_da", "rur_2D_log10_q2_vs_log10_x_da" ) ;






   RooUnfoldResponse* rur_log10_q2_gen_vs_obs_dnn = new RooUnfoldResponse( nbins_obs, log10(q2min), log10(q2max), nbins_gen, log10(q2min), log10(q2max), "rur_log10_q2_gen_vs_obs_dnn", "rur_log10_q2_gen_vs_obs_dnn" ) ;


   TH2F* h_log10_q2_gen_vs_obs_dnn    = new TH2F( "h_log10_q2_gen_vs_obs_dnn", "h_log10_q2_gen_vs_obs_dnn", nbins_obs, log10(q2min), log10(q2max),   nbins_gen, log10(q2min), log10(q2max) ) ;
   TH2F* h_log10_q2_gen_vs_obs_e      = new TH2F( "h_log10_q2_gen_vs_obs_e", "h_log10_q2_gen_vs_obs_e", nbins_obs, log10(q2min), log10(q2max),   nbins_gen, log10(q2min), log10(q2max) ) ;
   TH2F* h_log10_q2_gen_vs_obs_isigma = new TH2F( "h_log10_q2_gen_vs_obs_isigma", "h_log10_q2_gen_vs_obs_isigma", nbins_obs, log10(q2min), log10(q2max),   nbins_gen, log10(q2min), log10(q2max) ) ;
   TH2F* h_log10_q2_gen_vs_obs_da     = new TH2F( "h_log10_q2_gen_vs_obs_da", "h_log10_q2_gen_vs_obs_da", nbins_obs, log10(q2min), log10(q2max),   nbins_gen, log10(q2min), log10(q2max) ) ;


   TH2F* h_log10_x_gen_vs_obs_dnn    = new TH2F( "h_log10_x_gen_vs_obs_dnn", "h_log10_x_gen_vs_obs_dnn", nbins_obs, log10(xmin), log10(xmax),   nbins_gen, log10(xmin), log10(xmax) ) ;
   TH2F* h_log10_x_gen_vs_obs_e      = new TH2F( "h_log10_x_gen_vs_obs_e", "h_log10_x_gen_vs_obs_e", nbins_obs, log10(xmin), log10(xmax),   nbins_gen, log10(xmin), log10(xmax) ) ;
   TH2F* h_log10_x_gen_vs_obs_isigma = new TH2F( "h_log10_x_gen_vs_obs_isigma", "h_log10_x_gen_vs_obs_isigma", nbins_obs, log10(xmin), log10(xmax),   nbins_gen, log10(xmin), log10(xmax) ) ;
   TH2F* h_log10_x_gen_vs_obs_da     = new TH2F( "h_log10_x_gen_vs_obs_da", "h_log10_x_gen_vs_obs_da", nbins_obs, log10(xmin), log10(xmax),   nbins_gen, log10(xmin), log10(xmax) ) ;

   RooUnfoldResponse* rur_log10_x_gen_vs_obs_dnn = new RooUnfoldResponse( nbins_obs, log10(xmin), log10(xmax),  nbins_gen, log10(xmin), log10(xmax),  "rur_log10_x_gen_vs_obs_dnn", "rur_log10_x_gen_vs_obs_dnn" ) ;
   RooUnfoldResponse* rur_log10_x_gen_vs_obs_e = new RooUnfoldResponse( nbins_obs, log10(xmin), log10(xmax),  nbins_gen, log10(xmin), log10(xmax),  "rur_log10_x_gen_vs_obs_e", "rur_log10_x_gen_vs_obs_e" ) ;
   RooUnfoldResponse* rur_log10_x_gen_vs_obs_isigma = new RooUnfoldResponse( nbins_obs, log10(xmin), log10(xmax),  nbins_gen, log10(xmin), log10(xmax),  "rur_log10_x_gen_vs_obs_isigma", "rur_log10_x_gen_vs_obs_isigma" ) ;
   RooUnfoldResponse* rur_log10_x_gen_vs_obs_da = new RooUnfoldResponse( nbins_obs, log10(xmin), log10(xmax),  nbins_gen, log10(xmin), log10(xmax),  "rur_log10_x_gen_vs_obs_da", "rur_log10_x_gen_vs_obs_da" ) ;



   TH2F* h_log10_y_gen_vs_obs_dnn    = new TH2F( "h_log10_y_gen_vs_obs_dnn", "h_log10_y_gen_vs_obs_dnn", nbins_obs, log10(ymin), log10(ymax),   nbins_gen, log10(ymin), log10(ymax) ) ;
   TH2F* h_log10_y_gen_vs_obs_e      = new TH2F( "h_log10_y_gen_vs_obs_e", "h_log10_y_gen_vs_obs_e", nbins_obs, log10(ymin), log10(ymax),   nbins_gen, log10(ymin), log10(ymax) ) ;
   TH2F* h_log10_y_gen_vs_obs_isigma = new TH2F( "h_log10_y_gen_vs_obs_isigma", "h_log10_y_gen_vs_obs_isigma", nbins_obs, log10(ymin), log10(ymax),   nbins_gen, log10(ymin), log10(ymax) ) ;
   TH2F* h_log10_y_gen_vs_obs_da     = new TH2F( "h_log10_y_gen_vs_obs_da", "h_log10_y_gen_vs_obs_da", nbins_obs, log10(ymin), log10(ymax),   nbins_gen, log10(ymin), log10(ymax) ) ;

   RooUnfoldResponse* rur_log10_y_gen_vs_obs_dnn = new RooUnfoldResponse( nbins_obs, log10(ymin), log10(ymax),  nbins_gen, log10(ymin), log10(ymax),  "rur_log10_y_gen_vs_obs_dnn", "rur_log10_y_gen_vs_obs_dnn" ) ;
   RooUnfoldResponse* rur_log10_y_gen_vs_obs_e = new RooUnfoldResponse( nbins_obs, log10(ymin), log10(ymax),  nbins_gen, log10(ymin), log10(ymax),  "rur_log10_y_gen_vs_obs_e", "rur_log10_y_gen_vs_obs_e" ) ;
   RooUnfoldResponse* rur_log10_y_gen_vs_obs_isigma = new RooUnfoldResponse( nbins_obs, log10(ymin), log10(ymax),  nbins_gen, log10(ymin), log10(ymax),  "rur_log10_y_gen_vs_obs_isigma", "rur_log10_y_gen_vs_obs_isigma" ) ;
   RooUnfoldResponse* rur_log10_y_gen_vs_obs_da = new RooUnfoldResponse( nbins_obs, log10(ymin), log10(ymax),  nbins_gen, log10(ymin), log10(ymax),  "rur_log10_y_gen_vs_obs_da", "rur_log10_y_gen_vs_obs_da" ) ;



   TH2F* h_x_gen_vs_obs_dnn    = new TH2F( "h_x_gen_vs_obs_dnn", "h_x_gen_vs_obs_dnn", nbins_obs, (xmin), (xmax),   nbins_gen, (xmin), (xmax) ) ;
   TH2F* h_x_gen_vs_obs_e      = new TH2F( "h_x_gen_vs_obs_e", "h_x_gen_vs_obs_e", nbins_obs, (xmin), (xmax),   nbins_gen, (xmin), (xmax) ) ;
   TH2F* h_x_gen_vs_obs_isigma = new TH2F( "h_x_gen_vs_obs_isigma", "h_x_gen_vs_obs_isigma", nbins_obs, (xmin), (xmax),   nbins_gen, (xmin), (xmax) ) ;
   TH2F* h_x_gen_vs_obs_da     = new TH2F( "h_x_gen_vs_obs_da", "h_x_gen_vs_obs_da", nbins_obs, (xmin), (xmax),   nbins_gen, (xmin), (xmax) ) ;


   h_log10_q2_gen_vs_obs_dnn    -> SetXTitle( "Obs log10(Q2), DNN" ) ;
   h_log10_q2_gen_vs_obs_e      -> SetXTitle( "Obs log10(Q2), e" ) ;
   h_log10_q2_gen_vs_obs_isigma -> SetXTitle( "Obs log10(Q2), ISigma" ) ;
   h_log10_q2_gen_vs_obs_da     -> SetXTitle( "Obs log10(Q2), DA" ) ;

   h_log10_q2_gen_vs_obs_dnn    -> SetYTitle( "Gen log10(Q2), DNN" ) ;
   h_log10_q2_gen_vs_obs_e      -> SetYTitle( "Gen log10(Q2), e" ) ;
   h_log10_q2_gen_vs_obs_isigma -> SetYTitle( "Gen log10(Q2), ISigma" ) ;
   h_log10_q2_gen_vs_obs_da     -> SetYTitle( "Gen log10(Q2), DA" ) ;


   h_log10_x_gen_vs_obs_dnn    -> SetXTitle( "Obs log10(x), DNN" ) ;
   h_log10_x_gen_vs_obs_e      -> SetXTitle( "Obs log10(x), e" ) ;
   h_log10_x_gen_vs_obs_isigma -> SetXTitle( "Obs log10(x), ISigma" ) ;
   h_log10_x_gen_vs_obs_da     -> SetXTitle( "Obs log10(x), DA" ) ;

   h_log10_x_gen_vs_obs_dnn    -> SetYTitle( "Gen log10(x), DNN" ) ;
   h_log10_x_gen_vs_obs_e      -> SetYTitle( "Gen log10(x), e" ) ;
   h_log10_x_gen_vs_obs_isigma -> SetYTitle( "Gen log10(x), ISigma" ) ;
   h_log10_x_gen_vs_obs_da     -> SetYTitle( "Gen log10(x), DA" ) ;


   h_log10_y_gen_vs_obs_dnn    -> SetXTitle( "Obs log10(y), DNN" ) ;
   h_log10_y_gen_vs_obs_e      -> SetXTitle( "Obs log10(y), e" ) ;
   h_log10_y_gen_vs_obs_isigma -> SetXTitle( "Obs log10(y), ISigma" ) ;
   h_log10_y_gen_vs_obs_da     -> SetXTitle( "Obs log10(y), DA" ) ;

   h_log10_y_gen_vs_obs_dnn    -> SetYTitle( "Gen log10(y), DNN" ) ;
   h_log10_y_gen_vs_obs_e      -> SetYTitle( "Gen log10(y), e" ) ;
   h_log10_y_gen_vs_obs_isigma -> SetYTitle( "Gen log10(y), ISigma" ) ;
   h_log10_y_gen_vs_obs_da     -> SetYTitle( "Gen log10(y), DA" ) ;


   h_x_gen_vs_obs_dnn    -> SetXTitle( "Obs x, DNN" ) ;
   h_x_gen_vs_obs_e      -> SetXTitle( "Obs x, e" ) ;
   h_x_gen_vs_obs_isigma -> SetXTitle( "Obs x, ISigma" ) ;
   h_x_gen_vs_obs_da     -> SetXTitle( "Obs x, DA" ) ;

   h_x_gen_vs_obs_dnn    -> SetYTitle( "Gen x, DNN" ) ;
   h_x_gen_vs_obs_e      -> SetYTitle( "Gen x, e" ) ;
   h_x_gen_vs_obs_isigma -> SetYTitle( "Gen x, ISigma" ) ;
   h_x_gen_vs_obs_da     -> SetYTitle( "Gen x, DA" ) ;



   if ( last_event > 0 ) {
      nentries = last_event ;
   }

   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=first_event; jentry<nentries;jentry++) {

      int ei = jentry ;

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if ( !verbose && ei%100 == 0 ) {
         printf(" --- Event: %7d / %lld    %6.3f\r", ei, nentries, (1.*ei)/(1.*nentries) ) ;
         fflush(stdout) ;
      }


      if ( from_tlv_gen_Q2 < 220 ) continue ;
////  if ( Empz < 45 ) continue ; // unnecessary.  already in DNN ntuple maker
////  if ( Empz > 65 ) continue ;

      float log10_gen_Q2 = log10(from_tlv_gen_Q2) ;
      float log10_gen_x = log10(from_tlv_gen_x) ;
      float log10_gen_y = log10(from_tlv_gen_y) ;

      h_log10_q2_gen_vs_obs_dnn -> Fill( log10(dnn_Q2), log10_gen_Q2, wgt ) ;
      h_log10_q2_gen_vs_obs_e -> Fill( log10(obs_Q2_e), log10_gen_Q2, wgt ) ;
      h_log10_q2_gen_vs_obs_isigma -> Fill( log10(obs_Q2_ISigma), log10_gen_Q2, wgt ) ;
      h_log10_q2_gen_vs_obs_da -> Fill( log10(obs_Q2_DA), log10_gen_Q2, wgt ) ;

      h_log10_x_gen_vs_obs_dnn -> Fill( log10(dnn_x), log10_gen_x, wgt ) ;
      h_log10_x_gen_vs_obs_e -> Fill( log10(obs_x_e), log10_gen_x, wgt ) ;
      h_log10_x_gen_vs_obs_isigma -> Fill( log10(obs_x_ISigma), log10_gen_x, wgt ) ;
      h_log10_x_gen_vs_obs_da -> Fill( log10(obs_x_DA), log10_gen_x, wgt ) ;

      h_log10_y_gen_vs_obs_dnn -> Fill( log10(dnn_y), log10_gen_y, wgt ) ;
      h_log10_y_gen_vs_obs_e -> Fill( log10(obs_y_e), log10_gen_y, wgt ) ;
      h_log10_y_gen_vs_obs_isigma -> Fill( log10(obs_y_ISigma), log10_gen_y, wgt ) ;
      h_log10_y_gen_vs_obs_da -> Fill( log10(obs_y_DA), log10_gen_y, wgt ) ;

      h_x_gen_vs_obs_dnn -> Fill( (dnn_x), from_tlv_gen_x, wgt ) ;
      h_x_gen_vs_obs_e -> Fill( (obs_x_e), from_tlv_gen_x, wgt ) ;
      h_x_gen_vs_obs_isigma -> Fill( (obs_x_ISigma), from_tlv_gen_x, wgt ) ;
      h_x_gen_vs_obs_da -> Fill( (obs_x_DA), from_tlv_gen_x, wgt ) ;

      rur_log10_q2_gen_vs_obs_dnn -> Fill( log10(dnn_Q2), log10_gen_Q2, wgt ) ;

      rur_log10_x_gen_vs_obs_dnn -> Fill( log10(dnn_x), log10_gen_x, wgt ) ;
      rur_log10_x_gen_vs_obs_e -> Fill( log10(obs_x_e), log10_gen_x, wgt ) ;
      rur_log10_x_gen_vs_obs_isigma -> Fill( log10(obs_x_ISigma), log10_gen_x, wgt ) ;
      rur_log10_x_gen_vs_obs_da -> Fill( log10(obs_x_DA), log10_gen_x, wgt ) ;

      rur_log10_y_gen_vs_obs_dnn -> Fill( log10(dnn_y), log10_gen_y, wgt ) ;
      rur_log10_y_gen_vs_obs_e -> Fill( log10(obs_y_e), log10_gen_y, wgt ) ;
      rur_log10_y_gen_vs_obs_isigma -> Fill( log10(obs_y_ISigma), log10_gen_y, wgt ) ;
      rur_log10_y_gen_vs_obs_da -> Fill( log10(obs_y_DA), log10_gen_y, wgt ) ;

      rur_2D_log10_q2_vs_log10_x_dnn    -> Fill( log10(dnn_x), log10(dnn_Q2),   log10_gen_x, log10_gen_Q2,  wgt ) ;
      rur_2D_log10_q2_vs_log10_x_e      -> Fill( log10(obs_x_e), log10(obs_Q2_e),   log10_gen_x, log10_gen_Q2,  wgt ) ;
      rur_2D_log10_q2_vs_log10_x_isigma -> Fill( log10(obs_x_ISigma), log10(obs_Q2_ISigma),   log10_gen_x, log10_gen_Q2,  wgt ) ;
      rur_2D_log10_q2_vs_log10_x_da     -> Fill( log10(obs_x_DA), log10(obs_Q2_DA),   log10_gen_x, log10_gen_Q2,  wgt ) ;


   } // jentry

   printf("\n\n Done. \n\n") ;

   gDirectory -> ls() ;

   char save_file[1000] ;
   if ( strlen( out_file ) > 0 ) {
      sprintf( save_file, "%s", out_file ) ;
   } else {
      sprintf( save_file, "unfold-hists-input-nbins_gen%03d_obs%03d.root", nbins_gen, nbins_obs ) ;
   }

   saveHist( save_file, "*" ) ;

   TFile tf( save_file, "update" ) ;

   tf.WriteTObject( rur_log10_q2_gen_vs_obs_dnn, rur_log10_q2_gen_vs_obs_dnn->GetName() ) ;

   tf.WriteTObject( rur_log10_x_gen_vs_obs_dnn, rur_log10_x_gen_vs_obs_dnn->GetName() ) ;
   tf.WriteTObject( rur_log10_x_gen_vs_obs_e, rur_log10_x_gen_vs_obs_e->GetName() ) ;
   tf.WriteTObject( rur_log10_x_gen_vs_obs_isigma, rur_log10_x_gen_vs_obs_isigma->GetName() ) ;
   tf.WriteTObject( rur_log10_x_gen_vs_obs_da, rur_log10_x_gen_vs_obs_da->GetName() ) ;

   tf.WriteTObject( rur_log10_y_gen_vs_obs_dnn, rur_log10_y_gen_vs_obs_dnn->GetName() ) ;
   tf.WriteTObject( rur_log10_y_gen_vs_obs_e, rur_log10_y_gen_vs_obs_e->GetName() ) ;
   tf.WriteTObject( rur_log10_y_gen_vs_obs_isigma, rur_log10_y_gen_vs_obs_isigma->GetName() ) ;
   tf.WriteTObject( rur_log10_y_gen_vs_obs_da, rur_log10_y_gen_vs_obs_da->GetName() ) ;

   tf.WriteTObject( rur_2D_log10_q2_vs_log10_x_dnn, rur_2D_log10_q2_vs_log10_x_dnn->GetName() ) ;
   tf.WriteTObject( rur_2D_log10_q2_vs_log10_x_e, rur_2D_log10_q2_vs_log10_x_e->GetName() ) ;
   tf.WriteTObject( rur_2D_log10_q2_vs_log10_x_isigma, rur_2D_log10_q2_vs_log10_x_isigma->GetName() ) ;
   tf.WriteTObject( rur_2D_log10_q2_vs_log10_x_da, rur_2D_log10_q2_vs_log10_x_da->GetName() ) ;


   tf.Close() ;

} // Loop








