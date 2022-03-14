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

   gDirectory -> Delete( "h*" ) ;

   printf("\n\n Ntuple has %lld entries\n\n", nentries ) ;


  //-- make edge bins contain a fraction grf of a bin width of the good range
  //   That is, within the highest bin, grf * bw of that bin is within the good range and (1-grf)*bw is above the good range
  //   Do this for the gen binning.
  //
  //   This seems to help the unfolding a lot.  It avoids bias in the edge bins and heavy correlations.
  //

   //float grf = 0.2 ;
   float grf = 0.05 ;
   //float grf = 1.0 ; // no edge buffer bin



   float good_range_log10_xmin = -2.5  ;
   float good_range_log10_xmax = 0.0 ;



   //////////float good_range_log10_ymin = -2.3  ;
   //////////float good_range_log10_ymax = -0.2  ;

   ////////float good_range_log10_q2min = 2.342 ; // log10(220) = 2.342
   ////////float good_range_log10_q2max = 4 ;




   float good_range_log10_ymin = -2.0  ;
   float good_range_log10_ymax = 0.0  ;

   float good_range_log10_q2min = 2.3979400 ; // log10(250) = 2.3979400
   float good_range_log10_q2max = 4 ;





   float log10_xmin = 0.5 * ( good_range_log10_xmax + good_range_log10_xmin - ( good_range_log10_xmax - good_range_log10_xmin ) / ( 1. - 2.*(1.-grf)/(1.*nbins_gen) ) ) ;
   float log10_xmax = 0.5 * ( good_range_log10_xmax + good_range_log10_xmin + ( good_range_log10_xmax - good_range_log10_xmin ) / ( 1. - 2.*(1.-grf)/(1.*nbins_gen) ) ) ;

   float log10_ymin = 0.5 * ( good_range_log10_ymax + good_range_log10_ymin - ( good_range_log10_ymax - good_range_log10_ymin ) / ( 1. - 2.*(1.-grf)/(1.*nbins_gen) ) ) ;
   float log10_ymax = 0.5 * ( good_range_log10_ymax + good_range_log10_ymin + ( good_range_log10_ymax - good_range_log10_ymin ) / ( 1. - 2.*(1.-grf)/(1.*nbins_gen) ) ) ;

   float log10_q2min = 0.5 * ( good_range_log10_q2max + good_range_log10_q2min - ( good_range_log10_q2max - good_range_log10_q2min ) / ( 1. - 2.*(1.-grf)/(1.*nbins_gen) ) ) ;
   float log10_q2max = 0.5 * ( good_range_log10_q2max + good_range_log10_q2min + ( good_range_log10_q2max - good_range_log10_q2min ) / ( 1. - 2.*(1.-grf)/(1.*nbins_gen) ) ) ;

   printf("\n\n") ;

   printf("   log10  x :  good range  ( %9.3f , %9.3f ) ,   full range ( %9.3f , %9.5f )\n",
        good_range_log10_xmin  , good_range_log10_xmax  , log10_xmin, log10_xmax ) ;

   printf("   log10  y :  good range  ( %9.3f , %9.3f ) ,   full range ( %9.3f , %9.5f )\n",
        good_range_log10_ymin  , good_range_log10_ymax  , log10_ymin, log10_ymax ) ;

   printf("   log10 q2 :  good range  ( %9.3f , %9.3f ) ,   full range ( %9.3f , %9.5f )\n",
        good_range_log10_q2min  , good_range_log10_q2max  , log10_q2min, log10_q2max ) ;

   printf("\n\n") ;






   TH2F* h_dummy_log10_q2_vs_log10_x_gen = new TH2F( "h_dummy_log10_q2_vs_log10_x_gen", "Gen Template for 2D unfolding, log10 Q2 vs log10 x",
        nbins_gen, log10_xmin, log10_xmax,    nbins_gen, log10_q2min, log10_q2max ) ;

   TH2F* h_dummy_log10_q2_vs_log10_x_obs = new TH2F( "h_dummy_log10_q2_vs_log10_x_obs", "Obs Template for 2D unfolding, log10 Q2 vs log10 x",
        nbins_obs, log10_xmin, log10_xmax,    nbins_obs, log10_q2min, log10_q2max ) ;



   TH2F* h_dummy_log10_q2_vs_log10_y_gen = new TH2F( "h_dummy_log10_q2_vs_log10_y_gen", "Gen Template for 2D unfolding, log10 Q2 vs log10 y",
        nbins_gen, log10_ymin, log10_ymax,    nbins_gen, log10_q2min, log10_q2max ) ;

   TH2F* h_dummy_log10_q2_vs_log10_y_obs = new TH2F( "h_dummy_log10_q2_vs_log10_y_obs", "Obs Template for 2D unfolding, log10 Q2 vs log10 y",
        nbins_obs, log10_ymin, log10_ymax,    nbins_obs, log10_q2min, log10_q2max ) ;





   RooUnfoldResponse* rur_2D_log10_q2_vs_log10_x_dnn = new RooUnfoldResponse( h_dummy_log10_q2_vs_log10_x_obs, h_dummy_log10_q2_vs_log10_x_gen,
         "rur_2D_log10_q2_vs_log10_x_dnn", "rur_2D_log10_q2_vs_log10_x_dnn" ) ;
   RooUnfoldResponse* rur_2D_log10_q2_vs_log10_x_e   = new RooUnfoldResponse( h_dummy_log10_q2_vs_log10_x_obs, h_dummy_log10_q2_vs_log10_x_gen,
         "rur_2D_log10_q2_vs_log10_x_e", "rur_2D_log10_q2_vs_log10_x_e" ) ;
   RooUnfoldResponse* rur_2D_log10_q2_vs_log10_x_isigma = new RooUnfoldResponse( h_dummy_log10_q2_vs_log10_x_obs, h_dummy_log10_q2_vs_log10_x_gen,
         "rur_2D_log10_q2_vs_log10_x_isigma", "rur_2D_log10_q2_vs_log10_x_isigma" ) ;
   RooUnfoldResponse* rur_2D_log10_q2_vs_log10_x_esigma = new RooUnfoldResponse( h_dummy_log10_q2_vs_log10_x_obs, h_dummy_log10_q2_vs_log10_x_gen,
         "rur_2D_log10_q2_vs_log10_x_esigma", "rur_2D_log10_q2_vs_log10_x_esigma" ) ;
   RooUnfoldResponse* rur_2D_log10_q2_vs_log10_x_da = new RooUnfoldResponse( h_dummy_log10_q2_vs_log10_x_obs, h_dummy_log10_q2_vs_log10_x_gen,
         "rur_2D_log10_q2_vs_log10_x_da", "rur_2D_log10_q2_vs_log10_x_da" ) ;



   RooUnfoldResponse* rur_2D_log10_q2_vs_log10_y_dnn = new RooUnfoldResponse( h_dummy_log10_q2_vs_log10_y_obs, h_dummy_log10_q2_vs_log10_y_gen,
         "rur_2D_log10_q2_vs_log10_y_dnn", "rur_2D_log10_q2_vs_log10_y_dnn" ) ;
   RooUnfoldResponse* rur_2D_log10_q2_vs_log10_y_e   = new RooUnfoldResponse( h_dummy_log10_q2_vs_log10_y_obs, h_dummy_log10_q2_vs_log10_y_gen,
         "rur_2D_log10_q2_vs_log10_y_e", "rur_2D_log10_q2_vs_log10_y_e" ) ;
   RooUnfoldResponse* rur_2D_log10_q2_vs_log10_y_isigma = new RooUnfoldResponse( h_dummy_log10_q2_vs_log10_y_obs, h_dummy_log10_q2_vs_log10_y_gen,
         "rur_2D_log10_q2_vs_log10_y_isigma", "rur_2D_log10_q2_vs_log10_y_isigma" ) ;
   RooUnfoldResponse* rur_2D_log10_q2_vs_log10_y_esigma = new RooUnfoldResponse( h_dummy_log10_q2_vs_log10_y_obs, h_dummy_log10_q2_vs_log10_y_gen,
         "rur_2D_log10_q2_vs_log10_y_esigma", "rur_2D_log10_q2_vs_log10_y_esigma" ) ;
   RooUnfoldResponse* rur_2D_log10_q2_vs_log10_y_da = new RooUnfoldResponse( h_dummy_log10_q2_vs_log10_y_obs, h_dummy_log10_q2_vs_log10_y_gen,
         "rur_2D_log10_q2_vs_log10_y_da", "rur_2D_log10_q2_vs_log10_y_da" ) ;







   RooUnfoldResponse* rur_log10_q2_gen_vs_obs_dnn = new RooUnfoldResponse( nbins_obs, log10_q2min, log10_q2max, nbins_gen, log10_q2min, log10_q2max, "rur_log10_q2_gen_vs_obs_dnn", "rur_log10_q2_gen_vs_obs_dnn" ) ;


   TH2F* h_log10_q2_gen_vs_obs_dnn    = new TH2F( "h_log10_q2_gen_vs_obs_dnn", "h_log10_q2_gen_vs_obs_dnn", nbins_obs, log10_q2min, log10_q2max,   nbins_gen, log10_q2min, log10_q2max ) ;
   TH2F* h_log10_q2_gen_vs_obs_e      = new TH2F( "h_log10_q2_gen_vs_obs_e", "h_log10_q2_gen_vs_obs_e", nbins_obs, log10_q2min, log10_q2max,   nbins_gen, log10_q2min, log10_q2max ) ;
   TH2F* h_log10_q2_gen_vs_obs_isigma = new TH2F( "h_log10_q2_gen_vs_obs_isigma", "h_log10_q2_gen_vs_obs_isigma", nbins_obs, log10_q2min, log10_q2max,   nbins_gen, log10_q2min, log10_q2max ) ;
   TH2F* h_log10_q2_gen_vs_obs_esigma = new TH2F( "h_log10_q2_gen_vs_obs_esigma", "h_log10_q2_gen_vs_obs_esigma", nbins_obs, log10_q2min, log10_q2max,   nbins_gen, log10_q2min, log10_q2max ) ;
   TH2F* h_log10_q2_gen_vs_obs_da     = new TH2F( "h_log10_q2_gen_vs_obs_da", "h_log10_q2_gen_vs_obs_da", nbins_obs, log10_q2min, log10_q2max,   nbins_gen, log10_q2min, log10_q2max ) ;


   TH2F* h_log10_x_gen_vs_obs_dnn    = new TH2F( "h_log10_x_gen_vs_obs_dnn", "h_log10_x_gen_vs_obs_dnn", nbins_obs, log10_xmin, log10_xmax,   nbins_gen, log10_xmin, log10_xmax ) ;
   TH2F* h_log10_x_gen_vs_obs_e      = new TH2F( "h_log10_x_gen_vs_obs_e", "h_log10_x_gen_vs_obs_e", nbins_obs, log10_xmin, log10_xmax,   nbins_gen, log10_xmin, log10_xmax ) ;
   TH2F* h_log10_x_gen_vs_obs_isigma = new TH2F( "h_log10_x_gen_vs_obs_isigma", "h_log10_x_gen_vs_obs_isigma", nbins_obs, log10_xmin, log10_xmax,   nbins_gen, log10_xmin, log10_xmax ) ;
   TH2F* h_log10_x_gen_vs_obs_esigma = new TH2F( "h_log10_x_gen_vs_obs_esigma", "h_log10_x_gen_vs_obs_esigma", nbins_obs, log10_xmin, log10_xmax,   nbins_gen, log10_xmin, log10_xmax ) ;
   TH2F* h_log10_x_gen_vs_obs_da     = new TH2F( "h_log10_x_gen_vs_obs_da", "h_log10_x_gen_vs_obs_da", nbins_obs, log10_xmin, log10_xmax,   nbins_gen, log10_xmin, log10_xmax ) ;

   RooUnfoldResponse* rur_log10_x_gen_vs_obs_dnn = new RooUnfoldResponse( nbins_obs, log10_xmin, log10_xmax,  nbins_gen, log10_xmin, log10_xmax,  "rur_log10_x_gen_vs_obs_dnn", "rur_log10_x_gen_vs_obs_dnn" ) ;
   RooUnfoldResponse* rur_log10_x_gen_vs_obs_e = new RooUnfoldResponse( nbins_obs, log10_xmin, log10_xmax,  nbins_gen, log10_xmin, log10_xmax,  "rur_log10_x_gen_vs_obs_e", "rur_log10_x_gen_vs_obs_e" ) ;
   RooUnfoldResponse* rur_log10_x_gen_vs_obs_isigma = new RooUnfoldResponse( nbins_obs, log10_xmin, log10_xmax,  nbins_gen, log10_xmin, log10_xmax,  "rur_log10_x_gen_vs_obs_isigma", "rur_log10_x_gen_vs_obs_isigma" ) ;
   RooUnfoldResponse* rur_log10_x_gen_vs_obs_esigma = new RooUnfoldResponse( nbins_obs, log10_xmin, log10_xmax,  nbins_gen, log10_xmin, log10_xmax,  "rur_log10_x_gen_vs_obs_esigma", "rur_log10_x_gen_vs_obs_esigma" ) ;
   RooUnfoldResponse* rur_log10_x_gen_vs_obs_da = new RooUnfoldResponse( nbins_obs, log10_xmin, log10_xmax,  nbins_gen, log10_xmin, log10_xmax,  "rur_log10_x_gen_vs_obs_da", "rur_log10_x_gen_vs_obs_da" ) ;



   TH2F* h_log10_y_gen_vs_obs_dnn    = new TH2F( "h_log10_y_gen_vs_obs_dnn", "h_log10_y_gen_vs_obs_dnn", nbins_obs, log10_ymin, log10_ymax,   nbins_gen, log10_ymin, log10_ymax ) ;
   TH2F* h_log10_y_gen_vs_obs_e      = new TH2F( "h_log10_y_gen_vs_obs_e", "h_log10_y_gen_vs_obs_e", nbins_obs, log10_ymin, log10_ymax,   nbins_gen, log10_ymin, log10_ymax ) ;
   TH2F* h_log10_y_gen_vs_obs_isigma = new TH2F( "h_log10_y_gen_vs_obs_isigma", "h_log10_y_gen_vs_obs_isigma", nbins_obs, log10_ymin, log10_ymax,   nbins_gen, log10_ymin, log10_ymax ) ;
   TH2F* h_log10_y_gen_vs_obs_esigma = new TH2F( "h_log10_y_gen_vs_obs_esigma", "h_log10_y_gen_vs_obs_esigma", nbins_obs, log10_ymin, log10_ymax,   nbins_gen, log10_ymin, log10_ymax ) ;
   TH2F* h_log10_y_gen_vs_obs_da     = new TH2F( "h_log10_y_gen_vs_obs_da", "h_log10_y_gen_vs_obs_da", nbins_obs, log10_ymin, log10_ymax,   nbins_gen, log10_ymin, log10_ymax ) ;

   RooUnfoldResponse* rur_log10_y_gen_vs_obs_dnn = new RooUnfoldResponse( nbins_obs, log10_ymin, log10_ymax,  nbins_gen, log10_ymin, log10_ymax,  "rur_log10_y_gen_vs_obs_dnn", "rur_log10_y_gen_vs_obs_dnn" ) ;
   RooUnfoldResponse* rur_log10_y_gen_vs_obs_e = new RooUnfoldResponse( nbins_obs, log10_ymin, log10_ymax,  nbins_gen, log10_ymin, log10_ymax,  "rur_log10_y_gen_vs_obs_e", "rur_log10_y_gen_vs_obs_e" ) ;
   RooUnfoldResponse* rur_log10_y_gen_vs_obs_isigma = new RooUnfoldResponse( nbins_obs, log10_ymin, log10_ymax,  nbins_gen, log10_ymin, log10_ymax,  "rur_log10_y_gen_vs_obs_isigma", "rur_log10_y_gen_vs_obs_isigma" ) ;
   RooUnfoldResponse* rur_log10_y_gen_vs_obs_esigma = new RooUnfoldResponse( nbins_obs, log10_ymin, log10_ymax,  nbins_gen, log10_ymin, log10_ymax,  "rur_log10_y_gen_vs_obs_esigma", "rur_log10_y_gen_vs_obs_esigma" ) ;
   RooUnfoldResponse* rur_log10_y_gen_vs_obs_da = new RooUnfoldResponse( nbins_obs, log10_ymin, log10_ymax,  nbins_gen, log10_ymin, log10_ymax,  "rur_log10_y_gen_vs_obs_da", "rur_log10_y_gen_vs_obs_da" ) ;


   float xmin = pow( 10, log10_xmin ) ;
   float xmax = pow( 10, log10_xmax ) ;

   TH2F* h_x_gen_vs_obs_dnn    = new TH2F( "h_x_gen_vs_obs_dnn", "h_x_gen_vs_obs_dnn", nbins_obs, (xmin), (xmax),   nbins_gen, (xmin), (xmax) ) ;
   TH2F* h_x_gen_vs_obs_e      = new TH2F( "h_x_gen_vs_obs_e", "h_x_gen_vs_obs_e", nbins_obs, (xmin), (xmax),   nbins_gen, (xmin), (xmax) ) ;
   TH2F* h_x_gen_vs_obs_isigma = new TH2F( "h_x_gen_vs_obs_isigma", "h_x_gen_vs_obs_isigma", nbins_obs, (xmin), (xmax),   nbins_gen, (xmin), (xmax) ) ;
   TH2F* h_x_gen_vs_obs_esigma = new TH2F( "h_x_gen_vs_obs_esigma", "h_x_gen_vs_obs_esigma", nbins_obs, (xmin), (xmax),   nbins_gen, (xmin), (xmax) ) ;
   TH2F* h_x_gen_vs_obs_da     = new TH2F( "h_x_gen_vs_obs_da", "h_x_gen_vs_obs_da", nbins_obs, (xmin), (xmax),   nbins_gen, (xmin), (xmax) ) ;


   h_log10_q2_gen_vs_obs_dnn    -> SetXTitle( "Obs log10(Q2), DNN" ) ;
   h_log10_q2_gen_vs_obs_e      -> SetXTitle( "Obs log10(Q2), e" ) ;
   h_log10_q2_gen_vs_obs_isigma -> SetXTitle( "Obs log10(Q2), ISigma" ) ;
   h_log10_q2_gen_vs_obs_isigma -> SetXTitle( "Obs log10(Q2), eSigma" ) ;
   h_log10_q2_gen_vs_obs_da     -> SetXTitle( "Obs log10(Q2), DA" ) ;

   h_log10_q2_gen_vs_obs_dnn    -> SetYTitle( "Gen log10(Q2)" ) ;
   h_log10_q2_gen_vs_obs_e      -> SetYTitle( "Gen log10(Q2)" ) ;
   h_log10_q2_gen_vs_obs_isigma -> SetYTitle( "Gen log10(Q2)" ) ;
   h_log10_q2_gen_vs_obs_esigma -> SetYTitle( "Gen log10(Q2)" ) ;
   h_log10_q2_gen_vs_obs_da     -> SetYTitle( "Gen log10(Q2)" ) ;


   h_log10_x_gen_vs_obs_dnn    -> SetXTitle( "Obs log10(x), DNN" ) ;
   h_log10_x_gen_vs_obs_e      -> SetXTitle( "Obs log10(x), e" ) ;
   h_log10_x_gen_vs_obs_isigma -> SetXTitle( "Obs log10(x), ISigma" ) ;
   h_log10_x_gen_vs_obs_esigma -> SetXTitle( "Obs log10(x), eSigma" ) ;
   h_log10_x_gen_vs_obs_da     -> SetXTitle( "Obs log10(x), DA" ) ;

   h_log10_x_gen_vs_obs_dnn    -> SetYTitle( "Gen log10(x)" ) ;
   h_log10_x_gen_vs_obs_e      -> SetYTitle( "Gen log10(x)" ) ;
   h_log10_x_gen_vs_obs_isigma -> SetYTitle( "Gen log10(x)" ) ;
   h_log10_x_gen_vs_obs_esigma -> SetYTitle( "Gen log10(x)" ) ;
   h_log10_x_gen_vs_obs_da     -> SetYTitle( "Gen log10(x)" ) ;


   h_log10_y_gen_vs_obs_dnn    -> SetXTitle( "Obs log10(y), DNN" ) ;
   h_log10_y_gen_vs_obs_e      -> SetXTitle( "Obs log10(y), e" ) ;
   h_log10_y_gen_vs_obs_isigma -> SetXTitle( "Obs log10(y), ISigma" ) ;
   h_log10_y_gen_vs_obs_esigma -> SetXTitle( "Obs log10(y), eSigma" ) ;
   h_log10_y_gen_vs_obs_da     -> SetXTitle( "Obs log10(y), DA" ) ;

   h_log10_y_gen_vs_obs_dnn    -> SetYTitle( "Gen log10(y)" ) ;
   h_log10_y_gen_vs_obs_e      -> SetYTitle( "Gen log10(y)" ) ;
   h_log10_y_gen_vs_obs_isigma -> SetYTitle( "Gen log10(y)" ) ;
   h_log10_y_gen_vs_obs_esigma -> SetYTitle( "Gen log10(y)" ) ;
   h_log10_y_gen_vs_obs_da     -> SetYTitle( "Gen log10(y)" ) ;


   h_x_gen_vs_obs_dnn    -> SetXTitle( "Obs x, DNN" ) ;
   h_x_gen_vs_obs_e      -> SetXTitle( "Obs x, e" ) ;
   h_x_gen_vs_obs_isigma -> SetXTitle( "Obs x, ISigma" ) ;
   h_x_gen_vs_obs_esigma -> SetXTitle( "Obs x, eSigma" ) ;
   h_x_gen_vs_obs_da     -> SetXTitle( "Obs x, DA" ) ;

   h_x_gen_vs_obs_dnn    -> SetYTitle( "Gen x" ) ;
   h_x_gen_vs_obs_e      -> SetYTitle( "Gen x" ) ;
   h_x_gen_vs_obs_isigma -> SetYTitle( "Gen x" ) ;
   h_x_gen_vs_obs_esigma -> SetYTitle( "Gen x" ) ;
   h_x_gen_vs_obs_da     -> SetYTitle( "Gen x" ) ;



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


      if ( log10_gen_Q2 < good_range_log10_q2min ) continue ;
      if ( log10_gen_Q2 > good_range_log10_q2max ) continue ;

      if ( log10_gen_y < good_range_log10_ymin ) continue ;
      if ( log10_gen_y > good_range_log10_ymax ) continue ;




      h_log10_q2_gen_vs_obs_dnn -> Fill( log10(dnn_Q2), log10_gen_Q2, wgt ) ;
      h_log10_q2_gen_vs_obs_e -> Fill( log10(obs_Q2_e), log10_gen_Q2, wgt ) ;
      h_log10_q2_gen_vs_obs_isigma -> Fill( log10(obs_Q2_ISigma), log10_gen_Q2, wgt ) ;
      h_log10_q2_gen_vs_obs_esigma -> Fill( log10(obs_Q2_eSigma), log10_gen_Q2, wgt ) ;
      h_log10_q2_gen_vs_obs_da -> Fill( log10(obs_Q2_DA), log10_gen_Q2, wgt ) ;

      h_log10_x_gen_vs_obs_dnn -> Fill( log10(dnn_x), log10_gen_x, wgt ) ;
      h_log10_x_gen_vs_obs_e -> Fill( log10(obs_x_e), log10_gen_x, wgt ) ;
      h_log10_x_gen_vs_obs_isigma -> Fill( log10(obs_x_ISigma), log10_gen_x, wgt ) ;
      h_log10_x_gen_vs_obs_esigma -> Fill( log10(obs_x_eSigma), log10_gen_x, wgt ) ;
      h_log10_x_gen_vs_obs_da -> Fill( log10(obs_x_DA), log10_gen_x, wgt ) ;

      h_log10_y_gen_vs_obs_dnn -> Fill( log10(dnn_y), log10_gen_y, wgt ) ;
      h_log10_y_gen_vs_obs_e -> Fill( log10(obs_y_e), log10_gen_y, wgt ) ;
      h_log10_y_gen_vs_obs_isigma -> Fill( log10(obs_y_ISigma), log10_gen_y, wgt ) ;
      h_log10_y_gen_vs_obs_esigma -> Fill( log10(obs_y_eSigma), log10_gen_y, wgt ) ;
      h_log10_y_gen_vs_obs_da -> Fill( log10(obs_y_DA), log10_gen_y, wgt ) ;

      h_x_gen_vs_obs_dnn -> Fill( (dnn_x), from_tlv_gen_x, wgt ) ;
      h_x_gen_vs_obs_e -> Fill( (obs_x_e), from_tlv_gen_x, wgt ) ;
      h_x_gen_vs_obs_isigma -> Fill( (obs_x_ISigma), from_tlv_gen_x, wgt ) ;
      h_x_gen_vs_obs_esigma -> Fill( (obs_x_eSigma), from_tlv_gen_x, wgt ) ;
      h_x_gen_vs_obs_da -> Fill( (obs_x_DA), from_tlv_gen_x, wgt ) ;

      rur_log10_q2_gen_vs_obs_dnn -> Fill( log10(dnn_Q2), log10_gen_Q2, wgt ) ;

      rur_log10_x_gen_vs_obs_dnn -> Fill( log10(dnn_x), log10_gen_x, wgt ) ;
      rur_log10_x_gen_vs_obs_e -> Fill( log10(obs_x_e), log10_gen_x, wgt ) ;
      rur_log10_x_gen_vs_obs_isigma -> Fill( log10(obs_x_ISigma), log10_gen_x, wgt ) ;
      rur_log10_x_gen_vs_obs_esigma -> Fill( log10(obs_x_eSigma), log10_gen_x, wgt ) ;
      rur_log10_x_gen_vs_obs_da -> Fill( log10(obs_x_DA), log10_gen_x, wgt ) ;

      rur_log10_y_gen_vs_obs_dnn -> Fill( log10(dnn_y), log10_gen_y, wgt ) ;
      rur_log10_y_gen_vs_obs_e -> Fill( log10(obs_y_e), log10_gen_y, wgt ) ;
      rur_log10_y_gen_vs_obs_isigma -> Fill( log10(obs_y_ISigma), log10_gen_y, wgt ) ;
      rur_log10_y_gen_vs_obs_esigma -> Fill( log10(obs_y_eSigma), log10_gen_y, wgt ) ;
      rur_log10_y_gen_vs_obs_da -> Fill( log10(obs_y_DA), log10_gen_y, wgt ) ;

      rur_2D_log10_q2_vs_log10_x_dnn    -> Fill( log10(dnn_x), log10(dnn_Q2),   log10_gen_x, log10_gen_Q2,  wgt ) ;
      rur_2D_log10_q2_vs_log10_x_e      -> Fill( log10(obs_x_e), log10(obs_Q2_e),   log10_gen_x, log10_gen_Q2,  wgt ) ;
      rur_2D_log10_q2_vs_log10_x_isigma -> Fill( log10(obs_x_ISigma), log10(obs_Q2_ISigma),   log10_gen_x, log10_gen_Q2,  wgt ) ;
      rur_2D_log10_q2_vs_log10_x_esigma -> Fill( log10(obs_x_eSigma), log10(obs_Q2_eSigma),   log10_gen_x, log10_gen_Q2,  wgt ) ;
      rur_2D_log10_q2_vs_log10_x_da     -> Fill( log10(obs_x_DA), log10(obs_Q2_DA),   log10_gen_x, log10_gen_Q2,  wgt ) ;



      rur_2D_log10_q2_vs_log10_y_dnn    -> Fill( log10(dnn_y), log10(dnn_Q2),   log10_gen_y, log10_gen_Q2,  wgt ) ;
      rur_2D_log10_q2_vs_log10_y_e      -> Fill( log10(obs_y_e), log10(obs_Q2_e),   log10_gen_y, log10_gen_Q2,  wgt ) ;
      rur_2D_log10_q2_vs_log10_y_isigma -> Fill( log10(obs_y_ISigma), log10(obs_Q2_ISigma),   log10_gen_y, log10_gen_Q2,  wgt ) ;
      rur_2D_log10_q2_vs_log10_y_esigma -> Fill( log10(obs_y_eSigma), log10(obs_Q2_eSigma),   log10_gen_y, log10_gen_Q2,  wgt ) ;
      rur_2D_log10_q2_vs_log10_y_da     -> Fill( log10(obs_y_DA), log10(obs_Q2_DA),   log10_gen_y, log10_gen_Q2,  wgt ) ;






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
   tf.WriteTObject( rur_log10_x_gen_vs_obs_esigma, rur_log10_x_gen_vs_obs_esigma->GetName() ) ;
   tf.WriteTObject( rur_log10_x_gen_vs_obs_da, rur_log10_x_gen_vs_obs_da->GetName() ) ;

   tf.WriteTObject( rur_log10_y_gen_vs_obs_dnn, rur_log10_y_gen_vs_obs_dnn->GetName() ) ;
   tf.WriteTObject( rur_log10_y_gen_vs_obs_e, rur_log10_y_gen_vs_obs_e->GetName() ) ;
   tf.WriteTObject( rur_log10_y_gen_vs_obs_isigma, rur_log10_y_gen_vs_obs_isigma->GetName() ) ;
   tf.WriteTObject( rur_log10_y_gen_vs_obs_esigma, rur_log10_y_gen_vs_obs_esigma->GetName() ) ;
   tf.WriteTObject( rur_log10_y_gen_vs_obs_da, rur_log10_y_gen_vs_obs_da->GetName() ) ;

   tf.WriteTObject( rur_2D_log10_q2_vs_log10_x_dnn, rur_2D_log10_q2_vs_log10_x_dnn->GetName() ) ;
   tf.WriteTObject( rur_2D_log10_q2_vs_log10_x_e, rur_2D_log10_q2_vs_log10_x_e->GetName() ) ;
   tf.WriteTObject( rur_2D_log10_q2_vs_log10_x_isigma, rur_2D_log10_q2_vs_log10_x_isigma->GetName() ) ;
   tf.WriteTObject( rur_2D_log10_q2_vs_log10_x_esigma, rur_2D_log10_q2_vs_log10_x_esigma->GetName() ) ;
   tf.WriteTObject( rur_2D_log10_q2_vs_log10_x_da, rur_2D_log10_q2_vs_log10_x_da->GetName() ) ;

   tf.WriteTObject( rur_2D_log10_q2_vs_log10_y_dnn, rur_2D_log10_q2_vs_log10_y_dnn->GetName() ) ;
   tf.WriteTObject( rur_2D_log10_q2_vs_log10_y_e, rur_2D_log10_q2_vs_log10_y_e->GetName() ) ;
   tf.WriteTObject( rur_2D_log10_q2_vs_log10_y_isigma, rur_2D_log10_q2_vs_log10_y_isigma->GetName() ) ;
   tf.WriteTObject( rur_2D_log10_q2_vs_log10_y_esigma, rur_2D_log10_q2_vs_log10_y_esigma->GetName() ) ;
   tf.WriteTObject( rur_2D_log10_q2_vs_log10_y_da, rur_2D_log10_q2_vs_log10_y_da->GetName() ) ;


   tf.Close() ;

} // Loop








