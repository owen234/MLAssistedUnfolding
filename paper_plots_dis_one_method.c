

//
//  Following the example of tutorials/unfold/testUnfold1.c
//


//
//  This makes one plot per canvas and saves each plot to its own pdf file.
//
//

#include "histio.c"

#ifndef helpers
#define helpers

TH2F* get_hist( const char* hname ) {
   TH2F* hp = (TH2F*) gDirectory -> FindObject( hname ) ;
   if ( hp == 0x0 ) { printf("\n\n *** can't find %s\n\n", hname ) ; gSystem->Exit(-1) ; }
   return hp ;
}

#ifndef MyUnf
#define MyUnf
class MyTUnfoldDensity : public TUnfoldDensity {
   using TUnfoldDensity::TUnfoldDensity;
public:
   inline const TMatrixDSparse *MyGetE   (void) const { return GetE(); }
   inline const TMatrixDSparse *MyGetEinv(void) const { return GetEinv(); }
   /// folded back result
   inline const TMatrixDSparse *MyGetAx(void) const { return fA; }
   inline const TMatrixDSparse *MyGetDXDY(void) const { return GetDXDY(); }
   //   TMatrixDSparse *fDXDYVyy = MultiplyMSparseMSparse(fDXDY,fVyy)
   //   fVxx = MultiplyMSparseMSparseTranspVector(fDXDYVyy,fDXDY,0);
   // posterior response matrix
   inline const TMatrixDSparse *GetPosteriorResponseMatrix(void) const { return MultiplyMSparseMSparse(GetDXDY(),fA); } 
   TH2D GetPosteriorResponseHist(double zmin=0) const {
      TH1 *output = this->GetOutput("temp");
      TH2D PP("posterior response matrix","Posterior response matrix;Gen;Obs",
              output->GetNbinsX(),output->GetXaxis()->GetXbins()->GetArray(),
              output->GetNbinsX(),output->GetXaxis()->GetXbins()->GetArray() );
      const TMatrixDSparse *MatP = GetPosteriorResponseMatrix();
      for ( int ii = 0 ; ii<PP.GetNbinsX() ;ii++ ) {
         for ( int jj = 0 ; jj<PP.GetNbinsX() ;jj++ ) {
            double zval =(*MatP)[ii][jj];
            if ( zval != 0 && zval < zmin) zval = 0;//max(zval,zmin);
            PP.SetBinContent(ii+1,jj+1,zval);
         }
      }
      PP.SetMinimum(zmin);
      delete output;
      return PP;//MultiplyMSparseMSparse(GetDXDY(),fA);
   } 
};
#endif



double Calc_r(const TMatrixDSparse* AA) {
   const TMatrixDSparse& A = (*AA);
      //TMatrixT (Int_t nrows, Int_t ncols)
      const int NN = AA->GetNcols();
      TMatrixD oneT (  1, NN ); // j
      TMatrixD ones ( NN,  1 ); // jT
      TMatrixD numT (  1, NN ); // r
      TMatrixD nums ( NN,  1 ); // rT
      TMatrixD num2T(  1, NN ); // r2
      TMatrixD num2s( NN,  1 ); // r2T
      for ( int ii = 0 ; ii<NN; ii++ ) {
         oneT [0][ii] = 1;
         ones [ii][0] = 1;
         numT [0][ii] = ii+1;
         nums [ii][0] = ii+1;
         num2T[0][ii] = pow(ii+1,2);
         num2s[ii][0] = pow(ii+1,2);
      }
      double n   = (oneT  * A * ones )[0][0];
      double Sx  = (numT  * A * ones )[0][0];
      double Sy  = (oneT  * A * nums )[0][0];
      double Sx2 = (num2T * A * ones )[0][0];
      double Sy2 = (oneT  * A * num2s)[0][0];
      double Sxy = (numT  * A * nums )[0][0];

      double num    = n * Sxy - Sx*Sy;
      double denom2 = (n*Sx2 - Sx*Sx ) * (n*Sy2 - Sy*Sy);
      double denom  = sqrt(denom2);
      double r      = num/denom;
      return r;
}

double Calc_r(const TMatrixD& AA) {
   TMatrixDSparse BB(AA);
   return Calc_r(&BB);
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

   void paper_plots_dis_one_method(
                      const char* hist_name_a = "h_log10_x_gen_vs_obs_dnn",
                      const char* method_name = "DNN",
                      const char* var_name = "x",
                      int ngen = 1e5,
                      const char* input_file = "example-input-nbins_gen010_obs020.root"
                      ) {

      char htitle[1000] ;
      char fname[1000] ;

      float label_size = 0.050 ;
      float title_size = 0.055 ;
      float label_offset_y = 0.015 ;
      float label_offset_x = 0.015 ;
      float label_offset_z = 0.010 ;
      float title_offset_x = 1.1 ;
      float title_offset_y = 1.57 ;
      float title_x = 0.05 ;
      float title_y = 0.94 ;


      gSystem -> Exec( "mkdir -p paper-plots" ) ;

      TRandom3* tran = new TRandom3(1249) ;

      gDirectory -> Delete( "h*" ) ;
      gStyle -> SetPalette( kBird ) ;

      ////gStyle -> SetTitleW( 1.00 ) ;
      ////gStyle -> SetTitleH( 0.06 ) ;

      TText* tt_title = new TText() ;
      tt_title -> SetTextSize( 0.05 ) ;

      loadHist( input_file ) ;

      printf("\n\n") ;
      gDirectory -> ls() ;
      printf("\n\n") ;

      TH2F* h_in_gen_vs_obs_a = get_hist( hist_name_a ) ;

      TH1* h_obs_source_a = h_in_gen_vs_obs_a -> ProjectionX( "h_obs_source_a" ) ;
      TH1* h_gen_source_a = h_in_gen_vs_obs_a -> ProjectionY( "h_gen_source_a" ) ;


      TH2F* h_normalized_response = (TH2F*) h_in_gen_vs_obs_a -> Clone( "h_normalized_response" ) ;
      for ( int ybi=1; ybi<=h_in_gen_vs_obs_a->GetNbinsX(); ybi++ ) {
         float row_sum = 0. ;
         for ( int xbi=1; xbi<=h_in_gen_vs_obs_a->GetNbinsX(); xbi++ ) {
            row_sum += h_in_gen_vs_obs_a -> GetBinContent( xbi, ybi ) ;
         } // xbi
         for ( int xbi=1; xbi<=h_in_gen_vs_obs_a->GetNbinsX(); xbi++ ) {
            if ( row_sum > 0 ) {
               h_normalized_response -> SetBinContent( xbi, ybi, ( h_in_gen_vs_obs_a -> GetBinContent( xbi, ybi ) ) / row_sum ) ;
               h_normalized_response -> SetBinError( xbi, ybi, 0. ) ;
            } else {
               h_normalized_response -> SetBinContent( xbi, ybi, 0. ) ;
               h_normalized_response -> SetBinError( xbi, ybi, 0. ) ;
            }
         } // xbi
      } // ybi


      TH1* h_obs_random_a = (TH1*) h_obs_source_a -> Clone( "h_obs_random_a" ) ;
      h_obs_random_a->Reset() ;
      h_obs_random_a->FillRandom( h_obs_source_a, ngen, tran ) ;




     //-- Add underflows and overflows to fake data!
      float underflow_mean ;
      float overflow_mean ;
      int nbins_obs = h_obs_source_a -> GetNbinsX() ;
      float nobs_source ;
      float nobs_integral ;
      float nobs_underflow ;
      float nobs_overflow ;



      nobs_integral = h_obs_source_a -> Integral() ;
      nobs_underflow = h_obs_source_a -> GetBinContent(0) ;
      nobs_overflow = h_obs_source_a -> GetBinContent( nbins_obs+1 ) ;

      printf("  h_obs_source_a :  integral %9.5f, underflow %9.5f, overflow %9.5f\n", nobs_integral, nobs_underflow, nobs_overflow ) ;

      nobs_source = nobs_integral ;
      nobs_source += nobs_underflow ;
      nobs_source += nobs_overflow ;

      underflow_mean = ngen * ( nobs_underflow / nobs_integral ) ; // ngen only fills inside of histogram (integral)
      overflow_mean = ngen * ( nobs_overflow / nobs_integral ) ;

      h_obs_random_a -> SetBinContent( 0, tran->Poisson( underflow_mean ) ) ;
      h_obs_random_a -> SetBinContent( nbins_obs+1, tran->Poisson( overflow_mean ) ) ;






      MyTUnfoldDensity unfold_a( h_in_gen_vs_obs_a, TUnfold::kHistMapOutputVert ) ;

      int return_status_a = unfold_a.SetInput( h_obs_random_a ) ;
      printf("  Return status for SetInput A: %d\n", return_status_a ) ;




      //========================================================================
      // the unfolding is done here
      //
      // scan L curve and find best point
      Int_t nScan=30;
      // use automatic L-curve scan: start with taumin=taumax=0.0
      Double_t tauMin=0.0;
      Double_t tauMax=0.0;
      Int_t iBest;
      TSpline *logTauX_a,*logTauY_a;
      TSpline *logTauX_b,*logTauY_b;
      TGraph *lCurve_a;
      TGraph *lCurve_b;

      // if required, report Info messages (for debugging the L-curve scan)
      Int_t oldinfo=gErrorIgnoreLevel;
      gErrorIgnoreLevel=kInfo;

      // this method scans the parameter tau and finds the kink in the L curve
      // finally, the unfolding is done for the best choice of tau
      iBest=unfold_a.ScanLcurve( nScan, tauMin, tauMax, &lCurve_a, &logTauX_a, &logTauY_a );

      // if required, switch to previous log-level
      gErrorIgnoreLevel=oldinfo;


      //==========================================================================
      // print some results
      //
      std::cout<<"A tau="<<unfold_a.GetTau()<<"\n";
      std::cout<<"A chi**2="<<unfold_a.GetChi2A()<<"+"<<unfold_a.GetChi2L() <<" / "<<unfold_a.GetNdf()<<"\n";
      std::cout<<"A rho_avg="<<unfold_a.GetRhoAvg()<<"\trho_max / "<<unfold_a.GetRhoMax()<<"\n";

      
      //==========================================================================
      // retrieve results into histograms

      // get unfolded distribution
      TH1 *histMunfold_a = unfold_a.GetOutput("Unfolded_a");
      sprintf( htitle, "Unfolded distribution, %s", method_name ) ;
      histMunfold_a -> SetTitle( htitle ) ;
      {
         TString xtitle(histMunfold_a->GetXaxis()->GetTitle());
         xtitle.ReplaceAll("TAU1B","#tau_{ 1}^{ b}");
         histMunfold_a->GetXaxis()->SetTitle(xtitle);
      }
      // get unfolding result, folded back
      TH1 *histMdetFold_a = unfold_a.GetFoldedOutput("FoldedBack_a");


      // get total error matrix:
      //   migration matrix uncorrelated and correlated systematic errors
      //   added in quadrature to the data statistical errors
      TH2 *histEmatTotal_a=unfold_a.GetEmatrixTotal("EmatTotal_a");
      histEmatTotal_a -> SetTitle( "Covariance matrix, a" ) ;


      //-- calculate the correlation coefficients matrix
      TH2* correlation_matrix_a = (TH2*) histEmatTotal_a -> Clone( "correlation_matrix_a" ) ;
      sprintf( htitle, "Correlation coefficients, %s", method_name ) ;
      correlation_matrix_a -> SetTitle( htitle ) ;
      {
         TString xtitle(correlation_matrix_a->GetXaxis()->GetTitle());
         xtitle.ReplaceAll("TAU1B","#tau_{ 1}^{ b}");
         correlation_matrix_a->GetXaxis()->SetTitle(xtitle);
      }
      
      for ( int xbi=1; xbi<=histEmatTotal_a->GetNbinsX(); xbi++ ) {
         for ( int ybi=1; ybi<=histEmatTotal_a->GetNbinsY(); ybi++ ) {

            float rho = 1. ;
            float sigma_x2, sigma_y2 ;

            sigma_x2 = histEmatTotal_a->GetBinContent( xbi, xbi ) ;
            sigma_y2 = histEmatTotal_a->GetBinContent( ybi, ybi ) ;
            if ( sigma_x2 > 0 && sigma_y2 > 0 ) {
               rho = histEmatTotal_a -> GetBinContent( xbi, ybi )  / sqrt( sigma_x2 * sigma_y2 ) ;
            }
            correlation_matrix_a -> SetBinContent( xbi, ybi, rho ) ;


         } // ybi
      } // xbi

      correlation_matrix_a -> SetMinimum( -1. ) ;

      correlation_matrix_a -> SetMaximum(  1. ) ;



      // create data histogram with the total errors
      int nGen = histMunfold_a -> GetNbinsX() ;
      float xminGen = histMunfold_a -> GetXaxis() -> GetXmin() ;
      float xmaxGen = histMunfold_a -> GetXaxis() -> GetXmax() ;

      TH1D *histTotalError_a = new TH1D("TotalError_a","TotalError_a",nGen,xminGen,xmaxGen);
      for(Int_t bin=1;bin<=nGen;bin++) {
        histTotalError_a->SetBinContent(bin,histMunfold_a->GetBinContent(bin));
        histTotalError_a->SetBinError(bin,TMath::Sqrt(histEmatTotal_a->GetBinContent(bin,bin)));
      }


      // get global correlation coefficients
      // for this calculation one has to specify whether the
      // underflow/overflow bins are included or not
      // default: include all bins
      // here: exclude underflow and overflow bins
      TH2 *gHistInvEMatrix;
      TH1 *histRhoi_a = unfold_a.GetRhoItotal("rho_I_a",
                                        0, // use default title
                                        0, // all distributions
                                        "*[UO]", // discard underflow and overflow bins on all axes
                                        kTRUE, // use original binning
                                        &gHistInvEMatrix // store inverse of error matrix
                                        );

      histRhoi_a -> SetTitle( "Global correlation, a" ) ;


      gStyle -> SetOptStat(0) ;

      h_obs_source_a -> SetLineColor(4) ;
      h_gen_source_a -> SetLineColor(2) ;

      TH1* h_gen_compare_a = (TH1*) h_gen_source_a -> Clone( "h_gen_compare_a" ) ;
      h_gen_compare_a -> Scale( ( histMunfold_a -> Integral() )/( h_gen_source_a -> Integral() ) ) ;




      TH1* h_unfold_err_a = (TH1*) histMunfold_a -> Clone( "h_unfold_err_a" ) ;


      h_unfold_err_a -> SetFillColor( kBlue-9 ) ;
      h_unfold_err_a -> SetFillStyle( 1001 ) ;


      for ( int bi=1; bi<=histMunfold_a->GetNbinsX(); bi++ ) {
         h_unfold_err_a -> SetBinContent( bi, 0. ) ;
      }




      printf("\n\n") ;
      printf(" Gen bin edges:\n") ;
      for ( int bi=1; bi<=(histMunfold_a->GetNbinsX()+1); bi++ ) {
         float lowEdgeLog10 = histMunfold_a->GetBinLowEdge( bi ) ;
         float lowEdge = pow( 10., lowEdgeLog10 ) ;
         printf(" %5.3f  (%7.5f) |", lowEdgeLog10, lowEdge ) ;
      } // bi
      printf("\n\n") ;






      histRhoi_a->SetMaximum(1.1) ;


      //gStyle -> SetPadRightMargin(0.18) ;
      //gStyle -> SetPadLeftMargin(0.18) ;
      gStyle -> SetPadRightMargin(0.15) ;
      gStyle -> SetPadLeftMargin(0.20) ;
      gStyle -> SetPadBottomMargin(0.18) ;
      gStyle -> SetTitleBorderSize(0) ;
      gStyle -> SetTitleY(0.975) ;



      int can_width = 600 ;
      int can_height = 600 ;

      int can_spacing = 20 ;


     //-----

      TCanvas* can1 = (TCanvas*) gDirectory -> FindObject( "can1" ) ;
      if ( can1 == 0x0 ) { printf( "Making can1\n") ; can1 = new TCanvas( "can1", "", can_spacing, can_spacing, can_width, can_height ) ; }
      can1 -> Clear() ;
      can1 -> cd() ;


      //////////////h_in_gen_vs_obs_a -> SetTitleOffset( 1.2, "x" ) ;
      //////////////h_in_gen_vs_obs_a -> SetTitleOffset( 1.6, "y" ) ;
      //////////////h_in_gen_vs_obs_a -> SetTitleSize( 0.045, "x" ) ;
      //////////////h_in_gen_vs_obs_a -> SetTitleSize( 0.045, "y" ) ;

      h_in_gen_vs_obs_a -> SetTitleOffset( title_offset_x, "x" ) ;
      h_in_gen_vs_obs_a -> SetTitleOffset( title_offset_y, "y" ) ;
      h_in_gen_vs_obs_a -> SetTitleSize( title_size, "x" ) ;
      h_in_gen_vs_obs_a -> SetTitleSize( title_size, "y" ) ;
      h_in_gen_vs_obs_a -> SetLabelSize( label_size, "x" ) ;
      h_in_gen_vs_obs_a -> SetLabelSize( label_size, "y" ) ;
      h_in_gen_vs_obs_a -> SetLabelOffset( label_offset_x, "x" ) ;
      h_in_gen_vs_obs_a -> SetLabelOffset( label_offset_y, "y" ) ;


      sprintf( htitle, "Response matrix, %s", method_name ) ;
      h_in_gen_vs_obs_a -> SetTitle( htitle ) ;
      {
         TString xtitle(h_in_gen_vs_obs_a->GetXaxis()->GetTitle());
         xtitle.ReplaceAll("TAU1B","#tau_{ 1}^{ b}");
         h_in_gen_vs_obs_a->GetXaxis()->SetTitle(xtitle);
         TString ytitle(h_in_gen_vs_obs_a->GetYaxis()->GetTitle());
         ytitle.ReplaceAll("TAU1B","#tau_{ 1}^{ b}");
         h_in_gen_vs_obs_a->GetYaxis()->SetTitle(ytitle);
      }
      if ( strcmp( var_name, "y" ) == 0 ) { h_in_gen_vs_obs_a -> SetNdivisions( 605 ) ; }



      h_in_gen_vs_obs_a -> Draw("colz") ;
      TExec* change_hist_palette = new TExec( "change_hist_palette", "Setup2DhistPalette();" );
      change_hist_palette->Draw() ;
      h_in_gen_vs_obs_a -> Draw("colz same") ;
      h_in_gen_vs_obs_a -> Draw("axis same") ;

      h_in_gen_vs_obs_a -> GetZaxis() -> SetLabelSize( label_size ) ;
      h_in_gen_vs_obs_a -> GetZaxis() -> SetLabelOffset( label_offset_z ) ;

      h_in_gen_vs_obs_a -> SetTitle( "" ) ;
      tt_title -> DrawTextNDC( title_x, title_y, htitle ) ;

      can1 -> Update() ; can1 -> Draw() ; gSystem -> ProcessEvents() ;

      sprintf( fname, "paper-plots/dis-response-%s-%s.pdf", var_name, method_name ) ;
      can1 -> SaveAs( fname ) ;
      sprintf( fname, "paper-plots/dis-response-%s-%s.png", var_name, method_name ) ;
      can1 -> SaveAs( fname ) ;


     //-----

      gStyle -> SetPadRightMargin(0.05) ;

      TCanvas* can2 = (TCanvas*) gDirectory -> FindObject( "can2" ) ;
      if ( can2 == 0x0 ) { printf( "Making can2\n") ; can2 = new TCanvas( "can2", "", 2*can_spacing + can_width , can_spacing, can_width, can_height ) ; }
      can2 -> Clear() ;
      can2 -> cd() ;

      /////////histMunfold_a -> SetTitleOffset( 1.2, "x" ) ;
      /////////histMunfold_a -> SetTitleOffset( 1.8, "y" ) ;
      /////////histMunfold_a -> SetTitleSize( 0.045, "x" ) ;
      /////////histMunfold_a -> SetTitleSize( 0.045, "y" ) ;

      histMunfold_a -> SetTitleOffset( title_offset_x, "x" ) ;
      histMunfold_a -> SetTitleOffset( 1.2*title_offset_y, "y" ) ;
      histMunfold_a -> SetTitleSize( title_size, "x" ) ;
      histMunfold_a -> SetTitleSize( title_size, "y" ) ;
      histMunfold_a -> SetLabelSize( label_size, "x" ) ;
      histMunfold_a -> SetLabelSize( label_size, "y" ) ;
      histMunfold_a -> SetLabelOffset( label_offset_x, "x" ) ;
      histMunfold_a -> SetLabelOffset( label_offset_y, "y" ) ;


      histMunfold_a -> SetYTitle( "Events" ) ;

      if ( strcmp( var_name, "x" ) == 0 ) { histMunfold_a -> SetXTitle( "log10(x)" ) ; }
      if ( strcmp( var_name, "y" ) == 0 ) { histMunfold_a -> SetXTitle( "log10(y)" ) ; }
      if ( strcmp( var_name, "y" ) == 0 ) { histMunfold_a -> SetNdivisions( 605 ) ; }

      histMunfold_a -> SetLineWidth(3) ;
      h_gen_compare_a -> SetLineWidth(3) ;


      histMunfold_a -> Draw() ;
      h_gen_compare_a -> Draw("same hist") ;
      histMunfold_a -> Draw( "same" ) ;

      TLegend* legend = new TLegend( 0.60, 0.40, 0.90, 0.55 ) ;
      legend -> AddEntry( histMunfold_a, "Unfolded" ) ;
      legend -> AddEntry( h_gen_compare_a, "Gen" ) ;
      legend -> Draw() ;

      sprintf( htitle, "%s", histMunfold_a -> GetTitle() ) ;
      histMunfold_a -> SetTitle( "" ) ;
      tt_title -> DrawTextNDC( title_x, title_y, htitle ) ;

      can2 -> Update() ; can2 -> Draw() ; gSystem -> ProcessEvents() ;

      sprintf( fname, "paper-plots/dis-unfolded-%s-%s.pdf", var_name, method_name ) ;
      can2 -> SaveAs( fname ) ;
      sprintf( fname, "paper-plots/dis-unfolded-%s-%s.png", var_name, method_name ) ;
      can2 -> SaveAs( fname ) ;

     //-----

      gStyle -> SetPadRightMargin(0.18) ;

      TCanvas* can3 = (TCanvas*) gDirectory -> FindObject( "can3" ) ;
      if ( can3 == 0x0 ) { printf( "Making can3\n") ; can3 = new TCanvas( "can3", "", 3*can_spacing + 2*can_width , can_spacing, can_width, can_height ) ; }
      can3 -> Clear() ;
      can3 -> cd() ;

      /////////correlation_matrix_a -> SetTitleOffset( 1.2, "x" ) ;
      /////////correlation_matrix_a -> SetTitleOffset( 1.6, "y" ) ;
      /////////correlation_matrix_a -> SetTitleSize( 0.045, "x" ) ;
      /////////correlation_matrix_a -> SetTitleSize( 0.045, "y" ) ;

      correlation_matrix_a -> SetTitleOffset( title_offset_x, "x" ) ;
      correlation_matrix_a -> SetTitleOffset( title_offset_y, "y" ) ;
      correlation_matrix_a -> SetTitleSize( title_size, "x" ) ;
      correlation_matrix_a -> SetTitleSize( title_size, "y" ) ;
      correlation_matrix_a -> SetLabelSize( label_size, "x" ) ;
      correlation_matrix_a -> SetLabelSize( label_size, "y" ) ;
      correlation_matrix_a -> SetLabelOffset( label_offset_x, "x" ) ;
      correlation_matrix_a -> SetLabelOffset( label_offset_y, "y" ) ;

      if ( strcmp( var_name, "x" ) == 0 ) { correlation_matrix_a -> SetXTitle( "log10(x)" ) ; correlation_matrix_a -> SetYTitle( "log10(x)" ) ; }
      if ( strcmp( var_name, "y" ) == 0 ) { correlation_matrix_a -> SetXTitle( "log10(y)" ) ; correlation_matrix_a -> SetYTitle( "log10(y)" ) ; }
      if ( strcmp( var_name, "y" ) == 0 ) { correlation_matrix_a -> SetNdivisions( 605 ) ; }


      correlation_matrix_a -> Draw( "colz" ) ;
      TExec* change_cor_palette = new TExec( "change_cor_palette", "SetupCorrelationPalette();" );
      change_cor_palette->Draw() ;
      correlation_matrix_a -> Draw( "colz same" ) ;
      correlation_matrix_a -> Draw( "axis same" ) ;

      correlation_matrix_a -> GetZaxis() -> SetLabelSize( label_size ) ;
      correlation_matrix_a -> GetZaxis() -> SetLabelOffset( label_offset_z ) ;

      sprintf( htitle, "%s", correlation_matrix_a -> GetTitle() ) ;
      correlation_matrix_a -> SetTitle( "" ) ;
      tt_title -> DrawTextNDC( title_x, title_y, htitle ) ;

      can3 -> Update() ; can3 -> Draw() ; gSystem -> ProcessEvents() ;

      sprintf( fname, "paper-plots/dis-cormat-%s-%s.pdf", var_name, method_name ) ;
      can3 -> SaveAs( fname ) ;
      sprintf( fname, "paper-plots/dis-cormat-%s-%s.png", var_name, method_name ) ;
      can3 -> SaveAs( fname ) ;

     //-----

      TCanvas* can4 = (TCanvas*) gDirectory -> FindObject( "can4" ) ;
      if ( can4 == 0x0 ) { printf( "Making can4\n") ; can4 = new TCanvas( "can4", "", 4*can_spacing + 3*can_width , can_spacing, can_width, can_height ) ; }
      can4 -> Clear() ;
      can4 -> cd() ;

      /////////h_normalized_response -> SetTitleOffset( 1.2, "x" ) ;
      /////////h_normalized_response -> SetTitleOffset( 1.6, "y" ) ;
      /////////h_normalized_response -> SetTitleSize( 0.045, "x" ) ;
      /////////h_normalized_response -> SetTitleSize( 0.045, "y" ) ;

      h_normalized_response -> SetTitleOffset( title_offset_x, "x" ) ;
      h_normalized_response -> SetTitleOffset( title_offset_y, "y" ) ;
      h_normalized_response -> SetTitleSize( title_size, "x" ) ;
      h_normalized_response -> SetTitleSize( title_size, "y" ) ;
      h_normalized_response -> SetLabelSize( label_size, "x" ) ;
      h_normalized_response -> SetLabelSize( label_size, "y" ) ;
      h_normalized_response -> SetLabelOffset( label_offset_x, "x" ) ;
      h_normalized_response -> SetLabelOffset( label_offset_y, "y" ) ;



      sprintf( htitle, "Normalized response matrix, %s", method_name ) ;
      ////h_normalized_response -> SetTitle( htitle ) ;
      h_normalized_response -> SetTitle( "" ) ;
      {
         TString xtitle(h_normalized_response->GetXaxis()->GetTitle());
         xtitle.ReplaceAll("TAU1B","#tau_{ 1}^{ b}");
         h_normalized_response->GetXaxis()->SetTitle(xtitle);
      }
      if ( strcmp( var_name, "y" ) == 0 ) { h_normalized_response -> SetNdivisions( 605 ) ; }

      h_normalized_response -> Draw("colz") ;
      TExec* change_hist_palette2 = new TExec( "change_hist_palette2", "Setup2DhistPalette();" );
      change_hist_palette2->Draw() ;
      h_normalized_response -> Draw("colz same") ;
      h_normalized_response -> Draw("axis same") ;

      h_normalized_response -> GetZaxis() -> SetLabelSize( label_size ) ;
      h_normalized_response -> GetZaxis() -> SetLabelOffset( label_offset_z ) ;

      tt_title -> DrawTextNDC( title_x, title_y, htitle ) ;

      can4 -> Update() ; can4 -> Draw() ; gSystem -> ProcessEvents() ;

      sprintf( fname, "paper-plots/dis-normalized-response-%s-%s.pdf", var_name, method_name ) ;
      can4 -> SaveAs( fname ) ;
      sprintf( fname, "paper-plots/dis-normalized-response-%s-%s.png", var_name, method_name ) ;
      can4 -> SaveAs( fname ) ;


      //==========================================================================
      // 
      // calculate posterior response matrix [arXiv:2203.09579]
      // P = M A
      // M = E A^T V^-1   with E = (ATWA+L)^-1
      
      const TMatrixDSparse *PosteriorResponseMatrix = unfold_a.GetPosteriorResponseMatrix();
      // PosteriorResponseMatrix->Print();

      
      TH2D PosteriorResponseHist = unfold_a.GetPosteriorResponseHist(1e-5);
      PosteriorResponseHist.SetXTitle( h_normalized_response->GetXaxis()->GetTitle()); 
      PosteriorResponseHist.SetYTitle( h_normalized_response->GetYaxis()->GetTitle()); 

      //TMatrixT (Int_t nrows, Int_t ncols)
      double r = Calc_r(PosteriorResponseMatrix);
      int NN = PosteriorResponseMatrix->GetNrows();
      
      
      //PosteriorResponseHist.Print("all");
      double trace = 0;
      for ( int ii = 0 ; ii<PosteriorResponseHist.GetNbinsX() ; ii++ ) trace+=PosteriorResponseHist.GetBinContent(ii+1,ii+1);
      PosteriorResponseHist.SetTitle( Form("Posterior response matrix, %s - Trace/n=%.3f", method_name, trace/NN ) ) ;

      //PosteriorResponseHist.SetTitle( Form("#splitline{Posterior response matrix, %s}{#scale[0.8]{Trace/n=%.2f,  1-r(A)=%.2e}}", method_name, trace/NN, 1-r ) ) ;
      // c.SetTopMargin(0.16);
      // c.SetBottomMargin(0.14);

      TCanvas c("c","c",800,800);
      PosteriorResponseHist.Draw("colz");
      c.SetLogz();
      c.SaveAs(Form("paper-plots/dis-posterior-response-%s-%s.pdf",var_name, method_name));
      c.SaveAs(Form("paper-plots/dis-posterior-response-%s-%s.png",var_name, method_name));

      // reset
      gStyle->SetTitleY(0.975);
            
     //-----


      

      printf("\n\n\n") ;

      printf(" cut and paste for this:\n\n") ;

      printf("     paper_plots_dis_one_method(\"%s\",\"%s\",\"%s\",%d,\"%s\")\n", hist_name_a, method_name, var_name, ngen, input_file ) ;


      printf("\n\n\n") ;

   }















void test_calc_r() {

      {
         TMatrixD B(3,3);
         B[0][0] = 6;
         B[1][0] = 1;
         B[2][0] = 0;
         B[0][1] = 1;
         B[1][1] = 5;
         B[2][1] = 2;
         B[0][2] = 1;
         B[1][2] = 3;
         B[2][2] = 6;
         B.Print();
         double rB = Calc_r(B);
         cout<<"rB: "<<rB<<endl;
      }

      {
         TMatrixD D(5,5);
         D[0][0] = 2;
         D[1][0] = 1;

         D[0][1] = 1;
         D[1][1] = 3;
         D[2][1] = 2;

         D[1][2] = 2;
         D[2][2] = 3;
         D[3][2] = 4;

         D[2][3] = 1;
         D[3][3] = 2;
         D[4][3] = 3;

         D[3][4] = 1;
         D[4][4] = 1;

         D.Print();
         double rD = Calc_r(D);
         cout<<"rD: "<<rD<<endl;
      }
         exit(0);

}
