
#include "draw_obs_in_gen_slices.c"

   void run_doigs( const char* hname = "h_log10_x_gen_vs_obs_e", int last_slice_bin = 18, int first_slice_bin = 2 ) {

      TH2F* hp = (TH2F*) gDirectory -> FindObject( hname ) ;
      draw_obs_in_gen_slices( hp, last_slice_bin, first_slice_bin ) ;

   }
