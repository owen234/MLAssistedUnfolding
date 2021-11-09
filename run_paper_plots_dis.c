
#include "paper_plots_dis_one_method.c"
#include "paper_plots_dis_compare3.c"

   void run_paper_plots_dis(
            const char* input_file = "example-input-nbins_gen010_obs020.root",
            int nevts = 100000
          ) {

      char answ ;

    //---

      paper_plots_dis_one_method("h_log10_x_gen_vs_obs_dnn","DNN","x",nevts,input_file) ;
      answ = getchar() ; if ( answ == 'q' ) return ;

      paper_plots_dis_one_method("h_log10_x_gen_vs_obs_esigma","Sigma","x",nevts,input_file) ;
      answ = getchar() ; if ( answ == 'q' ) return ;

      paper_plots_dis_one_method("h_log10_x_gen_vs_obs_e","electron","x",nevts,input_file) ;
      answ = getchar() ; if ( answ == 'q' ) return ;

    //---

      paper_plots_dis_one_method("h_log10_y_gen_vs_obs_dnn","DNN","y",nevts,input_file) ;
      answ = getchar() ; if ( answ == 'q' ) return ;

      paper_plots_dis_one_method("h_log10_y_gen_vs_obs_esigma","Sigma","y",nevts,input_file) ;
      answ = getchar() ; if ( answ == 'q' ) return ;

      paper_plots_dis_one_method("h_log10_y_gen_vs_obs_e","electron","y",nevts,input_file) ;
      answ = getchar() ; if ( answ == 'q' ) return ;

    //---

      paper_plots_dis_compare3("h_log10_x_gen_vs_obs_dnn","h_log10_x_gen_vs_obs_esigma","h_log10_x_gen_vs_obs_e","x",100000,"example-input-nbins_gen010_obs020.root") ;
      answ = getchar() ; if ( answ == 'q' ) return ;

      paper_plots_dis_compare3("h_log10_y_gen_vs_obs_dnn","h_log10_y_gen_vs_obs_esigma","h_log10_y_gen_vs_obs_e","y",100000,"example-input-nbins_gen010_obs020.root") ;
      answ = getchar() ; if ( answ == 'q' ) return ;




   }
