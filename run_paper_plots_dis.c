
#include "paper_plots_dis_one_method.c"
#include "paper_plots_dis_compare3.c"
#include "paper_plots_dis_gen_syst.c"
#include "paper_plots_dis_gen_response_comp.c"

   void run_paper_plots_dis(
            const char* input_file = "paper-plots-input-1D-nbins_gen010_obs020.root",
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

      TCanvas* can ;

      can = (TCanvas*) gDirectory -> FindObject( "can1" ) ; delete can ;
      can = (TCanvas*) gDirectory -> FindObject( "can2" ) ; delete can ;
      can = (TCanvas*) gDirectory -> FindObject( "can3" ) ; delete can ;
      can = (TCanvas*) gDirectory -> FindObject( "can4" ) ; delete can ;




    //---

      paper_plots_dis_compare3("h_log10_x_gen_vs_obs_dnn","h_log10_x_gen_vs_obs_esigma","h_log10_x_gen_vs_obs_e","x",100000, input_file ) ;
      answ = getchar() ; if ( answ == 'q' ) return ;

      paper_plots_dis_compare3("h_log10_y_gen_vs_obs_dnn","h_log10_y_gen_vs_obs_esigma","h_log10_y_gen_vs_obs_e","y",100000, input_file, 1500 ) ;
      answ = getchar() ; if ( answ == 'q' ) return ;

      can = (TCanvas*) gDirectory -> FindObject( "can1" ) ; delete can ;
      can = (TCanvas*) gDirectory -> FindObject( "can2" ) ; delete can ;
      can = (TCanvas*) gDirectory -> FindObject( "can3" ) ; delete can ;




    //---

      paper_plots_dis_gen_syst( "x", 1, 60000, "paper-plots-input-1D-nbins_gen010_obs020.root", "paper-plots-input-1D-nbins_gen010_obs020-django.root"   ) ;
      answ = getchar() ; if ( answ == 'q' ) return ;

      paper_plots_dis_gen_syst( "y", 1, 60000, "paper-plots-input-1D-nbins_gen010_obs020.root", "paper-plots-input-1D-nbins_gen010_obs020-django.root"   ) ;
      answ = getchar() ; if ( answ == 'q' ) return ;


      can = (TCanvas*) gDirectory -> FindObject( "can1" ) ; delete can ;
      can = (TCanvas*) gDirectory -> FindObject( "can2" ) ; delete can ;
      can = (TCanvas*) gDirectory -> FindObject( "can3" ) ; delete can ;
      can = (TCanvas*) gDirectory -> FindObject( "can4" ) ; delete can ;


    //---

      paper_plots_dis_gen_response_comp("h_log10_x_gen_vs_obs_dnn","DNN","x", -1.00000,"paper-plots-input-1D-nbins_gen010_obs020.root","paper-plots-input-1D-nbins_gen010_obs020-django.root","Rapgap","Djangoh") ;
      answ = getchar() ; if ( answ == 'q' ) return ;

      paper_plots_dis_gen_response_comp("h_log10_x_gen_vs_obs_e","electron","x", -1.00000,"paper-plots-input-1D-nbins_gen010_obs020.root","paper-plots-input-1D-nbins_gen010_obs020-django.root","Rapgap","Djangoh") ;
      answ = getchar() ; if ( answ == 'q' ) return ;

      paper_plots_dis_gen_response_comp("h_log10_x_gen_vs_obs_esigma","Sigma","x", -1.00000,"paper-plots-input-1D-nbins_gen010_obs020.root","paper-plots-input-1D-nbins_gen010_obs020-django.root","Rapgap","Djangoh") ;
      answ = getchar() ; if ( answ == 'q' ) return ;


      paper_plots_dis_gen_response_comp("h_log10_y_gen_vs_obs_dnn","DNN","y", -1.00000,"paper-plots-input-1D-nbins_gen010_obs020.root","paper-plots-input-1D-nbins_gen010_obs020-django.root","Rapgap","Djangoh") ;
      answ = getchar() ; if ( answ == 'q' ) return ;

      paper_plots_dis_gen_response_comp("h_log10_y_gen_vs_obs_e","electron","y", -1.00000,"paper-plots-input-1D-nbins_gen010_obs020.root","paper-plots-input-1D-nbins_gen010_obs020-django.root","Rapgap","Djangoh") ;
      answ = getchar() ; if ( answ == 'q' ) return ;

      paper_plots_dis_gen_response_comp("h_log10_y_gen_vs_obs_esigma","Sigma","y", -1.00000,"paper-plots-input-1D-nbins_gen010_obs020.root","paper-plots-input-1D-nbins_gen010_obs020-django.root","Rapgap","Djangoh") ;
      answ = getchar() ; if ( answ == 'q' ) return ;



   }








