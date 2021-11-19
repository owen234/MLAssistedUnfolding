//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Nov 19 00:28:50 2021 by ROOT version 6.24/06
// from TTree dnnout/
// found on file: dnn-output-h1-v2-Rapgap.root
//////////////////////////////////////////////////////////

#ifndef dnnout_tau1b_h
#define dnnout_tau1b_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class dnnout_tau1b {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Char_t          has_isr;
   Char_t          has_fsr;
   Float_t         tower_sum_40;
   Long64_t        n_towers_40;
   Float_t         eta_pho_closest_to_ebeam;
   Double_t        e_pho_closest_to_ebeam;
   Float_t         phi_pho_closest_to_ebeam;
   Float_t         obs_e_e;
   Float_t         obs_e_phi;
   Float_t         obs_e_theta;
   Float_t         obs_hfs_e;
   Float_t         obs_hfs_pz;
   Float_t         obs_hfs_pt;
   Float_t         obs_hfs_phi;
   Float_t         obs_DeltaPhi;
   Float_t         obs_e_trk_e;
   Float_t         beam_e_e;
   Float_t         beam_p_e;
   Float_t         gen_HFS_Sigma;
   Float_t         gen_HFS_T;
   Float_t         gen_eUncomb_E;
   Float_t         gen_eUncomb_theta;
   Float_t         gen_eRecomb_E;
   Float_t         gen_eRecomb_theta;
   Float_t         gen_tau1b;
   Float_t         obs_tau1bs;
   Float_t         obs_tau1be;
   Float_t         obs_tau1bda;
   Float_t         obs_cHFSs_M;
   Float_t         obs_cHFSs_pt;
   Float_t         obs_cHFSs_theta;
   Float_t         obs_cHFSs_phi;
   Float_t         obs_cHFSs_N;
   Float_t         obs_cHFSs_dRavg;
   Float_t         obs_cHFSs_dR2avg;
   Float_t         obs_cHFSs_Empz;
   Float_t         obs_cHFSs_Eppz;
   Float_t         obs_cHFSe_Empz;
   Float_t         obs_cHFSe_Eppz;
   Float_t         obs_cHFSe_keE0;
   Float_t         obs_cHFSs_keE0;
   Float_t         obs_cHFSs_kesE0;
   Float_t         wgt;
   Float_t         obs_e_pz;
   Float_t         obs_e_pt;
   Float_t         obs_hfs_Sigma;
   Float_t         obs_hfs_T;
   Float_t         obs_e_Sigma;
   Float_t         obs_e_tantheta;
   Float_t         obs_e_Pt2;
   Float_t         obs_hfs_tangamma;
   Float_t         obs_hfs_Empz;
   Float_t         obs_e_Empz;
   Float_t         obs_event_Empz;
   Float_t         rot_pt1;
   Float_t         rot_pt2;
   Float_t         rot_Empz1;
   Float_t         rot_Empz2;
   Double_t        e_ecal_over_trk_ratio;
   Double_t        dphi_pho_closest_to_ebeam;
   Bool_t          has_norad;
   Float_t         obs_kine_ye;
   Float_t         obs_kine_yda;
   Float_t         obs_kine_yh;
   Float_t         obs_kine_ys;
   Float_t         obs_kine_Q2e;
   Float_t         obs_kine_Q2da;
   Float_t         obs_kine_Q2s;
   Float_t         obs_kine_Q2ida;
   Float_t         obs_kine_xis;
   Float_t         gen_e_Sigma;
   Float_t         gen_e_tantheta;
   Float_t         gen_e_Pt2;
   Float_t         gen_HFS_tangamma;
   Float_t         gen_kine_ye;
   Float_t         gen_kine_yda;
   Float_t         gen_kine_ys;
   Float_t         gen_kine_Q2e;
   Float_t         gen_kine_Q2da;
   Float_t         gen_kine_Q2s;
   Float_t         gen_kine_xis;
   Float_t         gen_log_x;
   Float_t         gen_log_y;
   Float_t         gen_log_Q2;
   Float_t         obs_ptbal;
   Float_t         obs_pzbal;
   Float_t         obs_hfs_theta;
   Float_t         dnn_x;
   Float_t         dnn_Q2;
   Float_t         dnn_y;
   Float_t         dnn_tau1b;
   Float_t         obs_kine_xda;
   Float_t         obs_kine_xida;
   Float_t         obs_kine_xe;

   // List of branches
   TBranch        *b_has_isr;   //!
   TBranch        *b_has_fsr;   //!
   TBranch        *b_tower_sum_40;   //!
   TBranch        *b_n_towers_40;   //!
   TBranch        *b_eta_pho_closest_to_ebeam;   //!
   TBranch        *b_e_pho_closest_to_ebeam;   //!
   TBranch        *b_phi_pho_closest_to_ebeam;   //!
   TBranch        *b_obs_e_e;   //!
   TBranch        *b_obs_e_phi;   //!
   TBranch        *b_obs_e_theta;   //!
   TBranch        *b_obs_hfs_e;   //!
   TBranch        *b_obs_hfs_pz;   //!
   TBranch        *b_obs_hfs_pt;   //!
   TBranch        *b_obs_hfs_phi;   //!
   TBranch        *b_obs_DeltaPhi;   //!
   TBranch        *b_obs_e_trk_e;   //!
   TBranch        *b_beam_e_e;   //!
   TBranch        *b_beam_p_e;   //!
   TBranch        *b_gen_HFS_Sigma;   //!
   TBranch        *b_gen_HFS_T;   //!
   TBranch        *b_gen_eUncomb_E;   //!
   TBranch        *b_gen_eUncomb_theta;   //!
   TBranch        *b_gen_eRecomb_E;   //!
   TBranch        *b_gen_eRecomb_theta;   //!
   TBranch        *b_gen_tau1b;   //!
   TBranch        *b_obs_tau1bs;   //!
   TBranch        *b_obs_tau1be;   //!
   TBranch        *b_obs_tau1bda;   //!
   TBranch        *b_obs_cHFSs_M;   //!
   TBranch        *b_obs_cHFSs_pt;   //!
   TBranch        *b_obs_cHFSs_theta;   //!
   TBranch        *b_obs_cHFSs_phi;   //!
   TBranch        *b_obs_cHFSs_N;   //!
   TBranch        *b_obs_cHFSs_dRavg;   //!
   TBranch        *b_obs_cHFSs_dR2avg;   //!
   TBranch        *b_obs_cHFSs_Empz;   //!
   TBranch        *b_obs_cHFSs_Eppz;   //!
   TBranch        *b_obs_cHFSe_Empz;   //!
   TBranch        *b_obs_cHFSe_Eppz;   //!
   TBranch        *b_obs_cHFSe_keE0;   //!
   TBranch        *b_obs_cHFSs_keE0;   //!
   TBranch        *b_obs_cHFSs_kesE0;   //!
   TBranch        *b_wgt;   //!
   TBranch        *b_obs_e_pz;   //!
   TBranch        *b_obs_e_pt;   //!
   TBranch        *b_obs_hfs_Sigma;   //!
   TBranch        *b_obs_hfs_T;   //!
   TBranch        *b_obs_e_Sigma;   //!
   TBranch        *b_obs_e_tantheta;   //!
   TBranch        *b_obs_e_Pt2;   //!
   TBranch        *b_obs_hfs_tangamma;   //!
   TBranch        *b_obs_hfs_Empz;   //!
   TBranch        *b_obs_e_Empz;   //!
   TBranch        *b_obs_event_Empz;   //!
   TBranch        *b_rot_pt1;   //!
   TBranch        *b_rot_pt2;   //!
   TBranch        *b_rot_Empz1;   //!
   TBranch        *b_rot_Empz2;   //!
   TBranch        *b_e_ecal_over_trk_ratio;   //!
   TBranch        *b_dphi_pho_closest_to_ebeam;   //!
   TBranch        *b_has_norad;   //!
   TBranch        *b_obs_kine_ye;   //!
   TBranch        *b_obs_kine_yda;   //!
   TBranch        *b_obs_kine_yh;   //!
   TBranch        *b_obs_kine_ys;   //!
   TBranch        *b_obs_kine_Q2e;   //!
   TBranch        *b_obs_kine_Q2da;   //!
   TBranch        *b_obs_kine_Q2s;   //!
   TBranch        *b_obs_kine_Q2ida;   //!
   TBranch        *b_obs_kine_xis;   //!
   TBranch        *b_gen_e_Sigma;   //!
   TBranch        *b_gen_e_tantheta;   //!
   TBranch        *b_gen_e_Pt2;   //!
   TBranch        *b_gen_HFS_tangamma;   //!
   TBranch        *b_gen_kine_ye;   //!
   TBranch        *b_gen_kine_yda;   //!
   TBranch        *b_gen_kine_ys;   //!
   TBranch        *b_gen_kine_Q2e;   //!
   TBranch        *b_gen_kine_Q2da;   //!
   TBranch        *b_gen_kine_Q2s;   //!
   TBranch        *b_gen_kine_xis;   //!
   TBranch        *b_gen_log_x;   //!
   TBranch        *b_gen_log_y;   //!
   TBranch        *b_gen_log_Q2;   //!
   TBranch        *b_obs_ptbal;   //!
   TBranch        *b_obs_pzbal;   //!
   TBranch        *b_obs_hfs_theta;   //!
   TBranch        *b_dnn_x;   //!
   TBranch        *b_dnn_Q2;   //!
   TBranch        *b_dnn_y;   //!
   TBranch        *b_dnn_tau1b;   //!
   TBranch        *b_obs_kine_xda;   //!
   TBranch        *b_obs_kine_xida;   //!
   TBranch        *b_obs_kine_xe;   //!

   dnnout_tau1b(const char* input_file_pattern);
   virtual ~dnnout_tau1b();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop( int nbins_gen=20, int nbins_obs=35, bool verbose=false, int last_event=-1, int first_event=0, const char* out_file = "" );
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef dnnout_tau1b_cxx
dnnout_tau1b::dnnout_tau1b(const char* input_file_pattern) : fChain(0) 
{
   TChain* chain = new TChain("dnnout") ;
   int nfiles = chain -> Add( input_file_pattern ) ;
   if ( nfiles <= 0 ) { printf("\n\n *** No files match %s\n\n", input_file_pattern ) ; gSystem -> Exit(-1) ; }
   int nevts = chain -> GetEntries() ;
   if ( nevts <= 0 ) { printf("\n\n *** No events in file(s) matching %s\n\n", input_file_pattern ) ; gSystem -> Exit(-1) ; }

   printf("\n\n Found %d file(s) matching %s.  Nevts = %d\n\n", nfiles, input_file_pattern, nevts ) ;

   TTree* tree = chain;
   Init(tree);


// // if parameter tree is not specified (or zero), connect the file
// // used to generate this class and read the Tree.
//    if (tree == 0) {
//       TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("dnn-output-h1-v2-Rapgap.root");
//       if (!f || !f->IsOpen()) {
//          f = new TFile("dnn-output-h1-v2-Rapgap.root");
//       }
//       f->GetObject("dnnout",tree);

//    }
//    Init(tree);
}

dnnout_tau1b::~dnnout_tau1b()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t dnnout_tau1b::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t dnnout_tau1b::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void dnnout_tau1b::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("has_isr", &has_isr, &b_has_isr);
   fChain->SetBranchAddress("has_fsr", &has_fsr, &b_has_fsr);
   fChain->SetBranchAddress("tower_sum_40", &tower_sum_40, &b_tower_sum_40);
   fChain->SetBranchAddress("n_towers_40", &n_towers_40, &b_n_towers_40);
   fChain->SetBranchAddress("eta_pho_closest_to_ebeam", &eta_pho_closest_to_ebeam, &b_eta_pho_closest_to_ebeam);
   fChain->SetBranchAddress("e_pho_closest_to_ebeam", &e_pho_closest_to_ebeam, &b_e_pho_closest_to_ebeam);
   fChain->SetBranchAddress("phi_pho_closest_to_ebeam", &phi_pho_closest_to_ebeam, &b_phi_pho_closest_to_ebeam);
   fChain->SetBranchAddress("obs_e_e", &obs_e_e, &b_obs_e_e);
   fChain->SetBranchAddress("obs_e_phi", &obs_e_phi, &b_obs_e_phi);
   fChain->SetBranchAddress("obs_e_theta", &obs_e_theta, &b_obs_e_theta);
   fChain->SetBranchAddress("obs_hfs_e", &obs_hfs_e, &b_obs_hfs_e);
   fChain->SetBranchAddress("obs_hfs_pz", &obs_hfs_pz, &b_obs_hfs_pz);
   fChain->SetBranchAddress("obs_hfs_pt", &obs_hfs_pt, &b_obs_hfs_pt);
   fChain->SetBranchAddress("obs_hfs_phi", &obs_hfs_phi, &b_obs_hfs_phi);
   fChain->SetBranchAddress("obs_DeltaPhi", &obs_DeltaPhi, &b_obs_DeltaPhi);
   fChain->SetBranchAddress("obs_e_trk_e", &obs_e_trk_e, &b_obs_e_trk_e);
   fChain->SetBranchAddress("beam_e_e", &beam_e_e, &b_beam_e_e);
   fChain->SetBranchAddress("beam_p_e", &beam_p_e, &b_beam_p_e);
   fChain->SetBranchAddress("gen_HFS_Sigma", &gen_HFS_Sigma, &b_gen_HFS_Sigma);
   fChain->SetBranchAddress("gen_HFS_T", &gen_HFS_T, &b_gen_HFS_T);
   fChain->SetBranchAddress("gen_eUncomb_E", &gen_eUncomb_E, &b_gen_eUncomb_E);
   fChain->SetBranchAddress("gen_eUncomb_theta", &gen_eUncomb_theta, &b_gen_eUncomb_theta);
   fChain->SetBranchAddress("gen_eRecomb_E", &gen_eRecomb_E, &b_gen_eRecomb_E);
   fChain->SetBranchAddress("gen_eRecomb_theta", &gen_eRecomb_theta, &b_gen_eRecomb_theta);
   fChain->SetBranchAddress("gen_tau1b", &gen_tau1b, &b_gen_tau1b);
   fChain->SetBranchAddress("obs_tau1bs", &obs_tau1bs, &b_obs_tau1bs);
   fChain->SetBranchAddress("obs_tau1be", &obs_tau1be, &b_obs_tau1be);
   fChain->SetBranchAddress("obs_tau1bda", &obs_tau1bda, &b_obs_tau1bda);
   fChain->SetBranchAddress("obs_cHFSs_M", &obs_cHFSs_M, &b_obs_cHFSs_M);
   fChain->SetBranchAddress("obs_cHFSs_pt", &obs_cHFSs_pt, &b_obs_cHFSs_pt);
   fChain->SetBranchAddress("obs_cHFSs_theta", &obs_cHFSs_theta, &b_obs_cHFSs_theta);
   fChain->SetBranchAddress("obs_cHFSs_phi", &obs_cHFSs_phi, &b_obs_cHFSs_phi);
   fChain->SetBranchAddress("obs_cHFSs_N", &obs_cHFSs_N, &b_obs_cHFSs_N);
   fChain->SetBranchAddress("obs_cHFSs_dRavg", &obs_cHFSs_dRavg, &b_obs_cHFSs_dRavg);
   fChain->SetBranchAddress("obs_cHFSs_dR2avg", &obs_cHFSs_dR2avg, &b_obs_cHFSs_dR2avg);
   fChain->SetBranchAddress("obs_cHFSs_Empz", &obs_cHFSs_Empz, &b_obs_cHFSs_Empz);
   fChain->SetBranchAddress("obs_cHFSs_Eppz", &obs_cHFSs_Eppz, &b_obs_cHFSs_Eppz);
   fChain->SetBranchAddress("obs_cHFSe_Empz", &obs_cHFSe_Empz, &b_obs_cHFSe_Empz);
   fChain->SetBranchAddress("obs_cHFSe_Eppz", &obs_cHFSe_Eppz, &b_obs_cHFSe_Eppz);
   fChain->SetBranchAddress("obs_cHFSe_keE0", &obs_cHFSe_keE0, &b_obs_cHFSe_keE0);
   fChain->SetBranchAddress("obs_cHFSs_keE0", &obs_cHFSs_keE0, &b_obs_cHFSs_keE0);
   fChain->SetBranchAddress("obs_cHFSs_kesE0", &obs_cHFSs_kesE0, &b_obs_cHFSs_kesE0);
   fChain->SetBranchAddress("wgt", &wgt, &b_wgt);
   fChain->SetBranchAddress("obs_e_pz", &obs_e_pz, &b_obs_e_pz);
   fChain->SetBranchAddress("obs_e_pt", &obs_e_pt, &b_obs_e_pt);
   fChain->SetBranchAddress("obs_hfs_Sigma", &obs_hfs_Sigma, &b_obs_hfs_Sigma);
   fChain->SetBranchAddress("obs_hfs_T", &obs_hfs_T, &b_obs_hfs_T);
   fChain->SetBranchAddress("obs_e_Sigma", &obs_e_Sigma, &b_obs_e_Sigma);
   fChain->SetBranchAddress("obs_e_tantheta", &obs_e_tantheta, &b_obs_e_tantheta);
   fChain->SetBranchAddress("obs_e_Pt2", &obs_e_Pt2, &b_obs_e_Pt2);
   fChain->SetBranchAddress("obs_hfs_tangamma", &obs_hfs_tangamma, &b_obs_hfs_tangamma);
   fChain->SetBranchAddress("obs_hfs_Empz", &obs_hfs_Empz, &b_obs_hfs_Empz);
   fChain->SetBranchAddress("obs_e_Empz", &obs_e_Empz, &b_obs_e_Empz);
   fChain->SetBranchAddress("obs_event_Empz", &obs_event_Empz, &b_obs_event_Empz);
   fChain->SetBranchAddress("rot_pt1", &rot_pt1, &b_rot_pt1);
   fChain->SetBranchAddress("rot_pt2", &rot_pt2, &b_rot_pt2);
   fChain->SetBranchAddress("rot_Empz1", &rot_Empz1, &b_rot_Empz1);
   fChain->SetBranchAddress("rot_Empz2", &rot_Empz2, &b_rot_Empz2);
   fChain->SetBranchAddress("e_ecal_over_trk_ratio", &e_ecal_over_trk_ratio, &b_e_ecal_over_trk_ratio);
   fChain->SetBranchAddress("dphi_pho_closest_to_ebeam", &dphi_pho_closest_to_ebeam, &b_dphi_pho_closest_to_ebeam);
   fChain->SetBranchAddress("has_norad", &has_norad, &b_has_norad);
   fChain->SetBranchAddress("obs_kine_ye", &obs_kine_ye, &b_obs_kine_ye);
   fChain->SetBranchAddress("obs_kine_yda", &obs_kine_yda, &b_obs_kine_yda);
   fChain->SetBranchAddress("obs_kine_yh", &obs_kine_yh, &b_obs_kine_yh);
   fChain->SetBranchAddress("obs_kine_ys", &obs_kine_ys, &b_obs_kine_ys);
   fChain->SetBranchAddress("obs_kine_Q2e", &obs_kine_Q2e, &b_obs_kine_Q2e);
   fChain->SetBranchAddress("obs_kine_Q2da", &obs_kine_Q2da, &b_obs_kine_Q2da);
   fChain->SetBranchAddress("obs_kine_Q2s", &obs_kine_Q2s, &b_obs_kine_Q2s);
   fChain->SetBranchAddress("obs_kine_Q2ida", &obs_kine_Q2ida, &b_obs_kine_Q2ida);
   fChain->SetBranchAddress("obs_kine_xis", &obs_kine_xis, &b_obs_kine_xis);
   fChain->SetBranchAddress("gen_e_Sigma", &gen_e_Sigma, &b_gen_e_Sigma);
   fChain->SetBranchAddress("gen_e_tantheta", &gen_e_tantheta, &b_gen_e_tantheta);
   fChain->SetBranchAddress("gen_e_Pt2", &gen_e_Pt2, &b_gen_e_Pt2);
   fChain->SetBranchAddress("gen_HFS_tangamma", &gen_HFS_tangamma, &b_gen_HFS_tangamma);
   fChain->SetBranchAddress("gen_kine_ye", &gen_kine_ye, &b_gen_kine_ye);
   fChain->SetBranchAddress("gen_kine_yda", &gen_kine_yda, &b_gen_kine_yda);
   fChain->SetBranchAddress("gen_kine_ys", &gen_kine_ys, &b_gen_kine_ys);
   fChain->SetBranchAddress("gen_kine_Q2e", &gen_kine_Q2e, &b_gen_kine_Q2e);
   fChain->SetBranchAddress("gen_kine_Q2da", &gen_kine_Q2da, &b_gen_kine_Q2da);
   fChain->SetBranchAddress("gen_kine_Q2s", &gen_kine_Q2s, &b_gen_kine_Q2s);
   fChain->SetBranchAddress("gen_kine_xis", &gen_kine_xis, &b_gen_kine_xis);
   fChain->SetBranchAddress("gen_log_x", &gen_log_x, &b_gen_log_x);
   fChain->SetBranchAddress("gen_log_y", &gen_log_y, &b_gen_log_y);
   fChain->SetBranchAddress("gen_log_Q2", &gen_log_Q2, &b_gen_log_Q2);
   fChain->SetBranchAddress("obs_ptbal", &obs_ptbal, &b_obs_ptbal);
   fChain->SetBranchAddress("obs_pzbal", &obs_pzbal, &b_obs_pzbal);
   fChain->SetBranchAddress("obs_hfs_theta", &obs_hfs_theta, &b_obs_hfs_theta);
   fChain->SetBranchAddress("dnn_x", &dnn_x, &b_dnn_x);
   fChain->SetBranchAddress("dnn_Q2", &dnn_Q2, &b_dnn_Q2);
   fChain->SetBranchAddress("dnn_y", &dnn_y, &b_dnn_y);
   fChain->SetBranchAddress("dnn_tau1b", &dnn_tau1b, &b_dnn_tau1b);
   fChain->SetBranchAddress("obs_kine_xda", &obs_kine_xda, &b_obs_kine_xda);
   fChain->SetBranchAddress("obs_kine_xida", &obs_kine_xida, &b_obs_kine_xida);
   fChain->SetBranchAddress("obs_kine_xe", &obs_kine_xe, &b_obs_kine_xe);
   Notify();
}

Bool_t dnnout_tau1b::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void dnnout_tau1b::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t dnnout_tau1b::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef dnnout_tau1b_cxx
