#ifndef GLOBAL_H
#include "TROOT.h"
#include "TFile.h"
#include "TLine.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TMacro.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"
#include "TMath.h"
#include <math.h>
#include <TLorentzVector.h>
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TText.h"
#include "TStyle.h"
#include "TGObject.h"
#include "TObject.h"
#include "TSystem.h"
#include "TMinuit.h"
#include <TRint.h>
#include <stdio.h>
#include <dlfcn.h>
#define GLOBAL_H

#endif

extern Int_t npart,segment,sector,indtype,n_P,n_PIp,n_PIm;
extern Short_t pmt_hit;
extern Int_t PdHit_EL,PdHit_PIp,PdHit_P,PdHit_PIm;

extern Float_t m_proton,m_pip,beta;
 extern   Float_t delta_mom_p_skor; 
extern Float_t LiveTime,inclusive,elastic,P_EL,th_EL,z_EL,z_P,z_PIp,z_PIm,ph_EL,ECT,nphe,theta_cc,ph_cc;

extern Float_t sigma;
extern Float_t ph_P,th_P,P_P,beta_P,beta_PIp,beta_PIm;

extern Float_t P_time,PIp_time,PIm_time;
extern Float_t P_dist,PIp_dist,PIm_dist;
extern Float_t P_PIp,ph_PIp,th_PIp,Nphe_PIp;
extern Float_t P_PIm,ph_PIm,th_PIm;
extern Float_t px_fermi,py_fermi,pz_fermi;
extern TLorentzVector  P4_inprot_miss;



extern Float_t sc_x,sc_y,sc_x_p,sc_y_p,sc_x_pip,sc_y_pip,sc_x_pim,sc_y_pim;

extern Double_t theta_PIm_cm,theta_PIp_cm,theta_P_cm,phi_PIm_cm,phi_P_cm,phi_PIp_cm,alpha_PPIp_piPIm, alpha_PIpPIm_pipf,alpha_PPIm_piPIp; 
extern TLorentzVector P4_EL,P4_ELP_reg,P4_PP_reg,P4_PIp_reg,P4_PIm_reg,P4_P;
extern Double_t inv_m_pip_pim,inv_m_pip_p,inv_m_pim_p;

extern Float_t fract_integ[3][6][18];


extern Short_t ph_cc_match;
extern TH1F *h_sim_mom_corr_test;

extern TH1F *h_pim_mis_main_top[5];
extern TH1F *h_pim_mis_main_top_sim[5];

extern TH1F *h_pim_mis_all_reg[5];
extern TH1F *h_pim_mis_all_reg_sim[5];

extern TH1F *h_0_mis_all_reg[5];
extern TH1F *h_0_mis_all_reg_sim[5];

extern TH1F *h_mom_all_reg[5];
extern TH1F *h_mom_all_reg_sim[5];


extern TH1F *h_1d_rc_0425;
extern TH1F *h_1d_rc_0475;
extern TH1F *h_1d_rc_0525;
extern TH1F *h_1d_rc_0575;
extern TH1F *h_1d_rc_0625;
extern TH1F *h_1d_rc_0675;
extern TH1F *h_1d_rc_0725;
extern TH1F *h_1d_rc_0775;
extern TH1F *h_1d_rc_0825;
extern TH1F *h_1d_rc_0875;
extern TH1F *h_1d_rc_0925;
extern TH1F *h_1d_rc_0975;

extern TH1F *h_z_P;
extern TH1F *h_z_PIp;
extern TH1F *h_z_PIm;

extern TH1F *h_z_P_sim;
extern TH1F *h_z_PIp_sim;
extern TH1F *h_z_PIm_sim;

extern TH1F *h_inv_NP[5];
extern TH1F *h_inv_NPIp[5];
extern TH1F *h_inv_NPIm[5];




extern TH2F *hist_ectot_sector1_sim,*hist_ectot_sector2_sim,*hist_ectot_sector3_sim,*hist_ectot_sector4_sim,*hist_ectot_sector5_sim,*hist_ectot_sector6_sim;

extern TH2F *th_cc_vs_seg_1_sim,*th_cc_vs_seg_2_sim,*th_cc_vs_seg_3_sim,*th_cc_vs_seg_4_sim,*th_cc_vs_seg_5_sim,*th_cc_vs_seg_6_sim;

extern TH2F *hist_ectot_sector1,*hist_ectot_sector2,*hist_ectot_sector3,*hist_ectot_sector4,*hist_ectot_sector5,*hist_ectot_sector6;

extern TH2F *th_cc_vs_seg_1,*th_cc_vs_seg_2,*th_cc_vs_seg_3,*th_cc_vs_seg_4,*th_cc_vs_seg_5,*th_cc_vs_seg_6;

extern TH2F *W_2pi_selection,*W_2pi_selection_sim;

extern TH1F *hist_PIm_miss,*hist_PIm_miss_en,*hist_PIp_miss,*hist_PIp_miss_d,*hist_PIp_miss_d_bef, *hist_PIp_miss_en,*hist_P_miss,*hist_P_miss_en, *h_miss_mass_0;


extern TH1F *hist_PIm_miss_all_reg_1, *hist_PIm_miss_all_reg_2;
extern TH1F *hist_PIm_miss_all_reg_1_sim, *hist_PIm_miss_all_reg_2_sim;

extern TH1F *hist_miss_en_0,*h_miss_mom_0,*h_miss_mom_0_sim,*h_miss_mom_0_d; 

extern TH1F *h_miss_mom_0_nocut, *h_miss_mom_0_cut_on0, *h_miss_mom_0_cut_onpim, *h_miss_mom_0_nocut_sim, *h_miss_mom_0_cut_on0_sim, *h_miss_mom_0_cut_onpim_sim;

extern TH1F *h_PIp_miss_d_sim,*h_PIp_miss_d_bef_sim,*h_PIm_miss_sim;
extern TH1F *h_PIp_miss_sim, *h_miss_mass_0_sim, *h_miss_en_0_sim, *h_PIm_miss_en_sim, *h_PIp_miss_en_sim;
extern TH1F *hist_w_hadr_all_reg, *hist_w_el_all_reg;
extern Float_t W,Q2;
extern TH2F *ph_vs_th_1,*ph_vs_th_2,*ph_vs_th_3,*ph_vs_th_4,*ph_vs_th_5,*ph_vs_th_6;
extern TH2F *ph_vs_th_1pe[17],*ph_vs_th_2pe[17],*ph_vs_th_3pe[17],*ph_vs_th_4pe[17],*ph_vs_th_5pe[17],*ph_vs_th_6pe[17];
extern TH2F *h_mm_0_vs_npart, *h_mm_pim_vs_npart, *h_mm_pip_vs_npart;
extern TH2F *h_mm_0_vs_npart_sim, *h_mm_pim_vs_npart_sim, *h_mm_pip_vs_npart_sim;

extern TH2F *ph_vs_th_el_sim[6][7];
extern TH2F *ph_vs_th_p_sim[6];
extern TH2F *ph_vs_th_pip_sim[6];
extern TH2F *ph_vs_th_pim_sim[6][5];

extern TH2F *h_dc_y_vs_x_el;
extern TH1F *h_z_corr1_data, *h_z_corr2_data,*h_z_corr1_sim, *h_z_corr2_sim;

extern TH2F  *h_cc_nphe_total_s1,*h_cc_nphe_total_s2,*h_cc_nphe_total_s3,*h_cc_nphe_total_s4,*h_cc_nphe_total_s5,*h_cc_nphe_total_s6;

extern TH2F  *h_cc_nphe_final_s1,*h_cc_nphe_final_s2,*h_cc_nphe_final_s3,*h_cc_nphe_final_s4,*h_cc_nphe_final_s5,*h_cc_nphe_final_s6;


extern TH2F *ph_vs_th_1pe_fid[17],*ph_vs_th_2pe_fid[17],*ph_vs_th_3pe_fid[17],*ph_vs_th_4pe_fid[17],*ph_vs_th_5pe_fid[17],*ph_vs_th_6pe_fid[17];

extern TH2F *hist_sector1,*hist_sector2,*hist_sector3,*hist_sector4,*hist_sector5,*hist_sector6;
extern TH2F *hist_nphe_sector1,*hist_nphe_sector2,*hist_nphe_sector3,*hist_nphe_sector4,*hist_nphe_sector5,*hist_nphe_sector6;
extern TH1F *nphe_sector1,*nphe_sector2,*nphe_sector3,*nphe_sector4,*nphe_sector5,*nphe_sector6;
extern TH1F *nphe_sector1_after,*nphe_sector2_after,*nphe_sector3_after,*nphe_sector4_after,*nphe_sector5_after,*nphe_sector6_after;
extern TH2F  *ph_vs_th_p_1,*ph_vs_th_p_2,*ph_vs_th_p_3,*ph_vs_th_p_4,*ph_vs_th_p_5,*ph_vs_th_p_6;
extern TH2F  *ph_vs_th_p_1_w,*ph_vs_th_p_2_w,*ph_vs_th_p_3_w,*ph_vs_th_p_4_w,*ph_vs_th_p_5_w,*ph_vs_th_p_6_w;
extern TH2F  *ph_th_p_1[6],*ph_th_p_2[6],*ph_th_p_3[6],*ph_th_p_4[6],*ph_th_p_5[6],*ph_th_p_6[6];
extern TH1F *hist_z_el_1, *hist_z_el_2, *hist_z_el_3, *hist_z_el_4, *hist_z_el_5, *hist_z_el_6;
extern TH1F *hist_z_el_1_empty, *hist_z_el_2_empty, *hist_z_el_3_empty, *hist_z_el_4_empty, *hist_z_el_5_empty, *hist_z_el_6_empty;
extern TH1F *hist_z_el_1_sim_1, *hist_z_el_2_sim_1, *hist_z_el_3_sim_1, *hist_z_el_4_sim_1, *hist_z_el_5_sim_1, *hist_z_el_6_sim_1;
extern TH1F *hist_z_el_1_sim_2, *hist_z_el_2_sim_2, *hist_z_el_3_sim_2, *hist_z_el_4_sim_2, *hist_z_el_5_sim_2, *hist_z_el_6_sim_2;



extern TH2F  *ph_th_pip_1[6],*ph_th_pip_2[6],*ph_th_pip_3[6],*ph_th_pip_4[6],*ph_th_pip_5[6],*ph_th_pip_6[6];

extern TH2F  *ph_vs_th_pip_1,*ph_vs_th_pip_2,*ph_vs_th_pip_3,*ph_vs_th_pip_4,*ph_vs_th_pip_5,*ph_vs_th_pip_6;

extern TH2F  *norm_nphe_s1,*norm_nphe_s2,*norm_nphe_s3,*norm_nphe_s4,*norm_nphe_s5,*norm_nphe_s6;

extern TH2F  *avrg_nphe_sector1,*avrg_nphe_sector2,*avrg_nphe_sector3,*avrg_nphe_sector4,*avrg_nphe_sector5,*avrg_nphe_sector6;
extern TFile *outFile;
extern TH1F *hist_ltime,*hist_ltime_1d,*hist_n_incl,*hist_n_incl_1d,*hist_n_elast,*hist_n_elast_1d;
 
extern TH2F  *W_2pi_fid_p;

 extern TH2F *beta_vs_p_p[6][48];
 extern TH2F *beta_vs_p_pim[6][48];
  extern TH2F *beta_vs_p_pip[6][48];
  
  extern TH2F *beta_vs_p_p_sim[6][48];
 extern TH2F *beta_vs_p_pim_sim[6][48];
  extern TH2F *beta_vs_p_pip_sim[6][48]; 
  
  
  extern TH2F *time_pip[6][48];
  extern TH2F *time_pim[6][48];
  extern TH2F *time_p[6][48];
  
  extern TH2F *time_pip_sim[6][48];
  extern TH2F *time_pim_sim[6][48];
  extern TH2F *time_p_sim[6][48];  
  
  
  extern TH2F *ph_vs_th_pim[6][15];

  //photoelectrons
  extern TH1F *ph_el_left[6][20];
  extern TH1F *ph_el_both[6][20];
  extern TH1F *ph_el_right[6][20];
  
  extern TH2F *ph_th_pim_all_p[6];

extern TH1F *h_inv_m_pip_pim,*h_inv_m_pip_p,*h_inv_m_pim_p; 
extern TH1F *h_inv_m_pip_pim_sim,*h_inv_m_pip_p_sim,*h_inv_m_pim_p_sim; 

extern TH1F *h_inv_m_pip_pim_bin, *h_inv_m_pip_p_bin, *h_inv_m_pim_p_bin, *h_theta_PIm_cm, *h_theta_PIp_cm, *h_theta_P_cm, *h_phi_PIm_cm, *h_phi_PIp_cm, *h_phi_P_cm, *h_alpha_PIpPIm_pipf, *h_alpha_PPIp_piPIm, *h_alpha_PPIm_piPIp;




extern THnSparseD *h_5dim_excl_1[12][21];
extern THnSparseD *h_5dim_excl_2[12][21];
extern THnSparseD *h_5dim_excl_3[12][21];


extern THnSparseD *h_5dim_excl_1_sim[12][21];
extern THnSparseD *h_5dim_excl_2_sim[12][21];
extern THnSparseD *h_5dim_excl_3_sim[12][21];

extern THnSparseD *h_5dim_excl_1_sim_evt[12][21];
extern THnSparseD *h_5dim_excl_2_sim_evt[12][21];
extern THnSparseD *h_5dim_excl_3_sim_evt[12][21];

extern THnSparseD *h_5dim_pim_1[12][21];
extern THnSparseD *h_5dim_pim_2[12][21];
extern THnSparseD *h_5dim_pim_3[12][21];

extern THnSparseD *h_5dim_pim_1_sim[12][21];
extern THnSparseD *h_5dim_pim_2_sim[12][21];
extern THnSparseD *h_5dim_pim_3_sim[12][21];

extern THnSparseD *h_5dim_pim_1_sim_evt[12][21];
extern THnSparseD *h_5dim_pim_2_sim_evt[12][21];
extern THnSparseD *h_5dim_pim_3_sim_evt[12][21];

extern THnSparseD *h_5dim_1_sim_gen[12][21];
extern THnSparseD *h_5dim_2_sim_gen[12][21];
extern THnSparseD *h_5dim_3_sim_gen[12][21];

extern THnSparseD *h_5dim_1_sim_gen_evt[12][21];
extern THnSparseD *h_5dim_2_sim_gen_evt[12][21];
extern THnSparseD *h_5dim_3_sim_gen_evt[12][21];




extern TH1D *h_cos_th;

extern TH1D *h_w_int[12];

extern TH2F *th_vs_p_e_1[6], *th_vs_p_e_2[6], *th_vs_p_p_1[6],*th_vs_p_p_2[6],*th_vs_p_pip_1[6],*th_vs_p_pip_2[6],*th_vs_p_pim_1[6],*th_vs_p_pim_2[6];

extern TH2F *th_vs_p_e_1_sim[6], *th_vs_p_e_2_sim[6], *th_vs_p_p_1_sim[6],*th_vs_p_p_2_sim[6],*th_vs_p_pip_1_sim[6],*th_vs_p_pip_2_sim[6],*th_vs_p_pim_1_sim[6],*th_vs_p_pim_2_sim[6];

int global();
