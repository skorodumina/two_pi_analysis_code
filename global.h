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

extern Float_t LiveTime,inclusive,elastic,P_EL,th_EL,z_EL,z_P,z_PIp,z_PIm,ph_EL,ECT,nphe,theta_cc,ph_cc;
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

extern TH1F *h_mixed_prod_excl;
extern TH1F *h_mixed_prod_pimmiss;

extern TH1F *h_mixed_prod_excl_sim;
extern TH1F *h_mixed_prod_pimmiss_sim;

extern TH1F *h_pim_mis_all_reg_th_dep_pim[5][12];
extern TH1F *h_pim_mis_all_reg_th_dep_pim_sim[5][12];

extern TH1F *h_mom_all_reg_th_dep_pim[5][12];
extern TH1F *h_mom_all_reg_th_dep_pim_sim[5][12];

extern TH1F *h_pim_mis_all_reg_th_dep_pip[5][12];
extern TH1F *h_pim_mis_all_reg_th_dep_pip_sim[5][12];

extern TH1F *h_mom_all_reg_th_dep_pip[5][12];
extern TH1F *h_mom_all_reg_th_dep_pip_sim[5][12];

extern TH1F *h_pim_mis_all_reg_th_dep_pr[5][12];
extern TH1F *h_pim_mis_all_reg_th_dep_pr_sim[5][12];

extern TH1F *h_pim_mis_th_dep_pr[5][12];
extern TH1F *h_pim_mis_th_dep_pr_sim[5][12];

extern TH1F *h_mom_all_reg_th_dep_pr[5][12];
extern TH1F *h_mom_all_reg_th_dep_pr_sim[5][12];


extern TH1F *h_pim_mis_th_dep[5][12][8];
extern TH1F *h_pim_mis_th_dep_sim[5][12][8];


extern TH1F *h_pim_mis_all_reg[5];
extern TH1F *h_pim_mis_all_reg_sim[5];

extern TH1F *h_pip_mis_all_reg[5];
extern TH1F *h_pip_mis_all_reg_sim[5];

extern TH1F *h_pr_mis_all_reg[5];
extern TH1F *h_pr_mis_all_reg_sim[5];


extern TH1F *h_mom_all_reg[5];
extern TH1F *h_mom_all_reg_sim[5];


extern TH1F *h_test;
extern TH1F *h_test_sim;
extern TH2F *h_test_2d_1,*h_test_2d_2,*h_test_2d_3,*h_test_2d_4,*h_test_2d_5,*h_test_2d_6,*h_test_2d_7,*h_test_2d_8,*h_test_2d_9;

extern TH2F *h_test_2d_1_sim,*h_test_2d_2_sim,*h_test_2d_3_sim,*h_test_2d_4_sim,*h_test_2d_5_sim,*h_test_2d_6_sim,*h_test_2d_7_sim,*h_test_2d_8_sim,*h_test_2d_9_sim;

extern TH1F *h_z_P;
extern TH1F *h_z_PIp;
extern TH1F *h_z_PIm;

extern TH1F *h_z_P_sim;
extern TH1F *h_z_PIp_sim;
extern TH1F *h_z_PIm_sim;
extern TH1F *h_pim_mis_fermi_nocut_1[4][21];
extern TH1F *h_pim_mis_fermi_nocut_sim_1[4][21];
extern TH1F *h_pim_mis_fermi_momcut_1[4][5];
extern TH1F *h_pim_mis_fermi_momcut_sim_1[4][5];
extern TH1F *h_mis_mom_fermi_1[4][5];
extern TH1F *h_mis_mom_fermi_sim_1[4][5];
extern TH1F *h_mis_mom_fermi_mmas_cut_1[4][5];
extern TH1F *h_mis_mom_fermi_mmas_cut_sim_1[4][5];



extern TH1F *h_pim_mis_fermi_nocut_2[2][19];
extern TH1F *h_pim_mis_fermi_nocut_sim_2[2][19];
extern TH1F *h_pim_mis_fermi_momcut_2[2][5];
extern TH1F *h_pim_mis_fermi_momcut_sim_2[2][5];
extern TH1F *h_mis_mom_fermi_2[2][5];
extern TH1F *h_mis_mom_fermi_sim_2[2][5];
extern TH1F *h_mis_mom_fermi_mmas_cut_2[2][5];
extern TH1F *h_mis_mom_fermi_mmas_cut_sim_2[2][5];


extern TH1F *h_pim_mis_fermi_nocut_3[2][17];
extern TH1F *h_pim_mis_fermi_nocut_sim_3[2][17];
extern TH1F *h_pim_mis_fermi_momcut_3[2][4];
extern TH1F *h_pim_mis_fermi_momcut_sim_3[2][4];
extern TH1F *h_mis_mom_fermi_3[2][4];
extern TH1F *h_mis_mom_fermi_sim_3[2][4];
extern TH1F *h_mis_mom_fermi_mmas_cut_3[2][4];
extern TH1F *h_mis_mom_fermi_mmas_cut_sim_3[2][4];


extern TH1F *h_pim_mis_fermi_nocut_4[14];
extern TH1F *h_pim_mis_fermi_nocut_sim_4[14];
extern TH1F *h_pim_mis_fermi_momcut_4[4];
extern TH1F *h_pim_mis_fermi_momcut_sim_4[4];
extern TH1F *h_mis_mom_fermi_4[4];
extern TH1F *h_mis_mom_fermi_sim_4[4];
extern TH1F *h_mis_mom_fermi_mmas_cut_4[4];
extern TH1F *h_mis_mom_fermi_mmas_cut_sim_4[4];

extern TH1F *h_pim_mis_fermi_nocut_5[13];
extern TH1F *h_pim_mis_fermi_nocut_sim_5[13];
extern TH1F *h_pim_mis_fermi_momcut_5[3];
extern TH1F *h_pim_mis_fermi_momcut_sim_5[3];
extern TH1F *h_mis_mom_fermi_5[3];
extern TH1F *h_mis_mom_fermi_sim_5[3];
extern TH1F *h_mis_mom_fermi_mmas_cut_5[3];
extern TH1F *h_mis_mom_fermi_mmas_cut_sim_5[3];

extern TH1F *h_pim_mis_fermi_nocut_6[12];
extern TH1F *h_pim_mis_fermi_nocut_sim_6[12];
extern TH1F *h_pim_mis_fermi_momcut_6[3];
extern TH1F *h_pim_mis_fermi_momcut_sim_6[3];
extern TH1F *h_mis_mom_fermi_6[3];
extern TH1F *h_mis_mom_fermi_sim_6[3];
extern TH1F *h_mis_mom_fermi_mmas_cut_6[3];
extern TH1F *h_mis_mom_fermi_mmas_cut_sim_6[3];


extern TH1F *h_pim_mis_fermi_nocut_7[10];
extern TH1F *h_pim_mis_fermi_nocut_sim_7[10];
extern TH1F *h_pim_mis_fermi_momcut_7[3];
extern TH1F *h_pim_mis_fermi_momcut_sim_7[3];
extern TH1F *h_mis_mom_fermi_7[3];
extern TH1F *h_mis_mom_fermi_sim_7[3];
extern TH1F *h_mis_mom_fermi_mmas_cut_7[3];
extern TH1F *h_mis_mom_fermi_mmas_cut_sim_7[3];

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


extern TH1F *h_1d_rc_0425_evt;
extern TH1F *h_1d_rc_0475_evt;
extern TH1F *h_1d_rc_0525_evt;
extern TH1F *h_1d_rc_0575_evt;
extern TH1F *h_1d_rc_0625_evt;
extern TH1F *h_1d_rc_0675_evt;
extern TH1F *h_1d_rc_0725_evt;
extern TH1F *h_1d_rc_0775_evt;
extern TH1F *h_1d_rc_0825_evt;
extern TH1F *h_1d_rc_0875_evt;
extern TH1F *h_1d_rc_0925_evt;
extern TH1F *h_1d_rc_0975_evt;



extern TH2F *hist_ectot_sector1_sim,*hist_ectot_sector2_sim,*hist_ectot_sector3_sim,*hist_ectot_sector4_sim,*hist_ectot_sector5_sim,*hist_ectot_sector6_sim;
extern TH2F *hist_ectot_sector1,*hist_ectot_sector2,*hist_ectot_sector3,*hist_ectot_sector4,*hist_ectot_sector5,*hist_ectot_sector6;
extern TH2F *W_2pi_selection,*W_2pi_selection_sim;

extern TH1F *hist_PIm_miss,*hist_PIm_miss_en,*hist_PIp_miss,*hist_PIp_miss_d,*hist_PIp_miss_d_bef, *hist_PIp_miss_en,*hist_P_miss,*hist_P_miss_en, *h_miss_mass_0;
extern TH1F *h_miss_mass_0_d, *h_miss_mass_0_d_sim;
extern TH1F *h_miss_mass_0_d_mmcut, *h_miss_mass_0_d_sim_mmcut;
extern TH1F *h_inprot_miss,*h_inprot_miss_en;
extern TH1F *h_inprot_miss_sim,*h_inprot_miss_en_sim;

extern TH1F *hist_PIm_miss_all_reg_1, *hist_PIm_miss_all_reg_2;
extern TH1F *hist_PIm_miss_all_reg_1_sim, *hist_PIm_miss_all_reg_2_sim;

extern TH1F *hist_miss_en_0,*h_miss_mom_0,*h_miss_mom_0_sim,*h_miss_mom_0_d,*h_miss_mom_0_d_sim; 
extern TH1F *h_miss_mom_0_d_mmcut,*h_miss_mom_0_d_sim_mmcut; 
extern TH1F *h_miss_mom_0_nocut, *h_miss_mom_0_cut_on0, *h_miss_mom_0_cut_onpim, *h_miss_mom_0_nocut_sim, *h_miss_mom_0_cut_on0_sim, *h_miss_mom_0_cut_onpim_sim;

extern TH1F *h_PIp_miss_d_sim,*h_PIp_miss_d_bef_sim,*h_PIm_miss_sim;
extern TH1F *h_PIp_miss_sim, *h_miss_mass_0_sim, *h_miss_en_0_sim, *h_PIm_miss_en_sim, *h_PIp_miss_en_sim;
extern TH1F *hist_w_hadr_all_reg, *hist_w_el_all_reg, *hist_w_sim_new_1dim,*hist_w_sim_old_1dim ;
extern Float_t W,Q2;
extern TH2F *ph_vs_th_1,*ph_vs_th_2,*ph_vs_th_3,*ph_vs_th_4,*ph_vs_th_5,*ph_vs_th_6;
extern TH2F *ph_vs_th_1pe[17],*ph_vs_th_2pe[17],*ph_vs_th_3pe[17],*ph_vs_th_4pe[17],*ph_vs_th_5pe[17],*ph_vs_th_6pe[17];
extern TH2F *h_mm_0_vs_npart, *h_mm_pim_vs_npart, *h_mm_pip_vs_npart;
extern TH2F *h_mm_0_vs_npart_sim, *h_mm_pim_vs_npart_sim, *h_mm_pip_vs_npart_sim;

extern TH2F *ph_vs_th_el_sim[6][7];
extern TH2F *ph_vs_th_p_sim[6];
extern TH2F *ph_vs_th_pip_sim[6];
extern TH2F *ph_vs_th_pim_sim[6][5];

extern TH2F  *h_cc_nphe_total_s1,*h_cc_nphe_total_s2,*h_cc_nphe_total_s3,*h_cc_nphe_total_s4,*h_cc_nphe_total_s5,*h_cc_nphe_total_s6;

extern TH2F  *h_cc_nphe_final_s1,*h_cc_nphe_final_s2,*h_cc_nphe_final_s3,*h_cc_nphe_final_s4,*h_cc_nphe_final_s5,*h_cc_nphe_final_s6;


extern TH2F *ph_vs_th_1pe_fid[17],*ph_vs_th_2pe_fid[17],*ph_vs_th_3pe_fid[17],*ph_vs_th_4pe_fid[17],*ph_vs_th_5pe_fid[17],*ph_vs_th_6pe_fid[17];

extern TH2F *hist_sector1,*hist_sector2,*hist_sector3,*hist_sector4,*hist_sector5,*hist_sector6;
extern TH2F *hist_nphe_sector1,*hist_nphe_sector2,*hist_nphe_sector3,*hist_nphe_sector4,*hist_nphe_sector5,*hist_nphe_sector6;
extern TH1F *nphe_sector1,*nphe_sector2,*nphe_sector3,*nphe_sector4,*nphe_sector5,*nphe_sector6;
extern TH1F *nphe_sector1_after,*nphe_sector2_after,*nphe_sector3_after,*nphe_sector4_after,*nphe_sector5_after,*nphe_sector6_after;
extern TH2F  *ph_vs_th_p_1,*ph_vs_th_p_2,*ph_vs_th_p_3,*ph_vs_th_p_4,*ph_vs_th_p_5,*ph_vs_th_p_6;
extern TH2F  *ph_vs_th_p_1_w,*ph_vs_th_p_2_w,*ph_vs_th_p_3_w,*ph_vs_th_p_4_w,*ph_vs_th_p_5_w,*ph_vs_th_p_6_w;
extern TH2F  *ph_th_p_1[5],*ph_th_p_2[5],*ph_th_p_3[5],*ph_th_p_4[5],*ph_th_p_5[5],*ph_th_p_6[5];
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
extern TH2F  *h_delta_w_vs_w_old_data;
extern TH2F  *h_delta_w_vs_w_old_data_mmcut;
extern TH2F  *h_delta_w_vs_w_old_gen;
extern TH2F  *h_delta_w_vs_w_old_rec;

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
  
  
  extern TH2F *ph_vs_th_pim[6][10];

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
