#include "TROOT.h"
#include "TFile.h"
#include "TLine.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TMacro.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"
#include "TMath.h"
#include <math.h>
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TText.h"
#include "TStyle.h"
#include "TGObject.h"
#include "TObject.h"
#include "TSystem.h"
#include "TMinuit.h"
#include <TRint.h>
#include <stdio.h>
#include <dlfcn.h>
#include "global.h"
#include <stdio.h>
#include <dlfcn.h>
#include <sstream>
#include <TLorentzVector.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring> 
#include <cstdlib>

 using namespace std; 
 
Int_t npart,segment,sector,indtype,n_P,n_PIp,n_PIm;
Short_t pmt_hit;
Int_t PdHit_EL,PdHit_PIp,PdHit_P,PdHit_PIm; 
 
 Float_t m_proton,m_pip,beta;
 
 Float_t LiveTime,inclusive,elastic,P_EL,th_EL,z_EL,z_P,z_PIp,z_PIm,ph_EL,ECT,nphe,theta_cc,ph_cc;
 Float_t sigma;

Float_t W,Q2;
Float_t ph_P,th_P,P_P,beta_P,beta_PIp,beta_PIm;
//Float_t beta_P_time,beta_PIp_time,beta_PIm_time;
Float_t P_time,PIp_time,PIm_time;
Float_t P_dist,PIp_dist,PIm_dist;
Float_t P_PIp,ph_PIp,th_PIp,Nphe_PIp;
Float_t P_PIm,ph_PIm,th_PIm;
 Float_t px_fermi,py_fermi,pz_fermi;
TLorentzVector  P4_inprot_miss;

Float_t sc_x,sc_y,sc_x_p,sc_y_p,sc_x_pip,sc_y_pip,sc_x_pim,sc_y_pim;
Double_t theta_PIm_cm,theta_PIp_cm,theta_P_cm,phi_PIm_cm,phi_P_cm,phi_PIp_cm,alpha_PPIp_piPIm, alpha_PIpPIm_pipf,alpha_PPIm_piPIp; 
TLorentzVector P4_EL,P4_ELP_reg,P4_PP_reg,P4_PIp_reg,P4_PIm_reg,P4_P;
Double_t inv_m_pip_pim,inv_m_pip_p,inv_m_pim_p; 

Float_t fract_integ[3][6][18]; 
Short_t ph_cc_match;
//booking electron histograms


/*TH2F *hist_sector1 = new TH2F("cc_sector1","cc_sector1",200, -5., 55.,200., -25., 25.);
TH2F *hist_sector2 = new TH2F("cc_sector2","cc_sector2",200, -5., 55.,200., -25., 25.);
TH2F *hist_sector3 = new TH2F("cc_sector3","cc_sector3",200, -5., 55.,200., -25., 25.);
TH2F *hist_sector4 = new TH2F("cc_sector4","cc_sector4",200, -5., 55.,200., -25., 25.);
TH2F *hist_sector5 = new TH2F("cc_sector5","cc_sector5",200, -5., 55.,200., -25., 25.);
TH2F *hist_sector6 = new TH2F("cc_sector6","cc_sector6",200, -5., 55.,200., -25., 25.);*/

TH1F *nphe_sector1 = new TH1F("nphe_sector1","nphe_sector1",501, -1., 500.);
TH1F *nphe_sector1_after = new TH1F("nphe_sector1_after","nphe_sector1_after",501, -1., 500.);
TH1F *nphe_sector2 = new TH1F("nphe_sector2","nphe_sector2",501, -1., 500.);
TH1F *nphe_sector2_after = new TH1F("nphe_sector2_after","nphe_sector2_after",501, -1., 500.);
TH1F *nphe_sector3 = new TH1F("nphe_sector3","nphe_sector3",501, -1., 500.);
TH1F *nphe_sector3_after = new TH1F("nphe_sector3_after","nphe_sector3_after",501, -1., 500.);
TH1F *nphe_sector4 = new TH1F("nphe_sector4","nphe_sector4",501, -1., 500.);
TH1F *nphe_sector4_after = new TH1F("nphe_sector4_after","nphe_sector4_after",501, -1., 500.);
TH1F *nphe_sector5 = new TH1F("nphe_sector5","nphe_sector5",501, -1., 500.);
TH1F *nphe_sector5_after = new TH1F("nphe_sector5_after","nphe_sector5_after",501, -1., 500.);
TH1F *nphe_sector6 = new TH1F("nphe_sector6","nphe_sector6",501, -1., 500.);
TH1F *nphe_sector6_after = new TH1F("nphe_sector6_after","nphe_sector6_after",501, -1., 500.);





/*TH2F *hist_nphe_sector1 = new TH2F("cc_nphe_sector1","cc_nphe_sector1",200, -5., 55.,200., -25., 25.);
TH2F *hist_nphe_sector2 = new TH2F("cc_nphe_sector2","cc_nphe_sector2",200, -5., 55.,200., -25., 25.);
TH2F *hist_nphe_sector3 = new TH2F("cc_nphe_sector3","cc_nphe_sector3",200, -5., 55.,200., -25., 25.);
TH2F *hist_nphe_sector4 = new TH2F("cc_nphe_sector4","cc_nphe_sector4",200, -5., 55.,200., -25., 25.);
TH2F *hist_nphe_sector5 = new TH2F("cc_nphe_sector5","cc_nphe_sector5",200, -5., 55.,200., -25., 25.);
TH2F *hist_nphe_sector6 = new TH2F("cc_nphe_sector6","cc_nphe_sector6",200, -5., 55.,200., -25., 25.);*/

TH2F  *avrg_nphe_sector1,*avrg_nphe_sector2,*avrg_nphe_sector3,*avrg_nphe_sector4,*avrg_nphe_sector5,*avrg_nphe_sector6;

TH2F  *norm_nphe_s1,*norm_nphe_s2,*norm_nphe_s3,*norm_nphe_s4,*norm_nphe_s5,*norm_nphe_s6;



//TH1F *h_pim_mis_fermi_nocut[22];
//TH1F *h_pim_mis_fermi_nocut_sim[22];
//TH1F *h_pim_mis_fermi_momcut[22];
//TH1F *h_pim_mis_fermi_momcut_sim[22];
//TH1F *h_mis_mom_fermi[22];
//TH1F *h_mis_mom_fermi_sim[22];
//TH1F *h_mis_mom_fermi_mmas_cut[22];
//TH1F *h_mis_mom_fermi_mmas_cut_sim[22];

TH1F *h_inv_NP[5];
TH1F *h_inv_NPIp[5];
TH1F *h_inv_NPIm[5];

/*TH1F *h_pim_mis_fermi_nocut_1[4][21];
TH1F *h_pim_mis_fermi_nocut_sim_1[4][21];
TH1F *h_pim_mis_fermi_momcut_1[4][5];
TH1F *h_pim_mis_fermi_momcut_sim_1[4][5];
TH1F *h_mis_mom_fermi_1[4][5];
TH1F *h_mis_mom_fermi_sim_1[4][5];
TH1F *h_mis_mom_fermi_mmas_cut_1[4][5];
TH1F *h_mis_mom_fermi_mmas_cut_sim_1[4][5];



TH1F *h_pim_mis_fermi_nocut_2[2][19];
TH1F *h_pim_mis_fermi_nocut_sim_2[2][19];
TH1F *h_pim_mis_fermi_momcut_2[2][5];
TH1F *h_pim_mis_fermi_momcut_sim_2[2][5];
TH1F *h_mis_mom_fermi_2[2][5];
TH1F *h_mis_mom_fermi_sim_2[2][5];
TH1F *h_mis_mom_fermi_mmas_cut_2[2][5];
TH1F *h_mis_mom_fermi_mmas_cut_sim_2[2][5];



TH1F *h_pim_mis_fermi_nocut_3[2][17];
TH1F *h_pim_mis_fermi_nocut_sim_3[2][17];
TH1F *h_pim_mis_fermi_momcut_3[2][4];
TH1F *h_pim_mis_fermi_momcut_sim_3[2][4];
TH1F *h_mis_mom_fermi_3[2][4];
TH1F *h_mis_mom_fermi_sim_3[2][4];
TH1F *h_mis_mom_fermi_mmas_cut_3[2][4];
TH1F *h_mis_mom_fermi_mmas_cut_sim_3[2][4];


TH1F *h_pim_mis_fermi_nocut_4[14];
TH1F *h_pim_mis_fermi_nocut_sim_4[14];
TH1F *h_pim_mis_fermi_momcut_4[4];
TH1F *h_pim_mis_fermi_momcut_sim_4[4];
TH1F *h_mis_mom_fermi_4[4];
TH1F *h_mis_mom_fermi_sim_4[4];
TH1F *h_mis_mom_fermi_mmas_cut_4[4];
TH1F *h_mis_mom_fermi_mmas_cut_sim_4[4];


TH1F *h_pim_mis_fermi_nocut_5[13];
TH1F *h_pim_mis_fermi_nocut_sim_5[13];
TH1F *h_pim_mis_fermi_momcut_5[3];
TH1F *h_pim_mis_fermi_momcut_sim_5[3];
TH1F *h_mis_mom_fermi_5[3];
TH1F *h_mis_mom_fermi_sim_5[3];
TH1F *h_mis_mom_fermi_mmas_cut_5[3];
TH1F *h_mis_mom_fermi_mmas_cut_sim_5[3];


TH1F *h_pim_mis_fermi_nocut_6[12];
TH1F *h_pim_mis_fermi_nocut_sim_6[12];
TH1F *h_pim_mis_fermi_momcut_6[3];
TH1F *h_pim_mis_fermi_momcut_sim_6[3];
TH1F *h_mis_mom_fermi_6[3];
TH1F *h_mis_mom_fermi_sim_6[3];
TH1F *h_mis_mom_fermi_mmas_cut_6[3];
TH1F *h_mis_mom_fermi_mmas_cut_sim_6[3];


TH1F *h_pim_mis_fermi_nocut_7[10];
TH1F *h_pim_mis_fermi_nocut_sim_7[10];
TH1F *h_pim_mis_fermi_momcut_7[3];
TH1F *h_pim_mis_fermi_momcut_sim_7[3];
TH1F *h_mis_mom_fermi_7[3];
TH1F *h_mis_mom_fermi_sim_7[3];
TH1F *h_mis_mom_fermi_mmas_cut_7[3];
TH1F *h_mis_mom_fermi_mmas_cut_sim_7[3];
*/


//TH1F *h_pim_mis_all_reg_th_dep_pim[5][12];
//TH1F *h_pim_mis_all_reg_th_dep_pim_sim[5][12];

//TH1F *h_mom_all_reg_th_dep_pim[5][12];
//TH1F *h_mom_all_reg_th_dep_pim_sim[5][12];


//TH1F *h_pim_mis_all_reg_th_dep_pip[5][12];
//TH1F *h_pim_mis_all_reg_th_dep_pip_sim[5][12];

//TH1F *h_mom_all_reg_th_dep_pip[5][12];
//TH1F *h_mom_all_reg_th_dep_pip_sim[5][12];


//TH1F *h_pim_mis_all_reg_th_dep_pr[5][12];
//TH1F *h_pim_mis_all_reg_th_dep_pr_sim[5][12];

//TH1F *h_pim_mis_th_dep_pr[5][12];
//TH1F *h_pim_mis_th_dep_pr_sim[5][12];


//TH1F *h_pim_mis_th_dep[5][12][8];
//TH1F *h_pim_mis_th_dep_sim[5][12][8];



//TH1F *h_mom_all_reg_th_dep_pr[5][12];
//TH1F *h_mom_all_reg_th_dep_pr_sim[5][12];

TH1F *h_pim_mis_all_reg[5];
TH1F *h_pim_mis_all_reg_sim[5];

TH1F *h_pip_mis_all_reg[5];
TH1F *h_pip_mis_all_reg_sim[5];

TH1F *h_pr_mis_all_reg[5];
TH1F *h_pr_mis_all_reg_sim[5];



TH1F *h_mom_all_reg[5];
TH1F *h_mom_all_reg_sim[5];




TH1F *h_sim_mom_corr_test = new TH1F ("h_sim_mom_corr_test","h_sim_mom_corr_test",200,-0.05,0.05);


TH1F *h_z_P = new TH1F ("h_z_p","h_z_p",200,-5.,5.);
TH1F *h_z_PIp = new TH1F ("h_z_pip","h_z_pip",200,-5.,5.);
TH1F *h_z_PIm = new TH1F ("h_z_pim","h_z_pim",200,-5.,5.);

TH1F *h_z_P_sim = new TH1F ("h_z_p_sim","h_z_p_sim",200,-5.,5.);
TH1F *h_z_PIp_sim = new TH1F ("h_z_pip_sim","h_z_pip_sim",200,-5.,5.);
TH1F *h_z_PIm_sim = new TH1F ("h_z_pim_sim","h_z_pim_sim",200,-5.,5.);













TH2F *h_cc_nphe_total_s1 = new TH2F("h_cc_nphe_total_s1","h_cc_nphe_total_s1",200, -5., 55.,200., -25., 25.);
TH2F *h_cc_nphe_total_s2 = new TH2F("h_cc_nphe_total_s2","h_cc_nphe_total_s2",200, -5., 55.,200., -25., 25.);
TH2F *h_cc_nphe_total_s3 = new TH2F("h_cc_nphe_total_s3","h_cc_nphe_total_s3",200, -5., 55.,200., -25., 25.);
TH2F *h_cc_nphe_total_s4 = new TH2F("h_cc_nphe_total_s4","h_cc_nphe_total_s4",200, -5., 55.,200., -25., 25.);
TH2F *h_cc_nphe_total_s5 = new TH2F("h_cc_nphe_total_s5","h_cc_nphe_total_s5",200, -5., 55.,200., -25., 25.);
TH2F *h_cc_nphe_total_s6 = new TH2F("h_cc_nphe_total_s6","h_cc_nphe_total_s6",200, -5., 55.,200., -25., 25.);


TH2F *h_cc_nphe_final_s1 = new TH2F("h_cc_nphe_final_s1","h_cc_nphe_final_s1",200, -5., 55.,200., -25., 25.);
TH2F *h_cc_nphe_final_s2 = new TH2F("h_cc_nphe_final_s2","h_cc_nphe_final_s2",200, -5., 55.,200., -25., 25.);
TH2F *h_cc_nphe_final_s3 = new TH2F("h_cc_nphe_final_s3","h_cc_nphe_final_s3",200, -5., 55.,200., -25., 25.);
TH2F *h_cc_nphe_final_s4 = new TH2F("h_cc_nphe_final_s4","h_cc_nphe_final_s4",200, -5., 55.,200., -25., 25.);
TH2F *h_cc_nphe_final_s5 = new TH2F("h_cc_nphe_final_s5","h_cc_nphe_final_s5",200, -5., 55.,200., -25., 25.);
TH2F *h_cc_nphe_final_s6 = new TH2F("h_cc_nphe_final_s6","h_cc_nphe_final_s6",200, -5., 55.,200., -25., 25.);




TH2F *hist_ectot_sector1 = new TH2F("ectot_sector1","ectot_sector1",220, -0.005, 2.195,80, -0.005, 0.795);
TH2F *hist_ectot_sector2 = new TH2F("ectot_sector2","ectot_sector2",220, -0.005, 2.195,80, -0.005, 0.795);
TH2F *hist_ectot_sector3 = new TH2F("ectot_sector3","ectot_sector3",220, -0.005, 2.195,80, -0.005, 0.795);
TH2F *hist_ectot_sector4 = new TH2F("ectot_sector4","ectot_sector4",220, -0.005, 2.195,80, -0.005, 0.795);
TH2F *hist_ectot_sector5 = new TH2F("ectot_sector5","ectot_sector5",220, -0.005, 2.195,80, -0.005, 0.795);
TH2F *hist_ectot_sector6 = new TH2F("ectot_sector6","ectot_sector6",220, -0.005, 2.195,80, -0.005, 0.795);

TH2F *hist_ectot_sector1_sim = new TH2F("ectot_sector1_sim","ectot_sector1_sim",440, -0.005, 2.195,160, -0.005, 0.795);
TH2F *hist_ectot_sector2_sim = new TH2F("ectot_sector2_sim","ectot_sector2_sim",440, -0.005, 2.195,160, -0.005, 0.795);
TH2F *hist_ectot_sector3_sim = new TH2F("ectot_sector3_sim","ectot_sector3_sim",440, -0.005, 2.195,160, -0.005, 0.795);
TH2F *hist_ectot_sector4_sim = new TH2F("ectot_sector4_sim","ectot_sector4_sim",440, -0.005, 2.195,160, -0.005, 0.795);
TH2F *hist_ectot_sector5_sim = new TH2F("ectot_sector5_sim","ectot_sector5_sim",440, -0.005, 2.195,160, -0.005, 0.795);
TH2F *hist_ectot_sector6_sim = new TH2F("ectot_sector6_sim","ectot_sector6_sim",440, -0.005, 2.195,160, -0.005, 0.795);


TH2F *th_cc_vs_seg_1 = new TH2F("th_cc_vs_seg_1","th_cc_vs_seg_1",20, -0.5, 19.5,400, 0., 50.);
TH2F *th_cc_vs_seg_2 = new TH2F("th_cc_vs_seg_2","th_cc_vs_seg_2",20, -0.5, 19.5,400, 0., 50.);
TH2F *th_cc_vs_seg_3 = new TH2F("th_cc_vs_seg_3","th_cc_vs_seg_3",20, -0.5, 19.5,400, 0., 50.);
TH2F *th_cc_vs_seg_4 = new TH2F("th_cc_vs_seg_4","th_cc_vs_seg_4",20, -0.5, 19.5,400, 0., 50.);
TH2F *th_cc_vs_seg_5 = new TH2F("th_cc_vs_seg_5","th_cc_vs_seg_5",20, -0.5, 19.5,400, 0., 50.);
TH2F *th_cc_vs_seg_6 = new TH2F("th_cc_vs_seg_6","th_cc_vs_seg_6",20, -0.5, 19.5,400, 0., 50.);

TH2F *th_cc_vs_seg_1_sim = new TH2F("th_cc_vs_seg_1_sim","th_cc_vs_seg_1_sim",20, -0.5, 19.5,400, 0.,50.);
TH2F *th_cc_vs_seg_2_sim = new TH2F("th_cc_vs_seg_2_sim","th_cc_vs_seg_2_sim",20, -0.5, 19.5,400, 0.,50.);
TH2F *th_cc_vs_seg_3_sim = new TH2F("th_cc_vs_seg_3_sim","th_cc_vs_seg_3_sim",20, -0.5, 19.5,400, 0.,50.);
TH2F *th_cc_vs_seg_4_sim = new TH2F("th_cc_vs_seg_4_sim","th_cc_vs_seg_4_sim",20, -0.5, 19.5,400, 0.,50.);
TH2F *th_cc_vs_seg_5_sim = new TH2F("th_cc_vs_seg_5_sim","th_cc_vs_seg_5_sim",20, -0.5, 19.5,400, 0.,50.);
TH2F *th_cc_vs_seg_6_sim = new TH2F("th_cc_vs_seg_6_sim","th_cc_vs_seg_6_sim",20, -0.5, 19.5,400, 0.,50.);




TH2F *W_2pi_selection = new TH2F("W_2pi_selection","W_2pi_selection",200, 1.2, 1.9,200, 0.2, 1.4);
TH2F *W_2pi_selection_sim = new TH2F("W_2pi_selection_sim","W_2pi_selection_sim",200, 1.2, 1.9,200, 0.2, 1.4);





TH1F *h_1d_rc_0425 = new TH1F("h_1d_rc_0425","h_1d_rc_0425",210, 1.3,1.825);
TH1F *h_1d_rc_0475 = new TH1F("h_1d_rc_0475","h_1d_rc_0475",210, 1.3,1.825);
TH1F *h_1d_rc_0525 = new TH1F("h_1d_rc_0525","h_1d_rc_0525",200, 1.3,1.8);
TH1F *h_1d_rc_0575 = new TH1F("h_1d_rc_0575","h_1d_rc_0575",200, 1.3,1.8);
TH1F *h_1d_rc_0625 = new TH1F("h_1d_rc_0625","h_1d_rc_0625",190, 1.3,1.775);
TH1F *h_1d_rc_0675 = new TH1F("h_1d_rc_0675","h_1d_rc_0675",180, 1.3,1.75);
TH1F *h_1d_rc_0725 = new TH1F("h_1d_rc_0725","h_1d_rc_0725",170, 1.3,1.725);
TH1F *h_1d_rc_0775 = new TH1F("h_1d_rc_0775","h_1d_rc_0775",160, 1.3,1.7);
TH1F *h_1d_rc_0825 = new TH1F("h_1d_rc_0825","h_1d_rc_0825",140, 1.3,1.65);
TH1F *h_1d_rc_0875 = new TH1F("h_1d_rc_0875","h_1d_rc_0875",130, 1.3,1.625);
TH1F *h_1d_rc_0925 = new TH1F("h_1d_rc_0925","h_1d_rc_0925",120, 1.3,1.6);
TH1F *h_1d_rc_0975 = new TH1F("h_1d_rc_0975","h_1d_rc_0975",100, 1.3,1.55);
/*
TH1F *h_1d_rc_0425_evt = new TH1F("h_1d_rc_0425_evt","h_1d_rc_0425_evt",21, 1.3,1.825);
TH1F *h_1d_rc_0475_evt = new TH1F("h_1d_rc_0475_evt","h_1d_rc_0475_evt",21, 1.3,1.825);
TH1F *h_1d_rc_0525_evt = new TH1F("h_1d_rc_0525_evt","h_1d_rc_0525_evt",20, 1.3,1.8);
TH1F *h_1d_rc_0575_evt = new TH1F("h_1d_rc_0575_evt","h_1d_rc_0575_evt",20, 1.3,1.8);
TH1F *h_1d_rc_0625_evt = new TH1F("h_1d_rc_0625_evt","h_1d_rc_0625_evt",19, 1.3,1.775);
TH1F *h_1d_rc_0675_evt = new TH1F("h_1d_rc_0675_evt","h_1d_rc_0675_evt",18, 1.3,1.75);
TH1F *h_1d_rc_0725_evt = new TH1F("h_1d_rc_0725_evt","h_1d_rc_0725_evt",17, 1.3,1.725);
TH1F *h_1d_rc_0775_evt = new TH1F("h_1d_rc_0775_evt","h_1d_rc_0775_evt",16, 1.3,1.7);
TH1F *h_1d_rc_0825_evt = new TH1F("h_1d_rc_0825_evt","h_1d_rc_0825_evt",14, 1.3,1.65);
TH1F *h_1d_rc_0875_evt = new TH1F("h_1d_rc_0875_evt","h_1d_rc_0875_evt",13, 1.3,1.625);
TH1F *h_1d_rc_0925_evt = new TH1F("h_1d_rc_0925_evt","h_1d_rc_0925_evt",12, 1.3,1.6);
TH1F *h_1d_rc_0975_evt = new TH1F("h_1d_rc_0975_evt","h_1d_rc_0975_evt",10, 1.3,1.55);*/

TH1F *h_inprot_miss = new TH1F("h_inprot_miss","h_inprot_miss",800, -2., 2.);
TH1F *h_inprot_miss_en = new TH1F("h_inprot_miss_en","h_inprot_miss_en",800, -2., 2.);

TH1F *h_inprot_miss_sim = new TH1F("h_inprot_miss_sim","h_inprot_miss_sim",800, -2., 2.);
TH1F *h_inprot_miss_en_sim = new TH1F("h_inprot_miss_en_sim","h_inprot_miss_en_sim",800, -2., 2.);

TH1F *hist_PIm_miss = new TH1F("hist_PIm_miss","hist_PIm_miss",400, -0.5, 0.5);
TH1F *h_PIm_miss_sim = new TH1F("h_PIm_miss_sim","h_PIm_miss_sim",400, -0.5, 0.5);

TH1F *hist_w_hadr_all_reg = new TH1F("hist_w_hadr_all_reg","hist_w_hadr_all_reg",210., 1.2, 1.9);
TH1F *hist_w_el_all_reg = new TH1F("hist_w_el_all_reg","hist_w_el_all_reg",210., 1.2, 1.9);
TH1F *hist_w_sim_old_1dim = new TH1F("hist_w_sim_old_1dim","hist_w_sim_old_1dim",210., 1.2, 1.9);
TH1F *hist_w_sim_new_1dim = new TH1F("hist_w_sim_new_1dim","hist_w_sim_new_1dim",210., 1.2, 1.9);


TH1F *hist_PIm_miss_all_reg_1 = new TH1F("hist_PIm_miss_all_reg_1","hist_PIm_miss_all_reg_1",400, -0.2, 0.2);
TH1F *hist_PIm_miss_all_reg_2 = new TH1F("hist_PIm_miss_all_reg_2","hist_PIm_miss_all_reg_2",400, -0.2, 0.2);

TH1F *hist_PIm_miss_all_reg_1_sim = new TH1F("hist_PIm_miss_all_reg_1_sim","hist_PIm_miss_all_reg_1_sim",400, -0.2, 0.2);
TH1F *hist_PIm_miss_all_reg_2_sim = new TH1F("hist_PIm_miss_all_reg_2_sim","hist_PIm_miss_all_reg_2_sim",400, -0.2, 0.2);


TH1F *h_PIp_miss_sim = new TH1F("h_PIp_miss_sim","h_PIp_miss_sim",500, -0.7, 0.3);
TH1F *h_miss_mass_0_sim = new TH1F("h_miss_mass_0_sim","h_miss_mass_0_sim",400, -0.2, 0.2);
TH1F *h_miss_en_0_sim = new TH1F("h_miss_en_0_sim","h_miss_en_0_sim",400, -2, 2);


TH1F *h_miss_mass_0 = new TH1F("h_miss_mass_0","h_miss_mass_0",400, -0.2, 0.2);
TH1F *h_PIm_miss_en_sim = new TH1F("h_PIm_miss_en_sim","h_PIm_miss_en_sim",400, -2, 2);
TH1F *h_PIp_miss_en_sim = new TH1F("h_PIp_miss_en_sim","h_PIp_miss_en_sim",400, -2, 2);

TH1F *hist_PIm_miss_en = new TH1F("hist_PIm_miss_en","hist_PIm_miss_en",400, -2, 2);

TH1F *h_miss_mass_0_d = new TH1F("h_miss_mass_0_d","h_miss_mass_0_d",2000, -0., 2.);
TH1F *h_miss_mass_0_d_sim = new TH1F("h_miss_mass_0_d_sim","h_miss_mass_0_d_sim",2000, -0., 2.);


TH1F *h_miss_mass_0_d_mmcut = new TH1F("h_miss_mass_0_d_mmcut","h_miss_mass_0_d_mmcut",2000, -0., 2.);
TH1F *h_miss_mass_0_d_sim_mmcut = new TH1F("h_miss_mass_0_d_sim_mmcut","h_miss_mass_0_d_sim_mmcut",2000, -0., 2.);

TH1F *h_miss_mom_0 = new TH1F("h_miss_mom_0","h_miss_mom_0",400, -0.1, 0.9);
TH1F *h_miss_mom_0_sim = new TH1F("h_miss_mom_0_sim","h_miss_mom_0_sim",400, -0.1, 0.9);


TH1F *h_miss_mom_0_nocut = new TH1F("h_miss_mom_0_nocut","h_miss_mom_0_nocut",400, -0.1, 0.9);
TH1F *h_miss_mom_0_cut_on0 = new TH1F("h_miss_mom_0_cut_on0","h_miss_mom_0_cut_on0",400, -0.1, 0.9);
TH1F *h_miss_mom_0_cut_onpim = new TH1F("h_miss_mom_0_cut_onpim","h_miss_mom_0_cut_onpim",400, -0.1, 0.9);

TH1F *h_miss_mom_0_nocut_sim = new TH1F("h_miss_mom_0_nocut_sim","h_miss_mom_0_nocut_sim",400, -0.1, 0.9);
TH1F *h_miss_mom_0_cut_on0_sim = new TH1F("h_miss_mom_0_cut_on0_sim","h_miss_mom_0_cut_on0_sim",400, -0.1, 0.9);
TH1F *h_miss_mom_0_cut_onpim_sim = new TH1F("h_miss_mom_0_cut_onpim_sim","h_miss_mom_0_cut_onpim_sim",400, -0.1, 0.9);


TH1F *h_miss_mom_0_d = new TH1F("h_miss_mom_0_d","h_miss_mom_0_d",400, -0.1, 0.9);
TH1F *h_miss_mom_0_d_sim = new TH1F("h_miss_mom_0_d_sim","h_miss_mom_0_d_sim",400, -0.1, 0.9);

TH1F *h_miss_mom_0_d_mmcut = new TH1F("h_miss_mom_0_d_mmcut","h_miss_mom_0_d_mmcut",400, -0.1, 0.9);
TH1F *h_miss_mom_0_d_sim_mmcut = new TH1F("h_miss_mom_0_d_sim_mmcut","h_miss_mom_0_d_sim_mmcut",400, -0.1, 0.9);
//extern 


TH1F *hist_PIp_miss = new TH1F("hist_PIp_miss","hist_PIp_miss",500, -0.7, 0.3);
TH1F *hist_PIp_miss_d = new TH1F("hist_PIp_miss_d","hist_PIp_miss_d",500, -1., 3.);
TH1F *hist_PIp_miss_d_bef = new TH1F("hist_PIp_miss_d_bef","hist_PIp_miss_d_bef",500, -1., 3.);
TH1F *h_PIp_miss_d_sim = new TH1F("h_PIp_miss_d_sim","h_PIp_miss_d_sim",500, -1., 3.);
TH1F *h_PIp_miss_d_bef_sim = new TH1F("h_PIp_miss_d_bef_sim","h_PIp_miss_d_bef_sim",500, -1., 3.);

TH1F *hist_PIp_miss_en = new TH1F("hist_PIp_miss_en","hist_PIp_miss_en",400, -2, 2);
TH1F *hist_P_miss = new TH1F("hist_P_miss","hist_P_miss",400, -0.5, 1.5);
TH1F *hist_P_miss_en = new TH1F("hist_P_miss_en","hist_P_miss_en",400, -2, 2);
TH1F *hist_miss_en_0 = new TH1F("hist_miss_en_0","hist_miss_en_0",400, -2, 2);
 

TH2F *h_mm_0_vs_npart = new TH2F("mm_0_vs_npart","mm_0_vs_npart",200, 0,10,300, -0.4, 0.2);
TH2F *h_mm_pim_vs_npart = new TH2F("mm_pim_vs_npart","mm_pim_vs_npart",200, 0,10,300, -0.4, 0.2);
TH2F *h_mm_pip_vs_npart = new TH2F("mm_pip_vs_npart","mm_pip_vs_npart",200, 0,10,300, -0.4, 0.2);

TH2F *h_mm_0_vs_npart_sim = new TH2F("mm_0_vs_npart_sim","mm_0_vs_npart_sim",200, 0,10,300, -0.4, 0.2);
TH2F *h_mm_pim_vs_npart_sim = new TH2F("mm_pim_vs_npart_sim","mm_pim_vs_npart_sim",200, 0,10,300, -0.4, 0.2);
TH2F *h_mm_pip_vs_npart_sim = new TH2F("mm_pip_vs_npart_sim","mm_pip_vs_npart_sim",200, 0,10,300, -0.4, 0.2);



TH1F *hist_z_el_1 = new TH1F("hist_z_el_1","hist_z_el_1",400, -10.,10.);
TH1F *hist_z_el_2 = new TH1F("hist_z_el_2","hist_z_el_2",400, -10.,10.);
TH1F *hist_z_el_3 = new TH1F("hist_z_el_3","hist_z_el_3",400, -10.,10.);
TH1F *hist_z_el_4 = new TH1F("hist_z_el_4","hist_z_el_4",400, -10.,10.);
TH1F *hist_z_el_5 = new TH1F("hist_z_el_5","hist_z_el_5",400, -10.,10.);
TH1F *hist_z_el_6 = new TH1F("hist_z_el_6","hist_z_el_6",400, -10.,10.);


TH1F *hist_z_el_1_empty = new TH1F("hist_z_el_1_empty","hist_z_el_1_empty",400, -10.,10.);
TH1F *hist_z_el_2_empty = new TH1F("hist_z_el_2_empty","hist_z_el_2_empty",400, -10.,10.);
TH1F *hist_z_el_3_empty = new TH1F("hist_z_el_3_empty","hist_z_el_3_empty",400, -10.,10.);
TH1F *hist_z_el_4_empty = new TH1F("hist_z_el_4_empty","hist_z_el_4_empty",400, -10.,10.);
TH1F *hist_z_el_5_empty = new TH1F("hist_z_el_5_empty","hist_z_el_5_empty",400, -10.,10.);
TH1F *hist_z_el_6_empty = new TH1F("hist_z_el_6_empty","hist_z_el_6_empty",400, -10.,10.);


TH1F *hist_z_el_1_sim_1 = new TH1F("hist_z_el_1_sim_1","hist_z_el_1_sim_1",400, -10.,10.);
TH1F *hist_z_el_2_sim_1 = new TH1F("hist_z_el_2_sim_1","hist_z_el_2_sim_1",400, -10.,10.);
TH1F *hist_z_el_3_sim_1 = new TH1F("hist_z_el_3_sim_1","hist_z_el_3_sim_1",400, -10.,10.);
TH1F *hist_z_el_4_sim_1 = new TH1F("hist_z_el_4_sim_1","hist_z_el_4_sim_1",400, -10.,10.);
TH1F *hist_z_el_5_sim_1 = new TH1F("hist_z_el_5_sim_1","hist_z_el_5_sim_1",400, -10.,10.);
TH1F *hist_z_el_6_sim_1 = new TH1F("hist_z_el_6_sim_1","hist_z_el_6_sim_1",400, -10.,10.);

TH1F *hist_z_el_1_sim_2 = new TH1F("hist_z_el_1_sim_2","hist_z_el_1_sim_2",400, -10.,10.);
TH1F *hist_z_el_2_sim_2 = new TH1F("hist_z_el_2_sim_2","hist_z_el_2_sim_2",400, -10.,10.);
TH1F *hist_z_el_3_sim_2 = new TH1F("hist_z_el_3_sim_2","hist_z_el_3_sim_2",400, -10.,10.);
TH1F *hist_z_el_4_sim_2 = new TH1F("hist_z_el_4_sim_2","hist_z_el_4_sim_2",400, -10.,10.);
TH1F *hist_z_el_5_sim_2 = new TH1F("hist_z_el_5_sim_2","hist_z_el_5_sim_2",400, -10.,10.);
TH1F *hist_z_el_6_sim_2 = new TH1F("hist_z_el_6_sim_2","hist_z_el_6_sim_2",400, -10.,10.);




//extern TH1F *inv_m_pip_pim_bin, *inv_m_pip_p_bin, *inv_m_pim_p_bin, *theta_PIm_cm, *theta_PIp_cm, *theta_P_cm, *phi_PIm_cm, *phi_PIp_cm, *phi_P_cm, *alpha_PIpPIm_pipf, *alpha_PPIp_piPIm, *alpha_PPIm_piPIp;

TH1F *h_inv_m_pip_pim_bin =  new TH1F("inv_m_pip_pim_bin","inv_m_pip_pim_bin",12.,0.228,0.69);
TH1F *h_inv_m_pip_p_bin =  new TH1F("inv_m_pip_p_bin","inv_m_pip_p_bin",12.,1.,1.485);
TH1F *h_inv_m_pim_p_bin =  new TH1F("inv_m_pim_p_bin","inv_m_pim_p_bin",12.,1.,1.485);
TH1F *h_theta_PIm_cm = new TH1F("theta_PIm_cm","theta_PIm_cm",10.,0.,180.);
TH1F *h_theta_PIp_cm = new TH1F("theta_PIp_cm","theta_PIp_cm",10.,0.,180.);
TH1F *h_theta_P_cm = new TH1F("theta_P_cm","theta_P_cm",10.,0.,180.);
TH1F *h_phi_PIm_cm = new TH1F("phi_PIm_cm","phi_PIm_cm",8.,0.,360.);
TH1F *h_phi_PIp_cm = new TH1F("phi_PIp_cm","phi_PIp_cm",8.,0.,360.);
TH1F *h_phi_P_cm = new TH1F("phi_P_cm","phi_P_cm",8.,0.,360.);
TH1F *h_alpha_PIpPIm_pipf = new TH1F("alpha_PIpPIm_pipf","alpha_PIpPIm_pipf",8.,0.,360.);
TH1F *h_alpha_PPIp_piPIm = new TH1F("alpha_PPIp_piPIm","alpha_PPIp_piPIm",8.,0.,360.);
TH1F *h_alpha_PPIm_piPIp = new TH1F("alpha_PPIm_piPIp","alpha_PPIm_piPIp",8.,0.,360.);

static const Int_t ndims = 5;
  Int_t bins[5];
  
  
  Double_t W_bin [21];
  
   Double_t xmin[5] = {(0.938272 + 0.13957), (0.13957 + 0.13957),0.,0.,0.};
   Double_t xmax_1[21];
   Double_t xmax_2[21];
   
   Double_t xmax[5];






TH1F *h_inv_m_pip_pim = new TH1F("inv_m_pip_pim","inv_m_pip_pim",210, 0.,2.1);
TH1F *h_inv_m_pip_p = new TH1F("inv_m_pip_p","inv_m_pip_p",210, 0.,2.1);
TH1F *h_inv_m_pim_p = new TH1F("inv_m_pim_p","inv_m_pim_p",210, 0.,2.1);

TH1F *h_inv_m_pip_pim_sim = new TH1F("inv_m_pip_pim_sim","inv_m_pip_pim_sim",110, 0.,2.1);
TH1F *h_inv_m_pip_p_sim = new TH1F("inv_m_pip_p_sim","inv_m_pip_p_sim",110, 0.,2.1);
TH1F *h_inv_m_pim_p_sim = new TH1F("inv_m_pim_p_sim","inv_m_pim_p_sim",110, 0.,2.1);


TH2F *ph_vs_th_1 = new TH2F("ph_vs_th_1","ph_vs_th_1",200, 0., 60.,150., -35., 35.);
TH2F *ph_vs_th_2 = new TH2F("ph_vs_th_2","ph_vs_th_2",200, 0., 60.,150., -35., 35.);
TH2F *ph_vs_th_3 = new TH2F("ph_vs_th_3","ph_vs_th_3",200, 0., 60.,150., -35., 35.);
TH2F *ph_vs_th_4 = new TH2F("ph_vs_th_4","ph_vs_th_4",200, 0., 60.,150., -35., 35.);
TH2F *ph_vs_th_5 = new TH2F("ph_vs_th_5","ph_vs_th_5",200, 0., 60.,150., -35., 35.);
TH2F *ph_vs_th_6 = new TH2F("ph_vs_th_6","ph_vs_th_6",200, 0., 60.,150., -35., 35.);

TH1D *h_cos_th;

ostringstream qqq;
Short_t i,j,k,p;
TH2F *ph_vs_th_1pe[17],*ph_vs_th_2pe[17],*ph_vs_th_3pe[17],*ph_vs_th_4pe[17],*ph_vs_th_5pe[17],*ph_vs_th_6pe[17];
TH2F *ph_vs_th_1pe_fid[17],*ph_vs_th_2pe_fid[17],*ph_vs_th_3pe_fid[17],*ph_vs_th_4pe_fid[17],*ph_vs_th_5pe_fid[17],*ph_vs_th_6pe_fid[17];
TH2F *beta_vs_p_p[6][48];
TH2F *beta_vs_p_pip[6][48];
TH2F *beta_vs_p_pim[6][48]; 

TH2F *beta_vs_p_p_sim[6][48];
TH2F *beta_vs_p_pip_sim[6][48];
TH2F *beta_vs_p_pim_sim[6][48]; 


TH2F *ph_vs_th_pim[6][10];
TH2F *ph_th_pim_all_p[6]; 

TH2F *time_pip[6][48];
TH2F *time_pim[6][48];
TH2F *time_p[6][48];

TH2F *time_pip_sim[6][48];
TH2F *time_pim_sim[6][48];
TH2F *time_p_sim[6][48];

TH2F *ph_vs_th_el_sim[6][7];
TH2F *ph_vs_th_p_sim[6];
TH2F *ph_vs_th_pip_sim[6];
TH2F *ph_vs_th_pim_sim[6][5];



//photoelectrons
TH1F *ph_el_left[6][20];
TH1F *ph_el_both[6][20];
TH1F *ph_el_right[6][20];

/*
THnSparseD *h_5dim_1[12][21];
THnSparseD *h_5dim_2[12][21];
THnSparseD *h_5dim_3[12][21];


THnSparseD *h_5dim_1_empty[12][21];
THnSparseD *h_5dim_2_empty[12][21];
THnSparseD *h_5dim_3_empty[12][21];

THnSparseD *h_5dim_1_sim_1[12][21];
THnSparseD *h_5dim_2_sim_1[12][21];
THnSparseD *h_5dim_3_sim_1[12][21];

THnSparseD *h_5dim_1_sim_2[12][21];
THnSparseD *h_5dim_2_sim_2[12][21];
THnSparseD *h_5dim_3_sim_2[12][21];
*/
THnSparseD *h_5dim_excl_1[12][21];
THnSparseD *h_5dim_excl_2[12][21];
THnSparseD *h_5dim_excl_3[12][21];


THnSparseD *h_5dim_excl_1_sim[12][21];
THnSparseD *h_5dim_excl_2_sim[12][21];
THnSparseD *h_5dim_excl_3_sim[12][21];

THnSparseD *h_5dim_excl_1_sim_evt[12][21];
THnSparseD *h_5dim_excl_2_sim_evt[12][21];
THnSparseD *h_5dim_excl_3_sim_evt[12][21];

THnSparseD *h_5dim_pim_1[12][21];
THnSparseD *h_5dim_pim_2[12][21];
THnSparseD *h_5dim_pim_3[12][21];

THnSparseD *h_5dim_pim_1_sim[12][21];
THnSparseD *h_5dim_pim_2_sim[12][21];
THnSparseD *h_5dim_pim_3_sim[12][21];

THnSparseD *h_5dim_pim_1_sim_evt[12][21];
THnSparseD *h_5dim_pim_2_sim_evt[12][21];
THnSparseD *h_5dim_pim_3_sim_evt[12][21];

THnSparseD *h_5dim_1_sim_gen[12][21];
THnSparseD *h_5dim_2_sim_gen[12][21];
THnSparseD *h_5dim_3_sim_gen[12][21];

THnSparseD *h_5dim_1_sim_gen_evt[12][21];
THnSparseD *h_5dim_2_sim_gen_evt[12][21];
THnSparseD *h_5dim_3_sim_gen_evt[12][21];





TH1D *h_w_int[12];

TH2F *th_vs_p_e_1[6], *th_vs_p_e_2[6], *th_vs_p_p_1[6],*th_vs_p_p_2[6],*th_vs_p_pip_1[6],*th_vs_p_pip_2[6],*th_vs_p_pim_1[6],*th_vs_p_pim_2[6];
TH2F *th_vs_p_e_1_sim[6], *th_vs_p_e_2_sim[6], *th_vs_p_p_1_sim[6],*th_vs_p_p_2_sim[6],*th_vs_p_pip_1_sim[6],*th_vs_p_pip_2_sim[6],*th_vs_p_pim_1_sim[6],*th_vs_p_pim_2_sim[6];
// booking proton histograms

TH2F *ph_vs_th_p_1 = new TH2F("ph_vs_th_p_1","ph_vs_th_p_1",220, 0., 110.,150., -35., 35.);
TH2F *ph_vs_th_p_2 = new TH2F("ph_vs_th_p_2","ph_vs_th_p_2",220, 0., 110.,150., -35., 35.);
TH2F *ph_vs_th_p_3 = new TH2F("ph_vs_th_p_3","ph_vs_th_p_3",220, 0., 110.,150., -35., 35.);
TH2F *ph_vs_th_p_4 = new TH2F("ph_vs_th_p_4","ph_vs_th_p_4",220, 0., 110.,150., -35., 35.);
TH2F *ph_vs_th_p_5 = new TH2F("ph_vs_th_p_5","ph_vs_th_p_5",220, 0., 110.,150., -35., 35.);
TH2F *ph_vs_th_p_6 = new TH2F("ph_vs_th_p_6","ph_vs_th_p_6",220, 0., 110.,150., -35., 35.);

TH2F *ph_vs_th_p_1_w = new TH2F("ph_vs_th_p_1_w","ph_vs_th_p_1_w",220, 0., 110.,150., -35., 35.);
TH2F *ph_vs_th_p_2_w = new TH2F("ph_vs_th_p_2_w","ph_vs_th_p_2_w",220, 0., 110.,150., -35., 35.);
TH2F *ph_vs_th_p_3_w = new TH2F("ph_vs_th_p_3_w","ph_vs_th_p_3_w",220, 0., 110.,150., -35., 35.);
TH2F *ph_vs_th_p_4_w = new TH2F("ph_vs_th_p_4_w","ph_vs_th_p_4_w",220, 0., 110.,150., -35., 35.);
TH2F *ph_vs_th_p_5_w = new TH2F("ph_vs_th_p_5_w","ph_vs_th_p_5_w",220, 0., 110.,150., -35., 35.);
TH2F *ph_vs_th_p_6_w = new TH2F("ph_vs_th_p_6_w","ph_vs_th_p_6_w",220, 0., 110.,150., -35., 35.);


TH2F *ph_vs_th_pip_1 = new TH2F("ph_vs_th_pip_1","ph_vs_th_pip_1",220, 0., 110.,150., -35., 35.);
TH2F *ph_vs_th_pip_2 = new TH2F("ph_vs_th_pip_2","ph_vs_th_pip_2",220, 0., 110.,150., -35., 35.);
TH2F *ph_vs_th_pip_3 = new TH2F("ph_vs_th_pip_3","ph_vs_th_pip_3",220, 0., 110.,150., -35., 35.);
TH2F *ph_vs_th_pip_4 = new TH2F("ph_vs_th_pip_4","ph_vs_th_pip_4",220, 0., 110.,150., -35., 35.);
TH2F *ph_vs_th_pip_5 = new TH2F("ph_vs_th_pip_5","ph_vs_th_pip_5",220, 0., 110.,150., -35., 35.);
TH2F *ph_vs_th_pip_6 = new TH2F("ph_vs_th_pip_6","ph_vs_th_pip_6",220, 0., 110.,150., -35., 35.);

TH2F  *h_delta_w_vs_w_old_data = new TH2F("h_delta_w_vs_w_old_data","h_delta_w_vs_w_old_data",180, 1.3, 1.9,300., -0.3,0.3);
TH2F  *h_delta_w_vs_w_old_data_mmcut = new TH2F("h_delta_w_vs_w_old_data_mmcut","h_delta_w_vs_w_old_data_mmcut",180, 1.3, 1.9,300., -0.3,0.3);
TH2F  *h_delta_w_vs_w_old_gen = new TH2F("h_delta_w_vs_w_old_gen","h_delta_w_vs_w_old_gen",180, 1.3, 1.9,300., -0.3,0.3);
TH2F  *h_delta_w_vs_w_old_rec = new TH2F("h_delta_w_vs_w_old_rec","h_delta_w_vs_w_old_rec",180, 1.3, 1.9,300., -0.3,0.3);

TH2F  *ph_th_p_1[5],*ph_th_p_2[5],*ph_th_p_3[5],*ph_th_p_4[5],*ph_th_p_5[5],*ph_th_p_6[5];

TH2F  *ph_th_pip_1[6],*ph_th_pip_2[6],*ph_th_pip_3[6],*ph_th_pip_4[6],*ph_th_pip_5[6],*ph_th_pip_6[6];

//output file

TFile *outFile;

// bookking quality check histograms

TH1F *hist_ltime = new TH1F("ltime","ltime",35000, -.5, 34999.5);
TH1F *hist_ltime_1d = new TH1F("ltime_1d","ltime_1d",1100, 0., 1.1);
TH1F *hist_n_incl = new TH1F("n_incl","n_incl",35000, -.5, 34999.5);
TH1F *hist_n_incl_1d = new TH1F("n_incl_1d","n_incl_1d",1000, 0., 300000.);
TH1F *hist_n_elast = new TH1F("n_elast","n_elast",35000, -.5, 34999.5);
TH1F *hist_n_elast_1d = new TH1F("n_elast_1d","n_elast_1d",1000, 0., 80000.);



int global() {

for (i=0; i<5; i++) {

qqq << "h_inv_NP_" <<i;
h_inv_NP[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, 0.938+0.939-0.05,1.4+0.1*i);
qqq.str("");

qqq << "h_inv_NPIp_" <<i;
h_inv_NPIp[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400,  0.140+0.939-0.05,1.4+0.1*i);
qqq.str("");

qqq << "h_inv_NPIm_" <<i;
h_inv_NPIm[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400,  0.140+0.939-0.05,1.4+0.1*i);
qqq.str("");


qqq << "h_pim_mis_all_reg_" <<i;
h_pim_mis_all_reg[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.4, 0.4);
qqq.str("");

qqq << "h_pim_mis_all_reg_sim_" <<i;
h_pim_mis_all_reg_sim[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.4, 0.4);
qqq.str("");

qqq << "h_pip_mis_all_reg_" <<i;
h_pip_mis_all_reg[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.4, 0.4);
qqq.str("");

qqq << "h_pip_mis_all_reg_sim_" <<i;
h_pip_mis_all_reg_sim[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.4, 0.4);
qqq.str("");

qqq << "h_pr_mis_all_reg_" <<i;
h_pr_mis_all_reg[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, 0.5, 1.3);
qqq.str("");

qqq << "h_pr_mis_all_reg_sim_" <<i;
h_pr_mis_all_reg_sim[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, 0.5, 1.3);
qqq.str("");



qqq << "h_mom_all_reg_" <<i;
h_mom_all_reg[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");

qqq << "h_mom_all_reg_sim_" <<i;
h_mom_all_reg_sim[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");
};

/*
for (j=0; j<12; j++) {


for (k=0; k<8; k++) {

qqq << "h_pim_mis_th_dep_arr_" <<i<<"_"<<j<<"_"<<k;
h_pim_mis_th_dep[i][j][k] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.05, 0.65);
qqq.str("");

qqq << "h_pim_mis_th_dep_arr_sim_" <<i<<"_"<<j<<"_"<<k;
h_pim_mis_th_dep_sim[i][j][k] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.05, 0.65);
qqq.str("");


};
};
for (j=0; j<12; j++) {

qqq << "h_pim_mis_all_reg_th_dep_pim_" <<i<<"_"<<j;
h_pim_mis_all_reg_th_dep_pim[i][j] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -1., 1.);
qqq.str("");

qqq << "h_pim_mis_all_reg_th_dep_pim_sim_" <<i<<"_"<<j;
h_pim_mis_all_reg_th_dep_pim_sim[i][j] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -1., 1.);
qqq.str("");

qqq << "h_mom_all_reg_th_dep_pim_" <<i<<"_"<<j;
h_mom_all_reg_th_dep_pim[i][j] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");

qqq << "h_mom_all_reg_th_dep_pim_sim_" <<i<<"_"<<j;
h_mom_all_reg_th_dep_pim_sim[i][j] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");


qqq << "h_pim_mis_all_reg_th_dep_pip_" <<i<<"_"<<j;
h_pim_mis_all_reg_th_dep_pip[i][j] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.5, 1.1);
qqq.str("");

qqq << "h_pim_mis_all_reg_th_dep_pip_sim_" <<i<<"_"<<j;
h_pim_mis_all_reg_th_dep_pip_sim[i][j] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.5, 1.1);
qqq.str("");

qqq << "h_mom_all_reg_th_dep_pip_" <<i<<"_"<<j;
h_mom_all_reg_th_dep_pip[i][j] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");

qqq << "h_mom_all_reg_th_dep_pip_sim_" <<i<<"_"<<j;
h_mom_all_reg_th_dep_pip_sim[i][j] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");





qqq << "h_mom_all_reg_th_dep_pr_" <<i<<"_"<<j;
h_mom_all_reg_th_dep_pr[i][j] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");

qqq << "h_mom_all_reg_th_dep_pr_sim_" <<i<<"_"<<j;
h_mom_all_reg_th_dep_pr_sim[i][j] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");



};
for (j=0; j<12; j++) {

qqq << "h_pim_mis_all_reg_th_dep_pr_" <<i<<"_"<<j;
h_pim_mis_all_reg_th_dep_pr[i][j] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.5, 1.1);
qqq.str("");

qqq << "h_pim_mis_all_reg_th_dep_pr_sim_" <<i<<"_"<<j;
h_pim_mis_all_reg_th_dep_pr_sim[i][j] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.5,1.1);
qqq.str("");


qqq << "h_pim_mis_th_dep_pr_" <<i<<"_"<<j;
h_pim_mis_th_dep_pr[i][j] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -1., 1.);
qqq.str("");

qqq << "h_pim_mis_th_dep_pr_sim_" <<i<<"_"<<j;
h_pim_mis_th_dep_pr_sim[i][j] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -1.,1.);
qqq.str("");


};
};*/


/*
for (j=0; j<4; j++) {
for (i=0; i<21; i++) {

qqq << "h_pim_mis_fermi_nocut_1_" <<j<<"_"<< i;
h_pim_mis_fermi_nocut_1[j][i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.05, 0.65);
qqq.str("");


qqq << "h_pim_mis_fermi_nocut_sim_1_" <<j<<"_"<< i;
h_pim_mis_fermi_nocut_sim_1[j][i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.05, 0.65);
qqq.str("");
};
for (i=0; i<5; i++) {
qqq << "h_pim_mis_fermi_momcut_1_" <<j<<"_"<< i;
h_pim_mis_fermi_momcut_1[j][i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.5, 0.5);
qqq.str("");


qqq << "h_pim_mis_fermi_momcut_sim_1_" <<j<<"_"<< i;
h_pim_mis_fermi_momcut_sim_1[j][i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.5, 0.5);
qqq.str("");


qqq << "h_mis_mom_fermi_1_" <<j<<"_"<< i;
h_mis_mom_fermi_1[j][i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");


qqq << "h_mis_mom_fermi_sim_1_" <<j<<"_"<< i;
h_mis_mom_fermi_sim_1[j][i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");


qqq << "h_mis_mom_fermi_mmas_cut_1_" <<j<<"_"<< i;
h_mis_mom_fermi_mmas_cut_1[j][i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");

qqq << "h_mis_mom_fermi_mmas_cut_sim_1_" <<j<<"_"<< i;
h_mis_mom_fermi_mmas_cut_sim_1[j][i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");

};
};


for (j=0; j<2; j++) {
for (i=0; i<19; i++) {

qqq << "h_pim_mis_fermi_nocut_2_" <<j<<"_"<< i;
h_pim_mis_fermi_nocut_2[j][i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.05, 0.65);
qqq.str("");


qqq << "h_pim_mis_fermi_nocut_sim_2_" <<j<<"_"<< i;
h_pim_mis_fermi_nocut_sim_2[j][i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.05, 0.65);
qqq.str("");


qqq << "h_pim_mis_fermi_momcut_2_" <<j<<"_"<< i;
h_pim_mis_fermi_momcut_2[j][i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -2.1, 0.5);
qqq.str("");


qqq << "h_pim_mis_fermi_momcut_sim_2_" <<j<<"_"<< i;
h_pim_mis_fermi_momcut_sim_2[j][i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -2.1, 0.5);
qqq.str("");


qqq << "h_mis_mom_fermi_2_" <<j<<"_"<< i;
h_mis_mom_fermi_2[j][i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");


qqq << "h_mis_mom_fermi_sim_2_" <<j<<"_"<< i;
h_mis_mom_fermi_sim_2[j][i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");


qqq << "h_mis_mom_fermi_mmas_cut_2_" <<j<<"_"<< i;
h_mis_mom_fermi_mmas_cut_2[j][i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");

qqq << "h_mis_mom_fermi_mmas_cut_sim_2_" <<j<<"_"<< i;
h_mis_mom_fermi_mmas_cut_sim_2[j][i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");

};
};



for (j=0; j<2; j++) {
for (i=0; i<17; i++) {

qqq << "h_pim_mis_fermi_nocut_3_" <<j<<"_"<< i;
h_pim_mis_fermi_nocut_3[j][i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.05, 0.65);
qqq.str("");


qqq << "h_pim_mis_fermi_nocut_sim_3_" <<j<<"_"<< i;
h_pim_mis_fermi_nocut_sim_3[j][i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.05, 0.65);
qqq.str("");


qqq << "h_pim_mis_fermi_momcut_3_" <<j<<"_"<< i;
h_pim_mis_fermi_momcut_3[j][i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -2.1, 0.5);
qqq.str("");


qqq << "h_pim_mis_fermi_momcut_sim_3_" <<j<<"_"<< i;
h_pim_mis_fermi_momcut_sim_3[j][i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -2.1, 0.5);
qqq.str("");


qqq << "h_mis_mom_fermi_3_" <<j<<"_"<< i;
h_mis_mom_fermi_3[j][i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");


qqq << "h_mis_mom_fermi_sim_3_" <<j<<"_"<< i;
h_mis_mom_fermi_sim_3[j][i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");


qqq << "h_mis_mom_fermi_mmas_cut_3_" <<j<<"_"<< i;
h_mis_mom_fermi_mmas_cut_3[j][i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");

qqq << "h_mis_mom_fermi_mmas_cut_sim_3_" <<j<<"_"<< i;
h_mis_mom_fermi_mmas_cut_sim_3[j][i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");

};
};




for (i=0; i<14; i++) {

qqq << "h_pim_mis_fermi_nocut_4_" << i;
h_pim_mis_fermi_nocut_4[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.05, 0.65);
qqq.str("");


qqq << "h_pim_mis_fermi_nocut_sim_4_" << i;
h_pim_mis_fermi_nocut_sim_4[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.05, 0.65);
qqq.str("");


qqq << "h_pim_mis_fermi_momcut_4_" << i;
h_pim_mis_fermi_momcut_4[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -2.1, 0.5);
qqq.str("");


qqq << "h_pim_mis_fermi_momcut_sim_4_" << i;
h_pim_mis_fermi_momcut_sim_4[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -2.1, 0.5);
qqq.str("");


qqq << "h_mis_mom_fermi_4_"<< i;
h_mis_mom_fermi_4[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");


qqq << "h_mis_mom_fermi_sim_4_"<< i;
h_mis_mom_fermi_sim_4[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");


qqq << "h_mis_mom_fermi_mmas_cut_4_"<< i;
h_mis_mom_fermi_mmas_cut_4[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");

qqq << "h_mis_mom_fermi_mmas_cut_sim_4_" << i;
h_mis_mom_fermi_mmas_cut_sim_4[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");

};


for (i=0; i<13; i++) {

qqq << "h_pim_mis_fermi_nocut_5_" << i;
h_pim_mis_fermi_nocut_5[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.05, 0.65);
qqq.str("");


qqq << "h_pim_mis_fermi_nocut_sim_5_" << i;
h_pim_mis_fermi_nocut_sim_5[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.05, 0.65);
qqq.str("");


qqq << "h_pim_mis_fermi_momcut_5_" << i;
h_pim_mis_fermi_momcut_5[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -2.1, 0.5);
qqq.str("");


qqq << "h_pim_mis_fermi_momcut_sim_5_" << i;
h_pim_mis_fermi_momcut_sim_5[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -2.1, 0.5);
qqq.str("");


qqq << "h_mis_mom_fermi_5_"<< i;
h_mis_mom_fermi_5[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");


qqq << "h_mis_mom_fermi_sim_5_"<< i;
h_mis_mom_fermi_sim_5[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");


qqq << "h_mis_mom_fermi_mmas_cut_5_"<< i;
h_mis_mom_fermi_mmas_cut_5[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");

qqq << "h_mis_mom_fermi_mmas_cut_sim_5_" << i;
h_mis_mom_fermi_mmas_cut_sim_5[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");





};


for (i=0; i<12; i++) {

qqq << "h_pim_mis_fermi_nocut_6_" << i;
h_pim_mis_fermi_nocut_6[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.05, 0.65);
qqq.str("");


qqq << "h_pim_mis_fermi_nocut_sim_6_" << i;
h_pim_mis_fermi_nocut_sim_6[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.05, 0.65);
qqq.str("");


qqq << "h_pim_mis_fermi_momcut_6_" << i;
h_pim_mis_fermi_momcut_6[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -2.1, 0.5);
qqq.str("");


qqq << "h_pim_mis_fermi_momcut_sim_6_" << i;
h_pim_mis_fermi_momcut_sim_6[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -2.1, 0.5);
qqq.str("");


qqq << "h_mis_mom_fermi_6_"<< i;
h_mis_mom_fermi_6[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");


qqq << "h_mis_mom_fermi_sim_6_"<< i;
h_mis_mom_fermi_sim_6[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");


qqq << "h_mis_mom_fermi_mmas_cut_6_"<< i;
h_mis_mom_fermi_mmas_cut_6[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");

qqq << "h_mis_mom_fermi_mmas_cut_sim_6_" << i;
h_mis_mom_fermi_mmas_cut_sim_6[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");

};

for (i=0; i<10; i++) {

qqq << "h_pim_mis_fermi_nocut_7_" << i;
h_pim_mis_fermi_nocut_7[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.05, 0.65);
qqq.str("");


qqq << "h_pim_mis_fermi_nocut_sim_7_" << i;
h_pim_mis_fermi_nocut_sim_7[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.05, 0.65);
qqq.str("");


qqq << "h_pim_mis_fermi_momcut_7_" << i;
h_pim_mis_fermi_momcut_7[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -2.1, 0.5);
qqq.str("");


qqq << "h_pim_mis_fermi_momcut_sim_7_" << i;
h_pim_mis_fermi_momcut_sim_7[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -2.1, 0.5);
qqq.str("");


qqq << "h_mis_mom_fermi_7_"<< i;
h_mis_mom_fermi_7[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");


qqq << "h_mis_mom_fermi_sim_7_"<< i;
h_mis_mom_fermi_sim_7[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");


qqq << "h_mis_mom_fermi_mmas_cut_7_"<< i;
h_mis_mom_fermi_mmas_cut_7[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");

qqq << "h_mis_mom_fermi_mmas_cut_sim_7_" << i;
h_mis_mom_fermi_mmas_cut_sim_7[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),400, -0.1, 1.9);
qqq.str("");

};*/


// loops for electron histograms

for (i=0; i<17; i++) {
qqq << "ph_vs_th_1pe[" << i << "]";
ph_vs_th_1pe[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),200, 0., 60.,150., -35., 35.);
qqq.str("");
};


for (i=0; i<17; i++) {
qqq << "ph_vs_th_2pe[" <<i << "]";
ph_vs_th_2pe[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),200, 0., 60.,150., -35., 35.);
qqq.str("");
};


for (i=0; i<17; i++) {
qqq << "ph_vs_th_3pe[" << i << "]";
ph_vs_th_3pe[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),200, 0., 60.,150., -35., 35.);
qqq.str("");
};


for (i=0; i<17; i++) {
qqq << "ph_vs_th_4pe[" << i << "]";
ph_vs_th_4pe[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),200, 0., 60.,150., -35., 35.);
qqq.str("");
};


for (i=0; i<17; i++) {
qqq << "ph_vs_th_5pe[" << i << "]";
ph_vs_th_5pe[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),200, 0., 60.,150., -35., 35.);
qqq.str("");
};


for (i=0; i<17; i++) {
qqq << "ph_vs_th_6pe[" << i << "]";
ph_vs_th_6pe[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),200, 0., 60.,150., -35., 35.);
qqq.str("");
};


for (i=0; i<17; i++) {
qqq << "ph_vs_th_1pe_fid[" << i << "]";
ph_vs_th_1pe_fid[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),200, 0., 60.,150., -35., 35.);
qqq.str("");
};


for (i=0; i<17; i++) {
qqq << "ph_vs_th_2pe_fid[" << i << "]";
ph_vs_th_2pe_fid[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),200, 0., 60.,150., -35., 35.);
qqq.str("");
};


for (i=0; i<17; i++) {
qqq << "ph_vs_th_3pe_fid[" << i << "]";
ph_vs_th_3pe_fid[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),200, 0., 60.,150., -35., 35.);
qqq.str("");
};


for (i=0; i<17; i++) {
qqq << "ph_vs_th_4pe_fid[" << i << "]";
ph_vs_th_4pe_fid[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),200, 0., 60.,150., -35., 35.);
qqq.str("");
};


for (i=0; i<17; i++) {
qqq << "ph_vs_th_5pe_fid[" << i << "]";
ph_vs_th_5pe_fid[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),200, 0., 60.,150., -35., 35.);
qqq.str("");
};


for (i=0; i<17; i++) {
qqq << "ph_vs_th_6pe_fid[" << i << "]";
ph_vs_th_6pe_fid[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),200, 0., 60.,150., -35., 35.);
qqq.str("");
};

//loops for proton histograms


for (i=0; i<5; i++) {
qqq << "ph_th_p_1[" << i << "]";
ph_th_p_1[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),220, 0., 110.,150., -35., 35.);
qqq.str("");
};

for (i=0; i<5; i++) {
qqq << "ph_th_p_2[" << i << "]";
ph_th_p_2[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),220, 0., 110.,150., -35., 35.);
qqq.str("");
};

for (i=0; i<5; i++) {
qqq << "ph_th_p_3[" << i << "]";
ph_th_p_3[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),220, 0., 110.,150., -35., 35.);
qqq.str("");
};

for (i=0; i<5; i++) {
qqq << "ph_th_p_4[" << i << "]";
ph_th_p_4[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),220, 0., 110.,150., -35., 35.);
qqq.str("");
};

for (i=0; i<5; i++) {
qqq << "ph_th_p_5[" << i << "]";
ph_th_p_5[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),220, 0., 110.,150., -35., 35.);
qqq.str("");
};

for (i=0; i<5; i++) {
qqq << "ph_th_p_6[" << i << "]";
ph_th_p_6[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),220, 0., 110.,150., -35., 35.);
qqq.str("");
};

//pip fid

for (i=0; i<6; i++) {
qqq << "ph_th_pip_1[" << i << "]";
ph_th_pip_1[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),220, 0., 110.,150., -35., 35.);
qqq.str("");
};

for (i=0; i<6; i++) {
qqq << "ph_th_pip_2[" << i << "]";
ph_th_pip_2[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),220, 0., 110.,150., -35., 35.);
qqq.str("");
};

for (i=0; i<6; i++) {
qqq << "ph_th_pip_3[" << i << "]";
ph_th_pip_3[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),220, 0., 110.,150., -35., 35.);
qqq.str("");
};

for (i=0; i<6; i++) {
qqq << "ph_th_pip_4[" << i << "]";
ph_th_pip_4[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),220, 0., 110.,150., -35., 35.);
qqq.str("");
};

for (i=0; i<6; i++) {
qqq << "ph_th_pip_5[" << i << "]";
ph_th_pip_5[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),220, 0., 110.,150., -35., 35.);
qqq.str("");
};

for (i=0; i<6; i++) {
qqq << "ph_th_pip_6[" << i << "]";
ph_th_pip_6[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),220, 0., 110.,150., -35., 35.);
qqq.str("");
};




for (i=0; i<6; i++) {
for(j=0; j<7; j++){
qqq << "ph_vs_th_el_sim_" << i+1 << "_" << j+1;
ph_vs_th_el_sim[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),240, 0., 60.,210., -35., 35.);
qqq.str("");
};
//for(j=0; j<5; j++){
qqq << "ph_vs_th_p_sim_" << i+1 << "_" << j+1;
ph_vs_th_p_sim[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),330, 0., 110.,210., -35., 35.);
qqq.str("");
//};
//for(j=0; j<6; j++){
qqq << "ph_vs_th_pip_sim_" << i+1 << "_" << j+1;
ph_vs_th_pip_sim[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),330, 0., 110.,210., -35., 35.);
qqq.str("");
//};
for(j=0; j<5; j++){
qqq << "ph_vs_th_pim_sim_" << i+1 << "_" << j+1;
ph_vs_th_pim_sim[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),360, 0., 120.,210., -35., 35.);
qqq.str("");

};
};




//beta_vs_p_pos[48][6]
for (i=0; i<6; i++) {
for(j=0; j<48; j++){
qqq << "beta_vs_p_p_" << i+1 << "_" << j+1;
beta_vs_p_p[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.0,300, 0., 1.5);
qqq.str("");
qqq << "beta_vs_p_pip_" << i+1 << "_" << j+1;
beta_vs_p_pip[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.0,300, 0., 1.5);
qqq.str("");
qqq << "beta_vs_p_pim_" << i+1 << "_" << j+1;
beta_vs_p_pim[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.0,300, 0., 1.5);
qqq.str("");


qqq << "time_p_" << i+1 << "_" << j+1;
time_p[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.0,600, -30., 30.);
qqq.str("");
qqq << "time_pip_" << i+1 << "_" << j+1;
time_pip[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.0,600, -30., 30.);
qqq.str("");
qqq << "time_pim_" << i+1 << "_" << j+1;
time_pim[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.0,600, -30., 30.);
qqq.str("");


qqq << "beta_vs_p_p_sim_" << i+1 << "_" << j+1;
beta_vs_p_p_sim[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.0,240, 0., 1.2);
qqq.str("");
qqq << "beta_vs_p_pip_sim_" << i+1 << "_" << j+1;
beta_vs_p_pip_sim[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.0,240, 0., 1.2);
qqq.str("");
qqq << "beta_vs_p_pim_sim_" << i+1 << "_" << j+1;
beta_vs_p_pim_sim[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.0,240, 0., 1.2);
qqq.str("");


qqq << "time_p_sim_" << i+1 << "_" << j+1;
time_p_sim[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.0,600, -30., 30.);
qqq.str("");
qqq << "time_pip_sim_" << i+1 << "_" << j+1;
time_pip_sim[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.0,600, -30., 30.);
qqq.str("");
qqq << "time_pim_sim_" << i+1 << "_" << j+1;
time_pim_sim[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.0,600, -30., 30.);
qqq.str("");







};
};



//TH2F *ph_vs_th_pim[6][10]
//TH2F *ph_th_pim_all_p[6];
for (i=0; i<6; i++) {
qqq << "ph_th_pim_all_p[" << i+1;
ph_th_pim_all_p[i]=new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0., 120.,150., -35., 35.);
qqq.str("");
for(j=0; j<10; j++){
qqq << "ph_vs_th_pim[" << i+1 << "][" << j+1 <<"]" ;
ph_vs_th_pim[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0., 120.,150., -35., 35.);
qqq.str("");
};
};


/*for (i=0; i<6; i++) {
for(j=0; j<57; j++){
qqq << "y_sc_vs_x_sc_p" << i+1 << "_" << j+1;
y_sc_vs_x_sc_p[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),600, -300.,300.,500, 0.,500.);
qqq.str("");
qqq << "y_sc_vs_x_sc_pim" << i+1 << "_" << j+1;
y_sc_vs_x_sc_pim[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),600, -300.,300.,500, 0.,500.);
qqq.str("");
qqq << "y_sc_vs_x_sc_pip" << i+1 << "_" << j+1;
y_sc_vs_x_sc_pip[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),600, -300.,300.,500, 0.,500.);
qqq.str("");
qqq << "y_sc_vs_x_sc_el" << i+1 << "_" << j+1;
y_sc_vs_x_sc_el[i][j] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),600, -300.,300.,500, 0.,500.);
qqq.str("");

};
};*/

 //Int_t bins[5] = {12, 12, 10, 8, 8};
  
  
 /* Double_t W_bin [20];
  for (i=0; i<20; i++) W_bin[i] = 1.325+0.025*i;
  
   Double_t xmin[5] = {1., 0.288,0.,0.,0.};
   Double_t xmax_1[5];
   Double_t xmax_1[5];
   for (i=0; i<20; i++) {
   xmax_1[i] = W_bin[i] - 0.13957 - 0.05;
   xmax_2[i] = W_bin[i] - 0.938272 - 0.05;
   };
   Double_t xmax[5] = {1.485,0.69, 180.,360.,360.};*/
   
   





/*
for (i=0; i<9; i++) {


W_bin[i] = 1.425+0.025*i;//right edge
xmax_1[i] = W_bin[i] - 0.13957 + 0.05;
xmax_2[i] = W_bin[i] - 0.938272 + 0.05;

xmax[0] = xmax_1[i]; 
xmax[1] = xmax_2[i];

bins[0]=12;
bins[1]=12;

qqq << "h_5dim_excl_1"<< "_w_" << 1000*(1.4125+0.025*i);
h_5dim_excl_1[i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_excl_2"<<"_w_" << 1000*(1.4125+0.025*i);
h_5dim_excl_2[i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_excl_3"<<"_w_" << 1000*(1.4125+0.025*i);
h_5dim_excl_3[i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");


qqq << "h_5dim_excl_1_empty" << "_w_" << 1000*(1.4125+0.025*i);
h_5dim_excl_1_empty[i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_excl_2_empty" << "_w_" << 1000*(1.4125+0.025*i);
h_5dim_excl_2_empty[i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_excl_3_empty" << "_w_" << 1000*(1.4125+0.025*i);
h_5dim_excl_3_empty[i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");


qqq << "h_5dim_1_sim_excl_1" << "_w_" << 1000*(1.4125+0.025*i);
h_5dim_1_sim_excl_1[i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_2_sim_excl_1"<< "_w_" << 1000*(1.4125+0.025*i);
h_5dim_2_sim_excl_1[i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_3_sim_excl_1" << "_w_" << 1000*(1.4125+0.025*i);
h_5dim_3_sim_excl_1[i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");


qqq << "h_5dim_1_sim_excl_2"<< "_w_" << 1000*(1.4125+0.025*i);
h_5dim_1_sim_excl_2[i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_2_sim_excl_2"<< "_w_" << 1000*(1.4125+0.025*i);
h_5dim_2_sim_excl_2[i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_3_sim_excl_2"<< "_w_" << 1000*(1.4125+0.025*i);
h_5dim_3_sim_excl_2[i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");


};


*/

xmax[2] = 180.;
xmax[3] = 360.;
xmax[4] = 360.;


for (i=0; i<21; i++) {


W_bin[i] = 1.3125+0.025*i;//center of the bin

if ((i>=0)&&(i<=1)){
bins[0]=8;
bins[1]=8;
bins[2]=6;
bins[3]=5;
bins[4]=5;
};
if ((i>=2)&&(i<=3)){
bins[0]=10;
bins[1]=10;
bins[2]=8;
bins[3]=5;
bins[4]=6;

};

if ((i>=4)&&(i<=6)){
bins[0]=12;
bins[1]=12;
bins[2]=10;
bins[3]=5;
bins[4]=8;
};

if ((i>=7)&&(i<=20)){
bins[0]=12;
bins[1]=12;
bins[2]=10;
bins[3]=6;
bins[4]=8;
};

//cout << i<<" bins0= " << bins[0]<<  " bins1= " << bins[1]<< " bins2= " << bins[2]<< " bins3= " << bins[3]<< " bins4= " << bins[4]<< "\n";

//xmax_1[i] = W_bin[i] - 0.13957 + 0.05;
//xmax_2[i] = W_bin[i] - 0.938272 + 0.05;

xmax_1[i] = W_bin[i] - 0.13957 + (W_bin[i]- 0.13957 - 0.13957 - 0.938272)/(bins[0]-1);
xmax_2[i] = W_bin[i] - 0.938272 +(W_bin[i]- 0.13957 - 0.13957 - 0.938272)/(bins[1]-1);


xmax[0] = xmax_1[i]; 
xmax[1] = xmax_2[i];





for(j=0; j<12; j++){



qqq << "h_5dim_excl_1_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 10000*(1.3125+0.025*i);
h_5dim_excl_1[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_excl_2_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 10000*(1.3125+0.025*i);
h_5dim_excl_2[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_excl_3_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 10000*(1.3125+0.025*i);
h_5dim_excl_3[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");



qqq << "h_5dim_excl_1_sim_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 10000*(1.3125+0.025*i);
h_5dim_excl_1_sim[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_excl_2_sim_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 10000*(1.3125+0.025*i);
h_5dim_excl_2_sim[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_excl_3_sim_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 10000*(1.3125+0.025*i);
h_5dim_excl_3_sim[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_excl_1_sim_evt_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 10000*(1.3125+0.025*i);
h_5dim_excl_1_sim_evt[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_excl_2_sim_evt_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 10000*(1.3125+0.025*i);
h_5dim_excl_2_sim_evt[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_excl_3_sim_evt_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 10000*(1.3125+0.025*i);
h_5dim_excl_3_sim_evt[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");


qqq << "h_5dim_pim_1_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 10000*(1.3125+0.025*i);
h_5dim_pim_1[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_pim_2_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 10000*(1.3125+0.025*i);
h_5dim_pim_2[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_pim_3_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 10000*(1.3125+0.025*i);
h_5dim_pim_3[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");



qqq << "h_5dim_pim_1_sim_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 10000*(1.3125+0.025*i);
h_5dim_pim_1_sim[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_pim_2_sim_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 10000*(1.3125+0.025*i);
h_5dim_pim_2_sim[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_pim_3_sim_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 10000*(1.3125+0.025*i);
h_5dim_pim_3_sim[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");


qqq << "h_5dim_pim_1_sim_evt_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 10000*(1.3125+0.025*i);
h_5dim_pim_1_sim_evt[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_pim_2_sim_evt_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 10000*(1.3125+0.025*i);
h_5dim_pim_2_sim_evt[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_pim_3_sim_evt_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 10000*(1.3125+0.025*i);
h_5dim_pim_3_sim_evt[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");


qqq << "h_5dim_1_sim_gen_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 10000*(1.3125+0.025*i);
h_5dim_1_sim_gen[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_2_sim_gen_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 10000*(1.3125+0.025*i);
h_5dim_2_sim_gen[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_3_sim_gen_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 10000*(1.3125+0.025*i);
h_5dim_3_sim_gen[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");


qqq << "h_5dim_1_sim_gen_evt_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 10000*(1.3125+0.025*i);
h_5dim_1_sim_gen_evt[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_2_sim_gen_evt_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 10000*(1.3125+0.025*i);
h_5dim_2_sim_gen_evt[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");

qqq << "h_5dim_3_sim_gen_evt_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 10000*(1.3125+0.025*i);
h_5dim_3_sim_gen_evt[j][i] = new THnSparseD(qqq.str().c_str(),qqq.str().c_str(),ndims,bins,xmin,xmax);
qqq.str("");



};
};

/*
for(j=0; j<12; j++){
qqq << "w_int_" << j+1;
h_w_int[j] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),21, 1.3,1.825);
qqq.str("");
};
*/
for (i=0; i<6; i++) {
for(j=0; j<20; j++){
qqq << "sector_" << i+1 << "_left_" << j+1 ;
ph_el_left[i][j] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),500, 0.,500);
qqq.str("");
};
};

for (i=0; i<6; i++) {
for(j=0; j<20; j++){
qqq << "sector_" << i+1 << "_both_" << j+1 ;
ph_el_both[i][j] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),500, 0.,500);
qqq.str("");
};
};

for (i=0; i<6; i++) {
for(j=0; j<20; j++){
qqq << "sector_" << i+1 << "_right_" << j+1 ;
ph_el_right[i][j] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),500, 0.,500);
qqq.str("");
};
};


for (i=0; i<6; i++) {
qqq << "th_vs_p_e_1_" << i+1;
th_vs_p_e_1[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),440, 0.,2.2,180,0.,60.);
qqq.str("");
qqq << "th_vs_p_e_2_" << i+1;
th_vs_p_e_2[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),440, 0.,2.2,180,0.,60.);
qqq.str("");
qqq << "th_vs_p_p_1_" << i+1;
th_vs_p_p_1[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.,360,0.,180.);
qqq.str("");
qqq << "th_vs_p_p_2_" << i+1;
th_vs_p_p_2[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.,360,0.,180.);
qqq.str("");
qqq << "th_vs_p_pip_1_" << i+1;
th_vs_p_pip_1[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.,360,0.,180.);
qqq.str("");
qqq << "th_vs_p_pip_2_" << i+1;
th_vs_p_pip_2[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.,360,0.,180.);
qqq.str("");
qqq << "th_vs_p_pim_1_" << i+1;
th_vs_p_pim_1[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.,360,0.,180.);
qqq.str("");
qqq << "th_vs_p_pim_2_" << i+1;
th_vs_p_pim_2[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.,360,0.,180.);
qqq.str("");

qqq << "th_vs_p_e_1_sim_" << i+1;
th_vs_p_e_1_sim[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),440, 0.,2.2,180,0.,60.);
qqq.str("");
qqq << "th_vs_p_e_2_sim_" << i+1;
th_vs_p_e_2_sim[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),440, 0.,2.2,180,0.,60.);
qqq.str("");
qqq << "th_vs_p_p_1_sim_" << i+1;
th_vs_p_p_1_sim[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.,360,0.,180.);
qqq.str("");
qqq << "th_vs_p_p_2_sim_" << i+1;
th_vs_p_p_2_sim[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.,360,0.,180.);
qqq.str("");
qqq << "th_vs_p_pip_1_sim_" << i+1;
th_vs_p_pip_1_sim[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.,360,0.,180.);
qqq.str("");
qqq << "th_vs_p_pip_2_sim_" << i+1;
th_vs_p_pip_2_sim[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.,360,0.,180.);
qqq.str("");
qqq << "th_vs_p_pim_1_sim_" << i+1;
th_vs_p_pim_1_sim[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.,360,0.,180.);
qqq.str("");
qqq << "th_vs_p_pim_2_sim_" << i+1;
th_vs_p_pim_2_sim[i] = new TH2F(qqq.str().c_str(),qqq.str().c_str(),400, 0.,2.,360,0.,180.);
qqq.str("");



};


};
