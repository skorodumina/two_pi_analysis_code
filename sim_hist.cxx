#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include <TLorentzVector.h>
#include "cuts_data.h"
#include "sim_hist.h"
#include "global.h"
#include <iostream>



using namespace std;


 void sim_hist() {
 Float_t beta_nom_pip,beta_nom_pim,beta_nom_p;
 beta_nom_pip = P_PIp/sqrt(m_pip*m_pip+P_PIp*P_PIp);
beta_nom_pim = P_PIm/sqrt(m_pip*m_pip+P_PIm*P_PIm);
beta_nom_p = P_P/sqrt(m_proton*m_proton+P_P*P_P);
 
 /*if (sector ==1) {
if (pmt_hit == -1) ph_el_left[0][segment]->Fill(nphe,1.);
if (pmt_hit == 0) ph_el_both[0][segment]->Fill(nphe,1.);
if (pmt_hit == 1) ph_el_right[0][segment]->Fill(nphe,1.);
 };
if (sector ==2)  {
if (pmt_hit == -1) ph_el_left[1][segment]->Fill(nphe,1.);
if (pmt_hit == 0) ph_el_both[1][segment]->Fill(nphe,1.);
if (pmt_hit == 1) ph_el_right[1][segment]->Fill(nphe,1.);

};

if (sector ==3) {
 if (pmt_hit == -1) ph_el_left[2][segment]->Fill(nphe,1.);
if (pmt_hit == 0) ph_el_both[2][segment]->Fill(nphe,1.);
if (pmt_hit == 1) ph_el_right[2][segment]->Fill(nphe,1.);
 
 };
if (sector ==4)  {
if (pmt_hit == -1) ph_el_left[3][segment]->Fill(nphe,1.);
if (pmt_hit == 0) ph_el_both[3][segment]->Fill(nphe,1.);
if (pmt_hit == 1) ph_el_right[3][segment]->Fill(nphe,1.);

};
if (sector ==5) {
if (pmt_hit == -1) ph_el_left[4][segment]->Fill(nphe,1.);
if (pmt_hit == 0) ph_el_both[4][segment]->Fill(nphe,1.);
if (pmt_hit == 1) ph_el_right[4][segment]->Fill(nphe,1.);

 };
if (sector ==6) {
if (pmt_hit == -1) ph_el_left[5][segment]->Fill(nphe,1.);
if (pmt_hit == 0) ph_el_both[5][segment]->Fill(nphe,1.);
if (pmt_hit == 1) ph_el_right[5][segment]->Fill(nphe,1.);

 };*/



/*
if (sector ==1) {
h_cc_nphe_total_s1->Fill(theta_cc,ph_cc,1.);
if (nphe > 70) h_cc_nphe_final_s1->Fill(theta_cc,ph_cc,1.);
 };
if (sector ==2)  {
h_cc_nphe_total_s2->Fill(theta_cc,ph_cc,1.);
if (nphe > 70) h_cc_nphe_final_s2->Fill(theta_cc,ph_cc,1.);

};

if (sector ==3) {
h_cc_nphe_total_s3->Fill(theta_cc,ph_cc,1.);
if (nphe > 70) h_cc_nphe_final_s3->Fill(theta_cc,ph_cc,1.); 
 };
if (sector ==4)  {
h_cc_nphe_total_s4->Fill(theta_cc,ph_cc,1.);
if (nphe > 70) h_cc_nphe_final_s4->Fill(theta_cc,ph_cc,1.);
};
if (sector ==5) {
h_cc_nphe_total_s5->Fill(theta_cc,ph_cc,1.);
if (nphe > 70) h_cc_nphe_final_s5->Fill(theta_cc,ph_cc,1.);
 };
if (sector ==6) {
h_cc_nphe_total_s6->Fill(theta_cc,ph_cc,1.);
if (nphe > 70) h_cc_nphe_final_s6->Fill(theta_cc,ph_cc,1.);
 };
*/


//cout << P_dist << "  "<< beta_nom_p << " "<< beta_P<<"\n";


if ((ph_P >= 330)&& (ph_P <= 360)) {

if (indtype==1) beta_vs_p_p_sim[0][PdHit_P-1]->Fill(P_P,beta_P,1);
if (indtype==1) time_p_sim[0][PdHit_P-1]->Fill(P_P,P_dist/30.*(1./beta_nom_p-1./beta_P),1.);
if (indtype==1) th_vs_p_p_1_sim[0] -> Fill(P_P,th_P,1);
 };//end of the first part of the first sector
 
 if ((ph_P >= 0) && (ph_P <= 30)) {

if (indtype==1) beta_vs_p_p_sim[0][PdHit_P-1]->Fill(P_P,beta_P,1);
if (indtype==1) time_p_sim[0][PdHit_P-1]->Fill(P_P,P_dist/30.*(1./beta_nom_p-1./beta_P),1.);
if (indtype==1) th_vs_p_p_1_sim[0] -> Fill(P_P,th_P,1);
 };//end of the second part of the first sector
  
 if ((ph_P >= 30) && (ph_P <=90)){ 

if (indtype==1) beta_vs_p_p_sim[1][PdHit_P-1]->Fill(P_P,beta_P,1);
if (indtype==1) time_p_sim[1][PdHit_P-1]->Fill(P_P,P_dist/30.*(1./beta_nom_p-1./beta_P),1.);
if (indtype==1) th_vs_p_p_1_sim[1] -> Fill(P_P,th_P,1);

 };//end of the sector2
 
 if ((ph_P >=90) && (ph_P <=150)) {

if (indtype==1) beta_vs_p_p_sim[2][PdHit_P-1]->Fill(P_P,beta_P,1);
if (indtype==1) time_p_sim[2][PdHit_P-1]->Fill(P_P,P_dist/30.*(1./beta_nom_p-1./beta_P),1.);
if (indtype==1) th_vs_p_p_1_sim[2] -> Fill(P_P,th_P,1);
 };//end of the sector3
 
 if ((ph_P >= 150) && (ph_P <= 210)) {

if (indtype==1) beta_vs_p_p_sim[3][PdHit_P-1]->Fill(P_P,beta_P,1);
if (indtype==1) time_p_sim[3][PdHit_P-1]->Fill(P_P,P_dist/30.*(1./beta_nom_p-1./beta_P),1.);
if (indtype==1) th_vs_p_p_1_sim[3] -> Fill(P_P,th_P,1);
};//end of the sector4

if ((ph_P >= 210) && (ph_P <=270)) {

if (indtype==1) beta_vs_p_p_sim[4][PdHit_P-1]->Fill(P_P,beta_P,1);
if (indtype==1) time_p_sim[4][PdHit_P-1]->Fill(P_P,P_dist/30.*(1./beta_nom_p-1./beta_P),1.);
if (indtype==1) th_vs_p_p_1_sim[4] -> Fill(P_P,th_P,1);
 };//end of the sector5
 
 if ((ph_P >= 270) && (ph_P <=330)) {
 
if (indtype==1) beta_vs_p_p_sim[5][PdHit_P-1]->Fill(P_P,beta_P,1);
if (indtype==1) time_p_sim[5][PdHit_P-1]->Fill(P_P,P_dist/30.*(1./beta_nom_p-1./beta_P),1.);
if (indtype==1) th_vs_p_p_1_sim[5] -> Fill(P_P,th_P,1);
 };//end of the sector6


//-------------------------------------
if ((ph_PIp >= 330)&& (ph_PIp <= 360)) {

  if (indtype==1) beta_vs_p_pip_sim[0][PdHit_PIp-1]->Fill(P_PIp,beta_PIp,1); 
  if (indtype==1) time_pip_sim[0][PdHit_PIp-1]->Fill(P_PIp,PIp_dist*(1./beta_nom_pip-1/beta_PIp)/30.,1.);
  if (indtype==1) th_vs_p_pip_1_sim[0]->Fill(P_PIp,th_PIp,1.);
 };//end of the first part of the first sector
 
 if ((ph_PIp >= 0) && (ph_PIp <= 30)) {

  if (indtype==1) beta_vs_p_pip_sim[0][PdHit_PIp-1]->Fill(P_PIp,beta_PIp,1); 
  if (indtype==1) time_pip_sim[0][PdHit_PIp-1]->Fill(P_PIp,PIp_dist*(1./beta_nom_pip-1/beta_PIp)/30.,1.);
  if (indtype==1) th_vs_p_pip_1_sim[0]->Fill(P_PIp,th_PIp,1.);
 };//end of the second part of the first sector
  
 if ((ph_PIp >= 30) && (ph_PIp <=90)){ 

 if (indtype==1) beta_vs_p_pip_sim[1][PdHit_PIp-1]->Fill(P_PIp,beta_PIp,1);
  if (indtype==1) time_pip_sim[1][PdHit_PIp-1]->Fill(P_PIp,PIp_dist*(1./beta_nom_pip-1/beta_PIp)/30.,1.);
  if (indtype==1) th_vs_p_pip_1_sim[1]->Fill(P_PIp,th_PIp,1.);
 };//end of the sector2
 
 if ((ph_PIp >=90) && (ph_PIp <=150)) {

 if (indtype==1) beta_vs_p_pip_sim[2][PdHit_PIp-1]->Fill(P_PIp,beta_PIp,1);
   if (indtype==1) time_pip_sim[2][PdHit_PIp-1]->Fill(P_PIp,PIp_dist*(1./beta_nom_pip-1/beta_PIp)/30.,1.);
   if (indtype==1) th_vs_p_pip_1_sim[2]->Fill(P_PIp,th_PIp,1.);
 };//end of the sector3
 
 if ((ph_PIp >= 150) && (ph_PIp <= 210)) {

 if (indtype==1) beta_vs_p_pip_sim[3][PdHit_PIp-1]->Fill(P_PIp,beta_PIp,1);
 if (indtype==1) time_pip_sim[3][PdHit_PIp-1]->Fill(P_PIp,PIp_dist*(1./beta_nom_pip-1/beta_PIp)/30.,1.);
 if (indtype==1) th_vs_p_pip_1_sim[3]->Fill(P_PIp,th_PIp,1.);
};//end of the sector4

if ((ph_PIp >= 210) && (ph_PIp <=270)) {

if (indtype==1) beta_vs_p_pip_sim[4][PdHit_PIp-1]->Fill(P_PIp,beta_PIp,1);
if (indtype==1) time_pip_sim[4][PdHit_PIp-1]->Fill(P_PIp,PIp_dist*(1./beta_nom_pip-1/beta_PIp)/30.,1.);
if (indtype==1) th_vs_p_pip_1_sim[4]->Fill(P_PIp,th_PIp,1.);
 };//end of the sector5
 
 if ((ph_PIp >= 270) && (ph_PIp <=330)) {
 

if (indtype==1) beta_vs_p_pip_sim[5][PdHit_PIp-1]->Fill(P_PIp,beta_PIp,1);
if (indtype==1) time_pip_sim[5][PdHit_PIp-1]->Fill(P_PIp,PIp_dist*(1./beta_nom_pip-1/beta_PIp)/30.,1.);
if (indtype==1) th_vs_p_pip_1_sim[5]->Fill(P_PIp,th_PIp,1.);
 };//end of the sector6

//-----------------------------------
if ((ph_PIm >= 330)&& (ph_PIm <= 360)) {

 if (indtype==1) beta_vs_p_pim_sim[0][PdHit_PIm-1]->Fill(P_PIm,beta_PIm,1);
 if (indtype==1) time_pim_sim[0][PdHit_PIm-1]->Fill(P_PIm,PIm_dist/30.*(1./beta_nom_pim-1./beta_PIm),1.);
 if (indtype==1) th_vs_p_pim_1_sim[0]->Fill(P_PIm,th_PIm,1.);
 };//end of the first part of the first sector
 
 if ((ph_PIm >= 0) && (ph_PIm <= 30)) {

 if (indtype==1) beta_vs_p_pim_sim[0][PdHit_PIm-1]->Fill(P_PIm,beta_PIm,1);
 if (indtype==1) time_pim_sim[0][PdHit_PIm-1]->Fill(P_PIm,PIm_dist/30.*(1./beta_nom_pim-1./beta_PIm),1.);
 if (indtype==1) th_vs_p_pim_1_sim[0]->Fill(P_PIm,th_PIm,1.);
 };//end of the second part of the first sector
  
 if ((ph_PIm >= 30) && (ph_PIm <=90)){ 

if (indtype==1) beta_vs_p_pim_sim[1][PdHit_PIm-1]->Fill(P_PIm,beta_PIm,1);
if (indtype==1) time_pim_sim[1][PdHit_PIm-1]->Fill(P_PIm,PIm_dist/30.*(1./beta_nom_pim-1./beta_PIm),1.);
if (indtype==1) th_vs_p_pim_1_sim[1]->Fill(P_PIm,th_PIm,1.);
 };//end of the sector2
 
 if ((ph_PIm >=90) && (ph_PIm <=150)) {

if (indtype==1) beta_vs_p_pim_sim[2][PdHit_PIm-1]->Fill(P_PIm,beta_PIm,1);
if (indtype==1) time_pim_sim[2][PdHit_PIm-1]->Fill(P_PIm,PIm_dist/30.*(1./beta_nom_pim-1./beta_PIm),1.);
if (indtype==1) th_vs_p_pim_1_sim[2]->Fill(P_PIm,th_PIm,1.);
 };//end of the sector3
 
 if ((ph_PIm >= 150) && (ph_PIm <= 210)) {

if (indtype==1) beta_vs_p_pim_sim[3][PdHit_PIm-1]->Fill(P_PIm,beta_PIm,1);
if (indtype==1) time_pim_sim[3][PdHit_PIm-1]->Fill(P_PIm,PIm_dist/30.*(1./beta_nom_pim-1./beta_PIm),1.);
if (indtype==1) th_vs_p_pim_1_sim[3]->Fill(P_PIm,th_PIm,1.);
};//end of the sector4

if ((ph_PIm >= 210) && (ph_PIm <=270)) {

if (indtype==1) beta_vs_p_pim_sim[4][PdHit_PIm-1]->Fill(P_PIm,beta_PIm,1);
if (indtype==1) time_pim_sim[4][PdHit_PIm-1]->Fill(P_PIm,PIm_dist/30.*(1./beta_nom_pim-1./beta_PIm),1.);
if (indtype==1) th_vs_p_pim_1_sim[4]->Fill(P_PIm,th_PIm,1.);
 };//end of the sector5
 
 if ((ph_PIm >= 270) && (ph_PIm <=330)) {
 
if (indtype==1) beta_vs_p_pim_sim[5][PdHit_PIm-1]->Fill(P_PIm,beta_PIm,1);
if (indtype==1) time_pim_sim[5][PdHit_PIm-1]->Fill(P_PIm,PIm_dist/30.*(1./beta_nom_pim-1./beta_PIm),1.);
if (indtype==1) th_vs_p_pim_1_sim[5]->Fill(P_PIm,th_PIm,1.);
 };//end of the sector6





/*if (sector ==1) {
 hist_z_el_1->Fill(z_EL,1.);
th_vs_p_e_1[0]->Fill(P_EL,th_EL,1.);
 };
if (sector ==2)  {

th_vs_p_e_1[1]->Fill(P_EL,th_EL,1.);
hist_z_el_2->Fill(z_EL,1.);
};

if (sector ==3) {
 hist_z_el_3->Fill(z_EL,1.);
th_vs_p_e_1[2]->Fill(P_EL,th_EL,1.);
 };
if (sector ==4)  {
th_vs_p_e_1[3]->Fill(P_EL,th_EL,1.);
hist_z_el_4->Fill(z_EL,1.);
};
if (sector ==5) {
th_vs_p_e_1[4]->Fill(P_EL,th_EL,1.);
 hist_z_el_5->Fill(z_EL,1.);
 };
if (sector ==6) {
th_vs_p_e_1[5]->Fill(P_EL,th_EL,1.);
 hist_z_el_6->Fill(z_EL,1.);
 };*/
 

 };
