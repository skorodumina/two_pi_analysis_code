#include <TMath.h>
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include "cuts_sim.h"
#include "global.h"
#include <iostream>



using namespace std;



 bool cuts_sim::Electron_cuts_sim(){
   
   
   

   
   bool cuts_sim;
   Float_t th_min,th_max,par1,par2,par3,fid_a,fid_b;
   Short_t i;
  
  
   Float_t ph_el_arr_sim[3][6][18] = {{{1000.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1000.},
                             {1000.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1000.},
		             {1000.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1000.},
	                     {1000.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1000.},
		             {1000.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1000.},
			     {1000.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1000.}},
			     
			     {{1000.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1000.},
                             {1000.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1000.},
		             {1000.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1000.},
	                     {1000.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1000.},
		             {1000.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1000.},
			     {1000.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1000.}},
			     
			     {{1000.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1000.},
                             {1000.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1000.},
		             {1000.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1000.},
	                     {1000.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1000.},
		             {1000.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1000.},
			     {1000.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1000.}}};

	
Float_t th_vs_seg_cc_arr[2][6][18] = {{{0.,13.0577,14.8198,16.9854,18.8497,20.7722, 23.1202, 24.7399, 26.3693, 28.2569, 31.6481, 33.7044, 37.5222, 40.1851, 42.4846, 45.,45.,0.},
                             {0., 12.3663, 14.5647, 16.3118, 18.7031, 20.8298, 22.5443, 24.8247, 26.3913, 28.6392, 31.5305, 33.566, 37.5775, 40.1132, 42.5322, 45., 45., 0.},
                             {0., 12.76, 14.6586, 16.6821, 18.6648, 20.5538, 22.5573, 24.962, 26.4137, 28.3029, 31.6033, 33.7665, 37.4934, 40.0381, 42.3285, 44.8663, 45.,0.},
                             {0., 13.0839, 14.7751, 16.7626, 18.6895, 20.9305, 23.1239, 24.7276, 26.5469, 28.1259, 31.2627, 33.6781, 37.3355, 39.6646, 42.5703, 45., 45., 0.},
                             {0., 12.6631, 16.052, 17.1266, 18.7471, 20.7836, 23.0648, 25.0117, 27.5244, 29.9016, 32.6042, 34.7115, 36.9986, 39.1818, 42.0926, 45., 45., 0.},
                             {0., 12.7948, 14.7214, 16.7988, 18.5596, 20.5496, 22.8541, 24.4748, 26.3951, 28.4766, 31.5388, 33.7486, 37.3584, 40.1278, 42.5831, 45., 45., 0.}},
			     
			     
			    {{0., 11., 11., 11.6061, 13.2822, 15.2113, 16.8047, 19.1317, 21.1941, 22.8422, 23.7821, 26.5811, 27.7998, 31.2606, 34.2189, 36.2558, 39.5832, 0.},
                             {0., 11., 11., 11.6132, 12.6424, 14.3677, 16.8573, 18.49, 20.9482, 22.5705, 23.7656, 26.4239, 27.6032, 30.7262, 33.8834, 36.3639, 38.9472, 0.},
		             {0., 11., 11., 11.643, 13.0375, 14.9574, 16.9836, 18.6356, 21.1218, 22.8459, 23.7062, 26.2906, 27.8548, 31.3241, 34.3427, 36.7193, 39.2374, 0.},
	                     {0., 11., 11., 11.7798, 13.2277, 14.8726, 16.6899, 19.1786, 21.0299, 22.9571, 23.9678, 26.1826, 27.4338, 31.2199, 33.6087, 36.1697, 39.178, 0.},
		             {0., 11., 11., 12.0142, 13.4243, 15.159, 16.8523, 19.0246, 20.7604, 22.5512, 24.7555, 26.934, 29.5094, 30.8662, 35.0397, 36.5057, 39.7693, 0.},
			     {0., 11., 11., 11.628, 13.1929, 15.0126, 16.3127, 19.0631, 20.872, 22.7261, 23.8125, 26.4752, 27.82, 30.8766, 34.2017, 36.4635, 40.1278, 0.}}};
			     
 cuts_sim = false; 

//calorimeter threshold cut + manually remove the first and last cc segments
if ((P_EL > 0.461)&&(segment!=0)&&(segment!=17)) {
     
   th_min= (11.7398+8.21504/(0.433327*P_EL+0.158076));
   par1=0.85+1.1*P_EL;
   par2=-62.8-30.*P_EL; 
   par3 = 0.0047*P_EL + 0.0079;

   fid_a = 41.3*pow((sin((th_EL-th_min)*par3)),(par1+par2/th_EL+1485./th_EL/th_EL))+1.; 
   fid_b = -41.3*pow((sin((th_EL-th_min)*par3)),(par1+par2/th_EL+1485./th_EL/th_EL))-1.; 
     
   th_max = 76.8617 -76.537*P_EL + 77.9387*P_EL*P_EL-28.389*P_EL*P_EL*P_EL;
  
   

	switch (sector) {
case 1 : 
if (indtype ==2) hist_z_el_1_sim_2->Fill(z_EL,1.);
/*if ((ph_EL >= 330) && (ph_EL <= 360)){
if ((P_EL < 1.79999) && (P_EL > 0.4)){
if (indtype==2) ph_vs_th_el_sim[0][int((P_EL*100-40)/20)]->Fill(th_EL,ph_EL-360,1);
};
};
if ((ph_EL >= 0) && (ph_EL <= 30)) {
if ((P_EL < 1.75999) && (P_EL > 0.4)){
if (indtype==2) ph_vs_th_el_sim[0][int((P_EL*100-40)/20)]->Fill(th_EL,ph_EL,1);
};
};*/

// ectot vs p cut
if ((ECT/P_EL > -0.0724295*P_EL*P_EL+0.213588*P_EL-0.0148421) && (ECT/P_EL < 0.0136346*P_EL*P_EL-0.0537425*P_EL+0.360596)){

	if (indtype==1) hist_ectot_sector1_sim->Fill(P_EL,ECT/P_EL,sigma); 

 //~fid cuts in Cherenkov plane
if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {

//phi_cc matching
if (((ph_cc>0)&&(pmt_hit==1))||((ph_cc<0)&&(pmt_hit==-1))||(pmt_hit==0)){

//th_cc vs seg cut
if ((theta_cc >th_vs_seg_cc_arr[1][0][segment])&&(theta_cc <th_vs_seg_cc_arr[0][0][segment])){

	if (indtype==1) th_cc_vs_seg_1_sim->Fill(segment+1,theta_cc,sigma);

// geometrical cut on number of photoelectrons
if (norm_nphe_s1->GetBinContent(norm_nphe_s1->GetXaxis()->FindBin(theta_cc),norm_nphe_s1->GetYaxis()->FindBin(ph_cc))>0.7){

//nphe cut that removes the 1 and 18 segments
if (nphe > ph_el_arr_sim[pmt_hit+1][0][segment]){

//vertex cut
if ((z_EL>-2.65) && (z_EL<1.85)){
	
if ((ph_EL >= 330) && (ph_EL <= 360)){

//fiducial cut		
if ((th_EL > th_min) && (th_EL < th_max)&& (ph_EL > fid_b+360) && (ph_EL < fid_a+360)){

	if ((P_EL < 1.79999) && (P_EL > 0.4)){
	if (indtype==1) ph_vs_th_el_sim[0][int((P_EL*100-40)/20)]->Fill(th_EL,ph_EL-360,1);
	};

	if ((W>1.3)&&(W<1.825)&&(Q2>0.4)&&(Q2<1.)&&(n_PIp==1)&&(indtype==1)) th_vs_p_e_2_sim[0]->Fill(P_EL,th_EL,sigma);
//	 hist_z_el_1_sim_1->Fill(z_EL,sigma);

   cuts_sim = true;
   
 }; //fiducial
 }; //second part of sector 1
   
if ((ph_EL >= 0) && (ph_EL <= 30)) {

//fiducial cut
if ((th_EL > th_min) && (th_EL < th_max)&& (ph_EL > fid_b) && (ph_EL < fid_a)){

	if ((P_EL < 1.75999) && (P_EL > 0.4)){
	if (indtype==1) ph_vs_th_el_sim[0][int((P_EL*100-40)/20)]->Fill(th_EL,ph_EL,1);
	};

	if ((W>1.3)&&(W<1.825)&&(Q2>0.4)&&(Q2<1.)&&(n_PIp==1)&&(indtype==1)) th_vs_p_e_2_sim[0]->Fill(P_EL,th_EL,sigma); 
//	if ((indtype==1)) hist_z_el_1_sim_1->Fill(z_EL,sigma);
	
   cuts_sim = true; 

 }; //fiducial
 }; //first part of sector 1

 };//vertex cut
   };//nphe cut that removes the 1 and 18 segments
   };//geometrical cut on number of photoelectrons
   };//th_cc vs seg cut   
   };//phi_cc matching
   };//~fid cuts in Cherenkov plane
   };//~fid cuts in Cherenkov plane
   };//~fid cuts in Cherenkov plane
   };//~fid cuts in Cherenkov plane
   
 }; // ectot vs p cut
 
 break;
 
case 2 : 
if (indtype ==2) hist_z_el_2_sim_2->Fill(z_EL,1.); 
//if ((P_EL < 1.75999) && (P_EL > 0.4)){
//if (indtype==2) ph_vs_th_el_sim[1][int((P_EL*100-40)/20)]->Fill(th_EL,ph_EL-60,1);
//};

// ectot vs p cut
if ((ECT/P_EL > -0.0514781*P_EL*P_EL+0.170698*P_EL+0.00710391) && (ECT/P_EL < 0.012651*P_EL*P_EL-0.0528591*P_EL+0.35928)) {

	if (indtype==1) hist_ectot_sector2_sim->Fill(P_EL,ECT/P_EL,sigma);
	  
//~fid cuts in Cherenkov plane
if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {

//phi_cc matching
if (((ph_cc>0)&&(pmt_hit==1))||((ph_cc<0)&&(pmt_hit==-1))||(pmt_hit==0)){

//th_cc vs seg cut
if ((theta_cc >th_vs_seg_cc_arr[1][1][segment])&&(theta_cc <th_vs_seg_cc_arr[0][1][segment])){	

	if (indtype==1) th_cc_vs_seg_2_sim->Fill(segment+1,theta_cc,sigma);

// geometrical cut on number of photoelectrons
if (norm_nphe_s2->GetBinContent(norm_nphe_s2->GetXaxis()->FindBin(theta_cc),norm_nphe_s2->GetYaxis()->FindBin(ph_cc))>0.65){

//nphe cut that removes the 1 and 18 segments
if (nphe > ph_el_arr_sim[pmt_hit+1][1][segment]){


//fiducial cut
if ((th_EL > th_min) && (th_EL < th_max)&& (ph_EL > fid_b+60) && (ph_EL < fid_a+60)&&(PdHit_EL!=16)){

	if ((P_EL < 1.75999) && (P_EL > 0.4)){
	if (indtype==1) ph_vs_th_el_sim[1][int((P_EL*100-40)/20)]->Fill(th_EL,ph_EL-60,1);
	};

	if ((W>1.3)&&(W<1.825)&&(Q2>0.4)&&(Q2<1.)&&(n_PIp==1)&&(indtype==1)) th_vs_p_e_2_sim[1]->Fill(P_EL,th_EL,sigma);

//th_vs_p
if ((th_EL > (11.7398+8.21504/(0.433327*(P_EL+0.1)+0.158076))+18.3)||(th_EL < (11.7398+8.21504/(0.433327*(P_EL+0.1)+0.158076))+16.)){
//if ((th_EL > (11.7398+8.21504/(0.433327*(P_EL+0.1)+0.125)) +3.05)||(th_EL < (11.7398+8.21504/(0.433327*(P_EL+0.1)+0.14)) +2.3)){

//if (th_EL > (11.7398+8.5/(0.433327*(P_EL+0.1)+0.15)) +3.3){
//if ((th_EL > (11.7398+8.5/(0.433327*(P_EL+0.1)+0.15)) +6.)||(th_EL < (11.7398+8.21504/(0.433327*(P_EL+0.1)+0.14)) +4.7)){


//vertex cut
if ((z_EL>-2.35) && (z_EL<2.15)){

//	if ((indtype==1)) hist_z_el_2_sim_1->Fill(z_EL,sigma);

   cuts_sim = true; 
   
   };//vertex cut
   };//th_vs_p
 //  };//th_vs_p   
 //  };//th_vs_p   
   };//fiducial 
   };//nphe cut that removes the 1 and 18 segments
   };//geometrical cut on number of photoelectrons
   };//th_cc vs seg cut   
   };//phi_cc matching
   };//~fid cuts in Cherenkov plane
   };//~fid cuts in Cherenkov plane
   };//~fid cuts in Cherenkov plane
   };//~fid cuts in Cherenkov plane

 }; // ectot vs p cut

 break; 
 
case 3 : 
  
if (indtype ==2) hist_z_el_3_sim_2->Fill(z_EL,1.); 
//if ((P_EL < 1.75999) && (P_EL > 0.4)){
//if (indtype==2) ph_vs_th_el_sim[2][int((P_EL*100-40)/20)]->Fill(th_EL,ph_EL-120,1);
//};

// ectot vs p cut
if ((ECT/P_EL > -0.0659335*P_EL*P_EL+0.199308*P_EL-0.00592097) && (ECT/P_EL < 0.0415148*P_EL*P_EL-0.111881*P_EL+0.387163)) { 

	if (indtype==1) hist_ectot_sector3_sim->Fill(P_EL,ECT/P_EL,sigma);
	
//~fid cuts in Cherenkov plane 
if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {  

//phi_cc matching
if (((ph_cc>0)&&(pmt_hit==1))||((ph_cc<0)&&(pmt_hit==-1))||(pmt_hit==0)){

//th_cc vs seg cut
if ((theta_cc >th_vs_seg_cc_arr[1][2][segment])&&(theta_cc <th_vs_seg_cc_arr[0][2][segment])){	

	if (indtype==1) th_cc_vs_seg_3_sim->Fill(segment+1,theta_cc,sigma);

// geometrical cut on number of photoelectrons
if (norm_nphe_s3->GetBinContent(norm_nphe_s3->GetXaxis()->FindBin(theta_cc),norm_nphe_s3->GetYaxis()->FindBin(ph_cc))>0.7){

//nphe cut that removes the 1 and 18 segments
if (nphe > ph_el_arr_sim[pmt_hit+1][2][segment]){

//fiducial cut
if ((th_EL > th_min) && (th_EL < th_max)&& (ph_EL > fid_b+120) && (ph_EL < fid_a+120)&&(PdHit_EL!=44)){

	if ((P_EL < 1.75999) && (P_EL > 0.4)){
	if (indtype==1) ph_vs_th_el_sim[2][int((P_EL*100-40)/20)]->Fill(th_EL,ph_EL-120,1);
	};

	if ((W>1.3)&&(W<1.825)&&(Q2>0.4)&&(Q2<1.)&&(n_PIp==1)&&(indtype==1)) th_vs_p_e_2_sim[2]->Fill(P_EL,th_EL,sigma);
 
 
//vertex cut
if ((z_EL>-2.25) && (z_EL< 2.25)){ 

//	if ((indtype==1)) hist_z_el_3_sim_1->Fill(z_EL,sigma); 
 
   cuts_sim = true; 

   };//vertex cut
   };//fiducial 
   };//nphe cut that removes the 1 and 18 segments
   };//geometrical cut on number of photoelectrons
   };//th_cc vs seg cut   
   };//phi_cc matching
   };//~fid cuts in Cherenkov plane 
   };//~fid cuts in Cherenkov plane 
   };//~fid cuts in Cherenkov plane 
   };//~fid cuts in Cherenkov plane 
   
 }; // ectot vs p cut

 break;  
 
case 4 : 

//if ((P_EL < 1.75999) && (P_EL > 0.4)){
//if (indtype==2) ph_vs_th_el_sim[3][int((P_EL*100-40)/20)]->Fill(th_EL,ph_EL-180,1);
//};
if (indtype ==2) hist_z_el_4_sim_2->Fill(z_EL,1.); 

// ectot vs p cut  
if ((ECT/P_EL > -0.00634486*P_EL*P_EL+0.0672081*P_EL+0.061873) && (ECT/P_EL < -0.0504573*P_EL*P_EL+0.0867658*P_EL+0.287867)) {

	if (indtype==1) hist_ectot_sector4_sim->Fill(P_EL,ECT/P_EL,sigma); 

//~fid cuts in Cherenkov plane
if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {

//phi_cc matching
if (((ph_cc>0)&&(pmt_hit==1))||((ph_cc<0)&&(pmt_hit==-1))||(pmt_hit==0)){

//th_cc vs seg cut
if ((theta_cc >th_vs_seg_cc_arr[1][3][segment])&&(theta_cc <th_vs_seg_cc_arr[0][3][segment])){

	if (indtype==1) th_cc_vs_seg_4_sim->Fill(segment+1,theta_cc,sigma);

// geometrical cut on number of photoelectrons
if (norm_nphe_s4->GetBinContent(norm_nphe_s4->GetXaxis()->FindBin(theta_cc),norm_nphe_s4->GetYaxis()->FindBin(ph_cc))>0.65){

//nphe cut that removes the 1 and 18 segments
if (nphe > ph_el_arr_sim[pmt_hit+1][3][segment]){

//fiducial cut
if ((th_EL > th_min) && (th_EL < th_max)&& (ph_EL > fid_b+180) && (ph_EL < fid_a+180)){


	if ((P_EL < 1.75999) && (P_EL > 0.4)){
	if (indtype==1) ph_vs_th_el_sim[3][int((P_EL*100-40)/20)]->Fill(th_EL,ph_EL-180,1);
	};  

//th_vs_p 
//if ((th_EL > (11.7398+5./(0.433327*(P_EL+0.1)+0.015)) +8.)||(th_EL < (11.7398+12./(0.433327*(P_EL+0.1)+0.24)) -0.4)){   	
//if ((th_EL > (11.7398+5./(0.433327*(P_EL+0.1)+0.04)) +15.6)||(th_EL < (11.7398+5./(0.433327*(P_EL+0.1)+0.04)) +14.1)){	 


	if ((W>1.3)&&(W<1.825)&&(Q2>0.4)&&(Q2<1.)&&(n_PIp==1)&&(indtype==1)) th_vs_p_e_2_sim[3]->Fill(P_EL,th_EL,sigma);

//vertex cut
if ((z_EL>-2.45) && (z_EL< 2.05)){
	
//	if ((indtype==1)) hist_z_el_4_sim_1->Fill(z_EL,sigma);	
	
cuts_sim = true; 

   };//vertex cut
//   };//th_vs_p
//   };//th_vs_p
   }; //fiducial 
   };//nphe cut that removes the 1 and 18 segments
   };//geometrical cut on number of photoelectrons
   };//th_cc vs seg cut   
   };//phi_cc matching
   };//~fid cuts in Cherenkov plane
   };//~fid cuts in Cherenkov plane
   };//~fid cuts in Cherenkov plane
   };//~fid cuts in Cherenkov plane
      
     
 }; // ectot vs p cut

 break;  
 
case 5 : 
if (indtype ==2)  hist_z_el_5_sim_2->Fill(z_EL,1.); 
//if ((P_EL < 1.75999) && (P_EL > 0.4)){
//if (indtype==2) ph_vs_th_el_sim[4][int((P_EL*100-40)/20)]->Fill(th_EL,ph_EL-240,1);
//}; 

// ectot vs p cut
if ((ECT/P_EL > -0.0578053*P_EL*P_EL+0.182369*P_EL+0.00157396) && (ECT/P_EL < -0.000865715*P_EL*P_EL-0.0193174*P_EL+0.340688)){

	if (indtype==1) hist_ectot_sector5_sim->Fill(P_EL,ECT/P_EL,sigma);  
   
//~fid cuts in Cherenkov plane
if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {

//phi_cc matching
if (((ph_cc>0)&&(pmt_hit==1))||((ph_cc<0)&&(pmt_hit==-1))||(pmt_hit==0)){

//th_cc vs seg cut
if ((theta_cc >th_vs_seg_cc_arr[1][4][segment])&&(theta_cc <th_vs_seg_cc_arr[0][4][segment])){	

	if (indtype==1) th_cc_vs_seg_5_sim->Fill(segment+1,theta_cc,sigma);
	
 // geometrical cut on number of photoelectrons
if (norm_nphe_s5->GetBinContent(norm_nphe_s5->GetXaxis()->FindBin(theta_cc),norm_nphe_s5->GetYaxis()->FindBin(ph_cc))>0.8){ 

//nphe cut that removes the 1 and 18 segments
if (nphe > ph_el_arr_sim[pmt_hit+1][4][segment]){ 

//fiducial cut
if ((th_EL > th_min) && (th_EL < th_max)&& (ph_EL > fid_b+240) && (ph_EL < fid_a+240)&&(PdHit_EL!=17)){

	if ((P_EL < 1.75999) && (P_EL > 0.4)){
	if (indtype==1) ph_vs_th_el_sim[4][int((P_EL*100-40)/20)]->Fill(th_EL,ph_EL-240,1);
	};

	if ((W>1.3)&&(W<1.825)&&(Q2>0.4)&&(Q2<1.)&&(n_PIp==1)&&(indtype==1)) th_vs_p_e_2_sim[4]->Fill(P_EL,th_EL,sigma);
	
//th_vs_p
if ((th_EL > (11.7398+8.21504/(0.433327*(P_EL+0.1)+0.158076))+19.9)||(th_EL < (11.7398+8.21504/(0.433327*(P_EL+0.1)+0.158076))+17.6)){
if ((th_EL > (11.7398+4.3/(0.433327*(P_EL+0.1)+0.158076))+15.8)||(th_EL < (11.7398+4.5/(0.433327*(P_EL+0.1)+0.158076))+13.5)){





//vertex cut
if ((z_EL>-2.9) && (z_EL< 1.6)){

//	if ((indtype==1)) hist_z_el_5_sim_1->Fill(z_EL,sigma);

   cuts_sim = true; 
  
}; //vertex cut  
};//th_vs_p
};//th_vs_p
}; //fiducial     
 


   };//nphe cut that removes the 1 and 18 segments
   };//geometrical cut on number of photoelectrons
   };//th_cc vs seg cut   
   };//phi_cc matching
   };//~fid cuts in Cherenkov plane
   };//~fid cuts in Cherenkov plane
   };//~fid cuts in Cherenkov plane
   };//~fid cuts in Cherenkov plane

 }; // ectot vs p cut

 break;   
 
case 6 : 
if (indtype ==2) hist_z_el_6_sim_2->Fill(z_EL,1.);
//if ((P_EL < 1.75999) && (P_EL > 0.4)){
//if (indtype==2) ph_vs_th_el_sim[5][int((P_EL*100-40)/20)]->Fill(th_EL,ph_EL-300,1);
//}; 
 
// ectot vs p cut
 if ((ECT/P_EL > -0.0514107*P_EL*P_EL+0.17164*P_EL+0.00253672) && (ECT/P_EL < 0.0123293*P_EL*P_EL-0.053898*P_EL+0.361742)) {
 
	if (indtype==1) hist_ectot_sector6_sim->Fill(P_EL,ECT/P_EL,sigma);  

//~fid cuts in Cherenkov plane
if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {

//phi_cc matching
if (((ph_cc>0)&&(pmt_hit==1))||((ph_cc<0)&&(pmt_hit==-1))||(pmt_hit==0)){

//th_cc vs seg cut 
if ((theta_cc >th_vs_seg_cc_arr[1][5][segment])&&(theta_cc <th_vs_seg_cc_arr[0][5][segment])){ 	

 	if (indtype==1) th_cc_vs_seg_6_sim->Fill(segment+1,theta_cc,sigma);

 // geometrical cut on number of photoelectrons
if (norm_nphe_s6->GetBinContent(norm_nphe_s6->GetXaxis()->FindBin(theta_cc),norm_nphe_s6->GetYaxis()->FindBin(ph_cc))>0.8){ 

//nphe cut that removes the 1 and 18 segments
if (nphe > ph_el_arr_sim[pmt_hit+1][5][segment]){

//fiducial cut
if ((th_EL > th_min) && (th_EL < th_max)&& (ph_EL > fid_b+300) && (ph_EL < fid_a+300)){

	if ((P_EL < 1.75999) && (P_EL > 0.4)){
	if (indtype==1) ph_vs_th_el_sim[5][int((P_EL*100-40)/20)]->Fill(th_EL,ph_EL-300,1);
	};


//th_vs_p 
//if ((th_EL > (11.7398+5./(0.433327*(P_EL+0.1)+0.015)) +8.3)||(th_EL < (11.7398+12./(0.433327*(P_EL+0.1)+0.27)) +0.3)){


	if ((W>1.3)&&(W<1.825)&&(Q2>0.4)&&(Q2<1.)&&(n_PIp==1)&&(indtype==1)) th_vs_p_e_2_sim[5]->Fill(P_EL,th_EL,sigma);


//vertex cut
if ((z_EL>-3.) && (z_EL<1.5)){

//	if ((indtype==1)) hist_z_el_6_sim_1->Fill(z_EL,sigma);

   cuts_sim = true; 

   }; //vertex cut
//   };//th_vs_p   
   }; //fiducial 
   };//nphe cut that removes the 1 and 18 segments
   };//geometrical cut on number of photoelectrons  
   };//th_cc vs seg cut   
   };//phi_cc matching
   };//~fid cuts in Cherenkov plane
   };//~fid cuts in Cherenkov plane
   };//~fid cuts in Cherenkov plane
   };//~fid cuts in Cherenkov plane

 }; // ectot vs p cut
 break;      
 
   
   }; // end of switch
   };  // end of calorimeter threshold cut
  
   
   return cuts_sim;
   
   };
   
   
   //////////////////////
   
      bool cuts_sim::Proton_cuts_sim(){
       
   bool cuts_sim;
   Float_t m_p,p_fid_a,p_fid_b, th_min, th_max;
   m_p=0.938272; 
   
   th_min = 12.;
   th_max = 60.;
  // th_max2 = 94.0079-98.7448*P_P + 84.5143*P_P*P_P-32.0668*P_P*P_P*P_P;   
   
  p_fid_a = 25.*(1-exp(-1.*0.12*(th_P-10.)))-3.;
  p_fid_b = -25.*(1-exp(-1.*0.12*(th_P-10.)))+3.; 
  
   cuts_sim = false; 
   
//Proton id beta_vs_p cut
if((n_P == 1)&& (beta_P < 0.958174*P_P/sqrt(m_p*m_p+0.921203*P_P*P_P-0.0289732)+ 0.0383292) && (beta_P > 0.97261*P_P/sqrt(m_p*m_p +0.89435*P_P*P_P -0.199588) -0.0775468)&&(PdHit_P !=48)){

//----------PROTON ENERGY LOSS----------------------------
P_P = P_P + delta_mom_p_skor;
//--------------------------------------------------------

//minimum proton momentum & mimimum proton beta cuts
if ((P_P > 0.235)&&(beta_P > 0.19)){

//vertex cut
if ((z_P > -4.8) && (z_P < 4.)){

if ((ph_P >= 330)&& (ph_P <= 360)) {

//fiducial sector1
if ((th_P > th_min)&&(th_P < th_max)&&(ph_P > p_fid_b+360) && (ph_P < p_fid_a+360)){

	if ((P_P > 0.2) && (P_P < 1.2)&&(indtype==1)) ph_vs_th_p_sim[0]->Fill(th_P,ph_P-360,1);

	if ((W>1.3)&&(n_PIp==1)&&(indtype==1)) th_vs_p_p_2_sim[0]->Fill(P_P,th_P,sigma);	
//th_vs_p1 
if (!((th_P<pow((P_P-0.03-0.304992),(0.0758186))*91.5643-48.2057 + 2. )&&!(th_P<pow((P_P-0.08-0.304992),(0.0758186))*91.5643-48.2057 +1.5))){
//th_vs_p2 
if (!((th_P <(304.23*(P_P+0.2)*(P_P+0.2)*(P_P+0.2) -255.798*(P_P+0.2)*(P_P+0.2)+497.462*(P_P+0.2) +38.0385)*exp(-1.6*(P_P+0.2))-104.)&&(th_P>(304.23*(P_P-0.12)*(P_P-0.12)*(P_P-0.12) -255.798*(P_P-0.12)*(P_P-0.12)+497.462*(P_P-0.12) +38.0385)*exp(-1.6*(P_P-0.12))-103.))){
	
//	if ((indtype==1)) hist_z_el_1_sim_1->Fill(z_P,sigma);

cuts_sim = true;

};//th_vs_p2
};//th_vs_p1
};//end of fiducial cut for the first part of the first sector

 
};//end of the first part of the first sector
 
if ((ph_P >= 0) && (ph_P <= 30)) {

//fiducial sector1 
if ((th_P > th_min)&&(th_P < th_max)&&(ph_P > p_fid_b) && (ph_P < p_fid_a)){

	if ((P_P > 0.2) && (P_P < 1.2)&&(indtype==1)) ph_vs_th_p_sim[0]->Fill(th_P,ph_P,1);
	
	if ((W>1.3)&&(n_PIp==1)&&(indtype==1)) th_vs_p_p_2_sim[0]->Fill(P_P,th_P,sigma);	
//th_vs_p1 
if (!((th_P<pow((P_P-0.03-0.304992),(0.0758186))*91.5643-48.2057 + 2. )&&!(th_P<pow((P_P-0.08-0.304992),(0.0758186))*91.5643-48.2057 +1.5))){
//th_vs_p2 
if (!((th_P <(304.23*(P_P+0.2)*(P_P+0.2)*(P_P+0.2) -255.798*(P_P+0.2)*(P_P+0.2)+497.462*(P_P+0.2) +38.0385)*exp(-1.6*(P_P+0.2))-104.)&&(th_P>(304.23*(P_P-0.12)*(P_P-0.12)*(P_P-0.12) -255.798*(P_P-0.12)*(P_P-0.12)+497.462*(P_P-0.12) +38.0385)*exp(-1.6*(P_P-0.12))-103.))){

//	if ((indtype==1)) hist_z_el_1_sim_1->Fill(z_P,sigma);

cuts_sim = true;

};//th_vs_p2
};//th_vs_p1
};//end of fiducial cut for the second part of the first sector

};//end of the second part of the first sector
 
 
if ((ph_P >= 30) && (ph_P <=90)&&(PdHit_P !=16)){ 

//fiducial sector2
if ((th_P > th_min)&&(th_P < th_max)&&(ph_P > p_fid_b+60) && (ph_P < p_fid_a+60)){

	if ((P_P > 0.2) && (P_P < 1.2)&&(indtype==1)) ph_vs_th_p_sim[1]->Fill(th_P,ph_P-60,1); 

	if ((W>1.3)&&(n_PIp==1)&&(indtype==1)) th_vs_p_p_2_sim[1]->Fill(P_P,th_P,sigma); 
	
//th_vs_p1
if (!((th_P<pow((P_P-0.03-0.304992),(0.0758186))*91.5643-48.2057 + 2. )&&!(th_P<pow((P_P-0.08-0.304992),(0.0758186))*91.5643-48.2057 +1.5))){
//th_vs_p2
if (!((th_P <pow((P_P-0.12-0.304992),(0.0758186))*91.5643-48.2057 +1.1)&&!(th_P<pow((P_P-0.16-0.304992),(0.0758186))*91.5643-48.2057))){
//th_vs_p3
if (!((th_P <26.5087*(P_P+0.04)*(P_P+0.04)*(P_P+0.04) -116.557*(P_P+0.04)*(P_P+0.04)+ 175.167*(P_P+0.04)-61.7717+1.5)&&(th_P>26.5087*(P_P-0.03)*(P_P-0.03)*(P_P-0.03) -116.557*(P_P-0.03)*(P_P-0.03)+ 175.167*(P_P-0.03)-61.7717-1.8))){
 
//	if ((indtype==1)) hist_z_el_2_sim_1->Fill(z_P,sigma);

cuts_sim = true;

};//th_vs_p3
};//th_vs_p2
};//th_vs_p1
};//end of the fiducial cut for sector2
};//end of the sector2

 
if ((ph_P >=90) && (ph_P <=150)&&(PdHit_P !=44)) {

//fiducial sector3
if ((th_P > th_min)&&(th_P < th_max)&&(ph_P > p_fid_b+120) && (ph_P < p_fid_a+120)){

	if ((P_P > 0.2) && (P_P < 1.2)&&(indtype==1)) ph_vs_th_p_sim[2]->Fill(th_P,ph_P-120,1);

	if ((W>1.3)&&(n_PIp==1)&&(indtype==1)) th_vs_p_p_2_sim[2]->Fill(P_P,th_P,sigma);
		
//th_vs_p1	
if (!((th_P<pow((P_P-0.03-0.304992),(0.0758186))*91.5643-48.2057 + 2. )&&!(th_P<pow((P_P-0.08-0.304992),(0.0758186))*91.5643-48.2057 +1.5))){
//th_vs_p2
if (!((th_P <pow((P_P-0.12-0.304992),(0.0758186))*91.5643-48.2057 +1.1)&&!(th_P<pow((P_P-0.16-0.304992),(0.0758186))*91.5643-48.2057))){	

//	if ((indtype==1)) hist_z_el_3_sim_1->Fill(z_P,sigma);

cuts_sim = true;

};//th_vs_p2
};//th_vs_p1
};//end of the fiducial cut for sector3
};//end of the sector3

 
if ((ph_P >= 150) && (ph_P <= 210)) {

//fiducial sector4
if ((th_P > th_min)&&(th_P < th_max)&&(ph_P > p_fid_b+180) && (ph_P < p_fid_a+180)){

	if ((P_P > 0.2) && (P_P < 1.2)&&(indtype==1)) ph_vs_th_p_sim[3]->Fill(th_P,ph_P-180,1);

	if ((W>1.3)&&(n_PIp==1)&&(indtype==1)) th_vs_p_p_2_sim[3]->Fill(P_P,th_P,sigma);
	
//th_vs_p1	
if (!((th_P<pow((P_P-0.03-0.304992),(0.0758186))*91.5643-48.2057 + 2. )&&!(th_P<pow((P_P-0.08-0.304992),(0.0758186))*91.5643-48.2057 +1.5))){
//th_vs_p2
if (!((th_P <pow((P_P-0.12-0.304992),(0.0758186))*91.5643-48.2057 +1.1)&&!(th_P<pow((P_P-0.16-0.304992),(0.0758186))*91.5643-48.2057))){	

//	if ((indtype==1)) hist_z_el_4_sim_1->Fill(z_P,sigma);

cuts_sim = true;

};//th_vs_p2
};//th_vs_p1
};//end of the fiducial cut for sector4
};//end of the sector4


if ((ph_P >= 210) && (ph_P <=270)&&(PdHit_P !=17)) {

//fiducial sector5
if ((th_P > th_min)&&(th_P < th_max)&&(ph_P > p_fid_b+240) && (ph_P < p_fid_a+240)){

	if ((P_P > 0.2) && (P_P < 1.2)&&(indtype==1)) ph_vs_th_p_sim[4]->Fill(th_P,ph_P-240,1);

	if ((W>1.3)&&(n_PIp==1)&&(indtype==1)) th_vs_p_p_2_sim[4]->Fill(P_P,th_P,sigma);
	
//th_vs_p1
if (!((th_P<pow((P_P-0.03-0.304992),(0.0758186))*91.5643-48.2057 + 2. )&&!(th_P<pow((P_P-0.08-0.304992),(0.0758186))*91.5643-48.2057 +1.5))){
//th_vs_p2
if (!((th_P <pow((P_P -0.321436),(0.0704348))*88.0419-46.9342-0.5)&&!(th_P<pow((P_P -0.12 -0.321436),(0.0704348))*88.0419-46.9342-2.5))){
//th_vs_p3
if (!((th_P <31.2482*(P_P+0.045)*(P_P+0.045)*(P_P+0.045) -135.817*(P_P+0.045)*(P_P+0.045)+ 198.038*(P_P+0.045)-66.968+1.5)&&(th_P>31.2482*(P_P-0.01)*(P_P-0.01)*(P_P-0.01) -135.817*(P_P-0.01)*(P_P-0.01)+ 198.038*(P_P-0.01)-66.968-2.5))){

//	if ((indtype==1)) hist_z_el_5_sim_1->Fill(z_P,sigma);

cuts_sim = true;

};//th_vs_p3
};//th_vs_p2
};//th_vs_p1
};//end of the fiducial cut for sector5
};//end of the sector5

 
if ((ph_P >= 270) && (ph_P <=330)) {

//fiducial sector6
if ((th_P > th_min)&&(th_P < th_max)&&(ph_P > p_fid_b+300) && (ph_P < p_fid_a+300)){

	if ((P_P > 0.2) && (P_P < 1.2)&&(indtype==1)) ph_vs_th_p_sim[5]->Fill(th_P,ph_P-300,1);

	if ((W>1.3)&&(n_PIp==1)&&(indtype==1)) th_vs_p_p_2_sim[5]->Fill(P_P,th_P,sigma);

//th_vs_p1
if (!((th_P<pow((P_P-0.03-0.304992),(0.0758186))*91.5643-48.2057 + 2. )&&!(th_P<pow((P_P-0.08-0.304992),(0.0758186))*91.5643-48.2057 +1.5))){
//th_vs_p2
if (!((th_P <pow((P_P-0.12-0.304992),(0.0758186))*91.5643-48.2057 +1.1)&&!(th_P<pow((P_P-0.16-0.304992),(0.0758186))*91.5643-48.2057))){

//	if ((indtype==1)) hist_z_el_6_sim_1->Fill(z_P,sigma);

cuts_sim = true;

};//th_vs_p2
};//th_vs_p1
};//end of the fiducial cut for sector6
};//end of the sector6
 
};//vertex cut 
};// end of minimum proton momentum & mimimum proton beta cuts 
};//end proton id beta vs p cut
 

return cuts_sim;
};  
    
    
    ///////////////////////////////
    
    bool cuts_sim::PIp_cuts_sim(){
       
   bool cuts_sim;
   Float_t m_pip,pip_fid_a,pip_fid_b,th_min,th_max;
   m_pip = 0.13957;

   th_min = 12.;
   th_max = 120.;
//    th_max2 = 222.947-505.444*P_PIp + 534.834*P_PIp*P_PIp-223.395*P_PIp*P_PIp*P_PIp;
    
  pip_fid_a = 25.*(1-exp(-1.*0.12*(th_PIp-10.)))-3.;
  pip_fid_b = -25.*(1-exp(-1.*0.12*(th_PIp-10.)))+3.;
  
   cuts_sim = false; 
 
     
if((n_PIp == 1)&&(beta_PIp < 0.838148*P_PIp/sqrt(m_pip*m_pip+0.793572*P_PIp*P_PIp -0.00291347)+ 0.105962) && (beta_PIp > 0.819388*P_PIp/sqrt(m_pip*m_pip +0.600785*P_PIp*P_PIp -0.00885673) -0.107038)&&(PdHit_PIp !=48)){ 
 
//vertex cut
if ((z_PIp>-4.8) && (z_PIp<4.)){
     
if ((ph_PIp >= 330) && (ph_PIp <=360)){


//fiducial sector1  
if ((th_PIp > th_min)&&(th_PIp < th_max)&&(ph_PIp > pip_fid_b+360) && (ph_PIp < pip_fid_a+360)){

	if ((P_PIp > 0) && (P_PIp < 1.2)&&(indtype==1)) ph_vs_th_pip_sim[0]->Fill(th_PIp,ph_PIp-360,1); 
	if (indtype==1) th_vs_p_pip_2_sim[0]->Fill(P_PIp,th_PIp,sigma);

//th_vs_p1
if (!((th_PIp<(304.23*(P_PIp+0.15)*(P_PIp+0.15)*(P_PIp+0.15) -255.798*(P_PIp+0.15)*(P_PIp+0.15)+497.462*(P_PIp+0.15) +38.0385)*exp(-1.85*(P_PIp+0.15)) +5.5 )&&!(th_PIp<(304.23*(P_PIp+0.15)*(P_PIp+0.15)*(P_PIp+0.15) -255.798*(P_PIp+0.15)*(P_PIp+0.15)+497.462*(P_PIp+0.15) +38.0385)*exp(-1.85*(P_PIp+0.15)) -3))){
//th_vs_p2
if (!((th_PIp<(pow((P_PIp-0.0575818),( 0.075643))*238.248-115.039)*exp(-0.5*P_PIp)-0.1)&&!(th_PIp<(pow((P_PIp-0.103718),(0.0703664))*252.822-133.024)*exp(-0.5*P_PIp)+0.1))){
//th_vs_p3
if (!((th_PIp<(304.23*(P_PIp+0.25)*(P_PIp+0.25)*(P_PIp+0.25) -255.798*(P_PIp+0.25)*(P_PIp+0.25)+497.462*(P_PIp+0.25) +38.0385)*exp(-1.6*(P_PIp+0.25))-104)&&!(th_PIp<(304.23*(P_PIp-0.12)*(P_PIp-0.12)*(P_PIp-0.12) -255.798*(P_PIp-0.12)*(P_PIp-0.12)+497.462*(P_PIp-0.12) +38.0385)*exp(-1.6*(P_PIp-0.12))-103))){

//	if ((indtype==1)) hist_z_el_1_sim_1->Fill(z_PIp,sigma);

cuts_sim = true;
   
};//th_vs_p3
};//th_vs_p2
};//th_vs_p1
};//end of the fiducial cut for first part of sector1
}; //end of the first part of sector1
  
if ((ph_PIp >= 0) && (ph_PIp <=30)){
  
//fiducial sector1  
if ((th_PIp > th_min)&&(th_PIp < th_max)&&(ph_PIp > pip_fid_b) && (ph_PIp < pip_fid_a)){

	if ((P_PIp > 0) && (P_PIp < 1.2)&&(indtype==1)) ph_vs_th_pip_sim[0]->Fill(th_PIp,ph_PIp,1);
	if (indtype==1) th_vs_p_pip_2_sim[0]->Fill(P_PIp,th_PIp,sigma);
	
//th_vs_p1
if (!((th_PIp<(304.23*(P_PIp+0.15)*(P_PIp+0.15)*(P_PIp+0.15) -255.798*(P_PIp+0.15)*(P_PIp+0.15)+497.462*(P_PIp+0.15) +38.0385)*exp(-1.85*(P_PIp+0.15)) +5.5 )&&!(th_PIp<(304.23*(P_PIp+0.15)*(P_PIp+0.15)*(P_PIp+0.15) -255.798*(P_PIp+0.15)*(P_PIp+0.15)+497.462*(P_PIp+0.15) +38.0385)*exp(-1.85*(P_PIp+0.15)) -3))){
//th_vs_p2
if (!((th_PIp<(pow((P_PIp-0.0575818),( 0.075643))*238.248-115.039)*exp(-0.5*P_PIp)-0.1)&&!(th_PIp<(pow((P_PIp-0.103718),(0.0703664))*252.822-133.024)*exp(-0.5*P_PIp)+0.1))){
//th_vs_p3
if (!((th_PIp<(304.23*(P_PIp+0.25)*(P_PIp+0.25)*(P_PIp+0.25) -255.798*(P_PIp+0.25)*(P_PIp+0.25)+497.462*(P_PIp+0.25) +38.0385)*exp(-1.6*(P_PIp+0.25))-104)&&!(th_PIp<(304.23*(P_PIp-0.12)*(P_PIp-0.12)*(P_PIp-0.12) -255.798*(P_PIp-0.12)*(P_PIp-0.12)+497.462*(P_PIp-0.12) +38.0385)*exp(-1.6*(P_PIp-0.12))-103))){

// 	if ((indtype==1)) hist_z_el_1_sim_1->Fill(z_PIp,sigma);
  
cuts_sim = true;
};//th_vs_p3
};//th_vs_p2
};//th_vs_p1
};//end of the fiducial cut for second part of sector1
}; //end of the second part of sector1
  
  
if ((ph_PIp >= 30) && (ph_PIp <=90)&&(PdHit_PIp !=16)) {
  
//fiducial sector2
if ((th_PIp > th_min)&&(th_PIp < th_max)&&(ph_PIp > pip_fid_b+60) && (ph_PIp < pip_fid_a+60)){

	if ((P_PIp > 0) && (P_PIp < 1.2)&&(indtype==1)) ph_vs_th_pip_sim[1]->Fill(th_PIp,ph_PIp-60,1);
	if (indtype==1) th_vs_p_pip_2_sim[1]->Fill(P_PIp,th_PIp,sigma);  
	
//th_vs_p1
if (!((th_PIp<pow((P_PIp-0.415068),(0.226449))*48.7564 + 2.79478-1.)&&!(th_PIp<pow((P_PIp-0.449975),( 0.315164 ))*36.608 +  9.74262-1.))){
//th_vs_p2
if (!((th_PIp<(387.289*P_PIp*P_PIp*P_PIp -758.466*P_PIp*P_PIp+ 842.881*P_PIp-299.953+15.)*exp(-2*P_PIp))&&!(th_PIp<(387.289*(P_PIp+0.03)*(P_PIp+0.03)*(P_PIp+0.03) -758.466*(P_PIp+0.03)*(P_PIp+0.03)+ 842.881*(P_PIp+0.03)-299.953-15.)*exp(-2*(P_PIp+0.03))-1.5))){

// 	if ((indtype==1)) hist_z_el_2_sim_1->Fill(z_PIp,sigma);

cuts_sim = true;
};//th_vs_p2
};//th_vs_p1
};//end of the fiducial cut for sector2
};//end of the sector2
  
  
if ((ph_PIp >= 90) && (ph_PIp <=150)&&(PdHit_PIp !=44)) {

//fiducial sector3
if ((th_PIp > th_min)&&(th_PIp < th_max)&&(ph_PIp > pip_fid_b+120) && (ph_PIp < pip_fid_a+120)){

	if ((P_PIp > 0) && (P_PIp < 1.2)&&(indtype==1)) ph_vs_th_pip_sim[2]->Fill(th_PIp,ph_PIp-120,1);
 	if (indtype==1) th_vs_p_pip_2_sim[2]->Fill(P_PIp,th_PIp,sigma);
	
//th_vs_p1
if (!((th_PIp<(10000*P_PIp*P_PIp*P_PIp-3607.41*P_PIp*P_PIp+ 1725.72*P_PIp-10.6776)*exp(-4.7*P_PIp))&&!(th_PIp<(10000*P_PIp*P_PIp*P_PIp-4505.62*P_PIp*P_PIp+  2056.24  *P_PIp -77.4077 +5.)*exp(-4.7*P_PIp)))){
//th_vs_p2
if (!((th_PIp<pow((P_PIp-0.416536),(0.108376))*67.4593-21.4374)&&!(th_PIp<pow((P_PIp -0.454898),( 0.289291))* 35.7267+6.65908+1.5))){

// 	if ((indtype==1)) hist_z_el_3_sim_1->Fill(z_PIp,sigma);

cuts_sim = true;
};//th_vs_p2
};//th_vs_p1
};//end of the fiducial cut for sector3
};//end of the sector3
  
if ((ph_PIp >= 150) && (ph_PIp <=210)){
 
//fiducial sector4
if ((th_PIp > th_min)&&(th_PIp < th_max)&&(ph_PIp > pip_fid_b+180) && (ph_PIp < pip_fid_a+180)){
 
	if ((P_PIp > 0) && (P_PIp < 1.2)&&(indtype==1)) ph_vs_th_pip_sim[3]->Fill(th_PIp,ph_PIp-180,1); 
	if (indtype==1) th_vs_p_pip_2_sim[3]->Fill(P_PIp,th_PIp,sigma);  	
	 
//th_vs_p1	
if (!((th_PIp<(304.23*(P_PIp+0.165)*(P_PIp+0.165)*(P_PIp+0.165) -255.798*(P_PIp+0.165)*(P_PIp+0.165)+497.462*(P_PIp+0.165) +38.0385)*exp(-1.85*(P_PIp+0.165)) +5.)&&!(th_PIp<(304.23*(P_PIp+0.18)*(P_PIp+0.18)*(P_PIp+0.18) -255.798*(P_PIp+0.18)*(P_PIp+0.18)+497.462*(P_PIp+0.18) +38.0385)*exp(-1.85*(P_PIp+0.18)) - 1.))){
//th_vs_p2
if (!((th_PIp<(1600*(P_PIp+0.03)*(P_PIp+0.03)*(P_PIp+0.03) -1068.36*(P_PIp+0.03)*(P_PIp+0.03)+ 775.016*(P_PIp+0.03)-1.13034)*exp(-2.75*(P_PIp+0.03)))&&!(th_PIp<(pow((P_PIp-0.103718),(0.0703664))*252.822-133.024)*exp(-0.45*P_PIp)-7.))){
//th_vs_p3
if (!((th_PIp<pow((P_PIp-0.416536),(0.108376))*67.4593-21.4374)&&!(th_PIp<pow((P_PIp -0.454898),( 0.289291))* 35.7267+6.65908+1.5))){

// 	if ((indtype==1)) hist_z_el_4_sim_1->Fill(z_PIp,sigma);

cuts_sim = true;
};//th_vs_p3
};//th_vs_p2
};//th_vs_p1
};//end of the fiducial cut for sector4
}; //end of the sector4
  
if ((ph_PIp >= 210) && (ph_PIp <=270)&&(PdHit_PIp !=17)) {

//fiducial sector5
if ((th_PIp > th_min)&&(th_PIp < th_max)&&(ph_PIp > pip_fid_b+240) && (ph_PIp < pip_fid_a+240)){

	if ((P_PIp > 0) && (P_PIp < 1.2)&&(indtype==1)) ph_vs_th_pip_sim[4]->Fill(th_PIp,ph_PIp-240,1);
	if (indtype==1) th_vs_p_pip_2_sim[4]->Fill(P_PIp,th_PIp,sigma);  
		
//th_vs_p1
if (!((th_PIp<(525.498*(P_PIp+0.03)*(P_PIp+0.03)*(P_PIp+0.03) -1284.98*(P_PIp+0.03)*(P_PIp+0.03)+1460.67*(P_PIp+0.03)-499.999)*exp(-1.94*(P_PIp+0.03)))&&!(th_PIp<(525.498*(P_PIp-0.02)*(P_PIp-0.02)*(P_PIp-0.02) -1284.98*(P_PIp-0.02)*(P_PIp-0.02)+1460.67*(P_PIp-0.02)-499.999)*exp(-1.94*(P_PIp-0.02)) -4.7))){
//th_vs_p2
if (!((th_PIp<(304.23*(P_PIp)*(P_PIp)*(P_PIp) -255.798*(P_PIp)*(P_PIp)+497.462*(P_PIp) +38.0385)*exp(-1.85*(P_PIp)))&&!(th_PIp<(304.23*(P_PIp+0.03)*(P_PIp+0.03)*(P_PIp+0.03) -255.798*(P_PIp+0.03)*(P_PIp+0.03)+497.462*(P_PIp+0.03) +38.0385)*exp(-1.85*(P_PIp+0.03)) -11.))){
//th_vs_p3
if (!((th_PIp<pow((P_PIp-0.304992),(0.0758186))*91.5643-48.2057 - 1.)&&!(th_PIp<pow((P_PIp -0.36848),( 0.0864219))*70.4769  -34.9998+ 1.5))){

// 	if ((indtype==1)) hist_z_el_5_sim_1->Fill(z_PIp,sigma);
	
cuts_sim = true;
};//th_vs_p3
};//th_vs_p2
};//th_vs_p1
};//end of the fiducial cut for sector5
};//end of the sector5
  
  
if ((ph_PIp >= 270) && (ph_PIp <=330)){
 
//fiducial sector6
if ((th_PIp > th_min)&&(th_PIp < th_max)&&(ph_PIp > pip_fid_b+300) && (ph_PIp < pip_fid_a+300)){
 
	if ((P_PIp > 0) && (P_PIp < 1.2)&&(indtype==1)) ph_vs_th_pip_sim[5]->Fill(th_PIp,ph_PIp-300,1); 
	if (indtype==1) th_vs_p_pip_2_sim[5]->Fill(P_PIp,th_PIp,sigma); 
		 
//th_vs_p1
if (!((th_PIp<pow((P_PIp-0.05-0.0942469),( 0.0582707))*114.358-50 -0.5)&&!(th_PIp<pow((P_PIp-0.05-0.126994),( 0.0706829))* 110.073-50+2.))){
//th_vs_p2
if (!((th_PIp<pow((P_PIp-0.416536),(0.108376))*67.4593-21.4374-1.)&&!(th_PIp<pow((P_PIp -0.454098),(0.0912936))*58.2946-20.4843+1.5))){

// 	if ((indtype==1)) hist_z_el_6_sim_1->Fill(z_PIp,sigma);

cuts_sim = true;
};//th_vs_p2
};//th_vs_p1
};//end of the fiducial cut for sector6
};//end of the sector6
  
  
  
 };//vertex cut  
  };//end of pi+ id beta vs p cut
    return cuts_sim;
   };
   
   
   
   //////////////////////////////////////////
   
   
   
  bool cuts_sim::PIm_cuts_sim(){
       
   bool cuts_sim;
   Float_t m_pim,th_min,th_max,th_edge,par1,par2,par3,pim_fid_a,pim_fid_b,th_min_lowp, a_shrink;
    
   m_pim = 0.13957; 

    th_min=(11.09+8./(0.472*(P_PIm-0.03)+0.117))-1.;
    th_min_lowp = 30.+5.24894e-05/(5.71075e-05*(P_PIm+0.004)*(P_PIm+0.004)+4.44089e-16)+3.;


   th_max = 140.;
// if (P_PIm<0.23) th_max = 140.;  
  //if (P_PIm>0.23) th_max = 238.213-537.237*P_PIm + 546.985*P_PIm*P_PIm-212.017*P_PIm*P_PIm*P_PIm;
  
  th_edge = th_max- 5.2/3.*(th_max-th_min)/2.;
  
  par1=0.61+1.18*P_PIm;
  par2=-59.2-35.3*P_PIm; 

  a_shrink = 1.;
  
/*  if ((Q2>0.45)&&(Q2<0.50)&&(W>1.3)&&(W<1.325)) a_shrink = 2.5;
  if ((Q2>0.50)&&(Q2<0.55)&&(((W>1.3)&&(W<1.325))||((W>1.35)&&(W<1.375)))) a_shrink = 2.5;
  if ((Q2>0.55)&&(Q2<0.60)&&(W>1.325)&&(W<1.4)) a_shrink = 2.5;
  if ((Q2>0.60)&&(Q2<0.65)&&(W>1.35)&&(W<1.375)) a_shrink = 2.5;
  if ((Q2>0.65)&&(Q2<0.70)&&(((W>1.3)&&(W<1.325))||((W>1.35)&&(W<1.4)))) a_shrink = 2.5;
  if ((Q2>0.70)&&(Q2<0.75)&&(W>1.35)&&(W<1.425)) a_shrink = 2.5;
  if ((Q2>0.75)&&(Q2<0.80)&&(W>1.35)&&(W<1.4)) a_shrink = 2.5;
  if ((Q2>0.80)&&(Q2<0.85)&&(((W>1.3)&&(W<1.325))||((W>1.35)&&(W<1.4)))) a_shrink = 2.5;
  if ((Q2>0.85)&&(Q2<0.90)&&(W>1.35)&&(W<1.4)) a_shrink = 2.5;
  if ((Q2>0.90)&&(Q2<0.95)&&(W>1.35)&&(W<1.375)) a_shrink = 2.5;
  if ((Q2>0.95)&&(Q2<1.00)&&(((W>1.325)&&(W<1.375))||((W>1.4)&&(W<1.425)))) a_shrink = 2.5;
  */
  
  par3=17.2*P_PIm-11.9*P_PIm*P_PIm-a_shrink;  
      
   pim_fid_a=+23.5*pow((sin((th_PIm-th_min)*0.015)),(par1+par2/th_PIm+1400./th_PIm/th_PIm))+par3;
   pim_fid_b=-23.5*pow((sin((th_PIm-th_min)*0.015)),(par1+par2/th_PIm+1400./th_PIm/th_PIm))-par3;
    
  if (th_PIm > th_edge) pim_fid_a=+23.5*pow((sin((th_edge-th_min)*0.015)),(par1+par2/th_edge+1400./th_edge/th_edge))+par3;
  if (th_PIm > th_edge) pim_fid_b=-23.5*pow((sin((th_edge-th_min)*0.015)),(par1+par2/th_edge+1400./th_edge/th_edge))-par3;   

	   
cuts_sim = false; 
  
//pi- id beta vs p cut        
if((n_PIm == 1)&& (beta_PIm < 0.779984*P_PIm/sqrt(m_pip*m_pip+0.700052*P_PIm*P_PIm-0.00506451 )+ 0.105877) && (beta_PIm > 0.854914*P_PIm/sqrt(m_pip*m_pip +0.629965*P_PIm*P_PIm -0.00757789)-0.113261 )&&(PdHit_PIm!=48)){    

//vertex cut
if ((z_PIm>-4.8) && (z_PIm<4.)){

//th_min for low momenta
if (((P_PIm<0.3)&&(th_PIm > th_min_lowp))||(P_PIm>0.3)){ 
 
if ((ph_PIm >= 330) && (ph_PIm <=360)){
 
//fiducial cut sector 1  
if ((th_PIm > th_min)&&(th_PIm < th_max)&&(ph_PIm > pim_fid_b+360) && (ph_PIm < pim_fid_a+360)){

	if ((P_PIm > 0) && (P_PIm < 2.)&&(indtype==1)) ph_vs_th_pim_sim[0][int(P_PIm/0.4)]->Fill(th_PIm,ph_PIm-360,1);

	if (indtype==1) th_vs_p_pim_2_sim[0]->Fill(P_PIm,th_PIm,sigma);
	
//th_vs_p
if (th_PIm < (11.09+8./(0.472*(P_PIm+0.3)+0.117))+77.){	 

 	if ((indtype==1)) hist_z_el_1_sim_1->Fill(z_PIm,sigma);	

 cuts_sim = true;
};//th_vs_p
};//end of the fiducial cut for first part of sector1
}; //end of the first part of sector1
  
if ((ph_PIm >= 0) && (ph_PIm <=30)){

//fiducial cut sector 1  
if ((th_PIm > th_min)&&(th_PIm < th_max)&&(ph_PIm > pim_fid_b) && (ph_PIm < pim_fid_a)){

	if ((P_PIm > 0) && (P_PIm < 2.)&&(indtype==1)) ph_vs_th_pim_sim[0][int(P_PIm/0.4)]->Fill(th_PIm,ph_PIm,1);

	if (indtype==1) th_vs_p_pim_2_sim[0]->Fill(P_PIm,th_PIm,sigma);
	
//th_vs_p
if (th_PIm < (11.09+8./(0.472*(P_PIm+0.3)+0.117))+77.){	 


 	if ((indtype==1)) hist_z_el_1_sim_1->Fill(z_PIm,sigma);
		 
cuts_sim = true;
};//th_vs_p
};//end of the fiducial cut for second part of sector1
};//end of the second part of sector1
 
if ((ph_PIm >= 30) && (ph_PIm <=90)&&(PdHit_PIm!=16)) {
    
//fiducial cut sector 2 
if ((th_PIm > th_min)&&(th_PIm < th_max)&&(ph_PIm > pim_fid_b+60) && (ph_PIm < pim_fid_a+60)){

	if ((P_PIm > 0) && (P_PIm < 2.)&&(indtype==1)) ph_vs_th_pim_sim[1][int(P_PIm/0.4)]->Fill(th_PIm,ph_PIm-60,1);

	if (indtype==1) th_vs_p_pim_2_sim[1]->Fill(P_PIm,th_PIm,sigma);	
//th_vs_p
if ((th_PIm < 36.152+3.69909e-05/(5.40783e-06*P_PIm+1.81169e-07)-2.5) || (th_PIm > 36.152+3.69909e-05/(5.40783e-06*P_PIm+1.81169e-07)+2.)){



 	if ((indtype==1)) hist_z_el_2_sim_1->Fill(z_PIm,sigma);

cuts_sim = true;
};//th_vs_p
};//end of the fiducial cut for sector2
};//end of the sector2
    
if ((ph_PIm >= 90) && (ph_PIm <=150)&&(PdHit_PIm!=44)) {

//fiducial cut sector 3
if ((th_PIm > th_min)&&(th_PIm < th_max)&&(ph_PIm > pim_fid_b+120) && (ph_PIm < pim_fid_a+120)){

	if ((P_PIm > 0) && (P_PIm < 2.)&&(indtype==1)) ph_vs_th_pim_sim[2][int(P_PIm/0.4)]->Fill(th_PIm,ph_PIm-120,1); 

	if (indtype==1) th_vs_p_pim_2_sim[2]->Fill(P_PIm,th_PIm,sigma);	
	
//th_vs_p
if (th_PIm <(11.09+8./(0.472*(P_PIm+0.2)+0.117))+85.){ 
if ((th_PIm <36.152+3.69909e-05/(5.40783e-06*P_PIm+1.5e-07)-8.5)||(th_PIm >36.152+3.69909e-05/(5.40783e-06*P_PIm+0.2e-07)-5.5)){ 
 
 
 	if ((indtype==1)) hist_z_el_3_sim_1->Fill(z_PIm,sigma);
	
cuts_sim = true;
};//th_vs_p
};//th_vs_p
};//end of the fiducial cut for sector3
};  //end of the sector3
  
if ((ph_PIm >= 150) && (ph_PIm <=210)){
 
//fiducial cut sector 4 
if ((th_PIm > th_min)&&(th_PIm < th_max)&&(ph_PIm > pim_fid_b+180) && (ph_PIm < pim_fid_a+180)){

	if ((P_PIm > 0) && (P_PIm < 2.)&&(indtype==1)) ph_vs_th_pim_sim[3][int(P_PIm/0.4)]->Fill(th_PIm,ph_PIm-180,1);

	if (indtype==1) th_vs_p_pim_2_sim[3]->Fill(P_PIm,th_PIm,sigma);	
//th_vs_p  
if ((th_PIm < 36.152+3.69909e-05/(5.40783e-06*(P_PIm+0.31)+1.81169e-07)+62.)||((th_PIm>36.152+3.69909e-05/(5.40783e-06*(P_PIm+0.31)+1.81169e-07)+66.)&&(th_PIm <36.152+3.69909e-05/(5.40783e-06*(P_PIm+0.25)+1.81169e-07)+80.))){


 	if ((indtype==1)) hist_z_el_4_sim_1->Fill(z_PIm,sigma);

cuts_sim = true;
};//th_vs_p
};//end of the fiducial cut for sector4
};//end of the sector4
  
if ((ph_PIm >= 210) && (ph_PIm <=270)&&(PdHit_PIm!=17)) {

//fiducial cut sector 5  
if ((th_PIm > th_min)&&(th_PIm < th_max)&&(ph_PIm > pim_fid_b+240) && (ph_PIm < pim_fid_a+240)){

	if ((P_PIm > 0) && (P_PIm < 2.)&&(indtype==1)) ph_vs_th_pim_sim[4][int(P_PIm/0.4)]->Fill(th_PIm,ph_PIm-240,1);
	
	if (indtype==1) th_vs_p_pim_2_sim[4]->Fill(P_PIm,th_PIm,sigma);
//th_vs_p
if ((th_PIm<36.152+3.69909e-05/(5.40783e-06*(P_PIm)+1.81169e-07)-1.)||((th_PIm>36.152+3.69909e-05/(5.40783e-06*(P_PIm-0.01)+1.81169e-07)+3)&&(th_PIm<36.152+3.69909e-05/(5.40783e-06*(P_PIm+0.1)+1.81169e-07)+63.))){ 


 	if ((indtype==1)) hist_z_el_5_sim_1->Fill(z_PIm,sigma);

cuts_sim = true;
};//th_vs_p
};//end of the fiducial cut for sector5
}; //end of the sector5
    
 if ((ph_PIm >= 270) && (ph_PIm <=330)){
 
//fiducial cut sector 6  
if ((th_PIm > th_min)&&(th_PIm < th_max)&&(ph_PIm > pim_fid_b+300) && (ph_PIm < pim_fid_a+300)){
 
	if ((P_PIm > 0) && (P_PIm < 2.)&&(indtype==1)) ph_vs_th_pim_sim[5][int(P_PIm/0.4)]->Fill(th_PIm,ph_PIm-300,1);

	if (indtype==1) th_vs_p_pim_2_sim[5]->Fill(P_PIm,th_PIm,sigma);	
		
//th_vs_p
//if ((th_PIm <36.152+3.69909e-05/(5.40783e-06*(P_PIm+0.15)+1.81169e-07)+71.) || (th_PIm > 36.152+3.69909e-05/(5.40783e-06*(P_PIm+0.1)+1.81169e-07)+72.5 )){  
if (th_PIm <36.152+3.69909e-05/(5.40783e-06*(P_PIm+0.15)+1.81169e-07)+71.){  



 	if ((indtype==1)) hist_z_el_6_sim_1->Fill(z_PIm,sigma);

cuts_sim = true;
};//th_vs_p
};//end of the fiducial cut for sector6
}; //end of the sector6
  
   
  
  };//th_min for low momenta
   };//vertex cut
   };//end of pi- id beta vs p cut 
    return cuts_sim;
   };   
