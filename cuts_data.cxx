#include <TMath.h>
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include "cuts_data.h"
#include "global.h"
#include <iostream>



using namespace std;

ostringstream qqq1;


   bool cuts_data::Electron_cuts_data(){
   


   
   bool cuts_data;
   Float_t th_min,th_max,par1,par2,par3,fid_a,fid_b;
   Short_t i;
  
 

   
   Float_t ph_el_arr[3][6][18] = {{{1000.,22.,20.,20.,25.,23.,26.,22.,25.,27.,27.,29.,29.,29.,32.,30.,29.,1000.},
                             {1000.,27.,20.,18.,22.,20.,23.,21.,23.,23.,24.,24.,26.,31.,30.,27.,25.,1000.},
                             {1000.,25.,23.,23.,25.,25.,22.,24.,31.,28.,27.,30.,32.,36.,49.,41.,40.,1000.},
                             {1000.,18.,23.,20.,20.,20.,25.,27.,24.,24.,27.,29.,22.,26.,43.,44.,30.,1000.},
                             {1000.,25.,25.,28.,26.,25.,26.,22.,23.,27.,30.,30.,45.,31.,35.,45.,35.,1000.},
                             {1000.,26.,22.,25.,25.,30.,26.,27.,25.,35.,28.,37.,30.,30.,40.,40.,40.,1000.}},
			     
			     {{1000.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1000.},
                             {1000.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1000.},
		             {1000.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1000.},
	                     {1000.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1000.},
		             {1000.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1000.},
			     {1000.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1000.}},
			     
			     {{1000.,22.,20.,18.,21.,20.,20.,22.,25.,25.,29.,25.,25.,30.,30.,30.,25.,1000.},
                             {1000.,25.,23.,18.,19.,23.,24.,25.,29.,28.,26.,23.,26.,44.,39.,38.,35.,1000.},
			     {1000.,24.,22.,25.,28.,23.,23.,25.,24.,30.,33.,28.,28.,34.,35.,36.,34.,1000.},
			     {1000.,25.,23.,23.,23.,20.,25.,31.,24.,26.,28.,28.,23.,25.,32.,24.,32.,1000.},
			     {1000.,23.,21.,21.,23.,45.,45.,30.,25.,31.,30.,35.,47.,30.,45.,55.,40.,1000.},
			     {1000.,25.,22.,22.,23.,35.,33.,23.,26.,30.,40.,35.,27.,40.,75.,40.,40.,1000.}}};
			     
			     
     

 Float_t th_vs_seg_cc_arr[2][6][18] = {{{0., 14.6443, 16.0725, 17.9887, 19.8438, 21.8559, 24.0962, 25.7115, 27.3194, 29.4738, 32.1473, 34.4818, 37.6353, 40.3543, 42.726, 45.3277, 45.4403, 0.},
                             {0., 14.1306, 15.7974, 17.539, 19.6698, 21.751, 23.4988, 25.692, 27.3278, 29.6457, 32.138, 34.2722, 37.7139, 40.2945, 42.666, 45.0813, 45.5352, 0.},
                             {0., 14.3257, 15.8933, 17.7586, 19.7038, 21.5607, 23.4988, 25.8161, 27.3858, 29.3724, 32.1998, 34.3168, 37.6387, 40.2826, 42.5966, 44.9048, 45.4617, 0.},
                             {0., 14.6212, 16.0671, 17.8442, 19.7198, 21.9049, 24.0132, 25.6748, 27.4177, 29.1068, 31.7789, 34.4718, 37.4024, 39.9344, 42.6817, 45.2509, 45.6462, 0.},
                             {0., 14.6693, 17.0764, 18.3016, 19.7753, 21.8377, 23.9542, 25.8572, 28.4107, 30.5323, 33.1446, 35.2059, 37.4483, 39.1344, 42.3635, 45.238, 45.4197, 0.},
                             {0., 14.3506, 16.0391, 17.8922, 19.5839, 21.663, 23.8015, 25.4501, 27.3368, 29.6662, 32.0932, 34.436, 37.5931, 40.2924, 42.8242, 45.5393, 45.8316, 0.}},
			     
			     
			    {{0., 10.054, 11.4847, 13.2483, 14.854, 16.7384, 18.3456, 20.5912, 22.6226, 24.1823, 25.4867, 28.0661, 29.7871, 33.066, 35.8994, 38.209, 41.2385, 0.},
                             {0., 10.4103, 11.5865, 13.0199, 14.2593, 15.9883, 18.4077, 20.0705, 22.3853, 24.0439, 25.3731, 27.9493, 29.6573, 32.5285, 35.6867, 38.3033, 40.9372, 0.},
		             {0., 10.36, 11.53, 13.1973, 14.583, 16.505, 18.5437, 20.2397, 22.5306, 24.2456, 25.3623, 27.9259, 29.9135, 33.098, 36.0244, 38.5221, 41.0999, 0.},
	                     {0., 10.2087, 11.6431, 13.3227, 14.8077, 16.4797, 18.3064, 20.6596, 22.5267, 24.3291, 25.653, 27.7081, 29.4662, 32.9525, 35.4147, 38.0797, 40.7784, 0.},
		             {0., 10.4149, 11.4033, 13.3655, 14.9656, 16.6926, 18.4285, 20.5829, 22.3409, 24.2592, 26.507, 28.5396, 31.1582, 32.6398, 36.7438, 38.4544, 41.3839, 0.},
			     {0., 10.3843, 11.486, 13.1869, 14.8111, 16.5356, 17.894, 20.5821, 22.3433, 24.1164, 25.4934, 28.0616, 29.7482, 32.738, 35.9267, 38.5048, 40.8638, 0.}}};	
//for (i=0;i<18;i++){
//cout << ph_el[0][5][i] << "\n";
//};

 cuts_data = false; 
   
if ((LiveTime > 0.8) && (LiveTime <0.9) && (inclusive > 85000) &&(inclusive < 100000) && (elastic > 24000) && (elastic < 28000)){
 //if ((LiveTime > 0.8) && (LiveTime <0.9) && (inclusive > 195000) &&(inclusive < 210000) && (elastic > 50000) && (elastic < 58000)){
 
//calorimeter threshold cut + manually remove the first and last cc segments
if ((P_EL > 0.461)&&(segment!=0)&&(segment!=17)) {

   th_min= (11.7398+8.21504/(0.433327*P_EL+0.158076));
   par1=0.85+1.1*P_EL;
   par2=-62.8-30.*P_EL; 
   par3 = 0.0047*P_EL + 0.0079;

   fid_a = 41.3*pow((sin((th_EL-th_min)*par3)),(par1+par2/th_EL+1485./th_EL/th_EL)); 
   fid_b = -41.3*pow((sin((th_EL-th_min)*par3)),(par1+par2/th_EL+1485./th_EL/th_EL)); 
     
   th_max = 76.8617 -76.537*P_EL + 77.9387*P_EL*P_EL-28.389*P_EL*P_EL*P_EL;
  

	switch (sector) {
	
case 1 : 
 
// ectot vs p cut
if ((ECT/P_EL > -0.0557144*P_EL*P_EL+0.180989*P_EL+0.0630868) && (ECT/P_EL < 0.0218734*P_EL*P_EL-0.0607997*P_EL+0.395324)) {
 
	hist_ectot_sector1->Fill(P_EL,ECT/P_EL,1.); 
 
//~fid cuts in Cherenkov plane
if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {

//phi_cc matching
if (((ph_cc>0)&&(pmt_hit==1))||((ph_cc<0)&&(pmt_hit==-1))||(pmt_hit==0)){

//th_cc vs seg cut	
if ((theta_cc > th_vs_seg_cc_arr[1][0][segment])&&(theta_cc <th_vs_seg_cc_arr[0][0][segment])){

	
	th_cc_vs_seg_1->Fill(segment+1,theta_cc,1.);	
		
	h_cc_nphe_total_s1->Fill(theta_cc, ph_cc,1.);
	if (nphe > 50) h_cc_nphe_final_s1->Fill(theta_cc, ph_cc,1.);


	nphe_sector1->Fill(nphe,1.);
// geometrical cut on number of photoelectrons
if (norm_nphe_s1->GetBinContent(norm_nphe_s1->GetXaxis()->FindBin(theta_cc),norm_nphe_s1->GetYaxis()->FindBin(ph_cc))>0.7){
//if (norm_nphe_s1->GetBinContent(int((theta_cc+5.)*200./60.)+1,int((ph_cc+25.)*200./50.)+1) > 0.8) {

	nphe_sector1_after->Fill(nphe,1.);

//nphe cut after poisson fit
if (nphe > ph_el_arr[pmt_hit+1][0][segment]){

	if (pmt_hit == -1) ph_el_left[0][segment]->Fill(nphe,1.);
	if (pmt_hit == 0) ph_el_both[0][segment]->Fill(nphe,1.);
	if (pmt_hit == 1) ph_el_right[0][segment]->Fill(nphe,1.);

//vertex cut
if ((z_EL>-2.65) && (z_EL<1.85)){

if ((ph_EL >= 330) && (ph_EL <= 360)){

	ph_vs_th_1 -> Fill(th_EL,ph_EL-360,1.);

	if ((P_EL < 1.75999) && (P_EL > 0.4)){
	ph_vs_th_1pe[int((P_EL*100-40)/8)]->Fill(th_EL,ph_EL-360,1);
	};

//fiducial cut 
if ((th_EL > th_min) && (th_EL < th_max) && (ph_EL > fid_b+360) && (ph_EL < fid_a+360)){
	
	if ((P_EL < 1.75999) && (P_EL > 0.4)){
	ph_vs_th_1pe_fid[int((P_EL*100-40)/8)] ->Fill(th_EL,ph_EL-360,1.);
	};

	if ((W>1.3)&&(W<1.825)&&(Q2>0.4)&&(Q2<1.)&&(n_PIp==1)&&(n_P==1)) th_vs_p_e_2[0]->Fill(P_EL,th_EL,1.);

//	hist_z_el_1->Fill(z_EL,1.);

   cuts_data = true;
      
   }; //fiducial
   }; //second part of sector 1
   
if ((ph_EL >= 0) && (ph_EL <= 30)) {
  
	ph_vs_th_1 -> Fill(th_EL,ph_EL,1.); 
	if ((P_EL < 1.75999) &&(P_EL > 0.4)){
	ph_vs_th_1pe[int((P_EL*100-40)/8)] ->Fill(th_EL,ph_EL,1.);
	};

//fiducial cut 	
if ((th_EL > th_min) && (th_EL < th_max) && (ph_EL > fid_b) && (ph_EL < fid_a)){

	if ((P_EL < 1.75999) &&(P_EL > 0.4)){
	ph_vs_th_1pe_fid[int((P_EL*100-40)/8)] -> Fill(th_EL,ph_EL,1.);
	};
 
	if ((W>1.3)&&(W<1.825)&&(Q2>0.4)&&(Q2<1.)&&(n_PIp==1)) th_vs_p_e_2[0]->Fill(P_EL,th_EL,1.); 

//	hist_z_el_1->Fill(z_EL,1.);

   cuts_data = true; 
   	
   }; //fiducial
   }; //first part of sector 1
 };//vertex cut
 };//nphe cut after poisson fit
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

 
// ectot vs p cut
if ((ECT/P_EL > -0.0637287*P_EL*P_EL+0.199291*P_EL+0.0617286) && (ECT/P_EL < 0.0130875*P_EL*P_EL-0.0380366*P_EL+0.436297)) {

	hist_ectot_sector2->Fill(P_EL,ECT/P_EL,1.); 
  
//~fid cuts in Cherenkov plane
 if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {

//phi_cc matching
if (((ph_cc>0)&&(pmt_hit==1))||((ph_cc<0)&&(pmt_hit==-1))||(pmt_hit==0)){

//th_cc vs seg cut
if ((theta_cc >th_vs_seg_cc_arr[1][1][segment])&&(theta_cc <th_vs_seg_cc_arr[0][1][segment])){
	
		
	th_cc_vs_seg_2->Fill(segment+1,theta_cc,1.);
		
	h_cc_nphe_total_s2->Fill(theta_cc, ph_cc,1.);
	if (nphe > 50) h_cc_nphe_final_s2->Fill(theta_cc, ph_cc,1.);

	nphe_sector2->Fill(nphe,1.);
// geometrical cut on number of photoelectrons
if (norm_nphe_s2->GetBinContent(norm_nphe_s2->GetXaxis()->FindBin(theta_cc),norm_nphe_s2->GetYaxis()->FindBin(ph_cc))>0.65){

	nphe_sector2_after->Fill(nphe,1.);


//nphe cut after poisson fit	
if (nphe > ph_el_arr[pmt_hit+1][1][segment]){

	if (pmt_hit == -1) ph_el_left[1][segment]->Fill(nphe,1.);
	if (pmt_hit == 0) ph_el_both[1][segment]->Fill(nphe,1.);
	if (pmt_hit == 1) ph_el_right[1][segment]->Fill(nphe,1.);
		
	ph_vs_th_2 -> Fill(th_EL,ph_EL-60,1.);

	if ((P_EL < 1.75999) &&(P_EL > 0.4))  {
	ph_vs_th_2pe[int((P_EL*100-40)/8)] ->Fill(th_EL,ph_EL-60,1.);
	}; 

//fiducial cut
if ((th_EL > th_min) && (th_EL < th_max) && (ph_EL > fid_b+60) && (ph_EL < fid_a+60)&&(PdHit_EL!=16)){
	if ((P_EL < 1.75999) &&(P_EL > 0.4))  {
	ph_vs_th_2pe_fid[int((P_EL*100-40)/8)] ->Fill(th_EL,ph_EL-60,1.);
	};

	if ((W>1.3)&&(W<1.825)&&(Q2>0.4)&&(Q2<1.)&&(n_PIp==1)) th_vs_p_e_2[1]->Fill(P_EL,th_EL,1.);
	
//th_vs_p    	
if ((th_EL > (11.7398+8.21504/(0.433327*(P_EL+0.1)+0.158076))+18.3)||(th_EL < (11.7398+8.21504/(0.433327*(P_EL+0.1)+0.158076))+16.)){


//vertex cut
if ((z_EL>-2.35) && (z_EL<2.15)){

//	hist_z_el_2->Fill(z_EL,1.);

   cuts_data = true; 
    
    };//vertex cut
    };//th_vs_p     
   };//fiducial 
   };//nphe cut after poisson fit
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


 // ectot vs p cut 
if ((ECT/P_EL > -0.0450008*P_EL*P_EL+0.146354*P_EL+0.0768507) && (ECT/P_EL < 0.020841*P_EL*P_EL-0.0494426*P_EL+0.377883)){ 
	hist_ectot_sector3->Fill(P_EL,ECT/P_EL,1.);  
  
//~fid cuts in Cherenkov plane
if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {

//phi_cc matching
if (((ph_cc>0)&&(pmt_hit==1))||((ph_cc<0)&&(pmt_hit==-1))||(pmt_hit==0)){
	
//th_cc vs seg cut
if ((theta_cc > th_vs_seg_cc_arr[1][2][segment])&&(theta_cc < th_vs_seg_cc_arr[0][2][segment])){	
	
	
	th_cc_vs_seg_3->Fill(segment+1,theta_cc,1.);
		
	h_cc_nphe_total_s3->Fill(theta_cc, ph_cc,1.);
	if (nphe > 50) h_cc_nphe_final_s3->Fill(theta_cc, ph_cc,1.);


	nphe_sector3->Fill(nphe,1.); 
// geometrical cut on number of photoelectrons
if (norm_nphe_s3->GetBinContent(norm_nphe_s3->GetXaxis()->FindBin(theta_cc),norm_nphe_s3->GetYaxis()->FindBin(ph_cc))>0.7){

	nphe_sector3_after->Fill(nphe,1.);

//nphe cut after poisson fit	
if (nphe > ph_el_arr[pmt_hit+1][2][segment]){

	if (pmt_hit == -1) ph_el_left[2][segment]->Fill(nphe,1.);
	if (pmt_hit == 0) ph_el_both[2][segment]->Fill(nphe,1.);
	if (pmt_hit == 1) ph_el_right[2][segment]->Fill(nphe,1.);
	
	ph_vs_th_3 -> Fill(th_EL,ph_EL-120,1.);

	if ((P_EL < 1.75999) &&(P_EL > 0.4))  {
	ph_vs_th_3pe[int((P_EL*100-40)/8)] ->Fill(th_EL,ph_EL-120,1.);
	};
//fiducial cut
if ((th_EL > th_min) && (th_EL < th_max) && (ph_EL > fid_b+120) && (ph_EL < fid_a+120)&&(PdHit_EL!=44)){

	if ((P_EL < 1.75999) &&(P_EL > 0.4))  {
	ph_vs_th_3pe_fid[int((P_EL*100-40)/8)] ->Fill(th_EL,ph_EL-120,1.);
	}; 

	if ((W>1.3)&&(W<1.825)&&(Q2>0.4)&&(Q2<1.)&&(n_PIp==1)) th_vs_p_e_2[2]->Fill(P_EL,th_EL,1.);


//vertex cut
if ((z_EL>-2.25) && (z_EL< 2.25)){

//	hist_z_el_3->Fill(z_EL,1.);

  cuts_data = true; 

   };//vertex cut
   };//fiducial 
   };//nphe cut after poisson fit
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

 
// ectot vs p cut 
if ((ECT/P_EL > -0.0509612*P_EL*P_EL+0.163016*P_EL+0.0787492) && (ECT/P_EL < 0.0135551*P_EL*P_EL-0.042061*P_EL+0.387975)){
  
	hist_ectot_sector4->Fill(P_EL,ECT/P_EL,1.);

//~fid cuts in Cherenkov plane
if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {

//phi_cc matching
if (((ph_cc>0)&&(pmt_hit==1))||((ph_cc<0)&&(pmt_hit==-1))||(pmt_hit==0)){

//th_cc vs seg cut
if ((theta_cc >th_vs_seg_cc_arr[1][3][segment])&&(theta_cc <th_vs_seg_cc_arr[0][3][segment])){

	
	th_cc_vs_seg_4->Fill(segment+1,theta_cc,1.);
		
	h_cc_nphe_total_s4->Fill(theta_cc, ph_cc,1.);
	if (nphe > 50) h_cc_nphe_final_s4->Fill(theta_cc, ph_cc,1.);

	
	nphe_sector4->Fill(nphe,1.);
// geometrical cut on number of photoelectrons
if (norm_nphe_s4->GetBinContent(norm_nphe_s4->GetXaxis()->FindBin(theta_cc),norm_nphe_s4->GetYaxis()->FindBin(ph_cc))>0.65){

	nphe_sector4_after->Fill(nphe,1.);


//nphe cut after poisson fit 	
if (nphe > ph_el_arr[pmt_hit+1][3][segment]){

	if (pmt_hit == -1) ph_el_left[3][segment]->Fill(nphe,1.);
	if (pmt_hit == 0) ph_el_both[3][segment]->Fill(nphe,1.);
	if (pmt_hit == 1) ph_el_right[3][segment]->Fill(nphe,1.);

	ph_vs_th_4 -> Fill(th_EL,ph_EL-180,1);

	if ((P_EL < 1.75999) &&(P_EL > 0.4))  {
	ph_vs_th_4pe[int((P_EL*100-40)/8)] ->Fill(th_EL,ph_EL-180,1.);
		};
//fiducial cut
if ((th_EL > th_min) && (th_EL < th_max) && (ph_EL > fid_b+180) && (ph_EL < fid_a+180)){

	if ((P_EL < 1.75999) &&(P_EL > 0.4))  {
	ph_vs_th_4pe_fid[int((P_EL*100-40)/8)] ->Fill(th_EL,ph_EL-180,1.);
	};

	if ((W>1.3)&&(W<1.825)&&(Q2>0.4)&&(Q2<1.)&&(n_PIp==1)) th_vs_p_e_2[3]->Fill(P_EL,th_EL,1.); 


//vertex cut
if ((z_EL>-2.45) && (z_EL< 2.05)){

//	hist_z_el_4->Fill(z_EL,1.);

   cuts_data = true; 
   };//vertex cut 
   };//fiducial 
   };//nphe cut after poisson fit 
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

  
// ectot vs p cut 
if ((ECT/P_EL > -0.0584054*P_EL*P_EL+0.18591*P_EL+0.0577907) && (ECT/P_EL < 0.0180595*P_EL*P_EL-0.052512*P_EL+0.404114)&&(PdHit_EL != 17)){

	hist_ectot_sector5->Fill(P_EL,ECT/P_EL,1.);   

//~fid cuts in Cherenkov plane
if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {
	
//phi_cc matching
if (((ph_cc>0)&&(pmt_hit==1))||((ph_cc<0)&&(pmt_hit==-1))||(pmt_hit==0)){
	
//th_cc vs seg cut
if ((theta_cc >th_vs_seg_cc_arr[1][4][segment])&&(theta_cc <th_vs_seg_cc_arr[0][4][segment])){	
	
	
	th_cc_vs_seg_5->Fill(segment+1,theta_cc,1.);
			
	h_cc_nphe_total_s5->Fill(theta_cc, ph_cc,1.);
	if (nphe > 50) h_cc_nphe_final_s5->Fill(theta_cc, ph_cc,1.);
	
	nphe_sector5->Fill(nphe,1.);
// geometrical cut on number of photoelectrons
if (norm_nphe_s5->GetBinContent(norm_nphe_s5->GetXaxis()->FindBin(theta_cc),norm_nphe_s5->GetYaxis()->FindBin(ph_cc))>0.8){

	nphe_sector5_after->Fill(nphe,1.);


//nphe cut after poisson fit 	
if (nphe > ph_el_arr[pmt_hit+1][4][segment]){

	if (pmt_hit == -1) ph_el_left[4][segment]->Fill(nphe,1.);
	if (pmt_hit == 0) ph_el_both[4][segment]->Fill(nphe,1.);
	if (pmt_hit == 1) ph_el_right[4][segment]->Fill(nphe,1.);
	
	ph_vs_th_5 -> Fill(th_EL,ph_EL-240,1);
  
	if ((P_EL < 1.75999) &&(P_EL > 0.4))  {
	ph_vs_th_5pe[int((P_EL*100-40)/8)] ->Fill(th_EL,ph_EL-240,1.);
		};
//fiducial cut
if ((th_EL > th_min) && (th_EL < th_max) && (ph_EL > fid_b+240) && (ph_EL < fid_a+240)&&(PdHit_EL!=17)){

	if ((W>1.3)&&(W<1.825)&&(Q2>0.4)&&(Q2<1.)&&(n_PIp==1)) th_vs_p_e_2[4]->Fill(P_EL,th_EL,1.);
 
//th_vs_p
if ((th_EL > (11.7398+8.21504/(0.433327*(P_EL+0.1)+0.158076))+19.9)||(th_EL < (11.7398+8.21504/(0.433327*(P_EL+0.1)+0.158076))+17.6)){

	
	
	if ((P_EL < 1.75999) &&(P_EL > 0.4))  {
	ph_vs_th_5pe_fid[int((P_EL*100-40)/8)] ->Fill(th_EL,ph_EL-240,1.);
	};


//vertex cut
if ((z_EL>-2.9) && (z_EL< 1.6)){
  
//	hist_z_el_5->Fill(z_EL,1.);  
	
   cuts_data = true; 
   
   };//vertex cut
   };//th_vs_p
   };//fiducial 
   };//nphe cut after poisson fit    
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


// ectot vs p cut  
if ((ECT/P_EL > -0.03943*P_EL*P_EL+0.149493*P_EL+0.0789496) && (ECT/P_EL < 0.0123265*P_EL*P_EL-0.0290198*P_EL+0.380917)){

	hist_ectot_sector6->Fill(P_EL,ECT/P_EL,1.);   

//~fid cuts in Cherenkov plane
if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {

//phi_cc matching
if (((ph_cc>0)&&(pmt_hit==1))||((ph_cc<0)&&(pmt_hit==-1))||(pmt_hit==0)){

//th_cc vs seg cut 
if ((theta_cc >th_vs_seg_cc_arr[1][5][segment])&&(theta_cc <th_vs_seg_cc_arr[0][5][segment])){


	th_cc_vs_seg_6->Fill(segment+1,theta_cc,1.);
	
	h_cc_nphe_total_s6->Fill(theta_cc, ph_cc,1.);
	if (nphe > 50) h_cc_nphe_final_s6->Fill(theta_cc, ph_cc,1.);

	nphe_sector6->Fill(nphe,1.);	
// geometrical cut on number of photoelectrons
if (norm_nphe_s6->GetBinContent(norm_nphe_s6->GetXaxis()->FindBin(theta_cc),norm_nphe_s6->GetYaxis()->FindBin(ph_cc))>0.8){

	nphe_sector6_after->Fill(nphe,1.);
	
//nphe cut after poisson fit	
if (nphe > ph_el_arr[pmt_hit+1][5][segment]){

	if (pmt_hit == -1) ph_el_left[5][segment]->Fill(nphe,1.);
	if (pmt_hit == 0) ph_el_both[5][segment]->Fill(nphe,1.);
	if (pmt_hit == 1) ph_el_right[5][segment]->Fill(nphe,1.);

	ph_vs_th_6 -> Fill(th_EL,ph_EL-300,1); 

	if ((P_EL < 1.75999) &&(P_EL > 0.4))  {
	ph_vs_th_6pe[int((P_EL*100-40)/8)] ->Fill(th_EL,ph_EL-300,1.);
	};

//fiducial cut	
if ((th_EL > th_min) && (th_EL < th_max) && (ph_EL > fid_b+300) && (ph_EL < fid_a+300)){

	if ((P_EL < 1.75999) &&(P_EL > 0.4))  {
	ph_vs_th_6pe_fid[int((P_EL*100-40)/8)] ->Fill(th_EL,ph_EL-300,1.);
	};

	if ((W>1.3)&&(W<1.825)&&(Q2>0.4)&&(Q2<1.)&&(n_PIp==1)) th_vs_p_e_2[5]->Fill(P_EL,th_EL,1.); 



//vertex cut
if ((z_EL>-3.) && (z_EL<1.5)){

//	hist_z_el_6->Fill(z_EL,1.);
	
   cuts_data = true; 
   
   };//vertex  
   };//fiducial
   };//nphe cut after poisson fit
   };//geometrical cut on number of photoelectrons
   };//th_cc vs seg cut      
   };//phi_cc matching
   };//~fid cuts in Cherenkov plane
   };//~fid cuts in Cherenkov plane
   };//~fid cuts in Cherenkov plane
   };//~fid cuts in Cherenkov plane   
 
 };//ectot vs p cut
 break;      
 
   
   };//end of switch
  };//end of calorimeter threshold cut
  };//end of quality check cuts
   

   
   
   return cuts_data;
   
   };




///////////////////////////////////////////////////////////////////////////////////////////////

   bool cuts_data::Proton_cuts_data(){
       
   bool cuts_data;
   Float_t m_p,p_fid_a,p_fid_b,th_min,th_max;	
    
   m_p=0.938272;   

   th_min = 12.;
   th_max = 60.;
  p_fid_a = 25.*(1-exp(-1.*0.12*(th_P-10.)))-2.5;
  p_fid_b = -25.*(1-exp(-1.*0.12*(th_P-10.)))+2.5; 

   cuts_data = false; 

//Proton id beta_vs_p cut
if((n_P == 1)&& (beta_P < 0.958174*P_P/sqrt(m_p*m_p+0.921203*P_P*P_P-0.0289732)+ 0.0383292) && (beta_P > 0.97261*P_P/sqrt(m_p*m_p +0.89435*P_P*P_P -0.199588) -0.0775468)&&(PdHit_P !=48)){

//----------PROTON ENERGY LOSS----------------------------
P_P = P_P + delta_mom_p_skor;
//--------------------------------------------------------

//minimum proton momentum & mimimum proton beta cuts
if ((P_P > 0.235)&&(beta_P > 0.19)){

//vertex cut
if ((z_P > -4.8) && (z_P < 4.)){

if ((ph_P >= 330)&& (ph_P <= 360) ) {

//fiducial sector1
if ((th_P > th_min)&&(th_P < th_max)&&(ph_P > p_fid_b+360) && (ph_P < p_fid_a+360)){
  
	if ((W>1.3)&&(P_P > 0.2) && (P_P < 1.4)&&(n_PIp==1)) ph_th_p_1[int((P_P-0.2)/0.2)]->Fill(th_P,ph_P-360,1);
	if ((W>1.3)&&(P_P > 0.2) && (P_P < 1.4)&&(n_PIp==1)) ph_vs_th_p_1_w -> Fill(th_P,ph_P-360,1);  

	if ((W>1.3)&&(n_PIp==1)) th_vs_p_p_2[0]->Fill(P_P,th_P,1.);

//th_vs_p1 
if (!((th_P<pow((P_P-0.03-0.304992),(0.0758186))*91.5643-48.2057 + 2. )&&!(th_P<pow((P_P-0.08-0.304992),(0.0758186))*91.5643-48.2057 +1.5))){
//th_vs_p2 
if (!((th_P <(304.23*(P_P+0.2)*(P_P+0.2)*(P_P+0.2) -255.798*(P_P+0.2)*(P_P+0.2)+497.462*(P_P+0.2) +38.0385)*exp(-1.6*(P_P+0.2))-104.)&&(th_P>(304.23*(P_P-0.12)*(P_P-0.12)*(P_P-0.12) -255.798*(P_P-0.12)*(P_P-0.12)+497.462*(P_P-0.12) +38.0385)*exp(-1.6*(P_P-0.12))-103.))){

//	hist_z_el_1->Fill(z_P,1.);
	
cuts_data = true;

};//th_vs_p2
};//th_vs_p1
};//end of fiducial cut for the first part of the first sector
  

	ph_vs_th_p_1 -> Fill(th_P,ph_P-360,1);
 
};//end of the first part of the first sector
 
if ((ph_P >= 0) && (ph_P <= 30)) {

//fiducial sector1
if ((th_P > th_min)&&(th_P < th_max)&&(ph_P > p_fid_b) && (ph_P < p_fid_a)){

 	if ((W>1.3)&&(P_P > 0.2) && (P_P < 1.4)&&(n_PIp==1)) ph_th_p_1[int((P_P-0.2)/0.2)]->Fill(th_P,ph_P,1);
	if ((W>1.3)&&(P_P > 0.2) && (P_P < 1.4)&&(n_PIp==1)) ph_vs_th_p_1_w -> Fill(th_P,ph_P,1);

	if ((W>1.3)&&(n_PIp==1)) th_vs_p_p_2[0]->Fill(P_P,th_P,1.); 

//th_vs_p1 
if (!((th_P<pow((P_P-0.03-0.304992),(0.0758186))*91.5643-48.2057 + 2. )&&!(th_P<pow((P_P-0.08-0.304992),(0.0758186))*91.5643-48.2057 +1.5))){
//th_vs_p2 
if (!((th_P <(304.23*(P_P+0.2)*(P_P+0.2)*(P_P+0.2) -255.798*(P_P+0.2)*(P_P+0.2)+497.462*(P_P+0.2) +38.0385)*exp(-1.6*(P_P+0.2))-104.)&&(th_P>(304.23*(P_P-0.12)*(P_P-0.12)*(P_P-0.12) -255.798*(P_P-0.12)*(P_P-0.12)+497.462*(P_P-0.12) +38.0385)*exp(-1.6*(P_P-0.12))-103.))){
 
//	hist_z_el_1->Fill(z_P,1.);

cuts_data = true;

};//th_vs_p2
};//th_vs_p1
};//end of fiducial cut for the second part of the first sector

	ph_vs_th_p_1 -> Fill(th_P,ph_P,1);
 
};//end of the second part of the first sector
 
 
 
if ((ph_P >= 30) && (ph_P <=90)&&(PdHit_P!=16)){ 

//fiducial sector2
if ((th_P > th_min)&&(th_P < th_max)&&(ph_P > p_fid_b+60) && (ph_P < p_fid_a+60)){
 
	if ((W>1.3)&&(P_P > 0.2) && (P_P < 1.4)&&(n_PIp==1)) ph_th_p_2[int((P_P-0.2)/0.2)]->Fill(th_P,ph_P-60,1);
	if ((W>1.3)&&(P_P > 0.2) && (P_P < 1.4)&&(n_PIp==1)) ph_vs_th_p_2_w -> Fill(th_P,ph_P-60,1); 
 
	if ((W>1.3)&&(n_PIp==1)) th_vs_p_p_2[1]->Fill(P_P,th_P,1.); 
	
//th_vs_p1
if (!((th_P<pow((P_P-0.03-0.304992),(0.0758186))*91.5643-48.2057 + 2. )&&!(th_P<pow((P_P-0.08-0.304992),(0.0758186))*91.5643-48.2057 +1.5))){
//th_vs_p2
if (!((th_P <pow((P_P-0.12-0.304992),(0.0758186))*91.5643-48.2057 +1.1)&&!(th_P<pow((P_P-0.16-0.304992),(0.0758186))*91.5643-48.2057))){
//th_vs_p3
if (!((th_P <26.5087*(P_P+0.04)*(P_P+0.04)*(P_P+0.04) -116.557*(P_P+0.04)*(P_P+0.04)+ 175.167*(P_P+0.04)-61.7717+1.5)&&(th_P>26.5087*(P_P-0.03)*(P_P-0.03)*(P_P-0.03) -116.557*(P_P-0.03)*(P_P-0.03)+ 175.167*(P_P-0.03)-61.7717-1.8))){

//	hist_z_el_2->Fill(z_P,1.); 

cuts_data = true;

};//th_vs_p3
};//th_vs_p2
};//th_vs_p1
};//end of the fiducial cut for sector2

 
	ph_vs_th_p_2 -> Fill(th_P,ph_P-60,1);
};//end of the sector2

 
if ((ph_P >=90) && (ph_P <=150)&&(PdHit_P!=44)) {

//fiducial sector3
if ((th_P > th_min)&&(th_P < th_max)&&(ph_P > p_fid_b+120) && (ph_P < p_fid_a+120)){

	if ((W>1.3)&&(P_P > 0.2) && (P_P < 1.4)&&(n_PIp==1)) ph_th_p_3[int((P_P-0.2)/0.2)]->Fill(th_P,ph_P-120,1);
	if ((W>1.3)&&(P_P > 0.2) && (P_P < 1.4)&&(n_PIp==1)) ph_vs_th_p_3_w -> Fill(th_P,ph_P-120,1); 

	if ((W>1.3)&&(n_PIp==1)) th_vs_p_p_2[2]->Fill(P_P,th_P,1.);

//th_vs_p1	
if (!((th_P<pow((P_P-0.03-0.304992),(0.0758186))*91.5643-48.2057 + 2. )&&!(th_P<pow((P_P-0.08-0.304992),(0.0758186))*91.5643-48.2057 +1.5))){

//th_vs_p2
if (!((th_P <pow((P_P-0.12-0.304992),(0.0758186))*91.5643-48.2057 +1.1)&&!(th_P<pow((P_P-0.16-0.304992),(0.0758186))*91.5643-48.2057))){

//	hist_z_el_3->Fill(z_P,1.); 

cuts_data = true;

};//th_vs_p2
};//th_vs_p1

};//end of the fiducial cut for sector3

	ph_vs_th_p_3 -> Fill(th_P,ph_P-120,1);
	
};//end of the sector3
 
if ((ph_P >= 150) && (ph_P <= 210)) {

//fiducial sector4
if ((th_P > th_min)&&(th_P < th_max)&&(ph_P > p_fid_b+180) && (ph_P < p_fid_a+180)){

	if ((W>1.3)&&(P_P > 0.2) && (P_P < 1.4)&&(n_PIp==1)) ph_th_p_4[int((P_P-0.2)/0.2)]->Fill(th_P,ph_P-180,1);
	if ((W>1.3)&&(P_P > 0.2) && (P_P < 1.4)&&(n_PIp==1)) ph_vs_th_p_4_w -> Fill(th_P,ph_P-180,1);

	if ((W>1.3)&&(n_PIp==1)) th_vs_p_p_2[3]->Fill(P_P,th_P,1.);

//th_vs_p1	
if (!((th_P<pow((P_P-0.03-0.304992),(0.0758186))*91.5643-48.2057 + 2. )&&!(th_P<pow((P_P-0.08-0.304992),(0.0758186))*91.5643-48.2057 +1.5))){
//th_vs_p2
if (!((th_P <pow((P_P-0.12-0.304992),(0.0758186))*91.5643-48.2057 +1.1)&&!(th_P<pow((P_P-0.16-0.304992),(0.0758186))*91.5643-48.2057))){

//	hist_z_el_4->Fill(z_P,1.);

cuts_data = true;

};//th_vs_p2
};//th_vs_p1

};//end of the fiducial cut for sector4

ph_vs_th_p_4 -> Fill(th_P,ph_P-180,1);

};//end of the sector4


if ((ph_P >= 210) && (ph_P <=270)&&(PdHit_P!=17)) {

//fiducial sector5
if ((th_P > th_min)&&(th_P < th_max)&&(ph_P > p_fid_b+240) && (ph_P < p_fid_a+240)){

	if ((W>1.3)&&(P_P > 0.2) && (P_P < 1.4)&&(n_PIp==1)) ph_th_p_5[int((P_P-0.2)/0.2)]->Fill(th_P,ph_P-240,1);
	if ((W>1.3)&&(P_P > 0.2) && (P_P < 1.4)&&(n_PIp==1)) ph_vs_th_p_5_w -> Fill(th_P,ph_P-240,1); 

	if ((W>1.3)&&(n_PIp==1)) th_vs_p_p_2[4]->Fill(P_P,th_P,1.);
	
//th_vs_p1
if (!((th_P<pow((P_P-0.03-0.304992),(0.0758186))*91.5643-48.2057 + 2. )&&!(th_P<pow((P_P-0.08-0.304992),(0.0758186))*91.5643-48.2057 +1.5))){
//th_vs_p2
if (!((th_P <pow((P_P -0.321436),(0.0704348))*88.0419-46.9342-0.5)&&!(th_P<pow((P_P -0.12 -0.321436),(0.0704348))*88.0419-46.9342-2.5))){
//th_vs_p3
if (!((th_P <31.2482*(P_P+0.045)*(P_P+0.045)*(P_P+0.045) -135.817*(P_P+0.045)*(P_P+0.045)+ 198.038*(P_P+0.045)-66.968+1.5)&&(th_P>31.2482*(P_P-0.01)*(P_P-0.01)*(P_P-0.01) -135.817*(P_P-0.01)*(P_P-0.01)+ 198.038*(P_P-0.01)-66.968-2.5))){

//	hist_z_el_5->Fill(z_P,1.);

 cuts_data = true;

};//th_vs_p3
};//th_vs_p2
};//th_vs_p1
};//end of the fiducial cut for sector5

//	ph_vs_th_p_5 -> Fill(th_P,ph_P-240,1);
 
 };//end of the sector5

 
if ((ph_P >= 270) && (ph_P <=330)) {

//fiducial sector6
if ((th_P > th_min)&&(th_P < th_max)&&(ph_P > p_fid_b+300) && (ph_P < p_fid_a+300)){

	if ((W>1.3)&&(P_P > 0.2) && (P_P < 1.4)&&(n_PIp==1)) ph_th_p_6[int((P_P-0.2)/0.2)]->Fill(th_P,ph_P-300,1);
	if ((W>1.3)&&(P_P > 0.2) && (P_P < 1.4)&&(n_PIp==1)) ph_vs_th_p_6_w -> Fill(th_P,ph_P-300,1); 

	if ((W>1.3)&&(n_PIp==1)) th_vs_p_p_2[5]->Fill(P_P,th_P,1.);

//th_vs_p1
if (!((th_P<pow((P_P-0.03-0.304992),(0.0758186))*91.5643-48.2057 + 2. )&&!(th_P<pow((P_P-0.08-0.304992),(0.0758186))*91.5643-48.2057 +1.5))){
//th_vs_p2
if (!((th_P <pow((P_P-0.12-0.304992),(0.0758186))*91.5643-48.2057 +1.1)&&!(th_P<pow((P_P-0.16-0.304992),(0.0758186))*91.5643-48.2057))){

//	hist_z_el_6->Fill(z_P,1.);

cuts_data = true;

};//th_vs_p2
};//th_vs_p1
};//end of the fiducial cut for sector6

	ph_vs_th_p_6 -> Fill(th_P,ph_P-300,1);
};//end of the sector6
 
};//vertex cut 
};// end of minimum proton momentum & mimimum proton beta cuts
};//end proton id beta vs p cut
 
return cuts_data;
};
    
    
/////////////////////////////
    
    bool cuts_data::PIp_cuts_data(){
       
   bool cuts_data;
   Float_t m_pip,pip_fid_a,pip_fid_b,th_min,th_max;
   m_pip = 0.13957;
   
   th_min = 12.;
   th_max = 120.;
  pip_fid_a = 25.*(1-exp(-1.*0.12*(th_PIp-10.)))-2.5;
  pip_fid_b = -25.*(1-exp(-1.*0.12*(th_PIp-10.)))+2.5;


   cuts_data = false; 

//pi+ id beta vs p cut   
if((n_PIp == 1)&&(beta_PIp < 0.838148*P_PIp/sqrt(m_pip*m_pip+0.793572*P_PIp*P_PIp -0.00291347)+ 0.105962) && (beta_PIp > 0.819388*P_PIp/sqrt(m_pip*m_pip +0.600785*P_PIp*P_PIp -0.00885673) -0.107038)&&(PdHit_PIp !=48)){
   
//vertex cut
if ((z_PIp>-4.8) && (z_PIp<4.)){
      
if ((ph_PIp >= 330) && (ph_PIp <=360)){
   
	ph_vs_th_pip_1 -> Fill(th_PIp,ph_PIp-360,1);

//fiducial sector 1 
if ((th_PIp > th_min)&&(th_PIp < th_max)&&(ph_PIp > pip_fid_b+360) && (ph_PIp < pip_fid_a+360)){
 
	if ((P_PIp > 0) && (P_PIp < 1.2)) ph_th_pip_1[int(P_PIp/0.2)]->Fill(th_PIp,ph_PIp-360,1); 

	th_vs_p_pip_2[0]->Fill(P_PIp,th_PIp,1.);
	
//th_vs_p1
if (!((th_PIp<(304.23*(P_PIp+0.15)*(P_PIp+0.15)*(P_PIp+0.15) -255.798*(P_PIp+0.15)*(P_PIp+0.15)+497.462*(P_PIp+0.15) +38.0385)*exp(-1.85*(P_PIp+0.15)) +5.5 )&&!(th_PIp<(304.23*(P_PIp+0.15)*(P_PIp+0.15)*(P_PIp+0.15) -255.798*(P_PIp+0.15)*(P_PIp+0.15)+497.462*(P_PIp+0.15) +38.0385)*exp(-1.85*(P_PIp+0.15)) -3))){
//th_vs_p2
if (!((th_PIp<(pow((P_PIp-0.0575818),( 0.075643))*238.248-115.039)*exp(-0.5*P_PIp)-0.1)&&!(th_PIp<(pow((P_PIp-0.103718),(0.0703664))*252.822-133.024)*exp(-0.5*P_PIp)+0.1))){
//th_vs_p3
if (!((th_PIp<(304.23*(P_PIp+0.25)*(P_PIp+0.25)*(P_PIp+0.25) -255.798*(P_PIp+0.25)*(P_PIp+0.25)+497.462*(P_PIp+0.25) +38.0385)*exp(-1.6*(P_PIp+0.25))-104)&&!(th_PIp<(304.23*(P_PIp-0.12)*(P_PIp-0.12)*(P_PIp-0.12) -255.798*(P_PIp-0.12)*(P_PIp-0.12)+497.462*(P_PIp-0.12) +38.0385)*exp(-1.6*(P_PIp-0.12))-103))){	

//	hist_z_el_1->Fill(z_PIp,1.);
		
   cuts_data = true;

};//th_vs_p3
};//th_vs_p2
};//th_vs_p1
};//end of the fiducial cut for first part of sector1
}; //end of the first part of sector1

  
if ((ph_PIp >= 0) && (ph_PIp <=30)){

	ph_vs_th_pip_1 -> Fill(th_PIp,ph_PIp,1);

//fiducial sector 1
if ((th_PIp > th_min)&&(th_PIp < th_max)&&(ph_PIp > pip_fid_b) && (ph_PIp < pip_fid_a)){
 
	if ((P_PIp > 0) && (P_PIp < 1.2)) ph_th_pip_1[int(P_PIp/0.2)]->Fill(th_PIp,ph_PIp,1); 
 
	th_vs_p_pip_2[0]->Fill(P_PIp,th_PIp,1.);
	 
//th_vs_p1
if (!((th_PIp<(304.23*(P_PIp+0.15)*(P_PIp+0.15)*(P_PIp+0.15) -255.798*(P_PIp+0.15)*(P_PIp+0.15)+497.462*(P_PIp+0.15) +38.0385)*exp(-1.85*(P_PIp+0.15)) +5.5 )&&!(th_PIp<(304.23*(P_PIp+0.15)*(P_PIp+0.15)*(P_PIp+0.15) -255.798*(P_PIp+0.15)*(P_PIp+0.15)+497.462*(P_PIp+0.15) +38.0385)*exp(-1.85*(P_PIp+0.15)) -3))){
//th_vs_p2
if (!((th_PIp<(pow((P_PIp-0.0575818),( 0.075643))*238.248-115.039)*exp(-0.5*P_PIp)-0.1)&&!(th_PIp<(pow((P_PIp-0.103718),(0.0703664))*252.822-133.024)*exp(-0.5*P_PIp)+0.1))){
//th_vs_p3
if (!((th_PIp<(304.23*(P_PIp+0.25)*(P_PIp+0.25)*(P_PIp+0.25) -255.798*(P_PIp+0.25)*(P_PIp+0.25)+497.462*(P_PIp+0.25) +38.0385)*exp(-1.6*(P_PIp+0.25))-104)&&!(th_PIp<(304.23*(P_PIp-0.12)*(P_PIp-0.12)*(P_PIp-0.12) -255.798*(P_PIp-0.12)*(P_PIp-0.12)+497.462*(P_PIp-0.12) +38.0385)*exp(-1.6*(P_PIp-0.12))-103))){

//	hist_z_el_1->Fill(z_PIp,1.); 

cuts_data = true;
};//th_vs_p3
};//th_vs_p2
};//th_vs_p1

};//end of the fiducial cut for second part of sector1
}; //end of the second part of sector1
  
  
if ((ph_PIp >= 30) && (ph_PIp <=90)&&(PdHit_PIp !=16)) {
  
	ph_vs_th_pip_2 -> Fill(th_PIp,ph_PIp-60,1);

//fiducial sector 2
if ((th_PIp > th_min)&&(th_PIp < th_max)&&(ph_PIp > pip_fid_b+60) && (ph_PIp < pip_fid_a+60)){

	if ((P_PIp > 0) && (P_PIp < 1.2)) ph_th_pip_2[int(P_PIp/0.2)]->Fill(th_PIp,ph_PIp-60,1);

	th_vs_p_pip_2[1]->Fill(P_PIp,th_PIp,1.);
	
//th_vs_p1
if (!((th_PIp<pow((P_PIp-0.415068),(0.226449))*48.7564 + 2.79478-1.)&&!(th_PIp<pow((P_PIp-0.449975),( 0.315164 ))*36.608 +  9.74262-1.))){
//th_vs_p2
if (!((th_PIp<(387.289*P_PIp*P_PIp*P_PIp -758.466*P_PIp*P_PIp+ 842.881*P_PIp-299.953+15.)*exp(-2*P_PIp))&&!(th_PIp<(387.289*(P_PIp+0.03)*(P_PIp+0.03)*(P_PIp+0.03) -758.466*(P_PIp+0.03)*(P_PIp+0.03)+ 842.881*(P_PIp+0.03)-299.953-15.)*exp(-2*(P_PIp+0.03))-1.5))){

//	hist_z_el_2->Fill(z_PIp,1.);
	
cuts_data = true;
};//th_vs_p2
};//th_vs_p1
};//end of the fiducial cut for sector2
};//end of the sector2
 
  
  
if ((ph_PIp >= 90) && (ph_PIp <=150)&&(PdHit_PIp!=44)) {

	ph_vs_th_pip_3 -> Fill(th_PIp,ph_PIp-120,1);

//fiducial sector 3
if ((th_PIp > th_min)&&(th_PIp < th_max)&&(ph_PIp > pip_fid_b+120) && (ph_PIp < pip_fid_a+120)){

	if ((P_PIp > 0) && (P_PIp < 1.2)) ph_th_pip_3[int(P_PIp/0.2)]->Fill(th_PIp,ph_PIp-120,1);

	th_vs_p_pip_2[2]->Fill(P_PIp,th_PIp,1.);

//th_vs_p1
if (!((th_PIp<(10000*P_PIp*P_PIp*P_PIp-3607.41*P_PIp*P_PIp+ 1725.72*P_PIp-10.6776)*exp(-4.7*P_PIp))&&!(th_PIp<(10000*P_PIp*P_PIp*P_PIp-4505.62*P_PIp*P_PIp+  2056.24  *P_PIp -77.4077 +5.)*exp(-4.7*P_PIp)))){
//th_vs_p2
if (!((th_PIp<pow((P_PIp-0.416536),(0.108376))*67.4593-21.4374)&&!(th_PIp<pow((P_PIp -0.454898),( 0.289291))* 35.7267+6.65908+1.5))){

//	hist_z_el_3->Fill(z_PIp,1.);
	
cuts_data = true;
};//th_vs_p2
};//th_vs_p1
};//end of the fiducial cut for sector3
};  //end of the sector3

  
if ((ph_PIp >= 150) && (ph_PIp <=210)){
 
	ph_vs_th_pip_4 -> Fill(th_PIp,ph_PIp-180,1);

//fiducial sector 4 
if ((th_PIp > th_min)&&(th_PIp < th_max)&&(ph_PIp > pip_fid_b+180) && (ph_PIp < pip_fid_a+180)){

	if ((P_PIp > 0) && (P_PIp < 1.2)) ph_th_pip_4[int(P_PIp/0.2)]->Fill(th_PIp,ph_PIp-180,1);

	th_vs_p_pip_2[3]->Fill(P_PIp,th_PIp,1.);
//th_vs_p1	
if (!((th_PIp<(304.23*(P_PIp+0.165)*(P_PIp+0.165)*(P_PIp+0.165) -255.798*(P_PIp+0.165)*(P_PIp+0.165)+497.462*(P_PIp+0.165) +38.0385)*exp(-1.85*(P_PIp+0.165)) +5.)&&!(th_PIp<(304.23*(P_PIp+0.18)*(P_PIp+0.18)*(P_PIp+0.18) -255.798*(P_PIp+0.18)*(P_PIp+0.18)+497.462*(P_PIp+0.18) +38.0385)*exp(-1.85*(P_PIp+0.18)) - 1.))){
//th_vs_p2
if (!((th_PIp<(1600*(P_PIp+0.03)*(P_PIp+0.03)*(P_PIp+0.03) -1068.36*(P_PIp+0.03)*(P_PIp+0.03)+ 775.016*(P_PIp+0.03)-1.13034)*exp(-2.75*(P_PIp+0.03)))&&!(th_PIp<(pow((P_PIp-0.103718),(0.0703664))*252.822-133.024)*exp(-0.45*P_PIp)-7.))){
//th_vs_p3
if (!((th_PIp<pow((P_PIp-0.416536),(0.108376))*67.4593-21.4374)&&!(th_PIp<pow((P_PIp -0.454898),( 0.289291))* 35.7267+6.65908+1.5))){

//	hist_z_el_4->Fill(z_PIp,1.);
	
cuts_data = true;
};//th_vs_p3
};//th_vs_p2
};//th_vs_p1
};//end of the fiducial cut for sector4
}; //end of the sector4
  
  
if ((ph_PIp >= 210) && (ph_PIp <=270)&&(PdHit_PIp!=17)) {

	ph_vs_th_pip_5 -> Fill(th_PIp,ph_PIp-240,1);

//fiducial sector 5 
if ((th_PIp > th_min)&&(th_PIp < th_max)&&(ph_PIp > pip_fid_b+240) && (ph_PIp < pip_fid_a+240)){

	if ((P_PIp > 0) && (P_PIp < 1.2)) ph_th_pip_5[int(P_PIp/0.2)]->Fill(th_PIp,ph_PIp-240,1);

	th_vs_p_pip_2[4]->Fill(P_PIp,th_PIp,1.);
	
//th_vs_p1
if (!((th_PIp<(525.498*(P_PIp+0.03)*(P_PIp+0.03)*(P_PIp+0.03) -1284.98*(P_PIp+0.03)*(P_PIp+0.03)+1460.67*(P_PIp+0.03)-499.999)*exp(-1.94*(P_PIp+0.03)))&&!(th_PIp<(525.498*(P_PIp-0.02)*(P_PIp-0.02)*(P_PIp-0.02) -1284.98*(P_PIp-0.02)*(P_PIp-0.02)+1460.67*(P_PIp-0.02)-499.999)*exp(-1.94*(P_PIp-0.02)) -4.7))){
//th_vs_p2
if (!((th_PIp<(304.23*(P_PIp)*(P_PIp)*(P_PIp) -255.798*(P_PIp)*(P_PIp)+497.462*(P_PIp) +38.0385)*exp(-1.85*(P_PIp)))&&!(th_PIp<(304.23*(P_PIp+0.03)*(P_PIp+0.03)*(P_PIp+0.03) -255.798*(P_PIp+0.03)*(P_PIp+0.03)+497.462*(P_PIp+0.03) +38.0385)*exp(-1.85*(P_PIp+0.03)) -11.))){
//th_vs_p3
if (!((th_PIp<pow((P_PIp-0.304992),(0.0758186))*91.5643-48.2057 - 1.)&&!(th_PIp<pow((P_PIp -0.36848),( 0.0864219))*70.4769  -34.9998+ 1.5))){

//	hist_z_el_5->Fill(z_PIp,1.);
	
cuts_data = true;
};//th_vs_p3
};//th_vs_p2
};//th_vs_p1
};//end of the fiducial cut for sector5
}; //end of the sector5
  
  
if ((ph_PIp >= 270) && (ph_PIp <=330)){

	ph_vs_th_pip_6 -> Fill(th_PIp,ph_PIp-300,1);

//fiducial sector 6 
if ((th_PIp > th_min)&&(th_PIp < th_max)&&(ph_PIp > pip_fid_b+300) && (ph_PIp < pip_fid_a+300)){

	if ((P_PIp > 0) && (P_PIp < 1.2)) ph_th_pip_6[int(P_PIp/0.2)]->Fill(th_PIp,ph_PIp-300,1);

	th_vs_p_pip_2[5]->Fill(P_PIp,th_PIp,1.);
		
//th_vs_p1
if (!((th_PIp<pow((P_PIp-0.05-0.0942469),( 0.0582707))*114.358-50 -0.5)&&!(th_PIp<pow((P_PIp-0.05-0.126994),( 0.0706829))* 110.073-50+2.))){
//th_vs_p2
if (!((th_PIp<pow((P_PIp-0.416536),(0.108376))*67.4593-21.4374-1.)&&!(th_PIp<pow((P_PIp -0.454098),(0.0912936))*58.2946-20.4843+1.5))){

//	hist_z_el_6->Fill(z_PIp,1.);
		
cuts_data = true;
};//th_vs_p2
};//th_vs_p1
};//end of the fiducial cut for sector6
}; //end of the sector6
  
}; //vertex cut
};//end of pi+ id beta vs p cut
return cuts_data;
};
   
   
   /////////////////////////////////////

  
    bool cuts_data::PIm_cuts_data(){
       
   bool cuts_data;
   Float_t m_pim,th_min,th_max,th_edge,par1,par2,par3,pim_fid_a,pim_fid_b, th_min_lowp;
   
   m_pim = 0.13957; 
   
  th_min=(11.09+8./(0.472*(P_PIm-0.03)+0.117))-1.;
  th_min_lowp = 30.+5.24894e-05/(5.71075e-05*(P_PIm+0.004)*(P_PIm+0.004)+4.44089e-16)+3.;
   
  if (P_PIm<0.31) th_max = 140.;  
  if (P_PIm>0.31) th_max =  254.914-506.82*P_PIm + 424.853*P_PIm*P_PIm-129.339*P_PIm*P_PIm*P_PIm;
  
  th_edge = th_max- 4./3.*(th_max-th_min)/2.;
    
  par1=0.61+1.18*P_PIm;
  par2=-59.2-35.3*P_PIm; 
  par3=17.2*P_PIm-11.9*P_PIm*P_PIm-1.;
   
   pim_fid_a=+23.5*pow((sin((th_PIm-th_min)*0.015)),(par1+par2/th_PIm+1400./th_PIm/th_PIm))+par3;
   pim_fid_b=-23.5*pow((sin((th_PIm-th_min)*0.015)),(par1+par2/th_PIm+1400./th_PIm/th_PIm))-par3;
    
  if (th_PIm > th_edge) pim_fid_a=+23.5*pow((sin((th_edge-th_min)*0.015)),(par1+par2/th_edge+1400./th_edge/th_edge))+par3;
  if (th_PIm > th_edge) pim_fid_b=-23.5*pow((sin((th_edge-th_min)*0.015)),(par1+par2/th_edge+1400./th_edge/th_edge))-par3; 

   cuts_data = false; 

//pi- id beta vs p cut   
if((n_PIm == 1)&& (beta_PIm < 0.779984*P_PIm/sqrt(m_pip*m_pip+0.700052*P_PIm*P_PIm-0.00506451 )+ 0.105877) && (beta_PIm > 0.854914*P_PIm/sqrt(m_pip*m_pip +0.629965*P_PIm*P_PIm -0.00757789)-0.113261 )&&(PdHit_PIm!=48)){ 

//vertex cut
if ((z_PIm>-4.8) && (z_PIm<4.)){

//th_min for low momenta  
if (((P_PIm<0.3)&&(th_PIm > th_min_lowp))||(P_PIm>0.3)){
   
if ((ph_PIm >= 330) && (ph_PIm <=360)){

	ph_th_pim_all_p[0] -> Fill(th_PIm,ph_PIm-360,1);

	if ((P_PIm > 0.1) && (P_PIm < 1.6)) ph_vs_th_pim[0][int((P_PIm-0.1)/0.1)]->Fill(th_PIm,ph_PIm-360,1);	
	
//fiducial cut sector 1
if ((th_PIm > th_min)&&(th_PIm < th_max)&&(ph_PIm > pim_fid_b+360) && (ph_PIm < pim_fid_a+360)){

	th_vs_p_pim_2[0]->Fill(P_PIm,th_PIm,1.);
		
//th_vs_p
if (th_PIm < (11.09+8./(0.472*(P_PIm+0.3)+0.117))+77.){		

	hist_z_el_1->Fill(z_PIm,1.);
	
cuts_data = true;
};//th_vs_p
};//end of the fiducial cut for first part of sector1
}; //end of the first part of sector1

if ((ph_PIm >= 0) && (ph_PIm <=30)){

	ph_th_pim_all_p[0] -> Fill(th_PIm,ph_PIm,1);
 
	if ((P_PIm > 0.1) && (P_PIm < 1.6)) ph_vs_th_pim[0][int((P_PIm-0.1)/0.1)]->Fill(th_PIm,ph_PIm,1);
		
//fiducial cut sector 1 
if ((th_PIm > th_min)&&(th_PIm < th_max)&&(ph_PIm > pim_fid_b) && (ph_PIm < pim_fid_a)){
	
	th_vs_p_pim_2[0]->Fill(P_PIm,th_PIm,1.);
	
//th_vs_p	
if (th_PIm < (11.09+8./(0.472*(P_PIm+0.3)+0.117))+77.){		

	hist_z_el_1->Fill(z_PIm,1.);
	
cuts_data = true;
};//th_vs_p
};//end of the fiducial cut for second part of sector1
};//end of the second part of sector1
 

  
if ((ph_PIm >= 30) && (ph_PIm <=90)&&(PdHit_PIm!=16)) {
    
	ph_th_pim_all_p[1] -> Fill(th_PIm,ph_PIm-60,1);

	if ((P_PIm > 0.1) && (P_PIm < 1.6)) ph_vs_th_pim[1][int((P_PIm-0.1)/0.1)]->Fill(th_PIm,ph_PIm-60,1);
	
//fiducial cut sector 2 
if ((th_PIm > th_min)&&(th_PIm < th_max)&&(ph_PIm > pim_fid_b+60) && (ph_PIm < pim_fid_a+60)){


	th_vs_p_pim_2[1]->Fill(P_PIm,th_PIm,1.);
	
//th_vs_p
if ((th_PIm < 36.152+3.69909e-05/(5.40783e-06*P_PIm+1.81169e-07)-2.5) || (th_PIm > 36.152+3.69909e-05/(5.40783e-06*P_PIm+1.81169e-07)+2.)){ 

	hist_z_el_2->Fill(z_PIm,1.);
	
cuts_data = true;
};//th_vs_p
};//end of the fiducial cut for sector2
};//end of the sector2

    
if ((ph_PIm >= 90) && (ph_PIm <=150)&&(PdHit_PIm!=44)) {

	ph_th_pim_all_p[2] -> Fill(th_PIm,ph_PIm-120,1);
 

	if ((P_PIm > 0.1) && (P_PIm < 1.6)) ph_vs_th_pim[2][int((P_PIm-0.1)/0.1)]->Fill(th_PIm,ph_PIm-120,1); 

//fiducial cut sector 3 
if ((th_PIm > th_min)&&(th_PIm < th_max)&&(ph_PIm > pim_fid_b+120) && (ph_PIm < pim_fid_a+120)){

	th_vs_p_pim_2[2]->Fill(P_PIm,th_PIm,1.);

//th_vs_p
if (th_PIm <(11.09+8./(0.472*(P_PIm+0.2)+0.117))+85.){ 

	hist_z_el_3->Fill(z_PIm,1.);
	
cuts_data = true;
};//th_vs_p
};//end of the fiducial cut for sector3
};//end of the sector3
  
  
if ((ph_PIm >= 150) && (ph_PIm <=210)){

	ph_th_pim_all_p[3] -> Fill(th_PIm,ph_PIm-180,1);

	if ((P_PIm > 0.1) && (P_PIm < 1.6)) ph_vs_th_pim[3][int((P_PIm-0.1)/0.1)]->Fill(th_PIm,ph_PIm-180,1);	
	
//fiducial cut sector 4
if ((th_PIm > th_min)&&(th_PIm < th_max)&&(ph_PIm > pim_fid_b+180) && (ph_PIm < pim_fid_a+180)){

	th_vs_p_pim_2[3]->Fill(P_PIm,th_PIm,1.);

//th_vs_p
if ((th_PIm < 36.152+3.69909e-05/(5.40783e-06*(P_PIm+0.31)+1.81169e-07)+62.)||((th_PIm>36.152+3.69909e-05/(5.40783e-06*(P_PIm+0.31)+1.81169e-07)+66.)&&(th_PIm <36.152+3.69909e-05/(5.40783e-06*(P_PIm+0.25)+1.81169e-07)+80.))){

	hist_z_el_4->Fill(z_PIm,1.);

cuts_data = true;
};//th_vs_p
};//end of the fiducial cut for sector4
}; //end of the sector4

 
if ((ph_PIm >= 210) && (ph_PIm <=270)&&(PdHit_PIm!=17)) {

	ph_th_pim_all_p[4] -> Fill(th_PIm,ph_PIm-240,1);

	if ((P_PIm > 0.1) && (P_PIm < 1.6)) ph_vs_th_pim[4][int((P_PIm-0.1)/0.1)]->Fill(th_PIm,ph_PIm-240,1);
	
//fiducial cut sector 5 
if ((th_PIm > th_min)&&(th_PIm < th_max)&&(ph_PIm > pim_fid_b+240) && (ph_PIm < pim_fid_a+240)){

	th_vs_p_pim_2[4]->Fill(P_PIm,th_PIm,1.);

//th_vs_p
if ((th_PIm<36.152+3.69909e-05/(5.40783e-06*(P_PIm)+1.81169e-07)-1.)||((th_PIm>36.152+3.69909e-05/(5.40783e-06*(P_PIm-0.01)+1.81169e-07)+3)&&(th_PIm<36.152+3.69909e-05/(5.40783e-06*(P_PIm+0.1)+1.81169e-07)+63.))){ 

	hist_z_el_5->Fill(z_PIm,1.);
	
cuts_data = true;

};//th_vs_p
};//end of the fiducial cut for sector5
}; //end of the sector5

    
if ((ph_PIm >= 270) && (ph_PIm <=330)){
 
	ph_th_pim_all_p[5] -> Fill(th_PIm,ph_PIm-300,1);
 
	if ((P_PIm > 0.1) && (P_PIm < 1.6)) ph_vs_th_pim[5][int((P_PIm-0.1)/0.1)]->Fill(th_PIm,ph_PIm-300,1); 
	
//fiducial cut sector 6 
if ((th_PIm > th_min)&&(th_PIm < th_max)&&(ph_PIm > pim_fid_b+300) && (ph_PIm < pim_fid_a+300)){

	th_vs_p_pim_2[5]->Fill(P_PIm,th_PIm,1.);

//th_vs_p
if ((th_PIm <36.152+3.69909e-05/(5.40783e-06*(P_PIm+0.15)+1.81169e-07)+71.) || (th_PIm > 36.152+3.69909e-05/(5.40783e-06*(P_PIm+0.1)+1.81169e-07)+72.5 )){  

	hist_z_el_6->Fill(z_PIm,1.);

cuts_data = true;

};//th_vs_p
};//end of the fiducial cut for sector6
}; //end of the sector6
  
};//th_min for low momenta
};//vertex cut
};//end of pi- id beta vs p cut 
return cuts_data;
};
  
    
   
   
