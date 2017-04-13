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
#include "beta_func_data.h"
#include "global.h"
#include <iostream>



using namespace std;


void beta_func_data(float& beta_PIp) {
 
 Float_t beta_nom_pip,beta_nom_pim,beta_nom_p,delta_t_pip,delta_t_pim,delta_t_p;
 Float_t p_fid_a_1, p_fid_b_1;
 Float_t pip_fid_a_1,pip_fid_b_1;
 Float_t th_min_1,par1_1,par2_1, pim_fid_a_1,pim_fid_b_1;
 cuts_data particle_ID_data;
 
 
   
beta_nom_pip = P_PIp/sqrt(m_pip*m_pip+P_PIp*P_PIp);
beta_nom_pim = P_PIm/sqrt(m_pip*m_pip+P_PIm*P_PIm);
beta_nom_p = P_P/sqrt(m_proton*m_proton+P_P*P_P);
delta_t_pip = PIp_dist*(1./beta_nom_pip-1/beta_PIp)/30.;
delta_t_pim = PIm_dist*(1./beta_nom_pim-1/beta_PIm)/30.;
delta_t_p = P_dist*(1./beta_nom_p-1/beta_P)/30.; 
 
 
 

pip_fid_a_1 = 24.*(1-exp(-1.*0.08*(th_PIp-9.)));
pip_fid_b_1 = -1.*25.*(1-exp(-1.*0.1*(th_PIp-10.)));
  
if ((ph_PIp >= 330) && (ph_PIp <=360)){
//    if ((ph_PIp > pip_fid_b_1+360) && (ph_PIp < pip_fid_a_1+360)){
    
//time_pip[0][PdHit_PIp-1]->Fill(P_PIp,delta_t_pip,1.); 

 if ((PdHit_PIp==36)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 0.39)*30./PIp_dist);
 
 if ((PdHit_PIp==41)&&(delta_t_pip < -2.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip + 3.8)*30./PIp_dist);
if ((PdHit_PIp==41)&&(delta_t_pip > -2.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip + 0.2)*30./PIp_dist);
 
 if ((PdHit_PIp==42)&&(delta_t_pip < 1.21)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip - 0.23)*30./PIp_dist);
if ((PdHit_PIp==42)&&(delta_t_pip > 1.21)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip - 2.19)*30./PIp_dist);

if ((PdHit_PIp==45)&&(delta_t_pip > -1.71)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip + 0.23)*30./PIp_dist);
if ((PdHit_PIp==45)&&(delta_t_pip < -1.71)&&(delta_t_pip > -3.945)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip + 3.19)*30./PIp_dist);
if ((PdHit_PIp==45)&&(delta_t_pip < -3.945)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip + 4.7)*30./PIp_dist);

if ((PdHit_PIp==45)&&(delta_t_pip > 2.25)&&(P_PIp>0.15)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip - 3.05)*30./PIp_dist);//

if ((PdHit_PIp==46)&&(delta_t_pip > 1.03)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip - 2.04)*30./PIp_dist);
if ((PdHit_PIp==46)&&(delta_t_pip < 1.03)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip - 0.02)*30./PIp_dist);

if ((PdHit_PIp==47)&&(delta_t_pip > -2.815)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 0.17)*30./PIp_dist);
if ((PdHit_PIp==47)&&(delta_t_pip < -2.815)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 5.47)*30./PIp_dist);

if ((PdHit_PIp==48)&&(delta_t_pip > 2.56)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip - 4.32)*30./PIp_dist);//
if ((PdHit_PIp==48)&&(delta_t_pip < 2.56)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip - 0.72)*30./PIp_dist);//



//   };//end of the fiducial cut for first part of sector1
  }; //end of the first part of sector1
  
  if ((ph_PIp >= 0) && (ph_PIp <=30)){
// if ((ph_PIp > pip_fid_b_1) && (ph_PIp < pip_fid_a_1)){ 

//time_pip[0][PdHit_PIp-1]->Fill(P_PIp,delta_t_pip,1.);

if ((PdHit_PIp==36)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 0.39)*30./PIp_dist);


if ((PdHit_PIp==41)&&(delta_t_pip < -2.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip + 3.8)*30./PIp_dist);
if ((PdHit_PIp==41)&&(delta_t_pip > -2.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip + 0.2)*30./PIp_dist);

if ((PdHit_PIp==42)&&(delta_t_pip < 1.21)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 0.23)*30./PIp_dist);
if ((PdHit_PIp==42)&&(delta_t_pip > 1.21)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 2.19)*30./PIp_dist);

if ((PdHit_PIp==45)&&(delta_t_pip > -1.71)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 0.23)*30./PIp_dist);
if ((PdHit_PIp==45)&&(delta_t_pip < -1.71)&&(delta_t_pip > -3.945)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 3.19)*30./PIp_dist);
if ((PdHit_PIp==45)&&(delta_t_pip < -3.945)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 4.7)*30./PIp_dist);

if ((PdHit_PIp==45)&&(delta_t_pip > 2.25)&&(P_PIp>0.15)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip - 3.05)*30./PIp_dist);//

if ((PdHit_PIp==46)&&(delta_t_pip > 1.03)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip - 2.04)*30./PIp_dist);
if ((PdHit_PIp==46)&&(delta_t_pip < 1.03)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip - 0.02)*30./PIp_dist);

if ((PdHit_PIp==47)&&(delta_t_pip > -2.815)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 0.17)*30./PIp_dist);
if ((PdHit_PIp==47)&&(delta_t_pip < -2.815)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 5.47)*30./PIp_dist);

if ((PdHit_PIp==48)&&(delta_t_pip > 2.56)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip - 4.32)*30./PIp_dist);//
if ((PdHit_PIp==48)&&(delta_t_pip < 2.56)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip - 0.72)*30./PIp_dist);//

//   };//end of the fiducial cut for second part of sector1
   }; //end of the second part of sector1
  
  
  if ((ph_PIp >= 30) && (ph_PIp <=90)) {
//  if ((ph_PIp > pip_fid_b_1+60) && (ph_PIp < pip_fid_a_1+60)){

//time_pip[1][PdHit_PIp-1]->Fill(P_PIp,delta_t_pip,1.);

if ((PdHit_PIp==24)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 1.13)*30./PIp_dist);

if ((PdHit_PIp==27)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 0.50)*30./PIp_dist);

if ((PdHit_PIp==29)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 0.41)*30./PIp_dist);

if ((PdHit_PIp==36)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 0.44)*30./PIp_dist);

if ((PdHit_PIp==37)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 0.29)*30./PIp_dist);

if ((PdHit_PIp==40)&&(delta_t_pip < -1.8)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip + 2.9)*30./PIp_dist);//
if ((PdHit_PIp==40)&&(delta_t_pip > 2.2)&&(delta_t_pip < 6.4)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip - 2.75)*30./PIp_dist);//
if ((PdHit_PIp==40)&&(delta_t_pip > 6.4)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip - 7.13)*30./PIp_dist);//


if ((PdHit_PIp==41)&&(delta_t_pip > -0.63)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 1.74)*30./PIp_dist);
if ((PdHit_PIp==41)&&(delta_t_pip < -0.63)&&(delta_t_pip > -3.535)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 3.)*30./PIp_dist);
if ((PdHit_PIp==41)&&(delta_t_pip > -5.7)&&(delta_t_pip < -3.535)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 4.06)*30./PIp_dist);
if ((PdHit_PIp==41)&&(delta_t_pip < -5.7)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 7.34)*30./PIp_dist);


if ((PdHit_PIp==42)&&(delta_t_pip > -2.515)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 0.28)*30./PIp_dist);
if ((PdHit_PIp==42)&&(delta_t_pip < -2.515)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 4.75)*30./PIp_dist);


if ((PdHit_PIp==43)&&(delta_t_pip > -5.875)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 4.82)*30./PIp_dist);
if ((PdHit_PIp==43)&&(delta_t_pip < -5.875)&&(delta_t_pip > -7.57)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 6.93)*30./PIp_dist);
if ((PdHit_PIp==43)&&(delta_t_pip < -7.57)&&(delta_t_pip > -9.51)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 8.21)*30./PIp_dist);
if ((PdHit_PIp==43)&&(delta_t_pip < -9.51)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 10.81)*30./PIp_dist);


if ((PdHit_PIp==44)&&(delta_t_pip > -5.28)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 3.92)*30./PIp_dist);
if ((PdHit_PIp==44)&&(delta_t_pip < -5.28)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 6.64)*30./PIp_dist);

if ((PdHit_PIp==45)&&(delta_t_pip > -2.495)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 0.25)*30./PIp_dist);
if ((PdHit_PIp==45)&&(delta_t_pip < -2.495)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 5.24)*30./PIp_dist);

if ((PdHit_PIp==46)&&(delta_t_pip > -4.43)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 0.82)*30./PIp_dist);
if ((PdHit_PIp==46)&&(delta_t_pip < -4.43)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 8.04)*30./PIp_dist);

if ((PdHit_PIp==47)&&(delta_t_pip < 0.795)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 0.14)*30./PIp_dist);
if ((PdHit_PIp==47)&&(delta_t_pip > 0.795)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 1.73)*30./PIp_dist);

if ((PdHit_PIp==48)&&(delta_t_pip < -2.44)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 4.3)*30./PIp_dist);//
// };//end of the fiducial cut for sector2
};//end of the sector2

  
if ((ph_PIp >= 90) && (ph_PIp <=150)) {
//if ((ph_PIp > pip_fid_b_1+120) && (ph_PIp < pip_fid_a_1+120)){

//time_pip[2][PdHit_PIp-1]->Fill(P_PIp,delta_t_pip,1.);

//if (delta_t_pip < 2.055) beta_new_vs_p_pip[2][PdHit_PIp-1]->Fill(P_PIp,beta_PIp,1.);


if ((PdHit_PIp==25)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 3.34)*30./PIp_dist);

if ((PdHit_PIp==35)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 0.38)*30./PIp_dist);

if ((PdHit_PIp==38)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 0.52)*30./PIp_dist);


if ((PdHit_PIp==40)&&(delta_t_pip < 2.055)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip-0.13)*30./PIp_dist);
//1./(1./beta_nom_pip-(delta_t_pip- 0.13)*30./PIp_dist);

if ((PdHit_PIp==40)&&(delta_t_pip > 2.055)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip-3.98)*30./PIp_dist);

//beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 3.98)*30./PIp_dist);

if ((PdHit_PIp==41)&&(delta_t_pip < 0.445)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 0.03)*30./PIp_dist);
if ((PdHit_PIp==41)&&(delta_t_pip > 0.445)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 0.92)*30./PIp_dist);


if ((PdHit_PIp==42)&&(delta_t_pip < 3.615)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 1.88)*30./PIp_dist);
if ((PdHit_PIp==42)&&(delta_t_pip > 3.615)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 5.35)*30./PIp_dist);

if ((PdHit_PIp==43)&&(delta_t_pip < 0.305)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 0.08)*30./PIp_dist);
if ((PdHit_PIp==43)&&(delta_t_pip > 0.305)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 0.53)*30./PIp_dist);
if ((PdHit_PIp==43)&&(delta_t_pip > 3.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 3.93)*30./PIp_dist);//
if ((PdHit_PIp==43)&&(delta_t_pip < -2.44)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 3.34)*30./PIp_dist);//



if ((PdHit_PIp==44)&&(delta_t_pip < 9.32)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 8.12)*30./PIp_dist);
if ((PdHit_PIp==44)&&(delta_t_pip > 9.32)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 10.52)*30./PIp_dist);

if ((PdHit_PIp==46)&&(delta_t_pip < 2.045)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 0.68)*30./PIp_dist);
if ((PdHit_PIp==46)&&(delta_t_pip > 2.045)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 3.41)*30./PIp_dist);

if ((PdHit_PIp==47)&&(delta_t_pip > -0.985)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 0.03)*30./PIp_dist);
if ((PdHit_PIp==47)&&(delta_t_pip < -0.985)&&(delta_t_pip > -3.245)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 1.94)*30./PIp_dist);
if ((PdHit_PIp==47)&&(delta_t_pip < -3.245)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 4.55)*30./PIp_dist);

if ((PdHit_PIp==48)&&(delta_t_pip < -2.)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 3.2)*30./PIp_dist);//

// };//end of the fiducial cut for sector3
};  //end of the sector3
  
 if ((ph_PIp >= 150) && (ph_PIp <=210)){
// if ((ph_PIp > pip_fid_b_1+180) && (ph_PIp < pip_fid_a_1+180)){

//time_pip[3][PdHit_PIp-1]->Fill(P_PIp,delta_t_pip,1.);

if ((PdHit_PIp==23)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip-0.4)*30./PIp_dist);

if ((PdHit_PIp==25)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip-0.4)*30./PIp_dist);

if ((PdHit_PIp==27)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip-0.47)*30./PIp_dist);

if ((PdHit_PIp==37)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip-0.41)*30./PIp_dist);

if ((PdHit_PIp==38)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip-0.54)*30./PIp_dist);

if ((PdHit_PIp==39)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 4.23)*30./PIp_dist);

if ((PdHit_PIp==40)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 0.21)*30./PIp_dist);

if ((PdHit_PIp==41)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 0.61)*30./PIp_dist);

if ((PdHit_PIp==42)&&(delta_t_pip > -1.46)) {
beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 0.05)*30./PIp_dist);
};
if ((PdHit_PIp==42)&&(delta_t_pip < -1.46)) {
beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 2.97)*30./PIp_dist);
};

if ((PdHit_PIp==43)&&(delta_t_pip > -0.18)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 1.03)*30./PIp_dist);
if ((PdHit_PIp==43)&&(delta_t_pip < -0.18)&&(delta_t_pip > -2.2)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 1.39)*30./PIp_dist);
if ((PdHit_PIp==43)&&(delta_t_pip < -2.2)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 3.01)*30./PIp_dist);


if ((PdHit_PIp==46)&&(delta_t_pip > -0.61)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 0.1)*30./PIp_dist);
if ((PdHit_PIp==46)&&(delta_t_pip < -0.61)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 1.12)*30./PIp_dist);

if ((PdHit_PIp==48)&&(delta_t_pip > -0.1)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 0.57)*30./PIp_dist);//
if ((PdHit_PIp==48)&&(delta_t_pip < -0.1)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 0.93)*30./PIp_dist);//

// };//end of the fiducial cut for sector4
  }; //end of the sector4
  
if ((ph_PIp >= 210) && (ph_PIp <=270)) {
//if ((ph_PIp > pip_fid_b_1+240) && (ph_PIp < pip_fid_a_1+240)){

//time_pip[4][PdHit_PIp-1]->Fill(P_PIp,delta_t_pip,1.);

if ((PdHit_PIp==37)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 0.43)*30./PIp_dist);

if ((PdHit_PIp==40)&&(delta_t_pip > -3.035)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 1.44)*30./PIp_dist);
if ((PdHit_PIp==40)&&(delta_t_pip < -3.035)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 4.63)*30./PIp_dist);

if ((PdHit_PIp==41)&&(delta_t_pip > 0.2421)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 0.48)*30./PIp_dist);
if ((PdHit_PIp==41)&&(delta_t_pip < 0.2421)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 0.0042)*30./PIp_dist);

if ((PdHit_PIp==42)&&(delta_t_pip < -2.31)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 3.96)*30./PIp_dist);
if ((PdHit_PIp==42)&&(delta_t_pip < 0.68)&&(delta_t_pip > -2.31)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 0.66)*30./PIp_dist);
if ((PdHit_PIp==42)&&(delta_t_pip > 0.68)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 2.02)*30./PIp_dist);

if ((PdHit_PIp==43)&&(delta_t_pip < -1.42)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 2.033)*30./PIp_dist);
if ((PdHit_PIp==43)&&(delta_t_pip > -1.42)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 0.51)*30./PIp_dist);

if ((PdHit_PIp==44)&&(delta_t_pip < -0.06)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 0.9)*30./PIp_dist);
if ((PdHit_PIp==44)&&(delta_t_pip > -0.06)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 0.78)*30./PIp_dist);

if ((PdHit_PIp==46)&&(delta_t_pip < 1.145)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 0.08)*30./PIp_dist);
if ((PdHit_PIp==46)&&(delta_t_pip > 1.145)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 2.21)*30./PIp_dist);

if ((PdHit_PIp==47)&&(delta_t_pip < 3.175)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 2.08)*30./PIp_dist);
if ((PdHit_PIp==47)&&(delta_t_pip > 3.175)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 4.27)*30./PIp_dist);

if ((PdHit_PIp==48)&&(delta_t_pip < -2.9)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 4.1)*30./PIp_dist);//
if ((PdHit_PIp==48)&&(delta_t_pip < -4.7)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 5.8)*30./PIp_dist);//
if ((PdHit_PIp==48)&&(delta_t_pip < -6.4)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 7.2)*30./PIp_dist);//


// };//end of the fiducial cut for sector5
}; //end of the sector5
  
  
 if ((ph_PIp >= 270) && (ph_PIp <=330)){
// if ((ph_PIp > pip_fid_b_1+300) && (ph_PIp < pip_fid_a_1+300)){

//time_pip[5][PdHit_PIp-1]->Fill(P_PIp,PIp_dist*(1./beta_nom_pip-1/beta_PIp)/30.,1.);

if ((PdHit_PIp==31)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 1.51)*30./PIp_dist);

if ((PdHit_PIp==25)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 0.38)*30./PIp_dist);

if ((PdHit_PIp==40)&&(delta_t_pip > -1.92)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 0.03)*30./PIp_dist);
if ((PdHit_PIp==40)&&(delta_t_pip < -1.92)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 3.88)*30./PIp_dist);

if ((PdHit_PIp==43)&&(delta_t_pip > 0.955)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 1.85)*30./PIp_dist);
if ((PdHit_PIp==43)&&(delta_t_pip < 0.955)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 0.06)*30./PIp_dist);


if ((PdHit_PIp==44)&&(delta_t_pip > -0.62)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip- 0.01)*30./PIp_dist);
if ((PdHit_PIp==44)&&(delta_t_pip < -0.62)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 1.25)*30./PIp_dist);


if ((PdHit_PIp==45)&&(delta_t_pip > -1.82)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 0.02)*30./PIp_dist);
if ((PdHit_PIp==45)&&(delta_t_pip < -1.82)&&(delta_t_pip > -4.325)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 3.62)*30./PIp_dist);
if ((PdHit_PIp==45)&&(delta_t_pip < -4.325)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 5.03)*30./PIp_dist);

if ((PdHit_PIp==47)&&(delta_t_pip > -1.845)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 0.83)*30./PIp_dist);
if ((PdHit_PIp==47)&&(delta_t_pip < -1.845)&&(delta_t_pip > -4.05)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 2.86)*30./PIp_dist);
if ((PdHit_PIp==47)&&(delta_t_pip < -4.05)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip+ 5.24)*30./PIp_dist);

if ((PdHit_PIp==48)&&(delta_t_pip < 1.25)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip - 0.15)*30./PIp_dist);//
if ((PdHit_PIp==48)&&(delta_t_pip > 1.25)) beta_PIp = 1./(1./beta_nom_pip-(delta_t_pip - 2.75)*30./PIp_dist);//

// };//end of the fiducial cut for sector6
  }; //end of the sector6





th_min_1=(11.09+8./(0.472*P_PIm+0.117));
  par1_1=0.705+1.1*P_PIm;
  par2_1=-63.2-29.3*P_PIm;       
   pim_fid_a_1=30.5*pow((sin((th_PIm-th_min_1)*0.01745)),(par1_1+par2_1/th_PIm+1530./th_PIm/th_PIm))-1;   pim_fid_b_1=-30.5*pow((sin((th_PIm-th_min_1)*0.01745)),(par1_1+par2_1/th_PIm+1530./th_PIm/th_PIm))+1;
   
 if ((ph_PIm >= 330) && (ph_PIm <=360)){
//  if ((ph_PIm > pim_fid_b_1+360) && (ph_PIm < pim_fid_a_1+360)){

// time_pim[0][PdHit_PIm-1]->Fill(P_PIm,PIm_dist/30.*(1./beta_nom_pim-1./beta_PIm),1.); 
  
 if ((PdHit_PIm==36)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 0.39)*30./PIm_dist);
 
 if ((PdHit_PIm==41)&&(delta_t_pim < -2.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim + 3.8)*30./PIm_dist);
if ((PdHit_PIm==41)&&(delta_t_pim > -2.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim + 0.2)*30./PIm_dist);

 if ((PdHit_PIm==42)&&(delta_t_pim < 1.21)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim - 0.23)*30./PIm_dist);
if ((PdHit_PIm==42)&&(delta_t_pim > 1.21)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim - 2.19)*30./PIm_dist);

if ((PdHit_PIm==45)&&(delta_t_pim > -1.71)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim + 0.23)*30./PIm_dist);
if ((PdHit_PIm==45)&&(delta_t_pim < -1.71)&&(delta_t_pim > -3.945)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim + 3.19)*30./PIm_dist);
if ((PdHit_PIm==45)&&(delta_t_pim < -3.945)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim + 4.7)*30./PIm_dist);

if ((PdHit_PIm==46)&&(delta_t_pim > 1.03)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim - 2.04)*30./PIm_dist);
if ((PdHit_PIm==46)&&(delta_t_pim < 1.03)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim - 0.02)*30./PIm_dist);

if ((PdHit_PIm==47)&&(delta_t_pim > -2.815)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim + 0.17)*30./PIm_dist);
if ((PdHit_PIm==47)&&(delta_t_pim < -2.815)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim + 5.47)*30./PIm_dist); 


 //  };//end of the fiducial cut for first part of sector1
  }; //end of the first part of sector1
  
  if ((ph_PIm >= 0) && (ph_PIm <=30)){
//   if ((ph_PIm > pim_fid_b_1) && (ph_PIm < pim_fid_a_1)){
   
 //time_pim[0][PdHit_PIm-1]->Fill(P_PIm,PIm_dist/30.*(1./beta_nom_pim-1./beta_PIm),1.); 
 
 if ((PdHit_PIm==36)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 0.39)*30./PIm_dist);
 
 if ((PdHit_PIm==41)&&(delta_t_pim < -2.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim + 3.8)*30./PIm_dist);
if ((PdHit_PIm==41)&&(delta_t_pim > -2.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim + 0.2)*30./PIm_dist);

 if ((PdHit_PIm==42)&&(delta_t_pim < 1.21)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim - 0.23)*30./PIm_dist);
if ((PdHit_PIm==42)&&(delta_t_pim > 1.21)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim - 2.19)*30./PIm_dist);

if ((PdHit_PIm==45)&&(delta_t_pim > -1.71)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim + 0.23)*30./PIm_dist);
if ((PdHit_PIm==45)&&(delta_t_pim < -1.71)&&(delta_t_pim > -3.945)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim + 3.19)*30./PIm_dist);
if ((PdHit_PIm==45)&&(delta_t_pim < -3.945)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim + 4.7)*30./PIm_dist);

if ((PdHit_PIm==46)&&(delta_t_pim > 1.03)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim - 2.04)*30./PIm_dist);
if ((PdHit_PIm==46)&&(delta_t_pim < 1.03)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim - 0.02)*30./PIm_dist);

if ((PdHit_PIm==47)&&(delta_t_pim > -2.815)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim + 0.17)*30./PIm_dist);
if ((PdHit_PIm==47)&&(delta_t_pim < -2.815)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim + 5.47)*30./PIm_dist); 
 
  
//    };//end of the fiducial cut for second part of sector1
   };//end of the second part of sector1
  
    if ((ph_PIm >= 30) && (ph_PIm <=90)) {
//if ((ph_PIm > pim_fid_b_1+60) && (ph_PIm < pim_fid_a_1+60)){

//time_pim[1][PdHit_PIm-1]->Fill(P_PIm,PIm_dist/30.*(1./beta_nom_pim-1./beta_PIm),1.);

if ((PdHit_PIm==24)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 1.13)*30./PIm_dist);

if ((PdHit_PIm==27)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 0.50)*30./PIm_dist);

if ((PdHit_PIm==29)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 0.41)*30./PIm_dist);

if ((PdHit_PIm==36)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 0.44)*30./PIm_dist);

if ((PdHit_PIm==37)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 0.29)*30./PIm_dist);

if ((PdHit_PIm==40)&&(delta_t_pim < -1.8)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim + 2.9)*30./PIm_dist);//
if ((PdHit_PIm==40)&&(delta_t_pim > 2.2)&&(delta_t_pim < 6.4)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim - 2.75)*30./PIm_dist);//
if ((PdHit_PIm==40)&&(delta_t_pim > 6.4)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim - 7.13)*30./PIm_dist);//

if ((PdHit_PIm==41)&&(delta_t_pim > -0.63)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 1.74)*30./PIm_dist);
if ((PdHit_PIm==41)&&(delta_t_pim < -0.63)&&(delta_t_pim > -3.535)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 3.)*30./PIm_dist);
if ((PdHit_PIm==41)&&(delta_t_pim > -5.7)&&(delta_t_pim < -3.535)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 4.06)*30./PIm_dist);
if ((PdHit_PIm==41)&&(delta_t_pim < -5.7)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 7.34)*30./PIm_dist);


if ((PdHit_PIm==42)&&(delta_t_pim > -2.515)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 0.28)*30./PIm_dist);
if ((PdHit_PIm==42)&&(delta_t_pim < -2.515)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 4.75)*30./PIm_dist);


if ((PdHit_PIm==43)&&(delta_t_pim > -5.875)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 4.82)*30./PIm_dist);
if ((PdHit_PIm==43)&&(delta_t_pim < -5.875)&&(delta_t_pim > -7.57)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 6.93)*30./PIm_dist);
if ((PdHit_PIm==43)&&(delta_t_pim < -7.57)&&(delta_t_pim > -9.51)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 8.21)*30./PIm_dist);
if ((PdHit_PIm==43)&&(delta_t_pim < -9.51)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 10.81)*30./PIm_dist);

if ((PdHit_PIm==44)&&(delta_t_pim > -5.28)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 3.92)*30./PIm_dist);
if ((PdHit_PIm==44)&&(delta_t_pim < -5.28)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 6.64)*30./PIm_dist);

if ((PdHit_PIm==45)&&(delta_t_pim > -2.495)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 0.25)*30./PIm_dist);
if ((PdHit_PIm==45)&&(delta_t_pim < -2.495)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 5.24)*30./PIm_dist);

if ((PdHit_PIm==46)&&(delta_t_pim > -4.43)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 0.82)*30./PIm_dist);
if ((PdHit_PIm==46)&&(delta_t_pim < -4.43)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 8.04)*30./PIm_dist);

if ((PdHit_PIm==47)&&(delta_t_pim < 0.795)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 0.14)*30./PIm_dist);
if ((PdHit_PIm==47)&&(delta_t_pim > 0.795)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 1.73)*30./PIm_dist);
  
// };//end of the fiducial cut for sector2
};//end of the sector2
    
if ((ph_PIm >= 90) && (ph_PIm <=150)) {
// if ((ph_PIm > pim_fid_b_1+120) && (ph_PIm < pim_fid_a_1+120)){

//time_pim[2][PdHit_PIm-1]->Fill(P_PIm,PIm_dist/30.*(1./beta_nom_pim-1./beta_PIm),1.); 

 if ((PdHit_PIm==25)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 3.34)*30./PIm_dist);

if ((PdHit_PIm==35)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 0.38)*30./PIm_dist);

if ((PdHit_PIm==38)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 0.52)*30./PIm_dist);


if ((PdHit_PIm==40)&&(delta_t_pim < 2.055)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim-0.13)*30./PIm_dist);
//1./(1./beta_nom_pim-(delta_t_pim- 0.13)*30./PIm_dist);

if ((PdHit_PIm==40)&&(delta_t_pim > 2.055)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim-3.98)*30./PIm_dist);

//beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 3.98)*30./PIm_dist);

if ((PdHit_PIm==41)&&(delta_t_pim < 0.445)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 0.03)*30./PIm_dist);
if ((PdHit_PIm==41)&&(delta_t_pim > 0.445)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 0.92)*30./PIm_dist);


if ((PdHit_PIm==42)&&(delta_t_pim < 3.615)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 1.88)*30./PIm_dist);
if ((PdHit_PIm==42)&&(delta_t_pim > 3.615)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 5.35)*30./PIm_dist);

if ((PdHit_PIm==43)&&(delta_t_pim < 0.305)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 0.08)*30./PIm_dist);
if ((PdHit_PIm==43)&&(delta_t_pim > 0.305)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 0.53)*30./PIm_dist);
if ((PdHit_PIm==43)&&(delta_t_pim > 3.)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 3.93)*30./PIm_dist);//
if ((PdHit_PIm==43)&&(delta_t_pim < -2.44)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 3.34)*30./PIm_dist);//



if ((PdHit_PIm==44)&&(delta_t_pim < 9.32)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 8.12)*30./PIm_dist);
if ((PdHit_PIm==44)&&(delta_t_pim > 9.32)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 10.52)*30./PIm_dist);


if ((PdHit_PIm==46)&&(delta_t_pim < 2.045)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 0.68)*30./PIm_dist);
if ((PdHit_PIm==46)&&(delta_t_pim > 2.045)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 3.41)*30./PIm_dist);

if ((PdHit_PIm==47)&&(delta_t_pim > -0.985)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 0.03)*30./PIm_dist);
if ((PdHit_PIm==47)&&(delta_t_pim < -0.985)&&(delta_t_pim > -3.245)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 1.94)*30./PIm_dist);
if ((PdHit_PIm==47)&&(delta_t_pim < -3.245)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 4.55)*30./PIm_dist);


// };//end of the fiducial cut for sector3
};  //end of the sector3
  
 if ((ph_PIm >= 150) && (ph_PIm <=210)){
//  if ((ph_PIm > pim_fid_b_1+180) && (ph_PIm < pim_fid_a_1+180)){

//time_pim[3][PdHit_PIm-1]->Fill(P_PIm,PIm_dist/30.*(1./beta_nom_pim-1./beta_PIm),1.);

if ((PdHit_PIm==23)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim-0.4)*30./PIm_dist);

if ((PdHit_PIm==25)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim-0.4)*30./PIm_dist);

if ((PdHit_PIm==27)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim-0.47)*30./PIm_dist);

if ((PdHit_PIm==37)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim-0.41)*30./PIm_dist);

if ((PdHit_PIm==38)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim-0.54)*30./PIm_dist);

if ((PdHit_PIm==39)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 4.23)*30./PIm_dist);

if ((PdHit_PIm==40)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 0.21)*30./PIm_dist);

if ((PdHit_PIm==41)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 0.61)*30./PIm_dist);

if ((PdHit_PIm==42)&&(delta_t_pim > -1.46)) {
beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 0.05)*30./PIm_dist);
};
if ((PdHit_PIm==42)&&(delta_t_pim < -1.46)) {
beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 2.97)*30./PIm_dist);
};

if ((PdHit_PIm==43)&&(delta_t_pim > -0.18)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 1.03)*30./PIm_dist);
if ((PdHit_PIm==43)&&(delta_t_pim < -0.18)&&(delta_t_pim > -2.2)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 1.39)*30./PIm_dist);
if ((PdHit_PIm==43)&&(delta_t_pim < -2.2)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 3.01)*30./PIm_dist);


if ((PdHit_PIm==46)&&(delta_t_pim > -0.61)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 0.1)*30./PIm_dist);
if ((PdHit_PIm==46)&&(delta_t_pim < -0.61)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 1.12)*30./PIm_dist);



// };//end of the fiducial cut for sector4
 }; //end of the sector4
  
if ((ph_PIm >= 210) && (ph_PIm <=270)) {
//if ((ph_PIm > pim_fid_b_1+240) && (ph_PIm < pim_fid_a_1+240)){

//time_pim[4][PdHit_PIm-1]->Fill(P_PIm,PIm_dist/30.*(1./beta_nom_pim-1./beta_PIm),1.);

if ((PdHit_PIm==37)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 0.43)*30./PIm_dist);

if ((PdHit_PIm==40)&&(delta_t_pim > -3.035)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 1.44)*30./PIm_dist);
if ((PdHit_PIm==40)&&(delta_t_pim < -3.035)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 4.63)*30./PIm_dist);

if ((PdHit_PIm==41)&&(delta_t_pim > 0.2421)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 0.48)*30./PIm_dist);
if ((PdHit_PIm==41)&&(delta_t_pim < 0.2421)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 0.0042)*30./PIm_dist);

if ((PdHit_PIm==42)&&(delta_t_pim < -2.31)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 3.96)*30./PIm_dist);
if ((PdHit_PIm==42)&&(delta_t_pim < 0.68)&&(delta_t_pim > -2.31)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 0.66)*30./PIm_dist);
if ((PdHit_PIm==42)&&(delta_t_pim > 0.68)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 2.02)*30./PIm_dist);

if ((PdHit_PIm==43)&&(delta_t_pim < -1.42)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 2.033)*30./PIm_dist);
if ((PdHit_PIm==43)&&(delta_t_pim > -1.42)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 0.51)*30./PIm_dist);

if ((PdHit_PIm==44)&&(delta_t_pim < -0.06)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 0.9)*30./PIm_dist);
if ((PdHit_PIm==44)&&(delta_t_pim > -0.06)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 0.78)*30./PIm_dist);

if ((PdHit_PIm==46)&&(delta_t_pim < 1.145)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 0.08)*30./PIm_dist);
if ((PdHit_PIm==46)&&(delta_t_pim > 1.145)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 2.21)*30./PIm_dist);

if ((PdHit_PIm==47)&&(delta_t_pim < 3.175)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 2.08)*30./PIm_dist);
if ((PdHit_PIm==47)&&(delta_t_pim > 3.175)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 4.27)*30./PIm_dist);



// };//end of the fiducial cut for sector5
}; //end of the sector5
    
 if ((ph_PIm >= 270) && (ph_PIm <=330)){
//  if ((ph_PIm > pim_fid_b_1+300) && (ph_PIm < pim_fid_a_1+300)){

//time_pim[5][PdHit_PIm-1]->Fill(P_PIm,PIm_dist/30.*(1./beta_nom_pim-1./beta_PIm),1.);
  
if ((PdHit_PIm==31)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 1.51)*30./PIm_dist);

if ((PdHit_PIm==25)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 0.38)*30./PIm_dist);

if ((PdHit_PIm==40)&&(delta_t_pim > -1.92)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 0.03)*30./PIm_dist);
if ((PdHit_PIm==40)&&(delta_t_pim < -1.92)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 3.88)*30./PIm_dist);

if ((PdHit_PIm==43)&&(delta_t_pim > 0.955)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 1.85)*30./PIm_dist);
if ((PdHit_PIm==43)&&(delta_t_pim < 0.955)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 0.06)*30./PIm_dist);


if ((PdHit_PIm==44)&&(delta_t_pim > -0.62)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim- 0.01)*30./PIm_dist);
if ((PdHit_PIm==44)&&(delta_t_pim < -0.62)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 1.25)*30./PIm_dist);

if ((PdHit_PIm==45)&&(delta_t_pim > -1.82)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 0.02)*30./PIm_dist);
if ((PdHit_PIm==45)&&(delta_t_pim < -1.82)&&(delta_t_pim > -4.325)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 3.62)*30./PIm_dist);
if ((PdHit_PIm==45)&&(delta_t_pim < -4.325)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 5.03)*30./PIm_dist);

if ((PdHit_PIm==47)&&(delta_t_pim > -1.845)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 0.83)*30./PIm_dist);
if ((PdHit_PIm==47)&&(delta_t_pim < -1.845)&&(delta_t_pim > -4.05)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 2.86)*30./PIm_dist);
if ((PdHit_PIm==47)&&(delta_t_pim < -4.05)) beta_PIm = 1./(1./beta_nom_pim-(delta_t_pim+ 5.24)*30./PIm_dist); 


// };//end of the fiducial cut for sector6
 }; //end of the sector6
  
 
   



p_fid_a_1 = 24.*(1-exp(-1.*0.08*(th_P-9.)));
   p_fid_b_1 = -1.*25.*(1-exp(-1.*0.1*(th_P-10.)));
if ((ph_P >= 330)&& (ph_P <= 360) ) {
//  if ((ph_P > p_fid_b_1+360) && (ph_P < p_fid_a_1+360)){
  
// time_p[0][PdHit_P-1]->Fill(P_P,P_dist/30.*(1./beta_nom_p-1./beta_P),1.);  
 
// if ((PdHit_P==36)&&(delta_t_p > 0.)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 0.39)*30./P_dist);
 
 //if ((PdHit_P==41)&&(delta_t_p < -2.)) beta_P = 1./(1./beta_nom_p-(delta_t_p + 3.8)*30./P_dist);
//if ((PdHit_P==41)&&(delta_t_p > -2.)) beta_P = 1./(1./beta_nom_p-(delta_t_p + 0.2)*30./P_dist);
 
 if ((PdHit_P==42)&&(delta_t_p < 1.21)) beta_P = 1./(1./beta_nom_p-(delta_t_p - 0.23)*30./P_dist);
if ((PdHit_P==42)&&(delta_t_p > 1.21)) beta_P = 1./(1./beta_nom_p-(delta_t_p - 2.19)*30./P_dist);

if ((PdHit_P==45)&&(delta_t_p > -1.71)) beta_P = 1./(1./beta_nom_p-(delta_t_p + 0.23)*30./P_dist);
if ((PdHit_P==45)&&(delta_t_p < -1.71)&&(delta_t_p > -3.945)) beta_P = 1./(1./beta_nom_p-(delta_t_p + 3.19)*30./P_dist);
if ((PdHit_P==45)&&(delta_t_p < -3.945)) beta_P = 1./(1./beta_nom_p-(delta_t_p + 4.7)*30./P_dist);

//if ((PdHit_P==46)&&(delta_t_p > 1.03)) beta_P = 1./(1./beta_nom_p-(delta_t_p - 2.04)*30./P_dist);
//if ((PdHit_P==46)&&(delta_t_p < 1.03)) beta_P = 1./(1./beta_nom_p-(delta_t_p - 0.02)*30./P_dist);

if ((PdHit_P==47)&&(delta_t_p > -2.815)) beta_P = 1./(1./beta_nom_p-(delta_t_p + 0.17)*30./P_dist);
if ((PdHit_P==47)&&(delta_t_p < -2.815)) beta_P = 1./(1./beta_nom_p-(delta_t_p + 5.47)*30./P_dist); 



//};//end of fiducial cut for the first part of the first sector
 };//end of the first part of the first sector
 
 if ((ph_P >= 0) && (ph_P <= 30)) {
//if ((ph_P > p_fid_b_1) && (ph_P < p_fid_a_1)){

 //time_p[0][PdHit_P-1]->Fill(P_P,P_dist/30.*(1./beta_nom_p-1./beta_P),1.);
 
// if ((PdHit_P==36)&&(delta_t_p > 0.)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 0.39)*30./P_dist);
 
// if ((PdHit_P==41)&&(delta_t_p < -2.)) beta_P = 1./(1./beta_nom_p-(delta_t_p + 3.8)*30./P_dist);
//if ((PdHit_P==41)&&(delta_t_p > -2.)) beta_P = 1./(1./beta_nom_p-(delta_t_p + 0.2)*30./P_dist);

 if ((PdHit_P==42)&&(delta_t_p < 1.21)) beta_P = 1./(1./beta_nom_p-(delta_t_p - 0.23)*30./P_dist);
if ((PdHit_P==42)&&(delta_t_p > 1.21)) beta_P = 1./(1./beta_nom_p-(delta_t_p - 2.19)*30./P_dist);

if ((PdHit_P==45)&&(delta_t_p > -1.71)) beta_P = 1./(1./beta_nom_p-(delta_t_p + 0.23)*30./P_dist);
if ((PdHit_P==45)&&(delta_t_p < -1.71)&&(delta_t_p > -3.945)) beta_P = 1./(1./beta_nom_p-(delta_t_p + 3.19)*30./P_dist);
if ((PdHit_P==45)&&(delta_t_p < -3.945)) beta_P = 1./(1./beta_nom_p-(delta_t_p + 4.7)*30./P_dist);

//if ((PdHit_P==46)&&(delta_t_p > 1.03)) beta_P = 1./(1./beta_nom_p-(delta_t_p - 2.04)*30./P_dist);
//if ((PdHit_P==46)&&(delta_t_p < 1.03)) beta_P = 1./(1./beta_nom_p-(delta_t_p - 0.02)*30./P_dist);

if ((PdHit_P==47)&&(delta_t_p > -2.815)) beta_P = 1./(1./beta_nom_p-(delta_t_p + 0.17)*30./P_dist);
if ((PdHit_P==47)&&(delta_t_p < -2.815)) beta_P = 1./(1./beta_nom_p-(delta_t_p + 5.47)*30./P_dist); 



//};//end of fiducial cut for the second part of the first sector
 };//end of the second part of the first sector
 
 
 
 if ((ph_P >= 30) && (ph_P <=90)){ 
// if ((ph_P > p_fid_b_1+60) && (ph_P < p_fid_a_1+60)){
  
//time_p[1][PdHit_P-1]->Fill(P_P,P_dist/30.*(1./beta_nom_p-1./beta_P),1.);
 
if ((PdHit_P==24)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 1.13)*30./P_dist);

//if ((PdHit_P==27)&&(delta_t_p > 0.)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 0.50)*30./P_dist);

//if ((PdHit_P==29)&&(delta_t_p > 0.)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 0.41)*30./P_dist);

//if ((PdHit_P==36)&&(delta_t_p > 0.)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 0.44)*30./P_dist);

//if ((PdHit_P==37)&&(delta_t_p > 0.)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 0.29)*30./P_dist);

if ((PdHit_P==40)&&(delta_t_p < -1.8)) beta_P = 1./(1./beta_nom_p-(delta_t_p + 2.9)*30./P_dist);//
if ((PdHit_P==40)&&(delta_t_p > 2.2)&&(delta_t_p < 6.4)) beta_P = 1./(1./beta_nom_p-(delta_t_p - 2.75)*30./P_dist);//
if ((PdHit_P==40)&&(delta_t_p > 6.4)) beta_P = 1./(1./beta_nom_p-(delta_t_p - 7.13)*30./P_dist);//


if ((PdHit_P==41)&&(delta_t_p > -0.63)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 1.74)*30./P_dist);
if ((PdHit_P==41)&&(delta_t_p < -0.63)&&(delta_t_p > -3.535)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 3.)*30./P_dist);
if ((PdHit_P==41)&&(delta_t_p > -5.7)&&(delta_t_p < -3.535)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 4.06)*30./P_dist);
if ((PdHit_P==41)&&(delta_t_p < -5.7)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 7.34)*30./P_dist);

if ((PdHit_P==42)&&(delta_t_p > -2.515)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 0.28)*30./P_dist);
if ((PdHit_P==42)&&(delta_t_p < -2.515)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 4.75)*30./P_dist);


if ((PdHit_P==43)&&(delta_t_p > -5.875)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 4.82)*30./P_dist);
if ((PdHit_P==43)&&(delta_t_p < -5.875)&&(delta_t_p > -7.57)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 6.93)*30./P_dist);
if ((PdHit_P==43)&&(delta_t_p < -7.57)&&(delta_t_p > -9.51)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 8.21)*30./P_dist);
if ((PdHit_P==43)&&(delta_t_p < -9.51)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 10.81)*30./P_dist);

if ((PdHit_P==44)&&(delta_t_p > -5.28)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 3.92)*30./P_dist);
if ((PdHit_P==44)&&(delta_t_p < -5.28)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 6.64)*30./P_dist);

if ((PdHit_P==45)&&(delta_t_p > -2.495)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 0.25)*30./P_dist);
if ((PdHit_P==45)&&(delta_t_p < -2.495)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 5.24)*30./P_dist);

if ((PdHit_P==46)&&(delta_t_p > -4.43)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 0.82)*30./P_dist);
if ((PdHit_P==46)&&(delta_t_p < -4.43)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 8.04)*30./P_dist);


//};//end of the fiducial cut for sector2
 };//end of the sector2
 
 if ((ph_P >=90) && (ph_P <=150)) {
//if ((ph_P > p_fid_b_1+120) && (ph_P < p_fid_a_1+120)){

//time_p[2][PdHit_P-1]->Fill(P_P,P_dist/30.*(1./beta_nom_p-1./beta_P),1.);

 if ((PdHit_P==25)&&(delta_t_p < 0.)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 3.34)*30./P_dist);

//if ((PdHit_P==35)&&(delta_t_p > 0.)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 0.38)*30./P_dist);

//if ((PdHit_P==38)&&(delta_t_p > 0.)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 0.52)*30./P_dist);


if ((PdHit_P==40)&&(delta_t_p < 2.055)) beta_P = 1./(1./beta_nom_p-(delta_t_p-0.13)*30./P_dist);
//1./(1./beta_nom_p-(delta_t_p- 0.13)*30./P_dist);

if ((PdHit_P==40)&&(delta_t_p > 2.055)) beta_P = 1./(1./beta_nom_p-(delta_t_p-3.98)*30./P_dist);

//beta_P = 1./(1./beta_nom_p-(delta_t_p- 3.98)*30./P_dist);

//if ((PdHit_P==41)&&(delta_t_p < 0.445)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 0.03)*30./P_dist);
//if ((PdHit_P==41)&&(delta_t_p > 0.445)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 0.92)*30./P_dist);


if ((PdHit_P==42)&&(delta_t_p < 3.615)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 1.88)*30./P_dist);
if ((PdHit_P==42)&&(delta_t_p > 3.615)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 5.35)*30./P_dist);

//if ((PdHit_P==43)&&(delta_t_p < 0.305)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 0.08)*30./P_dist);
if ((PdHit_P==43)&&(delta_t_p > 2.3)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 3.2)*30./P_dist);

if ((PdHit_P==44)&&(delta_t_p < 9.32)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 8.12)*30./P_dist);
if ((PdHit_P==44)&&(delta_t_p > 9.32)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 10.52)*30./P_dist);


if ((PdHit_P==46)&&(delta_t_p < 2.045)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 0.68)*30./P_dist);
if ((PdHit_P==46)&&(delta_t_p > 2.045)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 3.41)*30./P_dist);

//if ((PdHit_P==47)&&(delta_t_p > -0.985)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 0.03)*30./P_dist);
//if ((PdHit_P==47)&&(delta_t_p < -0.985)&&(delta_t_p > -3.245)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 1.94)*30./P_dist);
//if ((PdHit_P==47)&&(delta_t_p < -3.245)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 4.55)*30./P_dist);


//};//end of the fiducial cut for sector3
 };//end of the sector3
 
 if ((ph_P >= 150) && (ph_P <= 210)) {
//if ((ph_P > p_fid_b_1+180) && (ph_P < p_fid_a_1+180)){

//time_p[3][PdHit_P-1]->Fill(P_P,P_dist/30.*(1./beta_nom_p-1./beta_P),1.);

//if ((PdHit_P==23)&&(delta_t_p >0.)) beta_P = 1./(1./beta_nom_p-(delta_t_p-0.4)*30./P_dist);

//if ((PdHit_P==25)&&(delta_t_p >0.)) beta_P = 1./(1./beta_nom_p-(delta_t_p-0.4)*30./P_dist);

//if ((PdHit_P==27)&&(delta_t_p >0.)) beta_P = 1./(1./beta_nom_p-(delta_t_p-0.47)*30./P_dist);

//if ((PdHit_P==37)&&(delta_t_p >0.)) beta_P = 1./(1./beta_nom_p-(delta_t_p-0.41)*30./P_dist);

//if ((PdHit_P==38)&&(delta_t_p >0.)) beta_P = 1./(1./beta_nom_p-(delta_t_p-0.54)*30./P_dist);

if ((PdHit_P==39)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 4.23)*30./P_dist);

//if ((PdHit_P==40)&&(delta_t_p >0.)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 0.21)*30./P_dist);
//
//if ((PdHit_P==41)&&(delta_t_p <0.)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 0.61)*30./P_dist);

if ((PdHit_P==42)&&(delta_t_p > -1.46)) {
beta_P = 1./(1./beta_nom_p-(delta_t_p- 0.05)*30./P_dist);
};
if ((PdHit_P==42)&&(delta_t_p < -1.46)) {
beta_P = 1./(1./beta_nom_p-(delta_t_p+ 2.97)*30./P_dist);
};

if ((PdHit_P==43)&&(delta_t_p > -0.18)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 1.03)*30./P_dist);
if ((PdHit_P==43)&&(delta_t_p < -0.18)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 0.8)*30./P_dist);
//if ((PdHit_P==43)&&(delta_t_p < -2.2)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 3.01)*30./P_dist);


//if ((PdHit_P==46)&&(delta_t_p > -0.61)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 0.1)*30./P_dist);
//if ((PdHit_P==46)&&(delta_t_p < -0.61)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 1.12)*30./P_dist);



//};//end of the fiducial cut for sector4
};//end of the sector4

if ((ph_P >= 210) && (ph_P <=270)) {
// if ((ph_P > p_fid_b_1+240) && (ph_P < p_fid_a_1+240)){

 //time_p[4][PdHit_P-1]->Fill(P_P,P_dist/30.*(1./beta_nom_p-1./beta_P),1.);
 
//if ((PdHit_P==37)&&(delta_t_p < 0.)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 0.43)*30./P_dist);
 
if ((PdHit_P==40)&&(delta_t_p > -3.035)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 1.44)*30./P_dist);
if ((PdHit_P==40)&&(delta_t_p < -3.035)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 4.63)*30./P_dist);

//if ((PdHit_P==41)&&(delta_t_p > 0.2421)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 0.48)*30./P_dist);
//if ((PdHit_P==41)&&(delta_t_p < 0.2421)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 0.0042)*30./P_dist);

if ((PdHit_P==42)&&(delta_t_p < -2.31)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 3.96)*30./P_dist);
if ((PdHit_P==42)&&(delta_t_p < 0.68)&&(delta_t_p > -2.31)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 0.66)*30./P_dist);
if ((PdHit_P==42)&&(delta_t_p > 0.68)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 2.02)*30./P_dist);

if ((PdHit_P==43)&&(delta_t_p < -1.42)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 2.033)*30./P_dist);
if ((PdHit_P==43)&&(delta_t_p > -1.42)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 0.51)*30./P_dist);

//if ((PdHit_P==44)&&(delta_t_p < -0.06)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 0.9)*30./P_dist);
//if ((PdHit_P==44)&&(delta_t_p > -0.06)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 0.78)*30./P_dist);

if ((PdHit_P==46)&&(delta_t_p < 1.145)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 0.08)*30./P_dist);
if ((PdHit_P==46)&&(delta_t_p > 1.145)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 2.21)*30./P_dist);

//if ((PdHit_P==47)&&(delta_t_p < 3.175)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 2.08)*30./P_dist);
//if ((PdHit_P==47)&&(delta_t_p > 3.175)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 4.27)*30./P_dist);


 };//end of the sector5
// };
 if ((ph_P >= 270) && (ph_P <=330)) {
// if ((ph_P > p_fid_b_1+300) && (ph_P < p_fid_a_1+300)){

//time_p[5][PdHit_P-1]->Fill(P_P,P_dist/30.*(1./beta_nom_p-1./beta_P),1.);
   
if ((PdHit_P==31)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 1.51)*30./P_dist);

//if ((PdHit_P==25)&&(delta_t_p > 0.)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 0.38)*30./P_dist);

if ((PdHit_P==40)&&(delta_t_p > -1.92)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 0.03)*30./P_dist);
if ((PdHit_P==40)&&(delta_t_p < -1.92)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 3.88)*30./P_dist);

////if ((PdHit_P==43)&&(delta_t_p > 0.955)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 1.85)*30./P_dist);
//if ((PdHit_P==43)&&(delta_t_p < 0.955)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 0.06)*30./P_dist);

//if ((PdHit_P==44)&&(delta_t_p > -0.62)) beta_P = 1./(1./beta_nom_p-(delta_t_p- 0.01)*30./P_dist);
//if ((PdHit_P==44)&&(delta_t_p < -0.62)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 1.25)*30./P_dist);

if ((PdHit_P==45)&&(delta_t_p > -1.82)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 0.02)*30./P_dist);
if ((PdHit_P==45)&&(delta_t_p < -1.82)&&(delta_t_p > -4.325)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 3.62)*30./P_dist);
if ((PdHit_P==45)&&(delta_t_p < -4.325)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 5.03)*30./P_dist);

if ((PdHit_P==47)&&(delta_t_p > -1.845)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 0.83)*30./P_dist);
if ((PdHit_P==47)&&(delta_t_p < -1.845)&&(delta_t_p > -4.05)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 2.86)*30./P_dist);
if ((PdHit_P==47)&&(delta_t_p < -4.05)) beta_P = 1./(1./beta_nom_p-(delta_t_p+ 5.24)*30./P_dist); 

// };//end of the fiducial cut for sector6
 };//end of the sector6
  
 
 };
 
 
 //cout << "qqqq";
 
 
 
 
 
 
 
 
 
 
   
   
   
   
   
   
   
   
