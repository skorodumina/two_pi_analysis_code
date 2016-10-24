#include <TMath.h>
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include "cuts_empty.h"
#include "global.h"
#include <iostream>




using namespace std;

//ostringstream qqq1;



 bool cuts_empty::Electron_cuts_empty(){
   
   
   

   
   bool cuts_empty;
   Float_t th_min,par1,par2,fid_a,fid_b,a,b;
   Short_t i;
   Float_t ph_el_arr[3][6][18] = {{{1000.,1000.,22.,24.,25.,28.,35.,30.,30.,34.,26.,30.,32.,27.,36.,31.,28.,1000.},
                             {1000.,1000.,25.,16.,17.,22.,23.,21.,22.,24.,21.,31.,31.,31.,32.,29.,29.,1000.},
                             {1000.,1000.,20.,23.,24.,23.,19.,23.,31.,28.,27.,28.,32.,34.,35.,37.,37.,1000.},
                             {1000.,1000.,24.,20.,24.,19.,28.,26.,24.,27.,25.,29.,22.,23.,34.,37.,27.,1000.},
                             {1000.,1000.,31.,27.,25.,28.,29.,23.,22.,28.,35.,29.,51.,29.,39.,44.,36.,1000.},
                             {1000.,1000.,28.,28.,22.,22.,30.,30.,28.,46.,27.,39.,38.,29.,41.,42.,36.,1000.}},
			     {{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                             {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
		             {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
	                     {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
		             {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
			     {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}},
			     
			     {{1000.,1000.,18.,20.,24.,22.,22.,27.,26.,28.,32.,26.,23.,30.,28.,29.,25.,1000.},
                             {1000.,1000.,31.,20.,22.,33.,27.,26.,35.,35.,23.,24.,30.,45.,40.,40.,46.,1000.},
			     {1000.,1000.,23.,26.,28.,23.,23.,28.,24.,29.,25.,27.,29.,35.,37.,29.,36.,1000.},
			     {1000.,1000.,20.,21.,22.,24.,25.,34.,26.,27.,24.,27.,20.,28.,32.,27.,34.,1000.},
			     {1000.,1000.,27.,26.,27.,51.,40.,35.,31.,33.,33.,46.,40.,37.,48.,44.,31.,1000.},
			     {1000.,1000.,31.,21.,28.,33.,31.,28.,32.,38.,38.,36.,29.,35.,90.,35.,36.,1000.}}};


 cuts_empty = false; 
   
   
  
     if (P_EL > 0.461) {
     
  th_min=(9.5+17./(P_EL+0.2));
  par1=0.85+1.1*P_EL;
  par2=-62.8-30.*P_EL;       
   fid_a=37.3*pow((sin((th_EL-th_min)*0.01745)),(par1+par2/th_EL+1525./th_EL/th_EL));
     
   fid_b=-37.3*pow((sin((th_EL-th_min)*0.01745)),(par1+par2/th_EL+1525./th_EL/th_EL));  
   a = fid_a;
   b = fid_b; 
   
//   P_EL = floor((P_EL*1000) + 0.5)/1000;


       switch (sector) {
case 1 : 


 
 if (((ECT/P_EL) > (-0.0664988*P_EL*P_EL+0.208672*P_EL+0.0450621)) && ((ECT/P_EL) < (0.0235467*P_EL*P_EL-0.0656095*P_EL)+0.400512)) {
 
// if (nphe > 25.) {


if ((ph_EL >= 330) && (ph_EL <= 360)){


if ((th_EL > th_min) && (th_EL < 50) && (ph_EL > fid_b+360) && (ph_EL < fid_a+360)){

// hist_z_el_1_empty->Fill(z_EL,1.);
 //~fid cuts in Cherenkov plane
 if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {

//cut on average number of photoelectrons
//if (avrg_nphe_sector1->GetBinContent(int((theta_cc+5.)*200./60.),int((ph_cc+30.)*200./60.)) > 70.) {
//nphe_sector1_after->Fill(nphe,1.);
if ((theta_cc >1.64786+2.16008*segment-0.00433375*segment*segment)&&(theta_cc <11.3593+1.38551*segment+0.0451710*segment*segment)){

if (nphe > ph_el_arr[pmt_hit+1][0][segment]){

//if (nphe > 60){
   cuts_empty = true;
   };
   }; 
  // };
   
   };
   };
   };
   };
  
   
   }; //fiducial
   }; //second part of sector 1
   
  if ((ph_EL >= 0) && (ph_EL <= 30)) {


if ((th_EL > th_min) && (th_EL < 50) && (ph_EL > fid_b) && (ph_EL < fid_a)){

 
//hist_z_el_1_empty->Fill(z_EL,1.);
//~fid cuts in Cherenkov plane
if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {
//cut on average number of photoelectrons
//if (avrg_nphe_sector1->GetBinContent(int((theta_cc+5.)*200./60.),int((ph_cc+30.)*200./60.)) > 70.) {
//nphe_sector1_after->Fill(nphe,1.);
if ((theta_cc >1.64786+2.16008*segment-0.00433375*segment*segment)&&(theta_cc <11.3593+1.38551*segment+0.0451710*segment*segment)){

if (nphe > ph_el_arr[pmt_hit+1][0][segment]){
//if (nphe > 60){
   cuts_empty = true; 
   };
  // };
  }; 
   };
   };
   };
   };

 }; //fiducial
   }; //first part of sector 1


// }; //nphe cut
//??? 

     
 
 }; // ectot vs p cut
 
 break;
 
case 2 : 
 
 if (((ECT/P_EL) > (-0.081058*P_EL*P_EL+0.240978*P_EL+0.0362218)) && ((ECT/P_EL) < (0.0211696*P_EL*P_EL-0.058201*P_EL)+0.44936)) {
  
// if (nphe > 25.) {

 if ((th_EL > th_min) && (th_EL < 50) && (ph_EL > fid_b+60) && (ph_EL < fid_a+60)){

// hist_z_el_2_empty->Fill(z_EL,1.);
 
 //~fid cuts in Cherenkov plane
  if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {

//cut on average number of photoelectrons
//if (avrg_nphe_sector2->GetBinContent(int((theta_cc+5.)*200./60.),int((ph_cc+30.)*200./60.)) > 70.) {
//nphe_sector2_after->Fill(nphe,1.);
if ((theta_cc >1.95003+2.00182*segment+0.00390572*segment*segment)&&(theta_cc <11.1869+1.37368*segment+0.047233*segment*segment)){

if (nphe > ph_el_arr[pmt_hit+1][1][segment]){
///if (nphe > 60){
   cuts_empty = true; 
   };
   };
 //  };
   
   };
   };
   };
   }; 
   
   }; //fiducial 

 //}; //nphe cut
 
 
  
 
 }; // ectot vs p cut

 break; 
 
case 3 : 
 
if (((ECT/P_EL) > (-0.059223*P_EL*P_EL+0.181728*P_EL+0.0558898)) && ((ECT/P_EL) < (0.0212801*P_EL*P_EL-0.0537573*P_EL)+0.384177)) { 
  
// if (nphe > 25.) {

if ((th_EL > th_min) && (th_EL < 50) && (ph_EL > fid_b+120) && (ph_EL < fid_a+120)){
//hist_z_el_3_empty->Fill(z_EL,1.);
 //~fid cuts in Cherenkov plane
  if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {

//cut on average number of photoelectrons
//if (avrg_nphe_sector3->GetBinContent(int((theta_cc+5.)*200./60.),int((ph_cc+30.)*200./60.)) > 70.) {
//nphe_sector3_after->Fill(nphe,1.);
if ((theta_cc >1.9837+2.0442*segment+0.0022649*segment*segment)&&(theta_cc <11.2718+1.35251*segment+0.047262*segment*segment)){

if (nphe > ph_el_arr[pmt_hit+1][2][segment]){

//if (nphe > 60){
   cuts_empty = true; 
   };
   };
 //  };
   
   };
   };
   };
   }; 
   
   }; //fiducial  

// }; //nphe cut
 
 
 
 }; // ectot vs p cut

 break;  
 
case 4 : 
 
 if (((ECT/P_EL) > (-0.0622566*P_EL*P_EL+0.191642*P_EL+0.0607197)) && ((ECT/P_EL) < (0.0140287*P_EL*P_EL-0.0442112*P_EL)+0.392145)) {
  
// if (nphe > 25.) {

if ((th_EL > th_min) && (th_EL < 50) && (ph_EL > fid_b+180) && (ph_EL < fid_a+180)){
//hist_z_el_4_empty->Fill(z_EL,1.);
//~fid cuts in Cherenkov plane
  if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {
//cut on average number of photoelectrons
//if (avrg_nphe_sector4->GetBinContent(int((theta_cc+5.)*200./60.),int((ph_cc+30.)*200./60.)) > 70.) {
//nphe_sector4_after->Fill(nphe,1.);
if ((theta_cc >1.36283+2.27581*segment-0.0137087*segment*segment)&&(theta_cc <11.7047+1.27222*segment+0.0517623*segment*segment)){
if (nphe > ph_el_arr[pmt_hit+1][3][segment]){

//if (nphe > 60){
   cuts_empty = true; 
   };
  }; 
   //};
   
   };
   };
   };
   };
   
   
     }; //fiducial 
     
        

 //}; //nphe cut
 
  
 }; // ectot vs p cut

 break;  
 
case 5 : 
 
 if (((ECT/P_EL) > (-0.072118*P_EL*P_EL+0.219742*P_EL+0.0372605)) && ((ECT/P_EL) < (0.0180325*P_EL*P_EL-0.0538702*P_EL)+0.409009)&&(PdHit_EL != 17)) {
   
 //if (nphe > 25.) {

if ((th_EL > th_min) && (th_EL < 50) && (ph_EL > fid_b+240) && (ph_EL < fid_a+240)){
//hist_z_el_5_empty->Fill(z_EL,1.);
//~fid cuts in Cherenkov plane
if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {
//cut on average number of photoelectrons
//if (avrg_nphe_sector5->GetBinContent(int((theta_cc+5.)*200./60.),int((ph_cc+30.)*200./60.)) > 70.) {
//nphe_sector5_after->Fill(nphe,1.);
if ((theta_cc >2.17772+1.95055*segment+0.00993131*segment*segment)&&(theta_cc <11.9184+1.34684*segment+0.0471248*segment*segment)){

if (nphe > ph_el_arr[pmt_hit+1][4][segment]){
//if (nphe > 60){
   cuts_empty = true; 
   };
   };
 //  };
   
   };
   };
   };
   };
   
     }; //fiducial     

 //}; //nphe cuts
 
 
 }; // ectot vs p cut

 break;   
 
case 6 : 

 
 if (((ECT/P_EL) > (-0.0504534*P_EL*P_EL+0.176674*P_EL+0.0627479)) && ((ECT/P_EL) < (0.0112167*P_EL*P_EL-0.0273361*P_EL)+0.382789)) {
  
 //if (nphe > 25.) {

if ((th_EL > th_min) && (th_EL < 50) && (ph_EL > fid_b+300) && (ph_EL < fid_a+300)){
//hist_z_el_6_empty->Fill(z_EL,1.);
//~fid cuts in Cherenkov plane
if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {

//cut on average number of photoelectrons
//if (avrg_nphe_sector6->GetBinContent(int((theta_cc+5.)*200./60.),int((ph_cc+30.)*200./60.)) > 70.) {
//nphe_sector6_after->Fill(nphe,1.);
if ((theta_cc >2.18831+2.01682*segment+0.00288938*segment*segment)&&(theta_cc <11.5751+1.25302*segment+0.0536114*segment*segment)){

if (nphe > ph_el_arr[pmt_hit+1][5][segment]){

//if (nphe > 60){
   cuts_empty = true; 
   };
   };
  // };
   
   };
   };
   };
   };
   
      }; //fiducial   

 //}; //nphe cut
     
 
 }; // ectot vs p cut
 break;      
 
   
   }; // end of switch
   };  // end of calorimeter threshold cut
  
   

   
   
   return cuts_empty;
   
   };
   
   
   
   
   
   
   
   
   ////////////////////////////////////////////
      bool cuts_empty::Proton_cuts_empty(){
       
   bool cuts_empty;
   Float_t m_p,p_fid_a,p_fid_b;
   m_p=0.938272;   
   p_fid_a = 24.*(1-exp(-1.*0.08*(th_P-9.)));
   p_fid_b = -1.*25.*(1-exp(-1.*0.1*(th_P-10.)));
   cuts_empty = false; 
   



if ((n_P == 1)&& (beta_P < 0.9496*P_P/sqrt(m_p*m_p+0.9497*P_P*P_P-0.06649) + 0.04136) && (beta_P > 1.045*P_P/sqrt(m_p*m_p+0.896*P_P*P_P - 0.2) - 0.139)){

//&& (beta_P < 0.9675*P_P/sqrt(m_p*m_p+0.9386*P_P*P_P-0.1723) + 0.0063) && (beta_P > 0.9408*P_P/sqrt(m_p*m_p+0.7455*P_P*P_P - 0.2544) - 0.1126)

if ((ph_P >= 330)&& (ph_P <= 360)&&(PdHit_P !=48) ) {

//&&( PdHit_P != 42)&&(PdHit_P != 45)&&(PdHit_P != 47)&&(PdHit_P !=48)
// &&( PdHit_P != 42)&&(PdHit_P != 45)&&(PdHit_P != 47)&&(PdHit_P !=48)
//if ((W >= 1.2) && (npart >= 3) && (npart <=4)) {

  if ((ph_P > p_fid_b+360) && (ph_P < p_fid_a+360)){

cuts_empty = true;
};//end of fiducial cut for the first part of the first sector
  
 //};//end of W-cut
 
 };//end of the first part of the first sector
 
 if ((ph_P >= 0) && (ph_P <= 30)&&(PdHit_P !=48)) {
 
 //&&( PdHit_P != 42)&&(PdHit_P != 45)&&(PdHit_P != 47)&&(PdHit_P !=48)
 //&&( PdHit_P != 42)&&(PdHit_P != 45)&&(PdHit_P != 47)&&(PdHit_P !=48)
//if ((W >= 1.2) && (npart >= 3) && (npart <=4)) {
if ((ph_P > p_fid_b) && (ph_P < p_fid_a)){


cuts_empty = true;
};//end of fiducial cut for the second part of the first sector
//};//end of W-cut
  };//end of the second part of the first sector
 
 
 
 if ((ph_P >= 30) && (ph_P <=90)&&(PdHit_P!=48)){ 
 //&&(PdHit_P!=24)&&(PdHit_P!=41)&&(PdHit_P!=42)&&(PdHit_P!=43)&&(PdHit_P!=44)&&(PdHit_P!=45)&&(PdHit_P!=46)&&(PdHit_P!=48)
 //&&(PdHit_P!=45)&&(PdHit_P!=46)&&(PdHit_P!=47)&&(PdHit_P!=41)&&(PdHit_P!=42)&&(PdHit_P!=43)&&(PdHit_P!=44)&&(PdHit_P!=48)
//if ((W >= 1.2) && (npart >= 3) && (npart <=4)) {
 if ((ph_P > p_fid_b+60) && (ph_P < p_fid_a+60)){

cuts_empty = true;
};//end of the fiducial cut for sector2
//};//end of W-cut
 
 };//end of the sector2
 
 if ((ph_P >=90) && (ph_P <=150)&&(PdHit_P!=44)&&(PdHit_P!=48)) {
 
 //&&(PdHit_P!=25)&&(PdHit_P!=40)&&(PdHit_P!=42)&&(PdHit_P!=44)&&(PdHit_P!=46)&&(PdHit_P!=48)
 //&&(PdHit_P!=40)&&(PdHit_P!=41)&&(PdHit_P!=42)&&(PdHit_P!=43)&&(PdHit_P!=46)&&(PdHit_P!=47)&&(PdHit_P!=44)&&(PdHit_P!=25)&&(PdHit_P!=48)
//if ((W >= 1.2) && (npart >= 3) && (npart <=4)){

if ((ph_P > p_fid_b+120) && (ph_P < p_fid_a+120)){

cuts_empty = true;
};//end of the fiducial cut for sector3
//}; //end of W-cut
  };//end of the sector3
 
 if ((ph_P >= 150) && (ph_P <= 210)&&(PdHit_P!=48)) {
 
//&&(PdHit_P!=42)&&(PdHit_P!=39)&&(PdHit_P!=48) 
 //&&(PdHit_P!=42)&&(PdHit_P!=43)&&(PdHit_P!=46)&&(PdHit_P!=39)&&(PdHit_P!=48)
//if ((W >= 1.2) && (npart >= 3) && (npart <=4)){
if ((ph_P > p_fid_b+180) && (ph_P < p_fid_a+180)){

cuts_empty = true;
};//end of the fiducial cut for sector4
//}; //end of W-cut

};//end of the sector4

if ((ph_P >= 210) && (ph_P <=270)&&(PdHit_P!=17)&&(PdHit_P!=48)) {

//&&(PdHit_P!=40)&&(PdHit_P!=42)&&(PdHit_P!=43)&&(PdHit_P!=17)&&(PdHit_P!=48)
//&&(PdHit_P!=40)&&(PdHit_P!=42)&&(PdHit_P!=44)&&(PdHit_P!=46)&&(PdHit_P!=47)&&(PdHit_P!=17)&&(PdHit_P!=48)
//if ((W >= 1.2) && (npart >= 3) && (npart <=4)) {
 if ((ph_P > p_fid_b+240) && (ph_P < p_fid_a+240)){

 cuts_empty = true;
 };//end of the fiducial cut for sector5
// };//end of W-cut
 
 };//end of the sector5
 
 if ((ph_P >= 270) && (ph_P <=330)&&(PdHit_P!=48)) {
 
 //&&(PdHit_P!=40)&&(PdHit_P!=45)&&(PdHit_P!=31)&&(PdHit_P!=48)&&(PdHit_P!=47)
 //&&(PdHit_P!=40)&&(PdHit_P!=43)&&(PdHit_P!=44)&&(PdHit_P!=45)&&(PdHit_P!=31)&&(PdHit_P!=48)&&(PdHit_P!=47)
//if ((W >= 1.2) && (npart >= 3) && (npart <=4)){
if ((ph_P > p_fid_b+300) && (ph_P < p_fid_a+300)){

cuts_empty = true;
 };//end of the fiducial cut for sector6
// };//end of W-cut
 };//end of the sector6
 
 
 
 };
 

    return cuts_empty;
    };
    
    
    
    /////////////////////
    
    
   bool cuts_empty::PIp_cuts_empty(){
       
   bool cuts_empty;
   Float_t m_pip,pip_fid_a,pip_fid_b,beta_nom_pip;
   m_pip = 0.13957;
  pip_fid_a = 24.*(1-exp(-1.*0.08*(th_PIp-9.)));
  pip_fid_b = -1.*25.*(1-exp(-1.*0.1*(th_PIp-10.)));
  beta_nom_pip = P_PIp/sqrt(m_pip*m_pip+P_PIp*P_PIp);
   cuts_empty = false; 
     //cout << th_PIp<< " rgdgdf "<<ph_PIp<<" iiiiii "<<beta_PIp<<" riuthy "<<P_PIp<<"\n";
   if ((n_PIp == 1)&& (PIp_dist*(1./beta_nom_pip-1/beta_PIp)/30. < 0.0001769/(P_PIp*P_PIp*P_PIp*P_PIp+0.0001471)+0.8465)&&(PIp_dist*(1./beta_nom_pip-1/beta_PIp)/30. > -0.0002121/(P_PIp*P_PIp*P_PIp*P_PIp+5.685e-05)-0.8411)) {
   
   
   //(beta_PIp < 0.8*P_PIp/sqrt(m_pip*m_pip+0.9108*P_PIp*P_PIp-0.001768) + 0.2) && (beta_PIp > 1.054*P_PIp/sqrt(m_pip*m_pip+0.7001*P_PIp*P_PIp - 0.006497) - 0.2999)
   //&& (beta_PIp < 0.8*P_PIp/sqrt(m_pip*m_pip+0.91*P_PIp*P_PIp-0.0034) + 0.2) && (beta_PIp > 1.05*P_PIp/sqrt(m_pip*m_pip+0.7*P_PIp*P_PIp - 0.0056) - 0.297)
   
   //&& (beta_PIp < 0.8*P_PIp/sqrt(m_pip*m_pip+0.91*P_PIp*P_PIp-0.0034) + 0.2) && (beta_PIp > 1.05*P_PIp/sqrt(m_pip*m_pip+0.7*P_PIp*P_PIp - 0.0056) - 0.297)
   
     
   if ((ph_PIp >= 330) && (ph_PIp <=360)&&(PdHit_PIp !=48)){
   
   //&&(PdHit_PIp != 42)&&(PdHit_PIp != 45)&&(PdHit_PIp != 46)&&(PdHit_PIp != 47)&&(PdHit_PIp !=48)
   //&&(PdHit_PIp != 42)&&(PdHit_PIp != 45)&&(PdHit_PIp != 47)&&(PdHit_PIp !=48)
   
   if ((ph_PIp > pip_fid_b+360) && (ph_PIp < pip_fid_a+360)){
 
   cuts_empty = true;
   
   };//end of the fiducial cut for first part of sector1
  }; //end of the first part of sector1
  
  if ((ph_PIp >= 0) && (ph_PIp <=30)&&(PdHit_PIp !=48)){
  
  //&&(PdHit_PIp != 42)&&(PdHit_PIp != 45)&&(PdHit_PIp != 46)&&(PdHit_PIp != 47)&&(PdHit_PIp !=48)
  //&&( PdHit_PIp != 42)&&(PdHit_PIp != 45)&&(PdHit_PIp != 47)&&(PdHit_PIp !=48)

 if ((ph_PIp > pip_fid_b) && (ph_PIp < pip_fid_a)){

cuts_empty = true;
   };//end of the fiducial cut for second part of sector1
   }; //end of the second part of sector1
  
  
  if ((ph_PIp >= 30) && (ph_PIp <=90)&&(PdHit_PIp !=48)) {
  
 //&&(PdHit_PIp!=24)&&(PdHit_PIp!=45)&&(PdHit_PIp!=46)&&(PdHit_PIp!=47)&&(PdHit_PIp!=41)&&(PdHit_PIp!=42)&&(PdHit_PIp!=43)&&(PdHit_PIp!=44)&&(PdHit_PIp!=48) //&&(PdHit_PIp!=45)&&(PdHit_PIp!=46)&&(PdHit_PIp!=47)&&(PdHit_PIp!=41)&&(PdHit_PIp!=42)&&(PdHit_PIp!=43)&&(PdHit_PIp!=44)&&(PdHit_PIp!=48)

if ((ph_PIp > pip_fid_b+60) && (ph_PIp < pip_fid_a+60)){

cuts_empty = true;
 };//end of the fiducial cut for sector2
};//end of the sector2
  
  
if ((ph_PIp >= 90) && (ph_PIp <=150)&&(PdHit_PIp!=44)&&(PdHit_PIp!=48)) {

//&&(PdHit_PIp!=40)&&(PdHit_PIp!=41)&&(PdHit_PIp!=42)&&(PdHit_PIp!=43)&&(PdHit_PIp!=46)&&(PdHit_PIp!=47)&&(PdHit_PIp!=44)&&(PdHit_PIp!=25)&&(PdHit_PIp!=48)
//&&(PdHit_PIp!=40)&&(PdHit_PIp!=41)&&(PdHit_PIp!=42)&&(PdHit_PIp!=43)&&(PdHit_PIp!=46)&&(PdHit_PIp!=47)&&(PdHit_PIp!=44)&&(PdHit_PIp!=25)&&(PdHit_PIp!=48)

if ((ph_PIp > pip_fid_b+120) && (ph_PIp < pip_fid_a+120)){

cuts_empty = true;
 };//end of the fiducial cut for sector3
};  //end of the sector3
  
 if ((ph_PIp >= 150) && (ph_PIp <=210)&&(PdHit_PIp!=48)){
 
 //&&(PdHit_PIp!=42)&&(PdHit_PIp!=43)&&(PdHit_PIp!=46)&&(PdHit_PIp!=39)&&(PdHit_PIp!=48)
 //&&(PdHit_PIp!=42)&&(PdHit_PIp!=43)&&(PdHit_PIp!=46)&&(PdHit_PIp!=39)&&(PdHit_PIp!=48)
 
 if ((ph_PIp > pip_fid_b+180) && (ph_PIp < pip_fid_a+180)){

cuts_empty = true;
 };//end of the fiducial cut for sector4
 
 }; //end of the sector4
  
if ((ph_PIp >= 210) && (ph_PIp <=270)&&(PdHit_PIp!=17)&&(PdHit_PIp!=48)) {

//&&(PdHit_PIp!=40)&&(PdHit_PIp!=42)&&(PdHit_PIp!=44)&&(PdHit_PIp!=46)&&(PdHit_PIp!=47)&&(PdHit_PIp!=17)&&(PdHit_PIp!=48)
//&&(PdHit_PIp!=40)&&(PdHit_PIp!=42)&&(PdHit_PIp!=44)&&(PdHit_PIp!=46)&&(PdHit_PIp!=47)&&(PdHit_PIp!=17)&&(PdHit_PIp!=48)

if ((ph_PIp > pip_fid_b+240) && (ph_PIp < pip_fid_a+240)){


cuts_empty = true;
 };//end of the fiducial cut for sector5

}; //end of the sector5
  
  
 if ((ph_PIp >= 270) && (ph_PIp <=330)&&(PdHit_PIp!=48)){
 
 //&&(PdHit_PIp!=40)&&(PdHit_PIp!=44)&&(PdHit_PIp!=45)&&(PdHit_PIp!=31)&&(PdHit_PIp!=48)&&(PdHit_PIp!=47)
 //&&(PdHit_PIp!=40)&&(PdHit_PIp!=43)&&(PdHit_PIp!=44)&&(PdHit_PIp!=45)&&(PdHit_PIp!=31)&&(PdHit_PIp!=48)&&(PdHit_PIp!=47)

 if ((ph_PIp > pip_fid_b+300) && (ph_PIp < pip_fid_a+300)){

cuts_empty = true;
 };//end of the fiducial cut for sector6
 
 }; //end of the sector6
  
  
  
  
   };
    return cuts_empty;
   };
   
   ////////////////////
      
   bool cuts_empty::PIm_cuts_empty(){
       
   bool cuts_empty;
   Float_t m_pim,th_min,par1,par2,pim_fid_a,pim_fid_b,beta_nom_pim;
   
    
   m_pim = 0.13957; 
  th_min=(11.09+8./(0.472*P_PIm+0.117));
  par1=0.705+1.1*P_PIm;
  par2=-63.2-29.3*P_PIm;       
   pim_fid_a=30.5*pow((sin((th_PIm-th_min)*0.01745)),(par1+par2/th_PIm+1530./th_PIm/th_PIm))-1;
    beta_nom_pim = P_PIm/sqrt(m_pim*m_pim+P_PIm*P_PIm); 
   pim_fid_b=-30.5*pow((sin((th_PIm-th_min)*0.01745)),(par1+par2/th_PIm+1530./th_PIm/th_PIm))+1; 
   cuts_empty = false; 
  
     
   
 
   if ((n_PIm == 1)&& (PIm_dist*(1./beta_nom_pim-1/beta_PIm)/30. < 0.001284/(P_PIm*P_PIm*P_PIm*P_PIm+0.00206)+0.6335)&&(PIm_dist*(1./beta_nom_pim-1/beta_PIm)/30. > -0.002811/(P_PIm*P_PIm*P_PIm*P_PIm+0.002855)-0.5364)) {
   
   //&& (beta_PIm < 0.8*P_PIm/sqrt(m_pim*m_pim+1.2*P_PIm*P_PIm+0.008756) + 0.02963) && (beta_PIm > 0.7686*P_PIm/sqrt(m_pim*m_pim+0.4188*P_PIm*P_PIm - 0.01057)-0.2111)
   
   //&& (beta_PIm < 0.2338*P_PIm/sqrt(m_pim*m_pim+0.0518*P_PIm*P_PIm-0.0187) + 0.0001728) && (beta_PIm > 0.9659*P_PIm/sqrt(m_pim*m_pim+0.9729*P_PIm*P_PIm + 0.008634)-0.0003043)
   //&& (beta_PIm < 0.1717*P_PIm/sqrt(m_pim*m_pim+0.028*P_PIm*P_PIm-0.02) + 0.00023) && (beta_PIm > 0.112*P_PIm/sqrt(m_pim*m_pim+0.015*P_PIm*P_PIm - 0.02))
   
   
 if ((ph_PIm >= 330) && (ph_PIm <=360)&&(PdHit_PIm!=48)){
 
 //&&(PdHit_PIm!=42)&&(PdHit_PIm!=45)&&(PdHit_PIm!=46)&&(PdHit_PIm!=47)&&(PdHit_PIm!=48)
 //&&(PdHit_PIm!=42)&&(PdHit_PIm!=45)&&(PdHit_PIm!=46)&&(PdHit_PIm!=47)&&(PdHit_PIm!=48)



  if ((ph_PIm > pim_fid_b+360) && (ph_PIm < pim_fid_a+360)){
 
   cuts_empty = true;
   
   };//end of the fiducial cut for first part of sector1
  }; //end of the first part of sector1
  
  if ((ph_PIm >= 0) && (ph_PIm <=30)&&(PdHit_PIm!=48)){
  
  //&&(PdHit_PIm!=42)&&(PdHit_PIm!=45)&&(PdHit_PIm!=46)&&(PdHit_PIm!=47)&&(PdHit_PIm!=48)
  //&&(PdHit_PIm!=42)&&(PdHit_PIm!=45)&&(PdHit_PIm!=46)&&(PdHit_PIm!=47)&&(PdHit_PIm!=48)

 
 if ((ph_PIm > pim_fid_b) && (ph_PIm < pim_fid_a)){

cuts_empty = true;
   };//end of the fiducial cut for second part of sector1
   };//end of the second part of sector1
  
    if ((ph_PIm >= 30) && (ph_PIm <=90)&&(PdHit_PIm!=48)) {
    
 //  &&(PdHit_PIm!=16)&&(PdHit_PIm!=24)&&(PdHit_PIm!=41)&&(PdHit_PIm!=42)&&(PdHit_PIm!=43)&&(PdHit_PIm!=44)&&(PdHit_PIm!=45)&&(PdHit_PIm!=46)&&(PdHit_PIm!=47)&&(PdHit_PIm!=48) //&&(PdHit_PIm!=16)&&(PdHit_PIm!=24)&&(PdHit_PIm!=41)&&(PdHit_PIm!=42)&&(PdHit_PIm!=43)&&(PdHit_PIm!=44)&&(PdHit_PIm!=45)&&(PdHit_PIm!=46)&&(PdHit_PIm!=47)&&(PdHit_PIm!=48)


if ((ph_PIm > pim_fid_b+60) && (ph_PIm < pim_fid_a+60)){

cuts_empty = true;
 };//end of the fiducial cut for sector2
};//end of the sector2
    
if ((ph_PIm >= 90) && (ph_PIm <=150)&&(PdHit_PIm!=44)&&(PdHit_PIm!=48)) {

//&&(PdHit_PIm!=25)&&(PdHit_PIm!=40)&&(PdHit_PIm!=41)&&(PdHit_PIm!=42)&&(PdHit_PIm!=43)&&(PdHit_PIm!=44)&&(PdHit_PIm!=46)&&(PdHit_PIm!=47)&&(PdHit_PIm!=48)
//&&(PdHit_PIm!=25)&&(PdHit_PIm!=40)&&(PdHit_PIm!=41)&&(PdHit_PIm!=42)&&(PdHit_PIm!=43)&&(PdHit_PIm!=44)&&(PdHit_PIm!=46)&&(PdHit_PIm!=47)&&(PdHit_PIm!=48)

 
 if ((ph_PIm > pim_fid_b+120) && (ph_PIm < pim_fid_a+120)){

cuts_empty = true;
 };//end of the fiducial cut for sector3
};  //end of the sector3
  
 if ((ph_PIm >= 150) && (ph_PIm <=210)&&(PdHit_PIm!=48)){
 
 //&&(PdHit_PIm!=39)&&(PdHit_PIm!=42)&&(PdHit_PIm!=43)&&(PdHit_PIm!=46)&&(PdHit_PIm!=47)&&(PdHit_PIm!=48)
 //&&(PdHit_PIm!=39)&&(PdHit_PIm!=42)&&(PdHit_PIm!=43)&&(PdHit_PIm!=46)&&(PdHit_PIm!=47)&&(PdHit_PIm!=48)
 
 
 if ((ph_PIm > pim_fid_b+180) && (ph_PIm < pim_fid_a+180)){

cuts_empty = true;
 };//end of the fiducial cut for sector4
 }; //end of the sector4
  
if ((ph_PIm >= 210) && (ph_PIm <=270)&&(PdHit_PIm!=17)&&(PdHit_PIm!=48)) {

//&&(PdHit_PIm!=17)&&(PdHit_PIm!=40)&&(PdHit_PIm!=42)&&(PdHit_PIm!=43)&&(PdHit_PIm!=44)&&(PdHit_PIm!=46)&&(PdHit_PIm!=47)&&(PdHit_PIm!=48)
//&&(PdHit_PIm!=17)&&(PdHit_PIm!=40)&&(PdHit_PIm!=42)&&(PdHit_PIm!=43)&&(PdHit_PIm!=44)&&(PdHit_PIm!=46)&&(PdHit_PIm!=47)&&(PdHit_PIm!=48)

if ((ph_PIm > pim_fid_b+240) && (ph_PIm < pim_fid_a+240)){

cuts_empty = true;
 };//end of the fiducial cut for sector5
}; //end of the sector5
    
 if ((ph_PIm >= 270) && (ph_PIm <=330)&&(PdHit_PIm!=48)){
 
 //&&(PdHit_PIm!=31)&&(PdHit_PIm!=40)&&(PdHit_PIm!=43)&&(PdHit_PIm!=44)&&(PdHit_PIm!=45)&&(PdHit_PIm!=46)&&(PdHit_PIm!=47)&&(PdHit_PIm!=48) //&&(PdHit_PIm!=31)&&(PdHit_PIm!=40)&&(PdHit_PIm!=43)&&(PdHit_PIm!=44)&&(PdHit_PIm!=45)&&(PdHit_PIm!=46)&&(PdHit_PIm!=47)&&(PdHit_PIm!=48)

 
 if ((ph_PIm > pim_fid_b+300) && (ph_PIm < pim_fid_a+300)){

cuts_empty = true;
 };//end of the fiducial cut for sector6
 }; //end of the sector6
  
   
  
  
   };
    return cuts_empty;
   }; 
   
   
