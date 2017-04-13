#include "TROOT.h"
#include "TFile.h"
#include "TLine.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "THnSparse.h"
#include "TH3.h"
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
#include <TGClient.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <TGButtonGroup.h>
#include <RQ_OBJECT.h>
#include <TGNumberEntry.h>
#include <TGProgressBar.h>
#include <TGLabel.h>
#include <stdio.h>
#include <dlfcn.h>
#include "MyMainFrame.h"
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <TGFileDialog.h>
#include <GuiTypes.h>
#include <TGDoubleSlider.h>
#include <TGComboBox.h>
#include <TLeaf.h>
#include <TBranch.h>
#include <TLorentzVector.h>
#include <TError.h> 
#include <auto_ptr.h>
#ifndef __CINT__
#include <cstdlib>
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
#include <RooLandau.h> 
#include <RooNumConvPdf.h>
#include <RooDataHist.h>
#include "RooBinning.h"
#include <sys/types.h>
#include <wait.h>
#include <unistd.h>
#include <cstring>
#include <getopt.h>
#include <cstdlib>
#include "global.h"
#include "beta_func_data.h"
#include "beta_func_empty.h"
#include "rot_boost_cmsyst.h"
#include "fermi_bonn.h"
#include "data_hist.h"
#include "sim_hist.h"
#include "output.h"



 using namespace std; 
 

#define _USE_MATH_DEFINES


    void MyMainFrame::MainFrame(UChar_t flag, Float_t E_beam, Short_t nfiles, Short_t nfiles_empty, Short_t nfiles_sim, string inp_files[],string inp_files_empty[], string inp_files_sim[], string outfile_in) { 


	inpfile_inp = inp_files[0];
	outfile_inp = outfile_in;
	
	E0 = E_beam;

        n_files = nfiles;
	n_files_empty = nfiles_empty;
	n_files_sim = nfiles_sim;

//sozdaem massiv strok (novie)
        file = new string[nfiles];
	file_empty = new string[nfiles_empty];
	file_sim = new string[nfiles_sim];

//prisvaevaem novim massivam strok massivi strok iz input faila 
        file = inp_files;
	file_empty = inp_files_empty;
	file_sim = inp_files_sim;
	data_sim = flag;
	

DoDraw();


   
}

void MyMainFrame::DoDraw() {

 
 
       t20tot21(); 
   
}








    void MyMainFrame::t20tot21() { 
     gROOT->SetBatch(true);
gROOT->ProcessLine( "gErrorIgnoreLevel = kWarning; " );
  
    ostringstream qqq;
    ostringstream qqq1;
    
	Float_t Pgen,Prec; 
    Long64_t j;

 Short_t m, ti;
 Long64_t i,nstart,nstop,n_incl,n_elast, k_long;
 
 
 Float_t arr_pim_mis_cuts_min[12][12]={{0.,0.,0.,0.,0.,0.,-0.15,-0.25,-0.25,0.,0.,0.},{0.,0.,0.,0.,0.,0.,-0.25,-0.2,0.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,-0.3,-0.25,0.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,-0.2,-0.2,-0.2,0.,0.,0.},{0.,0.,0.,0.,0.,-0.2,-0.07,0.,0.,0.,0.,0},{0.,0.,0.,0.,-0.1,-0.1,0.,0.,0.,0.,0.,0.},{0.,0.,0.,0.,-0.1,-0.1,-0.1,0.,0.,0.,0.,0.},{0.,0.,-0.25,-0.3,-0.12,-0.15,0.,0.,0.,0.,0.,0.},{0.,0.,-0.2,-0.1,0.,-0.1,-0.1,0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},{0.,0.,0.,-0.2,0.,0.,0.,0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}};
 
 
 Float_t arr_pim_mis_cuts_max[12][12]={{0.,0.,0.,0.,0.,0.,0.15,0.4,0.3,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.25,0.35,0.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.25,0.35,0.3,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.15,0.17,0.25,0.,0.,0.},{0.,0.,0.,0.,0.,0.2,0.2,0.,0.,0.,0.,0},{0.,0.,0.,0.,0.25,0.1,0.,0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.2,0.1,0.25,0.,0.,0.,0.,0.},{0.,0.,0.15,0.2,0.25,0.2,0.,0.,0.,0.,0.,0.},{0.,0.,0.2,0.2,0.,0.08,0.,0.,0.,0.,0.,0.},{0.,0.,0.,0.25,0.25,0.25,0.,0.,0.,0.,0.,0.},{0.,0.3,0.3,0.3,0.3,0.,0.,0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}};
 
 
/* Float_t ARR_MMC_PIM_MISS[5][8][8] = {{{1.,1.,1.,1.,1.,1.,1.,1.},{1.,1.,1.,1.,1.,0.,0.,0.},{1.,1.,1.,1.,0.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.,0.,0.,0.,1.,1.,1.,1.},{0.,0.,0.,1.,1.,1.,1.,1.},{1.,1.,1.,1.,1.,1.,1.,1.},{1.,1.,1.,1.,1.,1.,1.,0.},{1.,1.,1.,1.,1.,1.,0.,0.},{1.,1.,1.,1.,0.,0.,0.,0.},{0.,1.,1.,0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.,0.,0.,0.,0.,1.,1.,1.},{0.,0.,0.,0.,0.,1.,1.,1.},{0.,0.,0.,0.,1.,1.,1.,1.},{0.,1.,1.,1.,1.,1.,1.,1.},{1.,1.,1.,1.,1.,1.,1.,0.},{1.,1.,1.,1.,1.,1.,0.,0.},{1.,1.,1.,1.,0.,0.,0.,0.},{1.,1.,0.,0.,0.,0.,0.,0.}},{{0.,0.,0.,0.,1.,1.,1.,1.},{0.,0.,0.,0.,1.,1.,1.,1.},{0.,0.,0.,0.,1.,1.,1.,1.},{0.,0.,0.,1.,1.,1.,1.,0.},{0.,0.,0.,1.,1.,1.,1.,0.},{0.,0.,1.,1.,1.,1.,0.,0.},{0.,0.,1.,1.,1.,0.,0.,0.},{1.,1.,1.,0.,0.,0.,0.,0.}},{{0.,0.,0.,1.,1.,1.,0.,0.},{0.,0.,0.,1.,1.,1.,0.,0.},{0.,0.,0.,1.,1.,1.,0.,0.},{0.,0.,0.,1.,1.,1.,0.,0.},{0.,0.,1.,1.,1.,1.,0.,0.},{0.,1.,1.,1.,1.,0.,0.,0.},{1.,1.,1.,1.,0.,0.,0.,0.},{1.,1.,1.,1.,0.,0.,0.,0.}}};*/

Float_t ARR_MMC_PIM_MISS[5][8][8] = {{{1.,1.,1.,1.,1.,1.,1.,1.},{1.,1.,1.,1.,1.,1.,1.,1.},{1.,1.,1.,1.,1.,1.,1.,1.},{0.,0.,0.,0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.,0.,1.,1.,1.,1.,1.,1.},{1.,1.,1.,1.,1.,1.,1.,1.},{1.,1.,1.,1.,1.,1.,1.,1.},{1.,1.,1.,1.,1.,1.,1.,1.},{1.,1.,1.,1.,1.,1.,1.,0.},{1.,1.,1.,1.,1.,1.,0.,0.},{1.,1.,1.,1.,1.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.,0.,0.,1.,1.,1.,1.,1.},{0.,0.,0.,0.,1.,1.,1.,1.},{0.,0.,1.,1.,1.,1.,1.,1.},{1.,1.,1.,1.,1.,1.,1.,1.},{1.,1.,1.,1.,1.,1.,1.,0.},{1.,1.,1.,1.,1.,1.,0.,0.},{1.,1.,1.,1.,1.,0.,0.,0.},{1.,1.,0.,0.,0.,0.,0.,0.}},{{0.,0.,1.,1.,1.,1.,1.,1.},{0.,0.,0.,1.,1.,1.,1.,1.},{0.,1.,1.,1.,1.,1.,1.,1.},{1.,1.,1.,1.,1.,1.,1.,0.},{1.,1.,1.,1.,1.,1.,1.,0.},{1.,1.,1.,1.,1.,1.,0.,0.},{1.,1.,1.,1.,1.,1.,0.,0.},{1.,1.,1.,1.,1.,0.,0.,0.}},{{1.,1.,1.,1.,1.,1.,0.,0.},{1.,1.,1.,1.,1.,1.,0.,0.},{1.,1.,1.,1.,1.,1.,0.,0.},{1.,1.,1.,1.,1.,1.,0.,0.},{1.,1.,1.,1.,1.,1.,0.,0.},{1.,1.,1.,1.,1.,0.,0.,0.},{1.,1.,1.,1.,1.,0.,0.,0.},{1.,1.,1.,1.,0.,0.,0.,0.}}};

//Float_t MMcut_pim_miss[21] = {0.2,0.2,0.2,0.25,0.25,0.25,0.25,0.22,0.22,0.22,0.22,0.27,0.3,0.27,0.3,0.25,0.27,0.27,0.27,0.27,0.27};
Float_t MMcut_pim_miss[21] = {0.2,0.22,0.25,0.26,0.26,0.3,0.27,0.25,0.25,0.24,0.24,0.27,0.3,0.27,0.3,0.27,0.3,0.3,0.3,0.3,0.3};

Float_t fract_fsi_corr[12][21] = {{1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.},{1., 1., 1., 1., 1., 1., 1., 0.949666, 0.946763, 0.92825, 0.92762, 0.942042, 0.969058, 0.939906, 0.972951, 0.920794, 0.938765, 0.952647, 0.943091, 0.948361, 0.963557},{1., 1., 1., 1., 1., 1., 1., 0.993911, 0.893549, 0.970045, 0.961446, 0.859731, 0.911005, 0.949457, 0.897599, 0.987884, 0.926785, 0.961528, 0.962106, 0.943108, 1.},{1., 1., 1., 1., 1., 1., 1., 0.994626, 0.92888, 0.974134, 0.938982, 0.912881, 0.886699, 0.95578, 0.953117, 0.980564, 0.935179, 0.955753, 0.971683, 0.966479, 1.},{1., 1., 1., 1., 1., 1., 1., 0.976887, 0.915209, 0.920231, 0.971177, 0.938058, 0.90531, 0.963204, 0.936656, 0.988965, 0.929378, 0.957025, 0.987681, 1., 1.},{1., 1., 1., 1., 1., 1., 1., 0.946144, 0.960353, 0.98249, 0.9808, 0.970148, 0.945595, 0.969558, 0.963168, 0.979624, 0.95783, 0.953134, 1., 1., 1.},{1., 1., 1., 1., 1., 1., 1., 0.999138, 0.989533, 0.970615, 0.99736, 0.93694, 0.938233, 0.948618, 0.958399, 0.977104, 1., 1., 1., 1., 1.},{1., 1., 1., 1., 1., 1., 1., 0.995786, 0.941906, 0.954921, 0.991707, 0.943651, 0.936837, 0.950524, 0.941602, 0.983574, 1., 1., 1., 1., 1.},{1., 1., 1., 1., 1., 1., 1., 0.995372, 0.992961, 0.996496, 0.992185, 0.956822, 0.943974, 0.88652, 1., 1., 1., 1., 1., 1., 1.},{1., 1., 1., 1., 1., 1., 1., 0.991166, 0.998852, 0.99667, 0.981678, 0.9876, 1., 1., 1., 1., 1., 1., 1., 1., 1.},{1., 1., 1., 1., 1., 1., 1., 0.9561, 0.976256, 0.982884, 0.985384, 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.},{1., 1., 1., 1., 1., 1., 1., 0.986958, 0.998352, 0.981677, 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.}} ;
TH1I *hist_adc_off[12]; 

 TFile *finp;
 Int_t block_total = 0;
 Int_t block_last = 0;
 Float_t Qfull = 0.;
  Float_t Qfull_empty = 0.;
  Float_t Qfull_sim = 0.; 
 
  Int_t block_curr = 0;
 Int_t block_prev = 0;
  Short_t m_old=0;
  
  
  
  
   Float_t P_EL_old,th_EL_old,Q2_old,P_EL_new,th_EL_new,W_new,Q2_new;
 bool selection;
 bool selection_pim_miss, selection_0_miss;
 bool selection_pim_miss_sim, selection_0_miss_sim;
 bool selection_pim_miss_empt, selection_0_miss_empt; 
 

 
  
	UChar_t pdhit;
	
	
	Float_t * p; 
        p = new Float_t [20];
	
        Float_t m_deutron,NpheCC_EL,Nphe_pip,Nphe_pim,ECtot_EL;
	Float_t p_miss_1,p_miss_2,p_miss_3,p_miss_4,p_miss_sqr;
	Float_t beta_nom_pip,beta_nom_pim,beta_nom_p,delta_t_pip,delta_t_pim,delta_t_p;
	Float_t p_fid_a_1, p_fid_b_1;
	Float_t  x_EL,y_EL;
	Float_t th_PIm_miss,ph_PIm_miss;
	Float_t  ECin_EL,ECout_EL;
	Int_t block, block_tot;
        Long64_t gpart,k,last_i,last_k;
	Float_t  q_l,Qdiff,Qcurr,Qprev,Qtotal,deltaQ;
	Int_t sc_part_local;
	Float_t sc_pd_local,sc_sect_local;
	Float_t sc_z,fid_a,fid_b,a,b;
	Float_t nx,ny,nz,par1,par2,th_min;
	Float_t sx,sy,sz,px,py,pz;
	Float_t delta_mom_p_ye, delta_p_el_sim; 
	Float_t E_gamma,E_p_gamma_lab,P_p_gamma_lab,beta1,gamma;
	Float_t pip_fid_a_1,pip_fid_b_1;
	Float_t th_min_1,par1_1,par2_1, pim_fid_a_1,pim_fid_b_1;
	Float_t W_old;
	Float_t pf_x,pf_y,pf_z,pxel_new,pyel_new,pzel_new;
	Float_t th_ph_pim, th_ph_pip,th_ph_pr;
	Float_t dc_x_el, dc_y_el;
	
	
   	Double_t integ, err, err1, old_bin_cont, new_bin_cont;
	
	Double_t Var1[5],Var2[5],Var3[5]; 
	Double_t Var_1[5],Var_2[5],Var_3[5]; 
	
 	bool cut_fiduch;
	bool bool_el_id_data, bool_proton_id_data, bool_pip_id_data,bool_pim_id_data;
	bool bool_el_id_sim, bool_proton_id_sim, bool_pip_id_sim,bool_pim_id_sim;
	bool bool_el_id_empt, bool_proton_id_empt, bool_pip_id_empt,bool_pim_id_empt;

	
	m_proton = 0.938272;
	//m_proton = 0.93957;
	m_deutron = 1.875612;
	m_pip = 0.13957;
	  
  
	TLorentzVector P4_D,P4_PIm_miss,P4_PIp_miss,P4_PIp_miss_d,P4_P_miss,P4_PIm_miss_0,P4_Pini_ferm_dat;
	TLorentzVector P4_PP_rot,P4_PP_rot_1,P4_PP_rot_2,P4_PP_rot_3,P4_PP_rot_2_boost;
	TLorentzVector  P4_PP_cor, P4_PIp_cor, P4_PIm_cor;
	TLorentzVector  P4_miss_0,P4_miss_0_d,P4_PIm_miss_d;
	TLorentzVector  P4_miss_0_en_comp,P4_PIm_miss_en_comp;
	TLorentzVector  P4_ELP_for_miss_en_comp;
	
	TVector3 V3_dir_gamma,P3_PP_rot,uz,ux;
		
	Float_t th_gamma, phi_gamma;
		 
	P4_EL.SetXYZT(0,0,2.039,2.039);
	P4_P.SetXYZT(0,0,0,m_proton);
	P4_D.SetXYZT(0,0,0,m_deutron);
	 
        Float_t p0_simomcor, p1_simomcor, p2_simomcor;
	


global();

 Float_t mms_all_reg_arr[2][5][20] = {{{-0.0,-0.2,-0.0,-0.08,-0.12,-0.12,0.,-0.12,-0.04,0.,-0.0,0.,-0.0,0.,0.,0.,0.,0.,0.,0.},
                             {-0.0,-0.08,-0.12,-0.05,-0.04,-0.0,0.,-0.05,-0.05,0.,-0.0,0.06,-0.12,1000.,1000.,1000.,1000.,1000.,1000.,1000.},
                             {1000.,-0.08,-0.1,-0.06,-0.05,-0.14,-0.08,-0.12,-0.09,-0.05,-0.08,-0.12,-0.12,-0.08,1000.,1000.,1000.,1000.,1000.,1000.},
                             {-0.04,-0.04,-0.15,-0.17,-0.08,-0.08,-0.08,-0.04,-0.08,-0.08,-0.04,-0.04,-0.16,-0.0,-0.08,-0.0,-0.12,-0.0,1000.,1000.},
                             {0.,0.,-0.2,-0.07,-0.07,-0.1,-0.12,-0.0,-0.08,-0.15,-0.08,-0.07,-0.13,-0.13,-0.05,-0.13,1000.,1000.,1000.,1000.}},
			     {{0.0,0.27,0.27,0.27,0.27,0.27,0.27,0.27,0.27,0.0,0.0,0.0,0.0,0.,0.,0.,0.,0.,0.,0.},
                             {0.17,0.18,0.22,0.18,0.2,0.18,0.19,0.3,0.3,0.3,0.25,0.25,0.25,-1000.,-1000.,-1000.,-1000.,-1000.,-1000.,-1000.},
	                     {-1000.,0.19,0.2,0.25,0.24,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,-1000.,-1000.,-1000.,-1000.,-1000.,-1000.},
		             {0.23,0.19,0.24,0.23,0.18,0.2,0.21,0.19,0.2,0.18,0.18,0.18,0.23,0.19,0.23,0.19,0.23,0.15,-1000.,-1000.},
			     {0.24,0.2,0.28,0.24,0.2,0.2,0.2,0.2,0.2,0.22,0.2,0.25,0.22,0.25,0.25,0.5,-1000.,-1000.,-1000.,-1000.}}};

//        TFile *nphefile = new TFile("norm_nphe_15jan17_gt_50.root","READ");
//	TFile *nphefile = new TFile("norm_nphe_17jan17_indiv_seg.root","READ");
//	 TFile *nphefile = new TFile("norm_nphe_26jan17_aft_ph_th_match.root","READ");
	 TFile *nphefile = new TFile("norm_nphe_27jan17_aft_cc_match_gt50.root","READ");
	 
        norm_nphe_s1 = (TH2F*)nphefile->Get("h_cc_nphe_final_s1");
	norm_nphe_s2 = (TH2F*)nphefile->Get("h_cc_nphe_final_s2");
	norm_nphe_s3 = (TH2F*)nphefile->Get("h_cc_nphe_final_s3");
	norm_nphe_s4 = (TH2F*)nphefile->Get("h_cc_nphe_final_s4");
	norm_nphe_s5 = (TH2F*)nphefile->Get("h_cc_nphe_final_s5");
	norm_nphe_s6 = (TH2F*)nphefile->Get("h_cc_nphe_final_s6");	

        //cout << " bin content1 = " << norm_nphe_s1->GetBinContent(100,100) << "\n";
	
      
	

      // cout << " bin content2 = " << avrg_nphe_sector1->GetBinContent(100,100) << "\n";
cout << "qqq1 \n";

ifstream input("phel_integr_fract_22jan17.txt");


if(input.is_open()){
i=0;
    while(!input.eof()){
          string line1,t_str, e_str,r_str,fr_str;
	   Int_t t,e,r;
	   Double_t fr;
           getline(input,line1); //read number
	   if (line1.length() != 0){ 
              t_str= line1.substr(0,line1.find(","));
            t = atof(t_str.c_str());
		   
	    e_str = line1.substr(t_str.length()+1,line1.substr(t_str.length()+1).find(","));
            e = atof(e_str.c_str());
	    	    
	    r_str = line1.substr(t_str.length()+e_str.length()+2, line1.substr(t_str.length()+e_str.length()+2).find(","));
            r = atof(r_str.c_str());
	    
	    fr_str = line1.substr(t_str.length()+e_str.length()+r_str.length()+3);
	    fr = atof(fr_str.c_str());
	    	    	    
	 //  cout << t<< "   " << e << "   " << r << "   " << fr <<" \n";
	    fract_integ[t][e][r] = fr;
	    i=i+1;
	    	    };
	    
    };
};

input.close();
//for(k=0; k<3; k++){
//for(i=0; i<6; i++){
//for(j=0; j<18; j++){

//cout << k << "," << i << "," <<j << "," << fract_integ[k][i][j]<< "\n";
//};
//};
//};

ph_cc_match = 1000;
  for (m=1; m<=n_files; m++) {

//sozdaem fail s imenem, vzyatim iz masiva strok  
  finp = new TFile(file[m-1].c_str()); 

  
  
 cout << "Processing file " << m << "\n"; 
 
 
 //cout << " bin content3 = " << norm_nphe_s1->GetBinContent(100,100) << "\n";
//berem derevo iz faila    
  TTree *t21 = (TTree*)finp->Get("t21");
  
   TBranch *br_indtype = t21->GetBranch("indtype");
    TBranch *br_npart = t21->GetBranch("npart");
    TBranch *br_p_el = t21->GetBranch("P_EL");
    TBranch *br_block = t21->GetBranch("block");
    TBranch *br_deltaQ = t21->GetBranch("deltaQ");
    TBranch *br_LiveTime = t21->GetBranch("LiveTime");
    TBranch *br_n_incl = t21->GetBranch("n_incl");
    TBranch *br_n_elast = t21->GetBranch("n_elast");
    TBranch *br_ph_EL = t21->GetBranch("ph_EL");
    TBranch *br_th_EL = t21->GetBranch("th_EL");
    TBranch *br_W = t21->GetBranch("W");
    TBranch *br_Q2 = t21->GetBranch("Q2");
    TBranch *br_nphe_el = t21->GetBranch("NpheCC_EL");
    TBranch *br_ectot_el = t21->GetBranch("ECtot_EL");
    TBranch *br_ecin_el = t21->GetBranch("ECin_EL");
    TBranch *br_ecout_el = t21->GetBranch("ECout_EL");
    TBranch *br_x_EL = t21->GetBranch("x_EL");
    TBranch *br_y_EL = t21->GetBranch("y_EL");
    TBranch *br_z_EL = t21->GetBranch("z_EL");
    TBranch *br_z_P = t21->GetBranch("dc_z_P");
    TBranch *br_z_PIp = t21->GetBranch("dc_z_PIp");
    TBranch *br_z_PIm = t21->GetBranch("dc_z_PIm");
    TBranch *br_PdHit_EL = t21->GetBranch("PdHit_EL");
    TBranch *br_PdHit_PIp = t21->GetBranch("PdHit_PIp");
    TBranch *br_PdHit_PIm = t21->GetBranch("PdHit_PIm");
    TBranch *br_PdHit_P = t21->GetBranch("PdHit_P");
    TBranch *br_sc_x = t21->GetBranch("sc_x");
    TBranch *br_sc_y = t21->GetBranch("sc_y");
    TBranch *br_sc_z = t21->GetBranch("sc_z");
    TBranch *br_pmt_hit = t21->GetBranch("pmt_hit");
    TBranch *br_segment = t21->GetBranch("segment");
    TBranch *br_theta_cc = t21->GetBranch("theta_cc");
    TBranch *br_ph_cc = t21->GetBranch("ph_cc");
    TBranch *br_sector = t21->GetBranch("sector");
    TBranch *br_n_PIp = t21->GetBranch("n_PIp");
    TBranch *br_n_PIm = t21->GetBranch("n_PIm");
    TBranch *br_n_P = t21->GetBranch("n_P");
    TBranch *br_p_pip = t21->GetBranch("P_PIp");
    TBranch *br_p_pim = t21->GetBranch("P_PIm");
    TBranch *br_p_p = t21->GetBranch("P_P");
    TBranch *br_th_PIp = t21->GetBranch("th_PIp");
    TBranch *br_th_PIm = t21->GetBranch("th_PIm");
    TBranch *br_th_P = t21->GetBranch("th_P");
    TBranch *br_ph_PIp = t21->GetBranch("ph_PIp");
    TBranch *br_ph_PIm = t21->GetBranch("ph_PIm");
    TBranch *br_ph_P = t21->GetBranch("ph_P");
    TBranch *br_beta_PIm = t21->GetBranch("beta_PIm");
    TBranch *br_beta_PIp = t21->GetBranch("beta_PIp");
    TBranch *br_beta_P = t21->GetBranch("beta_P");
    TBranch *br_nphe_pip = t21->GetBranch("NpheCC_PIp");
    TBranch *br_PIp_time = t21->GetBranch("PIp_time");
    TBranch *br_PIp_dist = t21->GetBranch("PIp_dist");
    TBranch *br_P_time = t21->GetBranch("P_time");
    TBranch *br_P_dist = t21->GetBranch("P_dist");
    TBranch *br_PIm_time = t21->GetBranch("PIm_time");
    TBranch *br_PIm_dist = t21->GetBranch("PIm_dist");
    TBranch *br_beta_P_time = t21->GetBranch("beta_P_time");
    TBranch *br_beta_PIp_time = t21->GetBranch("beta_PIp_time");
    TBranch *br_beta_PIm_time = t21->GetBranch("beta_PIm_time"); 
    
    TBranch *br_dc_x_EL = t21->GetBranch("dc_x_EL"); 
    TBranch *br_dc_y_EL = t21->GetBranch("dc_y_EL"); 
    
  Bool_t adc_cut_switch,tdc_cut_switch;
  
  Qdiff = 0.;
  Qcurr = 0.;
  Qprev = 0.;
  Qtotal = 0.;
  k = 0;
  block = 0;
  last_k = 0;
  nstart = 0;
  nstop = 0;
  n_incl = 0;
  n_elast = 0;
  

//tsikl po sobitiyam. ih odinakovoe kol-vo v kazdoi peremennoi  
  for (i=0; i<br_sector->GetEntries(); i++) { 
  
  Qprev = Qcurr;

br_segment->GetEntry(i);
br_pmt_hit->GetEntry(i);
 br_deltaQ->GetEntry(i);
  br_n_incl->GetEntry(i);
  br_n_elast->GetEntry(i);  
  br_LiveTime->GetEntry(i);
  br_block->GetEntry(i);  
  br_W->GetEntry(i);
  br_Q2->GetEntry(i); 
  br_npart->GetEntry(i);
  br_n_PIm->GetEntry(i);
  br_n_PIp->GetEntry(i);
  br_n_P->GetEntry(i);
  br_sector->GetEntry(i);
  br_LiveTime->GetEntry(i);
  br_block->GetEntry(i);  
  br_ph_cc->GetEntry(i);
  br_theta_cc->GetEntry(i);  
  br_p_el->GetEntry(i);  
  br_p_p->GetEntry(i); 
  br_p_pim->GetEntry(i);
  br_p_pip->GetEntry(i); 
  br_ectot_el->GetEntry(i); 
  br_ecin_el->GetEntry(i); 
  br_ecout_el->GetEntry(i);
  br_nphe_el->GetEntry(i);
  br_nphe_pip->GetEntry(i);
  br_th_EL->GetEntry(i);
  br_ph_EL->GetEntry(i); 
  br_th_P->GetEntry(i);
  br_ph_P->GetEntry(i);
  br_th_PIp->GetEntry(i);
  br_ph_PIp->GetEntry(i);
  br_th_PIm->GetEntry(i);
  br_ph_PIm->GetEntry(i);
  br_beta_PIp->GetEntry(i);
  br_beta_PIm->GetEntry(i);
  br_beta_P->GetEntry(i);
  br_PdHit_EL->GetEntry(i);
  br_PdHit_PIp->GetEntry(i);
  br_PdHit_PIm->GetEntry(i);
  br_PdHit_P->GetEntry(i);
  br_z_EL->GetEntry(i);
    br_z_P->GetEntry(i);
    br_z_PIp->GetEntry(i);
    br_z_PIm->GetEntry(i);
  br_W->GetEntry(i);
  br_Q2->GetEntry(i);
  br_indtype->GetEntry(i);
  
  
  br_PIp_time->GetEntry(i);
  br_PIm_time->GetEntry(i);
  br_P_time->GetEntry(i); 
  br_PIp_dist->GetEntry(i);
  br_PIm_dist->GetEntry(i);
  br_P_dist->GetEntry(i);
  br_beta_P_time->GetEntry(i);
  br_beta_PIp_time->GetEntry(i);
  br_beta_PIm_time->GetEntry(i);
  
  
  br_dc_x_EL->GetEntry(i);
  br_dc_y_EL->GetEntry(i);
  
  P_EL = br_p_el->GetLeaf("P_EL")->GetValue();
  P_P = br_p_p->GetLeaf("P_P")->GetValue();
  P_PIp = br_p_pip->GetLeaf("P_PIp")->GetValue();
  P_PIm = br_p_pim->GetLeaf("P_PIm")->GetValue();
  th_EL = br_th_EL->GetLeaf("th_EL")->GetValue();
  ph_EL = br_ph_EL->GetLeaf("ph_EL")->GetValue();
  th_P = br_th_P->GetLeaf("th_P")->GetValue();
  ph_P = br_ph_P->GetLeaf("ph_P")->GetValue();
  th_PIp = br_th_PIp->GetLeaf("th_PIp")->GetValue();
  ph_PIp = br_ph_PIp->GetLeaf("ph_PIp")->GetValue();
  npart = br_npart->GetLeaf("npart")->GetValue();
  n_P = br_n_P->GetLeaf("n_P")->GetValue();
  n_PIp = br_n_PIp->GetLeaf("n_PIp")->GetValue();
  Nphe_pip = br_nphe_pip->GetLeaf("NpheCC_PIp")->GetValue();
  PdHit_EL = br_PdHit_EL->GetLeaf("PdHit_EL")->GetValue();
  PdHit_PIp = br_PdHit_PIp->GetLeaf("PdHit_PIp")->GetValue();
  PdHit_PIm = br_PdHit_PIm->GetLeaf("PdHit_PIm")->GetValue();
  PdHit_P = br_PdHit_P->GetLeaf("PdHit_P")->GetValue();
//  beta_PIp = br_beta_PIp->GetLeaf("beta_PIp")->GetValue();
//  beta_PIm = br_beta_PIm->GetLeaf("beta_PIm")->GetValue();
//  beta_P = br_beta_P->GetLeaf("beta_P")->GetValue();
  n_PIm = br_n_PIm->GetLeaf("n_PIm")->GetValue();
  ph_PIm = br_ph_PIm->GetLeaf("ph_PIm")->GetValue();
  th_PIm = br_th_PIm->GetLeaf("th_PIm")->GetValue();
  z_EL = br_z_EL->GetLeaf("z_EL")->GetValue();
  z_P = br_z_P->GetLeaf("dc_z_P")->GetValue();
  z_PIp = br_z_PIp->GetLeaf("dc_z_PIp")->GetValue();
  z_PIm = br_z_PIm->GetLeaf("dc_z_PIm")->GetValue();
  W = br_W->GetLeaf("W")->GetValue();
  Q2 = br_Q2->GetLeaf("Q2")->GetValue();
  segment = br_segment->GetLeaf("segment")->GetValue();
  pmt_hit = br_pmt_hit->GetLeaf("pmt_hit")->GetValue();

  PIp_time = br_PIp_time->GetLeaf("PIp_time")->GetValue();
  PIm_time = br_PIm_time->GetLeaf("PIm_time")->GetValue();
  P_time = br_P_time->GetLeaf("P_time")->GetValue();
  PIp_dist = br_PIp_dist->GetLeaf("PIp_dist")->GetValue();
  PIm_dist = br_PIm_dist->GetLeaf("PIm_dist")->GetValue();
  P_dist = br_P_dist->GetLeaf("P_dist")->GetValue();
  indtype = br_indtype->GetLeaf("indtype")->GetValue();
  
  dc_x_el = br_dc_x_EL->GetLeaf("dc_x_EL")->GetValue();
   dc_y_el = br_dc_y_EL->GetLeaf("dc_y_EL")->GetValue();
  
//  cout << beta_P<<" "<<br_beta_P_time->GetLeaf("beta_P_time")->GetValue()<<" \n";

  beta_P = br_beta_P_time->GetLeaf("beta_P_time")->GetValue();
  beta_PIp = br_beta_PIp_time->GetLeaf("beta_PIp_time")->GetValue();
  beta_PIm = br_beta_PIm_time->GetLeaf("beta_PIm_time")->GetValue();
 
  if (br_ectot_el->GetLeaf("ECtot_EL")->GetValue() > (br_ecout_el->GetLeaf("ECout_EL")->GetValue()+br_ecin_el->GetLeaf("ECin_EL")->GetValue())) { 
          ECT  = br_ectot_el->GetLeaf("ECtot_EL")->GetValue();  
         } else {
          ECT  = (br_ecout_el->GetLeaf("ECout_EL")->GetValue()+br_ecin_el->GetLeaf("ECin_EL")->GetValue());    
          };
	  
	  //quatity check   
	     
	     block_curr = br_block->GetLeaf("block")->GetValue(); 

    if ((m_old != 0) && (m != m_old))block_total=block_total+block_prev;
    
    m_old=m;


   LiveTime=br_LiveTime->GetLeaf("LiveTime")->GetValue();
   
   inclusive = br_n_incl->GetLeaf("n_incl")->GetValue()/br_deltaQ->GetLeaf("deltaQ")->GetValue();
   elastic =  br_n_elast->GetLeaf("n_elast")->GetValue()/br_deltaQ->GetLeaf("deltaQ")->GetValue(); 

   if (block_curr != block_prev) {

   hist_ltime->Fill(block_total+block_curr,LiveTime);
   
   hist_ltime_1d->Fill(LiveTime);
   hist_n_incl->Fill(block_total+block_curr,inclusive);
hist_n_incl_1d->Fill(inclusive);  
 hist_n_elast->Fill(block_total+block_curr,elastic);
hist_n_elast_1d->Fill(elastic); 

 
if ((LiveTime > 0.8) && (LiveTime <0.9) && (inclusive > 85000) &&(inclusive < 100000) && (elastic > 24000) && (elastic < 28000)){
 
 Qfull = Qfull + br_deltaQ->GetLeaf("deltaQ")->GetValue();
 };
   
  
   
   block_prev = block_curr;
   };   
 

//----------------PROTON ENERGY LOSS Ye unfold-------------
delta_mom_p_ye = corrfunc.correct_energy_theta_pf(P_P, th_P);
P_P = P_P - delta_mom_p_ye;
//----------PROTON ENERGY LOSS----------------------------
delta_mom_p_skor = corrfunc.corr_pr_mom_skor(P_P, th_P);
//--------------------------------------------------------

//MOM CORR
/*
P_EL_new = corrfunc.correct_pel_e1_2039_2250_feb09(P_EL,th_EL,ph_EL);
th_EL_new = corrfunc.correct_thel_e1_2039_2250_feb09(P_EL,th_EL,ph_EL);
P_EL_old = P_EL;
P_EL = P_EL_new;
th_EL_old = th_EL;
th_EL = th_EL_new;



W_new = pow(float((E0+m_proton-P_EL)),2);
W_new = W_new -pow(P_EL*sin(th_EL*M_PI/180.)*cos(ph_EL*M_PI/180.),2);
W_new = W_new -pow(P_EL*sin(th_EL*M_PI/180.)*sin(ph_EL*M_PI/180.),2);
W_new = W_new -pow(float(E0 - P_EL*cos(th_EL*M_PI/180.)),2);
W_new = sqrt(W_new);
W = W_new;  

Q2_new = pow(float(E0-P_EL),2);
Q2_new = Q2_new -pow(P_EL*sin(th_EL*M_PI/180.)*cos(ph_EL*M_PI/180.),2);
Q2_new = Q2_new -pow(P_EL*sin(th_EL*M_PI/180.)*sin(ph_EL*M_PI/180.),2);
Q2_new = Q2_new -pow(float(E0 - P_EL*cos(th_EL*M_PI/180.)),2);
Q2 = -Q2_new;

*/   
 
sector =  int(br_sector->GetLeaf("sector")->GetValue());
nphe = br_nphe_el->GetLeaf("NpheCC_EL")->GetValue();
theta_cc = br_theta_cc->GetLeaf("theta_cc")->GetValue();
ph_cc = br_ph_cc->GetLeaf("ph_cc")->GetValue();  

  P4_EL.SetXYZT(0,0,2.039,2.039);
	P4_P.SetXYZT(0,0,0,m_proton);

//Missing mass pim
/*p_miss_1 = -P_EL*cos(ph_EL*M_PI/180.)*sin(th_EL*M_PI/180.)-P_P*cos(ph_P*M_PI/180.)*sin(th_P*M_PI/180.)-P_PIp*cos(ph_PIp*M_PI/180.)*sin(th_PIp*M_PI/180.);
p_miss_2 = -P_EL*sin(ph_EL*M_PI/180.)*sin(th_EL*M_PI/180.)-P_P*sin(ph_P*M_PI/180.)*sin(th_P*M_PI/180.)-P_PIp*sin(ph_PIp*M_PI/180.)*sin(th_PIp*M_PI/180.);
p_miss_3 = E0-P_EL*cos(th_EL*M_PI/180.)-P_P*cos(th_P*M_PI/180.)-P_PIp*cos(th_PIp*M_PI/180.);
p_miss_4 = E0+m_proton-P_EL-sqrt(m_proton*m_proton+P_P*P_P)-sqrt(m_pip*m_pip+P_PIp*P_PIp);

p_miss_sqr = p_miss_4*p_miss_4-p_miss_3*p_miss_3-p_miss_2*p_miss_2-p_miss_1*p_miss_1;*/

P4_ELP_reg.SetXYZT(P_EL*cos(ph_EL*M_PI/180.)*sin(th_EL*M_PI/180.),P_EL*sin(ph_EL*M_PI/180.)*sin(th_EL*M_PI/180.),P_EL*cos(th_EL*M_PI/180.),P_EL);

P4_PP_reg.SetXYZT(P_P*cos(ph_P*M_PI/180.)*sin(th_P*M_PI/180.),P_P*sin(ph_P*M_PI/180.)*sin(th_P*M_PI/180.),P_P*cos(th_P*M_PI/180.),sqrt(m_proton*m_proton+P_P*P_P));

P4_PIp_reg.SetXYZT(P_PIp*cos(ph_PIp*M_PI/180.)*sin(th_PIp*M_PI/180.),P_PIp*sin(ph_PIp*M_PI/180.)*sin(th_PIp*M_PI/180.),P_PIp*cos(th_PIp*M_PI/180.),sqrt(m_pip*m_pip+P_PIp*P_PIp));


P4_PIm_reg.SetXYZT(P_PIm*cos(ph_PIm*M_PI/180.)*sin(th_PIm*M_PI/180.),P_PIm*sin(ph_PIm*M_PI/180.)*sin(th_PIm*M_PI/180.),P_PIm*cos(th_PIm*M_PI/180.),sqrt(m_pip*m_pip+P_PIm*P_PIm));




//cout << th_ph_pr<<"\n";

P4_PIm_miss = P4_EL + P4_P - P4_ELP_reg - P4_PP_reg - P4_PIp_reg;
//P4_PIm_miss = P4_EL + P4_P - P4_ELP_reg - P4_PP_reg - P4_PIp_reg;
P4_PIp_miss = P4_EL + P4_P - P4_ELP_reg - P4_PP_reg - P4_PIm_reg;

P4_PIp_miss_d = P4_EL + P4_D - P4_ELP_reg - P4_PP_reg - P4_PIm_reg;
P4_PIm_miss_d = P4_EL + P4_D - P4_ELP_reg - P4_PP_reg - P4_PIp_reg;

P4_P_miss = P4_EL + P4_P - P4_ELP_reg - P4_PIp_reg - P4_PIm_reg;
P4_miss_0 = P4_EL + P4_P - P4_ELP_reg - P4_PP_reg - P4_PIm_reg -P4_PIp_reg;

P4_miss_0_d = P4_EL + P4_D - P4_ELP_reg - P4_PP_reg - P4_PIm_reg -P4_PIp_reg;


//cout<< P4_P[0]<<" "<< P4_P[1]<<" "<< P4_P[2]<<" "<< P4_P[3]<<" \n";

//P4_miss_N = 

//W_old = W;
//W = sqrt((P4_PIp_reg+P4_PIm_reg+P4_PP_reg).Mag2());
//W=sqrt((P4_P+P4_EL-P4_ELP_reg).Mag2());





th_PIm_miss = (180./M_PI)*acos(P4_PIm_miss[2]/sqrt(P4_PIm_miss[0]*P4_PIm_miss[0]+P4_PIm_miss[1]*P4_PIm_miss[1]+P4_PIm_miss[2]*P4_PIm_miss[2]));

/*
if(P4_PIm_miss[0] == 0.) {
if(P4_PIm_miss[1] > 0.) ph_PIm_miss = 90.;
if(P4_PIm_miss[1] < 0.) ph_PIm_miss = 270.;

};


if ((P4_PIm_miss[1]>0)&&(P4_PIm_miss[0]>0)&&(P4_PIm_miss[0] != 0.)) {
ph_PIm_miss = (180./M_PI)*atan(P4_PIm_miss[1]/P4_PIm_miss[0]);
};


if ((P4_PIm_miss[1]>0)&&(P4_PIm_miss[0]<0)&&(P4_PIm_miss[0] != 0.)) {
ph_PIm_miss = 180. - (180./M_PI)*atan(abs(P4_PIm_miss[1]/P4_PIm_miss[0]));
};

if ((P4_PIm_miss[1]<0)&&(P4_PIm_miss[0]<0)&&(P4_PIm_miss[0] != 0.)) {
ph_PIm_miss = 180. + (180./M_PI)*atan(abs(P4_PIm_miss[1]/P4_PIm_miss[0]));
};

if ((P4_PIm_miss[1]<0)&&(P4_PIm_miss[0]>0)&&(P4_PIm_miss[0] != 0.)) {
ph_PIm_miss = 360. - (180./M_PI)*atan(abs(P4_PIm_miss[1]/P4_PIm_miss[0]));
};*/


if (P4_PIm_miss[0] != 0.) {
ph_PIm_miss = (180./M_PI)*atan(P4_PIm_miss[1]/P4_PIm_miss[0]);
}
 else {
if(P4_PIm_miss[1] > 0.) ph_PIm_miss = 90.;
if(P4_PIm_miss[1] < 0.) ph_PIm_miss = 270.;
};
				   

if ((P4_PIm_miss[0] < 0.) && (P4_PIm_miss[1] > 0)) ph_PIm_miss = ph_PIm_miss+180.;
if (( P4_PIm_miss[0]< 0.) && ( P4_PIm_miss[1]< 0)) ph_PIm_miss = ph_PIm_miss+180.;
if ((P4_PIm_miss[0] > 0.) && ( P4_PIm_miss[1]< 0)) ph_PIm_miss = ph_PIm_miss+360.;

//E_miss = E0-m_proton-P_EL-sqrt(m_proton^2+P_P^2)-sqrt(m_pip^2+P_PIp^2)-sqrt(m_pip^2+P_PIm^2);

//cout << "theta_cc = " << theta_cc << " ph_cc = " << ph_cc << " bin content = " << norm_nphe_s1->GetBinContent(100,100) << "\n";


//if (PdHit_PIp==46) cout << "1 "<< beta_PIp << "\n";
//beta_func_data();
//if (PdHit_PIp==46) cout << "2 "<< beta_PIp << "\n";
//cout << n_PIp<<"\n";

bool_el_id_data=particle_ID_data.Electron_cuts_data(); 




selection = false;
selection_pim_miss = false;
selection_0_miss = false;

if (bool_el_id_data) {
//vertex difference cut
if ((abs(z_EL - z_P)<5.)&&(abs(z_EL - z_PIp)<5.)&&(abs(z_P - z_PIp)<5.)){

//cout << z_EL-z_P<<" "<< abs(z_EL-z_P)<<"\n";

//if ((bool_proton_id_data)&&(bool_pip_id_data)) W_2pi_selection->Fill(W,Q2,1.);
h_dc_y_vs_x_el->Fill(dc_x_el,dc_y_el,1.);



beta_func_data(beta_PIp);

if (particle_ID_data.Proton_cuts_data()) data_hist();

if (particle_ID_data.PIp_cuts_data()) h_z_corr1_data->Fill(z_PIp);
if ((particle_ID_data.PIp_cuts_data())&&(particle_ID_data.PIm_cuts_data())) h_z_corr2_data->Fill(z_PIm-z_PIp);



bool_proton_id_data=particle_ID_data.Proton_cuts_data();
bool_pip_id_data=particle_ID_data.PIp_cuts_data();
bool_pim_id_data=particle_ID_data.PIm_cuts_data();

if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)&&(npart >=3)){
if ((bool_proton_id_data)&&(bool_pim_id_data)&&(!bool_pip_id_data)){

hist_PIp_miss_d_bef-> Fill(P4_PIp_miss_d.Mag2(),1.);
if (sqrt((P4_PIp_miss+P4_PIm_reg)*(P4_PIp_miss+P4_PIm_reg)) > 2*m_pip - 0.05) {
if (sqrt((P4_PIp_miss+P4_PIm_reg)*(P4_PIp_miss+P4_PIm_reg)) < W - m_proton + 0.05) {
if (sqrt((P4_PIp_miss+P4_PP_reg)*(P4_PIp_miss+P4_PP_reg)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIp_miss+P4_PP_reg)*(P4_PIp_miss+P4_PP_reg)) < W - m_pip + 0.05) {
if (sqrt((P4_PIm_reg+P4_PP_reg)*(P4_PIm_reg+P4_PP_reg)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIm_reg+P4_PP_reg)*(P4_PIm_reg+P4_PP_reg)) < W - m_pip + 0.05) {

h_mm_pip_vs_npart->Fill(npart,P4_PIp_miss.Mag2(),1.);
//hist_PIp_miss_en->Fill(P4_PIp_miss[3],1.);

if (P4_PIp_miss[3] > 0.2) {



//if (P4_PIp_miss.Mag2() > m_pip*m_pip - 4*0.0302){
hist_PIp_miss_en->Fill(P4_PIp_miss[3],1.);
hist_PIp_miss-> Fill(P4_PIp_miss.Mag2(),1.);
hist_PIp_miss_d-> Fill(P4_PIp_miss_d.Mag2(),1.);
//hist_PIp_miss_cor-> Fill(P4_PIp_cor.Mag2(),1.);

selection = false;
//P4_PIp_reg = P4_PIp_miss;
//};
};
};
};
};
};
};
};


};
};


if ((bool_pip_id_data)&&(bool_pim_id_data)&&(!bool_proton_id_data)){
if (sqrt((P4_PIp_reg+P4_PIm_reg)*(P4_PIp_reg+P4_PIm_reg)) > 2*m_pip - 0.05) {
if (sqrt((P4_PIp_reg+P4_PIm_reg)*(P4_PIp_reg+P4_PIm_reg)) < W - m_proton + 0.05) {
if (sqrt((P4_PIp_reg+P4_P_miss)*(P4_PIp_reg+P4_P_miss)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIp_reg+P4_P_miss)*(P4_PIp_reg+P4_P_miss)) < W - m_pip + 0.05) {
if (sqrt((P4_PIm_reg+P4_P_miss)*(P4_PIm_reg+P4_P_miss)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIm_reg+P4_P_miss)*(P4_PIm_reg+P4_P_miss)) < W - m_pip + 0.05) {
if (P4_P_miss[3] > m_proton  - 0.05) {
hist_P_miss-> Fill(P4_P_miss.Mag2(),1.);
hist_P_miss_en->Fill(P4_P_miss[3],1.);
if (P4_P_miss.Mag2() > m_proton*m_proton - 4*0.103 ){
selection = false;
//P4_PP_reg = P4_P_miss;
};
};
};
};
};
};
};
};

};





if ((W > 1.3)&&(W < 1.8125)&&(Q2 > 0.4)&&(Q2 < 1.)&&(npart>=4)){
//if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)){
if ((bool_proton_id_data)&&(bool_pip_id_data)&&(bool_pim_id_data)){

//excl top vertex difference cut 
if ((abs(z_EL - z_PIm)<5.)&&(abs(z_P - z_PIm)<5.)&&(abs(z_PIp - z_PIm)<5.)){


if (sqrt((P4_PIp_reg+P4_PIm_reg)*(P4_PIp_reg+P4_PIm_reg)) > 2*m_pip - 0.01) {
if (sqrt((P4_PIp_reg+P4_PIm_reg)*(P4_PIp_reg+P4_PIm_reg)) < W - m_proton + 0.01) {
if (sqrt((P4_PIp_reg+P4_PP_reg)*(P4_PIp_reg+P4_PP_reg)) > m_pip + m_proton  - 0.01) {
if (sqrt((P4_PIp_reg+P4_PP_reg)*(P4_PIp_reg+P4_PP_reg)) < W - m_pip + 0.01) {
if (sqrt((P4_PIm_reg+P4_PP_reg)*(P4_PIm_reg+P4_PP_reg)) > m_pip + m_proton  - 0.01) {
if (sqrt((P4_PIm_reg+P4_PP_reg)*(P4_PIm_reg+P4_PP_reg)) < W - m_pip + 0.01) {

//P4_Pini_ferm_dat.SetXYZT(P4_miss_0[0],P4_miss_0[1],P4_miss_0[2],);


P4_inprot_miss =-(P4_EL - P4_ELP_reg - P4_PP_reg - P4_PIm_reg -P4_PIp_reg);
h_mm_0_vs_npart->Fill(npart,P4_miss_0.Mag2(),1.);

th_ph_pim = 180./M_PI*acos((P4_PIm_reg.Vect().Dot((P4_EL-P4_ELP_reg).Vect()))/(((P4_EL-P4_ELP_reg).Vect()).Mag())/((P4_PIm_reg.Vect()).Mag()));
th_ph_pip = 180./M_PI*acos((P4_PIp_reg.Vect().Dot((P4_EL-P4_ELP_reg).Vect()))/(((P4_EL-P4_ELP_reg).Vect()).Mag())/((P4_PIp_reg.Vect()).Mag()));
th_ph_pr = 180./M_PI*acos((P4_PP_reg.Vect().Dot((P4_EL-P4_ELP_reg).Vect()))/(((P4_EL-P4_ELP_reg).Vect()).Mag())/((P4_PP_reg.Vect()).Mag()));


//cut on missing energy

//if ((P4_miss_0[3] >  -0.05)) {

	h_z_P->Fill(z_P,fract_integ[pmt_hit+1][sector-1][segment]);
	h_z_PIp->Fill(z_PIp,fract_integ[pmt_hit+1][sector-1][segment]);
	h_z_PIm->Fill(z_PIm,fract_integ[pmt_hit+1][sector-1][segment]);

	if ((W>1.3)&&(W<1.8)&&(Q2>0.45)&&(Q2<1.)){
	
	if (((P4_miss_0.Vect()).Mag()<0.2)) {
	
	if ((P4_miss_0.Mag2()>-0.02)&&(P4_miss_0.Mag2()<0.001)){
	h_pim_mis_all_reg[int((W-1.3)/0.1)]-> Fill(P4_PIm_miss.Mag2(),1.);
	hist_miss_en_0->Fill(P4_miss_0[3],1.);
	
	};
	h_0_mis_all_reg[int((W-1.3)/0.1)]-> Fill(P4_miss_0.Mag2(),1.);
	
	};
	h_mom_all_reg[int((W-1.3)/0.1)]-> Fill((P4_miss_0.Vect()).Mag(),1.);
	
	};


/*
if ((W>1.3)&&(W<1.8)&&(Q2>0.4)&&(Q2<0.6)) {
h_pim_mis_fermi_nocut_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.1)]->Fill(P4_PIm_miss.Mag2(),fract_integ[pmt_hit+1][sector-1][segment]);
h_mis_mom_fermi_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.1)]->Fill((P4_miss_0.Vect()).Mag(),fract_integ[pmt_hit+1][sector-1][segment]);
};
if ((W>1.3)&&(W<1.775)&&(Q2>0.6)&&(Q2<0.7)) {
h_pim_mis_fermi_nocut_2[int((Q2-0.6)/0.05)][int((W-1.3)/0.1)]->Fill(P4_PIm_miss.Mag2(),fract_integ[pmt_hit+1][sector-1][segment]);
h_mis_mom_fermi_2[int((Q2-0.6)/0.05)][int((W-1.3)/0.1)]->Fill((P4_miss_0.Vect()).Mag(),fract_integ[pmt_hit+1][sector-1][segment]);
};
if ((W>1.3)&&(W<1.7)&&(Q2>0.7)&&(Q2<0.8)) {
h_pim_mis_fermi_nocut_3[int((Q2-0.7)/0.05)][int((W-1.3)/0.1)]->Fill(P4_PIm_miss.Mag2(),fract_integ[pmt_hit+1][sector-1][segment]);
h_mis_mom_fermi_3[int((Q2-0.7)/0.05)][int((W-1.3)/0.1)]->Fill((P4_miss_0.Vect()).Mag(),fract_integ[pmt_hit+1][sector-1][segment]);
};
if ((W>1.3)&&(W<1.65)&&(Q2>0.8)&&(Q2<0.85)) {
h_pim_mis_fermi_nocut_4[int((W-1.3)/0.1)]->Fill(P4_PIm_miss.Mag2(),fract_integ[pmt_hit+1][sector-1][segment]);
h_mis_mom_fermi_4[int((W-1.3)/0.1)]->Fill((P4_miss_0.Vect()).Mag(),fract_integ[pmt_hit+1][sector-1][segment]);
};
if ((W>1.3)&&(W<1.6)&&(Q2>0.85)&&(Q2<0.9)) {
h_pim_mis_fermi_nocut_5[int((W-1.3)/0.1)]->Fill(P4_PIm_miss.Mag2(),fract_integ[pmt_hit+1][sector-1][segment]);
h_mis_mom_fermi_5[int((W-1.3)/0.1)]->Fill((P4_miss_0.Vect()).Mag(),fract_integ[pmt_hit+1][sector-1][segment]);
};
if ((W>1.3)&&(W<1.575)&&(Q2>0.9)&&(Q2<0.95)) {
h_pim_mis_fermi_nocut_6[int((W-1.3)/0.1)]->Fill(P4_PIm_miss.Mag2(),fract_integ[pmt_hit+1][sector-1][segment]);
h_mis_mom_fermi_6[int((W-1.3)/0.1)]->Fill((P4_miss_0.Vect()).Mag(),fract_integ[pmt_hit+1][sector-1][segment]);
};
if ((W>1.3)&&(W<1.55)&&(Q2>0.95)&&(Q2<1.0)) {
h_pim_mis_fermi_nocut_7[int((W-1.3)/0.1)]->Fill(P4_PIm_miss.Mag2(),fract_integ[pmt_hit+1][sector-1][segment]);
h_mis_mom_fermi_7[int((W-1.3)/0.1)]->Fill((P4_miss_0.Vect()).Mag(),fract_integ[pmt_hit+1][sector-1][segment]);
};
*/

	h_miss_mom_0_nocut->Fill((P4_miss_0.Vect()).Mag(),1.);

	if ((P4_PIm_miss.Mag2()>0)&&(P4_PIm_miss.Mag2()<0.05)){
	h_miss_mom_0_cut_onpim->Fill((P4_miss_0.Vect()).Mag(),1.);
	};

	if ((W>1.3)&&(W<1.8)&&(Q2>0.4)&&(Q2<0.6)){
	h_inv_NP[int((W-1.3)/0.1)]->Fill(sqrt((P4_PP_reg+P4_miss_0_d ).Mag2()),1.);
	h_inv_NPIp[int((W-1.3)/0.1)]->Fill(sqrt((P4_PIp_reg+P4_miss_0_d ).Mag2()),1.);
	h_inv_NPIm[int((W-1.3)/0.1)]->Fill(sqrt((P4_PIm_reg+P4_miss_0_d ).Mag2()),1.);
	};

	hist_PIm_miss_all_reg_1-> Fill(P4_PIm_miss.Mag2(),1.);

//cut on missing momentum
if ((P4_miss_0.Vect()).Mag() <0.2){ 

	h_miss_mom_0->Fill((P4_miss_0.Vect()).Mag(),1.);

//cut on missing mass of pim
if ((P4_PIm_miss.Mag2()>-0.15)&&(P4_PIm_miss.Mag2()<0.15)){

	h_miss_mass_0-> Fill(P4_miss_0.Mag2(),1.);
	hist_PIm_miss_all_reg_2-> Fill(P4_PIm_miss.Mag2(),1.);

	h_miss_mom_0_cut_on0->Fill((P4_miss_0.Vect()).Mag(),1.);

	hist_w_hadr_all_reg->Fill(W,1.);
	hist_w_el_all_reg->Fill(W,1.);


	hist_PIm_miss_en->Fill(P4_PIm_miss[3],1.);

//cut on missing mass of 0
if ((P4_miss_0.Mag2()>-0.02)&&(P4_miss_0.Mag2()<0.001)){

selection_0_miss = true;
P4_PIm_reg = P4_PIm_miss;


//};
};
};
};
};
};
};
};
};
};
};
};
};


//if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)&&(npart>=3)){
if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)&&(npart>=3)){
if ((bool_pip_id_data)&&(bool_proton_id_data)&&(n_PIm == 0)){


if (sqrt((P4_PIp_reg+P4_PIm_miss)*(P4_PIp_reg+P4_PIm_miss)) > 2*m_pip ) {
if (sqrt((P4_PIp_reg+P4_PIm_miss)*(P4_PIp_reg+P4_PIm_miss)) < W - m_proton) {
if (sqrt((P4_PIp_reg+P4_PP_reg)*(P4_PIp_reg+P4_PP_reg)) > m_pip + m_proton ) {
if (sqrt((P4_PIp_reg+P4_PP_reg)*(P4_PIp_reg+P4_PP_reg)) < W - m_pip ) {
if (sqrt((P4_PIm_miss+P4_PP_reg)*(P4_PIm_miss+P4_PP_reg)) > m_pip + m_proton) {
if (sqrt((P4_PIm_miss+P4_PP_reg)*(P4_PIm_miss+P4_PP_reg)) < W - m_pip) {

h_mm_pim_vs_npart->Fill(npart,P4_PIm_miss.Mag2(),1.);
//hist_PIm_miss_en->Fill(P4_PIm_miss[3],1.);

th_ph_pip = 180./M_PI*acos((P4_PIp_reg.Vect().Dot((P4_EL-P4_ELP_reg).Vect()))/(((P4_EL-P4_ELP_reg).Vect()).Mag())/((P4_PIp_reg.Vect()).Mag()));
th_ph_pr = 180./M_PI*acos((P4_PP_reg.Vect().Dot((P4_EL-P4_ELP_reg).Vect()))/(((P4_EL-P4_ELP_reg).Vect()).Mag())/((P4_PP_reg.Vect()).Mag()));

if((P4_PIm_miss[3] > m_pip )) {



hist_PIm_miss-> Fill(P4_PIm_miss.Mag2(),fract_integ[pmt_hit+1][sector-1][segment]);

	if ((W>1.3)&&(W<1.8)&&(Q2>0.45)&&(Q2<0.5)){
		
	h_pim_mis_main_top[int((W-1.3)/0.1)]-> Fill(P4_PIm_miss.Mag2(),1.);
		
	};

//if ((P_P>0.25)&&(P_P<0.5)&&(P4_PIm_miss.Mag()>0))||((P_P>0.5)&&(P_P<0.75)&&)


//if ((th_PIp>10.)&&(th_PIp<130.)&&(th_P>10.)&&(th_P<50.)&&(P_P > 0.25)&&(P_P < 1.5)){
//if (((W>1.3)&&(W<1.5))||((W>1.5)&&(W<1.825)&&(ARR_MMC_PIM_MISS[int((P_P-0.25)/0.25)][int((th_PIp-10.)/15.)][int((th_P-10.)/5.)]==1))){
/*
if ((W>1.3)&&(W<1.825)&&(Q2>0.4)&&(Q2<0.6)) {
h_pim_mis_fermi_nocut_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(sqrt(abs(P4_PIm_miss.Mag2())),fract_integ[pmt_hit+1][sector-1][segment]);
};

if ((W>1.3)&&(W<1.775)&&(Q2>0.6)&&(Q2<0.7)) {
h_pim_mis_fermi_nocut_2[int((Q2-0.6)/0.05)][int((W-1.3)/0.025)]->Fill(sqrt(abs(P4_PIm_miss.Mag2())),fract_integ[pmt_hit+1][sector-1][segment]);
};

if ((W>1.3)&&(W<1.725)&&(Q2>0.7)&&(Q2<0.8)) {
h_pim_mis_fermi_nocut_3[int((Q2-0.7)/0.05)][int((W-1.3)/0.025)]->Fill(sqrt(abs(P4_PIm_miss.Mag2())),fract_integ[pmt_hit+1][sector-1][segment]);
};

if ((W>1.3)&&(W<1.65)&&(Q2>0.8)&&(Q2<0.85)) {
h_pim_mis_fermi_nocut_4[int((W-1.3)/0.025)]->Fill(sqrt(abs(P4_PIm_miss.Mag2())),fract_integ[pmt_hit+1][sector-1][segment]);
};

if ((W>1.3)&&(W<1.625)&&(Q2>0.85)&&(Q2<0.9)) {
h_pim_mis_fermi_nocut_5[int((W-1.3)/0.025)]->Fill(sqrt(abs(P4_PIm_miss.Mag2())),fract_integ[pmt_hit+1][sector-1][segment]);
};

if ((W>1.3)&&(W<1.6)&&(Q2>0.9)&&(Q2<0.95)) {
h_pim_mis_fermi_nocut_6[int((W-1.3)/0.025)]->Fill(sqrt(abs(P4_PIm_miss.Mag2())),fract_integ[pmt_hit+1][sector-1][segment]);
};

if ((W>1.3)&&(W<1.55)&&(Q2>0.95)&&(Q2<1.0)) {
h_pim_mis_fermi_nocut_7[int((W-1.3)/0.025)]->Fill(sqrt(abs(P4_PIm_miss.Mag2())),fract_integ[pmt_hit+1][sector-1][segment]);
};*/

//};
//};
//};


//};
//if ((P4_PIm_miss.Mag2()> -0.0236)&&(P4_PIm_miss.Mag2()< 0.063)){
if (sqrt(abs(P4_PIm_miss.Mag2())) < MMcut_pim_miss[int((W-1.3)/0.025)]){
W_2pi_selection->Fill(W,Q2,1.);
selection_pim_miss = true;
P4_PIm_reg = P4_PIm_miss;
};
};

};
};
};
};
};
};
};
//};
//};
//};

};

};//konets ifa vertex difference cut 
};//konets ifa el cutov




if ((selection_pim_miss)||(selection_0_miss)) {



//data_hist();

if((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)){
rot_boost_cmsyst();
Var_1[0] = inv_m_pip_p;
Var_1[1] = inv_m_pip_pim;
Var_1[2] = theta_P_cm;
Var_1[3] = phi_P_cm;
Var_1[4] = alpha_PIpPIm_pipf;

Var_2[0] = inv_m_pip_p;
Var_2[1] = inv_m_pip_pim;
Var_2[2] = theta_PIm_cm;
Var_2[3] = phi_PIm_cm;
Var_2[4] = alpha_PPIp_piPIm;

Var_3[0] = inv_m_pim_p;
Var_3[1] = inv_m_pip_pim;
Var_3[2] = theta_PIp_cm;
Var_3[3] = phi_PIp_cm;
Var_3[4] = alpha_PPIm_piPIp;


if((selection_pim_miss)){


h_5dim_pim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_1,fract_integ[pmt_hit+1][sector-1][segment]);
h_5dim_pim_2[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_2,fract_integ[pmt_hit+1][sector-1][segment]);
h_5dim_pim_3[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_3,fract_integ[pmt_hit+1][sector-1][segment]); 
};

if ((selection_0_miss)) {

h_5dim_excl_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_1,fract_integ[pmt_hit+1][sector-1][segment]);
h_5dim_excl_2[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_2,fract_integ[pmt_hit+1][sector-1][segment]);
h_5dim_excl_3[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_3,fract_integ[pmt_hit+1][sector-1][segment]); 

};
};



 


};//end selection??
//};



    };// konets tsicla po sobitiyam (po i)
  
  
  t21->Delete();

 finp->Close();  
 
     };//konets tsicla po failam (po m)
    
    
 /////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////
 /////////////////////////////////////////////////////////////// 
 
 block_curr = 0;
 block_prev = 0;
 m_old=0;
  
 for (m=1; m<=n_files_empty; m++) {
  
  finp = new TFile(file_empty[m-1].c_str()); 
  
  
  
 cout << "Processing file with empty target " << m << "\n"; 
  
  
  // cout << " bin content3 = " << avrg_nphe_sector1->GetBinContent(100,100) << "\n";
  
  TTree *t21 = (TTree*)finp->Get("t21");
  
   TBranch *br_indtype = t21->GetBranch("indtype");
    TBranch *br_npart = t21->GetBranch("npart");
    TBranch *br_p_el = t21->GetBranch("P_EL");
    TBranch *br_block = t21->GetBranch("block");
    TBranch *br_deltaQ = t21->GetBranch("deltaQ");
    TBranch *br_LiveTime = t21->GetBranch("LiveTime");
    TBranch *br_n_incl = t21->GetBranch("n_incl");
    TBranch *br_n_elast = t21->GetBranch("n_elast");
    TBranch *br_ph_EL = t21->GetBranch("ph_EL");
    TBranch *br_th_EL = t21->GetBranch("th_EL");
    TBranch *br_W = t21->GetBranch("W");
    TBranch *br_Q2 = t21->GetBranch("Q2");
    TBranch *br_nphe_el = t21->GetBranch("NpheCC_EL");
    TBranch *br_ectot_el = t21->GetBranch("ECtot_EL");
    TBranch *br_ecin_el = t21->GetBranch("ECin_EL");
    TBranch *br_ecout_el = t21->GetBranch("ECout_EL");
    TBranch *br_x_EL = t21->GetBranch("x_EL");
    TBranch *br_y_EL = t21->GetBranch("y_EL");
    TBranch *br_z_EL = t21->GetBranch("z_EL");
    TBranch *br_z_P = t21->GetBranch("z_P");
    TBranch *br_z_PIp = t21->GetBranch("z_PIp");
    TBranch *br_z_PIm = t21->GetBranch("z_PIm");
    TBranch *br_PdHit_EL = t21->GetBranch("PdHit_EL");
    TBranch *br_PdHit_PIp = t21->GetBranch("PdHit_PIp");
    TBranch *br_PdHit_PIm = t21->GetBranch("PdHit_PIm");
    TBranch *br_PdHit_P = t21->GetBranch("PdHit_P");
    TBranch *br_sc_x = t21->GetBranch("sc_x");
    TBranch *br_sc_y = t21->GetBranch("sc_y");
    TBranch *br_sc_z = t21->GetBranch("sc_z");
    TBranch *br_pmt_hit = t21->GetBranch("pmt_hit");
    TBranch *br_segment = t21->GetBranch("segment");
    TBranch *br_theta_cc = t21->GetBranch("theta_cc");
    TBranch *br_ph_cc = t21->GetBranch("ph_cc");
    TBranch *br_sector = t21->GetBranch("sector");
    TBranch *br_n_PIp = t21->GetBranch("n_PIp");
    TBranch *br_n_PIm = t21->GetBranch("n_PIm");
    TBranch *br_n_P = t21->GetBranch("n_P");
    TBranch *br_p_pip = t21->GetBranch("P_PIp");
    TBranch *br_p_pim = t21->GetBranch("P_PIm");
    TBranch *br_p_p = t21->GetBranch("P_P");
    TBranch *br_th_PIp = t21->GetBranch("th_PIp");
    TBranch *br_th_PIm = t21->GetBranch("th_PIm");
    TBranch *br_th_P = t21->GetBranch("th_P");
    TBranch *br_ph_PIp = t21->GetBranch("ph_PIp");
    TBranch *br_ph_PIm = t21->GetBranch("ph_PIm");
    TBranch *br_ph_P = t21->GetBranch("ph_P");
    TBranch *br_beta_PIm = t21->GetBranch("beta_PIm");
    TBranch *br_beta_PIp = t21->GetBranch("beta_PIp");
    TBranch *br_beta_P = t21->GetBranch("beta_P");
    TBranch *br_nphe_pip = t21->GetBranch("NpheCC_PIp");
    TBranch *br_PIp_time = t21->GetBranch("PIp_time");
    TBranch *br_PIp_dist = t21->GetBranch("PIp_dist");
    TBranch *br_P_time = t21->GetBranch("P_time");
    TBranch *br_P_dist = t21->GetBranch("P_dist");
    TBranch *br_PIm_time = t21->GetBranch("PIm_time");
    TBranch *br_PIm_dist = t21->GetBranch("PIm_dist");  
    
   
   

  Bool_t adc_cut_switch,tdc_cut_switch;
  
  Qdiff = 0.;
  Qcurr = 0.;
  Qprev = 0.;
  Qtotal = 0.;
  k = 0;
  block = 0;
  last_k = 0;
  nstart = 0;
  nstop = 0;
  n_incl = 0;
  n_elast = 0;
  
  
  for (i=0; i<br_sector->GetEntries(); i++) { 
  
  Qprev = Qcurr;

  br_segment->GetEntry(i);
  br_pmt_hit->GetEntry(i);
  br_deltaQ->GetEntry(i);
  br_n_incl->GetEntry(i);
  br_n_elast->GetEntry(i);  
  br_LiveTime->GetEntry(i);
  br_block->GetEntry(i);  
  br_W->GetEntry(i);
  br_Q2->GetEntry(i); 
  br_npart->GetEntry(i);
  br_n_PIm->GetEntry(i);
  br_n_PIp->GetEntry(i);
  br_n_P->GetEntry(i);
  br_sector->GetEntry(i);
  br_LiveTime->GetEntry(i);
  br_block->GetEntry(i);  
  br_ph_cc->GetEntry(i);
  br_theta_cc->GetEntry(i);  
  br_p_el->GetEntry(i);  
  br_p_p->GetEntry(i); 
  br_p_pim->GetEntry(i);
  br_p_pip->GetEntry(i); 
  br_ectot_el->GetEntry(i); 
  br_ecin_el->GetEntry(i); 
  br_ecout_el->GetEntry(i);
  br_nphe_el->GetEntry(i);
  br_nphe_pip->GetEntry(i);
  br_th_EL->GetEntry(i);
  br_ph_EL->GetEntry(i); 
  br_th_P->GetEntry(i);
  br_ph_P->GetEntry(i);
  br_th_PIp->GetEntry(i);
  br_ph_PIp->GetEntry(i);
  br_th_PIm->GetEntry(i);
  br_ph_PIm->GetEntry(i);
  br_beta_PIp->GetEntry(i);
  br_beta_PIm->GetEntry(i);
  br_beta_P->GetEntry(i);
  br_PdHit_EL->GetEntry(i);
  br_PdHit_PIp->GetEntry(i);
  br_PdHit_PIm->GetEntry(i);
  br_PdHit_P->GetEntry(i);
  br_z_EL->GetEntry(i);
   br_z_P->GetEntry(i);
  br_W->GetEntry(i);
  br_Q2->GetEntry(i);
  br_indtype->GetEntry(i);
  br_PIp_time->GetEntry(i);
  br_PIm_time->GetEntry(i);
  br_P_time->GetEntry(i); 
  br_PIp_dist->GetEntry(i);
  br_PIm_dist->GetEntry(i);
  br_P_dist->GetEntry(i);
  
  
  
  P_EL = br_p_el->GetLeaf("P_EL")->GetValue();
  P_P = br_p_p->GetLeaf("P_P")->GetValue();
  P_PIp = br_p_pip->GetLeaf("P_PIp")->GetValue();
  P_PIm = br_p_pim->GetLeaf("P_PIm")->GetValue();
  th_EL = br_th_EL->GetLeaf("th_EL")->GetValue();
  ph_EL = br_ph_EL->GetLeaf("ph_EL")->GetValue();
  th_P = br_th_P->GetLeaf("th_P")->GetValue();
  ph_P = br_ph_P->GetLeaf("ph_P")->GetValue();
  th_PIp = br_th_PIp->GetLeaf("th_PIp")->GetValue();
  ph_PIp = br_ph_PIp->GetLeaf("ph_PIp")->GetValue();
  npart = br_npart->GetLeaf("npart")->GetValue();
  n_P = br_n_P->GetLeaf("n_P")->GetValue();
  n_PIp = br_n_PIp->GetLeaf("n_PIp")->GetValue();
  Nphe_pip = br_nphe_pip->GetLeaf("NpheCC_PIp")->GetValue();
  PdHit_EL = br_PdHit_EL->GetLeaf("PdHit_EL")->GetValue();
  PdHit_PIp = br_PdHit_PIp->GetLeaf("PdHit_PIp")->GetValue();
  PdHit_PIm = br_PdHit_PIm->GetLeaf("PdHit_PIm")->GetValue();
  PdHit_P = br_PdHit_P->GetLeaf("PdHit_P")->GetValue();
  beta_PIp = br_beta_PIp->GetLeaf("beta_PIp")->GetValue();
  beta_PIm = br_beta_PIm->GetLeaf("beta_PIm")->GetValue();
  beta_P = br_beta_P->GetLeaf("beta_P")->GetValue();
  n_PIm = br_n_PIm->GetLeaf("n_PIm")->GetValue();
  ph_PIm = br_ph_PIm->GetLeaf("ph_PIm")->GetValue();
  th_PIm = br_th_PIm->GetLeaf("th_PIm")->GetValue();
  z_EL = br_z_EL->GetLeaf("z_EL")->GetValue();
   z_P = br_z_P->GetLeaf("z_P")->GetValue();
   z_PIp = br_z_PIp->GetLeaf("z_PIp")->GetValue();
   z_PIm = br_z_PIm->GetLeaf("z_PIm")->GetValue(); 
  W = br_W->GetLeaf("W")->GetValue();
  Q2 = br_Q2->GetLeaf("Q2")->GetValue();
  segment = br_segment->GetLeaf("segment")->GetValue();
  pmt_hit = br_pmt_hit->GetLeaf("pmt_hit")->GetValue();

  PIp_time = br_PIp_time->GetLeaf("PIp_time")->GetValue();
  PIm_time = br_PIm_time->GetLeaf("PIm_time")->GetValue();
  P_time = br_P_time->GetLeaf("P_time")->GetValue();
  PIp_dist = br_PIp_dist->GetLeaf("PIp_dist")->GetValue();
  PIm_dist = br_PIm_dist->GetLeaf("PIm_dist")->GetValue();
  P_dist = br_P_dist->GetLeaf("P_dist")->GetValue();
  indtype = br_indtype->GetLeaf("indtype")->GetValue();
  
  
  if (br_ectot_el->GetLeaf("ECtot_EL")->GetValue() > (br_ecout_el->GetLeaf("ECout_EL")->GetValue()+br_ecin_el->GetLeaf("ECin_EL")->GetValue())) { 
          ECT  = br_ectot_el->GetLeaf("ECtot_EL")->GetValue();  
         } else {
          ECT  = (br_ecout_el->GetLeaf("ECout_EL")->GetValue()+br_ecin_el->GetLeaf("ECin_EL")->GetValue());    
          };
	  
	  //quatity check   
	     
	     block_curr = br_block->GetLeaf("block")->GetValue(); 

    if ((m_old != 0) && (m != m_old))block_total=block_total+block_prev;
    
    m_old=m;


   LiveTime=br_LiveTime->GetLeaf("LiveTime")->GetValue();
   
   inclusive = br_n_incl->GetLeaf("n_incl")->GetValue()/br_deltaQ->GetLeaf("deltaQ")->GetValue();
   elastic =  br_n_elast->GetLeaf("n_elast")->GetValue()/br_deltaQ->GetLeaf("deltaQ")->GetValue(); 

   if (block_curr != block_prev) {


 
 Qfull_empty = Qfull_empty + br_deltaQ->GetLeaf("deltaQ")->GetValue();
 
   
  
   
   block_prev = block_curr;
   };   
	  
  
 
sector =  int(br_sector->GetLeaf("sector")->GetValue());
nphe = br_nphe_el->GetLeaf("NpheCC_EL")->GetValue();
theta_cc = br_theta_cc->GetLeaf("theta_cc")->GetValue();
ph_cc = br_ph_cc->GetLeaf("ph_cc")->GetValue();  
  


P4_ELP_reg.SetXYZT(P_EL*cos(ph_EL*M_PI/180.)*sin(th_EL*M_PI/180.),P_EL*sin(ph_EL*M_PI/180.)*sin(th_EL*M_PI/180.),P_EL*cos(th_EL*M_PI/180.),P_EL);

P4_PP_reg.SetXYZT(P_P*cos(ph_P*M_PI/180.)*sin(th_P*M_PI/180.),P_P*sin(ph_P*M_PI/180.)*sin(th_P*M_PI/180.),P_P*cos(th_P*M_PI/180.),sqrt(m_proton*m_proton+P_P*P_P));

P4_PIp_reg.SetXYZT(P_PIp*cos(ph_PIp*M_PI/180.)*sin(th_PIp*M_PI/180.),P_PIp*sin(ph_PIp*M_PI/180.)*sin(th_PIp*M_PI/180.),P_PIp*cos(th_PIp*M_PI/180.),sqrt(m_pip*m_pip+P_PIp*P_PIp));


P4_PIm_reg.SetXYZT(P_PIm*cos(ph_PIm*M_PI/180.)*sin(th_PIm*M_PI/180.),P_PIm*sin(ph_PIm*M_PI/180.)*sin(th_PIm*M_PI/180.),P_PIm*cos(th_PIm*M_PI/180.),sqrt(m_pip*m_pip+P_PIm*P_PIm));


P4_PIm_miss = P4_EL + P4_P - P4_ELP_reg - P4_PP_reg - P4_PIp_reg;
P4_PIp_miss = P4_EL + P4_P - P4_ELP_reg - P4_PP_reg - P4_PIm_reg;
P4_PIp_miss_d = P4_EL + P4_D - P4_ELP_reg - P4_PP_reg - P4_PIm_reg;
P4_P_miss = P4_EL + P4_P - P4_ELP_reg - P4_PIp_reg - P4_PIm_reg;
P4_miss_0 = P4_EL + P4_P - P4_ELP_reg - P4_PP_reg - P4_PIm_reg -P4_PIp_reg;

//W_old = W;
//W = sqrt((P4_PIp_reg+P4_PIm_reg+P4_PP_reg).Mag2());

th_PIm_miss = (180./M_PI)*acos(P4_PIm_miss[2]/sqrt(P4_PIm_miss[0]*P4_PIm_miss[0]+P4_PIm_miss[1]*P4_PIm_miss[1]+P4_PIm_miss[2]*P4_PIm_miss[2]));


if (P4_PIm_miss[0] != 0.) {
ph_PIm_miss = (180./M_PI)*atan(P4_PIm_miss[1]/P4_PIm_miss[0]);
}
 else {
if(P4_PIm_miss[1] > 0.) ph_PIm_miss = 90.;
if(P4_PIm_miss[1] < 0.) ph_PIm_miss = 270.;
};
				   

if ((P4_PIm_miss[0] < 0.) && (P4_PIm_miss[1] > 0)) ph_PIm_miss = ph_PIm_miss+180.;
if (( P4_PIm_miss[0]< 0.) && ( P4_PIm_miss[1]< 0)) ph_PIm_miss = ph_PIm_miss+180.;
if ((P4_PIm_miss[0] > 0.) && ( P4_PIm_miss[1]< 0)) ph_PIm_miss = ph_PIm_miss+360.;

beta_func_empty();

bool_el_id_empt=particle_ID_empty.Electron_cuts_empty(); 
bool_proton_id_empt=particle_ID_empty.Proton_cuts_empty();
bool_pip_id_empt=particle_ID_empty.PIp_cuts_empty();
bool_pim_id_empt=particle_ID_empty.PIm_cuts_empty();


selection = false;
selection_pim_miss_empt = false;
selection_0_miss_empt = false;

if (bool_el_id_empt) {


if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)&&(npart >=3)){
if ((bool_proton_id_empt)&&(bool_pim_id_empt)&&(!bool_pip_id_empt)){

if (sqrt((P4_PIp_miss+P4_PIm_reg)*(P4_PIp_miss+P4_PIm_reg)) > 2*m_pip - 0.05) {
if (sqrt((P4_PIp_miss+P4_PIm_reg)*(P4_PIp_miss+P4_PIm_reg)) < W - m_proton + 0.05) {
if (sqrt((P4_PIp_miss+P4_PP_reg)*(P4_PIp_miss+P4_PP_reg)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIp_miss+P4_PP_reg)*(P4_PIp_miss+P4_PP_reg)) < W - m_pip + 0.05) {
if (sqrt((P4_PIm_reg+P4_PP_reg)*(P4_PIm_reg+P4_PP_reg)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIm_reg+P4_PP_reg)*(P4_PIm_reg+P4_PP_reg)) < W - m_pip + 0.05) {

if (P4_PIp_miss[3] > m_pip  - 0.05) {

//if (P4_PIp_miss.Mag2() > m_pip*m_pip - 4*0.0302){
selection = false;
//P4_PIp_reg = P4_PIp_miss;
//};
};
};
};
};
};
};
};
};


};

if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)&&(npart >=3)){
if ((bool_pip_id_empt)&&(bool_pim_id_empt)&&(!bool_proton_id_empt)){

if (sqrt((P4_PIp_reg+P4_PIm_reg)*(P4_PIp_reg+P4_PIm_reg)) > 2*m_pip - 0.05) {
if (sqrt((P4_PIp_reg+P4_PIm_reg)*(P4_PIp_reg+P4_PIm_reg)) < W - m_proton + 0.05) {
if (sqrt((P4_PIp_reg+P4_P_miss)*(P4_PIp_reg+P4_P_miss)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIp_reg+P4_P_miss)*(P4_PIp_reg+P4_P_miss)) < W - m_pip + 0.05) {
if (sqrt((P4_PIm_reg+P4_P_miss)*(P4_PIm_reg+P4_P_miss)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIm_reg+P4_P_miss)*(P4_PIm_reg+P4_P_miss)) < W - m_pip + 0.05) {
if (P4_P_miss[3] > m_proton  - 0.05) {

//if (P4_P_miss.Mag2() >m_proton*m_proton - 4*0.103 ){
selection = false;
//P4_PP_reg = P4_P_miss;
//};
};
};
};
};
};
};
};
};
};
if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)&&(npart>=4)){
//if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)){
if ((bool_pip_id_empt)&&(bool_proton_id_empt)&&(bool_pim_id_empt)){

if (sector ==1)  hist_z_el_1_empty->Fill(z_EL-0.15,1.); 
if (sector ==2)  hist_z_el_2_empty->Fill(z_EL-0.15,1.);
if (sector ==3)  hist_z_el_3_empty->Fill(z_EL-0.15,1.);
if (sector ==4)  hist_z_el_4_empty->Fill(z_EL-0.15,1.);
if (sector ==5)  hist_z_el_5_empty->Fill(z_EL-0.15,1.);
if (sector ==6)  hist_z_el_6_empty->Fill(z_EL-0.15,1.);



if (sqrt((P4_PIp_reg+P4_PIm_reg)*(P4_PIp_reg+P4_PIm_reg)) > 2*m_pip - 0.05) {
if (sqrt((P4_PIp_reg+P4_PIm_reg)*(P4_PIp_reg+P4_PIm_reg)) < W - m_proton + 0.05) {
if (sqrt((P4_PIp_reg+P4_PP_reg)*(P4_PIp_reg+P4_PP_reg)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIp_reg+P4_PP_reg)*(P4_PIp_reg+P4_PP_reg)) < W - m_pip + 0.05) {
if (sqrt((P4_PIm_reg+P4_PP_reg)*(P4_PIm_reg+P4_PP_reg)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIm_reg+P4_PP_reg)*(P4_PIm_reg+P4_PP_reg)) < W - m_pip + 0.05) {

if ((P4_miss_0[3] >  - 0.05)) {

if ((P4_miss_0.Mag2()>-0.05)&&(P4_miss_0.Mag2()<0.01)){
selection_0_miss_empt = true;

};
};
};
};
};
};
};
};
};
};


if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)&&(npart >=3)){
//if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)){
if ((bool_pip_id_empt)&&(bool_proton_id_empt)&&(!bool_pim_id_empt)){
if (sqrt((P4_PIp_reg+P4_PIm_miss)*(P4_PIp_reg+P4_PIm_miss)) > 2*m_pip - 0.05) {
if (sqrt((P4_PIp_reg+P4_PIm_miss)*(P4_PIp_reg+P4_PIm_miss)) < W - m_proton + 0.05) {
if (sqrt((P4_PIp_reg+P4_PP_reg)*(P4_PIp_reg+P4_PP_reg)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIp_reg+P4_PP_reg)*(P4_PIp_reg+P4_PP_reg)) < W - m_pip + 0.05) {
if (sqrt((P4_PIm_miss+P4_PP_reg)*(P4_PIm_miss+P4_PP_reg)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIm_miss+P4_PP_reg)*(P4_PIm_miss+P4_PP_reg)) < W - m_pip + 0.05) {
if (P4_PIm_miss[3] > m_pip  - 0.05) {


if ((P4_PIm_miss.Mag2()> -0.01)&&(P4_PIm_miss.Mag2()< 0.05)){
selection_pim_miss_empt = true;
P4_PIm_reg = P4_PIm_miss;
};
};
};
};
};
};
};
};
};
};





if ((selection_pim_miss_empt)||(selection_0_miss_empt)) {

rot_boost_cmsyst();
/*
if((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)){
Var_1[0] = inv_m_pip_p;
Var_1[1] = inv_m_pip_pim;
Var_1[2] = theta_P_cm;
Var_1[3] = phi_P_cm;
Var_1[4] = alpha_PIpPIm_pipf;

Var_2[0] = inv_m_pip_p;
Var_2[1] = inv_m_pip_pim;
Var_2[2] = theta_PIm_cm;
Var_2[3] = phi_PIm_cm;
Var_2[4] = alpha_PPIp_piPIm;

Var_3[0] = inv_m_pim_p;
Var_3[1] = inv_m_pip_pim;
Var_3[2] = theta_PIp_cm;
Var_3[3] = phi_PIp_cm;
Var_3[4] = alpha_PPIm_piPIp;

//cout <<"ereee" << pmt_hit+1 << "  "<<sector-1 << "  "<< segment <<"   " << fract_integ[pmt_hit+1][sector-1][segment] << "\n";

//Q_empty = 464.797 micro C
h_5dim_1_empty[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_1,fract_integ[pmt_hit+1][sector-1][segment]);
h_5dim_2_empty[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_2,fract_integ[pmt_hit+1][sector-1][segment]);
h_5dim_3_empty[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_3,fract_integ[pmt_hit+1][sector-1][segment]); 

if ((selection_0_miss_empt)) {
if ((W > 1.4)&&(W < 1.625)&&(Q2 > 0.4)&&(Q2 < 0.9)) {

h_5dim_excl_1_empty[int((W-1.4)/0.025)]->Fill(Var_1,fract_integ[pmt_hit+1][sector-1][segment]);
h_5dim_excl_2_empty[int((W-1.4)/0.025)]->Fill(Var_2,fract_integ[pmt_hit+1][sector-1][segment]);
h_5dim_excl_3_empty[int((W-1.4)/0.025)]->Fill(Var_3,fract_integ[pmt_hit+1][sector-1][segment]);
};
};
};*/
};






 }; //konec ifa electronnih cutov
    };// konets tsicla po sobitiyam (po i)
  
  
  t21->Delete();

 finp->Close();  
 
     };//konets tsicla po failam (po m)  
  
  
  
  
  
  
//////////////////////////////////////////-sim  
 block_curr = 0;
 block_prev = 0;
 m_old=0;
 
for (m=1; m<=n_files_sim; m++) {
  
  finp = new TFile(file_sim[m-1].c_str()); 
  
  
  
 cout << "Processing file with simulation " << m << "\n"; 
  
  
  // cout << " bin content3 = " << avrg_nphe_sector1->GetBinContent(100,100) << "\n";
    
  TTree *t21 = (TTree*)finp->Get("t21");
  
    TBranch *br_indtype = t21->GetBranch("indtype");
    TBranch *br_npart = t21->GetBranch("npart");
    TBranch *br_p_el = t21->GetBranch("P_EL");
    TBranch *br_block = t21->GetBranch("block");
    TBranch *br_deltaQ = t21->GetBranch("deltaQ");
    TBranch *br_LiveTime = t21->GetBranch("LiveTime");
    TBranch *br_n_incl = t21->GetBranch("n_incl");
    TBranch *br_n_elast = t21->GetBranch("n_elast");
    TBranch *br_ph_EL = t21->GetBranch("ph_EL");
    TBranch *br_th_EL = t21->GetBranch("th_EL");
    TBranch *br_W = t21->GetBranch("W");
    TBranch *br_Q2 = t21->GetBranch("Q2");
    TBranch *br_nphe_el = t21->GetBranch("NpheCC_EL");
    TBranch *br_ectot_el = t21->GetBranch("ECtot_EL");
    TBranch *br_ecin_el = t21->GetBranch("ECin_EL");
    TBranch *br_ecout_el = t21->GetBranch("ECout_EL");
    TBranch *br_x_EL = t21->GetBranch("x_EL");
    TBranch *br_y_EL = t21->GetBranch("y_EL");
    TBranch *br_z_EL = t21->GetBranch("z_EL");
    TBranch *br_z_P = t21->GetBranch("z_P");
    TBranch *br_z_PIp = t21->GetBranch("z_PIp");
    TBranch *br_z_PIm = t21->GetBranch("z_PIm");
    TBranch *br_PdHit_EL = t21->GetBranch("PdHit_EL");
    TBranch *br_PdHit_PIp = t21->GetBranch("PdHit_PIp");
    TBranch *br_PdHit_PIm = t21->GetBranch("PdHit_PIm");
    TBranch *br_PdHit_P = t21->GetBranch("PdHit_P");
    TBranch *br_sc_x = t21->GetBranch("sc_x");
    TBranch *br_sc_y = t21->GetBranch("sc_y");
    TBranch *br_sc_z = t21->GetBranch("sc_z");
    TBranch *br_pmt_hit = t21->GetBranch("pmt_hit");
    TBranch *br_segment = t21->GetBranch("segment");
    TBranch *br_theta_cc = t21->GetBranch("theta_cc");
    TBranch *br_ph_cc = t21->GetBranch("ph_cc");
    TBranch *br_sector = t21->GetBranch("sector");
    TBranch *br_n_PIp = t21->GetBranch("n_PIp");
    TBranch *br_n_PIm = t21->GetBranch("n_PIm");
    TBranch *br_n_P = t21->GetBranch("n_P");
    TBranch *br_p_pip = t21->GetBranch("P_PIp");
    TBranch *br_p_pim = t21->GetBranch("P_PIm");
    TBranch *br_p_p = t21->GetBranch("P_P");
    TBranch *br_th_PIp = t21->GetBranch("th_PIp");
    TBranch *br_th_PIm = t21->GetBranch("th_PIm");
    TBranch *br_th_P = t21->GetBranch("th_P");
    TBranch *br_ph_PIp = t21->GetBranch("ph_PIp");
    TBranch *br_ph_PIm = t21->GetBranch("ph_PIm");
    TBranch *br_ph_P = t21->GetBranch("ph_P");
    TBranch *br_beta_PIm = t21->GetBranch("beta_PIm");
    TBranch *br_beta_PIp = t21->GetBranch("beta_PIp");
    TBranch *br_beta_P = t21->GetBranch("beta_P");
    TBranch *br_nphe_pip = t21->GetBranch("NpheCC_PIp");
    TBranch *br_PIp_time = t21->GetBranch("PIp_time");
    TBranch *br_PIp_dist = t21->GetBranch("PIp_dist");
    TBranch *br_P_time = t21->GetBranch("P_time");
    TBranch *br_P_dist = t21->GetBranch("P_dist");
    TBranch *br_PIm_time = t21->GetBranch("PIm_time");
    TBranch *br_PIm_dist = t21->GetBranch("PIm_dist"); 
     TBranch *br_sigma = t21->GetBranch("sigma");  
//    TBranch *br_pf_x = t21->GetBranch("pf_x");  
//    TBranch *br_pf_y = t21->GetBranch("pf_y"); 
 //   TBranch *br_pf_z = t21->GetBranch("pf_z");
    TBranch *br_beta_PIm_time = t21->GetBranch("beta_PIm_time");
    TBranch *br_beta_PIp_time = t21->GetBranch("beta_PIp_time");
    TBranch *br_beta_P_time = t21->GetBranch("beta_P_time"); 
  

  Bool_t adc_cut_switch,tdc_cut_switch;
  
  Qdiff = 0.;
  Qcurr = 0.;
  Qprev = 0.;
  Qtotal = 0.;
  k = 0;
  block = 0;
  last_k = 0;
  nstart = 0;
  nstop = 0;
  n_incl = 0;
  n_elast = 0;
  
  
  for (i=0; i<br_sector->GetEntries(); i++) { 
  
  Qprev = Qcurr;

br_segment->GetEntry(i);
br_pmt_hit->GetEntry(i);
 br_deltaQ->GetEntry(i);
  br_n_incl->GetEntry(i);
  br_n_elast->GetEntry(i);  
  br_LiveTime->GetEntry(i);
  br_block->GetEntry(i);  
  br_W->GetEntry(i);
  br_Q2->GetEntry(i); 
  br_npart->GetEntry(i);
  br_n_PIm->GetEntry(i);
  br_n_PIp->GetEntry(i);
  br_n_P->GetEntry(i);
  br_sector->GetEntry(i);
  br_LiveTime->GetEntry(i);
  br_block->GetEntry(i);  
  br_ph_cc->GetEntry(i);
  br_theta_cc->GetEntry(i);  
  br_p_el->GetEntry(i);  
  br_p_p->GetEntry(i); 
  br_p_pim->GetEntry(i);
  br_p_pip->GetEntry(i); 
  br_ectot_el->GetEntry(i); 
  br_ecin_el->GetEntry(i); 
  br_ecout_el->GetEntry(i);
  br_nphe_el->GetEntry(i);
  br_nphe_pip->GetEntry(i);
  br_th_EL->GetEntry(i);
  br_ph_EL->GetEntry(i); 
  br_th_P->GetEntry(i);
  br_ph_P->GetEntry(i);
  br_th_PIp->GetEntry(i);
  br_ph_PIp->GetEntry(i);
  br_th_PIm->GetEntry(i);
  br_ph_PIm->GetEntry(i);
  br_beta_PIp->GetEntry(i);
  br_beta_PIm->GetEntry(i);
  br_beta_P->GetEntry(i);
  br_PdHit_EL->GetEntry(i);
  br_PdHit_PIp->GetEntry(i);
  br_PdHit_PIm->GetEntry(i);
  br_PdHit_P->GetEntry(i);
  br_z_EL->GetEntry(i);
  br_z_P->GetEntry(i);
  br_z_PIp->GetEntry(i);
  br_z_PIm->GetEntry(i);
  br_W->GetEntry(i);
  br_Q2->GetEntry(i);
  br_indtype->GetEntry(i);
  
  
  br_PIp_time->GetEntry(i);
  br_PIm_time->GetEntry(i);
  br_P_time->GetEntry(i); 
  br_PIp_dist->GetEntry(i);
  br_PIm_dist->GetEntry(i);
  br_P_dist->GetEntry(i);
  br_sigma->GetEntry(i);
  br_beta_PIp_time->GetEntry(i);
  br_beta_PIm_time->GetEntry(i);
  br_beta_P_time->GetEntry(i);
   
   
//   br_pf_x->GetEntry(i);
 // br_pf_y->GetEntry(i);
//  br_pf_z->GetEntry(i);
  
  P_EL = br_p_el->GetLeaf("P_EL")->GetValue();
  P_P = br_p_p->GetLeaf("P_P")->GetValue();
  P_PIp = br_p_pip->GetLeaf("P_PIp")->GetValue();
  P_PIm = br_p_pim->GetLeaf("P_PIm")->GetValue();
  th_EL = br_th_EL->GetLeaf("th_EL")->GetValue();
  ph_EL = br_ph_EL->GetLeaf("ph_EL")->GetValue();
  th_P = br_th_P->GetLeaf("th_P")->GetValue();
  ph_P = br_ph_P->GetLeaf("ph_P")->GetValue();
  th_PIp = br_th_PIp->GetLeaf("th_PIp")->GetValue();
  ph_PIp = br_ph_PIp->GetLeaf("ph_PIp")->GetValue();
  npart = br_npart->GetLeaf("npart")->GetValue();
  n_P = br_n_P->GetLeaf("n_P")->GetValue();
  n_PIp = br_n_PIp->GetLeaf("n_PIp")->GetValue();
  Nphe_pip = br_nphe_pip->GetLeaf("NpheCC_PIp")->GetValue();
  PdHit_EL = br_PdHit_EL->GetLeaf("PdHit_EL")->GetValue();
  PdHit_PIp = br_PdHit_PIp->GetLeaf("PdHit_PIp")->GetValue();
  PdHit_PIm = br_PdHit_PIm->GetLeaf("PdHit_PIm")->GetValue();
  PdHit_P = br_PdHit_P->GetLeaf("PdHit_P")->GetValue();
//  beta_PIp = br_beta_PIp->GetLeaf("beta_PIp")->GetValue();
//  beta_PIm = br_beta_PIm->GetLeaf("beta_PIm")->GetValue();
//  beta_P = br_beta_P->GetLeaf("beta_P")->GetValue();
  n_PIm = br_n_PIm->GetLeaf("n_PIm")->GetValue();
  ph_PIm = br_ph_PIm->GetLeaf("ph_PIm")->GetValue();
  th_PIm = br_th_PIm->GetLeaf("th_PIm")->GetValue();
  z_EL = br_z_EL->GetLeaf("z_EL")->GetValue();
  z_P = br_z_P->GetLeaf("z_P")->GetValue();
  z_PIp = br_z_PIp->GetLeaf("z_PIp")->GetValue();
  z_PIm = br_z_PIm->GetLeaf("z_PIm")->GetValue();
  W = br_W->GetLeaf("W")->GetValue();
  Q2 = br_Q2->GetLeaf("Q2")->GetValue();
  segment = br_segment->GetLeaf("segment")->GetValue();
  pmt_hit = br_pmt_hit->GetLeaf("pmt_hit")->GetValue();
  
 
  PIp_time = br_PIp_time->GetLeaf("PIp_time")->GetValue();
  PIm_time = br_PIm_time->GetLeaf("PIm_time")->GetValue();
  P_time = br_P_time->GetLeaf("P_time")->GetValue();
  PIp_dist = br_PIp_dist->GetLeaf("PIp_dist")->GetValue();
  PIm_dist = br_PIm_dist->GetLeaf("PIm_dist")->GetValue();
  P_dist = br_P_dist->GetLeaf("P_dist")->GetValue();
  indtype = br_indtype->GetLeaf("indtype")->GetValue();
  sigma = br_sigma->GetLeaf("sigma")->GetValue();
//  pf_x = br_pf_x->GetLeaf("pf_x")->GetValue();
 //  pf_y = br_pf_y->GetLeaf("pf_y")->GetValue();
 //  pf_z = br_pf_z->GetLeaf("pf_z")->GetValue();
  beta_P = br_beta_P_time->GetLeaf("beta_P_time")->GetValue();
  beta_PIp = br_beta_PIp_time->GetLeaf("beta_PIp_time")->GetValue();
  beta_PIm = br_beta_PIm_time->GetLeaf("beta_PIm_time")->GetValue();

 
//----------ELECTRON MOMENTUM CORRECTION FOR SIM---------- 
  if (indtype==1){
  delta_p_el_sim = corrfunc.corr_el_mom_sim(P_EL, th_EL);
//  P_EL = P_EL + delta_p_el_sim;
  };
//--------------------------------------------------------

//----------PROTON ENERGY LOSS Ye unfold------------------
delta_mom_p_ye = corrfunc.correct_energy_theta_pf(P_P, th_P);
P_P = P_P - delta_mom_p_ye;
//----------PROTON ENERGY LOSS----------------------------
delta_mom_p_skor = corrfunc.corr_pr_mom_skor(P_P, th_P);
//--------------------------------------------------------

  
  if (br_ectot_el->GetLeaf("ECtot_EL")->GetValue() > (br_ecout_el->GetLeaf("ECout_EL")->GetValue()+br_ecin_el->GetLeaf("ECin_EL")->GetValue())) { 
          ECT  = br_ectot_el->GetLeaf("ECtot_EL")->GetValue();  
         } else {
          ECT  = (br_ecout_el->GetLeaf("ECout_EL")->GetValue()+br_ecin_el->GetLeaf("ECin_EL")->GetValue());    
          };
	  
	  //quatity check   
	     
	     block_curr = br_block->GetLeaf("block")->GetValue(); 

    if ((m_old != 0) && (m != m_old))block_total=block_total+block_prev;
    
    m_old=m;


   LiveTime=br_LiveTime->GetLeaf("LiveTime")->GetValue();
   
   inclusive = br_n_incl->GetLeaf("n_incl")->GetValue()/br_deltaQ->GetLeaf("deltaQ")->GetValue();
   elastic =  br_n_elast->GetLeaf("n_elast")->GetValue()/br_deltaQ->GetLeaf("deltaQ")->GetValue(); 

   if (block_curr != block_prev) {


 
 Qfull_sim = Qfull_sim + br_deltaQ->GetLeaf("deltaQ")->GetValue();
 
   
  
   
   block_prev = block_curr;
   };   
	  
  
 
sector =  int(br_sector->GetLeaf("sector")->GetValue());
nphe = br_nphe_el->GetLeaf("NpheCC_EL")->GetValue();
theta_cc = br_theta_cc->GetLeaf("theta_cc")->GetValue();
ph_cc = br_ph_cc->GetLeaf("ph_cc")->GetValue();  
  
//fermi_bonn();
//cout << px_fermi << "  "<<py_fermi<< " "<< pz_fermi<< "\n";
P4_P.SetXYZT(0,0,0,m_proton);

P4_EL.SetXYZT(0,0,2.039,2.039);
	
//P4_P.SetXYZT(px_fermi,py_fermi,pz_fermi,sqrt(m_proton*m_proton+px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi));
//cout <<"segment = " << segment <<"\n";
//cout <<" " << pf_x <<"  "<<pf_y<<" "<<pf_z<<"\n";

//pxel_new = P_EL*cos(ph_EL*M_PI/180.)*sin(th_EL*M_PI/180.)-pf_x;
//pyel_new = P_EL*sin(ph_EL*M_PI/180.)*sin(th_EL*M_PI/180.)-pf_y;
//pzel_new = P_EL*cos(th_EL*M_PI/180.)-pf_z;

//P4_ELP_reg.SetXYZT(P_EL*cos(ph_EL*M_PI/180.)*sin(th_EL*M_PI/180.),P_EL*sin(ph_EL*M_PI/180.)*sin(th_EL*M_PI/180.),P_EL*cos(th_EL*M_PI/180.),sqrt(pxel_new*pxel_new+pyel_new*pyel_new+pzel_new*pzel_new));
//P4_ELP_reg.SetXYZT(P_EL*cos(ph_EL*M_PI/180.)*sin(th_EL*M_PI/180.),P_EL*sin(ph_EL*M_PI/180.)*sin(th_EL*M_PI/180.),P_EL*cos(th_EL*M_PI/180.),sqrt(P_EL*P_EL+pf_x*pf_x+pf_y*pf_y+pf_z*pf_z));
//P4_ELP_for_miss_en_comp.SetXYZT(P_EL*cos(ph_EL*M_PI/180.)*sin(th_EL*M_PI/180.),P_EL*sin(ph_EL*M_PI/180.)*sin(th_EL*M_PI/180.),P_EL*cos(th_EL*M_PI/180.),sqrt(pxel_new*pxel_new+pyel_new*pyel_new+pzel_new*pzel_new));
//cout << th_PIp<<" hh\n";
P4_ELP_reg.SetXYZT(P_EL*cos(ph_EL*M_PI/180.)*sin(th_EL*M_PI/180.),P_EL*sin(ph_EL*M_PI/180.)*sin(th_EL*M_PI/180.),P_EL*cos(th_EL*M_PI/180.),P_EL);

P4_PP_reg.SetXYZT(P_P*cos(ph_P*M_PI/180.)*sin(th_P*M_PI/180.),P_P*sin(ph_P*M_PI/180.)*sin(th_P*M_PI/180.),P_P*cos(th_P*M_PI/180.),sqrt(m_proton*m_proton+P_P*P_P));

P4_PIp_reg.SetXYZT(P_PIp*cos(ph_PIp*M_PI/180.)*sin(th_PIp*M_PI/180.),P_PIp*sin(ph_PIp*M_PI/180.)*sin(th_PIp*M_PI/180.),P_PIp*cos(th_PIp*M_PI/180.),sqrt(m_pip*m_pip+P_PIp*P_PIp));


P4_PIm_reg.SetXYZT(P_PIm*cos(ph_PIm*M_PI/180.)*sin(th_PIm*M_PI/180.),P_PIm*sin(ph_PIm*M_PI/180.)*sin(th_PIm*M_PI/180.),P_PIm*cos(th_PIm*M_PI/180.),sqrt(m_pip*m_pip+P_PIm*P_PIm));






P4_PIm_miss = P4_EL + P4_P - P4_ELP_reg - P4_PP_reg - P4_PIp_reg;
//P4_PIm_miss_en_comp = P4_EL + P4_P - P4_ELP_for_miss_en_comp - P4_PP_reg - P4_PIp_reg;
P4_PIp_miss = P4_EL + P4_P - P4_ELP_reg - P4_PP_reg - P4_PIm_reg;
P4_PIp_miss_d = P4_EL + P4_D - P4_ELP_reg - P4_PP_reg - P4_PIm_reg;
P4_P_miss = P4_EL + P4_P - P4_ELP_reg - P4_PIp_reg - P4_PIm_reg;
P4_miss_0 = P4_EL + P4_P - P4_ELP_reg - P4_PP_reg - P4_PIm_reg -P4_PIp_reg;
//P4_miss_0_en_comp = P4_EL + P4_P - P4_ELP_for_miss_en_comp - P4_PP_reg - P4_PIm_reg -P4_PIp_reg;
P4_miss_0_d = P4_EL + P4_D - P4_ELP_reg - P4_PP_reg - P4_PIm_reg -P4_PIp_reg;



//if (indtype==2) cout << P4_PIm_miss.Mag()<<"\n";
W_old = W;
//W = sqrt((P4_PIp_reg+P4_PIm_reg+P4_PP_reg).Mag2());


//if ((P4_PIp_reg+P4_PIm_reg+P4_PP_reg).Mag2()>0) W=sqrt((P4_PIp_reg+P4_PIm_reg+P4_PP_reg).Mag2());
//if ((P4_PIp_reg+P4_PIm_reg+P4_PP_reg).Mag2()<=0) W=0;


//if ((P4_P+P4_EL-P4_ELP_reg).Mag2()>0) W=sqrt((P4_P+P4_EL-P4_ELP_reg).Mag2());
//if ((P4_P+P4_EL-P4_ELP_reg).Mag2()<=0) W=0;

th_PIm_miss = (180./M_PI)*acos(P4_PIm_miss[2]/sqrt(P4_PIm_miss[0]*P4_PIm_miss[0]+P4_PIm_miss[1]*P4_PIm_miss[1]+P4_PIm_miss[2]*P4_PIm_miss[2]));
/*if(P4_PIm_miss[0] == 0.) {
if(P4_PIm_miss[1] > 0.) ph_PIm_miss = 90.;
if(P4_PIm_miss[1] < 0.) ph_PIm_miss = 270.;

};


if ((P4_PIm_miss[1]>0)&&(P4_PIm_miss[0]>0)&&(P4_PIm_miss[0] != 0.)) {
ph_PIm_miss = (180./M_PI)*atan(P4_PIm_miss[1]/P4_PIm_miss[0]);
};


if ((P4_PIm_miss[1]>0)&&(P4_PIm_miss[0]<0)&&(P4_PIm_miss[0] != 0.)) {
ph_PIm_miss = 180. - (180./M_PI)*atan(abs(P4_PIm_miss[1]/P4_PIm_miss[0]));
};

if ((P4_PIm_miss[1]<0)&&(P4_PIm_miss[0]<0)&&(P4_PIm_miss[0] != 0.)) {
ph_PIm_miss = 180. + (180./M_PI)*atan(abs(P4_PIm_miss[1]/P4_PIm_miss[0]));
};

if ((P4_PIm_miss[1]<0)&&(P4_PIm_miss[0]>0)&&(P4_PIm_miss[0] != 0.)) {
ph_PIm_miss = 360. - (180./M_PI)*atan(abs(P4_PIm_miss[1]/P4_PIm_miss[0]));
};*/

if (P4_PIm_miss[0] != 0.) {
ph_PIm_miss = (180./M_PI)*atan(P4_PIm_miss[1]/P4_PIm_miss[0]);
}
 else {
if(P4_PIm_miss[1] > 0.) ph_PIm_miss = 90.;
if(P4_PIm_miss[1] < 0.) ph_PIm_miss = 270.;
};
				   

if ((P4_PIm_miss[0] < 0.) && (P4_PIm_miss[1] > 0)) ph_PIm_miss = ph_PIm_miss+180.;
if (( P4_PIm_miss[0]< 0.) && ( P4_PIm_miss[1]< 0)) ph_PIm_miss = ph_PIm_miss+180.;
if ((P4_PIm_miss[0] > 0.) && ( P4_PIm_miss[1]< 0)) ph_PIm_miss = ph_PIm_miss+360.;


//if (indtype==2) h_PIm_miss_sim-> Fill(P4_PIm_miss.Mag2(),1.);
if((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)){
if (indtype==2){
rot_boost_cmsyst();

Var_1[0] = inv_m_pip_p;
Var_1[1] = inv_m_pip_pim;
Var_1[2] = theta_P_cm;
Var_1[3] = phi_P_cm;
Var_1[4] = alpha_PIpPIm_pipf;

Var_2[0] = inv_m_pip_p;
Var_2[1] = inv_m_pip_pim;
Var_2[2] = theta_PIm_cm;
Var_2[3] = phi_PIm_cm;
Var_2[4] = alpha_PPIp_piPIm;

Var_3[0] = inv_m_pim_p;
Var_3[1] = inv_m_pip_pim;
Var_3[2] = theta_PIp_cm;
Var_3[3] = phi_PIp_cm;
Var_3[4] = alpha_PPIm_piPIp;

h_5dim_1_sim_gen[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_1,sigma);
h_5dim_2_sim_gen[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_2,sigma);
h_5dim_3_sim_gen[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_3,sigma); 

h_5dim_1_sim_gen_evt[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_1,sigma*sigma);
h_5dim_2_sim_gen_evt[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_2,sigma*sigma);
h_5dim_3_sim_gen_evt[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_3,sigma*sigma);






if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 0.45)) h_1d_rc_0425->Fill(W,sigma);
if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.45)&&(Q2 < 0.5)) h_1d_rc_0475->Fill(W,sigma);
if ((W > 1.3)&&(W < 1.8)&&(Q2 > 0.5)&&(Q2 < 0.55)) h_1d_rc_0525->Fill(W,sigma);
if ((W > 1.3)&&(W < 1.8)&&(Q2 > 0.55)&&(Q2 < 0.6)) h_1d_rc_0575->Fill(W,sigma);
if ((W > 1.3)&&(W < 1.775)&&(Q2 > 0.6)&&(Q2 < 0.65)) h_1d_rc_0625->Fill(W,sigma);
if ((W > 1.3)&&(W < 1.75)&&(Q2 > 0.65)&&(Q2 < 0.7)) h_1d_rc_0675->Fill(W,sigma);
if ((W > 1.3)&&(W < 1.725)&&(Q2 > 0.7)&&(Q2 < 0.75)) h_1d_rc_0725->Fill(W,sigma);
if ((W > 1.3)&&(W < 1.7)&&(Q2 > 0.75)&&(Q2 < 0.8)) h_1d_rc_0775->Fill(W,sigma);
if ((W > 1.3)&&(W < 1.65)&&(Q2 > 0.8)&&(Q2 < 0.85)) h_1d_rc_0825->Fill(W,sigma);
if ((W > 1.3)&&(W < 1.625)&&(Q2 > 0.85)&&(Q2 < 0.9)) h_1d_rc_0875->Fill(W,sigma);
if ((W > 1.3)&&(W < 1.6)&&(Q2 > 0.9)&&(Q2 < 0.95)) h_1d_rc_0925->Fill(W,sigma);
if ((W > 1.3)&&(W < 1.55)&&(Q2 > 0.95)&&(Q2 < 1.)) h_1d_rc_0975->Fill(W,sigma);

/*
if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 0.45)) h_1d_rc_0425_evt->Fill(W,1.);
if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.45)&&(Q2 < 0.5)) h_1d_rc_0475_evt->Fill(W,1.);
if ((W > 1.3)&&(W < 1.8)&&(Q2 > 0.5)&&(Q2 < 0.55)) h_1d_rc_0525_evt->Fill(W,1.);
if ((W > 1.3)&&(W < 1.8)&&(Q2 > 0.55)&&(Q2 < 0.6)) h_1d_rc_0575_evt->Fill(W,1.);
if ((W > 1.3)&&(W < 1.775)&&(Q2 > 0.6)&&(Q2 < 0.65)) h_1d_rc_0625_evt->Fill(W,1.);
if ((W > 1.3)&&(W < 1.75)&&(Q2 > 0.65)&&(Q2 < 0.7)) h_1d_rc_0675_evt->Fill(W,1.);
if ((W > 1.3)&&(W < 1.725)&&(Q2 > 0.7)&&(Q2 < 0.75)) h_1d_rc_0725_evt->Fill(W,1.);
if ((W > 1.3)&&(W < 1.7)&&(Q2 > 0.75)&&(Q2 < 0.8)) h_1d_rc_0775_evt->Fill(W,1.);
if ((W > 1.3)&&(W < 1.65)&&(Q2 > 0.8)&&(Q2 < 0.85)) h_1d_rc_0825_evt->Fill(W,1.);
if ((W > 1.3)&&(W < 1.625)&&(Q2 > 0.85)&&(Q2 < 0.9)) h_1d_rc_0875_evt->Fill(W,1.);
if ((W > 1.3)&&(W < 1.6)&&(Q2 > 0.9)&&(Q2 < 0.95)) h_1d_rc_0925_evt->Fill(W,1.);
if ((W > 1.3)&&(W < 1.55)&&(Q2 > 0.95)&&(Q2 < 1.)) h_1d_rc_0975_evt->Fill(W,1.);
*/



};
};


bool_el_id_sim=particle_ID_sim.Electron_cuts_sim();
bool_proton_id_sim=particle_ID_sim.Proton_cuts_sim();
bool_pip_id_sim=particle_ID_sim.PIp_cuts_sim();
bool_pim_id_sim=particle_ID_sim.PIm_cuts_sim();

selection = false;
selection_pim_miss_sim = false;
selection_0_miss_sim = false;

if (bool_el_id_sim) {

//vertex difference cut
if ((abs(z_EL - z_P)<5.)&&(abs(z_EL - z_PIp)<5.)&&(abs(z_P - z_PIp)<5.)){

//if ((bool_proton_id_sim)&&(bool_pip_id_sim)) W_2pi_selection_sim->Fill(W,Q2,sigma);
sim_hist();

if (particle_ID_sim.PIp_cuts_sim()) h_z_corr1_sim->Fill(z_PIp);
if ((particle_ID_sim.PIp_cuts_sim())&&(particle_ID_sim.PIm_cuts_sim())) h_z_corr2_sim->Fill(z_PIm-z_PIp);
if ((bool_proton_id_sim)&&(bool_pim_id_sim)&&(!bool_pip_id_sim)){

//if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)&&(npart>=3)) {
if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)) {

h_PIp_miss_d_bef_sim-> Fill(P4_PIp_miss_d.Mag2(),1.);
if (sqrt((P4_PIp_miss+P4_PIm_reg)*(P4_PIp_miss+P4_PIm_reg)) > 2*m_pip - 0.05) {
if (sqrt((P4_PIp_miss+P4_PIm_reg)*(P4_PIp_miss+P4_PIm_reg)) < W - m_proton + 0.05) {
if (sqrt((P4_PIp_miss+P4_PP_reg)*(P4_PIp_miss+P4_PP_reg)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIp_miss+P4_PP_reg)*(P4_PIp_miss+P4_PP_reg)) < W - m_pip + 0.05) {
if (sqrt((P4_PIm_reg+P4_PP_reg)*(P4_PIm_reg+P4_PP_reg)) > m_pip + m_proton  - 0.05) {
if (sqrt((P4_PIm_reg+P4_PP_reg)*(P4_PIm_reg+P4_PP_reg)) < W - m_pip + 0.05) {

h_mm_pip_vs_npart_sim->Fill(npart,P4_PIp_miss.Mag2(),1.);
 h_PIp_miss_en_sim->Fill(P4_PIp_miss[3],1.);

if (P4_PIp_miss[3] > 0.15) {

h_PIp_miss_sim-> Fill(P4_PIp_miss.Mag2(),1.);
h_PIp_miss_d_sim-> Fill(P4_PIp_miss_d.Mag2(),1.);
selection = false;
//P4_PIp_reg = P4_PIp_miss;
//};
};
};
};
};
};
};
};
};


};



if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)) {
if ((bool_pip_id_sim)&&(bool_proton_id_sim)&&(bool_pim_id_sim)&&(npart>=4)){
//excl top vertex diffrence cut 
if ((abs(z_EL - z_PIm)<5.)&&(abs(z_P - z_PIm)<5.)&&(abs(z_PIp - z_PIm)<5.)){
//cout << P4_PIm_reg[0] << " "<<P4_PIm_reg[1] <<" "<<P4_PIm_reg[2] << " "<<P4_PIm_reg[3] << " reg \n";
//cout << P4_PIm_miss[0] << " "<<P4_PIm_miss[1] <<" " <<P4_PIm_miss[2] << " "<<P4_PIm_miss[3] << " miss  \n";


if (sqrt((P4_PIp_reg+P4_PIm_reg)*(P4_PIp_reg+P4_PIm_reg)) > 2*m_pip - 0.01) {
if (sqrt((P4_PIp_reg+P4_PIm_reg)*(P4_PIp_reg+P4_PIm_reg)) < W - m_proton + 0.01) {
if (sqrt((P4_PIp_reg+P4_PP_reg)*(P4_PIp_reg+P4_PP_reg)) > m_pip + m_proton  - 0.01) {
if (sqrt((P4_PIp_reg+P4_PP_reg)*(P4_PIp_reg+P4_PP_reg)) < W - m_pip + 0.01) {
if (sqrt((P4_PIm_reg+P4_PP_reg)*(P4_PIm_reg+P4_PP_reg)) > m_pip + m_proton  - 0.01) {
if (sqrt((P4_PIm_reg+P4_PP_reg)*(P4_PIm_reg+P4_PP_reg)) < W - m_pip + 0.01) {


h_mm_0_vs_npart_sim->Fill(npart,P4_miss_0.Mag2(),1.);


th_ph_pim = 180./M_PI*acos((P4_PIm_reg.Vect().Dot((P4_EL-P4_ELP_reg).Vect()))/(((P4_EL-P4_ELP_reg).Vect()).Mag())/((P4_PIm_reg.Vect()).Mag()));
th_ph_pip = 180./M_PI*acos((P4_PIp_reg.Vect().Dot((P4_EL-P4_ELP_reg).Vect()))/(((P4_EL-P4_ELP_reg).Vect()).Mag())/((P4_PIp_reg.Vect()).Mag()));
th_ph_pr = 180./M_PI*acos((P4_PP_reg.Vect().Dot((P4_EL-P4_ELP_reg).Vect()))/(((P4_EL-P4_ELP_reg).Vect()).Mag())/((P4_PP_reg.Vect()).Mag()));

//cut on missing energy of 0
//if ((P4_miss_0[3] >  -0.05)) {

	h_z_P_sim ->Fill(z_P,sigma);
	h_z_PIp_sim->Fill(z_PIp,sigma);
	h_z_PIm_sim->Fill(z_PIm,sigma);
 
	if ((W>1.3)&&(W<1.8)&&(Q2>0.45)&&(Q2<1.)){
	
	if (((P4_miss_0.Vect()).Mag()<0.2)) {
	
	if ((P4_miss_0.Mag2()>-0.02)&&(P4_miss_0.Mag2()<0.001)){
	h_pim_mis_all_reg_sim[int((W-1.3)/0.1)]-> Fill(P4_PIm_miss.Mag2(),sigma);
	h_miss_en_0_sim->Fill(P4_miss_0[3],sigma);
	};
	h_0_mis_all_reg_sim[int((W-1.3)/0.1)]-> Fill(P4_miss_0.Mag2(),sigma);
	
	};
	h_mom_all_reg_sim[int((W-1.3)/0.1)]-> Fill((P4_miss_0.Vect()).Mag(),sigma);		
	};

/*
if ((W > 1.3)&&(W < 1.8)&&(Q2>0.4)&&(Q2<0.6)) {
h_pim_mis_fermi_nocut_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.1)]->Fill(P4_PIm_miss.Mag2(),sigma);
h_mis_mom_fermi_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.1)]->Fill((P4_miss_0.Vect()).Mag(),sigma);
};
if ((W > 1.3)&&(W < 1.775)&&(Q2>0.6)&&(Q2<0.7)) {
h_pim_mis_fermi_nocut_sim_2[int((Q2-0.6)/0.05)][int((W-1.3)/0.1)]->Fill(P4_PIm_miss.Mag2(),sigma);
h_mis_mom_fermi_sim_2[int((Q2-0.6)/0.05)][int((W-1.3)/0.1)]->Fill((P4_miss_0.Vect()).Mag(),sigma);
};
if ((W > 1.3)&&(W < 1.7)&&(Q2>0.7)&&(Q2<0.8)) {
h_pim_mis_fermi_nocut_sim_3[int((Q2-0.7)/0.05)][int((W-1.3)/0.1)]->Fill(P4_PIm_miss.Mag2(),sigma);
h_mis_mom_fermi_sim_3[int((Q2-0.7)/0.05)][int((W-1.3)/0.1)]->Fill((P4_miss_0.Vect()).Mag(),sigma);
};
if ((W > 1.3)&&(W < 1.65)&&(Q2>0.8)&&(Q2<0.85)) {
h_pim_mis_fermi_nocut_sim_4[int((W-1.3)/0.1)]->Fill(P4_PIm_miss.Mag2(),sigma);
h_mis_mom_fermi_sim_4[int((W-1.3)/0.1)]->Fill((P4_miss_0.Vect()).Mag(),sigma);
};
if ((W > 1.3)&&(W < 1.6)&&(Q2>0.85)&&(Q2<0.9)) {
h_pim_mis_fermi_nocut_sim_5[int((W-1.3)/0.1)]->Fill(P4_PIm_miss.Mag2(),sigma);
h_mis_mom_fermi_sim_5[int((W-1.3)/0.1)]->Fill((P4_miss_0.Vect()).Mag(),sigma);
};
if ((W > 1.3)&&(W < 1.575)&&(Q2>0.9)&&(Q2<0.95)) {
h_pim_mis_fermi_nocut_sim_6[int((W-1.3)/0.1)]->Fill(P4_PIm_miss.Mag2(),sigma);
h_mis_mom_fermi_sim_6[int((W-1.3)/0.1)]->Fill((P4_miss_0.Vect()).Mag(),sigma);
};
if ((W > 1.3)&&(W < 1.55)&&(Q2>0.95)&&(Q2<1.0)) {
h_pim_mis_fermi_nocut_sim_7[int((W-1.3)/0.1)]->Fill(P4_PIm_miss.Mag2(),sigma);
h_mis_mom_fermi_sim_7[int((W-1.3)/0.1)]->Fill((P4_miss_0.Vect()).Mag(),sigma);
};
*/

	h_miss_mom_0_nocut_sim->Fill((P4_miss_0.Vect()).Mag(),1.);

	if ((P4_PIm_miss.Mag2()>0)&&(P4_PIm_miss.Mag2()<0.05)){
	h_miss_mom_0_cut_onpim_sim->Fill((P4_miss_0.Vect()).Mag(),1.);
	};

	hist_PIm_miss_all_reg_1_sim-> Fill(P4_PIm_miss.Mag2(),sigma);

//cut on missing momentum
if ((P4_miss_0.Vect()).Mag() <0.2){

//cut on missing mass of 0
if ((P4_miss_0.Mag2()>-0.02)&&(P4_miss_0.Mag2()<0.001)){ 

//cut on missing mass of pim
if ((P4_PIm_miss.Mag2()>-0.15)&&(P4_PIm_miss.Mag2()<0.15)){

	hist_PIm_miss_all_reg_2_sim-> Fill(P4_PIm_miss.Mag2(),sigma);
	h_miss_mass_0_sim-> Fill(P4_miss_0.Mag2(),1.);
	h_miss_mom_0_cut_on0_sim->Fill((P4_miss_0.Vect()).Mag(),1.);


selection_0_miss_sim = true;
P4_PIm_reg = P4_PIm_miss;

//};
};
};
};
};
};
};
};
};
};
};
};
};

//if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)&&(npart>=3)){
if ((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)&&(npart>=3)){
if ((bool_pip_id_sim)&&(bool_proton_id_sim)&&(n_PIm == 0)){


if (sqrt((P4_PIp_reg+P4_PIm_miss)*(P4_PIp_reg+P4_PIm_miss)) > 2*m_pip) {
if (sqrt((P4_PIp_reg+P4_PIm_miss)*(P4_PIp_reg+P4_PIm_miss)) < W - m_proton) {
if (sqrt((P4_PIp_reg+P4_PP_reg)*(P4_PIp_reg+P4_PP_reg)) > m_pip + m_proton) {
if (sqrt((P4_PIp_reg+P4_PP_reg)*(P4_PIp_reg+P4_PP_reg)) < W - m_pip) {
if (sqrt((P4_PIm_miss+P4_PP_reg)*(P4_PIm_miss+P4_PP_reg)) > m_pip + m_proton) {
if (sqrt((P4_PIm_miss+P4_PP_reg)*(P4_PIm_miss+P4_PP_reg)) < W - m_pip) {

h_mm_pim_vs_npart_sim->Fill(npart,P4_PIm_miss.Mag2(),1.);
h_PIm_miss_en_sim->Fill(P4_PIm_miss[3],1.);

//if (P4_PIm_miss_en_comp[3] > m_pip  - 0.05) {
th_ph_pip = 180./M_PI*acos((P4_PIp_reg.Vect().Dot((P4_EL-P4_ELP_reg).Vect()))/(((P4_EL-P4_ELP_reg).Vect()).Mag())/((P4_PIp_reg.Vect()).Mag()));
th_ph_pr = 180./M_PI*acos((P4_PP_reg.Vect().Dot((P4_EL-P4_ELP_reg).Vect()))/(((P4_EL-P4_ELP_reg).Vect()).Mag())/((P4_PP_reg.Vect()).Mag()));
if (P4_PIm_miss[3] > m_pip) {

	if ((W>1.3)&&(W<1.8)&&(Q2>0.45)&&(Q2<0.5)){
		
	h_pim_mis_main_top_sim[int((W-1.3)/0.1)]-> Fill(P4_PIm_miss.Mag2(),1.);
		
	};
	
//if ((W>1.3)&&(W<1.8)&&(th_P<40.)&&(th_PIp<90.)&&(P_P>0.4)&&(th_ph_pip<85.)&&(th_ph_pr<35.)){
//if ((z_P>-2.5)&&(z_P<2.)&&(z_PIp>-2.5)&&(z_PIp<2.)&&((z_P-z_EL)>-2.5)&&((z_P-z_EL)<2.5)&&((z_PIp-z_EL)>-2.5)&&((z_PIp-z_EL)<2.5)&&((z_P-z_PIp)>-2.5)&&((z_P-z_PIp)<2.5)){
//if ((W<1.6)||((W>1.6)&&(P4_PIm_miss.Mag()<0.25))){
//if ((th_PIp>10.)&&(th_PIp<130.)&&(th_P>10.)&&(th_P<50.)&&(P_P > 0.25)&&(P_P < 1.5)){
//if (((W>1.3)&&(W<1.5))||((W>1.5)&&(W<1.825)&&(ARR_MMC_PIM_MISS[int((P_P-0.25)/0.25)][int((th_PIp-10.)/15.)][int((th_P-10.)/5.)]==1))){
/*
if ((W > 1.3)&&(W < 1.825)&&(Q2>0.4)&&(Q2<0.6)) {
h_pim_mis_fermi_nocut_sim_1[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(sqrt(abs(P4_PIm_miss.Mag2())),sigma);
};

if ((W > 1.3)&&(W < 1.775)&&(Q2>0.6)&&(Q2<0.7)) {
h_pim_mis_fermi_nocut_sim_2[int((Q2-0.6)/0.05)][int((W-1.3)/0.025)]->Fill(sqrt(abs(P4_PIm_miss.Mag2())),sigma);
};

if ((W > 1.3)&&(W < 1.725)&&(Q2>0.7)&&(Q2<0.8)) {
h_pim_mis_fermi_nocut_sim_3[int((Q2-0.7)/0.05)][int((W-1.3)/0.025)]->Fill(sqrt(abs(P4_PIm_miss.Mag2())),sigma);
};

if ((W > 1.3)&&(W < 1.65)&&(Q2>0.8)&&(Q2<0.85)) {
h_pim_mis_fermi_nocut_sim_4[int((W-1.3)/0.025)]->Fill(sqrt(abs(P4_PIm_miss.Mag2())),sigma);
};

if ((W > 1.3)&&(W < 1.625)&&(Q2>0.85)&&(Q2<0.9)) {
h_pim_mis_fermi_nocut_sim_5[int((W-1.3)/0.025)]->Fill(sqrt(abs(P4_PIm_miss.Mag2())),sigma);
};

if ((W > 1.3)&&(W < 1.6)&&(Q2>0.9)&&(Q2<0.95)) {
h_pim_mis_fermi_nocut_sim_6[int((W-1.3)/0.025)]->Fill(sqrt(abs(P4_PIm_miss.Mag2())),sigma);
};

if ((W > 1.3)&&(W < 1.55)&&(Q2>0.95)&&(Q2<1.0)) {
h_pim_mis_fermi_nocut_sim_7[int((W-1.3)/0.025)]->Fill(sqrt(abs(P4_PIm_miss.Mag2())),sigma);
};*/

//};
//};
//};


	h_PIm_miss_sim-> Fill(P4_PIm_miss.Mag2(),sigma);


if (sqrt(abs(P4_PIm_miss.Mag2())) < MMcut_pim_miss[int((W-1.3)/0.025)]){
//if ((P4_PIm_miss.Mag2()>-0.0272)&&(P4_PIm_miss.Mag2()<0.068)){ 
W_2pi_selection_sim->Fill(W,Q2,sigma);
selection_pim_miss_sim= true;
P4_PIm_reg = P4_PIm_miss;
};
};
};
};
};
};
};
};
};
};

//};

};//vertex difference cut
};//konets ifa electronnix cutov


if (indtype==2) Pgen=P_EL;
if (indtype==1) Prec=P_EL;
//cout << Pgen<<" "<<Prec<<"\n";
if (indtype==1) h_sim_mom_corr_test->Fill(Pgen-Prec,1.);

if ((selection_pim_miss_sim)||(selection_0_miss_sim)) {



if((W > 1.3)&&(W < 1.825)&&(Q2 > 0.4)&&(Q2 < 1.)){
if (indtype==1){
rot_boost_cmsyst();
Var_1[0] = inv_m_pip_p;
Var_1[1] = inv_m_pip_pim;
Var_1[2] = theta_P_cm;
Var_1[3] = phi_P_cm;
Var_1[4] = alpha_PIpPIm_pipf;

Var_2[0] = inv_m_pip_p;
Var_2[1] = inv_m_pip_pim;
Var_2[2] = theta_PIm_cm;
Var_2[3] = phi_PIm_cm;
Var_2[4] = alpha_PPIp_piPIm;

Var_3[0] = inv_m_pim_p;
Var_3[1] = inv_m_pip_pim;
Var_3[2] = theta_PIp_cm;
Var_3[3] = phi_PIp_cm;
Var_3[4] = alpha_PPIm_piPIp;



if ((selection_pim_miss_sim)){


h_5dim_pim_1_sim[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_1,sigma);
h_5dim_pim_2_sim[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_2,sigma);
h_5dim_pim_3_sim[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_3,sigma); 

h_5dim_pim_1_sim_evt[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_1,sigma*sigma);
h_5dim_pim_2_sim_evt[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_2,sigma*sigma);
h_5dim_pim_3_sim_evt[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_3,sigma*sigma); 
};
if ((selection_0_miss_sim)){

h_5dim_excl_1_sim[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_1,sigma);
h_5dim_excl_2_sim[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_2,sigma);
h_5dim_excl_3_sim[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_3,sigma);

h_5dim_excl_1_sim_evt[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_1,sigma*sigma);
h_5dim_excl_2_sim_evt[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_2,sigma*sigma);
h_5dim_excl_3_sim_evt[int((Q2-0.4)/0.05)][int((W-1.3)/0.025)]->Fill(Var_3,sigma*sigma);

};
};
};





}; //selection


 //}; //konec ifa electronnih cutov
    };// konets tsicla po sobitiyam (po i)
  
  
  t21->Delete();

 finp->Close();  
 
     };//konets tsicla po failam (po m)  
     
   /*    
  for (i=0; i<9;i++) {      
//ERRORS

for (k=1; k<=h_5dim_excl_1[i]->GetNbins(); k++) {
h_5dim_excl_1[i]->SetBinError(k,sqrt(h_5dim_excl_1[i]->GetBinContent(k)));
};
for (k=1; k<=h_5dim_excl_2[i]->GetNbins(); k++) {
h_5dim_excl_2[i]->SetBinError(k,sqrt(h_5dim_excl_2[i]->GetBinContent(k)));
};
for (k=1; k<=h_5dim_excl_3[i]->GetNbins(); k++) {
h_5dim_excl_3[i]->SetBinError(k,sqrt(h_5dim_excl_3[i]->GetBinContent(k)));
};

for (k=1; k<=h_5dim_excl_1_empty[i]->GetNbins(); k++) {
h_5dim_excl_1_empty[i]->SetBinError(k,sqrt(h_5dim_excl_1_empty[i]->GetBinContent(k)));
};
for (k=1; k<=h_5dim_excl_2_empty[i]->GetNbins(); k++) {
h_5dim_excl_2_empty[i]->SetBinError(k,sqrt(h_5dim_excl_2_empty[i]->GetBinContent(k)));
};
for (k=1; k<=h_5dim_excl_3_empty[i]->GetNbins(); k++) {
h_5dim_excl_3_empty[i]->SetBinError(k,sqrt(h_5dim_excl_3_empty[i]->GetBinContent(k)));
};
    
  };  
    
  */
     
     outFile = new TFile(outfile_inp.c_str(),"recreate");

h_cos_th = new TH1D("h_cos_th","h_cos_th",10,0.,180.);
     
Double_t temp;    
for (j=1; j<=10; j++) {
//temp = cos((h_cos_th->GetBinLowEdge(j))*M_PI/180.)-cos(M_PI/180.*(h_cos_th->GetBinLowEdge(j)+h_cos_th->GetBinWidth(j)));
//temp = cos(((h_cos_th->GetBinLowEdge(j))+h_cos_th->GetBinWidth(j)/2.)*M_PI/180.);
temp = sin(((h_cos_th->GetBinLowEdge(j))+h_cos_th->GetBinWidth(j)/2.)*M_PI/180.);
h_cos_th->SetBinContent(j,temp);
//h_cos_th->SetBinError(j,0.);
};     
     
    
       
//VIRTUAL PHOTON FLUX CALCULATION
Double_t W_bin[21];
Double_t Q2_bin[12];
Double_t omega[12][21],en_elp[12][21],th_elp[12][21],epsilon[12][21],flux[12][21],factor,L_0;
Long64_t Entr_sim_rec;
THnSparseD *ha[12][21];

for (j=0; j<12; j++) {
Q2_bin[j] = 0.425+0.05*j;
 for (i=0; i<21;i++) {
 W_bin[i] = 1.3125+0.025*i;
omega[j][i] = (W_bin[i]*W_bin[i] + Q2_bin[j] - m_proton*m_proton)/2./m_proton ;
en_elp[j][i] = 2.039 - omega[j][i];
th_elp[j][i]  = 2*asin(sqrt(Q2_bin[j]/4./2.039/en_elp[j][i]));

epsilon[j][i] = 1/(1. + 2.*(1. + omega[j][i]*omega[j][i]/Q2_bin[j])*(tan(th_elp[j][i]/2.))*(tan(th_elp[j][i]/2.)));
flux[j][i] = (omega[j][i]-Q2_bin[j]/2./m_proton)/137.;

flux[j][i]= flux[j][i]/2./(M_PI)/2.039/Q2_bin[j]/(1-epsilon[j][i]);
flux[j][i] = flux[j][i]*W_bin[i]/2.039/m_proton; 



//ERRORS

/*
for (k=1; k<=h_5dim_1[j][i]->GetNbins(); k++) {
h_5dim_1[j][i]->SetBinError(k,sqrt(h_5dim_1[j][i]->GetBinContent(k)));
};
for (k=1; k<=h_5dim_2[j][i]->GetNbins(); k++) {
h_5dim_2[j][i]->SetBinError(k,sqrt(h_5dim_2[j][i]->GetBinContent(k)));
};
for (k=1; k<=h_5dim_3[j][i]->GetNbins(); k++) {
h_5dim_3[j][i]->SetBinError(k,sqrt(h_5dim_3[j][i]->GetBinContent(k)));
};

for (k=1; k<=h_5dim_1_empty[j][i]->GetNbins(); k++) {
h_5dim_1_empty[j][i]->SetBinError(k,sqrt(h_5dim_1_empty[j][i]->GetBinContent(k)));
};
for (k=1; k<=h_5dim_2_empty[j][i]->GetNbins(); k++) {
h_5dim_2_empty[j][i]->SetBinError(k,sqrt(h_5dim_2_empty[j][i]->GetBinContent(k)));
};
for (k=1; k<=h_5dim_3_empty[j][i]->GetNbins(); k++) {
h_5dim_3_empty[j][i]->SetBinError(k,sqrt(h_5dim_3_empty[j][i]->GetBinContent(k)));
};*/


/*
for (k=1; k<=h_5dim_1_sim_1[j][i]->GetNbins(); k++) {
h_5dim_1_sim_1[j][i]->SetBinError(k,sqrt(h_5dim_1_sim_1[j][i]->GetBinContent(k)));
};
for (k=1; k<=h_5dim_2_sim_1[j][i]->GetNbins(); k++) {
h_5dim_2_sim_1[j][i]->SetBinError(k,sqrt(h_5dim_2_sim_1[j][i]->GetBinContent(k)));
};
for (k=1; k<=h_5dim_3_sim_1[j][i]->GetNbins(); k++) {
h_5dim_3_sim_1[j][i]->SetBinError(k,sqrt(h_5dim_3_sim_1[j][i]->GetBinContent(k)));
};

for (k=1; k<=h_5dim_1_sim_2[j][i]->GetNbins(); k++) {
h_5dim_1_sim_2[j][i]->SetBinError(k,sqrt(h_5dim_1_sim_2[j][i]->GetBinContent(k)));
};
for (k=1; k<=h_5dim_2_sim_2[j][i]->GetNbins(); k++) {
h_5dim_2_sim_2[j][i]->SetBinError(k,sqrt(h_5dim_2_sim_2[j][i]->GetBinContent(k)));
};
for (k=1; k<=h_5dim_3_sim_2[j][i]->GetNbins(); k++) {
h_5dim_3_sim_2[j][i]->SetBinError(k,sqrt(h_5dim_3_sim_2[j][i]->GetBinContent(k)));
};
*/

//cout << "j = " << j << " i = " << i << " flux = " << flux[j][i] << "\n";
factor = 1./flux[j][i];
L_0 = 0.63E12;
//factor = factor/L_0;


//h_5dim_1[j][i]->Scale(factor);
//h_5dim_2[j][i]->Scale(factor);
//h_5dim_3[j][i]->Scale(factor);

//h_5dim_1_empty[j][i]->Scale(factor);
//h_5dim_2_empty[j][i]->Scale(factor);
//h_5dim_3_empty[j][i]->Scale(factor);


 };
  };


gStyle->SetTitleSize(0.08,"t");
gStyle->SetOptStat("e");
gStyle->SetStatY(0.88); 
gStyle->SetTitleY(0.99);
gStyle->SetTitleX(0.445);
 gStyle->SetErrorX(0);
//gStyle->SetError(0); 
/*
 for (j=0; j<12; j++) {
 for (i=0; i<21;i++) {


h1prj_inv_m_pip_p[j][i] = h_5dim_1[j][i]->Projection(0,"0E"); 
h1prj_inv_m_pip_p[j][i]->SetMarkerStyle(20);
h1prj_inv_m_pip_p[j][i]->SetOption("e1pX0");
qqq << "m_{inv. #pi^{+}p}, GeV ";
h1prj_inv_m_pip_p[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h1prj_inv_m_pip_p_1_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h1prj_inv_m_pip_p[j][i]->SetName(qqq.str().c_str());




qqq.str("");

h1prj_inv_m_pip_pim[j][i] = h_5dim_1[j][i]->Projection(1,"0E"); 
h1prj_inv_m_pip_pim[j][i]->SetMarkerStyle(20);
h1prj_inv_m_pip_pim[j][i]->SetOption("e1pX0");
qqq << "m_{inv. #pi^{+}#pi^{-}}, GeV ";
h1prj_inv_m_pip_pim[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h1prj_inv_m_pip_pim_1_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h1prj_inv_m_pip_pim[j][i]->SetName(qqq.str().c_str());
qqq.str("");


h1prj_th_P[j][i] = h_5dim_1[j][i]->Projection(2,"0E"); 
h1prj_th_P[j][i] ->Divide(h_cos_th);
h1prj_th_P[j][i] ->SetEntries(h1prj_inv_m_pip_pim[j][i]->GetEntries());
h1prj_th_P[j][i]->SetMarkerStyle(20);
h1prj_th_P[j][i]->SetOption("e1pX0");
h1prj_th_P[j][i]->SetOption("X0");
qqq << "#theta_{p'} in c.m., deg ";
h1prj_th_P[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h1prj_th_P_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h1prj_th_P[j][i]->SetName(qqq.str().c_str());
qqq.str("");


h1prj_phi_P[j][i] = h_5dim_1[j][i]->Projection(3,"0E"); 
h1prj_phi_P[j][i]->SetMarkerStyle(20);
h1prj_phi_P[j][i]->SetOption("e1pX0");
qqq << "#phi_{p'} in c.m., deg";
h1prj_phi_P[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h1prj_phi_P_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h1prj_phi_P[j][i]->SetName(qqq.str().c_str());
qqq.str("");

h1prj_alpha_PIpPIm_pipf[j][i] = h_5dim_1[j][i]->Projection(4,"0E");
h1prj_alpha_PIpPIm_pipf[j][i]->SetMarkerStyle(20);
h1prj_alpha_PIpPIm_pipf[j][i]->SetOption("e1pX0");
qqq << "#alpha_{#pi^{+}#pi^{-}_pp'}, deg ";
h1prj_alpha_PIpPIm_pipf[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h1prj_alpha_PIpPIm_pipf_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h1prj_alpha_PIpPIm_pipf[j][i]->SetName(qqq.str().c_str());
qqq.str("");

/////////
h2prj_inv_m_pip_p[j][i] = h_5dim_2[j][i]->Projection(0,"0E"); 
h2prj_inv_m_pip_p[j][i]->SetMarkerStyle(20);
h2prj_inv_m_pip_p[j][i]->SetOption("e1pX0");
qqq << "m_{inv. #pi^{+}p}, GeV ";
h2prj_inv_m_pip_p[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h2prj_inv_m_pip_p_2_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h2prj_inv_m_pip_p[j][i]->SetName(qqq.str().c_str());
qqq.str("");

h2prj_inv_m_pip_pim[j][i] = h_5dim_2[j][i]->Projection(1,"0E"); 
h2prj_inv_m_pip_pim[j][i]->SetMarkerStyle(20);
h2prj_inv_m_pip_pim[j][i]->SetOption("e1pX0");
qqq << "m_{inv. #pi^{+}#pi^{-}}, GeV ";
h2prj_inv_m_pip_pim[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h2prj_inv_m_pip_pim_2_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h2prj_inv_m_pip_pim[j][i]->SetName(qqq.str().c_str());
qqq.str("");


h2prj_th_PIm[j][i]= h_5dim_2[j][i]->Projection(2,"0E");
h2prj_th_PIm[j][i] ->Divide(h_cos_th);
h2prj_th_PIm[j][i]->SetEntries(h1prj_inv_m_pip_pim[j][i]->GetEntries());
h2prj_th_PIm[j][i]->SetMarkerStyle(20); 
h2prj_th_PIm[j][i]->SetOption("e1pX0");
qqq << "#theta_{#pi^{-}} in c.m., deg";
h2prj_th_PIm[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h2prj_th_PIm_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h2prj_th_PIm[j][i]->SetName(qqq.str().c_str());
qqq.str("");


h2prj_phi_PIm[j][i]= h_5dim_2[j][i]->Projection(3,"0E"); 
h2prj_phi_PIm[j][i]->SetMarkerStyle(20);
h2prj_phi_PIm[j][i]->SetOption("e1pX0");
qqq << "#phi_{#pi^{-}} in c.m., deg";
h2prj_phi_PIm[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h2prj_phi_PIm_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h2prj_phi_PIm[j][i]->SetName(qqq.str().c_str());
qqq.str("");

h2prj_alpha_PPIp_piPIm[j][i]= h_5dim_2[j][i]->Projection(4,"0E");
h2prj_alpha_PPIp_piPIm[j][i]->SetMarkerStyle(20);
h2prj_alpha_PPIp_piPIm[j][i]->SetOption("e1pX0");
qqq << "#alpha_{p'#pi^{+}_p#pi^{-}}, deg";
h2prj_alpha_PPIp_piPIm[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h2prj_alpha_PPIp_piPIm_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h2prj_alpha_PPIp_piPIm[j][i]->SetName(qqq.str().c_str());
qqq.str("");


//////////
h3prj_inv_m_pim_p[j][i] = h_5dim_3[j][i]->Projection(0,"0E"); 
h3prj_inv_m_pim_p[j][i]->SetMarkerStyle(20);
h3prj_inv_m_pim_p[j][i]->SetOption("e1pX0");
qqq << "m_{inv. #pi^{-}p}, GeV ";
h3prj_inv_m_pim_p[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h3prj_inv_m_pim_p_3_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h3prj_inv_m_pim_p[j][i]->SetName(qqq.str().c_str());
qqq.str("");

h3prj_inv_m_pip_pim[j][i] = h_5dim_3[j][i]->Projection(1,"0E"); 
h3prj_inv_m_pip_pim[j][i]->SetMarkerStyle(20);
h3prj_inv_m_pip_pim[j][i]->SetOption("e1pX0");
qqq << "m_{inv. #pi^{+}#pi^{-}}, GeV ";
h3prj_inv_m_pip_pim[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h3prj_inv_m_pip_pim_2"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h3prj_inv_m_pip_pim[j][i]->SetName(qqq.str().c_str());
qqq.str("");

h3prj_th_PIp[j][i] = h_5dim_3[j][i]->Projection(2,"0E"); 
h3prj_th_PIp[j][i] ->Divide(h_cos_th);
h3prj_th_PIp[j][i] ->SetEntries(h1prj_inv_m_pip_pim[j][i]->GetEntries());
h3prj_th_PIp[j][i]->SetMarkerStyle(20);
h3prj_th_PIp[j][i]->SetOption("e1pX0");
qqq << "#theta_{#pi^{+}} in c.m., deg ";
h3prj_th_PIp[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h3prj_th_PIp_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h3prj_th_PIp[j][i]->SetName(qqq.str().c_str());
qqq.str("");

h3prj_phi_PIp[j][i] = h_5dim_3[j][i]->Projection(3,"0E"); 
h3prj_phi_PIp[j][i]->SetMarkerStyle(20);
h3prj_phi_PIp[j][i]->SetOption("e1pX0");
qqq << "#phi_{#pi^{+}} in c.m., deg";
h3prj_phi_PIp[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h3prj_phi_PIp_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h3prj_phi_PIp[j][i]->SetName(qqq.str().c_str());
qqq.str("");

h3prj_alpha_PPIm_piPIp[j][i] = h_5dim_3[j][i]->Projection(4,"0E");
h3prj_alpha_PPIm_piPIp[j][i]->SetMarkerStyle(20);
h3prj_alpha_PPIm_piPIp[j][i]->SetOption("e1pX0");
qqq << "#alpha_{p'#pi^{-}_p#pi^{+}}, deg";
h3prj_alpha_PPIm_piPIp[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h3prj_alpha_PPIm_piPIp_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h3prj_alpha_PPIm_piPIp[j][i]->SetName(qqq.str().c_str());
qqq.str("");



////////////////////////sim



h1prj_inv_m_pip_p_sim[j][i] = h_5dim_1_sim_1[j][i]->Projection(0,"0E"); 
h1prj_inv_m_pip_p_sim[j][i]->SetMarkerStyle(20);
h1prj_inv_m_pip_p_sim[j][i]->SetOption("e1pX0");
qqq << "m_{inv. #pi^{+}p}, GeV, sim rec";
h1prj_inv_m_pip_p_sim[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h1prj_inv_m_pip_p_1_sim_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h1prj_inv_m_pip_p_sim[j][i]->SetName(qqq.str().c_str());




qqq.str("");

h1prj_inv_m_pip_pim_sim[j][i] = h_5dim_1_sim_1[j][i]->Projection(1,"0E"); 
h1prj_inv_m_pip_pim_sim[j][i]->SetMarkerStyle(20);
h1prj_inv_m_pip_pim_sim[j][i]->SetOption("e1pX0");
qqq << "m_{inv. #pi^{+}#pi^{-}}, GeV, sim rec";
h1prj_inv_m_pip_pim_sim[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h1prj_inv_m_pip_pim_1_sim_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h1prj_inv_m_pip_pim_sim[j][i]->SetName(qqq.str().c_str());
qqq.str("");

h1prj_th_P_sim[j][i] = h_5dim_1_sim_1[j][i]->Projection(2,"0E"); 
h1prj_th_P_sim[j][i] ->Divide(h_cos_th);
h1prj_th_P_sim[j][i] ->SetEntries(h1prj_inv_m_pip_pim_sim[j][i]->GetEntries());
h1prj_th_P_sim[j][i]->SetMarkerStyle(20);
h1prj_th_P_sim[j][i]->SetOption("e1pX0");
h1prj_th_P_sim[j][i]->SetOption("X0");
qqq << "#theta_{p'} in c.m., deg, sim rec ";
h1prj_th_P_sim[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h1prj_th_P_sim_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h1prj_th_P_sim[j][i]->SetName(qqq.str().c_str());
qqq.str("");

h1prj_phi_P_sim[j][i] = h_5dim_1_sim_1[j][i]->Projection(3,"0E"); 
h1prj_phi_P_sim[j][i]->SetMarkerStyle(20);
h1prj_phi_P_sim[j][i]->SetOption("e1pX0");
qqq << "#phi_{p'} in c.m., deg, sim rec";
h1prj_phi_P_sim[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h1prj_phi_P_sim_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h1prj_phi_P_sim[j][i]->SetName(qqq.str().c_str());
qqq.str("");

h1prj_alpha_PIpPIm_pipf_sim[j][i] = h_5dim_1_sim_1[j][i]->Projection(4,"0E");
h1prj_alpha_PIpPIm_pipf_sim[j][i]->SetMarkerStyle(20);
h1prj_alpha_PIpPIm_pipf_sim[j][i]->SetOption("e1pX0");
qqq << "#alpha_{#pi^{+}#pi^{-}_pp'}, deg, sim rec ";
h1prj_alpha_PIpPIm_pipf_sim[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h1prj_alpha_PIpPIm_pipf_sim_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h1prj_alpha_PIpPIm_pipf_sim[j][i]->SetName(qqq.str().c_str());
qqq.str("");

/////////
h2prj_inv_m_pip_p_sim[j][i] = h_5dim_2_sim_1[j][i]->Projection(0,"0E"); 
h2prj_inv_m_pip_p_sim[j][i]->SetMarkerStyle(20);
h2prj_inv_m_pip_p_sim[j][i]->SetOption("e1pX0");
qqq << "m_{inv. #pi^{+}p}, GeV, sim rec ";
h2prj_inv_m_pip_p_sim[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h2prj_inv_m_pip_p_2_sim_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h2prj_inv_m_pip_p_sim[j][i]->SetName(qqq.str().c_str());
qqq.str("");

h2prj_inv_m_pip_pim_sim[j][i] = h_5dim_2_sim_1[j][i]->Projection(1,"0E"); 
h2prj_inv_m_pip_pim_sim[j][i]->SetMarkerStyle(20);
h2prj_inv_m_pip_pim_sim[j][i]->SetOption("e1pX0");
qqq << "m_{inv. #pi^{+}#pi^{-}}, GeV, sim rec ";
h2prj_inv_m_pip_pim_sim[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h2prj_inv_m_pip_pim_2_sim_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h2prj_inv_m_pip_pim_sim[j][i]->SetName(qqq.str().c_str());
qqq.str("");


h2prj_th_PIm_sim[j][i]= h_5dim_2_sim_1[j][i]->Projection(2,"0E");
h2prj_th_PIm_sim[j][i] ->Divide(h_cos_th);
h2prj_th_PIm_sim[j][i]->SetEntries(h1prj_inv_m_pip_pim_sim[j][i]->GetEntries());
h2prj_th_PIm_sim[j][i]->SetMarkerStyle(20); 
h2prj_th_PIm_sim[j][i]->SetOption("e1pX0");
qqq << "#theta_{#pi^{-}} in c.m., deg, sim rec";
h2prj_th_PIm_sim[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h2prj_th_PIm_sim_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h2prj_th_PIm_sim[j][i]->SetName(qqq.str().c_str());
qqq.str("");


h2prj_phi_PIm_sim[j][i]= h_5dim_2_sim_1[j][i]->Projection(3,"0E"); 
h2prj_phi_PIm_sim[j][i]->SetMarkerStyle(20);
h2prj_phi_PIm_sim[j][i]->SetOption("e1pX0");
qqq << "#phi_{#pi^{-}} in c.m., deg, sim rec";
h2prj_phi_PIm_sim[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h2prj_phi_PIm_sim_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h2prj_phi_PIm_sim[j][i]->SetName(qqq.str().c_str());
qqq.str("");

h2prj_alpha_PPIp_piPIm_sim[j][i]= h_5dim_2_sim_1[j][i]->Projection(4,"0E");
h2prj_alpha_PPIp_piPIm_sim[j][i]->SetMarkerStyle(20);
h2prj_alpha_PPIp_piPIm_sim[j][i]->SetOption("e1pX0");
qqq << "#alpha_{p'#pi^{+}_p#pi^{-}}, deg, sim rec";
h2prj_alpha_PPIp_piPIm_sim[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h2prj_alpha_PPIp_piPIm_sim_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h2prj_alpha_PPIp_piPIm_sim[j][i]->SetName(qqq.str().c_str());
qqq.str("");


//////////
h3prj_inv_m_pim_p_sim[j][i] = h_5dim_3_sim_1[j][i]->Projection(0,"0E"); 
h3prj_inv_m_pim_p_sim[j][i]->SetMarkerStyle(20);
h3prj_inv_m_pim_p_sim[j][i]->SetOption("e1pX0");
qqq << "m_{inv. #pi^{-}p}, GeV, sim rec ";
h3prj_inv_m_pim_p_sim[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h3prj_inv_m_pim_p_3_sim_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h3prj_inv_m_pim_p_sim[j][i]->SetName(qqq.str().c_str());
qqq.str("");

h3prj_inv_m_pip_pim_sim[j][i] = h_5dim_3_sim_1[j][i]->Projection(1,"0E"); 
h3prj_inv_m_pip_pim_sim[j][i]->SetMarkerStyle(20);
h3prj_inv_m_pip_pim_sim[j][i]->SetOption("e1pX0");
qqq << "m_{inv. #pi^{+}#pi^{-}}, GeV, sim rec ";
h3prj_inv_m_pip_pim_sim[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h3prj_inv_m_pip_pim_2_sim_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h3prj_inv_m_pip_pim_sim[j][i]->SetName(qqq.str().c_str());
qqq.str("");

h3prj_th_PIp_sim[j][i] = h_5dim_3_sim_1[j][i]->Projection(2,"0E"); 
h3prj_th_PIp_sim[j][i] ->Divide(h_cos_th);
h3prj_th_PIp_sim[j][i] ->SetEntries(h1prj_inv_m_pip_pim_sim[j][i]->GetEntries());
h3prj_th_PIp_sim[j][i]->SetMarkerStyle(20);
h3prj_th_PIp_sim[j][i]->SetOption("e1pX0");
qqq << "#theta_{#pi^{+}} in c.m., deg, sim rec ";
h3prj_th_PIp_sim[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h3prj_th_PIp_sim_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h3prj_th_PIp_sim[j][i]->SetName(qqq.str().c_str());
qqq.str("");

h3prj_phi_PIp_sim[j][i] = h_5dim_3_sim_1[j][i]->Projection(3,"0E"); 
h3prj_phi_PIp_sim[j][i]->SetMarkerStyle(20);
h3prj_phi_PIp_sim[j][i]->SetOption("e1pX0");
qqq << "#phi_{#pi^{+}} in c.m., deg, sim rec";
h3prj_phi_PIp_sim[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h3prj_phi_PIp_sim_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h3prj_phi_PIp_sim[j][i]->SetName(qqq.str().c_str());
qqq.str("");

h3prj_alpha_PPIm_piPIp_sim[j][i] = h_5dim_3_sim_1[j][i]->Projection(4,"0E");
h3prj_alpha_PPIm_piPIp_sim[j][i]->SetMarkerStyle(20);
h3prj_alpha_PPIm_piPIp_sim[j][i]->SetOption("e1pX0");
qqq << "#alpha_{p'#pi^{-}_p#pi^{+}}, deg, sim rec";
h3prj_alpha_PPIm_piPIp_sim[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h3prj_alpha_PPIm_piPIp_sim_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h3prj_alpha_PPIm_piPIp_sim[j][i]->SetName(qqq.str().c_str());
qqq.str("");

/////////////////////////sim gen


h1prj_inv_m_pip_p_sim_2[j][i] = h_5dim_1_sim_2[j][i]->Projection(0,"0E"); 
h1prj_inv_m_pip_p_sim_2[j][i]->SetMarkerStyle(20);
h1prj_inv_m_pip_p_sim_2[j][i]->SetOption("e1pX0");
qqq << "m_{inv. #pi^{+}p}, GeV, sim gen ";
h1prj_inv_m_pip_p_sim_2[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h1prj_inv_m_pip_p_1_sim_2_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h1prj_inv_m_pip_p_sim_2[j][i]->SetName(qqq.str().c_str());




qqq.str("");

h1prj_inv_m_pip_pim_sim_2[j][i] = h_5dim_1_sim_2[j][i]->Projection(1,"0E"); 
h1prj_inv_m_pip_pim_sim_2[j][i]->SetMarkerStyle(20);
h1prj_inv_m_pip_pim_sim_2[j][i]->SetOption("e1pX0");
qqq << "m_{inv. #pi^{+}#pi^{-}}, GeV, sim gen ";
h1prj_inv_m_pip_pim_sim_2[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h1prj_inv_m_pip_pim_1_sim_2_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h1prj_inv_m_pip_pim_sim_2[j][i]->SetName(qqq.str().c_str());
qqq.str("");

h1prj_th_P_sim_2[j][i] = h_5dim_1_sim_2[j][i]->Projection(2,"0E"); 
//h1prj_th_P_sim_2[j][i] ->Divide(h_cos_th);
//h1prj_th_P_sim_2[j][i] ->SetEntries(h1prj_inv_m_pip_pim_sim_2[j][i]->GetEntries());
h1prj_th_P_sim_2[j][i]->SetMarkerStyle(20);
h1prj_th_P_sim_2[j][i]->SetOption("e1pX0");
h1prj_th_P_sim_2[j][i]->SetOption("X0");
qqq << "#theta_{p'} in c.m., deg, sim gen ";
h1prj_th_P_sim_2[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h1prj_th_P_sim_2_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h1prj_th_P_sim_2[j][i]->SetName(qqq.str().c_str());
qqq.str("");

h1prj_phi_P_sim_2[j][i] = h_5dim_1_sim_2[j][i]->Projection(3,"0E"); 
h1prj_phi_P_sim_2[j][i]->SetMarkerStyle(20);
h1prj_phi_P_sim_2[j][i]->SetOption("e1pX0");
qqq << "#phi_{p'} in c.m., deg, sim gen";
h1prj_phi_P_sim_2[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h1prj_phi_P_sim_2_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h1prj_phi_P_sim_2[j][i]->SetName(qqq.str().c_str());
qqq.str("");

h1prj_alpha_PIpPIm_pipf_sim_2[j][i] = h_5dim_1_sim_2[j][i]->Projection(4,"0E");
h1prj_alpha_PIpPIm_pipf_sim_2[j][i]->SetMarkerStyle(20);
h1prj_alpha_PIpPIm_pipf_sim_2[j][i]->SetOption("e1pX0");
qqq << "#alpha_{#pi^{+}#pi^{-}_pp'}, deg, sim gen ";
h1prj_alpha_PIpPIm_pipf_sim_2[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h1prj_alpha_PIpPIm_pipf_sim_2_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h1prj_alpha_PIpPIm_pipf_sim_2[j][i]->SetName(qqq.str().c_str());
qqq.str("");

/////////
h2prj_inv_m_pip_p_sim_2[j][i] = h_5dim_2_sim_2[j][i]->Projection(0,"0E"); 
h2prj_inv_m_pip_p_sim_2[j][i]->SetMarkerStyle(20);
h2prj_inv_m_pip_p_sim_2[j][i]->SetOption("e1pX0");
qqq << "m_{inv. #pi^{+}p}, GeV, sim gen ";
h2prj_inv_m_pip_p_sim_2[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h2prj_inv_m_pip_p_2_sim_2_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h2prj_inv_m_pip_p_sim_2[j][i]->SetName(qqq.str().c_str());
qqq.str("");

h2prj_inv_m_pip_pim_sim_2[j][i] = h_5dim_2_sim_2[j][i]->Projection(1,"0E"); 
h2prj_inv_m_pip_pim_sim_2[j][i]->SetMarkerStyle(20);
h2prj_inv_m_pip_pim_sim_2[j][i]->SetOption("e1pX0");
qqq << "m_{inv. #pi^{+}#pi^{-}}, GeV, sim gen ";
h2prj_inv_m_pip_pim_sim_2[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h2prj_inv_m_pip_pim_2_sim_2_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h2prj_inv_m_pip_pim_sim_2[j][i]->SetName(qqq.str().c_str());
qqq.str("");


h2prj_th_PIm_sim_2[j][i]= h_5dim_2_sim_2[j][i]->Projection(2,"0E");
//h2prj_th_PIm_sim_2[j][i] ->Divide(h_cos_th);
//h2prj_th_PIm_sim_2[j][i]->SetEntries(h2prj_inv_m_pip_pim_sim_2[j][i]->GetEntries());
h2prj_th_PIm_sim_2[j][i]->SetMarkerStyle(20); 
h2prj_th_PIm_sim_2[j][i]->SetOption("e1pX0");
qqq << "#theta_{#pi^{-}} in c.m., deg, sim gen";
h2prj_th_PIm_sim_2[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h2prj_th_PIm_sim_2_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h2prj_th_PIm_sim_2[j][i]->SetName(qqq.str().c_str());
qqq.str("");


h2prj_phi_PIm_sim_2[j][i]= h_5dim_2_sim_2[j][i]->Projection(3,"0E"); 
h2prj_phi_PIm_sim_2[j][i]->SetMarkerStyle(20);
h2prj_phi_PIm_sim_2[j][i]->SetOption("e1pX0");
qqq << "#phi_{#pi^{-}} in c.m., deg, sim gen";
h2prj_phi_PIm_sim_2[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h2prj_phi_PIm_sim_2_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h2prj_phi_PIm_sim_2[j][i]->SetName(qqq.str().c_str());
qqq.str("");

h2prj_alpha_PPIp_piPIm_sim_2[j][i]= h_5dim_2_sim_2[j][i]->Projection(4,"0E");
h2prj_alpha_PPIp_piPIm_sim_2[j][i]->SetMarkerStyle(20);
h2prj_alpha_PPIp_piPIm_sim_2[j][i]->SetOption("e1pX0");
qqq << "#alpha_{p'#pi^{+}_p#pi^{-}}, deg, sim gen";
h2prj_alpha_PPIp_piPIm_sim_2[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h2prj_alpha_PPIp_piPIm_sim_2_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h2prj_alpha_PPIp_piPIm_sim_2[j][i]->SetName(qqq.str().c_str());
qqq.str("");


//////////
h3prj_inv_m_pim_p_sim_2[j][i] = h_5dim_3_sim_2[j][i]->Projection(0,"0E"); 
h3prj_inv_m_pim_p_sim_2[j][i]->SetMarkerStyle(20);
h3prj_inv_m_pim_p_sim_2[j][i]->SetOption("e1pX0");
qqq << "m_{inv. #pi^{-}p}, GeV, sim gen ";
h3prj_inv_m_pim_p_sim_2[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h3prj_inv_m_pim_p_3_sim_2_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h3prj_inv_m_pim_p_sim_2[j][i]->SetName(qqq.str().c_str());
qqq.str("");

h3prj_inv_m_pip_pim_sim_2[j][i] = h_5dim_3_sim_2[j][i]->Projection(1,"0E"); 
h3prj_inv_m_pip_pim_sim_2[j][i]->SetMarkerStyle(20);
h3prj_inv_m_pip_pim_sim_2[j][i]->SetOption("e1pX0");
qqq << "m_{inv. #pi^{+}#pi^{-}}, GeV, sim gen ";
h3prj_inv_m_pip_pim_sim_2[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h3prj_inv_m_pip_pim_2_sim_2_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h3prj_inv_m_pip_pim_sim_2[j][i]->SetName(qqq.str().c_str());
qqq.str("");

h3prj_th_PIp_sim_2[j][i] = h_5dim_3_sim_2[j][i]->Projection(2,"0E"); 
//h3prj_th_PIp_sim_2[j][i] ->Divide(h_cos_th);
//h3prj_th_PIp_sim_2[j][i] ->SetEntries(h3prj_inv_m_pip_pim_sim_2[j][i]->GetEntries());
h3prj_th_PIp_sim_2[j][i]->SetMarkerStyle(20);
h3prj_th_PIp_sim_2[j][i]->SetOption("e1pX0");
qqq << "#theta_{#pi^{+}} in c.m., deg, sim gen ";
h3prj_th_PIp_sim_2[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h3prj_th_PIp_sim_2_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h3prj_th_PIp_sim_2[j][i]->SetName(qqq.str().c_str());
qqq.str("");

h3prj_phi_PIp_sim_2[j][i] = h_5dim_3_sim_2[j][i]->Projection(3,"0E"); 
h3prj_phi_PIp_sim_2[j][i]->SetMarkerStyle(20);
h3prj_phi_PIp_sim_2[j][i]->SetOption("e1pX0");
qqq << "#phi_{#pi^{+}} in c.m., deg, sim gen";
h3prj_phi_PIp_sim_2[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h3prj_phi_PIp_sim_2_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h3prj_phi_PIp_sim_2[j][i]->SetName(qqq.str().c_str());
qqq.str("");

h3prj_alpha_PPIm_piPIp_sim_2[j][i] = h_5dim_3_sim_2[j][i]->Projection(4,"0E");
h3prj_alpha_PPIm_piPIp_sim_2[j][i]->SetMarkerStyle(20);
h3prj_alpha_PPIm_piPIp_sim_2[j][i]->SetOption("e1pX0");
qqq << "#alpha_{p'#pi^{-}_p#pi^{+}}, deg, sim gen";
h3prj_alpha_PPIm_piPIp_sim_2[j][i]->SetTitle(qqq.str().c_str());
qqq.str("");
qqq << "h3prj_alpha_PPIm_piPIp_sim_2_"<<"q2_" << 1000*(0.425+0.05*j) << "_w_" << 1000*(1.312+0.025*i);
h3prj_alpha_PPIm_piPIp_sim_2[j][i]->SetName(qqq.str().c_str());
qqq.str("");
*/

/*for (k=1;k<=10;k++){
old_bin_cont=h1prj_th_P[j][i]->GetBinContent(k);
new_bin_cont=old_bin_cont/sin((9+18*(k-1))*M_PI/180.);
h1prj_th_P[j][i]->SetBinContent(k,new_bin_cont);

old_bin_cont=h2prj_th_PIm[j][i]->GetBinContent(k);
new_bin_cont=old_bin_cont/sin((9+18*(k-1))*M_PI/180.);
h2prj_th_PIm[j][i]->SetBinContent(k,new_bin_cont);


old_bin_cont=h3prj_th_PIp[j][i]->GetBinContent(k);
new_bin_cont=old_bin_cont/sin((9+18*(k-1))*M_PI/180.);
h3prj_th_PIp[j][i]->SetBinContent(k,new_bin_cont);


old_bin_cont=h1prj_th_P_sim[j][i]->GetBinContent(k);
new_bin_cont=old_bin_cont/sin((9+18*(k-1))*M_PI/180.);
h1prj_th_P_sim[j][i]->SetBinContent(k,new_bin_cont);


old_bin_cont=h2prj_th_PIm_sim[j][i]->GetBinContent(k);
new_bin_cont=old_bin_cont/sin((9+18*(k-1))*M_PI/180.);
h2prj_th_PIm_sim[j][i]->SetBinContent(k,new_bin_cont);


old_bin_cont=h3prj_th_PIp_sim[j][i]->GetBinContent(k);
new_bin_cont=old_bin_cont/sin((9+18*(k-1))*M_PI/180.);
h3prj_th_PIp_sim[j][i]->SetBinContent(k,new_bin_cont);



//cout << 9+18*(k-1) << " sin = " << sin((9+18*(k-1))*M_PI/180.) << "\n";



};






};
};
*/




/*

 for (j=0; j<12; j++) {
 for (i=0; i<21;i++) {
 
 //h2prj_phi_PIm[j][i]->Sumw2();
//integ= h2prj_alpha_PPIp_piPIm[j][i]->IntegralAndError(0,8,err);
integ= h2prj_phi_PIm[j][i]->IntegralAndError(0,8,err);
err1 =err;
//integ= h2prj_alpha_PPIp_piPIm[j][i]->DoIntegral(0.,360.,err);
   h_w_int[j]->Fill(1.3125+0.025*i,integ);
 //  h_w_int[j]->SetBinContent(i,integ);
 //cout << W_bin[i] <<"\n";
 //if (j == 0)cout << "err1 = " << err1 << " err_new = " << h_w_int[j]->GetBinError(i+1) << "\n";
h_w_int[j]->SetBinError(i+1,err1);
h_w_int[j]->SetMarkerStyle(20);
h_w_int[j]->SetTitle("  ");
h_w_int[j]->SetOption("e1pX0");


};
};

 
 */

 
 
 
 
 
 
 nphefile->Close(); 
  output();
     
 	outFile->Close();

	/*hist_sector1->Delete();
	hist_nphe_sector1->Delete();
	hist_sector2->Delete();
	hist_nphe_sector2->Delete();
	hist_sector3->Delete();
	hist_nphe_sector3->Delete();
	hist_sector4->Delete();
	hist_nphe_sector4->Delete();
	hist_sector5->Delete();
	hist_nphe_sector5->Delete();
	hist_sector6->Delete();
	hist_nphe_sector6->Delete();*/
	
	cout << "Qfull = " << Qfull << "\n";
	cout << "Qfull_empty = " << Qfull_empty << "\n";
	cout << "Qfull_sim = " << Qfull_sim << "\n";
	
    };

    
     

















