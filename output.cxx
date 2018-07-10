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
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring> 
#include <cstdlib>
#include "global.h"
#include "macro.h"

 using namespace std; 


Int_t ti,io,jo,tj;

ostringstream qqq3;
 macro macros; 
 TMacro *crs_plot;
int output() {
TFile *file1 = TFile::Open("out_radcorr.root","RECREATE");
file1->cd();

h_1d_rc_0425->Write("", TObject::kOverwrite);
h_1d_rc_0475->Write("", TObject::kOverwrite);
h_1d_rc_0525->Write("", TObject::kOverwrite);
h_1d_rc_0575->Write("", TObject::kOverwrite);
h_1d_rc_0625->Write("", TObject::kOverwrite);
h_1d_rc_0675->Write("", TObject::kOverwrite);
h_1d_rc_0725->Write("", TObject::kOverwrite);
h_1d_rc_0775->Write("", TObject::kOverwrite);
h_1d_rc_0825->Write("", TObject::kOverwrite);
h_1d_rc_0875->Write("", TObject::kOverwrite);
h_1d_rc_0925->Write("", TObject::kOverwrite);
h_1d_rc_0975->Write("", TObject::kOverwrite);

h_1d_rc_0425_evt->Write("", TObject::kOverwrite);
h_1d_rc_0475_evt->Write("", TObject::kOverwrite);
h_1d_rc_0525_evt->Write("", TObject::kOverwrite);
h_1d_rc_0575_evt->Write("", TObject::kOverwrite);
h_1d_rc_0625_evt->Write("", TObject::kOverwrite);
h_1d_rc_0675_evt->Write("", TObject::kOverwrite);
h_1d_rc_0725_evt->Write("", TObject::kOverwrite);
h_1d_rc_0775_evt->Write("", TObject::kOverwrite);
h_1d_rc_0825_evt->Write("", TObject::kOverwrite);
h_1d_rc_0875_evt->Write("", TObject::kOverwrite);
h_1d_rc_0925_evt->Write("", TObject::kOverwrite);
h_1d_rc_0975_evt->Write("", TObject::kOverwrite);



    outFile->cd();
     outFile->mkdir("sector1");
     outFile->mkdir("sector2");
     outFile->mkdir("sector3");
     outFile->mkdir("sector4");
     outFile->mkdir("sector5");
     outFile->mkdir("sector6");
     
     outFile->mkdir("sector1_p_fid");
     outFile->mkdir("sector2_p_fid");
     outFile->mkdir("sector3_p_fid");
     outFile->mkdir("sector4_p_fid");
     outFile->mkdir("sector5_p_fid");
     outFile->mkdir("sector6_p_fid");
     
     outFile->mkdir("sector1_pip_fid");
     outFile->mkdir("sector2_pip_fid");
     outFile->mkdir("sector3_pip_fid");
     outFile->mkdir("sector4_pip_fid");
     outFile->mkdir("sector5_pip_fid");
     outFile->mkdir("sector6_pip_fid");  
     
       outFile->mkdir("sector1_pim_fid");
     outFile->mkdir("sector2_pim_fid");
     outFile->mkdir("sector3_pim_fid");
     outFile->mkdir("sector4_pim_fid");
     outFile->mkdir("sector5_pim_fid");
     outFile->mkdir("sector6_pim_fid");   
     
     outFile->mkdir("s1_el_fid_sim");
     outFile->mkdir("s2_el_fid_sim");
     outFile->mkdir("s3_el_fid_sim");
     outFile->mkdir("s4_el_fid_sim");
     outFile->mkdir("s5_el_fid_sim");
     outFile->mkdir("s6_el_fid_sim");
     
     outFile->mkdir("s1_p_fid_sim");
     outFile->mkdir("s2_p_fid_sim");
     outFile->mkdir("s3_p_fid_sim");
     outFile->mkdir("s4_p_fid_sim");
     outFile->mkdir("s5_p_fid_sim");
     outFile->mkdir("s6_p_fid_sim");
     
     outFile->mkdir("s1_pip_fid_sim");
     outFile->mkdir("s2_pip_fid_sim");
     outFile->mkdir("s3_pip_fid_sim");
     outFile->mkdir("s4_pip_fid_sim");
     outFile->mkdir("s5_pip_fid_sim");
     outFile->mkdir("s6_pip_fid_sim");  
     
       outFile->mkdir("s1_pim_fid_sim");
     outFile->mkdir("s2_pim_fid_sim");
     outFile->mkdir("s3_pim_fid_sim");
     outFile->mkdir("s4_pim_fid_sim");
     outFile->mkdir("s5_pim_fid_sim");
     outFile->mkdir("s6_pim_fid_sim");
     
     
     	
 
     outFile->mkdir("s1_beta_vs_p");
     outFile->mkdir("s2_beta_vs_p");
     outFile->mkdir("s3_beta_vs_p");
     outFile->mkdir("s4_beta_vs_p");
     outFile->mkdir("s5_beta_vs_p");
     outFile->mkdir("s6_beta_vs_p"); 
     
     outFile->mkdir("s1_beta_vs_pip");
     outFile->mkdir("s2_beta_vs_pip");
     outFile->mkdir("s3_beta_vs_pip");
     outFile->mkdir("s4_beta_vs_pip");
     outFile->mkdir("s5_beta_vs_pip");
     outFile->mkdir("s6_beta_vs_pip"); 
     

      outFile->mkdir("s1_beta_vs_p_pim");
     outFile->mkdir("s2_beta_vs_p_pim");
     outFile->mkdir("s3_beta_vs_p_pim");
     outFile->mkdir("s4_beta_vs_p_pim");
     outFile->mkdir("s5_beta_vs_p_pim");
     outFile->mkdir("s6_beta_vs_p_pim");
     outFile->mkdir("___________");
        outFile->mkdir("s1_beta_vs_p_sim");
     outFile->mkdir("s2_beta_vs_p_sim");
     outFile->mkdir("s3_beta_vs_p_sim");
     outFile->mkdir("s4_beta_vs_p_sim");
     outFile->mkdir("s5_beta_vs_p_sim");
     outFile->mkdir("s6_beta_vs_p_sim"); 
     
     outFile->mkdir("s1_beta_vs_pip_sim");
     outFile->mkdir("s2_beta_vs_pip_sim");
     outFile->mkdir("s3_beta_vs_pip_sim");
     outFile->mkdir("s4_beta_vs_pip_sim");
     outFile->mkdir("s5_beta_vs_pip_sim");
     outFile->mkdir("s6_beta_vs_pip_sim"); 
     

      outFile->mkdir("s1_beta_vs_p_pim_sim");
     outFile->mkdir("s2_beta_vs_p_pim_sim");
     outFile->mkdir("s3_beta_vs_p_pim_sim");
     outFile->mkdir("s4_beta_vs_p_pim_sim");
     outFile->mkdir("s5_beta_vs_p_pim_sim");
     outFile->mkdir("s6_beta_vs_p_pim_sim");  
     
     

  
     
     outFile->mkdir("s1_time_p");
     outFile->mkdir("s2_time_p");
     outFile->mkdir("s3_time_p");
     outFile->mkdir("s4_time_p");
     outFile->mkdir("s5_time_p");
     outFile->mkdir("s6_time_p");
     
     outFile->mkdir("s1_time_pip");
     outFile->mkdir("s2_time_pip");
     outFile->mkdir("s3_time_pip");
     outFile->mkdir("s4_time_pip");
     outFile->mkdir("s5_time_pip");
     outFile->mkdir("s6_time_pip");
     
     outFile->mkdir("s1_time_pim");
     outFile->mkdir("s2_time_pim");
     outFile->mkdir("s3_time_pim");
     outFile->mkdir("s4_time_pim");
     outFile->mkdir("s5_time_pim");
     outFile->mkdir("s6_time_pim");
     
     outFile->mkdir("____________");
     
      outFile->mkdir("s1_time_p_sim");
     outFile->mkdir("s2_time_p_sim");
     outFile->mkdir("s3_time_p_sim");
     outFile->mkdir("s4_time_p_sim");
     outFile->mkdir("s5_time_p_sim");
     outFile->mkdir("s6_time_p_sim");
     
     outFile->mkdir("s1_time_pip_sim");
     outFile->mkdir("s2_time_pip_sim");
     outFile->mkdir("s3_time_pip_sim");
     outFile->mkdir("s4_time_pip_sim");
     outFile->mkdir("s5_time_pip_sim");
     outFile->mkdir("s6_time_pip_sim");
     
     outFile->mkdir("s1_time_pim_sim");
     outFile->mkdir("s2_time_pim_sim");
     outFile->mkdir("s3_time_pim_sim");
     outFile->mkdir("s4_time_pim_sim");
     outFile->mkdir("s5_time_pim_sim");
     outFile->mkdir("s6_time_pim_sim");  
     
     
     
     
    outFile->mkdir("th_vs_p");  
    outFile->mkdir("th_vs_p_sim");  
  
    
    
    
 //   outFile->mkdir("5dim_hist_data");    
 //   outFile->mkdir("5dim_hist_sim");
 //    outFile->mkdir("5dim_excl_hist_data"); 
  //  outFile->mkdir("5dim_excl_hist_sim");
      
       outFile->mkdir("main_top_mm_pim_0425"); 
       outFile->mkdir("main_top_mm_pim_0475"); 
       outFile->mkdir("main_top_mm_pim_0525"); 
       outFile->mkdir("main_top_mm_pim_0575"); 
       outFile->mkdir("main_top_mm_pim_0625"); 
       outFile->mkdir("main_top_mm_pim_0675"); 
       outFile->mkdir("main_top_mm_pim_0725"); 
       outFile->mkdir("main_top_mm_pim_0775"); 
       outFile->mkdir("main_top_mm_pim_0825"); 
       outFile->mkdir("main_top_mm_pim_0875"); 
       outFile->mkdir("main_top_mm_pim_0925"); 
       outFile->mkdir("main_top_mm_pim_0975"); 
       
        outFile->mkdir("main_top_mm_pim_noq2dep");
       outFile->mkdir("excl_top_plots");
       
        outFile->mkdir("main_top_plots");
 /*      outFile->mkdir("miss_mass_pim_th_dep_pim"); 
       outFile->mkdir("miss_mass_pim_th_dep_pip"); 
       outFile->mkdir("miss_mass_pim_th_dep_pr"); 
       outFile->mkdir("miss_mass_pim_top_th_dep_pr"); 
       
       outFile->mkdir("miss_mass_pim_top_th_dep_arr"); */
       
         outFile->mkdir("inv_m_fsi"); 
       
      
  TDirectory *dir, *dir1, *dir2;  
    
for(jo=0; jo<12; jo++){
qqq3.str("");    
  qqq3 <<"q2_" << 0.425+0.05*jo;    
outFile->mkdir(qqq3.str().c_str());  
outFile->cd(qqq3.str().c_str());
dir = outFile->GetDirectory(qqq3.str().c_str());
dir->cd();
for (io=0; io<21; io++) {
qqq3.str("");    
  qqq3 <<"w_" << 1.3125+0.025*io ;    
dir->mkdir(qqq3.str().c_str());

dir->cd(qqq3.str().c_str());
qqq3.str(""); 
qqq3 <<"q2_" << 0.425+0.05*jo<<"/"<<"w_" << 1.3125+0.025*io ; 
dir1 = outFile->GetDirectory(qqq3.str().c_str());
dir1->cd();

h_5dim_excl_1[jo][io]->Write("", TObject::kOverwrite);
h_5dim_excl_2[jo][io]->Write("", TObject::kOverwrite);
h_5dim_excl_3[jo][io]->Write("", TObject::kOverwrite);


h_5dim_excl_1_sim[jo][io]->Write("", TObject::kOverwrite);
h_5dim_excl_2_sim[jo][io]->Write("", TObject::kOverwrite);
h_5dim_excl_3_sim[jo][io]->Write("", TObject::kOverwrite);

h_5dim_excl_1_sim_evt[jo][io]->Write("", TObject::kOverwrite);
h_5dim_excl_2_sim_evt[jo][io]->Write("", TObject::kOverwrite);
h_5dim_excl_3_sim_evt[jo][io]->Write("", TObject::kOverwrite);

h_5dim_pim_1[jo][io]->Write("", TObject::kOverwrite);
h_5dim_pim_2[jo][io]->Write("", TObject::kOverwrite);
h_5dim_pim_3[jo][io]->Write("", TObject::kOverwrite);

h_5dim_pim_1_sim[jo][io]->Write("", TObject::kOverwrite);
h_5dim_pim_2_sim[jo][io]->Write("", TObject::kOverwrite);
h_5dim_pim_3_sim[jo][io]->Write("", TObject::kOverwrite);

h_5dim_pim_1_sim_evt[jo][io]->Write("", TObject::kOverwrite);
h_5dim_pim_2_sim_evt[jo][io]->Write("", TObject::kOverwrite);
h_5dim_pim_3_sim_evt[jo][io]->Write("", TObject::kOverwrite);

h_5dim_1_sim_gen[jo][io]->Write("", TObject::kOverwrite);
h_5dim_2_sim_gen[jo][io]->Write("", TObject::kOverwrite);
h_5dim_3_sim_gen[jo][io]->Write("", TObject::kOverwrite);

h_5dim_1_sim_gen_evt[jo][io]->Write("", TObject::kOverwrite);
h_5dim_2_sim_gen_evt[jo][io]->Write("", TObject::kOverwrite);
h_5dim_3_sim_gen_evt[jo][io]->Write("", TObject::kOverwrite);


 
qqq3.str("");  

   };   
    }; 
   
    
    
/* for(jo=0; jo<12; jo++){
qqq3.str("");    
  qqq3 <<"sim_q2_" << 1000*(0.425+0.05*jo);    
outFile->mkdir(qqq3.str().c_str());  
outFile->cd(qqq3.str().c_str());
dir = outFile->GetDirectory(qqq3.str().c_str());
dir->cd();
for (io=0; io<21; io++) {
qqq3.str("");    
  qqq3 <<"sim_w_" << 1000*(1.312+0.025*io) ;    
dir->mkdir(qqq3.str().c_str());

dir->cd(qqq3.str().c_str());
qqq3.str(""); 
qqq3 <<"sim_q2_" << 1000*(0.425+0.05*jo)<<"/"<<"sim_w_" << 1000*(1.312+0.025*io) ; 
dir1 = outFile->GetDirectory(qqq3.str().c_str());
dir1->cd();

h1prj_inv_m_pip_p_sim[jo][io]->SetDirectory(dir1); 
h1prj_inv_m_pip_pim_sim[jo][io]->SetDirectory(dir1); 
h1prj_th_P_sim[jo][io]->SetDirectory(dir1);
h1prj_phi_P_sim[jo][io]->SetDirectory(dir1); 
h1prj_alpha_PIpPIm_pipf_sim[jo][io]->SetDirectory(dir1);

h2prj_inv_m_pip_p_sim[jo][io]->SetDirectory(dir1); 
h2prj_inv_m_pip_pim_sim[jo][io]->SetDirectory(dir1);
h2prj_th_PIm_sim[jo][io]->SetDirectory(dir1);
h2prj_phi_PIm_sim[jo][io] ->SetDirectory(dir1);
h2prj_alpha_PPIp_piPIm_sim[jo][io] ->SetDirectory(dir1);

h3prj_inv_m_pim_p_sim[jo][io]->SetDirectory(dir1); 
h3prj_inv_m_pip_pim_sim[jo][io]->SetDirectory(dir1);
h3prj_th_PIp_sim[jo][io]  ->SetDirectory(dir1);
h3prj_phi_PIp_sim[jo][io] ->SetDirectory(dir1); 
h3prj_alpha_PPIm_piPIp_sim[jo][io] ->SetDirectory(dir1);
     
qqq3.str("");  

 crs_plot = macros.dif_crs_plot(jo,io);
 //crs_plot->Exec();
 crs_plot->Write("", TObject::kOverwrite);
  macros.dif_crs_plot_write(jo,io);
   };   
    }; */
       
    
    
    
    
  
//photoelectrons

 for(jo=0; jo<6; jo++){
qqq3.str("");    
  qqq3 <<"ph_el_sector_" << jo+1;    
outFile->mkdir(qqq3.str().c_str());  
outFile->cd(qqq3.str().c_str());
dir = outFile->GetDirectory(qqq3.str().c_str());
dir->cd();
dir->mkdir("left");

qqq3.str("");    
  qqq3 <<"ph_el_sector_" << jo+1 << "/left";
 dir1 = outFile->GetDirectory(qqq3.str().c_str());
dir1->cd(); 
for (io=0; io<20; io++){
ph_el_left[jo][io]->SetDirectory(dir1); 
};
dir->mkdir("both");
qqq3.str("");    
  qqq3 <<"ph_el_sector_" << jo+1 << "/both";
 dir1 = outFile->GetDirectory(qqq3.str().c_str());
dir1->cd(); 
for (io=0; io<20; io++){
ph_el_both[jo][io]->SetDirectory(dir1);
}
dir->mkdir("right");
qqq3.str("");    
  qqq3 <<"ph_el_sector_" << jo+1 << "/right";
 dir1 = outFile->GetDirectory(qqq3.str().c_str());
dir1->cd(); 
for (io=0; io<20; io++){
ph_el_right[jo][io]->SetDirectory(dir1);
}


outFile->cd(); 
   
};  

h_sim_mom_corr_test->Write("", TObject::kOverwrite); 
h_z_P->Write("", TObject::kOverwrite); 
h_z_PIp->Write("", TObject::kOverwrite); 
h_z_PIm ->Write("", TObject::kOverwrite);

h_z_P_sim ->Write("", TObject::kOverwrite);
h_z_PIp_sim->Write("", TObject::kOverwrite);
h_z_PIm_sim->Write("", TObject::kOverwrite);

h_dc_y_vs_x_el->Write("", TObject::kOverwrite);

h_z_corr1_data->Write("", TObject::kOverwrite);
h_z_corr2_data->Write("", TObject::kOverwrite);
h_z_corr1_sim->Write("", TObject::kOverwrite);
h_z_corr2_sim->Write("", TObject::kOverwrite);





      h_cos_th->Write("", TObject::kOverwrite);
      
      hist_ectot_sector1->Write("", TObject::kOverwrite); 
      hist_ectot_sector2->Write("", TObject::kOverwrite);
      hist_ectot_sector3->Write("", TObject::kOverwrite);
      hist_ectot_sector4->Write("", TObject::kOverwrite);
      hist_ectot_sector5->Write("", TObject::kOverwrite);
      hist_ectot_sector6->Write("", TObject::kOverwrite);
      hist_ectot_sector1_sim->Write("", TObject::kOverwrite); 
      hist_ectot_sector2_sim->Write("", TObject::kOverwrite);
      hist_ectot_sector3_sim->Write("", TObject::kOverwrite);
      hist_ectot_sector4_sim->Write("", TObject::kOverwrite);
      hist_ectot_sector5_sim->Write("", TObject::kOverwrite);
      hist_ectot_sector6_sim->Write("", TObject::kOverwrite);
      h_mm_0_vs_npart->Write("", TObject::kOverwrite);
      h_mm_pim_vs_npart->Write("", TObject::kOverwrite);
      h_mm_pip_vs_npart->Write("", TObject::kOverwrite);
      h_mm_0_vs_npart_sim->Write("", TObject::kOverwrite);
      h_mm_pim_vs_npart_sim->Write("", TObject::kOverwrite);
      h_mm_pip_vs_npart_sim->Write("", TObject::kOverwrite);
      
    	
      hist_ltime->Write("", TObject::kOverwrite); 
      hist_n_incl->Write("", TObject::kOverwrite);
      hist_n_elast->Write("", TObject::kOverwrite); 
      hist_ltime_1d->Write("", TObject::kOverwrite); 
      hist_n_incl_1d->Write("", TObject::kOverwrite);
      hist_n_elast_1d->Write("", TObject::kOverwrite); 
            
     // hist_sector1->Write("", TObject::kOverwrite);
    //  hist_nphe_sector1->Write("", TObject::kOverwrite);
      

      
      W_2pi_selection->Write("", TObject::kOverwrite);
      W_2pi_selection_sim->Write("", TObject::kOverwrite);
      
      h_w_1dim->Write("", TObject::kOverwrite);
      h_w_1dim_sim->Write("", TObject::kOverwrite);
      
      
      hist_w_hadr_all_reg->Write("", TObject::kOverwrite);
      hist_w_el_all_reg->Write("", TObject::kOverwrite);
    
                    
       h_miss_mass_0->Write("", TObject::kOverwrite);
       h_miss_mass_0_sim->Write("", TObject::kOverwrite); 
       hist_miss_en_0->Write("", TObject::kOverwrite);
       h_miss_en_0_sim->Write("", TObject::kOverwrite);
      
       
       h_miss_mom_0_nocut->Write("", TObject::kOverwrite);
       h_miss_mom_0_cut_on0->Write("", TObject::kOverwrite);
       h_miss_mom_0_cut_onpim->Write("", TObject::kOverwrite);
       h_miss_mom_0_nocut_sim->Write("", TObject::kOverwrite);
       h_miss_mom_0_cut_on0_sim->Write("", TObject::kOverwrite);
       h_miss_mom_0_cut_onpim_sim ->Write("", TObject::kOverwrite);  
       
       
       
       h_miss_mom_0->Write("", TObject::kOverwrite);
       h_miss_mom_0_sim->Write("", TObject::kOverwrite);
     
       
      
       hist_PIm_miss->Write("", TObject::kOverwrite);
       h_PIm_miss_sim->Write("", TObject::kOverwrite);
       hist_PIm_miss_en->Write("", TObject::kOverwrite);
       h_PIm_miss_en_sim->Write("", TObject::kOverwrite);
       
       hist_PIp_miss->Write("", TObject::kOverwrite);
       h_PIp_miss_sim->Write("", TObject::kOverwrite); 
       hist_PIp_miss_en->Write("", TObject::kOverwrite);
       h_PIp_miss_en_sim->Write("", TObject::kOverwrite);
	
	
	hist_PIp_miss_d_bef->Write("", TObject::kOverwrite);
	h_PIp_miss_d_bef_sim->Write("", TObject::kOverwrite);
	hist_PIp_miss_d->Write("", TObject::kOverwrite);
	h_PIp_miss_d_sim->Write("", TObject::kOverwrite);
	


    
   hist_PIm_miss_all_reg_1->Write("", TObject::kOverwrite);
   hist_PIm_miss_all_reg_2->Write("", TObject::kOverwrite);
   hist_PIm_miss_all_reg_1_sim->Write("", TObject::kOverwrite);
   hist_PIm_miss_all_reg_2_sim->Write("", TObject::kOverwrite);
      
      
      
       hist_P_miss->Write("", TObject::kOverwrite);
       hist_P_miss_en->Write("", TObject::kOverwrite);


h_inv_m_pip_pim ->Write("", TObject::kOverwrite);
h_inv_m_pip_p->Write("", TObject::kOverwrite); 
h_inv_m_pim_p->Write("", TObject::kOverwrite); 
 
h_inv_m_pip_pim_sim ->Write("", TObject::kOverwrite);
h_inv_m_pip_p_sim->Write("", TObject::kOverwrite); 
h_inv_m_pim_p_sim->Write("", TObject::kOverwrite); 


hist_z_el_1->Write("", TObject::kOverwrite); 
hist_z_el_2->Write("", TObject::kOverwrite); 
hist_z_el_3->Write("", TObject::kOverwrite); 
hist_z_el_4->Write("", TObject::kOverwrite); 
hist_z_el_5->Write("", TObject::kOverwrite); 
hist_z_el_6->Write("", TObject::kOverwrite); 

hist_z_el_1_empty->Write("", TObject::kOverwrite); 
hist_z_el_2_empty->Write("", TObject::kOverwrite); 
hist_z_el_3_empty->Write("", TObject::kOverwrite); 
hist_z_el_4_empty->Write("", TObject::kOverwrite); 
hist_z_el_5_empty->Write("", TObject::kOverwrite); 
hist_z_el_6_empty->Write("", TObject::kOverwrite); 


hist_z_el_1_sim_1->Write("", TObject::kOverwrite); 
hist_z_el_2_sim_1->Write("", TObject::kOverwrite); 
hist_z_el_3_sim_1->Write("", TObject::kOverwrite); 
hist_z_el_4_sim_1->Write("", TObject::kOverwrite); 
hist_z_el_5_sim_1->Write("", TObject::kOverwrite); 
hist_z_el_6_sim_1->Write("", TObject::kOverwrite); 


hist_z_el_1_sim_2->Write("", TObject::kOverwrite); 
hist_z_el_2_sim_2->Write("", TObject::kOverwrite); 
hist_z_el_3_sim_2->Write("", TObject::kOverwrite); 
hist_z_el_4_sim_2->Write("", TObject::kOverwrite); 
hist_z_el_5_sim_2->Write("", TObject::kOverwrite); 
hist_z_el_6_sim_2->Write("", TObject::kOverwrite); 







h_inv_m_pip_pim_bin ->  Write("", TObject::kOverwrite);
h_inv_m_pip_p_bin ->  Write("", TObject::kOverwrite);
h_inv_m_pim_p_bin ->  Write("", TObject::kOverwrite);
h_theta_PIm_cm -> Write("", TObject::kOverwrite);
h_theta_PIp_cm -> Write("", TObject::kOverwrite);
h_theta_P_cm -> Write("", TObject::kOverwrite);
h_phi_PIm_cm -> Write("", TObject::kOverwrite);
h_phi_PIp_cm -> Write("", TObject::kOverwrite);
h_phi_P_cm -> Write("", TObject::kOverwrite);
h_alpha_PIpPIm_pipf -> Write("", TObject::kOverwrite);
h_alpha_PPIp_piPIm -> Write("", TObject::kOverwrite);
h_alpha_PPIm_piPIp -> Write("", TObject::kOverwrite);


h_cc_nphe_total_s1 -> Write("", TObject::kOverwrite);
h_cc_nphe_total_s2 -> Write("", TObject::kOverwrite);
h_cc_nphe_total_s3 -> Write("", TObject::kOverwrite);
h_cc_nphe_total_s4 -> Write("", TObject::kOverwrite);
h_cc_nphe_total_s5 -> Write("", TObject::kOverwrite);
h_cc_nphe_total_s6 -> Write("", TObject::kOverwrite);


h_cc_nphe_final_s1 -> Write("", TObject::kOverwrite);
h_cc_nphe_final_s2 -> Write("", TObject::kOverwrite);
h_cc_nphe_final_s3 -> Write("", TObject::kOverwrite);
h_cc_nphe_final_s4 -> Write("", TObject::kOverwrite);
h_cc_nphe_final_s5 -> Write("", TObject::kOverwrite);
h_cc_nphe_final_s6 -> Write("", TObject::kOverwrite);


th_cc_vs_seg_1  -> Write("", TObject::kOverwrite);
th_cc_vs_seg_2  -> Write("", TObject::kOverwrite);
th_cc_vs_seg_3  -> Write("", TObject::kOverwrite);
th_cc_vs_seg_4  -> Write("", TObject::kOverwrite);
th_cc_vs_seg_5  -> Write("", TObject::kOverwrite);
th_cc_vs_seg_6  -> Write("", TObject::kOverwrite);

th_cc_vs_seg_1_sim -> Write("", TObject::kOverwrite);
th_cc_vs_seg_2_sim -> Write("", TObject::kOverwrite);
th_cc_vs_seg_3_sim -> Write("", TObject::kOverwrite);
th_cc_vs_seg_4_sim -> Write("", TObject::kOverwrite);
th_cc_vs_seg_5_sim -> Write("", TObject::kOverwrite);
th_cc_vs_seg_6_sim -> Write("", TObject::kOverwrite);







/*norm_nphe_s1 -> Write("", TObject::kOverwrite);
norm_nphe_s2 -> Write("", TObject::kOverwrite);
norm_nphe_s3 -> Write("", TObject::kOverwrite);
norm_nphe_s4 -> Write("", TObject::kOverwrite);
norm_nphe_s5 -> Write("", TObject::kOverwrite);
norm_nphe_s6 -> Write("", TObject::kOverwrite);*/


    //  hist_sector1->Write("", TObject::kOverwrite);
    //  hist_nphe_sector1->Write("", TObject::kOverwrite);
	
     // avrg_nphe_sector1->Write("", TObject::kOverwrite);
      
     // hist_sector2->Write("", TObject::kOverwrite);
     // hist_nphe_sector2->Write("", TObject::kOverwrite);
      
     
     // avrg_nphe_sector2->Write("", TObject::kOverwrite);
            
     // hist_sector3->Write("", TObject::kOverwrite);
     // hist_nphe_sector3->Write("", TObject::kOverwrite); 
      
 

     // avrg_nphe_sector3->Write("", TObject::kOverwrite);     
      
      
    //  hist_sector4->Write("", TObject::kOverwrite);
     // hist_nphe_sector4->Write("", TObject::kOverwrite); 
      
      

     // avrg_nphe_sector4->Write("", TObject::kOverwrite);     
         
        
     // hist_sector5->Write("", TObject::kOverwrite);
     // hist_nphe_sector5->Write("", TObject::kOverwrite); 
      

      //avrg_nphe_sector5->Write("", TObject::kOverwrite);     
      
      
      
      
      //hist_sector6->Write("", TObject::kOverwrite);  
     // hist_nphe_sector6->Write("", TObject::kOverwrite);    
      

     // avrg_nphe_sector6->Write("", TObject::kOverwrite);     
      
        
      nphe_sector1->Write("", TObject::kOverwrite); 
      nphe_sector1_after->Write("", TObject::kOverwrite);  
      nphe_sector2->Write("", TObject::kOverwrite); 
      nphe_sector2_after->Write("", TObject::kOverwrite);        
      nphe_sector3->Write("", TObject::kOverwrite); 
      nphe_sector3_after->Write("", TObject::kOverwrite);  
      nphe_sector4->Write("", TObject::kOverwrite); 
      nphe_sector4_after->Write("", TObject::kOverwrite);   
      nphe_sector5->Write("", TObject::kOverwrite); 
      nphe_sector5_after->Write("", TObject::kOverwrite);  
      nphe_sector6->Write("", TObject::kOverwrite); 
      nphe_sector6_after->Write("", TObject::kOverwrite); 
      
     
      //z_el_1->Write("", TObject::kOverwrite);
      //z_el_2->Write("", TObject::kOverwrite);
      //z_el_3->Write("", TObject::kOverwrite);
      //z_el_4->Write("", TObject::kOverwrite);
      //z_el_5->Write("", TObject::kOverwrite);
      //z_el_6->Write("", TObject::kOverwrite);
     
//th_vs_p_e_1->Write("", TObject::kOverwrite);  
//th_vs_p_e_2->Write("", TObject::kOverwrite); 
//th_vs_p_p_1->Write("", TObject::kOverwrite);  
//th_vs_p_p_2->Write("", TObject::kOverwrite);  
//th_vs_p_pip_1->Write("", TObject::kOverwrite);  
//th_vs_p_pip_2->Write("", TObject::kOverwrite); 
  /*    
    outFile->cd("5dim_hist_data");
   for (jo=0; jo<12; jo++){  
    for (io=0; io<21; io++){   
     
      h_5dim_1[jo][io]->Write("", TObject::kOverwrite); 
        h_5dim_2[jo][io]->Write("", TObject::kOverwrite); 
        h_5dim_3[jo][io]->Write("", TObject::kOverwrite); 
       h_5dim_1_empty[jo][io]->Write("", TObject::kOverwrite); 
        h_5dim_2_empty[jo][io]->Write("", TObject::kOverwrite); 
        h_5dim_3_empty[jo][io]->Write("", TObject::kOverwrite);
           
      };
      };
    
     outFile->cd("5dim_hist_sim");
   for (jo=0; jo<12; jo++){  
    for (io=0; io<21; io++){   
     
      h_5dim_1_sim_1[jo][io]->Write("", TObject::kOverwrite); 
        h_5dim_2_sim_1[jo][io]->Write("", TObject::kOverwrite); 
        h_5dim_3_sim_1[jo][io]->Write("", TObject::kOverwrite); 
       h_5dim_1_sim_2[jo][io]->Write("", TObject::kOverwrite); 
        h_5dim_2_sim_2[jo][io]->Write("", TObject::kOverwrite); 
        h_5dim_3_sim_2[jo][io]->Write("", TObject::kOverwrite);
           
      };
      };   
      
  outFile->cd("5dim_excl_hist_data");
    for (io=0; io<9; io++){   
     
      h_5dim_excl_1[io]->Write("", TObject::kOverwrite); 
        h_5dim_excl_2[io]->Write("", TObject::kOverwrite); 
        h_5dim_excl_3[io]->Write("", TObject::kOverwrite); 
       h_5dim_excl_1_empty[io]->Write("", TObject::kOverwrite); 
        h_5dim_excl_2_empty[io]->Write("", TObject::kOverwrite); 
        h_5dim_excl_3_empty[io]->Write("", TObject::kOverwrite);
           
     };
    
     outFile->cd("5dim_excl_hist_sim");
    for (io=0; io<9; io++){   
     
      h_5dim_1_sim_excl_1[io]->Write("", TObject::kOverwrite); 
        h_5dim_2_sim_excl_1[io]->Write("", TObject::kOverwrite); 
        h_5dim_3_sim_excl_1[io]->Write("", TObject::kOverwrite); 
       h_5dim_1_sim_excl_2[io]->Write("", TObject::kOverwrite); 
        h_5dim_2_sim_excl_2[io]->Write("", TObject::kOverwrite); 
        h_5dim_3_sim_excl_2[io]->Write("", TObject::kOverwrite);
           
       };    */    
      
              
 outFile->cd("sector1");
 for (ti=0; ti<17; ti++){
  ph_vs_th_1pe[ti]->Write("", TObject::kOverwrite); 
  ph_vs_th_1pe_fid[ti]->Write("", TObject::kOverwrite);
  };
 ph_vs_th_1->Write("", TObject::kOverwrite);
 
  outFile->cd("sector2");
  for (ti=0; ti<17; ti++){
  ph_vs_th_2pe[ti]->Write("", TObject::kOverwrite); 
  ph_vs_th_2pe_fid[ti]->Write("", TObject::kOverwrite);
  };
 ph_vs_th_2->Write("", TObject::kOverwrite);
 
  outFile->cd("sector3");
  for (ti=0; ti<17; ti++){
  ph_vs_th_3pe[ti]->Write("", TObject::kOverwrite); 
  ph_vs_th_3pe_fid[ti]->Write("", TObject::kOverwrite);
  };
   ph_vs_th_3->Write("", TObject::kOverwrite);
 
  outFile->cd("sector4");
  for (ti=0; ti<17; ti++){
  ph_vs_th_4pe[ti]->Write("", TObject::kOverwrite); 
  ph_vs_th_4pe_fid[ti]->Write("", TObject::kOverwrite);
  };
 ph_vs_th_4->Write("", TObject::kOverwrite);
 
  outFile->cd("sector5");
  for (ti=0; ti<17; ti++){
  ph_vs_th_5pe[ti]->Write("", TObject::kOverwrite); 
  ph_vs_th_5pe_fid[ti]->Write("", TObject::kOverwrite);
  };
 ph_vs_th_5->Write("", TObject::kOverwrite);
 
  outFile->cd("sector6");
  for (ti=0; ti<17; ti++){
  ph_vs_th_6pe[ti]->Write("", TObject::kOverwrite); 
  ph_vs_th_6pe_fid[ti]->Write("", TObject::kOverwrite);
  };
 ph_vs_th_6->Write("", TObject::kOverwrite);
 
 ////kkkk
 outFile->cd("sector1_p_fid");
 for (ti=0; ti<6; ti++){
 ph_th_p_1[ti]->Write("", TObject::kOverwrite); 
 };
 ph_vs_th_p_1_w->Write("", TObject::kOverwrite);
 ph_vs_th_p_1->Write("", TObject::kOverwrite);
 
 outFile->cd("sector2_p_fid");
 for (ti=0; ti<6; ti++){
 ph_th_p_2[ti]->Write("", TObject::kOverwrite); 
 };
 ph_vs_th_p_2_w->Write("", TObject::kOverwrite);
 ph_vs_th_p_2->Write("", TObject::kOverwrite);
 
 outFile->cd("sector3_p_fid");
 for (ti=0; ti<6; ti++){
 ph_th_p_3[ti]->Write("", TObject::kOverwrite); 
 };
 ph_vs_th_p_3_w->Write("", TObject::kOverwrite);
 ph_vs_th_p_3->Write("", TObject::kOverwrite);
 
 outFile->cd("sector4_p_fid");
 for (ti=0; ti<6; ti++){
 ph_th_p_4[ti]->Write("", TObject::kOverwrite); 
 };
 ph_vs_th_p_4_w->Write("", TObject::kOverwrite);
 ph_vs_th_p_4->Write("", TObject::kOverwrite);
 
 outFile->cd("sector5_p_fid");
 for (ti=0; ti<6; ti++){
 ph_th_p_5[ti]->Write("", TObject::kOverwrite); 
 };
 ph_vs_th_p_5_w->Write("", TObject::kOverwrite);
 ph_vs_th_p_5->Write("", TObject::kOverwrite);
 
 outFile->cd("sector6_p_fid");
 for (ti=0; ti<6; ti++){
 ph_th_p_6[ti]->Write("", TObject::kOverwrite); 
 };
 ph_vs_th_p_6_w->Write("", TObject::kOverwrite);
 ph_vs_th_p_6->Write("", TObject::kOverwrite);
 
 /////uuuu
 
 outFile->cd("sector1_pip_fid");
 for (ti=0; ti<6; ti++){
 ph_th_pip_1[ti]->Write("", TObject::kOverwrite); 
 };
 ph_vs_th_pip_1->Write("", TObject::kOverwrite);
 outFile->cd("sector2_pip_fid");
 for (ti=0; ti<6; ti++){
 ph_th_pip_2[ti]->Write("", TObject::kOverwrite); 
 };
 ph_vs_th_pip_2->Write("", TObject::kOverwrite);
 outFile->cd("sector3_pip_fid");
 for (ti=0; ti<6; ti++){
 ph_th_pip_3[ti]->Write("", TObject::kOverwrite); 
 };
 ph_vs_th_pip_3->Write("", TObject::kOverwrite);
 outFile->cd("sector4_pip_fid"); 
 for (ti=0; ti<6; ti++){
 ph_th_pip_4[ti]->Write("", TObject::kOverwrite); 
 };
 ph_vs_th_pip_4->Write("", TObject::kOverwrite);
  outFile->cd("sector5_pip_fid");
  for (ti=0; ti<6; ti++){
 ph_th_pip_5[ti]->Write("", TObject::kOverwrite); 
 };
  ph_vs_th_pip_5->Write("", TObject::kOverwrite);
 outFile->cd("sector6_pip_fid");
 for (ti=0; ti<6; ti++){
 ph_th_pip_6[ti]->Write("", TObject::kOverwrite); 
 };
 ph_vs_th_pip_6->Write("", TObject::kOverwrite);
 
 ////pppp
 outFile->cd("sector1_pim_fid");
 for (ti=0; ti<15; ti++){
 ph_vs_th_pim[0][ti]->Write("", TObject::kOverwrite); 
 };
 ph_th_pim_all_p[0]->Write("", TObject::kOverwrite);
 
 outFile->cd("sector2_pim_fid");
  for (ti=0; ti<10; ti++){
ph_vs_th_pim[1][ti]->Write("", TObject::kOverwrite); 
 };
 ph_th_pim_all_p[1]->Write("", TObject::kOverwrite);
 
 outFile->cd("sector3_pim_fid");
  for (ti=0; ti<10; ti++){
 ph_vs_th_pim[2][ti]->Write("", TObject::kOverwrite); 
 };
 ph_th_pim_all_p[2]->Write("", TObject::kOverwrite);
 
 outFile->cd("sector4_pim_fid"); 
  for (ti=0; ti<10; ti++){
 ph_vs_th_pim[3][ti]->Write("", TObject::kOverwrite); 
 };
 ph_th_pim_all_p[3]->Write("", TObject::kOverwrite);
 
  outFile->cd("sector5_pim_fid");
   for (ti=0; ti<10; ti++){
ph_vs_th_pim[4][ti]->Write("", TObject::kOverwrite); 
 };
 ph_th_pim_all_p[4]->Write("", TObject::kOverwrite);
 
 outFile->cd("sector6_pim_fid");
  for (ti=0; ti<10; ti++){
 ph_vs_th_pim[5][ti]->Write("", TObject::kOverwrite); 
 };
 ph_th_pim_all_p[5]->Write("", TObject::kOverwrite);
 
 ///////////////
outFile->cd("s1_el_fid_sim");
 for (ti=0; ti<7; ti++){
 ph_vs_th_el_sim[0][ti]->Write("", TObject::kOverwrite); 
 };
 
 outFile->cd("s2_el_fid_sim");
  for (ti=0; ti<7; ti++){
ph_vs_th_el_sim[1][ti]->Write("", TObject::kOverwrite); 
 };

 outFile->cd("s3_el_fid_sim");
  for (ti=0; ti<7; ti++){
 ph_vs_th_el_sim[2][ti]->Write("", TObject::kOverwrite); 
 };
  
 outFile->cd("s4_el_fid_sim"); 
  for (ti=0; ti<7; ti++){
 ph_vs_th_el_sim[3][ti]->Write("", TObject::kOverwrite); 
 };
  
  outFile->cd("s5_el_fid_sim");
   for (ti=0; ti<7; ti++){
ph_vs_th_el_sim[4][ti]->Write("", TObject::kOverwrite); 
 };
  
 outFile->cd("s6_el_fid_sim");
  for (ti=0; ti<7; ti++){
 ph_vs_th_el_sim[5][ti]->Write("", TObject::kOverwrite); 
 };
  
 //////////
 
outFile->cd("s1_p_fid_sim");
 //for (ti=0; ti<5; ti++){
 ph_vs_th_p_sim[0]->Write("", TObject::kOverwrite); 
// };
 
 outFile->cd("s2_p_fid_sim");
//  for (ti=0; ti<5; ti++){
ph_vs_th_p_sim[1]->Write("", TObject::kOverwrite); 
// };

 outFile->cd("s3_p_fid_sim");
//  for (ti=0; ti<5; ti++){
 ph_vs_th_p_sim[2]->Write("", TObject::kOverwrite); 
// };
  
 outFile->cd("s4_p_fid_sim"); 
//  for (ti=0; ti<5; ti++){
 ph_vs_th_p_sim[3]->Write("", TObject::kOverwrite); 
// };
  
  outFile->cd("s5_p_fid_sim");
//   for (ti=0; ti<5; ti++){
ph_vs_th_p_sim[4]->Write("", TObject::kOverwrite); 
// };
  
 outFile->cd("s6_p_fid_sim");
//  for (ti=0; ti<5; ti++){
 ph_vs_th_p_sim[5]->Write("", TObject::kOverwrite); 
// };  
 ///
 outFile->cd("s1_pip_fid_sim");
// for (ti=0; ti<6; ti++){
 ph_vs_th_pip_sim[0]->Write("", TObject::kOverwrite); 
// };
 
 outFile->cd("s2_pip_fid_sim");
//  for (ti=0; ti<6; ti++){
ph_vs_th_pip_sim[1]->Write("", TObject::kOverwrite); 
// };

 outFile->cd("s3_pip_fid_sim");
//  for (ti=0; ti<6; ti++){
 ph_vs_th_pip_sim[2]->Write("", TObject::kOverwrite); 
// };
  
 outFile->cd("s4_pip_fid_sim"); 
//  for (ti=0; ti<6; ti++){
 ph_vs_th_pip_sim[3]->Write("", TObject::kOverwrite); 
// };
  
  outFile->cd("s5_pip_fid_sim");
//   for (ti=0; ti<6; ti++){
ph_vs_th_pip_sim[4]->Write("", TObject::kOverwrite); 
// };
  
 outFile->cd("s6_pip_fid_sim");
//  for (ti=0; ti<6; ti++){
ph_vs_th_pip_sim[5]->Write("", TObject::kOverwrite); 
// };
//
 outFile->cd("s1_pim_fid_sim");
 for (ti=0; ti<5; ti++){
 ph_vs_th_pim_sim[0][ti]->Write("", TObject::kOverwrite); 
 };
 
 outFile->cd("s2_pim_fid_sim");
  for (ti=0; ti<5; ti++){
ph_vs_th_pim_sim[1][ti]->Write("", TObject::kOverwrite); 
 };

 outFile->cd("s3_pim_fid_sim");
  for (ti=0; ti<5; ti++){
 ph_vs_th_pim_sim[2][ti]->Write("", TObject::kOverwrite); 
 };
  
 outFile->cd("s4_pim_fid_sim"); 
  for (ti=0; ti<5; ti++){
 ph_vs_th_pim_sim[3][ti]->Write("", TObject::kOverwrite); 
 };
  
  outFile->cd("s5_pim_fid_sim");
   for (ti=0; ti<5; ti++){
ph_vs_th_pim_sim[4][ti]->Write("", TObject::kOverwrite); 
 };
  
 outFile->cd("s6_pim_fid_sim");
  for (ti=0; ti<5; ti++){
 ph_vs_th_pim_sim[5][ti]->Write("", TObject::kOverwrite); 
 };
 
 
 
 
 
 /////iiii
  outFile->cd("s1_beta_vs_p");
  for (ti=0; ti<48; ti++){
 beta_vs_p_p[0][ti]->Write("", TObject::kOverwrite); 
 };
  outFile->cd("s2_beta_vs_p");
 for (ti=0; ti<48; ti++){
 beta_vs_p_p[1][ti]->Write("", TObject::kOverwrite);
 }; 
  outFile->cd("s3_beta_vs_p");
 for (ti=0; ti<48; ti++){
 beta_vs_p_p[2][ti]->Write("", TObject::kOverwrite);
 }; 
  outFile->cd("s4_beta_vs_p"); 
 for (ti=0; ti<48; ti++){
 beta_vs_p_p[3][ti]->Write("", TObject::kOverwrite);
 }; 
  outFile->cd("s5_beta_vs_p");
  for (ti=0; ti<48; ti++){
 beta_vs_p_p[4][ti]->Write("", TObject::kOverwrite); 
 };
  outFile->cd("s6_beta_vs_p");
  for (ti=0; ti<48; ti++){
 beta_vs_p_p[5][ti]->Write("", TObject::kOverwrite); 
 };
 
 
/////iiii
  outFile->cd("s1_beta_vs_pip");
  for (ti=0; ti<48; ti++){
 beta_vs_p_pip[0][ti]->Write("", TObject::kOverwrite); 
 };
  outFile->cd("s2_beta_vs_pip");
 for (ti=0; ti<48; ti++){
 beta_vs_p_pip[1][ti]->Write("", TObject::kOverwrite); 
 }; 
  outFile->cd("s3_beta_vs_pip");
 for (ti=0; ti<48; ti++){
 beta_vs_p_pip[2][ti]->Write("", TObject::kOverwrite); 
 }; 
  outFile->cd("s4_beta_vs_pip"); 
 for (ti=0; ti<48; ti++){
 beta_vs_p_pip[3][ti]->Write("", TObject::kOverwrite); 
 }; 
  outFile->cd("s5_beta_vs_pip");
  for (ti=0; ti<48; ti++){
 beta_vs_p_pip[4][ti]->Write("", TObject::kOverwrite); 
 };
  outFile->cd("s6_beta_vs_pip");
  for (ti=0; ti<48; ti++){
 beta_vs_p_pip[5][ti]->Write("", TObject::kOverwrite); 
 }; 
 
 
 
 ///oooo
  outFile->cd("s1_beta_vs_p_pim");
  for (ti=0; ti<48; ti++){
 beta_vs_p_pim[0][ti]->Write("", TObject::kOverwrite); 
 };
  outFile->cd("s2_beta_vs_p_pim");
 for (ti=0; ti<48; ti++){
 beta_vs_p_pim[1][ti]->Write("", TObject::kOverwrite); 
 }; 
  outFile->cd("s3_beta_vs_p_pim");
 for (ti=0; ti<48; ti++){
 beta_vs_p_pim[2][ti]->Write("", TObject::kOverwrite); 
 }; 
  outFile->cd("s4_beta_vs_p_pim"); 
 for (ti=0; ti<48; ti++){
 beta_vs_p_pim[3][ti]->Write("", TObject::kOverwrite); 
 }; 
  outFile->cd("s5_beta_vs_p_pim");
  for (ti=0; ti<48; ti++){
 beta_vs_p_pim[4][ti]->Write("", TObject::kOverwrite); 
 };
  outFile->cd("s6_beta_vs_p_pim");
  for (ti=0; ti<48; ti++){
 beta_vs_p_pim[5][ti]->Write("", TObject::kOverwrite); 
 };
 
outFile->cd("____________");
 /////iiii
 
  outFile->cd("s1_beta_vs_p_sim");
  for (ti=0; ti<48; ti++){
 beta_vs_p_p_sim[0][ti]->Write("", TObject::kOverwrite); 
 };
  outFile->cd("s2_beta_vs_p_sim");
 for (ti=0; ti<48; ti++){
 beta_vs_p_p_sim[1][ti]->Write("", TObject::kOverwrite);
 }; 
  outFile->cd("s3_beta_vs_p_sim");
 for (ti=0; ti<48; ti++){
 beta_vs_p_p_sim[2][ti]->Write("", TObject::kOverwrite);
 }; 
  outFile->cd("s4_beta_vs_p_sim"); 
 for (ti=0; ti<48; ti++){
 beta_vs_p_p_sim[3][ti]->Write("", TObject::kOverwrite);
 }; 
  outFile->cd("s5_beta_vs_p_sim");
  for (ti=0; ti<48; ti++){
 beta_vs_p_p_sim[4][ti]->Write("", TObject::kOverwrite); 
 };
  outFile->cd("s6_beta_vs_p_sim");
  for (ti=0; ti<48; ti++){
 beta_vs_p_p_sim[5][ti]->Write("", TObject::kOverwrite); 
 };
 
 
/////iiii
  outFile->cd("s1_beta_vs_pip_sim");
  for (ti=0; ti<48; ti++){
 beta_vs_p_pip_sim[0][ti]->Write("", TObject::kOverwrite); 
 };
  outFile->cd("s2_beta_vs_pip_sim");
 for (ti=0; ti<48; ti++){
 beta_vs_p_pip_sim[1][ti]->Write("", TObject::kOverwrite); 
 }; 
  outFile->cd("s3_beta_vs_pip_sim");
 for (ti=0; ti<48; ti++){
 beta_vs_p_pip_sim[2][ti]->Write("", TObject::kOverwrite); 
 }; 
  outFile->cd("s4_beta_vs_pip_sim"); 
 for (ti=0; ti<48; ti++){
 beta_vs_p_pip_sim[3][ti]->Write("", TObject::kOverwrite); 
 }; 
  outFile->cd("s5_beta_vs_pip_sim");
  for (ti=0; ti<48; ti++){
 beta_vs_p_pip_sim[4][ti]->Write("", TObject::kOverwrite); 
 };
  outFile->cd("s6_beta_vs_pip_sim");
  for (ti=0; ti<48; ti++){
 beta_vs_p_pip_sim[5][ti]->Write("", TObject::kOverwrite); 
 }; 
 
 
 
 ///oooo
  outFile->cd("s1_beta_vs_p_pim_sim");
  for (ti=0; ti<48; ti++){
 beta_vs_p_pim_sim[0][ti]->Write("", TObject::kOverwrite); 
 };
  outFile->cd("s2_beta_vs_p_pim_sim");
 for (ti=0; ti<48; ti++){
 beta_vs_p_pim_sim[1][ti]->Write("", TObject::kOverwrite); 
 }; 
  outFile->cd("s3_beta_vs_p_pim_sim");
 for (ti=0; ti<48; ti++){
 beta_vs_p_pim_sim[2][ti]->Write("", TObject::kOverwrite); 
 }; 
  outFile->cd("s4_beta_vs_p_pim_sim"); 
 for (ti=0; ti<48; ti++){
 beta_vs_p_pim_sim[3][ti]->Write("", TObject::kOverwrite); 
 }; 
  outFile->cd("s5_beta_vs_p_pim_sim");
  for (ti=0; ti<48; ti++){
 beta_vs_p_pim_sim[4][ti]->Write("", TObject::kOverwrite); 
 };
  outFile->cd("s6_beta_vs_p_pim_sim");
  for (ti=0; ti<48; ti++){
 beta_vs_p_pim_sim[5][ti]->Write("", TObject::kOverwrite); 
 };
 
 
///oooo
  outFile->cd("s1_time_p");
  for (ti=0; ti<48; ti++){
 time_p[0][ti]->Write("", TObject::kOverwrite); 
 };
outFile->cd("s2_time_p");
  for (ti=0; ti<48; ti++){
 time_p[1][ti]->Write("", TObject::kOverwrite); 
 };
outFile->cd("s3_time_p");
  for (ti=0; ti<48; ti++){
 time_p[2][ti]->Write("", TObject::kOverwrite); 
 };
outFile->cd("s4_time_p");
  for (ti=0; ti<48; ti++){
 time_p[3][ti]->Write("", TObject::kOverwrite); 
 };
outFile->cd("s5_time_p");
  for (ti=0; ti<48; ti++){
 time_p[4][ti]->Write("", TObject::kOverwrite); 
 };
outFile->cd("s6_time_p");
  for (ti=0; ti<48; ti++){
 time_p[5][ti]->Write("", TObject::kOverwrite); 
 };

///oooo
  outFile->cd("s1_time_pip");
  for (ti=0; ti<48; ti++){
 time_pip[0][ti]->Write("", TObject::kOverwrite); 
 };
outFile->cd("s2_time_pip");
  for (ti=0; ti<48; ti++){
 time_pip[1][ti]->Write("", TObject::kOverwrite); 
 };
outFile->cd("s3_time_pip");
  for (ti=0; ti<48; ti++){
 time_pip[2][ti]->Write("", TObject::kOverwrite); 
 };
outFile->cd("s4_time_pip");
  for (ti=0; ti<48; ti++){
 time_pip[3][ti]->Write("", TObject::kOverwrite); 
 };
outFile->cd("s5_time_pip");
  for (ti=0; ti<48; ti++){
 time_pip[4][ti]->Write("", TObject::kOverwrite); 
 };
outFile->cd("s6_time_pip");
  for (ti=0; ti<48; ti++){
 time_pip[5][ti]->Write("", TObject::kOverwrite); 
 };

///ooom
  outFile->cd("s1_time_pim");
  for (ti=0; ti<48; ti++){
 time_pim[0][ti]->Write("", TObject::kOverwrite); 
 };
outFile->cd("s2_time_pim");
  for (ti=0; ti<48; ti++){
 time_pim[1][ti]->Write("", TObject::kOverwrite); 
 };
outFile->cd("s3_time_pim");
  for (ti=0; ti<48; ti++){
 time_pim[2][ti]->Write("", TObject::kOverwrite); 
 };
outFile->cd("s4_time_pim");
  for (ti=0; ti<48; ti++){
 time_pim[3][ti]->Write("", TObject::kOverwrite); 
 };
outFile->cd("s5_time_pim");
  for (ti=0; ti<48; ti++){
 time_pim[4][ti]->Write("", TObject::kOverwrite); 
 };
outFile->cd("s6_time_pim");
  for (ti=0; ti<48; ti++){
 time_pim[5][ti]->Write("", TObject::kOverwrite); 
 };


///oooo
  outFile->cd("s1_time_p_sim");
  for (ti=0; ti<48; ti++){
 time_p_sim[0][ti]->Write("", TObject::kOverwrite); 
 };
outFile->cd("s2_time_p_sim");
  for (ti=0; ti<48; ti++){
 time_p_sim[1][ti]->Write("", TObject::kOverwrite); 
 };
outFile->cd("s3_time_p_sim");
  for (ti=0; ti<48; ti++){
 time_p_sim[2][ti]->Write("", TObject::kOverwrite); 
 };
outFile->cd("s4_time_p_sim");
  for (ti=0; ti<48; ti++){
 time_p_sim[3][ti]->Write("", TObject::kOverwrite); 
 };
outFile->cd("s5_time_p_sim");
  for (ti=0; ti<48; ti++){
 time_p_sim[4][ti]->Write("", TObject::kOverwrite); 
 };
outFile->cd("s6_time_p_sim");
  for (ti=0; ti<48; ti++){
 time_p_sim[5][ti]->Write("", TObject::kOverwrite); 
 };

///oooo
  outFile->cd("s1_time_pip_sim");
  for (ti=0; ti<48; ti++){
 time_pip_sim[0][ti]->Write("", TObject::kOverwrite); 
 };
outFile->cd("s2_time_pip_sim");
  for (ti=0; ti<48; ti++){
 time_pip_sim[1][ti]->Write("", TObject::kOverwrite); 
 };
outFile->cd("s3_time_pip_sim");
  for (ti=0; ti<48; ti++){
 time_pip_sim[2][ti]->Write("", TObject::kOverwrite); 
 };
outFile->cd("s4_time_pip_sim");
  for (ti=0; ti<48; ti++){
 time_pip_sim[3][ti]->Write("", TObject::kOverwrite); 
 };
outFile->cd("s5_time_pip_sim");
  for (ti=0; ti<48; ti++){
 time_pip_sim[4][ti]->Write("", TObject::kOverwrite); 
 };
outFile->cd("s6_time_pip_sim");
  for (ti=0; ti<48; ti++){
 time_pip_sim[5][ti]->Write("", TObject::kOverwrite); 
 };

///ooom
  outFile->cd("s1_time_pim_sim");
  for (ti=0; ti<48; ti++){
 time_pim_sim[0][ti]->Write("", TObject::kOverwrite); 
 };
outFile->cd("s2_time_pim_sim");
  for (ti=0; ti<48; ti++){
 time_pim_sim[1][ti]->Write("", TObject::kOverwrite); 
 };
outFile->cd("s3_time_pim_sim");
  for (ti=0; ti<48; ti++){
 time_pim_sim[2][ti]->Write("", TObject::kOverwrite); 
 };
outFile->cd("s4_time_pim_sim");
  for (ti=0; ti<48; ti++){
 time_pim_sim[3][ti]->Write("", TObject::kOverwrite); 
 };
outFile->cd("s5_time_pim_sim");
  for (ti=0; ti<48; ti++){
 time_pim_sim[4][ti]->Write("", TObject::kOverwrite); 
 };
outFile->cd("s6_time_pim_sim");
  for (ti=0; ti<48; ti++){
 time_pim_sim[5][ti]->Write("", TObject::kOverwrite); 
 };



  //outFile->cd("y_sc_vs_x_sc");
/*  for (jo=0; jo<6; jo++){
  for (ti=0; ti<57; ti++){
  
outFile->cd("y_sc_vs_x_sc_p");
 y_sc_vs_x_sc_p[jo][ti]->Write("", TObject::kOverwrite); 
outFile->cd("y_sc_vs_x_sc_pim");
 y_sc_vs_x_sc_pim[jo][ti]->Write("", TObject::kOverwrite);  
 outFile->cd("y_sc_vs_x_sc_pip");
 y_sc_vs_x_sc_pip[jo][ti]->Write("", TObject::kOverwrite); 
 outFile->cd("y_sc_vs_x_sc_el");
 y_sc_vs_x_sc_el[jo][ti]->Write("", TObject::kOverwrite); 
 };
  };*/
 
 
///oooo
  outFile->cd("th_vs_p");
  for (ti=0; ti<6; ti++){ 
//th_vs_p_e_1[ti]->Write("", TObject::kOverwrite);  
th_vs_p_e_2[ti]->Write("", TObject::kOverwrite); 
//th_vs_p_p_1[ti]->Write("", TObject::kOverwrite);  
th_vs_p_p_2[ti]->Write("", TObject::kOverwrite);  
//th_vs_p_pip_1[ti]->Write("", TObject::kOverwrite);  
th_vs_p_pip_2[ti]->Write("", TObject::kOverwrite);
//th_vs_p_pim_1[ti]->Write("", TObject::kOverwrite);  
th_vs_p_pim_2[ti]->Write("", TObject::kOverwrite); 
 
 
 }; 
 outFile->cd("th_vs_p_sim");
for (ti=0; ti<6; ti++){ 
//th_vs_p_e_1_sim[ti]->Write("", TObject::kOverwrite);  
th_vs_p_e_2_sim[ti]->Write("", TObject::kOverwrite); 
//th_vs_p_p_1_sim[ti]->Write("", TObject::kOverwrite);  
th_vs_p_p_2_sim[ti]->Write("", TObject::kOverwrite);  
//th_vs_p_pip_1_sim[ti]->Write("", TObject::kOverwrite);  
th_vs_p_pip_2_sim[ti]->Write("", TObject::kOverwrite);
//th_vs_p_pim_1_sim[ti]->Write("", TObject::kOverwrite);  
th_vs_p_pim_2_sim[ti]->Write("", TObject::kOverwrite); 
 
 
 };  

  outFile->cd("main_top_mm_pim_0425");
 
for (ti=0; ti<21; ti++){
h_maintop_pimismas_1[0][ti]->Write("", TObject::kOverwrite); 
h_maintop_pimismas_sim_1[0][ti]->Write("", TObject::kOverwrite); 

};


  outFile->cd("main_top_mm_pim_0475");
 
for (ti=0; ti<21; ti++){
h_maintop_pimismas_1[1][ti]->Write("", TObject::kOverwrite); 
h_maintop_pimismas_sim_1[1][ti]->Write("", TObject::kOverwrite); 
};




 outFile->cd("main_top_mm_pim_0525");
 
for (ti=0; ti<21; ti++){
h_maintop_pimismas_1[2][ti]->Write("", TObject::kOverwrite); 
h_maintop_pimismas_sim_1[2][ti]->Write("", TObject::kOverwrite); 
};



outFile->cd("main_top_mm_pim_0575");
 
for (ti=0; ti<21; ti++){
h_maintop_pimismas_1[3][ti]->Write("", TObject::kOverwrite); 
h_maintop_pimismas_sim_1[3][ti]->Write("", TObject::kOverwrite); 
};


outFile->cd("main_top_mm_pim_0625");
 
for (ti=0; ti<19; ti++){
h_maintop_pimismas_2[0][ti]->Write("", TObject::kOverwrite); 
h_maintop_pimismas_sim_2[0][ti]->Write("", TObject::kOverwrite); 
};


outFile->cd("main_top_mm_pim_0675");
 
for (ti=0; ti<19; ti++){
h_maintop_pimismas_2[1][ti]->Write("", TObject::kOverwrite); 
h_maintop_pimismas_sim_2[1][ti]->Write("", TObject::kOverwrite); 
};


outFile->cd("main_top_mm_pim_0725");
 
for (ti=0; ti<17; ti++){
h_maintop_pimismas_3[0][ti]->Write("", TObject::kOverwrite); 
h_maintop_pimismas_sim_3[0][ti]->Write("", TObject::kOverwrite); 
};


outFile->cd("main_top_mm_pim_0775");
 
for (ti=0; ti<17; ti++){
h_maintop_pimismas_3[1][ti]->Write("", TObject::kOverwrite); 
h_maintop_pimismas_sim_3[1][ti]->Write("", TObject::kOverwrite); 
};
  
  outFile->cd("main_top_mm_pim_0825");
 
for (ti=0; ti<14; ti++){
h_maintop_pimismas_4[ti]->Write("", TObject::kOverwrite); 
h_maintop_pimismas_sim_4[ti]->Write("", TObject::kOverwrite); 
};


 outFile->cd("main_top_mm_pim_0875");
 
for (ti=0; ti<13; ti++){
h_maintop_pimismas_5[ti]->Write("", TObject::kOverwrite); 
h_maintop_pimismas_sim_5[ti]->Write("", TObject::kOverwrite); 
};


 outFile->cd("main_top_mm_pim_0925");
 
for (ti=0; ti<12; ti++){
h_maintop_pimismas_6[ti]->Write("", TObject::kOverwrite); 
h_maintop_pimismas_sim_6[ti]->Write("", TObject::kOverwrite); 
};


 outFile->cd("main_top_mm_pim_0975");
 
for (ti=0; ti<10; ti++){
h_maintop_pimismas_7[ti]->Write("", TObject::kOverwrite); 
h_maintop_pimismas_sim_7[ti]->Write("", TObject::kOverwrite); 
};

 outFile->cd("main_top_mm_pim_noq2dep");
 
for (ti=0; ti<21; ti++){
h_maintop_pimismas_noq2dep[ti]->Write("", TObject::kOverwrite);
h_maintop_pimismas_noq2dep_sim[ti]->Write("", TObject::kOverwrite);
}; 
 


outFile->cd("excl_top_plots");

for (tj=0; tj<5; tj++){
 h_pim_mis_all_reg[tj]->Write("", TObject::kOverwrite);
  h_pim_mis_all_reg_sim[tj]->Write("", TObject::kOverwrite); 
  
 h_0_mis_all_reg[tj]->Write("", TObject::kOverwrite);
  h_0_mis_all_reg_sim[tj]->Write("", TObject::kOverwrite); 
    
   h_mom_all_reg[tj]->Write("", TObject::kOverwrite);
  h_mom_all_reg_sim[tj]->Write("", TObject::kOverwrite);
  
 
};

outFile->cd("main_top_plots");

for (tj=0; tj<5; tj++){
 h_pim_mis_main_top[tj]->Write("", TObject::kOverwrite);
  h_pim_mis_main_top_sim[tj]->Write("", TObject::kOverwrite); 
   
};
/*outFile->cd("miss_mass_pim_th_dep_pim");

for (tj=0; tj<5; tj++){

  
for (ti=0; ti<12; ti++){
 h_pim_mis_all_reg_th_dep_pim[tj][ti]->Write("", TObject::kOverwrite);
  h_pim_mis_all_reg_th_dep_pim_sim[tj][ti]->Write("", TObject::kOverwrite); 
 };
 };
 
 for (tj=0; tj<5; tj++){
for (ti=0; ti<12; ti++){
 h_mom_all_reg_th_dep_pim[tj][ti]->Write("", TObject::kOverwrite);
  h_mom_all_reg_th_dep_pim_sim[tj][ti]->Write("", TObject::kOverwrite); 
 };
 };
 

outFile->cd("miss_mass_pim_th_dep_pip");

for (tj=0; tj<5; tj++){
for (ti=0; ti<12; ti++){
 h_pim_mis_all_reg_th_dep_pip[tj][ti]->Write("", TObject::kOverwrite);
  h_pim_mis_all_reg_th_dep_pip_sim[tj][ti]->Write("", TObject::kOverwrite); 
 };
 };
 
 for (tj=0; tj<5; tj++){
for (ti=0; ti<12; ti++){
 h_mom_all_reg_th_dep_pip[tj][ti]->Write("", TObject::kOverwrite);
  h_mom_all_reg_th_dep_pip_sim[tj][ti]->Write("", TObject::kOverwrite); 
 };
 }; 
 
outFile->cd("miss_mass_pim_th_dep_pr");

for (tj=0; tj<5; tj++){
for (ti=0; ti<12; ti++){
 h_pim_mis_all_reg_th_dep_pr[tj][ti]->Write("", TObject::kOverwrite);
  h_pim_mis_all_reg_th_dep_pr_sim[tj][ti]->Write("", TObject::kOverwrite); 
 };
 };
 
 for (tj=0; tj<5; tj++){
for (ti=0; ti<12; ti++){
 h_mom_all_reg_th_dep_pr[tj][ti]->Write("", TObject::kOverwrite);
  h_mom_all_reg_th_dep_pr_sim[tj][ti]->Write("", TObject::kOverwrite); 
 };
 };  
 
 
 outFile->cd("miss_mass_pim_top_th_dep_pr");

for (tj=0; tj<5; tj++){
for (ti=0; ti<50; ti++){
 h_pim_mis_th_dep_pr[tj][ti]->Write("", TObject::kOverwrite);
  h_pim_mis_th_dep_pr_sim[tj][ti]->Write("", TObject::kOverwrite); 
 };
 };
 
 
   outFile->cd("miss_mass_pim_top_th_dep_arr");
 for (ti=0; ti<5; ti++) {
 for (io=0; io<12; io++) {

for (jo=0; jo<8; jo++) {

h_pim_mis_th_dep[ti][io][jo]->Write("", TObject::kOverwrite);
h_pim_mis_th_dep_sim[ti][io][jo]->Write("", TObject::kOverwrite);



};
};
};*/

 outFile->cd("inv_m_fsi");
 for (ti=0; ti<5; ti++) {
 
 h_inv_NP[ti]->Write("", TObject::kOverwrite);
h_inv_NPIp[ti]->Write("", TObject::kOverwrite);
h_inv_NPIm[ti]->Write("", TObject::kOverwrite);
 
 
 };


 
 
 
 
// TMacro *W_dep;
// W_dep = macros.W_plot();
// W_dep->Write("", TObject::kOverwrite);
// macros.W_plot_write();
 
    outFile->Write();
    file1->Write();
    file1->Close();	
      nphe_sector1->Delete();
      nphe_sector1_after->Delete();
      nphe_sector2->Delete();
      nphe_sector2_after->Delete();     
      nphe_sector3->Delete();
      nphe_sector3_after->Delete();
      nphe_sector4->Delete();
      nphe_sector4_after->Delete();      
      nphe_sector5->Delete();
      nphe_sector5_after->Delete();
      nphe_sector6->Delete();
      nphe_sector6_after->Delete();     
};
