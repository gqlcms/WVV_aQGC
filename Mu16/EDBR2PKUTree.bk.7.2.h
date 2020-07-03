//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jan  6 02:39:46 2016 by ROOT version 5.34/32
// from TTree EDBRCandidates/EDBR Candidates
// found on file: BulkGravWW750.root
//////////////////////////////////////////////////////////

#ifndef EDBR2PKUTree_h
#define EDBR2PKUTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>

//#include "DataFormats/Math/interface/deltaR.h"

#include <iostream>
#include <fstream>
using namespace std;
// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.
   const Int_t kMaxpassFilter_HBHE = 1;
   const Int_t kMaxpassFilter_HBHEIso = 1;
   const Int_t kMaxpassFilter_GlobalHalo = 1;
   const Int_t kMaxpassFilter_ECALDeadCell = 1;
   const Int_t kMaxpassFilter_GoodVtx = 1;
   const Int_t kMaxpassFilter_EEBadSc = 1;
   const Int_t kMaxpassFilter_badMuon = 1;
   const Int_t kMaxpassFilter_badChargedHadron = 1;
const Int_t kMaxpassecalBadCalibFilterUpdate = 1;

class EDBR2PKUTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   int num_of_lepton, num_of_AK8;
   // Declaration of leaf types

   Int_t           ls;
   Int_t           run;
   Int_t           nLooseEle;
   Int_t           nLooseMu;
//   Int_t           njets;
   Int_t           nbtag,nbtag_deep,num_bJet_deep;
   Int_t           num_bJet,num_sj1,num_sj2;
   Int_t           num_bJet_loose;
   Int_t           num_bJet_tight,num_bJet_loose_deep,num_bJet_tight_deep;
    Int_t           Nj8,Nj8_2,Nj4,status_1,status_2, status_3,numak4outsideak8,numak4insidephi,Nj4ex,Nj4in;
//   Double_t         jet2_pt;
//   Double_t         jet2_btag;
//   Double_t         jet3_pt;
//   Double_t         jet3_btag;

    double Mj_min, Mj_max,pt_min, pt_max, t21_max, t31_max, t41_max,  t21_min, t31_min, t41_min, PTw2divPTw1, PTw3divPTw1, PTw3divPTw2,t32_max, t42_max,t43_max,t32_min, t42_min,t43_min,PTw[3];
    double m_Jlvj,m_JJlvj,phi_max,phi_min,eta_max,eta_min,m_Jlvjj,m_JJlvjj,m_JJJlv,m_JJJlvj,m_JJJlvjj;
    double Mj_mid,pt_mid,eta_mid,phi_mid,t21_mid, t31_mid, t41_mid,t32_mid, t42_mid,t43_mid;
    
    double jetAK8puppi_dnnTop, jetAK8puppi_dnnW,jetAK8puppi_dnnH4q,jetAK8puppi_dnnTop_2, jetAK8puppi_dnnW_2,jetAK8puppi_dnnH4q_2,jetAK8puppi_dnnTop_3, jetAK8puppi_dnnW_3,jetAK8puppi_dnnH4q_3; //DeepAK8
    double jetAK8puppi_dnnDecorrTop, jetAK8puppi_dnnDecorrW,jetAK8puppi_dnnDecorrH4q,jetAK8puppi_dnnDecorrTop_2, jetAK8puppi_dnnDecorrW_2, jetAK8puppi_dnnDecorrH4q_2,jetAK8puppi_dnnDecorrTop_3, jetAK8puppi_dnnDecorrW_3, jetAK8puppi_dnnDecorrH4q_3; //Decorrelated DeepAK8
    double jetAK8puppi_dnnZ,jetAK8puppi_dnnZbb,jetAK8puppi_dnnHbb,jetAK8puppi_dnnZ_2,jetAK8puppi_dnnZbb_2,jetAK8puppi_dnnHbb_2,jetAK8puppi_dnnZ_3,jetAK8puppi_dnnZbb_3,jetAK8puppi_dnnHbb_3;
    double jetAK8puppi_dnnDecorrZ,jetAK8puppi_dnnDecorrZbb,jetAK8puppi_dnnDecorrHbb,jetAK8puppi_dnnDecorrZ_2,jetAK8puppi_dnnDecorrZbb_2,jetAK8puppi_dnnDecorrHbb_2,jetAK8puppi_dnnDecorrZ_3,jetAK8puppi_dnnDecorrZbb_3,jetAK8puppi_dnnDecorrHbb_3;
    double jetAK8puppi_dnnDecorrbb,jetAK8puppi_dnnDecorrcc,jetAK8puppi_dnnDecorrbbnog,jetAK8puppi_dnnDecorrccnog,jetAK8puppi_dnnDecorrbb_2,jetAK8puppi_dnnDecorrcc_2,jetAK8puppi_dnnDecorrbbnog_2,jetAK8puppi_dnnDecorrccnog_2,jetAK8puppi_dnnDecorrbb_3,jetAK8puppi_dnnDecorrcc_3,jetAK8puppi_dnnDecorrbbnog_3,jetAK8puppi_dnnDecorrccnog_3;
    double jetAK8puppi_dnnDecorrqcd,jetAK8puppi_dnnDecorrtop,jetAK8puppi_dnnDecorrw,jetAK8puppi_dnnDecorrz,jetAK8puppi_dnnDecorrzbb,jetAK8puppi_dnnDecorrhbb,jetAK8puppi_dnnDecorrh4q,jetAK8puppi_dnnDecorrqcd_2,jetAK8puppi_dnnDecorrtop_2,jetAK8puppi_dnnDecorrw_2,jetAK8puppi_dnnDecorrz_2,jetAK8puppi_dnnDecorrzbb_2,jetAK8puppi_dnnDecorrhbb_2,jetAK8puppi_dnnDecorrh4q_2,jetAK8puppi_dnnDecorrqcd_3,jetAK8puppi_dnnDecorrtop_3,jetAK8puppi_dnnDecorrw_3,jetAK8puppi_dnnDecorrz_3,jetAK8puppi_dnnDecorrzbb_3,jetAK8puppi_dnnDecorrhbb_3,jetAK8puppi_dnnDecorrh4q_3;

    double dnnt_max,dnnw_max,dnnh_max,dnndect_max,dnndecw_max,dnndech_max;
    double dnnt_mid,dnnw_mid,dnnh_mid,dnndect_mid,dnndecw_mid,dnndech_mid;
    double dnnt_min,dnnw_min,dnnh_min,dnndect_min,dnndecw_min,dnndech_min;
    double dnnz_max,dnnzbb_max,dnnhbb_max,dnndecz_max,dnndeczbb_max,dnndechbb_max,dnndecbb_max,dnndecbbnog_max,dnndeccc_max,dnndecccnog_max;
    double dnnz_mid,dnnzbb_mid,dnnhbb_mid,dnndecz_mid,dnndeczbb_mid,dnndechbb_mid,dnndecbb_mid,dnndecbbnog_mid,dnndeccc_mid,dnndecccnog_mid;
    double dnnz_min,dnnzbb_min,dnnhbb_min,dnndecz_min,dnndeczbb_min,dnndechbb_min,dnndecbb_min,dnndecbbnog_min,dnndeccc_min,dnndecccnog_min;
    double dnndecrawqcd_max,dnndecrawt_max,dnndecraww_max,dnndecrawz_max,dnndecrawzbb_max,dnndecrawhbb_max,dnndecrawh4q_max;
    double dnndecrawqcd_mid,dnndecrawt_mid,dnndecraww_mid,dnndecrawz_mid,dnndecrawzbb_mid,dnndecrawhbb_mid,dnndecrawh4q_mid;
    double dnndecrawqcd_min,dnndecrawt_min,dnndecraww_min,dnndecrawz_min,dnndecrawzbb_min,dnndecrawhbb_min,dnndecrawh4q_min;
    
    double ptgenwl[5],etagenwl[5],phigenwl[5],massgenwl[5],taggenwl[5];
    double ptgenzl[5],etagenzl[5],phigenzl[5],massgenzl[5],taggenzl[5];
    double ptgengl[10],etagengl[10],phigengl[10],egengl[10];
    double ptgenwf[5],etagenwf[5],phigenwf[5],massgenwf[5];
    double ptgenzf[5],etagenzf[5],phigenzf[5],massgenzf[5];
    double ptgengf[10],etagengf[10],phigengf[10],egengf[10];
    int matchedt[3],matchedtorg[3],matchedw[3],unmatched[3];
    int matchedft[3],matchedfw[3],unmatchedf[3];
    int matchedz[3],matchedg[3],matchedr[3],matchedrorg[3];
    double gent_w_tag,genantit_w_tag,mothergengf[10],mmothergengf[10],mmothergenq1f[5],mmothergenq2f[5],mmothergenq3f[5],mmothergenq4f[5],mmothergenq5f[5];
    double dRWW_R,dRWb,dRqq_W1,dRqq_W2;
    double mindRqq,maxdRqq,mindRqq_all,maxdRqq_all;
    double gent_b_pt,gent_b_phi,gent_b_eta,gent_b_mass;
    double genantit_b_pt,genantit_b_phi,genantit_b_eta,genantit_b_mass;
    double gent_w_pt,gent_w_phi,gent_w_eta,gent_w_mass;
    double genantit_w_pt,genantit_w_phi,genantit_w_eta,genantit_w_mass;
    double gent_w_q1_pt,gent_w_q1_phi,gent_w_q1_eta,gent_w_q1_e,gent_w_q1_pdg;
    double genantit_w_q1_pt,genantit_w_q1_phi,genantit_w_q1_eta,genantit_w_q1_e,genantit_w_q1_pdg;
    double gent_w_q2_pt,gent_w_q2_phi,gent_w_q2_eta,gent_w_q2_e,gent_w_q2_pdg;
    double genantit_w_q2_pt,genantit_w_q2_phi,genantit_w_q2_eta,genantit_w_q2_e,genantit_w_q2_pdg;
    double genw_q1_pt[5],genw_q1_eta[5],genw_q1_phi[5],genw_q1_e[5],genw_q1_pdg[5];
    double genw_q2_pt[5],genw_q2_eta[5],genw_q2_phi[5],genw_q2_e[5],genw_q2_pdg[5];
    double ptgenq1l[5],etagenq1l[5],phigenq1l[5],egenq1l[5];
    double ptgenq1f[5],etagenq1f[5],phigenq1f[5],egenq1f[5];
    double ptgenq2l[5],etagenq2l[5],phigenq2l[5],egenq2l[5];
    double ptgenq2f[5],etagenq2f[5],phigenq2f[5],egenq2f[5];
    double ptgenq3l[5],etagenq3l[5],phigenq3l[5],egenq3l[5];
    double ptgenq3f[5],etagenq3f[5],phigenq3f[5],egenq3f[5];
    double ptgenq4l[5],etagenq4l[5],phigenq4l[5],egenq4l[5];
    double ptgenq4f[5],etagenq4f[5],phigenq4f[5],egenq4f[5];
    double ptgenq5l[5],etagenq5l[5],phigenq5l[5],egenq5l[5];
    double ptgenq5f[5],etagenq5f[5],phigenq5f[5],egenq5f[5];
    double mothergenq1f[5],mothergenq2f[5],mothergenq3f[5],mothergenq4f[5],mothergenq5f[5];
    int matchingt[3],matchingw[3],matchingg[3],matchingz[3],matchingq[3],unmatching[3],matchingr[3],matchingtg[9][3],matchingtq[9][3];
    
    int matched0p4t[3],matched0p4torg[3],matched0p4w[3],unmatched0p4[3];
    int matched0p4ft[3],matched0p4fw[3],unmatched0p4f[3];
    int matched0p4z[3],matched0p4g[3],matched0p4r[3],matched0p4rorg[3];

    int matched0p8t[3],matched0p8torg[3],matched0p8w[3],unmatched0p8[3];
    int matched0p8ft[3],matched0p8fw[3],unmatched0p8f[3];
    int matched0p8z[3],matched0p8g[3],matched0p8r[3],matched0p8rorg[3];

    double mtbtvlep,deltaRbtvlep;
    double mtbmvlep,deltaRbmvlep,mbmvlep,deltaRantitop_wlepbm,deltaRtop_wlepbm;
    double mtblvlep,deltaRblvlep,mblvlep,deltaRantitop_wlepbl,deltaRtop_wlepbl;
    int matchnewr[4][3],nqqcr[4];
    int matchnewt[4][3],nqqct[4];
    int matchneww[4][3],nqqcw[4];
    Double_t pfMET;
    Double_t pfMETPhi;
    Double_t l_pt,l_pt_2;
    Double_t mtVlepnew;
    Double_t MTVlep;
    Double_t l_eta;
    Double_t l_phi;
    Double_t jet_pt;
    Double_t jet_pt_puppi;
    Double_t jet_eta;
    Double_t jet_phi;
    Double_t jet_mass_pruned;
    Double_t jet_mass_puppi;
    Double_t jet_mass_puppi_un;
    Double_t jet_mass_puppi_corr;
    Double_t jet_mass_puppi_casea=-99;
    Double_t jet_mass_puppi_casea_1=-99;

    Double_t jet_mass_puppi_caseb=-99;

    double pweight[446];
    TBranch        *b_pweight;
double jetAK8puppi_ptJEC_new,jetAK8puppi_ptJEC_JEC_up,jetAK8puppi_ptJEC_JEC_down,jetAK8puppi_ptJEC_JER_up,jetAK8puppi_ptJEC_JER_down;
double jetAK8puppi_ptJEC_2_new,jetAK8puppi_ptJEC_2_JEC_up,jetAK8puppi_ptJEC_2_JEC_down,jetAK8puppi_ptJEC_2_JER_up,jetAK8puppi_ptJEC_2_JER_down;
double jetAK8puppi_ptJEC_3_new,jetAK8puppi_ptJEC_3_JEC_up,jetAK8puppi_ptJEC_3_JEC_down,jetAK8puppi_ptJEC_3_JER_up,jetAK8puppi_ptJEC_3_JER_down;
double jetAK8puppi_e,jetAK8puppi_e_new,jetAK8puppi_e_JEC_up,jetAK8puppi_e_JEC_down,jetAK8puppi_e_JER_up,jetAK8puppi_e_JER_down;
double jetAK8puppi_e_2,jetAK8puppi_e_2_new,jetAK8puppi_e_2_JEC_up,jetAK8puppi_e_2_JEC_down,jetAK8puppi_e_2_JER_up,jetAK8puppi_e_2_JER_down;
double jetAK8puppi_e_3,jetAK8puppi_e_3_new,jetAK8puppi_e_3_JEC_up,jetAK8puppi_e_3_JEC_down,jetAK8puppi_e_3_JER_up,jetAK8puppi_e_3_JER_down;
double ptVlepJEC_JEC_up,yVlepJEC_JEC_up,phiVlepJEC_JEC_up,massVlepJEC_JEC_up,mtVlepJEC_JEC_up;
double ptVlepJEC_new,yVlepJEC_new,phiVlepJEC_new,massVlepJEC_new,mtVlepJEC_new;
double ptVlepJEC_JEC_down,yVlepJEC_JEC_down,phiVlepJEC_JEC_down,massVlepJEC_JEC_down,mtVlepJEC_JEC_down;
double ptVlepJEC_JER_up,yVlepJEC_JER_up,phiVlepJEC_JER_up,massVlepJEC_JER_up,mtVlepJEC_JER_up;
double ptVlepJEC_JER_down,yVlepJEC_JER_down,phiVlepJEC_JER_down,massVlepJEC_JER_down,mtVlepJEC_JER_down;
double MET_et_JEC_up,MET_et_JEC_down,MET_et_JER_up,MET_et_JER_down;
double MET_et_new;
double MET_phi_JEC_up,MET_phi_JEC_down,MET_phi_JER_up,MET_phi_JER_down;
double MET_phi_new;
    TBranch        *b_jetAK8puppi_ptJEC_new;   //!
    TBranch        *b_jetAK8puppi_ptJEC_JEC_up;   //!
    TBranch        *b_jetAK8puppi_ptJEC_JEC_down;   //!
    TBranch        *b_jetAK8puppi_ptJEC_JER_up;   //!
    TBranch        *b_jetAK8puppi_ptJEC_JER_down;   //!
    TBranch        *b_jetAK8puppi_ptJEC_2_new;   //!
    TBranch        *b_jetAK8puppi_ptJEC_2_JEC_up;   //!
    TBranch        *b_jetAK8puppi_ptJEC_2_JEC_down;   //!
    TBranch        *b_jetAK8puppi_ptJEC_2_JER_up;   //!
    TBranch        *b_jetAK8puppi_ptJEC_2_JER_down;   //!
    TBranch        *b_jetAK8puppi_ptJEC_3_new;   //!
    TBranch        *b_jetAK8puppi_ptJEC_3_JEC_up;   //!
    TBranch        *b_jetAK8puppi_ptJEC_3_JEC_down;   //!
    TBranch        *b_jetAK8puppi_ptJEC_3_JER_up;   //!
    TBranch        *b_jetAK8puppi_ptJEC_3_JER_down;   //!
    TBranch        *b_jetAK8puppi_e;   //!
    TBranch        *b_jetAK8puppi_e_new;   //!
    TBranch        *b_jetAK8puppi_e_JEC_up;   //!
    TBranch        *b_jetAK8puppi_e_JEC_down;   //!
    TBranch        *b_jetAK8puppi_e_JER_up;   //!
    TBranch        *b_jetAK8puppi_e_JER_down;   //!
    TBranch        *b_jetAK8puppi_e_2;   //!
    TBranch        *b_jetAK8puppi_e_2_new;   //!
    TBranch        *b_jetAK8puppi_e_2_JEC_up;   //!
    TBranch        *b_jetAK8puppi_e_2_JEC_down;   //!
    TBranch        *b_jetAK8puppi_e_2_JER_up;   //!
    TBranch        *b_jetAK8puppi_e_2_JER_down;   //!
    TBranch        *b_jetAK8puppi_e_3;   //!
    TBranch        *b_jetAK8puppi_e_3_new;   //!
    TBranch        *b_jetAK8puppi_e_3_JEC_up;   //!
    TBranch        *b_jetAK8puppi_e_3_JEC_down;   //!
    TBranch        *b_jetAK8puppi_e_3_JER_up;   //!
    TBranch        *b_jetAK8puppi_e_3_JER_down;   //!
    TBranch        *b_ptVlepJEC_JEC_up;   //!
    TBranch        *b_yVlepJEC_JEC_up;   //!
    TBranch        *b_phiVlepJEC_JEC_up;   //!
    TBranch        *b_massVlepJEC_JEC_up;   //!
    TBranch        *b_mtVlepJEC_JEC_up;   //!

    TBranch        *b_ptVlepJEC_new;   //!
    TBranch        *b_yVlepJEC_new;   //!
    TBranch        *b_phiVlepJEC_new;   //!
    TBranch        *b_massVlepJEC_new;   //!
    TBranch        *b_mtVlepJEC_new;   //!
    TBranch        *b_MET_phi_new;   //!

    TBranch        *b_ptVlepJEC_JEC_down;   //!
    TBranch        *b_yVlepJEC_JEC_down;   //!
    TBranch        *b_phiVlepJEC_JEC_down;   //!
    TBranch        *b_massVlepJEC_JEC_down;   //!
    TBranch        *b_mtVlepJEC_JEC_down;   //!
    TBranch        *b_ptVlepJEC_JER_up;   //!
    TBranch        *b_yVlepJEC_JER_up;   //!
    TBranch        *b_phiVlepJEC_JER_up;   //!
    TBranch        *b_massVlepJEC_JER_up;   //!
    TBranch        *b_mtVlepJEC_JER_up;   //!
    TBranch        *b_ptVlepJEC_JER_down;   //!
    TBranch        *b_yVlepJEC_JER_down;   //!
    TBranch        *b_phiVlepJEC_JER_down;   //!
    TBranch        *b_massVlepJEC_JER_down;   //!
    TBranch        *b_mtVlepJEC_JER_down;   //!
    TBranch        *b_MET_et_JEC_up;   //!
    TBranch        *b_MET_et_new;   //!
    TBranch        *b_MET_et_JEC_down;   //!
    TBranch        *b_MET_et_JER_up;   //!
    TBranch        *b_MET_et_JER_down;   //!
    TBranch        *b_MET_phi_JEC_up;   //!
    TBranch        *b_MET_phi_JEC_down;   //!
    TBranch        *b_MET_phi_JER_up;   //!
    TBranch        *b_MET_phi_JER_down;   //!


//   Double_t jet_mass_softdrop;
   Double_t W_pt;
   Double_t W_eta;
   Double_t W_phi;
   Double_t m_JJlv;
    double m_jlv,m_Jlv;
//   Double_t M_ww;
   Double_t fjet2_pt;
   Double_t fjet2_btag;
   Double_t fjet3_pt;
   Double_t fjet3_btag; 
   //JEC
//   Double_t corr_AK8;
//   Double_t corr;
   Double_t METraw_et;
   Double_t METraw_phi;
   Double_t METraw_sumEt;
   Double_t MET_et;
   Double_t MET_phi;
   Double_t MET_sumEt;
   Int_t CategoryID;
   Double_t isMatch;
   Double_t        weight,weight_tautopcorr1,weight_tautopcorr2,scale_center,scale_down,scale_up,weight_deep,weight_deep_topcorr1,weight_deep_topcorr2,weight_deep_matchedcorr1,weight_deep_matchedcorr2,weight_deep_tig_matchedcorr1,weight_deep_tig_matchedcorr2,weight_deep_los_matchedcorr1,weight_deep_los_matchedcorr2,weight_deep_lumi,weight_deep_lumi_corr1,weight_deep_lumi_corr2;
double weight_deep_lumi_nosft60,weight_deep_lumi_nosfq60,weight_deep_lumi_nosfw60,weight_deep_lumi_nosft100,weight_deep_lumi_nosfq100,weight_deep_lumi_nosfw100,weight_deep_lumi_nosft150,weight_deep_lumi_nosfq150,weight_deep_lumi_nosft200,weight_deep_lumi_nosfq200;
    double L1prefiring,L1prefiringup,L1prefiringdown;
    double weight_deep_prefire,weight_deep_prefireup,weight_deep_prefiredown;
    double weight_tig_deep,weight_tig_deep_topcorr1,weight_tig_deep_topcorr2;
    double weight_los_deep,weight_los_deep_topcorr1,weight_los_deep_topcorr2;
   Double_t        weight_nobtag;
   Double_t        IDweight;
   Double_t        IDweightISO;
   Double_t        IDweighttrk;
   Double_t        ToppTweight;
   Double_t        trigger_eff;
   Double_t        btagweight_center;
   Double_t        btagweight_up;
   Double_t        btagweight_down;
    Double_t        btagweight_center_deep;
    Double_t        btagweight_up_deep;
    Double_t        btagweight_down_deep;
    Double_t        btagweight_center_tig_deep;
    Double_t        btagweight_up_tig_deep;
    Double_t        btagweight_down_tig_deep;
    Double_t        btagweight_center_los_deep;
    Double_t        btagweight_up_los_deep;
    Double_t        btagweight_down_los_deep;
   Int_t           event;
   Int_t           nVtx;
   Double_t        ptVlep;
   Double_t        yVlep;
   Double_t        phiVlep;
   Double_t        massVlep;
   Double_t        mtVlep;
//   Double_t        massVhad_gen;
   Int_t           lep;
   Int_t           channel;
   Double_t        ptlep1;
   Double_t        ptlep2;
   Double_t        etalep1;
   Double_t        etalep2;
   Double_t        philep1;
   Double_t        philep2;
   Double_t        met;
   Double_t        metPhi;
   Double_t        theWeight;
   Double_t        nump;
   Double_t        numm;
   Double_t        npT;
   Double_t        npIT;
   Int_t           nBX;
   Double_t        triggerWeight;
   Double_t        lumiWeight;
   Double_t        pileupWeight;
   Double_t        delPhilepmet_m;
   Double_t        deltaRlepjet,deltaRjet1ak4[8]={-99,-99,-99,-99,-99,-99,-99,-99},deltaRjet1ak4_1[8]={-99,-99,-99,-99,-99,-99,-99,-99},deltaRphiak4[8]={-99,-99,-99,-99,-99,-99,-99,-99};
    
    double gen_tau_pt, gen_tau_eta, gen_tau_phi, gen_tau_e;
    double gen_tau_pt_2, gen_tau_eta_2, gen_tau_phi_2, gen_tau_e_2;
    double gen_tau_pt_3, gen_tau_eta_3, gen_tau_phi_3, gen_tau_e_3;
    
    double pttau[4],etatau[4],phitau[4],etau[4],pdgidtau[4];
    double pttau_2[4],etatau_2[4],phitau_2[4],etau_2[4],pdgidtau_2[4];
    double pttau_3[4],etatau_3[4],phitau_3[4],etau_3[4],pdgidtau_3[4];
    
    double ptq[3],etaq[3],phiq[3],eq[3],pdgidq[3];
    double ptq_2[3],etaq_2[3],phiq_2[3],eq_2[3],pdgidq_2[3];
    double ptq_3[3],etaq_3[3],phiq_3[3],eq_3[3],pdgidq_3[3];

    double gen_nele_pt, gen_nele_eta, gen_nele_phi, gen_nele_e;
    double gen_nele_pt_2, gen_nele_eta_2, gen_nele_phi_2, gen_nele_e_2;
    double gen_nmu_pt, gen_nmu_eta, gen_nmu_phi, gen_nmu_e;
    double gen_nmu_pt_2, gen_nmu_eta_2, gen_nmu_phi_2, gen_nmu_e_2;
    double gen_nele_pt_3, gen_nele_eta_3, gen_nele_phi_3, gen_nele_e_3;
    double gen_nmu_pt_3, gen_nmu_eta_3, gen_nmu_phi_3, gen_nmu_e_3;
    double gen_ntau_pt, gen_ntau_eta, gen_ntau_phi, gen_ntau_e;
    double gen_ntau_pt_2, gen_ntau_eta_2, gen_ntau_phi_2, gen_ntau_e_2;
    double gen_ntau_pt_3, gen_ntau_eta_3, gen_ntau_phi_3, gen_ntau_e_3;

    TLorentzVector *ak8sj11,*ak8sj12,*ak8sj13,*ak8sj14,*ak8sj15;//*puppi_softdropj1;
    TLorentzVector *ak8sj21,*ak8sj22,*ak8sj23,*ak8sj24,*ak8sj25;//*puppi_softdropj2;
    TLorentzVector vak8sj11,vak8sj12,vak8sj13,vak8sj14,vak8sj15;//,vpuppi_softdropj1;
    TLorentzVector vak8sj21,vak8sj22,vak8sj23,vak8sj24,vak8sj25;//,vpuppi_softdropj2;
    double ak4expt1,ak4exeta1,ak4exphi1,ak4exe1;
    double ak4expt2,ak4exeta2,ak4exphi2,ak4exe2;
    
    double deltaRjet1ak4p_dR0p8,deltaRjet2ak4p_dR0p8,deltaRjet1ak4p_dR0p4,deltaRak4mu_dR0p8,deltaRak4mu,deltaphijet1ak4,deltaphijet1ak4p_dR0p8,deltaRjet1inak4p_dR0p8,deltaphijet1inak4p_dR0p8,deltaRinak4mu_dR0p8,ak4inpt1_dR0p8,deltaPhiinak4mu_dR0p8;
    
Int_t Nj4ex_dR0p8,Nj4in_dR0p8;
    double ak4expt1_uncorr_dR0p8,ak4expt1_dR0p8,ak4exeta1_dR0p8,ak4exphi1_dR0p8,ak4exe1_dR0p8;
    double ak4expt2_dR0p8,ak4exeta2_dR0p8,ak4exphi2_dR0p8,ak4exe2_dR0p8;

    Int_t Nj4ex_dR1p0,Nj4in_dR1p0;
    double ak4expt1_dR1p0,ak4exeta1_dR1p0,ak4exphi1_dR1p0,ak4exe1_dR1p0;
    double ak4expt2_dR1p0,ak4exeta2_dR1p0,ak4exphi2_dR1p0,ak4exe2_dR1p0;
    
    Int_t Nj4ex_dR1p2,Nj4in_dR1p2;
    double ak4expt1_dR1p2,ak4exeta1_dR1p2,ak4exphi1_dR1p2,ak4exe1_dR1p2;
    double ak4expt2_dR1p2,ak4exeta2_dR1p2,ak4exphi2_dR1p2,ak4exe2_dR1p2;
    
    Int_t Nj4ex_dR0p5,Nj4in_dR0p5;
    double ak4expt1_dR0p5,ak4exeta1_dR0p5,ak4exphi1_dR0p5,ak4exe1_dR0p5;
    double ak4expt2_dR0p5,ak4exeta2_dR0p5,ak4exphi2_dR0p5,ak4exe2_dR0p5;
    
    Double_t        deltaRjet1jet2;
    Double_t deltaetajet1jet2,deltaetajet1lep,deltaetajet2lep;
   Double_t        delPhijetmet;
   Double_t        delPhijetlep;
    Double_t        deltaRlepjet_2;
    Double_t        delPhijetmet_2;
    Double_t        delPhijetlep_2;
    Double_t delPhilepmet,delPhiWjet1,delPhiWjet2,delPhiWjetclose,delPhiWjetfar;
    double dR_min_WW,dPhi_min_WW,dR_max_WW,dPhi_max_WW,dEta_max_WW,dEta_min_WW;
    double dR_j1W,dR_j2W,dR_j3W,dR_j1j2,dR_j1j3,dR_j2j3;
    double dPhi_j1W,dPhi_j2W,dPhi_j3W,dPhi_j1j2,dPhi_j1j3,dPhi_j2j3;
    double dEta_j1W,dEta_j2W,dEta_j3W,dEta_j1j2,dEta_j1j3,dEta_j2j3;
    double dR_min_WW_mass,dPhi_min_WW_mass,dR_max_WW_mass,dPhi_max_WW_mass,dEta_max_WW_mass,dEta_min_WW_mass;
    double dR_jmaxW,dR_jmidW,dR_jminW,dR_jmaxjmid,dR_jmaxjmin,dR_jmidjmin;
    double dPhi_jmaxW,dPhi_jmidW,dPhi_jminW,dPhi_jmaxjmid,dPhi_jmaxjmin,dPhi_jmidjmin;
    double dEta_jmaxW,dEta_jmidW,dEta_jminW,dEta_jmaxjmid,dEta_jmaxjmin,dEta_jmidjmin;

   Int_t           vbftag;
   Bool_t          IDLoose;
   Bool_t          IDTight;
    Bool_t          IDLoose_2;
    Bool_t          IDTight_2;
    Bool_t          IDLoose_3;
    Bool_t          IDTight_3;
   //Bool_t          isHighPt;
   //Bool_t          isHEEP;
   Double_t        trackIso;
   Double_t        muchaiso;
   Double_t        muneuiso;
   Double_t        muphoiso;
   Double_t        muPU;
   Double_t        muisolation;
   //Double_t        METraw_et;
   //Double_t        METraw_phi;
   //Double_t        METraw_sumEt;
   //Double_t        MET_et;
   //Double_t        MET_phi;
   //Double_t        MET_sumEt;
   Double_t        jetAK8puppi_pt1[4];
   Double_t        jetAK8puppi_eta1[4];
   Double_t        jetAK8puppi_mass1[4];
   Double_t        ptVlepJEC;
   Double_t        yVlepJEC;
   Double_t        phiVlepJEC;
   Double_t        massVlepJEC;
   Double_t        jetAK8puppi_sdJEC;
   Double_t        jetAK8puppi_sd;
   Double_t        jetAK8puppi_tau21;
   Double_t        jetAK8puppi_tau42;
   Double_t        jetAK8puppi_ptJEC;
   Double_t        jetAK8puppi_eta;
   Double_t        jetAK8puppi_phi;
   Double_t        jetAK8puppi_sdcorr;
   
	Double_t        jetAK8puppi_sdJEC_2;
   Double_t        jetAK8puppi_sd_2;
   Double_t        jetAK8puppi_tau21_2;
   Double_t        jetAK8puppi_tau42_2;
   Double_t        jetAK8puppi_ptJEC_2;
   Double_t        jetAK8puppi_eta_2;
   Double_t        jetAK8puppi_phi_2;
   Double_t        jetAK8puppi_sdcorr_2;
    Double_t        jetAK8puppi_sdJEC_3;
    Double_t        jetAK8puppi_sd_3;
    Double_t        jetAK8puppi_tau21_3;
    Double_t        jetAK8puppi_tau42_3;
    Double_t        jetAK8puppi_ptJEC_3;
    Double_t        jetAK8puppi_eta_3;
    Double_t        jetAK8puppi_phi_3;
    Double_t        jetAK8puppi_sdcorr_3;

    Double_t jetAK8puppi_tau1,jetAK8puppi_tau2,jetAK8puppi_tau3,jetAK8puppi_tau4,jetAK8puppi_tau1_2,jetAK8puppi_tau2_2,jetAK8puppi_tau3_2,jetAK8puppi_tau4_2,jetAK8puppi_tau1_3,jetAK8puppi_tau2_3,jetAK8puppi_tau3_3,jetAK8puppi_tau4_3;
Double_t        jetAK8puppi_tau31;
    Double_t        jetAK8puppi_tau32;
    Double_t        jetAK8puppi_tau43;
    Double_t        jetAK8puppi_tau41;
    double t1_max,t2_max,t3_max,t4_max;
    double t1_mid,t2_mid,t3_mid,t4_mid;
    double t1_min,t2_min,t3_min,t4_min;

    Double_t        jetAK8puppi_tau31_2;
    Double_t        jetAK8puppi_tau32_2;
    Double_t        jetAK8puppi_tau43_2;
    Double_t        jetAK8puppi_tau41_2;
    Double_t        jetAK8puppi_tau31_3;
    Double_t        jetAK8puppi_tau32_3;
    Double_t        jetAK8puppi_tau43_3;
    Double_t        jetAK8puppi_tau41_3;

Double_t MassVV[3];
Double_t        jet_mass_puppi_2;
Double_t        jet_mass_puppi_un_2;
Double_t        jet_mass_puppi_corr_2;
Double_t        jet_pt_puppi_2;
    
    
    Double_t        jet_mass_puppi_3;
    Double_t        jet_mass_puppi_un_3;
    Double_t        jet_mass_puppi_corr_3;
    Double_t        jet_pt_puppi_3;

Double_t massww[3];

	Double_t        candMasspuppiJEC;
    Double_t        candMasspuppiJEC_JEC_up;
    Double_t        candMasspuppiJEC_JEC_down;
    Double_t        candMasspuppiJEC_JER_up;
    Double_t        candMasspuppiJEC_JER_down;
    Double_t m_jlv_JEC_up,m_jlv_JEC_down,m_jlv_JER_up,m_jlv_JER_down;
    Double_t m_Jlv_JEC_up,m_Jlv_JEC_down,m_Jlv_JER_up,m_Jlv_JER_down;
    Double_t m_JJlv_JEC_up,m_JJlv_JEC_down,m_JJlv_JER_up,m_JJlv_JER_down;
    double m_Jlv_new,m_JJlv_new,m_jlv_new,candMasspuppiJEC_new;

//   Double_t        candMasspuppicorr;
//   Double_t        sdropJEC;
   Double_t        mtVlepJEC;
//   Double_t        delPhilepmetJEC;
//   Double_t        delPhijetmetJEC;
//   Double_t        delPhijetlepJEC;
   Int_t           HLT_Ele1;
   Int_t           HLT_Ele2;
   Int_t           HLT_Ele3;
   Int_t           HLT_Ele4;
   Int_t           HLT_Ele5;
   Int_t           HLT_Ele6;
   Int_t           HLT_Ele7;
   Int_t           HLT_Ele8;
   Int_t           HLT_Mu1;
   Int_t           HLT_Mu2;
   Int_t           HLT_Mu3;
   Int_t           HLT_Mu4;
   Int_t           HLT_Mu5;
   Int_t           HLT_Mu6;
   Int_t           HLT_Mu7;
   Int_t           HLT_Mu8;
   Int_t           HLT_Mu9;
   Int_t           HLT_Mu10;
   Int_t           HLT_Mu11;
   Int_t           HLT_Mu12;
   Bool_t          passFilter_HBHE;
   Bool_t          passFilter_HBHEIso;
   Bool_t          passFilter_GlobalHalo;
   Bool_t          passFilter_ECALDeadCell;
   Bool_t          passFilter_GoodVtx;
   Bool_t          passFilter_EEBadSc;
   Bool_t          passFilter_badMuon;
   Bool_t          passFilter_badChargedHadron;
    Bool_t          passecalBadCalibFilterUpdate;

   Int_t           ak4jet_hf[8];
   Int_t           ak4jet_pf[8];
   Double_t        ak4jet_pt[8];
    double          ak4jet_deepcsvb[8],ak4jet_deepcsvbb[8];
   Double_t        ak4jet_pt_uncorr[8];
   Double_t        ak4jet_eta[8];
   Double_t        ak4jet_phi[8];
   Double_t        ak4jet_e[8];
   Double_t        ak4jet_dr[8];
   Double_t        ak4jet_csv[8];
   Double_t        ak4jet_icsv[8];
   Double_t        deltaRAK4AK8_new[8];
   Double_t        ak4jet_IDLoose[8];
   Double_t        ak4jet_IDTight[8];

   Double_t        gen_gra_m;
   Double_t        gen_gra_pt,gen_gra_eta,gen_gra_phi;
    double gen_rad_m, gen_rad_pt, gen_rad_eta,gen_rad_phi;

   Double_t        gen_ele_pt;
   Double_t        gen_ele_eta;
   Double_t        gen_ele_phi;
   Double_t        gen_ele_e;
    Double_t        gen_ele_pt_2;
    Double_t        gen_ele_eta_2;
    Double_t        gen_ele_phi_2;
    Double_t        gen_ele_e_2;
    Double_t        gen_ele_pt_3;
    Double_t        gen_ele_eta_3;
    Double_t        gen_ele_phi_3;
    Double_t        gen_ele_e_3;

   Double_t        gen_mu_pt;
   Double_t        gen_mu_eta;
   Double_t        gen_mu_phi;
   Double_t        gen_mu_e;
    Double_t        gen_mu_pt_2;
    Double_t        gen_mu_eta_2;
    Double_t        gen_mu_phi_2;
    Double_t        gen_mu_e_2;
    Double_t        gen_mu_pt_3;
    Double_t        gen_mu_eta_3;
    Double_t        gen_mu_phi_3;
    Double_t        gen_mu_e_3;


   Double_t        genmatch_ele_pt;
   Double_t        genmatch_ele_eta;
   Double_t        genmatch_ele_phi;
   Double_t        genmatch_ele_e;
   Double_t        genmatch_ele_dr;
   Double_t        gentop_pt;
   Double_t        gentop_eta;
   Double_t        gentop_phi;
   Double_t        gentop_mass;
   Double_t        genantitop_pt;
   Double_t        genantitop_eta;
   Double_t        genantitop_phi;
   Double_t        genantitop_mass;
   Double_t        genmatch_mu_pt;
   Double_t        genmatch_mu_eta;
   Double_t        genmatch_mu_phi;
   Double_t        genmatch_mu_e;
   Double_t        genmatch_mu_dr;
   Double_t        etaGenVlep;
    Double_t        ptGenVlep;

   Double_t        phiGenVlep;
   Double_t        massGenVlep;
   Double_t        ptGenVhad;
   Double_t        etaGenVhad;
   Double_t        phiGenVhad;
   Double_t        massGenVhad;
    Double_t ptGenV_2, etaGenV_2, phiGenV_2, massGenV_2,ptGenV_3, etaGenV_3, phiGenV_3, massGenV_3;
    Double_t ptGenVhad_2, etaGenVhad_2, phiGenVhad_2, massGenVhad_2,ptGenVhad_3, etaGenVhad_3, phiGenVhad_3, massGenVhad_3;
    Double_t        ST;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_ls;   //!
   TBranch        *b_event;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_theWeight;   //!
   TBranch        *b_nump;   //!
   TBranch        *b_numm;   //!
   TBranch        *b_npT;   //!
   TBranch        *b_nLooseEle;   //!
   TBranch        *b_nLooseMu;   //!
   TBranch        *b_lep;   //!
   TBranch        *b_ptlep1;   //!
   TBranch        *b_etalep1;   //!
   TBranch        *b_philep1;   //!
   TBranch        *b_trackIso;   //!
   TBranch        *b_muisolation;   //add
   TBranch        *b_muchaiso;   //add
   TBranch        *b_muneuiso;   //add
   TBranch        *b_MET_et;   //!
   TBranch        *b_MET_phi;   //!
   TBranch        *b_ptVlepJEC;   //!
   TBranch        *b_yVlepJEC;   //!
   TBranch        *b_phiVlepJEC;   //!
   TBranch        *b_mtVlepJEC;   //!
   TBranch        *b_massVlepJEC;   //!
   TBranch        *b_jetAK8puppi_sdJEC;   //!
   TBranch        *b_jetAK8puppi_sd;   //!
   TBranch        *b_jetAK8puppi_tau21;   //!
   TBranch        *b_jetAK8puppi_tau42;   //!
   TBranch        *b_jetAK8puppi_ptJEC;   //!
   TBranch        *b_jetAK8puppi_eta;   //!
   TBranch        *b_jetAK8puppi_phi;   //!
   TBranch        *b_jetAK8puppi_sdcorr;   //!

    TBranch        *b_num_of_AK8;   //!
    TBranch        *b_num_of_lepton;   //!



    TBranch        *b_ptgenwl;   //!
    TBranch        *b_etagenwl;   //!
    TBranch        *b_phigenwl;   //!
    TBranch        *b_massgenwl;   //!
    TBranch        *b_taggenwl;   //!
    TBranch        *b_ptgenzl;   //!
    TBranch        *b_etagenzl;   //!
    TBranch        *b_phigenzl;   //!
    TBranch        *b_massgenzl;   //!
    TBranch        *b_taggenzl;   //!
    TBranch        *b_ptgengl;   //!
    TBranch        *b_etagengl;   //!
    TBranch        *b_phigengl;   //!
    TBranch        *b_egengl;   //!
    
    TBranch        *b_ptgenwf;   //!
    TBranch        *b_etagenwf;   //!
    TBranch        *b_phigenwf;   //!
    TBranch        *b_massgenwf;   //!
    TBranch        *b_ptgenzf;   //!
    TBranch        *b_etagenzf;   //!
    TBranch        *b_phigenzf;   //!
    TBranch        *b_massgenzf;   //!
    TBranch        *b_ptgengf;   //!
    TBranch        *b_etagengf;   //!
    TBranch        *b_phigengf;   //!
    TBranch        *b_egengf;   //!

    TBranch        *b_gent_b_pt;   //!
    TBranch        *b_gent_b_phi;   //!
    TBranch        *b_gent_b_eta;   //!
    TBranch        *b_gent_b_mass;   //!
    TBranch        *b_genantit_b_pt;   //!
    TBranch        *b_genantit_b_phi;   //!
    TBranch        *b_genantit_b_eta;   //!
    TBranch        *b_genantit_b_mass;   //!
    TBranch        *b_gent_w_pt;   //!
    TBranch        *b_gent_w_phi;   //!
    TBranch        *b_gent_w_eta;   //!
    TBranch        *b_gent_w_mass;   //!
    TBranch        *b_genantit_w_pt;   //!
    TBranch        *b_genantit_w_phi;   //!
    TBranch        *b_genantit_w_eta;   //!
    TBranch        *b_genantit_w_mass;   //!
    TBranch        *b_gent_w_q1_pt;   //!
    TBranch        *b_gent_w_q1_phi;   //!
    TBranch        *b_gent_w_q1_eta;   //!
    TBranch        *b_gent_w_q1_e;   //!
    TBranch        *b_gent_w_q1_pdg;   //!
    TBranch        *b_genantit_w_q1_pt;   //!
    TBranch        *b_genantit_w_q1_phi;   //!
    TBranch        *b_genantit_w_q1_eta;   //!
    TBranch        *b_genantit_w_q1_e;   //!
    TBranch        *b_genantit_w_q1_pdg;   //!
    TBranch        *b_gent_w_q2_pt;   //!
    TBranch        *b_gent_w_q2_phi;   //!
    TBranch        *b_gent_w_q2_eta;   //!
    TBranch        *b_gent_w_q2_e;   //!
    TBranch        *b_gent_w_q2_pdg;   //!
    TBranch        *b_genantit_w_q2_pt;   //!
    TBranch        *b_genantit_w_q2_phi;   //!
    TBranch        *b_genantit_w_q2_eta;   //!
    TBranch        *b_genantit_w_q2_e;   //!
    TBranch        *b_genantit_w_q2_pdg;   //!
    TBranch        *b_ptgenq1l;   //!
    TBranch        *b_etagenq1l;   //!
    TBranch        *b_phigenq1l;   //!
    TBranch        *b_egenq1l;   //!
    TBranch        *b_ptgenq1f;   //!
    TBranch        *b_etagenq1f;   //!
    TBranch        *b_phigenq1f;   //!
    TBranch        *b_egenq1f;   //!
    TBranch        *b_ptgenq2l;   //!
    TBranch        *b_etagenq2l;   //!
    TBranch        *b_phigenq2l;   //!
    TBranch        *b_egenq2l;   //!
    TBranch        *b_ptgenq2f;   //!
    TBranch        *b_etagenq2f;   //!
    TBranch        *b_phigenq2f;   //!
    TBranch        *b_egenq2f;   //!
    TBranch        *b_ptgenq3l;   //!
    TBranch        *b_etagenq3l;   //!
    TBranch        *b_phigenq3l;   //!
    TBranch        *b_egenq3l;   //!
    TBranch        *b_ptgenq3f;   //!
    TBranch        *b_etagenq3f;   //!
    TBranch        *b_phigenq3f;   //!
    TBranch        *b_egenq3f;   //!
    TBranch        *b_ptgenq4l;   //!
    TBranch        *b_etagenq4l;   //!
    TBranch        *b_phigenq4l;   //!
    TBranch        *b_egenq4l;   //!
    TBranch        *b_ptgenq4f;   //!
    TBranch        *b_etagenq4f;   //!
    TBranch        *b_phigenq4f;   //!
    TBranch        *b_egenq4f;   //!
    TBranch        *b_ptgenq5l;   //!
    TBranch        *b_etagenq5l;   //!
    TBranch        *b_phigenq5l;   //!
    TBranch        *b_egenq5l;   //!
    TBranch        *b_ptgenq5f;   //!
    TBranch        *b_etagenq5f;   //!
    TBranch        *b_phigenq5f;   //!
    TBranch        *b_egenq5f;   //!
    TBranch        *b_mothergenq1f;   //!
    TBranch        *b_mothergenq2f;   //!
    TBranch        *b_mothergenq3f;   //!
    TBranch        *b_mothergenq4f;   //!
    TBranch        *b_mothergenq5f;   //!
    TBranch        *b_gent_w_tag;   //!
    TBranch        *b_genantit_w_tag;   //!
    TBranch        *b_mothergengf;   //!
    TBranch        *b_mmothergengf;   //!
    TBranch        *b_mmothergenq1f;   //!
    TBranch        *b_mmothergenq2f;   //!
    TBranch        *b_mmothergenq3f;   //!
    TBranch        *b_mmothergenq4f;   //!
    TBranch        *b_mmothergenq5f;   //!

    TBranch        *b_genw_q1_pt;   //!
    TBranch        *b_genw_q1_eta;   //!
    TBranch        *b_genw_q1_phi;   //!
    TBranch        *b_genw_q1_e;   //!
    TBranch        *b_genw_q1_pdg;   //!
    TBranch        *b_genw_q2_pt;   //!
    TBranch        *b_genw_q2_eta;   //!
    TBranch        *b_genw_q2_phi;   //!
    TBranch        *b_genw_q2_e;   //!
    TBranch        *b_genw_q2_pdg;   //!
    
    TBranch        *b_L1prefiring;   //!
    TBranch        *b_L1prefiringup;   //!
    TBranch        *b_L1prefiringdown;   //!

   TBranch        *b_jetAK8puppi_sdJEC_2;   //!
   TBranch        *b_jetAK8puppi_sd_2;   //!
   TBranch        *b_jetAK8puppi_tau21_2;   //!
   TBranch        *b_jetAK8puppi_tau42_2;   //!
   TBranch        *b_jetAK8puppi_ptJEC_2;   //!
   TBranch        *b_jetAK8puppi_eta_2;   //!
   TBranch        *b_jetAK8puppi_phi_2;   //!
   TBranch        *b_jetAK8puppi_sdcorr_2;   //!
    TBranch        *b_jetAK8puppi_sdJEC_3;   //!
    TBranch        *b_jetAK8puppi_sd_3;   //!
    TBranch        *b_jetAK8puppi_tau21_3;   //!
    TBranch        *b_jetAK8puppi_tau42_3;   //!
    TBranch        *b_jetAK8puppi_ptJEC_3;   //!
    TBranch        *b_jetAK8puppi_eta_3;   //!
    TBranch        *b_jetAK8puppi_phi_3;   //!
    TBranch        *b_jetAK8puppi_sdcorr_3;   //!
    
    TBranch        *b_jetAK8puppi_dnnTop;   //!
    TBranch        *b_jetAK8puppi_dnnW;   //!
    TBranch        *b_jetAK8puppi_dnnH4q;   //!
    TBranch        *b_jetAK8puppi_dnnTop_2;   //!
    TBranch        *b_jetAK8puppi_dnnW_2;   //!
    TBranch        *b_jetAK8puppi_dnnH4q_2;   //!
    TBranch        *b_jetAK8puppi_dnnTop_3;   //!
    TBranch        *b_jetAK8puppi_dnnW_3;   //!
    TBranch        *b_jetAK8puppi_dnnH4q_3;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrTop;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrW;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrH4q;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrTop_2;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrW_2;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrH4q_2;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrTop_3;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrW_3;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrH4q_3;   //!

    TBranch        *b_jetAK8puppi_dnnZ;   //!
    TBranch        *b_jetAK8puppi_dnnZbb;   //!
    TBranch        *b_jetAK8puppi_dnnHbb;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrZ;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrZbb;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrHbb;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrbb;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrbbnog;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrcc;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrccnog;   //!
    TBranch        *b_jetAK8puppi_dnnZ_2;   //!
    TBranch        *b_jetAK8puppi_dnnZbb_2;   //!
    TBranch        *b_jetAK8puppi_dnnHbb_2;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrZ_2;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrZbb_2;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrHbb_2;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrbb_2;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrbbnog_2;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrcc_2;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrccnog_2;   //!
    TBranch        *b_jetAK8puppi_dnnZ_3;   //!
    TBranch        *b_jetAK8puppi_dnnZbb_3;   //!
    TBranch        *b_jetAK8puppi_dnnHbb_3;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrZ_3;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrZbb_3;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrHbb_3;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrbb_3;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrbbnog_3;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrcc_3;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrccnog_3;   //!

    TBranch        *b_jetAK8puppi_dnnDecorrqcd;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrtop;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrw;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrz;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrzbb;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrhbb;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrh4q;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrqcd_2;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrtop_2;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrw_2;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrz_2;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrzbb_2;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrhbb_2;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrh4q_2;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrqcd_3;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrtop_3;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrw_3;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrz_3;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrzbb_3;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrhbb_3;   //!
    TBranch        *b_jetAK8puppi_dnnDecorrh4q_3;   //!

   TBranch        *b_massww;   //!
   

TBranch        *b_IDLoose;   //!
    TBranch        *b_IDLoose_2;   //!
    TBranch        *b_IDLoose_3;   //!
   TBranch        *b_delPhilepmet;   //!
   TBranch        *b_deltaRlepjet;   //!
   TBranch        *b_delPhijetmet;   //!
   TBranch        *b_delPhijetlep;   //!
    TBranch        *b_deltaRlepjet_2;   //!
    TBranch        *b_delPhijetmet_2;   //!
    TBranch        *b_delPhijetlep_2;   //!

   TBranch        *b_vbftag;   //add
   TBranch        *b_candMasspuppiJEC;   //!
    TBranch        *b_m_jlv;   //!
   TBranch        *b_candMasspuppiJEC_new;   //!
    TBranch        *b_m_jlv_new;   //!
    TBranch        *b_candMasspuppiJEC_JEC_up;   //!
    TBranch        *b_candMasspuppiJEC_JEC_down;   //!
    TBranch        *b_candMasspuppiJEC_JER_up;   //!
    TBranch        *b_candMasspuppiJEC_JER_down;   //!
    TBranch        *b_m_jlv_JEC_up;   //!
    TBranch        *b_m_jlv_JEC_down;   //!
    TBranch        *b_m_jlv_JER_up;   //!
    TBranch        *b_m_jlv_JER_down;   //!


   TBranch        *b_gentop_pt; //add old
   TBranch        *b_gentop_eta; //add old
   TBranch        *b_genantitop_pt; //add old
   TBranch        *b_genantitop_eta; //add old
   TBranch        *b_HLT_Ele1;
   TBranch        *b_HLT_Ele2;
   TBranch        *b_HLT_Ele3;   //!
   TBranch        *b_HLT_Ele4;
   TBranch        *b_HLT_Ele5;
   TBranch        *b_HLT_Ele6;   //!
   TBranch        *b_HLT_Ele7;
   TBranch        *b_HLT_Ele8;
   TBranch        *b_HLT_Mu1;
   TBranch        *b_HLT_Mu2;   //!
   TBranch        *b_HLT_Mu3;
   TBranch        *b_HLT_Mu4;   //!
   TBranch        *b_HLT_Mu5;
   TBranch        *b_HLT_Mu6;   //!
   TBranch        *b_HLT_Mu7;
   TBranch        *b_HLT_Mu8;   //!
   TBranch        *b_HLT_Mu9;
   TBranch        *b_HLT_Mu10;   //!
   TBranch        *b_HLT_Mu11;
   TBranch        *b_HLT_Mu12;   //!
   TBranch        *b_ak4jet_hf;
   TBranch        *b_ak4jet_pf;
   TBranch        *b_ak4jet_pt;   //!
   TBranch        *b_ak4jet_pt_uncorr;   //!
   TBranch        *b_ak4jet_eta;   //!
   TBranch        *b_ak4jet_phi;   //!
   TBranch        *b_ak4jet_e;   //!
   TBranch        *b_ak4jet_dr;   //!
   TBranch        *b_ak4jet_csv;   //!
   TBranch        *b_ak4jet_icsv;   //!
    TBranch        *b_ak4jet_deepcsvb;   //!
    TBranch        *b_ak4jet_deepcsvbb;   //!
   TBranch        *b_ak4jet_IDLoose;   //!
   TBranch        *b_ak4jet_IDTight;   //!
    TBranch        *b_status_1;   //!
    TBranch        *b_status_2;   //!
    TBranch        *b_status_3;   //!

    TBranch        *b_pttau;   //!
    TBranch        *b_etatau;   //!
    TBranch        *b_phitau;   //!
    TBranch        *b_etau;   //!
    TBranch        *b_pdgidtau;   //!

    TBranch        *b_pttau_2;   //!
    TBranch        *b_etatau_2;   //!
    TBranch        *b_phitau_2;   //!
    TBranch        *b_etau_2;   //!
    TBranch        *b_pdgidtau_2;   //!
    
    TBranch        *b_pttau_3;   //!
    TBranch        *b_etatau_3;   //!
    TBranch        *b_phitau_3;   //!
    TBranch        *b_etau_3;   //!
    TBranch        *b_pdgidtau_3;   //!

    TBranch        *b_ptq;   //!
    TBranch        *b_etaq;   //!
    TBranch        *b_phiq;   //!
    TBranch        *b_eq;   //!
    TBranch        *b_pdgidq;   //!
    
    TBranch        *b_ptq_2;   //!
    TBranch        *b_etaq_2;   //!
    TBranch        *b_phiq_2;   //!
    TBranch        *b_eq_2;   //!
    TBranch        *b_pdgidq_2;   //!
    
    TBranch        *b_ptq_3;   //!
    TBranch        *b_etaq_3;   //!
    TBranch        *b_phiq_3;   //!
    TBranch        *b_eq_3;   //!
    TBranch        *b_pdgidq_3;   //!


   TBranch        *b_passFilter_HBHE_;   //!
   TBranch        *b_passFilter_HBHEIso_;   //!
   TBranch        *b_passFilter_GlobalHalo_;   //!
   TBranch        *b_passFilter_ECALDeadCell_;   //!
   TBranch        *b_passFilter_GoodVtx_;   //!
   TBranch        *b_passFilter_EEBadSc_;   //!
   TBranch        *b_passFilter_badMuon_;   //!
   TBranch        *b_passFilter_badChargedHadron_;   //!
    TBranch        *b_passecalBadCalibFilterUpdate_;   //!

   TBranch        *b_muphoiso;  //add
   TBranch        *b_muPU;  //add
   TBranch        *b_IDTight;   //!
    TBranch        *b_IDTight_2;   //!
    TBranch        *b_IDTight_3;   //!
   TBranch        *b_vbfeta;   //add
   TBranch        *b_vbfmjj;   //add
   TBranch        *b_nj1;   //add
   TBranch        *b_nj2;   //add
    TBranch        *b_vbfeta_2;   //add
    TBranch        *b_vbfmjj_2;   //add

   //TBranch        *b_isHighPt;   //!
   //TBranch        *b_isHEEP;   //!
   TBranch        *b_ptVlep;   //!
   TBranch        *b_yVlep;   //!
   TBranch        *b_phiVlep;   //!
   TBranch        *b_massVlep;   //!
   TBranch        *b_mtVlep;   //! 
    TBranch        *b_jetAK8puppi_tau1;   //!
    TBranch        *b_jetAK8puppi_tau2;   //!
    TBranch        *b_jetAK8puppi_tau3;   //!
    TBranch        *b_jetAK8puppi_tau4;   //!
    TBranch        *b_jetAK8puppi_tau1_2;   //!
    TBranch        *b_jetAK8puppi_tau2_2;   //!
    TBranch        *b_jetAK8puppi_tau3_2;   //!
    TBranch        *b_jetAK8puppi_tau4_2;   //!
    TBranch        *b_jetAK8puppi_tau1_3;   //!
    TBranch        *b_jetAK8puppi_tau2_3;   //!
    TBranch        *b_jetAK8puppi_tau3_3;   //!
    TBranch        *b_jetAK8puppi_tau4_3;   //!

   TBranch        *b_candMass;   //!
   TBranch        *b_channel;   //!
   TBranch        *b_ptlep2;   //!
   TBranch        *b_etalep2;   //!
   TBranch        *b_philep2;   //!
   TBranch        *b_met;   //!
   TBranch        *b_metPhi;   //!
   TBranch        *b_npIT;   //!
   TBranch        *b_nBX;   //!
   TBranch        *b_triggerWeight;   //!
   TBranch        *b_lumiWeight;   //!
   TBranch        *b_pileupWeight;   //!
   TBranch        *b_METraw_et;   //!
   TBranch        *b_METraw_phi;   //!
   TBranch        *b_METraw_sumEt;   //!
   TBranch        *b_MET_sumEt;   //!
   TBranch        *b_jetAK8puppi_pt1;   //!
   TBranch        *b_jetAK8puppi_eta1;   //!
   TBranch        *b_jetAK8puppi_mass1;   //!
   TBranch        *b_gen_gra_m;   //!
   TBranch        *b_gen_gra_pt;   //!
    TBranch        *b_gen_gra_eta;   //!
    TBranch        *b_gen_gra_phi;   //!

    TBranch        *b_gen_rad_m;   //!
    TBranch        *b_gen_rad_pt;   //!
    TBranch        *b_gen_rad_eta;   //!
    TBranch        *b_gen_rad_phi;   //!

   TBranch        *b_gen_ele_pt;   //!
   TBranch        *b_gen_ele_eta;   //!
   TBranch        *b_gen_ele_phi;   //!
   TBranch        *b_gen_ele_e;   //!
    TBranch        *b_gen_ele_pt_2;   //!
    TBranch        *b_gen_ele_eta_2;   //!
    TBranch        *b_gen_ele_phi_2;   //!
    TBranch        *b_gen_ele_e_2;   //!
    TBranch        *b_gen_ele_pt_3;   //!
    TBranch        *b_gen_ele_eta_3;   //!
    TBranch        *b_gen_ele_phi_3;   //!
    TBranch        *b_gen_ele_e_3;   //!

   TBranch        *b_gen_mu_pt;   //!
   TBranch        *b_gen_mu_eta;   //!
   TBranch        *b_gen_mu_phi;   //!
   TBranch        *b_gen_mu_e;   //!
    TBranch        *b_gen_mu_pt_2;   //!
    TBranch        *b_gen_mu_eta_2;   //!
    TBranch        *b_gen_mu_phi_2;   //!
    TBranch        *b_gen_mu_e_2;   //!
    TBranch        *b_gen_mu_pt_3;   //!
    TBranch        *b_gen_mu_eta_3;   //!
    TBranch        *b_gen_mu_phi_3;   //!
    TBranch        *b_gen_mu_e_3;   //!

    TBranch        *b_gen_tau_pt;   //!
    TBranch        *b_gen_tau_eta;   //!
    TBranch        *b_gen_tau_phi;   //!
    TBranch        *b_gen_tau_e;   //!
    TBranch        *b_gen_tau_pt_2;   //!
    TBranch        *b_gen_tau_eta_2;   //!
    TBranch        *b_gen_tau_phi_2;   //!
    TBranch        *b_gen_tau_e_2;   //!
    TBranch        *b_gen_tau_pt_3;   //!
    TBranch        *b_gen_tau_eta_3;   //!
    TBranch        *b_gen_tau_phi_3;   //!
    TBranch        *b_gen_tau_e_3;   //!

    TBranch        *b_gen_nele_pt;   //!
    TBranch        *b_gen_nele_eta;   //!
    TBranch        *b_gen_nele_phi;   //!
    TBranch        *b_gen_nele_e;   //!
    TBranch        *b_gen_nele_pt_2;   //!
    TBranch        *b_gen_nele_eta_2;   //!
    TBranch        *b_gen_nele_phi_2;   //!
    TBranch        *b_gen_nele_e_2;   //!
    TBranch        *b_gen_nele_pt_3;   //!
    TBranch        *b_gen_nele_eta_3;   //!
    TBranch        *b_gen_nele_phi_3;   //!
    TBranch        *b_gen_nele_e_3;   //!
    
    TBranch        *b_gen_nmu_pt;   //!
    TBranch        *b_gen_nmu_eta;   //!
    TBranch        *b_gen_nmu_phi;   //!
    TBranch        *b_gen_nmu_e;   //!
    TBranch        *b_gen_nmu_pt_2;   //!
    TBranch        *b_gen_nmu_eta_2;   //!
    TBranch        *b_gen_nmu_phi_2;   //!
    TBranch        *b_gen_nmu_e_2;   //!
    TBranch        *b_gen_nmu_pt_3;   //!
    TBranch        *b_gen_nmu_eta_3;   //!
    TBranch        *b_gen_nmu_phi_3;   //!
    TBranch        *b_gen_nmu_e_3;   //!
    
    TBranch        *b_gen_ntau_pt;   //!
    TBranch        *b_gen_ntau_eta;   //!
    TBranch        *b_gen_ntau_phi;   //!
    TBranch        *b_gen_ntau_e;   //!
    TBranch        *b_gen_ntau_pt_2;   //!
    TBranch        *b_gen_ntau_eta_2;   //!
    TBranch        *b_gen_ntau_phi_2;   //!
    TBranch        *b_gen_ntau_e_2;   //!
    TBranch        *b_gen_ntau_pt_3;   //!
    TBranch        *b_gen_ntau_eta_3;   //!
    TBranch        *b_gen_ntau_phi_3;   //!
    TBranch        *b_gen_ntau_e_3;   //!
    TBranch        *b_ak8sj11;   //!
    TBranch        *b_ak8sj12;   //!
    TBranch        *b_ak8sj13;   //!
    TBranch        *b_ak8sj14;   //!
    TBranch        *b_ak8sj15;   //!
    TBranch        *b_ak8sj21;   //!
    TBranch        *b_ak8sj22;   //!
    TBranch        *b_ak8sj23;   //!
    TBranch        *b_ak8sj24;   //!
    TBranch        *b_ak8sj25;   //!
    //TBranch        *b_puppi_softdropj1;   //!
    //TBranch        *b_puppi_softdropj2;   //!

   TBranch        *b_genmatch_ele_pt;   //!
   TBranch        *b_genmatch_ele_eta;   //!
   TBranch        *b_genmatch_ele_phi;   //!
   TBranch        *b_genmatch_ele_e;   //!
   TBranch        *b_genmatch_ele_dr;   //!
   TBranch        *b_genmatch_mu_pt;   //!
   TBranch        *b_genmatch_mu_eta;   //!
   TBranch        *b_genmatch_mu_phi;   //!
   TBranch        *b_genmatch_mu_e;   //!
   TBranch        *b_genmatch_mu_dr;   //!
   TBranch        *b_gentop_phi;   //add old
   TBranch        *b_gentop_mass;   //add old
   TBranch        *b_genantitop_phi;   //add old
   TBranch        *b_genantitop_mass;   //add old
   TBranch        *b_ptGenVlep;   //!
   TBranch        *b_etaGenVlep;   //!
   TBranch        *b_phiGenVlep;   //!
   TBranch        *b_massGenVlep;   //!
    
    TBranch        *b_ptGenV_2;   //!
    TBranch        *b_etaGenV_2;   //!
    TBranch        *b_phiGenV_2;   //!
    TBranch        *b_massGenV_2;   //!
    TBranch        *b_ptGenV_3;   //!
    TBranch        *b_etaGenV_3;   //!
    TBranch        *b_phiGenV_3;   //!
    TBranch        *b_massGenV_3;   //!

   TBranch        *b_ptGenVhad;   //!
   TBranch        *b_etaGenVhad;   //!
   TBranch        *b_phiGenVhad;   //!
   TBranch        *b_massGenVhad;   //!
    TBranch        *b_ptGenVhad_2;   //!
    TBranch        *b_etaGenVhad_2;   //!
    TBranch        *b_phiGenVhad_2;   //!
    TBranch        *b_massGenVhad_2;   //!
    TBranch        *b_ptGenVhad_3;   //!
    TBranch        *b_etaGenVhad_3;   //!
    TBranch        *b_phiGenVhad_3;   //!
    TBranch        *b_massGenVhad_3;   //!

/*
   TBranch        *b_njets;   //!
   TBranch        *b_nbtag;   //!
   TBranch        *b_jet2_pt;   //!
   TBranch        *b_jet2_btag;   //!
   TBranch        *b_jet3_pt;   //!
   TBranch        *b_jet3_btag;   //!
   TBranch        *b_massVhad_gen;   //!
   TBranch        *b_candMasspuppicorr;   //!
   TBranch        *b_delPhilepmetJEC;   //!
   TBranch        *b_delPhijetmetJEC;   //!
   TBranch        *b_delPhijetlepJEC;   //!
*/

   TString m_dataset;
   EDBR2PKUTree(TTree *tree=0, TString dataset="",Int_t IsData=0);

   virtual ~EDBR2PKUTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree,Int_t IsData);
   virtual void     Loop(TString channelname, Double_t XS, TTree *treew, Int_t IsData,Int_t IsMET);// channelname= "mu" or "el"
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     endJob() ;

   private:
   TTree *ExTree;
   TFile *fout; 

};

#endif

#ifdef EDBR2PKUTree_cxx
EDBR2PKUTree::EDBR2PKUTree(TTree *tree, TString dataset,Int_t IsData) : fChain(0)
{
//// if parameter tree is not specified (or zero), connect the file
//// used to generate this class and read the Tree.
//   if (tree == 0) {
//      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("BulkGravWW750.root");
//      if (!f || !f->IsOpen()) {
//         f = new TFile("BulkGravWW750.root");
//      }
//      TDirectory * dir = (TDirectory*)f->Get("BulkGravWW750.root:/treeDumper");
//      dir->GetObject("EDBRCandidates",tree);
//
//   }
   m_dataset=dataset;
   Init(tree,IsData);
}

EDBR2PKUTree::~EDBR2PKUTree()
{
   if (!fChain) return;
   //delete fChain->GetCurrentFile();
}

Int_t EDBR2PKUTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t EDBR2PKUTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void EDBR2PKUTree::Init(TTree *tree,Int_t IsData)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
    ak8sj11 = 0;
    ak8sj12 = 0;
    ak8sj13 = 0;
    ak8sj14 = 0;
    ak8sj15 = 0;
    ak8sj21 = 0;
    ak8sj22 = 0;
    ak8sj23 = 0;
    ak8sj24 = 0;
    ak8sj25 = 0;
    //puppi_softdropj1 = 0;
    //puppi_softdropj2 = 0;
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fout = new TFile(m_dataset, "RECREATE");
   ExTree = new TTree("PKUTree","PKUTree");

    //ExTree->Branch("CategoryID", &CategoryID, "CategoryID/I");
   ExTree->Branch("mtVlep", &mtVlep, "mtVlep/D");
   ExTree->Branch("mtVlepJEC", &mtVlepJEC, "mtVlepJEC/D");
    ExTree->Branch("mtbtvlep", &mtbtvlep, "mtbtvlep/D");
    
    ExTree->Branch("t1_max",&t1_max,"t1_max/D");
    ExTree->Branch("t2_max",&t2_max,"t2_max/D");
    ExTree->Branch("t3_max",&t3_max,"t3_max/D");
    ExTree->Branch("t4_max",&t4_max,"t4_max/D");
    ExTree->Branch("t1_mid",&t1_mid,"t1_mid/D");
    ExTree->Branch("t2_mid",&t2_mid,"t2_mid/D");
    ExTree->Branch("t3_mid",&t3_mid,"t3_mid/D");
    ExTree->Branch("t4_mid",&t4_mid,"t4_mid/D");
    ExTree->Branch("t1_min",&t1_min,"t1_min/D");
    ExTree->Branch("t2_min",&t2_min,"t2_min/D");
    ExTree->Branch("t3_min",&t3_min,"t3_min/D");
    ExTree->Branch("t4_min",&t4_min,"t4_min/D");
    ExTree->Branch("jetAK8puppi_tau1",&jetAK8puppi_tau1,"jetAK8puppi_tau1/D");
    ExTree->Branch("jetAK8puppi_tau2",&jetAK8puppi_tau2,"jetAK8puppi_tau2/D");
    ExTree->Branch("jetAK8puppi_tau3",&jetAK8puppi_tau3,"jetAK8puppi_tau3/D");
    ExTree->Branch("jetAK8puppi_tau4",&jetAK8puppi_tau4,"jetAK8puppi_tau4/D");
    ExTree->Branch("jetAK8puppi_tau1_2",&jetAK8puppi_tau1_2,"jetAK8puppi_tau1_2/D");
    ExTree->Branch("jetAK8puppi_tau2_2",&jetAK8puppi_tau2_2,"jetAK8puppi_tau2_2/D");
    ExTree->Branch("jetAK8puppi_tau3_2",&jetAK8puppi_tau3_2,"jetAK8puppi_tau3_2/D");
    ExTree->Branch("jetAK8puppi_tau4_2",&jetAK8puppi_tau4_2,"jetAK8puppi_tau4_2/D");
    ExTree->Branch("jetAK8puppi_tau1_3",&jetAK8puppi_tau1_3,"jetAK8puppi_tau1_3/D");
    ExTree->Branch("jetAK8puppi_tau2_3",&jetAK8puppi_tau2_3,"jetAK8puppi_tau2_3/D");
    ExTree->Branch("jetAK8puppi_tau3_3",&jetAK8puppi_tau3_3,"jetAK8puppi_tau3_3/D");
    ExTree->Branch("jetAK8puppi_tau4_3",&jetAK8puppi_tau4_3,"jetAK8puppi_tau4_3/D");

    ExTree->Branch("mtbmvlep", &mtbmvlep, "mtbmvlep/D");
    ExTree->Branch("deltaRbmvlep", &deltaRbmvlep, "deltaRbmvlep/D");
    ExTree->Branch("mbmvlep", &mbmvlep, "mbmvlep/D");
    ExTree->Branch("deltaRtop_wlepbm", &deltaRtop_wlepbm, "deltaRtop_wlepbm/D");
    ExTree->Branch("deltaRantitop_wlepbm", &deltaRantitop_wlepbm, "deltaRantitop_wlepbm/D");

    ExTree->Branch("mtblvlep", &mtblvlep, "mtblvlep/D");
    ExTree->Branch("deltaRblvlep", &deltaRblvlep, "deltaRblvlep/D");
    ExTree->Branch("mblvlep", &mblvlep, "mblvlep/D");
    ExTree->Branch("deltaRtop_wlepbl", &deltaRtop_wlepbl, "deltaRtop_wlepbl/D");
    ExTree->Branch("deltaRantitop_wlepbl", &deltaRantitop_wlepbl, "deltaRantitop_wlepbl/D");

//   ExTree->Branch("sdrop", &sdrop, "sdrop/D");
   ExTree->Branch("triggerWeight", &triggerWeight, "triggerWeight/D");
   ExTree->Branch("lumiWeight", &lumiWeight, "lumiWeight/D");
   ExTree->Branch("pileupWeight", &pileupWeight, "pileupWeight/D");

    ExTree->Branch("Nj8",&Nj8,"Nj8/I");
    ExTree->Branch("Nj8_2",&Nj8_2,"Nj8_2/I");

    ExTree->Branch("numak4outsideak8",&numak4outsideak8,"numak4outsideak8/I");
    ExTree->Branch("numak4insidephi",&numak4insidephi,"numak4insidephi/I");

    ExTree->Branch("Nj4",&Nj4,"Nj4/I");
    ExTree->Branch("Nj4ex",&Nj4ex,"Nj4ex/I");
    ExTree->Branch("Nj4in",&Nj4in,"Nj4in/I");
    //gen-level
    ExTree->Branch("dRWW_R",&dRWW_R,"dRWW_R/D");
    ExTree->Branch("dRWb",&dRWb,"dRWb/D");
    ExTree->Branch("dRqq_W1",&dRqq_W1,"dRqq_W1/D");
    ExTree->Branch("dRqq_W2",&dRqq_W2,"dRqq_W2/D");
    ExTree->Branch("mindRqq",&mindRqq,"mindRqq/D");
    ExTree->Branch("maxdRqq",&maxdRqq,"maxdRqq/D");
    ExTree->Branch("mindRqq_all",&mindRqq_all,"mindRqq_all/D");
    ExTree->Branch("maxdRqq_all",&maxdRqq_all,"maxdRqq_all/D");

   ExTree->Branch("lep",&lep,"lep/I");
   ExTree->Branch("channel",&channel,"channel/I");
   ExTree->Branch("run",&run,"run/I");
   ExTree->Branch("event",&event,"event/I");
   ExTree->Branch("lumi",&ls,"lumi/I");
   ExTree->Branch("nPV",&nVtx,"nPV/I");
   ExTree->Branch("pfMET",&pfMET,"pfMET/D");
   ExTree->Branch("pfMETPhi",&pfMETPhi,"pfMETPhi/D");
   ExTree->Branch("weight", &weight, "weight/D");
    ExTree->Branch("weight_deep_lumi", &weight_deep_lumi, "weight_deep_lumi/D");
    ExTree->Branch("weight_deep_lumi_corr1", &weight_deep_lumi_corr1, "weight_deep_lumi_corr1/D");
    ExTree->Branch("weight_deep_lumi_corr2", &weight_deep_lumi_corr2, "weight_deep_lumi_corr2/D");
        ExTree->Branch("weight_deep_lumi_nosft60",&weight_deep_lumi_nosft60,"weight_deep_lumi_nosft60/D");
        ExTree->Branch("weight_deep_lumi_nosfq60",&weight_deep_lumi_nosfq60,"weight_deep_lumi_nosfq60/D");
        ExTree->Branch("weight_deep_lumi_nosfw60",&weight_deep_lumi_nosfw60,"weight_deep_lumi_nosfw60/D");
        ExTree->Branch("weight_deep_lumi_nosft100",&weight_deep_lumi_nosft100,"weight_deep_lumi_nosft100/D");
        ExTree->Branch("weight_deep_lumi_nosfq100",&weight_deep_lumi_nosfq100,"weight_deep_lumi_nosfq100/D");
        ExTree->Branch("weight_deep_lumi_nosfw100",&weight_deep_lumi_nosfw100,"weight_deep_lumi_nosfw100/D");
        ExTree->Branch("weight_deep_lumi_nosft150",&weight_deep_lumi_nosft150,"weight_deep_lumi_nosft150/D");
        ExTree->Branch("weight_deep_lumi_nosfq150",&weight_deep_lumi_nosfq150,"weight_deep_lumi_nosfq150/D");
        ExTree->Branch("weight_deep_lumi_nosft200",&weight_deep_lumi_nosft200,"weight_deep_lumi_nosft200/D");
        ExTree->Branch("weight_deep_lumi_nosfq200",&weight_deep_lumi_nosfq200,"weight_deep_lumi_nosfq200/D");

    ExTree->Branch("weight_tautopcorr1", &weight_tautopcorr1, "weight_tautopcorr1/D");
    ExTree->Branch("weight_tautopcorr2", &weight_tautopcorr2, "weight_tautopcorr2/D");
    ExTree->Branch("weight_deep_topcorr1", &weight_deep_topcorr1, "weight_deep_topcorr1/D");
    ExTree->Branch("weight_deep_topcorr2", &weight_deep_topcorr2, "weight_deep_topcorr2/D");
    ExTree->Branch("weight_deep", &weight_deep, "weight_deep/D");
    ExTree->Branch("weight_tig_deep_topcorr1", &weight_tig_deep_topcorr1, "weight_tig_deep_topcorr1/D");
    ExTree->Branch("weight_tig_deep_topcorr2", &weight_tig_deep_topcorr2, "weight_tig_deep_topcorr2/D");
    ExTree->Branch("weight_tig_deep", &weight_tig_deep, "weight_tig_deep/D");
    ExTree->Branch("weight_los_deep_topcorr1", &weight_los_deep_topcorr1, "weight_los_deep_topcorr1/D");
    ExTree->Branch("weight_los_deep_topcorr2", &weight_los_deep_topcorr2, "weight_los_deep_topcorr2/D");
    ExTree->Branch("weight_los_deep", &weight_los_deep, "weight_los_deep/D");
    ExTree->Branch("weight_deep_matchedcorr1", &weight_deep_matchedcorr1, "weight_deep_matchedcorr1/D");
    ExTree->Branch("weight_deep_tig_matchedcorr1", &weight_deep_tig_matchedcorr1, "weight_deep_tig_matchedcorr1/D");
    ExTree->Branch("weight_deep_los_matchedcorr1", &weight_deep_los_matchedcorr1, "weight_deep_los_matchedcorr1/D");
    ExTree->Branch("weight_deep_matchedcorr2", &weight_deep_matchedcorr2, "weight_deep_matchedcorr2/D");
    ExTree->Branch("weight_deep_tig_matchedcorr2", &weight_deep_tig_matchedcorr2, "weight_deep_tig_matchedcorr2/D");
    ExTree->Branch("weight_deep_los_matchedcorr2", &weight_deep_los_matchedcorr2, "weight_deep_los_matchedcorr2/D");

    ExTree->Branch("scale_center", &scale_center, "scale_center/D");
    ExTree->Branch("scale_down", &scale_down, "scale_down/D");
    ExTree->Branch("scale_up", &scale_up, "scale_up/D");
   ExTree->Branch("weight_nobtag", &weight_nobtag, "weight_nobtag/D");
   ExTree->Branch("IDweight", &IDweight, "IDweight/D");
   ExTree->Branch("IDweightISO", &IDweightISO, "IDweightISO/D");
   ExTree->Branch("IDweighttrk", &IDweighttrk, "IDweighttrk/D");
   ExTree->Branch("ToppTweight", &ToppTweight, "ToppTweight/D");
   ExTree->Branch("trigger_eff", &trigger_eff, "trigger_eff/D");
   ExTree->Branch("btagweight_center", &btagweight_center, "btagweight_center/D");
   ExTree->Branch("btagweight_up", &btagweight_up, "btagweight_up/D");
   ExTree->Branch("btagweight_down", &btagweight_down, "btagweight_down/D");
   ExTree->Branch("isMatch", &isMatch, "isMatch/D");
   //ExTree->Branch("isHEEP", &isHEEP, "isHEEP/O");
   //ExTree->Branch("isHighPt", &isHighPt, "isHighPt/O");
   ExTree->Branch("passFilter_HBHE", &passFilter_HBHE, "passFilter_HBHE/O");
   ExTree->Branch("passFilter_GlobalHalo", &passFilter_GlobalHalo, "passFilter_GlobalHalo/O");
   ExTree->Branch("passFilter_HBHEIso", &passFilter_HBHEIso, "passFilter_HBHEIso/O");
   ExTree->Branch("passFilter_ECALDeadCell", &passFilter_ECALDeadCell, "passFilter_ECALDeadCell/O");
   ExTree->Branch("passFilter_GoodVtx", &passFilter_GoodVtx, "passFilter_GoodVtx/O");
   ExTree->Branch("passFilter_EEBadSc", &passFilter_EEBadSc, "passFilter_EEBadSc/O");
   ExTree->Branch("passFilter_badMuon", &passFilter_badMuon, "passFilter_badMuon/O");
   ExTree->Branch("passFilter_badChargedHadron", &passFilter_badChargedHadron, "passFilter_badChargedHadron/O");
    ExTree->Branch("passecalBadCalibFilterUpdate", &passecalBadCalibFilterUpdate, "passecalBadCalibFilterUpdate/O");

   ExTree->Branch("nLooseEle",&nLooseEle,"nLooseEle/I");
   ExTree->Branch("nLooseMu",&nLooseMu,"nLooseMu/I");
   ExTree->Branch("l_pt",&l_pt,"l_pt/D");
    ExTree->Branch("l_pt_2",&l_pt_2,"l_pt_2/D");

   ExTree->Branch("l_eta",&l_eta,"l_eta/D");
   ExTree->Branch("l_phi",&l_phi,"l_phi/D");
   ExTree->Branch("trackIso", &trackIso, "trackIso/D");
   ExTree->Branch("mtVlepnew",&mtVlepnew,"mtVlepnew/D");
   ExTree->Branch("MTVlep",&MTVlep,"MTVlep/D");
    ExTree->Branch("weight_deep_prefire",&weight_deep_prefire,"weight_deep_prefire/D");
    ExTree->Branch("weight_deep_prefireup",&weight_deep_prefireup,"weight_deep_prefireup/D");
    ExTree->Branch("weight_deep_prefiredown",&weight_deep_prefiredown,"weight_deep_prefiredown/D");

    if(IsData==10){
    ExTree->Branch("m_Jlv_new",&m_Jlv_new,"m_Jlv_new/D");
    ExTree->Branch("m_JJlv_new",&m_JJlv_new,"m_JJlv_new/D");
    ExTree->Branch("m_Jlv_JEC_up",&m_Jlv_JEC_up,"m_Jlv_JEC_up/D");
    ExTree->Branch("m_Jlv_JEC_down",&m_Jlv_JEC_down,"m_Jlv_JEC_down/D");
    ExTree->Branch("m_Jlv_JER_up",&m_Jlv_JER_up,"m_Jlv_JER_up/D");
    ExTree->Branch("m_Jlv_JER_down",&m_Jlv_JER_down,"m_Jlv_JER_down/D");
    ExTree->Branch("m_JJlv_JEC_up",&m_JJlv_JEC_up,"m_JJlv_JEC_up/D");
    ExTree->Branch("m_JJlv_JEC_down",&m_JJlv_JEC_down,"m_JJlv_JEC_down/D");
    ExTree->Branch("m_JJlv_JER_up",&m_JJlv_JER_up,"m_JJlv_JER_up/D");
    ExTree->Branch("m_JJlv_JER_down",&m_JJlv_JER_down,"m_JJlv_JER_down/D");

        ExTree->Branch("jetAK8puppi_ptJEC",&jetAK8puppi_ptJEC,"jetAK8puppi_ptJEC/D");
        ExTree->Branch("jetAK8puppi_ptJEC_new",&jetAK8puppi_ptJEC_new,"jetAK8puppi_ptJEC_new/D");
        ExTree->Branch("jetAK8puppi_ptJEC_JEC_up",&jetAK8puppi_ptJEC_JEC_up,"jetAK8puppi_ptJEC_JEC_up/D");
        ExTree->Branch("jetAK8puppi_ptJEC_JEC_down",&jetAK8puppi_ptJEC_JEC_down,"jetAK8puppi_ptJEC_JEC_down/D");
        ExTree->Branch("jetAK8puppi_ptJEC_JER_up",&jetAK8puppi_ptJEC_JER_up,"jetAK8puppi_ptJEC_JER_up/D");
        ExTree->Branch("jetAK8puppi_ptJEC_JER_down",&jetAK8puppi_ptJEC_JER_down,"jetAK8puppi_ptJEC_JER_down/D");
        ExTree->Branch("jetAK8puppi_ptJEC_2_new",&jetAK8puppi_ptJEC_2_new,"jetAK8puppi_ptJEC_2_new/D");
        ExTree->Branch("jetAK8puppi_ptJEC_2_JEC_up",&jetAK8puppi_ptJEC_2_JEC_up,"jetAK8puppi_ptJEC_2_JEC_up/D");
        ExTree->Branch("jetAK8puppi_ptJEC_2_JEC_down",&jetAK8puppi_ptJEC_2_JEC_down,"jetAK8puppi_ptJEC_2_JEC_down/D");
        ExTree->Branch("jetAK8puppi_ptJEC_2_JER_up",&jetAK8puppi_ptJEC_2_JER_up,"jetAK8puppi_ptJEC_2_JER_up/D");
        ExTree->Branch("jetAK8puppi_ptJEC_2_JER_down",&jetAK8puppi_ptJEC_2_JER_down,"jetAK8puppi_ptJEC_2_JER_down/D");
        ExTree->Branch("jetAK8puppi_ptJEC_3_new",&jetAK8puppi_ptJEC_3_new,"jetAK8puppi_ptJEC_3_new/D");
        ExTree->Branch("jetAK8puppi_ptJEC_3_JEC_up",&jetAK8puppi_ptJEC_3_JEC_up,"jetAK8puppi_ptJEC_3_JEC_up/D");
        ExTree->Branch("jetAK8puppi_ptJEC_3_JEC_down",&jetAK8puppi_ptJEC_3_JEC_down,"jetAK8puppi_ptJEC_3_JEC_down/D");
        ExTree->Branch("jetAK8puppi_ptJEC_3_JER_up",&jetAK8puppi_ptJEC_3_JER_up,"jetAK8puppi_ptJEC_3_JER_up/D");
        ExTree->Branch("jetAK8puppi_ptJEC_3_JER_down",&jetAK8puppi_ptJEC_3_JER_down,"jetAK8puppi_ptJEC_3_JER_down/D");
        ExTree->Branch("jetAK8puppi_e",&jetAK8puppi_e,"jetAK8puppi_e/D");
        ExTree->Branch("jetAK8puppi_e_new",&jetAK8puppi_e_new,"jetAK8puppi_e_new/D");
        ExTree->Branch("jetAK8puppi_e_JEC_up",&jetAK8puppi_e_JEC_up,"jetAK8puppi_e_JEC_up/D");
        ExTree->Branch("jetAK8puppi_e_JEC_down",&jetAK8puppi_e_JEC_down,"jetAK8puppi_e_JEC_down/D");
        ExTree->Branch("jetAK8puppi_e_JER_up",&jetAK8puppi_e_JER_up,"jetAK8puppi_e_JER_up/D");
        ExTree->Branch("jetAK8puppi_e_JER_down",&jetAK8puppi_e_JER_down,"jetAK8puppi_e_JER_down/D");
        ExTree->Branch("jetAK8puppi_e_2",&jetAK8puppi_e_2,"jetAK8puppi_e_2/D");
        ExTree->Branch("jetAK8puppi_e_2_new",&jetAK8puppi_e_2_new,"jetAK8puppi_e_2_new/D");
        ExTree->Branch("jetAK8puppi_e_2_JEC_up",&jetAK8puppi_e_2_JEC_up,"jetAK8puppi_e_2_JEC_up/D");
        ExTree->Branch("jetAK8puppi_e_2_JEC_down",&jetAK8puppi_e_2_JEC_down,"jetAK8puppi_e_2_JEC_down/D");
        ExTree->Branch("jetAK8puppi_e_2_JER_up",&jetAK8puppi_e_2_JER_up,"jetAK8puppi_e_2_JER_up/D");
        ExTree->Branch("jetAK8puppi_e_2_JER_down",&jetAK8puppi_e_2_JER_down,"jetAK8puppi_e_2_JER_down/D");
        ExTree->Branch("jetAK8puppi_e_3",&jetAK8puppi_e_3,"jetAK8puppi_e_3/D");
        ExTree->Branch("jetAK8puppi_e_3_new",&jetAK8puppi_e_3_new,"jetAK8puppi_e_3_new/D");
        ExTree->Branch("jetAK8puppi_e_3_JEC_up",&jetAK8puppi_e_3_JEC_up,"jetAK8puppi_e_3_JEC_up/D");
        ExTree->Branch("jetAK8puppi_e_3_JEC_down",&jetAK8puppi_e_3_JEC_down,"jetAK8puppi_e_3_JEC_down/D");
        ExTree->Branch("jetAK8puppi_e_3_JER_up",&jetAK8puppi_e_3_JER_up,"jetAK8puppi_e_3_JER_up/D");
        ExTree->Branch("jetAK8puppi_e_3_JER_down",&jetAK8puppi_e_3_JER_down,"jetAK8puppi_e_3_JER_down/D");
        ExTree->Branch("ptVlepJEC_JEC_up",&ptVlepJEC_JEC_up,"ptVlepJEC_JEC_up/D");
        ExTree->Branch("yVlepJEC_JEC_up",&yVlepJEC_JEC_up,"yVlepJEC_JEC_up/D");
        ExTree->Branch("phiVlepJEC_JEC_up",&phiVlepJEC_JEC_up,"phiVlepJEC_JEC_up/D");
        ExTree->Branch("massVlepJEC_JEC_up",&massVlepJEC_JEC_up,"massVlepJEC_JEC_up/D");
        ExTree->Branch("mtVlepJEC_JEC_up",&mtVlepJEC_JEC_up,"mtVlepJEC_JEC_up/D");

        ExTree->Branch("ptVlepJEC_new",&ptVlepJEC_new,"ptVlepJEC_new/D");
        ExTree->Branch("yVlepJEC_new",&yVlepJEC_new,"yVlepJEC_new/D");
        ExTree->Branch("phiVlepJEC_new",&phiVlepJEC_new,"phiVlepJEC_new/D");
        ExTree->Branch("massVlepJEC_new",&massVlepJEC_new,"massVlepJEC_new/D");
        ExTree->Branch("mtVlepJEC_new",&mtVlepJEC_new,"mtVlepJEC_new/D");
        ExTree->Branch("MET_phi_new",&MET_phi_new,"MET_phi_new/D");

        ExTree->Branch("ptVlepJEC_JEC_down",&ptVlepJEC_JEC_down,"ptVlepJEC_JEC_down/D");
        ExTree->Branch("yVlepJEC_JEC_down",&yVlepJEC_JEC_down,"yVlepJEC_JEC_down/D");
        ExTree->Branch("phiVlepJEC_JEC_down",&phiVlepJEC_JEC_down,"phiVlepJEC_JEC_down/D");
        ExTree->Branch("massVlepJEC_JEC_down",&massVlepJEC_JEC_down,"massVlepJEC_JEC_down/D");
        ExTree->Branch("mtVlepJEC_JEC_down",&mtVlepJEC_JEC_down,"mtVlepJEC_JEC_down/D");
        ExTree->Branch("ptVlepJEC_JER_up",&ptVlepJEC_JER_up,"ptVlepJEC_JER_up/D");
        ExTree->Branch("yVlepJEC_JER_up",&yVlepJEC_JER_up,"yVlepJEC_JER_up/D");
        ExTree->Branch("phiVlepJEC_JER_up",&phiVlepJEC_JER_up,"phiVlepJEC_JER_up/D");
        ExTree->Branch("massVlepJEC_JER_up",&massVlepJEC_JER_up,"massVlepJEC_JER_up/D");
        ExTree->Branch("mtVlepJEC_JER_up",&mtVlepJEC_JER_up,"mtVlepJEC_JER_up/D");
        ExTree->Branch("ptVlepJEC_JER_down",&ptVlepJEC_JER_down,"ptVlepJEC_JER_down/D");
        ExTree->Branch("yVlepJEC_JER_down",&yVlepJEC_JER_down,"yVlepJEC_JER_down/D");
        ExTree->Branch("phiVlepJEC_JER_down",&phiVlepJEC_JER_down,"phiVlepJEC_JER_down/D");
        ExTree->Branch("massVlepJEC_JER_down",&massVlepJEC_JER_down,"massVlepJEC_JER_down/D");
        ExTree->Branch("mtVlepJEC_JER_down",&mtVlepJEC_JER_down,"mtVlepJEC_JER_down/D");
        ExTree->Branch("MET_et_new",&MET_et_new,"MET_et_new/D");
        ExTree->Branch("MET_et_JEC_up",&MET_et_JEC_up,"MET_et_JEC_up/D");
        ExTree->Branch("MET_et_JEC_down",&MET_et_JEC_down,"MET_et_JEC_down/D");
        ExTree->Branch("MET_et_JER_up",&MET_et_JER_up,"MET_et_JER_up/D");
        ExTree->Branch("MET_et_JER_down",&MET_et_JER_down,"MET_et_JER_down/D");
        ExTree->Branch("MET_phi_JEC_up",&MET_phi_JEC_up,"MET_phi_JEC_up/D");
        ExTree->Branch("MET_phi_JEC_down",&MET_phi_JEC_down,"MET_phi_JEC_down/D");
        ExTree->Branch("MET_phi_JER_up",&MET_phi_JER_up,"MET_phi_JER_up/D");
        ExTree->Branch("MET_phi_JER_down",&MET_phi_JER_down,"MET_phi_JER_down/D");


    ExTree->Branch("gen_rad_m",&gen_rad_m,"gen_rad_m/D");
    ExTree->Branch("gen_rad_pt",&gen_rad_pt,"gen_rad_pt/D");
    ExTree->Branch("gen_rad_eta",&gen_rad_eta,"gen_rad_eta/D");
    ExTree->Branch("gen_rad_phi",&gen_rad_phi,"gen_rad_phi/D");

    ExTree->Branch("gen_gra_m",&gen_gra_m,"gen_gra_m/D");
    ExTree->Branch("gen_gra_pt",&gen_gra_pt,"gen_gra_pt/D");
    ExTree->Branch("gen_gra_eta",&gen_gra_eta,"gen_gra_eta/D");
    ExTree->Branch("gen_gra_phi",&gen_gra_phi,"gen_gra_phi/D");

    ExTree->Branch("gen_ele_e",&gen_ele_e,"gen_ele_e/D");
    ExTree->Branch("gen_ele_pt",&gen_ele_pt,"gen_ele_pt/D");
    ExTree->Branch("gen_ele_eta",&gen_ele_eta,"gen_ele_eta/D");
    ExTree->Branch("gen_ele_phi",&gen_ele_phi,"gen_ele_phi/D");
    ExTree->Branch("gen_ele_e_2",&gen_ele_e_2,"gen_ele_e_2/D");
    ExTree->Branch("gen_ele_pt_2",&gen_ele_pt_2,"gen_ele_pt_2/D");
    ExTree->Branch("gen_ele_eta_2",&gen_ele_eta_2,"gen_ele_eta_2/D");
    ExTree->Branch("gen_ele_phi_2",&gen_ele_phi_2,"gen_ele_phi_2/D");
    ExTree->Branch("gen_ele_e_3",&gen_ele_e_3,"gen_ele_e_3/D");
    ExTree->Branch("gen_ele_pt_3",&gen_ele_pt_3,"gen_ele_pt_3/D");
    ExTree->Branch("gen_ele_eta_3",&gen_ele_eta_3,"gen_ele_eta_3/D");
    ExTree->Branch("gen_ele_phi_3",&gen_ele_phi_3,"gen_ele_phi_3/D");

    ExTree->Branch("gen_mu_e",&gen_mu_e,"gen_mu_e/D");
    ExTree->Branch("gen_mu_pt",&gen_mu_pt,"gen_mu_pt/D");
    ExTree->Branch("gen_mu_eta",&gen_mu_eta,"gen_mu_eta/D");
    ExTree->Branch("gen_mu_phi",&gen_mu_phi,"gen_mu_phi/D");
    ExTree->Branch("gen_mu_e_2",&gen_mu_e_2,"gen_mu_e_2/D");
    ExTree->Branch("gen_mu_pt_2",&gen_mu_pt_2,"gen_mu_pt_2/D");
    ExTree->Branch("gen_mu_eta_2",&gen_mu_eta_2,"gen_mu_eta_2/D");
    ExTree->Branch("gen_mu_phi_2",&gen_mu_phi_2,"gen_mu_phi_2/D");
    ExTree->Branch("gen_mu_e_3",&gen_mu_e_3,"gen_mu_e_3/D");
    ExTree->Branch("gen_mu_pt_3",&gen_mu_pt_3,"gen_mu_pt_3/D");
    ExTree->Branch("gen_mu_eta_3",&gen_mu_eta_3,"gen_mu_eta_3/D");
    ExTree->Branch("gen_mu_phi_3",&gen_mu_phi_3,"gen_mu_phi_3/D");

    ExTree->Branch("gen_tau_e",&gen_tau_e,"gen_tau_e/D");
    ExTree->Branch("gen_tau_pt",&gen_tau_pt,"gen_tau_pt/D");
    ExTree->Branch("gen_tau_eta",&gen_tau_eta,"gen_tau_eta/D");
    ExTree->Branch("gen_tau_phi",&gen_tau_phi,"gen_tau_phi/D");
    ExTree->Branch("gen_tau_e_2",&gen_tau_e_2,"gen_tau_e_2/D");
    ExTree->Branch("gen_tau_pt_2",&gen_tau_pt_2,"gen_tau_pt_2/D");
    ExTree->Branch("gen_tau_eta_2",&gen_tau_eta_2,"gen_tau_eta_2/D");
    ExTree->Branch("gen_tau_phi_2",&gen_tau_phi_2,"gen_tau_phi_2/D");
    ExTree->Branch("gen_tau_e_3",&gen_tau_e_3,"gen_tau_e_3/D");
    ExTree->Branch("gen_tau_pt_3",&gen_tau_pt_3,"gen_tau_pt_3/D");
    ExTree->Branch("gen_tau_eta_3",&gen_tau_eta_3,"gen_tau_eta_3/D");
    ExTree->Branch("gen_tau_phi_3",&gen_tau_phi_3,"gen_tau_phi_3/D");

    ExTree->Branch("gen_nele_e",&gen_nele_e,"gen_nele_e/D");
    ExTree->Branch("gen_nele_pt",&gen_nele_pt,"gen_nele_pt/D");
    ExTree->Branch("gen_nele_eta",&gen_nele_eta,"gen_nele_eta/D");
    ExTree->Branch("gen_nele_phi",&gen_nele_phi,"gen_nele_phi/D");
    ExTree->Branch("gen_nele_e_2",&gen_nele_e_2,"gen_nele_e_2/D");
    ExTree->Branch("gen_nele_pt_2",&gen_nele_pt_2,"gen_nele_pt_2/D");
    ExTree->Branch("gen_nele_eta_2",&gen_nele_eta_2,"gen_nele_eta_2/D");
    ExTree->Branch("gen_nele_phi_2",&gen_nele_phi_2,"gen_nele_phi_2/D");
    ExTree->Branch("gen_nele_e_3",&gen_nele_e_3,"gen_nele_e_3/D");
    ExTree->Branch("gen_nele_pt_3",&gen_nele_pt_3,"gen_nele_pt_3/D");
    ExTree->Branch("gen_nele_eta_3",&gen_nele_eta_3,"gen_nele_eta_3/D");
    ExTree->Branch("gen_nele_phi_3",&gen_nele_phi_3,"gen_nele_phi_3/D");
    
    ExTree->Branch("gen_nmu_e",&gen_nmu_e,"gen_nmu_e/D");
    ExTree->Branch("gen_nmu_pt",&gen_nmu_pt,"gen_nmu_pt/D");
    ExTree->Branch("gen_nmu_eta",&gen_nmu_eta,"gen_nmu_eta/D");
    ExTree->Branch("gen_nmu_phi",&gen_nmu_phi,"gen_nmu_phi/D");
    ExTree->Branch("gen_nmu_e_2",&gen_nmu_e_2,"gen_nmu_e_2/D");
    ExTree->Branch("gen_nmu_pt_2",&gen_nmu_pt_2,"gen_nmu_pt_2/D");
    ExTree->Branch("gen_nmu_eta_2",&gen_nmu_eta_2,"gen_nmu_eta_2/D");
    ExTree->Branch("gen_nmu_phi_2",&gen_nmu_phi_2,"gen_nmu_phi_2/D");
    ExTree->Branch("gen_nmu_e_3",&gen_nmu_e_3,"gen_nmu_e_3/D");
    ExTree->Branch("gen_nmu_pt_3",&gen_nmu_pt_3,"gen_nmu_pt_3/D");
    ExTree->Branch("gen_nmu_eta_3",&gen_nmu_eta_3,"gen_nmu_eta_3/D");
    ExTree->Branch("gen_nmu_phi_3",&gen_nmu_phi_3,"gen_nmu_phi_3/D");
    
    ExTree->Branch("gen_ntau_e",&gen_ntau_e,"gen_ntau_e/D");
    ExTree->Branch("gen_ntau_pt",&gen_ntau_pt,"gen_ntau_pt/D");
    ExTree->Branch("gen_ntau_eta",&gen_ntau_eta,"gen_ntau_eta/D");
    ExTree->Branch("gen_ntau_phi",&gen_ntau_phi,"gen_ntau_phi/D");
    ExTree->Branch("gen_ntau_e_2",&gen_ntau_e_2,"gen_ntau_e_2/D");
    ExTree->Branch("gen_ntau_pt_2",&gen_ntau_pt_2,"gen_ntau_pt_2/D");
    ExTree->Branch("gen_ntau_eta_2",&gen_ntau_eta_2,"gen_ntau_eta_2/D");
    ExTree->Branch("gen_ntau_phi_2",&gen_ntau_phi_2,"gen_ntau_phi_2/D");
    ExTree->Branch("gen_ntau_e_3",&gen_ntau_e_3,"gen_ntau_e_3/D");
    ExTree->Branch("gen_ntau_pt_3",&gen_ntau_pt_3,"gen_ntau_pt_3/D");
    ExTree->Branch("gen_ntau_eta_3",&gen_ntau_eta_3,"gen_ntau_eta_3/D");
    ExTree->Branch("gen_ntau_phi_3",&gen_ntau_phi_3,"gen_ntau_phi_3/D");
        ExTree->Branch("pttau",pttau,"pttau[4]/D");
        ExTree->Branch("phitau",phitau,"phitau[4]/D");
        ExTree->Branch("etatau",etatau,"etatau[4]/D");
        ExTree->Branch("etau",etau,"etau[4]/D");
        ExTree->Branch("pdgidtau",pdgidtau,"pdgidtau[4]/D");
        ExTree->Branch("pttau_2",pttau_2,"pttau_2[4]/D");
        ExTree->Branch("phitau_2",phitau_2,"phitau_2[4]/D");
        ExTree->Branch("etatau_2",etatau_2,"etatau_2[4]/D");
        ExTree->Branch("etau_2",etau_2,"etau_2[4]/D");
        ExTree->Branch("pdgidtau_2",pdgidtau_2,"pdgidtau_2[4]/D");
        ExTree->Branch("pttau_3",pttau_3,"pttau_3[4]/D");
        ExTree->Branch("phitau_3",phitau_3,"phitau_3[4]/D");
        ExTree->Branch("etatau_3",etatau_3,"etatau_3[4]/D");
        ExTree->Branch("etau_3",etau_3,"etau_3[4]/D");
        ExTree->Branch("pdgidtau_3",pdgidtau_3,"pdgidtau_3[4]/D");
        
        ExTree->Branch("ptq",ptq,"ptq[3]/D");
        ExTree->Branch("phiq",phiq,"phiq[3]/D");
        ExTree->Branch("etaq",etaq,"etaq[3]/D");
        ExTree->Branch("eq",eq,"eq[3]/D");
        ExTree->Branch("pdgidq",pdgidq,"pdgidq[3]/D");
        ExTree->Branch("ptq_2",ptq_2,"ptq_2[3]/D");
        ExTree->Branch("phiq_2",phiq_2,"phiq_2[3]/D");
        ExTree->Branch("etaq_2",etaq_2,"etaq_2[3]/D");
        ExTree->Branch("eq_2",eq_2,"eq_2[3]/D");
        ExTree->Branch("pdgidq_2",pdgidq_2,"pdgidq_2[3]/D");
        ExTree->Branch("ptq_3",ptq_3,"ptq_3[3]/D");
        ExTree->Branch("phiq_3",phiq_3,"phiq_3[3]/D");
        ExTree->Branch("etaq_3",etaq_3,"etaq_3[3]/D");
        ExTree->Branch("eq_3",eq_3,"eq_3[3]/D");
        ExTree->Branch("pdgidq_3",pdgidq_3,"pdgidq_3[3]/D");

    }
    if(IsData>0){
        ExTree->Branch("matchedt",matchedt,"matchedt[3]/I");
        ExTree->Branch("matchedtorg",matchedtorg,"matchedtorg[3]/I");
        ExTree->Branch("matchedw",matchedw,"matchedw[3]/I");
        ExTree->Branch("unmatched",unmatched,"unmatched[3]/I");
        ExTree->Branch("matchedft",matchedft,"matchedft[3]/I");
        ExTree->Branch("matchedfw",matchedfw,"matchedfw[3]/I");
        ExTree->Branch("unmatchedf",unmatchedf,"unmatchedf[3]/I");
        ExTree->Branch("matchedz",matchedz,"matchedz[3]/I");
        ExTree->Branch("matchedg",matchedg,"matchedg[3]/I");
        ExTree->Branch("matchedr",matchedr,"matchedr[3]/I");
        ExTree->Branch("matchedrorg",matchedrorg,"matchedrorg[3]/I");

        ExTree->Branch("matched0p4t",matched0p4t,"matched0p4t[3]/I");
        ExTree->Branch("matched0p4torg",matched0p4torg,"matched0p4torg[3]/I");
        ExTree->Branch("matched0p4w",matched0p4w,"matched0p4w[3]/I");
        ExTree->Branch("unmatched0p4",unmatched0p4,"unmatched0p4[3]/I");
        ExTree->Branch("matched0p4ft",matched0p4ft,"matched0p4ft[3]/I");
        ExTree->Branch("matched0p4fw",matched0p4fw,"matched0p4fw[3]/I");
        ExTree->Branch("unmatched0p4f",unmatched0p4f,"unmatched0p4f[3]/I");
        ExTree->Branch("matched0p4z",matched0p4z,"matched0p4z[3]/I");
        ExTree->Branch("matched0p4g",matched0p4g,"matched0p4g[3]/I");
        ExTree->Branch("matched0p4r",matched0p4r,"matched0p4r[3]/I");
        ExTree->Branch("matched0p4rorg",matched0p4rorg,"matched0p4rorg[3]/I");

        ExTree->Branch("matched0p8t",matched0p8t,"matched0p8t[3]/I");
        ExTree->Branch("matched0p8torg",matched0p8torg,"matched0p8torg[3]/I");
        ExTree->Branch("matched0p8w",matched0p8w,"matched0p8w[3]/I");
        ExTree->Branch("unmatched0p8",unmatched0p8,"unmatched0p8[3]/I");
        ExTree->Branch("matched0p8ft",matched0p8ft,"matched0p8ft[3]/I");
        ExTree->Branch("matched0p8fw",matched0p8fw,"matched0p8fw[3]/I");
        ExTree->Branch("unmatched0p8f",unmatched0p8f,"unmatched0p8f[3]/I");
        ExTree->Branch("matched0p8z",matched0p8z,"matched0p8z[3]/I");
        ExTree->Branch("matched0p8g",matched0p8g,"matched0p8g[3]/I");
        ExTree->Branch("matched0p8r",matched0p8r,"matched0p8r[3]/I");
        ExTree->Branch("matched0p8rorg",matched0p8rorg,"matched0p8rorg[3]/I");
        ExTree->Branch("pweight",pweight,"pweight[446]/D");
        ExTree->Branch("taggenwl",taggenwl,"taggenwl[5]/D");
        ExTree->Branch("ptgenwl",ptgenwl,"ptgenwl[5]/D");
        ExTree->Branch("etagenwl",etagenwl,"etagenwl[5]/D");
        ExTree->Branch("phigenwl",phigenwl,"phigenwl[5]/D");
        ExTree->Branch("massgenwl",massgenwl,"massgenwl[5]/D");

/*
        ExTree->Branch("matched0p8t",&matched0p8t,"matched0p8t/I");
        ExTree->Branch("matched0p8w",&matched0p8w,"matched0p8w/I");
        ExTree->Branch("unmatched0p8",&unmatched0p8,"unmatched0p8/I");
        ExTree->Branch("matchedf0p8t",&matchedf0p8t,"matchedf0p8t/I");
        ExTree->Branch("matchedf0p8w",&matchedf0p8w,"matchedf0p8w/I");
        ExTree->Branch("unmatchedf0p8",&unmatchedf0p8,"unmatchedf0p8/I");*/
        
        /*ExTree->Branch("ptgenwl",ptgenwl,"ptgenwl[5]/D");
        ExTree->Branch("etagenwl",etagenwl,"etagenwl[5]/D");
        ExTree->Branch("phigenwl",phigenwl,"phigenwl[5]/D");
        ExTree->Branch("massgenwl",massgenwl,"massgenwl[5]/D");
        ExTree->Branch("taggenwl",taggenwl,"taggenwl[5]/D");
        ExTree->Branch("ptgenzl",ptgenzl,"ptgenzl[5]/D");
        ExTree->Branch("etagenzl",etagenzl,"etagenzl[5]/D");
        ExTree->Branch("phigenzl",phigenzl,"phigenzl[5]/D");
        ExTree->Branch("massgenzl",massgenzl,"massgenzl[5]/D");
        ExTree->Branch("taggenzl",taggenzl,"taggenzl[5]/D");
        ExTree->Branch("ptgengl",ptgengl,"ptgengl[10]/D");
        ExTree->Branch("etagengl",etagengl,"etagengl[10]/D");
        ExTree->Branch("phigengl",phigengl,"phigengl[10]/D");
        ExTree->Branch("egengl",egengl,"egengl[10]/D");
        
        ExTree->Branch("ptgenwf",ptgenwf,"ptgenwf[5]/D");
        ExTree->Branch("etagenwf",etagenwf,"etagenwf[5]/D");
        ExTree->Branch("phigenwf",phigenwf,"phigenwf[5]/D");
        ExTree->Branch("massgenwf",massgenwf,"massgenwf[5]/D");
        ExTree->Branch("ptgenzf",ptgenzf,"ptgenzf[5]/D");
        ExTree->Branch("etagenzf",etagenzf,"etagenzf[5]/D");
        ExTree->Branch("phigenzf",phigenzf,"phigenzf[5]/D");
        ExTree->Branch("massgenzf",massgenzf,"massgenzf[5]/D");
        ExTree->Branch("ptgengf",ptgengf,"ptgengf[10]/D");
        ExTree->Branch("etagengf",etagengf,"etagengf[10]/D");
        ExTree->Branch("phigengf",phigengf,"phigengf[10]/D");
        ExTree->Branch("egengf",egengf,"egengf[10]/D");
        
        ExTree->Branch("gentop_pt",&gentop_pt,"gentop_pt/D");
        ExTree->Branch("gentop_phi",&gentop_phi,"gentop_phi/D");
        ExTree->Branch("gentop_eta",&gentop_eta,"gentop_eta/D");
        ExTree->Branch("gentop_mass",&gentop_mass,"gentop_mass/D");
        ExTree->Branch("genantitop_pt",&genantitop_pt,"genantitop_pt/D");
        ExTree->Branch("genantitop_phi",&genantitop_phi,"genantitop_phi/D");
        ExTree->Branch("genantitop_eta",&genantitop_eta,"genantitop_eta/D");
        ExTree->Branch("genantitop_mass",&genantitop_mass,"genantitop_mass/D");

        ExTree->Branch("gent_b_pt",&gent_b_pt,"gent_b_pt/D");
        ExTree->Branch("gent_b_phi",&gent_b_phi,"gent_b_phi/D");
        ExTree->Branch("gent_b_eta",&gent_b_eta,"gent_b_eta/D");
        ExTree->Branch("gent_b_mass",&gent_b_mass,"gent_b_mass/D");
        ExTree->Branch("genantit_b_pt",&genantit_b_pt,"genantit_b_pt/D");
        ExTree->Branch("genantit_b_phi",&genantit_b_phi,"genantit_b_phi/D");
        ExTree->Branch("genantit_b_eta",&genantit_b_eta,"genantit_b_eta/D");
        ExTree->Branch("genantit_b_mass",&genantit_b_mass,"genantit_b_mass/D");
        ExTree->Branch("gent_w_pt",&gent_w_pt,"gent_w_pt/D");
        ExTree->Branch("gent_w_phi",&gent_w_phi,"gent_w_phi/D");
        ExTree->Branch("gent_w_eta",&gent_w_eta,"gent_w_eta/D");
        ExTree->Branch("gent_w_mass",&gent_w_mass,"gent_w_mass/D");
        ExTree->Branch("genantit_w_pt",&genantit_w_pt,"genantit_w_pt/D");
        ExTree->Branch("genantit_w_phi",&genantit_w_phi,"genantit_w_phi/D");
        ExTree->Branch("genantit_w_eta",&genantit_w_eta,"genantit_w_eta/D");
        ExTree->Branch("genantit_w_mass",&genantit_w_mass,"genantit_w_mass/D");
        ExTree->Branch("gent_w_q1_pt",&gent_w_q1_pt,"gent_w_q1_pt/D");
        ExTree->Branch("gent_w_q1_phi",&gent_w_q1_phi,"gent_w_q1_phi/D");
        ExTree->Branch("gent_w_q1_eta",&gent_w_q1_eta,"gent_w_q1_eta/D");
        ExTree->Branch("gent_w_q1_e",&gent_w_q1_e,"gent_w_q1_e/D");
        ExTree->Branch("gent_w_q1_pdg",&gent_w_q1_pdg,"gent_w_q1_pdg/D");
        ExTree->Branch("genantit_w_q1_pt",&genantit_w_q1_pt,"genantit_w_q1_pt/D");
        ExTree->Branch("genantit_w_q1_phi",&genantit_w_q1_phi,"genantit_w_q1_phi/D");
        ExTree->Branch("genantit_w_q1_eta",&genantit_w_q1_eta,"genantit_w_q1_eta/D");
        ExTree->Branch("genantit_w_q1_e",&genantit_w_q1_e,"genantit_w_q1_e/D");
        ExTree->Branch("genantit_w_q1_pdg",&genantit_w_q1_pdg,"genantit_w_q1_pdg/D");
        ExTree->Branch("gent_w_q2_pt",&gent_w_q2_pt,"gent_w_q2_pt/D");
        ExTree->Branch("gent_w_q2_phi",&gent_w_q2_phi,"gent_w_q2_phi/D");
        ExTree->Branch("gent_w_q2_eta",&gent_w_q2_eta,"gent_w_q2_eta/D");
        ExTree->Branch("gent_w_q2_e",&gent_w_q2_e,"gent_w_q2_e/D");
        ExTree->Branch("gent_w_q2_pdg",&gent_w_q2_pdg,"gent_w_q2_pdg/D");
        ExTree->Branch("genantit_w_q2_pt",&genantit_w_q2_pt,"genantit_w_q2_pt/D");
        ExTree->Branch("genantit_w_q2_phi",&genantit_w_q2_phi,"genantit_w_q2_phi/D");
        ExTree->Branch("genantit_w_q2_eta",&genantit_w_q2_eta,"genantit_w_q2_eta/D");
        ExTree->Branch("genantit_w_q2_e",&genantit_w_q2_e,"genantit_w_q2_e/D");
        ExTree->Branch("genantit_w_q2_pdg",&genantit_w_q2_pdg,"genantit_w_q2_pdg/D");
        ExTree->Branch("ptgenq1l",ptgenq1l,"ptgenq1l[5]/D");
        ExTree->Branch("etagenq1l",etagenq1l,"etagenq1l[5]/D");
        ExTree->Branch("phigenq1l",phigenq1l,"phigenq1l[5]/D");
        ExTree->Branch("egenq1l",egenq1l,"egenq1l[5]/D");
        ExTree->Branch("ptgenq1f",ptgenq1f,"ptgenq1f[5]/D");
        ExTree->Branch("etagenq1f",etagenq1f,"etagenq1f[5]/D");
        ExTree->Branch("phigenq1f",phigenq1f,"phigenq1f[5]/D");
        ExTree->Branch("egenq1f",egenq1f,"egenq1f[5]/D");
        ExTree->Branch("ptgenq2l",ptgenq2l,"ptgenq2l[5]/D");
        ExTree->Branch("etagenq2l",etagenq2l,"etagenq2l[5]/D");
        ExTree->Branch("phigenq2l",phigenq2l,"phigenq2l[5]/D");
        ExTree->Branch("egenq2l",egenq2l,"egenq2l[5]/D");
        ExTree->Branch("ptgenq2f",ptgenq2f,"ptgenq2f[5]/D");
        ExTree->Branch("etagenq2f",etagenq2f,"etagenq2f[5]/D");
        ExTree->Branch("phigenq2f",phigenq2f,"phigenq2f[5]/D");
        ExTree->Branch("egenq2f",egenq2f,"egenq2f[5]/D");
        ExTree->Branch("ptgenq3l",ptgenq3l,"ptgenq3l[5]/D");
        ExTree->Branch("etagenq3l",etagenq3l,"etagenq3l[5]/D");
        ExTree->Branch("phigenq3l",phigenq3l,"phigenq3l[5]/D");
        ExTree->Branch("egenq3l",egenq3l,"egenq3l[5]/D");
        ExTree->Branch("ptgenq3f",ptgenq3f,"ptgenq3f[5]/D");
        ExTree->Branch("etagenq3f",etagenq3f,"etagenq3f[5]/D");
        ExTree->Branch("phigenq3f",phigenq3f,"phigenq3f[5]/D");
        ExTree->Branch("egenq3f",egenq3f,"egenq3f[5]/D");
        ExTree->Branch("ptgenq4l",ptgenq4l,"ptgenq4l[5]/D");
        ExTree->Branch("etagenq4l",etagenq4l,"etagenq4l[5]/D");
        ExTree->Branch("phigenq4l",phigenq4l,"phigenq4l[5]/D");
        ExTree->Branch("egenq4l",egenq4l,"egenq4l[5]/D");
        ExTree->Branch("ptgenq4f",ptgenq4f,"ptgenq4f[5]/D");
        ExTree->Branch("etagenq4f",etagenq4f,"etagenq4f[5]/D");
        ExTree->Branch("phigenq4f",phigenq4f,"phigenq4f[5]/D");
        ExTree->Branch("egenq4f",egenq4f,"egenq4f[5]/D");
        ExTree->Branch("ptgenq5l",ptgenq5l,"ptgenq5l[5]/D");
        ExTree->Branch("etagenq5l",etagenq5l,"etagenq5l[5]/D");
        ExTree->Branch("phigenq5l",phigenq5l,"phigenq5l[5]/D");
        ExTree->Branch("egenq5l",egenq5l,"egenq5l[5]/D");
        ExTree->Branch("ptgenq5f",ptgenq5f,"ptgenq5f[5]/D");
        ExTree->Branch("etagenq5f",etagenq5f,"etagenq5f[5]/D");
        ExTree->Branch("phigenq5f",phigenq5f,"phigenq5f[5]/D");
        ExTree->Branch("egenq5f",egenq5f,"egenq5f[5]/D");
        ExTree->Branch("mothergenq1f",mothergenq1f,"mothergenq1f[5]/D");
        ExTree->Branch("mothergenq2f",mothergenq2f,"mothergenq2f[5]/D");
        ExTree->Branch("mothergenq3f",mothergenq3f,"mothergenq3f[5]/D");
        ExTree->Branch("mothergenq4f",mothergenq4f,"mothergenq4f[5]/D");
        ExTree->Branch("mothergenq5f",mothergenq5f,"mothergenq5f[5]/D");
        ExTree->Branch("gent_w_tag",&gent_w_tag,"gent_w_tag/D");
        ExTree->Branch("genantit_w_tag",&genantit_w_tag,"genantit_w_tag/D");
        ExTree->Branch("mothergengf",mothergengf,"mothergengf[10]/D");
        ExTree->Branch("mmothergengf",mmothergengf,"mmothergengf[10]/D");
        ExTree->Branch("mmothergenq1f",mmothergenq1f,"mmothergenq1f[5]/D");
        ExTree->Branch("mmothergenq2f",mmothergenq2f,"mmothergenq2f[5]/D");
        ExTree->Branch("mmothergenq3f",mmothergenq3f,"mmothergenq3f[5]/D");
        ExTree->Branch("mmothergenq4f",mmothergenq4f,"mmothergenq4f[5]/D");
        ExTree->Branch("mmothergenq5f",mmothergenq5f,"mmothergenq5f[5]/D");

        ExTree->Branch("genw_q1_pt",genw_q1_pt,"genw_q1_pt[5]/D");
        ExTree->Branch("genw_q1_eta",genw_q1_eta,"genw_q1_eta[5]/D");
        ExTree->Branch("genw_q1_phi",genw_q1_phi,"genw_q1_phi[5]/D");
        ExTree->Branch("genw_q1_e",genw_q1_e,"genw_q1_e[5]/D");
        ExTree->Branch("genw_q1_pdg",genw_q1_pdg,"genw_q1_pdg[5]/D");
        ExTree->Branch("genw_q2_pt",genw_q2_pt,"genw_q2_pt[5]/D");
        ExTree->Branch("genw_q2_eta",genw_q2_eta,"genw_q2_eta[5]/D");
        ExTree->Branch("genw_q2_phi",genw_q2_phi,"genw_q2_phi[5]/D");
        ExTree->Branch("genw_q2_e",genw_q2_e,"genw_q2_e[5]/D");
        ExTree->Branch("genw_q2_pdg",genw_q2_pdg,"genw_q2_pdg[5]/D");*/

        
        ExTree->Branch("matchingt",matchingt,"matchingt[3]/I");
        ExTree->Branch("matchingtg",matchingtg,"matchingtg[9][3]/I");
        ExTree->Branch("matchingtq",matchingtq,"matchingtq[9][3]/I");
        ExTree->Branch("matchingw",matchingw,"matchingw[3]/I");
        ExTree->Branch("matchingg",matchingg,"matchingg[3]/I");
        ExTree->Branch("matchingz",matchingz,"matchingz[3]/I");
        ExTree->Branch("matchingq",matchingq,"matchingq[3]/I");
        ExTree->Branch("unmatching",unmatching,"unmatching[3]/I");
        ExTree->Branch("matchingr",matchingr,"matchingr[3]/I");
        ExTree->Branch("matchnewr",matchnewr,"matchnewr[4][3]/I");
        ExTree->Branch("matchnewt",matchnewt,"matchnewt[4][3]/I");
        ExTree->Branch("matchneww",matchneww,"matchneww[4][3]/I");
        //ExTree->Branch("nqqcr",nqqcr,"nqqcr[4]/I");
        //ExTree->Branch("nqqct",nqqct,"nqqct[4]/I");
        //ExTree->Branch("nqqcw",nqqcw,"nqqcw[4]/I");

    }
    /*ExTree->Branch("vak8sj11",&vak8sj11);
    ExTree->Branch("vak8sj12",&vak8sj12);
    ExTree->Branch("vak8sj13",&vak8sj13);
    ExTree->Branch("vak8sj14",&vak8sj14);
    ExTree->Branch("vak8sj15",&vak8sj15);
    ExTree->Branch("vak8sj21",&vak8sj21);
    ExTree->Branch("vak8sj22",&vak8sj22);
    ExTree->Branch("vak8sj23",&vak8sj23);
    ExTree->Branch("vak8sj24",&vak8sj24);
    ExTree->Branch("vak8sj25",&vak8sj25);*/
    //ExTree->Branch("vpuppi_softdropj1",&vpuppi_softdropj1);
    //ExTree->Branch("vpuppi_softdropj2",&vpuppi_softdropj2);
    ExTree->Branch("ak4expt1",&ak4expt1,"ak4expt1/D");
    ExTree->Branch("ak4exphi1",&ak4exphi1,"ak4exphi1/D");
    ExTree->Branch("ak4exeta1",&ak4exeta1,"ak4exeta1/D");
    ExTree->Branch("ak4exe1",&ak4exe1,"ak4exe1/D");
    ExTree->Branch("ak4expt2",&ak4expt2,"ak4expt2/D");
    ExTree->Branch("ak4exphi2",&ak4exphi2,"ak4exphi2/D");
    ExTree->Branch("ak4exeta2",&ak4exeta2,"ak4exeta2/D");
    ExTree->Branch("ak4exe2",&ak4exe2,"ak4exe2/D");

    ExTree->Branch("deltaRjet1ak4p_dR0p8",&deltaRjet1ak4p_dR0p8,"deltaRjet1ak4p_dR0p8/D");
    ExTree->Branch("deltaRjet2ak4p_dR0p8",&deltaRjet2ak4p_dR0p8,"deltaRjet2ak4p_dR0p8/D");
    ExTree->Branch("deltaRjet1ak4p_dR0p4",&deltaRjet1ak4p_dR0p4,"deltaRjet1ak4p_dR0p4/D");
    ExTree->Branch("deltaRak4mu_dR0p8",&deltaRak4mu_dR0p8,"deltaRak4mu_dR0p8/D");
    ExTree->Branch("deltaRak4mu",&deltaRak4mu,"deltaRak4mu/D");
    ExTree->Branch("deltaphijet1ak4p_dR0p8",&deltaphijet1ak4p_dR0p8,"deltaphijet1ak4p_dR0p8/D");
    ExTree->Branch("deltaphijet1ak4",&deltaphijet1ak4,"deltaphijet1ak4/D");
    ExTree->Branch("ak4inpt1_dR0p8",&ak4inpt1_dR0p8,"ak4inpt1_dR0p8/D");

    ExTree->Branch("deltaRjet1inak4p_dR0p8",&deltaRjet1inak4p_dR0p8,"deltaRjet1inak4p_dR0p8/D");
    ExTree->Branch("deltaphijet1inak4p_dR0p8",&deltaphijet1inak4p_dR0p8,"deltaphijet1inak4p_dR0p8/D");
    ExTree->Branch("deltaRinak4mu_dR0p8",&deltaRinak4mu_dR0p8,"deltaRinak4mu_dR0p8/D");
    ExTree->Branch("deltaPhiinak4mu_dR0p8",&deltaPhiinak4mu_dR0p8,"deltaPhiinak4mu_dR0p8/D");

    ExTree->Branch("Nj4ex_dR0p5",&Nj4ex_dR0p5,"Nj4ex_dR0p5/I");
    ExTree->Branch("Nj4in_dR0p5",&Nj4in_dR0p5,"Nj4in_dR0p5/I");
    ExTree->Branch("ak4expt1_dR0p5",&ak4expt1_dR0p5,"ak4expt1_dR0p5/D");
    ExTree->Branch("ak4exphi1_dR0p5",&ak4exphi1_dR0p5,"ak4exphi1_dR0p5/D");
    ExTree->Branch("ak4exeta1_dR0p5",&ak4exeta1_dR0p5,"ak4exeta1_dR0p5/D");
    ExTree->Branch("ak4exe1_dR0p5",&ak4exe1_dR0p5,"ak4exe1_dR0p5/D");
    ExTree->Branch("ak4expt2_dR0p5",&ak4expt2_dR0p5,"ak4expt2_dR0p5/D");
    ExTree->Branch("ak4exphi2_dR0p5",&ak4exphi2_dR0p5,"ak4exphi2_dR0p5/D");
    ExTree->Branch("ak4exeta2_dR0p5",&ak4exeta2_dR0p5,"ak4exeta2_dR0p5/D");
    ExTree->Branch("ak4exe2_dR0p5",&ak4exe2_dR0p5,"ak4exe2_dR0p5/D");

    ExTree->Branch("Nj4ex_dR0p8",&Nj4ex_dR0p8,"Nj4ex_dR0p8/I");
    ExTree->Branch("Nj4in_dR0p8",&Nj4in_dR0p8,"Nj4in_dR0p8/I");
    ExTree->Branch("ak4expt1_dR0p8",&ak4expt1_dR0p8,"ak4expt1_dR0p8/D");
    ExTree->Branch("ak4exphi1_dR0p8",&ak4exphi1_dR0p8,"ak4exphi1_dR0p8/D");
    ExTree->Branch("ak4exeta1_dR0p8",&ak4exeta1_dR0p8,"ak4exeta1_dR0p8/D");
    ExTree->Branch("ak4exe1_dR0p8",&ak4exe1_dR0p8,"ak4exe1_dR0p8/D");
    ExTree->Branch("ak4expt2_dR0p8",&ak4expt2_dR0p8,"ak4expt2_dR0p8/D");
    ExTree->Branch("ak4exphi2_dR0p8",&ak4exphi2_dR0p8,"ak4exphi2_dR0p8/D");
    ExTree->Branch("ak4exeta2_dR0p8",&ak4exeta2_dR0p8,"ak4exeta2_dR0p8/D");
    ExTree->Branch("ak4exe2_dR0p8",&ak4exe2_dR0p8,"ak4exe2_dR0p8/D");
    ExTree->Branch("ak4expt1_uncorr_dR0p8",&ak4expt1_uncorr_dR0p8,"ak4expt1_uncorr_dR0p8/D");

    ExTree->Branch("Nj4ex_dR1p0",&Nj4ex_dR1p0,"Nj4ex_dR1p0/I");
    ExTree->Branch("Nj4in_dR1p0",&Nj4in_dR1p0,"Nj4in_dR1p0/I");
    ExTree->Branch("ak4expt1_dR1p0",&ak4expt1_dR1p0,"ak4expt1_dR1p0/D");
    ExTree->Branch("ak4exphi1_dR1p0",&ak4exphi1_dR1p0,"ak4exphi1_dR1p0/D");
    ExTree->Branch("ak4exeta1_dR1p0",&ak4exeta1_dR1p0,"ak4exeta1_dR1p0/D");
    ExTree->Branch("ak4exe1_dR1p0",&ak4exe1_dR1p0,"ak4exe1_dR1p0/D");
    ExTree->Branch("ak4expt2_dR1p0",&ak4expt2_dR1p0,"ak4expt2_dR1p0/D");
    ExTree->Branch("ak4exphi2_dR1p0",&ak4exphi2_dR1p0,"ak4exphi2_dR1p0/D");
    ExTree->Branch("ak4exeta2_dR1p0",&ak4exeta2_dR1p0,"ak4exeta2_dR1p0/D");
    ExTree->Branch("ak4exe2_dR1p0",&ak4exe2_dR1p0,"ak4exe2_dR1p0/D");

    ExTree->Branch("Nj4ex_dR1p2",&Nj4ex_dR1p2,"Nj4ex_dR1p2/I");
    ExTree->Branch("Nj4in_dR1p2",&Nj4in_dR1p2,"Nj4in_dR1p2/I");
    ExTree->Branch("ak4expt1_dR1p2",&ak4expt1_dR1p2,"ak4expt1_dR1p2/D");
    ExTree->Branch("ak4exphi1_dR1p2",&ak4exphi1_dR1p2,"ak4exphi1_dR1p2/D");
    ExTree->Branch("ak4exeta1_dR1p2",&ak4exeta1_dR1p2,"ak4exeta1_dR1p2/D");
    ExTree->Branch("ak4exe1_dR1p2",&ak4exe1_dR1p2,"ak4exe1_dR1p2/D");
    ExTree->Branch("ak4expt2_dR1p2",&ak4expt2_dR1p2,"ak4expt2_dR1p2/D");
    ExTree->Branch("ak4exphi2_dR1p2",&ak4exphi2_dR1p2,"ak4exphi2_dR1p2/D");
    ExTree->Branch("ak4exeta2_dR1p2",&ak4exeta2_dR1p2,"ak4exeta2_dR1p2/D");
    ExTree->Branch("ak4exe2_dR1p2",&ak4exe2_dR1p2,"ak4exe2_dR1p2/D");

  
    ExTree->Branch("jetAK8puppi_mass1",jetAK8puppi_mass1,"jetAK8puppi_mass1[4]/D");
    ExTree->Branch("jetAK8puppi_pt1",jetAK8puppi_pt1,"jetAK8puppi_pt1[4]/D");
    ExTree->Branch("jetAK8puppi_eta1",jetAK8puppi_eta1,"jetAK8puppi_eta1[4]/D");

    ExTree->Branch("ptGenVlep",&ptGenVlep,"ptGenVlep/D");
    ExTree->Branch("etaGenVlep",&etaGenVlep,"etaGenVlep/D");
    ExTree->Branch("phiGenVlep",&phiGenVlep,"phiGenVlep/D");
    ExTree->Branch("massGenVlep",&massGenVlep,"massGenVlep/D");

    ExTree->Branch("ptGenV_2",&ptGenV_2,"ptGenV_2/D");
    ExTree->Branch("etaGenV_2",&etaGenV_2,"etaGenV_2/D");
    ExTree->Branch("phiGenV_2",&phiGenV_2,"phiGenV_2/D");
    ExTree->Branch("massGenV_2",&massGenV_2,"massGenV_2/D");
    ExTree->Branch("ptGenV_3",&ptGenV_3,"ptGenV_3/D");
    ExTree->Branch("etaGenV_3",&etaGenV_3,"etaGenV_3/D");
    ExTree->Branch("phiGenV_3",&phiGenV_3,"phiGenV_3/D");
    ExTree->Branch("massGenV_3",&massGenV_3,"massGenV_3/D");

    ExTree->Branch("ptGenVhad",&ptGenVhad,"ptGenVhad/D");
    ExTree->Branch("etaGenVhad",&etaGenVhad,"etaGenVhad/D");
    ExTree->Branch("phiGenVhad",&phiGenVhad,"phiGenVhad/D");
    ExTree->Branch("massGenVhad",&massGenVhad,"massGenVhad/D");
    ExTree->Branch("ptGenVhad_2",&ptGenVhad_2,"ptGenVhad_2/D");
    ExTree->Branch("etaGenVhad_2",&etaGenVhad_2,"etaGenVhad_2/D");
    ExTree->Branch("phiGenVhad_2",&phiGenVhad_2,"phiGenVhad_2/D");
    ExTree->Branch("massGenVhad_2",&massGenVhad_2,"massGenVhad_2/D");
    ExTree->Branch("ptGenVhad_3",&ptGenVhad_3,"ptGenVhad_3/D");
    ExTree->Branch("etaGenVhad_3",&etaGenVhad_3,"etaGenVhad_3/D");
    ExTree->Branch("phiGenVhad_3",&phiGenVhad_3,"phiGenVhad_3/D");
    ExTree->Branch("massGenVhad_3",&massGenVhad_3,"massGenVhad_3/D");

   ExTree->Branch("jet_pt",&jet_pt,"jet_pt/D");
   ExTree->Branch("jet_eta",&jet_eta,"jet_eta/D");
   ExTree->Branch("jet_phi",&jet_phi,"jet_phi/D");
   ExTree->Branch("jet_mass_pruned",&jet_mass_pruned,"jet_mass_pruned/D");
   ExTree->Branch("jet_mass_puppi",&jet_mass_puppi,"jet_mass_puppi/D");
   ExTree->Branch("jet_mass_puppi_un",&jet_mass_puppi_un,"jet_mass_puppi_un/D");
   ExTree->Branch("jet_mass_puppi_corr",&jet_mass_puppi_corr,"jet_mass_puppi_corr/D");
    ExTree->Branch("jet_mass_puppi_casea",&jet_mass_puppi_casea,"jet_mass_puppi_casea/D");
    ExTree->Branch("jet_mass_puppi_casea_1",&jet_mass_puppi_casea_1,"jet_mass_puppi_casea_1/D");

    ExTree->Branch("jet_mass_puppi_caseb",&jet_mass_puppi_caseb,"jet_mass_puppi_caseb/D");

   ExTree->Branch("jet_pt_puppi",&jet_pt_puppi,"jet_pt_puppi/D");
   ExTree->Branch("jetAK8puppi_tau21",&jetAK8puppi_tau21,"jetAK8puppi_tau21/D");
//   ExTree->Branch("jet_mass_softdrop",&jet_mass_softdrop,"jet_mass_softdrop/D");
    ExTree->Branch("jetAK8puppi_tau31",&jetAK8puppi_tau31,"jetAK8puppi_tau31/D");
    ExTree->Branch("jetAK8puppi_tau32",&jetAK8puppi_tau32,"jetAK8puppi_tau32/D");
    ExTree->Branch("jetAK8puppi_tau43",&jetAK8puppi_tau43,"jetAK8puppi_tau43/D");
    ExTree->Branch("jetAK8puppi_tau41",&jetAK8puppi_tau41,"jetAK8puppi_tau41/D");

    ExTree->Branch("jetAK8puppi_tau31_2",&jetAK8puppi_tau31_2,"jetAK8puppi_tau31_2/D");
    ExTree->Branch("jetAK8puppi_tau32_2",&jetAK8puppi_tau32_2,"jetAK8puppi_tau32_2/D");
    ExTree->Branch("jetAK8puppi_tau43_2",&jetAK8puppi_tau43_2,"jetAK8puppi_tau43_2/D");
    ExTree->Branch("jetAK8puppi_tau41_2",&jetAK8puppi_tau41_2,"jetAK8puppi_tau41_2/D");
    
    ExTree->Branch("jetAK8puppi_tau31_3",&jetAK8puppi_tau31_3,"jetAK8puppi_tau31_3/D");
    ExTree->Branch("jetAK8puppi_tau32_3",&jetAK8puppi_tau32_3,"jetAK8puppi_tau32_3/D");
    ExTree->Branch("jetAK8puppi_tau43_3",&jetAK8puppi_tau43_3,"jetAK8puppi_tau43_3/D");
    ExTree->Branch("jetAK8puppi_tau41_3",&jetAK8puppi_tau41_3,"jetAK8puppi_tau41_3/D");

   ExTree->Branch("jetAK8puppi_eta",&jetAK8puppi_eta,"jetAK8puppi_eta/D");
   ExTree->Branch("jetAK8puppi_phi",&jetAK8puppi_phi,"jetAK8puppi_phi/D");
   ExTree->Branch("jetAK8puppi_sdcorr",&jetAK8puppi_sdcorr,"jetAK8puppi_sdcorr/D");
   ExTree->Branch("candMasspuppiJEC",&candMasspuppiJEC,"candMasspuppiJEC/D");
//   ExTree->Branch("candMasspuppicorr",&candMasspuppicorr,"candMasspuppicorr/D");


 ExTree->Branch("jet_mass_puppi_2",&jet_mass_puppi_2,"jet_mass_puppi_2/D");
   ExTree->Branch("jet_mass_puppi_un_2",&jet_mass_puppi_un_2,"jet_mass_puppi_un_2/D");
   ExTree->Branch("jet_mass_puppi_corr_2",&jet_mass_puppi_corr_2,"jet_mass_puppi_corr_2/D");
   ExTree->Branch("jet_pt_puppi_2",&jet_pt_puppi_2,"jet_pt_puppi_2/D");
   ExTree->Branch("jetAK8puppi_tau21_2",&jetAK8puppi_tau21_2,"jetAK8puppi_tau21_2/D");
   ExTree->Branch("jetAK8puppi_tau42",&jetAK8puppi_tau42,"jetAK8puppi_tau42/D");
   ExTree->Branch("jetAK8puppi_tau42_2",&jetAK8puppi_tau42_2,"jetAK8puppi_tau42_2/D");
    ExTree->Branch("jetAK8puppi_eta_2",&jetAK8puppi_eta_2,"jetAK8puppi_eta_2/D");
    ExTree->Branch("jetAK8puppi_phi_2",&jetAK8puppi_phi_2,"jetAK8puppi_phi_2/D");
    ExTree->Branch("jetAK8puppi_sdcorr_2",&jetAK8puppi_sdcorr_2,"jetAK8puppi_sdcorr_2/D");
    
    ExTree->Branch("jet_mass_puppi_3",&jet_mass_puppi_3,"jet_mass_puppi_3/D");
    ExTree->Branch("jet_mass_puppi_un_3",&jet_mass_puppi_un_3,"jet_mass_puppi_un_3/D");
    ExTree->Branch("jet_mass_puppi_corr_3",&jet_mass_puppi_corr_3,"jet_mass_puppi_corr_3/D");
    ExTree->Branch("jet_pt_puppi_3",&jet_pt_puppi_3,"jet_pt_puppi_3/D");
    ExTree->Branch("jetAK8puppi_tau21_3",&jetAK8puppi_tau21_3,"jetAK8puppi_tau21_3/D");
    ExTree->Branch("jetAK8puppi_tau42",&jetAK8puppi_tau42,"jetAK8puppi_tau42/D");
    ExTree->Branch("jetAK8puppi_tau42_3",&jetAK8puppi_tau42_3,"jetAK8puppi_tau42_3/D");
    ExTree->Branch("jetAK8puppi_eta_3",&jetAK8puppi_eta_3,"jetAK8puppi_eta_3/D");
    ExTree->Branch("jetAK8puppi_phi_3",&jetAK8puppi_phi_3,"jetAK8puppi_phi_3/D");
    ExTree->Branch("jetAK8puppi_sdcorr_3",&jetAK8puppi_sdcorr_3,"jetAK8puppi_sdcorr_3/D");

 ExTree->Branch("MassVV",&MassVV,"MassVV[3]/D");

   ExTree->Branch("W_pt",&W_pt,"W_pt/D");
   ExTree->Branch("W_eta",&W_eta,"W_eta/D");
   ExTree->Branch("W_phi",&W_phi,"W_phi/D");
   ExTree->Branch("m_JJlv",&m_JJlv,"m_JJlv/D");
    ExTree->Branch("m_Jlv",&m_Jlv,"m_Jlv/D");

//   ExTree->Branch("M_ww",&M_ww,"M_ww/D");

//   ExTree->Branch("njets",&njets,"njets/I");
   ExTree->Branch("nbtag",&nbtag,"nbtag/I");
   ExTree->Branch("num_bJet",&num_bJet,"num_bJet/I");
    ExTree->Branch("num_sj1",&num_sj1,"num_sj1/I");
    ExTree->Branch("num_sj2",&num_sj2,"num_sj2/I");



   ExTree->Branch("num_of_lepton",&num_of_lepton,"num_of_lepton/I");
   ExTree->Branch("num_of_AK8",&num_of_AK8,"num_of_AK8/I");


   ExTree->Branch("num_bJet_loose",&num_bJet_loose,"num_bJet_loose/I");
   ExTree->Branch("num_bJet_tight",&num_bJet_tight,"num_bJet_tight/I");
    ExTree->Branch("nbtag_deep",&nbtag_deep,"nbtag_deep/I");
    ExTree->Branch("num_bJet_deep",&num_bJet_deep,"num_bJet_deep/I");
    ExTree->Branch("num_bJet_loose_deep",&num_bJet_loose_deep,"num_bJet_loose_deep/I");
    ExTree->Branch("num_bJet_tight_deep",&num_bJet_tight_deep,"num_bJet_tight_deep/I");

   ExTree->Branch("deltaRAK4AK8_new",deltaRAK4AK8_new,"deltaRAK4AK8_new[8]/D");
    ExTree->Branch("status_1",&status_1,"status_1/I");
    ExTree->Branch("status_2",&status_2,"status_2/I");
    ExTree->Branch("status_3",&status_3,"status_3/I");
    ExTree->Branch("ak4jet_eta",ak4jet_eta,"ak4jet_eta[8]/D");
    ExTree->Branch("ak4jet_phi",ak4jet_phi,"ak4jet_phi[8]/D");

    ExTree->Branch("Mj_min",&Mj_min,"Mj_min/D");
    ExTree->Branch("Mj_max",&Mj_max,"Mj_max/D");
    ExTree->Branch("pt_min",&pt_min,"pt_min/D");
    ExTree->Branch("pt_max",&pt_max,"pt_max/D");
    ExTree->Branch("phi_min",&phi_min,"phi_min/D");
    ExTree->Branch("phi_max",&phi_max,"phi_max/D");
    ExTree->Branch("eta_min",&eta_min,"eta_min/D");
    ExTree->Branch("eta_max",&eta_max,"eta_max/D");
    ExTree->Branch("m_Jlvj",&m_Jlvj,"m_Jlvj/D");
    ExTree->Branch("m_JJlvj",&m_JJlvj,"m_JJlvj/D");
    ExTree->Branch("m_Jlvjj",&m_Jlvjj,"m_Jlvjj/D");
    ExTree->Branch("m_JJlvjj",&m_JJlvjj,"m_JJlvjj/D");
    
    ExTree->Branch("m_JJJlv",&m_JJJlv,"m_JJJlv/D");
    ExTree->Branch("m_JJJlvj",&m_JJJlvj,"m_JJJlvj/D");
    ExTree->Branch("m_JJJlvjj",&m_JJJlvjj,"m_JJJlvjj/D");

    ExTree->Branch("Mj_mid",&Mj_mid,"Mj_mid/D");
    ExTree->Branch("pt_mid",&pt_mid,"pt_mid/D");
    ExTree->Branch("phi_mid",&phi_mid,"phi_mid/D");
    ExTree->Branch("eta_mid",&eta_mid,"eta_mid/D");
    ExTree->Branch("t21_mid",&t21_mid,"t21_mid/D");
    ExTree->Branch("t31_mid",&t31_mid,"t31_mid/D");
    ExTree->Branch("t41_mid",&t41_mid,"t41_mid/D");
    ExTree->Branch("t32_mid",&t32_mid,"t32_mid/D");
    ExTree->Branch("t42_mid",&t42_mid,"t42_mid/D");
    ExTree->Branch("t43_mid",&t43_mid,"t43_mid/D");

    ExTree->Branch("t21_max",&t21_max,"t21_max/D");
    ExTree->Branch("t31_max",&t31_max,"t31_max/D");
    ExTree->Branch("t41_max",&t41_max,"t41_max/D");
    ExTree->Branch("t32_max",&t32_max,"t32_max/D");
    ExTree->Branch("t42_max",&t42_max,"t42_max/D");
    ExTree->Branch("t43_max",&t43_max,"t43_max/D");

    ExTree->Branch("t32_min",&t32_min,"t32_min/D");
    ExTree->Branch("t42_min",&t42_min,"t42_min/D");
    ExTree->Branch("t43_min",&t43_min,"t43_min/D");

    ExTree->Branch("t21_min",&t21_min,"t21_min/D");
    ExTree->Branch("t31_min",&t31_min,"t31_min/D");
    ExTree->Branch("t41_min",&t41_min,"t41_min/D");
    
    ExTree->Branch("jetAK8puppi_dnnTop",&jetAK8puppi_dnnTop,"jetAK8puppi_dnnTop/D");
    ExTree->Branch("jetAK8puppi_dnnW",&jetAK8puppi_dnnW,"jetAK8puppi_dnnW/D");
    ExTree->Branch("jetAK8puppi_dnnH4q",&jetAK8puppi_dnnH4q,"jetAK8puppi_dnnH4q/D");
    
    ExTree->Branch("jetAK8puppi_dnnTop_2",&jetAK8puppi_dnnTop_2,"jetAK8puppi_dnnTop_2/D");
    ExTree->Branch("jetAK8puppi_dnnW_2",&jetAK8puppi_dnnW_2,"jetAK8puppi_dnnW_2/D");
    ExTree->Branch("jetAK8puppi_dnnH4q_2",&jetAK8puppi_dnnH4q_2,"jetAK8puppi_dnnH4q_2/D");
    
    ExTree->Branch("jetAK8puppi_dnnTop_3",&jetAK8puppi_dnnTop_3,"jetAK8puppi_dnnTop_3/D");
    ExTree->Branch("jetAK8puppi_dnnW_3",&jetAK8puppi_dnnW_3,"jetAK8puppi_dnnW_3/D");
    ExTree->Branch("jetAK8puppi_dnnH4q_3",&jetAK8puppi_dnnH4q_3,"jetAK8puppi_dnnH4q_3/D");

    ExTree->Branch("jetAK8puppi_dnnDecorrTop",&jetAK8puppi_dnnDecorrTop,"jetAK8puppi_dnnDecorrTop/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrW",&jetAK8puppi_dnnDecorrW,"jetAK8puppi_dnnDecorrW/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrH4q",&jetAK8puppi_dnnDecorrH4q,"jetAK8puppi_dnnDecorrH4q/D");
    
    ExTree->Branch("jetAK8puppi_dnnDecorrTop_2",&jetAK8puppi_dnnDecorrTop_2,"jetAK8puppi_dnnDecorrTop_2/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrW_2",&jetAK8puppi_dnnDecorrW_2,"jetAK8puppi_dnnDecorrW_2/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrH4q_2",&jetAK8puppi_dnnDecorrH4q_2,"jetAK8puppi_dnnDecorrH4q_2/D");
    
    ExTree->Branch("jetAK8puppi_dnnDecorrTop_3",&jetAK8puppi_dnnDecorrTop_3,"jetAK8puppi_dnnDecorrTop_3/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrW_3",&jetAK8puppi_dnnDecorrW_3,"jetAK8puppi_dnnDecorrW_3/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrH4q_3",&jetAK8puppi_dnnDecorrH4q_3,"jetAK8puppi_dnnDecorrH4q_3/D");

    ExTree->Branch("dnnt_max",&dnnt_max,"dnnt_max/D");
    ExTree->Branch("dnnw_max",&dnnw_max,"dnnw_max/D");
    ExTree->Branch("dnnh_max",&dnnh_max,"dnnh_max/D");
    ExTree->Branch("dnndect_max",&dnndect_max,"dnndect_max/D");
    ExTree->Branch("dnndecw_max",&dnndecw_max,"dnndecw_max/D");
    ExTree->Branch("dnndech_max",&dnndech_max,"dnndech_max/D");

    ExTree->Branch("dnnt_mid",&dnnt_mid,"dnnt_mid/D");
    ExTree->Branch("dnnw_mid",&dnnw_mid,"dnnw_mid/D");
    ExTree->Branch("dnnh_mid",&dnnh_mid,"dnnh_mid/D");
    ExTree->Branch("dnndect_mid",&dnndect_mid,"dnndect_mid/D");
    ExTree->Branch("dnndecw_mid",&dnndecw_mid,"dnndecw_mid/D");
    ExTree->Branch("dnndech_mid",&dnndech_mid,"dnndech_mid/D");

    ExTree->Branch("dnnt_min",&dnnt_min,"dnnt_min/D");
    ExTree->Branch("dnnw_min",&dnnw_min,"dnnw_min/D");
    ExTree->Branch("dnnh_min",&dnnh_min,"dnnh_min/D");
    ExTree->Branch("dnndect_min",&dnndect_min,"dnndect_min/D");
    ExTree->Branch("dnndecw_min",&dnndecw_min,"dnndecw_min/D");
    ExTree->Branch("dnndech_min",&dnndech_min,"dnndech_min/D");

    /*ExTree->Branch("jetAK8puppi_dnnZ",&jetAK8puppi_dnnZ,"jetAK8puppi_dnnZ/D");
    ExTree->Branch("jetAK8puppi_dnnZbb",&jetAK8puppi_dnnZbb,"jetAK8puppi_dnnZbb/D");
    ExTree->Branch("jetAK8puppi_dnnHbb",&jetAK8puppi_dnnHbb,"jetAK8puppi_dnnHbb/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrZ",&jetAK8puppi_dnnDecorrZ,"jetAK8puppi_dnnDecorrZ/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrZbb",&jetAK8puppi_dnnDecorrZbb,"jetAK8puppi_dnnDecorrZbb/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrHbb",&jetAK8puppi_dnnDecorrHbb,"jetAK8puppi_dnnDecorrHbb/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrbb",&jetAK8puppi_dnnDecorrbb,"jetAK8puppi_dnnDecorrbb/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrbbnog",&jetAK8puppi_dnnDecorrbbnog,"jetAK8puppi_dnnDecorrbbnog/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrcc",&jetAK8puppi_dnnDecorrcc,"jetAK8puppi_dnnDecorrcc/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrccnog",&jetAK8puppi_dnnDecorrccnog,"jetAK8puppi_dnnDecorrccnog/D");
    ExTree->Branch("jetAK8puppi_dnnZ_2",&jetAK8puppi_dnnZ_2,"jetAK8puppi_dnnZ_2/D");
    ExTree->Branch("jetAK8puppi_dnnZbb_2",&jetAK8puppi_dnnZbb_2,"jetAK8puppi_dnnZbb_2/D");
    ExTree->Branch("jetAK8puppi_dnnHbb_2",&jetAK8puppi_dnnHbb_2,"jetAK8puppi_dnnHbb_2/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrZ_2",&jetAK8puppi_dnnDecorrZ_2,"jetAK8puppi_dnnDecorrZ_2/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrZbb_2",&jetAK8puppi_dnnDecorrZbb_2,"jetAK8puppi_dnnDecorrZbb_2/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrHbb_2",&jetAK8puppi_dnnDecorrHbb_2,"jetAK8puppi_dnnDecorrHbb_2/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrbb_2",&jetAK8puppi_dnnDecorrbb_2,"jetAK8puppi_dnnDecorrbb_2/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrbbnog_2",&jetAK8puppi_dnnDecorrbbnog_2,"jetAK8puppi_dnnDecorrbbnog_2/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrcc_2",&jetAK8puppi_dnnDecorrcc_2,"jetAK8puppi_dnnDecorrcc_2/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrccnog_2",&jetAK8puppi_dnnDecorrccnog_2,"jetAK8puppi_dnnDecorrccnog_2/D");
    ExTree->Branch("jetAK8puppi_dnnZ_3",&jetAK8puppi_dnnZ_3,"jetAK8puppi_dnnZ_3/D");
    ExTree->Branch("jetAK8puppi_dnnZbb_3",&jetAK8puppi_dnnZbb_3,"jetAK8puppi_dnnZbb_3/D");
    ExTree->Branch("jetAK8puppi_dnnHbb_3",&jetAK8puppi_dnnHbb_3,"jetAK8puppi_dnnHbb_3/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrZ_3",&jetAK8puppi_dnnDecorrZ_3,"jetAK8puppi_dnnDecorrZ_3/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrZbb_3",&jetAK8puppi_dnnDecorrZbb_3,"jetAK8puppi_dnnDecorrZbb_3/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrHbb_3",&jetAK8puppi_dnnDecorrHbb_3,"jetAK8puppi_dnnDecorrHbb_3/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrbb_3",&jetAK8puppi_dnnDecorrbb_3,"jetAK8puppi_dnnDecorrbb_3/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrbbnog_3",&jetAK8puppi_dnnDecorrbbnog_3,"jetAK8puppi_dnnDecorrbbnog_3/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrcc_3",&jetAK8puppi_dnnDecorrcc_3,"jetAK8puppi_dnnDecorrcc_3/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrccnog_3",&jetAK8puppi_dnnDecorrccnog_3,"jetAK8puppi_dnnDecorrccnog_3/D");
    */
    ExTree->Branch("jetAK8puppi_dnnDecorrqcd",&jetAK8puppi_dnnDecorrqcd,"jetAK8puppi_dnnDecorrqcd/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrtop",&jetAK8puppi_dnnDecorrtop,"jetAK8puppi_dnnDecorrtop/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrz",&jetAK8puppi_dnnDecorrz,"jetAK8puppi_dnnDecorrz/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrw",&jetAK8puppi_dnnDecorrw,"jetAK8puppi_dnnDecorrw/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrzbb",&jetAK8puppi_dnnDecorrzbb,"jetAK8puppi_dnnDecorrzbb/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrhbb",&jetAK8puppi_dnnDecorrhbb,"jetAK8puppi_dnnDecorrhbb/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrh4q",&jetAK8puppi_dnnDecorrh4q,"jetAK8puppi_dnnDecorrh4q/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrqcd_2",&jetAK8puppi_dnnDecorrqcd_2,"jetAK8puppi_dnnDecorrqcd_2/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrtop_2",&jetAK8puppi_dnnDecorrtop_2,"jetAK8puppi_dnnDecorrtop_2/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrz_2",&jetAK8puppi_dnnDecorrz_2,"jetAK8puppi_dnnDecorrz_2/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrw_2",&jetAK8puppi_dnnDecorrw_2,"jetAK8puppi_dnnDecorrw_2/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrzbb_2",&jetAK8puppi_dnnDecorrzbb_2,"jetAK8puppi_dnnDecorrzbb_2/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrhbb_2",&jetAK8puppi_dnnDecorrhbb_2,"jetAK8puppi_dnnDecorrhbb_2/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrh4q_2",&jetAK8puppi_dnnDecorrh4q_2,"jetAK8puppi_dnnDecorrh4q_2/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrqcd_3",&jetAK8puppi_dnnDecorrqcd_3,"jetAK8puppi_dnnDecorrqcd_3/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrtop_3",&jetAK8puppi_dnnDecorrtop_3,"jetAK8puppi_dnnDecorrtop_3/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrz_3",&jetAK8puppi_dnnDecorrz_3,"jetAK8puppi_dnnDecorrz_3/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrw_3",&jetAK8puppi_dnnDecorrw_3,"jetAK8puppi_dnnDecorrw_3/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrzbb_3",&jetAK8puppi_dnnDecorrzbb_3,"jetAK8puppi_dnnDecorrzbb_3/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrhbb_3",&jetAK8puppi_dnnDecorrhbb_3,"jetAK8puppi_dnnDecorrhbb_3/D");
    ExTree->Branch("jetAK8puppi_dnnDecorrh4q_3",&jetAK8puppi_dnnDecorrh4q_3,"jetAK8puppi_dnnDecorrh4q_3/D");


    /*ExTree->Branch("dnnz_max",&dnnz_max,"dnnz_max/D");
    ExTree->Branch("dnnzbb_max",&dnnzbb_max,"dnnzbb_max/D");
    ExTree->Branch("dnnhbb_max",&dnnhbb_max,"dnnhbb_max/D");
    ExTree->Branch("dnndecz_max",&dnndecz_max,"dnndecz_max/D");
    ExTree->Branch("dnndeczbb_max",&dnndeczbb_max,"dnndeczbb_max/D");
    ExTree->Branch("dnndechbb_max",&dnndechbb_max,"dnndechbb_max/D");
    ExTree->Branch("dnndecbb_max",&dnndecbb_max,"dnndecbb_max/D");
    ExTree->Branch("dnndecbbnog_max",&dnndecbbnog_max,"dnndecbbnog_max/D");
    ExTree->Branch("dnndeccc_max",&dnndeccc_max,"dnndeccc_max/D");
    ExTree->Branch("dnndecccnog_max",&dnndecccnog_max,"dnndecccnog_max/D");
    
    ExTree->Branch("dnnz_mid",&dnnz_mid,"dnnz_mid/D");
    ExTree->Branch("dnnzbb_mid",&dnnzbb_mid,"dnnzbb_mid/D");
    ExTree->Branch("dnnhbb_mid",&dnnhbb_mid,"dnnhbb_mid/D");
    ExTree->Branch("dnndecz_mid",&dnndecz_mid,"dnndecz_mid/D");
    ExTree->Branch("dnndeczbb_mid",&dnndeczbb_mid,"dnndeczbb_mid/D");
    ExTree->Branch("dnndechbb_mid",&dnndechbb_mid,"dnndechbb_mid/D");
    ExTree->Branch("dnndecbb_mid",&dnndecbb_mid,"dnndecbb_mid/D");
    ExTree->Branch("dnndecbbnog_mid",&dnndecbbnog_mid,"dnndecbbnog_mid/D");
    ExTree->Branch("dnndeccc_mid",&dnndeccc_mid,"dnndeccc_mid/D");
    ExTree->Branch("dnndecccnog_mid",&dnndecccnog_mid,"dnndecccnog_mid/D");
    
    ExTree->Branch("dnnz_min",&dnnz_min,"dnnz_min/D");
    ExTree->Branch("dnnzbb_min",&dnnzbb_min,"dnnzbb_min/D");
    ExTree->Branch("dnnhbb_min",&dnnhbb_min,"dnnhbb_min/D");
    ExTree->Branch("dnndecz_min",&dnndecz_min,"dnndecz_min/D");
    ExTree->Branch("dnndeczbb_min",&dnndeczbb_min,"dnndeczbb_min/D");
    ExTree->Branch("dnndechbb_min",&dnndechbb_min,"dnndechbb_min/D");
    ExTree->Branch("dnndecbb_min",&dnndecbb_min,"dnndecbb_min/D");
    ExTree->Branch("dnndecbbnog_min",&dnndecbbnog_min,"dnndecbbnog_min/D");
    ExTree->Branch("dnndeccc_min",&dnndeccc_min,"dnndeccc_min/D");
    ExTree->Branch("dnndecccnog_min",&dnndecccnog_min,"dnndecccnog_min/D");*/

    ExTree->Branch("dnndecrawqcd_max",&dnndecrawqcd_max,"dnndecrawqcd_max/D");
    ExTree->Branch("dnndecrawt_max",&dnndecrawt_max,"dnndecrawt_max/D");
    ExTree->Branch("dnndecraww_max",&dnndecraww_max,"dnndecraww_max/D");
    ExTree->Branch("dnndecrawz_max",&dnndecrawz_max,"dnndecrawz_max/D");
    ExTree->Branch("dnndecrawzbb_max",&dnndecrawzbb_max,"dnndecrawzbb_max/D");
    ExTree->Branch("dnndecrawhbb_max",&dnndecrawhbb_max,"dnndecrawhbb_max/D");
    ExTree->Branch("dnndecrawh4q_max",&dnndecrawh4q_max,"dnndecrawh4q_max/D");
    ExTree->Branch("dnndecrawqcd_mid",&dnndecrawqcd_mid,"dnndecrawqcd_mid/D");
    ExTree->Branch("dnndecrawt_mid",&dnndecrawt_mid,"dnndecrawt_mid/D");
    ExTree->Branch("dnndecraww_mid",&dnndecraww_mid,"dnndecraww_mid/D");
    ExTree->Branch("dnndecrawz_mid",&dnndecrawz_mid,"dnndecrawz_mid/D");
    ExTree->Branch("dnndecrawzbb_mid",&dnndecrawzbb_mid,"dnndecrawzbb_mid/D");
    ExTree->Branch("dnndecrawhbb_mid",&dnndecrawhbb_mid,"dnndecrawhbb_mid/D");
    ExTree->Branch("dnndecrawh4q_mid",&dnndecrawh4q_mid,"dnndecrawh4q_mid/D");
    ExTree->Branch("dnndecrawqcd_min",&dnndecrawqcd_min,"dnndecrawqcd_min/D");
    ExTree->Branch("dnndecrawt_min",&dnndecrawt_min,"dnndecrawt_min/D");
    ExTree->Branch("dnndecraww_min",&dnndecraww_min,"dnndecraww_min/D");
    ExTree->Branch("dnndecrawz_min",&dnndecrawz_min,"dnndecrawz_min/D");
    ExTree->Branch("dnndecrawzbb_min",&dnndecrawzbb_min,"dnndecrawzbb_min/D");
    ExTree->Branch("dnndecrawhbb_min",&dnndecrawhbb_min,"dnndecrawhbb_min/D");
    ExTree->Branch("dnndecrawh4q_min",&dnndecrawh4q_min,"dnndecrawh4q_min/D");
    
    ExTree->Branch("PTw2divPTw1",&PTw2divPTw1,"PTw2divPTw1/D");
    ExTree->Branch("PTw3divPTw1",&PTw3divPTw1,"PTw3divPTw1/D");
    ExTree->Branch("PTw3divPTw2",&PTw3divPTw2,"PTw3divPTw2/D");
    
    ExTree->Branch("HLT_Mu2",&HLT_Mu2,"HLT_Mu2/I");
    ExTree->Branch("HLT_Mu3",&HLT_Mu3,"HLT_Mu3/I");
    ExTree->Branch("HLT_Mu12",&HLT_Mu12,"HLT_Mu12/I");
    ExTree->Branch("HLT_Mu2",&HLT_Mu2,"HLT_Mu2/I");
    ExTree->Branch("HLT_Ele3",&HLT_Ele3,"HLT_Ele3/I");
    ExTree->Branch("HLT_Ele4",&HLT_Ele4,"HLT_Ele4/I");
    ExTree->Branch("HLT_Ele8",&HLT_Ele8,"HLT_Ele8/I");

//   ExTree->Branch("jet2_pt",&fjet2_pt,"jet2_pt/D");
//   ExTree->Branch("jet2_btag",&fjet2_btag,"jet2_btag/D");
//   ExTree->Branch("jet3_pt",&fjet3_pt,"jet3_pt/D");
//   ExTree->Branch("jet3_btag",&fjet3_btag,"jet3_btag/D");
//JEC

//   ExTree->Branch("corr_AK8",&corr_AK8,"corr_AK8/D");
//   ExTree->Branch("corr",&corr,"corr/D");
   ExTree->Branch("METraw_et",&METraw_et,"METraw_et/D");
   ExTree->Branch("METraw_phi",&METraw_phi,"METraw_phi/D");
   ExTree->Branch("METraw_sumEt",&METraw_sumEt,"METraw_sumEt/D");
   ExTree->Branch("MET_et",&MET_et,"MET_et/D");
   ExTree->Branch("met",&met,"met/D");
   ExTree->Branch("metPhi",&metPhi,"metPhi/D");
   ExTree->Branch("MET_phi",&MET_phi,"MET_phi/D");
   ExTree->Branch("MET_sumEt",&MET_sumEt,"MET_sumEt/D");
    ExTree->Branch("deltaRjet1jet2",&deltaRjet1jet2,"deltaRjet1jet2/D");
    ExTree->Branch("deltaetajet1jet2",&deltaetajet1jet2,"deltaetajet1jet2/D");
    ExTree->Branch("deltaetajet1lep",&deltaetajet1lep,"deltaetajet1lep/D");
    ExTree->Branch("deltaetajet2lep",&deltaetajet2lep,"deltaetajet2lep/D");
    ExTree->Branch("ak4jet_pt",&ak4jet_pt,"ak4jet_pt[8]/D");
    ExTree->Branch("ak4jet_pt_uncorr",&ak4jet_pt_uncorr,"ak4jet_pt_uncorr[8]/D");

    ExTree->Branch("delPhilepmet",&delPhilepmet,"delPhilepmet/D");
    ExTree->Branch("delPhilepmet_m",&delPhilepmet_m,"delPhilepmet_m/D");

    ExTree->Branch("delPhiWjet1",&delPhiWjet1,"delPhiWjet1/D");
    ExTree->Branch("delPhiWjet2",&delPhiWjet2,"delPhiWjet2/D");
    ExTree->Branch("delPhiWjetclose",&delPhiWjetclose,"delPhiWjetclose/D");
    ExTree->Branch("delPhiWjetfar",&delPhiWjetfar,"delPhiWjetfar/D");
    
    ExTree->Branch("dR_max_WW",&dR_max_WW,"dR_max_WW/D");
    ExTree->Branch("dR_min_WW",&dR_min_WW,"dR_min_WW/D");
    ExTree->Branch("dPhi_max_WW",&dPhi_max_WW,"dPhi_max_WW/D");
    ExTree->Branch("dPhi_min_WW",&dPhi_min_WW,"dPhi_min_WW/D");
    ExTree->Branch("dEta_max_WW",&dEta_max_WW,"dEta_max_WW/D");
    ExTree->Branch("dEta_min_WW",&dEta_min_WW,"dEta_min_WW/D");

    ExTree->Branch("dR_j1W",&dR_j1W,"dR_j1W/D");
    ExTree->Branch("dR_j2W",&dR_j2W,"dR_j2W/D");
    ExTree->Branch("dR_j3W",&dR_j3W,"dR_j3W/D");
    ExTree->Branch("dR_j1j2",&dR_j1j2,"dR_j1j2/D");
    ExTree->Branch("dR_j1j3",&dR_j1j3,"dR_j1j3/D");
    ExTree->Branch("dR_j2j3",&dR_j2j3,"dR_j2j3/D");
    ExTree->Branch("dPhi_j1W",&dPhi_j1W,"dPhi_j1W/D");
    ExTree->Branch("dPhi_j2W",&dPhi_j2W,"dPhi_j2W/D");
    ExTree->Branch("dPhi_j3W",&dPhi_j3W,"dPhi_j3W/D");
    ExTree->Branch("dPhi_j1j2",&dPhi_j1j2,"dPhi_j1j2/D");
    ExTree->Branch("dPhi_j1j3",&dPhi_j1j3,"dPhi_j1j3/D");
    ExTree->Branch("dPhi_j2j3",&dPhi_j2j3,"dPhi_j2j3/D");
    ExTree->Branch("dEta_j1W",&dEta_j1W,"dEta_j1W/D");
    ExTree->Branch("dEta_j2W",&dEta_j2W,"dEta_j2W/D");
    ExTree->Branch("dEta_j3W",&dEta_j3W,"dEta_j3W/D");
    ExTree->Branch("dEta_j1j2",&dEta_j1j2,"dEta_j1j2/D");
    ExTree->Branch("dEta_j1j3",&dEta_j1j3,"dEta_j1j3/D");
    ExTree->Branch("dEta_j2j3",&dEta_j2j3,"dEta_j2j3/D");

    ExTree->Branch("dR_max_WW_mass",&dR_max_WW_mass,"dR_max_WW_mass/D");
    ExTree->Branch("dR_min_WW_mass",&dR_min_WW_mass,"dR_min_WW_mass/D");
    ExTree->Branch("dPhi_max_WW_mass",&dPhi_max_WW_mass,"dPhi_max_WW_mass/D");
    ExTree->Branch("dPhi_min_WW_mass",&dPhi_min_WW_mass,"dPhi_min_WW_mass/D");
    ExTree->Branch("dEta_max_WW_mass",&dEta_max_WW_mass,"dEta_max_WW_mass/D");
    ExTree->Branch("dEta_min_WW_mass",&dEta_min_WW_mass,"dEta_min_WW_mass/D");
    
    ExTree->Branch("dR_jmaxW",&dR_jmaxW,"dR_jmaxW/D");
    ExTree->Branch("dR_jmidW",&dR_jmidW,"dR_jmidW/D");
    ExTree->Branch("dR_jminW",&dR_jminW,"dR_jminW/D");
    ExTree->Branch("dR_jmaxjmid",&dR_jmaxjmid,"dR_jmaxjmid/D");
    ExTree->Branch("dR_jmaxjmin",&dR_jmaxjmin,"dR_jmaxjmin/D");
    ExTree->Branch("dR_jmidjmin",&dR_jmidjmin,"dR_jmidjmin/D");
    ExTree->Branch("dPhi_jmaxW",&dPhi_jmaxW,"dPhi_jmaxW/D");
    ExTree->Branch("dPhi_jmidW",&dPhi_jmidW,"dPhi_jmidW/D");
    ExTree->Branch("dPhi_jminW",&dPhi_jminW,"dPhi_jminW/D");
    ExTree->Branch("dPhi_jmaxjmid",&dPhi_jmaxjmid,"dPhi_jmaxjmid/D");
    ExTree->Branch("dPhi_jmaxjmin",&dPhi_jmaxjmin,"dPhi_jmaxjmin/D");
    ExTree->Branch("dPhi_jmidjmin",&dPhi_jmidjmin,"dPhi_jmidjmin/D");
    ExTree->Branch("dEta_jmaxW",&dEta_jmaxW,"dEta_jmaxW/D");
    ExTree->Branch("dEta_jmidW",&dEta_jmidW,"dEta_jmidW/D");
    ExTree->Branch("dEta_jminW",&dEta_jminW,"dEta_jminW/D");
    ExTree->Branch("dEta_jmaxjmid",&dEta_jmaxjmid,"dEta_jmaxjmid/D");
    ExTree->Branch("dEta_jmaxjmin",&dEta_jmaxjmin,"dEta_jmaxjmin/D");
    ExTree->Branch("dEta_jmidjmin",&dEta_jmidjmin,"dEta_jmidjmin/D");

   ExTree->Branch("deltaRlepjet",&deltaRlepjet,"deltaRlepjet/D");
    ExTree->Branch("deltaRjet1ak4",deltaRjet1ak4,"deltaRjet1ak4[8]/D");
    ExTree->Branch("deltaRjet1ak4_1",deltaRjet1ak4_1,"deltaRjet1ak4_1[8]/D");
    ExTree->Branch("PTw",PTw,"PTw[3]/D");

    ExTree->Branch("deltaRphiak4",deltaRphiak4,"deltaRphiak4[8]/D");

   ExTree->Branch("delPhijetmet",&delPhijetmet,"delPhijetmet/D");
   ExTree->Branch("delPhijetlep",&delPhijetlep,"delPhijetlep/D");
    ExTree->Branch("deltaRlepjet_2",&deltaRlepjet_2,"deltaRlepjet_2/D");
    ExTree->Branch("delPhijetmet_2",&delPhijetmet_2,"delPhijetmet_2/D");
    ExTree->Branch("delPhijetlep_2",&delPhijetlep_2,"delPhijetlep_2/D");

   ExTree->Branch("vbftag",&vbftag,"vbftag/I");
   ExTree->Branch("muisolation",&muisolation,"muisolation/D");
    ExTree->Branch("ST",&ST,"ST/D");
    if(IsData>0) {
        //fChain->SetBranchAddress("pweight", pweight, &b_pweight);

    fChain->SetBranchAddress("m_jlv_new", &m_jlv_new, &b_m_jlv_new);




    fChain->SetBranchAddress("m_jlv_JEC_up", &m_jlv_JEC_up, &b_m_jlv_JEC_up);
    fChain->SetBranchAddress("m_jlv_JEC_down", &m_jlv_JEC_down, &b_m_jlv_JEC_down);
    fChain->SetBranchAddress("m_jlv_JER_up", &m_jlv_JER_up, &b_m_jlv_JER_up);
    fChain->SetBranchAddress("m_jlv_JER_down", &m_jlv_JER_down, &b_m_jlv_JER_down);

        fChain->SetBranchAddress("jetAK8puppi_ptJEC_new", &jetAK8puppi_ptJEC_new, &b_jetAK8puppi_ptJEC_new);
        fChain->SetBranchAddress("jetAK8puppi_ptJEC_JEC_up", &jetAK8puppi_ptJEC_JEC_up, &b_jetAK8puppi_ptJEC_JEC_up);
        fChain->SetBranchAddress("jetAK8puppi_ptJEC_JEC_down", &jetAK8puppi_ptJEC_JEC_down, &b_jetAK8puppi_ptJEC_JEC_down);
        fChain->SetBranchAddress("jetAK8puppi_ptJEC_JER_up", &jetAK8puppi_ptJEC_JER_up, &b_jetAK8puppi_ptJEC_JER_up);
        fChain->SetBranchAddress("jetAK8puppi_ptJEC_JER_down", &jetAK8puppi_ptJEC_JER_down, &b_jetAK8puppi_ptJEC_JER_down);
        fChain->SetBranchAddress("jetAK8puppi_ptJEC_2_new", &jetAK8puppi_ptJEC_2_new, &b_jetAK8puppi_ptJEC_2_new);
        fChain->SetBranchAddress("jetAK8puppi_ptJEC_2_JEC_up", &jetAK8puppi_ptJEC_2_JEC_up, &b_jetAK8puppi_ptJEC_2_JEC_up);
        fChain->SetBranchAddress("jetAK8puppi_ptJEC_2_JEC_down", &jetAK8puppi_ptJEC_2_JEC_down, &b_jetAK8puppi_ptJEC_2_JEC_down);
        fChain->SetBranchAddress("jetAK8puppi_ptJEC_2_JER_up", &jetAK8puppi_ptJEC_2_JER_up, &b_jetAK8puppi_ptJEC_2_JER_up);
        fChain->SetBranchAddress("jetAK8puppi_ptJEC_2_JER_down", &jetAK8puppi_ptJEC_2_JER_down, &b_jetAK8puppi_ptJEC_2_JER_down);
        fChain->SetBranchAddress("jetAK8puppi_ptJEC_3_new", &jetAK8puppi_ptJEC_3_new, &b_jetAK8puppi_ptJEC_3_new);
        fChain->SetBranchAddress("jetAK8puppi_ptJEC_3_JEC_up", &jetAK8puppi_ptJEC_3_JEC_up, &b_jetAK8puppi_ptJEC_3_JEC_up);
        fChain->SetBranchAddress("jetAK8puppi_ptJEC_3_JEC_down", &jetAK8puppi_ptJEC_3_JEC_down, &b_jetAK8puppi_ptJEC_3_JEC_down);
        fChain->SetBranchAddress("jetAK8puppi_ptJEC_3_JER_up", &jetAK8puppi_ptJEC_3_JER_up, &b_jetAK8puppi_ptJEC_3_JER_up);
        fChain->SetBranchAddress("jetAK8puppi_ptJEC_3_JER_down", &jetAK8puppi_ptJEC_3_JER_down, &b_jetAK8puppi_ptJEC_3_JER_down);
        fChain->SetBranchAddress("jetAK8puppi_e", &jetAK8puppi_e, &b_jetAK8puppi_e);
        fChain->SetBranchAddress("jetAK8puppi_e_new", &jetAK8puppi_e_new, &b_jetAK8puppi_e_new);
        fChain->SetBranchAddress("jetAK8puppi_e_JEC_up", &jetAK8puppi_e_JEC_up, &b_jetAK8puppi_e_JEC_up);
        fChain->SetBranchAddress("jetAK8puppi_e_JEC_down", &jetAK8puppi_e_JEC_down, &b_jetAK8puppi_e_JEC_down);
        fChain->SetBranchAddress("jetAK8puppi_e_JER_up", &jetAK8puppi_e_JER_up, &b_jetAK8puppi_e_JER_up);
        fChain->SetBranchAddress("jetAK8puppi_e_JER_down", &jetAK8puppi_e_JER_down, &b_jetAK8puppi_e_JER_down);
        fChain->SetBranchAddress("jetAK8puppi_e_2", &jetAK8puppi_e_2, &b_jetAK8puppi_e_2);
        fChain->SetBranchAddress("jetAK8puppi_e_2_new", &jetAK8puppi_e_2_new, &b_jetAK8puppi_e_2_new);
        fChain->SetBranchAddress("jetAK8puppi_e_2_JEC_up", &jetAK8puppi_e_2_JEC_up, &b_jetAK8puppi_e_2_JEC_up);
        fChain->SetBranchAddress("jetAK8puppi_e_2_JEC_down", &jetAK8puppi_e_2_JEC_down, &b_jetAK8puppi_e_2_JEC_down);
        fChain->SetBranchAddress("jetAK8puppi_e_2_JER_up", &jetAK8puppi_e_2_JER_up, &b_jetAK8puppi_e_2_JER_up);
        fChain->SetBranchAddress("jetAK8puppi_e_2_JER_down", &jetAK8puppi_e_2_JER_down, &b_jetAK8puppi_e_2_JER_down);
        fChain->SetBranchAddress("jetAK8puppi_e_3", &jetAK8puppi_e_3, &b_jetAK8puppi_e_3);
        fChain->SetBranchAddress("jetAK8puppi_e_3_new", &jetAK8puppi_e_3_new, &b_jetAK8puppi_e_3_new);
        fChain->SetBranchAddress("jetAK8puppi_e_3_JEC_up", &jetAK8puppi_e_3_JEC_up, &b_jetAK8puppi_e_3_JEC_up);
        fChain->SetBranchAddress("jetAK8puppi_e_3_JEC_down", &jetAK8puppi_e_3_JEC_down, &b_jetAK8puppi_e_3_JEC_down);
        fChain->SetBranchAddress("jetAK8puppi_e_3_JER_up", &jetAK8puppi_e_3_JER_up, &b_jetAK8puppi_e_3_JER_up);
        fChain->SetBranchAddress("jetAK8puppi_e_3_JER_down", &jetAK8puppi_e_3_JER_down, &b_jetAK8puppi_e_3_JER_down);
        fChain->SetBranchAddress("ptVlepJEC_JEC_up", &ptVlepJEC_JEC_up, &b_ptVlepJEC_JEC_up);
        fChain->SetBranchAddress("yVlepJEC_JEC_up", &yVlepJEC_JEC_up, &b_yVlepJEC_JEC_up);
        fChain->SetBranchAddress("phiVlepJEC_JEC_up", &phiVlepJEC_JEC_up, &b_phiVlepJEC_JEC_up);
        fChain->SetBranchAddress("massVlepJEC_JEC_up", &massVlepJEC_JEC_up, &b_massVlepJEC_JEC_up);
        fChain->SetBranchAddress("mtVlepJEC_JEC_up", &mtVlepJEC_JEC_up, &b_mtVlepJEC_JEC_up);

        fChain->SetBranchAddress("ptVlepJEC_new", &ptVlepJEC_new, &b_ptVlepJEC_new);
        fChain->SetBranchAddress("yVlepJEC_new", &yVlepJEC_new, &b_yVlepJEC_new);
        fChain->SetBranchAddress("phiVlepJEC_new", &phiVlepJEC_new, &b_phiVlepJEC_new);
        fChain->SetBranchAddress("massVlepJEC_new", &massVlepJEC_new, &b_massVlepJEC_new);
        fChain->SetBranchAddress("mtVlepJEC_new", &mtVlepJEC_new, &b_mtVlepJEC_new);
        fChain->SetBranchAddress("MET_phi_new", &MET_phi_new, &b_MET_phi_new);

        fChain->SetBranchAddress("ptVlepJEC_JEC_down", &ptVlepJEC_JEC_down, &b_ptVlepJEC_JEC_down);
        fChain->SetBranchAddress("yVlepJEC_JEC_down", &yVlepJEC_JEC_down, &b_yVlepJEC_JEC_down);
        fChain->SetBranchAddress("phiVlepJEC_JEC_down", &phiVlepJEC_JEC_down, &b_phiVlepJEC_JEC_down);
        fChain->SetBranchAddress("massVlepJEC_JEC_down", &massVlepJEC_JEC_down, &b_massVlepJEC_JEC_down);
        fChain->SetBranchAddress("mtVlepJEC_JEC_down", &mtVlepJEC_JEC_down, &b_mtVlepJEC_JEC_down);
        fChain->SetBranchAddress("ptVlepJEC_JER_up", &ptVlepJEC_JER_up, &b_ptVlepJEC_JER_up);
        fChain->SetBranchAddress("yVlepJEC_JER_up", &yVlepJEC_JER_up, &b_yVlepJEC_JER_up);
        fChain->SetBranchAddress("phiVlepJEC_JER_up", &phiVlepJEC_JER_up, &b_phiVlepJEC_JER_up);
        fChain->SetBranchAddress("massVlepJEC_JER_up", &massVlepJEC_JER_up, &b_massVlepJEC_JER_up);
        fChain->SetBranchAddress("mtVlepJEC_JER_up", &mtVlepJEC_JER_up, &b_mtVlepJEC_JER_up);
        fChain->SetBranchAddress("ptVlepJEC_JER_down", &ptVlepJEC_JER_down, &b_ptVlepJEC_JER_down);
        fChain->SetBranchAddress("yVlepJEC_JER_down", &yVlepJEC_JER_down, &b_yVlepJEC_JER_down);
        fChain->SetBranchAddress("phiVlepJEC_JER_down", &phiVlepJEC_JER_down, &b_phiVlepJEC_JER_down);
        fChain->SetBranchAddress("massVlepJEC_JER_down", &massVlepJEC_JER_down, &b_massVlepJEC_JER_down);
        fChain->SetBranchAddress("mtVlepJEC_JER_down", &mtVlepJEC_JER_down, &b_mtVlepJEC_JER_down);
        fChain->SetBranchAddress("MET_et_new", &MET_et_new, &b_MET_et_new);
        fChain->SetBranchAddress("MET_et_JEC_up", &MET_et_JEC_up, &b_MET_et_JEC_up);
        fChain->SetBranchAddress("MET_et_JEC_down", &MET_et_JEC_down, &b_MET_et_JEC_down);
        fChain->SetBranchAddress("MET_et_JER_up", &MET_et_JER_up, &b_MET_et_JER_up);
        fChain->SetBranchAddress("MET_et_JER_down", &MET_et_JER_down, &b_MET_et_JER_down);
        fChain->SetBranchAddress("MET_phi_JEC_up", &MET_phi_JEC_up, &b_MET_phi_JEC_up);
        fChain->SetBranchAddress("MET_phi_JEC_down", &MET_phi_JEC_down, &b_MET_phi_JEC_down);
        fChain->SetBranchAddress("MET_phi_JER_up", &MET_phi_JER_up, &b_MET_phi_JER_up);
        fChain->SetBranchAddress("MET_phi_JER_down", &MET_phi_JER_down, &b_MET_phi_JER_down);

        
        fChain->SetBranchAddress("num_of_AK8", num_of_AK8, &num_of_AK8);
        fChain->SetBranchAddress("num_of_lepton", num_of_lepton, &b_num_of_lepton);



        fChain->SetBranchAddress("ptgenwl", ptgenwl, &b_ptgenwl);
        fChain->SetBranchAddress("etagenwl", etagenwl, &b_etagenwl);
        fChain->SetBranchAddress("phigenwl", phigenwl, &b_phigenwl);
        fChain->SetBranchAddress("massgenwl", massgenwl, &b_massgenwl);
        fChain->SetBranchAddress("taggenwl", taggenwl, &b_taggenwl);
        fChain->SetBranchAddress("pweight", pweight, &b_pweight);
        fChain->SetBranchAddress("ptgenzl", ptgenzl, &b_ptgenzl);
        fChain->SetBranchAddress("etagenzl", etagenzl, &b_etagenzl);
        fChain->SetBranchAddress("phigenzl", phigenzl, &b_phigenzl);
        fChain->SetBranchAddress("massgenzl", massgenzl, &b_massgenzl);
        fChain->SetBranchAddress("taggenzl", taggenzl, &b_taggenzl);
        fChain->SetBranchAddress("ptgengl", ptgengl, &b_ptgengl);
        fChain->SetBranchAddress("etagengl", etagengl, &b_etagengl);
        fChain->SetBranchAddress("phigengl", phigengl, &b_phigengl);
        fChain->SetBranchAddress("egengl", egengl, &b_egengl);
        fChain->SetBranchAddress("ptgenwf", ptgenwf, &b_ptgenwf);
        fChain->SetBranchAddress("etagenwf", etagenwf, &b_etagenwf);
        fChain->SetBranchAddress("phigenwf", phigenwf, &b_phigenwf);
        fChain->SetBranchAddress("massgenwf", massgenwf, &b_massgenwf);
        fChain->SetBranchAddress("ptgenzf", ptgenzf, &b_ptgenzf);
        fChain->SetBranchAddress("etagenzf", etagenzf, &b_etagenzf);
        fChain->SetBranchAddress("phigenzf", phigenzf, &b_phigenzf);
        fChain->SetBranchAddress("massgenzf", massgenzf, &b_massgenzf);
        fChain->SetBranchAddress("ptgengf", ptgengf, &b_ptgengf);
        fChain->SetBranchAddress("etagengf", etagengf, &b_etagengf);
        fChain->SetBranchAddress("phigengf", phigengf, &b_phigengf);
        fChain->SetBranchAddress("egengf", egengf, &b_egengf);

        fChain->SetBranchAddress("gent_b_pt", &gent_b_pt, &b_gent_b_pt);
        fChain->SetBranchAddress("gent_b_phi", &gent_b_phi, &b_gent_b_phi);
        fChain->SetBranchAddress("gent_b_eta", &gent_b_eta, &b_gent_b_eta);
        fChain->SetBranchAddress("gent_b_mass", &gent_b_mass, &b_gent_b_mass);
        fChain->SetBranchAddress("genantit_b_pt", &genantit_b_pt, &b_genantit_b_pt);
        fChain->SetBranchAddress("genantit_b_phi", &genantit_b_phi, &b_genantit_b_phi);
        fChain->SetBranchAddress("genantit_b_eta", &genantit_b_eta, &b_genantit_b_eta);
        fChain->SetBranchAddress("genantit_b_mass", &genantit_b_mass, &b_genantit_b_mass);
        fChain->SetBranchAddress("gent_w_pt", &gent_w_pt, &b_gent_w_pt);
        fChain->SetBranchAddress("gent_w_phi", &gent_w_phi, &b_gent_w_phi);
        fChain->SetBranchAddress("gent_w_eta", &gent_w_eta, &b_gent_w_eta);
        fChain->SetBranchAddress("gent_w_mass", &gent_w_mass, &b_gent_w_mass);
        fChain->SetBranchAddress("genantit_w_pt", &genantit_w_pt, &b_genantit_w_pt);
        fChain->SetBranchAddress("genantit_w_phi", &genantit_w_phi, &b_genantit_w_phi);
        fChain->SetBranchAddress("genantit_w_eta", &genantit_w_eta, &b_genantit_w_eta);
        fChain->SetBranchAddress("genantit_w_mass", &genantit_w_mass, &b_genantit_w_mass);
        fChain->SetBranchAddress("gent_w_q1_pt", &gent_w_q1_pt, &b_gent_w_q1_pt);
        fChain->SetBranchAddress("gent_w_q1_phi", &gent_w_q1_phi, &b_gent_w_q1_phi);
        fChain->SetBranchAddress("gent_w_q1_eta", &gent_w_q1_eta, &b_gent_w_q1_eta);
        fChain->SetBranchAddress("gent_w_q1_e", &gent_w_q1_e, &b_gent_w_q1_e);
        fChain->SetBranchAddress("gent_w_q1_pdg", &gent_w_q1_pdg, &b_gent_w_q1_pdg);
        fChain->SetBranchAddress("genantit_w_q1_pt", &genantit_w_q1_pt, &b_genantit_w_q1_pt);
        fChain->SetBranchAddress("genantit_w_q1_phi", &genantit_w_q1_phi, &b_genantit_w_q1_phi);
        fChain->SetBranchAddress("genantit_w_q1_eta", &genantit_w_q1_eta, &b_genantit_w_q1_eta);
        fChain->SetBranchAddress("genantit_w_q1_e", &genantit_w_q1_e, &b_genantit_w_q1_e);
        fChain->SetBranchAddress("genantit_w_q1_pdg", &genantit_w_q1_pdg, &b_genantit_w_q1_pdg);
        fChain->SetBranchAddress("gent_w_q2_pt", &gent_w_q2_pt, &b_gent_w_q2_pt);
        fChain->SetBranchAddress("gent_w_q2_phi", &gent_w_q2_phi, &b_gent_w_q2_phi);
        fChain->SetBranchAddress("gent_w_q2_eta", &gent_w_q2_eta, &b_gent_w_q2_eta);
        fChain->SetBranchAddress("gent_w_q2_e", &gent_w_q2_e, &b_gent_w_q2_e);
        fChain->SetBranchAddress("gent_w_q2_pdg", &gent_w_q2_pdg, &b_gent_w_q2_pdg);
        fChain->SetBranchAddress("genantit_w_q2_pt", &genantit_w_q2_pt, &b_genantit_w_q2_pt);
        fChain->SetBranchAddress("genantit_w_q2_phi", &genantit_w_q2_phi, &b_genantit_w_q2_phi);
        fChain->SetBranchAddress("genantit_w_q2_eta", &genantit_w_q2_eta, &b_genantit_w_q2_eta);
        fChain->SetBranchAddress("genantit_w_q2_e", &genantit_w_q2_e, &b_genantit_w_q2_e);
        fChain->SetBranchAddress("genantit_w_q2_pdg", &genantit_w_q2_pdg, &b_genantit_w_q2_pdg);
        fChain->SetBranchAddress("ptgenq1l",ptgenq1l, &b_ptgenq1l);
        fChain->SetBranchAddress("etagenq1l",etagenq1l, &b_etagenq1l);
        fChain->SetBranchAddress("phigenq1l",phigenq1l, &b_phigenq1l);
        fChain->SetBranchAddress("egenq1l",egenq1l, &b_egenq1l);
        fChain->SetBranchAddress("ptgenq1f",ptgenq1f, &b_ptgenq1f);
        fChain->SetBranchAddress("etagenq1f",etagenq1f, &b_etagenq1f);
        fChain->SetBranchAddress("phigenq1f",phigenq1f, &b_phigenq1f);
        fChain->SetBranchAddress("egenq1f",egenq1f, &b_egenq1f);
        fChain->SetBranchAddress("ptgenq2l",ptgenq2l, &b_ptgenq2l);
        fChain->SetBranchAddress("etagenq2l",etagenq2l, &b_etagenq2l);
        fChain->SetBranchAddress("phigenq2l",phigenq2l, &b_phigenq2l);
        fChain->SetBranchAddress("egenq2l",egenq2l, &b_egenq2l);
        fChain->SetBranchAddress("ptgenq2f",ptgenq2f, &b_ptgenq2f);
        fChain->SetBranchAddress("etagenq2f",etagenq2f, &b_etagenq2f);
        fChain->SetBranchAddress("phigenq2f",phigenq2f, &b_phigenq2f);
        fChain->SetBranchAddress("egenq2f",egenq2f, &b_egenq2f);
        fChain->SetBranchAddress("ptgenq3l",ptgenq3l, &b_ptgenq3l);
        fChain->SetBranchAddress("etagenq3l",etagenq3l, &b_etagenq3l);
        fChain->SetBranchAddress("phigenq3l",phigenq3l, &b_phigenq3l);
        fChain->SetBranchAddress("egenq3l",egenq3l, &b_egenq3l);
        fChain->SetBranchAddress("ptgenq3f",ptgenq3f, &b_ptgenq3f);
        fChain->SetBranchAddress("etagenq3f",etagenq3f, &b_etagenq3f);
        fChain->SetBranchAddress("phigenq3f",phigenq3f, &b_phigenq3f);
        fChain->SetBranchAddress("egenq3f",egenq3f, &b_egenq3f);
        fChain->SetBranchAddress("ptgenq4l",ptgenq4l, &b_ptgenq4l);
        fChain->SetBranchAddress("etagenq4l",etagenq4l, &b_etagenq4l);
        fChain->SetBranchAddress("phigenq4l",phigenq4l, &b_phigenq4l);
        fChain->SetBranchAddress("egenq4l",egenq4l, &b_egenq4l);
        fChain->SetBranchAddress("ptgenq4f",ptgenq4f, &b_ptgenq4f);
        fChain->SetBranchAddress("etagenq4f",etagenq4f, &b_etagenq4f);
        fChain->SetBranchAddress("phigenq4f",phigenq4f, &b_phigenq4f);
        fChain->SetBranchAddress("egenq4f",egenq4f, &b_egenq4f);
        fChain->SetBranchAddress("ptgenq5l",ptgenq5l, &b_ptgenq5l);
        fChain->SetBranchAddress("etagenq5l",etagenq5l, &b_etagenq5l);
        fChain->SetBranchAddress("phigenq5l",phigenq5l, &b_phigenq5l);
        fChain->SetBranchAddress("egenq5l",egenq5l, &b_egenq5l);
        fChain->SetBranchAddress("ptgenq5f",ptgenq5f, &b_ptgenq5f);
        fChain->SetBranchAddress("etagenq5f",etagenq5f, &b_etagenq5f);
        fChain->SetBranchAddress("phigenq5f",phigenq5f, &b_phigenq5f);
        fChain->SetBranchAddress("egenq5f",egenq5f, &b_egenq5f);
        fChain->SetBranchAddress("mothergenq1f",mothergenq1f, &b_mothergenq1f);
        fChain->SetBranchAddress("mothergenq2f",mothergenq2f, &b_mothergenq2f);
        fChain->SetBranchAddress("mothergenq3f",mothergenq3f, &b_mothergenq3f);
        fChain->SetBranchAddress("mothergenq4f",mothergenq4f, &b_mothergenq4f);
        fChain->SetBranchAddress("mothergenq5f",mothergenq5f, &b_mothergenq5f);
        fChain->SetBranchAddress("genw_q1_pt",genw_q1_pt, &b_genw_q1_pt);
        fChain->SetBranchAddress("genw_q1_eta",genw_q1_eta, &b_genw_q1_eta);
        fChain->SetBranchAddress("genw_q1_phi",genw_q1_phi, &b_genw_q1_phi);
        fChain->SetBranchAddress("genw_q1_e",genw_q1_e, &b_genw_q1_e);
        fChain->SetBranchAddress("genw_q1_pdg",genw_q1_pdg, &b_genw_q1_pdg);
        fChain->SetBranchAddress("genw_q2_pt",genw_q2_pt, &b_genw_q2_pt);
        fChain->SetBranchAddress("genw_q2_eta",genw_q2_eta, &b_genw_q2_eta);
        fChain->SetBranchAddress("genw_q2_phi",genw_q2_phi, &b_genw_q2_phi);
        fChain->SetBranchAddress("genw_q2_e",genw_q2_e, &b_genw_q2_e);
        fChain->SetBranchAddress("genw_q2_pdg",genw_q2_pdg, &b_genw_q2_pdg);
        fChain->SetBranchAddress("gent_w_tag", &gent_w_tag, &b_gent_w_tag);
        fChain->SetBranchAddress("genantit_w_tag", &genantit_w_tag, &b_genantit_w_tag);
        fChain->SetBranchAddress("mothergengf",mothergengf, &b_mothergengf);
        fChain->SetBranchAddress("mmothergengf",mmothergengf, &b_mmothergengf);
        fChain->SetBranchAddress("mmothergenq1f",mmothergenq1f, &b_mmothergenq1f);
        fChain->SetBranchAddress("mmothergenq2f",mmothergenq2f, &b_mmothergenq2f);
        fChain->SetBranchAddress("mmothergenq3f",mmothergenq3f, &b_mmothergenq3f);
        fChain->SetBranchAddress("mmothergenq4f",mmothergenq4f, &b_mmothergenq4f);
        fChain->SetBranchAddress("mmothergenq5f",mmothergenq5f, &b_mmothergenq5f);

    }
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("ls", &ls, &b_ls);
   fChain->SetBranchAddress("nLooseEle", &nLooseEle, &b_nLooseEle);
   fChain->SetBranchAddress("nLooseMu", &nLooseMu, &b_nLooseMu);
    fChain->SetBranchAddress("L1prefiring", &L1prefiring, &b_L1prefiring);
    fChain->SetBranchAddress("L1prefiringup", &L1prefiringup, &b_L1prefiringup);
    fChain->SetBranchAddress("L1prefiringdown", &L1prefiringdown, &b_L1prefiringdown);

//   fChain->SetBranchAddress("njets", &njets, &b_njets);
//   fChain->SetBranchAddress("nbtag", &nbtag, &b_nbtag);
//   fChain->SetBranchAddress("jet2_pt", &jet2_pt, &b_jet2_pt);
//   fChain->SetBranchAddress("jet2_btag", &jet2_btag, &b_jet2_btag);
//   fChain->SetBranchAddress("jet3_pt", &jet3_pt, &b_jet3_pt);
//   fChain->SetBranchAddress("jet3_btag", &jet3_btag, &b_jet3_btag);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("ptVlep", &ptVlep, &b_ptVlep);
   fChain->SetBranchAddress("yVlep", &yVlep, &b_yVlep);
   fChain->SetBranchAddress("phiVlep", &phiVlep, &b_phiVlep);
   fChain->SetBranchAddress("massVlep", &massVlep, &b_massVlep);
   fChain->SetBranchAddress("mtVlep", &mtVlep, &b_mtVlep);
//   fChain->SetBranchAddress("massVhad_gen", &massVhad_gen, &b_massVhad_gen);
    fChain->SetBranchAddress("jetAK8puppi_tau1", &jetAK8puppi_tau1, &b_jetAK8puppi_tau1);
    fChain->SetBranchAddress("jetAK8puppi_tau2", &jetAK8puppi_tau2, &b_jetAK8puppi_tau2);
    fChain->SetBranchAddress("jetAK8puppi_tau3", &jetAK8puppi_tau3, &b_jetAK8puppi_tau3);
    fChain->SetBranchAddress("jetAK8puppi_tau4", &jetAK8puppi_tau4, &b_jetAK8puppi_tau4);
    fChain->SetBranchAddress("jetAK8puppi_tau1_2", &jetAK8puppi_tau1_2, &b_jetAK8puppi_tau1_2);
    fChain->SetBranchAddress("jetAK8puppi_tau2_2", &jetAK8puppi_tau2_2, &b_jetAK8puppi_tau2_2);
    fChain->SetBranchAddress("jetAK8puppi_tau3_2", &jetAK8puppi_tau3_2, &b_jetAK8puppi_tau3_2);
    fChain->SetBranchAddress("jetAK8puppi_tau4_2", &jetAK8puppi_tau4_2, &b_jetAK8puppi_tau4_2);
    fChain->SetBranchAddress("jetAK8puppi_tau1_3", &jetAK8puppi_tau1_3, &b_jetAK8puppi_tau1_3);
    fChain->SetBranchAddress("jetAK8puppi_tau2_3", &jetAK8puppi_tau2_3, &b_jetAK8puppi_tau2_3);
    fChain->SetBranchAddress("jetAK8puppi_tau3_3", &jetAK8puppi_tau3_3, &b_jetAK8puppi_tau3_3);
    fChain->SetBranchAddress("jetAK8puppi_tau4_3", &jetAK8puppi_tau4_3, &b_jetAK8puppi_tau4_3);

   fChain->SetBranchAddress("lep", &lep, &b_lep);
   fChain->SetBranchAddress("channel", &channel, &b_channel);
   fChain->SetBranchAddress("ptlep1", &ptlep1, &b_ptlep1);
   fChain->SetBranchAddress("ptlep2", &ptlep2, &b_ptlep2);
   fChain->SetBranchAddress("etalep1", &etalep1, &b_etalep1);
   fChain->SetBranchAddress("etalep2", &etalep2, &b_etalep2);
   fChain->SetBranchAddress("philep1", &philep1, &b_philep1);
   fChain->SetBranchAddress("philep2", &philep2, &b_philep2);
   fChain->SetBranchAddress("met", &met, &b_met);
   fChain->SetBranchAddress("metPhi", &metPhi, &b_metPhi);
   fChain->SetBranchAddress("theWeight", &theWeight, &b_theWeight);
   fChain->SetBranchAddress("nump", &nump, &b_nump);
   fChain->SetBranchAddress("numm", &numm, &b_numm);
   fChain->SetBranchAddress("npT", &npT, &b_npT);
   fChain->SetBranchAddress("npIT", &npIT, &b_npIT);
   fChain->SetBranchAddress("nBX", &nBX, &b_nBX);



   fChain->SetBranchAddress("delPhilepmet", &delPhilepmet, &b_delPhilepmet);
   fChain->SetBranchAddress("deltaRlepjet", &deltaRlepjet, &b_deltaRlepjet);
   fChain->SetBranchAddress("delPhijetmet", &delPhijetmet, &b_delPhijetmet);
   fChain->SetBranchAddress("delPhijetlep", &delPhijetlep, &b_delPhijetlep);
    fChain->SetBranchAddress("deltaRlepjet_2", &deltaRlepjet_2, &b_deltaRlepjet_2);
    fChain->SetBranchAddress("delPhijetmet_2", &delPhijetmet_2, &b_delPhijetmet_2);
    fChain->SetBranchAddress("delPhijetlep_2", &delPhijetlep_2, &b_delPhijetlep_2);

   fChain->SetBranchAddress("vbftag", &vbftag, &b_vbftag);
   fChain->SetBranchAddress("IDLoose", &IDLoose, &b_IDLoose);
   fChain->SetBranchAddress("IDTight", &IDTight, &b_IDTight);
    fChain->SetBranchAddress("IDLoose_2", &IDLoose_2, &b_IDLoose_2);
    fChain->SetBranchAddress("IDTight_2", &IDTight_2, &b_IDTight_2);
    fChain->SetBranchAddress("IDLoose_3", &IDLoose_3, &b_IDLoose_3);
    fChain->SetBranchAddress("IDTight_3", &IDTight_3, &b_IDTight_3);

   //fChain->SetBranchAddress("isHighPt", &isHighPt, &b_isHighPt);
   //fChain->SetBranchAddress("isHEEP", &isHEEP, &b_isHEEP);
   fChain->SetBranchAddress("trackIso", &trackIso, &b_trackIso);
   fChain->SetBranchAddress("muchaiso", &muchaiso, &b_muchaiso);
   fChain->SetBranchAddress("muneuiso", &muneuiso, &b_muneuiso);
   fChain->SetBranchAddress("muphoiso", &muphoiso, &b_muphoiso);
   fChain->SetBranchAddress("muPU", &muPU, &b_muPU);
   fChain->SetBranchAddress("muisolation", &muisolation, &b_muisolation);
   fChain->SetBranchAddress("METraw_et", &METraw_et, &b_METraw_et);
   fChain->SetBranchAddress("METraw_phi", &METraw_phi, &b_METraw_phi);
   fChain->SetBranchAddress("METraw_sumEt", &METraw_sumEt, &b_METraw_sumEt);
   fChain->SetBranchAddress("MET_et", &MET_et, &b_MET_et);
   fChain->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
   fChain->SetBranchAddress("MET_sumEt", &MET_sumEt, &b_MET_sumEt);
   fChain->SetBranchAddress("jetAK8puppi_pt1", jetAK8puppi_pt1, &b_jetAK8puppi_pt1);
   fChain->SetBranchAddress("jetAK8puppi_eta1", jetAK8puppi_eta1, &b_jetAK8puppi_eta1);
   fChain->SetBranchAddress("jetAK8puppi_mass1", jetAK8puppi_mass1, &b_jetAK8puppi_mass1);
   fChain->SetBranchAddress("ptVlepJEC", &ptVlepJEC, &b_ptVlepJEC);
   fChain->SetBranchAddress("yVlepJEC", &yVlepJEC, &b_yVlepJEC);
   fChain->SetBranchAddress("phiVlepJEC", &phiVlepJEC, &b_phiVlepJEC);
   fChain->SetBranchAddress("massVlepJEC", &massVlepJEC, &b_massVlepJEC);

   fChain->SetBranchAddress("jetAK8puppi_sd", &jetAK8puppi_sd, &b_jetAK8puppi_sd);
   fChain->SetBranchAddress("jetAK8puppi_tau21", &jetAK8puppi_tau21, &b_jetAK8puppi_tau21);
   fChain->SetBranchAddress("jetAK8puppi_tau42", &jetAK8puppi_tau42, &b_jetAK8puppi_tau42);
   fChain->SetBranchAddress("jetAK8puppi_ptJEC", &jetAK8puppi_ptJEC, &b_jetAK8puppi_ptJEC);
   fChain->SetBranchAddress("jetAK8puppi_eta", &jetAK8puppi_eta, &b_jetAK8puppi_eta);
   fChain->SetBranchAddress("jetAK8puppi_phi", &jetAK8puppi_phi, &b_jetAK8puppi_phi);


    fChain->SetBranchAddress("m_jlv", &m_jlv, &b_m_jlv);
   fChain->SetBranchAddress("candMasspuppiJEC", &candMasspuppiJEC, &b_candMasspuppiJEC);

   fChain->SetBranchAddress("jetAK8puppi_sd_2", &jetAK8puppi_sd_2, &b_jetAK8puppi_sd_2);
   fChain->SetBranchAddress("jetAK8puppi_tau21_2", &jetAK8puppi_tau21_2, &b_jetAK8puppi_tau21_2);
   fChain->SetBranchAddress("jetAK8puppi_tau42_2", &jetAK8puppi_tau42_2, &b_jetAK8puppi_tau42_2);
   fChain->SetBranchAddress("jetAK8puppi_ptJEC_2", &jetAK8puppi_ptJEC_2, &b_jetAK8puppi_ptJEC_2);
   fChain->SetBranchAddress("jetAK8puppi_eta_2", &jetAK8puppi_eta_2, &b_jetAK8puppi_eta_2);
   fChain->SetBranchAddress("jetAK8puppi_phi_2", &jetAK8puppi_phi_2, &b_jetAK8puppi_phi_2);

    

    fChain->SetBranchAddress("jetAK8puppi_sd_3", &jetAK8puppi_sd_3, &b_jetAK8puppi_sd_3);
    fChain->SetBranchAddress("jetAK8puppi_tau21_3", &jetAK8puppi_tau21_3, &b_jetAK8puppi_tau21_3);
    fChain->SetBranchAddress("jetAK8puppi_tau42_3", &jetAK8puppi_tau42_3, &b_jetAK8puppi_tau42_3);
    fChain->SetBranchAddress("jetAK8puppi_ptJEC_3", &jetAK8puppi_ptJEC_3, &b_jetAK8puppi_ptJEC_3);
    fChain->SetBranchAddress("jetAK8puppi_eta_3", &jetAK8puppi_eta_3, &b_jetAK8puppi_eta_3);
    fChain->SetBranchAddress("jetAK8puppi_phi_3", &jetAK8puppi_phi_3, &b_jetAK8puppi_phi_3);


fChain->SetBranchAddress("massww",&massww);


//   fChain->SetBranchAddress("sdropJEC", &sdropJEC, &b_sdropJEC);
   fChain->SetBranchAddress("mtVlepJEC", &mtVlepJEC, &b_mtVlepJEC);
//   fChain->SetBranchAddress("delPhilepmetJEC", &delPhilepmetJEC, &b_delPhilepmetJEC);
//   fChain->SetBranchAddress("delPhijetmetJEC", &delPhijetmetJEC, &b_delPhijetmetJEC);
//   fChain->SetBranchAddress("delPhijetlepJEC", &delPhijetlepJEC, &b_delPhijetlepJEC);
   fChain->SetBranchAddress("HLT_Ele1", &HLT_Ele1, &b_HLT_Ele1);
   fChain->SetBranchAddress("HLT_Ele2", &HLT_Ele2, &b_HLT_Ele2);
   fChain->SetBranchAddress("HLT_Ele3", &HLT_Ele3, &b_HLT_Ele3);
   fChain->SetBranchAddress("HLT_Ele4", &HLT_Ele4, &b_HLT_Ele4);
   fChain->SetBranchAddress("HLT_Ele5", &HLT_Ele5, &b_HLT_Ele5);
   fChain->SetBranchAddress("HLT_Ele6", &HLT_Ele6, &b_HLT_Ele6);
   fChain->SetBranchAddress("HLT_Ele7", &HLT_Ele7, &b_HLT_Ele7);
   fChain->SetBranchAddress("HLT_Ele8", &HLT_Ele8, &b_HLT_Ele8);
   fChain->SetBranchAddress("HLT_Mu1", &HLT_Mu1, &b_HLT_Mu1);
   fChain->SetBranchAddress("HLT_Mu2", &HLT_Mu2, &b_HLT_Mu2);
   fChain->SetBranchAddress("HLT_Mu3", &HLT_Mu3, &b_HLT_Mu3);
   fChain->SetBranchAddress("HLT_Mu4", &HLT_Mu4, &b_HLT_Mu4);
   fChain->SetBranchAddress("HLT_Mu5", &HLT_Mu5, &b_HLT_Mu5);
   fChain->SetBranchAddress("HLT_Mu6", &HLT_Mu6, &b_HLT_Mu6);
   fChain->SetBranchAddress("HLT_Mu7", &HLT_Mu7, &b_HLT_Mu7);
   fChain->SetBranchAddress("HLT_Mu8", &HLT_Mu8, &b_HLT_Mu8);
   fChain->SetBranchAddress("HLT_Mu9", &HLT_Mu9, &b_HLT_Mu9);
   fChain->SetBranchAddress("HLT_Mu10", &HLT_Mu10, &b_HLT_Mu10);
   fChain->SetBranchAddress("HLT_Mu11", &HLT_Mu11, &b_HLT_Mu11);
   fChain->SetBranchAddress("HLT_Mu12", &HLT_Mu12, &b_HLT_Mu12);
   
   fChain->SetBranchAddress("passFilter_HBHE", &passFilter_HBHE, &b_passFilter_HBHE_);
   fChain->SetBranchAddress("passFilter_HBHEIso", &passFilter_HBHEIso, &b_passFilter_HBHEIso_);
   fChain->SetBranchAddress("passFilter_GlobalHalo", &passFilter_GlobalHalo, &b_passFilter_GlobalHalo_);
   fChain->SetBranchAddress("passFilter_ECALDeadCell", &passFilter_ECALDeadCell, &b_passFilter_ECALDeadCell_);
   fChain->SetBranchAddress("passFilter_GoodVtx", &passFilter_GoodVtx, &b_passFilter_GoodVtx_);
   fChain->SetBranchAddress("passFilter_EEBadSc", &passFilter_EEBadSc, &b_passFilter_EEBadSc_);
   fChain->SetBranchAddress("passFilter_badMuon", &passFilter_badMuon, &b_passFilter_badMuon_);
   fChain->SetBranchAddress("passFilter_badChargedHadron", &passFilter_badChargedHadron, &b_passFilter_badChargedHadron_);
    fChain->SetBranchAddress("passecalBadCalibFilterUpdate", &passecalBadCalibFilterUpdate, &b_passecalBadCalibFilterUpdate_);
fChain->SetBranchAddress("ak4jet_hf", ak4jet_hf, &b_ak4jet_hf);
   fChain->SetBranchAddress("ak4jet_pf", ak4jet_pf, &b_ak4jet_pf);
   fChain->SetBranchAddress("ak4jet_pt", ak4jet_pt, &b_ak4jet_pt);
   fChain->SetBranchAddress("ak4jet_pt_uncorr", ak4jet_pt_uncorr, &b_ak4jet_pt_uncorr);
   fChain->SetBranchAddress("ak4jet_eta", ak4jet_eta, &b_ak4jet_eta);
   fChain->SetBranchAddress("ak4jet_phi", ak4jet_phi, &b_ak4jet_phi);
   fChain->SetBranchAddress("ak4jet_e", ak4jet_e, &b_ak4jet_e);
   fChain->SetBranchAddress("ak4jet_dr", ak4jet_dr, &b_ak4jet_dr);
   fChain->SetBranchAddress("ak4jet_csv", ak4jet_csv, &b_ak4jet_csv);
   fChain->SetBranchAddress("ak4jet_icsv", ak4jet_icsv, &b_ak4jet_icsv);
    fChain->SetBranchAddress("ak4jet_deepcsvb", ak4jet_deepcsvb, &b_ak4jet_deepcsvb);
    fChain->SetBranchAddress("ak4jet_deepcsvbb", ak4jet_deepcsvbb, &b_ak4jet_deepcsvbb);
   fChain->SetBranchAddress("ak4jet_IDLoose", ak4jet_IDLoose, &b_ak4jet_IDLoose);
   fChain->SetBranchAddress("ak4jet_IDTight", ak4jet_IDTight, &b_ak4jet_IDTight);
    
    //fChain->SetBranchAddress("ak8sj11", &ak8sj11, &b_ak8sj11);
    //fChain->SetBranchAddress("ak8sj12", &ak8sj12, &b_ak8sj12);
    //fChain->SetBranchAddress("ak8sj13", &ak8sj13, &b_ak8sj13);
    //fChain->SetBranchAddress("ak8sj14", &ak8sj14, &b_ak8sj14);
    //fChain->SetBranchAddress("ak8sj15", &ak8sj15, &b_ak8sj15);
    //fChain->SetBranchAddress("ak8sj21", &ak8sj21, &b_ak8sj21);
    //fChain->SetBranchAddress("ak8sj22", &ak8sj22, &b_ak8sj22);
    //fChain->SetBranchAddress("ak8sj23", &ak8sj23, &b_ak8sj23);
    //fChain->SetBranchAddress("ak8sj24", &ak8sj24, &b_ak8sj24);
    //fChain->SetBranchAddress("ak8sj25", &ak8sj25, &b_ak8sj25);
    fChain->SetBranchAddress("jetAK8puppi_dnnTop", &jetAK8puppi_dnnTop, &b_jetAK8puppi_dnnTop);
    fChain->SetBranchAddress("jetAK8puppi_dnnW", &jetAK8puppi_dnnW  , &b_jetAK8puppi_dnnW  );
    fChain->SetBranchAddress("jetAK8puppi_dnnH4q", &jetAK8puppi_dnnH4q, &b_jetAK8puppi_dnnH4q);
    fChain->SetBranchAddress("jetAK8puppi_dnnTop_2", &jetAK8puppi_dnnTop_2, &b_jetAK8puppi_dnnTop_2);
    fChain->SetBranchAddress("jetAK8puppi_dnnW_2", &jetAK8puppi_dnnW_2  , &b_jetAK8puppi_dnnW_2  );
    fChain->SetBranchAddress("jetAK8puppi_dnnH4q_2", &jetAK8puppi_dnnH4q_2, &b_jetAK8puppi_dnnH4q_2);
    fChain->SetBranchAddress("jetAK8puppi_dnnTop_3", &jetAK8puppi_dnnTop_3, &b_jetAK8puppi_dnnTop_3);
    fChain->SetBranchAddress("jetAK8puppi_dnnW_3", &jetAK8puppi_dnnW_3  , &b_jetAK8puppi_dnnW_3  );
    fChain->SetBranchAddress("jetAK8puppi_dnnH4q_3", &jetAK8puppi_dnnH4q_3, &b_jetAK8puppi_dnnH4q_3);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrTop", &jetAK8puppi_dnnDecorrTop, &b_jetAK8puppi_dnnDecorrTop);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrW", &jetAK8puppi_dnnDecorrW  , &b_jetAK8puppi_dnnDecorrW  );
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrH4q", &jetAK8puppi_dnnDecorrH4q, &b_jetAK8puppi_dnnDecorrH4q);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrTop_2", &jetAK8puppi_dnnDecorrTop_2, &b_jetAK8puppi_dnnDecorrTop_2);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrW_2", &jetAK8puppi_dnnDecorrW_2  , &b_jetAK8puppi_dnnDecorrW_2  );
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrH4q_2", &jetAK8puppi_dnnDecorrH4q_2, &b_jetAK8puppi_dnnDecorrH4q_2);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrTop_3", &jetAK8puppi_dnnDecorrTop_3, &b_jetAK8puppi_dnnDecorrTop_3);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrW_3", &jetAK8puppi_dnnDecorrW_3  , &b_jetAK8puppi_dnnDecorrW_3  );
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrH4q_3", &jetAK8puppi_dnnDecorrH4q_3, &b_jetAK8puppi_dnnDecorrH4q_3);
    
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrbb", &jetAK8puppi_dnnDecorrbb, &b_jetAK8puppi_dnnDecorrbb);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrbbnog", &jetAK8puppi_dnnDecorrbbnog, &b_jetAK8puppi_dnnDecorrbbnog);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrcc", &jetAK8puppi_dnnDecorrcc, &b_jetAK8puppi_dnnDecorrcc);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrccnog", &jetAK8puppi_dnnDecorrccnog, &b_jetAK8puppi_dnnDecorrccnog);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrbb_2", &jetAK8puppi_dnnDecorrbb_2, &b_jetAK8puppi_dnnDecorrbb_2);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrbbnog_2", &jetAK8puppi_dnnDecorrbbnog_2, &b_jetAK8puppi_dnnDecorrbbnog_2);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrcc_2", &jetAK8puppi_dnnDecorrcc_2, &b_jetAK8puppi_dnnDecorrcc_2);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrccnog_2", &jetAK8puppi_dnnDecorrccnog_2, &b_jetAK8puppi_dnnDecorrccnog_2);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrbb_3", &jetAK8puppi_dnnDecorrbb_3, &b_jetAK8puppi_dnnDecorrbb_3);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrbbnog_3", &jetAK8puppi_dnnDecorrbbnog_3, &b_jetAK8puppi_dnnDecorrbbnog_3);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrcc_3", &jetAK8puppi_dnnDecorrcc_3, &b_jetAK8puppi_dnnDecorrcc_3);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrccnog_3", &jetAK8puppi_dnnDecorrccnog_3, &b_jetAK8puppi_dnnDecorrccnog_3);

    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrqcd", &jetAK8puppi_dnnDecorrqcd, &b_jetAK8puppi_dnnDecorrqcd);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrtop", &jetAK8puppi_dnnDecorrtop, &b_jetAK8puppi_dnnDecorrtop);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrz", &jetAK8puppi_dnnDecorrz, &b_jetAK8puppi_dnnDecorrz);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrw", &jetAK8puppi_dnnDecorrw, &b_jetAK8puppi_dnnDecorrw);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrzbb", &jetAK8puppi_dnnDecorrzbb, &b_jetAK8puppi_dnnDecorrzbb);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrhbb", &jetAK8puppi_dnnDecorrhbb, &b_jetAK8puppi_dnnDecorrhbb);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrh4q", &jetAK8puppi_dnnDecorrh4q, &b_jetAK8puppi_dnnDecorrh4q);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrqcd_2", &jetAK8puppi_dnnDecorrqcd_2, &b_jetAK8puppi_dnnDecorrqcd_2);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrtop_2", &jetAK8puppi_dnnDecorrtop_2, &b_jetAK8puppi_dnnDecorrtop_2);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrz_2", &jetAK8puppi_dnnDecorrz_2, &b_jetAK8puppi_dnnDecorrz_2);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrw_2", &jetAK8puppi_dnnDecorrw_2, &b_jetAK8puppi_dnnDecorrw_2);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrzbb_2", &jetAK8puppi_dnnDecorrzbb_2, &b_jetAK8puppi_dnnDecorrzbb_2);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrhbb_2", &jetAK8puppi_dnnDecorrhbb_2, &b_jetAK8puppi_dnnDecorrhbb);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrh4q_2", &jetAK8puppi_dnnDecorrh4q_2, &b_jetAK8puppi_dnnDecorrh4q_2);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrqcd_3", &jetAK8puppi_dnnDecorrqcd_3, &b_jetAK8puppi_dnnDecorrqcd_3);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrtop_3", &jetAK8puppi_dnnDecorrtop_3, &b_jetAK8puppi_dnnDecorrtop_3);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrz_3", &jetAK8puppi_dnnDecorrz_3, &b_jetAK8puppi_dnnDecorrz_3);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrw_3", &jetAK8puppi_dnnDecorrw_3, &b_jetAK8puppi_dnnDecorrw_3);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrzbb_3", &jetAK8puppi_dnnDecorrzbb_3, &b_jetAK8puppi_dnnDecorrzbb_3);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrhbb_3", &jetAK8puppi_dnnDecorrhbb_3, &b_jetAK8puppi_dnnDecorrhbb);
    fChain->SetBranchAddress("jetAK8puppi_dnnDecorrh4q_3", &jetAK8puppi_dnnDecorrh4q_3, &b_jetAK8puppi_dnnDecorrh4q_3);

    if(IsData>0){



        










   fChain->SetBranchAddress("gen_ele_pt", &gen_ele_pt, &b_gen_ele_pt);
   fChain->SetBranchAddress("gen_ele_eta", &gen_ele_eta, &b_gen_ele_eta);
   fChain->SetBranchAddress("gen_ele_phi", &gen_ele_phi, &b_gen_ele_phi);
   fChain->SetBranchAddress("gen_ele_e", &gen_ele_e, &b_gen_ele_e);
    fChain->SetBranchAddress("gen_ele_pt_2", &gen_ele_pt_2, &b_gen_ele_pt_2);
    fChain->SetBranchAddress("gen_ele_eta_2", &gen_ele_eta_2, &b_gen_ele_eta_2);
    fChain->SetBranchAddress("gen_ele_phi_2", &gen_ele_phi_2, &b_gen_ele_phi_2);
    fChain->SetBranchAddress("gen_ele_e_2", &gen_ele_e_2, &b_gen_ele_e_2);
    fChain->SetBranchAddress("gen_ele_pt_3", &gen_ele_pt_3, &b_gen_ele_pt_3);
    fChain->SetBranchAddress("gen_ele_eta_3", &gen_ele_eta_3, &b_gen_ele_eta_3);
    fChain->SetBranchAddress("gen_ele_phi_3", &gen_ele_phi_3, &b_gen_ele_phi_3);
    fChain->SetBranchAddress("gen_ele_e_3", &gen_ele_e_3, &b_gen_ele_e_3);

    fChain->SetBranchAddress("gen_mu_pt", &gen_mu_pt, &b_gen_mu_pt);
    fChain->SetBranchAddress("gen_mu_eta", &gen_mu_eta, &b_gen_mu_eta);
    fChain->SetBranchAddress("gen_mu_phi", &gen_mu_phi, &b_gen_mu_phi);
    fChain->SetBranchAddress("gen_mu_e", &gen_mu_e, &b_gen_mu_e);
    fChain->SetBranchAddress("gen_mu_pt_2", &gen_mu_pt_2, &b_gen_mu_pt_2);
    fChain->SetBranchAddress("gen_mu_eta_2", &gen_mu_eta_2, &b_gen_mu_eta_2);
    fChain->SetBranchAddress("gen_mu_phi_2", &gen_mu_phi_2, &b_gen_mu_phi_2);
    fChain->SetBranchAddress("gen_mu_e_2", &gen_mu_e_2, &b_gen_mu_e_2);
    fChain->SetBranchAddress("gen_mu_pt_3", &gen_mu_pt_3, &b_gen_mu_pt_3);
    fChain->SetBranchAddress("gen_mu_eta_3", &gen_mu_eta_3, &b_gen_mu_eta_3);
    fChain->SetBranchAddress("gen_mu_phi_3", &gen_mu_phi_3, &b_gen_mu_phi_3);
    fChain->SetBranchAddress("gen_mu_e_3", &gen_mu_e_3, &b_gen_mu_e_3);


























    












    












    //fChain->SetBranchAddress("puppi_softdropj1", &puppi_softdropj1, &b_puppi_softdropj1);
    //fChain->SetBranchAddress("puppi_softdropj2", &puppi_softdropj2, &b_puppi_softdropj2);

    fChain->SetBranchAddress("pttau", pttau, &b_pttau);
    fChain->SetBranchAddress("phitau", phitau, &b_phitau);
    fChain->SetBranchAddress("etatau", etatau, &b_etatau);
    fChain->SetBranchAddress("etau", etau, &b_etau);
    fChain->SetBranchAddress("pdgidtau", pdgidtau, &b_pdgidtau);
    fChain->SetBranchAddress("pttau_2", pttau_2, &b_pttau_2);
    fChain->SetBranchAddress("phitau_2", phitau_2, &b_phitau_2);
    fChain->SetBranchAddress("etatau_2", etatau_2, &b_etatau_2);
    fChain->SetBranchAddress("etau_2", etau_2, &b_etau_2);
    fChain->SetBranchAddress("pdgidtau_2", pdgidtau_2, &b_pdgidtau_2);
    fChain->SetBranchAddress("pttau_3", pttau_3, &b_pttau_3);
    fChain->SetBranchAddress("phitau_3", phitau_3, &b_phitau_3);
    fChain->SetBranchAddress("etatau_3", etatau_3, &b_etatau_3);
    fChain->SetBranchAddress("etau_3", etau_3, &b_etau_3);
    fChain->SetBranchAddress("pdgidtau_3", pdgidtau_3, &b_pdgidtau_3);

    fChain->SetBranchAddress("ptq", ptq, &b_ptq);
    fChain->SetBranchAddress("phiq", phiq, &b_phiq);
    fChain->SetBranchAddress("etaq", etaq, &b_etaq);
    fChain->SetBranchAddress("eq", eq, &b_eq);
    fChain->SetBranchAddress("pdgidq", pdgidq, &b_pdgidq);
    fChain->SetBranchAddress("ptq_2", ptq_2, &b_ptq_2);
    fChain->SetBranchAddress("phiq_2", phiq_2, &b_phiq_2);
    fChain->SetBranchAddress("etaq_2", etaq_2, &b_etaq_2);
    fChain->SetBranchAddress("eq_2", eq_2, &b_eq_2);
    fChain->SetBranchAddress("pdgidq_2", pdgidq_2, &b_pdgidq_2);
    fChain->SetBranchAddress("ptq_3", ptq_3, &b_ptq_3);
    fChain->SetBranchAddress("phiq_3", phiq_3, &b_phiq_3);
    fChain->SetBranchAddress("etaq_3", etaq_3, &b_etaq_3);
    fChain->SetBranchAddress("eq_3", eq_3, &b_eq_3);
    fChain->SetBranchAddress("pdgidq_3", pdgidq_3, &b_pdgidq_3);



   fChain->SetBranchAddress("genmatch_ele_pt", &genmatch_ele_pt, &b_genmatch_ele_pt);
   fChain->SetBranchAddress("genmatch_ele_eta", &genmatch_ele_eta, &b_genmatch_ele_eta);
   fChain->SetBranchAddress("genmatch_ele_phi", &genmatch_ele_phi, &b_genmatch_ele_phi);
   fChain->SetBranchAddress("genmatch_ele_e", &genmatch_ele_e, &b_genmatch_ele_e);
   fChain->SetBranchAddress("genmatch_ele_dr", &genmatch_ele_dr, &b_genmatch_ele_dr);
   fChain->SetBranchAddress("genmatch_mu_pt", &genmatch_mu_pt, &b_genmatch_mu_pt);
   fChain->SetBranchAddress("genmatch_mu_eta", &genmatch_mu_eta, &b_genmatch_mu_eta);
   fChain->SetBranchAddress("genmatch_mu_phi", &genmatch_mu_phi, &b_genmatch_mu_phi);
   fChain->SetBranchAddress("genmatch_mu_e", &genmatch_mu_e, &b_genmatch_mu_e);
   fChain->SetBranchAddress("genmatch_mu_dr", &genmatch_mu_dr, &b_genmatch_mu_dr);
   fChain->SetBranchAddress("gentop_pt", &gentop_pt, &b_gentop_pt);
   fChain->SetBranchAddress("gentop_eta", &gentop_eta, &b_gentop_eta);
   fChain->SetBranchAddress("gentop_phi", &gentop_phi, &b_gentop_phi);
   fChain->SetBranchAddress("gentop_mass", &gentop_mass, &b_gentop_mass);
   fChain->SetBranchAddress("genantitop_pt", &genantitop_pt, &b_genantitop_pt);
   fChain->SetBranchAddress("genantitop_eta", &genantitop_eta, &b_genantitop_eta);
   fChain->SetBranchAddress("genantitop_phi", &genantitop_phi, &b_genantitop_phi);
   fChain->SetBranchAddress("genantitop_mass", &genantitop_mass, &b_genantitop_mass);
















    








    }

   Notify();
}

Bool_t EDBR2PKUTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void EDBR2PKUTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

void EDBR2PKUTree::endJob() {
   fout->cd();
   ExTree->Write();
   fout->Write();
   fout->Close();
}

Int_t EDBR2PKUTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef EDBR2PKUTree_cxx
