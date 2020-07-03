#define EDBR2PKUTree_cxx
#include "EDBR2PKUTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <unistd.h>
#include <fcntl.h>

#define Pi 3.141593
//#include "BTagCalibrationStandalone.h"
//
vector<Double_t> generate_weights(TH1* data_npu_estimated, Int_t isForSynch){
	// see SimGeneral/MixingModule/python/mix_2015_25ns_Startup_PoissonOOTPU_cfi.pyy; copy and paste from there:
    //https://github.com/cms-sw/cmssw/blob/master/SimGeneral/MixingModule/python/mix_2017_25ns_WinterMC_PUScenarioV1_PoissonOOTPU_cfi.py

	const Double_t npu_probs[75] = {
1.78653e-05 ,2.56602e-05 ,5.27857e-05 ,8.88954e-05 ,0.000109362 ,0.000140973 ,0.000240998 ,0.00071209 ,0.00130121 ,0.00245255 ,0.00502589 ,0.00919534 ,0.0146697 ,0.0204126 ,0.0267586 ,0.0337697 ,0.0401478 ,0.0450159 ,0.0490577 ,0.0524855 ,0.0548159 ,0.0559937 ,0.0554468 ,0.0537687 ,0.0512055 ,0.0476713 ,0.0435312 ,0.0393107 ,0.0349812 ,0.0307413 ,0.0272425 ,0.0237115 ,0.0208329 ,0.0182459 ,0.0160712 ,0.0142498 ,0.012804 ,0.011571 ,0.010547 ,0.00959489 ,0.00891718 ,0.00829292 ,0.0076195 ,0.0069806 ,0.0062025 ,0.00546581 ,0.00484127 ,0.00407168 ,0.00337681 ,0.00269893 ,0.00212473 ,0.00160208 ,0.00117884 ,0.000859662 ,0.000569085 ,0.000365431 ,0.000243565 ,0.00015688 ,9.88128e-05 ,6.53783e-05 ,3.73924e-05 ,2.61382e-05 ,2.0307e-05 ,1.73032e-05 ,1.435e-05 ,1.36486e-05 ,1.35555e-05 ,1.37491e-05 ,1.34255e-05 ,1.33987e-05 ,1.34061e-05 ,1.34211e-05 ,1.34177e-05 ,1.32959e-05 ,1.33287e-05};
	if (isForSynch==0) { //OFFICIAL RECIPE
		vector<Double_t> result(75);
		Double_t s = 0.0;
		for(Int_t npu=0; npu<75; ++npu){
			Double_t npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));
			result[npu] = npu_estimated / npu_probs[npu];
			s += npu_estimated;
		}
		// normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
		for(Int_t npu=0; npu<75; ++npu){
			result[npu] /= s;
		}
		return result;
	}
	else { //THIS IS FOR THE SYNCH ONLY. THIS IS NOT THE OFFICIAL RECIPE!
		vector<Double_t> result(60);
		for(Int_t npu=0; npu<60; ++npu){
			if (data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu))==NULL)
			  result[npu] = 0.;
			else {
				Double_t npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));            
				result[npu] = npu_estimated;
			}
		}
		return result;
	}

}

double SumPt(double ptphi[5][2],int n){
    double pt=-99;
    double x1=0,x2=0;
    for(int i=0;i<n;i++){x1+=ptphi[i][0]*TMath::Cos(ptphi[i][1]);x2+=ptphi[i][0]*TMath::Sin(ptphi[i][1]);}
    //cout<<"  sum pt "<<ptphi[0][0]<<"   "<<x1<<"   "<<x2<<endl;
    if(n>=1)pt=sqrt(x1*x1+x2*x2);
    return pt;
}
int endpoint(int nl,int np,int p1or2){
    int p1,p2;
    if(np==5){
        if(nl==0){p1=0;p2=1;}
        if(nl==1){p1=0;p2=2;}
        if(nl==2){p1=0;p2=3;}
        if(nl==3){p1=0;p2=4;}
        if(nl==4){p1=1;p2=2;}
        if(nl==5){p1=1;p2=3;}
        if(nl==6){p1=1;p2=4;}
        if(nl==7){p1=2;p2=3;}
        if(nl==8){p1=2;p2=4;}
        if(nl==9){p1=3;p2=4;}
        }
    if(np==4){
        if(nl==0){p1=0;p2=1;}
        if(nl==1){p1=0;p2=2;}
        if(nl==2){p1=0;p2=3;}
        if(nl==3){p1=1;p2=2;}
        if(nl==4){p1=1;p2=3;}
        if(nl==5){p1=2;p2=3;}
    }
    if(np==3){
        if(nl==0){p1=0;p2=1;}
        if(nl==1){p1=0;p2=2;}
        if(nl==2){p1=1;p2=2;}
    }
    if(np==2){
        if(nl==0){p1=0;p2=1;}
    }
    if(p1or2==2) p1=p2;
    return p1;
}
Double_t bsv_deep (Int_t cud, Double_t x ) // (b tag middle workpoint) cud=1,2,3 for central,up,down; x for pt
{
    double result=1.0;
    if (cud==1) {  //central
        result=0.653526*((1.+(0.220245*x))/(1.+(0.14383*x)));
    }
    else if (cud==2) { //up
        if(x<20) {result=1;}
        else if(x<30)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))+0.043795019388198853;}
        else if(x<50)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))+0.015845479443669319;}
        else if(x<70)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))+0.014174085110425949;}
        else if(x<100)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))+0.013200919143855572;}
        else if(x<140)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))+0.012912030331790447;}
        else if(x<200)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))+0.019475525245070457;}
        else if(x<300)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))+0.01628459244966507;}
        else if(x<600)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))+0.034840557724237442;}
        else if(x<1000)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))+0.049875054508447647;}
        else {result=1;}
    }
    else if (cud==3) { //down
        if(x<20) {result=1;}
        else if(x<30)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))-0.043795019388198853;}
        else if(x<50)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))-0.015845479443669319;}
        else if(x<70)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))-0.014174085110425949;}
        else if(x<100)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))-0.013200919143855572;}
        else if(x<140)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))-0.012912030331790447;}
        else if(x<200)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))-0.019475525245070457;}
        else if(x<300)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))-0.01628459244966507;}
        else if(x<600)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))-0.034840557724237442;}
        else if(x<1000)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))-0.049875054508447647;}
        else {result=1;}
    }
    return result;
}
Double_t csv_deep (Int_t cud, Double_t x ) // (b tag middle workpoint) cud=1,2,3 for central,up,down; x for pt
{
    double result=1.0;
    if (cud==1) {  //central
        result=0.653526*((1.+(0.220245*x))/(1.+(0.14383*x)));
    }
    else if (cud==2) { //up
        if(x<20) {result=1;}
        else if(x<30)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))+0.13138505816459656;}
        else if(x<50)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))+0.047536440193653107;}
        else if(x<70)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))+0.042522255331277847;}
        else if(x<100)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))+0.039602756500244141;}
        else if(x<140)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))+0.038736090064048767;}
        else if(x<200)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))+0.058426573872566223;}
        else if(x<300)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))+0.048853777348995209;}
        else if(x<600)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))+0.10452167689800262;}
        else if(x<1000)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))+0.14962516725063324;}
        else {result=1;}
    }
    else if (cud==3) { //down
        if(x<20) {result=1;}
        else if(x<30)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))-0.13138505816459656;}
        else if(x<50)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))-0.047536440193653107;}
        else if(x<70)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))-0.042522255331277847;}
        else if(x<100)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))-0.039602756500244141;}
        else if(x<140)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))-0.038736090064048767;}
        else if(x<200)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))-0.058426573872566223;}
        else if(x<300)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))-0.048853777348995209;}
        else if(x<600)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))-0.10452167689800262;}
        else if(x<1000)  {result=(0.653526*((1.+(0.220245*x))/(1.+(0.14383*x))))-0.14962516725063324;}
        else {result=1;}
    }
    return result;
}
Double_t lsv_deep (Int_t cud, Double_t x ) // (b tag middle workpoint) cud=1,2,3 for central,up,down; x for pt
{
    double result=1.0;
    if (cud==1) {  //central
        result=1.09286+-0.00052597*x+1.88225e-06*x*x+-1.27417e-09*x*x*x;
    }
    else if (cud==2) { //up
        if(x<20) {result=1;}
        else if(x<1000)  {result=(1.09286+-0.00052597*x+1.88225e-06*x*x+-1.27417e-09*x*x*x)*(1+(0.101915+0.000192134*x+-1.94974e-07*x*x));}
        else {result=1;}
    }
    else if (cud==3) { //down
        if(x<20) {result=1;}
        else if(x<1000)  {result=(1.09286+-0.00052597*x+1.88225e-06*x*x+-1.27417e-09*x*x*x)*(1-(0.101915+0.000192134*x+-1.94974e-07*x*x));}
        else {result=1;}
    }
    return result;
}
Double_t bsv_los_deep (Int_t cud, Double_t x ) // (b tag loose workpoint) cud=1,2,3 for central,up,down; x for pt
{
    double result=1.0;
    if (cud==1) {  //central
        result=0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x)));
    }
    else if (cud==2) { //up
        if(x<20) {result=1;}
        else if(x<30)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))+0.039079215377569199;}
        else if(x<50)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))+0.014356007799506187;}
        else if(x<70)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))+0.012236535549163818;}
        else if(x<100)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))+0.011882896535098553;}
        else if(x<140)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))+0.011785224080085754;}
        else if(x<200)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))+0.012215510942041874;}
        else if(x<300)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))+0.014544816687703133;}
        else if(x<600)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))+0.026683652773499489;}
        else if(x<1000)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))+0.055047694593667984;}
        else {result=1;}
    }
    else if (cud==3) { //down
        if(x<20) {result=1;}
        else if(x<30)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))-0.039079215377569199;}
        else if(x<50)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))-0.014356007799506187;}
        else if(x<70)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))-0.012236535549163818;}
        else if(x<100)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))-0.011882896535098553;}
        else if(x<140)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))-0.011785224080085754;}
        else if(x<200)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))-0.012215510942041874;}
        else if(x<300)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))-0.014544816687703133;}
        else if(x<600)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))-0.026683652773499489;}
        else if(x<1000)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))-0.055047694593667984;}
        else {result=1;}
    }
    return result;
}
Double_t csv_los_deep (Int_t cud, Double_t x ) // (b tag loose workpoint) cud=1,2,3 for central,up,down; x for pt
{
    double result=1.0;
    if (cud==1) {  //central
        result=0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x)));
    }
    else if (cud==2) { //up
        if(x<20) {result=1;}
        else if(x<30)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))+0.097698040306568146;}
        else if(x<50)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))+0.035890020430088043;}
        else if(x<70)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))+0.030591338872909546;}
        else if(x<100)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))+0.029707241803407669;}
        else if(x<140)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))+0.029463060200214386;}
        else if(x<200)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))+0.030538776889443398;}
        else if(x<300)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))+0.036362040787935257;}
        else if(x<600)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))+0.066709131002426147;}
        else if(x<1000)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))+0.13761924207210541;}
        else {result=1;}
    }
    else if (cud==3) { //down
        if(x<20) {result=1;}
        else if(x<30)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))-0.097698040306568146;}
        else if(x<50)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))-0.035890020430088043;}
        else if(x<70)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))-0.030591338872909546;}
        else if(x<100)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))-0.029707241803407669;}
        else if(x<140)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))-0.029463060200214386;}
        else if(x<200)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))-0.030538776889443398;}
        else if(x<300)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))-0.036362040787935257;}
        else if(x<600)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))-0.066709131002426147;}
        else if(x<1000)  {result=(0.933791*((1.+(0.0115268*x))/(1.+(0.0103699*x))))-0.13761924207210541;}
        else {result=1;}
    }
    return result;
}
Double_t lsv_los_deep (Int_t cud, Double_t x ) // (b tag loose workpoint) cud=1,2,3 for central,up,down; x for pt
{
    double result=1.0;
    if (cud==1) {  //central
        result=1.06337+-0.000276004*x+1.25504e-06*x*x+-8.9312e-10*x*x*x;
    }
    else if (cud==2) { //up
        if(x<20) {result=1;}
        else if(x<1000)  {result=(1.06337+-0.000276004*x+1.25504e-06*x*x+-8.9312e-10*x*x*x)*(1+(0.0421943+5.30087e-05*x+-6.87049e-08*x*x));}
        else {result=1;}
    }
    else if (cud==3) { //down
        if(x<20) {result=1;}
        else if(x<1000)  {result=(1.06337+-0.000276004*x+1.25504e-06*x*x+-8.9312e-10*x*x*x)*(1-(0.0421943+5.30087e-05*x+-6.87049e-08*x*x));}
        else {result=1;}
    }
    return result;
}
Double_t bsv_tig_deep (Int_t cud, Double_t x ) // (b tag tight workpoint) cud=1,2,3 for central,up,down; x for pt
{
    double result=1.0;
    if (cud==1) {  //central
        result=0.495614*((1.+(0.371043*x))/(1.+(0.184225*x)));
    }
    else if (cud==2) { //up
        if(x<20) {result=1;}
        else if(x<30)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))+0.064842469990253448;}
        else if(x<50)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))+0.015950050204992294;}
        else if(x<70)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))+0.014914392493665218;}
        else if(x<100)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))+0.014344207942485809;}
        else if(x<140)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))+0.014204045757651329;}
        else if(x<200)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))+0.014788268133997917;}
        else if(x<300)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))+0.020897330716252327;}
        else if(x<600)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))+0.03883708268404007;}
        else if(x<1000)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))+0.053850434720516205;}
        else {result=1;}
    }
    else if (cud==3) { //down
        if(x<20) {result=1;}
        else if(x<30)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))-0.064842469990253448;}
        else if(x<50)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))-0.015950050204992294;}
        else if(x<70)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))-0.014914392493665218;}
        else if(x<100)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))-0.014344207942485809;}
        else if(x<140)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))-0.014204045757651329;}
        else if(x<200)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))-0.014788268133997917;}
        else if(x<300)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))-0.020897330716252327;}
        else if(x<600)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))-0.03883708268404007;}
        else if(x<1000)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))-0.053850434720516205;}
        else {result=1;}
    }
    return result;
}
Double_t csv_tig_deep (Int_t cud, Double_t x ) // (b tag tight workpoint) cud=1,2,3 for central,up,down; x for pt
{
    double result=1.0;
    if (cud==1) {  //central
        result=0.495614*((1.+(0.371043*x))/(1.+(0.184225*x)));
    }
    else if (cud==2) { //up
        if(x<20) {result=1;}
        else if(x<30)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))+0.22694864869117737;}
        else if(x<50)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))+0.055825173854827881;}
        else if(x<70)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))+0.052200373262166977;}
        else if(x<100)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))+0.050204727798700333;}
        else if(x<140)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))+0.049714159220457077;}
        else if(x<200)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))+0.051758937537670135;}
        else if(x<300)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))+0.073140658438205719;}
        else if(x<600)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))+0.13592979311943054;}
        else if(x<1000)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))+0.18847651779651642;}
        else {result=1;}
    }
    else if (cud==3) { //down
        if(x<20) {result=1;}
        else if(x<30)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))-0.22694864869117737;}
        else if(x<50)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))-0.055825173854827881;}
        else if(x<70)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))-0.052200373262166977;}
        else if(x<100)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))-0.050204727798700333;}
        else if(x<140)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))-0.049714159220457077;}
        else if(x<200)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))-0.051758937537670135;}
        else if(x<300)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))-0.073140658438205719;}
        else if(x<600)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))-0.13592979311943054;}
        else if(x<1000)  {result=(0.495614*((1.+(0.371043*x))/(1.+(0.184225*x))))-0.18847651779651642;}
        else {result=1;}
    }
    return result;
}
Double_t lsv_tig_deep (Int_t cud, Double_t x ) // (b tag tight workpoint) cud=1,2,3 for central,up,down; x for pt
{
    double result=1.0;
    if (cud==1) {  //central
        result=0.963784+0.423962/sqrt(x);
    }
    else if (cud==2) { //up
        if(x<20) {result=1;}
        else if(x<1000)  {result=(0.963784+0.423962/sqrt(x))*(1+(0.21914+0.000183013*x+-1.26742e-07*x*x));}
        else {result=1;}
    }
    else if (cud==3) { //down
        if(x<20) {result=1;}
        else if(x<1000)  {result=(0.963784+0.423962/sqrt(x))*(1-(0.21914+0.000183013*x+-1.26742e-07*x*x));}
        else {result=1;}
    }
    return result;
}

Double_t bsv (Int_t cud, Double_t x ) // (b tag middle workpoint) cud=1,2,3 for central,up,down; x for pt
{
    double result=1.0;
    if (cud==1) {  //central
        result=0.561694*((1.+(0.31439*x))/(1.+(0.17756*x)));
    }
    else if (cud==2) { //up
        if(x<20) {result=1;}
        else if(x<30)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.040213499218225479;}
        else if(x<50)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.014046305790543556;}
        else if(x<70)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.012372690252959728;}
        else if(x<100)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.012274007312953472;}
        else if(x<140)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.011465910822153091;}
        else if(x<200)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.012079551815986633;}
        else if(x<300)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.014995276927947998;}
        else if(x<600)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.021414462476968765;}
        else if(x<1000)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.032377112656831741;}
        else {result=1;}
    }
    else if (cud==3) { //down
        if(x<20) {result=1;}
        else if(x<30)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.040213499218225479;}
        else if(x<50)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.014046305790543556;}
        else if(x<70)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.012372690252959728;}
        else if(x<100)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.012274007312953472;}
        else if(x<140)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.011465910822153091;}
        else if(x<200)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.012079551815986633;}
        else if(x<300)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.014995276927947998;}
        else if(x<600)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.021414462476968765;}
        else if(x<1000)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.032377112656831741;}
        else {result=1;}
    }
    return result;
}
Double_t csv (Int_t cud, Double_t x ) // (b tag middle workpoint) cud=1,2,3 for central,up,down; x for pt
{
    double result=1.0;
    if (cud==1) {  //central
        result=0.561694*((1.+(0.31439*x))/(1.+(0.17756*x)));
    }
    else if (cud==2) { //up
        if(x<20) {result=1;}
        else if(x<30)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.12064050137996674;}
        else if(x<50)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.042138919234275818;}
        else if(x<70)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.03711806982755661;}
        else if(x<100)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.036822021007537842;}
        else if(x<140)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.034397732466459274;}
        else if(x<200)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.0362386554479599;}
        else if(x<300)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.044985830783843994;}
        else if(x<600)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.064243391156196594;}
        else if(x<1000)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))+0.097131341695785522;}
        else {result=1;}
    }
    else if (cud==3) { //down
        if(x<20) {result=1;}
        else if(x<30)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.12064050137996674;}
        else if(x<50)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.042138919234275818;}
        else if(x<70)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.03711806982755661;}
        else if(x<100)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.036822021007537842;}
        else if(x<140)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.034397732466459274;}
        else if(x<200)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.0362386554479599;}
        else if(x<300)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.044985830783843994;}
        else if(x<600)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.064243391156196594;}
        else if(x<1000)  {result=(0.561694*((1.+(0.31439*x))/(1.+(0.17756*x))))-0.097131341695785522;}
        else {result=1;}
    }
    return result;
}
Double_t lsv (Int_t cud, Double_t x ) // (b tag middle workpoint) cud=1,2,3 for central,up,down; x for pt
{
    double result=1.0;
    if (cud==1) {  //central
        result=1.0589+0.000382569*x+-2.4252e-07*x*x+2.20966e-10*x*x*x;
    }
    else if (cud==2) { //up
        if(x<20) {result=1;}
        else if(x<1000)  {result=(1.0589+0.000382569*x+-2.4252e-07*x*x+2.20966e-10*x*x*x)*(1+(0.100485+3.95509e-05*x+-4.90326e-08*x*x));}
        else {result=1;}
    }
    else if (cud==3) { //down
        if(x<20) {result=1;}
        else if(x<1000)  {result=(1.0589+0.000382569*x+-2.4252e-07*x*x+2.20966e-10*x*x*x)*(1-(0.100485+3.95509e-05*x+-4.90326e-08*x*x));}
        else {result=1;}
    }
    return result;
}
Double_t bsv_los (Int_t cud, Double_t x ) // (b tag loose workpoint) cud=1,2,3 for central,up,down; x for pt
{
    double result=1.0;
    if (cud==1) {  //central
        result=0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x)));
    }
    else if (cud==2) { //up
        if(x<20) {result=1;}
        else if(x<30)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.025381835177540779;}
        else if(x<50)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.012564006261527538;}
        else if(x<70)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.011564776301383972;}
        else if(x<100)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.011248723603785038;}
        else if(x<140)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.010811596177518368;}
        else if(x<200)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.010882497765123844;}
        else if(x<300)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.013456921093165874;}
        else if(x<600)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.017094610258936882;}
        else if(x<1000)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.02186630479991436;}
        else {result=1;}
    }
    else if (cud==3) { //down
        if(x<20) {result=1;}
        else if(x<30)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.025381835177540779;}
        else if(x<50)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.012564006261527538;}
        else if(x<70)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.011564776301383972;}
        else if(x<100)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.011248723603785038;}
        else if(x<140)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.010811596177518368;}
        else if(x<200)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.010882497765123844;}
        else if(x<300)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.013456921093165874;}
        else if(x<600)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.017094610258936882;}
        else if(x<1000)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.02186630479991436;}
        else {result=1;}
    }
    return result;
}
Double_t csv_los (Int_t cud, Double_t x ) // (b tag loose workpoint) cud=1,2,3 for central,up,down; x for pt
{
    double result=1.0;
    if (cud==1) {  //central
        result=0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x)));
    }
    else if (cud==2) { //up
        if(x<20) {result=1;}
        else if(x<30)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.063454590737819672;}
        else if(x<50)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.031410016119480133;}
        else if(x<70)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.02891194075345993;}
        else if(x<100)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.028121808543801308;}
        else if(x<140)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.027028990909457207;}
        else if(x<200)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.027206243947148323;}
        else if(x<300)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.033642303198575974;}
        else if(x<600)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.04273652657866478;}
        else if(x<1000)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))+0.054665762931108475;}
        else {result=1;}
    }
    else if (cud==3) { //down
        if(x<20) {result=1;}
        else if(x<30)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.063454590737819672;}
        else if(x<50)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.031410016119480133;}
        else if(x<70)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.02891194075345993;}
        else if(x<100)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.028121808543801308;}
        else if(x<140)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.027028990909457207;}
        else if(x<200)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.027206243947148323;}
        else if(x<300)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.033642303198575974;}
        else if(x<600)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.04273652657866478;}
        else if(x<1000)  {result=(0.887973*((1.+(0.0523821*x))/(1.+(0.0460876*x))))-0.054665762931108475;}
        else {result=1;}
    }
    return result;
}
Double_t lsv_los (Int_t cud, Double_t x ) // (b tag loose workpoint) cud=1,2,3 for central,up,down; x for pt
{
    double result=1.0;
    if (cud==1) {  //central
        result=1.13904+-0.000594946*x+1.97303e-06*x*x+-1.38194e-09*x*x*x;
    }
    else if (cud==2) { //up
        if(x<20) {result=1;}
        else if(x<1000)  {result=(1.13904+-0.000594946*x+1.97303e-06*x*x+-1.38194e-09*x*x*x)*(1+(0.0996438+-8.33354e-05*x+4.74359e-08*x*x));}
        else {result=1;}
    }
    else if (cud==3) { //down
        if(x<20) {result=1;}
        else if(x<1000)  {result=(1.13904+-0.000594946*x+1.97303e-06*x*x+-1.38194e-09*x*x*x)*(1-(0.0996438+-8.33354e-05*x+4.74359e-08*x*x));}
        else {result=1;}
    }
    return result;
}
Double_t bsv_tig (Int_t cud, Double_t x ) // (b tag tight workpoint) cud=1,2,3 for central,up,down; x for pt
{
    double result=1.0;
    if (cud==1) {  //central
        result=0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x)));
    }
    else if (cud==2) { //up
        if(x<20) {result=1;}
        else if(x<30)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.033732704818248749;}
        else if(x<50)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.01562843844294548;}
        else if(x<70)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.013530348427593708;}
        else if(x<100)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.013609844259917736;}
        else if(x<140)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.013236711733043194;}
        else if(x<200)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.013806583359837532;}
        else if(x<300)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.019633084535598755;}
        else if(x<600)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.030928170308470726;}
        else if(x<1000)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.052857179194688797;}
        else {result=1;}
    }
    else if (cud==3) { //down
        if(x<20) {result=1;}
        else if(x<30)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.033732704818248749;}
        else if(x<50)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.01562843844294548;}
        else if(x<70)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.013530348427593708;}
        else if(x<100)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.013609844259917736;}
        else if(x<140)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.013236711733043194;}
        else if(x<200)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.013806583359837532;}
        else if(x<300)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.019633084535598755;}
        else if(x<600)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.030928170308470726;}
        else if(x<1000)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.052857179194688797;}
        else {result=1;}
    }
    return result;
}
Double_t csv_tig (Int_t cud, Double_t x ) // (b tag tight workpoint) cud=1,2,3 for central,up,down; x for pt
{
    double result=1.0;
    if (cud==1) {  //central
        result=0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x)));
    }
    else if (cud==2) { //up
        if(x<20) {result=1;}
        else if(x<30)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.11806446313858032;}
        else if(x<50)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.054699532687664032;}
        else if(x<70)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.047356218099594116;}
        else if(x<100)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.047634456306695938;}
        else if(x<140)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.04632849246263504;}
        else if(x<200)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.048323042690753937;}
        else if(x<300)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.068715795874595642;}
        else if(x<600)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.10824859887361526;}
        else if(x<1000)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))+0.18500012159347534;}
        else {result=1;}
    }
    else if (cud==3) { //down
        if(x<20) {result=1;}
        else if(x<30)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.11806446313858032;}
        else if(x<50)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.054699532687664032;}
        else if(x<70)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.047356218099594116;}
        else if(x<100)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.047634456306695938;}
        else if(x<140)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.04632849246263504;}
        else if(x<200)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.048323042690753937;}
        else if(x<300)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.068715795874595642;}
        else if(x<600)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.10824859887361526;}
        else if(x<1000)  {result=(0.817647*((1.+(0.038703*x))/(1.+(0.0312388*x))))-0.18500012159347534;}
        else {result=1;}
    }
    return result;
}
Double_t lsv_tig (Int_t cud, Double_t x ) // (b tag tight workpoint) cud=1,2,3 for central,up,down; x for pt
{
    double result=1.0;
    if (cud==1) {  //central
        result=0.971945+163.215/(x*x)+0.000517836*x;
    }
    else if (cud==2) { //up
        if(x<20) {result=1;}
        else if(x<1000)  {result=(0.971945+163.215/(x*x)+0.000517836*x)*(1+(0.291298+-0.000222983*x+1.69699e-07*x*x));}
        else {result=1;}
    }
    else if (cud==3) { //down
        if(x<20) {result=1;}
        else if(x<1000)  {result=(0.971945+163.215/(x*x)+0.000517836*x)*(1-(0.291298+-0.000222983*x+1.69699e-07*x*x));}
        else {result=1;}
    }
    return result;
}

void EDBR2PKUTree::Loop(TString channelname, Double_t XS, TTree *treew, Int_t IsData,Int_t IsMET) {
    cout<<endl<<"IsMET  "<<IsMET<<"   "<<IsData<<endl;
	std::vector<Double_t> weights_pu1; //these are made with our recipe
	std::vector<Double_t> weights_pu2; //these are made with the official recipe
	TFile* pileupFile1 = TFile::Open("MyDataPileupHistogram.root");//pileupDataRun2016BH_63mb_80X.root");  
	TH1F* pileupHisto1 = (TH1F*)pileupFile1->Get("pileup");  
	weights_pu1 = generate_weights(pileupHisto1,0);
	pileupFile1->Close();

	//  TFile* pileupFile2 = TFile::Open("puweights.root");  
	TFile* pileupFile2 = TFile::Open("PUxSynch.root");  
	TH1F *pileupHisto2 = (TH1F*)pileupFile2->Get("puweights");
	weights_pu2 = generate_weights(pileupHisto2,1);
	pileupFile2->Close();

	//TFile * input1 = new TFile ("puweights.root");
	//TH1F* hR1= (TH1F*)input1->Get("puweights");
	//zixu
	TFile * input1 = new TFile ("puweight.root");	
	TH1F* hR1= (TH1F*)input1->Get("h2");
	//TFile * input1 = new TFile ("test_mu.root");
	//TH1F* hR1= (TH1F*)input1->Get("hRatio"); //"pileup");//hRatio");


	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();

	Double_t n_deltaRlepjet = 0; 
	Double_t n_delPhijetlep = 0; 
	Double_t ntau = 0;
	Double_t number_qq = 0; 
	Double_t nmassVhad = 0; 
	Double_t nptVlepJEC = 0;
	Double_t nID_e = 0;
	Double_t npt_e = 0;
	Double_t nmet_e = 0; 
	Double_t nnum_bJet_e = 0; 
	Double_t n_delPhijetmet = 0; 

	Double_t nID_mu = 0;
	Double_t npt_mu = 0;
	Double_t nmet_mu = 0; 
	Double_t nnum_bJet_mu = 0; 
	//Double_t nbtb_mu = 0; 

	Double_t nptVhad = 0;
	Double_t yields = 0;
	//TLorentzVector jetV, genjetV;
	//some constants inside this analysis
	Double_t pi_2=1.57079632679;
    Double_t theWeightw;
    //treew->GetBranchAdrress("theWeight",&theWeightw);
	Long64_t npp = treew->GetEntries("theWeight>0.");
	Long64_t nmm = treew->GetEntries("theWeight<0.");
	cout<<"npp="<<npp<<" nmm="<<nmm<<" totaleventnumber="<<(treew->GetEntries())<<"   number of events passed= "<<nentries<<endl;

	Double_t nn;
	Double_t eff_and_pu_Weight;
	Double_t eff_and_pu_Weight1;
	Double_t Identical_lumiWeight = XS;//All the events inside a sample are same lumiweight
	//Double_t Identical_lumiWeight = XS/totaleventnumber;//All the events inside a sample are same lumiweight

    TFile * elereco = new TFile ("EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root");
    TH2D* sfelereco= (TH2D*) elereco->Get("EGamma_SF2D");

    TFile * elehlt = new TFile ("egammaEffi.txt_EGM2D.root");
    TH2D* sfelehltbar= (TH2D*) elehlt->Get("SF_TH2F_Barrel");
    TH2D* sfelehltend= (TH2D*) elehlt->Get("SF_TH2F_EndCap");

    TFile * muhlt = new TFile ("EfficienciesAndSF_RunBtoF.root");
    TH2D* sfmuhlt= (TH2D*) muhlt->Get("Mu50_OR_TkMu50_PtEtaBins/efficienciesDATA/abseta_pt_DATA");
    TH2D* sfmuhlt2= (TH2D*) muhlt->Get("Mu50_OR_TkMu50_PtEtaBins/efficienciesMC/abseta_pt_MC");

    TFile * muhlt1 = new TFile ("EfficienciesAndSF_Period4.root");
    TH2D* sfmuhlt3= (TH2D*) muhlt1->Get("Mu50_OR_TkMu50_PtEtaBins/efficienciesDATA/abseta_pt_DATA");
    TH2D* sfmuhlt4= (TH2D*) muhlt1->Get("Mu50_OR_TkMu50_PtEtaBins/efficienciesMC/abseta_pt_MC");

    TFile * methlt = new TFile ("myTriggerScaleFactors_2016.root");
    TH2D* sfmethlt= (TH2D*) methlt->Get("MET_MU_DATA");
    TH2D* sfmethlt1= (TH2D*) methlt->Get("MET_MU_MC");

    TFile * muid = new TFile ("RunGH_SF_ID.root");
    TH2D* sfmuid= (TH2D*) muid->Get("NUM_HighPtID_DEN_genTracks_eta_pair_newTuneP_probe_pt");

    TFile * muid1 = new TFile ("RunBCDEF_SF_ID-2.root");
    TH2D* sfmuid1= (TH2D*) muid1->Get("NUM_HighPtID_DEN_genTracks_eta_pair_newTuneP_probe_pt");

    TFile * muiso = new TFile ("RunBCDEF_SF_ISO-2.root");
    TH2D* sfmuiso= (TH2D*) muiso->Get("NUM_LooseRelTkIso_DEN_HighPtIDandIPCut_eta_pair_newTuneP_probe_pt");

    TFile * muiso1 = new TFile ("RunGH_SF_ISO.root");
    TH2D* sfmuiso1= (TH2D*) muiso1->Get("NUM_LooseRelTkIso_DEN_HighPtIDandIPCut_eta_pair_newTuneP_probe_pt");

//double dnn_corr_t60[20] = {0.995885951944,0.983347452525,0.975759095113,0.953894714958,0.959454451025,0.960101121283,1.00515144059,1.00780247321,1.00122648613,0.999014002774,1.02867582857,1.03145755377,1.01231623338,1.0179995388,1.11338432762,1.12434265093,1.20040172725,1.23214781346,1.10853261871,1.18987881019};
//double dnn_corr_w60[20] = {2.39458356376,2.10675676757,1.96389723144,1.83305260353,1.62227004119,1.35159175913,1.28347995393,1.12865755511,1.06173304965,1.11626583343,1.10663444106,1.05993269819,1.03361665259,0.973676940489,0.954494833793,0.963835929878,0.854030933747,0.751693936006,0.727429712752,0.550473104372};
//double dnn_corr_q60[20] = {0.572375796481,0.708645299286,0.77648989791,0.840909818904,0.938448382313,0.971713004382,1.03147428396,1.05720560765,1.15566459444,1.14361862677,1.14984326299,1.16409295556,1.23784332694,1.26328103249,1.25892474891,1.20810615469,1.25082136728,1.31198884853,1.21616014376,1.4531834221};
//double dnn60[21] = {0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0};
//double dnn_corr_t100[15] = {0.847154349946,0.861446597227,0.875788794595,0.903514289215,0.935834360913,0.912597817545,0.916403678881,0.989538106113,0.998313251219,1.0488928671,1.07979175572,1.13121365835,1.07244126924,1.05710685554,0.9882153811};
//double dnn_corr_w100[15] = {3.28680454623,2.92532785121,2.59195693247,2.19822178861,1.8928291334,1.69629303853,1.84852637016,1.11132708885,0.717401303109,0.654081382093,0.624780735049,0.657329768143,0.708462529573,0.676771377827,0.459211205967};
//double dnn_corr_q100[15] = {0.708805080084,0.787588888102,0.859639275672,0.939090503188,0.995054549756,1.05411856889,1.01579702529,1.15344721429,1.19651054091,1.19572801079,1.1935752368,1.236249576,1.23675216132,1.20064677386,1.26795999152};
//double dnn100[16] = {0.0,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.85,0.9,0.95,1.0};
//double dnn_corr_t150[14] = {1.05488971447,1.03686633172,1.02207216476,1.00909077151,1.05208103692,1.05088622997,1.05614235922,1.04932323179,1.01335667257,1.01843056033,1.01381950436,1.01904926622,0.895666339885,0.86830177928};
//double dnn_corr_q150[14] = {0.734651023281,0.878344946553,0.996293498269,1.09978945966,1.09283498655,1.10644838559,0.992094502248,1.07406306746,1.09594019193,1.02955610664,1.03473230636,1.02886157044,1.16736664901,1.19808508574};
//double dnn150[15] = {0.0,0.05,0.1,0.15,0.2,0.25,0.4,0.55,0.6,0.75,0.8,0.85,0.9,0.95,1.0};
//double dnn_corr_t200[4] = {0.892895478407,0.868785588271,0.927177903317,0.996584114548};
//double dnn_corr_q200[4] = {0.911922472216,0.989819137309,1.24754288345,1.21823830571};
//double dnn200[5] = {0.0,0.25,0.5,0.75,1.0};

double dnn_corr_t60[20] = {1.17632061448,1.16495811349,1.01810830074,0.966049196558,0.92420465487,0.911264861325,0.988900385026,0.9813525651,0.992416197462,0.983912990526,0.965200560156,0.958772731415,0.962461915697,0.967269065268,1.11356707886,1.02358424442,1.37415244161,1.21025178237,1.26051558367,1.41164967254};
double dnn_corr_w60[20] = {1.85821858885,1.34905683327,0.908316834659,0.699749783686,0.329113071317,0.987950896932,0.835575087979,1.03205969712,0.891734031096,0.990452521359,0.919696889432,0.983140795981,0.841388664052,0.812643985365,0.698568003698,0.757988210937,0.651463495899,0.737621065734,0.718820005326,0.489286047378};
double dnn_corr_q60[20] = {0.394536303807,0.578529365642,0.753942891599,0.834897792016,0.972789065137,0.924580619938,1.02973873497,1.00997355902,1.18149830842,1.16898018904,1.26564063313,1.25471537646,1.48257454199,1.48976864361,1.63419529423,1.76031648635,1.72778472568,1.96525626849,2.10094382271,2.24586230036};
double dnn60[21] = {0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0};
double dnn_corr_t100[15] = {1.0606607231,1.0826085179,1.17347785441,1.01982391517,1.06953541688,1.01386610667,0.849017039436,1.02662222472,0.899843102587,0.786123810158,0.894546804748,1.02779348751,0.8754927327,0.640298731453,0.639906664811};
double dnn_corr_w100[15] = {1.77421989463,1.35924770238,1.12356916013,0.392676033627,0.284038409702,0.246824657173,0.233686248196,-0.590573717697,0.48150504793,0.739377051037,0.654824294156,0.92610489067,1.05063080233,1.47795790806,1.31692510236};
double dnn_corr_q100[15] = {0.687266918148,0.769754709091,0.78448292235,1.00855273501,1.01213629046,1.04314586015,1.1135755598,1.22252149154,1.21300091816,1.31665008822,1.31161309561,1.23505206691,1.23689640143,1.2559544879,1.27809530557};
double dnn100[16] = {0.0,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.85,0.9,0.95,1.0};
double dnn_corr_t150[14] = {0.939347855856,0.938784626029,0.917074494858,0.914195871734,1.02700373303,1.00642665098,1.09494793887,1.06643657764,0.941273856315,0.850004598804,0.846970304408,0.998947860073,0.921089010587,0.982692873036};
double dnn_corr_q150[14] = {0.884607323807,0.887239712064,0.988707162068,1.00216109036,0.987035366601,1.12873509025,0.820557825051,1.01399365087,1.06037711928,1.64653085285,1.64826171567,1.56156864401,1.60598193115,1.5708410313};
double dnn150[15] = {0.0,0.05,0.1,0.15,0.2,0.25,0.4,0.55,0.6,0.75,0.8,0.85,0.9,0.95,1.0};
double dnn_corr_t200[4] = {1.10171188733,1.22285375206,0.539648903836,1.13016576977};
double dnn_corr_q200[4] = {1.03701429685,0.784239673081,1.26670788311,1.09459657083};
double dnn200[5] = {0.0,0.25,0.5,0.75,1.0};


	int bin60=sizeof(dnn_corr_t60)/sizeof(dnn_corr_t60[0]);
	int bin100=sizeof(dnn_corr_t100)/sizeof(dnn_corr_t100[0]);
	int bin150=sizeof(dnn_corr_t150)/sizeof(dnn_corr_t150[0]);
	int bin200=sizeof(dnn_corr_t200)/sizeof(dnn_corr_t200[0]);
        Double_t massdnn[5]={60,100,150,200,1000000};

	TH1D *h_t60=new TH1D("h_t60","",bin60,0,1);
	TH1D *h_q60=new TH1D("h_q60","",bin60,0,1);
	TH1D *h_w60=new TH1D("h_w60","",bin60,0,1);
	TH1D *h_t100=new TH1D("h_t100","",bin100,0,1);
	TH1D *h_q100=new TH1D("h_q100","",bin100,0,1);
	TH1D *h_w100=new TH1D("h_w100","",bin100,0,1);
	TH1D *h_t150=new TH1D("h_t150","",bin150,0,1);
	TH1D *h_q150=new TH1D("h_q150","",bin150,0,1);
	TH1D *h_t200=new TH1D("h_t200","",bin200,0,1);
	TH1D *h_q200=new TH1D("h_q200","",bin200,0,1);
	h_t60->SetBins(bin60,dnn60);
	h_q60->SetBins(bin60,dnn60);
	h_w60->SetBins(bin60,dnn60);
	h_t100->SetBins(bin100,dnn100);
	h_q100->SetBins(bin100,dnn100);
	h_w100->SetBins(bin100,dnn100);
	h_t150->SetBins(bin150,dnn150);
	h_q150->SetBins(bin150,dnn150);
	h_t200->SetBins(bin200,dnn200);
	h_q200->SetBins(bin200,dnn200);
	for(int i=0;i<bin60;i++){
	h_t60->SetBinContent(i+1,dnn_corr_t60[i]);
	h_q60->SetBinContent(i+1,dnn_corr_q60[i]);
	h_w60->SetBinContent(i+1,dnn_corr_w60[i]);
	}
	for(int i=0;i<bin100;i++){
	h_t100->SetBinContent(i+1,dnn_corr_t100[i]);
	h_q100->SetBinContent(i+1,dnn_corr_q100[i]);
	h_w100->SetBinContent(i+1,dnn_corr_w100[i]);
	}
	for(int i=0;i<bin150;i++){
	h_t150->SetBinContent(i+1,dnn_corr_t150[i]);
	h_q150->SetBinContent(i+1,dnn_corr_q150[i]);
	}
	for(int i=0;i<bin200;i++){
	h_t200->SetBinContent(i+1,dnn_corr_t200[i]);
	h_q200->SetBinContent(i+1,dnn_corr_q200[i]);
	}

	Long64_t nbytes = 0, nb = 0;
	//for (Long64_t jentry=0; jentry<10;jentry++)
	for (Long64_t jentry=0; jentry<nentries;jentry++)
	{
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
                   if (jentry%1000000==0)
                      {std::cout<<jentry<<std::endl;}
		nb = fChain->GetEntry(jentry);
        //cout<<ientry<<"    "<<theWeight<<endl;
		nbytes += nb;
		pfMET             = Double_t(met);
		pfMETPhi          = Double_t(metPhi);
		l_pt              = Double_t(ptlep1);
        l_pt_2              = Double_t(ptlep2);

		l_eta             = Double_t(etalep1);
		l_phi             = Double_t(philep1);
                jet_mass_puppi    = Double_t(jetAK8puppi_sd);
                jet_mass_puppi_un = Double_t(jetAK8puppi_sd);
                jet_pt_puppi      = Double_t(jetAK8puppi_ptJEC);
//		jet_mass_softdrop = Double_t(sdropJEC);
		W_pt              = Double_t(ptVlepJEC);
		W_eta             = Double_t(yVlep);
		W_phi             = Double_t(phiVlep);
                m_JJlv             = Double_t(candMasspuppiJEC);
        m_Jlv_new=-99;
        m_JJlv_new=-99;
        m_Jlv_JEC_up=-99;
        m_Jlv_JEC_down=-99;
        m_Jlv_JER_up=-99;
        m_Jlv_JER_down=-99;
        m_JJlv_JEC_up=-99;
        m_JJlv_JEC_down=-99;
        m_JJlv_JER_up=-99;
        m_JJlv_JER_down=-99;
        triggerWeight=1;
        if (candMasspuppiJEC_new>0) m_JJlv_new             = Double_t(candMasspuppiJEC_new);
        if (candMasspuppiJEC_JEC_up>0) m_JJlv_JEC_up             = Double_t(candMasspuppiJEC_JEC_up);
        if (candMasspuppiJEC_JEC_down>0) m_JJlv_JEC_down             = Double_t(candMasspuppiJEC_JEC_down);
        if (candMasspuppiJEC_JER_up>0) m_JJlv_JER_up             = Double_t(candMasspuppiJEC_JER_up);
        if (candMasspuppiJEC_JER_down>0) m_JJlv_JER_down             = Double_t(candMasspuppiJEC_JER_down);

        m_Jlv = m_jlv;
        if (m_jlv_new>0) m_Jlv_new = m_jlv_new;
        if (m_jlv_JEC_up>0) m_Jlv_JEC_up = m_jlv_JEC_up;
        if (m_jlv_JEC_down>0) m_Jlv_JEC_down = m_jlv_JEC_down;
        if (m_jlv_JER_up>0) m_Jlv_JER_up = m_jlv_JER_up;
        if (m_jlv_JER_down>0) m_Jlv_JER_down = m_jlv_JER_down;

		jet_mass_puppi_2    = Double_t(jetAK8puppi_sd_2);
                jet_mass_puppi_un_2 = Double_t(jetAK8puppi_sd_2);
                jet_pt_puppi_2      = Double_t(jetAK8puppi_ptJEC_2);
        
        jet_mass_puppi_3    = Double_t(jetAK8puppi_sd_3);
        jet_mass_puppi_un_3 = Double_t(jetAK8puppi_sd_3);
        jet_pt_puppi_3      = Double_t(jetAK8puppi_ptJEC_3);

        jetAK8puppi_tau31=jetAK8puppi_tau3/jetAK8puppi_tau1;
        jetAK8puppi_tau32=jetAK8puppi_tau3/jetAK8puppi_tau2;
        jetAK8puppi_tau43=jetAK8puppi_tau4/jetAK8puppi_tau3;
        jetAK8puppi_tau41=jetAK8puppi_tau4/jetAK8puppi_tau1;
        jetAK8puppi_tau42=jetAK8puppi_tau4/jetAK8puppi_tau2;

        jetAK8puppi_tau31_2=jetAK8puppi_tau3_2/jetAK8puppi_tau1_2;
        jetAK8puppi_tau32_2=jetAK8puppi_tau3_2/jetAK8puppi_tau2_2;
        jetAK8puppi_tau43_2=jetAK8puppi_tau4_2/jetAK8puppi_tau3_2;
        jetAK8puppi_tau41_2=jetAK8puppi_tau4_2/jetAK8puppi_tau1_2;
        jetAK8puppi_tau42_2=jetAK8puppi_tau4_2/jetAK8puppi_tau2_2;

        jetAK8puppi_tau31_3=jetAK8puppi_tau3_3/jetAK8puppi_tau1_3;
        jetAK8puppi_tau32_3=jetAK8puppi_tau3_3/jetAK8puppi_tau2_3;
        jetAK8puppi_tau43_3=jetAK8puppi_tau4_3/jetAK8puppi_tau3_3;
        jetAK8puppi_tau41_3=jetAK8puppi_tau4_3/jetAK8puppi_tau1_3;
        jetAK8puppi_tau42_3=jetAK8puppi_tau4_3/jetAK8puppi_tau2_3;
 
        for(Int_t iii=0;iii<3;iii++)
MassVV[iii] =massww[iii];
//cout<<jet_mass_puppi_un_2<<endl;
//cout<<jetAK8puppi_sd_2<<endl;
//		
//		fjet2_pt          = Double_t(jet2_pt);
//		fjet2_btag        = Double_t(jet2_btag);
//		fjet3_pt          = Double_t(jet3_pt);
//		fjet3_btag        = Double_t(jet3_btag);
                mtVlepnew         = Double_t(sqrt(2*ptlep1*met*(1.0-cos(philep1-metPhi))));
                MTVlep            = Double_t(sqrt(2*ptlep1*MET_et*(1.0-cos(philep1-MET_phi))));
                //puppi+softdrop recorrected by Thea's "JEC"
                Double_t gencorrect=1.0;
                Double_t recocorrect_0eta1p3=1.0;
                Double_t recocorrect_1p3eta2p5=1.0;
        jet_mass_puppi_corr=-99;jet_mass_puppi_corr_2=-99;jet_mass_puppi_corr_3=-99;
//                gencorrect=1.0-0.321*pow(jet_pt_puppi*0.0354,-1.1);
//                recocorrect_0eta1p3=1.09-1.69e-04*jet_pt_puppi+3.34e-07*pow(jet_pt_puppi,2)-2.47e-10*pow(jet_pt_puppi,3)+7.8e-14*pow(jet_pt_puppi,4)-8.83e-18*pow(jet_pt_puppi,5);
//                recocorrect_1p3eta2p5=1.3-7.76e-04*jet_pt_puppi+1.11e-06*pow(jet_pt_puppi,2)-6.79e-10*pow(jet_pt_puppi,3)+1.87e-13*pow(jet_pt_puppi,4)-1.9e-17*pow(jet_pt_puppi,5);
                gencorrect=1.006-1.062*pow(jet_pt_puppi*0.08,-1.2);
                recocorrect_0eta1p3=1.093-1.501e-04*jet_pt_puppi+3.449e-07*pow(jet_pt_puppi,2)-2.681e-10*pow(jet_pt_puppi,3)+8.674e-14*pow(jet_pt_puppi,4)-1.001e-17*pow(jet_pt_puppi,5);
                recocorrect_1p3eta2p5=1.272-5.72e-04*jet_pt_puppi+8.37e-07*pow(jet_pt_puppi,2)-5.204e-10*pow(jet_pt_puppi,3)+1.454e-13*pow(jet_pt_puppi,4)-1.504e-17*pow(jet_pt_puppi,5);
                if (fabs(jetAK8puppi_eta)<=1.3){jet_mass_puppi_corr=jet_mass_puppi_un*gencorrect*recocorrect_0eta1p3;}
                else if (fabs(jetAK8puppi_eta)<2.5 && fabs(jetAK8puppi_eta)>1.3){jet_mass_puppi_corr=jet_mass_puppi_un*gencorrect*recocorrect_1p3eta2p5;}


Double_t gencorrect_2=1.0;
                Double_t recocorrect_0eta1p3_2=1.0;
                Double_t recocorrect_1p3eta2p5_2=1.0;
gencorrect_2=1.006-1.062*pow(jet_pt_puppi_2*0.08,-1.2);
                recocorrect_0eta1p3_2=1.093-1.501e-04*jet_pt_puppi_2+3.449e-07*pow(jet_pt_puppi_2,2)-2.681e-10*pow(jet_pt_puppi_2,3)+8.674e-14*pow(jet_pt_puppi_2,4)-1.001e-17*pow(jet_pt_puppi_2,5);
                recocorrect_1p3eta2p5_2=1.272-5.72e-04*jet_pt_puppi_2+8.37e-07*pow(jet_pt_puppi_2,2)-5.204e-10*pow(jet_pt_puppi_2,3)+1.454e-13*pow(jet_pt_puppi_2,4)-1.504e-17*pow(jet_pt_puppi_2,5);
                if (fabs(jetAK8puppi_eta_2)<=1.3){jet_mass_puppi_corr_2=jet_mass_puppi_un_2*gencorrect_2*recocorrect_0eta1p3_2;}
                else if (fabs(jetAK8puppi_eta_2)<2.5 && fabs(jetAK8puppi_eta_2)>1.3){jet_mass_puppi_corr_2=jet_mass_puppi_un_2*gencorrect_2*recocorrect_1p3eta2p5_2;}

        Double_t gencorrect_3=1.0;
        Double_t recocorrect_0eta1p3_3=1.0;
        Double_t recocorrect_1p3eta2p5_3=1.0;
        gencorrect_3=1.006-1.062*pow(jet_pt_puppi_3*0.08,-1.2);
        recocorrect_0eta1p3_3=1.093-1.501e-04*jet_pt_puppi_3+3.449e-07*pow(jet_pt_puppi_3,2)-2.681e-10*pow(jet_pt_puppi_3,3)+8.674e-14*pow(jet_pt_puppi_3,4)-1.001e-17*pow(jet_pt_puppi_3,5);
        recocorrect_1p3eta2p5_3=1.272-5.72e-04*jet_pt_puppi_3+8.37e-07*pow(jet_pt_puppi_3,2)-5.204e-10*pow(jet_pt_puppi_3,3)+1.454e-13*pow(jet_pt_puppi_3,4)-1.504e-17*pow(jet_pt_puppi_3,5);
        if (fabs(jetAK8puppi_eta_3)<=1.3){jet_mass_puppi_corr_3=jet_mass_puppi_un_3*gencorrect_3*recocorrect_0eta1p3_3;}
        else if (fabs(jetAK8puppi_eta_3)<2.5 && fabs(jetAK8puppi_eta_3)>1.3){jet_mass_puppi_corr_3=jet_mass_puppi_un_3*gencorrect_3*recocorrect_1p3eta2p5_3;}
        jet_mass_puppi_corr=jetAK8puppi_sd;
        jet_mass_puppi_corr_2=jetAK8puppi_sd_2;
        jet_mass_puppi_corr_3=jetAK8puppi_sd_3;
		// GEN-RECO match
		//deltaRleplep = deltaR(etalep1,philep1,etalep2,philep2);
		//deltaRWlepGen = deltaR(etaGenVlep, phiGenVlep, yVlep, phiVlep);
		//jetV.SetPtEtaPhiM(ptVhad, yVhad, phiVhad, massVhad);
		//genjetV.SetPtEtaPhiM(ptGenVhad, etaGenVhad, phiGenVhad, massGenVhad);
		//deltaRWhadGen = deltaR(etaGenVhad, phiGenVhad, yVhad, phiVhad);
		//Double_t deltaRWhadGen = sqrt(pow(etaGenVhad-yVhad,2) + pow(phiGenVhad-phiVhad,2));
		Double_t deltaRWhadGen = sqrt(pow(etaGenVhad-jetAK8puppi_eta,2) + pow(phiGenVhad-jetAK8puppi_phi,2));

		//Weight Calculation
        if(IsData>0 && npT<0) continue;
		Int_t bin = hR1->FindBin(npT);
		pileupWeight = hR1->GetBinContent(bin);		

		eff_and_pu_Weight = 0;
		eff_and_pu_Weight1 = 0;
		if(IsData>0) {
			if(npT < weights_pu1.size()){
				eff_and_pu_Weight = weights_pu1[npT];
			}
			if(npT < weights_pu2.size()){
				eff_and_pu_Weight1 = weights_pu2[npT];
			}
		}
/*
      trigger_eff=1.0;

      TFile * input_trigger = new TFile ("MuonEff-TH2D.root");
      input_trigger->cd("");
      TH2D* HLTeff= (TH2D*) input_trigger->Get("2Dh");

      Double_t mu_pt=ptlep1;
      if (mu_pt>500) {mu_pt=499.0;}
      int mubin=HLTeff->FindBin(fabs(etalep1),mu_pt);
      trigger_eff=HLTeff->GetBinContent(mubin); 
      input_trigger->Close();

*/
//      std::cout<<trigger_eff<<std::endl;
        if(IsData>0){
    if(channelname=="el"){
        trigger_eff=1.0;
        Double_t El_MC=1.0;
        Double_t El_Data=1.0;
        Double_t MET_MC=1.0;
        Double_t MET_Data=1.0;
        //https://indico.cern.ch/event/787353/contributions/3417731/attachments/1842466/3127748/190821__JaesungKim__EGammaHLTMeeting__EGammaHLTEfficiency_for_WRAnalysis_update.pdf
        if(ptlep1>=sfelehltbar->GetYaxis()->GetBinLowEdge(1)&&ptlep1<=sfelehltbar->GetYaxis()->GetBinUpEdge(sfelehltbar->GetNbinsY())){
        if (fabs(etalep1)<1.444) {
            int elebinhltbar=sfelehltbar->FindBin(etalep1,ptlep1);
            trigger_eff=sfelehltbar->GetBinContent(elebinhltbar);}
        if (fabs(etalep1)>1.566) {
            int elebinhltend=sfelehltend->FindBin(etalep1,ptlep1);
            trigger_eff=sfelehltend->GetBinContent(elebinhltend);}}

        IDweight=1.0;
        /*
         if (ptlep1<20.)
         {if (fabs(etalep1)<0.8){IDweight=0.992611;}
         else if (fabs(etalep1)<1.44){IDweight=1.11702;}
         else if (fabs(etalep1)<1.56){IDweight=0.825893;}
         else if (fabs(etalep1)<2.0){IDweight=0.486322;}
         else if (fabs(etalep1)>=2.0){IDweight=0.795031;}
         }
         else if (ptlep1<30.)
         {if (fabs(etalep1)<0.8){IDweight=0.970909;}
         else if (fabs(etalep1)<1.44){IDweight=0.949343;}
         else if (fabs(etalep1)<1.56){IDweight=0.812865;}
         else if (fabs(etalep1)<2.0){IDweight=0.604607;}
         else if (fabs(etalep1)>=2.0){IDweight=0.794971;}
         }
         else if (ptlep1<40.)
         {if (fabs(etalep1)<0.8){IDweight=0.975867;}
         else if (fabs(etalep1)<1.44){IDweight=0.975232;}
         else if (fabs(etalep1)<1.56){IDweight=0.819502;}
         else if (fabs(etalep1)<2.0){IDweight=0.694099;}
         else if (fabs(etalep1)>=2.0){IDweight=0.853312;}
         }
         else if (ptlep1<50.)
         {if (fabs(etalep1)<0.8){IDweight=0.974359;}
         else if (fabs(etalep1)<1.44){IDweight=0.972603;}
         else if (fabs(etalep1)<1.56){IDweight=0.819805;}
         else if (fabs(etalep1)<2.0){IDweight=0.740027;}
         else if (fabs(etalep1)>=2.0){IDweight=0.861777;}
         }
         */
        //https://indico.cern.ch/event/831669/contributions/3485543/attachments/1871797/3084930/ApprovalSlides_EE_v3.pdf
        //https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRunIIRecommendations
        if (ptlep1>30.)
        {if (fabs(etalep1)<1.444){IDweight=1*0.983;}
        else if (fabs(etalep1)>=1.566){IDweight=1*0.991;}
        }
        
        IDweightISO=1.0;
        
        IDweighttrk=1.0;//Electron Reconstruction Scale Factor
        if(ptlep1>=sfelereco->GetYaxis()->GetBinLowEdge(1)&&ptlep1<=sfelereco->GetYaxis()->GetBinUpEdge(sfelereco->GetNbinsY())){
        Double_t ele_pt=ptlep1;
        if (ele_pt>500) {ele_pt=499.0;}
        int elebin=sfelereco->FindBin(etalep1,ele_pt);
            IDweighttrk=sfelereco->GetBinContent(elebin);}

        /*if (etalep1<-2.45){IDweighttrk=1.317604;}
        else if(etalep1<-2.4){IDweighttrk=1.11378;}
        else if(etalep1<-2.3){IDweighttrk=1.024625;}
        else if(etalep1<-2.2){IDweighttrk=1.013641;}
        else if(etalep1<-2.){IDweighttrk=1.007277;}
        else if(etalep1<-1.8){IDweighttrk=0.994817;}
        else if(etalep1<-1.63){IDweighttrk=0.994786;}
        else if(etalep1<-1.566){IDweighttrk=0.991632;}
        else if(etalep1<-1.444){IDweighttrk=0.963129;}
        else if(etalep1<-1.2){IDweighttrk=0.989701;}
        else if(etalep1<-1.){IDweighttrk=0.985656;}
        else if(etalep1<-0.6){IDweighttrk=0.981595;}
        else if(etalep1<-0.4){IDweighttrk=0.984678;}
        else if(etalep1<-0.2){IDweighttrk=0.981614;}
        else if(etalep1<0.0){IDweighttrk=0.980433;}
        else if(etalep1<0.2){IDweighttrk=0.984552;}
        else if(etalep1<0.4){IDweighttrk=0.98876;}
        else if(etalep1<0.6){IDweighttrk=0.987743;}
        else if(etalep1<1.0){IDweighttrk=0.987743;}
        else if(etalep1<1.2){IDweighttrk=0.987743;}
        else if(etalep1<1.444){IDweighttrk=0.987679;}
        else if(etalep1<1.566){IDweighttrk=0.967598;}
        else if(etalep1<1.63){IDweighttrk=0.989627;}
        else if(etalep1<1.8){IDweighttrk=0.992761;}
        else if(etalep1<2.0){IDweighttrk=0.991761;}
        else if(etalep1<2.2){IDweighttrk=0.99794;}
        else if(etalep1<2.3){IDweighttrk=1.001037;}
        else if(etalep1<2.4){IDweighttrk=0.989507;}
        else if(etalep1<2.45){IDweighttrk=0.970519;}
        else if(etalep1<2.5){IDweighttrk=0.906667;}*/
        
        
        /*
         TFile * input_ID = new TFile ("elesf.root");
         input_ID->cd("");
         TH2D* hele= (TH2D*) input_ID->Get("elesf");
         
         Double_t elept=ptlep1;
         if (elept>200) {elept=199.0;}
         int elebin0=hele->FindBin(etalep1,elept);
         IDweight=hele->GetBinContent(elebin0);
         input_ID->Close();
         */
        //cout << "pileupWeight:"<<pileupWeight<< " eff_and_pu_Weight:" << eff_and_pu_Weight << " eff_and_pu_Weight1:" << eff_and_pu_Weight1 << endl;
        if(theWeight>0) nn=1;
        else nn= -1;
        if(npp>0) lumiWeight=Identical_lumiWeight/(npp-nmm);
        else lumiWeight=Identical_lumiWeight/nentries;
        weight_nobtag=lumiWeight*triggerWeight*eff_and_pu_Weight*nn*trigger_eff;
        //weight=lumiWeight*triggerWeight*pileupWeight*nn;
        if (IsData>1 &&IsData<9 ) weight_nobtag = weight_nobtag*1.21;
        
        //lumiWeight=Identical_lumiWeight;
        //if(npp>0) weight=lumiWeight*triggerWeight*pileupWeight/(npp-nmm)*nn*0.04024;//0.00559;
        //else weight=lumiWeight*triggerWeight*pileupWeight/nentries*0.04024;//0.00559;
        if ( IsData==0 ) weight_nobtag=1;
        //Weight Calculation Done
        

    }
        if(channelname=="mu"){
            trigger_eff=1.0;
            double effmud1=1,effmud2=1,effmetd=1,effmumc1=1,effmumc2=1,effmetmc=1;
            //https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults
            if(ptlep1>=sfmuhlt->GetYaxis()->GetBinLowEdge(1)&&ptlep1<=sfmuhlt->GetYaxis()->GetBinUpEdge(sfmuhlt->GetNbinsY())){
            double tmp_eta=fabs(etalep1);
            if(tmp_eta>2.4)tmp_eta=2.39;
            double tmp_pthlt=ptlep1;
            if(tmp_pthlt>1200)tmp_pthlt=1199;
            int mubinhlt=sfmuhlt->FindBin(tmp_eta,tmp_pthlt);
                effmud1=sfmuhlt->GetBinContent(mubinhlt);
                effmumc1=sfmuhlt2->GetBinContent(mubinhlt);
                effmud2=sfmuhlt3->GetBinContent(mubinhlt);
                effmumc2=sfmuhlt4->GetBinContent(mubinhlt);
            }
            if(MET_et>=sfmethlt->GetXaxis()->GetBinLowEdge(1)&&MET_et<=sfmethlt->GetXaxis()->GetBinUpEdge(sfmethlt->GetNbinsX())){
                effmetd=sfmethlt->GetBinContent(sfmethlt->FindBin(MET_et));
                effmetmc=sfmethlt1->GetBinContent(sfmethlt1->FindBin(MET_et));
            }
            double lumibf=5.746010293+2.572903489+4.242291557+4.025228049+3.104509119;
            double lumigh=7.575579377+8.650628335;
            trigger_eff=((lumibf*effmud1+lumigh*effmud2)/(lumibf+lumigh)+effmetd-(lumibf*effmud1+lumigh*effmud2)/(lumibf+lumigh)*effmetd)/((lumibf*effmumc1+lumigh*effmumc2)/(lumibf+lumigh)+effmetmc-(lumibf*effmumc1+lumigh*effmumc2)/(lumibf+lumigh)*effmetmc);
            //cout<<effmud1<<"   "<<effmumc1<<"   "<<effmud2<<"   "<<effmumc2<<"   "<<effmetd<<"   "<<effmetmc<<"   "<<trigger_eff<<"   "<<endl;
        
		if(theWeight>0) nn=1;
		else nn= -1;
		if(npp>0) lumiWeight=Identical_lumiWeight/(npp-nmm);
		else lumiWeight=Identical_lumiWeight/nentries;
		weight_nobtag=lumiWeight*triggerWeight*eff_and_pu_Weight*nn*trigger_eff;
		//weight=lumiWeight*triggerWeight*pileupWeight*nn;
		if (IsData>1 &&IsData<9 ) weight_nobtag = weight_nobtag*1.21;

		//lumiWeight=Identical_lumiWeight;
		//if(npp>0) weight=lumiWeight*triggerWeight*pileupWeight/(npp-nmm)*nn*0.04024;//0.00559;
		//else weight=lumiWeight*triggerWeight*pileupWeight/nentries*0.04024;//0.00559;
		if ( IsData==0 ) weight_nobtag=1;
		//Weight Calculation Done


      IDweight=1.0;
            if(ptlep1>=sfmuid->GetYaxis()->GetBinLowEdge(1)&&ptlep1<=sfmuid->GetYaxis()->GetBinUpEdge(sfmuid->GetNbinsY())){
            double tmp_ptlep1=ptlep1;
            if(tmp_ptlep1>120)tmp_ptlep1=119.;
                double tmp_etaid=etalep1;
                if(tmp_etaid>2.4)tmp_etaid=2.39;
            int mubinid=sfmuid->FindBin(tmp_etaid,tmp_ptlep1);
                IDweight=sfmuid->GetBinContent(mubinid);}

                /*if (ptlep1>=55 && ptlep1<60)
                   {if (fabs(etalep1)<0.9){IDweight=0.982688;}
                    else if (fabs(etalep1)<1.2){IDweight=0.976116;}
                    else if (fabs(etalep1)<2.1){IDweight=0.986235;}
                    else if (fabs(etalep1)<2.4){IDweight=0.966154;}
                   }
                else if (ptlep1>=60)
                   {if (fabs(etalep1)<0.9){IDweight=0.99385;}
                    else if (fabs(etalep1)<1.2){IDweight=0.976965;}
                    else if (fabs(etalep1)<2.1){IDweight=0.990001;}
                    else if (fabs(etalep1)<2.4){IDweight=0.967312;}
                   }*/

      IDweightISO=1.0;
            if(ptlep1>=sfmuiso->GetYaxis()->GetBinLowEdge(1)&&ptlep1<=sfmuiso->GetYaxis()->GetBinUpEdge(sfmuiso->GetNbinsY())){
            double tmp_ptlep1=ptlep1;
            if(tmp_ptlep1>120)tmp_ptlep1=119.;
                double tmp_etaiso=etalep1;
                if(tmp_etaiso>2.4)tmp_etaiso=2.39;

            int mubiniso=sfmuiso->FindBin(tmp_etaiso,tmp_ptlep1);
                IDweightISO=sfmuiso->GetBinContent(mubiniso);}

                /*if (ptlep1>=55 && ptlep1<60)
                   {if (fabs(etalep1)<0.9){IDweightISO=0.997968;}
                    else if (fabs(etalep1)<1.2){IDweightISO=0.998721;}
                    else if (fabs(etalep1)<2.1){IDweightISO=1.000081;}
                    else if (fabs(etalep1)<2.4){IDweightISO=1.00056;}
                   }
                else if (ptlep1>=60)
                   {if (fabs(etalep1)<0.9){IDweightISO=0.998726;}
                    else if (fabs(etalep1)<1.2){IDweightISO=0.999314;}
                    else if (fabs(etalep1)<2.1){IDweightISO=0.999228;}
                    else if (fabs(etalep1)<2.4){IDweightISO=1.002153;}
                   }*/

      IDweighttrk=1.0;
/*
        if (etalep1<-2.1){IDweighttrk=9.82399009186853522e-01;}
        else if(etalep1<-1.6){IDweighttrk=9.91746789037933008e-01;}
        else if(etalep1<-1.1){IDweighttrk=9.95944961092376846e-01;}
        else if(etalep1<-0.6){IDweighttrk=9.93413142541369476e-01;}
        else if(etalep1<0.0){IDweighttrk=9.91460688530866996e-01;}
        else if(etalep1<0.6){IDweighttrk=9.94680143661991423e-01;}
        else if(etalep1<1.1){IDweighttrk=9.96666389348924819e-01;}
        else if(etalep1<1.6){IDweighttrk=9.94933892427240618e-01;}
        else if(etalep1<2.1){IDweighttrk=9.91186607207322878e-01;}
        else if(etalep1<2.4){IDweighttrk=9.76811919457875155e-01;}
*/
        }
        ToppTweight=1.0;
        double a_top=0.0615;
        double b_top=- 0.0005;
        if (gentop_pt<0 || genantitop_pt<0){ToppTweight=1.0;}
        else if(gentop_pt>0 && genantitop_pt>0){ToppTweight=sqrt(exp(a_top+b_top*gentop_pt)*exp(a_top+b_top*genantitop_pt));}
        else if(gentop_pt>0 && genantitop_pt<0){ToppTweight=sqrt(exp(a_top+b_top*gentop_pt));}
        else if(gentop_pt<0 && genantitop_pt>0){ToppTweight=sqrt(exp(a_top+b_top*genantitop_pt));}

        }
        if(IsData==0){
            trigger_eff=1.0;
            if(theWeight>0) nn=1;
            else nn= -1;
            if(npp>0) lumiWeight=Identical_lumiWeight/(npp-nmm);
            else lumiWeight=Identical_lumiWeight/nentries;
            weight_nobtag=lumiWeight*triggerWeight*eff_and_pu_Weight*nn*trigger_eff;
            //weight=lumiWeight*triggerWeight*pileupWeight*nn;
            if (IsData>1 &&IsData<9 ) weight_nobtag = weight_nobtag*1.21;
            if ( IsData==0 ) weight_nobtag=1;
            IDweight=1.0;
            IDweightISO=1.0;
            IDweighttrk=1.0;
            ToppTweight=1.0;
        }


        
//--BSF----------------------------
      btagweight_center=1.0;
      btagweight_up=1.0;
      btagweight_down=1.0;
        btagweight_center_deep=1.0;
        btagweight_up_deep=1.0;
        btagweight_down_deep=1.0;
        btagweight_center_tig_deep=1.0;
        btagweight_up_tig_deep=1.0;
        btagweight_down_tig_deep=1.0;
        btagweight_center_los_deep=1.0;
        btagweight_up_los_deep=1.0;
        btagweight_down_los_deep=1.0;

        double deltaRAK4AK8_new_2[8]={0};
        Nj8_2=1,Nj8=0;
        
        if(jet_pt_puppi_2>0){
            Nj8_2=2;}
        
        Nj4=0;
        for(Int_t ii=0;ii<8;ii++){
            if(ak4jet_pt[ii]>0){
                Nj4++;
            }
        }
        for(Int_t ii=0;ii<4;ii++){
            if(jetAK8puppi_pt1[ii]>0&&abs(jetAK8puppi_eta1[ii])<2.4){
                Nj8++;
            }
        }
        if(jet_mass_puppi<-100||jet_mass_puppi_2<-100||jet_mass_puppi_3<-100) continue;

        if(IsData>0){
      double  bweight=1.0, dweight=1.0;
      double  bweightup=1.0;
      double  bweightdown=1.0;
            double  bweight_deep=1.0, dweight_deep=1.0;
            double  bweightup_deep=1.0;
            double  bweightdown_deep=1.0;
            double  bweight_tig_deep=1.0, dweight_tig_deep=1.0;
            double  bweightup_tig_deep=1.0;
            double  bweightdown_tig_deep=1.0;
            double  bweight_los_deep=1.0, dweight_los_deep=1.0;
            double  bweightup_los_deep=1.0;
            double  bweightdown_los_deep=1.0;

      double  beff=1.0, ceff=1.0, leff=1.0;

        
      for(Int_t i=0; i<8; i++)  {
       deltaRAK4AK8_new[i]=0.;
       deltaRAK4AK8_new[i]=sqrt(pow(fabs(ak4jet_eta[i]-jetAK8puppi_eta),2)+pow(TMath::Min(fabs(ak4jet_phi[i]-jetAK8puppi_phi),2*Pi-fabs(ak4jet_phi[i]-jetAK8puppi_phi)),2));
          if(Nj8==2) deltaRAK4AK8_new_2[i]=sqrt(pow(fabs(ak4jet_eta[i]-jetAK8puppi_eta_2),2)+pow(TMath::Min(fabs(ak4jet_phi[i]-jetAK8puppi_phi_2),2*Pi-fabs(ak4jet_phi[i]-jetAK8puppi_phi_2)),2));

       if(ak4jet_pt[i]>30 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8_new[i]>=0.8 && (Nj8==1||deltaRAK4AK8_new_2[i]>=0.8)) {
         if(abs(ak4jet_hf[i])==5 && ak4jet_pt[i]<3000.) 
          { 
            double jet_sf    = bsv(1,ak4jet_pt[i]); 
            double jet_sfu    = bsv(2,ak4jet_pt[i]);
            double jet_sfd    = bsv(3,ak4jet_pt[i]);
              double jet_sf_deep    = bsv_deep(1,ak4jet_pt[i]);
              double jet_sfu_deep    = bsv_deep(2,ak4jet_pt[i]);
              double jet_sfd_deep    = bsv_deep(3,ak4jet_pt[i]);
              double jet_sf_tig_deep    = bsv_tig_deep(1,ak4jet_pt[i]);
              double jet_sfu_tig_deep    = bsv_tig_deep(2,ak4jet_pt[i]);
              double jet_sfd_tig_deep    = bsv_tig_deep(3,ak4jet_pt[i]);
              double jet_sf_los_deep    = bsv_los_deep(1,ak4jet_pt[i]);
              double jet_sfu_los_deep    = bsv_los_deep(2,ak4jet_pt[i]);
              double jet_sfd_los_deep    = bsv_los_deep(3,ak4jet_pt[i]);

//            int bbin=hb->FindBin(ak4jet_pt[i],ak4jet_eta[i]);
//            beff=hb->GetBinContent(bbin);
              //b jet middle wp csv_2
              if (ak4jet_pt[i]>=30.0 && ak4jet_pt[i]<50.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.647332650923;}
              else if (ak4jet_eta[i]<-0.8){beff=0.574428968825;}
              else if (ak4jet_eta[i]>-0.8){beff=0.580733478893;}
              }
              else if (ak4jet_pt[i]>=50.0 && ak4jet_pt[i]<70.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.681764680183;}
              else if (ak4jet_eta[i]<-0.8){beff=0.607928437101;}
              else if (ak4jet_eta[i]>-0.8){beff=0.615291984235;}
              }
              else if (ak4jet_pt[i]>=70.0 && ak4jet_pt[i]<100.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.695375820567;}
              else if (ak4jet_eta[i]<-0.8){beff=0.622144210455;}
              else if (ak4jet_eta[i]>-0.8){beff=0.629384525268;}
              }
              else if (ak4jet_pt[i]>=100.0 && ak4jet_pt[i]<140.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.696089410422;}
              else if (ak4jet_eta[i]<-0.8){beff=0.619753472781;}
              else if (ak4jet_eta[i]>-0.8){beff=0.627676830787;}
              }
              else if (ak4jet_pt[i]>=140.0 && ak4jet_pt[i]<200.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.680199371866;}
              else if (ak4jet_eta[i]<-0.8){beff=0.611069273396;}
              else if (ak4jet_eta[i]>-0.8){beff=0.616045887084;}
              }
              else if (ak4jet_pt[i]>=200.0 && ak4jet_pt[i]<300.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.639681418349;}
              else if (ak4jet_eta[i]<-0.8){beff=0.563957042283;}
              else if (ak4jet_eta[i]>-0.8){beff=0.568304340723;}
              }
              else if (ak4jet_pt[i]>=300.0 && ak4jet_pt[i]<600.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.560957890044;}
              else if (ak4jet_eta[i]<-0.8){beff=0.523701843813;}
              else if (ak4jet_eta[i]>-0.8){beff=0.526011078231;}
              }
              else if (ak4jet_pt[i]>=600.0 && ak4jet_pt[i]<1000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.46878094322;}
              else if (ak4jet_eta[i]<-0.8){beff=0.47157152604;}
              else if (ak4jet_eta[i]>-0.8){beff=0.471656125556;}
              }
              else if (ak4jet_pt[i]>=1000.0 && ak4jet_pt[i]<3000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.41957313899;}
              else if (ak4jet_eta[i]<-0.8){beff=0.450354609929;}
              else if (ak4jet_eta[i]>-0.8){beff=0.425266903915;}
              }
              if(ak4jet_icsv[i]>0.8838) {
                  bweight=bweight*beff*jet_sf; dweight=dweight*beff;
                  bweightup=bweightup*beff*jet_sfu;
                  bweightdown=bweightdown*beff*jet_sfd;
              }
              else {
                  bweight=bweight*(1-beff*jet_sf); dweight=dweight*(1-beff);
                  bweightup=bweightup*(1-beff*jet_sfu);
                  bweightdown=bweightdown*(1-beff*jet_sfd);
              }
              //b jet middle wp DeepCSV
              if (ak4jet_pt[i]>=30.0 && ak4jet_pt[i]<50.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.69289054503;}
              else if (ak4jet_eta[i]<-0.8){beff=0.626876762927;}
              else if (ak4jet_eta[i]>-0.8){beff=0.637414613848;}
              }
              else if (ak4jet_pt[i]>=50.0 && ak4jet_pt[i]<70.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.726955160895;}
              else if (ak4jet_eta[i]<-0.8){beff=0.663040270851;}
              else if (ak4jet_eta[i]>-0.8){beff=0.672930844849;}
              }
              else if (ak4jet_pt[i]>=70.0 && ak4jet_pt[i]<100.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.739178943248;}
              else if (ak4jet_eta[i]<-0.8){beff=0.671524512238;}
              else if (ak4jet_eta[i]>-0.8){beff=0.681053968587;}
              }
              else if (ak4jet_pt[i]>=100.0 && ak4jet_pt[i]<140.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.742580432054;}
              else if (ak4jet_eta[i]<-0.8){beff=0.670239757334;}
              else if (ak4jet_eta[i]>-0.8){beff=0.677429556976;}
              }
              else if (ak4jet_pt[i]>=140.0 && ak4jet_pt[i]<200.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.735916249616;}
              else if (ak4jet_eta[i]<-0.8){beff=0.662461650933;}
              else if (ak4jet_eta[i]>-0.8){beff=0.667615935214;}
              }
              else if (ak4jet_pt[i]>=200.0 && ak4jet_pt[i]<300.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.706052668703;}
              else if (ak4jet_eta[i]<-0.8){beff=0.628799388706;}
              else if (ak4jet_eta[i]>-0.8){beff=0.631470724493;}
              }
              else if (ak4jet_pt[i]>=300.0 && ak4jet_pt[i]<600.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.652342361457;}
              else if (ak4jet_eta[i]<-0.8){beff=0.587623735846;}
              else if (ak4jet_eta[i]>-0.8){beff=0.588928672311;}
              }
              else if (ak4jet_pt[i]>=600.0 && ak4jet_pt[i]<1000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.561566159168;}
              else if (ak4jet_eta[i]<-0.8){beff=0.537073003542;}
              else if (ak4jet_eta[i]>-0.8){beff=0.532968340871;}
              }
              else if (ak4jet_pt[i]>=1000.0 && ak4jet_pt[i]<3000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.498954521694;}
              else if (ak4jet_eta[i]<-0.8){beff=0.479573712256;}
              else if (ak4jet_eta[i]>-0.8){beff=0.467857142857;}
              }
              if(ak4jet_icsv[i]>0.6321) {
                  bweight_deep=bweight_deep*beff*jet_sf_deep; dweight_deep=dweight_deep*beff;
                  bweightup_deep=bweightup_deep*beff*jet_sfu_deep;
                  bweightdown_deep=bweightdown_deep*beff*jet_sfd_deep;
              }
              else {
                  bweight_deep=bweight_deep*(1-beff*jet_sf_deep); dweight_deep=dweight_deep*(1-beff);
                  bweightup_deep=bweightup_deep*(1-beff*jet_sfu_deep);
                  bweightdown_deep=bweightdown_deep*(1-beff*jet_sfd_deep);
              }
              
              
              //b jet loose wp DeepCSV
              if (ak4jet_pt[i]>=30.0 && ak4jet_pt[i]<50.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.841210764847;}
              else if (ak4jet_eta[i]<-0.8){beff=0.810264618953;}
              else if (ak4jet_eta[i]>-0.8){beff=0.816790149527;}
              }
              else if (ak4jet_pt[i]>=50.0 && ak4jet_pt[i]<70.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.864728995719;}
              else if (ak4jet_eta[i]<-0.8){beff=0.834669348882;}
              else if (ak4jet_eta[i]>-0.8){beff=0.839980838846;}
              }
              else if (ak4jet_pt[i]>=70.0 && ak4jet_pt[i]<100.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.876020129011;}
              else if (ak4jet_eta[i]<-0.8){beff=0.844805459205;}
              else if (ak4jet_eta[i]>-0.8){beff=0.849756794514;}
              }
              else if (ak4jet_pt[i]>=100.0 && ak4jet_pt[i]<140.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.882550883747;}
              else if (ak4jet_eta[i]<-0.8){beff=0.851453339837;}
              else if (ak4jet_eta[i]>-0.8){beff=0.854276310286;}
              }
              else if (ak4jet_pt[i]>=140.0 && ak4jet_pt[i]<200.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.884411924905;}
              else if (ak4jet_eta[i]<-0.8){beff=0.857190397309;}
              else if (ak4jet_eta[i]>-0.8){beff=0.858229372025;}
              }
              else if (ak4jet_pt[i]>=200.0 && ak4jet_pt[i]<300.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.877699871512;}
              else if (ak4jet_eta[i]<-0.8){beff=0.853839246023;}
              else if (ak4jet_eta[i]>-0.8){beff=0.854308470831;}
              }
              else if (ak4jet_pt[i]>=300.0 && ak4jet_pt[i]<600.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.866537339179;}
              else if (ak4jet_eta[i]<-0.8){beff=0.855038623192;}
              else if (ak4jet_eta[i]>-0.8){beff=0.854764213243;}
              }
              else if (ak4jet_pt[i]>=600.0 && ak4jet_pt[i]<1000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.851237560602;}
              else if (ak4jet_eta[i]<-0.8){beff=0.866103050383;}
              else if (ak4jet_eta[i]>-0.8){beff=0.859676267555;}
              }
              else if (ak4jet_pt[i]>=1000.0 && ak4jet_pt[i]<3000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.861474124412;}
              else if (ak4jet_eta[i]<-0.8){beff=0.852575488455;}
              else if (ak4jet_eta[i]>-0.8){beff=0.857142857143;}
              }
              if(ak4jet_icsv[i]>0.2217) {
                  bweight_los_deep=bweight_los_deep*beff*jet_sf_los_deep; dweight_los_deep=dweight_los_deep*beff;
                  bweightup_los_deep=bweightup_los_deep*beff*jet_sfu_los_deep;
                  bweightdown_los_deep=bweightdown_los_deep*beff*jet_sfd_los_deep;
              }
              else {
                  bweight_los_deep=bweight_los_deep*(1-beff*jet_sf_los_deep); dweight_los_deep=dweight_los_deep*(1-beff);
                  bweightup_los_deep=bweightup_los_deep*(1-beff*jet_sfu_los_deep);
                  bweightdown_los_deep=bweightdown_los_deep*(1-beff*jet_sfd_los_deep);
              }
              
              
              //b jet tight wp DeepCSV
              if (ak4jet_pt[i]>=30.0 && ak4jet_pt[i]<50.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.520869700243;}
              else if (ak4jet_eta[i]<-0.8){beff=0.43937276358;}
              else if (ak4jet_eta[i]>-0.8){beff=0.452455595379;}
              }
              else if (ak4jet_pt[i]>=50.0 && ak4jet_pt[i]<70.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.563580454628;}
              else if (ak4jet_eta[i]<-0.8){beff=0.483288908117;}
              else if (ak4jet_eta[i]>-0.8){beff=0.495775694296;}
              }
              else if (ak4jet_pt[i]>=70.0 && ak4jet_pt[i]<100.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.57455874948;}
              else if (ak4jet_eta[i]<-0.8){beff=0.492246553058;}
              else if (ak4jet_eta[i]>-0.8){beff=0.503952520111;}
              }
              else if (ak4jet_pt[i]>=100.0 && ak4jet_pt[i]<140.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.571409970752;}
              else if (ak4jet_eta[i]<-0.8){beff=0.486007686416;}
              else if (ak4jet_eta[i]>-0.8){beff=0.495542760314;}
              }
              else if (ak4jet_pt[i]>=140.0 && ak4jet_pt[i]<200.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.554259448319;}
              else if (ak4jet_eta[i]<-0.8){beff=0.467166648355;}
              else if (ak4jet_eta[i]>-0.8){beff=0.474168752141;}
              }
              else if (ak4jet_pt[i]>=200.0 && ak4jet_pt[i]<300.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.501582653719;}
              else if (ak4jet_eta[i]<-0.8){beff=0.409965311191;}
              else if (ak4jet_eta[i]>-0.8){beff=0.415229321835;}
              }
              else if (ak4jet_pt[i]>=300.0 && ak4jet_pt[i]<600.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.411013411297;}
              else if (ak4jet_eta[i]<-0.8){beff=0.337036158045;}
              else if (ak4jet_eta[i]>-0.8){beff=0.339300846649;}
              }
              else if (ak4jet_pt[i]>=600.0 && ak4jet_pt[i]<1000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.292875166567;}
              else if (ak4jet_eta[i]<-0.8){beff=0.261624585856;}
              else if (ak4jet_eta[i]>-0.8){beff=0.259223994287;}
              }
              else if (ak4jet_pt[i]>=1000.0 && ak4jet_pt[i]<3000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.2232096184;}
              else if (ak4jet_eta[i]<-0.8){beff=0.241563055062;}
              else if (ak4jet_eta[i]>-0.8){beff=0.180357142857;}
              }
              if(ak4jet_icsv[i]>0.8953) {
                  bweight_tig_deep=bweight_tig_deep*beff*jet_sf_tig_deep; dweight_tig_deep=dweight_tig_deep*beff;
                  bweightup_tig_deep=bweightup_tig_deep*beff*jet_sfu_tig_deep;
                  bweightdown_tig_deep=bweightdown_tig_deep*beff*jet_sfd_tig_deep;
              }
              else {
                  bweight_tig_deep=bweight_tig_deep*(1-beff*jet_sf_tig_deep); dweight_tig_deep=dweight_tig_deep*(1-beff);
                  bweightup_tig_deep=bweightup_tig_deep*(1-beff*jet_sfu_tig_deep);
                  bweightdown_tig_deep=bweightdown_tig_deep*(1-beff*jet_sfd_tig_deep);
              }
              
              

          }
        else if(abs(ak4jet_hf[i])==4 && ak4jet_pt[i]<3000.)
          {   
            double jet_sf    = csv(1,ak4jet_pt[i]);
            double jet_sfu    = csv(2,ak4jet_pt[i]);
            double jet_sfd    = csv(3,ak4jet_pt[i]);
              double jet_sf_deep    = csv_deep(1,ak4jet_pt[i]);
              double jet_sfu_deep    = csv_deep(2,ak4jet_pt[i]);
              double jet_sfd_deep    = csv_deep(3,ak4jet_pt[i]);
              double jet_sf_tig_deep    = csv_tig_deep(1,ak4jet_pt[i]);
              double jet_sfu_tig_deep    = csv_tig_deep(2,ak4jet_pt[i]);
              double jet_sfd_tig_deep    = csv_tig_deep(3,ak4jet_pt[i]);
              double jet_sf_los_deep    = csv_los_deep(1,ak4jet_pt[i]);
              double jet_sfu_los_deep    = csv_los_deep(2,ak4jet_pt[i]);
              double jet_sfd_los_deep    = csv_los_deep(3,ak4jet_pt[i]);

//            int cbin=hc->FindBin(ak4jet_pt[i],ak4jet_eta[i]);
//            ceff=hc->GetBinContent(cbin);
              //c jet middle wp csv_2
              if (ak4jet_pt[i]>=30.0 && ak4jet_pt[i]<50.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.134641720495;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.121166477571;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.121321670268;}
              }
              else if (ak4jet_pt[i]>=50.0 && ak4jet_pt[i]<70.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.131091961461;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.117347838189;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.118499893002;}
              }
              else if (ak4jet_pt[i]>=70.0 && ak4jet_pt[i]<100.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.136123314875;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.124065366633;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.125565528143;}
              }
              else if (ak4jet_pt[i]>=100.0 && ak4jet_pt[i]<140.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.143276930452;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.130265186364;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.132288574989;}
              }
              else if (ak4jet_pt[i]>=140.0 && ak4jet_pt[i]<200.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.141473999343;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.137933649736;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.139164307147;}
              }
              else if (ak4jet_pt[i]>=200.0 && ak4jet_pt[i]<300.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.134456602805;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.129486943226;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.130881506146;}
              }
              else if (ak4jet_pt[i]>=300.0 && ak4jet_pt[i]<600.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.12118518329;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.139912735705;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.136867817858;}
              }
              else if (ak4jet_pt[i]>=600.0 && ak4jet_pt[i]<1000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.125290422511;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.152911249293;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.144447536877;}
              }
              else if (ak4jet_pt[i]>=1000.0 && ak4jet_pt[i]<3000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.120165745856;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.174418604651;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.156756756757;}
              }
              if(ak4jet_icsv[i]>0.8838) {
                  bweight=bweight*ceff*jet_sf; dweight=dweight*ceff;
                  bweightup=bweightup*ceff*jet_sfu;
                  bweightdown=bweightdown*ceff*jet_sfd;
              }
              else {
                  bweight=bweight*(1-ceff*jet_sf); dweight=dweight*(1-ceff);
                  bweightup=bweightup*(1-ceff*jet_sfu);
                  bweightdown=bweightdown*(1-ceff*jet_sfd);
              }
              //c jet middle wp DeepCSV
              if (ak4jet_pt[i]>=30.0 && ak4jet_pt[i]<50.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.129222981653;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.12950337414;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.131247230564;}
              }
              else if (ak4jet_pt[i]>=50.0 && ak4jet_pt[i]<70.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.117511658516;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.122105267909;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.123202524755;}
              }
              else if (ak4jet_pt[i]>=70.0 && ak4jet_pt[i]<100.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.11615914238;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.118675833028;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.119520303458;}
              }
              else if (ak4jet_pt[i]>=100.0 && ak4jet_pt[i]<140.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.121660554992;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.121737542076;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.121667035755;}
              }
              else if (ak4jet_pt[i]>=140.0 && ak4jet_pt[i]<200.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.129610177752;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.129728045772;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.12830101997;}
              }
              else if (ak4jet_pt[i]>=200.0 && ak4jet_pt[i]<300.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.135430225509;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.130401636845;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.128649540396;}
              }
              else if (ak4jet_pt[i]>=300.0 && ak4jet_pt[i]<600.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.144650004058;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.141602206632;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.13753036073;}
              }
              else if (ak4jet_pt[i]>=600.0 && ak4jet_pt[i]<1000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.152528207271;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.15451577802;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.160846851801;}
              }
              else if (ak4jet_pt[i]>=1000.0 && ak4jet_pt[i]<3000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.138200782269;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.154255319149;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.178378378378;}
              }
              if(ak4jet_icsv[i]>0.6321) {
                  bweight_deep=bweight_deep*ceff*jet_sf_deep; dweight_deep=dweight_deep*ceff;
                  bweightup_deep=bweightup_deep*ceff*jet_sfu_deep;
                  bweightdown_deep=bweightdown_deep*ceff*jet_sfd_deep;
              }
              else {
                  bweight_deep=bweight_deep*(1-ceff*jet_sf_deep); dweight_deep=dweight_deep*(1-ceff);
                  bweightup_deep=bweightup_deep*(1-ceff*jet_sfu_deep);
                  bweightdown_deep=bweightdown_deep*(1-ceff*jet_sfd_deep);
              }
              
              
              //c jet loose wp DeepCSV
              if (ak4jet_pt[i]>=30.0 && ak4jet_pt[i]<50.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.407532419829;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.414756479434;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.41552030184;}
              }
              else if (ak4jet_pt[i]>=50.0 && ak4jet_pt[i]<70.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.396198586996;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.406929414522;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.40737402914;}
              }
              else if (ak4jet_pt[i]>=70.0 && ak4jet_pt[i]<100.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.395266574669;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.405871869675;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.406265413537;}
              }
              else if (ak4jet_pt[i]>=100.0 && ak4jet_pt[i]<140.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.404117752108;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.413086131753;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.411639878628;}
              }
              else if (ak4jet_pt[i]>=140.0 && ak4jet_pt[i]<200.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.417210082833;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.429569432325;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.425682017993;}
              }
              else if (ak4jet_pt[i]>=200.0 && ak4jet_pt[i]<300.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.42772089264;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.442946517278;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.436942165497;}
              }
              else if (ak4jet_pt[i]>=300.0 && ak4jet_pt[i]<600.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.464082567948;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.496432212029;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.486341419533;}
              }
              else if (ak4jet_pt[i]>=600.0 && ak4jet_pt[i]<1000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.504972837443;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.55359085963;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.554577948859;}
              }
              else if (ak4jet_pt[i]>=1000.0 && ak4jet_pt[i]<3000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.552803129074;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.622340425532;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.578378378378;}
              }
              if(ak4jet_icsv[i]>0.2217) {
                  bweight_los_deep=bweight_los_deep*ceff*jet_sf_los_deep; dweight_los_deep=dweight_los_deep*ceff;
                  bweightup_los_deep=bweightup_los_deep*ceff*jet_sfu_los_deep;
                  bweightdown_los_deep=bweightdown_los_deep*ceff*jet_sfd_los_deep;
              }
              else {
                  bweight_los_deep=bweight_los_deep*(1-ceff*jet_sf_los_deep); dweight_los_deep=dweight_los_deep*(1-ceff);
                  bweightup_los_deep=bweightup_los_deep*(1-ceff*jet_sfu_los_deep);
                  bweightdown_los_deep=bweightdown_los_deep*(1-ceff*jet_sfd_los_deep);
              }
              
              
              //c jet tight wp DeepCSV
              if (ak4jet_pt[i]>=30.0 && ak4jet_pt[i]<50.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.0235693322031;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.0245835119498;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.0251726314782;}
              }
              else if (ak4jet_pt[i]>=50.0 && ak4jet_pt[i]<70.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.0221002199997;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.0236972600643;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.0241362178695;}
              }
              else if (ak4jet_pt[i]>=70.0 && ak4jet_pt[i]<100.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.0227875353133;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.0234240379638;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.0238317547109;}
              }
              else if (ak4jet_pt[i]>=100.0 && ak4jet_pt[i]<140.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.0248622843829;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.0242885304979;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.0248959753838;}
              }
              else if (ak4jet_pt[i]>=140.0 && ak4jet_pt[i]<200.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.0273376534458;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.0269927593883;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.0266163842966;}
              }
              else if (ak4jet_pt[i]>=200.0 && ak4jet_pt[i]<300.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.0290333494669;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.0261314787673;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.0257881920605;}
              }
              else if (ak4jet_pt[i]>=300.0 && ak4jet_pt[i]<600.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.0297079776972;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.028077591893;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.0269123511949;}
              }
              else if (ak4jet_pt[i]>=600.0 && ak4jet_pt[i]<1000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.0308399498537;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.0329162132753;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.0299697552928;}
              }
              else if (ak4jet_pt[i]>=1000.0 && ak4jet_pt[i]<3000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.0286831812256;}
              else if (ak4jet_eta[i]<-0.8){ceff=0.0212765957447;}
              else if (ak4jet_eta[i]>-0.8){ceff=0.0540540540541;}
              }
              if(ak4jet_icsv[i]>0.8953) {
                  bweight_tig_deep=bweight_tig_deep*ceff*jet_sf_tig_deep; dweight_tig_deep=dweight_tig_deep*ceff;
                  bweightup_tig_deep=bweightup_tig_deep*ceff*jet_sfu_tig_deep;
                  bweightdown_tig_deep=bweightdown_tig_deep*ceff*jet_sfd_tig_deep;
              }
              else {
                  bweight_tig_deep=bweight_tig_deep*(1-ceff*jet_sf_tig_deep); dweight_tig_deep=dweight_tig_deep*(1-ceff);
                  bweightup_tig_deep=bweightup_tig_deep*(1-ceff*jet_sfu_tig_deep);
                  bweightdown_tig_deep=bweightdown_tig_deep*(1-ceff*jet_sfd_tig_deep);
              }

          }
        else if (abs(ak4jet_hf[i])==0 && ak4jet_pt[i]<3000.)
          {  
            double jet_sf    = lsv(1,ak4jet_pt[i]);
            double jet_sfu    = lsv(2,ak4jet_pt[i]);
            double jet_sfd    = lsv(3,ak4jet_pt[i]);
              double jet_sf_deep    = lsv_deep(1,ak4jet_pt[i]);
              double jet_sfu_deep    = lsv_deep(2,ak4jet_pt[i]);
              double jet_sfd_deep    = lsv_deep(3,ak4jet_pt[i]);
              double jet_sf_tig_deep    = lsv_tig_deep(1,ak4jet_pt[i]);
              double jet_sfu_tig_deep    = lsv_tig_deep(2,ak4jet_pt[i]);
              double jet_sfd_tig_deep    = lsv_tig_deep(3,ak4jet_pt[i]);
              double jet_sf_los_deep    = lsv_los_deep(1,ak4jet_pt[i]);
              double jet_sfu_los_deep    = lsv_los_deep(2,ak4jet_pt[i]);
              double jet_sfd_los_deep    = lsv_los_deep(3,ak4jet_pt[i]);

//            int lbin=hl->FindBin(ak4jet_pt[i],ak4jet_eta[i]);
//            leff=hl->GetBinContent(lbin);
              //l jet middle wp csv_2
              if (ak4jet_pt[i]>=20.0 && ak4jet_pt[i]<200.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.00883804539947;}
              else if (ak4jet_eta[i]<-0.8){leff=0.0126624606715;}
              else if (ak4jet_eta[i]>-0.8){leff=0.0118134413393;}
              }
              else if (ak4jet_pt[i]>=200.0 && ak4jet_pt[i]<500.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.0108378571471;}
              else if (ak4jet_eta[i]<-0.8){leff=0.0209253471333;}
              else if (ak4jet_eta[i]>-0.8){leff=0.0202266537294;}
              }
              else if (ak4jet_pt[i]>=500.0 && ak4jet_pt[i]<1000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.0194805814586;}
              else if (ak4jet_eta[i]<-0.8){leff=0.0362205015431;}
              else if (ak4jet_eta[i]>-0.8){leff=0.0329325280287;}
              }
              else if (ak4jet_pt[i]>=1000.0 && ak4jet_pt[i]<3000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.0339058666072;}
              else if (ak4jet_eta[i]<-0.8){leff=0.0450112528132;}
              else if (ak4jet_eta[i]>-0.8){leff=0.0479452054795;}
              }
              if(ak4jet_icsv[i]>0.8838) {
                  bweight=bweight*leff*jet_sf; dweight=dweight*leff;
                  bweightup=bweightup*leff*jet_sfu;
                  bweightdown=bweightdown*leff*jet_sfd;
              }
              else {
                  bweight=bweight*(1-leff*jet_sf); dweight=dweight*(1-leff);
                  bweightup=bweightup*(1-leff*jet_sfu);
                  bweightdown=bweightdown*(1-leff*jet_sfd);
              }
              //l jet middle wp DeepCSV
              if (ak4jet_pt[i]>=20.0 && ak4jet_pt[i]<200.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.011375461494;}
              else if (ak4jet_eta[i]<-0.8){leff=0.0166255497284;}
              else if (ak4jet_eta[i]>-0.8){leff=0.0163500521057;}
              }
              else if (ak4jet_pt[i]>=200.0 && ak4jet_pt[i]<500.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.0136983493182;}
              else if (ak4jet_eta[i]<-0.8){leff=0.0207888628102;}
              else if (ak4jet_eta[i]>-0.8){leff=0.0196059086432;}
              }
              else if (ak4jet_pt[i]>=500.0 && ak4jet_pt[i]<1000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.0297503594554;}
              else if (ak4jet_eta[i]<-0.8){leff=0.0386133774486;}
              else if (ak4jet_eta[i]>-0.8){leff=0.0353440504808;}
              }
              else if (ak4jet_pt[i]>=1000.0 && ak4jet_pt[i]<3000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.0398640996602;}
              else if (ak4jet_eta[i]<-0.8){leff=0.0518292682927;}
              else if (ak4jet_eta[i]>-0.8){leff=0.0468509984639;}
              }
              if(ak4jet_icsv[i]>0.6321) {
                  bweight_deep=bweight_deep*leff*jet_sf_deep; dweight_deep=dweight_deep*leff;
                  bweightup_deep=bweightup_deep*leff*jet_sfu_deep;
                  bweightdown_deep=bweightdown_deep*leff*jet_sfd_deep;
              }
              else {
                  bweight_deep=bweight_deep*(1-leff*jet_sf_deep); dweight_deep=dweight_deep*(1-leff);
                  bweightup_deep=bweightup_deep*(1-leff*jet_sfu_deep);
                  bweightdown_deep=bweightdown_deep*(1-leff*jet_sfd_deep);
              }
              
              
              //l jet loose wp DeepCSV
              if (ak4jet_pt[i]>=20.0 && ak4jet_pt[i]<200.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.0955787532571;}
              else if (ak4jet_eta[i]<-0.8){leff=0.148789610047;}
              else if (ak4jet_eta[i]>-0.8){leff=0.145664194575;}
              }
              else if (ak4jet_pt[i]>=200.0 && ak4jet_pt[i]<500.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.113664749688;}
              else if (ak4jet_eta[i]<-0.8){leff=0.180917356463;}
              else if (ak4jet_eta[i]>-0.8){leff=0.173081118715;}
              }
              else if (ak4jet_pt[i]>=500.0 && ak4jet_pt[i]<1000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.216589950689;}
              else if (ak4jet_eta[i]<-0.8){leff=0.2958792345;}
              else if (ak4jet_eta[i]>-0.8){leff=0.286283052885;}
              }
              else if (ak4jet_pt[i]>=1000.0 && ak4jet_pt[i]<3000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.334314835787;}
              else if (ak4jet_eta[i]<-0.8){leff=0.383384146341;}
              else if (ak4jet_eta[i]>-0.8){leff=0.357142857143;}
              }
              if(ak4jet_icsv[i]>0.2217) {
                  bweight_los_deep=bweight_los_deep*leff*jet_sf_los_deep; dweight_los_deep=dweight_los_deep*leff;
                  bweightup_los_deep=bweightup_los_deep*leff*jet_sfu_los_deep;
                  bweightdown_los_deep=bweightdown_los_deep*leff*jet_sfd_los_deep;
              }
              else {
                  bweight_los_deep=bweight_los_deep*(1-leff*jet_sf_los_deep); dweight_los_deep=dweight_los_deep*(1-leff);
                  bweightup_los_deep=bweightup_los_deep*(1-leff*jet_sfu_los_deep);
                  bweightdown_los_deep=bweightdown_los_deep*(1-leff*jet_sfd_los_deep);
              }
              
              
              //l jet tight wp DeepCSV
              if (ak4jet_pt[i]>=20.0 && ak4jet_pt[i]<200.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.00115701280759;}
              else if (ak4jet_eta[i]<-0.8){leff=0.00153811907691;}
              else if (ak4jet_eta[i]>-0.8){leff=0.0015196417698;}
              }
              else if (ak4jet_pt[i]>=200.0 && ak4jet_pt[i]<500.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.00174898076121;}
              else if (ak4jet_eta[i]<-0.8){leff=0.00228349093204;}
              else if (ak4jet_eta[i]>-0.8){leff=0.00217531498754;}
              }
              else if (ak4jet_pt[i]>=500.0 && ak4jet_pt[i]<1000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.00353182989318;}
              else if (ak4jet_eta[i]<-0.8){leff=0.00475617550852;}
              else if (ak4jet_eta[i]>-0.8){leff=0.00368088942308;}
              }
              else if (ak4jet_pt[i]>=1000.0 && ak4jet_pt[i]<3000.0)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.00724801812005;}
              else if (ak4jet_eta[i]<-0.8){leff=0.00533536585366;}
              else if (ak4jet_eta[i]>-0.8){leff=0.00614439324117;}
              }
              if(ak4jet_icsv[i]>0.8953) {
                  bweight_tig_deep=bweight_tig_deep*leff*jet_sf_tig_deep; dweight_tig_deep=dweight_tig_deep*leff;
                  bweightup_tig_deep=bweightup_tig_deep*leff*jet_sfu_tig_deep;
                  bweightdown_tig_deep=bweightdown_tig_deep*leff*jet_sfd_tig_deep;
              }
              else {
                  bweight_tig_deep=bweight_tig_deep*(1-leff*jet_sf_tig_deep); dweight_tig_deep=dweight_tig_deep*(1-leff);
                  bweightup_tig_deep=bweightup_tig_deep*(1-leff*jet_sfu_tig_deep);
                  bweightdown_tig_deep=bweightdown_tig_deep*(1-leff*jet_sfd_tig_deep);
              }

          }
      //std::cout<<" yy1 "<<abs(ak4jet_hf[i])<<" "<<dweight<<" "<<bweight<<std::endl;
        } 
      //std::cout<<" yy2 "<<abs(ak4jet_hf[i])<<" "<<ak4jet_pt[i]<<" "<<fabs(ak4jet_eta[i])<<" "<<ak4jet_IDLoose[i]<<" "<<deltaRAK4AK8[i]<<std::endl;
     }
//     input->Close();
//     std::cout<<" xx "<<bweight/dweight<<" "<<bweightup/dweight<<" "<<bweightdown/dweight<<std::endl;
        

     btagweight_center=bweight/dweight;
     btagweight_up=bweightup/dweight;
     btagweight_down=bweightdown/dweight;
     btagweight_center_deep=bweight_deep/dweight_deep;
     btagweight_up_deep=bweightup_deep/dweight_deep;
     btagweight_down_deep=bweightdown_deep/dweight_deep;

            btagweight_center_tig_deep=bweight_tig_deep/dweight_tig_deep;
            btagweight_up_tig_deep=bweightup_tig_deep/dweight_tig_deep;
            btagweight_down_tig_deep=bweightdown_tig_deep/dweight_tig_deep;
            btagweight_center_los_deep=bweight_los_deep/dweight_los_deep;
            btagweight_up_los_deep=bweightup_los_deep/dweight_los_deep;
            btagweight_down_los_deep=bweightdown_los_deep/dweight_los_deep;

//     std::cout<<btagweight_center<<std::endl;
 
     
//--BSF----------------------------
        }
                if(theWeight>0) nn=1;
                else nn= -1;
                if(npp>0) lumiWeight=Identical_lumiWeight/(npp-nmm);
                else lumiWeight=Identical_lumiWeight/nentries;
//                std::cout<<IDweight<<"        "<<btagweight_center<<"   "<<trigger_eff<<std::endl;
                weight=lumiWeight*triggerWeight*eff_and_pu_Weight*nn*trigger_eff*IDweight*IDweightISO*IDweighttrk*ToppTweight;
        //if(weight==0) cout << "  lumiWeight  "<<lumiWeight<<"  triggerWeight   "<<triggerWeight<< " eff_and_pu_Weight:" << eff_and_pu_Weight << "  trigger_eff  "<<trigger_eff<<"  IDweight   "<<IDweight<<"  IDweightISO   "<<IDweightISO<<"  IDweighttrk  "<<IDweighttrk<<"  ToppTweight   "<<ToppTweight<< endl;
        //cout <<"gentop_pt  "<<gentop_pt<< "  etalep1  "<<etalep1<<"  ptlep1   "<<ptlep1<< " eff_and_pu_Weight:" << eff_and_pu_Weight << "  trigger_eff  "<<trigger_eff<<"  IDweight   "<<IDweight<<"  IDweightISO   "<<IDweightISO<<"  IDweighttrk  "<<IDweighttrk<<"  ToppTweight   "<<ToppTweight<< endl;

        //https://github.com/jmhogan/GenHTweight/blob/master/WJetsToLNuSFs.txt
                if (IsData==2 ) weight = weight*1.21*0.998056;
                if (IsData==3 ) weight = weight*1.21*0.978569;
                if (IsData==4 ) weight = weight*1.21*0.928054;
                if (IsData==5 ) weight = weight*1.21*0.856705;
                if (IsData==6 ) weight = weight*1.21*0.757463;
                if (IsData==7 ) weight = weight*1.21*0.608292;
                if (IsData==8 ) weight = weight*1.21*0.454246;
                if ( IsData==0 ) weight=1;
//final weight------
        scale_center=1;scale_down=1;scale_up=1;
        Double_t pdf_center=1,pdf_down=1,pdf_up=1;
        
        /*if(IsData!=11&&IsData>0) {
            scale_center=pweight[0];
            double div[8];
            for(Int_t i=0;i<8;i++) div[i]=(pweight[i+1]-pweight[0])/pweight[0];
            Int_t *indexpw=new Int_t[8];
            TMath::Sort(8,div,indexpw,1);
            scale_down=scale_center+div[indexpw[7]];scale_up=scale_center+div[indexpw[0]];
            pdf_center=pweight[110];
             double mpdf=0,pdf=0;
            for(int j=110;j<211;j++)    mpdf+=pweight[j];
            for(int j=110;j<211;j++)
                pdf+=pow((pweight[j]-mpdf/101),2)/100;
             pdf_down=pdf_center-mpdf;pdf_up=pdf_center+mpdf;
         
        }*/
        Double_t tau_corr1=1,tau_corr2=1;
        /*
        if((IsData==9||IsData==11)&&(access("TauCorrection_top_mu.root",R_OK)==0)){
        TFile * input_tau_top_corr = new TFile ("TauCorrection_top_mu.root");
        input_tau_top_corr->cd("");
        TH2D* tau_corr= (TH2D*) input_tau_top_corr->Get("top_tau_corr_use_ST_0");
        Double_t jet_mass1=jet_mass_puppi_corr,jet_mass2=jet_mass_puppi_corr_2;
        if (jet_mass1>=50&&jet_mass1<100) {tau_corr1=tau_corr->GetBinContent(tau_corr->GetXaxis()->FindBin(jet_mass1),tau_corr->GetYaxis()->FindBin(jetAK8puppi_tau21));}
        if (jet_mass1>=100) {if(jet_mass1>=200) jet_mass1=210.0;
            tau_corr1=tau_corr->GetBinContent(tau_corr->GetXaxis()->FindBin(jet_mass1),tau_corr->GetYaxis()->FindBin(jetAK8puppi_tau42));}
        if (jet_mass2>=50&&jet_mass2<100) {tau_corr2=tau_corr->GetBinContent(tau_corr->GetXaxis()->FindBin(jet_mass2),tau_corr->GetYaxis()->FindBin(jetAK8puppi_tau21_2));}
        if (jet_mass2>=100) {if(jet_mass2>=200) jet_mass2=210.0;
            tau_corr2=tau_corr->GetBinContent(tau_corr->GetXaxis()->FindBin(jet_mass2),tau_corr->GetYaxis()->FindBin(jetAK8puppi_tau42_2));}
            //cout<<jet_mass1<<"  "<<tau_corr->GetXaxis()->FindBin(jet_mass1)<<"   "<<tau_corr->GetYaxis()->FindBin(jetAK8puppi_tau42))<<endl;
        input_tau_top_corr->Close();
        }*/
        /*
         Double_t corr_tau_top[8][8]={0.88368171179,1.08142263609,0.952048715331,1.01450142262,1.01740012028,0.992085999379,0.890395001515,1.85377913532,
            1.17163189105,1.12633473475,1.00520190811,0.910652497628,1.05384696569,0.896949536522,0.981140063552,1.00667824318,
            0.986215256915,1.0762916331,0.868096775328,1.24566676532,0.897194362079,0.822301332844,1.22386486371,2.63218300907,
            0.916857375252,0.984275971648,1.13991616446,0.957040265942,0.904987000019,1.05658482069,0.843908811455,1.44428185639,
            1.02959773418,1.05697462236,1.00287218061,1.00780327896,0.993209893351,0.776684418614,0.959931632303,0.828227969764,
            0.379736559518,1.31582250636,1.14827881494,1.03337524737,1.07237138028,0.995844680494,0.881326166252,0.738593723698,
            1.0789446713,1.00926262552,1.00093511771,1.05236614269,0.995073140143,1.05023448361,0.82395499251,0.568156815962,
            3.43926566453,1.10303959852,1.00197029873,0.996030514591,0.990453941042,1.0512168889,0.865403960178,0.563547976891};*/
        Double_t corr_tau_top[8][8]={0};
    if( channelname=="mu" ){
        Double_t corr_tau_top_tmp[8][8]={1.03402672187,1.13558084027,1.00360681557,1.00246285777,0.99502594773,0.941244969132,0.873699375981,1.31334378138,
            1.27402792202,1.20715713692,1.02959119179,0.899276775573,1.02657444441,0.884478365492,0.841185385275,0.728595008922,
            1.05702477048,1.12832937166,0.895318658766,1.23145431738,0.842625831911,0.750238436967,1.09420467154,1.91003197061,
            0.980603593355,1.02947129926,1.129511571,0.949527312819,0.868362761076,0.951283339336,0.743157927249,1.12065620372,
            1.08967172869,1.11645005347,0.993568975319,0.993556032941,0.955014915484,0.727630260733,0.855622724783,0.737436065466,
            0.402541005911,1.38553543957,1.17279430632,1.0355029671,1.04503566577,0.991063708437,0.888307528307,0.762422249715,
            1.14159534938,1.05231203514,1.02810296348,1.05506976435,0.952612431097,0.983597829862,0.811889597808,0.543743762695,
            3.38517432139,1.16865887847,1.04707550932,0.967030785779,0.970556964857,1.0061406186,0.870983443008,0.532607396537};
        for(int i=0;i<8;i++){for(int j=0;j<8;j++) corr_tau_top[i][j]=corr_tau_top_tmp[i][j];}
    }
        if( channelname=="el" ){
            Double_t corr_tau_top_tmp[8][8]={1.41025975583,1.15858307886,1.0452690859,1.07341169412,0.971667001311,0.761359522858,0.968412543795,1.11905262524,
                1.06615963681,0.991362674834,1.0832466723,1.00584033907,0.921137288222,0.911435773944,1.16787220624,1.16606543631,
                1.01413792196,1.08350052267,1.05526182615,0.884223867405,0.999574049646,0.956885484244,0.993222583722,0.624462847374,
                1.07504405837,1.10128884231,0.96863454652,1.03623323744,0.985294528506,0.740692055289,0.724261698936,0.285626769713,
                0.991639809344,1.08779087537,0.999179928681,0.989316527303,0.979710236647,0.868862648216,0.812155977336,0.0,
                0.904798527337,1.07531968905,1.06744039711,1.02579027671,0.961549277631,1.07560869706,0.936769114354,0.791346031794,
                0.670408713531,1.06802236901,1.08949541824,0.964228451804,0.905098361023,1.05957081004,0.962326062293,0.414373981633,
                0.0,1.36172365544,0.862197044607,1.09692614851,1.06289307319,0.68202643517,1.10112065357,1.73629657039};
            for(int i=0;i<8;i++){for(int j=0;j<8;j++) corr_tau_top[i][j]=corr_tau_top_tmp[i][j];}
        }

        if(IsData==1||IsData==9||IsData==11||IsData==10){
        Double_t jet_mass1=jet_mass_puppi_corr,jet_mass2=jet_mass_puppi_corr_2;
            Double_t taubin[9]={0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1};
            Double_t massbin[9]={50,60,70,80,90,100,150,200,1000000};
            int ix=0;
            for(ix=0;ix<8;ix++){
                if(jet_mass1<massbin[0]) {ix=-1;break;}
                if(jet_mass1>=massbin[ix]&&jet_mass1<massbin[ix+1]) break;
                }
            int iy=0;
            for(iy=0;iy<8;iy++){
                if(jet_mass1<100&&(jetAK8puppi_tau21<0||jetAK8puppi_tau21>1)) {iy=-1;break;}
                if(jet_mass1>=100&&(jetAK8puppi_tau42<0||jetAK8puppi_tau42>1)) {iy=-1;break;}
                if(jet_mass1<100&&(jetAK8puppi_tau21>=taubin[iy]&&jetAK8puppi_tau21<taubin[iy+1])) break;
                if(jet_mass1>=100&&(jetAK8puppi_tau42>=taubin[iy]&&jetAK8puppi_tau42<taubin[iy+1])) break;
            }
            //cout<<ix<<iy<<endl;
            if(ix>=0&&iy>=0) tau_corr1=corr_tau_top[ix][iy];
            int ix2=0;
            for(ix2=0;ix2<8;ix2++){
                if(jet_mass2<massbin[0]) {ix2=-1;break;}
                if(jet_mass2>=massbin[ix2]&&jet_mass2<massbin[ix2+1]) break;
            }
            int iy2=0;
            for(iy2=0;iy2<8;iy2++){
                if(jet_mass2<100&&(jetAK8puppi_tau21_2<0||jetAK8puppi_tau21_2>1)) {iy2=-1;break;}
                if(jet_mass2>=100&&(jetAK8puppi_tau42_2<0||jetAK8puppi_tau42_2>1)) {iy2=-1;break;}
                if(jet_mass2<100&&(jetAK8puppi_tau21_2>=taubin[iy2]&&jetAK8puppi_tau21_2<taubin[iy2+1])) break;
                if(jet_mass2>=100&&(jetAK8puppi_tau42_2>=taubin[iy2]&&jetAK8puppi_tau42_2<taubin[iy2+1])) break;
            }
            if(ix2>=0&&iy2>=0) tau_corr2=corr_tau_top[ix2][iy2];
            //if(ix==1||ix2==1) cout<<jet_mass1<<"  "<<jet_mass2<<"  "<<ix<<"   "<<iy<<"   "<<ix2<<"   "<<iy2<<"  "<<tau_corr1<<"  "<<tau_corr2<<endl;
        }
        if(tau_corr1==0) tau_corr1=1;
        if(tau_corr2==0) tau_corr2=1;
        //cout<<tau_corr1<<endl;
        double weight_tmp=weight;
        weight=btagweight_center*weight;
        weight_tautopcorr1=weight*tau_corr1*tau_corr2;
        weight_tautopcorr2=weight*tau_corr1*tau_corr2*tau_corr1*tau_corr2;
        
        Double_t corr_tau_wjets[8][8]={0};
        if( channelname=="mu" ){
            /*Double_t corr_tau_wjets_tmp[8][8]={0.890338365029,1.36951339498,1.21829324738,1.04545353353,0.995498767509,0.805068336795,0.722901182608,0.832207091212,
             1.16017385313,1.18981503432,1.02978323511,1.13913303627,0.921533456764,0.86787873397,0.725135022092,0.70267716147,
             1.11936289195,1.12168078333,1.14756670293,0.978534224507,0.955173074638,0.871852454955,0.625612240463,0.704241418756,
             1.07420405138,1.13400570899,1.06806934222,1.03140715048,0.965049849752,0.715746556747,0.747355563966,0.317885258035,
             1.09893144453,1.15117667382,1.07886079746,1.05295843998,0.874291679157,0.814226043194,0.478201367609,0.513121488775,
             4.06075521949,1.02619590385,1.0907372046,0.99956946094,0.966399210913,1.03230481645,0.980526429312,0.919910712067,
             1.15555399072,0.987914427408,0.983878342707,1.00971148859,0.906918387498,1.07184436283,1.0881906704,0.838550936066,
             0.758809700766,1.21311137728,0.902559232848,0.848916725301,0.912020239288,1.15903946763,1.22424273191,1.08247703796};*/
            Double_t corr_tau_wjets_tmp[8][8]={0.878093740102,1.35156804567,1.21240539342,1.04662781583,0.994611210851,0.810461996416,0.7266833886,0.815895651623,
                1.09906271047,1.14509245011,1.02713684287,1.14865862365,0.926744533372,0.875857159651,0.726745853445,0.706692247635,
                1.08368794285,1.07212760852,1.16308231664,0.967811047954,0.973180942193,0.89205332303,0.63844271766,0.681519372528,
                1.04284648379,1.11675705123,1.05408987696,1.04076524132,0.983030478416,0.723921597151,0.767649447996,0.326313172076,
                1.07258152449,1.11993999688,1.07836471575,1.0553622067,0.885342143572,0.832897953632,0.492360188901,0.525759746067,
                4.26894599699,0.957789410667,1.07579486091,0.989669333318,0.968463342718,1.03636625932,0.986605654291,0.924427560405,
                1.06961238475,0.961331819599,0.963633250292,0.999552852707,0.925391876451,1.08318790875,1.11442288616,0.869397442857,
                0.738869452843,1.15491886733,0.89828967094,0.846040089658,0.913241660166,1.17428972598,1.24037277721,1.10446103014};
            for(int i=0;i<8;i++){for(int j=0;j<8;j++) corr_tau_wjets[i][j]=corr_tau_wjets_tmp[i][j];}
        }
        if( channelname=="el" ){
            Double_t corr_tau_wjets_tmp[8][8]={1.17292641723,1.10891790676,1.12095389311,0.973386080896,1.05078999779,0.905922788641,0.836077714046,0.858081094185,
                1.28701282261,1.14259715244,1.14058794551,1.05418747543,0.996798379208,0.817133813464,0.666248461507,0.787187559044,
                1.13685253873,1.10369052406,1.10811962364,1.1022258448,0.938749995542,0.750294624097,0.699666712499,0.71074161622,
                0.97332603372,1.03820385269,1.16188763707,1.02336152663,0.879555321097,0.895620009769,0.67726050421,0.178077617303,
                0.559131233724,0.651974135154,0.924824246403,1.20979144196,1.61383950661,1.71306880469,2.04605542934,0.679655808092,
                0.859778527345,1.08807037348,0.997113874532,0.991518733497,1.05233073551,1.06841187914,0.911995001616,0.754012650151,
                1.96017343947,0.935163911645,0.920617472488,0.931746105761,1.06198942512,1.01536409014,1.11621738285,1.326943875,
                1.90352970457,0.919726446215,0.934034407982,0.951325045954,0.89204333586,1.26855244267,0.984714168594,1.11868970913};
            for(int i=0;i<8;i++){for(int j=0;j<8;j++) corr_tau_wjets[i][j]=corr_tau_wjets_tmp[i][j];}
        }
        
        if((IsData>=2&&IsData<=8)){
            Double_t jet_mass1=jet_mass_puppi_corr,jet_mass2=jet_mass_puppi_corr_2;
            Double_t taubin[9]={0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1};
            Double_t massbin[9]={50,60,70,80,90,100,150,200,1000000};
            int ix=0;
            for(ix=0;ix<8;ix++){
                if(jet_mass1<massbin[0]) {ix=-1;break;}
                if(jet_mass1>=massbin[ix]&&jet_mass1<massbin[ix+1]) break;
            }
            int iy=0;
            for(iy=0;iy<8;iy++){
                if(jet_mass1<100&&(jetAK8puppi_tau21<0||jetAK8puppi_tau21>1)) {iy=-1;break;}
                if(jet_mass1>=100&&(jetAK8puppi_tau42<0||jetAK8puppi_tau42>1)) {iy=-1;break;}
                if(jet_mass1<100&&(jetAK8puppi_tau21>=taubin[iy]&&jetAK8puppi_tau21<taubin[iy+1])) break;
                if(jet_mass1>=100&&(jetAK8puppi_tau42>=taubin[iy]&&jetAK8puppi_tau42<taubin[iy+1])) break;
            }
            //cout<<ix<<iy<<endl;
            if(ix>=0&&iy>=0) tau_corr1=corr_tau_wjets[ix][iy];
            int ix2=0;
            for(ix2=0;ix2<8;ix2++){
                if(jet_mass2<massbin[0]) {ix2=-1;break;}
                if(jet_mass2>=massbin[ix2]&&jet_mass2<massbin[ix2+1]) break;
            }
            int iy2=0;
            for(iy2=0;iy2<8;iy2++){
                if(jet_mass2<100&&(jetAK8puppi_tau21_2<0||jetAK8puppi_tau21_2>1)) {iy2=-1;break;}
                if(jet_mass2>=100&&(jetAK8puppi_tau42_2<0||jetAK8puppi_tau42_2>1)) {iy2=-1;break;}
                if(jet_mass2<100&&(jetAK8puppi_tau21_2>=taubin[iy2]&&jetAK8puppi_tau21_2<taubin[iy2+1])) break;
                if(jet_mass2>=100&&(jetAK8puppi_tau42_2>=taubin[iy2]&&jetAK8puppi_tau42_2<taubin[iy2+1])) break;
            }
            if(ix2>=0&&iy2>=0) tau_corr2=corr_tau_wjets[ix2][iy2];
            
            
            if(tau_corr1==0) tau_corr1=1;
            if(tau_corr2==0) tau_corr2=1;
            //cout<<jet_mass1<<" a "<<jet_mass_puppi_corr<<"  "<<ix<<"  "<<jetAK8puppi_tau21<<"  "<<iy<<"  "<<tau_corr1<<endl;
            //cout<<jet_mass2<<" b "<<jet_mass_puppi_corr_2<<"  "<<ix2<<"  "<<jetAK8puppi_tau21_2<<"  "<<iy2<<"  "<<tau_corr2<<endl;
            //cout<<tau_corr1<<endl;
            weight_tautopcorr1=weight*tau_corr1*tau_corr2;
            weight_tautopcorr2=weight*tau_corr1*tau_corr2*tau_corr1*tau_corr2;
        }

        ///***** corr deepak8 ******///
        Double_t deepak8_corr1=1,deepak8_corr2=1;
        /*
         if((IsData==9||IsData==11)&&(access("DeepAK8Correction_top_lep.root",R_OK)==0)){
         TFile * input_deepak8_top_corr = new TFile ("deepak8Correction_top_mu.root");
         input_deepak8_top_corr->cd("");
         TH2D* deepak8_corr= (TH2D*) input_deepak8_top_corr->Get("top_deepak8_corr_use_ST_0");
         Double_t jet_mass1=jet_mass_puppi_corr,jet_mass2=jet_mass_puppi_corr_2;
         if (jet_mass1>=50&&jet_mass1<100) {deepak8_corr1=deepak8_corr->GetBinContent(deepak8_corr->GetXaxis()->FindBin(jet_mass1),deepak8_corr->GetYaxis()->FindBin(jetAK8puppi_dnnDecorrW));}
         if (jet_mass1>=100) {if(jet_mass1>=200) jet_mass1=210.0;
         deepak8_corr1=deepak8_corr->GetBinContent(deepak8_corr->GetXaxis()->FindBin(jet_mass1),deepak8_corr->GetYaxis()->FindBin(((jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q)/(jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q+jetAK8puppi_dnnDecorrqcd))));}
         if (jet_mass2>=50&&jet_mass2<100) {deepak8_corr2=deepak8_corr->GetBinContent(deepak8_corr->GetXaxis()->FindBin(jet_mass2),deepak8_corr->GetYaxis()->FindBin(jetAK8puppi_dnnDecorrW_2));}
         if (jet_mass2>=100) {if(jet_mass2>=200) jet_mass2=210.0;
         deepak8_corr2=deepak8_corr->GetBinContent(deepak8_corr->GetXaxis()->FindBin(jet_mass2),deepak8_corr->GetYaxis()->FindBin(((jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2)/(jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2+jetAK8puppi_dnnDecorrqcd_2))));}
         //cout<<jet_mass1<<"  "<<deepak8_corr->GetXaxis()->FindBin(jet_mass1)<<"   "<<deepak8_corr->GetYaxis()->FindBin(((jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q)/(jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q+jetAK8puppi_dnnDecorrqcd))))<<endl;
         input_deepak8_top_corr->Close();
         }*/
        Double_t corr_deepak8_top[4][12]={0};
        if( channelname=="mu"||channelname=="el" ){
            Double_t corr_deepak8_top_tmp[4][12]={0.644693981106,0.813884479014,0.9581094444,1.10369465581,1.06317264272,1.11860474771,1.01341532195,1.11330373755,1.16134793033,1.00208593929,1.15777162098,0.971844322721,
                0.709693803338,0.796987633814,1.09747903963,1.10616039205,0.981750115444,0.880359074993,1.12636126636,1.13281039152,1.21878327287,1.1731795458,1.17219072371,1.00616642582,
                0.694304423163,0.93189847275,1.06798882602,0.917641891085,1.01675453078,0.987941134099,1.00122560345,1.02471007355,1.0754614055,1.0754614055,1.05078356126,1.05078356126,
                0.947338757927,0.947338757927,1.04751871879,1.04751871879,0.906131230598,0.906131230598,1.04013842851,1.04013842851,1.0890765813,1.0890765813,1.0890765813,1.0890765813};
            for(int i=0;i<4;i++){for(int j=0;j<12;j++) corr_deepak8_top[i][j]=corr_deepak8_top_tmp[i][j];}
        }
        if(IsData==1||IsData==9||IsData==11||IsData==10){
            Double_t jet_mass1=jet_mass_puppi_corr,jet_mass2=jet_mass_puppi_corr_2;
            Double_t deepak8bin[13]={0,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
            Double_t massbin[5]={60,100,150,200,1000000};
            int ixx=0;
            for(ixx=0;ixx<4;ixx++){
                if(jet_mass1<massbin[0]) {ixx=-1;break;}
                if(jet_mass1>=massbin[ixx]&&jet_mass1<massbin[ixx+1]) break;
            }
            int iyy=0;
            for(iyy=0;iyy<12;iyy++){
                if(jet_mass1<100&&(jetAK8puppi_dnnDecorrW<0||jetAK8puppi_dnnDecorrW>1)) {iyy=-1;break;}
                if(jet_mass1>=100&&(((jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q)/(jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q+jetAK8puppi_dnnDecorrqcd))<0||((jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q)/(jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q+jetAK8puppi_dnnDecorrqcd))>1)) {iyy=-1;break;}
                if(jet_mass1<100&&(jetAK8puppi_dnnDecorrW>=deepak8bin[iyy]&&jetAK8puppi_dnnDecorrW<deepak8bin[iyy+1])) break;
                if(jet_mass1>=100&&(((jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q)/(jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q+jetAK8puppi_dnnDecorrqcd))>=deepak8bin[iyy]&&((jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q)/(jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q+jetAK8puppi_dnnDecorrqcd))<deepak8bin[iyy+1])) break;
            }
            //cout<<ixx<<iyy<<endl;
            if(ixx>=0&&iyy>=0) deepak8_corr1=corr_deepak8_top[ixx][iyy];
            int ixx2=0;
            for(ixx2=0;ixx2<4;ixx2++){
                if(jet_mass2<massbin[0]) {ixx2=-1;break;}
                if(jet_mass2>=massbin[ixx2]&&jet_mass2<massbin[ixx2+1]) break;
            }
            int iyy2=0;
            for(iyy2=0;iyy2<12;iyy2++){
                if(jet_mass2<100&&(jetAK8puppi_dnnDecorrW_2<0||jetAK8puppi_dnnDecorrW_2>1)) {iyy2=-1;break;}
                if(jet_mass2>=100&&(((jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2)/(jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2+jetAK8puppi_dnnDecorrqcd_2))<0||((jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2)/(jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2+jetAK8puppi_dnnDecorrqcd_2))>1)) {iyy2=-1;break;}
                if(jet_mass2<100&&(jetAK8puppi_dnnDecorrW_2>=deepak8bin[iyy2]&&jetAK8puppi_dnnDecorrW_2<deepak8bin[iyy2+1])) break;
                if(jet_mass2>=100&&(((jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2)/(jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2+jetAK8puppi_dnnDecorrqcd_2))>=deepak8bin[iyy2]&&((jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2)/(jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2+jetAK8puppi_dnnDecorrqcd_2))<deepak8bin[iyy2+1])) break;
            }
            if(ixx2>=0&&iyy2>=0) deepak8_corr2=corr_deepak8_top[ixx2][iyy2];
            //if(ixx==4&&(iyy==4)) cout<<jet_mass1<<"  "<<jetAK8puppi_dnnDecorrW<<"  "<<(((jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q)/(jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q+jetAK8puppi_dnnDecorrqcd)))<<"  "<<jet_mass2<<"  "<<ixx<<"   "<<iyy<<"   "<<ixx2<<"   "<<iyy2<<"  "<<deepak8_corr1<<"  "<<deepak8_corr2<<endl;
        
        if(deepak8_corr1==0) deepak8_corr1=1;
        if(deepak8_corr2==0) deepak8_corr2=1;

        weight_deep=btagweight_center_deep*weight_tmp;
        weight_deep_topcorr1=weight_deep*deepak8_corr1*deepak8_corr2;
        weight_deep_topcorr2=weight_deep*deepak8_corr1*deepak8_corr2*deepak8_corr1*deepak8_corr2;

            weight_tig_deep=btagweight_center_tig_deep*weight_tmp;
            weight_tig_deep_topcorr1=weight_tig_deep*deepak8_corr1*deepak8_corr2;
            weight_tig_deep_topcorr2=weight_tig_deep*deepak8_corr1*deepak8_corr2*deepak8_corr1*deepak8_corr2;
            weight_los_deep=btagweight_center_los_deep*weight_tmp;
            weight_los_deep_topcorr1=weight_los_deep*deepak8_corr1*deepak8_corr2;
            weight_los_deep_topcorr2=weight_los_deep*deepak8_corr1*deepak8_corr2*deepak8_corr1*deepak8_corr2;

        //if(jet_mass_puppi_corr>=100&&jet_mass_puppi_corr<150&&((jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q)/(jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q+jetAK8puppi_dnnDecorrqcd))>=0.8) cout<<deepak8_corr1<<"  "<<deepak8_corr2<<"  "<<(weight_deep_topcorr1/weight_deep)<<endl;
        }
        deepak8_corr1=1;deepak8_corr2=1;
        Double_t corr_deepak8_wjets[4][12]={0};
        if( channelname=="mu" ||channelname=="el" ){
            Double_t corr_deepak8_wjets_tmp[4][12]={0.654254918739,0.815210932655,0.910038032827,0.932839478126,1.05400461733,1.11882859013,1.23655146368,1.18968825811,1.19391549836,1.29924949799,1.19039910678,1.18338655006,
                0.700083721023,0.88425728413,0.874207706088,0.951375830523,1.07075137866,1.11385148772,1.13544437581,1.13784196743,1.09417735717,1.0740804446,1.24843160445,1.26607139424,
                0.912517340367,0.99157165686,1.0305264157,1.03556005263,0.984675213002,0.966808237809,0.989273515549,1.0195136825,1.01241770892,1.01241770892,1.05230576083,1.05230576083,
                1.03387166559,1.03387166559,0.865032759479,0.865032759479,0.984061741661,0.984061741661,1.00605199371,1.00605199371,1.11677164941,1.11677164941,1.11677164941,1.11677164941};
            for(int i=0;i<4;i++){for(int j=0;j<12;j++) corr_deepak8_wjets[i][j]=corr_deepak8_wjets_tmp[i][j];}
        }

        if((IsData>=2&&IsData<=8)){
            Double_t jet_mass1=jet_mass_puppi_corr,jet_mass2=jet_mass_puppi_corr_2;
            Double_t deepak8bin[13]={0,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
            Double_t massbin[5]={60,100,150,200,1000000};
            int ixx=0;
            for(ixx=0;ixx<4;ixx++){
                if(jet_mass1<massbin[0]) {ixx=-1;break;}
                if(jet_mass1>=massbin[ixx]&&jet_mass1<massbin[ixx+1]) break;
            }
            int iyy=0;
            for(iyy=0;iyy<12;iyy++){
                if(jet_mass1<100&&(jetAK8puppi_dnnDecorrW<0||jetAK8puppi_dnnDecorrW>1)) {iyy=-1;break;}
                if(jet_mass1>=100&&(((jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q)/(jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q+jetAK8puppi_dnnDecorrqcd))<0||((jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q)/(jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q+jetAK8puppi_dnnDecorrqcd))>1)) {iyy=-1;break;}
                if(jet_mass1<100&&(jetAK8puppi_dnnDecorrW>=deepak8bin[iyy]&&jetAK8puppi_dnnDecorrW<deepak8bin[iyy+1])) break;
                if(jet_mass1>=100&&(((jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q)/(jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q+jetAK8puppi_dnnDecorrqcd))>=deepak8bin[iyy]&&((jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q)/(jetAK8puppi_dnnDecorrw+jetAK8puppi_dnnDecorrh4q+jetAK8puppi_dnnDecorrqcd))<deepak8bin[iyy+1])) break;
            }
            //cout<<ixx<<iyy<<endl;
            if(ixx>=0&&iyy>=0) deepak8_corr1=corr_deepak8_wjets[ixx][iyy];
            int ixx2=0;
            for(ixx2=0;ixx2<4;ixx2++){
                if(jet_mass2<massbin[0]) {ixx2=-1;break;}
                if(jet_mass2>=massbin[ixx2]&&jet_mass2<massbin[ixx2+1]) break;
            }
            int iyy2=0;
            for(iyy2=0;iyy2<12;iyy2++){
                if(jet_mass2<100&&(jetAK8puppi_dnnDecorrW_2<0||jetAK8puppi_dnnDecorrW_2>1)) {iyy2=-1;break;}
                if(jet_mass2>=100&&(((jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2)/(jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2+jetAK8puppi_dnnDecorrqcd_2))<0||((jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2)/(jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2+jetAK8puppi_dnnDecorrqcd_2))>1)) {iyy2=-1;break;}
                if(jet_mass2<100&&(jetAK8puppi_dnnDecorrW_2>=deepak8bin[iyy2]&&jetAK8puppi_dnnDecorrW_2<deepak8bin[iyy2+1])) break;
                if(jet_mass2>=100&&(((jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2)/(jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2+jetAK8puppi_dnnDecorrqcd_2))>=deepak8bin[iyy2]&&((jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2)/(jetAK8puppi_dnnDecorrw_2+jetAK8puppi_dnnDecorrh4q_2+jetAK8puppi_dnnDecorrqcd_2))<deepak8bin[iyy2+1])) break;
            }
            if(ixx2>=0&&iyy2>=0) deepak8_corr2=corr_deepak8_wjets[ixx2][iyy2];
            
        
        if(deepak8_corr1==0) deepak8_corr1=1;
        if(deepak8_corr2==0) deepak8_corr2=1;
        weight_deep=btagweight_center_deep*weight_tmp;
        weight_deep_topcorr1=weight_deep*deepak8_corr1*deepak8_corr2;
        weight_deep_topcorr2=weight_deep*deepak8_corr1*deepak8_corr2*deepak8_corr1*deepak8_corr2;
            weight_tig_deep=btagweight_center_tig_deep*weight_tmp;
            weight_tig_deep_topcorr1=weight_tig_deep*deepak8_corr1*deepak8_corr2;
            weight_tig_deep_topcorr2=weight_tig_deep*deepak8_corr1*deepak8_corr2*deepak8_corr1*deepak8_corr2;
            weight_los_deep=btagweight_center_los_deep*weight_tmp;
            weight_los_deep_topcorr1=weight_los_deep*deepak8_corr1*deepak8_corr2;
            weight_los_deep_topcorr2=weight_los_deep*deepak8_corr1*deepak8_corr2*deepak8_corr1*deepak8_corr2;

        }
        if (IsData==0 ) {weight_deep=btagweight_center_deep*weight_tmp;weight_deep_topcorr1=weight_deep;weight_deep_topcorr2=weight_deep;
        weight_tig_deep=btagweight_center_tig_deep*weight_tmp;weight_tig_deep_topcorr1=weight_tig_deep;weight_tig_deep_topcorr2=weight_tig_deep;
            weight_los_deep=btagweight_center_los_deep*weight_tmp;weight_los_deep_topcorr1=weight_los_deep;weight_los_deep_topcorr2=weight_los_deep;}
        weight_deep_prefire=weight_deep*L1prefiring;
        weight_deep_prefireup=weight_deep*L1prefiringup;
        weight_deep_prefiredown=weight_deep*L1prefiringdown;
        weight_deep_lumi=weight_deep*L1prefiring*35.92;
        if (IsData==0 ){
            weight_deep_prefire=weight_deep;
            weight_deep_prefireup=weight_deep;
            weight_deep_prefiredown=weight_deep;
            weight_deep_lumi=weight_deep;
        }
        //number of bjet calculation
		num_bJet=0.;
		num_bJet_loose=0.;
		num_bJet_tight=0.;
        num_bJet_deep=0.;
        num_bJet_loose_deep=0.;
        num_bJet_tight_deep=0.;

		for(Int_t i=0; i<8; i++)  {
                        deltaRAK4AK8_new[i]=0.;
                        deltaRAK4AK8_new[i]=sqrt(pow(fabs(ak4jet_eta[i]-jetAK8puppi_eta),2)+pow(TMath::Min(fabs(ak4jet_phi[i]-jetAK8puppi_phi),2*Pi-fabs(ak4jet_phi[i]-jetAK8puppi_phi)),2));
            if(Nj8==2) deltaRAK4AK8_new_2[i]=sqrt(pow(fabs(ak4jet_eta[i]-jetAK8puppi_eta_2),2)+pow(TMath::Min(fabs(ak4jet_phi[i]-jetAK8puppi_phi_2),2*Pi-fabs(ak4jet_phi[i]-jetAK8puppi_phi_2)),2));

			if(ak4jet_pt[i]>30 && ak4jet_icsv[i]>0.8838 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8_new[i]>=0.8 && (Nj8==1 || deltaRAK4AK8_new_2[i]>=0.8)) {num_bJet=num_bJet+1;}
			if(ak4jet_pt[i]>30 && ak4jet_icsv[i]>0.5426 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8_new[i]>=0.8 && (Nj8==1 || deltaRAK4AK8_new_2[i]>=0.8)) {num_bJet_loose=num_bJet_loose+1;}
			if(ak4jet_pt[i]>30 && ak4jet_icsv[i]>0.9693 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8_new[i]>=0.8 && (Nj8==1 || deltaRAK4AK8_new_2[i]>=0.8)) {num_bJet_tight=num_bJet_tight+1;}
            if(ak4jet_pt[i]>30 && (ak4jet_deepcsvb[i]+ak4jet_deepcsvbb[i])>0.6321 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8_new[i]>=0.8 && (Nj8==1 || deltaRAK4AK8_new_2[i]>=0.8)) {num_bJet_deep=num_bJet_deep+1;}
            if(ak4jet_pt[i]>30 && (ak4jet_deepcsvb[i]+ak4jet_deepcsvbb[i])>0.2217 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8_new[i]>=0.8 && (Nj8==1 || deltaRAK4AK8_new_2[i]>=0.8)) {num_bJet_loose_deep=num_bJet_loose_deep+1;}
            if(ak4jet_pt[i]>30 && (ak4jet_deepcsvb[i]+ak4jet_deepcsvbb[i])>0.8953 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8_new[i]>=0.8 && (Nj8==1 || deltaRAK4AK8_new_2[i]>=0.8)) {num_bJet_tight_deep=num_bJet_tight_deep+1;}

		}
		nbtag=num_bJet;
        nbtag_deep=num_bJet_deep;

		//number of bjet calculation Done

		Int_t nLooseLep=nLooseEle+nLooseMu;//the tight Lep included
        
         if(Nj8==1)
        {jetAK8puppi_tau31_2=-99;jetAK8puppi_tau32_2=-99;jetAK8puppi_tau43_2=-99;jetAK8puppi_tau41_2=-99;jetAK8puppi_tau42_2=-99;jetAK8puppi_tau31_3=-99;jetAK8puppi_tau32_3=-99;jetAK8puppi_tau43_3=-99;jetAK8puppi_tau41_3=-99;jetAK8puppi_tau42_3=-99;}
        if(Nj8==2)
        {jetAK8puppi_tau32_3=-99;jetAK8puppi_tau43_3=-99;jetAK8puppi_tau41_3=-99;jetAK8puppi_tau42_3=-99;}

        deltaRjet1jet2=sqrt(pow(fabs(jetAK8puppi_eta_2-jetAK8puppi_eta),2)+pow(TMath::Min(fabs(jetAK8puppi_phi_2-jetAK8puppi_phi),2*Pi-fabs(jetAK8puppi_phi_2-jetAK8puppi_phi)),2));
        deltaetajet1jet2=jetAK8puppi_eta-jetAK8puppi_eta_2;
        deltaetajet1lep=jetAK8puppi_eta-etalep1;
        deltaetajet2lep=jetAK8puppi_eta_2-etalep1;
        delPhilepmet_m=TMath::Min(fabs(philep1-MET_phi),2*Pi-fabs(philep1-MET_phi));
        delPhiWjet1=TMath::Min(fabs(W_phi-jetAK8puppi_phi),2*Pi-fabs(W_phi-jetAK8puppi_phi));
        delPhiWjet2=TMath::Min(fabs(W_phi-jetAK8puppi_phi_2),2*Pi-fabs(W_phi-jetAK8puppi_phi_2));
        delPhiWjetclose=(fabs(delPhiWjet1)<fabs(delPhiWjet2))?delPhiWjet1:delPhiWjet2;
        delPhiWjetfar=(fabs(delPhiWjet1)>fabs(delPhiWjet2))?delPhiWjet1:delPhiWjet2;
        if(Nj8==2)
            ST=W_pt+jet_pt_puppi+jet_pt_puppi_2;
        if(Nj8==1)
            ST=W_pt+jet_pt_puppi;
        if(Nj8>=3)
            ST=W_pt+jet_pt_puppi+jet_pt_puppi_2+jet_pt_puppi_3;

        for(int i=0;i<8;i++){
            deltaRjet1ak4[i]=-99;deltaRjet1ak4_1[i]=-99;deltaRphiak4[i]=-99;}
        for(Int_t i=0;i<8;i++){
        deltaRjet1ak4[i]=sqrt(pow(fabs(ak4jet_eta[i]-jetAK8puppi_eta),2)+pow(TMath::Min(fabs(ak4jet_phi[i]-jetAK8puppi_phi),2*Pi-fabs(ak4jet_phi[i]-jetAK8puppi_phi)),2));
        deltaRjet1ak4_1[i]=deltaR(ak4jet_eta[i],ak4jet_phi[i],jetAK8puppi_eta,jetAK8puppi_phi);

        }
        Double_t tmp_deltaRjet1ak4=999,tmp_deltaRjet1ak4_1=999;
        Int_t tmp_mini=-99,tmp_mini_1=-99;
        for(Int_t i=0;i<8;i++){
            if((tmp_deltaRjet1ak4>=deltaRjet1ak4[i])&&deltaRjet1ak4[i]<1.6&&deltaRjet1ak4[i]>0.8&&deltaRjet1ak4[i]>0&&ak4jet_pt[i]>0)
            {tmp_mini=i;
                tmp_deltaRjet1ak4=deltaRjet1ak4[i];}
            if((tmp_deltaRjet1ak4_1>=deltaRjet1ak4_1[i])&&deltaRjet1ak4_1[i]<1.6&&deltaRjet1ak4_1[i]>0.8&&ak4jet_pt[i]>0)
            {tmp_mini_1=i;
                tmp_deltaRjet1ak4_1=deltaRjet1ak4[i];}


        }
        jet_mass_puppi_casea=-99;        jet_mass_puppi_casea_1=-99;

        //if(tmp_mini<0) cout<<tmp_mini<<"  "<<jet_mass_puppi_casea<<endl;
        if(tmp_mini>0){
            TLorentzVector ak8jet1,ak4closejet1;
            ak8jet1.SetPtEtaPhiM(jet_pt_puppi,jetAK8puppi_eta,jetAK8puppi_phi,jet_mass_puppi_corr);
            ak4closejet1.SetPtEtaPhiE(ak4jet_pt[tmp_mini],ak4jet_eta[tmp_mini],ak4jet_phi[tmp_mini],ak4jet_e[tmp_mini]);
            jet_mass_puppi_casea = (ak8jet1+ak4closejet1).Mag();

        }
        //cout<<tmp_mini<<"  "<<jet_mass_puppi_casea<<"  "<<jet_mass_puppi_corr<<endl;
        numak4outsideak8=0;
        Double_t tmpak4jet_pt[8],tmpak4jet_eta[8],tmpak4jet_phi[8],tmpak4jet_e[8];
        Double_t tmpdRjet1ak4[8]={999,999,999,999,999,999,999,999};
        for(Int_t i=0;i<8;i++){
            if(deltaRjet1ak4_1[i]<1.6&&deltaRjet1ak4_1[i]>0.8&&ak4jet_pt[i]>0)
            {
                tmpak4jet_pt[numak4outsideak8]=ak4jet_pt[i];
                tmpak4jet_eta[numak4outsideak8]=ak4jet_eta[i];
                tmpak4jet_phi[numak4outsideak8]=ak4jet_phi[i];
                tmpak4jet_e[numak4outsideak8]=ak4jet_e[i];
                tmpdRjet1ak4[numak4outsideak8]=deltaRjet1ak4_1[i];
            numak4outsideak8++;
            }
        }
        if(numak4outsideak8==0) jet_mass_puppi_casea_1=-99;

        if(numak4outsideak8==1)
        {   TLorentzVector ak8jet1_1,ak4closejet1_1;
            ak8jet1_1.SetPtEtaPhiM(jet_pt_puppi,jetAK8puppi_eta,jetAK8puppi_phi,jet_mass_puppi_corr);
            ak4closejet1_1.SetPtEtaPhiE(tmpak4jet_pt[0],tmpak4jet_eta[0],tmpak4jet_phi[0],tmpak4jet_e[0]);
            jet_mass_puppi_casea_1 = (ak8jet1_1+ak4closejet1_1).Mag();
}
        if(numak4outsideak8>=2){
            Int_t *indexx=new Int_t[8];
            TMath::Sort(8,tmpdRjet1ak4,indexx,0);
            TLorentzVector ak8jet1_1,ak4closejet1_1,ak4closejet1_2;
            ak8jet1_1.SetPtEtaPhiM(jet_pt_puppi,jetAK8puppi_eta,jetAK8puppi_phi,jet_mass_puppi_corr);
            ak4closejet1_1.SetPtEtaPhiE(ak4jet_pt[indexx[0]],ak4jet_eta[indexx[0]],ak4jet_phi[indexx[0]],ak4jet_e[indexx[0]]);
            ak4closejet1_2.SetPtEtaPhiE(ak4jet_pt[indexx[1]],ak4jet_eta[indexx[1]],ak4jet_phi[indexx[1]],ak4jet_e[indexx[1]]);
            
            jet_mass_puppi_casea_1 = (ak8jet1_1+ak4closejet1_1+ak4closejet1_2).Mag();
            
        }
        //cout<<numak4outsideak8<<"    "<<jet_mass_puppi_casea_1<<endl;
        /*if(Nj4>=3){
            TLorentzVector ak8jet1_1,ak4closejet1_1,ak4closejet1_2;
            ak8jet1_1.SetPtEtaPhiM(jet_pt_puppi,jetAK8puppi_eta,jetAK8puppi_phi,jet_mass_puppi_corr);
            ak4closejet1_1.SetPtEtaPhiE(ak4jet_pt[1],ak4jet_eta[1],ak4jet_phi[1],ak4jet_e[1]);
            ak4closejet1_2.SetPtEtaPhiE(ak4jet_pt[2],ak4jet_eta[2],ak4jet_phi[2],ak4jet_e[2]);

            jet_mass_puppi_casea_1 = (ak8jet1_1+ak4closejet1_1+ak4closejet1_2).Mag();
            
        }
        if(Nj4==2){
            TLorentzVector ak8jet1_1,ak4closejet1_1;
            ak8jet1_1.SetPtEtaPhiM(jet_pt_puppi,jetAK8puppi_eta,jetAK8puppi_phi,jet_mass_puppi_corr);
            ak4closejet1_1.SetPtEtaPhiE(ak4jet_pt[1],ak4jet_eta[1],ak4jet_phi[1],ak4jet_e[1]);
            jet_mass_puppi_casea_1 = (ak8jet1_1+ak4closejet1_1).Mag();
            
        }*/


 
        for(Int_t i=0;i<8;i++){
            //deltaRphiak4[i]=sqrt(pow(fabs(-ak4jet_eta[i]-jetAK8puppi_eta),2)+pow(TMath::Min(fabs(ak4jet_phi[i]-(-Pi+jetAK8puppi_phi)),2*Pi-fabs(Pi+ak4jet_phi[i]-jetAK8puppi_phi)),2));
            deltaRphiak4[i]=deltaR(-jetAK8puppi_eta,((jetAK8puppi_phi>0)?(-Pi+jetAK8puppi_phi):(Pi+jetAK8puppi_phi)),ak4jet_eta[i],ak4jet_phi[i]);
        }
        numak4insidephi=0;
        Double_t tmp1ak4jet_pt[8],tmp1ak4jet_eta[8],tmp1ak4jet_phi[8],tmp1ak4jet_e[8];
        Double_t tmpdRphiak4[8]={-999,-999,-999,-999,-999,-999,-999,-999};
        for(Int_t i=0;i<8;i++){
            if(deltaRphiak4[i]<1.6&&ak4jet_pt[i]>0)
            {
                tmp1ak4jet_pt[numak4outsideak8]=ak4jet_pt[i];
                tmp1ak4jet_eta[numak4outsideak8]=ak4jet_eta[i];
                tmp1ak4jet_phi[numak4outsideak8]=ak4jet_phi[i];
                tmp1ak4jet_e[numak4outsideak8]=ak4jet_e[i];
                tmpdRphiak4[numak4outsideak8]=tmpdRphiak4[i];
                numak4insidephi++;
            }
        }
        Int_t *indexa=new Int_t[8];
        TMath::Sort(8,deltaRphiak4,indexa,0);

        //if((status_1==4 &&(status_2!=4||status_3!=4) )&&(! (status_2==4&&status_3==4)))
        //cout<<jetAK8puppi_phi<<"   "<<((jetAK8puppi_phi>0)?(-Pi+jetAK8puppi_phi):(Pi+jetAK8puppi_phi))<<"  "<<deltaRphiak4[indexa[0]]<<" "<<deltaR(jetAK8puppi_eta,jetAK8puppi_phi,W_eta,W_phi)<<" "<<deltaR(jetAK8puppi_eta,jetAK8puppi_phi,ak4jet_eta[0],ak4jet_phi[0])<<" "<<numak4insidephi<<endl;
        if(numak4insidephi==0)
        {                    jet_mass_puppi_caseb = massVlepJEC;
        }
       
        if(numak4insidephi==1)
        {   TLorentzVector phireco_1,ak4closephi_1;
            phireco_1.SetPtEtaPhiM(W_pt,W_eta,W_phi,massVlepJEC);
            ak4closephi_1.SetPtEtaPhiE(tmp1ak4jet_pt[0],tmp1ak4jet_eta[0],tmp1ak4jet_phi[0],tmp1ak4jet_e[0]);
            jet_mass_puppi_caseb = (phireco_1+ak4closephi_1).Mag();
        }
        if(numak4insidephi>=2){
            Int_t *index=new Int_t[8];
            TMath::Sort(8,tmpdRphiak4,index,1);
            TLorentzVector phireco_1,ak4closephi_1,ak4closephi_2;
            phireco_1.SetPtEtaPhiM(W_pt,W_eta,W_phi,massVlepJEC);
            ak4closephi_1.SetPtEtaPhiE(tmp1ak4jet_pt[index[0]],tmp1ak4jet_eta[index[0]],tmp1ak4jet_phi[index[0]],tmp1ak4jet_e[index[0]]);
            ak4closephi_2.SetPtEtaPhiE(tmp1ak4jet_pt[index[1]],tmp1ak4jet_eta[index[1]],tmp1ak4jet_phi[index[1]],tmp1ak4jet_e[index[1]]);

            jet_mass_puppi_caseb = (phireco_1+ak4closephi_1+ak4closephi_2).Mag();
        }
        
        Mj_min=-99; Mj_max=-99;pt_min=-99; pt_max=-99; t21_max=-99; t31_max=-99; t41_max=-99;  t21_min=-99; t31_min=-99; t41_min=-99; PTw2divPTw1=-99; PTw3divPTw1=-99; PTw3divPTw2=-99;phi_max=-99;phi_min=-99;eta_max=-99;eta_min=-99;t32_max=-99;t42_max=-99;t43_max=-99;t32_min=-99;t42_min=-99;t43_min=-99; t1_max=-99;t2_max=-99;t3_max=-99;t4_max=-99;t1_mid=-99;t2_mid=-99;t3_mid=-99;t4_mid=-99;t1_min=-99;t2_min=-99;t3_min=-99;t4_min=-99;
        Mj_mid=-99;pt_mid=-99; t21_mid=-99; t31_mid=-99; t41_mid=-99;phi_mid=-99;eta_mid=-99;t32_mid=-99;t42_mid=-99;t43_mid=-99;
        dnnt_max=-99;dnnw_max=-99;dnnh_max=-99;dnndect_max=-99;dnndecw_max=-99;dnndech_max=-99;
        dnnt_mid=-99;dnnw_mid=-99;dnnh_mid=-99;dnndect_mid=-99;dnndecw_mid=-99;dnndech_mid=-99;
        dnnt_min=-99;dnnw_min=-99;dnnh_min=-99;dnndect_min=-99;dnndecw_min=-99;dnndech_min=-99;
        dnnz_max=-99;dnnzbb_max=-99;dnnhbb_max=-99;dnndecz_max=-99;dnndeczbb_max=-99;dnndechbb_max=-99;dnndecbb_max=-99;dnndecbbnog_max=-99;dnndeccc_max=-99;dnndecccnog_max=-99;
        dnnz_mid=-99;dnnzbb_mid=-99;dnnhbb_mid=-99;dnndecz_mid=-99;dnndeczbb_mid=-99;dnndechbb_mid=-99;dnndecbb_mid=-99;dnndecbbnog_mid=-99;dnndeccc_mid=-99;dnndecccnog_mid=-99;
        dnnz_min=-99;dnnzbb_min=-99;dnnhbb_min=-99;dnndecz_min=-99;dnndeczbb_min=-99;dnndechbb_min=-99;dnndecbb_min=-99;dnndecbbnog_min=-99;dnndeccc_min=-99;dnndecccnog_min=-99;
        dnndecrawqcd_max=-99;dnndecraww_max=-99;dnndecrawz_max=-99;dnndecrawzbb_max=-99;dnndecrawhbb_max=-99;dnndecrawt_max=-99;dnndecrawh4q_max=-99;
        dnndecrawqcd_mid=-99;dnndecraww_mid=-99;dnndecrawz_mid=-99;dnndecrawzbb_mid=-99;dnndecrawhbb_mid=-99;dnndecrawt_mid=-99;dnndecrawh4q_mid=-99;
        dnndecrawqcd_min=-99;dnndecraww_min=-99;dnndecrawz_min=-99;dnndecrawzbb_min=-99;dnndecrawhbb_min=-99;dnndecrawt_min=-99;dnndecrawh4q_min=-99;
        Double_t Mj_sort[3]={jet_mass_puppi_corr,jet_mass_puppi_corr_2,jet_mass_puppi_corr_3},
            pt_sort[3]={jet_pt_puppi,jet_pt_puppi_2,jet_pt_puppi_3},
            eta_sort[3]={jetAK8puppi_eta,jetAK8puppi_eta_2,jetAK8puppi_eta_3},
            phi_sort[3]={jetAK8puppi_phi,jetAK8puppi_phi_2,jetAK8puppi_phi_3},
        t1_sort[3]={jetAK8puppi_tau1,jetAK8puppi_tau1_2,jetAK8puppi_tau1_3},
        t2_sort[3]={jetAK8puppi_tau2,jetAK8puppi_tau2_2,jetAK8puppi_tau2_3},
        t3_sort[3]={jetAK8puppi_tau3,jetAK8puppi_tau3_2,jetAK8puppi_tau3_3},
        t4_sort[3]={jetAK8puppi_tau4,jetAK8puppi_tau4_2,jetAK8puppi_tau4_3},
            t21_sort[3]={jetAK8puppi_tau21,jetAK8puppi_tau21_2,jetAK8puppi_tau21_3},
            t31_sort[3]={jetAK8puppi_tau31,jetAK8puppi_tau31_2,jetAK8puppi_tau31_3},
            t41_sort[3]={jetAK8puppi_tau41,jetAK8puppi_tau41_2,jetAK8puppi_tau41_3},
            t32_sort[3]={jetAK8puppi_tau32,jetAK8puppi_tau32_2,jetAK8puppi_tau32_3},
            t42_sort[3]={jetAK8puppi_tau42,jetAK8puppi_tau42_2,jetAK8puppi_tau42_3},
            t43_sort[3]={jetAK8puppi_tau43,jetAK8puppi_tau43_2,jetAK8puppi_tau43_3},
            dnnt_sort[3]={jetAK8puppi_dnnTop,jetAK8puppi_dnnTop_2,jetAK8puppi_dnnTop_3},
            dnnw_sort[3]={jetAK8puppi_dnnW,jetAK8puppi_dnnW_2,jetAK8puppi_dnnW_3},
            dnnh_sort[3]={jetAK8puppi_dnnH4q,jetAK8puppi_dnnH4q_2,jetAK8puppi_dnnH4q_3},
            dnndect_sort[3]={jetAK8puppi_dnnDecorrTop,jetAK8puppi_dnnDecorrTop_2,jetAK8puppi_dnnDecorrTop_3},
            dnndecw_sort[3]={jetAK8puppi_dnnDecorrW,jetAK8puppi_dnnDecorrW_2,jetAK8puppi_dnnDecorrW_3},
            dnndech_sort[3]={jetAK8puppi_dnnDecorrH4q,jetAK8puppi_dnnDecorrH4q_2,jetAK8puppi_dnnDecorrH4q_3},
        dnnz_sort[3]={jetAK8puppi_dnnZ,jetAK8puppi_dnnZ_2,jetAK8puppi_dnnZ_3},
        dnnzbb_sort[3]={jetAK8puppi_dnnZbb,jetAK8puppi_dnnZbb_2,jetAK8puppi_dnnZbb_3},
        dnnhbb_sort[3]={jetAK8puppi_dnnHbb,jetAK8puppi_dnnHbb_2,jetAK8puppi_dnnHbb_3},
        dnndecz_sort[3]={jetAK8puppi_dnnDecorrZ,jetAK8puppi_dnnDecorrZ_2,jetAK8puppi_dnnDecorrZ_3},
        dnndeczbb_sort[3]={jetAK8puppi_dnnDecorrZbb,jetAK8puppi_dnnDecorrZbb_2,jetAK8puppi_dnnDecorrZbb_3},
        dnndechbb_sort[3]={jetAK8puppi_dnnDecorrHbb,jetAK8puppi_dnnDecorrHbb_2,jetAK8puppi_dnnDecorrHbb_3},
        dnndecbb_sort[3]={jetAK8puppi_dnnDecorrbb,jetAK8puppi_dnnDecorrbb_2,jetAK8puppi_dnnDecorrbb_3},
        dnndecbbnog_sort[3]={jetAK8puppi_dnnDecorrbbnog,jetAK8puppi_dnnDecorrbbnog_2,jetAK8puppi_dnnDecorrbbnog_3},
        dnndeccc_sort[3]={jetAK8puppi_dnnDecorrcc,jetAK8puppi_dnnDecorrcc_2,jetAK8puppi_dnnDecorrcc_3},
        dnndecccnog_sort[3]={jetAK8puppi_dnnDecorrccnog,jetAK8puppi_dnnDecorrccnog_2,jetAK8puppi_dnnDecorrccnog_3},
        dnndecrawqcd_sort[3]={jetAK8puppi_dnnDecorrqcd,jetAK8puppi_dnnDecorrqcd_2,jetAK8puppi_dnnDecorrqcd_3},
        dnndecrawt_sort[3]={jetAK8puppi_dnnDecorrtop,jetAK8puppi_dnnDecorrtop_2,jetAK8puppi_dnnDecorrtop_3},
        dnndecrawz_sort[3]={jetAK8puppi_dnnDecorrz,jetAK8puppi_dnnDecorrz_2,jetAK8puppi_dnnDecorrz_3},
        dnndecraww_sort[3]={jetAK8puppi_dnnDecorrw,jetAK8puppi_dnnDecorrw_2,jetAK8puppi_dnnDecorrw_3},
        dnndecrawzbb_sort[3]={jetAK8puppi_dnnDecorrzbb,jetAK8puppi_dnnDecorrzbb_2,jetAK8puppi_dnnDecorrzbb_3},
        dnndecrawhbb_sort[3]={jetAK8puppi_dnnDecorrhbb,jetAK8puppi_dnnDecorrhbb_2,jetAK8puppi_dnnDecorrhbb_3},
        dnndecrawh4q_sort[3]={jetAK8puppi_dnnDecorrh4q,jetAK8puppi_dnnDecorrh4q_2,jetAK8puppi_dnnDecorrh4q_3}
        ;
        
        if(Nj8==1){
            Mj_max=Mj_sort[0];pt_max=pt_sort[0];eta_max=eta_sort[0];phi_max=phi_sort[0];t21_max=t21_sort[0];t31_max=t31_sort[0];t41_max=t41_sort[0];t32_max=t32_sort[0];t42_max=t42_sort[0];t43_max=t43_sort[0];
            Mj_min=Mj_sort[0];pt_min=pt_sort[0];eta_min=eta_sort[0];phi_min=phi_sort[0];t21_min=t21_sort[0];t31_min=t31_sort[0];t41_min=t41_sort[0];t32_min=t32_sort[0];t42_min=t42_sort[0];t43_min=t43_sort[0];
            dnnt_max=dnnt_sort[0];dnnw_max=dnnw_sort[0];dnnh_max=dnnh_sort[0];dnndect_max=dnndect_sort[0];dnndecw_max=dnndecw_sort[0];dnndech_max=dnndech_sort[0];
            dnnt_min=dnnt_sort[0];dnnw_min=dnnw_sort[0];dnnh_min=dnnh_sort[0];dnndect_min=dnndect_sort[0];dnndecw_min=dnndecw_sort[0];dnndech_min=dnndech_sort[0];
            dnnz_max=dnnz_sort[0];dnnzbb_max=dnnzbb_sort[0];dnnhbb_max=dnnhbb_sort[0];dnndecz_max=dnndecz_sort[0];dnndeczbb_max=dnndeczbb_sort[0];dnndechbb_max=dnndechbb_sort[0];dnndecbb_max=dnndecbb_sort[0];dnndecbbnog_max=dnndecbbnog_sort[0];dnndeccc_max=dnndeccc_sort[0];dnndecccnog_max=dnndecccnog_sort[0];
            dnnz_min=dnnz_sort[0];dnnzbb_min=dnnzbb_sort[0];dnnhbb_min=dnnhbb_sort[0];dnndecz_min=dnndecz_sort[0];dnndeczbb_min=dnndeczbb_sort[0];dnndechbb_min=dnndechbb_sort[0];dnndecbb_min=dnndecbb_sort[0];dnndecbbnog_min=dnndecbbnog_sort[0];dnndeccc_min=dnndeccc_sort[0];dnndecccnog_min=dnndecccnog_sort[0];
            dnndecrawqcd_max=dnndecrawqcd_sort[0];dnndecrawt_max=dnndecrawt_sort[0];dnndecraww_max=dnndecraww_sort[0];dnndecrawz_max=dnndecrawz_sort[0];dnndecrawzbb_max=dnndecrawzbb_sort[0];dnndecrawhbb_max=dnndecrawhbb_sort[0];dnndecrawh4q_max=dnndecrawh4q_sort[0];
            dnndecrawqcd_min=dnndecrawqcd_sort[0];dnndecrawt_min=dnndecrawt_sort[0];dnndecraww_min=dnndecraww_sort[0];dnndecrawz_min=dnndecrawz_sort[0];dnndecrawzbb_min=dnndecrawzbb_sort[0];dnndecrawhbb_min=dnndecrawhbb_sort[0];dnndecrawh4q_min=dnndecrawh4q_sort[0];
            t1_max=t1_sort[0];t2_max=t2_sort[0];t3_max=t3_sort[0];t4_max=t4_sort[0];
            t1_min=t1_sort[0];t2_min=t2_sort[0];t3_min=t3_sort[0];t4_min=t4_sort[0];
        }
        if(Nj8==2){
            Int_t *indexmass=new Int_t[3];TMath::Sort(2,Mj_sort,indexmass,1);
            Mj_max=Mj_sort[indexmass[0]];pt_max=pt_sort[indexmass[0]];eta_max=eta_sort[indexmass[0]];phi_max=phi_sort[indexmass[0]];t21_max=t21_sort[indexmass[0]];t31_max=t31_sort[indexmass[0]];t41_max=t41_sort[indexmass[0]];t32_max=t32_sort[indexmass[0]];t42_max=t42_sort[indexmass[0]];t43_max=t43_sort[indexmass[0]];
            Mj_min=Mj_sort[indexmass[1]];pt_min=pt_sort[indexmass[1]];eta_min=eta_sort[indexmass[1]];phi_min=phi_sort[indexmass[1]];t21_min=t21_sort[indexmass[1]];t31_min=t31_sort[indexmass[1]];t41_min=t41_sort[indexmass[1]];t32_min=t32_sort[indexmass[1]];t42_min=t42_sort[indexmass[1]];t43_min=t43_sort[indexmass[1]];
            dnnt_max=dnnt_sort[indexmass[0]];dnnw_max=dnnw_sort[indexmass[0]];dnnh_max=dnnh_sort[indexmass[0]];dnndect_max=dnndect_sort[indexmass[0]];dnndecw_max=dnndecw_sort[indexmass[0]];dnndech_max=dnndech_sort[indexmass[0]];
            dnnt_min=dnnt_sort[indexmass[1]];dnnw_min=dnnw_sort[indexmass[1]];dnnh_min=dnnh_sort[indexmass[1]];dnndect_min=dnndect_sort[indexmass[1]];dnndecw_min=dnndecw_sort[indexmass[1]];dnndech_min=dnndech_sort[indexmass[1]];
            dnnz_max=dnnz_sort[indexmass[0]];dnnzbb_max=dnnzbb_sort[indexmass[0]];dnnhbb_max=dnnhbb_sort[indexmass[0]];dnndecz_max=dnndecz_sort[indexmass[0]];dnndeczbb_max=dnndeczbb_sort[indexmass[0]];dnndechbb_max=dnndechbb_sort[indexmass[0]];dnndecbb_max=dnndecbb_sort[indexmass[0]];dnndecbbnog_max=dnndecbbnog_sort[indexmass[0]];dnndeccc_max=dnndeccc_sort[indexmass[0]];dnndecccnog_max=dnndecccnog_sort[indexmass[0]];
            dnnz_min=dnnz_sort[indexmass[1]];dnnzbb_min=dnnzbb_sort[indexmass[1]];dnnhbb_min=dnnhbb_sort[indexmass[1]];dnndecz_min=dnndecz_sort[indexmass[1]];dnndeczbb_min=dnndeczbb_sort[indexmass[1]];dnndechbb_min=dnndechbb_sort[indexmass[1]];dnndecbb_min=dnndecbb_sort[indexmass[1]];dnndecbbnog_min=dnndecbbnog_sort[indexmass[1]];dnndeccc_min=dnndeccc_sort[indexmass[1]];dnndecccnog_min=dnndecccnog_sort[indexmass[1]];
            dnndecrawqcd_max=dnndecrawqcd_sort[indexmass[0]];dnndecrawt_max=dnndecrawt_sort[indexmass[0]];dnndecraww_max=dnndecraww_sort[indexmass[0]];dnndecrawz_max=dnndecrawz_sort[indexmass[0]];dnndecrawzbb_max=dnndecrawzbb_sort[indexmass[0]];dnndecrawhbb_max=dnndecrawhbb_sort[indexmass[0]];dnndecrawh4q_max=dnndecrawh4q_sort[indexmass[0]];
            dnndecrawqcd_min=dnndecrawqcd_sort[indexmass[1]];dnndecrawt_min=dnndecrawt_sort[indexmass[1]];dnndecraww_min=dnndecraww_sort[indexmass[1]];dnndecrawz_min=dnndecrawz_sort[indexmass[1]];dnndecrawzbb_min=dnndecrawzbb_sort[indexmass[1]];dnndecrawhbb_min=dnndecrawhbb_sort[indexmass[1]];dnndecrawh4q_min=dnndecrawh4q_sort[indexmass[1]];
            t1_max=t1_sort[indexmass[0]];t2_max=t2_sort[indexmass[0]];t3_max=t3_sort[indexmass[0]];t4_max=t4_sort[indexmass[0]];
            t1_min=t1_sort[indexmass[1]];t2_min=t2_sort[indexmass[1]];t3_min=t3_sort[indexmass[1]];t4_min=t4_sort[indexmass[1]];
        }
        if(Nj8>=3){
            Int_t *indexmass=new Int_t[3];TMath::Sort(3,Mj_sort,indexmass,1);
            Mj_max=Mj_sort[indexmass[0]];pt_max=pt_sort[indexmass[0]];eta_max=eta_sort[indexmass[0]];phi_max=phi_sort[indexmass[0]];t21_max=t21_sort[indexmass[0]];t31_max=t31_sort[indexmass[0]];t41_max=t41_sort[indexmass[0]];t32_max=t32_sort[indexmass[0]];t42_max=t42_sort[indexmass[0]];t43_max=t43_sort[indexmass[0]];
            Mj_mid=Mj_sort[indexmass[1]];pt_mid=pt_sort[indexmass[1]];eta_mid=eta_sort[indexmass[1]];phi_mid=phi_sort[indexmass[1]];t21_mid=t21_sort[indexmass[1]];t31_mid=t31_sort[indexmass[1]];t41_mid=t41_sort[indexmass[1]];t32_mid=t32_sort[indexmass[1]];t42_mid=t42_sort[indexmass[1]];t43_mid=t43_sort[indexmass[1]];
            Mj_min=Mj_sort[indexmass[2]];pt_min=pt_sort[indexmass[2]];eta_min=eta_sort[indexmass[2]];phi_min=phi_sort[indexmass[2]];t21_min=t21_sort[indexmass[2]];t31_min=t31_sort[indexmass[2]];t41_min=t41_sort[indexmass[2]];t32_min=t32_sort[indexmass[2]];t42_min=t42_sort[indexmass[2]];t43_min=t43_sort[indexmass[2]];
            dnnt_max=dnnt_sort[indexmass[0]];dnnw_max=dnnw_sort[indexmass[0]];dnnh_max=dnnh_sort[indexmass[0]];dnndect_max=dnndect_sort[indexmass[0]];dnndecw_max=dnndecw_sort[indexmass[0]];dnndech_max=dnndech_sort[indexmass[0]];
            dnnt_mid=dnnt_sort[indexmass[1]];dnnw_mid=dnnw_sort[indexmass[1]];dnnh_mid=dnnh_sort[indexmass[1]];dnndect_mid=dnndect_sort[indexmass[1]];dnndecw_mid=dnndecw_sort[indexmass[1]];dnndech_mid=dnndech_sort[indexmass[1]];
            dnnt_min=dnnt_sort[indexmass[2]];dnnw_min=dnnw_sort[indexmass[2]];dnnh_min=dnnh_sort[indexmass[2]];dnndect_min=dnndect_sort[indexmass[2]];dnndecw_min=dnndecw_sort[indexmass[2]];dnndech_min=dnndech_sort[indexmass[2]];
            dnnz_max=dnnz_sort[indexmass[0]];dnnzbb_max=dnnzbb_sort[indexmass[0]];dnnhbb_max=dnnhbb_sort[indexmass[0]];dnndecz_max=dnndecz_sort[indexmass[0]];dnndeczbb_max=dnndeczbb_sort[indexmass[0]];dnndechbb_max=dnndechbb_sort[indexmass[0]];dnndecbb_max=dnndecbb_sort[indexmass[0]];dnndecbbnog_max=dnndecbbnog_sort[indexmass[0]];dnndeccc_max=dnndeccc_sort[indexmass[0]];dnndecccnog_max=dnndecccnog_sort[indexmass[0]];
            dnnz_mid=dnnz_sort[indexmass[1]];dnnzbb_mid=dnnzbb_sort[indexmass[1]];dnnhbb_mid=dnnhbb_sort[indexmass[1]];dnndecz_mid=dnndecz_sort[indexmass[1]];dnndeczbb_mid=dnndeczbb_sort[indexmass[1]];dnndechbb_mid=dnndechbb_sort[indexmass[1]];dnndecbb_mid=dnndecbb_sort[indexmass[1]];dnndecbbnog_mid=dnndecbbnog_sort[indexmass[1]];dnndeccc_mid=dnndeccc_sort[indexmass[1]];dnndecccnog_mid=dnndecccnog_sort[indexmass[1]];
            dnnz_min=dnnz_sort[indexmass[2]];dnnzbb_min=dnnzbb_sort[indexmass[2]];dnnhbb_min=dnnhbb_sort[indexmass[2]];dnndecz_min=dnndecz_sort[indexmass[2]];dnndeczbb_min=dnndeczbb_sort[indexmass[2]];dnndechbb_min=dnndechbb_sort[indexmass[2]];dnndecbb_min=dnndecbb_sort[indexmass[2]];dnndecbbnog_min=dnndecbbnog_sort[indexmass[2]];dnndeccc_min=dnndeccc_sort[indexmass[2]];dnndecccnog_min=dnndecccnog_sort[indexmass[2]];
            dnndecrawqcd_max=dnndecrawqcd_sort[indexmass[0]];dnndecrawt_max=dnndecrawt_sort[indexmass[0]];dnndecraww_max=dnndecraww_sort[indexmass[0]];dnndecrawz_max=dnndecrawz_sort[indexmass[0]];dnndecrawzbb_max=dnndecrawzbb_sort[indexmass[0]];dnndecrawhbb_max=dnndecrawhbb_sort[indexmass[0]];dnndecrawh4q_max=dnndecrawh4q_sort[indexmass[0]];
            dnndecrawqcd_mid=dnndecrawqcd_sort[indexmass[1]];dnndecrawt_mid=dnndecrawt_sort[indexmass[1]];dnndecraww_mid=dnndecraww_sort[indexmass[1]];dnndecrawz_mid=dnndecrawz_sort[indexmass[1]];dnndecrawzbb_mid=dnndecrawzbb_sort[indexmass[1]];dnndecrawhbb_mid=dnndecrawhbb_sort[indexmass[1]];dnndecrawh4q_mid=dnndecrawh4q_sort[indexmass[1]];
            dnndecrawqcd_min=dnndecrawqcd_sort[indexmass[2]];dnndecrawt_min=dnndecrawt_sort[indexmass[2]];dnndecraww_min=dnndecraww_sort[indexmass[2]];dnndecrawz_min=dnndecrawz_sort[indexmass[2]];dnndecrawzbb_min=dnndecrawzbb_sort[indexmass[2]];dnndecrawhbb_min=dnndecrawhbb_sort[indexmass[2]];dnndecrawh4q_min=dnndecrawh4q_sort[indexmass[2]];
            t1_max=t1_sort[indexmass[0]];t2_max=t2_sort[indexmass[0]];t3_max=t3_sort[indexmass[0]];t4_max=t4_sort[indexmass[0]];
            t1_mid=t1_sort[indexmass[1]];t2_mid=t2_sort[indexmass[1]];t3_mid=t3_sort[indexmass[1]];t4_mid=t4_sort[indexmass[1]];
            t1_min=t1_sort[indexmass[2]];t2_min=t2_sort[indexmass[2]];t3_min=t3_sort[indexmass[2]];t4_min=t4_sort[indexmass[2]];
        }

  
        Int_t *indexx11=new Int_t[3];
        double ptwdivw[3]={-99,-99,-99};
        ptwdivw[0]=W_pt;ptwdivw[1]=jet_pt_puppi;ptwdivw[2]=jet_pt_puppi_2;
        TMath::Sort(3,ptwdivw,indexx11,1);
        PTw2divPTw1=ptwdivw[indexx11[1]]/ptwdivw[indexx11[0]];
        PTw3divPTw1=ptwdivw[indexx11[2]]/ptwdivw[indexx11[0]];
        PTw3divPTw2=ptwdivw[indexx11[2]]/ptwdivw[indexx11[1]];
        PTw[0]=ptwdivw[indexx11[0]];PTw[1]=ptwdivw[indexx11[1]];PTw[2]=ptwdivw[indexx11[2]];
        if(Nj8==1)
        {PTw3divPTw1=-99;PTw3divPTw2=-99;}
        
        double deltaRjet1ak4p[8]={-99,-99,-99,-99,-99,-99,-99,-99};
        double deltaRjet2ak4p[8]={-99,-99,-99,-99,-99,-99,-99,-99};
        double deltaRjet3ak4p[8]={-99,-99,-99,-99,-99,-99,-99,-99};

        TLorentzVector tempak4,ak8jet1_2;
        ak8jet1_2.SetPtEtaPhiM(jet_pt_puppi,jetAK8puppi_eta,jetAK8puppi_phi,jet_mass_puppi_corr);
        
        Nj4ex=0;
        Nj4in=0;
            ak4expt1=-99,ak4exeta1=-99,ak4exphi1=-99,ak4exe1=-99;
            ak4expt2=-99,ak4exeta2=-99,ak4exphi2=-99,ak4exe2=-99;
            
            Nj4ex_dR0p5=0;
            Nj4in_dR0p5=0;
            ak4expt1_dR0p5=-99,ak4exeta1_dR0p5=-99,ak4exphi1_dR0p5=-99,ak4exe1_dR0p5=-99;
            ak4expt2_dR0p5=-99,ak4exeta2_dR0p5=-99,ak4exphi2_dR0p5=-99,ak4exe2_dR0p5=-99;

        
        Nj4ex_dR0p8=0;
        Nj4in_dR0p8=0;
        ak4expt1_dR0p8=-99,ak4exeta1_dR0p8=-99,ak4exphi1_dR0p8=-99,ak4exe1_dR0p8=-99;
        ak4expt2_dR0p8=-99,ak4exeta2_dR0p8=-99,ak4exphi2_dR0p8=-99,ak4exe2_dR0p8=-99;
        
        Nj4ex_dR1p0=0;
        Nj4in_dR1p0=0;
        ak4expt1_dR1p0=-99,ak4exeta1_dR1p0=-99,ak4exphi1_dR1p0=-99,ak4exe1_dR1p0=-99;
        ak4expt2_dR1p0=-99,ak4exeta2_dR1p0=-99,ak4exphi2_dR1p0=-99,ak4exe2_dR1p0=-99;
        
        Nj4ex_dR1p2=0;
        Nj4in_dR1p2=0;
        ak4expt1_dR1p2=-99,ak4exeta1_dR1p2=-99,ak4exphi1_dR1p2=-99,ak4exe1_dR1p2=-99;
        ak4expt2_dR1p2=-99,ak4exeta2_dR1p2=-99,ak4exphi2_dR1p2=-99,ak4exe2_dR1p2=-99;

            double dR0p8=0.8;
        
        for(Int_t ii=0; ii<8; ii++)  {
            if(ak4jet_pt[ii]>30 && fabs(ak4jet_eta[ii])<2.4 && ak4jet_IDLoose[ii]>0){
                deltaRjet1ak4p[ii]=sqrt(pow(fabs(ak4jet_eta[ii]-jetAK8puppi_eta),2)+pow(TMath::Min(fabs(ak4jet_phi[ii]-jetAK8puppi_phi),2*Pi-fabs(ak4jet_phi[ii]-jetAK8puppi_phi)),2));
                if(Nj8>=2)
                deltaRjet2ak4p[ii]=sqrt(pow(fabs(ak4jet_eta[ii]-jetAK8puppi_eta_2),2)+pow(TMath::Min(fabs(ak4jet_phi[ii]-jetAK8puppi_phi_2),2*Pi-fabs(ak4jet_phi[ii]-jetAK8puppi_phi_2)),2));
                if(Nj8>=3)
                deltaRjet3ak4p[ii]=sqrt(pow(fabs(ak4jet_eta[ii]-jetAK8puppi_eta_3),2)+pow(TMath::Min(fabs(ak4jet_phi[ii]-jetAK8puppi_phi_3),2*Pi-fabs(ak4jet_phi[ii]-jetAK8puppi_phi_3)),2));
                

            if((Nj8==1&&(ak4jet_pt[ii]>0)&&deltaRjet1ak4p[ii]<=dR0p8)||(Nj8==2&&(ak4jet_pt[ii]>0)&&(deltaRjet1ak4p[ii]<=dR0p8||deltaRjet2ak4p[ii]<=dR0p8))||(Nj8>=3&&(ak4jet_pt[ii]>0)&&(deltaRjet1ak4p[ii]<=dR0p8||deltaRjet2ak4p[ii]<=dR0p8||deltaRjet3ak4p[ii]<=dR0p8))) {
                ++Nj4in_dR0p8;
                if(Nj4in_dR0p8==1&&(ak4jet_pt[ii]>0)) {
                    //ak4expt1_dR0p8=ak4jet_pt[ii];ak4exeta1_dR0p8=ak4jet_eta[ii];ak4exphi1_dR0p8=ak4jet_phi[ii];ak4exe1_dR0p8=ak4jet_e[ii];ak4expt1_uncorr_dR0p8=ak4jet_pt_uncorr[ii];
                    ak4inpt1_dR0p8=ak4jet_pt[ii];
                    deltaRjet1inak4p_dR0p8=deltaRjet1ak4p[ii];
                    deltaphijet1inak4p_dR0p8=deltaPhi(ak4jet_phi[ii],jetAK8puppi_phi);
                    deltaRinak4mu_dR0p8=deltaR(l_eta,l_phi,ak4jet_eta[ii],ak4jet_phi[ii]);
                    deltaPhiinak4mu_dR0p8=deltaPhi(l_phi,ak4jet_phi[ii]);

            }

            }
            if((Nj8==1&&(ak4jet_pt[ii]>0)&&deltaRjet1ak4p[ii]>dR0p8)||(Nj8==2&&(ak4jet_pt[ii]>0)&&deltaRjet1ak4p[ii]>dR0p8&&deltaRjet2ak4p[ii]>dR0p8)||(Nj8>=3&&(ak4jet_pt[ii]>0)&&deltaRjet1ak4p[ii]>dR0p8&&deltaRjet2ak4p[ii]>dR0p8&&deltaRjet3ak4p[ii]>dR0p8))
            {++Nj4ex_dR0p8;
            if(Nj4ex_dR0p8==1&&(ak4jet_pt[ii]>0)) {
                    ak4expt1_dR0p8=ak4jet_pt[ii];ak4exeta1_dR0p8=ak4jet_eta[ii];ak4exphi1_dR0p8=ak4jet_phi[ii];ak4exe1_dR0p8=ak4jet_e[ii];ak4expt1_uncorr_dR0p8=ak4jet_pt_uncorr[ii];
                deltaRjet1ak4p_dR0p8=deltaRjet1ak4p[ii];
                deltaphijet1ak4p_dR0p8=deltaPhi(ak4exphi1_dR0p8,jetAK8puppi_phi);
                deltaRak4mu_dR0p8=deltaR(l_eta,l_phi,ak4exeta1_dR0p8,ak4exphi1_dR0p8);
                //cout<<" deltaRjet1ak4  "<<deltaR(jetAK8puppi_eta,jetAK8puppi_phi,ak4jet_eta[ii],ak4jet_phi[ii])<<"  deltaRak4lep  "<<deltaR(l_eta,l_phi,ak4jet_eta[ii],ak4jet_phi[ii])<<"  deltaRjet1lep    "<<deltaR(jetAK8puppi_eta,jetAK8puppi_phi,l_eta,l_phi)<<endl;
                /*if(ak4expt1_uncorr_dR0p8>100){
                    cout<<"deltaRjet1ak4p_dR0p8 "<<deltaRjet1ak4p_dR0p8<<"  "<<jet_pt_puppi<<endl;
                for(int i=0;i<numsj1;i++)
                    cout<<"deltaRak4sj "<<i<<" "<<deltaRak4sj[i][ii]<<endl;
                for(int i=0;i<numsj2;i++)
                    cout<<"deltaRak4sj2 "<<i<<" "<<deltaRak4sj2[i][ii]<<endl;
                
                }*/
                }
            if(Nj4ex_dR0p8==2&&(ak4jet_pt[ii]>0)) {
                    ak4expt2_dR0p8=ak4jet_pt[ii];ak4exeta2_dR0p8=ak4jet_eta[ii];ak4exphi2_dR0p8=ak4jet_phi[ii];ak4exe2_dR0p8=ak4jet_e[ii];
                deltaRjet2ak4p_dR0p8=deltaRjet2ak4p[ii];

                }
            }
            double dR1p0=1.0;
            
            if(Nj8==1&&(ak4jet_pt[ii]>0)&&deltaRjet1ak4p[ii]>dR1p0) ++Nj4ex_dR1p0;
            if(Nj8==1&&(ak4jet_pt[ii]>0)&&deltaRjet1ak4p[ii]<dR1p0) ++Nj4in_dR1p0;
            if(Nj8==2&&(ak4jet_pt[ii]>0)&&deltaRjet1ak4p[ii]>dR1p0&&deltaRjet2ak4p[ii]>dR1p0) ++Nj4ex_dR1p0;
            if(Nj8==2&&(ak4jet_pt[ii]>0)&&(deltaRjet1ak4p[ii]<dR1p0||deltaRjet2ak4p[ii]<dR1p0)) ++Nj4in_dR1p0;
            
            if(Nj4ex_dR1p0==1&&(ak4jet_pt[ii]>0)) {
                ak4expt1_dR1p0=ak4jet_pt[ii];ak4exeta1_dR1p0=ak4jet_eta[ii];ak4exphi1_dR1p0=ak4jet_phi[ii];ak4exe1_dR1p0=ak4jet_e[ii];
            }
            if(Nj4ex_dR1p0==2&&(ak4jet_pt[ii]>0)) {
                ak4expt2_dR1p0=ak4jet_pt[ii];ak4exeta2_dR1p0=ak4jet_eta[ii];ak4exphi2_dR1p0=ak4jet_phi[ii];ak4exe2_dR1p0=ak4jet_e[ii];
            }

            double dR1p2=1.2;
            
            if(Nj8==1&&(ak4jet_pt[ii]>0)&&deltaRjet1ak4p[ii]>dR1p2) ++Nj4ex_dR1p2;
            if(Nj8==1&&(ak4jet_pt[ii]>0)&&deltaRjet1ak4p[ii]<dR1p2) ++Nj4in_dR1p2;
            if(Nj8==2&&(ak4jet_pt[ii]>0)&&deltaRjet1ak4p[ii]<dR1p2&&deltaRjet2ak4p[ii]<dR1p2) ++Nj4ex_dR1p2;
            if(Nj8==2&&(ak4jet_pt[ii]>0)&&(deltaRjet1ak4p[ii]>dR1p2||deltaRjet2ak4p[ii]>dR1p2)) ++Nj4in_dR1p2;
            
            if(Nj4ex_dR1p2==1&&(ak4jet_pt[ii]>0)) {
                ak4expt1_dR1p2=ak4jet_pt[ii];ak4exeta1_dR1p2=ak4jet_eta[ii];ak4exphi1_dR1p2=ak4jet_phi[ii];ak4exe1_dR1p2=ak4jet_e[ii];
            }
            if(Nj4ex_dR1p2==2&&(ak4jet_pt[ii]>0)) {
                ak4expt2_dR1p2=ak4jet_pt[ii];ak4exeta2_dR1p2=ak4jet_eta[ii];ak4exphi2_dR1p2=ak4jet_phi[ii];ak4exe2_dR1p2=ak4jet_e[ii];
            }


        }
        //jet_mass_puppi_caseb = ak8jet1_2.Mag();
        }
        
        
        m_Jlvj=-99;m_JJlvj=-99;m_Jlvjj=-99;m_JJlvjj=-99;m_JJJlv=-99;m_JJJlvj=-99;m_JJJlvjj=-99;
        TLorentzVector ak8jet1_3,ak8jet2_3,ak8jet3_3,ak4exjet1,ak4exjet2,Wlv;
        ak8jet1_3.SetPtEtaPhiM(jet_pt_puppi,jetAK8puppi_eta,jetAK8puppi_phi,jet_mass_puppi_corr);
        ak8jet2_3.SetPtEtaPhiM(jet_pt_puppi_2,jetAK8puppi_eta_2,jetAK8puppi_phi_2,jet_mass_puppi_corr_2);
        ak8jet3_3.SetPtEtaPhiM(jet_pt_puppi_3,jetAK8puppi_eta_3,jetAK8puppi_phi_3,jet_mass_puppi_corr_3);
        Wlv.SetPtEtaPhiM(W_pt,W_eta,W_phi,massVlepJEC);
        m_JJJlv=(ak8jet1_3+ak8jet2_3+ak8jet3_3+Wlv).Mag();
        if(Nj4ex_dR0p8==0){
            m_Jlvj=m_jlv;m_Jlvjj=m_jlv;
            if(Nj8>=2)      {m_JJlvj=m_JJlv;m_JJlvjj=m_JJlv;}
            if(Nj8>=3)      {m_JJJlvj=m_JJJlv;m_JJJlvjj=m_JJJlv;}
        }
        if(Nj4ex_dR0p8==1){
            ak4exjet1.SetPtEtaPhiE(ak4expt1_dR0p8,ak4exeta1_dR0p8,ak4exphi1_dR0p8,ak4exe1_dR0p8);
            m_Jlvj = (ak8jet1_3+Wlv+ak4exjet1).Mag(); m_Jlvjj=m_Jlvj;
            if(Nj8==2)  {m_JJlvj=(ak8jet1_3+ak8jet2_3+Wlv+ak4exjet1).Mag();m_JJlvjj=m_JJlvj;}
            if(Nj8>=3)  {m_JJJlvj=(ak8jet1_3+ak8jet2_3+ak8jet3_3+Wlv+ak4exjet1).Mag();m_JJJlvjj=m_JJlvj;}
        }
        if(Nj4ex_dR0p8>=2){
            ak4exjet1.SetPtEtaPhiE(ak4expt1_dR0p8,ak4exeta1_dR0p8,ak4exphi1_dR0p8,ak4exe1_dR0p8);
            ak4exjet2.SetPtEtaPhiE(ak4expt2_dR0p8,ak4exeta2_dR0p8,ak4exphi2_dR0p8,ak4exe2_dR0p8);
            m_Jlvj = (ak8jet1_3+Wlv+ak4exjet1).Mag(); m_Jlvjj=(ak8jet1_3+Wlv+ak4exjet1+ak4exjet2).Mag();
            if(Nj8==2)  {m_JJlvj=(ak8jet1_3+ak8jet2_3+Wlv+ak4exjet1).Mag();m_JJlvjj=(ak8jet1_3+ak8jet2_3+Wlv+ak4exjet1+ak4exjet2).Mag();}
            if(Nj8>=3)  {m_JJJlvj=(ak8jet1_3+ak8jet2_3+ak8jet3_3+Wlv+ak4exjet1).Mag();m_JJJlvjj=(ak8jet1_3+ak8jet2_3+ak8jet3_3+Wlv+ak4exjet1+ak4exjet2).Mag();}
        }

        dR_min_WW=-99;dPhi_min_WW=-99;dR_max_WW=-99;dPhi_max_WW=-99;dEta_max_WW=-99;dEta_max_WW=-99;
        dR_j1W=-99;dR_j2W=-99;dR_j3W=-99;dR_j1j2=-99;dR_j1j3=-99;dR_j2j3=-99;
        dPhi_j1W=-99;dPhi_j2W=-99;dPhi_j3W=-99;dPhi_j1j2=-99;dPhi_j1j3=-99;dPhi_j2j3=-99;
        dEta_j1W=-99;dEta_j2W=-99;dEta_j3W=-99;dEta_j1j2=-99;dEta_j1j3=-99;dEta_j2j3=-99;
        dR_min_WW_mass=-99;dPhi_min_WW_mass=-99;dR_max_WW_mass=-99;dPhi_max_WW_mass=-99;dEta_max_WW_mass=-99;dEta_max_WW_mass=-99;
        dR_jmaxW=-99;dR_jmidW=-99;dR_jminW=-99;dR_jmaxjmid=-99;dR_jmaxjmin=-99;dR_jmidjmin=-99;
        dPhi_jmaxW=-99;dPhi_jmidW=-99;dPhi_jminW=-99;dPhi_jmaxjmid=-99;dPhi_jmaxjmin=-99;dPhi_jmidjmin=-99;
        dEta_jmaxW=-99;dEta_jmidW=-99;dEta_jminW=-99;dEta_jmaxjmid=-99;dEta_jmaxjmin=-99;dEta_jmidjmin=-99;

        double dRWW[6],dPhiWW[6],dEtaWW[6];
        dRWW[0]=deltaR(W_eta,W_phi,jetAK8puppi_eta,jetAK8puppi_phi);
        dRWW[1]=deltaR(W_eta,W_phi,jetAK8puppi_eta_2,jetAK8puppi_phi_2);
        dRWW[3]=deltaR(W_eta,W_phi,jetAK8puppi_eta_3,jetAK8puppi_phi_3);
        dRWW[2]=deltaR(jetAK8puppi_eta,jetAK8puppi_phi,jetAK8puppi_eta_2,jetAK8puppi_phi_2);
        dRWW[4]=deltaR(jetAK8puppi_eta,jetAK8puppi_phi,jetAK8puppi_eta_3,jetAK8puppi_phi_3);
        dRWW[5]=deltaR(jetAK8puppi_eta_2,jetAK8puppi_phi_2,jetAK8puppi_eta_3,jetAK8puppi_phi_3);
        
        dPhiWW[0]=deltaPhi(W_phi,jetAK8puppi_phi);
        dPhiWW[1]=deltaPhi(W_phi,jetAK8puppi_phi_2);
        dPhiWW[3]=deltaPhi(W_phi,jetAK8puppi_phi_3);
        dPhiWW[2]=deltaPhi(jetAK8puppi_phi,jetAK8puppi_phi_2);
        dPhiWW[4]=deltaPhi(jetAK8puppi_phi,jetAK8puppi_phi_3);
        dPhiWW[5]=deltaPhi(jetAK8puppi_phi_2,jetAK8puppi_phi_3);

        dEtaWW[0]=(W_eta-jetAK8puppi_eta);
        dEtaWW[1]=(W_eta-jetAK8puppi_eta_2);
        dEtaWW[3]=(W_eta-jetAK8puppi_eta_3);
        dEtaWW[2]=(jetAK8puppi_eta-jetAK8puppi_eta_2);
        dEtaWW[4]=(jetAK8puppi_eta-jetAK8puppi_eta_3);
        dEtaWW[5]=(jetAK8puppi_eta_2-jetAK8puppi_eta_3);
        dR_j1W=dRWW[0];dR_j2W=dRWW[1];dR_j3W=dRWW[3];dR_j1j2=dRWW[2];dR_j1j3=dRWW[4];dR_j2j3=dRWW[5];
        dPhi_j1W=dPhiWW[0];dPhi_j2W=dPhiWW[1];dPhi_j3W=dPhiWW[3];dPhi_j1j2=dPhiWW[2];dPhi_j1j3=dPhiWW[4];dPhi_j2j3=dPhiWW[5];
        dEta_j1W=dEtaWW[0];dEta_j2W=dEtaWW[1];dEta_j3W=dEtaWW[3];dEta_j1j2=dEtaWW[2];dEta_j1j3=dEtaWW[4];dEta_j2j3=dEtaWW[5];
        if(Nj8==1) {dR_min_WW=dRWW[0];dR_max_WW=dRWW[0];
        dPhi_min_WW=dPhiWW[0];dPhi_max_WW=dPhiWW[0];
        dEta_min_WW=dEtaWW[0];dEta_max_WW=dEtaWW[0];}
        if(Nj8==2) {
            Int_t *indexxwwr=new Int_t[6];
            TMath::Sort(3,dRWW,indexxwwr,1);
            dR_min_WW=dRWW[indexxwwr[2]];dR_max_WW=dRWW[indexxwwr[0]];
            Int_t *indexxwwPhi=new Int_t[6];
            TMath::Sort(3,dPhiWW,indexxwwPhi,1);
            dPhi_min_WW=dPhiWW[indexxwwPhi[2]];dPhi_max_WW=dPhiWW[indexxwwPhi[0]];
            Int_t *indexxwwEta=new Int_t[6];
            TMath::Sort(3,dEtaWW,indexxwwEta,1);
            dEta_min_WW=dEtaWW[indexxwwEta[2]];dEta_max_WW=dEtaWW[indexxwwEta[0]];
        }
    if(Nj8>=3) {
        Int_t *indexxwwr=new Int_t[6];
        TMath::Sort(6,dRWW,indexxwwr,1);
        dR_min_WW=dRWW[indexxwwr[5]];dR_max_WW=dRWW[indexxwwr[0]];
        Int_t *indexxwwPhi=new Int_t[6];
        TMath::Sort(6,dPhiWW,indexxwwPhi,1);
        dPhi_min_WW=dPhiWW[indexxwwPhi[5]];dPhi_max_WW=dPhiWW[indexxwwPhi[0]];
        Int_t *indexxwwEta=new Int_t[6];
        TMath::Sort(6,dEtaWW,indexxwwEta,1);
        dEta_min_WW=dEtaWW[indexxwwEta[5]];dEta_max_WW=dEtaWW[indexxwwEta[0]];
    }
        //cout<<"NJ  " <<Nj8<<dR_min_WW<<"    "<<dR_max_WW<<"    "<<dR_j1W<<"    "<<dR_j2W<<"    "<<dR_j3W<<"    "<<dR_j1j2<<"    "<<dR_j1j3<<"    "<<dR_j2j3<<"    "<<endl;
    double dRWW_mass[6],dPhiWW_mass[6],dEtaWW_mass[6];
    dRWW_mass[0]=deltaR(W_eta,W_phi,eta_max,phi_max);
    dRWW_mass[3]=deltaR(W_eta,W_phi,eta_mid,phi_mid);
    dRWW_mass[1]=deltaR(W_eta,W_phi,eta_min,phi_min);
    dRWW_mass[4]=deltaR(eta_max,phi_max,eta_mid,phi_mid);
    dRWW_mass[2]=deltaR(eta_max,phi_max,eta_min,phi_min);
    dRWW_mass[5]=deltaR(eta_mid,phi_mid,eta_min,phi_min);
    
    dPhiWW_mass[0]=deltaPhi(W_phi,phi_max);
    dPhiWW_mass[3]=deltaPhi(W_phi,phi_mid);
    dPhiWW_mass[1]=deltaPhi(W_phi,phi_min);
    dPhiWW_mass[4]=deltaPhi(phi_max,phi_mid);
    dPhiWW_mass[2]=deltaPhi(phi_max,phi_min);
    dPhiWW_mass[5]=deltaPhi(phi_mid,phi_min);
    
    dEtaWW_mass[0]=(W_eta-eta_max);
    dEtaWW_mass[3]=(W_eta-eta_mid);
    dEtaWW_mass[2]=(W_eta-eta_min);
    dEtaWW_mass[4]=(eta_max-eta_mid);
    dEtaWW_mass[2]=(eta_max-eta_min);
    dEtaWW_mass[5]=(eta_mid-eta_min);
    dR_jmaxW=dRWW_mass[0];dR_jmidW=dRWW_mass[3];dR_jminW=dRWW_mass[1];dR_jmaxjmid=dRWW_mass[4];dR_jmaxjmin=dRWW_mass[2];dR_jmidjmin=dRWW_mass[5];
    dPhi_jmaxW=dPhiWW_mass[0];dPhi_jmidW=dPhiWW_mass[3];dPhi_jminW=dPhiWW_mass[2];dPhi_jmaxjmid=dPhiWW_mass[4];dPhi_jmaxjmin=dPhiWW_mass[2];dPhi_jmidjmin=dPhiWW_mass[5];
    dEta_jmaxW=dEtaWW_mass[0];dEta_jmidW=dEtaWW_mass[3];dEta_jminW=dEtaWW_mass[1];dEta_jmaxjmid=dEtaWW_mass[4];dEta_jmaxjmin=dEtaWW_mass[2];dEta_jmidjmin=dEtaWW_mass[5];
    if(Nj8==1) {dR_min_WW_mass=dRWW_mass[0];dR_max_WW_mass=dRWW_mass[0];
        dPhi_min_WW_mass=dPhiWW_mass[0];dPhi_max_WW_mass=dPhiWW_mass[0];
        dEta_min_WW_mass=dEtaWW_mass[0];dEta_max_WW_mass=dEtaWW_mass[0];}
if(Nj8==2) {
    Int_t *indexxwwr_mass=new Int_t[6];
    TMath::Sort(3,dRWW_mass,indexxwwr_mass,1);
    dR_min_WW_mass=dRWW_mass[indexxwwr_mass[2]];dR_max_WW_mass=dRWW_mass[indexxwwr_mass[0]];
    Int_t *indexxwwPhi_mass=new Int_t[6];
    TMath::Sort(3,dPhiWW_mass,indexxwwPhi_mass,1);
    dPhi_min_WW_mass=dPhiWW_mass[indexxwwPhi_mass[2]];dPhi_max_WW_mass=dPhiWW_mass[indexxwwPhi_mass[0]];
    Int_t *indexxwwEta_mass=new Int_t[6];
    TMath::Sort(3,dEtaWW_mass,indexxwwEta_mass,1);
    dEta_min_WW_mass=dEtaWW_mass[indexxwwEta_mass[2]];dEta_max_WW_mass=dEtaWW_mass[indexxwwEta_mass[0]];
}
if(Nj8>=3) {
    Int_t *indexxwwr_mass=new Int_t[6];
    TMath::Sort(6,dRWW_mass,indexxwwr_mass,1);
    dR_min_WW_mass=dRWW_mass[indexxwwr_mass[5]];dR_max_WW_mass=dRWW_mass[indexxwwr_mass[0]];
    Int_t *indexxwwPhi_mass=new Int_t[6];
    TMath::Sort(6,dPhiWW_mass,indexxwwPhi_mass,1);
    dPhi_min_WW_mass=dPhiWW_mass[indexxwwPhi_mass[5]];dPhi_max_WW_mass=dPhiWW_mass[indexxwwPhi_mass[0]];
    Int_t *indexxwwEta_mass=new Int_t[6];
    TMath::Sort(6,dEtaWW_mass,indexxwwEta_mass,1);
    dEta_min_WW_mass=dEtaWW_mass[indexxwwEta_mass[5]];dEta_max_WW_mass=dEtaWW_mass[indexxwwEta_mass[0]];
}
        for(int i=0;i<3;i++){
            matchedt[i]=-99;matchedw[i]=-99;unmatched[i]=-99;matchedz[i]=-99;matchedg[i]=-99;matchedr[i]=-99;
            matchedft[i]=-99;matchedfw[i]=-99;unmatchedf[i]=-99;
            matched0p4t[i]=-99;matched0p4w[i]=-99;unmatched0p4[i]=-99;matched0p4z[i]=-99;matched0p4g[i]=-99;matched0p4r[i]=-99;
            matched0p4ft[i]=-99;matched0p4fw[i]=-99;unmatched0p4f[i]=-99;
            matched0p8t[i]=-99;matched0p8w[i]=-99;unmatched0p8[i]=-99;matched0p8z[i]=-99;matched0p8g[i]=-99;matched0p8r[i]=-99;
            matched0p8ft[i]=-99;matched0p8fw[i]=-99;unmatched0p8f[i]=-99;}
        double deltaRtj[3]={99,99,99},deltaRantitj[3]={99,99,99},deltaRwj[5][3],deltaRwfj[5][3],deltaRzj[5][3],deltaRgj[10][3],deltaRzfj[5][3],deltaRgfj[10][3],deltaRrj[3]={99,99,99};
        if (gentop_pt>0){
            deltaRtj[0] = deltaR(gentop_eta,gentop_phi,eta_max,phi_max);
            deltaRtj[1] = deltaR(gentop_eta,gentop_phi,eta_mid,phi_mid);
            deltaRtj[2] = deltaR(gentop_eta,gentop_phi,eta_min,phi_min);}
        if(IsData==10&&gen_rad_pt>0){
            deltaRrj[0] = deltaR(gen_rad_eta,gen_rad_phi,eta_max,phi_max);
            deltaRrj[1] = deltaR(gen_rad_eta,gen_rad_phi,eta_mid,phi_mid);
            deltaRrj[2] = deltaR(gen_rad_eta,gen_rad_phi,eta_min,phi_min);}
        if (genantitop_pt>0){
            deltaRantitj[0] = deltaR(genantitop_eta,genantitop_phi,eta_max,phi_max);
            deltaRantitj[1] = deltaR(genantitop_eta,genantitop_phi,eta_mid,phi_mid);
            deltaRantitj[2] = deltaR(genantitop_eta,genantitop_phi,eta_min,phi_min);}
        for(int i=0;i<5;i++){
            for (int j=0;j<3;j++){
                deltaRwj[i][j]=99;
                deltaRwfj[i][j]=99;
                deltaRzj[i][j]=99;
                deltaRzfj[i][j]=99;
                }
            if(ptgenwl[i]>0){
                deltaRwj[i][0] = deltaR(etagenwl[i],phigenwl[i],eta_max,phi_max);
                deltaRwj[i][1] = deltaR(etagenwl[i],phigenwl[i],eta_mid,phi_mid);
                deltaRwj[i][2] = deltaR(etagenwl[i],phigenwl[i],eta_min,phi_min);}
            if(ptgenwf[i]>0){
                deltaRwfj[i][0] = deltaR(etagenwf[i],phigenwf[i],eta_max,phi_max);
                deltaRwfj[i][1] = deltaR(etagenwf[i],phigenwf[i],eta_mid,phi_mid);
                deltaRwfj[i][2] = deltaR(etagenwf[i],phigenwf[i],eta_min,phi_min);}
            if(ptgenzl[i]>0){
                deltaRzj[i][0] = deltaR(etagenzl[i],phigenzl[i],eta_max,phi_max);
                deltaRzj[i][1] = deltaR(etagenzl[i],phigenzl[i],eta_mid,phi_mid);
                deltaRzj[i][2] = deltaR(etagenzl[i],phigenzl[i],eta_min,phi_min);}
        }
        for(int i=0;i<10;i++){
            for (int j=0;j<3;j++){
                deltaRgj[i][j]=99;
                deltaRgfj[i][j]=99;}
            if(ptgengl[i]>0){
                deltaRgj[i][0] = deltaR(etagengl[i],phigengl[i],eta_max,phi_max);
                deltaRgj[i][1] = deltaR(etagengl[i],phigengl[i],eta_mid,phi_mid);
                deltaRgj[i][2] = deltaR(etagengl[i],phigengl[i],eta_min,phi_min);}
        }
        
        double deltaRwj1[5],deltaRwj2[5],deltaRwj3[5];
        for(int i=0;i<5;i++){
            deltaRwj1[i]=deltaRwj[i][0];
            deltaRwj2[i]=deltaRwj[i][1];
            deltaRwj3[i]=deltaRwj[i][2];
            if(Nj8==2) deltaRwj2[i]=deltaRwj[i][2];
        }
        Int_t *mindrwj1=new Int_t[5];
        TMath::Sort(5,deltaRwj1,mindrwj1,0);
        Int_t *mindrwj2=new Int_t[5];
        TMath::Sort(5,deltaRwj2,mindrwj2,0);
        Int_t *mindrwj3=new Int_t[5];
        TMath::Sort(5,deltaRwj3,mindrwj3,0);
        
        double rmatched=0.6;
        double masstorw=100;//110->100
        if(Nj8==1){
        if (deltaRtj[0]<rmatched||deltaRantitj[0]<rmatched) {matchedt[0]=1;matchedft[0]=1;}
        if (deltaRrj[0]<rmatched) {matchedr[0]=1;}
        for (int i=0;i<5;i++){
            if (deltaRwj[i][0]<rmatched) matchedw[0]=1;
            if (deltaRwfj[i][0]<rmatched) matchedfw[0]=1;
            if (deltaRzj[i][0]<rmatched) matchedz[0]=1;
        }
        for (int i=0;i<10;i++){if (deltaRgj[i][0]<rmatched) matchedg[0]=1;}

        if(matchedt[0]<=0&&matchedr[0]<=0&&matchedw[0]<=0 && matchedz[0]<=0&&matchedg[0]<=0)  unmatched[0]=1;
            matchedtorg[0]=matchedt[0];matchedrorg[0]=matchedr[0];
            if(IsData>0 && IsData!=10&&jet_mass_puppi_corr<masstorw && min(deltaRtj[0],deltaRantitj[0])> deltaRwj1[mindrwj1[0]]&&matchedt[0]==1) matchedt[0]=0;
            if(IsData==10&&jet_mass_puppi_corr<masstorw && deltaRrj[0]> deltaRwj1[mindrwj1[0]]&&matchedr[0]==1) matchedr[0]=0;
        if(matchedft[0]<=0&&matchedfw[0]<=0)  unmatchedf[0]=1;
        }
        if(Nj8==2){
            if (deltaRtj[0]<rmatched||deltaRantitj[0]<rmatched) {matchedt[0]=1;matchedft[0]=1;}
            if (deltaRrj[0]<rmatched) {matchedr[0]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][0]<rmatched) matchedw[0]=1;
                if (deltaRwfj[i][0]<rmatched) matchedfw[0]=1;
                if (deltaRzj[i][0]<rmatched) matchedz[0]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[i][0]<rmatched) matchedg[0]=1;}
            if(matchedt[0]<=0&&matchedr[0]<=0&&matchedw[0]<=0 && matchedz[0]<=0&&matchedg[0]<=0)  unmatched[0]=1;
            matchedtorg[0]=matchedt[0];matchedrorg[0]=matchedr[0];
            if(IsData>0 && IsData!=10&&Mj_max<masstorw && min(deltaRtj[0],deltaRantitj[0])> deltaRwj1[mindrwj1[0]]&&matchedt[0]==1) matchedt[0]=0;
            if(IsData==10&&Mj_max<masstorw && deltaRrj[0]> deltaRwj1[mindrwj1[0]]&&matchedr[0]==1) matchedr[0]=0;

            if(matchedft[0]<=0&&matchedfw[0]<=0)  unmatchedf[0]=1;
            
            if (deltaRtj[2]<rmatched||deltaRantitj[2]<rmatched) {matchedt[1]=1;matchedft[1]=1;}
            if (deltaRrj[2]<rmatched) {matchedr[1]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][2]<rmatched) matchedw[1]=1;
                if (deltaRwfj[i][2]<rmatched) matchedfw[1]=1;
                if (deltaRzj[i][2]<rmatched) matchedz[1]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[i][2]<rmatched) matchedg[2]=1;}

            if(matchedt[1]<=0&&matchedr[1]<=0&&matchedw[1]<=0 && matchedz[1]<=0&&matchedg[1]<=0)  unmatched[1]=1;
            matchedtorg[1]=matchedt[1];matchedrorg[1]=matchedr[1];
            if(IsData>0 && IsData!=10&&Mj_mid<masstorw && min(deltaRtj[2],deltaRantitj[2])> deltaRwj2[mindrwj2[0]]&&matchedt[1]==1) matchedt[1]=0;
            if(IsData==10&&Mj_mid<masstorw && deltaRrj[2]> deltaRwj2[mindrwj2[0]]&&matchedr[1]==1) matchedr[1]=0;

            if(matchedft[1]<=0&&matchedfw[1]<=0)  unmatchedf[1]=1;
        }
        if(Nj8==3){
            if (deltaRtj[0]<rmatched||deltaRantitj[0]<rmatched) {matchedt[0]=1;matchedft[0]=1;}
            if (deltaRrj[0]<rmatched) {matchedr[0]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][0]<rmatched) matchedw[0]=1;
                if (deltaRwfj[i][0]<rmatched) matchedfw[0]=1;
                if (deltaRzj[i][0]<rmatched) matchedz[0]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[i][0]<rmatched) matchedg[0]=1;}
            if(matchedt[0]<=0&&matchedr[0]<=0&&matchedw[0]<=0 && matchedz[0]<=0&&matchedg[0]<=0)  unmatched[0]=1;
            matchedtorg[0]=matchedt[0];matchedrorg[0]=matchedr[0];
            if(IsData>0 && IsData!=10&&Mj_max<masstorw && min(deltaRtj[0],deltaRantitj[0])> deltaRwj1[mindrwj1[0]]&&matchedt[0]==1) matchedt[0]=0;
            if(IsData==10&&Mj_max<masstorw && deltaRrj[0]> deltaRwj1[mindrwj1[0]]&&matchedr[0]==1) matchedr[0]=0;
            if(matchedft[0]<=0&&matchedfw[0]<=0)  unmatchedf[0]=1;
            
            if (deltaRtj[1]<rmatched||deltaRantitj[1]<rmatched) {matchedt[1]=1;matchedft[1]=1;}
            if (deltaRrj[1]<rmatched) {matchedr[1]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][1]<rmatched) matchedw[1]=1;
                if (deltaRwfj[i][1]<rmatched) matchedfw[1]=1;
                if (deltaRzj[i][1]<rmatched) matchedz[1]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[i][1]<rmatched) matchedg[1]=1;}
            if(matchedt[1]<=0&&matchedr[1]<=0&&matchedw[1]<=0 && matchedz[1]<=0&&matchedg[1]<=0)  unmatched[1]=1;
            matchedtorg[1]=matchedt[1];matchedrorg[1]=matchedr[1];
            if(IsData>0 && IsData!=10&&Mj_mid<masstorw && min(deltaRtj[1],deltaRantitj[1])> deltaRwj2[mindrwj2[0]]&&matchedt[1]==1) matchedt[1]=0;
            if(IsData==10&&Mj_mid<masstorw && deltaRrj[1]> deltaRwj2[mindrwj2[0]]&&matchedr[1]==1) matchedr[1]=0;
            if(matchedft[1]<=0&&matchedfw[1]<=0)  unmatchedf[1]=1;
            
            if (deltaRtj[2]<rmatched||deltaRantitj[2]<rmatched) {matchedt[2]=1;matchedft[2]=1;}
            if (deltaRrj[2]<rmatched) {matchedr[2]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][2]<rmatched) matchedw[2]=1;
                if (deltaRwfj[i][2]<rmatched) matchedfw[2]=1;
                if (deltaRzj[i][2]<rmatched) matchedz[2]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[2][0]<rmatched) matchedg[2]=1;}
            if(matchedt[2]<=0&&matchedr[2]<=0&&matchedw[2]<=0 && matchedz[2]<=0&&matchedg[2]<=0)  unmatched[2]=1;
            matchedtorg[2]=matchedt[2];matchedrorg[2]=matchedr[2];
            if(IsData>0 && IsData!=10&&Mj_min<masstorw && min(deltaRtj[2],deltaRantitj[2])> deltaRwj3[mindrwj3[0]]&&matchedt[2]==1) matchedt[2]=0;
            if(IsData==10&&Mj_min<masstorw && deltaRrj[2]> deltaRwj3[mindrwj3[0]]&&matchedr[2]==1) matchedr[2]=0;
            if(matchedft[2]<=0&&matchedfw[2]<=0)  unmatchedf[2]=1;
        }
        
        double rmatched0p4=0.4;
        if(Nj8==1){
            if (deltaRtj[0]<rmatched0p4||deltaRantitj[0]<rmatched0p4) {matched0p4t[0]=1;matched0p4ft[0]=1;}
            if (deltaRrj[0]<rmatched0p4) {matched0p4r[0]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][0]<rmatched0p4) matched0p4w[0]=1;
                if (deltaRwfj[i][0]<rmatched0p4) matched0p4fw[0]=1;
                if (deltaRzj[i][0]<rmatched0p4) matched0p4z[0]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[i][0]<rmatched0p4) matched0p4g[0]=1;}

            if(matched0p4t[0]<=0&&matched0p4r[0]<=0&&matched0p4w[0]<=0 && matched0p4z[0]<=0&&matched0p4g[0]<=0)  unmatched0p4[0]=1;
            matched0p4torg[0]=matched0p4t[0];matched0p4rorg[0]=matched0p4r[0];
            if(IsData>0 && IsData!=10&&jet_mass_puppi_corr<masstorw && min(deltaRtj[0],deltaRantitj[0])> deltaRwj1[mindrwj1[0]]&&matched0p4t[0]==1) matched0p4t[0]=0;
            if(IsData==10&&jet_mass_puppi_corr<masstorw && deltaRrj[0]> deltaRwj1[mindrwj1[0]]&&matched0p4r[0]==1) matched0p4r[0]=0;
            if(matched0p4ft[0]<=0&&matched0p4fw[0]<=0)  unmatched0p4f[0]=1;
        }
        if(Nj8==2){
            if (deltaRtj[0]<rmatched0p4||deltaRantitj[0]<rmatched0p4) {matched0p4t[0]=1;matched0p4ft[0]=1;}
            if (deltaRrj[0]<rmatched0p4) {matched0p4r[0]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][0]<rmatched0p4) matched0p4w[0]=1;
                if (deltaRwfj[i][0]<rmatched0p4) matched0p4fw[0]=1;
                if (deltaRzj[i][0]<rmatched0p4) matched0p4z[0]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[i][0]<rmatched0p4) matched0p4g[0]=1;}

            if(matched0p4t[0]<=0&&matched0p4r[0]<=0&&matched0p4w[0]<=0 && matched0p4z[0]<=0&&matched0p4g[0]<=0)  unmatched0p4[0]=1;
            matched0p4torg[0]=matched0p4t[0];matched0p4rorg[0]=matched0p4r[0];
            if(IsData>0 && IsData!=10&&Mj_max<masstorw && min(deltaRtj[0],deltaRantitj[0])> deltaRwj1[mindrwj1[0]]&&matched0p4t[0]==1) matched0p4t[0]=0;
            if(IsData==10&&Mj_max<masstorw && deltaRrj[0]> deltaRwj1[mindrwj1[0]]&&matched0p4r[0]==1) matched0p4r[0]=0;
            
            if(matched0p4ft[0]<=0&&matched0p4fw[0]<=0)  unmatched0p4f[0]=1;
            
            if (deltaRtj[2]<rmatched0p4||deltaRantitj[2]<rmatched0p4) {matched0p4t[1]=1;matched0p4ft[1]=1;}
            if (deltaRrj[2]<rmatched0p4) {matched0p4r[1]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][2]<rmatched0p4) matched0p4w[1]=1;
                if (deltaRwfj[i][2]<rmatched0p4) matched0p4fw[1]=1;
                if (deltaRzj[i][2]<rmatched0p4) matched0p4z[1]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[i][2]<rmatched0p4) matched0p4g[1]=1;}

            if(matched0p4t[1]<=0&&matched0p4r[1]<=0&&matched0p4w[1]<=0 && matched0p4z[1]<=0&&matched0p4g[1]<=0)  unmatched0p4[1]=1;
            matched0p4torg[1]=matched0p4t[1];matched0p4rorg[1]=matched0p4r[1];
            if(IsData>0 && IsData!=10&&Mj_mid<masstorw && min(deltaRtj[2],deltaRantitj[2])> deltaRwj2[mindrwj2[0]]&&matched0p4t[1]==1) matched0p4t[1]=0;
            if(IsData==10&&Mj_mid<masstorw && deltaRrj[2]> deltaRwj2[mindrwj2[0]]&&matched0p4r[1]==1) matched0p4r[1]=0;
            
            if(matched0p4ft[1]<=0&&matched0p4fw[1]<=0)  unmatched0p4f[1]=1;
        }
        if(Nj8==3){
            if (deltaRtj[0]<rmatched0p4||deltaRantitj[0]<rmatched0p4) {matched0p4t[0]=1;matched0p4ft[0]=1;}
            if (deltaRrj[0]<rmatched0p4) {matched0p4r[0]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][0]<rmatched0p4) matched0p4w[0]=1;
                if (deltaRwfj[i][0]<rmatched0p4) matched0p4fw[0]=1;
                if (deltaRzj[i][0]<rmatched0p4) matched0p4z[0]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[i][0]<rmatched0p4) matched0p4g[0]=1;}

            if(matched0p4t[0]<=0&&matched0p4r[0]<=0&&matched0p4w[0]<=0 && matched0p4z[0]<=0&&matched0p4g[0]<=0)  unmatched0p4[0]=1;
            matched0p4torg[0]=matched0p4t[0];matched0p4rorg[0]=matched0p4r[0];
            if(IsData>0 && IsData!=10&&Mj_max<masstorw && min(deltaRtj[0],deltaRantitj[0])> deltaRwj1[mindrwj1[0]]&&matched0p4t[0]==1) matched0p4t[0]=0;
            if(IsData==10&&Mj_max<masstorw && deltaRrj[0]> deltaRwj1[mindrwj1[0]]&&matched0p4r[0]==1) matched0p4r[0]=0;
            if(matched0p4ft[0]<=0&&matched0p4fw[0]<=0)  unmatched0p4f[0]=1;
            
            if (deltaRtj[1]<rmatched0p4||deltaRantitj[1]<rmatched0p4) {matched0p4t[1]=1;matched0p4ft[1]=1;}
            if (deltaRrj[1]<rmatched0p4) {matched0p4r[1]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][1]<rmatched0p4) matched0p4w[1]=1;
                if (deltaRwfj[i][1]<rmatched0p4) matched0p4fw[1]=1;
                if (deltaRzj[i][1]<rmatched0p4) matched0p4z[1]=1;
                if (deltaRgj[i][1]<rmatched0p4) matched0p4g[1]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[i][1]<rmatched0p4) matched0p4g[1]=1;}

            if(matched0p4t[1]<=0&&matched0p4r[1]<=0&&matched0p4w[1]<=0 && matched0p4z[1]<=0&&matched0p4g[1]<=0)  unmatched0p4[1]=1;
            matched0p4torg[1]=matched0p4t[1];matched0p4rorg[1]=matched0p4r[1];
            if(IsData>0 && IsData!=10&&Mj_mid<masstorw && min(deltaRtj[1],deltaRantitj[1])> deltaRwj2[mindrwj2[0]]&&matched0p4t[1]==1) matched0p4t[1]=0;
            if(IsData==10&&Mj_mid<masstorw && deltaRrj[1]> deltaRwj2[mindrwj2[0]]&&matched0p4r[1]==1) matched0p4r[1]=0;
            if(matched0p4ft[1]<=0&&matched0p4fw[1]<=0)  unmatched0p4f[1]=1;
            
            if (deltaRtj[2]<rmatched0p4||deltaRantitj[2]<rmatched0p4) {matched0p4t[2]=1;matched0p4ft[2]=1;}
            if (deltaRrj[2]<rmatched0p4) {matched0p4r[2]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][2]<rmatched0p4) matched0p4w[2]=1;
                if (deltaRwfj[i][2]<rmatched0p4) matched0p4fw[2]=1;
                if (deltaRzj[i][2]<rmatched0p4) matched0p4z[2]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[i][2]<rmatched0p4) matched0p4g[2]=1;}

            if(matched0p4t[2]<=0&&matched0p4r[2]<=0&&matched0p4w[2]<=0 && matched0p4z[2]<=0&&matched0p4g[2]<=0)  unmatched0p4[2]=1;
            matched0p4torg[2]=matched0p4t[2];matched0p4rorg[2]=matched0p4r[2];
            if(IsData>0 && IsData!=10&&Mj_min<masstorw && min(deltaRtj[2],deltaRantitj[2])> deltaRwj3[mindrwj3[0]]&&matched0p4t[2]==1) matched0p4t[2]=0;
            if(IsData==10&&Mj_min<masstorw && deltaRrj[2]> deltaRwj3[mindrwj3[0]]&&matched0p4r[2]==1) matched0p4r[2]=0;
            if(matched0p4ft[2]<=0&&matched0p4fw[2]<=0)  unmatched0p4f[2]=1;
        }

        double rmatched0p8=0.8;
        if(Nj8==1){
            if (deltaRtj[0]<rmatched0p8||deltaRantitj[0]<rmatched0p8) {matched0p8t[0]=1;matched0p8ft[0]=1;}
            if (deltaRrj[0]<rmatched0p8) {matched0p8r[0]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][0]<rmatched0p8) matched0p8w[0]=1;
                if (deltaRwfj[i][0]<rmatched0p8) matched0p8fw[0]=1;
                if (deltaRzj[i][0]<rmatched0p8) matched0p8z[0]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[i][0]<rmatched0p8) matched0p8g[0]=1;}

            if(matched0p8t[0]<=0&&matched0p8r[0]<=0&&matched0p8w[0]<=0 && matched0p8z[0]<=0&&matched0p8g[0]<=0)  unmatched0p8[0]=1;
            matched0p8torg[0]=matched0p8t[0];matched0p8rorg[0]=matched0p8r[0];
            if(IsData>0 && IsData!=10&&jet_mass_puppi_corr<masstorw && min(deltaRtj[0],deltaRantitj[0])> deltaRwj1[mindrwj1[0]]&&matched0p8t[0]==1) matched0p8t[0]=0;
            if(IsData==10&&jet_mass_puppi_corr<masstorw && deltaRrj[0]> deltaRwj1[mindrwj1[0]]&&matched0p8r[0]==1) matched0p8r[0]=0;
            if(matched0p8ft[0]<=0&&matched0p8fw[0]<=0)  unmatched0p8f[0]=1;
        }
        if(Nj8==2){
            if (deltaRtj[0]<rmatched0p8||deltaRantitj[0]<rmatched0p8) {matched0p8t[0]=1;matched0p8ft[0]=1;}
            if (deltaRrj[0]<rmatched0p8) {matched0p8r[0]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][0]<rmatched0p8) matched0p8w[0]=1;
                if (deltaRwfj[i][0]<rmatched0p8) matched0p8fw[0]=1;
                if (deltaRzj[i][0]<rmatched0p8) matched0p8z[0]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[i][0]<rmatched0p8) matched0p8g[0]=1;}

            if(matched0p8t[0]<=0&&matched0p8r[0]<=0&&matched0p8w[0]<=0 && matched0p8z[0]<=0&&matched0p8g[0]<=0)  unmatched0p8[0]=1;
            matched0p8torg[0]=matched0p8t[0];matched0p8rorg[0]=matched0p8r[0];
            if(IsData>0 && IsData!=10&&Mj_max<masstorw && min(deltaRtj[0],deltaRantitj[0])> deltaRwj1[mindrwj1[0]]&&matched0p8t[0]==1) matched0p8t[0]=0;
            if(IsData==10&&Mj_max<masstorw && deltaRrj[0]> deltaRwj1[mindrwj1[0]]&&matched0p8r[0]==1) matched0p8r[0]=0;
            
            if(matched0p8ft[0]<=0&&matched0p8fw[0]<=0)  unmatched0p8f[0]=1;
            
            if (deltaRtj[2]<rmatched0p8||deltaRantitj[2]<rmatched0p8) {matched0p8t[1]=1;matched0p8ft[1]=1;}
            if (deltaRrj[2]<rmatched0p8) {matched0p8r[1]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][2]<rmatched0p8) matched0p8w[1]=1;
                if (deltaRwfj[i][2]<rmatched0p8) matched0p8fw[1]=1;
                if (deltaRzj[i][2]<rmatched0p8) matched0p8z[1]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[i][2]<rmatched0p8) matched0p8g[1]=1;}

            if(matched0p8t[1]<=0&&matched0p8r[1]<=0&&matched0p8w[1]<=0 && matched0p8z[1]<=0&&matched0p8g[1]<=0)  unmatched0p8[1]=1;
            matched0p8torg[1]=matched0p8t[1];matched0p8rorg[1]=matched0p8r[1];
            if(IsData>0 && IsData!=10&&Mj_mid<masstorw && min(deltaRtj[2],deltaRantitj[2])> deltaRwj2[mindrwj2[0]]&&matched0p8t[1]==1) matched0p8t[1]=0;
            if(IsData==10&&Mj_mid<masstorw && deltaRrj[2]> deltaRwj2[mindrwj2[0]]&&matched0p8r[1]==1) matched0p8r[1]=0;
            
            if(matched0p8ft[1]<=0&&matched0p8fw[1]<=0)  unmatched0p8f[1]=1;
        }
        if(Nj8==3){
            if (deltaRtj[0]<rmatched0p8||deltaRantitj[0]<rmatched0p8) {matched0p8t[0]=1;matched0p8ft[0]=1;}
            if (deltaRrj[0]<rmatched0p8) {matched0p8r[0]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][0]<rmatched0p8) matched0p8w[0]=1;
                if (deltaRwfj[i][0]<rmatched0p8) matched0p8fw[0]=1;
                if (deltaRzj[i][0]<rmatched0p8) matched0p8z[0]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[i][0]<rmatched0p8) matched0p8g[0]=1;}

            if(matched0p8t[0]<=0&&matched0p8r[0]<=0&&matched0p8w[0]<=0 && matched0p8z[0]<=0&&matched0p8g[0]<=0)  unmatched0p8[0]=1;
            matched0p8torg[0]=matched0p8t[0];matched0p8rorg[0]=matched0p8r[0];
            if(IsData>0 && IsData!=10&&Mj_max<masstorw && min(deltaRtj[0],deltaRantitj[0])> deltaRwj1[mindrwj1[0]]&&matched0p8t[0]==1) matched0p8t[0]=0;
            if(IsData==10&&Mj_max<masstorw && deltaRrj[0]> deltaRwj1[mindrwj1[0]]&&matched0p8r[0]==1) matched0p8r[0]=0;
            if(matched0p8ft[0]<=0&&matched0p8fw[0]<=0)  unmatched0p8f[0]=1;
            
            if (deltaRtj[1]<rmatched0p8||deltaRantitj[1]<rmatched0p8) {matched0p8t[1]=1;matched0p8ft[1]=1;}
            if (deltaRrj[1]<rmatched0p8) {matched0p8r[1]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][1]<rmatched0p8) matched0p8w[1]=1;
                if (deltaRwfj[i][1]<rmatched0p8) matched0p8fw[1]=1;
                if (deltaRzj[i][1]<rmatched0p8) matched0p8z[1]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[i][1]<rmatched0p8) matched0p8g[1]=1;}

            if(matched0p8t[1]<=0&&matched0p8r[1]<=0&&matched0p8w[1]<=0 && matched0p8z[1]<=0&&matched0p8g[1]<=0)  unmatched0p8[1]=1;
            matched0p8torg[1]=matched0p8t[1];matched0p8rorg[1]=matched0p8r[1];
            if(IsData>0 && IsData!=10&&Mj_mid<masstorw && min(deltaRtj[1],deltaRantitj[1])> deltaRwj2[mindrwj2[0]]&&matched0p8t[1]==1) matched0p8t[1]=0;
            if(IsData==10&&Mj_mid<masstorw && deltaRrj[1]> deltaRwj2[mindrwj2[0]]&&matched0p8r[1]==1) matched0p8r[1]=0;
            if(matched0p8ft[1]<=0&&matched0p8fw[1]<=0)  unmatched0p8f[1]=1;
            
            if (deltaRtj[2]<rmatched0p8||deltaRantitj[2]<rmatched0p8) {matched0p8t[2]=1;matched0p8ft[2]=1;}
            if (deltaRrj[2]<rmatched0p8) {matched0p8r[2]=1;}
            for (int i=0;i<5;i++){
                if (deltaRwj[i][2]<rmatched0p8) matched0p8w[2]=1;
                if (deltaRwfj[i][2]<rmatched0p8) matched0p8fw[2]=1;
                if (deltaRzj[i][2]<rmatched0p8) matched0p8z[2]=1;
            }
            for (int i=0;i<10;i++){if (deltaRgj[i][2]<rmatched0p8) matched0p8g[2]=1;}

            if(matched0p8t[2]<=0&&matched0p8r[2]<=0&&matched0p8w[2]<=0 && matched0p8z[2]<=0&&matched0p8g[2]<=0)  unmatched0p8[2]=1;
            matched0p8torg[2]=matched0p8t[2];matched0p8rorg[2]=matched0p8r[2];
            if(IsData>0 && IsData!=10&&Mj_min<masstorw && min(deltaRtj[2],deltaRantitj[2])> deltaRwj3[mindrwj3[0]]&&matched0p8t[2]==1) matched0p8t[2]=0;
            if(IsData==10&&Mj_min<masstorw && deltaRrj[2]> deltaRwj3[mindrwj3[0]]&&matched0p8r[2]==1) matched0p8r[2]=0;
            if(matched0p8ft[2]<=0&&matched0p8fw[2]<=0)  unmatched0p8f[2]=1;
        }

        TLorentzVector wlep,bjet;
        wlep.SetPtEtaPhiM(W_pt,W_eta,W_phi,massVlepJEC);
        mtbtvlep=-99;        mtbmvlep=-99;deltaRbtvlep=-99;deltaRbmvlep=-99;deltaRblvlep=-99;mtblvlep=-99;mbmvlep=-99;
        deltaRantitop_wlepbm=99;deltaRtop_wlepbm=99;deltaRantitop_wlepbl=99;deltaRtop_wlepbl=99;mblvlep=-99;
        if(num_bJet_tight_deep==1){
        for(Int_t i=0; i<8; i++)  {
            deltaRAK4AK8_new[i]=0.;
            deltaRAK4AK8_new[i]=sqrt(pow(fabs(ak4jet_eta[i]-jetAK8puppi_eta),2)+pow(TMath::Min(fabs(ak4jet_phi[i]-jetAK8puppi_phi),2*Pi-fabs(ak4jet_phi[i]-jetAK8puppi_phi)),2));
            if(Nj8==2) deltaRAK4AK8_new_2[i]=sqrt(pow(fabs(ak4jet_eta[i]-jetAK8puppi_eta_2),2)+pow(TMath::Min(fabs(ak4jet_phi[i]-jetAK8puppi_phi_2),2*Pi-fabs(ak4jet_phi[i]-jetAK8puppi_phi_2)),2));
            if(ak4jet_pt[i]>30 && (ak4jet_deepcsvb[i]+ak4jet_deepcsvbb[i])>0.8958 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8_new[i]>=0.8 && (Nj8==1 || deltaRAK4AK8_new_2[i]>=0.8)) {bjet.SetPtEtaPhiE(ak4jet_pt[i],ak4jet_eta[i],ak4jet_phi[i],ak4jet_e[i]);
                mtbtvlep=Double_t(sqrt(massVlepJEC*massVlepJEC+bjet.Mag()*bjet.Mag()+2*ak4jet_pt[i]*W_pt*(1.0-cos(ak4jet_phi[i]-W_phi))));
                deltaRbtvlep=deltaR(W_eta,W_phi,ak4jet_eta[i],ak4jet_phi[i]);
                //cout<<"hhhh "<<massVlepJEC<<"   "<<bjet.Mag()<<"   "<<(wlep+bjet).Mt()<<"   "<<mtbtvlep<<"   "<<Double_t(sqrt(massVlepJEC*massVlepJEC+bjet.Mag()*bjet.Mag()+2*ak4jet_pt[i]*W_pt*(1.0-cos(ak4jet_phi[i]-W_phi))))<<"   "<<deltaRbtvlep<<endl;
            }
            
        }}
        if(num_bJet_deep>=1){
            double wlepak4[8]={99,99,99,99,99,99,99,99};
            int iclosewlepak4=0;
            for(Int_t i=0; i<8; i++)  {
                deltaRAK4AK8_new[i]=0.;
                deltaRAK4AK8_new[i]=sqrt(pow(fabs(ak4jet_eta[i]-jetAK8puppi_eta),2)+pow(TMath::Min(fabs(ak4jet_phi[i]-jetAK8puppi_phi),2*Pi-fabs(ak4jet_phi[i]-jetAK8puppi_phi)),2));
                if(Nj8==2) deltaRAK4AK8_new_2[i]=sqrt(pow(fabs(ak4jet_eta[i]-jetAK8puppi_eta_2),2)+pow(TMath::Min(fabs(ak4jet_phi[i]-jetAK8puppi_phi_2),2*Pi-fabs(ak4jet_phi[i]-jetAK8puppi_phi_2)),2));
                if(ak4jet_pt[i]>30 && (ak4jet_deepcsvb[i]+ak4jet_deepcsvbb[i])>0.6324 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8_new[i]>=0.8 && (Nj8==1 || deltaRAK4AK8_new_2[i]>=0.8)) {
                    wlepak4[i]=deltaR(W_eta,W_phi,ak4jet_eta[i],ak4jet_phi[i]);
                    //cout<<"hhhhm "<<massVlepJEC<<"   "<<bjet.Mag()<<"   "<<(wlep+bjet).Mt()<<"   "<<mtbmvlep<<"   "<<Double_t(sqrt(massVlepJEC*massVlepJEC+bjet.Mag()*bjet.Mag()+2*ak4jet_pt[i]*W_pt*(1.0-cos(ak4jet_phi[i]-W_phi))))<<"   "<<deltaRbmvlep<<endl;
                }
            }
            Int_t *indexx_wlepak4=new Int_t[8];
            TMath::Sort(8,wlepak4,indexx_wlepak4,0);
            bjet.SetPtEtaPhiE(ak4jet_pt[indexx_wlepak4[0]],ak4jet_eta[indexx_wlepak4[0]],ak4jet_phi[indexx_wlepak4[0]],ak4jet_e[indexx_wlepak4[0]]);
            mtbmvlep=Double_t(sqrt(massVlepJEC*massVlepJEC+bjet.Mag()*bjet.Mag()+2*ak4jet_pt[indexx_wlepak4[0]]*W_pt*(1.0-cos(ak4jet_phi[indexx_wlepak4[0]]-W_phi))));
            deltaRbmvlep=deltaR(W_eta,W_phi,ak4jet_eta[indexx_wlepak4[0]],ak4jet_phi[indexx_wlepak4[0]]);
            TLorentzVector wlepak4v;wlepak4v=wlep+bjet;
            mbmvlep=wlepak4v.Mag();
            deltaRtop_wlepbm=99;deltaRantitop_wlepbm=99;
            if (gentop_pt>0) deltaRtop_wlepbm = deltaR(gentop_eta,gentop_phi,wlepak4v.Eta(),wlepak4v.Phi());
            if (genantitop_pt >0) deltaRantitop_wlepbm = deltaR(genantitop_eta,genantitop_phi,wlepak4v.Eta(),wlepak4v.Phi());
        }
        if(num_bJet_loose_deep>=1){
            double wlepak4[8]={99,99,99,99,99,99,99,99};
            int iclosewlepak4=0;
            for(Int_t i=0; i<8; i++)  {
                deltaRAK4AK8_new[i]=0.;
                deltaRAK4AK8_new[i]=sqrt(pow(fabs(ak4jet_eta[i]-jetAK8puppi_eta),2)+pow(TMath::Min(fabs(ak4jet_phi[i]-jetAK8puppi_phi),2*Pi-fabs(ak4jet_phi[i]-jetAK8puppi_phi)),2));
                if(Nj8==2) deltaRAK4AK8_new_2[i]=sqrt(pow(fabs(ak4jet_eta[i]-jetAK8puppi_eta_2),2)+pow(TMath::Min(fabs(ak4jet_phi[i]-jetAK8puppi_phi_2),2*Pi-fabs(ak4jet_phi[i]-jetAK8puppi_phi_2)),2));
                if(ak4jet_pt[i]>30 && (ak4jet_deepcsvb[i]+ak4jet_deepcsvbb[i])>0.2219 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8_new[i]>=0.8 && (Nj8==1 || deltaRAK4AK8_new_2[i]>=0.8)) {
                    wlepak4[i]=deltaR(W_eta,W_phi,ak4jet_eta[i],ak4jet_phi[i]);
                    //cout<<"hhhhm "<<massVlepJEC<<"   "<<bjet.Mag()<<"   "<<(wlep+bjet).Mt()<<"   "<<mtblvlep<<"   "<<Double_t(sqrt(massVlepJEC*massVlepJEC+bjet.Mag()*bjet.Mag()+2*ak4jet_pt[i]*W_pt*(1.0-cos(ak4jet_phi[i]-W_phi))))<<"   "<<deltaRblvlep<<endl;
                }
            }
            Int_t *indexx_wlepak4=new Int_t[8];
            TMath::Sort(8,wlepak4,indexx_wlepak4,0);
            bjet.SetPtEtaPhiE(ak4jet_pt[indexx_wlepak4[0]],ak4jet_eta[indexx_wlepak4[0]],ak4jet_phi[indexx_wlepak4[0]],ak4jet_e[indexx_wlepak4[0]]);
            mtblvlep=Double_t(sqrt(massVlepJEC*massVlepJEC+bjet.Mag()*bjet.Mag()+2*ak4jet_pt[indexx_wlepak4[0]]*W_pt*(1.0-cos(ak4jet_phi[indexx_wlepak4[0]]-W_phi))));
            deltaRblvlep=deltaR(W_eta,W_phi,ak4jet_eta[indexx_wlepak4[0]],ak4jet_phi[indexx_wlepak4[0]]);
            TLorentzVector wlepak4v;wlepak4v=wlep+bjet;
            mblvlep=wlepak4v.Mag();
            deltaRtop_wlepbl=99;deltaRantitop_wlepbl=99;
            if (gentop_pt>0) deltaRtop_wlepbl = deltaR(gentop_eta,gentop_phi,wlepak4v.Eta(),wlepak4v.Phi());
            if (genantitop_pt >0) deltaRantitop_wlepbl = deltaR(genantitop_eta,genantitop_phi,wlepak4v.Eta(),wlepak4v.Phi());
        }
        
        //sqrt(pow(fabs(gentop_eta-genantitop_eta),2)+pow(TMath::Min(fabs(gentop_phi-genantitop_phi),2*Pi-fabs(gentop_phi-genantitop_phi)),2))
        
        for(int i=0;i<3;i++){
            matchingt[i]=-99;matchingw[i]=-99;unmatching[i]=-99;matchingz[i]=-99;matchingg[i]=-99;matchingr[i]=-99;matchingq[i]=-99;
            for (int j=0;j<9;j++){
                matchingtg[j][i]=-99;matchingtq[j][i]=-99;}
        }
        double dRtj[3]={99,99,99},dRantitj[3]={99,99,99},dRwj[5][3],dRzj[5][3],dRgj[10][3],dRrj[3]={99,99,99},dRtg[10],dRq1j[5][3],dRq2j[5][3],dRq3j[5][3],dRq4j[5][3],dRq5j[5][3],dRtq1[5],dRtq2[5],dRtq3[5],dRtq4[5],dRtq5[5];
        double dRw_q1j[5][3],dRw_q2j[5][3],dRw2j[3]={99,99,99},dRw3j[3]={99,99,99},matchingtmp[3]={-99,-99,-99};
        if (gentop_pt>0){
            dRtj[0] = deltaR(gentop_eta,gentop_phi,eta_max,phi_max);
            if(Nj8==3)dRtj[1] = deltaR(gentop_eta,gentop_phi,eta_mid,phi_mid);
            if(Nj8==2)dRtj[1] = deltaR(gentop_eta,gentop_phi,eta_min,phi_min);
            dRtj[2] = deltaR(gentop_eta,gentop_phi,eta_min,phi_min);}
        if(IsData==10&&gen_rad_pt>0){
            dRrj[0] = deltaR(gen_rad_eta,gen_rad_phi,eta_max,phi_max);
            if(Nj8==3) dRrj[1] = deltaR(gen_rad_eta,gen_rad_phi,eta_mid,phi_mid);
            if(Nj8==2) dRrj[1] = deltaR(gen_rad_eta,gen_rad_phi,eta_min,phi_min);
            dRrj[2] = deltaR(gen_rad_eta,gen_rad_phi,eta_min,phi_min);}
        if (genantitop_pt>0){
            dRantitj[0] = deltaR(genantitop_eta,genantitop_phi,eta_max,phi_max);
            if(Nj8==3) dRantitj[1] = deltaR(genantitop_eta,genantitop_phi,eta_mid,phi_mid);
            if(Nj8==2) dRantitj[1] = deltaR(genantitop_eta,genantitop_phi,eta_min,phi_min);
            dRantitj[2] = deltaR(genantitop_eta,genantitop_phi,eta_min,phi_min);}
        for(int i=0;i<5;i++){
            dRtq1[i]=99;dRtq2[i]=99;dRtq3[i]=99;dRtq4[i]=99;dRtq5[i]=99;
            for (int j=0;j<3;j++){
                dRwj[i][j]=99;
                dRzj[i][j]=99;
                dRq1j[i][j]=99;dRq2j[i][j]=99;dRq3j[i][j]=99;dRq4j[i][j]=99;dRq5j[i][j]=99;
                dRw_q1j[i][j]=99;dRw_q2j[i][j]=99;
            }
            if(ptgenwl[i]>0){
                dRwj[i][0] = deltaR(etagenwl[i],phigenwl[i],eta_max,phi_max);
                if(Nj8==3) dRwj[i][1] = deltaR(etagenwl[i],phigenwl[i],eta_mid,phi_mid);
                if(Nj8==2) dRwj[i][1] = deltaR(etagenwl[i],phigenwl[i],eta_min,phi_min);
                dRwj[i][2] = deltaR(etagenwl[i],phigenwl[i],eta_min,phi_min);}
            if(ptgenzl[i]>0){
                dRzj[i][0] = deltaR(etagenzl[i],phigenzl[i],eta_max,phi_max);
                if(Nj8==3) dRzj[i][1] = deltaR(etagenzl[i],phigenzl[i],eta_mid,phi_mid);
                if(Nj8==2) dRzj[i][1] = deltaR(etagenzl[i],phigenzl[i],eta_min,phi_min);
                dRzj[i][2] = deltaR(etagenzl[i],phigenzl[i],eta_min,phi_min);}
            if(ptgenq1l[i]>0){
                dRq1j[i][0] = deltaR(etagenq1l[i],phigenq1l[i],eta_max,phi_max);
                if(Nj8==3) dRq1j[i][1] = deltaR(etagenq1l[i],phigenq1l[i],eta_mid,phi_mid);
                if(Nj8==2) dRq1j[i][1] = deltaR(etagenq1l[i],phigenq1l[i],eta_min,phi_min);
                dRq1j[i][2] = deltaR(etagenq1l[i],phigenq1l[i],eta_min,phi_min);}
            if(ptgenq2l[i]>0){
                dRq2j[i][0] = deltaR(etagenq2l[i],phigenq2l[i],eta_max,phi_max);
                if(Nj8==3) dRq2j[i][1] = deltaR(etagenq2l[i],phigenq2l[i],eta_mid,phi_mid);
                if(Nj8==2) dRq2j[i][1] = deltaR(etagenq2l[i],phigenq2l[i],eta_min,phi_min);
                dRq2j[i][2] = deltaR(etagenq2l[i],phigenq2l[i],eta_min,phi_min);}
            if(ptgenq3l[i]>0){
                dRq3j[i][0] = deltaR(etagenq3l[i],phigenq3l[i],eta_max,phi_max);
                if(Nj8==3) dRq3j[i][1] = deltaR(etagenq3l[i],phigenq3l[i],eta_mid,phi_mid);
                if(Nj8==2) dRq3j[i][1] = deltaR(etagenq3l[i],phigenq3l[i],eta_min,phi_min);
                dRq3j[i][2] = deltaR(etagenq3l[i],phigenq3l[i],eta_min,phi_min);}
            if(ptgenq4l[i]>0){
                dRq4j[i][0] = deltaR(etagenq4l[i],phigenq4l[i],eta_max,phi_max);
                if(Nj8==3) dRq4j[i][1] = deltaR(etagenq4l[i],phigenq4l[i],eta_mid,phi_mid);
                if(Nj8==2) dRq4j[i][1] = deltaR(etagenq4l[i],phigenq4l[i],eta_min,phi_min);
                dRq4j[i][2] = deltaR(etagenq4l[i],phigenq4l[i],eta_min,phi_min);}
            if(ptgenq5l[i]>0){
                dRq5j[i][0] = deltaR(etagenq5l[i],phigenq5l[i],eta_max,phi_max);
                if(Nj8==3) dRq5j[i][1] = deltaR(etagenq5l[i],phigenq5l[i],eta_mid,phi_mid);
                if(Nj8==2) dRq5j[i][1] = deltaR(etagenq5l[i],phigenq5l[i],eta_min,phi_min);
                dRq5j[i][2] = deltaR(etagenq5l[i],phigenq5l[i],eta_min,phi_min);}

        }
        for(int i=0;i<10;i++){
            dRtg[i]=99;
            for (int j=0;j<3;j++){dRgj[i][j]=99;}
            if(ptgengl[i]>0){
                dRgj[i][0] = deltaR(etagengl[i],phigengl[i],eta_max,phi_max);
                if(Nj8==3) dRgj[i][1] = deltaR(etagengl[i],phigengl[i],eta_mid,phi_mid);
                if(Nj8==2) dRgj[i][1] = deltaR(etagengl[i],phigengl[i],eta_min,phi_min);
                dRgj[i][2] = deltaR(etagengl[i],phigengl[i],eta_min,phi_min);}
        }
        
        double dRt_bj[3]={99,99,99},dRt_wq1j[3]={99,99,99},dRt_wq2j[3]={99,99,99},dRt_wj[3]={99,99,99},dRt_wq1q2[3]={99,99,99};
        double tmpbpt=-99,tmpbm=-99,tmpbeta=-99,tmpbphi=-99,tmpq1pt=-99,tmpq1eta=-99,tmpq1phi=-99,tmpq2pt=-99,tmpq2eta=-99,tmpq2phi=-99,tmpwpt=-99,tmpweta=-99,tmpwphi=-99,tmptpt=-99,tmpteta=-99,tmptphi=-99,tmpt_wtag=-99;
        double dRr_q1j[3]={99,99,99},dRr_q2j[3]={99,99,99},dRr_q3j[3]={99,99,99},dRr_q4j[3]={99,99,99};
        double dRt1_bj[3]={99,99,99},dRt2_bj[3]={99,99,99},dRt1_wq1j[3]={99,99,99},dRt1_wq2j[3]={99,99,99},dRt2_wq1j[3]={99,99,99},dRt2_wq2j[3]={99,99,99},dRt_wb[3]={99,99,99},dRt1_wj[3]={99,99,99},dRt2_wj[3]={99,99,99},dRt1_wb[3]={99,99,99},dRt2_wb[3]={99,99,99};
        double rmatching=0.6,rbmatching=0.8;
        TLorentzVector t_b,matchedgv,matchedq1v,matchedq2v,matchedq3v,matchedq4v,matchedq5v;
        dRWW_R=-99;dRWb=-99;dRqq_W1=-99;dRqq_W2=-99;mindRqq=-99;maxdRqq=-99;
        double mindRqq_all=-99,maxdRqq_all=-99;
        //if (jentry==100000) break;
        //cout<<jentry<<endl;
        for(int i=0;i<Nj8;i++){
            //if (Nj8==2) break;
            if(gentop_pt>0||genantitop_pt>0){
                if (dRtj[i]<rmatching) matchingtmp[i]=-1;
                if (matchingtmp[i]!=-1&&dRantitj[i]<rmatching) matchingtmp[i]=-2;
                if (matchingtmp[i]==-1&&dRantitj[i]<rmatching&&dRtj[i]<dRantitj[i]) matchingtmp[i]=-1;
                if (matchingtmp[i]==-1&&dRantitj[i]<rmatching&&dRtj[i]>dRantitj[i]) matchingtmp[i]=-2;
                //cout<<i<<" jet  "<<Nj8<<"   "<<matchingtmp[i]<<endl;
                if(matchingtmp[i]==-1||matchingtmp[i]==-2){
                    if(matchingtmp[i]==-1){tmptpt=gentop_pt;tmpteta=gentop_eta;tmptphi=gentop_phi;tmpbpt=gent_b_pt;tmpbeta=gent_b_eta;tmpbphi=gent_b_phi;tmpwpt=gent_w_pt;tmpweta=gent_w_eta;tmpwphi=gent_w_phi;tmpq1pt=gent_w_q1_pt;tmpq1eta=gent_w_q1_eta;tmpq1phi=gent_w_q1_phi;tmpq2pt=gent_w_q2_pt;tmpq2eta=gent_w_q2_eta;tmpq2phi=gent_w_q2_phi;tmpbm=gent_b_mass;tmpt_wtag=gent_w_tag;}
                    if(matchingtmp[i]==-2){tmptpt=genantitop_pt;tmpteta=genantitop_eta;tmptphi=genantitop_phi;tmpbpt=genantit_b_pt;tmpbeta=genantit_b_eta;tmpbphi=genantit_b_phi;tmpwpt=genantit_w_pt;tmpweta=genantit_w_eta;tmpwphi=genantit_w_phi;tmpq1pt=genantit_w_q1_pt;tmpq1eta=genantit_w_q1_eta;tmpq1phi=genantit_w_q1_phi;tmpq2pt=genantit_w_q2_pt;tmpq2eta=genantit_w_q2_eta;tmpq2phi=genantit_w_q2_phi;tmpbm=genantit_b_mass;tmpt_wtag=genantit_w_tag;}
                        if(tmpbpt<0) matchingt[i]=-1;
                        if(tmpbpt>0){
                            if(i==0) dRt_bj[i]=deltaR(tmpbeta,tmpbphi,eta_max,phi_max);
                            if(i==1&&Nj8!=2) dRt_bj[i]=deltaR(tmpbeta,tmpbphi,eta_mid,phi_mid);
                            if(i==1&&Nj8==2) dRt_bj[i]=deltaR(tmpbeta,tmpbphi,eta_min,phi_min);
                            if(i==2) dRt_bj[i]=deltaR(tmpbeta,tmpbphi,eta_min,phi_min);
                            if(i==0) dRt_wq1j[i]=deltaR(tmpq1eta,tmpq1phi,eta_max,phi_max);
                            if(i==1&&Nj8!=2) dRt_wq1j[i]=deltaR(tmpq1eta,tmpq1phi,eta_mid,phi_mid);
                            if(i==1&&Nj8==2) dRt_wq1j[i]=deltaR(tmpq1eta,tmpq1phi,eta_min,phi_min);
                            if(i==2) dRt_wq1j[i]=deltaR(tmpq1eta,tmpq1phi,eta_min,phi_min);
                            if(i==0) dRt_wq2j[i]=deltaR(tmpq2eta,tmpq2phi,eta_max,phi_max);
                            if(i==1&&Nj8!=2) dRt_wq2j[i]=deltaR(tmpq2eta,tmpq2phi,eta_mid,phi_mid);
                            if(i==1&&Nj8==2) dRt_wq2j[i]=deltaR(tmpq2eta,tmpq2phi,eta_min,phi_min);
                            if(i==2) dRt_wq2j[i]=deltaR(tmpq2eta,tmpq2phi,eta_min,phi_min);
                            if(i==0) dRt_wj[i]=deltaR(tmpweta,tmpwphi,eta_max,phi_max);
                            if(i==1&&Nj8!=2) dRt_wj[i]=deltaR(tmpweta,tmpwphi,eta_mid,phi_mid);
                            if(i==1&&Nj8==2) dRt_wj[i]=deltaR(tmpweta,tmpwphi,eta_min,phi_min);
                            if(i==2) dRt_wj[i]=deltaR(tmpweta,tmpwphi,eta_min,phi_min);
                            dRt_wb[i]=deltaR(tmpweta,tmpwphi,tmpbeta,tmpbphi);
                            if(i==0&&Mj_max>100){
                                dRWb=deltaR(tmpweta,tmpwphi,tmpbeta,tmpbphi);
                                dRqq_W1=deltaR(tmpq1eta,tmpq1phi,tmpq2eta,tmpq2phi);
                            }
                            //if(dRt_bj[i]<rbmatching&&dRt_wq1j[i]<rbmatching&&dRt_wq2j[i]<rbmatching&&dRt_wj[i]<rmatching&&tmpt_wtag==4)
                            if(dRt_bj[i]<rbmatching&&dRt_wq1j[i]<rbmatching&&dRt_wq2j[i]<rbmatching&&tmpt_wtag==4){
                                matchingt[i]=3;
                                for(int j=0;j<10;j++){
                                    dRtg[j]=deltaR(tmpteta,tmptphi,etagengl[j],phigengl[j]);
                                    if(dRtg[j]<rbmatching&&dRgj[j][i]<rbmatching) matchingt[i]=4;
                                    t_b.SetPtEtaPhiM(tmpbpt,tmpbeta,tmpbphi,tmpbm);
                                    matchedgv.SetPtEtaPhiE(ptgengl[j],etagengl[j],phigengl[j],egengl[j]);
                                    if(matchingt[i]==4&&(t_b+matchedgv).Mag()>40) matchingtg[0][i]=4;
                                    if(matchingt[i]==4&&(t_b+matchedgv).Mag()>60) matchingtg[1][i]=4;
                                    if(matchingt[i]==4&&(t_b+matchedgv).Mag()>40&&(t_b+matchedgv).Mag()<120) matchingtg[2][i]=4;
                                    if(matchingt[i]==4&&(t_b+matchedgv).Mag()>60&&(t_b+matchedgv).Mag()<100) matchingtg[3][i]=4;
                                    if(matchingt[i]==4&&ptgengl[j]>0.1*tmptpt) matchingtg[4][i]=4;
                                    if(matchingt[i]==4&&ptgengl[j]>0.15*tmptpt) matchingtg[5][i]=4;
                                    if(matchingt[i]==4&&ptgengl[j]>0.2*tmptpt) matchingtg[6][i]=4;
                                    if(matchingt[i]==4&&ptgengl[j]>0.25*tmptpt) matchingtg[7][i]=4;
                                    if(matchingt[i]==4&&(t_b+matchedgv).Pt()>0.2*tmpwpt) matchingtg[8][i]=4;
                                }
                                for(int j=0;j<5;j++){
                                    dRtq1[j]=deltaR(tmpteta,tmptphi,etagenq1l[j],phigenq1l[j]);
                                    dRtq2[j]=deltaR(tmpteta,tmptphi,etagenq2l[j],phigenq2l[j]);
                                    dRtq3[j]=deltaR(tmpteta,tmptphi,etagenq3l[j],phigenq3l[j]);
                                    dRtq4[j]=deltaR(tmpteta,tmptphi,etagenq4l[j],phigenq4l[j]);
                                    dRtq5[j]=deltaR(tmpteta,tmptphi,etagenq5l[j],phigenq5l[j]);
                                    if(matchingt[i]!=4&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching))) matchingt[i]=6;
                                    if(matchingt[i]==4&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching))) matchingt[i]=46;
                                    t_b.SetPtEtaPhiM(tmpbpt,tmpbeta,tmpbphi,tmpbm);
                                    matchedq1v.SetPtEtaPhiE(ptgenq1l[j],etagenq1l[j],phigenq1l[j],egenq1l[j]);
                                    matchedq2v.SetPtEtaPhiE(ptgenq2l[j],etagenq2l[j],phigenq2l[j],egenq2l[j]);
                                    matchedq3v.SetPtEtaPhiE(ptgenq3l[j],etagenq3l[j],phigenq3l[j],egenq3l[j]);
                                    matchedq4v.SetPtEtaPhiE(ptgenq4l[j],etagenq4l[j],phigenq4l[j],egenq4l[j]);
                                    matchedq5v.SetPtEtaPhiE(ptgenq5l[j],etagenq5l[j],phigenq5l[j],egenq5l[j]);
                                    if((matchingt[i]==6||matchingt[i]==46)&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&(t_b+matchedq1v).Mag()>40)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&(t_b+matchedq2v).Mag()>40)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&(t_b+matchedq3v).Mag()>40)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&(t_b+matchedq4v).Mag()>40)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&(t_b+matchedq5v).Mag()>40))) matchingtq[0][i]=6;
                                    if((matchingt[i]==6||matchingt[i]==46)&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&(t_b+matchedq1v).Mag()>60)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&(t_b+matchedq2v).Mag()>60)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&(t_b+matchedq3v).Mag()>60)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&(t_b+matchedq4v).Mag()>60)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&(t_b+matchedq5v).Mag()>60))) matchingtq[1][i]=6;
                                    if((matchingt[i]==6||matchingt[i]==46)&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&(t_b+matchedq1v).Mag()>40&&(t_b+matchedq1v).Mag()<120)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&(t_b+matchedq2v).Mag()>40&&(t_b+matchedq2v).Mag()<120)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&(t_b+matchedq3v).Mag()>40&&(t_b+matchedq3v).Mag()<120)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&(t_b+matchedq4v).Mag()>40&&(t_b+matchedq4v).Mag()<120)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&(t_b+matchedq5v).Mag()>40&&(t_b+matchedq5v).Mag()<120))) matchingtq[2][i]=6;
                                    if((matchingt[i]==6||matchingt[i]==46)&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&(t_b+matchedq1v).Mag()>60&&(t_b+matchedq1v).Mag()<100)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&(t_b+matchedq2v).Mag()>60&&(t_b+matchedq2v).Mag()<100)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&(t_b+matchedq3v).Mag()>60&&(t_b+matchedq3v).Mag()<100)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&(t_b+matchedq4v).Mag()>60&&(t_b+matchedq4v).Mag()<100)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&(t_b+matchedq5v).Mag()>60&&(t_b+matchedq5v).Mag()<100))) matchingtq[3][i]=6;
                                    if((matchingt[i]==6||matchingt[i]==46)&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&ptgenq1l[j]>0.1*tmptpt)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&ptgenq2l[j]>0.1*tmptpt)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&ptgenq3l[j]>0.1*tmptpt)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&ptgenq4l[j]>0.1*tmptpt)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&ptgenq5l[j]>0.1*tmptpt))) matchingtq[4][i]=6;
                                    if((matchingt[i]==6||matchingt[i]==46)&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&ptgenq1l[j]>0.15*tmptpt)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&ptgenq2l[j]>0.15*tmptpt)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&ptgenq3l[j]>0.15*tmptpt)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&ptgenq4l[j]>0.15*tmptpt)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&ptgenq5l[j]>0.15*tmptpt))) matchingtq[5][i]=6;
                                    if((matchingt[i]==6||matchingt[i]==46)&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&ptgenq1l[j]>0.2*tmptpt)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&ptgenq2l[j]>0.2*tmptpt)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&ptgenq3l[j]>0.2*tmptpt)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&ptgenq4l[j]>0.2*tmptpt)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&ptgenq5l[j]>0.2*tmptpt))) matchingtq[6][i]=6;
                                    if((matchingt[i]==6||matchingt[i]==46)&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&ptgenq1l[j]>0.25*tmptpt)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&ptgenq2l[j]>0.25*tmptpt)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&ptgenq3l[j]>0.25*tmptpt)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&ptgenq4l[j]>0.25*tmptpt)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&ptgenq5l[j]>0.25*tmptpt))) matchingtq[7][i]=6;
                                    if((matchingt[i]==6||matchingt[i]==46)&&(((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&(t_b+matchedq1v).Pt()>0.2*tmpwpt)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&(t_b+matchedq2v).Pt()>0.2*tmpwpt)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&(t_b+matchedq3v).Pt()>0.2*tmpwpt)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&(t_b+matchedq4v).Pt()>0.2*tmpwpt)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&(t_b+matchedq5v).Pt()>0.2*tmpwpt)))) matchingtq[8][i]=6;
                                }
                            }
                            if(dRt_bj[i]<rbmatching&&((dRt_wq1j[i]<rbmatching&&dRt_wq2j[i]>rbmatching)||(dRt_wq1j[i]>rbmatching&&dRt_wq2j[i]<rbmatching))&&tmpt_wtag==4){
                                matchingt[i]=7;
                            }
                            //cout<<tmpt_wtag<<" w tag"<<endl;
                            if(dRt_bj[i]<rbmatching&&dRt_wq1j[i]>rbmatching&&dRt_wq2j[i]>rbmatching&&tmpt_wtag==4){
                                matchingt[i]=-3;
                                if(matchingq[i]!=-99) matchingq[i]+=1;
                                if(matchingq[i]==-99) matchingq[i]=1;}
                            if(dRt_bj[i]<rbmatching&&tmpt_wtag!=4&&dRt_wq1j[i]<rbmatching){
                                //cout<<"   end"<<endl;
                                //cout<<"wlnu  "<<dRt_wq1j[i]<<endl;
                                matchingt[i]=1;
                                for(int j=0;j<10;j++){
                                    dRtg[j]=deltaR(tmpteta,tmptphi,etagengl[j],phigengl[j]);
                                    if(dRtg[j]<rbmatching&&dRgj[j][i]<rbmatching) matchingt[i]=2;
                                    t_b.SetPtEtaPhiM(tmpbpt,tmpbeta,tmpbphi,tmpbm);
                                    matchedgv.SetPtEtaPhiE(ptgengl[j],etagengl[j],phigengl[j],egengl[j]);
                                    if(matchingt[i]==2&&(t_b+matchedgv).Mag()>40) matchingtg[0][i]=2;
                                    if(matchingt[i]==2&&(t_b+matchedgv).Mag()>60) matchingtg[1][i]=2;
                                    if(matchingt[i]==2&&(t_b+matchedgv).Mag()>40&&(t_b+matchedgv).Mag()<120) matchingtg[2][i]=2;
                                    if(matchingt[i]==2&&(t_b+matchedgv).Mag()>60&&(t_b+matchedgv).Mag()<100) matchingtg[3][i]=2;
                                    if(matchingt[i]==2&&ptgengl[j]>0.1*tmptpt) matchingtg[4][i]=2;
                                    if(matchingt[i]==2&&ptgengl[j]>0.15*tmptpt) matchingtg[5][i]=2;
                                    if(matchingt[i]==2&&ptgengl[j]>0.2*tmptpt) matchingtg[6][i]=2;
                                    if(matchingt[i]==2&&ptgengl[j]>0.25*tmptpt) matchingtg[7][i]=2;
                                    if(matchingt[i]==2&&(t_b+matchedgv).Pt()>0.2*tmpwpt) matchingtg[8][i]=2;
                                }
                                for(int j=0;j<5;j++){
                                    dRtq1[j]=deltaR(tmpteta,tmptphi,etagenq1l[j],phigenq1l[j]);
                                    dRtq2[j]=deltaR(tmpteta,tmptphi,etagenq2l[j],phigenq2l[j]);
                                    dRtq3[j]=deltaR(tmpteta,tmptphi,etagenq3l[j],phigenq3l[j]);
                                    dRtq4[j]=deltaR(tmpteta,tmptphi,etagenq4l[j],phigenq4l[j]);
                                    dRtq5[j]=deltaR(tmpteta,tmptphi,etagenq5l[j],phigenq5l[j]);
                                    if(matchingt[i]!=2&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching))) matchingt[i]=5;
                                    if(matchingt[i]==2&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching))) matchingt[i]=25;
                                    t_b.SetPtEtaPhiM(tmpbpt,tmpbeta,tmpbphi,tmpbm);
                                    matchedq1v.SetPtEtaPhiE(ptgenq1l[j],etagenq1l[j],phigenq1l[j],egenq1l[j]);
                                    matchedq2v.SetPtEtaPhiE(ptgenq2l[j],etagenq2l[j],phigenq2l[j],egenq2l[j]);
                                    matchedq3v.SetPtEtaPhiE(ptgenq3l[j],etagenq3l[j],phigenq3l[j],egenq3l[j]);
                                    matchedq4v.SetPtEtaPhiE(ptgenq4l[j],etagenq4l[j],phigenq4l[j],egenq4l[j]);
                                    matchedq5v.SetPtEtaPhiE(ptgenq5l[j],etagenq5l[j],phigenq5l[j],egenq5l[j]);
                                    if((matchingt[i]==5||matchingt[i]==25)&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&(t_b+matchedq1v).Mag()>40)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&(t_b+matchedq2v).Mag()>40)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&(t_b+matchedq3v).Mag()>40)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&(t_b+matchedq4v).Mag()>40)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&(t_b+matchedq5v).Mag()>40))) matchingtq[0][i]=5;
                                    if((matchingt[i]==5||matchingt[i]==25)&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&(t_b+matchedq1v).Mag()>60)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&(t_b+matchedq2v).Mag()>60)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&(t_b+matchedq3v).Mag()>60)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&(t_b+matchedq4v).Mag()>60)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&(t_b+matchedq5v).Mag()>60))) matchingtq[1][i]=5;
                                    if((matchingt[i]==5||matchingt[i]==25)&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&(t_b+matchedq1v).Mag()>40&&(t_b+matchedq1v).Mag()<120)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&(t_b+matchedq2v).Mag()>40&&(t_b+matchedq2v).Mag()<120)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&(t_b+matchedq3v).Mag()>40&&(t_b+matchedq3v).Mag()<120)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&(t_b+matchedq4v).Mag()>40&&(t_b+matchedq4v).Mag()<120)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&(t_b+matchedq5v).Mag()>40&&(t_b+matchedq5v).Mag()<120))) matchingtq[2][i]=5;
                                    if((matchingt[i]==5||matchingt[i]==25)&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&(t_b+matchedq1v).Mag()>60&&(t_b+matchedq1v).Mag()<100)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&(t_b+matchedq2v).Mag()>60&&(t_b+matchedq2v).Mag()<100)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&(t_b+matchedq3v).Mag()>60&&(t_b+matchedq3v).Mag()<100)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&(t_b+matchedq4v).Mag()>60&&(t_b+matchedq4v).Mag()<100)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&(t_b+matchedq5v).Mag()>60&&(t_b+matchedq5v).Mag()<100))) matchingtq[3][i]=5;
                                    if((matchingt[i]==5||matchingt[i]==25)&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&ptgenq1l[j]>0.1*tmptpt)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&ptgenq2l[j]>0.1*tmptpt)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&ptgenq3l[j]>0.1*tmptpt)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&ptgenq4l[j]>0.1*tmptpt)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&ptgenq5l[j]>0.1*tmptpt))) matchingtq[4][i]=5;
                                    if((matchingt[i]==5||matchingt[i]==25)&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&ptgenq1l[j]>0.15*tmptpt)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&ptgenq2l[j]>0.15*tmptpt)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&ptgenq3l[j]>0.15*tmptpt)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&ptgenq4l[j]>0.15*tmptpt)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&ptgenq5l[j]>0.15*tmptpt))) matchingtq[5][i]=5;
                                    if((matchingt[i]==5||matchingt[i]==25)&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&ptgenq1l[j]>0.2*tmptpt)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&ptgenq2l[j]>0.2*tmptpt)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&ptgenq3l[j]>0.2*tmptpt)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&ptgenq4l[j]>0.2*tmptpt)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&ptgenq5l[j]>0.2*tmptpt))) matchingtq[6][i]=5;
                                    if((matchingt[i]==5||matchingt[i]==25)&&((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&ptgenq1l[j]>0.25*tmptpt)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&ptgenq2l[j]>0.25*tmptpt)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&ptgenq3l[j]>0.25*tmptpt)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&ptgenq4l[j]>0.25*tmptpt)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&ptgenq5l[j]>0.25*tmptpt))) matchingtq[7][i]=5;
                                    if((matchingt[i]==5||matchingt[i]==25)&&(((dRtq1[j]<rbmatching&&dRq1j[j][i]<rbmatching&&(t_b+matchedq1v).Pt()>0.2*tmpwpt)||(dRtq2[j]<rbmatching&&dRq2j[j][i]<rbmatching&&(t_b+matchedq2v).Pt()>0.2*tmpwpt)||(dRtq3[j]<rbmatching&&dRq3j[j][i]<rbmatching&&(t_b+matchedq3v).Pt()>0.2*tmpwpt)||(dRtq4[j]<rbmatching&&dRq4j[j][i]<rbmatching&&(t_b+matchedq4v).Pt()>0.2*tmpwpt)||(dRtq5[j]<rbmatching&&dRq5j[j][i]<rbmatching&&(t_b+matchedq5v).Pt()>0.2*tmpwpt)))) matchingtq[8][i]=5;

                                }
                            }
                            dRt_wq1q2[i]=deltaR(tmpq1eta,tmpq1phi,tmpq2eta,tmpq2phi);
                            if(dRt_bj[i]>rbmatching&&dRt_wq1j[i]<rbmatching&&dRt_wq2j[i]<rbmatching&&dRt_wj[i]<rmatching&&tmpt_wtag==4){
                                matchingt[i]=-2;
                                if(matchingw[i]!=-99) matchingw[i]+=1;
                                if(matchingw[i]==-99) matchingw[i]=1;
                                for(int j=0;j<10;j++){
                                    if(dRgj[j][i]<rbmatching&&ptgengl[j]>0.1*tmpwpt) matchingw[i]=100;
                                    if(dRgj[j][i]<rbmatching&&ptgengl[j]>0.15*tmpwpt) matchingw[i]=101;
                                }
                                for(int j=0;j<5;j++){
                                    if((dRq1j[j][i]<rbmatching&&ptgenq1l[j]<0.1*tmpwpt)||(dRq2j[j][i]<rbmatching&&ptgenq2l[j]<0.1*tmpwpt)||(dRq3j[j][i]<rbmatching&&ptgenq3l[j]<0.1*tmpwpt)||(dRq4j[j][i]<rbmatching&&ptgenq4l[j]<0.1*tmpwpt)||(dRq5j[j][i]<rbmatching&&ptgenq5l[j]<0.1*tmpwpt)) matchingw[i]=102;
                                    if((dRq1j[j][i]<rbmatching&&ptgenq1l[j]>0.15*tmpwpt)||(dRq2j[j][i]<rbmatching&&ptgenq2l[j]>0.15*tmpwpt)||(dRq3j[j][i]<rbmatching&&ptgenq3l[j]>0.15*tmpwpt)||(dRq4j[j][i]<rbmatching&&ptgenq4l[j]>0.15*tmpwpt)||(dRq5j[j][i]<rbmatching&&ptgenq5l[j]>0.15*tmpwpt)) matchingw[i]=103;
                                }

                            }
                            if(dRt_bj[i]>rbmatching&&((dRt_wq1j[i]>rbmatching&&dRt_wq2j[i]<rbmatching)||(dRt_wq1j[i]<rbmatching&&dRt_wq2j[i]>rbmatching))&&tmpt_wtag==4){
                                matchingt[i]=-4;
                                if(matchingq[i]!=-99) matchingq[i]+=1;
                                if(matchingq[i]==-99) matchingq[i]=1;
                            }
                        }
                }
                if(!(matchingtmp[i]==-1||matchingtmp[i]==-2)){
                    if(i==0) dRt1_bj[i]=deltaR(gent_b_eta,gent_b_phi,eta_max,phi_max);
                    if(i==1&&Nj8!=2) dRt1_bj[i]=deltaR(gent_b_eta,gent_b_phi,eta_mid,phi_mid);
                    if(i==1&&Nj8==2) dRt1_bj[i]=deltaR(gent_b_eta,gent_b_phi,eta_min,phi_min);
                    if(i==2) dRt1_bj[i]=deltaR(gent_b_eta,gent_b_phi,eta_min,phi_min);
                    if(i==0) dRt2_bj[i]=deltaR(genantit_b_eta,genantit_b_phi,eta_max,phi_max);
                    if(i==1&&Nj8!=2) dRt2_bj[i]=deltaR(genantit_b_eta,genantit_b_phi,eta_mid,phi_mid);
                    if(i==1&&Nj8==2) dRt2_bj[i]=deltaR(genantit_b_eta,genantit_b_phi,eta_min,phi_min);
                    if(i==2) dRt2_bj[i]=deltaR(genantit_b_eta,genantit_b_phi,eta_min,phi_min);
                    
                    if(i==0) dRt1_wq1j[i]=deltaR(gent_w_q1_eta,gent_w_q1_phi,eta_max,phi_max);
                    if(i==1&&Nj8!=2) dRt1_wq1j[i]=deltaR(gent_w_q1_eta,gent_w_q1_phi,eta_mid,phi_mid);
                    if(i==1&&Nj8==2) dRt1_wq1j[i]=deltaR(gent_w_q1_eta,gent_w_q1_phi,eta_min,phi_min);
                    if(i==2) dRt1_wq1j[i]=deltaR(gent_w_q1_eta,gent_w_q1_phi,eta_min,phi_min);
                    if(i==0) dRt2_wq1j[i]=deltaR(genantit_w_q1_eta,genantit_w_q1_phi,eta_max,phi_max);
                    if(i==1&&Nj8!=2) dRt2_wq1j[i]=deltaR(genantit_w_q1_eta,genantit_w_q1_phi,eta_mid,phi_mid);
                    if(i==1&&Nj8==2) dRt2_wq1j[i]=deltaR(genantit_w_q1_eta,genantit_w_q1_phi,eta_min,phi_min);
                    if(i==2) dRt2_wq1j[i]=deltaR(genantit_w_q1_eta,genantit_w_q1_phi,eta_min,phi_min);
                    if(i==0) dRt1_wq2j[i]=deltaR(gent_w_q2_eta,gent_w_q2_phi,eta_max,phi_max);
                    if(i==1&&Nj8!=2) dRt1_wq2j[i]=deltaR(gent_w_q2_eta,gent_w_q2_phi,eta_mid,phi_mid);
                    if(i==1&&Nj8==2) dRt1_wq2j[i]=deltaR(gent_w_q2_eta,gent_w_q2_phi,eta_min,phi_min);
                    if(i==2) dRt1_wq2j[i]=deltaR(gent_w_q2_eta,gent_w_q2_phi,eta_min,phi_min);
                    if(i==0) dRt2_wq2j[i]=deltaR(genantit_w_q2_eta,genantit_w_q2_phi,eta_max,phi_max);
                    if(i==1&&Nj8!=2) dRt2_wq2j[i]=deltaR(genantit_w_q2_eta,genantit_w_q2_phi,eta_mid,phi_mid);
                    if(i==1&&Nj8==2) dRt2_wq2j[i]=deltaR(genantit_w_q2_eta,genantit_w_q2_phi,eta_min,phi_min);
                    if(i==2) dRt2_wq2j[i]=deltaR(genantit_w_q2_eta,genantit_w_q2_phi,eta_min,phi_min);
                    
                    if(i==0) dRt1_wj[i]=deltaR(gent_w_eta,gent_w_phi,eta_max,phi_max);
                    if(i==1&&Nj8!=2) dRt1_wj[i]=deltaR(gent_w_eta,gent_w_phi,eta_mid,phi_mid);
                    if(i==1&&Nj8==2) dRt1_wj[i]=deltaR(gent_w_eta,gent_w_phi,eta_min,phi_min);
                    if(i==2) dRt1_wj[i]=deltaR(gent_w_eta,gent_w_phi,eta_min,phi_min);
                    if(i==0) dRt2_wj[i]=deltaR(genantit_w_eta,genantit_w_phi,eta_max,phi_max);
                    if(i==1&&Nj8!=2) dRt2_wj[i]=deltaR(genantit_w_eta,genantit_w_phi,eta_mid,phi_mid);
                    if(i==1&&Nj8==2) dRt2_wj[i]=deltaR(genantit_w_eta,genantit_w_phi,eta_min,phi_min);
                    if(i==2) dRt2_wj[i]=deltaR(genantit_w_eta,genantit_w_phi,eta_min,phi_min);
                    dRt1_wb[i]=deltaR(gent_w_eta,gent_w_phi,gent_b_eta,gent_b_phi);
                    dRt2_wb[i]=deltaR(genantit_w_eta,genantit_w_phi,genantit_b_eta,genantit_b_phi);
                    
                    if((dRt1_wq1j[i]<rbmatching&&dRt1_wq2j[i]<rbmatching&&dRt1_bj[i]>rbmatching&&dRt1_wj[i]<rmatching&&gent_w_tag==4)) {
                        matchingt[i]=-5;
                        if(matchingw[i]!=-99) matchingw[i]+=1;
                        if(matchingw[i]==-99) matchingw[i]=1;
                        for(int j=0;j<10;j++){
                            if(dRgj[j][i]<rbmatching&&ptgengl[j]>0.1*gent_w_pt) matchingw[i]=100;
                            if(dRgj[j][i]<rbmatching&&ptgengl[j]>0.15*gent_w_pt) matchingw[i]=101;
                        }
                        for(int j=0;j<5;j++){
                            if((dRq1j[j][i]<rbmatching&&ptgenq1l[j]<0.1*gent_w_pt)||(dRq2j[j][i]<rbmatching&&ptgenq2l[j]<0.1*gent_w_pt)||(dRq3j[j][i]<rbmatching&&ptgenq3l[j]<0.1*gent_w_pt)||(dRq4j[j][i]<rbmatching&&ptgenq4l[j]<0.1*gent_w_pt)||(dRq5j[j][i]<rbmatching&&ptgenq5l[j]<0.1*gent_w_pt)) matchingw[i]=102;
                            if((dRq1j[j][i]<rbmatching&&ptgenq1l[j]>0.15*gent_w_pt)||(dRq2j[j][i]<rbmatching&&ptgenq2l[j]>0.15*gent_w_pt)||(dRq3j[j][i]<rbmatching&&ptgenq3l[j]>0.15*gent_w_pt)||(dRq4j[j][i]<rbmatching&&ptgenq4l[j]>0.15*gent_w_pt)||(dRq5j[j][i]<rbmatching&&ptgenq5l[j]>0.15*gent_w_pt)) matchingw[i]=103;
                        }

                    }
                    if((dRt2_wq1j[i]<rbmatching&&dRt2_wq2j[i]<rbmatching&&dRt2_bj[i]>rbmatching&&dRt2_wj[i]<rmatching&&genantit_w_tag==4)) {
                        matchingt[i]=-6;
                        if(matchingw[i]!=-99) matchingw[i]+=1;
                        if(matchingw[i]==-99) matchingw[i]=1;
                        for(int j=0;j<10;j++){
                            if(dRgj[j][i]<rbmatching&&ptgengl[j]>0.1*genantit_w_pt) matchingw[i]=100;
                            if(dRgj[j][i]<rbmatching&&ptgengl[j]>0.15*genantit_w_pt) matchingw[i]=101;
                        }
                        for(int j=0;j<5;j++){
                            if((dRq1j[j][i]<rbmatching&&ptgenq1l[j]<0.1*genantit_w_pt)||(dRq2j[j][i]<rbmatching&&ptgenq2l[j]<0.1*genantit_w_pt)||(dRq3j[j][i]<rbmatching&&ptgenq3l[j]<0.1*genantit_w_pt)||(dRq4j[j][i]<rbmatching&&ptgenq4l[j]<0.1*genantit_w_pt)||(dRq5j[j][i]<rbmatching&&ptgenq5l[j]<0.1*genantit_w_pt)) matchingw[i]=102;
                            if((dRq1j[j][i]<rbmatching&&ptgenq1l[j]>0.15*genantit_w_pt)||(dRq2j[j][i]<rbmatching&&ptgenq2l[j]>0.15*genantit_w_pt)||(dRq3j[j][i]<rbmatching&&ptgenq3l[j]>0.15*genantit_w_pt)||(dRq4j[j][i]<rbmatching&&ptgenq4l[j]>0.15*genantit_w_pt)||(dRq5j[j][i]<rbmatching&&ptgenq5l[j]>0.15*genantit_w_pt)) matchingw[i]=103;
                        }

                    }

                    if((dRt1_wq1j[i]<rbmatching&&dRt1_wq2j[i]>rbmatching&&dRt1_bj[i]>rbmatching&&gent_w_tag==4)||(dRt1_wq1j[i]>rbmatching&&dRt1_wq2j[i]<rbmatching&&dRt1_bj[i]>rbmatching&&gent_w_tag==4)||(dRt1_wq1j[i]>rbmatching&&dRt1_wq2j[i]>rbmatching&&dRt1_bj[i]<rbmatching&&gent_w_tag==4)) {
                        matchingt[i]=-7;
                        if(matchingq[i]!=-99) matchingq[i]+=1;
                        if(matchingq[i]==-99) matchingq[i]=1;
                    }
                    if((dRt2_wq1j[i]<rbmatching&&dRt2_wq2j[i]>rbmatching&&dRt2_bj[i]>rbmatching&&genantit_w_tag==4)||(dRt2_wq1j[i]>rbmatching&&dRt2_wq2j[i]<rbmatching&&dRt2_bj[i]>rbmatching&&genantit_w_tag==4)||(dRt2_wq1j[i]>rbmatching&&dRt2_wq2j[i]>rbmatching&&dRt2_bj[i]<rbmatching&&genantit_w_tag==4)) {
                        matchingt[i]=-8;
                        if(matchingq[i]!=-99) matchingq[i]+=1;
                        if(matchingq[i]==-99) matchingq[i]=1;
                    }

                }
                for(int j=0;j<5;j++){
                    if (dRwj[j][i]<rmatching&&ptgenwf[j]!=gent_w_pt&&ptgenwf[j]!=genantit_w_pt){
                        if(i==0) dRw_q1j[j][i]=deltaR(genw_q1_eta[j],genw_q1_phi[j],eta_max,phi_max);
                        if(i==1&&Nj8!=2) dRw_q1j[j][i]=deltaR(genw_q1_eta[j],genw_q1_phi[j],eta_mid,phi_mid);
                        if(i==1&&Nj8==2) dRw_q1j[j][i]=deltaR(genw_q1_eta[j],genw_q1_phi[j],eta_min,phi_min);
                        if(i==2) dRw_q1j[j][i]=deltaR(genw_q1_eta[j],genw_q1_phi[j],eta_min,phi_min);
                        if(i==0) dRw_q2j[j][i]=deltaR(genw_q2_eta[j],genw_q2_phi[j],eta_max,phi_max);
                        if(i==1&&Nj8!=2) dRw_q2j[j][i]=deltaR(genw_q2_eta[j],genw_q2_phi[j],eta_mid,phi_mid);
                        if(i==1&&Nj8==2) dRw_q2j[j][i]=deltaR(genw_q2_eta[j],genw_q2_phi[j],eta_min,phi_min);
                        if(i==2) dRw_q2j[j][i]=deltaR(genw_q2_eta[j],genw_q2_phi[j],eta_min,phi_min);
                        if(dRw_q1j[j][i]<rbmatching&&dRw_q2j[j][i]<rbmatching&&taggenwl[j]==4) {
                            if(matchingw[i]!=-99) matchingw[i]+=1;
                            if(matchingw[i]==-99) matchingw[i]=1;
                            for(int j=0;j<10;j++){
                                if(dRgj[j][i]<rbmatching&&ptgengl[j]>0.1*gent_w_pt) matchingw[i]=100;
                                if(dRgj[j][i]<rbmatching&&ptgengl[j]>0.15*gent_w_pt) matchingw[i]=101;
                            }
                            for(int j=0;j<5;j++){
                                if((dRq1j[j][i]<rbmatching&&ptgenq1l[j]<0.1*gent_w_pt)||(dRq2j[j][i]<rbmatching&&ptgenq2l[j]<0.1*gent_w_pt)||(dRq3j[j][i]<rbmatching&&ptgenq3l[j]<0.1*gent_w_pt)||(dRq4j[j][i]<rbmatching&&ptgenq4l[j]<0.1*gent_w_pt)||(dRq5j[j][i]<rbmatching&&ptgenq5l[j]<0.1*gent_w_pt)) matchingw[i]=102;
                                if((dRq1j[j][i]<rbmatching&&ptgenq1l[j]>0.15*gent_w_pt)||(dRq2j[j][i]<rbmatching&&ptgenq2l[j]>0.15*gent_w_pt)||(dRq3j[j][i]<rbmatching&&ptgenq3l[j]>0.15*gent_w_pt)||(dRq4j[j][i]<rbmatching&&ptgenq4l[j]>0.15*gent_w_pt)||(dRq5j[j][i]<rbmatching&&ptgenq5l[j]>0.15*gent_w_pt)) matchingw[i]=103;
                            }

                        }
                        if(((dRw_q1j[j][i]<rbmatching&&dRw_q2j[j][i]>rbmatching)||(dRw_q1j[j][i]>rbmatching&&dRw_q2j[j][i]<rbmatching))&&taggenwl[j]==4) {
                            if(matchingq[i]!=-99) matchingq[i]+=1;
                            if(matchingq[i]==-99) matchingq[i]=1;
                        }
                    }
                }

                    for(int j=0;j<5;j++){
                        if (dRzj[j][i]<rmatching){
                            if(matchingz[i]!=-99) matchingz[i]+=1;
                            if(matchingz[i]==-99) matchingz[i]=1;
                        }
                    }
                    for(int j=0;j<10;j++){
                        if (dRgj[j][i]<rmatching){
                            if(matchingg[i]!=-99) matchingg[i]+=1;
                            if(matchingg[i]==-99) matchingg[i]=1;
                        }
                    }
                    for(int j=0;j<5;j++){
                        if ((dRq1j[j][i]<rmatching)||(dRq2j[j][i]<rmatching)||(dRq3j[j][i]<rmatching)||(dRq4j[j][i]<rmatching)||(dRq5j[j][i]<rmatching)){
                            if(matchingq[i]!=-99) matchingq[i]+=1;
                            if(matchingq[i]==-99) matchingq[i]=1;
                        }
                    }
                //if((matchingtmp[i]==-1||matchingtmp[i]==-2)) cout<<i<<"  matched t "<<dRtj[i]<<"   "<<dRantitj[i]<<"  "<<dRt_bj[i]<<"   "<<dRt_wq1j[i]<<"   "<<dRt_wq2j[i]<<"   "<<matchingt[i]<<"   "<<matchingw[i]<<"   "<<matchingq[i]<<endl;

                //if(!(matchingtmp[i]==-1||matchingtmp[i]==-2)) cout<<i<<" matched no t  "<<dRtj[i]<<"   "<<dRantitj[i]<<"  "<<dRt1_bj[i]<<"   "<<dRt1_wq1j[i]<<"   "<<dRt1_wq2j[i]<<"  "<<dRt2_bj[i]<<"   "<<dRt2_wq1j[i]<<"   "<<dRt2_wq2j[i]<<"   "<<matchingt[i]<<"   "<<matchingw[i]<<"   "<<matchingq[i]<<"   "<<matchingg[i]<<endl;
                
                if(matchingt[i]<0&&matchingw[i]<0&&matchingg[i]<0&&matchingz[i]<0&&matchingq[i]<0&&matchingr[i]<0) unmatching[i]=1;
            }
            if(gentop_pt<0&&genantitop_pt<0){
                for(int j=0;j<5;j++){
                    if (dRwj[j][i]<rmatching){
                        if(i==0) dRw_q1j[j][i]=deltaR(genw_q1_eta[j],genw_q1_phi[j],eta_max,phi_max);
                        if(i==1&&Nj8!=2) dRw_q1j[j][i]=deltaR(genw_q1_eta[j],genw_q1_phi[j],eta_mid,phi_mid);
                        if(i==1&&Nj8==2) dRw_q1j[j][i]=deltaR(genw_q1_eta[j],genw_q1_phi[j],eta_min,phi_min);
                        if(i==2) dRw_q1j[j][i]=deltaR(genw_q1_eta[j],genw_q1_phi[j],eta_min,phi_min);
                        if(i==0) dRw_q2j[j][i]=deltaR(genw_q2_eta[j],genw_q2_phi[j],eta_max,phi_max);
                        if(i==1&&Nj8!=2) dRw_q2j[j][i]=deltaR(genw_q2_eta[j],genw_q2_phi[j],eta_mid,phi_mid);
                        if(i==1&&Nj8==2) dRw_q2j[j][i]=deltaR(genw_q2_eta[j],genw_q2_phi[j],eta_min,phi_min);
                        if(i==2) dRw_q2j[j][i]=deltaR(genw_q2_eta[j],genw_q2_phi[j],eta_min,phi_min);
                        if(dRw_q1j[j][i]<rbmatching&&dRw_q2j[j][i]<rbmatching&&taggenwl[j]==4) {
                            if(matchingw[i]!=-99) matchingw[i]+=1;
                            if(matchingw[i]==-99) matchingw[i]=1;
                            for(int j=0;j<10;j++){
                                if(dRgj[j][i]<rbmatching&&ptgengl[j]>0.1*genantit_w_pt) matchingw[i]=100;
                                if(dRgj[j][i]<rbmatching&&ptgengl[j]>0.15*genantit_w_pt) matchingw[i]=101;
                            }
                            for(int j=0;j<5;j++){
                                if((dRq1j[j][i]<rbmatching&&ptgenq1l[j]<0.1*genantit_w_pt)||(dRq2j[j][i]<rbmatching&&ptgenq2l[j]<0.1*genantit_w_pt)||(dRq3j[j][i]<rbmatching&&ptgenq3l[j]<0.1*genantit_w_pt)||(dRq4j[j][i]<rbmatching&&ptgenq4l[j]<0.1*genantit_w_pt)||(dRq5j[j][i]<rbmatching&&ptgenq5l[j]<0.1*genantit_w_pt)) matchingw[i]=102;
                                if((dRq1j[j][i]<rbmatching&&ptgenq1l[j]>0.15*genantit_w_pt)||(dRq2j[j][i]<rbmatching&&ptgenq2l[j]>0.15*genantit_w_pt)||(dRq3j[j][i]<rbmatching&&ptgenq3l[j]>0.15*genantit_w_pt)||(dRq4j[j][i]<rbmatching&&ptgenq4l[j]>0.15*genantit_w_pt)||(dRq5j[j][i]<rbmatching&&ptgenq5l[j]>0.15*genantit_w_pt)) matchingw[i]=103;
                            }

                        }
                        if(((dRw_q1j[j][i]<rbmatching&&dRw_q2j[j][i]>rbmatching)||(dRw_q1j[j][i]>rbmatching&&dRw_q2j[j][i]<rbmatching))&&taggenwl[j]==4) {
                            if(matchingq[i]!=-99) matchingq[i]+=1;
                            if(matchingq[i]==-99) matchingq[i]=1;
                        }
                    }
                }
                if(IsData==10){
                    double dRw2w3=99;
                    /*if(status_2==1) dRqq_W1=deltaR(gen_ele_eta_2,gen_ele_phi_2,gen_nele_eta_2,gen_nele_phi_2);
                    if(status_2==2) dRqq_W1=deltaR(gen_mu_eta_2,gen_mu_phi_2,gen_nmu_eta_2,gen_nmu_phi_2);
                    if(status_2==3) dRqq_W1=deltaR(gen_tau_eta_2,gen_tau_phi_2,gen_ntau_eta_2,gen_ntau_phi_2);
                    if(status_2==4) dRqq_W1=deltaR(etaq_2[0],phiq_2[0],etaq_2[1],phiq_2[1]);
                    if(status_3==1) dRqq_W2=deltaR(gen_ele_eta_3,gen_ele_phi_3,gen_nele_eta_3,gen_nele_phi_3);
                    if(status_3==2) dRqq_W2=deltaR(gen_mu_eta_3,gen_mu_phi_3,gen_nmu_eta_3,gen_nmu_phi_3);
                    if(status_3==3) dRqq_W2=deltaR(gen_tau_eta_3,gen_tau_phi_3,gen_ntau_eta_3,gen_ntau_phi_3);
                    if(status_3==4) dRqq_W2=deltaR(etaq_3[0],phiq_3[0],etaq_3[1],phiq_3[1]);
                    mindRqq=TMath::Min(dRqq_W1,dRqq_W2);
                    maxdRqq=TMath::Max(dRqq_W1,dRqq_W2);*/

                    if(status_2==4&&status_3==4) {
                        dRw2w3=deltaR(etaGenVhad_2,phiGenVhad_2,etaGenVhad_3,phiGenVhad_3);
                        if(i==0) dRw2j[i]=deltaR(etaGenVhad_2,phiGenVhad_2,eta_max,phi_max);
                        if(i==1&&Nj8!=2) dRw2j[i]=deltaR(etaGenVhad_2,phiGenVhad_2,eta_mid,phi_mid);
                        if(i==1&&Nj8==2) dRw2j[i]=deltaR(etaGenVhad_2,phiGenVhad_2,eta_min,phi_min);
                        if(i==2) dRw2j[i]=deltaR(etaGenVhad_2,phiGenVhad_2,eta_min,phi_min);
                        if(i==0) dRw3j[i]=deltaR(etaGenVhad_3,phiGenVhad_3,eta_max,phi_max);
                        if(i==1&&Nj8!=2) dRw3j[i]=deltaR(etaGenVhad_3,phiGenVhad_3,eta_mid,phi_mid);
                        if(i==1&&Nj8==2) dRw3j[i]=deltaR(etaGenVhad_3,phiGenVhad_3,eta_min,phi_min);
                        if(i==2) dRw3j[i]=deltaR(etaGenVhad_3,phiGenVhad_3,eta_min,phi_min);
                        dRWW_R=deltaR(etaGenVhad_2,phiGenVhad_2,etaGenVhad_3,phiGenVhad_3);
                    }
                    if(status_2==4&&status_3!=4) {
                        dRw2w3=deltaR(etaGenVhad_2,phiGenVhad_2,etaGenV_3,phiGenV_3);
                        if(i==0) dRw2j[i]=deltaR(etaGenVhad_2,phiGenVhad_2,eta_max,phi_max);
                        if(i==1&&Nj8!=2) dRw2j[i]=deltaR(etaGenVhad_2,phiGenVhad_2,eta_mid,phi_mid);
                        if(i==1&&Nj8==2) dRw2j[i]=deltaR(etaGenVhad_2,phiGenVhad_2,eta_min,phi_min);
                        if(i==2) dRw2j[i]=deltaR(etaGenVhad_2,phiGenVhad_2,eta_min,phi_min);
                        if(i==0) dRw3j[i]=deltaR(etaGenV_3,phiGenV_3,eta_max,phi_max);
                        if(i==1&&Nj8!=2) dRw3j[i]=deltaR(etaGenV_3,phiGenV_3,eta_mid,phi_mid);
                        if(i==1&&Nj8==2) dRw3j[i]=deltaR(etaGenV_3,phiGenV_3,eta_min,phi_min);
                        if(i==2) dRw3j[i]=deltaR(etaGenV_3,phiGenV_3,eta_min,phi_min);
                        dRWW_R=deltaR(etaGenVhad_2,phiGenVhad_2,etaGenV_3,phiGenV_3);
                    }
                    if(status_2!=4&&status_3==4) {
                        dRw2w3=deltaR(etaGenV_2,phiGenV_2,etaGenVhad_3,phiGenVhad_3);
                        if(i==0) dRw2j[i]=deltaR(etaGenV_2,phiGenV_2,eta_max,phi_max);
                        if(i==1&&Nj8!=2) dRw2j[i]=deltaR(etaGenV_2,phiGenV_2,eta_mid,phi_mid);
                        if(i==1&&Nj8==2) dRw2j[i]=deltaR(etaGenV_2,phiGenV_2,eta_min,phi_min);
                        if(i==2) dRw2j[i]=deltaR(etaGenV_2,phiGenV_2,eta_min,phi_min);
                        if(i==0) dRw3j[i]=deltaR(etaGenVhad_3,phiGenVhad_3,eta_max,phi_max);
                        if(i==1&&Nj8!=2) dRw3j[i]=deltaR(etaGenVhad_3,phiGenVhad_3,eta_mid,phi_mid);
                        if(i==1&&Nj8==2) dRw3j[i]=deltaR(etaGenVhad_3,phiGenVhad_3,eta_min,phi_min);
                        if(i==2) dRw3j[i]=deltaR(etaGenVhad_3,phiGenVhad_3,eta_min,phi_min);
                        dRWW_R=deltaR(etaGenV_2,phiGenV_2,etaGenVhad_3,phiGenVhad_3);
                    }

                    //if(dRrj[i]<rmatching&&dRw2j[i]<rmatching&&dRw3j[i]<rmatching&&dRw2w3<1.6){
                    if(dRrj[i]<rmatching&&dRw2j[i]<rmatching&&dRw3j[i]<rmatching){
                    //if(dRrj[i]<rmatching){
                        if(i==0) dRr_q1j[i]=deltaR(etaq_2[0],phiq_2[0],eta_max,phi_max);
                        if(i==1&&Nj8!=2) dRr_q1j[i]=deltaR(etaq_2[0],phiq_2[0],eta_mid,phi_mid);
                        if(i==1&&Nj8==2) dRr_q1j[i]=deltaR(etaq_2[0],phiq_2[0],eta_min,phi_min);
                        if(i==2) dRr_q1j[i]=deltaR(etaq_2[0],phiq_2[0],eta_min,phi_min);
                        if(i==0) dRr_q2j[i]=deltaR(etaq_2[1],phiq_2[1],eta_max,phi_max);
                        if(i==1&&Nj8!=2) dRr_q2j[i]=deltaR(etaq_2[1],phiq_2[1],eta_mid,phi_mid);
                        if(i==1&&Nj8==2) dRr_q2j[i]=deltaR(etaq_2[1],phiq_2[1],eta_min,phi_min);
                        if(i==2) dRr_q2j[i]=deltaR(etaq_2[1],phiq_2[1],eta_min,phi_min);
                        if(i==0) dRr_q3j[i]=deltaR(etaq_3[0],phiq_3[0],eta_max,phi_max);
                        if(i==1&&Nj8!=2) dRr_q3j[i]=deltaR(etaq_3[0],phiq_3[0],eta_mid,phi_mid);
                        if(i==1&&Nj8==2) dRr_q3j[i]=deltaR(etaq_3[0],phiq_3[0],eta_min,phi_min);
                        if(i==2) dRr_q3j[i]=deltaR(etaq_3[0],phiq_3[0],eta_min,phi_min);
                        if(i==0) dRr_q4j[i]=deltaR(etaq_3[1],phiq_3[1],eta_max,phi_max);
                        if(i==1&&Nj8!=2) dRr_q4j[i]=deltaR(etaq_3[1],phiq_3[1],eta_mid,phi_mid);
                        if(i==1&&Nj8==2) dRr_q4j[i]=deltaR(etaq_3[1],phiq_3[1],eta_min,phi_min);
                        if(i==2) dRr_q4j[i]=deltaR(etaq_3[1],phiq_3[1],eta_min,phi_min);
                        double dRr_qj[4]={99,99,99,99};
                        dRr_qj[0]=dRr_q1j[i];dRr_qj[1]=dRr_q2j[i];dRr_qj[2]=dRr_q3j[i];dRr_qj[3]=dRr_q4j[i];
                        Int_t *idRr_qj=new Int_t[4];
                        TMath::Sort(4,dRr_qj,idRr_qj,0);
                        //cout<<dRr_qj[idRr_qj[0]]<<"  "<<dRr_qj[idRr_qj[1]]<<"  "<<dRr_qj[idRr_qj[2]]<<"  "<<dRr_qj[idRr_qj[3]]<<endl;
                        if(dRr_qj[idRr_qj[3]]<rbmatching) {
                        matchingr[i]=4;
                        for(int j=0;j<10;j++){
                            if(dRgj[j][i]<rbmatching&&ptgengl[j]>0.1*gen_rad_pt) matchingr[i]=5;
                            if(dRgj[j][i]<rbmatching&&ptgengl[j]>0.15*gen_rad_pt) matchingr[i]=6;
                        }
                        for(int j=0;j<5;j++){
                            if((dRq1j[j][i]<rbmatching&&ptgenq1l[j]<0.1*gen_rad_pt)||(dRq2j[j][i]<rbmatching&&ptgenq2l[j]<0.1*gen_rad_pt)||(dRq3j[j][i]<rbmatching&&ptgenq3l[j]<0.1*gen_rad_pt)||(dRq4j[j][i]<rbmatching&&ptgenq4l[j]<0.1*gen_rad_pt)||(dRq5j[j][i]<rbmatching&&ptgenq5l[j]<0.1*gen_rad_pt)) matchingr[i]=7;
                            if((dRq1j[j][i]<rbmatching&&ptgenq1l[j]>0.15*gen_rad_pt)||(dRq2j[j][i]<rbmatching&&ptgenq2l[j]>0.15*gen_rad_pt)||(dRq3j[j][i]<rbmatching&&ptgenq3l[j]>0.15*gen_rad_pt)||(dRq4j[j][i]<rbmatching&&ptgenq4l[j]>0.15*gen_rad_pt)||(dRq5j[j][i]<rbmatching&&ptgenq5l[j]>0.15*gen_rad_pt)) matchingr[i]=8;
                        }
                        }
                        if(dRr_qj[idRr_qj[2]]<rbmatching&&dRr_qj[idRr_qj[3]]>rbmatching) matchingr[i]=3;
                        if(dRr_qj[idRr_qj[1]]<rbmatching&&dRr_qj[idRr_qj[2]]>rbmatching) matchingr[i]=2;
                        if(dRr_qj[idRr_qj[0]]<rbmatching&&dRr_qj[idRr_qj[1]]>rbmatching){
                            if(status_2==4&&status_3==4) matchingr[i]=-1;
                            if(!(status_2==4&&status_3==4)) matchingr[i]=1;
                        }
                    }
                }
                for(int j=0;j<5;j++){
                    if (dRzj[j][i]<rmatching){
                            if(matchingz[i]!=-99) matchingz[i]+=1;
                            if(matchingz[i]==-99) matchingz[i]=1;
                    }
                }
                for(int j=0;j<10;j++){
                    if (dRgj[j][i]<rmatching){
                        if(matchingg[i]!=-99) matchingg[i]+=1;
                        if(matchingg[i]==-99) matchingg[i]=1;
                    }
                }
                for(int j=0;j<5;j++){
                    if ((dRq1j[j][i]<rmatching)||(dRq2j[j][i]<rmatching)||(dRq3j[j][i]<rmatching)||(dRq4j[j][i]<rmatching)||(dRq5j[j][i]<rmatching)){
                        if(matchingq[i]!=-99) matchingq[i]+=1;
                        if(matchingq[i]==-99) matchingq[i]=1;
                    }
                }
            }
            //cout<<matchingr[i]<<" cc  "<<matchingw[i]<<"   "<<matchingq[i]<<"   "<<matchingg[i]<<endl;
            if(matchingt[i]<0&&matchingw[i]<0&&matchingg[i]<0&&matchingz[i]<0&&matchingq[i]<0&&matchingr[i]<0) unmatching[i]=1;
        }
        
	double dnncorr=1;//useless
	double dnncorr_all=1;        
	double dnncorr_nosft60=1;
        double dnncorr_nosfq60=1;
        double dnncorr_nosfw60=1;
        double dnncorr_nosft100=1;
        double dnncorr_nosfq100=1;
        double dnncorr_nosfw100=1;
        double dnncorr_nosft150=1;
        double dnncorr_nosfq150=1;
        double dnncorr_nosft200=1;
        double dnncorr_nosfq200=1;
	
	double dnn_corrt60[2]={1,1},dnn_corrq60[2]={1,1},dnn_corrw60[2]={1,1},dnn_corrt100[2]={1,1},dnn_corrq100[2]={1,1},dnn_corrw100[2]={1,1},dnn_corrt150[2]={1,1},dnn_corrq150[2]={1,1},dnn_corrt200[2]={1,1},dnn_corrq200[2]={1,1};
	
       if(IsData!=0){
            Double_t jet_mass[2]={Mj_max,Mj_min};
	    for(int inj8=0;inj8<Nj8&&inj8<3&&inj8<2;inj8++){
            	int ixd1=0;
	    	if(inj8==1&&Nj8==3) jet_mass[inj8]=Mj_mid;
            	double dnntagd1=dnndecw_max;
		if(inj8==1&&Nj8==2) dnntagd1=dnndecw_min;
		if(inj8==1&&Nj8==3) dnntagd1=dnndecw_mid;
            	if(jet_mass[inj8]>=100) dnntagd1=(dnndecraww_max+dnndecrawh4q_max)/(dnndecrawqcd_max+dnndecrawh4q_max+dnndecraww_max);
            	if(jet_mass[inj8]>=100&&inj8==1&&Nj8==2) dnntagd1=(dnndecraww_min+dnndecrawh4q_min)/(dnndecrawqcd_min+dnndecrawh4q_min+dnndecraww_min);
            	if(jet_mass[inj8]>=100&&inj8==1&&Nj8==3) dnntagd1=(dnndecraww_mid+dnndecrawh4q_mid)/(dnndecrawqcd_mid+dnndecrawh4q_mid+dnndecraww_mid);
            	for(ixd1=0;ixd1<4;ixd1++){
            	    if(jet_mass[inj8]<massdnn[0]) {ixd1=-1;break;}
            	    if(jet_mass[inj8]>=massdnn[ixd1]&&jet_mass[inj8]<massdnn[ixd1+1]) break;
            	}
            	bool tmatching1=(matchingt[inj8]>0)||(matchingr[inj8]==3 &&( (status_2==4&&status_3==4)))||(matchingr[inj8]>=4);
		bool wmatching1=(matchingr[inj8]<=0 &&matchingt[inj8]<=0 && matchingw[inj8]>0)||(matchingr[inj8]==2 &&(! (status_2==4&&status_3==4)));
		//bool wmatching1=(matchingr[inj8]<=0 &&matchingt[inj8]<=0 && matchingw[inj8]>0);
		bool qmatching1=matchingr[inj8]<=0 &&matchingt[inj8]<=0 && matchingw[inj8]<=0 &&(matchingg[inj8]>0 || matchingq[inj8]>0);
                if(dnntagd1>=0 && dnntagd1<=1){
		if(dnntagd1==1)dnntagd1=0.99;
            	if(ixd1==0){
			if(tmatching1) dnn_corrt60[inj8]=h_t60->GetBinContent(h_t60->FindBin(dnntagd1));
			if(wmatching1) dnn_corrw60[inj8]=h_w60->GetBinContent(h_w60->FindBin(dnntagd1));
			if(qmatching1) dnn_corrq60[inj8]=h_q60->GetBinContent(h_q60->FindBin(dnntagd1));
		}
            	if(ixd1==1){
			if(tmatching1) dnn_corrt100[inj8]=h_t100->GetBinContent(h_t100->FindBin(dnntagd1));
			if(wmatching1) dnn_corrw100[inj8]=h_w100->GetBinContent(h_w100->FindBin(dnntagd1));
			if(qmatching1) dnn_corrq100[inj8]=h_q100->GetBinContent(h_q100->FindBin(dnntagd1));
		}
            	if(ixd1==2){
			if(wmatching1&&IsData==10) dnn_corrw100[inj8]=h_w100->GetBinContent(h_w100->FindBin(dnntagd1));
			if(tmatching1) dnn_corrt150[inj8]=h_t150->GetBinContent(h_t150->FindBin(dnntagd1));
			if(qmatching1) dnn_corrq150[inj8]=h_q150->GetBinContent(h_q150->FindBin(dnntagd1));
		}
            	if(ixd1==3){
			if(wmatching1&&IsData==10) dnn_corrw100[inj8]=h_w100->GetBinContent(h_w100->FindBin(dnntagd1));
			if(tmatching1) dnn_corrt200[inj8]=h_t200->GetBinContent(h_t200->FindBin(dnntagd1));
			if(qmatching1) dnn_corrq200[inj8]=h_q200->GetBinContent(h_q200->FindBin(dnntagd1));
		}}
		dnncorr_all*=dnn_corrt60[inj8]*dnn_corrw60[inj8]*dnn_corrq60[inj8]*dnn_corrt100[inj8]*dnn_corrq100[inj8]*dnn_corrw100[inj8]*dnn_corrt150[inj8]*dnn_corrq150[inj8]*dnn_corrt200[inj8]*dnn_corrq200[inj8];
		dnncorr_nosft60*=dnn_corrt60[inj8]*dnn_corrw60[inj8]*dnn_corrq60[inj8]*dnn_corrt100[inj8]*dnn_corrq100[inj8]*dnn_corrw100[inj8]*dnn_corrt150[inj8]*dnn_corrq150[inj8]*dnn_corrt200[inj8]*dnn_corrq200[inj8]/dnn_corrt60[inj8];
		dnncorr_nosfq60*=dnn_corrt60[inj8]*dnn_corrw60[inj8]*dnn_corrq60[inj8]*dnn_corrt100[inj8]*dnn_corrq100[inj8]*dnn_corrw100[inj8]*dnn_corrt150[inj8]*dnn_corrq150[inj8]*dnn_corrt200[inj8]*dnn_corrq200[inj8]/dnn_corrq60[inj8];
		dnncorr_nosfw60*=dnn_corrt60[inj8]*dnn_corrw60[inj8]*dnn_corrq60[inj8]*dnn_corrt100[inj8]*dnn_corrq100[inj8]*dnn_corrw100[inj8]*dnn_corrt150[inj8]*dnn_corrq150[inj8]*dnn_corrt200[inj8]*dnn_corrq200[inj8]/dnn_corrw60[inj8];
		dnncorr_nosft100*=dnn_corrt60[inj8]*dnn_corrw60[inj8]*dnn_corrq60[inj8]*dnn_corrt100[inj8]*dnn_corrq100[inj8]*dnn_corrw100[inj8]*dnn_corrt150[inj8]*dnn_corrq150[inj8]*dnn_corrt200[inj8]*dnn_corrq200[inj8]/dnn_corrt100[inj8];
		dnncorr_nosfq100*=dnn_corrt60[inj8]*dnn_corrw60[inj8]*dnn_corrq60[inj8]*dnn_corrt100[inj8]*dnn_corrq100[inj8]*dnn_corrw100[inj8]*dnn_corrt150[inj8]*dnn_corrq150[inj8]*dnn_corrt200[inj8]*dnn_corrq200[inj8]/dnn_corrq100[inj8];
		dnncorr_nosfw100*=dnn_corrt60[inj8]*dnn_corrw60[inj8]*dnn_corrq60[inj8]*dnn_corrt100[inj8]*dnn_corrq100[inj8]*dnn_corrw100[inj8]*dnn_corrt150[inj8]*dnn_corrq150[inj8]*dnn_corrt200[inj8]*dnn_corrq200[inj8]/dnn_corrw100[inj8];
		dnncorr_nosft150*=dnn_corrt60[inj8]*dnn_corrw60[inj8]*dnn_corrq60[inj8]*dnn_corrt100[inj8]*dnn_corrq100[inj8]*dnn_corrw100[inj8]*dnn_corrt150[inj8]*dnn_corrq150[inj8]*dnn_corrt200[inj8]*dnn_corrq200[inj8]/dnn_corrt150[inj8];
		dnncorr_nosfq150*=dnn_corrt60[inj8]*dnn_corrw60[inj8]*dnn_corrq60[inj8]*dnn_corrt100[inj8]*dnn_corrq100[inj8]*dnn_corrw100[inj8]*dnn_corrt150[inj8]*dnn_corrq150[inj8]*dnn_corrt200[inj8]*dnn_corrq200[inj8]/dnn_corrq150[inj8];
		dnncorr_nosft200*=dnn_corrt60[inj8]*dnn_corrw60[inj8]*dnn_corrq60[inj8]*dnn_corrt100[inj8]*dnn_corrq100[inj8]*dnn_corrw100[inj8]*dnn_corrt150[inj8]*dnn_corrq150[inj8]*dnn_corrt200[inj8]*dnn_corrq200[inj8]/dnn_corrt200[inj8];
		dnncorr_nosfq200*=dnn_corrt60[inj8]*dnn_corrw60[inj8]*dnn_corrq60[inj8]*dnn_corrt100[inj8]*dnn_corrq100[inj8]*dnn_corrw100[inj8]*dnn_corrt150[inj8]*dnn_corrq150[inj8]*dnn_corrt200[inj8]*dnn_corrq200[inj8]/dnn_corrq200[inj8];
		//if(inj8==0) cout<<" massdnn[ixd1]  "<<massdnn[ixd1]<<" Mj_max  "<<Mj_max<<" dnndecw_max "<<dnndecw_max<<" wh "<<(dnndecraww_max+dnndecrawh4q_max)/(dnndecrawqcd_max+dnndecrawh4q_max+dnndecraww_max)<<" dnntagd1 "<<dnntagd1<<" wmatching1  "<<wmatching1<<" tmatching1 "<<tmatching1<<" qmatching1 "<<qmatching1<<" sfall "<<dnncorr_all<<endl;
		//if(inj8==1) cout<<" massdnn[ixd1]  "<<massdnn[ixd1]<<" Mj_max  "<<Mj_max<<" dnndecw_max "<<dnndecw_max<<" wh "<<(dnndecraww_max+dnndecrawh4q_max)/(dnndecrawqcd_max+dnndecrawh4q_max+dnndecraww_max)<<" dnntagd1 "<<dnntagd1<<" wmatching1  "<<wmatching1<<" tmatching1 "<<tmatching1<<" qmatching1 "<<qmatching1<<" sfall "<<dnncorr_all<<endl;
	     }
	}
	//cout<<" inj8  "<<Nj8<<" dnncorr_all "<<dnncorr_all<<endl;
	//	if(jentry>=200) break;
        weight_deep_matchedcorr1=weight_deep*dnncorr;
        weight_deep_matchedcorr2=weight_deep*dnncorr*dnncorr;
        weight_deep_tig_matchedcorr1=weight_tig_deep*dnncorr;
        weight_deep_tig_matchedcorr2=weight_tig_deep*dnncorr*dnncorr;
        weight_deep_los_matchedcorr1=weight_los_deep*dnncorr;
        weight_deep_los_matchedcorr2=weight_los_deep*dnncorr*dnncorr;
	//cout<<dnncorr_all<<" dnncorr_all "<<dnncorr_nosft60<<"  "<<dnncorr_nosfq60<<"  "<<dnncorr_nosfw60<<"  "<<dnncorr_nosft100<<"  "<<dnncorr_nosfq100<<"  "<<dnncorr_nosfw100<<"  "<<dnncorr_nosft150<<"  "<<dnncorr_nosfq150<<"  "<<dnncorr_nosft200<<"  "<<dnncorr_nosfq200<<"  "<<endl;
        weight_deep_lumi_corr1=weight_deep_lumi*dnncorr_all;
	weight_deep_lumi_nosft60=weight_deep_lumi*dnncorr_nosft60;
	weight_deep_lumi_nosfq60=weight_deep_lumi*dnncorr_nosfq60;
	weight_deep_lumi_nosfw60=weight_deep_lumi*dnncorr_nosfw60;
	weight_deep_lumi_nosft100=weight_deep_lumi*dnncorr_nosft100;
	weight_deep_lumi_nosfq100=weight_deep_lumi*dnncorr_nosfq100;
	weight_deep_lumi_nosfw100=weight_deep_lumi*dnncorr_nosfw100;
	weight_deep_lumi_nosft150=weight_deep_lumi*dnncorr_nosft150;
	weight_deep_lumi_nosfq150=weight_deep_lumi*dnncorr_nosfq150;
	weight_deep_lumi_nosft200=weight_deep_lumi*dnncorr_nosft200;
	weight_deep_lumi_nosfq200=weight_deep_lumi*dnncorr_nosfq200;
        weight_deep_lumi_corr2=weight_deep_lumi*dnncorr_all*dnncorr_all;

        double etaq1=-99,etaq2=-99,etaq3=-99,etaq4=-99,phiq1=-99,phiq2=-99,phiq3=-99,phiq4=-99,ptq1=-99,ptq2=-99,ptq3=-99,ptq4=-99;
        double q4v_R[4][3]={-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99};
        for(int i=0;i<4;i++){
            nqqcr[i]=0;nqqct[i]=0;nqqcw[i]=0;
            for(int j=0;j<3;j++)    {matchnewr[i][j]=-99;matchnewt[i][j]=-99;matchneww[i][j]=-99;}
        }
        double ignoreprong=0.15;
        if(IsData==10){
            if(status_2==1){etaq1=gen_ele_eta_2;phiq1=gen_ele_phi_2;ptq1=gen_ele_pt_2;etaq2=gen_nele_eta_2;phiq2=gen_nele_phi_2;ptq2=gen_nele_phi_2;}
            if(status_2==2){etaq1=gen_mu_eta_2;phiq1=gen_mu_phi_2;etaq2=gen_nmu_eta_2;phiq2=gen_nmu_phi_2;}
            if(status_2==3){etaq1=gen_tau_eta_2;phiq1=gen_tau_phi_2;etaq2=gen_ntau_eta_2;phiq2=gen_ntau_phi_2;}
            if(status_2==4){etaq1=etaq_2[0];phiq1=phiq_2[0];etaq2=etaq_2[1];phiq2=phiq_2[1];}
            if(status_3==1){etaq3=gen_ele_eta_3;phiq3=gen_ele_phi_3;etaq4=gen_nele_eta_3;phiq4=gen_nele_phi_3;}
            if(status_3==2){etaq3=gen_mu_eta_3;phiq3=gen_mu_phi_3;etaq4=gen_nmu_eta_3;phiq4=gen_nmu_phi_3;}
            if(status_3==3){etaq3=gen_tau_eta_3;phiq3=gen_tau_phi_3;etaq4=gen_ntau_eta_3;phiq4=gen_ntau_phi_3;}
            if(status_3==4){etaq3=etaq_3[0];phiq3=phiq_3[0];etaq4=etaq_3[1];phiq4=phiq_3[1];}
            
            if(status_2==1){q4v_R[0][0]=gen_ele_pt_2;q4v_R[0][1]=gen_ele_eta_2;q4v_R[0][2]=gen_ele_phi_2;q4v_R[1][0]=gen_nele_pt_2;q4v_R[1][1]=gen_nele_eta_2;q4v_R[1][2]=gen_nele_phi_2;}
            if(status_2==2){q4v_R[0][0]=gen_mu_pt_2;q4v_R[0][1]=gen_mu_eta_2;q4v_R[0][2]=gen_mu_phi_2;q4v_R[1][0]=gen_nmu_pt_2;q4v_R[1][1]=gen_nmu_eta_2;q4v_R[1][2]=gen_nmu_phi_2;}
            if(status_2==3){q4v_R[0][0]=gen_tau_pt_2;q4v_R[0][1]=gen_tau_eta_2;q4v_R[0][2]=gen_tau_phi_2;q4v_R[1][0]=gen_ntau_pt_2;q4v_R[1][1]=gen_ntau_eta_2;q4v_R[1][2]=gen_ntau_phi_2;}
            if(status_2==4){q4v_R[0][0]=ptq_2[0];q4v_R[0][1]=etaq_2[0];q4v_R[0][2]=phiq_2[0];q4v_R[1][0]=ptq_2[1];q4v_R[1][1]=etaq_2[1];q4v_R[1][2]=phiq_2[1];}
            if(status_3==1){q4v_R[2][0]=gen_ele_pt_3;q4v_R[2][1]=gen_ele_eta_3;q4v_R[2][2]=gen_ele_phi_3;q4v_R[3][0]=gen_nele_pt_3;q4v_R[3][1]=gen_nele_eta_3;q4v_R[3][2]=gen_nele_phi_3;}
            if(status_3==2){q4v_R[2][0]=gen_mu_pt_3;q4v_R[2][1]=gen_mu_eta_3;q4v_R[2][2]=gen_mu_phi_3;q4v_R[3][0]=gen_nmu_pt_3;q4v_R[3][1]=gen_nmu_eta_3;q4v_R[3][2]=gen_nmu_phi_3;}
            if(status_3==3){q4v_R[2][0]=gen_tau_pt_3;q4v_R[2][1]=gen_tau_eta_3;q4v_R[2][2]=gen_tau_phi_3;q4v_R[3][0]=gen_ntau_pt_3;q4v_R[3][1]=gen_ntau_eta_3;q4v_R[3][2]=gen_ntau_phi_3;}
            if(status_3==4){q4v_R[2][0]=ptq_3[0];q4v_R[2][1]=etaq_3[0];q4v_R[2][2]=phiq_3[0];q4v_R[3][0]=ptq_3[1];q4v_R[3][1]=etaq_3[1];q4v_R[3][2]=phiq_3[1];}

            dRqq_W1=deltaR(q4v_R[0][1],q4v_R[0][2],q4v_R[1][1],q4v_R[1][2]);
            dRqq_W2=deltaR(q4v_R[2][1],q4v_R[2][2],q4v_R[3][1],q4v_R[3][2]);
            mindRqq=TMath::Min(dRqq_W1,dRqq_W2);
            maxdRqq=TMath::Max(dRqq_W1,dRqq_W2);
            
             double dRqq_all[6]={deltaR(q4v_R[0][1],q4v_R[0][2],q4v_R[1][1],q4v_R[1][2]),deltaR(q4v_R[0][1],q4v_R[0][2],q4v_R[2][1],q4v_R[2][2]),deltaR(q4v_R[0][1],q4v_R[0][2],q4v_R[3][1],q4v_R[3][2]),deltaR(q4v_R[1][1],q4v_R[1][2],q4v_R[2][1],q4v_R[2][2]),deltaR(q4v_R[1][1],q4v_R[1][2],q4v_R[3][1],q4v_R[3][2]),deltaR(q4v_R[2][1],q4v_R[2][2],q4v_R[3][1],q4v_R[3][2])};
             Int_t *index_all=new Int_t[6];
             TMath::Sort(6,dRqq_all,index_all,0);
             mindRqq_all=dRqq_all[index_all[0]];
             maxdRqq_all=dRqq_all[index_all[5]];

            double extraqgr[5][3]={0};
            int insideqgr=0;
            for(int j=0;j<10;j++){
                if(dRgj[j][0]<rbmatching&&ptgengl[j]>0.15*gen_rad_pt&&ptgengl[j]>extraqgr[0][0]) {extraqgr[0][0]=ptgengl[j];extraqgr[0][1]=etagengl[j];extraqgr[0][2]=phigengl[j];}
            }
            for(int j=0;j<5;j++){
                if((dRq1j[j][0]<rbmatching&&ptgenq1l[j]>0.15*gen_rad_pt&&ptgenq1l[j]>extraqgr[0][0])){extraqgr[0][0]=ptgenq1l[j];extraqgr[0][1]=etagenq1l[j];extraqgr[0][2]=phigenq1l[j];}
                if((dRq2j[j][0]<rbmatching&&ptgenq2l[j]>0.15*gen_rad_pt&&ptgenq2l[j]>extraqgr[0][0])){extraqgr[0][0]=ptgenq2l[j];extraqgr[0][1]=etagenq2l[j];extraqgr[0][2]=phigenq2l[j];}
                if((dRq3j[j][0]<rbmatching&&ptgenq3l[j]>0.15*gen_rad_pt&&ptgenq3l[j]>extraqgr[0][0])){extraqgr[0][0]=ptgenq3l[j];extraqgr[0][1]=etagenq3l[j];extraqgr[0][2]=phigenq3l[j];}
                if((dRq4j[j][0]<rbmatching&&ptgenq4l[j]>0.15*gen_rad_pt&&ptgenq4l[j]>extraqgr[0][0])){extraqgr[0][0]=ptgenq4l[j];extraqgr[0][1]=etagenq4l[j];extraqgr[0][2]=phigenq4l[j];}
                if((dRq5j[j][0]<rbmatching&&ptgenq5l[j]>0.15*gen_rad_pt&&ptgenq5l[j]>extraqgr[0][0])){extraqgr[0][0]=ptgenq5l[j];extraqgr[0][1]=etagenq5l[j];extraqgr[0][2]=phigenq5l[j];}
            }
            dRr_q1j[0]=deltaR(q4v_R[0][1],q4v_R[0][2],eta_max,phi_max);
            dRr_q2j[0]=deltaR(q4v_R[1][1],q4v_R[1][2],eta_max,phi_max);
            dRr_q3j[0]=deltaR(q4v_R[2][1],q4v_R[2][2],eta_max,phi_max);
            dRr_q4j[0]=deltaR(q4v_R[3][1],q4v_R[3][2],eta_max,phi_max);
            if(status_2==4&&status_3==4) {
                dRw2j[0]=deltaR(etaGenVhad_2,phiGenVhad_2,eta_max,phi_max);
                dRw3j[0]=deltaR(etaGenVhad_3,phiGenVhad_3,eta_max,phi_max);
            }
            if(status_2==4&&status_3!=4) {
                dRw2j[0]=deltaR(etaGenVhad_2,phiGenVhad_2,eta_max,phi_max);
                dRw3j[0]=deltaR(etaGenV_3,phiGenV_3,eta_max,phi_max);
            }
            if(status_2!=4&&status_3==4) {
                dRw2j[0]=deltaR(etaGenV_2,phiGenV_2,eta_max,phi_max);
                dRw3j[0]=deltaR(etaGenVhad_3,phiGenVhad_3,eta_max,phi_max);
            }
            int tmpinr=0;
            //if(dRrj[0]<rmatching){
            if(dRrj[0]<rmatching&&dRw2j[0]<rmatching&&dRw3j[0]<rmatching){
                if(extraqgr[0][0]>0) tmpinr++;
                //if(tmpinr==1){cout<<" checkgpt "<<extraqgr[0][0]<<endl;}
                if(dRr_q1j[0]<rbmatching){extraqgr[tmpinr][0]=q4v_R[0][0];extraqgr[tmpinr][1]=q4v_R[0][1];extraqgr[tmpinr][2]=q4v_R[0][2];tmpinr++;}
                if(dRr_q2j[0]<rbmatching&&status_2==4){extraqgr[tmpinr][0]=q4v_R[1][0];extraqgr[tmpinr][1]=q4v_R[1][1];extraqgr[tmpinr][2]=q4v_R[1][2];tmpinr++;}
                if(dRr_q3j[0]<rbmatching){extraqgr[tmpinr][0]=q4v_R[2][0];extraqgr[tmpinr][1]=q4v_R[2][1];extraqgr[tmpinr][2]=q4v_R[2][2];tmpinr++;}
                if(dRr_q4j[0]<rbmatching&&status_3==4){extraqgr[tmpinr][0]=q4v_R[3][0];extraqgr[tmpinr][1]=q4v_R[3][1];extraqgr[tmpinr][2]=q4v_R[3][2];tmpinr++;}
                /*for(int i=0;i<tmpinr;i++){
                    cout<<"check inr pt"<<i<<"    "<<extraqgr[i][0];}
                cout<<endl;*/
            }
            for(int k=0;k<4;k++) {matchnewr[k][0]=tmpinr;matchnewr[k][2]=tmpinr;
                if(extraqgr[0][0]>0)matchnewr[k][2]=tmpinr-1;
            }
            
            double PTtar=gen_rad_pt;
            const int maxgroupr=tmpinr;
            int pgroupr[4][2*maxgroupr],subgroupr[4][maxgroupr][maxgroupr],nsubgroupr[4][maxgroupr],flagw1[4][2*maxgroupr],flagalone[4][maxgroupr];
            double ptgroupr[4][maxgroupr];
            for(int i=0;i<4;i++){
                for(int j=0;j<maxgroupr;j++){
                    ptgroupr[i][j]=0;nsubgroupr[i][j]=-99;flagalone[i][j]=-99;
                    for(int k=0;k<maxgroupr;k++) subgroupr[i][j][k]=-99;}
                for(int j=0;j<2*maxgroupr;j++) {pgroupr[i][j]=-99;flagw1[i][j]=-99;}
            }
            double drqqin5[2*tmpinr];
            for(int i=0;i<2*tmpinr;i++){drqqin5[i]=-99;}
            int idrqqin5=0;
            for(int m=0;m<tmpinr;m++){
                for(int n=m+1;n<tmpinr;n++){
                    drqqin5[idrqqin5]=deltaR(extraqgr[m][1],extraqgr[m][2],extraqgr[n][1],extraqgr[n][2]);
                    //cout<<m<<"   "<<n<<"   "<<extraqgr[m][1]<<"   "<<extraqgr[m][2]<<"   "<<extraqgr[n][1]<<"   "<<extraqgr[n][2]<<"   "<<drqqin5[idrqqin5]<<endl;
                    idrqqin5++;
                }
            }
            int i1sub,i2sub;
            double drwp[4]={0.1,0.15,0.2,0.25};
            int end1,end2;
            for(int i=0;i<4;i++){
                i1sub=0;i2sub=0;
            for(int j=0;j<2*tmpinr;j++){
                if(0<drqqin5[j]&&drqqin5[j]<drwp[i]){
                    end1=endpoint(j,tmpinr,1);end2=endpoint(j,tmpinr,2);
                    //cout<<" check  "<<j<<"   "<<end1<<"   "<<end2<<"   "<<flagw1[i][end1]<<"   "<<flagw1[i][end2]<<" i1sub  "<<i1sub<<"   "<<i2sub<<"   "<<endl;
                    if(flagw1[i][end1]==-99&&flagw1[i][end2]==-99) {if(subgroupr[i][i1sub][0]!=-99){i1sub++;i2sub=0;}subgroupr[i][i1sub][i2sub]=end1;subgroupr[i][i1sub][i2sub+1]=end2;
                        flagw1[i][end1]=1;flagw1[i][end2]=1;i2sub+=2;
                        for(int k=0;k<2*tmpinr;k++){
                            if(0<drqqin5[k]&&drqqin5[k]<drwp[i]){
                                end1=endpoint(k,tmpinr,1);end2=endpoint(k,tmpinr,2);
                                if((flagw1[i][end1]==-99&&flagw1[i][end2]==1)||(flagw1[i][end1]==1&&flagw1[i][end2]==-99)){
                                    int endtmp=end1;if(flagw1[i][end1]==1) endtmp=end2;
                                    subgroupr[i][i1sub][i2sub]=endtmp;i2sub++;flagw1[i][endtmp]=1;
                                    for(int m=0;m<2*tmpinr;m++){
                                        if(0<drqqin5[m]&&drqqin5[m]<drwp[i]){
                                            end1=endpoint(m,tmpinr,1);end2=endpoint(m,tmpinr,2);
                                            if((flagw1[i][end1]==-99&&flagw1[i][end2]==1)||(flagw1[i][end1]==1&&flagw1[i][end2]==-99)){
                                                int endtmp=end1;if(flagw1[i][end1]==1) endtmp=end2;
                                                subgroupr[i][i1sub][i2sub]=endtmp;i2sub++;flagw1[i][endtmp]=1;
                                                for(int n=0;n<2*tmpinr;n++){
                                                    if(0<drqqin5[n]&&drqqin5[n]<drwp[i]){
                                                        end1=endpoint(n,tmpinr,1);end2=endpoint(n,tmpinr,2);
                                                        if((flagw1[i][end1]==-99&&flagw1[i][end2]==1)||(flagw1[i][end1]==1&&flagw1[i][end2]==-99)){
                                                            int endtmp=end1;if(flagw1[i][end1]==1) endtmp=end2;
                                                            subgroupr[i][i1sub][i2sub]=endtmp;i2sub++;flagw1[i][endtmp]=1;
                                                        }}}
                                            }}}
                                }
                            }
                        }
                    }
                }//end
            }
                int totalgroup=0;
                for(int j=0;j<tmpinr;j++){if(subgroupr[i][j][0]!=-99)totalgroup++;}
                const int ntotal=totalgroup;
                int flagw2[ntotal][tmpinr];
                for(int j=0;j<tmpinr;j++){
                    for(int k=0;k<ntotal;k++){
                        flagw2[k][j]=-99;
                    for(int m=0;m<tmpinr;m++){
                        if(subgroupr[i][k][m]==j){flagw2[k][j]=1;}
                    }}
                    int flagw3=-99;
                    for(int k=0;k<ntotal;k++){
                        if(flagw2[k][j]==1) flagw3=1;}
                    if(flagw3==-99){subgroupr[i][totalgroup][0]=j;totalgroup++;}
                     }
            }
            for(int i=0;i<4;i++){
                for(int j=0;j<tmpinr;j++){
                    int igroup=0;
                    double tmpsubg[5][2]={0};
                    for(int k=0;k<tmpinr;k++){
                        if(subgroupr[i][j][k]!=-99) {tmpsubg[igroup][0]=extraqgr[subgroupr[i][j][k]][0];tmpsubg[igroup][1]=extraqgr[subgroupr[i][j][k]][2];igroup++;}
                    }
                    ptgroupr[i][j]=SumPt(tmpsubg,igroup);
                    if(ptgroupr[i][j]/PTtar>ignoreprong) {if(matchnewr[i][1]!=-99) matchnewr[i][1]++;if(matchnewr[i][1]==-99) matchnewr[i][1]=1;}
                    /*cout<<"groupr"<<j;
                    for(int k=0;k<tmpinr;k++)cout<<"   "<<subgroupr[i][j][k];
                    cout<<"   "<<ptgroupr[i][j]<<"    "<<ptgroupr[i][j]/PTtar<<"    "<<matchnewr[i][1]<<endl;*/
            }
            }
        }
        
        if(IsData>0&&(gentop_pt>0||genantitop_pt>0)){
            
            double extraqgt[4][3]={0};
            if (dRtj[0]<rmatching) matchingtmp[0]=-1;
            if (matchingtmp[0]!=-1&&dRantitj[0]<rmatching) matchingtmp[0]=-2;
            if (matchingtmp[0]==-1&&dRantitj[0]<rmatching&&dRtj[0]<dRantitj[0]) matchingtmp[0]=-1;
            if (matchingtmp[0]==-1&&dRantitj[0]<rmatching&&dRtj[0]>dRantitj[0]) matchingtmp[0]=-2;
            //cout<<i<<" jet  "<<Nj8<<"   "<<matchingtmp[0]<<endl;
            if(matchingtmp[0]==-1||matchingtmp[0]==-2){
                if(matchingtmp[0]==-1){tmptpt=gentop_pt;tmpteta=gentop_eta;tmptphi=gentop_phi;tmpbpt=gent_b_pt;tmpbeta=gent_b_eta;tmpbphi=gent_b_phi;tmpwpt=gent_w_pt;tmpweta=gent_w_eta;tmpwphi=gent_w_phi;tmpq1pt=gent_w_q1_pt;tmpq1eta=gent_w_q1_eta;tmpq1phi=gent_w_q1_phi;tmpq2pt=gent_w_q2_pt;tmpq2eta=gent_w_q2_eta;tmpq2phi=gent_w_q2_phi;tmpbm=gent_b_mass;tmpt_wtag=gent_w_tag;}
                if(matchingtmp[0]==-2){tmptpt=genantitop_pt;tmpteta=genantitop_eta;tmptphi=genantitop_phi;tmpbpt=genantit_b_pt;tmpbeta=genantit_b_eta;tmpbphi=genantit_b_phi;tmpwpt=genantit_w_pt;tmpweta=genantit_w_eta;tmpwphi=genantit_w_phi;tmpq1pt=genantit_w_q1_pt;tmpq1eta=genantit_w_q1_eta;tmpq1phi=genantit_w_q1_phi;tmpq2pt=genantit_w_q2_pt;tmpq2eta=genantit_w_q2_eta;tmpq2phi=genantit_w_q2_phi;tmpbm=genantit_b_mass;tmpt_wtag=genantit_w_tag;}

            int insideqgr=0;
            for(int j=0;j<10;j++){
                if(dRgj[j][0]<rbmatching&&ptgengl[j]>0.15*tmptpt&&ptgengl[j]>extraqgt[0][0]) {extraqgt[0][0]=ptgengl[j];extraqgt[0][1]=etagengl[j];extraqgt[0][2]=phigengl[j];}
            }
            for(int j=0;j<5;j++){
                if((dRq1j[j][0]<rbmatching&&ptgenq1l[j]>0.15*tmptpt&&ptgenq1l[j]>extraqgt[0][0])){extraqgt[0][0]=ptgenq1l[j];extraqgt[0][1]=etagenq1l[j];extraqgt[0][2]=phigenq1l[j];}
                if((dRq2j[j][0]<rbmatching&&ptgenq2l[j]>0.15*tmptpt&&ptgenq2l[j]>extraqgt[0][0])){extraqgt[0][0]=ptgenq2l[j];extraqgt[0][1]=etagenq2l[j];extraqgt[0][2]=phigenq2l[j];}
                if((dRq3j[j][0]<rbmatching&&ptgenq3l[j]>0.15*tmptpt&&ptgenq3l[j]>extraqgt[0][0])){extraqgt[0][0]=ptgenq3l[j];extraqgt[0][1]=etagenq3l[j];extraqgt[0][2]=phigenq3l[j];}
                if((dRq4j[j][0]<rbmatching&&ptgenq4l[j]>0.15*tmptpt&&ptgenq4l[j]>extraqgt[0][0])){extraqgt[0][0]=ptgenq4l[j];extraqgt[0][1]=etagenq4l[j];extraqgt[0][2]=phigenq4l[j];}
                if((dRq5j[j][0]<rbmatching&&ptgenq5l[j]>0.15*tmptpt&&ptgenq5l[j]>extraqgt[0][0])){extraqgt[0][0]=ptgenq5l[j];extraqgt[0][1]=etagenq5l[j];extraqgt[0][2]=phigenq5l[j];}
            }
            dRt_bj[0]=deltaR(tmpbeta,tmpbphi,eta_max,phi_max);
                if(dRt_bj[0]<rbmatching){
            dRt_wq1j[0]=deltaR(tmpq1eta,tmpq1phi,eta_max,phi_max);
            dRt_wq2j[0]=deltaR(tmpq2eta,tmpq2phi,eta_max,phi_max);
            int tmpinr=0;
                if(extraqgt[0][0]>0) tmpinr++;
                if(dRt_bj[0]<rbmatching){extraqgt[tmpinr][0]=tmpbpt;extraqgt[tmpinr][1]=tmpbeta;extraqgt[tmpinr][2]=tmpbphi;tmpinr++;}
                if(dRt_wq1j[0]<rbmatching){extraqgt[tmpinr][0]=tmpq1pt;extraqgt[tmpinr][1]=tmpq1eta;extraqgt[tmpinr][2]=tmpq1phi;tmpinr++;}
                if(dRt_wq2j[0]<rbmatching&&tmpt_wtag==4){extraqgt[tmpinr][0]=tmpq2pt;extraqgt[tmpinr][1]=tmpq2eta;extraqgt[tmpinr][2]=tmpq2phi;tmpinr++;}
                int thorl=(tmpt_wtag==4)?1:-1;

            for(int k=0;k<4;k++) {matchnewt[k][0]=tmpinr;matchnewt[k][2]=tmpinr*thorl;
                if(extraqgt[0][0]>0)matchnewt[k][2]=(tmpinr-1)*thorl;
            }
                
                double PTtar=tmptpt;
                const int maxgroupr=tmpinr;
                int pgroupr[4][2*maxgroupr],subgroupr[4][maxgroupr][maxgroupr],nsubgroupr[4][maxgroupr],flagw1[4][2*maxgroupr],flagalone[4][maxgroupr];
                double ptgroupr[4][maxgroupr];
                for(int i=0;i<4;i++){
                    for(int j=0;j<maxgroupr;j++){
                        ptgroupr[i][j]=0;nsubgroupr[i][j]=-99;flagalone[i][j]=-99;
                        for(int k=0;k<maxgroupr;k++) subgroupr[i][j][k]=-99;}
                    for(int j=0;j<2*maxgroupr;j++) {pgroupr[i][j]=-99;flagw1[i][j]=-99;}
                }
                double drqqin5[2*tmpinr];
                for(int i=0;i<2*tmpinr;i++){drqqin5[i]=-99;}
                int idrqqin5=0;
                for(int m=0;m<tmpinr;m++){
                    for(int n=m+1;n<tmpinr;n++){
                        drqqin5[idrqqin5]=deltaR(extraqgt[m][1],extraqgt[m][2],extraqgt[n][1],extraqgt[n][2]);
                        //cout<<m<<"   "<<n<<"   "<<extraqgt[m][1]<<"   "<<extraqgt[m][2]<<"   "<<extraqgt[n][1]<<"   "<<extraqgt[n][2]<<"   "<<drqqin5[idrqqin5]<<endl;
                        idrqqin5++;
                    }
                }
                int i1sub,i2sub;
                double drwp[4]={0.1,0.15,0.2,0.25};
                int end1,end2;
                for(int i=0;i<4;i++){
                    i1sub=0;i2sub=0;
                    for(int j=0;j<2*tmpinr;j++){
                        if(0<drqqin5[j]&&drqqin5[j]<drwp[i]){
                            end1=endpoint(j,tmpinr,1);end2=endpoint(j,tmpinr,2);
                            //cout<<" check  "<<j<<"   "<<end1<<"   "<<end2<<"   "<<flagw1[i][end1]<<"   "<<flagw1[i][end2]<<" i1sub  "<<i1sub<<"   "<<i2sub<<"   "<<endl;
                            if(flagw1[i][end1]==-99&&flagw1[i][end2]==-99) {if(subgroupr[i][i1sub][0]!=-99){i1sub++;i2sub=0;}subgroupr[i][i1sub][i2sub]=end1;subgroupr[i][i1sub][i2sub+1]=end2;
                                flagw1[i][end1]=1;flagw1[i][end2]=1;i2sub+=2;
                                for(int k=0;k<2*tmpinr;k++){
                                    if(0<drqqin5[k]&&drqqin5[k]<drwp[i]){
                                        end1=endpoint(k,tmpinr,1);end2=endpoint(k,tmpinr,2);
                                        if((flagw1[i][end1]==-99&&flagw1[i][end2]==1)||(flagw1[i][end1]==1&&flagw1[i][end2]==-99)){
                                            int endtmp=end1;if(flagw1[i][end1]==1) endtmp=end2;
                                            subgroupr[i][i1sub][i2sub]=endtmp;i2sub++;flagw1[i][endtmp]=1;
                                            for(int m=0;m<2*tmpinr;m++){
                                                if(0<drqqin5[m]&&drqqin5[m]<drwp[i]){
                                                    end1=endpoint(m,tmpinr,1);end2=endpoint(m,tmpinr,2);
                                                    if((flagw1[i][end1]==-99&&flagw1[i][end2]==1)||(flagw1[i][end1]==1&&flagw1[i][end2]==-99)){
                                                        int endtmp=end1;if(flagw1[i][end1]==1) endtmp=end2;
                                                        subgroupr[i][i1sub][i2sub]=endtmp;i2sub++;flagw1[i][endtmp]=1;
                                                        for(int n=0;n<2*tmpinr;n++){
                                                            if(0<drqqin5[n]&&drqqin5[n]<drwp[i]){
                                                                end1=endpoint(n,tmpinr,1);end2=endpoint(n,tmpinr,2);
                                                                if((flagw1[i][end1]==-99&&flagw1[i][end2]==1)||(flagw1[i][end1]==1&&flagw1[i][end2]==-99)){
                                                                    int endtmp=end1;if(flagw1[i][end1]==1) endtmp=end2;
                                                                    subgroupr[i][i1sub][i2sub]=endtmp;i2sub++;flagw1[i][endtmp]=1;
                                                                }}}
                                                    }}}
                                        }
                                    }
                                }
                            }
                        }//end
                    }
                    int totalgroup=0;
                    for(int j=0;j<tmpinr;j++){if(subgroupr[i][j][0]!=-99)totalgroup++;}
                    const int ntotal=totalgroup;
                    int flagw2[ntotal][tmpinr];
                    for(int j=0;j<tmpinr;j++){
                        for(int k=0;k<ntotal;k++){
                            flagw2[k][j]=-99;
                            for(int m=0;m<tmpinr;m++){
                                if(subgroupr[i][k][m]==j){flagw2[k][j]=1;}
                            }}
                        int flagw3=-99;
                        for(int k=0;k<ntotal;k++){
                            if(flagw2[k][j]==1) flagw3=1;}
                        if(flagw3==-99){subgroupr[i][totalgroup][0]=j;totalgroup++;}
                    }
                }
                for(int i=0;i<4;i++){
                    for(int j=0;j<tmpinr;j++){
                        int igroup=0;
                        double tmpsubg[5][2]={0};
                        for(int k=0;k<tmpinr;k++){
                            if(subgroupr[i][j][k]!=-99) {tmpsubg[igroup][0]=extraqgt[subgroupr[i][j][k]][0];tmpsubg[igroup][1]=extraqgt[subgroupr[i][j][k]][2];igroup++;}
                        }
                        ptgroupr[i][j]=SumPt(tmpsubg,igroup);
                        if(ptgroupr[i][j]/PTtar>ignoreprong) {if(matchnewt[i][1]!=-99) matchnewt[i][1]++;if(matchnewt[i][1]==-99) matchnewt[i][1]=1;}
                        /*cout<<"groupt"<<j;
                         for(int k=0;k<tmpinr;k++)cout<<"   "<<subgroupr[i][j][k];
                         cout<<"   "<<ptgroupr[i][j]<<"    "<<ptgroupr[i][j]/PTtar<<"    "<<matchnewt[i][1]<<endl;*/
                    }
                    matchnewt[i][1]=matchnewt[i][1]*thorl;
                }
                }
        }
        }
        if(IsData>0&&(gentop_pt>0||genantitop_pt>0)){
            
            double extraqgw[3][3]={0};
            if (dRtj[0]<rmatching) matchingtmp[0]=-1;
            if (matchingtmp[0]!=-1&&dRantitj[0]<rmatching) matchingtmp[0]=-2;
            if (matchingtmp[0]==-1&&dRantitj[0]<rmatching&&dRtj[0]<dRantitj[0]) matchingtmp[0]=-1;
            if (matchingtmp[0]==-1&&dRantitj[0]<rmatching&&dRtj[0]>dRantitj[0]) matchingtmp[0]=-2;
            //cout<<i<<" jet  "<<Nj8<<"   "<<matchingtmp[0]<<endl;
            if(matchingtmp[0]==-1||matchingtmp[0]==-2){
                if(matchingtmp[0]==-1){tmptpt=gentop_pt;tmpteta=gentop_eta;tmptphi=gentop_phi;tmpbpt=gent_b_pt;tmpbeta=gent_b_eta;tmpbphi=gent_b_phi;tmpwpt=gent_w_pt;tmpweta=gent_w_eta;tmpwphi=gent_w_phi;tmpq1pt=gent_w_q1_pt;tmpq1eta=gent_w_q1_eta;tmpq1phi=gent_w_q1_phi;tmpq2pt=gent_w_q2_pt;tmpq2eta=gent_w_q2_eta;tmpq2phi=gent_w_q2_phi;tmpbm=gent_b_mass;tmpt_wtag=gent_w_tag;}
                if(matchingtmp[0]==-2){tmptpt=genantitop_pt;tmpteta=genantitop_eta;tmptphi=genantitop_phi;tmpbpt=genantit_b_pt;tmpbeta=genantit_b_eta;tmpbphi=genantit_b_phi;tmpwpt=genantit_w_pt;tmpweta=genantit_w_eta;tmpwphi=genantit_w_phi;tmpq1pt=genantit_w_q1_pt;tmpq1eta=genantit_w_q1_eta;tmpq1phi=genantit_w_q1_phi;tmpq2pt=genantit_w_q2_pt;tmpq2eta=genantit_w_q2_eta;tmpq2phi=genantit_w_q2_phi;tmpbm=genantit_b_mass;tmpt_wtag=genantit_w_tag;}
                
                int insideqgr=0;
                for(int j=0;j<10;j++){
                    if(dRgj[j][0]<rbmatching&&ptgengl[j]>0.15*tmpwpt&&ptgengl[j]>extraqgw[0][0]) {extraqgw[0][0]=ptgengl[j];extraqgw[0][1]=etagengl[j];extraqgw[0][2]=phigengl[j];}
                }
                for(int j=0;j<5;j++){
                    if((dRq1j[j][0]<rbmatching&&ptgenq1l[j]>0.15*tmpwpt&&ptgenq1l[j]>extraqgw[0][0])){extraqgw[0][0]=ptgenq1l[j];extraqgw[0][1]=etagenq1l[j];extraqgw[0][2]=phigenq1l[j];}
                    if((dRq2j[j][0]<rbmatching&&ptgenq2l[j]>0.15*tmpwpt&&ptgenq2l[j]>extraqgw[0][0])){extraqgw[0][0]=ptgenq2l[j];extraqgw[0][1]=etagenq2l[j];extraqgw[0][2]=phigenq2l[j];}
                    if((dRq3j[j][0]<rbmatching&&ptgenq3l[j]>0.15*tmpwpt&&ptgenq3l[j]>extraqgw[0][0])){extraqgw[0][0]=ptgenq3l[j];extraqgw[0][1]=etagenq3l[j];extraqgw[0][2]=phigenq3l[j];}
                    if((dRq4j[j][0]<rbmatching&&ptgenq4l[j]>0.15*tmpwpt&&ptgenq4l[j]>extraqgw[0][0])){extraqgw[0][0]=ptgenq4l[j];extraqgw[0][1]=etagenq4l[j];extraqgw[0][2]=phigenq4l[j];}
                    if((dRq5j[j][0]<rbmatching&&ptgenq5l[j]>0.15*tmpwpt&&ptgenq5l[j]>extraqgw[0][0])){extraqgw[0][0]=ptgenq5l[j];extraqgw[0][1]=etagenq5l[j];extraqgw[0][2]=phigenq5l[j];}
                }
                dRt_bj[0]=deltaR(tmpbeta,tmpbphi,eta_max,phi_max);
                dRt_wq1j[0]=deltaR(tmpq1eta,tmpq1phi,eta_max,phi_max);
                dRt_wq2j[0]=deltaR(tmpq2eta,tmpq2phi,eta_max,phi_max);
                if(dRt_bj[0]>rbmatching&&tmpt_wtag==4){
                int tmpinr=0;
                if(extraqgw[0][0]>0) tmpinr++;
                if(dRt_wq1j[0]<rbmatching){extraqgw[tmpinr][0]=tmpq1pt;extraqgw[tmpinr][1]=tmpq1eta;extraqgw[tmpinr][2]=tmpq1phi;tmpinr++;}
                if(dRt_wq2j[0]<rbmatching&&tmpt_wtag==4){extraqgw[tmpinr][0]=tmpq2pt;extraqgw[tmpinr][1]=tmpq2eta;extraqgw[tmpinr][2]=tmpq2phi;tmpinr++;}
                
                for(int k=0;k<4;k++) {matchneww[k][0]=tmpinr;matchneww[k][2]=tmpinr;
                    if(extraqgw[0][0]>0)matchneww[k][2]=tmpinr-1;
                }
                int thorl=(tmpt_wtag==4)?1:-1;
                
                    double PTtar=tmpwpt;
                    const int maxgroupr=tmpinr;
                    int pgroupr[4][2*maxgroupr],subgroupr[4][maxgroupr][maxgroupr],nsubgroupr[4][maxgroupr],flagw1[4][2*maxgroupr],flagalone[4][maxgroupr];
                    double ptgroupr[4][maxgroupr];
                    for(int i=0;i<4;i++){
                        for(int j=0;j<maxgroupr;j++){
                            ptgroupr[i][j]=0;nsubgroupr[i][j]=-99;flagalone[i][j]=-99;
                            for(int k=0;k<maxgroupr;k++) subgroupr[i][j][k]=-99;}
                        for(int j=0;j<2*maxgroupr;j++) {pgroupr[i][j]=-99;flagw1[i][j]=-99;}
                    }
                    double drqqin5[2*tmpinr];
                    for(int i=0;i<2*tmpinr;i++){drqqin5[i]=-99;}
                    int idrqqin5=0;
                    for(int m=0;m<tmpinr;m++){
                        for(int n=m+1;n<tmpinr;n++){
                            drqqin5[idrqqin5]=deltaR(extraqgw[m][1],extraqgw[m][2],extraqgw[n][1],extraqgw[n][2]);
                            //cout<<m<<"   "<<n<<"   "<<extraqgw[m][1]<<"   "<<extraqgw[m][2]<<"   "<<extraqgw[n][1]<<"   "<<extraqgw[n][2]<<"   "<<drqqin5[idrqqin5]<<endl;
                            idrqqin5++;
                        }
                    }
                    int i1sub,i2sub;
                    double drwp[4]={0.1,0.15,0.2,0.25};
                    int end1,end2;
                    for(int i=0;i<4;i++){
                        i1sub=0;i2sub=0;
                        for(int j=0;j<2*tmpinr;j++){
                            if(0<drqqin5[j]&&drqqin5[j]<drwp[i]){
                                end1=endpoint(j,tmpinr,1);end2=endpoint(j,tmpinr,2);
                                //cout<<" check  "<<j<<"   "<<end1<<"   "<<end2<<"   "<<flagw1[i][end1]<<"   "<<flagw1[i][end2]<<" i1sub  "<<i1sub<<"   "<<i2sub<<"   "<<endl;
                                if(flagw1[i][end1]==-99&&flagw1[i][end2]==-99) {if(subgroupr[i][i1sub][0]!=-99){i1sub++;i2sub=0;}subgroupr[i][i1sub][i2sub]=end1;subgroupr[i][i1sub][i2sub+1]=end2;
                                    flagw1[i][end1]=1;flagw1[i][end2]=1;i2sub+=2;
                                    for(int k=0;k<2*tmpinr;k++){
                                        if(0<drqqin5[k]&&drqqin5[k]<drwp[i]){
                                            end1=endpoint(k,tmpinr,1);end2=endpoint(k,tmpinr,2);
                                            if((flagw1[i][end1]==-99&&flagw1[i][end2]==1)||(flagw1[i][end1]==1&&flagw1[i][end2]==-99)){
                                                int endtmp=end1;if(flagw1[i][end1]==1) endtmp=end2;
                                                subgroupr[i][i1sub][i2sub]=endtmp;i2sub++;flagw1[i][endtmp]=1;
                                                for(int m=0;m<2*tmpinr;m++){
                                                    if(0<drqqin5[m]&&drqqin5[m]<drwp[i]){
                                                        end1=endpoint(m,tmpinr,1);end2=endpoint(m,tmpinr,2);
                                                        if((flagw1[i][end1]==-99&&flagw1[i][end2]==1)||(flagw1[i][end1]==1&&flagw1[i][end2]==-99)){
                                                            int endtmp=end1;if(flagw1[i][end1]==1) endtmp=end2;
                                                            subgroupr[i][i1sub][i2sub]=endtmp;i2sub++;flagw1[i][endtmp]=1;
                                                            for(int n=0;n<2*tmpinr;n++){
                                                                if(0<drqqin5[n]&&drqqin5[n]<drwp[i]){
                                                                    end1=endpoint(n,tmpinr,1);end2=endpoint(n,tmpinr,2);
                                                                    if((flagw1[i][end1]==-99&&flagw1[i][end2]==1)||(flagw1[i][end1]==1&&flagw1[i][end2]==-99)){
                                                                        int endtmp=end1;if(flagw1[i][end1]==1) endtmp=end2;
                                                                        subgroupr[i][i1sub][i2sub]=endtmp;i2sub++;flagw1[i][endtmp]=1;
                                                                    }}}
                                                        }}}
                                            }
                                        }
                                    }
                                }
                            }//end
                        }
                        int totalgroup=0;
                        for(int j=0;j<tmpinr;j++){if(subgroupr[i][j][0]!=-99)totalgroup++;}
                        const int ntotal=totalgroup;
                        int flagw2[ntotal][tmpinr];
                        for(int j=0;j<tmpinr;j++){
                            for(int k=0;k<ntotal;k++){
                                flagw2[k][j]=-99;
                                for(int m=0;m<tmpinr;m++){
                                    if(subgroupr[i][k][m]==j){flagw2[k][j]=1;}
                                }}
                            int flagw3=-99;
                            for(int k=0;k<ntotal;k++){
                                if(flagw2[k][j]==1) flagw3=1;}
                            if(flagw3==-99){subgroupr[i][totalgroup][0]=j;totalgroup++;}
                        }
                    }
                    for(int i=0;i<4;i++){
                        for(int j=0;j<tmpinr;j++){
                            int igroup=0;
                            double tmpsubg[5][2]={0};
                            for(int k=0;k<tmpinr;k++){
                                if(subgroupr[i][j][k]!=-99) {tmpsubg[igroup][0]=extraqgw[subgroupr[i][j][k]][0];tmpsubg[igroup][1]=extraqgw[subgroupr[i][j][k]][2];igroup++;}
                            }
                            ptgroupr[i][j]=SumPt(tmpsubg,igroup);
                            if(ptgroupr[i][j]/PTtar>ignoreprong) {if(matchneww[i][1]!=-99) matchneww[i][1]++;if(matchneww[i][1]==-99) matchneww[i][1]=1;}
                            /*cout<<"group"<<j;
                             for(int k=0;k<tmpinr;k++)cout<<"   "<<subgroupr[i][j][k];
                             cout<<"   "<<ptgroupr[i][j]<<"    "<<ptgroupr[i][j]/PTtar<<"    "<<matchneww[i][1]<<endl;*/
                        }
                    }

            }
        }
        }
        if(IsData>0&&(gentop_pt<0&&genantitop_pt<0)){
            
            double extraqgw[3][3]={0};
            double allw[5][9],dRqq_wall[5];
            for(int i=0;i<5;i++){
                for(int j=0;j<9;j++) allw[i][j]=-99;
                if(taggenwl[i]==4){
                allw[i][0]=ptgenwl[i];allw[i][1]=etagenwl[i];allw[i][2]=phigenwl[i];
                allw[i][3]=genw_q1_pt[i];allw[i][4]=genw_q1_eta[i];allw[i][5]=genw_q1_phi[i];
                    allw[i][6]=genw_q2_pt[i];allw[i][7]=genw_q2_eta[i];allw[i][8]=genw_q2_phi[i];}
                dRqq_wall[i]=deltaR(allw[i][1],allw[i][1],eta_max,phi_max);
            }
            Int_t *index_wall=new Int_t[5];
            TMath::Sort(5,dRqq_wall,index_wall,0);

                int insideqgr=0;
                for(int j=0;j<10;j++){
                    if(dRgj[j][0]<rbmatching&&ptgengl[j]>0.15*allw[index_wall[0]][0]&&ptgengl[j]>extraqgw[0][0]) {extraqgw[0][0]=ptgengl[j];extraqgw[0][1]=etagengl[j];extraqgw[0][2]=phigengl[j];}
                }
                for(int j=0;j<5;j++){
                    if((dRq1j[j][0]<rbmatching&&ptgenq1l[j]>0.15*allw[index_wall[0]][0]&&ptgenq1l[j]>extraqgw[0][0])){extraqgw[0][0]=ptgenq1l[j];extraqgw[0][1]=etagenq1l[j];extraqgw[0][2]=phigenq1l[j];}
                    if((dRq2j[j][0]<rbmatching&&ptgenq2l[j]>0.15*allw[index_wall[0]][0]&&ptgenq2l[j]>extraqgw[0][0])){extraqgw[0][0]=ptgenq2l[j];extraqgw[0][1]=etagenq2l[j];extraqgw[0][2]=phigenq2l[j];}
                    if((dRq3j[j][0]<rbmatching&&ptgenq3l[j]>0.15*allw[index_wall[0]][0]&&ptgenq3l[j]>extraqgw[0][0])){extraqgw[0][0]=ptgenq3l[j];extraqgw[0][1]=etagenq3l[j];extraqgw[0][2]=phigenq3l[j];}
                    if((dRq4j[j][0]<rbmatching&&ptgenq4l[j]>0.15*allw[index_wall[0]][0]&&ptgenq4l[j]>extraqgw[0][0])){extraqgw[0][0]=ptgenq4l[j];extraqgw[0][1]=etagenq4l[j];extraqgw[0][2]=phigenq4l[j];}
                    if((dRq5j[j][0]<rbmatching&&ptgenq5l[j]>0.15*allw[index_wall[0]][0]&&ptgenq5l[j]>extraqgw[0][0])){extraqgw[0][0]=ptgenq5l[j];extraqgw[0][1]=etagenq5l[j];extraqgw[0][2]=phigenq5l[j];}
                }
            double drwq1j=deltaR(allw[index_wall[0]][4],allw[index_wall[0]][5],eta_max,phi_max);
            double drwq2j=deltaR(allw[index_wall[0]][7],allw[index_wall[0]][8],eta_max,phi_max);

                    int tmpinr=0;
                    if(extraqgw[0][0]>0) tmpinr++;
                    if(drwq1j<rbmatching){extraqgw[tmpinr][0]=allw[index_wall[0]][3];extraqgw[tmpinr][1]=allw[index_wall[0]][4];extraqgw[tmpinr][2]=allw[index_wall[0]][5];tmpinr++;}
                    if(drwq2j<rbmatching){extraqgw[tmpinr][0]=allw[index_wall[0]][6];extraqgw[tmpinr][1]=allw[index_wall[0]][7];extraqgw[tmpinr][2]=allw[index_wall[0]][8];tmpinr++;}
                    
                    for(int k=0;k<4;k++) {matchneww[k][0]=tmpinr;matchneww[k][2]=tmpinr;
                        if(extraqgw[0][0]>0)matchneww[k][2]=tmpinr-1;
                    }
                    int thorl=1;
                    
            double PTtar=allw[index_wall[0]][0];
            const int maxgroupr=tmpinr;
            int pgroupr[4][2*maxgroupr],subgroupr[4][maxgroupr][maxgroupr],nsubgroupr[4][maxgroupr],flagw1[4][2*maxgroupr],flagalone[4][maxgroupr];
            double ptgroupr[4][maxgroupr];
            for(int i=0;i<4;i++){
                for(int j=0;j<maxgroupr;j++){
                    ptgroupr[i][j]=0;nsubgroupr[i][j]=-99;flagalone[i][j]=-99;
                    for(int k=0;k<maxgroupr;k++) subgroupr[i][j][k]=-99;}
                for(int j=0;j<2*maxgroupr;j++) {pgroupr[i][j]=-99;flagw1[i][j]=-99;}
            }
            double drqqin5[2*tmpinr];
            for(int i=0;i<2*tmpinr;i++){drqqin5[i]=-99;}
            int idrqqin5=0;
            for(int m=0;m<tmpinr;m++){
                for(int n=m+1;n<tmpinr;n++){
                    drqqin5[idrqqin5]=deltaR(extraqgw[m][1],extraqgw[m][2],extraqgw[n][1],extraqgw[n][2]);
                    //cout<<m<<"   "<<n<<"   "<<extraqgw[m][1]<<"   "<<extraqgw[m][2]<<"   "<<extraqgw[n][1]<<"   "<<extraqgw[n][2]<<"   "<<drqqin5[idrqqin5]<<endl;
                    idrqqin5++;
                }
            }
            int i1sub,i2sub;
            double drwp[4]={0.1,0.15,0.2,0.25};
            int end1,end2;
            for(int i=0;i<4;i++){
                i1sub=0;i2sub=0;
                for(int j=0;j<2*tmpinr;j++){
                    if(0<drqqin5[j]&&drqqin5[j]<drwp[i]){
                        end1=endpoint(j,tmpinr,1);end2=endpoint(j,tmpinr,2);
                        //cout<<" check  "<<j<<"   "<<end1<<"   "<<end2<<"   "<<flagw1[i][end1]<<"   "<<flagw1[i][end2]<<" i1sub  "<<i1sub<<"   "<<i2sub<<"   "<<endl;
                        if(flagw1[i][end1]==-99&&flagw1[i][end2]==-99) {if(subgroupr[i][i1sub][0]!=-99){i1sub++;i2sub=0;}subgroupr[i][i1sub][i2sub]=end1;subgroupr[i][i1sub][i2sub+1]=end2;
                            flagw1[i][end1]=1;flagw1[i][end2]=1;i2sub+=2;
                            for(int k=0;k<2*tmpinr;k++){
                                if(0<drqqin5[k]&&drqqin5[k]<drwp[i]){
                                    end1=endpoint(k,tmpinr,1);end2=endpoint(k,tmpinr,2);
                                    if((flagw1[i][end1]==-99&&flagw1[i][end2]==1)||(flagw1[i][end1]==1&&flagw1[i][end2]==-99)){
                                        int endtmp=end1;if(flagw1[i][end1]==1) endtmp=end2;
                                        subgroupr[i][i1sub][i2sub]=endtmp;i2sub++;flagw1[i][endtmp]=1;
                                        for(int m=0;m<2*tmpinr;m++){
                                            if(0<drqqin5[m]&&drqqin5[m]<drwp[i]){
                                                end1=endpoint(m,tmpinr,1);end2=endpoint(m,tmpinr,2);
                                                if((flagw1[i][end1]==-99&&flagw1[i][end2]==1)||(flagw1[i][end1]==1&&flagw1[i][end2]==-99)){
                                                    int endtmp=end1;if(flagw1[i][end1]==1) endtmp=end2;
                                                    subgroupr[i][i1sub][i2sub]=endtmp;i2sub++;flagw1[i][endtmp]=1;
                                                    for(int n=0;n<2*tmpinr;n++){
                                                        if(0<drqqin5[n]&&drqqin5[n]<drwp[i]){
                                                            end1=endpoint(n,tmpinr,1);end2=endpoint(n,tmpinr,2);
                                                            if((flagw1[i][end1]==-99&&flagw1[i][end2]==1)||(flagw1[i][end1]==1&&flagw1[i][end2]==-99)){
                                                                int endtmp=end1;if(flagw1[i][end1]==1) endtmp=end2;
                                                                subgroupr[i][i1sub][i2sub]=endtmp;i2sub++;flagw1[i][endtmp]=1;
                                                            }}}
                                                }}}
                                    }
                                }
                            }
                        }
                    }//end
                }
                int totalgroup=0;
                for(int j=0;j<tmpinr;j++){if(subgroupr[i][j][0]!=-99)totalgroup++;}
                const int ntotal=totalgroup;
                int flagw2[ntotal][tmpinr];
                for(int j=0;j<tmpinr;j++){
                    for(int k=0;k<ntotal;k++){
                        flagw2[k][j]=-99;
                        for(int m=0;m<tmpinr;m++){
                            if(subgroupr[i][k][m]==j){flagw2[k][j]=1;}
                        }}
                    int flagw3=-99;
                    for(int k=0;k<ntotal;k++){
                        if(flagw2[k][j]==1) flagw3=1;}
                    if(flagw3==-99){subgroupr[i][totalgroup][0]=j;totalgroup++;}
                }
            }
            for(int i=0;i<4;i++){
                for(int j=0;j<tmpinr;j++){
                    int igroup=0;
                    double tmpsubg[5][2]={0};
                    for(int k=0;k<tmpinr;k++){
                        if(subgroupr[i][j][k]!=-99) {tmpsubg[igroup][0]=extraqgw[subgroupr[i][j][k]][0];tmpsubg[igroup][1]=extraqgw[subgroupr[i][j][k]][2];igroup++;}
                    }
                    ptgroupr[i][j]=SumPt(tmpsubg,igroup);
                    if(ptgroupr[i][j]/PTtar>ignoreprong) {if(matchneww[i][1]!=-99) matchneww[i][1]++;if(matchneww[i][1]==-99) matchneww[i][1]=1;}
                    /*cout<<"group"<<j;
                     for(int k=0;k<tmpinr;k++)cout<<"   "<<subgroupr[i][j][k];
                     cout<<"   "<<ptgroupr[i][j]<<"    "<<ptgroupr[i][j]/PTtar<<"    "<<matchneww[i][1]<<endl;*/
                }
            }
                    }
                    
        //if(jentry==10000) break;
        
		Double_t isAnaHP=1.;
		Double_t isAnaLP=1.;
		Double_t isAnaNP=1.;
		Double_t isTTBarControl=1.;
		Int_t tmp_categoryID_channel=0;
		double pt_cut=200;
		if( channelname=="el" ){
		tmp_categoryID_channel=-1;// -1 for el; 1 for mu
            if(IsData>0){
			if ( lep==11 && (HLT_Ele4>0 || HLT_Ele5>0 || HLT_Ele8>0) && Mj_max > 40 && nLooseLep==1 && ptlep1>55 && fabs(etalep1)<2.5 && MET_et>80 && ptVlepJEC > 200. && jet_pt_puppi>200 && fabs(jetAK8puppi_eta)<2.4 && ((IDLoose>0&&Nj8==1)||(IDLoose>0&&IDLoose_2>0&&Nj8==2&& fabs(jetAK8puppi_eta_2)<2.4&&jet_pt_puppi_2>pt_cut)||(IDLoose>0&&IDLoose_2>0&&IDLoose_3>0&&Nj8>=3&& fabs(jetAK8puppi_eta_2)<2.4&& fabs(jetAK8puppi_eta_3)<2.4&&jet_pt_puppi_2>pt_cut&&jet_pt_puppi_3>pt_cut))&&passFilter_HBHE>0  &&  passFilter_GlobalHalo>0 && passFilter_HBHEIso>0 && passFilter_ECALDeadCell>0 && passFilter_GoodVtx>0 && passFilter_badMuon>0>0)// &&  passFilter_EEBadSc>0
                ExTree->Fill();}
            if(IsData==0&&IsMET==0){
                if ( lep==11 && (HLT_Ele4>0 || HLT_Ele5>0 || HLT_Ele8>0) && Mj_max > 40 && nLooseLep==1 && ptlep1>55 && fabs(etalep1)<2.5 && MET_et>80 && ptVlepJEC > 200. && jet_pt_puppi>200 && fabs(jetAK8puppi_eta)<2.4 && ((IDLoose>0&&Nj8==1)||(IDLoose>0&&IDLoose_2>0&&Nj8==2&& fabs(jetAK8puppi_eta_2)<2.4&&jet_pt_puppi_2>pt_cut)||(IDLoose>0&&IDLoose_2>0&&IDLoose_3>0&&Nj8>=3&& fabs(jetAK8puppi_eta_2)<2.4&& fabs(jetAK8puppi_eta_3)<2.4&&jet_pt_puppi_2>pt_cut&&jet_pt_puppi_3>pt_cut))&&passFilter_HBHE>0  &&  passFilter_GlobalHalo>0 && passFilter_HBHEIso>0 && passFilter_ECALDeadCell>0 && passFilter_GoodVtx>0 && passFilter_badMuon>0>0&&  passFilter_EEBadSc>0)// &&  passFilter_EEBadSc>0
                    ExTree->Fill();}
            /*if(IsData==0&&IsMET==1){
                if ( lep==11 && (HLT_Mu12>0 && HLT_Mu2<1 && HLT_Mu3<1 && HLT_Ele3<1 && HLT_Ele4<1 && HLT_Ele8<1) && Mj_max > 40 && nLooseLep==1 && ptlep1>55 && fabs(etalep1)<2.5 && MET_et>80 && ptVlepJEC > 200. && jet_pt_puppi>200 && fabs(jetAK8puppi_eta)<2.4 && ((IDLoose>0&&Nj8==1)||(IDLoose>0&&IDLoose_2>0&&Nj8==2&& fabs(jetAK8puppi_eta_2)<2.4)||(IDLoose>0&&IDLoose_2>0&&IDLoose_3>0&&Nj8>=3&& fabs(jetAK8puppi_eta_2)<2.4&& fabs(jetAK8puppi_eta_3)<2.4))&&passFilter_HBHE>0  &&  passFilter_GlobalHalo>0 && passFilter_HBHEIso>0 && passFilter_ECALDeadCell>0 && passFilter_GoodVtx>0 && passFilter_badMuon>0&&  passFilter_EEBadSc>0)// &&  passFilter_EEBadSc>0
                    ExTree->Fill();}*/
        }
		else if( channelname=="mu" ){
		tmp_categoryID_channel=1;// -1 for el; 1 for mu
            if(IsData>0){
            ExTree->Fill();
        //    if(lep==13 &&Nj8==2 ){    
        //    if(Nj8==2 ){    
       //     if(lep==13 && (HLT_Mu2>0 || HLT_Mu3>0|| HLT_Mu12>0) && trackIso/ptlep1<0.1 && muisolation<0.05 && fabs(etalep1)<2.4 && Mj_max > 40 && nLooseLep==1 && ptlep1>55 && MET_et>40 && ptVlepJEC>200 && jet_pt_puppi>200 && fabs(jetAK8puppi_eta)<2.4 && ((IDLoose>0&&Nj8==1)||(IDLoose>0&&IDLoose_2>0&&Nj8==2&& fabs(jetAK8puppi_eta_2)<2.4&&jet_pt_puppi_2>pt_cut)||(IDLoose>0&&IDLoose_2>0&&IDLoose_3>0&&Nj8>=3&& fabs(jetAK8puppi_eta_2)<2.4&& fabs(jetAK8puppi_eta_3)<2.4&&jet_pt_puppi_2>pt_cut&&jet_pt_puppi_3>pt_cut))&&passFilter_HBHE>0  &&  passFilter_GlobalHalo>0 && passFilter_HBHEIso>0 && passFilter_ECALDeadCell>0 && passFilter_GoodVtx>0 && passFilter_badMuon>0>0)
        //        ExTree->Fill();}
}
             if(IsData==0&&IsMET==0){
            if(lep==13 &&Nj8==2 )
       //         if(lep==13 && (HLT_Mu2>0 || HLT_Mu3>0) && trackIso/ptlep1<0.1 && muisolation<0.05 && fabs(etalep1)<2.4 && Mj_max > 40 && nLooseLep==1 && ptlep1>55 && MET_et>40 && ptVlepJEC>200 && jet_pt_puppi>200 && fabs(jetAK8puppi_eta)<2.4 && ((IDLoose>0&&Nj8==1)||(IDLoose>0&&IDLoose_2>0&&Nj8==2&& fabs(jetAK8puppi_eta_2)<2.4&&jet_pt_puppi_2>pt_cut)||(IDLoose>0&&IDLoose_2>0&&IDLoose_3>0&&Nj8>=3&& fabs(jetAK8puppi_eta_2)<2.4&& fabs(jetAK8puppi_eta_3)<2.4&&jet_pt_puppi_2>pt_cut&&jet_pt_puppi_3>pt_cut))&&passFilter_HBHE>0  &&  passFilter_GlobalHalo>0 && passFilter_HBHEIso>0 && passFilter_ECALDeadCell>0 && passFilter_GoodVtx>0 && passFilter_badMuon>0>0&&  passFilter_EEBadSc>0)
                    ExTree->Fill();}
            if(IsData==0&&IsMET==1){
            if(lep==13 &&Nj8==2 )
       //         if(lep==13 && (HLT_Mu2<1&&HLT_Mu3<1&& HLT_Mu12>0&&HLT_Ele4<1&& HLT_Ele5<1&& HLT_Ele8<1) && trackIso/ptlep1<0.1 && muisolation<0.05 && fabs(etalep1)<2.4 && Mj_max > 40 && nLooseLep==1 && ptlep1>55 && MET_et>40 && ptVlepJEC>200 && jet_pt_puppi>200 && fabs(jetAK8puppi_eta)<2.4 && ((IDLoose>0&&Nj8==1)||(IDLoose>0&&IDLoose_2>0&&Nj8==2&& fabs(jetAK8puppi_eta_2)<2.4&&jet_pt_puppi_2>pt_cut)||(IDLoose>0&&IDLoose_2>0&&IDLoose_3>0&&Nj8>=3&& fabs(jetAK8puppi_eta_2)<2.4&& fabs(jetAK8puppi_eta_3)<2.4&&jet_pt_puppi_2>pt_cut&&jet_pt_puppi_3>pt_cut))&&passFilter_HBHE>0  &&  passFilter_GlobalHalo>0 && passFilter_HBHEIso>0 && passFilter_ECALDeadCell>0 && passFilter_GoodVtx>0 && passFilter_badMuon>0>0&&  passFilter_EEBadSc>0)

                //if(lep==13 && (HLT_Mu2<1&&HLT_Mu3<1&& HLT_Mu12>0) && trackIso/ptlep1<0.1 && muisolation<0.05 && fabs(etalep1)<2.4 && Mj_max > 40 && nLooseLep==1 && ptlep1>55 && MET_et>40 && ptVlepJEC>200 && jet_pt_puppi>200 && fabs(jetAK8puppi_eta)<2.4 && ((IDLoose>0&&Nj8==1)||(IDLoose>0&&IDLoose_2>0&&Nj8==2&& fabs(jetAK8puppi_eta_2)<2.4)||(IDLoose>0&&IDLoose_2>0&&IDLoose_3>0&&Nj8>=3&& fabs(jetAK8puppi_eta_2)<2.4&& fabs(jetAK8puppi_eta_3)<2.4))&&passFilter_HBHE>0  &&  passFilter_GlobalHalo>0 && passFilter_HBHEIso>0 && passFilter_ECALDeadCell>0 && passFilter_GoodVtx>0 && passFilter_badMuon>0>0&&  passFilter_EEBadSc>0)
                    ExTree->Fill();}

		}else{
			cout<<"We don't know channelname:"<<channelname<<endl;
		}


//cout<<nLooseMu<<"      "<<nLooseEle<<"     "<<endl;
	}
    elereco->Close();
    elehlt->Close();
    muhlt->Close();
    methlt->Close();
    muid->Close();
    muiso->Close();
    muhlt1->Close();
    muid1->Close();
    muiso1->Close();

}

/*
 2019/08/15:
 add dR(b/2q,j)<0.6:t matched
 add dR(2q,j)&dR(b,j)<0.6:W matched
 judge t or W; jet mass<110 ->100
*/
