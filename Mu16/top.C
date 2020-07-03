void top(){
	double ToppTweight=1;
	double gentop_pt=250;
	double genantitop_pt=250;
        double a_top=0.0615;
        double b_top=- 0.0005;
        if (gentop_pt<0 || genantitop_pt<0){ToppTweight=1.0;}
        else if(gentop_pt>0 && genantitop_pt>0){ToppTweight=sqrt(exp(a_top+b_top*gentop_pt)*exp(a_top+b_top*genantitop_pt));}
        else if(gentop_pt>0 && genantitop_pt<0){ToppTweight=sqrt(exp(a_top+b_top*gentop_pt));}
        else if(gentop_pt<0 && genantitop_pt>0){ToppTweight=sqrt(exp(a_top+b_top*genantitop_pt));}
	cout<<ToppTweight<<endl;
}
