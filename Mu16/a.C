void a(){
Double_t b[3]={2,1,0},
                     c[3]={3,6,5};
Int_t *indexmass=new Int_t[3];
                TMath::Sort(2,b,indexmass,1);
                cout<<indexmass[1]<<" "<<b[indexmass[0]]<<endl;
                cout<<c[indexmass[0]]<<endl;
}


