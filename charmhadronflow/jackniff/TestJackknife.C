void TestJackknife(const Int_t split=12,const Int_t ngroups=6){

  Double_t average=500.;
  Double_t variance=5.;
  
  Double_t reduceSampleVariance=variance*TMath::Sqrt(split);

  Double_t values[split];
  Double_t normalAverage=0;
  for(Int_t i=0;i<split;i++){
    values[i]=gRandom->Gaus(average,reduceSampleVariance);
    normalAverage+=values[i];
  }
  normalAverage/=split;
  Printf("normalAverage is: %f, ref. variance is %f",normalAverage,variance);

  TH1D *hAverages=new TH1D("hAverages","hAverages",400,average-40*variance,average+40*variance);
  
  // now jacknife
  Int_t minInt=(Int_t)TMath::Power(2.,ngroups)-1;
  Int_t maxInt=(Int_t)TMath::Power(2.,split)-1;

  for(Int_t k=minInt;k<=maxInt;k++){
    Double_t avTrial=0;
    Int_t sumbits=0;
    for(Int_t jbit=0;jbit<=split;jbit++){
      if( k & (1<<jbit)){
	avTrial+=values[jbit];
	sumbits++;
      }
    }
    if(sumbits==ngroups){
      avTrial/=ngroups;
      hAverages->Fill(avTrial);
    }
  }

  TCanvas *cJK=new TCanvas("cJK","cJK",800,800);
  cJK->cd();
  hAverages->Draw();
}
