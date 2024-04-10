void GenerateSamples()
{

 ifstream infile;
 infile.open("Combination_12_6.txt");
 if(!infile) cout<<"error"<<endl;

 string str;
 Int_t t1;
 Int_t a[924][6];
 Int_t*p=&a[0][0];
 while(infile>>t1)             
 {
   *p=t1;
    p++;
 }
 infile.close();

 //TFile *file1 = TFile::Open("../SubSamples/AnalysisResults_fixed_1.root","READ");
 //TFile *file2 = TFile::Open("../SubSamples/AnalysisResults_fixed_2.root","READ"); 

 TString ifName;
 TFileMerger merger;

 for(Int_t k=1;k<=924;++k)
 {
  for(Int_t j=0;j<6;++j)
  {
   ifName = Form("/Users/zbiao/noprompt_study/DmesonAnalysis/filterdata/JK12/AnalysisResul%d.root",a[k-1][j]);
   merger.AddFile(ifName.Data(),kFALSE); 
  }
  //merger.OutputFile(Form("./Samples/Samples_%d.root",k)); 
  merger.OutputFile(Form("jack_output/12_6/Samples_%d.root",k));
  merger.Merge(); 
  merger.ClearObjectNames();
  cout <<"Merge " << k << endl;
 }

 



}

