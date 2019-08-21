#if !defined (__CINT__) || defined (__CLING__)

#include <Riostream.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TParameter.h>

#include <AliRDHFCutsDstoKKpi.h>
#include <AliRDHFCutsDplustoKpipi.h>

#include "CutsDstoKKpi_PbPb2018_010_Raa.C"
#include "CutsDstoKKpi_PbPb2018_3050_Raa.C"
#include "CutsDstoKKpi_PbPb2018_6080_Raa.C"

#include "CutsDplustoKpipi_PbPb2018_010_Raa.C"
#include "CutsDplustoKpipi_PbPb2018_3050_Raa.C"
#include "CutsDplustoKpipi_PbPb2018_6080_Raa.C"

#endif

enum centclasses {k010, k3050, k6080};

void makeCutsTreeCreator(bool fIncludeDs=true, bool fIncludeDplus=false, int cent=k010, bool fIsMC=false)
{
    if(!fIncludeDs && !fIncludeDplus) {
        std::cerr << "You have to enable at least a meson species! Exit " << std::endl;
        return;
    }

    AliRDHFCutsDstoKKpi* cutsDsFilt = NULL;
    AliRDHFCutsDstoKKpi* cutsDsCent = NULL;

    AliRDHFCutsDplustoKpipi* cutsDplusFilt = NULL;
    AliRDHFCutsDplustoKpipi* cutsDplusCent = NULL;

    switch(cent) {
        case k010:
        {
            if(fIncludeDs) {
                cutsDsFilt = MakeFileForCutsDs010_FiltTreeCreator2018(false, -1.0, fIsMC);
                cutsDsCent = MakeFileForCutsDs010_Central2018(true, 8.0, fIsMC);
            }
            if(fIncludeDplus) {
                cutsDplusFilt = MakeFileForCutsDplus010_FiltTreeCreator2018(false, -1.0, fIsMC);
                cutsDplusCent = MakeFileForCutsDplus010_Central2018(true, 3.0, fIsMC);
            }

            break;
        }
        case k3050:
        {
            if(fIncludeDs) {
                cutsDsFilt = MakeFileForCutsDs3050_FiltTreeCreator2018(false, -1.0, fIsMC);
                cutsDsCent = MakeFileForCutsDs3050_Central2018(true, 8.0, fIsMC);
            }
            if(fIncludeDplus) {
                cutsDplusFilt = MakeFileForCutsDplus3050_FiltTreeCreator2018(false, -1.0, fIsMC);
                cutsDplusCent = MakeFileForCutsDplus3050_Central2018(false, -1.0, fIsMC);
            }

            break;
        }
        case k6080:
        {
            if(fIncludeDs) {
                cutsDsFilt = MakeFileForCutsDs6080_FiltTreeCreator2018(false, -1.0, fIsMC);
                cutsDsCent = MakeFileForCutsDs6080_Central2018(true, 8.0, fIsMC);
            }
            if(fIncludeDplus) {
                cutsDplusFilt = MakeFileForCutsDplus3050_FiltTreeCreator2018(false, -1.0, fIsMC);
                cutsDplusCent = MakeFileForCutsDplus3050_Central2018(false, -1.0, fIsMC);
            }

            break;
        }

    }

    TString triggername = "kINT7_kCentral";
    if(fIsMC) triggername = "kMB";
    TString mesonname = "";
    if(fIncludeDs) mesonname += "Ds";
    if(fIncludeDplus) mesonname += "Dplus";

    TFile fout(Form("%sCuts_PbPb2018_Central_%s.root",mesonname.Data(), triggername.Data()),"recreate");
    fout.cd();
    if(fIncludeDs) {
        cutsDsFilt->Write("DstoKKpiFilteringCuts");
        cutsDsCent->Write("DstoKKpiAnalysisCuts");
    }
    if(fIncludeDplus) {
        cutsDplusFilt->Write("DplustoKpopiFilteringCuts");
        cutsDplusCent->Write("DplustoKpipiAnalysisCuts");
    }
    fout.Close();    
}
