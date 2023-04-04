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

void makeCutsTreeCreator(bool fIncludeDs=true, bool fIncludeDplus=false, int cent=k010, bool fIsMC=false, double ptDsmin=2., double ptDsmax=50., int preselType=kTightQM)
{
    if(!fIncludeDs && !fIncludeDplus) {
        std::cerr << "You have to enable at least a meson species! Exit " << std::endl;
        return;
    }

    AliRDHFCutsDstoKKpi* cutsDsFilt = NULL;
    AliRDHFCutsDstoKKpi* cutsDsCent = NULL;

    AliRDHFCutsDplustoKpipi* cutsDplusFilt = NULL;
    AliRDHFCutsDplustoKpipi* cutsDplusCent = NULL;

    TString centname = "";
    TString triggername = "";
    switch(cent) {
        case k010:
        {
            centname="010";
            triggername="kINT7_kCentral";
            if(fIncludeDs) {
                cutsDsFilt = MakeFileForCutsDs010_FiltTreeCreator2018QM(fIsMC, ptDsmin, ptDsmax, preselType);
                cutsDsCent = MakeFileForCutsDs010_Central2018(true, 8.0, fIsMC, 0, ptDsmin, ptDsmax);
            }
            if(fIncludeDplus) {
                cutsDplusFilt = MakeFileForCutsDplus010_FiltTreeCreator2018(false, -1.0, fIsMC);
                cutsDplusCent = MakeFileForCutsDplus010_Central2018(true, 3.0, fIsMC, 0);
            }

            break;
        }
        case k3050:
        {
            centname="3050";
            triggername="kINT7_kSemiCentral";
            if(fIncludeDs) {
                cutsDsFilt = MakeFileForCutsDs3050_Filt2018_Pass3(fIsMC, ptDsmin, ptDsmax);
                cutsDsCent = MakeFileForCutsDs3050_Central2018_Pass3(true, 8.0, fIsMC, 0, ptDsmin, ptDsmax);
            }
            if(fIncludeDplus) {
                cutsDplusFilt = MakeFileForCutsDplus3050_FiltTreeCreator2018(false, -1.0, fIsMC);
                cutsDplusCent = MakeFileForCutsDplus3050_Central2018(false, -1.0, fIsMC, 0);
            }

            break;
        }
        case k6080:
        {
            centname="6080";
            triggername="kINT7";
            if(fIncludeDs) {
                cutsDsFilt = MakeFileForCutsDs6080_FiltTreeCreator2018(false, -1.0, fIsMC, ptDsmin, ptDsmax);
                cutsDsCent = MakeFileForCutsDs6080_Central2018(true, 8.0, fIsMC, 0, ptDsmin, ptDsmax);
            }
            if(fIncludeDplus) {
                cutsDplusFilt = MakeFileForCutsDplus3050_FiltTreeCreator2018(false, -1.0, fIsMC);
                cutsDplusCent = MakeFileForCutsDplus3050_Central2018(false, -1.0, fIsMC, 0);
            }

            break;
        }
        default:
        {
            std::cerr << "Centrality not implemented! Exit " << std::endl;
            return;
        }

    }

    if(fIsMC) triggername = "kMB";
    TString mesonname = "";
    if(fIncludeDs) mesonname += "Ds";
    if(fIncludeDplus) mesonname += "Dplus";
    TString preselname = "";
    if(cent==k010)
    {
        if(preselType == kLooseQM)
            preselname = "_loosepresel";
        else if(preselType == kMediumQM)
            preselname = "_mediumpresel";
        else if(preselType == kTightQM)
            preselname = "_tightpresel";
        else
            preselname = "_undefinedpresel";
    }

    TFile fout(Form("%sCuts_treecreator_PbPb2018_%s_%s_ptDs%0.f_%0.f%s.root",mesonname.Data(), centname.Data(), triggername.Data(), ptDsmin, ptDsmax, preselname.Data()),"recreate");
    fout.cd();
    if(fIncludeDs) {
        cutsDsFilt->SetName("DstoKKpiFilteringCuts");
        cutsDsCent->SetName("DstoKKpiAnalysisCuts");
        cutsDsFilt->Write("DstoKKpiFilteringCuts");
        cutsDsCent->Write("DstoKKpiAnalysisCuts");
    }
    if(fIncludeDplus) {
        cutsDplusFilt->SetName("DplustoKpipiFilteringCuts");
        cutsDplusCent->SetName("DplustoKpipiAnalysisCuts");
        cutsDplusFilt->Write("DplustoKpipiFilteringCuts");
        cutsDplusCent->Write("DplustoKpipiAnalysisCuts");
    }
    fout.Close();
}
