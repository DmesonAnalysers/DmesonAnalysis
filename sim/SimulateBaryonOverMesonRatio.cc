#if defined(__tum_batchfarm__)
R__ADD_INCLUDE_PATH(/home/software/pythia/pythia8307/include/)
R__LOAD_LIBRARY(/home/software/pythia/pythia8307/lib/libpythia8.so)
#endif

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TFile.h>
#include <TH1.h>

#include "Pythia8/Pythia.h"
#include "TNtuple.h"

#endif

using namespace Pythia8;

namespace {
enum tunes { kMonash = 0, kCRMode0, kCRMode2, kCRMode3 };
enum processes { kSoftQCD = 0, kHardQCD };
}  // namespace

//__________________________________________________________________________________________________
template <typename T>
bool IsFromBeauty(int &particleId, T &pythia);

//__________________________________________________________________________________________________
void SimulateBaryonOverMesonRatio(int nEvents = 50000, int charmedBaryonPDG = 4122, int charmedMesonPDG = 421,
                                  int strangeBaryonPDG = 3122, int strangeMesonPDG = 310, int tune = kCRMode2,
                                  int process = kSoftQCD, float energy = 13000,
                                  std::string oFileNameBase = "pythia_npLc_over_npD0") {
    const int nBinsPerGeV = 1000;

    //__________________________________________________________
    // create and configure pythia generator

    Pythia pythia;

    if (process == kSoftQCD) {
        pythia.readString("SoftQCD:all = on");
    } else if (process == kHardQCD) {
        pythia.readString("HardQCD:hardccbar = on");
        pythia.readString("HardQCD:hardbbbar = on");
    }

    // set tune
    if (tune == kMonash) {
        pythia.readString(Form("Tune:pp = 14"));
    } else if (tune == kCRMode0) {
        pythia.readString(Form("Tune:pp = 14"));
        pythia.readString("ColourReconnection:mode = 1");
        pythia.readString("ColourReconnection:allowDoubleJunRem = off");
        pythia.readString("ColourReconnection:m0 = 2.9");
        pythia.readString("ColourReconnection:allowJunctions = on");
        pythia.readString("ColourReconnection:junctionCorrection = 1.43");
        pythia.readString("ColourReconnection:timeDilationMode = 0");
        pythia.readString("StringPT:sigma = 0.335");
        pythia.readString("StringZ:aLund = 0.36");
        pythia.readString("StringZ:bLund = 0.56");
        pythia.readString("StringFlav:probQQtoQ = 0.078");
        pythia.readString("StringFlav:ProbStoUD = 0.2");
        pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.readString("MultiPartonInteractions:pT0Ref = 2.12");
        pythia.readString("BeamRemnants:remnantMode = 1");
        pythia.readString("BeamRemnants:saturation = 5");
    } else if (tune == kCRMode2) {
        pythia.readString(Form("Tune:pp = 14"));
        pythia.readString("ColourReconnection:mode = 1");
        pythia.readString("ColourReconnection:allowDoubleJunRem = off");
        pythia.readString("ColourReconnection:m0 = 0.3");
        pythia.readString("ColourReconnection:allowJunctions = on");
        pythia.readString("ColourReconnection:junctionCorrection = 1.20");
        pythia.readString("ColourReconnection:timeDilationMode = 2");
        pythia.readString("ColourReconnection:timeDilationPar = 0.18");
        pythia.readString("StringPT:sigma = 0.335");
        pythia.readString("StringZ:aLund = 0.36");
        pythia.readString("StringZ:bLund = 0.56");
        pythia.readString("StringFlav:probQQtoQ = 0.078");
        pythia.readString("StringFlav:ProbStoUD = 0.2");
        pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.readString("MultiPartonInteractions:pT0Ref = 2.15");
        pythia.readString("BeamRemnants:remnantMode = 1");
        pythia.readString("BeamRemnants:saturation = 5");
    } else if (tune == kCRMode3) {
        pythia.readString(Form("Tune:pp = 14"));
        pythia.readString("ColourReconnection:mode = 1");
        pythia.readString("ColourReconnection:allowDoubleJunRem = off");
        pythia.readString("ColourReconnection:m0 = 0.3");
        pythia.readString("ColourReconnection:allowJunctions = on");
        pythia.readString("ColourReconnection:junctionCorrection = 1.15");
        pythia.readString("ColourReconnection:timeDilationMode = 3");
        pythia.readString("ColourReconnection:timeDilationPar = 0.073");
        pythia.readString("StringPT:sigma = 0.335");
        pythia.readString("StringZ:aLund = 0.36");
        pythia.readString("StringZ:bLund = 0.56");
        pythia.readString("StringFlav:probQQtoQ = 0.078");
        pythia.readString("StringFlav:ProbStoUD = 0.2");
        pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.readString("MultiPartonInteractions:pT0Ref = 2.05");
        pythia.readString("BeamRemnants:remnantMode = 1");
        pythia.readString("BeamRemnants:saturation = 5");
    }

    // init
    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed %d", 0));
    pythia.settings.mode("Beams:idA", 2212);
    pythia.settings.mode("Beams:idB", 2212);
    pythia.settings.parm("Beams:eCM", energy);
    pythia.init();

    //__________________________________________________________
    // define outputs
    std::string tuneStr;
    if (tune == kMonash)
        tuneStr = "Monash";
    else if (tune == kCRMode0)
        tuneStr = "CRMode0";
    else if (tune == kCRMode2)
        tuneStr = "CRMode2";
    else if (tune == kCRMode3)
        tuneStr = "CRMode3";

    std::string processStr;
    if (process == kSoftQCD)
        processStr = "SoftQCD";
    else if (process == kHardQCD)
        processStr = "hardQCD";

    auto oFile = new TFile(Form("%s_tune-%s_process-%s.root", oFileNameBase.data(), tuneStr.data(), processStr.data()),
                           "recreate");
    std::map<std::string, TH1F *> hCharmedBaryon;
    std::map<std::string, TH1F *> hCharmedMeson;

    double maxPt = 24.;
    std::string charmedBaryonLatex;
    if (charmedBaryonPDG == 4122) {
        charmedBaryonLatex = "#Lambda_{c}^{+}";
    } else {
        printf("Error: particle %d is not supported. Exit!", charmedBaryonPDG);
        return;
    }

    std::string charmedMesonLatex;
    if (charmedMesonPDG == 421) {
        charmedMesonLatex = "D^{0}";
    } else {
        printf("Error: particle %d is not supported. Exit!", charmedMesonPDG);
        return;
    }

    std::string strangeBaryonLatex;
    if (strangeBaryonPDG == 3122) {
        strangeBaryonLatex = "#Lambda^{0}";
    } else {
        printf("Error: particle %d is not supported. Exit!", strangeBaryonPDG);
        return;
    }

    std::string strangeMesonLatex;
    if (strangeMesonPDG == 310) {
        strangeMesonLatex = "K^{0}_{S}";
    } else {
        printf("Error: particle %d is not supported. Exit!", strangeMesonPDG);
        return;
    }

    // baryon histograms
    int nBins = nBinsPerGeV * maxPt;
    auto title = Form(";#it{p}_{T}(%s) (GeV/#it{c});Counts", charmedBaryonLatex.c_str());
    hCharmedBaryon.insert({"prompt", new TH1F("hCharmedBaryonPrompt", title, nBins, 0., maxPt)});

    title = Form(";#it{p}_{T}(%s) (GeV/#it{c});Counts", charmedBaryonLatex.c_str());
    hCharmedBaryon.insert({"nonprompt", new TH1F("hCharmedBaryonNonPrompt", title, nBins, 0., maxPt)});

    title = Form(";#it{p}_{T}(%s) (GeV/#it{c});Counts", strangeBaryonLatex.c_str());
    TH1F *hStrangeBaryon = new TH1F("hStrangeBaryon", title, nBins, 0., maxPt);

    // meson histograms
    title = Form(";#it{p}_{T}(%s) (GeV/#it{c});Counts", charmedMesonLatex.c_str());
    hCharmedMeson.insert({"prompt", new TH1F("hCharmedMesonPrompt", title, nBins, 0., maxPt)});

    title = Form(";#it{p}_{T}(%s) (GeV/#it{c});Counts", charmedMesonLatex.c_str());
    hCharmedMeson.insert({"nonprompt", new TH1F("hCharmedMesonNonPrompt", title, nBins, 0., maxPt)});

    title = Form(";#it{p}_{T}(%s) (GeV/#it{c});Counts", strangeMesonLatex.c_str());
    TH1F *hStrangeMeson = new TH1F("hStrangeMeson", title, nBins, 0., maxPt);

    //__________________________________________________________
    // perform the simulation

    for (int iEvent = 0; iEvent < nEvents; iEvent++) {
        pythia.next();
        for (int iPart = 3; iPart < pythia.event.size(); iPart++) {
            auto &part = pythia.event[iPart];

            if (abs(part.y()) > 0.5) continue;
            int absPdg = std::abs(part.id());

            if (charmedBaryonPDG == absPdg) {
                if (IsFromBeauty(iPart, pythia))
                    hCharmedBaryon["nonprompt"]->Fill(part.pT());
                else
                    hCharmedBaryon["prompt"]->Fill(part.pT());
            } else if (charmedMesonPDG == absPdg) {
                if (IsFromBeauty(iPart, pythia))
                    hCharmedMeson["nonprompt"]->Fill(part.pT());
                else
                    hCharmedMeson["prompt"]->Fill(part.pT());
            } else if (strangeBaryonPDG == absPdg) {
                double rProd = std::sqrt(part.xProd() * part.xProd() + part.yProd() * part.yProd() +
                                         part.zProd() * part.zProd());  // units: mm
                if (rProd < 1)                                          // exclude weak decays
                    hStrangeBaryon->Fill(part.pT());
            } else if (strangeMesonPDG == absPdg) {
                double rProd = std::sqrt(part.xProd() * part.xProd() + part.yProd() * part.yProd() +
                                         part.zProd() * part.zProd());  // units: mm
                if (rProd < 1)                                          // exclude weak decays
                    hStrangeMeson->Fill(part.pT());
            }
        }
    }
    hCharmedBaryon["prompt"]->Write();
    hCharmedBaryon["nonprompt"]->Write();
    hCharmedMeson["prompt"]->Write();
    hCharmedMeson["nonprompt"]->Write();
    hStrangeBaryon->Write();
    hStrangeMeson->Write();

    oFile->Close();
}

//__________________________________________________________________________________________________
template <typename T>
bool IsFromBeauty(int &particleIndex, T &pythia) {
    std::vector<std::vector<int>> indecesToTest{};
    indecesToTest.push_back({particleIndex});
    auto depth = 0;

    while (indecesToTest[depth].size() > 0) {
        std::vector<int> idsTmp{};
        for (auto &id : indecesToTest[depth]) {
            auto mothers = pythia.event[id].motherList();
            for (auto &mom : mothers) {
                idsTmp.push_back(mom);
                int absPdgMom = std::abs(pythia.event[mom].id());
                if (absPdgMom == 5 ||                   // b quark
                    absPdgMom / 100 == 5 ||             // b mesons
                    absPdgMom / 1000 == 5 ||            // b baryons
                    (absPdgMom - 10000) / 100 == 5 ||   // bbbar resonances
                    (absPdgMom - 20000) / 100 == 5 ||   // bbbar resonances
                    (absPdgMom - 30000) / 100 == 5 ||   // bbbar resonances
                    (absPdgMom - 100000) / 100 == 5 ||  // bbbar resonances
                    (absPdgMom - 200000) / 100 == 5 ||  // bbbar resonances
                    (absPdgMom - 300000) / 100 == 5     // bbbar resonances
                ) {
                    return true;
                }
            }
        }
        indecesToTest.push_back(idsTmp);
        depth++;
    }
    return false;
}
