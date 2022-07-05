#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TDatabasePDG.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THn.h>
#include <TMath.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TSpline.h>
#include <TSystem.h>

#include <array>
#include <deque>
#include <map>
#include <string>
#include <vector>

#include "Math/GenVector/Boost.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Pythia8/Pythia.h"

#endif

using namespace Pythia8;

namespace pythiasettings {
enum tunes { kMonash = 0, kCRMode0, kCRMode2, kCRMode3 };

enum processes { kSoftQCD = 0, kHardQCD };

std::array<int, 2> baryonPDG{411, 421};  // D+, D0
std::array<int, 2> mesonPDG{413, 423};   // D*+, D*0

// Lc

}  // namespace

//__________________________________________________________________________________________________
template <typename T>
bool CheckDauAcc(T &ptDauArray, T &etaDauArray, T &eDauArray, T &pdgDauArray, double ptMin, double etaMin,
                 double etaMax);
template <typename T, typename T2>
bool IsFromBeauty(T &mothers, T2 &pythia);

//__________________________________________________________________________________________________
void SimulateBaryonOverMesonRatio(int nEvents = 1000, int baryonPDG = 4122, int mesonPDG = 421, int tune = pythiasettings::kCRMode2,
                                  int process = pythiasettings::kSoftQCD, float energy = 13000, int seed = 42,
                                  std::string outFileNameRoot = "pythia_npLcnpD0", int threadlabel = 0) {
    //__________________________________________________________
    // create and configure pythia generator
    ROOT::EnableImplicitMT(16);

    outFileNameRoot = std::string(Form("%s_%d_t%d.root", outFileNameRoot.c_str(), nEvents, threadlabel));
    Pythia pythia;

    if (process == pythiasettings::kSoftQCD) {
        pythia.readString("SoftQCD:all = on");
    } else if (process == pythiasettings::kHardQCD) {
        pythia.readString("HardQCD:hardccbar = on");
        pythia.readString("HardQCD:hardbbbar = on");
    }

    // set tune
    if (tune == pythiasettings::kMonash) {
        pythia.readString(Form("Tune:pp = 14"));
    } else if (tune == pythiasettings::kCRMode0) {
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
    } else if (tune == pythiasettings::kCRMode2) {
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
    } else if (tune == pythiasettings::kCRMode3) {
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

    // keep only interesting decays, to be reweighted a posteriori
    // baryon decay
    // if (baryonPDG == 4122) {
    //     pythia.readString("4122:onMode = off");
    //     pythia.readString("4122:onIfMatch = 2212 321 211");
    //     pythia.readString("4122:onIfMatch = 2212 310"); // todo add pKpi
    // }

    // // meson decay
    // if (mesonPDG == 421) {
    //     pythia.readString("421:onMode = off");
    //     pythia.readString("421:onIfMatch = 321 321 211");
    // }

    // pythia.readString("4122:onMode = off");
    // pythia.readString("411:onIfMatch = 211 211 321");
    // pythia.readString("413:onIfMatch = 211 421");
    // pythia.readString("423:onIfMatch = 22 421");  // for simplicity, let's consider only D0* --> D0 gamma

    // init
    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed %d", seed));
    pythia.settings.mode("Beams:idA", 2212);
    pythia.settings.mode("Beams:idB", 2212);
    pythia.settings.parm("Beams:eCM", energy);
    pythia.init();

    //__________________________________________________________
    // define outputs
    auto oFile = new TFile(outFileNameRoot.data(), "recreate");
    std::map<std::string, TH1F *> hBaryon;
    std::map<std::string, TH1F *> hMeson;

    double maxPt = 24.;
    auto hCiccio = new TH1F("hCiccio", "Counts;#it{p}_{T}^{Lc} (GeV/#it{c});", 100, 0., maxPt);
    hBaryon.insert({"prompt", new TH1F("hBaryonPrompt", "Counts;#it{p}_{T}^{Lc} (GeV/#it{c});", 100, 0., maxPt)});
    hBaryon.insert({"nonprompt", new TH1F("hBaryonNonPrompt", "Counts;#it{p}_{T}^{Lc} (GeV/#it{c});", 100, 0., maxPt)});

    hMeson.insert({"prompt", new TH1F("hMesonPrompt", "Counts;#it{p}_{T}^{Lc} (GeV/#it{c});", 100, 0., maxPt)});
    hMeson.insert({"nonprompt", new TH1F("hMesonNonPrompt", "Counts;#it{p}_{T}^{Lc} (GeV/#it{c});", 100, 0., maxPt)});

    int nBins[4] = {80, 80, 400, 2};
    double mins[4] = {-4., -4., 0., 0.5};
    double maxs[4] = {4., 4., 2., 2.5};

    //__________________________________________________________
    // perform the simulation
    std::vector<ROOT::Math::PxPyPzMVector> partDmeson{};
    std::vector<ROOT::Math::PxPyPzMVector> partDstar{};
    std::deque<std::vector<ROOT::Math::PxPyPzMVector>> partBufferDmeson{};
    std::deque<std::vector<ROOT::Math::PxPyPzMVector>> partBufferDstar{};
    std::vector<int> pdgDmeson{};
    std::vector<int> pdgDstar{};
    std::vector<float> effDmeson{};
    std::vector<float> effDstar{};
    std::vector<int> motherDmeson{};
    std::vector<int> idxDstar{};
    std::vector<float> yDmeson{};
    std::vector<float> yDstar{};
    std::deque<std::vector<int>> pdgBufferDmeson{};
    std::deque<std::vector<int>> pdgBufferDstar{};
    std::deque<std::vector<float>> effBufferDmeson{};
    std::deque<std::vector<float>> effBufferDstar{};
    std::vector<bool> isMesonInAcceptance{};
    std::vector<bool> isBaryonInAcceptance{};

    for (int iEvent = 0; iEvent < nEvents; iEvent++) {
        // continue;
        pythia.next();
        for (int iPart = 3; iPart < pythia.event.size(); iPart++) {
            int pdg = pythia.event[iPart].id();
            int absPdg = std::abs(pdg);

            bool isBaryonInEvent = baryonPDG == absPdg;
            bool isMesonInEvent = mesonPDG == absPdg;
            // bool isBaryonInEvent = std::find(baryonPDG.begin(), baryonPDG.end(), absPdg) != baryonPDG.end();
            // bool isMesonInEvent = std::find(mesonPDG.begin(), mesonPDG.end(), absPdg) != mesonPDG.end();

            if (!isMesonInEvent && !isBaryonInEvent) continue;
            // if (isBaryonInEvent && pythia.event[iPart].pT() < 1) continue;

            std::vector<int> mothers = pythia.event[iPart].motherList();

            // todo: check decay kinem of daughters
            auto dauList = pythia.event[iPart].daughterList();
            // std::vector<double> ptDau{}, etaDau{}, pdgDau{}, eDau{};
            // for (auto &dau : dauList) {
            //     auto absPdgDau = std::abs(pythia.event[dau].id());
            //     if ((absPdg != 423 && (absPdgDau == 211 || absPdgDau == 321)) ||
            //         (absPdg == 423 && (absPdgDau == 211 || absPdgDau == 321 || absPdgDau == 22))) {
            //         ptDau.push_back(std::sqrt(pythia.event[dau].px() * pythia.event[dau].px() +
            //                                   pythia.event[dau].py() * pythia.event[dau].py() +
            //                                   pythia.event[dau].pz() * pythia.event[dau].pz()));
            //         etaDau.push_back(pythia.event[dau].eta());
            //         pdgDau.push_back(std::abs(pythia.event[dau].id()));
            //         eDau.push_back(pythia.event[dau].e());
            //     }
            // }
            // if (etaDau.size() == 0) continue;

            // todo: check  eta
            // double etaMaxDau = 0.;
            // for (auto &eta : etaDau) {
            //     if (std::abs(eta) > std::abs(etaMaxDau)) etaMaxDau = eta;
            // }

            if (isMesonInEvent) {
                bool isFromB = IsFromBeauty(mothers, pythia);
                // if (!CheckDauAcc(ptDau, etaDau, eDau, pdgDau, 0.05, -4., 4.)) continue;
                // isMesonInAcceptance.push_back(CheckDauAcc(ptDau, etaDau, eDau, pdgDau, 0.05, -0.8, 0.8));

                auto ptMeson = pythia.event[dauList[0]].pT();
                // printf("ptmeson: %f    ---    ", ptMeson);
                // // auto yMeson = abs(pythia.event[dauList[0]].y()); // todo: need?
                // // if (ptMeson < 1 || yMeson > 4) continue; // todo: need??
                // for (auto dau : dauList) {
                //     printf("%d ", pythia.event[dau].id());
                // }
                // printf("\n");
                isFromB ? hMeson["nonprompt"]->Fill(ptMeson) : hMeson["prompt"]->Fill(ptMeson);
                // hCiccio->Fill(ptMeson);
                // break;

            } else if (isBaryonInEvent) {
                bool isFromB = IsFromBeauty(mothers, pythia);

                auto ptBaryon = pythia.event[dauList[0]].pT();
                // printf("ptbaryon: %f    ---    ", ptBaryon);
                // //     // if (!CheckDauAcc(ptDau, etaDau, eDau, pdgDau, 0.05, -4., 4.)) continue;
                // //     // isMesonInAcceptance.push_back(CheckDauAcc(ptDau, etaDau, eDau, pdgDau, 0.05, -0.8, 0.8));
                // for (auto dau : dauList) {
                //     printf("%d ", pythia.event[dau].id());
                // }
                // printf("\n");
                //     // printf("ptBaryon: %f", ptBaryon);
                //     // auto yMeson = abs(pythia.event[dauList[0]].y()); // todo: need?
                //     // if (ptMeson < 1 || yMeson > 4) continue; // todo: need??

                isFromB ? hBaryon["nonprompt"]->Fill(ptBaryon) : hBaryon["prompt"]->Fill(ptBaryon);
            }
            // break;
        }

        // partDstar.clear();
        // pdgDstar.clear();
        // idxDstar.clear();
        // effDstar.clear();
        // yDstar.clear();
        // isMesonInAcceptance.clear();

        // partDmeson.clear();
        // pdgDmeson.clear();
        // motherDmeson.clear();
        // effDmeson.clear();
        // yDmeson.clear();
        // isBaryonInAcceptance.clear();
    }

    // return;
    // save root output file
    // for (auto &pdgDstar : mesonPDG) {
    //     for (auto &pdgDmeson : baryonPDG) {
    hBaryon["prompt"]->Write();
    hBaryon["nonprompt"]->Write();
    hMeson["prompt"]->Write();
    hMeson["nonprompt"]->Write();
    // //     }
    // hCiccio->Write();
    // }

    std::map<std::string, TH1F *> hBaryonOverMeson;
    hBaryonOverMeson.insert({"prompt", (TH1F *)hBaryon["prompt"]->Clone()});
    hBaryonOverMeson["prompt"]->SetName("hBaryonOverMesonPrompt");
    hBaryonOverMeson["prompt"]->Sumw2();
    hBaryonOverMeson["prompt"]->Divide(hMeson["prompt"]);

    hBaryonOverMeson.insert({"nonprompt", (TH1F *)hBaryon["nonprompt"]->Clone()});
    hBaryonOverMeson["nonprompt"]->SetName("hBaryonOverMesonNonPrompt");
    hBaryonOverMeson["nonprompt"]->Sumw2();
    hBaryonOverMeson["nonprompt"]->Divide(hMeson["nonprompt"]);


    hBaryonOverMeson["prompt"]->Write();
    hBaryonOverMeson["nonprompt"]->Write();

    oFile->Close();
}

//__________________________________________________________________________________________________
template <typename T>
bool CheckDauAcc(T &ptDauArray, T &etaDauArray, T &eDauArray, T &pdgDauArray, double ptMin, double etaMin,
                 double etaMax) {
    for (size_t iDau = 0; iDau < ptDauArray.size(); iDau++) {
        if (pdgDauArray[iDau] != 22)  // no photons
        {
            if (ptDauArray[iDau] < ptMin || etaDauArray[iDau] < etaMin ||
                etaDauArray[iDau] > etaMax)  // pT>50 MeV and |eta|<4
                return false;
        } else  // photons
        {
            if (eDauArray[iDau] < 0.4 || etaDauArray[iDau] > 4. ||
                etaDauArray[iDau] < -1.6)  // E>400 MeV and -1.6<eta<4
                return false;
        }
    }

    return true;
}

//__________________________________________________________________________________________________
template <typename T, typename T2>
bool IsFromBeauty(T &mothers, T2 &pythia) {
    for (auto &mom : mothers) {
        int absPdgMom = std::abs(pythia.event[mom].id());
        if (absPdgMom == 5 || absPdgMom / 100 == 5 || absPdgMom / 1000 == 5 || (absPdgMom - 10000) / 100 == 5 ||
            (absPdgMom - 20000) / 100 == 5 || (absPdgMom - 30000) / 100 == 5 || (absPdgMom - 100000) / 100 == 5 ||
            (absPdgMom - 200000) / 100 == 5 || (absPdgMom - 300000) / 100 == 5) {  // remove beauty feed-down
            return true;
        }
    }
    return false;
}