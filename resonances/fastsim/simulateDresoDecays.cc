R__LOAD_LIBRARY(EvtGen)
R__ADD_INCLUDE_PATH($EVTGEN_ROOT/include)

#include <array>
#include <string>
#include <vector>

#include <TDatabasePDG.h>
#include <TFile.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TF1.h>

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"

#include "Pythia8/Pythia.h"
#include "EvtGen/EvtGen.hh"
#include "Pythia8Plugins/EvtGen.h"

using namespace Pythia8;

namespace {
enum decayer { kPythia8 = 0, kEvtGen };
}

//__________________________________________________________________________________________________
void SimulateDresoDecays(int nEvents = 10000, int pdgReso = 435,
                         int decayer = kEvtGen, int seed = 42, std::string trigger = "MB");
Particle GetParticle(int pdgCode = 435, TF1* fPtShape = nullptr);
TF1* GetPtShape(std::string inFileName = "HEPDataDvsMult.root", std::string trigger = "MB");
double PowLaw(double* pt, double* pars);

//__________________________________________________________________________________________________
void simulateDresoDecays(int nEvents, int pdgReso, int decayer, int seed, std::string trigger) {
  //__________________________________________________________
  // create and configure pythia generator

  Pythia pythia;
  pythia.readString("HardQCD:hardccbar = on");
  pythia.readString("HardQCD:hardbbbar = on");

  EvtGenDecays *evtgen = nullptr;
  std::string decName = "";
  if (decayer == kPythia8) {
    decName = "Pythia8";
    // keep only interesting decays
    pythia.readString("435:onMode = off");
    pythia.readString("435:onIfMatch = 311 411");
    pythia.readString("10433:onMode = off");
    pythia.readString("10433:onIfMatch = 311 413");
    // switch off D-meson decays
    pythia.readString("421:onMode = off");
    pythia.readString("413:onMode = off");
    pythia.readString("411:onMode = off");
    pythia.readString("310:onMode = off");
    pythia.readString("311:onMode = off");
    pythia.readString("311:onIfMatch = 310");
  } else if (decayer == kEvtGen) {
    // keep only interesting decays
    decName = "EvtGen";
    evtgen = new EvtGenDecays(&pythia, "./DRESODECAYS.DEC", "./evt.pdl");
  } else {
    std::cerr << "ERROR: Invalid decayer selected! Exit." << std::endl;
    return;
  }

  // init
  pythia.readString("Random:setSeed = on");
  pythia.readString(Form("Random:seed %d", seed));
  pythia.init();
  gRandom->SetSeed(seed);

  //__________________________________________________________
  // define outputs
  std::string outFileNameRoot = Form("Decays_%d_%s_decayer_%s.root", pdgReso, decName.data(), trigger.data());
  auto tupleDreso = new TNtuple(
      "tupleDreso", "tupleDreso",
      "mass_reso:pt_reso:y_reso:phi_reso:pt_d:y_d:phi_d:pt_v0:y_v0:phi_v0");

  //__________________________________________________________
  // define reso mass limit
  double massLimit = -1.;
  switch(pdgReso) {
    case 435:
    {
      massLimit = TDatabasePDG::Instance()->GetParticle(411)->Mass() + TDatabasePDG::Instance()->GetParticle(310)->Mass();
      break;
    }
    case 10433:
    {
      massLimit = TDatabasePDG::Instance()->GetParticle(413)->Mass() + TDatabasePDG::Instance()->GetParticle(310)->Mass();
      break;
    }
    default:
    {
      std::cerr << "ERROR: PDG code of resonance " << pdgReso << " not implemented! Exit." << std::endl;
      break;
      return;
    }
  }
  // define pT shape
  auto fPtShape = GetPtShape("HEPDataDvsMult.root", trigger);

  //__________________________________________________________
  // perform the simulation
  for (auto iEvent{1}; iEvent < nEvents; ++iEvent) {

    // reset pythia event and put resonance only
    pythia.event.reset();
    pythia.event.append(GetParticle(pdgReso, fPtShape));
    if (pythia.event[1].m() < massLimit) {
      continue;
    }
    int idPart = pythia.event[1].id();
    pythia.particleData.mayDecay(idPart, true);

    if (evtgen) {
      evtgen->decay();
    } else {
      pythia.moreDecays();
    }

    std::array<float, 4> arrDreso{};
    std::array<float, 3> arrDdau{};
    std::array<float, 3> arrV0dau{};
    bool isResoFilled = false, isDmesonFilled = false, isV0Filled = false;
    for (auto iPart{1}; iPart < pythia.event.size(); ++iPart) {
      int pdg = pythia.event[iPart].id();
      int absPdg = std::abs(pdg);

      float mass = pythia.event[iPart].m();
      float pT = pythia.event[iPart].pT();
      float phi = pythia.event[iPart].phi();
      float y = pythia.event[iPart].y();

      if (!isResoFilled && absPdg == pdgReso) {
        arrDreso = std::array<float, 4>{mass, pT, y, phi};
        isResoFilled = true;
      }
      else if (pdgReso == 435) {
        if (!isDmesonFilled && absPdg == 411) {
          arrDdau = std::array<float, 3>{pT, y, phi};
          isDmesonFilled = true;
        }
        else if (!isV0Filled && absPdg == 310) {
          arrV0dau = std::array<float, 3>{pT, y, phi};
          isV0Filled = true;
        }
      }
      else if (pdgReso == 10433) {
        if (!isDmesonFilled && absPdg == 413) {
          arrDdau = std::array<float, 3>{pT, y, phi};
          isDmesonFilled = true;
        }
        else if (!isV0Filled && absPdg == 310) {
          arrV0dau = std::array<float, 3>{pT, y, phi};
          isV0Filled = true;
        }
      }
    }

    if (!isResoFilled || !isDmesonFilled || !isV0Filled) {
      continue;
    }

    std::array<float, 10> array4tuple{};
    std::copy(arrDreso.begin(), arrDreso.end(), array4tuple.begin());
    std::copy(arrDdau.begin(), arrDdau.end(), array4tuple.begin() + arrDreso.size());
    std::copy(arrV0dau.begin(), arrV0dau.end(), array4tuple.begin() + arrDreso.size() + arrDdau.size());
    tupleDreso->Fill(array4tuple.data());
  }

  // save root output file
  TFile outFile(outFileNameRoot.data(), "recreate");
  fPtShape->Write();
  tupleDreso->Write();
  outFile.Close();
}

//__________________________________________________________
Particle GetParticle(int pdg, TF1* fPtShape)
{

  double mass = gRandom->BreitWigner(TDatabasePDG::Instance()->GetParticle(pdg)->Mass(), TDatabasePDG::Instance()->GetParticle(pdg)->Width());
  double phi = gRandom->Uniform(2 * TMath::Pi());
  double pt = fPtShape ? fPtShape->GetRandom() : gRandom->Uniform(50.);
  double mt = TMath::Sqrt(mass * mass + pt * pt);
  double y = -1.5 + gRandom->Uniform(3.);

  auto fourMom = ROOT::Math::PxPyPzMVector(pt * TMath::Cos(phi), pt * TMath::Sin(phi), mt * TMath::SinH(y), mass);

  Particle part;
  part.id(pdg);
  part.status(11);
  part.xProd(0.);
  part.yProd(0.);
  part.zProd(0.);
  part.tProd(0.);
  part.m(fourMom.M());
  part.e(fourMom.E());
  part.px(fourMom.Px());
  part.py(fourMom.Py());
  part.pz(fourMom.Pz());

  return part;
}

//__________________________________________________________
TF1* GetPtShape(std::string inFileName, std::string trigger) {

  auto inFile = TFile::Open(inFileName.data());
  TGraphErrors* gPtShape = nullptr;
  if (trigger == "MB") {
    gPtShape = (TGraphErrors*)inFile->Get("Table sup3b/Graph1D_y1");
  }
  else if (trigger == "HM") {
    gPtShape = (TGraphErrors*)inFile->Get("Table sup3b/Graph1D_y4");
  }
  else {
    std::cerr << "WARNING: Trigger " << trigger << " not supported for pT shape! Return nullptr." << std::endl;
    return nullptr;
  }
  inFile->Close();

  TF1* fPtShape = new TF1("fPtShape", PowLaw, 0., 24., 4);
  fPtShape->SetName("fPtShape");
  fPtShape->SetParameters(gPtShape->Integral(), 2.6, 2.8, 2.);
  gPtShape->Fit(fPtShape, "IME0Q", "", 2., 24.);

  return fPtShape;
}

//__________________________________________________________
double PowLaw(double *pt, double *pars)
{
    return pars[0] * pt[0] / TMath::Power((1 + TMath::Power(pt[0] / pars[1], pars[3])), pars[2]);
}
