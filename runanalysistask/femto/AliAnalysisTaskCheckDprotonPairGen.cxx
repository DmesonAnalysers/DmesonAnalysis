#include <vector>

#include <TChain.h>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliVertexingHFUtils.h"
#include <TH1F.h>
#include "AliAODInputHandler.h"
#include "AliAODMCHeader.h"
#include "AliLog.h"
#include "AliAnalysisTaskCheckDprotonPairGen.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskCheckDprotonPairGen);
/// \endcond

//______________________________________________________________________________
AliAnalysisTaskCheckDprotonPairGen::AliAnalysisTaskCheckDprotonPairGen() : AliAnalysisTaskSE("DprotonMCChecks")
{
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}

//___________________________________________________________________________
AliAnalysisTaskCheckDprotonPairGen::~AliAnalysisTaskCheckDprotonPairGen()
{
    delete fOutput;
}

//___________________________________________________________________________
void AliAnalysisTaskCheckDprotonPairGen::UserCreateOutputObjects()
{
    /// create output histos

    fOutput = new TList();
    fOutput->SetOwner();
    fOutput->SetName("OutputHistos");

    fHistNEvents = new TH1F("fHistNEvents", "Number of processed events", 3, -0.5, 2.5);
    fHistNEvents->GetXaxis()->SetBinLabel(1, "Read events");
    fHistNEvents->GetXaxis()->SetBinLabel(2, "Selected SE events");
    fHistNEvents->GetXaxis()->SetBinLabel(3, "ME events");
    fHistNEvents->SetMinimum(0);
    fOutput->Add(fHistNEvents);

    fHistDplusProtonPairs[0] = new TH1F("fHistDplusProtonPairsDPrompt", ";#it{k}* (GeV/#it{c});", 1000, 0., 1.);
    fHistDplusProtonPairs[1] = new TH1F("fHistDplusProtonPairsDFromB", ";#it{k}* (GeV/#it{c});", 1000, 0., 1.);
    fHistDplusProtonPairs[2] = new TH1F("fHistDplusProtonPairsFromSameB", ";#it{k}* (GeV/#it{c});", 1000, 0., 1.);
    fHistDplusProtonPairs[3] = new TH1F("fHistDplusProtonPairsDFromDstar", ";#it{k}* (GeV/#it{c});", 1000, 0., 1.);
    fHistDplusAntiProtonPairs[0] = new TH1F("fHistDplusAntiProtonPairsDPrompt", ";#it{k}* (GeV/#it{c});", 1000, 0., 1.);
    fHistDplusAntiProtonPairs[1] = new TH1F("fHistDplusAntiProtonPairsDFromB", ";#it{k}* (GeV/#it{c});", 1000, 0., 1.);
    fHistDplusAntiProtonPairs[2] = new TH1F("fHistDplusAntiProtonPairsFromSameB", ";#it{k}* (GeV/#it{c});", 1000, 0., 1.);
    fHistDplusAntiProtonPairs[3] = new TH1F("fHistDplusAntiProtonPairsDFromDstar", ";#it{k}* (GeV/#it{c});", 1000, 0., 1.);
    fHistDplusProtonPairsME[0] = new TH1F("fHistDplusProtonPairsMEDPrompt", ";#it{k}* (GeV/#it{c});", 1000, 0., 1.);
    fHistDplusProtonPairsME[1] = new TH1F("fHistDplusProtonPairsMEDFromB", ";#it{k}* (GeV/#it{c});", 1000, 0., 1.);
    fHistDplusProtonPairsME[2] = new TH1F("fHistDplusProtonPairsMEFromSameB", ";#it{k}* (GeV/#it{c});", 1000, 0., 1.);
    fHistDplusProtonPairsME[3] = new TH1F("fHistDplusProtonPairsMEDFromDstar", ";#it{k}* (GeV/#it{c});", 1000, 0., 1.);
    fHistDplusAntiProtonPairsME[0] = new TH1F("fHistDplusAntiProtonPairsMEDPrompt", ";#it{k}* (GeV/#it{c});", 1000, 0., 1.);
    fHistDplusAntiProtonPairsME[1] = new TH1F("fHistDplusAntiProtonPairsMEDFromB", ";#it{k}* (GeV/#it{c});", 1000, 0., 1.);
    fHistDplusAntiProtonPairsME[2] = new TH1F("fHistDplusAntiProtonPairsMEFromSameB", ";#it{k}* (GeV/#it{c});", 1000, 0., 1.);
    fHistDplusAntiProtonPairsME[3] = new TH1F("fHistDplusAntiProtonPairsMEDFromDstar", ";#it{k}* (GeV/#it{c});", 1000, 0., 1.);
    for (int iHist = 0; iHist < 4; iHist++)
    {
        fOutput->Add(fHistDplusProtonPairs[iHist]);
        fOutput->Add(fHistDplusAntiProtonPairs[iHist]);
        fOutput->Add(fHistDplusProtonPairsME[iHist]);
        fOutput->Add(fHistDplusAntiProtonPairsME[iHist]);
    }

    PostData(1, fOutput);
}
//______________________________________________________________________________
void AliAnalysisTaskCheckDprotonPairGen::UserExec(Option_t *)
{
    fAOD = dynamic_cast<AliAODEvent *>(InputEvent());
    if (!fAOD)
    {
        AliWarning("Bad event");
        return;
    }

    fHistNEvents->Fill(0);

    TClonesArray *arrayMC = dynamic_cast<TClonesArray *>(fAOD->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
    if (!arrayMC)
    {
        AliWarning("MC particles branch not found");
        return;
    }
    AliAODMCHeader *mcHeader = dynamic_cast<AliAODMCHeader *>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    if (!mcHeader)
    {
        printf("MC header branch not found");
        return;
    }
    if (TMath::Abs(mcHeader->GetVtxZ()) > 10)
        return;

    fHistNEvents->Fill(1);

    // vectors for SE
    std::vector<std::array<double, 4> > DpVec{};
    std::vector<std::array<double, 4> > pVec{};
    std::vector<int> motherIdxDpVec{};
    std::vector<int> motherIdxpVec{};
    std::vector<int> originMotherDpVec{};
    std::vector<int> originMotherpVec{};
    std::vector<int> pdgMotherDpVec{};
    std::vector<int> pdgMotherpVec{};

    int nParticles = arrayMC->GetEntriesFast();
    for (int iPart = 0; iPart < nParticles; iPart++)
    {
        AliAODMCParticle *mcPart = dynamic_cast<AliAODMCParticle *>(arrayMC->At(iPart));
        int pdg = mcPart->PdgCode();
        int absPdg = TMath::Abs(pdg);
        if (absPdg != 411 && absPdg != 2212)
            continue;

        double pt = mcPart->Pt();
        double eta = mcPart->Eta();
        double rapid = mcPart->Y();

        bool searchUpToQuark = true ? absPdg == 411 : false;
        int origin = AliVertexingHFUtils::CheckOrigin(arrayMC, mcPart, searchUpToQuark);
        if (origin < 4 && absPdg == 411)
            continue;

        std::array<double, 4> part = {static_cast<double>(mcPart->GetPdgCode()), mcPart->Px(), mcPart->Py(), mcPart->Pz()};
        if(absPdg == 411)
        {
            DpVec.push_back(part);
            originMotherDpVec.push_back(origin);
        }
        else
        {
            pVec.push_back(part);
            originMotherpVec.push_back(origin);
        }

        if(origin == 5 || (origin == 4 && absPdg == 411))
        {
            bool isNotFromDstar = true;
            int mother = mcPart->GetMother();
            while (mother >= 0)
            {
                AliAODMCParticle *mcGranma = dynamic_cast<AliAODMCParticle *>(arrayMC->At(mother));
                if (mcGranma)
                {
                    int pdgGranma = mcGranma->GetPdgCode();
                    int abspdgGranma = TMath::Abs(pdgGranma);
                    if(origin == 5)
                    {
                        if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000))
                        {
                            if(absPdg == 411)
                            {
                                motherIdxDpVec.push_back(mother);
                                pdgMotherDpVec.push_back(pdgGranma);
                            }
                            else
                            {
                                motherIdxpVec.push_back(mother);
                                pdgMotherpVec.push_back(pdgGranma);
                            }
                            break; // stop at first one
                        }
                    }
                    else if (origin == 4 && absPdg == 411)
                    {
                        if(abspdgGranma == 413)
                        {
                            motherIdxDpVec.push_back(mother);
                            pdgMotherDpVec.push_back(pdgGranma);
                            isNotFromDstar = false;
                            break; // stop at first one
                        }
                    }
                    mother = mcGranma->GetMother();
                }
            }
            if (origin == 4 && absPdg == 411 && isNotFromDstar)
            {
                motherIdxDpVec.push_back(-1);
                pdgMotherDpVec.push_back(0);
            }
        }
        else            
        {
            if(absPdg == 411)
            {
                motherIdxDpVec.push_back(-1);
                pdgMotherDpVec.push_back(0);
            }
            else
            {
                motherIdxpVec.push_back(-1);
                pdgMotherpVec.push_back(0);
            }
        }
    }

    // fill buffers
    fDVecBuffer.push_back(DpVec);
    fpVecBuffer.push_back(pVec);
    fMotherIdxDVecBuffer.push_back(motherIdxDpVec);
    fMotherIdxpVecBuffer.push_back(motherIdxpVec);
    fOriginMotherDVecBuffer.push_back(originMotherDpVec);
    fOriginMotherpVecBuffer.push_back(originMotherpVec);
    fPdgMotherDVecBuffer.push_back(pdgMotherDpVec);
    fPdgMotherpVecBuffer.push_back(pdgMotherpVec);
    if (fDVecBuffer.size() > 10) { //buffer full, let's kill the first entry
        fDVecBuffer.pop_front();
        fpVecBuffer.pop_front();
        fMotherIdxDVecBuffer.pop_front();
        fMotherIdxpVecBuffer.pop_front();
        fOriginMotherDVecBuffer.pop_front();
        fOriginMotherpVecBuffer.pop_front();
        fPdgMotherDVecBuffer.pop_front();
        fPdgMotherpVecBuffer.pop_front();
    } 

    // loop on same event
    for(size_t iD=0; iD<DpVec.size(); iD++)
    {
        for(size_t ip=0; ip<pVec.size(); ip++)
        {
            double kStar = ComputeKstar(DpVec[iD], pVec[ip]);
            if((DpVec[iD][0] > 0 && pVec[iD][ip] > 0) || 
               (DpVec[iD][0] < 0 && pVec[iD][ip] < 0))
            {
                if(originMotherDpVec[iD] == 4)
                {
                    fHistDplusProtonPairs[0]->Fill(kStar);
                    if(abs(pdgMotherDpVec[iD]) == 413)
                        fHistDplusProtonPairs[3]->Fill(kStar);
                }
                else
                {
                    if(originMotherpVec[ip] == 5 && motherIdxpVec[ip] == motherIdxpVec[iD])
                        fHistDplusProtonPairs[2]->Fill(kStar);
                    else
                        fHistDplusProtonPairs[1]->Fill(kStar);
                }
            }
            else
            {
                if(originMotherDpVec[iD] == 4)
                {
                    fHistDplusAntiProtonPairs[0]->Fill(kStar);
                    if(abs(pdgMotherDpVec[iD]) == 413)
                        fHistDplusAntiProtonPairs[3]->Fill(kStar);
                }
                else
                {
                    if(originMotherpVec[ip] == 5 && motherIdxpVec[ip] == motherIdxpVec[iD])
                        fHistDplusAntiProtonPairs[2]->Fill(kStar);
                    else
                        fHistDplusAntiProtonPairs[1]->Fill(kStar);
                }
            }
        }
    }

    // loop on mixed event
    if(fpVecBuffer.size() < 10) // to avoid repetitions
        return;

    // mix last D-meson vector with all previous proton vectors
    for(size_t iME=0; iME<fpVecBuffer.size(); iME++)
    {
        fHistNEvents->Fill(2);
        for(size_t iD=0; iD<fDVecBuffer[9].size(); iD++)
        {
            for(size_t ip=0; ip<fpVecBuffer[iME].size(); ip++)
            {
                double kStar = ComputeKstar(fDVecBuffer[9][iD], fpVecBuffer[iME][ip]);
                if((fDVecBuffer[9][iD][0] > 0 && fpVecBuffer[iME][ip][0] > 0) || 
                   (fDVecBuffer[9][iD][0] < 0 && fpVecBuffer[iME][ip][0] < 0))
                {
                    if(fOriginMotherDVecBuffer[9][iD] == 4)
                    {
                        fHistDplusProtonPairsME[0]->Fill(kStar);
                        if(abs(fPdgMotherDVecBuffer[9][iD]) == 413)
                            fHistDplusProtonPairsME[3]->Fill(kStar);
                    }
                    else
                        fHistDplusProtonPairsME[1]->Fill(kStar);
                }
                else
                {
                    if(fOriginMotherDVecBuffer[9][iD] == 4)
                    {
                        fHistDplusAntiProtonPairsME[0]->Fill(kStar);
                        if(abs(fPdgMotherDVecBuffer[9][iD]) == 413)
                            fHistDplusAntiProtonPairsME[3]->Fill(kStar);
                    }
                    else
                        fHistDplusAntiProtonPairsME[1]->Fill(kStar);
                }
            }
        }
    }

    PostData(1, fOutput);
}

//______________________________________________________________________________
void AliAnalysisTaskCheckDprotonPairGen::Terminate(Option_t * /*option*/)
{
    /// Terminate analysis
    fOutput = dynamic_cast<TList *>(GetOutputData(1));
    if (!fOutput)
    {
        AliWarning("fOutput not available");
        return;
    }

    return;
}

double AliAnalysisTaskCheckDprotonPairGen::ComputeKstar(std::array<double, 4> part1, std::array<double, 4> part2)
{
    TLorentzVector SPtrack, TPProng, trackSum, SPtrackCMS, TPProngCMS;
    SPtrack.SetXYZM(part1[1], part1[2], part1[3],
                    TDatabasePDG::Instance()->GetParticle(abs(static_cast<int>(part1[0])))->Mass());
    TPProng.SetXYZM(part2[1], part2[2], part2[3],
                    TDatabasePDG::Instance()->GetParticle(abs(static_cast<int>(part2[0])))->Mass());
    trackSum = SPtrack + TPProng;

    double beta = trackSum.Beta();
    double betax = beta * cos(trackSum.Phi()) * sin(trackSum.Theta());
    double betay = beta * sin(trackSum.Phi()) * sin(trackSum.Theta());
    double betaz = beta * cos(trackSum.Theta());

    SPtrackCMS = SPtrack;
    TPProngCMS = TPProng;

    SPtrackCMS.Boost(-betax, -betay, -betaz);
    TPProngCMS.Boost(-betax, -betay, -betaz);

    TLorentzVector trackRelK;

    trackRelK = SPtrackCMS - TPProngCMS;
    double kStar = 0.5 * trackRelK.P();
    return kStar;
}