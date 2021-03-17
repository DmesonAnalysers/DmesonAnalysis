#ifndef ALIANALYSISTASKCHECKDPROTONPAIRGEN_H
#define ALIANALYSISTASKCHECKDPROTONPAIRGEN_H

#include <array>
#include <deque>

#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskCheckDprotonPairGen : public AliAnalysisTaskSE
{

public:
    AliAnalysisTaskCheckDprotonPairGen();
    virtual ~AliAnalysisTaskCheckDprotonPairGen();

    virtual void UserExec(Option_t *option);
    virtual void UserCreateOutputObjects();
    virtual void Terminate(Option_t *option);

    double ComputeKstar(std::array<double, 4> part1, std::array<double, 4> part2);

private:
    AliAnalysisTaskCheckDprotonPairGen(const AliAnalysisTaskCheckDprotonPairGen &source);
    AliAnalysisTaskCheckDprotonPairGen &operator=(const AliAnalysisTaskCheckDprotonPairGen &source);

    TList *fOutput = nullptr;                                                    //!<! list of output histos
    TH1F *fHistNEvents = nullptr;                                                //!<! histo with N of events
    TH1F *fHistDplusProtonPairs[4] = {nullptr, nullptr, nullptr, nullptr};       //!<! histos with D+p pairs vs k*
    TH1F *fHistDplusAntiProtonPairs[4] = {nullptr, nullptr, nullptr, nullptr};   //!<! histos with D+antip pairs vs k*
    TH1F *fHistDplusProtonPairsME[4] = {nullptr, nullptr, nullptr, nullptr};     //!<! histos with D+p pairs vs k* for ME
    TH1F *fHistDplusAntiProtonPairsME[4] = {nullptr, nullptr, nullptr, nullptr}; //!<! histos with D+antip pairs vs k* for ME
    AliAODEvent *fAOD = nullptr;                                                 /// pointer to current event

    // buffers for ME
    std::deque<std::vector<std::array<double, 4> > > fDVecBuffer;                ///buffer for D+
    std::deque<std::vector<std::array<double, 4> > > fpVecBuffer;                ///buffer for p
    std::deque<std::vector<int> > fMotherIdxDVecBuffer;                          ///buffer for D+ mother index
    std::deque<std::vector<int> > fMotherIdxpVecBuffer;                          ///buffer for p mother index
    std::deque<std::vector<int> > fOriginMotherDVecBuffer;                       ///buffer for D+ origin
    std::deque<std::vector<int> > fOriginMotherpVecBuffer;                       ///buffer for p origin
    std::deque<std::vector<int> > fPdgMotherDVecBuffer;                          ///buffer for D mother pdg
    std::deque<std::vector<int> > fPdgMotherpVecBuffer;                          ///buffer for p mother pdg

    /// \cond CLASSIMP
    ClassDef(AliAnalysisTaskCheckDprotonPairGen, 1);
    /// \endcond
};

#endif
