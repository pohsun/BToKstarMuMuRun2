#ifndef COMMON_BDECAY_SELECT_H
#define COMMON_BDECAY_SELECT_H

#define A_LARGE_NUMBER 9999
#include "BPHRecoSelect.h"
#include "BPHMomentumSelect.h"
#include "BPHVertexSelect.h"

class MuonPtSelect: public BPHRecoSelect {
    public:
        MuonPtSelect(double minPt=0.7, double maxPt=A_LARGE_NUMBER):
            minPtCut(minPt),
            maxPtCut(maxPt){}
        virtual bool accept(const reco::Candidate& cand) const {
            const pat::Muon* p = reinterpret_cast<const pat::Muon*>(&cand);
            if (p == 0) return false;
            return (p->p4().pt() > minPtCut && p->p4().pt() < maxPtCut);
        }
    private:
        double minPtCut, maxPtCut;
}

class MuonEtaSelect: public BPHRecoSelect {
    public:
        MuonEtaSelect(double minEta=0, double maxEta=2.4):
            minEtaCut(minEta),
            maxEtaCut(maxEta){}
        virtual bool accept(const reco::Candidate& cand) const {
            const pat::Muon* p = reinterpret_cast<const pat::Muon*>(&cand);
            if (p == 0) return false;
            return (fabs(p->p4().eta()) > minEtaCut && fabs(p->p4().eta()) < maxEtaCut);
        }
    private:
        double minEtaCut, maxEtaCut;
}

class TrackPtSelect: public BPHRecoSelect {
    public:
        TrackPtSelect(double minPt=0.7, double maxPt=A_LARGE_NUMBER):
            minPtCut(minPt),
            maxPtCut(maxPt){}
        virtual bool accept(const reco::Candidate& cand) const {
            return (cand.p4().pt() > minPtCut && cand.p4().pt() < maxPtCut);
        }
    private:
        double minPtCut, maxPtCut;
}

class TrackEtaSelect: public BPHRecoSelect {
    public:
        TrackEtaSelect(double minEta=0, double maxEta=2.4):
            minEtaCut(minEta),
            maxEtaCut(maxEta){}
        virtual bool accept(const reco::Candidate& cand) const {
            return (fabs(cand.p4().eta()) > minEtaCut && fabs(cand.p4().eta()) < maxEtaCut);
        }
    private:
        double minEtaCut, maxEtaCut;
}

class MassSelect: public BPHMomentumSelect {
    public:
        MassSelect(double minMass, double maxMass=A_LARGE_NUMBER):
            minMassCut(minMass),
            maxMassCut(maxMass){}
        virtual bool accept(const BPHDecayMomentum& cand) const {
            double mass = cand.composite().mass();
            return ( (mass > minMassCut) && (mass < maxMassCut) );
        }
    private:
        double minMassCut, maxMassCut;
}

class VtxProbSelect: public BPHVertexSelect {
    public:
        VtxProbSelect(double minProb):
            minProbCut(minProb){}
        virtual bool accept(const BPHDecayVertex& cand) const {
            const reco::Vertex& v = cand.vertex();
            if (v.isFake())     return false;
            if (!v.isValid())   return false;
            return ( TMath::Prob( v.chi2(), lround(v.ndof())) > minProbCut );

        }
    private:
        double minProbCut;
}

#endif
