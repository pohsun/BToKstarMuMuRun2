#ifndef MUON_BRANCH_PLANTER_H
#define MUON_BRANCh_PLANTER_H

#include "BphAna/BToKstarMuMuRun2/interface/AbstractBranchPlanter.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#define MAX_MUONS 64

class MuonBranchPlanter : public AbstractBranchPlanter{
    public:
        MuonBranchPlanter();
        ~MuonBranchPlanter();

        // cuttage, plant, reap

        // setters, getters, and adders
        bool setBranchAddress(TTree*,const char prefix[]="MuonInfo");

        // services
        virtual const std::string getClassName()        {return "MuonBranchPlanter";}
        virtual const std::string getMotherClassName()  {return "AbstractBranchPlanter";}
        virtual void print();

        // Branch variables defined as PUBLIC for convenient accessment.
        unsigned int nMu;
        std::vector<int>        charge;
        std::vector<double>     pt;
        std::vector<double>     eta;
        std::vector<double>     phi;
            // Trajectory qualities
        std::vector<int>        nPixHit;
        std::vector<int>        nTrkHit;
        std::vector<int>        nPixLayer;
        std::vector<int>        nTrkLayer;
        std::vector<double>     dxy;
        std::vector<double>     dz;
        std::vector<double>     normChi2;
        std::vector<double>     dcaBS;
        std::vector<double>     dcaBSErr;

        // Following branches to be filled by other Planter
            // Trigger and MC info
        std::vector<bool>       isL1Matched;
        std::vector<bool>       isMCMatched;

    private:

        // Variables and objects not in branch


        // cuttage
        bool isHandlesOk();
        void regBranches(const char prefix[]="MuonInfo");
        //void regHists(const char prefix[]="MuonInfo");

        // plant
        void clearVariables();
        void buildBranches(const edm::Event&, const edm::EventSetup&);
        void l1Match();//no L3 muon in the skim
        void calcDcaBS(reco::TrackRef,double&,double&);

};

#endif
