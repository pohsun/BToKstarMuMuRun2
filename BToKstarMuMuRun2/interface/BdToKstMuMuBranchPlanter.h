#ifndef BD_BRANCH_PLANTER_H
#define BD_BRANCH_PLANTER_H

#include "BphAna/BToKstarMuMuRun2/interface/AbstractBranchPlanter.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"

#define MAX_MUONS 64

class BdToKstMuMuBranchPlanter : public AbstractBranchPlanter{
    public:
        BdToKstMuMuBranchPlanter();
        ~BdToKstMuMuBranchPlanter();

        // cuttage, plant, reap

        // setters, getters, and adders
        bool setBranchAddress(TTree*,const char prefix[]="BdInfo");

        // services
        virtual const std::string getClassName()        {return "BdToKstMuMuBranchPlanter";}
        virtual const std::string getMotherClassName()  {return "AbstractBranchPlanter";}
        virtual void print();

        // Branch variables defined as PUBLIC for convenient accessment.
        //unsigned int nKs;
        //std::vector<double>     ptKs;
        //std::vector<double>     etaKs;
        //std::vector<double>     phiKs;
        //std::vector<double>     xKs;
        //std::vector<double>     yKs;
        //std::vector<double>     zKs;
        //std::vector<double>     normChi2Ks;

        unsigned int nKst;
        std::vector<double>     ptKst;
        std::vector<double>     etaKst;
        std::vector<double>     phiKst;
        std::vector<double>     xKst;
        std::vector<double>     yKst;
        std::vector<double>     zKst;
        std::vector<double>     normChi2Kst;

        unsigned int nBd;
        std::vector<double>     ptBd;
        std::vector<double>     etaBd;
        std::vector<double>     phiBd;
        std::vector<double>     normChi2Bd;
            // Trajectory qualities
        std::vector<double>     dxy;
        std::vector<double>     dz;
        std::vector<double>     dcaBS;
        std::vector<double>     dcaBSErr;

            // Trigger and MC info
        std::vector<bool>       bIsGenMatched;

    private:

        // Variables and objects not in branch


        // cuttage
        bool isHandlesOk();
        void regBranches(const char prefix[]="BdInfo");

        // plant
        void clearVariables();
        void buildBranches(const edm::Event&);
        void calLS( double,double,double,
                    double,double,double,
                    double,double,double,
                    double,double,double,
                    double,double,double,
                    double,double,double,
                    double*,double*);
        void calCosAlpha(   double,double,double,
                            double,double,double,
                            double,double,double,
                            double,double,double,
                            double,double,double,
                            double,double,double,
                            double*,double*);
        void calCtau(RefCountedKinematicTree, double&, double&);
        bool calClosestApproachTracks(const reco::TransientTrack, const reco::TransientTrack, double&, double&, double&);
        bool hasGoodDcaBs(const reco::TransientTrack, double&, double&, double&);

};

#endif
