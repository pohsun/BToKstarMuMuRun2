#ifndef EVENT_BRANCH_PLANTER_H
#define EVENT_BRANCh_PLANTER_H

#include "BphAna/BToKstarMuMuRun2/interface/AbstractBranchPlanter.h"
#include "BphAna/BToKstarMuMuRun2/interface/CommonBDecayAssistant.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

 // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideHighLevelTrigger#Access_to_the_HLT_configuration

class EvtBranchPlanter : public AbstractBranchPlanter{
    public:
        EvtBranchPlanter();
        ~EvtBranchPlanter();

        // setters, getters, and adders
        bool setBranchAddress(TTree*,const char prefix[]="EvtInfo");

        // services
        virtual const std::string getClassName()        {return "EvtBranchPlanter";}
        virtual const std::string getMotherClassName()  {return "AbstractBranchPlanter";}
        virtual void print();

        // Branch variables defined as PUBLIC for convenient accessment.
        unsigned int lumiBlock = 0;
        unsigned int nPV       = 0;
        std::vector<std::string> *triggeredNames = 0;
        std::vector<int>         *triggeredPrescales = 0;

    private:
        // Variables and objects not in branch

        // cuttage
        bool isHandlesOk();
        void regBranches(const char prefix[]="EvtInfo");
        //void regHists(const char prefix[]="EvtInfo");

        // plant
        void clearVariables();
        void buildBranches(const edm::Event&, const edm::EventSetup&);
        void hltReport(const edm::Event&);

};
#endif
