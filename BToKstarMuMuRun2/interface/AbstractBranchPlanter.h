#ifndef ABSTRACT_BRANCH_PLANTER_H
#define ABSTRACT_BRANCH_PLANTER_H

#include "BphAna/BToKstarMuMuRun2/interface/CommonBDecayAssistant.h"

#include <map>
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

struct HistArgs{
    const char      name [128];
    const char      title[128];
    unsigned int    nBins;
    double          x_min;
    double          x_max;
};

class AbstractBranchPlanter{
    /*
     * cuttage, plant, and reap. That's how branches grow up.
     * cuttage      : Put branches to a known tree.
     *                Skip this step if you don't want to keep the tree.
     * plant        : Fill branch contents.
     * reap         : Build index and post-processing.
     *
     * In derived class, at least three functions to be defined:
     *      regBranches
     *      clearVariables
     *      buildBranches
     *
     * It's better to define following
     *      setBranchAddress                # Retrieving the branches
     *      writePSetDescription            # Write PSetDescription in EDAnalyzer::fillDescription
     *      getClassName
     *      getMotherClassName
     *
     *
     * It's also useful to cooperate with other trees with
     * TTree::GetEntryWithIndex(ParentTree.runNo,ParentTree.evtNo)
     * Ref: https://root.cern.ch/root/html/tutorials/tree/treefriend.C.html
     */

    public:
        AbstractBranchPlanter();
        ~AbstractBranchPlanter();

        void cuttage(TTree*, const char prefix[] = "AbsInfo");
        void plant(const edm::Event&, const edm::EventSetup&);
        void reap(const char prefix[] = "AbsInfo");

        // setters, getters, and adders
        bool    setHistFile(TFile*);
        bool    setBranchAddress(TTree*, const char prefix[] = "AbsInfo");
        void    linkAssistant(CommonBDecayAssistant*);

        TTree*                      getTree() const {return tree_;}
        TFile*                      getHistFile() const {return file_;}
        TH1D*                       getHist(std::string) const;
        std::vector<std::string>    getHistList() const;
        void                        getEntryWithIndex(unsigned int, unsigned int);

        bool                        addHist(const TH1D&);
        bool                        addHist(const char*, const char*, unsigned int, double, double);
        bool                        addHist(HistArgs&);

        // services
        virtual const std::string getClassName() const {return "AbstractBranchPlanter";}
        virtual const std::string getMotherClassName() const {return "";}
        virtual void print() const;
        //static void writePSetDescription(edm::ParameterSetDescription&);

        // Variables defined as PUBLIC in daughter class for convenient accessment.
        unsigned int runNo = 0;
        unsigned int evtNo = 0;

    protected:
        // member data
        TFile                *file_     =0;
        TTree                *tree_     =0;
        std::map<std::string,TH1D*> histMap_;
        CommonBDecayAssistant *assistant_ = 0;

        // cuttage
        // regBranches is replaced in derived class. keyword: 'runtime-binding'
        virtual bool isHandlesOk();
        virtual void regBranches(const char prefix[]="AbsInfo");
        //void regHists(const char prefix[]="AbsInfo");

        // plant, these pure virtual functions to be defined.
        virtual void clearVariables() = 0;
        virtual void buildBranches(const edm::Event&, const edm::EventSetup&) = 0;

    private:
};

#endif
