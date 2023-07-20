#ifndef BASESELECTOR_H
#define BASESELECTOR_H

#include <TSelector.h>
#include <TTree.h>
#include <TTreeReader.h>

#include <exception>
#include <iostream>
#include <vector>

#include "analysis_suite/skim/interface/CommonEnums.h"
#include "analysis_suite/skim/interface/Electron.h"
#include "analysis_suite/skim/interface/Jet.h"
#include "analysis_suite/skim/interface/Muon.h"
#include "analysis_suite/skim/interface/Met.h"
#include "analysis_suite/skim/interface/GenParticle.h"
#include "analysis_suite/skim/interface/ScaleFactors.h"
#include "analysis_suite/skim/interface/Variable.h"
#include "analysis_suite/skim/interface/CommonStructs.h"
#include "analysis_suite/skim/interface/Bar.h"

enum class Channel;
enum class Subchannel;

class BaseSelector : public TSelector {
public:
    BaseSelector(TTree* /*tree*/ = 0) {}
    virtual ~BaseSelector() {}
    virtual Int_t Version() const { return 2; }

    virtual void Init(TTree* tree);
    virtual Bool_t Process(Long64_t entry);
    virtual void SlaveTerminate();

    /**
     * @brief Set output file used by python code
     *
     * @param outfile_ Output TFile to put trees and info
     **/
    void setOutputFile(TFile* outfile_) { outfile = outfile_; }

    /**
     * @brief
     **/
    std::vector<TreeInfo> getTrees() {
        std::vector<TreeInfo> tree_vec;
        for (auto& [chan, tree]: trees) {
            tree_vec.push_back(tree);
        }
        return tree_vec;
    }
    TDirectory* getOutdir() { return outdir; }

    ClassDef(BaseSelector, 0);

protected:
    size_t numSystematics() { return systematics_.size()*2 - 1; }
    virtual void clearParticles();
    virtual void clearOutputs(){};
    void setupSyst(size_t systNum);

    void createTree(std::string name, Channel chan)
    {
        trees.emplace(chan, TreeInfo(name, outdir, fOutput));
        SetupOutTreeBranches(trees.at(chan).tree);
    }

    void setupTrigger(Subchannel sChan, Dataset dataset, std::vector<std::string> trigs = {}) {
        trig_cuts.setup_channel(sChan, dataset, fReader, trigs);
    }

    void fillCutFlow(Channel chan, CutInfo& cuts) {
        if (trees.find(chan) != trees.end())
            trees.at(chan).fillCutFlow(cuts, *weight);
    }

    size_t nLeps(Level level) { return muon.size(level) + elec.size(level); }

    // To be filled by Child class
    virtual void ApplyScaleFactors(){};
    virtual void FillValues(const Bitmap& event_bitmap) {};
    virtual void setOtherGoodParticles(size_t syst) {};
    virtual bool getCutFlow() { return true; }
    virtual void SetupOutTreeBranches(TTree* tree);

    // Protected Variables
    TTreeReader fReader;
    std::string groupName_;
    TFile* outfile;
    Long64_t stop_point = -1;

    std::map<Channel, TreeInfo> trees;

    TRVariable<Float_t> genWeight;
    TRVariable<Int_t> PV_npvs;

    Year year_;
    bool isMC_ = true;
    Subchannel subChannel_;
    Dataset dataset_ = Dataset::None;

    // Current weight and all weights
    float* weight;
    std::vector<Float_t> o_weight;

    // Current channel and all channels
    Channel* currentChannel_;
    std::vector<Channel> o_channels;

    // Systs that change stuff
    std::vector<size_t> syst_to_index = { 0 };
    std::vector<size_t> novel_systs = { 0 };
    size_t n_novel_systs = 1;

    std::vector<Bool_t> o_pass_event;

    OutValue<UInt_t> run, lumiblock;
    OutValue<ULong64_t> event;

    VarCollection<Bool_t> metfilters;

    ScaleFactors sfMaker;
    Muon muon;
    Electron elec;
    Jet jet;
    GenParticle rGen;
    GenJet rGenJet;
    Met met;

    MET_Type met_type = MET_Type::PF;

    TriggerInfo trig_cuts;

private:
    void SetupEvent(size_t systNum);
    void setupSystematicInfo();

    std::vector<Systematic> systematics_ = { Systematic::Nominal };
    TDirectory* outdir;
    Bar bar;
    std::vector<std::pair<Systematic, eVar>> syst_var_pair;
};

#endif
