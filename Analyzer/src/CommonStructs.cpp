#include "analysis_suite/Analyzer/interface/CommonStructs.h"
#include "analysis_suite/Analyzer/interface/CommonFuncs.h"


TreeInfo::TreeInfo(std::string name, std::set<Channel> goodChannels_, TDirectory* outdir, TSelectorList* fOutput) :
    goodChannels(goodChannels_)
{
    tree = createObject<TTree>(fOutput, name);
    tree->SetDirectory(outdir);
    std::string treeName = tree->GetName();
    cutflow = createObject<TH1F>(fOutput, "cutFlow_"+treeName, 15, 0, 15);
    cutflow_ind = createObject<TH1F>(fOutput, "cutFlow_ind_"+treeName, 15, 0, 15);
}

void TreeInfo::fillCutFlow(CutInfo& cuts, Channel chan, float weight) {
    // Setup cutflow histogram
    if (!initialize_axis) {
        initialize_axis = true;
        for (size_t i = 0; i < cuts.size(); ++i) {
            cutflow->GetXaxis()->SetBinLabel(i+1, cuts.name(i).c_str());
            cutflow_ind->GetXaxis()->SetBinLabel(i+1, cuts.name(i).c_str());
        }
        cutflow->GetXaxis()->SetBinLabel(cuts.size()+1, "PassChannel");
        cutflow_ind->GetXaxis()->SetBinLabel(cuts.size()+1, "PassChannel");
    }

    bool passAll = true;
    for (auto& [cutName, truth] : cuts.cuts) {
        passAll &= truth;
        if (truth)
            cutflow_ind->Fill(cutName.c_str(), weight);
        if (passAll)
            cutflow->Fill(cutName.c_str(), weight);
    }
    if (contains(chan)) {
        cutflow_ind->Fill("PassChannel", weight);
        if (passAll)
            cutflow->Fill("PassChannel", weight);
    }
}
