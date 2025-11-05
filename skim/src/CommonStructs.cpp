#include "analysis_suite/skim/interface/CommonStructs.h"
#include "analysis_suite/skim/interface/CommonFuncs.h"


TreeInfo::TreeInfo(std::string name, TDirectory* outdir, TSelectorList* fOutput)
{
    tree = new TTree(name.c_str(), name.c_str());
    // tree->SetDirectory(outdir);
    fOutput->Add(tree);
}

void TreeInfo::fillCutFlow(CutInfo& cuts, float weight) {
    return;
}
