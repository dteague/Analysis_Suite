#include "analysis_suite/skim/interface/BaseSelector.h"
#include "analysis_suite/skim/interface/ThreeTop.h"
#include "analysis_suite/skim/interface/FakeRate.h"
#include "analysis_suite/skim/interface/Nonprompt_Closure.h"
#include "analysis_suite/skim/interface/Lepton_MisId_Closure.h"
#include "analysis_suite/skim/interface/BEfficiency.h"
#include "analysis_suite/skim/interface/DY_test.h"
#include "analysis_suite/skim/interface/SingleLep_Trigger.h"
#include "analysis_suite/skim/interface/TwoLepton.h"
#include "analysis_suite/skim/interface/TagProbe.h"
#include "analysis_suite/skim/interface/trigger_eff.h"
#include "analysis_suite/skim/interface/BScale.h"

namespace {
namespace {
    BaseSelector pBaseSelector;
    ThreeTop pThreeTop;
    FakeRate pFakeRate;
    Nonprompt_Closure pNonprompt_Closure;
    Closure_MisId pClosure_MisId;
    BEfficiency pBEfficiency;
    DY_test pDY_test;
    SingleLep_Trigger pSingleLep_Trigger;
    TwoLepton pTwoLepton;
    TagProbe pTagProbe;
    Trigger_Eff pTrigger_Eff;
    BScale pBScale;
} // namespace
} // namespace
