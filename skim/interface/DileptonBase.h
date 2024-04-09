#ifndef DILEPTONBASE_H_
#define DILEPTONBASE_H_

#include "analysis_suite/skim/interface/BaseSelector.h"

#include <nlohmann/json.hpp>

enum class Subchannel {
    MM,
    EM,
    EE,
    Single_E,
    Single_M,
    None,
};

class DileptonBase : public BaseSelector {
public:
    virtual void Init(TTree* tree) override;
    void setup_trigger();
    bool getTriggerValue();
    float getLeadPt();
    void setSubChannel();
    bool isSameSign(Level level);

    nlohmann::json dilep_trigs;
};





#endif // DILEPTONBASE_H_
