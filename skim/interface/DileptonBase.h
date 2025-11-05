#ifndef DILEPTONBASE_H_
#define DILEPTONBASE_H_

#include "analysis_suite/skim/interface/BaseSelector.h"

#include <nlohmann/json.hpp>

class DileptonBase : public BaseSelector {
public:
    virtual void Init(TTree* tree) override;
    void setup_trigger();
    float getLeadPt();
    void setSubChannel();
    bool isSameSign(Level level);
};





#endif // DILEPTONBASE_H_
