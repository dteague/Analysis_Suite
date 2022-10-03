#include "analysis_suite/skim/interface/Systematic.h"


size_t SystematicWeights::nSyst = 0;
Year SystematicWeights::year_ = Year::yrDefault;
std::string SystematicWeights::scaleDir_ = "";
bool SystematicWeights::isMC = true;
eVar SystematicWeights::currentVar = eVar::Nominal;
Systematic SystematicWeights::currentSyst = Systematic::Nominal;
