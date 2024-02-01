#include "analysis_suite/skim/interface/CommonFuncs.h"
#include <DataFormats/Math/interface/deltaPhi.h>

#include <math.h>


float deltaR(float eta1, float eta2, float phi1, float phi2)
{
    return pow(eta1-eta2, 2) + pow(deltaPhi(phi1, phi2), 2);
}
