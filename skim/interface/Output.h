#ifndef __OUTPUT_H_
#define __OUTPUT_H_
// Cannot use inherited classes unfortunately
// To get around this, helper functions and macros reduce overhead
// of the main particle filling needs

// Used to fill basic kinematic variables
#define PARTICLE_VARS                           \
    std::vector<Float_t> pt;                    \
    std::vector<Float_t> eta;                   \
    std::vector<Float_t> phi;                   \
    std::vector<Float_t> mass;                  \
    std::vector<Int_t> syst_bitMap;

// Used to setup basic clearing functionality
template <class T>
void generic_clear(T* output)
{
    output->pt.clear();
    output->eta.clear();
    output->phi.clear();
    output->mass.clear();
    output->syst_bitMap.clear();
}


////////////////////
// Output Classes //
////////////////////


struct ParticleOut {
    PARTICLE_VARS;
    void clear() {
        generic_clear(this);
    }
};

struct JetOut {
    PARTICLE_VARS;
    std::vector<std::vector<Float_t>> jet_shifts;
    std::vector<Float_t> discriminator;
    std::vector<Int_t> pass_btag;
    void clear() {
        generic_clear(this);
        jet_shifts.clear();
        discriminator.clear();
        pass_btag.clear();
    }
};


struct TopOut {
    PARTICLE_VARS;
    std::vector<Float_t> disc;
    void clear() {
        generic_clear(this);
        disc.clear();
    }
};

struct LeptonOut {
    PARTICLE_VARS;
    std::vector<Bool_t> flip;
    void clear() {
        generic_clear(this);
        flip.clear();
    }
};

struct ElectronOut {
    PARTICLE_VARS;
    std::vector<std::pair<Float_t, Float_t>> dEScale;
    std::vector<std::pair<Float_t, Float_t>> dESigma;
    void clear() {
        generic_clear(this);
        dEScale.clear();
        dESigma.clear();
    }
};

struct ElectronOut_Endcap {
    PARTICLE_VARS;
    std::vector<Float_t> dxy;
    std::vector<Float_t> dz;
    std::vector<Bool_t> convVeto;
    std::vector<UChar_t> lostHits;
    std::vector<Float_t> sip3d;
    std::vector<Float_t> mvaTTH;
    std::vector<UChar_t> truth;
    void clear() {
        generic_clear(this);
        dxy.clear();
        dz.clear();
        convVeto.clear();
        lostHits.clear();
        sip3d.clear();
        mvaTTH.clear();
        truth.clear();
    }
};

struct LeptonOut_small {
    std::vector<Float_t> rawPt;
    std::vector<Float_t> ptRatio;
    std::vector<Float_t> mvaTTH;
    std::vector<Int_t> syst_bitMap;

    void clear() {
        rawPt.clear();
        ptRatio.clear();
        mvaTTH.clear();
        syst_bitMap.clear();
    }
};


struct LeptonOut_Fake {
    PARTICLE_VARS;
    std::vector<Float_t> rawPt;
    std::vector<Float_t> ptRatio;
    std::vector<Float_t> ptRatio2;
    std::vector<Float_t> mvaTTH;
    std::vector<Float_t> iso;
    std::vector<Float_t> jet_btag;
    std::vector<Int_t> jet_flav;
    void clear() {
        generic_clear(this);
        rawPt.clear();
        ptRatio.clear();
        ptRatio2.clear();
        mvaTTH.clear();
        iso.clear();
        jet_btag.clear();
        jet_flav.clear();
    }
};

struct ElectronOut_Fake {
    PARTICLE_VARS;
    std::vector<Float_t> rawPt;
    std::vector<Float_t> ptRatio;
    std::vector<Float_t> ptRatio2;
    std::vector<Float_t> old_mvaTTH;
    std::vector<Float_t> mvaTTH;
    std::vector<Float_t> iso;
    std::vector<Float_t> jet_btag;
    std::vector<Int_t> jet_flav;
    void clear() {
        generic_clear(this);
        rawPt.clear();
        ptRatio.clear();
        ptRatio2.clear();
        old_mvaTTH.clear();
        mvaTTH.clear();
        iso.clear();
        jet_btag.clear();
        jet_flav.clear();
    }
};


struct BEffOut {
    PARTICLE_VARS;
    std::vector<Int_t> flavor;
    std::vector<Int_t> pass_btag;
    void clear() {
        generic_clear(this);
        flavor.clear();
        pass_btag.clear();
    }
};

#endif // __OUTPUT_H_
