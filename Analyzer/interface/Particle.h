#ifndef __PARTICLE_H_
#define __PARTICLE_H_

#include "DataFormats/Math/interface/LorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <vector>

#include "analysis_suite/Analyzer/interface/CommonEnums.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>
    LorentzVector;
typedef std::vector<std::vector<size_t>> PartList;

template <class T>
using TTRArray = TTreeReaderArray<T>;

class ParticleOut;
class BJetOut;
class TopOut;

class Particle {
public:
    Particle(){};
    virtual ~Particle(){};

    void setup(std::string name, TTreeReader& fReader);
    virtual void setGoodParticles(size_t syst);

    virtual float getScaleFactor() { return 1.0; };

    size_t size() const { return (m_pt) ? m_pt->GetSize() : 0; }
    size_t size(Level level) const { return list(level).size(); }
    Float_t pt(size_t idx) const { return m_pt->At(idx); }
    Float_t eta(size_t idx) const { return m_eta->At(idx); }
    Float_t phi(size_t idx) const { return m_phi->At(idx); }
    Float_t mass(size_t idx) const { return m_mass->At(idx); }

    const std::vector<size_t>& list(Level level) const { return *m_partList.at(level); };
    const std::vector<size_t>& list(Level level, size_t syst) const { return m_partArray.at(level).at(syst); };
    const std::vector<size_t>& bitmap(Level level) const { return m_bitArray.at(level); };

    void moveLevel(Level level_start, Level level_end);
    virtual void clear();

    static size_t nSyst;
    static TFile* f_scale_factors;
    static Year year_;
    static std::string yearStr_;
    static std::string scaleDir_;

protected:
    TTRArray<Float_t>* m_pt;
    TTRArray<Float_t>* m_eta;
    TTRArray<Float_t>* m_phi;
    TTRArray<Float_t>* m_mass;

    std::unordered_map<Level, std::vector<size_t>*> m_partList;
    std::unordered_map<Level, PartList> m_partArray;
    std::unordered_map<Level, std::vector<size_t>> m_bitArray;

    void setup_map(Level level);
    void fill_bitmap();

    template <class T>
    void setSF(T*& sf, std::string name)
    {
        sf = static_cast<T*>(f_scale_factors->Get((yearStr_ + "/" + name).c_str()));
    }

    template <class T, class... Args>
    float getWeight(T* hist, Args... args)
    {
        return hist->GetBinContent(hist->FindBin(args...));
    }
};

#endif // __PARTICLE_H_
