#ifndef __PARTICLE_H_
#define __PARTICLE_H_

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "Math/GenVector/VectorUtil.h"
#include "analysis_suite/Analyzer/interface/Systematic.h"
#include "analysis_suite/Analyzer/interface/Variable.h"

#include <TTreeReader.h>

#include <vector>

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> LorentzVector;
typedef math::RhoEtaPhiVector Vector3D;
typedef std::vector<std::vector<size_t>> PartList;

class Particle : public SystematicWeights {
public:
    Particle() {};
    virtual ~Particle() {};

    virtual void setup(std::string name, TTreeReader& fReader) {};

    virtual size_t _size() const {return 0;};
    virtual Float_t _pt(size_t idx) const { return -1;};
    virtual Float_t _eta(size_t idx) const { return -1;};
    virtual Float_t _phi(size_t idx) const { return -1;};
    virtual Float_t _mass(size_t idx) const { return -1;};

    size_t idx(Level level, size_t i) const;

    size_t size() const { return _size(); };
    size_t size(Level level) const { return list(level).size(); }
    Float_t pt(size_t i) const { return _pt(i); }
    Float_t pt(Level level, size_t i) const { return pt(idx(level, i)); }
    Float_t eta(size_t i) const { return _eta(i); }
    Float_t eta(Level level, size_t i) const { return eta(idx(level, i)); }
    Float_t phi( size_t i) const { return _phi(i); }
    Float_t phi(Level level, size_t i) const { return phi(idx(level, i)); }
    Float_t mass(size_t i) const { return _mass(i); }
    Float_t mass(Level level, size_t i) const { return mass(idx(level, i)); }

    LorentzVector p4(size_t idx) const {return LorentzVector(pt(idx), eta(idx), phi(idx), mass(idx));}
    LorentzVector p4(Level level, size_t i) const { return p4(idx(level, i)); }
    Vector3D p3(size_t idx) const { return Vector3D(pt(idx), eta(idx), phi(idx)); }
    Vector3D p3(Level level, size_t i) const { return p3(idx(level, i)); }

    const std::vector<size_t>& list(Level level) const { return *m_partList.at(level); };
    const std::vector<size_t>& list(Level level, size_t syst) const { return m_partArray.at(level).at(syst); };
    virtual void clear();

    virtual void setupGoodLists() {std::cout << "SHOULDN'T BE CALLED" << std::endl;}
    virtual void setupGoodLists(Particle&) {std::cout << "SHOULDN'T BE CALLED (with jet)" << std::endl;}
    virtual void setupGoodLists(Particle&, Particle&) {std::cout << "SHOULDN'T BE CALLED (with jet)" << std::endl;}

    template <class... Args>
    void setGoodParticles(size_t syst, Args&&... args);

    virtual float getScaleFactor() { return 1.0; };
    const std::vector<size_t>& bitmap(Level level) const { return m_bitArray.at(level); };

    void moveLevel(Level level_start, Level level_end);
    void xorLevel(Level big, Level small, Level target);


protected:
    virtual void setup_map(Level level);

    std::unordered_map<Level, std::vector<size_t>*> m_partList;
    std::unordered_map<Level, PartList> m_partArray;
    std::unordered_map<Level, std::vector<size_t>> m_bitArray;
};

template <class... Args>
void Particle::setGoodParticles(size_t syst, Args&&... args)
{
    // Setup variables
    for (auto& [key, plist] : m_partArray) {
        m_partList[key] = &plist[syst];
        m_bitArray[key].assign(size(), 0);
    }

    // Virutal class specific list making
    setupGoodLists(std::forward<Args>(args)...);

    // Fill the bitmap
    for (const auto& [key, plist] : m_partArray) {
        for (size_t syst = 0; syst < nSyst; ++syst) {
            for (auto idx : plist[syst]) {
                m_bitArray[key][idx] += 1 << syst;
            }
        }
    }
}


#endif // __PARTICLE_H_
