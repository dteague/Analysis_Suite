#ifndef COMMONFUNCS_H_
#define COMMONFUNCS_H_

#include <unordered_map>
#include <cmath>
#include <string>
#include <stdexcept>

template <class A, class B>
A get_by_val(std::unordered_map<A, B> map, B val_) {
    for (auto& [key, val] : map) {
        if (val == val_)
            return key;
    }
    throw std::out_of_range("Could not find value in map");
}


template <class T, class List, class... Args>
T* createObject(List list, std::string name, Args... args)
{
    T* obj = new T(name.c_str(), name.c_str(), args...);
    list->Add(obj);
    return obj;
}

float deltaR(float eta1, float eta2, float phi1, float phi2);

#endif // COMMONFUNCS_H_
