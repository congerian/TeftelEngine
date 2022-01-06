#pragma once


#include "../public/util.hxx"

namespace TFN{

class Particle {

    public:

    FP mass;
    FP effectiveRadius;

    FP vx;
    FP vy;

    FP x;
    FP y;

    Particle(  
        FP _x = 0.0, FP _y = 0.0, FP _vx = 0.0, FP _vy = 0.0,
        const FP& _mass = 1, const FP& _effectiveRadius = 1) :
        
        x{_x}, y{_y}, vx{_vx}, vy{_vy}, mass{_mass}, effectiveRadius{_effectiveRadius}
    {}

    Particle& operator= (const Particle& other){
        x = other.x;
        y = other.y;
        vx = other.vx;
        vy = other.vy;
        mass = other.mass;
        effectiveRadius = other.effectiveRadius;
        return *this;
    }
};

}

