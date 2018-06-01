#ifndef GMGInterpolator_H
#define GMGInterpolator_H
#include "PW.h"
#include "petscvec.h"
class GMGInterpolator{
    public:
    virtual void interpolate(PW<Vec> coarse, PW<Vec> fine)=0;
};
#endif
