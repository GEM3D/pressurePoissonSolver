#ifndef GMGRestrictor_H
#define GMGRestrictor_H
#include "PW.h"
#include "petscvec.h"
class GMGRestrictor{
    public:
    virtual void restrict(PW<Vec> coarse, PW<Vec> fine)=0;
};
#endif
