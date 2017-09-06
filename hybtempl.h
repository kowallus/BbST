#ifndef BBST_HYBTEMPL_H
#define BBST_HYBTEMPL_H

#include "common.h"

class RMQAPI {
public:
    virtual t_array_size rmq(const t_array_size &begIdx, const t_array_size &endIdx);
};

class RMQCounter: public RMQAPI {
private:
    uint64_t counter = 0;
public:
    t_array_size rmq(const t_array_size &begIdx, const t_array_size &endIdx) {
        counter++;
        return MAX_T_ARRAYSIZE;
    }

    void resetCounter() {
        counter = 0;
    }

    uint64_t getRMQCount() {
        return counter;
    }
};



#endif //BBST_HYBTEMPL_H
