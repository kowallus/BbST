#ifndef BBST_HYBTEMPL_H
#define BBST_HYBTEMPL_H

#include "common.h"
#include "includes/RMQRMM64.h"

class RMQAPI {
public:
    virtual t_array_size rmq(const t_array_size &begIdx, const t_array_size &endIdx);
    virtual size_t memUsageInBytes();
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

    size_t memUsageInBytes() {
        return 0;
    }
};

class FNRMQBP: public RMQAPI {
private:
    RMQRMM64 *rmqImpl;
public:

    FNRMQBP(const t_value* valuesArray, const t_array_size n) {
        rmqImpl = new RMQRMM64((t_value*) valuesArray, n);
    }

    ~FNRMQBP() {
        delete rmqImpl;
    }

    t_array_size rmq(const t_array_size &begIdx, const t_array_size &endIdx) {
        return rmqImpl->queryRMQ(begIdx, endIdx);
    }

    size_t memUsageInBytes() {
        return rmqImpl->getSize();
    }
};



#endif //BBST_HYBTEMPL_H
