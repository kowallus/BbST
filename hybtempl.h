#ifndef BBST_HYBTEMPL_H
#define BBST_HYBTEMPL_H

#include "common.h"

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

class RMQCounterDecorator: public RMQAPI {
private:
    RMQAPI* coreRMQ;
    uint64_t counter = 0;
    double sumLogM = 0;
public:

    RMQCounterDecorator(RMQAPI* coreRMQ):
            coreRMQ(coreRMQ) {}

    t_array_size rmq(const t_array_size &begIdx, const t_array_size &endIdx) {
        counter++;
        sumLogM += log2(endIdx - begIdx + 1);
//        cout << sumLogM << "\t" << endIdx - begIdx << endl;
        return coreRMQ->rmq(begIdx, endIdx);
    }

    void resetCounter() {
        counter = 0;
    }

    uint64_t getRMQCount() {
        return counter;
    }

    double getLogAvaregedM() {
        return pow(2, sumLogM / counter);
    }

    size_t memUsageInBytes() {
        return coreRMQ->memUsageInBytes();
    }
};


#endif //BBST_HYBTEMPL_H
