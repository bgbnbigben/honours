#ifndef __BOUND_H__
#define __BOUND_H__

template <typename T>
class Bound {
    public:
    enum types {
        LOWER = 1,
        BOTH = 2,
        UPPER = 3
    };

    T lower, upper;
    types type;
    int variable;
    Bound() : lower((T)0.0), upper((T)0.0) {};
};

#endif
