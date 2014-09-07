#ifndef __VARIABLE_H__
#define __VARIABLE_H__
enum class VariableType : unsigned char {
    Continuous,
    Discrete
};

union VariableValue {
    double d;
    long long ll;
};

class Variable {
    protected:
        VariableValue _val;
        VariableType _type;

    public:

        Variable() : _type(VariableType::Continuous) {_val.d = 0.5;}

        void val(double v) {
            if (_type == VariableType::Discrete) 
                throw new std::invalid_argument("Can't call with double");
            _val.d = v;
        }
        
        void val(long long v) {
            if (_type == VariableType::Continuous) 
                throw new std::invalid_argument("Can't call with ll");
            _val.ll = v;
        }

        operator double() {
            return _val.d;
        }

        operator long long() {
            return _val.ll;
        }
};

class VariableContainer : protected Variable {
    public:
        VariableValue rangeStart;
        VariableValue rangeEnd;
        VariableContainer() : _type(VariableType::Continuous) {rangeStart.d = 0.0, rangeEnd.d = 1.0, _val.d = 0.5;}

        void setType(const VariableType& vt) { _type = vt; }
        VariableType getType() { return _type; }
};

#endif // __VARIABLE_H__
