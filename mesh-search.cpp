#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <utilities/Variable.h>

const double EPS = 1e-10;


//class Function;

class Problem {
    private:
        std::vector<Variable> _meshCentre;
        int _constrictions;
        struct window {
            VariableValue left;
            VariableValue right;
        };

    protected:
        std::vector<Variable> _variables;
        double _constrictFactor;
        double _lengthTolerance;
        double _sideWidth;

    public:

        //Function* f;

        Problem(std::vector<VariableContainer>& v, double constrictionFactor, double sideWidth, double length) : _variables(v), _constrictFactor(constrictionFactor), _constrictions(0), _sideWidth(sideWidth), _lengthTolerance(length) {
            _meshCentre.resize(_variables.size());
            for (int i = 0; i < _variables.size(); i++)
                if (_variables[i].type == VariableType::Discrete) 
                    _meshCentre[i].val((long long)_variables[i]);
                else
                    _meshCentre[i].val(_variables[i]);
        };

        void drive() {
            VariableValue left, right;
            std::vector<VariableValue> currentPoint(_meshCentre.size());
            std::vector<double> carry(_variables.size());
            std::vector<window> windows(_variables.size());
            long long nodesInCurrentRow;
            long long nodesPerRow;
            long long numNodes;

            /* set up some "helper" functions. They're not on the class because
             * wtf why would they be there? They're not at all useful to
             * anything else
             */
            auto nodesInRow = [&] (int i) {
                return (_variables[i].type == VariableType::Discrete ?
                            std::min(nodesPerRow,
                                _variables[i].rangeEnd.ll - _variables[i].rangeStart.ll + 1ll) :
                            nodesPerRow);
            };
            auto getWindowForVariable = [&] (int idx) {
                nodesInCurrentRow = nodesInRow(idx);
                window w;
                if (_variables[idx].type == VariableType::Continuous) {
                    w.left.d = std::max(_variables[idx].rangeStart.d, (double)_meshCentre[idx] - _sideWidth/2.0);
                    w.right.d = std::min(_variables[idx].rangeEnd.d, (double)_meshCentre[idx] + _sideWidth/2.0);
                } else if (_variables[idx].type == VariableType::Discrete) {
                    w.left.ll = std::max(_variables[idx].rangeStart.ll, (long long)_meshCentre[idx] - nodesInCurrentRow / 2);
                    w.right.ll = std::min(_variables[idx].rangeEnd.ll, (long long)_meshCentre[idx] + nodesInCurrentRow / 2);
                }
                return w;
            };
            auto incrementVariable = [&](int i) {
                if (_variables[i].type == VariableType::Discrete) {
                    currentPoint[i].ll += (long long)((windows[i].right.ll - windows[i].left.ll) / (double) (nodesInCurrentRow - 1));
                    // For discrete vars, add with carry to deal with possible
                    // windows-size problems (right - left > numNodes).
                    // This logic could be horribly horribly fucked.
                    double _ = (double)currentPoint[i].ll;
                    _ += carry[i];
                    // I don't want to round() this, I want to floor + deal with
                    // potential double-precision problems
                    if (floor(_ + EPS) != floor(_)) {
                        carry[i] = 0.0;
                        currentPoint[i].ll = (long long)ceil(_);
                    } else {
                        carry[i] = _ - floor(_);
                        currentPoint[i].ll = (long long)floor(_);
                    }
                } else {
                    currentPoint[i].d += ((windows[i].right.d - windows[i].left.d) / (double) (nodesInCurrentRow - 1));
                }
            };

            for (int z = 0; z < 3; z++, _constrictions++) {
                std::fill(carry.begin(), carry.end(), 0.0);
                nodesPerRow = pow((1/_constrictFactor), _constrictions) + 1;
                numNodes = 1;

                for (int i = 0; i < _variables.size(); i++) {
                    windows[i] = getWindowForVariable(i);
                    numNodes *= nodesInCurrentRow;
                    currentPoint[i] = windows[i].left;
                }

                std::cout << numNodes << std::endl;

                // This is stupid, but I'm doing it anyway to fix some
                // stupidity:
                VariableValue _ = currentPoint[0];
                nodesInCurrentRow = nodesInRow(0);
                incrementVariable(0);

                if (currentPoint[0] - _ <= _lengthTolerance) return;

                if (_variables[0].type == VariableType::Continuous)
                    currentPoint[0].d = _.d - (currentPoint[0].d - _.d);
                else
                    currentPoint[0].ll = _.ll - (currentPoint[0].ll - _.ll);

                for (long long i = 0; i < numNodes; i++) {
                    nodesInCurrentRow = nodesInRow(0);
                    incrementVariable(0);
                    for (int j = 0; j < _variables.size(); j++) {
                        //nodesInCurrentRow = nodesInRow(j); 
                        // TODO hack
                        if ((_variables[j].type == VariableType::Continuous && currentPoint[j].d > windows[j].right.d) ||
                            (_variables[j].type == VariableType::Discrete && currentPoint[j].ll > windows[j].right.ll)) {
                            currentPoint[j] = windows[j].left;
                            if (j < _variables.size() - 1) {
                                nodesInCurrentRow = nodesInRow(j+1);
                                incrementVariable(j+1);
                            }
                        } else break;
                    }

                    std::cout << "(";
                    for (int j = 0; j < _variables.size(); j++)
                        std::cout << (_variables[j].type == VariableType::Discrete ? currentPoint[j].ll : currentPoint[j].d) << (j == _variables.size() - 1 ? "" : ",\t\t");
                    std::cout << ")\n";
                }

            }
        }

};

int main() {
    std::vector<VariableContainer> variables(3);
    variables[1].setType(VariableType::Discrete);
    variables[1].rangeStart.ll = 0;
    variables[1].rangeEnd.ll = 3;
    variables[1].val(1ll);
    variables[2].setType(VariableType::Discrete);
    variables[2].rangeStart.ll = 0;
    variables[2].rangeEnd.ll = 10;
    variables[2].val(5ll);

    Problem problem(variables, .5, .5, .01);

    problem.drive();
    return 0;
}
