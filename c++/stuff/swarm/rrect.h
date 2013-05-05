#ifndef RRECT_H
#define RRECT_H

#include <vector>

struct RRect {
    RRect(const std::vector<double>& point, bool rect = true) {
        min = new double[point.size()];
        max = new double[point.size()];
        for (auto i = 0u; i < point.size(); i++) {
            min[i] = point[i];
            max[i] = point[i];
            if (true || rect) {
                min[i] -= 5e-40;
                max[i] += 5e-40;
            }
        }
    }

    ~RRect() {
        delete[] min;
        delete[] max;
    }


    double* min;
    double* max;
}; 

#endif
