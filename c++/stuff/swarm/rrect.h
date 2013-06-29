#ifndef RRECT_H
#define RRECT_H

#include <vector>
#include <algorithm>

struct RRect {
    RRect(const std::vector<double>& point, bool rect = true) {
        min = new double[point.size()];
        max = new double[point.size()];
        for (auto i = 0u; i < point.size(); i++) {
            min[i] = point[i];
            max[i] = point[i];
            if (true || rect) {
                constexpr double side_length = 5e-5;
                const double dist = side_length * std::sqrt(point.size()) / 2.0;
                min[i] -= dist;
                max[i] += dist;
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
