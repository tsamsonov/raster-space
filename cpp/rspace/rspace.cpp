// [[Rcpp::plugins(cpp14)]]

#include<vector>
#include<cmath>
#include<iostream>
#include<pybind11/pybind11.h>
#include<pybind11/numpy.h>
#include<pybind11/stl.h>

namespace py = pybind11;
using Shifts = std::vector<std::pair<int, int>>;

// Bresenham's algorithm by X
Shifts get_shifts_byX(int x0, int y0, int x1, int y1) {
    int dx = x1 - x0;
    int dy = y1 - y0;
    int yi = 1;
    if (dy < 0) {
        yi = -1;
        dy = -dy;
    }

    int D = (2 * dy) - dx;
    int y = y0;

    Shifts s;
    for (auto x = x0; x <= x1; x++) {
        s.emplace_back(x, y);
        if (D > 0) {
            y += yi;
            D += 2 * (dy - dx);
        } else {
            D += 2 * dy;
        }
    }

    return s;
}

// Bresenham's algorithm by Y
Shifts get_shifts_byY(int x0, int y0, int x1, int y1) {
    int dx = x1 - x0;
    int dy = y1 - y0;
    int xi = 1;
    if (dx < 0) {
        xi = -1;
        dx = -dx;
    }

    int D = (2 * dx) - dy;
    int x = x0;

    Shifts s;
    for (auto y = y0; y <= y1; y++) {
        s.emplace_back(x, y);
        if (D > 0) {
            x += xi;
            D += 2 * (dx - dy);
        } else {
            D += 2 * dx;
        }
    }

    return s;
}

// Bresenham's algorithm in any direction
Shifts get_shifts(int x0, int y0, int x1, int y1) {
    if (abs(y1 - y0) < abs(x1 - x0)) {
        if (x0 > x1) {
            auto s = get_shifts_byX(x1, y1, x0, y0);
            std::reverse(std::begin(s), std::end(s));
            return s;
        } else {
            return get_shifts_byX(x0, y0, x1, y1);
        }
    } else {
        if (y0 > y1) {
            auto s = get_shifts_byY(x1, y1, x0, y0);
            std::reverse(std::begin(s), std::end(s));
            return s;
        } else {
            return get_shifts_byY(x0, y0, x1, y1);
        }
    }
}

std::pair<int, int> operator +(const std::pair<int, int>& x, const std::pair<int, int>& y) {
    return std::make_pair(x.first + y.first, x.second + y.second);
}

double dist(int dx, int dy) {
    return sqrt(dx * dx + dy * dy);
}

bool is_within(int i, int j, int imax, int jmax) {
    return (i >= 0) and (j >= 0) and (i < imax) and (j < jmax);
}

// Region length calculation
std::vector<std::vector<double>> estimate_length(py::array_t<double> cobst, py::array_t<double> cwidth,
                                    double cellsize, double nodata, int ndirs, double radius) {

    std::vector<Shifts> shifts;
    auto r = int(ceil(radius/cellsize));
    auto angle = 2 * M_PI / ndirs;
    double a = M_PI_2;
    int i, j;

    for (auto k = 0; k < ndirs; k++) {
        i = -r * cos(a);
        j = r * sin(a);
        shifts.push_back(std::move(get_shifts(0, 0, j, i)));
        a -= angle;
    }

    auto buf_obst = cobst.request();
    auto buf_width = cwidth.request();

//    if (buf_obst.size != buf_length.size)
//        throw std::runtime_error("Input shapes must match");

    auto nrow = buf_obst.shape[0];
    auto ncol = buf_obst.shape[1];

    auto *obst = (double *) buf_obst.ptr;
    auto *width = (double *) buf_width.ptr;

    std::vector<std::vector<double>> output(4, std::vector<double>(nrow * ncol));

    int ndirs2 = ndirs / 2;

    int ik, jk;

    std::vector<double> d1(ndirs2); // reusable vector;
    std::vector<double> d2(ndirs2); // reusable vector;

    std::vector<double> d12(ndirs2); // reusable vector;
    std::vector<double> w12(ndirs2); // reusable vector;
    std::vector<double> n12(ndirs2); // reusable vector;

    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {

            for (auto k = 0; k < ndirs2; k++) {
                n12[k] = 1;
                w12[k] = width[i * ncol + j];

                d1[k] = r;
                for (auto s: shifts[k]) {
                    ik = i + s.first;
                    jk = j + s.second;
                    if (is_within(ik, jk, nrow, ncol)) {
                        if (obst[ik * ncol + jk] == 1) {
                            d1[k] = dist(s.first, s.second);
                            break;
                        } else {
                            w12[k] += width[ik * ncol + jk];
                            n12[k]++;
                        }
                    } else {
                        d1[k] = dist(s.first, s.second);
                        break;
                    }

                }

                d2[k] = r;
                for (auto s: shifts[ndirs2 + k]) {
                    ik = i + s.first;
                    jk = j + s.second;
                    if (is_within(ik, jk, nrow, ncol)) {
                        if (obst[ik * ncol + jk] == 1) {
                            d2[k] = dist(s.first, s.second);
                            break;
                        } else {
                            w12[k] += width[ik * ncol + jk];
                            n12[k]++;
                        }
                    } else {
                        d2[k] = dist(s.first, s.second);
                        break;
                    }
                }

                d12[k] = d1[k] + d2[k];
            }

//            length[i * ncol + j] = *std::max_element(d12.begin(), d12.end());
//            length[i * ncol + j] = std::distance(d12.begin(), std::max_element(d12.begin(), d12.end()));

            auto max_iterator = std::max_element(d12.begin(), d12.end());
            auto max_position = std::distance(d12.begin(), max_iterator);
            auto max_length = *max_iterator;

//            auto max_iterator = std::max_element(w12.begin(), w12.end());
//            auto max_position = std::distance(w12.begin(), max_iterator);
//            auto max_length = d12[max_position];

            output[0][i * ncol + j] = max_length * cellsize;
            output[1][i * ncol + j] = 180 * angle * max_position / M_PI;
            output[2][i * ncol + j] = max_length > 0 ? d1[max_position] / max_length : 1;
            output[3][i * ncol + j] = w12[max_position] / n12[max_position];

        }
    }

    return output;
}

// Region width calculation
py::array_t<double> estimate_width(py::array_t<double> ceuc, py::array_t<double> cwidth, double cellsize, double nodata) {


    auto buf_euc = ceuc.request(), buf_width = cwidth.request();

    if (buf_euc.size != buf_width.size)
        throw std::runtime_error("Input shapes must match");

    auto nrow = buf_euc.shape[0];
    auto ncol = buf_euc.shape[1];

    auto *euc = (double *) buf_euc.ptr,
         *width = (double *) buf_width.ptr;

    for (int i = 0; i < nrow; i++) {

        for (int j = 0; j < ncol; j++) {

            auto radius = euc[i * ncol + j];

            if (radius == nodata) {
                width[i * ncol + j] = nodata;
                continue;
            }

            auto w = int(ceil(radius / cellsize));

            for (int k = -w; k <= w; ++k) {
                for (int l = -w; l <= w; ++l) {

                    if (k * k + l * l > w * w)
                        continue;

                    auto ik = i + k;
                    auto jl = j + l;

                    if (ik < 0 || ik >= nrow || jl < 0 || jl >= ncol)
                        continue;

                    if (width[ik * ncol + jl] < 2 * radius)
                        width[ik * ncol + jl] = 2 * radius;
                }
            }
        }
    }

    return cwidth;
}

PYBIND11_MODULE(rspace, m) {
    m.doc() = R"pbdoc(
        C++ plugin for region width estimation
        -----------------------

        .. currentmodule:: rspace

        .. autosummary::
           :toctree: _generate

           estimate_width, estimate_length

    )pbdoc";

    m.def("estimate_width", &estimate_width, R"pbdoc(
        Estimate local region width
    )pbdoc");

    m.def("estimate_length", &estimate_length, R"pbdoc(
        Estimate local region length
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
