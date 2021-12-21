// [[Rcpp::plugins(cpp14)]]

#include<vector>
#include<cmath>
#include<iostream>
#include<fstream>
#include<queue>
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
        shifts.push_back(get_shifts(0, 0, j, i));
        a -= angle;
    }

    auto buf_obst = cobst.request();
    auto buf_width = cwidth.request();

    auto nrow = buf_obst.shape[0];
    auto ncol = buf_obst.shape[1];

    auto *obst = (double *) buf_obst.ptr;
    auto *width = (double *) buf_width.ptr;

    std::vector<std::vector<double>> output(5, std::vector<double>(nrow * ncol));

    int ndirs2 = ndirs / 2;

    int ik, jk;

    // distances
    std::vector<double> d1(ndirs2); // reusable vector;
    std::vector<double> d2(ndirs2); // reusable vector;
    std::vector<double> d12(ndirs2); // reusable vector;

    // obstacle heights
    std::vector<double> cos1(ndirs2); // reusable vector;
    std::vector<double> cos2(ndirs2); // reusable vector;

    // widths and counts
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
                        if (obst[ik * ncol + jk] > 0) {
                            d1[k] = dist(s.first, s.second);
                            cos1[k] = ((ik == i) and (jk == j)) ? 1 : pow(cos(atan2(obst[ik * ncol + jk], d1[k] * cellsize)), 2);
                            break;
                        } else {
                            w12[k] += width[ik * ncol + jk];
                            n12[k]++;
                        }
                    } else {
                        d1[k] = dist(s.first, s.second);
                        cos1[k] = pow(cos(atan2(obst[ik * ncol + jk], d1[k] * cellsize)), 2);
                        break;
                    }

                }

                d2[k] = r;
                for (auto s: shifts[ndirs2 + k]) {
                    ik = i + s.first;
                    jk = j + s.second;
                    if (is_within(ik, jk, nrow, ncol)) {
                        if (obst[ik * ncol + jk] > 0) {
                            d2[k] = dist(s.first, s.second);
                            cos2[k] = ((ik == i) and (jk == j)) ? 1 : pow(cos(atan2(obst[ik * ncol + jk], d2[k] * cellsize)), 2);
                            break;
                        } else {
                            w12[k] += width[ik * ncol + jk];
                            n12[k]++;
                        }
                    } else {
                        d2[k] = dist(s.first, s.second);
                        cos2[k] = pow(cos(atan2(obst[ik * ncol + jk], d2[k] * cellsize)), 2);
                        break;
                    }
                }

                d12[k] = d1[k] + d2[k];
            }

            auto max_iterator = std::max_element(d12.begin(), d12.end());
            auto max_position = std::distance(d12.begin(), max_iterator);
            auto max_length = *max_iterator;

            output[0][i * ncol + j] = max_length * cellsize;
            output[1][i * ncol + j] = 180 * angle * max_position / M_PI;
            output[2][i * ncol + j] = max_length > 0 ? d1[max_position] / max_length : 1;
            output[3][i * ncol + j] = w12[max_position] / n12[max_position];
            output[4][i * ncol + j] = (std::accumulate(cos1.begin(), cos1.end(), 0.0) + std::accumulate(cos2.begin(), cos2.end(), 0.0)) / ndirs;

        }
    }

    return output;
}

// Region width calculation
std::vector<std::vector<double>> estimate_width(py::array_t<double> ceuc, py::array_t<double> cheight,
                                                double cellsize, double nodata) {


    auto buf_euc = ceuc.request(), buf_height = cheight.request();

    if (buf_euc.size != buf_height.size)
        throw std::runtime_error("Input shapes must match");

    auto nrow = buf_euc.shape[0];
    auto ncol = buf_euc.shape[1];

    auto *euc = (double *) buf_euc.ptr,
         *height = (double *) buf_height.ptr;

    std::vector<std::vector<double>> output(4, std::vector<double>(nrow * ncol));

    for (int i = 0; i < nrow; i++) {

        for (int j = 0; j < ncol; j++) {

            auto radius = euc[i * ncol + j];

            if (radius == nodata) {
                output[0][i * ncol + j] = nodata;
                output[1][i * ncol + j] = nodata;
                output[2][i * ncol + j] = nodata;
                output[3][i * ncol + j] = nodata;
                continue;
            }

            auto w = int(ceil(radius / cellsize));

            std::vector<int> covered_idx;
            double height_sum = 0;
            double height_max = 0;

            for (int k = -w; k <= w; ++k) {
                for (int l = -w; l <= w; ++l) {

                    if (k * k + l * l > w * w)
                        continue;

                    auto ik = i + k;
                    auto jl = j + l;

                    if (ik < 0 || ik >= nrow || jl < 0 || jl >= ncol)
                        continue;

                    auto idx = ik * ncol + jl;

                    if (output[0][idx] < 2 * radius) {
                        output[0][idx] = 2 * radius;
                        covered_idx.push_back(idx);
                        height_sum += height[idx];
                    }

                    height_max = (height[idx] > height_max) ? height[idx] : height_max;

                    if (covered_idx.size() > 0) {
                        auto mean_height = height_sum / covered_idx.size();
                        for (auto cidx: covered_idx) {
                            output[1][cidx] = height_max;
                            output[2][cidx] = mean_height;
                            output[3][cidx] = 0.5 * mean_height / radius;
                        }
                    }
                }
            }
        }
    }

    return output;
}

py::array_t<double> estimate_space(py::array_t<double> bld_euc, py::array_t<double> bld_alloc,
                                   py::array_t<double> veg_euc, py::array_t<double> veg_alloc,
                                   double cellsize, double nodata,
                                   int ndirs = 360, double radius = 0) {

    auto bld_euc_buf = bld_euc.request(),
         veg_euc_buf = bld_alloc.request(),
         bld_alloc_buf = bld_alloc.request(),
         veg_alloc_buf = veg_alloc.request();

    if (bld_euc_buf.size != veg_euc_buf.size)
        throw std::runtime_error("Input shapes must match");

    auto nrow = bld_euc_buf.shape[0];
    auto ncol = bld_euc_buf.shape[1];

    auto nelem = nrow * ncol;

    auto *beuc = (double *) bld_euc_buf.ptr,
         *bheight = (double *) bld_alloc_buf.ptr,
         *veuc = (double *) veg_euc_buf.ptr,
         *vheight = (double *) veg_alloc_buf.ptr;

    // out array that will store all desired information
    py::array_t<double> output = py::array_t<double>(std::vector<ptrdiff_t>{8, nrow, ncol});
    auto outbuf = output.request();
    auto* out = (double *) outbuf.ptr;

    for (int i = 0; i < nrow; i++) {

        for (int j = 0; j < ncol; j++) {

            auto radius = beuc[i * ncol + j];

            if (radius == nodata) {
                out[0 * nelem + i * ncol + j] = nodata;
                out[1 * nelem + i * ncol + j] = nodata;
                out[2 * nelem + i * ncol + j] = nodata;
                continue;
            }

            auto w = int(ceil(radius / cellsize));

            std::vector<int> covered_idx;
            double height_sum = 0;

            for (int k = -w; k <= w; ++k) {
                for (int l = -w; l <= w; ++l) {

                    if (k * k + l * l > w * w)
                        continue;

                    auto ik = i + k;
                    auto jl = j + l;

                    if (ik < 0 || ik >= nrow || jl < 0 || jl >= ncol)
                        continue;

                    auto idx = ik * ncol + jl;

                    if (out[0 * nelem + idx] < 2 * radius) {
                        out[0 * nelem + idx] = 2 * radius;
                        covered_idx.push_back(idx);
                        height_sum += bheight[idx];
                    }

                    if (height_sum > 0) {
                        auto mean_height = height_sum / covered_idx.size();
                        for (auto cidx: covered_idx) {
                            out[1 * nelem + cidx] = mean_height;
                            out[2 * nelem + cidx] = 0.5 * mean_height / radius;
                        }
                    }
                }
            }
        }
    }

//    return out;

    std::vector<Shifts> shifts;
    auto r = int(ceil(radius/cellsize));
    auto angle = 2 * M_PI / ndirs;
    double a = M_PI_2;
    int i, j;

    for (auto k = 0; k < ndirs; k++) {
        i = -r * cos(a);
        j = r * sin(a);
        shifts.push_back(get_shifts(0, 0, j, i));
        a -= angle;
    }

//    const auto& width = out[0];

    int ndirs2 = ndirs / 2;

    int ik, jk;

    // distances
    std::vector<double> d1(ndirs2); // reusable vector;
    std::vector<double> d2(ndirs2); // reusable vector;
    std::vector<double> d12(ndirs2); // reusable vector;

    // obstacle heights
    std::vector<double> cos1(ndirs2); // reusable vector;
    std::vector<double> cos2(ndirs2); // reusable vector;

    // widths and counts
    std::vector<double> w12(ndirs2); // reusable vector;
    std::vector<double> n12(ndirs2); // reusable vector;

    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            for (auto k = 0; k < ndirs2; k++) {
                n12[k] = 1;
                w12[k] = out[i * ncol + j];

                d1[k] = r;
                for (auto s: shifts[k]) {
                    ik = i + s.first;
                    jk = j + s.second;
                    if (is_within(ik, jk, nrow, ncol)) {
                        if (beuc[ik * ncol + jk] < cellsize) {
                            d1[k] = dist(s.first, s.second);
                            cos1[k] = ((ik == i) and (jk == j)) ? 1 : pow(cos(atan2(bheight[ik * ncol + jk], d1[k] * cellsize)), 2);
                            break;
                        } else {
                            w12[k] += out[ik * ncol + jk];
                            n12[k]++;
                        }
                    } else {
                        d1[k] = dist(s.first, s.second);
                        cos1[k] = pow(cos(atan2(bheight[ik * ncol + jk], d1[k] * cellsize)), 2);
                        break;
                    }

                }

                d2[k] = r;
                for (auto s: shifts[ndirs2 + k]) {
                    ik = i + s.first;
                    jk = j + s.second;
                    if (is_within(ik, jk, nrow, ncol)) {
                        if (beuc[ik * ncol + jk] < cellsize) {
                            d2[k] = dist(s.first, s.second);
                            cos2[k] = ((ik == i) and (jk == j)) ? 1 : pow(cos(atan2(bheight[ik * ncol + jk], d2[k] * cellsize)), 2);
                            break;
                        } else {
                            w12[k] += out[ik * ncol + jk];
                            n12[k]++;
                        }
                    } else {
                        d2[k] = dist(s.first, s.second);
                        cos2[k] = pow(cos(atan2(bheight[ik * ncol + jk], d2[k] * cellsize)), 2);
                        break;
                    }
                }

                d12[k] = d1[k] + d2[k];
            }

            auto max_iterator = std::max_element(d12.begin(), d12.end());
            auto max_position = std::distance(d12.begin(), max_iterator);
            auto max_length = *max_iterator;

            out[3 * nelem + i * ncol + j] = max_length * cellsize;
            out[4 * nelem + i * ncol + j] = 180 * angle * max_position / M_PI;
            out[5 * nelem + i * ncol + j] = max_length > 0 ? d1[max_position] / max_length : 1;
            out[6 * nelem + i * ncol + j] = w12[max_position] / n12[max_position];
            out[7 * nelem + i * ncol + j] = (std::accumulate(cos1.begin(), cos1.end(), 0.0) + std::accumulate(cos2.begin(), cos2.end(), 0.0)) / ndirs;
        }
    }

    return output;
}

void estimate_space2(py::array_t<double> ceuc, py::array_t<double> ceuc_full, py::array_t<double> cheight, py::array_t<double> cheight_full,
                     py::array_t<double> cwidth, py::array_t<double> chgt, py::array_t<double> chw, py::array_t<double> cbuilt,
                     py::array_t<double> clength, py::array_t<double> cdir, py::array_t<double> cpos, py::array_t<double> cdom,
                     py::array_t<double> cawidth, py::array_t<double> cahgt, py::array_t<double> cahw, py::array_t<double> cabuilt,
                     py::array_t<double> csvf,
                     double cellsize, double nodata, int ndirs, double search_radius, double min_hw = 0.1) {

//    std::ofstream myfile;
//    myfile.open ("/Volumes/Work/__UCLIM/Kosheleva/log.txt");

    auto buf_euc = ceuc.request(), buf_euc_full = ceuc_full.request(), buf_height = cheight.request(), buf_height_full = cheight_full.request(),
         buf_width = cwidth.request(), buf_hgt = chgt.request(), buf_hw = chw.request(), buf_built = cbuilt.request(),
         buf_length = clength.request(), buf_dir = cdir.request(), buf_pos = cpos.request(), buf_dom = cdom.request(),
         buf_awidth = cawidth.request(), buf_ahgt = cahgt.request(), buf_ahw = cahw.request(), buf_abuilt = cabuilt.request(),
         buf_svf = csvf.request();

    if (buf_euc.size != buf_height.size)
        throw std::runtime_error("Input shapes must match");

    auto nrow = buf_euc.shape[0];
    auto ncol = buf_euc.shape[1];
    auto nelem = nrow * ncol;

    auto *euc = (double *) buf_euc.ptr,
         *euc_full = (double *) buf_euc_full.ptr,
         *height = (double *) buf_height.ptr,
         *height_full = (double *) buf_height_full.ptr,
         *width = (double *) buf_width.ptr,
         *hgt = (double *) buf_hgt.ptr,
         *hw = (double *) buf_hw.ptr,
         *built = (double *) buf_built.ptr,
         *length = (double *) buf_length.ptr,
         *dir = (double *) buf_dir.ptr,
         *pos = (double *) buf_pos.ptr,
         *dom = (double *) buf_dom.ptr,
         *awidth = (double *) buf_awidth.ptr,
         *ahgt = (double *) buf_ahgt.ptr,
         *ahw = (double *) buf_ahw.ptr,
         *abuilt = (double *) buf_abuilt.ptr,
         *svf = (double *) buf_svf.ptr;

    int idx, i, j, ik, jk, jl, n;
    double radius, w, height_sum, built_sum, mean_height, mean_hw, min_hw_val;
    double *height_ptr, *euc_ptr;

    auto cmp = [](std::pair<int,double> left, std::pair<int,double> right) { return left.second < right.second; };
    std::priority_queue<std::pair<int,double>, std::vector<std::pair<int,double>>, decltype(cmp)> q(cmp);

    std::vector<bool> valid(nelem, true);

    std::vector<double> euc_derived(nelem, 0);
    std::vector<double> height_derived(nelem, 0);

    for (auto built_type: {1, 0}) {

        if (built_type == 1) {
            height_ptr = height;
            euc_ptr = euc;
            min_hw_val = min_hw;
        } else {
            height_ptr = height_full;
            euc_ptr = euc_full;
            min_hw_val = 0;
            std::fill(valid.begin(), valid.end(), true);
        }

        for (auto i = 0; i < nelem; ++i) {
            if (euc_ptr[i] > 0)
                q.emplace(std::pair<int, double>(i, euc_ptr[i]));
        }

        while (!q.empty()) {

            auto cell = q.top();
            q.pop();

//            if (width[cell.first] != nodata)
//                continue;

            idx = cell.first;
            radius = cell.second;

            i = floor(idx / ncol);
            j = idx - i * ncol;

            w = int(ceil(radius / cellsize));

            std::vector<int> covered_idx;
            height_sum = 0;
            built_sum = 0;
            n = 0;

            for (int k = -w; k <= w; ++k) {
                for (int l = -w; l <= w; ++l) {

                    if (k * k + l * l > w * w)
                        continue;

                    ik = i + k;
                    jl = j + l;

                    if (!is_within(ik, jl, nrow, ncol))
                        continue;

                    auto ikdx = ik * ncol + jl;

//                    if ((width[ikdx] < 2 * radius) && valid[ikdx]) {
                        covered_idx.push_back(ikdx);
//                    }

                    height_sum += height_ptr[ikdx];
                    n++;

                }
            }

            if (height_sum > 0) {
                mean_height = height_sum / n;
                mean_hw = 0.5 * mean_height / radius;
                if (mean_hw >= min_hw_val) {
                    for (auto cidx: covered_idx) {
                        if ((width[cidx] < 2 * radius) && valid[cidx]) {
                            width[cidx] = 2 * radius;
                            dom[cidx] = idx;
                            hgt[cidx] = mean_height;
                            hw[cidx] = mean_hw;
                            built[cidx] = built_type;
                            euc_derived[cidx] = euc_ptr[cidx];
                            height_derived[cidx] = height_ptr[cidx];
                        }
                    } 
                } else {
                    for (auto cidx: covered_idx) {
                        if (width[cidx] < 2 * radius) {
                            valid[cidx] = false;
                        }
                    }
                }
            }
        }
    }

    std::vector<Shifts> shifts;
    auto r = int(ceil(search_radius/cellsize));
    auto angle = 2 * M_PI / ndirs;
    double a = M_PI_2;

    for (auto k = 0; k < ndirs; k++) {
        i = -r * cos(a);
        j = r * sin(a);
        shifts.push_back(get_shifts(0, 0, j, i));
        a -= angle;
    }

    int ndirs2 = ndirs / 2;


    // distances
    std::vector<double> d1(ndirs2); // reusable vector;
    std::vector<double> d2(ndirs2); // reusable vector;
    std::vector<double> d12(ndirs2); // reusable vector;

    // obstacle heights
    std::vector<double> cos1(ndirs2); // reusable vector;
    std::vector<double> cos2(ndirs2); // reusable vector;

    // widths and counts
    std::vector<double> w12(ndirs2); // reusable vector;
    std::vector<double> h12(ndirs2); // reusable vector;
    std::vector<double> hw12(ndirs2); // reusable vector;
    std::vector<double> bld12(ndirs2); // reusable vector;
    std::vector<double> n12(ndirs2); // reusable vector;

    for (int i = 0; i < nrow; i++) {
//        myfile << i << std::endl;
        for (int j = 0; j < ncol; j++) {
//            myfile << j << std::endl;
            idx = i * ncol + j;
            if (width[idx] > 0) {
                for (auto k = 0; k < ndirs2; k++) {
                    n12[k] = 1;
                    w12[k] = width[idx];
                    h12[k] = hgt[idx];
                    hw12[k] = hw[idx];
                    bld12[k] = built[idx];
    //                    myfile << "w12[k] = " << w12[k] << std::endl;
                    d1[k] = r;
                    for (auto s: shifts[k]) {
                        ik = i + s.first;
                        jk = j + s.second;
                        if (is_within(ik, jk, nrow, ncol)) {
                            if (euc_derived[ik * ncol + jk] < cellsize) {
                                d1[k] = dist(s.first, s.second);
                                cos1[k] = ((ik == i) and (jk == j)) ? 1 : pow(cos(atan2(height_derived[ik * ncol + jk], d1[k] * cellsize)), 2);
    //                                myfile << "cos1[k] = " << cos1[k] << std::endl;
                                break;
                            } else {
                                w12[k] += width[ik * ncol + jk];
                                h12[k] += hgt[ik * ncol + jk];
                                hw12[k] += hw[ik * ncol + jk];
                                bld12[k] += built[ik * ncol + jk];
    //                                myfile << "plus1 w12[k] = " << w12[k] << std::endl;
                                n12[k]++;
                            }
                        } else {
                            d1[k] = dist(s.first, s.second);
                            cos1[k] = 1;
    //                            cos1[k] = pow(cos(atan2(height[ik * ncol + jk], d1[k] * cellsize)), 2);
    //                            myfile << "alt1 cos1[k] = " << cos1[k] << std::endl;
                            break;
                        }

                    }

                    d2[k] = r;
                    for (auto s: shifts[ndirs2 + k]) {
                        ik = i + s.first;
                        jk = j + s.second;
                        if (is_within(ik, jk, nrow, ncol)) {
                            if (euc_derived[ik * ncol + jk] < cellsize) {
                                d2[k] = dist(s.first, s.second);
                                cos2[k] = ((ik == i) and (jk == j)) ? 1 : pow(cos(atan2(height_derived[ik * ncol + jk], d2[k] * cellsize)), 2);
    //                                myfile << "cos2[k] = " << cos2[k] << std::endl;
                                break;
                            } else {
                                w12[k] += width[ik * ncol + jk];
                                h12[k] += hgt[ik * ncol + jk];
                                hw12[k] += hw[ik * ncol + jk];
                                bld12[k] += built[ik * ncol + jk];
    //                                myfile << "plus2 w12[k] = " << w12[k] << std::endl;
                                n12[k]++;
                            }
                        } else {
                            d2[k] = dist(s.first, s.second);
                            cos1[k] = 1;
    //                            cos2[k] = pow(cos(atan2(height[ik * ncol + jk], d2[k] * cellsize)), 2);
    //                            myfile << "alt2 cos2[k] = " << cos2[k] << std::endl;
                            break;
                        }
                    }

                    d12[k] = d1[k] + d2[k];
    //                    myfile << "sum d12[k] = " << d12[k] << std::endl;

                }

                //                myfile << "iterators start" << std::endl;

                auto max_iterator = std::max_element(d12.begin(), d12.end());
                auto max_position = std::distance(d12.begin(), max_iterator);
                auto max_length = *max_iterator;

                length[idx] = max_length * cellsize;
                dir[idx] = 180 * angle * max_position / M_PI;
                pos[idx] = max_length > 0 ? d1[max_position] / max_length : 1;
                awidth[idx] = w12[max_position] / n12[max_position];
                ahgt[idx] = h12[max_position] / n12[max_position];
                ahw[idx] = hw12[max_position] / n12[max_position];
                abuilt[idx] = bld12[max_position] / n12[max_position];
                svf[idx] = (std::accumulate(cos1.begin(), cos1.end(), 0.0) + std::accumulate(cos2.begin(), cos2.end(), 0.0)) / ndirs;
    //                myfile << "iterators stop" << std::endl;
            }
        }
    }

//    myfile << "Ended computations";
//    myfile.close();

    return;
}

PYBIND11_MODULE(rspace, m) {
    m.doc() = R"pbdoc(
        C++ plugin for region width estimation
        -----------------------

        .. currentmodule:: rspace

        .. autosummary::
           :toctree: _generate

           estimate_width, estimate_length, estimate_space

    )pbdoc";

    m.def("estimate_width", &estimate_width, R"pbdoc(
        Estimate local region width
    )pbdoc");

    m.def("estimate_length", &estimate_length, R"pbdoc(
        Estimate local region length
    )pbdoc");

    m.def("estimate_space", &estimate_space, R"pbdoc(
        Estimate local region space
    )pbdoc");

    m.def("estimate_space2", &estimate_space2, R"pbdoc(
        Estimate local region space
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
