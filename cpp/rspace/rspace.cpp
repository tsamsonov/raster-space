// [[Rcpp::plugins(cpp14)]]

#include<atomic>
#include<vector>
#include<cmath>
#include<iostream>
#include<fstream>
#include<queue>
#include<chrono>
#include<thread>
#include<functional>
#include<pybind11/pybind11.h>
#include<pybind11/numpy.h>
#include<pybind11/stl.h>
#include<pybind11/functional.h>

namespace py = pybind11;
using Shifts = std::vector<std::pair<int, int>>;

// Order of variables in resulting buffer
const int WIDTH = 0, HGT = 1, HW = 2, BUILT = 3, OFFSET = 4, SVF = 5, DOM = 6,
          LENGTH = 7, DIR = 8, POS = 9, AWIDTH = 10, AHGT = 11, AHW = 12, ABUILT = 13, AOFFSET = 14, ASVF = 15;

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
//    return get_shifts_byX(0, 0, 0, 0);
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

void calculate_length_params(double *out, const double *euc, const std::vector<Shifts>& shifts, const std::vector<int>& main_dir,
                             double cellsize, int nrow, int ncol, int first, int last) {
    int nelem = nrow * ncol;
    int i, j, ik, jk, ikdx;

    auto ndirs = shifts.size();
    auto ndirs2 = shifts.size() / 2;
    auto angle = 2 * M_PI / ndirs;

    const int WIDTH = 0, HGT = 1, HW = 2, BUILT = 3, OFFSET = 4, SVF = 5, DOM = 6,
            LENGTH = 7, DIR = 8, POS = 9, AWIDTH = 10, AHGT = 11, AHW = 12, ABUILT = 13, AOFFSET = 14, ASVF = 15;

    // offsets
    long width = WIDTH*nelem, hgt = HGT*nelem, hw = HW*nelem, built = BUILT*nelem, offset = OFFSET*nelem,
         svf = SVF*nelem, dom = DOM*nelem, length = LENGTH*nelem, dir = DIR*nelem, pos = POS*nelem,
         awidth = AWIDTH*nelem, ahgt = AHGT*nelem, ahw = AHW*nelem, abuilt = ABUILT*nelem, aoffset = AOFFSET*nelem, asvf = ASVF*nelem;

    for (int idx = first; idx < last; idx++) {
        i = idx / ncol;
        j = idx % ncol;

        if (euc[idx] > 0) {

            double d[] {0, 0};

            double w12 = 0;
            double h12 = 0;
            double hw12 = 0;
            double bld12 = 0;
            double n12 = 0;
            double svf12 = 0;
            double off12 = 0;

            for (auto k: {0, 1}) {
                for (auto s: shifts[main_dir[idx] + k*ndirs2]) {
                    ik = i + s.first;
                    jk = j + s.second;
                    if (is_within(ik, jk, nrow, ncol)) {
                        ikdx = ik * ncol + jk;
                        if (euc[ikdx] == 0) {
                            d[k] = dist(s.first, s.second);
                            break;
                        } else {
                            w12 += out[width + ikdx];
                            h12 += out[hgt + ikdx];
                            hw12 += out[hw + ikdx];
                            bld12 += out[built + ikdx];
                            off12 += out[offset + ikdx];
                            svf12 += out[svf + ikdx];
                            n12++;
                        }
                    } else {
                        d[k] = dist(s.first, s.second);
                        break;
                    }
                }
            }

            auto d01 = d[0] + d[1];

            out[length + idx] = d01 * cellsize;
            out[dir + idx] = 180 * angle * main_dir[idx] / M_PI;
            out[pos + idx] = d01 > 0 ? d[0] / d01 : 1;
            out[awidth + idx] = w12 / n12;
            out[ahgt + idx] = h12 / n12;
            out[ahw + idx] = hw12 / n12;
            out[abuilt + idx] = bld12 / n12;
            out[aoffset + idx] = off12 / n12;
            out[asvf + idx] = svf12 / n12;
        }
    }
}

void find_main_dir(double cellsize, int dir1, int dir2, long nrow, long ncol, const double *euc, const double *hgt,
                   const std::vector<Shifts>& shifts, std::vector<int> &main_dir, std::vector<double> &main_len,
                   std::vector<double> &cossum, std::ofstream &log) {

    auto ndirs = shifts.size();
    auto ndirs2 = ndirs / 2;
    auto nelem = nrow * ncol;
    std::vector<double> len(nelem, 0);
    std::vector<bool> newdir(nelem, true);

    for (auto dk = dir1; dk < dir2; dk++) {
        auto sh = std::vector<Shifts> {shifts[dk], shifts[ndirs2 + dk]};
        std::vector<int> sidx {0, 0};
        std::vector<double> h {0, 0};

        for (auto idx = 0; idx < nelem; idx++) {
            if (len[idx] == 0 and euc[idx] > 0) {
                auto i = idx / ncol;
                auto j = idx % ncol;
                sidx[0] = 0;
                sidx[1] = 0;
                h[0] = 0;
                h[1] = 0;

                for (auto k: {0, 1}){
                    auto sk = 0;
                    for (auto s: sh[k]) {
                        auto ik = i + s.first;
                        auto jk = j + s.second;
                        if (is_within(ik, jk, nrow, ncol)) {
                            auto ikdx = ik * ncol + jk;
                            if (euc[ikdx] == 0) {
                                sidx[k] = sk;
                                h[k] = hgt[ikdx];
                                break;
                            }
                        } else {
                            sidx[k] = sk;
                            h[k] = 0;
                            break;
                        }
                        sk++;
                    }
                }

                double d = dist(sh[0][sidx[0]].first, sh[0][sidx[0]].second) +
                           dist(sh[1][sidx[1]].first, sh[1][sidx[1]].second);

                for (auto k: {0, 1}){
                    auto s = sh[k];
                    for (auto sk = 0; sk < sidx[k]; sk++) {
                        auto ik = i + s[sk].first;
                        auto jk = j + s[sk].second;
                        auto ikdx = ik * ncol + jk;
                        len[ikdx] = d;

                        if (newdir[ikdx]) {
                            cossum[ikdx] += pow(cos(atan2(h[k], dist(s[sidx[k]].first - s[sk].first, s[sidx[k]].second - s[sk].second) * cellsize)), 2);
                            newdir[ikdx] = false;
                        }

//                    cossum[ikdx] += 1;
                    }
                }
            }
        }

        for (auto i = 0; i < nelem; ++i) {
            if (len[i] > main_len[i]) {
                main_len[i] = len[i];
                main_dir[i] = dk;
            }
        }

        std::fill(len.begin(), len.end(), 0);
        std::fill(newdir.begin(), newdir.end(), true);
    }
}

py::array_t<double> estimate_space(py::array_t<double> bld_euc, py::array_t<double> bld_alloc,
                                   py::array_t<double> veg_euc, py::array_t<double> veg_alloc,
                                   double cellsize, double nodata,
                                   unsigned int ndirs, double search_radius, unsigned int nthreads,
                                   py::function &feedback, py::function &progress) {

    std::ofstream log;
    std::remove("/Volumes/Work/__UCLIM/Kosheleva/log.txt"); // delete file
    log.open("/Volumes/Work/__UCLIM/Kosheleva/log.txt");

    auto bld_euc_buf = bld_euc.request(),
         bld_alloc_buf = bld_alloc.request(),
         veg_euc_buf = veg_euc.request(),
         veg_alloc_buf = veg_alloc.request();

    if (bld_euc_buf.size != veg_euc_buf.size)
        throw std::runtime_error("Input shapes must match");

    auto nrow = bld_euc_buf.shape[0];
    auto ncol = bld_euc_buf.shape[1];
    auto nelem = nrow * ncol;
    auto frac = nelem / 100;

    auto *beuc = (double *) bld_euc_buf.ptr,
         *bheight = (double *) bld_alloc_buf.ptr,
         *veuc = (double *) veg_euc_buf.ptr,
         *vheight = (double *) veg_alloc_buf.ptr;

    // out array that will store all desired information
    py::array_t<double> output = py::array_t<double>(std::vector<ptrdiff_t>{16, nrow, ncol});

    auto outbuf = output.request();
    auto *out = (double *) outbuf.ptr;

    // offsets
    long width = WIDTH*nelem, hgt = HGT*nelem, hw = HW*nelem, built = BUILT*nelem, offset = OFFSET*nelem,
         svf = SVF*nelem, dom = DOM*nelem, length = LENGTH*nelem, dir = DIR*nelem, pos = POS*nelem,
         awidth = AWIDTH*nelem, ahgt = AHGT*nelem, ahw = AHW*nelem, abuilt = ABUILT*nelem, aoffset = AOFFSET*nelem, asvf = ASVF*nelem;

    int idx, ikdx, i, j, ik, jk, jl, n;
    double radius, w, height_sum, built_sum, mean_hw, mean_height;
    double *height_ptr, *euc_ptr;

    euc_ptr = beuc;
    height_ptr = bheight;

    // the pixels are considered in the order of decreasing distance
    auto cmp = [](std::pair<int,double> left, std::pair<int,double> right) { return left.second < right.second; };
    std::priority_queue<std::pair<int,double>, std::vector<std::pair<int,double>>, decltype(cmp)> q(cmp);

    std::vector<bool> valid(nelem, true);
    std::vector<double> euc_derived(nelem, 0);
    std::vector<double> height_derived(nelem, 0);

    for (auto i = 0; i < nelem; ++i) {
        out[offset + i] = euc_ptr[i];
        if (euc_ptr[i] > 0)
            q.emplace(std::pair<int, double>(i, euc_ptr[i]));
    }

    feedback("CALCULATING WIDTH-BASED PARAMETERS...");

    // WIDTH COMPUTATION
    while (!q.empty()) {

        auto cell = q.top();
        q.pop();

        idx = cell.first;
        radius = cell.second;

//        if (out[width + idx] > 0)
//            continue;

        i = floor(idx / ncol);
        j = idx - i * ncol;

        w = int(floor(radius / cellsize));

        height_sum = 0;
        built_sum = 0;
        n = 0;

        std::vector<int> covered_idx;

        for (int k = -w; k <= w; ++k) {
            for (int l = -w; l <= w; ++l) {

                if (k * k + l * l >= w * w)
                    continue;

                ik = i + k;
                jl = j + l;

                if (!is_within(ik, jl, nrow, ncol))
                    continue;

                ikdx = ik * ncol + jl;

                if (out[width + ikdx] < 2 * radius)
                    covered_idx.push_back(ikdx);

                height_sum += height_ptr[ikdx];
                n++;
            }
        }

        if (!covered_idx.empty()) {
            mean_height = height_sum / n;
            mean_hw = 0.5 * mean_height / radius;
            for (auto cidx: covered_idx) {
                out[width + cidx] = 2 * radius;
                out[dom + cidx] = idx;
                out[hgt + cidx] = mean_height;
                out[hw + cidx] = mean_hw;
                out[built + cidx] = 1;
                euc_derived[cidx] = euc_ptr[cidx];
                height_derived[cidx] = height_ptr[cidx];
            }
        }

//        std::vector<int>().swap(covered_idx);

//        cur++;
//        if (cur % frac == 0) {
//            progress(100 * cur / nelem);
//        }

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


    auto num_threads = nthreads > 0 ? std::min(nthreads, ndirs/2) : std::min(std::thread::hardware_concurrency(), ndirs/2);

    std::vector<std::thread> threads;
    threads.reserve(num_threads);

    auto th_size = ndirs / (2 * num_threads);

    feedback("FINDING MAIN DIRECTIONS...");

    std::vector<int> main_dir(nelem, 0);
    std::vector<double> main_len(nelem, 0.0);
    std::vector<double> cossum(nelem, 0.0);

    if (num_threads > 1) {

        feedback("Using " + std::to_string(num_threads) + " threads with " + std::to_string(th_size) + " directions each");

        std::vector<std::vector<int>> main_dirs(num_threads, std::vector<int>(nelem, 0));
        std::vector<std::vector<double>> main_lens(num_threads, std::vector<double>(nelem, 0.0));
        std::vector<std::vector<double>> cossums(num_threads, std::vector<double>(nelem, 0.0));

        i = 0;
        while (i < num_threads-1) {
            threads.emplace_back(std::thread(find_main_dir, cellsize, th_size * i, th_size * (i + 1),
                                             nrow, ncol, euc_ptr, height_ptr, std::ref(shifts), std::ref(main_dirs[i]),
                                             std::ref(main_lens[i]), std::ref(cossums[i]), std::ref(log)));
            i++;
        }

        // last thread
        threads.emplace_back(std::thread(find_main_dir, cellsize, th_size * i, ndirs / 2,
                                         nrow, ncol, euc_ptr, height_ptr, std::ref(shifts), std::ref(main_dirs[i]),
                                         std::ref(main_lens[i]), std::ref(cossums[i]), std::ref(log)));

        for (auto &th : threads) {
            th.join();
        }

        for (auto i = 0; i < num_threads; ++i) {
            for (auto j = 0; j < nelem; ++j) {
                if (main_lens[i][j] > main_len[j]) {
                    main_len[j] = main_lens[i][j];
                    main_dir[j] = main_dirs[i][j];
                }
                cossum[j] += cossums[i][j];
            }
        }

        threads.clear();

    } else {
        try {
            find_main_dir(cellsize, 0, ndirs / 2, nrow, ncol, euc_ptr, height_ptr, shifts, main_dir, main_len, cossum, log);
        } catch (std::exception &e) {
            py::print(e.what());
        }
    }

    // calculate svf from sum
    for (int k = 0; k < cossum.size(); ++k) {
        out[svf + k] = 2 * cossum[k] / ndirs;
    }

    feedback("CALCULATING LENGTH-BASED PARAMETERS...");

    num_threads = nthreads > 0 ? nthreads : std::thread::hardware_concurrency();

    if (num_threads > 1) {
        th_size = nelem / num_threads;
        feedback("Using " + std::to_string(num_threads) + " threads with " + std::to_string(th_size) + " pixels each");

        i = 0;
        while (i < num_threads-1) {
            threads.emplace_back(
                    std::thread(calculate_length_params, out, euc_ptr,
                                std::ref(shifts), std::ref(main_dir),
                                cellsize,  nrow, ncol, th_size * i, th_size * (i + 1)));
            i++;
        }

        // last thread
        threads.emplace_back(
                std::thread(calculate_length_params, out, euc_ptr,
                            std::ref(shifts), std::ref(main_dir),
                            cellsize, nrow, ncol, th_size * i, nelem));

        for (auto &th : threads) {
            th.join();
        }
    } else {
        calculate_length_params(out, euc_ptr, shifts, main_dir, cellsize, nrow, ncol, 0, nelem);
    }

    log.close();

    return output;
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

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
