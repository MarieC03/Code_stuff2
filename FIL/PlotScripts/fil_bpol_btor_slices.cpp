// Portable note:
// This code should be compiled with g++ or icpx, not with legacy icpc/h5c++ wrappers
// that are tied to Intel classic + GCC13/glibc headers.
// Example:
//   g++ -O3 -std=c++17 fil_bpol_btor_slices_portable.cpp -I${HDF5_DIR}/include -L${HDF5_DIR}/lib -lhdf5 -o fil_bpol_btor_slices
// or
//   icpx -O3 -std=c++17 fil_bpol_btor_slices_portable.cpp -I${HDF5_DIR}/include -L${HDF5_DIR}/lib -lhdf5 -o fil_bpol_btor_slices

#if defined(__INTEL_COMPILER) && !defined(__INTEL_LLVM_COMPILER)
#error "This source should not be compiled with Intel classic icpc. Use g++ or icpx instead; the _Float32/_Float64 errors come from an icpc <-> GCC13/glibc incompatibility, not from this code."
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <hdf5.h>

namespace fs = std::filesystem;

static constexpr double RSCALE_KM = 1.47662504;
static constexpr double CU_TO_MS = 4.925490947e-3;

struct Options {
    fs::path root;
    fs::path outdir = "derived_bfield_slices_cpp";
    std::string plane = "xz";
    double coord0 = 0.0;
    bool mirror_z = false;
    int itmin = std::numeric_limits<int>::min();
    int itmax = std::numeric_limits<int>::max();
    std::optional<double> tmin_ms;
    std::optional<double> tmax_ms;
    std::optional<double> tmerg_code;
    std::optional<int> requested_rl;
    bool write_abs_btor = true;
};

struct DatasetRef {
    fs::path file;
    std::string dset_name;
    std::string group;
    std::string var;
    int iteration = -1;
    int tl = -1;
    int rl = -1;
    int c = -1;
    double time = std::numeric_limits<double>::quiet_NaN();
    std::array<double, 3> origin{};
    std::array<double, 3> delta{};
    std::array<int, 3> shape_xyz{}; // x,y,z order
    std::array<double, 3> corner{};
};

struct VarAliases {
    std::string canonical;
    std::vector<std::string> aliases;
};

struct Mesh2D {
    std::vector<double> x_km;
    std::vector<double> y_km;
    std::size_t nx = 0;
    std::size_t ny = 0;
};

struct Field2D {
    std::size_t nx = 0;
    std::size_t ny = 0;
    std::vector<float> a;
    Field2D() = default;
    Field2D(std::size_t ny_, std::size_t nx_)
        : nx(nx_), ny(ny_), a(nx_ * ny_, std::numeric_limits<float>::quiet_NaN()) {}
    float &operator()(std::size_t j, std::size_t i) { return a[j * nx + i]; }
    const float &operator()(std::size_t j, std::size_t i) const { return a[j * nx + i]; }
};

struct SliceFields {
    Field2D bx, by, bz;
    Field2D b0, b2, W;
    Field2D bpol, btor;
};

[[noreturn]] void die(const std::string &msg) {
    throw std::runtime_error(msg);
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) out.push_back(item);
    return out;
}

bool starts_with(const std::string &s, const std::string &prefix) {
    return s.rfind(prefix, 0) == 0;
}

bool ends_with(const std::string &s, const std::string &suffix) {
    return s.size() >= suffix.size() && s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

Options parse_args(int argc, char **argv) {
    Options opt;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto need = [&](const std::string &name) -> std::string {
            if (i + 1 >= argc) die("Missing value after " + name);
            return argv[++i];
        };
        if (a == "--root") opt.root = need(a);
        else if (a == "--outdir") opt.outdir = need(a);
        else if (a == "--plane") opt.plane = need(a);
        else if (a == "--coord0") opt.coord0 = std::stod(need(a));
        else if (a == "--mirror-z") opt.mirror_z = true;
        else if (a == "--itmin") opt.itmin = std::stoi(need(a));
        else if (a == "--itmax") opt.itmax = std::stoi(need(a));
        else if (a == "--tmin") opt.tmin_ms = std::stod(need(a));
        else if (a == "--tmax") opt.tmax_ms = std::stod(need(a));
        else if (a == "--tmerg-code") opt.tmerg_code = std::stod(need(a));
        else if (a == "--rl") opt.requested_rl = std::stoi(need(a));
        else if (a == "--signed-btor") opt.write_abs_btor = false;
        else if (a == "-h" || a == "--help") {
            std::cout << "Usage: " << argv[0] << " --root SIMDIR [options]\n"
                      << "Options:\n"
                      << "  --outdir DIR\n"
                      << "  --plane {xy,xz,yz}\n"
                      << "  --coord0 VALUE      fixed coordinate in code units for the slice\n"
                      << "  --mirror-z          reflect z>=0 data to z<0 before slicing in xz/yz\n"
                      << "  --itmin N --itmax N\n"
                      << "  --tmin MS --tmax MS  filter on t or (t-tmerg) in ms\n"
                      << "  --tmerg-code VALUE   merger time in code units\n"
                      << "  --rl N               use only one rl; default: use all rl with fine overwriting coarse\n"
                      << "  --signed-btor        keep sign of b_tor instead of abs(b_tor)\n";
            std::exit(0);
        } else die("Unknown argument: " + a);
    }
    if (opt.root.empty()) die("You must pass --root");
    if (!(opt.plane == "xy" || opt.plane == "xz" || opt.plane == "yz")) die("--plane must be xy, xz, or yz");
    return opt;
}

std::vector<fs::path> candidate_data_dirs(const fs::path &root) {
    std::vector<fs::path> dirs;
    if (!fs::exists(root)) die("Root does not exist: " + root.string());
    auto natural_key = [](const fs::path &p) {
        const std::string s = p.string();
        std::vector<std::string> out;
        std::string cur;
        bool in_digits = false;
        for (char ch : s) {
            bool d = std::isdigit(static_cast<unsigned char>(ch));
            if (cur.empty()) {
                cur.push_back(ch);
                in_digits = d;
            } else if (d == in_digits) {
                cur.push_back(ch);
            } else {
                out.push_back(cur);
                cur = std::string(1, ch);
                in_digits = d;
            }
        }
        if (!cur.empty()) out.push_back(cur);
        return out;
    };
    if (root.filename() == "data_hdf5_3D") {
        dirs.push_back(root);
        return dirs;
    }
    if (fs::is_directory(root)) {
        if (starts_with(root.filename().string(), "output")) {
            fs::path d = root / "data_hdf5_3D";
            if (fs::exists(d)) dirs.push_back(d);
            return dirs;
        }
        for (const auto &entry : fs::directory_iterator(root)) {
            if (!entry.is_directory()) continue;
            const auto name = entry.path().filename().string();
            if (!starts_with(name, "output")) continue;
            fs::path d = entry.path() / "data_hdf5_3D";
            if (fs::exists(d)) dirs.push_back(d);
        }
        std::sort(dirs.begin(), dirs.end(), [&](const fs::path &a, const fs::path &b) {
            return natural_key(a) < natural_key(b);
        });
    }
    return dirs;
}

std::vector<VarAliases> required_vars() {
    return {
        {"rho", {"rho", "rho_b"}},
        {"Bx", {"Bvec[0]", "Bx", "bvecx"}},
        {"By", {"Bvec[1]", "By", "bvecy"}},
        {"Bz", {"Bvec[2]", "Bz", "bvecz"}},
        {"vx", {"vel[0]", "vx"}},
        {"vy", {"vel[1]", "vy"}},
        {"vz", {"vel[2]", "vz"}},
        {"alp", {"alp", "alpha"}},
        {"betax", {"betax"}},
        {"betay", {"betay"}},
        {"betaz", {"betaz"}},
        {"gxx", {"gxx"}},
        {"gxy", {"gxy"}},
        {"gxz", {"gxz"}},
        {"gyy", {"gyy"}},
        {"gyz", {"gyz"}},
        {"gzz", {"gzz"}},
        {"W", {"w_lorentz", "W"}}, // optional in practice
    };
}

bool filename_matches_alias(const std::string &filename, const std::string &alias) {
    if (!ends_with(filename, ".h5")) return false;
    return starts_with(filename, alias + ".xyz") || filename == alias + ".h5";
}

std::optional<std::array<double, 3>> read_attr_double3(hid_t obj, const char *name) {
    if (H5Aexists(obj, name) <= 0) return std::nullopt;
    hid_t attr = H5Aopen(obj, name, H5P_DEFAULT);
    if (attr < 0) return std::nullopt;
    std::array<double, 3> val{};
    herr_t status = H5Aread(attr, H5T_NATIVE_DOUBLE, val.data());
    H5Aclose(attr);
    if (status < 0) return std::nullopt;
    return val;
}

std::optional<double> read_attr_double_scalar(hid_t obj, const char *name) {
    if (H5Aexists(obj, name) <= 0) return std::nullopt;
    hid_t attr = H5Aopen(obj, name, H5P_DEFAULT);
    if (attr < 0) return std::nullopt;
    double val = std::numeric_limits<double>::quiet_NaN();
    herr_t status = H5Aread(attr, H5T_NATIVE_DOUBLE, &val);
    H5Aclose(attr);
    if (status < 0) return std::nullopt;
    return val;
}

std::vector<std::string> list_root_datasets(hid_t file) {
    std::vector<std::string> names;
    H5G_info_t info{};
    if (H5Gget_info(file, &info) < 0) return names;
    for (hsize_t i = 0; i < info.nlinks; ++i) {
        ssize_t n = H5Lget_name_by_idx(file, ".", H5_INDEX_NAME, H5_ITER_INC, i, nullptr, 0, H5P_DEFAULT);
        if (n < 0) continue;
        std::string name(n + 1, '\0');
        H5Lget_name_by_idx(file, ".", H5_INDEX_NAME, H5_ITER_INC, i, name.data(), name.size(), H5P_DEFAULT);
        name.pop_back();
        names.push_back(name);
    }
    return names;
}

std::optional<DatasetRef> parse_dataset_ref(const fs::path &file_path, const std::string &dset_name) {
    static const std::regex re(R"((\S*)::(\S*) it=(\d+) tl=(\d+) rl=(\d+)(?: c=(\d+))?)");
    std::smatch m;
    if (!std::regex_match(dset_name, m, re)) return std::nullopt;

    hid_t file = H5Fopen(file_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) return std::nullopt;
    hid_t dset = H5Dopen2(file, dset_name.c_str(), H5P_DEFAULT);
    if (dset < 0) {
        H5Fclose(file);
        return std::nullopt;
    }
    hid_t space = H5Dget_space(dset);
    int ndims = H5Sget_simple_extent_ndims(space);
    if (ndims != 3) {
        H5Sclose(space); H5Dclose(dset); H5Fclose(file);
        return std::nullopt;
    }
    hsize_t dims[3]{};
    H5Sget_simple_extent_dims(space, dims, nullptr);

    auto origin_opt = read_attr_double3(dset, "origin");
    auto delta_opt = read_attr_double3(dset, "delta");
    auto time_opt = read_attr_double_scalar(dset, "time");
    if (!origin_opt || !delta_opt) {
        H5Sclose(space); H5Dclose(dset); H5Fclose(file);
        return std::nullopt;
    }

    DatasetRef ref;
    ref.file = file_path;
    ref.dset_name = dset_name;
    ref.group = m[1];
    ref.var = m[2];
    ref.iteration = std::stoi(m[3]);
    ref.tl = std::stoi(m[4]);
    ref.rl = std::stoi(m[5]);
    ref.c = m[6].matched ? std::stoi(m[6]) : -1;
    ref.time = time_opt.value_or(std::numeric_limits<double>::quiet_NaN());
    ref.origin = *origin_opt;
    ref.delta = *delta_opt;
    ref.shape_xyz = {static_cast<int>(dims[2]), static_cast<int>(dims[1]), static_cast<int>(dims[0])};
    for (int a = 0; a < 3; ++a) ref.corner[a] = ref.origin[a] + ref.delta[a] * (ref.shape_xyz[a] - 1);

    H5Sclose(space); H5Dclose(dset); H5Fclose(file);
    return ref;
}

using RLIndex = std::map<int, DatasetRef>; // keyed by component c
using VarIndex = std::map<std::string, std::map<int, std::map<int, RLIndex>>>;

VarIndex build_index(const Options &opt) {
    const auto dirs = candidate_data_dirs(opt.root);
    if (dirs.empty()) die("No output*/data_hdf5_3D directories found under " + opt.root.string());

    VarIndex index;
    const auto vars = required_vars();

    for (const auto &dir : dirs) {
        for (const auto &entry : fs::directory_iterator(dir)) {
            if (!entry.is_regular_file()) continue;
            const auto fname = entry.path().filename().string();
            std::optional<std::string> canonical;
            for (const auto &v : vars) {
                for (const auto &alias : v.aliases) {
                    if (filename_matches_alias(fname, alias)) {
                        canonical = v.canonical;
                        break;
                    }
                }
                if (canonical) break;
            }
            if (!canonical) continue;

            hid_t file = H5Fopen(entry.path().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
            if (file < 0) continue;
            for (const auto &dname : list_root_datasets(file)) {
                auto ref = parse_dataset_ref(entry.path(), dname);
                if (!ref) continue;
                if (ref->tl != 0) continue;
                const int comp = (ref->c >= 0 ? ref->c : 0);
                // Because candidate_data_dirs() is naturally sorted, later outputs overwrite earlier ones.
                index[*canonical][ref->iteration][ref->rl][comp] = *ref;
            }
            H5Fclose(file);
        }
    }

    for (const auto &v : vars) {
        if (v.canonical == "W") continue; // optional
        if (!index.count(v.canonical)) die("Missing required variable file(s) for " + v.canonical);
    }
    return index;
}

std::set<int> intersection_iterations(const VarIndex &index, bool include_W_optional = false) {
    std::set<int> common;
    bool first = true;
    for (const auto &kv : index) {
        if (kv.first == "W" && !include_W_optional) continue;
        std::set<int> its;
        for (const auto &it_kv : kv.second) its.insert(it_kv.first);
        if (first) {
            common = std::move(its);
            first = false;
        } else {
            std::set<int> tmp;
            std::set_intersection(common.begin(), common.end(), its.begin(), its.end(), std::inserter(tmp, tmp.begin()));
            common = std::move(tmp);
        }
    }
    return common;
}

std::map<int, double> iteration_to_time_ms(const VarIndex &index, const std::string &reference_var = "rho") {
    std::map<int, double> out;
    auto itv = index.find(reference_var);
    if (itv == index.end()) return out;
    for (const auto &[it, rl_map] : itv->second) {
        for (const auto &[rl, refs] : rl_map) {
            if (!refs.empty()) {
                out[it] = refs.begin()->second.time * CU_TO_MS;
                break;
            }
        }
    }
    return out;
}

bool pass_filters(const Options &opt, int it, double time_ms) {
    if (it < opt.itmin || it > opt.itmax) return false;
    double tplot = time_ms;
    if (opt.tmerg_code) tplot -= (*opt.tmerg_code) * CU_TO_MS;
    if (opt.tmin_ms && tplot < *opt.tmin_ms) return false;
    if (opt.tmax_ms && tplot > *opt.tmax_ms) return false;
    return true;
}

std::vector<int> choose_iterations(const Options &opt, const VarIndex &index) {
    const auto common = intersection_iterations(index, false);
    const auto times = iteration_to_time_ms(index);
    std::vector<int> out;
    for (int it : common) {
        auto tt = times.find(it);
        if (tt == times.end()) continue;
        if (pass_filters(opt, it, tt->second)) out.push_back(it);
    }
    if (out.empty()) die("No iterations remain after intersecting variables and applying filters.");
    return out;
}

std::vector<int> available_rls_for_iteration(const VarIndex &index, int iteration, bool include_W = false) {
    std::set<int> common;
    bool first = true;
    for (const auto &kv : index) {
        if (kv.first == "W" && !include_W) continue;
        auto it = kv.second.find(iteration);
        if (it == kv.second.end()) return {};
        std::set<int> rls;
        for (const auto &rlkv : it->second) rls.insert(rlkv.first);
        if (first) {
            common = std::move(rls);
            first = false;
        } else {
            std::set<int> tmp;
            std::set_intersection(common.begin(), common.end(), rls.begin(), rls.end(), std::inserter(tmp, tmp.begin()));
            common = std::move(tmp);
        }
    }
    return std::vector<int>(common.begin(), common.end());
}

std::vector<DatasetRef> collect_refs(const VarIndex &index, const std::string &var, int iteration, const Options &opt) {
    std::vector<DatasetRef> refs;
    auto iv = index.find(var);
    if (iv == index.end()) return refs;
    auto it = iv->second.find(iteration);
    if (it == iv->second.end()) return refs;
    for (const auto &[rl, comp_map] : it->second) {
        if (opt.requested_rl && rl != *opt.requested_rl) continue;
        for (const auto &[comp, ref] : comp_map) refs.push_back(ref);
    }
    std::sort(refs.begin(), refs.end(), [](const DatasetRef &a, const DatasetRef &b) {
        if (a.rl != b.rl) return a.rl < b.rl;
        if (a.c != b.c) return a.c < b.c;
        return a.file.string() < b.file.string();
    });
    return refs;
}

std::pair<int, std::pair<int,int>> choose_plane_axes(const std::string &plane) {
    if (plane == "xy") return {2, {0,1}};
    if (plane == "xz") return {1, {0,2}};
    if (plane == "yz") return {0, {1,2}};
    die("Unsupported plane");
}

Mesh2D build_mesh_for_iteration(const VarIndex &index, int iteration, const Options &opt) {
    auto plane_info = choose_plane_axes(opt.plane);
    const int dir_idx = plane_info.first;
    const int a0 = plane_info.second.first;
    const int a1 = plane_info.second.second;

    bool have = false;
    double xmin = 0.0, xmax = 0.0, ymin = 0.0, ymax = 0.0;
    double dx_min = std::numeric_limits<double>::max();
    double dy_min = std::numeric_limits<double>::max();

    for (const auto &kv : index) {
        if (kv.first == "W") continue;
        for (const auto &ref : collect_refs(index, kv.first, iteration, opt)) {
            double lo_dir = std::min(ref.origin[dir_idx], ref.corner[dir_idx]);
            double hi_dir = std::max(ref.origin[dir_idx], ref.corner[dir_idx]);
            const double tol = std::max(1e-12, 0.51 * std::abs(ref.delta[dir_idx]));
            bool intersects = (opt.coord0 >= lo_dir - tol && opt.coord0 <= hi_dir + tol);
            if (opt.mirror_z && (opt.plane == "xz" || opt.plane == "yz") && ref.origin[2] >= -1e-12 && ref.corner[2] > 0.0) {
                intersects = true;
            }
            if (!intersects) continue;

            double x0 = std::min(ref.origin[a0], ref.corner[a0]) * RSCALE_KM;
            double x1 = std::max(ref.origin[a0], ref.corner[a0]) * RSCALE_KM;
            double y0 = std::min(ref.origin[a1], ref.corner[a1]) * RSCALE_KM;
            double y1 = std::max(ref.origin[a1], ref.corner[a1]) * RSCALE_KM;
            if (!have) {
                xmin = x0; xmax = x1; ymin = y0; ymax = y1; have = true;
            } else {
                xmin = std::min(xmin, x0); xmax = std::max(xmax, x1);
                ymin = std::min(ymin, y0); ymax = std::max(ymax, y1);
            }
            dx_min = std::min(dx_min, std::abs(ref.delta[a0]) * RSCALE_KM);
            dy_min = std::min(dy_min, std::abs(ref.delta[a1]) * RSCALE_KM);
        }
    }
    if (!have) die("Could not determine a slice mesh for iteration " + std::to_string(iteration));

    auto make_axis = [](double lo, double hi, double d) {
        std::vector<double> v;
        const long n = static_cast<long>(std::llround((hi - lo) / d)) + 1;
        v.reserve(std::max<long>(n,1));
        for (long i = 0; i < n; ++i) v.push_back(lo + d * static_cast<double>(i));
        return v;
    };

    Mesh2D m;
    m.x_km = make_axis(xmin, xmax, dx_min);
    m.y_km = make_axis(ymin, ymax, dy_min);
    m.nx = m.x_km.size();
    m.ny = m.y_km.size();
    return m;
}

std::vector<float> read_full_dataset_float(const DatasetRef &ref) {
    hid_t file = H5Fopen(ref.file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) die("Failed to open file " + ref.file.string());
    hid_t dset = H5Dopen2(file, ref.dset_name.c_str(), H5P_DEFAULT);
    if (dset < 0) { H5Fclose(file); die("Failed to open dataset " + ref.dset_name); }
    std::size_t n = static_cast<std::size_t>(ref.shape_xyz[0]) * ref.shape_xyz[1] * ref.shape_xyz[2];
    std::vector<float> data(n);
    if (H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data()) < 0) {
        H5Dclose(dset); H5Fclose(file);
        die("Failed reading dataset " + ref.dset_name);
    }
    H5Dclose(dset); H5Fclose(file);
    return data;
}

inline std::size_t idx3_zyx(int iz, int iy, int ix, int ny, int nx) {
    return (static_cast<std::size_t>(iz) * ny + iy) * nx + ix;
}

struct TileSlice {
    std::vector<double> x_km;
    std::vector<double> y_km;
    std::vector<float> vals;
    int nx = 0;
    int ny = 0;
};

bool is_odd_under_mirror_z(const std::string &var) {
    static const std::set<std::string> odd = {"vel[2]", "vz", "betaz", "gxz", "gyz", "Bvec[0]", "Bvec[1]", "Bx", "By"};
    return odd.count(var) > 0;
}

TileSlice extract_tile_slice(const DatasetRef &ref, const std::vector<float> &data, const Options &opt) {
    auto plane_info = choose_plane_axes(opt.plane);
    const int dir_idx = plane_info.first;
    const int a0 = plane_info.second.first;
    const int a1 = plane_info.second.second;

    std::vector<double> x(ref.shape_xyz[0]), y(ref.shape_xyz[1]), z(ref.shape_xyz[2]);
    for (int i = 0; i < ref.shape_xyz[0]; ++i) x[i] = ref.origin[0] + ref.delta[0] * i;
    for (int j = 0; j < ref.shape_xyz[1]; ++j) y[j] = ref.origin[1] + ref.delta[1] * j;
    for (int k = 0; k < ref.shape_xyz[2]; ++k) z[k] = ref.origin[2] + ref.delta[2] * k;

    std::vector<float> ext = data;
    int nz_ext = ref.shape_xyz[2];
    std::vector<double> zext = z;
    if (opt.mirror_z && (opt.plane == "xz" || opt.plane == "yz") && ref.origin[2] >= -1e-12 && ref.corner[2] > 0.0 && ref.shape_xyz[2] > 1) {
        const int nx = ref.shape_xyz[0], ny = ref.shape_xyz[1], nz = ref.shape_xyz[2];
        const int nz_neg = nz - 1;
        nz_ext = nz + nz_neg;
        std::vector<float> merged(static_cast<std::size_t>(nz_ext) * ny * nx);
        std::vector<double> zneg(nz_neg);
        for (int k = 0; k < nz_neg; ++k) zneg[k] = -z[nz - 1 - k - 1];
        for (int k = 0; k < nz_neg; ++k) {
            const int src_k = nz - 1 - k - 1;
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    float v = data[idx3_zyx(src_k, j, i, ny, nx)];
                    if (is_odd_under_mirror_z(ref.var)) v = -v;
                    merged[idx3_zyx(k, j, i, ny, nx)] = v;
                }
            }
        }
        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    merged[idx3_zyx(k + nz_neg, j, i, ny, nx)] = data[idx3_zyx(k, j, i, ny, nx)];
                }
            }
        }
        zext = zneg;
        zext.insert(zext.end(), z.begin(), z.end());
        ext.swap(merged);
    }

    auto closest_index = [](const std::vector<double> &v, double x0) {
        int best = 0;
        double dbest = std::abs(v[0] - x0);
        for (int i = 1; i < static_cast<int>(v.size()); ++i) {
            double d = std::abs(v[i] - x0);
            if (d < dbest) { dbest = d; best = i; }
        }
        return best;
    };

    TileSlice out;
    const int nx = ref.shape_xyz[0], ny = ref.shape_xyz[1];
    if (opt.plane == "xy") {
        int k = closest_index(zext, opt.coord0);
        out.nx = nx; out.ny = ny;
        out.x_km.resize(nx); out.y_km.resize(ny); out.vals.resize(static_cast<std::size_t>(nx) * ny);
        for (int i = 0; i < nx; ++i) out.x_km[i] = x[i] * RSCALE_KM;
        for (int j = 0; j < ny; ++j) out.y_km[j] = y[j] * RSCALE_KM;
        for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) out.vals[static_cast<std::size_t>(j) * nx + i] = ext[idx3_zyx(k, j, i, ny, nx)];
    } else if (opt.plane == "xz") {
        int j = closest_index(y, opt.coord0);
        out.nx = nx; out.ny = nz_ext;
        out.x_km.resize(nx); out.y_km.resize(nz_ext); out.vals.resize(static_cast<std::size_t>(nx) * nz_ext);
        for (int i = 0; i < nx; ++i) out.x_km[i] = x[i] * RSCALE_KM;
        for (int k = 0; k < nz_ext; ++k) out.y_km[k] = zext[k] * RSCALE_KM;
        for (int k = 0; k < nz_ext; ++k) for (int i = 0; i < nx; ++i) out.vals[static_cast<std::size_t>(k) * nx + i] = ext[idx3_zyx(k, j, i, ny, nx)];
    } else {
        int i = closest_index(x, opt.coord0);
        out.nx = ny; out.ny = nz_ext;
        out.x_km.resize(ny); out.y_km.resize(nz_ext); out.vals.resize(static_cast<std::size_t>(ny) * nz_ext);
        for (int j = 0; j < ny; ++j) out.x_km[j] = y[j] * RSCALE_KM;
        for (int k = 0; k < nz_ext; ++k) out.y_km[k] = zext[k] * RSCALE_KM;
        for (int k = 0; k < nz_ext; ++k) for (int j = 0; j < ny; ++j) out.vals[static_cast<std::size_t>(k) * ny + j] = ext[idx3_zyx(k, j, i, ny, nx)];
    }
    return out;
}

std::vector<double> centers_to_edges(const std::vector<double> &v) {
    if (v.empty()) return {};
    if (v.size() == 1) return {v[0] - 0.5, v[0] + 0.5};
    std::vector<double> e(v.size() + 1);
    e[0] = v[0] - 0.5 * (v[1] - v[0]);
    for (std::size_t i = 1; i < v.size(); ++i) e[i] = 0.5 * (v[i-1] + v[i]);
    e.back() = v.back() + 0.5 * (v.back() - v[v.size()-2]);
    return e;
}

void paste_regular_tile(Field2D &canvas, const Mesh2D &mesh, const TileSlice &tile) {
    auto ex = centers_to_edges(mesh.x_km);
    auto ey = centers_to_edges(mesh.y_km);
    auto tx = centers_to_edges(tile.x_km);
    auto ty = centers_to_edges(tile.y_km);
    for (int j = 0; j < tile.ny; ++j) {
        double ylo = std::min(ty[j], ty[j+1]);
        double yhi = std::max(ty[j], ty[j+1]);
        auto gy0 = static_cast<int>(std::max<std::ptrdiff_t>(0, std::upper_bound(ey.begin(), ey.end(), ylo) - ey.begin() - 1));
        auto gy1 = static_cast<int>(std::min<std::ptrdiff_t>(canvas.ny, std::lower_bound(ey.begin(), ey.end(), yhi) - ey.begin()));
        if (gy1 <= gy0) continue;
        for (int i = 0; i < tile.nx; ++i) {
            float v = tile.vals[static_cast<std::size_t>(j) * tile.nx + i];
            if (!std::isfinite(v)) continue;
            double xlo = std::min(tx[i], tx[i+1]);
            double xhi = std::max(tx[i], tx[i+1]);
            auto gx0 = static_cast<int>(std::max<std::ptrdiff_t>(0, std::upper_bound(ex.begin(), ex.end(), xlo) - ex.begin() - 1));
            auto gx1 = static_cast<int>(std::min<std::ptrdiff_t>(canvas.nx, std::lower_bound(ex.begin(), ex.end(), xhi) - ex.begin()));
            if (gx1 <= gx0) continue;
            for (int jj = gy0; jj < gy1; ++jj) {
                for (int ii = gx0; ii < gx1; ++ii) canvas(jj, ii) = v;
            }
        }
    }
}

Field2D join_var_to_mesh(const VarIndex &index, const std::string &var, int iteration, const Mesh2D &mesh, const Options &opt) {
    Field2D out(mesh.ny, mesh.nx);
    auto refs = collect_refs(index, var, iteration, opt);
    // Coarse-to-fine overwrite: each component is unique after deduplication.
    std::sort(refs.begin(), refs.end(), [](const DatasetRef &a, const DatasetRef &b) {
        if (a.rl != b.rl) return a.rl < b.rl;
        if (a.c != b.c) return a.c < b.c;
        return a.file.string() < b.file.string();
    });
    for (const auto &ref : refs) {
        auto data = read_full_dataset_float(ref);
        auto tile = extract_tile_slice(ref, data, opt);
        paste_regular_tile(out, mesh, tile);
    }
    return out;
}

inline double spatial_dot(double ax, double ay, double az,
                          double bx, double by, double bz,
                          double gxx, double gxy, double gxz,
                          double gyy, double gyz, double gzz) {
    return gxx * ax * bx + gyy * ay * by + gzz * az * bz
         + gxy * (ax * by + ay * bx)
         + gxz * (ax * bz + az * bx)
         + gyz * (ay * bz + az * by);
}

void normalize_basis(double rx, double ry, double rz,
                     double gxx, double gxy, double gxz,
                     double gyy, double gyz, double gzz,
                     double &ux, double &uy, double &uz) {
    double n2 = spatial_dot(rx, ry, rz, rx, ry, rz, gxx, gxy, gxz, gyy, gyz, gzz);
    if (!(n2 > 1e-28) || !std::isfinite(n2)) {
        ux = uy = uz = std::numeric_limits<double>::quiet_NaN();
        return;
    }
    double invn = 1.0 / std::sqrt(n2);
    ux = rx * invn; uy = ry * invn; uz = rz * invn;
}

SliceFields compute_derived(const Mesh2D &mesh,
                            const std::unordered_map<std::string, Field2D> &f,
                            const Options &opt) {
    SliceFields s;
    s.bx = Field2D(mesh.ny, mesh.nx);
    s.by = Field2D(mesh.ny, mesh.nx);
    s.bz = Field2D(mesh.ny, mesh.nx);
    s.b0 = Field2D(mesh.ny, mesh.nx);
    s.b2 = Field2D(mesh.ny, mesh.nx);
    s.W  = Field2D(mesh.ny, mesh.nx);
    s.bpol = Field2D(mesh.ny, mesh.nx);
    s.btor = Field2D(mesh.ny, mesh.nx);

    for (std::size_t j = 0; j < mesh.ny; ++j) {
        for (std::size_t i = 0; i < mesh.nx; ++i) {
            auto get = [&](const std::string &name) -> double { return f.at(name)(j, i); };
            const double Bx = get("Bx"), By = get("By"), Bz = get("Bz");
            const double vx = get("vx"), vy = get("vy"), vz = get("vz");
            const double alp = get("alp");
            const double betax = get("betax"), betay = get("betay"), betaz = get("betaz");
            const double gxx = get("gxx"), gxy = get("gxy"), gxz = get("gxz");
            const double gyy = get("gyy"), gyz = get("gyz"), gzz = get("gzz");
            double W = std::numeric_limits<double>::quiet_NaN();
            if (f.count("W") && std::isfinite(f.at("W")(j, i))) {
                W = f.at("W")(j, i);
            } else {
                double vsq = spatial_dot(vx, vy, vz, vx, vy, vz, gxx, gxy, gxz, gyy, gyz, gzz);
                if (std::isfinite(vsq)) {
                    vsq = std::clamp(vsq, 0.0, 1.0 - 1e-12);
                    W = 1.0 / std::sqrt(1.0 - vsq);
                }
            }
            if (!(std::isfinite(Bx) && std::isfinite(By) && std::isfinite(Bz) && std::isfinite(vx) && std::isfinite(vy) && std::isfinite(vz) &&
                  std::isfinite(alp) && std::isfinite(W) && std::abs(alp) > 1e-14 && std::isfinite(gxx) && std::isfinite(gxy) && std::isfinite(gxz) &&
                  std::isfinite(gyy) && std::isfinite(gyz) && std::isfinite(gzz))) {
                continue;
            }

            double vdotB = spatial_dot(Bx, By, Bz, vx, vy, vz, gxx, gxy, gxz, gyy, gyz, gzz);
            double b0 = W * vdotB / alp;
            double ux = W * (vx - betax / alp);
            double uy = W * (vy - betay / alp);
            double uz = W * (vz - betaz / alp);
            double bx = (Bx + alp * b0 * ux) / W;
            double by = (By + alp * b0 * uy) / W;
            double bz = (Bz + alp * b0 * uz) / W;
            double b2 = spatial_dot(bx, by, bz, bx, by, bz, gxx, gxy, gxz, gyy, gyz, gzz) - (alp * alp) * (b0 * b0);

            double Xc = 0.0, Yc = 0.0, Zc = 0.0;
            double xcode = mesh.x_km[i] / RSCALE_KM;
            double ycode = mesh.y_km[j] / RSCALE_KM;
            if (opt.plane == "xy") {
                Xc = xcode; Yc = ycode; Zc = opt.coord0;
            } else if (opt.plane == "xz") {
                Xc = xcode; Yc = opt.coord0; Zc = ycode;
            } else {
                Xc = opt.coord0; Yc = xcode; Zc = ycode;
            }

            double r_cyl = std::sqrt(Xc * Xc + Yc * Yc);
            if (!(r_cyl > 1e-14)) {
                s.bx(j,i) = static_cast<float>(bx);
                s.by(j,i) = static_cast<float>(by);
                s.bz(j,i) = static_cast<float>(bz);
                s.b0(j,i) = static_cast<float>(b0);
                s.b2(j,i) = static_cast<float>(b2);
                s.W(j,i)  = static_cast<float>(W);
                continue;
            }

            double eRx, eRy, eRz;
            normalize_basis(Xc, Yc, 0.0, gxx, gxy, gxz, gyy, gyz, gzz, eRx, eRy, eRz);
            double ephix, ephiy, ephiz;
            normalize_basis(-Yc, Xc, 0.0, gxx, gxy, gxz, gyy, gyz, gzz, ephix, ephiy, ephiz);
            double z_dot_eR   = spatial_dot(0.0, 0.0, 1.0, eRx, eRy, eRz, gxx, gxy, gxz, gyy, gyz, gzz);
            double z_dot_ephi = spatial_dot(0.0, 0.0, 1.0, ephix, ephiy, ephiz, gxx, gxy, gxz, gyy, gyz, gzz);
            double ezx, ezy, ezz;
            normalize_basis(-z_dot_eR * eRx - z_dot_ephi * ephix,
                            -z_dot_eR * eRy - z_dot_ephi * ephiy,
                             1.0 - z_dot_eR * eRz - z_dot_ephi * ephiz,
                            gxx, gxy, gxz, gyy, gyz, gzz, ezx, ezy, ezz);
            if (!(std::isfinite(eRx) && std::isfinite(ephix) && std::isfinite(ezx))) continue;

            double bR = spatial_dot(bx, by, bz, eRx, eRy, eRz, gxx, gxy, gxz, gyy, gyz, gzz);
            double bphi = spatial_dot(bx, by, bz, ephix, ephiy, ephiz, gxx, gxy, gxz, gyy, gyz, gzz);
            double bzproj = spatial_dot(bx, by, bz, ezx, ezy, ezz, gxx, gxy, gxz, gyy, gyz, gzz);
            double bpol = std::sqrt(std::max(0.0, bR * bR + bzproj * bzproj));

            s.bx(j,i) = static_cast<float>(bx);
            s.by(j,i) = static_cast<float>(by);
            s.bz(j,i) = static_cast<float>(bz);
            s.b0(j,i) = static_cast<float>(b0);
            s.b2(j,i) = static_cast<float>(b2);
            s.W(j,i)  = static_cast<float>(W);
            s.bpol(j,i) = static_cast<float>(bpol);
            s.btor(j,i) = static_cast<float>(opt.write_abs_btor ? std::abs(bphi) : bphi);
        }
    }
    return s;
}

void write_string_attr(hid_t obj, const char *name, const std::string &value) {
    hid_t type = H5Tcopy(H5T_C_S1);
    H5Tset_size(type, value.size());
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t attr = H5Acreate2(obj, name, type, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, type, value.data());
    H5Aclose(attr); H5Sclose(space); H5Tclose(type);
}

void write_double_attr(hid_t obj, const char *name, double value) {
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t attr = H5Acreate2(obj, name, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_DOUBLE, &value);
    H5Aclose(attr); H5Sclose(space);
}

void write_int_attr(hid_t obj, const char *name, int value) {
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t attr = H5Acreate2(obj, name, H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_INT, &value);
    H5Aclose(attr); H5Sclose(space);
}

void write_vec1d(hid_t file, const std::string &name, const std::vector<double> &v) {
    hsize_t dims[1] = {v.size()};
    hid_t space = H5Screate_simple(1, dims, nullptr);
    hid_t ds = H5Dcreate2(file, name.c_str(), H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());
    H5Dclose(ds); H5Sclose(space);
}

void write_field2d(hid_t file, const std::string &name, const Field2D &f) {
    hsize_t dims[2] = {f.ny, f.nx};
    hid_t space = H5Screate_simple(2, dims, nullptr);
    hid_t ds = H5Dcreate2(file, name.c_str(), H5T_NATIVE_FLOAT, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(ds, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, f.a.data());
    H5Dclose(ds); H5Sclose(space);
}

void write_output(const fs::path &path, int iteration, double time_ms, const Mesh2D &mesh, const SliceFields &s, const Options &opt) {
    fs::create_directories(path.parent_path());
    hid_t file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) die("Could not create output file " + path.string());

    write_vec1d(file, "axis0_km", mesh.x_km);
    write_vec1d(file, "axis1_km", mesh.y_km);
    write_field2d(file, "bpol", s.bpol);
    write_field2d(file, "btor", s.btor);
    write_field2d(file, "bx_comov", s.bx);
    write_field2d(file, "by_comov", s.by);
    write_field2d(file, "bz_comov", s.bz);
    write_field2d(file, "b0_comov", s.b0);
    write_field2d(file, "b2_comov", s.b2);
    write_field2d(file, "W", s.W);

    write_int_attr(file, "iteration", iteration);
    write_double_attr(file, "time_ms", time_ms);
    write_double_attr(file, "coord0_code", opt.coord0);
    write_string_attr(file, "plane", opt.plane);
    write_string_attr(file, "axis0_name", opt.plane == "yz" ? "y" : "x");
    write_string_attr(file, "axis1_name", opt.plane == "xy" ? "y" : "z");
    write_int_attr(file, "mirror_z", opt.mirror_z ? 1 : 0);
    write_string_attr(file, "btor_convention", opt.write_abs_btor ? "abs(bphi)" : "signed bphi");
    write_string_attr(file, "formula_bpol", "sqrt(b_R^2 + b_zproj^2) using orthonormal cylindrical basis from spatial metric");
    write_string_attr(file, "formula_btor", opt.write_abs_btor ? "abs(b_phi) using orthonormal cylindrical basis from spatial metric" : "b_phi using orthonormal cylindrical basis from spatial metric");

    H5Fclose(file);
}

int main(int argc, char **argv) {
    try {
        Options opt = parse_args(argc, argv);
        VarIndex index = build_index(opt);
        auto iterations = choose_iterations(opt, index);

        std::cout << "Found " << iterations.size() << " eligible iteration(s).\n";
        auto times = iteration_to_time_ms(index);

        for (int it : iterations) {
            auto available = available_rls_for_iteration(index, it, false);
            if (available.empty()) {
                std::cerr << "[skip] iteration " << it << ": no common rl across required variables\n";
                continue;
            }
            if (opt.requested_rl && std::find(available.begin(), available.end(), *opt.requested_rl) == available.end()) {
                std::cerr << "[skip] iteration " << it << ": requested rl=" << *opt.requested_rl << " not available\n";
                continue;
            }

            Mesh2D mesh = build_mesh_for_iteration(index, it, opt);
            std::unordered_map<std::string, Field2D> fields;
            for (const auto &name : {std::string("rho"), std::string("Bx"), std::string("By"), std::string("Bz"),
                                     std::string("vx"), std::string("vy"), std::string("vz"),
                                     std::string("alp"), std::string("betax"), std::string("betay"), std::string("betaz"),
                                     std::string("gxx"), std::string("gxy"), std::string("gxz"), std::string("gyy"), std::string("gyz"), std::string("gzz")}) {
                fields.emplace(name, join_var_to_mesh(index, name, it, mesh, opt));
            }
            if (index.count("W")) fields.emplace("W", join_var_to_mesh(index, "W", it, mesh, opt));

            SliceFields derived = compute_derived(mesh, fields, opt);
            double tms = times[it];
            double tplot = opt.tmerg_code ? (tms - (*opt.tmerg_code) * CU_TO_MS) : tms;

            std::ostringstream fname;
            fname << "bfield_slice_" << opt.plane << "_it" << std::setw(8) << std::setfill('0') << it << ".h5";
            fs::path out = opt.outdir / opt.plane / fname.str();
            write_output(out, it, tms, mesh, derived, opt);

            std::cout << "[ok] it=" << it
                      << " rls=" << (opt.requested_rl ? std::to_string(*opt.requested_rl) : std::string("all"))
                      << " mesh=" << mesh.nx << "x" << mesh.ny
                      << " t_ms=" << tms
                      << " t_plot_ms=" << tplot
                      << " -> " << out << "\n";
        }
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
