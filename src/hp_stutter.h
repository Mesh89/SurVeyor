#ifndef HP_STUTTER_H
#define HP_STUTTER_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

const int MIN_READS_FOR_HP_STUTTER_RATE = 100;

// Empirical HP stutter rates for one library; unobserved bins remain zero.
struct hp_stutter_profile_t {
    double rates[4][31][7] = {}; // A/C/G/T, HP length, delta + 3.
    long long reads[4][31] = {};

    // Construct a library profile from TSV file fname; unobserved bins remain zero.
    // Each delta row repeats the read total for its base and HP length.
    explicit hp_stutter_profile_t(const std::string& fname) {
        std::ifstream fin(fname);
        if (!fin) throw std::runtime_error("Unable to open " + fname + ".");
        std::string line;
        if (!std::getline(fin, line) || line != "BASE\tL\tDELTA\tREADS\tRATE") throw std::runtime_error("Invalid HP stutter profile header in " + fname + ".");
        const std::string bases = "ACGT";
        while (std::getline(fin, line)) {
            if (line.empty()) continue;
            std::istringstream row(line);
            std::string base, extra;
            int hp_len, delta;
            long long read_count;
            double stutter_rate;
            if (!(row >> base >> hp_len >> delta >> read_count >> stutter_rate) || row >> extra || base.size() != 1 || bases.find(base) == std::string::npos || hp_len < 5 || hp_len > 30 || delta < -3 || delta > 3 || read_count <= 0 || !std::isfinite(stutter_rate) || stutter_rate < 0 || stutter_rate > 1) throw std::runtime_error("Invalid HP stutter profile row in " + fname + ".");
            int base_idx = bases.find(base);
            rates[base_idx][hp_len][delta + 3] = stutter_rate;
            reads[base_idx][hp_len] = read_count;
        }
        if (fin.bad()) throw std::runtime_error("Unable to read HP stutter profile " + fname + ".");
    }

    // Return the empirical stutter rate for HP base b (A/C/G/T) and true length L (5..30).
    // d is observed minus true length (-3..3); +/-3 bins include larger changes.
    // Prefer the base alone, then pool A/T or C/G, then try shorter lengths; assume no stutter if unsupported.
    double rate(char b, int L, int d) const {
        int base = std::string("ACGT").find(b);
        int partner = 3 - base;
        for (int LL = L; LL >= 5; LL--) {
            if (reads[base][LL] >= MIN_READS_FOR_HP_STUTTER_RATE) return rates[base][LL][d + 3];
            long long n = reads[base][LL] + reads[partner][LL];
            if (n >= MIN_READS_FOR_HP_STUTTER_RATE) return (rates[base][LL][d + 3] * reads[base][LL] + rates[partner][LL][d + 3] * reads[partner][LL]) / n;
        }
        return d == 0 ? 1 : 0; // No data: no stutter.
    }

    // Return P(observed length o | true length t, HP base b), floored at 1e-5; b is A/C/G/T.
    // Profile lookup lengths are clamped to 5..30; +/-3 bins are spread over geometric tails.
    double K(int o, int t, char b) const {
        // Clamp the profile bin while measuring the length change from the actual true length.
        int tt = std::min(std::max(t, 5), 30);
        int d = o - t;
        double value;
        if (std::abs(d) <= 2) {
            // Use measured rates directly for unchanged lengths and changes of one or two bases.
            value = rate(b, tt, d);
        } else {
            // Split the +/-3 bin across distances 3, 4, ... with weights 1/2, 1/4, ... .
            value = rate(b, tt, d > 0 ? 3 : -3) * std::pow(0.5, std::abs(d) - 2);
        }
        // Keep unseen or very rare outcomes usable in likelihood calculations.
        return std::max(value, 1e-5);
    }
};

// Return the length in nonempty allele-length vector G closest to observed HP length o.
// Resolve equal distances toward the shorter length, independent of G's ordering.
inline int nearest(int o, const std::vector<int>& G) {
    int best = G.front();
    int best_dist = std::abs(o - best);
    for (int g : G) {
        int dist = std::abs(o - g);
        if (dist < best_dist || (dist == best_dist && g < best)) {
            best = g;
            best_dist = dist;
        }
    }
    return best;
}

// Return sorted unique candidate lengths from allele lengths G and histogram h (length -> read count).
// Add observed lengths with at least two reads among the eight largest bins; count ties favor shorter lengths.
// Note: allele lengths not present in G but present in T are NOT reported as extra candidates. They are only used in subsequent calculations
inline std::vector<int> candidate_true_lengths(const std::map<int,int>& h, const std::vector<int>& G) {
    std::vector<std::pair<int,int>> observed(h.begin(), h.end());
    std::sort(observed.begin(), observed.end(), [](const std::pair<int,int>& a, const std::pair<int,int>& b) { return a.second != b.second ? a.second > b.second : a.first < b.first; });
    std::vector<int> T = G;
    for (size_t i = 0; i < observed.size() && i < 8; i++) {
        if (observed[i].second >= 2) T.push_back(observed[i].first);
    }
    std::sort(T.begin(), T.end());
    T.erase(std::unique(T.begin(), T.end()), T.end());
    return T;
}

// Return unconstrained least-squares weights for the passive columns of row-major M and observations hv.
// Inactive weights are zero; the selected columns are linearly independent.
inline std::vector<double> hp_stutter_least_squares(const std::vector<std::vector<double>>& M, const std::vector<double>& hv, const std::vector<bool>& passive) {
    int rows = M.size(), cols = passive.size();
    std::vector<int> selected;
    for (int j = 0; j < cols; j++) if (passive[j]) selected.push_back(j);
    int count = selected.size();
    std::vector<std::vector<double>> qr(rows, std::vector<double>(count));
    std::vector<double> rhs = hv, weights(cols, 0);
    for (int i = 0; i < rows; i++) for (int j = 0; j < count; j++) qr[i][j] = M[i][selected[j]];

    // Householder reflections solve each passive fit without forming normal equations.
    for (int j = 0; j < count; j++) {
        double norm = 0;
        for (int i = j; i < rows; i++) norm = std::hypot(norm, qr[i][j]);
        double diagonal = -std::copysign(norm, qr[j][j]);
        qr[j][j] -= diagonal;
        double norm_sq = 0;
        for (int i = j; i < rows; i++) norm_sq += qr[i][j] * qr[i][j];
        for (int k = j + 1; k < count; k++) {
            double projection = 0;
            for (int i = j; i < rows; i++) projection += qr[i][j] * qr[i][k];
            projection *= 2 / norm_sq;
            for (int i = j; i < rows; i++) qr[i][k] -= projection * qr[i][j];
        }
        double projection = 0;
        for (int i = j; i < rows; i++) projection += qr[i][j] * rhs[i];
        projection *= 2 / norm_sq;
        for (int i = j; i < rows; i++) rhs[i] -= projection * qr[i][j];
        qr[j][j] = diagonal;
    }
    for (int j = count - 1; j >= 0; j--) {
        double value = rhs[j];
        for (int k = j + 1; k < count; k++) value -= qr[j][k] * weights[selected[k]];
        weights[selected[j]] = value / qr[j][j];
    }
    return weights;
}

// Return nonnegative weights minimising ||M w - hv||_2 for full-column-rank row-major M and observations hv.
inline std::vector<double> hp_stutter_nnls(const std::vector<std::vector<double>>& M, const std::vector<double>& hv) {
    int rows = M.size(), cols = M[0].size();
    std::vector<double> weights(cols, 0), gradient(cols, 0);
    std::vector<bool> passive(cols, false);
    auto update_gradient = [&]() {
        std::fill(gradient.begin(), gradient.end(), 0);
        for (int i = 0; i < rows; i++) {
            double residual = hv[i];
            for (int j = 0; j < cols; j++) residual -= M[i][j] * weights[j];
            for (int j = 0; j < cols; j++) gradient[j] += M[i][j] * residual;
        }
    };
    double matrix_norm = 0, observations_norm_sq = 0;
    for (int j = 0; j < cols; j++) {
        double column_sum = 0;
        for (int i = 0; i < rows; i++) column_sum += std::abs(M[i][j]);
        matrix_norm = std::max(matrix_norm, column_sum);
    }
    for (double value : hv) observations_norm_sq += value * value;
    double tolerance = 10 * std::numeric_limits<double>::epsilon() * rows * matrix_norm * std::sqrt(observations_norm_sq);
    update_gradient();

    // Lawson-Hanson: admit the column that most reduces the residual, then refit the passive set.
    while (true) {
        int entering = -1;
        for (int j = 0; j < cols; j++) {
            if (!passive[j] && gradient[j] > tolerance && (entering == -1 || gradient[j] > gradient[entering])) entering = j;
        }
        if (entering == -1) break;
        passive[entering] = true;
        while (true) {
            std::vector<double> trial = hp_stutter_least_squares(M, hv, passive);
            double step = 1;
            int blocking = -1;
            for (int j = 0; j < cols; j++) {
                if (!passive[j] || trial[j] > 0) continue;
                double fraction = weights[j] == 0 ? 0 : weights[j] / (weights[j] - trial[j]);
                if (blocking == -1 || fraction < step) {
                    step = fraction;
                    blocking = j;
                }
            }
            if (blocking == -1) {
                weights = std::move(trial);
                break;
            }
            // Move only as far as the first zero weight, then remove that column from the passive set.
            for (int j = 0; j < cols; j++) weights[j] += step * (trial[j] - weights[j]);
            weights[blocking] = 0;
            passive[blocking] = false;
            for (int j = 0; j < cols; j++) {
                if (passive[j] && weights[j] <= 0) {
                    weights[j] = 0;
                    passive[j] = false;
                }
            }
        }
        update_gradient();
    }
    return weights;
}

// Return the genuine-read fraction for allele length x in G, using histogram h and the profile for HP base b.
// An allele bin contains reads nearest to x over G; return NaN if that bin has no reads.
inline double hp_stutter_alpha(const hp_stutter_profile_t& profile, char b, const std::map<int,int>& h, const std::vector<int>& G, int x) {
    double den = 0;
    for (const auto& observation : h) {
        if (observation.second > 0 && nearest(observation.first, G) == x) den += observation.second;
    }
    if (den == 0) return std::numeric_limits<double>::quiet_NaN();
    std::vector<int> T = candidate_true_lengths(h, G);
    int lo = std::max(0, std::min(T.front(), h.begin()->first) - 3);
    int hi = std::max(T.back(), h.rbegin()->first) + 3;
    int rows = hi - lo + 1, cols = T.size();
    std::vector<std::vector<double>> M(rows, std::vector<double>(cols));
    std::vector<double> hv(rows, 0);
    for (const auto& observation : h) hv[observation.first - lo] = observation.second;
    // Normalise over the full window, including rows with no observed reads.
    for (int j = 0; j < cols; j++) {
        double total = 0;
        for (int i = 0; i < rows; i++) total += M[i][j] = profile.K(lo + i, T[j], b);
        for (int i = 0; i < rows; i++) M[i][j] /= total;
    }
    std::vector<double> weights = hp_stutter_nnls(M, hv);
    int column = std::lower_bound(T.begin(), T.end(), x) - T.begin();
    double num = 0;
    for (int i = 0; i < rows; i++) {
        if (hv[i] <= 0 || nearest(lo + i, G) != x) continue;
        double predicted = 0;
        for (int j = 0; j < cols; j++) predicted += M[i][j] * weights[j];
        if (predicted == 0) continue;
        num += hv[i] * M[i][column] * weights[column] / predicted;
    }
    return std::min(1.0, num / den);
}

// Return the best diploid log-likelihood among genotypes containing x minus the best gt excluding x, over candidate lengths T.
// Use histogram h and raw profile probabilities for HP base b; clamp to +/-1e4, or return 1e4 if no pair excludes x.
inline double hp_stutter_run_score(const hp_stutter_profile_t& profile, char b, const std::map<int,int>& h, const std::vector<int>& T, int x) {
    double best_with = -std::numeric_limits<double>::infinity();
    double best_without = -std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < T.size(); i++) {
        for (size_t j = i; j < T.size(); j++) {
            double LL = 0;
            for (const auto& observation : h) {
                if (observation.second > 0) LL += observation.second * std::log(0.5 * profile.K(observation.first, T[i], b) + 0.5 * profile.K(observation.first, T[j], b));
            }
            if (T[i] == x || T[j] == x) best_with = std::max(best_with, LL);
            else best_without = std::max(best_without, LL);
        }
    }
    if (best_without == -std::numeric_limits<double>::infinity()) return 1e4;
    return std::max(-1e4, std::min(1e4, best_with - best_without));
}

#endif
