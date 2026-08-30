#ifndef HP_MISMATCH_RATE_THRESHOLDS_H
#define HP_MISMATCH_RATE_THRESHOLDS_H

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

static const char HP_MISMATCH_RATE_THRESHOLDS_FILENAME[] = "hp_3p_mismatch_rate_thresholds.tsv";

struct hp_mismatch_rate_thresholds_t {
    std::array<std::vector<double>, 4> thresholds_by_base;

    explicit hp_mismatch_rate_thresholds_t(const std::string& fname) { load(fname); }

    void load(const std::string& fname) {
        std::ifstream fin(fname);
        if (!fin) throw std::runtime_error("Unable to open " + fname + ".");

        std::string hp_len_header, a_header, c_header, g_header, t_header;
        if (!(fin >> hp_len_header >> a_header >> c_header >> g_header >> t_header) || hp_len_header != "HP_LEN" || a_header != "A" || c_header != "C" || g_header != "G" || t_header != "T") {
            throw std::runtime_error("Invalid header in " + fname + ".");
        }

        for (std::vector<double>& thresholds : thresholds_by_base) thresholds.clear();
        int hp_len;
        double thresholds[4];
        int expected_hp_len = 0;
        while (fin >> hp_len) {
            if (!(fin >> thresholds[0] >> thresholds[1] >> thresholds[2] >> thresholds[3])) throw std::runtime_error("Invalid data in " + fname + ".");
            if (hp_len != expected_hp_len) throw std::runtime_error("Invalid HP length in " + fname + ".");
            for (int hp_base_idx = 0; hp_base_idx < 4; hp_base_idx++) {
                if (!std::isfinite(thresholds[hp_base_idx]) || thresholds[hp_base_idx] < 0 || thresholds[hp_base_idx] > 1) throw std::runtime_error("Invalid threshold in " + fname + ".");
                thresholds_by_base[hp_base_idx].push_back(thresholds[hp_base_idx]);
            }
            expected_hp_len++;
        }
        if (!fin.eof()) throw std::runtime_error("Invalid data in " + fname + ".");
        if (expected_hp_len == 0) throw std::runtime_error("No HP mismatch-rate thresholds in " + fname + ".");
    }

    double get_threshold(int hp_len, char hp_base) const {
        if (hp_len < 0) throw std::invalid_argument("HP length cannot be negative.");
        int hp_base_idx = get_base_idx(hp_base);
        if (hp_base_idx < 0) throw std::invalid_argument("HP base must be A, C, G, or T.");
        const std::vector<double>& thresholds = thresholds_by_base[hp_base_idx];
        return thresholds[std::min<size_t>(hp_len, thresholds.size() - 1)];
    }

private:
    static int get_base_idx(char base) {
        if (base == 'A' || base == 'a') return 0;
        if (base == 'C' || base == 'c') return 1;
        if (base == 'G' || base == 'g') return 2;
        if (base == 'T' || base == 't') return 3;
        return -1;
    }
};

#endif // HP_MISMATCH_RATE_THRESHOLDS_H
