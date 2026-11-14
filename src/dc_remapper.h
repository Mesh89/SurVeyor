#ifndef DC_REMAPPER_H_
#define DC_REMAPPER_H_

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "types.h"
#include "sam_utils.h"
#include "clustering_utils.h"
#include <memory>
#include <set>

struct insertion_cluster_t {

    hts_pos_t start, end;
    std::shared_ptr<cluster_t> cluster = NULL;
    std::vector<bam1_t*> semi_mapped_reads;
    std::vector<bam1_t*> semi_mapped_mates;
    std::shared_ptr<consensus_t> clip_consensus = NULL;

    insertion_cluster_t() : start(0), end(0) {
        cluster = std::make_shared<cluster_t>();
    }

    insertion_cluster_t(std::shared_ptr<cluster_t> cluster) : cluster(cluster), start(cluster->la_start), end(cluster->la_end) {
        std::vector<std::shared_ptr<bam1_t>> reads;
        reads.swap(cluster->reads);
        for (std::shared_ptr<bam1_t> read : reads) {
            add_stable_read(read);
        }
    }

    void add_stable_read(std::shared_ptr<bam1_t> read) {
        cluster->reads.push_back(read);
        start = std::min(start, read->core.pos);
        end = std::max(end, bam_endpos(read.get()));
        cluster->la_max_mapq = std::max(cluster->la_max_mapq, read->core.qual);
        cluster->ra_max_mapq = std::max(cluster->ra_max_mapq, read->core.qual);
    }
    void add_semi_mapped_reads(bam1_t* read) {
        if (clip_consensus) {
            for (auto clipped_read : clip_consensus->clipped_reads) {
                if (get_read_key(read) == get_read_key(clipped_read.get())) return;
            }
        }
    	semi_mapped_reads.push_back(read);
        start = std::min(start, read->core.mpos);
        end = std::max(end, get_mate_endpos(read));
    }
    void add_clip_cluster(std::shared_ptr<consensus_t> clip_consensus) {
        if (!clip_consensus) return;
        if (!this->clip_consensus || this->clip_consensus->reads() < clip_consensus->reads()) {
            this->clip_consensus = clip_consensus;
            if (clip_consensus->left_clipped) start = clip_consensus->anchor_start();
            else start = std::min(start, clip_consensus->anchor_start());
            if (!clip_consensus->left_clipped) end = clip_consensus->anchor_end();
            else end = std::max(end, clip_consensus->anchor_end());
        }
    }

    void deallocate_reads() {
        cluster->reads.clear();
        for (bam1_t* read : semi_mapped_reads) {
        	bam_destroy1(read);
        }
        semi_mapped_reads.clear();
        for (bam1_t* read : semi_mapped_mates) {
        	bam_destroy1(read);
        }
        semi_mapped_mates.clear();
    }

    int size() {
        return cluster->reads.size() + semi_mapped_reads.size() + (clip_consensus ? 1 : 0);
    }

    bool empty() {
    	return this->size() == 0;
    }
};

size_t count_assembly_reads(std::shared_ptr<insertion_cluster_t> r_cluster, std::shared_ptr<insertion_cluster_t> l_cluster) {
    std::set<std::pair<std::string, int>> reads;
    for (auto cluster : {r_cluster, l_cluster}) {
        if (cluster->clip_consensus) {
            for (auto read : cluster->clip_consensus->clipped_reads) reads.insert(get_read_key(read.get()));
        }
        for (auto read : cluster->cluster->reads) {
            reads.insert(get_read_key(read.get()));
            reads.insert(get_read_key(read.get(), true));
        }
        for (bam1_t* read : cluster->semi_mapped_reads) reads.insert(get_read_key(read));
    }
    return reads.size();
}

struct assembly_read_filter_t {
    std::set<std::pair<std::string, int>> clipped_reads, bam_reads, used_reads;

    assembly_read_filter_t(std::shared_ptr<insertion_cluster_t> r_cluster, std::shared_ptr<insertion_cluster_t> l_cluster, bool include_semi_mapped_mates = false) {
        for (auto cluster : {r_cluster, l_cluster}) {
            if (cluster->clip_consensus) {
                for (auto read : cluster->clip_consensus->clipped_reads) clipped_reads.insert(get_read_key(read.get()));
            }
            for (auto read : cluster->cluster->reads) bam_reads.insert(get_read_key(read.get()));
            for (bam1_t* read : cluster->semi_mapped_reads) bam_reads.insert(get_read_key(read));
            if (include_semi_mapped_mates) {
                for (bam1_t* read : cluster->semi_mapped_mates) bam_reads.insert(get_read_key(read));
            }
        }
        bam_reads.insert(clipped_reads.begin(), clipped_reads.end());
    }

    bool add_read(bam1_t* read, bool clipped = false) {
        auto key = get_read_key(read);
        if (!clipped && clipped_reads.count(key)) return false;
        return used_reads.insert(key).second;
    }

    bool add_mate(bam1_t* read) {
        auto key = get_read_key(read, true);
        return !bam_reads.count(key) && used_reads.insert(key).second;
    }
};

std::string get_mate_seq(bam1_t* read, std::unordered_map<std::string, std::string>& mateseqs) {
    std::string qname = get_mate_lookup_qname(read);
    if (!mateseqs.count(qname)) {
    	std::cerr << "Warning: mate not found for " << qname << std::endl;
    	return "";
    }
    return mateseqs[qname];
}
std::string get_mate_qual(bam1_t* read, std::unordered_map<std::string, std::string>& matequals) {
    std::string qname = get_mate_lookup_qname(read);
    if (!matequals.count(qname)) {
    	std::cerr << "Warning: mate not found for " << qname << std::endl;
    	return "";
    }
    return matequals[qname];
}


#endif /* DC_REMAPPER_H_ */
