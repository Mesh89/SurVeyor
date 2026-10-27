#ifndef DC_REMAPPER_H_
#define DC_REMAPPER_H_

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "types.h"
#include "sam_utils.h"
#include "clustering_utils.h"
#include <memory>

struct insertion_cluster_t {

    hts_pos_t start, end;
    std::shared_ptr<cluster_t> cluster = NULL;
    std::vector<bam1_t*> semi_mapped_reads;
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
            bool sm_left_clipped = is_left_clipped(read, 1);
            bool sm_right_clipped = is_right_clipped(read, 1);
            if (sm_left_clipped && clip_consensus->left_clipped && read->core.pos == clip_consensus->breakpoint && get_left_clip_size(read) <= clip_consensus->clip_len) {
                return;
            } else if (sm_right_clipped && !clip_consensus->left_clipped && bam_endpos(read) == clip_consensus->breakpoint && get_right_clip_size(read) <= clip_consensus->clip_len) {
                return;
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
    }

    int size() {
        return cluster->reads.size() + semi_mapped_reads.size() + (clip_consensus ? 1 : 0);
    }

    bool empty() {
    	return this->size() == 0;
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
