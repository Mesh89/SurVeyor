#ifndef CLUSTERING_UTILS_H
#define CLUSTERING_UTILS_H

#include "htslib/sam.h"
#include "sam_utils.h"
#include "types.h"
#include <cstdint>
#include <iostream>
#include <memory>
#include <tuple>
#include <list>
#include <queue>
#include <map>

struct cluster_t {
	int id = -1;
	bool la_rev, ra_rev;
	hts_pos_t la_start, la_end, ra_start, ra_end;
	hts_pos_t la_start_unclipped, la_end_unclipped, ra_start_unclipped, ra_end_unclipped;
	int count = 1, la_confident_count = 0, ra_confident_count = 0;
	uint8_t la_max_mapq, ra_max_mapq;
	bool used = false;

	// furthermost sequence for left and right-anchor, respectively, along the direction of the anchor
	// meaning, if l/ra_rev is true, then it is the leftmost sequence, otherwise it is the rightmost sequence
	std::string la_furthermost_seq, ra_furthermost_seq;
	
	std::vector<std::shared_ptr<bam1_t>> reads; // reads that are part of this cluster

	cluster_t() : 
		la_start(INT32_MAX), la_end(0), la_rev(false), la_start_unclipped(INT32_MAX), la_end_unclipped(0),
		ra_start(INT32_MAX), ra_end(0), ra_rev(false), ra_start_unclipped(INT32_MAX), ra_end_unclipped(0),
		la_max_mapq(0), ra_max_mapq(0) {}

	cluster_t(bam1_t* read, int high_confidence_mapq) : 
			la_start(read->core.pos), la_end(bam_endpos(read)), la_rev(bam_is_rev(read)), 
			la_start_unclipped(get_unclipped_start(read)), la_end_unclipped(get_unclipped_end(read)),
			ra_start(read->core.mpos), ra_end(get_mate_endpos(read)), ra_rev(bam_is_mrev(read)),
			ra_start_unclipped(get_mate_unclipped_start(read)), ra_end_unclipped(get_mate_unclipped_end(read)),
			la_max_mapq(read->core.qual), ra_max_mapq(get_mq(read)) {
		
		if (read->core.qual >= high_confidence_mapq) {
			la_confident_count = 1;
		}
		if ((uint8_t) get_mq(read) >= high_confidence_mapq) {
			ra_confident_count = 1;
		}
		la_furthermost_seq = get_sequence(read);
		ra_furthermost_seq = bam_get_qname(read);
	}

	static hts_pos_t distance(const std::shared_ptr<cluster_t>& c1, const std::shared_ptr<cluster_t>& c2) {
		hts_pos_t la_dist = std::max(c1->la_end, c2->la_end) - std::min(c1->la_start, c2->la_start);
		hts_pos_t ra_dist = std::max(c1->ra_end, c2->ra_end) - std::min(c1->ra_start, c2->ra_start);
		return std::max(la_dist, ra_dist);
	}

	static bool can_merge(const std::shared_ptr<cluster_t>& c1, const std::shared_ptr<cluster_t>& c2, int max_distance) {
		return c1->la_rev == c2->la_rev && c1->ra_rev == c2->ra_rev && distance(c1, c2) <= max_distance;
	}

	static std::shared_ptr<cluster_t> merge(const std::shared_ptr<cluster_t>& c1, const std::shared_ptr<cluster_t>& c2) {
		std::shared_ptr<cluster_t> merged = std::make_shared<cluster_t>();
		merged->la_start = std::min(c1->la_start, c2->la_start);
		merged->la_end = std::max(c1->la_end, c2->la_end);
		merged->la_rev = c1->la_rev;
		merged->ra_start = std::min(c1->ra_start, c2->ra_start);
		merged->ra_end = std::max(c1->ra_end, c2->ra_end);
		merged->ra_rev = c1->ra_rev;
		merged->count = c1->count + c2->count;
		merged->la_confident_count += c1->la_confident_count + c2->la_confident_count;
		merged->ra_confident_count += c1->ra_confident_count + c2->ra_confident_count;
		merged->la_max_mapq = std::max(c1->la_max_mapq, c2->la_max_mapq);
		merged->ra_max_mapq = std::max(c1->ra_max_mapq, c2->ra_max_mapq);
		merged->la_start_unclipped = std::min(c1->la_start_unclipped, c2->la_start_unclipped);
		merged->la_end_unclipped = std::max(c1->la_end_unclipped, c2->la_end_unclipped);
		merged->ra_start_unclipped = std::min(c1->ra_start_unclipped, c2->ra_start_unclipped);
		merged->ra_end_unclipped = std::max(c1->ra_end_unclipped, c2->ra_end_unclipped);

		// set leftmost_lseq
		if (c1->la_rev) {
			if (c1->la_start_unclipped <= c2->la_start_unclipped) merged->la_furthermost_seq = c1->la_furthermost_seq;
			else merged->la_furthermost_seq = c2->la_furthermost_seq;
		} else {
			if (c1->la_end_unclipped >= c2->la_end_unclipped) merged->la_furthermost_seq = c1->la_furthermost_seq;
			else merged->la_furthermost_seq = c2->la_furthermost_seq;
		}

		// set rightmost_lseq
		if (c1->ra_rev) {
			if (c1->ra_start_unclipped <= c2->ra_start_unclipped) merged->ra_furthermost_seq = c1->ra_furthermost_seq;
			else merged->ra_furthermost_seq = c2->ra_furthermost_seq;
		} else {
			if (c1->ra_end_unclipped >= c2->ra_end_unclipped) merged->ra_furthermost_seq = c1->ra_furthermost_seq;
			else merged->ra_furthermost_seq = c2->ra_furthermost_seq;
		}

		return merged;
	}
};

bool operator < (const cluster_t& c1, const cluster_t& c2) {
	return std::make_tuple(c1.la_start, c1.la_end, c1.la_rev, c1.ra_start, c1.ra_end, c1.ra_rev) < std::make_tuple(c2.la_start, c2.la_end, c2.la_rev, c2.ra_start, c2.ra_end, c2.ra_rev);
}
bool operator == (const cluster_t& c1, const cluster_t& c2) {
	return c1.la_start == c2.la_start && c1.la_end == c2.la_end && c1.ra_start == c2.ra_start && c1.ra_end == c2.ra_end && c1.la_rev == c2.la_rev && c1.ra_rev == c2.ra_rev;
}
bool operator != (const cluster_t& c1, const cluster_t& c2) {
	return !(c1 == c2);
}

struct cc_distance_t {
    int distance;
    std::shared_ptr<cluster_t> c1, c2;

    cc_distance_t(const std::shared_ptr<cluster_t>& c1, const std::shared_ptr<cluster_t>& c2) : distance(cluster_t::distance(c1, c2)), c1(c1), c2(c2) {}
};
bool operator < (const cc_distance_t& ccd1, const cc_distance_t& ccd2) { // reverse op for priority queue
    return ccd1.distance > ccd2.distance;
}

struct union_find_t {

	int* parents;
	int* sizes;
	int n;

	union_find_t(int n) : n(n) {
		parents = new int[n];
		sizes = new int[n];
		for (int i = 0; i < n; i++) {
			parents[i] = i;
			sizes[i] = 1;
		}
	}

	int find(int i) {
		int root = i;
		while (root != parents[root]) {
			root = parents[root];
		}
		while (i != root) {
			int newp = parents[i];
			parents[i] = root;
			i = newp;
		}
		return root;
	}
	int merge(int x, int y) {
		int i = find(x);
		int j = find(y);
		if (i == j) {
			return i;
		}

		if (sizes[i] < sizes[j]) {
			parents[i] = j;
			sizes[j] += sizes[i];
			return j;
		} else {
			parents[j] = i;
			sizes[i] += sizes[j];
			return i;
		}
	}

	~union_find_t() {
		delete[] parents;
		delete[] sizes;
	}
};

void remove_cluster_from_mm(std::multimap<int, std::shared_ptr<cluster_t>>& mm, std::shared_ptr<cluster_t> c, int pos) {
    auto bounds = mm.equal_range(pos);
    for (auto it = bounds.first; it != bounds.second; it++) {
        if (it->second == c) {
            mm.erase(it);
            break;
        }
    }
}
void remove_cluster_from_mm(std::multimap<int, std::shared_ptr<cluster_t>>& mm, std::shared_ptr<cluster_t> c) {
    remove_cluster_from_mm(mm, c, c->la_start);
    remove_cluster_from_mm(mm, c, c->la_end);
}

void cluster_clusters(std::vector<std::shared_ptr<cluster_t>>& clusters, std::vector<std::shared_ptr<bam1_t>>& reads, int max_distance, int max_cluster_size,
	 bool assign_reads_to_clusters = false) {

	if (clusters.empty()) return;
	
	if (assign_reads_to_clusters) { // sort clusters and reads consistently
		if (clusters.size() != reads.size()) {
			throw std::runtime_error("Error: number of clusters and reads do not match");
		}

		std::vector<std::pair<std::shared_ptr<cluster_t>, std::shared_ptr<bam1_t>>> cluster_reads;
		for (int i = 0; i < clusters.size(); i++) {
			cluster_reads.push_back(std::make_pair(clusters[i], reads[i]));
		}
		std::sort(cluster_reads.begin(), cluster_reads.end(), [](std::pair<std::shared_ptr<cluster_t>, std::shared_ptr<bam1_t>>& cr1, std::pair<std::shared_ptr<cluster_t>, std::shared_ptr<bam1_t>>& cr2) {
			return  std::make_tuple(cr1.first->la_start, cr1.first->la_end, cr1.first->ra_start, cr1.first->ra_end) < 
					std::make_tuple(cr2.first->la_start, cr2.first->la_end, cr2.first->ra_start, cr2.first->ra_end);
		});

		for (int i = 0; i < cluster_reads.size(); i++) {
			clusters[i] = cluster_reads[i].first;
			reads[i] = cluster_reads[i].second;
		}
		for (int i = 0; i < clusters.size(); i++) clusters[i]->id = i;
	} else {
		std::sort(clusters.begin(), clusters.end(), [](std::shared_ptr<cluster_t> c1, std::shared_ptr<cluster_t> c2) {
			return  std::make_tuple(c1->la_start, c1->la_end, c1->ra_start, c1->ra_end) < 
					std::make_tuple(c2->la_start, c2->la_end, c2->ra_start, c2->ra_end);
		});
		for (int i = 0; i < clusters.size(); i++) clusters[i]->id = i;
	}

	int n_reads = clusters.size();
	union_find_t* uf = new union_find_t(n_reads);

	int curr = 0;
    for (int i = 1; i < clusters.size(); i++) {
        if (*clusters[curr] == *clusters[i]) {
            std::shared_ptr<cluster_t> merged = cluster_t::merge(clusters[curr], clusters[i]);
            int parent_id = uf->merge(curr, i);
            merged->id = parent_id;
            clusters[parent_id] = merged;
            clusters[curr+i-parent_id] = NULL;
        } else {
            curr = i;
        }
    }

    std::priority_queue<cc_distance_t> pq;
	for (int i = 0; i < clusters.size(); i++) {
		if (clusters[i] == NULL || clusters[i]->used) continue;

		std::vector<cc_distance_t> ccps;
		for (int j = i+1; j < clusters.size(); j++) {
			if (clusters[j] == NULL || clusters[j]->used) continue;
			if (clusters[j]->la_start-clusters[i]->la_start > max_distance) break;
			if (cluster_t::can_merge(clusters[i], clusters[j], max_distance)) {
				ccps.push_back(cc_distance_t(clusters[i], clusters[j]));
			}
		}
		if (ccps.size() <= max_cluster_size) {
			for (cc_distance_t& ccp : ccps) pq.push(ccp);
		} else { // very large cluster - probably abnormal region - mark all involved pairs to be discarded
			clusters[i]->used = true;
			for (cc_distance_t& ccp : ccps) {
				ccp.c2->used = true;
			}
		}
	}

	std::multimap<int, std::shared_ptr<cluster_t>> clusters_map;
    for (std::shared_ptr<cluster_t> c : clusters) {
        if (c == NULL || c->used) continue;
        clusters_map.insert(std::make_pair(c->la_start, c));
        clusters_map.insert(std::make_pair(c->la_end, c));
    }

	while (!pq.empty()) {
        cc_distance_t ccd = pq.top();
        pq.pop();

        if (ccd.c1->used || ccd.c2->used) continue;

        auto new_cluster = std::shared_ptr<cluster_t>(cluster_t::merge(ccd.c1, ccd.c2));
        int parent_id = uf->merge(ccd.c1->id, ccd.c2->id);
        new_cluster->id = parent_id;
        clusters[parent_id] = new_cluster;
        clusters[ccd.c1->id+ccd.c2->id-parent_id] = NULL;

        ccd.c1->used = true;
        remove_cluster_from_mm(clusters_map, ccd.c1);
        ccd.c2->used = true;
        remove_cluster_from_mm(clusters_map, ccd.c2);

        auto end = clusters_map.upper_bound(new_cluster->la_end + max_distance);
        for (auto map_it = clusters_map.lower_bound(new_cluster->la_start - max_distance);
             map_it != end; map_it++) {
            if (!map_it->second->used && cluster_t::can_merge(new_cluster, map_it->second, max_distance)) {
                pq.push(cc_distance_t(new_cluster, map_it->second));
            }
        }
        clusters_map.insert(std::make_pair(new_cluster->la_start, new_cluster));
        clusters_map.insert(std::make_pair(new_cluster->la_end, new_cluster));
	}

	for (int i = 0; i < clusters.size(); i++) {
		if (clusters[i] != NULL && clusters[i]->used) {
			clusters[i] = NULL;
		}
	}

	if (assign_reads_to_clusters) {
		// for each set of reads, make a cluster-vector
		for (int i = 0; i < n_reads; i++) {
			int parent_id = uf->find(i);
			if (clusters[parent_id] != NULL) {
				clusters[parent_id]->reads.push_back(reads[i]);
			}
		}
	}
	delete uf;
}

#endif // CLUSTERING_UTILS_H
