#ifndef SURVEYOR_SSW_SCORING_IDENTITY_H
#define SURVEYOR_SSW_SCORING_IDENTITY_H

#include <map>
#include <mutex>
#include "exact_cache.h"

namespace surveyor_cache {

inline uint64_t ssw_scoring_identity(const int8_t* matrix, int matrix_size, const int8_t* translation, int translation_size, int target_size) {
	key_builder_t key;
	key.add(matrix_size);
	key.add_bytes(reinterpret_cast<const char*>(matrix), matrix_size*matrix_size);
	key.add_bytes(reinterpret_cast<const char*>(translation), translation_size);
	key.add(target_size);
	static std::map<std::string, uint64_t> identities;
	static std::mutex mutex;
	std::lock_guard<std::mutex> lock(mutex);
	return identities.emplace(key.take(), identities.size()+1).first->second;
}

}

#endif
