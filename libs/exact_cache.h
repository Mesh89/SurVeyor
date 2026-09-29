#ifndef SURVEYOR_EXACT_CACHE_H
#define SURVEYOR_EXACT_CACHE_H

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>

namespace surveyor_cache {

inline uint64_t& block_generation() {
	static thread_local uint64_t generation = 0;
	return generation;
}

inline void start_block() {
	++block_generation();
}

inline size_t cache_limit(const char* env_name, size_t default_limit) {
	const char* value = std::getenv(env_name);
	if (!value || !*value) return default_limit;
	size_t parsed = 0;
	for (const char* c = value; *c; c++) {
		if (*c < '0' || *c > '9') return default_limit;
		size_t digit = *c - '0';
		if (parsed > (std::numeric_limits<size_t>::max() - digit) / 10) return default_limit;
		parsed = parsed * 10 + digit;
	}
	return parsed;
}

class key_builder_t {
public:
	template<typename T>
	void add(const T& value) {
		static_assert(std::is_trivially_copyable<T>::value, "cache keys require trivially copyable values");
		key.append(reinterpret_cast<const char*>(&value), sizeof(value));
	}

	void add_bytes(const char* data, size_t length) {
		uint64_t stored_length = length;
		add(stored_length);
		if (length) key.append(data, length);
	}

	void add_string(const std::string& value) {
		add_bytes(value.data(), value.size());
	}

	std::string take() {
		return std::move(key);
	}

private:
	std::string key;
};

template<typename Value>
class exact_cache_t {
public:
	exact_cache_t(const char* env_name, size_t default_limit) : limit(cache_limit(env_name, default_limit)), generation(block_generation()) {}

	const Value* find(const std::string& key) {
		sync_block();
		if (limit == 0) return NULL;
		typename std::unordered_map<std::string, Value>::const_iterator it = entries.find(key);
		if (it == entries.end()) return NULL;
		return &it->second;
	}

	void store(std::string key, Value value) {
		sync_block();
		if (limit == 0) return;
		if (entries.size() >= limit) entries.clear();
		entries.emplace(std::move(key), std::move(value));
	}

private:
	void sync_block() {
		if (generation == block_generation()) return;
		entries.clear();
		generation = block_generation();
	}

	size_t limit;
	uint64_t generation;
	std::unordered_map<std::string, Value> entries;
};

}

#endif
