#ifndef DYNAMIC_INTERVAL_TREE_H
#define DYNAMIC_INTERVAL_TREE_H

#include <algorithm>
#include <cstdint>
#include <limits>
#include <vector>

template<typename Coord, typename T>
class DynamicIntervalTree {
    struct Node {
        Coord start, end, max_end;
        uint64_t key_id, priority;
        T value;
        Node* left;
        Node* right;

        Node(Coord start, Coord end, uint64_t key_id, uint64_t priority, T value)
            : start(start), end(end), max_end(end), key_id(key_id), priority(priority), value(value),
              left(nullptr), right(nullptr) {}
    };

    Node* root = nullptr;
    uint64_t next_key_id = 0;

    static uint64_t mix64(uint64_t x) {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }

    static bool key_less(Coord start1, Coord end1, uint64_t key_id1,
                         Coord start2, Coord end2, uint64_t key_id2) {
        if (start1 != start2) return start1 < start2;
        if (end1 != end2) return end1 < end2;
        return key_id1 < key_id2;
    }

    static Coord max_end(Node* node) {
        return node == nullptr ? std::numeric_limits<Coord>::min() : node->max_end;
    }

    static void update(Node* node) {
        if (node == nullptr) return;
        node->max_end = std::max(node->end, std::max(max_end(node->left), max_end(node->right)));
    }

    static Node* rotate_right(Node* node) {
        Node* left = node->left;
        node->left = left->right;
        left->right = node;
        update(node);
        update(left);
        return left;
    }

    static Node* rotate_left(Node* node) {
        Node* right = node->right;
        node->right = right->left;
        right->left = node;
        update(node);
        update(right);
        return right;
    }

    static Node* insert(Node* root, Node* node) {
        if (root == nullptr) return node;

        if (key_less(node->start, node->end, node->key_id, root->start, root->end, root->key_id)) {
            root->left = insert(root->left, node);
            if (root->left->priority > root->priority) {
                root = rotate_right(root);
            }
        } else {
            root->right = insert(root->right, node);
            if (root->right->priority > root->priority) {
                root = rotate_left(root);
            }
        }

        update(root);
        return root;
    }

    static void query(Node* node, Coord ql, Coord qr, std::vector<T>& out) {
        if (node == nullptr || node->max_end < ql) return;

        if (node->left != nullptr && node->left->max_end >= ql) {
            query(node->left, ql, qr, out);
        }

        if (node->start <= qr && node->end >= ql) {
            out.push_back(node->value);
        }

        if (node->start <= qr) {
            query(node->right, ql, qr, out);
        }
    }

    static void destroy(Node* node) {
        if (node == nullptr) return;
        destroy(node->left);
        destroy(node->right);
        delete node;
    }

public:
    ~DynamicIntervalTree() {
        destroy(root);
    }

    DynamicIntervalTree() = default;

    DynamicIntervalTree(const DynamicIntervalTree&) = delete;
    DynamicIntervalTree& operator=(const DynamicIntervalTree&) = delete;

    void insert(Coord start, Coord end, T value) {
        uint64_t key_id = next_key_id++;
        root = insert(root, new Node(start, end, key_id, mix64(key_id), value));
    }

    std::vector<T> query(Coord start, Coord end) const {
        std::vector<T> out;
        query(root, start, end, out);
        return out;
    }
};

#endif /* DYNAMIC_INTERVAL_TREE_H */
