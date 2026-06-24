#include <algorithm>
#include <vector>

struct SegTree {

    typedef u_int8_t T; // type of the segment tree values
    
    struct Node { T mx = 0, lazy = 0; };
    int n;
    std::vector<Node> st;

    SegTree(int n) : n(n), st(4*n) {}

    void apply(int p, T v) {
        st[p].mx   += v;
        st[p].lazy += v;
    }
    void push(int p) {             // propagate one level down
        if (st[p].lazy) {
            apply(p<<1, st[p].lazy);
            apply(p<<1|1, st[p].lazy);
            st[p].lazy = 0;
        }
    }

    void range_add(int p, int l, int r, int ql, int qr, T val) {
        if (qr < l || r < ql) return;
        if (ql <= l && r <= qr) { apply(p,val); return; }
        push(p);
        int mid = (l+r)/2;
        range_add(p<<1,   l,   mid, ql, qr, val);
        range_add(p<<1|1, mid+1,r,  ql, qr, val);
        st[p].mx = std::max(st[p<<1].mx, st[p<<1|1].mx);
    }
    // External wrapper: inclusive indices
    void add(int ql, int qr, T val){ range_add(1,0,n-1,ql,qr,val); }

    bool range_has_ge(int p,int l,int r,int ql,int qr,T m) {
        if (qr < l || r < ql || st[p].mx < m) return false;
        if (ql <= l && r <= qr) return true;
        push(p);
        int mid = (l+r)/2;
        return range_has_ge(p<<1,l,mid,ql,qr,m) ||
               range_has_ge(p<<1|1,mid+1,r,ql,qr,m);
    }
    bool any_ge(int ql,int qr,T m){
        return range_has_ge(1,0,n-1,ql,qr,m);
    }
};
