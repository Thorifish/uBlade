#define WEIGHTED 1
#define OPENMP
#include <cmath>
#include "ligra.h"
#include "index_map.h"
#include "bucket.h"
#include "DS_nonest.h"

int tt = 0;

struct Visit_F {
  array_imap<uintE> dists;
  uintE* distances;
  uintE* prime;
  uintE batch_size;
  uintE k;
  uintE delta;
  uintE* idx;
  uintE id;

  Visit_F(uintE* _idx, array_imap<uintE>& _dists, uintE* _distances, uintE* _prime, uintE _delta, uintE _id, uintE _batch_size, uintE _k) : idx(_idx), dists(_dists), distances(_distances), prime(_prime), delta(_delta), id(_id), batch_size(_batch_size), k(_k) {}

  inline Maybe<uintE> update(const uintE& s, const uintE& d, const intE& w) {    
    uintE oval;
    uintE dist, n_dist;
    uintE ww = w / 1000;
    uintE pro = (w % 1000);
    bool up = false;
    uintE re;
    float p = float(pro) / 1000;
    if (pro == 0) {
      ww--;
      p = 1.00;
    }

    float _r1;
    float _r2;
    float _r3 = HybridTaus(LCGStep(s + d,19993,19997) + LCGStep(s + d,30001,30007), s + d, s + d, 3);
    int pas_num = int(p * batch_size);
    float res = p * batch_size - pas_num;
    _r1 = (_r3*197) - int(_r3*197);
    _r2 = (_r3*301) - int(_r3*301);
    _r3 = (_r3*19) - int(_r3*19);
    for (int i=0;i<batch_size;i++) {
    if (distances[s*batch_size + i] != INT_E_MAX) {
      uintE re_dist = dists.s[d];
      oval = distances[d*batch_size + i];
      dist = oval | TOP_BIT;
      n_dist = (distances[s*batch_size + i] | TOP_BIT) + ww;
        if (n_dist < dist) {
          //c++;
          uintE B = idx[k*batch_size + i] / batch_size;
          float r1 = _r1;
          float r2 = _r2;
          float r3 = _r3;
          int pass_num = pas_num;
          r1 = r1 * (B + 1) - int(r1 * (B + 1));
          r2 = r2 * (B + 1) - int(r2 * (B + 1));
          r3 = r3 * (B + 1) - int(r3 * (B + 1));
          uintE str = int(r1 * batch_size) % batch_size;
          uintE step = int(r2 * 7) % 7;
          uintE ii = idx[k*batch_size + i] % batch_size;
          bool pass = false;
          if (ii < str) ii = ii + batch_size;
          if (r3 < res) pass_num++;
          while ((((ii - str) / prime[step]) < pass_num)) {
            if (((ii - str) % prime[step]) == 0) pass = true;
            ii = ii + batch_size;
          }
          if (pass) {
            if (!(re_dist & TOP_BIT)) { // First visitor
            //if (n_dist < (re_dist | TOP_BIT)) {
              uintE ndist = n_dist;
              uintE nbkt = (((n_dist & INT_E_MAX)) / delta);
              uintE pbkt = (((re_dist & INT_E_MAX)) / delta);
              if ((nbkt > pbkt) && (pbkt >= id)) ndist = re_dist | TOP_BIT; 
              if (CAS(&(dists[d]), re_dist, ndist)) { 
                CAS(&(distances[d*batch_size + i]), oval , n_dist);
                up = true;
                re = re_dist;
              }
            }
            writeMin(&(dists[d]), n_dist);
            CAS(&(distances[d*batch_size + i]), oval , n_dist);
          }
        }
    }
    }
    if (up) return Maybe<uintE>(re);
    return Maybe<uintE>();
  }

  inline Maybe<uintE> updateAtomic(const uintE& s, const uintE& d, const intE& w) {
    uintE oval;
    uintE dist, n_dist;
    uintE ww = w / 1000;
    uintE pro = (w % 1000);
    bool up = false;
    uintE re;
    float p = float(pro) / 1000;
    if (pro == 0) {
      ww--;
      p = 1.00;
    }

    float _r1;
    float _r2;
    float _r3 = HybridTaus(LCGStep(s + d,19993,19997) + LCGStep(s + d,30001,30007), s + d, s + d, 3);
    int pas_num = int(p * batch_size);
    float res = p * batch_size - pas_num;
    _r1 = (_r3*197) - int(_r3*197);
    _r2 = (_r3*301) - int(_r3*301);
    _r3 = (_r3*19) - int(_r3*19);
    for (int i=0;i<batch_size;i++) {
    if (distances[s*batch_size + i] != INT_E_MAX) {
      uintE re_dist = dists.s[d];
      oval = distances[d*batch_size + i];
      dist = oval | TOP_BIT;
      n_dist = (distances[s*batch_size + i] | TOP_BIT) + ww;
        if (n_dist < dist) {
          //c++;
          uintE B = idx[k*batch_size + i] / batch_size;
          float r1 = _r1;
          float r2 = _r2;
          float r3 = _r3;
          int pass_num = pas_num;
          r1 = r1 * (B + 1) - int(r1 * (B + 1));
          r2 = r2 * (B + 1) - int(r2 * (B + 1));
          r3 = r3 * (B + 1) - int(r3 * (B + 1));
          uintE str = int(r1 * batch_size) % batch_size;
          uintE step = int(r2 * 7) % 7;
          uintE ii = idx[k*batch_size + i] % batch_size;
          bool pass = false;
          if (ii < str) ii = ii + batch_size;
          if (r3 < res) pass_num++;
          while ((((ii - str) / prime[step]) < pass_num)) {
            if (((ii - str) % prime[step]) == 0) pass = true;
            ii = ii + batch_size;
          }
          if (pass) {
            if (!(re_dist & TOP_BIT)) { // First visitor
            //if (n_dist < (re_dist | TOP_BIT)) {
              uintE ndist = n_dist;
              uintE nbkt = (((n_dist & INT_E_MAX)) / delta);
              uintE pbkt = (((re_dist & INT_E_MAX)) / delta);
              if ((nbkt > pbkt) && (pbkt >= id)) ndist = re_dist | TOP_BIT; 
              if (CAS(&(dists[d]), re_dist, ndist)) { 
                CAS(&(distances[d*batch_size + i]), oval , n_dist);
                up = true;
                re = re_dist;
              }
            }
            writeMin(&(dists[d]), n_dist);
            CAS(&(distances[d*batch_size + i]), oval , n_dist);
          }
        }
        //}
    }
    }
    if (up) return Maybe<uintE>(re);
    //}
    return Maybe<uintE>();
  }

  inline bool cond(const uintE& d) const { return true; }
};

template <class vertex>
void DeltaStepping(graph<vertex>& G, uintE* di, uintE* reli, uintE* idx, uintE src, uintE delta, uintE k, uintE* prime, uintE batch_size, size_t num_buckets=128) {
  auto V = G.V; size_t n = G.n, m = G.m;
  auto dists = array_imap<uintE>(n, [&] (size_t i) { return INT_E_MAX; });
  dists[src] = 0;

  //uintE batch_size = 50;
  uintE* distances = newA(uintE,batch_size*n);
  for (int i=0;i<n;i++) {
    for (int j=0;j<batch_size;j++) {
      distances[i*batch_size + j] = INT_E_MAX;
    }
  }
  for (int j=0;j<batch_size;j++) {
    distances[src*batch_size + j] = 0;
  }


  auto get_bkt = [&] (const uintE& dist) -> const uintE {
    return (dist >= INT_E_MAX) ? UINT_E_MAX : (dist / delta); };
  auto get_ring = [&] (const size_t& v) -> const uintE {
    auto d = dists[v];
    return (d >= INT_E_MAX) ? UINT_E_MAX : (d / delta); };
  auto b = make_buckets(n, get_ring, increasing, num_buckets);

  auto apply_f = [&] (const uintE v, uintE& oldDist) -> void {
    uintE newDist;
    for (int i=0;i<batch_size;i++) {
      newDist = distances[v*batch_size + i] & VAL_MASK;
      distances[v*batch_size + i] = newDist;
    } 
    newDist = dists.s[v] & VAL_MASK;
    dists.s[v] = newDist; // Remove the TOP_BIT in the distance.
    // Compute the previous bucket and new bucket for the vertex.
    uintE prev_bkt = get_bkt(oldDist), new_bkt = get_bkt(newDist);
    bucket_dest dest = b.get_bucket(prev_bkt, new_bkt);
    oldDist = dest; // write back
  };

  auto bkt = b.next_bucket();

  while (bkt.id != b.null_bkt) {
    uintE id = bkt.id;
    auto active = bkt.identifiers;

    // The outvisitedput of the edgeMap is a vertexSubsetData<uintE> where the value
    // stored with each vertex is its original distance in this round
  
    auto res = edgeMapData<uintE>(G, active, Visit_F(idx,dists,distances,prime,delta,id,batch_size,k), G.m/20, sparse_no_filter | dense_forward);
    vertexMap(res, apply_f);
    
    if (res.dense()) {
        b.update_buckets(res.get_fn_repr(), n);
    } else {
        b.update_buckets(res.get_fn_repr(), res.size());
    }
    res.del(); active.del();
    bkt = b.next_bucket();
  }
  uintE temp;
  for (int i = 0;i<batch_size;i++) {
    temp = 0;
    while (di[temp]!= distances[tt*batch_size + i]) {
      temp++;
      if (di[temp] == 0) di[temp] = distances[tt*batch_size + i];
    }
    reli[temp]++;
  }
  free(distances);
}

template <class vertex>
void Compute(graph<vertex>& GA, commandLine P, uintE src, uintE tar) {
  uintE delta = P.getOptionLongValue("-delta",1);
  size_t num_buckets = P.getOptionLongValue("-nb", 128);
  uintE batch_size = P.getOptionLongValue("-batch",1);
  if (num_buckets != (1 << pbbs::log2_up(num_buckets))) {
    cout << "Please specify a number of buckets that is a power of two" << endl;
    exit(-1);
  }

  uintE n = GA.n;
  uintE total_sample = P.getOptionLongValue("-sample",100);
  
  //bool ord = P.getOptionValue("-order");
  bool ord = true;

  long sampleSize = 1;
  uintE* starts = newA(uintE,1);
  tt = tar;

  uintE* prime = newA(uintE,7);
  prime[0] = 1;
  prime[1] = 3;
  prime[2] = 7;
  prime[3] = 11;
  prime[4] = 13;
  prime[5] = 17;
  prime[6] = 19;

  uintE* dis = newA(uintE,n);
  uintE* ord_count = newA(uintE,total_sample);
  uintE* ord_idx = newA(uintE,total_sample);
  for (int i=0;i<total_sample;i++) {
    ord_idx[i] = i;
    ord_count[i] = 0;
  }
  if (ord) {
  starts[0] = src;
  vertexSubset front(n,sampleSize,starts);
  vertexSubset FF = front;
  for (int i = 0;i<total_sample;i++) {
    FF = front;
    for (int j = 0;j<n;j++) {
      dis[j] = INT_E_MAX / 2;
    }
    dis[src] = 0;
    for (int j = 0;j<2;j++) {
      vertexSubset output = edgeMap(GA,FF,ORDER(i,dis,prime,batch_size));
      FF = output;
      ord_count[i] = ord_count[i] + FF.size();
    }
  }
  }

  int temp = 0;

  if (ord) {
  for (int i=0;i<total_sample;i++) {
    for (int j=0;j<total_sample-i-1;j++) {
      if (ord_count[j] > ord_count[j+1]) {
        temp = ord_count[j];
        ord_count[j] = ord_count[j + 1];
        ord_count[j + 1] = temp;
        temp = ord_idx[j];
        ord_idx[j] = ord_idx[j + 1];
        ord_idx[j + 1] = temp;
      }
    }
  }
  }

  uintE* di = newA(uintE,total_sample + 1);
  uintE* reli = newA(uintE,total_sample + 1);
  for (int i=0;i<total_sample + 1;i++) {
    reli[i] = 0;
    di[i] = 0;
  }
  di[0] = INT_E_MAX;

  for (int i=0;i<total_sample/batch_size;i++) {
    DeltaStepping(GA, di, reli, ord_idx, src, delta, i, prime, batch_size, num_buckets);
    for (int j = 0;j<n;j++) {
    }
  }
 
  uintE cc = 0;
  uintE rel = 0;
  uintE total_rel = 0;
  for (int i=1;i<total_sample + 1;i++) {
    total_rel = total_rel + reli[i]; 
    if (reli[i] > rel) {
      rel = reli[i];
      cc = di[i];
    } else if (reli[i] == rel) {
      if (di[i] < cc) {
        cc = di[i];
      }
    }
  }
  printf("most probable distance: %d, reliability: %d\n",cc,rel);
  free(di);
  free(reli);
}
