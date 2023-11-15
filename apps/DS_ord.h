#include "ligra.h"

constexpr uintE TOP_BIT = ((uintE)INT_E_MAX) + 1;
constexpr uintE VAL_MASK = INT_E_MAX;


unsigned TausStep(unsigned z, int S1, int S2, int S3, unsigned M) {
	unsigned b=(((z << S1) ^ z) >> S2);
	return (((z & M) << S3) ^ b);
}

unsigned LCGStep(unsigned z, unsigned A, unsigned C) {
	return (A*z+C);
}

float HybridTaus(unsigned z1, unsigned z2, unsigned z3, unsigned z4 ) {
  return 2.3283064365387e-10 * ( TausStep(z1, 13, 19, 12, 4294967294UL) ^ TausStep(z2, 2, 25, 4, 4294967288UL) ^ TausStep(z3, 3, 11, 17, 4294967280UL) ^ LCGStep(z4, 1664525, 1013904223UL) );
}

struct BFS_S {
  bool* visit;
  BFS_S(bool* _visit) :
    visit(_visit)
  {}
  inline bool update (uintE s, uintE d, intE edgeLen){
    int up = 0;
    if (!visit[d]) {
      up = 1;
      visit[d] = true;
    }
    return up;
  }

  inline bool updateAtomic (uintE s, uintE d, intE edgeLen){
    int up = 0;
    if (!visit[d]) {
      up = 1;
      visit[d] = true;
    }
    return up;
  }
  inline bool cond (uintE d) { return cond_true(d); }
};

////////////
struct BFS_T {
  bool* visit;
  bool* visited;
  bool* uped;
  uintE* shortest;
  uintE *longest;
  BFS_T(bool* _visited, bool* _visit, bool* _uped, uintE* _shortest, uintE *_longest) :
    visited(_visited), visit(_visit), uped(_uped), shortest(_shortest), longest(_longest)
  {}
  inline bool update (uintE s, uintE d, intE edgeLen){
    int up = 0;
    uintE len = edgeLen / 1000;
    if ((edgeLen % 1000) == 0) len--;
    if ((visit[d]) && (shortest[s] + len < shortest[d])) {
      shortest[d] = shortest[s] + len;
      visited[d] = true;
      if (!uped[d]) {
        uped[d] = true;
        up = 1;
      }
    }
    //printf("%d,%ld\n",len,longest);
    if (len > *longest) {
      *longest = len;
    }
    return up;
  }

  inline bool updateAtomic (uintE s, uintE d, intE edgeLen){
    int up = 0;
    uintE len = edgeLen / 1000;
    if ((edgeLen % 1000) == 0) len--;
    if ((visit[d]) && (shortest[s] + len < shortest[d])) {
      uintE toWrite = shortest[s] + len;
      uintE olddis = shortest[d];
      if (CAS(&(shortest[d]), olddis, toWrite)) {
        visited[d] = true;
        if (!uped[d]) {
          uped[d] = true;
          up = 1;
        }
      }
    }
    if (len > *longest) {
      uintE toWrite = len;
      uintE olddis = *longest;
      CAS(longest, olddis, toWrite);
    }
    return up;
  }
  inline bool cond (uintE d) { return cond_true(d); }
};

///////////////////////

struct ORDER {
  bool* visited;
  int idx;
  uintE* dis;
  uintE* prime;
  uintE batch_size;
  ORDER(bool* _visited, int _idx, uintE* _dis, uintE* _prime, uintE _batch_size) :
    visited(_visited), idx(_idx), dis(_dis), prime(_prime), batch_size(_batch_size)
  {}
  inline bool update (uintE s, uintE d, intE edgeLen){
    //printf("ni\n");
    if (visited[d]) {
    int up = 0;
    uintE ww = edgeLen / 1000;
    uintE pro = edgeLen % 1000;
    if (pro == 0) {
      pro = 1000;
      ww--;
    }
    float pp = float(pro) / 1000;
    float p = HybridTaus(LCGStep(s,19993,19997) + LCGStep(d,30001,30007), s, d, idx*idx);
    p = (p*19) - int(p*19);
    if ((p<pp) && (dis[s] + ww < dis[d])) {
      dis[d] = dis[s] + ww;
      up = 1;
    }
    return up;
    }
    return 0;
  }

  inline bool updateAtomic (uintE s, uintE d, intE edgeLen){
    if (visited[d]) {
    int up = 0;
    uintE ww = edgeLen / 1000;
    uintE pro = edgeLen % 1000;
    if (pro == 0) {
      pro = 1000;
      ww--;
    }
    float pp = float(pro) / 1000;

    int k = idx / batch_size;
    float _r1 = HybridTaus(LCGStep(s,193,197) + LCGStep(d,301,307), s + d, s + d, 1);
    float _r2 = HybridTaus(LCGStep(s,3,7) + LCGStep(d,1,7), s, d, 2);
    float _r3 = HybridTaus(LCGStep(s + d,19993,19997) + LCGStep(s + d,30001,30007), s + d, s + d, 3);
    int pass_num = int(pp * batch_size);
    float res = pp * batch_size - pass_num;
    //printf("%.3f,%.3f,%.3f\n",r1,r2,r3);
    _r1 = (_r1*23) - int(_r1*23);
    _r2 = (_r2*11) - int(_r2*11);
    _r3 = (_r3*19) - int(_r3*19);
    _r3 = _r3*(k + 1) - int(_r3*(k + 1));
    float r1 = _r1;
    float r2 = _r2;
    float r3 = _r3;
    if (r3 < res) pass_num++;
    for (int i=0;i<k + 1;i++) {
      r1 = r1 - int(r1);
      r1 = r1 + _r1;
      r2 = r2 - int(r2);
      r2 = r2 + _r2;
    }
    //if (pass_num == 20) printf("%.3f\n",p);
    //printf("%d\n",pass_num);
    //printf("%d,%d,%.3f\n",s,d,p);
    uintE str = int(r1 * batch_size) % batch_size;
    uintE step = int(r2 * 7) % 7;
    //int c = 0;

    bool pass = false;
    uintE id = idx % batch_size;
    if (id < str) id = id + batch_size;
    while ((((id - str) / prime[step]) < pass_num)) {
      if (((id - str) % prime[step]) == 0) pass = true;
      id = id + batch_size;
    }
    
    
    if ((pass) && (dis[s] + ww < dis[d])) {
      uintE temp = dis[d];
      if (CAS(&(dis[d]), temp, dis[s] + ww)) {
        up = 1;
      }
    }
    return up;
    }
    return 0;
  }
  inline bool cond (uintE d) { return cond_true(d); }
};






struct Visit_F2 {
  array_imap<uintE> dists;
  uintE* distances;
  uintE batch_size;
  uintE k;
  uintE tar;
  uintE delta;
  uintE* share;
  uintE* parent;
  bool* mark;
  uintE* idx;
  uintE id;

  bool* visited;
  uintE* shortest;
  Visit_F2(uintE* _idx, array_imap<uintE>& _dists, uintE* _distances, uintE* _parent, uintE* _share, bool* _mark, uintE _delta, uintE _id, uintE _batch_size, uintE _k, uintE _tar, bool* _visited, uintE* _shortest) : idx(_idx), dists(_dists), distances(_distances), parent(_parent), share(_share), mark(_mark), delta(_delta), id(_id), batch_size(_batch_size), k(_k), tar(_tar), visited(_visited), shortest(_shortest) {}

  inline Maybe<uintE> update(const uintE& s, const uintE& d, const intE& w) {
    if (visited[d]) {
      printf("hi\n");
      exit(0);
    uintE oval;
    uintE dist, n_dist;
    uintE ww = w / 1000;
    uintE pro = (w % 1000);
    float p = float(pro) / 1000;
    if (pro == 0) {
      ww--;
      p = 1.00;
    }
    for (int i=0;i<batch_size;i++) {
      if (distances[s*batch_size + i] != INT_E_MAX) {
        uintE est = distances[tar*batch_size + i] & INT_E_MAX;
        //uintE oval = dists.s[d];
        uintE re_dist = dists[d];
        oval = distances[d*batch_size + i];
        dist = oval | TOP_BIT;
        //n_dist = (dists.s[s] | TOP_BIT) + w;
        n_dist = (distances[s*batch_size + i] | TOP_BIT) + ww;
        if (((n_dist + shortest[d]) & INT_E_MAX) < est) {
          float pp = HybridTaus(LCGStep(s + d,19993,19997) + LCGStep(s + d,30001,30007), s + d, s + d, i*i);
          pp = (pp*19) - int(pp*19);
          if ((n_dist < dist) && (pp < p)) {
            if (!(re_dist & TOP_BIT)) { // First visitor
              dists[d] = n_dist;
              distances[d*batch_size + i] = n_dist;
              return Maybe<uintE>(re_dist);
            }
            dists[d] = n_dist;
            distances[d*batch_size + i] = n_dist;
          }
        }
      }
    }
    }
    return Maybe<uintE>();
  }

  inline Maybe<uintE> updateAtomic(const uintE& s, const uintE& d, const intE& w) {
    //uintE oval = dists.s[d];
    //uintE dist = oval | TOP_BIT;
    //uintE n_dist = (dists.s[s] | TOP_BIT) + w;
    //if ((s == 2684915) && (idx[k*batch_size + 0] == 38) && (d == 1299452)) printf("yes\n"); 
    if (visited[d]) {
    uintE oval;
    uintE dist, n_dist;
    uintE ww = w / 1000;
    uintE pro = (w % 1000);
    bool up = false;
    uintE re;
    float p = float(pro) / 1000;
    bool mar = false;
    if (pro == 0) {
      ww--;
      p = 1.00;
    }
    //printf("%d\n",ww);
    //if (n_dist < dist) {
    //  if (!(oval & TOP_BIT) && CAS(&(dists[d]), oval, n_dist)) { // First visitor
    //    return Maybe<uintE>(oval);
    //  }
    //  writeMin(&(dists[d]), n_dist);
    //}
    for (int i=0;i<batch_size;i++) {
      if (distances[s*batch_size + i] != INT_E_MAX) {
        if (!mark[s]) {
          share[omp_get_thread_num()]++;
          mar = true;
          //printf("%d,%d\n",s,omp_get_thread_num());
        }
        uintE est = distances[tar*batch_size + i] & INT_E_MAX;
        uintE re_dist = dists.s[d];
        uintE ii = idx[(k*batch_size + i)];
        oval = distances[d*batch_size + i];
        dist = oval | TOP_BIT;
        n_dist = (distances[s*batch_size + i] | TOP_BIT) + ww;
        //if (idx[(k*batch_size + i)] == 66) {
        //    if ((s == 1008823) && (d == 1692325)) printf("->%d:%d,%d,%d,%d\n",d,distances[s*batch_size + i],shortest[d],est,n_dist);
        //}
        if (((n_dist + shortest[d]) & INT_E_MAX) < est) {
          float pp = HybridTaus(LCGStep(s + d,19993,19997) + LCGStep(s + d,30001,30007), s + d, s + d, ii*ii);
          pp = (pp*19) - int(pp*19);
          //pp = 0.03;
          //pp = 0;
          //int ra = rand() % 1000;
          //float pp = float(ra) / 1000;
          if ((n_dist < dist) && (pp < p)) {
            if (!(re_dist & TOP_BIT)) { // First visitor
              //if (n_dist < (re_dist | TOP_BIT)) {
                uintE ndist = n_dist;
                uintE nbkt = (((n_dist & INT_E_MAX) + shortest[d]) / delta);
                uintE pbkt = (((re_dist & INT_E_MAX) + shortest[d]) / delta);
                if ((nbkt > pbkt) && (pbkt >= id)) ndist = re_dist | TOP_BIT; 
                if (CAS(&(dists[d]), re_dist, ndist)) { 
                  CAS(&(distances[d*batch_size + i]), oval , n_dist);
                  //printf("%d->%d\n",re_dist,dists[d]);
                  up = true;
                  re = re_dist;
                  //return Maybe<uintE>(re_dist);
                }
              //} else {
                //if (CAS(&(dists[d]), re_dist, (re_dist | TOP_BIT))) { 
                //  CAS(&(distances[d*batch_size + i]), oval , n_dist);
                //  up = true;
                //  re = re_dist;
                //}
              //}
            }
            writeMin(&(dists[d]), n_dist);
            CAS(&(distances[d*batch_size + i]), oval , n_dist);
            //if (idx[(k*batch_size + i)] == 87) {
            //  if ((s == 1582490) && (d == 2684915)) printf("%d->%d:%d, %d+%d bkt:%d\n",s,d,distances[d*batch_size + i] & INT_E_MAX,distances[s*batch_size + i] & INT_E_MAX,ww,(dists[d])&INT_E_MAX);
            //  uintE ori = parent[d];
            //  CAS(&(parent[d]), ori , s);
            //}
            //printf("%d,%d,%d,%d\n",oval,n_dist,dist,distances[d*batch_size + i]);
          }
        }
      }
    }
    if (mar) {
      mark[s] = true;
      //printf("%d\n",s);
    }
    if (up) return Maybe<uintE>(re);
    }
    return Maybe<uintE>();
  }

  inline bool cond(const uintE& d) const { return true; }
};

struct nei {
	uintE v;
	nei* next;
};


struct MyWorld {
  nei* neighbor;

  MyWorld() {
    neighbor = NULL;
  }

  bool Is_in(uintE i) {
    nei* temp = new(nei);
    temp = neighbor;
    while (temp != NULL) {
      if (temp->v == i) {
        return true;
      }
      temp = temp->next;
    }
    return false;
  }

  void push(uintE i) {
    nei* temp = new(nei);
    temp->v = i;
    temp->next = neighbor;
    neighbor = temp; 
  }
};

