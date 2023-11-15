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
  int idx;
  uintE* dis;
  uintE* prime;
  uintE batch_size;
  ORDER(int _idx, uintE* _dis, uintE* _prime, uintE _batch_size) :
    idx(_idx), dis(_dis), prime(_prime), batch_size(_batch_size)
  {}
  inline bool update (uintE s, uintE d, intE edgeLen){
    //printf("ni\n");
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

  inline bool updateAtomic (uintE s, uintE d, intE edgeLen){
    int up = 0;
    uintE ww = edgeLen / 1000;
    uintE pro = edgeLen % 1000;
    if (pro == 0) {
      pro = 1000;
      ww--;
    }
    float pp = float(pro) / 1000;

    int k = idx / batch_size;
    float _r1;
    float _r2;
    float _r3 = HybridTaus(LCGStep(s + d,19993,19997) + LCGStep(s + d,30001,30007), s + d, s + d, 3);
    int pass_num = int(pp * batch_size);
    float res = pp * batch_size - pass_num;
    //printf("%.3f,%.3f,%.3f\n",r1,r2,r3);
    _r1 = (_r3*197) - int(_r3*197);
    _r2 = (_r3*301) - int(_r3*301);
    _r3 = (_r3*19) - int(_r3*19);
    _r3 = _r3*(k + 1) - int(_r3*(k + 1));
    float r1 = _r1;
    float r2 = _r2;
    float r3 = _r3;
    if (r3 < res) pass_num++;
    r1 = r1 * (k + 1) - int(r1 * (k + 1));
    r2 = r2 * (k + 1) - int(r2 * (k + 1));
          
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
  inline bool cond (uintE d) { return cond_true(d); }
};

struct ORDER_UW {
  int idx;
  uintE* dis;
  uintE* prime;
  uintE batch_size;
  ORDER_UW(int _idx, uintE* _dis, uintE* _prime, uintE _batch_size) :
    idx(_idx), dis(_dis), prime(_prime), batch_size(_batch_size)
  {}
  inline bool update (uintE s, uintE d, intE edgeLen){
    //printf("ni\n");
    int up = 0;
    uintE ww = 1;
    uintE pro = edgeLen % 1000;
    if (pro == 0) {
      pro = 1000;
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

  inline bool updateAtomic (uintE s, uintE d, intE edgeLen){
    int up = 0;
    uintE ww = 1;
    uintE pro = edgeLen % 1000;
    if (pro == 0) {
      pro = 1000;
    }
    float pp = float(pro) / 1000;

    int k = idx / batch_size;
    float _r1;
    float _r2;
    float _r3 = HybridTaus(LCGStep(s + d,19993,19997) + LCGStep(s + d,30001,30007), s + d, s + d, 3);
    int pass_num = int(pp * batch_size);
    float res = pp * batch_size - pass_num;
    //printf("%.3f,%.3f,%.3f\n",r1,r2,r3);
    _r1 = (_r3*197) - int(_r3*197);
    _r2 = (_r3*301) - int(_r3*301);
    _r3 = (_r3*19) - int(_r3*19);
    _r3 = _r3*(k + 1) - int(_r3*(k + 1));
    float r1 = _r1;
    float r2 = _r2;
    float r3 = _r3;
    if (r3 < res) pass_num++;
    r1 = r1 * (k + 1) - int(r1 * (k + 1));
    r2 = r2 * (k + 1) - int(r2 * (k + 1));
          
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
  inline bool cond (uintE d) { return cond_true(d); }
};
