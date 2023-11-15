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


struct Vert {
	int dis;
  int cons;
	struct Vert* next;
};

bool Insert(Vert** d,uintE dis, uintE cost) {
    Vert* temp = new(Vert);
    temp = *d;
    if ((*d == NULL) || (temp->dis > dis)) {
      Vert* sert = new(Vert);
      sert->dis = dis;
      sert->cons = cost;
      sert->next = temp;
      *d = sert;
      return true;      
    }
    while (temp != NULL) {
        if ((temp->dis <= dis) && (temp->cons <= cost)) return false;
        if ((temp->dis >= dis) && (temp->cons >= cost)) {
          temp->dis = dis;
          temp->cons = cost;
          return true;
        }
        if ((temp->next == NULL) || ((dis < temp->next->dis) && (cost > temp->next->cons))) {
            Vert* sert = new(Vert);
            sert->dis = dis;
            sert->cons = cost;
            sert->next = temp->next;
            temp->next = sert;
            return true;
        }
        temp = temp->next;
    }
}
struct adj {
	int degree;
	Vert* neighbor;
};

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
    uintE len = edgeLen / 100000;
    if ((edgeLen % 100000) == 0) len--;
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
    uintE len = edgeLen / 100000;
    if ((edgeLen % 100000) == 0) len--;
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
  ORDER(bool* _visited, int _idx, uintE* _dis) :
    visited(_visited), idx(_idx), dis(_dis)
  {}
  inline bool update (uintE s, uintE d, intE edgeLen){
    //printf("ni\n");
    if (visited[d]) {
    int up = 0;
    uintE ww = edgeLen / 100000;
    uintE pro = edgeLen % 100;
    if (pro == 0) {
      pro = 100;
      ww--;
    }
    float pp = float(pro) / 100;
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
    uintE ww = edgeLen / 100000;
    uintE pro = edgeLen % 100;
    if (pro == 0) {
      pro = 100;
      ww--;
    }
    float pp = float(pro) / 100;
    float p = HybridTaus(LCGStep(s,19993,19997) + LCGStep(d,30001,30007), s, d, idx*idx);
    p = (p*19) - int(p*19);
    if ((p<pp) && (dis[s] + ww < dis[d])) {
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
