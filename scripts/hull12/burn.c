#include <stdio.h>
#include <string.h>
#include "ico_data.h"
/* count edge subsets D of the icosahedral graph, invariant under automorphism g, with max degree(D) <= 2, |D|<=MAXD */
#define MAXD 13
int ncyc; int cyc[30][30]; int clen[30];
long long cnt[MAXD+1];
int deg[12];
void dfs(int c,int size){
  if(c==ncyc){ cnt[size]++; return; }
  dfs(c+1,size);
  if(size+clen[c]>MAXD) return;
  int ok=1;
  for(int k=0;k<clen[c];k++){int e=cyc[c][k];deg[EU[e]]++;deg[EV[e]]++;}
  for(int k=0;k<clen[c];k++){int e=cyc[c][k];if(deg[EU[e]]>2||deg[EV[e]]>2)ok=0;}
  if(ok) dfs(c+1,size+clen[c]);
  for(int k=0;k<clen[c];k++){int e=cyc[c][k];deg[EU[e]]--;deg[EV[e]]--;}
}
int main(){
  double tot[MAXD+1]; long long lab[MAXD+1]; memset(tot,0,sizeof tot);
  for(int a=0;a<NA;a++){
    int seen[30]={0}; ncyc=0;
    for(int e=0;e<30;e++) if(!seen[e]){ int f=e; clen[ncyc]=0; while(!seen[f]){seen[f]=1; cyc[ncyc][clen[ncyc]++]=f; f=PERM[a][f];} ncyc++; }
    memset(cnt,0,sizeof cnt); memset(deg,0,sizeof deg); dfs(0,0);
    for(int s=0;s<=MAXD;s++){ tot[s]+=cnt[s]; }
    int isid=1; for(int e=0;e<30;e++) if(PERM[a][e]!=e) isid=0;
    if(isid) for(int s=0;s<=MAXD;s++) lab[s]=cnt[s];
  }
  for(int s=0;s<=MAXD;s++) printf("deleted %d kept %d labeled %lld orbits %.1f\n",s,30-s,lab[s],tot[s]/NA);
}
