/******************************************************************
 * dense_ds_allH.cpp                 (C++14, portable, all H=2-6)
 * ----------------------------------------------------------------
 * Exact densest subgraph for   H = 2,3,4,5,6   (edge / h-clique).
 * Writes every result block into  tempoutput.txt  (truncated).
 *
 * Designed for small-to-moderate graphs (≤ 7 000 vertices) such as
 * the As-733 benchmark.  For larger graphs a core-based algorithm
 * is required.
 *
 * Build :  g++ -O2 -std=c++14 dense_ds_allH.cpp -o dense_ds
 * Usage :  dense_ds  graph.txt
 *          dense_ds  --all  folder_with_txt_or_gml
 ******************************************************************/
#include <bits/stdc++.h>
using namespace std;

/* ---------------- Dinic -------------------------------------- */
struct Dinic{
    struct E{int to,rev;double cap;};
    vector<vector<E>> G; vector<int> lvl,it;
    Dinic(int n=0){init(n);}
    void init(int n){G.assign(n,{}); lvl.resize(n); it.resize(n);}
    void addEdge(int u,int v,double c){
        G[u].push_back({v,(int)G[v].size(),c});
        G[v].push_back({u,(int)G[u].size()-1,0});
    }
    bool bfs(int s,int t){
        fill(lvl.begin(),lvl.end(),-1); queue<int>q; lvl[s]=0; q.push(s);
        while(!q.empty()){
            int v=q.front(); q.pop();
            for(auto&e:G[v]) if(e.cap>1e-9 && lvl[e.to]<0){
                lvl[e.to]=lvl[v]+1; q.push(e.to);
            }
        } return lvl[t]>=0;
    }
    double dfs(int v,int t,double f){
        if(v==t) return f;
        for(int &i=it[v]; i<(int)G[v].size(); ++i){
            auto &e=G[v][i];
            if(e.cap>1e-9 && lvl[v]<lvl[e.to]){
                double r=dfs(e.to,t,min(f,e.cap));
                if(r>1e-9){ e.cap-=r; G[e.to][e.rev].cap+=r; return r; }
            }
        } return 0;
    }
    double maxflow(int s,int t){
        const double INF=1e100; double F=0,aug;
        while(bfs(s,t)){ fill(it.begin(),it.end(),0);
            while((aug=dfs(s,t,INF))>1e-9) F+=aug; }
        return F;
    }
    vector<char> reach(int s){
        vector<char> vis(G.size(),0); queue<int>q; q.push(s); vis[s]=1;
        while(!q.empty()){ int v=q.front(); q.pop();
            for(auto&e:G[v]) if(e.cap>1e-9 && !vis[e.to]){
                vis[e.to]=1; q.push(e.to); } }
        return vis;
    }
};

/* ---------- naive (h-1)-clique enumeration (small graphs) ---- */
static void enumCliquesRec(const vector<vector<int>>&adj,
                           vector<int>&cand,vector<int>&cur,int k,
                           vector<vector<int>>&out){
    if((int)cur.size()==k){ out.push_back(cur); return; }
    while(!cand.empty()){
        int v=cand.back(); cand.pop_back();
        vector<int> next;
        for(int u:cand)
            if(find(adj[v].begin(),adj[v].end(),u)!=adj[v].end())
                next.push_back(u);
        cur.push_back(v);
        enumCliquesRec(adj,next,cur,k,out);
        cur.pop_back();
    }
}
static void enumCliques(const vector<vector<int>>&adj,int k,
                        vector<vector<int>>&out){
    /* enumerate all k-cliques, k ≤ 5 (because H ≤ 6) */
    vector<int> cand(adj.size()); iota(cand.begin(),cand.end(),0);
    vector<int> cur; enumCliquesRec(adj,cand,cur,k,out);
}

/* ---------- build flow for given H --------------------------- */
static void build_flow(int H,const vector<vector<int>>&adj,double a,Dinic&D)
{
    int n=adj.size(); int S=0, T=1;                    // temporary
    if(H==2){
        D.init(n+2);
        for(int u=0;u<n;++u){
            D.addEdge(S,u+2,(double)adj[u].size());
            D.addEdge(u+2,T,2.0*a);
            for(int v:adj[u]) if(u<v){
                D.addEdge(u+2,v+2,1.0); D.addEdge(v+2,u+2,1.0);
            }
        }
        return;
    }
    /* enumerate (H-1)-cliques */
    vector<vector<int>> inst;
    enumCliques(adj,H-1,inst);

    int m = inst.size();
    int baseV = 2;             // vertices start at 2
    int baseC = 2 + n;         // clique-instance nodes
    D.init(2 + n + m);         // S,T already counted

    for(int v=0;v<n;++v){
        D.addEdge(S, baseV+v, (double)adj[v].size());
        D.addEdge(baseV+v, T, H * a);                  // alpha*|VΨ|
    }
    for(int i=0;i<m;++i){
        int cnode = baseC + i;
        for(int v:inst[i]){
            D.addEdge(cnode, baseV+v, std::numeric_limits<double>::infinity());
            D.addEdge(baseV+v, cnode, 1.0);
        }
    }
}

/* ---------- exact densest subgraph for H --------------------- */
static vector<int> exactDS(int H,const vector<vector<int>>&adj,double&alpha){
    int n=adj.size(); double lo=0,hi=0;
    for(auto&nb:adj) if(nb.size()>hi) hi=nb.size();

    vector<int> best;
    auto feas=[&](double a)->bool{
        Dinic D; build_flow(H,adj,a,D);
        D.maxflow(0,1);
        auto vis=D.reach(0); bool any=false;
        for(int v=0;v<n;++v) if(vis[2+v]){ any=true; break; }
        if(any){ best.clear(); for(int v=0;v<n;++v) if(vis[2+v]) best.push_back(v); }
        return any;
    };
    while(hi-lo>1e-6){ double mid=0.5*(lo+hi); (feas(mid)?lo:hi)=mid; }
    alpha=lo; return best;
}

/* ---------- read edge list ----------------------------------- */
struct Graph{ vector<vector<int>> adj; long long m=0; };
static Graph readEL(const string&fn){
    ifstream fin(fn.c_str()); if(!fin) throw runtime_error("open "+fn);
    string line; size_t Nd=0; vector<pair<int,int>> E;
    while(getline(fin,line)){
        if(line.empty()) continue;
        if(line[0]=='#'){ size_t p=line.find("Nodes:"); if(p!=string::npos)
            Nd=(size_t)atoll(line.substr(p+6).c_str()); continue; }
        istringstream iss(line); int u,v; if(iss>>u>>v) E.push_back({u,v});
    }
    size_t n=Nd; for(auto&e:E) n=max(n,(size_t)max(e.first,e.second)+1);
    Graph G; G.adj.assign(n,{});
    for(auto&e:E){int u=e.first,v=e.second;if(u==v)continue;
        G.adj[u].push_back(v); G.adj[v].push_back(u); ++G.m;}
    return G;
}

/* ---------- recursive dir walk (C++14) ----------------------- */
#ifdef _WIN32
#include <windows.h>
static inline void low(string&s){for(char&c:s)c=(char)tolower(c);}
static void walk(const string&dir,vector<string>&out){
    string pat=dir.back()=='\\'||dir.back()=='/'?dir:dir+"\\";
    pat+="*"; WIN32_FIND_DATAA ffd; HANDLE h=FindFirstFileA(pat.c_str(),&ffd);
    if(h==INVALID_HANDLE_VALUE) return;
    do{ string n=ffd.cFileName; if(n=="."||n=="..") continue;
        string full=dir+"\\"+n;
        if(ffd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) walk(full,out);
        else{ string ext=n.substr(n.find_last_of('.')+1); low(ext);
              if(ext=="txt"||ext=="gml") out.push_back(full);}
    }while(FindNextFileA(h,&ffd)); FindClose(h);
}
#else
#include <dirent.h>
#include <sys/stat.h>
#include <fnmatch.h>
static void walk(const string&dir,vector<string>&out){
    DIR*dp=opendir(dir.c_str()); if(!dp) return; dirent*de;
    while((de=readdir(dp))){ string n=de->d_name;
        if(n=="."||n=="..") continue; string full=dir+"/"+n;
        struct stat st; if(stat(full.c_str(),&st)<0) continue;
        if(S_ISDIR(st.st_mode)) walk(full,out);
        else if(!fnmatch("*.txt",n.c_str(),FNM_CASEFOLD)||
                !fnmatch("*.gml",n.c_str(),FNM_CASEFOLD))
            out.push_back(full);
    } closedir(dp);
}
#endif
static void listFiles(const string&dir,vector<string>&out){
    out.clear(); walk(dir,out); sort(out.begin(),out.end());
}

/* ---------- write one block ---------------------------------- */
static void block(ofstream&f,int H,const string&fn,double a,
                  const vector<int>&sub,long long N,long long M,double ms){
    f<<"---------------------------------------------\n";
    f<<"Nodes: "<<N<<" Edges: "<<M<<"\n";
    f<<"H: "<<H<<"\n";
    f.setf(ios::fixed); f<<setprecision(6);
    f<<"Alpha: "<<a<<"\n";
    f<<"Nodes in CDS: "<<sub.size()<<"\nCDS:";
    for(int v:sub) f<<" "<<v;
    f<<"\nTime taken: "<<(long long)(ms+0.5)<<" ms\n";
}

/* ---------------- MAIN --------------------------------------- */
int main(int argc,char**argv){
    ios::sync_with_stdio(false);
    bool batch=false; string tgt;
    if(argc>=3 && (string(argv[1])=="--all"||string(argv[1])=="-all"))
        batch=true,tgt=argv[2];
    else if(argc>=2) tgt=argv[1];
    else{ cerr<<"usage: "<<argv[0]<<" file | --all dir\n"; return 1; }

    vector<string> files; if(batch){ listFiles(tgt,files);
        if(files.empty()){ cerr<<"no txt/gml\n"; return 1; }} else files.push_back(tgt);

    ofstream fout("tempoutput.txt",ios::trunc);
    if(!fout){ cerr<<"cannot write tempoutput.txt\n"; return 1;}

    for(const string&fn:files){
        Graph G=readEL(fn);
        for(int H=2;H<=4;++H){
            double alpha=0; clock_t t0=clock();
            vector<int> cds=exactDS(H,G.adj,alpha);
            clock_t t1=clock(); double ms=1000.0*(t1-t0)/CLOCKS_PER_SEC;
            block(fout,H,fn,alpha,cds,(long long)G.adj.size(),G.m,ms);
        }
    }
    fout.close(); return 0;
}