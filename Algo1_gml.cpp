
#include <bits/stdc++.h>
#include <windows.h>
using namespace std;

/* ----------------------------------------------------------------
 *  1.  Helper: recursive .gml discovery under <rootDir>
 * ----------------------------------------------------------------*/
static void listGMLFiles(const string &root, vector<string> &out)
{
    WIN32_FIND_DATAA ffd;
    HANDLE h = FindFirstFileA((root + "\\*").c_str(), &ffd);
    if (h == INVALID_HANDLE_VALUE)
        return;

    do
    {
        string name = ffd.cFileName;
        if (name == "." || name == "..")
            continue;
        string path = root + "\\" + name;
        if (ffd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)
            listGMLFiles(path, out); // recurse
        else if (path.size() >= 4 &&
                 _stricmp(path.c_str() + path.size() - 4, ".gml") == 0)
            out.push_back(path);
    } while (FindNextFileA(h, &ffd));
    FindClose(h);
    sort(out.begin(), out.end());
}

/* ----------------------------------------------------------------
 *  2.  Very tolerant GML parser  (keeps undirected edges)
 * ----------------------------------------------------------------*/
static inline string trim(const string &s)
{
    auto st = s.find_first_not_of(" \t\r\n");
    auto en = s.find_last_not_of(" \t\r\n");
    return (st == string::npos) ? string() : s.substr(st, en - st + 1);
}

static bool parseGML(const string &path,
                     vector<pair<int, int>> &edges,
                     int &V)
{
    ifstream fin(path);
    if (!fin)
    {
        cerr << "Cannot open " << path << '\n';
        return false;
    }

    V = 0; // largest id seen + 1
    int curSrc = -1;
    bool inNode = false, inEdge = false;
    string line;
    while (getline(fin, line))
    {
        string s = trim(line);
        if (!inNode && !inEdge)
        {
            if (s == "node")
                inNode = true;
            else if (s == "edge")
                inEdge = true;
            continue;
        }
        if (s == "[")
            continue; // begin block
        if (s == "]")
        {
            inNode = inEdge = false;
            continue;
        }

        istringstream iss(s);
        string key;
        iss >> key;

        if (inNode && key == "id")
        {
            int id;
            if (iss >> id)
                V = max(V, id + 1);
        }
        else if (inEdge)
        {
            if (key == "source")
                iss >> curSrc;
            else if (key == "target")
            {
                int tgt;
                if (iss >> tgt)
                {
                    edges.emplace_back(curSrc, tgt);
                    V = max(V, max(curSrc, tgt) + 1);
                }
            }
        }
    }
    return true;
}

/* ----------------------------------------------------------------
 *  3.  Dinic max-flow (unchanged)
 * ----------------------------------------------------------------*/
struct Dinic
{
    struct Edge
    {
        int to, rev;
        double cap;
    };
    vector<vector<Edge>> G;
    vector<int> lvl, it;

    void init(int n)
    {
        G.assign(n, {});
        lvl.assign(n, -1);
        it.assign(n, 0);
    }
    void add_edge(int u, int v, double c)
    {
        G[u].push_back({v, (int)G[v].size(), c});
        G[v].push_back({u, (int)G[u].size() - 1, 0.0});
    }
    bool bfs(int s, int t)
    {
        fill(lvl.begin(), lvl.end(), -1);
        queue<int> q;
        lvl[s] = 0;
        q.push(s);
        while (!q.empty())
        {
            int u = q.front();
            q.pop();
            for (auto &e : G[u])
                if (e.cap > 1e-9 && lvl[e.to] < 0)
                {
                    lvl[e.to] = lvl[u] + 1;
                    q.push(e.to);
                }
        }
        return lvl[t] >= 0;
    }
    double dfs(int u, int t, double f)
    {
        if (u == t)
            return f;
        for (int &i = it[u]; i < (int)G[u].size(); ++i)
        {
            auto &e = G[u][i];
            if (e.cap > 1e-9 && lvl[e.to] == lvl[u] + 1)
            {
                double got = dfs(e.to, t, min(f, e.cap));
                if (got > 1e-9)
                {
                    e.cap -= got;
                    G[e.to][e.rev].cap += got;
                    return got;
                }
            }
        }
        return 0.0;
    }
    double max_flow(int s, int t)
    {
        double flow = 0, f;
        while (bfs(s, t))
        {
            fill(it.begin(), it.end(), 0);
            while ((f = dfs(s, t, 1e100)) > 1e-9)
                flow += f;
        }
        return flow;
    }
    vector<char> reachable(int s)
    {
        int n = G.size();
        vector<char> vis(n, 0);
        queue<int> q;
        vis[s] = 1;
        q.push(s);
        while (!q.empty())
        {
            int u = q.front();
            q.pop();
            for (auto &e : G[u])
                if (e.cap > 1e-9 && !vis[e.to])
                {
                    vis[e.to] = 1;
                    q.push(e.to);
                }
        }
        return vis;
    }
};

/* ----------------------------------------------------------------
 *  4.  Fast triangle list (H = 3) and naïve clique list (H ≥ 4)
 *      –– exactly the same as in your original file
 * ----------------------------------------------------------------*/
static void listTriangles(const vector<vector<int>> &adj,
                          vector<array<int, 3>> &tri)
{
    int n = adj.size();
    vector<int> ord(n);
    iota(ord.begin(), ord.end(), 0);
    sort(ord.begin(), ord.end(),
         [&](int a, int b)
         {
             if (adj[a].size() != adj[b].size())
                 return adj[a].size() < adj[b].size();
             return a < b;
         });
    vector<int> pos(n);
    for (int i = 0; i < n; ++i)
        pos[ord[i]] = i;
    vector<vector<int>> fwd(n);
    for (int u = 0; u < n; ++u)
        for (int v : adj[u])
            if (pos[u] < pos[v])
                fwd[u].push_back(v);
    vector<char> mark(n, 0);
    for (int u = 0; u < n; ++u)
    {
        for (int v : fwd[u])
            mark[v] = 1;
        for (int v : fwd[u])
            for (int w : fwd[v])
                if (mark[w])
                    tri.push_back({u, v, w});
        for (int v : fwd[u])
            mark[v] = 0;
    }
}
static void enumCliqueRec(const vector<vector<int>> &adj,
                          vector<int> &cand, vector<int> &cur, int k,
                          vector<vector<int>> &out)
{
    if ((int)cur.size() == k)
    {
        out.push_back(cur);
        return;
    }
    while (!cand.empty())
    {
        int v = cand.back();
        cand.pop_back();
        vector<int> nxt;
        for (int u : cand)
            if (binary_search(adj[v].begin(), adj[v].end(), u))
                nxt.push_back(u);
        cur.push_back(v);
        enumCliqueRec(adj, nxt, cur, k, out);
        cur.pop_back();
    }
}
static void listCliques(const vector<vector<int>> &adj, int k,
                        vector<vector<int>> &out)
{
    vector<int> cand(adj.size());
    iota(cand.begin(), cand.end(), 0);
    vector<int> cur;
    enumCliqueRec(adj, cand, cur, k, out);
}

/* ----------------------------------------------------------------
 *  5.  Build flow network and exact densest-H-subgraph
 *      (identical to original “Fixed.cpp”)
 * ----------------------------------------------------------------*/
static void build_flow(int H, const vector<vector<int>> &adj, double a, Dinic &D)
{
    int n = adj.size();
    if (H == 2)
    { // edge-density
        D.init(n + 2);
        for (int u = 0; u < n; ++u)
        {
            D.add_edge(0, 2 + u, adj[u].size());
            D.add_edge(2 + u, 1, 2 * a);
            for (int v : adj[u])
                if (u < v)
                {
                    D.add_edge(2 + u, 2 + v, 1);
                    D.add_edge(2 + v, 2 + u, 1);
                }
        }
        return;
    }
    vector<vector<int>> inst;
    if (H == 3)
    {
        vector<array<int, 3>> tri;
        listTriangles(adj, tri);
        inst.reserve(tri.size());
        for (auto &t : tri)
            inst.push_back({t[0], t[1], t[2]});
    }
    else
    {
        listCliques(adj, H - 1, inst);
    }
    int m = inst.size();
    D.init(2 + n + m);
    for (int u = 0; u < n; ++u)
    {
        D.add_edge(0, 2 + u, adj[u].size());
        D.add_edge(2 + u, 1, H * a);
    }
    for (int i = 0; i < m; ++i)
    {
        int c = 2 + n + i;
        for (int v : inst[i])
        {
            D.add_edge(c, 2 + v, 1e100);
            D.add_edge(2 + v, c, 1);
        }
    }
}
static vector<int> exactDS(int H, const vector<vector<int>> &adj, double &alpha)
{
    int n = adj.size();
    double lo = 0, hi = 0;
    for (auto &v : adj)
        hi = max(hi, (double)v.size());
    vector<int> best;
    auto feas = [&](double a)
    {
        Dinic D;
        build_flow(H, adj, a, D);
        D.max_flow(0, 1);
        auto vis = D.reachable(0);
        vector<int> cur;
        for (int i = 0; i < n; ++i)
            if (vis[2 + i])
                cur.push_back(i);
        if (!cur.empty())
            best = move(cur);
        return !best.empty();
    };
    const double tol = 1e-6;
    while (hi - lo > tol)
    {
        double mid = (lo + hi) / 2;
        if (feas(mid))
            lo = mid;
        else
            hi = mid;
    }
    alpha = lo;
    return best;
}

/* ----------------------------------------------------------------
 *  6.  Write block exactly like “Graph1”
 * ----------------------------------------------------------------*/
static void writeBlock(ostream &f, int H, int N, long long M, double alpha,
                       const vector<int> &cds, long long ms)
{
    f << "---------------------------------------------\n";
    f << "Nodes: " << N << " Edges: " << M << "\n";
    f << "H: " << H << "\n";
    f << fixed << setprecision(6) << "Alpha: " << alpha << "\n";
    f << "Nodes in CDS: " << cds.size() << "\nCDS:";
    for (int v : cds)
        f << " " << v;
    f << "\nTime taken: " << ms << " ms\n";
}

/* ----------------------------------------------------------------
 *  7.  Main
 * ----------------------------------------------------------------*/
int main(int argc, char **argv)
{
    ios::sync_with_stdio(false);

    if (argc != 3 && !(argc == 4 && (string(argv[1]) == "--all" || string(argv[1]) == "-all")))
    {
        cerr << "Usage:\n  " << argv[0] << "  graph.gml  H\n"
             << "  " << argv[0] << "  --all  folder  H\n";
        return 1;
    }

    bool batch = (string(argv[1]) == "--all" || string(argv[1]) == "-all");
    string target = batch ? argv[2] : argv[1];
    int H = batch ? atoi(argv[3]) : atoi(argv[2]);

    vector<string> files;
    if (batch)
    {
        listGMLFiles(target, files);
        if (files.empty())
        {
            cerr << "No .gml files under \"" << target << "\"\n";
            return 1;
        }
    }
    else
        files.push_back(target);

    ofstream fout("results_temp.txt", ios::trunc);
    if (!fout)
    {
        cerr << "Cannot create results_temp.txt\n";
        return 1;
    }

    for (auto &path : files)
    {
        vector<pair<int, int>> edges;
        int N = 0;
        if (!parseGML(path, edges, N))
            continue;

        vector<vector<int>> adj(N);
        long long M = 0;
        for (auto &e : edges)
        {
            int u = e.first, v = e.second;
            if (u >= 0 && v >= 0 && u != v && u < N && v < N)
            {
                adj[u].push_back(v);
                adj[v].push_back(u);
                ++M;
            }
        }
        M /= 2; // undirected

        auto t0 = chrono::high_resolution_clock::now();
        double alpha = 0;
        auto cds = exactDS(H, adj, alpha);
        auto t1 = chrono::high_resolution_clock::now();
        long long ms = chrono::duration_cast<chrono::milliseconds>(t1 - t0).count();

        writeBlock(fout, H, N, M, alpha, cds, ms);
    }
    return 0;
}