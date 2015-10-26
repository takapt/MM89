#ifndef LOCAL
#define NDEBUG
#endif

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <climits>
#include <cfloat>
#include <ctime>
#include <cassert>
#include <map>
#include <utility>
#include <set>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>
#include <sstream>
#include <complex>
#include <stack>
#include <queue>
#include <numeric>
#include <list>
#include <iomanip>
#include <fstream>
#include <bitset>

using namespace std;

#define foreach(it, c) for (__typeof__((c).begin()) it=(c).begin(); it != (c).end(); ++it)
template <typename T> void print_container(ostream& os, const T& c) { const char* _s = " "; if (!c.empty()) { __typeof__(c.begin()) last = --c.end(); foreach (it, c) { os << *it; if (it != last) os << _s; } } }
template <typename T> ostream& operator<<(ostream& os, const vector<T>& c) { print_container(os, c); return os; }
template <typename T> ostream& operator<<(ostream& os, const set<T>& c) { print_container(os, c); return os; }
template <typename T> ostream& operator<<(ostream& os, const multiset<T>& c) { print_container(os, c); return os; }
template <typename T> ostream& operator<<(ostream& os, const deque<T>& c) { print_container(os, c); return os; }
template <typename T, typename U> ostream& operator<<(ostream& os, const map<T, U>& c) { print_container(os, c); return os; }
template <typename T, typename U> ostream& operator<<(ostream& os, const pair<T, U>& p) { os << "(" << p.first << ", " << p.second << ")"; return os; }

template <typename T> void print(T a, int n, const string& split = " ") { for (int i = 0; i < n; i++) { cerr << a[i]; if (i + 1 != n) cerr << split; } cerr << endl; }
template <typename T> void print2d(T a, int w, int h, int width = -1, int br = 0) { for (int i = 0; i < h; ++i) { for (int j = 0; j < w; ++j) { if (width != -1) cerr.width(width); cerr << a[i][j] << ' '; } cerr << endl; } while (br--) cerr << endl; }
template <typename T> void input(T& a, int n) { for (int i = 0; i < n; ++i) cin >> a[i]; }
#define dump(v) (cerr << #v << ": " << v << endl)
// #define dump(v)

#define rep(i, n) for (int i = 0; i < (int)(n); ++i)
#define erep(i, n) for (int i = 0; i <= (int)(n); ++i)
#define all(a) (a).begin(), (a).end()
#define rall(a) (a).rbegin(), (a).rend()
#define clr(a, x) memset(a, x, sizeof(a))
#define sz(a) ((int)(a).size())
#define mp(a, b) make_pair(a, b)
#define ten(n) ((long long)(1e##n))

template <typename T, typename U> void upmin(T& a, const U& b) { a = min<T>(a, b); }
template <typename T, typename U> void upmax(T& a, const U& b) { a = max<T>(a, b); }
template <typename T> void uniq(T& a) { sort(a.begin(), a.end()); a.erase(unique(a.begin(), a.end()), a.end()); }
template <class T> string to_s(const T& a) { ostringstream os; os << a; return os.str(); }
template <class T> T to_T(const string& s) { istringstream is(s); T res; is >> res; return res; }
bool in_rect(int x, int y, int w, int h) { return 0 <= x && x < w && 0 <= y && y < h; }

typedef long long ll;
typedef pair<int, int> pint;
typedef unsigned long long ull;

const int DX[] = { 0, 1, 0, -1 };
const int DY[] = { 1, 0, -1, 0 };




ull rdtsc()
{
#ifdef __amd64
    ull a, d;
    __asm__ volatile ("rdtsc" : "=a" (a), "=d" (d));
    return (d<<32) | a;
#else
    ull x;
    __asm__ volatile ("rdtsc" : "=A" (x));
    return x;
#endif
}
#ifdef LOCAL
const double CYCLES_PER_SEC = 3.30198e9;
#else
const double CYCLES_PER_SEC = 2.5e9;
#endif
double get_absolute_sec()
{
    return (double)rdtsc() / CYCLES_PER_SEC;
}
#ifdef _MSC_VER
#include <Windows.h>
    double get_ms() { return (double)GetTickCount64() / 1000; }
#else
#include <sys/time.h>
    double get_ms() { struct timeval t; gettimeofday(&t, NULL); return (double)t.tv_sec * 1000 + (double)t.tv_usec / 1000; }
#endif

#ifndef LOCAL
#define USE_RDTSC
#endif
class Timer
{
private:
    double start_time;
    double elapsed;

#ifdef USE_RDTSC
    double get_sec() { return get_absolute_sec(); }
#else
    double get_sec() { return get_ms() / 1000; }
#endif

public:
    Timer() {}

    void start() { start_time = get_sec(); }
    double get_elapsed() { return elapsed = get_sec() - start_time; }
};
Timer g_timer;
#ifdef LOCAL
#define USE_TIMER
#ifdef USE_TIMER
const double G_TL_SEC = 9.5;
#else
const double G_TL_SEC = 1e9;
#endif
#else
#define USE_TIMER
const double G_TL_SEC = 9.8;
#endif


struct Pos
{
    int x, y;
    Pos(int x, int y)
        : x(x), y(y)
    {
    }
    Pos()
        : x(0), y(0)
    {
    }

    bool operator==(const Pos& other) const
    {
        return x == other.x && y == other.y;
    }
    bool operator !=(const Pos& other) const
    {
        return x != other.x || y != other.y;
    }

    void operator+=(const Pos& other)
    {
        x += other.x;
        y += other.y;
    }
    void operator-=(const Pos& other)
    {
        x -= other.x;
        y -= other.y;
    }

    Pos operator+(const Pos& other) const
    {
        Pos res = *this;
        res += other;
        return res;
    }
    Pos operator-(const Pos& other) const
    {
        Pos res = *this;
        res -= other;
        return res;
    }
    Pos operator*(int a) const
    {
        return Pos(x * a, y * a);
    }

    bool operator<(const Pos& other) const
    {
        if (x != other.x)
            return x < other.x;
        else
            return y < other.y;
    }

    int sq_dist(const Pos& p) const
    {
        int dx = x - p.x;
        int dy = y - p.y;
        return dx * dx + dy * dy;
    }
    int dist(const Pos& p) const
    {
        return abs(p.x - x) + abs(p.y - y);
    }

    Pos next(int dir) const
    {
        return Pos(x + DX[dir], y + DY[dir]);
    }

    void move(int dir)
    {
        x += DX[dir];
        y += DY[dir];
    }

    int dir(const Pos& to) const
    {
        rep(dir, 4)
        {
            if (next(dir) == to)
                return dir;
        }
        assert(false);
        return -1;
    }

    uint pack() const
    {
        return (y << 7) | x;
    }
};
Pos operator*(int a, const Pos& pos)
{
    return pos * a;
}
ostream& operator<<(ostream& os, const Pos& pos)
{
    os << "(" << pos.x << ", " << pos.y << ")";
    return os;
}

class Random
{
private:
    unsigned int  x, y, z, w;
public:
    Random(unsigned int x
             , unsigned int y
             , unsigned int z
             , unsigned int w)
        : x(x), y(y), z(z), w(w) { }
    Random()
        : x(123456789), y(362436069), z(521288629), w(88675123) { }
    Random(unsigned int seed)
        : x(123456789), y(362436069), z(521288629), w(seed) { }

    unsigned int next()
    {
        unsigned int t = x ^ (x << 11);
        x = y;
        y = z;
        z = w;
        return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
    }

    int next_int() { return next(); }

    // [0, upper)
    int next_int(int upper) { return next() % upper; }

    // [low, high)
    int next_int(int low, int high) { return next_int(high - low) + low; }

    double next_double(double upper) { return upper * next() / UINT_MAX; }
    double next_double(double low, double high) { return next_double(high - low) + low; }

    template <typename T>
    int select(const vector<T>& ratio)
    {
        T sum = accumulate(ratio.begin(), ratio.end(), (T)0);
        double v = next_double(sum) + 1e-6;
        for (int i = 0; i < (int)ratio.size(); ++i)
        {
            v -= ratio[i];
            if (v <= 0)
                return i;
        }
        return (int)ratio.size() - 1;
    }
};
Random g_rand;

const char* CELL_STR = "SLURE.";
const int STRAIGHT = 0;
const int LEFT = 1;
const int UTURN = 2;
const int RIGHT = 3;
const int EVERY = 4;
const int OUT = 5;
int to_dir_diff(char c)
{
    int d = strchr(CELL_STR, c) - CELL_STR;
    assert(0 <= d && d < 6);
    return d;
}
char to_char(int dir_diff)
{
    assert(0 <= dir_diff && dir_diff < 6);
    return CELL_STR[dir_diff];
}

class BoolBoard
{
public:
    BoolBoard(int w, int h) :
        w(w), h(h), f{}
    {
    }

    BoolBoard(){}

    bool get(int x, int y) const
    {
        assert(in_rect(x, y, w, h));
        return f[y][x];
    }

    void set(int x, int y, bool v)
    {
        assert(in_rect(x, y, w, h));
        f[y][x] = v;
    }

    int count() const
    {
        int c = 0;
//         rep(y, h)
//             c += f[y].count();
        rep(y, h) rep(x, w)
            c += f[y][x];
        return c;
    }

private:
    int w, h;
//     bitset<80> f[80];
    bool f[80][80];
};

int calls;
int waste_calls;
class Maze
{
public:
    Maze(const vector<string>& maze) :
        w(maze[0].size()), h(maze.size())
    {
        rep(y, h) rep(x, w)
            set(x, y, to_dir_diff(maze[y][x]));
    }
    Maze(){}

    int width() const { return w; }
    int height() const { return h; }

    int get(int x, int y) const
    {
        assert(in_rect(x, y, w, h));
        return a[y][x];
    }

    void set(int x, int y, int d)
    {
        assert(in_rect(x, y, w, h));
        assert(0 <= d && d < 6);
        assert(!outside(x, y));
        a[y][x] = d;
    }

    bool outside(int x, int y) const
    {
        assert(in_rect(x, y, w, h));
        return get(x, y) == OUT;
    }

    int count_changes(const Maze& ori_maze) const
    {
        int c = 0;
        rep(y, h) rep(x, w)
            c += get(x, y) != ori_maze.get(x, y);
        return c;
    }

    int count_covers() const
    {
        BoolBoard f = search_covered(list_borders());
        int c = 0;
        rep(y, h) rep(x, w)
            c += f.get(x, y);
        return c;
    }

    bool border(int x, int y) const
    {
        if (!outside(x, y))
            return false;

        rep(dir, 4)
        {
            int nx = x + DX[dir], ny = y + DY[dir];
            if (in_rect(nx, ny, w, h) && !outside(nx, ny))
                return true;
        }
        return false;
    }

    vector<Pos> list_borders() const
    {
        vector<Pos> borders;
        rep(y, h) rep(x, w)
            if (border(x, y))
                borders.push_back(Pos(x, y));
        return borders;
    }

    BoolBoard search_covered(const vector<Pos>& borders) const
    {
//         calls = 0;
//         waste_calls = 0;
        BoolBoard visited(w, h);
        BoolBoard covered(w, h);
        for (auto& start : borders)
        {
            assert(border(start.x, start.y));
            rep(dir, 4)
            {
                int nx = start.x + DX[dir], ny = start.y + DY[dir];
                if (in_rect(nx, ny, w, h) && !outside(nx, ny))
                    search(nx, ny, dir, visited, covered);
            }
        }
        return covered;
    }

    vector<vector<Pos>> list_paths() const
    {
        auto borders = list_borders();
        BoolBoard visited(w, h);
        vector<vector<Pos>> paths;
        vector<Pos> path;
        for (auto& start : borders)
        {
            assert(border(start.x, start.y));
            path.push_back(start);
            rep(dir, 4)
            {
                int nx = start.x + DX[dir], ny = start.y + DY[dir];
                if (in_rect(nx, ny, w, h) && !outside(nx, ny))
                    search_path(nx, ny, dir, visited, path, paths);
            }
            path.pop_back();
        }
        return paths;
    }

private:
    bool search(int x, int y, int dir, BoolBoard& visited, BoolBoard& covered) const
    {
        assert(in_rect(x, y, w, h));
        assert(!visited.get(x, y));

        if (outside(x, y))
            return true;

        visited.set(x, y, true);

        bool any_exit = false;
        if (get(x, y) == EVERY)
        {
            rep(ndir, 4)
            {
                int nx = x + DX[ndir], ny = y + DY[ndir];
                assert(in_rect(nx, ny, w, h));
                if (!visited.get(nx, ny))
                    any_exit |= search(nx, ny, ndir, visited, covered);
            }
        }
        else
        {
            assert(0 <= get(x, y) && get(x, y) < 4);
            int ndir = (dir + get(x, y)) % 4;
            int nx = x + DX[ndir], ny = y + DY[ndir];
            assert(in_rect(nx, ny, w, h));
            if (!visited.get(nx, ny))
                any_exit |= search(nx, ny, ndir, visited, covered);
        }

        visited.set(x, y, false);

        if (any_exit)
            covered.set(x, y, true);

//         ++calls;
//         if (!any_exit)
//             ++waste_calls;

        return any_exit;
    }

    bool search_path(int x, int y, int dir, BoolBoard& visited, vector<Pos>& path, vector<vector<Pos>>& paths) const
    {
        assert(in_rect(x, y, w, h));
        assert(!visited.get(x, y));

        if (outside(x, y))
        {
            assert(path.size() >= 2);
            paths.push_back(path);
            paths.back().push_back(Pos(x, y));
            return true;
        }

        visited.set(x, y, true);
        path.push_back(Pos(x, y));

        bool any_exit = false;
        if (get(x, y) == EVERY)
        {
            rep(ndir, 4)
            {
                int nx = x + DX[ndir], ny = y + DY[ndir];
                assert(in_rect(nx, ny, w, h));
                if (!visited.get(nx, ny))
                    any_exit |= search_path(nx, ny, ndir, visited, path, paths);
            }
        }
        else
        {
            assert(0 <= get(x, y) && get(x, y) < 4);
            int ndir = (dir + get(x, y)) % 4;
            int nx = x + DX[ndir], ny = y + DY[ndir];
            assert(in_rect(nx, ny, w, h));
            if (!visited.get(nx, ny))
                any_exit |= search_path(nx, ny, ndir, visited, path, paths);
        }

        visited.set(x, y, false);
        path.pop_back();

        return any_exit;
    }

    int w, h;
    int a[80][80];
};


#ifdef LOCAL
int test_case_seed;

// #define VIS
#ifdef VIS
#include <sys/stat.h>
#include "cairo.h"

string get_dir_path(string filepath)
{
    assert(filepath.size() > 0);
    while (filepath[filepath.size() - 1] != '/')
        filepath.erase(filepath.begin() + (filepath.size() - 1));
    return filepath;
}

string make_filename(int index)
{
    char filename[256];
    sprintf(filename, "image/%d/%d.png", test_case_seed, index);
    return filename;
}

struct DrawPaths
{
    vector<vector<Pos>> paths;
    double r, g, b;
    DrawPaths(const vector<vector<Pos>>& paths, double r, double g, double b) :
        paths(paths), r(r), g(g), b(b)
    {
    }
};
void save_image(const string& filename, const Maze& maze, const Maze& prev_maze, const Maze& start_maze, const vector<DrawPaths>& draw_paths)
{
    const auto borders = maze.list_borders();
    const auto covered = maze.search_covered(borders);
    const auto prev_covered = prev_maze.search_covered(borders);

    const int CELL_SIZE = 30;
    const int OFFSET = 5;
    const int w = maze.width(), h = maze.height();
    const int surface_w = w * CELL_SIZE + 2 * OFFSET;
    const int surface_h = h * CELL_SIZE + 2 * OFFSET;

    cairo_surface_t* surface = cairo_image_surface_create(CAIRO_FORMAT_RGB24, surface_w, surface_h);
    cairo_t* cr = cairo_create(surface);

    cairo_set_source_rgb(cr, 1, 1, 1);
    cairo_rectangle(cr, 0, 0, surface_w, surface_h);
    cairo_fill(cr);

    rep(y, h) rep(x, w)
    {
        if (!maze.outside(x, y))
        {
            double r = 1, g = 1, b = 1;
            if (covered.get(x, y))
                r = g = b = 0.8;
            if (covered.get(x, y) && !prev_covered.get(x, y))
            {
                r = 1;
                g = b = 0.4;
            }
            else if (!covered.get(x, y) && prev_covered.get(x, y))
            {
                r = 0.2;
                g = b = 0.8;
            }

            if (maze.get(x, y) != prev_maze.get(x, y))
            {
                r = 0.8;
                g = b = 0.2;
            }
//             if (maze.get(x, y) != start_maze.get(x, y))
//             {
//                 g -= 0.25;
//                 b -= 0.25;
//             }
            cairo_set_source_rgb(cr, r, g, b);
            cairo_rectangle(cr, OFFSET + x * CELL_SIZE, OFFSET + y * CELL_SIZE, CELL_SIZE, CELL_SIZE);
            cairo_fill(cr);
        }
    }

    rep(y, h + 1)
    {
        cairo_move_to(cr, OFFSET + 0, OFFSET + y * CELL_SIZE);
        cairo_line_to(cr, OFFSET + w * CELL_SIZE, OFFSET + y * CELL_SIZE);
    }
    rep(x, w + 1)
    {
        cairo_move_to(cr, OFFSET + x * CELL_SIZE, OFFSET + 0);
        cairo_line_to(cr, OFFSET + x * CELL_SIZE, OFFSET + h * CELL_SIZE);
    }
    cairo_set_source_rgb(cr, 0, 0, 0);
    cairo_set_line_width(cr, 2);
    cairo_stroke(cr);

    rep(y, h) for (int x = 1; x < w; ++x)
    {
        if (maze.outside(x - 1, y) != maze.outside(x, y))
        {
            cairo_move_to(cr, OFFSET + x * CELL_SIZE, OFFSET + y * CELL_SIZE);
            cairo_line_to(cr, OFFSET + x * CELL_SIZE, OFFSET + (y + 1) * CELL_SIZE);
        }
    }
    rep(x, w) for (int y = 1; y < h; ++y)
    {
        if (maze.outside(x, y - 1) != maze.outside(x, y))
        {
            cairo_move_to(cr, OFFSET + x * CELL_SIZE, OFFSET + y * CELL_SIZE);
            cairo_line_to(cr, OFFSET + (x + 1) * CELL_SIZE, OFFSET + y * CELL_SIZE);
        }
    }
    cairo_set_line_width(cr, 4);
    cairo_stroke(cr);

    for (auto& draw_path : draw_paths)
    {
        cairo_set_source_rgb(cr, draw_path.r, draw_path.g, draw_path.b);
        for (auto& path : draw_path.paths)
        {
            assert(path.size() >= 3);

            cairo_rectangle(cr, OFFSET + path[0].x * CELL_SIZE + CELL_SIZE * 0.4, OFFSET + path[0].y * CELL_SIZE + CELL_SIZE * 0.4, 0.2 * CELL_SIZE, 0.2 * CELL_SIZE);
            rep(i, (int)path.size() - 1)
            {
                int x = path[i].x, y = path[i].y;
                int nx = path[i + 1].x, ny = path[i + 1].y;
                cairo_move_to(cr, OFFSET + x * CELL_SIZE + 0.5 * CELL_SIZE, OFFSET + y * CELL_SIZE + 0.5 * CELL_SIZE);
                cairo_line_to(cr, OFFSET + nx * CELL_SIZE + 0.5 * CELL_SIZE, OFFSET + ny * CELL_SIZE + 0.5 * CELL_SIZE);
            }
            cairo_set_line_width(cr, 4);
            cairo_stroke(cr);
        }
    }

    cairo_set_source_rgb(cr, 0, 0, 0);
    cairo_set_font_size(cr, CELL_SIZE * 0.7);
    rep(y, h) rep(x, w)
    {
        if (!maze.outside(x, y))
        {
            cairo_move_to(cr, OFFSET + x * CELL_SIZE + 0.25 * CELL_SIZE, OFFSET + (y + 1) * CELL_SIZE - 0.2 * CELL_SIZE);
            char buf[] = {to_char(maze.get(x, y)), '\0'};
            cairo_show_text(cr, buf);
        }
    }


    mkdir(get_dir_path(filename).c_str(), 0777);
    cairo_surface_write_to_png(surface, filename.c_str());

    cairo_destroy(cr);
    cairo_surface_destroy(surface);
}
#endif
#endif

double sum = 0;
Timer prof_timer;
class Solver
{
public:
    Solver(const Maze& start_maze, const int max_changes) :
        start_maze(start_maze), max_changes(max_changes), borders(start_maze.list_borders()),
        w(start_maze.width()), h(start_maze.height())
    {
        rep(y, h) rep(x, w)
            dist_to_outside[y][x] = -1;

        queue<Pos> q;
        for (auto& p : borders)
        {
            dist_to_outside[p.y][p.x] = 0;
            q.push(p);
        }
        while (!q.empty())
        {
            Pos p = q.front();
            q.pop();

            rep(dir, 4)
            {
                Pos next = p.next(dir);
                if (in_rect(next.x, next.y, w, h) && !start_maze.outside(next.x, next.y) && dist_to_outside[next.y][next.x] == -1)
                {
                    dist_to_outside[next.y][next.x] = dist_to_outside[p.y][p.x] + 1;
                    q.push(next);
                }
            }
        }
    }

    struct State
    {
        int dir;

        int new_covers;
        int orig_to_change;
        int change_to_change;
        int change_to_orig;
        int everys;
        int cover_breaks;

        Pos pos;
        State* prev;

        State(int dir, int new_covers, int orig_to_change, int change_to_change, int change_to_orig, int everys, const Pos& pos, State* prev) :
            dir(dir),
            new_covers(new_covers),
            orig_to_change(orig_to_change),
            change_to_change(change_to_change),
            change_to_orig(change_to_orig),
            everys(everys),
            cover_breaks(0),
            pos(pos),
            prev(prev)
        {
        }
        State(){}

        State next_state(const Pos& next, int dir)
        {
            State s = *this;
            s.pos = next;
            s.dir = dir;
            s.prev = this;
            return s;
        }

        int score;
        int eval() const
        {
            return score;
        }

        vector<Pos> make_rev_path() const
        {
            vector<Pos> rev_path;
            for (const State* s = this; s != nullptr; s = s->prev)
                rev_path.push_back(s->pos);
            return rev_path;
        }
        void make_rev_path(vector<Pos>& rev_path) const
        {
            rev_path.clear();
            for (const State* s = this; s != nullptr; s = s->prev)
                rev_path.push_back(s->pos);
        }

        vector<Pos> make_path() const
        {
            auto path = make_rev_path();
            reverse(all(path));
            return path;
        }
        void make_path(vector<Pos>& path) const
        {
            path.clear();
            make_rev_path(path);
            reverse(all(path));
        }

        bool on_path(const Pos& p) const
        {
            for (const State* s = this; s->prev != nullptr; s = s->prev)
                if (p == s->pos)
                    return true;
            return false;
        }

        bool operator<(const State& other) const
        {
            return eval() > other.eval();
        }
    };
    struct SearchNewPathResult
    {
        Maze maze;
        State state;
        int changes;
        vector<Pos> path;

#ifdef VIS
        vector<DrawPaths> draw_paths;
#endif
    };
    SearchNewPathResult search_new_path(const Maze& init_maze, const double progress) const
    {
        const auto init_covered = init_maze.search_covered(borders);
        const int init_covers = init_covered.count();
        int init_changes = 0;
        int inside = 0;
        rep(y, h) rep(x, w)
        {
            if (start_maze.get(x, y) != init_maze.get(x, y))
                ++init_changes;
            if (!start_maze.outside(x, y))
                ++inside;
        }

        const int MAX_COEF_DEPTH = 3;

        const bool USE_BREAK = max_changes >= 810;

        const int BEAM_WIDTH = USE_BREAK ? 10 : 50;
        const double COEF_DEPTH = MAX_COEF_DEPTH * (1 - progress);
        const int MAX_DEPTH = COEF_DEPTH * (init_maze.width() + init_maze.height()) / 2;
        static vector<State> states[MAX_COEF_DEPTH * 80 + 10];
        rep(i, MAX_DEPTH)
            states[i].clear();

        const Pos start_pos = borders[g_rand.next_int(borders.size())];
        State start_state;
        start_state.pos = start_pos;
        start_state.prev = nullptr;
        rep(dir, 4)
        {
            int nx = start_pos.x + DX[dir];
            int ny = start_pos.y + DY[dir];
            if (in_rect(nx, ny, w, h) && !init_maze.outside(nx, ny))
            {
                states[0].push_back(State(dir, !init_covered.get(nx, ny), 0, 0, 0, 0, Pos(nx, ny), &start_state));
            }
        }

        auto eval = [&](const State& state)
        {
            if ((double)max_changes / inside < 0.3)
            {
                int score = 0;
                score += 100 * state.new_covers;
                score -= 200 * state.orig_to_change;

                if (USE_BREAK)
                    score -= 200 * state.cover_breaks;// * (1 - progress);
                else
                    score -= 200 * state.change_to_change;
                return score;
            }
            else
            {
                int score = 0;
                score += 100 * state.new_covers;
                score -= 100 * state.orig_to_change;

                //             if (USE_BREAK)
                score -= 200 * state.cover_breaks;// * (1 - progress);
                //             else
                //                 score -= 200 * state.change_to_change;
                return score;
            }
        };
//         auto eval = [&](const State& state)
//         {
//             int score = 0;
//             score += 100 * state.new_covers;
//             score -= 100 * state.orig_to_change;
//
// //             if (USE_BREAK)
//                 score -= 200 * state.cover_breaks;// * (1 - progress);
// //             else
// //                 score -= 200 * state.change_to_change;
//             return score;
//         };

        static vector<State> complete_states;
        complete_states.clear();

        rep(depth, MAX_DEPTH - 1)
        {
            if (g_timer.get_elapsed() > G_TL_SEC)
                break;

            sort(all(states[depth]));
            if (states[depth].size() > BEAM_WIDTH)
                states[depth].erase(states[depth].begin() + BEAM_WIDTH, states[depth].end());

            for (auto& state : states[depth])
            {
                const int cur_x = state.pos.x;
                const int cur_y = state.pos.y;

                static BoolBoard on_path(80, 80);
                static vector<Pos> rev_path;
                state.make_rev_path(rev_path);
                for (auto& p : rev_path)
                    on_path.set(p.x, p.y, true);

                static vector<int> dirs;
                dirs.clear();
                rep(ndir, 4)
                {
                    const int nx = cur_x + DX[ndir];
                    const int ny = cur_y + DY[ndir];
                    if (((state.dir - ndir + 4) % 4 != 2 && dist_to_outside[ny][nx] <= MAX_DEPTH - depth - 1) &&
                            (init_maze.outside(nx, ny) || !on_path.get(nx, ny)))
                    {
                        dirs.push_back(ndir);
                    }
                }
                for (int ndir : dirs)
                {
                    const int nx = cur_x + DX[ndir];
                    const int ny = cur_y + DY[ndir];
                    const Pos npos(nx, ny);
                    assert(in_rect(nx, ny, w, h));

                    if (init_maze.outside(nx, ny))
                    {
                        if (state.new_covers >= 1)
                            complete_states.push_back(state.next_state(npos, ndir));
                        continue;
                    }

                    {
                        int cell_to = (ndir - state.dir + 4) % 4;
                        bool need_change = init_maze.get(cur_x, cur_y) != EVERY && cell_to != init_maze.get(cur_x, cur_y);

                        State nstate = state.next_state(npos, ndir);

                        if (!init_covered.get(nx, ny))
                            ++nstate.new_covers;

                        if (init_maze.get(cur_x, cur_y) == EVERY)
                            ++nstate.everys;

                        if (need_change)
                        {
                            if (cell_to == start_maze.get(cur_x, cur_y))
                                ++nstate.change_to_orig;
                            else if (init_maze.get(cur_x, cur_y) == start_maze.get(cur_x, cur_y))
                                ++nstate.orig_to_change;
                            else
                                ++nstate.change_to_change;

                            if (init_covered.get(cur_x, cur_y))
                                ++nstate.cover_breaks;
                        }

                        nstate.score = eval(nstate);

                        states[depth + 1].push_back(nstate);
                    }
                }
                for (auto& p : rev_path)
                    on_path.set(p.x, p.y, false);
            }
        }

        SearchNewPathResult best_result;
        best_result.maze = init_maze;
        best_result.changes = w * h;
        int best_covers = init_covers;
        for (auto& state : complete_states)
        {
            if (g_timer.get_elapsed() > G_TL_SEC)
                break;

            static vector<Pos> path;
            state.make_path(path);
            assert(path.size() >= 2);

            Maze maze = init_maze;
            for (int i = 1; i < (int)path.size() - 1; ++i)
            {
                int cur_dir = path[i - 1].dir(path[i]);
                int need_cell = (path[i].dir(path[i + 1]) - cur_dir + 4) % 4;
                int x = path[i].x, y = path[i].y;
                if (maze.get(x, y) != EVERY)
                    maze.set(x, y, need_cell);
            }

            auto covered = maze.search_covered(borders);
            int covers = covered.count();
            int changes = 0;
            rep(y, h) rep(x, w)
            {
                if (!maze.outside(x, y))
                {
                    if (!covered.get(x, y))
                        maze.set(x, y, start_maze.get(x, y));

                    changes += start_maze.get(x, y) != maze.get(x, y);
                    if (changes > max_changes)
                        goto END;
                }
            }
END:

//             if (best_covers > covers && changes > max_changes)
//             {
//                 fprintf(stderr, "%4d -> %4d, (%4d %4d) %4d %4d %4d %4d\n", best_covers, covers, changes, max_changes, state.new_covers, state.change_to_change, state.orig_to_change, state.change_to_orig);
//             }

            if (changes <= max_changes)
            {
                auto sc = [&](int covers, int changes)
                {
                    return covers - changes * changes * (1 - progress) / max_changes;
                };
//                 if (covers > best_covers)
                if (sc(covers, changes) > sc(best_covers, best_result.changes))
                {
                    best_covers = covers;
                    best_result.maze = maze;
                    best_result.state = state;
                    best_result.changes = changes;

//                     fprintf(stderr, "%4d %4d %4d %4d\n", state.new_covers, state.change_to_change, state.orig_to_change, state.change_to_orig);
                }
            }
        }
#ifdef VIS
        if (best_covers > init_covers)
        {
            //             fprintf(stderr, "%4d -> %4d, (%4d %4d) (%4d %4d) %4d %4d %4d %4d\n", init_covers, best_covers, best_result.changes, max_changes, (int)path.size(), MAX_DEPTH, state.new_covers, state.change_to_change, state.orig_to_change, state.change_to_orig);

            auto& state = best_result.state;
            auto path = state.make_path();

            auto init_paths = init_maze.list_paths();
            auto best_paths = best_result.maze.list_paths();

            set<vector<Pos>> init_path_set(all(init_paths)), best_path_set(all(best_paths));

            vector<vector<Pos>> common_paths, remove_paths, new_paths;
            for (auto& path : init_paths)
            {
                if (best_path_set.count(path))
                    common_paths.push_back(path);
                else
                    remove_paths.push_back(path);
            }
            for (auto& path : best_paths)
                if (!init_path_set.count(path))
                    new_paths.push_back(path);

            best_result.draw_paths = {
                DrawPaths(common_paths, 0, 0, 0),
                DrawPaths(remove_paths, 0, 1, 0),
                DrawPaths(new_paths, 0, 1, 1),
                DrawPaths({path}, 1, 0, 0),
            };
        }
#endif

        return best_result;
    }


    struct EndState
    {
        int new_covers;
        int orig_to_change;
        int change_to_change;
        int change_to_orig;
        int everys;

        deque<Pos> path;

        template <typename T>
        EndState(const T& other) :
            new_covers(other.new_covers),
            orig_to_change(other.orig_to_change),
            change_to_change(other.change_to_change),
            change_to_orig(other.change_to_orig),
            everys(other.everys)
        {
        }

        EndState(int new_covers, int orig_to_change, int change_to_change, int change_to_orig, int everys, const deque<Pos>&& path) :
            new_covers(new_covers),
            orig_to_change(orig_to_change),
            change_to_change(change_to_change),
            change_to_orig(change_to_orig),
            everys(everys),
            path(path)
        {
        }

        EndState(){}

        int eval() const
        {
            int score = 0;
            score += new_covers;
            //                 score += change_to_orig;
            score -= 2 * orig_to_change;
            score -= 2 * change_to_change;
            return score;
        }

        bool operator<(const EndState& other) const
        {
            return eval() > other.eval();
        }
    };
    struct EndSearchNewPathResult
    {
        Maze maze;
        EndState state;
    };
    EndSearchNewPathResult search_new_path_end(const Maze& init_maze) const
    {
        const auto init_covered = init_maze.search_covered(borders);
        const int init_covers = init_covered.count();
        int init_changes = 0;
        rep(y, h) rep(x, w)
        {
            if (start_maze.get(x, y) != init_maze.get(x, y))
                ++init_changes;
        }


        const int MAX_COEF_DEPTH = 4;

        const int BEAM_WIDTH = 50;
        const int COEF_DEPTH = g_timer.get_elapsed() < G_TL_SEC * 0.5 ? 4 : 2;
        const int MAX_DEPTH = COEF_DEPTH * max(init_maze.width(), init_maze.height());
        static vector<EndState> states[MAX_COEF_DEPTH * 80 + 10];
        rep(i, MAX_DEPTH)
            states[i].clear();

        vector<Pos> not_covered;
        rep(y, h) rep(x, w)
            if (!init_maze.outside(x, y) && !init_covered.get(x, y))
                not_covered.push_back(Pos(x, y));
        const Pos start_pos = not_covered[g_rand.next_int(not_covered.size())];
        rep(dir, 4)
        {
            int x = start_pos.x;
            int y = start_pos.y;
            int nx = start_pos.x + DX[dir];
            int ny = start_pos.y + DY[dir];
            assert(in_rect(nx, ny, w, h));

            int new_covers = (!init_maze.outside(x, y) && !init_covered.get(x, y))
               + (!init_maze.outside(nx, ny) && !init_covered.get(nx, ny));
            states[0].push_back(EndState(new_covers, 0, 0, 0, 0, {start_pos, Pos(nx, ny)}));
        }

        if (false)
        {
        const Pos start_pos = borders[g_rand.next_int(borders.size())];
        rep(dir, 4)
        {
            int nx = start_pos.x + DX[dir];
            int ny = start_pos.y + DY[dir];
            if (in_rect(nx, ny, w, h) && !init_maze.outside(nx, ny))
            {
//                 states[0].push_back(EndState(!init_covered.get(nx, ny), 0, 0, 0, 0, {start_pos, Pos(nx, ny)}));
                states[0].push_back(EndState(!init_covered.get(nx, ny), 0, 0, 0, 0, {Pos(nx, ny), start_pos}));
            }
        }
        }

        static vector<EndState> complete_states;
        complete_states.clear();

        rep(depth, MAX_DEPTH - 1)
        {
            if (g_timer.get_elapsed() > G_TL_SEC)
                break;

            sort(all(states[depth]));
            if (states[depth].size() > BEAM_WIDTH)
                states[depth].erase(states[depth].begin() + BEAM_WIDTH, states[depth].end());

            for (auto& state : states[depth])
            {
                const auto& path = state.path;

                struct Move
                {
                    Move(const Pos& pos, const EndState& state) :
                        pos(pos),
                        new_covers(state.new_covers),
                        orig_to_change(state.orig_to_change),
                        change_to_change(state.change_to_change),
                        change_to_orig(state.change_to_orig),
                        everys(state.everys)
                    {
                    }

                    Pos pos;

                    int new_covers;
                    int orig_to_change;
                    int change_to_change;
                    int change_to_orig;
                    int everys;

                    int eval() const
                    {
                        int score = 0;
                        score += new_covers;
                        //                 score += change_to_orig;
                        score -= 2 * orig_to_change;
                        score -= 2 * change_to_change;
                        return score;
                    }
                };


                auto generate_moves = [&](int cur_x, int cur_y, int cur_dir, int rev_sign, int fix_dir)
                {
                    assert(in_rect(cur_x, cur_y, w, h));

                    vector<int> dirs;
                    rep(ndir, 4)
                    {
                        const int nx = cur_x + DX[ndir];
                        const int ny = cur_y + DY[ndir];
                        if ((ndir + fix_dir + 2 + 4) % 4 != cur_dir &&
                            (init_maze.outside(nx, ny) || find(path.begin(), path.end(), Pos(nx, ny)) == path.end()))
                        {
                            dirs.push_back(ndir);
                        }
                    }

                    vector<Move> moves;
                    for (int ndir : dirs)
                    {
                        const int nx = cur_x + DX[ndir];
                        const int ny = cur_y + DY[ndir];
                        assert(in_rect(nx, ny, w, h));

                        if (init_maze.outside(nx, ny))
                        {
                            if (state.new_covers >= 1)
                                moves.push_back(Move(Pos(nx, ny), state));
                            continue;
                        }

                        {
                            int cell_to = (rev_sign * (ndir - cur_dir) + fix_dir + 4) % 4;
                            bool need_change = init_maze.get(cur_x, cur_y) != EVERY && cell_to != init_maze.get(cur_x, cur_y);

                            Move move(Pos(nx, ny), state);

                            if (!init_covered.get(nx, ny))
                                ++move.new_covers;

                            if (init_maze.get(cur_x, cur_y) == EVERY)
                                ++move.everys;

                            if (need_change)
                            {
                                if (cell_to == start_maze.get(cur_x, cur_y))
                                    ++move.change_to_orig;
                                else if (init_maze.get(cur_x, cur_y) == start_maze.get(cur_x, cur_y))
                                    ++move.orig_to_change;
                                else
                                    ++move.change_to_change;
                            }

                            moves.push_back(move);
                        }
                    }
                    return moves;
                };
                vector<Move> front_moves;
                if (!start_maze.outside(path.back().x, path.back().y))
                    front_moves = generate_moves(path.back().x, path.back().y, (path.end() - 2)->dir(path.back()), +1, 0);

                vector<Move> back_moves;
                if (!start_maze.outside(path.front().x, path.front().y))
                    back_moves = generate_moves(path.front().x, path.front().y, path.front().dir(path[1]), -1, 2);

                auto eval_moves = [](const vector<Move>& moves)
                {
                    int score = -1;
                    for (auto& m : moves)
                        upmax(score, m.eval());
                    return score;
                };

                if (eval_moves(front_moves) >= eval_moves(back_moves))
                {
                    for (auto& move : front_moves)
                    {
                        EndState nstate(move);
                        nstate.path = state.path;
                        nstate.path.push_back(move.pos);
                        if (start_maze.outside(nstate.path.front().x, nstate.path.front().y) && start_maze.outside(nstate.path.back().x, nstate.path.back().y))
                            complete_states.push_back(nstate);
                        else
                            states[depth + 1].push_back(nstate);
                    }
                }
                else
                {
                    for (auto& move : back_moves)
                    {
                        EndState nstate(move);
                        nstate.path = state.path;
                        nstate.path.push_front(move.pos);
                        if (start_maze.outside(nstate.path.front().x, nstate.path.front().y) && start_maze.outside(nstate.path.back().x, nstate.path.back().y))
                            complete_states.push_back(nstate);
                        else
                            states[depth + 1].push_back(nstate);
                    }
                }
            }
        }

        EndSearchNewPathResult best_result;
        best_result.maze = init_maze;
        int best_covers = init_covers;
        for (auto& state : complete_states)
        {
            if (g_timer.get_elapsed() > G_TL_SEC)
                break;

            auto& path = state.path;
            assert(path.size() >= 2);

            Maze maze = init_maze;
            for (int i = 1; i < (int)path.size() - 1; ++i)
            {
                int cur_dir = path[i - 1].dir(path[i]);
                int need_cell = (path[i].dir(path[i + 1]) - cur_dir + 4) % 4;
                int x = path[i].x, y = path[i].y;
                if (maze.get(x, y) != EVERY)
                    maze.set(x, y, need_cell);
            }

            auto covered = maze.search_covered(borders);
            int covers = covered.count();
            int changes = 0;
            rep(y, h) rep(x, w)
            {
                if (!maze.outside(x, y) && !covered.get(x, y))
                    maze.set(x, y, start_maze.get(x, y));

                if (start_maze.get(x, y) != maze.get(x, y))
                    ++changes;

                if (changes > max_changes)
                    goto END;
            }
END:

            if (changes <= max_changes)
            {
                double score = state.eval();
                if (covers > best_covers)
                {
                    best_covers = covers;
                    best_result.maze = maze;
                    best_result.state = state;
                }
            }
        }

        return best_result;
    }



    Maze solve()
    {
#ifdef VIS
        save_image(make_filename(-1), start_maze, start_maze, start_maze, {DrawPaths({}, 0, 0, 0)});
#endif

        Maze best_maze;
        int best_covers = -1;

        const int MAX_MAZES = min(8, 1 * 80 * 80 / (w * h));
        vector<Maze> current_mazes(MAX_MAZES, start_maze);


#ifndef USE_TIMER
        const int MAX_TRIES = 568;
#endif

        double last_halving_progress = 0;
        int try_i;
        for (try_i = 0;
#ifdef USE_TIMER
                g_timer.get_elapsed() < G_TL_SEC;
#else
                try_i < MAX_TRIES;
#endif
                ++try_i)
        {
#ifdef USE_TIMER
            double progress = g_timer.get_elapsed() / G_TL_SEC;
#else
            double progress = (double)try_i / MAX_TRIES;
#endif

            if (current_mazes.size() > 1 && progress - last_halving_progress > 0.2)
            {
                vector<pair<int, int>> order;
                rep(i, current_mazes.size())
                    order.push_back(make_pair(current_mazes[i].search_covered(borders).count(), i));
                sort(rall(order));
//                 dump(order);

                const int take = current_mazes.size() / 2;
                vector<Maze> next_mazes(take);
                rep(i, take)
                    next_mazes[i] = current_mazes[order[i].second];
                current_mazes = std::move(next_mazes);

                last_halving_progress = progress;

//                 dump(current_mazes.size());
            }

            for (auto& current_maze : current_mazes)
            {
                SearchNewPathResult result = search_new_path(current_maze, progress);
                auto& next_maze = result.maze;
                auto next_coverd = next_maze.search_covered(borders);
                int next_covers = next_coverd.count();
                auto sc = [&](int covers, int changes)
                {
                    return covers - changes * changes * (1 - progress) / max_changes;
                };
//                 if (next_covers > current_maze.search_covered(borders).count())
                if (sc(next_covers, next_maze.count_changes(start_maze)) > sc(current_maze.count_covers(), current_maze.count_changes(start_maze)))
                {
//                     fprintf(stderr, "%6d %6d\n", calls, waste_calls);
                    //                 fprintf(stderr, "%4d (%4.2f): %4d -> %4d\n", try_i, g_timer.get_elapsed(), current_covers, next_covers);
#ifdef VIS
                    save_image(make_filename(try_i), next_maze, current_maze, start_maze, result.draw_paths);
#endif
                    current_maze = next_maze;
                }
                if (next_covers > best_covers)
                {
                    best_covers = next_covers;
                    best_maze = next_maze;
                }
            }
        }
        dump(try_i);

        best_maze = use_up(best_maze);
        int covers = best_maze.count_covers();
        for (auto& maze : current_mazes)
        {
            maze = use_up(maze);
            int covers = maze.search_covered(borders).count();
            if (covers > best_covers)
            {
                best_covers = covers;
                best_maze = maze;
            }
        }
        return best_maze;
    }

    Maze use_up(Maze maze) const
    {
        auto covered = maze.search_covered(borders);

        int rem_changes = max_changes;
        rep(y, h) rep(x, w)
            if (maze.get(x, y) != start_maze.get(x, y))
                --rem_changes;

        for (auto& p : borders)
        {
            rep(dir, 4)
            {
                int x = p.x + DX[dir];
                int y = p.y + DY[dir];
                if (rem_changes > 0 && in_rect(x, y, w, h) && !maze.outside(x, y) && !covered.get(x, y))
                {
                    --rem_changes;
                    maze.set(x, y, UTURN);
                    covered.set(x, y, true);
                }
            }
        }
        assert(rem_changes >= 0);

        return maze;
    }

private:
    const Maze& start_maze;
    const int max_changes;
    const vector<Pos> borders;
    const int w, h;

    int dist_to_outside[80][80];
};

class MazeFixing
{
public:
    vector<string> improve(const vector<string>& maze_, const int f)
    {
        g_timer.start();

        Maze start_maze(maze_);
        Solver solver(start_maze, f);
        Maze final_maze = solver.solve();

#ifdef VIS
//         save_image(make_filename(114514), final_maze, start_maze, final_maze.list_paths());
#endif

        vector<string> res;
        rep(y, start_maze.height()) rep(x, start_maze.width())
        {
            if (start_maze.get(x, y) != final_maze.get(x, y))
            {
                assert(!start_maze.outside(x, y));
                res.push_back(make_change(x, y, to_char(final_maze.get(x, y))));
            }
        }
        dump(start_maze.list_borders().size());
        dump(res.size());
        assert(res.size() <= f);

        fprintf(stderr, "time: %.3f\n", g_timer.get_elapsed());

//         dump(sum / g_timer.get_elapsed());

        return res;
    }

    string make_change(int x, int y, char v)
    {
        char buf[128];
        sprintf(buf, "%d %d %c", y, x, v);
        return buf;
    }
};


#ifdef LOCAL
void gen_input()
{
    int seed;
    cin >> seed;
    int h;
    cin >> h;
    vector<string> maze(h);
    input(maze, h);
    int f;
    cin >> f;

    ofstream os("input");
    os << seed << endl;
    os << h << endl;
    for (auto& s : maze)
        os << s << endl;
    os << endl;
    os << f << endl;
    os.close();
}
int main()
{
    cin >> test_case_seed;

    int h;
    cin >> h;
    vector<string> maze(h);
    input(maze, h);
    int f;
    cin >> f;

    auto ret = MazeFixing().improve(maze, f);
    cout << ret.size() << endl;
    for (auto& s : ret)
        cout << s << endl;
    cout.flush();
}
#endif
