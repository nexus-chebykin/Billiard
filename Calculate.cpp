#include <iostream>
#include <vector>
#include <algorithm>
#include <deque>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <stack>
#include <tuple>
#include <set>

#define all(x) (x).begin(), (x).end()

using namespace std;

typedef pair<int, int> pii;

const int N = 10; // Максимальное количество отрезков!
const int segmentPrecision = 1000;
const int anglePrecision = 1000;
const double pi = 3.141592653589793;
const double INF = 1e9;

struct Point {
    double x = 0, y = 0;

    Point() = default;
    Point(const double& a, const double& b) {
        x = a;
        y = b;
    }
};

inline double dist(Point fst, Point oth) // Расстояние между двумя точками
{
    return hypot(fst.x - oth.x, fst.y - oth.y);
}

inline pair<pair<double, double>, double> line_from_points(Point F, Point S) // Прямая по 2 точкам
{
    pair<double, double> tmp = {F.y - S.y, S.x - F.x};
    return {tmp, F.x * S.y - F.y * S.x};
}

Point cross_lines(pair<pair<double, double>, double> F,
                  pair<pair<double, double>, double> S) //Пересекает две прямые (если параллельные => (INF, INF))
{
    double A1 = F.first.first;
    double B1 = F.first.second;
    double C1 = F.second;
    double A2 = S.first.first;
    double B2 = S.first.second;
    double C2 = S.second;
    if (A1 * B2 - A2 * B1 == 0)
        return Point(INF, INF);
    return Point((C2 * B1 - C1 * B2) / (A1 * B2 - A2 * B1), (C2 * A1 - C1 * A2) / (B1 * A2 - A1 * B2));
}

double sp(Point F, Point S) // Скалярное произведение
{
    return F.x * S.x + F.y * S.y;
}

double vp(Point F, Point S) // Модуль векторного произведения
{
    return F.x * S.y - F.y * S.x;
}

Point intersect(double dx,
                double dy,
                Point st,
                pair<Point, Point> seg) // Пересекает луч из st с вектором направления (dx, dy) с отрезком seg
{
    Point second_beam_point = Point(st.x + dx, st.y + dy);
    auto beam = line_from_points(st, second_beam_point);
    auto seg_line = line_from_points(seg.first, seg.second);
    Point intersection = cross_lines(beam, seg_line);
    if (intersection.x >= INF) // Параллельные луч и отрезок
        return Point(INF, INF);
    // Проверка, что точка пересечения лежит внутри отрезка
    Point v1 = Point(intersection.x - seg.first.x, intersection.y - seg.first.y);
    Point v2 = Point(intersection.x - seg.second.x, intersection.y - seg.second.y);
    if (sp(v1, v2) > 0)
        return Point(INF, INF);
    // Проверка, что точка пересечения лежит в нужной стороне от начала луча
    v1 = Point(dx, dy);
    v2 = Point(intersection.x - st.x, intersection.y - st.y);
    if (sp(v1, v2) <= 0)
        return Point(INF, INF);
    return intersection;
}

pair<pair<double, double>, double> line_from_vector(Point V, Point F) // Прямая по точке и направляющему вектору
{
    pair<pair<double, double>, double> line = {{-V.y, V.x}, 0};
    line.second = -(line.first.first * F.x + line.first.second * F.y);
    return line;
}

Point reflect(Point F, pair<pair<double, double>, double> line) //Отражает точку F относительно прямой line
{
    double A = line.first.first;
    double B = line.first.second;
    double C = line.second;
    pair<pair<double, double>, double> perp_line = line_from_vector(Point(A, B), F);
    Point inter = cross_lines(line, perp_line);
    return Point(2 * inter.x - F.x, 2 * inter.y - F.y);
}

void normalize(Point& F) // Нормализует вектор F
{
    double tmp = hypot(F.x, F.y);
    F.x /= tmp;
    F.y /= tmp;
}

pair<Point, Point> segm[N];
vector<pair<pii, int>> g[N][segmentPrecision + 1][anglePrecision];
pair<double, double> angles[anglePrecision]; //sin _ cos
double segment_length[N];

int was[N][segmentPrecision + 1][anglePrecision];
vector<pair<pii, int>> good_vertex;
set<pair<pii, int>> cycle_vertex;

pair<pii, int> par[N][segmentPrecision + 1][anglePrecision];
int dst[N][segmentPrecision + 1][anglePrecision];
int best_cycle_length = INF;
vector<pair<pii, int>> best_cycle, exact_length_cycle;

Point kth_point_on_segm(int k, int ind) //Возвращает к-ую точку на отрезке
{
    double rati = (double) k / segmentPrecision;
    return Point(segm[ind].first.x * (1 - rati) + segm[ind].second.x * rati,
                 segm[ind].first.y * (1 - rati) + segm[ind].second.y * rati);
}

void dfs(int segm_id, int point, int angle) {
    if (good_vertex.size() >= 10000)
        return;
    was[segm_id][point][angle] = 1;
    bool f = false;
    for (pair<pii, int> nex: g[segm_id][point][angle]) {
        if (was[nex.first.first][nex.first.second][nex.second] == 0)
            dfs(nex.first.first, nex.first.second, nex.second);
        else if (was[nex.first.first][nex.first.second][nex.second] == 1)
            f = true;
    }
    if (f)
        good_vertex.push_back({{segm_id, point}, angle});
    was[segm_id][point][angle] = 2;
}

void dfs_non_rec(int segm_id, int point, int angle) {
    vector<pair<pii, pii>> s;
    s.push_back({{segm_id, point}, {angle, 0}});
    was[segm_id][point][angle] = 1;
    while (!s.empty()) {
        int seg = s.back().first.first;
        int p = s.back().first.second;
        int a = s.back().second.first;
        int ind = s.back().second.second;
        s.pop_back();
        if (ind < g[seg][p][a].size())
            s.push_back({{seg, p}, {a, ind + 1}});
        else {
            was[seg][p][a] = 2;
            continue;
        }
        auto& nex = g[seg][p][a][ind];
        if (was[nex.first.first][nex.first.second][nex.second] == 0) {
            was[nex.first.first][nex.first.second][nex.second] = 1;
            s.push_back({{nex.first.first, nex.first.second}, {nex.second, 0}});
        } else if (was[nex.first.first][nex.first.second][nex.second] == 1)
            cycle_vertex.insert({{seg, p}, a});
    }
}

void bfs(int segm_id, int point, int angle) {
    deque<pair<pii, int>> dq;
    dq.push_back({{segm_id, point}, angle});
    int iter = 0;
    int max_dist = 0;
    vector<pair<pii, int>> used_vertex;
    dst[segm_id][point][angle] = 0;
    pair<pii, int> final_vert = {{-1, -1}, -1};
    while (max_dist < best_cycle_length && dq.size() > 0) {
        pair<pii, int> cur_v = dq.front();
        dq.pop_front();
        for (auto el: g[cur_v.first.first][cur_v.first.second][cur_v.second]) {
            if (el == make_pair(make_pair(segm_id, point), angle)) {
                final_vert = cur_v;
                break;
            }
            if (dst[el.first.first][el.first.second][el.second]
                > dst[cur_v.first.first][cur_v.first.second][cur_v.second] + 1) {
                dst[el.first.first][el.first.second][el.second] = max_dist =
                    dst[cur_v.first.first][cur_v.first.second][cur_v.second] + 1;
                used_vertex.push_back(el);
                par[el.first.first][el.first.second][el.second] = cur_v;
                dq.push_back(el);
            }
        }
    }
    dst[segm_id][point][angle] = INF;
    for (auto el: used_vertex)
        dst[el.first.first][el.first.second][el.second] = INF;
    vector<pair<pii, int>> cycle;
    if (final_vert.first.first != -1) {
        while (final_vert != make_pair(make_pair(segm_id, point), angle)) {
            cycle.push_back(final_vert);
            final_vert = par[final_vert.first.first][final_vert.first.second][final_vert.second];
        }
        cycle.push_back(final_vert);
        if (cycle.size() < best_cycle_length) {
            best_cycle_length = cycle.size();
            best_cycle = cycle;
        }
    }
}

void bfs_exact_length(int len, int segm_id, int point, int angle) {
    deque<pair<pii, int>> dq;
    dq.push_back({{segm_id, point}, angle});
    int iter = 0;
    int max_dist = 0;
    vector<pair<pii, int>> used_vertex;
    dst[segm_id][point][angle] = 0;
    pair<pii, int> final_vert = {{-1, -1}, -1};
    while (max_dist <= len && dq.size() > 0) {
        pair<pii, int> cur_v = dq.front();
        dq.pop_front();
        for (auto el: g[cur_v.first.first][cur_v.first.second][cur_v.second]) {
            if (el == make_pair(make_pair(segm_id, point), angle)) {
                final_vert = cur_v;
                break;
            }
            if (dst[el.first.first][el.first.second][el.second]
                > dst[cur_v.first.first][cur_v.first.second][cur_v.second] + 1) {
                dst[el.first.first][el.first.second][el.second] = max_dist =
                    dst[cur_v.first.first][cur_v.first.second][cur_v.second] + 1;
                used_vertex.push_back(el);
                par[el.first.first][el.first.second][el.second] = cur_v;
                dq.push_back(el);
            }
        }
    }
    dst[segm_id][point][angle] = INF;
    for (auto el: used_vertex)
        dst[el.first.first][el.first.second][el.second] = INF;
    vector<pair<pii, int>> cycle;
    if (final_vert.first.first != -1) {
        while (final_vert != make_pair(make_pair(segm_id, point), angle)) {
            cycle.push_back(final_vert);
            final_vert = par[final_vert.first.first][final_vert.first.second][final_vert.second];
        }
        cycle.push_back(final_vert);
        if (cycle.size() == len) {
            //best_cycle_length = cycle.size();
            exact_length_cycle = cycle;
        }
    }
}

signed main() {
    cout << "Reading walls...\n";

    ifstream cin("../my_test.txt");
    int n;
    cin >> n;
    cout << setprecision(7) << fixed;
    // Предподсчет координат векторов углов (все углы в программе рассматриваются, будто отложены от начала координат)
    for (int i = 0; i < anglePrecision; ++i) {
        double tmp = 2 * pi * i / anglePrecision;
        angles[i] = {sin(tmp), cos(tmp)};
    }
    // Ввод отрезков
    for (int i = 0; i < n; ++i)
        cin >> segm[i].first.x >> segm[i].first.y >> segm[i].second.x >> segm[i].second.y;

    //Корректен ли был ввод
    cout << n << ' ' << segm[0].first.x << endl;

    // Длины отрезков
    for (int i = 0; i < n; ++i)
        segment_length[i] = dist(segm[i].first, segm[i].second);
    double min_length = (*min_element(segment_length, segment_length + n)) / 10;

    cout << "Building graph...\n";
    for (int currentSegment = 0; currentSegment < n; ++currentSegment) {
        Point cur = Point(segm[currentSegment].first.x, segm[currentSegment].first.y);
        const double step_x = (segm[currentSegment].second.x - segm[currentSegment].first.x) / segmentPrecision;
        const double step_y = (segm[currentSegment].second.y - segm[currentSegment].first.y) / segmentPrecision;
        for (int pt = 0; pt <= segmentPrecision; ++pt) {
            if (min(dist(cur, segm[currentSegment].first), dist(cur, segm[currentSegment].second)) <= min_length)
                goto Label;
            for (int currentAngle = 0; currentAngle < anglePrecision; ++currentAngle) {
                // Поиск первого пересекаемого отрезка
                pair<double, int> intersected_segm = {INF, -1};
                for (int otherSegment = 0; otherSegment < n; ++otherSegment) {
                    if (otherSegment != currentSegment) {
                        Point inter =
                            intersect(angles[currentAngle].second, angles[currentAngle].first, cur, segm[otherSegment]);
                        if (inter.x != INF) {
                            double distanceToIntersection = dist(cur, inter);
                            if (intersected_segm.first > distanceToIntersection)
                                intersected_segm = {distanceToIntersection, otherSegment};
                        }
                    }
                }
                if (intersected_segm.second != -1) // Нашелся отрезок, об который ударится мячик, работаем с ним
                {
                    Point inter = intersect(angles[currentAngle].second,
                                            angles[currentAngle].first,
                                            cur,
                                            segm[intersected_segm.second]);
                    double tmp = dist(segm[intersected_segm.second].first,
                                      inter); // Находим расстояние от начала пересеченного отрезка до точки пересечения, чтобы узнать, какие точки разбиения ближайшие
                    tmp *= segmentPrecision / segment_length[intersected_segm.second];
                    int fragm_point = floor(tmp); // Ближайшие точки разбиения теперь fragm_point и fragm_point + 1

                    // Рассмотрим случай, если точка разбиения, в которую попадет мячик - fragm_point
                    /*double rati = (double)fragm_point / segmentPrecision;
                    // Находим координаты этой точки разбиения
                    Point other_segm_point = Point(segm[intersected_segm.second].first.x * (1 - rati) + segm[intersected_segm.second].second.x * rati,
                                                    segm[intersected_segm.second].first.y * (1 - rati) + segm[intersected_segm.second].second.y * rati);
                    */
                    Point other_segm_point = kth_point_on_segm(fragm_point, intersected_segm.second);
                    auto segment_line =
                        line_from_points(segm[intersected_segm.second].first, segm[intersected_segm.second].second);
                    // Прямая через рассматриваемую точку разбиения, перпендикулярная к отрезку
                    auto reflecting_line =
                        line_from_vector(Point(segment_line.first.first, segment_line.first.second), other_segm_point);
                    Point refl_point = reflect(cur, reflecting_line);
                    Point new_trajectory = Point(refl_point.x - other_segm_point.x, refl_point.y - other_segm_point.y);
                    normalize(new_trajectory);
                    // Поиск угла отражения
                    double angle;
                    if (new_trajectory.x >= 0) {
                        if (new_trajectory.y >= 0)
                            angle = asin(new_trajectory.y) / (2 * pi);
                        else
                            angle = 1 + asin(new_trajectory.y) / (2 * pi);
                    } else {
                        if (new_trajectory.y >= 0)
                            angle = acos(new_trajectory.x) / (2 * pi);
                        else
                            angle = 0.5 - asin(new_trajectory.y) / (2 * pi);
                    }
                    angle *= anglePrecision;
                    int angle_number = floor(angle);
                    // Проверка что нет пары ребер a -> b, b -> a
                    pair<pii, int> to1 = {{intersected_segm.second, fragm_point}, angle_number};
                    pair<pii, int> to2 = {{intersected_segm.second, fragm_point}, (angle_number + 1) % anglePrecision};
                    pair<pii, int> from = {{currentSegment, pt}, currentAngle};
                    if (find(all(g[to1.first.first][to1.first.second][to1.second]), from)
                        == g[to1.first.first][to1.first.second][to1.second].end())
                        g[currentSegment][pt][currentAngle].push_back(to1);
                    if (find(all(g[to2.first.first][to2.first.second][to2.second]), from)
                        == g[to2.first.first][to2.first.second][to2.second].end())
                        g[currentSegment][pt][currentAngle].push_back(to2);
                    if (fragm_point != segmentPrecision) {
                        // Рассматриваем случай fragm_point + 1, если он есть! (этот if)
                        ++fragm_point;
                        /*double rati = (double)fragm_point / segmentPrecision;
                        // Находим координаты этой точки разбиения
                        Point other_segm_point = Point(segm[intersected_segm.second].first.x * (1 - rati) + segm[intersected_segm.second].second.x * rati,
                                                        segm[intersected_segm.second].first.y * (1 - rati) + segm[intersected_segm.second].second.y * rati);
                        */
                        Point other_segm_point = kth_point_on_segm(fragm_point, intersected_segm.second);
                        // Прямая через рассматриваемую точку разбиения, перпендикулярная к отрезку
                        auto reflecting_line =
                            line_from_vector(Point(segment_line.first.first, segment_line.first.second),
                                             other_segm_point);
                        Point refl_point = reflect(cur, reflecting_line);
                        Point new_trajectory =
                            Point(refl_point.x - other_segm_point.x, refl_point.y - other_segm_point.y);
                        normalize(new_trajectory);
                        // Поиск угла отражения
                        double angle;
                        if (new_trajectory.x >= 0) {
                            if (new_trajectory.y >= 0)
                                angle = asin(new_trajectory.y) / (2 * pi);
                            else
                                angle = 1 + asin(new_trajectory.y) / (2 * pi);
                        } else {
                            if (new_trajectory.y >= 0)
                                angle = acos(new_trajectory.x) / (2 * pi);
                            else
                                angle = 0.5 - asin(new_trajectory.y) / (2 * pi);
                        }
                        angle *= anglePrecision;
                        int angle_number = floor(angle);
                        // Проверка что нет пары ребер a -> b, b -> a
                        pair<pii, int> to1 = {{intersected_segm.second, fragm_point}, angle_number};
                        pair<pii, int>
                            to2 = {{intersected_segm.second, fragm_point}, (angle_number + 1) % anglePrecision};
                        pair<pii, int> from = {{currentSegment, pt}, currentAngle};
                        if (find(all(g[to1.first.first][to1.first.second][to1.second]), from)
                            == g[to1.first.first][to1.first.second][to1.second].end())
                            g[currentSegment][pt][currentAngle].push_back(to1);
                        if (find(all(g[to2.first.first][to2.first.second][to2.second]), from)
                            == g[to2.first.first][to2.first.second][to2.second].end())
                            g[currentSegment][pt][currentAngle].push_back(to2);
                    }
                }
            }
            Label:
            cur.x += step_x;
            cur.y += step_y;
        }
    }

    cout << "Running DFS...\n";
    cout << n << endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < segmentPrecision + 1; ++j)
            for (int k = 0; k < anglePrecision; ++k)
                if (was[i][j][k] == 0)
                    dfs_non_rec(i, j, k);
    }

    cout << "cycle_vertex.size " << cycle_vertex.size() << endl;

    cout << "Running BFS for shortest cycle...\n";
    for (auto el: cycle_vertex)
        bfs(el.first.first, el.first.second, el.second);
    ofstream fout("../pyinput.txt");
    fout << setprecision(7) << fixed;

    fout << n << endl;
    for (int i = 0; i < n; ++i)
        fout << segm[i].first.x << ' ' << segm[i].first.y << ' ' << segm[i].second.x << ' ' << segm[i].second.y << endl;
    if (!best_cycle.empty()) {
        fout << best_cycle.size() << endl;
        for (auto el: best_cycle) {
            Point tmp = kth_point_on_segm(el.first.second, el.first.first);
            fout << tmp.x << " " << tmp.y << ' ' << el.first.first + 1 << endl;
        }
    } else
        fout << "No cycle found :(\n";
    fout.close();
    system("py -3 ../Visualize.py");
    return 0;
}

