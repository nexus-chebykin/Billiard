import random



def big_enough():
    global ans
    if (len(ans) >= 300):
        ans = ans[:-1] + "})"
        print(ans)
        ans = "Execute({"


def add_a_point(a, b, ch):
    global ans
    global last_point
    ans += '"' + ch + str(last_point[ch]) + '=' + str((a, b)) + '",'
    point_to_name[(a, b)] = ch + str(last_point[ch])
    last_point[ch] += 1
    big_enough()


def add_a_segment(a, b, c, d, ch):
    global ans
    if (a, b) not in point_to_name:
        add_a_point(a, b, 'A')
    if (c, d) not in point_to_name:
        add_a_point(c, d, 'A')
    ans += '"' + ch + str(last_point[ch]) + '=Segment(' + \
        point_to_name[(a, b)] + ',' + point_to_name[(c, d)] + ')",'
    last_point[ch] += 1
    big_enough()


def perp(i, x, y):
    global ans
    ans += '"P' + str(last_point["P"]) + '=PerpendicularLine(' + \
        point_to_name[(x, y)] + ',S' + str(i) + ')",'
    last_point["P"] += 1
    big_enough()


def add_angle(i):
    global ans
    if (i == 1):
        ans += '"L' + str(last_point["L"]) + \
            '=Angle(B' + str(i + 1) + ',' + 'B' + str(i) + ',' + 'X' + str(i) + ')",'
        last_point['L'] += 1
        ans += '"L' + str(last_point["L"]) + \
            '=Angle(B' + str(k) + ',' + 'B' + str(i) + ',' + 'X' + str(i) + ')",'
        last_point['L'] += 1
    elif i == k:
        ans += '"L' + str(last_point["L"]) + \
            '=Angle(B' + str(1) + ',' + 'B' + str(i) + ',' + 'X' + str(i) + ')",'
        last_point['L'] += 1
        ans += '"L' + str(last_point["L"]) + \
            '=Angle(B' + str(i - 1) + ',' + 'B' + str(i) + ',' + 'X' + str(i) + ')",'
        last_point['L'] += 1
    else:
        ans += '"L' + str(last_point["L"]) + \
            '=Angle(B' + str(i + 1) + ',' + 'B' + str(i) + ',' + 'X' + str(i) + ')",'
        last_point['L'] += 1
        ans += '"L' + str(last_point["L"]) + \
            '=Angle(B' + str(i - 1) + ',' + 'B' + str(i) + ',' + 'X' + str(i) + ')",'
        last_point['L'] += 1
    big_enough()


def hide(i):
    global ans
    ans += '"SetConditionToShowObject(P' + str(i) + ',sl==' + str(i) + ')",'
    big_enough()

def hide3(i):
    global ans
    ans += '"SetConditionToShowObject(X' + str(i) + ',sl==' + str(i) + ')",'
    big_enough()

def hide2(i):
    global ans
    ans += '"SetConditionToShowObject(L' + \
        str(2 * i - 1) + ',sl==' + str(i) + ')",'
    ans += '"SetConditionToShowObject(L' + \
        str(2 * i) + ',sl==' + str(i) + ')",'
    big_enough()

def add_point_on_line(i):
    global ans
    ans += '"X' + str(last_point['X']) + '=Point(P' + str(i) + ')",'
    last_point['X'] += 1
    big_enough() 


n = int(input('Количество отрезков:\n'))
ans = 'Execute({'
last_point = dict()
last_point['A'] = 1  # Вершины многоугольника
last_point['B'] = 1  # Вершины цикла
last_point['S'] = 1  # Название сторон многоугольника
last_point['C'] = 1  # Название сторон цикла
last_point['P'] = 1  # Название перпендикуляров
last_point['L'] = 1  # Название углов
last_point['X'] = 1  # Название точек на перпендикулярах
point_to_name = dict()
segments = []
for i in range(n):
    a, b, c, d = map(float, input().split())
    segments.append((a, b, c, d))
    add_a_segment(a, b, c, d, "S")
k = int(input('Количество вершин цикла:\n'))
cycle = []
for i in range(k):
    a, b, t = map(float, input().split())
    add_a_point(a, b, 'B')
    t = int(t)
    cycle.append((a, b, t))
    perp(t, a, b)
ans += '"sl=Slider(1,' + str(k) + ',1)",'
for i in range(k):
    add_a_segment(cycle[i][0], cycle[i][1], cycle[(i + 1) %
                                                  k][0], cycle[(i + 1) % k][1], 'C')
for i in range(1, k + 1):
    hide(i)
for i in range(1, k + 1):
    add_point_on_line(i)
for i in range(1, k + 1):
    add_angle(i)
for i in range(1, k + 1):
    hide2(i)
    hide3(i)
ans = ans[:-1] + '})'
print(repr(ans))