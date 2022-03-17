# -*- coding: utf-8 -*-
import os
import random
import sys
import webbrowser

sys.stdin = open("../pyinput.txt")
sys.stdout = open("../OPENME.txt", 'w')


BROWSER_PATH = r'C:\Program Files\Mozilla Firefox\firefox.exe'


def big_enough():
    global ans
    return


def add_a_point(a, b, ch):
    global ans
    global last_point
    ans += f"{ch}{last_point[ch]}={(a, b)}\n"
    point_to_name[(a, b)] = ch + str(last_point[ch])
    last_point[ch] += 1
    big_enough()


def add_a_segment(a, b, c, d, ch):
    global ans
    if (a, b) not in point_to_name:
        add_a_point(a, b, 'A')
    if (c, d) not in point_to_name:
        add_a_point(c, d, 'A')
    ans += f"{ch}{last_point[ch]}=Segment({point_to_name[(a, b)]},{point_to_name[(c, d)]})\n"
    last_point[ch] += 1
    big_enough()


def perp(i, x, y):
    global ans
    ans += f"P{last_point['P']}=PerpendicularLine({point_to_name[(x, y)]},S{i})\n"
    last_point["P"] += 1
    big_enough()


def add_angle(i):
    global ans
    if (i == 1):
        ans += f"L{last_point['L']}=Angle(B{i + 1}, B{i}, X{i})\n"
        last_point['L'] += 1
        ans += f"L{last_point['L']}=Angle(B{k}, B{i}, X{i})\n"
        last_point['L'] += 1
    elif i == k:
        ans += f"L{last_point['L']}=Angle(B{1}, B{i}, X{i})\n"
        last_point['L'] += 1
        ans += f"L{last_point['L']}=Angle(B{i - 1}, B{i}, X{i})\n"
        last_point['L'] += 1
    else:
        ans += f"L{last_point['L']}=Angle(B{i + 1}, B{i}, X{i})\n"
        last_point['L'] += 1
        ans += f"L{last_point['L']}=Angle(B{i - 1}, B{i}, X{i})\n"
        last_point['L'] += 1
    big_enough()


def hide(i):
    global ans
    ans += f"SetConditionToShowObject(P{i},sl=={i})\n"
    big_enough()


def hide3(i):
    global ans
    ans += f"SetConditionToShowObject(X{i},sl=={str(i)})\n"
    big_enough()


def hide2(i):
    global ans
    ans += f"SetConditionToShowObject(L{2 * i - 1},sl=={i})\n"
    ans += f"SetConditionToShowObject(L{2 * i},sl=={i})\n"
    big_enough()


def add_point_on_line(i):
    global ans
    ans += f"X{last_point['X']}=Point(P{i})\n"
    last_point['X'] += 1
    big_enough()

# 'Количество отрезков'
n = int(input())
ans = ''
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
# Вершин в цикле
k = int(input())
cycle = []
for i in range(k):
    a, b, t = map(float, input().split())
    add_a_point(a, b, 'B')
    t = int(t)
    cycle.append((a, b, t))
    perp(t, a, b)
ans += f"sl=Slider(1,{k},1)\n"
for i in range(k):
    add_a_segment(cycle[i][0], cycle[i][1], cycle[(i + 1) %
                                                  k][0], cycle[(i + 1) % k][1],
                  'C')
for i in range(1, k + 1):
    hide(i)
for i in range(1, k + 1):
    add_point_on_line(i)
for i in range(1, k + 1):
    add_angle(i)
for i in range(1, k + 1):
    hide2(i)
    hide3(i)
print(ans,end='')
webbrowser.register('Firefox', None, webbrowser.BackgroundBrowser(BROWSER_PATH))
webbrowser.get(using='Firefox').open(f'file://{os.path.realpath("../show.html")}')
