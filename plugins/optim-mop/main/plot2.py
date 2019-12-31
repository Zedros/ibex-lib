import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animation
import os
import math
import random
import subprocess
import multiprocessing
import traceback
import argparse


class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __lt__(self, other):
        if self.x < other.x:
            return True
        elif self.x == other.x:
            return self.y < other.y
        return False


class Box:
    def __init__(self, box_string):
        box_string = box_string.replace('inf', "math.inf")
        box_string = box_string.replace('nan', "0")
        box_var = eval(box_string)
        self.createBox(box_var['id'], box_var['pts'][0], box_var['pts'][1],
                       box_var['diam_x'], box_var['diam_y'], box_var['pA'],
                       box_var['pB'])

    def __repr__(self):
        return str((self.x, self.y))

    def __lt__(self, other):
        if self.x < other.x:
            return True
        elif self.x == other.x:
            return self.y < other.y
        return False

    def createBox(self, id_box, x, y,  diam_x, diam_y, pA, pB):
        self.id_box = id_box
        self.x = x
        self.y = y
        self.diam_x = diam_x
        self.diam_y = diam_y
        self.pA = pA
        self.pB = pB


class LB:
    def __init__(self, lb_string):
        lb_string = lb_string.replace('inf', 'math.inf')
        lb_string = lb_string.replace('nan', '0')
        lb_var = eval(lb_string)
        self.key = lb_var['id']
        self.x = lb_var['pts'][0]
        self.y = lb_var['pts'][1]

    def __lt__(self, other):
        if isinstance(other, LB):
            if self.x < other.x:
                return True
            elif self.x == other.x:
                return self.y > other.y
            else:
                return False
        return NotImplemented


class UB:
    def __init__(self, ub_string):
        ub_string = ub_string.replace('inf', 'math.inf')
        ub_string = ub_string.replace('nan', '0')
        ub_var = eval(ub_string)
        self.key = ub_string
        self.x = ub_var['pts'][0]
        self.y = ub_var['pts'][1]


class DataDict:
    box_set = {}
    patch_set = {}
    UB_set = {}
    LB_set = {}
    pending_box = None

    def get_ub(self):
        ub_list_x = []
        ub_list_y = []
        ub_list = []
        for x in self.UB_set.values():
            ub_list.append(Point(x.x, x.y))
        ub_list = sorted(ub_list)
        i = 0
        if ub_list:
            ub_list_x.append(ub_list[0].x)
            ub_list_y.append(ub_list[0].y)
            for point in ub_list[1:]:
                ub_list_x.append(point.x)
                ub_list_y.append(ub_list[i].y)
                ub_list_x.append(point.x)
                ub_list_y.append(point.y)
                i = i + 1
            box = sorted(list(self.box_set.values()))[-1]
            ub_list_x.append(box.x + box.diam_x)
            ub_list_y.append(ub_list_y[-1])
        return ub_list_x, ub_list_y

    def get_lb(self):
        lb_list_x = []
        lb_list_y = []
        box_list = []
        for x in self.LB_set.values():
            box_list.append(Point(x.x, x.y))
        for x in self.box_set.values():
            box_list.append(Point(x.x, x.y))
        box_list = sorted(box_list)
        i = 0
        for box in box_list[1:]:
            if box.y > box_list[i].y:
                box_list.remove(box)
            else:
                i = i + 1
        i = 0
        if box_list:
            lb_list_x.append(box_list[0].x)
            lb_list_y.append(box_list[0].y)
            for box in box_list[1:]:
                if box.x != box_list[i].x:
                    lb_list_x.append(box.x)
                    lb_list_y.append(box_list[i].y)
                lb_list_x.append(box.x)
                lb_list_y.append(box.y)
                i = i + 1
            box = sorted(list(self.box_set.values()))[-1]
            lb_list_x.append(box.x + box.diam_x)
            lb_list_y.append(lb_list_y[-1])
        return lb_list_x, lb_list_y

    def get_cy_lines(self):
        line_list = []
        for box in self.box_set.values():
            if(box.pA[0] < box.pB[0]) and (box.pA[1] > box.pB[1]):
                line = plt.Line2D((box.pA[0], box.pB[0]),
                                  (box.pA[1], box.pB[1]),
                                  lw=0.5, markeredgecolor='black')
                line_list.append(line)
        return line_list

    def append_lb(self, lb):
        self.LB_set[lb.key] = lb

    def append_ub(self, ub):
        self.UB_set[ub.key] = ub

    def remove_ub(self, key):
        key = key.replace('inf', 'math.inf')
        key = key.replace('nan', '0')
        del self.UB_set[key]

    def append_box(self, box):
        self.box_set[box.id_box] = box
        self.patch_set[box.id_box] = plt.Rectangle(
            (box.x, box.y),
            box.diam_x, box.diam_y,
            fill=False,
            edgecolor='black',
            linestyle='solid',
            lw=0.1,
            label="Box"
        )

    def get_patches(self):
        patch_list = []
        for x in self.patch_set.values():
            patch_list.append(x)
        return patch_list

    def get_box(self, id_box):
        return self.box_set[id_box], self.patch_set[id_box]

    def remove_box(self, id_box):
        if self.pending_box is not None and id_box != self.pending_box:
            del self.box_set[self.pending_box]
            del self.patch_set[self.pending_box]
            self.pending_box = id_box
        else:
            self.pending_box = id_box


def ibex(q, args):
    params = args.ibex.split(" ")
    proc = subprocess.Popen(
        ["__build__/plugins/optim-mop/ibexmop", "--plot"] + params,
        stdin=None, stdout=subprocess.PIPE
    )
    state = "run"
    if "--help" in params or "-h" in params:
        print(proc.stdout.read().decode())
    else:
        if 'variable' in proc.stdout.readline().decode():
            q.put('start')
        else:
            q.put('stop')
            state = 'stop'
        while state != "stop":
            for line in proc.stdout:
                if 'add' in line.decode() or 'del' in line.decode():
                    q.put(line.decode())
            proc.stdout.close()
            proc.wait()
            break


def plot():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # ax.set_aspect(1)
    line, = plt.plot([], [], 'b-', label='Upperbound')
    line2, = plt.plot([], [], 'r-', label='Lowerbound')
    plt.xlabel('Funcion Objetivo 1')
    plt.ylabel('Funcion Objetivo 2')
    plt.title('Plotter')
    # plt.setp(line, color='b', linewidth=2.0, label='UB', marker='-')
    # line.set_data([], [])
    plt.xlim(-300, 300)
    plt.ylim(-300, 300)
    return fig, line, ax, line2


def updateplot(num, q, l, data, ax, patch_list, l2, line_list):
    if pause:
        patch_list = [ax.add_patch(x) for x in data.get_patches()]
        return tuple(patch_list) + (l, l2,) + tuple(line_list)
    else:
        try:
            result = q.get_nowait()
            if result != 'Q':
                # ax.clear()
                if 'add:' in result:
                    box = Box(result.split("add: ")[1])
                    data.append_box(box)
                if 'del:' in result:
                    id_box = eval(result.split("del: ")[1])['id']
                    data.remove_box(id_box)
                if 'add ub:' in result:
                    ub = UB(result.split("add ub: ")[1])
                    data.append_ub(ub)
                if 'del ub:' in result:
                    key = result.split("del ub: ")[1]
                    data.remove_ub(key)
                if 'add lb:' in result:
                    lb = LB(result.split("add lb: ")[1])
                    data.append_lb(lb)
                l.set_data(data.get_ub())
                patch_list = [ax.add_patch(x) for x in data.get_patches()]
                line_list = [ax.add_line(x) for x in data.get_cy_lines()]
                l2.set_data(data.get_lb())
                # plt.legend(handles=[l, l2, patch_list[0]])
                return tuple(patch_list) + (l, ) + (l2, ) + tuple(line_list)
            else:
                q.close()
        except Exception as e:
            # print("Esto es un error:")
            # print(e)
            return tuple(patch_list) + (l, ) + (l2, ) + tuple(line_list)


pause = False


def onClick(event, line, line2):
    global pause
    if(event.button == 2):
        if min(line.get_xdata()[:-1]) < min(line2.get_xdata()):
            x_min = min(line.get_xdata())
        else:
            x_min = min(line2.get_xdata())
        if max(line.get_xdata()[:-1]) > max(line2.get_xdata()):
            x_max = max(line.get_xdata())
        else:
            x_max = max(line2.get_xdata())
        y_min = min(line2.get_ydata())
        y_max = max(line.get_ydata())
        x_prom = (x_max + x_min)/2
        y_prom = (y_max + y_min)/2
        if x_max - x_min > y_max - y_min:
            x_a = x_min - (x_max - x_min)*.1
            x_b = x_max + (x_max - x_min)*.1
            y_a = x_a - (x_prom - y_prom)
            y_b = x_b - (x_prom - y_prom)
        else:
            y_a = y_min - (y_max - y_min)*.1
            y_b = y_max + (y_max - y_min)*.1
            x_a = y_a - (y_prom - x_prom)
            x_b = y_b - (y_prom - x_prom)
        plt.xlim(x_a,
                 x_b)
        plt.ylim(y_a,
                 y_b)
    if(event.button == 3):
        pause ^= True


def main():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-ib", "--ibex",
                        help="Ibex arguments",
                        type=str, required=True)
    parser.add_argument("-s", "--speed",
                        help="Speed of the plotter",
                        type=int, default=100)
    args = parser.parse_args()
    try:
        q = multiprocessing.Queue()
        ibex_proc = multiprocessing.Process(None, ibex, args=(q, args))
        ibex_proc.start()
        if "-h" not in args.ibex and "--h" not in args.ibex:
            while True:
                result = q.get()
                if result == 'start':
                    break
                elif result == 'stop':
                    ibex_prox.terminate()
                    return 0
            fig, line, ax, line2 = plot()
            patch_list = []
            line_list = []
            fig.canvas.mpl_connect('button_press_event',
                                   lambda event: onClick(event, line, line2))
            data = DataDict()
            line_ani = animation.FuncAnimation(
                fig, updateplot, q.qsize(),
                fargs=(q, line, data, ax, patch_list, line2, line_list),
                interval=args.speed, blit=True, repeat=False
            )
            plt.show()
        ibex_proc.join()
        print("Done")
        ibex_proc.terminate()
    except Exception:
        ibex_proc.terminate()


if __name__ == '__main__':
    main()
