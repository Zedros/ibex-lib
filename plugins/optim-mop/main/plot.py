import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animation
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import math
import os
import random
from stat import *
import multiprocessing
from tkinter import *
import time
from itertools import islice
# window=Tk()

def main():
    global filedate
    filedate = None
    q = multiprocessing.Queue()

    #Create and start the simulation process
    simulate=multiprocessing.Process(None,simulation,args=(q,))
    simulate.start()

    #Create the base plot
    plot()

    #Call a function to update the plot when there is new data
    updateplot(q)

    plt.show()
    print('Done')

def plot():
    global line,ax1,canvas
    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)

    # plt.ion()
    # canvas = FigureCanvasTkAgg(fig, master=window)
    # canvas.show()
    # canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
    # canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    line, = ax1.plot([], [])

def updateplot(q):
    try:       #Try to check if there is data in the queue
        result=q.get_nowait()

        if result !='Q':
            # print(result)
            LB, UB, LB2 = result
            UBx = []
            UBy = []
            LBx = []
            LBy = []
            ax1.clear()

            if len(LB) > 0:
                ax1.add_patch(patches.Rectangle(
                    (LB[0]['pts'][0], LB[0]['pts'][1]),
                    LB[0]['diam_x'], LB[0]['diam_y'],
                    fill=False,
                    edgecolor='red',
                    linestyle='solid',
                    lw=0.1
                    ))

            for p in [patches.Rectangle(
                        (lb['pts'][0], lb['pts'][1]),
                        lb['diam_x'], lb['diam_y'],
                        fill=False,
                        edgecolor='black',
                        linestyle='solid',
                        lw=0.1
                        ) for lb in islice(LB, 1, len(LB))]:
                ax1.add_patch(p)
            for ub in UB:
                UBx.append(ub[0])
                UBy.append(ub[1])
            for lb in LB2:
                LBx.append(lb[0])
                LBy.append(lb[1])
            for lb in LB:
            	if (lb['pA'][0] < lb['pB'][0]) and (lb['pA'][1] > lb['pB'][1]):
            		line = plt.Line2D((lb['pA'][0], lb['pB'][0]),(lb['pA'][1], lb['pB'][1]),lw=0.5,markeredgecolor='black')
            		ax1.add_line(line)
            ax1.plot()
            
            plt.plot(UBx, UBy, '-g', lw=0.5)
            plt.plot(LBx, LBy, '-b', lw=0.5)
            #ax1.set_xlim([0.0,0.35])
            #ax1.set_ylim([0.7,1.0])
            plt.pause(1)
            updateplot(q)
             # print(result)
                 #here get crazy with the plotting, you have access to all the global variables that you defined in the plot function, and have the data that the simulation sent.
             # line.set_ydata([1,result,10])
             # ax1.draw_artist(line)
             # plt.draw()
             # plt.pause(0.1)
             # updateplot(q)
             # window.after(500,updateplot,q)
        else:
             print('done')
    except:
        # print("empty")
        plt.pause(1)
        updateplot(q)
        # window.after(500,updateplot,q)

def simulation(q):
    global filedate
    # print(filedate)
    while True:
        st = os.stat('output2.txt')
        # print(st[ST_MTIME])
        if(st[ST_MTIME] == filedate):
            time.sleep(1)
            pass
        else:
            filedate = st[ST_MTIME]
            f = open("output2.txt")
            reader = f.read()
            reader = reader.replace('inf', "math.inf")
            reader = reader.replace('nan', "0")
            LB = {} #eval(reader.split("\n")[0])
            UB = eval(reader.split("\n")[0])
            LB2 = {}
            # LB2 = eval(reader.split("\n")[2])
            q.put((LB, UB, LB2))

    # iterations = range(100)
    # for i in iterations:
    #     if not i % 10:
    #         time.sleep(1)
    #             #here send any data you want to send to the other process, can be any pickable object
    #         q.put(random.randint(1,10))
    # q.put('Q')


if __name__ == '__main__':
    main()
