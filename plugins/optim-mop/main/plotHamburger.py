import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animation
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib import collections  as mc
import math
import os
import random
from stat import *
import multiprocessing
from tkinter import *
import time
from itertools import islice
# window=Tk()
filedate = None

def main():
    global filedate
    # filedate = None
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
            LB, UB, LB2, LB3 = result
            UBx = []
            UBy = []
            LBx = []
            LBy = []

            ax1.clear()

            for ub in UB:
                UBx.append(ub[0])
                UBy.append(ub[1])

            # add line upperbound
            lines = []
            #p = [(-100, 100), (100, -100)]
            if len(LB)>0:
                lines.append(LB)
                lc = mc.LineCollection(lines, colors='blue', linewidths=0.5)
                ax1.add_collection(lc)

            # add lines function
            lines = []
            # p = [(100,100),(20,60),(-100,-100)]
            if len(LB2)>0:
                for i in range(1,len(LB2)):
                    arr1 = []
                    arr1.append(LB2[i-1])
                    arr1.append(LB2[i])
                    lines.append(arr1)
                lc = mc.LineCollection(lines, colors='green', linewidths=0.5)
                ax1.add_collection(lc)

            # add lines function
            lines2 = []
            # p = [(100,100),(20,60),(-100,-100)]
            if len(LB3)>0:
                for i in range(1,len(LB3)):
                    arr1 = []
                    arr1.append(LB3[i-1])
                    arr1.append(LB3[i])
                    lines2.append(arr1)
                lc = mc.LineCollection(lines2, colors='brown', linewidths=0.5)
                ax1.add_collection(lc)


            ax1.plot()
            plt.plot(UBx, UBy, 'r-', markersize=3)
            # plt.plot(LBx, LBy, '-b', lw=0.5)
            #ax1.set_xlim([0.0,0.35])
            #ax1.set_ylim([0.7,1.0])
            #plt.gca().set_aspect('equal', adjustable='box')
            plt.axis('scaled')
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
        st = os.stat('outputH.txt')
        # print(st[ST_MTIME])
        if(st[ST_MTIME] == filedate):
            time.sleep(1)
            pass
        else:
            filedate = st[ST_MTIME]
            f = open("outputH.txt")
            reader = f.read()
            reader = reader.replace('inf', "math.inf")
            reader = reader.replace('nan', "0")
            UB = eval(reader.split("\n")[0])
            LB = eval(reader.split("\n")[1])
            LB2 = eval(reader.split("\n")[2])
            LB3 = eval(reader.split("\n")[3])

            # LB2 = eval(reader.split("\n")[2])
            q.put((LB, UB, LB2, LB3))

    # iterations = range(100)
    # for i in iterations:
    #     if not i % 10:
    #         time.sleep(1)
    #             #here send any data you want to send to the other process, can be any pickable object
    #         q.put(random.randint(1,10))
    # q.put('Q')


if __name__ == '__main__':
    main()
