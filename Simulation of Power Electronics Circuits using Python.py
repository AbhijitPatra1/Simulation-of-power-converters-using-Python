"""
Created on Fri Apr 15 10:50:32 2022

@author: Abhijt Patra
        
"""



import numpy as np
import tkinter as tk
import matplotlib.pyplot as plt
import math

from pylab import *
from scipy.optimize import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from tkinter import ttk

np.seterr(divide='ignore', invalid='ignore')

win = tk.Tk()
win.title("Simulation of Power Converters")
win.minsize(300, 570)

a_title = tk.Label(win, text="SIMULATION OF POWER CONVERTERS", fg="blue", bg="yellow", font="Helvetica 16 bold italic")

a_title.grid(column=0, row=0, columnspan=3)

simpara = tk.LabelFrame(win, text='Simulation Parameters', fg="brown", font="Verdana 8 bold")

simpara.grid(column=0, row=1, sticky=tk.W + tk.N, padx=10, pady=10)

contypes = tk.LabelFrame(win, text='Select One in each group, Simulate and Plot(full)', fg="brown",
                         font="Verdana 8 bold")

contypes.grid(column=0, row=2, sticky=tk.W + tk.N, padx=10, pady=10)

numres = tk.LabelFrame(win, text='Numerical Results', fg="brown", font="Verdana 8 bold")

numres.grid(column=1, row=1, sticky=tk.W + tk.N, padx=10, pady=10)

plttype = tk.LabelFrame(win, text='Select One for Single Plot', fg="brown", font="Verdana 8 bold")

plttype.grid(column=1, row=2, sticky=tk.W + tk.N, padx=10, pady=10)

res = tk.LabelFrame(win, text='Plot Area')

res.grid(column=2, row=1, rowspan=3)

a_input = tk.Label(simpara, text="RMS Value of Sine Wave/DC Voltage", font="Verdana 10 bold", fg="red")

a_input.grid(column=0, row=0, sticky=tk.W)

peamp = tk.Entry(simpara, width=10)

peamp.insert(0, "230")

peamp.grid(column=1, row=0)

a_supply = tk.Label(simpara, text="Supply Frequency in Hz", font="Verdana 10 bold", fg="red", justify=tk.LEFT)

a_supply.grid(column=0, row=1, sticky=tk.W)

freq = tk.Entry(simpara, width=10)

freq.insert(0, "50")

freq.grid(column=1, row=1)

a_time = tk.Label(simpara, text="Simulation Time in secs", font="Verdana 10 bold", fg="red")

a_time.grid(column=0, row=2, sticky=tk.W)

simtime = tk.Entry(simpara, width=10)

simtime.insert(0, "0.06")

simtime.grid(column=1, row=2)

a_freq = tk.Label(simpara, text="Switching Frequency in Hz", font="Verdana 10 bold", fg="magenta")

a_freq.grid(column=0, row=3, sticky=tk.W)

sf = tk.Entry(simpara, width=10)

sf.insert(0, "1000")

sf.grid(column=1, row=3)

a_mi = tk.Label(simpara, text="Modulation Index/Duty Cycle", font="Verdana 10 bold", fg="magenta", justify=tk.LEFT)

a_mi.grid(column=0, row=4, sticky=tk.W)

m = tk.Entry(simpara, width=10)

m.insert(0, "0.8")

m.grid(column=1, row=4)

a_fire = tk.Label(simpara, text="Firing Angle in Degrees", font="Verdana 10 bold", fg="magenta")

a_fire.grid(column=0, row=5, sticky=tk.W)

alpha = tk.Entry(simpara, width=10)

alpha.insert(0, "45")

alpha.grid(column=1, row=5)

tconvt = tk.Label(contypes, text="Type of Converter", font="Verdana 10 bold", fg="green").grid(column=0, row=1,
                                                                                               sticky=tk.W)

convtype = tk.IntVar()

rad1 = tk.Radiobutton(contypes, text='Half Wave Controlled Rectifier', value=1, variable=convtype,
                      font="Verdana 10 bold", fg="blue", )

rad2 = tk.Radiobutton(contypes, text='Full Wave Controlled Rectifier', value=2, variable=convtype,
                      font="Verdana 10 bold", fg="blue")

rad3 = tk.Radiobutton(contypes, text='Pulse Width Modulation Inverter', value=3, variable=convtype,
                      font="Verdana 10 bold", fg="blue")

rad4 = tk.Radiobutton(contypes, text='AC Voltage Regulator', value=4, variable=convtype, font="Verdana 10 bold",
                      fg="blue")

rad5 = tk.Radiobutton(contypes, text='DC Chopper', value=5, variable=convtype, font="Verdana 10 bold", fg="blue")

rad1.grid(column=0, row=2, sticky=tk.W)

rad2.grid(column=0, row=3, sticky=tk.W)

rad3.grid(column=0, row=4, sticky=tk.W)

rad4.grid(column=0, row=5, sticky=tk.W)

rad5.grid(column=0, row=6, sticky=tk.W)

convtype.set(1)

nphases = tk.Label(contypes, text="Number of Phases", font="Verdana 10 bold", fg="green").grid(column=0, row=7,
                                                                                               sticky=tk.W)

oneplt = tk.IntVar()

ipplt = tk.Radiobutton(plttype, text='Input Voltage', value=1, variable=oneplt, font="Verdana 10 bold", fg="blue", )

lvplt = tk.Radiobutton(plttype, text='Load Voltage', value=2, variable=oneplt, font="Verdana 10 bold", fg="blue", )

lcplt = tk.Radiobutton(plttype, text='Load Current', value=3, variable=oneplt, font="Verdana 10 bold", fg="blue", )

ilvplt = tk.Radiobutton(plttype, text='Input and Load Voltage', value=4, variable=oneplt, font="Verdana 10 bold",
                        fg="blue", )

ipplt.grid(column=0, row=1, sticky=tk.W)

lvplt.grid(column=0, row=2, sticky=tk.W)

lcplt.grid(column=0, row=3, sticky=tk.W)

ilvplt.grid(column=0, row=4, sticky=tk.W)

oneplt.set(1)

ivc = tk.IntVar()
lvc = tk.IntVar()
lcc = tk.IntVar()
ac = tk.IntVar()
rv = tk.IntVar()
cv = tk.IntVar()

tk.Checkbutton(plttype, text="Input Voltage", variable=ivc, font="Verdana 10 bold", fg="blue", ).grid(column=0, row=6,
                                                                                                      sticky=tk.W)
tk.Checkbutton(plttype, text="Load Voltage", variable=lvc, font="Verdana 10 bold", fg="blue", ).grid(column=0, row=7,
                                                                                                     sticky=tk.W)
tk.Checkbutton(plttype, text="Load Current", variable=lcc, font="Verdana 10 bold", fg="blue", ).grid(column=0, row=8,
                                                                                                     sticky=tk.W)
tk.Checkbutton(plttype, text="Angle(WT)", variable=ac, font="Verdana 10 bold", fg="blue", ).grid(column=0, row=9,
                                                                                                 sticky=tk.W)
tk.Checkbutton(plttype, text="Reference Waves", variable=rv, font="Verdana 10 bold", fg="blue", ).grid(column=0, row=10,
                                                                                                       sticky=tk.W)
tk.Checkbutton(plttype, text="Carrier Waves", variable=cv, font="Verdana 10 bold", fg="blue", ).grid(column=0, row=11,
                                                                                                     sticky=tk.W)


def sel():
    global l
    l = int(loadchoice.get())
    if l == 1:
        lvalue.config(state='disabled')
    else:
        lvalue.config(state='normal')


loadchoice = tk.IntVar()

phasechoice = tk.IntVar()

phase1 = tk.Radiobutton(contypes, text='Single Phase', value=1, variable=phasechoice, font="Verdana 10 bold", fg="blue")

phase3 = tk.Radiobutton(contypes, text='Three Phase', value=2, variable=phasechoice, font="Verdana 10 bold", fg="blue")

phase1.grid(column=0, row=8, sticky=tk.W)

phase3.grid(column=0, row=9, sticky=tk.W)

phasechoice.set(1)

ltype = tk.Label(contypes, text="Type of Load", font="Verdana 10 bold", fg="green").grid(column=0, row=10, sticky=tk.W)

load1 = tk.Radiobutton(contypes, text='R Load', value=1, variable=loadchoice, font="Verdana 10 bold", fg="blue",
                       command=sel)

load2 = tk.Radiobutton(contypes, text='RL Load', value=2, variable=loadchoice, font="Verdana 10 bold", fg="blue",
                       command=sel)

load1.grid(column=0, row=11, sticky=tk.W)

load2.grid(column=0, row=12, sticky=tk.W)

loadchoice.set(2)

r_input = tk.Label(contypes, text="R Value in Ohms", font="Verdana 10 bold", fg="red")

r_input.place(x=130, y=250)

rvalue = tk.Entry(contypes, width=10)

rvalue.insert(0, "100")

rvalue.grid(column=2, row=11)

l_input = tk.Label(contypes, text="L Value (H) ", font="Verdana 10 bold", fg="red")

l_input.place(x=130, y=275)

lvalue = tk.Entry(contypes, width=10)

lvalue.insert(0, "0.1")

lvalue.grid(column=2, row=12)

lvnv = tk.Label(numres, text="Average Load Voltage :", font="Verdana 10 bold", fg="red")

lvnv.grid(column=0, row=0, sticky=tk.W)

lcnv = tk.Label(numres, text="Average Load Current :", font="Verdana 10 bold", fg="red")

lcnv.grid(column=0, row=1, sticky=tk.W)

rlvnv = tk.Label(numres, text="RMS Load Voltage :", font="Verdana 10 bold", fg="red")

rlvnv.grid(column=0, row=2, sticky=tk.NE)

rlcnv = tk.Label(numres, text="RMS Load Current :", font="Verdana 10 bold", fg="red")

rlcnv.grid(column=0, row=3, sticky=tk.NE)

lpow = tk.Label(numres, text="Load Power :", font="Verdana 10 bold", fg="red")

lpow.grid(column=0, row=4, sticky=tk.NE)

ea = tk.Label(numres, text="Extinction Angle :", font="Verdana 10 bold", fg="red")

ea.grid(column=0, row=5, sticky=tk.NE)

evavl = tk.StringVar()
lvavl = tk.StringVar()
lcavl = tk.StringVar()
rmsvl = tk.StringVar()
rmscl = tk.StringVar()
lpowl = tk.StringVar()


def func(b):
    k = int(loadchoice.get())
    a1 = math.radians(float(alpha.get()))
    fr = int(freq.get())
    L1 = float(lvalue.get())
    R1 = float(rvalue.get())
    if k == 1:
        langle = (math.atan(2 * math.pi * fr * 0 / R1))
    else:
        langle = (math.atan(2 * math.pi * fr * L1 / R1))

    r1 = (a1 - b) / math.tan(langle)
    return math.sin(langle - b) - (math.sin(langle - a1) * math.exp(r1))


def simulate():
    # global s, l,pch, time, amplitude1, amplitude2, amplitude3, amplitude4, amplitude5, amplitude6, lcurrent, lvav, lcav, rmsv, rmsc,lpowv
    global s, l, pch, time, VR, VY, VB, LV, LC, ANG, R1, R2, C1, P1, P2, lcurrent, lvav, lcav, rmsv, rmsc, lpowv
    s = int(convtype.get())
    l = int(loadchoice.get())
    tt = float(simtime.get())
    V = float(peamp.get())
    f1 = int(freq.get())
    a = float(alpha.get())
    L = float(lvalue.get())
    R = float(rvalue.get())
    m1 = float(m.get())
    f2 = int(sf.get())
    pch = int(phasechoice.get())
    time = np.arange(0, tt, 0.00001);
    g = len(time)
    beta = (fsolve(func, 3.14))

    if l == 1:
        langle = (math.atan(2 * math.pi * f1 * 0 / R))
    else:
        langle = (math.atan(2 * math.pi * f1 * L / R))

    if s == 1 and pch == 1:
        P = 2 / f1
        ANG = (360 * (P - abs((time * 2 % P) - P))) / P
        VR = V * math.sqrt(2) * (np.sin(2 * math.pi * f1 * time))
        LV = V * math.sqrt(2) * (np.sin(2 * math.pi * f1 * time))
        Z = math.sqrt((R * R) + (2 * math.pi * f1 * 2 * math.pi * f1 * L * L))
        ad = math.radians(a)
        x1 = np.sin((2 * math.pi * f1 * time) - langle)
        x2 = math.sin(langle - ad)
        e1 = -((2 * math.pi * f1 * time) - ad) / math.tan(langle)
        x3 = np.exp(e1)
        LC = (V / Z) * ((x2 * x3) + x1)
        if l == 1:
            lvav = ((1 + cos(ad)) / (2 * math.pi)) * (V * sqrt(2))
            lcav = lvav / R
            rmsv = math.sqrt((math.pi - ad) + ((math.sin(2 * ad)) / 2)) * (V * sqrt(2) / (2 * math.sqrt(math.pi)))
            rmsc = rmsv / R
            lpowv = (rmsc * rmsc) * R
        else:
            lvav = (((np.cos(ad) - np.cos(beta))) / (2 * math.pi)) * (V * math.sqrt(2))
            lcav = lvav / R
            rmsv = math.sqrt(((beta - ad) - (((math.sin(2 * beta)) - (math.sin(2 * ad))) / 2)) * (
                        V * math.sqrt(2) * V * math.sqrt(2) / (4 * math.pi)))
            rmsc = rmsv / math.sqrt((R * R) + (2 * math.pi * f1 * 2 * math.pi * f1 * L * L))
            lpowv = (rmsc * rmsc) * R
        i = 0
        while i < g:
            if ANG[i] >= a and ANG[i] <= math.degrees(beta):
                LV[i] = LV[i]
                LC[i] = LC[i]
            else:
                LV[i] = 0
                LC[i] = 0
            i += 1

        lav = tk.Label(numres, textvariable=lvavl, font="Verdana 10 bold", fg="green")

        lav.grid(column=1, row=0, sticky=tk.W)

        lvavl.set(np.round(lvav, 2))

        lcurrav = tk.Label(numres, textvariable=lcavl, font="Verdana 10 bold", fg="green")

        lcurrav.grid(column=1, row=1, sticky=tk.W)

        lcavl.set(np.round(lcav, 2))

        rmsvv = tk.Label(numres, textvariable=rmsvl, font="Verdana 10 bold", fg="green")

        rmsvv.grid(column=1, row=2, sticky=tk.W)

        rmsvl.set(np.round(rmsv, 2))

        rmscv = tk.Label(numres, textvariable=rmscl, font="Verdana 10 bold", fg="green")

        rmscv.grid(column=1, row=3, sticky=tk.W)

        rmscl.set(np.round(rmsc, 2))

        lpow1 = tk.Label(numres, textvariable=lpowl, font="Verdana 10 bold", fg="green")

        lpow1.grid(column=1, row=4, sticky=tk.W)

        lpowl.set(np.round(lpowv, 0))

        eav = tk.Label(numres, textvariable=evavl, font="Verdana 10 bold", fg="green")

        eav.grid(column=1, row=5, sticky=tk.W)

        evavl.set(np.round(math.degrees(beta), 0))

        tk.messagebox.showinfo("Simulation", "Simulation Completed")

    elif s == 1 and pch == 2:
        P = 2 / f1
        ANG = (360 * (P - abs((time * 2 % P) - P))) / P
        VR = V * math.sqrt(2) * (np.sin(2 * math.pi * f1 * time))
        VY = V * sqrt(2) * (np.sin(2 * math.pi * f1 * time - (2 * math.pi / 3)))
        VB = V * sqrt(2) * (np.sin(2 * math.pi * f1 * time + (2 * math.pi / 3)))
        LV = 0 * (np.sin(2 * math.pi * f1 * time + (2 * math.pi / 3)))

        if (a + 150) > 180:
            a1 = 180
        else:
            a1 = a + 150

        if (a + 270) > 300 and (a + 270) <= 360:
            a2 = 300
            a3 = a + 270
            a5 = 0
            a6 = (a2 - (a + 150)) - (360 - (a + 270))
        elif (a + 270) <= 300:
            a2 = a + 270
            a3 = a2
            a5 = 0
            a6 = a + 30
        else:
            a2 = 300
            a3 = 360
            a5 = (a + 270) - 360
            a6 = (a2 - (a + 150)) - (360 - (a + 270))
        i = 0
        while i < g:
            if ANG[i] >= a + 30 and ANG[i] <= a1:
                LV[i] = VR[i]
            elif ANG[i] >= a + 150 and ANG[i] <= a2:
                LV[i] = VY[i]
            elif ANG[i] >= a3 and ANG[i] <= 360:
                LV[i] = VB[i]
            elif ANG[i] >= a5 and ANG[i] <= a6:
                LV[i] = VB[i]
            else:
                LV[i] = 0

            i += 1

        tk.messagebox.showinfo("Simulation", "Simulation Completed")
    elif s == 2 and pch == 1:
        beta1 = math.degrees(beta)
        beta2 = beta1 - 180
        if beta1 >= (180 + a):
            beta2 = a
        else:
            beta2 = beta2

        print(math.degrees(beta))
        print(beta2)
        P = 2 / f1
        ANG = (360 * (P - abs((time * 2 % P) - P))) / P
        VR = V * sqrt(2) * np.sin(2 * math.pi * f1 * time)
        LV = V * sqrt(2) * (np.sin(2 * math.pi * f1 * time))
        Z = math.sqrt((R * R) + (2 * math.pi * f1 * 2 * math.pi * f1 * L * L))
        ad = math.radians(a)
        x1 = sin((2 * math.pi * f1 * time) - langle)
        x2 = math.sin(langle - ad)
        e1 = -((2 * math.pi * f1 * time) - ad) / math.tan(langle)
        x3 = np.exp(e1)
        LC = (V / Z) * ((x2 * x3) + x1)
        if l == 1:
            lvav = ((1 + cos(ad)) / (math.pi)) * (V * sqrt(2))
            lcav = lvav / R
            rmsv = (math.sqrt(0.5 - (ad / (2 * math.pi)) + (math.sin(2 * ad) / (4 * math.pi)))) * (V * sqrt(2))
            rmsc = rmsv / R
            lpowv = (rmsc * rmsc) * R
        else:
            lvav = (((cos(ad) - cos(beta))) / (math.pi)) * (V * sqrt(2))
            lcav = lvav / R
            rmsv = math.sqrt(((beta - ad) - (((math.sin(2 * beta)) - (math.sin(2 * ad))) / 2)) * (
                        V * sqrt(2) * V * sqrt(2) / (2 * math.pi)))
            rmsc = rmsv / math.sqrt((R * R) + (2 * math.pi * f1 * 2 * math.pi * f1 * L * L))
            lpowv = (rmsc * rmsc) * R
        i = 0
        while i < g:
            if ANG[i] >= 0 and ANG[i] <= beta2:
                LV[i] = -LV[i]
                LC[i] = -LC[i]
            elif ANG[i] >= a and ANG[i] <= 180 + beta2:
                LV[i] = LV[i]
                LC[i] = LC[i]
            elif ANG[i] >= (180 + a) and ANG[i] <= 360 + beta2:
                LV[i] = -LV[i]
                LC[i] = -LC[i]
            else:
                LV[i] = 0
                LC[i] = 0
            i += 1

        lav = tk.Label(numres, textvariable=lvavl, font="Verdana 10 bold", fg="green")

        lav.grid(column=1, row=0, sticky=tk.W)

        lvavl.set(np.round(lvav, 2))

        lcurrav = tk.Label(numres, textvariable=lcavl, font="Verdana 10 bold", fg="green")

        lcurrav.grid(column=1, row=1, sticky=tk.W)

        lcavl.set(np.round(lcav, 2))

        rmsvv = tk.Label(numres, textvariable=rmsvl, font="Verdana 10 bold", fg="green")

        rmsvv.grid(column=1, row=2, sticky=tk.W)

        rmsvl.set(np.round(rmsv, 2))

        rmscv = tk.Label(numres, textvariable=rmscl, font="Verdana 10 bold", fg="green")

        rmscv.grid(column=1, row=3, sticky=tk.W)

        rmscl.set(np.round(rmsc, 2))

        lpow1 = tk.Label(numres, textvariable=lpowl, font="Verdana 10 bold", fg="green")

        lpow1.grid(column=1, row=4, sticky=tk.W)

        lpowl.set(np.round(lpowv, 0))

        eav = tk.Label(numres, textvariable=evavl, font="Verdana 10 bold", fg="green")

        eav.grid(column=1, row=5, sticky=tk.W)

        evavl.set(np.round(math.degrees(beta), 0))

        tk.messagebox.showinfo("Simulation", "Simulation Completed")

    elif s == 3:
        P = 1 / f2
        C1 = (1 * (P - abs((time * 2 % P) - P))) / P  # C1
        R1 = m1 * np.sin(2 * math.pi * f1 * time)  # R1
        R2 = m1 * np.sin(2 * math.pi * f1 * time - math.pi)  # R2
        P1 = np.sin(2 * math.pi * f1 * time - math.pi)  # P1
        P2 = np.sin(2 * math.pi * f1 * time - math.pi)  # P2
        LV = np.sin(2 * math.pi * f1 * time - math.pi)  # LV
        VR = V
        i = 0
        while i < g:
            if C1[i] > R1[i]:
                P1[i] = 0
                LV[i] = 0
            elif C1[i] < R1[i]:
                P1[i] = 1
                LV[i] = V
            if C1[i] > R2[i]:
                P2[i] = 0
            elif C1[i] < R2[i]:
                P2[i] = 1

            if C1[i] > R2[i] and R2[i] > 0:
                LV[i] = 0
            elif C1[i] < R2[i] and R2[i] > 0:
                LV[i] = -V

            i += 1
        tk.messagebox.showinfo("Simulation", "Simulation Completed")
    elif s == 4:
        beta1 = math.degrees(beta)
        beta2 = beta1 - 180
        if beta1 >= (180 + a):
            beta2 = a
        else:
            beta2 = beta2

        print(math.degrees(beta))
        print(beta2)
        P = 2 / f1
        ANG = (360 * (P - abs((time * 2 % P) - P))) / P  # ANG
        VR = V * sqrt(2) * np.sin(2 * math.pi * f1 * time)  # VR
        LV = V * sqrt(2) * (np.sin(2 * math.pi * f1 * time))  # LV
        Z = math.sqrt((R * R) + (2 * math.pi * f1 * 2 * math.pi * f1 * L * L))
        ad = math.radians(a)
        x1 = sin((2 * math.pi * f1 * time) - langle)
        x2 = math.sin(langle - ad)
        e1 = -((2 * math.pi * f1 * time) - ad) / math.tan(langle)
        x3 = np.exp(e1)
        LC = (V / Z) * ((x2 * x3) + x1)

        i = 0
        while i < g:
            if ANG[i] >= 0 and ANG[i] <= beta2:
                LV[i] = LV[i]
                LC[i] = LC[i]
            elif ANG[i] >= a and ANG[i] <= 180 + beta2:
                LV[i] = LV[i]
                LC[i] = LC[i]
            elif ANG[i] >= (180 + a) and ANG[i] <= 360 + beta2:
                LV[i] = LV[i]
                LC[i] = LC[i]
            else:
                LV[i] = 0
                LC[i] = 0
            i += 1

        tk.messagebox.showinfo("Simulation", "Simulation Completed")

    elif s == 5:
        P = 1 / f2
        C1 = (1 * (P - abs((time * 2 % P) - P))) / P  # C1

        LV = np.sin(2 * math.pi * f1 * time - math.pi)  # LV

        i = 0
        while i < g:
            if C1[i] > m1:
                LV[i] = 0
            else:
                LV[i] = V

            i += 1
        tk.messagebox.showinfo("Simulation", "Simulation Completed")


def plt1():
    if s == 1 and pch == 1:
        fig = Figure(figsize=(6, 5), dpi=100)
        a1 = fig.add_subplot(211)
        a1.plot(time, VR, label="INPUT VOLTAGE")
        a1.plot(time, LV, label="LOAD VOLTAGE")
        a1.set_title('SINGLE PHASE CONTROLLED RECTIFIER-HALF CONTROLLER')
        a1.set_ylabel('Volts')
        a1.legend(bbox_to_anchor=(1.05, 1), loc='upper right', borderaxespad=0.)
        a1.grid(b=True, which='major', color='#666666', linestyle=':')

        a3 = fig.add_subplot(212)
        a3.plot(time, LC, label="Load Current")
        a3.set_xlabel('Time in ms')
        a3.set_ylabel('Current in Amps')
        a3.grid(b=True, which='major', color='#666666', linestyle=':')

        canvas = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(fig, master=win)
        canvas.draw()
        canvas.get_tk_widget().grid(row=1, column=2, rowspan=3)
        toolbarFrame = tk.Frame(win)
        toolbarFrame.grid(row=5, column=2)
        toolbar = matplotlib.backends.backend_tkagg.NavigationToolbar2Tk(canvas, toolbarFrame)
        toolbar.update()

    elif s == 1 and pch == 2:
        fig5 = Figure(figsize=(6, 5), dpi=100)
        a1 = fig5.add_subplot(311)
        a1.plot(time, VR, label="R Phase")
        a1.plot(time, VY, label="Y Phase")
        a1.plot(time, VB, label="B Phase")
        a1.plot(time, LV, label="Load Voltage")
        a1.set_title('THREE PHASE CONTROLLED RECTIFIER-HALF CONTROLLER')
        a1.set_ylabel('Volts')
        # a1.legend(bbox_to_anchor=(1.05, 1), loc='upper right', borderaxespad=0.)
        a1.grid(b=True, which='major', color='#666666', linestyle=':')

        a3 = fig5.add_subplot(312)
        a3.plot(time, LV, label="LOAD VOLTAGE")
        a3.set_xlabel('Time in ms')
        a3.set_ylabel('LOAD VOLTAGE')
        a3.grid(b=True, which='major', color='#666666', linestyle=':')

        a4 = fig5.add_subplot(313)
        a4.plot(time, ANG, label="ANGLE")
        a4.set_xlabel('Time in ms')
        a4.set_ylabel('ANGLE')
        a4.grid(b=True, which='major', color='#666666', linestyle=':')

        canvas = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(fig5, master=win)
        canvas.draw()
        canvas.get_tk_widget().grid(row=1, column=2, rowspan=3)
        toolbarFrame = tk.Frame(win)
        toolbarFrame.grid(row=5, column=2)
        toolbar = matplotlib.backends.backend_tkagg.NavigationToolbar2Tk(canvas, toolbarFrame)
        toolbar.update()

    elif s == 2 and pch == 1:
        fig1 = Figure(figsize=(6, 5), dpi=100)
        # t = np.arange(0, 3, .01)
        a1 = fig1.add_subplot(211)
        a1.plot(time, VR, label="INPUT VOLTAGE")
        a1.plot(time, LV, label="LOAD VOLTAGE")
        a1.set_title('SINGLE PHASE CONTROLLED RECTIFIER-FULL CONTROLLER ')
        a1.set_ylabel('Volts')
        a1.legend(bbox_to_anchor=(1.05, 1), loc='upper right', borderaxespad=0.)
        a1.grid(b=True, which='major', color='#666666', linestyle=':')

        a3 = fig1.add_subplot(212)
        a3.plot(time, LC, label="Load Current")
        a3.set_xlabel('Time in ms')
        a3.set_ylabel('Current in Amps')
        a3.grid(b=True, which='major', color='#666666', linestyle=':')

        canvas = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(fig1, master=win)  # A tk.DrawingArea.
        canvas.draw()
        canvas.get_tk_widget().grid(row=1, column=2, rowspan=3)
        toolbarFrame = tk.Frame(win)
        toolbarFrame.grid(row=5, column=2)
        toolbar = matplotlib.backends.backend_tkagg.NavigationToolbar2Tk(canvas, toolbarFrame)
        toolbar.update()
    elif s == 3:
        fig2 = Figure(figsize=(6, 5), dpi=100)
        a1 = fig2.add_subplot(311)
        a1.plot(time, C1, label="Carrier")
        a1.plot(time, R1, label="Reference 1")
        a1.plot(time, R2, label="Reference 2")
        a1.set_title('PWM GENERATION FOR SINGLE PHASE INVERTER')
        a1.set_ylabel('Volts')
        a1.legend(bbox_to_anchor=(1.05, 1), loc='upper right', borderaxespad=0.)
        a1.grid(b=True, which='major', color='#666666', linestyle=':')

        a3 = fig2.add_subplot(312)
        a3.plot(time, P1, label="PWM FOR SWITCHES S1 AND S2")
        a3.set_ylabel('Volts')
        a3.grid(b=True, which='major', color='#666666', linestyle=':')
        a3.legend(bbox_to_anchor=(1.05, 1), loc='upper right', borderaxespad=0.)


        a4 = fig2.add_subplot(313)
        a4.plot(time, LV, label="LOAD VOLTAGE")
        a4.legend(bbox_to_anchor=(1.05, 1), loc='upper right', borderaxespad=0.)
        a4.grid(b=True, which='major', color='#666666', linestyle=':')
        a4.set_ylabel('Volts')
        a4.set_xlabel('Time in ms')

        canvas = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(fig2, master=win)  # A tk.DrawingArea.
        canvas.draw()
        canvas.get_tk_widget().grid(row=1, column=2, rowspan=3)
        toolbarFrame = tk.Frame(win)
        toolbarFrame.grid(row=5, column=2)
        toolbar = matplotlib.backends.backend_tkagg.NavigationToolbar2Tk(canvas, toolbarFrame)
        toolbar.update()
    elif s == 4:
        fig3 = Figure(figsize=(6, 5), dpi=100)
        # t = np.arange(0, 3, .01)
        a1 = fig3.add_subplot(211)
        a1.plot(time, VR, label="INPUT VOLTAGE")
        a1.plot(time, LV, label="LOAD VOLTAGE")
        a1.set_title('SINGLE PHASE AC VOLTAGE REGULATOR')
        a1.set_ylabel('Volts')
        a1.legend(bbox_to_anchor=(1.05, 1), loc='upper right', borderaxespad=0.)
        a1.grid(b=True, which='major', color='#666666', linestyle=':')

        a3 = fig3.add_subplot(212)
        a3.plot(time, LC, label="Load Current")
        a3.set_xlabel('Time in ms')
        a3.set_ylabel('Current in Amps')
        a3.grid(b=True, which='major', color='#666666', linestyle=':')

        canvas = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(fig3, master=win)  # A tk.DrawingArea.
        canvas.draw()
        canvas.get_tk_widget().grid(row=1, column=2, rowspan=3)
        toolbarFrame = tk.Frame(win)
        toolbarFrame.grid(row=5, column=2)
        toolbar = matplotlib.backends.backend_tkagg.NavigationToolbar2Tk(canvas, toolbarFrame)
        toolbar.update()
    elif s == 5:
        fig4 = Figure(figsize=(6, 5), dpi=100)
        a1 = fig4.add_subplot(211)
        a1.plot(time, LV, label="Load Voltage")
        a1.set_title('LOAD VOLTAGE DC-DC BUCK CONVERTER')
        a1.set_ylabel('Volts')
        a1.grid(b=True, which='major', color='#666666', linestyle=':')

        a3 = fig4.add_subplot(212)
        a3.plot(time, C1, label="CARRIER SIGNAL")
        a3.set_ylabel('Volts')
        a3.grid(b=True, which='major', color='#666666', linestyle=':')
        a3.legend(bbox_to_anchor=(1.05, 1), loc='upper right', borderaxespad=0.)

        canvas = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(fig4, master=win)  # A tk.DrawingArea.
        canvas.draw()
        canvas.get_tk_widget().grid(row=1, column=2, rowspan=3)
        toolbarFrame = tk.Frame(win)
        toolbarFrame.grid(row=5, column=2)
        toolbar = matplotlib.backends.backend_tkagg.NavigationToolbar2Tk(canvas, toolbarFrame)
        toolbar.update()


def fplt():
    FP = int(oneplt.get())
    if FP == 1:
        fig = Figure(figsize=(6, 5), dpi=100)
        p1 = fig.add_subplot(111)
        p1.plot(time, VR, label="INPUT VOLTAGE")
        p1.set_title('SINGLE PHASE CONTROLLED RECTIFIER-HALF CONTROLLER(INPUT VOLTAGE)')
        p1.set_ylabel('Volts')
        p1.legend(bbox_to_anchor=(1.05, 1), loc='upper right', borderaxespad=0.)
        p1.grid(b=True, which='major', color='#666666', linestyle=':')

        canvas = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(fig, master=win)
        canvas.draw()
        canvas.get_tk_widget().grid(row=1, column=2, rowspan=3)
        toolbarFrame = tk.Frame(win)
        toolbarFrame.grid(row=5, column=2)
        toolbar = matplotlib.backends.backend_tkagg.NavigationToolbar2Tk(canvas, toolbarFrame)
        toolbar.update()
    elif FP == 2:
        fig = Figure(figsize=(6, 5), dpi=100)
        p2 = fig.add_subplot(111)
        p2.plot(time, LV, label="LOAD VOLTAGE")
        p2.set_title('SINGLE PHASE CONTROLLED RECTIFIER-HALF CONTROLLER(LOAD VOLTAGE)')
        p2.set_ylabel('Volts')
        p2.legend(bbox_to_anchor=(1.05, 1), loc='upper right', borderaxespad=0.)
        p2.grid(b=True, which='major', color='#666666', linestyle=':')

        canvas = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(fig, master=win)
        canvas.draw()
        canvas.get_tk_widget().grid(row=1, column=2, rowspan=3)
        toolbarFrame = tk.Frame(win)
        toolbarFrame.grid(row=5, column=2)
        toolbar = matplotlib.backends.backend_tkagg.NavigationToolbar2Tk(canvas, toolbarFrame)
        toolbar.update()
    elif FP == 3:
        fig = Figure(figsize=(6, 5), dpi=100)
        p3 = fig.add_subplot(111)
        p3.plot(time, LC, label="LOAD CURRENT")
        p3.set_title('SINGLE PHASE CONTROLLED RECTIFIER-HALF CONTROLLER(LOAD CURRENT)')
        p3.set_ylabel('Amps')
        p3.legend(bbox_to_anchor=(1.05, 1), loc='upper right', borderaxespad=0.)
        p3.grid(b=True, which='major', color='#666666', linestyle=':')

        canvas = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(fig, master=win)
        canvas.draw()
        canvas.get_tk_widget().grid(row=1, column=2, rowspan=3)
        toolbarFrame = tk.Frame(win)
        toolbarFrame.grid(row=5, column=2)
        toolbar = matplotlib.backends.backend_tkagg.NavigationToolbar2Tk(canvas, toolbarFrame)
        toolbar.update()

    elif FP == 4:
        fig = Figure(figsize=(6, 5), dpi=100)
        p4 = fig.add_subplot(111)
        p4.plot(time, VR, label="INPUT VOLTAGE")
        p4.plot(time, LV, label="LOAD VOLTAGE")
        p4.set_title('SINGLE PHASE CONTROLLED RECTIFIER-HALF CONTROLLER')
        p4.set_ylabel('Volts')
        p4.legend(bbox_to_anchor=(1.05, 1), loc='upper right', borderaxespad=0.)
        p4.grid(b=True, which='major', color='#666666', linestyle=':')

        canvas = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(fig, master=win)
        canvas.draw()
        canvas.get_tk_widget().grid(row=1, column=2, rowspan=3)
        toolbarFrame = tk.Frame(win)
        toolbarFrame.grid(row=5, column=2)
        toolbar = matplotlib.backends.backend_tkagg.NavigationToolbar2Tk(canvas, toolbarFrame)
        toolbar.update()


def Splt1():
    global fig2
    fig2 = Figure(figsize=(6, 5), dpi=100)

    s1 = fig2.add_subplot(311)

    if (ivc.get() == 1):
        s1.plot(time, VR, label="INPUT VOLTAGE")
        s1.set_ylabel('Volts')
    if (lvc.get() == 1):
        s1.plot(time, LV, label="LOAD VOLTAGE")
        s1.set_ylabel('Volts')
    if (lcc.get() == 1):
        s1.plot(time, LC, label="LOAD CURRENT")
        s1.set_ylabel('Amps')
    if (ac.get() == 1):
        s1.plot(time, ANG, label="ANGLE")
        s1.set_ylabel('Degrees')
    if (rv.get() == 1):
        s1.plot(time, R1, label="Reference 1")
        s1.plot(time, R2, label="Reference 2")
        s1.set_ylabel('Volts')
    if (cv.get() == 1):
        s1.plot(time, C1, label="Carrier")
        s1.set_ylabel('Volts')

    s1.set_title('PWM GENERATION FOR SINGLE PHASE INVERTER')

    s1.legend(bbox_to_anchor=(1.05, 1), loc='upper right', borderaxespad=0.)
    s1.grid(b=True, which='major', color='#666666', linestyle=':')

    canvas = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(fig2, master=win)
    canvas.draw()
    canvas.get_tk_widget().grid(row=1, column=2, rowspan=3)
    toolbarFrame = tk.Frame(win)
    toolbarFrame.grid(row=5, column=2)
    toolbar = matplotlib.backends.backend_tkagg.NavigationToolbar2Tk(canvas, toolbarFrame)
    toolbar.update()


def Splt2():
    # fig2 = Figure(figsize=(6, 5), dpi=100)

    s2 = fig2.add_subplot(312)

    if (ivc.get() == 1):
        s2.plot(time, VR, label="INPUT VOLTAGE")
        s2.set_ylabel('Volts')
    if (lvc.get() == 1):
        s2.plot(time, LV, label="LOAD VOLTAGE")
        s2.set_ylabel('Volts')
    if (lcc.get() == 1):
        s2.plot(time, LC, label="LOAD CURRENT")
        s2.set_ylabel('Amps')
    if (ac.get() == 1):
        s2.plot(time, ANG, label="ANGLE")
        s2.set_ylabel('Degrees')
    if (rv.get() == 1):
        s2.plot(time, R1, label="Reference 1")
        s2.plot(time, R2, label="Reference 2")
        s2.set_ylabel('Volts')
    if (cv.get() == 1):
        s2.plot(time, C1, label="Carrier")
        s2.set_ylabel('Volts')

    # s2.set_title('PWM GENERATION FOR SINGLE PHASE INVERTER')

    s2.legend(bbox_to_anchor=(1.05, 1), loc='upper right', borderaxespad=0.)
    s2.grid(b=True, which='major', color='#666666', linestyle=':')

    canvas = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(fig2, master=win)
    canvas.draw()
    canvas.get_tk_widget().grid(row=1, column=2, rowspan=3)
    toolbarFrame = tk.Frame(win)
    toolbarFrame.grid(row=5, column=2)
    toolbar = matplotlib.backends.backend_tkagg.NavigationToolbar2Tk(canvas, toolbarFrame)
    toolbar.update()


def Splt3():
    # fig2 = Figure(figsize=(6, 5), dpi=100)

    s3 = fig2.add_subplot(313)

    if (ivc.get() == 1):
        s3.plot(time, VR, label="INPUT VOLTAGE")
        s3.set_ylabel('Volts')
    if (lvc.get() == 1):
        s3.plot(time, LV, label="LOAD VOLTAGE")
        s3.set_ylabel('Volts')
    if (lcc.get() == 1):
        s3.plot(time, LC, label="LOAD CURRENT")
        s3.set_ylabel('Amps')
    if (ac.get() == 1):
        s3.plot(time, ANG, label="ANGLE")
        s3.set_ylabel('Degrees')
    if (rv.get() == 1):
        s3.plot(time, R1, label="Reference 1")
        s3.plot(time, R2, label="Reference 2")
        s3.set_ylabel('Volts')
    if (cv.get() == 1):
        s3.plot(time, C1, label="Carrier")
        s3.set_ylabel('Volts')

    # s3.set_title('PWM GENERATION FOR SINGLE PHASE INVERTER')

    s3.legend(bbox_to_anchor=(1.05, 1), loc='upper right', borderaxespad=0.)
    s3.grid(b=True, which='major', color='#666666', linestyle=':')

    canvas = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(fig2, master=win)
    canvas.draw()
    canvas.get_tk_widget().grid(row=1, column=2, rowspan=3)
    toolbarFrame = tk.Frame(win)
    toolbarFrame.grid(row=5, column=2)
    toolbar = matplotlib.backends.backend_tkagg.NavigationToolbar2Tk(canvas, toolbarFrame)
    toolbar.update()

btn = tk.Button(win, text="SIMULATE", command=simulate, font="Verdana 10 bold")

btn.place(x=25, y=550, anchor="w")

btn1 = tk.Button(win, text="FULL PLOT", command=plt1, font="Verdana 10 bold")

btn1.place(x=150, y=550, anchor="w")

btn2 = tk.Button(win, text="S.PLOT 1", command=Splt1, font="Verdana 10 bold")

btn2.place(x=275, y=550, anchor="w")

btn3 = tk.Button(win, text="S.PLOT 2", command=Splt2, font="Verdana 10 bold")

btn3.place(x=375, y=550, anchor="w")

btn4 = tk.Button(win, text="S.PLOT 3", command=Splt3, font="Verdana 10 bold")

btn4.place(x=475, y=550, anchor="w")

fullplt = tk.Button(plttype, text="SINGLE PLOT", command=fplt, font="Verdana 10 bold")

fullplt.grid(column=0, row=5, sticky=tk.N, padx=10, pady=5)

win.mainloop()