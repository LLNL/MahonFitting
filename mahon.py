"""
Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at the Lawrence Livermore National Laboratory.
Written by Reto Trappitsch, trappitsch1@llnl.gov

LLNL-CODE-745740 All rights reserved. This file is part of MahonFitting v1.0

Please also read this link - Our Disclaimer (https://github.com/LLNL/MahonFitting/blob/master/DISCLAIMER) and
GNU General Public License (https://github.com/LLNL/MahonFitting/blob/master/LICENSE).

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
License (as published by the Free Software Foundation) version 2, dated June 1991.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the
Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
"""

# calculate linear regression w/ error calculation according Mahon 1996 (New York regression)
try:
    import Tkinter as tk
    import tkFileDialog as filedialog
    import tkMessageBox as messagebox
except ImportError:
    import tkinter as tk
    from tkinter import filedialog
    from tkinter import messagebox
import numpy as np

import matplotlib
matplotlib.use('TkAgg')

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure


class Mahon:
    def __init__(self, master):
        self.master = master
        master.title("Mahon linear regression")

        # Geometry and stuff
        rows, cols = 9, 3
        BW, BH, pad = 120, 30, 15
        winsize = ((1.5 + cols / 2) * pad + cols * BW,
                   (1.5 + rows / 2) * pad + rows * BH,
                   10, 10)  # width, height, xoffset, yoffset
        self.rows, self.cols = rows, cols
        self.BW, self.BH, self.pad = BW, BH, pad
        self.winsize = winsize

        # set master size
        master.geometry("%dx%d+%d+%d" % winsize)

        # fixed intercept
        self.afx = None

        # some variables to have for later
        self.slope = None
        self.yinter = None
        self.slopeunc = None
        self.yinterunc = None
        self.xinterunc = None
        self.mswd = None
        self.fname = None

        # some variables to remember
        self.xbar = None
        self.ybar = None

        # create variables for holding the data
        self.xdat, self.xunc = None, None
        self.ydat, self.yunc = None, None
        self.p = None
        # create variables for holding the backup
        self._xdat, self._xunc = None, None
        self._ydat, self._yunc = None, None
        self._p = None

        # make load general text and load csv file
        self.loadtxt_button = tk.Button(master, text='Load txt', command=self.loadtxt)
        self.loadtxt_button.place(x=self.x(1), y=self.y(1), width=BW, height=BH)
        self.loadcsv_button = tk.Button(master, text='Load csv', command=self.loadcsv)
        self.loadcsv_button.place(x=self.x(1), y=self.y(2), width=BW, height=BH)

        # make the GUI
        self.calculate_button = tk.Button(master, text='Calculate', command=self.calculate)
        self.calculate_button.place(x=self.x(2), y=self.y(1), width=BW, height=BH)
        self.plot_button = tk.Button(master, text='Plot', command=self.plotdata)
        self.plot_button.place(x=self.x(3), y=self.y(1), width=BW, height=BH)

        # correlation value
        self.p_corr_var = tk.IntVar()
        self.p_corr_box = tk.Checkbutton(master, text='Calculate with correlation values',
                                         variable=self.p_corr_var)
        self.p_corr_box.place(x=self.x(2), y=self.y(2), width=2*BW+pad, height=BH, anchor=tk.NW)

        # how many sigma uncertainties are in the file?
        self.sig_label = tk.Label(master, text='Uncertainties:', justify=tk.LEFT)
        self.sig_label.place(x=self.x(1), y=self.y(3), width=BW, height=BH, anchor=tk.NW)
        self.sigvar = tk.IntVar()
        self.sig_rad1 = tk.Radiobutton(master, text='1 sigma', variable=self.sigvar, value=1)
        self.sig_rad1.place(x=self.x(2), y=self.y(3), width=BW, height=BH, anchor=tk.NW)
        self.sig_rad2 = tk.Radiobutton(master, text='2 sigma', variable=self.sigvar, value=2)
        self.sig_rad2.place(x=self.x(3), y=self.y(3), width=BW, height=BH, anchor=tk.NW)
        self.sig_rad1.select()

        # Force intercept through zero toggle
        self.force_int_var = tk.IntVar()
        self.force_int_label = tk.Label(master, text='Intercept:', justify=tk.LEFT)
        self.force_int_label.place(x=self.x(1), y=self.y(4), width=BW, height=BH)
        self.force_int_box = tk.Checkbutton(master, text='Fix x=0 to y =',
                                            variable=self.force_int_var, command=self.enablefixed)
        self.force_int_box.place(x=self.x(2), y=self.y(4), width=1.5 * BW + pad, height=BH)
        self.force_int_entry = tk.Entry(master, state=tk.DISABLED)
        self.force_int_entry.place(x=self.x(3)+BW/2., y=self.y(4), width=BW/2., height=BH)

        # make the slope and intercept labels
        self.a_label = tk.Label(master, text='slope:', justify=tk.LEFT)
        self.a_label.place(x=self.x(1), y = self.y(5), width=BW, height=BH, anchor=tk.NW)
        self.b_label = tk.Label(master, text='y-intercept:', justify=tk.LEFT)
        self.b_label.place(x=self.x(1), y=self.y(6), width=BW, height=BH, anchor=tk.NW)
        self.mswd_label = tk.Label(master, text='x-intercept:', justify=tk.LEFT)
        self.mswd_label.place(x=self.x(1), y=self.y(7), width=BW, height=BH, anchor=tk.NW)
        self.mswd_label = tk.Label(master, text='MSWD:', justify=tk.LEFT)
        self.mswd_label.place(x=self.x(1), y=self.y(8), width=BW, height=BH, anchor=tk.NW)
        # text boxes for output
        self.a_text = tk.Entry(master)
        self.a_text.place(x=self.x(2), y=self.y(5), width=2*BW+pad/2, height=BH)
        self.b_text = tk.Entry(master)
        self.b_text.place(x=self.x(2), y=self.y(6), width=2*BW+pad/2, height=BH)
        self.x_text = tk.Entry(master)
        self.x_text.place(x=self.x(2), y=self.y(7), width=2*BW+pad/2, height=BH)
        self.mswd_text = tk.Entry(master)
        self.mswd_text.place(x=self.x(2), y=self.y(8), width=2*BW+pad/2, height=BH)

        # Help button
        self.help_button = tk.Button(master, text='Help', command=self.get_help)
        self.help_button.place(x=self.x(2), y=self.y(9), width=BW, height=BH)
        # Quit button
        self.quit_button = tk.Button(master, text='Quit', command=self.quit_app)
        self.quit_button.place(x=self.x(3), y=self.y(9), width=BW, height=BH)

    def x(self, col):
        # Calculates x coordinate of gui element
        return self.pad + 0.5 * (col - 1) * self.pad + (col - 1) * self.BW

    def y(self, row):
        # Calculates y coordinate of gui element
        return self.pad + 0.5 * (row - 1) * self.pad + (row - 1) * self.BH

    def enablefixed(self):
        if self.force_int_var.get():
            self.force_int_entry.config(state='normal')
        else:
            self.force_int_entry.config(state='disabled')

    def loadfile(self, fname, splitter=None):
        f = open(fname, 'r')
        datain = []
        for line in f:
            datain.append(line)
        f.close()

        # check for mac file format as exported from excel
        if len(datain) == 1:
            tmp = []
            datsplit = datain[0].split('\r')
            for it in range(len(datsplit)):
                tmp.append(datsplit[it])
            datain = tmp

        # put filename into self
        self.fname = fname

        # create data
        data = []
        xdat, ydat, xunc, yunc, p = [], [], [], [], []

        for it in range(len(datain)):
            # strip newline characters
            datain[it].strip()
            datain[it].rstrip()
            if splitter is None:
                data.append(datain[it].split())
            else:
                data.append(datain[it].split(splitter))

        for it in range(len(data)):
            if len(data[it]) == 4:
                # catch header issues
                if data[it][0][0] == '#':
                    pass
                else:
                    try:
                        xdat.append(float(data[it][0]))
                        xunc.append(float(data[it][1]))
                        ydat.append(float(data[it][2]))
                        yunc.append(float(data[it][3]))
                        # append a zero for p
                        p.append(0.)
                    except ValueError:
                        # probably a header, not a number for sure
                        continue
            if len(data[it]) > 4:
                # catch header issues
                if data[it][0][0] == '#':
                    pass
                else:
                    self.p_corr_box.select()
                    try:
                        xdat.append(float(data[it][0]))
                        xunc.append(float(data[it][1]))
                        ydat.append(float(data[it][2]))
                        yunc.append(float(data[it][3]))
                        if data[it][4] != '':
                            p.append(float(data[it][4]))
                        else:   # if no correlation is given but there is still more than 4 columns for some reason
                            p.append(0.)
                    except ValueError:
                        # probably a header, not a number for sure
                        continue

        # create backup
        self._xdat, self._ydat, self._xunc, self._yunc = np.array(xdat), np.array(ydat), np.array(xunc), np.array(yunc)
        self._p = p

    def loadtxt(self):
        options = dict()
        options['defaultextension'] = '.txt'
        options['filetypes'] = [('txt files', '.txt'), ('all files', '.*')]
        fname = filedialog.askopenfilename(**options)
        if fname == '':
            return
        self.loadfile(fname)

    def loadcsv(self):
        options = dict()
        options['defaultextension'] = '.csv'
        options['filetypes'] = [('csv files', '.csv'), ('all files', '.*')]
        fname = filedialog.askopenfilename(**options)
        if fname == '':
            return
        self.loadfile(fname, splitter=',')

    def calculate(self):
        if self._xdat is None or self._ydat is None or self._xunc is None or self._yunc is None:
            messagebox.showerror('Error', 'No input data loaded')
            return

        # grab the data from the read in values and calculate the appropriate uncertainties
        self.xdat = np.array(self._xdat)
        self.ydat = np.array(self._ydat)
        self.xunc = np.array(self._xunc / self.sigvar.get())
        self.yunc = np.array(self._yunc / self.sigvar.get())

        # get the correlation array if it was chosen, otherwise get None
        if self.p_corr_var.get():
            self.p = np.array(self._p)
        else:
            self.p = np.zeros(len(self._p))

        # clear entry
        self.a_text.delete(0, 'end')
        self.b_text.delete(0, 'end')
        self.x_text.delete(0, 'end')
        self.mswd_text.delete(0, 'end')

        # fixed intercept checking:
        if self.force_int_var.get():
            try:
                self.afx = float(self.force_int_entry.get())
            except ValueError:
                messagebox.showerror('Error', 'Please define the fixed intercept or untoggle that specific option.')
                return
            # add a point the the input data that has 1e12 times smaller errors than the smallest error in the system
            xdat = np.zeros(len(self.xdat) + 1)
            ydat = np.zeros(len(self.ydat) + 1)
            xunc = np.zeros(len(self.xunc) + 1)
            yunc = np.zeros(len(self.yunc) + 1)
            p = np.zeros(len(self.p) + 1)
            # find the minimum error to add
            errintercept = np.min(np.array([np.min(self.xunc), np.min(self.yunc)])) / 1.e18
            # now add the point to the new array and then add all the existing data
            xdat[0] = 0.
            ydat[0] = self.afx
            xunc[0] = errintercept
            yunc[0] = errintercept
            p[0] = 0.
            for it in range(len(self.xdat)):
                xdat[it+1] = self.xdat[it]
                ydat[it+1] = self.ydat[it]
                xunc[it+1] = self.xunc[it]
                yunc[it+1] = self.yunc[it]
                p[it+1] = self.p[it]
            # now write back
            self.xdat = xdat
            self.ydat = ydat
            self.xunc = xunc
            self.yunc = yunc
            self.p = p
            # now calculate
            self.calcparams()
            self.calcunc()
            self.calcunc(calcxintunc=True)
        else:
            # run the calculation with the loaded x and y data
            self.calcparams()
            self.calcunc()
            self.calcunc(calcxintunc=True)

        # calculate the MSWD
        self.calcmswd()

        # update output
        slopestring = '%1.6e +/- %1.6e' % (self.slope, float(self.sigvar.get()) * self.slopeunc)
        if self.force_int_var.get():
            yinterstring = '%1.6e' % (self.yinter)
        else:
            yinterstring = '%1.6e +/- %1.6e' % (self.yinter, float(self.sigvar.get()) * self.yinterunc)
        xinterstring = '%1.6e +/- %1.6e' % (self.xinter, float(self.sigvar.get()) * self.xinterunc)
        mswdstring = '%.2f' % (self.mswd)
        self.a_text.insert(tk.END, slopestring)
        self.b_text.insert(tk.END, yinterstring)
        self.x_text.insert(tk.END, xinterstring)
        self.mswd_text.insert(tk.END, mswdstring)

        # write an output file
        outname = self.fname[0:len(self.fname) - 4] + '_out.txt'
        f = open(outname, 'w')
        f.writelines('All uncertainties are ' + str(self.sigvar.get()) + ' sigma!\n')
        # write header
        f.writelines('Input data\nxdata\txdata_unc\tydata\tydata_unc\tCorrelation_values\n')
        for it in range(len(self._xdat)):
            f.writelines(str(self._xdat[it]) + '\t' + str(float(self.sigvar.get()) * self._xunc[it]) + '\t' +
                         str(self._ydat[it]) + '\t' + str(float(self.sigvar.get()) * self._yunc[it]) + '\t' +
                         str(self._p[it]) + '\n')
        # now write the slope and stuff out
        f.writelines('\nLinear Regression after Mahon (1996):\n=====================================\n\n')
        # used correlation values or not?
        if self.p_corr_var.get():
            f.writelines('Regression calculation considering the correlation values\n')
        else:
            f.writelines('Regression calculation without considering the correlation values\n')
        # parameters
        slopestring = '%1.6e +/- %1.6e' % (self.slope, float(self.sigvar.get()) * self.slopeunc)
        if self.force_int_var.get():
            yinterstring = '%1.6e' % (self.yinter)
        else:
            yinterstring = '%1.6e +/- %1.6e' % (self.yinter, float(self.sigvar.get()) * self.yinterunc)
        xinterstring = '%1.6e +/- %1.6e' % (self.xinter, float(self.sigvar.get()) * self.xinterunc)
        mswdstring = '%.2f' % (self.mswd)
        # write to file
        f.writelines('All uncertainties are ' + str(self.sigvar.get()) + ' sigma:\n')
        if self.force_int_var.get():
            f.writelines('Intercept at x=0 fixed at y=' + str(self.afx) + '\n')
        f.writelines('Slope:       ' + slopestring + '\n')
        f.writelines('Y-Intercept: ' + yinterstring + '\n')
        f.writelines('X-Intercept: ' + xinterstring + '\n')
        f.writelines('MSWD:        ' + mswdstring)
        # flush and close
        f.flush()
        f.close()

    def plotdata(self):
        # subroutines for plotting function
        def xzeroplt():
            # x axis
            if defaultxlim[0] < 0. < defaultxlim[1]:
                a.set_xlim(defaultxlim)
            elif defaultxlim[0] > 0.:
                a.set_xlim([0., defaultxlim[1]])
            else:
                a.set_xlim([defaultxlim[0], 0.])

            # y axis
            if defaultylim[0] < self.yinter:
                minylim = defaultylim[0]
            else:
                minylim = self.yinter
            if defaultylim[1] > self.yinter:
                maxylim = defaultylim[1]
            else:
                maxylim = self.yinter
            a.set_ylim([minylim, maxylim])
            canvas.show()

        def resetxlim():
            a.set_xlim(defaultxlim)
            a.set_ylim(defaultylim)
            canvas.show()

        # check if input available
        if self.xdat is None or self.ydat is None or self.xunc is None or self.yunc is None:
            messagebox.showerror('Error', 'No input data loaded')
            return

        xdat, ydat, xunc, yunc = self._xdat, self._ydat, float(self.sigvar.get()) * self._xunc, \
                                 float(self.sigvar.get()) * self._yunc

        pltwin = tk.Toplevel()
        pltwin.geometry("800x600")

        f = Figure(figsize=(5, 4), dpi=100)
        a = f.add_subplot(111)

        a.errorbar(xdat, ydat, xerr=xunc, yerr=yunc, fmt='o', color='b', mfc='b', capsize=0, markersize=9)


        # default limits
        defaultxlim = a.get_xlim()
        defaultylim = a.get_ylim()

        # find min max in x and y, if x < 0
        if np.min(xdat - xunc) < 0 < np.max(xdat + xunc):
            xlims = np.array([np.min(xdat - xunc), np.max(xdat + xunc)])
        elif np.min(xdat - xunc) > 0:
            xlims = np.array([0., np.max(xdat + xunc)])
        else:
            xlims = np.array([np.min(xdat-xunc), 0.])

        # calculate fitline
        regy = self.slope * xlims + self.yinter

        # draw the fit line through the plot
        a.plot(xlims, regy, '-', linewidth=2, color='k')

        # set default limits
        a.set_xlim(defaultxlim)
        a.set_ylim(defaultylim)

        # put the parameters in the title:
        slopestring = '%1.3e +/- %1.3e' % (self.slope, float(self.sigvar.get()) * self.slopeunc)
        yinterstring = '%1.3e +/- %1.3e' % (self.yinter, float(self.sigvar.get()) * self.yinterunc)
        mswdstring = '%.2f' % (self.mswd)

        a.set_title('a * x + b / MSWD: ' + mswdstring + '\na = ' + slopestring + ' / b = ' + yinterstring)
        a.set_xlabel('x')
        a.set_ylabel('y')

        # a tk.DrawingArea
        canvas = FigureCanvasTkAgg(f, master=pltwin)
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        toolbar = NavigationToolbar2TkAgg(canvas, pltwin)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        # linlog toggle control
        tk.Button(pltwin, text='Plot to x=0', command=xzeroplt).pack(side=tk.LEFT)
        tk.Button(pltwin, text='Reset x limits', command=resetxlim).pack(side=tk.LEFT)

    def get_help(self):
        hlp_txt = 'This program makes a linear regression including errors in x and y ' \
                  'according to Mahon (1996). The linear regression is defined as\n\n' \
                  'y = a*x + b \n\n To use it, load a file first. You can ' \
                  'load either a text file with spaces between the columns (Load txt)' \
                  ', or a comma separated file (Load csv). Then hit calculate and the ' \
                  'slope and y-intercept will be calculated, including the appropriate ' \
                  'uncertainties.\n\n' \
                  'Choose if your file has 1 sigma or 2 sigma uncertainties using the ' \
                  'appropriate radio buttons. The program will give back the uncertainties ' \
                  'in the chosen format.\n\n' \
                  'Your file has to be structured as following: \n' \
                  '1st column: x values\n' \
                  '2nd column: uncertainties of x values\n' \
                  '3rd column: y values\n' \
                  '4th column: uncertainties of y values\n' \
                  '5th column: Correlation value \n' \
                  '   (optional: assumed 0 if not provided)\n\n' \
                  'LLNL-CODE-745740\n'
        messagebox.showinfo(title='Help', message=hlp_txt)

    def quit_app(self):
        self.master.quit()
        self.master.destroy()

    def calcparams(self):
        """
        Calculate the parameters, both intercepts as well as the slope of the regression
        :return:
        """
        bcalc = 1   # initial guess for slope
        bold = 0    # to compare later

        # read in selfs
        xdat = self.xdat
        xunc = self.xunc
        ydat = self.ydat
        yunc = self.yunc

        # run thorough the while loop
        whilecounter = 0
        whilecountermax = 1e5
        while np.abs((bold-bcalc) / bcalc) > 1.e-10 and whilecounter < whilecountermax:
            whilecounter += 1
            # prep for while loop, start before this line and compare bold to bcalc
            bold = bcalc
            # calculate xbar and ybar
            xbar = 0
            ybar = 0
            weightsum = 0
            for it in range(len(xdat)):
                wi = self.calc_wi(xunc[it], yunc[it], bcalc, self.p[it])
                xbar += xdat[it] * wi
                ybar += ydat[it] * wi
                weightsum += wi
            xbar /= weightsum
            ybar /= weightsum

            # now calculate b
            btop = 0   # top sum
            bbot = 0   # bottom sum

            for it in range(len(xdat)):
                xi = xdat[it]
                yi = ydat[it]
                sxi = xunc[it]
                syi = yunc[it]
                pit = self.p[it]
                wi = self.calc_wi(sxi, syi, bcalc, pit)
                ui = xi - xbar
                vi = yi - ybar
                # add to sums
                btop += wi**2. * vi * (ui * syi**2. + bcalc * vi * sxi**2. - pit * vi * sxi * syi)
                bbot += wi**2. * ui * (ui * syi**2. + bcalc * vi * sxi**2. - bcalc * pit * ui * sxi * syi)

            # new bcalc
            bcalc = btop / bbot

        # error message if whilecounter timed out
        if whilecounter == whilecountermax:
            messagebox.showwarning('Convergence warning', 'Warning! Your calculation might not have converged ' +
                                                          'properly. The difference between the last calculated '
                                                          'slope  and the current slope is: ' +
                                                          str(np.abs((bold-bcalc) / bcalc)) + ' You can ignore this ' +
                                                          'Message if this is an acceptable convergence for you.')

        # now that slope is determined, calculate the y intercept
        self.yinter = ybar - bcalc * xbar

        # now done, so write back slope
        self.slope = bcalc

        # calculate x intercept
        self.xinter = -self.yinter / self.slope

        # write back xbar and ybar
        self.xbar = xbar
        self.ybar = ybar

    def calcunc(self, calcxintunc=False):
        """
        Calculates the uncertainty of the slope and y
        :param calcxintunc: If it needs to calculate the x uncertainty, then set this to true
        :return:
        """
        if calcxintunc:
            # read in selfs
            # since this is for x uncertainty, simply switch it x and y.
            xdat = self.ydat
            xunc = self.yunc
            ydat = self.xdat
            yunc = self.xunc
            xbar = self.ybar
            ybar = self.xbar
            b = 1. / self.slope
        else:
            # read in selfs
            xdat = self.xdat
            xunc = self.xunc
            ydat = self.ydat
            yunc = self.yunc
            xbar = self.xbar
            ybar = self.ybar
            b = self.slope

        # let us first calculate the derivatives
        # dell theta / dell b (dthdb) calculation
        sum1 = 0.
        sum2 = 0.
        for it in range(len(xdat)):
            xi = xdat[it]
            yi = ydat[it]
            sxi = xunc[it]
            syi = yunc[it]
            pit = self.p[it]
            wi = self.calc_wi(xunc[it], yunc[it], b, pit)
            ui = xi - xbar
            vi = yi - ybar
            sxyi = pit * sxi * syi
            sum1 += wi**2. * (2 * b * (ui * vi * sxi**2. - ui**2. * sxyi) + (ui**2. * syi**2. - vi**2 * sxi**2.))
            sum2 += wi**3. * (sxyi - b * sxi**2.) * (b**2. * (ui * vi * sxi**2 - ui**2 * sxyi) +
                                                     b * (ui**2 * syi**2 - vi**2 * sxi**2) -
                                                     (ui * vi * syi**2 - vi**2 * sxyi))
        dthdb = sum1 + 4. * sum2

        # calculate the sum of all weights
        wksum = 0.
        for it in range(len(xdat)):
            wksum += self.calc_wi(xunc[it], yunc[it], b, self.p[it])

        # now calculate sigasq and sigbsq
        sigasq = 0.
        sigbsq = 0.
        for it in range(len(xdat)):
            sxi = xunc[it]
            syi = yunc[it]
            pit = self.p[it]
            wi = self.calc_wi(sxi, syi, b, pit)
            sxyi = pit * sxi * syi

            # calculate dell theta / dell xi and dell theta / dell yi
            dthdxi = 0.
            dthdyi = 0.
            for jt in range(len(xdat)):
                xj = xdat[jt]
                yj = ydat[jt]
                sxj = xunc[jt]
                syj = yunc[jt]
                pjt = self.p[jt]
                wj = self.calc_wi(sxj, syj, b, pjt)
                uj = xj - xbar
                vj = yj - ybar
                sxyj = pjt * sxj * syj
                # add to dthdxi and dthdyi
                dthdxi += wj**2. * (kron(it, jt) - wi / wksum) * (b**2 * (vj * sxj**2 - 2 * uj * sxyj) +
                                                                  2 * b * uj * syj**2 - vj * syj**2)
                # correct equation! not equal to equation 21 in Mahon (1996)
                dthdyi += wj ** 2. * (kron(it, jt) - wi / wksum) * (b ** 2 * uj * sxj ** 2 + 2 * vj * sxyj -
                                                                    2 * b * vj * sxj**2. - uj * syj ** 2)

            # now calculate dell a / dell xi and dell a / dell yi
            dadxi = -b * wi / wksum - xbar * dthdxi / dthdb
            dadyi = wi / wksum - xbar * dthdyi / dthdb

            # now finally add to sigasq and sigbsq
            sigbsq += dthdxi**2. * sxi**2. + dthdyi**2. * syi**2. + 2 * sxyi * dthdxi * dthdyi
            sigasq += dadxi**2. * sxi**2. + dadyi**2. * syi**2. + 2 * sxyi * dadxi * dadyi

        # now divide sigbsq
        sigbsq /= dthdb**2.

        # now write slope uncertainty and y intercept uncertainty back to class
        if calcxintunc:
            self.xinterunc = np.sqrt(sigasq)
        else:
            self.yinterunc = np.sqrt(sigasq)
            self.slopeunc = np.sqrt(sigbsq)

    def calcmswd(self):
        xdat, ydat, xunc, yunc = self.xdat, self.ydat, self.xunc, self.yunc
        mswd = 0.
        for it in range(len(xdat)):
            xi = xdat[it]
            yi = ydat[it]
            sxi = xunc[it]
            syi = yunc[it]
            pit = self.p[it]
            wi = self.calc_wi(sxi, syi, self.slope, pit)
            mswd += wi * (yi - self.slope * xi - self.yinter)**2.

        # now divide by degrees of freedom minus 2, since 2 fixed parameters
        mswd /= (len(xdat) - 2.)
        self.mswd = mswd

    def calc_wi(self, sx, sy, b, p):
        return 1. / (sy**2 + b**2 * sx**2 - 2 * b * p * sx * sy)


def kron(i, j):
    # calculates Kronecker delta
    if i == j:
        return 1.
    else:
        return 0.


# run GUI
root = tk.Tk()
my_gui = Mahon(root)
root.mainloop()
