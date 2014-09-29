# -*- coding: utf-8 -*-
"""
Created on Tue Aug 12 15:47:32 2014

@author: Andrew
"""

from automatedAFManalysis.afm_sql_code.afmsqlcode import SQLConnection, SQLData, FeatureInfo
from scipy import signal
import numpy as np
from PyQt4 import QtCore, QtGui
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from automatedAFManalysis.dataviewer.autoanalysisgui import Ui_MainWindow
from automatedAFManalysis.dataviewer.dataviewer import Database, InputWindow, ConnectDialog
import csv
from multiprocessing import freeze_support

class Errbox(QtGui.QMessageBox):
    '''
    Simple error box. The text can be set with errbox.setText().
    '''

    def __init__(self, parent=None):
        QtGui.QMessageBox.__init__(self)
        self.setStandardButtons(QtGui.QMessageBox.Ok)
        self.setDefaultButton(QtGui.QMessageBox.Ok)
        self.setWindowTitle('ERROR')

class MainWindow(QtGui.QMainWindow, Ui_MainWindow):

    def __init__(self, parent=None, f=QtCore.Qt.WindowFlags()):
        QtGui.QMainWindow.__init__(self, parent, f)
        self.setupUi(self)

class Viewer(MainWindow, InputWindow):

    '''
    This class describes the GUI and related methods for the program.  Since it
    borrows heavily from the design of the GUI of the dataviewer program, it
    inherits from that programs InputWindow class.

    Member Variables/Objects:•

        Private:
            _cell_flag: Serves no purpose in this program, vestigial from the
                        InputWindow class.  Always set to true.
            _tethers: Serves no purpose in this program, vestigial from the
                      InputWindow class.  Always set to true.
            _database: The MLDatabase class object for the application, handles
                       most SQL interactions
            _analyze_task: The instance of the AnalyzeTask class used to analyze the
                        data associated with the current selected parameters
            _tscatplot: The scatter plot figure for the tethers
            _thistplot: The histogram figure for the tethers
            _utscatplot: The scatter plot figure for the untethered bonds
            _uthistplot: The histogram figure for the untethered bonds
            cancelsig: Signal used to cancel the _analyze_task
            _tforces: Calculated tether forces
            _tslopes: Calculated tether slopes
            _utforces: Calcuated untethered bond forces
            _utslopes: Calculated untethered bond slopes

    Methods:

        Private:
            __init__: Builds window, sets initial values
            _bind_callbacks: Binds the callbacks for the UI
            _cancelcallback: Callback that runs when you press the cancel
                             button, kills execution of the _analyze_task
            _plot_stuff: Creates plots upon receiving a signal from the
                         _analyze_task
            _onprogress: Updates progress bar on signal from the _analyze_task
            _error: Creates an error dialog
            _plot_callback: Callback that runs when user presses the plot button
            _save_callback: Callback for the save menu action
            _connect_callback: Callback for the connect menu action
    '''

    cancelsig = QtCore.pyqtSignal(bool)

    def __init__(self, parent=None):
        '''
        The __init__ method creates the window and initializes the display,
        including populating the various menus with valid options.
        '''

        #Intialize various members, setup window
        QtGui.QDialog.__init__(self, parent)
        MainWindow.__init__(self)
        self.setWindowTitle('Single Bond Automatic Analysis')
        self._cell_flag = True
        self._tethers = True
        self._analyze_task = None
        self._tforces = []
        self._tslopes = []
        self._utforces = []
        self._utslopes = []

        #Connect controls to callbacks
        self._bind_callbacks()
        self.plot_button.setEnabled(False)
        #Initialize progressbar
        self.progressBar.setValue(0)

        #Setup plots
        self._tscatplot = MyMplCanvas(self.tscattframe,
                                   width=5, height=8, dpi=100)
        self._thistplot = MyMplCanvas(self.thistframe,
                                   width=5, height=8, dpi=100)
        self._utscatplot = MyMplCanvas(self.utscattframe,
                                   width=5, height=8, dpi=100)
        self._uthistplot = MyMplCanvas(self.uthistframe,
                                   width=5, height=8, dpi=100)
        self._plotstuff({'tetherforces':[], 'tetherslopes':[],
                      'untetheredforces':[], 'untetheredslopes':[],
                      'truetetherforces':[], 'truetetherslopes':[],
                      'trueuntetheredforces':[], 'trueuntetheredslopes':[],
                      'confusionmatrix':[[0, 0, 0], [0, 0, 0], [0, 0, 0]]})

        #Hides some vestigial interface options from InputWindow
        self.holdbox.hide()
        self.groupBox.hide()
        self.attachments.hide()

        #Setup confusionmatrix tableview
        self.confusionmatrix.setHorizontalScrollBarPolicy(1)
        self.confusionmatrix.setVerticalScrollBarPolicy(1)


    def _bind_callbacks(self):
        '''
        The _bind_callbacks method connects the various controls to the methods
        that they run.
        '''

        InputWindow._bind_callbacks(self)
        QtCore.QObject.connect(self.cancelbutton, QtCore.SIGNAL('clicked()'),
                               self._cancelcallback)

    def _cancelcallback(self):
        '''
        The _cancelcallback method runs when the user clicks the cancel button.
        It stops all current fitting/classification tasks in the _analyze_task
        object.
        '''

        #Disconnect _analyze_task from progressbar
        self._analyze_task.notifyProgress.disconnect()
        #Send abort signal
        self.cancelsig.emit(False)
        #Reset progressbar/pushbuttons
        self._onprogress(0)
        self.cancelbutton.setEnabled(False)
        self.plot_button.setEnabled(True)

    def _plot_callback(self):
        '''
        The _plot_callback method runs when the plot button is pushed.  It starts
        the _analyze_task and connects it to the matplotlib graphs/progressbar.
        '''

        try:
            del self._analyze_task
        except:
            pass
        try:
            self.plot_button.setEnabled(False)
            self.cancelbutton.setEnabled(True)
            self.progressBar.setRange(0, self._database.totalcycs)
            #Creates _analyze_task and starts it
            self._analyze_task = AnalyzeTask(self._database.list_of_ids, self, self._conn)
            self._analyze_task.start()
            #Connects various signals between the viewer and the _analyze_task
            self._analyze_task.notifyProgress.connect(self._onprogress)
            self._analyze_task.datasig.connect(self._plotstuff)
            self._analyze_task.errsig.connect(self.error)
            self.cancelsig.connect(self._analyze_task.abort)
        except Exception as errtext:
            err_box = Errbox(self)
            err_box.setText(str(errtext))
            err_box.exec_()

    def _plotstuff(self, data):
        '''
        The _plotstuff method plots the data generated by the _analyze_task, it
        runs when it gets a signal from the _analyze_task.

        Inputs:
            data: Dictionary containing:
                     {'tetherforces':list, 'tetherslopes':list,
                      'untetheredforces':list, 'untetheredslopes':list,
                      'truetetherforces':list, 'truetetherslopes':list,
                      'trueuntetheredforces':list, 'trueuntetheredslopes':list,
                      'confusionmatrix':3x3 list of lists}
        '''

        self.cancelbutton.setEnabled(False)
        self.plot_button.setEnabled(True)

        forces = [data['tetherforces'], data['untetheredforces'],
                  data['truetetherforces'], data['trueuntetheredforces']]
        slopes = [data['tetherslopes'], data['untetheredslopes'],
                  data['truetetherslopes'], data['trueuntetheredslopes']]

        self._tforces = forces[0]
        self._tslopes = slopes[0]
        self._utforces = forces[1]
        self._utslopes = slopes[1]
        confusion = data['confusionmatrix']

        hist_plot_names = ['Histogram of Calculated Tether Forces',
                           'Histogram of Calculated Untethered Forces',
                           'Histogram of True Tether Forces',
                           'Histogram of True Untethered Forces']
        scat_plot_names = ['Scatter Plot of Calculated Tether Forces',
                           'Scatter Plot of Calculated Untethered Forces',
                           'Scatter Plot of True Tether Forces',
                           'Scatter Plot of True Untethered Forces']
        hist_plot_list = [self._thistplot.calc, self._uthistplot.calc,
                          self._thistplot.true, self._uthistplot.true]
        scat_plot_list = [self._tscatplot.calc, self._utscatplot.calc,
                          self._tscatplot.true, self._utscatplot.true]
        figlist = [self._thistplot, self._uthistplot,
                   self._tscatplot, self._utscatplot]

        #Plot the histograms
        binlist = range(0, 125, 5)
        for i, graph in enumerate(hist_plot_list):
            try:
                graph.axes.cla()
                graph.axes.hist(forces[i], normed=True, bins=binlist)
            #If no data, don't plot
            except:
                pass
            #Limits, labels
            finally:
                graph.set_xlim(0, 100)
                graph.set_xlabel('Force (pN)', fontsize=12)
                graph.set_title(hist_plot_names[i])

        #Plot scatterplots
        for i, graph in enumerate(scat_plot_list):
            try:
                graph.axes.cla()
                graph.axes.plot(slopes[i], forces[i], '.')
            #If no data, don't plot
            except:
                pass
            #Limits, labels
            finally:
                graph.set_xlim(100, 10000)
                graph.set_ylim(0, 100)
                graph.set_xscale('log')
                graph.set_xlabel('Loading Rate (pN/s)', fontsize=12)
                graph.set_ylabel('Force (pN)', fontsize=12)
                graph.set_title(scat_plot_names[i])
        #Draw 'em
        for fig in figlist:
            fig.draw()

        #Displays the confusion matrix and sizes the tableview appropriately
        tabmodel = Mytablemodel()
        self.confusionmatrix.setModel(tabmodel)
        tabmodel.setdat(confusion)
        self.confusionmatrix.resizeColumnsToContents()
        w = 0
        for i in xrange(0, 4):
            w = w + self.confusionmatrix.columnWidth(i)
        h = 0
        for i in xrange(0, 3):
            h = h + self.confusionmatrix.rowHeight(i)
        h = h + self.confusionmatrix.rowHeight(0)
        self.confusionmatrix.setFixedWidth(w)
        self.confusionmatrix.setFixedHeight(h)

    def _save_callback(self):
        '''
        The save_data method saves the extracted slopes and forces to a CSV
        file when the user clicks the save menu option
        '''

        #Save tether data:
        #Prompt user for save location:
        fname = QtGui.QFileDialog.getSaveFileName(self, 'Save Tethered Data Location',
                                                  'c:\\', '*.txt')
        fname = str(fname)
        #If save location provided, save
        if fname != '':
            with open(fname, 'wb') as dataout:
                dataout = csv.writer(dataout, delimiter=',')
                dataout.writerow(['slopes (pN/s)', 'forces (pN)'])
                for i in xrange(0, len(self._tforces)):
                    dataout.writerow([self._tslopes[i], self._tforces[i]])

        #Save untethered data:
        fname = QtGui.QFileDialog.getSaveFileName(self, 'Save Untethered Data Location',
                                                  'c:\\', '*.txt')
        fname = str(fname)
        if fname != '':
            with open(fname, 'wb') as dataout:
                dataout = csv.writer(dataout, delimiter=',')
                dataout.writerow(['slopes (pN/s)', 'forces (pN)'])
                for i in xrange(0, len(self._utforces)):
                    dataout.writerow([self._utslopes[i], self._utforces[i]])

    def _onprogress(self, i):
        '''
        The _onprogress method updates the progress bar. It updates on signals
        from the _analyze_task.

        Inputs:
            i: Progress of the _analyze_task
        '''

        if i > 0 and i < self._database.totalcycs-1:
            self.proglabel.setText('Analyzing...')
        elif i == self._database.totalcycs-1:
            self.proglabel.setText('Done.')
        self.progressBar.setValue(i)

    def _connect_callback(self):
        '''
        The _connect_callback runs when the user picks the connect menu option.
        It connects to an SQL database with experimental data locations and
        then retrieves valid parameters for various menu options.
        '''

        connectdialog = ConnectDialog()
        try:
            self._conn = SQLConnection(connectdialog.username, connectdialog.password, connectdialog.host, connectdialog.port)
            self._database = MLDatabase(self._conn)
            self._init_tables()
            self.plot_button.setEnabled(True)
        except Exception as e:
            self.error(str(e))
            return


    def error(self, errtext):
        '''
        Creates generic error dialog.
        Inputs:
            errtext: The error text
        '''
        err_box = Errbox(self)
        err_box.setText(errtext)
        err_box.exec_()
        return

class MLDatabase(Database):
    '''
    The MLDatabase object handles most of the interactions with the SQL
    database, aside from retrieving individual cycle time series.  It inherits
    from the Database class used in the dataviewer program.

    Objects:
        Public:
            list_of_ids: List of all raw data file ids that correspond to the
                         current list of parameters
            totalcycs: Total number of data cycles that correspond to the
                       current list of parameters

    Methods:
        Public:
            get_valid_params: Returns remaining valid parameter choices for the
                              experimental data given a choice of alrady chosen
                              parameters
            update_list: Updates the objects list_of_ids, returns
                                 total number of cycles for selected parameters

    '''

    def __init__(self, conn):
        '''
        The __init__ method runs the Database __init__ method.
        '''
        self.list_of_ids = []
        self.totalcycs = 0
        Database.__init__(self, conn)
        
    def get_valid_params(self, tuple_in, cell_flag, tetherflag):
        '''
        The get_valid_params method returns the valid parameter choices given a
        set of chosen parameters.  For example, if you select a given date
        you'll only see buffers and pipette objects that were used on that day.

        Inputs:
            tuple_in: Tuple of lists containing the various constraints for
                      each selection parameter
                      (pipette object list, cantilever object list,
                       buffer list, date list, loading rate list)
            cell_flag: True for cells, False for beads
            tetherflag: True for tethered, false for untethered

        Outputs:
            valid_params: List of valid parameters
        '''

        names = ['pipette_obj', 'cantilever_obj',
                 'buffer', 'date', 'loading_rate']
        valid_params = []
        #Loop over selected parameters to build their tables
        for i in xrange(0, len(names)):
            #We build constraints with knockout set to i (we don't want the
            #case where selecting a date means other dates are invalid choices
            #for example, date choices should not constrain valid date choices)
            conditions = self._build_conditions(tuple_in, i,
                                                cell_flag, tetherflag,
                                                buildjoin=False)
            select = ''.join(['SELECT DISTINCT ', names[i],
                         ' FROM afm_data.filelist'])
            statement = ''.join([select, conditions])
            #Build the valid parameters out:
            #ANY is always a valid choice
            temp = ['ANY']

            results = self._conn.execute(statement)
            for params in results:
                temp.append(params[0])
            temp = [temp]
            valid_params.append(temp)
        #Converts loading rates to floats
        for params in xrange(len(valid_params[4][0])):
            if valid_params[4][0][params] == 'ANY':
                pass
            else:
                valid_params[4][0][params] = float(str(valid_params[4][0][params]))

        return valid_params
    def update_list(self, tuple_in, cell_flag, tetherflag):
        '''
        The update_list method updates the list_of_ids.

        Inputs:
            tuple_in: Tuple of lists containing user selected selection
                      parameters
                      (pipette object values, cantilever object values,
                      buffer values, date values, loading rate values)
            cell_flag: Flag indicating whether bead or cell data is to be
                       retrieved, True for cells
            tether_flag: Flag indicating whether tethered or untethered data
                         is to be retrieved, True for tethered

        Outputs:
            A tuple consisting of:
                totalcycs*2: The total number of cycles meeting the user
                             supplied selection criteria, there is two of them
                             because the calling methods expect two returns
        '''

        #Get features:
        #SQL query statement made of two parts, a and b, where a is selection
        #information and b is constraints/conditions
        select = ('SELECT id, number_of_cycles FROM afm_data.filelist ')

        conditions = self._build_conditions(tuple_in, 9, cell_flag, tetherflag, buildjoin=False)
        statement = ''.join([select, conditions])
        result = self._conn.execute(statement)
        self.list_of_ids = []
        listoftotcycs = []
        for rows in result:
            self.list_of_ids.append(int(rows[0]))
            listoftotcycs.append(int(rows[1]))
        self.totalcycs = sum(listoftotcycs)

        #Calling functions expect two returns
        return (self.totalcycs, self.totalcycs)


class MyMplCanvas(FigureCanvas):
    '''
    Matplotlib canvas object
    '''

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        '''
        Sets up the canvas, creates a figure, sets various parameters.
        '''

        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.fig.patch.set_alpha(1)
        self.fig.patch.set_facecolor([.95, 0.95, 0.95])
        self.fig.subplots_adjust(left=0.18, top=0.95, bottom=0.07, hspace=.5)
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self,
                                   QtGui.QSizePolicy.Expanding,
                                   QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.calc = self.fig.add_subplot(211)
        self.true = self.fig.add_subplot(212)

class EmittingSQLDat(SQLData):
    '''
    The EmittingSQLDat class handles retrieval, classification and certain
    processing tasks of individual cycle time series.  It fetches its data from
    a location specified by the SQL database. It also emits progress updates
    via the signal object given to it.

    Member Variables/Objects:

        Public:
            emitter: Signal linked to a progress bar
            prog: Integer counter, increments up with each classified cycle

    Member Methods:

        Private:
            __init__: Sets the emitter and prog members to initial values, runs
                      SQLData.__init__

        Public:
            get_next_cycle: Gets the next valid (non-junk) cycle, emits a
                            progress signal with each classified cycle
    '''

    def __init__(self, fileid, method_flag, conn, emitter, prog, thresh_val=25,
                 win=21, parent=None, cyclist=None):
        self.emitter = emitter
        self.prog = prog
        cycs = cyclist
        SQLData.__init__(self, fileid, method_flag, conn, thresh_val=25,
                         win=21, parent=None, cyclist=cycs)
    def get_next_cycle(self):
        '''
        The get_next_cycle method gets the next cycle that passes the cycle
        selection criteria. It also moves the current cycle to that cycle and
        emits an update signal to the viewers progress bar.

        Output:
            result: None if no cycle found, else a 2d numpy array with the
                    force vs time data
        '''

        self.prog = self.prog + 1
        self.emitter.emit(self.prog)
        return SQLData.get_next_cycle(self)

class Mytablemodel(QtCore.QAbstractTableModel):
    '''
    The table model used for bead data.
    '''

    def __init__(self):
        QtCore.QAbstractTableModel.__init__(self)
        self.initdata = []
        self.rows = 3
        self.columns = 4
        self.header = ["", "Junk", "Tethered", "Untethered"]

    def rowCount(self, parent=QtCore.QModelIndex()):
        return self.rows

    def columnCount(self, parent=QtCore.QModelIndex()):
        return self.columns

    def data(self, index, role=QtCore.Qt.DisplayRole):
        if not index.isValid():
            return None
        elif role != QtCore.Qt.DisplayRole:
            return None
        return self.initdata[index.column()][index.row()]

    def headerData(self, col, orientation, role):
        if orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
            return self.header[col]
        return None

    def setdat(self, data):
        self.initdata = [["Junk", "Tethered", "Untethered"],
                         data[0], data[1], data[2]]

class AnalyzeTask(QtCore.QThread):
    '''
    The AnalyzeTask class handles the fitting, feature extraction and classification
    of a given list of raw data file ids whose time series locations are
    stored in the SQL database.

    Member Objects/Variables:

        Public:
            errsig: Signal used to send an error message to the caller
            notifyProgress: Signal used to send progress updates to the caller
            datasig: Signal used to send finished slope/force data as well as
                     the confusion matrix to the caller

        Private:
            _prog: Cycles analyzed by the instance so far, updated after each
                  raw data file analyzed
            _conn: SQL database connection
            _listofids: List of file ids of raw data files to be analyzed
            _totalanalyzedcycs: Total cycles that have been analyzed by hand
                               previously, used to build confusion matrix
            _abortf: When set to true kills analysis
            _data: EmittingSQLdata object used by the task
            _cycdat: The force time data for the current cycle from the _data
                     object
            _smoothparam: Gaussian window size for _smooth method
            _output: The output, a dictionary of form:
                     {'tetherforces':list, 'tetherslopes':list,
                      'untetheredforces':list, 'untetheredslopes':list,
                      'truetetherforces':list, 'truetetherslopes':list,
                      'trueuntetheredforces':list, 'trueuntetheredslopes':list,
                      'confusionmatrix':3x3 list of lists}
                     sent out of task via datasig signal

    Methods:

        Public:
            abort: Kills the currently running analysis
            run: Runs the analysis

        Private:
            _get_true_data: Gets the ground truth (hand analyzed) tethered and
                            untethered forces and slopes
            _smooth: Smooths data with a gaussian window, data must be a 1D
                     array
            _find_fitting_range: Find the fit interval of a tethered cycle
            _find_crossover_force_loc: Finds the location of the crossover of a
                                   tethered cycle
            _analyze_tethers: Analyzes single bond cycles classified as
                              tethered
            _analyze_untethered: Analyzes single bond cycles classified as
                                 untethered
            _get_features: Gets the slope and force of a single bond cycle
            _calc_rsq: Calculates the R sq value of a linear fit over a given
                       interval
            _update_confusion_mat: Updates the values of some confusion matrix
                                   elements when a single bond cycle is
                                   analyzed
            _generate_all_features: Analyzes the data and gets the slopes,
                                    forces and confusion matrix
    '''

    errsig = QtCore.pyqtSignal(str)
    notifyProgress = QtCore.pyqtSignal(int)
    datasig = QtCore.pyqtSignal(dict)

    def __init__(self, ids, viewer, conn):
        '''
        The __init__ method sets initial values for instance objects/variables.

        Inputs:
            ids: A list of raw data file id numbers to be analyzed by the
                 class instance
            viewer: The main application window/GUI
        '''

        self._prog = 0
        self._conn = conn
        self._listofids = ids
        self._totalanalyzedcycs = 0
        self._abortf = False
        self._output = {'tetherforces':[], 'tetherslopes':[],
                      'untetheredforces':[], 'untetheredslopes':[],
                      'truetetherforces':None, 'truetetherslopes':None,
                      'trueuntetheredforces':None, 'trueuntetheredslopes':None,
                      'confusionmatrix':None}
        self._data = None
        self._cycdat = None
        self._smoothparam = 1
        QtCore.QThread.__init__(self)

    def abort(self, signal):
        '''
        The abort method kills currently running analysis.  It is triggered by
        a signal from the Viewer.
        '''
        self._abortf = True

    def run(self):
        '''
        The run method runs when .start() is called.  It conducts all the
        analysis of the cycles associated with the raw data files in _listofids.
        '''

        try:

            self._get_true_data()
            self._generate_all_features()
            if self._output['tetherforces'] != None:
                try:
                    self.datasig.emit(self._output)
                except Exception as errtext:
                    self.errsig.emit(str(errtext))
            else:
                self._data.pool.terminate()
                del self._data
        except Exception as errtext:
            self.errsig.emit(str(errtext))

    def _get_true_data(self):
        '''
        The _get_true_data method retrieves the ground truth (analyzed by hand)
        numbers for the forces and slopes for a given _listofids.
        '''

        #Gets true tether slopes and forces
        s = ('SELECT forces, slope FROM afm_data.cell_data '
             'INNER JOIN afm_data.filelist '
             'ON filelist.id = file_id '
             'WHERE file_id=%s '
             'AND tethered=\'Y\';')

        self._output['truetetherforces'] = []
        self._output['truetetherslopes'] = []
        for ids in self._listofids:
            for rows in self._conn.execute(s%ids):
                self._output['truetetherforces'].append(float(rows[0]))
                self._output['truetetherslopes'].append(float(rows[1]))

        #Gets true untethered slopes and forces
        s = ('SELECT forces, slope FROM afm_data.cell_data '
             'INNER JOIN afm_data.filelist '
             'ON filelist.id = file_id '
             'WHERE file_id=%s '
             'AND tethered=\'N\';')

        self._output['trueuntetheredforces'] = []
        self._output['trueuntetheredslopes'] = []
        for ids in self._listofids:
            for rows in self._conn.execute(s%ids):
                self._output['trueuntetheredforces'].append(float(rows[0]))
                self._output['trueuntetheredslopes'].append(float(rows[1]))

    def _smooth(self, ydata):
        '''
        Creates smoothed ydata using convolution with a gaussian window.

        Inputs:
            ydata: The input data (1d Numpy array)
            gwin: The size of the smoothing Gaussian window

        Output:
            smoothed_y: The smoothed array, same length as input
        '''

        gwin =  self._smoothparam
        halfgwin = int(np.floor(gwin / 2))
        #Makes window
        window = (signal.gaussian(gwin, std=np.floor(gwin/6)) /
                  np.sum(signal.gaussian(gwin, std=np.floor(gwin/6))))
        #Smoothes force data
        smoothed_y = np.convolve(ydata, window, 'valid')
        #Pads data so it is the same length as the input data
        smoothed_y = np.hstack((np.ones((halfgwin)) * smoothed_y[0],
                                smoothed_y, np.ones((halfgwin)) *
                                smoothed_y[smoothed_y.size-1]))
        return smoothed_y

    def _find_fitting_range(self, graphchunk, dt):
        '''
        The _find_fitting_range method finds the range over which to fit a line
        to a tether force curve.  It does this by finding the location of the
        maximum slope of the input data graphchunk and then setting the fit
        interval to a certain range to the left and right of that.

        Inputs:
            graphchunk: The input data, a chunk of an experimental cycle, a
                        2d numpy array of force and time
            dt: The time interval between datapoints in graphchunk

        Outputs:
            fitwinst: The start of the fitting interval
            fitwinen: The end of the fitting interval
            maxloc: The location of the maximum slope
        '''

        #Find the gradient of the data
        slope = np.gradient(graphchunk[1, :], dt)
        #The half size of the fitting window, at 100 kHz 2000 seems pretty good
        halfwinsize = int(2000*0.00001/dt)
        #Initial gaussian window size to smooth the gradient with, we don't
        #the gradient shold be slowly changing so setting it large works well
        self._smoothparam = np.min([int(8001*0.00001/dt), graphchunk[1, :].size-5])
        if self._smoothparam % 2 == 0:
            self._smoothparam = self._smoothparam + 1
        #If the max slope is right up against the left or right endpoint, do
        #a little less smoothing
        flag = False
        while flag == False:
            smoothslope = self._smooth(slope)
            if self._smoothparam > 1000:
                if (smoothslope[0] > 0.95*np.max(smoothslope) or
                    smoothslope[-1] > 0.95*np.max(smoothslope)):
                    self._smoothparam = self._smoothparam-500
                else:
                    flag = True
            #Don't let the smoothing parameter get too small though
            else:
                self._smoothparam = 2001
                flag = True
        #Smooth the gradient and find the location of its maximum
        smoothslope = self._smooth(slope)
        maxloc = np.argmax(smoothslope)
        #fitwinen and fitwinst are maxloc+-halfiwnsize unless that is outside
        #the range of the data, then we set it to the corresponding endpoint of
        #the data
        fitwinen = maxloc + halfwinsize
        if fitwinen >= smoothslope.size:
            fitwinen = smoothslope.size-1
        fitwinst = maxloc - halfwinsize
        if fitwinst < 0:
            fitwinst = 0

        return (fitwinst, fitwinen, maxloc)

    def _find_crossover_force_loc(self, graphchunk, maxloc, line, dt):
        '''
        The _find_crossover_force method finds the location of the crossover
        event within a graphchunk given the location of its maximum, the
        parameters of a linear fit and the time between each datapoint

        Inputs:
            graphchunk: The input data, a chunk of an experimental cycle, a
                        2d numpy array of force and time
            maxloc: Location of the maximum gradient
            line: List containing the fit parameters (slope and yint) of the
               linear fit of the data.
            dt: The time interval between datapoints in graphchunk

        Outputs:
            crossoverloc: The location of crossover
        '''

        #Set initial gaussian window size to smooth the data with
        self._smoothparam = np.min([int(2500*0.00001/dt), graphchunk[1, :].size-5])
        if self._smoothparam % 2 == 0:
            self._smoothparam = self._smoothparam + 1
        #Look for the first place the fit residual is less than -2.5 pN,
        #starting at maxloc. If this ends up being within the _smoothparam
        #(the size of the Gaussian window) of the end of the file decrease
        #_smoothparam, resmooth and do it again.
        flag = False
        while flag == False:
            smoothedresid = (self._smooth(graphchunk[1, :]) -
                             (line[1]+graphchunk[0, :]*line[0]))
            crossoverloc = maxloc
            try:
                while smoothedresid[crossoverloc] > -2.5:
                    crossoverloc = crossoverloc + 1
            except IndexError:
                crossoverloc = smoothedresid.size - 1
            if crossoverloc > smoothedresid.size - self._smoothparam:
                self._smoothparam = self._smoothparam-500
                if self._smoothparam < 0:
                    self._smoothparam = 1
                    flag = True
            else:
                flag = True
        return crossoverloc

    def _analyze_tethers(self, dis_ind, zer):
        '''
        The _analyze_tethers method analyzes identified tether cycles and
        extracts the crossover force and the slope of the linear portion of the
        graph.

        Inputs:
            dis_ind: 1D Numpy array with the indices of the discontinuities in
                     the current cycle
            zer: 1D Numpy array with the average values of the force
                 immediately after the discontinuities (these are zero points)

        Outputs:
            A 2 element tuple containing the slope of the fitted line and the
            value of the crossover force.  If for whatever reason analysis
            fails then this returns (False, False). Slope is in pN/s, force
            is in pN.
        '''

        #If the cycle has more than one discontinuity we break it into chunks
        #and take the chunk that begins with the second to last discontinuity
        #and ends with the last one
        if dis_ind.size > 1:
            st = dis_ind[-2]
            en = dis_ind[-1]
        #If it has one discontinuity the chunk starts at the beginning of the
        #data and runs to the discontinuity
        elif dis_ind.size == 1:
            st = 0
            en = dis_ind[0]
        #With no discontinuities we assume the classifier screwed up and return
        #False values
        elif dis_ind.size == 0:
            return (False, False)
        #The piece of the cycle we are working with
        graphchunk = self._cycdat[:, st:en]
        #Time between data points
        dt = graphchunk[0, 2]-graphchunk[0, 1]
        #Finds the fit interval for the linear fit
        (fitwinst, fitwinen, maxloc) = self._find_fitting_range(graphchunk, dt)
        #Fits the data
        try:
            line = np.polyfit(graphchunk[0, fitwinst:fitwinen],
                              graphchunk[1, fitwinst:fitwinen], 1)
        except TypeError:
            return (False, False)
        #Very low slopes for the linear fit are symptomatic of bad fits
        if line[0] < 100:
            return (False, False)
        #Get the crossover location
        cossoverloc = self._find_crossover_force_loc(graphchunk, maxloc, line, dt)
        #Calculate the R squared of the fit, if it's good, return the slope and
        #the value of the fit at crossover (minus the location of zero force)
        if self._calc_rsq(graphchunk, fitwinst, cossoverloc, line) > 0.95:
            return (line[0], line[1] + line[0]*graphchunk[0, cossoverloc] - zer[-1])
        #If the fit is bad return False values
        else:
            return (False, False)

    def _analyze_untethered(self, force_abs, dis_ind, slopes, y_int, zer):
        '''
        The _analyze_untethered method analyzes identified untethered cycles
        and extracts the slope of the linear fit and the value of the rupture
        force.

        Inputs:
            force_abs: 1D numpy array with the absolute values of the linear
                       fits from fin_peaks at the discontinuities
            dis_ind: 1D Numpy array with the indices of the discontinuities in
                     the current cycle
            slopes: 1D Numpy array with the slopes in pN/index of the linear
                    fits from fin_peaks at the discontinuities
            y_int: 1D Numpy array of the Y intercept values of the linear fits
                   from fin_peaks at the discontinuities
            zer: 1D Numpy array with the average values of the force
                 immediately after the discontinuities (these are zero points)

        Outputs:
            A 2 element tuple containing the slope of the fitted line and the
            value of the crossover force.  If for whatever reason analysis
            fails then this returns (False, False). Slope is in pN/s, force
            is in pN.
        '''

        #With no discontinuities, we assume the classifier screwed up

        if dis_ind.size == 0:
            return (False, False)
        else:
            #If the R squared of the linear fit is good enough, convert the
            #slope to pN/s (from pN/index) and return it and the force - the
            #zero point
            dt = self._cycdat[0, 2] - self._cycdat[0, 1]
            self._smoothparam = 161
            if (self._calc_rsq(self._cycdat, (dis_ind[-1] - 500),
                               dis_ind[-1] - 80, [slopes[-1]/dt,
                               y_int[-1]]) > 0.5):
                return (slopes[-1] / dt, force_abs[-1]-zer[-1])
            #Otherwise return False values
            else:
                return (False, False)

    def _get_features(self):
        '''
        The _get_features method extracts values for slope and force from the
        current cycle.

        Output:
            A 2 element tuple containing a slope in pN/s and a force in pN
        '''

        if self._cycdat == None:
            return None

        #Find discontinuities in the cycle
        (force_abs, dis_ind, slopes, y_int, zer, smoothed_y) = self._data.fin_peaks(4.5, 161, 500, 1000)

        if self._data.tethered:
            return self._analyze_tethers(dis_ind, zer)

        else:
            return self._analyze_untethered(force_abs, dis_ind,
                                            slopes, y_int, zer)

    def _calc_rsq(self, graphchunk, fitwinst, fitwinen, line):
        '''
        The _calc_rsq method calculates the R squared valued of a fit

        Inputs:
            graphchunk: 2D array containing the data that was fit
            fitwinst: Starting index of the fit region to evaluate
            i: Ending index of the fit region to evaluate
            a: 2 element Numpy array containing the linear fit parameters
        '''

        if fitwinst > self._smoothparam / 2:
            pass
        else:
            fitwinst = int(self._smoothparam / 2)+1
        sm = self._smooth(graphchunk[1, :])[fitwinst:fitwinen]
        ybar = np.mean(sm)
        sstot = np.sum((sm - ybar) ** 2)
        ssres = np.sum((sm - line[0] * self._cycdat[0, fitwinst:fitwinen] -
                        line[1]) ** 2)
        return 1 - ssres/sstot

    def _update_confusion_mat(self, ids, confusion):
        '''
        The _update_confusion_mat method takes in a confusion matrix and
        updates some of its elements based on the current data cycle and the
        ground truth values stored in the SQL database.

        Inputs:
            ids: The file_id of the current cycle
            confusion: A confusion matrix (3x3 list of lists)

        Output:
            confusion: Updated confusion matrix (3x3 list of lists)
        '''

        #Get the ground truth for this cycle:
        s = ('SELECT tethered FROM afm_data.cell_data '
             'WHERE file_id=%s '
             'AND cycle=%s;')
        handres = []
        for rows in self._conn.execute(s%(ids, self._data.cyc)):
            handres.append(rows[0])

        #Update the confusion matrix:
        #Current cycle not in database (ground truth is junk)
        if len(handres) == 0:
            if self._data.tethered:
                confusion[1, 0] = confusion[1, 0] + 1
            else:
                confusion[2, 0] = confusion[2, 0] + 1
        #Current cycle is in database as a tether
        elif handres[0] == 'Y':
            if self._data.tethered:
                confusion[1, 1] = confusion[1, 1] + 1
            else:
                confusion[2, 1] = confusion[2, 1] + 1
        #Current cycle is in databse as untethered attachment
        else:
            if self._data.tethered:
                confusion[1, 2] = confusion[1, 2] + 1
            else:
                confusion[2, 2] = confusion[2, 2] + 1
        return confusion

    def _generate_all_features(self):
        '''
        The _generate_all_features method analyzes all the data for a given
        list of file ids, classifying each cycle and returning force and slope
        values for those that are classified as either tethers or untethered
        single bonds.
        '''

        confusion = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])

        #Cycling through all the values of _listofids
        for ids in self._listofids:
            #Current data object, corresponding to 1 raw data file
            self._data = EmittingSQLDat(ids, 4, self._conn, self.notifyProgress, self._prog)
            #Current cycle from the data object, 2d array of force, time data,
            #only retrieves cycles the classifier has marked as single bonds
            self._cycdat = self._data.get_next_cycle()
            #Determine whether the current raw data file has been analyzed and
            #how many total cycles it has
            s = ('SELECT analyzed, number_of_cycles FROM afm_data.filelist WHERE id = %s;')
            for rows in self._conn.execute(s%ids):
                analyzed = rows[0]
                cycs = int(rows[1])
            #If analyzed add nymber of cycles to the total number of analyzed
            #cycles
            if analyzed == 'Y':
                self._totalanalyzedcycs = self._totalanalyzedcycs + cycs

            while self._cycdat != None:
                try:
                    #Get slope and force for the current cycle
                    (slope, force) = self._get_features()

                    if slope == False:
                        pass
                    #Update lists of slopes and forces
                    else:
                        if self._data.tethered:
                            self._output['tetherslopes'].append(slope)
                            self._output['tetherforces'].append(force)
                        else:
                            self._output['untetheredslopes'].append(slope)
                            self._output['untetheredforces'].append(force)
                        #If this raw file has been analyzed, update the
                        #confusion matrix
                        if analyzed == 'Y':
                            confusion = self._update_confusion_mat(ids, confusion)
                except:
                    pass
                #Get the next cycle
                finally:
                    try:
                        self._cycdat = self._data.get_next_cycle()
                    except:
                        self._cycdat = self._data.get_next_cycle()
                #If abort signal received stop, return None
                if self._abortf == True:
                    self._output['tetherforces'] = None
                    return
            #When done with the _data object kill it's pool of workers, update
            #the progress value and then delete the object
            self._data.pool.terminate()
            self._prog = self._data.prog
            del self._data
        #Set the remaining values of the confusion matrix (any entry where
        #the classifier predicted junk still has not been accounted for)
        confusion[0, 0] = (self._totalanalyzedcycs -
                           len(self._output['truetetherforces']) -
                           len(self._output['trueuntetheredforces']) -
                           confusion[1, 0] - confusion[2, 0])
        confusion[0, 1] = (len(self._output['truetetherforces']) -
                           confusion[1, 1] - confusion[2, 1])
        confusion[0, 2] = (len(self._output['trueuntetheredforces']) -
                           confusion[1, 2] - confusion[2, 2])
        confusion = confusion.tolist()

        self._output['confusionmatrix'] = confusion
        return

#Application entry point
if __name__ == '__main__':
    freeze_support()
    import sys
    APP = QtGui.QApplication(sys.argv)
    APP_WIN = Viewer()
    APP_WIN.show()
    sys.exit(APP.exec_())
