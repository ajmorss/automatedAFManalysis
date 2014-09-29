# -*- coding: utf-8 -*-
"""
Created on Tue Jun 10 15:39:35 2014

@author: Andrew
"""

from automatedAFManalysis.dataviewer.legacy_interface import Ui_MainWindow
from PyQt4 import QtGui, QtCore
import numpy as np
import automatedAFManalysis.dataviewer.connectdialog as cd
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
import csv
from automatedAFManalysis.afm_sql_code.afmsqlcode import SQLConnection

class MainWindow(QtGui.QMainWindow, Ui_MainWindow):

    def __init__(self, parent=None, f=QtCore.Qt.WindowFlags()):
        QtGui.QMainWindow.__init__(self, parent, f)
        self.setupUi(self)

class ConnectDialog(QtGui.QDialog, cd.Ui_Dialog):
    '''
    This class defines the connect dialog box.
    '''
    def __init__(self,parent=None):
        QtGui.QDialog.__init__(self,parent)
        self.setupUi(self)
        self.okB.clicked.connect(self._ok_callback)
        self.setWindowTitle('Connect To SQL Database')
        self.exec_()
    
    def _ok_callback(self):
        self.username = self.usernameLE.text()
        self.password = self.passwordLE.text()
        self.host = self.hostnameLE.text()
        self.port = self.portLE.text()
        self.close()

class MyMplCanvas(FigureCanvas):
    '''
    Matplotlib canvas object
    '''

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        '''
        Sets up the canvas, creates a figure, adds subplots, sets various
        parameters.
        '''

        fig = Figure(figsize=(width, height), dpi=dpi)
        fig.patch.set_alpha(0)
        self.histogram = fig.add_subplot(211)
        self.scatterplot = fig.add_subplot(212)
        fig.subplots_adjust(left=0.18, top=0.95, bottom=0.07, hspace=.5)
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self,
                                   QtGui.QSizePolicy.Expanding,
                                   QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

class InputWindow(MainWindow):
    '''
    The window for the dataviewer

    Member Variables/Objects:

        Private:
            _cell_flag: Flag that indicates bead/cell data, True for cells
            _tethers: Flag that indicates tethered/untethred data, True for
                      tethers
            _forces_arr_list: List of force arrays (1d numpy), used in plotting
            _slopes_arr_list: List of slope arrays (1d numpy), used in plotting
            _database: The Database class object for the application
            _pipette_object_data: Data for the pipette_object tableview
            _cantilever_object_data: Data for the cantilever_object tableview
            _buffer_data: Data for the buffer tableview
            _date_data: Data for the date tableview
            _loading_rate_data: Data for the loading rate tableview

    Methods:

        Private:
            __init__: Builds window, sets initial values and builds a plot
            _bind_callbacks: Binds the callbacks for the UI
            _init_tables: Sets up the tables
            _cell_radio_onclick: Callback for the cell radio button
            _bead_radio_onclick: Callback for the bead radio button
            _radio_onclick: Callback for tethered/untethered radio buttons
            _hold_onclick: Callback for the hold checkbox
            _plot_callback: Builds plots
            _onclick: Method used in various callbacks for user interface
                      actions
            _update_tables: Updates the tables
            _get_selected: Get user selected values from the tables
            _get_rows: Get user selected rows from the tables
            _save_callback: Callback for the save menu action
            _exit_callback: Exit menu action callback
    '''

    def __init__(self, parent=None):
        '''
        The __init__ method builds the window, sets various initial values
        and builds a first plot.
        '''

        self._cell_flag = True
        self._tethers = True
        self.pipette_object_data = None
        self.cantilever_object_data = None
        self.buffer_data = None
        self.date_data = None
        self.loading_rate_data = None
        self.forces_arr_list = []
        self.slopes_arr_list = []
        QtGui.QDialog.__init__(self, parent)
        self.setupUi(self)
        self.matplot = MyMplCanvas(self.matplotframe,
                                   width=5, height=8, dpi=100)
        NavigationToolbar(self.matplot, self.toolframe, coordinates=True)
        self._database = Database()
        self._bind_callbacks()
        self._init_tables()
        self._plot_callback()

    def _connect_callback(self):
        connectdialog = ConnectDialog()
        self._conn = SQLConnection(connectdialog.username, connectdialog.password, connectdialog.host, connectdialog.port)
        self._database = Database(self._conn)
        self._init_tables()

    def _bind_callbacks(self):
        '''
        The _bind_callbacks method...binds the callbacks
        '''

        QtCore.QObject.connect(self.tethered, QtCore.SIGNAL('clicked()'),
                               self._radio_onclick)
        QtCore.QObject.connect(self.holdbox, QtCore.SIGNAL('clicked()'),
                               self._hold_onclick)
        QtCore.QObject.connect(self.celldat, QtCore.SIGNAL('clicked()'),
                               self._cell_radio_onclick)
        QtCore.QObject.connect(self.beaddat, QtCore.SIGNAL('clicked()'),
                               self._bead_radio_onclick)
        QtCore.QObject.connect(self.untethered, QtCore.SIGNAL('clicked()'),
                               self._radio_onclick)
        QtCore.QObject.connect(self.plot_button, QtCore.SIGNAL('clicked()'),
                               self._plot_callback)
        QtCore.QObject.connect(self.actionExit, QtCore.SIGNAL('triggered()'),
                               self._exit_callback)
        QtCore.QObject.connect(self.actionSave, QtCore.SIGNAL('triggered()'),
                               self._save_callback)
        QtCore.QObject.connect(self.actionConnect, QtCore.SIGNAL('triggered()'),
                               self._connect_callback)
        tablist = [self.pipette_object, self.cantilever_object, self.date,
                   self.loading_rate, self.buffer]
        for table in tablist:
            table.clicked.connect(self._onclick)

    def _init_tables(self):
        '''
        The _init_tables method initializes the tables when the program is
        opened.
        '''

        #Start by setting all table data to 'ANY'
        self.pipette_object_data = [['ANY']]
        self.cantilever_object_data = [['ANY']]
        self.buffer_data = [['ANY']]
        self.date_data = [['ANY']]
        self.loading_rate_data = [['ANY']]
        self._update_tables()
        #Make sure the 'ANY' option is actually selected for each table
        tablist = [self.pipette_object, self.cantilever_object, self.date,
                   self.loading_rate, self.buffer]
        for table in tablist:
            table.selectRow(0)
            table.verticalHeader().setVisible(False)
            table.horizontalHeader().setVisible(False)
        #Set initial attachments and total cycle number labels
        cycnumbers = self._database.update_list(self._get_selected(),
                                                     self._cell_flag,
                                                     self._tethers)
        self.attachments.setText(str(cycnumbers[0]))
        self.total_cycles.setText(str(cycnumbers[1]))
        #Get all the valid choices for the various tables
        vld_params = self._database.get_valid_params(self._get_selected(),
                                                   self._cell_flag,
                                                   self._tethers)
        (self.pipette_object_data, self.cantilever_object_data,
         self.buffer_data, self.date_data, self.loading_rate_data) = vld_params
        #Run _onclick and resize the columns
        self._onclick()
        for table in tablist:
            table.resizeColumnsToContents()

    def _cell_radio_onclick(self):
        '''
        The _cell_radio_onclick method defines the callback for the cells radio
        button.  It enables the tethered and untethered radio buttons, sets
        the cell_flag to True and runs the _onclick method.
        '''

        self.untethered.setEnabled(True)
        self.tethered.setEnabled(True)
        self.cell_flag = True
        self._onclick()

    def _bead_radio_onclick(self):
        '''
        The _bead_radio_onclick method defines the callback for the beads radio
        button.  It disables the tethered and untethered radio buttons, sets
        the cell_flag to False and runs the _onclick method.
        '''

        self.untethered.setEnabled(False)
        self.tethered.setEnabled(False)
        self.cell_flag = False
        self._onclick()

    def _radio_onclick(self):
        '''
        The _radio_onclick method defines the callback for the tethered
        and untethered radio buttons
        '''

        if self.tethered.isChecked():
            self.tethers = True
        elif self.untethered.isChecked():
            self.tethers = False
        self._onclick()

    def _hold_onclick(self):
        '''
        The _hold_onclick callback defines the callback for the hold check box.
        If the box is unchecked it clears the force and slope array lists.
        '''

        if self.holdbox.isChecked():
            pass
        else:
            self.forces_arr_list = []
            self.slopes_arr_list = []

    def _plot_callback(self):
        '''
        The _plot_callback method builds the plots for the matplotlib widget.
        '''

        hist = self.matplot.histogram
        scpl = self.matplot.scatterplot
        #Clear previous plots
        scpl.axes.cla()
        hist.axes.cla()
        #If hold check box isn't checked, clear the force and slopes array
        #lists
        if self.holdbox.isChecked():
            pass
        else:
            self.forces_arr_list = []
            self.slopes_arr_list = []
        #Get forces and slopes for the current user selected paramaters
        forces = [x[3] for x in self.database.list_of_features]
        forces = np.array(forces)
        slopes = [x[4] for x in self.database.list_of_features]
        slopes = np.array(slopes)
        self.forces_arr_list.append(forces)
        self.slopes_arr_list.append(slopes)
        #If holdbox checked plot everything, set alpha = 0.3 for histogram
        if self.holdbox.isChecked():
            scpl.axes.hold(True)
            hist.axes.hold(True)
            scpl.axes.set_xscale('log')
            for i, forces in enumerate(self.forces_arr_list):
                scpl.axes.plot(self.slopes_arr_list[i],
                          forces, '.', ms=3)
                hist.axes.hist(forces,
                               bins=range(0, int(np.max(forces))+5, 5),
                               alpha=0.3, normed=True)
        #Else just plot current stuff, alpha = 1
        else:
            scpl.axes.plot(slopes, forces, '.', ms=3)
            scpl.axes.set_xscale('log')
            hist.axes.hist(forces, normed=True,
                      bins=range(0, int(np.max(forces))+5, 5))
        #Set titles appropriately
        if self.tethered.isChecked() and self.cell_flag:
            hist.set_title('Histogram of Crossover Forces', fontsize=14)
            scpl.set_title(('Scatter Plot of Crossover Forces \n'
                            'vs. Loading Rate'), fontsize=14)
        else:
            hist.set_title('Histogram of Rupture Forces', fontsize=14)
            scpl.set_title(('Scatter Plot of Rupture Forces \n'
                            'vs. Loading Rate'), fontsize=14)
        #Set axis parameters and draw it
        hist.set_xlabel('Force (pN)', fontsize=12)
        hist.set_ylabel('Frequency', fontsize=12)
        scpl.set_xlabel('Loading Rate (pN/s)', fontsize=12)
        scpl.set_ylabel('Force (pN)', fontsize=12)
        hist.tick_params(axis='both', which='major', labelsize=10)
        scpl.tick_params(axis='both', which='major', labelsize=10)
        hist.set_xlim(0,100)
        scpl.set_xlim(100,10000)
        scpl.set_ylim(0,100)
        self.matplot.draw()

    def _onclick(self):
        '''
        The _onclick methed is used as a response to various user clicks.
        It updates updates the tables with valid parameter options for the
        selected rows and updates the attachment and total cycle text labels.
        '''

        #Get valid parameter options for user selected rows:
        selected = self._get_selected()
        vld_params = self._database.get_valid_params(selected,
                                             self._cell_flag,
                                             self._tethers)
        (self.pipette_object_data, self.cantilever_object_data,
         self.buffer_data, self.date_data, self.loading_rate_data) = vld_params
        #Update the tables with the new parameter options:
        self._update_tables()
        #Restore user row selections:
        tablist = [self.pipette_object, self.cantilever_object, self.buffer,
                   self.date, self.loading_rate]
        datlist = [self.pipette_object_data, self.cantilever_object_data,
                   self.buffer_data, self.date_data, self.loading_rate_data]
        for i in xrange(0, len(tablist)):
            for j in selected[i][0]:
                tablist[i].setSelectionMode(QtGui.QAbstractItemView.MultiSelection)
                try:
                    row = datlist[i][0].index(j)
                    tablist[i].selectRow(row)
                except ValueError:
                    pass
                tablist[i].setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
        #Update the attachment and total cycle text labels:
        cycnumbers = self._database.update_list(self._get_selected(),
                                                     self._cell_flag,
                                                     self._tethers)
        self.attachments.setText(str(cycnumbers[0]))
        self.total_cycles.setText(str(cycnumbers[1]))

    def _update_tables(self):
        '''
        The _update_tables method updates the various tables using the current
        valid parameters for the user selected rows.
        '''

        tablist = [self.pipette_object, self.cantilever_object, self.buffer,
                   self.date, self.loading_rate]
        datlist = [self.pipette_object_data, self.cantilever_object_data,
                   self.buffer_data, self.date_data, self.loading_rate_data]

        for i in xrange(0, len(tablist)):
            tabmodel = Mytablemodel()
            tablist[i].setModel(tabmodel)
            tabmodel.update_table(datlist[i])

    def _get_selected(self):
        '''
        The _get_selected method gets the user selected values for the various
        parameters.

        Output:
            selected: User selected values for the various tables (tuple of
                      lists)
                      (pipette table values, cantilever object values,
                      buffer values, date values, loading rate values)
        '''

        #Get user selected rows:
        rows = self._get_rows()
        #Use those rows to get the values associated with them:
        datlist = [self.pipette_object_data, self.cantilever_object_data,
                   self.buffer_data, self.date_data, self.loading_rate_data]
        selected = []
        #Iterate over each table:
        for i in xrange(0, len(rows)):
            temp = []
            #Iterate over each selected row within a table:
            for j in rows[i]:
                temp.append(datlist[i][0][j])
            temp = [temp]
            selected.append(temp)
        selected = tuple(selected)

        return selected

    def _get_rows(self):
        '''
        The _get_rows method gets the user selected rows for the various
        tables.

        Output:
            table_rows: A tuple of lists containing all the user selected rows.
                        (pipette table rows, cantilever object rows, buffer
                         rows, date rows, loading rate rows)
        '''

        table_rows = []
        tablist = [self.pipette_object, self.cantilever_object, self.buffer,
                   self.date, self.loading_rate]
        #Iterate over tables
        for tables in tablist:
            temp = []
            #Iterate over selected rows within tables
            for sel_rows in tables.selectionModel().selectedRows():
                temp.append(sel_rows.row())
            table_rows.append(temp)
        table_rows = tuple(table_rows)

        return table_rows

    def _save_callback(self):
        '''
        The _save_callback method is the callback for the save menu action, it
        runs the database save_data method.
        '''
        self._database.save_data(self, self.cell_flag)

    def _exit_callback(self):
        '''
        The _exit_callback is the callback for the exit menu option.
        '''

        self._database.close()
        self.close()

class Mytablemodel(QtCore.QAbstractTableModel):
    '''
    The table model used for bead data.
    '''
    def __init__(self):
        QtCore.QAbstractTableModel.__init__(self)
        self.initdata = []
        self.rows = 0
        self.columns = 1

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

    def setdat(self, data):
        self.rows = len(data[0])
        self.initdata = data

    def update_table(self, dat):
        '''
        The update_table method updates the table with tether data.

        Inputs:
            dat: A list of lists to input to the table.
        '''

        self.setdat(dat)

class Database(object):
    '''
    The Database class defines all interactions with the SQL database.

    Member Variables/Objects:

        Private:
            _conn: The SQLConnection used by the class

        Public:
            list_of_features: Currently extracted features from the database.
                              A list of lists, where each entry is the features
                              for a given cycle.

    Methods:

        Private:
            __init__: Establishes database connection
            _build_conditions: Builds the conditions/constraints used by the
                               update_list and get_valid_params methods

        Public:
            update_list: Updates the list_of_features object for a
                                 given set of constraints, returns numbers
                                 for total attachments and total cycles
            get_valid_params: Returns the valid parameter choices given a set of
                       chosen parameter
            save_data: Saves list_of_features to a CSV
            close: Closes the SQLConnection object
    '''

    def __init__(self, conn):
        '''
        The __init__ method establishes a connection to the database
        '''

        self._conn = conn
        self.list_of_features = []

    def update_list(self, tuple_in, cell_flag, tetherflag):
        '''
        The update_list method updates the list_of_features.

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
                The length of the list_of_features variable (total attachments)
                tot: The total number of cycles meeting the user supplied
                     selection criteria
        '''

        #Get features:
        #SQL query statement made of two parts, a and b, where a is selection
        #information and b is constraints/conditions
        if cell_flag:
            select = ('SELECT full_file_name, number, date, forces, slope, '
                 'start_time, cyt_break_time, end_time, integral, '
                 'cell_num FROM afm_data.filelist')
        else:
            select = ('SELECT full_file_name, number, date, forces, slope, time '
                 'FROM afm_data.filelist')
        conditions = self._build_conditions(tuple_in, 9, cell_flag, tetherflag)
        statement = ''.join([select, conditions])
        result = self._conn.execute(statement)
        self.list_of_features = []
        for cycle in result:
            self.list_of_features.append(cycle)
        #Get total number of cycles (regardless of attachment data):
        #SQL query statement made of two parts, a and b, where a is selection
        #information and b is constraints/conditions
        if cell_flag:
            select = ('SELECT number_of_cycles FROM afm_data.filelist '
                 'INNER JOIN '
                 '(SELECT DISTINCT file_id FROM afm_data.cell_data) S '
                 'ON afm_data.filelist.id = s.file_id')
        else:
            select = ('SELECT number_of_cycles FROM afm_data.filelist '
                 'INNER JOIN '
                 '(SELECT DISTINCT file_id FROM afm_data.bead_data) S '
                 'ON afm_data.filelist.id = s.file_id')
        conditions = self._build_conditions(tuple_in, 9, cell_flag,
                                   tetherflag, buildjoin=False)
        statement = ''.join([select, conditions])
        tot = 0
        self._conn.execute(statement)
        for cycles in self._conn.execute(statement):
            tot = tot + cycles[0]
        #Return number of attachments and cycles
        return (len(self.list_of_features), tot)

    def _build_conditions(self, tuple_in, knockout,
                          cell_flag, tetherflag, buildjoin=True):
        '''
        The _build_conditions method builds the constraints/conditions used by
        the update_list and get_valid_params methods

        Inputs:
            tuple_in: Tuple of lists containing the various constraints for
                      each selection parameter
                      (pipette object list, cantilever object list,
                       buffer list, date list, loading rate list)
            knockout: Integer that indicates a parameter from the tuple_in over
                      which we don't wish to apply constraints, used by
                      get_valid_params
            cell_flag: True for cells, False for beads
            tetherflag: True for tethered, false for untethered

        Parameters:
            buildjoin: True builds a join statement and a tether constraint,
                       False does not

        Output:
            conditions: The output constraints, string SQL query fragment
        '''

        def build_single_condition(column, tuple_in, whereandflag):
            '''
            build_single_condition builds the constraints for a single paramter
            for the _build_conditions method

            Inputs:
                column: Which element of the names table we are building a
                        constraint for where names is a list of the various
                        selection parameters
                        0 = pipette object
                        1 = cantilever object
                        2 = buffer
                        3 = date
                        4 = loading rate
                tuple_in: Tuple of lists containing the various constraints for
                          each selection parameter
                          (pipette object list, cantilever object list,
                           buffer list, date list, loading rate list)
                whereandflag: Whether we need a WHERE or an AND (need where
                              for first parameter in a query, AND otherwise).
                              True for WHERE.

            Outputs:
                condition: String with all the constraints for the given
                           parameter
            '''

            names = ['pipette_obj', 'cantilever_obj',
                     'buffer', 'date', 'loading_rate']
            condition = ''
            if whereandflag:
                whereand = ' WHERE ('
            else:
                whereand = ' AND ('

            #Loop over all the constraints for the given parameter
            for cond in xrange(0, len(tuple_in[column][0])):
                if cond == 0:
                    #Special case for loading rate (not a string)
                    if names[column] == 'loading_rate':
                        temp = ''.join([whereand, names[column], ' = ',
                                        str(tuple_in[column][0][cond])])
                    else:
                        temp = ''.join([whereand, names[column], ' = \'',
                                        str(tuple_in[column][0][cond]), '\''])
                else:
                    #Special case for loading rate (not a string)
                    if names[column] == 'loading_rate':
                        temp = ''.join([' OR ', names[column], ' = ',
                                        str(tuple_in[column][0][cond])])
                    else:
                        temp = ''.join([' OR ', names[column], ' = \'',
                                        str(tuple_in[column][0][cond]), '\''])
                #temp is a single constraint, all constraints build the
                #statement for this parameter
                condition = ''.join([condition, temp])
            #All the constraints added, toss a parenthesis at the end
            condition = ''.join([condition, ')'])
            return condition

        def build_join(conditions, cell_flag, tetherflag):
            '''
            Builds the JOIN statement for the _build_conditions method

            Inputs:
                conditions: Conditions in
                cell_flag: True for cells, False for beads
                tetherflag: True for tethered, false for untethered

            Outputs
                conditions: Conditions out, with an INNER JOIN
            '''

            #Need to add a constraint on tethers if it's cell data
            if cell_flag:
                #If no conditions then we need a WHERE
                if len(conditions) == 0:
                    whereand = ' WHERE '
                else:
                    whereand = ' AND '
                #Tethers Y or N
                if tetherflag:
                    yorn = '\'Y\''
                else:
                    yorn = '\'N\''
                #Set the table name for the INNER JOIN and the tether
                #constraint
                jointab = 'cell_data'
                cellbit = ''.join([' ', whereand, 'tethered = ', yorn])
            #If not cells, use the bead_data table and set the tether
            #constraint to ''
            else:
                jointab = 'bead_data'
                cellbit = ''
            #Build final conditions out and return
            conditions = ''.join([' INNER JOIN afm_data.', jointab,
                                  ' ON filelist.id=', jointab, '.file_id',
                                  conditions, cellbit, ';'])
            return conditions

        #METHOD ENTRY POINT
        conditions = ''
        #Loop over the parameters
        for column in xrange(0, len(tuple_in)):
            condition = ''
            #No constraints on the knockout
            if column == knockout:
                pass
            #No conditions yet need a WHERE
            elif len(conditions) == len(''):
                #No constraints on 'ANY'
                if 'ANY' in tuple_in[column][0]:
                    pass
                #This shouldn't happen but just in case
                elif len(tuple_in[column][0]) == 0:
                    pass
                else:
                    condition = build_single_condition(column, tuple_in, True)
            #Already conditions need an AND
            else:
                if 'ANY' in tuple_in[column][0]:
                    pass
                elif len(tuple_in[column][0]) == 0:
                    pass
                else:
                    condition = build_single_condition(column, tuple_in, False)
            conditions = ''.join([conditions, condition])
        #If we need an INNER JOIN on fileid and a constraint on tethers:
        if buildjoin:
            conditions = build_join(conditions, cell_flag, tetherflag)

        return conditions

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
            stuff: List of valid parameters
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
                                                cell_flag, tetherflag)
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
        for params in xrange(len(valid_params[4][0])):
            if valid_params[4][0][params] == 'ANY':
                pass
            else:
                valid_params[4][0][params] = float(str(valid_params[4][0][params]))

        return valid_params

    def save_data(self, viewer, cell_flag):
        '''
        The save_data method saves the extracted database features to a CSV
        file

        Inputs:
            viewer: The data viewer window
            cell_flag: True for cells, False for beads
        '''

        #Prompt user for save location:
        fname = QtGui.QFileDialog.getSaveFileName(viewer, 'Open File',
                                                  'c:\\data', '*.txt')
        fname = str(fname)
        #If save location provided, save
        if fname != '':
            with open(fname, 'wb') as dataout:
                dataout = csv.writer(dataout, delimiter=',')
                #Write cell data:
                if cell_flag:
                    dataout.writerow(['full_file_name', 'number', 'date',
                                      'forces', 'slope', 'start_time',
                                      'cyt_break_time', 'end_time', 'integral',
                                      'cell_num'])
                    for cycle in self.list_of_features:
                        dataout.writerow(cycle)
                #Write bead data:
                else:
                    dataout.writerow(['full_file_name', 'number', 'date',
                                      'forces', 'slope', 'time'])
                    for cycle in self.list_of_features:
                        dataout.writerow(cycle)

    def close(self):
        '''
        The close method closes the connection to the database.
        '''

        self._conn.close()

if __name__ == '__main__':
    import sys
    APP = QtGui.QApplication(sys.argv)
    APP_WIN = InputWindow()
    APP_WIN.show()
    sys.exit(APP.exec_())
