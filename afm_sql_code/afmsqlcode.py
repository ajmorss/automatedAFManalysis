# -*- coding: utf-8 -*-
"""
Created on Mon Jun 16 11:10:09 2014

@author: Andrew
"""
import sqlalchemy as sqla
from PyQt4 import QtGui
from multiprocessing import Pool
import numpy as np
import csv
import cPickle as pickle
from scipy import signal

class SQLConnection(object):
    '''
    Connection to the SQL database where all the afm results/file info/file
    locations are stored.
    '''
    def __init__(self, username, password, host, port):
        try:
            #ENTER DATABASE INFO HERE
            connection = 'mysql+mysqldb://%s:%s@%s:%s/afm_data'
            connection = connection % (username, password, host, port)
            self.engine = sqla.create_engine(connection)
            self.metadata = sqla.MetaData()
            self.filelist = sqla.Table('filelist', self.metadata,
                                       sqla.Column('id', sqla.Integer, primary_key=True),
                                       sqla.Column('date', sqla.String(6)),
                                       sqla.Column('number', sqla.Integer),
                                       sqla.Column('full_file_name', sqla.String(300)),
                                       sqla.Column('k', sqla.Float),
                                       sqla.Column('OLS', sqla.Float),
                                       sqla.Column('loading_rate', sqla.Numeric(10,3)),
                                       sqla.Column('touch_time', sqla.Numeric(10,3)),
                                       sqla.Column('approach_rate', sqla.Numeric(10,3)),
                                       sqla.Column('touch_strength', sqla.Numeric(10,3)),
                                       sqla.Column('number_of_cycles',sqla.Integer),
                                       sqla.Column('pipette_obj', sqla.String(300)),
                                       sqla.Column('cantilever_obj', sqla.String(300)),
                                       sqla.Column('buffer', sqla.String(300)),
                                       sqla.Column('analyzed', sqla.String(1)),
                                       sqla.Column('pressure',sqla.Float),
                                       sqla.Column('cell_num',sqla.Integer))
            self.dat = sqla.Table('dat', self.metadata,
                                  sqla.Column('id', sqla.Integer, primary_key=True),
                                  sqla.Column('file_id',  None, sqla.ForeignKey('filelist.id')),
                                  sqla.Column('cycle', sqla.Integer),
                                  sqla.Column('piezo_state', sqla.Integer),
                                  sqla.Column('location', sqla.String(300)))
            self.cell_data = sqla.Table('cell_data', self.metadata,
                                        sqla.Column('id', sqla.Integer, primary_key=True),
                                        sqla.Column('file_id',  None, sqla.ForeignKey('filelist.id')),
                                        sqla.Column('cycle', sqla.Integer),
                                        sqla.Column('forces', sqla.Float),
                                        sqla.Column('slope', sqla.Float),
                                        sqla.Column('zero_position', sqla.Float),
                                        sqla.Column('start_time', sqla.Float),
                                        sqla.Column('cyt_break_time', sqla.Float),
                                        sqla.Column('end_time', sqla.Float),
                                        sqla.Column('integral', sqla.Float),
                                        sqla.Column('tethered', sqla.String(3)))
            self.metadata.create_all(self.engine)
            self.conn = self.engine.connect()
        except Exception as e:
            mb = QtGui.QMessageBox(parent=None)
            mb.setWindowTitle('ERR')
            mb.setText(str(e))
            mb.exec_()    

    def close(self):
        self.conn.close()
    
    def execute(self, s):
        return self.conn.execute(s)

class FeatureInfo(object):
    '''
    Class used by the _build_classification_features method of the SQLData class.
    '''
    def __init__(self, features):
        self.featuremeans = []
        self.featurestds = []
        for f in xrange(0, features.shape[1]):
            self.featuremeans.append(np.mean(features[:,f]))
            self.featurestds.append(np.std(features[:,f]))

def imp1(location):
    '''
    The imp1 function is the target function of the worker pools used by the
    SQLData class.
    
    Inputs:
        location: Location on the hard drive of the piezo state 3 data
    
    Outputs:
        chdat: The force vs time data of the given cycle.
    '''
    
    try:
        with open(location) as datain:
            #Reads data in
            datain = csv.reader(datain, delimiter=',')
            datax = []
            datay = []
            for rows in datain:
                datax.append(float(rows[0]))
                datay.append(float(rows[1]))
            data = np.vstack((datax, datay))
    
            #Generates a normalized rectangular window
            win = np.ones((int(21)))/21
    
            #Creates the smoothed data by convolution with the rectangular window
            tempx = np.convolve(data[0, :], win, 'valid')
            tempy = np.convolve(data[1, :], win, 'valid')
            dataav = np.vstack((tempx, tempy))
    
            #Flag to control while loop
            reading = True
            #Current read position
            index = 0
            #Runs through the smoothed cycle until it finds a value greater
            #than zero
            while reading:
                if dataav[1, index] > 0:
                    zero = index
                    reading = False
                #If for whatever reason no crossing is found, entire cycle is
                #kept
                if index == dataav[1, :].size - 1:
                    zero = 1
                    reading = False
                index = index + 1
    
            #Chops the data using the zero crossing point
            tempx = (data[0,
                     int(zero):data[0, :].size] -
                     data[0, int(zero)])
            tempy = data[1,
                    int(zero):data[1, :].size]
            chdat = np.vstack((tempx, tempy))
    except:
        chdat = None
    return chdat

class SQLData(object):
    '''
    The SQLData class retrieves force vs time information from the hard drive
    via the SQL Database.  It largely replaces the older Dat class for the
    purposes of interacting with forde_dat_improved.py.

    After initilized the class the class keeps a current position within the
    data set.  Calling get_next_cycle moves forward and returns the next cycle
    while get_previous_cycle moves back and returns the previous cycle.    
    
    The class uses the SQLConnection class to interface with the database and
    retrieve file locations.  The associated function imp1 assumes that the
    retrieved data is a comma delimited two column file, with the first column
    containing the time and the second the calibrated force.
    
    The SQL database is meant to be the afm_data database.  It should have a
    table called 'filelist' and another called 'dat'.
    
    Member Variables/Objects:
        Private:
            _classifier: Classifier object (scikit learn SVM) for automagic
                         selection
            _cycs: Total number of cycles in the selected data set
            _cyc_list: List of cycles that meet the selection criteria, note
                       that this list is built "as you go" by calling
                       get_next_cycle().
            _cyc_list_pos: Current position within the cycle list (for example
                           if the current cycle is 5 and the cycle list is 
                           [0, 2, 4, 5, 6, 7] then this is 3).
            _datares: The results of the import worker pools. Calling
                      self._datares[self.cyc].get() returns current cycle.
            _err_suppress: Flag to determine whether or not further output of 
                           cycle not found error message box is suppressed. 
            _imported_list: Imported cycle list used with selection method 3
            _method_flag: The selection method, takes the following values:
                          0: all cycles returned
                          1: thresholding max values
                          2: automagic selection
                          3: imported list
            _parent: Parent window that the cycle not found error box will
                     associate with
            _thresh_val: Value to threshold at if _method_flag is 1
            _win: Window size to smooth over before thresholding if 
                  _method_flag is 1
        
        Public:
            cyc: The current cycle
            fileid: The id of the current data sets raw file in the SQL 
                    database
            ps0_locations: A list of the locations of the piezo state 0 files
            ps1_locations: A list of the locations of the piezo state 1 files
            ps2_locations: A list of the locations of the piezo state 2 files
            ps3_locations: A list of the locations of the piezo state 3 files
            pool: Worker pool (multiprocessing package) for file imports
    
    Methods:
        Private:
            __init__: Initializes variables, retrieves file locations and
                      starts worker pools
            _bin_dat: Bins data
            _build_classification_features: Builds features used by the automagic
                                   classifier 
            _test_cyc: Tests to see whether the current cycle meets the
                       selection criteria
            
        Public:
            get_next_cycle: Returns next cycle that passes the selection 
                            criteria and sets the current cycle to that cycle
            get_previous_cycle: Gets the last cycle that passed the selection
                                criteria and sets the current cycle and cycle
                                list position to that cycle.
            fin_peaks: Finds discontinuities in the current cycle
    '''    
    
    def __init__(self, fileid, method_flag, conn, thresh_val=25, win=21, parent = None, cyclist = None):
        '''
        The __init__ method retrieves the file locations from the sql database
        and sends them to the worker pools to be imported. It also sets initial
        values for various member variables.
        '''

        #Initialize various variables
        self.cyc = -1
        self._cyc_list = []  
        self._cyc_list_pos = -1
        self.fileid = fileid   
        self.initial_zero = 0
        self._method_flag = method_flag
        self._imported_list = cyclist
        self._err_suppress = False        
        self._parent = parent
        self._thresh_val = thresh_val
        self._win = win
        if method_flag == 2:
            self._classifier = pickle.load( open( "machine_learning\\best_bond_classifier.p", "rb" ) )
        if method_flag == 4:
            self._classifier = pickle.load( open( "machine_learning\\best_bond_classifier.p", "rb" ) )
        self.tethered = None
        #Get file locations from the database
        s = sqla.select([conn.dat.c.location]).where((conn.filelist.c.id == fileid) & (conn.dat.c.file_id == conn.filelist.c.id) & (conn.dat.c.piezo_state == 0)).order_by(conn.dat.c.cycle)
        self.ps0_locations = []
        for rows in conn.execute(s):
            self.ps0_locations.append(rows[0])
        self.ps1_locations = []       
        s = sqla.select([conn.dat.c.location]).where((conn.filelist.c.id == fileid) & (conn.dat.c.file_id == conn.filelist.c.id) & (conn.dat.c.piezo_state == 1)).order_by(conn.dat.c.cycle)
        for rows in conn.execute(s):
            self.ps1_locations.append(rows[0])
        self.ps2_locations = []
        s = sqla.select([conn.dat.c.location]).where((conn.filelist.c.id == fileid) & (conn.dat.c.file_id == conn.filelist.c.id) & (conn.dat.c.piezo_state == 2)).order_by(conn.dat.c.cycle)
        for rows in conn.execute(s):
            self.ps2_locations.append(rows[0])
        self.ps3_locations = []
        s = sqla.select([conn.dat.c.location]).where((conn.filelist.c.id == fileid) & (conn.dat.c.file_id == conn.filelist.c.id) & (conn.dat.c.piezo_state == 3)).order_by(conn.dat.c.cycle)
        for rows in conn.execute(s):
            self.ps3_locations.append(rows[0])
        s = sqla.select([conn.filelist.c.loading_rate, conn.filelist.c.k]).where((conn.filelist.c.id == fileid))
        for rows in conn.execute(s):
            self.retrate = float(rows[0])/float(rows[1])

        #Get total number of cycle locations retrieved
        self._cycs = len(self.ps3_locations)

        #Send the locations to the worker pools to be imported
        self.pool = Pool(processes=10)
        self._datares=[]
        for i in xrange(0, self._cycs):
            self._datares.append(self.pool.apply_async(imp1, [self.ps3_locations[i]]))

    def get_next_cycle(self):
        '''
        The get_next_cycle method gets the next cycle that passes the cycle
        selection criteria. It also moves the current cycle to that cycle.
        
        Output:
            result: None if no cycle found, else a 2d numpy array with the
                    force vs time data
        '''
        
        #If we're not at the end of the cycle list, just move forward within
        #the list
        if self._cyc_list_pos + 1 < len(self._cyc_list):
            self._cyc_list_pos = self._cyc_list_pos + 1
            self.cyc = self._cyc_list[self._cyc_list_pos]
            self._initial_zero()
            result = self._datares[self.cyc].get()
        #If we still have more cycles before the end of the data set increment
        #cyc by 1 and test it, if it passes return it, if not then run the
        #method again
        elif self.cyc + 1 < self._cycs:
            self.cyc = self.cyc + 1
            #No result means error:
            if self._datares[self.cyc].get() != None:
                #If passes test append current cycle to the list, set the list
                #position to the end of the list and return the cycle information
                if self._test_cyc():
                    self._cyc_list.append(self.cyc)
                    self._cyc_list_pos = len(self._cyc_list) - 1
                    self._initial_zero()
                    result = self._datares[self.cyc].get()
                #Otherwise, we do reucursion
                else:
                    result = self.get_next_cycle()
            else:
                self._initial_zero()
                result = self.get_next_cycle()
        #If we're at the end of the data set there are no more cycles, return
        #None
        else:
            result = None
        return result

    def get_previous_cycle(self):
        '''
        The get_previous_cycle method goes to the previous cycle in _cyc_list.
        It also moves the current cycle and cycle list position to that cycle.

        Output:
            result: None if at beginning already, else a 2d numpy array with 
                    the force vs time data
        '''
        
        if self._cyc_list_pos - 1 > -1:
            self._cyc_list_pos = self._cyc_list_pos - 1
            self.cyc = self._cyc_list[self._cyc_list_pos]
            self._initial_zero()
            return self._datares[self.cyc].get()
        #If already at the beginning of the list return None
        else:
            return None

    def _initial_zero(self):
        fileloc = self.ps0_locations[self.cyc]
        forces = []
        with open(fileloc) as datin:
            datin = csv.reader(datin)
            for rows in datin:
                forces.append(float(rows[1]))
        self.initial_zero = sum(forces)/len(forces)
        
    def _test_cyc(self):
        '''
        The _test_cyc method determines whether or not the current cycle meets
        the selection criteria.
        
        Output:
            Boolean True for pass, False for fail
        '''
        
        #If method is "return all" then always passes
        if self._method_flag == 0:
            return True
        #Thresholding
        if self._method_flag == 1:
            self._initial_zero()
            if self.initial_zero < -40:
                return False
            #Creates a rectangular window
            window = np.ones((int(self._win)))/self._win
            #Creates the smoothed data by convolution with the rectangular window
            tempx = np.convolve(self._datares[self.cyc].get()[0, :], window, 'valid')
            tempy = np.convolve(self._datares[self.cyc].get()[1, :], window, 'valid')
            dataav = np.vstack((tempx, tempy))
            return np.max(dataav[1, :]) > self._thresh_val
        #Automagic
        if self._method_flag == 2 or self._method_flag == 4:
            self._initial_zero()
            if self.initial_zero < -40:
                return False
            #Get the features to pass to the classifier
            features = self._build_classification_features()
            features = np.reshape(features, (1,features.size))
            #Return classifier results
            if self._classifier.predict(features) == 0:
                return False
            if self._classifier.predict(features) == 1:
                self.tethered = True
                return True
            if self._classifier.predict(features) == 2:
                self.tethered = False
                return True

        #Imported list
        if self._method_flag == 3:
            if self.cyc in self._imported_list:
                return True
            else:
                return False

    def fin_peaks(self, det_thresh, win, sl_win, zer_len):
        '''
        fin_peaks locates discontinuities in a cycle, then uses them to
        generate values for rupture forces, slopes etc. It does this by looking
        for large values of the gradient/derivative of the cycle, where "large"
        is defined in terms of the standard deviation.

        Inputs:
            det_thresh: How big, in units of standard deviations, a derivative
                        has to be to be recognized as large (try starting with
                        4.5)
            win: Size of the window to smooth the data over before looking for
                 discontinuties

            sl_win: Size of the window over which to calculate the slope of the
                    linear fit of the data before a discontinuity

            zer_len: Size of the window over which to calculate the location of
                     "zero force" (the value of the foce immediately after the
                     break)

        Outputs:
            The method outputs a tuple of format:

            force_abs, dis_ind, slopes, y_int, zer, smoothed_y)

            Each member of the tuple is a numpy array which refers to
            all discontinuities detected by the method

            force_abs: The absolute value of the rupture forces as measured by
                       the AFM, relative to the zero defined for the cycle. The
                       location of the zero forces immediately after the breaks
                       (zer) must be subtracted off from this to get the true
                       value of the rupture force.  This force is calculated by
                       using a linear fit and finding the value of that fit at
                       the time of rupture.
            dis_ind: The locations of the indices of the discontinuties
            slopes: The slope of a linear fit of the data immediately before
                    the rupture, *this slope is in units of pN/indices (in
                    other words if your data was sampled at 10 hz then this
                    slope will have to be multiplied by 10 to put it in pN/s)*
            y_int: Y intercept value of the linear fit of the data immediately
                   before the rupture, in pN
            zer: The value of the force reading from the QPD immediately after
                 the rupture (zero location)
            smoothed_y: The smoothed y data used to find the discontinuties
        '''

        #This section smoothes the input data by convolution with a gaussian
        #window
        #Window Size
        gwin = win
        halfgwin = int(np.floor(gwin / 2))
        #Makes window
        window = (signal.gaussian(gwin, std=np.floor(gwin/6)) /
                 np.sum(signal.gaussian(gwin, std=np.floor(gwin/6))))
        #Smoothes force data
        smoothed_y = np.convolve(self._datares[self.cyc].get()[1, :], window, 'valid')
        #Pads data so it is the same length as the input data
        smoothed_y = np.hstack((np.ones((halfgwin)) * smoothed_y[0],
                                smoothed_y, np.ones((halfgwin)) *
                                smoothed_y[smoothed_y.size-1]))
        #Creates an 'x' array of index values of the force data
        smoothed_x = np.arange(smoothed_y.size)

        #Takes the gradient of the smoothed data
        grad_y = np.gradient(smoothed_y)

        #This section smoothes the gradient data with a hamming window
        #Creates the window
        window = signal.hamming(21)/np.sum(signal.hamming(21))
        halfwin = int(np.floor(21 / 2))
        #Smoothes the gradient
        grad_y_s = np.convolve(grad_y, window, 'valid')
        #Pads the smoothed data so it is the same length as the input data
        grad_y_s = np.hstack((np.zeros(halfwin), grad_y_s,
                            np.zeros(halfwin)))
        #Creates an 'x' array of the index values of the gradient data
        grad_x = np.arange(np.size(grad_y))

        #Creates a list of large gradients
        try:
            large_grads_x = grad_x[-grad_y_s > det_thresh * np.std(grad_y_s)]
            large_grads_y = grad_y_s[-grad_y_s > det_thresh * np.std(grad_y_s)]
        except ValueError:
            return (np.array([]), np.array([]), np.array([]), np.array([]),
                    np.array([]), np.array([]))

        if large_grads_y.size != 0:
            #First we find the approximate locations of discontinuities based
            #on the smoothed gradient list
            #If there is more than one entry in large_grads we must remove
            #duplicate entries for single discontinuities (strings of large
            #gradients must be reduced to single entries)
            if large_grads_y.size > 1:
                #Breaks large grads up into "connected blocks" (blocks of
                #large grads seperated by less than 10 indices).
                large_grad_diff = np.diff(large_grads_x)
                large_grad_diff_index = np.arange(large_grad_diff.size)
                #The + 1 adjusts for an offset from the diff function
                blocks = np.hstack((0, large_grad_diff_index[large_grad_diff >
                                    10] + 1, large_grads_y.size))
                #List of discontinuity locations
                dis_ind = []

                #Finds the location of the largest (most negative) gradient in
                #each block.  a is the location within each block of the
                #largest gradient.
                for index in xrange(0, blocks.size-1):
                    #If a block consists of only one entry, a=0
                    if blocks[index]+1 == blocks[index+1]:
                        max_ind = 0
                    #Else a is the location within the block of the largest
                    #gradient (most negative)
                    else:
                        max_ind = np.argmin(large_grads_y[blocks[index]:
                                            blocks[index+1]])
                    #Append the location of the largest gradient in the
                    #block to dis_ind (index at location within the block
                    #+ start of block)
                    dis_ind.append(large_grads_x[max_ind + blocks[index]])

            #If only 1 entry in large_grads, then dis_ind is simply the index
            #of that entry
            else:
                dis_ind = large_grads_x

            #Uses the unsmoothed gradient list to find exact location of
            #various discontinuities based on the approximate locations found
            #above
            for index in xrange(0, len(dis_ind)):
                disc = dis_ind[index]
                max_ind = np.argmin(grad_y[disc - 25:
                                    disc + 25])
                dis_ind[index] = max_ind + grad_x[disc - 25]

            #Various outputs
            slopes = []
            y_int = []
            zer = []
            force_abs = []
            #For each entry in dis_ind, finds slope, y_int, zer and force_abs
            #and appends them to the list
            for index in xrange(0, len(dis_ind)):
                disc = dis_ind[index]
                if disc > sl_win:
                    line_fit_x = smoothed_x[disc - sl_win:
                                            disc]
                    line_fit_y = np.array(self._datares[self.cyc].get()[1, disc -
                                sl_win:disc])
                    line = np.polyfit(line_fit_x, line_fit_y, 1)
                else:
                    line_fit_x = smoothed_x[0:disc]
                    line_fit_y = np.array(self._datares[self.cyc].get()[1, 0:disc])
                    line = np.polyfit(line_fit_x, line_fit_y, 1)

                slopes.append(line[0])
                y_int.append(line[1])
                force_abs.append(line[0] * disc + line[1])

                if dis_ind[index] + zer_len + 5 < len(self._datares[self.cyc].get()[1, :]):
                    zer.append(np.mean(self._datares[self.cyc].get()[1,
                                       disc:
                                       disc + zer_len]))
                else:
                    zer.append(np.mean(self._datares[self.cyc].get()[1,
                                       disc:
                                       len(self._datares[self.cyc].get()[1, :])]))

        #If large_grads is empty, return empty lists
        else:
            dis_ind = []
            force_abs = []
            slopes = []
            y_int = []
            zer = []

        #Convert lists to numpy arrays
        force_abs = np.array(force_abs)
        dis_ind = np.array(dis_ind)
        slopes = np.array(slopes)
        y_int = np.array(y_int)
        zer = np.array(zer)

        return force_abs, dis_ind, slopes, y_int, zer, smoothed_y

    def _bindat(self, data, newdatlen):
        '''
        The _bindat method bins data.
        
        Inputs:
            data: Data to be binned, 1d numpy array
            newdatlen: Desired final data length
        
        Outputs:
            datbinned: Binned data, a 1d numpy array
        '''

        datlen = data.size
        if datlen>=newdatlen:
            bins = np.linspace(0, datlen-1, newdatlen+1)
            datbinned = np.zeros(newdatlen)
            for x in xrange(0, datbinned.size):
                datbinned[int(x)] = np.mean(data[int(bins[x]):int(bins[x+1])])
            
        else:
            datbinned = np.hstack((data,np.zeros(newdatlen-datlen,)))
            
        return datbinned
        
    def _build_classification_features(self):
        '''
        The _build_classification_features method creates the features that are passed
        to the automagic classifier.
        
        Output:
            features: 1d numpy array of features
        '''
        
        if self._method_flag ==2:
            featureinfo = pickle.load( open( "machine_learning\\feature_info_for_bond_classifier.p", "rb" ) )   
        elif self._method_flag == 4:
            featureinfo = pickle.load( open( "machine_learning\\feature_info_for_bond_classifier.p", "rb" ) )   
        peaks = self.fin_peaks(4.5, 201, 5000, 50)
        if peaks[0].size > 0:
            finbreak = peaks[1][peaks[1].size - 1]
            endp = int(finbreak*1.25)
        else:
            endp = self._datares[self.cyc].get()[1, :].size - 1
        samples_dat = self._bindat(self._datares[self.cyc].get()[1, 0:endp], 512)
        ffts = (np.fft.fft(samples_dat[:]))*np.conj((np.fft.fft(samples_dat[:])))
        ffts = ffts[0:ffts.size/2]
        ffts = np.real(ffts)
        samples_dat = self._bindat(samples_dat, samples_dat.size/2)
        ffts = self._bindat(ffts, ffts.size/2)
    
        grad_y = np.gradient(samples_dat)
        if self._method_flag == 2:
            large_grads = grad_y[-grad_y > 2.5 * np.std(grad_y)]
            disconts = large_grads.size  
        elif self._method_flag == 4:
            disconts = peaks[0].size  

        largestdrop = np.min(grad_y)
        slopes = np.max(grad_y)
        
        integrals = np.trapz(samples_dat)
        if self._method_flag == 4:
            integrals = 0
#            dxp=self._datares[self.cyc].get()[0, 1]-self._datares[self.cyc].get()[0, 0]
#            integrals = np.trapz(self._datares[self.cyc].get()[1, 0:endp], dx=dxp)    

        means = np.mean(samples_dat)

        maxs = np.max(samples_dat)
    
        fft_peak_no = ffts[ffts>2.5*np.std(ffts)].size
        
        binned = self._bindat(samples_dat, 16)

        features = [disconts, means, maxs, slopes, largestdrop, integrals, fft_peak_no]
        features = np.transpose(np.array(features))
        features = np.hstack((features, binned))
        
        if self._method_flag == 4:
            ffts = self._bindat(ffts, 16)
            features = np.hstack((features, ffts))

        for f in xrange(0, features.size):
            if featureinfo.featurestds[f] == 0:
                features[f] = (features[f]-featureinfo.featuremeans[f])/1
            else:
                features[f] = (features[f]-featureinfo.featuremeans[f])/featureinfo.featurestds[f]
        return features
