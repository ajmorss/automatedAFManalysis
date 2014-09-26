# -*- coding: utf-8 -*-

from automatedAFManalysis.afm_sql_code.afmsqlcode import SQLConnection, SQLData
import numpy as np
import multiprocessing
import csv

class SQLDataForFeatures(SQLData):
    '''
    The SQLDataForFeatures class is a modified version of the SQLData class
    that facilitates the creation and saving of a classified feature set from
    AFM force time series for the purposes of supervised learning.
    
    Methods:
        Public:
            build_write_all_features: Builds classification features for all 
                                      cycles in the given raw data file then 
                                      writes them to tethersamples.txt
        Private:
            _build_classification_features: Build the classification features
                                            for the current cycle
    '''

    def build_write_all_features(self):
        '''
        Builds classification features for all cycles in the given raw data 
        file then writes them to tethersamples.txt
        '''

        while self.get_next_cycle() != None:
            #Get the features for the current cycle
            features = self._build_classification_features()
            #Get the ground truth classification
            conn = SQLConnection()
            s = ('SELECT tethered FROM cell_data '
                 'WHERE file_id=%s '
                 'AND cycle=%s;')
            tethered = []
            for rows in conn.execute((s % (self.fileid, self.cyc))):
                tethered.append(rows[0])
            if len(tethered) == 0:
                classification = 0
            elif tethered[0] == 'Y':
                classification = 1
            elif tethered[0] == 'N':
                classification = 2
            #Tack the classification and the unique id of the cycle (so we can
            #retrieve it later if we want)
            features = np.hstack((features, classification))
            s = ('select id from dat '
                 'WHERE file_id=%s '
                 'AND cycle=%s '
                 'AND piezo_state=3;')
            datid = []
            for rows in conn.execute((s % (self.fileid, self.cyc))):
                datid.append(int(rows[0]))            
            features = np.hstack((features, datid[0]))
            conn.close()
            #Write the cycles features to file
            features = features.tolist()
            with open('bondsamples.txt','ab') as datout:
                datout=csv.writer(datout)
                datout.writerow(features)

    def _build_classification_features(self):
        '''
        The _build_classification_features method creates the features that are
        passed to the classifier.
        
        Output:
            features: 1d numpy array of features for the current cycle
        '''

        #Preprocesses data (shortens it, bins it so the length is 512)       
        peaks = self.fin_peaks(4.5, 201, 5000, 50)
        if peaks[0].size > 0:
            finbreak = peaks[1][peaks[1].size - 1]
            endp = int(finbreak*1.25)
        else:
            endp = self._datares[self.cyc].get()[1, :].size - 1
        samples_dat = self._bindat(self._datares[self.cyc].get()[1, 0:endp], 512)
        #Get FFT, square it, take first half
        ffts = (np.fft.fft(samples_dat[:]))*np.conj((np.fft.fft(samples_dat[:])))
        ffts = ffts[0:ffts.size/2]
        ffts = np.real(ffts)
        #Bin the data a little more
        samples_dat = self._bindat(samples_dat, samples_dat.size/2)
        ffts = self._bindat(ffts, ffts.size/2)

        #First two features: largest discontinuity (drop/negative slope) and
        #largest positive slope
        grad_y = np.gradient(samples_dat)    
        largestdrop = np.min(grad_y)
        slopes = np.max(grad_y)

        #Next feature, total discontinuities
        disconts = peaks[0].size

        #TO DO: Rerun and retrain with actual integral values, currently all 0
        #dxp=self._datares[self.cyc].get()[0, 1]-self._datares[self.cyc].get()[0, 0]
        #integrals = np.trapz(self._datares[self.cyc].get()[1, 0:endp], dx=dxp)
        integrals = 0

        #Next features, the mean force of the data and the max of the data
        means = np.mean(samples_dat)
        maxs = np.max(samples_dat)

        #Next feature, the number of peaks in the power spectrum    
        fft_peak_no = ffts[ffts>2.5*np.std(ffts)].size

        #The last features are the the data itself and the power spectrum,
        #binned to a final length of 16 points each        
        binned = self._bindat(samples_dat, 16)
        ffts = self._bindat(ffts, 16)

        #Build and return features        
        features = [disconts, means, maxs, slopes, largestdrop, integrals, fft_peak_no]
        features = np.transpose(np.array(features))
        features = np.hstack((features, binned))
        features = np.hstack((features, ffts))
        return features

def write_all_features(idnum):
    '''
    Worker task, writes all the features for a given raw data file
    '''

    dataobj = SQLDataForFeatures(idnum, 0)
    dataobj.build_write_all_features()
    dataobj.pool.terminate()
    del(dataobj)
    return True


if __name__ == '__main__':
    #Get the list of all raw files that training/testing data will be pulled from
    conn = SQLConnection()
    s = ('SELECT id FROM afm_data.filelist '
         'WHERE id>374 '
         'AND id<614 '
         'AND analyzed=\'Y\' '
         'AND pipette_obj!=\'ECADBEAD (SNAP)\';')
    listofids=[]
    for rows in conn.execute(s):
        listofids.append(int(rows[0]))
    conn.close()
    print len(listofids)      

    #Break the complete list of file ids up into chunks    
    chunksize = float(5)
    chunkedids = []
    last = 0.0  
    while last < len(listofids):
      chunkedids.append(listofids[int(last):int(last + chunksize)])
      last += chunksize
    
    #Go through each entry in chunkeds and send its elements out to worker
    #processes, each process take a single raw data file, finds its features
    #and writes them
    prog=0
    for i in chunkedids:
        jobs = []
        for j in i:
            p = multiprocessing.Process(target=write_all_features, args=(j,))
            jobs.append(p)
            p.start()
        for p in jobs:
            p.join()
        prog = prog+1
        print prog