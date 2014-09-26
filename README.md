automatedAFManalysis
====================

Automated analysis of single bond cell/afm adhesion experiments

This program automatically analyzes data taken in AFM single cell force spectroscopy experiments. A screen shot of the typical
output is included.

Overview:

The program uses a support vector classifier to automatically go through single experiment force vs time series and pick out those
experiments (cycles) that possess a single bond.  It classifies these experiments as either tethered (the adhesion force between 
the cell and the cantilever pulls out a nano tube from the membrane) or untethered (the external bond between the surface protein 
and the cantilever breaks before a tube can be pulled). Once a cycle is classified it is automatically analyzed to extract a slope
(loading rate in pN/s) and a force.  In the case of tethered experiments the slope it reports is the greatest slope and the force
is the "crossover" force, the force where the cytoskeletal linkage is broken.  This is the adhesive force at the point where the
nanotube begins to be pulled out.  In the case of untethered experiments the slope is that right before the rupture of the bond
and the force is the bond rupture force.

In addition to analyzing the data the program will also compare its output (both of the classifier and the force/slope data) with
that of any "ground truth" (hand analyzed) data that maybe present in the SQL database with which it is associated.

Using the Program:

In order to run the program it will require a connection to an SQL database containing file locations for force time results.
You will have to go to afm_sql_code\afmsqlcode.py and enter your database information where it says #ENTER DATABASE INFO HERE.
A picture of the interface can be seen in interface.png

The left side of the interface contains lists with several experimental parameters.  These are loaded from the database.
Selecting an option or set of options from the lists will remove options from other lists if they are no longer valid.  For
example if you took data with cell type 1 on the 3rd of April and data with cell type 2 on the 4th, then selecting the 3rd from
the date column will remove cell type 2 from the pipette object list.

Once a set of experimental parameters are picked the plot button will run the analysis.  The results and any ground truth data
will then be displayed on the graphs on the right.  The results can then be saved via the menu option.

Classification:

The script feature_creation.py takes force vs. time data from individual experiments and extracts features from them that are
used for machine learning.  The features used include:

        mean force, max force, total number of discontinuities (bond ruptures), largest slope, largest drop (bond rupture),
        number of peaks in the power spectrum, the force time data (binned into 16 bins), power spectrum (binned into 16 bins)

The script also loads in the labels for the data from the database as well as the file id for the experiment and writes them
to a csv file.

The script bond_learning.py then feads the labeled data into scikitlearn SVC object and uses cross validation to create a
classifier.

Data Structure:

Data is assumed to be taken from the experimental apparatus and saved in large initial raw data files containing dozens or
hundreds of experiments all with the same parameters.  After this it is also assumed that they have been broken up into smaller,
calibrated force time csv's.  Each experiment is assumed to have 4 associated files, one contains data of the approach to the
cantilever, one contains data for the compression of the cell against the cantilever, on contains a holding or 'touch' phase,
and finally one contains the retraction data.

  SQL Database:
  
  In order to make the program work you will have to go to afm_sql_code\afmsqlcode.py and enter your database information where
  it says #ENTER DATABASE INFO HERE.
  
  The database is assumed to have the following tables:
  
  filelist: A table listing raw data files each with hundreds of experiments.  It contains data like the buffer the experiments
            were done in, the object being tested, the loading rate, the date, etc.
  
  dat: A table containing the locations of the individual calibrated force time files.
  
  cell_data: A table containing the extracted force, time etc.
  
  The exact columns and database structure the program expects are laid out in afmsqlcode.py
