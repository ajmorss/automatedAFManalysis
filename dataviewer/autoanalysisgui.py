# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'autmagicanalysisgui.ui'
#
# Created: Fri Aug 29 10:17:25 2014
#      by: PyQt4 UI code generator 4.9.6
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(1203, 931)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.groupBox = QtGui.QGroupBox(self.centralwidget)
        self.groupBox.setGeometry(QtCore.QRect(130, 10, 181, 80))
        self.groupBox.setObjectName(_fromUtf8("groupBox"))
        self.untethered = QtGui.QRadioButton(self.groupBox)
        self.untethered.setGeometry(QtCore.QRect(20, 50, 151, 17))
        self.untethered.setObjectName(_fromUtf8("untethered"))
        self.tethered = QtGui.QRadioButton(self.groupBox)
        self.tethered.setGeometry(QtCore.QRect(20, 20, 131, 17))
        self.tethered.setChecked(True)
        self.tethered.setObjectName(_fromUtf8("tethered"))
        self.groupBox_2 = QtGui.QGroupBox(self.centralwidget)
        self.groupBox_2.setGeometry(QtCore.QRect(10, 100, 631, 621))
        self.groupBox_2.setObjectName(_fromUtf8("groupBox_2"))
        self.label = QtGui.QLabel(self.groupBox_2)
        self.label.setGeometry(QtCore.QRect(20, 20, 81, 16))
        self.label.setObjectName(_fromUtf8("label"))
        self.label_2 = QtGui.QLabel(self.groupBox_2)
        self.label_2.setGeometry(QtCore.QRect(20, 330, 91, 16))
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.label_3 = QtGui.QLabel(self.groupBox_2)
        self.label_3.setGeometry(QtCore.QRect(330, 20, 91, 16))
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.label_4 = QtGui.QLabel(self.groupBox_2)
        self.label_4.setGeometry(QtCore.QRect(330, 330, 81, 16))
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.label_5 = QtGui.QLabel(self.groupBox_2)
        self.label_5.setGeometry(QtCore.QRect(450, 330, 111, 16))
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.pipette_object = QtGui.QTableView(self.groupBox_2)
        self.pipette_object.setGeometry(QtCore.QRect(20, 40, 271, 271))
        self.pipette_object.setObjectName(_fromUtf8("pipette_object"))
        self.cantilever_object = QtGui.QTableView(self.groupBox_2)
        self.cantilever_object.setGeometry(QtCore.QRect(330, 40, 271, 271))
        self.cantilever_object.setObjectName(_fromUtf8("cantilever_object"))
        self.buffer = QtGui.QTableView(self.groupBox_2)
        self.buffer.setGeometry(QtCore.QRect(20, 350, 271, 251))
        self.buffer.setObjectName(_fromUtf8("buffer"))
        self.date = QtGui.QTableView(self.groupBox_2)
        self.date.setGeometry(QtCore.QRect(330, 350, 91, 251))
        self.date.setObjectName(_fromUtf8("date"))
        self.loading_rate = QtGui.QTableView(self.groupBox_2)
        self.loading_rate.setGeometry(QtCore.QRect(450, 350, 111, 251))
        self.loading_rate.setObjectName(_fromUtf8("loading_rate"))
        self.plot_button = QtGui.QPushButton(self.centralwidget)
        self.plot_button.setGeometry(QtCore.QRect(30, 790, 75, 23))
        self.plot_button.setObjectName(_fromUtf8("plot_button"))
        self.toolframe = QtGui.QFrame(self.centralwidget)
        self.toolframe.setGeometry(QtCore.QRect(830, 830, 311, 41))
        self.toolframe.setFrameShape(QtGui.QFrame.StyledPanel)
        self.toolframe.setFrameShadow(QtGui.QFrame.Raised)
        self.toolframe.setObjectName(_fromUtf8("toolframe"))
        self.label_7 = QtGui.QLabel(self.centralwidget)
        self.label_7.setGeometry(QtCore.QRect(30, 760, 71, 16))
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.total_cycles = QtGui.QLabel(self.centralwidget)
        self.total_cycles.setGeometry(QtCore.QRect(120, 760, 46, 16))
        self.total_cycles.setObjectName(_fromUtf8("total_cycles"))
        self.groupBox_3 = QtGui.QGroupBox(self.centralwidget)
        self.groupBox_3.setGeometry(QtCore.QRect(30, 10, 91, 80))
        self.groupBox_3.setObjectName(_fromUtf8("groupBox_3"))
        self.celldat = QtGui.QRadioButton(self.groupBox_3)
        self.celldat.setGeometry(QtCore.QRect(10, 20, 82, 17))
        self.celldat.setChecked(True)
        self.celldat.setObjectName(_fromUtf8("celldat"))
        self.beaddat = QtGui.QRadioButton(self.groupBox_3)
        self.beaddat.setEnabled(False)
        self.beaddat.setGeometry(QtCore.QRect(10, 50, 82, 17))
        self.beaddat.setObjectName(_fromUtf8("beaddat"))
        self.tabWidget = QtGui.QTabWidget(self.centralwidget)
        self.tabWidget.setGeometry(QtCore.QRect(650, 10, 541, 861))
        self.tabWidget.setAutoFillBackground(True)
        self.tabWidget.setStyleSheet(_fromUtf8("\n"
"\n"
"    QTabWidget>QWidget>QWidget{background:rgba(95%,95%, 95%, 100%);}\n"
""))
        self.tabWidget.setObjectName(_fromUtf8("tabWidget"))
        self.tscatter = QtGui.QWidget()
        self.tscatter.setObjectName(_fromUtf8("tscatter"))
        self.tscattframe = QtGui.QFrame(self.tscatter)
        self.tscattframe.setGeometry(QtCore.QRect(10, 10, 511, 811))
        self.tscattframe.setFrameShape(QtGui.QFrame.StyledPanel)
        self.tscattframe.setFrameShadow(QtGui.QFrame.Raised)
        self.tscattframe.setObjectName(_fromUtf8("tscattframe"))
        self.tabWidget.addTab(self.tscatter, _fromUtf8(""))
        self.thist = QtGui.QWidget()
        self.thist.setObjectName(_fromUtf8("thist"))
        self.thistframe = QtGui.QFrame(self.thist)
        self.thistframe.setGeometry(QtCore.QRect(10, 10, 511, 811))
        self.thistframe.setFrameShape(QtGui.QFrame.StyledPanel)
        self.thistframe.setFrameShadow(QtGui.QFrame.Raised)
        self.thistframe.setObjectName(_fromUtf8("thistframe"))
        self.tabWidget.addTab(self.thist, _fromUtf8(""))
        self.utscatter = QtGui.QWidget()
        self.utscatter.setObjectName(_fromUtf8("utscatter"))
        self.utscattframe = QtGui.QFrame(self.utscatter)
        self.utscattframe.setGeometry(QtCore.QRect(10, 10, 511, 811))
        self.utscattframe.setFrameShape(QtGui.QFrame.StyledPanel)
        self.utscattframe.setFrameShadow(QtGui.QFrame.Raised)
        self.utscattframe.setObjectName(_fromUtf8("utscattframe"))
        self.tabWidget.addTab(self.utscatter, _fromUtf8(""))
        self.uthist = QtGui.QWidget()
        self.uthist.setObjectName(_fromUtf8("uthist"))
        self.uthistframe = QtGui.QFrame(self.uthist)
        self.uthistframe.setGeometry(QtCore.QRect(10, 10, 511, 811))
        self.uthistframe.setFrameShape(QtGui.QFrame.StyledPanel)
        self.uthistframe.setFrameShadow(QtGui.QFrame.Raised)
        self.uthistframe.setObjectName(_fromUtf8("uthistframe"))
        self.tabWidget.addTab(self.uthist, _fromUtf8(""))
        self.confusionmattab = QtGui.QWidget()
        self.confusionmattab.setObjectName(_fromUtf8("confusionmattab"))
        self.confusionmatrix = QtGui.QTableView(self.confusionmattab)
        self.confusionmatrix.setGeometry(QtCore.QRect(140, 300, 256, 192))
        self.confusionmatrix.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.confusionmatrix.setObjectName(_fromUtf8("confusionmatrix"))
        self.tabWidget.addTab(self.confusionmattab, _fromUtf8(""))
        self.attachments = QtGui.QLabel(self.centralwidget)
        self.attachments.setGeometry(QtCore.QRect(120, 740, 46, 13))
        self.attachments.setObjectName(_fromUtf8("attachments"))
        self.holdbox = QtGui.QCheckBox(self.centralwidget)
        self.holdbox.setGeometry(QtCore.QRect(330, 750, 70, 17))
        self.holdbox.setObjectName(_fromUtf8("holdbox"))
        self.progressBar = QtGui.QProgressBar(self.centralwidget)
        self.progressBar.setGeometry(QtCore.QRect(30, 820, 261, 23))
        self.progressBar.setProperty("value", 100)
        self.progressBar.setObjectName(_fromUtf8("progressBar"))
        self.cancelbutton = QtGui.QPushButton(self.centralwidget)
        self.cancelbutton.setEnabled(False)
        self.cancelbutton.setGeometry(QtCore.QRect(120, 790, 75, 23))
        self.cancelbutton.setObjectName(_fromUtf8("cancelbutton"))
        self.proglabel = QtGui.QLabel(self.centralwidget)
        self.proglabel.setGeometry(QtCore.QRect(30, 850, 251, 21))
        self.proglabel.setText(_fromUtf8(""))
        self.proglabel.setAlignment(QtCore.Qt.AlignCenter)
        self.proglabel.setObjectName(_fromUtf8("proglabel"))
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1203, 21))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        self.menuFile = QtGui.QMenu(self.menubar)
        self.menuFile.setObjectName(_fromUtf8("menuFile"))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)
        self.actionSave = QtGui.QAction(MainWindow)
        self.actionSave.setObjectName(_fromUtf8("actionSave"))
        self.actionExit = QtGui.QAction(MainWindow)
        self.actionExit.setObjectName(_fromUtf8("actionExit"))
        self.menuFile.addAction(self.actionSave)
        self.menuFile.addAction(self.actionExit)
        self.menubar.addAction(self.menuFile.menuAction())

        self.retranslateUi(MainWindow)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow", None))
        self.groupBox.setTitle(_translate("MainWindow", "Give me:", None))
        self.untethered.setText(_translate("MainWindow", "Untethered Attachments", None))
        self.tethered.setText(_translate("MainWindow", "Tethered Attachments", None))
        self.groupBox_2.setTitle(_translate("MainWindow", "With Parameters:", None))
        self.label.setText(_translate("MainWindow", "Pipette Object", None))
        self.label_2.setText(_translate("MainWindow", "Buffer", None))
        self.label_3.setText(_translate("MainWindow", "Cantilever Object", None))
        self.label_4.setText(_translate("MainWindow", "Dates", None))
        self.label_5.setText(_translate("MainWindow", "Nominal Loading Rate", None))
        self.plot_button.setText(_translate("MainWindow", "Plot!", None))
        self.label_7.setText(_translate("MainWindow", "Total Cycles:", None))
        self.total_cycles.setText(_translate("MainWindow", "0", None))
        self.groupBox_3.setTitle(_translate("MainWindow", "Give me:", None))
        self.celldat.setText(_translate("MainWindow", "Cell Data", None))
        self.beaddat.setText(_translate("MainWindow", "Bead Data", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tscatter), _translate("MainWindow", "Tether Scatter Plots", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.thist), _translate("MainWindow", "Tether Histograms", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.utscatter), _translate("MainWindow", "Untethered Scatter Plots", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.uthist), _translate("MainWindow", "Untethered Histograms", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.confusionmattab), _translate("MainWindow", "Confusion Matrix", None))
        self.attachments.setText(_translate("MainWindow", "TextLabel", None))
        self.holdbox.setText(_translate("MainWindow", "CheckBox", None))
        self.cancelbutton.setText(_translate("MainWindow", "Cancel", None))
        self.menuFile.setTitle(_translate("MainWindow", "File", None))
        self.actionSave.setText(_translate("MainWindow", "Save", None))
        self.actionExit.setText(_translate("MainWindow", "Exit", None))

