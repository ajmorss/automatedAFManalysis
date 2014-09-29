# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'connectdialog.ui'
#
# Created: Fri Sep 26 13:39:41 2014
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

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName(_fromUtf8("Dialog"))
        Dialog.resize(292, 183)
        self.label = QtGui.QLabel(Dialog)
        self.label.setGeometry(QtCore.QRect(10, 10, 61, 16))
        self.label.setObjectName(_fromUtf8("label"))
        self.label_2 = QtGui.QLabel(Dialog)
        self.label_2.setGeometry(QtCore.QRect(150, 10, 46, 13))
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.label_3 = QtGui.QLabel(Dialog)
        self.label_3.setGeometry(QtCore.QRect(10, 80, 61, 16))
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.label_4 = QtGui.QLabel(Dialog)
        self.label_4.setGeometry(QtCore.QRect(150, 80, 46, 13))
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.hostnameLE = QtGui.QLineEdit(Dialog)
        self.hostnameLE.setGeometry(QtCore.QRect(10, 30, 113, 20))
        self.hostnameLE.setObjectName(_fromUtf8("hostnameLE"))
        self.portLE = QtGui.QLineEdit(Dialog)
        self.portLE.setGeometry(QtCore.QRect(150, 30, 113, 20))
        self.portLE.setObjectName(_fromUtf8("portLE"))
        self.usernameLE = QtGui.QLineEdit(Dialog)
        self.usernameLE.setGeometry(QtCore.QRect(10, 100, 113, 20))
        self.usernameLE.setObjectName(_fromUtf8("usernameLE"))
        self.passwordLE = QtGui.QLineEdit(Dialog)
        self.passwordLE.setGeometry(QtCore.QRect(150, 100, 113, 20))
        self.passwordLE.setEchoMode(QtGui.QLineEdit.Password)
        self.passwordLE.setObjectName(_fromUtf8("passwordLE"))
        self.okB = QtGui.QPushButton(Dialog)
        self.okB.setGeometry(QtCore.QRect(190, 150, 75, 23))
        self.okB.setObjectName(_fromUtf8("okB"))

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(_translate("Dialog", "Dialog", None))
        self.label.setText(_translate("Dialog", "Host Name", None))
        self.label_2.setText(_translate("Dialog", "Port", None))
        self.label_3.setText(_translate("Dialog", "Username", None))
        self.label_4.setText(_translate("Dialog", "Password", None))
        self.okB.setText(_translate("Dialog", "PushButton", None))

