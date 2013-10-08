#! /bin/python 
# read and display the MIT-BIH ecg data

__author__ = "qile68@163.com"

import os,sys,struct
import matplotlib.pyplot as plt
import numpy as np



FILE_DIR = "./MIT-BIH"

# annotation about ecg see 
# http://www.physionet.org/physiobank/annotations.shtml
# http://www.physionet.org/physiotools/wfdb/lib/ecgcodes.h

ANNOT_LIST = ["NOT",	# NOTQRS	0	/* not-QRS (not a getann/putann code) */
			  "N",		# NORMAL	1	/* normal beat */
			  "L",		# LBBB	2	/* left bundle branch block beat */
			  "R",		# RBBB	3	/* right bundle branch block beat */
			  "a",		# ABERR	4	/* aberrated atrial premature beat */
			  "V",		# PVC	5	/* premature ventricular contraction */
			  "F",		# FUSION	6	/* fusion of ventricular and normal beat */
			  "J",		# NPC	7	/* nodal (junctional) premature beat *
			  "A",		# APC	8	/* atrial premature contraction */
			  "S",		# SVPB	9	/* premature or ectopic supraventricular beat */
			  "E",		# VESC	10	/* ventricular escape beat */
			  "j",		# NESC	11	/* nodal (junctional) escape beat */
			  "/",		# PACE	12	/* paced beat */
			  "Q",		# UNKNOWN	13	/* unclassifiable beat */
			  "~",		# NOISE	14	/* signal quality change */
			  "|",		# ARFCT	16	/* isolated QRS-like artifact */
			  "s",		# STCH	18	/* ST change */
			  "T",		# TCH	19	/* T-wave change */
			  "*",		# SYSTOLE	20	/* systole */
			  "D",		# DIASTOLE 21	/* diastole */
			  "\"",		# NOTE	22	/* comment annotation */
			  "=",		# MEASURE 23	/* measurement annotation */
			  "p",		# PWAVE	24	/* P-wave peak */
			  "B",		# BBB	25	/* left or right bundle branch block */
			  "^",		# PACESP	26	/* non-conducted pacer spike */
			  "t",		# TWAVE	27	/* T-wave peak */
			  "+",		# RHYTHM	28	/* rhythm change */
			  "u",		# UWAVE	29	/* U-wave peak */
			  "?",		# LEARN	30	/* learning */
			  "!",		# FLWAV	31	/* ventricular flutter wave */
			  "[",		# VFON	32	/* start of ventricular flutter/fibrillation */
			  "]",		# VFOFF	33	/* end of ventricular flutter/fibrillation */
			  "e",		# AESC	34	/* atrial escape beat */
			  "n",		# SVESC	35	/* supraventricular escape beat */
			  "@",		# LINK    36	/* link to external data (aux contains URL) */
			  "x",		# NAPC	37	/* non-conducted P-wave (blocked APB) */
			  "f",		# PFUS	38	/* fusion of paced and normal beat */
			  "(",		# WFON	39	/* waveform onset */
			  "`",		# PQ	WFON	/* PQ junction (beginning of QRS) */
			  ")",		# WFOFF	40	/* waveform end */
			  "'",		# JPT	WFOFF	/* J point (end of QRS) */
			  "r",		# RONT	41	/* R-on-T premature ventricular contraction */
			  "MAX",		# ACMAX	49	/* value of largest valid annot code (must be < 50) */			  
			  ""]
class EcgData():
	filename = ""
	fs = 0			# sample rate
	num = 0			# sample number
	format = "212"  # data format saved
	chn = {}		# ecg lead name
	gain = {}		# channel gain
	data = {}		# ecg data array
	def __init__(self,name = ""):
		self.filename = name
		self.readHea()
	
	def resetProperty(self):
		self.filename = ""
		self.fs = 0
		self.num = 0
		self.format = 0
		self.chn = {}
		self.gain = {}
	
	def getProperty(self):
		return [self.filename,self.fs,self.num,self.format,self.chn,self.gain]
	
	def readHea(self):
		if self.filename == "":
			return
		dir = FILE_DIR + os.sep + self.filename + ".hea"
		fp = open(dir,"r")
		line = fp.readline()
		vals = line.split(" ")
		if len(vals) < 4:
			return
		n = vals[1]
		self.fs = int(vals[2])
		self.num = int(vals[3])
		for i in range(int(n)):
			line = fp.readline()
			vals = line.split(" ")
			self.chn[i] = vals[8]
			self.gain[i] = vals[2]
			self.format = vals[1]
			self.data[i] = []
		fp.close()

	def __opcode(self,code,time):
		if code < 50:
			return True,ANNOT_LIST[code]
		elif code == 59:		
			return False,4
		elif 59 < code < 63:
			return False,0
		elif code == 63:
			if time % 2 == 1:
				return False,time + 1
			else:
				return False,time
		else:
			return False,0
			
		
	def readAtr(self):
		if self.filename == "":
			return
		dir = FILE_DIR + os.sep + self.filename + ".atr"
		fp = open(dir,"rb")
		point = [0]
		annotation = []
		while True:
			v = struct.unpack("BB",fp.read(2))
			if v[0] == 0 and v[1] == 0:
				break
			code = v[1] >> 2
			time = ((v[1] & 0x03) << 8) + v[0]
			ret,val = self.__opcode(code,time)
			if ret:
				point.append(point[-1] + time)
				annotation.append(val)
			else:
				x = fp.read(val)
		return point[1:],annotation
	
	def readDat(self):
		if self.filename == "" or self.format != "212":
			return
		dir = FILE_DIR + os.sep + self.filename + ".dat"
		fp = open(dir,"rb") 							# must add "b" under windows for read as binary
		for i in range(self.num):		
			v = struct.unpack("BBB",fp.read(3))			# read 3 bytes(as unsigned char) for two channel value once a time
			ch1 = ((v[1] & 0x0f) << 8) + v[0]			# example: v : E3 33 F3
			ch2 = ((v[1] & 0xf0) << 4) + v[2]			# 			ch1: 3E3	ch2: 3F3	
			self.data[0].append(ch1)
			self.data[1].append(ch2)
		fp.close()
			
	def setFile(self,name = ""):
		if name == "":
			return False
		self.resetProperty()
		self.filename = name
		self.readHea()

	def getDat(self,ch = 0,start = 0,n = 500):
		if ch >= len(self.data):
			return []
		if n <= 0:  # get all data
			return self.data[ch]
		num = len(self.data[ch])
		if start > num:
			return []
		if start + n > num:
			n = num - start
		return self.data[ch][start:start + n]
	
	def getFs(self):
		return self.fs
		
	# plot the ecg during "dur" seconds from the "start" time
	def plot(self,start = 0,dur = 10):	
		n = self.fs * dur
		if self.fs * start + n > self.num:
			n = self.num - self.fs*start
			dur = n / self.fs
		x = np.linspace(start,start + dur,n)
		# plot channel 1
		plt.subplot(211)
		plt.plot(x,self.data[0][:n],"r",label=self.chn[0])
		plt.xlabel("Time(s)")
		plt.ylabel("Volt(mv)")
		plt.legend()
		# plot channel 2
		plt.subplot(212)
		plt.plot(x,self.data[1][:n],"r",label=self.chn[1])
		plt.xlabel("Time(s)")
		plt.ylabel("Volt(mv)")
		plt.legend()
		plt.show()
		
if __name__ == "__main__":
	# test EcgData
	ecgData = EcgData("100")
	ecgData.readDat()
	ecgData.plot(0,10)
