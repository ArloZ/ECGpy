#! /bin/python 
# read and display the MIT-BIH ecg data

__author__ = "qile68@163.com"

import os,sys,struct
import matplotlib.pyplot as plt
import numpy as np



FILE_DIR = "../MIT-BIH"
#FILE_DIR = "E:\Workspace\labplatform\ecg\MIT-BIH"

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
	_filename = ""
	_fs = 0			# sample rate
	_num = 0			# sample number
	_format = "212"  # data format saved
	_chn = {}		# ecg lead name
	_gain = {}		# channel gain
	_data = {}		# ecg data array
	_atr = []		# data annotations 
	def __init__(self,name = ""):
		self._filename = name
		self.read_hea()
	
	def reset_property(self):
		self._filename = ""
		self._fs = 0
		self._num = 0
		self._format = 0
		self._chn = {}
		self._gain = {}
		self._atr = []
	
	def get_property(self):
		return [self._filename,self._fs,self._num,self._format,self._chn,self._gain]
	
	def read_hea(self):
		if self._filename == "":
			return
		dir = FILE_DIR + os.sep + self._filename + ".hea"
		fp = open(dir,"r")
		line = fp.readline()
		vals = line.split(" ")
		if len(vals) < 4:
			return
		n = vals[1]
		self._fs = int(vals[2])
		self._num = int(vals[3])
		for i in range(int(n)):
			line = fp.readline()
			vals = line.split(" ")
			self._chn[i] = vals[8]
			self._gain[i] = vals[2]
			self._format = vals[1]
			self._data[i] = []
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
			
		
	def read_atr(self):
		if self._filename == "":
			return
		dir = FILE_DIR + os.sep + self._filename + ".atr"
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
		self._atr = [point[1:],annotation]
	
	def read_dat(self):
		if self._filename == "" or self._format != "212":
			return
		dir = FILE_DIR + os.sep + self._filename + ".dat"
		fp = open(dir,"rb") 							# must add "b" under windows for read as binary
		for i in range(self._num):		
			v = struct.unpack("BBB",fp.read(3))			# read 3 bytes(as unsigned char) for two channel value once a time
			ch1 = ((v[1] & 0x0f) << 8) + v[0]			# example: v : E3 33 F3
			ch2 = ((v[1] & 0xf0) << 4) + v[2]			# 			ch1: 3E3	ch2: 3F3	
			self._data[0].append(ch1)
			self._data[1].append(ch2)
		fp.close()
			
	def set_file(self,name = ""):
		if name == "":
			return False
		self.reset_property()
		self._filename = name
		self.read_hea()
	
	def get_total_peak(self):
		"""
		get the total peaks of this ecg data
		"""
		annot = self._atr[1]
		annd = {}
		for ann in annot:
			if ann in annd:
				annd[ann] += 1
			else:
				annd[ann] = 1
		total = 0
		other = ["?","~","s","p","t","u","(",")","[","]"]
		for (k,v) in annd.items():
			if k in other:
				continue
			total += v
		return total
		
	def get_atr(self,start = 0,seconds = 10):
		"""
		get annotations of the ecg data
		input : start means the annotations data where to start
				seconds means the length of the signal in seconds
						if seconds is less than zero will return all the data
		"""
		if seconds <= 0:
			return self._atr
		start = start * self._fs
		n = seconds * self._fs + start
		if start > max(self._atr[0]):
			return []
		point,annot = self._atr
		begin = 0
		end = -1
		for i in range(len(point)):
			if point[i] >= start:
				begin = i
				for j in range(len(point)-i):
					if point[i + j] > n:
						end = i+j
						break
				break
		if end > 0:
			return [point[begin:end],annot[begin:end]]
		else:
			return [point[begin:],annot[begin:]]

	def get_dat(self,ch = 0,start = 0,seconds = 10):
		"""
		get ecgdata which signal data is a list
		input : start means the data where to start
				seconds means the length of the signal in seconds
						if seconds is less than zero will return all the data
		"""
		if ch >= len(self._data):
			return []
		if seconds <= 0:  # get all data
			return self._data[ch]
		n = seconds * self._fs
		start = start * self._fs
		num = len(self._data[ch])
		if start > num:
			return []
		if start + n > num:
			n = num - start
		return self._data[ch][start:start + n]
	
	def get_fs(self):
		return self._fs
		
	# plot the ecg during "dur" seconds from the "start" time
	def plot(self,start = 0,dur = 10):	
		n = self._fs * dur
		if self._fs * start + n > self._num:
			n = self._num - self._fs*start
			dur = n / self._fs
		x = np.linspace(start,start + dur,n)
		# plot channel 1
		plt.subplot(211)
		plt.plot(x,self._data[0][:n],"r",label=self._chn[0])
		plt.xlabel("Time(s)")
		plt.ylabel("Volt(mv)")
		plt.legend()
		# plot channel 2
		plt.subplot(212)
		plt.plot(x,self._data[1][:n],"r",label=self._chn[1])
		plt.xlabel("Time(s)")
		plt.ylabel("Volt(mv)")
		plt.legend()
		plt.show()
		
if __name__ == "__main__":
	# test EcgData
	ecgData = EcgData("100")
	ecgData.read_dat()
	ecgData.plot(0,10)
