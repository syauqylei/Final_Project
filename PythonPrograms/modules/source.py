"""
this program contains
Source function that is used in the research

"""
from math import sin ,exp,pi

def source_1(f,t):#zahradnik et al 1993
	return sin(2*pi*f*t)*exp(-4*(pi*f*t)**2/16)


def source_2(f,t):#Yang et al 2012
	return 5.76*f*f*(1.0-16.0*(0.6*f*t-1.0)**2)*exp(-8.0*(0.6*f*t-1)**2)
