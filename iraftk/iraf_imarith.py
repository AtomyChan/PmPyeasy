#! /usr/bin/python

from pyraf import iraf 
from iraf import imarith

def imarith_iraf(operand1, operator, operand2, result):
	'''
	INPUTS:
		operand1: operand image or list of images or constant
		operator: oerator to be applied to the operands, allowed operators are '+','-','*','/','min','max'

		result: resultant image
	'''


	imarith.unlearn()
	imarith(operand1=operand1, op=operator, operand2=operand2, result=result)


if __name__ == "__main__":
	import sys

	

	opd1 = sys.argv[1]

	if opd1 == '--help':
		print "usage: thisfile operand1 operator operand2 output"
		sys.exit()

	op = sys.argv[2]
	opd2 = sys.argv[3]

	result = sys.argv[4]

	imarith_iraf(opd1, op, opd2, result)

