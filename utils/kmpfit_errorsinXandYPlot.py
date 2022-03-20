#!/usr/bin/env python
#------------------------------------------------------------
# Purpose: Program to straight line parameters
#          to data with errors in both coordinates. Compare 
#          the results with SciPy's ODR routine.
# Vog, 27 Nov, 2011
#------------------------------------------------------------
import numpy
from matplotlib.pyplot import figure, show, rc
from numpy.random import normal
from kmpfit import kmpfit
from scipy.odr import Data, Model, ODR, RealData, odr_stop


def kmp_linear_fitting(x,errx,y,erry):

	def model(p, x):
	   # Model is staight line: y = a + b*x
	   a, b = p
	   return a + b*x
	
	def residuals(p, data):
	   # Residuals function for effective variance
	   a, b = p
	   x, y, ex, ey = data
	   w = ey*ey + b*b*ex*ex
	   wi = numpy.sqrt(numpy.where(w==0.0, 0.0, 1.0/(w)))
	   d = wi*(y-model(p,x))
	   return d
	
	def residuals2(p, data):
	   # Minimum distance formula with expression for x_model
	   a, b = p
	   x, y, ex, ey = data
	   wx = 1/(ex*ex)
	   wy = 1/(ey*ey)
	   df = b
	   xd = x + (wy*(y-model(p,x))*df)/(wx+wy*df*df)
	   yd = model(p,xd)
	   D = numpy.sqrt( wx*(x-xd)**2+wy*(y-yd)**2 )
	   return D
	
	
	beta0 = [0,0]
	print "\n========== Results SciPy's ODR ============"
	linear = Model(model)
	mydata = RealData(x, y, sx=errx, sy=erry)
	myodr = ODR(mydata, linear, beta0=beta0, maxit=5000)
	myoutput = myodr.run()
	print "Fitted parameters:      ", myoutput.beta
	print "Covariance errors:      ", numpy.sqrt(myoutput.cov_beta.diagonal())
	print "Standard errors:        ", myoutput.sd_beta
	print "Minimum (reduced)chi^2: ", myoutput.res_var
	beta = myoutput.beta
	
	# Prepare fit routine
	fitobj = kmpfit.Fitter(residuals=residuals, data=(x, y, errx, erry))
	try:
	   fitobj.fit(params0=beta0)
	except Exception, mes:
	   print "Something wrong with fit: ", mes
	   raise SystemExit
	
	print "\n\n======== Results kmpfit: w1 = ey*ey + b*b*ex*ex ========="
	print "Params:                 ", fitobj.params
	print "Covariance errors:      ", fitobj.xerror
	print "Standard errors         ", fitobj.stderr
	print "Chi^2 min:              ", fitobj.chi2_min
	print "Reduced Chi^2:          ", fitobj.rchi2_min
	print "Message:                ", fitobj.message
	
	fitobj2 = kmpfit.Fitter(residuals=residuals2, data=(x, y, errx, erry))
	try:
	   fitobj2.fit(params0=beta0)
	except Exception, mes:
	   print "Something wrong with fit: ", mes
	   raise SystemExit
	
	print "\n\n======== Results kmpfit: r = ex*ex/(ey*ey), xd = (x-a*r+y*b*r)/(1+r) ========="
	print "Params:                 ", fitobj2.params
	print "Covariance errors:      ", fitobj2.xerror
	print "Standard errors         ", fitobj2.stderr
	print "Chi^2 min:              ", fitobj2.chi2_min
	print "Reduced Chi^2:          ", fitobj2.rchi2_min
	print "Message:                ", fitobj2.message
	
	
	t = "\nTHE WILLAMSON APPROACH"
	print t, "\n", "="*len(t)
	# Step 1: Get a and b for a, b with standard weighted least squares calculation
	
	def lingres(xa, ya, w):
	   # Return a, b for the relation y = a + b*x
	   # given data in xa, ya and weights in w
	   sum   =  w.sum()
	   sumX  = (w*xa).sum()
	   sumY  = (w*ya).sum()
	   sumX2 = (w*xa*xa).sum()
	   sumY2 = (w*ya*ya).sum()
	   sumXY = (w*xa*ya).sum()
	   delta = sum * sumX2 - sumX * sumX 
	   a = (sumX2*sumY - sumX*sumXY) / delta
	   b = (sumXY*sum - sumX*sumY) / delta
	   return a, b
	
	w = numpy.where(erry==0.0, 0.0, 1.0/(erry*erry))
	a,b = lingres(x, y, w)
	a_y = a; b_y = b       # Williamson initial Parameters
	
	ui = errx**2
	vi = erry**2
	
	n = 0 
	cont = True
	while cont:
	   # Step 2: Use this slope to find weighting for each point
	   wi = (vi+b*b*ui)**-1
	
	   # Step 3: Calcu;ate weighted avarages
	   w_sum = wi.sum()
	   x_av = (wi*x).sum() / w_sum
	   x_diff = x - x_av
	   y_av = (wi*y).sum() / w_sum
	   y_diff = y - y_av
	
	   # Step 4: Calculate the 'improvement' vector zi
	   zi = wi*(vi*x_diff + b*ui*y_diff)
	   b_will = (wi*zi*y_diff).sum()/ (wi*zi*x_diff).sum()
	   cont = abs(b-b_will) > 1e-12 and n < 100
	   n += 1
	   b = b_will
	
	# Step 5: Repeat steps 2-4 until convergence
	
	# Step 6: Calculate 'a' using the averages of a and y
	a_will = y_av - b_will*x_av    # Improved parameters
	
	# Step 7: The variances
	wi = (vi+b_will*b_will*ui)**-1
	w_sum = wi.sum()
	
	z_av = (wi*zi).sum() / w_sum
	zi2 = zi - z_av
	Q =1.0/(wi*(x_diff*y_diff/b_will + 4*zi2*(zi-x_diff))).sum()
	sigb2 = Q*Q * (wi*wi*(x_diff**2*vi+y_diff**2*ui)).sum()
	siga2 = 1.0/w_sum + 2*(x_av+2*z_av)*z_av*Q + (x_av+2*z_av)**2*sigb2
	siga = numpy.sqrt(siga2)
	sigb = numpy.sqrt(sigb2)
	
	print "Williamson Fitted A, B: ", a_will, b_will
	print "Parameter errors: ", siga, sigb
	
	# Some plotting
	rc('font', size=9)
	rc('legend', fontsize=8)
	fig = figure(1)
	frame = fig.add_subplot(1,1,1, aspect=1, adjustable='datalim')
	frame.errorbar(x, y, xerr=errx, yerr=erry,  fmt='bo')
	# Plot first fit
	frame.plot(x, model(beta,x), '-y', lw=4, label="ODR", alpha=0.6)
	frame.plot(x, model(fitobj.params,x), 'c', ls='--', lw=2, label="kmpfit")
	frame.plot(x, model(fitobj2.params,x), '#ffaa00', label="kmpfit correct")
	frame.plot(x, model((a_will,b_will),x), 'g', label="Williamson")
	frame.set_xlabel("X")
	frame.set_ylabel("Y")
	frame.set_title("Weights in both coordinates. Model: $y=a+bx$")
	frame.grid(True)
	leg = frame.legend(loc=1)
	show()

	return fitobj2.params
