SimSpec.iraf.task(bias4k=SimSpec.os.path.dirname(SimSpec.param['precleanscript'])+'/bias4k.cl')
SimSpec.iraf.bias4k(SimSpec.WorkDir+'*_SQ???.fits', suffix='_bcor')
SimSpec.iraf.imdelete(SimSpec.WorkDir+'*_bcor?.fits')
for f in SimSpec.glob.glob(SimSpec.WorkDir+'*_bcor.fits'):
    SimSpec.os.system("mv -f '"+f+"' '"+f.replace("_bcor.fits",".fits")+"'")


SimSpec.iraf.hedit.unlearn()
#SimSpec.iraf.hedit(SimSpec.WorkDir+'*_SQ???.fits', 'EPOCH', "(EQUINOX)", add='yes', addonly='yes', verify='no')
SimSpec.iraf.hedit(SimSpec.WorkDir+'*_SQ???.fits', 'ccdsec,datasec', delete='yes', verify='no', update='yes')
SimSpec.iraf.hedit(SimSpec.WorkDir+'*_SQ???.fits', 'UT', value='(@"TIME-OBS")', add='yes', verify='no')
SimSpec.iraf.hedit(SimSpec.WorkDir+'*_SQ???.fits', 'ST', value='(@"LST")', add='yes', verify='no')
SimSpec.iraf.hedit(SimSpec.WorkDir+'*_SQ???.fits', 'airmass', value='(@"SECZ")', add='yes', verify='no')

#SimSpec.iraf.setairmass.unlearn()
#SimSpec.iraf.setairmass(SimSpec.WorkDir+'*_SQ???.fits', equinox="equinox")

#SimSpec.iraf.ccdproc(SimSpec.WorkDir+'*_SQ???.fits', output="", ccdtype="",  fixpix='no', oversca='no', trim='yes', zerocor='no', darkcor='no', flatcor='no', trimsec='[100:200,100:200]', zero='bias', flat='nflat')