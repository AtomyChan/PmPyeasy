SimSpec.iraf.hedit.unlearn()
SimSpec.iraf.hedit(SimSpec.WorkDir+'*_SQ???.fits', 'OBSERVAT', '(@"SITENAME")', add='yes', addonly='yes', verify='no')

SimSpec.iraf.setjd.unlearn()
SimSpec.iraf.setjd(SimSpec.WorkDir+'*_SQ???.fits', time="UT-TIME", epoch="EQUINOX")

#SimSpec.iraf.hedit(SimSpec.WorkDir+'*_SQ???.fits', 'ccdsec,datasec', delete='yes', verify='no', update='yes')
#SimSpec.iraf.hedit(SimSpec.WorkDir+'*_SQ???.fits', 'EPOCH', "(EQUINOX)", add='yes', addonly='yes', verify='no')

SimSpec.iraf.ccdproc.unlearn()
#SimSpec.iraf.ccdproc(SimSpec.WorkDir+'*_SQ???.fits', output="", ccdtype="",  fixpix='no', oversca='no', trim='yes', zerocor='no', darkcor='no', flatcor='no', trimsec='[100:200,100:200]', zero='bias', flat='nflat')

SimSpec.iraf.ccdproc(SimSpec.WorkDir+'*_SQ???.fits', output="", ccdtype="",  fixpix='no', oversca='yes', trim='no', zerocor='no', darkcor='no', flatcor='no', biassec="image", readaxis="column", trimsec='[*,860:3130]', zero='bias', flat='nflat')

#try:
#   SimSpec.iraf.ccdproc(SimSpec.WorkDir+'*_SQ???.fits', output="", ccdtype="",  fixpix='no', oversca='no', trim='yes', zerocor='no', darkcor='no', flatcor='no', biassec="image", readaxis="column", trimsec='[*,1508:2786]', zero='bias', flat='nflat')
#except:
#   pass

SimSpec.iraf.ccdproc.unlearn()
for fl in SimSpec.glob.glob(SimSpec.WorkDir+'*_SQ???.fits'):
    y=int(SimSpec.iraf.hselect(fl,'naxis'+('2' if SimSpec.param['dispaxis'].lower()=='x' else '1'),'yes',Stdout=1)[0])
    if y>SimSpec.param['cropheight']:
        dl=(y-SimSpec.param['cropheight'])/2
        tsec="%d:%d" % (dl+1,dl+SimSpec.param['cropheight'])
        if SimSpec.param['dispaxis'].lower()=='x':
            tsec='[*,'+tsec+']'
        else:
            tsec='['+tsec+',*]'
    else:
        tsec='[*,*]'
    SimSpec.iraf.ccdproc(fl, output="", ccdtype="",  fixpix='no', oversca='no', trim='yes', zerocor='no', darkcor='no', flatcor='no', biassec="image", readaxis="column", trimsec=tsec, zero='bias', flat='nflat')

import time
print ('Waiting a bit to allow the buggy IRAF to settle down. Sorry for that!!')
time.sleep(15)

SimSpec.iraf.hedit(SimSpec.WorkDir+'*_SQ???.fits', 'ccdsec,datasec', delete='yes', verify='no', update='yes')
