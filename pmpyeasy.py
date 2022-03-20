
class PmPyeasy:
    import ConfigParser
    import os
    import numpy as np
    import glob
    from time import time
    
    __name__='PmPyeasy'
    __ver__='2.4.4'
    __date__='1 Apr 2019'
    __author__='Subhash Bose'
    __description__='SIMplified SPECtrum reduction Interface'
    
    try:
        ProDir=os.path.dirname(os.path.realpath(__file__))
    except:
        ProDir=os.getcwd()
    
    ProDir=os.path.join(ProDir,'','')
    ConfDir=ProDir+'modules/'

    #if not os.path.exists(ConfDir):
    #    ConfDir = os.path.join(os.path.realpath(ProDir + '../.'+__name__+'Modules/'),'')

    
    def __init__(self, inifile='setup.conf', Object=False, DataDir='', WorkDir=''):

        from pyraf import iraf
        self.iraf=iraf

        self.iraf.imred(_doprint=0)
        self.iraf.ccdred(_doprint=0)
        self.iraf.crutil(_doprint=0)
        # self.iraf.nproto(_doprint=0)
        # self.iraf.stsdas(_doprint=0)
        self.iraf.twodspec(_doprint=0)
        self.iraf.apextract(_doprint=0)
        self.iraf.onedspec(_doprint=0)
        self.iraf.ccdred.instrument = self.ProDir + 'extern/blank.dat'
        # self.iraf.task(ss_lacos=ProDir+'extern/lacos_spec.cl')

        #print inifile
        self.param={}
        self.BiasFiles=[]
        self.FlatFiles=[]
        self.StdFiles=[]
        self.ArcFiles=[]
        self.ObjectFiles=[]
        self.Files={}
        self.inifile=inifile
        if not self.os.path.isfile(inifile):
            raise self.PmPyeasyError("Unable to find configuration file: "+inifile)
        self.iniconf = self.ConfigParser.RawConfigParser(allow_no_value=True)
        self.iniconf.read(inifile)
        self.Caseiniconf = self.ConfigParser.RawConfigParser(allow_no_value=True)
        self.Caseiniconf.optionxform = str
        self.Caseiniconf.read(inifile)
        if DataDir=='':
            self.DataDir=self.os.path.dirname(self.inifile)
        else:
            self.DataDir=DataDir
        self.DataDir=self.os.path.join(self.DataDir,'')
        if WorkDir=='':
            self.WorkDir=self.os.path.join(self.DataDir,'ImgReduction')
        else:
            self.WorkDir=WorkDir
        self.WorkDir=self.os.path.join(self.WorkDir,'')
        self.parseSetup()
        self.DbDir=self.param['fileprefix']+'db'
        
        #self.StdDir=self.os.path.join(self.os.path.realpath(self.param['stddir']),'')
        
        if Object:
            self.parseObject(Object)
        else:
            self.ObjectGrp=False
            self.StdName=False

    class PmPyeasyError(Exception):
        __module__ = Exception.__module__

    def parseObject(self, Object):
        self.ObjectGrp=Object
        
        self.parseFlatBias(Object)
        
        try:
            self.StdGrp=self.iniconf.get(Object,'std')
        except:
            self.StdGrp='std'
        self.StdFiles=self.parseGrp(self.StdGrp)
        self.Files['StdSQ']=self.WorkDir+self.param['fileprefix']+self.StdGrp+'_SQ'
        self.Files['Std2D']=self.WorkDir+self.param['fileprefix']+self.StdGrp+'.fits'
        self.Files['Std1D']=self.WorkDir+self.param['fileprefix']+self.StdGrp+'.0001.fits'
        self.Files['StdWcal']=self.WorkDir+self.param['fileprefix']+self.StdGrp+'.0001.w.fits'
        self.Files['StdFcal']=self.WorkDir+self.param['fileprefix']+self.StdGrp+'.0001.f.fits'
        self.Files['StdSens']=self.WorkDir+self.param['fileprefix']+self.StdGrp+'_sens.fits'
        self.Files['StdStdcal']=self.WorkDir+self.param['fileprefix']+self.StdGrp+'_stdcal'
        
        if self.StdFiles:
            try:
                self.StdName=self.iniconf.get(self.StdGrp,'name')
            except:
                raise self.PmPyeasyError('No name defined for "'+self.StdGrp+'" group standard.')
        
        try:
            self.ArcGrp=self.iniconf.get(Object,'arc')
        except:
            self.ArcGrp='arc'
        self.ArcFiles=self.parseGrp(self.ArcGrp)
        self.Files['ArcSQ']=self.WorkDir+self.param['fileprefix']+self.ArcGrp+'_SQ'
        self.Files['Arc2D']=self.WorkDir+self.param['fileprefix']+self.ArcGrp+'.fits'
        self.Files['Arc1D']=self.WorkDir+self.param['fileprefix']+self.ArcGrp+'.0001.fits'
        
        self.ObjectFiles=self.parseGrp(self.ObjectGrp)
        self.Files['ObjectSQ']=self.WorkDir+self.param['fileprefix']+self.ObjectGrp+'_SQ'
        self.Files['Object2D']=self.WorkDir+self.param['fileprefix']+self.ObjectGrp+'.fits'
        self.Files['Object1D']=self.WorkDir+self.param['fileprefix']+self.ObjectGrp+'.0001.fits'
        self.Files['ObjectWcal']=self.WorkDir+self.param['fileprefix']+self.ObjectGrp+'.0001.w.fits'
        self.Files['ObjectFcal']=self.WorkDir+self.param['fileprefix']+self.ObjectGrp+'.0001.f.fits'
        if not self.ObjectFiles:
            raise self.PmPyeasyError('No files listed for the "'+Object+'" group object.')
    
    def parseFlatBias(self, GrpName):
        if GrpName!=self.ObjectGrp:
            self.parseFlatBias(self.ObjectGrp)
            
        try:
            self.BiasGrp=self.iniconf.get(GrpName,'bias')
        except:
            if GrpName==self.ObjectGrp:
                self.BiasGrp='bias'                
        self.BiasFiles=self.parseGrp(self.BiasGrp)
        self.Files['BiasSQ']=self.WorkDir+self.param['fileprefix']+self.BiasGrp+'_SQ'
        self.Files['Bias']=self.WorkDir+self.param['fileprefix']+self.BiasGrp+'.fits'
        
        try:
            self.FlatGrp=self.iniconf.get(GrpName,'flat')
        except:
            if GrpName==self.ObjectGrp:
                self.FlatGrp='flat'
        self.FlatFiles=self.parseGrp(self.FlatGrp)
        self.Files['FlatSQ']=self.WorkDir+self.param['fileprefix']+self.FlatGrp+'_SQ'
        self.Files['Flat']=self.WorkDir+self.param['fileprefix']+self.FlatGrp+'.fits'
        self.Files['FlatNorm']=self.WorkDir+self.param['fileprefix']+self.FlatGrp+'_norm.fits'
        
    
    def parseGrp(self,Grp):
        lis=[]
        try:
            for k,l in self.Caseiniconf.items(Grp):
                if k.lower() not in ['arc','bias','flat','std','name']:
                    lis=lis+[self.os.path.join(self.DataDir,k)]
        except:
            #lis=False
            pass
        return lis
    
    def parseParam(self,k,v):
        if k in ['ref-arc','ref-id','precleanscript','postcleanscript']:
            v=self.os.path.join(self.os.path.dirname(self.inifile), v)
        elif k in ['extinction', 'line-list','stddir']:
            if '{prodir}' in v:
                v=v.replace('{prodir}',self.ProDir)
            if not '$' in v:
                v=self.os.path.realpath(self.os.path.join(self.os.path.dirname(self.inifile), v))
                if k=='stddir':
                    v=self.os.path.join(v,'')
        if str(v).lower()=='true':
            v=True
        elif str(v).lower()=='false':
            v=False
        else:
            try:
                v=float(v)
            except:
                v=v.split(';')[0]
        self.param[k]=v
    
    def parseSetup(self):
        try:
            self.Instrument=self.iniconf.get('setup','instrument')
        except:
            self.Instrument='default'
        if self.inifile!=self.ConfDir+'default/setup.conf':
            if self.os.path.exists(self.ConfDir+self.Instrument+'/setup.conf'):
                self.param.update(PmPyeasy(self.ConfDir+self.Instrument+'/setup.conf').param)
            else:
                raise self.PmPyeasyError('Module '+self.Instrument+' does not exist')
        for k,v in self.iniconf.items('setup'):
            self.parseParam(k,v)
        #print self.inifile
    
    def ListObjects(self):
        lis=[]
        for l in self.iniconf.sections():
	    if ('bias' not in l.lower()) and ('flat' not in l.lower()) and ('setup' not in l.lower()):
                lis=lis+[l]
        return(lis)
    
    def DeleteOldFiles(self,files,timestamp='now'):
        if timestamp:
            if str(timestamp).lower()=='now':
                timestamp=self.time()
            for f in self.glob.glob(files):
                if self.os.stat(f).st_mtime<timestamp:
                    self.os.remove(f)
    
    def isInputNew(self, inp, out, delout=False):
        if self.os.path.isfile(out):
            if self.os.path.isfile(inp) and self.os.stat(inp).st_ctime>self.os.stat(out).st_ctime+0.1:
                if delout:
                    self.os.remove(out)
                return True
            else:
                return False
        else:
            return True
    
    def LsStdDir(self):
        v=self.param['stddir']
        if '$' in v:
            irdir=v.split('$')[0]+'$'
            uxdir='/'+str(self.iraf.path(irdir,Stdout=1)[0]).split('!/')[-1]
            v=v.replace(irdir,uxdir)
        v=self.os.path.join(v,'')
        lis=[]
        for k in self.glob.glob(v+'*.dat'):
            lis.append(k.replace(v,'').replace('.dat',''))
        if self.os.path.isfile(v+'names.men'):
            with open(v+'names.men') as f:
                for row in f:
                    row=row.strip()
                    if row=='' or row[0]=='#':
                        continue
                    lis.append(row.split()[0])
        return(sorted(lis, key=lambda s: s.lower()))
    
    def GetFname(self,name):
        res=self.os.path.basename(name).rsplit('.', 1)[0]
        #if len(res)==0:
        #    res=name
        return(res)
    
    def CombineBiasFlat(self,DelOlder=False):

    	from collections import OrderedDict
        self.iraf.imdelete.unlearn()
        self.iraf.imcopy.unlearn()
        
        self.DeleteOldFiles(self.Files['BiasSQ']+'???.fits',DelOlder)
        self.DeleteOldFiles(self.Files['FlatSQ']+'???.fits',DelOlder)
        self.DeleteOldFiles(self.Files['Bias'],DelOlder)
        self.DeleteOldFiles(self.Files['Flat'],DelOlder)
        self.DeleteOldFiles(self.Files['FlatNorm'],DelOlder)
            
        if self.param['biascor'] and self.BiasFiles and not self.os.path.isfile(self.Files['Bias']):
            self.iraf.imdelete(self.Files['BiasSQ']+'???.fits')
            for i in range(1,len(self.BiasFiles)+1):
                self.iraf.imcopy(self.BiasFiles[i-1]+self.param['fitsext'], self.Files['BiasSQ']+("%03d" % i)+'.fits')
            self.PreCleanExec(self.Files['BiasSQ']+'???.fits')
            self.iraf.zerocombine.unlearn()
            self.iraf.zerocombine(self.Files['BiasSQ']+'???.fits', output=self.Files['Bias'], combine='median', reject='ccdclip', ccdtype="", delete='yes' if self.param['deletecalsingles'] else 'no', rdnoise=self.param['readnoise'], gain=self.param['gain'])

        if self.param['flatcor'] and self.FlatFiles and not self.os.path.isfile(self.Files['Flat']):
            self.iraf.imdelete(self.Files['FlatSQ']+'???.fits')
            for i in range(1,len(self.FlatFiles)+1):
                self.iraf.imcopy(self.FlatFiles[i-1]+self.param['fitsext'], self.Files['FlatSQ']+("%03d" % i)+'.fits')
            self.PreCleanExec(self.Files['FlatSQ']+'???.fits')
            self.iraf.ccdproc.unlearn()
            self.iraf.ccdproc(self.Files['FlatSQ']+'???.fits', output="", ccdtype="", fixpix='no', overscan='no', trim='no', zerocor='yes' if self.param['biascor'] else 'no', darkcor='no', flatcor='no',trimsec='', zero=self.Files['Bias'], biassec='[*,*]')
            self.iraf.flatcombine.unlearn()
            self.iraf.flatcombine(self.Files['FlatSQ']+'???.fits', output=self.Files['Flat'], combine='median', reject='avsigclip', ccdtype="", process='no', subsets='no', delete='yes' if self.param['deletecalsingles'] else 'no', scale='mode', rdnoise=self.param['readnoise'], gain=self.param['gain'])
            self.iraf.imreplace.unlearn()
            self.iraf.imreplace(self.Files['Flat'], value=1, upper=0)


	    self.iraf.imstatistics.unlearn()	
	    imgstatfields = "image,npix,mode, midpt,mean,stddev,min,max"
	    imgstatsec = self.param['nflat.statsec']
	    imgstatout = self.iraf.imstatistics(self.Files['Flat']+imgstatsec, fields = imgstatfields, lower = 'INDEF', upper = 'INDEF',nclip = 0,lsigma = 3,usigma =3,binwidth = 1, Stdout=1)
	    imgstatvalues = imgstatout[1].split()
	    stat_retdict = OrderedDict()
	    for i,field in enumerate(imgstatfields.split(',')):
		stat_retdict[field] = imgstatvalues[i]	

	    self.iraf.imarith.unlearn()
	    if self.param['nflat.method'] in ['mean', 'max', 'mode']:
		operand2 = stat_retdict[self.param['nflat.method']]
	    else:
		raise self.PmPyeasyError('flat normalization method %s not valid'%self.param['nflat.method'])

	    print self.Files['Flat']
	    print operand2
	    print self.Files['FlatNorm']
	    self.iraf.imarith(operand1=self.Files['Flat'], op='/', operand2=operand2, result=self.Files['FlatNorm'])



    def CosmicClean(self,img):
        from extern import cosmics
        def subf():
            array, header = cosmics.fromfits(img)
            c = cosmics.cosmicsimage(array, gain=self.param['gain'],
                    readnoise=self.param['readnoise'],
                    sigclip = self.param['lacos.sigclip'], sigfrac = self.param['lacos.sigfrac'],
                    objlim = self.param['lacos.objlim'], satlevel=self.param['ap.saturation'],
                    skyOrder = int(self.param['lacos.skyorder']), objectOrder = int(self.param['lacos.objorder']))
            c.run(maxiter = int(self.param['lacos.niter']))
            header.set('COSMIC', 'cleaned')
            cosmics.tofits(img+".cos.fits", c.cleanarray, header)
        subf()
        if self.os.path.isfile(img+".cos.fits"):
            self.os.system('mv -f "'+img+'.cos.fits" "'+img+'"')
        
    def PreCleanExec(self,files):
        self.iraf.hedit(files, fields="Red-Soft", value=self.__name__+' v'+self.__ver__+' ', add='yes', addonly='yes', ver='no', show='no');
        self.iraf.hedit(files, fields="SSModule", value=str(self.param['modulever'])+' ', add='yes', addonly='yes', ver='no', show='no');
        if 'reducer' in self.param and self.param['reducer']!='':
            self.iraf.hedit(files, fields="Reducer", value=self.param['reducer'], add='yes', addonly='yes', ver='no', show='no');
        execfile(self.param['precleanscript'],{},{'PmPyeasy':self})

    def Clean(self,DelOlder=False):
        self.iraf.imdelete.unlearn()
        self.iraf.imcopy.unlearn()
        
        if not self.os.path.exists(self.WorkDir):
            self.os.makedirs(self.WorkDir)
        

        self.DeleteOldFiles(self.Files['ArcSQ']+'???.fits',DelOlder)
        self.DeleteOldFiles(self.Files['Arc2D'],DelOlder)
        self.DeleteOldFiles(self.Files['StdSQ']+'???.fits',DelOlder)
        self.DeleteOldFiles(self.Files['Std2D'],DelOlder)
        self.DeleteOldFiles(self.Files['ObjectSQ']+'???.fits',DelOlder)
        self.DeleteOldFiles(self.Files['Object2D'],DelOlder)
        
        self.iraf.imcombine.unlearn()
        
        #if self.param['cosmiccor']:
            #self.iraf.ss_lacos.unlearn()
            #self.iraf.ss_lacos.gain    = self.param['gain']            #gain (electrons/ADU)
            #self.iraf.ss_lacos.readn   = self.param['readnoise']            #read noise (electrons)
            #self.iraf.ss_lacos.xorder  = self.param['lacos.xorder']              #order of object fit (0=no fit)
            #self.iraf.ss_lacos.yorder  = self.param['lacos.yorder']             #order of sky line fit (0=no fit)
            #self.iraf.ss_lacos.sigclip = self.param['lacos.sigclip']            #detection limit for cosmic rays (sigma)
            #self.iraf.ss_lacos.sigfrac = self.param['lacos.sigfrac']            #fractional detection limit for neighbouring pixels
            #self.iraf.ss_lacos.objlim  = self.param['lacos.objlim']             #contrast limit between CR and underlying object
            #self.iraf.ss_lacos.niter   = self.param['lacos.niter']              #maximum number of iterations
        
        
        self.iraf.ccdproc.unlearn()
        if self.ArcFiles and not self.os.path.isfile(self.Files['Arc2D']):
	    print "Clean Arc files..."
            self.iraf.imdelete(self.Files['ArcSQ']+'???.fits')
            self.parseFlatBias(self.ArcGrp)
            self.CombineBiasFlat(DelOlder)
            for i in range(1,len(self.ArcFiles)+1):
                self.iraf.imcopy(self.ArcFiles[i-1]+self.param['fitsext'], self.Files['ArcSQ']+("%03d" % i)+'.fits')
            self.PreCleanExec(self.Files['ArcSQ']+'???.fits')
            self.iraf.ccdproc(self.Files['ArcSQ']+'???.fits', output="", ccdtype="", fixpix='no', oversca='no', trim='yes' if self.param['trimsec']!='[*,*]' else 'no', zerocor='yes' if self.param['biascor'] else 'no', darkcor='no', flatcor='yes' if self.param['flatcor'] else 'no', trimsec=self.param['trimsec'], zero=self.Files['Bias'], flat=self.Files['FlatNorm'])
            print self.Files['ArcSQ']+'???.fits'    
            self.iraf.imcombine(self.Files['ArcSQ']+'???.fits', output=self.Files['Arc2D'], comb=self.param['combineobjtype'] )
            if self.param['deleteobjsingles']:
                self.iraf.imdelete(self.Files['ArcSQ']+'???.fits')
        
        
        if self.StdFiles and not self.os.path.isfile(self.Files['Std2D']):
	    print "Clean Std files..."
            self.iraf.imdelete(self.Files['StdSQ']+'???.fits')
            self.parseFlatBias(self.StdGrp)
            self.CombineBiasFlat(DelOlder)
            for i in range(1,len(self.StdFiles)+1):
                self.iraf.imcopy(self.StdFiles[i-1]+self.param['fitsext'], self.Files['StdSQ']+("%03d" % i)+'.fits')
            self.PreCleanExec(self.Files['StdSQ']+'???.fits')
            self.iraf.ccdproc(self.Files['StdSQ']+'???.fits', output="", ccdtype="", fixpix='no', oversca='no', trim='yes' if self.param['trimsec']!='[*,*]' else 'no', zerocor='yes' if self.param['biascor'] else 'no', darkcor='no', flatcor='yes' if self.param['flatcor'] else 'no', trimsec=self.param['trimsec'], zero=self.Files['Bias'], flat=self.Files['FlatNorm'])
            
            if self.param['cosmiccor']:
                for i in range(1,len(self.StdFiles)+1):
                    self.CosmicClean(self.Files['StdSQ']+("%03d" % i)+'.fits')
        
            self.iraf.imcombine(self.Files['StdSQ']+'???.fits', output=self.Files['Std2D'], comb=self.param['combineobjtype'] )
            if self.param['deleteobjsingles']:
                self.iraf.imdelete(self.Files['StdSQ']+'???.fits')
        
        
        if not self.os.path.isfile(self.Files['Object2D']):
            self.iraf.imdelete(self.Files['ObjectSQ']+'???.fits')
            self.parseFlatBias(self.ObjectGrp)
            self.CombineBiasFlat(DelOlder)
            for i in range(1,len(self.ObjectFiles)+1):
                self.iraf.imcopy(self.ObjectFiles[i-1]+self.param['fitsext'], self.Files['ObjectSQ']+("%03d" % i)+'.fits')
            self.PreCleanExec(self.Files['ObjectSQ']+'???.fits')
	    
            self.iraf.ccdproc(self.Files['ObjectSQ']+'???.fits', output="", ccdtype="", fixpix='no', oversca='no', trim='no', zerocor='yes' if self.param['biascor'] else 'no', darkcor='no', flatcor='yes' if self.param['flatcor'] else 'no', trimsec='', zero=self.Files['Bias'], flat=self.Files['FlatNorm'], minrepl=0.1)
            

	    print self.param['trimsec']
            self.iraf.ccdproc(self.Files['ObjectSQ']+'???.fits', output="", ccdtype="", fixpix='no', oversca='no', trim='yes' if self.param['trimsec']!='[*,*]' else 'no', zerocor='no', darkcor='no', flatcor='no', trimsec=self.param['trimsec'], zero='', flat='', minrepl=0.1)

	    print "Clean Obj image done!"
            if self.param['cosmiccor']:
                for i in range(1,len(self.ObjectFiles)+1):
                    self.CosmicClean(self.Files['ObjectSQ']+("%03d" % i)+'.fits')
            
            self.iraf.imcombine(self.Files['ObjectSQ']+'???.fits', output=self.Files['Object2D'], comb=self.param['combineobjtype'] )
            if self.param['deleteobjsingles']:
                self.iraf.imdelete(self.Files['ObjectSQ']+'???.fits')
        execfile(self.param['postcleanscript'],{},{'PmPyeasy':self})
    
    def ApExtract(self, spectype='object',DelOlder=False):
        self.iraf.apextract.unlearn()
        self.iraf.apextract.database=self.WorkDir+self.DbDir
        self.iraf.apextract.dispaxi   = 1 if self.param['dispaxis'].lower()=='x' else 2 
        self.iraf.apall.unlearn()
        self.iraf.apall.nfind         = 1              #Number of apertures to be found automatically
        self.iraf.apall.output        = ""             #List of output spectra
        self.iraf.apall.apertures     = ""             #Apertures
        self.iraf.apall.format        = "onedspec"     #Extracted spectra format
        self.iraf.apall.references    = ""             #List of aperture reference images
        self.iraf.apall.profiles      = ""             #List of aperture profile images\n
        self.iraf.apall.interactive   = 'yes'            #Run task interactively?
        self.iraf.apall.find          = 'yes'            #Find apertures?
        self.iraf.apall.recenter      = 'no'             #Recenter apertures?
        self.iraf.apall.resize        = 'no'             #Resize apertures?
        self.iraf.apall.edit          = 'yes'            #Edit apertures?
        self.iraf.apall.trace         = 'yes'            #Trace apertures?
        self.iraf.apall.fittrace      = 'yes'            #Fit the traced points interactively?
        self.iraf.apall.extract       = 'yes'            #Extract spectra?
        self.iraf.apall.extras        = 'yes'            #Extract sky, sigma, etc.?
        self.iraf.apall.review        = 'yes'            #Review extractions?\n
        self.iraf.apall.line          = self.param['ap.line']           #Dispersion line
        self.iraf.apall.nsum          = self.param['ap.nsum']             #Number of dispersion lines to sum or median\n\n# DEFAULT APERTURE PARAMETERS\n
        self.iraf.apall.lower         = self.param['ap.lower']           #Lower aperture limit relative to center
        self.iraf.apall.upper         = self.param['ap.upper']            #Upper aperture limit relative to center
        self.iraf.apall.apidtable     = ""             #Aperture ID table (optional)\n\n# DEFAULT BACKGROUND PARAMETERS\n
        self.iraf.apall.b_function    = self.param['ap.b_function']    #Background function
        self.iraf.apall.b_order       = self.param['ap.b_order']              #Background function order
        self.iraf.apall.b_sample      = self.param['ap.b_sample']    #Background sample regions
        self.iraf.apall.b_naverage    = -4           #Background average or median
        self.iraf.apall.b_niterate    = 0              #Background rejection iterations
        self.iraf.apall.b_low_reject  = 3.             #Background lower rejection sigma
        self.iraf.apall.b_high_rejec  = 3.             #Background upper rejection sigma
        self.iraf.apall.b_grow        = 0.             #Background rejection growing radius\n\n# APERTURE CENTERING PARAMETERS\n
        self.iraf.apall.width         = self.param['ap.width']             #Profile centering width
        self.iraf.apall.radius        = self.param['ap.radius']             #Profile centering radius
        self.iraf.apall.threshold     = 0.             #Detection threshold for profile centering\n\n# AUTOMATIC FINDING AND ORDERING PARAMETERS\n
        self.iraf.apall.minsep        = 5.             #Minimum separation between spectra
        self.iraf.apall.maxsep        = 1000.          #Maximum separation between spectra
        self.iraf.apall.order         = "increasing"   #Order of apertures\n\n# RECENTERING PARAMETERS\n
        self.iraf.apall.aprecenter    = ""             #Apertures for recentering calculation
        self.iraf.apall.npeaks        = 'INDEF'          #Select brightest peaks
        self.iraf.apall.shift         = 'yes'            #Use average shift instead of recentering?\n\n# RESIZING PARAMETERS\n
        self.iraf.apall.llimit        = 'INDEF'          #Lower aperture limit relative to center
        self.iraf.apall.ulimit        = 'INDEF'          #Upper aperture limit relative to center
        self.iraf.apall.ylevel        = 0.1            #Fraction of peak or intensity for automatic width
        self.iraf.apall.peak          = 'yes'            #Is ylevel a fraction of the peak?
        self.iraf.apall.bkg           = 'yes'            #Subtract background in automatic width?
        self.iraf.apall.r_grow        = 0.             #Grow limits by this factor
        self.iraf.apall.avglimits     = 'no'             #Average limits over all apertures?\n\n# TRACING PARAMETERS\n
        self.iraf.apall.t_nsum        = 12             #Number of dispersion lines to sum
        self.iraf.apall.t_step        = 8              #Tracing step
        self.iraf.apall.t_nlost       = 10             #Number of consecutive times profile is lost before quitting 3 default give large fro low signal
        self.iraf.apall.t_function    = "chebyshev"    #Trace fitting function
        self.iraf.apall.t_order       = self.param['ap.t_order']              #Trace fitting function order
        self.iraf.apall.t_sample      = "*"            #Trace sample regions
        self.iraf.apall.t_naverage    = -10            #Trace average or median give large -ve number for median
        self.iraf.apall.t_niterate    = 0              #Trace rejection iterations
        self.iraf.apall.t_low_reject  = 3.             #Trace lower rejection sigma
        self.iraf.apall.t_high_rejec  = 3.             #Trace upper rejection sigma
        self.iraf.apall.t_grow        = 0.             #Trace rejection growing radius\n\n# EXTRACTION PARAMETERS\n
        self.iraf.apall.background    = "fit"          #Background to subtract
        self.iraf.apall.skybox        = 1              #Box car smoothing length for sky
        self.iraf.apall.weights       = "variance"     #Extraction weights (none|variance)
        self.iraf.apall.pfit          = "fit1d"        #Profile fitting type (fit1d|fit2d)
        self.iraf.apall.clean         = 'yes'            #Detect and replace bad pixels?
        self.iraf.apall.saturation    = self.param['ap.saturation']         #Saturation level
        self.iraf.apall.readnoise     = self.param['readnoise']          #Read out noise sigma (photons)
        self.iraf.apall.gain          = self.param['gain']         #Photon gain (photons/data number)
        self.iraf.apall.lsigma        = 3.             #Lower rejection threshold
        self.iraf.apall.usigma        = 3.             #Upper rejection threshold
        self.iraf.apall.nsubaps       = 1              #Number of subapertures per aperture
        
        if spectype.lower()=='object':
            self.DeleteOldFiles(self.Files['Object1D'],DelOlder)
            if self.isInputNew(self.Files['Object2D'], self.Files['Object1D'], delout=True):
                self.iraf.apall(self.Files['Object2D'])
        elif spectype.lower()=='standard':
            self.DeleteOldFiles(self.Files['Std1D'],DelOlder)
            if self.isInputNew(self.Files['Std2D'],self.Files['Std1D'], delout=True):
                self.iraf.apall(self.Files['Std2D'])
        elif spectype.lower()=='arc':
            self.DeleteOldFiles(self.Files['Arc1D'],DelOlder)
            if self.isInputNew(self.Files['Arc2D'],self.Files['Arc1D'], delout=True):
                self.iraf.apall(self.Files['Arc2D'], ref=self.Files['Object2D'], extras='no', interac='no', recen='no', resize='no', trace='no', line='INDEF', background="none")
        else:
            raise self.PmPyeasyError('Invalid SpecType '+spectype)
    
    def AutoId(self,DelOlder=False):
        self.DeleteOldFiles(self.WorkDir+self.DbDir+'/id'+self.GetFname(self.Files['Arc1D']),DelOlder)
        if self.isInputNew( self.Files['Arc1D'], self.WorkDir+self.DbDir+'/id'+self.GetFname(self.Files['Arc1D']), delout=True):
            refcopy=self.WorkDir+'refarc_00.0001.fits'
            idcopy=self.WorkDir+self.DbDir+'/idrefarc_00.0001'
            import shutil
            shutil.copy2(self.param['ref-arc'], refcopy)
            shutil.copy2(self.param['ref-id'], idcopy)
            self.os.system("sed -i -e 's/"+self.GetFname(self.param['ref-arc'])+"/"+self.GetFname(refcopy)+"/g' "+idcopy)

            self.iraf.reidentify.unlearn()
            cwd=self.os.getcwd()
            self.os.chdir( self.WorkDir )
            try:
                self.iraf.reidentify(reference=self.os.path.basename(refcopy),images=self.os.path.basename(self.Files['Arc1D']),coordli=self.param['line-list'],cradius=self.param['id.cradius'],thresho=0,interactive='no',shift='INDEF',search='INDEF',refit='yes',database=self.DbDir)
            finally:
                self.os.chdir( cwd )
            
            self.os.remove(refcopy)
            self.os.remove(idcopy)
            #self.DeleteOldFiles(refcopy,self.time()+10)
            #self.DeleteOldFiles(idcopy,self.time()+10)
            return(True)
        return(False)
        
    def MannId(self):
        self.iraf.identify.unlearn()
        cwd=self.os.getcwd()
        #if self.os.path.isfile(self.WorkDir + self.DbDir + '/id' + self.GetFname(self.Files['Arc1D'])) \
        #        and raw_input("Start over fresh and discard old line identification? [y/N] ").strip().lower()=='y':
        #    self.os.remove(self.WorkDir + self.DbDir + '/id' + self.GetFname(self.Files['Arc1D']))
        self.os.chdir( self.WorkDir )
        try:
            self.iraf.identify(images=self.os.path.basename(self.Files['Arc1D']),coordlist=self.param['line-list'],cradius=self.param['id.cradius'], fwidth=self.param['id.fwidth'],thresho=0, function='spline3',order=4,database=self.DbDir)
        finally:
            self.os.chdir( cwd )
    
    def DispCor(self,DelOlder=False):
        self.iraf.hedit.unlearn()
        self.iraf.dispcor.unlearn()
        
        self.DeleteOldFiles(self.Files['StdWcal'],DelOlder)
        if self.StdFiles and ( self.isInputNew( self.Files['Std1D'], self.Files['StdWcal'], delout=True) or self.isInputNew( self.Files['Arc1D'], self.Files['StdWcal'], delout=True) ):
            self.iraf.hedit(self.Files['Std1D'], fields="REFSPEC1", value=self.os.path.basename(self.Files['Arc1D']), add='yes', addonly='yes', ver='no', show='no')
            cwd=self.os.getcwd()
            self.os.chdir( self.WorkDir )
            try:
                self.iraf.dispcor(self.os.path.basename(self.Files['Std1D']), self.os.path.basename(self.Files['StdWcal']), w1=self.param['dispc.w1'], w2=self.param['dispc.w2'], confirm='yes' if str(self.param['dispc.w1']).lower()=='indef' and str(self.param['dispc.w2']).lower()=='indef' else 'no', database=self.DbDir)
            finally:
                self.os.chdir( cwd )
        
        self.DeleteOldFiles(self.Files['ObjectWcal'],DelOlder)
        if self.isInputNew( self.Files['Object1D'], self.Files['ObjectWcal'], delout=True) or self.isInputNew( self.Files['Arc1D'], self.Files['ObjectWcal'], delout=True):
            from astropy.io import fits
            try:
                with fits.open(self.Files['StdWcal']) as hdul:
                    header=hdul[0].header
                w1=header['CRVAL1']
                w2=w1+(float(header['NAXIS1'])-1)*float(header['CD1_1'])
            except:
                w1=self.param['dispc.w1']
                w2=self.param['dispc.w2']

            self.iraf.hedit(self.Files['Object1D'], fields="REFSPEC1", value=self.os.path.basename(self.Files['Arc1D']), add='yes', addonly='yes', ver='no', show='no')
            cwd=self.os.getcwd()
            self.os.chdir( self.WorkDir )
            try:
                self.iraf.dispcor(self.os.path.basename(self.Files['Object1D']), self.os.path.basename(self.Files['ObjectWcal']), w1=w1, w2=w2, confirm='no', database=self.DbDir)
            finally:
                self.os.chdir( cwd )

    def FluxCal(self,DelOlder=False):
        if self.StdFiles:
            self.DeleteOldFiles(self.Files['StdSens'],DelOlder)
            if self.isInputNew( self.Files['StdWcal'], self.Files['StdSens'], delout=True):
                self.DeleteOldFiles(self.Files['StdStdcal'])
                self.iraf.standard.unlearn()
                self.iraf.sensfunc.unlearn()
                cwd=self.os.getcwd()
                self.os.chdir( self.WorkDir )
                try:
                    self.iraf.standard(self.os.path.basename(self.Files['StdWcal']), self.os.path.basename(self.Files['StdStdcal']), extinct=self.param['extinction'], caldir=self.param['stddir'], star_name=self.StdName, answer="no")

                    self.iraf.sensfunc(self.os.path.basename(self.Files['StdStdcal']), self.os.path.basename(self.Files['StdSens']), ignoreaps="yes", extinct=self.param['extinction'])
                finally:
                    self.os.chdir( cwd )
        
        for stype in ['Std','Object']:
            self.iraf.calibrate.unlearn()
            self.DeleteOldFiles(self.Files[stype+'Fcal'],DelOlder)
            if self.isInputNew( self.Files[stype+'Wcal'], self.Files[stype+'Fcal'], delout=True) or  self.isInputNew( self.Files['StdSens'], self.Files[stype+'Fcal'], delout=True):
                cwd=self.os.getcwd()
                self.os.chdir( self.WorkDir )
                try:
                    self.iraf.calibrate(self.os.path.basename(self.Files[stype+'Wcal']), self.os.path.basename(self.Files[stype+'Fcal']), extinction=self.param['extinction'],sensitivity=self.os.path.basename(self.Files['StdSens']), ignoreaps="yes")
                finally:
                    self.os.chdir( cwd )
    
    def SkyLineCor(self,stype='Object'):
        if self.param['skylinecor'].lower()=='oi':
            rest_wl=5577.3387
        elif self.param['skylinecor'].lower()=='hgi':
            rest_wl=4358.34
        else:
            print('Invalid skyline for correction: '+self.param['skylinecor'])
            return
        
        f=self.Files[stype+'Fcal']
        from astropy.io import fits
        v=fits.getdata(f, ext=0)[2][0]
        with fits.open(f) as hdul:
            header=hdul[0].header
        import numpy as np
        w=np.arange(0,float(header['NAXIS1']))*float(header['CD1_1'])+header['CRVAL1']
        tmp=(w>rest_wl-30) & (w<rest_wl+30)
        v=v[tmp]
        w=w[tmp]
        cen=w[v==np.max(v)][0]
        
        
        shift=rest_wl-cen;

        print("Appling Shift: %.3f Ang" % shift)
        cwd=self.os.getcwd()
        self.os.chdir( self.WorkDir )
        try:
            self.iraf.specshift(self.os.path.basename(f), shift=shift, verbose='yes')
        finally:
            self.os.chdir( cwd )

        ShOld=self.iraf.hselect(f,"SB_SHIFT",'yes',Stdout=1)[0]

        try:
            ShOld=float(ShOld);
        except:
            ShOld=0.0


        shift=ShOld+shift
        
        self.param['skyshift']=shift

        self.iraf.hedit(f, fields="SB_SHIFT", value=shift, add='yes', addonly='yes', ver='no', show='no');

        print ("Net Shift applied = "+str(shift))

class SSutil:
    
    @staticmethod
    def ListInstrument():
        lis=[]
        for ini in PmPyeasy.glob.glob(PmPyeasy.ConfDir+'*/setup.conf'):
            inst=PmPyeasy.os.path.basename(PmPyeasy.os.path.dirname(ini))
            if inst.lower()=='default':
                continue
            tmp = PmPyeasy.ConfigParser.RawConfigParser(allow_no_value=True)
            tmp.read(ini)
            try:
                desc=tmp.get('setup','description')
            except:
                desc=''
            lis=lis+[(inst,desc)]
        return (sorted(lis))
    
    @staticmethod
    def ListDefaultStds(conf=PmPyeasy.ConfDir+'default/setup.conf'):
        return(PmPyeasy(conf).LsStdDir())
    
    @staticmethod
    def ObsLog(inst='default',lsinp='*.fits', datadir='', outlogfile='obslogfile'):
        if type(lsinp)==str:
            lsinp=[lsinp]
        try:
            p=PmPyeasy(PmPyeasy.ConfDir+inst+'/setup.conf').param
        except PmPyeasy.PmPyeasyError as zz:
            print "Invalid instrument name"
            return
        import sys, glob, os
        from astropy.io import fits
	from astropy.table import Table
        from collections import OrderedDict

        index=p['oblg.fitindex']
        keywords=p['oblg.keys']

        head=keywords.split(',');
        #print ' \t'.join([' File    ']+head);
        split=index.split('[',1)
        
        if (len(split)>1):
            index='['+split[1]
            indexn=index[1:-1]
            try:
                indexn=int(indexn)
            except: pass
        else:
            index=''
            indexn=0

	logdict = OrderedDict()
	logdict['File'] = []
	for headkey in head:
		logdict[headkey] = []
	
        for lsinp_ech in lsinp:
            for ls in sorted(glob.glob(os.path.join(datadir,lsinp_ech))):
                row=[ls+index];
		logdict['File'].append(ls+index)
                try:
                    hdulist = fits.open(ls, ignore_missing_end=True)[indexn]
                except:
                    print("Failed reading fits file "+ls)
                    continue
                for key in head:
                    try:
                        if isinstance(hdulist.header[key],float):
                            val=repr(hdulist.header[key])
                        else:
                            val=str(hdulist.header[key])
                    except:
                        val='---' #"\033[31m---\033[39m"
			
		    logdict[key].append(val)

	logtable = Table(logdict.values(), names=logdict.keys()) 
	logtable.write(outlogfile, format='ascii.fixed_width')          

    @staticmethod
    def GenConf(inst='default',ofile='setup.conf'):
        try:
            p=PmPyeasy(PmPyeasy.ConfDir+inst+'/setup.conf').param
        except PmPyeasy.PmPyeasyError as zz:
            print "Invalid instrument name"
            return
        iniout = PmPyeasy.ConfigParser.RawConfigParser(allow_no_value=True)
        iniout.optionxform = str
        iniout.add_section('setup')
        iniout.set('setup','; '+PmPyeasy.__name__+' template file for reduction setup.')
        iniout.set('setup','; Generated by '+PmPyeasy.__name__+' (lib) v'+PmPyeasy.__ver__)
        iniout.set('setup','; with '+inst+' Module v'+str(p['modulever']))
        iniout.set('setup','Instrument',inst)
        for req_par in [k.strip() for k in p['_required_user_pars'].split(',') if k.strip()!='']:
            iniout.set('setup', req_par, str(p[req_par.lower()])+'   ;Please ensure this parameter is correctly set')
        iniout.set('setup','')
        if p['biascor']:
            iniout.add_section('bias')
            iniout.set('bias','')
        if p['flatcor']:
            iniout.add_section('flat')
            iniout.set('flat','')
        iniout.add_section('arc')
        iniout.set('arc','')
        iniout.add_section('std')
        iniout.set('std','; Use "pmpyeasy -liststds" to get list of standard names available by default')
        iniout.set('std','Name','')
        iniout.set('std','')
        iniout.add_section('object')
        
        if not PmPyeasy.os.path.isfile(ofile):
            with open(ofile, 'wb') as configfile:
                iniout.write(configfile)
            print "Created file: "+ofile
        else:
            print "Config file already exist"
        
    @staticmethod
    def RunSNID(sfile,SnidExe):
        from pyraf import iraf
        fname=PmPyeasy.os.path.basename(sfile).rsplit('.', 1)[0]
        iraf.wspectext(sfile+'[*,1]' ,PmPyeasy.os.path.join(PmPyeasy.os.path.dirname(sfile),fname+'.txt'), header='no')
        cwd=PmPyeasy.os.getcwd()
        PmPyeasy.os.chdir( PmPyeasy.os.path.dirname(sfile) )
        try:
            PmPyeasy.os.system(SnidExe+" '"+fname+'.txt'+"'")
        except Exception as e:
            echo ("Failed to run SNID")
        finally:
            PmPyeasy.os.chdir( cwd )
    
    @staticmethod
    def DeleteUnusedFiles(inifile='setup.conf',datadir=''):
        SS=PmPyeasy(inifile,DataDir=datadir)
        safekeep=[]        
        for sec in SS.iniconf.sections():
            if sec.lower() not in ['setup']:
                safekeep=safekeep+SS.parseGrp(sec)
        print("Warning!!!! Following files will be deleted")
        print("Files:")
        lis=[]
        for l in SS.glob.glob(SS.os.path.join(datadir,'*.fits')):
            if l not in safekeep:
                print (l+' -- '+SS.iraf.hselect(l+SS.param['oblg.fitindex'],SS.param['oblg.keys'],'yes',Stdout=1)[0])
                lis=lis+[l]
        if raw_input("Are you really sure to delete these "+str(len(lis)) +" files? [y/n]: ").lower()=='y':
            if raw_input("By what percent you are sure?: ")=='100%':
                for l in lis:
                    try:
                        SS.os.remove(l)
                    except:
                        print ("Failed to delete "+l)
                print("Deleted")
            else:
                print 'Makeup your mind and come back when you are "100%" sure :)'
        else:
            print 'Nothing deleted.'
    
    @staticmethod
    def InstallModule(modpackage):
        import tempfile, zipfile,shutil
        mdir=PmPyeasy.ConfDir

        if not PmPyeasy.os.path.exists( mdir ):
            PmPyeasy.os.mkdir(mdir)

        if modpackage.lower()[:8] == 'https://' or modpackage.lower()[:7] == 'http://':
            import requests
            from StringIO import StringIO
            print "URL supplied, Downloading module package"
            modpackage=StringIO(requests.get(modpackage).content)
        elif modpackage.lower()[-6:] == '.ssmod':
            print "Local package file specified"
        else:
            repo='http://astro.subhashbose.com/pmpyeasy/modules/'
            print "Assuming "+PmPyeasy.__name__+" module name specified"
            print "getting module info and package file from the repository "+repo
            import re, requests
            from StringIO import StringIO
            index=requests.get(repo+'index.ssmod').content
            b = re.findall('<a href="(.*)">(.*)</a></td><td>(.*)</td><td>(.*)</td>', index)
            index=dict([[i[1].strip().lower(), [i[0].strip(), i[2], i[3]]] for i in b])
            try:
                modpackage=repo+index[modpackage.lower()][0]
            except Exception as e:
                raise PmPyeasy.PmPyeasyError("Specified package name not found in repository!! Ensure the package name match exactly listed in "+repo)
            modpackage = StringIO(requests.get(modpackage).content)

        tdir=tempfile.mkdtemp(prefix=PmPyeasy.__name__+'_')

        zf = zipfile.ZipFile(modpackage, 'r')
        zf.extractall(tdir)
        zf.close()

        if not PmPyeasy.os.access(mdir, PmPyeasy.os.W_OK):
            raise PmPyeasy.PmPyeasyError("Insufficient file permission!! You can not install/upgrade Modules, try as superuser")

        for mod in sorted(next(PmPyeasy.os.walk(tdir))[1]):
            target_fullpath=PmPyeasy.os.path.join(mdir,mod)
            if PmPyeasy.os.path.exists( target_fullpath ):
                if raw_input(mod+": module already exist. Do you want to replace with new one? [y/N] : ").lower()=='y':
                    shutil.rmtree(target_fullpath)
                elif raw_input("Do you want to rename the existing module to a new name, and then proceed with installation? [y/N] : ").lower()=='y':
                    new = raw_input("Enter the new name (old is "+mod+") : ")
                    shutil.move(target_fullpath, PmPyeasy.os.path.join(mdir,new))
            if not PmPyeasy.os.path.exists( target_fullpath ):
                shutil.move(PmPyeasy.os.path.join(tdir,mod),  target_fullpath)
                print ('Successfully installed module '+mod+'!!')
            else:
                print ('Skipped installation of module '+mod+'!')
        shutil.rmtree(tdir)

    @staticmethod
    def PackModule(moduledir,modulefile=''):
        if not PmPyeasy.os.path.isdir(moduledir):
            raise PmPyeasy.PmPyeasyError('Module packing error! %s is not a directory' % moduledir)
        if modulefile=='' or PmPyeasy.os.path.isdir(modulefile):
            iniconf = PmPyeasy.ConfigParser.RawConfigParser(allow_no_value=True)
            if len(iniconf.read(PmPyeasy.os.path.join(moduledir,'setup.conf')))==0:
                raise PmPyeasy.PmPyeasyError('Module packing error! %s is not a valid module' % moduledir)
            try:
                ver=str(iniconf.get('setup', 'modulever'))
            except:
                ver='1.0'
            modulefile=PmPyeasy.os.path.join(modulefile,PmPyeasy.os.path.basename(PmPyeasy.os.path.join(moduledir,'')[0:-1])+'-v'+ver+'.ssmod')
        import zipfile
        zf=zipfile.ZipFile(modulefile, 'w', zipfile.ZIP_DEFLATED)
        for root, dirs, files in PmPyeasy.os.walk(moduledir):
            for file in files:
                zf.write(PmPyeasy.os.path.join(root, file),
                         arcname=PmPyeasy.os.path.relpath(PmPyeasy.os.path.join(root, file) ,PmPyeasy.os.path.join(moduledir, '..')) )
        zf.close()
        if PmPyeasy.os.path.isfile(modulefile):
            print "Created: "+modulefile
            return modulefile
        else:
            print "Failed creating: " + modulefile
            return None


        
