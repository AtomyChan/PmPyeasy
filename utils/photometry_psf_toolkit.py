#! /usr/bin/env python

def prepare_pm(pm_in,pm_out,input_image,output_obj,fwhm,skyvalue,):
    
    fid_in = open(pm_in,'r')
    fid_out = open(pm_out,'wt')
    for line in fid_in.readlines():
        if line.split()[0]== 'FWHM':
            if fwhm == -99.99:
                fwhm = 5.0
            line = 'FWHM = ' + str(fwhm) + '    Approx FWHM of objects (pixels) along major axis.\n'
        
        if line.split()[0]== 'SKY':            
            line = 'SKY = ' + str(skyvalue) +   '   Approximate mean sky value in data numbers.\n'
            
        if line.split()[0]== 'IMAGE_IN':            
            line = 'IMAGE_IN = ' + input_image + '  Input image name.\n'
        
        if line.split()[0]== 'OBJECTS_OUT':            
            line = 'OBJECTS_OUT = ' + output_obj +' \n'
            
        fid_out.write(line)
        
    fid_in.close()
    fid_out.close()
    
    
def prepare_pm2(pm_out,input_image,output_obj,fwhm,skyvalue):
    fid_out = open(pm_out,'wt')
    fid_out.write('SKY = ' + str(skyvalue) +'\n')
    fid_out.write('FWHM = '+ str(fwhm) +'\n')
    fid_out.write('EPERDN = 1.4\n')
    fid_out.write('RDNOISE = 13.5\n')
    fid_out.write('TOP = 48000.0\n')
    fid_out.write('IMAGE_IN = '+input_image+'\n')
    fid_out.write('OBJECTS_OUT = '+output_obj+'\n')
    fid_out.write('OBJTYPE_OUT = COMPLETE \n')
    fid_out.write('PARAMS_DEFAULT = paramdefault_fortran \n')
    fid_out.close()