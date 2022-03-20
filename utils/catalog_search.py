#! /usr/bin/python

#Here I put some catalog search function for convenient catalog search
#First created by Ping Chen on 2016-11-11

#supported methods:
#1. local APASS
#2. IRSA catalogs (API and wget)
#3. General Catalog Access from MAST (API and wget)
#3. astroquery.vizier
#4. astroquery.mast


#asrtoquery.mast.Catalogs The Catalogs class provides access to a subset of the astronomical catalogs stored at MAST. The catalogs currently available through this interface are:
#The Hubble Source Catalog (HSC)
#The GALEX Catalog (V2 and V3)
#The Gaia (DR1 and DR2) and TGAS Catalogs
#The TESS Input Catalog (TIC)
#The TESS Candidate Target List (CTL)
#The Disk Detective Catalog
#PanSTARRS (DR1, DR2)

import os
import numpy as np
import wget
from collections import OrderedDict
from astropy.table import Table

from astropy.table import Table, vstack, Column
from astroquery.vizier import Vizier
import astropy.units as u
import astropy.coordinates as coord
from astroquery.mast import Catalogs

def query_local_APASS(ra,dec,distance,output_file='/home/asassn/ASASSN/photometry_pipeline/txtfiles/apass_example_out.tbl', local_database_dir='/home/asassn/SNeData/APASS/', verbose=1):
	'''
	Simle search of local APASS database; this will return entries within square region defined by ra,dec and distance

	INPUTS:
	    	ra: the right ascension of the center of search region in degree unit
	    	dec: the declination of the center of search region in degree unit
	    	distance: the half width of the square in degree unit
		output_file: the filename of output catalog
		local_database_dir: the path to the local APASS data

	PRODUCTS:
		The search results will be store in the output file
	'''
	#define the search squrare
	ramin = ra - distance
	ramax = ra + distance
	decmin = dec - distance
	decmax = dec + distance

	if verbose:
		print "RA search range: [%s, %s]"%(ramin,ramax)
		print "Dec search range: [%s, %s]"%(decmin, decmax)

	fileindex = __APASS_get_region_index(decmin, decmax) #first get the search region index
	print "fileindex involved in the search:", fileindex
	filenames = __APASS_get_region_file(local_database_dir,fileindex) #then get the related filename of APASS data
	print "filenames involved in the search:", filenames

	if os.path.isfile(output_file):
		os.remove(output_file)

	raw = []
	for key in filenames.keys():
		std_ref =  os.path.join(local_database_dir,filenames[key])
		fid = open(std_ref)
		line= fid.readline()
		line= fid.readline()
		while line:
			records = line.split()
			try:
				ra = float(records[1])
				dec = float(records[3])
			except:
				line = fid.readline()
				continue
			if ra>ramin and ra<ramax and dec>decmin and dec<decmax:
				#print ra,dec
				raw.append(records)
			line  = fid.readline()
	fid.close()

	colnames = ['fieldname','RAJ2000','e_RAJ2000','DEJ2000','e_DEJ2000','nobs','mobs','Vmag','B-V','Bmag','g_mag','r_mag','i_mag','e_Vmag','e_B-V','e_Bmag','e_g_mag','e_r_mag','e_i_mag']
	want = Table(np.array(raw), names=colnames)

	if output_file is not None:
		want.write(output_file, format='ascii.csv')

	return want


def __APASS_get_region_index(decl, dech):
	'''
	get the corresponding file index given dec range
	'''
	segments = np.arange(-90,90,5)
	index = segments[(int(np.floor(decl/5.0))+18):(int(np.floor(dech/5.0))+18+1)]
	return index


def __APASS_get_region_file(stdref_database_dir,fileindex):
	'''
	prepare
	'''
	std_ref_files = os.listdir(stdref_database_dir)
	filenames = {}

	for blg in fileindex:
		if blg>0:
			blg_str = str(blg)
			if blg_str == '5':
				blg_str = '05'
			prefix = 'zp'+ blg_str
		elif blg == 0:
			prefix = 'zp00'
		else:
			blg_str = str(blg)
			if str(blg) == '-5':
				blg_str = '-05'
			prefix = 'z'+blg_str.replace('-','m')
		filename = [ref for ref in std_ref_files if prefix in ref]
		filenames[prefix] = filename[0]
	return filenames


def query_local_ATLAS_Refcat2(ra,dec,distance, local_database_master_dir='/home/asassn/SNeData/ATLAS_Refcat2/', outfile=None, verbose=0):
	'''
	Simle search of local ATLAS Refcat2 database; this will return entries within square region defined by ra,dec and distance

	INPUTS:
	    	ra: the right ascension of the center of search region in degree unit
	    	dec: the declination of the center of search region in degree unit
	    	distance: the half width of the square in degree unit
		local_database_master_dir: the path to the local ATLAS-refcat2 data

	PRODUCTS:
		astropy.table
	'''
	#define the search squrare
	ramin = ra - distance
	ramax = ra + distance
	decmin = dec - distance
	decmax = dec + distance
	filenames = __Refcat2_get_region_filenames(ramin, ramax, decmin, decmax) #first search files within the relative region
	print "files involved in the search:", filenames
	units_converter = __Refcat2_prepare_colnames_and_unit_converter()
	colnames = units_converter.keys()
	if not os.path.exists(local_database_master_dir):
		raise ValueError('The expected directory %s containing the database not exists'%local_database_master_dir)

	raw = Table()
	subdirs = ['00_m_16', '16_m_17', '17_m_18', '18_m_19']
	for subdir in subdirs:
		datadir = os.path.join(local_database_master_dir, subdir)
		if not os.path.exists(datadir):
			print "The current subdir %s not available yet"%subdir
			continue

		for i,filename in enumerate(filenames):
			datafile = os.path.join(datadir, filename)
			if not os.path.exists(datafile):
				print 'expected file %s under %s not found'%(filename, subdir)
				continue
			temp = Table.read(datafile, format='ascii')
			for col in ['col25', 'col29', 'col33', 'col37']: #solve the contributor indicator
				temp.replace_column(col, Column(map(lambda x: bin(int(str(x), 16))[2:].zfill(8), temp[col].data), name=col))
			raw = vstack((raw, temp))

	for colnew, colold in zip(units_converter.keys(), raw.colnames):
		convertfactor = units_converter[colnew]
		if convertfactor is not None:
			decimals = int(np.log10(1.0/convertfactor))
			raw[colold] = np.around(raw[colold]*convertfactor, decimals=decimals)

		raw.rename_column(colold, colnew)
	if verbose:
		print raw
	want = raw[(raw['RAJ2000']>ramin)*(raw['RAJ2000']<ramax)*(raw['DEJ2000']>decmin)*(raw['DEJ2000']<decmax)]
	if outfile is not None:
		want.write(outfile, format='ascii.csv')

	return want

def __Refcat2_get_region_filenames(ramin,ramax,decmin, decmax):
	'''
	Get RA, Dec index for files
	'''
	RA_indexs  = [str(raindex).rjust(3,'0') for raindex in np.arange(int(np.floor(ramin)), int(np.floor(ramax))+1)]
	Dec_indexs = [str(decindex) for decindex in np.arange(int(np.floor(decmin)), int(np.floor(decmax))+1)]
	for i,decindex in enumerate(Dec_indexs):
		if decindex[0] == '-':
			decindex = decindex[1:]
			sign = '-'
		else:
			sign = '+'
		decindex = sign+decindex.rjust(2,'0')
		Dec_indexs[i] = decindex
	filenames = []
	for raindex in RA_indexs:
		for decindex in Dec_indexs:
			filenames.append(raindex+decindex+'.rc2')

	return filenames



def __Refcat2_prepare_colnames_and_unit_converter():
	'''
	See Table 4 of https://ui.adsabs.harvard.edu/abs/2018ApJ...867..105T/abstract for the detailed information on column names and units
	'''
	units_converter = OrderedDict()
	units_converter['RAJ2000'] = 1e-8  #column 1
	units_converter['DEJ2000'] = 1e-8
	units_converter['plx'] = 1e-2 # parallax from Gaia DR2; convert to milliarsecond
	units_converter['e_plx'] = 1e-2 #convert to milliarsecond
	units_converter['pmra'] = 1e-2	#proper motion from Gaia DR2
	units_converter['e_pmra'] = 1e-2
	units_converter['pmdec'] = 1e-2
	units_converter['e_pmdec'] = 1e-2
	units_converter['Gmag'] = 1e-3  #Gaia DR2 G magnitude
	units_converter['e_Gmag'] = 1e-3    #column 10
	units_converter['BPmag'] = 1e-3 #Gaia G_{BP} magnitude
	units_converter['e_BPmag'] = 1e-3
	units_converter['RPmag'] = 1e-3
	units_converter['e_RPmag'] = 1e-3
	units_converter['Teff'] = None  #Gaia stellar effective temperature
	units_converter['AGaia'] = 1e-3 #Gaia estimate of G-band extinction
	units_converter['dupvar'] = None #Gaia flags as constant(0), variable(1), or NOT_AVAILABLE(2)+4*DUPLICATE
	units_converter['Ag'] = 1e-3    #SFD estimate of g-band extinction
	units_converter['rp1'] = 0.1    #radius where cumulative G flux exceeds 0.1x this star
	units_converter['r1'] = 0.1     #column 20; radius where cumulative G flux exceeds 1x this star
	units_converter['r10'] = 0.1    #radius where cumulative G flux exceeds 10x this star
	units_converter['gmag'] = 1e-3  #Pan-STARRS g_{P1} magnitude
	units_converter['e_gmag'] = 1e-3#Pan-STARRS g_{P1} magnitude uncertainty
	units_converter['gchi'] = 1e-2  #chi2/dof for contributor to g
	units_converter['gcontrib'] = None #bitmap of contributing catalog to g
	units_converter['rmag'] = 1e-3  #Pan-STARRS r_{P1} magnitude
	units_converter['e_rmag'] = 1e-3#Pan-STARRS r_{P1} magnitude uncertainty
	units_converter['rchi'] = 1e-2  #chi2/dof for contributor to r
	units_converter['rcontrib'] = None #bitmap of contributing catalog to r
	units_converter['imag'] = 1e-3  #Column 30; Pan-STARRS i_{P1} magnitude
	units_converter['e_imag'] = 1e-3#Pan-STARRS i_{P1} magnitude uncertainty
	units_converter['ichi'] = 1e-2  #chi2/dof for contributor to i
	units_converter['icontrib'] = None #bitmap of contributing catalog to i
	units_converter['zmag'] = 1e-3  #Pan-STARRS z_{P1} magnitude
	units_converter['e_zmag'] = 1e-3#Pan-STARRS z_{P1} magnitude uncertainty
	units_converter['zchi'] = 1e-2  #chi2/dof for contributor to z
	units_converter['zcontrib'] = None #bitmap of contributing catalog to z
	units_converter['nstat'] = None  #count of griz deweighted outliers
	units_converter['Jmag'] = 1e-3  #2mass J magnitude
	units_converter['e_Jmag'] = 1e-3#Column 40; 2mass J magnitude uncertainty
	units_converter['Hmag'] = 1e-3  #2mass H magnitude
	units_converter['e_Hmag'] = 1e-3#2mass H magnitude uncertainty
	units_converter['Kmag'] = 1e-3  #2mass K magnitude
	units_converter['e_Kmag'] = 1e-3#2mass K magnitude uncertainty

	return units_converter

def query_VO_SCS(RA,Dec,SR,table='fp_psc',out_format='csv', output_file = '/home/asassn/ASASSN/photometry_pipeline/txtfiles/voscs_example_out.tbl'):
	'''
	IRSA catalogs can be searched via the VO Simple Cone Search.
	Please refer to the following link for details:
	https://irsa.ipac.caltech.edu/docs/vo_scs.html

	The SCS API will return all of the available columns.

	The base URL for the Simple Cone Search is:
	http://irsa.ipac.caltech.edu/SCS?

	INPUTS:
		paramters	values					default description
		table:		fp_psc, iraspsc, etc.			fp_psc	IRSA catalog string identifier "catname" (see below).
		RA:		0.0-360.0				NA	Right Ascension in deg (ICRS, but see Note below).
		Dec:		-90.0-90.0				NA	Declination in deg (ICRS, but see Note below).
		SR:		0.0 < SR < 0.75				NA	Cone search radius in deg.
		format:		votable, ipac_table, csv, tsv, fits	csv	Format of the output table (optional).

		Some Popular Catalogs
		catname for "table" parameter	Description
		wise_allwise_p3as_psd		AllWISE Source Catalog
		fp_psc				2MASS Point Source Catalog
		glimpse_s07			GLIMPSE I Spring 07 Catalog (Spitzer)
		cosmos_phot			COSMOS Photometry Catalog
		iraspsc				IRAS Point Source Catalog

	UUTPUT:


	'''

	url_template = "http://irsa.ipac.caltech.edu/SCS?table=which_table&RA=RA_value_degree&DEC=DEC_value_degree&SR=SR_value_degree&format=out_format"
	url = url_template.replace('which_table',table)
	url = url.replace('RA_value_degree', str(RA))
	url = url.replace('DEC_value_degree', str(Dec))
	url = url.replace('SR_value_degree', str(SR))
	url = url.replace('out_format',out_format)

	if os.path.isfile(output_file):
		os.remove(output_file)

	wget.download(url, out=output_file)


def generic_query_Vizier(ra_deg, dec_deg, radius_deg, catalog, outfile=None, allrows=True, allcolumns=False):
	'''
	Generic query through Vizier with catalog name provided.
	INPUTS:
		ra_deg: RA in degrees
		dec_deg: Declination in degrees
		radius_deg: field radius in degrees

	Catalogs:
	WISE All-Sky Data Release (Cutri+ 2012): II/311/wise
	VISTA Hemisphere Survey band-merged multi-waveband catalogues (VHS) DR4.1: II/359
	VISTA Variables in the Via Lactea Survey (VVV). Catalog Data Release 2 (VVV CAT): II/348
	The VLT Survey Telescope ATLAS. (2015): II/350
	'''
	if allrows:
		Vizier.ROW_LIMIT = -1   #default 50
	if allcolumns:
		Vizier.columns = ['**'] #default ['*']
	field = coord.SkyCoord(ra=ra_deg, dec=dec_deg, unit=(u.deg, u.deg), frame='fk5')
	result =  Vizier.query_region(field, radius=radius_deg*u.deg, catalog=catalog)
	if len(result) == 0:
		print "No catalog found..."
		data = Table()
	else:
		if len(result) > 1:
			print "More than one catalog found as listed below. Data from the first catalog will be returned."
			print result
		data = result[0]

	return data


def SDSS_DR12_query_Vizier(ra_deg, dec_deg, radius_deg, outfile=None, allrows=True, allcolumns=False):
	'''
	https://cdsarc.unistra.fr/ftp/V/147/ReadMe
	'''
	catalog="V/147"
	data = generic_query_Vizier(ra_deg, dec_deg, radius_deg, catalog, outfile=outfile, allrows=allrows, allcolumns=allcolumns)
	if outfile is None:
		outfile = 'SDSS_DR12_RA%s_Dec%s_table_Vizier.csv'%(str(ra_deg, dec_deg))
	data.write(outfile, format='ascii.csv')

def panstarrs_query_Vizier(ra_deg, dec_deg, radius_deg, outfile=None, allrows=True, allcolumns=False):
	"""
	Query PanSTARRS @ VizieR using astroquery.vizier
	https://cdsarc.unistra.fr/viz-bin/ReadMe/II/349?format=html&tex=true
	"""
	catalog="II/349/ps1"
	if dec_deg < -30:
		raise ValueError("requested region not covered by Pan-STARRS")
	data = generic_query_Vizier(ra_deg, dec_deg, radius_deg, catalog, outfile=outfile, allrows=allrows, allcolumns=allcolumns)
	if outfile is None:
		outfile = 'PanSTARRS_DR1_RA%s_Dec%s_table_Vizier.csv'%(str(ra_deg, dec_deg))
	data.write(outfile, format='ascii.csv')

def apass_query_Vizier(ra_deg, dec_deg, radius_deg, outfile=None, allrows=True, allcolumns=False):
	"""
	Query APASS @ VizieR using astroquery.vizier
	https://cdsarc.unistra.fr/viz-bin/ReadMe/II/336?format=html&tex=true
	"""
	catalog="II/336/apass9"
	data = generic_query_Vizier(ra_deg, dec_deg, radius_deg, catalog, outfile=outfile, allrows=allrows, allcolumns=allcolumns)
	if outfile is None:
		outfile = 'APASS_DR9_RA%s_Dec%s_table_Vizier.csv'%(str(ra_deg, dec_deg))
	data.write(outfile, format='ascii.csv')

def twomass_query_Vizier(ra_deg, dec_deg, radius_deg, outfile=None, allrows=True, allcolumns=False):
	"""
	Query 2MASS @ VizieR using astroquery.vizier
	https://cdsarc.unistra.fr/viz-bin/ReadMe/II/246?format=html&tex=true
	"""
	catalog="II/246/out"
	data = generic_query_Vizier(ra_deg, dec_deg, radius_deg, catalog, outfile=outfile, allrows=allrows, allcolumns=allcolumns)
	if outfile is None:
		outfile = 'TwoMASS_RA%s_Dec%s_table_Vizier.csv'%(str(ra_deg, dec_deg))
	data.write(outfile, format='ascii.csv')


def gaia2_query_Vizier(ra_deg, dec_deg, radius_deg, outfile=None, allrows=True, allcolumns=False):
	"""
	Query Gaia DR2 @ VizieR using astroquery.vizier
	https://cdsarc.unistra.fr/viz-bin/ReadMe/I/345?format=html&tex=true
	"""
	catalog="I/345/gaia2"
	data = generic_query_Vizier(ra_deg, dec_deg, radius_deg, catalog, outfile=outfile, allrows=allrows, allcolumns=allcolumns)
	if outfile is None:
		outfile = 'Gaia2_RA%s_Dec%s_table_Vizier.csv'%(str(ra_deg, dec_deg))
	data.write(outfile, format='ascii.csv')

def UKIDSS_query_Vizier(ra_deg, dec_deg, radius_deg, outfile=None, allrows=True, allcolumns=False):
	'''
	UKIDSS is a set of five surveys:
	Large Area Survey (LAS);
	Galactic Plane Survey (GPS);
	Galactic Clusters Survey (GCS);
	Deep Extragalactic Survey (DXS);
	Ultra Deep Survey (UDS).
	See below for more information:
	http://www.ukidss.org/surveys/surveys.html

	UKIDSS DR9 contains data from LAS, GCS and DXS. (II/319)
	https://cdsarc.unistra.fr/viz-bin/ReadMe/II/319?format=html&tex=true

	UKIDSS DR6 contains data from GPS. (II/316/gps6)
	https://cdsarc.unistra.fr/viz-bin/ReadMe/II/316?format=html&tex=true
	'''
	catalog1 = 'II/319'
	data1 = generic_query_Vizier(ra_deg, dec_deg, radius_deg, catalog1, outfile=outfile, allrows=allrows, allcolumns=allcolumns)

	catalog2 = 'II/316/gps6'
	data2 = generic_query_Vizier(ra_deg, dec_deg, radius_deg, catalog1, outfile=outfile, allrows=allrows, allcolumns=allcolumns)

	if len(data1)>0:
		data = data1
	elif len(data2)>0:
		data = data2
	else:
		raise ValueError("requested region not covered by UKIDSS")

	if outfile is None:
		outfile = 'UKIDSS_RA%s_Dec%s_table_Vizier.csv'%(str(ra_deg, dec_deg))
	data.write(outfile, format='ascii.csv')

def gaia_query_mast_Catalogs(ra_deg, dec_deg, radius_deg, outfile=None, version=2):
	'''
	An optional version parameter allows you to select which version you want, the default is the highest version.
	version=2 equals DR2
	'''
	catalog_data = Catalogs.query_region("%s %s"%(ra_deg, dec_deg), radius=radius_deg, catalog="Gaia", version=version)
	if outfile is None:
		outfile = 'Gaia%s_RA%s_Dec%s_table_MAST.csv'%(version,str(ra_deg), str(dec_deg))
	catalog.write(outfile, format='ascii.csv')

def ps1_query_mast_Catalogs(ra_deg, dec_deg, radius_deg, outfile=None, data_release='dr2', table='mean'):
	'''
	The PanSTARRS Catalog has multiple data releases as well as multiple queryable tables.
	An optional data release parameter allows you to select which data release is desired,
	with the default being the latest version (dr2). The table to query is a required parameter.

	DR 1 Mean, and Stack: 'mean', 'stack'
	DR 2 Mean, Stack, and Detection: 'mean', 'stack', 'detection', 'forced_mean'

	See details below:
	https://catalogs.mast.stsci.edu/docs/panstarrs.html
	'''
	catalog_data = Catalogs.query_region("%s %s"%(ra_deg, dec_deg), radius=radius_deg, catalog="Panstarrs", data_release=data_release, table=table)
	if outfile is None:
		outfile = 'PS1_%s_RA%s_Dec%s_table_MAST.csv'%(data_release,str(ra_deg), str(dec_deg))
	catalog.write(outfile, format='ascii.csv')


def query_General_MAST(RA, Dec, SR, FORMAT=None, catalog=None, filename=None, maxobj=None, magFaintEnd = None, magBrightEnd = None, mindet = None):
	'''
	General Catalog Access : http://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?Parameters...

	Required Parameter List
	1 of the following 3 queries - VO ConeSearch, BoxSearch, IDsearch

	RA=ra(deg) &DEC=dec(deg) &SR=search radius(deg)
	BBOX=raMin(deg),decMin(deg),raMax(deg),decMax(deg)
	ID=catID
	Optional Parameters

	FORMAT= VOTABLE(default) | HTML | KML | CSV | TSV | JSON | TEXT(limited set of catalogs)
	CATALOG=GSC23(default) | GSC11 | GSC12 | USNOB | SDSS | FIRST | 2MASS | IRAS | GALEX | GAIA | TGAS | WISE
	| CAOM_OBSCORE | CAOM_OBSPOINTING | PS1V3OBJECTS | PS1V3DETECTIONS
	FILENAME=outputname (directs output to file)
	MAXOBJ=n (limits number of entries returned by brightest magnitude)
	MAGRANGE=bright,faint (limits number of entries returned by limits)
	MINDET=n (minimum numbr of detections PanSTARRS only)
	e.g.

	http://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?RA=83.633083&DEC=22.0145&SR=.01&FORMAT=html&CAT=PS1V3OBJECTS&MINDET=25&MAXOBJ=5
	http://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?CAT=PS1V3OBJECTS&ID=134410836341049171&FORMAT=HTML
	http://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?ra=10&dec=20&sr=0.05&format=html
	http://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?bbox=10,20,10.5,21&format=csv&catalog=gsc11
	http://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?id=NBQI004317

	More details:
	http://gsss.stsci.edu/Software/WebServices.htm


	NOTE: Currently only Cone search implemented!

	INPUTS:
		RA:  right ascension, degree
		Dec: declination, degree
		SR:  search radius, degree
		FORMAT: default VOTABLE
		catalog: default: GSC23
		filename: output filename
		maxobj: limits number of entries returned by brightest magnitude
		magFaintEnd:
		magBrightEnd:
		mindet: minimum numbr of detections PanSTARRS only
	'''

	#check RA, Dec and SR input
	RA_deg_str = str(RA)
	Dec_deg_str = str(Dec)
	SR_deg_str = str(SR)

	RA_deg_yn  = RA_deg_str.replace('.','',1).isdigit()

	Dec_temp = Dec_deg_str.replace('.','',1)
	if Dec_temp[0] == '-':
		Dec_temp = Dec_temp[1:]

	Dec_deg_yn = Dec_temp.isdigit()
	SR_deg_yn  = SR_deg_str.replace('.','',1).isdigit()

	if (not RA_deg_yn) or (not Dec_deg_yn) or (not SR_deg_yn):
		raise ValueError('RA, Dec and SR must be in the format of degree')


	URL_origin = 'http://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?RA=RA_input&DEC=Dec_input&SR=search_radius'

	url = URL_origin.replace('RA_input', RA_deg_str)
	url = url.replace('Dec_input', Dec_deg_str)
	url = url.replace('search_radius', SR_deg_str)

	if FORMAT is not None:
		url = url + '&FORMAT=%s'%FORMAT

	if catalog is not None:
		url = url + '&CATALOG=%s'%catalog

#	if filename is not None:
#		if os.path.isfile(filename):
#			os.remove(filename)
#		url = url + '&FILENAME=%s'%filename

	if maxobj is not None:
		url = url + '&MAXOBJ=%s'%maxobj

	if magFaintEnd is not None and magBrightEnd is not None:
		url = url + '&MAGRANGE=%s,%s'%(magFaintEnd, magBrightEnd)

	if catalog == 'PS1V3OBJECTS' and mindet is not None:
		url = url + '&MINDET=%s'%mindet


	print "The catalog will be downloaded from the address: %s"%url

	if filename is not None:
		if os.path.isfile(filename):
			os.remove(filename)

	wget.download(url, out=filename)


	outlines = open(filename).readlines()[1:]
	#print outlines

	fid = open(filename,'wt')
	for line in outlines:
		fid.write(line)

	fid.close()



if __name__ == '__main__':
	'''
	use the Simple Cone Search (SCS) API to search and download 2MASS catalog
	'''
	import optparse
	parser = optparse.OptionParser()

	def_RA = ''
        parser.add_option('-a', '--ra','--RA', dest="RA", type="string", default=def_RA,   help = "Target Right Ascension in degree [%s]" % def_RA)

        def_DEC = ''
        parser.add_option('-d', '--dec','--DEC', dest="DEC", type="string", default=def_DEC, help = "Target Declination in degrees [%s]" % def_DEC)

        def_SR = 0.15
        parser.add_option('-r', '--radius', dest="search_radius", type="float", default=def_SR,   help = "cone search radius in degree [%s]" % def_SR)

	def_catalog = 'fp_psc'
	parser.add_option('-c', '--catalog', dest='search_catalog', type='string', default=def_catalog, help="the catalog in interest, avaiable ones: apass, fp_psc, glimpse_s07, cosmos_phot, iraspsc [%s]"%def_catalog)

	def_degree = False
        parser.add_option('-D','--degree', dest="radec_degree",action="store_true", default=def_degree, help="whether the RA and Dec inputs are in degree[%s]"%def_degree)

	options, remainder = parser.parse_args()

	RA_str = options.RA
	Dec_str = options.DEC
	SR = options.search_radius
	catalog = options.search_catalog
	radec_degree = options.radec_degree

	if RA_str == '' or Dec_str == '':
		raise IOError( "RA and Dec are required" )

	if not radec_degree:
		from common import radec_format_transformation
		RA, Dec = radec_format_transformation(RA_str, Dec_str)
	else:
		RA = float(RA_str)
		Dec = float(Dec_str)



	valid_catalogs = ['apass', 'fp_psc', 'glimpse_s07', 'cosmos_phot', 'iraspsc', 'GSC23','GSC11','GSC12','USNOB','SDSS','FIRST','2MASS','IRAS','GALEX','GAIA','TGAS','WISE','CAOM_OBSCORE','CAOM_OBSPOINTING','PS1V3OBJECTS', 'PS1V3DETECTIONS']

	IRSA_catalogs = ['fp_psc', 'glimpse_s07', 'cosmos_phot', 'iraspsc']
	MAST_catalogs = ['GSC23','GSC11','GSC12','USNOB','SDSS','FIRST','2MASS','IRAS','GALEX','GAIA','TGAS','WISE','CAOM_OBSCORE','CAOM_OBSPOINTING','PS1V3OBJECTS', 'PS1V3DETECTIONS']

	if catalog not in valid_catalogs:
		raise ValueError("the required catalog not supported yet")

	if catalog == 'apass':
		query_local_APASS(RA, Dec, SR, )
	elif catalog in IRSA_catalogs:
		query_VO_SCS(RA,Dec, SR, table=catalog)
	else:
		FORMAT = 'csv'
		filename = 'VO_MAST_' + catalog + '_example.csv'
		query_General_MAST(RA, Dec, SR, FORMAT=FORMAT, catalog=catalog, filename=filename)
