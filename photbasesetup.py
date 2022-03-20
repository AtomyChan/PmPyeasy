from __future__ import print_function
import os
from astropy.table import Table
import ConfigParser

class basesetup:
	'''
	set up directories, load telescoep info, load pre-defined targets info
	'''
	def __init_attributes__(self):
		self.upload_data_dir = "" #remote server address:folder
		self.save_data_local_dir = "" #photometry result archival folder

		self.data_repository_master_dir = "/Volumes/MyPassport/images" #for pipeline set up
		self.workplace_master_dir = "/Volumes/MyPassport/photometry"
		self.stdcat_dir = "/mnt/md0/stdcat" 
		#self.data_repository_master_dir = "/home/asassn/ASASSN/repository" #for pipeline set up
		#self.workplace_master_dir = "/home/asassn/ASASSN/"

		self.base_dir     = os.path.dirname(os.path.realpath(__file__))
		#self.stdcat_dir = os.path.join(self.base_dir, 'stdcat')
		self.conffile_dir = os.path.join(self.base_dir, 'conf')
		self.extern_dir = os.path.join(self.base_dir, 'extern')

		#------------------ set below properly; if not provided, then try global system command ------
		self.diapl_dir   = os.path.join(self.extern_dir, 'DIAPL/DIAPL_BIN/') #the executable binary folder of diapl
		self.ccdproc_dir = os.path.join(self.extern_dir, 'CCDProc')
		self.dophot_C_dir= os.path.join(self.extern_dir, 'DoPHOT_C-master')
		self.dophot_fortran_dir = os.path.join(self.extern_dir, 'dophot')
		self.fitsh_dir   = os.path.join(self.extern_dir, 'fitsh-0.9.4/src')
		self.hotpants_dir= os.path.join(self.extern_dir, 'hotpants-master')

		self.first_targets_info = None
		self.first_targets_file = os.path.join(self.conffile_dir, 'targets_info_lookup_table_01.txt')
		self.second_targets_info = None
		self.second_targets_file = os.path.join(self.conffile_dir, 'targets_info_lookup_table_02.txt')
		self.telescopes_info_file = os.path.join(self.conffile_dir,'telescopes.info')
		self.telescopes_available = self.__get_telescopes_names(self.telescopes_info_file)
		self.telescopes_info = {}
		self.registerd_photometry_info_file = os.path.join(self.conffile_dir,'photometry_objects.info')
		self.registerd_photometry_objects = self.__get_registerd_photometry_objects(self.registerd_photometry_info_file)
		self.registerd_photometry_objects_info = {}
		self.phot_ret = {}

	def __init__(self):
		'''
		init the class
		'''
		self.__init_attributes__()
		for telescope in self.telescopes_available:
			self.telescopes_info[telescope] = self.__ConfigParser_get_options_given_section(self.telescopes_info_file,telescope)

		for photometry_object in self.registerd_photometry_objects:
			self.registerd_photometry_objects_info[photometry_object] = self.__ConfigParser_get_options_given_section(self.registerd_photometry_info_file,photometry_object)

		try:
			self.first_targets_info = Table.read(self.first_targets_file, format='ascii.fixed_width')
		except:
			print("Failded to load the first target info table %s"%self.first_targets_file)

		try:
			self.second_targets_info = Table.read(self.second_targets_file, format='ascii.fixed_width')
		except:
			print("Failed to load the second target info table %s"%self.second_targets_file)

	def __get_registerd_photometry_objects(self,info_file):
		objects_names = self.__get_configfile_sections(info_file)
		return objects_names


	def __get_telescopes_names(self,info_file):
		telescopes_codes = self.__get_configfile_sections(info_file)
		return telescopes_codes

	def __get_configfile_sections(self,info_file):
                Config = ConfigParser.ConfigParser()
                Config.read(info_file)
	        sections = Config.sections()

		return sections


	def __ConfigParser_get_options_given_section(self,info_file,section):
		'''
		Output:
			info_file: the input file
			section: the config section name
		'''
		Config = ConfigParser.ConfigParser()
		Config.read(info_file)

		info_dict = {}
		options = Config.options(section)
		for option in options:
		    try:
		    	info_dict[option] = Config.get(section, option)
		    	if info_dict[option] == -1:
		        	print("Skip: %s" % option)
		    except:
		        print("Exception on %s!" % option)
		        info_dict[option] = None
		return info_dict
