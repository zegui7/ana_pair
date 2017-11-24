#!/usr/bin/env python

"""
ana_pair.py - a method to analyse the interface of related protein complexes.

usage:
example:

ana_pair.py can use as input a family of protein complexes, a MD trajectory or
a PDB ensemble and as the backend of the AnaPair webserver, implemented in shiny
at ____________.

Apart from itself, it requires the "scripts" folder, which includes:
	#R - "fluctuations.R" to calculate fluctuation values using normal mode
	analysis
	#VMD - "get_salt_bridges.tcl" to retrieve salt bridges
	#VMD - "get_sasa_hbonds.tcl" to retrieve SASA and H-bonds
	#VMD - "frame2pdb.tcl" to get individual PDB structures from a MD trajectory

"""

__author__ = "J.G. Almeida"
__email__ = "josegcpa@ebi.ac.uk"

USAGE = __doc__.format(__author__, __email__)

import time
import os
import shutil
import argparse
import sys
from glob import glob,iglob
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
import toolz
import requests
from bs4 import BeautifulSoup
import threading
import queue

#Disable warnings. Can't remember if this works tbh
requests.packages.urllib3.disable_warnings()

#To make sure the script is using the correct file/folder separator
if 'win' in sys.platform:
	SYS_SEP = '\\'
else:
	SYS_SEP = '/'

#Folder locations for webdrivers. Chrome driver works best for debugging
DRIVER_PATH = r'D:\MEGA_Backup\Google Drive\MastersThesis\Z_OTHER\DRIVERS'
CHROME_DRIVER = DRIVER_PATH + SYS_SEP + 'chromedriver.exe'
PHANTOM_DRIVER = DRIVER_PATH + SYS_SEP + 'phantomjs.exe'

#PATH location for external programs
VMD_PATH = 'vmd'
R_PATH = 'Rscript'
CLUSTALO_PATH = r'clustal-omega-1.2.2-win64\clustalo.exe'
SCRIPT_PATH = r'scripts'

#URLs required along the script
COCOMAPS_URL = 'https://www.molnac.unisa.it/BioTools/cocomaps/'
INTERPROSURF_URL = 'http://curie.utmb.edu/usercomplex.html'
CONSURF_URL = 'http://consurf.tau.ac.il/cgi-bin/consurf_2016.cgi'
CONSURF_RESULT_URL = 'http://consurf.tau.ac.il/results/'

RES_CODES = [('CYS', 'C'), ('ASP', 'D'), ('SER', 'S'), ('GLN', 'Q'),
             ('LYS', 'K'), ('ILE', 'I'), ('PRO', 'P'), ('THR', 'T'),
             ('PHE', 'F'), ('ASN', 'N'), ('GLY', 'G'), ('HIS', 'H'),
             ('LEU', 'L'), ('ARG', 'R'), ('TRP', 'W'), ('ALA', 'A'),
             ('VAL', 'V'), ('GLU', 'E'), ('TYR', 'Y'), ('MET', 'M'),
             ('MSE', 'M'), ('SOC', 'C'), ('CYX', 'C'), ('HID', 'H'),
			 ('HIE', 'H'), ('HIP', 'H')]

RES_CODES_DICT = dict(RES_CODES)

#Identifier format: ResnameResnoChain
#Types of data: res,res,value res,value monomer,value complex,value

#TEMP CHECKLIST
#COCOMAPS - CHECK
#INTERPROSURF - CHECK
#HB - CHECK
#HB COCOMAPS - CHECK
#SB - CHECK
#SASA - CHECK
#FLUCTUATIONS - CHECK
#CONSURF - CHECK

#MD TRAJECTORY ANALYSIS -

#NMR ENSEMBLE -

#threading-related functions

def thread_maker(func):
	def thread_sub(*k,**kw):
		thread = threading.Thread(target = func, args = k, kwargs = kw)
		thread.daemon = True
		thread.start()
		return thread
	return thread_sub

def get_paths():
	paths = []
	with open('paths','r') as o:
		lines = o.readlines()
		for line in lines:
			data = line.split('=')
			paths.append(data[1])
	CHROME_DRIVER = paths[0]
	PHANTOM_DRIVER = paths[1]
	VMD_PATH = paths[2]
	R_PATH = paths[3]
	CLUSTALO_PATH = paths[4]
	SCRIPT_PATH = paths[5]
	return paths

#Classes

class Complex():

	"""
	This will define the main .pdb class. It will be used to fetch scores from
	CoCoMaps, InterProSurf, HB, SB, SASA and fluctuations, and use outputs from
	Consurf and EVFold.
	threading is being used to run all processes parallelly to speed up the
	whole thing.
	"""

	def __init__(self,pdb,download_path = 'temp_data',
				 inter_path = 'inter',output_path = 'output'):

		self.pdb = pdb
		self.download_path = download_path
		self.inter_path = inter_path
		self.output_path = output_path
		self.identifier = self.get_id()
		self.thread_dict = {}

		self.folder = SYS_SEP.join(self.pdb.split(SYS_SEP)[:-1])
		self.atoms = toolz.pipe(self.pdb,self.open_lines,self.get_atoms)
		self.chains = self.get_chains()
		self.pdb2fasta(True)

		self.thread_dict["CoCoMaps"] = self.cocomaps_wraper()
		self.thread_dict["InterProSurf"] = self.interprosurf_wraper()
		self.thread_dict["Fluctuations"] = self.fluctuation_wraper()
		self.thread_dict["Structural"] = self.vmd_wraper()
		self.thread_waiter()

		print("Data for " + self.identifier + " has been retrieved.")

	def __repr__(self):
		return "<Complex(" + self.identifier + ') object>'

	def get_id(self):
		temp = self.pdb.split(SYS_SEP)[-1]
		identifier = temp.split('.')[0]
		return identifier

	def open_lines(self,file):
		with open(file,'r') as o:
			lines = o.readlines()
			return lines

	def open_write(self,file_name,data):
		with open(file_name,'w') as o:
			o.write(data)

	def write_lines(self,lines,file_name):
		with open(file_name,'w') as o:
			for line in lines:
				data = ','.join([str(x) for x in line]) + '\n'
				o.write(data)

	#threading-related scripts

	def thread_waiter(self):
		ready = False
		while ready == False:
			for thread in self.thread_dict:
				self.thread_dict[thread]
				if not self.thread_dict[thread].isAlive():
					print(self.identifier,"OUTPUT -",thread,
						  "data has been retrieved.")
					self.thread_dict[thread].handled = True
					self.thread_dict[thread].join()
				else:
					self.thread_dict[thread].handled = False
			self.thread_dict = {thread: self.thread_dict[thread]
								for thread in self.thread_dict
								if self.thread_dict[thread].handled == False}
			if len(self.thread_dict) == 0:
				ready = True
			time.sleep(5)

	# PDB-file related scripts

	def get_atoms(self,lines):
		atom_lines = []
		for line in lines:
			if line[0:4] == 'ATOM':
				atom_lines.append(line)
		return atom_lines

	def get_chains(self):
		chains = set()
		for line in self.atoms:
			chain = line[21]
			chains.add(chain)
		chains = list(chains)
		return sorted(chains)

	def pdb2fasta(self,pdbsplit = False):
		fasta_ids = ['_'.join([self.identifier,chain]) for chain in self.chains]
		self.fasta_dict = {identifier: '' for identifier in fasta_ids}
		if pdbsplit == True:
			struc_ids = ['_'.join([self.identifier,chain])
						 for chain in self.chains]
			struc_dict = {identifier: '' for identifier in struc_ids}
		temp = []
		counter = 0
		for atom in self.atoms:
			resname = atom[17:20].strip()
			resno = atom[22:26].strip()
			chain = atom[21].strip()
			identifier = resname + resno + chain
			chain_id = '_'.join([self.identifier,chain])
			if identifier not in temp:
				self.fasta_dict[chain_id] += RES_CODES_DICT[resname]
				temp.append(identifier)
				counter += 1
				if counter == 75:
					self.fasta_dict[chain_id] += '\n'
					counter = 0
			if pdbsplit == True:
				struc_dict[chain_id] += atom +'\n'
		if pdbsplit == True:
			for pdb in struc_dict:
				file_name = self.inter_path + SYS_SEP + pdb + '.pdb'
				with open(file_name,'w') as o:
					o.write(struc_dict[pdb])

		self.fasta_files = []
		for fasta in self.fasta_dict:
			file_name = self.inter_path + SYS_SEP + fasta + '.fasta'
			with open(file_name,'w') as o:
				o.write('>' + fasta + '\n' + self.fasta_dict[fasta])
			self.fasta_files.append(file_name)

	#Webdriver related functions

	def start_driver(self,drive_func):
		if drive_func == 'phantom':
			return webdriver.PhantomJS(PHANTOM_DRIVER)
		elif drive_func == 'chrome':
			return webdriver.Chrome(CHROME_DRIVER)

	#CoCoMaps related functions

	def cocomaps_digger(self):
		params = {'pdb_file':self.pdb,
				  'local_chain_molecule1':self.chains[0],
				  'local_chain_molecule2':self.chains[1]}
		new_names = ["interaction_overview","minimum_distances",
					 "hbonds","asa_all","asa_mol1","asa_mol2"]

		driver = self.start_driver("phantom")
		driver.get(COCOMAPS_URL)
		assert "CoCoMaps" in driver.title

		driver.find_elements_by_name("type_pdb")[1].click()
		for param in params:
			pdb = driver.find_element_by_name(param)
			pdb.send_keys(params[param])
		driver.find_element_by_name("submit").click()
		time.sleep(5)
		download_list = driver.find_elements_by_class_name("downloadTable")[3:]

		for i in range(0,len(download_list)):
			curr_download = download_list[i]
			url_name = curr_download.get_attribute("href")
			new_name = self.download_path + SYS_SEP + 'COCOMAPS_'
			new_name += '_'.join([self.identifier,new_names[i]]) + '.txt'
			r = requests.get(url_name,verify = False)
			with open(new_name, "wb") as table:
				table.write(r.content)

		dist_element = driver.find_element_by_xpath('//a[@target ="_table"]')
		dist_element.click()
		time.sleep(5)
		while len(driver.window_handles) == 1:
			time.sleep(1)
		driver.switch_to_window(driver.window_handles[1])
		download_dist = driver.find_element_by_class_name("downloadTable")
		table_url = download_dist.get_attribute("href")
		r = requests.get(table_url,verify = False)
		table_name = self.download_path + SYS_SEP + 'COCOMAPS_'
		table_name += self.identifier + "_dist_table.txt"
		with open(table_name, "wb") as table:
			table.write(r.content)
		driver.quit()

	def process_cocomaps(self):
		comp_path = self.download_path + SYS_SEP + 'COCOMAPS_'
		comp_path += self.identifier + '_asa_all.txt'
		residue_path = self.download_path + SYS_SEP + 'COCOMAPS_'
		residue_path += self.identifier + '_asa_mol*'
		hbond_path = self.download_path + SYS_SEP + 'COCOMAPS_'
		hbond_path +=self.identifier + '_hbonds.txt'

		complex_data = []
		residue_data = []
		hbond_data = []

		comp_data = self.open_lines(comp_path)
		for line in comp_data:
			temp = line.split()
			value = temp[-1]
			identifier = '_'.join(temp[:-1])
			content = [identifier,value]
			complex_data.append(content)

		residue_glob = glob(residue_path)
		for globed in residue_glob:
			res_data = self.open_lines(globed)
			for line in res_data:
				temp = line.split()
				value = temp[-1]
				identifier = ''.join(temp[0:3])
				content = [identifier,value]
				residue_data.append(content)

		hb_data = self.open_lines(hbond_path)
		for line in hb_data:
			temp = line.split()
			value = temp[-1]
			id1 = ''.join([temp[6],temp[5],temp[4]])
			id2 = ''.join([temp[6],temp[5],temp[4]])
			content = [identifier,value]
			hbond_data.append(content)

		return complex_data,residue_data,hbond_data

	#InterProSurf-related scripts

	def interprosurf_digger(self):
		complex_col = ["chain","polar_area_energy","apolar_area_energy",
					   "total_area_energy","no_surf_atoms","no_buried_atoms"]
		all_data = []
		complex_data = []
		residue_data = []
		chains = set()
		driver = self.start_driver("phantom")
		driver.get(INTERPROSURF_URL)
		assert "InterProSurf" in driver.title

		name = driver.find_element_by_name("name")
		email = driver.find_element_by_name("email")
		file = driver.find_element_by_name("PDBfile")

		name.send_keys(self.identifier)
		email.send_keys('jgalmeida7@hotmail.com')
		file.send_keys(self.pdb)

		driver.find_element_by_xpath("//input[@type='submit']").click()
		time.sleep(5)
		driver.switch_to_window(driver.window_handles[1])
		all_html = driver.page_source
		soup = BeautifulSoup(all_html,"html.parser")
		for mess in soup.find_all('table'):
			bowl = mess.get_text()
			data = bowl.split('\n')
			data = [x.strip() for x in data if x.strip() != '']
			all_data.append(data)
		for data in all_data[4:]:
			if len(data[-2]) == 1:
				chains.add(data[-2])
				residue_data.append(data)
		chains = sorted(list(chains))
		nchains = len(chains)
		for i in range(0,nchains):
			new_data = [chains[i]] + [x.split('=')[-1].strip()
						for x in all_data[i][1:]]
			complex_data.append(new_data)
		complex_full = [''.join(chains)] + [x.split('=')[-1].strip()
						for x in all_data[nchains][1:]]
		complex_data.append(complex_full)
		residue_col = all_data[nchains + 1]
		residue_name = self.download_path + SYS_SEP + 'INTERPROSURF_'
		residue_name += self.identifier + '_residues.csv'
		complex_name = self.download_path + SYS_SEP + 'INTERPROSURF_'
		complex_name += self.identifier + '_complex.csv'
		residue_write = ','.join(residue_col) + '\n'
		residue_write +='\n'.join([','.join(data) for data in residue_data])
		complex_write = ','.join(complex_col) + '\n'
		complex_write += '\n'.join([','.join(data) for data in complex_data])
		self.open_write(residue_name,residue_write)
		self.open_write(complex_name,complex_write)
		driver.quit()

	def process_interprosurf(self):
		comp_path = self.download_path + SYS_SEP + 'INTERPROSURF_'
		comp_path += self.identifier + '_complex.csv'
		res_path = self.download_path + SYS_SEP + 'INTERPROSURF_'
		res_path += self.identifier + '_residues.csv'
		comp_data = self.open_lines(comp_path)
		complex_data = []
		residue_data = []
		header = comp_data[0].split(',')
		monomer_a = comp_data[1].split(',')
		monomer_b = comp_data[2].split(',')
		comp = comp_data[3].split(',')
		for i in range(0,len(header)):
			temp = [header[i],monomer_a[i],monomer_b[i],comp[i]]
			complex_data.append(temp)
		res_data = self.open_lines(res_path)
		for res in res_data[1:]:
			pre_temp = res.split(',')
			identifier = pre_temp[4] + pre_temp[1] + pre_temp[5]
			residue_data.append([identifier,pre_temp[6]])
		return complex_data,residue_data

	#VMD-related scripts

	def vmd_path_normalizer(self):
		self.vmd_path = self.pdb.replace('\\','/')

	def run_vmd(self,script):
		os.system(' '.join([VMD_PATH,'-dispdev none -nt -e',script,'-args',
							self.vmd_path, '> nul']))

	def process_saltbr(self):
		data = []
		path = SYS_SEP.join([self.folder,"saltbr_output",
							 self.identifier,"*.dat"])
		all_saltbr = glob(path)
		for saltbr in all_saltbr:
			x = saltbr[:-4]
			spl = x.replace('chain','')
			spl = spl.split('-')
			saltbr_a = spl[-2].split('_')
			saltbr_b = spl[-1].split('_')
			if saltbr_a[1] != saltbr_b[1]:
				in_data = [''.join(saltbr_a),''.join(saltbr_a)]
				with open(saltbr,'r') as o:
					content = o.read()
					content = in_data + [str(content.split()[1])]
					data.append(content)
		return data

	def process_sasa(self):
		data = []
		path = SYS_SEP.join([self.folder,"sasa_hb",
							 "sasa_" + self.identifier + ".log"])
		rl = self.open_lines(path)
		for line in rl:
			spl = line.split()
			content = [spl[2] + spl[0] + spl[1],spl[3]]
			data.append(content)
		return data

	def process_hbonds(self):
		data = []
		path = SYS_SEP.join([self.folder,"sasa_hb",
							 "hb_" + self.identifier + ".log"])

		def get_id(string):
			temp = string.split('-')
			final = temp[1] + temp[0][-1]
			return final

		rl = self.open_lines(path)
		for line in rl[2:]:
			spl = line.split()
			id1 = get_id(spl[0])
			id2 = get_id(spl[1])
			content = [id1,id2,1]
			data.append(content)
		return data

	#R-related functions

	def run_r(self,script,arguments):
		arguments_ready = ' '.join(arguments)
		print(' '.join([R_PATH,script,arguments_ready,'> nul']))
		os.system(' '.join([R_PATH,script,arguments_ready,'> nul']))


	#Wrapers

	@thread_maker
	def cocomaps_wraper(self):
		self.cocomaps_digger()
		self.cocomaps_comp,self.cocomaps_res,self.cocomaps_hbonds = self.process_cocomaps()

		self.write_lines(self.cocomaps_comp,self.inter_path + SYS_SEP +
						 'cocomaps_comp_' + self.identifier + '.csv')
		self.write_lines(self.cocomaps_res,self.inter_path + SYS_SEP +
						 'cocomaps_res_' + self.identifier + '.csv')
		self.write_lines(self.cocomaps_hbonds,self.inter_path + SYS_SEP +
						 'cocomaps_hbonds_' + self.identifier + '.csv')

	@thread_maker
	def vmd_wraper(self):
		self.vmd_path_normalizer()
		self.run_vmd(SCRIPT_PATH + SYS_SEP + "get_salt_bridges.tcl")
		self.saltbr = self.process_saltbr()

		self.run_vmd(SCRIPT_PATH + SYS_SEP + "get_sasa_hbonds.tcl")
		self.sasa = self.process_sasa()
		self.hbonds = self.process_hbonds()

		self.write_lines(self.saltbr,self.inter_path + SYS_SEP +
						 'vmd_saltbr_' + self.identifier + '.csv')
		self.write_lines(self.sasa,self.inter_path + SYS_SEP +
						 'vmd_sasa_' + self.identifier + '.csv')
		self.write_lines(self.hbonds,self.inter_path + SYS_SEP +
						 'vmd_hbonds_' + self.identifier + '.csv')

	@thread_maker
	def interprosurf_wraper(self):
		self.interprosurf_digger()
		self.interprosurf_comp,self.interprosurf_res = self.process_interprosurf()

		self.write_lines(self.interprosurf_comp,self.inter_path + SYS_SEP +
						 'interprosurf_comp_' + self.identifier + '.csv')
		self.write_lines(self.interprosurf_res,self.inter_path + SYS_SEP +
						 'interprosurf_res_' + self.identifier + '.csv')

	@thread_maker
	def fluctuation_wraper(self):
		self.run_r(SCRIPT_PATH + SYS_SEP + "fluctuations.R",
				   [self.pdb,self.identifier,self.inter_path + SYS_SEP +
				    "fluctuations_" + self.identifier + '.csv'])

class Complex_Consurf():
	#Consurf-related scripts

	def __init__(self,pdb,inter,chain,identifier):
		self.pdb = pdb
		self.inter_path = inter
		self.chain = chain
		self.identifier = identifier
		self.thread_consurf = self.consurf_wraper()
		self.thread_consurf.join()

	def write_lines(self,lines,file_name):
		with open(file_name,'w') as o:
			for line in lines:
				data = ','.join([str(x) for x in line]) + '\n'
				o.write(data)

	#Consurf-related scripts

	def consurf_digger(self,chain):
		consurf_data = {'Run_Number':'NONE',
		'DNA_AA':'AA',
		'Run_Number':'NONE',
		'PDB_yes_no':'yes',
		'MAX_FILE_SIZE':'2000000',
		'MSA_yes_no':'no',
		'FASTA_text':'',
		'MSAprogram':'MAFFT',
		'proteins_DB':'UNIREF90',
		'MAX_NUM_HOMOL':'150',
		'MAX_REDUNDANCY':'95',
		'MIN_IDENTITY':'35',
		'ITERATIONS':'3',
		'E_VALUE':'0.001',
		'Homolog_search_algorithm':'BLAST',
		'ALGORITHM':'Bayes',
		'SUB_MATRIX':'BEST',
		'user_email':'',
		'send_user_mail':'no',
		'submitForm':'Submit',
		'PDB_chain':chain
		}

		consurf_files = {'pdb_FILE': open(self.pdb),
		'PDB_chain': chain,
		'JOB_TITLE': self.identifier + '_' + chain
		}

		r = requests.post(CONSURF_URL, data = consurf_data,
						  files = consurf_files)

		job_id = r.url.split('/')[4]
		result_url = CONSURF_RESULT_URL + job_id + '/consurf.grades'

		while r.status_code != requests.codes.ok:
			time.sleep(1)

		print("Consurf results will be stored in",r.url,'for',self.identifier)
		return result_url

	def consurf_retriever(self,url):
		min_counter = 0
		output = []
		ready = False
		while ready == False:
			r = requests.get(url)
			time.sleep(60)
			if r.status_code == requests.codes.ok:
				if 'Sorry.' in r.text:
					ready = False
				else:
					ready = True
				min_counter += 1
				if min_counter == 240:
					self.consurf_res = 'time_out'
		for line in r.text.split('\n'):
			data = line.split()
			if len(data) > 0 and data[0].isdigit() == True:
				if data[2] != '-':
					entry = [data[2].replace(':',''), data[3]]
					output.append(entry)

		write_lines(output,'consurf_' + self.inter_path + SYS_SEP +
						   self.identifier + '.csv')

	def consurf_wraper(self):
		output = self.consurf_digger(self.chain)
		thread = self.consurf_retriever(output)
		thread.join()

class AllComplexes():

	"""
	This class is used to compile all Complex classes in the same object.
	ana_type defines whether this should be a family analysis, a trajectory
	analysis or an ensemble analysis. The difference between the three are based
	on the user	input:
	- "family" analyses a family of different complexes between related proteins
	i.e. GPCR-partner complexes;
	- "trajectory" analyses a MD trajectory;
	- "ensemble" analyses a PDB file with several different models i.e. NMR
	structures
	"""

	def __init__(self,download_path = 'temp_data',inter_path = 'inter',
				 output_path = 'output',ana_type = "family",**kw):
		self.download_path = download_path
		self.inter_path = inter_path
		self.output_path = output_path
		self.folder = kw["folder"]
		self.ana_type = ana_type

		self.complexes = {}
		self.thread_dict = {}
		self.consurf_submitted = []
		self.submitted = []

		print("Downloads will be stored in " + self.download_path + '.')
		print("Intermediate files will be stored in " + self.inter_path + '.')
		print("Output files will be stored stored in " + self.output_path + '.')

		self.create_folders()
		if ana_type == "family":
			self.consurf_thread_dict = {}
			self.monomer_1_partners = {}
			self.monomer_2_partners = {}
			self.get_pdbs()
			self.get_combinations()
			self.get_complexes_family()
			self.align_wraper()
			while len(self.consurf_thread_dict) > 0:
				self.consurf_thread_dict = self.check_thread_dict(self.consurf_thread_dict)
				time.sleep(60)

		elif ana_type == "trajectory":
			self.coordinates = self.folder + SYS_SEP + kw["coordinates"]
			self.trajectory = self.folder + SYS_SEP + kw["trajectory"]
			self.split_trajectory()
			self.get_complexes_trajectory()

		elif ana_type == "ensemble":
			pass

	def open_write(self,file_name,to_write):
		with open(file_name,'w') as o:
			o.write(to_write)

	def grab_chains(self,file):
		with open(file,'r') as o:
			lines = o.readlines()
			chains = set(line[21] for line in lines if line[0:4] == 'ATOM')
		return chains

	def create_folders(self):
		folders = [self.download_path,self.inter_path,self.output_path]
		for folder in folders:
			if os.path.isdir(folder) == False:
				os.makedirs(folder)

	def get_pdbs(self):
		all_pdbs = glob(self.folder + SYS_SEP + '*.pdb')
		self.all_pdbs = [os.path.abspath(x) for x in all_pdbs]

	def get_combinations(self):
		def add_mon(mon_dict,id1,id2):
			if id1 not in mon_dict:
				mon_dict[id1] = [id2]
			else:
				mon_dict[id1].append(id2)
		for pdb in self.all_pdbs:
			pdb_id = pdb.split(SYS_SEP)[-1][:-4]
			id1,id2 = pdb_id.split('-')
			add_mon(self.monomer_1_partners,id1,id2)
			add_mon(self.monomer_2_partners,id2,id1)

	@thread_maker
	def get_complex(self,file):
		pdb = file.replace('.pdb','')
		pdb = pdb.split(SYS_SEP)[-1]
		sep_ids = pdb.split('-')
		self.complexes[pdb] = [file,(sep_ids[0],sep_ids[1]),Complex(file)]

	@thread_maker
	def get_consurf(self,file,chain,identifier):
		Complex_Consurf(file,self.inter_path,chain,identifier)

	def check_thread_dict(self,thread_dict):
		for thread in thread_dict:
			if not thread_dict[thread].isAlive():
				thread_dict[thread].handled = True
				thread_dict[thread].join()
			else:
				thread_dict[thread].handled = False
		thread_dict = {thread: thread_dict[thread]
					   for thread in thread_dict
					   if thread_dict[thread].handled == False}
		return thread_dict

	def get_complexes_family(self): #wraper for the family method
		#This submits to everything to Consurf
		for pdb in self.all_pdbs:
			id1,id2 = pdb[:-4].split(SYS_SEP)[-1].split('-')
			if id1 not in self.consurf_thread_dict:
				chain = list(self.grab_chains(pdb))[0]
				self.consurf_thread_dict[id1] = self.get_consurf(pdb,chain,id1)
			if id2 not in self.consurf_thread_dict:
				chain = list(self.grab_chains(pdb))[1]
				self.consurf_thread_dict[id2] = self.get_consurf(pdb,chain,id2)

		#This calculates the rest of the features
		while len(self.complexes) != len(self.all_pdbs):
			self.thread_dict = self.check_thread_dict(self.thread_dict)
			for pdb in self.all_pdbs:
				if len(self.thread_dict) < 2 and pdb not in self.submitted:
					self.thread_dict[pdb] = self.get_complex(pdb)
					self.submitted.append(pdb)
			time.sleep(1)

	def combine_fastas(self):
		self.multi_fasta = []
		fasta_a = {}
		fasta_b = {}
		fasta_a_write = ''
		fasta_b_write = ''
		for pdb in self.complexes:
			id1,id2 = pdb.split('-')
			fasta_dict = self.complexes[pdb][2].fasta_dict
			sorted_keys = sorted(fasta_dict.keys())
			if id1 not in fasta_a:
				fasta_a[id1] = ['>',id1,'\n',fasta_dict[sorted_keys[0]],'\n']
				fasta_a[id1] = ''.join(fasta_a[id1])
			if id2 not in fasta_b:
				fasta_b[id2] = ['>',id2,'\n',fasta_dict[sorted_keys[1]],'\n']
				fasta_b[id2] = ''.join(fasta_b[id2])
		for fasta in fasta_a:
			fasta_a_write += fasta_a[fasta]
		for fasta in fasta_b:
			fasta_b_write += fasta_b[fasta]

		self.open_write(self.inter_path + SYS_SEP + 'monomer_1.fasta',
						fasta_a_write)
		self.open_write(self.inter_path + SYS_SEP + 'monomer_2.fasta',
						fasta_b_write)

	def run_clustalo(self,in_file,out_file):
		os.system(' '.join([CLUSTALO_PATH,'--force -i',in_file,'-o',out_file]))

	def align_wraper(self):
		print("Aligning...")
		self.combine_fastas()
		self.run_clustalo(self.inter_path + SYS_SEP + 'monomer_1.fasta',
						  self.inter_path + SYS_SEP + 'monomer_1_aln.fasta')
		self.run_clustalo(self.inter_path + SYS_SEP + 'monomer_2.fasta',
						  self.inter_path + SYS_SEP + 'monomer_2_aln.fasta')
		print("Monomer pairs have been aligned.")

	#trajectory-related scripts

	def split_trajectory(self):
		def run_vmd(script,arguments):
			os.system(' '.join([VMD_PATH,'-dispdev none -nt -e',script,'-args',
								arguments, '> nul']))
		def get_id():
			return self.coordinates[:-4].split(SYS_SEP)[-1]

		arguments = ' '.join([
		self.coordinates,
		self.trajectory,
		get_id(),
		self.inter_path
		])

		run_vmd(SYS_SEP.join([SCRIPT_PATH,"frame2pdb.tcl"]),arguments)

	def get_complexes_trajectory(self):
		def run_r(script,arguments):
			arguments_ready = ' '.join(arguments)
			os.system(' '.join([R_PATH,script,arguments_ready,'> nul']))

		def process_chains():
			with open("chains","r") as o:
				lines = o.readlines()[1:]
			lines = [x.replace('"','') for x in lines]
			lines = [x.split()[1] for x in lines]
			return lines

		def correct_chains(file):
			with open(file,'r') as o:
				lines = o.readlines()
				o.close()
			new_lines = [lines[0]]
			for index,line in enumerate(lines[1:]):
				if index > len(chains) - 1:
					chain_index = len(chains) - 1
				else: chain_index = index
				if chains[chain_index] != 'X':
					line = line[:21] + chains[chain_index] + line[22:]
					new_lines.append(line)
			return ''.join(new_lines)

		self.pdb_path = SYS_SEP.join([os.getcwd(),self.inter_path,
									  'all_frames','*pdb'])
		self.all_pdbs = glob(self.pdb_path)

		#this is essentially to add chains to unformatted files coming from
		#molecular dynamics.

		run_r(SCRIPT_PATH + SYS_SEP + "get_chains.R",
			  [SYS_SEP.join([os.getcwd(),self.inter_path,"all_frames"])])
		chains = process_chains()
		consurf_submitted = False
		while len(self.complexes) != len(self.all_pdbs):
			self.thread_dict = self.check_thread_dict(self.thread_dict)
			for pdb in self.all_pdbs:
				new_file = correct_chains(pdb)
				if new_file[0:4] == 'CRYS':
					with open(pdb,'w') as o: o.write(new_file)
				if consurf_submitted == False:
					id1,id2 = pdb[:-4].split(SYS_SEP)[-1].split('_')[0:2]
					chains = list(self.grab_chains(pdb))
					self.get_consurf(pdb,chains[0],id1)
					self.get_consurf(pdb,chains[1],id2)
					consurf_submitted = True
				if len(self.thread_dict) < 1 and pdb not in self.submitted:
					self.thread_dict[pdb] = self.get_complex(pdb)
					self.submitted.append(pdb)
			time.sleep(1) 

q = queue.Queue()
#AllComplexes(download_path = 'temp_data',inter_path = 'inter', output_path = 'output',
#			 ana_type = "family",folder = "test_family")
AllComplexes(download_path = 'temp_data',inter_path = 'inter',
			 output_path = 'output',ana_type = "trajectory",folder = "test_trajectory",
			 coordinates = "A21-R1.gro",trajectory = "A21-R1.xtc")
