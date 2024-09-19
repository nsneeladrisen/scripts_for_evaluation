import pickle
import os
import sys

'''
This script identifies the interface residues based on pairwise distance between atoms
Each residue is labelled 1/0 as interface/non-interface
PLEASE PROVIDE THE DIRECTORY CONTAINING THE FILES AS INPUT
'''

'''
this function reads files
'''
def read_file(files_dir,file):
	with open(files_dir+file,'r') as filehandler:
		lines=filehandler.readlines()
	return lines

'''
this function returns a list of interface residues for each chain
'''

def identify_interface_residue(files_dir,file):
	interface_res={}
	lines_dict=read_file(files_dir,file)
	for line_dict in lines_dict:
		l=line_dict.split()
		part1,part2=l[0].split(':'),l[1].split(':')
		chain1,chain2=part1[3],part2[3]
		res1,res2=part1[2],part2[2]
		if chain1 not in interface_res:
			interface_res[chain1]=set()
		if chain2 not in interface_res:
			interface_res[chain2]=set()
		if chain1!=chain2:
			interface_res[chain1].add(res1)
			interface_res[chain2].add(res2)
	return interface_res

'''
this function creates a dictionary containing if a residue is interface or not
'''

def interface_noninterface_label(files_dir,file,interface_res):
	i_ni_label={}
	lines_pdb=read_file(files_dir,file)
	for line_pdb in lines_pdb:
		if line_pdb[0:4]=='ATOM':
			chain=line_pdb[21]
			res_no=line_pdb[22:27].strip()
			if chain not in i_ni_label:
				i_ni_label[chain]={}
			if res_no in interface_res[chain] and res_no not in i_ni_label[chain]:
				i_ni_label[chain][res_no]=1
			elif res_no not in interface_res[chain] and res_no not in i_ni_label[chain]:
				i_ni_label[chain][res_no]=0
	return i_ni_label
			
'''
this function calls other functions to create a final dictionary of all the pdb files
'''
def labels(files_dir):
	labels_dict={}
	if '/' not in files_dir:
		files_dir+='/'
	try:
		os.listdir(files_dir)
	except:
		print('Input directory not available')
		sys.exit(1)
	for file in os.listdir(files_dir):
		if file.endswith('.dist'):
			print(file)
			interface_res=identify_interface_residue(files_dir,file)
			int_nonint_label=interface_noninterface_label(files_dir,file[:-5],interface_res)
			labels_dict[file[:-5]]=int_nonint_label
	with open(files_dir+'interface_non_interface_label.pickle','wb') as filehandler:
		pickle.dump(labels_dict,filehandler)
	return labels_dict

if __name__=="__main__":
	labels(sys.argv[1])
