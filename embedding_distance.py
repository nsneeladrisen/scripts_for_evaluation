import decimal
import numpy as np
import torch
import os
import pickle
from sklearn.metrics.pairwise import euclidean_distances

'''
This script calculates the embedding distance between the sequences with unknown products to that of known products
The distance are only calculated for binding site residues
'''


pretty_print = lambda x: np.format_float_positional(x, trim="-")

'''
This function creates a dictionary of the binding site residue number for each uniprot ID
'''
def sdp(sdp_table):
	dict_uniprot_bs={}
	with open(sdp_table) as filehandle:
		lines_sdp=filehandle.readlines()
		for line_sdp in lines_sdp:
			line_sdp=line_sdp.split('\t')
			dict_uniprot_bs[line_sdp[0]]=line_sdp[1][:-1].split(',')
	return dict_uniprot_bs

'''
This function creates an array for the embeddings for sequences with known products
'''
def embedding_known_product_sequences(embedding_known_dir):
	embedding_known=[]#array of all the embeddings of sequences with known products
	for file in os.listdir(embedding_known_dir):
		if file.endswith('.pt'):
			embedding_known.append(torch.load(embedding_known_dir+file))
	return embedding_known

'''
This function calculates the embedding distance between the corresponding binding site residues of sequences with known products to that of unknown products
'''
def embedding_distance(embedding_unknown,embedding_known,sdp_unknown,sdp_known):
	dict_embedding={}
	for j in range(0,len(embedding_known)):
		sdp_j_embedding=0.0
		sdp_i_embedding=0.0
		count=0
		sum_sdp_embedding=0.0
		for value in range(0,len(sdp_unknown[embedding_unknown['label']])):
			sdp1_known=sdp_known[embedding_known[j]['label']][value]
			sdp1_unknown=sdp_unknown[embedding_unknown['label']][value]
			if sdp1_known !='NA' and sdp1_unknown != 'NA':
				#print(sdp1_known,sdp1_unknown)
				sdp_j_embedding=embedding_known[j]['representations'][33][int(sdp1_known)-1]
				sdp_i_embedding=embedding_unknown['representations'][33][int(sdp1_unknown)-1]
				#print(sdp_j_embedding,sdp_i_embedding)
				sum_sdp_embedding+=(float(euclidean_distances([np.array(sdp_i_embedding)],[np.array(sdp_j_embedding)])))
				count+=1
		dict_embedding[embedding_known[j]['label']]=float(sum_sdp_embedding/count)
	return dict_embedding

'''
This function calls the other functions and sends individual embedding files for sequences with unknown product to calculate embedding distances
'''
def calculate_distance():
	embedding_known=embedding_known_product_sequences('known_product_embedding/')
	sdp_known=sdp('sdp_table_known.tsv')
	sdp_unknown=sdp('sdp_table_unknown.tsv')
	embedding_dict={}
	for file in os.listdir('unknown_seq_embedding/'):
		embedding_unknown=torch.load('unknown_seq_embedding/'+file)
		embedding_dict[embedding_unknown['label']]=embedding_distance(embedding_unknown,embedding_known,sdp_unknown,sdp_known)
		print(file)
	with open('eucledian_distance.pickle', 'wb') as filehandle:
		pickle.dump(embedding_dict,filehandle)

if __name__=="__main__":
	calculate_distance()
