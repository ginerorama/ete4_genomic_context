## ETE Context v0.1 2021 
# Joaquin Giner
# 
# Requirements
# -----------------------
# Python packages [Biopython]
# Binaries [MAFFT,Fasttree]
# -----------------------
#

import sys,os,re
from collections import Counter
from Prosite import *
from subprocess import Popen, PIPE
import argparse
import time

start_time = time.time()







def query_protein_parser(query_protein):
	
	"""
	transform query fasta protein into a list of two element query_protein_list, 
	first contains fasta header and second contain fasta sequence 
	["header","sequence"]
	"-","_","\t" and " " characters are removed
	
	"""

	protein_sequence = [] #list to join protein seq of fasta with "\n"
	query_protein_list = []	
	with open(query_protein,"r") as query_file:
		lines = query_file.readlines()
		for element in lines:
			if element.startswith(">"):
				query_protein_name = element.rstrip().replace(">","").replace(" ","_")
			else: 
				if element != "\n":
					clean_element = element.rstrip().replace("\t","").replace(" ","")
					protein_sequence.append(clean_element)

	query_protein_list = [query_protein_name, "".join(protein_sequence)]


	return query_protein_list



def make_list_and_dict(query_protein_list):

	"""
	Parseamos utilizando biopython un gbk para extraer primero
	una lista (list_of_genes) con todos los locus_tag ordenados
	segun su presencia en el genoma y por otro lado un diccionario
	(dict_of_genes) en el que almacenamos para cada locus_tag su 
	strand y su product (función).
	Fasta_dict contain the proteome retrieved from the gbff in fasta format
	"""
	
	print("parsing gbffs with biopython")
	from Bio import SeqIO

	list_of_genes = [] # list of genes for each genomes ordered as in the gbk file
	dict_of_genes = {} # key:locus_tag ; strand, product, genome_name

	fasta_for_clustering = open("DB.fasta","w") #fasta with all protein sequence used as input for mmseqs2
	genbank_directory = os.getcwd()+"/"+gbff_folder

	dir_list = os.listdir(genbank_directory)
	fasta_dict = {}

	for file in dir_list:
		if file.endswith(".gbff") or file.endswith(".gbk"): 
			with open(genbank_directory+file, "r") as input_handle:
				for record in SeqIO.parse(input_handle, "genbank"):
					name = record.description.split(",")[0].replace("chromosome","") #genome name
					for feature in record.features:
						if feature.type == 'CDS':
							locus_tag = feature.qualifiers['locus_tag'][0] 
							#contains all genes from a gbff ordered
							list_of_genes.append(locus_tag)

							strand = feature.location.strand

							try:
								product = feature.qualifiers['product'][0] 
							except:
								product = "NA"

							try:
								gene = feature.qualifiers['gene'][0] 
							except:
								gene = ""

							dict_of_genes[locus_tag]=[strand,product,name,gene]

							try:
								sequence = feature.qualifiers['translation'][0]
								fasta_for_clustering.write(">"+locus_tag+"\n"+sequence+"\n")
								fasta_dict[locus_tag]=sequence
							except:
								#This try except is due to the fact that in gbff pseudogenes
								#print("was not possible to obtain the "+locus_tag+" sequence")
								continue
	
	locus_tag,sequence = query_protein_list					
	fasta_for_clustering.write(">"+locus_tag+"\n"+sequence+"\n")
	fasta_for_clustering.close()


	return 	dict_of_genes,list_of_genes,fasta_dict		


def make_mmseq_cluster():

	"""
	Require mmseqs2.
	Call mmseqs2 to perform a fast protein clustering using easy-cluster option
	citing mmseqs2 documentation:
	easy-cluster by default clusters the entries of a FASTA/FASTQ file using a 
	cascaded clustering algorithm
	
	"""

	print("making MMseqs2 clustering")

	mmseq_proceso = Popen(['mmseqs','easy-cluster','DB.fasta','clusterRes','tmp'],stdout=PIPE,stderr=PIPE)
	error = mmseq_proceso.stderr.read()
	done = mmseq_proceso.stdout.read()

	if not error:
		print("MMseqs2 clustering done!")
	else:
		print("error in mmseq_cluster step!!\n")




def parser_DB_cluster(query_protein_list):

	"""
	parsea el archivo clusterRes_cluster.tsv, donde se encuentran todos
	los cluster que ha encontrado mmseq, y devuelve un diccionario(clusters_dict)
	que contiene clave(representative)-valor(protein members) y tambien devuelve
	el cluster del gen query que hemos introducido

	"""

	locus_tag,sequence = query_protein_list

	MMseq_db = open("clusterRes_cluster.tsv","r")

	clusters_dict = {} #representative of the cluster(Key);member_list(values)

	for line in MMseq_db:
		fields = line.rstrip().split("\t")
		representative = fields[0]
		member = fields[1]

		if representative not in clusters_dict:
			clusters_dict[representative]=[member]
		else:
			member_list = clusters_dict[representative]
			member_list.append(member)			
			clusters_dict[representative]=member_list

	count = 0
	for k,v in clusters_dict.items():
		count+=1
		if locus_tag in k or locus_tag in v:
			#print(k,v)
			query_cluster = v
			query_cluster.remove(locus_tag)
	
	print("MMseqs2 # cluster "+str(count))


	return clusters_dict,query_cluster



def get_genomic_neighbourhood_list(query_cluster,list_of_genes,nside):

	"""
	retrieve neighbours genes for all the genes presented in the query cluster
	nside control the number of neighbour genes to get up/downstream of anchor gene(query protein)
	The output is a dictionary (Key: anchor) : neighbour genes list [-2,-1,0,+1,+2]

	"""

	neighbourhood_dict = {}

	for protein in query_cluster:
		genomic_neighbourhood_list = []
		gene_position = list_of_genes.index(protein)
		for n in range(-nside,nside+1,1):
			pos = gene_position + n
			genomic_neighbourhood_list.append(list_of_genes[pos])

		neighbourhood_dict[protein] = genomic_neighbourhood_list


	return neighbourhood_dict




def calculate_clusters(cluster_dict,neighbourhood_dict,dict_of_genes,neigh_protein_sequences):

	"""
	Compute gene neighbour clusters
	------------------------------- 
	Calculate all uniq cluster, assign a numeric ID for each cluster and 
	annotate each neighbour gene with its numeric ID cluster in 
	genomic_cluster_dict

	"""
	
	neighbour_cluster_dict = {} # key (neigh gene) : cluster representative
	uniq_cluster_list = [] 		# contains all the uniq clusters
	genomic_cluster_dict = {}	# contain the numeric ID for every cluster key(cluster): int ID
	genomic_context_dict = {}

	# genero un diccionaro que me alamecene para cada protein center de cada genoma
	# a analizar su representative de su cluster
	for protein_center,neighbours_list in neighbourhood_dict.items():
		for gene in neighbours_list:
			neighbour_cluster_dict[gene]="Pseudo"	
			for representative,list_of_members in cluster_dict.items():
				if gene in list_of_members:
					neighbour_cluster_dict[gene]= representative
							
	
	# create a list with all the uniqes clusters
	
	uniq_cluster_list = []
	for cluster in neighbour_cluster_dict.values():
		if cluster not in uniq_cluster_list:
			uniq_cluster_list.append(cluster)
	
	
	
	#print(uniq_cluster_list)
	# assign an uniq (int) identificator (ID) for each representative. 
	# Cluster number is stored in synteny_cluster_dict
	cluster_number = 0
	for cluster in uniq_cluster_list:
		if cluster == "Pseudo": # Pseudogenes have assigned exclusively with cluster 0
			genomic_cluster_dict[cluster]=0
		else:	
			cluster_number+=1	
			genomic_cluster_dict[cluster]=cluster_number
	

	#Creo un diccionario que contiene todos los ortologos a la proteina query y sus contextos genomicos
	#en los diferentes genomas
	for gene,cluster_representative in neighbour_cluster_dict.items():
		strand,product,genome,gene_symbol = dict_of_genes[gene]
		cluster_number = genomic_cluster_dict[cluster_representative]
		sequence = neigh_protein_sequences[gene]
		genomic_context_dict[gene]= cluster_number,strand,product,genome,gene_symbol,sequence

	#clean comma character for each field in the dictionary
	clean_genomic_context_dict = clean_comma(genomic_context_dict)

	return clean_genomic_context_dict 


def calculate_conservation(genomic_context_dict,neighbourhood_dict):

	conservation_count_list = []
	count_of_anchors = 0
	conservation_dict = {}

	# generate a list with all clusters present in the neighbours genes
	# avoiding those repeated in the same anchor
	for anchor, neighbours_list in neighbourhood_dict.items():
		count_of_anchors +=1
		clusters_per_anchor= []
		for gene in neighbours_list:
			cluster = genomic_context_dict[gene][0]

			if cluster not in clusters_per_anchor:
				clusters_per_anchor.append(cluster) 
		#print(clusters_per_anchor)
		conservation_count_list += clusters_per_anchor

	# create a dictionary with cluster and their % of conservation 
	cluster_counts = Counter(conservation_count_list)
	for cluster,count in cluster_counts.items():
		conservation_dict[cluster] = count/count_of_anchors*100


	for k,v in genomic_context_dict.items():
		cluster = v[0]
		conservation_score = conservation_dict[cluster]
		v = v + [conservation_score]
		genomic_context_dict[k]=v
	
	#print(genomic_context_dict)

	return genomic_context_dict


def neighbour_fasta_dict(fasta_dict,neighbourhood_dict):
	"""
	create a fasta dict containning all neighbours protein sequences
	"""

	neigh_protein_sequences = {}
	for anchor_gene,genomic_neighbourhood_list in neighbourhood_dict.items():
		for gene in genomic_neighbourhood_list:
			try:
				sequence = fasta_dict[gene]
				neigh_protein_sequences[gene] = sequence
			except:
				neigh_protein_sequences[gene] = "NA"
	
	return neigh_protein_sequences 


def search_prosite_domains(neigh_protein_sequences, prosite = "False"):
	
	"""
    Recorre el diccionario que contiene todas las secuencias proteícas
	de los neighbours y busca los patrones de Prosite usando finditer.
	Usa las funciones del module Prosite
    """

	prosite_hash = {}

	import re

	
	#if prosite False this function return a domain_dict with only "No domain"
	if prosite == "False":
		prosite_hash = {}
		domain_list = "No domain"
		for protein,sequence in neigh_protein_sequences.items():
			prosite_hash[protein]=domain_list
	
	#if prosite True this function return a domain_dict with RE domains
	else:
	
		dict_pattern = create_prosite_dict()
				
		for protein,sequence in neigh_protein_sequences.items():	
			domain_list = []
			for pattern in list(dict_pattern.keys()):
				ungapped_seq = sequence.replace("-", "")
				matches = re.finditer(pattern, ungapped_seq)       
				for match in matches:
					# FORMAT: 'PS00538', 'CHEMOTAXIS_TRANSDUC_1', 'Bacterial chemotaxis sensory transducers signature.'
					accesion,name,description = dict_pattern[pattern]
					protein_match =accesion+"|"+name+"|"+description+"|"+str(match.start())+"|"+str(match.end())
					#print(protein_match)
					domain_list.append(protein_match)
			if domain_list == []:
				domain_list = "No domain"	
			prosite_hash[protein]=domain_list

	return prosite_hash


def genomic_context_output(neighbourhood_dict,genomic_context_dict,dict_of_genes,prosite_hash):
	
	"""
	Control strand orientation of Anchor genes making all anchor genes have positive orientation
	and change concordantly the orientation of the genomic context. 
	
	Join all prosite list result using a "@" in case prosite = True

	and make up output result to be saved in "gcontext.csv"
	"""


	# output file
	output = open("gcontext.csv","w")
	output.write("#ORF,cluster_ID,strand,description,gene_symbol,genome_ID,domain,sequence\n")

	for anchor_gene,genomic_neighbourhood_list in neighbourhood_dict.items():
		strand_order = 1
		anchor_strand = dict_of_genes[anchor_gene][0] # get anchor stand
		if anchor_strand == -1: # force anchor gene to have a positive strand positive orientation
			genomic_neighbourhood_list = genomic_neighbourhood_list[::-1]
			strand_order = -1
		for gene in genomic_neighbourhood_list:
			cluster_number,strand,description,genome,gene_symbol,sequence,conservation_score = genomic_context_dict[gene]
			strand = strand_order * strand #multiplicamos por -1 para que los strands cambien
			domains = prosite_hash[gene]
			if domains != "No domain":
				domains = "@".join(prosite_hash[gene])
			if cluster_number == 0: #its Pseudo gene
				description = "Pseudogen; "+description
				domains = "No domain"
			output.write("{},{},{},{},{},{},{},{},{}\n".format(gene,cluster_number,strand,description,gene_symbol,genome,domains,sequence,conservation_score))

	output.close()


def clean_comma(genomic_context_dict):

	"""
	Sice the output is in csv format is a good idea to scan and remove any comma
	presents in the strings to plot 
	"""

	for gene in genomic_context_dict.keys():
		clean_list = []
		for element in genomic_context_dict[gene]:
			if type(element) == str: 
				element = element.replace(",","")
			clean_list.append(element)
		genomic_context_dict[gene]=clean_list
	
	return genomic_context_dict

def mafft(neighbourhood_dict,fasta_dict):

	print("mafft alignment")

	output = open("anchor_cluster.fasta","w")

	for anchor_gene in neighbourhood_dict.keys():
		fasta = fasta_dict[anchor_gene]
		output.write(">"+anchor_gene+"\n"+fasta+"\n")
	output.close()

	mafft_process = Popen(['mafft-linsi anchor_cluster.fasta > anchor_cluster_mafft.fasta'],stdout=PIPE,stderr=PIPE,shell=True)
	
	error = mafft_process.stderr.read()
	done = mafft_process.stdout.read()

	if error:
		print(error)

def Fasttree():
	fasttree = Popen(['fasttree anchor_cluster_mafft.fasta > anchor_cluster_mafft.nwx'],stdout=PIPE,stderr=PIPE,shell=True)
	error = fasttree.stderr.read()
	done = fasttree.stdout.read()
	
	if error:
		print(error)


def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""

    # from whichcraft import which
    from shutil import which

    return which(name) is not None


def check_binaries(binaries_required):

	""" Check if binary is installed or in the path, if not sys.exit """

	for binary in binaries_required:
		
		if is_tool(binary) == False:
			sys.exit("\nError: ETE context requires "+ binary + " please install it or put in the path")
		

def check_biopython():

	try: 
		import Bio 
	except:
		sys.exit("\nError: ETE context requires "+ "Biopython" + " please install it or put in the path" )	


## Argument control with argparse
parser = argparse.ArgumentParser()
parser.add_argument("query_protein",
					help="query protein in fasta format")
parser.add_argument("gbff_folder",
					help="path to folder containing genomic files in .gbff")
parser.add_argument("-p","--prosite", 
					help="annotate proteins with prosites REs (default = False)", 
					action="store_true")
parser.add_argument("-n","--nside", 
					help="window of neighbours genes visualizated up/downstream (integer value)", 
					type=int)

args = parser.parse_args()

prosite = "True" if args.prosite else "False"
nside = args.nside if args.nside else 3

if args.query_protein:
	query_protein = args.query_protein

if args.gbff_folder:
	gbff_folder = args.gbff_folder


## Check for binaries and python packages
binaries_required = [ "fasttree",
					"mafft"]

check_binaries(binaries_required)
check_biopython()



query_protein_list = query_protein_parser(query_protein)
dict_of_genes,list_of_genes,fasta_dict	= make_list_and_dict(query_protein_list)
make_mmseq_cluster()
clusters_dict,query_cluster = parser_DB_cluster(query_protein_list)
neighbourhood_dict = get_genomic_neighbourhood_list(query_cluster,list_of_genes,nside)
neigh_protein_sequences = neighbour_fasta_dict(fasta_dict,neighbourhood_dict)
genomic_context_dict = calculate_clusters(clusters_dict,neighbourhood_dict,dict_of_genes,neigh_protein_sequences)
genomic_context_dict = calculate_conservation(genomic_context_dict,neighbourhood_dict)
prosite_hash = search_prosite_domains(neigh_protein_sequences,prosite)
genomic_context_output(neighbourhood_dict,genomic_context_dict,dict_of_genes,prosite_hash)
mafft(neighbourhood_dict,fasta_dict)
Fasttree()



print("--- %s seconds ---" % (time.time() - start_time))

