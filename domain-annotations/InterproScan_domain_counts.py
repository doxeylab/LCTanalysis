#!/bin/python3

######################################################################
#
#  Copyright Michael Mansfield, BSc
#
#  Parses certain annotations out of an InterproScan output file and
#  returns a file containing unique domain annotations with lists of
#  proteins that have that domain architecture.
#
#                         CAVEAT EMPTOR
#  Basically, this is not really meant to be used by anyone except
#  for the analysis specific to Orrell et al. 2019 (in review).
#  With some specific changes it might be used for other protein
#  families.
#
#  If you have found this script through other means, I recommend
#  looking elsewhere to solve your problem, or else tweet me and we
#  can work it out together: @Mansfi3ld
#
#####################################################################

import sys
import operator
from collections import OrderedDict
from operator import itemgetter

# All of the files and paths are hard-coded. You will have to 
# change each of these if you want to replicate this work.
base_directory = "/disk2/mmansfield/LCT/Round4/revisions/collection-analysis-stuff/domain-annotations"
fasta_file = base_directory + "/PSI-BLAST-Subset.fa"
interpro_output = base_directory + "/PSI-BLAST-Subset-Interproscan.tsv"
# Minimum protein length. Excludes proteins below a length cutoff from output.
minimum_length = 0

###############################################################
# The following lines attach annotation identifiers to
# text descriptions of their domains. However, some Gene3d
# models do not have good annotations, so I had to look these
# up manually and try to work backwards from the structures that
# constitute the model in order to annotate them.
#
# Not very sophisticated or future-proof, so please beware.
###############################################################

# I recycle some names here. I'm not proud of myself.
Gene3d_dict = {}
with open(base_directory + "/Gene3D_model_to_family_map.csv") as gene3d_file:
	for line in gene3d_file:
		line = line.strip()
		line = line.replace('"', '')
		pdbid = line.split(",")[0]
		gene3did = line.split(",")[1]
		desc = line.split(',')[2]
		other = line.split(',')[3]
		if gene3did not in Gene3d_dict:
			Gene3d_dict[gene3did] = [ desc, other, [pdbid] ]
		else:
			Gene3d_dict[gene3did][2].append(pdbid)

PfamA_dict = {}
with open(base_directory + "/Pfam-A.clans.tsv") as pfam_file:
	for line in pfam_file:
		line = line.strip()
		split_line = line.split('\t')
		PfamA_dict[split_line[0]] = split_line[1:]

CATH_dict = {}
with open(base_directory + "/CATH_names.txt") as cath_file:
	for line in cath_file:
		line = line.strip()
		if line.startswith("#"):
			pass
		else:
			split_line = list(filter(None, line.split(" ")))
			if split_line[0] not in CATH_dict:
				CATH_dict[split_line[0]] = [split_line[1], " ".join(split_line[2:]).split(":")[1] ]

# The piece de resistance: a manually-derived list of annotations
# based on me just looking at the Gene3d website.
# Do as I say, not as I do!

manual_gene3d_annots = {'1.10.10.1780':'LCT-like superfamily', '1.10.274.80':'C. novyi alpha toxin-associated superfamily', '1.10.287.130':'histidine kinase superfamily', '1.10.3730.30':'C. novyi alpha toxin-associated superfamily','1.20.140.180':'Dermonecrotic/RTX toxin, membrane localization domain', '1.20.58.1190':'TcdA/TcdB N-terminal helical domain', '1.25.40.20':"ankyrin repeat superfamily", '3.30.450.20':'PAS domain superfamily', '3.30.565.10':'Histidine kinase-like ATPase, C-terminal domain', '3.30.70.1720':'Adenylyl cyclase superfamily', '3.30.70.330':'giant nucleotide-binding superfamily', '3.40.50.11050':'peptidase C60 superfamily', '3.40.50.11550':'Region of Pasteurella dermonecrotic toxin', '3.90.550.20':'Glucosyltransferase/PaTox superfamily', '2.10.270.10':'Cholin binding', '3.10.670.10':'SseI-like cysteine protease', '3.30.160.280':'Heme-binding hemoglobin serine protease' }

###############################################################
# The following lines read in the Interproscan files. It makes
# some assumptions about the structure of the fasta headers, as
# well as the interproscan output.
###############################################################

protein_dict1 = {}
with open(fasta_file) as i:
	for line in i:
		line = line.strip()
		if line.startswith(">"):
			protein_id = line.split(">")[1].split(" ")[0]
			metadata = " ".join(line.split(" ")[1:])
			protein_dict1[protein_id] = [metadata]
		elif line.startswith("-"):
			pass
		else:
			protein_dict1[protein_id].append(len(line))

protein_dict = {}
for p in protein_dict1:
	if protein_dict1[p][1] >= minimum_length:
		protein_dict[p] = protein_dict1[p]

# I used this script at one point to generate protein domain diagrams, so
# the dictionary here is assigned values to HEX colour values.
# The current version only cares about the keys so this could be a list
# but it works so ¯\_(ツ)_/¯
annotations_I_care_about = {"Gene3D":"#ff6db6", "Pfam":"#074751", "CDD":"#db6d00", "ProSiteProfiles":"#b66dff", "TIGRFAM":"#ffff6d", "CDD":"#db6d00"}

#Interproscan reader.
with open(interpro_output) as i2:
	for line in i2:
		line = line.strip()
		if line.startswith("#"):
			pass
		elif len(line) < 1:
			pass
		else:
			# The filter is to remove empty columns.
			# (to avoid lists looking like this: ["A", "B", "", "", "", "X")
			split_line = list(filter(None, line.split("\t")))
			protein = split_line[0]
			# Basically, for every annotation, make a blank one, but replace
			# it with information parsed from the interproscan output.
			protein_dict_info = ["", "", "", 0, 0, ""]
			if protein in protein_dict:
				annotation_type = split_line[3]
				if annotation_type in annotations_I_care_about:
					if annotation_type == 'Pfam':
						model = split_line[4].split(":")[-1]
						desc = split_line[5]
						start = int(split_line[6])
						stop = int(split_line[7])
						evalue = split_line[8]
						protein_dict_info = [annotation_type, model, desc, start, stop, evalue]
					elif annotation_type == "Gene3D":
						model = split_line[4].split(":")[-1]
						gene3d_entry = Gene3d_dict[model]
						simple_desc = gene3d_entry[0]
						start = int(split_line[5])
						stop = int(split_line[6])
						evalue = split_line[7]
						if simple_desc == "":
							if model in manual_gene3d_annots:
								simple_desc = manual_gene3d_annots[model]
							else:
								print("Somehow there is no manual or automatic Gene3d annotation for the following: %s" % model)
								print("Look it up and add it to the manual_3d_dict object above...")
								print("The offending line:")
								print(split_line)
								sys.exit()
						protein_dict_info = [annotation_type, model, simple_desc, start, stop, evalue]
					elif annotation_type == 'ProSiteProfiles':
						model = split_line[4]
						desc = split_line[-1]
						start = int(split_line[6])
						stop = int(split_line[7])
						evalue = split_line[8]
						protein_dict_info = [annotation_type, model, desc, start, stop, evalue]
					elif annotation_type == 'TIGRFAM':
						model = split_line[4]
						desc = split_line[5]
						start = int(split_line[6])
						stop = int(split_line[7])
						protein_dict_info = [annotation_type, model, desc, start, stop, evalue]
					elif annotation_type == "CDD":
						model = split_line[4]
						desc = split_line[5]
						start = int(split_line[6])
						stop = int(split_line[7])
						evalue = split_line[8]
						protein_dict_info = [annotation_type, model, desc, start, stop, evalue]
					else:
						model = split_line[4]
						desc = split_line[5]
						start = int(split_line[6])
						stop = int(split_line[7])
						evalue = split_line[8]
						protein_dict_info = [annotation_type, model, desc, start, stop, evalue]
					protein_length = protein_dict[protein][1]
					#print(split_line)
					#print(protein_dict_info)
					if int(protein_dict_info[4]) > protein_length:
						print("Somehow the lengths of proteins got confused. Parsing error maybe? Check the following line:")
						print(split_line)
						print(protein_dict_info)
						sys.exit()
					if protein in protein_dict:
						if protein_dict_info == ["", "", "", 0,0, ""]:
							print("Yeah, something's wrong here, the following line is being annotated to the protein %s incorrectly" % protein)
							print(split_line)
							print(protein_dict_info)
							sys.exit()
						else:
							#print(split_line)
							#print(protein_dict_info)
							protein_dict[protein].append(protein_dict_info)
				else:
					pass

###############################################################
# Now that the annotations are all translated into usable
# text annotations, the dictionary must be reversed. 
# Currently, the dict is structured dict[protein] = [annotations]
# But what we want is dict[Domain_Architecture] = list(proteins)
###############################################################

unique_domain_architectures = {}
for protein in protein_dict:
	fasta_head = protein_dict[protein][0]
	protein_length = protein_dict[protein][1]
	domain_arch_redundant = []
	domain_arch_nonredundant = []
	for domain in sorted( protein_dict[protein][2:], key=itemgetter(3)):
		domain_formatted = domain[2] +  " (" + domain[1] + ")"
		domain_arch_redundant.append( domain_formatted )
		if domain_formatted  not in domain_arch_nonredundant:
			domain_arch_nonredundant.append(domain_formatted)
	#print(domain_arch_redundant)
	#domain_arch_nonredundant = list(set(list(domain_arch_redundant)))
	#print(domain_arch_nonredundant)
	domain_arch_str = ','.join(domain_arch_nonredundant)

	if domain_arch_str not in unique_domain_architectures:
		unique_domain_architectures[domain_arch_str] = [1, protein]
	else:
		unique_domain_architectures[domain_arch_str][0] += 1
		unique_domain_architectures[domain_arch_str].append(protein)

for architecture in unique_domain_architectures:
	# Finally, prints the annotations. Yay.
	# Prints to stdout because I do everything at the command line :)
	# Obviously, change to with open('output.tsv', 'w') as output: etc.
	# if desired
	print(architecture + "\t" + str(unique_domain_architectures[architecture][0]) + "\t" + ', '.join(unique_domain_architectures[architecture][1:]))

