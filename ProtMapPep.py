
"""
ProtMapPep.py

Created by Pathmanaban Ramasamy on 18th Sep 2018
Copyright (c) 2018 Pathmanaban Ramasamy. All rights reserved.

"""
import sys, getopt
import csv
import os
import argparse
import difflib
from Bio import SeqIO
from collections import defaultdict


def UPseq():
	records = list(SeqIO.parse("uniprot_sprot.fasta", "fasta"))
	return records
UPseqlis=UPseq()				
	
def UPseqfetch(UPseqlis,accsn,pep):

	for record in UPseqlis:
			recid=record.id.split('|')[1]
			if str(recid)==str(accsn).strip():
				
				upseq=record.seq
				
				newpos1= [n for n in xrange(len(upseq)) if upseq.find(pep, n) == n]
				if newpos1:
					newpos1=[x+1 for x in newpos1]
					
					return ("No",newpos1)
				
				
def extract_mod(configfile):
	with open(configfile,'r') as infilex:
		next(infilex,None)
		modific=[]
		for line in infilex:
			line=line.split(',')
			mod=line[0].split('=')[1].lower()
			if mod not in modific:
				modific.append(mod)
			
		
			
		return modific
		
	



def map_and_locate():
	
	
	parser = argparse.ArgumentParser(description='Pepmap ')
	parser.add_argument('-s','--sp', help='Specify the species name(as in database): case sensitive', required=True)
	parser.add_argument('-i','--ifile', help='Specify the infile name with extension', required=True)
	parser.add_argument('-o','--ofile', help='Specify the outfile name with extension', required=True)
	parser.add_argument('-p','--pep',help='Column index of peptide in infile',required=True)
	parser.add_argument('-m','--mod',help='column index of modified peptide in infile',required=True)
	parser.add_argument('-a','--acc',help='column index of protein accession in infile',required=True)
	parser.add_argument('-c','--conf',help='Config file used for searching',required=True)
	args = vars(parser.parse_args())
	
	if len(args)==7:
		infile=args['ifile']
		outfile=args['ofile']
		propepindx=int(args['pep'])
		pepindx=int(args['mod'])
		accession=int(args['acc']) 
		species='_'+args['sp']
		ptmfile=args['conf']
		modifi=extract_mod(ptmfile)
		
		with open(infile,'r') as infile1,open(outfile,'w') as outfile:
				writer1 = csv.writer(outfile,delimiter='\t',quoting=csv.QUOTE_NONE,escapechar=' ')
				
				next(infile1,None)
				decoy=['_crap','Random_']
				an,bn=0,0
				for line1 in infile1:
					an+=1
					
					line=line1.split('\t')
					protpep=line[propepindx]
					pep=line[pepindx]
					
					if species in line[accession] and '_crap' not in line[accession] and 'Random_' not in line[accession]:
						
						if '|' in line[accession]:
							bn+=1
							accsn=line[accession].split('|')[1].strip()
							
							
							pep_start1=UPseqfetch(UPseqlis,accsn,protpep)
							if pep_start1:
								pep_start=pep_start1[1]
								mismatch=pep_start1[0]
						
								
								if pep_start:
									
									for ev_p_start in pep_start:
										ev_p_start=int(ev_p_start)
										peplis=[]
										modlis=[]
							
										if len(line[propepindx]) != len(line[pepindx]):
												
												for s in modifi:
												
													if s in pep.lower() and s not in modlis:
														
														modlis.append(s)
												
												for mod in modlis:		
													mod_pos_pep = [n for n in xrange(len(pep)) if pep.lower().find(mod, n) == n]
													
													for x in mod_pos_pep:
														peplis.append(mod+'__'+str(x+ev_p_start))
													
													peplis=sorted(peplis, key=lambda x: int(x.split('__')[1]))
													
												if len(peplis) >1:
													peplis1=[x.split('__')[0] for x in peplis]
													newpos=[]
													lis2=[len(x.split('__')[0]) for x in peplis]
													lis1=[int(x.split('__')[1]) for x in peplis]
													for elem in lis1:
														elem1=elem-sum(lis2[:lis1.index(elem)])
														newpos.append(elem1)
													
													peplis2=[]
													for a,b in enumerate(peplis1):
														
														elem2=b+'__'+str(newpos[a])
														peplis2.append(elem2)
													
													writer1.writerow([line1.strip('\n'),ev_p_start,','.join(peplis2),mismatch])
												elif len(peplis)==1:
													writer1.writerow([line1.rstrip('\n'),ev_p_start,','.join(peplis),mismatch])
										else:
											writer1.writerow([line1.strip('\n'),ev_p_start,"Unmodified",mismatch])
		print "\n"									
		print "Results are wirtten to: ", outfile
		print "\n"
		print "*****************************************************"
		print "              Thanks for using Pepmap                "
		print "*****************************************************"
	else:
		print "Some argument missing"
		sys.exit() 
		

if __name__ == "__main__":
	print "\n"
	print "************************************************"
	print "              Pepmap version 1.0                "
	print "************************************************"
	print "\n"
	map_and_locate()
	


