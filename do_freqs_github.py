#!/usr/bin/python

# Download all human variations in vcf format:
# ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/

# in .vcf file, header row starts with "#CHROM"
# Col: #CHROM (ex: chr7)
# Col: POS (ex: 33204597)
# Col: ID (ex: rs4491797)
# Col: REF (ex: G)
# Col: ALT (ex: T)
# Col: FILTER (ex: PASS)


import re
import csv
import os
import string

def openfile(filename, cols):
	raw = list(csv.reader(open(filename, 'rU'), delimiter='\t'))
	if len(cols) != 0:
		headers = raw[0]
		raw = raw[1:]

		count = 0
		inds = []
		for a in headers:
			for b in cols:
				if a == b:
					inds.append(count)

			count += 1
	else:
		headers = []
		inds = []

	return raw, headers, inds


def savefile(filename, out, headers_out):

				if headers_out != []:
								out.insert(0,headers_out)

				with open(filename, 'w') as a_file:
								for result in out:
												result = '\t'.join(result)
												a_file.write(result + '\n')





def gen_vcf(allele_list, full_data, outfile):
	cols = ['#CHROM', 'POS', 'REF', 'ALT']
#	raw = openfile(os.path.expanduser(full_data), [])

#	count = 0
#	row_ind = 0
#	for a in raw:
#		for b in a:
#			if len(b) > 0:
	#			print b
#				if b[0] == '#CHROM':
		#			print b, count
#					row_ind = count
#			count += 1

#	top_info = raw[0][:row_ind]
#	headers = raw[0][row_ind]
#	raw = raw[0][row_ind:]

#	count = 0
#	for a in headers:
#		if a == '#CHROM':
#			chrom_ind = count
#		if a == 'POS':
#			pos_ind = count
#		if a == 'REF':
#			ref_ind = count
#		if a == 'ALT':
#			alt_ind = count
#		count += 1

#	out = top_info
#	out.append(headers)

	chrom_ind = 0
	pos_ind = 1
	rs_ind = 2
	ref_ind = 3
	alt_ind = 4


	def strip_newline(line_in):
		new_line = []
		for x in line_in:
			if '\n' not in x:
				new_line.append(x)
			else:
				x_split = x.split('\n')
				new_line.append(x_split[0])
		return new_line

	out = []
	found = 0

#	for a in raw:
	count = 0
	with open(os.path.expanduser(full_data)) as f:
		for line in f:
			count += 1
#			print 'Trying:', count, '(Found:', found, ')'

#			a = list(line)
			a = line.split('\t')
			#print a
			if len(a) > 2:
				print 'Trying:', count, a[chrom_ind], a[pos_ind], '(Found:', found, ')'
			if len(a) > 0:
				if '#' in a[0]:
					new_line = strip_newline(a)
					out.append(new_line)

				else:
					chrom_num = int(a[chrom_ind].split('chr')[-1])
					if chrom_num > 10:
						break
					
					
	#		print a[chrom_ind], a[pos_ind], a[ref_ind], a[alt_ind]
			for b in allele_list:
#				if a[chrom_ind] == b[1] or a[chrom_ind] == b[1].split('chr')[-1]:
#				if a[chrom_ind] == '6':
				if len(a) > 2:
					for i in range(0,len(b[2])):
						pos = b[3][i]
						chg_split = b[4][i].split('>')
						ref = chg_split[0]
						alt = chg_split[1]
						if pos == a[pos_ind] and ref == a[ref_ind].upper() and alt == a[alt_ind].upper():
							new_line = strip_newline(a)
							out.append(new_line)
							print ''
							print b[0], 'found:'	
							print new_line
							found += 1
						if a[rs_ind] == b[2][i]:
							new_line = strip_newline(a)
							out.append(new_line)
							print ''
							print b[0], 'found:'	
							print new_line
							found += 1				


	out2 = out
	out = []
	for a in out2:
		if a not in out:
			out.append(a)
							
						
							
							
	savefile(outfile, out, [])





def extract_af(allele_list, vcf_in):
	# Header = 'INFO'

	vcf_raw = openfile(vcf_in, [])

	count = 0
	row_ind = 0
	for a in vcf_raw:
		for b in a:
			if len(b) > 0:
				if b[0] == '#CHROM':
					row_ind = count
			count += 1

	top_info = vcf_raw[0][:row_ind]
	headers = vcf_raw[0][row_ind]
	vcf_data = vcf_raw[0][row_ind+1:]

	pos_ind = 1
	id_ind = 2
	ref_ind = 3
	alt_ind = 4
	info_ind = 7
	start_ind = 9


	seen = []
	for a in allele_list:
		for snp in a[2]:
			found = False
			for b in vcf_data:
				if snp == b[id_ind]:
					found = True
			if found == False and snp not in seen:
				print 'Missing SNP:', snp, '('+a[0]+')'
				seen.append(snp)


	# get population freqs for multi-SNP alleles of DR
	# each individual has diploid readout: e.g. 0|1
	# 0 = REF (not SNP)
	# 1 = ALT (SNP)
	# for DR, consider both homozygous (1|1) and hetero (0|1 and 1|0)
	# if an individual has either case, they "have" the SNP
	# then group all individuals that have all appropriate SNPs
	# the allele freq can then be calculated

	all_ins, h, i = openfile('individuals_ALL.txt', [])
	afr_ins, h, i = openfile('individuals_AFR.txt', [])
	amr_ins, h, i = openfile('individuals_AMR.txt', [])
	eur_ins, h, i = openfile('individuals_EUR.txt', [])
	eas_ins, h, i = openfile('individuals_EAS.txt', [])
	sas_ins, h, i = openfile('individuals_SAS.txt', [])

	indiv_heads = ['Individual', 'Population']
	for a in allele_list:
		indiv_heads.append(a[0] + ' hetero_1')
		indiv_heads.append(a[0] + ' hetero_2')
		indiv_heads.append(a[0] + ' homo')
		indiv_heads.append(a[0] + ' total')

	individs = []
	for i in range(start_ind,len(headers)):
		entry = ['0']*len(indiv_heads)
		entry[0] = headers[i]

		pop = ''		
		if [entry[0]] in afr_ins:
			pop = 'AFR'
		if [entry[0]] in amr_ins:
			pop = 'AMR'
		if [entry[0]] in eur_ins:
			pop = 'EUR'
		if [entry[0]] in eas_ins:
			pop = 'EAS'
		if [entry[0]] in sas_ins:
			pop = 'SAS'
		entry[1] = pop

		start_ind = 2
		for x in allele_list:
			snp_pass = False
			homo = 0
			hetero_1 = 0
			hetero_2 = 0

			snps = x[2]
			for snp in snps:
				for y in vcf_data:
					if y[2] == snp:
						#print snp, y[i]
						#if y[i] != '0|0':
						#	snp_pass += 1
						if y[i] == '1|1':
							homo += 1
						elif y[i] == '1|0':
							hetero_1 += 1
						elif y[i] == '0|1':
							hetero_2 += 1
			#print entry[0], snps, snp_pass

			if homo == len(snps):
				snp_pass = True
				entry[start_ind+2] = '1'
			if hetero_1 == len(snps):
				snp_pass = True
				entry[start_ind] = '1'
			if hetero_2 == len(snps):
				snp_pass = True
				entry[start_ind+1] = '1'

			if snp_pass == True:
				entry[start_ind+3] = '1'
			#else:
			#	entry[start_ind] = '0'
			start_ind += 4

		#print entry
		individs.append(entry)
	out_file = 'DRO_allele_freqs_individuals.txt'
	savefile(out_file, individs, indiv_heads)
	individs = individs[1:]


	headers_out = [
			'Allele', 'SNPs',
			'ALL Allele freq',
			'AFR Allele freq',
			'AMR Allele freq',
			'EAS Allele freq',
			'EUR Allele freq',
			'SAS Allele freq'
			] 
	out = []

	for a in allele_list:
		entry = ['']*len(headers_out)
		entry[0] = a[0]

		snp_str = ''
		for x in a[2]:
			snp_str += x + ';'
		snp_str = snp_str[:-1]
		entry[1] = snp_str

		for i in range(0,len(indiv_heads)):
			if a[0] + ' total' == indiv_heads[i]:
				ind = i
		# ALL
		pass_num = 0
		total = 0
		for b in individs:
			if b[ind] == '1':
				pass_num += 1
			total += 1
		freq = float(pass_num) / float(total)
		print entry[0], snp_str, 'ALL:', pass_num, '/', total, '('+str(freq)+'%)'
		entry[2] = str(freq)

		start_ind = 3
		pops = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
		for pop in pops:
			pass_num = 0
			total = 0
			for b in individs:
				if b[1] == pop:
					if b[ind] == '1':
						pass_num += 1
					total += 1

			freq = float(pass_num) / float(total)
			print entry[0], snp_str, pop+':', pass_num, '/', total, '('+str(freq)+'%)'
			entry[start_ind] = str(freq)
			start_ind += 1
		out.append(entry)

	out_file = 'DRO_allele_freqs.txt'
	savefile(out_file, out, headers_out)


	headers_out = ['SNP', 'Allele', 'Position (chr6)', 'Ref>Alt (nt)',
			'Ref>Alt (AA)',
			'AF',
			'EAS_AF', 'AMR_AF', 'AFR_AF', 'EUR_AF',
			'SAS_AF',
			'Full info'
			]
	out = []

	for a in vcf_data:
		entry = ['']*len(headers_out)
		entry[0] = a[id_ind]
		entry[2] = a[pos_ind]
		entry[3] = a[ref_ind] + '>' + a[alt_ind]

		info = a[info_ind]
		for i in info.split(';'):
			if i[:3] == 'AF=':
				entry[5] = i.split('AF=')[-1]
#				print entry[0], entry[5]			
			if 'EAS_AF' in i:
				entry[6] = i.split('AF=')[-1]
			if 'AMR_AF' in i:
				entry[7] = i.split('AF=')[-1]
			if 'AFR_AF' in i:
				entry[8] = i.split('AF=')[-1]
			if 'EUR_AF' in i:
				entry[9] = i.split('AF=')[-1]
			if 'SAS_AF' in i:
				entry[10] = i.split('AF=')[-1]

		entry[11] = info

		allele = 'n.d.'
		aa = 'n.d.'
		for b in allele_list:
			#print entry[0], b
			for c in range(0,len( b[2])):
				if b[2][c] == entry[0]:
					aa = b[5][c]
					allele = b[0]
		entry[1] = allele
		entry[4] = aa
		#print entry
		
		out.append(entry)

	#savefile(out_file, out, headers_out)	





def split_multiallelic(in_file):
	raw, headers, inds = openfile(in_file, [])
#	raw = raw[0]

	count = 0
	for a in raw:
#		print a
		if len(a) > 0:
			if a[0] == '#CHROM':
				head_ind = count

		count += 1

	top_info = raw[:head_ind]
	headers = raw[head_ind]
	raw = raw[head_ind+1:]

	count = 0
	for a in headers:
		if a == 'INFO':
			info_ind = count
		if a == 'ID':
			id_ind = count
		if a == 'REF':
			ref_ind = count
		if a == 'ALT':
			alt_ind = count
		if a == 'FORMAT':
			format_ind = count


		count += 1

	out = top_info
	out.append(headers)

	for a in raw:
		if 'MULTI_ALLELIC' in a[info_ind]:
			entry_1 = ['']*len(headers)
			entry_2 = ['']*len(headers)

			entry_1[0] = a[0]
			entry_2[0] = a[0]

			entry_1[1] = a[1]
			entry_2[1] = a[1]

			entry_1[2] = a[2] + '.1'
			entry_2[2] = a[2] + '.2'

			entry_1[3] = a[3]
			entry_2[3] = a[3]

			entry_1[4] = a[4].split(',')[0]
			entry_2[4] = a[4].split(',')[1]

			entry_1[5] = a[5]
			entry_2[5] = a[5]

			entry_1[6] = a[6]
			entry_2[6] = a[6]

			info_1 = ''
			info_2 = ''
#			print a[info_ind]

			for x in a[info_ind].split(';'):
				comma_split = x.split(',')
#				print comma_split
				if len(comma_split) > 1:
					front_sp = comma_split[0].split('=')
					front = front_sp[0] + '='
					info_1 += front + front_sp[-1] + ';'
					info_2 += front + comma_split[1] + ';'
				else:
					info_1 += x + ';'
					info_2 += x + ';'

			info_1 = info_1[:-1]
			info_2 = info_2[:-1]

			entry_1[7] = info_1
			entry_2[7] = info_2

			entry_1[8] = a[8]
			entry_2[8] = a[8]

			for x in range(9,len(a[9:])):
				geno_1a = ''
				geno_1b = ''
				geno_2a = ''
				geno_2b = ''
				geno_sp = a[x].split('|')

				if geno_sp[0] == '0':
					geno_1a = '0'
					geno_2a = '0'
				elif geno_sp[0] == '1':
					geno_1a = '1'
					geno_2a = '0'
				elif geno_sp[0] == '2':
					geno_1a = '0'
					geno_2a = '1'

				if geno_sp[1] == '0':
					geno_1b = '0'
					geno_2b = '0'
				elif geno_sp[1] == '1':
					geno_1b = '1'
					geno_2b = '0'
				elif geno_sp[1] == '2':
					geno_1b = '0'
					geno_2b = '1'

				geno_1 = geno_1a + '|' + geno_1b
				geno_2 = geno_2a + '|' + geno_2b

				#print a[x], geno_1, geno_2
				entry_1[x] = geno_1
				entry_2[x] = geno_2

			out.append(entry_1)
			out.append(entry_2)

		else:
			out.append(a)

	out_file = '1_alleles_multisplit.vcf'
	savefile(out_file, out, [])

#split_multiallelic('1_DM_alleles_raw_vcf.vcf')




def get_allele_freqs(freq_file, blocks):
	cols = ['Population']
	raw, headers, inds = openfile(freq_file, cols)

	for i in range(0,len(blocks)):
		#block_let = letters[i]
		block = blocks[i]
		neg_block = block[0][0].split('*')[0] + '*0101'
		headers.append(neg_block + ' hetero_1')
		headers.append(neg_block + ' hetero_2')
		headers.append(neg_block + ' homo')
		headers.append(neg_block + ' total')

	blocks2 = []

	raw2 = raw
	raw = []
	for a in raw2:
		entry = ['0']*len(headers)
		for i in range(0,len(a)):
			entry[i] = a[i]
		
		for i in range(0,len(blocks)):
			#block_let = letters[i]
			block = blocks[i]
			neg_block = block[0][0].split('*')[0] + '*0101'

			for j in range(0,len(headers)):
				if neg_block + ' hetero_1' == headers[j]:
					neg_het_1_ind = j
				if neg_block + ' hetero_2' == headers[j]:
					neg_het_2_ind = j
				if neg_block + ' homo' == headers[j]:
					neg_hom_ind = j
				if neg_block + ' total' == headers[j]:
					neg_tot_ind = j
			block2 = []
			for x in block:
				block2.append(x)
			block2.append([neg_block, ''])
			block2 = sorted(block2)
			if block2 not in blocks2:
				blocks2.append(block2)

			found_het_1 = False
			found_het_2 = False
			found_hom = False
			het_1_num = 0
			het_2_num = 0
			for c in block:
				allele = c[0]
				het_ind = ''
				hom_ind = ''
				tot_ind = ''
				print allele
				for j in range(0,len(headers)):
					if allele + ' hetero_1' == headers[j]:
						het_1_ind = j
					if allele + ' hetero_2' == headers[j]:
						het_2_ind = j
					if allele + ' homo' == headers[j]:
						hom_ind = j
					if allele + ' total' == headers[j]:
						tot_ind = j
				if a[het_1_ind] != '0':
					found_het_1 = True
					het_1_num += 1
				if a[het_2_ind] != '0':
					found_het_2 = True
					het_2_num += 1
				if a[hom_ind] != '0':
					found_hom = True
			
			if found_het_1 == False and found_het_2 == False and found_hom == False:
				entry[neg_hom_ind] = '1'
				entry[neg_tot_ind] = '1'
			if found_het_1 == True and found_het_2 == False and found_hom == False:
				entry[neg_het_2_ind] = '1'
				entry[neg_tot_ind] = '1'
			if found_het_1 == False and found_het_2 == True and found_hom == False:
				entry[neg_het_1_ind] = '1'
				entry[neg_tot_ind] = '1'
	
			#if found_het_1 == True and found_het_2 == False and found_hom == False:
			#	if het_1_num == 1:
			#		entry[neg_het_1_ind] = '1'
			#		entry[neg_tot_ind] = '1'
			#elif found_het_1 == False and found_het_2 == True and found_hom == False:
			#	if het_2_num == 1:
			#		entry[neg_het_2_ind] = '1'
			#		entry[neg_tot_ind] = '1'
			#elif found_het_1 == False and found_het_2 == False and found_hom == True:
			#	entry[neg_hom_ind] = '1'
			#	entry[neg_tot_ind] = '1'
				
		raw.append(entry)


	unique_pops = []
	for a in raw:
		unique_pops.append(a[inds[0]])
	unique_pops = list(set(unique_pops))
	unique_pops = sorted(unique_pops)
	unique_pops.insert(0,'ALL')

	out = []
	out.append([''])
	letters = [i for i in string.ascii_uppercase]

	for pop in unique_pops:
		out.append([pop])

		if pop != 'ALL':
			raw_to_use = []
			for a in raw:
				if a[inds[0]] == pop:
					raw_to_use.append(a)
		else:
			raw_to_use = raw

		for i in range(0,len(blocks2)):
			print blocks2[i]
			#print i, letters
			block_let = letters[i]
			block = blocks2[i]

			out.append(['Block ' + block_let])
			
			for a in block:
				allele = a[0]
				het_1_ind = ''
				het_2_ind = ''
				hom_ind = ''
				tot_ind = ''
				for j in range(0,len(headers)):
					if allele + ' hetero_1' == headers[j]:
						het_1_ind = j
					if allele + ' hetero_2' == headers[j]:
						het_2_ind = j
					if allele + ' homo' == headers[j]:
						hom_ind = j
					if allele + ' total' == headers[j]:
						tot_ind = j

				het_1_tot = 0
				het_2_tot = 0
				hom_tot = 0
				tot_tot = 0
				for b in raw_to_use:
					if b[het_1_ind] != '0':
						het_1_tot += 1	
					if b[het_2_ind] != '0':
						het_2_tot += 1	
					if b[hom_ind] != '0':
						hom_tot += 1	
					if b[tot_ind] != '0':
						tot_tot += 1	
				het_1_per = float(het_1_tot) / float(len(raw_to_use))			
				het_2_per = float(het_2_tot) / float(len(raw_to_use))			
				hom_per = float(hom_tot) / float(len(raw_to_use))			
				tot_per = float(tot_tot) / float(len(raw_to_use))			

	
				out.append(['  Individuals with ' + allele])
				out.append(['  Heterozygous (0|1): ' + str(het_1_tot) + ' (' + str(het_1_per) + ')'])
				out.append(['  Heterozygous (1|0): ' + str(het_2_tot) + ' (' + str(het_2_per) + ')'])
				out.append(['  Homozygous (1|1): ' + str(hom_tot) + ' (' + str(hom_per) + ')'])
				out.append(['  Total: ' + str(tot_tot) + ' (' + str(tot_per) + ')'])

				out.append([''])

			
			out.append(['** Block combinations'])
			out.append([''])

			block_combs = []
			for j in block:
				for h in block:
					comb = [j[0], h[0]]
					comb = sorted(comb)
					if comb not in block_combs:
						block_combs.append(comb)

			for comb in block_combs:
				all_1 = comb[0]
				all_2 = comb[1]
				#print '*', comb
					
				het_1_ind_1 = ''
				het_2_ind_1 = ''
				hom_ind_1 = ''
				tot_ind_1 = ''
				het_1_ind_2 = ''
				het_2_ind_2 = ''
				hom_ind_2 = ''
				tot_ind_2 = ''
				for j in range(0,len(headers)):
						if all_1 + ' hetero_1' == headers[j]:
							het_1_ind_1 = j
						if all_1 + ' hetero_2' == headers[j]:
							het_2_ind_1 = j
						if all_1 + ' homo' == headers[j]:
							hom_ind_1 = j
						if all_1 + ' total' == headers[j]:
							tot_ind_1 = j
						if all_2 + ' hetero_1' == headers[j]:
							het_1_ind_2 = j
						if all_2 + ' hetero_2' == headers[j]:
							het_2_ind_2 = j
						if all_2 + ' homo' == headers[j]:
							hom_ind_2 = j
						if all_2 + ' total' == headers[j]:
							tot_ind_2 = j
				c_11 = 0
				c_01 = 0
				c_10 = 0
				c_00 = 0	
				for b in raw_to_use:
					test_1 = False
					test_2 = False

					if b[hom_ind_1] != '0' or b[het_1_ind_1] != '0':
						test_1 = True
					if b[hom_ind_2] != '0' or b[het_2_ind_2] != '0':
						test_2 = True

					if test_1 == True and test_2 == True:
						c_11 += 1
					elif test_1 == False and test_2 == True:
						c_01 += 1
					elif test_1 == True and test_2 == False:
						c_10 += 1
					elif test_1 == False and test_2 == False:
						c_00 += 1

				c_11_per = float(c_11) / float(len(raw_to_use))
				c_01_per = float(c_01) / float(len(raw_to_use))
				c_10_per = float(c_10) / float(len(raw_to_use))
				c_00_per = float(c_00) / float(len(raw_to_use))


				out.append(['  ' + all_1 + '-' + all_2])
					
				out.append(['  Individuals with (' + all_1 + '|' + all_2 +') (1|1): ' + str(c_11) + ' (' + str(c_11_per) + ')'])	
				out.append(['  Individuals with (' + all_1 + '|' + all_2 +') (0|1): ' + str(c_01) + ' (' + str(c_01_per) + ')'])	
				out.append(['  Individuals with (' + all_1 + '|' + all_2 +') (1|0): ' + str(c_10) + ' (' + str(c_10_per) + ')'])	
				out.append(['  Individuals with (' + all_1 + '|' + all_2 +') (0|0): ' + str(c_00) + ' (' + str(c_00_per) + ')'])	

				out.append([''])


		#for i in range(0,len(blocks2)):
		if len(blocks2) == 2:
			#print blocks2[i]
			#print i, letters
			#block_let = letters[i]
			#block = blocks2[i]


			combs = []
			if len(blocks2) == 2:
				for e in blocks2[0]:
					for f in blocks2[1]:
						comb = [e[0], f[0]]
						#print comb	
						comb = sorted(comb)
						if comb not in combs:
							combs.append(comb)

				for comb in combs:
					all_1 = comb[0]
					all_2 = comb[1]
					
					het_1_ind_1 = ''
					het_2_ind_1 = ''
					hom_ind_1 = ''
					tot_ind_1 = ''
					het_1_ind_2 = ''
					het_2_ind_2 = ''
					hom_ind_2 = ''
					tot_ind_2 = ''
					for j in range(0,len(headers)):
						if all_1 + ' hetero_1' == headers[j]:
							het_1_ind_1 = j
						if all_1 + ' hetero_2' == headers[j]:
							het_2_ind_1 = j
						if all_1 + ' homo' == headers[j]:
							hom_ind_1 = j
						if all_1 + ' total' == headers[j]:
							tot_ind_1 = j
						if all_2 + ' hetero_1' == headers[j]:
							het_1_ind_2 = j
						if all_2 + ' hetero_2' == headers[j]:
							het_2_ind_2 = j
						if all_2 + ' homo' == headers[j]:
							hom_ind_2 = j
						if all_2 + ' total' == headers[j]:
							tot_ind_2 = j

					het_1_tot_1 = 0
					het_2_tot_1 = 0
					hom_tot_1 = 0
					tot_tot_1 = 0
					het_1_tot_2 = 0
					het_2_tot_2 = 0
					hom_tot_2 = 0
					tot_tot_2 = 0
					hom_1_2 = 0
					hom_het_1_2 = 0
					hom_het_mix_1_2 = 0
					hom_1_hom_2 = 0
					hom_1_het_2 = 0
					het_1_hom_2 = 0
					het_1_het_2 = 0
					tot_1_tot_2 = 0
					tot_1_het_2 = 0
					tot_1_hom_2 = 0
					het_1_tot_2 = 0
					hom_1_tot_2 = 0


					c_11_11 = 0
					c_01_10 = 0
					c_10_01 = 0
					c_11_01 = 0
					c_11_10 = 0
					c_01_11 = 0
					c_10_11 = 0
					c_00_00 = 0

					for b in raw_to_use:
						if b[hom_ind_1] != '0' and b[hom_ind_2] != '0':
							c_11_11 += 1
						if b[tot_ind_1] == '0' and b[tot_ind_2] == '0':
							c_00_00 += 1
						if b[het_2_ind_1] != '0' and b[het_1_ind_2] != '0':
							c_01_10 += 1
						if b[het_1_ind_1] != '0' and b[het_2_ind_2] != '0':
							c_10_01 += 1
						if b[hom_ind_1] != '0' and b[het_2_ind_2] != '0':
							c_11_01 += 1
						if b[hom_ind_1] != '0' and b[het_1_ind_2] != '0':
							c_11_10 += 1
						if b[het_2_ind_1] != '0' and b[hom_ind_2] != '0':
							c_01_11 += 1
						if b[het_1_ind_1] != '0' and b[hom_ind_2] != '0':
							c_10_11 += 1


						if b[het_1_ind_1] != '0':
							het_1_tot_1 += 1	
						if b[het_2_ind_1] != '0':
							het_2_tot_1 += 1	
						if b[hom_ind_1] != '0':
							hom_tot_1 += 1	
						if b[tot_ind_1] != '0':
							tot_tot_1 += 1	
						if b[het_1_ind_2] != '0':
							het_1_tot_2 += 1	
						if b[het_2_ind_2] != '0':
							het_2_tot_2 += 1	
						if b[hom_ind_2] != '0':
							hom_tot_2 += 1	
						if b[tot_ind_2] != '0':
							tot_tot_2 += 1

						if b[hom_ind_1] != '0' and b[hom_ind_2] != '0':
							hom_1_hom_2 += 1

						if b[het_1_ind_1] != '0' and b[het_2_ind_2] != '0' and b[het_2_ind_1] == '0' and b[het_1_ind_2] == '0':
							het_1_het_2 += 1
						elif b[het_2_ind_1] != '0' and b[het_1_ind_2] != '0' and b[het_1_ind_1] == '0' and b[het_2_ind_2] == '0':
							het_1_het_2 += 1

						if b[het_1_ind_1] != '0' and b[hom_ind_2] != '0':
							het_1_hom_2 += 1
						elif b[het_2_ind_1] != '0' and b[hom_ind_2] != '0':
							het_1_hom_2 += 1

						if b[hom_ind_1] != '0' and b[het_1_ind_2] != '0':
							hom_1_het_2 += 1
						elif b[hom_ind_1] != '0' and b[het_2_ind_2] != '0':
							hom_1_het_2 += 1

						if b[tot_ind_1] != '0' and b[tot_ind_2] != '0':
							tot_1_tot_2 += 1

						#if b[tot_ind_1] != '0' and b[het_ind_2] != '0':
						#	tot_1_het_2 += 1
						#if b[tot_ind_1] != '0' and b[hom_ind_2] != '0':
						#	tot_1_hom_2 += 1
						#if b[het_ind_1] != '0' and b[tot_ind_2] != '0':
						#	het_1_tot_2 += 1
						#if b[hom_ind_1] != '0' and b[tot_ind_2] != '0':
						#	hom_1_tot_2 += 1

						if b[hom_ind_1] != '0' and b[hom_ind_2] != '0':
							hom_1_2 += 1
						if b[tot_ind_1] != '0' and b[tot_ind_2] != '0':
							hom_het_1_2 += 1

						#if '*0101' in all_1 and '*0101' not in all_2:
						#	if b[het_ind_1] != '0' and b[tot_ind_2] != '0' and b[hom_ind_1] == '0':
						#		hom_het_mix_1_2 += 1	
						#if '*0101' in all_2 and '*0101' not in all_1:
						#	if b[het_ind_2] != '0' and b[tot_ind_1] != '0' and b[hom_ind_2] == '0':
						#		hom_het_mix_1_2 += 1	
						#if '*0101' in all_2 and '*0101' in all_1:
						#	if b[het_ind_2] != '0' and b[het_ind_1] != '0' and b[hom_ind_2] == '0' and b[hom_ind_1] == '0':
						#		hom_het_mix_1_2 += 1	


					hom_per = float(hom_1_2) / float(len(raw_to_use))			
					hom_het_per = float(hom_het_1_2) / float(len(raw_to_use))			
					hom_het_mix_per = float(hom_het_mix_1_2) / float(len(raw_to_use))			

					hom_1_hom_2_per = float(hom_1_hom_2) / float(len(raw_to_use))
					hom_1_het_2_per = float(hom_1_het_2) / float(len(raw_to_use))
					het_1_hom_2_per = float(het_1_hom_2) / float(len(raw_to_use))
					het_1_het_2_per = float(het_1_het_2) / float(len(raw_to_use))

					tot_1_tot_2_per = float(tot_1_tot_2) / float(len(raw_to_use))
					tot_1_het_2_per = float(tot_1_het_2) / float(len(raw_to_use))
					tot_1_hom_2_per = float(tot_1_hom_2) / float(len(raw_to_use))
					het_1_tot_2_per = float(het_1_tot_2) / float(len(raw_to_use))
					hom_1_tot_2_per = float(hom_1_tot_2) / float(len(raw_to_use))


					c_11_11_per = float(c_11_11) / float(len(raw_to_use))
					c_00_00_per = float(c_00_00) / float(len(raw_to_use))
					c_01_10_per = float(c_01_10) / float(len(raw_to_use))
					c_10_01_per = float(c_10_01) / float(len(raw_to_use))
					c_11_01_per = float(c_11_01) / float(len(raw_to_use))
					c_11_10_per = float(c_11_10) / float(len(raw_to_use))
					c_01_11_per = float(c_01_11) / float(len(raw_to_use))
					c_10_11_per = float(c_10_11) / float(len(raw_to_use))




					out.append(['' + all_1 + '-' + all_2])
					
					out.append(['    Individuals with ' + all_1 + '(1|1) + ' + all_2 + '(1|1): ' + str(c_11_11) + ' (' + str(c_11_11_per) + ')'])	
					out.append(['    Individuals with ' + all_1 + '(0|0) + ' + all_2 + '(0|0): ' + str(c_00_00) + ' (' + str(c_00_00_per) + ')'])	
					out.append(['    Individuals with ' + all_1 + '(1|0) + ' + all_2 + '(0|1): ' + str(c_10_01) + ' (' + str(c_10_01_per) + ')'])	
					out.append(['    Individuals with ' + all_1 + '(0|1) + ' + all_2 + '(1|0): ' + str(c_01_10) + ' (' + str(c_01_10_per) + ')'])	
					out.append(['    Individuals with ' + all_1 + '(1|1) + ' + all_2 + '(0|1): ' + str(c_11_01) + ' (' + str(c_11_01_per) + ')'])	
					out.append(['    Individuals with ' + all_1 + '(1|1) + ' + all_2 + '(1|0): ' + str(c_11_10) + ' (' + str(c_11_10_per) + ')'])	
					out.append(['    Individuals with ' + all_1 + '(0|1) + ' + all_2 + '(1|1): ' + str(c_01_11) + ' (' + str(c_01_11_per) + ')'])	
					out.append(['    Individuals with ' + all_1 + '(1|0) + ' + all_2 + '(1|1): ' + str(c_10_11) + ' (' + str(c_10_11_per) + ')'])	



					#out.append(['    Individuals with ' + all_1 + '(1|1) + ' + all_2 + '(1|1): ' + str(hom_1_hom_2) + ' (' + str(hom_1_hom_2_per) + ')'])
					#out.append(['    Individuals with ' + all_1 + '(homo) + ' + all_2 + '(hetero): ' + str(hom_1_het_2) + ' (' + str(hom_1_het_2_per) + ')'])
					#out.append(['    Individuals with ' + all_1 + '(hetero) + ' + all_2 + '(homo): ' + str(het_1_hom_2) + ' (' + str(het_1_hom_2_per) + ')'])
					#out.append(['    Individuals with ' + all_1 + '(hetero) + ' + all_2 + '(hetero): ' + str(het_1_het_2) + ' (' + str(het_1_het_2_per) + ')'])
					#out.append(['    Individuals with ' + all_1 + '(homo/hetero) + ' + all_2 + '(homo/hetero): ' + str(tot_1_tot_2) + ' (' + str(tot_1_tot_2_per) + ')'])
					#out.append(['    Individuals with ' + all_1 + '(homo/hetero) + ' + all_2 + '(homo): ' + str(tot_1_hom_2) + ' (' + str(tot_1_hom_2_per) + ')'])
					#out.append(['    Individuals with ' + all_1 + '(homo/hetero) + ' + all_2 + '(hetero): ' + str(tot_1_het_2) + ' (' + str(tot_1_het_2_per) + ')'])
					#out.append(['    Individuals with ' + all_1 + '(homo) + ' + all_2 + '(homo/hetero): ' + str(hom_1_tot_2) + ' (' + str(hom_1_tot_2_per) + ')'])
					#out.append(['    Individuals with ' + all_1 + '(hetero) + ' + all_2 + '(homo/hetero): ' + str(het_1_tot_2) + ' (' + str(het_1_tot_2_per) + ')'])


					#out.append(['   Individuals with ' + all_1 + '(homo.) + ' + all_2 + '(homo.): ' + str(hom_1_2) + ' (' + str(hom_per) + ')'])
					#out.append(['   Individuals with ' + all_1 + '(homo/hetero.) + ' + all_2 + '(homo/hetero.): ' + str(hom_het_1_2) + ' (' + str(hom_het_per) + ')'])

					#if '*0101' in all_1 and '*0101' not in all_2:
					#	out.append(['   Individuals with ' + all_1 + '(hetero.) + ' + all_2 + '(homo/hetero.): ' + str(hom_het_mix_1_2) + ' (' + str(hom_het_mix_per) + ')'])
					#if '*0101' in all_2 and '*0101' not in all_1:
					#	out.append(['   Individuals with ' + all_1 + '(homo/hetero.) + ' + all_2 + '(hetero.): ' + str(hom_het_mix_1_2) + ' (' + str(hom_het_mix_per) + ')'])
					#if '*0101' in all_2 and '*0101' in all_1:
					#	out.append(['   Individuals with ' + all_1 + '(hetero.) + ' + all_2 + '(hetero.): ' + str(hom_het_mix_1_2) + ' (' + str(hom_het_mix_per) + ')'])

					out.append([''])


		out.append([''])
		out.append(['#################################'])
		out.append([''])


	out_file = freq_file.split('.txt')[0] + '_allele_freq_totals.txt'
	savefile(out_file, out, [])


def get_pos(rs_list):
	file_1 = 'DOA_exac_ENSG00000204252_2019_01_22_13_18_34_man.txt'
	cols = ['Chrom', 'Position', 'RSID', 'Reference', 'Alternate', 'Consequence', 'Annotation']
	raw1, headers1, inds1 = openfile(file_1, cols)

	file_2 = 'DOB_exac_ENSG00000241106_2019_01_22_13_18_14_man.txt'
	cols = ['Chrom', 'Position', 'RSID', 'Reference', 'Alternate', 'Consequence', 'Annotation']
	raw2, headers2, inds2 = openfile(file_2, cols)

	out = []
	headers_out = ['Name', 'RSID', 'Chrom', 'Position', 'Change', 'Consequence', 'Annotation']

	raw_both = []
	for a in raw1:
		entry = ['']*len(a)
		for i in range(0,len(a)):
			entry[i] = a[i]
		raw_both.append(entry)	

	for a in raw2:
		entry = ['']*len(a)
		for i in range(0,len(a)):
			entry[i] = a[i]
		raw_both.append(entry)	

	allele_list = []

	for a in rs_list:
		entry = ['']*len(headers_out)
		name = a[0]
		rs = a[1]

		chrom = ''
		pos = ''
		change = ''
		cons = ''
		annot = ''

		for b in raw_both:
			if b[inds1[2]] == rs:
				chrom = b[inds1[0]]
				pos = b[inds1[1]]
				change = b[inds1[3]] + '>' + b[inds1[4]]
				cons = b[inds1[5]]
				annot = b[inds1[6]]

		entry[0] = name
		entry[1] = rs
		entry[2] = chrom
		entry[3] = pos
		entry[4] = change
		entry[5] = cons
		entry[6] = annot

		out.append(entry)
		print entry

		allele_ent = [name, 'chr'+chrom, [rs], [pos], [change], [cons]]
		allele_list.append(allele_ent)

	out_file = '0_DO_rs_list_pos.txt'
	savefile(out_file, out, headers_out)


	return allele_list


def extract_af_2(allele_list, vcf_in, file_out_2):
	# Header = 'INFO'

	vcf_raw = openfile(vcf_in, [])

	count = 0
	row_ind = 0
	for a in vcf_raw:
		for b in a:
			if len(b) > 0:
				if b[0] == '#CHROM':
					row_ind = count
			count += 1

	top_info = vcf_raw[0][:row_ind]
	headers = vcf_raw[0][row_ind]
	vcf_data = vcf_raw[0][row_ind+1:]

	pos_ind = 1
	id_ind = 2
	ref_ind = 3
	alt_ind = 4
	info_ind = 7
	start_ind = 9


	seen = []
	for a in allele_list:
		for snp in a[2]:
			found = False
			for b in vcf_data:
				if snp == b[id_ind]:
					found = True
			if found == False and snp not in seen:
				print 'Missing SNP:', snp, '('+a[0]+')'
				seen.append(snp)


	# get population freqs for multi-SNP alleles of DR
	# each individual has diploid readout: e.g. 0|1
	# 0 = REF (not SNP)
	# 1 = ALT (SNP)
	# for DR, consider both homozygous (1|1) and hetero (0|1 and 1|0)
	# if an individual has either case, they "have" the SNP
	# then group all individuals that have all appropriate SNPs
	# the allele freq can then be calculated

	all_ins, h, i = openfile('individuals_ALL.txt', [])
	afr_ins, h, i = openfile('individuals_AFR.txt', [])
	amr_ins, h, i = openfile('individuals_AMR.txt', [])
	eur_ins, h, i = openfile('individuals_EUR.txt', [])
	eas_ins, h, i = openfile('individuals_EAS.txt', [])
	sas_ins, h, i = openfile('individuals_SAS.txt', [])

	indiv_heads = ['Individual', 'Population']
	for a in allele_list:
		indiv_heads.append(a[0] + ' (1|1)')
		indiv_heads.append(a[0] + ' (1|0)')
		indiv_heads.append(a[0] + ' (0|1)')
		indiv_heads.append(a[0] + ' (0|0)')
		#indiv_heads.append(a[0] + ' total')

	individs = []
	for i in range(start_ind,len(headers)):
		entry = ['0']*len(indiv_heads)
		entry[0] = headers[i]

		pop = ''		
		if [entry[0]] in afr_ins:
			pop = 'AFR'
		if [entry[0]] in amr_ins:
			pop = 'AMR'
		if [entry[0]] in eur_ins:
			pop = 'EUR'
		if [entry[0]] in eas_ins:
			pop = 'EAS'
		if [entry[0]] in sas_ins:
			pop = 'SAS'
		entry[1] = pop

		start_ind = 2
		for x in allele_list:
			snp_pass = False
			#homo = 0
			#hetero_1 = 0
			#hetero_2 = 0
			c_11 = 0
			c_10 = 0
			c_01 = 0
			c_00 = 0

			snps = x[2]
			for snp in snps:
				for y in vcf_data:
					if y[2] == snp:
						snp_pass = True
						#print snp, y[i]
						#if y[i] != '0|0':
						#	snp_pass += 1
						if y[i] == '1|1':
							#homo += 1
							c_11 += 1
						elif y[i] == '1|0':
							#hetero_1 += 1
							c_10 += 1
						elif y[i] == '0|1':
							#hetero_2 += 1
							c_01 += 1
						elif y[i] == '0|0':
							#hetero_2 += 1
							c_00 += 1

			if c_11 == len(snps):
				entry[start_ind] = '1'
			elif c_10 == len(snps):
				entry[start_ind+1] = '1'
			elif c_01 == len(snps):
				entry[start_ind+2] = '1'
			elif c_00 == len(snps):
				entry[start_ind+3] = '1'
			if snp_pass == False:
				entry[start_ind+3] = '1'

			start_ind += 4

		individs.append(entry)

	savefile(file_out_2, individs, indiv_heads)






def all_neg_allele(freq_file, blocks):
	cols = ['Population']
	raw, headers, inds = openfile(freq_file, cols)

	for i in range(0,len(blocks)):
		neg_block = ''
		block = blocks[i]
		#neg_block = block[0][0].split('*')[0] + '*0101'
		#headers.append(neg_block + ' hetero_1')
		#headers.append(neg_block + ' hetero_2')
		for a in block:
			if a[1] == 'neg':
				neg_block = a[0]

		headers.append(neg_block + ' (1|1)')
		headers.append(neg_block + ' (1|0)')
		headers.append(neg_block + ' (0|1)')
		headers.append(neg_block + ' (0|0)')


	out = []

	for a in raw:
		entry = ['0']*len(headers)
		for i in range(0,len(a)):
			entry[i] = a[i]


		for i in range(0,len(blocks)):
			for b in blocks[i]:
				if b[1] == 'neg':
					neg_allele = b[0]
			
			num_10 = 0
			num_01 = 0
			num_11 = 0
			num_00 = 0
			total_num = 0

			for b in blocks[i]:
				if b[1] != 'neg':
					allele = b[0]
					ind_11 = ''
					ind_10 = ''
					ind_01 = ''
					ind_00 = ''

					for j in range(0,len(headers)):
						if allele + ' (1|1)' == headers[j]:
							ind_11 = j
						if allele + ' (1|0)' == headers[j]:
							ind_10 = j
						if allele + ' (0|1)' == headers[j]:
							ind_01 = j
						if allele + ' (0|0)' == headers[j]:
							ind_00 = j
						if neg_allele + ' (1|1)' == headers[j]:
							neg_ind_11 = j
						if neg_allele + ' (1|0)' == headers[j]:
							neg_ind_10 = j
						if neg_allele + ' (0|1)' == headers[j]:
							neg_ind_01 = j
						if neg_allele + ' (0|0)' == headers[j]:
							neg_ind_00 = j

					if a[ind_11] != '0':
						num_11 += 1
					if a[ind_10] != '0':
						num_10 += 1
					if a[ind_01] != '0':
						num_01 += 1
					if a[ind_00] != '0':
						num_00 += 1
					total_num += 1

			# if all other SNPs in block are 0|0, neg gets 1|1
			if num_00 == total_num:
				entry[neg_ind_11] = '1'

			# if other SNPs have 1|0 and no 0|1s, neg gets 0|1
			if num_10 > 0 and num_01 == 0:
				entry[neg_ind_01] = '1'

			# if other SNPs have 0|1 and no 1|0s, neg gets 1|0
			if num_01 > 0 and num_10 == 0:
				entry[neg_ind_10] = '1'
			
			# if other SNPs have 1|1 or 1|0+0|1, neg gets 0|0
			if num_11 > 0:
				entry[neg_ind_00] = '1'
			elif num_01 > 0 and num_10 > 0:
				entry[neg_ind_00] = '1'

		out.append(entry)

	out_file = freq_file.split('.txt')[0] + '_neg.txt'
	savefile(out_file, out, headers)



def get_inds(headers, allele):
	ind_00 = ''
	ind_01 = ''
	ind_10 = ''
	ind_11 = ''

	for j in range(0,len(headers)):
		if allele + ' (1|1)' == headers[j]:
			ind_11 = j
		if allele + ' (0|1)' == headers[j]:
			ind_01 = j
		if allele + ' (1|0)' == headers[j]:
			ind_10 = j
		if allele + ' (0|0)' == headers[j]:
			ind_00 = j

	return ind_00, ind_10, ind_01, ind_11


def allele_nums_2(freq_file, blocks):
	cols = ['Population']
	raw, headers, inds = openfile(freq_file, cols)

	unique_pops = []
	for a in raw:
		unique_pops.append(a[inds[0]])
	unique_pops = list(set(unique_pops))
	unique_pops = sorted(unique_pops)
	unique_pops.insert(0,'ALL')

	out = []
	out.append([''])
	letters = [i for i in string.ascii_uppercase]


	freq_all_out = []
	freq_hom_out = []
	freq_heads = ['Haplotype']
	for pop in unique_pops:
		freq_heads.append(pop)

	freq_temp = []
	unique_combs = []

	for pop in unique_pops:
		out.append(['Population: ' + pop])
		out.append([''])

		if pop != 'ALL':
			raw_to_use = []
			for a in raw:
				if a[inds[0]] == pop:
					raw_to_use.append(a)
		else:
			raw_to_use = raw

		for i in range(0,len(blocks)):
			block_let = letters[i]
			block = blocks[i]

			out.append(['Block ' + block_let])

			block_alls = []
			for a in block:
				block_alls.append(a[0])
			block_alls = sorted(block_alls)

			for allele in block_alls:
				ind_00, ind_10, ind_01, ind_11 = get_inds(headers, allele)

				tot_00 = 0
				tot_10 = 0
				tot_01 = 0
				tot_11 = 0
				for b in raw_to_use:
					if b[ind_00] != '0':
						tot_00 += 1
					if b[ind_10] != '0':
						tot_10 += 1
					if b[ind_01] != '0':
						tot_01 += 1
					if b[ind_11] != '0':
						tot_11 += 1
				per_00 = float(tot_00) / float(len(raw_to_use))
				per_10 = float(tot_10) / float(len(raw_to_use))
				per_01 = float(tot_01) / float(len(raw_to_use))
				per_11 = float(tot_11) / float(len(raw_to_use))

				out.append(['  ' + allele])
				out.append(['  Individuals with ' + allele + '(1|1): ' + str(tot_11) + ' (' + str(per_11) + ')'])	
				out.append(['  Individuals with ' + allele + '(1|0): ' + str(tot_10) + ' (' + str(per_10) + ')'])	
				out.append(['  Individuals with ' + allele + '(0|1): ' + str(tot_01) + ' (' + str(per_01) + ')'])	
				out.append(['  Individuals with ' + allele + '(0|0): ' + str(tot_00) + ' (' + str(per_00) + ')'])	
				
				tot = tot_11 + tot_10 + tot_01 + tot_00
				tot_per = float(tot) / float(len(raw_to_use))
				out.append(['  Total: ' + str(tot) + ' (' + str(tot_per) + ')'])

				out.append([''])

			out.append(['Block Combinations'])
			out.append([''])

			block_combs = []
			for h in block_alls:
				for j in block_alls:
					comb = [h, j]
					comb = sorted(comb)
					if comb not in block_combs:
						block_combs.append(comb)
			for comb in block_combs:
				allele_1 = comb[0]
				allele_2 = comb[1]
				
				ind_00_1, ind_10_1, ind_01_1, ind_11_1 = get_inds(headers, allele_1)
				ind_00_2, ind_10_2, ind_01_2, ind_11_2 = get_inds(headers, allele_2)

				c_00 = 0
				c_10 = 0
				c_01 = 0
				c_11 = 0

				for b in raw_to_use:
					left_1 = False
					right_1 = False

					if b[ind_10_1] != '0' or b[ind_11_1] != '0':
						left_1 = True
					if b[ind_01_2] != '0' or b[ind_11_2] != '0':
						right_1 = True

					if left_1 == True and right_1 == True:
						c_11 += 1
					if left_1 == True and right_1 == False:
						c_10 += 1
					if left_1 == False and right_1 == True:
						c_01 += 1
					if left_1 == False and right_1 == False:
						c_00 += 1

				per_00 = float(c_00) / float(len(raw_to_use))
				per_10 = float(c_10) / float(len(raw_to_use))
				per_01 = float(c_01) / float(len(raw_to_use))
				per_11 = float(c_11) / float(len(raw_to_use))

				out.append(['  ' + allele_1 + '/' + allele_2])
				out.append(['  Individuals with (' + allele_1 + '|' + allele_2 + ') (1|1): ' + str(c_11) + ' (' + str(per_11) + ')'])	
				out.append(['  Individuals with (' + allele_1 + '|' + allele_2 + ') (1|0): ' + str(c_10) + ' (' + str(per_10) + ')'])	
				out.append(['  Individuals with (' + allele_1 + '|' + allele_2 + ') (0|1): ' + str(c_01) + ' (' + str(per_01) + ')'])	
				out.append(['  Individuals with (' + allele_1 + '|' + allele_2 + ') (0|0): ' + str(c_00) + ' (' + str(per_00) + ')'])	
				
				tot = c_11 + c_10 + c_01 + c_00
				tot_per = float(tot) / float(len(raw_to_use))
				out.append(['  Total: ' + str(tot) + ' (' + str(tot_per) + ')'])
				
				out.append([''])



		out.append(['Allele Combinations'])

		if len(blocks) == 2:
			combs = []
			for j in blocks[0]:
				for h in blocks[1]:
					comb = [j[0],h[0]]
					comb = sorted(comb)
					if comb not in combs:
						combs.append(comb)
	
			for comb in combs:
				allele_1 = comb[0]
				allele_2 = comb[1]	
				ind_00_1, ind_10_1, ind_01_1, ind_11_1 = get_inds(headers, allele_1)
				ind_00_2, ind_10_2, ind_01_2, ind_11_2 = get_inds(headers, allele_2)


				c_11_11 = 0
				c_11_10 = 0
				c_11_01 = 0
				c_11_00 = 0
				
				c_10_11 = 0
				c_10_10 = 0
				c_10_01 = 0
				c_10_00 = 0
				
				c_01_11 = 0
				c_01_10 = 0
				c_01_01 = 0
				c_01_00 = 0

				c_00_11 = 0
				c_00_10 = 0
				c_00_01 = 0
				c_00_00 = 0

				for b in raw_to_use:
					if b[ind_11_1] != '0' and b[ind_11_2] != '0':
						c_11_11 += 1
					if b[ind_11_1] != '0' and b[ind_10_2] != '0':
						c_11_10 += 1
					if b[ind_11_1] != '0' and b[ind_01_2] != '0':
						c_11_01 += 1
					if b[ind_11_1] != '0' and b[ind_00_2] != '0':
						c_11_00 += 1

					if b[ind_10_1] != '0' and b[ind_11_2] != '0':
						c_10_11 += 1
					if b[ind_10_1] != '0' and b[ind_10_2] != '0':
						c_10_10 += 1
					if b[ind_10_1] != '0' and b[ind_01_2] != '0':
						c_10_01 += 1
					if b[ind_10_1] != '0' and b[ind_00_2] != '0':
						c_10_00 += 1

					if b[ind_01_1] != '0' and b[ind_11_2] != '0':
						c_01_11 += 1
					if b[ind_01_1] != '0' and b[ind_10_2] != '0':
						c_01_10 += 1
					if b[ind_01_1] != '0' and b[ind_01_2] != '0':
						c_01_01 += 1
					if b[ind_01_1] != '0' and b[ind_00_2] != '0':
						c_01_00 += 1

					if b[ind_00_1] != '0' and b[ind_11_2] != '0':
						c_00_11 += 1
					if b[ind_00_1] != '0' and b[ind_10_2] != '0':
						c_00_10 += 1
					if b[ind_00_1] != '0' and b[ind_01_2] != '0':
						c_00_01 += 1
					if b[ind_00_1] != '0' and b[ind_00_2] != '0':
						c_00_00 += 1

				per_11_11 = float(c_11_11) / float(len(raw_to_use))
				per_11_10 = float(c_11_10) / float(len(raw_to_use))
				per_11_01 = float(c_11_01) / float(len(raw_to_use))
				per_11_00 = float(c_11_00) / float(len(raw_to_use))

				per_10_11 = float(c_10_11) / float(len(raw_to_use))
				per_10_10 = float(c_10_10) / float(len(raw_to_use))
				per_10_01 = float(c_10_01) / float(len(raw_to_use))
				per_10_00 = float(c_10_00) / float(len(raw_to_use))
		
				per_01_11 = float(c_01_11) / float(len(raw_to_use))
				per_01_10 = float(c_01_10) / float(len(raw_to_use))
				per_01_01 = float(c_01_01) / float(len(raw_to_use))
				per_01_00 = float(c_01_00) / float(len(raw_to_use))

				per_00_11 = float(c_00_11) / float(len(raw_to_use))
				per_00_10 = float(c_00_10) / float(len(raw_to_use))
				per_00_01 = float(c_00_01) / float(len(raw_to_use))
				per_00_00 = float(c_00_00) / float(len(raw_to_use))

				out.append(['  ' + allele_1 + '/' + allele_2])
				out.append(['  Individuals with ' + allele_1 + '(1|1)' + allele_2 + '(1|1): ' + str(c_11_11) + ' (' + str(per_11_11) + ')'])	
				out.append(['  Individuals with ' + allele_1 + '(1|1)' + allele_2 + '(1|0): ' + str(c_11_10) + ' (' + str(per_11_10) + ')'])	
				out.append(['  Individuals with ' + allele_1 + '(1|1)' + allele_2 + '(0|1): ' + str(c_11_01) + ' (' + str(per_11_01) + ')'])	
				out.append(['  Individuals with ' + allele_1 + '(1|1)' + allele_2 + '(0|0): ' + str(c_11_00) + ' (' + str(per_11_00) + ')'])	

				out.append(['  Individuals with ' + allele_1 + '(1|0)' + allele_2 + '(1|1): ' + str(c_10_11) + ' (' + str(per_10_11) + ')'])	
				out.append(['  Individuals with ' + allele_1 + '(1|0)' + allele_2 + '(1|0): ' + str(c_10_10) + ' (' + str(per_10_10) + ')'])	
				out.append(['  Individuals with ' + allele_1 + '(1|0)' + allele_2 + '(0|1): ' + str(c_10_01) + ' (' + str(per_10_01) + ')'])	
				out.append(['  Individuals with ' + allele_1 + '(1|0)' + allele_2 + '(0|0): ' + str(c_10_00) + ' (' + str(per_10_00) + ')'])	
		
				out.append(['  Individuals with ' + allele_1 + '(0|1)' + allele_2 + '(1|1): ' + str(c_01_11) + ' (' + str(per_01_11) + ')'])	
				out.append(['  Individuals with ' + allele_1 + '(0|1)' + allele_2 + '(1|0): ' + str(c_01_10) + ' (' + str(per_01_10) + ')'])	
				out.append(['  Individuals with ' + allele_1 + '(0|1)' + allele_2 + '(0|1): ' + str(c_01_01) + ' (' + str(per_01_01) + ')'])	
				out.append(['  Individuals with ' + allele_1 + '(0|1)' + allele_2 + '(0|0): ' + str(c_01_00) + ' (' + str(per_01_00) + ')'])	
				
				out.append(['  Individuals with ' + allele_1 + '(0|0)' + allele_2 + '(1|1): ' + str(c_00_11) + ' (' + str(per_00_11) + ')'])	
				out.append(['  Individuals with ' + allele_1 + '(0|0)' + allele_2 + '(1|0): ' + str(c_00_10) + ' (' + str(per_00_10) + ')'])	
				out.append(['  Individuals with ' + allele_1 + '(0|0)' + allele_2 + '(0|1): ' + str(c_00_01) + ' (' + str(per_00_01) + ')'])	
				out.append(['  Individuals with ' + allele_1 + '(0|0)' + allele_2 + '(0|0): ' + str(c_00_00) + ' (' + str(per_00_00) + ')'])	

				tot = 0
				tot += c_11_11
				tot += c_11_10
				tot += c_11_01
				tot += c_11_00
				tot += c_10_11
				tot += c_10_10
				tot += c_10_01
				tot += c_10_00
				tot += c_01_11
				tot += c_01_10
				tot += c_01_01
				tot += c_01_00
				tot += c_00_11
				tot += c_00_10
				tot += c_00_01
				tot += c_00_00
				tot_per = float(tot) / float(len(raw_to_use))
				out.append(['  Total: ' + str(tot) + ' (' + str(tot_per) + ')'])

				out.append([''])

				freq_all = 0
				freq_all += c_11_11
				freq_all += c_11_10
				freq_all += c_11_01
				freq_all += c_10_11
				freq_all += c_10_10
				freq_all += c_01_11
				freq_all += c_01_01
				freq_all_per = float(freq_all) / float(len(raw_to_use))

				freq_hom = 0
				freq_hom += c_11_11

				freq_hom_per = float(freq_hom) / float(len(raw_to_use))

				comb_str = allele_1 + '/' + allele_2
				comb_1 = [allele_1, allele_2]
				comb_1 = sorted(comb)
				if comb_1 not in unique_combs:
					unique_combs.append(comb_1)

				freq_temp.append([pop, comb_str, str(freq_all_per), str(freq_hom_per)])
				

		out.append([''])
		out.append(['###############################'])
		out.append([''])


	freq_all_out = []
	freq_hom_out = []

	unique_combs = sorted(unique_combs)

	for comb in unique_combs:
		comb_str = comb[0] + '/' + comb[1]
		ent_all = ['0.0']*len(freq_heads)
		ent_hom = ['0.0']*len(freq_heads)
		ent_all[0] = comb_str
		ent_hom[0] = comb_str

		pop_ind = 1
		for pop in unique_pops:
			for a in freq_temp:
				if a[0] == pop and a[1] == comb_str:
					ent_all[pop_ind] = a[2]
					ent_hom[pop_ind] = a[3]
			pop_ind += 1

		freq_all_out.append(ent_all)
		freq_hom_out.append(ent_hom)


	savefile(freq_file.split('.txt')[0] + '_stats.txt', out, [])

	out_file = freq_file.split('.txt')[0] + '_freq_all_table.txt'
	savefile(out_file, freq_all_out, freq_heads)
	out_file = freq_file.split('.txt')[0] + '_freq_hom_table.txt'
	savefile(out_file, freq_hom_out, freq_heads)



# 1. generate allele list:

# get positions/chromosome from https://www.ncbi.nlm.nih.gov/projects/SNP/ with rs code
# get position from GRCh37.p13

# 19/01/22: HLA-DO
# 3 alphas, 4 betas
# Look for all combinations

rs_list = [
	# Assigned alleles
	['DOA*0103', 'rs41542323'],
	['DOA*0102', 'rs11575906'],
	['DOA*0104N', 'rs41541116'],
	['DOB*0102', 'rs2071554'],
	['DOB*0105', 'rs11575907'],
	['DOB*0104', 'rs2070121'],
	['DOB*0103', 'rs2621330'],

	['DOA*0105', 'rs10947368'],
	]

allele_list = get_pos(rs_list)


block_1 = [
	['DOA*0101', 'neg'],
	['DOA*0103', 'rs41542323'],
	['DOA*0102', 'rs11575906'],
	['DOA*0104N', 'rs41541116'],
	['DOA*0105', 'rs10947368'],
	]

block_2 = [
	['DOB*0101', 'neg'],
	['DOB*0102', 'rs2071554'],
	['DOB*0105', 'rs11575907'],
	['DOB*0104', 'rs2070121'],
	['DOB*0103', 'rs2621330'],
	]

blocks = [block_1, block_2]


# 2. download raw .vcf data from 1000 Genomes Project:
# specify the range of positions you want (e.g. position 31000000 - 33000000)
# tabix -fh ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 6:31000000-33000000 > 0_raw_31000000_33000000.vcf

# For DO:
# 31000000 - 34000000


# 3. run "gen_vcf" with allele_list and the raw .vcf file
# specify output file (e.g. '1_DR_alleles_raw.vcf')
in_file = '0_raw_31000000_34000000.vcf'
out_file = '1_DO_alleles_raw_small.vcf'
#gen_vcf(allele_list, in_file, out_file)

# extract_af will give you allele frequencies
vcf_in = out_file
out_file = '1_DO_allele_freqs_small.txt'
#extract_af(allele_list, vcf_in)

freq_file = 'DRO_allele_freqs_individuals.txt'
#get_allele_freqs(freq_file, blocks)


file_out = 'DO_indiv_freqs_table.txt'
extract_af_2(allele_list, vcf_in, file_out)
all_neg_allele(file_out, blocks)

freq_file = file_out.split('.txt')[0] + '_neg.txt'
allele_nums_2(freq_file, blocks)


# 4. run vcftools for each population (and 'ALL')
# this creates a file ready for plink
# vcftools --vcf 1_DO_alleles_raw_small.vcf --keep individuals_AFR.txt --out 2_DO_alleles_small_AFR --plink-tped

pops = ['ALL', 'AFR', 'EAS', 'AMR', 'EUR', 'SAS']
#for pop in pops:
#	cmd = 'vcftools --vcf 1_DO_alleles_raw_small.vcf --keep individuals_' + pop + '.txt --out 2_DO_alleles_' + pop + '_small --plink-tped'
#	os.system(cmd)

#	cmd = './plink --tfile 2_DO_alleles_' + pop + '_small --recodeHV --out 3_DO_alleles_' + pop + '_small_forHaploview'
#	os.system(cmd)



# 5. generate pedigree files with plink
# ./plink --tfile 2_DO_alleles_AFR --recodeHV --out 3_DR_alleles_AFR_forHaploview


# 6. run HaploView.jar and open the .ped files
# define blocks (alpha and beta)



