import os, sys
usage = '<file with loci of interest><MAGE annotations><MAGE COG file><PFAM pfam_scan tab output><TIGRFAM hmmscan tab output><output name>'

argv = sys.argv[1:]

if len(argv) == 6:
	cogseed = argv.pop(0)
	newproteinmageinfo = argv.pop(0)
	magecog = argv.pop(0)
	pfamscan = argv.pop(0)
	tigrfam = argv.pop(0)
	output = argv.pop(0)
else:
	sys.exit(usage)

mage_anno = {}
seed_dict = {}
output_dict = {}
with open (cogseed,'r') as c:
	for ln in c:
		otherinfo = ln.split('\t')
		label = otherinfo.pop(0)
		seed_dict[label] = otherinfo
with open (magecog,'r') as mc:
	for ln in mc:
		label = ln.split('\t')[0]
		if 'Label' in ln:
			continue
		cog = ln.split('\t')[9]
		if label not in mage_anno:
			mage_anno[label] = {}
		if 'cog' not in mage_anno[label]:
			mage_anno[label]['cog'] = []
		mage_anno[label]['cog'].append(cog)

with open (pfamscan,'r') as pf:
	for ln in pf:
		if '#' in ln or ln == '\n':
			continue
		else:
			loc = ln.split()[0].split('|')[0]
			pfam = ln.split()[6]
			if loc not in mage_anno:
				mage_anno[loc] = {}
			if 'pfam' not in mage_anno[loc]:
				mage_anno[loc]['pfam'] = []
			mage_anno[loc]['pfam'].append(pfam)	


with open (tigrfam,'r') as tf:
	for ln in tf:
		if '#' in ln or ln == '\n':
			continue
		else:
			loc = ln.split()[2].split('|')[0]
			tfam = ln.split()[0]
			if loc not in mage_anno:
				mage_anno[loc] = {}
			if 'tfam' not in mage_anno[loc]:
				mage_anno[loc]['tfam'] = []
			mage_anno[loc]['tfam'].append(tfam)	

with open (newproteinmageinfo,'r') as np:
	linesBefore = list()
	post5 = [12,"NA"]
	for ln in np:
		label = ln.split('\t')[0]
		if label == 'Label':
			continue
		#ln =ln.decode('utf-16','ignore')
		frame = ln.split('\t')[2]
		begin = ln.split('\t')[3]
		end = ln.split('\t')[4]
		length = ln.split('\t')[5]
		product = ln.split('\t')[10]
		meta = [label, frame, begin, end, length, product]
		if label not in mage_anno:
			mage_anno[label] = {}
		mage_anno[label]['meta'] = meta
		mjoin = '\t'.join(mage_anno[label]['meta'])
		if 'pfam' in mage_anno[label]:
			pjoin = ','.join(mage_anno[label]['pfam'])
		else:
			pjoin = 'NA'
		if 'tfam' in mage_anno[label]:
			tjoin = ','.join(mage_anno[label]['tfam'])
		else:
			tjoin = 'NA'
		if 'cog' in mage_anno[label]:
			cjoin = ','.join(mage_anno[label]['cog'])
		else:
			cjoin = 'NA'
		fullmeta = '{0}\t{1}\t{2}\t{3}\n'.format(mjoin,pjoin,tjoin,cjoin)
		linesBefore.append(fullmeta) #ln
		if len(linesBefore) > 6:
			linesBefore.pop(0)
		if post5[0] < 12:
			output_dict[post5[1]][post5[0]] = fullmeta#ln.rstrip() 
			post5[0] = post5[0] + 1
		if label in seed_dict:
			output_dict[label] = {}
			post5 = [6,label]
			for i in range(len(linesBefore)):
				output_dict[label][i] = linesBefore[i] #.rstrip()

with open(output +'_full','w') as out:
	out.write('Locus\tFrame\tBegin\tEnd\tLength\tProduct\tPFAM\tTIGRFAM\tCOG\n')
	for i in mage_anno:
		mjoin = '\t'.join(mage_anno[i]['meta'])
		if 'pfam' in mage_anno[i]:
			pjoin = ','.join(mage_anno[i]['pfam'])
		else:
			pjoin = 'NA'
		if 'tfam' in mage_anno[i]:
			tjoin = ','.join(mage_anno[i]['tfam'])
		else:
			tjoin = 'NA'
		if 'cog' in mage_anno[i]:
			cjoin = ','.join(mage_anno[i]['cog'])
		else:
			cjoin = 'NA'
		out.write('{0}\t{1}\t{2}\t{3}\n'.format(mjoin,pjoin,tjoin,cjoin))

with open(output +'_seed','w') as out:
	out.write('Locus\tFrame\tBegin\tEnd\tLength\tProduct\tPFAM\tTIGRFAM\tCOG\tSeedDistance\n')
	for i in output_dict:
		out.write('\n')
		for x in output_dict[i].keys(): #may need to sort?
			out.write(output_dict[i][x].rstrip()+'\t'+i+'__'+str(x)+'\n')


