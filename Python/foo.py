filename = raw_input("what is the file to parse: ")

f1 = open(filename,'r')
f2 = open(filename+".vector.out",'w')

positions = {}
data = {}


for line in f1:
	a = line.rstrip().split('\t')
	b = a[1].split(', ')
	muts = []
	for x in b:
		if x.find('/') != -1: x = x[:-2]
		muts.append(x)

	data[a[0]] = muts
	
	for var in b:
		if var.find('/') != -1:
			var = var[:-2]
		if not positions.has_key(var[1:]):
			positions[int(var[1:])] = var


pos_sorted = positions.keys()
pos_sorted.sort()
head = ""
for x in pos_sorted:
	head = head + positions[x]+'\t'
head.rstrip()

f2.write('Accession'+'\t'+head+'\n')
for acc in data:
	f2.write(acc+'\t')
	for pos in pos_sorted:
		if positions[pos] in data[acc]: f2.write('1')
		else: f2.write('0')
		f2.write('\t')
	f2.write('\n')

f1.close()
f2.close()
