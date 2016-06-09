import sys
from pprint import pprint

def main(inputname):
	results = open(inputname).readlines()

	vocab = {}
	for line in results[2:]:
		try:
			vocab[int(line)] += 1 
		except:
			vocab[int(line)] = 1 
	pprint (vocab)

if __name__ == "__main__":
	
	if len(sys.argv) > 2:
		print "Too many arguments"
		sys.exit(1)

	if len(sys.argv) < 2:
		print "Missing arguments"
		sys.exit(1)

	inputname = sys.argv[1]

	main(inputname)