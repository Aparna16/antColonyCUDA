import sys
from pprint import pprint

def main(inputname):
	results = open(inputname).readlines()

	cities = {}
	# skip first two lines 
	# for each city, verifies its occurrence 
	for line in results[2:]:
		try:
			cities[int(line)] += 1 
		except:
			cities[int(line)] = 1 
	
	for city in cities.keys():
		if cities[city] > 1:
			print "Solucao encontrada nao eh valida!"
			return
	
	print "Solucao valida!"

if __name__ == "__main__":
	
	if len(sys.argv) > 2:
		print "Too many arguments"
		sys.exit(1)

	if len(sys.argv) < 2:
		print "Missing arguments"
		sys.exit(1)

	inputname = sys.argv[1]

	main(inputname)