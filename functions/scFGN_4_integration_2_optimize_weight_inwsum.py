#!/usr/bin/python3 -u

import builtins
import sys
from os import execv
import argparse
from math import log
from numpy import arange
from operator import itemgetter
from math import floor

def print(*objects, **kwargs):
	kwargs['flush'] = True
	return builtins.print(*objects, **kwargs)

parser = argparse.ArgumentParser(
	prog="scNet_integration_2_find_best_weight.py",
	description="""Find best weight for network integration.
(Python-ported script of optimize_D_inWS.pl)"""
)
parser.add_argument(
	"gold_standard",
	type=str,
	metavar="GOLDSTANDARD",
	help="Directory of gold standard file."
)
parser.add_argument(
	"jointable",
	type=str,
	metavar="JOINTABLE",
	help="Directory of joined network file."
)
parser.add_argument(
	"report",
	type=str,
	metavar="REPORT",
	help="Directory of output report file."
)
parser.add_argument(
	"-c", "--coverage",
	type=int,
	default=70,
	help="Coverage limit. Default: Maximum coverage(the program will find maximum coverage and restart)"
)
parser.add_argument(
	"-s", "--wstart",
	type=int,
	default=0,
	help="Start point of the weight. Default: 0(naive sum)"
)
parser.add_argument(
	"-e", "--wend",
	type=int,
	default=51,
	help="""End point of the weight. 
			Weight larger than this value will be considered
			as Positive Infinity. Default: 51"""
)
parser.add_argument(
	"-i", "--winterval",
	type=int,
	default=1,
	help="Interval value for the weight search. Default: 1"
)
parser.add_argument(
	"-S", "--llsstart",
	type=float,
	default=1.0,
	help="Start point of the LLS cutoff. Default: 1.0"
)
parser.add_argument(
	"-E", "--llsend",
	type=float,
	default=1.5,
	help="End point of the LLS cutoff. Default: 1.5"
)
parser.add_argument(
	"-I", "--llsinterval",
	type=float,
	default=0.5,
	help="Interval value for the LLS cutoff search. Default: 0.5"
)
parser.add_argument(
	"-r", "--resolution",
	type=int,
	default=1,
	help="Resolution of the cumulative LLS calculation (1~10). Default: 1"
)

def logging(*args, **kwargs):
	"""
	Print given string to stdout and log file.
	Parameters are preserved.
	"""
	print(*args, **kwargs)
	kwargs["file"] = output
	print(*args, **kwargs)

	return None

def weighted_sum(w, values):
	"""
	Take weight and sorted list of scores(by descending order),
	Return a single final score.
	"""
	if w == 0: #  naive sum
		return sum(values)
	elif w == args.wend + 1: #  positive infinity
		return values[0] #  equivalent with max(values)
	else:
		wsum = values[0]
		for d in range(len(values[1:])):
			wsum += values[d + 1] / (w * (d + 1))
		return wsum

args = parser.parse_args(sys.argv[1:])

output = open(args.report, "w")

print("\nReading gold standard...", end="\t")

gs = {tuple(sorted([a,b])) for a,b in (line.strip().split("\t") for line in open(args.gold_standard))}
gs_genes = set([i for k in gs for i in k])

total_annot_pair = (len(gs_genes) * (len(gs_genes) - 1)) // 2
neg_pair = total_annot_pair - len(gs) #  len(gs) = pos_pair
prior_odds = len(gs) / neg_pair

prior_pos_prob = len(gs) / total_annot_pair
prior_neg_prob = neg_pair / total_annot_pair

print("Done.\n")

logging("Total number of annotated genes: %s" % len(gs_genes))
logging("Total number of annotated pairs: %s" % total_annot_pair)
logging("Total number of positive pairs: %s" % len(gs))
logging("Total number of negative pairs: %s" % neg_pair)
logging("Prior odds: %s" % prior_odds)

print("\nReading join table...", end="\t")

jointable = dict()
with open(args.jointable) as J:
	for line in J:
		elements = line.strip().split("\t")
		jointable[tuple(sorted(elements[:2]))] = sorted([float(k) for k in elements[2:] if k != "NA"], reverse=True)
		#  Sort the values at this point & discard 'NA' values to reduce memory usage & running time.
print("Done.")

logging("\nTest start for:")
logging("Weight range %s ~ %s (with interval %s) & positive infinity" % (args.wstart, args.wend, args.winterval))
logging("LLS cutoff range %s ~ %s (with interval %s)" % (args.llsstart, args.llsend, args.llsinterval))
logging("Coverage limit %s (with resolution %s)" % (args.coverage, args.resolution))

lls_range = arange(args.llsstart, args.llsend + args.llsinterval, args.llsinterval) #  range() function doesn't take float values as arguments
for cutoff in [round(k, len(str(args.llsinterval).split(".")[1])) for k in lls_range]:
	#  Trim jointable values first. This will reduce memory usage & running time.
	cut_jointable = dict()
	for pair in jointable:
		cut_value = [k for k in jointable[pair] if k >= cutoff]
		if cut_value:
			cut_jointable[pair] = cut_value
	
	target_coverage = args.coverage
	flag_continue = True
	for weight in range(args.wstart, args.wend + 2, args.winterval): # end point is (weight_end + 1), for positive infinity
		if weight == 0:
			logging("\nTest for weight = %s, LLS cutoff = %s:" % ("0(Naive_sum)", cutoff))
		elif weight == args.wend + 1:
			logging("\nTest for weight = %s, LLS cutoff = %s:" % ("Positive infinity", cutoff))
		else:
			logging("\nTest for weight = %s, LLS cutoff = %s:" % (weight, cutoff))
		sum_result = [(pair,weighted_sum(weight, cut_jointable[pair])) for pair in cut_jointable]
#		sum_result.sort(reverse=True, key=lambda x:x[1]) #  lambda is slower than operator.itemgetter
		sum_result.sort(reverse=True, key=itemgetter(1))
		if len(sum_result) == 0: # Bug fix for 'target coverage 0.0 infinite loop' shlee 20211019
			print(f"Discarding tests for LLS cutoff {cutoff} (No wsum results produced)")
			logging("Sum of cumulative LLS = 0.0")
			continue
		while True:
			unique_genes = set()
			bin_pos_pair = 0
			bin_neg_pair = 0
			posterior_odds = 0.0
			sum_lls = 0.0 # equivalent with 'area' variable in legacy code
			lls = 0.0
			prev_lls = 100.0
			coverage = 0.0
			datapoint = args.resolution
			gsp = 0
			gsn = 0
#			for pair, _ in sum_result:
			for idx, (pair, _) in enumerate(sum_result):
				#  Pair between a GS gene and a non-GS gene will not be considered.
				if pair[0] in gs_genes and pair[1] in gs_genes:
					gsp += 1
					unique_genes.add(pair[0])
					unique_genes.add(pair[1])
					if pair in gs:
						bin_pos_pair += 1
					else:
						bin_neg_pair += 1
				else:
					gsn += 1
				coverage = (len(unique_genes) / len(gs_genes)) * 100
				if coverage < datapoint:
					continue
				else:
					posterior_odds = (bin_pos_pair + 2 * prior_pos_prob) / (bin_neg_pair + 2 * prior_neg_prob)
					lls = log(posterior_odds / prior_odds)
					if lls > prev_lls:
						logging("Abnormal plot shape (%s --> %s), test discarded." % (prev_lls, lls))
						flag_continue = False
						break
					sum_lls += lls
				if datapoint >= target_coverage:
					logging("Sum of cumulative LLS = %s" % sum_lls)
					flag_continue = False
					break
				else:
					datapoint += args.resolution
					continue
			else:
				# logging("The given network doesn't cover %s%% of the gold standard genes with this conditions." % target_coverage)
				# logging("Maximum coverage at current LLS cutoff (%s): %s" % (cutoff, coverage))
				# logging("Reset target coverage to maximum coverage %s and recalculate: " % floor(coverage))
				print("The given network doesn't cover %s%% of the gold standard genes with this conditions." % target_coverage)
				print("Maximum coverage at current LLS cutoff (%s): %s" % (cutoff, coverage))
				print("Reset target coverage to maximum coverage %s and recalculate: " % floor(coverage))
				target_coverage = floor(coverage)
			if not flag_continue:
				break
with open('final_coverage.txt','w') as f:
	f.write("final coverage = "+str(target_coverage))

logging("\nEntire test completed.")
output.close()
