#!/usr/bin/python3 -u

import builtins
import sys
import argparse
from math import log
from numpy import arange
from operator import itemgetter

def print(*objects, **kwargs):
	kwargs['flush'] = True
	return builtins.print(*objects, **kwargs)

parser = argparse.ArgumentParser(
	prog="scNet_integration_4_summation.py",
	description="""Make concensus table with given jointable & weight & LLS cutoff.
(Python-ported script of concensus_by*.pl)"""
)
parser.add_argument(
	"jointable",
	type=str,
	metavar="JOINTABLE",
	help="Directory of joined table file."
)
parser.add_argument(
	"concensus",
	type=str,
	metavar="CONCENSUS",
	help="Directory of output concensus table file."
)
parser.add_argument(
	"method",
	type=str,
	metavar="METHOD",
	choices=["max", "wsum", "naivesum"],
	help="Summation method. weight = 1: naive Bayesian summation. weight = positive infinity: 'max'. else: weighted summation"
)
parser.add_argument(
	"cutoff",
	type=float,
	metavar="CUTOFF",
	help="LLS cutoff value."
)
parser.add_argument(
	"-w", "--weight",
	type=int,
	metavar="WEIGHT",
	default=0,
	help="Weight value for weighted summation method. Default: 0 (will raise error)"
)

args = parser.parse_args(sys.argv[1:])

def weighted_sum(w, values):
	"""
	Take weight and sorted list of scores(by descending order),
	Return a single final score.
	"""
	if w == 1: #  naive sum
		return sum(values)
	elif w == "Inf": #  positive infinity
		return values[0] #  equivalent with max(values)
	else:
		wsum = values[0]
		for d in range(len(values[1:])):
			wsum += values[d + 1] / (w * (d + 1))
		return wsum

output = open(args.concensus, "w")

if args.method == "max":
	weight = "Inf"
elif args.method == "naivesum":
	weight = 1
else:
	weight = args.weight

jointable = list()
with open(args.jointable) as J:
  for line in J:
    values = sorted([float(k) for k in line.strip().split("\t")[2:] if k != "NA" and float(k) >= args.cutoff], reverse=True)
    if values:
      jointable.append((line.strip(), weighted_sum(weight, values)))
      
for linecontent, wsum_value in sorted(jointable, key=lambda x:x[1], reverse=True):
  print(linecontent, wsum_value, sep="\t", file=output)
      

output.close()
