#!/usr/bin/python3 -u

import sys
import os

# failed = os.popen("grep 'Maximum coverage at current LLS cutoff' %s" % sys.argv[1]).read()
# if failed:
# 	from math import floor
# 	print("Optimization failed. Try with coverage %s" % floor(float(failed.split(":")[-1])))
# 	exit()

result = list()
with open(sys.argv[1]) as REPORT:
	current_weight = ""
	current_thr = ""
	for line in REPORT:
		if line.startswith("Test for"):
			elements = line.strip().split(" ")
			current_weight = elements[4][:-1]
			current_thr = elements[8][:-1]
			if current_weight == "Positiv":
				current_weight = "Positive infinity"
				current_thr = elements[9][:-1]
			else:
				current_thr = elements[8][:-1]
			continue
		elif line.startswith("Sum of cumulative"):
			elements = line.strip().split(" ")
			result.append((current_weight, current_thr, float(elements[5])))
		else:
			pass

result.sort(key=lambda x: x[2], reverse=True)
max_auc = result[0][2]

to_print = [[k for k in result if k[2] == max_auc][0]]
for k in to_print:
	print("Best result: Weight %s, Threshold %s, AUC %s" % (k[0], k[1], k[2]))
