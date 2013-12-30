#! /usr/bin/env python

import re
import random
import sys

InputFile = sys.argv[1]
InFile = open(InputFile)
Lines = InFile.readlines()

for Line in Lines:
	Random_number1 = random.randint(10000,99999)
	Random_number2 = random.randint(10000,99999)
	MyRe = r"(\>gi)\|12345\|(emb)\|A12345\|( \w.+)"
	MySub = r'\1|' + str(Random_number1) + r'|\2|' + str(Random_number2) + r'|\3'
	Result = re.sub(MyRe, MySub, Line)
	print Result,
InFile.close()
