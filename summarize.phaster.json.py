#!/usr/bin/env python


import json
import os
import sys
from argparse import ArgumentParser


def parseArgs():
	parser = ArgumentParser(description='Extracts summary data from a '
		'JavaScript Object Notation (JSON) file returned by PHASTER\'s '
		'(PHAge Search Tool Enhanced Release) application programming '
		'interface (API)', add_help=False)
	req = parser.add_argument_group('Required')
	req.add_argument('-i', '--infile', required=True, metavar='FILE',
		help='input JSON file')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-o', '--outfile', required=False, metavar='FILE',
		default=None, help='[stdout]')
	return parser.parse_args()

def main():
	opt = parseArgs()
	ifh = os.path.abspath(os.path.expanduser(opt.infile))

	with open(ifh) as json_file:
		json_data = json.load(json_file)
		if opt.outfile is not None:
			ofh = os.path.abspath(os.path.expanduser(opt.outfile))
			ofh.write(json_data[u'summary'])
		else:
			print(json_data[u'summary']) 

if __name__ == '__main__':
	main()