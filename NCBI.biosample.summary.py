#!/usr/bin/env python3


import gzip
import json
import os
import shutil
import sys
import xml.etree.ElementTree as ET
from argparse import ArgumentParser
from tempfile import mkdtemp

def parseArgs():
	parser = ArgumentParser(description='Creates a JSON database of '
		'BioSample information from NCBI\'s FTP site or from a BioSample '
		'XML file', add_help=False)
	req = parser.add_argument_group('Required')
	req.add_argument('-o', '--output', required=True, metavar='FILE',
		help='output JSON database file')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-s', '--biosample-data', metavar='FILE',
		help='XML local file of biosample data (Schema v2.0 '
		'format as NCBI has specified in https://www.ncbi.nlm.nih.gov/'
		'biosample/docs/submission/validation-service/) rather than directly '
		'fetching from NCBI\'s FTP site; especially useful if FTP is '
		'blocked or creating custom database')
	opt.add_argument('-u', '--url', metavar='STR', type=str, default=None,
		help='web address to a file of biosample data; useful for bypassing '
		'NCBI access (e.g., EMBL-EBI)')
	return parser.parse_args()

def decompress_file(infile, outdir):
	uncompressed_file = os.path.basename(infile).rstrip('.gz')
	outfile = os.path.join(outdir, uncompressed_file)
	with gzip.open(infile, 'rb') as ifh, open(outfile, 'wb') as ofh:
		shutil.copyfileobj(ifh, ofh)
	return outfile

def main():
	opt = parseArgs()

	# Prepare or fetch summary file of biosample data
	if opt.biosample_data:
		summary_file = os.path.realpath(os.path.expanduser(
			opt.biosample_data))
		if summary_file.endswith('.gz'):
			summary_file = decompress_file(summary_file, tmp)
		# TO-DO: verify XML file given conforms to NCBI BioSample (use their
		#        tester?)
	else:
		tfh = os.path.join(tmp, 'biosample_set.xml.gz')
		if opt.url is not None:
			url = opt.url
		else:
			url = 'ftp://ftp.ncbi.nlm.nih.gov/biosample/biosample_set.xml.gz'
		os.system('curl --max-time 300 --silent --fail --show-error '
			'{} -o {}'.format(url, tfh))
		if os.path.exists(tfh) and os.path.getsize(tfh) > 0:
			if tfh.endswith('.gz'):
				summary_file = decompress_file(tfh, tmp)
			else:
				summary_file = tfh
		else:
			sys.stderr.write('ERROR: unable to fetch biosample data file '
				'from NCBI\'s FTP site\n')
			shutil.rmtree(tmp)
			sys.exit(1)

	# Read in XML data
	# NOTE: primary tags from:
	# 	xml_data = ET.iterparse(f, events=('start', 'end'))
	# 	tags = []
	# 	for _, elem in xml_data:
	# 		tags.append(elem.tag)
	# 	tags = list(set(tags))
	group_data = {}
	empty_values = ('null', '', None, 'missing', 'not determined')
	xml_data = ET.iterparse(summary_file, events=('start', 'end'))
	if not opt.biosample_data:
		shutil.rmtree(tmp)
	for event, elem in xml_data:
		# All biosample data are within a single 'BioSampleSet' tag
		if event == 'end' and elem.tag == 'BioSample':
			# d = dict.fromkeys(tags)
			d = {}
			# Extract primary info from 'Id' tags within the 'Ids' tag
			# adds at most 3 key-value pairs to the sample dict
			#    * selects all child elements
			#    . select the current node
			#   // selects all subelements
			biosample_accn = elem.findtext(
				'.//Id[@db="BioSample"][@is_primary="1"]')
			if biosample_accn is None:
				break
			for el in elem.iterfind('.//Id'):
				id_attrs = el.attrib
				if 'db_label' in id_attrs and \
				id_attrs['db_label'] == 'Sample name':
					d['Sample_Name'] = el.text
				elif 'SRA' == id_attrs['db']:
					d['SRA'] = el.text

			# Extract info within the 'Description' tag
			d['Description'] = elem.findtext(
				'.//Description/Comment/Paragraph')
			d['Title'] = elem.findtext('.//Description/Title')
			d['Organism'] = elem.findtext('.//OrganismName')
			organism = elem.find('.//Organism')
			d['Taxonomy_Name'] = organism.get('taxonomy_name')
			d['TaxID'] = int(organism.get('taxonomy_id'))

			# Extract info within the 'Owner' tag
			first = elem.findtext('.//Owner/Contacts/Contact/Name/First')
			middle = elem.findtext('.//Owner/Contacts/Contact/Name/Middle','')
			last = elem.findtext('.//Owner/Contacts/Contact/Name/Last')
			if any((first, middle, last)):
				d['Owner'] = ('{}, {} {}'.format(last, first, middle)).strip()
			contact = elem.find('.//Owner/Contacts/Contact')
			if contact is not None:
				d['Contact'] = contact.attrib['email']

			# Extract info from the 'Model' tag
			d['Model'] = elem.findtext('.//Models/Model')

			# Extract info from the 'Package' tag
			d['Package'] = elem.findtext('.//Package')

			# Extract info from the 'Attributes' tag
			for el in elem.iterfind('.//Attribute'):
				attribute_attrs = el.attrib
				if 'harmonized_name' in attribute_attrs:
					# 415 unique harmonized_name keys exist according to
					# https://www.ncbi.nlm.nih.gov/biosample/docs/attributes/
					# which is why JSON is used as database rather than sqlite
					if el.text not in empty_values:
						d[attribute_attrs['harmonized_name']] = el.text
				elif 'Alias' == attribute_attrs['attribute_name']:
					d['Alias'] = el.text
				elif 'INSDC center name' == attribute_attrs['attribute_name']:
					d['Center'] = el.text
				# NOTE: non-harmonized data below might be problematic; skip
				# elif 'attribute_name' in attribute_attrs:
				# 	d[attribute_attrs['attribute_name']] = el.text

			# Extract info from the 'Links/Link' tags
			for el in elem.iterfind('.//Links/Link'):
				link_attrs = el.attrib
				if link_attrs['type'] == 'entrez' and \
				link_attrs['target'] == 'bioproject':
					d['BioProject'] = int(el.text)
				elif link_attrs['type'] == 'entrez' and \
				link_attrs['target'] == 'pubmed':
					d['PubMed'] = int(el.text)
				elif link_attrs['type'] == 'url':
					d['URL'] = el.text

			# Extract info from the 'Status' tags (uses self-closing tags)
			for el in elem.iterfind('.//Status'):
				status_attrs = el.attrib
				d['Status'] = status_attrs['status']
				d['Status_date'] = status_attrs['when'].split('T')[0]
			d = {k: v for k, v in d.items() if v not in empty_values}
			group_data.update({biosample_accn: d})

	# Output data as JSON
	with open(os.path.realpath(os.path.expanduser(opt.output)), 'w') as o:
		json.dump(group_data, o)

if __name__ == '__main__':
	main()
