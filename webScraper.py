#!/usr/bin/python3

from lxml import html
from concurrent import futures
from itertools import chain
import requests, re, progressbar
import pandas as pd

#############
# Variables #
#############

# URL of the pages to explore
pages 		= ['http://www.thegoodscentscompany.com/allproc-{}.html'.format(i) for i in range(1,13)]
# Xpath to get the <a> tag containing the url
xpath 		= '//td[@class="llstw8"]/a[@onclick]'
# Regex to get the url from the <a> tag
regex 		= r'openMainWindow\([\'"]?([^\'" >]+)[\'"]?\);.*'
# Order of appearance of the properties in the Html page: property names can be different from the webpage
html_order 	= ['Name','InChIKey','SMILES','CAS Number','Molecular Weight','Odor Type','Odor Strength','Odor Description','Taste Description']
# Order wanted for the CSV file
csv_order 	= ['CAS Number','Name','SMILES','InChIKey','Molecular Weight','Odor Type','Odor Strength','Odor Description','Taste Description']
# Xpath to retrieve the data: Property names must match here
xpath1 		= ['//*[@class="radw8" and ./text()="{}: "]/../*[@class="radw7"]/text() | '.format(i) for i in ['Std.InChIKey','SMILES','Molecular Weight']]
xpath_cas 	= '//*[@class="radw8" and ./text()="CAS Number: "]/../*[@class="radw11"]/a/text()'
xpath_name 	= '//title/text()'
xpath2 		= ['//*[@class="radw4" and ./text()="{}: "]/../*[@class="radw11"]/text() | '.format(i) for i in ['Odor Type','Odor Strength','Odor Description','Taste Description']]
# Concatenate the above xpaths
xpath3 		= ''.join(xpath1 + xpath2 + [xpath_cas,' | ', xpath_name])
# Name of the output file
outputFile 	= 'TheGoodScents.csv'

###########
# Functions
###########

def generateProgessBar(array, string):
	'''Generates a progress bar on the terminal'''
	widgets = [string ," - [", progressbar.ETA(format='Remaining:  %(eta)s'), "] ", progressbar.Bar(), " ", progressbar.Percentage()]
	return progressbar.ProgressBar(widgets=widgets, max_value=len(array))

def getLinks(page):
	'''Retrieves each compound's <a> tag, returns 1 list of strings per webpage given'''
	# Parse the webpage as an HTML document, and save it in tree
	scrapped = requests.get(page)
	tree = html.fromstring(scrapped.content)
	# Retrieve the link to each compound using an XPath query.
	# Here, we are looking for the <a> tag contained in a <td class="llstw8">
	# It is possible to use the Inspector tool of Chrome to automatically build complex Xpath queries:
	# //*[@id="alltableList1"]/tbody/tr[5]/td[3]/a
	# Here, looking at the source code of the page, we get the following XPath:
	compound_links = tree.xpath(xpath)
	return compound_links

def regexSearch(item):
	'''Get the URL contained in each <a> tag'''
	# Search in the attribute onclick for the regex. Ex of value of the attribute:
	# openMainWindow('http://www.thegoodscentscompany.com/data/rw1247381.html');return false;
	# [^\'" >]+ will match anything that is not a ' " space or > character
	url = re.search(regex, item.attrib['onclick'])
	if url:
		return url.group(1)
	else:
		return None

def dataFetch(url):
	'''Go to each URL and collect data'''
	scrapped = requests.get(url)
	tree = html.fromstring(scrapped.content)
	content = tree.xpath(xpath3)
	d = {}
	for i in range(len(html_order)):
		try:
			# Only keep the Name, remove other characters produced by HTML
			if i==0:
				name = content[i].replace('\xa0',';').split(';')[0]
				d[html_order[i]] = name
			else:
				# Sometimes, only the name, CAS and Mw are available,
				# which breaks the way we assign properties to their value
				# This is easy to detect as it will put the Mw in the SMILES cell
				if re.match(r'\d', content[2]):
					# If only digits are found in the cell supposed to contained SMILES
					d = {'Name': name, 'CAS Number': content[1], 'Molecular Weight': content[2]}
					break
				else:
					d[html_order[i]] = content[i]
		except IndexError:
			# If the content doesn't exist
			d[html_order[i]] = ''
	return d

def exThreadSubmit(array, function, string):
	'''Submit calculations to multiple threads.
	Applies a function to each element of the array given, and calls the progress bar with a specific description.'''
	# uses a pool of threads to execute calls asynchronously
	with futures.ThreadPoolExecutor() as executor:
		jobs = []
		results = []
		# Submit jobs
		for item in array:
			job = executor.submit(function, item)
			jobs.append(job)
		pbar = generateProgessBar(jobs, string)
		# Get results (with progress bar) as they are completed
		for job in pbar(futures.as_completed(jobs)):
			# If result is not None
			if job.result():
				results.append(job.result())
	return results

########
# Main #
########

# Retrieve each compound's <a> tag, returns 1 list of tags per webpage given
list_of_items_lists = exThreadSubmit(pages, getLinks, 'Scraping HTML pages')
# Convert to a flat list
aTag_list = list(chain.from_iterable(list_of_items_lists))
print("Found",len(aTag_list), "compounds")
# Get the URL contained in each <a> tag
url_list = exThreadSubmit(aTag_list, regexSearch, 'Fetching URL')
print("Fetched", len(url_list), "URL")
# Go to each URL and collect data
data = exThreadSubmit(url_list, dataFetch, 'Collecting data')
# Create a dataframe with the data
df = pd.DataFrame(data)
print("Collected data on", len(df), "compounds")
# Reorganize the columns
df = df[csv_order]
# Write dataframe to CSV file
df.to_csv(outputFile, sep=';', index=False)