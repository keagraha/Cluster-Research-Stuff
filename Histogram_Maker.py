"""This program will be used to filter and make histograms using the Matcha data.
The user will input where they want cuts to be for: redshift, richness, and
the ratio of core_temp/r500_core_cropped_temp. The program will output multiple
labelled histograms."""

import numpy as np
import matplotlib.pyplot as plt
import pandas
from astropy.table import Table

data = Table.read("y3a2-6.4.22+2_merged_and_filteredforealsies.fits")
pdata = data.copy() 
pdata.remove_columns(['Obsids','serendipitous_dist', '500_kiloparsecs_fractional_source_exposures', 
	'500_kiloparsecs_fractional_background_exposures']) 
	#These columns are multidimensional so they can't be converted to pandas
df = pdata.to_pandas() #Converting table to pandas dataframe for easier use

def histograms(rss = 0.45, rns = 100, trs = 0.9, histos_made = True):
	"""Makes histograms using various splits in redshift, richness, and the ratio.

	Parameters
	----------

	rss: double, in [0, inf) (default = 0.45)
		Redshift Split; where the redshift column is split to differentiate
		between low and high redshift
	rns: double, in [0, inf) (default = 100)
		Richness Split; where the richness column is split to differentiate
		between low and high richness
	trs: double, in [0, inf) (default = 0.9)
		Ratio Split; where the ratio between core_temp/r500_core_cropped_temp
		is split to differentiate between cool core clusters and non cool core
		clusters
	histos_made: boolean (default = True)
		Determines if histograms are made or not, just to check size_string w/o
		having to close a bunch of graphs
	"""

	#Assertions
	assert rss >= 0, "Redshift Split is equal to or greater than 0."
	assert rns >= 0, "Richness Split is equal to or greater than 0"
	assert trs >= 0, "The ratio is nonnegative"

	#Filters
	high_redshift = df['Redshift'] >= rss
	low_redshift = df['Redshift'] < rss
	high_richness = df['lambda'] >= rns
	low_richness = df['lambda'] < rns
	non_cool_core = (df['core_temperature']/df['r500_core_cropped_temperature']) >= 0.9
	cool_core = (df['core_temperature']/df['r500_core_cropped_temperature']) < 0.9

	#Sizes of each sample
	size = len(df)
	high_redshift_size = len(df['Redshift'][high_redshift]) 
	low_redshift_size = size - high_redshift_size
	high_richness_size = len(df['lambda'][high_richness]) 
	low_richness_size = size - high_richness_size
	non_cool_core_size = len(df['Redshift'][non_cool_core]) 
	cool_core_size = size - non_cool_core_size

	#Some Info Based on Splits
	size_string = f"""\
	The number of clusters in the inputted table is {size}.
	------------------------------------------------------
	The number of clusters with high redshift is {high_redshift_size}.
	The number of clusters with low redshift is {low_redshift_size}.
	------------------------------------------------------
	The number of clusters with high richness is {high_richness_size}.
	The number of clusters with low richness is {low_richness_size}.
	------------------------------------------------------
	The number of clusters with cool cores is {cool_core_size}.
	The number of clusters with non-cool cores is {non_cool_core_size}.
	"""
	print(size_string)

	#Histograms
	if(histos_made == True):
		fig1, ax = plt.subplots()
		ax.hist((df['core_temperature']/df['r500_core_cropped_temperature']))
		ax.set_xlabel("Ratio of Core Temp to r500 Core Cropped Temp")
		ax.set_ylabel("Frequency")
		ax.set_title("Histogram of Ratio of Core Temp to r500 Core Cropped Temp")
		plt.show()

		fig2, ay = plt.subplots()
		ay.hist((df['core_temperature'][high_redshift]/df['r500_core_cropped_temperature'][high_redshift]),
			label = 'Redshift > %.2f' %(rss), color = 'red', alpha = 0.25)
		ay.hist((df['core_temperature'][low_redshift]/df['r500_core_cropped_temperature'][low_redshift]),
			label = 'Redshift < %.2f' %(rss), color = 'blue', alpha = 0.25)
		ay.set_xlabel("Ratio of Core Temp to r500 Core Cropped Temp")
		ay.set_ylabel("Frequency")
		ay.set_title("Histogram of Ratio of Core Temp to r500 Core Cropped Temp with Bins of Redshift")
		ay.legend()
		plt.show()

		fig3, az = plt.subplots()
		az.hist((df['core_temperature'][high_richness]/df['r500_core_cropped_temperature'][high_richness]),
			label = 'Richness > %.2f' %(rns), color = 'red', alpha = 0.25)
		az.hist((df['core_temperature'][low_richness]/df['r500_core_cropped_temperature'][low_richness]),
			label = 'Richness < %.2f' %(rns), color = 'blue', alpha = 0.25)
		az.set_xlabel("Ratio of Core Temp to r500 Core Cropped Temp")
		az.set_ylabel("Frequency")
		az.set_title("Histogram of Ratio of Core Temp to r500 Core Cropped Temp with Bins of Richness")
		az.legend()
		plt.show()

		fig4, ax = plt.subplots()
		ax.hist((df['core_temperature'][cool_core]/df['r500_core_cropped_temperature'][cool_core]))
		ax.set_xlabel("Ratio of Core Temp to r500 Core Cropped Temp")
		ax.set_ylabel("Frequency")
		ax.set_title("Histogram of Ratio of Core Temp to r500 Core Cropped Temp Less Than %.2f" %(trs))
		plt.show()

		fig5, ay = plt.subplots()
		ay.hist((df['core_temperature'][cool_core][high_redshift]/df['r500_core_cropped_temperature'][cool_core][high_redshift]),
			label = 'Redshift > %.2f' %(rss), color = 'red', alpha = 0.25)
		ay.hist((df['core_temperature'][cool_core][low_redshift]/df['r500_core_cropped_temperature'][cool_core][low_redshift]),
			label = 'Redshift < %.2f' %(rss), color = 'blue', alpha = 0.25)
		ay.set_xlabel("Ratio of Core Temp to r500 Core Cropped Temp")
		ay.set_ylabel("Frequency")
		ay.set_title("Histogram of Cool Core Clusters with Bins of Redshift")
		ay.legend()
		plt.show()

		fig6, az = plt.subplots()
		az.hist((df['core_temperature'][cool_core][high_richness]/df['r500_core_cropped_temperature'][cool_core][high_richness]),
			label = 'Richness > %.2f' %(rns), color = 'red', alpha = 0.25)
		az.hist((df['core_temperature'][cool_core][low_richness]/df['r500_core_cropped_temperature'][cool_core][low_richness]),
			label = 'Richness < %.2f' %(rns), color = 'blue', alpha = 0.25)
		az.set_xlabel("Ratio of Core Temp to r500 Core Cropped Temp")
		az.set_ylabel("Frequency")
		az.set_title("Histogram of Cool Core Clusters with Bins of Richness")
		az.legend()
		plt.show()

		fig7, ax = plt.subplots()
		ax.hist((df['core_temperature'][non_cool_core]/df['r500_core_cropped_temperature'][non_cool_core]))
		ax.set_xlabel("Ratio of Core Temp to r500 Core Cropped Temp")
		ax.set_ylabel("Frequency")
		ax.set_title("Histogram of Ratio of Core Temp to r500 Core Cropped Temp Greater Than %.2f" %(trs))
		plt.show()

		fig8, ay = plt.subplots()
		ay.hist((df['core_temperature'][non_cool_core][high_redshift]/df['r500_core_cropped_temperature'][non_cool_core][high_redshift]),
			label = 'Redshift > %.2f' %(rss), color = 'red', alpha = 0.25)
		ay.hist((df['core_temperature'][non_cool_core][low_redshift]/df['r500_core_cropped_temperature'][non_cool_core][low_redshift]),
			label = 'Redshift < %.2f' %(rss), color = 'blue', alpha = 0.25)
		ay.set_xlabel("Ratio of Core Temp to r500 Core Cropped Temp")
		ay.set_ylabel("Frequency")
		ay.set_title("Histogram of Non-Cool Core Clusters with Bins of Redshift")
		ay.legend()
		plt.show()

		fig9, az = plt.subplots()
		az.hist((df['core_temperature'][non_cool_core][high_richness]/df['r500_core_cropped_temperature'][non_cool_core][high_richness]),
			label = 'Richness > %.2f' %(rns), color = 'red', alpha = 0.25)
		az.hist((df['core_temperature'][non_cool_core][low_richness]/df['r500_core_cropped_temperature'][non_cool_core][low_richness]),
			label = 'Richness < %.2f' %(rns), color = 'blue', alpha = 0.25)
		az.set_xlabel("Ratio of Core Temp to r500 Core Cropped Temp")
		az.set_ylabel("Frequency")
		az.set_title("Histogram of Non-Cool Core Clusters with Bins of Richness")
		az.legend()
		plt.show()