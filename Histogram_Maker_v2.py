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

def histograms(rss_1 = 0.29, rss_2 = 0.44 , rns_1 = 91, rns_2 = 114.7 , trs = 0.7, histos_made = True):
	"""Makes histograms using various splits in redshift, richness, and the ratio.

	Parameters
	----------

	rss_1: double, in [0, inf) (default = 0.29)
		Redshift Split 1; where the redshift column is split to differentiate
		between low and medium redshift; the default for this and rss_2 split
		the data set into equally sized groups. 
	rss_2: double, in [0, inf) (default = 0.44)
		Redshift Split 2; where the redshift column is split to differentiate
		between medium and high redshift; the default for this and rss_1 split
		the data set into equally sized groups. 
	rns_1: double, in [0, inf) (default = 91)
		Richness Split 1; where the richness column is split to differentiate
		between low and medium richness; the default for this and rns_2 split
		the data set into equally sized groups. 
	rns_1: double, in [0, inf) (default = 114.7)
		Richness Split 2; where the richness column is split to differentiate
		between medium and high richness, the default for this and rns_1 split
		the data set into equally sized groups. 
	trs: double, in [0, inf) (default = 0.7)
		Ratio Split; where the ratio between core_temp/r500_core_cropped_temp
		is split to differentiate between cool core clusters and non cool core
		clusters
	histos_made: boolean (default = True)
		Determines if histograms are made or not, just to check size_string w/o
		having to close a bunch of graphs
	"""

	#Assertions
	assert rss_1 >= 0, "Redshift Splits must be equal to or greater than 0."
	assert rss_2 >= 0, "Redshift Splits must be equal to or greater than 0."
	assert rns_1 >= 0, "Richness Splits must be equal to or greater than 0."
	assert rns_2 >= 0, "Richness Splits must be equal to or greater than 0."
	assert trs >= 0, "The ratio must be nonnegative."

	#Filters
	high_redshift = df['Redshift'] >= rss_2
	medium_redshift = ((df['Redshift'] < rss_2) & (df['Redshift'] >= rss_1))
	low_redshift = df['Redshift'] < rss_1
	high_richness = df['lambda'] >= rns_2
	medium_richness = ((df['lambda'] < rns_2) & (df['lambda'] >= rns_1))
	low_richness = df['lambda'] < rns_1
	non_cool_core = (df['core_temperature']/df['r500_core_cropped_temperature']) >= trs
	cool_core = (df['core_temperature']/df['r500_core_cropped_temperature']) < trs

	#Sizes of each sample
	size = len(df)
	high_redshift_size = len(df['Redshift'][high_redshift]) 
	medium_redshift_size = len(df['Redshift'][medium_redshift])
	low_redshift_size = size - high_redshift_size - medium_redshift_size
	high_richness_size = len(df['lambda'][high_richness])
	medium_richness_size = len(df['lambda'][medium_richness]) 
	low_richness_size = size - high_richness_size - medium_richness_size
	non_cool_core_size = len(df['Redshift'][non_cool_core]) 
	cool_core_size = size - non_cool_core_size

	#Some Info Based on Splits
	size_string = f"""\
	The number of clusters in the inputted table is {size}.
	------------------------------------------------------
	The number of clusters with high redshift is {high_redshift_size}.
	The number of clusters with medium redshift is {medium_redshift_size}.
	The number of clusters with low redshift is {low_redshift_size}.
	------------------------------------------------------
	The number of clusters with high richness is {high_richness_size}.
	The number of clusters with medium richness is {medium_richness_size}.
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
			label = 'Redshift >= %.2f' %(rss_2), color = 'red', alpha = 0.25)
		ay.hist((df['core_temperature'][medium_redshift]/df['r500_core_cropped_temperature'][medium_redshift]),
			label = '%.2f =< Redshift < %.2f' %(rss_1, rss_2), color = 'green', alpha = 0.25)
		ay.hist((df['core_temperature'][low_redshift]/df['r500_core_cropped_temperature'][low_redshift]),
			label = 'Redshift < %.2f' %(rss_1), color = 'blue', alpha = 0.25)
		ay.set_xlabel("Ratio of Core Temp to r500 Core Cropped Temp")
		ay.set_ylabel("Frequency")
		ay.set_title("Histogram of Ratio of Core Temp to r500 Core Cropped Temp with Bins of Redshift")
		ay.legend()
		plt.show()

		fig3, az = plt.subplots()
		az.hist((df['core_temperature'][high_richness]/df['r500_core_cropped_temperature'][high_richness]),
			label = 'Richness >= %.2f' %(rns_2), color = 'red', alpha = 0.25)
		az.hist((df['core_temperature'][medium_richness]/df['r500_core_cropped_temperature'][medium_richness]),
			label = '%.2f =< Richness < %.2f' %(rns_1, rns_2), color = 'green', alpha = 0.25)
		az.hist((df['core_temperature'][low_richness]/df['r500_core_cropped_temperature'][low_richness]),
			label = 'Richness < %.2f' %(rns_1), color = 'blue', alpha = 0.25)
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
			label = 'Redshift >= %.2f' %(rss_2), color = 'red', alpha = 0.25)
		ay.hist((df['core_temperature'][cool_core][medium_redshift]/df['r500_core_cropped_temperature'][cool_core][medium_redshift]),
			label = '%.2f =< Redshift < %.2f' %(rss_1, rss_2), color = 'green', alpha = 0.25)
		ay.hist((df['core_temperature'][cool_core][low_redshift]/df['r500_core_cropped_temperature'][cool_core][low_redshift]),
			label = 'Redshift < %.2f' %(rss_1), color = 'blue', alpha = 0.25)
		ay.set_xlabel("Ratio of Core Temp to r500 Core Cropped Temp")
		ay.set_ylabel("Frequency")
		ay.set_title("Histogram of Cool Core Clusters with Bins of Redshift")
		ay.legend()
		plt.show()

		fig6, az = plt.subplots()
		az.hist((df['core_temperature'][cool_core][high_richness]/df['r500_core_cropped_temperature'][cool_core][high_richness]),
			label = 'Richness >= %.2f' %(rns_2), color = 'red', alpha = 0.25)
		az.hist((df['core_temperature'][cool_core][medium_richness]/df['r500_core_cropped_temperature'][cool_core][medium_richness]),
			label = '%.2f =< Richness < %.2f' %(rns_1, rns_2), color = 'green', alpha = 0.25)
		az.hist((df['core_temperature'][cool_core][low_richness]/df['r500_core_cropped_temperature'][cool_core][low_richness]),
			label = 'Richness < %.2f' %(rns_1), color = 'blue', alpha = 0.25)
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
			label = 'Redshift >= %.2f' %(rss_2), color = 'red', alpha = 0.25)
		ay.hist((df['core_temperature'][non_cool_core][medium_redshift]/df['r500_core_cropped_temperature'][non_cool_core][medium_redshift]),
			label = '%.2f =< Redshift < %.2f' %(rss_1, rss_2), color = 'green', alpha = 0.25)	
		ay.hist((df['core_temperature'][non_cool_core][low_redshift]/df['r500_core_cropped_temperature'][non_cool_core][low_redshift]),
			label = 'Redshift < %.2f' %(rss_1), color = 'blue', alpha = 0.25)
		ay.set_xlabel("Ratio of Core Temp to r500 Core Cropped Temp")
		ay.set_ylabel("Frequency")
		ay.set_title("Histogram of Non-Cool Core Clusters with Bins of Redshift")
		ay.legend()
		plt.show()

		fig9, az = plt.subplots()
		az.hist((df['core_temperature'][non_cool_core][high_richness]/df['r500_core_cropped_temperature'][non_cool_core][high_richness]),
			label = 'Richness >= %.2f' %(rns_2), color = 'red', alpha = 0.25)
		az.hist((df['core_temperature'][non_cool_core][medium_richness]/df['r500_core_cropped_temperature'][non_cool_core][medium_richness]),
			label = '%.2f =< Richness < %.2f' %(rns_1, rns_2), color = 'green', alpha = 0.25)
		az.hist((df['core_temperature'][non_cool_core][low_richness]/df['r500_core_cropped_temperature'][non_cool_core][low_richness]),
			label = 'Richness < %.2f' %(rns_1), color = 'blue', alpha = 0.25)
		az.set_xlabel("Ratio of Core Temp to r500 Core Cropped Temp")
		az.set_ylabel("Frequency")
		az.set_title("Histogram of Non-Cool Core Clusters with Bins of Richness")
		az.legend()
		plt.show()

