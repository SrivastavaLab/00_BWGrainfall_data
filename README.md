## Purpose of this folder

To supply a single place in which we can organize all the data for the drought experiment. The hope is to concentrate all the effort here, to make everything both easier and less error-prone.

There are two files that do all the work. For simplicity they are not in a folder. They are:

`01_accessing_data.R` : this file downloads all the data straight from Dropbox. It also performs several checks (usually introduced _after_ we found a problem, to make sure that we find it early if it reappears!). Finally, it downloads taxonomic information from the BWG database, and accesses the latest BWG traits. it then combines all this information with the main dataset.

`02_aggregate_functional_groups` : this file aggregates some of the data at the level of the bromeliad: ibutton data and hydrological data. It also summarizes the invertebrate community by totals of functional groups and taxonomy.