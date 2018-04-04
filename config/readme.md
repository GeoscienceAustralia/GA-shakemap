# Local changes to config files

* Changes made in model.conf are:

	source_network = au
	vs30file = <INSTALL_DIR>/data/vs30/asscm_wii_vs30.grd
	gmpe = stable_continental_au
	max_range = 150

* Updated gmpe_sets.conf to include “stable_continental_au”

* Modules.conf was changed to add:

	SomervilleEtAl2009NonCratonic
	Allen2012

* Updated products.conf to include:

	topography = <INSTALL_DIR>/data/mapping/topography/au_gebco.nc
