SPECIFIC TO THE VESSEL CLASS
=================================

Many of these specificities are
quarter specific, meaning that 
data can be realod at the beginnig of each quarter
during the simulation
when a discrete event (beginnig oof a quarter)
is detected.

******************************
*** STL multimap containers
******************************


vesselsspe_harbours_quarterXX
index of pt_graph (node) visited by each vessel

vesselsspe_freq_harbours_quarterXX
frequency of numbers of visits by each vessel


vesselsspe_fgrounds_quarterXX
index of pt_graph (node) visited by each vessel when fishing

vesselsspe_freq_fgrounds_quarterXX
frequency of numbers of visits by each vessel on this fishing area

vesselsspe_percent_tacs_per_pop_semesterXX
Vessel ID | prop | pop

vesselsspe_betas_semesterXX
Vessel ID | prop | pop





******************************
*** STL multimap containers PER VESSEL
******************************

XXXX_gshape_cpue_per_stk_on_nodes_quarterXX
pt_graph | pop code | param value


XXXX_gscale_cpue_per_stk_on_nodes_quarterXX
pt_graph | pop code | param value


XXXX_cpue_per_stk_on_nodes_quarterXX
pt_graph | pop code | value

XXXX_freq_possible_metiers_quarterXX
pt_graph | metier code | prop


******************************
*** STL map containers
******************************

XXXX_possible_metiers_quarterXX
pt_graph | metier code




******************************
*** Other format
******************************

vesselsspe_features_quarterXX
   c('VE_REF', 'speed', 'fuelconsrate', 'length', 'KW', 'carrying_capacity_model_nls', 
   'tank_capacity_model_nls', 'nb_pings_per_trip', 'shape', 'scale', 'av.trip.duration')

vesselsspe_features_quarterXX_subset
as a trick to subset for some vessels...