
INPUT DATA SPECIFIC TO THE 'POPULATION' CLASS
============================================

STL multimap containers are used for storing data in 3 dimensions


XX give either the pop index
or the number of the concerned biolsce


for info:
pop_names_baltic_only.data
list the names of the population and coding

a population can be described
with numbers-at-age 
(numbers of individuals in each age class)
or with numbers-at-size group 
(numbers of individuals in each class of body length)


*********************************
***To build STL multimap in c++:
*********************************

init_weight_per_szgroup_biolsceXX
stock | weight | szgroup

init_proprecru_per_szgroup_biolsceXX
stock | prop | szgroup

init_pops_per_szgroup_biolsceXX
stock | numbers | szgroup

init_maturity_per_szgroup_biolsceXX
stock | prop | szgroup

init_M_per_szgroup_biolsceXX
stock | rate | szgroup

init_fecundity_per_szgroup_biolsceXX
stock | prop | szgroup

comcat_per_szgroup_done_by_hand
stock | code | szgroup


*********************************
***To build STL map in c++:
*********************************
avaiXX_betas_se_semesterXX
stock | avai


XXspe_stecf_oth_land_per_month_per_node_semesterXX
pt_graph (node) | landings this pop

percent_landings_from_simulated_vessels

*********************************
***Other formats
*********************************

XXspe_SSB_R_parameters_biolsceXX
2 double numbers for 2 parameters

XXspe_initial_tac
1 double number

XXspe_fbar_amin_amax_ftarget_Fpercent_TACpercent
5 double numbers

XXspe_size_transition_matrix_biolsceXX
a matrix of dim szgroup | szgroup

XXspe_percent_szgroup_per_age_biolsceXX
a matrix of dim szgroup | nbages

XXspe_percent_age_per_szgroup_biolsceXX
a matrix of dim szgroup | nbages

XXctrysspe_relative_stability_semesterXX
few double numbers

multiplier_for_biolsce

*********************************
***Uncommon format 
*********************************

just 5 size groups for optimisation
the_selected_szgroups.dat



