using SparseArrays
using LinearAlgebra
using GLMakie
#using AbstractPlotting
#using AbstractPlotting.MakieLayout
#using Plots
using Colors
using BSON
using OffsetArrays 

    struct phreeqc_output_struct # archived form returned by run_phreeqc, interp_sparse 
        driving_parameters::Vector{Float32}
        solute_concentrations::Vector{Float32}
        mineral_saturation_indices::Vector{Float32}
        mineral_omegas::Vector{Float32}
        mineral_reaction_extents::Vector{Float32}
        CO2_uptake_extent::Float32
    end
    struct phreeqc_box_struct # box-specific expanded form 
        driving_parameters 
        pore_saturation 
        solute_concentrations
        mineral_saturation_indices
        mineral_omegas
        mineral_reaction_extents
        mineral_porewater_reaction_rates
        mineral_meters_for_residual_array
        CO2_uptake_extent
    end
    mutable struct box_transect_model_struct 
        timepoints
        parameter_timeseries
        box_surface_elevations
        layer_thickness_timeseries
        mineral_fraction_timeseries
        mineral_reaction_rates_timeseries
        mineral_saturation_state_timeseries
        mineral_summary_timeseries
        component_summary_timeseries 
        tracer_summary_timeseries 
        box_solute_timeseries
        box_soil_CO2_timeseries
        river_solute_timeseries
        regolith_age_timeseries
        erosion_source_timeseries
        CO2_uptake_rate_timeseries 
        #alkalinity_production_bulk_rate_timeseries
        #Be_10_timeseries
    end 

    # Solutes 
    solute_names = ["Ca","Mg","Na","K","Si","Al","Li","Li7","Be9","Be10","Alk","C","pCO2","pH" ]
    n_solutes = length(solute_names)
    Ca_2plus_solute,Mg_2plus_solute,Na_plus_solute,K_plus_solute,Si_4plus_solute,Al_solute,
        Li_solute,Li7_solute,Be9_solute,Be10_solute,Alk_solute,Tot_CO2_solute,pCO2_solute,pH_solute = 
        1,2,3,4,5,6,7,8,9,10,11,12,13,14
    typical_river_solute_concentrations =
    #   ["Ca","Mg","Na","K","Si","Al","Li","Li7","Be9","Be10","Alk","Tot_CO2","pCO2"]
        #[585.,245.,240.,44.,145., NaN, 0.1, -20.,NaN,   NaN,  2000.,    2200.] .* # umol / L Maybeck 2003
        [335.0, 130.0, 223.0, 44.0, 173.0, NaN, 0.1, -20.0, NaN, NaN, 2000.0, 2200.0,3.e6,8.e3] .* # Stallard 1983 
        1.e-3 # mol / m3 except pCO2 which is ppm 

    # Chemical components
    chemical_component_names = ["CaO","MgO","Na2O","K2O","SiO2","Al2O3","CO2","FeO","Fe2O3","H2O"]
    n_chemical_components = length(chemical_component_names)
    CaO_component,MgO_component,
        Na2O_component,K2O_component,SiO2_component,
        Al2O3_component,CO2_component,FeO_component,Fe2O3_component,
        H2O_component = 1,2,3,4,5,6,7,8,9,10
    n_reaction_components = Al2O3_component
    n_solute_aligned_components = SiO2_component
    chemical_component_molwts =
        [56.0, 40.0, 62.0, 94.0, 60.0, 102.0, 44.0, 72.0, 160.0, 18.0]
    solute_component_stoiciometry = fill(1.,n_solutes)
    solute_component_stoiciometry[Na_plus_solute:K_plus_solute] .= 2.
    solute_component_stoiciometry[Al_solute] = 2.

    parameter_names = [["temp","soil CO2"] ; chemical_component_names[1:7]]
    phreeqc_temperatures = [1.]
    phreeqc_temperature_step = 2.
    while phreeqc_temperatures[end] < 32.
        append!(phreeqc_temperatures, 
            phreeqc_temperatures[end] + phreeqc_temperature_step)
    end 
    n_phreeqc_temperatures = size(phreeqc_temperatures)[1]

    #n_points_per_decade = 10
    #omega_n_points_per_decade = 6

    rct_factor_lores = 10. ^ ( 1. / 6. )    # sqrt(sqrt(sqrt(sqrt(sqrt(sqrt(10))))))
    rct_factor_hires = 10. ^ ( 1. / 10. ) 
    phreeqc_pCO2s = [300.] #
    while phreeqc_pCO2s[end] < 300000. 
        append!(phreeqc_pCO2s,phreeqc_pCO2s[end] * rct_factor_hires)
    end
    #,1000.,3000.,10000.,30000.]
    n_phreeqc_pCO2s = size(phreeqc_pCO2s)[1]

    phreeqc_mineral_reaction_constant_grid = [0.,1.e-18]
    while phreeqc_mineral_reaction_constant_grid[end] < 1.e-4
        append!(phreeqc_mineral_reaction_constant_grid, 
            phreeqc_mineral_reaction_constant_grid[end] * rct_factor_hires)
    end 
    n_phreeqc_mineral_reaction_constant_grid = 
        size(phreeqc_mineral_reaction_constant_grid)[1]
    #=
    phreeqc_SiO2_reaction_extents = [1.e-7]
    #for i_rct in 1:48
    while phreeqc_SiO2_reaction_extents[end] < 1.0

        #new_rct = phreeqc_SiO2_reaction_extents[end] * rxn_mult
        append!(phreeqc_SiO2_reaction_extents,
            phreeqc_SiO2_reaction_extents[end] * rct_factor_lores)
    end # to 1000 uM 
    n_phreeqc_SiO2_reaction_extents = size(phreeqc_SiO2_reaction_extents)[1]

    phreeqc_component_silica_ratio_ranges = [ 
        1., # Ca
        2.,  # Mg
        0.1, # Na, 1:1 dissolution with anorthite 
        0.5 / 3., # K, pure k-spar 
        1.,
        0.5] # Al, Ca-spar is the heavy one 
        #0. ] # CO2 in 

        
    phreeqc_silicate_cation_reaction_extents = [1.e-3] # , 0.1, 0.25, 0.5, 0.75, 1.0]# 
    #rct_factor = sqrt(sqrt(sqrt(sqrt(10.0))))
    #i_rct = 1
    while phreeqc_silicate_cation_reaction_extents[end] < 2.0
        append!(
            phreeqc_silicate_cation_reaction_extents,
            phreeqc_silicate_cation_reaction_extents[end] * rct_factor_lores)
    end 
    n_phreeqc_silicate_cation_reaction_extents = size(phreeqc_silicate_cation_reaction_extents)[1]

    phreeqc_CaO_reaction_extents = [[0.] ; phreeqc_silicate_cation_reaction_extents ]
    while phreeqc_CaO_reaction_extents[end] < 1.e3
        new_rct = phreeqc_CaO_reaction_extents[end] * 1.25
        append!(phreeqc_CaO_reaction_extents,new_rct)
    end 
    n_phreeqc_CaO_reaction_extents = size(phreeqc_CaO_reaction_extents)[1]
    
    phreeqc_MgO_reaction_extents = deepcopy(phreeqc_CaO_reaction_extents)
    n_phreeqc_MgO_reaction_extents = size(phreeqc_MgO_reaction_extents)[1]

    phreeqc_Na2O_reaction_extents = phreeqc_CaO_reaction_extents 
    n_phreeqc_Na2O_reaction_extents = size(phreeqc_Na2O_reaction_extents)[1]

    phreeqc_K2O_reaction_extents = [[0.] ; phreeqc_silicate_cation_reaction_extents ]
    n_phreeqc_Na2O_reaction_extents = size(phreeqc_Na2O_reaction_extents)[1]

    phreeqc_Al2O3_reaction_extents = [[0.] ; phreeqc_silicate_cation_reaction_extents ]
    n_phreeqc_Al2O3_reaction_extents = size(phreeqc_Al2O3_reaction_extents)[1]

    phreeqc_Calcite_max_reaction_rates = [[0.] ; phreeqc_SiO2_reaction_extents ]
    n_phreeqc_Calcite_max_reaction_rates = size(phreeqc_Calcite_max_reaction_rates)[1]

    phreeqc_Dolomite_max_reaction_rates = [[0.] ; phreeqc_SiO2_reaction_extents ]
    n_phreeqc_Dolomite_max_reaction_rates = size(phreeqc_Dolomite_max_reaction_rates)[1]


    #=omega_rct_factor = 10. ^ ( 1. / omega_n_points_per_decade )
    phreeqc_carbonate_omega_values = [NaN, 1.e-6]
    while phreeqc_carbonate_omega_values[end] < 1.
        append!(phreeqc_carbonate_omega_values, 
            phreeqc_carbonate_omega_values[end]* omega_rct_factor )
    end 
    n_phreeqc_carbonate_omega_values = size(phreeqc_carbonate_omega_values)[1]
    =#

    #phreeqc_CO2_alk_reaction_ratios = phreeqc_Al2O3_reaction_extents  
    #[0.,1.e-0.05]
    #for i_rct in 1:19
    #    new_rct = phreeqc_CO2_alk_reaction_ratios[end] + 
    #        phreeqc_CO2_alk_reaction_ratios[2]
    #    append!(phreeqc_CO2_alk_reaction_ratios,new_rct)
    #end 
    #n_phreeqc_CO2_alk_reaction_ratios = size(phreeqc_CO2_alk_reaction_ratios)[1]

    #phreeqc_carbonate_reaction_extents = [0., 1.e-7]
    #rxn_mult = sqrt(sqrt(10)) #sqrt(sqrt(10))
    #for i_rct in 1:24
    #    new_rct = phreeqc_carbonate_reaction_extents[end] * rxn_mult
    #    append!(phreeqc_carbonate_reaction_extents,new_rct)
    #end 
    #n_phreeqc_carbonate_reaction_extents = size(phreeqc_carbonate_reaction_extents)[1]

    #phreeqc_CO2_reaction_extents = deepcopy(phreeqc_carbonate_reaction_extents)#[0.,1.e-7]
    #n_phreeqc_CO2_reaction_extents = size(phreeqc_CO2_reaction_extents)[1]

    
    n_phreeqc_reactants = n_reaction_components
    =#

    # Kitchen Sink Minerals 

        mineral_names = ["Anorthite", "Albite", "K-Feldspar", # feldspars
            "Annite", "Phlogopite", # components of biotite, Fe++/Mg 
            "Enstatite", "Ferrosilite", "Wollastonite", # components of pyroxenes/amphiboles, Mg/Fe++/Ca 
            "Forsterite", "Fayalite", # components of olivine, Mg/Fe++
            "Quartz", "Hematite", # Fe+++ 
            "Montmor-Ca", "Montmor-Mg", "Montmor-Na", "Montmor-K", # secondary clays, contain Mg
            "Beidellite-Ca", "Beidellite-Mg", "Beidellite-Na", "Beidellite-K", # no Mg 
            #"Saponite-Ca", "Saponite-Mg", "Saponite-Na", "Saponite-K", # no Mg 
            "Illite", "Muscovite", "Pyrophyllite", "Laumontite",
            "Kaolinite", "Gibbsite", "SiO2(am)", #"Celadonite",
            "Calcite", "Dolomite"] # , "Cristobalite(alpha)"] # , "SiO2(am)"]
        Anorthite_mineral, Albite_mineral, K_Feldspar_mineral,
        Annite_mineral, Phlogopite_mineral,
        Enstatite_mineral, Ferrosilite_mineral, Wollastonite_mineral,
        Forsterite_mineral, Fayalite_mineral,
        Quartz_mineral, Hematite_mineral,
        Montmor_Ca_mineral, Montmor_Mg_mineral, Montmor_Na_mineral, Montmor_K_mineral, 
        Beidellite_Ca_mineral, Beidellite_Mg_mineral, Beidellite_Na_mineral, Beidellite_K_mineral,
        #Saponite_Ca_mineral, Saponite_Mg_mineral, Saponite_Na_mineral, Saponite_K_mineral,
        Illite_mineral, Muscovite_mineral, Pyrophyllite_mineral, Laumontite_mineral, 
        Kaolinite_mineral, Gibbsite_mineral, SiO2_am_mineral, #Celadonite_mineral,
        Calcite_mineral, Dolomite_mineral = #, Cristobalite_mineral =
                1, 2, 3,
                4, 5, 
                6, 7, 8, 
                9, 10,  
                11, 12, 
                13, 14, 15, 16,
                17, 18, 19, 20,
                21, 22, 23, 24,
                25, 26, 27,
                28, 29 
        n_minerals = length(mineral_names)
        n_primary_minerals = Hematite_mineral
        phreeqc_dissolving_minerals = [
            Anorthite_mineral, Albite_mineral, K_Feldspar_mineral,
            Annite_mineral, Phlogopite_mineral,
            Enstatite_mineral, Ferrosilite_mineral, Wollastonite_mineral,
            Forsterite_mineral, Fayalite_mineral,
            Calcite_mineral, Dolomite_mineral]
        n_phreeqc_dissolving_minerals = size(phreeqc_dissolving_minerals)[1]
        primary_minerals = [
            Anorthite_mineral, Albite_mineral, K_Feldspar_mineral,
            Annite_mineral, Phlogopite_mineral,
            Enstatite_mineral, Ferrosilite_mineral, Wollastonite_mineral,
            Forsterite_mineral, Fayalite_mineral,
            Quartz_mineral, Hematite_mineral]
        authigenic_minerals = [
            #Hematite_mineral, 
            Montmor_Ca_mineral, Montmor_Mg_mineral, 
            Montmor_Na_mineral, Montmor_K_mineral,
            Beidellite_Ca_mineral, Beidellite_Mg_mineral,
            Beidellite_Na_mineral, Beidellite_K_mineral,
            Illite_mineral, Muscovite_mineral, Pyrophyllite_mineral, Laumontite_mineral, 
            Kaolinite_mineral, Gibbsite_mineral,
            Calcite_mineral, SiO2_am_mineral]
        secondary_minerals = [
            Montmor_Ca_mineral, Montmor_Mg_mineral, 
            Montmor_Na_mineral, Montmor_K_mineral, 
            Beidellite_Ca_mineral, Beidellite_Mg_mineral, 
            Beidellite_Na_mineral, Beidellite_K_mineral,
            #Saponite_Ca_mineral, Saponite_Mg_mineral, 
            #Saponite_Na_mineral, Saponite_K_mineral,
            Illite_mineral, Muscovite_mineral, Pyrophyllite_mineral, Laumontite_mineral, 
            Kaolinite_mineral, Gibbsite_mineral, SiO2_am_mineral, 
            Calcite_mineral, Dolomite_mineral ] #, Cristobalite_mineral]
        secondary_silicate_minerals = [
            Montmor_Ca_mineral, Montmor_Mg_mineral,
            Montmor_Na_mineral, Montmor_K_mineral,
            Beidellite_Ca_mineral, Beidellite_Mg_mineral,
            Beidellite_Na_mineral, Beidellite_K_mineral,
            #Saponite_Ca_mineral, Saponite_Mg_mineral,
            #Saponite_Na_mineral, Saponite_K_mineral,
            Illite_mineral, Muscovite_mineral, Pyrophyllite_mineral, Laumontite_mineral, 
            Kaolinite_mineral, Gibbsite_mineral,
            SiO2_am_mineral] # , Celadonite_mineral] #,
            #Cristobalite_mineral]
        carbonate_minerals = [Calcite_mineral, Dolomite_mineral]
        pure_phase_secondary_silicate_minerals = [
            Muscovite_mineral, Pyrophyllite_mineral, Laumontite_mineral,
            Kaolinite_mineral, Gibbsite_mineral,
            SiO2_am_mineral]
        secondary_solid_solution_minerals = 
            [Beidellite_Ca_mineral, Beidellite_Mg_mineral,
            Beidellite_Na_mineral, Beidellite_K_mineral, 
            Montmor_Ca_mineral, Montmor_Mg_mineral,
            Montmor_Na_mineral, Montmor_K_mineral,
            Illite_mineral]
        cation_rich_clays = [ Montmor_Ca_mineral, Montmor_Mg_mineral, 
            Montmor_Na_mineral, Montmor_K_mineral, 
            Beidellite_Ca_mineral, Beidellite_Mg_mineral, 
            Beidellite_Na_mineral, Beidellite_K_mineral,
            #Saponite_Ca_mineral, Saponite_Mg_mineral, 
            #Saponite_Na_mineral, Saponite_K_mineral,
            Illite_mineral, Muscovite_mineral, Pyrophyllite_mineral, Laumontite_mineral, 
        ]
        Plagioclase_group, K_Feldspar_group, Mica_primary_group, 
            Pyroxene_group, Olivine_group, Refractory_group,
            Smectite_group, Mica_secondary_group, Laterite_group,
            Carbonate_group = 
            1,2,3,4,5,6,7,8,9,10
        mineral_groups = [ 
            Plagioclase_group, Plagioclase_group, 
            K_Feldspar_group,   
            Mica_primary_group, Mica_primary_group,
            Pyroxene_group, Pyroxene_group, Pyroxene_group, 
            Olivine_group, Olivine_group,
            Refractory_group, Refractory_group, 
            Smectite_group, Smectite_group, Smectite_group, Smectite_group,
            Smectite_group, Smectite_group, Smectite_group, Smectite_group,
            Smectite_group, Smectite_group, Laterite_group, Laterite_group,
            Laterite_group, Laterite_group, Laterite_group, 
            Carbonate_group, Carbonate_group ]
        mineral_group_names = [
            "Plagioclase",
            "K_Feldspar",
            "Mica_primary",
            "Pyroxene",
            "Olivine",
            "Refractory",
            "Smectites",
            "Mica_secondary",
            "Laterites",
            "Carbonates"]
        n_mineral_groups = 10
        #= solid_solution_names = ["Smectites"] 
        n_solid_solutions = 1
        Smectite_solid_solution = 1
        solid_solution_sets = fill([],n_solid_solutions)
        solid_solution_sets[Smectite_solid_solution] =
            [Beidellite_Ca_mineral, Beidellite_Mg_mineral,
                Beidellite_Na_mineral, Beidellite_K_mineral, 
                #Montmor_Ca_mineral, Montmor_Mg_mineral,
                #Montmor_Na_mineral, Montmor_K_mineral,
                Illite_mineral] =#

        # from wen 2022 
            
        mica_dissolution_solid_rate_constant = 2.56e-13
        pyroxene_dissolution_solid_rate_constant = 1.0e-10 # "Enstatite"
        olivine_dissolution_solid_rate_constant = 1.2e-12 # "Forsterite"
        n_phreeqc_dissolving_phases = 8
        area_to_volume_ratio = 1.e2
        lasaga_rates_pH5 = area_to_volume_ratio .* # mol/m2 s 
              [ 5.6e-9, # Anorthite
                5.6e-9, # Lasaga said 1.19e-11 but plag is solid sol w anorthite, # Albite
                1.67e-12, # K-Feldspar
                mica_dissolution_solid_rate_constant, # Annite biotite with phlogopite: K, Fe2+ 
                mica_dissolution_solid_rate_constant, # Phlogopite: K, Mg 
                pyroxene_dissolution_solid_rate_constant, # Enstatite orthopyroxene with ferrosilite: MgSiO3
                pyroxene_dissolution_solid_rate_constant, # Ferrosilite: FeSiO3 
                pyroxene_dissolution_solid_rate_constant, # Wollastonite olivine with forsterite and fayalite: CaSiO3
                olivine_dissolution_solid_rate_constant, # Forsterite: Mg2SiO4
                olivine_dissolution_solid_rate_constant, # Fayalite: Fe2SiO4 
                0.0, # Quartz
                0.0, # Hematite: Fe2O3 
                0.0, # Montmor_Ca_mineral, 
                0.0, # Montmor_Mg_mineral, 
                0.0, # Montmor_Na_mineral, 
                0.0, # Montmor_K_mineral,
                0.0, # Beidellite_Ca_mineral, 
                0.0, # Beidellite_Mg_mineral, 
                0.0, # Beidellite_Na_mineral, 
                0.0, # Beidellite_K_mineral,
                0.0, # Illite_mineral, 
                0.0, # Muscovite_mineral
                0.0, # Pyrophyllite_mineral
                0.0, # Laumontite_mineral, 
                0.0, # Kaolinite_mineral, 
                0.0, # Gibbsite_mineral, 
                0.0, # SiO2_am_mineral
                3.e-6, # Calcite_mineral, (not from Lasaga)
                3.e-7] # Dolomite_mineral], (not from Lasaga)
        #= mineral_dissolution_pH_rate_orders = # yr-1
            [   0.6, # Anorthite
                0.6, # Albite
                1.0, # K-Feldspar
                0.6, # Annite biotite with phlogopite: K, Fe2+ 
                0.6, # Phlogopite: K, Mg 
                0.25, # Enstatite orthopyroxene with ferrosilite: MgSiO3
                0.25, # Ferrosilite: FeSiO3 
                0.25, # Wollastonite olivine with forsterite and fayalite: CaSiO3
                0.34, # Forsterite: Mg2SiO4
                0.34, # Fayalite: Fe2SiO4 
                0.0, # Quartz
                0.0, # Hematite: Fe2O3 
                0.0, # Montmor_Ca_mineral, 
                0.0, # Montmor_Mg_mineral, 
                0.0, # Montmor_Na_mineral, 
                0.0, # Montmor_K_mineral,
                0.0, # Beidellite_Ca_mineral, 
                0.0, # Beidellite_Mg_mineral, 
                0.0, # Beidellite_Na_mineral, 
                0.0, # Beidellite_K_mineral,
                0.0, # Illite_mineral, 
                0.0, # Muscovite_mineral
                0.0, # Pyrophyllite_mineral
                0.0, # Laumontite_mineral, 
                0.0, # Kaolinite_mineral, 
                0.0, # Gibbsite_mineral, 
                0.0, # SiO2_am_mineral
                0.6, # Calcite_mineral, 
                0.6] # Dolomite_mineral]
     
        mineral_neutral_dissolution_rate_constants = fill(0.,n_minerals)

        
        for i_mineral in 1:n_minerals # calculate mineral_neutral_dissolution_rate_constants
            if lasaga_rates_pH5[i_mineral] > 0.
                mineral_neutral_dissolution_rate_constants[i_mineral] = 
                    lasaga_rates_pH5[i_mineral] * # mol/m2 s
                    area_to_volume_ratio * # mol/m3 pw s
                    1. / 1.e2^mineral_dissolution_pH_rate_orders[i_mineral] # to pH 7
            end 
        end =#
            #=[   3.0e-5, # Anorthite
                6.0e-8, # Albite
                1.4e-9, # K-Feldspar
                mica_dissolution_solid_rate_constant, # Annite biotite with phlogopite: K, Fe2+ 
                mica_dissolution_solid_rate_constant, # Phlogopite: K, Mg 
                pyroxene_dissolution_solid_rate_constant, # Enstatite orthopyroxene with ferrosilite and Wollastonite: MgSiO3
                pyroxene_dissolution_solid_rate_constant, # Ferrosilite: FeSiO3 
                pyroxene_dissolution_solid_rate_constant, # Wollastonite: CaSiO3
                olivine_dissolution_solid_rate_constant, # Forsterite olivine: Mg2SiO4
                olivine_dissolution_solid_rate_constant, # Fayalite: Fe2SiO4 
                3.0e-9, # Quartz
                0.0e-6, # Hematite: Fe2O3 
                0.0e-6, # Montmor_Ca_mineral, 
                0.0e-6, # Montmor_Mg_mineral, 
                0.0e-6, # Montmor_Na_mineral, 
                0.0e-6, # Montmor_K_mineral,
                0.0e-6, # Beidellite_Ca_mineral, 
                0.0e-6, # Beidellite_Mg_mineral, 
                0.0e-6, # Beidellite_Na_mineral, 
                0.0e-6, # Beidellite_K_mineral,
                0.0e-6, # Illite_mineral, 
                0.0e-6, # Muscovite_mineral
                0.0e-6, # Pyrophyllite_mineral
                0.0e-6, # Laumontite_mineral, 
                0.0e+0, # Kaolinite_mineral, 
                0.0e+0, # Gibbsite_mineral, 
                0.0e+0, # SiO2_am_mineral
                1.0e-4, # Calcite_mineral, 
                1.0e-4] # Dolomite_mineral]=#

        #= mineral_dissolution_pH_start_points = # yr-1
            [   7., # Anorthite
                7., # Albite
                7., # K-Feldspar
                7., # Annite biotite with phlogopite: K, Fe2+ 
                7., # Phlogopite: K, Mg 
                7., # Enstatite orthopyroxene with ferrosilite: MgSiO3
                7., # Ferrosilite: FeSiO3 
                7., # Wollastonite olivine with forsterite and fayalite: CaSiO3
                7., # Forsterite: Mg2SiO4
                7., # Fayalite: Fe2SiO4 
                0.0, # Quartz
                0.0, # Hematite: Fe2O3 
                0.0, # Montmor_Ca_mineral, 
                0.0, # Montmor_Mg_mineral, 
                0.0, # Montmor_Na_mineral, 
                0.0, # Montmor_K_mineral,
                0.0, # Beidellite_Ca_mineral, 
                0.0, # Beidellite_Mg_mineral, 
                0.0, # Beidellite_Na_mineral, 
                0.0, # Beidellite_K_mineral,
                0.0, # Illite_mineral, 
                0.0, # Muscovite_mineral
                0.0, # Pyrophyllite_mineral
                0.0, # Laumontite_mineral, 
                0.0, # Kaolinite_mineral, 
                0.0, # Gibbsite_mineral, 
                0.0, # SiO2_am_mineral
                5.5, # Calcite_mineral, 
                5.5] # Dolomite_mineral]
        =#
        #mineral_dissolution_porewater_resat_constants = fill(0.,n_minerals)
        #mineral_dissolution_porewater_resat_constants[Calcite_mineral:Dolomite_mineral] .=
        #    1.e-1 # yr-1 

        mineral_relative_Li_affinity = fill(0.0, n_minerals)
        mineral_relative_Li_affinity[Beidellite_Mg_mineral] = 1.0
        mineral_relative_Li_affinity[Kaolinite_mineral] = 0.2
        mineral_Li_fractionation = fill(0.0, n_minerals)
        mineral_Li_fractionation[Beidellite_Mg_mineral] = 30.0 # solid - dissolved 
        mineral_Li_fractionation[Kaolinite_mineral] = 0.0

        rho_continent_crust = 2.7e6 # g / m3 

        mineral_component_stoic = fill(0.0, n_minerals, n_chemical_components)
        mineral_component_stoic[Anorthite_mineral, :] =
            [1.0, 0.0, 0.0, 0.0, 2.0, 1.0, 0.0, 0.0, 0.0, 4.0]
        #    CaO,  MgO, Na2O, K2O,SiO2,Al2O3,CO2,FeO,Fe2O3,H2O
        mineral_component_stoic[Albite_mineral, :] =
            [0.0, 0.0, 0.5, 0.0, 3.0, 0.5, 0.0, 0.0, 0.0, 2.0]
        mineral_component_stoic[K_Feldspar_mineral, :] =
            [0.0, 0.0, 0.0, 0.5, 3.0, 0.5, 0.0, 0.0, 0.0, 2.0]
        mineral_component_stoic[Annite_mineral, :] =
            [0.0, 0.0, 0.0, 0.5, 3.0, 0.5, 0.0, 3.0, 0.0, 6.0]
        mineral_component_stoic[Phlogopite_mineral, :] =
            [0.0, 3.0, 0.0, 0.5, 3.0, 0.5, 0.0, 0.0, 0.0, 6.0]
        mineral_component_stoic[Enstatite_mineral, :] =
            [0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0]
        mineral_component_stoic[Ferrosilite_mineral, :] =
            [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0]
        #    CaO,  MgO, Na2O, K2O,SiO2,Al2O3,CO2,FeO,Fe2O3,H2O
        mineral_component_stoic[Wollastonite_mineral, :] =
            [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0]
        mineral_component_stoic[Forsterite_mineral, :] =
            [0.0, 2.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2.0]
        mineral_component_stoic[Fayalite_mineral, :] =
            [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 2.0]
        mineral_component_stoic[Quartz_mineral, :] =
            [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        mineral_component_stoic[Hematite_mineral, :] =
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
        mineral_component_stoic[Montmor_Ca_mineral, :] =
            [0.165, 0.33, 0.0, 0.0, 4.0, 1.67 / 2, 0.0, 0.0, 0.0, 4.0]
        mineral_component_stoic[Montmor_Mg_mineral, :] =
            [0.0, 0.495, 0.0, 0.0, 4.0, 1.67 / 2.0, 0.0, 0.0, 0.0, 4.0]
        mineral_component_stoic[Montmor_Na_mineral, :] =
            [0.0, 0.33, 0.33 / 2.0, 0.0, 4.0, 1.67 / 2.0, 0.0, 0.0, 0.0, 4.0]
        mineral_component_stoic[Montmor_K_mineral, :] =
            [0.0, 0.33, 0.0, 0.33 / 2, 4.0, 1.67 / 2.0, 0.0, 0.0, 0.0, 4.0]
        mineral_component_stoic[Beidellite_Ca_mineral, :] =
            [0.165, 0.0, 0.0, 0.0, 3.67, 2.33 / 2, 0.0, 0.0, 0.0, 4.66]
        mineral_component_stoic[Beidellite_Mg_mineral, :] =
            [0.0, 0.165, 0.0, 0.0, 3.67, 2.33 / 2, 0.0, 0.0, 0.0, 4.66]
        mineral_component_stoic[Beidellite_Na_mineral, :] =
            [0.0, 0.0, 0.33 / 2, 0.0, 3.67, 2.33 / 2, 0.0, 0.0, 0.0, 4.66]
        mineral_component_stoic[Beidellite_K_mineral, :] =
            [0.0, 0.0, 0.0, 0.165, 3.67, 2.33 / 2, 0.0, 0.0, 0.0, 4.66]
        mineral_component_stoic[Illite_mineral, :] =
            [0.0, 0.25, 0.0, 0.6 / 2, 3.5, 2.3 / 2, 0.0, 0.0, 0.0, 5.0]
        mineral_component_stoic[Muscovite_mineral, :] =
            [0.0, 0.0, 0.0, 0.5, 3.0, 1.5, 0.0, 0.0, 0.0, 6.0]
        mineral_component_stoic[Pyrophyllite_mineral, :] =
            [0.0, 0.0, 0.0, 0.0, 4.0, 1.0, 0.0, 0.0, 0.0, 4.0]
        mineral_component_stoic[Laumontite_mineral, :] =
            [1.0, 0.0, 0.0, 0.0, 4.0, 1.0, 0.0, 0.0, 0.0, 8.0]
        mineral_component_stoic[Kaolinite_mineral, :] =
            [0.0, 0.0, 0.0, 0.0, 2.0, 1.0, 0.0, 0.0, 0.0, 5.0]
        mineral_component_stoic[Gibbsite_mineral, :] =
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 3.0]
        mineral_component_stoic[SiO2_am_mineral, :] =
            [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        #mineral_component_stoic[Celadonite_mineral, :] =
        #    [0.0, 1.0, 0.0, 0.5, 4.0, 0.5, 0.0, 0.0, 0.0, 4.0]
        mineral_component_stoic[Calcite_mineral, :] =
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
        mineral_component_stoic[Dolomite_mineral, :] =
            [1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0]
        #mineral_component_stoic[Cristobalite_mineral, :] =
        #    [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]

        mineral_molwts = fill(0.0, n_minerals)
        mineral_mol_per_m3 = fill(0.0, n_minerals)
        for i_mineral in 1:n_minerals
            for i_component in 1:n_chemical_components
                mineral_molwts[i_mineral] +=
                    mineral_component_stoic[i_mineral, i_component] *
                    chemical_component_molwts[i_component]
                mineral_mol_per_m3[i_mineral] =
                    rho_continent_crust * # g / m3 
                    1.0 / mineral_molwts[i_mineral] # mol / g -> mol / m3 
            end
        end






    #mineral_dissolution_solid_rate_constants .= 0.0
    #mineral_dissolution_solid_rate_constants[Phlogopite_mineral] = 
    #    1.e-4       

    Li_mol_per_m3 = rho_continent_crust / 6.
    #=
    mineral_component_volume_stoic = fill(0.,n_minerals,n_chemical_components)
    # min vol * stoic = component vol ; min vol = min mol * molwt 
    for i_mineral in 1:n_minerals
        for i_component in 1:n_chemical_components
            mineral_component_volume_stoic[i_mineral,i_component] = 
                mineral_component_stoic[i_mineral, i_component] * # mol component / mol mineral
                chemical_component_molwts[i_component] * # g component / mol mineral 
                1. / mineral_molwts[i_mineral] # g component / g mineral proportional to volume 
        end 
    end 
    =#

    #"Anorthite", "Albite", "K-Feldspar","Enstatite", "Forsterite", "Quartz",
    basement_mineral_fraction_endmembers = fill(0.0, 3, n_minerals)
    granite_endmember, basalt_endmember, ultramafic_endmember = 1,2,3 
    mica_mineral_fractions = [0.7,0.3]
    pyroxene_mineral_fractions = [0.4,0.4,0.2]
    olivine_mineral_fractions = [0.8,0.2]
    basement_mineral_fraction_endmembers[granite_endmember, 1:n_primary_minerals] =
        [0.18, # Ca-spar
        0.18,  # Na spar
        0.26,   # K spar 
        0.0239 * mica_mineral_fractions[1],  # Annite carries FeO 
        0.0239 * mica_mineral_fractions[2], # Phlogopite has MgO 
        0., # Enstatite 
        0., # Ferrosilite 
        0., # Wollastonite 
        0., # Forsterite
        0., # Fayalite
        0., # Quartz
        0.0122 ] # Hematite carries Fe2O3 
    basement_mineral_fraction_endmembers[basalt_endmember, 1:n_primary_minerals] =
        [0.2, # Ca-spar
        0.1,  # Na spar 
        0.0,   # K spar
        0.0,  # Annite # biotite 
        0.0, # Phlogopite 
        0.5 * pyroxene_mineral_fractions[1], # Enstatite # pyroxenes
        0.5 * pyroxene_mineral_fractions[2], # Ferrosilite 
        0.5 * pyroxene_mineral_fractions[3], # Wollastonite 
        0.1 * olivine_mineral_fractions[1], # Forsterite # olivine 
        0.1 * olivine_mineral_fractions[2], # Fayalite
        0.05, # Quartz
        0.0] # Hematite 
    basement_mineral_fraction_endmembers[ultramafic_endmember, 1:n_primary_minerals] =
        [0.0, # Ca-spar
        0.0,  # Na spar 
        0.0,   # K spar
        0.0,  # Annite # biotite
        0.0, # Phlogopite 
        0.1 * pyroxene_mineral_fractions[1], # Enstatite # pyroxenes
        0.1 * pyroxene_mineral_fractions[2], # Ferrosilite 
        0.1 * pyroxene_mineral_fractions[3], # Wollastonite 
        0.9 * olivine_mineral_fractions[1], # Forsterite # olivine 90% 
        0.9 * olivine_mineral_fractions[2], # Fayalite
        0.0, # Quartz
        0.0] # Hematite 
    for i_endmember in 1:3 
        fraction_tot = 0.
        basement_mineral_fraction_endmembers[i_endmember,Quartz_mineral] = 0.
        for i_mineral in 1:n_primary_minerals 
            fraction_tot += basement_mineral_fraction_endmembers[i_endmember,i_mineral]
        end 
        basement_mineral_fraction_endmembers[i_endmember,Quartz_mineral] = 1. - fraction_tot 
    end 
    depleted_mineral_fractions = fill(0.,n_minerals)
    depleted_mineral_fractions[Quartz_mineral] = 1.
    depleted_carbonated_mineral_fractions = fill(0.,n_minerals)
    depleted_carbonated_mineral_fractions[Quartz_mineral] = 0.5
    depleted_carbonated_mineral_fractions[Dolomite_mineral] = 1.e-5
    depleted_carbonated_mineral_fractions[Calcite_mineral] = 1. - 
        depleted_carbonated_mineral_fractions[Quartz_mineral] - 
        depleted_carbonated_mineral_fractions[Dolomite_mineral]
    
    #basement_mineral_fractions = depleted_basement_mineral_fractions

    # Tracers 
    n_tracers = 2; Li_tracer, Li7_tracer = 1, 2
    tracer_names = ["Li","Li7"]
    lithium_bearing_minerals = [Phlogopite_mineral,Beidellite_Mg_mineral,Kaolinite_mineral]
    lithium_bearing_primary_minerals = [Phlogopite_mineral]
    lithium_bearing_secondary_minerals = [Beidellite_Mg_mineral,Kaolinite_mineral]
    dissolved_Li_2_Mg_ratio = 1.2E-3 # 100 uM -> 120 nM 
    basement_mineral_tracer_volume_fractions = fill(0.,n_minerals,n_tracers)
    basement_mineral_tracer_volume_fractions[Phlogopite_mineral,Li_tracer:Li7_tracer] .= 
        1000. * # mg / Kg, from Zhang 2021  
        1.e-6 #* # g Li / g mineral  

    #Be_10_influx, Be_10_inventory, Be_10_concentration = 1, 2, 3
    # E12 atoms / m2 yr, E12 atoms,E12 atoms / m3 

    n_layers = 2
    regolith_layer, saprolite_layer, bedrock_layer = 1, 2, 3
    layer_names = ["Regolith","Saprolite","Bedrock"]

    box_width = 4.e5

    phreeqc_input_file_name = "/Users/archer/Synched/papers/gridplates/mineral_solubilities/auto.in"
    phreeqc_output_file_name = "/Users/archer/Synched/papers/gridplates/mineral_solubilities/auto.out"
    exec_command = "/Users/archer/Synched/papers/gridplates/mineral_solubilities/run_auto_phreeqc.x"
    #phreeqc_input_file_name = "/Volumes/RAMDisk/auto.in"
    #phreeqc_output_file_name = "/Volumes/RAMDisk/auto.out"
    #exec_command = "/Volumes/RAMDisk/run_auto_phreeqc.x"
  
    minimum_phreeqc_reaction_extent = 1.e-7


    #function save_results_block(file_name,results_block)
    #    BSON.@save 

    temperature_driver_offset = 0.0
    atm_CO2 = 300.0
    soil_CO2_factor = 1.25 
    
    basement_mineral_basalt_fraction = 0.01
    basement_calcite_fraction = 0.25
    basement_dolomite_fraction = 0.05
    precip_minus_evap = 1.5 # m / yr 
    groundwater_fraction = 1.0
    #erosion_rate_parameter = 1. 
    saprolite_thickness = 30.0
    aquifer_recharge_fraction = 0.25
    saprolite_erosion_constant = 3.e-3 # m / yr  
    saprolite_erosion_efold_depth = 3.0 # meter 
    land_regolith_flow_parameter = 3.e2
    mineral_reaction_rate_parameter = 1.0
    #secondary_mineral_redissolution_rate_parameter = 1.0
    n_floodplain_boxes = 4

    #= 
    memories_of_phreeqc = Array{phreeqc_output_struct,1}(undef,0)
    phreeqc_memory_pointers = Dict()

    BSON.@save "phreeqc_memory_pointers_6pd_2.bson" phreeqc_memory_pointers
    n_phreeqc_memory_files = Int(floor(size(memories_of_phreeqc)[1] / 1e6 )) +1
    for i_file in 1:n_phreeqc_memory_files
        filename = "memories_6pd_2." * string(i_file) * ".bson"
        begin_pointer = Int((i_file - 1) * 1e6 + 1)
        end_pointer = begin_pointer + 1000000 - 1
        if end_pointer > size(memories_of_phreeqc)[1]
            end_pointer = size(memories_of_phreeqc)[1]
        end 
        memory_block = memories_of_phreeqc[begin_pointer:end_pointer]
        BSON.@save filename memory_block
        println("wrote ", filename)
    end 
    
    memories_of_phreeqc = Array{phreeqc_output_struct,1}(undef,0)
    phreeqc_memory_pointers = Dict()

    BSON.@load "phreeqc_memory_pointers_6pd_2.bson" phreeqc_memory_pointers
    n_phreeqc_memory_files = length(phreeqc_memory_pointers)
    n_phreeqc_memory_files = Int(floor(n_phreeqc_memory_files / 1e6 )) + 1
    memories_of_phreeqc = Array{phreeqc_output_struct,1}(undef,0)
    for i_file in 1:n_phreeqc_memory_files
        filename = "memories_6pd_2." * string(i_file) * ".bson"
        println("reading ", filename)
        BSON.@load filename memory_block
        append!(memories_of_phreeqc, memory_block)
    end  
    =#

  
    if @isdefined(memories_of_phreeqc) == false
        fuckme = "hello"
        global memories_of_phreeqc = Array{phreeqc_output_struct,1}(undef,0)
        global phreeqc_memory_pointers = Dict()
        #cd("../mineral_solubilities/april_13")
        #phreeqc_memory_pointers, memories_of_phreeqc = 
        #    read_phreeqc_memories( "april_11_39" )
    end 


    n_memories_of_phreeqc_filled = size(memories_of_phreeqc)[1]
    println("loaded run bank ", n_memories_of_phreeqc_filled)
    #n_interpolated_phreeqc_runs = 0
    #n_total_phreeqc_runs = 0

model_duration = 1.e8
    float_year = 1.5e8
    spike_year = 1.5e8
    spike_duration = 2.e4 
    model_time_stages = [float_year, spike_year, spike_duration,model_duration]
    #regolith_transportation_depth = 2.0

    function mudscape(parameter_list, model_time_stages)
        results_block =
            mudscape(
                parameter_list[1], 
                parameter_list[2], parameter_list[3],
                parameter_list[4], parameter_list[5],
                parameter_list[6], parameter_list[7], parameter_list[8],
                parameter_list[9], parameter_list[10],
                parameter_list[11], parameter_list[12],
                parameter_list[13], parameter_list[14],
                model_time_stages)
        return results_block
end

enable_interpolated_phreeqc = true    
enable_check_interpolated_phreeqc = false 
check_interpolated_phreeqc_tolerance = 1.e-1
enable_multithread_phreeqc = true

#enable_default_solid_solutions = false
#enable_closed_CO2_system = true

# mudscape model 
function mudscape(
    temperature_driver_offset, 
    atm_CO2, soil_CO2_factor, 
    precip_minus_evap, 
    basement_mineral_basalt_fraction, 
    basement_calcite_fraction, basement_dolomite_fraction,
    saprolite_thickness, aquifer_recharge_fraction, 
    saprolite_erosion_constant, saprolite_erosion_efold_depth, 
    land_regolith_flow_parameter, mineral_reaction_rate_parameter, 
    n_floodplain_boxes,
    model_time_stages,
    restart_file_name = "", restart_time_point = 0,
    acid_added = 0., acidified_box_list = [])

# Setup 
    #global enable_interpolated_phreeqc = false

    enable_bedrock_erosion = true
    enable_saprolite_erosion = true
    enable_left_coastal_runoff = true            
    enable_right_coastal_runoff = true
    enable_interbox_transport = true     
    enable_primary_mineral_dissolution = true    
    enable_secondary_mineral_redissolution = false  
    enable_restart_from_file = false   
    #enable_saprolite_to_regolith_transport = true 
    #enable_growing_saprolite = true 
    #enable_lithium_dissolution = false  
    enable_lithium_reprecipitation = true 
    enable_depleted_mineral_initial_contition = false 
    enable_depleted_silicate_mineral_initial_contition = true 

    enable_saprolite_mineral_diagnostics = false  
    enable_mineral_diagnostics = false     
    enable_tracer_diagnostics = false 
    enable_component_diagnostics = false 
    enable_river_chemistry_diagnostics = false 
    enable_bulk_regolith_transport_diagnostics = false  

    float_year = model_time_stages[1]
    spike_year = model_time_stages[2]
    spike_duration = model_time_stages[3]
    model_duration = model_time_stages[4]

    main_time_step = 5.e3
    spike_time_step = 100.
    initial_regolith_thickness = 8.0
    coastal_abrasion_factor = 3.e1

    #baby_step_era = true 
    #reset_domain_now = false 
    #time_step_growth_rate = 1.
    #time_step_maximum = 1.e4
    time_steps = [main_time_step]
    time_points = [0.]
    CO2_eras = ["spinup"]
    spiked_already = false 
    while time_points[end] < model_duration
        time_step = time_steps[end]
        time_now = time_points[end]
        previous_CO2_era = CO2_eras[end]
        CO2_era = previous_CO2_era
        if time_points[end] >= float_year && previous_CO2_era == "spinup"
            CO2_era = "float_init"
        end 
        if previous_CO2_era == "float_init"
            CO2_era = "float"
        end 
        if time_points[end] >= spike_year &&
            spiked_already == false 
            spiked_already = true 
            CO2_era = "spike" 
        end 
        if previous_CO2_era == "spike"
            CO2_era = "float"
        end 
        if time_points[end] >= spike_year + spike_duration
            time_step = main_time_step 
        end 
        #if time_points[end] >= baby_step_duration
        #    time_step *= time_step_growth_rate
        #    time_step = min(time_step, time_step_maximum)
        #end
        #if time_points[end] >= spike_year + spike_duration
        #    time_step = baby_time_step
        #end 
        time_now += time_step
        append!(time_steps, time_step)
        append!(time_points, time_now)
        push!(CO2_eras, CO2_era)
    end 
    n_steps = size(time_steps)[1]
    #n_baby_steps = n_steps 

    box_surface_elevations = [0.,5000.,10000.,5000.] # ,1000.,900.,800.,700.,600.,500.,400.,300.,200.,100.,0.]
    elevation_increment = 1000. / n_floodplain_boxes 
    if n_floodplain_boxes > 0 
        elevation = 1000. 
        for i_floodplain_box in 1:Int(n_floodplain_boxes)
            append!(box_surface_elevations,elevation)
            elevation -= elevation_increment
        end 
    end 
    append!(box_surface_elevations, 0)
    n_boxes = size(box_surface_elevations)[1]
    box_temperatures = []
    for i_box in 1:n_boxes 
        append!(box_temperatures, 
            15. + # 20. - ( box_surface_elevations[i_box] / 10000. ) * 10. +
            temperature_driver_offset )
    end 
    box_base_temperatures = deepcopy(box_temperatures)

    river_watershed_fractions = fill(0.0, 2, n_boxes)
    river_watershed_fractions[1, 1:2] = [1.0, 1.0]
    river_watershed_fractions[2, :] = 1.0 .- river_watershed_fractions[1, :]
    boxes_in_play = []
    for i_box in 1:n_boxes 
        append!(boxes_in_play,i_box)
    end 
    n_boxes_in_play = n_boxes 

    box_coordinates = fill(0.,n_boxes)
  
    box_layer_thicknesses = fill(0.,n_layers,n_boxes) 
    box_layer_thicknesses[regolith_layer, :] .= initial_regolith_thickness
    box_layer_thicknesses[saprolite_layer, :] .= saprolite_thickness
    box_fluid_flushing_rates = fill(0.0, n_layers, n_boxes) # meters / year
    box_fluid_flushing_ages = fill(0.0, n_layers, n_boxes)
    box_pore_saturations = fill(0.0, n_layers, n_boxes)
    #aquifer_recharge_fraction = 0.25
    groundwater_fraction = 1.
    shallow_flow_time, aquifer_flow_time = 2. / 12. , 30. 

    box_erosion_rates = fill(0.,regolith_layer:bedrock_layer,n_boxes)
    regolith_transportation_depth = 2. # meters 

    box_porosities = fill(0.5,n_layers,n_boxes)
    box_solid_volumes = box_layer_thicknesses .* (1. .- box_porosities)
    box_mineral_fractions = fill(0.0, n_minerals, n_layers+1, n_boxes)
    box_mineral_volumes = fill(0.0, n_minerals, n_layers, n_boxes)
    basement_mineral_fractions = fill(0.0, n_minerals)
    basement_mineral_fractions[:] =
        basement_mineral_fraction_endmembers[granite_endmember, :] .* 
            (1.0 - basement_mineral_basalt_fraction - 
            basement_calcite_fraction - basement_dolomite_fraction) +
        basement_mineral_fraction_endmembers[basalt_endmember, :] .* 
            basement_mineral_basalt_fraction
    basement_mineral_fractions[Calcite_mineral] = 
        basement_calcite_fraction
    basement_mineral_fractions[Dolomite_mineral] =
        basement_dolomite_fraction
   
    for i_box in 1:n_boxes
        for i_layer in regolith_layer:saprolite_layer
            #if enable_depleted_mineral_initial_contition 
            #    box_mineral_fractions[:, i_layer, i_box] = 
            #        depleted_mineral_fractions[:]
            #elseif enable_depleted_silicate_mineral_initial_contition
                box_mineral_fractions[:, i_layer, i_box] = 
                    basement_mineral_fractions
                #    depleted_carbonated_mineral_fractions[:] .* 0.95 .+
                #basement_mineral_fraction_endmembers[granite_endmember, :] .* 0.04 .+
                #basement_mineral_fraction_endmembers[basalt_endmember, :] .* 0.01
            #else
            #    box_mineral_fractions[:, i_layer, i_box] = 
            #        basement_mineral_fractions[:]
            #end 
            box_mineral_volumes[:,i_layer, i_box] = update_mineral_volumes(
                box_mineral_fractions[:, i_layer, i_box],
                box_layer_thicknesses[i_layer, i_box],
                box_porosities[i_layer, i_box])
        end
        box_mineral_fractions[:, bedrock_layer, i_box] = 
            basement_mineral_fractions[:]
    end 

    soil_water_content = 0.3
    box_soil_CO2_values = []
    for i_box in 1:n_boxes
        log_soil_CO2 = log(atm_CO2 / 1.e6) / log(10.0) +
            exp(-3.0 * soil_water_content +
                -0.25 / soil_water_content) /
            (0.09 + exp(-0.34 * box_temperatures[i_box]))
        soil_CO2 = 10^log_soil_CO2 * 1.e6
        soil_CO2 *= soil_CO2_factor
        append!(box_soil_CO2_values,soil_CO2)
    end 
    box_precipitation_rates = fill(1.,n_boxes) 
    base_precipitation_rates = deepcopy(box_precipitation_rates)
    box_soil_water_contents = fill(0.3,n_boxes)
    base_soil_water_contents = deepcopy(box_soil_water_contents)
    box_solute_concentrations = fill(0.,n_solutes,n_layers,n_boxes)
    box_mineral_saturation_indices = fill(0.,n_minerals,n_layers,n_boxes)
    box_mineral_porewater_reaction_rates = fill(0.0, n_minerals, n_layers, n_boxes)
    box_mineral_reaction_extents_from_phreeqc = fill(0.0, n_minerals, n_layers, n_boxes)
    box_mineral_meters_for_residual_array = fill(0.0, n_minerals, n_layers, n_boxes)
    #box_mineral_dissolution_solid_rate_constants = fill(NaN, n_minerals, n_layers, n_boxes)
    box_CO2_uptake_bulk_rates = fill(0.0, n_layers, n_boxes)
    #box_alkalinity_production_bulk_rates = fill(0.0, n_layers, n_boxes)

    timepoints = fill(0.,n_steps+1)
    parameter_timeseries = fill(0.,14,n_steps+1)
    layer_thickness_timeseries = fill(0.,n_layers,n_boxes,n_steps+1)
    layer_thickness_timeseries[:,:,1] = box_layer_thicknesses 
    mineral_fraction_timeseries = fill(0.0, n_minerals, n_layers+1, n_boxes, n_steps + 1)
    for i_box in 1:n_boxes
        for i_mineral in 1:n_minerals
            mineral_fraction_timeseries[i_mineral,bedrock_layer,i_box,:] .= 
                basement_mineral_fractions[i_mineral]
        end 
    end
    mineral_reaction_rates_timeseries = fill(0.0, n_minerals, n_layers, n_boxes, n_steps + 1)
    mineral_saturation_state_timeseries = fill(0.0, n_minerals, n_layers, n_boxes, n_steps + 1)
    mineral_summary_timeseries = fill(0.,n_minerals,6,n_steps + 1) 
    component_summary_timeseries = fill(0.,n_chemical_components,6,n_steps + 1)
    # erosion, reaction, runoff, inventory, change inv, balance 
    tracer_summary_timeseries = fill(0.0, 11, n_steps + 1)
    # regolith C, sap C, reg alk, sap alk, tot C, Li src l/r, runoff l/r, % l/r  
    box_solute_timeseries = fill(0.0, n_solutes, n_layers, n_boxes, n_steps + 1)
    box_soil_CO2_timeseries = fill(0.0, n_boxes, n_steps + 1)
    box_temperature_timeseries = fill(0.0, n_boxes, n_steps + 1)
    box_precipitation_rate_timeseries = fill(0.0, n_boxes, n_steps + 1)
    river_solute_timeseries = fill(0.0, 2, n_solutes, n_steps + 1)
    regolith_age_timeseries = fill(0.0, n_boxes, n_steps + 1)
    erosion_source_timeseries = fill(0.0, regolith_layer:bedrock_layer, n_boxes, n_steps + 1)
    CO2_uptake_rate_timeseries = fill(0.0, n_layers, n_boxes, n_steps + 1)
    #alkalinity_production_bulk_rate_timeseries = fill(0.0, n_layers, n_boxes, n_steps + 1)
    #Be_10_timeseries = fill(0.0, 3, n_boxes, n_steps + 1)
    transport_matrix_all_boxes = fill(0.0, n_boxes, n_boxes)
    residual_matrix = fill(0.,n_boxes)

    solid_components_total = fill(0.,n_chemical_components) # m2 
    for i_box in 1:n_boxes
        for i_layer in 1:n_layers
            for i_mineral in 1:n_minerals
                solid_components_total += box_mineral_volumes[i_mineral,i_layer,i_box] .*
                    box_width .*
                    mineral_component_stoic[i_mineral,:]
            end 
            #=box_mineral_dissolution_solid_rate_constants[1:n_primary_minerals,i_layer,i_box] = 
                mineral_neutral_dissolution_rate_constants[1:n_primary_minerals] .*
                mineral_reaction_rate_parameter 
            for i_mineral in secondary_minerals
                box_mineral_dissolution_solid_rate_constants[i_mineral,i_layer,i_box] = 
                    mineral_neutral_dissolution_rate_constants[i_mineral] *
                    secondary_mineral_redissolution_rate_parameter
            end =#
        end 
    end 
    #=if enable_primary_mineral_dissolution == false
        for i_mineral in 2:n_minerals
            box_mineral_dissolution_solid_rate_constants[i_mineral] = 0.
        end 
    end 
    if enable_secondary_mineral_redissolution == false
        for i_box in 1:n_boxes
            for i_layer in 1:n_layers 
                for i_mineral in cation_rich_clays
                    box_mineral_dissolution_solid_rate_constants[i_mineral,i_layer,i_box] = 0.
                end 
            end 
        end
    end =#

    # Lithium setup 
    box_dissolved_del_Li7 = fill(0.,n_layers,n_boxes)   
    box_mineral_tracer_volume_fractions = fill(0.,n_minerals,n_tracers,n_layers,n_boxes)
    box_mineral_tracer_mole_fractions = fill(0.,n_minerals,n_tracers,n_layers,n_boxes)
    box_mineral_tracer_volumes = fill(0.,n_minerals,n_tracers,n_layers,n_boxes)
    box_mineral_tracer_volume_fluxes = fill(0.,n_minerals,n_tracers,n_layers,n_boxes)
    for i_box in 1:n_boxes
        for i_layer in 1:n_layers
            for i_mineral in lithium_bearing_primary_minerals
                for i_tracer in [Li_tracer,Li7_tracer]
                    box_mineral_tracer_volume_fractions[i_mineral,i_tracer,i_layer,i_box] = 
                        basement_mineral_tracer_volume_fractions[i_mineral,i_tracer]
                    box_mineral_tracer_mole_fractions[i_mineral,i_tracer,i_layer,i_box] = 
                        box_mineral_tracer_volume_fractions[i_mineral,i_tracer,i_layer,i_box] * # g Li / g min 
                        1. / 6. * # mol Li / g min 
                        mineral_molwts[i_mineral] # mol Li / mol min 
                    box_mineral_tracer_volumes[i_mineral,i_tracer,i_layer,i_box] = 
                        box_mineral_volumes[i_mineral,i_layer,i_box] * # m3 mineral / m2 
                        box_mineral_tracer_volume_fractions[i_mineral,i_tracer,i_layer,i_box] # m3 tracer / m3 mineral 
                        
                end 
            end 
        end 
    end 

    # erosional age
    box_erosional_age = fill(0.,n_layers,n_boxes)

    if restart_file_name !== "" 
        box_layer_thicknesses, box_mineral_fractions, box_mineral_volumes = 
            initialize_mudscape_from_results_file(file_name,time_point)
    end 
    
    
     
# time step 
    base_CO2_uptake_rate_recorded = false
    spike_C_ppm = 1250.
    base_CO2_uptake = NaN 
    base_atm_CO2 = atm_CO2 
    carbon_cycle_equivalent_buffer = 0.
    #base_T = surface_temperature
    base_LRFP = land_regolith_flow_parameter
    base_precip_minus_evap = precip_minus_evap
    base_SEC = saprolite_erosion_constant
    base_MRRP = mineral_reaction_rate_parameter
    time_of_spike = NaN
    #new_output_blocks = Array{phreeqc_box_struct,1}(undef,n_boxes)
    i_step = 1

    for i_step in 1:n_steps  # n_steps 
    # setup
        #n_total_phreeqc_runs = 0
        #n_interpolated_phreeqc_runs = 0
        time_step = time_steps[i_step]
        time_now = time_points[i_step]
        #if time_now > 1.e7
        #    global enable_interpolated_phreeqc = true
        #end 
        CO2_era = CO2_eras[i_step]
        #time_in_spike = NaN
        #if time_now >= baby_step_duration
        #    baby_step_era = false
        #    reset_domain_now = true 
        #end 
        if CO2_era == "float_init" # drawdown in eq with buffered ocean 
            base_CO2_uptake = tracer_summary_timeseries[5, i_step] 
            carbon_cycle_equivalent_buffer = 0.
            CO2_era = "float"
            println("starting to float CO2")
        end 
        if CO2_era == "spike"
            carbon_cycle_equivalent_buffer += spike_C_ppm
            time_of_spike = time_now 
            #CO2_era = "float"
            println("spiking CO2")
            #atm_CO2_to_liquid = spike_C_ppm * 0.25
            #atm_CO2_to_CaCO3 = spike_C_ppm * 0.15
        end 
        if CO2_era == "float" || CO2_era == "spike"
            airborne_fraction = 0.1 
            if time_of_spike == time_of_spike
                time_since_spike = time_now - time_of_spike
                airborne_fraction = 
                    0.25 * exp(-time_since_spike / 1000.0) +
                    0.15 * exp(-time_since_spike / 10000.0) +
                    0.1
            end 
            CO2_drawdown =
                (tracer_summary_timeseries[5, i_step] - base_CO2_uptake) *  
                0.1 / base_CO2_uptake * # Gton / yr
                time_step * # Gton C 
                0.5 # ppt from Gton 
            carbon_cycle_equivalent_buffer -= CO2_drawdown 
            atm_CO2 = base_atm_CO2 + 
                carbon_cycle_equivalent_buffer *
                airborne_fraction
            println("time, CO2, drawdown ", [time_now,atm_CO2, carbon_cycle_equivalent_buffer,
                tracer_summary_timeseries[5, i_step],
                tracer_summary_timeseries[5, i_step] - base_CO2_uptake])
            println("atm CO2 ", atm_CO2)
        end
        #if CO2_era == "spike" || CO2_era == "float"
        temperature_offset =
            4.0 * log(atm_CO2 / base_atm_CO2) / log(2.0)
        box_temperatures = box_base_temperatures .+ temperature_offset
        #surface_temperature = base_T + temperature_offset
        box_precipitation_rates = base_precipitation_rates .* 
            (1.0 + (temperature_offset * 0.02))

        #high_T_saturation = 0.6 + precip_minus_evap * 1.4
        for i_box in 1:n_boxes
            box_soil_water_contents[i_box] = base_soil_water_contents[i_box] * 
                box_precipitation_rates[i_box] / 
                base_precipitation_rates[i_box]
            log_soil_CO2 = log(atm_CO2 / 1.e6) / log(10.0) +
                exp(-3.0 * box_soil_water_contents[i_box] +
                    -0.25 / box_soil_water_contents[i_box]) /
                (0.09 + exp(-0.34 * box_temperatures[i_box]))
            box_soil_CO2_values[i_box] = 10^log_soil_CO2 * 1.e6 *
                soil_CO2_factor
        end 
         
        land_regolith_flow_parameter = base_LRFP *
            precip_minus_evap / base_precip_minus_evap
        saprolite_erosion_constant = base_SEC *
            precip_minus_evap / base_precip_minus_evap
        mineral_reaction_rate_parameter = base_MRRP *
            precip_minus_evap / base_precip_minus_evap
        #println("CO2 values ", [atm_CO2,soil_CO2])

        previous_mineral_volumes = deepcopy(box_mineral_volumes)
        previous_layer_thicknesses = deepcopy(box_layer_thicknesses)

        max_land_regolith_flow_parameter = land_regolith_flow_parameter * # meters / year, * H -> m2 / yr 
            regolith_transportation_depth * # dE
            1.0 / box_width * # m3 / m yr  
            time_step * # m3 / m  
            1.0 / box_width

        #if time_now >= spin_up_phase
        #    time_step = accelerated_time_step
        #end 
        #time_now += time_step 

        box_pore_volumes = box_porosities .* box_layer_thicknesses
        box_solid_volumes = (1.0 .- box_porosities) .* box_layer_thicknesses

        box_fluid_flushing_ages[regolith_layer,:] .= shallow_flow_time
        box_fluid_flushing_ages[saprolite_layer,:] .= aquifer_flow_time

        box_fluid_flushing_rates[regolith_layer,:] .= 
            precip_minus_evap * 
            groundwater_fraction * 
            (1. - aquifer_recharge_fraction) 
        box_fluid_flushing_rates[saprolite_layer,:] .= 
            precip_minus_evap * 
            groundwater_fraction * 
            aquifer_recharge_fraction
        
        #=
        for i_box in 1:n_boxes 
            regolith_relative_flow = box_layer_thicknesses[regolith_layer,i_box]
            saprolite_relative_flow = box_layer_thicknesses[saprolite_layer,i_box] *
                saprolite_relative_permeability
            box_fluid_flushing_rates[regolith_layer, i_box] = 
                precip_minus_evap * groundwater_fraction * 
                regolith_relative_flow / 
                ( regolith_relative_flow + saprolite_relative_flow )
            box_fluid_flushing_rates[saprolite_layer, i_box] = 
                precip_minus_evap * groundwater_fraction - 
                box_fluid_flushing_rates[regolith_layer, i_box]
            box_fluid_flushing_ages[:, i_box] = 
                box_pore_volumes[:, i_box] ./ 
                box_fluid_flushing_rates[:, i_box]
        end =#
        
            #box_fluid_flushing_rates[regolith_layer, i_box] = 
            #    saprolite_exchange_lifetime_factor * 
            #    box_pore_volumes[regolith_layer, i_box] *
            #    precip_minus_evap * groundwater_fraction *
            #    1. / ( box_pore_volumes[saprolite_layer, i_box] + 
            #        saprolite_exchange_lifetime_factor * 
            #        box_pore_volumes[regolith_layer, i_box] )
            #box_fluid_flushing_rates[saprolite_layer, i_box] = 
            #    precip_minus_evap * groundwater_fraction - 
            #    box_fluid_flushing_rates[regolith_layer, i_box]
            #box_fluid_flushing_ages[:, i_box] = 
            #    box_pore_volumes[:, i_box] ./ 
            #    box_fluid_flushing_rates[:, i_box]
            #=mean_flushing_age = ( box_pore_volumes[regolith_layer, i_box] +
                box_pore_volumes[saprolite_layer, i_box] ) / 
                (precip_minus_evap * groundwater_fraction) 
            saprolite_age_relative_to_regolith = 
                mean_flushing_age * 
                saprolite_exchange_lifetime_factor
            box_fluid_flushing_rates[saprolite_layer, i_box] = 
                box_pore_volumes[saprolite_layer, i_box] / 
                saprolite_age_relative_to_regolith
            box_fluid_flushing_rates[regolith_layer, i_box] = 
                precip_minus_evap * groundwater_fraction #- 
                #box_fluid_flushing_rates[saprolite_layer, i_box]
            box_fluid_flushing_ages[regolith_layer, i_box] = 
                box_pore_volumes[regolith_layer, i_box] / # m3 / m2 
                box_fluid_flushing_rates[regolith_layer, i_box] # m3 / m2 yr -> yr  
            box_fluid_flushing_ages[saprolite_layer, i_box] = 
                box_pore_volumes[saprolite_layer, i_box] / 
                box_fluid_flushing_rates[saprolite_layer, i_box] + 
                box_fluid_flushing_ages[regolith_layer, i_box]=#
        #end 

        #=box_fluid_flushing_ages[regolith_layer, :] =
            box_pore_volumes[regolith_layer, :] ./
            (precip_minus_evap * groundwater_fraction)

        box_fluid_flushing_timescales[saprolite_layer, :] .=     
            saprolite_exchange_lifetime_factor
        box_fluid_flushing_rates[saprolite_layer, :] .= 
            box_pore_volumes[saprolite_layer,:] ./
            box_fluid_flushing_timescales[saprolite_layer, :]
        box_fluid_flushing_rates[regolith_layer, :] .= 
            precip_minus_evap * groundwater_fraction 
        box_fluid_flushing_rates[regolith_layer, :] -= 
            box_fluid_flushing_rates[saprolite_layer, :]
        box_fluid_flushing_timescales[regolith_layer, :] = # tau = vol / flux 
            box_pore_volumes[regolith_layer, :] ./      
            box_fluid_flushing_rates[regolith_layer, :]=#

        for i_box in 1:n_boxes # update erosion fluxes 
            if enable_bedrock_erosion
                #if time_now <= baby_step_duration #|| 
                    #box_layer_thicknesses[regolith_layer, i_box] > 3.
                box_erosion_rates[saprolite_layer,i_box] = 
                    saprolite_erosion_constant *
                    #erosion_rate_parameter * 
                    exp( - box_layer_thicknesses[regolith_layer,i_box] / 
                        saprolite_erosion_efold_depth )
                    # m3 solid / m2 year 
                box_erosion_rates[bedrock_layer,i_box] = 
                    box_erosion_rates[saprolite_layer,i_box] +
                    ( saprolite_thickness - 
                    box_layer_thicknesses[saprolite_layer, i_box] ) /
                    time_step 
            end 
        end
        #println("box check ", box_layer_thicknesses[:,1])
        
        saved_lu_transport_all_boxes = 0 # place_holder
        saved_lu_transport_interior_boxes = 0
        #transport_matrix .= 0.
 
    # derivation
        # saprolite point-by-point 
            # saprolite erosion rate is specified, 
            # bedrock erosion resets thickness to the spec value,
            # but prior to any reaction.  slop is dumped into a 
            # mass-conserving small deviation of sap thickness. 
            # compute new box inventories without erosion of saprolite
            # then impose sap erosion, compute new thickness 
            # Inv' = Inv + E_basement F_basement dt - K Inv' dt
            # Inv' ( 1 + K dt ) = Inv + E_basement F_basement dt
            # Inv' = (Inv + E_basement F_basement dt) / (1 + K dt)
            # 
            # compute new fractions in the saprolite 
            # run phreeqc, get rates of secondary mineral formation / dissolution 
            # set inventories of secondary minerals, recompute fractions 
            # E_saprolite = new sap inv - old sap inv (const sap thickness)
        # regolith transportation matrix 
            # baby step period. fluxes proportional to thickness throughout the domain.  
            # diffuse each mineral independently. 
            # elevation gradient is based on initial thicknesses to avoid nonlinearity
            # M' = mineral_solid_volume end of step, T = box_layer_thickness, beginning of step 
            # D T1 dt/dx2 (E1+T1-E2-T2) = m mineral / m2 yr 
            # fraction m = M (vol min) / T / (1-por)
            # case 1, E1 > E2, flow from 1 -> 2 is regulated by T1, as evaluated at 1
            # M'1 = M1 + Es1 dt Fms - K dt M'1 - D (1/(1-por)) dt/dx2 (E1+T1-E2-T2) M'1  
            # M'2 = M2 + Es2 dt Fms - K dt M'2 + D (1/(1-por)) dt/dx2 (E1+T1-E2-T2) M'1 - D (1/(1-por)) dt/dx2 (E2+T2) M'2
            # grouping terms 
            # M'1 ( 1 + K dt + D (1/(1-por)) dt/dx2 (E1+T1-E2-T2) ) + M'2 (O.)                         = M1 + Es1 dt Fms 
            # M'1 (-D (1/(1-por)) dt/dx2 (E1+T1-E2-T2))            + M'2 (1 + K dt + D (1/(1-por)) dt/dx2 (E2+T2) = M2 + Es2 dt Fms 

            # case 2, if E1 < E2, flow from 2 -> 1 is regulated by T2, as evaluated at 1 
            # M'1 = M1 + Es1 dt Fms - K dt M'1 - D (1/(1-por)) dt / dx2 (E1+T1-E2-T2) M'2 
            # M'2 = M2 + Es2 dt Fms - K dt M'2 + D (1/(1-por)) dt / dx2 (E1+T1-E2-T2) M'2 - D (1/(1-por)) dt / dx2 (E2+T2) M'2
            # grouping terms
            # M'1 ( 1 + K dt ) + M'2 ( D (1/(1-por)) dt/dx2 (E1+T1-E2-T2) )                           = T1 + Es1 dt Fb 
            # M'1 (0.)         + M'2 (1 + K dt - D (1/(1-por)) dt/dx2 (E1+T1-E2-T2) + D  (1/(1-por)) dt/dx2 E2) = T2 + Es2 dt Fb 

            # interior transport period. fluxes independent of regolith thickness.   
            # composite calculation neglect change in volume, figure that the feedback to transport would
            # be minor (mostly driven by bedrock elevation gradients). T=total thickness, M=mineral thickness
            # assume boundary condition influx from point 0 
            # T'1 = T1 + Es1 dt - D Tmax dt/dx2 (E1+T'1-E2-T'2) + D T0 dt/dx2 (E0+T0-E1-T1)
            # T'2 = T2 + Es2 dt + D Tmax dt/dx2 (E1+T'1-E2-T'2) - D Tmax dt/dx2 (E2+T'2)

            # T'1 ( 1 + D Tmax dt/dx2 ) + T'2 ( - D Tmax dt/dx2 ) = 
            #                  T1 + Es1 (1-por) dt/dx2 - D Tmax dt/dx2 (E1-E2) + D T0 dt/dx2 (E0+T0-E1-T1)
            # T'1 ( - D Tmax dt/dx2 ) + T'2 ( 1 + D Tmax dt/dx2 + D Tmax dt/dx2 ) = 
            #                  T2 + Es2 dt + D Tmax dt/dx2 (E1-E2) - D Tmax dt/dx2 E2
        
            # tot_flux = D dt/dx (E1+T'1-E2-T'2)
            # tot_runoff = D dt/dx (E2+T'2) / 2 

            # Upstream differencing. Case 1 flow from 1 -> 2 
            # M'1 = M1 + Es1 dt Fbm - Km dt M'1 - tflux1 M'1/T'1 / dx + tflux0 M0/T0
            # M'2 = M2 + Es1 dt Fbm - Km dt M'2 + tflux1 M'1/T'1 / dx - trunoff * M'2 / T'2
            
            # M'1 ( 1 + Km dt + tflux1 / T'1 / dx ) + M'2 ( 0 ) = 
            #                  M1 + Es1 dt Fbm + tflux0 M0/T0 
            # M'1 ( - tflux1 / T'1 / dx ) + M'2 ( 1 + Km dt + trunoff / T'2 ) =
            #                  M2 + Es2 dt Fbm 

            # Case 2 flow from 2 -> 1, eval at 1, tflux1 < 0
            # M'1 = M1 + Es1 dt Fbm - Km dt M'1 - tflux1 M'2/T'2 / dx + tflux0 M0/T0
            # M'2 = M2 + Es1 dt Fbm - Km dt M'2 + tflux1 M'2/T'2 / dx - trunoff * M'2 / T'2
            
            # M'1 ( 1 + Km dt ) + M'2 ( tflux1 / T'1 / dx ) = 
            #                  M1 + Es1 dt Fbm + tflux0 M0/T0 
            # M'1 ( 0 ) + M'2 ( 1 + Km dt - tflux1 / T'1 / dx + trunoff / T'2 ) =
            #                  M2 + Es2 dt Fbm 

            # check to see if the sum of the new mineral thickness is close enough to T', no holes 

            # for all primaries, find the dissolved mass fraction -> kaolinite or whatever 
            # composite reaction volume change comp_dV = Sum primary(K dt T'1 Fdiss)
            # T'1 = T1 + Es1 dt - dt(Km M'1 Fdm + Kn N'1 Fdn) - D dt / dx2 (E1+T1-E2-T2)
            #                           oops where do these come from?
            # T'2 = T2 + Es2 dt - Sum primary(K dt T'2 Fdiss) + D dt / dx2 (E1+T1-E2-T2) - D dt / dx2 (E2+T2)

            # what if T'1 < Tmin and T'2 > Tmin 
            # T'1 = T1 + Es1 dt Fs - K dt T'1 - D dt / dx2 (E1+T1-E2-T2) T'1  
            # T'2 = T2 + Es2 dt Fs - K dt T'2 + D dt / dx2 (E1+T1-E2-T2) T'1 - D dt / dx2 Tmin (E2+T2) 
            # T'1 ( 1 + K dt + D dt/dx2(E1+T1-E2-T2) ) + T'2 (O.)                         = T1 + Es1 dt Fs 
            # T'1 (-D dt/dx2 (E1+T1-E2-T2))            + T'2 (1 + K dt + D dt/dx2 (E2+T2) = T2 + Es2 dt Fs 



            # simultaneous treatment of all minerals  
            # T'm1 = Tm1 + Es1 dt Fsm - K dt T'm1 - D Tmax dt / dx2 (E1+T'm1+T'n2-E2-T'm2-T'n2)
            # T'n1 = Tn1 + Es1 dt Fsn - K dt T'n1 - D Tmax dt / dx2 (E1+T'm1+T'n2-E2-T'm2-T'n2)
            # T'm2 = Tm2 + Es2 dt Fms - K dt T'm2 + D Tmax dt / dx2 (E1+T'm1+T'n2-E2-T'm2-T'n2)
            # T'n1 oops fails because you need the mineral fractions at the end of the step for transport, nonlinear 

            # case 3, both T's higher than the max. applied value, so no longer dependent on T
            # T'1 = T1 + Es1 dt Fs - K dt T'1 - D dt / dx2 (E1+T'1-E2-T'2) Tmax  
            # T'2 = T2 + Es2 dt Fs - K dt T'2 + D dt / dx2 (E1+T'1-E2-T'2) Tmax  
            # grouping 
            # T'1 ( 1 + K dt + D Tmax dt/dx2 ) + T'2 ( - D Tmax dt/dx2 ) = T1 + Es1 dt Fs - D Tmax dt/dx2 (E1-E2)
            # T'1 ( - D Tmax dt/dx2 ) + T'2 ( 1 + K dt + D Tmax dt/dx2 ) = T2 + Es2 dt Fs + D Tmax dt/dx2 (E1-E2)
            # oops fails because other minerals affect the elevation hence the transport 


        # regolith reacton solver after the fact 

    # entrain bedrock into saprolite 
        for i_box in 1:n_boxes  
            for i_mineral in 1:n_minerals 
                new_thickness =
                    ( box_mineral_volumes[i_mineral,saprolite_layer,i_box] + 
                    box_erosion_rates[bedrock_layer,i_box] * 
                        basement_mineral_fractions[i_mineral] * 
                        time_step ) #
                #= if i_mineral in primary_minerals && enable_primary_mineral_dissolution 
                    new_thickness /=
                        ( 1. + 
                        box_mineral_dissolution_solid_rate_constants[i_mineral,saprolite_layer,i_box] *
                            time_step )
                end =#
                box_mineral_volumes[i_mineral,saprolite_layer,i_box] = 
                    new_thickness
            end 
            box_mineral_fractions[:,saprolite_layer,i_box],
                box_solid_volumes[saprolite_layer,i_box],
                box_layer_thicknesses[saprolite_layer, i_box] =
                update_mineral_fractions(box_mineral_volumes[:,saprolite_layer,i_box],
                    box_porosities[saprolite_layer,i_box])
        end 
    # run phreeqc on saprolite   
        if enable_interpolated_phreeqc
            for i_box in 1:n_boxes
                #println("saprolite box = ", i_box)
                output_block = # , n_new_total_phreeqc_runs, n_new_interpolated_phreeqc_runs = 
                #output_block_interp, n_new_total_phreeqc_runs, n_new_interpolated_phreeqc_runs = 
                phreeqc_box( 
                    box_mineral_fractions[:, saprolite_layer, i_box], 
                    box_layer_thicknesses[saprolite_layer,i_box],
                    box_porosities[saprolite_layer,i_box], 
                    box_fluid_flushing_rates[saprolite_layer,i_box], 
                    box_fluid_flushing_ages[saprolite_layer,i_box], 
                    box_temperatures[i_box], box_soil_CO2_values[i_box],
                    time_step, 1 )
                #println("Alk ", output_block.solute_concentrations[Alk_solute])
                #n_total_phreeqc_runs += n_new_total_phreeqc_runs
                #n_interpolated_phreeqc_runs += n_new_interpolated_phreeqc_runs
                box_pore_saturations[saprolite_layer,i_box] = 
                    output_block.pore_saturation
                box_solute_concentrations[:,saprolite_layer,i_box] = 
                    output_block.solute_concentrations
                box_mineral_saturation_indices[:,saprolite_layer,i_box] = 
                    output_block.mineral_saturation_indices
                box_mineral_reaction_extents_from_phreeqc[:, saprolite_layer, i_box] = 
                    output_block.mineral_reaction_extents
                box_mineral_porewater_reaction_rates[:,saprolite_layer,i_box] = 
                    output_block.mineral_porewater_reaction_rates
                box_mineral_meters_for_residual_array[:,saprolite_layer,i_box] = 
                    output_block.mineral_meters_for_residual_array
                box_CO2_uptake_bulk_rates[saprolite_layer, i_box] = 
                    output_block.CO2_uptake_extent
                #box_alkalinity_production_bulk_rates[saprolite_layer, i_box] =
                #    output_block.summary[2]
            end
        else # direct, not interpolated 
            Threads.@threads for i_box in 1:n_boxes
                #println()
                #println("saprolite box = ", i_box)
                thread_number = Threads.threadid()
                #=box_solute_concentrations[:,saprolite_layer,i_box], 
                    box_mineral_saturation_indices[:, saprolite_layer, i_box],
                    box_mineral_reaction_extents_from_phreeqc[:, saprolite_layer, i_box],
                    box_mineral_porewater_reaction_rates[:,saprolite_layer,i_box],
                    box_mineral_meters_for_residual_array[:,saprolite_layer,i_box],
                    #box_mineral_dissolution_solid_rate_constants[:,saprolite_layer,i_box],
                    box_CO2_uptake_bulk_rates[saprolite_layer, i_box],
                    box_alkalinity_production_bulk_rates[saprolite_layer, i_box] =# # eq / year 
                output_block = # , n_new_total_phreeqc_runs, n_new_interpolated_phreeqc_runs = 
                phreeqc_box( 
                    box_mineral_fractions[:, saprolite_layer, i_box], 
                    box_layer_thicknesses[saprolite_layer,i_box],
                    box_porosities[saprolite_layer,i_box], 
                    box_fluid_flushing_rates[saprolite_layer,i_box], 
                    box_fluid_flushing_ages[saprolite_layer,i_box], 
                    box_temperatures[i_box], box_soil_CO2_values[i_box],
                    time_step,thread_number)
                #n_total_phreeqc_runs += n_new_total_phreeqc_runs
                #n_interpolated_phreeqc_runs += n_new_interpolated_phreeqc_runs
                #end    
                #for i_box in 1:n_boxes
                box_pore_saturations[saprolite_layer,i_box] = 
                    output_block.pore_saturation
                box_solute_concentrations[:,saprolite_layer,i_box] = 
                    output_block.solute_concentrations
                box_mineral_saturation_indices[:,saprolite_layer,i_box] = 
                    output_block.mineral_saturation_indices
                box_mineral_reaction_extents_from_phreeqc[:, saprolite_layer, i_box] = 
                    output_block.mineral_reaction_extents
                box_mineral_porewater_reaction_rates[:,saprolite_layer,i_box] = 
                    output_block.mineral_porewater_reaction_rates
                box_mineral_meters_for_residual_array[:,saprolite_layer,i_box] = 
                    output_block.mineral_meters_for_residual_array
                box_CO2_uptake_bulk_rates[saprolite_layer, i_box] = 
                    output_block.CO2_uptake_extent
                #box_alkalinity_production_bulk_rates[saprolite_layer, i_box] =
                #    output_block.summary[2]
                #if new_output_blocks[i_box].origin == "direct"
                #    match_set = matching_driving_parameter_set(
                #        new_output_blocks[i_box].driving_parameters)
                #    if (match_set == match_set ) == false
                #        #println("saving saprolite ", i_box, " ", n_memories_of_phreeqc_filled)
                #        n_memories_of_phreeqc_filled += 1
                #        memories_of_phreeqc[n_memories_of_phreeqc_filled] = 
                #            new_output_blocks[i_box]
                #    end 
                #end 
            end
            #n_total_phreeqc_runs += n_boxes
        end 
    # react saprolite minerals 
        for i_box in 1:n_boxes 
            for i_mineral in 1:n_minerals 
                box_mineral_volumes[i_mineral,saprolite_layer,i_box] += 
                    box_mineral_meters_for_residual_array[i_mineral, saprolite_layer, i_box]
            end 
            box_mineral_fractions[:,saprolite_layer,i_box],
            box_solid_volumes[saprolite_layer, i_box],
            box_layer_thicknesses[saprolite_layer, i_box] =
                update_mineral_fractions(box_mineral_volumes[:,saprolite_layer,i_box],
                    box_porosities[saprolite_layer,i_box])
        end 
    # saprolite to regolith transport
        if enable_saprolite_erosion 
            for i_box in 1:n_boxes
                #box_erosion_rates[regolith_layer,i_box] = 
                #    ( box_layer_thicknesses[saprolite_layer,i_box] - 
                #    previous_layer_thicknesses[saprolite_layer,i_box] ) * 
                #    ( 1. - box_porosities[saprolite_layer,i_box])  / 
                #    time_step
                for i_mineral in 1:n_minerals
                    box_mineral_volumes[i_mineral,saprolite_layer,i_box] -= 
                        box_erosion_rates[saprolite_layer, i_box] *
                        box_mineral_fractions[i_mineral,saprolite_layer,i_box] *
                        time_step 
                end 
                box_mineral_fractions[:,saprolite_layer,i_box],
                box_solid_volumes[saprolite_layer, i_box],
                box_layer_thicknesses[saprolite_layer, i_box] =
                    update_mineral_fractions(box_mineral_volumes[:,saprolite_layer,i_box],
                        box_porosities[saprolite_layer,i_box])
                #box_mineral_fractions[:, regolith_layer, i_box],
                #box_solid_volumes[regolith_layer, i_box],
                #box_layer_thicknesses[regolith_layer, i_box] =
                #    update_mineral_fractions(box_mineral_volumes[:, regolith_layer, i_box],
                #        box_porosities[regolith_layer, i_box])
            end 
        end 
    # saprolite solids diagnostics 
        if enable_saprolite_mineral_diagnostics
            println(); println("Saprolite Mineral Diagnostics")
            bedrock_erosive_mineral_sources = fill(0.,n_minerals) # m3 / yr 
            saprolite_erosive_mineral_sinks = fill(0.,n_minerals)
            suspended_runoff_mineral_sinks = fill(0.0, n_minerals)
            old_mineral_inventories = fill(0.0, n_minerals)
            mineral_inventories = fill(0.0, n_minerals)
            mineral_inventory_changes = fill(0.0, n_minerals)
            mineral_reaction_rates = fill(0.0, n_minerals)
            for i_box in 1:n_boxes 
                for i_mineral in 1:n_minerals
                    bedrock_erosive_mineral_sources[i_mineral] += 
                        box_erosion_rates[bedrock_layer, i_box] * #  m3 solid / m2 year
                        basement_mineral_fractions[i_mineral] * # m / yr 
                        box_width * box_width # m3 solid / yr 
                    saprolite_erosive_mineral_sinks[i_mineral] += 
                        box_erosion_rates[saprolite_layer, i_box] * #  m3 solid / m2 year
                        box_mineral_fractions[i_mineral,saprolite_layer,i_box] * # m / yr 
                        box_width * box_width # m3 solid / yr 
                    old_mineral_inventories[i_mineral] += 
                        previous_mineral_volumes[i_mineral,saprolite_layer,i_box] * # meters solid 
                        box_width * box_width # m3 solid 
                    mineral_inventories[i_mineral] += 
                        box_mineral_volumes[i_mineral,saprolite_layer,i_box] *
                        box_width * box_width
                    mineral_reaction_rates[i_mineral] += # m3 / yr 
                        box_mineral_meters_for_residual_array[i_mineral,saprolite_layer,i_box] *
                        box_width * box_width / time_step 
                end     
            end 
            mineral_inventory_changes = 
                ( mineral_inventories .- 
                old_mineral_inventories ) ./ time_step 
            for i_mineral in 1:n_minerals 
                balance = bedrock_erosive_mineral_sources[i_mineral] -
                    saprolite_erosive_mineral_sinks[i_mineral] +
                    mineral_reaction_rates[i_mineral] -
                    suspended_runoff_mineral_sinks[i_mineral] -
                    mineral_inventory_changes[i_mineral]
                println(mineral_names[i_mineral],
                    " BEr,SEr,Rx,Rp,Chg,Bal ", [
                    bedrock_erosive_mineral_sources[i_mineral],
                    saprolite_erosive_mineral_sinks[i_mineral],
                    mineral_reaction_rates[i_mineral],
                    suspended_runoff_mineral_sinks[i_mineral],
                    mineral_inventory_changes[i_mineral],
                    balance])
            end 
            dissolving_solids_component_fluxes = fill(0.0, n_chemical_components)
            dissolving_solids_solute_fluxes = fill(0.0, n_chemical_components)
            for i_component in 1:n_chemical_components
                for i_mineral in 1:n_minerals
                    dissolving_solids_component_fluxes[i_component] +=
                        mineral_reaction_rates[i_mineral] * # m3 / year
                        mineral_mol_per_m3[i_mineral] * # mol mineral / yr 
                        mineral_component_stoic[i_mineral, i_component] # mol comp / yr 
                    #=if dissolving_solids_component_fluxes[i_component] == 
                        dissolving_solids_component_fluxes[i_component]
                    else 
                        error(i_component," ", i_mineral )
                    end
                    println(mineral_names[i_mineral], " ", 
                        chemical_component_names[i_component]," ", 
                        dissolving_solids_component_fluxes[i_component])=#
                end 
            end 
            runoff_flow = 0. #precip_minus_evap * (1. - groundwater_fraction) *
                #box_width * box_width 
            for i_box in 1:n_boxes
                runoff_flow +=
                    box_fluid_flushing_rates[saprolite_layer, i_box] *
                    box_width * box_width
            end 
            dissolving_solids_runoff_concentrations = fill(0.,n_solutes)
            for i_comp in 1:n_solute_aligned_components
                i_solute = i_comp
                dissolving_solids_solute_fluxes[i_solute] = # mol / yr 
                    - dissolving_solids_component_fluxes[i_comp] *
                    solute_component_stoiciometry[i_solute]
                dissolving_solids_runoff_concentrations[i_solute] =
                    dissolving_solids_solute_fluxes[i_solute] * # mol / yr 
                    1 / runoff_flow # yr / m3 -> mol / m3 
            end
            dissolved_runoff_solute_fluxes = fill(0.0, n_solutes)
            for i_box in 1:n_boxes
                for i_solute in 1:n_solutes
                    dissolved_runoff_solute_fluxes[i_solute] +=
                        box_solute_concentrations[i_solute, saprolite_layer, i_box] * # mol / m3 pw
                        box_fluid_flushing_rates[saprolite_layer,i_box] * # m / yr -> mol / m2 yr 
                        box_width * box_width
                end
            end 
            if enable_saprolite_mineral_diagnostics
                println()
                println("Saprolite solutes")
                for i_solute in 1:Si_4plus_solute
                    balance = dissolved_runoff_solute_fluxes[i_solute] / 
                        dissolving_solids_solute_fluxes[i_solute]
                        
                    println(solute_names[i_solute] * " Dsolid,Dsolute,bal ",
                        [dissolving_solids_solute_fluxes[i_solute], 
                        dissolved_runoff_solute_fluxes[i_solute],
                        balance ])
                end 
                println()
            end 
        end 
    # regolith transportation matrix
        regolith_transport_fluxes = fill(0.0, 0:n_boxes)
        #box_shallow_soil_boundary_influx = fill(0.,n_boxes)
        transport_matrix_all_boxes .= 0.
        #transport_matrix_bulk = 0 
        for i_box in 1:n_boxes # construct a no-reaction transport_matrix_all_boxes
            transport_matrix_all_boxes[i_box, i_box] += 1.0
            if enable_interbox_transport
                if i_box < n_boxes # proxy for connected to the right
                    i_neighbor_box = i_box + 1
                    transported = land_regolith_flow_parameter * # meters / year, * H -> m2 / yr 
                        (box_surface_elevations[i_box] + 
                        box_layer_thicknesses[regolith_layer,i_box] -
                        box_surface_elevations[i_neighbor_box] -
                        box_layer_thicknesses[regolith_layer, i_neighbor_box]
                        ) * # dE
                        1. / box_width * # m3 / m yr  
                        time_step * # m3 / m  
                        1. / box_width * # m3 / m2 yr; dM = D H_lim dt dE / dx2 
                        (1. / (1. - box_porosities[regolith_layer,i_box]))
                    if transported > 0. # self is higher than neighbor, controls flux 
                        # matrix indices are [equation, term in equation] 
                        transport_matrix_all_boxes[i_box, i_box] += transported
                        transport_matrix_all_boxes[i_neighbor_box, i_box] -= transported
                    else # neighbor is higher 
                        transport_matrix_all_boxes[i_box, i_neighbor_box] += transported
                        transport_matrix_all_boxes[i_neighbor_box, i_neighbor_box] -= transported
                    end 
                end 
            end
        end
        if enable_right_coastal_runoff
            transport_matrix_all_boxes[n_boxes, n_boxes] +=
                land_regolith_flow_parameter * # m / yr 
                coastal_abrasion_factor * 
                (box_surface_elevations[n_boxes] + box_layer_thicknesses[regolith_layer,n_boxes]) * # dE 
                1. / box_width * # m3 solid / m interface yr 
                time_step * # m3 solid / m interface 
                1.0 / box_width * # m3 / m2 yr; dM = D H_lim dt dE / dx2 
                (1.0 / (1.0 - box_porosities[regolith_layer, n_boxes]))# meters 
        end 
        if enable_left_coastal_runoff
            transport_matrix_all_boxes[1, 1] +=
                land_regolith_flow_parameter * 
                coastal_abrasion_factor * # m / yr 
                (box_surface_elevations[1] + box_layer_thicknesses[regolith_layer,1]) * # dE 
                1. / box_width * # m3 solid / m interface yr 
                time_step * # m3 solid / m interface 
                1.0 / box_width * # m3 / m2 yr; dM = D H_lim dt dE / dx2 
                (1.0 / (1.0 - box_porosities[regolith_layer, 1]))# meters 
                # * mineral_solid_volumes gives change in mineral_solid_volumes 
        end 
    # baby step era
        #if baby_step_era
            #println("baby step")
            # calculate fluxes that were imposed based on elevations at the beginning of the time step 
            # for consistency 
        if enable_interbox_transport
            for i_box in 1:n_boxes-1 # fluxes calculated at beginning because thats how they are applied in solver
                i_neighbor_box = i_box + 1
                regolith_transport_fluxes[i_box] =
                    land_regolith_flow_parameter * # meters / year, * H -> m2 / yr 
                    (box_surface_elevations[i_box] +
                    box_layer_thicknesses[regolith_layer, i_box] -
                    box_surface_elevations[i_neighbor_box] -
                    box_layer_thicknesses[regolith_layer, i_neighbor_box]
                    ) *
                    #box_layer_thicknesses[regolith_layer, i_box] * # m3 / yr 
                    1.0 / (box_width * box_width) #* # m3 / m2 yr; dM = D H_lim dt dE / dx2 
                    #(1.0 / (1.0 - box_porosities[regolith_layer, i_box]))
                    # flux of bulk, not just solid phase 
            end
        end 
        if enable_right_coastal_runoff
            regolith_transport_fluxes[n_boxes] =
                land_regolith_flow_parameter * # meters bulk / year, * H -> m2 bulk / yr 
                coastal_abrasion_factor * 
                (box_surface_elevations[n_boxes] +
                box_layer_thicknesses[regolith_layer, n_boxes]
                ) *
                1.0 / (box_width * box_width) #* # m3 / m2 yr; dM = D H_lim dt dE / dx2 
        end 
        if enable_left_coastal_runoff
            regolith_transport_fluxes[0] =
                - land_regolith_flow_parameter * # meters bulk / year, * H -> m2 bulk / yr 
                coastal_abrasion_factor * 
                (box_surface_elevations[1] +
                    box_layer_thicknesses[regolith_layer, 1]
                ) *
                1.0 / (box_width * box_width) #* # m3 / m2 yr; dM = D H_lim dt dE / dx2 
        end 

        # regolith transport all minerals considering primary mineral dissolution only
        # secondary mineral reaction fluxes will be dealt with afterward 
        saved_lu_transport_all_boxes = lu(transport_matrix_all_boxes)
        #original_box_mineral_volumes = deepcopy(box_mineral_volumes)
        for i_mineral in 1:n_minerals # add mineral-specific term, compute residual, run it
            residual_matrix = fill(0.,n_boxes_in_play)
            reacting_transport_matrix = deepcopy(transport_matrix_all_boxes)
            need_reacting_transport_matrix = false
            for i_box in 1:n_boxes # add mineral-specific reaction reaction terms to matrix 
                residual_matrix[i_box] = 
                    box_mineral_volumes[i_mineral,regolith_layer,i_box] +
                    box_erosion_rates[saprolite_layer, i_box] *
                    box_mineral_fractions[i_mineral,saprolite_layer,i_box] *
                    time_step              
                #=if i_mineral in primary_minerals && 
                    box_mineral_dissolution_solid_rate_constants[i_mineral] > 0
                    need_reacting_transport_matrix = true
                    reacting_transport_matrix[i_box, i_box] +=
                        box_mineral_dissolution_solid_rate_constants[i_mineral] *
                        time_step 
                end=#
            end
            lu_transport = deepcopy(saved_lu_transport_all_boxes)
            if need_reacting_transport_matrix 
                lu_transport = lu(reacting_transport_matrix)
            end
            new_volume_array = lu_transport \ residual_matrix
            for i_box in 1:n_boxes # get new layer thicknesses 
                box_mineral_volumes[i_mineral, regolith_layer, i_box] =
                    new_volume_array[i_box]
            end 
        end 
        for i_box in 1:n_boxes
            box_mineral_fractions[:, regolith_layer, i_box],
            box_solid_volumes[regolith_layer, i_box],
            box_layer_thicknesses[regolith_layer, i_box] =
                update_mineral_fractions(box_mineral_volumes[:, regolith_layer, i_box],
                    box_porosities[regolith_layer, i_box])
        end 
    # long step era 

        #end 
    # run phreeqc on regolith boxes 
        if enable_interpolated_phreeqc
            for i_box in boxes_in_play 
                #println("regolith box number ", i_box)
                output_block = # , n_new_total_phreeqc_runs, n_new_interpolated_phreeqc_runs = 
                phreeqc_box( 
                    box_mineral_fractions[:, regolith_layer, i_box],
                    box_layer_thicknesses[regolith_layer,i_box],
                    box_porosities[regolith_layer,i_box], 
                    box_fluid_flushing_rates[regolith_layer,i_box], 
                    box_fluid_flushing_ages[regolith_layer,i_box], 
                    box_temperatures[i_box],box_soil_CO2_values[i_box],
                    #enable_interpolated_phreeqc_lattice, # enable_default_solid_solutions,
                    time_step,1) # ,thread_number)
                #n_total_phreeqc_runs += n_new_total_phreeqc_runs
                #n_interpolated_phreeqc_runs += n_new_interpolated_phreeqc_runs
                box_pore_saturations[regolith_layer,i_box] = 
                    output_block.pore_saturation
                box_solute_concentrations[:,regolith_layer,i_box] = 
                    output_block.solute_concentrations
                box_mineral_saturation_indices[:,regolith_layer,i_box] = 
                    output_block.mineral_saturation_indices
                box_mineral_reaction_extents_from_phreeqc[:, regolith_layer, i_box] = 
                    output_block.mineral_reaction_extents
                box_mineral_porewater_reaction_rates[:,regolith_layer,i_box] = 
                    output_block.mineral_porewater_reaction_rates
                box_mineral_meters_for_residual_array[:,regolith_layer,i_box] = 
                    output_block.mineral_meters_for_residual_array
                box_CO2_uptake_bulk_rates[regolith_layer, i_box] = 
                    output_block.CO2_uptake_extent
                #box_alkalinity_production_bulk_rates[regolith_layer, i_box] =
                #    output_block.summary[2]
            end # direct, not interpolated 
        else 
            Threads.@threads for i_box in 1:n_boxes
                #println()
                #println("regolith box = ", i_box)
                thread_number = Threads.threadid()
                output_block = # , n_new_total_phreeqc_runs, n_new_interpolated_phreeqc_runs = 
                phreeqc_box( 
                    box_mineral_fractions[:, regolith_layer, i_box], 
                    box_layer_thicknesses[regolith_layer,i_box],
                    box_porosities[regolith_layer,i_box], 
                    box_fluid_flushing_rates[regolith_layer,i_box], 
                    box_fluid_flushing_ages[regolith_layer,i_box], 
                    box_temperatures[i_box], box_soil_CO2_values[i_box],
                    time_step,thread_number)
                box_pore_saturations[regolith_layer,i_box] = 
                    output_block.pore_saturation
                box_solute_concentrations[:,regolith_layer,i_box] = 
                    output_block.solute_concentrations
                box_mineral_saturation_indices[:,regolith_layer,i_box] = 
                    output_block.mineral_saturation_indices
                box_mineral_reaction_extents_from_phreeqc[:, regolith_layer, i_box] = 
                    output_block.mineral_reaction_extents
                box_mineral_porewater_reaction_rates[:,regolith_layer,i_box] = 
                    output_block.mineral_porewater_reaction_rates
                box_mineral_meters_for_residual_array[:,regolith_layer,i_box] = 
                    output_block.mineral_meters_for_residual_array
                box_CO2_uptake_bulk_rates[regolith_layer, i_box] = 
                    output_block.CO2_uptake_extent
            end 
        end 
    # react regolith secondary minerals 
        for i_mineral in 1:n_minerals # add mineral-specific term, compute residual, run it
            #transport_matrix .= 0.0
            residual_matrix = fill(0.0,n_boxes_in_play)
            for (i_pos, i_box) in enumerate(boxes_in_play)
                residual_matrix[i_pos] =
                    box_mineral_meters_for_residual_array[i_mineral, regolith_layer, i_box]
            end 
            #lu_transport = deepcopy(saved_lu_transport)
            #if need_reacting_transport_matrix
            #    lu_transport = lu(reacting_transport_matrix)
            #end
            saved_lu_transport = saved_lu_transport_all_boxes
            #if baby_step_era == false
            #    saved_lu_transport = saved_lu_transport_interior_boxes
            #end 
            new_volume_array = saved_lu_transport \ residual_matrix
            for (i_pos, i_box) in enumerate(boxes_in_play) # get new layer thicknesses 
                if box_mineral_volumes[i_mineral, regolith_layer, i_box] + new_volume_array[i_pos] < 0.
                    println("catching negative ", mineral_names[i_mineral], " in box ", i_box)
                    box_mineral_volumes[i_mineral, regolith_layer, i_box] = 1.e-6
                else
                    box_mineral_volumes[i_mineral, regolith_layer, i_box] +=
                        new_volume_array[i_pos]
                end
            end
        end
        for i_box in boxes_in_play
            #= for i_mineral in secondary_minerals 
                box_mineral_volumes[i_mineral,regolith_layer,i_box] += 
                    box_mineral_meters_for_residual_array[i_mineral,regolith_layer,i_box]
            end =# 
            box_mineral_fractions[:,regolith_layer,i_box],
            box_solid_volumes[regolith_layer, i_box],
            box_layer_thicknesses[regolith_layer, i_box] =
                update_mineral_fractions(box_mineral_volumes[:,regolith_layer,i_box],
                    box_porosities[regolith_layer,i_box])
            #= for i_mineral in 1:n_minerals # carbonate_minerals
                if box_mineral_fractions[i_mineral, regolith_layer, i_box] < 0.0
                    error("too much loss of " * mineral_names[i_mineral] * " ",  
                        i_box)
                end 
            end =#
        end  
    # mineral diagnostics  
        erosive_mineral_sources = fill(0.,n_minerals) # m3 / yr 
        suspended_runoff_mineral_sinks = fill(0.0, n_minerals)
        old_mineral_inventories = fill(0.0, n_minerals)
        mineral_inventories = fill(0.0, n_minerals)
        mineral_inventory_changes = fill(0.0, n_minerals)
        mineral_reaction_rates = fill(0.0, n_minerals)
        mineral_balances = fill(0.0, n_minerals)
        thin_soil_boundary_influx = fill(0., n_minerals)
        for i_box in boxes_in_play 
            for i_mineral in 1:n_minerals
                erosive_mineral_sources[i_mineral] += 
                    box_erosion_rates[bedrock_layer, i_box] * #  m3 solid / m2 year
                    basement_mineral_fractions[i_mineral] * # m / yr 
                    box_width * box_width # m3 solid / yr 
                for i_layer in 1:n_layers 
                    old_mineral_inventories[i_mineral] += 
                        previous_mineral_volumes[i_mineral,i_layer,i_box] * # meters solid 
                        #1. / ( 1. - box_porosity[i_layer,i_box] ) * # meters bulk
                        box_width * box_width # m3 solid 
                    mineral_inventories[i_mineral] += 
                        box_mineral_volumes[i_mineral,i_layer,i_box] *
                        box_width * box_width
                    mineral_reaction_rates[i_mineral] +=
                        #box_mineral_porewater_reaction_rates[i_mineral,i_layer,i_box] * # mol/m3 pw yr 
                        #box_porosities[i_layer,i_box] * # mol/m3 bulk yr 
                        #box_layer_thicknesses[i_layer,i_box] * # mol/m2 yr 
                        #box_width * box_width * # mol / yr 
                        #1. / mineral_mol_per_m3[i_mineral]
                        box_mineral_meters_for_residual_array[i_mineral,i_layer,i_box] * # m3 solid / step 
                        box_width * box_width *
                        1. / time_step # m3 solid / m2 year 
                end
            end     
        end 
        mineral_inventory_changes = ( mineral_inventories .- old_mineral_inventories ) ./ time_step 
        #if enable_coastal_runoff 
            #if baby_step_era 
                box_erosion_rates[regolith_layer,:] .= 0.
                for i_box in 1:n_boxes-1
                    if regolith_transport_fluxes[i_box] !== 0. 
                        hilltop_box = i_box # flow to the right default 
                        depositing_box = i_box + 1 
                        if regolith_transport_fluxes[i_box] < 0
                            # flow to the left
                            hilltop_box = i_box + 1
                            depositing_box = i_box 
                        end 
                        box_erosion_rates[regolith_layer, hilltop_box] += # want m3 solid / m2 year 
                            land_regolith_flow_parameter * # meters / year, * H -> m2 / yr 
                            (box_surface_elevations[hilltop_box] +
                            box_layer_thicknesses[regolith_layer, hilltop_box] -
                            box_surface_elevations[depositing_box] -
                            box_layer_thicknesses[regolith_layer, depositing_box]
                            ) *
                            box_layer_thicknesses[regolith_layer, hilltop_box] * # m3 bulk / m interface yr 
                            ( 1. - box_porosities[regolith_layer, hilltop_box]) * # m3 solid / m yr 
                            1.0 / ( box_width * box_width ) 
                    end 
                end 
                for i_flux in [0,n_boxes]
                    if regolith_transport_fluxes[i_flux] !== 0.
                        hilltop_box = i_flux # flow to the right default 
                        #depositing_box = i_flux + 1 
                        if regolith_transport_fluxes[i_flux] < 0
                            # flow to the left
                            hilltop_box = i_flux + 1
                            #depositing_box = i_flux 
                        end 
                        box_erosion_rates[regolith_layer, hilltop_box] += # want m3 solid / m2 year 
                            land_regolith_flow_parameter * # meters / year, * H -> m2 / yr 
                            (box_surface_elevations[hilltop_box] +
                            box_layer_thicknesses[regolith_layer, hilltop_box]
                            ) *
                            box_layer_thicknesses[regolith_layer, hilltop_box] * # m3 bulk / m interface yr 
                            ( 1. - box_porosities[regolith_layer, hilltop_box]) * # m3 solid / m yr 
                            1.0 / ( box_width * box_width ) 
                    end 
                end 
   
                for i_mineral in 1:n_minerals
                    suspended_runoff_mineral_sinks[i_mineral] = 0.
                    if enable_left_coastal_runoff 
                        suspended_runoff_mineral_sinks[i_mineral] += 
                            regolith_transport_fluxes[0] * # m3 bulk 
                            (1. / (1. - box_porosities[regolith_layer,1])) *
                            box_mineral_volumes[i_mineral, regolith_layer, 1] *
                            box_width * box_width 
                    end
                    if enable_right_coastal_runoff
                        suspended_runoff_mineral_sinks[i_mineral] += 
                            - regolith_transport_fluxes[n_boxes] * # m3 bulk 
                            (1.0 / (1.0 - box_porosities[regolith_layer, n_boxes])) *
                            box_mineral_volumes[i_mineral, regolith_layer, n_boxes] *
                            box_width * box_width
                    end 
                end
            #= else 
                for i_mineral in 1:n_minerals
                    suspended_runoff_mineral_sinks[i_mineral] = # 
                        - land_regolith_flow_parameter * # meters / year, * H -> m2 / yr 
                        regolith_transportation_depth * # dE
                        ( box_surface_elevations[n_boxes] + 
                        box_layer_thicknesses[regolith_layer,n_boxes] ) * # bulk fluxes 
                        (1.0 - box_porosities[regolith_layer, n_boxes]) *
                        box_mineral_fractions[i_mineral, regolith_layer, n_boxes]

                        #regolith_transport_fluxes[n_boxes] * # m3 solid / m 
                        #box_mineral_fractions[i_mineral, regolith_layer, n_boxes]
                end 
            end =#
        #end 
        for (i_pos,i_box) in enumerate(boxes_in_play)
            if i_pos == 1 && i_box > 1 # proxy for next to an excluded thin-soil cell 
                #i_excluded_box = i_box - 1
                for i_mineral in 1:n_minerals
                    thin_soil_boundary_influx[i_mineral] += 
                        box_shallow_soil_boundary_influx[i_box] *   
                        ( 1. - box_porosities[regolith_layer,i_box]) * 
                        box_mineral_fractions[i_mineral,regolith_layer,i_box] *
                        box_width * box_width / time_step
                end 
            end 
        end   
        for i_mineral in 1:n_minerals
            mineral_balances[i_mineral] = erosive_mineral_sources[i_mineral] +
                mineral_reaction_rates[i_mineral] +
                suspended_runoff_mineral_sinks[i_mineral] -
                mineral_inventory_changes[i_mineral]
        end 

        if enable_mineral_diagnostics
            println()
            println("Mineral Diagnostics")
            for i_mineral in 1:n_minerals 
                #if baby_step_era
                     println(mineral_names[i_mineral],
                        " Er,Rx,Rp,Chg,Bal ", [
                        erosive_mineral_sources[i_mineral],
                        mineral_reaction_rates[i_mineral],
                        suspended_runoff_mineral_sinks[i_mineral],
                        mineral_inventory_changes[i_mineral],
                        mineral_balances[i_mineral]])
                #=else 
                    balance = erosive_mineral_sources[i_mineral] +
                        thin_soil_boundary_influx[i_mineral] +
                        mineral_reaction_rates[i_mineral] +
                        suspended_runoff_mineral_sinks[i_mineral] -
                        mineral_inventory_changes[i_mineral]
                    println(mineral_names[i_mineral],
                        " Er,Bf,Rx,Rp,Chg,Bal ", [
                        erosive_mineral_sources[i_mineral],
                        thin_soil_boundary_influx[i_mineral],
                        mineral_reaction_rates[i_mineral],
                        suspended_runoff_mineral_sinks[i_mineral],
                        mineral_inventory_changes[i_mineral],
                        balance]) 
                end =#
            end 
        end 
   
    # convert minerals to moles of components 
        erosive_component_sources = fill(0.0, n_chemical_components) # mobile_layer_transport / yr 
        dissolving_solids_component_fluxes = fill(0.0, n_chemical_components)
        suspended_runoff_component_fluxes = fill(0.0, n_chemical_components)
        thin_soil_boundary_component_influxes = fill(0.0, n_chemical_components)
        component_inventories = fill(0.0, n_chemical_components)
        component_inventory_changes = fill(0.0, n_chemical_components)
        component_balances = fill(0.0, n_chemical_components)
        for i_component in 1:n_chemical_components
            for i_mineral in 1:n_minerals
                erosive_component_sources[i_component] += 
                    erosive_mineral_sources[i_mineral] * # m3 / year
                    mineral_mol_per_m3[i_mineral] * # mol mineral / year 
                    mineral_component_stoic[i_mineral,i_component] # moles of component /m2 year
                dissolving_solids_component_fluxes[i_component] +=
                    mineral_reaction_rates[i_mineral] * # m3 / year
                    mineral_mol_per_m3[i_mineral] *
                    mineral_component_stoic[i_mineral, i_component] # moles / yr 
                suspended_runoff_component_fluxes[i_component] +=
                    suspended_runoff_mineral_sinks[i_mineral] * # m3 / year
                    mineral_mol_per_m3[i_mineral] *
                    mineral_component_stoic[i_mineral, i_component]
                thin_soil_boundary_component_influxes[i_component] +=
                    thin_soil_boundary_influx[i_mineral] * 
                    mineral_mol_per_m3[i_mineral] *
                    mineral_component_stoic[i_mineral, i_component]
                component_inventories[i_component] +=
                    mineral_inventories[i_mineral] *
                    mineral_component_stoic[i_mineral, i_component]
                component_inventory_changes[i_component] +=
                    mineral_inventory_changes[i_mineral] * # m3 / year
                    mineral_mol_per_m3[i_mineral] *
                    mineral_component_stoic[i_mineral, i_component]
                component_balances[i_component] =
                    erosive_component_sources[i_component] +
                    thin_soil_boundary_component_influxes[i_component] -
                    component_inventory_changes[i_component] +
                    dissolving_solids_component_fluxes[i_component] +
                    suspended_runoff_component_fluxes[i_component]
            end 
        end 
    # component diagnostics 
        if enable_component_diagnostics 
            println(); println("Component Diagnostics")
            for i_component in 1:6 
                        #solute_balance =
                #    dissolving_solids_component_fluxes[i_component] /
                #    dissolved_runoff_component_fluxes[i_component]
                println(chemical_component_names[i_component], " E,Bf,Ds,Rp,dI,SlBal ",
                    [erosive_component_sources[i_component], 
                    thin_soil_boundary_component_influxes[i_component],
                    dissolving_solids_component_fluxes[i_component],
                    suspended_runoff_component_fluxes[i_component],
                    #dissolved_runoff_component_fluxes[i_component],
                    component_inventory_changes[i_component],
                    solid_balance])
            end  

        end 

    # lithium tracer model 
        previous_mineral_tracer_volumes = deepcopy(box_mineral_tracer_volumes)
        # primmary minerals just set to bedrock fraction
            for i_box in 1:n_boxes 
                for i_layer in 1:n_layers 
                    for i_mineral in lithium_bearing_primary_minerals
                        for i_tracer in 1:n_tracers
                            box_mineral_tracer_volumes[i_mineral,i_tracer,i_layer,i_box] =
                                box_mineral_volumes[i_mineral,i_layer,i_box] *
                                basement_mineral_tracer_volume_fractions[i_mineral,i_tracer] 
                        end 
                    end 
                end 
            end 
        # saprolite compute and apply secondary reactions 
            i_layer = saprolite_layer 
            for i_box in 1:n_boxes
                box_solute_concentrations[Li_solute,i_layer,i_box],
                    box_solute_concentrations[Li7_solute,i_layer,i_box],
                    box_dissolved_del_Li7[i_layer,i_box],
                    box_mineral_tracer_volume_fluxes[:,Li_tracer:Li7_tracer,i_layer,i_box] = 
                    derive_lithium_porewater_isotopes(
                        box_solute_concentrations[Mg_2plus_solute,i_layer,i_box], 
                        box_fluid_flushing_rates[i_layer,i_box],
                        box_mineral_meters_for_residual_array[:,i_layer,i_box], 
                        box_mineral_tracer_volume_fractions[:,:,i_layer,i_box],
                        enable_lithium_reprecipitation, time_step)
                box_mineral_tracer_volumes[:,:,i_layer,i_box] += 
                    box_mineral_tracer_volume_fluxes[:,:,i_layer,i_box]
                for i_tracer in 1:n_tracers
                    box_mineral_tracer_volume_fractions[:,i_tracer,i_layer,i_box],
                        box_mineral_tracer_mole_fractions[:,i_tracer,i_layer,i_box] =
                        update_tracer_fractions(
                            box_mineral_volumes[:,i_layer,i_box],
                            box_mineral_tracer_volumes[:,i_tracer,i_layer,i_box])
                end 
            end 
        # saprolite erosion to regolith
            for i_box in 1:n_boxes
                for i_mineral in lithium_bearing_secondary_minerals
                    for i_tracer in [Li_tracer,Li7_tracer]
                        sap_erosion_flux = # moves existing tracer but not newly precipitated 
                            box_erosion_rates[saprolite_layer, i_box] * # meters / year 
                            box_mineral_fractions[i_mineral,saprolite_layer,i_box] *
                            box_mineral_tracer_volume_fractions[i_mineral,i_tracer,saprolite_layer,i_box] *
                            time_step 
                        box_mineral_tracer_volumes[i_mineral,i_tracer,saprolite_layer,i_box] -=
                            sap_erosion_flux
                        box_mineral_tracer_volumes[i_mineral,i_tracer,regolith_layer,i_box] +=
                            sap_erosion_flux
                    end 
                end 
            end 
        # regolith transport without reaction 
            for i_mineral in lithium_bearing_secondary_minerals
                for i_tracer in [Li_tracer, Li7_tracer]
                    residual_matrix = fill(0.0,n_boxes)
                    for i_box in 1:n_boxes 
                        residual_matrix[i_box] =
                            box_mineral_tracer_volumes[i_mineral, i_tracer, regolith_layer, i_box] #+
                    end 
                    new_thickness_array = saved_lu_transport_all_boxes \ residual_matrix
                    for i_box in 1:n_boxes
                        box_mineral_tracer_volumes[i_mineral, i_tracer, regolith_layer, i_box] =
                            new_thickness_array[i_box]
                    end # i_box 
                end # i_tracer 
            end # i_mineral 
        # update fractions 
            for i_box in 1:n_boxes
                for i_layer in 1:n_layers 
                    for i_tracer in 1:n_tracers
                        box_mineral_tracer_volume_fractions[:,i_tracer,i_layer,i_box],
                            box_mineral_tracer_mole_fractions[:,i_tracer,i_layer,i_box] =
                            update_tracer_fractions(
                                box_mineral_volumes[:,i_layer,i_box],
                                box_mineral_tracer_volumes[:,i_tracer,i_layer,i_box])
                    end 
                end 
            end 
        # regolith compute and apply secondary reactions 
            i_layer = regolith_layer 
            for i_box in 1:n_boxes
                box_solute_concentrations[Li_solute,i_layer,i_box],
                    box_solute_concentrations[Li7_solute,i_layer,i_box],
                    box_dissolved_del_Li7[i_layer,i_box],
                    box_mineral_tracer_volume_fluxes[:,:,i_layer,i_box] = 
                    derive_lithium_porewater_isotopes(
                        box_solute_concentrations[Mg_2plus_solute,i_layer,i_box], 
                        box_fluid_flushing_rates[i_layer,i_box],
                        box_mineral_meters_for_residual_array[:,i_layer,i_box], 
                        box_mineral_tracer_volume_fractions[:,:,i_layer,i_box],
                        enable_lithium_reprecipitation,time_step)
            end 
            for i_mineral in lithium_bearing_secondary_minerals # add mineral-specific term, compute residual, run it
                for i_tracer = 1:n_tracers 
                    residual_matrix = fill(0.0,n_boxes)
                    for i_box in 1:n_boxes # add mineral-specific reaction reaction terms to matrix 
                        residual_matrix[i_box] =
                            box_mineral_tracer_volume_fluxes[i_mineral, i_tracer, i_layer, i_box]
                    end
                    new_volume_array = saved_lu_transport_all_boxes \ residual_matrix
                    for i_box in 1:n_boxes # get new layer thicknesses 
                        box_mineral_tracer_volumes[i_mineral, i_tracer, i_layer, i_box] +=
                            new_volume_array[i_box]
                    end
                    for i_box in 1:n_boxes
                        box_mineral_tracer_volume_fractions[:,i_tracer,i_layer,i_box],
                        box_mineral_tracer_mole_fractions[:,i_tracer,i_layer,i_box] =
                        update_tracer_fractions(
                            box_mineral_volumes[:,i_layer,i_box],
                            box_mineral_tracer_volumes[:,i_tracer,i_layer,i_box])
                    end 
                end 
            end 
        #= react regolith secondary lithium
            for i_box in 1:n_boxes
                for i_mineral in lithium_bearing_secondary_minerals 

                    box_mineral_volumes[i_mineral,saprolite_layer,i_box] += 
                        box_mineral_meters_for_residual_array[i_mineral, saprolite_layer, i_box]
                end
            end 
            for i_box in 1:n_boxes
                for i_tracer in 1:n_tracers
                    box_mineral_tracer_volume_fractions[:,i_tracer,regolith_layer,i_box],
                        box_mineral_tracer_mole_fractions[:,i_tracer,regolith_layer,i_box] =
                        update_tracer_fractions(
                            box_mineral_volumes[:,regolith_layer,i_box],
                            box_mineral_tracer_volumes[:,i_tracer,regolith_layer,i_box])
                end 
            end =# 
        # tracer mineral-bound tracer diagnostics
            lithium_dissolution_fluxes = [0.,0.]
            lithium_runoff_fluxes = [0.,0.]
            lithium_runoff_fractions = [0.,0.]
            for i_tracer in 1:n_tracers 
                bedrock_erosive_mineral_tracer_sources = fill(0.,n_minerals) # m3 / yr 
                #saprolite_erosive_mineral_tracer_sinks = fill(0.,n_minerals)
                suspended_runoff_mineral_tracer_sinks = fill(0.0, n_minerals)
                old_mineral_tracer_inventories = fill(0.0, n_minerals)
                new_mineral_tracer_inventories = fill(0.0, n_minerals)
                mineral_tracer_inventory_changes = fill(0.0, n_minerals)
                mineral_tracer_reaction_rates = fill(0.0, n_minerals)
                tracer_dissolution_rate = 0.
                tracer_precipitation_rate = 0.
                tracer_runoff_rate = 0. 
                for i_box in 1:n_boxes 
                    for i_mineral in lithium_bearing_minerals
                        bedrock_erosive_mineral_tracer_sources[i_mineral] += 
                                box_erosion_rates[bedrock_layer, i_box] * #  m3 solid / m2 year
                                basement_mineral_fractions[i_mineral] * # m / yr 
                                basement_mineral_tracer_volume_fractions[i_mineral,i_tracer] * 
                                box_width * box_width # m3 solid / yr 
                        for i_layer in 1:n_layers
                            old_mineral_tracer_inventories[i_mineral] += 
                                previous_mineral_tracer_volumes[i_mineral,i_tracer,i_layer,i_box] * # meters solid 
                                box_width * box_width # m3 solid 
                            new_mineral_tracer_inventories[i_mineral] += 
                                box_mineral_tracer_volumes[i_mineral,i_tracer,i_layer,i_box] *
                                box_width * box_width
                            if box_mineral_meters_for_residual_array[i_mineral, i_layer, i_box] < 0 # dissolving src  
                                rate = 
                                    box_mineral_meters_for_residual_array[i_mineral,i_layer,i_box] *
                                    box_mineral_tracer_volume_fractions[i_mineral, i_tracer, i_layer, i_box] *
                                    box_width * box_width / time_step # m3 / year  
                                mineral_tracer_reaction_rates[i_mineral] += rate
                                tracer_dissolution_rate -= rate * # m3 / step 
                                    Li_mol_per_m3 # mol Li / year 
                                if i_tracer == 1
                                    for i_dir in 1:2
                                        lithium_dissolution_fluxes[i_dir] -= rate * # m3 / step 
                                        Li_mol_per_m3 *
                                        river_watershed_fractions[i_dir,i_box]
                                    end 
                                end 
                            else # precipitating mineral sink for solute Li 
                                rate = 
                                    box_mineral_tracer_volume_fluxes[i_mineral,i_tracer,i_layer,i_box] *
                                    box_width * box_width / time_step # m3 / year 
                                mineral_tracer_reaction_rates[i_mineral] += rate
                                tracer_precipitation_rate += rate * # m3 Li / year 
                                    Li_mol_per_m3 # mol Li / year 
                                #end 
                            end
                        end 
                    end    
                    i_solute = Li_solute + i_tracer - 1
                    for i_layer in 1:n_layers
                        flux = 
                            box_solute_concentrations[i_solute,i_layer,i_box] * # mol / m3 pw
                            box_fluid_flushing_rates[i_layer,i_box] * # m3 / m2 yr -> mol / m2 yr
                            box_width * box_width # mol Li / year 
                        tracer_runoff_rate += flux
                        if i_tracer == 1
                            for i_dir in 1:2
                                lithium_runoff_fluxes[i_dir] +=
                                    flux *
                                    river_watershed_fractions[i_dir,i_box]
                            end 
                        end 
                    end 
                end 
                lithium_runoff_fractions = lithium_runoff_fluxes ./ lithium_dissolution_fluxes
                mineral_tracer_inventory_changes = ( new_mineral_tracer_inventories .- old_mineral_tracer_inventories ) ./ time_step 
                if enable_left_coastal_runoff
                    for i_mineral in 1:n_minerals
                        flux = # 
                            land_regolith_flow_parameter * 
                            box_surface_elevations[1] *
                            #1. / box_width * 
                            box_mineral_volumes[i_mineral,1,1] * #* # m2 / m yr 
                            box_mineral_tracer_volume_fractions[i_mineral, i_tracer, i_layer, 1]
                            #1. / box_width 
                        suspended_runoff_mineral_tracer_sinks[i_mineral] = flux
                        #lithium_runoff_fluxes[1] += flux
                    end
                end
                if enable_right_coastal_runoff
                    for i_mineral in 1:n_minerals
                        
                        flux = 
                            land_regolith_flow_parameter *
                            box_surface_elevations[n_boxes] *
                            #1. / box_width * 
                            box_mineral_volumes[i_mineral, 1, n_boxes] * #* # m2 / m yr 
                            box_mineral_tracer_volume_fractions[i_mineral, i_tracer, i_layer, n_boxes]
                        #1. / box_width 
                        suspended_runoff_mineral_tracer_sinks[i_mineral] = flux
                        #lithium_runoff_fluxes[2] += flux
                    end
                end
                
                if enable_tracer_diagnostics 
                    println()
                    println(tracer_names[i_tracer], " Mineral Diagnostics")
                    for i_mineral in lithium_bearing_minerals 
                        balance = bedrock_erosive_mineral_tracer_sources[i_mineral] +
                            #saprolite_erosive_mineral_tracer_sinks[i_mineral] +
                            mineral_tracer_reaction_rates[i_mineral] -
                            suspended_runoff_mineral_tracer_sinks[i_mineral] -
                            mineral_tracer_inventory_changes[i_mineral]

                        println(mineral_names[i_mineral] * "-Li BEr,Rx,Rp,Chg,Bal ", [
                            bedrock_erosive_mineral_tracer_sources[i_mineral],
                            #saprolite_erosive_mineral_tracer_sinks[i_mineral],
                            mineral_tracer_reaction_rates[i_mineral],
                            suspended_runoff_mineral_tracer_sinks[i_mineral],
                            mineral_tracer_inventory_changes[i_mineral],
                            balance])
            
                        tracer_balance = tracer_dissolution_rate - 
                            tracer_precipitation_rate -
                            tracer_runoff_rate
                        println(tracer_names[i_tracer]," Diss,Pcp,Run,Bal ",
                            [tracer_dissolution_rate,
                            tracer_precipitation_rate,
                            tracer_runoff_rate,
                            tracer_balance])
                    end 
                end 
                    #println()
            end 

    # Erosional age
        residual_matrix = fill(0.,n_boxes)
        #reacting_transport_matrix = deepcopy(transport_matrix)
        for i_box in 1:n_boxes
            residual_matrix[i_box] = 
                ( box_erosional_age[regolith_layer,i_box] + time_step ) *
                previous_layer_thicknesses[regolith_layer, i_box]
        end 
        new_age_array = saved_lu_transport_all_boxes \ residual_matrix
        for i_box in 1:n_boxes # get new layer thicknesses 
            box_erosional_age[regolith_layer, i_box] =
                new_age_array[i_box] / 
                box_layer_thicknesses[regolith_layer, i_box]
        end 

    # river concentrations 
        river_concentrations_from_phreeqc = fill(0.0, 2, n_solutes)
        dissolved_runoff_solute_fluxes = fill(0.0, 2, n_solutes)
        dissolved_runoff_component_fluxes = fill(0.0, 2, n_chemical_components) # mol/yr 
        dissolving_solids_solute_fluxes = fill(0.0, n_solutes)
        river_concentrations_from_solids = fill(0.0, n_solutes)
        river_del_Li7_values = fill(0.,2) 
        river_flows = fill(0.,2)
      
        box_flow_volumes = fill(0.,n_layers,n_boxes)
        for i_box in 1:n_boxes
            for i_dir in 1:2 # add direct runoff first 
                river_flows[i_dir] += precip_minus_evap * (1.0 - groundwater_fraction) *
                    river_watershed_fractions[i_dir,i_box] * # m3 / m2 yr 
                    box_width * box_width 
            end 
            for i_layer in 1:n_layers
                box_flow_volumes[i_layer,i_box] =
                    box_fluid_flushing_rates[i_layer, i_box] *
                    box_width * box_width # m3 water / yr
                for i_dir in 1:2
                    river_flows[i_dir] += box_flow_volumes[i_layer,i_box] * 
                        river_watershed_fractions[i_dir,i_box]
                    for i_solute in 1:n_solutes
                        dissolved_runoff_solute_fluxes[i_dir,i_solute] +=
                            box_solute_concentrations[i_solute, i_layer, i_box] * # mol / m3 pw
                            box_flow_volumes[i_layer,i_box] * # mol / yr 
                            river_watershed_fractions[i_dir,i_box]
                    end
                    river_del_Li7_values[i_dir] += box_dissolved_del_Li7[i_layer, i_box] * 
                        box_flow_volumes[i_layer,i_box] * 
                        river_watershed_fractions[i_dir,i_box]
                end
            end
        end
        for i_dir in 1:2
            if river_flows[i_dir] > 0.
                river_del_Li7_values[i_dir,:] /= river_flows[i_dir]
            end
        end
        for i_component in 1:n_solute_aligned_components
            i_solute = i_component
            dissolved_runoff_component_fluxes[:,i_component] =
                dissolved_runoff_solute_fluxes[:,i_solute] ./
                solute_component_stoiciometry[i_solute]
        end
        river_temperatures = [ 
            box_temperatures[1], box_temperatures[end]
        ]
        for i_dir in 1:2
            river_concentrations_from_phreeqc[i_dir, 1:pCO2_solute-1] = 
                dissolved_runoff_solute_fluxes[i_dir, 1:pCO2_solute-1] ./
                river_flows[i_dir] # mol / m3 
            pCO2, pH = 
                solution_pCO2(
                    river_concentrations_from_phreeqc[i_dir, Alk_solute], 
                    river_concentrations_from_phreeqc[i_dir, Tot_CO2_solute],
                    river_temperatures[i_dir])
            river_concentrations_from_phreeqc[i_dir,pCO2_solute] = pCO2
            river_concentrations_from_phreeqc[i_dir,pH_solute] = pH
        end 
        
        dissolving_solids_solute_fluxes = # [n_solutes] 
            - dissolving_solids_component_fluxes[CaO_component:SiO2_component] .*
            solute_component_stoiciometry[CaO_component:SiO2_component]
        river_concentrations_from_solids[Ca_2plus_solute:Si_4plus_solute] =
            dissolving_solids_solute_fluxes ./ 
            (river_flows[1] + river_flows[2])
        #dissolving_solids_component_fluxes[CaO_component:SiO2_component] .* 
        #solute_component_stoiciometry[i_solute] ./ river_flow
    # Be-10 
        box_Be_fluxes = fill(0.,Be9_solute:Be10_solute,n_boxes)
        #Be_9, Be_10 = 1, 2 
        for i_box in 1:n_boxes
            for i_layer in 1:n_layers 
                for i_mineral in primary_minerals
                    box_Be_fluxes[Be9_solute, i_box] +=
                        # moles / m3 porewater year 
                        - box_mineral_porewater_reaction_rates[i_mineral,i_layer,i_box] * # mol / m3 pw yr 
                        mineral_molwts[i_mineral] * # grams mineral / m3 pw yr 
                        box_porosities[i_layer,i_box] * # g min / m3 bulk yr 
                        box_layer_thicknesses[i_layer,i_box] * # g mineral / m2 yr
                        2.5e-6 / # g Be9 / m2 bulk yr, from Wittmann 2015 
                        9. * # mol Be9 / m2 yr 
                        box_width * box_width # moles / year 
                end 
            end 
        end 
        rainfall_Be10_coefficient = 1.22E6 / # atoms / cm2 yr at Berkeley Monaghan 
            0.625 * # rainfall at Berkeley m / yr 
            1.e4 * # atoms / m2 year 
            box_width * box_width * # atoms / box year
            1 / 6.023e23 # moles / box yr 
        for i_box in 1:n_boxes
            box_Be_fluxes[Be10_solute,i_box] = 
                rainfall_Be10_coefficient * 
                box_fluid_flushing_rates[regolith_layer, i_box] # moles / year
        end 
        for i_box in 1:n_boxes
            for i_Be = Be9_solute:Be10_solute
                box_solute_concentrations[i_Be,regolith_layer,i_box] = 
                    box_Be_fluxes[i_Be, i_box] * # moles / yr 
                    1.0 / box_flow_volumes[i_layer, i_box]
                for i_dir in 1:2
                    river_concentrations_from_phreeqc[i_dir, i_Be] +=
                        box_Be_fluxes[i_Be, i_box] * # moles / yr 
                        #1. / box_flow_volumes[i_layer,i_box] * # y / m3 -> mol / m3  
                        river_watershed_fractions[i_dir,i_box] # mol / m3 
                end
            end 
        end 
        for i_dir in 1:2
            river_concentrations_from_phreeqc[i_dir, Be9_solute:Be10_solute] ./=
                river_flows 
        end 

        #= reacting_transport_matrix = deepcopy(transport_matrix_all_boxes)
        residual_matrix_Be10 = fill(0.,n_boxes)
        rainfall_Be10_coefficient = 1.22E-6 / # E12 atoms / cm2 yr at Berkeley Monaghan 
            0.625 * # scale for rainfall
            1.e4   # atoms / m2 year 
        for i_box in 1:n_boxes
            box_Be_10[Be_10_influx,i_box] = 
                rainfall_Be10_coefficient * 
                box_fluid_flushing_rates[regolith_layer, i_box] # E12 atoms / m2 year
            residual_matrix_Be10[i_box] =
                box_Be_10[Be_10_inventory, i_box] + # E12 atoms / m3 solid 
                box_Be_10[Be_10_influx,i_box] * # E12 atoms / m2 yr
                #1. / box_solid_volumes[regolith_layer,i_box] * # E12 atoms / m3 solid yr 
                time_step # E12 atoms / m2 
            reacting_transport_matrix[i_box,i_box] += 4.6e-7 * time_step # 1.5myr half life
        end
        new_Be10_array = lu(reacting_transport_matrix) \ residual_matrix_Be10
        for i_box in 1:n_boxes # get new layer thicknesses 
            box_Be_10[Be_10_inventory,i_box] =
                new_Be10_array[i_box]
            box_Be_10[Be_10_concentration,i_box] =
                box_Be_10[Be_10_inventory, i_box] / # E12 atoms 
                box_solid_volumes[regolith_layer,i_box] # E12 atoms / m3 solid 
        end =#

  
        # river chemistry diagnostics 
        if enable_river_chemistry_diagnostics
            println()
            println("Solutes")
            for i_solute in 1:5
                println(solute_names[i_solute] * " flux Ds Sl, conc Ds Sl ", [
                    dissolved_runoff_solute_fluxes[i_solute],
                    dissolving_solids_solute_fluxes[1, i_solute] + dissolving_solids_solute_fluxes[2, i_solute],
                    river_concentrations_from_phreeqc[2,i_solute],
                    river_concentrations_from_solids[i_solute]])
            end
        end 

    # wrap-up output 
       
     
        #println("saprolite ", box_layer_thicknesses[saprolite_layer, :])
        #= if enable_primary_mineral_dissolution &&
            mineral_dissolution_solid_rate_constants[Anorthite_mineral] > 0.

            println("high regolith primary " *
                base_cation_abundance_report(
                    primary_minerals, box_mineral_volumes[:, 1, 1]) * " secondary " * 
                base_cation_abundance_report(
                    secondary_minerals, box_mineral_volumes[:, 1, 1]) * " " *
                clay_maturity_report(box_mineral_fractions[:, 1, 1]))
            #= println("low regolith  " *
                    string(box_layer_thicknesses[1, 2]) * " primary " *
                    base_cation_abundance_report(
                        primary_minerals, box_mineral_volumes[:, 1, 2]) * " secondary " *
                    base_cation_abundance_report(
                        secondary_minerals, box_mineral_volumes[:, 1, 2])) =# 
            println("saprolite     primary " *
                base_cation_abundance_report(
                    primary_minerals, box_mineral_volumes[:, 2, 1]) * " secondary " *
                base_cation_abundance_report(
                    secondary_minerals, box_mineral_volumes[:, 2, 1]) * " " *
                clay_maturity_report(box_mineral_fractions[:, 2, 1]))
            #print_table(solute_names, river_concentrations, typical_river_solute_concentrations)
            println("river ca ", river_concentrations_from_phreeqc[:,Ca_2plus_solute])
            println("river d7Li ", river_del_Li7_values)
            #println("soil Be10 ", box_Be_10[Be_10_concentration,:])
            println("age ", box_erosional_age[regolith_layer,1])
        end =# 
        # returns 
        timepoints[i_step+1] = time_now + time_step
        parameter_timeseries[:, i_step+1] = [
            temperature_offset, atm_CO2, 
            soil_CO2_factor, 
            precip_minus_evap, 
            basement_mineral_basalt_fraction, basement_calcite_fraction,
            basement_dolomite_fraction, 
            saprolite_thickness, aquifer_recharge_fraction,
            saprolite_erosion_constant, saprolite_erosion_efold_depth,
            land_regolith_flow_parameter, mineral_reaction_rate_parameter, 
            n_floodplain_boxes]

        layer_thickness_timeseries[:, :, i_step+1] = box_layer_thicknesses
        mineral_fraction_timeseries[:, regolith_layer:saprolite_layer, :, i_step+1] = 
            box_mineral_fractions[:, regolith_layer:saprolite_layer, :]
        mineral_reaction_rates_timeseries[:, :, :, i_step+1] = box_mineral_porewater_reaction_rates
        mineral_saturation_state_timeseries[:, :, :, i_step+1] = box_mineral_saturation_indices
        box_solute_timeseries[:, :, :, i_step+1] = box_solute_concentrations
        box_soil_CO2_timeseries[:, i_step+1] = box_soil_CO2_values
        box_temperature_timeseries[:,i_step+1] = box_temperatures
        box_precipitation_rate_timeseries[:,i_step+1] = box_precipitation_rates
        river_solute_timeseries[:, :, i_step+1] = river_concentrations_from_phreeqc
        river_solute_timeseries[:,Li7_solute, i_step+1] = river_del_Li7_values
        regolith_age_timeseries[:, i_step+1] = box_erosional_age[regolith_layer, :]
        erosion_source_timeseries[:, :, i_step+1] = box_erosion_rates[:, :]
        CO2_uptake_rate_timeseries[:, :, i_step+1] = box_CO2_uptake_bulk_rates
        #alkalinity_production_bulk_rate_timeseries[:, :, i_step+1] = box_alkalinity_production_bulk_rates

        mineral_summary_timeseries[:, 1, i_step+1] =
            erosive_mineral_sources
        mineral_summary_timeseries[:, 2, i_step+1] =
            mineral_reaction_rates
        mineral_summary_timeseries[:, 3, i_step+1] =
            suspended_runoff_mineral_sinks
        mineral_summary_timeseries[:, 4, i_step+1] =
            mineral_inventories
        mineral_summary_timeseries[:, 5, i_step+1] =
            mineral_inventory_changes
        mineral_summary_timeseries[:, 6, i_step+1] =
            mineral_balances

        # erosion, reaction, runoff, inventory 
        component_summary_timeseries[:, 1, i_step+1] =
            erosive_component_sources
        component_summary_timeseries[:, 2, i_step+1] =
            dissolving_solids_component_fluxes
        component_summary_timeseries[:, 3, i_step+1] =
            suspended_runoff_component_fluxes
        component_summary_timeseries[:, 4, i_step+1] =
            component_inventories
        component_summary_timeseries[:, 5, i_step+1] =
            component_inventory_changes
        component_summary_timeseries[:, 6, i_step+1] =
            component_balances

        #carbon_flux_summary_timeseries[:,i_step+1] .= 0.
        for i_box in 1:n_boxes
            tracer_summary_timeseries[1,i_step+1] += 
                CO2_uptake_rate_timeseries[regolith_layer, i_box, i_step+1]
            tracer_summary_timeseries[2,i_step+1] += 
                CO2_uptake_rate_timeseries[saprolite_layer, i_box, i_step+1]
            #tracer_summary_timeseries[3,i_step+1] += 
            #    alkalinity_production_bulk_rate_timeseries[regolith_layer, i_box, i_step+1]
            #tracer_summary_timeseries[4, i_step+1] +=
            #    alkalinity_production_bulk_rate_timeseries[saprolite_layer, i_box, i_step+1]
        end
        for ii in 1:4
            tracer_summary_timeseries[5, i_step+1] +=
                tracer_summary_timeseries[ii, i_step+1]
        end 
        tracer_summary_timeseries[6:7, i_step+1] = lithium_dissolution_fluxes
        tracer_summary_timeseries[8:9, i_step+1] = lithium_runoff_fluxes
        tracer_summary_timeseries[10:11, i_step+1] = lithium_runoff_fractions
    # stdout 
        i_output_frequency = 1 # 0
        #if n_total_phreeqc_runs - n_interpolated_phreeqc_runs < 10 && 
        if  enable_interpolated_phreeqc
            i_output_frequency = 100
        end
        if mod(i_step, i_output_frequency) == 0 || i_step == 1
            println()
            println("year        ", time_now + time_step)
            println("CO2, uptake ", [atm_CO2,tracer_summary_timeseries[5, i_step+1]])
            #println(box_solute_concentrations[Alk_solute,saprolite_layer,1])

            #println(box_mineral_reaction_extents_from_phreeqc[:,saprolite_layer,1])
            for i_box in [3,n_boxes-2]
                println("Box ",i_box,"  thick  ", box_layer_thicknesses[1:2,i_box])
                println("      age     ", box_fluid_flushing_ages[:,i_box])
                println("      flow fc ", box_fluid_flushing_rates[:,i_box] ./ 
                    (box_fluid_flushing_rates[1,i_box] + box_fluid_flushing_rates[2,i_box]) )
                println("      sat fc  ", box_pore_saturations[:,i_box] )
                println("      pCO2    ", box_solute_concentrations[pCO2_solute,1:2,i_box])
                println("      alk     ", box_solute_concentrations[Alk_solute,1:2,i_box])
                println("      TCO2    ", box_solute_concentrations[Tot_CO2_solute,1:2,i_box])
                println("      pH      ", box_solute_concentrations[pH_solute,1:2,i_box]) 
                println("      Cal Om  ", box_mineral_saturation_indices[Calcite_mineral,1:2,i_box])
                println("      Dol Om  ", box_mineral_saturation_indices[Dolomite_mineral,1:2,i_box])
                depletion_factors = fill(0.,2)
                labels = ["Csp","Cal","Dol"]
                mineral_loop = [Anorthite_mineral,Calcite_mineral,Dolomite_mineral]
                for i_loop in 1:3
                    i_mineral = mineral_loop[i_loop]
                    println("          ",labels[i_loop]," ",box_mineral_fractions[i_mineral,1:3,i_box])
                end 
                for i_loop in 1:3
                    i_mineral = mineral_loop[i_loop]
                    for i_layer in 1:2
                        depletion_factors[i_layer] = ( 
                            box_mineral_fractions[i_mineral,i_layer,i_box] - 
                            box_mineral_fractions[i_mineral,bedrock_layer,i_box] ) / 
                            box_mineral_fractions[i_mineral,bedrock_layer,i_box]
                    end 
                    println("         f",labels[i_loop]," ",depletion_factors)
                end 
            end
            println("River pCO2  ", river_concentrations_from_phreeqc[:,pCO2_solute], [typical_river_solute_concentrations[pCO2_solute]])
            println("      alk   ", river_concentrations_from_phreeqc[:,Alk_solute], [typical_river_solute_concentrations[Alk_solute]])
            println("      TCO2  ", river_concentrations_from_phreeqc[:,Tot_CO2_solute], [typical_river_solute_concentrations[Tot_CO2_solute]])
            println("      pH    ", river_concentrations_from_phreeqc[:,pH_solute], [typical_river_solute_concentrations[pH_solute]])
            println("      Ca    ", river_concentrations_from_phreeqc[:,Ca_2plus_solute], [typical_river_solute_concentrations[Ca_2plus_solute]])
            println("       K    ", river_concentrations_from_phreeqc[:,K_plus_solute], [typical_river_solute_concentrations[K_plus_solute]])
    
             #println("     pw reac ", [box_mineral_porewater_reaction_rates[Dolomite_mineral, regolith_layer, :]])
            #println("     sl reac ", [box_mineral_meters_for_residual_array[Dolomite_mineral, regolith_layer, :]])
            println(clay_maturity_report(box_mineral_fractions[:,1,1]))
            #if enable_interpolated_phreeqc #&& 
            #    #n_total_phreeqc_runs - n_interpolated_phreeqc_runs > 1
            #        println("archived runs ", size(memories_of_phreeqc)[1], 
            #            " interpolation pre-filled for ", floor(100 * n_interpolated_phreeqc_runs / n_total_phreeqc_runs), "%")
            #end 
            #println("regolith thickness ", box_layer_thicknesses[regolith_layer, :])
            #println("saprolite thickness ", box_layer_thicknesses[saprolite_layer, :])
            #println("river ca ", river_concentrations_from_phreeqc[:, Ca_2plus_solute])
            #println("river d7Li ", river_del_Li7_values)
            #println("age ", box_erosional_age[regolith_layer, :])
            #println("Ca-spar ", [box_mineral_fractions[1, :, 1]])
            #println("B/K ", [box_mineral_fractions[Beidellite_Ca_mineral, 1:2, 1] ./
            #                 (box_mineral_fractions[Beidellite_Ca_mineral, 1:2, 1] .+
            #                  box_mineral_fractions[Kaolinite_mineral, 1:2, 1])])
        end
    end 
 
    
        results_block = box_transect_model_struct( 
            timepoints, 
            parameter_timeseries,
            box_surface_elevations,
            layer_thickness_timeseries,
            mineral_fraction_timeseries, 
            mineral_reaction_rates_timeseries,
            mineral_saturation_state_timeseries,
            mineral_summary_timeseries,
            component_summary_timeseries, 
            tracer_summary_timeseries,
            box_solute_timeseries, 
            box_soil_CO2_timeseries,
            river_solute_timeseries,
            regolith_age_timeseries,
            erosion_source_timeseries,
            CO2_uptake_rate_timeseries)#, 
            #alkalinity_production_bulk_rate_timeseries)#,
            #Be_10_timeseries )
        return results_block
    end
function run_and_plot_mudscape(temperature_driver_offset,
        atm_CO2, soil_CO2_factor,
        precip_minus_evap, groundwater_fraction,
        basement_mineral_basalt_fraction, basement_calcite_fraction, basement_dolomite_fraction,
        saprolite_thickness, aquifer_recharge_fraction,
        saprolite_erosion_constant, saprolite_erosion_efold_depth,
        land_regolith_flow_parameter, mineral_reaction_rate_parameter,
        n_floodplain_boxes,
        model_time_stages)

    run_name = "amazon"
    cd("/Users/archer/Synched/papers/gridplates/mineral_solubilities")
    if run_name in readdir() 
    else
        println("creating " * run_name)
        mkdir(run_name)
    end 
    #global enable_interpolated_phreeqc_lattice = false
    @time results_direct = mudscape(
        temperature_driver_offset,
        atm_CO2, soil_CO2_factor,
        precip_minus_evap, 
        basement_mineral_basalt_fraction, basement_calcite_fraction, basement_dolomite_fraction,
        saprolite_thickness, aquifer_recharge_fraction,
        saprolite_erosion_constant, saprolite_erosion_efold_depth,
        land_regolith_flow_parameter, mineral_reaction_rate_parameter,
        n_floodplain_boxes,
        model_time_stages)
    fuckme = "direct"
    file_name = "direct.bson" # run_name * ".bson"
    BSON.@save file_name results_direct
    #BSON.@load file_name results_direct
    #enable_interpolated_phreeqc_lattice = true
    @time results_interp = mudscape(
        temperature_driver_offset,
        atm_CO2, soil_CO2_factor,
        precip_minus_evap, 
        basement_mineral_basalt_fraction, basement_calcite_fraction, basement_dolomite_fraction,
        saprolite_thickness, aquifer_recharge_fraction,
        saprolite_erosion_constant, saprolite_erosion_efold_depth,
        land_regolith_flow_parameter, mineral_reaction_rate_parameter,
        n_floodplain_boxes,
        model_time_stages)
    fuckme = "interp"
    file_name = "interp.precip.bson" # run_name * ".bson"
    BSON.@save file_name results_interp

    f = Figure(resolution=(1000, 800))
    panel_for_plots = f[1,1] = GridLayout()
    ax = Axis(panel_for_plots[1, 1], xlabel="Myr", ylabel="Atm CO2", xlabelsize=20, ylabelsize=20)
    scatterlines!(ax, results_block.timepoints[100:end-1] ./ 1.e6 .- 2.0, atm_CO2_timeseries[100:end-1], 
        markersize=0, label="Gridplates/Phreeqc Model",color=:red)
    berner_curve = deepcopy(atm_CO2_timeseries)
    for i in 1:size(berner_curve)[1]
        if results_block.timepoints[i] > 2.e6
            berner_curve[i] = 400. + 60. * (exp(-(results_block.timepoints[i]-2.e6)/4.e5))
        end 
    end 
    scatterlines!(ax, results_block.timepoints[100:end-1] ./ 1.e6 .- 2.0, berner_curve[100:end-1],
        markersize=0, label="Berner Model",color=:black)
    axislegend(ax)
    Makie.save("CO2_spike.png", f)

    img_interp = plot_layout_regolith(results_interp, 3, "Interp")
    Makie.save("april_9.interp.1E5.png", img_interp)

    img_direct = plot_layout_regolith(results_direct, 20, "Direct")
    Makie.save("april_9.direct.1E5.png", img_direct)


    file_name = "bump_base.bson" # run_name * ".bson"
    BSON.@save file_name results_block
    #BSON.@load file_name results_block
    animate_plot_regolith(results_interp,
        run_name)

    #BSON.@load file_name results_block
    #plot(timepoints,layer_thickness_timeseries[1,4,:])
    #plot(timepoints, river_solute_timeseries[1, :])
    #plot(timepoints, mineral_fraction_timeseries[1,1,2,:])
    #plot(timepoints, river_solute_timeseries[Li7_solute, :])
    #plot(layer_thickness_timeseries[1, :, end])
    #plot(regolith_age_timeseries[:, end])

 
end
function mudscape_parameter_sensitivities()

    dir_name = "april_16"
    #model_end_time = 1.e6

    cd("/Users/archer/Synched/papers/gridplates/mineral_solubilities")
    if dir_name in readdir()
    else
        println("creating " * dir_name)
        mkdir(dir_name)
    end
    cd(dir_name)

    temperature_driver_offset = 0.
    atm_CO2 = 300. 
    soil_CO2_factor = 1.25
    precip_minus_evap = 1. # m / yr 
    #groundwater_fraction = 0.25
    basement_mineral_basalt_fraction = 0.05 
    basement_calcite_fraction = 0.2
    basement_dolomite_fraction = 0.02
    saprolite_thickness = 30. 
    aquifer_recharge_fraction = 0.25 # years
    saprolite_erosion_constant = 3.e-3 # m / yr  
    saprolite_erosion_efold_depth = 3.0
    land_regolith_flow_parameter = 3.e2 # m / yr 
    mineral_reaction_rate_parameter = 1. 
    #secondary_mineral_redissolution_rate_parameter = 1.0
    n_floodplain_boxes = 10. 
    #model_end_time = 1.e7 
    base_parameters = [
        temperature_driver_offset, 
        atm_CO2, 
        soil_CO2_factor,
        precip_minus_evap, 
        basement_mineral_basalt_fraction, 
        basement_calcite_fraction, 
        basement_dolomite_fraction,
        saprolite_thickness, 
        aquifer_recharge_fraction, 
        saprolite_erosion_constant, 
        saprolite_erosion_efold_depth,
        land_regolith_flow_parameter, 
        mineral_reaction_rate_parameter,
        n_floodplain_boxes]
    parm_tags = ["dT", "atmCO2", "SoilCO2fac",
        "Pcp", 
        "Basalt", "Calcite","Dolomite",
        "SapThk", "SapFrac", 
        "Eros", "Erexp",
        "Trans", "Diss",
        "Width"]
    parm_labels = ["Temperature","Atm pCO2","Soil CO2 factor",
        "Rainfall Runoff Rate",
        "Basalt Fraction of Bedrock",
        "Calcite Fraction of Bedrock","Dolomite Fraction of Bedrock",
        "Saprolite Thickness","Deep Flow Fraction",
        "Erosion Rate Scale","Erosion Rate E-folding Depth",
        "Regolith Transport Coefficient","Primary Dissolution Rate Scale",
        "Floodplain Width" ]
    #twiddled_parameter_list = 3:11 
    twiddled_parameter_values_list = [ 
        [-2.,0.,2.], #                                   1 temperature_driver_offset
        [100.,300.,1000.,3000.], #                       2 atm_CO2
        [0.5,1.0,2.0], #                                 3 soil_CO2_factor 
        [0.1, 0.25,0.5,1.0,1.5,2.0], #                   4 precip_minus_evap
        [0.0, 0.05, 0.1], #                              5 basement_mineral_basalt_fraction
        [0.0, 0.02, 0.04, 0.08, 0.1], #                  6 basement_calcite_fraction
        [0.0, 0.02, 0.05], #                             7 basement_dolomite_fraction
        [10.,30.0,100.0], #                              8 saprolite_thickness
        [0.1,0.15,0.25,0.4], #                           9 aquifer_recharge_fraction
        [1.e-3,3.e-3,1.e-2, 3.e-2, 1.e-1, 3.e-1, 1.], # 10 saprolite_erosion_constant
        [1.,3.,10.], #                                  11 saprolite_erosion_efold_depth 
        [10.,30.,100.,300.,500.,1000.], #               12 land_regolith_flow_parameter
        [0.1,0.3,1.,3.,10.],  #                         13 mineral_reaction_rate_parameter
        [1.,5.,10.,20.] ] #                             14 floodplain boxes 

    run_tag = "april_16"

    results_block = mudscape(base_parameters, model_time_stages)
    file_name_base = run_tag * ".base"
    file_name = "results." * file_name_base * ".bson"
    BSON.@save file_name results_block
    results_block_base = deepcopy(results_block)

    #trial_number = 1
    #file_name = 
    save_phreeqc_memories(file_name_base)
    reset_phreeqc_memories()
    for i_variable_twist in 1:14 # 4:12 # [7,8,9,10,11]# 1:11
        parameter_list = deepcopy(base_parameters)
        #file_name = "results_interp." * parm_tags[i_variable_twist] *
        #    "." * string(base_parameters[i_variable_twist]) * ".bson"
        #results_block = deepcopy(results_block_base)
        #println(file_name)
        #BSON.@save file_name results_block
        for parameter_tweak in twiddled_parameter_values_list[i_variable_twist]
            file_name_base = 
                run_tag * "." * parm_tags[i_variable_twist] *
                "." * string(parameter_tweak)
            file_name = "results." * file_name_base * ".bson"
            #println();println(file_name);println()
            results_block = deepcopy(results_block_base)
            if base_parameters[i_variable_twist] == parameter_tweak 
                println("copying base run to ", file_name)
            else 
                parameter_list[i_variable_twist] = parameter_tweak
                println("running new ", file_name)
                results_block = weather(parameter_list, model_time_stages)
                #id_tag = "v2" * string(trial_number)
            end 
            BSON.@save file_name results_block 
            save_phreeqc_memories(file_name_base)
            reset_phreeqc_memories()

            #= for i_box in 1:n_boxes
                for i_mineral in 1:n_minerals
                    basement_mineral_fractions = fill(0.0, n_minerals)
                    basement_mineral_fractions[1:n_primary_minerals] =
                        basement_mineral_fraction_endmembers[granite_endmember, :] .* 
                            (1.0 - basement_mineral_basalt_fraction) +
                        basement_mineral_fraction_endmembers[basalt_endmember, :] .* 
                            basement_mineral_basalt_fraction
                    results_block.mineral_fraction_timeseries[i_mineral,3,i_box,:] .=
                        basement_mineral_fractions[i_mineral]
                end 
            end =#
            #file_name = "results_1myr_b." * parm_tags[i_variable_twist] *
            #            "." * string(parameter_tweak) * ".bson"
            #BSON.@save file_name results_block 
        end
    end 

    for i_variable_twist in 1:14# 2 # animation of end states 
        img_number = 0
        for file in readdir()
            if file[1:3] == "fra"
                rm(file,force=true)
            end 
        end 
        for parameter_tweak in twiddled_parameter_values_list[i_variable_twist]
            file_name =
                "results.april_10." * parm_tags[i_variable_twist] *
                "." * string(parameter_tweak) * ".bson"
            BSON.@load file_name results_block
            img_number += 1
            plot_file_name = "frame." * lpad(img_number, 3, "0") * ".png"
            println(plot_file_name) 
            title = parm_labels[i_variable_twist] * " " * string(parameter_tweak)
            i_step = 201

            img = plot_layout_regolith(results_block, i_step, title)
            
            #filename = "fra." * lpad(img_number, 3, "0") * ".png"
            #rm(filename,force=true)
            Makie.save(plot_file_name, img)
        end   
        mp4_file = "regolith.april_10." * parm_tags[i_variable_twist]  * ".mp4"
        rm(mp4_file, force=true)
        run(`/usr/local/bin/ffmpeg -r 2 -f image2 -s 1920x1080 -i frame.%03d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p $mp4_file`)
    end       

    #for extraction_function in [tracer_5() ]
   
end 
function mudscape_parameter_short_term_sensitivities()
    restart_file_name = "results.bson"
    BSON.@load file_name results_block
    time_point = size(results_block.time_points)[1]
    temperature_driver_offset = 0.
    atm_CO2 = 300. 
    soil_CO2_factor = 1.25
    precip_minus_evap = 1. # m / yr 
    acid_added = 0.
    acidified_box = 0
    base_parameters = [
        temperature_driver_offset, 
        atm_CO2, 
        soil_CO2_factor,
        precip_minus_evap, 
        aquifer_recharge_fraction, 
        saprolite_erosion_constant, 
        saprolite_erosion_efold_depth, 
        land_regolith_flow_parameter, 
        mineral_reaction_rate_parameter, 
        acid_added, 
        acidified_box ]
    parm_tags = ["dT", 
        "atmCO2", 
        "SoilCO2fac",
        "P_minus_E", 
        "SapFrac", 
        "Eros", 
        "Erexp",
        "Trans", 
        "Diss",
        "Acid_added",
        "Acidified_box"
        ]
    parm_labels = ["Temperature",
        "Atm pCO2",
        "Soil CO2 factor",
        "Rainfall Runoff Rate",
        "Deep Flow Fraction",
        "Erosion Rate Scale",
        "Erosion Rate E-folding Depth",
        "Regolith Transport Coefficient",
        "Primary Dissolution Rate Scale",
        "Acid flux",
        "Acidified box" ]
    #twiddled_parameter_list = 3:11 
    twiddled_parameter_values_list = [ 
        [-2.,0.,2.], #                                   1 temperature_driver_offset
        [100.,300.,1000.,3000.], #                       2 atm_CO2
        [0.5,1.0,2.0], #                                 3 soil_CO2_factor 
        [0.1, 0.25,0.5,1.0,1.5,2.0], #                   4 precip_minus_evap
        [0.1,0.15,0.25,0.4], #                           9 aquifer_recharge_fraction
        [1.e-3,3.e-3,1.e-2, 3.e-2, 1.e-1, 3.e-1, 1.], # 10 saprolite_erosion_constant
        [1.,3.,10.], #                                  11 saprolite_erosion_efold_depth 
        [10.,30.,100.,300.,500.,1000.], #               12 land_regolith_flow_parameter
        [0.1,0.3,1.,3.,10.],  #                         13 mineral_reaction_rate_parameter
        [1.,5.,10.,20.] ] #                             14 floodplain boxes 



function plot_mudscape_parameter_sensitivies(external_function, i_step)
    f = Figure(resolution=(1800, 1200))
    panel_for_title = f[1, 1] = GridLayout()

    panel_for_plots = f[2, 1] = GridLayout()

    n_columns = 3
    i_row, i_column, i_plot = 1, 1, 1
    axes = []
    variable_names = ["blah"]
    plot_name = "blah"
    colors = [:black,:red,:orange,:blue]
    n_values = 1
    for i_variable_twist in 1:14 # [1,3,4,5,7,8,9,10,11,12,13,14]
        println(i_variable_twist)
        ax = Axis(panel_for_plots[i_row, i_column],
            xlabel=parm_labels[i_variable_twist],
            xlabelsize=20)
        push!(axes, ax)
        #linkyaxes!(ax,axes[1])
        if i_column < n_columns
            i_column += 1
        else
            i_row += 1
            i_column = 1
        end
        i_plot += 1
        results_blocks = []
        for parameter_tweak in twiddled_parameter_values_list[i_variable_twist]
            file_name =
                "results.april_10." * parm_tags[i_variable_twist] *
                "." * string(parameter_tweak) * ".bson"
            BSON.@load file_name results_block
            push!(results_blocks, results_block)
        end
        n_model_runs = size(results_blocks)[1]
        test1, plot_name, variable_names = external_function(results_blocks[1],i_step)
        n_values = size(variable_names)[1]
        if n_values == 1
            i = 1
            output_parameter_response = fill(0.0, n_model_runs)
            for (i_run, results_block) in enumerate(results_blocks)
                output_parameter_response[i_run],plot_name, variable_names =
                    external_function(results_block,i_step)
            end
            Makie.scatterlines!(ax,
                twiddled_parameter_values_list[i_variable_twist],
                output_parameter_response[:],label=variable_names[i],
                color=colors[i],markercolor=colors[i])
        else 
            output_parameter_response = fill(0.0, n_values, n_model_runs)
            for (i_run, results_block) in enumerate(results_blocks)
                output_parameter_response[:,i_run], plot_name, variable_names =
                    external_function(results_block,i_step)
            end
            for i in 1:n_values
                Makie.scatterlines!(ax,
                    twiddled_parameter_values_list[i_variable_twist],
                    output_parameter_response[i,:],label=variable_names[i],
                    color=colors[i],markercolor=colors[i])
            end 
        end 
        #Makie.scatterlines!(ax,
        #    twiddled_parameter_values_list[i_variable_twist],
        #    output_parameter_response[2, :])
        if size(axes)[1] > 1
            linkyaxes!(axes[1], axes[end])
        end
    end
    #linkyaxes!(axes[1],axes[2],axes[3])
    Label(panel_for_title[1, 1], plot_name, tellwidth=false, textsize=30)
    if n_values > 1
        axislegend(axes[1])
    end 
    return f
    #Makie.save("carbon_fluxes.april_10.png", f)
end
function plot_mudscape_CO2_spikes()
    results_blocks = []
    file_names = ["bump_everything.bson",  #"bump_T4.bson",
        "bump_pcp_1.08.bson", "bump_LD_1.08.bson", "bump_CO2_plus_250.bson",
        "bump_rxt_1.4.bson", "bump_SE_1.08.bson",
        "bump_T4.bson"]
    variable_names = ["All feedbacks", #"T + 4 C",
        "Precipitation + 8%", "Land Transport + 8%", "CO2 + 250 ppm",
        "Dissolution Rates + 40%",  "Erosion rate + 8%",
        "T + 4 C"]
    #=file_names = ["bump_everything_but_T.bson",#"bump_T4.bson",
        "bump_pcp_1.08.bson", "bump_LD_1.08.bson",
        "bump_rxt_1.4.bson", "bump_SE_1.08.bson"]
    variable_names = ["All feedbacks",#"T + 4 C",
        "Precipitation + 8%", "Land Transport + 8%",
        "Dissolution Rates + 40%", "Erosion rate + 8%"]=#
    colors = [:red,:blue,:orange,:black,:violet,:green,:pink]
    for file_name in file_names
        BSON.@load file_name results_block 
        push!(results_blocks,results_block)
    end 

    f = Figure(resolution=(1200, 1200))
    panel_for_title = f[1,1] = GridLayout()
    #title = "Doubled CO2 Effect on Carbon Uptake Rates"
    #title = "Doubled CO2 Effect on CaO Dissolution Rates"
    title = "Doubled CO2 Effect on Carbon Fluxes"
    Label(panel_for_title[1, 1], title, tellwidth=false, textsize=30)
    panel_for_plots = f[2,1] = GridLayout()

    ax1 = Axis(panel_for_plots[1,1],ylabel="mol / m2 yr",xlabel="kyr",title="CO2 uptake",
        xlabelfontsize=20,ylabelfontsize=20)
    timepoints = ( results_blocks[1].timepoints[10:end] .- 5.e4 ) ./ 1.e3
    scatterlines!(ax1,timepoints,
        results_blocks[1].tracer_summary_timeseries[5, 10:end],
        #results_blocks[1].component_summary_timeseries[CaO_component, 2, 10:end],
        #results_blocks[1].mineral_summary_timeseries[Anorthite_mineral,2,10:end],
        markersize=0,color=colors[1],label=variable_names[1],linewidth=2)
    for i_block in 2:size(results_blocks)[1]
        scatterlines!(ax1,timepoints,
            results_blocks[i_block].tracer_summary_timeseries[5,10:end],
            #results_blocks[i_block].component_summary_timeseries[CaO_component, 2, 10:end],
            #results_blocks[i_block].mineral_summary_timeseries[Anorthite_mineral, 2, 10:end],
            label=variable_names[i_block],color=colors[i_block],
            markersize=0,linewidth=2)
    end
    #scatterlines!(ax1, [0.], [results_blocks[i_block].tracer_summary_timeseries[5, 10] * 1.4],
    #    markersize = 20)
    #axislegend(ax1,fontsize=20)

    ax2 = Axis(panel_for_plots[1, 2], ylabel="mol / m2 yr", 
        xlabel="kyr", title="Ca dissolution",
        xlabelfontsize=20, ylabelfontsize=20)
    scatterlines!(ax2, timepoints,
        - results_blocks[1].component_summary_timeseries[CaO_component, 2, 10:end],
        #results_blocks[1].mineral_summary_timeseries[Anorthite_mineral,2,10:end],
        markersize=0, color=colors[1], label=variable_names[1], linewidth=2)
    for i_block in 2:size(results_blocks)[1]
        scatterlines!(ax2, timepoints,
            - results_blocks[i_block].component_summary_timeseries[CaO_component, 2, 10:end],
            #results_blocks[i_block].mineral_summary_timeseries[Anorthite_mineral, 2, 10:end],
            label=variable_names[i_block], color=colors[i_block],
            markersize=0, linewidth=2)
    end
    axislegend(ax2, position=:rb)
    ax3 = Axis(panel_for_plots[2, 1], ylabel="mol / m2 yr",
        xlabel="kyr", title="Mg dissolution",
        xlabelfontsize=20, ylabelfontsize=20)
    scatterlines!(ax3, timepoints,
        - results_blocks[1].component_summary_timeseries[MgO_component, 2, 10:end],#results_blocks[1].mineral_summary_timeseries[Anorthite_mineral,2,10:end],
        markersize=0, color=colors[1], label=variable_names[1], linewidth=2)
    for i_block in 2:size(results_blocks)[1]
        scatterlines!(ax3, timepoints,
            - results_blocks[i_block].component_summary_timeseries[MgO_component, 2, 10:end],
            #results_blocks[i_block].component_summary_timeseries[CaO_component, 2, 10:end],
            #results_blocks[i_block].mineral_summary_timeseries[Anorthite_mineral, 2, 10:end],
            label=variable_names[i_block], color=colors[i_block],
            markersize=0, linewidth=2)
    end

    ax4 = Axis(panel_for_plots[2, 2], ylabel="mol / m2 yr", xlabel="kyr", title="CO2 uptake",
        xlabelfontsize=20, ylabelfontsize=20)
    scatterlines!(ax4, timepoints,
        results_blocks[1].tracer_summary_timeseries[5, 10:end],
        #results_blocks[1].component_summary_timeseries[CaO_component, 2, 10:end],
        #results_blocks[1].mineral_summary_timeseries[Anorthite_mineral,2,10:end],
        markersize=0, color=colors[1], label="Gridplates", linewidth=2)
    #for i_block in 2:size(results_blocks)[1]
    #    scatterlines!(ax1, timepoints,
    #        results_blocks[i_block].tracer_summary_timeseries[5, 10:end],
            #results_blocks[i_block].component_summary_timeseries[CaO_component, 2, 10:end],
            #results_blocks[i_block].mineral_summary_timeseries[Anorthite_mineral, 2, 10:end],
    #        label=variable_names[i_block], color=colors[i_block],
    #        markersize=0, linewidth=2)
    #end
    scatterlines!(ax4, [0.], [results_blocks[i_block].tracer_summary_timeseries[5, 10] * 1.4],
        markersize = 20, label = "Berner")
    axislegend(ax4)

    #legendobject = Legend(panel_for_plots[2, 2],ax1)
    #legendobject.tellwidth = false
    #set_theme!(Theme(fontsize=20))
    f
    Makie.save("bumped.png",f)
end 
function plot_layout_regolith(results_block,
    i_step, plot_title )
    # setup 
        f = Figure(resolution=(1800, 1500))
        panel_for_title = f[1,1] = GridLayout()
        Label(panel_for_title[1, 1], plot_title, tellwidth=false, textsize=30)
        panel_for_plots = f[2,1] = GridLayout()
        #title_panel = f[4,1:3] = GridLayout()
        
        panel_for_space_plots = panel_for_plots[1, 1] = GridLayout()
        #panel_for_space_plots = f[1, 1] = GridLayout()
        panel_for_miscellany = panel_for_plots[1, 2] = GridLayout()
        panel_for_square_plots = panel_for_miscellany[1,1] = GridLayout()
        panel_for_timeseries_plots = panel_for_miscellany[2,1] = GridLayout()
        panel_for_origin_bar_chart_and_legends = panel_for_miscellany[3,1] = GridLayout()
        #panel_for_bar_chart = panel_for_origin_bar_chart_and_legends[1,3] = GridLayout()
        #panel_for_solid_legend = panel_for_origin_bar_chart_and_legends[1,2] = GridLayout()
        #panel_for_solute_legend = panel_for_origin_bar_chart_and_legends[1,3] = GridLayout()

        ax_topography_space = Axis(panel_for_space_plots[1,1], title="Topography",
            ylabel = "meters(regolith), m/100 (topo)", xlabelsvisible = false, ticklabels=false )
        ax_minerals_space = Axis(panel_for_space_plots[2, 1], title="Regolith Mineralology",
            ylabel = "cumulative %", xlabelsvisible=false, ticklabels=false)
        ax_weathering_rates_space = Axis(panel_for_space_plots[3, 1], title="Weathering Rates",
            xlabelsvisible=false, ticklabels=false)
        #ax_carbon_fluxes_space = Axis(panel_for_space_plots[3, 1], title="Carbon Fluxes",
        #    yaxisposition=:right, rightspinecolor=:red)
        ax_erosion_rates_space = Axis(panel_for_space_plots[4, 1], title="Erosion Rates",
            xlabelsvisible=false, ticklabels=false,
            ylabel="meters / Myr")
        ax_CIA_space = Axis(panel_for_space_plots[5, 1], title="Regolith Diagnostics",
            xlabelsvisible=false, ticklabels=false,
            ylabel="Ca Depletion %")
        ax_age_space = Axis(panel_for_space_plots[5, 1], 
            xlabelsvisible=false, ticklabels=false,
            ylabel="Age, kyr", yaxisposition=:right, rightspinecolor=:red)
        ax_pCO2_space = Axis(panel_for_space_plots[6, 1], title="Porewater pCO2",
            xlabelsvisible=false, ticklabels=false,
            ylabel="ppm", xlabel="Distance, km")
        Makie.xlims!(ax_CIA_space, 0.0, 4000.0)
    linkxaxes!(ax_topography_space, ax_minerals_space,
            ax_weathering_rates_space, ax_erosion_rates_space, 
            ax_CIA_space, ax_age_space, ax_pCO2_space)

        ax_willenbring = Axis(panel_for_square_plots[1,1],title="Denudation Rates",
            #xscale = log, yscale = log, 
            xlabel = "Log10 slope m/km", ylabel = "Log10 Erosion meters / Myr",
            xticklabelsize=20, yticklabelsize=20)
        ax_solute_plot = Axis(panel_for_square_plots[1, 2], title="River Solute Concentrations",
            xlabel = "umol / L global mean observed", ylabel = "umol / L model")
        timestamp = "Elapsed time " * 
            string(results_block.timepoints[i_step]) * " years"
        ax_solute_timeseries = Axis(panel_for_timeseries_plots[1,1],title=timestamp,
            xticksvisible=false, ylabel="River Conc mol / L")
        ax_river_CO2_timeseries = Axis(panel_for_timeseries_plots[2, 1], title="River pCO2",
            xticksvisible=false, ylabel="ppm")
        ax_d7Li_timeseries = Axis(panel_for_timeseries_plots[3,1],title="d7Li and d10Be in Rivers",
            xlabel = "Kyr",ylabel = "d7Li")
        ax_d10Be_timeseries = Axis(panel_for_timeseries_plots[3, 1],
            ylabel = "d10Be", yaxisposition = :right, rightspinecolor = :red)
        #hidedecorations!(ax_d10Be_timeseries)
        ax_origins_barchart = Axis(panel_for_origin_bar_chart_and_legends[1,1],
            title="Solute Source and Sink Minerals",
            xticks=(1:2:9, solute_names[1:5]))
        #ax_solid_legend = Axis(panel_for_origin_bar_chart_and_legends[1,2],title="Solid Colors")
        #ax_solute_legend = Axis(panel_for_origin_bar_chart_and_legends[1,3],title="Solute Colors")

    linkxaxes!(ax_solute_timeseries, ax_d7Li_timeseries, ax_d10Be_timeseries, ax_river_CO2_timeseries)

        n_boxes = size(results_block.box_surface_elevations)[1]
        n_layers = 2
        x_locs = [0.]
        for i_loc in 1:size(results_block.box_surface_elevations)[1]-1
            append!(x_locs,x_locs[end]+400.)
        end 
        timestamp = string(results_block.timepoints[i_step]) * " years"
    # Topography 
        scaled_regolith_surface = 
            ( results_block.box_surface_elevations[:] .+
            results_block.layer_thickness_timeseries[1,:,i_step] ) .*
            1.e-2 
        not_to_scale_regolith_base = scaled_regolith_surface .-
            results_block.layer_thickness_timeseries[1,:,i_step]
        Makie.scatterlines!(ax_topography_space, x_locs, 
            scaled_regolith_surface,
            markersize=0, color=:black)
        Makie.scatterlines!(ax_topography_space, x_locs,
            not_to_scale_regolith_base,
            markersize=0, color=:brown)
    # Regolith Mineralogy 
        #Makie.ylims!(ax_minerals, 0.0, 100.0)
        for i_layer in 1:2
            cumulative_fractions = fill(-100.0 * (i_layer-1), n_boxes)
            next_cumulative_fractions = fill(-100.0 * (i_layer-1), n_boxes)
            #cmap = mineral_group_colormap()
            for i_mineral_group in 1:n_mineral_groups
                plot_me_as_region = false
                for i_mineral in 1:n_minerals
                    if mineral_groups[i_mineral] == i_mineral_group
                        for i_box in 1:n_boxes
                            if results_block.mineral_fraction_timeseries[i_mineral, i_layer, i_box, i_step] > 0.0
                                plot_me_as_region = true 
                                next_cumulative_fractions[i_box] +=
                                    results_block.mineral_fraction_timeseries[i_mineral, i_layer, i_box, i_step] * 100.0
                            end
                        end
                    end
                end
                
                if plot_me_as_region
                    point_list = get_stripe_coords(cumulative_fractions, next_cumulative_fractions, x_locs)
                    try
                    Makie.poly!(ax_minerals_space, point_list, color=mineral_group_colormap()[i_mineral_group],
                            label=mineral_group_names[i_mineral_group])
                    catch 
                    end 
                else
                Makie.scatterlines!(ax_minerals_space, x_locs, next_cumulative_fractions,
                        markersize=0, color=mineral_group_colormap()[i_mineral_group],
                        label=mineral_group_names[i_mineral_group])
                end
                cumulative_fractions = deepcopy(next_cumulative_fractions)
            end
        end 
    # Weathering rates 
        for i_solute in 1:5
            box_solute_dissolution_rates = fill(0.0, n_boxes)
            for i_mineral in 1:n_minerals
                for i_box in 1:n_boxes
                    for i_layer in 1:n_layers
                        box_solute_dissolution_rates[i_box] -=
                            results_block.mineral_reaction_rates_timeseries[i_mineral, i_layer, i_box, i_step] *
                            mineral_component_stoic[i_mineral, i_solute]
                    end
                end
            end
            Makie.scatterlines!(ax_weathering_rates_space, x_locs, box_solute_dissolution_rates,
                color=solute_colormap()[i_solute], label=solute_names[i_solute],
                markersize=0)
        end
        Makie.scatterlines!(ax_weathering_rates_space, x_locs,
            results_block.CO2_uptake_rate_timeseries[1, :, i_step] .+
            results_block.CO2_uptake_rate_timeseries[2, :, i_step],
            color=:black, label="CaCO3 dep")#,
            #markersize=0)
        #Makie.scatterlines!(ax_weathering_rates_space, x_locs,
            #results_block.alkalinity_production_bulk_rate_timeseries[1, :, i_step] .+
            #results_block.alkalinity_production_bulk_rate_timeseries[2, :, i_step],
        #    color=:red, label="Alk")#,
            #markersize=0)
        Makie.axislegend(ax_weathering_rates_space)
    # Erosion rates 
        erosion_colors = [:green,:brown,:darkgrey]
        for i_layer in regolith_layer:bedrock_layer
            Makie.scatterlines!(ax_erosion_rates_space,x_locs,
                results_block.erosion_source_timeseries[i_layer,:,i_step] .* 1.e6,
                label=layer_names[i_layer],
                markersize=0, color=erosion_colors[i_layer])
        end 
        Makie.axislegend(ax_erosion_rates_space)
    # Chemical Index of Alteration 
        calcium_relative_depletions = fill(0.,n_boxes)
    
        for i_box in 1:n_boxes 
            calcium_content = 0.
            bedrock_calcium_content = 0.0
            #total_content = 0.
            for i_mineral in 1:n_minerals 
                calcium_content +=
                    results_block.mineral_fraction_timeseries[i_mineral, regolith_layer, i_box, i_step] *
                    mineral_component_stoic[i_mineral,CaO_component]
                bedrock_calcium_content +=
                    results_block.mineral_fraction_timeseries[i_mineral, bedrock_layer, i_box, i_step] *
                    mineral_component_stoic[i_mineral, CaO_component]
                #total_content += results_block.mineral_fraction_timeseries[i_mineral,regolith_layer,i_box,i_step]
            end 
            calcium_relative_depletions[i_box] = (bedrock_calcium_content - calcium_content) /
                bedrock_calcium_content * 100.0
        end 
        Makie.scatterlines!(ax_CIA_space, x_locs,
            calcium_relative_depletions,
            markersize=0,label="CIA")
        Makie.scatterlines!(ax_age_space, x_locs,
            results_block.regolith_age_timeseries[:,i_step] .* 1.e-3,
            color=:red, markersize=0, linestyle="-",label="Age")
    # porewater pCO2 
        Makie.scatterlines!(ax_pCO2_space, x_locs,
            results_block.box_solute_timeseries[pCO2_solute,1, :, i_step],
            color=:black, markersize=0, label="Regolith")
        Makie.scatterlines!(ax_pCO2_space, x_locs,
            results_block.box_solute_timeseries[pCO2_solute, 2, :, i_step],
            color=:red, markersize=0, label="Saprolite")
        Makie.axislegend(ax_pCO2_space)
        #Makie.axislegend(ax_willenbring, merge=true)
                    
        #=Makie.scatterlines!(ax_diagnostics, x_locs, 
            results_block.Be_10_timeseries[Be_10_concentration,:,i_step],
            color=:red, label="Soil 10-B E-7 atoms/gram",
            markersize=0)=#

        #Makie.scatterlines!(ax_carbon_fluxes_space, x_locs,
        #        results_block.alkalinity_production_bulk_rate_timeseries[2, :, i_step],
        #        color=:blue, label="Saprolite Alk Src",
        #        markersize=0, linestyle="-")
        #Makie.axislegend(ax_carbon_fluxes_space)
    # Willenbring plot 
        box_denudation_rates = fill(0.,n_boxes)
        for i_box in 1:n_boxes
            box_denudation_rates[i_box] = log10( # m bulk / Myr 
                results_block.erosion_source_timeseries[regolith_layer,i_box,i_step] * # m3 solid / m2 year 
                1. / 0.8 * # ( 1. - box_porosities[regolith_layer,i_box]) * # m3 bulk / m2 yr 
                1.e6 )# m / Myr 
        end 
        #=for i_box in 1:n_boxes 
            if results_block.Be_10_timeseries[Be_10_concentration,i_box, i_step] > 0.0
                deposition =
                    results_block.Be_10_timeseries[Be_10_influx,i_box,i_step] * # E12 atoms / m2 year 
                    1.e-4 # E12 atoms / cm2 yr 
                concentration = 
                    results_block.Be_10_timeseries[Be_10_concentration,i_box, i_step] * # E12 atoms / m3 solid 
                    1. / rho_continent_crust # E12 atoms / gram solid
                erosion_g_cm2yr = deposition / concentration # g / cm2 yr 
                erosion_mm_kyr = erosion_g_cm2yr / rho_continent_crust * # m3 solid / m2 yr 
                    1.e6 # mm / kyr 
                box_denudation_rates_from_Be_10[i_box] = erosion_mm_kyr #  Willinbring
            end 
        end =# 
        box_slopes = fill(0.0, n_boxes) # want m / km 
        #box_width = 4.e5
        for i_box in 2:n_boxes-1
            box_slopes[i_box] = 
                    (results_block.box_surface_elevations[i_box] +
                    results_block.layer_thickness_timeseries[1, i_box, i_step] -
                    results_block.box_surface_elevations[i_box+1] - 
                    results_block.layer_thickness_timeseries[1, i_box+1, i_step]) /
                    box_width *
                    1.e3
        end 
        for i_box in [1,n_boxes]
            box_slopes[i_box] = 
                ( results_block.box_surface_elevations[i_box] + 
                results_block.layer_thickness_timeseries[1,i_box,i_step] ) /
                box_width *
                1.e3 # meters / km 
        end 
        for i_box in 1:n_boxes
            box_slopes[i_box] = log10(abs(box_slopes[i_box]))
        end 
        #Makie.scatterlines!(ax_willenbring,
        #    Float32.(box_slopes),
        #    Float32.(box_denudation_rates),
        #    color=:red, label="Model",
        #    markersize=0)
        for i_box in 1:n_boxes
            if box_slopes[i_box] > 0 && box_denudation_rates[i_box] > 0
                scatter_point = Point2f(box_slopes[i_box],
                    box_denudation_rates[i_box])
                    #println(solute_point)
                Makie.scatter!(ax_willenbring, scatter_point,color=:black,
                    label = "model") # ,
                        #color=solute_colormap()[i_solute], label=solute_names[i_solute])
            end 
        end
        #Makie.scatterlines!(ax_willenbring, box_slopes, box_denudation_rates_from_Be_10,
        #    color=:red, label="Model", marker)#,#,
            #linestyle="nothing")
        willenbring_slopes = [1.,1.5,2.,2.5,3.]
        #while willenbring_slopes[end] < 900.
        #    next_slope = willenbring_slopes[end] * sqrt(10.)
        #    append!( willenbring_slopes,next_slope )
        #end 
        willenbring_denudation_rates = []
        for slope in willenbring_slopes 
            denudation_rate = log10(11.9 * exp( 0.0065 * 10. ^ slope ) )
            append!(willenbring_denudation_rates,denudation_rate)
        end 
        #willenbring_slopes .= log10(willenbring_slopes)
        #willenbring_denudation_rates .= log10(willenbring_denudation_rates)
        Makie.scatterlines!(ax_willenbring, 
            Float32.(willenbring_slopes),
            Float32.(willenbring_denudation_rates),
            color=:black, label="Willenbring",
            markersize=0)
        Makie.axislegend(ax_willenbring, merge=true)
    # Solute scatterplot 
        Makie.scatterlines!(ax_solute_plot, [0.0, 500.0], [0.0, 500.0],
            markersize=0)
        for i_solute in 1:Li_solute
            if i_solute != Al_solute
                
                obs_concs = [typical_river_solute_concentrations[i_solute] * 1.e3,
                    typical_river_solute_concentrations[i_solute] * 1.e3]
                model_concs = [results_block.river_solute_timeseries[1, i_solute, i_step] * 1.e3,
                    results_block.river_solute_timeseries[2, i_solute, i_step] * 1.e3]
                solute_name = solute_names[i_solute]
                if i_solute == Li_solute
                    obs_concs .*= 1.e3
                    model_concs .*= 1.e3
                    solute_name *= " * 100"
                end 
                Makie.scatterlines!(ax_solute_plot,
                    obs_concs,
                    model_concs,
                    color=solute_colormap()[i_solute], 
                    markercolor=solute_colormap()[i_solute],
                    strokecolor=solute_colormap()[i_solute],
                    markersize=20,
                    strokewidth=2,
                    label = solute_name)
                solute_point = Point2f(obs_concs[1],
                    model_concs[1])
                Makie.scatter!(ax_solute_plot, solute_point,
                    color=:white, strokecolor=solute_colormap()[i_solute],
                    markersize=20,
                    strokewidth=2)

                #solute_point = Point2f(typical_river_solute_concentrations[i_solute] * 1.e3,
                #    results_block.river_solute_timeseries[2, i_solute, i_step] * 1.e3)
                #Makie.scatter!(ax_solute_plot, solute_point,
                #    color=solute_colormap()[i_solute], 
                #    strokecolor=solute_colormap()[i_solute],
                #    strokewidth=2,
                #    label="Amazon")

            end
        end
        Makie.axislegend(ax_solute_plot,merge=true)
    # Solute_timeseries
        for i_solute in 1:7
            #println(i_solute)
            solute_timeseries = 
                results_block.river_solute_timeseries[:,i_solute,:] .* 1.e3
            label = solute_names[i_solute]
            if i_solute == Li_solute
                solute_timeseries[:,:] .*= 1.e3
                label = label * "*1E3"
            end
            if i_solute == Al_solute
                solute_timeseries[:, :] .*= 1.e3
                label = label * "*1E3"
            end
            Makie.scatterlines!(ax_solute_timeseries,
                results_block.timepoints ./ 1.e3,
                solute_timeseries[1, :],
                color=solute_colormap()[i_solute], markersize=0,# label=label * " Andean",
                linestyle=".")
            Makie.scatterlines!(ax_solute_timeseries,
                results_block.timepoints ./ 1.e3,
                solute_timeseries[2, :],
                color=solute_colormap()[i_solute], markersize=0, label=label)
        end
        current_point = Point2f(results_block.timepoints[i_step] / 1.e3, 0.)
        Makie.scatter!(ax_solute_timeseries, current_point, 
            color=:black,markersize=20)
        #Makie.axislegend(ax_solute_timeseries)#, merge=true)
    # River pCO2 timeseries 
        Makie.scatterlines!(ax_river_CO2_timeseries,
            results_block.timepoints ./ 1.e3,
            results_block.river_solute_timeseries[1, pCO2_solute, :],
            color=:black, markersize=0,
            label="Andean", linestyle="-")
        Makie.scatterlines!(ax_river_CO2_timeseries,
            results_block.timepoints ./ 1.e3,
            results_block.river_solute_timeseries[2, pCO2_solute, :],
            color=:black, markersize=0,
            label="Amazon")
            #=Makie.scatterlines!(ax_solute_timeseries,
                results_block.timepoints ./ 1.e3,
                solute_timeseries[1,:],
                color=solute_colormap()[i_solute], markersize=0, label=label * " Andean",
                linestyle = ".")
            current_point = Point2f(results_block.timepoints[i_step] / 1.e3, solute_timeseries[i_step])
            Makie.scatter!(ax_solute_timeseries, current_point, 
                color=solute_colormap()[i_solute])

            
            current_point = Point2f(results_block.timepoints[i_step] / 1.e3, solute_timeseries[i_step])
            Makie.scatter!(ax_solute_timeseries, current_point, color=solute_colormap()[i_solute])=#
        #end
        #text!(ax, temp_summary, position=(-50, -30), textsize=20)
        #Makie.text!(ax_solute_timeseries,"fuckme",position=(0,0))
    # River Isotope timeseries
        Makie.scatterlines!(ax_d7Li_timeseries,
            results_block.timepoints ./ 1.e3,
            results_block.river_solute_timeseries[1,Li7_solute, :],
            color=:black, markersize=0,
            label="Andean",linestyle="-")
        Makie.scatterlines!(ax_d7Li_timeseries,
            results_block.timepoints ./ 1.e3,
            results_block.river_solute_timeseries[2, Li7_solute, :],
            color=:black, markersize=0,
            label="Amazon")
        Makie.scatterlines!(ax_d10Be_timeseries,
            results_block.timepoints ./ 1.e3,
            results_block.river_solute_timeseries[1, Be10_solute, :] ./ 
            results_block.river_solute_timeseries[1, Be9_solute, :] .*
            1.E7,
            color=:red, markersize=0,
            label="Andean", linestyle="-")
        Makie.scatterlines!(ax_d10Be_timeseries,
            results_block.timepoints ./ 1.e3,
            results_block.river_solute_timeseries[2, Be10_solute, :] ./
            results_block.river_solute_timeseries[2, Be9_solute, :] .*
            1.E7,
            color=:red, markersize=0,
            label="Amazon")
        Makie.axislegend(ax_d7Li_timeseries)
    # Solute mineral origins barchart
        solute_mineral_sources = fill(0.0, n_solutes, n_mineral_groups)
        for i_solute in 1:5
            for i_mineral in 1:n_minerals
                i_mineral_group = mineral_groups[i_mineral]
                for i_layer in 1:n_layers
                    for i_box in 1:n_boxes
                        solute_mineral_sources[i_solute, i_mineral_group] -=
                            results_block.mineral_reaction_rates_timeseries[i_mineral, i_layer, i_box, i_step] *
                            mineral_component_stoic[i_mineral, i_solute]
                    end
                end
            end
        end
        #mineral_group_colormap()
        for i_mineral_group in 1:n_mineral_groups
            #println(i_mineral_group)
            barplot!(
                ax_origins_barchart,
                [1, 3, 5, 7, 9],
                solute_mineral_sources[1:5, i_mineral_group],
                color=mineral_group_colormap()[i_mineral_group],
                label=mineral_group_names[i_mineral_group]
                #title="Stacked bars")
                #labels=([solute_names[1:5]])
            )
        end
    # Legends 
        solid_legend = Legend(panel_for_origin_bar_chart_and_legends[1, 2], 
            ax_origins_barchart, "Solid Phases")
            solid_legend.tellheight = true 
        solute_legend = Legend(panel_for_origin_bar_chart_and_legends[1, 3],
            ax_weathering_rates_space, "Solutes")
        solute_legend.tellheight = true
    f
end 
function animate_plot_regolith(results_block, 
    run_name)
    cd("/Users/archer/Synched/papers/gridplates/mineral_solubilities/run_libraries/")# * run_name)
    file_name = "results.april_10.base.bson"
    BSON.@load file_name results_block

    n_steps = size(results_block.timepoints)[1]
    #rm("/*.png",force=true)
    #rm(directory_name, force=true)
    #mkdir(directory_name)
    #cd(directory_name)

    #rm(outfile)
    year_list = []
    for i_step in 2:n_steps
        year = results_block.timepoints[i_step]
        if year < 1.e5
            push!(year_list,year)
        elseif year < 1.e6 && mod(year,1.e5) == 0.
            push!(year_list,year)
        elseif year < 1.e7 && mod(year,1.e6) == 0.
            push!(year_list, year)
        elseif mod(year, 1.e7) == 0.0
            push!(year_list, year)
        end 
    end

    img_number = 1
    for i_step in 1:n_steps
        year = results_block.timepoints[i_step]
        if year in year_list
            img_number += 1
            timestamp = "Elapsed time " *
                string(results_block.timepoints[i_step]) * " years"
            println(timestamp)
            
            img = plot_layout_regolith(results_block, i_step, timestamp)
            filename = "frame." * lpad(img_number, 3, "0") * ".png"
            rm(filename,force=true)
            Makie.save(filename,img)
        end 
        #println(i_step)
    end
    last_filename = "frame." * lpad(img_number, 3, "0") * ".png"
    first_filename = "frame." * lpad(1, 3, "0") * ".png"
    cp(last_filename, first_filename,force=true)
    #cd("..")
    mp4_file = "regolith.transect." * run_name * ".mp4"
    rm(mp4_file, force=true)
    #rm()
    run(`/usr/local/bin/ffmpeg -r 2 -f image2 -s 1920x1080 -i frame.%03d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p $mp4_file`)
    # ffmpeg -r 2 -f image2 -s 1920x1080 -i frame.%03d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p regolith.transect.april_9.mp4
end
function initialize_mudscape_from_results_file(file_name,time_point)
    BSON.@load file_name results_block
    n_boxes = size(results_block.box_surface_elevations)[1]
    n_layers = 2
    box_layer_thicknesses = results_block.layer_thickness_timeseries[:, :, i_step]
    box_mineral_fractions = fill(0.0, n_minerals, n_layers+1, n_boxes)
    box_mineral_volumes = fill(0.0, n_minerals, n_layers, n_boxes)
    for i_box in 1:n_boxes
        for i_layer in 1:n_layers
            box_mineral_fractions[:, i_layer, i_box] = 
                results_block.mineral_fraction_timeseries[:, i_layer, i_box, end]
            box_mineral_volumes[:,i_layer,i_box] = 
                update_mineral_volumes(box_mineral_fractions[:,i_layer,i_box], 
                    box_layer_thicknesses[i_layer, i_box], box_porosities[i_layer, i_box])
        end 
    end 
    return box_layer_thicknesses, box_mineral_fractions, box_mineral_volumes
end 


        

# aquifer model (phreeqc 10 boxes) and eroding hilltop (1-point mudscape)
function run_aquifer_phreeqc(aquifer_basalt_fraction, 
    aquifer_calcite_fraction, aquifer_dolomite_fraction, 
    aquifer_inert_fraction, 
    soil_CO2, surface_temperature, 
    flushing_age, 
    acid_added, 
    n_cells)
    #=aquifer_basalt_fraction = 0.05
    aquifer_calcite_fraction = 0.04
    aquifer_inert_fraction = 0. 
    aquifer_dolomite_fraction = 0.02
    porosity = 0.5; layer_thickness = 1.
    soil_CO2 = 2.5e4; surface_temperature = 20.
    flushing_time = 30.
    acid_added = 0.; n_cells = 10 =#
    #n_shifts = n_cells
    aquifer_granite_fraction = 1. - aquifer_basalt_fraction - 
        aquifer_calcite_fraction - aquifer_dolomite_fraction - 
        aquifer_inert_fraction 
    aquifer_mineral_fractions = fill(0.0, n_minerals)
    aquifer_mineral_fractions[:] =
        basement_mineral_fraction_endmembers[granite_endmember, :] .* 
            aquifer_granite_fraction +
        basement_mineral_fraction_endmembers[basalt_endmember, :] .* 
            aquifer_basalt_fraction
    aquifer_mineral_fractions[Calcite_mineral] = 
        aquifer_calcite_fraction
    aquifer_mineral_fractions[Dolomite_mineral] =
        aquifer_dolomite_fraction
    aquifer_mineral_fractions[Quartz_mineral] = 
        aquifer_inert_fraction
    aquifer_mineral_reaction_rate_constants = fill(0.,n_minerals)
    for i_mineral in 1:n_minerals  
        aquifer_mineral_reaction_rate_constants[i_mineral] = # want meters / year 
            aquifer_mineral_fractions[i_mineral] * # m3 mineral / m3 solid 
            lasaga_rates_pH5[i_mineral] * # m3 mineral /m3 solid yr 
            flushing_time * # mol / m3 pw step 
            1.e-3 *
            1. / n_cells
    end 
    generate_phreeqc_input_file( # omegas and rxns 
        surface_temperature, soil_CO2,
        aquifer_mineral_reaction_rate_constants,
        n_cells, acid_added, 1)
    run_command = "/Users/archer/Synched/papers/gridplates/mineral_solubilities/run_auto_phreeqc_1.x"
    cmd = `$run_command`
    run(cmd)
    #aquifer_out_filename = "aquifer_elbow7_CO2_open.out"
    aquifer_out_filename = "auto_1.out" # aquifer_no_elbow1_CO2_open.out"
    solute_concentrations,
        mineral_saturation_indices,
        mineral_reaction_extents, 
        CO2_uptake_extents = 
    parse_aquifer_phreeqc_file(aquifer_out_filename, n_cells)
    plot_title = "50 mM CaO"
    fig = plot_aquifer(solute_concentrations,
                mineral_saturation_indices,
                mineral_reaction_extents, 
                CO2_uptake_extents, plot_title)
    Makie.save("april_10/aquifer.50mMCaO.png",fig)
    discharge_rates = fill(0.,n_cells)
    log_discharge_rates = fill(0.,n_cells)
    log_concs = fill(0.,n_solutes,n_cells)
    for i_cell in 1:n_cells
        discharge_rates[i_cell] = 1. * # m3 pore volume / m2 
            1. / ( flushing_age * i_cell ) * # m / yr 
            1.e3 / 365.         # mm / day 
        log_discharge_rates[i_cell] = log(discharge_rates[i_cell])
        for i_solute in 1:n_solutes
            log_concs[i_solute,i_cell] = log(solute_concentrations[i_solute,i_cell])
        end 
    end 
    Makie.plot(log_discharge_rates,log_concs[1,:])
  
    return solute_concentrations,
        mineral_saturation_indices,
        mineral_reaction_extents, 
        CO2_uptake_extents
end
function plot_aquifer(solute_concentrations,
                mineral_saturation_indices,
                mineral_reaction_extents, 
                CO2_uptake_extents, plot_title)
    fig = Figure(resolution=(1200, 1500))
    panel_for_title = fig[1,1] = GridLayout()
    Label(panel_for_title[1, 1], plot_title, tellwidth=false, textsize=30)
    panel_for_plots = fig[2,1] = GridLayout()
    left_panel_for_plots = panel_for_plots[1,1] = GridLayout()
    middle_slot_for_legends = panel_for_plots[1,2] = GridLayout()
    right_panel_for_plots = panel_for_plots[1,3] = GridLayout()
    panel_for_solutes_plot = left_panel_for_plots[1,1] = GridLayout()
    panel_for_omegas_plot = left_panel_for_plots[2,1] = GridLayout()
    panel_for_reaction_rate_plot = left_panel_for_plots[3,1] = GridLayout()
    panel_for_CO2_reaction_rate_plot = right_panel_for_plots[1,1] = GridLayout()
    panel_for_solute_comp = right_panel_for_plots[2,1] = GridLayout()
    panel_for_origins_barchart = right_panel_for_plots[3,1] = GridLayout()

    ax_solutes = Axis(panel_for_solutes_plot[1,1],
        title="Solute Concentrations")#,limits=(nothing,(0.,0.5)))
    ax_pH =  Axis(panel_for_solutes_plot[1,1],
        xlabelsvisible=false, ticklabels=false,
        ylabel="pH", yaxisposition=:right, rightspinecolor=:red,
        limits=(nothing,(6,10)))#,
        #limits=(nothing,(6.,9.5)))
    ax_reaction_rates = Axis(panel_for_reaction_rate_plot[1,1],
        title="Mineral Reactions")
    ax_omegas = Axis(panel_for_omegas_plot[1,1],
        title="Mineral Saturation Indices")
    ax_CO2_fluxes = Axis(panel_for_CO2_reaction_rate_plot[1,1],
        title="CO2 Uptake")#,limits=(nothing,(0.,0.3)))
    ax_solute_compare = Axis(panel_for_solute_comp[1,1],
        title="Solute Concentrations")
    ax_origins_barchart = Axis(panel_for_origins_barchart[1,1],
        title="Solute Source Minerals",
        xticks=(1:2:9, solute_names[1:5]),
        limits=(nothing,(-4.,4.)))

    x_locs = []
    for i_cell in 1:n_cells
        append!(x_locs,i_cell)
    end 
    for (index,i_solute) in enumerate([1,2,3,4,5,6,11,12])
        Makie.scatterlines!(ax_solutes,Float32.(x_locs),
            Float32.(solute_concentrations[i_solute,:]), 
            label=solute_names[i_solute],
            color=generic_colormap()[index], 
            markercolor=generic_colormap()[index],
            markersize=20)
    end
    middle_slot_for_legends[1, 1] = Legend(fig, ax_solutes)

    #Makie.axislegend(ax_solutes,position=:lt)

    Makie.scatterlines!(ax_pH,Float32.(x_locs),
        Float32.(solute_concentrations[pH_solute,:]), 
        label=solute_names[pH_solute],
        color=:black, 
        markercolor=:black,
        linestyle="-",
        linewidth=2,
        markersize=0)

    for (index, i_mineral) in enumerate([Anorthite_mineral,Albite_mineral,K_Feldspar_mineral,Calcite_mineral,Dolomite_mineral])
        Makie.scatterlines!(ax_omegas,Float32.(x_locs),
            Float32.(mineral_saturation_indices[i_mineral,:]), 
            label=mineral_names[i_mineral],
            color=generic_colormap()[index],
            #colormap=:rainbow,
            #colorrange=(1,n_minerals),
            #color=mineral_colormap()[i_mineral],
            markercolor=generic_colormap()[index],
            markersize=20)
    end
    #Makie.axislegend(ax_omegas,position=:lb)
    panel_for_two_legends_middle = middle_slot_for_legends[2, 1] = GridLayout()
    
    panel_for_two_legends_middle[1, 1] = Legend(fig, ax_omegas)


    #local_colormap = [colorant"red",colorant"blue",colorant"green"]
    index = 0
    for (index,i_mineral) in enumerate([Anorthite_mineral,Albite_mineral,K_Feldspar_mineral, Calcite_mineral,Dolomite_mineral])
        Makie.scatterlines!(ax_reaction_rates,Float32.(x_locs),
            Float32.(-mineral_reaction_extents[i_mineral,:] .* 1.e3), 
            label=mineral_names[i_mineral],
            #colormap=:rainbow,
            #colorrange=(1,n_minerals),
            color=generic_colormap()[index],
            #color=local_colormap[i_mineral],
            markercolor=generic_colormap()[index],
            markersize=20)
    end
    i_color = 0
    for i_mineral = secondary_silicate_minerals
        plot_me = false
        for i_cell in 1:n_cells
            if mineral_reaction_extents[i_mineral,i_cell] !== 0.
                plot_me = true
                i_color += 1
                if i_color > 12
                    i_color = 1
                end
            end
        end 
        if plot_me
            Makie.scatterlines!(ax_reaction_rates,Float32.(x_locs),
                Float32.(-mineral_reaction_extents[i_mineral,:]), 
                label=mineral_names[i_mineral],
                #colormap=:rainbow,
                #colorrange=(1,n_minerals),
                color=generic_colormap()[i_color],
                #color=mineral_colormap()[i_mineral],
                markercolor=generic_colormap()[i_color],
                linestyle="-",
                markersize=10)
        end
    end
    #Makie.axislegend(ax_reaction_rates,position=:lb)    
    panel_for_two_legends = middle_slot_for_legends[3, 1] = GridLayout()
    
    panel_for_two_legends[1,1] = Legend(fig, ax_reaction_rates)



    Makie.scatterlines!(ax_CO2_fluxes, 
        Float32.(x_locs),
        Float32.(CO2_uptake_extents), 
        #label="CO2 uptake",
        #colormap=:rainbow,
        #colorrange=(1,n_minerals),
        color=generic_colormap()[1],
        #color=mineral_colormap()[i_mineral],
        markercolor=generic_colormap()[1],
        #linestyle="-",
        markersize=20)
    #=Makie.scatterlines!(ax_CO2_fluxes, 
        Float32.(x_locs),
        Float32.(CO2_sinks_from_alkalinity), 
        label="CO2 sink from Alk release",
        #colormap=:rainbow,
        #colorrange=(1,n_minerals),
        color=generic_colormap()[2],
        #color=mineral_colormap()[i_mineral],
        markercolor=generic_colormap()[2],
        #linestyle="-",
        markersize=20)=#
    #Makie.axislegend(ax_CO2_fluxes,position=:lt)

    cell_picks = [1,10]
    for (index,i_solute) in enumerate([1,2,3,4,5,6])
        obs_concs = [typical_river_solute_concentrations[i_solute] * 1.e3,
            typical_river_solute_concentrations[i_solute] * 1.e3]
        model_concs = [solute_concentrations[i_solute, cell_picks[1]] * 1.e3,
            solute_concentrations[i_solute, cell_picks[2]] * 1.e3]
        solute_name = solute_names[i_solute]
        Makie.scatterlines!(ax_solute_compare,
            obs_concs,
            model_concs,
            color=generic_colormap()[index], 
            markercolor=generic_colormap()[index],
            #strokecolor=solute_colormap()[i_solute],
            markersize=20,
            #strokewidth=2,
            label = solute_name)
        solute_point = Point2f(obs_concs[1],
            model_concs[1])
        Makie.scatter!(ax_solute_compare, solute_point,
            color=:white, strokecolor=generic_colormap()[index],
            markersize=20,
            strokewidth=2)
    end
    Makie.scatterlines!(ax_solute_compare, [0.0, 500.0], [0.0, 500.0],
        markersize=0)
    #Makie.axislegend(ax_solute_compare,merge=true,position=:lt)

    for (index,i_solute) in enumerate([11,12])
        obs_concs = [typical_river_solute_concentrations[i_solute] * 1.e3,
            typical_river_solute_concentrations[i_solute] * 1.e3]
        model_concs = [solute_concentrations[i_solute, cell_picks[1]] * 1.e3,
            solute_concentrations[i_solute, cell_picks[2]] * 1.e3]
        scale_factor = 10
        solute_name = solute_names[i_solute] * "/" * string(scale_factor)
        Makie.scatterlines!(ax_solute_compare,
            obs_concs ./ scale_factor,
            model_concs ./ scale_factor,
            color=generic_colormap()[index+6], 
            markercolor=generic_colormap()[index+6],
            #strokecolor=solute_colormap()[i_solute],
            markersize=20,
            #strokewidth=2,
            label = solute_name)
        solute_point = Point2f(obs_concs[1] / scale_factor,
            model_concs[1] / scale_factor)
        Makie.scatter!(ax_solute_compare, solute_point,
            color=:white, strokecolor=generic_colormap()[index+6],
            markersize=20,
            strokewidth=2)
    end
    #Makie.scatterlines!(ax_CO2_compare, [0.0, 3000.0], [0.0, 3000.0],
    #    markersize=0)
    #Makie.axislegend(ax_solute_compare,merge=true,position=:lt)
    panel_for_two_legends_middle[2, 1] = Legend(fig, ax_solute_compare)

    solute_mineral_sources = fill(0.0, 5, n_mineral_groups)
    for i_solute in 1:5
        for i_mineral in 1:n_minerals
            i_mineral_group = mineral_groups[i_mineral]
            for i_cell in 1:n_cells
                solute_mineral_sources[i_solute, i_mineral_group] -=
                    mineral_reaction_extents[i_mineral, i_cell] *
                    mineral_component_stoic[i_mineral, i_solute]
            end
        end
    end
    #=for i_mineral_group in 2:n_mineral_groups
        for i_solute in 1:5
            solute_mineral_sources[i_solute, i_mineral_group] +=
                solute_mineral_sources[i_solute, i_mineral_group-1]
        end 
    end =#

    #mineral_group_colormap()
    cat_piece = [1,1,1,1,1]
    group_piece = [1, 3, 5, 7, 9]
    cat = []
    group = []
    height = []
    labels = []
    cum_solute_mineral_sources = fill(0.0, 5, n_mineral_groups)
    highest_positive = fill(0.,5)
    lowest_negative = fill(0.,5)
    #=for i_solute in 1:5
        if solute_mineral_sources[i_solute, i_mineral_group] > 0.
            highest_positive[i_solute] = 
                solute_mineral_sources[i_solute, i_mineral_group]
        else
            lowest_negative[i_solute] = 
                solute_mineral_sources[i_solute, i_mineral_group]
        end 
    end =#
    #        cum_solute_mineral_sources[i_solute,i_mineral_group] +=

    for i_mineral_group in 1:n_mineral_groups
        for i_solute in 1:5
            if solute_mineral_sources[i_solute, i_mineral_group] > 0.
                cum_solute_mineral_sources[i_solute, i_mineral_group] =
                    solute_mineral_sources[i_solute, i_mineral_group] +
                    highest_positive[i_solute]
                highest_positive[i_solute] = 
                    cum_solute_mineral_sources[i_solute, i_mineral_group]
            else
                cum_solute_mineral_sources[i_solute, i_mineral_group] =
                    solute_mineral_sources[i_solute, i_mineral_group] +
                    lowest_negative[i_solute]
                lowest_negative[i_solute] = 
                    cum_solute_mineral_sources[i_solute, i_mineral_group]
            end 
        end
    end
    for i_mineral_group = n_mineral_groups:-1:1
        #append!(cat,cat_piece)
        #append!(group,group_piece)
        #append!(height,solute_mineral_sources[1:5, i_mineral_group])
        #for i_solute in 1:5
        #    push!(labels,mineral_group_names[i_mineral_group])
        #end
        #cum_sources .+= solute_mineral_sources[1:5, i_mineral_group]

        barplot!(
            ax_origins_barchart,
            Float32.(group_piece),
            Float32.(cum_solute_mineral_sources[1:5,i_mineral_group] .* 1.e3),
            #color=mineral_group_colormap()[i_mineral_group],
            stack=Int.(cat_piece),
            color=mineral_group_colormap()[i_mineral_group],
            label=mineral_group_names[i_mineral_group])
        cat_piece .+= 1
    end 
        #println(i_mineral_group)
    #=barplot!(
        ax_origins_barchart,
        Float32.(group),
        Float32.(height),
        #color=mineral_group_colormap()[i_mineral_group],
        stack=Int.(cat),
        color=Int.(cat),
        label="fuckme"
    )=#

    panel_for_two_legends[2,1] = Legend(fig, ax_origins_barchart)
    #Makie.axislegend(ax_origins_barchart,merge=true,position=:lt)
    #solid_legend = Legend(panel_for_origin_bar_chart_and_legends[1, 2], 
    #    ax_origins_barchart, "Solid Phases")
    #    solid_legend.tellheight = true 
    return fig
end       
function aquifer_parameter_sensitivies()
   #base_driving_parameters = [aquifer_basalt_fraction,
    #    aquifer_calcite_fraction, aquifer_dolomite_fraction,
    #    porosity, layer_thickness, soil_CO2, soil_temperature, n_cells]
 
    sensitivity_run_names = ["CO2","Calcite","Basalt","AddMgO","AddCaO","AddCaOInert"]
    #=parameter_value_lists = [
        [1.e3,1.8e3,3.e3,6.e3,1.e4,1.8e4,3.e4,6e4,1.e5,1.8e5,3.e5,6.e5], # CO2
        [0.,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5], # Calcite
        [0.,0.01,0.018,0.03,0.06,0.1],
        [-1.e-3,-1.e-4,-1.e-5,-1.e-6,1.e-6,1.e-5,1.e-4,1.e-3]]=#
    parameter_value_lists = [
        [1.8e3,6.e3,1.e4,1.8e4,3.e4,6.e4,1.e5,1.8e5,3.e5], # CO2
        [0.,0.05,0.1,0.25,0.5], # Calcite
        [0.,0.05,0.1,0.15,0.2], # Basalt
        [-1.e-3,-3.e-4,-1.e-4,-3.e-5,-1.e-5,0.,1.e-5,3.e-5,1.e-4,3.e-4,1.e-3],
        [-1.e-3,-3.e-4,-1.e-4,-3.e-5,-1.e-5,0.,1.e-5,3.e-5,1.e-4,3.e-4,1.e-3],
        [-1.e-3,-3.e-4,-1.e-4,-3.e-5,-1.e-5,0.,1.e-5,3.e-5,1.e-4,3.e-4,1.e-3]]
    n_cells = 1
    CO2_uptake_rates = deepcopy(parameter_value_lists) .* 0.
    driver_values = deepcopy(parameter_value_lists) .* 0.
    for i_series in 1:size(sensitivity_run_names)[1]
        aquifer_basalt_fraction = 0.1
        aquifer_calcite_fraction = 0.1
        aquifer_dolomite_fraction = 0.01
        aquifer_dolomite_ratio = aquifer_dolomite_fraction / aquifer_calcite_fraction
        aquifer_inert_fraction = 0.
        porosity = 0.3; layer_thickness = 1.
        soil_CO2 = 3.e4; soil_temperature = 20. 
        reactant_added = 0.; reactant_name = ""
        # for animation
        image_number = 0
        file_list = readdir()
        for file = file_list 
            if file[1:3] == "img" 
                println("removing ", file)
                rm(file,force=true)
            end 
        end
        mp4_file = "aquifer_" * sensitivity_run_names[i_series] * ".mp4"
        #
        aquifer_out_filename = ""; plot_title = ""
        i_parameter_list = min(i_series,4)
        for (i_parm,parm_value) in enumerate(parameter_value_lists[i_parameter_list])
            if i_series == 1 
                soil_CO2 = parm_value
                aquifer_out_filename = "aquifer_CO2_" * string(parm_value) * ".out"
                plot_title = "Initial CO2 " * string(parm_value)
            elseif i_series == 2
                aquifer_dolomite_fraction = parm_value * 
                    aquifer_dolomite_ratio
                aquifer_calcite_fraction = parm_value
                aquifer_out_filename = "aquifer_calcite_" * string(parm_value) * ".out"
                plot_title = "Calcite Fraction " * string(parm_value)
            elseif i_series == 3 
                aquifer_basalt_fraction = parm_value
                aquifer_out_filename = "aquifer_basalt_" * string(parm_value) * ".out"
                plot_title = "Basalt Fraction " * string(parm_value)
            else  # 4 HCl/MgO, 5 HCl/CaO, 6 HCl/CaO inert
                reactant_name = "HCl(g)"
                reactant_added = parm_value
                #acid_added = parm_value
                if parm_value < 0
                    if i_series == 4
                        reactant_name = "Brucite"
                    else
                        reactant_name = "Portlandite"
                    end 
                    reactant_added = - parm_value / 2.
                end
                out_base_filename = "aquifer_react_"
                plot_title = "Adding " * string(reactant_added) * " eq of " * 
                    reactant_name
                if i_series == 6  
                    out_base_filename = "aquifer_inert_"
                    plot_title = plot_title * " Inert Bedrock"
                    aquifer_inert_fraction = 1.
                    aquifer_basalt_fraction, 
                        aquifer_calcite_fraction,
                        aquifer_dolomite_fraction = 0.,0.,0.
                end 
                aquifer_out_filename = out_base_filename * reactant_name * "_" *
                    string(reactant_added) * ".out"
                
                #append!(driving_parameters,reactant_added)
            end 
            println("running ", aquifer_out_filename)
            solute_concentrations,
                mineral_saturation_indices,
                mineral_reaction_extents, 
                CO2_uptake_extents =
            run_aquifer_phreeqc(
                aquifer_basalt_fraction, 
                aquifer_calcite_fraction,
                aquifer_dolomite_fraction, 
                aquifer_inert_fraction, 
                porosity, layer_thickness, 
                soil_CO2, soil_temperature, 
                reactant_name, reactant_added, 
                n_cells)
            cp("aquifer_phreeqc.out", aquifer_out_filename, force=true)

            #=solute_concentrations,
                mineral_saturation_indices,
                mineral_reaction_extents, 
                CO2_uptake_extents = 
            parse_aquifer_phreeqc_file(aquifer_out_filename, n_cells) =#

            #=fig = plot_aquifer(solute_concentrations,
                mineral_saturation_indices,
                mineral_reaction_extents, 
                CO2_uptake_extents,
                plot_title)
            image_number += 1
            image_file_name = "img." * lpad(image_number, 3, "0") * ".png"
            Makie.save(image_file_name,fig)
            println("saved ", image_file_name)=#

            CO2_uptake_rates[i_series][i_parm] = CO2_uptake_extents[n_cells]
            driver_values[i_series][i_parm] = parm_value
            #println(CO2_uptake_extents[10], " ", CO2_uptake_rates[i_series])
        end
        #=println("compiling ",mp4_file )
        rm(mp4_file, force=true)
        run(`ffmpeg -r 2 -f image2 -s 1920x1080 -i img.%03d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p $mp4_file`)
        =#
    end 
    #=i_mid = Int(( size(CO2_uptake_rates[4])[1] - 1 ) / 2 + 1 )
    CO2_uptake_rates[4][:] .-= CO2_uptake_rates[4][i_mid]
    CO2_uptake_rates[5][:] .-= CO2_uptake_rates[5][i_mid]
    CO2_uptake_rates[6][:] .-= CO2_uptake_rates[6][i_mid]=#
  
    f = Figure(resolution=(1800, 1500))
    panel_for_title = f[1,1] = GridLayout()
    plot_title = "CO2 uptake rates"
    Label(panel_for_title[1, 1], plot_title, tellwidth=false, textsize=30)
    panel_for_plots = f[2,1] = GridLayout()
    ax_CO2 = Axis(panel_for_plots[1,1], title="Soil CO2",
        ylabel = "uptake rate", xlabelsvisible = false, ticklabels=false )
    ax_calcite = Axis(panel_for_plots[1,2], title="Bedrock Calcite",
        ylabel = "uptake rate", xlabelsvisible = false, ticklabels=false )
    ax_basalt = Axis(panel_for_plots[2,1], title="Bedrock Basalt",
        ylabel = "uptake rate", xlabelsvisible = false, ticklabels=false )
    ax_react = Axis(panel_for_plots[2,2], title="Soil Treatment",
        ylabel = "uptake rate", xlabelsvisible = false, ticklabels=false )
    linkyaxes!(ax_CO2,ax_calcite,ax_basalt)
    Makie.scatterlines!(ax_CO2,Float32.(driver_values[1]),Float32.(CO2_uptake_rates[1]))
    Makie.scatterlines!(ax_calcite,Float32.(driver_values[2]),Float32.(CO2_uptake_rates[2]))
    Makie.scatterlines!(ax_basalt,Float32.(driver_values[3]),Float32.(CO2_uptake_rates[3]))
    Makie.scatterlines!(ax_react,Float32.(driver_values[4]),Float32.(CO2_uptake_rates[4]),
        markercolor=:red,label="Add MgO")
    Makie.scatterlines!(ax_react,Float32.(driver_values[5]),Float32.(CO2_uptake_rates[5]),
        markercolor=:green,label="Add CaO")
    Makie.scatterlines!(ax_react,Float32.(driver_values[6]),Float32.(CO2_uptake_rates[6]),
        markercolor=:blue, label = "Add CaO Inert Bedrock")
    Makie.axislegend(ax_react)
    f
    Makie.save("carbon_fluxes_lasaga.png",f)
    #return fig
end
function eroding_hilltop_model(runoff,layer_thickness,deep_exch,pCO2,t_scaling,file_name)
    # runoff = 1.; layer_thickness = 10.; deep_exch = 10.; pCO2 = 4000.; t_scaling = 1. 
    surface_temperature = 15.
    porosity = 0.2 
    
    layer_thicknesses = [layer_thickness, layer_thickness]
    layer_pore_volumes = layer_thicknesses .* porosity 
    runoff_volume = runoff # m / year 
    deep_aquifer_exchange_time = deep_exch
    fluid_flushing_timescales = fill(0.0, n_layers)
    runoff_volumes = fill(0.0, n_layers)
    #total_thickness = layer_thicknesses[regolith_layer] + layer_thicknesses[saprolite_layer]
    fluid_flushing_timescales[regolith_layer] = # years
        (layer_thicknesses[regolith_layer] * porosity) / # 
        runoff_volume
    fluid_flushing_timescales[saprolite_layer] = fluid_flushing_timescales[regolith_layer] + 
        deep_aquifer_exchange_time
    runoff_volumes[saprolite_layer] =
        layer_pore_volumes[saprolite_layer] /
        deep_aquifer_exchange_time
    runoff_volumes[regolith_layer] = runoff_volume - runoff_volumes[saprolite_layer]

    basement_mineral_fractions = 
        basement_mineral_fraction_endmembers[granite_endmember, :] .* 0.9 .+
        basement_mineral_fraction_endmembers[basalt_endmember, :] .* 0.1

    mineral_dissolution_solid_rate_constants =
        [0.00015, # Anorthite
            0.00001, # Albite
            0.0001, # K-Feldspar
            0.0, # Annite
            0.0, # Phlogopite
            0.0003, # Enstatite
            0.0, # Ferrosilite
            0.0000, # Wollastonite
            0.0000, # Forsterite
            0.0, # Fayalite
            0.0, # Quartz
            0.0] .* # Hematite
        t_scaling
    mineral_reaction_timescales = fill(NaN, n_minerals)
    for i_mineral in primary_minerals
        if mineral_dissolution_solid_rate_constants[i_mineral] != 0.0
            mineral_reaction_timescales[i_mineral] = 1.0 /
                mineral_dissolution_solid_rate_constants[i_mineral]
        end
    end

    erosion_rate = typical_river_solute_concentrations[Ca_2plus_solute] * # mol / m3 runoff
        runoff_volume * # m3 runoff / m2 year -> mol / m2 year
        mineral_molwts[Anorthite_mineral] * # g feldspar / m2 year
        1.0 / rho_continent_crust * # m3 feldspar / m2 year 
        1.0 / basement_mineral_fractions[Anorthite_mineral] # m3 solid / m2 year 

    layer_mineral_fractions = fill(0.,n_minerals,2)
    for i_mineral in primary_minerals 
        layer_mineral_fractions[i_mineral,saprolite_layer] =
            basement_mineral_fractions[i_mineral] 
        if mineral_reaction_timescales[i_mineral] == mineral_reaction_timescales[i_mineral]
            layer_mineral_fractions[i_mineral,saprolite_layer] *=
                erosion_rate / 
                (erosion_rate + 
                layer_thicknesses[saprolite_layer] / mineral_reaction_timescales[i_mineral])
        end
        layer_mineral_fractions[i_mineral, regolith_layer] =
            layer_mineral_fractions[i_mineral, saprolite_layer]
        if mineral_reaction_timescales[i_mineral] == mineral_reaction_timescales[i_mineral]
            layer_mineral_fractions[i_mineral, regolith_layer] *=
                erosion_rate /
                (erosion_rate +
                layer_thicknesses[regolith_layer] / mineral_reaction_timescales[i_mineral])
        end         
    end 

    mineral_reaction_extents = fill(0.,n_minerals,2)
    solute_concentrations = fill(0.,n_solutes,n_layers) # mol / m3 pw 
    mineral_saturation_indices = fill(0.,n_minerals,n_layers)
    mineral_reaction_extents = fill(0.,n_minerals+1,n_layers)
    for i_layer in 1:2
        for i_mineral in primary_minerals 
            mineral_reaction_extents[i_mineral,i_layer] = # want mol / L pw
                mineral_dissolution_solid_rate_constants[i_mineral] * #  
                layer_mineral_fractions[i_mineral,i_layer] * # m3 rxn / m3 solid yr 
                rho_continent_crust / 1.e3 * # g / L yr 
                1. / mineral_molwts[i_mineral] * # mol / L bulk yr 
                1. / porosity * # mol / L pw yr 
                fluid_flushing_timescales[i_layer] # mol / L pw 
        end 
        generate_phreeqc_input_file_from_minerals(
            surface_temperature, pCO2,
            - mineral_reaction_extents[:,i_layer])
        
        cmd = `$exec_command`
        run(cmd)
        solute_concentrations[:,i_layer], # mol / m3 pw 
            mineral_saturation_indices[:,i_layer],
            mineral_reaction_extents[:,i_layer] = # mol / m3 pw
            parse_phreeqc_output_file()
    end 

    f = open(file_name, "w")

    mean_solute_concentrations = solute_concentrations[:,regolith_layer] + 
        solute_concentrations[:, saprolite_layer] / 
        fluid_flushing_timescales[saprolite_layer]
        println(f,"      shallow solutes ")
    print_table(f, solute_names, solute_concentrations[:, 1],
        typical_river_solute_concentrations)
    println(f, "     deep ")
    print_table(f,solute_names, solute_concentrations[:,2],
        typical_river_solute_concentrations)
    println(f,"      combined ")
    river_concs = ( 
        solute_concentrations[:, 1] .* runoff_volume .+ 
        solute_concentrations[:, 2] * saprolite_runoff_volume ) ./
        runoff_volume
    print_table(f,solute_names, river_concs,
        typical_river_solute_concentrations)
    println(f)
    #=for i_layer in 1:2
        println(
            mineral_abundance_report(primary_minerals, layer_mineral_fractions[:, i_layer]) * " .... " *
            mineral_abundance_report(secondary_silicate_minerals, mineral_reaction_extents[:, i_layer]) * " .... " *
            base_cation_abundance_report(primary_minerals, layer_mineral_fractions[:, i_layer]))
    end =#
    println(f,"      primary minerals regolith")
    mineral_sorted_table(f,primary_minerals, layer_mineral_fractions[:,1])
    println(f,"             saprolite")
    mineral_sorted_table(f,primary_minerals, layer_mineral_fractions[:, 2])
    println(f,"               basement")
    mineral_sorted_table(f,primary_minerals, basement_mineral_fractions)
    println(f,"      secondary regolith")
    mineral_sorted_table(f,secondary_silicate_minerals, mineral_reaction_extents[:, 1])
    println(f,"            saprolite")
    mineral_sorted_table(f,secondary_silicate_minerals, mineral_reaction_extents[:, 2])
    mineral_reaction_rates = deepcopy(mineral_reaction_extents)
    mineral_deep_fraction = fill(0.,n_minerals)
    total_reaction = 0.
    println(f,"      deep rxn percent")
    for i_mineral in 1:n_minerals
        mineral_reaction_rates[i_mineral,:] = 
            mineral_reaction_extents[i_mineral,:] ./
            fluid_flushing_timescales[:]
        total_reaction =
            mineral_reaction_rates[i_mineral, 1] +
            mineral_reaction_rates[i_mineral, 2]
        if total_reaction !== 0.
            mineral_deep_fraction[i_mineral] = mineral_reaction_rates[i_mineral, 2] /
                total_reaction
            println(f,mineral_names[i_mineral], "  ", Int(floor(mineral_deep_fraction[i_mineral]*100)))
        end 
    end 
    println(f,"       shallow rxn mol/m3 yr")
    for i_mineral in 1:n_minerals
        if mineral_reaction_rates[i_mineral, 1] !== 0.0
            println(f, mineral_names[i_mineral], "  ",
                Int(floor(mineral_reaction_rates[i_mineral, 1] * 100.0)) / 100.0)
        end
    end
    println(f, "       deep ")
    for i_mineral in 1:n_minerals
        if mineral_reaction_rates[i_mineral, 2] !== 0.0
            println(f, mineral_names[i_mineral], "  ",
                Int(floor(mineral_reaction_rates[i_mineral, 2] * 100.0)) / 100.0)
        end
    end
    println(f,"       shallow rxn mol/m3 pw ")
    for i_mineral in 1:n_minerals
        if mineral_reaction_extents[i_mineral,1] !== 0.0
            println(f,mineral_names[i_mineral], "  ", 
                Int(floor(mineral_reaction_extents[i_mineral,1] * 100.)) / 100.) 
        end
    end
    println(f,"       deep ")
    for i_mineral in 1:n_minerals
        if mineral_reaction_extents[i_mineral, 2] !== 0.0
            println(f,mineral_names[i_mineral], "  ",
                Int(floor(mineral_reaction_extents[i_mineral, 2] * 100.0)) / 100.0)
        end
    end
    close(f)
end 

# phreeqc utilities
function phreeqc_box( mineral_fractions, 
    layer_thickness, 
    porosity, 
    flushing_rate, 
    flushing_age, 
    surface_temperature, soil_CO2,
    time_step, thread_number ) 
    #=
    i_layer = saprolite_layer; i_box = 1
    mineral_fractions = box_mineral_fractions[:,i_layer,i_box]
    layer_thickness = box_layer_thicknesses[i_layer,i_box]
    porosity = box_porosities[i_layer,i_box]
    flushing_rate = box_fluid_flushing_rates[i_layer,i_box]
    flushing_age = box_fluid_flushing_ages[i_layer,i_box]
    surface_temperature = box_temperatures[i_box]
    soil_CO2 = box_soil_CO2_values[i_box]
    thread_number = 1
    =#
    # returned variables
    mineral_reaction_constants_for_phreeqc = fill(0.,n_minerals)

    pore_volume = layer_thickness * porosity
    fluid_volume = flushing_rate * flushing_age 
    pore_saturation = fluid_volume / pore_volume

    #mineral_solid_reaction_rates = fill(0.0, n_minerals)
    #mineral_porewater_reaction_rates = fill(0.,n_minerals)
    #mineral_meters_for_residual_array = fill(0.,n_minerals)
    #n_new_total_phreeqc_runs, n_new_interpolated_phreeqc_runs = 0,0
 
    buffered_mineral_fraction = 1.e-2 # set to 0 to turn off 
    buffered_mineral_order = 2.5
    for i_mineral in phreeqc_dissolving_minerals  
        buffered_term = 
            min(max(1.e-6,mineral_fractions[i_mineral]), buffered_mineral_fraction) ^ 
            buffered_mineral_order
        linear_term = max(0.,mineral_fractions[i_mineral] - buffered_mineral_fraction)
        mineral_reaction_constants_for_phreeqc[i_mineral] = 
            #mineral_fractions[i_mineral] * # m3 mineral / m3 solid 
            (buffered_term + linear_term) * # standin for mineral_fractions
            lasaga_rates_pH5[i_mineral] * # mol mineral / m3 pw sec 
            flushing_age * # mol / m3 pw step 
            1.e-3  # mol / l pw sec
    end

    output_block = # , n_new_total_phreeqc_runs, n_new_interpolated_phreeqc_runs =
        run_or_interpolate_phreeqc(surface_temperature, soil_CO2, 
            mineral_reaction_constants_for_phreeqc, thread_number)
            # if interpolation is enabled, tread_number will be a dummy 
            # and threading will happen inside run_phreeqc_to_omega

    mineral_porewater_reaction_rates = fill(0.,n_minerals)
    mineral_meters_for_residual_array = fill(0.,n_minerals)
    #mineral_omegas == fill(0.,n_minerals)
    for i_mineral in 1:n_minerals 
        #mineral_omegas[i_mineral] = 
        #    10 ^ output_block.mineral_saturation_indices[i_mineral]
        mineral_porewater_reaction_rates[i_mineral] = # mol / m3 pw yr 
            mineral_porewater_rate_from_phreeqc_rate(
                output_block.mineral_reaction_extents[i_mineral], # mol / l pw 
                flushing_age)
        mineral_meters_for_residual_array[i_mineral] = 
            mineral_solid_rate_from_porewater_rate(
                mineral_porewater_reaction_rates[i_mineral],
                mineral_mol_per_m3[i_mineral], 
                porosity, 
                layer_thickness) * 
            time_step * pore_saturation
    end 
    box_output_block = phreeqc_box_struct(
        output_block.driving_parameters,
        pore_saturation,
        output_block.solute_concentrations, 
        output_block.mineral_saturation_indices,
        output_block.mineral_omegas, 
        output_block.mineral_reaction_extents, 
        mineral_porewater_reaction_rates,
        mineral_meters_for_residual_array,
        output_block.CO2_uptake_extent)
    return box_output_block # , n_new_total_phreeqc_runs, n_new_interpolated_phreeqc_runs
end 
function solution_pCO2(alkalinity,total_CO2,temperature)
    if total_CO2 == 0. 
        return 0.
    end 
    input_file_name = 
        "/Users/archer/Synched/papers/gridplates/mineral_solubilities/auto_CO2.in"
    rm(input_file_name, force=true)
    f = open(input_file_name, "w")
    println(f, "SOLUTION 1")
    println(f, "    pH 7.0")
    println(f, "    temp " * string(temperature))
    println(f, "REACTION")
    cao_reaction = alkalinity * 1.e-3 / 2. 
    cao_reaction = max(cao_reaction,0.)
    println(f,"   CaO " * string(cao_reaction))
    println(f,"   CO2 " * string(total_CO2 * 1.e-3))
    println(f,"end")
    close(f)
    run(`/Users/archer/Synched/papers/gridplates/mineral_solubilities/run_pCO2.x`)
    output_file_name = input_file_name[1:end-2] * "out"
    f = open(output_file_name)
    file_lines = readlines(f)
    iline = 1
    still_looking_for_header = true
    while still_looking_for_header
        file_line = file_lines[iline]
        if findlast("Beginning of batch-reaction calculations", file_line) == nothing
            iline += 1
        else
            still_looking_for_header = false
        end
    end
    still_looking_for_pH = true 
    pH = NaN
    while still_looking_for_pH
        file_line = file_lines[iline]
        if findlast("pH", file_line) == nothing
            iline += 1
        else
            words = split(file_line)
            pH = parse(Float64, words[3])
            still_looking_for_pH = false
        end
    end

    still_looking_for_header = true
    while still_looking_for_header
        file_line = file_lines[iline]
        if findlast("Saturation indices", file_line) == nothing
            iline += 1
        else
            still_looking_for_header = false
        end
    end
    still_looking_for_CO2 = true
    pCO2 = 0. 
    while still_looking_for_CO2
        #println(iline)
        file_line = file_lines[iline]
        if findlast("CO2(g)", file_line) == nothing
            iline += 1
        else
            words = split(file_line)
            pCO2 = 10^(parse(Float64, words[2])) * 1.e6
            #println("CO2 ", words[2])
            still_looking_for_CO2 = false 
        end
    end
    return pCO2, pH 
end 
function run_or_interpolate_phreeqc( 
    surface_temperature, soil_CO2, 
    mineral_reaction_constants_for_phreeqc,
    thread_number)

    #= 
    n_memories_of_phreeqc_filled = 0
    n_memories_filled = n_memories_of_phreeqc_filled
    phreeqc_memory_pointers = Dict() 
    memories_of_phreeqc = Array{phreeqc_output_struct,1}(undef,n_memories_of_phreeqc)
     =#
    #println("running ", surface_temperature, " ", soil_CO2, component_reaction_fluxes)
    #output_block = 0
    #n_total_phreeqc_runs, n_interpolated_phreeqc_runs = 0,0

    # dissolving phases 
    n_dissolving_phases = 8
    dissolving_phase_names = ["Anorthite_phase","Albite_phase","K_Feldspar_phase,",
        "Mica_phase","Pyroxene_phase","Olivine_phase","Calcite_phase","Dolomite_phase"]
    phase_reaction_rate_constants = fill(0., n_dissolving_phases)
    for i_phase in 1:3
        phase_reaction_rate_constants[i_phase] = 
            mineral_reaction_constants_for_phreeqc[ i_phase ]
    end 
    phase_reaction_rate_constants[4] = 
        mineral_reaction_constants_for_phreeqc[ Annite_mineral ] / mica_mineral_fractions[1]
    phase_reaction_rate_constants[5] = 
        mineral_reaction_constants_for_phreeqc[ Enstatite_mineral ] / pyroxene_mineral_fractions[1]
    phase_reaction_rate_constants[6] = 
        mineral_reaction_constants_for_phreeqc[ Forsterite_mineral ] / olivine_mineral_fractions[1]
    phase_reaction_rate_constants[7] = 
        mineral_reaction_constants_for_phreeqc[ Calcite_mineral ]
    phase_reaction_rate_constants[8] = 
        mineral_reaction_constants_for_phreeqc[ Dolomite_mineral ]

    exact_driving_parameters = Float32.([[surface_temperature, soil_CO2] ; 
        phase_reaction_rate_constants])
    #n_total_phreeqc_runs, n_interpolated_phreeqc_runs = 0,0
    if enable_interpolated_phreeqc == true # enable_interpolated_phreeqc 
        #global n_memories_of_phreeqc_filled
        #n_memories_filled = n_memories_of_phreeqc_filled
        
        #println("this works? ", n_memories_filled)

        solute_concentrations = fill(0.0, n_solutes)
        mineral_saturation_indices = fill(0.0, n_minerals)
        mineral_reaction_extents = fill(0.0, n_minerals) 
        mineral_omegas = fill(0.,n_minerals)
        CO2_uptake_extent = 0.

        index_sizes, index_range_lists, fraction_lists, interp_domains = 
            interpolation_values_from_driving_parameters(exact_driving_parameters)

        fill_fractions = fill(0.,
            index_sizes[1],
            index_sizes[2],
            index_sizes[3],
            index_sizes[4],
            index_sizes[5],
            index_sizes[6],
            index_sizes[7],
            index_sizes[8],
            index_sizes[9],
            index_sizes[10])
        
        #=index_range_lists[1][1]:index_range_lists[1][2],
            index_range_lists[2][1]:index_range_lists[2][2],
            index_range_lists[3][1]:index_range_lists[3][2],
            index_range_lists[4][1]:index_range_lists[4][2],
            index_range_lists[5][1]:index_range_lists[5][2],
            index_range_lists[6][1]:index_range_lists[6][2],
            index_range_lists[7][1]:index_range_lists[7][2],
            index_range_lists[8][1]:index_range_lists[8][2],
            index_range_lists[9][1]:index_range_lists[9][2] )=#

        for ii in 1:index_sizes[1]
            fill_fractions[ii, :, :, :, :, :, :, :, :, :] .= fraction_lists[1][ii]
        end
        for ii in 1:index_sizes[2]
            fill_fractions[:, ii, :, :, :, :, :, :, :, :] .*= fraction_lists[2][ii]
        end
        for ii in 1:index_sizes[3]
            fill_fractions[:, :, ii, :, :, :, :, :, :, :] .*= fraction_lists[3][ii]
        end
        for ii in 1:index_sizes[4]
            fill_fractions[:, :, :, ii, :, :, :, :, :, :] .*= fraction_lists[4][ii]
        end
        for ii in 1:index_sizes[5]
            fill_fractions[:, :, :, :, ii, :, :, :, :, :] .*= fraction_lists[5][ii]
        end
        for ii in 1:index_sizes[6]
            fill_fractions[:, :, :, :, :, ii, :, :, :, :] .*= fraction_lists[6][ii]
        end
        for ii in 1:index_sizes[7]
            fill_fractions[:, :, :, :, :, :, ii, :, :, :] .*= fraction_lists[7][ii]
        end
        for ii in 1:index_sizes[8]
            fill_fractions[:, :, :, :, :, :, :, ii, :, :] .*= fraction_lists[8][ii]
        end
        for ii in 1:index_sizes[9]
            fill_fractions[:, :, :, :, :, :, :, :, ii, :] .*= fraction_lists[9][ii]
        end
        for ii in 1:index_sizes[10]
            fill_fractions[:, :, :, :, :, :, :, :, :, ii] .*= fraction_lists[10][ii]
        end
       
        # iL1,iL2,iL3,iL4,iL5,iL6,iL7,iL8 = 1,1,1,1,1,1,1,1

       
        total_fill_fraction = 0.
        driving_parameters_to_run_list = [] # only needing new runs
        full_address_tag_list = []               # for all node points
        rerun_address_tag_list = []
        memory_pointer_list = []
        fill_fraction_list = []
        
        for iL1 in 1:index_sizes[1]
            for iL2 in 1:index_sizes[2]
                for iL3 in 1:index_sizes[3]
                    for iL4 in 1:index_sizes[4]
                        for iL5 in 1:index_sizes[5]
                            for iL6 in 1:index_sizes[6]
                                for iL7 in 1:index_sizes[7]
                                    for iL8 in 1:index_sizes[8]
                                        for iL9 in 1:index_sizes[9]
                                            for iL10 in 1:index_sizes[10]
                                                # iL1,iL2,iL3,iL4,iL5,iL6,iL7,iL8,iL9,iL10 = 1,1,1,1,1,1,1,1,1,1
                                                # iL1,iL2,iL3,iL4,iL5,iL6,iL7,iL8,iL9,iL10 = 2,2,2,2,2,2,2,2,2,2
                                                iG1 = index_range_lists[1][iL1] # T 
                                                iG2 = index_range_lists[2][iL2] # CO2 
                                                iG3 = index_range_lists[3][iL3] # Anorthite
                                                iG4 = index_range_lists[4][iL4] # Albite
                                                iG5 = index_range_lists[5][iL5] # K-spar
                                                iG6 = index_range_lists[6][iL6] # mica
                                                iG7 = index_range_lists[7][iL7] # pyroxene
                                                iG8 = index_range_lists[8][iL8] # olivine 
                                                iG9 = index_range_lists[9][iL9] # calcite 
                                                iG10 = index_range_lists[10][iL10] # calcite 

                                                fill_fraction = fill_fractions[iL1, iL2, iL3, iL4, iL5, iL6, iL7, iL8, iL9, iL10]
                                                total_fill_fraction += fill_fraction
                                                append!(fill_fraction_list,fill_fraction)

                                                                        #T,  CO3,Ca,Mg,  Na, K,  Si, Al, CO2
                                                address_tag = string(iG1) * "," * string(iG2) * "," * 
                                                    string(iG3) * "," * string(iG4) * "," * 
                                                    string(iG5) * "," * string(iG6) * "," * 
                                                    string(iG7) * "," * string(iG8) * "," *
                                                    string(iG9) * "," * string(iG10)
                                                append!(full_address_tag_list,[address_tag])
                                                #println(address_tag)
                                                try 
                                                    memory_pointer = phreeqc_memory_pointers[address_tag] 
                                                    #println("found it ",address_tag) #, " ", 
                                                        #memories_of_phreeqc[memory_pointer].solute_concentrations[K_plus_solute] )
                                                    #append!(memory_pointer_list,memory_pointer)
                                                catch 
                                                    #var_indices = [iG3, iG4, iG5, iG6, iG7, iG8]
                                                    #println("need new node ", address_tag)
                                                    append!(rerun_address_tag_list,[address_tag])
                                                #end
                                                    gridpoint_surface_temperature = phreeqc_temperatures[iG1]
                                                    gridpoint_soil_CO2 = phreeqc_pCO2s[iG2]
                                                    gridpoint_mineral_reaction_rate_constants = fill(0.0, n_phreeqc_dissolving_phases)
                                                    gridpoint_mineral_reaction_rate_constants[1] = 
                                                        phreeqc_mineral_reaction_constant_grid[iG3]
                                                    gridpoint_mineral_reaction_rate_constants[2] = 
                                                        phreeqc_mineral_reaction_constant_grid[iG4]
                                                    gridpoint_mineral_reaction_rate_constants[3] =  
                                                        phreeqc_mineral_reaction_constant_grid[iG5]
                                                    gridpoint_mineral_reaction_rate_constants[4] = 
                                                        phreeqc_mineral_reaction_constant_grid[iG6]
                                                    gridpoint_mineral_reaction_rate_constants[5] = 
                                                        phreeqc_mineral_reaction_constant_grid[iG7]
                                                    gridpoint_mineral_reaction_rate_constants[6] = 
                                                        phreeqc_mineral_reaction_constant_grid[iG8]
                                                    gridpoint_mineral_reaction_rate_constants[7] = 
                                                        phreeqc_mineral_reaction_constant_grid[iG9]
                                                    gridpoint_mineral_reaction_rate_constants[8] = 
                                                        phreeqc_mineral_reaction_constant_grid[iG10]

                                                    gridpoint_driving_parameters = [[gridpoint_surface_temperature, gridpoint_soil_CO2] ; 
                                                        gridpoint_mineral_reaction_rate_constants]
                                                    append!(driving_parameters_to_run_list,[gridpoint_driving_parameters])
                                                    #n_memories_filled += 1
                                                    #append!(memory_pointer_list,n_memories_filled)
                                                    
                                                    
                                                    #output_block = 
                                                    #    run_phreeqc(
                                                    #        surface_temperature, soil_CO2,
                                                    #        gridpoint_component_reaction_fluxes, 
                                                    #        1)
                                                    #println([iG1,iG2,iG3,iG4,iG5,iG6,iG7,iG8,iG9])
                                                    #println(driving_parameters)
                                                    #memory_pointer += 1 
                                                    #phreeqc_memories[address_tag] = memory_pointer
                                                    #n_memories_of_phreeqc_filled = n_memories_of_phreeqc_filled + 1
                                                    #n_memories_of_phreeqc_filled += 1
                                                    #memories_of_phreeqc[memory_pointer+1] =
                                                        #output_block
                                                    #phreeqc_memory_addresses[iG1,iG2,iG3,iG4,iG5,iG6,iG7,iG8,iG9] =
                                                    #    n_memories_of_phreeqc_filled
                                                    #println("filling ",memory_pointer+1,[iG1,iG2,iG3,iG4,iG5,iG6,iG7,iG8,iG9])
                                                #else
                                                #    println("using ",n_memories_of_phreeqc_filled,[iG1,iG2,iG3,iG4,iG5,iG6,iG7,iG8,iG9])
                                                end
                                            end 
                                        end
                                    end 
                                end
                            end
                        end
                    end
                end
            end
        end
   
        n_total_phreeqc_runs = size(full_address_tag_list)[1]
        n_threads = 1
        if enable_multithread_phreeqc
            n_threads = Threads.nthreads()
        end
        n_runs_to_do = size(driving_parameters_to_run_list)[1]
        n_interpolated_phreeqc_runs = n_total_phreeqc_runs - n_runs_to_do
        block_buffer = Array{phreeqc_output_struct,1}(undef,n_threads)

        n_threaded_blocks = Int(trunc(
            n_runs_to_do / n_threads ) ) + 1
        #if n_threaded_blocks > 1 && n_threads > 1
            #println("running ", n_threaded_blocks, " blocks of ", n_threads, 
            #    " filling to ", size(memories_of_phreeqc)[1] + n_runs_to_do)    
        #end 
        for i_threaded_block in 1:n_threaded_blocks
            if n_threaded_blocks > 1
                println("running ", i_threaded_block, " of ", n_threaded_blocks, 
                    " blocks filled ", size(memories_of_phreeqc)[1] + n_threads)
            end 
            n_threads_in_block = n_threads
            if i_threaded_block == n_threaded_blocks
                n_threads_in_block = n_runs_to_do - 
                    ( n_threaded_blocks - 1 ) * n_threads 
            end 
            if n_threads_in_block > 0
                #println("block ", i_threaded_block, " running ", n_threads_in_block)
                Threads.@threads for i_thread in 1:n_threads_in_block
                    i_list = ( i_threaded_block - 1 ) * n_threads + i_thread 
                    #println("block ", i_threaded_block, " running ", i_list)
                    gridpoint_driving_parameters = driving_parameters_to_run_list[i_list]
                    
                    #gridpoint_driving_parameters[9:10] = exact_driving_parameters[9:10]

                    block_buffer[i_thread] = 
                        run_phreeqc_direct(gridpoint_driving_parameters, i_thread)
                end 
                for i_thread in 1:n_threads_in_block
                    push!(memories_of_phreeqc, block_buffer[i_thread])
                    memory_pointer = size(memories_of_phreeqc)[1]
                    i_list = ( i_threaded_block - 1 ) * n_threads + i_thread 
                    address_tag = rerun_address_tag_list[i_list]
                    phreeqc_memory_pointers[address_tag] = memory_pointer
                end 
                #=
                for i_blow_up_soon in 1:100000000
                    try 
                        push!(memories_of_phreeqc, output_block)
                        phreeqc_memory_pointers[i_blow_up_soon] = i_blow_up_soon
                    catch 
                        println(length(phreeqc_memory_pointers))
                    end 
                end
                println(length(phreeqc_memory_pointers))
                println(sizeof(memories_of_phreeqc)/1e9)
                save_phreeqc_memories("testing")
                =#
            end
        end 

        #n_memories_of_phreeqc_filled = size(memories_of_phreeqc)[1]
        #println("filled slots ", n_memories_of_phreeqc_filled)

        for (i_list, address_tag) in enumerate(full_address_tag_list)
            memory_pointer = phreeqc_memory_pointers[address_tag]
            gridpoint_solute_concentrations = 
                #gridpoint_outblock.solute_concentrations 
                #phreeqc_memories[address_tag].solute_concentrations
                memories_of_phreeqc[memory_pointer].solute_concentrations
            gridpoint_mineral_saturation_indices = 
                memories_of_phreeqc[memory_pointer].mineral_saturation_indices
            gridpoint_mineral_omegas =  
                memories_of_phreeqc[memory_pointer].mineral_omegas
                #phreeqc_memories[address_tag].mineral_saturation_indices
            gridpoint_mineral_reaction_extents = 
                #gridpoint_outblock.mineral_reaction_extents 
                #phreeqc_memories[address_tag].mineral_reaction_extents
                memories_of_phreeqc[memory_pointer].mineral_reaction_extents
                #println("dolomite reaction, omega ", gridpoint_mineral_reaction_extents[Dolomite_mineral],
                #    " ", gridpoint_mineral_saturation_indices[Dolomite_mineral] )
            gridpoint_CO2_uptake_extent = 
                #gridpoint_outblock.summary 
                #phreeqc_memories[address_tag].summary
                memories_of_phreeqc[memory_pointer].CO2_uptake_extent
            solute_concentrations .+=
                gridpoint_solute_concentrations .*
                fill_fraction_list[i_list]
            mineral_saturation_indices .+=
                gridpoint_mineral_saturation_indices .*
                fill_fraction_list[i_list]
            mineral_omegas .+=
                gridpoint_mineral_omegas .*
                fill_fraction_list[i_list]
            mineral_reaction_extents .+=
                gridpoint_mineral_reaction_extents .*
                fill_fraction_list[i_list]
            CO2_uptake_extent += 
                gridpoint_CO2_uptake_extent *
                fill_fraction_list[i_list]
        end 
        
        output_block = phreeqc_output_struct(
            exact_driving_parameters,
            solute_concentrations, 
            mineral_saturation_indices,
            mineral_omegas,
            mineral_reaction_extents, 
            CO2_uptake_extent)

        if enable_check_interpolated_phreeqc 
            direct_output_block =
                run_phreeqc_direct(exact_driving_parameters, 1)
            misfits = block_results_misfit(output_block,direct_output_block)
            if misfits[0] > check_interpolated_phreeqc_tolerance
                #println("interplated, direct")
                compare_two_output_blocks(output_block,direct_output_block,
                    misfits, check_interpolated_phreeqc_tolerance)
                #println()
            else 
                println("interp check ok ", misfits)
                #println()
                #println("misfits = ", misfits)
            #println("misfits tot,c,r,o ", block_results_misfit(output_block,direct_output_block))
            end 
            #output_block = deepcopy(direct_output_block)
        end # enable_check_interpolated_phreeqc == false 
    else # enable_interpolated_phreeqc == false
        output_block =
            run_phreeqc_direct(exact_driving_parameters,thread_number)
        n_total_phreeqc_runs = 1
    end
    #println()

    return output_block
end 
function run_phreeqc_direct(driving_parameters, thread_number) 
    surface_temperature = driving_parameters[1]
    soil_CO2 = driving_parameters[2]
    mineral_reaction_constants_for_phreeqc = fill(0.,n_minerals)
    mineral_reaction_constants_for_phreeqc[Anorthite_mineral:K_Feldspar_mineral] = 
        driving_parameters[3:5]
    mineral_reaction_constants_for_phreeqc[Annite_mineral:Phlogopite_mineral] = 
        mica_mineral_fractions .* driving_parameters[6]
    mineral_reaction_constants_for_phreeqc[Enstatite_mineral:Wollastonite_mineral] = 
        pyroxene_mineral_fractions .* driving_parameters[7]
    mineral_reaction_constants_for_phreeqc[Forsterite_mineral:Fayalite_mineral] = 
        olivine_mineral_fractions .* driving_parameters[8]
    mineral_reaction_constants_for_phreeqc[Calcite_mineral:Dolomite_mineral] = 
        driving_parameters[9:10]
    #omega_unlimited_carbonate_dissolution_rates = 
    #    fill(0.,Calcite_mineral:Dolomite_mineral)
    #omega_unlimited_carbonate_dissolution_rates[Calcite_mineral:Dolomite_mineral] =
    #    driving_parameters[9:10]
    #carbonate_omegas = fill(0.,n_minerals)
    #carbonate_omegas[Calcite_mineral:Dolomite_mineral] = driving_parameters[end-1:end]
    output_block = run_phreeqc_direct( surface_temperature, soil_CO2,
        mineral_reaction_constants_for_phreeqc,
        thread_number )
    return output_block 
end
function run_phreeqc_direct( surface_temperature, soil_CO2,
    mineral_reaction_constants_for_phreeqc, 
    thread_number )
    output_block = run_phreeqc_direct( surface_temperature, soil_CO2,
        mineral_reaction_constants_for_phreeqc, 
        1, 0., 
        thread_number )
    return output_block
end 
function run_phreeqc_direct( surface_temperature, soil_CO2,
    mineral_reaction_constants_for_phreeqc, 
    n_cells, acid_added, thread_number )
    generate_phreeqc_input_file( # omegas and rxns 
        surface_temperature, soil_CO2,
        mineral_reaction_constants_for_phreeqc,
        n_cells, acid_added, thread_number)
    run_command = exec_command[1:end-2] * "_" * string(thread_number) * ".x"
    cmd = `$run_command`
    try 
        run(cmd)
    catch
        error("failed thread ", thread_number," ", surface_temperature, " ", soil_CO2,
            mineral_reaction_rate_constants) #,carbonate_omegas)
    end
    #output_block = # mol / m3 pw
    driving_parameters = Float32.([[surface_temperature, soil_CO2] ; mineral_reaction_constants_for_phreeqc])
    #try 
        #println(thread_number)
        #global solute_concentrations, 
        #    mineral_saturation_indices, 
        #    mineral_reaction_extents, 
        #    summary
    solute_concentrations, 
        mineral_saturation_indices, 
        mineral_reaction_extents, 
        CO2_uptake_extent = 
            parse_phreeqc_output_file(thread_number)
    mineral_omegas = fill(Float32(0.),n_minerals)
    for i_mineral in 1:n_minerals 
        mineral_omegas[i_mineral] = 10 ^ mineral_saturation_indices[i_mineral]
    end 
    
    #catch 
    #    error("failed ", driving_parameters,thread_number)
    #end 
    output_block = phreeqc_output_struct(        
        driving_parameters,
        solute_concentrations,
        mineral_saturation_indices,
        mineral_omegas,
        mineral_reaction_extents,
        CO2_uptake_extent)

    #push!(memories_of_phreeqc,output_block)
    #output_block.mineral_reaction_extents .+= mineral_reaction_extents_for_phreeqc
    #output_block.driving_parameters = 
    #    [[surface_temperature, soil_CO2] ; component_reaction_fluxes]
    return output_block
    #return solute_concentrations, 
    #    mineral_saturation_indices, 
    #    mineral_reaction_extents, 
    #    summary
end 
function generate_phreeqc_input_file( 
    surface_temperature, soil_CO2,
    mineral_reaction_rate_constants,
    n_cells, acid_added, thread_number)

    new_input_file_name = phreeqc_input_file_name
    if thread_number !== 0
        new_input_file_name = phreeqc_input_file_name[1:end-3] * "_" *
                              string(thread_number) * ".in"
    end
    rm(new_input_file_name, force=true)
    f = open(new_input_file_name, "w")
    println(f, "SOLUTION 0")
    println(f, "    pH 7.0")
    println(f, "    temp " * string(surface_temperature))
    println(f, "EQUILIBRIUM_PHASES 0")
    println(f, "    CO2(g) " * string(log10(soil_CO2 / 1.e6)))
    if acid_added !== 0.
        reactant_name = "HCl"; reactant_added = acid_added
        if acid_added < 0.
            reactant_name = "CaO"; reactant_added = - acid_added / 2.
        end 
        println(f,"REACTION 0")
        println(f,"   ", reactant_name, " 1. " )
        println(f,"  ", reactant_added)
    end
    println(f, "SAVE SOLUTION 0")
    println(f, "END")
    destination_cells = "1"
    if n_cells > 1
        destination_cells = destination_cells * "-" * string(n_cells)
    end 
    println(f, "SOLUTION_SPECIES") # magic incantation to soothe phreeqc
    println(f, "  H2O + 0.01e- = H2O-0.01")
    println(f, "  log_k   -9.0")
    println(f, "COPY SOLUTION 0 " * destination_cells)
    println(f, "END")
    println(f, "ADVECTION")
    println(f, "    -cells ", n_cells)
    println(f, "    -shifts ", n_cells)
    println(f, "    -time_step 1 year")
    println(f, "EQUILIBRIUM_PHASES " * destination_cells)
    for i_mineral in secondary_silicate_minerals # secondary_silicate_minerals
        println(f, "  " * mineral_names[i_mineral] * " 0. 0.")
    end
    println(f,"KINETICS " * destination_cells)
    for i_mineral in 1:n_minerals
        if mineral_reaction_rate_constants[i_mineral] > 0.
            println(f,mineral_names[i_mineral])
            println(f,"   -m 10.")
            println(f,"   -m0 10.")
            println(f,"   -parms ", 
                mineral_reaction_rate_constants[i_mineral]) #," ",
                #mineral_dissolution_pH_start_points[i_mineral]," ",
                #mineral_dissolution_pH_rate_orders[i_mineral] )
            println(f,"   -tol 1.e-8")
        end
    end 
    println(f, "RATES " * destination_cells)
    for i_mineral in carbonate_minerals
        println(f,mineral_names[i_mineral])
        println(f,"   -start")
        #println(f,"   10 pH_parm = PARM(2)")
        #println(f,"   20 pH = - LA(\"H+\")")
        println(f,"   30 aH = ACT(\"H+\")")
        println(f,"   40 aH_ref = 1.e-5")
        #println(f,"   50 IF pH > pH_parm THEN aH = aH_neut")
        println(f,"   60 pH5_rate = PARM(1)")
        println(f,"   70 rate_constant = pH5_rate * aH / aH_ref")
        println(f,"   80 moles = rate_constant * (1. - SR(\"", mineral_names[i_mineral], "\")) * TIME")
        if i_mineral !== Calcite_mineral
            println(f,"   90 IF ( SR(\"", mineral_names[i_mineral], "\") > 1. ) THEN moles = 0.")
        end 
        println(f,"   100 SAVE moles")
        println(f,"   -end")
    end 
    for i_mineral in primary_minerals 
        if mineral_reaction_rate_constants[i_mineral] > 0.
            println(f,mineral_names[i_mineral])
            println(f,"   -start")
            #println(f,"   10 pH_parm = PARM(2)")
            #println(f,"   20 pH = - LA(\"H+\")")
            println(f,"   30 aH = ACT(\"H+\")")
            println(f,"   40 aH_ref = 1.e-5")
            #println(f,"   50 IF pH > pH_parm THEN aH = aH_neut")
            println(f,"   60 pH5_rate = PARM(1)")
            println(f,"   70 rate_constant = pH5_rate * aH / aH_ref")
            if i_mineral == K_Feldspar_mineral
                println(f,"   75 rate_constant = rate_constant * (1. - SR(\"K-Feldspar\")) ")
            end 
            println(f,"   80 moles = rate_constant * TIME")
            println(f,"   90 IF ( SR(\"", mineral_names[i_mineral], "\") > 1. ) THEN moles = 0.")
            println(f,"   100 SAVE moles")
            println(f,"   -end")
        end 
    end
    close(f)
end
function parse_phreeqc_output_file(thread_number)
    n_cells = 1
    aquifer_out_filename = phreeqc_output_file_name
    if thread_number !== 0
        aquifer_out_filename =
            phreeqc_output_file_name[1:end-4] * "_" *
            string(thread_number) * ".out"
    end
    solute_concentrations, 
        mineral_saturation_indices, 
        mineral_reaction_extents, 
        CO2_uptake_extents = 
    parse_aquifer_phreeqc_file(aquifer_out_filename, n_cells)
    return solute_concentrations[:,1], 
        mineral_saturation_indices[:,1], 
        mineral_reaction_extents[:,1], 
        CO2_uptake_extents[1]
end 
function parse_aquifer_phreeqc_file(aquifer_out_filename, n_cells)
    solute_concentrations = fill(Float32(0.0),n_solutes,n_cells)
    mineral_saturation_indices = fill(Float32(0.0), n_minerals,n_cells )
    mineral_reaction_extents = fill(Float32(0.0), n_minerals,n_cells) 
    total_mineral_reaction_extents = fill(Float32(0.0), n_minerals)
    CO2_uptake_extents = fill(Float32(0.0), n_cells)

    f = open(aquifer_out_filename)
    file_lines = readlines(f)

    iline = 1
    iline = next_line_number(file_lines, iline, 
        "Beginning of batch-reaction calculations")
    initial_total_CO2, iline = 
        extract_value_from_next_line_number(file_lines, iline, 
            "Total CO2 (mol/kg)",5)    
    initial_total_CO2 *= 1.e3
    n_shifts = n_cells # - 1   
    for i_cell in 1:n_cells 
        #println(i_cell)
        iline = 1
        cell_name = "Cell " * string(i_cell) * "."
        for i_blob_of_output in 1:n_shifts
            iline = next_line_number(file_lines, iline, cell_name) + 1
        end
        one_cell_mineral_reaction_extents,iline = 
        extract_values_from_lines(file_lines,iline,
            "-Phase assemblage-",
            mineral_names,1,2)
        one_cell_secondary_mineral_reaction_extents,iline = 
        extract_values_from_lines(file_lines,iline,
            "-Solution composition-",
            mineral_names,1,7)
        mineral_reaction_extents[:,i_cell] += 
            one_cell_mineral_reaction_extents + 
            one_cell_secondary_mineral_reaction_extents
        total_mineral_reaction_extents[:] +=
            mineral_reaction_extents[:,i_cell]
        solute_concentrations[:,i_cell], iline = 
            extract_values_from_lines(file_lines,iline,
                "-Description of solution-",
                solute_names,1,2)
        solute_concentrations[pH_solute,i_cell],iline = 
            extract_value_from_next_line_number(file_lines, iline, 
                "pH",3)
        solute_concentrations[Alk_solute,i_cell],iline = 
            extract_value_from_next_line_number(file_lines, iline, 
                "Total alkalinity",5)
        final_total_CO2, iline = 
            extract_value_from_next_line_number(file_lines, iline, 
            "Total CO2 (mol/kg)",5)
        final_total_CO2 *= 1.e3
        pre_sat_indices_iline = next_line_number(file_lines, iline, 
            "---Saturation indices---")
        mineral_saturation_indices[:,i_cell], iline = 
            extract_values_from_lines(file_lines,pre_sat_indices_iline,
                "**For a gas",
                mineral_names,1,2) 
        CO2_saturation, iline = 
            extract_value_from_next_line_number(file_lines, pre_sat_indices_iline, 
            "CO2(g)",2)
        solute_concentrations[1:Tot_CO2_solute,i_cell] .*= # mol / L
            1.e3 # mol / m3
        solute_concentrations[pCO2_solute,i_cell] = 10^(CO2_saturation) * 1.e6
        #total_CO2_change = 
            #( solute_concentrations[Tot_CO2_solute, i_cell] -
            #initial_total_CO2 )
        CO2_uptake_extents[i_cell] = 
            solute_concentrations[Alk_solute, i_cell] / 2. - 
            ( solute_concentrations[Tot_CO2_solute, i_cell] -
            initial_total_CO2 ) -
            total_mineral_reaction_extents[Calcite_mineral] -
            2. * total_mineral_reaction_extents[Dolomite_mineral]
    end 
    close(f)
    return solute_concentrations,
        mineral_saturation_indices,
        mineral_reaction_extents, 
        CO2_uptake_extents
end 

# phreeqc interpolation 
function interpolation_values_from_driving_parameters( 
    driving_parameters) 

    interp_domains = [
        [phreeqc_temperatures, driving_parameters[1]], # 1
        [phreeqc_pCO2s, driving_parameters[2]]]
    for i_dissolving_mineral in 3:size(driving_parameters)[1]
        push!(interp_domains, 
            [ phreeqc_mineral_reaction_constant_grid, 
                driving_parameters[i_dissolving_mineral]])
    end
        
        #=,                   # 2
        [phreeqc_mineral_reaction_constant_grid,   # 3, CaO 
            driving_parameters[3]],
        [phreeqc_mineral_reaction_constant_grid,   # 3, CaO 
            driving_parameters[3]],
        [phreeqc_Na2O_reaction_extents,   # 5, Na2O 
            component_reaction_fluxes[Na2O_component] /
            component_reaction_fluxes[SiO2_component] /
            phreeqc_component_silica_ratio_ranges[Na2O_component]],
        [phreeqc_K2O_reaction_extents,   # 6, K2O
            component_reaction_fluxes[K2O_component] /
            component_reaction_fluxes[SiO2_component] /
            phreeqc_component_silica_ratio_ranges[K2O_component]],
        [phreeqc_SiO2_reaction_extents,              # 7, SiO2  
            component_reaction_fluxes[SiO2_component]],
        [phreeqc_Al2O3_reaction_extents,   # 8, Al2O3 
            component_reaction_fluxes[Al2O3_component] /
            component_reaction_fluxes[SiO2_component] /
            phreeqc_component_silica_ratio_ranges[Al2O3_component]],
        [phreeqc_Calcite_max_reaction_rates, 
            - omega_unlimited_carbonate_dissolution_rates[Calcite_mineral]],
        [phreeqc_Dolomite_max_reaction_rates, 
            - omega_unlimited_carbonate_dissolution_rates[Dolomite_mineral]]]
    =#

        #=[phreeqc_CO2_alk_reaction_ratios,  # 9, CO2 
            component_reaction_fluxes[CO2_component] / 
            ( 2 * component_reaction_fluxes[CaO_component] + 
              2 * component_reaction_fluxes[MgO_component] + 
              component_reaction_fluxes[Na2O_component] +
              component_reaction_fluxes[K2O_component] + 
              3. * component_reaction_fluxes[Al2O3_component])]
            ] =#

    index_sizes = []
    index_range_lists = []
    fraction_lists = []
    #parameter_list = []
    for i_domain in 1:size(interp_domains)[1]
        #println(i_domain)
        interp_domain = interp_domains[i_domain]
        index_size, indices, fractions = find_fill_fractions(interp_domain[1], interp_domain[2])
        push!(index_sizes, index_size)
        push!(index_range_lists, indices)
        push!(fraction_lists, fractions)
        #push!(parameter_list, interp_domain[2])
    end

    return index_sizes, index_range_lists, fraction_lists, interp_domains
end 
function reset_phreeqc_memories()
    global memories_of_phreeqc = Array{phreeqc_output_struct,1}(undef,0)
    global phreeqc_memory_pointers = Dict()
end
phreeqc_memory_file_blocksize = 100000
function save_phreeqc_memories(id_string)
    file_name = "phreeqc_memory_pointers." * id_string * ".bson"
    files_in_directory = readdir()
    if file_name in files_in_directory
        backup_file_name = file_name * ".bak"
        mv(file_name, backup_file_name, force=true)
        println("archiving ", backup_file_name)
    end 
    BSON.@save file_name phreeqc_memory_pointers
    println("saved " * file_name)
    n_phreeqc_memory_files = Int(floor(size(memories_of_phreeqc)[1] / phreeqc_memory_file_blocksize )) + 1
    clobbered_the_last_file_already = false 
    for i_file in n_phreeqc_memory_files:-1:1
        file_name = "memories." * id_string * "." * string(i_file) * ".bson"
        saving_this_file = false 
        if file_name in files_in_directory 
            if clobbered_the_last_file_already == false 
                backup_file_name = file_name * ".bak"
                mv(file_name, backup_file_name,force=true)
                println("archiving ", file_name)
                clobbered_the_last_file_already = true 
                saving_this_file = true 
            else
                println("skipping ", file_name)
            end 
        else 
            saving_this_file = true 
        end 
        if saving_this_file
            begin_pointer = Int((i_file - 1) * phreeqc_memory_file_blocksize + 1)
            end_pointer = begin_pointer + phreeqc_memory_file_blocksize - 1
            if end_pointer > size(memories_of_phreeqc)[1]
                end_pointer = size(memories_of_phreeqc)[1]
            end 
            memory_block = memories_of_phreeqc[begin_pointer:end_pointer]
            println("saving ", file_name)
            BSON.@save file_name memory_block
        end 
    end 
end
function read_phreeqc_memories(id_string)
    file_name = "phreeqc_memory_pointers." * id_string * ".bson"
    BSON.@load file_name phreeqc_memory_pointers
    println("loaded " * file_name)
    n_phreeqc_memory_files = Int(floor(length(phreeqc_memory_pointers) / phreeqc_memory_file_blocksize )) +1
    memories_of_phreeqc = Array{phreeqc_output_struct,1}(undef,0)
    for i_file in 1:n_phreeqc_memory_files
        file_name = "memories." * id_string * "." * string(i_file) * ".bson"
        begin_pointer = Int((i_file - 1) * phreeqc_memory_file_blocksize + 1)
        end_pointer = begin_pointer + phreeqc_memory_file_blocksize - 1
        if end_pointer > length(phreeqc_memory_pointers)
            end_pointer = length(phreeqc_memory_pointers)
        end 
        #memory_block = memories_of_phreeqc[begin_pointer:end_pointer]
        BSON.@load file_name memory_block
        append!(memories_of_phreeqc, memory_block)
        println("read ", file_name)
    end 
    return phreeqc_memory_pointers, memories_of_phreeqc
end
function block_parameter_misfit(driving_parameters,driving_parameters2)
    total_misfit = 0.
    n_driving_parameters = size(driving_parameters)[1]
    misfit_scales = fill(1.,n_driving_parameters)
    misfit_scales[1] = 10. 
    for i_parm in 1:size(driving_parameters)[1] 
        misfit = 
            ( driving_parameters[i_parm] - 
            driving_parameters2[i_parm] )
        mean =  
            ( driving_parameters[i_parm] + 
            driving_parameters2[i_parm] )
        if abs(mean) > 0.
            total_misfit += (misfit / mean)^2 * misfit_scales[i_parm]
        end 
    end 
    return total_misfit 
end 
function block_results_misfit(block1,block2)
    total_misfits = fill(0.,0:3) # solutes, rxns, omegas 
    for i_solute in 1:n_solutes
        val1 = block1.solute_concentrations[i_solute]
        val2 = block2.solute_concentrations[i_solute]
        misfit = val1 - val2  
        mean = val1 + val2 
        abs_min = min(abs(val1),abs(val2))
        if abs(mean) > 0. && abs_min > 0.
            total_misfits[1] += (misfit/mean)^2 
        end 
    end 
    for i_mineral in secondary_minerals 
        val1 = block1.mineral_reaction_extents[i_mineral]
        val2 = block2.mineral_reaction_extents[i_mineral]
        misfit = val1 - val2  
        mean = val1 + val2 
        abs_min = min(abs(val1),abs(val2))
        if abs(mean) > 0. && abs_min > 0.
            total_misfits[2] += (misfit/mean)^2 
        end 
    end 
    for i_mineral in 1:n_minerals
        val1 = block1.mineral_saturation_indices[i_mineral]
        val2 = block2.mineral_saturation_indices[i_mineral]
        misfit = val1 - val2  
        mean = val1 + val2 
        abs_min = min(abs(val1),abs(val2))
        if abs(mean) > 0. && abs_min > 0.
            total_misfits[3] += (misfit/mean)^2 
        end 
        #println(i_mineral,"  ",misfit," ",mean," ",total_misfits[3])
    end 
    for i_fit in 1:3
        total_misfits[0] += total_misfits[i_fit]
    end 
    return total_misfits 
end 
function compare_interpolated_with_direct(output_block,driving_parameters)
    direct_output_block =
        run_phreeqc(driving_parameters, 
            1)
    results_misfits = block_results_misfit(output_block,direct_output_block)
    #fraction_RMS = 0.
    #for i_frac in 1:size(fraction_vector)[1]
    #    fraction_RMS += fraction_vector[i_frac]^2
    #end 
    #println("comparing ", direct_output_block, results_misfits, fraction_vector, fraction_RMS)
    compare_two_output_blocks(output_block, direct_output_block)
    return results_misfits
end 
function find_interpolation_axis_sensitivities(
    surface_temperature, soil_CO2, component_reaction_fluxes)

    #= 
    surface_temperature = 15.5
    soil_CO2 = 1000. 
    component_reaction_fluxes = [ # mol / L pw 
        0.0021441526406411424, 
        0.0005074868084956037, 
        0.0011419680829870213, 
        0.0002943196018873364, 
        0.013104343729840947, 
        0.0034361418817879484, 
        0.01]
    =# 

    index_sizes, index_range_lists, fraction_lists, interp_domains =
        interpolation_values_from_component_fluxes(
            surface_temperature, soil_CO2, component_reaction_fluxes, 
            omega_unlimited_carbonate_dissolution_rates)

    #=
    for i_size in 1:size(index_sizes)[1]
        index_sizes[i_size] = 1
        index_range_lists[i_size] = [index_range_lists[i_size][1]]
        fraction_lists[i_size] = [1.0]
    end
    component_reaction_fluxes = component_fluxes_from_interpolation_values( 
        index_sizes, index_range_lists, fraction_lists )
    surface_temperature = interp_domains[1][1][index_range_lists[1][1]]
    soil_CO2 = interp_domains[2][1][index_range_lists[2][1]]
    =# 


    parameter_list = [[surface_temperature, soil_CO2] ; component_reaction_fluxes; 
        omega_unlimited_carbonate_dissolution_rates[Calcite_mineral:Dolomite_mineral]] 
    output_block_ideal = 
    run_phreeqc_direct(parameter_list[1], 
        parameter_list[2],
        component_reaction_fluxes,
        omega_unlimited_carbonate_dissolution_rates, #enable_default_solid_solutions, 
        1)
    #i_solute = pCO2_solute ; thread_number = 1
    #println("base ", output_block_ideal.solute_concentrations[i_solute])
    #
    test_var = output_block_ideal.mineral_saturation_indices[Dolomite_mineral]
    println("base ", test_var)

    #= 
    n_memories_of_phreeqc_filled = 0
    n_memories_filled = n_memories_of_phreeqc_filled
    phreeqc_memory_pointers = Dict() 
     =#
    output_block_interp = interpolate_sparse_phreeqc_lattice( # with carbonate omega 
        surface_temperature, soil_CO2, component_reaction_fluxes,
        omega_unlimited_carbonate_dissolution_rates,thread_number) # runs comparison to stdout

    for i_variable in 1:10
        axis_mean = 0. 
        for i_run in 1:index_sizes[i_variable]
            altered_parameter_list = deepcopy(parameter_list)
            altered_parameter_list[i_variable] = # one side of interp range or other 
                interp_domains[i_variable][1][index_range_lists[i_variable][i_run]]
            #if i_variable in [3,4,5,6,8]
            #    altered_parameter_list[i_variable] *= 
            #        parameter_list[SiO2_component+2] *
            #        phreeqc_component_silica_ratio_ranges[i_variable-2]
            #end
                
                #altered_component_reaction_fluxes = component_fluxes_from_interpolation_values(
            #    index_sizes, index_range_lists, fraction_lists)
            output_block = 
                run_phreeqc_direct(altered_parameter_list, # enable_default_solid_solutions,
                    1)
            test_var = output_block.mineral_saturation_indices[Dolomite_mineral]
  
            axis_mean += test_var * 
                fraction_lists[i_variable][i_run]
            println(parameter_names[i_variable]," ", 
                fraction_lists[i_variable][i_run]," ",
                altered_parameter_list[i_variable], " ", 
                parameter_list[i_variable], " ", 
                test_var, " ", # solute_concentrations[i_solute], " ", 
                index_range_lists[i_variable][i_run])
        end 
        println("mean ", axis_mean)
    end 
 
end

# functions to pass to the sensitivity plotter 
function beryllium_isotopes(results_block)
    n_boxes = size(results_block.mineral_fraction_timeseries)[3]
    result = results_block.river_solute_timeseries[:, Be10_solute, end] ./
        results_block.river_solute_timeseries[:, Be9_solute, end]
    plot_name = "Dissolved Beryllium Isotope Ratio"
    variable_names = ["Andean", "Amazon"]
    return result, plot_name, variable_names
end
function lithium_isotopes(results_block)
    n_boxes = size(results_block.mineral_fraction_timeseries)[3]
    result = results_block.river_solute_timeseries[:, Li7_solute, end] 
    plot_name = "Dissolved Lithium Isotope Ratio"
    variable_names = ["Andean","Amazon"]
    return result, plot_name, variable_names
end
function Li_dissolved_fraction(results_block)
    results = []
    variable_names = []
    for i_river in [1, 2]
            i_column = 10 + i_river - 1 
            result =
                results_block.tracer_summary_timeseries[i_column, end]
            push!(results, result)
            variable_name = "Andean "
            if i_river == 2
                variable_name = "Amazon "
            end
            variable_name = variable_name 
            push!(variable_names, variable_name)
    end
    plot_name = "Li Solid Phase Runoff Fraction"
    return results, plot_name, variable_names
end
function Li_budgets(results_block)
    results = []
    variable_base_names = ["Dissolution","Runoff"]
    variable_names = []
    for i_river in [1,2]
        for i_flux in [1,2]
            i_column = 6 + i_river - 1 + (i_flux-1) * 2
            result = 
                results_block.tracer_summary_timeseries[i_column, end]
                push!(results,result)
            variable_name = "Andean "
            if i_river == 2
                variable_name = "Amazon "
            end 
            variable_name = variable_name * variable_base_names[i_flux]
            push!(variable_names, variable_name)            
        end 
    end
    plot_name = "Li flux terms"
    return results, plot_name, variable_names
end
function floodplain_clay_richness(results_block)
    n_boxes = size(results_block.mineral_fraction_timeseries)[3]
    result = results_block.mineral_fraction_timeseries[Illite_mineral, 1,n_boxes, end] / 
        (results_block.mineral_fraction_timeseries[Illite_mineral, 1, n_boxes, end] + 
        results_block.mineral_fraction_timeseries[Kaolinite_mineral, 1, n_boxes, end])
    plot_name = "Illite / Illite + Kaolinite"
    variable_names = [""]
    return result, plot_name, variable_names
end
function CaO_reaction(results_block)
    result = results_block.component_summary_timeseries[1, 2, end]
    plot_name = "CaO reaction"
    variable_names = [""]
    return result, plot_name, variable_names
end
function base_cation_summary(results_block)
    results = []
    i_component = 1
    for i_flux in [1,2,3,5]
        result = results_block.component_summary_timeseries[i_component, i_flux, end]
        push!(results,result)
    end
    plot_name = "CaO flux terms"
    variable_names = ["Erosion","Reaction","Runoff","Growth"]
    return results, plot_name, variable_names
end
function feldspar_erosion(results_block)
    result = results_block.mineral_summary_timeseries[1,1, end]
    variable_name = "Feldspar erosion"
    return result, variable_name, [""]
end
function feldspar_reaction(results_block)
    result = results_block.mineral_summary_timeseries[1, 2, end]
    plot_name = "Feldspar reaction"
    variable_names = [""]
    return result, plot_name, variable_names
end
function extract_CO2_uptake_rate(results_block,i_step)
    results = [results_block.tracer_summary_timeseries[1,i_step], 
        results_block.tracer_summary_timeseries[2,i_step],
        results_block.tracer_summary_timeseries[1,i_step] +
        results_block.tracer_summary_timeseries[2,i_step]]
    step_time = results_block.timepoints[i_step]
    plot_name = "Total C Uptake After " * string(step_time) * " Years"
    variable_names = ["Regolith","Saprolite","Combined"]
    return results,plot_name, variable_names
end 

# colormaps
function mineral_group_colormap()
     cmap = [colorant"red", # Plag
        colorant"yellow",   # K-feldspar
        colorant"pink",    # Micas
        colorant"aqua",    # Pyroxene 
        colorant"blue",    # Olivine
        colorant"brown",    # Quartz and Fe
        colorant"orange",   # Smectites
        colorant"pink",  # Micas 2nd
        colorant"violet",   # Kaolinite+Gibbsite
        colorant"lightgrey" ] # Carbonates
    return cmap
end
function get_stripe_coords(cumulative_fractions, next_cumulative_fractions, x_locs)
    n_boxes = size(cumulative_fractions)[1]
    point_list = rand(Point2f, 2 * n_boxes + 1) .* 0.0 # fill(0,n_steps) #{Point{2, Float64}}#poly_xs = []; poly_ys = []
    i_place_in_array = 0
    for i_box in 1:n_boxes
        #x_loc = (i_box-1) * box_width 
        new_point = Point2f(x_locs[i_box], cumulative_fractions[i_box])
        i_place_in_array += 1
        point_list[i_place_in_array] += new_point#append!(poly_xs,i_box); append!(poly_ys,cumulative_fractions[i_box])
    end
    for i_box in n_boxes:-1:1
        #x_loc = (i_box - 1) * box_width
        new_point = Point2f(x_locs[i_box], next_cumulative_fractions[i_box])
        i_place_in_array += 1
        point_list[i_place_in_array] += new_point
    end
    point_list[end] = point_list[1]
    return point_list
end
function generic_colormap()
        cmap = [colorant"red", # Ca
        colorant"green",   # Mg
        colorant"violet",    # Na
        colorant"grey",    # K 
        colorant"blue",    # Si 
        colorant"yellow",  # Al
        colorant"black",   # Li
        colorant"orange",  # C 
        colorant"brown",   # pCO2
        colorant"tan",
        colorant"magenta",
        colorant"lightblue"]   # pH
    return cmap 
end 
function solute_colormap()
    cmap = [colorant"red", # Ca
        colorant"green",   # Mg
        colorant"violet",    # Na
        colorant"grey",    # K 
        colorant"blue",    # Si 
        colorant"yellow",  # Al
        colorant"black",   # Li
        colorant"black",   # Li7
        colorant"black",   # Alk
        colorant"orange",  # C 
        colorant"black",   # pCO2
        colorant"black"]   # pH
    #=for i_solute in 1:n_solutes
        if i_solute == Ca_2plus_solute
            push!(cmap, colorant"red")
        elseif i_solute == Mg_2plus_solute
            push!(cmap, colorant"green")
        elseif i_solute == Na_plus_solute 
            push!(cmap, colorant"grey")
        elseif i_solute == K_plus_solute
            push!(cmap, colorant"grey")
        elseif i_solute == Si_4plus_solute
            push!(cmap, colorant"blue")
        elseif i_solute == Li_solute
            push!(cmap, colorant"yellow")
        end
    end =#
    return cmap 
end 
function mineral_colormap()
    cmap = fill(colorant"grey",n_minerals)
    for i_mineral in 1:n_minerals
        if i_mineral == Anorthite_mineral
            cmap[i_mineral] = colorant"red"
        elseif i_mineral == Calci
        elseif i_mineral == Albite_mineral
            push!(cmap, colorant"white")
        elseif i_mineral == K_Feldspar_mineral
            push!(cmap, colorant"blue")
        elseif mineral_component_stoic[i_mineral,MgO_component] > 0
            push!(cmap, colorant"green")
        elseif i_mineral in primary_minerals 
            push!(cmap, colorant"grey")
        elseif i_mineral in cation_rich_clays
            push!(cmap, colorant"orange")
        elseif i_mineral in authigenic_minerals
            push!(cmap, colorant"yellow")
        end
    end 
    return cmap 
end  
function animate(plot_function, plot_type)
    cd(animation_directory)
    if plot_type in readdir()
    else
        mkdir(plot_type)
        println("creating ", plot_type)
    end
    cd(plot_type)
    starting_file_list = readdir()
    if length(starting_file_list) > 0
        rm(starting_file_list[end])
    end
    starting_file_list = readdir()
    image_number = 0
    ages = [animation_final_age]
    for age in animation_initial_age:-main_time_step*animation_n_step:animation_final_age
        main_time_step*animation_n_step:animation_final_age
        main_time_step*animation_n_step:animation_final_age
        push!(ages, age)
    end
    for age in ages
        image_number += 1
        image_file_name = "img." * lpad(image_number, 3, "0") * ".png"
        if image_file_name in starting_file_list
            println("already done ", age, " ", plot_type * "/" * image_file_name)
        else
            println("creating ", age, " ", plot_type * "/" * image_file_name)
            global world = read_world(age)
            scene = plot_function()
            Makie.save(image_file_name, scene)
        end
    end
    mp4_file = "../" * plot_type * "." * output_tag * ".mp4"
    println("compiling ", mp4_file)
    rm(mp4_file, force=true)
    run(`ffmpeg -r 2 -f image2 -s 1920x1080 -i img.%03d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p $mp4_file`)
    cp("img.001.png", "../" * plot_type * "." * output_tag * ".png", force=true)
    for imgfile in readdir()
        rm(imgfile)
    end
    cd(animation_directory)
    rm(plot_type, force=true)
    cd(code_base_directory)
end

# general utilities
function min_max_array(a)
    min_val = a[end]; max_val = a[end]
    min_pos = size(a)[1]; max_pos = min_pos
    for i in size(a)[1]-1:-1:1
        if a[i] <= min_val
            min_val = a[i]
            min_pos = i 
        end
        if a[i] >= max_val 
            max_val = a[i]
            max_pos = i
        end 
    end 
    return min_val,max_val, min_pos,max_pos
end 
function derive_lithium_porewater_isotopes(Mg_conc, flushing_rate,
    mineral_volume_step_reactions, mineral_tracer_volume_fractions,
    enable_lithium_reprecipitation, time_step)
    #= i_layer = 1
    Mg_conc = box_solute_concentrations[Mg_2plus_solute,i_layer,i_box] 
    flushing_rate = box_fluid_flushing_rates[i_layer,i_box] 
    mineral_volume_step_reactions = box_mineral_meters_for_residual_array[:,i_layer,i_box]
    mineral_tracer_volume_fractions = 
        box_mineral_tracer_volume_fractions[:,:,i_layer,i_box] =# 

    mineral_tracer_volume_fluxes = fill(0.,n_minerals,n_tracers)
    # porewater [Li] from [Mg] correlation 
    Li_conc = 
        Mg_conc * dissolved_Li_2_Mg_ratio 
    runoff_Li_flux = 
        Li_conc * flushing_rate # mol / m2 yr
    dissolving_Li_flux = 0.   #dissolving_Mg_flux = 0.
    for i_mineral in lithium_bearing_minerals # calculate dissolution 
        if mineral_volume_step_reactions[i_mineral] < 0. # dissolving
            #println("dissolving ", mineral_names[i_mineral])
            box_dissolving_Li_flux =
                - mineral_volume_step_reactions[i_mineral] * # m / step 
                mineral_tracer_volume_fractions[i_mineral,Li_tracer] * # m3 Li / m2 step 
                Li_mol_per_m3 * # mol Li / m2 step 
                1. / time_step # mol Li / yr 
            dissolving_Li_flux += box_dissolving_Li_flux
            mineral_tracer_volume_fluxes[i_mineral,Li_tracer] = 
                mineral_volume_step_reactions[i_mineral] * # m / step 
                mineral_tracer_volume_fractions[i_mineral, Li_tracer] # m3 Li / m2 step 
                #1.0 / time_step # m3 Li / yr 
        end 
    end
    precipitating_Li_rate_fractions = fill(0.,n_minerals)
    total_precipitating_Li_bearing_mineral_rates = 0.
    for i_mineral in lithium_bearing_secondary_minerals # calculate uptake fractions 
        if mineral_volume_step_reactions[i_mineral] > 0. # precipitating
            #println("precipitating " * mineral_names[i_mineral])
            potential_ride = 
                mineral_volume_step_reactions[i_mineral] *
                mineral_relative_Li_affinity[i_mineral]
            precipitating_Li_rate_fractions[i_mineral] = potential_ride
            total_precipitating_Li_bearing_mineral_rates += potential_ride
        end 
    end 
    if total_precipitating_Li_bearing_mineral_rates > 0
        precipitating_Li_rate_fractions ./= 
            total_precipitating_Li_bearing_mineral_rates
    end 
    if runoff_Li_flux > dissolving_Li_flux # all dissolved 
        #println("runoff Li exceeds dissolution sources ")
        Li_conc *= 
            dissolving_Li_flux /
            runoff_Li_flux
        runoff_Li_flux = dissolving_Li_flux
    end 
    #if enable_lithium_dissolution == false
    #    runoff_Li_flux = 0.
    #    box_solute_concentrations[Li_solute,i_layer,i_box] = 0.
    #end 
    if enable_lithium_reprecipitation == false
        precipitating_Li_rate_fractions .= 0.
        Li_conc *= 
            dissolving_Li_flux /
            runoff_Li_flux
        runoff_Li_flux = dissolving_Li_flux
    end 
    precipitating_Li_rates = precipitating_Li_rate_fractions .*
        (dissolving_Li_flux - runoff_Li_flux)
        # lithium mass balance in a box but what about saprolite erosion? 

    #box_mineral_tracer_dissolution_flux[:,Li_tracer,i_layer,i_box] .-=
    #    precipitating_Li_rates

    # compute fractionated porewater and tracer RHS for transport 
    # diss = runoff + pcp1 + pcp2
    # d_pcp1 = d_runoff + D_pcp1 
    # d_pcp2 = d_runoff + D_pcp2
    # diss * d_diss = runoff * d_runoff + pcp1 * (d_runoff + D_pcp1) + pcp2 * (d_runoff + D_pcp2)
    # d_runoff (runoff + pcp1 + pcp2) = diss * d_diss - pcp1 * D_pcp1 - pcp2 * D_pcp2
    # d_runoff = (diss * d_diss - pcp1 * D_pcp1 - pcp2 * D_pcp2) / (runoff + pcp1 + pcp2)
    
    dissolved_del_Li7 = 0. # d_diss of bedrock = 0
    denominator = runoff_Li_flux
    for i_mineral in lithium_bearing_secondary_minerals 
        dissolved_del_Li7 -= 
            precipitating_Li_rates[i_mineral] *
            mineral_Li_fractionation[i_mineral] 
        denominator += precipitating_Li_rates[i_mineral]
    end 
    if denominator > 0
        dissolved_del_Li7 /= denominator
    end 

    Li7_conc = 
        Li_conc * 
        ( 1. + dissolved_del_Li7 / 1.e3 )

    mineral_tracer_volume_fluxes = fill(0.,n_minerals,n_tracers)
    for i_mineral in lithium_bearing_secondary_minerals
        if mineral_volume_step_reactions[i_mineral] > 0 # precipitating 
            mineral_tracer_volume_fluxes[i_mineral,Li_tracer] = 
                ( dissolving_Li_flux - runoff_Li_flux ) * # mol Li / m2 yr 
                precipitating_Li_rate_fractions[i_mineral] *  # still mol Li / m2 year
                1 / Li_mol_per_m3  * # m3 Li / m2 yr 
                time_step # meters Li / step 
            mineral_del_Li7 = dissolved_del_Li7 +
                mineral_Li_fractionation[i_mineral]
            mineral_tracer_volume_fluxes[i_mineral,Li7_tracer] = 
                mineral_tracer_volume_fluxes[i_mineral,Li_tracer] * 
                ( 1. + mineral_del_Li7 / 1.e3 )
        else # dissolving 
            for i_tracer in [Li_tracer,Li7_tracer]
                mineral_tracer_volume_fluxes[i_mineral,i_tracer] = 
                    mineral_volume_step_reactions[i_mineral] *
                    mineral_tracer_volume_fractions[i_mineral,i_tracer]
            end 
        end 
    end 

    return Li_conc, Li7_conc, dissolved_del_Li7, mineral_tracer_volume_fluxes
end 
function update_mineral_volumes(mineral_fractions, layer_thickness, porosity)
    solid_thickness = layer_thickness * ( 1. - porosity )
    mineral_volumes = mineral_fractions .* solid_thickness
    return mineral_volumes
end 
function update_mineral_fractions(mineral_volumes, porosity)
    solid_volume = 0.
    for i_mineral in 1:n_minerals
        solid_volume += mineral_volumes[i_mineral]
    end 
    mineral_fractions = mineral_volumes ./ solid_volume
    layer_thickness = solid_volume / (1.0 - porosity)
    return mineral_fractions, solid_volume, layer_thickness
end 
function update_tracer_volumes(tracer_fractions, layer_thickness, porosity)
    solid_thickness = layer_thickness * ( 1. - porosity )
    mineral_volumes = mineral_fractions .* solid_thickness
    return mineral_volumes
end 
function update_tracer_fractions(mineral_volumes, tracer_volumes)
    #= tracer_volumes = box_mineral_tracer_volumes[:,1,2,1]
    mineral_volumes = box_mineral_volumes[:,2,1] =#
    volume_fractions = fill(0.,n_minerals)
    mole_fractions = fill(0.,n_minerals)
    for i_mineral in 1:n_minerals 
        if mineral_volumes[i_mineral] > 0
            volume_fractions[i_mineral] = 
                tracer_volumes[i_mineral] / 
                mineral_volumes[i_mineral]
            mole_fractions[i_mineral] = 
                volume_fractions[i_mineral] *
                mineral_molwts[i_mineral] ./ 6.
        end 
    end 
    return volume_fractions, mole_fractions
end
function mineral_fluxes_to_components(mineral_fluxes)
    component_flux = fill(0.,n_reaction_components)
    for (i_mineral, flux) in pairs(mineral_fluxes)
        #println(i_mineral," ",flux)
        for i_component in 1:n_reaction_components
            component_flux[i_component] +=
                flux * # mol / l pw 
                mineral_component_stoic[i_mineral, i_component]
        end 
        #println( )
    end 
    return component_flux
end 
function components_to_river_concentrations()
    surface_temperature = 15.
    soil_CO2 = 4000.
    component_reaction_flux = 
        [2.0e-4,8.0e-4,1.0e-4,1.5e-5,1.0e-4,1.0e-8,1.0e-8,1.0e-8,0.0,0.0]
    #      CaO,    MgO,  Na2O,   K2O,  SiO2, Al2O3,   FeO, Fe2O3,CO2,H2O
    generate_phreeqc_input_file_from_components(
        surface_temperature, soil_CO2, component_reaction_flux, true)
    cmd = `$exec_command`
    run(cmd)
    solute_concentrations, # mol / m3 pw 
    mineral_saturation_indices,
    mineral_reaction_extents = # mol / m3 pw
        parse_phreeqc_output_file()
    print_table(solute_names,solute_concentrations,
        typical_river_solute_concentrations)
end 

function base_cation_abundance_report(mineral_list, mineral_inventories)
    #= mineral_list = primary_minerals
    mineral_inventories = layer_mineral_fractions[:,1] =#
    base_cation_inventories = [0.0, 0.0]
    total_inventory = 0
    for i_mineral in mineral_list
        total_inventory += mineral_inventories[i_mineral]
        for i_component in 1:2
            base_cation_inventories[i_component] +=
                mineral_inventories[i_mineral] * # m3 / m2 year
                mineral_component_stoic[i_mineral, i_component] # moles of component /m2 year
        end
    end
    base_cation_fractions = base_cation_inventories ./ total_inventory
    report_string = "Ca/Mg " *
        string(Int(floor(base_cation_fractions[1] * 100))) * "/" *
        string(Int(floor(base_cation_fractions[2] * 100)))
    return report_string
end
function clay_maturity_report(mineral_fractions)
    # mineral_fractions = box_mineral_fractions[:,2,1]
    total_cation_rich_clays = 0.
    for i_mineral = cation_rich_clays
        total_cation_rich_clays += mineral_fractions[i_mineral]
    end 
    total_cation_poor_clays = 
        mineral_fractions[Kaolinite_mineral] +
        mineral_fractions[Gibbsite_mineral]
    total_clays = total_cation_rich_clays + total_cation_poor_clays
    if total_clays == 0.
        return ""
    end 
    rich_fraction = total_cation_rich_clays / 
        total_clays * 100.0
    #poor_fraction = 100 - rich_fraction 
    total_authigenic_fraction = 0.
    for i_mineral in authigenic_minerals
        total_authigenic_fraction += mineral_fractions[i_mineral]
    end 
    report_string = "clay richness " *
        string(Int(floor(rich_fraction))) * "% " * "authigenic " *
        string(Int(floor(total_authigenic_fraction * 100))) * "% "
    return report_string
end 
function print_table(name_list, value_list)
    n_lines = length(name_list)
    for i_line in 1:n_lines
        println(name_list[i_line], " ", value_list[i_line])
    end 
end
function print_table(f,name_list, value_list)
    n_lines = length(name_list)
    for i_line in 1:n_lines
        println(f,name_list[i_line], " ", value_list[i_line])
    end 
end
function print_table(name_list, value_list, ref_value_list)
    n_lines = length(name_list)
    for i_line in 1:n_lines
        if ref_value_list[i_line] == ref_value_list[i_line]
            println(name_list[i_line], " ", value_list[i_line], " ",ref_value_list[i_line])
        end
    end
end
function print_table(f,name_list, value_list, ref_value_list)
    n_lines = length(name_list)
    for i_line in 1:n_lines
        if ref_value_list[i_line] == ref_value_list[i_line]
            println(f,name_list[i_line], " ", value_list[i_line], " ",ref_value_list[i_line])
        end
    end
end
function vector_total(vector)
    total = 0.
    for i in 1:length(vector)
        total += vector[i]
    end 
    return total 
end
function mineral_table(mineral_field)
    #n_columns = size(mineral_field)[2]
    for i_mineral in 1:n_minerals
        println(mineral_names[i_mineral],"  ",mineral_field[i_mineral,:])
    end
end
function mineral_sorted_table(mineral_list,mineral_fractions)
    #= mineral_list = secondary_silicate_minerals
    mineral_fractions = mineral_reaction_extents[:, 1] =#
    #println(mineral_list)
    edited_fraction_list = []
    edited_id_list = []
    values_total = 0.
    for i_mineral in mineral_list
        append!(edited_fraction_list, mineral_fractions[i_mineral])
        append!(edited_id_list, i_mineral)
        values_total += mineral_fractions[i_mineral]
    end 
    permutation_vector = sortperm(edited_fraction_list, rev=true)
    for i_p in permutation_vector
        println(mineral_names[edited_id_list[i_p]], "  ",
            Int(floor(edited_fraction_list[i_p]/values_total *100.)))
    end
end
function mineral_sorted_table(f,mineral_list, mineral_fractions)
    #= mineral_list = secondary_silicate_minerals
    mineral_fractions = mineral_reaction_extents[:, 1] =#
    #println(mineral_list)
    edited_fraction_list = []
    edited_id_list = []
    values_total = 0.0
    for i_mineral in mineral_list
        append!(edited_fraction_list, mineral_fractions[i_mineral])
        append!(edited_id_list, i_mineral)
        values_total += mineral_fractions[i_mineral]
    end
    permutation_vector = sortperm(edited_fraction_list, rev=true)
    for i_p in permutation_vector
        println(f,mineral_names[edited_id_list[i_p]], "  ",
            Int(floor(edited_fraction_list[i_p] / values_total * 100.0)))
    end
end
function mineral_abundance_report(mineral_list,mineral_inventories)
    #= mineral_list = secondary_silicate_minerals
    mineral_inventories = mineral_reaction_extents[:, 2] =#
    most_abundant_minerals = [0,0]
    most_abundant_fractions = [0.,0.]
    total_inventory = 0.
    mineral_fractions = fill(0.,n_minerals)
    for i_mineral in mineral_list
        total_inventory += mineral_inventories[i_mineral]
        if mineral_inventories[i_mineral] > most_abundant_fractions[1]
            most_abundant_minerals[1] = i_mineral
            most_abundant_fractions[1] = mineral_inventories[i_mineral]
        end 
        #println(i_mineral, [most_abundant_minerals[1], most_abundant_fractions[1]])
    end 
    for i_mineral in mineral_list
        #println(i_mineral)
        mineral_fractions[i_mineral] =
            mineral_inventories[i_mineral] /
            total_inventory
        if i_mineral != most_abundant_minerals[1]
            if mineral_inventories[i_mineral] > most_abundant_fractions[2]
                most_abundant_minerals[2] = i_mineral
                most_abundant_fractions[2] = mineral_inventories[i_mineral]
            end 
        end 
        #println(i_mineral, [most_abundant_minerals[2], most_abundant_minerals[2]])
    end 
    report_string = mineral_names[most_abundant_minerals[1]] * "/" * 
        mineral_names[most_abundant_minerals[2]] * " " *
        string(Int(floor(mineral_fractions[most_abundant_minerals[1]]*100))) *
        "/" *
        string(Int(floor(mineral_fractions[most_abundant_minerals[2]]*100)))
    #println(report_string)
    return report_string
end
function next_line_number(file_lines, current_line, match_string)
    still_looking_for_match_string = true 
    match_line_number = 0
    iline = current_line
    while still_looking_for_match_string 
        file_line = file_lines[iline]
        if findlast(match_string,file_line) == nothing && 
            size(file_lines)[1] > iline
                iline += 1 
        else 
            still_looking_for_match_string = false 
            match_line_number = iline
        end 
    end 
    return match_line_number
end 
function extract_value_from_line(file_line,word_position)
    words = split(file_line)
    extracted_value = parse(Float64, words[word_position])
    return extracted_value
end 
function extract_value_from_next_line_number(
    file_lines, current_line, match_string, word_position)
    iline = next_line_number(file_lines, current_line, match_string)
    value = extract_value_from_line(file_lines[iline],word_position)
    return value, iline
end 
function extract_values_from_lines(file_lines,iline,
    exit_string, 
    variable_names,label_position,word_position)
    n_names = size(variable_names)[1]
    values_list = fill(0.,n_names)
    still_finding_values = true
    while still_finding_values
        iline += 1 
        file_line = file_lines[iline]
        if findlast(exit_string, file_line) !== nothing 
            still_finding_values = false 
            #println(still_finding_values)
        else
            words = split(file_line)
            if size(words)[1] >= word_position
                if words[label_position] in variable_names
                    value_index = findfirst(isequal(words[label_position]), variable_names)
                    values_list[value_index] = parse(Float64, words[word_position])
                    #println(iline," ",words[label_position])
                end 
                
            end
            #println(iline," ",file_line)
        end
    end
    return values_list, iline
end 
function mineral_phreeqc_rate_from_porewater_rate(
    mineral_porewater_reaction_rate,
    #layer_thickness,porosity,
    flushing_age)
    # phreeqc_rate = mol / l pore fluid at steady state with flushing 
    # porewater_rate = mol / m3 porewater yr 
    # positive means dissolution of a solid, src to pw 
    mineral_phreeqc_rate = # because phreeqc rate is a solid rate 
        - mineral_porewater_reaction_rate * # mol / m3 pw year 
        flushing_age * # mol / m3 pw
        1.e-3 # mol / l pw 
        #porosity * # mol / m3 bulk yr 
        #layer_thickness *  # mol / m2 yr 
    return mineral_phreeqc_rate
end 
function mineral_porewater_rate_from_phreeqc_rate(
    mineral_phreeqc_reaction_rate,
    #layer_thickness, porosity, 
    flushing_age)
    mineral_porewater_rate = - mineral_phreeqc_reaction_rate * # mol/l pw
        #1. / layer_thickness *
        #1. / porosity *
        1. / flushing_age * # mol / l yr 
        1.e3 # mol / m3 yr 
    return mineral_porewater_rate
end
function mineral_solid_rate_from_porewater_rate(
    mineral_porewater_rate, mineral_mol_conc, porosity, layer_thickness)
    # solid rate = meters mineral / yr 
    # sign switch, negative solid flux = positive porewater flux 
    mineral_solid_rate = - mineral_porewater_rate * # mol / m3 porewater yr 
        porosity * # mol mineral / m3 bulk yr 
        1. / mineral_mol_conc * # m3 mineral / m3 bulk yr 
        layer_thickness # m3 mineral / m2 yr = m / yr 
    return mineral_solid_rate
end 
function mineral_porewater_rate_from_solid_rate(
    mineral_solid_rate, mineral_mol_conc, porosity, layer_thickness)
    #= mineral_solid_rate = mineral_solid_reaction_rates[1]
    mineral_mol_conc = mineral_mol_per_m3[1] 
    =#
    mineral_porewater_rate = # mol / m3 porewater yr 
        - mineral_solid_rate * # m3 mineral / m2 yr 
        1. / layer_thickness * # m3 mineral / m3 bulk yr 
        mineral_mol_conc * # mol / m3 bulk yr 
        1. / porosity # mol / m3 pw yr 
    return mineral_porewater_rate
end
function mineral_phreeqc_rate_from_solid_rate(
    mineral_solid_rate, mineral_mol_conc, porosity, layer_thickness,flushing_age)
    mineral_porewater_rate = mineral_porewater_rate_from_solid_rate(
        mineral_solid_rate, mineral_mol_conc, porosity, layer_thickness)
    mineral_phreeqc_rate = mineral_phreeqc_rate_from_porewater_rate(
        mineral_porewater_rate,
        #layer_thickness,porosity,
        flushing_age)
    return mineral_phreeqc_rate
end 
function mineral_solid_rate_from_phreeqc_rate(
    mineral_phreeqc_rate, mineral_mol_conc, porosity, layer_thickness, flushing_age)
    mineral_porewater_rate = mineral_porewater_rate_from_phreeqc_rate(
        mineral_phreeqc_rate, flushing_age)
    mineral_solid_rate = mineral_solid_rate_from_porewater_rate(
        mineral_porewater_rate,
        mineral_mol_conc, porosity, layer_thickness)
    return mineral_solid_rate
end
function component_fluxes_from_interpolation_values( 
    index_sizes, index_range_lists, fraction_lists )

    component_reaction_fluxes = fill(0.0, n_reaction_components)
    for i_frac in 1:index_sizes[7] # SiO2
        component_reaction_fluxes[SiO2_component] += 
            phreeqc_SiO2_reaction_extents[index_range_lists[7][i_frac]] * 
            fraction_lists[7][i_frac]
    end 
    for i_frac in 1:index_sizes[3] # CaO
        component_reaction_fluxes[CaO_component] += 
            phreeqc_CaO_reaction_extents[index_range_lists[3][i_frac]] * 
            component_reaction_fluxes[SiO2_component] * 
            phreeqc_component_silica_ratio_ranges[CaO_component] *
            fraction_lists[3][i_frac]
    end 
    for i_frac in 1:index_sizes[4]
        component_reaction_fluxes[MgO_component] += 
            phreeqc_MgO_reaction_extents[index_range_lists[4][i_frac]] * 
            component_reaction_fluxes[SiO2_component] * 
            phreeqc_component_silica_ratio_ranges[MgO_component] *
            fraction_lists[4][i_frac]
    end 
    for i_frac in 1:index_sizes[5]
        component_reaction_fluxes[Na2O_component] += 
            phreeqc_Na2O_reaction_extents[index_range_lists[5][i_frac]] * 
            component_reaction_fluxes[SiO2_component] * 
            phreeqc_component_silica_ratio_ranges[Na2O_component] *
            fraction_lists[5][i_frac]
    end 
    for i_frac in 1:index_sizes[6]
        component_reaction_fluxes[K2O_component] += 
            phreeqc_K2O_reaction_extents[index_range_lists[6][i_frac]] * 
            component_reaction_fluxes[SiO2_component] * 
            phreeqc_component_silica_ratio_ranges[K2O_component] *
            fraction_lists[6][i_frac]
    end 
    for i_frac in 1:index_sizes[8]
        component_reaction_fluxes[Al2O3_component] += 
            phreeqc_Al2O3_reaction_extents[index_range_lists[8][i_frac]] * 
            component_reaction_fluxes[SiO2_component] * 
            phreeqc_component_silica_ratio_ranges[Al2O3_component] *
            fraction_lists[8][i_frac]
    end 
    for i_frac in 1:index_sizes[9]
        component_reaction_fluxes[CO2_component] += 
            phreeqc_CO2_alk_reaction_ratios[index_range_lists[9][i_frac]] * 
            ( 2 * component_reaction_fluxes[CaO_component] + 
              2 * component_reaction_fluxes[MgO_component] + 
              component_reaction_fluxes[Na2O_component] +
              component_reaction_fluxes[K2O_component] + 
              3. * component_reaction_fluxes[Al2O3_component]) * 
            fraction_lists[9][i_frac]
    end 





    #for i_comp in [1,2,3,4,6] 
    #    for i_frac in 1:index_sizes[i_comp+2]
    #        component_reaction_fluxes[i_comp] += 
    #            component_reaction_fluxes[SiO2_component] *
    #            phreeqc_silicate_cation_reaction_extents[index_range_lists[i_comp+2][i_frac]] * 
    #            phreeqc_component_silica_ratio_ranges[i_comp] * 
    #            fraction_lists[i_comp+2][i_frac]
    #    end 
    #end 
    #carbonate_omegas = fill(0.,Calcite_mineral:Dolomite_mineral)
    ##for i_carb in 1:2
    #    for i_frac in 1:index_sizes[i_carb+8]
    #        carbonate_omegas[Calcite_mineral+i_carb-1] += 
    #            phreeqc_carbonate_omegas[index_range_lists[i_carb+8][i_frac]] *
    #            fraction_lists[i_carb+8][i_frac]
    ##    end 
    #end 

    #for i_comp in [MgO_component,K2O_component,Al2O3_component] 
    #    for i_frac in 1:index_sizes[i_comp+2]
    #        component_reaction_fluxes[i_comp] += 
    #            component_reaction_fluxes[SiO2_component] *
    #            phreeqc_silicate_cation_reaction_extents[index_range_lists[i_comp+2]][i_frac] * 
    #            phreeqc_component_silica_ratio_ranges[i_comp] *
    #            fraction_lists[i_comp+2][i_frac]
    #    end 
    #end 
    #i_comp = CO2_component
    #=for i_frac in 1:index_sizes[CO2_component+2]
        component_reaction_fluxes[CO2_component] +=
            component_reaction_fluxes[SiO2_component] *
            phreeqc_carbonate_to_silicate_ratios[index_range_lists[CO2_component+2]][i_frac] *
            fraction_lists[CO2_component+2][i_frac]
        component_reaction_fluxes[CaO_component] +=
            component_reaction_fluxes[SiO2_component] *
            phreeqc_carbonate_to_silicate_ratios[index_range_lists[CO2_component+2]][i_frac] *
            fraction_lists[CO2_component+2][i_frac]
    end =#
    return component_reaction_fluxes # , carbonate_omegas
end 
function test_phreeqc_interpolation() 
    #  itemp, iCO2, ic1, ic2, ic3, ic4, ic5, ic6, icaco3 = 1, 2, 5, 5, 5, 5, 15, 5, 3
    component_reaction_fluxes = fill(0.,n_reaction_components)
    for i_component in [SiO2_component]
        component_reaction_fluxes[i_component] = #rand() * 
            phreeqc_SiO2_reaction_extents[9]
    end 
    for i_component in [CaO_component, MgO_component, 
        Na2O_component, K2O_component, Al2O3_component]
        component_reaction_fluxes[i_component] = #rand() *
            component_reaction_fluxes[SiO2_component] *
            phreeqc_silicate_cation_reaction_extents[3] *
            phreeqc_component_silica_ratio_ranges[i_component]
    end 
    #=alkalinity_reaction_flux = ( 
        2 * component_reaction_fluxes[CaO_component] + 
        2 * component_reaction_fluxes[MgO_component] + 
        component_reaction_fluxes[Na2O_component] +
        component_reaction_fluxes[K2O_component] + 
        3. * component_reaction_fluxes[Al2O3_component]
        )
    component_reaction_fluxes[CO2_component] = phreeqc_CO2_alk_reaction_ratios[1] * 
        alkalinity_reaction_flux=#
    #carbonate_omegas = fill(1., Calcite_mineral:Dolomite_mineral)
    #=for i_mineral in Calcite_mineral:Dolomite_mineral
        carbonate_omegas[i_mineral] = rand()
    end =#
    #carbonate_reaction_extents = fill(0., Calcite_mineral:Dolomite_mineral)
    #for i_mineral in Calcite_mineral:Dolomite_mineral
    #    carbonate_reaction_extents[i_mineral] = #rand() * 
    #        - phreeqc_carbonate_reaction_extents[12]
    #end 
    #carbonate_reaction_extents[Dolomite_mineral] /= 10.
    #for i_mineral in Calcite_mineral:Dolomite_mineral
    #    for i_component in 1:n_reaction_components
    #        component_reaction_fluxes[i_component] -=
    #            carbonate_reaction_extents[i_mineral] *
    #            mineral_component_stoic[i_mineral, i_component]
    #    end
    #end 
    #carbonate_omegas[Calcite_mineral] = rand()
    #carbonate_omegas[Dolomite_mineral] = rand()

    surface_temperature = 10.
    soil_CO2 = 5000.
    #=
    interpolated_output_block = 
        interpolate_phreeqc_lattice(
            surface_temperature, soil_CO2, 
            component_reaction_fluxes, 
            carbonate_omegas, 
            carbonate_reaction_extents,
            false) =#
    #fake_interpolated_output_block =
    #    interpolate_phreeqc_lattice(
    #        surface_temperature, soil_CO2,
    #        component_reaction_fluxes, 
    #        carbonate_omegas,
    #        carbonate_reaction_extents, 
    #        true) # , carbonate_omegas, true)

    interp_output_block = interpolate_sparse_phreeqc_lattice( # with carbonate omega 
        surface_temperature, soil_CO2, component_reaction_fluxes, 
        #n_memories_of_phreeqc_filled,
        true, true )
    direct_output_block =
        run_phreeqc(surface_temperature, soil_CO2,
            component_reaction_fluxes, 
            #carbonate_omegas,
            #carbonate_reaction_extents, 
            #enable_default_solid_solutions,
            1)
    compare_two_output_blocks(output_block_int, output_block,
        [1.,1.,1.],.5)
end 
function compare_two_output_blocks(block1, block2, misfits, tolerance)
    is_ok = true 
    if misfits[1] > tolerance
        println("Conc misfit ", misfits[1])
        for i_solute in 1:n_solutes
            println(solute_names[i_solute], 
                [block1.solute_concentrations[i_solute],
                block2.solute_concentrations[i_solute]])
        end 
        is_ok = false
        println()
    end 
    if misfits[2] > tolerance
        println("Rxn misfit ", misfits[2])
        for i_mineral in secondary_minerals
            println(mineral_names[i_mineral],[
                block1.mineral_reaction_extents[i_mineral],
                block2.mineral_reaction_extents[i_mineral],])
        end
        is_ok = false
        println()
    end 
    if misfits[3] > tolerance
        println("Omega misfit ", misfits[3])
        for i_mineral in secondary_minerals
            println(mineral_names[i_mineral], [
                #10^interpolated_output_block.mineral_saturation_indices[i_mineral],
                10^block1.mineral_saturation_indices[i_mineral],
                10^block2.mineral_saturation_indices[i_mineral]])
        end
        is_ok = false
        println()
    end 
    if is_ok == false
        println("")
    end 
end
function search_sorted_below(a, x)
    ijustbelow = NaN 
    for i in 1:length(a)
        if a[1] < a[end] # for example xcoords
            if a[i] <= x
                ijustbelow = i
            end
        else
            reverse_i = length(a) - i + 1
            if a[reverse_i] <= x
                ijustbelow = reverse_i
            end
        end
    end
    return ijustbelow
end
function find_fill_fractions(a, x) # returns [ilow,ihigh],[lowfrac,highfrac]
    # a = interp_domains[5][1]; x = interp_domains[5][2]
    len_a = size(a)[1]
    indices = []
    fractions = []
    index_below = 0
    if x == x 
        if a[1] == a[1]
            index_below = search_sorted_below(a, x)
            if index_below < 1
                error("extrapolating low",a,x)
            elseif index_below == len_a 
                error("extrapolating high",a,x)
            end 
        else
            index_below = search_sorted_below(a[2:end], x) + 1
            println("how did I get here?")
        end 
        push!(indices, index_below)
        #indices[2] = indices[1] + 1
        push!(fractions,1.) 
        if x > a[index_below]
            if index_below < len_a 
                frac_below = 
                    (a[index_below+1] - x) /
                    (a[index_below+1] - a[index_below])
                if frac_below < 1.0
                    push!(indices, index_below + 1)
                    fractions[1] = frac_below
                    push!(fractions, 1.0 - frac_below)
                end 
            end
        end 
    else # x = NaN?  
        push!(fractions,1.)
        push!(indices,1)
        #println("wassup now? ")
    end 
    return size(indices)[1], indices, fractions
end
