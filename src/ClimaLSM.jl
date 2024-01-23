module ClimaLSM
using UnPack
using DocStringExtensions
using Hyb_LSM
using ClimaCore
import ClimaCore: Fields, Spaces

include("shared_utilities/Parameters.jl")
import .Parameters as LSMP

include("shared_utilities/Regridder.jl")
include("shared_utilities/Domains.jl")
include("shared_utilities/FileReader.jl")
using .Domains
include("shared_utilities/utils.jl")
include("shared_utilities/models.jl")
include("shared_utilities/drivers.jl")
include("shared_utilities/boundary_conditions.jl")
include("shared_utilities/sources.jl")
include("shared_utilities/implicit_tendencies.jl")
include("shared_utilities/implicit_functions.jl")
include("standalone/Bucket/Bucket.jl")
export make_interactions_update_aux
"""
     AbstractLandModel{FT} <: AbstractModel{FT}

An abstract type for all land model types, which are used
to simulated multiple land surface components as
a single system. standalone component runs do not require
this interface and it should not be used for that purpose.

Many methods taking an argument of type `AbstractLandModel` are
extensions of functions defined for `AbstractModel`s.
There are default methods that apply for all `AbstractLandModel`s,
including `make_update_aux`, `make_exp_tendency`, `make_imp_tendency`,
`make_compute_exp_tendency`, `make_compute_imp_tendency`,
`initialize_prognostic`, `initialize_auxiliary`, `initialize`,
and `coordinates`.

Methods which dispatch on a specific type of AbstractLandModel
include any function involving interactions between components,
as these interactions depend on the components in the land model
and the versions of these component models being used.
"""
abstract type AbstractLandModel{FT} <: AbstractModel{FT} end

ClimaLSM.name(::AbstractLandModel) = :land

"""
    Domains.coordinates(model::AbstractLandModel)

Returns a NamedTuple of the unique set of coordinates for the LSM
`model`, where the unique set is taken over the coordinates of all
of the subcomponents.

For example, an LSM with a single layer snow model, multi-layer 
soil model, and canopy model would have a coordinate set corresponding
to the coordinates of the surface (snow), the subsurface coordinates (soil)
and the coordinates of the surface (canopy). This would return the coordinates
of the surface and subsurface. These are distinct because the subsurface
coordinates correspond to the centers of the layers, while the surface
corresponds to the top face of the domain.
"""
function Domains.coordinates(model::AbstractLandModel)
    components = land_components(model)
    coords_list = map(components) do (component)
        Domains.coordinates(getproperty(model, component))
    end
    unique_coords_list = merge(coords_list...)
    return unique_coords_list
end


function initialize_prognostic(
    model::AbstractLandModel{FT},
    coords::NamedTuple,
) where {FT}
    components = land_components(model)
    Y_state_list = map(components) do (component)
        submodel = getproperty(model, component)
        getproperty(initialize_prognostic(submodel, coords), component)
    end
    Y = ClimaCore.Fields.FieldVector(; NamedTuple{components}(Y_state_list)...)
    return Y
end

function initialize_auxiliary(
    model::AbstractLandModel{FT},
    coords::NamedTuple,
) where {FT}
    components = land_components(model)
    p_state_list = map(components) do (component)
        submodel = getproperty(model, component)
        getproperty(initialize_auxiliary(submodel, coords), component)
    end
    p_interactions = initialize_interactions(model, coords)
    p = (; p_interactions..., NamedTuple{components}(p_state_list)...)
    domains_list = map(components) do (component)
        submodel = getproperty(model, component)
        getproperty(submodel, :domain)
    end
    for domain in unique(domains_list)
        p = add_dss_buffer_to_aux(p, domain)
    end
    return p
end

"""
    initialize_interactions(land::AbstractLandModel) end

Initializes interaction variables, which are a type of auxiliary
variable, to empty objects of the correct type for the model.

Interaction variables are specified by `interaction_vars`, their types
by `interaction_types`, and their spaces by `interaction_spaces`.
This function should be called during `initialize_auxiliary` step.
"""
function initialize_interactions(land::AbstractLandModel, land_coords)
    vars = interaction_vars(land)
    types = interaction_types(land)
    domains = interaction_domain_names(land)
    interactions = map(zip(types, domains)) do (T, domain)
        zero_instance = zero(T)
        map(_ -> zero_instance, getproperty(land_coords, domain))
    end
    return NamedTuple{vars}(interactions)
end

function make_imp_tendency(land::AbstractLandModel)
    components = land_components(land)

    # If all component models are stepped explicitly, do nothing in imp_tendency!
    if all(c -> typeof(c) .<: AbstractExpModel, components)
        function do_nothing_imp_tendency!(dY, Y, p, t) end
        return do_nothing_imp_tendency!
    else
        compute_imp_tendency_list = map(
            x -> make_compute_imp_tendency(getproperty(land, x)),
            components,
        )
        update_aux! = make_update_aux(land)
        function imp_tendency!(dY, Y, p, t)
            update_aux!(p, Y, t)
            for f! in compute_imp_tendency_list
                f!(dY, Y, p, t)
            end
        end
        return imp_tendency!
    end
end

function make_exp_tendency(land::AbstractLandModel)
    components = land_components(land)
    compute_exp_tendency_list =
        map(x -> make_compute_exp_tendency(getproperty(land, x)), components)
    update_aux! = make_update_aux(land)
    function exp_tendency!(dY, Y, p, t)
        update_aux!(p, Y, t)
        for f! in compute_exp_tendency_list
            f!(dY, Y, p, t)
        end
    end
    return exp_tendency!
end

function make_set_initial_aux_state(land::AbstractLandModel)
    interactions_update_aux! = make_interactions_update_aux(land)
    components = land_components(land)
    # These functions also call update_aux
    set_initial_aux_function_list =
        map(x -> make_set_initial_aux_state(getproperty(land, x)), components)
    function set_initial_aux_state!(p, Y0, t0)
        for f! in set_initial_aux_function_list
            f!(p, Y0, t0)
        end
        interactions_update_aux!(p, Y0, t0) # this has to come last.
    end
    return set_initial_aux_state!
end

function make_update_aux(land::AbstractLandModel)
    interactions_update_aux! = make_interactions_update_aux(land)
    components = land_components(land)
    update_aux_function_list =
        map(x -> make_update_aux(getproperty(land, x)), components)
    function update_aux!(p, Y, t)
        for f! in update_aux_function_list
            f!(p, Y, t)
        end
        interactions_update_aux!(p, Y, t) # this has to come last.
    end
    return update_aux!
end

"""
    make_interactions_update_aux(land::AbstractLandModel) end

Makes and returns a function which updates the interaction variables,
which are a type of auxiliary variable.

The `update_aux!` function returned is evaluted during the right hand
side evaluation.

This is a stub which concrete types of LSMs extend.
"""
function make_interactions_update_aux(land::AbstractLandModel)
    function interactions_update_aux!(p, Y, t)
        nothing
    end
    return interactions_update_aux!
end

"""
    land_components(land::AbstractLandModel)

Returns the component names of the `land` model, by calling
`propertynames(land)`.
"""
land_components(land::AbstractLandModel) = propertynames(land)

function prognostic_vars(land::AbstractLandModel)
    components = land_components(land)
    prognostic_list = map(components) do model
        prognostic_vars(getproperty(land, model))
    end
    return NamedTuple{components}(prognostic_list)
end

function prognostic_types(land::AbstractLandModel)
    components = land_components(land)
    prognostic_list = map(components) do model
        prognostic_types(getproperty(land, model))
    end
    return NamedTuple{components}(prognostic_list)
end

function auxiliary_vars(land::AbstractLandModel)
    components = land_components(land)
    auxiliary_list = map(components) do model
        auxiliary_vars(getproperty(land, model))
    end
    return NamedTuple{(components..., :interactions)}((
        auxiliary_list...,
        interaction_vars(land),
    ))
end

function auxiliary_types(land::AbstractLandModel)
    components = land_components(land)
    auxiliary_list = map(components) do model
        auxiliary_types(getproperty(land, model))
    end
    return NamedTuple{(components..., :interactions)}((
        auxiliary_list...,
        interaction_types(land),
    ))
end


"""
   interaction_vars(m::AbstractLandModel)

Returns the interaction variable symbols for the model in the form of a tuple.
"""
interaction_vars(m::AbstractLandModel) = ()

"""
   interaction_types(m::AbstractLandModel)

Returns the shared interaction variable types for the model in the form of a tuple.
"""
interaction_types(m::AbstractLandModel) = ()

"""
   interaction_domain_names(m::AbstractLandModel)

Returns the interaction domain symbols in the form of a tuple e.g. :surface or :subsurface.

This is only required for variables shared between land submodels, and only needed
for multi-component models, not standalone components. Component-specific variables
should be listed as prognostic or auxiliary variables which do not require this to
initialize.
"""
interaction_domain_names(m::AbstractLandModel) = ()

# Methods extended by the LSM models we support
include("standalone/SurfaceWater/Pond.jl")
using .Pond
import .Pond: surface_runoff
include("standalone/Soil/Soil.jl")
using .Soil
import .Soil: soil_boundary_fluxes
import .Soil.Biogeochemistry: soil_temperature, soil_moisture
include("standalone/Snow/Snow.jl")
using .Snow
include("standalone/Vegetation/Canopy.jl")
using .Canopy
using .Canopy.PlantHydraulics
import .Canopy.PlantHydraulics: root_flux_per_ground_area!
### Concrete types of AbstractLandModels
### and associated methods
include("integrated/soil_energy_hydrology_biogeochemistry.jl")
include("integrated/pond_soil_model.jl")
include("integrated/soil_canopy_model.jl")

end
