"""
Spline loading and interpolation utilities for cross-sections.
"""

"""
    load_spline(filepath::String)

Load a spline from a JLD2 file.

The file should contain:
- `knots`: Array of knot points
- `values`: Array of values at knot points
- `order`: (optional) Spline order

Returns an interpolation object from Interpolations.jl.
"""
function load_spline(filepath::String)
    data = load(filepath)

    knots = Float64.(data["knots"])
    values = Float64.(data["values"])

    # Create linear interpolation (can upgrade to cubic if needed)
    return linear_interpolation(knots, values, extrapolation_bc=Line())
end

"""
    load_spline_2d(filepath::String)

Load a 2D spline from a JLD2 file.

The file should contain:
- `knots_x`: Array of knot points for first dimension
- `knots_y`: Array of knot points for second dimension
- `values`: 2D array of values

Returns a 2D interpolation object.
"""
function load_spline_2d(filepath::String)
    data = load(filepath)

    knots_x = Float64.(data["knots_x"])
    knots_y = Float64.(data["knots_y"])
    values = Float64.(data["values"])

    return linear_interpolation((knots_x, knots_y), values, extrapolation_bc=Line())
end

"""
    save_spline(filepath::String, knots::AbstractVector, values::AbstractVector)

Save a 1D spline to a JLD2 file.
"""
function save_spline(filepath::String, knots::AbstractVector, values::AbstractVector)
    jldsave(filepath; knots=collect(knots), values=collect(values))
end

"""
    save_spline_2d(filepath::String, knots_x, knots_y, values::AbstractMatrix)

Save a 2D spline to a JLD2 file.
"""
function save_spline_2d(filepath::String, knots_x, knots_y, values::AbstractMatrix)
    jldsave(filepath;
            knots_x=collect(knots_x),
            knots_y=collect(knots_y),
            values=collect(values))
end
