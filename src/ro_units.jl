module ROUnits

export lmh, lmhbar

using Unitful

@unit lmh "lmh" LMH 1.0u"L/m^2/hr" false # Liters per square meter per hour. Unit often used for flux.

@unit lmmhbar "lmhbar" LMHperBar 1.0u"L/m^2/hr/bar" false # Liters per square meter per hour per bar. Unit often used for permeability.

end # module ROUnits