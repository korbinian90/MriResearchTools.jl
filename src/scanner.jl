const SCANNER = Dict(   
    :B7T => Dict(
        :γ => 42.5756e6, # in [MHz/T]
        :SlewR => 200, # in [(T/m)/s] => 200 max
        :Glimit => 42e-3, # in [T/m] => 70 max
        :Texc => 2.1, # in [ms]
        :Tspoil => 3.0, # in [ms]
    ),
    :B3T => Dict(
        :γ => 42.5756e6, # in [MHz/T]
    ),
)
