"""
`init_subspec!(m::COTHANK)`

Initializes a model subspecification by overwriting parameters from
the original model object with new parameter objects. This function is
called from within the model constructor.
"""
function init_subspec!(m::COTHANK)

    if subspec(m) == "ss1"
        return
        error("This subspec is not defined.")
    end
end
