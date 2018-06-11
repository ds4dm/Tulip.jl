env = TulipEnv()

# check default value was specified for all parameters
for k in fieldnames(env)
    @test haskey(Tulip.TLP_DEFAULT_PARAM, k)
end

# check all parameters are initialized to default
# this also tests all the getters
for (k, v) in Tulip.TLP_DEFAULT_PARAM
    @test env[k] == Tulip.TLP_DEFAULT_PARAM[k]
    @test env[String(k)] == Tulip.TLP_DEFAULT_PARAM[k]
end

# test setters
d = Dict()