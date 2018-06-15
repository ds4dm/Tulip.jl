env = Tulip.TulipEnv()

# check default value was specified for all parameters
for k in fieldnames(env)
    @test haskey(Tulip.TLP_DEFAULT_PARAM, k)
end

# check all parameters are initialized to default
# this also tests all the getters
for (k, v) in Tulip.TLP_DEFAULT_PARAM
    @test env[k] == v
    @test env[String(k)] == v
end

# test setters
for (k, v) in Tulip.TLP_DEFAULT_PARAM
    env[k] = v
    env[String(k)] = v
end

# Errors should be caught
try env["UNDEFINED"] = 0.0
catch err
    @test isa(err, ErrorException)
end

try env[:output_level] = -1
catch err
    @test isa(err, AssertionError)
end

try env[:barrier_iter_max] = -1
catch err
    @test isa(err, AssertionError)
end

try env[:time_limit] = -1
catch err
    @test isa(err, AssertionError)
end

try env[:barrier_tol_feas] = -0.01
catch err
    @test isa(err, AssertionError)
end
try env[:barrier_tol_feas] = 1.01
catch err
    @test isa(err, AssertionError)
end

try env[:barrier_tol_opt] = -0.01
catch err
    @test isa(err, AssertionError)
end
try env[:barrier_tol_opt] = 1.01
catch err
    @test isa(err, AssertionError)
end

try env[:barrier_tol_conv] = -0.01
catch err
    @test isa(err, AssertionError)
end
try env[:barrier_tol_conv] = 1.01
catch err
    @test isa(err, AssertionError)
end