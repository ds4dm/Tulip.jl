env = Tulip.TulipEnv()

env_ = copy(env)
@test env_.verbose.val == env.verbose.val

# Test getters
v = env.verbose.val
@test env[:verbose] == v
@test env["verbose"] == v

# Test setters
env[:verbose] = 1
@test env.verbose.val == 1
env["verbose"] = 0
@test env.verbose.val == 0

Tulip.set_param!(env, verbose=1, time_limit=10.0)
@test env.verbose.val == 1
@test env.time_limit.val == 10.0

# Test reset
Tulip.reset!(env)
for s in fieldnames(Tulip.TulipEnv)
    v_ = env[s]
    @test v_ == Core.getfield(env, s).def_val
end