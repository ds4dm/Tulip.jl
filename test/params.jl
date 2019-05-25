@testset "Params" begin

p = Tulip.RealParam(:param, 0, 0, 10)
p.val = 1
@test Tulip.get_param_name(p) == :param
@test Tulip.get_param_type(p) == Int64
@test Tulip.get_param_value(p) == 1

Tulip.set_param_default!(p)
@test p.val == 0

Tulip.set_param_value!(p, 5)
@test p.val == 5

@test Tulip.test_param_value(p, 5)
@test !Tulip.test_param_value(p, 100)

end  # testset