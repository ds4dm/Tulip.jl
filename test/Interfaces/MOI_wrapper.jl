#  Copyright 2018-2019: Mathieu Tanneau
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

using Test

import MathOptInterface as MOI
import Tulip

@testset "Direct optimizer" begin
    model = Tulip.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    MOI.Test.runtests(
        model,
        MOI.Test.Config(
            Float64;
            atol = 1e-6,
            rtol = 1e-6,
            exclude = Any[MOI.ConstraintBasisStatus, MOI.VariableBasisStatus],
        ),
    )
end

@testset "MOI Bridged" begin
    model = MOI.Bridges.full_bridge_optimizer(Tulip.Optimizer(), Float64)
    MOI.set(model, MOI.Silent(), true)
    MOI.Test.runtests(
        model,
        MOI.Test.Config(
            Float64;
            atol = 1e-6,
            rtol = 1e-6,
            exclude = Any[MOI.ConstraintBasisStatus, MOI.VariableBasisStatus],
        ),
        exclude=[
            r"^test_conic_NormInfinityCone_INFEASIBLE$",
            r"^test_conic_NormOneCone_INFEASIBLE$",
        ],
    )
end

# Run the MOI tests with HSD and MPC algorithms
@testset "MOI Linear tests - $ipm" for ipm in [Tulip.HSD, Tulip.MPC]
    model = Tulip.Optimizer()
    model.inner.params.IPM.Factory = Tulip.Factory(ipm)
    MOI.set(model, MOI.Silent(), true)
    MOI.Test.runtests(
        model,
        MOI.Test.Config(
            Float64;
            atol = 1e-6,
            rtol = 1e-6,
            exclude = Any[MOI.ConstraintBasisStatus, MOI.VariableBasisStatus],
        ),
        include=["linear"],
    )
end

@testset "Cached optimizer" begin
    inner = MOI.Utilities.CachingOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        Tulip.Optimizer(),
    )
    model = MOI.Bridges.full_bridge_optimizer(inner, Float64)
    MOI.set(model, MOI.Silent(), true)
    MOI.Test.runtests(
        model,
        MOI.Test.Config(
            Float64;
            atol = 1e-6,
            rtol = 1e-6,
            exclude = Any[MOI.ConstraintBasisStatus, MOI.VariableBasisStatus],
        ),
        exclude=[
            r"^test_conic_NormInfinityCone_INFEASIBLE$",
            r"^test_conic_NormOneCone_INFEASIBLE$",
        ],
    )
end

@testset "test_attribute_TimeLimitSec" begin
    model = Tulip.Optimizer()
    @test MOI.supports(model, MOI.TimeLimitSec())
    @test MOI.get(model, MOI.TimeLimitSec()) === nothing
    MOI.set(model, MOI.TimeLimitSec(), 0.0)
    @test MOI.get(model, MOI.TimeLimitSec()) == 0.0
    MOI.set(model, MOI.TimeLimitSec(), nothing)
    @test MOI.get(model, MOI.TimeLimitSec()) === nothing
    MOI.set(model, MOI.TimeLimitSec(), 1.0)
    @test MOI.get(model, MOI.TimeLimitSec()) == 1.0
end
