using siRNATools
using Test, Dates, PyCall

@testset "Clean Values" begin
    @test siRNATools.clean_value("10", Int) == 10
    @test siRNATools.clean_value("10", Float64) == 10.0
    @test siRNATools.clean_value(10, String) == "10"
    @test siRNATools.clean_value(10.6, Int) == 10
    @test siRNATools.clean_value(43000, Date, siRNATools.int_to_date) isa Date
end

@testset "Pyxlsb Loaded" begin
    @test occursin("pyxlsb", siRNATools.GetPyxlsb().__path__[1])
end