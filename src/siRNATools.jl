module siRNATools
export SSF_DF

using PyCall, DataFrames, Dates, ProgressMeter

#ENV["PYTHON"] = "C:\\Users\\mwilson\\AppData\\Local\\Continuum\\anaconda3\\python.exe"

pb = pyimport("pyxlsb")

searchdir(path, key) = filter(x->occursin(key,x), readdir(path))

function int_to_date(x::Number) :: Date
    Date(1900) + Day(x)
end

function clean_value(x, type::DataType, f=eval::Function, g=eval::Function) :: type
    (x === nothing) && (x = 0)
    if !(x isa type)
        if x isa String
            try 
                x = parse(Float64, x)
            catch
                x = 0
            end
        elseif type == String
            try
                x = string(x)
            catch
                x = ""
            end
        end
    end
    try
        x = f(x)
    catch
    end
    try
        x = g(x)
    catch
    end
    (type == Int) ? (return floor(x)) : return x
end

function calc_batch_number(notebook::Int, page::Int, type::Int, position::Int) :: String
    (type == 12) ? "$notebook-$(page)M_$position" : "$notebook-$(page)$(string(Char(((position-1)%8)+65)))_$(Int(floor((position-1)/8)+1))"
end

function SSF_DF(NB::Int64) :: DataFrame
    path = "R:\\Chemistry\\siRNA\\Single Strands\\$NB\\"
    df = DataFrame(
    StrandID = String[],
    NoteBook = Int64[],
    Page = Int64[],
    SheetType = Int64[],
    Position = Int64[],
    Scale = Float64[],
    Needed = Float64[],
    NeedBy = Date[],
    Target = String[],
    Sequence = String[],
    )
    for file in searchdir(path, ".xlsb")
        sht = try
            pb.open_workbook(path * file).get_sheet("Sample Set Form") 
        catch 
            println(file, " could not be opened")
            continue
        end
        PG, ST = 0, 0
        for (index, row) in enumerate(sht.rows(sparse=true))
            if index == 3
                PG = row[3][3]
                (row[4][3] == "M") ? ST = 12 : ST = 48
            end
            if (row[3][3] != nothing) && (row[3][3] != "") && (index in 7:54)
                SCL = clean_value(row[4][3], Float64)
                NDD = clean_value(row[5][3], Float64)
                NDB = clean_value(row[8][3], Date, int_to_date)
                TRG = clean_value(row[14][3], String)
                SEQ = clean_value(row[15][3], String)
                push!(df, [row[3][3], NB, PG, ST, (index-6), SCL, NDD, NDB, TRG, SEQ])
            end
        end
    end
    df.BatchNumber = calc_batch_number.(df.NoteBook, df.Page, df.SheetType, df.Position)
    df
end

function CF_DF(NB::Int64) :: DataFrame
    path = "R:\\Chemistry\\siRNA\\Single Strands\\$NB\\"
    df = DataFrame(
    StrandID = String[],
    BatchNumber = String[],
    CleavageStartDate = Date[],
    Comments = String[],
    CrudeMS = Float64[],
    CrudePurity = Float64[],
    MW = Float64[],
    CrudeYield = Float64[],
    )
    for file in searchdir(path, ".xlsb")
        sht = try
            pb.open_workbook(path * file).get_sheet("Cleavage Form") 
        catch 
            println(file, " could not be opened")
            continue
        end
        PG, ST = 0, 0
        for (index, row) in enumerate(sht.rows(sparse=true))
            if (row[2][3] != nothing) && (row[2][3] != "") && (row[2][3][end] == 'S') &&(index in 7:54)
                SID = clean_value(row[2][3], String)
                BN = replace(clean_value(row[3][3], String), " crude" => "")
                CSD = clean_value(row[7][3], Date, int_to_date)
                COM = clean_value(row[9][3], String)
                CMS = clean_value(row[13][3], Float64)
                CP = clean_value(row[14][3], Float64)
                MW = clean_value(row[18][3], Float64)
                CY = clean_value(row[20][3], Float64)
                push!(df, [SID, BN, CSD, COM, CMS, CP, MW, CY])
            end
        end
    end
    df
end

end # module
