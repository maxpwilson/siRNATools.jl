try
catch
end

using DataFrames, Dates, ProgressMeter, Query, PyCall

searchdir(path, key) = filter(x->occursin(key,x), readdir(path))

function int_to_date(x::Number) :: Date
    d = Date(1900) + Day(floor(x)) 
    (d > Date(1900)) && (d = d - Day(2))
    d
end

function GetPyxlsb() :: PyObject
    try
        return pyimport("pyxlsb")
    catch
        try 
            run(`python -m pip install pyxlsb --quiet`);
            return pyimport("pyxlsb")
        catch
            error("Couldn't get pyxlsb")
        end
    end
end

function clean_value(x::Any, type::DataType, f::Function...) :: type
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
    for g in f
        try
            x = g(x)
        catch
        end
    end
    (type == Int) ? (return floor(x)) : return x
end
function clean_value(x::Any, type::DataType) :: type
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
    (type == Int) ? (return floor(x)) : return x
end

function calc_batch_number(notebook::Int, page::Int, type::Int, position::Int) :: String
    (type == 12) ? "$notebook-$(page)M_$position" : "$notebook-$(page)$(string(Char(((position-1)%8)+65)))_$(Int(floor((position-1)/8)+1))"
end

function SSF_DF(NB::Int64) :: DataFrame
    pb = GetPyxlsb()
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
    pb = GetPyxlsb()
    path = "R:\\Chemistry\\siRNA\\Single Strands\\$NB\\"
    df = DataFrame(
    StrandID = String[],
    BatchNumber = String[],
    CleavageStartDate = Date[],
    SynthesisComments = String[],
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

function DF_DF(NB::Int64) :: DataFrame
    pb = GetPyxlsb()
    path = "R:\\Chemistry\\siRNA\\Single Strands\\$NB\\"
    df = DataFrame(
        BatchNumber = String[],
        DeprotectionComments = String[]
    )
    for file in searchdir(path, ".xlsb")
        sht = try
            pb.open_workbook(path * file).get_sheet("Deprotection Form") 
        catch 
            println(file, " could not be opened")
            continue
        end
        for (index, row) in enumerate(sht.rows(sparse=true))
            if (row[2][3] != nothing) && (row[2][3] != "") && (row[2][3][end] == 'S') &&(index in 7:54)
                BN = replace(clean_value(row[3][3], String), " postTBDMS" => "")
                COM = clean_value(row[16][3], String)
                push!(df, [BN, COM])
            end
        end
    end
    df
end

function PF_DF(NB::Int64) :: DataFrame
    pb = GetPyxlsb()
    path = "R:\\Chemistry\\siRNA\\Single Strands\\$NB\\"
    df = DataFrame(
        StrandID = String[],
        BatchNumber = String[],
        Instrument = String[],
        Column = String[],
        Method = String[],
        PurStartDate = Date[],
        PrepDate = Date[],
        RetentionTime = Float64[],
        FractionsKept = String[],
        InjectionVolume = Float64[],
        PurificationComments = String[],
    )
    for file in searchdir(path, ".xlsb")
        sht = try
            pb.open_workbook(path * file).get_sheet("Purification Form") 
        catch 
            println(file, " could not be opened")
            continue
        end
        for (index, row) in enumerate(sht.rows(sparse=true))
            if (row[2][3] != nothing) && (row[2][3] != "") && (row[2][3][end] == 'S') &&(index in 7:54)
                SID = clean_value(row[2][3], String)
                BN = replace(clean_value(row[3][3], String), " purified" => "")
                INS = clean_value(row[5][3], String)
                COL = clean_value(row[6][3], String)
                MET = clean_value(row[7][3], String)
                PSD = clean_value(row[8][3], Date, int_to_date)
                PRD = clean_value(clean_value(row[9][3], Int), Date, int_to_date)
                RT = clean_value(row[10][3], Float64)
                FK = clean_value(row[11][3], String)
                IV = clean_value(row[12][3], Float64)
                COM = clean_value(row[13][3], String)
                push!(df, [SID, BN, INS, COL, MET, PSD, PRD, RT, FK, IV, COM])
            end
        end
    end
    df
end

function sht_df(NB::Int64) :: DataFrame
    df_pf = PF_DF(NB)
    df_ssf = SSF_DF(NB)
    df_cf = CF_DF(NB)
    df_df = DF_DF(NB)
    df = @from pf in df_pf begin
        @join ssf in df_ssf on pf.BatchNumber equals ssf.BatchNumber
        @join cf in df_cf on pf.BatchNumber equals cf.BatchNumber
        @join df in df_df on pf.BatchNumber equals df.BatchNumber
        @select {ssf.StrandID, pf.BatchNumber, ssf.NoteBook, ssf.Page, ssf.SheetType, ssf.Position, ssf.Scale, ssf.Needed, ssf.NeedBy, ssf.Sequence, ssf.Target,
            cf.CleavageStartDate, cf.SynthesisComments, cf.CrudeMS, cf.MW, cf.CrudePurity, cf.CrudeYield, df.DeprotectionComments,
            pf.Instrument, pf.Column, pf.Method, pf.PurificationComments, pf.PurStartDate, pf.PrepDate, pf.RetentionTime, pf.FractionsKept, pf.InjectionVolume}
        @collect DataFrame
    end
    df
end

function shts_df(NBs::Int64...) :: DataFrame
    df = DataFrame()
    for NB in NBs
        df = vcat(df, sht_df(NB))
    end
    df
end

function q_df() :: DataFrame
    pb = GetPyxlsb()
    file = "R:\\Chemistry\\siRNA\\Synthesis Q.xlsb"
    df = DataFrame(
        RequestedBy = String[],
        RequestedDate = Date[],
        EstimatedDate = Date[],
        InjectionDate = Date[],
        Target = String[],
        DuplexID = String[],
        SetDetail = String[],
        AntiSenseID = String[],
        ASSynthesisScale = String[],
        ASStatus = String[],
        SenseID = String[],
        SSSynthesisScale = String[],
        SSStatus = String[],
        Notes = String[],
        SpecialAmidites = String[],
        AmtNeeded = String[],
        DuplexStatus = String[],
        DateCompleted = Date[],
        AnnealedBy = String[],
        DuplexBatchID = String[],
        Concentration = String[],
        AmountPrepared = Float64[],
    )
    sht = try
        pb.open_workbook(file).get_sheet("DUPLEX QUEUE") 
    catch 
        println("Queue could not be opened")
        return df
    end
    for (index, row) in enumerate(sht.rows(sparse=true))
        if (row[6][3] != nothing) && (row[6][3] != "") 
            REQBY = clean_value(row[1][3], String)
            REQD = clean_value(row[2][3], Date, int_to_date)
            ESTD = clean_value(row[3][3], Date, int_to_date)
            INJD = clean_value(row[4][3], Date, int_to_date)
            TRGT = clean_value(row[5][3], String)
            DUPID = clean_value(row[6][3], String)
            SETDT = clean_value(row[7][3], String)
            ASID = clean_value(row[8][3], String)
            ASSS = clean_value(row[9][3], String)
            ASS = clean_value(row[10][3], String)
            SSID = clean_value(row[12][3], String)
            SSSS = clean_value(row[13][3], String)
            SSS = clean_value(row[14][3], String)
            NTS = clean_value(row[15][3], String)
            SPAM = clean_value(row[16][3], String)
            AMTN = clean_value(row[17][3], String)
            DUPS = clean_value(row[18][3], String)
            COMD = clean_value(row[19][3], Date, int_to_date)
            ANB = clean_value(row[20][3], String)
            DUBID = clean_value(row[21][3], String)
            CONC = clean_value(row[23][3], String)
            AMTP = clean_value(row[24][3], Float64)
            push!(df, [REQBY, REQD, ESTD, INJD, TRGT, DUPID, SETDT, ASID, ASSS, ASS, SSID, SSSS, SSS, NTS, SPAM, AMTN, DUPS, COMD, ANB, DUBID, CONC, AMTP])
        end
    end
    df
end

function InQueue() :: DataFrame
    df_q = q_df()
    df_q[(df_q.AnnealedBy .== "") .& (df_q.DateCompleted .== Date(1900, 1,1)) .& (df_q.RequestedDate .>= (today() - Month(7))) .& (df_q.DuplexStatus .!= "Canceled") .& (df_q.DuplexStatus .!= "Archived") .& (df_q.DuplexStatus .!= "Madison") .& (df_q.DuplexStatus .!= "Handed Off"), :]
end

function InProgress(status::String) :: DataFrame
    df_q = InQueue()
    df_q[(df_q.ASStatus .== status) .| (df_q.SSStatus .== status), :]
end