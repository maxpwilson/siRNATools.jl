
function OTA_Excelfile(dfs, species, name)
    @assert name[end-3:end] == "xlsx"
    xlsxwriter = pyimport("xlsxwriter")
    workbook = xlsxwriter.Workbook(name)
    cellformat = workbook.add_format(Dict("font_name" => "Courier New"))
    red = workbook.add_format(Dict("color" => "red", "font_name" => "Courier New"))
    black = workbook.add_format(Dict("color" => "black", "font_name" => "Courier New"))
    for x in 1:length(dfs)
        worksheet = workbook.add_worksheet(species[x])
        cols = []

        push!(cols, Dict("header" => "Acc#"))
        push!(cols, Dict("header" => "Gene ID"))
        push!(cols, Dict("header" => "Gene Symbol"))
        push!(cols, Dict("header" => "Description"))
        push!(cols, Dict("header" => "Region"))
        push!(cols, Dict("header" => "MM# 2-18"))
        push!(cols, Dict("header" => "AS 1-19 (5'->3')"))
        push!(cols, Dict("header" => "RC of off-target site in mRNA (5'->3')"))
        push!(cols, Dict("header" => "Mismatch Positions in AS"))
        push!(cols, Dict("header" => "Transcript Location"))
        widths = [18.86,11.86,17.14,60.14,17.57,13.14,28.14,53.14,34.13,27.86]
        for y in 1:size(dfs[x])[1]
            for z in 1:size(dfs[x])[2]
                if z != 8
                    worksheet.write(y, z-1, "$(dfs[x][y, z])")
                else
                    printer = []
                    for j in dfs[x][y,z]
                        if lowercase(j) == j
                            push!(printer, red)
                            push!(printer, string(uppercase(j)))
                        else
                            push!(printer, black)
                            push!(printer, string(j))
                        end
                    end
                    worksheet.write_rich_string(y, z-1, printer...)
                end
            end
        end
        worksheet.add_table(0,0,size(dfs[x])[1], size(dfs[x])[2]-1, Dict("columns" => cols))
        for i in 1:size(dfs[x])[2]
            worksheet.set_column(i-1, i-1, widths[i], cellformat)
        end
    end
    workbook.close()
end

function Counts_Excelfile(dfs, species, name)
    @assert name[end-3:end] == "xlsx"
    @assert length(unique(size.(dfs))) == 1
    @assert unique(names.(dfs))[1] == ["ID", "Pattern", "Zero", "One", "Two", "Three", "Four", "Score"]
    xlsxwriter = pyimport("xlsxwriter")
    workbook = xlsxwriter.Workbook(name)
    header_format = workbook.add_format(Dict("bg_color" => "white", "font_name" => "Calibri", "font_size" => 11, "align" => "center", "right" => 1, "top" => 1, "left" => 1))

    header_format_left = workbook.add_format(Dict("bg_color" => "white", "font_name" => "Calibri", "font_size" => 11, "align" => "center", "right" => 1, "top" => 1, "left" => 2))

    header_format_right = workbook.add_format(Dict("bg_color" => "white", "font_name" => "Calibri", "font_size" => 11, "align" => "center", "right" => 2, "top" => 1, "left" => 1))


    header_format_bold = workbook.add_format(Dict("bg_color" => "white", "font_name" => "Calibri", "font_size" => 11, "bold" => true, "align" => "center", "border" => 2))

    header_format_green = workbook.add_format(Dict("bg_color" => "#CCFFCC", "font_name" => "Calibri", "font_size" => 11, "align" => "top", "bottom" => 2, "top" => 1, "left" => 1, "right" => 1))
    header_format_green.set_left_color("#D8DBD8")
    header_format_green.set_right_color("#D8DBD8")
    header_format_green.set_top_color("#D8DBD8")
    header_format_green.set_align("center")

    header_format_green_right = workbook.add_format(Dict("bg_color" => "#CCFFCC", "font_name" => "Calibri", "font_size" => 11, "align" => "top", "bottom" => 2, "top" => 1, "left" => 1, "right" => 1))
    header_format_green_right.set_left_color("#D8DBD8")
    header_format_green_right.set_align("center")
    header_format_green_right.set_top_color("#D8DBD8")
    header_format_green_right.set_right(2)

    header_format_green_rotated = workbook.add_format(Dict("bg_color" => "#CCFFCC", "font_name" => "Calibri", "font_size" => 11, "align" => "top", "rotation" => 90, "bottom" => 2, "left" => 1, "right" => 1, "top" => 1))
    header_format_green_rotated.set_right_color("#D8DBD8")
    header_format_green_rotated.set_top_color("#D8DBD8")
    header_format_green_rotated.set_align("center")

    header_format_green_rotated_left = workbook.add_format(Dict("bg_color" => "#CCFFCC", "font_name" => "Calibri", "font_size" => 11, "align" => "top", "rotation" => 90, "bottom" => 2, "left" => 1, "right" => 1, "top" => 1))
    header_format_green_rotated_left.set_right_color("#D8DBD8")
    header_format_green_rotated_left.set_top_color("#D8DBD8")
    header_format_green_rotated_left.set_align("center")
    header_format_green_rotated_left.set_left(2)

    default_format = workbook.add_format(Dict("bg_color" => "white", "font_name" => "Arial", "font_size" => 10, "align" => "center", "border_color" => "#D8DBD8", "border" => 1))

    ot_format = workbook.add_format(Dict("bg_color" => "#FC9CCC", "font_name" => "Arial", "font_size" => 10, "align" => "center", "border" => 1, "border_color" => "#FC9CCC"))

    score_format_low = workbook.add_format(Dict("bg_color" => "#FFCC99", "font_name" => "Arial", "font_size" => 10, "align" => "center", "border" => 1, "border_color" => "#FFCC99"))

    score_format_mid = workbook.add_format(Dict("bg_color" => "#FFFF99", "font_name" => "Arial", "font_size" => 10, "align" => "center", "border" => 1, "border_color" => "#FFFF99"))

    score_format_high = workbook.add_format(Dict("bg_color" => "#99CC00", "font_name" => "Arial", "font_size" => 10, "align" => "center", "border" => 1, "border_color" => "#99CC00"))

    worksheet = workbook.add_worksheet("Counts")
    for i in 1:(size(dfs[1])[1])
        worksheet.write(i+3, 0, "$(dfs[1][i, :ID])", default_format)
        worksheet.write(i+3, 1, "$(dfs[1][i, :Pattern])", default_format)
    end
    worksheet.set_column(0, 0, 6.71)
    worksheet.set_column(1, 1, 30.57)
    worksheet.autofilter(3, 2, size(dfs[1])[1] + 3, ((length(dfs)-1) * 6) + 7)
    worksheet.set_row(3, 61)
    for i in 1:(length(dfs)*6)
        if mod(i, 6) == 1
            worksheet.set_column(i+1, i+1, 4.57)
            if i ==1
                worksheet.write(3, i+1, "Score", header_format_green_rotated_left)
            else
                worksheet.write(3, i+1, "Score", header_format_green_rotated)
            end
        elseif mod(i, 6) == 2
            worksheet.set_column(i+1, i+1, 2.71)
            worksheet.write(3, i+1, 0, header_format_green)
        elseif mod(i, 6) == 3
            worksheet.set_column(i+1, i+1, 2.71)
            worksheet.write(3, i+1, 1, header_format_green)
        elseif mod(i, 6) == 4
            worksheet.set_column(i+1, i+1, 3.57)
            worksheet.write(3, i+1, 2, header_format_green)
        elseif mod(i, 6) == 5
            worksheet.set_column(i+1, i+1, 5)
            worksheet.write(3, i+1, 3, header_format_green)
        elseif mod(i, 6) == 0
            worksheet.set_column(i+1, i+1, 5.57)
            if i == length(dfs) * 6
                worksheet.write(3, i+1, 4, header_format_green_right)
            else
                worksheet.write(3, i+1, 4, header_format_green)
            end
        end
    end
    worksheet.merge_range(0, 2, 0, ((length(dfs) - 1) * 6) + 7, "siRNA AS specificity (detailed) (GENES) RefSeq rel. $(VERSION)", header_format_bold)
    for i in 1:length(species)
        if i == 1
            worksheet.merge_range(1, ((i - 1) * 6) + 2, 1, ((i - 1) * 6) + 7, "$(species[i])", header_format_left)
            worksheet.merge_range(2, ((i - 1) * 6) + 2, 2, ((i - 1) * 6) + 7, "OT freq by # of mm", header_format_left)
        elseif i == length(species)
            worksheet.merge_range(1, ((i - 1) * 6) + 2, 1, ((i - 1) * 6) + 7, "$(species[i])", header_format_right)
            worksheet.merge_range(2, ((i - 1) * 6) + 2, 2, ((i - 1) * 6) + 7, "OT freq by # of mm", header_format_right)
        else
            worksheet.merge_range(1, ((i - 1) * 6) + 2, 1, ((i - 1) * 6) + 7, "$(species[i])", header_format)
            worksheet.merge_range(2, ((i - 1) * 6) + 2, 2, ((i - 1) * 6) + 7, "OT freq by # of mm", header_format)
        end
    end
    for x in 1:length(dfs)
        df = dfs[x]
        for i in 1:(size(df))[1]
            if df[i, :Score] > 3
                worksheet.write_number(i+3, (6 * (x - 1)) + 2, (df[i, :Score]), score_format_high)
            elseif df[i, :Score] == 3
                worksheet.write_number(i+3, (6 * (x - 1)) + 2, (df[i, :Score]), score_format_mid)
            elseif df[i, :Score] > 2
                worksheet.write_number(i+3, (6 * (x - 1)) + 2, (df[i, :Score]), score_format_low)
            else
                worksheet.write_number(i+3, (6 * (x - 1)) + 2, (df[i, :Score]), default_format)
            end
            if df[i, :Zero] == 0
                worksheet.write_number(i+3, (6 * (x - 1)) + 3, (df[i, :Zero]), default_format)
            else
                worksheet.write_number(i+3, (6 * (x - 1)) + 3, (df[i, :Zero]), ot_format)
            end
            if df[i, :One] == 0
                worksheet.write_number(i+3, (6 * (x - 1)) + 4, (df[i, :One]), default_format)
            else
                worksheet.write_number(i+3, (6 * (x - 1)) + 4, (df[i, :One]), ot_format)
            end
            if df[i, :Two] == 0
                worksheet.write_number(i+3, (6 * (x - 1)) + 5, (df[i, :Two]), default_format)
            else
                worksheet.write_number(i+3, (6 * (x - 1)) + 5, (df[i, :Two]), ot_format)
            end
            if df[i, :Three] == 0
                worksheet.write_number(i+3, (6 * (x - 1)) + 6, (df[i, :Three]), default_format)
            else
                worksheet.write_number(i+3, (6 * (x - 1)) + 6, (df[i, :Three]), ot_format)
            end
            if df[i, :Four] == 0
                worksheet.write_number(i+3, (6 * (x - 1)) + 7, (df[i, :Four]), default_format)
            else
                worksheet.write_number(i+3, (6 * (x - 1)) + 7, (df[i, :Four]), ot_format)
            end

        end
    end
    workbook.close()
end

function GeneralReport_Excelfile(report, transcripts, tSpecies, mSpecies, transcript, name, species)
        @assert name[end-3:end] == "xlsx"
        xlsxwriter = pyimport("xlsxwriter")
        workbook = xlsxwriter.Workbook("$PATH/Output_Files/$name")
        report_sheet = workbook.add_worksheet("$species $transcript")
        transcript_sheet = workbook.add_worksheet("Transcripts")

        ts_columns = []
        push!(ts_columns, Dict("header" => "Species"))
        push!(ts_columns, Dict("header" => "Transcript"))
        push!(ts_columns, Dict("header" => "BasePairs"))
        push!(ts_columns, Dict("header" => "Ensembl Match"))
        push!(ts_columns, Dict("header" => "Include?"))
        transcripts[!, :Include] .= "Yes"
        data = [Vector(transcripts[i, [:Species, :Transcript, :BasePairs, :EnsemblMatch, :Include]]) for i in 1:size(transcripts)[1]]
        transcript_sheet.add_table(1, 1, size(transcripts)[1]+1, length(ts_columns), Dict("data" => data, "columns" => ts_columns, "style" => "Table Style Light 1"))

        report_columns = []
        for name in names(report)
            push!(report_columns, Dict("header" => name))
        end
        data = [Vector(report[i, :]) for i in 1:size(report)[1]]

        report_sheet.add_table(1,1, size(report)[1]+1, size(report)[2], Dict("data" => data, "columns" => report_columns, "style" => "Table Style Light 1"))
        workbook.close()
end

function miRNA_ExcelFile(report, name)
    @assert name[end-3:end] == "xlsx"
    for col in names(report)
        report[ismissing.(report[col]),col] = ""
    end
    xlsxwriter = pyimport("xlsxwriter")
    workbook = xlsxwriter.Workbook("$PATH/Output_Files/$name")
    worksheet = workbook.add_worksheet("miRNA Matches")
    columns = []
    push!(columns, Dict("header" => "Family ID"))
    push!(columns, Dict("header" => "ID"))
    push!(columns, Dict("header" => "Organism"))
    push!(columns, Dict("header" => "Common Name"))
    push!(columns, Dict("header" => "Family"))
    push!(columns, Dict("header" => "Mature Sequence"))
    data = [Vector(report[i, :]) for i in 1:size(report)[1]]
    worksheet.add_table(0,0, size(report)[1], size(report)[2] - 1, Dict("data" => data, "columns" => columns, "style" => "Table Style Light 1"))
    widths = [21.71, 15, 28.57, 37.71, 15.86, 43.86]
    for i in 1:size(report)[2]
        worksheet.set_column(i-1, i-1, widths[i])
    end
    workbook.close()
end

function SNP_ExcelFile(report, name)
    @assert name[end-3:end] == "xlsx"
    for col in names(report)[3:end]
        report[ismissing.(report[!, col]),col] .= ""
    end
    xlsxwriter = pyimport("xlsxwriter")
    workbook = xlsxwriter.Workbook("$PATH/Output_Files/$name")
    worksheet = workbook.add_worksheet("SNP Matches")
    columns = []
    push!(columns, Dict("header" => "Sequence"))
    push!(columns, Dict("header" => "Pos"))
    push!(columns, Dict("header" => "ID"))
    push!(columns, Dict("header" => "SNP Pos"))
    push!(columns, Dict("header" => "MAF"))
    push!(columns, Dict("header" => "ID MAF>0"))
    push!(columns, Dict("header" => "SNP Pos MAF>0"))
    push!(columns, Dict("header" => "MAF>0"))
    push!(columns, Dict("header" => "IDs MAF>0.01"))
    push!(columns, Dict("header" => "SNP Pos MAF>0.01"))
    push!(columns, Dict("header" => "MAF>0.01"))
    data = [string.(Vector(report[i, :])) for i in 1:size(report)[1]]
    worksheet.add_table(0,0, size(report)[1], size(report)[2] - 1, Dict("data" => data, "columns" => columns, "style" => "Table Style Light 1"))
#    widths = [0,0,0,0,0,0,0,0,0,0,0]
#    for i in 1:size(report)[2]
#        worksheet.set_column(i-1, i-1, widths[i])
#    end
    workbook.close()
end
