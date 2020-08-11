
src = @__DIR__

if !("Specificity_Path.txt" in readdir(src))
    touch("$(src)/Specificity_Path.txt")
end

if !isdir(readline(open("$(src)/Specificity_Path.txt")))
    println("Current specified Path does not exist.  Update using siRNATools.Informatics.Update_Path(PATH).")
end

"""
    Update_Path(::String)

Updates default value for constant PATH.  All created or downloaded files are put into PATH directory or sub-directory created by that function.
Persistent on each computer, but needs to be updated after install before calculating specificity.
"""
function Update_Path(PATH::String)
    write("$(src)/Specificity_Path.txt", PATH);
end

function Load_Path()
    global PATH = readline(open("$(src)/Specificity_Path.txt"))
end
function Load_Version()
    try
        if !("RefSeq_version.txt" in readdir(PATH))
            touch("$(PATH)/RefSeq_version.txt")
        end
        global VERSION = readline(open("$(PATH)/RefSeq_version.txt"))
    catch
    end
end

Load_Path()
Load_Version()
