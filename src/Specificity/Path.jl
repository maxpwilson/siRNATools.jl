
src = @__DIR__

if !("Specificity_Path.txt" in readdir(src))
    touch("$(src)/Specificity_Path.txt")
end

if !isdir(readline(open("$(src)/Specificity_Path.txt")))
    println("Current specified Path does not exist.  Update using siRNATools.Specificity.Update_Path(PATH). Then re-building siRNATools ")
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

Load_Path()
