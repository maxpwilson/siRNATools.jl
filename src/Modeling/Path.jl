
src = @__DIR__

if !("Modeling_Path.txt" in readdir(src))
    touch("$(src)/Modeling_Path.txt")
end

if !isdir(readline(open("$(src)/Modeling_Path.txt")))
    println("Current specified Path does not exist.  Update using siRNATools.Modeling.Update_Path(PATH).")
end

"""
    Update_Path(::String)

Updates default value for constant PATH.  All created or downloaded files are put into PATH directory or sub-directory created by that function.
Persistent on each computer, but needs to be updated after install before calculating specificity.
"""
function Update_Path(PATH::String)
    write("$(src)/Modeling_Path.txt", PATH);
end

function Load_Path()
    global PATH = readline(open("$(src)/Modeling_Path.txt"))
end

Load_Path()
