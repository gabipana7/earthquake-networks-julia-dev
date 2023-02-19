using PyCall

### Julia wrapper on nemomap Python 
function motifs_discovery(inputName,queryName)

    # Define workind directory and path to nemopmap
    scriptdir = @__DIR__ 
    nemomapdir = scriptdir * "/src/nemomap"
    pushfirst!(PyVector(pyimport("sys")."path"), nemomapdir)

    # Include motifsdiscovery python code
    @pyinclude("./src/nemomap/motifsdiscovery.py")
    
    
    stats = py"getMotif"(inputName,queryName)

    # This code changes the output of original nemomap code to a more readable csv
    py"""

    import ast, os, csv

    fileMotif=open("./output.txt")
    linesMotif = fileMotif.readlines()

    fileMotif.close()
    os.remove('output.txt')

    # Properly evaluate the Lines to get the Lists
    motifNodes=[]
    for item in linesMotif:
        motifNodes.append(ast.literal_eval(item))

    newmotifs=[]
    for motifs in motifNodes:
        res = [eval(i) for i in motifs]
        newmotifs.append(res)

    with open("output.csv", "w", newline='') as f:
        wr = csv.writer(f)
        wr.writerows(newmotifs)   
    """

    return stats
end