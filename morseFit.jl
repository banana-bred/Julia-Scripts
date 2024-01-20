##################################################################
# fit a potential energy curve to a Morse potential
##################################################################

using DelimitedFiles: writedlm, readdlm
using LsqFit

global PROGRAM = basename(@__FILE__)

##################################################################
######################### FUNCTIONS ##############################
##################################################################

# IDEAS: FIT ONLY PART OF IT (HEAD, TAIL)
# IDEAS: FIT BOTH (HEAD, TAIL)

function usage()
  println()
  println(PROGRAM,": Fit a potential energy curve to a Lennard-Jones potential")
  println()
  println("  usage: ", PROGRAM, "[FILE] [OPTIONS]")
  println()
  return true
end # usage

function die(message::String)
    println("@@@@@@@@@@@@@@@@@@@@@@@@@")
    println("ERROR: ",message)
    println("@@@@@@@@@@@@@@@@@@@@@@@@@")
    usage()
    exit()
end # die

function Morse(r, c::AbstractArray)
    a  = c[1]
    r0 = c[2] # equilibrium distance
    D  = c[3] # depth
    diss = c[4] # dissociation limit
    0 in r ? die("At least one of the supplied internuclear distances is 0") :
        @. return D * (exp(-2 * a * (r - r0)) - 2 * exp(-a * (r - r0))) + diss
end # Morse

function main(
        infile::String,
        outfile::String,
        tailOnly::Bool,
        rCutoff::Float64,
        rMin::Float64,
        rMax::Float64,
        rStep::Float64,
        printParams::Bool
    )

    # -- read data into x and y arrays
    data = readdlm(infile)
    xdata = data[:,1]
    ydata = data[:,2]
    data = 0 # <-- clear data (it shouldn't be that big anyway)

    if tailOnly
         y = ydata[xdata .> rCutoff]
         x = xdata[xdata .> rCutoff]
     else
         x = xdata
         y = ydata
    end # if

    global rmin = x[argmin(y)]

    a = 0.18    # -- potential width
    r0 = rmin   # -- location of minimum
    D = 0.15    # -- well depth
    diss = -75  # -- dissociation limit
    c = [a, r0, D, diss]

    fit = curve_fit(Morse, x, y, c)

    if tailOnly
        rMin = last(x)
    end #if

    xfit = [r for r in rMin:rStep:rMax]
    yfit = Morse(xfit, fit.param)

    if tailOnly
        xfit = [xdata; xfit]
        yfit = [ydata; yfit]
    end #if

    writedlm(outfile, [xfit yfit])

    if printParams
        println("a, r0, D, dissociationLimit  = ",  fit.param)
    end # if

end

infile   = ARGS[1]
outfile  = ARGS[2]
tailOnly = parse(Bool, ARGS[3]) # -- T: only fit the tail of the potential
                                #    F: fit the whole thing
rCutoff = parse(Float64, ARGS[4]) # -- cutoff distance for fitting the tail. Values smaller than this are ignored
rMin    = parse(Float64, ARGS[5]) # -- smallest distance for theARGS[5] fit. Ignored if tailOnly is true
rMax    = parse(Float64, ARGS[6]) # -- largest  distance for theARGS[6] fit
rStep   = parse(Float64, ARGS[7]) # -- linear step size for the ARGS[7]fit
printParams = parse(Bool, ARGS[8])

main(infile, outfile, tailOnly, rCutoff, rMin, rMax, rStep, printParams)