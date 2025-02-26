module LagrangeInterp

using LinearAlgebra


function lag(x,j,Xdata)
    n = length(Xdata)
    l = 1
    for m = 1:n
        if m ≠ j
            l *= (x - Xdata[m])/(Xdata[j] - Xdata[m])
        end
    end
    return l
end

function interp(Xdata,Ydata)

    n = length(Xdata)

    return function (x)
        return sum( lag(x,j,Xdata)*Ydata[j] for j in 1:n  )
    end

end


end