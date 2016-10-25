module Helper
export part

"""
part{T}(arr::Array{T,1}, f)

takes in an array (arr) and boolean function (f) and
returns a tuple of arrays; the left containing the
elements of (a) that yielded f(a) == false, and the
right containint the remaining elements.
"""
function part{T}(arr::Array{T,1}, f)
  l = Array{T,1}()
  r = Array{T,1}()

  for a in arr
    if f(a)
      push!(r,a)
    else
      push!(l,a)
    end
  end

  (l,r)
end

end
