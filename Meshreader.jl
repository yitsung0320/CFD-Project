
# Meshreader2D
# file format is # of grids in i direction; # of

# read from txt input = filename string
using DelimitedFiles
function readfiles(str::String)
  coordinate = readdlm(str,Float64)
  ni = convert(Int64, coordinate[1,1])
  nj = convert(Int64, coordinate[1,2])

  return ni , nj , coordinate[2:end,:]
end

using LinearAlgebra

# mutable makes the obj able to be change
mutable struct cell
   point1::Array{Float64,1}
   point2::Array{Float64,1}
   point3::Array{Float64,1}
   point4::Array{Float64,1}
   cellcenter::Array{Float64,1}
   area::Float64

   vector12::Array{Float64,1}
   vector23::Array{Float64,1}
   vector34::Array{Float64,1}
   vector41::Array{Float64,1}

   wallmidpt12::Array{Float64,1}
   wallmidpt23::Array{Float64,1}
   wallmidpt34::Array{Float64,1}
   wallmidpt41::Array{Float64,1}

   length12:: Float64
   length23:: Float64
   length34:: Float64
   length41:: Float64

   out_norm12::Array{Float64,1}
   out_norm23::Array{Float64,1}
   out_norm34::Array{Float64,1}
   out_norm41::Array{Float64,1}
   # object construction function (can only initialize))
   function cell(p1::T,p2::T,p3::T,p4::T) where T<:Array{Float64,1}
    point1 = p1
    point2 = p2
    point3 = p3
    point4 = p4
    new(p1,p2,p3,p4)
   end
end

#  cell property calculat the cellcenter, outward_norm, area, srfacemidpt,
#  surface length
function cellproperty(coord,ni::Int64,nj::Int64)

 meshcell = Array{cell,2}(undef,nj-1,ni-1)
 w = Array{Float64,1}(undef,3)   # angular vector for judge the mesh roation

  for i = 1: ni-1
  for j = 1 : nj-1
   meshcell[j,i] = cell(coord[i + ni*(j-1),:],coord[i + 1+ni*(j-1),:],
                        coord[i + ni*j+1,:],coord[i+ni*j,:])
   meshcell[j,i].cellcenter = (coord[i + ni*(j-1),:] + coord[i + 1+ni*(j-1),:]
                              +coord[i + ni*j+1,:] + coord[i+ni*j,:])/4

   meshcell[j,i].vector12 = coord[i+1+ni*(j-1),:] - coord[i+ni*(j-1),:]
   meshcell[j,i].vector23 = coord[i+1+ni*j,:] - coord[i+1+ni*(j-1),:]
   meshcell[j,i].vector34 = coord[i+ni*j,:] - coord[i+1+ni*j,:]
   meshcell[j,i].vector41 = coord[i+ni*(j-1),:] - coord[i+ni*j,:]

   meshcell[j,i].length12 = norm(meshcell[j,i].vector12)
   meshcell[j,i].length23 = norm(meshcell[j,i].vector23)
   meshcell[j,i].length34 = norm(meshcell[j,i].vector34)
   meshcell[j,i].length41 = norm(meshcell[j,i].vector41)

   meshcell[j,i].wallmidpt12 = (coord[i+1+ni*(j-1),:] + coord[i+ni*(j-1),:])/2
   meshcell[j,i].wallmidpt23 = (coord[i+1+ni*j,:] + coord[i+1+ni*(j-1),:])/2
   meshcell[j,i].wallmidpt34 = (coord[i+ni*j,:] + coord[i+1+ni*j,:])/2
   meshcell[j,i].wallmidpt41 = (coord[i+ni*(j-1),:] + coord[i+ni*j,:])/2

   meshcell[j,i].area = (norm(cross([meshcell[j,i].vector12;0],
                                    [meshcell[j,i].vector23;0])/2)+
                          norm(cross([meshcell[j,i].vector34;0],
                                    [meshcell[j,i].vector41;0])/2))

  # Determination of the curl (for deciding normal)
   w = cross([meshcell[j,i].vector12;0],[meshcell[j,i].vector23;0])
   d = w[3]/abs(w[3]) #d = 1 ccw, d = -1 cw

   meshcell[j,i].out_norm12 = [d*meshcell[j,i].vector12[2],
                        -d*meshcell[j,i].vector12[1]]/meshcell[j,i].length12
   meshcell[j,i].out_norm23 = [d*meshcell[j,i].vector23[2],
                        -d*meshcell[j,i].vector23[1]]/meshcell[j,i].length23
   meshcell[j,i].out_norm34 = [d*meshcell[j,i].vector34[2],
                        -d*meshcell[j,i].vector34[1]]/meshcell[j,i].length34
   meshcell[j,i].out_norm41 = [d*meshcell[j,i].vector41[2],
                        -d*meshcell[j,i].vector41[1]]/meshcell[j,i].length41

   end
  end
 return meshcell
end  #end function
