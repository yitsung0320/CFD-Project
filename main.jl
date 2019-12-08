 # CFD projet main file

 # 2-D incompressible Navier Stokes Problem
 # Coordinate generate / mesh builder

include("Meshreader.jl")
using LinearAlgebra
#=
prompt = "Insert the value of the square domain L"
L = parse(Float64, input(prompt))
prompt = " Insert the up wall velocity U"
U = parse(Float64, input(prompt))
Lx,Ly = L,L
mu = parse(Float64, input(prompt))
=#
#=
Lx, Ly = 2*pi , 2*pi
nx, ny = 100,100
L = Lx
rho = 1.0
mu = 2*pi*1e-3
U = 1.0
nu = mu/rho
Re = rho*U*L/mu
=#
#
Lx, Ly = 1 , 1
nx, ny = 129,129
L = 1
rho = 1
U = 1
mu = 0.25*1e-2
nu = mu/rho
Re = U*L/nu
#
grid = zeros(Float64,nx*ny,2)
for j = 1 : ny
 for i = 1 : nx
   k = (j-1)*nx+i
   grid[k,1] = Lx/(nx-1)*(i-1)
   grid[k,2] = Ly/(ny-1)*(j-1)
 end
end

meshcell = cellproperty(grid,nx,ny)

dist_i = ones(Float64,ny-1,nx)
dist_j = ones(Float64,ny,nx-1)
sd_i = zeros(Float64,ny-1,nx)
sd_j = zeros(Float64,ny,nx-1)

for i = 2 : nx-1
  for j = 1 : ny-1
    dist_i[j,i] = norm(meshcell[j,i].cellcenter-meshcell[j,i-1].cellcenter)
  end
end

for j = 1:ny-1
    dist_i[j,1] = dist_i[j,2]
    dist_i[j,nx] = dist_i[j,nx-1]
end

for i =  1 : nx-1
  for j = 2 : ny-1
  dist_j[j,i] = norm(meshcell[j,i].cellcenter-meshcell[j-1,i].cellcenter)
  end
end

for i = 1:nx-1
  dist_j[1,i] = dist_j[2,i]
  dist_j[ny,i] = dist_j[ny-1,i]
end

# calculate the S/dist
for i = 1 : nx-1
  for j = 1 : ny-1
    sd_i[j,i] = 1.0/dist_i[j,i]*meshcell[j,i].length41
    sd_j[j,i] = 1.0/dist_j[j,i]*meshcell[j,i].length12
  end
end
for j = 1:ny-1
   sd_i[j,nx] = 1.0/dist_i[j,nx]*meshcell[j,nx-1].length23
end
for i = 1:nx-1
   sd_j[ny,i] = 1.0/dist_j[ny,i]*meshcell[ny-1,i].length34
end

# ========== Initial Condition Setting ====================
function init(ny::Int64,nx::Int64)
  # initialize the stored value .

 global  u = zeros(Float64,ny+1,nx+1)
 global  v = zeros(Float64,ny+1,nx+1)
 global  p = ones(Float64,ny+1,nx+1)
 global  p0 = ones(Float64,ny+1,nx+1) # store the previous iter result
 # store the u_s = u_star, v_s = v_star
 global  u_s = zeros(Float64,ny+1,nx+1)
 global  v_s = zeros(Float64,ny+1,nx+1)
 #=
 # tyler green init
 for j = 1:ny-1
  for i = 1:nx-1
  u[j+1,i+1] = -(cos(meshcell[j,i].cellcenter[1])*
                 sin(meshcell[j,i].cellcenter[2]))
  v[j+1,i+1] =  (sin(meshcell[j,i].cellcenter[1])*
                 cos(meshcell[j,i].cellcenter[2]))
  p[j+1,i+1] = 1.0 - 1/4*(cos(2*meshcell[j,i].cellcenter[1])
                       +cos(2*meshcell[j,i].cellcenter[2]))
   end
  end
  =#
end
# =========== Boundary Condition Setting ==================
function bcs(ny::Int64,nx::Int64)
 #=
 # ======= Tyler Green problem==========
 # Periodic Boundary condition
  for j = 2: ny
   global u[j,1] = u[j,nx]
   global u_s[j,1] = u_s[j,nx]
   global v[j,1] = v[j,nx]
   global v_s[j,1] = v_s[j,nx]
   global p[j,1] = p[j,nx]
   global p0[j,1] = p0[j,nx]

   global u[j,nx+1] = u[j,2]
   global u_s[j,nx+1] = u_s[j,2]
   global v[j,nx+1] = v[j,2]
   global v_s[j,nx+1] = v_s[j,2]
   global p[j,nx+1] = p[j,2]
   global p0[j,nx+1] = p0[j,2]
  end

  for i = 2: nx
   global u[1,i] = u[ny,i]
   global u_s[1,i] = u_s[ny,i]
   global v[1,i] = v[ny,i]
   global v_s[1,i] = v_s[ny,i]
   global p[1,i] = p[ny,i]
   global p0[1,i] = p0[ny,i]

   global u[ny+1,i] = u[2,i]
   global u_s[ny+1,i] = u_s[2,i]
   global v[ny+1,i] = v[2,i]
   global v_s[ny+1,i] = v_s[2,i]
   global p[ny+1,i] = p[2,i]
   global p0[ny+1,i] = p0[2,i]
  end
 =#
 #
 # ======  cavity problem ==============
 #  Dirichlet boundary Condition bc = 0
 #  left & right
   for j = 2:ny
      global u[j,1] = -1*u[j,2]
      global u_s[j,1] = -1*u_s[j,2]

      global u[j,nx+1] = -1*u[j,nx]
      global u_s[j,nx+1] = -1*u_s[j,nx]

      global v[j,1] = -1*v[j,2]
      global v_s[j,1] = -1*v_s[j,2]

      global v[j,nx+1] =-1*v[j,nx]
      global v_s[j,nx+1] = -1*v_s[j,nx]
   end
   # bottom & top
   for i = 2:nx
      global u[1,i] = -1*u[2,i]
      global u_s[1,i] =  -1*u_s[2,i]

      global u[ny+1,i] = U*2.0 - u[ny,i]
      global u_s[ny+1,i] = U*2.0 - u_s[ny,i]

      global v[1,i] = -1*v[2,i]
      global v_s[1,i] = -1*v_s[2,i]

      global v[ny+1,i] = -1*v[ny,i]
      global v_s[ny+1,i] = -1*v_s[ny,i]
   end

 # Neumann BC of pressure
 for j = 2:ny
    global p[j,1] = p[j,2]
    global p[j,nx+1] = p[j,nx]
 end
  for i = 2:nx
    global p[1,i] = p[2,i]
    global p[ny+1,i] = p[ny,i]
  end
  #

end
function pressure_correction()
    global nx2 = zeros(Float64,ny-1,nx-1)
    global ny2 = zeros(Float64,ny-1,nx-1)
    global nxy = zeros(Float64,ny-1,nx-1)

    for j = 1 : ny-1
      for i = 1 : nx-1
      nx2[j,i] =(meshcell[j,i].out_norm12[1]^2*meshcell[j,i].length12^2
               + meshcell[j,i].out_norm23[1]^2*meshcell[j,i].length23^2
               + meshcell[j,i].out_norm34[1]^2*meshcell[j,i].length34^2
               + meshcell[j,i].out_norm41[1]^2*meshcell[j,i].length41^2)

      ny2[j,i] =(meshcell[j,i].out_norm12[2]^2*meshcell[j,i].length12^2
               + meshcell[j,i].out_norm23[2]^2*meshcell[j,i].length23^2
               + meshcell[j,i].out_norm34[2]^2*meshcell[j,i].length34^2
               + meshcell[j,i].out_norm41[2]^2*meshcell[j,i].length41^2)

      nxy[j,i] =(meshcell[j,i].out_norm12[1]*meshcell[j,i].out_norm12[2]*meshcell[j,i].length12^2
               + meshcell[j,i].out_norm23[1]*meshcell[j,i].out_norm23[2]*meshcell[j,i].length23^2
               + meshcell[j,i].out_norm34[1]*meshcell[j,i].out_norm34[2]*meshcell[j,i].length34^2
               + meshcell[j,i].out_norm41[1]*meshcell[j,i].out_norm41[2]*meshcell[j,i].length41^2)
      end
    end
end
# ====== test section ======
init(ny,nx)
bcs(ny,nx)
pressure_correction()
# spatial(u,v,p,ny,nx)
# ====== Maain Loop ====================
t_min = 20.0
t_max = 25.0
t_step = 0.0005
#t_step = 0.0005 tyler green
t_end = convert(Int64, (t_max-t_min)/t_step)
t = zeros(t_end)
K = zeros(t_end)
iter = 0

Nu_1 = 0.0
Nv_1 = 0.0
Vu_1 = 0.0
Vv_1 = 0.0

vn_ss = zeros(ny-1,nx-1)

for it = 1 : t_end

  t[it] = t_step*it
  print("time = $(t_min+(it*t_step))",'\n')

  for j = 1 : ny-1
    for i = 1 : nx-1
     # surface u v derivative
     u12 = (u[j+1,i+1]+u[j,i+1])/2
     u23 = (u[j+1,i+1]+u[j+1,i+2])/2
     u34 = (u[j+1,i+1]+u[j+2,i+1])/2
     u41 = (u[j+1,i+1]+u[j+1,i])/2

     v12 = (v[j+1,i+1]+v[j,i+1])/2
     v23 = (v[j+1,i+1]+v[j+1,i+2])/2
     v34 = (v[j+1,i+1]+v[j+2,i+1])/2
     v41 = (v[j+1,i+1]+v[j+1,i])/2
    # ==============================
     #=
     vn12 = u12*meshcell[j,i].out_norm12[1] + v12*meshcell[j,i].out_norm12[2]
     vn23 = u23*meshcell[j,i].out_norm23[1] + v23*meshcell[j,i].out_norm23[2]
     vn34 = u34*meshcell[j,i].out_norm34[1] + v34*meshcell[j,i].out_norm34[2]
     vn41 = u41*meshcell[j,i].out_norm41[1] + v41*meshcell[j,i].out_norm41[2]
     =#
     vn12 = dot([u12;v12],meshcell[j,i].out_norm12)
     vn23 = dot([u23;v23],meshcell[j,i].out_norm23)
     vn34 = dot([u34;v34],meshcell[j,i].out_norm34)
     vn41 = dot([u41;v41],meshcell[j,i].out_norm41)

     u_vn_s = (u12*vn12*meshcell[j,i].length12
              +u23*vn23*meshcell[j,i].length23
              +u34*vn34*meshcell[j,i].length34
              +u41*vn41*meshcell[j,i].length41)
     v_vn_s = (v12*vn12*meshcell[j,i].length12
              +v23*vn23*meshcell[j,i].length23
              +v34*vn34*meshcell[j,i].length34
              +v41*vn41*meshcell[j,i].length41)

    Nu = -1/meshcell[j,i].area*u_vn_s
    Nv = -1/meshcell[j,i].area*v_vn_s

    # surface derivative flux

    u_n12_s = (u[j,i+1]-u[j+1,i+1])*sd_j[j,i]
    u_n23_s = (u[j+1,i+2]-u[j+1,i+1])*sd_i[j,i+1]
    u_n34_s = (u[j+2,i+1]-u[j+1,i+1])*sd_j[j+1,i]
    u_n41_s = (u[j+1,i]-u[j+1,i+1])*sd_i[j,i]

    v_n12_s = (v[j,i+1]-v[j+1,i+1])*sd_j[j,i]
    v_n23_s = (v[j+1,i+2]-v[j+1,i+1])*sd_i[j,i+1]
    v_n34_s = (v[j+2,i+1]-v[j+1,i+1])*sd_j[j+1,i]
    v_n41_s = (v[j+1,i]-v[j+1,i+1])*sd_i[j,i]

    Vu = 1/meshcell[j,i].area*nu*(u_n12_s+u_n23_s+u_n34_s+u_n41_s)
    Vv = 1/meshcell[j,i].area*nu*(v_n12_s+v_n23_s+v_n34_s+v_n41_s)

    # AB2
    #u_s = u_star
    #=
    if it == 1
      global Nu_1 = Nu
      global Nv_1 = Nv
      global Vu_1 = Vu
      global Vv_1 = Vv
    end
    =#
    u_s[j+1,i+1] = u[j+1,i+1] + t_step*(1.5*(Nu+Vu)-0.5*(Nu_1+Vu_1))
    v_s[j+1,i+1] = v[j+1,i+1] + t_step*(1.5*(Nv+Vv)-0.5*(Nv_1+Vv_1))

    Nu_1 = Nu
    Nv_1 = Nv
    Vu_1 = Vu
    Vv_1 = Vv

  end
 end
 bcs(nx,ny)

 # end of calculate the u*.v* value at cellcenter

 # Pressure calculation parts

  for j = 1:ny-1
    for i = 1:nx-1

    # u_s12 ,23 ,34 ,41 = u* on surface
    u_s12 = (u_s[j+1,i+1]+u_s[j,i+1])/2
    u_s23 = (u_s[j+1,i+1]+u_s[j+1,i+2])/2
    u_s34 = (u_s[j+1,i+1]+u_s[j+2,i+1])/2
    u_s41 = (u_s[j+1,i+1]+u_s[j+1,i])/2

    v_s12 = (v_s[j+1,i+1]+v_s[j,i+1])/2
    v_s23 = (v_s[j+1,i+1]+v_s[j+1,i+2])/2
    v_s34 = (v_s[j+1,i+1]+v_s[j+2,i+1])/2
    v_s41 = (v_s[j+1,i+1]+v_s[j+1,i])/2

    #  vn = unx + vny ======================
    vn_s12 = dot([u_s12;v_s12],meshcell[j,i].out_norm12)
    vn_s23 = dot([u_s23;v_s23],meshcell[j,i].out_norm23)
    vn_s34 = dot([u_s34;v_s34],meshcell[j,i].out_norm34)
    vn_s41 = dot([u_s41;v_s41],meshcell[j,i].out_norm41)

    vn_ss[j,i] = rho/t_step*(vn_s12*meshcell[j,i].length12
                             +vn_s23*meshcell[j,i].length23
                             +vn_s34*meshcell[j,i].length34
                             +vn_s41*meshcell[j,i].length41)
      end
    end

  # surface pressure derivative flux
  for iter = 1:10000
    res = 0.0
    for j = 1:ny-1
    for i = 1:nx-1
     # G_S
     p0[j+1,i+1] = p[j+1,i+1]
     p[j+1,i+1] = -1*((vn_ss[j,i]-p[j+1,i+2]*sd_i[j,i+1]-p[j+2,i+1]*sd_j[j+1,i]
                       -p[j+1,i]*sd_i[j,i]-p[j,i+1]*sd_j[j,i])/
                       (sd_i[j,i]+sd_i[j,i+1]+sd_j[j,i]+sd_j[j+1,i]))
     # SOR
      #w = 1.7 for tyler green
      w = 1.985
      p[j+1,i+1] = p0[j+1,i+1] + w*(p[j+1,i+1]-p0[j+1,i+1])

      res = res +  (p[j+1,i+1]-p0[j+1,i+1])^2
     end
    end
    bcs(ny,nx)

     if res^0.5 <= 0.0001
       print("res = $(res^0.5)",'\n')
       print("iter = $iter",'\n')
       print("iter finished",'\n')
       break
     end
 end

    # Pressure correction to link u_s to u at next time step
    for j = 1:ny-1
      for i = 1 : nx-1

      pnx12 = ((p[j,i+1]-p[j+1,i+1])*sd_j[j,i]*meshcell[j,i].out_norm12[1]
               *meshcell[j,i].length12)
      pnx23 = ((p[j+1,i+2]-p[j+1,i+1])*sd_i[j,i+1]*meshcell[j,i].out_norm23[1]
                *meshcell[j,i].length23)
      pnx34 = ((p[j+2,i+1]-p[j+1,i+1])*sd_j[j+1,i]*meshcell[j,i].out_norm34[1]
                *meshcell[j,i].length34)
      pnx41 = ((p[j+1,i]-p[j+1,i+1])*sd_i[j,i]*meshcell[j,i].out_norm41[1]
                *meshcell[j,i].length34)
      pnx = pnx12 + pnx23 + pnx34 + pnx41

      pny12 = ((p[j,i+1]-p[j+1,i+1])*sd_j[j,i]*meshcell[j,i].out_norm12[2]
               *meshcell[j,i].length12)
      pny23 = ((p[j+1,i+2]-p[j+1,i+1])*sd_i[j,i+1]*meshcell[j,i].out_norm23[2]
                *meshcell[j,i].length23)
      pny34 = ((p[j+2,i+1]-p[j+1,i+1])*sd_j[j+1,i]*meshcell[j,i].out_norm34[2]
                *meshcell[j,i].length34)
      pny41 = ((p[j+1,i]-p[j+1,i+1])*sd_i[j,i]*meshcell[j,i].out_norm41[2]
                *meshcell[j,i].length41)
      pny = pny12 + pny23 + pny34 + pny41

      (px,py) = [nx2[j,i] nxy[j,i];nxy[j,i] ny2[j,i]]\[pnx;pny]
      #print("px = $px")
      #print("py = $py")
      global u[j+1,i+1] = u_s[j+1,i+1]-t_step/rho*px
      global v[j+1,i+1] = v_s[j+1,i+1]-t_step/rho*py

      #K[it] = K[it] + (u[j+1,i+1]^2 + v[j+1,i+1]^2)/2*meshcell[j,i].area
      end
    end
    bcs(ny,nx)


end

#K = K/(4*pi^2)


# plot and figure

#using PyPlot
#t_max = 0
x_coord = zeros((ny-1),(nx-1))
y_coord = zeros((ny-1),(nx-1))
for j = 1:ny-1
  for i = 1:nx-1
  x_coord[j,i] = meshcell[j,i].cellcenter[1]
  y_coord[j,i] = meshcell[j,i].cellcenter[2]
  end
end
#=
# 2-d Tyler Green Porblem
#t_max = 0.0001
v_tyler = zeros((ny-1),(nx-1))
u_tyler = zeros((ny-1),(nx-1))
p_tyler = zeros((ny-1),(nx-1))
for j = 1:ny-1
  for i = 1:nx-1
   u_tyler[j,i] = -cos(x_coord[j,i])*sin(y_coord[j,i])*exp(-2*nu*t_max)
   v_tyler[j,i] =  sin(x_coord[j,i])*cos(y_coord[j,i])*exp(-2*nu*t_max)
   p_tyler[j,i] = 1.0 - 1/4*(cos(2*x_coord[j,i])+cos(2*y_coord[j,i]))*exp(-4*nu*t_max)
  end
end

K_T = zeros(t_end)
for i = 1:t_end
 K_T[i] = 1/4*exp(-4*nu*i*t_step)
end
=#

using Plots
pyplot()

contour(x_coord,y_coord,u[2:ny,2:nx])
#contour(x_coord,y_coord,u_tyler)
contour(x_coord,y_coord,v[2:ny,2:nx])
#contour(x_coord,y_coord,v_tyler)
contour(x_coord,y_coord,p[2:ny,2:nx])
#contour(x_coord,y_coord,p_tyler)

#=
for it = 1 : t_end
  t[it] = t_step*it
end
plot(t,K,label ="CFD",
     title = "Re = $Re KE decay",
     xlabel = "t", ylabel = "KE")
plot!(t,K_T,label = "Theoretical")

#png("re1000")
=#
#
# Data restart process
using DelimitedFiles
writedlm("p data of $t_max re400",p)
writedlm("u data of $t_max re400",u)
writedlm("v data of $t_max re400",v)


u = readdlm("u data of 10.0 re400")
v = readdlm("v data of 10.0 re400")
p = readdlm("p data of 10.0 re400")
