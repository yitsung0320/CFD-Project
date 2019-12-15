 # CFD projet main file

 # 2-D incompressible Navier Stokes Problem
 # Coordinate generate / mesh builder

include("Meshreader.jl")
using LinearAlgebra
using DelimitedFiles
using Plots
pyplot()


function tylergreen(x::Real)
  global Lx, Ly = 2*pi , 2*pi
  L = Lx
  global Re = x
  global nx, ny = 100,100
  global rho = 1.0
  global U = 1.0
  global mu = rho*U*L/Re
  global nu = mu/rho
  global problem_type = 1
  global w = 1.7
end

function cavity(x::Real)
  global Lx, Ly = 1 , 1
  L = Lx
  global Re = x
  global nx, ny = 129,129
  global rho = 1
  global U = 1
  global mu = rho*U*L/Re
  global nu = mu/rho
  global problem_type = 2
  global w = 1.985
end

function grid_generator!(nx::Int64,ny::Int64)
   global grid = zeros(Float64,nx*ny,2)
   for i = 1 : nx
    for j = 1 : ny
     k = (j-1)*nx+i
     grid[k,1] = Lx/(nx-1)*(i-1)
     end
   end
   for i = 1 : nx
    for j = 1 : ny
     k = (j-1)*nx+i
     grid[k,2] = Ly/(ny-1)*(j-1)
    end
   end
end
# ==  extract the needed variables based from geometry ====
function geo_factor(nx::Int64,ny::Int64)
  global dist_i = ones(Float64,ny-1,nx) # the x dir distance between cellcenter
  global dist_j = ones(Float64,ny,nx-1) # the y dir distance between cellcenter
  global sd_i = zeros(Float64,ny-1,nx) # surface length23 or 41/dist_i
  global sd_j = zeros(Float64,ny,nx-1) # surface length12 pr 34/dist_j
  global sum_sd = zeros(Float64,ny-1,nx-1) # sum of sd_i sd_j for each real cell

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
  for i = 1 : nx-1
  for j = 1 : ny-1
      sum_sd[j,i] =sd_i[j,i]+sd_i[j,i+1]+sd_j[j,i]+sd_j[j+1,i]
  end
  end
end

# ========== Initial Condition Setting ====================
function init(nx::Int64,ny::Int64,problem_type::Int64)
  # initialize the stored value .

 global  u = zeros(Float64,ny+1,nx+1)
 global  v = zeros(Float64,ny+1,nx+1)
 global  p = ones(Float64,ny+1,nx+1)
 # store the previous iter result
 global  p0 = ones(Float64,ny+1,nx+1)
 # store the u_s = u_star, v_s = v_star
 global  u_s = zeros(Float64,ny+1,nx+1)
 global  v_s = zeros(Float64,ny+1,nx+1)
 if problem_type == 1
  # tyler green init
   for i = 1:nx-1
   for j = 1:ny-1
    u[j+1,i+1] = -(cos(meshcell[j,i].cellcenter[1])*
                   sin(meshcell[j,i].cellcenter[2]))
    u_s[j+1,i+1] = u[j+1,i+1]
    v[j+1,i+1] =  (sin(meshcell[j,i].cellcenter[1])*
                   cos(meshcell[j,i].cellcenter[2]))
    v_s[j+1,i+1] = v[j+1,i+1]
    p[j+1,i+1] = 1.0 - 1/4*(cos(2*meshcell[j,i].cellcenter[1])
                          +cos(2*meshcell[j,i].cellcenter[2]))
    end
    end
  else
  end
end

# =========== Boundary Condition Setting ==================
function bcs(nx::Int64,ny::Int64,problem_type::Int64,)

 if problem_type == 1
 # ======= Tyler Green problem==========
 # Periodic Boundary condition
   for j = 2: ny
     u[j,1]   = u[j,nx]
     u_s[j,1] = u_s[j,nx]
     v[j,1]   = v[j,nx]
     v_s[j,1] = v_s[j,nx]
     p[j,1]   = p[j,nx]
     p0[j,1]  = p0[j,nx]

     u[j,nx+1]   = u[j,2]
     u_s[j,nx+1] = u_s[j,2]
     v[j,nx+1]   = v[j,2]
     v_s[j,nx+1] = v_s[j,2]
     p[j,nx+1]   = p[j,2]
     p0[j,nx+1]  = p0[j,2]
  end
   for i = 2: nx
     u[1,i]   = u[ny,i]
     u_s[1,i] = u_s[ny,i]
     v[1,i]   = v[ny,i]
     v_s[1,i] = v_s[ny,i]
     p[1,i]   = p[ny,i]
     p0[1,i]  = p0[ny,i]

     u[ny+1,i]   = u[2,i]
     u_s[ny+1,i] = u_s[2,i]
     v[ny+1,i]   = v[2,i]
     v_s[ny+1,i] = v_s[2,i]
     p[ny+1,i]   = p[2,i]
     p0[ny+1,i]  = p0[2,i]
  end
  else
 # ======  cavity problem ==============
 #  Dirichlet boundary Condition bc = 0
 #  left & right
   for j = 2:ny
       u[j,1] = -1*u[j,2]
       u_s[j,1] = -1*u_s[j,2]

       u[j,nx+1] = -1*u[j,nx]
       u_s[j,nx+1] = -1*u_s[j,nx]

       v[j,1] = -1*v[j,2]
       v_s[j,1] = -1*v_s[j,2]

       v[j,nx+1] =-1*v[j,nx]
       v_s[j,nx+1] = -1*v_s[j,nx]
   end
   # bottom & top
    for i = 2:nx
       u[1,i] = -1*u[2,i]
       u_s[1,i] =  -1*u_s[2,i]

       u[ny+1,i] = U*2.0 - u[ny,i]
       u_s[ny+1,i] = U*2.0 - u_s[ny,i]

       v[1,i] = -1*v[2,i]
       v_s[1,i] = -1*v_s[2,i]

       v[ny+1,i] = -1*v[ny,i]
       v_s[ny+1,i] = -1*v_s[ny,i]
   end
   # Neumann BC of pressure
   for j = 2:ny
       p[j,1] = p[j,2]
       p[j,nx+1] = p[j,nx]
   end
   for i = 2:nx
       p[1,i] = p[2,i]
       p[ny+1,i] = p[ny,i]
   end
 end #end if
end

# =========== pressure coorrection varaible ===============
function pressure_correction()
    global nx2 = zeros(Float64,ny-1,nx-1)
    global ny2 = zeros(Float64,ny-1,nx-1)
    global nxy = zeros(Float64,ny-1,nx-1)

    for i = 1 : nx-1
    for j = 1 : ny-1
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

# ====== precondition =================
#tylergreen(10)
cavity(1000)
grid_generator!(nx,ny)
meshcell = cellproperty(grid,nx,ny)
geo_factor(nx,ny)

init(nx,ny,problem_type)
bcs(nx,ny,problem_type)
pressure_correction()

# ====== Maain Loop  setting ====================
global t_min = 0.0
global t_max = 25.0
global t_step = 0.005
#t_step = 0.001 tyler green
global t_end = convert(Int64, (t_max-t_min)/t_step)
global t = zeros(t_end)
global K = zeros(t_end)

# ======== main code calculation stage function definition ==============
function v_prediction(Nu::Array{Float64,2},Nu_1::Array{Float64,2},
                      Nv::Array{Float64,2},Nv_1::Array{Float64,2},
                      Vu::Array{Float64,2},Vu_1::Array{Float64,2},
                      Vv::Array{Float64,2},Vv_1::Array{Float64,2})

   for i = 1 : nx-1
   for j = 1 : ny-1
     # ==== surface u v derivative ===================
     u12 = (u[j+1,i+1]+u[j,i+1])/2
     u23 = (u[j+1,i+1]+u[j+1,i+2])/2
     u34 = (u[j+1,i+1]+u[j+2,i+1])/2
     u41 = (u[j+1,i+1]+u[j+1,i])/2

     v12 = (v[j+1,i+1]+v[j,i+1])/2
     v23 = (v[j+1,i+1]+v[j+1,i+2])/2
     v34 = (v[j+1,i+1]+v[j+2,i+1])/2
     v41 = (v[j+1,i+1]+v[j+1,i])/2
     # ===== surface normal velocity u*nx+v*ny =======
     vn12 = u12*meshcell[j,i].out_norm12[1] + v12*meshcell[j,i].out_norm12[2]
     vn23 = u23*meshcell[j,i].out_norm23[1] + v23*meshcell[j,i].out_norm23[2]
     vn34 = u34*meshcell[j,i].out_norm34[1] + v34*meshcell[j,i].out_norm34[2]
     vn41 = u41*meshcell[j,i].out_norm41[1] + v41*meshcell[j,i].out_norm41[2]

     u_vn_s = (u12*vn12*meshcell[j,i].length12+u23*vn23*meshcell[j,i].length23
              +u34*vn34*meshcell[j,i].length34+u41*vn41*meshcell[j,i].length41)
     v_vn_s = (v12*vn12*meshcell[j,i].length12+v23*vn23*meshcell[j,i].length23
              +v34*vn34*meshcell[j,i].length34+v41*vn41*meshcell[j,i].length41)

     Nu[j,i] = -1.0/meshcell[j,i].area*u_vn_s
     Nv[j,i] = -1.0/meshcell[j,i].area*v_vn_s

     #  ====== surface derivative flux ======================
     u_n12_s = (u[j,i+1]-u[j+1,i+1])*sd_j[j,i]
     u_n23_s = (u[j+1,i+2]-u[j+1,i+1])*sd_i[j,i+1]
     u_n34_s = (u[j+2,i+1]-u[j+1,i+1])*sd_j[j+1,i]
     u_n41_s = (u[j+1,i]-u[j+1,i+1])*sd_i[j,i]

     v_n12_s = (v[j,i+1]-v[j+1,i+1])*sd_j[j,i]
     v_n23_s = (v[j+1,i+2]-v[j+1,i+1])*sd_i[j,i+1]
     v_n34_s = (v[j+2,i+1]-v[j+1,i+1])*sd_j[j+1,i]
     v_n41_s = (v[j+1,i]-v[j+1,i+1])*sd_i[j,i]

     Vu[j,i] = 1.0/meshcell[j,i].area*nu*(u_n12_s+u_n23_s+u_n34_s+u_n41_s)
     Vv[j,i] = 1.0/meshcell[j,i].area*nu*(v_n12_s+v_n23_s+v_n34_s+v_n41_s)
   end
   end
   #=
   #  ======== u_s = u_star AB2 method for nonlinear and vicous ============
   for i = 1:nx-1
   for j = 1:ny-1
     u_s[j+1,i+1] = u[j+1,i+1] + t_step*(1.5*(Nu[j,i]+Vu[j,i])-0.5*(Nu_1[j,i]+Vu_1[j,i]))
     v_s[j+1,i+1] = v[j+1,i+1] + t_step*(1.5*(Nv[j,i]+Vv[j,i])-0.5*(Nv_1[j,i]+Vv_1[j,i]))

     Nu_1[j,i] = Nu[j,i]
     Nv_1[j,i] = Nv[j,i]
     Vu_1[j,i] = Vu[j,i]
     Vv_1[j,i] = Vv[j,i]
   end
   end
    =#
    #
    #  ======== using Crank -Nicolson Method for viscous term ===============
    # G_S to solve iteration first and then SOR

     u_s0 = zeros(Float64,ny,nx) # used to store previous calculation
     v_s0 = zeros(Float64,ny,nx)

     # SOR setting gain
     w = 1.7
     for iter = 1:10000
        resu = 0.0
        resv = 0.0
      for i = 1:nx-1
      for j = 1:ny-1
       u_s0[j+1,i+1] = u_s[j+1,i+1]
       v_s0[j+1,i+1] = v_s[j+1,i+1]

       RHSu = u[j+1,i+1] + t_step*(1.5*Nu[j,i]-0.5*Nu_1[j,i]+0.5*Vu[j,i])
       RHSv = v[j+1,i+1] + t_step*(1.5*Nv[j,i]-0.5*Nv_1[j,i]+0.5*Vv[j,i])

        u_s[j+1,i+1] = ((RHSu + 0.5*t_step*nu/meshcell[j,i].area*
                    (sd_i[j,i+1]*u_s[j+1,i+2] + sd_i[j,i]*u_s[j+1,i]
                    +sd_j[j+1,i]*u_s[j+2,i+1] + sd_j[j,i]*u_s[j,i+1]))/
                    (1 + 0.5*t_step*nu/meshcell[j,i].area*sum_sd[j,i]))
        v_s[j+1,i+1] = ((RHSv + 0.5*t_step*nu/meshcell[j,i].area*
                    (sd_i[j,i+1]*v_s[j+1,i+2] + sd_i[j,i]*v_s[j+1,i]
                    +sd_j[j+1,i]*v_s[j+2,i+1] + sd_j[j,i]*v_s[j,i+1]))/
                    (1 + 0.5*t_step*nu/meshcell[j,i].area*sum_sd[j,i]))

       #  ======= SOR iteration =============
         u_s[j+1,i+1] = u_s0[j+1,i+1] + w*(u_s[j+1,i+1]-u_s0[j+1,i+1])
         v_s[j+1,i+1] = v_s0[j+1,i+1] + w*(v_s[j+1,i+1]-v_s0[j+1,i+1])

         resu = resu + (u_s[j+1,i+1]-u_s0[j+1,i+1])^2
         resv = resv + (v_s[j+1,i+1]-v_s0[j+1,i+1])^2
        end
        end

       if resu^0.5 <=1e-4 && resv^0.5<=1e-4
         print("resu = $(resu^0.5)",'\n')
         print("resv = $(resv^0.5)",'\n')
         print("iter = $iter",'\n')
         print("velocity prediction iter finished",'\n')
       break
      end
    end

     for i = 1:nx-1
     for j = 1:ny-1
        Nu_1[j,i] = Nu[j,i]
        Nv_1[j,i] = Nv[j,i]
      end
      end
     #
end
function p_calculation()

  vn_ss = zeros(Float64,ny-1,nx-1)

   for i = 1:nx-1
   for j = 1:ny-1
    # u_s12 ,23 ,34 ,41 = u* on surface
    u_s12 = (u_s[j+1,i+1]+u_s[j,i+1])/2
    u_s23 = (u_s[j+1,i+1]+u_s[j+1,i+2])/2
    u_s34 = (u_s[j+1,i+1]+u_s[j+2,i+1])/2
    u_s41 = (u_s[j+1,i+1]+u_s[j+1,i])/2

    v_s12 = (v_s[j+1,i+1]+v_s[j,i+1])/2
    v_s23 = (v_s[j+1,i+1]+v_s[j+1,i+2])/2
    v_s34 = (v_s[j+1,i+1]+v_s[j+2,i+1])/2
    v_s41 = (v_s[j+1,i+1]+v_s[j+1,i])/2

    #  vn_s = u_snx + v_sny ======================
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

    #print(vn_ss[3,3])
  # surface pressure derivative flux
   for iter = 1:10000
    res = 0.0
    for i = 1:nx-1
    for j = 1:ny-1
     # G_S
     p0[j+1,i+1] = p[j+1,i+1]
     p[j+1,i+1] = -1*((vn_ss[j,i]-p[j+1,i+2]*sd_i[j,i+1]-p[j+2,i+1]*sd_j[j+1,i]
                       -p[j+1,i]*sd_i[j,i]-p[j,i+1]*sd_j[j,i])/
                       (sd_i[j,i]+sd_i[j,i+1]+sd_j[j,i]+sd_j[j+1,i]))
     # SOR
      #w = 1.7 #for tyler green
      w = 1.985
      p[j+1,i+1] = p0[j+1,i+1] + w*(p[j+1,i+1]-p0[j+1,i+1])

      res = res +  (p[j+1,i+1]-p0[j+1,i+1])^2
    end
    end
    bcs(nx,ny,problem_type)

     if res^0.5 <= 0.0001
       print("resp = $(res^0.5)",'\n')
       print("iter = $iter",'\n')
       print("pressure iter finished",'\n')
       break
     end
   end
end
function p_correction(it::Int64)
   for i = 1 : nx-1
   for j = 1 : ny-1

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

    u[j+1,i+1] = u_s[j+1,i+1]-t_step/rho*px
    v[j+1,i+1] = v_s[j+1,i+1]-t_step/rho*py

    K[it] = K[it] + 1/(4*pi^2)*(u[j+1,i+1]^2+v[j+1,i+1]^2)/2*meshcell[j,i].area
    end
    end
end

#  =================== plot and figure =============================

function tylergreen_plot(t_max::Float64)

  x_coord = zeros((ny-1),(nx-1))
  y_coord = zeros((ny-1),(nx-1))
  for j = 1:ny-1
  for i = 1:nx-1
    x_coord[j,i] = meshcell[j,i].cellcenter[1]
    y_coord[j,i] = meshcell[j,i].cellcenter[2]
  end
  end
  # 2-d Tyler Green Porblem
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
     K_T[i] = 1/4*exp(-4*nu*(t_min+i*t_step))
   end
   plot(t,K,label ="CFD",
      title = "Re = $Re KE decay",
      xlabel = "t", ylabel = "KE")
   plot!(t,K_T,label = "Theoretical")
   png("Tyler green K comparison")

   contour(x_coord,y_coord,u[2:ny,2:nx],xlabel = "x",ylabel = "y")
   png(" u contour tyler at $t_max Re$Re")
   contour(x_coord,y_coord,u_tyler,xlabel = "x",ylabel = "y")
   png(" u contour tylergreen ana at $t_max Re$Re")
   contour(x_coord,y_coord,v[2:ny,2:nx],xlabel = "x",ylabel = "y")
   png(" v contour tyler at $t_max Re$Re")
   contour(x_coord,y_coord,v_tyler,xlabel = "x",ylabel = "y")
   png(" v contour tylergreen ana at $t_max Re$Re")
   contour(x_coord,y_coord,p[2:ny,2:nx],xlabel = "x",ylabel = "y",fill = true)
   png(" p contour tyler at $t_max Re$Re")
   contour(x_coord,y_coord,p_tyler,xlabel = "x",ylabel = "y",fill = true)
   png(" p contour tylergreen ana at $t_max Re$Re")

   x_1d = zeros((nx-1))
   y_1d = zeros((ny-1))
   for i = 1:nx-1
    x_1d[i] = x_coord[1,i]
   end
   for j = 1:ny-1
     y_1d[j] = y_coord[j,1]
   end
   PyPlot.figure()
   PyPlot.xlabel("x")
   PyPlot.ylabel("y")
   PyPlot.title("Streamline")
   PyPlot.streamplot(x_1d,y_1d,u[2:ny,2:nx],v[2:ny,2:nx],density=0.7)
   PyPlot.savefig("tyler green streamplot")
end
function tylergreen_result(t_max::Float64)

   writedlm("tyler p data of $t_max Re$Re",p)
   writedlm("tyler u data of $t_max Re$Re",u)
   writedlm("tyler v data of $t_max Re$Re",v)
end
#ylergreen_result()
#tylergreen_plot()
function cavity_result(t_max::Float64)

  # written the original calulate data
  writedlm("cavity p data of $t_max Re$Re",p)
  writedlm("cavity u data of $t_max Re$Re",u)
  writedlm("cavity v data of $t_max Re$Re",v)
  # for cavity comparison ,
  # the data points is set to be on the grids instead of cell
  x_coord = zeros(ny,nx)
  y_coord = zeros(ny,nx)
  u_c = zeros(ny,nx)
  v_c = zeros(ny,nx)
  p_c = zeros(ny,nx)

  for i = 1:nx-1
  for j = 1:ny-1
    x_coord[j,i] = meshcell[j,i].point1[1]
    y_coord[j,i] = meshcell[j,i].point1[2]
    x_coord[j+1,i] = meshcell[j,i].point4[1]
    y_coord[j+1,i] = meshcell[j,i].point4[2]
    x_coord[j+1,i+1] = meshcell[j,i].point3[1]
    y_coord[j+1,i+1] = meshcell[j,i].point3[2]
    x_coord[j,i+1] = meshcell[j,i].point2[1]
    y_coord[j,i+1] = meshcell[j,i].point2[2]
  end
  end

  # ======== assign corner value
  u[1,1] = u[2,2]
  v[1,1] = v[2,2]
  p[1,1] = p[2,2]

  u[ny+1,nx+1] = u[ny,nx]
  v[ny+1,nx+1] = v[ny,nx]
  p[ny+1,nx+1] = p[ny,nx]

  u[ny+1,1] = u[ny,2]
  v[ny+1,1] = v[ny,2]
  p[ny+1,1] = p[ny,2]

  u[1,nx+1] = u[2,nx]
  v[1,nx+1] = v[2,nx]
  p[1,nx+1] = p[2,nx]
  # ========= but the 4 corner points value is meaningless
  for i = 1:nx
  for j = 1:ny
   u_c[j,i] = 0.25*(u[j,i]+u[j,i+1]+u[j+1,i]+u[j+1,i+1])
   v_c[j,i] = 0.25*(v[j,i]+v[j,i+1]+v[j+1,i]+v[j+1,i+1])
   p_c[j,i] = 0.25*(p[j,i]+p[j,i+1]+p[j+1,i]+p[j+1,i+1])
  end
  end

  #writedlm("cavity p grid data of $t_max Re$Re",p_c)
  #writedlm("cavity u grid data of $t_max Re$Re",u_c)
  #writedlm("cavity v grid data of $t_max Re$Re",v_c)

  u_data = zeros(17)
  v_data = zeros(17)
  g_series = [129,126,125,124,123,110,95,80,65,59,37,23,14,10,9,8,1]
  for i = 1:17
    u_data[i] = u_c[g_series[i],65]
  end
  g_series = [129,125,124,123,122,117,111,104,65,31,30,21,13,11,10,9,1]
  for i = 1:17
    v_data[i] = v_c[65,g_series[i]]
  end
  writedlm("cavity u center data of $t_max Re$Re",u_data)
  writedlm("cavity v center data of $t_max Re$Re",v_data)

  curl = zeros((ny-1),(nx-1))
  x_coordc = zeros((ny-1),(nx-1))
  y_coordc = zeros((ny-1),(nx-1))

  for i = 1:nx-1
  for j = 1:ny-1
    x_coordc[j,i] = meshcell[j,i].cellcenter[1]
    y_coordc[j,i] = meshcell[j,i].cellcenter[2]
    u12 = (u[j+1,i+1]+u[j,i+1])/2
    u23 = (u[j+1,i+1]+u[j+1,i+2])/2
    u34 = (u[j+1,i+1]+u[j+2,i+1])/2
    u41 = (u[j+1,i+1]+u[j+1,i])/2

    v12 = (v[j+1,i+1]+v[j,i+1])/2
    v23 = (v[j+1,i+1]+v[j+1,i+2])/2
    v34 = (v[j+1,i+1]+v[j+2,i+1])/2
    v41 = (v[j+1,i+1]+v[j+1,i])/2

    curl[j,i] = (v23-v41)/meshcell[j,i].length12-(u34-u12)/meshcell[j,i].length23
   end
   end
    #writedlm("cavity curl data of $t_max Re$Re",curl)
end
function cavity_plot(t_max::Float64)
  x_coord = zeros((ny-1),(nx-1))
  y_coord = zeros((ny-1),(nx-1))
  for j = 1:ny-1
  for i = 1:nx-1
    x_coord[j,i] = meshcell[j,i].cellcenter[1]
    y_coord[j,i] = meshcell[j,i].cellcenter[2]
  end
  end
  contour(x_coord,y_coord,u[2:ny,2:nx],xlabel = "x",ylabel = "y")
  png(" cavity u contour tyler at $t_max Re $Re")
  contour(x_coord,y_coord,v[2:ny,2:nx],xlabel = "x",ylabel = "y")
  png(" cavity v contour tyler at $t_max Re $Re")
  contour(x_coord,y_coord,p[2:ny,2:nx],xlabel = "x",ylabel = "y",fill = true)
  png(" cavity p contour tyler at $t_max Re $Re")

  x_1d = zeros((nx-1))
  y_1d = zeros((ny-1))
  for i = 1:nx-1
   x_1d[i] = x_coord[1,i]
  end
  for j = 1:ny-1
    y_1d[j] = y_coord[j,1]
  end
  PyPlot.figure()
  PyPlot.xlabel("x")
  PyPlot.ylabel("y")
  PyPlot.title("Streamline")
  PyPlot.streamplot(x_1d,y_1d,u[2:ny,2:nx],v[2:ny,2:nx],density=0.7)
  PyPlot.savefig("Cavity Streamplot")
end
# ====== main loop function contains the whole three stage calculation ========
function main_loop(t_min::Float64,t_max::Float64,t_step::Float64,t_end::Int64)

  # pre define t-derivative slope
  Nu = zeros(Float64,ny-1,nx-1)
  Nv = zeros(Float64,ny-1,nx-1)
  Vu = zeros(Float64,ny-1,nx-1)
  Vv = zeros(Float64,ny-1,nx-1)

  Nu_1 = zeros(Float64,ny-1,nx-1)
  Nv_1 = zeros(Float64,ny-1,nx-1)
  Vu_1 = zeros(Float64,ny-1,nx-1)
  Vv_1 = zeros(Float64,ny-1,nx-1)

  for it = 1 : t_end

    t[it] = t_min + t_step*it
    print('\n',"==== time = $(t[it]), calculation start ====",'\n')

   # =========== velocity prediction for u*,v* =============================
    v_prediction(Nu,Nu_1,Nv,Nv_1,Vu,Vu_1,Vv,Vv_1)
    bcs(nx,ny,problem_type)
   # =========== Pressure calculation parts ================================
    p_calculation()
    bcs(nx,ny,problem_type)
   # =========== Pressure Correction =======================================
    p_correction(it)
    print("========= timestep finished =================== \n")
    bcs(nx,ny,problem_type)

   if t[it] == 0.5|| t[it]==5.0|| t[it]==10.0||t[it]== 15.0||t[it]==20.0||t[it]==25.0
     #tylergreen_result(t[it])
     #if t[it] == 1.0
     #tylergreen_plot(t[it])
     #end
     cavity_result(t[it])
     cavity_plot(t[it])
  end
  end
end
main_loop(t_min,t_max,t_step,t_end)


#  =================== plot and figure =============================


t_max = 5.0
#function  Data_extraction()
  u_g = zeros(ny,nx)
  v_g = zeros(ny,nx)
  p_g = zeros(ny,nx)
  u_g = readdlm("cavity u grid data of 5.0 Re100")
  v_g = readdlm("cavity v grid data of 5.0 Re100.txt")
  p_g = readdlm("cavity p grid data of $t_max Re$Re")



end
writedlm("p data of $t_max re$re",p)
writedlm("u data of $t_max re$re",u)
writedlm("v data of $t_max re$re",v)
=#
