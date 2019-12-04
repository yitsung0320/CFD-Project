# MAE 560 CFD
# Time Derivation Method
using LinearAlgebra

function LUsolver(A::Array{Float64,2},B::Array{Float64,1})

      d = size(A,1)
      y = Array{Float64,1}(undef,d)
      x = Array{Float64,1}(undef,d)
      L = lu(A).L
      U = lu(A).U

      # Ly = B
      y[1] = B[1]/L[1,1]
      for i = 2 : d
         yl = 0.0
         for j = 1 : i-1
            yl = yl+y[j]*L[i,j]
         end
        y[i] = (B[i]-yl)/L[i,i]
      end

      #Ux = y
      x[d] = y[d]/U[d,d]
      for i = 2 : d
         ux = 0.0
         for j = 1 : i-1
            ux = ux + x[d+1-j]*U[d+1-i,d+1-j]
         end
       x[d+1-i] = (y[d+1-i]-ux)/U[d+1-i,d+1-i]
      end
    return x
end
# =============== Runge Kutta 3 method ===================================
function Runge_Kutta3(t_step::Float64,f::Array{Float64,1},F_x::Array{Float64,1},
                      num_cells::Int64)
  #t_step = time step
  #F_x = Flux term

         k1 = Array{Float64,1}(undef,num_cells+2)
         k2 = Array{Float64,1}(undef,num_cells+2)
         k3 = Array{Float64,1}(undef,num_cells+2)
         f_t= Array{Float64,1}(undef,num_cells+2)

        for i = 1: num_cells+2
         k1[i] = F_x[i]
         k2[i] = F_x[i] + t_step/2*k1[i]
         k3[i] = F_x[i] - 1*t_step*k1[i]+2*t_step*k2[i]
         f_t[i] = f[i] + t_step/6*(k1[i]+4*k2[i]+k3[i])
        end
     return f_t
end
# =============== Crank Nicolson method ==================================
function Crank_Nicolson(t_step::Float64,f::Array{Float64,1},fn::Array{Float64,2}
                        ,fn1::Array{Float64,2},fd::Array{Float64,1},num_cells::Int64)
        # f = cell solution at t = n
        # fd = dump vector usually contains bcS
        A = -1.0*fn1*t_step/2.0
        B = fn*t_step/2.0
        fd = fd*t_step/2

        for i = 1:num_cells
         A[i,i] = A[i,i]+1
         B[i,i] = B[i,i]+1
        end
        B = B*f[2:num_cells+1] + fd

      # now ready to solve Ax = B
      # output f at  t = n+1
      f = LUsolver(A,B)

end
# =============== Adams Bashforth 2 method ================================
function Adams_Bashforth2(t_step::Float64,f::Array{Float64,2},
                          df_n::Array{Float64,2},df_n1::Array{Float64,2})
 f[:]= f[:] + t_step*(3/2*df_n[:]-1/2*df_n1[:])
 return f

end
