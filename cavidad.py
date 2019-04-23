"""
Programa que resuelve el problema de Navier-Stokes estacionario dado por

-\nabla \cdot \sigma(u,p) + (\nabla u)u = f   en   \Omega
                         \nabla \cdot u = 0   en   \Omega     
                                      u = u_D en   \partial\Omega

donde:

 \sigma(u,p) :=  2\nu(u) \varepsilon(u) - pI


Rodolfo Araya D
rodolfo.araya@udec.cl


Benchmark de la cavidad

"""
# Importando los modulos de FeniCS
from dolfin import *

# Opciones de optimizacion para C++ 
parameters["form_compiler"]["cpp_optimize"]       = True
parameters['form_compiler']['cpp_optimize_flags'] = '-O3 -ffast-math -march=native'
parameters["form_compiler"]["optimize"]           = True

#dimension del espacio
dim = 2

# Generando la malla donde se haran los calculos
#
nx = ny = 300 # si agrandas estos valores la malla sera mas fina y
              # la solucion mas precisa
mesh = UnitSquareMesh(nx, ny)

#Referencias de los lados para imponer las condiciones de frontera
def upper(x, on_boundary):
    return near(x[1], 1.0)

def wall(x, on_boundary):
    if on_boundary:
        if near(x[1], 1.0):
            return False
        else:
            return True
    else:
        return False
        
# Datos del problema
Re = 200.0                # Reynolds
nu = Constant(1.0/Re)     # Viscosidad 
f  = Constant((0.0, 0.0)) # lado derecho

# Espacios discretos
degree = 2
elemento = mesh.ufl_cell() # malla triangular
H    = VectorElement("Lagrange", elemento, degree)    # Polinomios de grado 2 para aproximar la velocidad
Q    = FiniteElement("Lagrange", elemento, degree-1)  # Polinomios de grado 2 para aproximar la presion
HxQ  = FunctionSpace(mesh, H*Q)

# Funciones test e incognitas
(v, q) = TestFunctions(HxQ)
up     = Function (HxQ)
(u,p)  = split(up) #separa up 
 
# Condiciones de frontera de tipo Dirichlet
u_upper = Constant((1.0,0.0))
u_wall  = Constant((0.0,0.0))

bc_upper = DirichletBC(HxQ.sub(0), u_upper, upper)
bc_wall  = DirichletBC(HxQ.sub(0), u_wall,  wall)

dirichlet=[bc_wall,bc_upper]

# Tensor de deformaciones:
def epsilon(r):
    return sym(grad(r))

# Funcional al cual se le debe encontrar el cero:
F = (2*nu*inner(epsilon(u), epsilon(v)) +  inner(grad(u)*u, v) - div(v)*p + q*div(u))*dx -inner(f,v)*dx

# Calculo de la solucion

# Primer metodo: Todo se hace de forma automatica
solve(F == 0, up, dirichlet)

#Segundo metodo: Un poco mas de control de lo que hace Fenics (metodo mas robusto)
#J = derivative(F, up )
#
#problem = NonlinearVariationalProblem(F, up , dirichlet, J)
#solver  = NonlinearVariationalSolver(problem)
#
## Parametros para resolver el problema no lineal:
#prm = solver.parameters                                 # se acorta el nombre para hacerlo mas facil de escribir
#prm["nonlinear_solver"]                  = "snes"       # Scalable Nonlinear Equations Solvers (SNES) del sistema PETSc
#prm["snes_solver"]["line_search"]        = "bt"         # Se usa un proceso de continuacion
#prm["snes_solver"]["linear_solver"]      = "mumps"      # Metodo para resolver los problemas lineales que aparecen
##prm["snes_solver"]["linear_solver"]     = "gmres"      # Metodo para resolver los problemas lineales que aparecen
##prm["snes_solver"]["preconditioner"]    = "hypre_amg"  # Precondicionador para resolver los problemas lineales que aparecen
##prm["snes_solver"]["krylov_solver"]["nonzero_initial_guess"] = False  # vector de partida del metodo iterativo
#prm["snes_solver"]["report"]             = True         # Informacion por pantalla
#
#solver.solve()

#recuperando las variables originales
(u, p) = up.split()

#generando los archivos de salida
u.rename('velocidad', 'velocidad del fluido')
p.rename('presion', 'presion del fluido')

velo = File('resultados/velocidad.pvd')
pres = File('resultados/presion.pvd')

velo << u
pres << p

import matplotlib.pyplot as plt
plot(u)
#plot(p)
plt.show()
