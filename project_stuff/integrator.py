# Graduate Computational Physics
# Gravity on a Mesh Project
# Integrator file

import numpy as np
import densities
import solver

# function to run one time step of the integrator
def integrate(positions, velocities, density,time_step,grid_size=32):
    '''
    Integrates the particles' positions and velocities forward by one time step
    using the gradient of the gravitational potential and the Verlet integration method

    :param positions: array size (num_particles, 3)
                    array of position values (x, y, z) for each particle
    :param velocities: array size (num_particles, 3)
                    array of velocity values (vx, vy, vz) for each particle
    :param potential: array size (32, 32, 32)
                    array of potential values at each point in 32x32x32 space
    :param time_step: float
                    size of whole time step
    :return: new_positions:
                    returns array of positions values (x, y, z) after one time step
    :return: new_velocities:
                    returns array of positions values (x, y, z) after one time step
    :return: new_potential:
                    returns array of potential values
    '''

    potential = solver.solve_poisson(density)

    # get values of force (3 components) on each particle at the current positions
    F_current = grav_force(potential, positions)

    # calculate v(t + half step)
    v_half = velocities + 0.5*time_step*F_current

    #calculate x(t + step)
    new_positions = positions + time_step*v_half

    #enforcing toroidal boundary conditions
    new_positions = new_positions % 32

    # calculate updated density with new positions
    new_density = densities.cic_density(new_positions)


    # update potential
    new_potential = solver.solve_poisson(new_density)




    # calculate new values of force from the new potential
    new_F = grav_force(new_potential, new_positions)
    # calculate v(t + step)
    new_velocities = v_half + time_step*new_F

    # return position, velocity, new_potential after a whole time step
    return new_positions, new_velocities, new_potential

def grav_force(potential, positions,grid_size=32,cell_size=1):
    '''
    Calculates the gravitational force on each particle using CIC interpolation (to keep consistency with density asssignment)
    Takes the gradient of the potential and returns the value of the potential gradient (force) nearest
    to this position (approximation since our potential is a discrete function)

    :param potential: array size (32, 32, 32)
                    array of potential values at each point in 32x32x32 space
    :param positions: array size (num_particles, 3)
                    array of position values for each particle (x, y, z)
    :return: F: array size (num_particles, 3)
                    returns array of force on each particle in each direction (Fx, Fy, Fz)
    '''

    # acceleration is the negative gradient of the potential (g_x,g_y,g_z)
    # assuming mass =1
    g = np.gradient(potential)
    g = [(-1* grad) for grad in g]
    #stacking so that the x,y,z components of the force can be accessed through the first index
    # [0,:,:,:] gives x component of force for each cell in the grid
    g = np.stack(g, axis=0)

    #(num_particles x 3) array to store acceleration for each particle
    #each row stores (g_x,g_y,g_z)
    g_p = np.zeros(np.shape(positions))

    #index for each particle
    p_num=0

    # populating acceleration array according to CIC interpolation
    for x_p in positions:

        # indices of the cell containing the particle (x_p,y_p,z_p)
        indices = np.int32(np.floor(x_p)) % grid_size

        i, j, k = indices[0], indices[1], indices[2]

        # indices of the next cell in each direction (i+1,j+1,k+1) with boundary conditions enforced
        ii, jj, kk = (i + 1) % grid_size, (j + 1) % grid_size, (k + 1) % grid_size

        # centers of the cells containing the particle (x_c,y_c,z_c)
        x_c = np.array([i + (cell_size / 2), j + (cell_size / 2), k + (cell_size / 2)])

        # distance of the particle from the cell center (|x_p-x_c|,|y_p-y_c|.|z_p-z_c|)
        #effectively gives us the fraction of the particle in neighboring cells
        d = np.abs(x_p - x_c)

        #fraction of particle in the parent cells
        t = np.abs(cell_size - d)

        # implementation of CIC interpolation for acceleration
        # (interpolated from the cells to which particle contributed density)
        # ensures absence of artificial self-forces and the third Newton’s law ("Writing a PM code")
        for index in range(3):
            g_p[p_num,index] = ((g[index,i,j,k]*(t[0]*t[1]*t[2]) + g[index,ii,j,k]*(d[0]*t[1]*t[2])+g[index,i,jj,k]*(t[0]*d[1]*t[2]) +g[index,ii,jj,k]*(d[0]*d[1]*t[2])) +g[index,i,j,kk]*(t[0]*t[1]*d[2])+g[index,ii,j,kk]*(d[0]*t[1]*d[2])  +g[index,i,jj,kk]*(t[0]*d[1]*d[2])+g[index,ii,jj,kk]*(d[0]*d[1]*d[2]))

        p_num += 1

    return g_p

# placeholder for laplace equation solver (will come from different source file)
#def solve_laplace(density):
#    potential = density
#    return potential

# old stuff from grav_force
'''
for i in range(np.size(xpos)):
        # Fx = dphi/dx
        F[i][0] = dphi[0][xpos[i], ypos[i], zpos[i]]
        # Fy = dphi/dy
        F[i][1] = dphi[1][xpos[i], ypos[i], zpos[i]]
        # Fz = dphi/dz
        F[i][2] = dphi[2][xpos[i], ypos[i], zpos[i]]
'''


