# from simple_test import run

# def u0(x):
#     #return 1.0*(x <= 0.5) + 2.0 * (x>0.5)
#     return 1 + np.sin(2*np.pi*x)

# N = 25
# x = np.linspace(0,1,N)
# u = u0(x)
# gw = 2
# T = 1.
# M = 100
# dt = T/M

# dx = x[1]-x[0]
# uGhost = np.zeros(N+2*gw)

# xGhost = np.zeros(N+2*gw)
# xGhost[0] = x[-1]
# xGhost[1:N+1] = x
# xGhost[N+1] = x[0]

# tend = np.zeros(N)

# # Input:
# # ui: u[i] for i the center point of the stencil
# # uj: u[j] for a single j, or array of values for all j in the stencil
# # xi: x[i] for the center point of the stencil
# # xj: x[j] (like uj). If array, it must be length(xj)==length(uj)
# # Output:
# # If uj is a number, returns g_i(x_j) 
# # If uj is an array, returns [g_i(x_j) for every j in the stencil],
# #                             which can be unpacked with *
# def g(ui,uj,xi,xj):
#     ### linear
# #    return uj - ui*np.exp(xj-xi)   

#     ### Burgers
#    return uj*uj*0.5 - 0.5*ui*ui*np.exp(2*(xj-xi))  


# for n in range(M):
#     print u
#     plt.plot(x,u)
#     plt.show()
#     uGhost[:gw] = u[-gw:]
#     uGhost[gw:N+gw] = u
#     uGhost[-gw:] = u[:gw]

#     for i in range(N):
#         iOff = i+gw
#         u_st = uGhost[iOff-gw:iOff+gw+1] # u at the stencil for ui
#         x_st = xGhost[iOff-gw:iOff+gw+1] # x at the stencil for ui
#         print (i, x_st)
#         Gr = wr.weno5_rec(*g(u[i],u_st,x[i],x_st))
#         Gl = wr.weno5_rec(*g(u[i],u_st[::-1],x[i],x_st[::-1]))
#         tend[i] = -(Gr - Gl)/dx;
#     print tend
#     u = u + dt*tend

# plt.plot(x,u)
# plt.show()

