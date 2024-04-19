import numpy as np
from scipy import special
import matplotlib.pyplot as plt
import sys, time

#
# Simulating Quantum Field Fluctuations in flat Minkowski space
#

def main1():
  # 
  # Quantum fluctuations using ground state wave functional 
  # in n-space with dimension d
  # 
  
  # Set mass for the field
  m = 1
  # Set size of the box
  L = 1
  # Set number of Fourier Coefficients
  n = 20  #20 - good middle
  
  
  
  """
  #
  # Plotting \phi_n Fourier Coefficient distribution.
  #
  """
  
  
  
  # Define omegas
  def omega_n(n):
    return np.sqrt(m**2+np.pi**2*n**2/L**2)
  
  # Define phi_n range of values
  phi_ns = np.linspace(-1.5,1.5,200)

  # For each n, have a probability distribution  
  def n_pdf(n):
    return np.sqrt(omega_n(n)/np.pi)*np.e**(-omega_n(n)*phi_ns**2) 
  
  
  ## Plot probability distributions of \phi_n for a few different scenarious
  

  # Fix m and L, vary n's
  plt.figure(1)
  plt.title(r"Probability Distributions of $\tilde\varphi_n$ (varying $n$)")
  plt.xlabel(r"$\tilde\varphi_n$")
  plt.plot(phi_ns, n_pdf(1), color="black", label="$n=1$")
  plt.plot(phi_ns, n_pdf(2), linestyle=":", color="red", label="$n=2$")
  plt.plot(phi_ns, n_pdf(3), linestyle="--", color="blue", label="$n=3$")
  plt.plot(phi_ns, n_pdf(4), linestyle="-.", color="green", label="$n=4$")
  plt.legend()
  plt.savefig(f"./n_pdf-m{m}-L{L}.pdf")
  
  # Fix n and m, vary L
  plt.figure(2)
  plt.title(r"Probability Distributions of $\tilde\varphi_n$ (varying $L$)")
  plt.xlabel(r"$\tilde\varphi_n$")
  
  L=1/2
  plt.plot(phi_ns, n_pdf(1), color="black", label=rf"$L={L}$")
  L=1
  plt.plot(phi_ns, n_pdf(1), linestyle=":", color="red", label=rf"$L={L}$")
  L=10
  plt.plot(phi_ns, n_pdf(1), linestyle="--", color="blue", label=rf"$L={L}$")
  L=100
  plt.plot(phi_ns, n_pdf(1), linestyle="-.", color="green", label=rf"$L={L}$")
  plt.legend()
  plt.savefig(f"./L_pdf-n{n}-m{m}.pdf")
  # Return to original L
  L=1
  
  # Fix n and L, vary m
  plt.figure(3)
  plt.title(r"Probability Distributions of $\tilde\varphi_n$ (varying $m$)")
  plt.xlabel(r"$\tilde\varphi_n$")
  
  m=0
  plt.plot(phi_ns, n_pdf(1), color="black", label=rf"$m={m}$")
  m=2
  plt.plot(phi_ns, n_pdf(1), linestyle=":", color="red", label=rf"$m={m}$")
  m=5
  plt.plot(phi_ns, n_pdf(1), linestyle="--", color="blue", label=rf"$m={m}$")
  m=10
  plt.plot(phi_ns, n_pdf(1), linestyle="-.", color="green", label=rf"$m={m}$")
  plt.legend()
  plt.savefig(f"./m_pdf-n{n}-L{L}.pdf")
  
  
  
  """
  #
  # Plotting \phi(x) in multiple dimensions using 
  # our distribution of Fourier Coefficients \phi_n
  #
  """
  
  
  
  # The Inverse Cumulative distribution function.
  # Used for sampling our distributed Fourier coefficients \phi_n
  def n_cdfinv(n,y):
    return special.erfinv(2*y-1)/np.sqrt(omega_n(n))
    
  # Check if sampling works!!!
  x=[]
  n=2
  for i in range(200):
    dd = np.random.random()
    x.append(n_cdfinv(n,dd))
  x=sorted(x)  
  
  plt.figure(2)
  plt.title(rf"Random Samplings of $\tilde\varphi_n$ for n={n}")
  plt.xlabel(r"$\tilde\varphi_n$")
  plt.hist(x,label=r"$\tilde\varphi_n$ samples", density=True)
  plt.plot(phi_ns, n_pdf(n),label=r"$\tilde\varphi_n$ distribution")
  plt.legend()
  plt.savefig("./n_samples.pdf")
  
  # Works!!!
  

  #
  # Plot \phi(x) in 2D and 3D
  #

   
  def phi_2d(x,y,t):
    # D=2
    s=0 # Init sum
    # Create range of n values for sum
    ns = np.linspace(1,n,n,dtype=int)
    # Iterate through each n, taking into account the dimension.
    for i in range(n): 
      n_1 = ns[i]
      for j in range(n): 
        n_2 = ns[j]
        # Caluclate n = \sqrt{n_1^2 + n_2^2}
        nn = np.sqrt(n_1**2+n_2**2)
        # Sample values of the Fourier Coeffecient given n
        dd = np.random.random()
        phi_n = n_cdfinv(nn,dd)
        # Calculate the sum for phi(x)
        s = s + (2/L)**(d/2)*phi_n*np.sin(np.pi/L*(n_1*x + n_2*y)-omega_n(nn)*t)
    
    return s
  
  
  def phi_3d(x,y,t): # With z=0
    # D=3
    s=0 # Init sum
    # Create range of n values for sum
    ns = np.linspace(1,n,n,dtype=int)
    # Iterate through each n, taking into account the dimension.
    for i in range(n): 
      n_1 = ns[i]
      for j in range(n): 
        n_2 = ns[j]
        for k in range(n): 
          n_3 = ns[k]
          # Caluclate n = \sqrt{n_1^2 + n_2^2 + n_3^3}
          nn = np.sqrt(n_1**2+n_2**2+n_3**2)
          # Sample values of the Fourier Coeffecient given n
          dd = np.random.random()
          phi_n = n_cdfinv(nn,dd)
          # Calculate the sum for phi(x)
          s = s + (2/L)**(d/2)*phi_n*np.sin(np.pi/L*(n_1*x + n_2*y + n_3*0)-omega_n(nn)*t)
            
    return s
  
  
  ## Plotting 
  
  
  # Set mass for the field
  m = 1
  # Set size of the box
  L = 1
  # Set number of Fourier Coefficients
  n = 20  #20 - good middle
  
  
  
  # Parameters
  x = np.linspace(-1/2,1/2,100)
  t = 0 # Work with single time slice for now
  d = 2

  # Plot in 2D (t=0)
  xs, ys = np.meshgrid(x, x)
  phis_2d = phi_2d(xs,ys,t) 
  
  plt.figure(4)
  plt.title(rf"$\varphi(\vec x)$ fluctuations in {d}D")
  plt.xlabel(r"$x$")
  plt.ylabel(r"$y$")
  heatmap_2D = plt.imshow(phis_2d, extent =[-1/2, 1/2, -1/2, 1/2], cmap="seismic")
  plt.colorbar(heatmap_2D)
  plt.savefig(f"./x{d}D_m{m}-L{L}-n{n}.pdf")
  plt.close()
  
  
  # Plot in 3D (t=0)
  
  # Parameters
  d  = 3
  
  xs, ys = np.meshgrid(x, x) #zs=0
  phis_3d = phi_3d(xs,ys,t) 
  
  plt.figure(5)
  plt.title(rf"$\varphi(\vec x)$ fluctuations in {d}D")
  plt.xlabel(r"$x$")
  plt.ylabel(r"$y$")
  heatmap_3D = plt.imshow(phis_3d, extent =[-1/2, 1/2, -1/2, 1/2], cmap="seismic")
  plt.colorbar(heatmap_3D)
  plt.savefig(f"./x{d}D_m{m}-L{L}-n{n}.pdf")
  plt.close()
  
  
  #
  # Plotting time varying field amplitudes movies! # Ignore this in submission for now?
  # 
  
  
  # NOTE!! That this is essentially the same, in the case of Minkowski spacetime, as recalculating \phi(\bold x,t) at different time intervals, as this would emulate the behaviour of \phi_n(t) as a time dependent fluctuating fourier coefficient?
  #
  
  ts = np.linspace(0,10)
  phis_2dt_0 = phi_2d(xs,ys,ts[0]) # initialize 
  d = 2
  
  # Plot in 2D (varying t)
  # Plot initial map
  fig = plt.figure(6)
  fig.suptitle(rf"$\varphi(\vec x)$ fluctuations in {d}D")
  fig.supxlabel(r"$x$")
  fig.supylabel(r"$y$")
  ax = fig.add_subplot(111)
  heatmap2Dt = ax.imshow(phis_2dt_0, extent=[-1/2, 1/2, -1/2, 1/2], cmap="seismic")
  # Show heat map
  fig.colorbar(heatmap2Dt)
  plt.show(block=False)
  
  # Update time and re-plot
  for i in range(1,len(ts)):
    t = ts[i]
    # time sleep maybe if too fast
#    time.sleep(0.01)
    # recalculate for at new time step
    phis_2dt_step = phi_2d(xs,ys,t)
    #replace imshow image contents
    heatmap2Dt.set_array(phis_2dt_step)
    # redraw the new plot
    fig.canvas.draw()
    fig.canvas.flush_events()
    
  
  
if __name__ == "__main__": main1()
