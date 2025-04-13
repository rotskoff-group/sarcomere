import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Wedge
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import math
import joblib
from joblib import Parallel, delayed
import os
import sys
import argparse

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

def plot_filaments(center,thetas,l,Lx,Ly,ax,color='k',color_spectrum=None,transparency=None,**kwargs):
    #check if l is a scalar or an array
    if np.isscalar(l):
        l = np.ones(center.shape[0])*l
    if color_spectrum is not None:
        color_spectrum = np.sqrt(color_spectrum)
        # Use a colormap to get color values from the color spectrum
        colormap = plt.cm.Blues  # Adjust colormap if desired
        custom_colormap = mcolors.LinearSegmentedColormap.from_list(
        "custom_blues", colormap(np.linspace(0.3, 1, 256))
        )  # Start from a light blue (0.3) instead of white (0.0)
        # Normalize color spectrum values between 0 and 1
        norm = plt.Normalize(0, 1)
        # Create ScalarMappable for color bar based on color spectrum
        sm = plt.cm.ScalarMappable(cmap=custom_colormap, norm=norm)
        sm.set_array([])  # Dummy array for the color bar
        
    for i in range(center.shape[0]):
        if l[i]<0.01:
            continue
        # Determine color based on the color spectrum
        if color_spectrum is not None:
            color = custom_colormap(norm(color_spectrum[i]))
        
        if transparency is not None:
            transparency[i] = np.clip(transparency[i],0,1)
            kwargs['alpha'] = 0.3 + 0.7 * transparency[i]
            
        # transparency is 0.3 where color_spectrum is 0 and 1 elsewhere
        if color_spectrum is not None and transparency is None:
            alpha = 0.3 + 0.7 * (color_spectrum[i] > 0)
            kwargs['alpha'] = alpha

        x = center[i,0]
        y = center[i,1]
        theta = thetas[i]
        l_i = l[i]
        delta_x = 0.5*l_i*np.cos(theta)
        delta_y = 0.5*l_i*np.sin(theta)
        ax.arrow(x-delta_x, y-delta_y, 2*delta_x, 2*delta_y,
                    head_width=0.1, head_length=0.1, fc=color, ec=color,**kwargs)
        #plot the periodic images
        if x+np.abs(delta_x)>Lx/2:
            ax.arrow(x-delta_x-Lx, y-delta_y, 2*delta_x, 2*delta_y,
                        head_width=0.1, head_length=0.1, fc=color, ec=color,**kwargs)
            if y+np.abs(delta_y)>Ly/2:
                ax.arrow(x-delta_x-Lx, y-delta_y-Ly, 2*delta_x, 2*delta_y,
                            head_width=0.1, head_length=0.1, fc=color, ec=color,**kwargs)
            elif y-np.abs(delta_y)<-Ly/2:
                ax.arrow(x-delta_x-Lx, y-delta_y+Ly, 2*delta_x, 2*delta_y,
                            head_width=0.1, head_length=0.1, fc=color, ec=color,**kwargs)
        elif x-np.abs(delta_x)<-Lx/2:
            ax.arrow(x-delta_x+Lx, y-delta_y, 2*delta_x, 2*delta_y,
                        head_width=0.1, head_length=0.1, fc=color, ec=color,**kwargs)
            if y+np.abs(delta_y)>Ly/2:
                ax.arrow(x-delta_x+Lx, y-delta_y-Ly, 2*delta_x, 2*delta_y,
                            head_width=0.1, head_length=0.1, fc=color, ec=color,**kwargs)
            elif y-np.abs(delta_y)<-Ly/2:
                ax.arrow(x-delta_x+Lx, y-delta_y+Ly, 2*delta_x, 2*delta_y,
                            head_width=0.1, head_length=0.1, fc=color, ec=color,**kwargs)    
        if y+np.abs(delta_y)>Ly/2:
            ax.arrow(x-delta_x, y-delta_y-Ly, 2*delta_x, 2*delta_y,
                        head_width=0.1, head_length=0.1, fc=color, ec=color,**kwargs)
            if x+np.abs(delta_x)>Lx/2:
                ax.arrow(x-delta_x-Lx, y-delta_y-Ly, 2*delta_x, 2*delta_y,
                            head_width=0.1, head_length=0.1, fc=color, ec=color,**kwargs)
            elif x-np.abs(delta_x)<-Lx/2:
                ax.arrow(x-delta_x+Lx, y-delta_y-Ly, 2*delta_x, 2*delta_y,
                            head_width=0.1, head_length=0.1, fc=color, ec=color,**kwargs)
        elif y-np.abs(delta_y)<-Ly/2:
            ax.arrow(x-delta_x, y-delta_y+Ly, 2*delta_x, 2*delta_y,
                        head_width=0.1, head_length=0.1, fc=color, ec=color,**kwargs)
            if x+np.abs(delta_x)>Lx/2:
                ax.arrow(x-delta_x-Lx, y-delta_y+Ly, 2*delta_x, 2*delta_y,
                            head_width=0.1, head_length=0.1, fc=color, ec=color,**kwargs)
            elif x-np.abs(delta_x)<-Lx/2:
                ax.arrow(x-delta_x+Lx, y-delta_y+Ly, 2*delta_x, 2*delta_y,
                            head_width=0.1, head_length=0.1, fc=color, ec=color,**kwargs)
    # Add color bar
    if color_spectrum is not None:
        cbar = plt.colorbar(sm, ax=ax, orientation='vertical', fraction=0.02, pad=0.04)
        cbar.set_label('Catch bond strength')

def plot_myosin(center,thetas,l,Lx,Ly,radius,ax,color='k'):
    for i in range(center.shape[0]):
        x = center[i,0]
        y = center[i,1]
        theta = thetas[i]
        delta_x = l/2 * np.cos(theta) - radius * np.sin(theta)
        delta_y = l/2 * np.sin(theta) + radius * np.cos(theta)
        #plot a rectangle to represent the myosin
        ax.add_patch(Rectangle((x-delta_x,y-delta_y),width=l,height=2*radius,angle=theta*180/math.pi,fill=True,facecolor=color,alpha=0.5))
        #plot the periodic images
        ax.add_patch(Rectangle((x-delta_x-Lx,y-delta_y),width=l,height=2*radius,angle=theta*180/math.pi,fill=True,facecolor=color,alpha=0.5))
        ax.add_patch(Rectangle((x-delta_x+Lx,y-delta_y),width=l,height=2*radius,angle=theta*180/math.pi,fill=True,facecolor=color,alpha=0.5))
        ax.add_patch(Rectangle((x-delta_x,y-delta_y-Ly),width=l,height=2*radius,angle=theta*180/math.pi,fill=True,facecolor=color,alpha=0.5))
        ax.add_patch(Rectangle((x-delta_x,y-delta_y+Ly),width=l,height=2*radius,angle=theta*180/math.pi,fill=True,facecolor=color,alpha=0.5))
        ax.add_patch(Rectangle((x-delta_x-Lx,y-delta_y-Ly),width=l,height=2*radius,angle=theta*180/math.pi,fill=True,facecolor=color,alpha=0.5))
        ax.add_patch(Rectangle((x-delta_x+Lx,y-delta_y-Ly),width=l,height=2*radius,angle=theta*180/math.pi,fill=True,facecolor=color,alpha=0.5))
        ax.add_patch(Rectangle((x-delta_x-Lx,y-delta_y+Ly),width=l,height=2*radius,angle=theta*180/math.pi,fill=True,facecolor=color,alpha=0.5))
        ax.add_patch(Rectangle((x-delta_x+Lx,y-delta_y+Ly),width=l,height=2*radius,angle=theta*180/math.pi,fill=True,facecolor=color,alpha=0.5))
        #draw spheres at the ends of the myosin
        end_pts = np.array([[x-l/2*np.cos(theta),y-l/2*np.sin(theta)],[x+l/2*np.cos(theta),y+l/2*np.sin(theta)]])
        theta_perp = np.array([[theta*180/math.pi+90,theta*180/math.pi-90],[theta*180/math.pi-90,theta*180/math.pi+90]])
        for pt,theta_p in zip(end_pts,theta_perp):
            #draw a circle with no edge to represent the myosin head
            ax.add_artist(Wedge((pt[0], pt[1]), radius, color=color, alpha=0.5, linewidth=0,theta1=theta_p[0],theta2=theta_p[1]))
            ax.add_artist(Wedge((pt[0]-Lx, pt[1]), radius, color=color, alpha=0.5, linewidth=0,theta1=theta_p[0],theta2=theta_p[1]))
            ax.add_artist(Wedge((pt[0]+Lx, pt[1]), radius, color=color, alpha=0.5, linewidth=0,theta1=theta_p[0],theta2=theta_p[1]))
            ax.add_artist(Wedge((pt[0], pt[1]-Ly), radius, color=color, alpha=0.5, linewidth=0,theta1=theta_p[0],theta2=theta_p[1]))
            ax.add_artist(Wedge((pt[0], pt[1]+Ly), radius, color=color, alpha=0.5, linewidth=0,theta1=theta_p[0],theta2=theta_p[1]))
            ax.add_artist(Wedge((pt[0]-Lx, pt[1]-Ly), radius, color=color, alpha=0.5, linewidth=0,theta1=theta_p[0],theta2=theta_p[1]))
            ax.add_artist(Wedge((pt[0]+Lx, pt[1]-Ly), radius, color=color, alpha=0.5, linewidth=0,theta1=theta_p[0],theta2=theta_p[1]))
            ax.add_artist(Wedge((pt[0]-Lx, pt[1]+Ly), radius, color=color, alpha=0.5, linewidth=0,theta1=theta_p[0],theta2=theta_p[1]))
            ax.add_artist(Wedge((pt[0]+Lx, pt[1]+Ly), radius, color=color, alpha=0.5, linewidth=0,theta1=theta_p[0],theta2=theta_p[1]))

           

def plot_system(frame,data,myosin_radius,actin_length,myosin_length, Lx,Ly):

        fig, ax = plt.subplots(dpi=500)
        actin_center = data["/actin/center"][frame]
        #color is based on cb_flag
        cb_strength = data["/actin/cb_strength"][frame].flatten()
        f_load = data["/actin/f_load"][frame].flatten()
        cb_strength = cb_strength*f_load
        binding_ratio = data["/actin/partial_binding_ratio"][frame].flatten()
        binding_ratio = np.clip(binding_ratio*3,0,1)
        crosslinking_ratio = data["/actin/crosslink_ratio"][frame].flatten()
        crosslinking_ratio = np.clip(crosslinking_ratio*3,0,1)
        actin_velocity = data["/actin/velocity"][frame]
        velocity_theta = np.arctan2(actin_velocity[:,1],actin_velocity[:,0])
        velocity_mag = np.linalg.norm(actin_velocity,axis=-1)/5
        plot_filaments(actin_center, data["/actin/theta"][frame][:,0],actin_length,Lx=Lx, Ly=Ly,ax=ax, color_spectrum=cb_strength)
        # actin_force = data["/actin/force"][frame]/10
        # force_theta = np.arctan2(actin_force[:,1],actin_force[:,0])
        # force_mag = np.linalg.norm(actin_force,axis=-1)
        # force_centers = actin_center + actin_force/2
        # plot_filaments(force_centers, force_theta,force_mag,Lx=Lx, Ly=Ly,ax=ax, color='grey',linestyle='dashed')
        # velocity_centers = actin_center + actin_velocity/2
        # plot_filaments(velocity_centers, velocity_theta,velocity_mag,Lx=Lx, Ly=Ly,ax=ax, color='pink',linestyle='dashed',alpha=0.8)

        actin_angular_force = data["/actin/angular_force"][frame][:,0]
        angular_theta = ((actin_angular_force>0)-0.5)*np.pi
        #plot_filaments(actin_center,angular_theta,np.abs(actin_angular_force),Lx=Lx, Ly=Ly,ax=ax, color='pink',linestyle='dotted')
        plot_myosin(data["/myosin/center"][frame], data["/myosin/theta"][frame][:,0],myosin_length,Lx=Lx, Ly=Ly,radius=myosin_radius,ax=ax,color='C1')

        myosin_force = data["/myosin/force"][frame]
        force_theta = np.arctan2(myosin_force[:,1],myosin_force[:,0])
        force_mag = np.linalg.norm(myosin_force,axis=-1)
        force_centers = data["/myosin/center"][frame] + myosin_force/2
        plot_filaments(force_centers, force_theta,force_mag,Lx=Lx, Ly=Ly,ax=ax,color='pink',linestyle='dashed')

        myosin_angular_force = data["/myosin/angular_force"][frame][:,0]
        myosin_theta = data["/myosin/theta"][frame][:,0]
        angular_theta = (((myosin_angular_force*myosin_theta)>0)-0.5)*np.pi
        #plot_filaments(data["/myosin/center"][frame],angular_theta,np.abs(myosin_angular_force),Lx=Lx, Ly=Ly,ax=ax, color='pink',linestyle='dotted')
        ax.set_xlim(-Lx/2, Lx/2)
        ax.set_ylim(-Ly/2, Ly/2)
        ax.set_aspect("equal")
        return fig,ax

def plot(ind, nframes, nworkers,**kwargs):
    frame_start = int(ind*(nframes)/nworkers)
    frame_end = max(int((ind+1)*(nframes)/nworkers),1)
    for frame in range(frame_start,frame_end):
        fig, ax = plot_system(frame=frame,**kwargs)
        plt.savefig(file_format.format(frame_dir,frame))
        plt.close(fig)


def hdf5_to_dict(hdf5_file):
    data_dict = {}
    def _hdf5_to_dict_recursive(group, group_name):
        for key in group.keys():
            if isinstance(group[key], h5py.Dataset):
                data_dict[f"{group_name}/{key}"] = np.array(group[key])
            elif isinstance(group[key], h5py.Group):
                _hdf5_to_dict_recursive(group[key], f"{group_name}/{key}")
    
    _hdf5_to_dict_recursive(hdf5_file, "")
    return data_dict



#write an argument parser
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--filename", type=str, default='data/traj.h5')
    parser.add_argument("--frame_dir", type=str, default="frames")
    parser.add_argument("--Lx", type=float, default=10)
    parser.add_argument("--Ly", type=float, default=10)
    parser.add_argument("--myosin_radius", type=float, default=0.2)
    parser.add_argument("--actin_length", type=float, default=1)
    parser.add_argument("--myosin_length", type=float, default=1.5)
    parser.add_argument("--print_frame",type=int,default=1000)
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    filename = args.filename
    frame_dir = args.frame_dir
    Lx = args.Lx
    Ly = args.Ly
    myosin_radius = args.myosin_radius
    actin_length = args.actin_length
    myosin_length = args.myosin_length

    traj = h5py.File(filename, 'r')
    data = hdf5_to_dict(traj)
    # energy = data["/energy/total_energy"]
    traj = h5py.File(filename, 'r')
    data = hdf5_to_dict(traj)
    nframes = data["/actin/center"].shape[0]
    print(f"Number of frames: {nframes}")
    if args.print_frame<nframes:
        actin_center = data["/actin/center"][args.print_frame]
        cb_strength = data["/actin/cb_strength"][args.print_frame]
        for i in range(actin_center.shape[0]):
            print(f"Actin filament {i}: center: {actin_center[i]}")
            if cb_strength[i]>0:
                print(f"Actin filament {i}: crossbridge strength: {cb_strength[i]}")
        myosin_center = data["/myosin/center"][args.print_frame]
        for i in range(myosin_center.shape[0]):
            print(f"Myosin filament {i}: center: {myosin_center[i]}")
            print(f"Myosin filament {i}: force: {data['/myosin/force'][args.print_frame][i]}")

    cpu_workers = joblib.cpu_count()

    n_digits = len(str(nframes))
    file_format = "{}/frame_{:0"+str(n_digits)+"d}.png"
    if not os.path.exists(frame_dir):
        os.mkdir(frame_dir)
    Parallel(n_jobs=cpu_workers)(delayed(plot)(i, nframes, cpu_workers,
        data=data, myosin_radius=myosin_radius, 
        actin_length=actin_length, myosin_length=myosin_length, Lx=Lx, Ly=Ly) for i in range(cpu_workers))
    
