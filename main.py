
# import pygame as pg
import numpy as np
import json, sys, os
import matplotlib.pyplot as plt

from raytracing import Sphere, Scene, Light


def plot_canvas(canvas):
    plt.figure(figsize=(8,8))
    plt.imshow(canvas,)
    plt.axis("off")
    plt.show()

def main():
    Cw = 400
    Ch = 400
    Vw = 1
    Vh = 1
    f  = 1

    spheres = [Sphere(radius=0.1,specular_coeff=1000), 
               Sphere(position=(2,0,7), color=(0,0,255),radius=1,specular_coeff=800),
               Sphere(position=(2,4,10), color=(112, 41, 99),radius=2.,specular_coeff=1000),
               Sphere(position=(2,-4,10), color=(49, 203, 164),radius=2.,specular_coeff=1000), 
               Sphere(color=(255,255,0), radius=5000, position=(5003,0,0))]
    lights  = [Light(intensity=0.2), Light(intensity=0.6, position=(0,0,0)), Light(intensity=0.2,direction=(-1,-1,1))]
    scene = Scene(Cw, Ch, Vw, Vh, f, spheres, lights, BG_COLOR=(255,255,255))

    O = np.zeros(3)
    for x in np.arange(-Cw//2, Cw//2):
        for y in np.arange(-Ch//2, Ch//2):
            xa,ya,za = scene.canvas_to_viewport(x,y)
            color = scene.trace_ray(O, np.array((xa,ya,za)))
            scene.put_pixel(x,y,color)
    
    plot_canvas(scene.canvas)

if __name__ == "__main__":
    main()