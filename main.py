
# import pygame as pg
import numpy as np
import json, sys, os, time
import matplotlib.pyplot as plt

from raytracing import Sphere, Scene, Light

def plot_canvas(canvas):
    plt.figure(figsize=(8,8))
    plt.imshow(canvas,interpolation="bicubic")
    plt.axis("off")
    plt.show()

def main():
    Cw = 400
    Ch = 400
    Vw = 1
    Vh = 1
    f  = 1

    start_time = time.time()

    spheres = [Sphere(radius=0.1,specular_coeff=10,reflectance=0.2), 
               Sphere(position=(3,0,7), color=(0,0,255),radius=1,specular_coeff=20,reflectance=0.3),
               Sphere(position=(2,4,10), color=(112, 41, 99),radius=2.,specular_coeff=1000,reflectance=0.4),
               Sphere(position=(2,-4,10), color=(49, 203, 164),radius=2.,specular_coeff=1000,reflectance=0.5), 
               Sphere(color=(255,255,0), radius=5000, position=(5003,0,0))]
    lights  = [Light(intensity=0.1), Light(intensity=0.6, position=(-4,0,0)), Light(intensity=0.3,direction=(-0.5,0.5,0.7))]
    scene = Scene(Cw, Ch, Vw, Vh, f, spheres, lights, BG_COLOR=(0,0,0))

    O = np.zeros(3)
    for x in np.arange(-Cw//2, Cw//2):
        for y in np.arange(-Ch//2, Ch//2):
            xa,ya,za = scene.canvas_to_viewport(x,y)
            color = scene.trace_ray(O, np.array((xa,ya,za)),recursion_depth=3)
            scene.put_pixel(x,y,color)
    
    end_time = time.time()
    print(f"Took: {end_time-start_time:.1f} seconds to render")

    plot_canvas(scene.canvas)

if __name__ == "__main__":
    main()
    