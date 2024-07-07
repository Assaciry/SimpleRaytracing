
import pygame as pg
import numpy as np
import json, sys, os
import matplotlib.pyplot as plt

from raytracing import Sphere, Scene


def plot_canvas(canvas):
    plt.figure(figsize=(8,8))
    plt.imshow(canvas)
    plt.axis("off")
    plt.show()

def main():
    Cw = 200
    Ch = 200
    Vw = 1
    Vh = 1
    f  = 1

    spheres = [Sphere(radius=0.1), Sphere(position=(2,0,4), color=(0,0,255),radius=1)]
    scene = Scene(Cw, Ch, Vw, Vh, f, spheres, BG_COLOR=(0,0,0))

    O = np.zeros(3)
    for x in np.arange(-Cw//2, Cw//2):
        for y in np.arange(-Ch//2, Ch//2):
            xa,ya,za = scene.canvas_to_viewport(x,y)
            color = scene.trace_ray(O, np.array((xa,ya,za)))
            scene.put_pixel(x,y,color)
    
    plot_canvas(scene.canvas)

if __name__ == "__main__":
    main()